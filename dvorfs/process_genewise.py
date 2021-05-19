#!/usr/bin/env python3
import sys, argparse, os
import pandas as pd
import numpy as np
from collections import defaultdict
from . import parse_genewise


# NOTE for sequence arithmetic: records use BED-style index

def target_overlap(h1,h2): # uses bed index
    overlap = min(h1['tend'],h2['tend']) - max(h1['tstart'],h2['tstart'])
    return overlap


def query_overlap(h1,h2): # uses bed index
    overlap = min(h1['qend'],h2['qend']) - max(h1['qstart'],h2['qstart'])
    return overlap


def correct_order(h1,h2):
    q_order = sorted([h1,h2], key=lambda x: x['qstart'])
    if h1['tstrand'] == '+':
        t_order = sorted([h1,h2], key=lambda x: x['tstart'])
    elif h1['tstrand'] == '-':
        t_order = sorted([h1,h2], key=lambda x: x['tend'], reverse=True)

    if q_order == t_order:
        return True
    else:
        return False


def simplify_record(r):
    return {k:v for k,v in r.items() if k not in ['target','tstrand','query']}


def combine(r1,r2):
    overlap = query_overlap(r1,r2)
    if overlap > 0:
        r1['aaseq'] = r1['aaseq'][:-overlap] + r2['aaseq']
        r1['naseq'] = r1['naseq'][:-overlap] + r2['naseq']
        r1['alpos'] = r1['alpos'][:-overlap] + r2['alpos']
        r1['matches'] = r1['matches'] + r2['matches'] - overlap
    else:
        r1['aaseq'] = r1['aaseq'] + ['x']*-overlap + r2['aaseq']
        r1['naseq'] = r1['naseq'] + ['nnn']*-overlap + r2['naseq']
        r1['alpos'] = (r1['alpos']
            + list(zip(range(r1['alpos'][-1][0]+1, r2['alpos'][0][0]),[0]*-overlap))
            + r2['alpos'])
        r1['matches'] = r1['matches'] + r2['matches']
    r1['qstart'] = min([r1['qstart'],r2['qstart']])
    r1['qend'] = max([r1['qend'],r2['qend']])
    r1['tstart'] = min([r1['tstart'],r2['tstart']])
    r1['tend'] = max([r1['tend'],r2['tend']])
    r1['bits'] = r1['bits'] + r2['bits']
    r1['overlapped'] += r2['overlapped']
    r1['hsps'].append(simplify_record(r2))
    return r1


def merge_hits(df, max_overlap=0, max_distance=0):
    merged = []

    for g_tup, g in df.groupby(['target','tstrand','query'], sort=False):
        target, tstrand, query = g_tup

        # remove lower scoring overlapping hits
        stack = g.sort_values('bits', ascending=False).to_dict(orient='records')
        no_overlaps = []
        while len(stack) > 0:
            r1 = stack.pop(0)
            r1['target'], r1['tstrand'], r1['query'] = g_tup
            r1['overlapped'] = 0
            new_stack = []
            for r2 in stack:
                if target_overlap(r1,r2) > 0:
                    # r2 is a worse hit and not useful, remove it
                    r1['overlapped'] += 1
                else:
                    # keep r2
                    new_stack.append(r2)
            stack = new_stack
            no_overlaps.append(r1)

        # sort by position and merge
        if tstrand == '+':
            merge_stack = sorted(no_overlaps, key=lambda x: x['tstart'])
        else:
            merge_stack = sorted(no_overlaps, key=lambda x: x['tend'], reverse=True)

        while len(merge_stack) > 0:
            r1 = merge_stack.pop(0)
            r1['hsps'] = [simplify_record(r1)]
            new_merge_stack = []
            for i, r2 in enumerate(merge_stack):
                if query_overlap(r1,r2) > max_overlap:
                    # don't merge but keep r2
                    new_merge_stack.append(r2)
                elif -1*target_overlap(r1,r2) > max_distance:
                    # don't merge but keep all following r2s
                    new_merge_stack += merge_stack[i:]
                    break
                elif not correct_order(r1,r2):
                    # query pieces in wrong order:
                    # don't merge but keep all following r2s
                    new_merge_stack += merge_stack[i:]
                    break
                else:
                    # no overlaps and not too far between r1 and r2
                    # so combine them and remove r2 from the stack
                    r1 = combine(r1,r2)

            merge_stack = new_merge_stack

            r1['no_hsps'] = len(r1['hsps'])
            merged.append(r1)


    all_merged_hits = pd.DataFrame.from_dict(merged)
    return all_merged_hits


def parse_mask_tsv(f):
    d = defaultdict(list)
    for l in f:
        k, *v = l.strip().split('\t')
        if len(v) < 2:
            d[k].append((0, 1000000))
        else:
            d[k].append((int(v[0]), int(v[1])))
    return dict(d)


def make_ali_array(hit_df):
    """
    Make a numpy array of the codon alignment with final row as codon reference row.
    hit_df must contain the columns: id, alpos, naseq.
    Outputs the array and a list of the hit ids representing the row-wise order of sequences in the
    alignment.
    """

    col_dict = defaultdict(lambda: np.full((len(hit_df)+1,1), '---', dtype=object))

    # Mark all the canoncial cols
    for c in range(hit_df['qstart'].min(), hit_df['qend'].max()+1):
        col_dict[(c,0)][-1,0] = 'N--'

    seq_order = []
    # Fill in the codon sequences in the cols
    for row, hit in enumerate(hit_df.itertuples()):
        seq_order.append(hit.id)
        for pos,s in zip(hit.alpos, hit.naseq):
            col_dict[pos][row,0] = s
            # if its non-canonical col label it
            if pos[1] > 0:
                col_dict[pos][-1,0] = 'n--'

    col_items = sorted(list(col_dict.items()), key=lambda i: i[0])
    cols = [ i[1] for i in col_items ]

    # pad each column to mod 3 length
    for c in cols:
        width = max(len(i[0]) for i in c)
        width += -width % 3
        for i in range(len(c)):
            c[i,0] += '-'*(width-len(c[i,0]))

    # remove empty cols
    cols = [ c for c in cols if not np.all(c[:-1] == '---') ]

    ali_arr = np.hstack(cols)

    return ali_arr, seq_order



def process_genewise(infile, fasta, windowed=False,
    merge=False, merge_distance=1000, merge_overlap=2,
    filter_type='no-overlap', hit_mask=None, bit_cutoff=15.0, length_cutoff=30,
    out_cols=[], make_alis=False):

    results = parse_genewise.parse(infile, fasta)
    df = pd.DataFrame.from_records(results, columns=parse_genewise.HSP._fields)
    df = df[df['bits'] > 0.0]

    if windowed:
        # parse the target seq id of window to calculate real coordinates on the contig and real contig name
        df['wstart'] = df['target'].apply(lambda x: int(x.split(':')[1].split('-')[0]))
        df['target'] = df['target'].apply(lambda x: x.split(':')[0])
        df['tstart'] = df['tstart'] + df['wstart']
        df['tend'] = df['tend'] + df['wstart']
        df = df.drop('wstart', axis=1)

    # filter / merge hits depending on filter mode
    if merge:
        hits = merge_hits(df, max_overlap=merge_overlap, max_distance=merge_distance)
    else:
        hits = df
        hits['no_hsps'] = 1

    if len(hits) < 1:
        filtered_hits = hits

    elif filter_type == 'all':
        filtered_hits = hits

    elif filter_type == 'best-per':
        filtered_hits = (hits
                         .groupby('target')
                         .apply(lambda x: x.nlargest(1,'bits'))
                         .sort_values('bits', ascending=False)
                         .reset_index(drop=True)
                       )

    elif filter_type == 'no-overlap':
        keep = []
        for g_tup, g in hits.groupby(['target','tstrand'], sort=False):
            stack = list(g.sort_values('bits', ascending=False).itertuples())
            while len(stack) > 0:
                r1 = stack.pop(0)
                new_stack = []
                for r2 in stack:
                    if min(r1.tend,r2.tend) - max(r1.tstart,r2.tstart) <= 0:
                        new_stack.append(r2)
                stack = new_stack
                keep.append(r1.Index)

        filtered_hits = hits.iloc[keep]

    # apply hit_mask
    if len(filtered_hits) < 1:
        pass
    elif hit_mask:
        with open(hit_mask) as f:
            hit_mask = parse_mask_tsv(f)
        def mask_hits(r):
            masked = False
            regions = hit_mask.get(r.query,[])
            for start,end in regions:
                if r.qstart >= start and r.qend <= end:
                    masked = True
            return not masked

        mask = filtered_hits.apply(mask_hits, axis=1)
        filtered_hits = filtered_hits[mask]

    # apply quality filters
    if len(filtered_hits) < 1:
        final_hits = filtered_hits
    else:
        mask = filtered_hits.apply(lambda r: r.matches >= length_cutoff and r.bits >= bit_cutoff, axis=1)
        final_hits = filtered_hits[mask].reset_index(drop=True)

    final_hits = (final_hits
            .sort_values(['query','bits'],ascending=[True,False])
            .reset_index(drop=True)
            .reset_index()
           )
    final_hits['id'] = final_hits['index'] + 1

    # finalise output format
    cols = ['id', 'bits', 'query', 'qstart', 'qend', 'target',
            'tstrand', 'tstart', 'tend', 'no_hsps', 'matches' ]

    cols += out_cols
    output_df = final_hits[cols]

    if make_alis:
        alis = []
        for name, hit_df in final_hits.groupby(['query']):
            alis.append( (name, *make_ali_array(hit_df)) )
    else:
        alis = None

    return output_df, alis



def main():
    parser = argparse.ArgumentParser(description=f"process_genewise.py is part of DVORFS")

    parser.add_argument('infile',
        type=argparse.FileType('r'),
        help="""Genewise output file. Genewise must have been run with the `-alb` argument.""")

    parser.add_argument('-s', '--fasta',
        type=argparse.FileType('r'), required=True,
        help="""Exact fasta file that genewise was run with.""")

    parser.add_argument('-w', '--windowed',
        action='store_true',
        help="""Input fasta to genewise is windowed. Parses the seq IDs to get real coordinates.""")

    parser.add_argument('-m', '--merge',
        action='store_true',
        help="""Hits will be merged before filtering.
                (Worse hits from same query are removed at overlaps.)""")

    parser.add_argument('-d', '--merge-distance',
        type=int, default=1000,
        help="""Maximum allowed distance between two hits in the subject sequence (in bp) for them to be merged.""")

    parser.add_argument('-o', '--merge-overlap',
        type=int, default=2,
        help="""Maximum number of positions overlapping in the query pHMM for two hits to be merged.""")

    parser.add_argument('-f', '--filter',
        choices=['all', 'no-overlap', 'best-per'], default='no-overlap',
        help="""
        all:        All hits are kept.
        no_overlap: Hits are removed if they are overlapped by a better hit from a
                    different query.
        best_per:   Only the highest scoring hit per contig is kept.""")

    parser.add_argument('-k','--hit-mask',
        type=argparse.FileType('r'),
        help="""TSV file with 3 columns: 1. name of query, 2. start position of masked region, 3. end position of masked region.""")

    parser.add_argument('-b', '--bit-cutoff',
        type=float, default=15.0,
        help="""Minimum bit score for hits to be kept after merging.""")

    parser.add_argument('-l', '--length-cutoff',
        type=int, default=30,
        help="""Minimum number of codons aligned to the query pHMM for hits to be kept after merging.""")

    parser.add_argument('-a', '--aaseq',
        action='store_true',
        help="""Output predicted AA sequence.""")

    parser.add_argument('-n', '--naseq',
        action='store_true',
        help="""Output nucleotide sequence with comma seperated codons/indels.""")

    parser.add_argument('--full',
        action='store_true')

    parser.add_argument('--out',
        type=argparse.FileType('w'), default=sys.stdout,
        help="""Specify output file. By default, output goes to STDOUT.""")

    parser.add_argument('--aliout',
        help="""Output an alignment fasta of hits for each HMM with any hits into specified dir.""")


    args = parser.parse_args()

    outcols = []
    if args.full:
        args.aaseq = args.naseq = True
    if args.aaseq:
        outcols.append('aaseq')
    if args.naseq:
        outcols.append('naseq')
    if args.full:
        outcols  += ['overlapped', 'hsps']

    make_alis = True if args.aliout else False
    hit_mask = args.hit_mask.name if args.hit_mask else None

    df, alis = process_genewise(args.infile.name, args.fasta.name, windowed=args.windowed,
        merge=args.merge, merge_distance=args.merge_distance, merge_overlap=args.merge_overlap,
        filter_type=args.filter, hit_mask=hit_mask, bit_cutoff=args.bit_cutoff, length_cutoff=args.length_cutoff,
        out_cols=outcols, make_alis=make_alis)

    if args.aaseq:
        df['aaseq'] = df['aaseq'].apply(lambda l:''.join(l))
    if args.naseq:
        df['naseq'] = df['naseq'].apply(lambda l:','.join(l))

    df.to_csv(args.out, sep='\t', index=False, float_format='%.2f')

    if args.aliout:
        os.makedirs(args.aliout, exist_ok=True)

        for name, ali, order in alis:
            with open(os.path.join(args.aliout, f'{name}.ali.fa'), 'w') as f:
                print(">CODONS", file=f)
                print(''.join(ali[-1]), file=f)

                for id, row in zip(order,ali[:-1]):
                    print(f">hitid-{id}", file=f)
                    print(''.join(row), file=f)
