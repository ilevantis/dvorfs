#!/usr/bin/env python3
import sys, argparse
import pandas as pd
from collections import defaultdict
import parse_genewise


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
        r1['matches'] = r1['matches'] + r2['matches'] - overlap
    else:
        r1['aaseq'] = r1['aaseq'] + ['x']*-overlap + r2['aaseq']
        r1['naseq'] = r1['naseq'] + ['nnn']*-overlap + r2['naseq']
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



def main(args):

    results = parse_genewise.parse(args.infile.name, args.fasta.name)
    df = pd.DataFrame.from_records(results, columns=parse_genewise.HSP._fields)
    df = df[df['bits'] > 0.0]

    if args.windowed:
        # parse the target seq id of window to calculate real coordinates on the contig and real contig name
        df['wstart'] = df['target'].apply(lambda x: int(x.split(':')[1].split('-')[0]))
        df['target'] = df['target'].apply(lambda x: x.split(':')[0])
        df['tstart'] = df['tstart'] + df['wstart']
        df['tend'] = df['tend'] + df['wstart']
        df = df.drop('wstart', axis=1)

    # filter / merge hits depending on filter mode
    if args.merge:
        hits = merge_hits(df, max_overlap=args.merge_overlap, max_distance=args.merge_distance)
    else:
        hits = df
        hits['no_hsps'] = 1

    if len(hits) < 1:
        filtered_hits = hits

    elif args.filter == 'all':
        filtered_hits = hits

    elif args.filter == 'best-per':
        filtered_hits = (hits
                         .groupby('target')
                         .apply(lambda x: x.nlargest(1,'bits'))
                         .sort_values('bits', ascending=False)
                         .reset_index(drop=True)
                       )

    elif args.filter == 'no-overlap':
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
    elif args.hit_mask:
        hit_mask = parse_mask_tsv(args.hit_mask)
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
        mask = filtered_hits.apply(lambda r: r.matches >= args.length_cutoff and r.bits >= args.bit_cutoff, axis=1)
        final_hits = filtered_hits[mask].reset_index(drop=True)

    # finalise output format
    cols = ['bits', 'query', 'qstart', 'qend', 'target',
            'tstrand', 'tstart', 'tend', 'no_hsps', 'matches' ]
    if args.full:
        args.aaseq = args.naseq = True
    if args.aaseq:
        final_hits['aaseq'] = final_hits['aaseq'].apply(lambda l:''.join(l))
        cols.append('aaseq')
    if args.naseq:
        final_hits['naseq'] = final_hits['naseq'].apply(lambda l:','.join(l))
        cols.append('naseq')
    if args.full:
        cols += ['overlapped', 'hsps']

    final_hits = (final_hits[cols]
            .sort_values(['target','bits'],ascending=[True,False])
            .reset_index(drop=True)
            .reset_index()
           )
    final_hits['index'] = final_hits['index'] + 1

    final_hits.to_csv(args.out, sep='\t', index=False, float_format='%.2f')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('infile',
        type=argparse.FileType('r'))

    parser.add_argument('-s', '--fasta',
        type=argparse.FileType('r'), required=True)

    parser.add_argument('-w', '--windowed',
        action='store_true',
        help="""Input fasta to genewise is windowed. Parses the seq IDs to get real coordinates.""")

    parser.add_argument('-m', '--merge',
        action='store_true',
        help="""Hits will be merged before filtering.
                (Worse hits from same query are removed at overlaps.)""")

    parser.add_argument('-d', '--merge-distance',
        type=int, default=1000)

    parser.add_argument('-o', '--merge-overlap',
        type=int, default=2)

    parser.add_argument('-f', '--filter',
        choices=['all', 'no-overlap', 'best-per'], default='no-overlap',
        help="""
        all:        All hits are kept.
        no_overlap: Hits are removed if they are overlapped by a better hit from a
                    different query.
        best_per:   Only the highest scoring hit per contig is kept.""")

    parser.add_argument('-k','--hit-mask',
        type=argparse.FileType('r'))

    parser.add_argument('-b', '--bit-cutoff',
        type=float, default=15.0)

    parser.add_argument('-l', '--length-cutoff',
        type=int, default=30)

    parser.add_argument('-a', '--aaseq',
        action='store_true',
        help="""Output predicted AA sequence.""")

    parser.add_argument('-n', '--naseq',
        action='store_true',
        help="""Output nucleotide sequence with comma seperated codons/indels.""")

    parser.add_argument('--full',
        action='store_true')

    parser.add_argument('--out',
        type=argparse.FileType('w'), default=sys.stdout)

    args = parser.parse_args()

    main(args)
