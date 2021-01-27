#!/usr/bin/env python3
import parse_hmmer3
import argparse, sys
from collections import defaultdict


def parse_fai(handle):
    with handle as f:
        lengths = { i.split('\t')[0] : int(i.split('\t')[1]) for i in f }
    return lengths


def fasta_lengths(handle):
    l = next(handle)
    try:
        while True:
            title = l[1:].rstrip().split(' ')[0]
            length = 0
            l = next(handle)
            while l[0] != ">":
                length += len(l.rstrip().replace(" ", "").replace("\r", ""))
                l = next(handle)
            yield title, length
    except StopIteration:
        yield title, length
        return

def hsp_filter(hit, e_cutoff):
    hit['hsps'] = [ hsp for hsp in hit['hsps'] if hsp['hsp_i_eval'] < e_cutoff ]
    return hit


def print_header(seq_lens, windowed=False):
    if windowed:
        lengths = defaultdict(lambda:0)
        for k,v in seq_lens.items():
            seq_id = ':'.join(k.split(':')[:-1])
            offset = int(k.split(':')[-1].split('-'))
            length = offset + v
            lengths[seq_id] = max(lengths[seq_id],length)
        lengths = dict(lengths)
    else:
        lengths = seq_lens

    print('##gff-version 3')
    for k,v in lengths.items():
        print('##sequence-region',k,1,v,sep=' ')


def to_nucpos(pos, frame, nuc_len):
    '''Gives the middle nucleotide of the codon'''
    if frame < 4: # forward frame
        return pos*3 + (frame - 1) % 3 - 1
    else: # reverse frame
        return nuc_len + 1 - (pos*3 + (frame - 1) % 3) + 1


def print_hits(hits, seq_lens, windowed=False):
    for query, hits in hits.items():
        for hit in hits:
            frame = int(hit['hit_name'][-1])
            strand = '+' if frame < 4 else '-'
            seq_id = hit['hit_name'][:-2]
            nuc_len = seq_lens[seq_id]
            for hsp in hit['hsps']:
                range = (to_nucpos(hsp['hsp_e_start'],frame, nuc_len),
                         to_nucpos(hsp['hsp_e_end'],frame, nuc_len))
                if strand == '+':
                    start, end = range
                else:
                    end, start = range
                start = 1 if start < 2 else start-1
                end =   nuc_len if end > nuc_len-1 else end+1

                if windowed:
                    seq_name = ':'.join(seq_id.split(':')[:-1])
                    offset = int(seq_id.split(':')[-1])
                    start += offset
                    end += offset
                else:
                    seq_name = seq_id

                line = [seq_name,
                        'hmmer',
                        'alignment',
                        start,
                        end,
                        hsp['hsp_i_eval'],
                        strand,
                        '.',
                        'Name='+query ]
                print(*line, sep='\t', flush=True)

def main(args):

    seq_lens = { i:length for i,length in fasta_lengths(args.fasta) }

    print_header(seq_lens,windowed=args.windowed)
    sys.stdout.flush()

    with args.hits_file as f:
        hits = { k : [ hsp_filter(hit, args.e_cutoff) for hit in v ] for k,v in parse_hmmer3.parse_hmmout(f) if v }

    hits = { k : [ hit for hit in v if hit['hsps'] ] for k,v in hits.items() }
    hits = { k : v for k,v in hits.items() if v }

    print_hits(hits, seq_lens, windowed=args.windowed)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("hits_file", type=argparse.FileType('r'))
    parser.add_argument("-f", "--fasta", type=argparse.FileType('r'), required=True)
    parser.add_argument('-e', '--e_cutoff', type=float, default=1e-2)
    parser.add_argument("-w", "--windowed", action='store_true')

    args = parser.parse_args()
    main(args)
