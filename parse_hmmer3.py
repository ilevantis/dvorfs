#!/usr/bin/env python3
import re
from ast import literal_eval
import numpy as np


def parse_val(s):
    try:
        return literal_eval(s)
    except:
        if s in ['nan', '-nan']:
            return np.nan
        elif s in ['inf']:
            return np.inf
        else:
            return s

def parse_score_line(line):
    split = line.split()
    col_vals = split[0:3] + split[8:9] + [' '.join(split[9:])]
    return { 'hit_eval' : float(col_vals[0]),
             'hit_score': float(col_vals[1]),
             'hit_bias' : float(col_vals[2]),
             'hit_name' : col_vals[3],
             'hit_desc' : col_vals[4]         }

_DOMAIN_COLS = ['hsp_num', 'hsp_score', 'hsp_bias', 'hsp_c_eval', 'hsp_i_eval', 'hsp_q_start', 'hsp_q_end', 'hsp_h_start', 'hsp_h_end', 'hsp_e_start', 'hsp_e_end', 'hsp_acc']
def parse_domain_line(line):
    split = line.split()
    col_vals = split[0:1] + split[2:8] + split[9:11] + split[12:14] + split[15:1]
    return { i[0]:parse_val(i[1]) for i in zip(_DOMAIN_COLS, col_vals) }


_RECORD_SPLIT_STRING = 'Internal pipeline statistics summary:\n-------------------------------------\n'
_QUERY_NAME_PATTERN = re.compile(r'Query:(.*)\[M=\d+\]\n')
_SCORES_TABLE_PATTERN = re.compile(r'Scores for complete sequences \(score includes all domains\):.*-----------\n(.*)\nDomain annotation for each sequence \(and alignments\):',
    re.DOTALL )
_HSP_PARENT_PATTERN = re.compile(r'^([^\s]*)\s+')
_HSP_NUM_PATTERN = re.compile(r'domain (\d+)\s')
_ALIGN_PATTERN = re.compile(r'  .*\s+\d+\s+([a-zA-Z\.]+)')

def parse_record(record):
    """take in a record as a list of lines, return a tuple:
       ('query_name','dict_of_hits')"""
    query = _QUERY_NAME_PATTERN.search(record[0]).group(1).strip()

    if  '   [No hits detected that satisfy reporting thresholds]\n' in record:
        vals = []

    else:
        hits = {}
        record = ''.join(record)
        results, stats = record.split(_RECORD_SPLIT_STRING)

        scores_table = _SCORES_TABLE_PATTERN.search(results).group(1).splitlines()
        for line in scores_table:
            if line != '  ------ inclusion threshold ------':
                hit = parse_score_line(line)
                hits[hit['hit_name']] = hit

        hits_hsps = results.split('>> ')[1:]
        for hsp in hits_hsps:
            hsp_dict = {}
            hsp_parent = _HSP_PARENT_PATTERN.search(hsp).group(1)
            hsp_data_parts = hsp.split('Alignments for each domain:')

            domain_table = ( s for s in hsp_data_parts[0].splitlines()[3:] if s.strip() )
            for line in domain_table:
                hsp_info = parse_domain_line(line)
                hsp_dict[hsp_info['hsp_num']] = hsp_info

            alignments = ( s.splitlines() for s in hsp_data_parts[1].split('==')[1:] )
            for alignment in alignments:
                hsp_num = literal_eval(_HSP_NUM_PATTERN.search(alignment[0]).group(1))
                first_line = 1
                for a in alignment[1:]:
                    if _ALIGN_PATTERN.search(a):
                        break
                    else:
                        first_line += 1
                align_span = _ALIGN_PATTERN.search(alignment[first_line]).span(1)
                hsp_dict[hsp_num]['alignstr_q'] = alignment[first_line][slice(*align_span)]
                hsp_dict[hsp_num]['alignstr_m'] = alignment[first_line+1][slice(*align_span)]
                hsp_dict[hsp_num]['alignstr_h'] = alignment[first_line+2][slice(*align_span)]
                hsp_dict[hsp_num]['alignstr_a'] = alignment[first_line+3][slice(*align_span)]

            hits[hsp_parent]['hsps'] = [v for v in hsp_dict.values()]

        vals = [v for v in hits.values()]

    return (query, vals)



def parse_hmmout(handle):
    record = []
    for line in handle:
        if line[0] in ['#', '\n']:
            pass
        elif line != '//\n':
            record.append(line)
        else:
            yield parse_record(record)
            record = []
