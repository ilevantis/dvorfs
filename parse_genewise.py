import re
import sys
from collections import namedtuple

codon_table = {'AAA': 'K', 'GAA': 'E', 'TAA': '*', 'CAA': 'Q',
               'AAG': 'K', 'GAG': 'E', 'TAG': '*', 'CAG': 'Q',
               'AAT': 'N', 'GAT': 'D', 'TAT': 'Y', 'CAT': 'H',
               'AAC': 'N', 'GAC': 'D', 'TAC': 'Y', 'CAC': 'H',
               'ATA': 'I', 'GTA': 'V', 'TTA': 'L', 'CTA': 'L',
               'ATG': 'M', 'GTG': 'V', 'TTG': 'L', 'CTG': 'L',
               'ATT': 'I', 'GTT': 'V', 'TTT': 'F', 'CTT': 'L',
               'ATC': 'I', 'GTC': 'V', 'TTC': 'F', 'CTC': 'L',
               'AGA': 'R', 'GGA': 'G', 'TGA': '*', 'CGA': 'R',
               'AGG': 'R', 'GGG': 'G', 'TGG': 'W', 'CGG': 'R',
               'AGT': 'S', 'GGT': 'G', 'TGT': 'C', 'CGT': 'R',
               'AGC': 'S', 'GGC': 'G', 'TGC': 'C', 'CGC': 'R',
               'ACA': 'T', 'GCA': 'A', 'TCA': 'S', 'CCA': 'P',
               'ACG': 'T', 'GCG': 'A', 'TCG': 'S', 'CCG': 'P',
               'ACT': 'T', 'GCT': 'A', 'TCT': 'S', 'CCT': 'P',
               'ACC': 'T', 'GCC': 'A', 'TCC': 'S', 'CCC': 'P' }

dna_complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N',
                  'Y':'R', 'R':'Y', 'W':'W', 'S':'S', 'K':'M',
                  'M':'K', 'D':'H', 'V':'B', 'H':'D', 'B':'V' }


def revcomp(dna_seq):
    return ''.join([dna_complement[n.upper()] for n in dna_seq[::-1]])


def translate_codon(triplet):
    try:
        aa = codon_table[triplet]
    except KeyError:
        aa = 'X'
    return aa


ptrns = { 'tab_header': r'Bits.*introns\n',
          'pep_header': r'>.*\.pep\n',
          'pretty_header': r'\s{21}Alignment.*\n',
          'hit_header': r'>Results for ([^ ]+) vs ([^ ]+) \((reverse|forward)\) \[\d+\]\n',
          'alb_line': r'([0-9\.-]+) \[(.*)\],\[(.*)\]\n' }
ptrns = {k:re.compile(v) for k,v in ptrns.items()}


def SimpleFastaParser(handle):
# adapted from biopython Bio.SeqIO.FastaIO SimpleFastaParser
    l = next(handle)
    try:
        while True:
            title = l[1:].rstrip().split(' ')[0]
            lines = []
            l = next(handle)
            while l[0] != ">":
                lines.append(l.rstrip())
                l = next(handle)
            yield title, "".join(lines).replace(" ", "").replace("\r", "")
    except StopIteration:
        yield title, "".join(lines).replace(" ", "").replace("\r", "")
        return

def get_fasta_seq(seqid, handle):
    try:
        seq = ''
        sfp = SimpleFastaParser(handle)
        while seq == '':
            i,s = next(sfp)
            if i == seqid:
                seq = s
    except StopIteration:
        print('Could not find seq id in fasta!', file=sys.stderr)
    return seq



AlState = namedtuple('AlState', ['score', 'qpos', 'qstate', 'tpos', 'tstate'])

def parse_alb_line(l):
    s, q, t = ptrns['alb_line'].fullmatch(l).groups()
    s = float(s)
    # add one to each position to get correct index
    qpos = tuple(int(i)+1 for i in q.split(' ')[0].split(':'))
    qstate = q.split(' ')[1][1:-1]
    tpos = tuple(int(i)+1 for i in t.split(' ')[0].split(':'))
    tstate = t.split(' ')[1][1:-1]
    return AlState(s,qpos,qstate,tpos,tstate)

HSP = namedtuple('HSP', ['query', 'qstart', 'qend',
                         'target', 'tstrand', 'tstart', 'tend',
                         'bits', 'matches', 'naseq', 'aaseq' ])


def parse(genewisefilename, fastafilename, low_mem=False):
    '''
    Generator that iterates through a genewise output file (genewise must be
    called with '-alb'). Also requires the fasta file that genewise was used on.
    Outputs an HSP namedtuple for each alignment in the genewise file.
    Coordinates are BED-style indexed.
    low_mem: does not load fastafile into memory (Takes much much longer!)
    '''
    handle = open(genewisefilename, 'r')
    l = next(handle)

    if not low_mem:
        with open(fastafilename) as f:
            rawseqs = {i:s for i,s in SimpleFastaParser(f)}
    else:
        prev_target = ''

    try:
        while True:
            if l[:8] == '>Results':
                query, target, direction = ptrns['hit_header'].fullmatch(l).groups()
                if direction == 'forward':
                    strand = '+'
                else:
                    strand = '-'

                if not low_mem :
                    rawseq = rawseqs[target]
                else:
                    if target != prev_target:
                        # get fasta sequence from fasta file
                        with open(fastafilename) as f:
                            rawseq = get_fasta_seq(target, f)
                        prev_target = target

                if strand == '-':
                    seq = revcomp(rawseq)
                else:
                    seq = rawseq.upper()

                l = next(handle)


            elif ptrns['alb_line'].fullmatch(l):
                a = parse_alb_line(l)
                aaseq = []
                naseq = []
                score = 0.0
                matches = 0
                tstart = a.tpos[0]
                qstart = a.qpos[0]

                while (    (a.tstate, a.qstate) != ('RANDOM_SEQUENCE', 'LOOP')
                       and (a.tstate, a.qstate) != ('END', 'END')
                      ):

                    if a.tpos[0] >= len(seq):
                        print("Fasta sequence '{}' shorter than expected!".format(target), file=sys.stderr)
                    elif (a.tstate, a.qstate) == ('CODON', 'MATCH_STATE'):
                        # codon
                        s = seq[a.tpos[0]:a.tpos[1]]
                        naseq.append(s)
                        aaseq.append(translate_codon(s))
                        score += a.score
                        matches += 1
                    elif (a.tstate, a.qstate) == ('CODON', 'INSERT_STATE'):
                        # target has insertion compared to query
                        s = seq[a.tpos[0]:a.tpos[1]]
                        naseq.append(s)
                        aaseq.append(translate_codon(s))
                        score += a.score
                    elif (a.tstate, a.qstate) == ('INSERT', 'DELETE_STATE'):
                        # target has deletion compared to query
                        score += a.score
                        matches += 1
                    elif (a.tstate, a.qstate) == ('SEQUENCE_DELETION', 'INSERT_STATE'):
                        # likely insertion frameshift mutation between two alignable codons
                        s = seq[a.tpos[0]:a.tpos[1]]
                        naseq.append(s)
                        aaseq.append('!')
                        score += a.score
                    elif (a.tstate, a.qstate) == ('SEQUENCE_INSERTION', 'INSERT_STATE'):
                        # likely insertion frameshift mutation between two alignable codons
                        s = seq[a.tpos[0]:a.tpos[1]]
                        naseq.append(s)
                        aaseq.append('!')
                        score += a.score
                    elif (a.tstate, a.qstate) == ('SEQUENCE_DELETION', 'MATCH_STATE'):
                        # likely deletion frameshift mutation that disrupted an alignable codon
                        s = seq[a.tpos[0]:a.tpos[1]]
                        naseq.append(s)
                        aaseq.append('?')
                        score += a.score
                        matches += 1
                    elif (a.tstate, a.qstate) == ('SEQUENCE_INSERTION', 'MATCH_STATE'):
                        # likely insertion frameshift and highly diverged codon (but don't know what codon)
                        s = seq[a.tpos[0]:a.tpos[1]]
                        naseq.append(s)
                        aaseq.append('&')
                        score += a.score
                        matches += 1
                    else:
                        #should not get here
                        print("{} vs {} ({})".format(query,target,direction), file=sys.stderr)
                        print(a, file=sys.stderr)

                    l = next(handle)
                    a = parse_alb_line(l)

                # outside while loop means we hit random_sequence or end state
                tend = a.tpos[0]
                qend = a.qpos[0]
                if strand == '-':
                    seqlen = len(seq)
                    tstart, tend = seqlen - tend, seqlen - tstart

                hsp = HSP(query, qstart, qend,
                          target, strand, tstart, tend,
                          score, matches, naseq, aaseq)

                if hsp.matches > 0: #ignore empty hsps caused by random seq -> end transition
                    yield hsp

                l = next(handle)

            else:
                l = next(handle)

    except StopIteration:
        pass




###############################################################################

# OLD CODE

def parse_tab_line(line):
    cols = ['bits', 'query', 'qstart', 'qend', 'target', 'tstart', 'tend', 'idels', 'introns']
    d = dict(zip(cols,line.strip().split()))
    for k,v in d.items():
        if k in ['bits', 'qstart', 'qend', 'tstart', 'tend', 'idels', 'introns']:
            d[k] = literal_eval(v)

    if d['tstart'] < d['tend']:
        d['tstrand'] = '+'
    else:
        d['tstart'], d['tend'] = d['tend'], d['tstart']
        d['tstrand'] = '-'

    return d

def parse_nonalb(handle):
    '''Coordinates are GFF-style indexed.'''
    from ast import literal_eval
    ptrns = { 'tab_header': r'Bits.*introns\n',
              'pep_header': r'>.*\.pep\n',
              'pretty_header': r'\s{21}Alignment.*\n' }
    ptrns = {k:re.compile(v) for k,v in ptrns.items()}
    summs, peps, qstrings, tstrings  = [], [], [], []

    l = next(handle)
    try:
        while True:
            if ptrns['tab_header'].fullmatch(l):
                l = next(handle)
                while l != '//\n':
                    d = parse_tab_line(l)
                    summs.append(d)
                    l = next(handle)

            elif ptrns['pep_header'].fullmatch(l):
                l = next(handle)
                pep = ''
                while l[0] not in ['>','/']:
                    pep += l.strip()
                    l = next(handle)
                peps.append(pep)

            elif ptrns['pretty_header'].fullmatch(l):
                l = next(handle)
                ali_block = []
                while not any([ptrns['pretty_header'].fullmatch(l), l == '//\n']):
                    ali_block.append(l)
                    l = next(handle)
                qstring = ''.join([s[21:].rstrip() for s in ali_block[2::8]])
                tstring = ''.join([s[21:].rstrip() for s in ali_block[4::8]])
                qstrings.append(qstring)
                tstrings.append(tstring)

            else:
                l = next(handle)

    except StopIteration:
        pass

    combine = lambda x: { **x[0], **{'pep':x[1]}, **{'qstring':x[2]}, **{'tstring':x[3]} }
    results = [ combine(i) for i in zip(summs, peps, qstrings, tstrings) ]
    return results
