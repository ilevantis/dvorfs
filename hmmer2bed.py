from collections import namedtuple

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


def to_nucpos(pos, frame, nuc_len):
    '''Gives the middle nucleotide of the codon (bed style index), given hmmer position (translated)'''
    if frame < 4: # forward frame
        return pos*3 + ((frame-1)%3) - 2
    else: # reverse frame
        return nuc_len + 1 - (pos*3 + ((frame-1)% 3))


BedRecord = namedtuple('BedRecord', ['chrom','start','end','name','score','strand'])

def hmmer2bed(nuc_fasta, hmmer_table, windowed=False):
    """create bed file from windowed nucleotide fasta and hmmer output table"""

    with open(nuc_fasta) as f:
        seq_lens = { i:length for i,length in fasta_lengths(f) }

    with open(hmmer_table) as f:
        rows = ( l.split() for l in f if l[0] != '#' )
        for row in rows:
            frame = int(row[0][-1])
            strand = '+' if frame < 4 else '-'

            seq_id = row[0][:-2]
            nuc_len = seq_lens[seq_id]
            e_start = int(row[19])
            e_end = int(row[20])
            range = (to_nucpos(e_start, frame, nuc_len),
                     to_nucpos(e_end, frame, nuc_len))
            if strand == '+':
                start, end = range
            else:
                end, start = range

            # expand bed coordinates to encompass whole start and end codons
            start -= 1
            end += 2

            if windowed:
                seq_name = ':'.join(seq_id.split(':')[:-1])
                offset = int(seq_id.split(':')[-1].split('-')[0])
                start += offset
                end += offset
            else:
                seq_name = seq_id

            i_e_val = float(row[12])
            qname = row[3]

            yield BedRecord(seq_name, start, end, qname, i_e_val, strand)
    return


def bed_slop(bed_record, slop, chrom_lens):
    """reimplements most basic function of bedtools slop"""
    start = max(0, bed_record.start - slop)
    end = min(chrom_lens[bed_record.chrom], bed_record.start + slop)
    return bed_record._replace(start=start, end=end)


def make_hmmer_presearchbed(nuc_fasta, hmmer_table, outbed, chrom_lens,
    slop=3000, windowed=False, e_cutoff=1e-2):
    """
    Creates a bed file with regions in which to search based on HMMER output table.
    Filters HMMER hits by evalue then applies slop to the envelope regions of the hits.
    Output bedfile is sorted
    """

    records = hmmer2bed(nuc_fasta, hmmer_table, windowed=windowed)
    records = ( r for r in records if r.score < e_cutoff ) # apply e cutoff filter
    records = ( bed_slop(r, slop, chrom_lens) for r in records ) # apply slop
    records = list(set(records)) # remove duplicate records
    # sort the bed records:
    records.sort(key=lambda r: r.start)
    records.sort(key=lambda r: r.chrom)

    with open(outbed, 'w') as f:
        for r in records:
            print(r.chrom, r.start, r.end, sep='\t', file=f)

    return outbed
