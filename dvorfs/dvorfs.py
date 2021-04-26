#!/usr/bin/env python3

__version__ = '1.0.0'

from subprocess import Popen, PIPE, DEVNULL
import sys, os, shutil, argparse
from os import path
from .hmmer2bed import make_hmmer_presearchbed
from .process_genewise import process_genewise

###

debug = False
workdir = ''

###

def run_pipeline(cmd_list, stdout):
    procs = []
    for i, cmd in enumerate(cmd_list):
        if i == 0:
            procs.append(Popen(cmd, stdin=None, stdout=PIPE))
        elif i == len(cmd_list)-1:
            procs.append(Popen(cmd, stdin=procs[-1].stdout, stdout=stdout, stderr=PIPE))
        else:
            procs.append(Popen(cmd, stdin=procs[-1].stdout, stdout=PIPE))


    for i,p in enumerate(procs[:-1]):
        if debug:
            print(' '.join(cmd_list[i]), file=sys.stderr)
        p.stdout.close()

    if debug:
        print(' '.join(cmd_list[-1]), file=sys.stderr)

    output, err = procs[-1].communicate()
    exit_code = procs[0].wait()
    if exit_code != 0: raise

    return output, err, exit_code


def clean_and_exit(keep_workdir=False):
    global workdir
    if not keep_workdir:
        print("Cleaning up working directory.", file=sys.stderr)
        shutil.rmtree(workdir)
    exit()



def hmmconversion(hmminput, hmm2db, hmm3db, type):

    if type == 'seed':
        cmd = ['hmmbuild2', hmm2db, hmminput]
        p1 = Popen(cmd, stdout=DEVNULL)
        cmd = ['hmmbuild', hmm3db, hmminput]
        p2 = Popen(cmd, stdout=DEVNULL)
        exit_codes = (p1.wait(), p2.wait())
        for ec in exit_codes:
            if ec != 0: raise

    elif type == 'hmm2':
        cmd = ['hmmconvert', hmminput]
        with open(hmm3db, 'w') as f:
            p = Popen(cmd, stdout=f)
            exit_code = p.wait()
            if exit_code != 0: raise
        shutil.copy(hmminput, hmm2db)

    elif type == 'hmm3':
        cmd = ['hmmconvert', '-2', hmminput]
        with open(hmm2db, 'w') as f:
            p = Popen(cmd, stdout=f)
            exit_code = p.wait()
            if exit_code != 0: raise
        shutil.copy(hmminput, hmm3db)

    else:
        raise

    return hmm2db, hmm3db



def hmmer_pre_search(dbfafai, dbfasta, hmm3db,
    threads=2, slop=3000, e_cutoff=1e-2):

    # window the input fasta and replace ':' with '$' in the fasta IDs
    # so transeq doesn't mess with the IDs
    cmds = [
        ['bedtools', 'makewindows', '-g', dbfafai, '-w', '20000', '-s', '15000'],
        ['bedtools', 'getfasta', '-fi', dbfasta, '-bed', '-'],
        ['sed', '/^>/ s/:/\$/'],
        ['fold', '-w', '80']
    ]

    window_fa_name = path.join(workdir,'pre_windows.fa')
    with open(window_fa_name,'w') as f:
        exit_info = run_pipeline(cmds, f)


    # translate the windows using EMBOSS tools transeq
    tr_fa_name = path.join(workdir,'tr.faa')
    cmd = ['transeq', '-frame', '6', '-alternative', window_fa_name, tr_fa_name]
    p = Popen(cmd)
    exit_code = p.wait()
    if exit_code != 0: raise


    # replace '$' with ':' again in both nucleotide and translated fasta
    cmd = ['sed', '-i', '/^>/ s/\$/:/', window_fa_name, tr_fa_name]
    p = Popen(cmd)
    exit_code = p.wait()
    if exit_code != 0: raise


    # run HMMER3 on the translated windows
    hmmer_out_name = path.join(workdir,'hmmer.out')
    cmd = ['hmmsearch', '--cpu', str(int(threads)), '--domtblout', hmmer_out_name, hmm3db, tr_fa_name]
    p = Popen(cmd, stdout=DEVNULL)
    exit_code = p.wait()
    if exit_code != 0: raise


    with open(dbfafai) as f:
        chrom_lens = { l.split('\t')[0]:int(l.split('\t')[1]) for l in f }

    # make the presearch bed based on HMMER output
    hmmer_bed_name = path.join(workdir,'hmmer.bed')
    make_hmmer_presearchbed(window_fa_name, hmmer_out_name, hmmer_bed_name, chrom_lens,
        slop=slop, windowed=True, e_cutoff=e_cutoff)

    return hmmer_bed_name



def bed2wdw_fa(bedfile, dbfasta):

    cmds = [
        ['bedtools', 'merge', '-d', '200', '-i', bedfile],
        ['bedtools', 'makewindows', '-b', '-', '-w', '20000', '-s', '15000'],
        ['gawk', 'BEGIN{FS=OFS="\t"}$3-$2>=90{print $0}'], # filter out windows < 30AA
        ['bedtools', 'getfasta', '-fi', dbfasta, '-bed', '-'],
        ['fold', '-w', '80']
    ]

    window_fa_name = path.join(workdir,'windows.fa')
    with open(window_fa_name,'w') as f:
        exit_info = run_pipeline(cmds, f)

    return window_fa_name


def genome2wdw_fa(dbfafai, dbfasta):

    cmds = [
        ['bedtools', 'makewindows', '-g', dbfafai, '-w', '20000', '-s', '15000'],
        ['gawk', 'BEGIN{FS=OFS="\t"}$3-$2>=90{print $0}'], # filter out windows < 30AA
        ['bedtools', 'getfasta', '-fi', dbfasta, '-bed', '-'],
        ['fold', '-w', '80']
    ]

    window_fa_name = path.join(workdir,'windows.fa')
    with open(window_fa_name,'w') as f:
        exit_info = run_pipeline(cmds, f)

    return window_fa_name



def run_genewise(sefasta, hmm2db, gw_out, threads=2, mem=4e6):

    cmd = ['genewisedb', '-pfam', hmm2db, '-dnadb', sefasta,
            '-alg', '333', '-aalg', '333L', '-init', 'local',
            '-gap', '6', '-ext', '0', '-subs', '1e-2', '-indel', '1e-4',
            '-cut', '15',
            '-pthread', '-pthr_no', str(int(threads)), '-kbyte', str(int(mem)),
            '-aln', '9999999', '-alb', '-quiet']

    # make an environment with the WISECONFIGDIR environment variable so genewise is happy
    # assumes we are running in a conda environment with wise2 installed through conda
    gw_env = os.environ.copy()
    gw_env['WISECONFIGDIR'] = path.join(os.environ['CONDA_PREFIX'], 'share/wise2/wisecfg/')

    gw_out_temp = path.join(workdir,'gw.out.temp')
    with open(gw_out_temp, 'w') as f:
        p = Popen(cmd, stdout=f, env=gw_env)
        exit_code = p.wait()
        if exit_code != 0: raise
    os.rename(gw_out_temp, gw_out)

    return gw_out



def process_genewise_output(sefasta, gw_out,
    make_alis=False, out_cols=[], **kwargs):

    # output file/dir names
    gw_tsv = path.join(workdir,'gw.tsv')
    gw_alidir = path.join(workdir,'gw_ali')

    df, alis = process_genewise(gw_out, sefasta,
        make_alis=make_alis, out_cols=out_cols, **kwargs)

    if 'aaseq' in out_cols:
        df['aaseq'] = df['aaseq'].apply(lambda l:''.join(l))
    if 'naseq' in out_cols:
        df['naseq'] = df['naseq'].apply(lambda l:','.join(l))

    df.to_csv(gw_tsv, sep='\t', index=False, float_format='%.2f')

    if make_alis:
        os.makedirs(gw_alidir, exist_ok=True)

        for name, ali, order in alis:
            with open(path.join(gw_alidir, f'{name}.ali.fa'), 'w') as f:
                print(">CODONS", file=f)
                print(''.join(ali[-1]), file=f)

                for id, row in zip(order,ali[:-1]):
                    print(f">hitid-{id}", file=f)
                    print(''.join(row), file=f)

    return gw_tsv, gw_alidir



###


def run_dvorfs(args):

    # set work directory for temp files
    global workdir
    if args.workdir:
        workdir = args.workdir
    else:
        workdir = path.join(args.outdir,f'{args.prefix}.tmp')
    os.makedirs(workdir, exist_ok=True)

    if not args.fai:
        args.fai = args.fasta + '.fai'

    # convert HMM files to both formats
    hmm2db = path.join(workdir,'hmmdb.hmm2')
    hmm3db = path.join(workdir,'hmmdb.hmm3')
    if not (path.isfile(hmm2db) and path.isfile(hmm3db)):
        print("Converting HMM databases...", file=sys.stderr)
        if args.seed:
            hmmconversion(args.seed, hmm2db, hmm3db, type='seed')
        elif args.hmm2:
            hmmconversion(args.hmm2, hmm2db, hmm3db, type='hmm2')
        elif args.hmm3:
            hmmconversion(args.hmm3, hmm2db, hmm3db, type='hmm3')
    else:
        print("Converted HMM databases already present, skipping step.", file=sys.stderr)

    if not args.full_search:
        if not args.bed:
            # do a presearch
            print("Running presearch using HMMER...", file=sys.stderr)
            search_bed = hmmer_pre_search(args.fai, args.fasta, hmm3db, threads=args.procs, slop=args.presearch_slop)
            shutil.copy(search_bed,path.join(args.outdir,f'{args.prefix}.presearch.bed'))

            if os.stat(search_bed).st_size == 0:
                print("Presearch found no potential hits. Finishing early.", file=sys.stderr)
                clean_and_exit(args.keep_workdir)

        else:
            # use supplied bed file
            search_bed = args.bed

        # chop search regions into 20kbp windows to stop GW from breaking
        windowed_fasta = bed2wdw_fa(search_bed, args.fasta)

    else:
        # chop whole genome into 20kbp windows to stop GW from breaking
        windowed_fasta = genome2wdw_fa(args.fai, args.fasta)

    # run GW
    gw_out = path.join(workdir,'gw.out')
    if not path.isfile(gw_out):
        # complete genewise output does not exist - run genewise
        print("Running GeneWise...", file=sys.stderr)
        run_genewise(windowed_fasta, hmm2db, gw_out, threads=args.procs, mem=4e6)
    else:
        print("GeneWise outfile already present, skipping step.", file=sys.stderr)

    # Process hits
    out_cols = ['aaseq']
    if args.full_tsv:
        args.nuc_tsv = True
    if args.nuc_tsv:
        out_cols.append('naseq')
    if args.full_tsv:
        out_cols  += ['overlapped', 'hsps']

    hit_mask = args.hit_mask.name if args.hit_mask else None

    print("Processing GeneWise hits...", file=sys.stderr)
    gw_tsv, gw_alidir = process_genewise_output(windowed_fasta, gw_out, windowed=True,
        merge=(not args.nomerge), merge_distance=args.merge_distance,
        filter_type=args.filter, hit_mask=hit_mask, bit_cutoff=args.bit_cutoff, length_cutoff=args.length_cutoff,
        out_cols=out_cols, make_alis=(not args.noali))

    shutil.copy(gw_tsv, path.join(args.outdir,f'{args.prefix}.hits.tsv'))

    if not args.noali:
        shutil.copytree(gw_alidir, path.join(args.outdir,f'{args.prefix}.alis'), dirs_exist_ok=True)


    print("DVORFS finsihed running.", file=sys.stderr)

    # clean up the work directory
    clean_and_exit(args.keep_workdir)






###

def main():
    parser = argparse.ArgumentParser(description=f"DVORFS v{__version__}")

    queryargs = parser.add_argument_group("Input target sequence")
    queryargs.add_argument('-f', '--fasta',
        type=str, required=True,
        help="""Input fasta file.""")

    queryargs.add_argument('-i', '--fai',
        type=str,
        help="""Input fasta index file (default: <input.fasta>.fai).""")

    hmmargs = parser.add_argument_group("Input query profiles")
    hmmargs = hmmargs.add_mutually_exclusive_group(required=True)
    hmmargs.add_argument('--hmm2', type=str,
    help="""Input query hmm file (HMMER2 file format).""")
    hmmargs.add_argument('--hmm3', type=str,
    help="""Input query hmm file (HMMER3 file format).""")
    hmmargs.add_argument('--seed', type=str,
    help="""Input query seed alignment file (stockholm or fasta format).""")

    presearchargs = parser.add_argument_group("Presearch")
    presearchargs.add_argument('--presearch-slop',
        type=int, default=3000,
        help="""Size of flanking regions next to presearch hits in which to search (default: 3000).""")

    presearchargs.add_argument('-b', '--bed',
        type=str,
        help="""Skip presearch step and use a bed file to limit search regions.""")

    presearchargs.add_argument('--full-search',
        action='store_true',
        help="""Skip presearch step and run genewise on the whole input fasta.
                WARNING: this can take a very long time.""")

    runargs = parser.add_argument_group("Runtime/output")
    runargs.add_argument('-p', '--procs',
        type=int, default=2,
        help="""Number of processor threads to use for running HMMER and GeneWise.""")

    runargs.add_argument('-o','--outdir',
        type=str, default='./',
        help="""Output directory. Defaults to current working directory.""")

    runargs.add_argument('--prefix',
        type=str, default='dvorfs',
        help="""Prefix for output files (default: "dvorfs").""")

    runargs.add_argument('-d','--workdir',
        type=str,
        help="""Directory in which DVORFS will save files during a run (default: <outdir>/dvorfs.tmp).""")

    runargs.add_argument('-k', '--keep-workdir',
        action='store_true',
        help="""Do not delete the temporary working directory after DVORFS has finished.""")

    runargs.add_argument('--nuc-tsv',
        action='store_true',
        help="""Nucleotide sequences with comma seperated codons are inlcuded in the output tsv.""")

    runargs.add_argument('--full-tsv',
        action='store_true',
        help="""Extra columns containing postprocessing details are included in the output tsv.""")

    runargs.add_argument('--noali',
        action='store_true',
        help="""Do not output explicit codon alignments of hits.""")

    postprocargs = parser.add_argument_group("Hit postprocessing and filtering")
    postprocargs.add_argument('--nomerge',
        action='store_true',
        help="""Adjacent GeneWise hits will not be merged.""")

    postprocargs.add_argument('--merge-distance',
        type=int, default=1000,
        help="""Maximum allowed distance between GeneWise hits to be merged.""")

    postprocargs.add_argument('--filter',
        choices=['all', 'no-overlap', 'best-per'], default='all',
        help="""
        all:        All hits are kept (default).
        no_overlap: Hits are removed if they are overlapped by a better hit from a
                    different query.
        best_per:   Only the highest scoring hit per contig is kept.""")

    postprocargs.add_argument('--hit-mask',
        type=argparse.FileType('r'),
        help="""Three column TSV file for filtering out hits based on query regions.""")

    postprocargs.add_argument('--bit-cutoff',
        type=float, default=15.0,
        help="""Bit score threshold for reported hits (default: 15.0).""")

    postprocargs.add_argument('--length-cutoff',
        type=int, default=30,
        help="""Length threshold (no. of codons) for reported hits (default: 30).""")

    parser.add_argument('-v','--version', action='version', version=f"DVORFS v{__version__}")


    args = parser.parse_args()
    run_dvorfs(args)
