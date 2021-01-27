#!/usr/bin/env python3
from subprocess import Popen, PIPE
import sys, os, shutil
from os import path

###

debug = False
script_dir = path.dirname(path.abspath(__file__))
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
    if exit_code != 0:
        raise

    return output, err, exit_code


def hmmconversion(hmminput, type):

    hmm2db = path.join(workdir,'hmmdb.hmm2')
    hmm3db = path.join(workdir,'hmmdb.hmm3')

    if type == 'seed':
        cmd = ['hmmbuild2', hmm2db, hmminput]
        p1 = Popen(cmd)
        cmd = ['hmmbuild', hmm3db, hmminput]
        p2 = Popen(cmd)
        exit_codes = (p1.wait(), p2.wait())

    elif type == 'hmm2':
        cmd = ['hmmconvert', hmminput]
        with open(hmm3db, 'w') as f:
            p = Popen(cmd, stdout=f)
            exit_code = p.wait()
        hmm2db = hmminput

    elif type == 'hmm3':
        cmd = ['hmmconvert', hmminput]
        with open(hmm2db, 'w') as f:
            p = Popen(cmd, stdout=f)
            exit_code = p.wait()
        hmm3db = hmminput

    else:
        raise

    return hmm2db, hmm3db


def hmmer_pre_search(dbfafai, dbfasta, hmm3db,
                        threads=2):

    cmds = [
        ['bedtools', 'makewindows', '-g', dbfafai, '-w', '20000', '-s', '15000'],
        ['bedtools', 'getfasta', '-fi', dbfasta, '-bed', '-'],
        ['sed', '/^>/ s/:/\$/'],
        ['fold', '-w', '80']
    ]

    window_fa_name = path.join(workdir,'windows.fa')
    with open(window_fa_name,'w') as f:
        exit_info = run_pipeline(cmds, f)


    tr_fa_name = path.join(workdir,'tr.faa')
    cmd = ['transeq', '-frame', '6', '-alternative', window_fa_name, tr_fa_name]
    p = Popen(cmd)
    exit_code = p.wait()
    cmd = ['sed', '-i', '/^>/ s/\$/:/', window_fa_name, tr_fa_name]
    p = Popen(cmd)
    exit_code = p.wait()

    hmmer_out_name = path.join(workdir,'hmmer.out')
    cmd = ['hmmsearch', '--notextw', '--cpu', str(int(threads)), '-o', hmmer_out_name, hmm3db, tr_fa_name]
    p = Popen(cmd)
    exit_code = p.wait()


    cmds = [
        [path.join(script_dir,'hmm2gff.py'), '-w', '-e', '1e5', '-f', window_fa_name,  hmmer_out_name],
        ['gawk', 'BEGIN{FS=OFS="\t"}/^#/{next}{print $1, $4 - 1, $5}'],
        ['sort', '-k', '1,1', '-k2,2n'],
        ['bedtools', 'slop', '-g', dbfafai, '-b', '3000']
    ]

    hmmer_bed_name = path.join(workdir,'hmmer.bed')
    with open(hmmer_bed_name, 'w') as f:
        exit_info = run_pipeline(cmds, f)

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


def run_genewise(sefasta, hmm2db, hit_mask=None,
                    threads=2, mem=4e6):

    cmd = ['genewisedb', '-pfam', hmm2db, '-dnadb', sefasta,
            '-alg', '333', '-aalg', '333L', '-init', 'local',
            '-gap', '6', '-ext', '0', '-subs', '1e-2', '-indel', '1e-4',
            '-cut', '15',
            '-pthread', '-pthr_no', str(int(threads)), '-kbyte', str(int(mem)),
            '-aln', '9999999', '-alb', '-quiet']

    gw_out = path.join(workdir,'gw.out')
    with open(gw_out, 'w') as f:
        p = Popen(cmd, stdout=f)
        exit_code = p.wait()

    cmd = [path.join(script_dir,'process_genewise.py'),
            gw_out, '--fasta', sefasta, '--windowed',
            '--merge', '--merge-distance', '1000',
            '--filter', 'no-overlap',
            '--bit-cutoff', '15.0', '--length-cutoff', '20',
            '--full']

    if hit_mask:
        cmd += ['--hit-mask', hit_mask]

    gw_tsv = path.join(workdir,'gw.tsv')
    with open(gw_tsv, 'w') as f:
        p = Popen(cmd, stdout=f)
        exit_code = p.wait()

    return gw_tsv


###


def main(args):

    # set work directory for temp files
    global workdir
    workdir = path.join(args.workdir,'dvorfs.tmp')
    os.mkdirs(workdir, exist_ok=True)

    # convert HMM files to both formats
    if args.seed:
        hmm2db, hmm3db = hmmconversion('seed', args.seed)
    elif args.hmm2:
        hmm2db, hmm3db = hmmconversion('hmm2', args.hmm2)
    elif args.hmm3:
        hmm2db, hmm3db = hmmconversion('hmm3', args.hmm3)


    if not args.bed:
        # do a presearch
        search_bed = hmmer_pre_search(args.fai, args.fasta, hmm3db, threads=args.procs)
        shutil.copy(search_bed,path.join(args.outdir,'dvorfs_presearch.bed'))
    else:
        # use supplied bed file
        search_bed = args.bed

    # chop into 20kbp windows to stop GW from breaking
    windowed_fasta = bed2wdw_fa(search_bed, args.fasta)

    # run GW and process hits
    gw_tsv = run_genewise(windowed_fasta, hmm2db, hit_mask=None, threads=args.procs, mem=4e6)
    shutil.copy(gw_tsv,path.join(args.outdir,'dvorfs_hits.tsv'))

    # clean up the work directory
    if not args.keep_workdir:
        shutil.rmtree(workdir)






###

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--fasta',
        type=str, required=True)

    parser.add_argument('-i', '--fai',
        type=str, required=True)

    hmmarg = parser.add_mutually_exclusive_group(required=True)
    hmmarg.add_argument('--hmm2', type=str)
    hmmarg.add_argument('--hmm3', type=str)
    hmmarg.add_argument('--seed', type=str)

    parser.add_argument('-b', '--bed',
        type=str)

    parser.add_argument('-p', '--procs',
        type=int, default=2)

    parser.add_argument('-o','--outdir',
        type=str, default='./')
    parser.add_argument('-d','--workdir',
        type=str, default='./')
    parser.add_argument('-k', '--keep-workdir',
        action='store_true')


    main(args)
