# -*- coding: utf-8 -*-
# Mudule contains function, which prints help message.

def print_help(version, last_update_date):
    # Function, which prints help message.
    #
    # :param version: program version;
    # :type version: str;
    # :param last_update_date: last update version of program;
    # :type last_update_date: str;

    print(' <<< preprocess16S >>>')
    print(' Version {}. {} edition.'.format(version, last_update_date))
    print('-'*10)
    print('# Basic usage is:')
    print('  ./preprocess16S.py --tasks task1,task2 -1 forward_R1_reads.fastq.gz \
-2 reverse_R2_reads.fastq.gz -o outdir/')
    print('\nSee details at https://github.com/masikol/preprocess16S')
    print('-'*10)
    print('# Options:')
    print('Print-and-exit options:')
    print('  -h (--help) -- show help message.')
    print('  -v (--version) -- show version.')
    print('## General:')
    print("""  --tasks -- comma-separated list of tasks to run.
    Permitted values: rm-crosstalks, ngmerge.
    Default: rm-crosstalks.""")
    print('  -1 (--R1) -- FASTQ file of forward reads.')
    print('  -2 (--R2) -- FASTQ file of reverse reads.')
    print('  -o (--outdir) -- output directory.')
    print("""  -z (--gzip-output) -- [0, 1]. 1 -- gzip output files after work is done.
    0 -- keep output files uncompressed.
    Default: 1.""")
    print('## Crosstalks detection:')
    print("""  -r (--primers) -- FASTA file, where primers sequences are stored (one line per sequence).
    Illumina V3-V4 primer sequences are used by default.""")
    print("""  -x (--threshold) -- threshold value used in crosstalks detection;
    Real number from 0 to 1. Default: 0.52.""")
    print("""  -s (--max-offset) -- maximum offset used in crosstalks detection;
    Integer > 0. Default: 2.""")
    print("""  -с (--cut-off-primers) [0, 1]. 0 -- keep primers, 1 -- cut them off.
    Default -- 1.""")
    print('## Read merging:')
    print("""  -m (--min-overlap) -- minimum overlap of the paired-end reads to be merged with NGmerge.
    Default: 20 nt.""")
    print("""  -p (--mismatch-frac) -- fraction of mismatches to allow in the overlapped region
    (a fraction of the overlap length).
    Default: 0.1.""")
    print("""  -t (--threads) <int> -- number of threads to launch.
    Default: 1.""")
    print("""  -q (--phred-offset) [33, 64] -- Phred quality offset.
    Default: 33.""")
    print("""  --ngmerge-path -- path to NGmerge executable.
    You can specify it if bundled NGmerge 0.3 is not suitable for you.""")
# end def print_help
