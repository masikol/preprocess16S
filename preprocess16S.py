#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = '5.0.a'
# Year, month, day
__last_update_date__ = '2021-XX-XX'
__author__ = 'Maxim Sikolenko'


import os
import sys
import argparse
import itertools


import src.runners
from src.platform import platf_depend_exit
from src.arguments import CrosstalksArguments
from src.crosstalks.get_primers_seqs import get_primers_seqs
from src.printlog import config_logging, printlog_info_time, getwt



def handle_args():
    parser = argparse.ArgumentParser(description="""
  preprocess16S (Version {}. {} edition. Author: {}).
  A tool for preprocessing 16S amplicon reads."""\
.format(__version__, __last_update_date__, __author__))

    parser.add_argument('-v', '--version', help='print version and exit',
        action='store_true')
    parser.add_argument('-1', '--forw-fpath', help='fastq file of forward reads (required)',
        required=True)
    parser.add_argument('-2', '--revr-fpath', help='fastq file of reverse reads (optional)')
    parser.add_argument('-p', '--primers-fpath',
        help='fasta file of primer sequences (default: Illimina V3-V4 amplicon primers)')
    parser.add_argument('-o', '--outdir', help='output directory',
        default=os.path.join(os.getcwd(), 'preprocess16S-outdir'))
    parser.add_argument('-x', '--threshold', help='percent identity threshold (default: 0.52)',
        type=float, default=0.52)
    parser.add_argument('-s', '--max-offset', help='maximum offset (default: 2)',
        type=int, default=2)
    parser.add_argument('-c', '--cut-off-primers',
        help='whether to cut off primers [0, 1] (default: 1, i.e. do cut them off)',
        type=bool, default=True)

    # Shameful ad hoc thing
    if '-v' in sys.argv[1:] or '--verbose' in sys.argv[1:]:
        print(__version__)
        platf_depend_exit(1)
    # end if

    args = parser.parse_args()

    return args
# end def handle_args


def main():
    args = handle_args()
    config_logging(args.outdir)

    crosstalks_args = CrosstalksArguments(
        tuple(
            itertools.takewhile(
                lambda x: not x is None,
                [args.forw_fpath, args.revr_fpath]
            )
        ),
        get_primers_seqs(args.primers_fpath),
        args.threshold,
        args.max_offset,
        args.cut_off_primers,
        args.outdir
    )
    src.runners.crosstalks_runner(crosstalks_args)

    printlog_info_time('End.')
# end def main


if __name__ == '__main__':
    main()
# end if
