# -*- coding: utf-8 -*-
# Module contains function, which prints status bar for processing fastq file(s).

import os
import sys

import src.fastq
from src.printlog import getwt


def run_status_bar(iteration_func, infpaths):
    # Function, which prints status bar for processing fastq file(s).
    # :param iteration_func: function which manipulates fastq data;
    # :param infpaths: collection of paths to input files;
    # :type infpaths: list<str>;

    # Configure initial parameters for status bar.
    nreads = src.fastq.count_reads(infpaths)
    bar_len = _get_bar_len()
    next_print_num = int(nreads * 0.01)
    inc_num = next_print_num
    i = 0
    sys.stdout.write('{} - [>{}] 0/{} (0%)'.format(getwt(), ' '*bar_len, nreads))
    sys.stdout.flush()

    # Start processing
    for fastq_records in src.fastq.fastq_generator(infpaths):

        iteration_func(fastq_records)

        i += 1

        # Update status bar
        if i > next_print_num:
            # Get new length of terminal window
            bar_len = _get_bar_len()
            done_ratio = i / nreads
            sys.stdout.write('\r{} - [{}>{}] {}/{} ({}%)'\
                .format(getwt(), '='*int(bar_len*done_ratio),
                ' '*int(bar_len*(1-done_ratio)), i, nreads, int(done_ratio*100)))
            sys.stdout.flush()
            next_print_num += inc_num
        # end if
    # end for

    # Update status bar last time
    bar_len = _get_bar_len()
    done_ratio = i / nreads
    sys.stdout.write('\r{} - [{}{}] {}/{} ({}%)\n'\
        .format(getwt(), '='*int(bar_len*done_ratio),
        ' '*int(bar_len*(1-done_ratio)), i, nreads, int(done_ratio*100)))
    sys.stdout.flush()
    next_print_num += inc_num
# end def run_status_bar


def _get_bar_len():
    # Function returns length of status bar.
    # Firstly it tries to obtain "adaptive" bar length from terminal width.
    # It might fail if output is redirected somewhere, e.g. with `tee`.
    # If obtaining "adaptive" length fails, bar length will be set to some low value.
    # Return type: int.
    try:
        bar_len = int(os.get_terminal_size().columns * 0.40)
    except OSError:
        # Default low value
        bar_len = 40
    # end try
    return bar_len
# end _get_bar_len
