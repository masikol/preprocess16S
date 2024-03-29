# -*- coding: utf-8 -*-
# Module contains runner functions for preprocess16S tasks.

import os
import re
import glob
import functools
import subprocess as sp

import src.fastq
import src.filesystem
import src.status_bar
import src.binary_statistics
import src.output_filenames as ofn
from src.platform import platf_depend_exit
import src.crosstalks.provide_crosstalk_pipe as pcp
from src.printlog import printlog_info_time, printlog_error, printlog_info


def crosstalks_runner(args):
    # Runner function for crosstalk detection task.
    #
    # :param args: arguments for crosstalk detection task;
    # :type args: CrosstalksArguments;
    #
    # Returns two collections:
    # 1. A collection of valid (non-crosstalk) paths.
    # 2. A collection of trash (crosstalk) paths.

    def write_positive(fastq_records, valid_files, statistics):
        # Function for writing non-crosstalk fastq records.
        src.fastq.write_fastq_records(fastq_records, valid_files)
        statistics.increment_positive()
    # end def crosstalk_write_positive


    def write_negative(fastq_records, trash_files, statistics):
        # Function for writing crosstalk fastq records.
        src.fastq.write_fastq_records(fastq_records, trash_files)
        statistics.increment_negative()
    # end def crosstalk_write_positive


    def crosstalk_iteration(fastq_records, primers, threshold, max_offset,
            write_positive, write_negative):

        not_crosstalk, fastq_records = crosstalk_pipe(
            primers, fastq_records,
            threshold, max_offset
        )

        if not_crosstalk:
            write_positive(fastq_records)
        else:
            write_negative(fastq_records)
        # end if
    # end def

    printlog_info_time('Start cross-talk detection.')

    # Get crosstalk pipe
    crosstalk_pipe = pcp.provide_crosstalk_pipe(args.cut_off_primers)

    # Cnofigure and open output files
    valid_fpaths, trash_fpaths = ofn.get_crosstalk_outfpaths(args.outdir, args.infpaths)
    valid_files = src.filesystem.open_files_for_appending(valid_fpaths)
    trash_files = src.filesystem.open_files_for_appending(trash_fpaths)

    # Create BinaryStatistics object
    crosstalk_statistics = src.binary_statistics.BinaryStatistics()

    # Connect non-crosstalk output files and statistics to `write_positive` function
    crosstalk_write_positive = functools.partial(
        write_positive,
        valid_files=valid_files,
        statistics=crosstalk_statistics
    )

    # Connect crosstalk output files and statistics to `write_negative` function
    crosstalk_write_negative = functools.partial(
        write_negative,
        trash_files=trash_files,
        statistics=crosstalk_statistics
    )

    # Connect arguments to `crosstalk_iteration` function
    #   for status_bar to be wrapped around it.
    loaded_crosstalk_iteration = functools.partial(
        crosstalk_iteration,
        primers=args.primers, threshold=args.threshold, max_offset=args.max_offset,
        write_positive=crosstalk_write_positive,
        write_negative=crosstalk_write_negative
    )

    # Proceed
    src.status_bar.run_status_bar(loaded_crosstalk_iteration, args.infpaths)

    # Close files
    for file_collection in (valid_files, trash_files):
        src.filesystem.close_files(file_collection)
    # end for

    printlog_info_time('{} ({}%) cross-talks detected.'\
        .format(crosstalk_statistics.negative_stat,
                crosstalk_statistics.get_negative_percents())
    )

    return valid_fpaths, trash_fpaths
# end def crosstalks_runner


def ngmerge_runner(args):
    # Runner function for NGmerge task.
    #
    # :param args: arguments for NGmerge task;
    # :type args: NGmergeArguments;
    #
    # Returns two collections:
    # 1. A collection of valid ("merged") paths.
    # 2. A collection of trash ("unmerged") paths.

    print()
    printlog_info_time('Running NGmerge..')

    # NGmerge puts result files into working directory --
    #   we will temporarily go to output directory
    old_dir = os.getcwd()
    os.chdir(args.outdir)

    # Conigure output files' names
    merged_basename, unmerged_prefix = ofn.get_ngmerge_outprefixes(args.infpaths[0])

    # Configure command
    ngmerge_cmd = '{} -1 {} -2 {} -o {} -f {} -n {} -v -m {} -p {} -q {}'\
        .format(args.ngmerge, args.infpaths[0], args.infpaths[1],
        merged_basename, unmerged_prefix, args.n_thr,
        args.min_overlap, args.mismatch_frac, args.phred_offset)
    printlog_info('Command: `{}`'.format(ngmerge_cmd))

    # Run NGmerge
    print('NGmerge is doing it\'s job silently...')
    pipe = sp.Popen(ngmerge_cmd, shell = True, stderr=sp.PIPE)
    stderr = pipe.communicate()[1].decode('utf-8') # run NGmerge

    if pipe.returncode != 0:
        # error
        printlog_error('Error running NGmerge.: {}'.format(stderr))
        platf_depend_exit(pipe.returncode)
    # end if

    # Parse merging statistics from NGmerge's stderr
    stderr = stderr.splitlines()[1:]
    reads_pattern = r'Fragments \(pairs of reads\) analyzed: ([0-9]+)'
    merged_pattern = r'Successfully stitched: ([0-9]+)'

    # Collect statistics
    try:
        reads_processed = int(re.search(reads_pattern, stderr[0]).group(1))
        merged_reads = int(re.search(merged_pattern, stderr[1]).group(1))
    except (ValueError, AttributeError) as err:
        printlog_error('Error 78 ({}). Please, contact the developer.'.format(err))
        platf_depend_exit(78)
    # end try

    os.chdir(old_dir) # return to old dir

    printlog_info_time('NGmerge merged {}/{} ({}%) read pairs.'\
        .format(merged_reads, reads_processed,
            round(merged_reads / reads_processed * 100, 2)))

    # Configure absolute paths to output files.
    merged_fpath = os.path.join(args.outdir, merged_basename)
    unmerged_fpaths = sorted(
        glob.glob(
            os.path.join(args.outdir, '{}*.fastq'.format(unmerged_prefix))
        )
    )

    # Oh yeah, first returned value must be a collection.
    return [merged_fpath], unmerged_fpaths
# end def ngmerge_runner
