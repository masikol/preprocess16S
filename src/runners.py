# -*- coding: utf-8 -*-

import functools

import src.fastq
import src.filesystem
import src.status_bar
import src.binary_statistics
from src.platform import platf_depend_exit
import src.crosstalks.provide_crosstalk_pipe as pcp
from src.output_filenames import get_crosstalk_outfpaths
from src.printlog import printlog_info_time, printlog_error


def crosstalks_runner(args):

    if len(args.primers) != len(args.infpaths):
        printlog_error('Error: {} input fastq file(s) provided,\
 and {} primer(s) provided. Exiting...'.format(len(args.infpaths), len(args.primers)))
        platf_depend_exit(1)
    # end if

    def write_positive(fastq_records, valid_files, statistics):
        src.fastq.write_fastq_records(fastq_records, valid_files)
        statistics.increment_positive()
    # end def crosstalk_write_positive


    def write_negative(fastq_records, trash_files, statistics):
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

    crosstalk_pipe = pcp.provide_crosstalk_pipe(args.cut_off_primers)

    valid_fpaths, trash_fpaths = get_crosstalk_outfpaths(args.outdir, args.infpaths)
    valid_files = src.filesystem.open_files_for_appending(valid_fpaths)
    trash_files = src.filesystem.open_files_for_appending(trash_fpaths)

    crosstalk_statistics = src.binary_statistics.BinaryStatistics()

    crosstalk_write_positive = functools.partial(
        write_positive,
        valid_files=valid_files,
        statistics=crosstalk_statistics
    )

    crosstalk_write_negative = functools.partial(
        write_negative,
        trash_files=trash_files,
        statistics=crosstalk_statistics
    )

    loaded_crosstalk_iteration = functools.partial(
        crosstalk_iteration,
        primers=args.primers, threshold=args.threshold, max_offset=args.max_offset,
        write_positive=crosstalk_write_positive,
        write_negative=crosstalk_write_negative
    )

    src.status_bar.run_status_bar(loaded_crosstalk_iteration, args.infpaths)

    for file_collection in (valid_files, trash_files):
        src.filesystem.close_files(file_collection)
    # end for

    printlog_info_time('{} ({}%) cross-talks detected.'\
        .format(crosstalk_statistics.negative_stat,
                crosstalk_statistics.get_negative_percents())
    )
# end def crosstalks_runner
