# -*- coding: utf-8 -*-

import sys

import src.compression
from src.printlog import printlog_error, printlog_info_time
from src.platform import platf_depend_exit

class FastqRecord:

    def __init__(self, read_name=None, seq=None, comment=None, qual_str=None):
        self.read_name = read_name
        self.seq = seq
        self.comment = comment
        self.qual_str = qual_str
    # end def __init__

    def update_record(self, read_name, seq, comment, qual_str):
        self.read_name = read_name
        self.seq = seq
        self.comment = comment
        self.qual_str = qual_str
    # end def update_record

    def validate_fastq(self):
        if self.read_name[0] != '@':
            return 'Invalid first line of fastq record: `{}`'\
                .format(self.read_name)
        elif len(self.seq) != len(self.qual_str):
            return 'Length of sequence of is not equal to length of quality string. Read: `{}`'\
                .format(self.read_name)
        else:
            return None
        # end if
    # end def validate_fastq

    def __repr__(self):
        return '\n<{}\n{}...\n{}\n{}...>\n'\
            .format(self.read_name, self.seq[:20], self.comment, self.qual_str[:20])
    # end def __repr__

    def __str__(self):
        return '{}\n{}\n{}\n{}\n'.format(self.read_name, self.seq, self.comment, self.qual_str)
# end class FasrqRecord


def fastq_generator(fq_fpaths):

    open_funcs = src.compression.provide_open_funcs(fq_fpaths)

    fq_files = list()
    fq_records = list()
    for fpath, open_func in zip(fq_fpaths, open_funcs):
        fq_files.append(open_func(fpath))
        fq_records.append(FastqRecord(None, None, None, None))
    # end for

    eof = False

    while not eof:

        for fq_record, fq_file in zip(fq_records, fq_files):
            fq_record.update_record(
                fq_file.readline().strip(),
                fq_file.readline().strip(),
                fq_file.readline().strip(),
                fq_file.readline().strip()
            )
        # end for

        if fq_records[0].read_name == '':
            eof = True
        else:
            for fq_record in fq_records:
                error_response = fq_record.validate_fastq()
                if not error_response is None:
                    printlog_error('Fastq error: {}'.format(error_response))
                    platf_depend_exit(1)
                # end if
            # end for
            yield fq_records
        # end if
    # end while

    for fq_file in fq_files:
        fq_file.close()
    # end for
# end def fastq_generator


def write_fastq_records(fq_records, outfiles):
    for fq_record, outfile in zip(fq_records, outfiles):
        outfile.write(str(fq_record))
    # end for
# end def write_fastq_records


def count_reads(fq_fpaths):

    printlog_info_time('Counting reads...')

    open_func = src.compression.provide_open_funcs(fq_fpaths)[0]

    with open_func(fq_fpaths[0]) as infile:
        nreads = sum(1 for _ in infile) // 4
    # end with

    printlog_info_time('{} reads.'.format(nreads))

    return nreads
# end def count_reads
