# -*- coding: utf-8 -*-
# Module contains functions for manipulating fastq. And a container class for fastq.

import src.compression
from src.platform import platf_depend_exit
from src.printlog import printlog_error, printlog_info_time


class FastqRecord:
    # Class represents a fastq record.
    # Fields:
    # :read_name: read name WITH preceding '@';
    # :seq: sequence;
    # :comment: comment string;
    # :qual_str: quality string;

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
        # If first line doesn't look like read name:
        if self.read_name[0] != '@':
            return 'Invalid first line of fastq record: `{}`'\
                .format(self.read_name)
        # If length of sequence and length of quality string doesn't equal:
        elif len(self.seq) != len(self.qual_str):
            return 'Length of sequence of is not equal to length of quality string. Read: `{}`'\
                .format(self.read_name)
        else:
            return None
        # end if
    # end def validate_fastq

    def __str__(self):
        return '{}\n{}\n{}\n{}\n'.format(self.read_name, self.seq, self.comment, self.qual_str)
# end class FastqRecord


def fastq_generator(fq_fpaths):
    # Function yields fastq records.
    # It does not create new FastqRecord object each time.
    # Instead it just updates extant object.
    # :param fq_fpaths: list ot paths to input fastq files;
    # :type fq_fpaths: list<str>, tuple<str>;
    # Yields list of FastqRecord-s, list<FastqRecord>.

    # Get open funtions for both files
    open_funcs = src.compression.provide_open_funcs(fq_fpaths)

    # Open input files and create FastqRecord objects for forward and reverse reads.
    fq_files = list()
    fq_records = list()
    for fpath, open_func in zip(fq_fpaths, open_funcs):
        fq_files.append(open_func(fpath))
        fq_records.append(FastqRecord(None, None, None, None))
    # end for

    eof = False

    while not eof:

        for fq_record, fq_file in zip(fq_records, fq_files):
            # Update FastqRecord
            fq_record.update_record(
                fq_file.readline().strip(),
                fq_file.readline().strip(),
                fq_file.readline().strip(),
                fq_file.readline().strip()
            )
        # end for

        if fq_records[0].read_name == '':
            eof = True # end of file
        else:
            # Validate fastq record(s)
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

    # Close input files.
    for fq_file in fq_files:
        fq_file.close()
    # end for
# end def fastq_generator


def write_fastq_records(fq_records, outfiles):
    # Function writes forward read to "forward" output file
    #   and reverse read to "reverse" output file.
    # :param fq_records: collection of fastq records to write;
    # :type fq_records: list<FastqRecord>;
    # :param outfiles: collection of output files;
    # :type outfiles: list<_io.TextIOWrapper>, list<gzip.GzipFile>, list<bz2.BZ2File>;
    for fq_record, outfile in zip(fq_records, outfiles):
        outfile.write(str(fq_record))
    # end for
# end def write_fastq_records


def count_reads(fq_fpaths):
    # Function counts reads in a fastq file.
    # :param fq_fpaths: list of paths to fastq files;
    # :type fq_fpaths: list<str>;
    # Returns number of reads in "forward" file
    #   (assuming that there are as many reads in "reverse" file).

    printlog_info_time('Counting reads...')

    # Get open function for "forward" file
    open_func = src.compression.provide_open_funcs(fq_fpaths)[0]

    # Count reads in "forward" file
    with open_func(fq_fpaths[0]) as infile:
        nreads = sum(1 for _ in infile) // 4
    # end with

    printlog_info_time('{} reads.'.format(nreads))

    return nreads
# end def count_reads
