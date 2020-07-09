# -*- coding: utf-8 -*-
# Module for reading and writing fastq files.

import sys
from src.printing import *


def read_fastq_pair(read_files, fmt_func):
    """
    Function reads pair of FASTQ records from two FASTQ files: R1 and R2,
       assumming that these FASTQ files are sorted.

    :param read_files: a dictionary of the following structure:
    {
        "R1": R1 file object,
        "R2": R2 file object
    };
    :type read_files: dict<str: _io.TextIOWrapper or 'gzip.GzipFile>;
    :param fmt_func: a function from 'FORMATTING_FUNCS' tuple;

    Returns a dictionary of the following structure:
    {
        "R1": R1_dict,
        "R1": R1_dict
    }
    'R1_dict' and 'R2_dict' are dictionaries of the following structure:
    {
        "seq_id": ID if the sequence,
        "seq": actually sequence,
        "opt_id": third line of FASTQ record,
        "qual_str": quality line
    }
    I.e. type of returned value is 'dict< dict<str, str> >'
    """

    if len(read_files) != 2 and len(read_files) != 1:
        print_error("You can only pass 1 or 2 files to the function 'read_fastq_pair'!\a")
        print("Contact the developer")
        sys.exit(1)
    # end if

    fastq_recs = dict()           # this dict should consist of two fastq-records: R1 and R2
    for key in read_files.keys():
        fastq_recs[key] = {                    #read all 4 lines of fastq-record
            "seq_id": fmt_func(read_files[key].readline()),
            "seq": fmt_func(read_files[key].readline()).upper(), # searching for cross-talks is case-dependent
            "opt_id": fmt_func(read_files[key].readline()),
            "qual_str": fmt_func(read_files[key].readline())
        }

        # If all lines from files are read
        if fastq_recs[key]["seq_id"] == "":
            return None
        # end if
    # end for
    return fastq_recs
# end def read_fastq_pair


def write_fastq_record(outfile, fastq_record):
    """
    Function writes FASTQ record to 'outfile'.

    :param outfile: file object that describes file to write in;
    :type outfile: _io.TextIOWrapper;
    :param fastq_record: dict of 4 elements. Elements are four corresponding lines of FASTQ file.
       Structure of this dictionary if following:
    {
        "seq_id": ID if the sequence,
        "seq": actually sequence,
        "opt_id": third line of FASTQ record,
        "qual_str": quality line
    }
    :type fastq_record: dict<str: str>
    """

    try:
        outfile.write(fastq_record["seq_id"] + '\n')
        outfile.write(fastq_record["seq"] + '\n')
        outfile.write(fastq_record["opt_id"] + '\n')
        outfile.write(fastq_record["qual_str"] + '\n')
        outfile.flush()
    except Exception as exc:
        print_error("\n error while writing to output file")
        print( str(exc) )
        sys.exit(1)
    # end try
# end def write_fastq_record
