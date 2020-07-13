# -*- coding: utf-8 -*-
# Module for reading and writing fastq files.

import sys
from src.printing import *
from src.filesystem import *


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


def fastq_read_packets(read_paths, reads_at_all, n_thr):
    """
    Function-generator for retrieving FASTQ records from PE files
        and distributing them evenly between 'n_thr' processes
        for further parallel processing.
        # end for
    :param read_paths: dictionary (dict<str: str> of the following structure:
    {
        "R1": path_to_file_with_forward_reads,
        "R2": path_to_file_with_reverse_reads
    }
    :param reads_at_all: number of read pairs in these files;
    :type reads_at_all: int;
    Yields lists of FASTQ-records (structure of these records is described in 'write_fastq_record' function).
    Returns None when end of file(s) is reached.
    """

    how_to_open = OPEN_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]
    fmt_func = FORMATTING_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]
    read_files = dict() # dictionary for file objects

    # Open files that contain reads meant to be processed.
    for key, path in read_paths.items():
        read_files[key] = how_to_open(path)
    # end for

    # Compute packet size (one packet -- one thread).
    pack_size = reads_at_all // n_thr
    if reads_at_all % n_thr != 0:
        pack_size += 1
    # end if

    for i in range(n_thr):
        packet = list()
        for j in range(pack_size):
            tmp = read_fastq_pair(read_files, fmt_func)
            # if the end of the line is reached
            if tmp is None:
                # Yield partial packet and return 'None' on next call
                yield packet
                return
            else:
                packet.append(tmp)
            # end if
        # end for
        yield packet # yield full packet
    # end for
# end def fastq_read_packets