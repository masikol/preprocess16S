#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Version 2.1;
# 09.09.2019 edition

# |===== Check python interpreter version. =====|

from sys import version_info as verinf

if verinf.major < 3:#{
    print( "\nYour python interpreter version is " + "%d.%d" % (verinf.major, verinf.minor) )
    print("\tPlease, use Python 3!\a")
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python 2.
    raw_input("Press ENTER to exit:")
    exit(1)
#}


from sys import stdout as sys_stdout
def printn(text):#{
    """
    Function prints text to the console without adding '\n' in the end of the line.
    Why not just to use 'print(text, end="")'?
    In order to display informative error message if Python 2.X is launched
        instead if awful error traceback.
    """
    sys_stdout.write(text)
#}


def print_error(text):#{
    """Function for printing pretty error messages"""
    print("\n   \a!! - ERROR: " + text + '\n')
#}

# |===== Stuff for dealing with time =====|

from time import time, strftime, localtime, gmtime
start_time = time()


def get_work_time():#{
    return strftime("%H:%M:%S", gmtime( time() - start_time))
#}

# |===========================================|


# |===== Handle CL arguments =====|

import os
import getopt
from sys import argv
import multiprocessing as mp

usage_msg = """
DESCRIPTION:\n
preprocess16S.py -- this script is designed for preprocessing reads from 16S regions of rDNA.\n
It works with Illumina pair-end reads.\n
It detects and removes reads, that came from other sample (aka cross-talks),
    relying on the information, whether there are PCR primer sequences in these reads or not.\n
More precisely, if required primer sequence is found at the 5'-end of a read, therefore,
    it is a read from 16S rDNA and we need it.\n
Moreover, it can cut these primers off and merge reads by using 'read_merging_16S' module.\n
Reads should be stored in files, both of which are of FASTQ format or both are gzipped (i.e. fastq.gz).
Sequences of required primers are retrieved from FASTA file (or fasta.gz).\n
Attention! This script cannot be executed by python interpreter version < 3.0!
----------------------------------------------------------\n
OPTIONS:\n
    -h (--help) --- show help message;\n
    -c (--cutoff) --- Flag option. If specified, primer sequences will be cut off;\n
    -m (--merge-reads) --- Flag option. If specified, reads will be merged together 
        with 'read_merging_16S' module;\n
    -q (--quality-plot) --- Flag option. If specified, a graph of read quality distribution
        will be plotted. Requires 'numpy' and 'matplotlib' to be installed;\n
    -p (--primers) --- FASTA file, in which primer sequences are stored;\n
    -1 (--R1) --- FASTQ file, in which forward reads are stored;\n
    -2 (--R2) --- FASTQ file, in which reverse reads are stored;\n
    -o (--outdir) --- directory, in which result files will be placed;\n
    -t (--threads) <int> --- number of threads to launch;
----------------------------------------------------------\n
EXAMPLES:\n
  1) Filter cross-talks from files 'forw_R1_reads.fastq.gz' and 'rev_R2_reads.fastq.gz' depending on
     primer sequenes in file 'V3V4_primers.fasta'. Cut primer sequences off and put result files in
     the directory named 'outdir':\n
       ./preprocess16S.py -p V3V4_primers.fasta -1 forw_R1_reads.fastq.gz -2 rev_R2_reads.fastq.gz -o outdir -c\n
  2) Filter cross-talks from files 'forw_R1_reads.fastq.gz' and 'rev_R2_reads.fastq.gz' depending on
     primer sequenes in file 'V3V4_primers.fasta'. Cut primer sequences off, merge reads,
     plot a quality graph and put result files in the directory named 'outdir'. Launch 4 threads to calculations:\n
       ./preprocess16S.py -p V3V4_primers.fasta -1 forw_R1_reads.fastq.gz -2 rev_R2_reads.fastq.gz -o outdir -cmq -t 4\n
"""

try:#{
    opts, args = getopt.getopt(argv[1:], "hcmqp:1:2:o:t:", 
        ["help", "cutoff", "merge-reads", "quality-plot", "primers=", "R1=", "R2=", "outdir=", "threads="])
#}
except getopt.GetoptError as opt_err:#{
    print( str(opt_err) )
    print(usage_msg)
    exit(2)
#}

def check_file_existance(path):#{
    if os.path.exists(path):#{
        return
    #}
    else:#{
        print_error("File '{}' does not exist!".format(path))
        exit(1)
    #}
#}

cutoff = False
merge_reads = False
quality_plot = False
primer_path = None
outdir_path = "{}{}preprocess16S_result_{}".format(os.getcwd(), os.sep, start_time).replace(" ", "_") # default path
read_paths = dict()
n_thr = 1

for opt, arg in opts:#{

    if opt in ("-h", "--help"):#{
        print(usage_msg)
        exit(0)
    #}

    elif opt in ("-c", "--cutoff"):#{
        cutoff = True
    #}

    elif opt in ("-m", "--merge-reads"):#{
        merge_reads = True
    #}

    elif opt in ("-q", "--quality-plot"):#{
        quality_plot = True
    #}

    elif opt in ("-p", "--primers"):#{
        check_file_existance(arg)
        primer_path = arg
    #{

    elif opt in ("-o", "--outdir"):#{
        outdir_path = os.path.abspath(arg)
    #}

    elif opt in ("-1", "--R1"):#{
        if not "-2" in argv and not "--R2" in argv:
            print_error("You should specify both forward and reverse reads!")
            exit(1)
        check_file_existance(arg)
        read_paths["R1"] = arg
    #{

    elif opt in ("-2", "--R2"):#{
        if not "-1" in argv and not "--R1" in argv:
            print_error("You should specify both forward and reverse reads!")
            exit(1)
        check_file_existance(arg)
        read_paths["R2"] = arg
    #}

    elif opt in ("-t", "--threads"):#{
        try:#{
            n_thr = int(arg)
            if n_thr < 1:
                raise ValueError
        #}
        except ValueError:#{
            print_error("number of threads must be positive integer number!")
            print(" And here is your value: '{}'".format(arg))
            exit(1)
        #}
        if n_thr > mp.cpu_count():#{
            print("\n  Warning! {} CPU threads are available".format(mp.cpu_count()))
            print("And you have specified {} threads to use.\n".format(n_thr))
            reply = input("""If this is just what you want, press ENTER
    or enter 'q' to exit:>>""")
            if reply == "":#{
                pass
            #}
            else:#{
                exit(0)
            #}
        #}
    #}
#}


print("\n |=== preprocess16S.py (version 2.1) ===|\n")

# Check packages needed for plotting
if quality_plot:#{

    print("\nChecking packages needed for plotting...")

    try:#{
        import numpy as np
    #}
    except ImportError as imperr:#{
        print_error("'numpy' package is not installed!")
        print( str(imperr) )
        print("Please install numpy (e.g. pip3 install numpy)")
        exit(1)
    #}

    try:#{
        import matplotlib.pyplot as plt
    #}
    except ImportError as imperr:#{
        print_error("'matplotlib' package is not installed!")
        print( str(imperr) )
        print("Please install matplotlib (e.g. pip3 install matplotlib)")
        exit(1)
    #}

    print("ok\n")
#}



# Check utilities for read merging
if merge_reads:#{

    pathdirs = os.environ["PATH"].split(os.pathsep)

    for utility in ("blastn", "blastdbcmd"):#{

        utility_found = False

        for directory in pathdirs:#{
            if os.path.exists(directory) and utility in os.listdir(directory):#{
                utility_found = True
                break
            #}
        #}

        if not utility_found:#{
            print_error("{} is not installed!".format(utility))
            print("Please install {}".format(utility))
            print("""If this error still occure although you have installed everything
    -- make sure that this program is added to PATH)""")
            exit(1)
        #}
    #}

    try:#{
        import read_merging_16S
    #}
    except ImportError as imperr:#{
        print_error("Module 'read_merging_16S' not found.")
        print_error( str(imperr) )
        print("Please, run 'install_read_merging_16S.sh' provided with 'preprocess_16S.py'.")
        print("For further info see README.md in repository.")
        exit(0)
    #}
#}


print( get_work_time() + " ({}) ".format(strftime("%d.%m.%Y %H:%M:%S", localtime(start_time))) + "- Start working\n")

from re import match, findall, search
from gzip import open as open_as_gzip
from gzip import GzipFile
from _io import TextIOWrapper

# |====== Some data used by this script =====|

MAX_SHIFT = 4
# Cross-talks are searched by shifting primer sequences across 5'-end of a read,
#   starting with zero-shift, i.e.:
# RRRRRRRRRRRRRRRR
# PPPPPPP
# And for example, if MAX_SHIFT is 3, edge cases will be:
# ---RRRRRRRRRRRRRR
# PPPPPPP
#   and:
# RRRRRRRRRRRRRR
# ---PPPPPPP,
# where R is nucleotide from read, P is nucleotide from primer and dash means gap

RECOGN_PERCENTAGE = 0.52
# E.g., RECOGN_PERCENTAGE == 0.70, then if less than 70% nucleotides from primer sequence match
# corresponding nucleotides from a read, this read will be considered as cross-talk.

MATCH_DICT = {
    "A": "ANRWMDHV",
    "G": "GNRSKBDV",
    "C": "CNYSMBHV",
    "T": "TNYWKBDH",
    "N": ""
}
# This dictionary is necessary beacause primers used in 16S amplicon researches are often redundant.
# Keys in this dictionary are nucleotides that read sequence may contain.
# Values are lists (I do not mean build-in python type, obviously) of symbols a redundant primer
#   sequence may contain and match corresponding key-nucleotide.
# E.g., if we have got A from read sequence and one of ANRWMDHV symbols is located at
#    corresponding position in primer sequence, this A-nucleotide will be considered as matching.
# 'N'-nucleotide in read sequence can not match anything.

is_gzipped = lambda f: True if f.endswith(".gz") else False

# Following tuples of functions are here in order to process '.fastq' and '.fastq.gz'
#    files using the same interface.

OPEN_FUNCS = (open, open_as_gzip)

FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format .fastq line
    lambda line: line.decode("utf-8").strip()   # format .fastq.gz line
)


# |========================= Functions =========================|

def close_files(*files):#{
    """
    Function closes all files passed to it, no matter whether they are single file objects or collenctions.

    :type of the argumets:
        'dict<str: _io.TextIOWrapper' or 'gzip.GzipFile>'' or
        'list<_io.TextIOWrapper or gzip.GzipFile>' or
        'tuple<_io.TextIOWrapper or gzip.GzipFile>' or
        '_io.TextIOWrapper' or
        'gzip.GzipFile'
    """
    for obj in files:#{

        if isinstance(obj, dict) or isinstance(obj, list) or isinstance(obj, tuple):#{
            for indx_or_key in obj:#{
                obj[indx_or_key].close()
            #}
        #}

        elif isinstance(obj, TextIOWrapper) or isinstance(obj, GzipFile):#{
            obj.close()
        #}

        else:#{
            print("""If you use function 'close_files', please, store file objects in 
    lists, tuples or dictionaries or pass file objects itself to the function.""")
            print("You have passed '{}'".format( type(obj) ))
            try:#{
                obj.close()
            #}
            except:#{
                print_error("Object passed to 'close_files' function can't be closed")
                exit(1)
            #}
        #}
    #}
#}


def read_fastq_pair(read_files, fmt_func):#{
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

    if len(read_files) != 2:#{
        print("You can only pass 2 files to the function 'read_fastq_pair'!\a")
        exit(1)
    #}

    fastq_recs = dict()           # this dict should consist of two fastq-records: R1 and R2
    for key in read_files.keys():#{
        fastq_recs[key] = {                    #read all 4 lines of fastq-record
            "seq_id": fmt_func(read_files[key].readline()),
            "seq": fmt_func(read_files[key].readline()).upper(), # searching for cross-talks is case-dependent
            "opt_id": fmt_func(read_files[key].readline()),
            "qual_str": fmt_func(read_files[key].readline())
        }

        # If all lines from files are read
        if fastq_recs[key]["seq_id"] is "":
            return None
    #}
    return fastq_recs
#}


def write_fastq_record(outfile, fastq_record):#{
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
    try:#{
        outfile.write(fastq_record["seq_id"] + '\n')
        outfile.write(fastq_record["seq"] + '\n')
        outfile.write(fastq_record["opt_id"] + '\n')
        outfile.write(fastq_record["qual_str"] + '\n')
    #}
    except Exception as exc:#{
        print("\nAn error occured while writing to outfile")
        print( str(exc) )
        exit(1)
    #}
#}


def find_primer(primers, fastq_rec):#{
    """
    This function figures out, whether a primer sequence is at 5'-end of a read passed to it.
    Moreover it cuts primer sequences off if there are any.

    :param primers: list of required primer sequences
    :type primers: list<str>
    :param fastq_rec: fastq-record containing read, which this function searches for required primer in
    :type fastq_rec: dict<str: str>

    Returns tuple of following values:
    0-element: True if primer sequence is found in read, otherwise returns False.
    1-element: read (as fastq_rec-like dict) with a primer sequence cut off if there were any (and if -c is specified)
        and intact read otherwise.
    """

    read = fastq_rec["seq"]

    for primer in primers:#{
        primer_len = len(primer)

        for shift in range(0, MAX_SHIFT + 1):#{
            score_1, score_2 = 0, 0

            for pos in range(0, primer_len - shift):#{
                # if match, 1(one) will be added to score, 0(zero) otherwise
                score_1 += primer[pos+shift] in MATCH_DICT[read[pos]]
            #}
            if score_1 / primer_len >= RECOGN_PERCENTAGE:#{
                if cutoff:#{
                    cutlen = len(primer) - shift
                    fastq_rec["seq"] = read[cutlen : ]
                    fastq_rec["qual_str"] = fastq_rec["qual_str"][cutlen : ]
                    return (True, fastq_rec)
                #}
                else:#{
                    return (True, fastq_rec)
                #}
            #}

            # If there is no primer in 0-position, loop below will do unnecessary work.
            # Let it go a bit further.
            shift += 1

            for pos in range(0, primer_len - shift):#{
                score_2 += primer[pos] in MATCH_DICT[read[pos + shift]]
            #}
            if score_2 / primer_len >= RECOGN_PERCENTAGE:#{
                if cutoff:#{
                    cutlen = len(primer) - shift
                    fastq_rec["seq"] = read[cutlen : ]
                    fastq_rec["qual_str"] = fastq_rec["qual_str"][cutlen : ]
                    return (True, fastq_rec)
                #}
                else:#{
                    return (True, fastq_rec)
                #}
            #}
        #}
    #}

    return (False, fastq_rec)
#}


# This is a decorator.
# All functions that process reads do some same operations: they open read files, open result files,
#   they count how many reads are already processed and they all files mentioned above.
# So why not to use one interface during multiple proceduers?
def progress_counter(process_func, read_paths, result_paths=None, stats=None, inc_percentage=0.05):#{

    def organizer():#{

        # Collect some info
        file_type = int(match(r".*\.gz$", read_paths["R1"]) is not None)
        how_to_open = OPEN_FUNCS[file_type]
        actual_format_func = FORMATTING_FUNCS[file_type]
        readfile_length = sum(1 for line in how_to_open(read_paths["R1"], 'r'))  # lengths of read files are equal
        read_pairs_num = int(readfile_length / 4)            # divizion by 4, because there are 4 lines per one fastq-record

        # Open files
        read_files = dict()
        result_files = dict()
        try:#{
            for key in read_paths.keys():#{
                read_files[key] = how_to_open(read_paths[key])
            #}
            if result_paths is not None:#{
                for key in result_paths.keys():#{
                    result_files[key] = open(result_paths[key], 'w')
                #}
            #}
        #}
        except OSError as oserror:#{
            print("Error while opening file", str(oserror))
            close_files(read_files, result_files)
            input("Press enter to exit:")
            exit(1)
        #}

        # Proceed
        inc_percentage = 0.01
        reads_processed = 0
        spaces = 50
        next_done_percentage = inc_percentage
        printn("[" + " "*50 + "]" + " 0%")
        while reads_processed < read_pairs_num:#{

            fastq_recs = read_fastq_pair(read_files, actual_format_func)

            # Do what you need with these reads
            if result_paths is not None:#{
                process_func(fastq_recs, result_files, stats)
            #}
            else:#{
                process_func(fastq_recs)    # it is None while calulating data for plotting
            #}

            reads_processed += 1
            if reads_processed / read_pairs_num >= next_done_percentage:#{
                count = round(next_done_percentage * 100)
                spaces = 50 - int(count/2)

                printn("\r[" + "="*int(count/2) + ">" + " "*spaces + "]" + " {}% ({}/{})"
                    .format(count, reads_processed, read_pairs_num))
                next_done_percentage += inc_percentage
            #}
        #}

        print("\r[" + "="*50 + "]" + " 100% ({}/{})\n\n"
            .format(reads_processed, read_pairs_num))
        close_files(read_files, result_files)
    #}

    return organizer
#}


def find_primer_organizer(fastq_recs, result_files, stats):#{

    primer_in_R1, fastq_recs["R1"] = find_primer(primers, fastq_recs["R1"])
    primer_in_R2, fastq_recs["R2"] = find_primer(rev_primers, fastq_recs["R2"])

    if primer_in_R1 or primer_in_R2:#{
        write_fastq_record(result_files["mR1"], fastq_recs["R1"])
        write_fastq_record(result_files["mR2"], fastq_recs["R2"])
        stats["match"] += 1
    #}
    else:#{
        write_fastq_record(result_files["trR1"], fastq_recs["R1"])
        write_fastq_record(result_files["trR2"], fastq_recs["R2"])
        stats["trash"] += 1
    #}
#}


# |===== Functions for parallel quality computation =====|

# There is mo need to store these functions in memory if only one thread will be launched
if n_thr > 1:#{

    def fastq_read_packets(read_paths, reads_at_all):#{
        """
        Function-generator for retrieving FASTQ records from PE files
            and distributing them evenly between 'n_thr' processes
            for further parallel processing.

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
        how_to_open = OPEN_FUNCS[ is_gzipped(read_paths["R1"]) ]
        fmt_func = FORMATTING_FUNCS[ is_gzipped(read_paths["R1"]) ]
        read_files = dict() # dictionary for file objects

        # Open files that contain reads meant to be processed.
        for key, path in read_paths.items():#{
            read_files[key] = how_to_open(path)
        #}

        # Compute packet size (one packet -- one thread).
        pack_size = reads_at_all // n_thr
        if reads_at_all % n_thr != 0: # in order not to import 'math'
            pack_size += 1

        for i in range(n_thr):#{
            packet = list()
            for j in range(pack_size):#{
                tmp = read_fastq_pair(read_files, fmt_func)
                if tmp is None:#{ if the end of the line is reached
                    # Yield partial packet and return 'None' on next call
                    yield packet
                    # Python 2 throws a syntax error on the atttmpt of returning 'None' from the generator explicitly.
                    # But single 'return' statement returns 'None' anyway, so here it is:
                    return # None
                #}
                else:#{
                    packet.append(tmp)
                #}
            #}
            yield packet # yield full packet
        #}
    #}

    def single_qual_calcer(data, reads_at_all):#{
        """
        Function that performs task meant to be done by one process while parallel quality calculation.

        :param data: list of FASTQ-records
            (structure of these records is described in 'write_fastq_record' function);
        :type data: list< dict<str: str> >;
        :param reads_at_all: total number of read pairs in input files;
        :type reads_at_all: int;

        Returns numpy.ndarray<int> performing quality distribution of reads.
        """

        # In 2019 sequencators do not read better that Q40.
        top_x_scale, step = 40.0, 0.5
        # average read quality
        X = np.arange(0, top_x_scale + step, step)
        # amount of reads with sertain average quality
        Y = np.zeros(int(top_x_scale / step), dtype=int)

        get_phred33 = lambda symb: ord(symb) - 33

        # Processes will print number of processed reads every 'delay' reads.
        delay, i = 1000, 0

        for fastq_recs in data:#{

            for rec in fastq_recs.values():#{

                qual_str = rec["qual_str"]

                qual_array = np.array(list( map(get_phred33, qual_str) )) # get qualities
                avg_qual = round(np.mean(qual_array), 2) # get mean quality
                min_indx = ( np.abs(X - avg_qual) ).argmin() # find index in Y to increment
                Y[min_indx] += 1
            #}

            if i == delay:#{ if next 1000 reads are processed
                with count_lock:#{ synchronized incrementing
                    counter.value += delay
                #}
                
                bar_len = 50 # length of status bar
                eqs = int( bar_len*(counter.value / reads_at_all) ) # number of '=' characters
                spaces = bar_len - eqs # number of ' ' characters
                with print_lock:#{
                    printn("\r[" + "="*eqs + '>' + ' '*spaces +"] {}% ({}/{})".format(int(counter.value/reads_at_all*100),
                        counter.value, reads_at_all))
                #}
                i = 0 # reset i
            #}
            i += 1
        #}

        return Y
    #}


    def proc_init(print_lock_buff, counter_buff, count_lock_buff):#{
        """
        Function that initializes global variables that all processes shoud have access to.
        This function is meant to be passed as 'initializer' argument to 'multiprocessing.Pool' function.

        :param print_lock_buff: lock that synchronizes printing to the console;
        :type print_lock_buff: multiprocessing.Lock;
        :param counter_buff: integer number representing number of processed reads;
        :type counter_buff: multiprocessing.Value;
        :param count_lock_buff: lock that synchronizes incrementing 'counter' variable;
        :type count_lock_buff: multiprocessing.Lock;
        """

        global print_lock
        print_lock = print_lock_buff

        global counter
        counter = counter_buff

        global count_lock
        count_lock = count_lock_buff
    #}


    from functools import reduce

    def parallel_qual(read_paths, n_thr):#{
        """
        Function launches parallel quality calculations.

        :param read_paths: dict of paths to read files. it's structure is described in 'fastq_read_packets' function;
        :type read_paths: dict<str: str>;
        :param n_thr: int;
        :type n_thr: int:
        """
        
        print_lock = mp.Lock() # lock that synchronizes printing to the console;
        count_lock = mp.Lock() # lock that synchronizes 'counter' variable incrementing;
        counter = mp.Value('i', 0) # integer number representing number of processed reads;

        # In 2019 sequencators do not read better that Q40.
        top_x_scale, step = 40.0, 0.5
        # Array of average read qualities
        X = np.arange(0, top_x_scale + step, step)

        # Count number of read pairs
        reads_at_all = int( sum(1 for line in how_to_open(read_paths["R1"])) / 4 )

        # Create pool of processes
        pool = mp.Pool(n_thr, initializer=proc_init, initargs=(print_lock, counter, count_lock))
        # Run parallel calculations
        Y = pool.starmap(single_qual_calcer, [(data, reads_at_all) for data in fastq_read_packets(read_paths, reads_at_all)])
        print("\r["+"="*50+"] 100% ({}/{})\n".format(reads_at_all, reads_at_all))

        def arr_sum(Y1, Y2): # function to perform 'functools.reduce' sum
            return Y1 + Y2

        return reduce(arr_sum, Y) # element-wise sum of all processes' results
    #}
#}

# |===== END of functions for parallel quality computation =====|



# |======================= Start calculating =======================|


# |===== Prepare files for read merging =====|

# === Search for primer file in working directory if it is not specified ===

if primer_path is None:#{ # if primer file is not specified by CL argument
    # search for primer file in current directory
    looks_like_primer_file = lambda f: False if match(r".*[pP]rimers.*\.(m)?fa(sta)?(\.gz)?$", f) is None else True

    try:#{
        primer_path = ( list(filter(looks_like_primer_file, os.listdir('.'))) )[0]
    #}
    except IndexError:#{
        print("No file considered as primer file is found.")
        print("Please, specify it with '-p' option")
        exit(1)
    #}
#}


# Variable named 'file_type' below will be 0(zero) if primers are in plain FASTA file
#   and 1(one) if they are in fasta.gz file
# If primers are in .gz file, it will be opened as .gz file and as standard text file otherwise.
file_type = int(match(r".*\.gz$", primer_path) is not None)
how_to_open = OPEN_FUNCS[file_type]
actual_format_func = FORMATTING_FUNCS[file_type]


# === Retrieve primers sequences from file ===

# Assuming that each primer sequence is written in one line 
primers = list()     # list of required primer sequences
primer_ids = list()  # list of their names
primer_file = None
try:#{
    primer_file = how_to_open(primer_path, 'r')
    for i, line in enumerate(primer_file):#{
        if i % 2 == 0:#{                                            # if line is sequense id
            primer_ids.append(actual_format_func(line))
        #}
        else:#{                                                     # if line is a sequence
            line = actual_format_func(line).upper()
            err_set = set(findall(r"[^ATGCRYSWKMBDHVN]", line))     # primer validation
            if len(err_set) != 0:#{
                print_error("There are some inappropriate symbols in your primers!")
                print("Here they are:")
                for err_symb in err_set:#{
                    print(" " + err_symb)
                #}
                print("See you later with correct data")
                input("Press enter to exit:")
                exit(1)
            #}
            primers.append(line)
        #}
    #}
#}
except OSError as oserror:#{
    print("Error while reading primer file.\n", str(oserror))
    input("Press enter to exit:")
    exit(1)
#}
finally:#{
    if primer_file is not None:
        primer_file.close()
#}

# The last line in fasta file can be a blank line.
# We do not need it.
try:#{
    primer_ids.remove('')
#}
except ValueError:#{
    pass
#}


# Reversed verfion of this list is used to speed up cross-talks searching:
#    when reverse read is processing, this reversed lisr is passed to "find_primer" funtion
#    since finding reverse primer in reverse read is more likely than finding a forward one.
rev_primers = reversed(primers)


# === Select read files if they are not specified ===

# If read files are not specified by CL arguments
if len(read_paths) == 0:#{
    # Search for read files in current directory.
    looks_like_forw_reads = lambda f: False if match(r".*R1.*\.f(ast)?q(\.gz)?$", smth) is None else True
    forw_reads = ( list(filter(looks_like_forw_reads, os.listdir('.'))) )[0]
    rev_reads = forw_reads.replace("R1", "R2")
    if not os.path.exists(rev_reads):#{
        print_error("File with reverse reads not found in current directory!")
        print("File with forward ones: '{}'".format(forw_reads))
        exit(1)
    #}
    read_paths = {"R1": forw_reads, "R2": rev_reads}
#}


# I need to keep names of read files in memory in order to name result files properly.
names = dict()
for key, path in read_paths.items():#{
    file_name_itself = os.path.basename(read_paths[key])
    names[key] = search(r"(.*)\.f(ast)?q", file_name_itself).group(1)
#}


# === Create output directory. ===
if not os.path.exists(outdir_path):#{
    try:#{
        os.mkdir(outdir_path)
    #}
    except OSError as oserror:#{
        print_error("Error while creating result directory")
        print( str(oserror) )
        close_files(read_files)
        exit(1)
    #}
#}



# === Create and open result files. ===

files_to_gzip = list()
result_files = dict()
# Keys description:
# 'm' -- matched (i.e. sequence with primer in it);
# 'tr' -- trash (i.e. sequence without primer in it);
# I need to keep these paths in memory in order to gzip corresponding files afterwards.
result_paths = {
    # We need trash anyway (trash without primers and, therefore, without 16S data):
    "mR1": "{}{}{}.16S.fastq".format(outdir_path, os.sep, names["R1"]),
    "mR2": "{}{}{}.16S.fastq".format(outdir_path, os.sep, names["R2"]),
    "trR1": "{}{}{}.trash.fastq".format(outdir_path, os.sep, names["R1"]),
    "trR2": "{}{}{}.trash.fastq".format(outdir_path, os.sep, names["R2"])
}

for path in result_paths.values():#{
    rm_content = open(path, 'w') # remove data that might be in these files -- it will be rewritten
    rm_content.close()
#}

files_to_gzip.extend(result_paths.values())

primer_stats = {
    "match": 0,           # number of read pairs with primers
    "trash": 0           # number of cross-talks
}

print("\nFollowing files will be processed:")
for i, path in enumerate(read_paths.values()):#{
    print("  {}. '{}'".format(i+1, os.path.abspath(path)))
#}
print("Number of theads: {}".format(n_thr))
if cutoff:
    print("Primer sequences will be cut off.")
if merge_reads:
    print("Reads will be merged together.")
if quality_plot:
    print("Quality graph will be plotted.")
print('-'*20+'\n')


# |===== Start the process of searching for cross-talks =====|

print("\n{} - Searching for cross-talks started".format(get_work_time()))
print("Proceeding...\n")

primer_task = progress_counter(find_primer_organizer, read_paths, result_paths, primer_stats)
primer_task()

print("{} - Searching for cross-talks is completed".format(get_work_time()))
print("""\n{} read pairs with primer sequences are found.
{} read pairs considered as cross-talks are found.""".format(primer_stats["match"], primer_stats["trash"]))
print('\n' + '~' * 50 + '\n')

# |===== The process of searching for cross-talks is completed =====|



# |===== Start the process of read merging =====|

if merge_reads:#{

    merge_result_files = read_merging_16S.merge_reads(result_paths["mR1"], result_paths["mR2"], 
        outdir_path=outdir_path, n_thr=n_thr)
    merging_stats = read_merging_16S.get_merging_stats()

    files_to_gzip.extend(merge_result_files.values())
#}


# |===== The process of merging reads is completed =====|


# |===== Prepare data for plotting =====|

if quality_plot:#{

    from math import log # for log scale in plot

    top_x_scale, step = 40.0, 0.5
    # average read quality
    X = np.arange(0, top_x_scale + step, step)
    # amount of reads with sertain average quality
    Y = np.zeros(int(top_x_scale / step), dtype=int)

    get_phred33 = lambda symb: ord(symb) - 33

    # This function will be used as "organizer"
    def add_data_for_qual_plot(fastq_recs):#{

        # quality_strings = list()
        for rec in fastq_recs.values():#{

            qual_str = rec["qual_str"]

            qual_array = np.array(list( map(get_phred33, qual_str) ))
            avg_qual = round(np.mean(qual_array), 2)
            min_indx = ( np.abs(X - avg_qual) ).argmin()
            Y[min_indx] += 1
            #}
        #}
    #}

    if merge_reads:#{
        data_plotting_paths = {
            "R1": merge_result_files["merg"]
        }
    #}
    else:#{
        data_plotting_paths = {
            "R1": result_paths["mR1"],
            "R2": result_paths["mR2"]
        }
    #}

    print("\n{} - Calculations for plotting started".format(get_work_time()))
    print("Proceeding...\n")


    if n_thr == 1:#{
        plotting_task = progress_counter(add_data_for_qual_plot, data_plotting_paths)
        plotting_task()
    #}
    else:#{
        Y = parallel_qual(data_plotting_paths, n_thr)
    #}

    print("\n{} - Calculations for plotting are completed".format(get_work_time()))
    print('\n' + '~' * 50 + '\n')


 # |===== Plot a graph =====|

    image_path = os.sep.join((outdir_path, "quality_plot.png"))

    fig, ax = plt.subplots()

    Y = np.array(list(map( lambda y: log(y) if y > 0 else 0 , Y)))

    ax.plot(X[:len(Y)], Y, 'r')
    ax.set(title="Quality per read distribution",
        xlabel="avg quality (Phred33)",
        ylabel="ln(amount of reads)")
    ax.grid()
    fig.savefig(image_path)

    # Open image as image via xdg-open in order not to pause executing main script while
    #   you are watching at pretty plot.

    img_open_util = "xdg-open"
    util_found = False
    for directory in os.environ["PATH"].split(os.pathsep):#{
        if os.path.isdir(directory) and img_open_util in os.listdir(directory):#{
            util_found = True
            break
        #}
    #}

    if util_found:
        os.system("{} {}".format(img_open_util, image_path))
    print("Plot image is here:\n  '{}'\n\n".format(image_path))
#}


# |===== The process of plotting is completed =====|


# Remove empty files
for file in files_to_gzip:#{
    if os.stat(file).st_size == 0:#{
        os.remove(file)
        print("'{}' is removed since it is empty".format(file))
    #}
#}


# Gzip result files
gzip_util = "gzip"
util_found = False
for directory in os.environ["PATH"].split(os.pathsep):#{
    if os.path.isdir(directory) and gzip_util in os.listdir(directory):#{
        util_found = True
        break
    #}
#}

print("\n{} - Gzipping result files...".format(get_work_time()))
for file in files_to_gzip:#{
    if os.path.exists(file):#{
        if util_found:#{
            os.system("{} {}".format(gzip_util, file))
        #}
        else:#{
            with open(file, 'r') as plain_file, open_as_gzip(file+".gz", 'wb') as gz_file:#{
                for line in plain_file:#{
                    gz_file.write(bytes(line, "utf-8"))
                #}
            #}
            os.unlink(file)
        #}
        print("\'{}\' is gzipped".format(file))
    #}
#}
print("{} - Gzipping is completed\n".format(get_work_time()))
print("Result files are placed in the following directory:\n  '{}'\n\n".format(os.path.abspath(outdir_path)))


# Create log file
start_time_fmt = strftime("%d_%m_%Y_%H_%M_%S", localtime(start_time))
with open("{}{}preprocess16S_{}.log".format(outdir_path, os.sep, start_time_fmt).replace(" ", "_"), 'w') as logfile:#{

    logfile.write("Script 'preprocess_16S.py' was ran at {}\n".format(start_time_fmt))
    end_time = strftime("%d_%m_%Y_%H_%M_%S", localtime(time()))
    logfile.write("Script completed it's job at {}\n\n".format(end_time))
    logfile.write("The man who ran it was searching for following primer sequences:\n\n")

    for i in range(len(primers)):#{
        logfile.write("{}\n".format(primer_ids[i]))
        logfile.write("{}\n".format(primers[i]))
    #}

    logfile.write("""\n{} read pairs with primer sequences have been found.
{} read pairs without primer sequences have been found.\n""".format(primer_stats["match"], primer_stats["trash"]))

    if cutoff:#{
        logfile.write("These primers were cut off.\n")
    #}

    if merge_reads:#{
        logfile.write("\n\tReads were merged\n\n")
        logfile.write("{} read pairs have been merged.\n".format(merging_stats[0]))
        logfile.write("{} read pairs haven't been merged.\n".format(merging_stats[1]))
        logfile.write("{} read pairs have been considered as too short for merging.\n".format(merging_stats[2]))
    #}

    logfile.write("\nResults are in the following directory:\n  '{}'\n".format(outdir_path))

    if quality_plot:#{
        logfile.write("\nQuality graph was plotted. Here it is: \n  '{}'\n".format(image_path))
    #}
#}


print('\n'+get_work_time() + " ({}) ".format(strftime("%d.%m.%Y %H:%M:%S", localtime(time()))) + "- Job is successfully completed!\n")
exit(0)

