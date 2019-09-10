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
    # Neither does 'input' in python2.
    raw_input("Press ENTER to exit:")
    exit(1)
#}

def print_error(text):#{
    "Function for printing error messages"
    print("\n   \a!! - ERROR: " + text + '\n')
#}

from datetime import datetime
start_time = datetime.now().strftime("%Y-%m-%d %H.%M.%S")


# |===== Handle CL arguments =====|

import os
import getopt
from sys import argv
import multiprocessing as mp

usage_msg = """
DESCRIPTION:\n
preprocess16S.py -- this script is designed for preprocessing reads from 16S regions of rDNA.\n
It works with Illumina pair-end reads. More precisely, it detects and removes reads,
    that came from other sample (aka cross-talks), relying on the information,
    whether there are PCR primer sequences in these reads or not. If required primer sequence is
    found in a read, therefore, it is a read from 16S rDNA and we need it.\n
Moreover, it can cut these primers off and merge reads by using 'read_merging_16S' module.\n
Reads should be stored in files, both of which are of FASTQ format or both are gzipped (i.e. fastq.gz).
Sequences of required primers are retrieved from FASTA file (or fasta.gz).\n
Attention! This script cannot be executed by python interpreter version < 3.0!
----------------------------------------------------------\n
OPTIONS:\n
    -c or --cutoff 
        cut primer sequences off;
    -m or --merge-reads
        if this option is specified, script will merge reads together;
    -q or --quality-plot
        plot a graph of read quality distribution;
    -p or --primers
        FASTA file, in which primer sequences are stored;
    -1 or --R1
        FASTQ file, in which forward reads are stored;
    -2 or --R2
        FASTQ file, in which reverse reads are stored;
    -o or --outdir
        directory, in which result files will be placed.
----------------------------------------------------------\n
EXAMPLES:\n
  1) Remove cross-talks in files 'forw_R1_reads.fastq.gz' and 'rev_R2_reads.fastq.gz' depending on
     primer sequenes in file 'V3V4_primers.fasta'. Cut primer sequences off and put result files in
     the directory named 'outdir'.\n
       ./preprocess16S.py -p V3V4_primers.fasta -1 forw_R1_reads.fastq.gz -2 rev_R2_reads.fastq.gz -o outdir -c
"""

try:#{
    opts, args = getopt.getopt(argv[1:], "hcmqp:1:2:o:t:", 
        ["help", "cutoff", "merge-reads", "quality-plot", "primers=", "R1=", "R2=", "outdir=",
        "V3-V4", "threads="])
#}
except getopt.GetoptError as opt_err:#{
    print(opt_err)
    print(usage_msg)
    exit(2)
#}

def check_file_existance(path):#{
    if os.path.exists(path):#{
        return
    #}
    else:#{
        print("\nFile '{}' does not exist!".format(path))
        exit(1)
    #}
#}

cutoff = False
merge_reads = False
V3V4 = False
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

    # elif opt == "--V3-V4":#{
    #     V3V4 = True
    # #}

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
            print("ATTENTION!\n\tYou should specify both forward and reverse reads!")
            exit(1)
        check_file_existance(arg)
        read_paths["R1"] = arg
    #{

    elif opt in ("-2", "--R2"):#{
        if not "-1" in argv and not "--R1" in argv:
            print("\nATTENTION!\n\tYou should specify both forward and reverse reads!")
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
            print(" Warning! - {} threads has been specified, while there are {} cores in your CPU.\n".format(n_thr,
                mp.cpu_count()))
        #}
    #}
#}


# # --V3-V4 option can't be specified without --merge-reads (-m)
# if V3V4 and not merge_reads:#{
#     print("\nError!\a")
#     print("'--V3-V4' option can't be specified without '--merge-reads' (-m) option")
#     print("""The reason for this is that checking for constant region 
#     between V3 and V4 variable regions is performed only while merging reads and 
#     only if your target sequences contain V3 and V4 regions 
#     and a constant region between them.\n""")
#     exit(1)
# #}

print("\n |=== preprocess16S.py (version 2.1) ===|\n")

# Check packages needed for plotting
if quality_plot:#{

    print("\nChecking packages needed for plotting...")
    imp_error_msg = """\tAttention!\nPKG package is not installed. Graph will not be plotted
If you want to use this feature, please install PKG (e.g. pip install PKG)"""

    try:#{
        import numpy as np
    #}
    except ImportError as imperr:#{
        print("\nError: {}".format(str(imperr)))
        print(imp_error_msg.replace("PKG", "numpy"))
        quality_plot = False
    #}

    try:#{
        import matplotlib.pyplot as plt
    #}
    except ImportError as imperr:#{
        print("\nError: {}".format(str(imperr)))
        print(imp_error_msg.replace("PKG", "matplotlib"))
        quality_plot = False
    #}

    if quality_plot:
        print("ok")
        
    del imp_error_msg
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
            print("\tAttention!\n{} is not found in your system. Reads will not be merged.".format(utility))
            print("If you want to use this feature, please install {}".format(utility))
            print("""If this error still occure although you have installed everything 
    -- make sure that this program is added to PATH)""")
            merge_reads = False
        #}
    #}

    try:#{
        import read_merging_16S
    #}
    except ImportError as imperr:#{
        print("\nError: {}".format(str(imperr)))
        print("\tModule 'read_merging_16S' not found. Reads will not be merged.")
        print("To use this feature, you need to install 'read_merging_16S'.")
        print("More precisely, you need to run 'install.sh' provided with 'preprocess_16S.py'.")
        print("For further info see README.md in repo.")
        merge_reads = False
    #}
#}


from re import match
from re import findall
from gzip import open as open_as_gzip
from gzip import GzipFile
from _io import TextIOWrapper
from array import array

# |====== Some data used by this script =====|

MAX_SHIFT = 4
# Primer sequences are searched with respect to shift.
# E.g., if MAX_SHIFT == 3, primer sequence will be searched from this collocation:
# ---RRRRRRRRRRRRRR
# PPPPPPP
#   to this collocation:
# RRRRRRRRRRRRRR
# ---PPPPPPP,
# where R is nucleotide from read, P is nucleotide from primer and dash means gap

RECOGN_PERCENTAGE = 0.52
# E.g., RECOGN_PERCENTAGE == 0.7, then if 70% or more nucleotides from primer sequence match
# corresponding nucleotides form read, this read will be considered as read with primer sequence in it
# and should be written to appropriate file.

MATCH_DICT = {
    "A": "ANRWMDHV",
    "G": "GNRSKBDV",
    "C": "CNYSMBHV",
    "T": "TNYWKBDH",
    "N": ""
}
# Key is nucleotide that read sequence may contain.
# Value is a list (i do not mean biuld-in python type, obviously) of symbols primer
# sequence may contain and match corresponding key-nucleotide.
# E.g., if we have got A from read sequence and one of ANRWMDHV symbols is located at
# corresponding position in primer sequence, this A-nucleotide will be considered as matched.
# N from read can not match any.
# Info is got from this table: https://www.bioinformatics.org/sms/iupac.html

OPEN_FUNCS = (open, open_as_gzip)

FORMATTING_FUNCS = (
    lambda line: line.strip().upper(),   # format .fastq line
    lambda line: line.decode("utf-8").strip().upper()   # format .fastq.gz line 
)
# Data from .fastq and .fastq.gz should be parsed in different way,
#   because data from .gz is read as bytes, not str.


# |========================= Functions =========================|

def close_files(*files):#{
    """
    Function closes all files passed to it, no matter whether they are single file objects or collenctions.

    :type of the argumets:
        dict<str: _io.TextIOWrapper or gzip.GzipFile>
        or list<_io.TextIOWrapper or gzip.GzipFile>
        or tuple<_io.TextIOWrapper or gzip.GzipFile>
        or _io.TextIOWrapper
        or gzip.GzipFile
    :return: void
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
            try:#{
                obj.close()
            #}
            except:#{
                print("Object passed to 'close_files' function can't be closed")
                exit(1)
            #}
        #}
    #}
#}


def write_fastq_record(outfile, fastq_record):#{
    """
    :param outfile: file, which data from fastq_record is written in
    :type outfile: _io.TextIOWrapper
    :param fastq_record: dict of 4 elements. Elements are four corresponding lines of fastq
    :type fastq_record: dict<str: str>
    :return: void
    """
    try:#{
        outfile.write(fastq_record["seq_id"] + '\n')
        outfile.write(fastq_record["seq"] + '\n')
        outfile.write(fastq_record["optional_id"] + '\n')
        outfile.write(fastq_record["quality_str"] + '\n')
    #}
    except Exception as exc:#{
        print("\nAn error occured while writing to outfile")
        print(str(exc))
        exit(1)
    #}
#}


def read_fastq_record(read_files, fmt_func):#{

    if len(read_files) != 1 and len(read_files) != 2:#{
        print("You can pass only 1 or 2 files to the function 'read_pair_of_reads'!\a")
        exit(1)
    #}

    fastq_recs = dict()           # this dict should consist of two fastq-records: from R1 and from R2
    for key in read_files.keys():#{
        fastq_recs[key] = {                    #read all 4 lines of fastq-record
            "seq_id": fmt_func(read_files[key].readline()),
            "seq": fmt_func(read_files[key].readline()),
            "optional_id": fmt_func(read_files[key].readline()),
            "quality_str": fmt_func(read_files[key].readline())
        }

        if fastq_recs[key]["seq_id"] is "":
            return None
    #}
    return fastq_recs
#}


def find_primer(primers, fastq_rec):#{
    """
    This function figures out, whether a primer sequence is in read passed to it.
    Pluralistically it cuts primer sequences off if there are any.

    :param primers: list of required primer sequences
    :type primers: list<str>
    :param fastq_rec: fastq-record containing read, which this function searches for required primer in
    :type fastq_rec: dict<str: str>
    :return: if primer sequence is found in read -- returns True, else returns False
    Moreover, returns a read with a primer sequence cut off if there were any (and if -c is specified)
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
                    fastq_rec["quality_str"] = fastq_rec["quality_str"][cutlen : ]
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
                    fastq_rec["quality_str"] = fastq_rec["quality_str"][cutlen : ]
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


def select_file_manually(message):#{
    """
    Function provides you with an interface of selecting a file by typing path to it to the command line.

    :param message: message shown on the screen offering to enter the path
    :type message: str
    :return: path to manually selected file
    :return type: str
    """
    print('\n' + '~' * 50 + '\n')

    while True:#{
        path = input(message)
        if path == 'q!':#{
            exit(0)
        #}
        if os.path.exists(path):#{
            print("\tok...")
            return path
        #}
        else:#{
            print("\tERROR\tThere is no file named \'{}\'\a".format(path))
        #}
    #}
#}



# This is a decorator.
# All functions that process reads do some same operations: they open read files, open result files, 
#   they count how many reads are already processed and they close all files mentioned above.
# This is why I use decorator here, although in some cases I didn't find the way how to implement it smartly enough.
# But the temptation was too strong to resist ;)
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
        print("Proceeding...\n")
        print("[" + " "*50 + "]" + " 0%", end="")
        while reads_processed < read_pairs_num:#{

            fastq_recs = read_fastq_record(read_files, actual_format_func)

            # for file_key in fastq_recs.keys():#{
            #     for field_key in fastq_recs[file_key].keys():#{
            #         fastq_recs[file_key][field_key] = actual_format_func(fastq_recs[file_key][field_key])
            #     #}
            # #}

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

                print("\r[" + "="*int(count/2) + ">" + " "*spaces + "]" + " {}% ({}/{})"
                    .format(count, reads_processed, read_pairs_num), end="")
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



# |======================= Start calculating =======================|


# |===== Prepare files for read merging =====|

# === Select primer file if it is not specified ===

if primer_path is None:#{ # if primer file is not specified by CL argument
    # search for primer file in current directory
    print('\n' + '~' * 50 + '\n')
    for smth in os.listdir('.'):#{
        if match(r".*[pP]rimers.*\.(m)?fa(sta)?$", smth) is not None:#{
            primer_path = smth
            break
        #}
    #}
    reply = None
    if primer_path is not None:#{
        with open(primer_path, 'r') as primer_file:#{
            file_content = primer_file.read()
        #}
        message = """File named '{}' is found.\n Here are primers stored in this file: 
    \n{} \nFile \'{}\' will be used as primer file.
    Do you want to select another file manually instead of it?
        Enter \'y\' to select other file,
        \'q!\' -- to exit
        or anything else (e.g. merely press ENTER) to continue:""".format(primer_path, file_content, primer_path)
    #}
    else:#{
        message = """No file considered as primer file is found.
    Do you want to select file manually?
    Enter \'y\' to select other file or anything else (e.g. merely press ENTER) to exit:"""
    #}
    reply = input(message)
    if reply == "q!":#{
        print("Exiting...")
        exit(0)
    #}
    if reply == 'y':#{
        message = "Enter path to .fa primer file (or \'q!\' to exit):"
        primer_path = select_file_manually(message)
    #}
    else:#{
        if primer_path is not None:#{
            pass
        #}
        else:#{
            exit(0)
        #}
    #}
    print('\n' + '~' * 50 + '\n')
#}


# Variable named 'file_type' below will be 0(zero) if primers are in .fa file
#   and 1(one) if they are in fa.gz file
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
            line = actual_format_func(line)       
            err_set = set(findall(r"[^ATGCRYSWKMBDHVN]", line))     # primer validation
            if len(err_set) != 0:#{
                print("! There are some inappropriate symbols in your primers. Here they are:")
                for err_symb in err_set:#{
                    print(err_symb)
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


rev_primers = reversed(primers)


# === Select read files if they are not specified ===

# If read files are not specified by CL arguments
if len(read_paths) == 0:#{
    # Search for read files in current directory.
    read_paths = dict()
    for smth in os.listdir('.'):#{
        if match(r".*_R1_.*\.fastq(\.gz)?$", smth) is not None and os.path.exists(smth.replace("_R1_", "_R2_")):#{
            read_paths["R1"] = smth
            read_paths["R2"] = smth.replace("_R1_", "_R2_")
            break
        #}
    #}
    reply = None
    if len(read_paths) == 0:#{
        message = """\nNo files considered as read files are found. 
        Do you want to select other files manually?
        Enter \'y\' to select other files or anything else (e.g. merely press ENTER) to exit:"""
    #}
    else:#{
        message = """Files named\n\t'{}',\n\t'{}'\n are found. 
        They will be used as read files. 
    Do you want to select other files manually instead of them?
        Enter \'y\' to select other files,
        \'q!\' -- to exit
        or anything else (e.g. merely press ENTER) to continue:""".format(read_paths["R1"], read_paths["R2"])
    #}
    reply = input(message)
    if reply == "q!":#{
        print("Exiting...")
        exit(0)
    #}
    if reply == "y":#{
        read_paths = dict()
        for R_12, forw_rev in zip(("R1", "R2"), ("forward", "reverse")):#{
            message = "Enter path to the .fastq file with {} reads (or \'q!\' to exit):".format(forw_rev)
            read_paths[R_12] = select_file_manually(message)
        #}
        # check if both of read files are of fastq format or both are gzipped
        check = 0
        for i in read_paths.keys():#{
            check += int(match(r".*\.gz$", read_paths[i]) is not None)
        #}
        if check == 1:#{
            print("""\n\tATTENTION!
    \tBoth of read files should be of fastq format or both be gzipped (i.e. fastq.gz).
    \tPlease, make sure this requirement is satisfied.""")
            input("Press enter to exit:")
            exit(1)
        #}
    #}
    else:#{
        if len(read_paths) == 2:#{
            pass
        #}
        else:#{
            exit(0)
        #}
    #}
    print('\n' + '~' * 50 + '\n')
#}

# I need to keep names of read files in memory in order to name result files properly.
names = dict()
# file_name_itself = read_paths["R1"][read_paths["R1"].rfind(os.sep)+1 :] # get rid of, e.g. '/../' in path
file_name_itself = os.path.basename(read_paths["R1"])
names["R1"] = file_name_itself[:file_name_itself.rfind(".fastq")]
# names["R1"] = match(r"(.*)\.f(ast)?q(\.gz)?$", file_name_itself).groups(0)[0]
file_name_itself = os.path.basename(read_paths["R2"])
# names["R2"] = match(r"(.*)\.f(ast)?q(\.gz)?$", file_name_itself).groups(0)[0]
names["R2"] = file_name_itself[:file_name_itself.rfind(".fastq")]


# === Create output directory. ===
if not os.path.exists(outdir_path):#{
    try:#{
        os.mkdir(outdir_path)
    #}
    except OSError as oserror:#{
        print("Error while creating result directory\n", str(oserror))
        close_files(read_files)
        input("Press enter to exit:")
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

files_to_gzip.extend(result_paths.values())

primer_stats = {
    "match": 0,           # number of read pairs with primers
    "trash": 0           # number of cross-talks
}

print("\nFollowing files will be processed:")
for i, path in enumerate(read_paths.values()):#{
    print("  {}. '{}'".format(i+1, os.path.abspath(path)))
#}
print() # just print blank line


is_gzipped = lambda f: True if f.endswith(".gz") else False


# |===== Start the process of searching for cross-talks =====|

print("\nSearching for cross-talks started")
primer_task = progress_counter(find_primer_organizer, read_paths, result_paths, primer_stats)

primer_task()


# parallel_uncrossing(primers, read_paths, result_paths, n_thr)


print("\nSearching for cross-talks is completed")
print("""\n{} read pairs with primer sequences are found.
{} read pairs considered as cross-talks are found.""".format(primer_stats["match"], primer_stats["trash"]))
print('\n' + '~' * 50 + '\n')

# |===== The process of searching for cross-talks is completed =====|



# |===== Start the process of read merging =====|

if merge_reads:#{

    merge_result_files = read_merging_16S.merge_reads(result_paths["mR1"], result_paths["mR2"], 
        outdir_path=outdir_path, V3V4=V3V4)
    merging_stats = read_merging_16S.get_merging_stats()

    files_to_gzip.extend(merge_result_files.values())
#}


# |===== The process of merging reads is completed =====|


# /=/=/=/=//=/=/=/=//=/=/=/=//=/=/=/=//=/=/=/=//=/=/=/=//=/=/=/=/


def fastq_read_packets(read_paths, reads_at_all):#{
    how_to_open = OPEN_FUNCS[ is_gzipped(read_paths["R1"]) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(read_paths["R1"]) ]
    read_files = dict()
    for key, path in read_paths.items():#{
        read_files[key] = how_to_open(path)
    #}

    pack_size = reads_at_all // n_thr
    if reads_at_all % n_thr != 0:
        pack_size += 1

    iters = reads_at_all // pack_size
    if reads_at_all % pack_size != 0:
        iters += 1


    for i in range(iters):#{
        packet = list()
        for j in range(pack_size):#{
            tmp = read_fastq_record(read_files, fmt_func)
            if tmp is None:
                yield packet
                return None
            else:
                packet.append(tmp)
        #}
        yield packet
    #}
#}

def single_qual_calcer(data, reads_at_all):#{

    # print("Process starts")

    top_x_scale, step = 40.0, 0.5
    # average read quality
    X = np.arange(0, top_x_scale + step, step)
    # amount of reads with sertain average quality
    Y = np.zeros(int(top_x_scale / step), dtype=int)

    get_phred33 = lambda symb: ord(symb) - 33

    # print(' -- data len = ' + str(len(data)))

    delay, i = 1000, 0

    for fastq_recs in data:#{

        for rec in fastq_recs.values():#{

            qual_str = rec["quality_str"]

            qual_array = np.array(list( map(get_phred33, qual_str) ))
            avg_qual = round(np.mean(qual_array), 2)
            min_indx = ( np.abs(X - avg_qual) ).argmin()
            Y[min_indx] += 1
        #}

        if i == delay:#{
            with count_lock:#{
                counter.value += delay
            #}
            
            bar_len = 50
            eqs = int( bar_len*(counter.value / reads_at_all) )
            spaces = bar_len - eqs
            arr = '>' if eqs < bar_len else ''
            with print_lock:#{
                print("\r[" + "="*eqs + arr + ' '*spaces +"] {}% ({}/{})".format(int(counter.value/reads_at_all*100),
                    counter.value, reads_at_all), end="")
            #}
            i = 0
        #}
        i += 1
    #}

    return Y
#}


def proc_init(print_lock_buff, counter_buff, count_lock_buff):#{

    global print_lock
    print_lock = print_lock_buff

    global counter
    counter = counter_buff

    global count_lock
    count_lock = count_lock_buff
#}

from functools import reduce
def parallel_qual(read_paths, n_thr):#{

    print("\nCalculations for plotting started")
    
    print_lock = mp.Lock()
    count_lock = mp.Lock()
    counter = mp.Value('i', 0)

    top_x_scale, step = 40.0, 0.5
    # average read quality
    X = np.arange(0, top_x_scale + step, step)

    print("Proceeding...")

    reads_at_all = int( sum(1 for line in how_to_open(read_paths["R1"])) / 4 )

    pool = mp.Pool(n_thr, initializer=proc_init, initargs=(print_lock, counter, count_lock))
    Y = pool.starmap(single_qual_calcer, [(data, reads_at_all) for data in fastq_read_packets(read_paths, reads_at_all)])
    print("\r["+"="*50+"] 100% ({}/{})\n".format(reads_at_all, reads_at_all))

    print("\nCalculations for plotting are completed")
    print('\n' + '~' * 50 + '\n')

    def arr_sum(Y1, Y2):
        return Y1 + Y2

    return reduce(arr_sum, Y)
#}


# /=/=/=/=//=/=/=/=//=/=/=/=//=/=/=/=//=/=/=/=//=/=/=/=//=/=/=/=/



# |===== Prepare data for plotting =====|

if quality_plot:#{

    from math import log # for log scale in plot

    img_open_util = "xdg-open"
    util_found = False
    for directory in os.environ["PATH"].split(os.pathsep):#{
        if os.path.isdir(directory) and img_open_util in os.listdir(directory):#{
            util_found = True
            break
        #}
    #}

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

            qual_str = rec["quality_str"]

            qual_array = np.array(list( map(get_phred33, qual_str) ))
            avg_qual = round(np.mean(qual_array), 2)
            min_indx = ( np.abs(X - avg_qual) ).argmin()
            Y[min_indx] += 1
            #}
        #}
    #}

    # Name without "__R1__" and "__R2__":
    # more_common_name = os.path.basename(read_paths["R1"])[: os.path.basename(read_paths["R1"]).find("_R1_")]

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

    if n_thr == 1:#{
        print("Single")
        print("\nCalculations for plotting started")
        plotting_task = progress_counter(add_data_for_qual_plot, data_plotting_paths)
        plotting_task()
        print("\nCalculations for plotting are completed")
        print('\n' + '~' * 50 + '\n')
    #}
    else:#{
        print("Multi")
        Y = parallel_qual(data_plotting_paths, n_thr)
    #}


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

# print("\nGzipping result files...")
# for file in files_to_gzip:#{
#     if os.path.exists(file):#{
#         if util_found:#{
#             os.system("{} {}".format(gzip_util, file))
#         #}
#         else:#{
#             with open(file, 'r') as plain_file, open_as_gzip(file+".gz", 'wb') as gz_file:#{
#                 for line in plain_file:#{
#                     gz_file.write(bytes(line, "utf-8"))
#                 #}
#             #}
#             os.unlink(file)
#         #}
#         print("\'{}\' is gzipped".format(file))
#     #}
# #}
# print("Gzipping is completed\n")
# print("Result files are placed in the following directory:\n\t'{}'\n\n".format(os.path.abspath(outdir_path)))


# Create log file
with open("{}{}preprocess16S_{}.log".format(outdir_path, os.sep, start_time).replace(" ", "_"), 'w') as logfile:#{

    logfile.write("Script 'preprocess_16S.py' was ran at {}\n".format(start_time.replace('.', ':')))
    end_time = datetime.now().strftime("%Y-%m-%d %H.%M.%S")
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

exit(0)

