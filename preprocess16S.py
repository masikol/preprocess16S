#!/usr/bin/env python3
"""This script preprocesses reads from 16S regions of rDNA. It works with Illumina pair-end reads.
More precisely, it detects and removes reads, that came from other sample (aka cross-talks), relying on the information,
    whether there are PCR primer sequences in these reads or not. If required primer sequence is
    found in a read, therefore, it is a read from 16S rDNA and we need it.
Moreover, it can cut these primers off and merge reads by using 'read_merging_16S' module.
Reads should be stored in files, both of which are of fastq format or both are gzipped (i.e. fastq.gz).
Sequences of required primers are retrieved from .fa or fa.gz file.

Attention! This script cannot be executed by python interpreter version < 3.0!

Usage:
    python preprocess16S.py [-p primer_file] [-1 forward_reads -2 reverse_reads] [-o output_dir]
Options:
    -c or --cutoff 
        cut primer sequences off;
    -m or --merge-reads
        all of a sudden, if this option is specified, script will merge reads together
    --V3-V4
        by specifying this option, more accurate read merging can be performed,
        but only if '--merge-reads' (-m) option is specified and 
        target sequences contain V3 and V4 regions and a constant region between them.
    -q or --quality-plot
        plot a graph of read quality distribution
    -p or --primers
        file, in which primer sequences are stored;
    -1 or --R1
        file, in which forward reads are stored;
    -2 or --R2
        file, in which reverse reads are stored;
    -o or --outdir
        directory, in which result files will be placed.

Last modified 14.06.2019
"""

# |===== Colors =====|

def get_colored_print(color):
    def colored_print(text):
        print(color + text + "\033[0m")
    return colored_print

YELLOW = "\033[1;33m"
RED = "\033[1;31m"
GREEN = "\033[1;32m"

print_red = get_colored_print(RED)
print_yellow = get_colored_print(YELLOW)
print_green = get_colored_print(GREEN)


# |===== Check python interpreter version. =====|

from sys import version_info
if (version_info.major + 0.1 * version_info.minor) < (3.0 - 1e-6):        # just-in-case 1e-6 substraction
    print_red("\t\nATTENTION!\nThis script cannot be executed by python interpreter version < 3.0!")
    print_red("\tYour python version: {}.{}".format(version_info.major, version_info.minor))
    print_red("Try 'python3 preprocess16S.py' or update your interpreter.")
    exit(0)


from datetime import datetime
now = datetime.now().strftime("%Y-%m-%d %H.%M.%S")


# |===== Handle CL arguments =====|

import os
import getopt
from sys import argv

usage_msg = """usage:
    python preprocess16S.py -p <primer_file.mfa> -1 <forward_reads.fastq> -2 <reverse_reads.fastq> [-o output_dir]"""

try:
    opts, args = getopt.getopt(argv[1:], "hcmqp:1:2:o:", 
        ["help", "cutoff", "merge-reads", "V3-V4", "quality-plot", "primers=", "R1=", "R2=", "outdir="])
except getopt.GetoptError as opt_err:
    print_red(opt_err)
    print(usage_msg)
    exit(2)

def check_file_existance(path):
    if os.path.exists(path):
        return
    else:
        print_red("\nFile '{}' does not exist!".format(path))
        exit(1)

cutoff = False
merge_reads = False
V3V4 = False
quality_plot = False
primer_path = None
outdir_path = "{}{}preprocess16S_result_{}".format(os.getcwd(), os.sep, now).replace(" ", "_") # default path
read_paths = dict()
for opt, arg in opts:
    if opt in ("-h", "--help"):
        print(usage_msg)
        exit(0)
    elif opt in ("-c", "--cutoff"):
        cutoff = True
    elif opt in ("-m", "--merge-reads"):
        merge_reads = True
    elif opt == "--V3-V4":
        V3V4 = True
    elif opt in ("-q", "--quality-plot"):
        quality_plot = True
    elif opt in ("-p", "--primers"):
        check_file_existance(arg)
        primer_path = arg
    elif opt in ("-o", "--outdir"):
        outdir_path = os.path.abspath(arg)
    elif opt in ("-1", "--R1"):
        if not "-2" in argv and not "--R2" in argv:
            print_red("ATTENTION!\n\tYou should specify both forward and reverse reads!")
            exit(1)
        check_file_existance(arg)
        read_paths["R1"] = arg
    elif opt in ("-2", "--R2"):
        if not "-1" in argv and not "--R1" in argv:
            print_red("\nATTENTION!\n\tYou should specify both forward and reverse reads!")
            exit(1)
        check_file_existance(arg)
        read_paths["R2"] = arg


# --V3-V4 option can't be specified without --merge-reads (-m)
if V3V4 and not merge_reads:
    print_red("\nError!\a")
    print_red("'--V3-V4' option can't be specified without '--merge-reads' (-m) option")
    print("""The reason for this is that checking for constant region 
    between V3 and V4 variable regions is performed only while merging reads and 
    only if your target sequences contain V3 and V4 regions 
    and a constant region between them.\n""")
    exit(1)


# Check packages needed for plotting
if quality_plot:
    print("\nChecking packages needed for plotting...")
    imp_error_msg = """\tAttention!\nPKG package is not installed. Graph will not be plotted
If you want to use this feature, please install PKG (e.g. pip install PKG)"""
    try:
        import numpy as np
    except ImportError as imperr:
        print_yellow("\nError: {}".format(str(imperr)))
        print_yellow(imp_error_msg.replace("PKG", "numpy"))
        quality_plot = False
    try:
        import matplotlib.pyplot as plt
    except ImportError as imperr:
        print_yellow("\nError: {}".format(str(imperr)))
        print_yellow(imp_error_msg.replace("PKG", "matplotlib"))
        quality_plot = False
    print_green("ok")
    del imp_error_msg
            

# Check utilities for read merging
if merge_reads:
    pathdirs = os.environ["PATH"].split(os.pathsep)
    for utility in ("fasta36", "blastn", "blastdbcmd"):
        utility_found = False
        for directory in pathdirs:
            if utility in os.listdir(directory):
                utility_found = True
                break
        if not utility_found:
            print_yellow("\tAttention!\n{} is not found in your system. Reads will not be merged.".format(utility))
            print_yellow("If you want to use this feature, please install {}".format(utility))
            print_yellow("""If this error still occure although you have installed everything 
    -- make sure that this program is added to PATH)""")
            merge_reads = False

    try:
        import read_merging_16S
    except ImportError as imperr:
        print_yellow("\nError: {}".format(str(imperr)))
        print_yellow("\tModule 'read_merging_16S' not found. Reads will not be merged.")
        print_yellow("To use this feature, you need to install 'read_merging_16S'.")
        print_yellow("More precisely, you need to run 'install.sh' provided with 'preprocess_16S.py'.")
        print_yellow("For further info see README.md in repo.")
        merge_reads = False


from re import match
from re import findall
from gzip import open as open_as_gzip
from gzip import GzipFile
from _io import TextIOWrapper
from collections import namedtuple

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

def close_files(*files):
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
    for obj in files:
        if isinstance(obj, dict) or isinstance(obj, list) or isinstance(obj, tuple):
            for indx_or_key in obj:
                obj[indx_or_key].close()
        elif isinstance(obj, TextIOWrapper) or isinstance(obj, GzipFile):
            obj.close()
        else:
            print_yellow("""If you use function 'close_files', please, store file objects in 
    lists, tuples or dictionaries or pass file objects itself to the function.""")
            try:
                obj.close()
            except:
                print_red("Object passed to 'close_files' function can't be closed")
                exit(1)


def write_fastq_record(outfile, fastq_record):
    """
    :param outfile: file, which data from fastq_record is written in
    :type outfile: _io.TextIOWrapper
    :param fastq_record: dict of 4 elements. Elements are four corresponding lines of fastq
    :type fastq_record: dict<str: str>
    :return: void
    """
    try:
        outfile.write(fastq_record["seq_id"] + '\n')
    except TypeError:
        print(fastq_record)
        exit(1)
    outfile.write(fastq_record["seq"] + '\n')
    outfile.write(fastq_record["optional_id"] + '\n')
    outfile.write(fastq_record["quality_str"] + '\n')


def read_fastq_record(read_files):

    if len(read_files) != 1 and len(read_files) != 2:
        print_red("You can pass only 1 or 2 files to the function 'read_pair_of_reads'!\a")
        exit(1)

    fastq_recs = dict()           # this dict should consist of two fastq-records: from R1 and from R2
    for key in read_files.keys():
        fastq_recs[key] = {                    #read all 4 lines of fastq-record
            "seq_id": read_files[key].readline(),
            "seq": read_files[key].readline(),
            "optional_id": read_files[key].readline(),
            "quality_str": read_files[key].readline()
        }
    return fastq_recs


def find_primer(primers, fastq_rec):
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

    for primer in primers:
        primer_len = len(primer)

        for shift in range(0, MAX_SHIFT + 1):
            score_1, score_2 = 0, 0

            for pos in range(0, primer_len - shift):
                # if match, 1(one) will be added to score, 0(zero) otherwise
                score_1 += int(primer[pos+shift] in MATCH_DICT[read[pos]])
            if score_1 / primer_len >= RECOGN_PERCENTAGE:
                if not cutoff:
                    return (True, fastq_rec)
                else:
                    cutlen = len(primer) - shift
                    fastq_rec["seq"] = read[cutlen : ]
                    fastq_rec["quality_str"] = fastq_rec["quality_str"][cutlen : ]
                    return (True, fastq_rec)

            # If there is no primer in 0-position, loop below will do unnecessary work.
            # Let it go a bit further. 
            shift += 1 

            for pos in range(0, primer_len - shift):
                score_2 += int(primer[pos] in MATCH_DICT[read[pos + shift]])
            if score_2 / primer_len >= RECOGN_PERCENTAGE:
                if not cutoff:
                    return (True, fastq_rec)
                else:
                    cutlen = len(primer) + shift
                    fastq_rec["seq"] = read[cutlen : ]
                    fastq_rec["quality_str"] = fastq_rec["quality_str"][cutlen : ]
                    return (True, fastq_rec)

    return (False, fastq_rec)


def select_file_manually(message):
    """
    Function provides you with an interface of selecting a file by typing path to it to the command line.

    :param message: message shown on the screen offering to enter the path
    :type message: str
    :return: path to manually selected file
    :return type: str
    """
    print('\n' + '~' * 50 + '\n')
    while True:
        path = input(message)
        if path == 'q!':
            exit(0)
        if os.path.exists(path):
            print_green("\tok...")
            return path
        else:
            print_red("\tERROR\tThere is no file named \'{}\'\a".format(path))



# This is a decorator.
# All functions that process reads do some same operations: they open read files, open result files, 
#   they count how many reads are already processed and they close all files mentioned above.
# This is why I use decorator here, although in some cases I didn't find the way how to implement it smartly enough.
# But the temptation was too strong to resist ;)
def progress_counter(process_func, read_paths, result_paths=None, stats=None, inc_percentage=0.05):

    def organizer():

        # Collect some info
        file_type = int(match(r".*\.gz$", read_paths["R1"]) is not None)
        how_to_open = OPEN_FUNCS[file_type]
        actual_format_func = FORMATTING_FUNCS[file_type]
        readfile_length = sum(1 for line in how_to_open(read_paths["R1"], 'r'))  # lengths of read files are equal
        read_pairs_num = int(readfile_length / 4)            # divizion by 4, because there are 4 lines per one fastq-record

        # Open files
        read_files = dict()
        result_files = dict()
        try:
            for key in read_paths.keys():
                read_files[key] = how_to_open(read_paths[key])
            if result_paths is not None:
                for key in result_paths.keys():
                    result_files[key] = open(result_paths[key], 'w')
        except OSError as oserror:
            print_red("Error while opening file", str(oserror))
            close_files(read_files, result_files)
            input("Press enter to exit:")
            exit(1)

        # Proceed
        inc_percentage = 0.01
        reads_processed = 0
        spaces = 50
        next_done_percentage = inc_percentage
        print("\nProceeding...\n\n")
        print("[" + " "*50 + "]" + "  0%\r", end="")
        while reads_processed < read_pairs_num:

            fastq_recs = read_fastq_record(read_files)

            for file_key in fastq_recs.keys():
                for field_key in fastq_recs[file_key].keys():
                    fastq_recs[file_key][field_key] = actual_format_func(fastq_recs[file_key][field_key])

            # Do what you need with these reads
            if result_paths is not None:
                process_func(fastq_recs, result_files, stats)
            else:
                process_func(fastq_recs)    # it is None while calulating data for plotting

            reads_processed += 1
            if reads_processed / read_pairs_num >= next_done_percentage:
                count = round(next_done_percentage * 100)
                spaces = 50 - int(count/2)

                print("[" + "="*int(count/2) + ">" + " "*spaces + "]" + "  {}% ({}/{} read pairs are processed)\r"
                    .format(count, reads_processed, read_pairs_num), end="")
                next_done_percentage += inc_percentage

        print("[" + "="*50 + "]" + "  100% ({}/{} read pairs are processed)\n\n"
            .format(reads_processed, read_pairs_num))
        close_files(read_files, result_files)

    return organizer


def find_primer_organizer(fastq_recs, result_files, stats):

    primer_in_R1, fastq_recs["R1"] = find_primer(primers, fastq_recs["R1"])
    primer_in_R2, fastq_recs["R2"] = find_primer(primers, fastq_recs["R2"])

    if primer_in_R1 or primer_in_R2:
        write_fastq_record(result_files["mR1"], fastq_recs["R1"])
        write_fastq_record(result_files["mR2"], fastq_recs["R2"])
        stats["match"] += 1
    else:
        write_fastq_record(result_files["trR1"], fastq_recs["R1"])
        write_fastq_record(result_files["trR2"], fastq_recs["R2"])
        stats["trash"] += 1



# |======================= Start calculating =======================|


# |===== Prepare files for read merging =====|

# === Select primer file if it is not specified ===

if primer_path is None:     # if primer file is not specified by CL argument
    # search for primer file in current directory
    print('\n' + '~' * 50 + '\n')
    for smth in os.listdir('.'):
        if match(r".*[pP]rimers.*\.(m)?fa(sta)?$", smth) is not None:
            primer_path = smth
            break
    reply = None
    if primer_path is not None:
        with open(primer_path, 'r') as primer_file:
            file_content = primer_file.read()
        message = """File named '{}' is found.\n Here are primers stored in this file: 
    \n{} \nFile \'{}\' will be used as primer file.
    Do you want to select another file manually instead of it?
        Enter \'y\' to select other file,
        \'q!\' -- to exit
        or anything else (e.g. merely press ENTER) to continue:""".format(primer_path, file_content, primer_path)
    else:
        message = """No file considered as primer file is found.
    Do you want to select file manually?
    Enter \'y\' to select other file or anything else (e.g. merely press ENTER) to exit:"""
    reply = input(message)
    if reply == "q!":
        print("Exiting...")
        exit(0)
    if reply == 'y':
        message = "Enter path to .fa primer file (or \'q!\' to exit):"
        primer_path = select_file_manually(message)
    else:
        if primer_path is not None:
            pass
        else:
            exit(0)
    print('\n' + '~' * 50 + '\n')


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
try:
    primer_file = how_to_open(primer_path, 'r')
    for i, line in enumerate(primer_file):
        if i % 2 == 0:                                              # if line is sequense id
            primer_ids.append(actual_format_func(line))
        else:                                                       # if line is a sequence
            line = actual_format_func(line)       
            err_set = set(findall(r"[^ATGCRYSWKMBDHVN]", line))     # primer validation
            if len(err_set) != 0:
                print_red("! There are some inappropriate symbols in your primers. Here they are:")
                for err_symb in err_set:
                    print(err_symb)
                print("See you later with correct data")
                input("Press enter to exit:")
                exit(1)
            primers.append(line)
except OSError as oserror:
    print_red("Error while reading primer file.\n", str(oserror))
    input("Press enter to exit:")
    exit(1)
finally:
    if primer_file is not None:
        primer_file.close()

# The last line in fasta file can be a blank line.
# We do not need it.
try:
    primer_ids.remove('')
except ValueError:
    pass


# === Select read files if they are not specified ===

# If read files are not specified by CL arguments
if len(read_paths) == 0:
    # Search for read files in current directory.
    read_paths = dict()
    for smth in os.listdir('.'):
        if match(r".*_R1_.*\.fastq(\.gz)?$", smth) is not None and os.path.exists(smth.replace("_R1_", "_R2_")):
            read_paths["R1"] = smth
            read_paths["R2"] = smth.replace("_R1_", "_R2_")
            break
    reply = None
    if len(read_paths) == 0:
        message = """\nNo files considered as read files are found. 
        Do you want to select other files manually?
        Enter \'y\' to select other files or anything else (e.g. merely press ENTER) to exit:"""
    else:
        message = """Files named\n\t'{}',\n\t'{}'\n are found. 
        They will be used as read files. 
    Do you want to select other files manually instead of them?
        Enter \'y\' to select other files,
        \'q!\' -- to exit
        or anything else (e.g. merely press ENTER) to continue:""".format(read_paths["R1"], read_paths["R2"])
    reply = input(message)
    if reply == "q!":
        print("Exiting...")
        exit(0)
    if reply == "y":
        read_paths = dict()
        for R_12, forw_rev in zip(("R1", "R2"), ("forward", "reverse")):
            message = "Enter path to the .fastq file with {} reads (or \'q!\' to exit):".format(forw_rev)
            read_paths[R_12] = select_file_manually(message)
        # check if both of read files are of fastq format or both are gzipped
        check = 0
        for i in read_paths.keys():
            check += int(match(r".*\.gz$", read_paths[i]) is not None)
        if check == 1:
            print_red("""\n\tATTENTION!
    \tBoth of read files should be of fastq format or both be gzipped (i.e. fastq.gz).
    \tPlease, make sure this requirement is satisfied.""")
            input("Press enter to exit:")
            exit(1)
    else:
        if len(read_paths) == 2:
            pass
        else:
            exit(0)
    print('\n' + '~' * 50 + '\n')

# I need to keep names of read files in memory in order to name result files properly.
names = dict()
file_name_itself = read_paths["R1"][read_paths["R1"].rfind(os.sep)+1 :] # get rid of, e.g. '/../' in path
names["R1"] = match(r"(.*)\.f(ast)?q(\.gz)?$", file_name_itself).groups(0)[0]
file_name_itself = read_paths["R2"][read_paths["R2"].rfind(os.sep)+1 :]
names["R2"] = match(r"(.*)\.f(ast)?q(\.gz)?$", file_name_itself).groups(0)[0]


# === Create output directory. ===
if not os.path.exists(outdir_path):
    try:
        os.mkdir(outdir_path)
    except OSError as oserror:
        print_red("Error while creating result directory\n", str(oserror))
        close_files(read_files)
        input("Press enter to exit:")
        exit(1)



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



# |===== Start the process of searching for cross-talks =====|

print_yellow("\nSearching for cross-talks started")
primer_task = progress_counter(find_primer_organizer, read_paths, result_paths, primer_stats)

primer_task()

print_green("\nSearching for cross-talks is completed")
print("""\n{} read pairs with primer sequences are found.
{} read pairs considered as cross-talks are found.""".format(primer_stats["match"], primer_stats["trash"]))
print('\n' + '~' * 50 + '\n')

# |===== The process of searching for cross-talks is completed =====|



# |===== Start the process of read merging =====|

if merge_reads:

    merge_result_files = read_merging_16S.merge_reads(result_paths["mR1"], result_paths["mR2"], 
        outdir_path=outdir_path, V3V4=V3V4)
    merging_stats = read_merging_16S.get_merging_stats()

    files_to_gzip.extend(merge_result_files.values())


# |===== The process of merging reads is completed =====|


# |===== Prepare data for plotting =====|

if quality_plot:

    top_x_scale, step = 40.0, 0.5
    # average read quality
    X = np.arange(0.0, top_x_scale + step, step)
    # amount of reads with sertain average quality
    Y = np.zeros(int(top_x_scale / step), dtype=int)

    # This function will be used as "organizer"
    def add_data_for_qual_plot(fastq_recs):

        quality_strings = list()
        for rec in fastq_recs.values():
            quality_strings.append(rec["quality_str"])

        for qual_str in quality_strings:

            qual_array = np.array( [ord(qual_str[i]) - 33 for i in range(len(qual_str))])
            avg_qual = round(np.mean(qual_array), 2)
            min_indx = (np.abs( X - avg_qual )).argmin()

            Y[min_indx] += 1

    # Name without "__R1__" and "__R2__":
    more_common_name = os.path.basename(read_paths["R1"])[: os.path.basename(read_paths["R1"]).find("_R1_")]

    if merge_reads:
        data_plotting_paths = {
            "R1": merge_result_files["merg"]
        }
    else:
        data_plotting_paths = {
            "R1": result_paths["mR1"],
            "R2": result_paths["mR2"]
        }

    print_yellow("\nCalculations for plotting started")
    plotting_task = progress_counter(add_data_for_qual_plot, data_plotting_paths)
    plotting_task()
    print_green("\nCalculations for plotting are completed")
    print('\n' + '~' * 50 + '\n')


 # |===== Plot a graph =====|    

    image_path = os.sep.join((outdir_path, "quality_plot.png"))

    fig, ax = plt.subplots()
    ax.plot(X[:len(Y)], Y, 'r')
    ax.set(title="Quality per read distribution", 
        xlabel="avg quality (Phred33)",
        ylabel="amount of reads")
    ax.grid()
    fig.savefig(image_path)

    # Open image as image via xdg-open in order not to pause executing main script while 
    #   you are watching at pretty plot.

    os.system("xdg-open {}".format(image_path))
    # plt.show()


# |===== The process of plotting is completed =====|


# Remove empty files
for file in files_to_gzip:
    if os.stat(file).st_size == 0:
        os.remove(file)
        print("'{}' is removed since it is empty".format(file))

# Gzip result files
print_yellow("\nGzipping result files...")
for file in files_to_gzip:
    if os.path.exists(file):
        os.system("gzip -f {}".format(file))
        print("\'{}\' is gzipped".format(file))
print_green("Gzipping is completed\n")
print_yellow("Result files are placed in the following directory:\n\t{}\n\n".format(os.path.abspath(outdir_path)))


# Create log file
with open("{}{}preprocess16S_{}.log".format(outdir_path, os.sep, now).replace(" ", "_"), 'w') as logfile:
    logfile.write("Script 'preprocess_16S.py' was ran on {}\n".format(now.replace('.', ':')))
    logfile.write("The man who ran it was searching for following primer sequences:\n\n")
    for i in range(len(primers)):
        logfile.write("{}\n".format(primer_ids[i]))
        logfile.write("{}\n".format(primers[i]))
    logfile.write("""\n{} read pairs with primer sequences have been found.
{} read pairs without primer sequences have been found.\n""".format(primer_stats["match"], primer_stats["trash"]))
    if cutoff:
        logfile.write("These primers were cut off.\n")

    if merge_reads:
        logfile.write("\n\tReads were merged\n\n")
        logfile.write("{} read pairs have been merged.\n".format(merging_stats[0]))
        logfile.write("{} read pairs haven't been merged.\n".format(merging_stats[1]))
        logfile.write("{} read pairs have been considered as too short for merging.\n".format(merging_stats[2]))

    if quality_plot:
        logfile.write("\nQuality graph was plotted. Here it is: \n\t'{}'\n".format(image_path))

exit(0)

