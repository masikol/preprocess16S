#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "4.0.a"
# Year, month, day
# __last_update_date__ = "2020-08-07"

# |===== Check python interpreter version. =====|

import sys

if sys.version_info.major < 3:
    print( "\nYour python interpreter version is " + "%d.%d" % (sys.version_info.major,
        sys.version_info.minor) )
    print("   Please, use Python 3.\a")
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    if sys.platform.startswith("win"):
        raw_input("Press ENTER to exit:")
    # end if
    sys.exit(1)
else:
    print("\nPython {}.{}.{}".format(sys.version_info.major, sys.version_info.minor, sys.version_info.micro))
# end if


# |===== Handle CL arguments =====|

import os
import re
import getopt
from src.printing import *

try:
    opts, args = getopt.getopt(sys.argv[1:], "hvkmqr:1:2:o:t:f:N:m:p:",
        ["help", "version", "keep-primers", "merge-reads", "quality-plot",
        "primers=", "R1=", "R2=", "outdir=", "threads=", "phred-offset",
        "ngmerge-path=", "num-N=", "min-overlap=", "mismatch-frac="])
except getopt.GetoptError as opt_err:
    print( str(opt_err) )
    print("See help ('-h' option)")
    exit(2)
# end try


# First search for information-providing options:

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    print("\npreprocess16S.py -- version {}; {} edition.\n".format(__version__, __last_update_date__))
    if "--help" in sys.argv[1:]:
        print("""  DESCRIPTION:\n
preprocess16S.py -- this script is designed for preprocessing reads from 16S regions of rDNA.\n
It works with Illumina pair-end reads.\n
It detects and removes reads, that came from other sample (aka cross-talks),
    relying on the information, whether there are PCR primer sequences in these reads or not.\n
More precisely, if required primer sequence is found at the 5'-end of a read, therefore,
    it is a read from 16S rDNA and we need it.\n
Moreover, it can trim these primer sequences and merge reads by using 'read_merging_16S' module.\n
Reads should be stored in files, both of which are of FASTQ format or both are gzipped (i.e. fastq.gz).
Sequences of required primers are retrieved from FASTA file (or fasta.gz).\n
Attention! This script cannot be executed by python interpreter version < 3.0!""")
        print("----------------------------------------------------------\n")
    # end if

    print("""  OPTIONS:\n
-h (--help) --- show help message.
   '-h' -- brief, '--help' -- full;\n
-v (--version) --- show version;\n
-k (--keep-primers) --- Flag option. If specified, primer sequences will not be trimmed;\n
-m (--merge-reads) --- Flag option. If specified, reads will be merged together 
    with 'read_merging_16S' module;\n
-q (--quality-plot) --- Flag option. If specified, a plot of read quality distribution
    will be created. Requires 'numpy' and 'matplotlib' to be installed;\n
-r (--primers) --- FASTA file, in which primer sequences are stored.
    Illumina V3-V4 primer sequences are used by default:
    https://support.illumina.com/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf
    See "Amplicon primers" section.\n
-1 (--R1) --- FASTQ file, in which forward reads are stored;\n
-2 (--R2) --- FASTQ file, in which reverse reads are stored;\n
-o (--outdir) --- directory, in which result files will be placed;\n
-t (--threads) <int> --- number of threads to launch;
-f (--phred-offset) [33, 64] --- Phred quality offset (default -- 33);\n""")

    if "--help" in sys.argv[1:]:
        print("----------------------------------------------------------\n")
        print("""  EXAMPLES:\n
1) Filter cross-talks from files 'forw_R1_reads.fastq.gz' and 'rev_R2_reads.fastq.gz' depending on
 default Illumina V3-V4 primer sequenes. Do not trim primer sequences and put result files in
 the directory named 'outdir':\n
   ./preprocess16S.py -1 forw_R1_reads.fastq.gz -2 rev_R2_reads.fastq.gz -o outdir -c\n
2) Filter cross-talks from files 'forw_R1_reads.fastq.gz' and 'rev_R2_reads.fastq.gz' depending on
 primer sequenes in file 'V3V4_primers.fasta'. Trim primer sequences, merge reads,
 create a quality plot and put result files in the directory named 'outdir'. Launch 4 threads to calculations:\n
   ./preprocess16S.py -p V3V4_primers.fasta -1 forw_R1_reads.fastq.gz -2 rev_R2_reads.fastq.gz -o outdir -cmq -t 4\n""")
    # end if
    exit(0)
# end if

if "-v" in sys.argv[1:] or "--version" in sys.argv[1:]:
    print(__version__)
    exit(0)
# end if


keep_primers = False
merge_reads = False
quality_plot = False
primer_path = None
outdir_path = "{}{}preprocess16S_result_{}".format(os.getcwd(), os.sep, start_time).replace(" ", "_") # default path
read_paths = dict()
n_thr = 1
phred_offset = 33

ngmerge = os.path.join(os.path.dirname(os.path.abspath(__file__)), "binaries", "NGmerge")
# Length of (consensus?) conservative region between V3 and V4 in bacteria:
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0007401
num_N = 35
min_overlap = 20 # as default in NGmerge
mismatch_frac = 0.1 # as default in NGmerge

for opt, arg in opts:

    if opt in ("-k", "--keep-primers"):
        keep_primers = True

    elif opt in ("-m", "--merge-reads"):
        merge_reads = True

    elif opt in ("-q", "--quality-plot"):
        quality_plot = True

    elif opt in ("-r", "--primers"):
        if not os.path.exists(arg):
            print_error("File '{}' does not exist!".format(arg))
            sys.exit(1)
        # end if
        primer_path = arg

    elif opt in ("-o", "--outdir"):
        outdir_path = os.path.abspath(arg)

    elif opt in ("-1", "--R1"):
        if not "-2" in sys.argv and not "--R2" in sys.argv:
            print_error("You should specify both forward and reverse reads!")
            exit(1)
        # end if
        if not os.path.exists(arg):
            print_error("File '{}' does not exist!".format(arg))
            sys.exit(1)
        # end if
        read_paths["R1"] = arg

    elif opt in ("-2", "--R2"):
        if not "-1" in sys.argv and not "--R1" in sys.argv:
            print_error("You should specify both forward and reverse reads!")
            exit(1)
        # end if
        if not os.path.exists(arg):
            print_error("File '{}' does not exist!".format(arg))
            sys.exit(1)
        # end if
        read_paths["R2"] = arg

    elif opt in ("-t", "--threads"):
        try:
            n_thr = int(arg)
            if n_thr < 1:
                raise ValueError
            # end if
        except ValueError:
            print_error("number of threads must be positive integer number!")
            print(" And here is your value: '{}'".format(arg))
            exit(1)
        # end try

        if n_thr > len(os.sched_getaffinity(0)):
            print("""\nWarning! You have specified {} threads to use
  although {} are available.""".format(n_thr, len(os.sched_getaffinity(0))))
            error = True
            while error:
                reply = input("""\nPress ENTER to switch to {} threads,
  or enter 'c' to continue with {} threads,
  or enter 'q' to exit:>>""".format(len(os.sched_getaffinity(0)), n_thr))
                if reply in ("", 'c', 'q'):
                    error = False
                    if reply == "":
                        n_thr = len(os.sched_getaffinity(0))
                        print("\nNumber of threads switched to {}\n".format(n_thr))
                    elif reply == 'c':
                        pass
                    elif reply == 'q':
                        sys.exit(0)
                    # end if
                else:
                    print("\nInvalid reply!\n")
                # end if
            # end while
        # end if

    elif opt in ("-f", "--phred-offset"):
        try:
            phred_offset = int(arg)
            if phred_offset != 33 and phred_offset != 64:
                raise ValueError
            # end if
        except ValueError:
            print("\nError: invalid Phred offset specified: {}".format(phred_offset))
            print("Available values: 33, 64.")
            exit(1)
        # end try

    elif opt == "--ngmerge-path":

        if os.path.isfile(arg):
            ngmerge = arg
        else:
            print_error("file '{}' does not exist!".format(arg))
            sys.exit(1)
        # end if

    elif opt in ("-N", "--num-N"):
        try:
            num_N = int(arg)
            if num_N < 0:
                raise ValueError
            # end if
        except ValueError:
            print("Invalid minimum number of Ns (-N option): '{}'".format(arg))
            print("It must be integer number > 0.")
            sys.exit(1)
        # end try

    elif opt in ("-m", "--min-overlap"):
        try:
            min_overlap = int(arg)
            if min_overlap < 0:
                raise ValueError
            # end if
        except ValueError:
            print("Invalid minimum overlap (-m option): '{}'".format(arg))
            print("It must be integer number > 0.")
            sys.exit(1)
        # end try

    elif opt in ("-p", "--mismatch-frac"):
        try:
            mismatch_frac = float(arg)
            if mismatch_frac < 0 or mismatch_frac > 1:
                raise ValueError
            # end if
        except ValueError:
            print("Invalid minimum fraction of mismatches in the overlap (-p option): '{}'".format(arg))
            print("It must be fraction of the overlap length -- from 0.0 to 1.0.")
            sys.exit(1)
        # end try
    # end if
# end for


# Check if NGmerge is executable
if not os.access(ngmerge, os.X_OK):
    print_error("NGmerge file is not executable: '{}'".format(ngmerge))
    print("Please, make it executable. You can do it in this way:")
    print(" chmod +x {}".format(ngmerge))
# end if


print("\n |=== preprocess16S.py (version {}) ===|\n".format(__version__))

# Check packages needed for plotting
if quality_plot:

    print("\nChecking packages needed for plotting...")

    try:
        if n_thr == 1:
            from src.quality_plot import *
        else:
            from src.parallel_quality_plot import *
        # end if
    except ImportError as imperr:
        print_error("module integrity is corrupted.")
        print( str(imperr) )
        sys.exit(1)
    # end try
    print("ok")
# end if


# Check utilities for read merging
if merge_reads:

    pathdirs = os.environ["PATH"].split(os.pathsep)

    for utility in ("blastn", "blastdbcmd"):

        utility_found = False

        for directory in pathdirs:
            if os.path.exists(directory) and utility in os.listdir(directory):
                utility_found = True
                break
            # end if
        # end for

        if not utility_found:
            print_error("{} is not installed!".format(utility))
            print("Please install {}".format(utility))
            print("""If this error still occure although you have installed everything
    -- make sure that this program is added to PATH)""")
            exit(1)
        # end if
    # end for

    try:
        import read_merging_16S
    except ImportError as imperr:
        print_error("module integrity is corrupted.")
        print( str(imperr) )
        exit(1)
    # end try
# end if


print( '\n'+get_work_time() + " ({}) ".format(strftime("%d.%m.%Y %H:%M:%S", localtime(start_time))) + "- Start working\n")


from src.fastq import *
from src.filesystem import *
from src.crosstalks import *
import math


# This is a decorator.
# All functions that process reads do some same operations: they open read files, open result files,
#   they count how many reads are already processed.
# So we will use one interface during multiple procedures.

def progress_counter(process_func, read_paths, result_paths=None, **kwargs):

    def organizer():

        # Collect some info
        file_type = get_archv_fmt_indx(read_paths["R1"])
        how_to_open = OPEN_FUNCS[file_type]
        actual_format_func = FORMATTING_FUNCS[file_type]
        readfile_length = sum(1 for line in how_to_open(read_paths["R1"]))  # lengths of read files are equal
        read_pairs_num = int(readfile_length / 4)            # divizion by 4, since there are 4 lines per one fastq-record

        # Open files
        read_files = open_files(read_paths, how_to_open)
        if not result_paths is None:
            result_files = open_files(result_paths, open, 'w')
        # end if

        # Proceed
        inc_percentage = 0.01
        reads_processed = 0
        # Progress bar shoud not be larger than terminal
        readnum_digits = math.ceil(math.log(read_pairs_num, 10)) # number of digits in number of read pairs
        # Number of spaces-equals to print. 11 is number of some 'technical' characters in progress bar, such as '['.
        width = min(50, os.get_terminal_size().columns - (11+2*readnum_digits))
        width_ratio = 100 / width
        spaces = width
        next_done_percentage = inc_percentage
        printn("[" + " "*width + "]" + " 0%")

        while reads_processed < read_pairs_num:

            fastq_recs = read_fastq_pair(read_files, actual_format_func)

            # Do what you need with these reads
            if result_paths is not None:
                process_func(fastq_recs, result_files, **kwargs)
            else:
                process_func(fastq_recs, **kwargs)    # result_files is None while calulating data for plotting
            # end if

            reads_processed += 1
            if reads_processed / read_pairs_num >= next_done_percentage:
                width = min(50, os.get_terminal_size().columns - (11+2*readnum_digits))
                perc_done = round(next_done_percentage * 100)
                spaces = width - int(perc_done / width_ratio)

                printn("\r[{}>{}] {}% ({}/{})".format("="*int(perc_done / width_ratio), " "*spaces,
                    perc_done, reads_processed, read_pairs_num))
                next_done_percentage += inc_percentage
            # end if
        # end while

        print("\r[" + "="*width + "]" + " 100% ({}/{})\n\n"
            .format(reads_processed, read_pairs_num))

        close_files(read_files)
        if not result_paths is None:
            close_files(result_files)
        # end if
    # end def organizer

    return organizer
# end def progress_counter


# |======================= Proceed =======================|

# === Retrieve primers sequences from file ===

primers, primer_ids = get_primers(sys.argv[1:], primer_path)


# === Select read files if they are not specified ===

# If read files are not specified by CL arguments
if len(read_paths) == 0:
    # Search for read files in current directory.
    looks_like_forw_reads = lambda f: False if re.match(r".*R1.*\.f(ast)?q(\.gz)?$", smth) is None else True
    forw_reads = ( list(filter(looks_like_forw_reads, os.listdir('.'))) )[0]
    rev_reads = forw_reads.replace("R1", "R2")
    if not os.path.exists(rev_reads):
        print_error("File with reverse reads not found in current directory!")
        print("File with forward ones: '{}'".format(forw_reads))
        exit(1)
    # end if
    read_paths = {"R1": forw_reads, "R2": rev_reads}
# end if


# I need to keep names of read files in memory in order to name result files properly.
names = dict()
for key, path in read_paths.items():
    file_name_itself = os.path.basename(read_paths[key])
    names[key] = re.search(r"(.*)\.f(ast)?q", file_name_itself).group(1)
# end for


# === Create output directory. ===
if not os.path.exists(outdir_path):
    try:
        os.makedirs( outdir_path )
    except OSError as oserror:
        print_error("Error while creating result directory")
        print( str(oserror) )
        close_files(read_files)
        exit(1)
    # end try
# end if


# Check if there are old files in the output directory.
# If so -- ask user premission to remove them or quite. There is not neerd to append new reads to old data.
n_old_files = len(os.listdir(outdir_path))
if n_old_files != 0:
    error = True
    while error:
        reply = input("""Output directory is not empty.
Remove old content or quite? [R/q] """)
        if reply.upper() == 'R' or reply == "":
            error = False
            import shutil
            print()
            try:
                for fpath in os.listdir(outdir_path):
                    fpath = os.path.join(outdir_path, fpath)
                    print("Removing '{}'".format(fpath))
                    if os.path.isfile(fpath) or os.path.islink(fpath):
                        os.unlink(fpath) # remove plain file file or link
                    elif os.path.isdir(fpath):
                        shutil.rmtree(fpath) # remove directory tree
                    else:
                        print_error("error 55")
                        print("Please, contact the developer.")
                        sys.exit(55)
                    # end if
                # end for
            except OSError as oserr:
                print_error( str(oserr) )
                sys.exit(1)
            # end try
        elif reply.upper() == 'Q':
            sys.exit(0)
        else:
            print("Invalid reply: '{}'".format(reply))
        # end if
    # end while
# end if

artif_dir = os.path.join(outdir_path, "putative_artifacts")

# === Create directory for trash. ===
if not os.path.exists(artif_dir):
    try:
        os.makedirs( artif_dir )
    except OSError as oserror:
        print_error("Error while creating result directory")
        print( str(oserror) )
        close_files(read_files)
        exit(1)
    # end try
# end if


# === Create and open result files. ===

files_to_gzip = list()
result_files = dict()
# Keys description:
# 'm' -- matched (i.e. sequence with primer in it);
# 'tr' -- trash (i.e. sequence without primer in it);
# I need to keep these paths in memory in order to gzip corresponding files afterwards.
result_paths = {
    # We need trash anyway (trash without primers and, therefore, without 16S data):
    "mR1": "{}{}{}.16S.fastq".format( outdir_path, os.sep, names["R1"]),
    "mR2": "{}{}{}.16S.fastq".format( outdir_path, os.sep, names["R2"]),
    "trR1": "{}{}{}.trash.fastq".format( artif_dir, os.sep, names["R1"]),
    "trR2": "{}{}{}.trash.fastq".format( artif_dir, os.sep, names["R2"])
}

files_to_gzip.extend(result_paths.values())

primer_stats = {
    "match": 0,           # number of read pairs with primers
    "trash": 0            # number of cross-talks
}

print("\nFollowing files will be processed:")
for i, path in enumerate(read_paths.values()):
    print("  {}. '{}'".format(i+1, os.path.abspath(path)))
# end for

print("Number of threads: {};".format(n_thr))
print("Phred offset: {};".format(phred_offset))
if keep_primers:
    print("Primer sequences will not be trimmed.")
# end if
if merge_reads:
    print("Reads will be merged together.")
# end if
if quality_plot:
    print("Quality plot will be created.")
# end if
print('-'*20+'\n')


# |===== Start the process of searching for cross-talks =====|

print("{} - Searching for cross-talks started".format(get_work_time()))
print("Proceeding...\n")

primer_task = progress_counter(find_primer_organizer, read_paths, result_paths,
    primers=primers, stats=primer_stats, keep_primers=keep_primers)
primer_task()
del primer_task

print("{} - Searching for cross-talks is completed".format(get_work_time()))
print("""{} read pairs have been processed.
{} read pairs with primer sequences have been detected.
{} cross-talk read pairs have been detected.""".format(primer_stats["match"] + primer_stats["trash"],
    primer_stats["match"], primer_stats["trash"]))

cr_talk_rate = round(100 * primer_stats["trash"] / (primer_stats["match"] + primer_stats["trash"]), 3)
print("I.e. cross-talk rate is {}%".format(cr_talk_rate))
print('\n' + '~' * 50)

# |===== The process of searching for cross-talks is completed =====|


# |===== Start the process of read merging =====|

if merge_reads:

    merge_result_files = read_merging_16S.merge_reads(result_paths["mR1"], result_paths["mR2"], 
        ngmerge=ngmerge, outdir_path=outdir_path, n_thr=n_thr, phred_offset=phred_offset,
        num_N=num_N, min_overlap=min_overlap, mismatch_frac=mismatch_frac)
    merging_stats = read_merging_16S.get_merging_stats()

    files_to_gzip.extend(merge_result_files.values())
# end if


# |===== The process of merging reads is completed =====|


# |===== Create a quality plot =====|

if quality_plot:

    # Function for getting Q value from Phred-encoded character:
    def substr_phred_offs(q_symb):
        return ord(q_symb) - phred_offset
    # end def substr_phred_offs

    if merge_reads:
        data_plotting_paths = {
            "R1": merge_result_files["merg"]
        }
    else:
        data_plotting_paths = {
            "R1": result_paths["mR1"],
            "R2": result_paths["mR2"]
        }
    # end if

    print("\n{} - Calculations for plotting started".format(get_work_time()))
    print("Proceeding...\n")

    if n_thr == 1:
        plotting_task = progress_counter(calc_qual_disrib, data_plotting_paths,
            substr_phred_offs=substr_phred_offs)
        plotting_task()
        del plotting_task
        image_path = create_plot(outdir_path, phred_offset)
    else:
        Y = parallel_qual(data_plotting_paths, n_thr, substr_phred_offs)
        image_path = create_plot(Y, outdir_path, phred_offset)
    # end if

    print("{} - Calculations for plotting are completed".format(get_work_time()))
    print('\n' + '~' * 50)
# end if


# |===== The process of plotting is completed =====|


# Remove empty files
for file in files_to_gzip:
    if os.stat(file).st_size == 0:
        os.remove(file)
        print("'{}' is removed since it is empty".format(file))
    # end if
# end for

# Gzip result files
gzip_util = "gzip"
util_found = False
for directory in os.environ["PATH"].split(os.pathsep):
    if os.path.isdir(directory) and gzip_util in os.listdir(directory):
        util_found = True
        break
    # end if
# end for

print("\n{} - Gzipping result files...".format(get_work_time()))
for file in files_to_gzip:
    if os.path.exists(file):
        if util_found:
            os.system("{} {}".format(gzip_util, file))
        else:
            with open(file, 'r') as plain_file, open_as_gzip(file+".gz", 'wb') as gz_file:
                for line in plain_file:
                    gz_file.write(bytes(line, "utf-8"))
                # end for
            # end with
            os.unlink(file)
        # end if
        print("\'{}\' is gzipped".format(file))
    # end if
# end for

print("{} - Gzipping is completed\n".format(get_work_time()))
print("Result files are placed in the following directory:\n  '{}'\n".format(os.path.abspath(outdir_path)))

# Create log file
start_time_fmt = strftime("%d_%m_%Y_%H_%M_%S", localtime(start_time))
with open("{}{}preprocess16S_{}.log".format(outdir_path, os.sep, start_time_fmt).replace(" ", "_"), 'w') as logfile:

    logfile.write("Script 'preprocess_16S.py' was run at {}\n".format(start_time_fmt))
    end_time = strftime("%d_%m_%Y_%H_%M_%S", localtime(time()))
    logfile.write("Script completed it's job at {}\n".format(end_time))

    logfile.write("Phred offset: {}.\n\n".format(phred_offset))
    logfile.write("Number of threads used: {}.\n".format(n_thr))

    logfile.write("Primer sequences:\n\n")

    for i in range(len(primers)):
        logfile.write("{}\n".format(primer_ids[i]))
        logfile.write("{}\n".format(primers[i]))
    # end for

    logfile.write("\nFollowing files have been processed:\n")
    logfile.write("  '{}'\n  '{}'\n\n".format(os.path.abspath(read_paths["R1"]), os.path.abspath(read_paths["R2"])))

    logfile.write("""{} read pairs have been processed.
{} read pairs with primer sequences have been detected.
{} cross-talk read pairs have been detected.\n""".format(primer_stats["match"] + primer_stats["trash"],
    primer_stats["match"], primer_stats["trash"]))

    logfile.write("I.e. cross-talk rate is {}%\n\n".format(cr_talk_rate))

    if not keep_primers:
        logfile.write("Primer sequences were not trimmed.\n")
    # end if

    if merge_reads:
        logfile.write("\n\tReads were merged\n\n")
        logfile.write("{} read pairs have been merged.\n".format(merging_stats[0]))
        logfile.write("{} read pairs haven't been merged.\n".format(merging_stats[1]))
    # end if

    logfile.write("\nResults are in the following directory:\n  '{}'\n".format(outdir_path))

    if quality_plot:
        logfile.write("\nQuality plot was created. Here it is: \n  '{}'\n".format(image_path))
    # end if

    logfile.flush()
# end with

print(get_work_time() + " ({}) ".format(strftime("%d.%m.%Y %H:%M:%S", localtime(time()))) + "- Job is successfully completed!\n")
sys.exit(0)
