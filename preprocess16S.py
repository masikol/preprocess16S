"""This script preprocesses reads from 16S regions of rDNA. It works with Illumina pair-end reads.
Excactly: it detects and removes reads, that came from other sample, relying on the information,
    whether there are PCR primer sequences in these reads or not. If required primer sequence is
    found in a read, therefore, it is a read from 16S rDNA and we need it.
Moreover, it can cut these primers off and merge reads by using 'read_merging_16S' module.
Reads should be stored in files, both of which are of fastq format or both are gzipped (i.e. fastq.gz).
Sequences of required primers are retrieved from .fa or fa.gz file.

Attention! This script cannot be executed by python interpreter version < 3.0!

Usage:
    python preprocess16S.py [--cutoff] [-p primer_file] [-1 forward_reads -2 reverse_reads] [-o output_dir]
Options:
    -c or --cutoff 
        cut primer sequences off;
    -m or --merge-reads
        all of a sudden, if this option is specified, script will merge reads together
    -q or --quality-plot
        plot a graph number of reads as a function of average read quality
    -p or --primers
        file, in which primer sequences are stored;
    -1 or --R1
        file, in which forward reads are stored;
    -2 or --R2
        file, in which reverse reads are stored;
    -o or --outdir
        directory, in which result files will be plased.

Last modified 26.04.2019
"""


# |===== Check python interpreter version. =====|

# It won't work on python2 because I use 'input()' function for data input by user.
from sys import version_info
if (version_info.major + 0.1 * version_info.minor) < (3.0 - 1e-6):        # just-in-case 1e-6 substraction
    print("\t\nATTENTION!\nThis script cannot be executed by python interpreter version < 3.0!")
    print("\tYour python version: {}.{}.".format(version_info.major, version_info.minor))
    print("Try '$ python3 preprocess16S.py' or update your interpreter.")
    exit(0)


from datetime import datetime
now = datetime.now().strftime("%Y-%m-%d %H.%M.%S")


# |===== Handle CL arguments =====|

import os
import getopt
from sys import argv
usage_msg = """usage:
    python preprocess16S.py [--cutoff] [--merge-reads] -p <primer_file> -1 <forward_reads> -2 <reverse_reads> [-o output_dir]"""

try:
    opts, args = getopt.getopt(argv[1:], "hcmqp:1:2:o:", ["help", "cutoff", "merge-reads", "quality-plot", "primers=", "R1=", "R2=", "outdir="])
except getopt.GetoptError as opt_err:
    print(opt_err)
    print(usage_msg)
    exit(2)

def check_file_existance(path):
    if os.path.exists(path):
        return
    else:
        print("\nFile '{}' does not exist!".format(path))
        exit(1)

cutoff = False
merge_reads = False
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
    elif opt in ("-q", "--quality-plot"):
        quality_plot = True
    elif opt in ("-p", "--primers"):
        check_file_existance(arg)
        primer_path = arg
    elif opt in ("-o", "--outdir"):
        outdir_path = os.path.abspath(arg)
    elif opt in ("-1", "--R1"):
        if not "-2" in argv and not "--R2" in argv:
            print("ATTENTION!\n\tYou should specify both forward and reverse reads!")
            exit(1)
        check_file_existance(arg)
        read_paths["R1"] = arg
    elif opt in ("-2", "--R2"):
        if not "-1" in argv and not "--R1" in argv:
            print("ATTENTION!\n\tYou should specify both forward and reverse reads!")
            exit(1)
        check_file_existance(arg)
        read_paths["R2"] = arg


from re import match
from re import findall
from gzip import open as open_as_gzip
from gzip import GzipFile
from _io import TextIOWrapper 

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

RECOGN_PERCENTAGE = 0.51
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
            print("""If you use function 'close_files', please, store file objects in 
    lists, tuples or dictionaries or pass file objects itself to the function.""")


def write_fastq_record(outfile, fastq_record):
    """
    :param outfile: file, which data from fastq_record is written in
    :type outfile: _io.TextIOWrapper
    :param fastq_record: dict of 4 elements. Elements are four corresponding lines of fastq
    :type fastq_record: dict<str: str>
    :return: void
    """
    outfile.write(fastq_record["seq_id"] + '\n')
    outfile.write(fastq_record["seq"] + '\n')
    outfile.write(fastq_record["optional_id"] + '\n')
    outfile.write(fastq_record["quality_str"] + '\n')


def find_primer(primers, read):
    """
    This function figures out, whether a primer sequence is in read passed to it.
    Pluralistically it cuts primer sequences off if there are any.

    :param primers: list of required primer sequences
    :type primers: list<str>
    :param read: read, which this function searches for required primer in
    :type read: str
    :return: if primer sequence is found in read -- returns True, else returns False
    Moreover, returns a read with a primer sequence cut off if there were any (and if -c is specified)
        and intact read otherwise.
    """

    for primer in primers:
        primer_len = len(primer)
        for shift in range(0, MAX_SHIFT + 1):
            score_1, score_2 = 0, 0
            for pos in range(0, primer_len - shift):
                # if match, 1(one) will be added to score, 0(zero) otherwise
                score_1 += int(primer[pos+shift] in MATCH_DICT[read[pos]])
            if score_1 / primer_len >= RECOGN_PERCENTAGE:
                if not cutoff:
                    return (True, read)
                return (True, read[len(primer) - shift : ])
            for pos in range(0, primer_len - shift):
                score_2 += int(primer[pos] in MATCH_DICT[read[pos + shift]])
            if score_2 / primer_len >= RECOGN_PERCENTAGE:
                if not cutoff:
                    return (True, read)
                return (True, read[len(primer) + shift : ])
    return (False, read)


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
            print("\tok...")
            return path
        else:
            print("\tERROR\tThere is no file named \'{}\'".format(path))


# If we do not want to merge reads, therefore, we do not need following objects to be stored in RAM:
if merge_reads:

    # I import it here in order not to make another if-statement
    import read_merging_16S

    # The followig dictionary represents some statistics of merging. 
    # Again, I define it here in order not to make another if-statement
    # Keys in this dictionary are in consistency with return codes of the function "handle_read_merging_result"
    merging_stats = {
        0: 0,           # number of merged reads
        1: 0,           # number of putative chimeras
        2: 0            # number of too short sequences
    }

    def handle_read_merging_result(merging_result, fastq_recs, result_files, read_files):
        """
        This function handles the result of read merging and writes sequences in corresponding files.

        :param merging_result: a tuple of strings (if reads are merged), or an integer, if they cannot be merged
        :param fastq_reqs: a dictionary of two fastq-records stored as dictionary of it's fields
        :type fastq_reads: dict<str: dict<str, str>>
        :type result_files: dict<str: _io.TextIOWrapper>
        :param read_files: this parameter is needed only for closing it if a fatal error occured
        :type read_files: dict<str: _io.TextIOWrapper>
        """

        # if reads are merged
        if merging_result == 0:

            merged_rec = {
                "seq_id": fastq_recs["R1"]["seq_id"],
                "seq": read_merging_16S.curr_merged_seq,
                "optional_id": '+',
                "quality_str": read_merging_16S.curr_merged_qual
            }
            if quality_plot:
                add_data_for_qual_plot(read_merging_16S.curr_merged_qual)
            write_fastq_record(result_files["merg"], merged_rec)
            return 0

        # if they probably are chimeras
        elif merging_result == 1:
            write_fastq_record(result_files["chR1"], fastq_recs["R1"])
            write_fastq_record(result_files["chR2"], fastq_recs["R2"])
            return 1

        # if resulting sequence is too short to distinguish taxa
        elif merging_result == 2:
            write_fastq_record(result_files["shrtR1"], fastq_recs["R1"])
            write_fastq_record(result_files["shrtR2"], fastq_recs["R2"])
            return 2

        # if unforseen situation occured in 'read_merging_16S'
        elif merging_result == 3:
            close_files(read_files, result_files)
            input("Press ENTER to exit:")
            exit(1)

        # if 'read_merging_16S' returnes something unexpected and undesigned
        else:
            print("ERROR!!!\n\tModule that was merging reads returned an unexpected and undesigned value.")
            try:
                print("\tHere it is: {}".format(str(merging_result)))
            except:
                print("\tUnfortunately, it cannot be printed. It is of {} type".format(type(merging_result)))
            print("It is my fault. Report it to me -- I will fix it.")
            close_files(read_files, result_files)
            input("Press ENTER to exit:")
            exit(1)


if quality_plot:

    try:
        import numpy as np
    except ImportError as imperr:
        print(str(immperr))
        print("Please, make sure that numpy is installed")
        close_files(read_files, result_files)
        exit(1)

    top_x_scale, step = 40.0, 0.5
    # average read quality
    X = np.arange(0.0, top_x_scale + step, step)
    # amount of reads with sertain average quality
    Y = np.zeros(int(top_x_scale / step), dtype=int)

    def add_data_for_qual_plot(*quality_strings):

        for qual_str in quality_strings:

            qual_array = np.array( [ord(qual_str[i]) - 33 for i in range(len(qual_str))])
            avg_qual = round(np.mean(qual_array), 2)
            min_indx = (np.abs( X - avg_qual )).argmin()

            Y[min_indx] += 1




# |======================= Start proceeding =======================|


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
                print("! There are some inappropriate symbols in your primers. Here they are:")
                for err_symb in err_set:
                    print(err_symb)
                print("See you later with correct data")
                input("Press enter to exit:")
                exit(1)
            primers.append(line)
except OSError as oserror:
    print("Error while reading primer file.\n", str(oserror))
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
            print("""\n\tATTENTION!
    \tBoth of read files should be of fastq format or both be gzipped (i.e. fastq.gz).
    \tPlease, make sure this requirement is satisfied.""")
            input("Press enter to exit:")
            exit(0)
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


# === Open read files ===

read_files = dict()
file_type = int(match(r".*\.gz$", read_paths["R1"]) is not None)
how_to_open = OPEN_FUNCS[file_type]
actual_format_func = FORMATTING_FUNCS[file_type]
print("Counting reads...")
readfile_length = sum(1 for line in how_to_open(read_paths["R1"], 'r'))  # lengths of read files are equal
print("""Done\nThere are {} reads in each file
(e.i. {} at all, since reads are pair-end)\n""".format(int(readfile_length/4), int(readfile_length/2)))
try:
    read_files["R1"] = how_to_open(read_paths["R1"])
    read_files["R2"] = how_to_open(read_paths["R2"])
except OSError as oserror:
    print("Error while opening one of .fastq read files.\n", str(oserror))
    close_files(read_files)
    input("Press enter to exit:")
    exit(1)


# === Create output directory. ===
if not os.path.exists(outdir_path):
    try:
        os.mkdir(outdir_path)
    except OSError as oserror:
        print("Error while creating result directory\n", str(oserror))
        close_files(read_files)
        input("Press enter to exit:")
        exit(1)

# If we want to merge reads -- create a directory for putative artifacts
if merge_reads:
    artif_dir = "{}{}putative_artifacts".format(outdir_path, os.sep)
    if not os.path.exists(artif_dir):
        try:
            os.mkdir(artif_dir)
        except OSError as oserror:
            print("Error while creating result directory\n", str(oserror))
            close_files(read_files)
            input("Press enter to exit:")
            exit(1)


# === Create and open result files. ===

result_files = dict()
# Keys description: 
# 'm' -- matched (i.e. sequence with primer in it); 
# 'tr' -- trash (i.e. sequence without primer in it);
# 'merg' -- merged sequences
# 'ch' -- putative chimeras
# 'shrt' -- those reads, that form too short sequence to distinguish taxa after merging 
# 'R1', 'R2' -- forward, reverse reads correspondingly;
# I need to keep these paths in memory in order to gzip corresponding files afterwards.
result_paths = {
    # We need trash anyway (trash without primers and, therefore, without 16S data):
    "trR1": "{}{}{}.trash.fastq".format(outdir_path, os.sep, names["R1"]),
    "trR2": "{}{}{}.trash.fastq".format(outdir_path, os.sep, names["R2"]),
}

if merge_reads:
    # Name without "__R1__" and "__R2__":
    more_common_name = names["R1"][: names["R1"].find("_R1_")]
    # File for merged sequences:
    result_paths["merg"] = "{}{}{}.16S.merged.fastq".format(outdir_path, os.sep, more_common_name)
    # Files for purativa chimeras:
    result_paths["chR1"] = "{}{}{}.chimera.fastq".format(artif_dir, os.sep, names["R1"])
    result_paths["chR2"] = "{}{}{}.chimera.fastq".format(artif_dir, os.sep, names["R2"])
    # Files for those reads, that form too short sequence after merging:
    result_paths["shrtR1"] = "{}{}{}.too_short.fastq".format(artif_dir, os.sep, names["R1"])
    result_paths["shrtR2"] = "{}{}{}.too_short.fastq".format(artif_dir, os.sep, names["R2"])
else:
    # Files for files, in which there are primer sequences:
    result_paths["mR1"] = "{}{}{}.16S.fastq".format(outdir_path, os.sep, names["R1"])
    result_paths["mR2"] = "{}{}{}.16S.fastq".format(outdir_path, os.sep, names["R2"])

try:
    for key in result_paths.keys():
        result_files[key] = open(result_paths[key], 'w')
except OSError as oserror:
    print("Error while opening one of result files", str(oserror))
    close_files(read_files, result_files)
    input("Press enter to exit:")
    exit(1)

# Some values for statistics:
m_count, tr_count = 0, 0
read_pairs_num = int(readfile_length / 4)            # divizion by 4, because there are 4 lines per one fastq-record
reads_processed, next_done_percentage = 0, 0.05


# |===== Start the process of searching for primer sequences in reads =====|

print("Proceeding...")
while reads_processed < read_pairs_num:
    try:
        fastq_recs = dict()           # this dict should consist of two fastq-records: from R1 and from R2
        for key in read_files.keys():
            fastq_recs[key] = { 	               #read all 4 lines of fastq-record
                "seq_id": read_files[key].readline(),
                "seq": read_files[key].readline(),
                "optional_id": read_files[key].readline(),
                "quality_str": read_files[key].readline()
            }
    except IOError as ioerror:
        print("Error while parsing one of the .fastq read files", str(ioerror))
        close_files(read_files, result_files)
        input("Press enter to exit:")
        exit(1)
    for file_key in fastq_recs.keys():
        for field_key in fastq_recs[file_key].keys():
            fastq_recs[file_key][field_key] = actual_format_func(fastq_recs[file_key][field_key])
        # As it has been mentioned, function find_primer() returns it's verdict about 
        #     whether there was any primer sequence in a read passed to it.
        # Moreover, it returns a read with a primer sequence cut off if there were any and intact read otherwise.
        # Therefore, usage is following albeit a bit awkward:
    primer_in_R1, fastq_recs["R1"]["seq"] = find_primer(primers, fastq_recs["R1"]["seq"])
    primer_in_R2, fastq_recs["R2"]["seq"] = find_primer(primers, fastq_recs["R2"]["seq"])
    try:
        if primer_in_R1 or primer_in_R2:
            if not merge_reads:
                write_fastq_record(result_files["mR1"], fastq_recs["R1"])
                write_fastq_record(result_files["mR2"], fastq_recs["R2"])
                if quality_plot:
                    add_data_for_qual_plot(fastq_recs["R1"]["quality_str"], fastq_recs["R2"]["quality_str"])
            else:
                merging_result = read_merging_16S.merge_reads(fastq_recs)
                reply = handle_read_merging_result(merging_result, fastq_recs, result_files, read_files)
                merging_stats[reply] += 1
            m_count += 1
        else:
            write_fastq_record(result_files["trR1"], fastq_recs["R1"])
            write_fastq_record(result_files["trR2"], fastq_recs["R2"])
            tr_count += 1
    except IOError as ioerror:
        print("Error while writing to one of the result files", str(ioerror))
        close_files(read_files, result_files)
        input("Press enter to exit:")
        exit(1)
    reads_processed += 1
    if reads_processed / read_pairs_num >= next_done_percentage:
        print("{}% of reads are processed\nProceeding...".format(round(next_done_percentage * 100)))
        next_done_percentage += 0.05
close_files(read_files, result_files)

print("100% of reads are processed")
print('\n' + '~' * 50 + '\n')
print("""{} read pairs with primer sequences are found.
{} read pairs without primer sequences are found.""".format(m_count, tr_count))

# Remove temporary and empty result files
if merge_reads:
    read_merging_16S.del_temp_files()
for file in result_paths.values():
    if os.stat(file).st_size == 0:
        os.remove(file)
        print("'{}' is removed since it is empty".format(file))

# We want to look at pretty plot while result files are gzipping
if quality_plot:
    gnuplot_script = "quality_plot.gp"
    data_file = "{}{}quality_data.tsv".format(outdir_path, os.sep)
    image_path = "{}{}quality_plot.png".format(outdir_path, os.sep)
    cmd_for_gnuplot = "gnuplot -c {} {} {}".format(gnuplot_script, data_file, image_path)
    with open(data_file, 'w') as qual_data_file:
        qual_data_file.write("Agv_read_quality,_Pherd33\tNumber-of-reads\n")
        i = 0
        while i < len(Y):
            qual_data_file.write("{}\t{}\n".format(X[i], Y[i]))
            i += 1
    os.system(cmd_for_gnuplot)
    os.system("xdg-open {}".format(image_path))

# Gzip result files
print("\nGzipping result files...")
for file in result_paths.values():
    if os.path.exists(file):
        os.system("gzip -f {}".format(file))
        print("\'{}\' is gzipped".format(file))
print("Done\n")
print("Result files are placed in the following directory:\n\t{}".format(os.path.abspath(outdir_path)))


# Create log file
with open("{}{}preprocess16S_{}.log".format(outdir_path, os.sep, now).replace(" ", "_"), 'w') as logfile:
    logfile.write("The script 'preprocess16S.py' was ran on {}\n".format(now.replace('.', ':')))
    logfile.write("The man who ran it was searching for following primer sequences:\n\n")
    for i in range(len(primers)):
        logfile.write("{}\n".format(primer_ids[i]))
        logfile.write("{}\n".format(primers[i]))
    logfile.write("""\n{} read pairs with primer sequences have been found.
{} read pairs without primer sequences have been found.\n""".format(m_count, tr_count))

    if merge_reads:
        logfile.write("\n\tReads were merged\n\n")
        logfile.write("{} read pairs have been merged.\n".format(merging_stats[0]))
        logfile.write("{} read pairs have been considered as chimeras.\n".format(merging_stats[1]))
        logfile.write("{} read pairs have been considered as too short for merging.\n".format(merging_stats[2]))

exit(0)
