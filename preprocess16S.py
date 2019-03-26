"""Attention! This script cannot be executed by python interpreter version < 3.6!
If your interpreter throws an unhandled SyntaxError -- make sure you are using python 3.6+.

This script processes reads from 16S regions of rDNA.
Excactly: it removes those reads, that came from other loci, relying on the information,
    whether there are primer sequences in these reads or not. If required primer sequence is
    found in a read, therefore, it is a read from 16S rDNA and we need it.
Moreover, it cuts these primers off.
Reads should be stored in files, both of which are of fastq format or both be gzipped (i.e. fastq.gz).
Sequences of required primers are retrieved from .fa or fa.gz file.
It is recommended to place file with primers and files with reads
    in the same directory with this .py file.
It is recommended to name file with primers as 'primers.fasta' or something like it.
Excactly: the program will find file with primers automatically, if it's name
    contains word "primers" and has .fa, .fasta or .mfa extention.
Result files named '...16S.fastq.gz' and '...trash.fastq.gz' will be
    placed in the directory nested in directory, where this .py file is located.
This result directory will be named preprocess16S_result... and so on according to time it was ran.
Last modified 25.03.2019
"""

from re import match
from re import findall
from gzip import open as open_as_gzip
import os
from datetime import datetime

now = datetime.now().strftime("%Y-%m-%d %H.%M.%S")

MAX_SHIFT = 4
# Primer sequences are searched with respect to shift.
# E.g., if MAX_SHIFT == 3, primer sequence will be searched from this collocation:
# ---RRRRRRRRRRRRRR
# PPPPPPP
# to this collocation:
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


def close_all_files(read_files: dict, result_files: dict):
    """
    :param read_files: files with reads
    :type read_files: dict<str: _io.TextIOWrapper or gzip.GzipFile>
    :param result_files: files where reads are written in
    :type result_files: dict<str: _io.TextIOWrapper>
    :return: void
    """
    for key in read_files.keys():
        read_files[key].close()
    for key in result_files.keys():
        result_files[key].close()


from _io import TextIOWrapper    # import just for type hint
def write_fastq_record(outfile: TextIOWrapper, fastq_record: dict):
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


# This function decides, whether a primer sequence is in read passed to it.
# Pluralistically it cuts primer sequences off if there are any.
def find_primer(primers: list, read: str) -> (bool, str):
    """
    :param primers: list of required primer sequences
    :type primers: list<str>
    :param read: read, which this function searches for required primer in
    :type read: str
    :return: if primer sequence is found in read -- returns True, else returns False
    Moreover, returns a read with a primer sequence cut off if there were any and intact read otherwise.
    """
    global MAX_SHIFT
    global RECOGN_PERCENTAGE
    global MATCH_DICT
    for primer in primers:
        primer_len = len(primer)
        for shift in range(0, MAX_SHIFT + 1):
            score_1, score_2 = 0, 0
            for pos in range(0, primer_len - shift):
                # if match, 1(one) will be added to score, 0(zero) otherwise
                score_1 += int(primer[pos+shift] in MATCH_DICT[read[pos]])
            if score_1 / primer_len >= RECOGN_PERCENTAGE:
                return (True, read[len(primer) - shift : ])
            for pos in range(0, primer_len - shift):
                score_2 += int(primer[pos] in MATCH_DICT[read[pos + shift]])
            if score_2 / primer_len >= RECOGN_PERCENTAGE:
                return (True, read[len(primer) + shift : ])
    return (False, read)


def select_file_manually(message: str) -> str:
    """
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
            print(f"\tERROR\tThere is no file named \'{path}\'")


# search for file with primers in current directory
print('\n' + '~' * 50 + '\n')
primer_path = None
for smth in os.listdir('.'):
    if match(r".*[pP]rimers.*\.(m)?fa(sta)?$", smth) is not None:
        primer_path = smth
        break
reply = None
if primer_path is not None:
    with open(primer_path, 'r') as primer_file:
        file_content = primer_file.read()
    message = f"""File named '{primer_path}' is found.\n Here are primers stored in this file: 
\n{file_content} \nFile \'{primer_path}\' will be used as file with primers.
Do you want to select another file manually instead of it?
    Enter \'y\' to select other file,
    \'q!\' -- to exit
    or anything else (e.g. merely press ENTER) to continue:"""
else:
    message = """No file considered as file with primers is found.
Do you want to select file manually?
Enter \'y\' to select other file or anything else (e.g. merely press ENTER) to exit:"""
reply = input(message)
if reply == "q!":
    print("Exiting...")
    exit(0)
if reply == 'y':
    message = "Enter path to .fa file with primers (or \'q!\' to exit):"
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


# Retrieve primers sequences from file 
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
    print("Error while reading file with primers.\n", repr(oserror))
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


# Search for files with reads in current directory
read_paths = list()
for smth in os.listdir('.'):
    if match(r".*_R1_.*\.fastq(\.gz)?$", smth) is not None and os.path.exists(smth.replace("_R1_", "_R2_")):
        read_paths.append(smth)
        read_paths.append(smth.replace("_R1_", "_R2_"))
        break
reply = None
if len(read_paths) == 0:
    message = """\nNo files considered as files with reads are found. 
    Do you want to select other files manually?
    Enter \'y\' to select other files or anything else (e.g. merely press ENTER) to exit:"""
else:
    message = f"""Files named\n\t'{read_paths[0]}',\n\t'{read_paths[1]}'\n are found. 
    They will be used as files with reads. 
Do you want to select other files manually instead of them?
    Enter \'y\' to select other files,
    \'q!\' -- to exit
    or anything else (e.g. merely press ENTER) to continue:"""
reply = input(message)
if reply == "q!":
    print("Exiting...")
    exit(0)
if reply == "y":
    read_paths = list()
    for forw_rev in ("forward", "reverse"):
        message = f"Enter path to the .fastq file with {forw_rev} reads (or \'q!\' to exit):"
        read_paths.append(select_file_manually(message))
    # check if both of files with reads are of fastq format or both are gzipped
    check = 0
    for i in range(len(read_paths)):
        check += int(match(r".*\.gz$", read_paths[i]) is not None)
    if check == 1:
        print("""\n\tATTENTION!
\tBoth of files with reads should be of fastq format or both be gzipped (i.e. fastq.gz).
\tPlease, make sure this requirement is satisfied.""")
        input("Press enter to exit:")
        exit(0)
else:
    if len(read_paths) == 2:
        pass
    else:
        input("Press enter to exit:")
        exit(0)
print('\n' + '~' * 50 + '\n')

# I need to keep names of files with reads in memory in order to name result files properly.
names = list()
for path in read_paths:
    names.append(match(r"(.*)\.f(ast)?q(\.gz)?$", path).groups(0)[0])


# open files with reads
read_files = dict()
file_type = int(match(r".*\.gz$", read_paths[0]) is not None)
how_to_open = OPEN_FUNCS[file_type]
actual_format_func = FORMATTING_FUNCS[file_type]
print("Counting reads...")
readfile_length = sum(1 for line in how_to_open(read_paths[0], 'r'))  # lengths of files with reads are equal
print(f"""Done\nThere are {int(readfile_length/4)} reads in each file
(e.i. {int(readfile_length/2)} at all, since reads are pair-end)\n""")
try:
    read_files["R1"] = how_to_open(read_paths[0])
    read_files["R2"] = how_to_open(read_paths[1])
except OSError as oserror:
    print("Error while opening one of .fastq files with reads.\n", repr(oserror))
    for k in read_files.keys():
        read_files[k].close()
    input("Press enter to exit:")
    exit(1)


# Create result directory
outdir_path = f"{os.getcwd()}{os.sep}preprocess16S_result_{now}".replace(" ", "_")
try:
    os.mkdir(outdir_path)
except OSError as oserror:
    print("Error while creating result directory\n", repr(oserror))
    input("Press enter to exit:")
    exit(1)


# Create and open result files.
result_files = dict()
# Keys description: 
# 'm' -- matched (i.e. sequence with primer in it); 
# 'tr' -- trash (i.e. sequence without primer in it);
# 'R1', 'R2' -- forward, reverse reads correspondingly;
# I need to keep these paths in memory in order to gzip corresponding files afterwards.
result_paths = {
    "mR1": f"{outdir_path}{os.sep}{names[0]}.16S.fastq",
    "trR1": f"{outdir_path}{os.sep}{names[0]}.trash.fastq",
    "mR2": f"{outdir_path}{os.sep}{names[1]}.16S.fastq",
    "trR2": f"{outdir_path}{os.sep}{names[1]}.trash.fastq"
}
try:
    result_files["mR1"] = open(result_paths["mR1"], 'w')
    result_files["trR1"] = open(result_paths["trR1"], 'w')
    result_files["mR2"] = open(result_paths["mR2"], 'w')
    result_files["trR2"] = open(result_paths["trR2"], 'w')
except OSError as oserror:
    print("Error while opening one of result files", repr(oserror))
    close_all_files(read_files, result_files)
    input("Press enter to exit:")
    exit(1)


m_count, tr_count = 0, 0
# There are 4 lines per record in fastq file and we have two files, therefore:
reads_at_all = readfile_length / 2            # reads_at_all = (readfile_length / 4) * 2
reads_proceeded, next_done_percentage = 0, 0.05
# Start the process of searching for primer sequences in reads
print("Proceeding...")
for line in range(0, readfile_length, 4):
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
        print("Error while parsing one of the .fastq files with reads", repr(ioerror))
        close_all_files(read_files, result_files)
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
            write_fastq_record(result_files["mR1"], fastq_recs["R1"])
            write_fastq_record(result_files["mR2"], fastq_recs["R2"])
            m_count += 2
        else:
            write_fastq_record(result_files["trR1"], fastq_recs["R1"])
            write_fastq_record(result_files["trR2"], fastq_recs["R2"])
            tr_count += 2
    except IOError as ioerror:
        print("Error while writing to one of the result files", repr(ioerror))
        close_all_files(read_files, result_files)
        input("Press enter to exit:")
        exit(1)
    reads_proceeded += 2
    if reads_proceeded / reads_at_all >= next_done_percentage:
        print(f"{round(next_done_percentage * 100)}% of reads are processed\nProceeding...")
        next_done_percentage += 0.05
close_all_files(read_files, result_files)

print("100% of reads are processed")
print('\n' + '~' * 50 + '\n')
print(f"""{m_count} reads with primer sequences are found.
{tr_count} reads without primer sequences are found.""")

#gzip result files
print("Gzipping result files...")
for key in result_files.keys():
        os.system(f"gzip {result_paths[key]}")
        print(f"\'{result_paths[key]}\' is gzipped")
print("Done\n")
print(f"Result files are placed in the following directory:\n\t{os.path.abspath(outdir_path)}")

# create log file
with open(f"{outdir_path}{os.sep}preprocess16S_{now}.log".replace(" ", "_"), 'w') as logfile:
    logfile.write(f"The script 'preprocess16S.py' was ran on {now.replace('.', ':')}\n")
    logfile.write("The man who ran it was searching for following primer sequences:\n\n")
    for i in range(len(primers)):
        logfile.write(f"{primer_ids[i]}\n")
        logfile.write(f"{primers[i]}\n")
    logfile.write(f"""\n{m_count} reads with primer sequences have been found.
{tr_count} reads without primer sequences have been found.\n""")
input("Press enter to exit:")
exit(0)
