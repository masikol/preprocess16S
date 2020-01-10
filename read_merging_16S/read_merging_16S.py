#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "3.3.0"
# Year, month, day
__last_update_date__ = "2020-01-10"

import os
from re import search
from re import match
from gzip import open as open_as_gzip
from bz2 import open as open_as_bz2
from subprocess import Popen as sp_Popen, PIPE as sp_PIPE
import multiprocessing as mp


# |===============================  Data  ===============================|

# |===== Stuff for dealing with time =====|

from time import time, strftime, localtime, gmtime
start_time = time()

def get_work_time():
    return strftime("%H:%M:%S", gmtime( time() - start_time))
# end def get_work_time


# |===========================================|

def get_archv_fmt_indx(fpath):
    if fpath.endswith(".gz"):
        return 1
    elif fpath.endswith(".bz2"):
        return 2
    else:
        return 0
    # end if
# end def get_archv_fmt_indx

# Following tuples of functions are here in order to process '.fastq', '.fastq.gz' and '.fastq.bz2'
#    files using the same interface.

_OPEN_FUNCS = (open, open_as_gzip, open_as_bz2)

_FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format .fastq line
    lambda line: line.decode("utf-8").strip(),   # format .fastq.gz line
    lambda line: line.decode("utf-8").strip(),   # format .fastq.bz2 line (same as previous, but for consistency let it be so)
)

_ncbi_fmt_db = "REPLACE_DB"

class SilvaDBNotInstalledError(Exception):
    pass
# end class SilvaDBNotInstalledError


if not os.path.exists(_ncbi_fmt_db + ".nhr"):
    raise SilvaDBNotInstalledError("Silva database is not installed!\n  Please, run 'install_read_merging_16S.sh'")
# end if


def print_error(text):
    """Function for printing error messages"""
    print("\n   \a!! - ERROR: " + text + '\n')
# end def print_error


from sys import stdout as sys_stdout
def printn(text):
    """
    Function prints text to the console without adding '\n' in the end of the line.
    Why not just to use 'print(text, end="")'?
    In order to display informative error message if Python 2.X is launched
        instead if awful error traceback.
    """
    sys_stdout.write(text)
# end def printn


QSTART, QEND, SSTART, SEND, EVALUE, SACC, SSTRAND = range(7)
_cmd_for_blastn = """blastn -db {} -penalty -1 -reward 2 -ungapped \
-outfmt "6 qstart qend sstart send evalue sacc sstrand" \
-task megablast -max_target_seqs 1""".format(_ncbi_fmt_db)

_draft_cmd_for_blastdbcmd = "blastdbcmd -db {} -entry REPLACE_ME".format(_ncbi_fmt_db)

_RC_DICT = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
    'H': 'D', 'V': 'B', 'U': 'A', 'N': 'N'
}

# === Functionality to getting reverse-complement sequences ===
_single_nucl_rc = lambda nucl: _RC_DICT[nucl]
_rc = lambda seq: "".join(map(_single_nucl_rc, seq[::-1]))

# The next constant is used for detecting random alignments in the center of reads.

# It is calculated in the following way:
#               250_000 * 1/(4 ^ 14) = 0.00093 or 0.093%, where
#   250_000 is the number of reads, 4 is number of possible nucleotides, 14 is the length of a sequence.

# Assuming that 250_000 reads is common for a 16S metagenomic research,
#    0.00093 can be interprated as an expected number of 250_000 randomly generated sequences
#    of the length of 14, that tottally match some fixed sequence of the same length.

# 0.093% is little enough to kepp 14.
_MAX_ALIGN_OFFSET = 14
# --------------------------------------------

# Aligners can miss short overlapping region, especially if there are some sequencing errors at the end of reads
_MIN_OVERLAP = 11


# Maximum credible gap length. Very long gap isn't probable enough.
# It's half of E. coli's constant region between V3 and V4.
_MAX_CRED_GAP_LEN = 39


_curr_merged_seq = ""
_curr_merged_qual = ""


# The following dictionary represents some statistics of merging.
# Keys in this dictionary are in consistency with return codes of the function "_merge_pair"
_merging_stats = None    # it in None in the beginning, because there is no statistics before merging


# Constants for naive "aligning"
_SEED_LEN = _MIN_OVERLAP
_INDCS = range(_SEED_LEN)

# 5'-FFFFFFFFFFFFFFFFFFF-3'  -- forward read
#         3'-RRRRRRRRRRRRRRRRRRRR-5'  -- reverse read
# Length of overlap region above is ~136 nt, because:
# 1. Common read length is 301 nt.
# 2. Distance between start of forward primer and start of reverse primer is ~466 nt.
# 3. 466 - 301 = 165; -- overhang
# 4. 301 - 165 = 136; -- overlap, assumming that the picture above is +- symmetric.
# Length of constant region might vary, so I'll add 24 untill pretty and great enough number: 160.
_MAX_SHIFT = 160

# I can't choose and explain this numbert sensibly by now.
_MIN_PIDENT = 0.75


# ===============================  Internal functions  ===============================


def _close_files(*files):
    """
    Function closes all files passed to it, no matter whether they are single file objects or collenctions.

    :type of the argumets:
        dict<str: _io.TextIOWrapper or gzip.GzipFile>
        or list<_io.TextIOWrapper or gzip.GzipFile>
        or tuple<_io.TextIOWrapper or gzip.GzipFile>
        or _io.TextIOWrapper
        or gzip.GzipFile
    """
    for obj in files:
        if isinstance(obj, dict) or isinstance(obj, list) or isinstance(obj, tuple):
            for indx_or_key in obj:
                obj[indx_or_key].close()
            # end for
        
        elif isinstance(obj, TextIOWrapper) or isinstance(obj, GzipFile):
            obj.close()
        
        else:
            print("""If you use function 'close_files', please, store file objects in 
    lists, tuples or dictionaries or pass file objects itself to the function.""")
            try:
                obj.close()
            except:
                print_error("Object passed to 'close_files' function can't be closed")
                exit(1)
            # end try
        # end if
    # end for
# end def _close_files


def _write_fastq_record(outfile, fastq_record):
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
    outfile.write(fastq_record["seq_id"] + '\n')
    outfile.write(fastq_record["seq"] + '\n')
    outfile.write(fastq_record["opt_id"] + '\n')
    outfile.write(fastq_record["qual_str"] + '\n')
    outfile.flush()
# end def _write_fastq_record


def _read_fastq_pair(read_files, fmt_func):
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

    if len(read_files) != 2:
        print_error("You can only pass 2 files to the function 'read_fastq_pair'!")
        exit(1)
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
# end def _read_fastq_pair


def _naive_align(fseq, rseq):
    """
    This function takes the tail of forward read and slides it through the reverse read.
    If significant similarity is found -- return value of shift.
    Otherwise rerurns -1.

    :param fseq: forward read;
    :type fseq: str;
    :param rseq: reverse read;
    :type rseq: str;
    :return type: int;
    """

    shift = 0

    for seed_len in range(_SEED_LEN, _SEED_LEN + _MAX_SHIFT):

        seed = fseq[len(fseq) - seed_len :]
        score = 0

        for pos in range(seed_len):
            if seed[pos] == rseq[pos]:
                score += 1
            # end if
        # end for

        if score / seed_len > _MIN_PIDENT:
            return shift
        else:
            shift += 1
        # end if
    # end for

    return -1
# end def _naive_align


def _merge_by_overlap(loffset, overl, fseq, fqual, rseq, rqual):
    """
    Function merges two reads according to their overlapping redion.
    Function leaves nucleotide with higher quality.

    :param loffset: number of nucleotides in forward read from it's beginning to start of the overlapping region;
    :type loffset: int;
    :param overl: number of nucleotides in the overlapping region;
    :type overl: int;
    :param fseq: forward read itself;
    :type fseq: str;
    :param fqual: quality string of the forward read;
    :type fqual: str;
    :param rseq: reverse read itself;
    :type rseq: str;
    :param rqual: quality string of the reverse read;
    :type rqual: str;

    Returns tuple of merged sequence itself and it's quality line;
    Return type: tuple<str>;
    """

    # 5'-end comes from forward read
    merged_seq = fseq[: loffset]
    merged_qual = fqual[: loffset]

    # "recover" overlapping region depending on quality
    for i, j in zip(range(loffset, loffset + overl), range(overl)):
        if ord(fqual[i]) >= ord(rqual[j]):
            qual, nucl = fqual[i], fseq[i]
        else:
            qual, nucl = rqual[j], rseq[j]
        # end if
        
        merged_seq += nucl
        merged_qual += qual
    # end for

    # 3'-end comes from reverse read
    merged_seq += rseq[overl :]
    merged_qual += rqual[overl :]

    return (merged_seq, merged_qual)
# end def _merge_by_overlap


class NoRefAlignError(Exception):
    pass
# end class NoRefAlignError


def _blast_and_align(fseq, f_id, rseq, r_id):
    """
    Function does the following:
        1. Blast forward read against the database.
        2. Find and get reference sequence, against which forward read 
           has aligned with the better score (e. i. the first in list).
        3. Align reverse read agaist this reference sequence.
        4. Return results of these alignments.
    :param rseq: reverse read;
    :type rseq: str;

    Returns tuple of reports about aligning forward and reverse reads against 'reference' sequence.
    0-element of this tuple is 'forward-against-reference' aligning report.
    1-element of this tuple is 'reverse-against-reference' aligning report.
    Reports are tuples of 'str'.
    """

    # Align forward read. Find the best hit.
    pipe = sp_Popen(_cmd_for_blastn, shell=True, stdout=sp_PIPE, stderr=sp_PIPE, stdin=sp_PIPE)

    pipe.stdin.write( bytes(f_id.replace('@', '>') + '\n', "utf-8") )
    pipe.stdin.write( bytes(fseq + '\n', "utf-8") )
    pipe.stdin.close()

    exit_code = pipe.wait()
    if exit_code != 0:
        print_error("error while aligning a sequence against local database")
        print(pipe.stderr.read().decode("utf-8"))
        print("Id if erroneous read:\n  '{}'".format(f_id))
        exit(exit_code)
    # end if

    faref_report = pipe.stdout.read().decode("utf-8").strip().split('\t')

    # If there is no significant similarity
    if faref_report[0] == '':
        raise NoRefAlignError()
    # end if

    # Retrieve reference sequence
    cmd_for_blastbdcmd = _draft_cmd_for_blastdbcmd.replace("REPLACE_ME", faref_report[SACC])
    pipe = sp_Popen(cmd_for_blastbdcmd, shell=True, stdout=sp_PIPE, stderr=sp_PIPE)
    stdout_stderr = pipe.communicate()
    
    exit_code = pipe.returncode
    if exit_code != 0:
        print_error("error while aligning a sequence against local database")
        print(stdout_stderr[1].decode("utf-8"))
        exit(exit_code)
    # end if

    sbjct_fasta_lines = stdout_stderr[0].decode("utf-8").splitlines()
    sbjct_seq_id = sbjct_fasta_lines[0].strip() # get seq id
    sbjct_seq = "".join(sbjct_fasta_lines[1:]) # get sequence itself

    # Turn sequence around if forward read has aligned to minus-strand
    if faref_report[SSTRAND] == "minus":
        sbjct_seq = _rc(sbjct_seq)
    # end if

    # Align reverse read against reference sequence
    raref_report = SW_align(rseq, sbjct_seq)
    if raref_report is None:
        raise NoRefAlignError()
    # end if

    return faref_report, raref_report
# end def _blast_and_align


def _handle_unforseen_case(f_id, fseq, r_id, rseq):
    """
    Handle unforseen case and provide man who uses the program with information 
        about how erroneous reads align one against another.
    :param f_id: id of the forward sequence
    :type f_id: str
    :param fseq: forward sequence itself
    :type fseq: str
    :param r_id: id of the reverse sequence
    :type r_id: str
    :param rseq: reverse sequence itself
    :type rseq: str
    """
    error_report = "error_report.txt"
    print("Read merging crashed. Please contact the developer")
    print("You can get some info about reads caused this crash in the following file:\n\t'{}'"
        .format(os.path.abspath(error_report)))

    with open(error_report, 'w') as errfile:
        errfile.write("Forward read:\n\n" + f_id + '\n')
        errfile.write(fseq + "\n\n" + '~' * 40 + '\n')
        errfile.write("Reverse read:\n\n" + r_id + '\n')
        errfile.write(rseq + "\n\n" + ('~' * 40 + '\n') * 2 + '\n')
        errfile.write("ALIGNMENT (FORWARD READ AGAINST REVERSE ONE):\n\n")
    # end with

    # ===  Provide man who uses the program with information about how erroneous reads align one against another  ===

    erralign = SW_align(fseq, rseq)
    with open(error_report, 'a') as errfile:
        errfile.write("Q_START={}; Q_END={}; S_START={}; S_END={}; PIDENT={}\n".format( erralign.q_start, erralign.q_end,
                                                            erralign.s_start, erralign.s_end, erral.get_pident() ))
    # end with
# end def _handle_unforseen_case


def _del_temp_files(put_artif_dir):
    """
    Delete temporary files used by functions of these module.
    """
    if put_artif_dir is not None:
        for file in os.listdir(put_artif_dir):
            if os.path.exists(os.sep.join([put_artif_dir, file])) and "roughly" in file:
                os.remove(os.sep.join( (put_artif_dir, file) ))
            # end if
        # end for
    # end if
# end def _del_temp_files


def _naive_merging(fastq_recs):
    """
    The first of the "kernel" functions in this module. Performs rough process of read merging.
    :param fastq_reqs: a dictionary of two fastq-records stored as dictionary of it's fields
    :type fastq_reads: dict<str: dict<str, str>>

    The description of this weird parameter:
    The outer dictionary contains two dictionaries accessable by the following keys: "R1" (forward read), "R2" (reverse read).
    Fields of each record should be accessable by the follosing keys:
    1) "seq_id" (ID of the read)
    2) "seq" (sequence itsef)
    3) "opt_id" (the third line, where '+' is usually written)
    4) "qual_str" (quality string in Phred33)

    # Return values:
    # 0 -- if reads can be merged;
    # 1 -- if reads can't be merged;
    # 2 -- too short
    """

    f_id = fastq_recs["R1"]["seq_id"]
    fseq = fastq_recs["R1"]["seq"]
    fqual = fastq_recs["R1"]["qual_str"]
    r_id = fastq_recs["R2"]["seq_id"]
    rseq = _rc(fastq_recs["R2"]["seq"])         # reverse-complement
    rqual = fastq_recs["R2"]["qual_str"][::-1]    # reverse


    # |==== Firstly try to find overlapping region with naive-but-relatively-fast algorithm ====|

    # If reads are too short -- IndexError will be raised and these reads will be considered too short.
    # We do not need these reads anyway.
    try:
        shift = _naive_align(fseq, rseq)
    except IndexError:
        return 2
    # end try

    if shift != -1:
        overl = shift + _SEED_LEN
        loffset = len(fseq) - overl

        merged_seq, merged_qual = _merge_by_overlap(loffset, overl, fseq, fqual, rseq, rqual)
        
        globals()["_curr_merged_seq"] = merged_seq
        globals()["_curr_merged_qual"] = merged_qual
        return 0
    else:
        return 1
    # end if
# end def _naive_merging


#      |===== Smith-Waterman algorithm implementation =====|

class SWAlignLensNotEq(Exception):
    """
    An exception meant to be raised when length of the first aligned 
    sequence the second's one are not equal after aligning.
    Is subclass of Exception.
    """
    pass
# end class SWAlignLensNotEq


class SWInvalidMove(Exception):
    """
    An exception meant to be raised when Smith-Waterman algorithm can't perform the next move.
    Is subclass of Exception.
    """
    pass
# end class SWInvalidMove


class SWAlignResult:
    """
    Class SWAlignResult is dedicated to perform result of Smith-Waterman aligning.

    :field q_align: aligned first (aka query) sequence (e.i. with gaps represented as '-');
    :type q_align: str;
    :field s_align: aligned seconf (aka subject) sequence (e.i. with gaps represented as '-');
    :type s_align: str;
    :field q_start: 1-based number of query sequence, in which alignment starts;
    :type q_start: int;
    :field q_end: 1-based number of query sequence, in which alignment ends;
    :type q_end: int;
    :field s_start: 1-based number of subject sequence, in which alignment starts;
    :type s_start: int;
    :field s_end: 1-based number of subject sequence, in which alignment ends;
    :type s_end: int;

    :method __init__: accepts: (q_align, s_align, q_start, q_end, s_start, s_end);
        After calling the __init__ methof all these arguments will
        become corresponding fields of class's instance
    :method get_align_len: returns length of alignment length of 'int';
    :method get_pident: returns identity of the alignment of 'float' [0, 1];
    """

    def __init__(self, q_align, s_align, q_start, q_end, s_start, s_end):
        """
        :param q_align: aligned first (aka query) sequence (e.i. with gaps represented as '-');
        :type q_align: str;
        :param s_align: aligned seconf (aka subject) sequence (e.i. with gaps represented as '-');
        :type s_align: str;
        :param q_start: 1-based number of query sequence, in which alignment starts;
        :type q_start: int;
        :param q_end: 1-based number of query sequence, in which alignment ends;
        :type q_end: int;
        :param s_start: 1-based number of subject sequence, in which alignment starts;
        :type s_start: int;
        :param s_end: 1-based number of subject sequence, in which alignment ends;
        :type s_end: int;
        """
        if len(q_align) != len(s_align):
            raise SWAlignLensNotEq("""Smith-Waterman algorithm: lengths of alignments aren't equal:
        # end if
        query alignment length = {}; subject alignment length = {}""".format( len(q_align), len(s_align) ))
        

        self.q_align = q_align
        self.s_align = s_align
        self.q_start = q_start
        self.q_end = q_end
        self.s_start = s_start
        self.s_end = s_end
    # end def __init__


    def get_align_len(self):
        """
        Function returns length of the alignment of 'int'.
        Raises SWAlignLensNotEq if lengths of alignments aren't equal.
        """
        if len(self.q_align) == len(self.s_align):
            return len(self.q_align)
        
        else:
            raise SWAlignLensNotEq("""Smith-Waterman algorithm: lengths of alignments aren't equal:
        # end if
        query alignment length = {}; subject alignment length = {}""".format( len(self.q_align), len(self.s_align) ))
    # end def get_align_len


    def get_pident(self):
        """
        Function returns identity of the alignment of 'float' [0, 1].
        Raises SWAlignLensNotEq if lengths of alignments aren't equal.
        """
        if len(self.q_align) != len(self.s_align):
            raise SWAlignLensNotEq("""Smith-Waterman algorithm: lengths of alignments aren't equal:
        # end if
        query alignment length = {}; subject alignment length = {}""".format( len(q_align), len(s_align) ))
        

        pident = 0
        for i in range(len(self.q_align)):
            if self.q_align[i] == '-' or self.s_align[i] == '-':
                continue
            # end if
            pident += self.q_align[i] == self.s_align[i]
        # end for
        return round( pident / len(self.q_align), 2 )
    # end def get_pident


    def __repr__(self):
        return "\nQS={}; QE={}; SS={}; SE={}; PID={}".format( self.q_start, self.q_end,
                                                            self.s_start, self.s_end, self.get_pident() )
    # end def __repr__
# end class SWAlignResult


from functools import partial
from array import array


def SW_get_new_score(up, left, diag,
                  matched, gap_penalty, match, mismatch):
    """
    Function is dedicated to fill an element of alignment matrix during
    forward step of Smith-Waterman algorithm.

    :param up: value of the upper element of matrix;
    :type up: int;
    :param left: value of the left element of matrix;
    :type left: int;
    :param diag:value of the diagonal element of matrix;
    :type diagonal: int;
    :param matched: logic value: is True if nucleotides are equal and False otherwise;
    :type matched: bool;
    :param gap_penalty: penalty for gap;
    :type gap_penalty: int;
    :param match: reward for matching;
    :type match: int;
    :param mismatch: penalty for mismatch;
    :type mismatch: int;
    """

    add_to_diag = match if matched else mismatch
    return max(0, diag+add_to_diag, left-gap_penalty, up-gap_penalty)
# end def SW_get_new_score


# Constants for traceback.
END, DIAG, UP, LEFT = range(4)


def SW_next_move(SWsm, i, j, match, mismatch, gap_penalty, matched):
    """
    Function is dedicated to perform a move during traceback of Smith-Waterman algorithm.

    :param SWsm: Smith-Waterman scoring matrix;
    :type SWsm: list<list<int>>;   meant to be changed to list<array.array<unsigned_int>>
    :param i: row index of scoring matrix;
    :type i: int;
    :param j: column index of scoring matrix;
    :type j: int;
    :param gap_penalty: penalty for gap;
    :type gap_penalty: int;
    :param matched: logic value: is True if nucleotides are equal and False otherwise;
    :type matched: bool;

    Returns 1 if move is diagonal, 2 if move is upper, 3 if move is left.

    Raises SWInvalidMove when Smith-Waterman algorithm can't perform the next move.
    """
    diag = SWsm[i - 1][j - 1]
    up   = SWsm[i - 1][  j  ]
    left = SWsm[  i  ][j - 1]

    if (diag + (mismatch, match)[matched]) == SWsm[i][j]:
        return DIAG if diag!=0 else END
    
    elif (up - gap_penalty) == SWsm[i][j]:
        return UP if up!=0 else END
    
    elif (left - gap_penalty) == SWsm[i][j]:
        return LEFT if up!=0 else END
    
    else:
        # Execution should not reach here.
        raise SWInvalidMove("""Smith-Waterman algorithm: invalid move during traceback:
            diag={}; up={}; left={}; matched={}; SWsm[i][j]={}""".format(diag, up, left, matched, SWsm[i][j]))
    # end if
# end def SW_next_move


def SW_align(jseq, iseq, gap_penalty=5, match=1, mismatch=-1):
    """
    Function performs Smith-Waterman aligning algorithm.
    Reward and all penalties are int for the sake of acceleration of this slow snake.
    All this algorithm is the first candidate to be rewritten in C as Python extention.

    :param jseq: 'query' sequence, which is loated to the top of scoring matrix, across the j index varying;
    :type jseq: str;
    :param iseq: 'subject' sequence, which is loated to the left of scoring matrix, across the i index varying;
    :type iseq: str;
    :param gap_penalty: penalty for gap;
    :type gap_penalty: int;
    :param match: reward for matching;
    :type match: int;
    :param mismatch: penalty for mismatch;
    :type mismatch: int;

    Returns instance of SWAlignResult performing results of alignment.
    """

    SWsm = [array('I', [0 for j in range(len(jseq)+1)]) for i in range(len(iseq)+1)]
    # SWsm = [[0 for j in range(len(jseq)+1)] for i in range(len(iseq)+1)]
    score = partial(SW_get_new_score, gap_penalty=gap_penalty, match=match, mismatch=mismatch)
    max_val, max_pos = 0, None

    for i in range(1, len(iseq)+1):
        for j in range(1, len(jseq)+1):
            SWsm[i][j] = score(up=SWsm[i-1][j], left=SWsm[i][j-1], diag=SWsm[i-1][j-1], matched=iseq[i-1] == jseq[j-1])

            if SWsm[i][j] > max_val:
                max_val = SWsm[i][j]
                max_pos = (i, j)
            # end if
        # end for
    # end for

    # It happens (no alignment)
    if max_pos is None:
        return None
    # end if

    aligned_jseq = ""
    aligned_iseq = ""
    i, j = max_pos
    move = SW_next_move(SWsm, i, j, match, mismatch,
                        gap_penalty, iseq[i-1] == jseq[j-1])

    while move != END:

        if move == DIAG:
            aligned_jseq += jseq[j-1]
            aligned_iseq += iseq[i-1]
            j -= 1
            i -= 1

        elif move == UP:
            aligned_jseq += '-'
            aligned_iseq += iseq[i-1]
            i -= 1

        elif move == LEFT:
            aligned_jseq += jseq[j-1]
            aligned_iseq += '-'
            j -= 1
        # end if

        move = SW_next_move(SWsm, i, j, match, mismatch,
                            gap_penalty, iseq[i-1] == jseq[j-1])
    # end while

    if jseq[j-1] == iseq[i-1]:
        aligned_jseq += jseq[j-1]
        aligned_iseq += iseq[i-1]
    # end if

    return SWAlignResult(aligned_jseq[::-1], aligned_iseq[::-1], j, max_pos[1], i, max_pos[0])
# end def SW_align


def _accurate_merging(fastq_recs):
    """
    The second of the "kernel" functions in this module. Performs accurate process of read merging.
    :param fastq_reqs: a dictionary of two fastq-records stored as dictionary of it's fields
    :type fastq_reads: dict<str: dict<str, str>>
    Description of this weird parameter:
    The outer dictionary contains two dictionaries accessable by the following keys: "R1" (forward read), "R2" (reverse read).
    Fields of each record should be accessable by the follosing keys:
    1) "seq_id" (ID of the read)
    2) "seq" (sequence itsef)
    3) "opt_id" (the third line, where '+' is usually written)
    4) "qual_str" (quality string in Phred33)
    # Return values:
    # 0 -- if reads can be merged;
    # 1 -- if reads can't be merged;
    # 2 -- too short
    # 3 -- fatal error, unforseen case
    """

    f_id = fastq_recs["R1"]["seq_id"]
    fseq = fastq_recs["R1"]["seq"]
    fqual = fastq_recs["R1"]["qual_str"]
    r_id = fastq_recs["R2"]["seq_id"]
    rseq = _rc(fastq_recs["R2"]["seq"])         # reverse-complement
    rqual = fastq_recs["R2"]["qual_str"][::-1]      # reverse

    if len(fseq) != len(fqual) or len(rseq) != len(rqual):
        print("\n\nInvalid FASTQ format!\a")
        print("Lengths of the sequence and the quality line are unequal")
        print("ID if this read: '{}'".format(f_id))
        exit(1)
    # end if

    far_report = SW_align(fseq, rseq)

    if far_report is None:
        return 1
    # end if

    # |==== Check how they have aligned ====|
    
    # Consider normal overlap as:
    # FFFFFFFF----
    # ----RRRRRRRR
    reads_normally_overlap = True   # naive assumption

    # Discard following situation:
    # --FFFFFF
    # RRRRRR--
    too_short = far_report.q_start < far_report.s_start
    too_short = too_short or ( far_report.q_end / len(fseq) ) < ( far_report.s_end / len(rseq) )

    reads_normally_overlap = reads_normally_overlap and not too_short

    # Catch randomly occured alignment in center of sequences,
    # e. i. reads don't align against one another in a proper way
    rand_align = (len(fseq) - far_report.q_end) > _MAX_ALIGN_OFFSET or far_report.s_start > _MAX_ALIGN_OFFSET


    # resulting conslusion
    reads_normally_overlap = reads_normally_overlap and not rand_align

    # If there is not enough identity, it is better to make sure -- it is better to blast
    if reads_normally_overlap and far_report.get_pident() < _MIN_PIDENT:
        reads_normally_overlap = False
        rand_align = True
    # end if


    # |==== Decide what to do according to "analysis" above. ====|

    # === Ok, reads overlap in the way they should do it. ===
    # FFFFFFFF----
    # ----RRRRRRRR
    if reads_normally_overlap:

        loffset = far_report.q_start -  far_report.s_start # zero-based index of overlap start
        overl = len(fseq) - loffset # length of the overlap

        merged_seq, merged_qual = _merge_by_overlap(loffset, overl, fseq, fqual, rseq, rqual)
        globals()["_curr_merged_seq"] = merged_seq
        globals()["_curr_merged_qual"] = merged_qual

        return 0

    # === Randomly occured alignment in center of sequences, need to blast ===
    elif rand_align:

        # === Blast forward read, align reverse read against the reference ===
        # "faref" means Forward [read] Against REFerence [sequence]
        # "raref" means Reverse [read] Against REFerence [sequence]
        try:
            faref_report, raref_report = _blast_and_align(fseq, f_id, rseq, r_id)
        except NoRefAlignError:
            return 1
        # end try

        if float(faref_report[EVALUE]) > 1e-2:
            return 1
        # end if

        # Calculate some "features".
        # All these "features" are in coordinates of the reference sequence
        forw_start = int(faref_report[SSTART]) - int(faref_report[QSTART])
        forw_end = forw_start + len(fseq)
        rev_start = raref_report.s_start - raref_report.q_start
        rev_end = rev_start + len(rseq)

        gap = forw_end < rev_start

        if forw_end > rev_end:
            return 2
        # end if

        # ===  Handle alignment results ====

        # Overlap region is long enough
        if not gap and forw_end - rev_start > _MIN_OVERLAP:
            return 1

        # Length of the overlapping region is short,
        # but reverse read alignes as they should do it, so let it be.
        elif not gap and forw_end - rev_start <= _MIN_OVERLAP:
            
            loffset = rev_start - forw_start
            overl = len(fseq) - loffset
            merged_seq, merged_qual = _merge_by_overlap(loffset, overl, fseq, fqual, rseq, rqual)
            globals()["_curr_merged_seq"] = merged_seq
            globals()["_curr_merged_qual"] = merged_qual
            return 0

        # Here we have a gap.
        elif gap:
            
            gap_len = rev_start - forw_end
            # Very long gap isn't probable enough.
            if gap_len > _MAX_CRED_GAP_LEN:
                return 1
            
            else:
                # Fill in the gap with 'N'.
                merged_seq = fseq + 'N' * (gap_len - 1) + rseq
                merged_qual = fqual + chr(33) * (gap_len - 1) + rqual   # Illumina uses Phred33

                globals()["_curr_merged_seq"] = merged_seq
                globals()["_curr_merged_qual"] = merged_qual
                return 0
            # end if
        # end if

        # === Unforseen case. This code should not be ran. But if it does-- report about it. ===
        _handle_unforseen_case(f_id, fseq, r_id, rseq)
        return 3

    # Too short sequence, even if is merged. We do not need it.
    elif too_short:
        # --FFFFFFF
        # RRRRRRR--
        return 2
    # end if

    # === Unforseen case. This code should not be ran. But if it'll do -- report about it. ===
    _handle_unforseen_case(f_id, fseq, r_id, rseq)
    return 3
# end def _accurate_merging


def _handle_merge_pair_result(merging_result, fastq_recs, result_files, accurate=False):
    """
    This function handles the result of read merging and writes sequences in corresponding files.

    :param merging_result: a tuple of strings (if reads are merged), or an integer, if they cannot be merged
    :param fastq_reqs: a dictionary of two fastq-records stored as dictionary of it's fields
    :type fastq_reads: dict<str: dict<str, str>>
    :type result_files: dict<str: _io.TextIOWrapper>
    :type read_files: dict<str: _io.TextIOWrapper>
    """

    # if reads are merged
    if merging_result == 0:

        merged_rec = {
            "seq_id": fastq_recs["R1"]["seq_id"],
            "seq": _curr_merged_seq,
            "opt_id": fastq_recs["R1"]["opt_id"],
            "qual_str": _curr_merged_qual
        }
        _write_fastq_record(result_files["merg"], merged_rec)
        _merging_stats[merging_result] += 1
        if accurate:
            _merging_stats[1] -= 1
        # end if
        return 0
    
    # can't merge reads
    elif merging_result == 1:
        _write_fastq_record(result_files["umR1"], fastq_recs["R1"])
        _write_fastq_record(result_files["umR2"], fastq_recs["R2"])
        if not accurate:
            _merging_stats[merging_result] += 1
        # end if
        
        return 1
    
    # if resulting sequence is too short to distinguish taxa
    elif merging_result == 2:
        _write_fastq_record(result_files["shrtR1"], fastq_recs["R1"])
        _write_fastq_record(result_files["shrtR2"], fastq_recs["R2"])
        _merging_stats[merging_result] += 1
        if accurate:
            _merging_stats[1] -= 1
        # end if
        
        return 2
    
    # if unforseen situation occured in 'read_merging_16S'
    elif merging_result == 3:
        _close_files(result_files)
        input("Press ENTER to exit:")
        exit(1)
    
    # if 'read_merging_16S' returnes something unexpected and undesigned
    else:
        print("ERROR!!!\n\tModule that was merging reads returned an unexpected and undesigned value.")
        try:
            print("\tHere it is: {}".format(str(merging_result)))
        
        except:
            print("\tUnfortunately, it cannot be printed. It is of {} type".format(type(merging_result)))
        # end try
        
        print("Please, contact the developer.")
        _close_files(read_files, result_files)
        input("Press ENTER to exit:")
        exit(1)
    # end if
# end def _handle_merge_pair_result


def _check_file_existance(path):

    if os.path.exists(path):
        return
    else:
        print("\nFile '{}' does not exist!".format(path))
        exit(1)
    # end if
# end def _check_file_existance


def _open_files(read_paths, result_paths, wmode='w'):
    """
    Function opens read files and result files.

    'read_paths' and 'result_paths' are dictionaries (dict<str: str>).
    :param wmode: mode of 'open' function;
    :type wmode: str;

    Returns ruple of dictionaries with file objects.
    These dictionaries have the same structure as dictionaries passed to this function.
    """

    how_to_open = _OPEN_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]

    read_files = dict()
    result_files = dict()
    try:
        for key in read_paths.keys():
            read_files[key] = how_to_open(read_paths[key])
        # end for
        for key in result_paths.keys():
            result_files[key] = open(result_paths[key], wmode)
        # end for
    except OSError as oserror:
        print("Error while opening file", str(oserror))
        _close_files(read_files, result_files)
        input("Press enter to exit:")
        exit(1)
    # end try
    return read_files, result_files
# end def _open_files


def _fastq_read_packets(read_paths, reads_at_all, n_thr):
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

    how_to_open = _OPEN_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]
    fmt_func = _FORMATTING_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]
    read_files = dict() # dictionary for file objects

    # Open files that contain reads meant to be processed.
    for key, path in read_paths.items():
        read_files[key] = how_to_open(path)
    # end for

    # Compute packet size (one packet -- one thread).
    pack_size = reads_at_all // n_thr
    if reads_at_all % n_thr != 0: # in order not to import 'math'
        pack_size += 1
    # end if

    for i in range(n_thr):
        packet = list()
        for j in range(pack_size):
            tmp = _read_fastq_pair(read_files, fmt_func)
            # if the end of the line is reached
            if tmp is None:
                # Yield partial packet and return 'None' on next call
                yield packet
                # Python 2 throws a syntax error on the attempt of returning 'None' from the generator explicitly.
                # But single 'return' statement returns 'None' anyway, so here it is:
                return # None
            
            else:
                packet.append(tmp)
            # end if
        # end for
        yield packet # yield full packet
    # end for
# end def _fastq_read_packets


def _proc_init(print_lock_buff, counter_buff, count_lock_buff, write_lock_buff,
    result_files_buff,_merging_stats_buff):
    """
    Function that initializes global variables that all processes shoud have access to.
    This function is meant to be passed as 'initializer' argument to 'multiprocessing.Pool' function.

    :param print_lock_buff: lock that synchronizes printing to the console;
    :type print_lock_buff: multiprocessing.Lock;
    :param counter_buff: integer number representing number of processed reads;
    :type counter_buff: multiprocessing.Value;
    :param count_lock_buff: lock that synchronizes incrementing 'counter' variable;
    :type count_lock_buff: multiprocessing.Lock;
    :param write_lock_buff: lock that synchronizes writing to result files;
    :type write_lock_buff: multiprocessing.Lock;
    :param result_paths_buff: dictionary of paths to result files;
    :type result_paths_buff: dict<str: str>;
    :param _merging_stats_buff: array of values in which merging result statistics will be stored;
    :type _merging_stats_buff: multiprocessing.Array<int>;
    """

    global print_lock
    print_lock = print_lock_buff

    global counter
    counter = counter_buff

    global count_lock
    count_lock = count_lock_buff

    global write_lock
    write_lock = write_lock_buff

    global result_files
    result_files = result_files_buff

    global _merging_stats
    _merging_stats = _merging_stats_buff
# end def _proc_init


def _single_merger(merging_function, data, reads_at_all, accurate, delay=5, max_unwr_size=10):
    """
    Function that performs task meant to be done by one process while parallel accurate read merging.

    :param data: list of FASTQ-records
        (structure of these records is described in 'write_fastq_record' function);
    :type data: list< dict<str: str> >;
    :param reads_at_all: total number of read pairs in input files;
    :type reads_at_all: int;
    :param accurate: flag that is True if accurate merging is performing;
    :type accurate: bool;
    :param delay: number of reads that will be processed silently.
        I.e. status bar will be updated every 'delay' reads;
    :type delay: int;
    :param max_unwr_size: procesed reads will be written to result files every time next 'max_unwr_size'
        reads are processed;
    :type max_unwr_size: int;
    """

    # Processes will print number of processed reads every 'delay' reads.
    i =  0
    bar_len = 50 # length of status bar

    # Processes will write processed reads every 'max_unwr_size' reads.
    j = 0
    tmp_fq_recs = list()
    tmp_merge_res_list = list()
    
    for fastq_recs in data:

        tmp_merge_res_list.append(merging_function(fastq_recs))
        tmp_fq_recs.append(fastq_recs)
        j += 1

        if j == max_unwr_size:

            with write_lock:
                for k in range(max_unwr_size):
                    _handle_merge_pair_result(tmp_merge_res_list[k], tmp_fq_recs[k], result_files, accurate=accurate)
                # end for
            # end with

            j = 0
            tmp_fq_recs.clear()
            tmp_merge_res_list.clear()
        # end if

        if i == delay:
            # synchronized incrementing
            with count_lock:
                counter.value += i
            # end with
            
            i = 0

            eqs = int( bar_len*(counter.value / reads_at_all) ) # number of '=' characters
            spaces = bar_len - eqs # number of ' ' characters
            with print_lock:
                printn("\r[" + "="*eqs + '>' + ' '*spaces +"] {}% ({}/{})".format(int(counter.value/reads_at_all*100),
                    counter.value, reads_at_all))
            # end with
        # end if
        i += 1
    # end for

    # A little 'tail' of reads often will remain -- we do not want to lose them:
    with write_lock:
        for k in range(len(tmp_merge_res_list)):
            _handle_merge_pair_result(tmp_merge_res_list[k], tmp_fq_recs[k], result_files, accurate=accurate)
        # end for
    # end with
# end def _single_merger


def _parallel_merging(merging_function, read_paths, result_paths, n_thr, accurate, delay=5, max_unwr_size=10):
    """
    Function launches parallel read merging.

    :param merging_function: function that will be applied to reads;
    :param read_paths: dict of paths to read files. it's structure is described in '_fastq_read_packets' function;
    :type read_paths: dict<str: str>;
    :param result_paths: dict of paths to read files;
    :type result_paths: dict<str: str>;
    :param n_thr: int;
    :type n_thr: int:
    :param accurate: flag that is True if accurate merging is performing;
    :type accurate: bool;
    :param delay: number of reads that will be processed silently.
        I.e. status bar will be updated every 'delay' reads;
    :type delay: int;
    :param max_unwr_size: procesed reads will be written to result files every time next 'max_unwr_size'
        reads are processed;
    :type max_unwr_size: int;
    """
    
    print_lock = mp.Lock() # lock that synchronizes printing to the console;
    write_lock = mp.Lock() # lock that synchronizes writing to result files;
    counter = mp.Value('i', 0) # integer number representing number of processed reads;
    count_lock = mp.Lock() # lock that synchronizes 'counter' variable incrementing;
    sync_merg_stats = mp.Array('i', 3)
    for i in range(len(_merging_stats.values())):
        sync_merg_stats[i] = _merging_stats[i]
    # end for

    how_to_open = _OPEN_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]

    reads_at_all = int( sum(1 for line in how_to_open(read_paths["R1"])) / 4 )

    # Open result files
    result_files = dict()
    for key, path in result_paths.items():
        result_files[key] = open(path, 'a')
    # end for

    pool = mp.Pool(n_thr, initializer=_proc_init, initargs=(print_lock, counter, count_lock, write_lock, result_files, sync_merg_stats))
    pool.starmap(_single_merger, [(merging_function, data, reads_at_all, accurate, delay, max_unwr_size)
        for data in _fastq_read_packets(read_paths, reads_at_all, n_thr)])
        # end for

    _close_files(result_files)

    globals()["_merging_stats"] = {
        0: sync_merg_stats[0],
        1: sync_merg_stats[1],
        2: sync_merg_stats[2]
    }

    print("\r["+"="*50+"] 100% ({}/{})\n".format(reads_at_all, reads_at_all))
# end def _parallel_merging


def _one_thread_merging(merging_function, read_paths, wmode, result_paths, accurate, delay=1):
    """
    Function launches one-thread merging.
    
    :param merging_function: function that will be applied to reads;
    :param read_paths: dict of paths to read files. it's structure is described in '_fastq_read_packets' function;
    :type read_paths: dict<str: str>;
    :param wmode: mode of 'open' function;
    :type wmode: str;
    :param result_paths: dict of paths to read files;
    :type result_paths: dict<str: str>;
    :param accurate: flag that is True if accurate merging is performing;
    :type accurate: bool;
    :param delay: number of reads that will be processed silently.
        I.e. status bar will be updated every 'delay' reads;
    :type delay: int;
    """

    # Collect some info
    how_to_open = _OPEN_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]
    actual_format_func = _FORMATTING_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]
    readfile_length = sum(1 for line in how_to_open(read_paths["R1"]))  # lengths of read files are equal
    read_pairs_num = int(readfile_length / 4) # divizion by 4, because there are 4 lines per one fastq-record

    # Open files
    read_files, result_files = _open_files(read_paths, result_paths, wmode=wmode)

    # Proceed
    reads_processed, i = 0, 0
    bar_len = 50
    while reads_processed < read_pairs_num:

        fastq_recs = _read_fastq_pair(read_files, actual_format_func)

        merging_result = merging_function(fastq_recs)
        _handle_merge_pair_result(merging_result, fastq_recs, result_files, accurate=accurate)
        reads_processed += 1
        i += 1

        if i == delay:
            eqs = int( bar_len*(reads_processed / read_pairs_num) ) # number of '=' characters
            spaces = bar_len - eqs # number of ' ' characters
            printn("\r[" + "="*eqs + '>' + ' '*spaces +"] {}% ({}/{})".format(int(reads_processed/read_pairs_num*100),
                reads_processed, read_pairs_num))
            i = 0
        # end if
    # end while
    print("\r[" + "="*50 + "]" + "  100% ({}/{})\n\n"
        .format(reads_processed, read_pairs_num))
    _close_files(read_files, result_files)
# end def _one_thread_merging


# ===============================  "Public" stuff  ===============================


class AccessStatsBeforeMergingError(Exception):
    "Exception raising on the attempt of accessing merging statisttics before merging itself"
    pass
# end class AccessStatsBeforeMergingError


def get_merging_stats():
    """
    Function returns a dict<int: int> of the following format:
    {
        0: merged_read_pairs_number,
        1: unmerged_read_pairs_number,
        2: too_short_reads_number
    }
    Function rises an AccessStatsBeforeMergingError on the attempt of 
        accessing merging statisttics before merging (e.i. before calling 'merge_reads' function).
    """
    if _merging_stats is not None:
        return _merging_stats
    else:
        raise AccessStatsBeforeMergingError("You cannot access merging statistics before merging is completed!\a")
    # end if
# end def get_merging_stats


def merge_reads(R1_path, R2_path,
    outdir_path="read_merging_result_{}".format(strftime("%d_%m_%Y_%H_%M_%S", localtime(start_time))).replace(' ', '_'), n_thr=1):
    """
    This is the function that you should actually call from the outer scope in order to merge reads
    (and 'get_merging_stats' after it, if you want).

    :param R1_path: path to a file with forward reads
    :type R1_path: str
    :param R2_path: path to a file with reverse reads
    :type R2_path: str
    :param outdir_path: path to a directory, in which result files will be stored
    :type outdir_path: str
    :param n_thr: number of execution threads;
    :type n_thr: int;

    Function returns a dict<str: str> of the following format:
    {   
        "merg": path to a file with successfully merged reads,
        "umR1": path to a file with forward unmerged reads,
        "umR2": path to a file with reverse unmerged reads,
        "shrtR1": path to a file with forward reads considered as too short,
        "shrtR2": path to a file with reverse reads considered as too short
    }
    """

    # Create output directory
    if not os.path.exists(outdir_path):
        try:
            os.mkdir(outdir_path)
        
        except OSError as oserror:
            print("Error while creating result directory\n", str(oserror))
            input("Press ENTER to exit:")
            exit(1)
        # end try
    # end if

    # Create a directory for putative artifacts
    artif_dir = "{}{}putative_artifacts".format(outdir_path, os.sep)
    if not os.path.exists(artif_dir):
        try:
            os.mkdir(artif_dir)
        
        except OSError as oserror:
            print("Error while creating result directory\n", str(oserror))
            input("Press ENTER to exit:")
            exit(1)
        # end try
    # end if

    # The followig dictionary represents some statistics of merging.
    # Keys in this dictionary are in consistency with return codes of the function "_merge_pair"
    globals()["_merging_stats"] = {
        0: 0,           # number of merged reads
        1: 0,           # number of unmerged reads
        2: 0            # number of too short sequences
    }

    # |==== Naive merging ===|

    read_paths = {
        "R1": R1_path,
        "R2": R2_path
    }

    names = dict()
    for key, path in read_paths.items():
        file_name_itself = os.path.basename(read_paths[key])
        names[key] = search(r"(.*)\.f(ast)?q", file_name_itself).group(1)
    # end for
    
    # Name without "__R1__" and "__R2__":
    more_common_name = names["R1"][: names["R1"].rfind("R1")].strip('_')

    result_paths = {
        # File for merged sequences:
        "merg": "{}{}{}.merged.fastq".format(outdir_path, os.sep, more_common_name),
        # Files for unmerged sequences:
        "umR1": "{}{}{}.roughly_unmerged.fastq".format(artif_dir, os.sep, names["R1"]),
        "umR2": "{}{}{}.roughly_unmerged.fastq".format(artif_dir, os.sep, names["R2"]),
        # Files for those reads, that form too short sequence after merging:
        "shrtR1": "{}{}{}.too_short.fastq".format(artif_dir, os.sep, names["R1"]),
        "shrtR2": "{}{}{}.too_short.fastq".format(artif_dir, os.sep, names["R2"])
    }

    print("\n{} - Read merging started".format(get_work_time()))
    print("\nProceeding...\n\n")
    printn("[" + " "*50 + "]" + "  0%")

    if n_thr == 1:
        _one_thread_merging(_naive_merging, read_paths, 'w', result_paths,
            accurate=False, delay=100)
    else:
        _parallel_merging(_naive_merging, read_paths, result_paths, n_thr,
            accurate=False, delay=100, max_unwr_size=100)
    # end if

    print("\n{} - 1-st step of read merging is completed".format(get_work_time()))
    print("""\n{} read pairs have been merged together
{} read pairs haven't been merged together
{} read pairs have been considered as too short"""
.format(_merging_stats[0], _merging_stats[1], _merging_stats[2]))
    print('\n' + '~' * 50 + '\n')
    
    _del_temp_files(None)


    # |==== Accurate merging ===|

    read_paths = {
        "R1": result_paths["umR1"],
        "R2": result_paths["umR2"]
    }

    result_paths["umR1"] = result_paths["umR1"].replace("roughly_", "")
    result_paths["umR2"] = result_paths["umR2"].replace("roughly_", "")

    print("""\nNow the program will merge the rest of reads
    with accurate-but-slow strategy""")
    # end with
    print("\tIt will take a while")
    print("\n{} - Proceeding...\n\n".format(get_work_time()))
    printn("[" + " "*50 + "]" + "  0%")

    if n_thr == 1:
        _one_thread_merging(_accurate_merging, read_paths, 'a', result_paths, accurate=True)
    else:
        _parallel_merging(_accurate_merging, read_paths, result_paths, n_thr, accurate=True, delay=n_thr)
    # end if

    print("\n{} - Read merging is completed".format(get_work_time()))
    print("\nFinally,")
    print("""  {} read pairs have been merged together
  {} read pairs haven't been merged together
  {} read pairs have been considered as too short"""
  .format(_merging_stats[0], _merging_stats[1], _merging_stats[2]))
    print('\n' + '~' * 50 + '\n')
    _del_temp_files(artif_dir)

    return result_paths
# end def merge_reads


# |==== 'read_merging_16S.py' can be used as script ====|

if __name__ == "__main__":

    from sys import version_info as verinf
    if verinf.major < 3:
        print( "\nYour python interpreter version is " + "%d.%d" % (verinf.major, verinf.minor) )
        print("\tPlease, use Python 3!\a")
        # In python 2 'raw_input' does the same thing as 'input' in python 3.
        # Neither does 'input' in python 2.
        raw_input("Press ENTER to exit:")
        exit(1)
    # end if

    import getopt
    from sys import argv

    usage_msg = """
Module "read_merging_16S" is designed for merging Illumina (MiSeq) pair-end reads from 16S-rDNA.

Attention! This script cannot be executed by python interpreter version < 3.0!

It can be used as script as well as be imorted as module from outer Python program
    and then used via calling 'merge_reads' function.

Usage:
    ./read_merging_16S.py -1 forward_reads -2 reverse_reads [-o output_dir]
Options:
    -h (--help) --- show help message;\n
    -v (--version) --- show version;\n
    -1 (--R1) --- FASTQ file, in which forward reads are stored;\n
    -2 (--R2) --- FASTQ file, in which reverse reads are stored;\n
    -o (--outdir) --- directory, in which result files will be placed;\n
    -t (--threads) <int> --- number of threads to launch;
"""

    try:
        opts, args = getopt.getopt(argv[1:], "hv1:2:o:t:", ["help", "version", "R1=", "R2=", "outdir=", "threads="])
    except getopt.GetoptError as opt_err:
        print(str(opt_err) + '\a')
        print("See help ('-h' option)")
        exit(2)
    # end try

    outdir_path = "{}{}read_merging_16S_result_{}".format(os.getcwd(), os.sep,
        strftime("%d_%m_%Y_%H_%M_%S", localtime(start_time))).replace(" ", "_") # default path
    read_paths = dict()
    n_thr = 1

    # First search for information-providing options:

    if "-h" in argv[1:] or "--help" in argv[1:]:
        print(usage_msg)
        exit(0)
    # end if

    if "-v" in argv[1:] or "--version" in argv[1:]:
        print(__version__)
        exit(0)
    # end if

    for opt, arg in opts:

        if opt in ("-o", "--outdir"):
            outdir_path = os.path.abspath(arg)

        elif opt in ("-1", "--R1"):
            if not "-2" in argv and not "--R2" in argv:
                print("\nATTENTION!\n\tYou should specify both forward and reverse reads!\a")
                exit(1)
            # end if
            _check_file_existance(arg)
            read_paths["R1"] = arg

        elif opt in ("-2", "--R2"):
            if not "-1" in argv and not "--R1" in argv:
                print("\nATTENTION!\n\tYou should specify both forward and reverse reads!\a")
                exit(1)
            # end if
            _check_file_existance(arg)
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
            
            if n_thr > mp.cpu_count():
                print("""\n  Warning! You have specified {} threads to use
        although {} are available.\n""".format(n_thr, mp.cpu_count()))
                reply = input("""If this is just what you want, press ENTER
        or enter 'q' to exit:>>""")
                if reply == "":
                    pass
                else:
                    exit(0)
                # end if
            # end if
        # end if
    # end for

    # Check utilities for read merging
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
            print("\tAttention!\n{} is not found in your system.\a".format(utility))
            print("If you want to use this feature, please install {}".format(utility))
            print("""If this error still occure although you have installed everything 
    -- make sure that this program is added to PATH)""")
        # end if
    # end for

    print("\n |=== read_merging_16S.py (version {}) ===|\n".format(__version__))

    print( get_work_time() + " ({}) ".format(strftime("%d.%m.%Y %H:%M:%S", localtime(start_time))) + "- Start working\n")

    print("\nFollowing files will be processed:")
    for i, path in enumerate(read_paths.values()):
        print("  {}. '{}'".format(i+1, os.path.abspath(path)))
    # end for
    print() # just print blank line

    print("\nResult files will be placed in the following directory:\n\t'{}'\n".format(outdir_path))

    # Proceed
    result_files = merge_reads(read_paths["R1"], read_paths["R2"], outdir_path=outdir_path, n_thr=n_thr)

    print("{} read pairs have been processed.".format(_merging_stats[0]+_merging_stats[1]+_merging_stats[2]))

    # Remove empty files
    for file in result_files.values():
        if os.stat(file).st_size == 0:
            os.remove(file)
            print("'{}' is removed since it is empty".format(file))
        # end if
    # end for

    # Gzip result files
    print("\nGzipping result files...")
    for file in result_files.values():
        if os.path.exists(file):
            with open(file, 'r') as plain_file, open_as_gzip(file+".gz", 'wb') as gz_file:
                for line in plain_file:
                    gz_file.write(bytes(line, "utf-8"))
                # end for
            # end with
            os.unlink(file)
            print("\'{}\' is gzipped".format(file))
        # end if
    # end for
    
    print("Gzipping is completed\n")
    print("Result files are placed in the following directory:\n\t'{}'\n".format(outdir_path))

    # Write log file
    with open("{}{}read_merging_16S_{}.log".format(outdir_path, os.sep,
            strftime("%d_%m_%Y_%H_%M_%S", localtime(start_time))).replace(" ", "_"), 'w') as logfile:
        logfile.write("Script 'read_merging_16S.py' was ran at {}\n".format(strftime("%d.%m.%Y.%H:%M:%S", localtime(start_time))))
        end_time = strftime("%d_%m_%Y_%H_%M_%S", localtime(time()))

        logfile.write("\nFollowing files have been processed:\n")
        logfile.write("  '{}'\n  '{}'\n\n".format(os.path.abspath(read_paths["R1"]), os.path.abspath(read_paths["R2"])))

        logfile.write("{} read pairs have been processed.\n".format(_merging_stats[0]+_merging_stats[1]+_merging_stats[2]))

        logfile.write("Script completed it's job at {}\n\n".format(end_time))
        logfile.write("{} read pairs have been merged.\n".format(_merging_stats[0]))
        logfile.write("{} read pairs haven't been merged together.\n".format(_merging_stats[1]))
        logfile.write("{} read pairs have been considered as too short for merging.\n".format(_merging_stats[2]))
    # end with

    print('\n'+get_work_time() + " ({}) ".format(strftime("%d.%m.%Y %H:%M:%S", localtime(time()))) + "- Job is successfully completed!\n")
# end if
