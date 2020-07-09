#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "4.0.a"
# Year, month, day
# __last_update_date__ = "2020-08-07"

import os
import re
import sys
import math

from gzip import open as open_as_gzip
from bz2 import open as open_as_bz2

from subprocess import Popen as sp_Popen, PIPE as sp_PIPE
import multiprocessing as mp

from functools import partial
from array import array

from src.printing import *
from src.fastq import *
from src.filesystem import *

_ncbi_fmt_db = "/mnt/1TB.Toshiba/sikol_tmp/SILVA_DB/SILVA_138_SSURef_NR99_tax_silva_trunc.fasta"

class SilvaDBNotInstalledError(Exception):
    pass
# end class SilvaDBNotInstalledError


if not os.path.exists(_ncbi_fmt_db + ".nhr"):
    raise SilvaDBNotInstalledError("Silva database is not installed!\n  Please, run 'install_read_merging_16S.sh'")
# end if


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
_MAX_ALIGN_OFFSET = 5
# --------------------------------------------

# According to https://support.illumina.com/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf
_INSERT_LEN = 550


# The following dictionary represents some statistics of merging.
# Keys in this dictionary are in consistency with return codes of the function "_merge_pair"
_merging_stats = None    # it in None in the beginning, because there is no statistics before merging


# ===============================  Internal functions  ===============================


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
        print("ID if erroneous read:\n  '{}'".format(f_id))
        sys.exit(exit_code)
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
        sys.exit(exit_code)
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


def _del_temp_files(put_artif_dir=None):
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
        query alignment length = {}; subject alignment length = {}""".format( len(q_align), len(s_align) ))
        # end if
        

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
        query alignment length = {}; subject alignment length = {}""".format( len(self.q_align), len(self.s_align) ))
        # end if
    # end def get_align_len


    def get_pident(self):
        """
        Function returns identity of the alignment of 'float' [0, 1].
        Raises SWAlignLensNotEq if lengths of alignments aren't equal.
        """
        if len(self.q_align) != len(self.s_align):
            raise SWAlignLensNotEq("""Smith-Waterman algorithm: lengths of alignments aren't equal:
        query alignment length = {}; subject alignment length = {}""".format( len(q_align), len(s_align) ))
        # end if
        

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


def _accurate_merging(fastq_recs, phred_offset, num_N, min_overlap, mismatch_frac):
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
    4) "qual_str" (quality string)
    # Return values:
    # 0 -- if reads can be merged;
    # 1 -- if reads can't be merged;
    # 2 -- dovetailed
    # 3 -- fatal error, unforseen case
    """

    f_id = fastq_recs["R1"]["seq_id"]
    fseq = fastq_recs["R1"]["seq"]
    fqual = fastq_recs["R1"]["qual_str"]
    r_id = fastq_recs["R2"]["seq_id"]
    rseq = _rc(fastq_recs["R2"]["seq"])         # reverse-complement
    rqual = fastq_recs["R2"]["qual_str"][::-1]      # reverse

    len_fseq = len(fseq)
    len_rseq = len(rseq)

    if len_fseq != len(fqual) or len_rseq != len(rqual):
        print("\n\nInvalid FASTQ format!\a")
        print("Lengths of the sequence and the quality line are unequal")
        print("ID if this read: '{}'".format(f_id))
        sys.exit(1)
    # end if

    if len_fseq + len_rseq + num_N < _INSERT_LEN:
        return 1, None
    # end if

    # === Blast forward read, align reverse read against the reference ===
    # "faref" means Forward [read] Against REFerence [sequence]
    # "raref" means Reverse [read] Against REFerence [sequence]
    try:
        faref_report, raref_report = _blast_and_align(fseq, f_id, rseq, r_id)
    except NoRefAlignError:
        return 1, None
    # end try

    if float(faref_report[EVALUE]) > 1e-2:
        return 1, None
    # end if

    # Calculate some "features".
    # All these "features" are in coordinates of the reference sequence
    forw_start = int(faref_report[SSTART]) - int(faref_report[QSTART])
    forw_end = forw_start + len_fseq
    rev_start = raref_report.s_start - raref_report.q_start
    rev_end = rev_start + len_rseq

    gap = forw_end < rev_start

    if forw_end > rev_end:
        # Dovetailed read pair
        return 1, None
    # end if

    # ===  Handle alignment results ====

    # Overlap region is long enough
    if not gap and forw_end - rev_start > min_overlap:
        return 1, None

    # Length of the overlapping region is short,
    # but reverse read alignes as they should do it, so let it be.
    elif not gap and forw_end - rev_start <= min_overlap:
        
        loffset = rev_start - forw_start
        overl = len_fseq - loffset
        merged_seq, merged_qual = _merge_by_overlap(loffset, overl, fseq, fqual, rseq, rqual)

        return 0, {"seq": merged_seq, "qual_str": merged_qual}

    # Here we have a gap.
    elif gap:
        
        gap_len = rev_start - forw_end
        # If gap is too long -- discard this pair.
        if gap_len > num_N:
            return 1, None
        
        else:
            # Fill in the gap with 'N'.
            merged_seq = fseq + 'N' * (gap_len - 1) + rseq
            merged_qual = fqual + chr(phred_offset+3) * (gap_len - 1) + rqual   # Illumina uses Phred33

            return 0, {"seq": merged_seq, "qual_str": merged_qual}
        # end if
    # end if

    # === Unforseen case. This code should not be ran. But if it does-- report about it. ===
    _handle_unforseen_case(f_id, fseq, r_id, rseq)
    return 3, None
# end def _accurate_merging


def _handle_merge_pair_result(merging_result, fastq_recs, result_files, merged_strs, accurate=False):
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

        if merged_strs is None:
            print("Fatal error 77. Please, contact the developer.")
            sys.exit(77)
        # end if

        merged_rec = {
            "seq_id": fastq_recs["R1"]["seq_id"],
            "seq": merged_strs["seq"],
            "opt_id": fastq_recs["R1"]["opt_id"],
            "qual_str": merged_strs["qual_str"]
        }
        write_fastq_record(result_files["merg"], merged_rec)
        _merging_stats[merging_result] += 1
        if accurate:
            _merging_stats[1] -= 1
        # end if
        return 0
    
    # can't merge reads
    elif merging_result == 1:
        write_fastq_record(result_files["umR1"], fastq_recs["R1"])
        write_fastq_record(result_files["umR2"], fastq_recs["R2"])
        if not accurate:
            _merging_stats[merging_result] += 1
        # end if
        
        return 1
    
    # if unforseen situation occured in 'read_merging_16S'
    elif merging_result == 3:
        close_files(result_files)
        input("Press ENTER to exit:")
        sys.exit(1)
    
    # if 'read_merging_16S' returnes something unexpected and undesigned
    else:
        print("ERROR!!!\n\tModule that was merging reads returned an unexpected and undesigned value.")
        print("\tHere it is: {}".format(str(merging_result)))
        
        print("Please, contact the developer.")
        close_files(read_files, result_files)
        input("Press ENTER to exit:")
        sys.exit(1)
    # end if
# end def _handle_merge_pair_result


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

    how_to_open = OPEN_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]
    fmt_func = FORMATTING_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]
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
            tmp = read_fastq_pair(read_files, fmt_func)
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


def _one_thread_merging(merging_function, read_paths, wmode,
    result_paths, accurate, num_N, min_overlap, mismatch_frac, delay=1, phred_offset=33):
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
    how_to_open = OPEN_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]
    actual_format_func = FORMATTING_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]
    readfile_length = sum(1 for line in how_to_open(read_paths["R1"]))  # lengths of read files are equal
    read_pairs_num = int(readfile_length / 4) # divizion by 4, because there are 4 lines per one fastq-record

    # Open files
    read_files = open_files(read_paths, how_to_open)
    result_files = open_files(result_paths, how_to_open, wmode)

    # Proceed
    reads_processed, i = 0, 0
    readnum_digits = math.ceil(math.log(read_pairs_num, 10))
    bar_len = min(50, os.get_terminal_size().columns - (11+2*readnum_digits))
    while reads_processed < read_pairs_num:

        fastq_recs = read_fastq_pair(read_files, actual_format_func)

        merging_result, merged_strs = merging_function(fastq_recs, phred_offset, num_N, min_overlap, mismatch_frac)
        _handle_merge_pair_result(merging_result, fastq_recs, result_files, merged_strs, accurate=accurate)
        reads_processed += 1
        i += 1

        if i == delay:
            bar_len = min(50, os.get_terminal_size().columns - (11+2*readnum_digits))
            eqs = int( bar_len*(reads_processed / read_pairs_num) ) # number of '=' characters
            spaces = bar_len - eqs # number of ' ' characters
            printn("\r[" + "="*eqs + '>' + ' '*spaces +"] {}% ({}/{})".format(int(reads_processed/read_pairs_num*100),
                reads_processed, read_pairs_num))
            i = 0
        # end if
    # end while
    print("\r[" + "="*bar_len + "]" + "  100% ({}/{})\n\n"
        .format(reads_processed, read_pairs_num))
    close_files(read_files, result_files)
# end def _one_thread_merging


def _parallel_merging(merging_function, read_paths, result_paths, n_thr,
    accurate, num_N, min_overlap, mismatch_frac, delay=5, max_unwr_size=10, phred_offset=33):
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

    how_to_open = OPEN_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]

    reads_at_all = int( sum(1 for line in how_to_open(read_paths["R1"])) / 4 )

    # Open result files
    result_files = dict()
    for key, path in result_paths.items():
        result_files[key] = open(path, 'a')
    # end for

    pool = mp.Pool(n_thr, initializer=_proc_init,
        initargs=(print_lock, counter, count_lock, write_lock, result_files, sync_merg_stats))
    pool.starmap(_single_merger,
        [(merging_function, data, reads_at_all, accurate, delay, max_unwr_size, phred_offset,
            num_N, min_overlap, mismatch_frac)
        for data in _fastq_read_packets(read_paths, reads_at_all, n_thr)])

    # Reaping zombies
    pool.close()
    pool.join()

    close_files(result_files)

    globals()["_merging_stats"] = {
        0: sync_merg_stats[0],
        1: sync_merg_stats[1]
    }

    print("\r["+"="*50+"] 100% ({}/{})\n".format(reads_at_all, reads_at_all))
# end def _parallel_merging



def _single_merger(merging_function, data, reads_at_all, accurate,
    delay, max_unwr_size, phred_offset,
    num_N, min_overlap, mismatch_frac):
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
    readnum_digits = math.ceil(math.log(reads_at_all, 10))
    bar_len = min(50, os.get_terminal_size().columns - (11+2*readnum_digits)) # length of status bar

    # Processes will write processed reads every 'max_unwr_size' reads.
    j = 0
    tmp_fq_recs = list()
    tmp_merge_res_list = list()

    EXT_CODE, SEQS = range(2)
    
    for fastq_recs in data:

        tmp_merge_res_list.append(merging_function(fastq_recs, phred_offset, num_N, min_overlap, mismatch_frac))
        tmp_fq_recs.append(fastq_recs)
        j += 1

        if j == max_unwr_size:

            with write_lock:
                for k in range(max_unwr_size):
                    _handle_merge_pair_result(tmp_merge_res_list[k][EXT_CODE], tmp_fq_recs[k],
                        result_files, tmp_merge_res_list[k][SEQS], accurate=accurate)
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

            bar_len = min(50, os.get_terminal_size().columns - (11+2*readnum_digits)) # length of status bar
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
            _handle_merge_pair_result(tmp_merge_res_list[k][EXT_CODE], tmp_fq_recs[k],
                result_files, tmp_merge_res_list[k][SEQS], accurate=accurate)
        # end for
    # end with
# end def _single_merger


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
        1: unmerged_read_pairs_number
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


def merge_reads(R1_path, R2_path, ngmerge,
    outdir_path="read_merging_result_{}".format(strftime("%d_%m_%Y_%H_%M_%S", localtime(start_time))).replace(' ', '_'),
    n_thr=1, phred_offset=33, num_N=35, min_overlap=20, mismatch_frac=0.1):
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
        "umR2": path to a file with reverse unmerged reads
    }
    """

    # Create a directory for putative artifacts
    artif_dir = "{}{}putative_artifacts".format(outdir_path, os.sep)
    if not os.path.exists(artif_dir):
        try:
            os.mkdir(artif_dir)
        
        except OSError as oserror:
            print("Error while creating result directory\n", str(oserror))
            sys.exit(1)
        # end try
    # end if

    # The followig dictionary represents some statistics of merging.
    # Keys in this dictionary are in consistency with return codes of the function "_merge_pair"
    globals()["_merging_stats"] = {
        0: 0,           # number of merged reads
        1: 0           # number of unmerged reads
    }

    # |==== Run NGmerge ===|

    read_paths = {
        "R1": R1_path,
        "R2": R2_path
    }

    names = dict()
    for key, path in read_paths.items():
        file_name_itself = os.path.basename(read_paths[key])
        names[key] = re.search(r"(.*)\.f(ast)?q", file_name_itself).group(1)
    # end for
    
    # Name without "__R1__" and "__R2__":
    more_common_name = names["R1"][: names["R1"].rfind("R1")].strip('_')

    result_paths = {
        # File for merged sequences:
        "merg": "{}{}{}.merged.fastq".format(outdir_path, os.sep, more_common_name),
        # Files for unmerged sequences:
        "umR1": "{}{}{}.unmerged.fastq".format(artif_dir, os.sep, names["R1"]),
        "umR2": "{}{}{}.unmerged.fastq".format(artif_dir, os.sep, names["R2"])
    }

    # Run NGmerge
    print("\n{} - Read merging started".format(get_work_time()))
    print("\nRunning NGmerge...\n")

    # NGmerge puts result files into working directory -- we will temporary go to output directory
    old_dir = os.getcwd()
    os.chdir(outdir_path)

    unmerged_prefix = "{}.unmerged".format(more_common_name)

    ngmerge_cmd = "{} -1 {} -2 {} -o {} -f {} -n {} -v -m {} -p {}".format(ngmerge, read_paths["R1"], read_paths["R2"],
        "{}.merged.fastq".format(more_common_name), unmerged_prefix, n_thr, min_overlap, mismatch_frac)
    print(ngmerge_cmd + '\n')
    print("NGmerge is doing it's job silently...")
    pipe = sp_Popen(ngmerge_cmd, shell = True, stderr=sp_PIPE)
    stderr = pipe.communicate()[1].decode("utf-8")

    if pipe.returncode != 0:
        # error
        print_error("error running NGmerge.")
        print(stderr)
        sys.exit(pipe.returncode)
    else:
        # Parse merging statistoct from NGmerge's stderr
        stderr = stderr.splitlines()[1:]
        reads_pattern = r"Fragments \(pairs of reads\) analyzed: ([0-9]+)"
        merged_pattern = r"Successfully stitched: ([0-9]+)"
        try:
            reads_processed = int(re.search(reads_pattern, stderr[0]).group(1))
            merged_reads = int(re.search(merged_pattern, stderr[1]).group(1))
            globals()["_merging_stats"][0] = merged_reads
            globals()["_merging_stats"][1] = reads_processed - merged_reads
        except AttributeError:
            print_error("error 78")
            print("Please, contact the developer.")
            sys.exit(78)
        # end try
    # end if

    # Move files with unmerged reads to directory 'putative_artifacts'
    os.rename("{}_1.fastq".format(unmerged_prefix), result_paths["umR1"])
    os.rename("{}_2.fastq".format(unmerged_prefix), result_paths["umR2"])
    os.chdir(old_dir) # returs to old dir

    print("\n{} - 1-st step of read merging is completed".format(get_work_time()))
    print("""\n{} read pairs have been merged together
{} read pairs haven't been merged together"""
.format(_merging_stats[0], _merging_stats[1]))
    print('\n' + '~' * 50 + '\n')

    _del_temp_files()


    # |==== Accurate merging ===|

    read_paths = {
        "R1": result_paths["umR1"],
        "R2": result_paths["umR2"]
    }

    print("""\nNow the program will merge the rest of reads
  filling gaps between them.""")
    # end with
    print("\tIt will take a while")
    print("\n{} - Proceeding...\n\n".format(get_work_time()))
    printn("[" + " "*50 + "]" + "  0%")

    if n_thr == 1:
        _one_thread_merging(_accurate_merging, read_paths, 'a', result_paths,
            True, num_N, min_overlap, mismatch_frac, phred_offset=phred_offset)
    else:
        _parallel_merging(_accurate_merging, read_paths, result_paths, n_thr,
            True, num_N, min_overlap, mismatch_frac, delay=n_thr, phred_offset=phred_offset)
    # end if

    print("\n{} - Read merging is completed".format(get_work_time()))
    print("\nFinally,")
    print("""  {} read pairs have been merged together
  {} read pairs haven't been merged together."""
  .format(_merging_stats[0], _merging_stats[1]))
    print('\n' + '~' * 50 + '\n')
    _del_temp_files(artif_dir)

    return result_paths
# end def merge_reads


# |==== 'read_merging_16S.py' can be used as script ====|

if __name__ == "__main__":

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

    import getopt

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
    -f (--phred-offset) [33, 64] --- Phred quality offset (default -- 33);
"""

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hv1:2:o:t:f:", ["help", "version", "R1=", "R2=",
            "outdir=", "threads=", "phred_offset=", "ngmerge-path=", "num-N=", "min-overlap=", "mismatch-frac="])
    except getopt.GetoptError as opt_err:
        print(str(opt_err) + '\a')
        print("See help ('-h' option)")
        sys.exit(2)
    # end try

    outdir_path = "{}{}read_merging_16S_result_{}".format(os.getcwd(), os.sep,
        strftime("%d_%m_%Y_%H_%M_%S", localtime(start_time))).replace(" ", "_") # default path
    read_paths = dict()
    n_thr = 1
    phred_offset = 33
    ngmerge = os.path.join(os.path.dirname(os.path.abspath(__file__)), "binaries", "NGmerge")
    # Length of (consensus?) conservative region between V3 and V4 in bacteria:
    # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0007401
    num_N = 35
    min_overlap = 20 # as default in NGmerge
    mismatch_frac = 0.1 # as default in NGmerge

    # First search for information-providing options:

    if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
        print(usage_msg)
        sys.exit(0)
    # end if

    if "-v" in sys.argv[1:] or "--version" in sys.argv[1:]:
        print(__version__)
        sys.exit(0)
    # end if

    for opt, arg in opts:

        if opt in ("-o", "--outdir"):
            outdir_path = os.path.abspath(arg)

        elif opt in ("-1", "--R1"):
            if not "-2" in sys.argv and not "--R2" in sys.argv:
                print("\nATTENTION!\n\tYou should specify both forward and reverse reads!\a")
                sys.exit(1)
            # end if
            if not os.path.exists(arg):
                print_error("File '{}' does not exist!".format(arg))
                sys.exit(1)
            # end if
            read_paths["R1"] = os.path.abspath(arg)

        elif opt in ("-2", "--R2"):
            if not "-1" in sys.argv and not "--R1" in sys.argv:
                print("\nATTENTION!\n\tYou should specify both forward and reverse reads!\a")
                sys.exit(1)
            # end if
            if not os.path.exists(arg):
                print_error("File '{}' does not exist!".format(arg))
                sys.exit(1)
            # end if
            read_paths["R2"] = os.path.abspath(arg)

        elif opt in ("-t", "--threads"):
            try:
                n_thr = int(arg)
                if n_thr < 1:
                    raise ValueError
                # end if
            except ValueError:
                print_error("number of threads must be positive integer number!")
                print(" And here is your value: '{}'".format(arg))
                sys.exit(1)
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
                sys.exit(1)
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
        # end if
    # end for

    # Check if NGmerge is executable
    if not os.access(ngmerge, os.X_OK):
        print_error("NGmerge file is not executable: '{}'".format(ngmerge))
        print("Please, make it executable. You can do it in this way:")
        print(" chmod +x {}".format(ngmerge))
    # end if

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

    # === Create output directory. ===
    if not os.path.exists(outdir_path):
        try:
            os.makedirs(outdir_path)
        except OSError as oserror:
            print_error("Error while creating result directory")
            print( str(oserror) )
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
                            print_error("error 77")
                            print("Please, contact the developer.")
                            sys.exit(77)
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

    # == Proceed ==
    result_files = merge_reads(read_paths["R1"], read_paths["R2"], ngmerge=ngmerge,
        outdir_path=outdir_path, n_thr=n_thr, phred_offset=phred_offset)

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
        logfile.write("{} read pairs have been considered as dovetailed.\n".format(_merging_stats[2]))
    # end with

    print('\n'+get_work_time() + " ({}) ".format(strftime("%d.%m.%Y %H:%M:%S", localtime(time()))) + "- Job is successfully completed!\n")
# end if
