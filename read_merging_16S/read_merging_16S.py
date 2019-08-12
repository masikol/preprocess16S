#!/usr/bin/env python3
"""
Module "read_merging_16S" is designed for merging Illumina (MiSeq) pair-end reads from 16S-rDNA.

Attention! This script cannot be executed by python interpreter version < 3.0!

It can be used as script as well as be imorted as module from outer Python program
    and then used via calling 'merge_reads' function.

Usage:
    python read_merging_16S.py -1 forward_reads -2 reverse_reads [-o output_dir]
Options:
    -h or --help
        output help message
    -1 or --R1
        file, in which forward reads are stored;
    -2 or --R2
        file, in which reverse reads are stored;
    --V3-V4
        by specifying this option, more accurate merging can be performed,
        but only if target sequences contain V3 and V4 regions and a constant region between them.
    -o or --outdir
        directory, in which result files will be placed.

Last modified 29.06.2019
"""

import os
from re import search
from re import match
from gzip import open as open_as_gzip


# |===== Colors =====|

def get_colored_print(color):#{

    def colored_print(text):#{
        print(color + text + "\033[0m")
    #}
    return colored_print
#}

YELLOW = "\033[1;33m"
RED = "\033[1;31m"
GREEN = "\033[1;32m"

print_red = get_colored_print(RED)
print_yellow = get_colored_print(YELLOW)
print_green = get_colored_print(GREEN)


# |===============================  Data   ===============================|

from datetime import datetime
_now = datetime.now().strftime("%Y-%m-%d %H.%M.%S")

_OPEN_FUNCS = (open, open_as_gzip)

# Data from .fastq and .fastq.gz should be parsed in different way,
#   because data from .gz is read as bytes, not str.
_FORMATTING_FUNCS = (
    lambda line: line.strip().upper(),   # format .fastq line
    lambda line: line.decode("utf-8").strip().upper()   # format .fastq.gz line 
)

# _blast_rep = "blastn_report.txt"
# _fasta_rep = "fasta36_report.txt"
_query = "query.fasta"
_sbjct = "subject.fasta"
_ncbi_fmt_db = "REPLACE_DB"
_constV3V4_path = "REPLACE_CONST"

class SilvaDBNotInstalledError(Exception):#{
    pass
#}

if not os.path.exists(_ncbi_fmt_db + ".nhr") or not os.path.exists(_constV3V4_path):#{
    raise SilvaDBNotInstalledError(RED + "Silva database is not installed!\n\tRun 'install_read_merging_16S.sh'" + "\033[0m")
#}


# blastn supports <= 4 threads
from multiprocessing import cpu_count
if cpu_count() >= 4:#{
    _num_cpus_blast = 4
#}
else:#{
    _num_cpus_blast = cpu_count()
#}


QSEQID, SSEQID, PIDENT, LENGTH, MISMATCH, GAPOPEN, QSTART, QEND, SSTART, SEND, EVALUE, BITSCORE, SACC, SSTRAND = range(14)
_cmd_for_blastn = """blastn -query {} -db {} -penalty -1 -reward 2 -ungapped \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sacc sstrand" \
-task megablast -max_target_seqs 1 -num_threads {}""".format(_query, _ncbi_fmt_db, _num_cpus_blast)

_draft_cmd_for_blastdbcmd = "blastdbcmd -db {} -entry REPLACE_ME -out {}".format(_ncbi_fmt_db, _sbjct)

_cmd_for_fasta = "fasta36 {} {} -n -f 20 -g 10 -m 8 -3 -r +5/-2".format(_query, _sbjct)


_RC_DICT = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'U': 'A',   # in this case it does not matter
    'N': 'N'    # just in case
}
_single_nucl_rc = lambda nucl: _RC_DICT[nucl]
_rc = lambda seq: "".join(map(_single_nucl_rc, seq[::-1]))

# We need just it's length stored in memory
_const_V3_V4_len = 0
with open(_constV3V4_path, 'r') as const_seq_file:#{
    const_seq_file.readline()
    for line in const_seq_file:#{
        _const_V3_V4_len += len(line.strip())
    #}
#}

# === These constants had been chosen empirically: ===

# This constant is used fir detecting random alignments in the center of reads:
_MAX_ALIGN_OFFSET = 40

# Aligner can miss short overlapping region, especially if there are some sequencing enrrors at the end of reads
_MIN_OVERLAP = 11


# Constant region should be in center of a merged sequence
_MAX_CONST_MID_OFFS = 70

# Constant region should align good enough, since it is constant
# (percent of identity is provided by fasta36, coverage is not and we need to calculate it manually,
# therefore coverage is if fraction form and identity -- in percent form)
_MIN_CONST_COV = 0.90
_MIN_CONST_IDENT = 80.0

# Maximum credible gap length. Very long gap isn't probable enough.
_MAX_CRED_GAP_LEN = 90


_curr_merged_seq = ""
_curr_merged_qual = ""


# The following dictionary represents some statistics of merging. 
# Keys in this dictionary are in consistency with return codes of the function "_merge_pair"
_merging_stats = None    # it in None in the beginning, because there is no statistics before merging


# Constants for naive "aligning"
_SEED_LEN = 11
_INDCS = range(_SEED_LEN)
_MAX_SHIFT = 150
_MIN_PIDENT = 0.90


# ===============================  Internal functions  ===============================


def _close_files(*files):#{
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
            for indx_or_key in obj:
                obj[indx_or_key].close()
        #}
        elif isinstance(obj, TextIOWrapper) or isinstance(obj, GzipFile):#{
            obj.close()
        #}
        else:#{
            print_yellow("""If you use function 'close_files', please, store file objects in 
    lists, tuples or dictionaries or pass file objects itself to the function.""")
            try:#{
                obj.close()
            #}
            except:#{}
                print_red("Object passed to 'close_files' function can't be closed")
                exit(1)
            #}
        #}
    #}
#}


def _write_fastq_record(outfile, fastq_record):#{
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
#}


def _read_fastq_record(read_files):#{

    if len(read_files) != 1 and len(read_files) != 2:#{
        print_red("You can pass only 1 or 2 files to the function 'read_pair_of_reads'!")
        exit(1)
    #}

    fastq_recs = dict()           # this dict should consist of two fastq-records: from R1 and from R2
    for key in read_files.keys():#{
        fastq_recs[key] = {                    #read all 4 lines of fastq-record
            "seq_id": read_files[key].readline(),
            "seq": read_files[key].readline(),
            "optional_id": read_files[key].readline(),
            "quality_str": read_files[key].readline()
        }
    #}
    return fastq_recs
#}


def _naive_align(fseq, rseq):#{
    """
    This function takes the tail of forward read and slides it through the reverse read.
    If sugnificant similarity is found -- return value of shift.
    Otherwise rerurn -1.
    :param fseq: forward read
    :type fseq: str
    :param rseq: reverse read
    :type rseq: str
    :return: int
    """    

    shift = 0

    for seed_len in range(_SEED_LEN, _SEED_LEN + _MAX_SHIFT):#{

        seed = fseq[len(fseq) - seed_len :]
        score = 0

        for pos in range(seed_len):#{
            if seed[pos] == rseq[pos]:
                score += 1
        #}

        if score / seed_len > _MIN_PIDENT:#{
            return shift
        #}
        else:#{
            shift += 1
        #}
    #}

    return -1
#}


def _al_ag_one_anoth(f_id, fseq, r_id, rseq):#{
    """
    Internal function.
    Align reads against one another with fasta36.
    :param f_id: id of the forward sequence
    :type f_id: str
    :param fseq: forward sequence itself
    :type fseq: str
    :param r_id: id of the reverse sequence
    :type r_id: str
    :param rseq: reverse sequence itself
    :type rseq: str
    :return: report generated by fasta36
    :return type: list<str> 
    """

    # form query file
    with open(_query, 'w') as query_file:#{
        query_file.write(f_id.replace('@', '>') + '\n')
        query_file.write(fseq)
    #}

    # form subject file
    with open(_sbjct, 'w') as sbjct_file:#{
        sbjct_file.write(r_id.replace('@', '>') + '\n')
        sbjct_file.write(rseq)
    #}

    # align forward read against reverse read
    pipe = os.popen(_cmd_for_fasta)
    far_report = pipe.readline().strip().split('\t')  # "far" means Forward [read] Against Reverse [read]
    pipe.close()

    return far_report
#}


def _merge_by_overlap(loffset, overl, fseq, fqual, rseq, rqual):#{
    """
    Internal function.
    Merge two reads according to their overlapping redion.
    Function leaves nucleotide with higher quality.
    :param loffset: number of nucleotides in forward read from it's beginning to start of the overlapping region
    :type loffset: int
    :param overl: number of nucleotides in the overlapping region
    :type overl: int
    :param fseq: forward read itself
    :type fseq: str
    :param fqual: quality string of the forward read
    :type fqual: str
    :param rseq: reverse read itself
    :type rseq: str
    :param rqual: quality string of the reverse read
    :type rqual: str
    """

    # beginning comes from forward read
    merged_seq = fseq[: loffset]
    merged_qual = fqual[: loffset]

    # "recover" overlapping region depending on quality 
    for i, j in zip(range(loffset, loffset + overl), range(overl)):#{
        if ord(fqual[i]) >= ord(rqual[j]):#{
            qual, nucl = fqual[i], fseq[i]
        #}
        else:#{
            qual, nucl = rqual[j], rseq[j]
        #}
        merged_seq += nucl
        merged_qual += qual
    #}

    # end comes from reverse read
    merged_seq += rseq[overl :]
    merged_qual += rqual[overl :]

    return (merged_seq, merged_qual)
#}


class NoRefAlignError(Exception):#{
    pass
#}


def _blast_and_align(rseq, r_id):#{
    """
    Internal function.
    Sequence of actions performed by this function:
    1. Blast forward read against the database.
    2. Find and get reference sequence, against which forward read has aligned with the better score (e. i. the first in list).
    3. Align reverse read agaist this reference sequence.
    4. Return results of these alignments.
    :param rseq: reverse read
    :type rseq: str
    """

    # Assume that query is still the same.
    # Blast forward read.
    pipe = os.popen(_cmd_for_blastn)   # query is still the same
    faref_report = pipe.readline().strip().split('\t')
    pipe.close()

    if faref_report[0] == '':#{
        raise NoRefAlignError()
    #}

    # get ref sequence
    cmd_for_blastbdcmd = _draft_cmd_for_blastdbcmd.replace("REPLACE_ME", faref_report[SACC])
    os.system(cmd_for_blastbdcmd)
    if faref_report[SSTRAND] == "minus":#{
        with open(_sbjct, 'r') as sbjct_file:#{
            sbjct_seq_id = sbjct_file.readline().strip()
            sbjct_seq = ""
            for line in sbjct_file:#{
                sbjct_seq += line.strip()
            #}
        #}
        sbjct_seq = _rc(sbjct_seq)
        with open(_sbjct, 'w') as sbjct_file:#{
            sbjct_file.write(sbjct_seq_id + '\n')
            sbjct_file.write(sbjct_seq)
        #}
    #}

    # align reverse read against the reference
    with open(_query, 'w') as query_file:#{
        query_file.write(r_id.replace('@', '>') + '\n')
        query_file.write(rseq)
    #}

    pipe = os.popen(_cmd_for_fasta)
    raref_report = pipe.readline().strip().split('\t')   # raref means Reverse [read] Against REFerence [sequence]
    pipe.close()

    return faref_report, raref_report
#}


def _search_for_constant_region(loffset, overl, fseq, fqual, rseq, rqual):#{
    """
    Internal function.
    This function tryes to "merge-by-overlap" forward and reverse reads and searches for
        the consensus sequence of constant region of 16S rRNA gene located between
        V3 and V4 variable regions.
    If there is a constant region there -- consider reads as properly merged and return them.
    Else -- do not merge them and return 1.
    :param loffset: number of nucleotides in forward read from it's beginning to start of the overlapping region
    :type loffset: int
    :param overl: number of nucleotides in the overlapping region
    :type overl: int
    :param fseq: forward read itself
    :type fseq: str
    :param fqual: quality string of the forward read
    :type fqual: str
    :param rseq: reverse read itself
    :type rseq: str
    :param rqual: quality string of the reverse read
    :type rqual: str
    :return: True if there is a constant region, otherwise False
    """

    # try to merge
    merged_seq, merged_qual = _merge_by_overlap(loffset, overl, fseq, fqual, rseq, rqual)

    # form query file
    with open(_sbjct, 'w') as sbjct_file:#{
        sbjct_file.write(">presumably_merged_read\n")
        sbjct_file.write(merged_seq)
    #}

    # form subject file
    os.system("cat {} > {}".format(_constV3V4_path, _query))

    # align
    pipe = os.popen(_cmd_for_fasta.replace(" -3", ""))  # in this case we need both strands
    pmer_report = pipe.readline().strip().split('\t')  # "pmer" means Presumably MERged [read]
    pipe.close()

    # check if there are a constant region
    enough_cov = float(pmer_report[LENGTH]) / _const_V3_V4_len > _MIN_CONST_COV
    enough_indent = float(pmer_report[PIDENT]) > _MIN_CONST_IDENT
    # constant region should be in the center of merged sequence:
    in_the_midl = int(pmer_report[SSTART]) > _MAX_CONST_MID_OFFS
    in_the_midr = int(pmer_report[SEND]) > len(merged_seq) - _MAX_CONST_MID_OFFS


    if enough_cov and enough_indent and in_the_midl and in_the_midr:#{
        _curr_merged_seq = merged_seq
        _curr_merged_qual = merged_qual
        return True
    #}
    else:#{
        # can't merge
        return False
    #}
#}


def _handle_unforseen_case(f_id, fseq, r_id, rseq):#{
    """
    Internal function.
    Handle unforseen case and Provide man who uses the program with information 
        about how erroneous reads align one against another.
    :param f_id: id of the forward sequence
    :type f_id: str
    :param fseq: forward sequence itself
    :type fseq: str
    :param r_id: id of the reverse sequence
    :type r_id: str
    :param rseq: reverse sequence itself
    :type rseq: str
    :return: void
    """
    error_report = "error_report.txt"
    print_red("Read merging crashed. It is my fault. Report to me -- I will fix it.")
    print_red("You can get some info about reads caused this crash in the following file:\n\t'{}'"
        .format(os.path.abspath(error_report)))

    with open(error_report, 'w') as errfile:#{
        errfile.write("Forward read:\n\n" + f_id + '\n')
        errfile.write(fseq + "\n\n" + '~' * 40 + '\n')
        errfile.write("Reverse read:\n\n" + r_id + '\n')
        errfile.write(rseq + "\n\n" + ('~' * 40 + '\n') * 2 + '\n')
        errfile.write("ALIGNMENT (FORWARD READ AGAINST ERVERSE ONE):\n\n")
    #}

    # ===  Provide man who uses the program with information about how erroneous reads align one against another  ===
    # form query file
    with open(_query, 'w') as query_file:#{
        query_file.write(f_id.replace('@', '>') + '\n')
        query_file.write(fseq)
    #}

    # form subject file
    with open(_sbjct, 'w') as sbjct_file:#{
        sbjct_file.write(r_id.replace('@', '>') + '\n')
        sbjct_file.write(rseq)
    #}

    # align forward read against reverse read and append result to error report
    os.system("fasta36 {} {} -n -f 20 -g 10 -3 >> {}".format(_query, _sbjct, error_report))
#}


def _del_temp_files(put_artif_dir):#{
    """
    Delete temporary files used by functions of these module.
    """
    for file in _query, _sbjct:#{
        if os.path.exists(file):#{
            os.remove(file)
        #}
    #}
    if put_artif_dir is not None:#{
        for file in os.listdir(put_artif_dir):#{
            if os.path.exists(os.sep.join([put_artif_dir, file])) and "roughly" in file:#{
                os.remove(os.sep.join([put_artif_dir, file]))
            #}
        #}
    #}
#}


def _naive_merging(fastq_recs):#{
    """
    The first of the "kernel" functions in this module. Performs rough process of read merging.
    :param fastq_reqs: a dictionary of two fastq-records stored as dictionary of it's fields
    :type fastq_reads: dict<str: dict<str, str>>

    The description of this weird parameter:
    The outer dictionary contains two dictionaries accessable by the following keys: "R1" (forward read), "R2" (reverse read).
    Fields of each record should be accessable by the follosing keys:
    1) "seq_id" (ID of the read)
    2) "seq" (sequence itsef)
    3) "optional_id" (the third line, where '+' is usually written)
    4) "quality_str" (quality string in Phred33)
    # Return values:
    # 0 -- if reads can be merged;
    # 1 -- if reads can't be merged;
    # 2 -- too short
    """

    f_id = fastq_recs["R1"]["seq_id"]
    fseq = fastq_recs["R1"]["seq"]
    fqual = fastq_recs["R1"]["quality_str"]
    r_id = fastq_recs["R2"]["seq_id"]
    rseq = _rc(fastq_recs["R2"]["seq"])         # reverse-complement
    rqual = fastq_recs["R2"]["quality_str"][::-1]      # reverse


    # |==== Firstly try to find overlapping region with naive-but-relatively-fast algorithm ====|

    # If reads are too short -- IndexError will be raised and these reads will be considered too short.
    # We do not need these reads anyway.
    try:#{
        shift = _naive_align(fseq, rseq)
    #}
    except IndexError:#{
        return 2
    #}

    if shift != -1:#{
        overl = shift + _SEED_LEN
        loffset = len(fseq) - overl

        merged_seq, merged_qual = _merge_by_overlap(loffset, overl, fseq, fqual, rseq, rqual)
        
        globals()["_curr_merged_seq"] = merged_seq
        globals()["_curr_merged_qual"] = merged_qual
        return 0
    #}
    else:#{
        return 1
    #}
#}


def _accurate_merging(fastq_recs, V3V4=False):#{}
    """
    The second of the "kernel" functions in this module. Performs accurate process of read merging.
    :param fastq_reqs: a dictionary of two fastq-records stored as dictionary of it's fields
    :type fastq_reads: dict<str: dict<str, str>>
    Description of this weird parameter:
    The outer dictionary contains two dictionaries accessable by the following keys: "R1" (forward read), "R2" (reverse read).
    Fields of each record should be accessable by the follosing keys:
    1) "seq_id" (ID of the read)
    2) "seq" (sequence itsef)
    3) "optional_id" (the third line, where '+' is usually written)
    4) "quality_str" (quality string in Phred33)
    # Return values:
    # 0 -- if reads can be merged;
    # 1 -- if reads can't be merged;
    # 2 -- too short
    # 3 -- fatal error, unforseen case
    """

    f_id = fastq_recs["R1"]["seq_id"]
    fseq = fastq_recs["R1"]["seq"]
    fqual = fastq_recs["R1"]["quality_str"]
    r_id = fastq_recs["R2"]["seq_id"]
    rseq = _rc(fastq_recs["R2"]["seq"])         # reverse-complement
    rqual = fastq_recs["R2"]["quality_str"][::-1]      # reverse


    far_report = _al_ag_one_anoth(f_id, fseq, r_id, rseq)


    # |==== Check how they have aligned ====|
    
    # Consider normal overlap as:
    # FFFFFFFF----
    # ----RRRRRRRR
    reads_normally_overlap = True   # naive assumption

    # Discard following situation:
    # --FFFFFF
    # RRRRRR--
    too_short = int(far_report[QSTART]) < int(far_report[SSTART])
    too_short = too_short or int(far_report[QEND]) / len(fseq) < int(far_report[SEND]) / len(rseq)

    reads_normally_overlap = reads_normally_overlap and not too_short

    # Catch randomly occured alignment in center of sequences, 
    # e. i. reads don't align against one another in a proper way
    rand_align = len(fseq) - int(far_report[QEND]) > _MAX_ALIGN_OFFSET or int(far_report[SSTART]) > _MAX_ALIGN_OFFSET
    # resulting conslusion
    reads_normally_overlap = reads_normally_overlap and not rand_align

    if reads_normally_overlap and float(far_report[PIDENT]) < 75.0:#{
        return 1
    #}


    # |==== Decide what to do according to "analysis" above. ====|

    # === Ok, reads overlap in the way they should do it. ===
    # FFFFFFFF----
    # ----RRRRRRRR
    if reads_normally_overlap:#{

        loffset = int(far_report[QSTART]) - int(far_report[SSTART])  # offset from the left
        # in the next line we have 1-based terms, hense I need to substract 1:
        overl = int(far_report[LENGTH]) + int(far_report[SSTART]) + (len(fseq) - int(far_report[QEND])) - 1

        merged_seq, merged_qual = _merge_by_overlap(loffset, overl, fseq, fqual, rseq, rqual)
        globals()["_curr_merged_seq"] = merged_seq
        globals()["_curr_merged_qual"] = merged_qual

        return 0
    #}

    # === Randomly occured alignment in center of sequences, need to blast ===
    elif rand_align:#{

        # === Blast forward read, align reverse read against the reference ===
        # "faref" means Forward [read] Against REFerence [sequence]
        # "raref" means Reverse [read] Against REFerence [sequence]
        try:#{
            faref_report, raref_report = _blast_and_align(rseq, r_id)
        #}
        except NoRefAlignError:#{
            return 1
        #}

        if float(faref_report[EVALUE]) > 0.05 or float(raref_report[EVALUE]) > 0.05:#{
            return 1
        #}


        # Calculate some "features".
        # All these "features" are in coordinates of the reference sequence
        forw_start = int(faref_report[SSTART]) - int(faref_report[QSTART])
        forw_end = forw_start + len(fseq)
        rev_start = int(raref_report[SSTART]) - int(raref_report[QSTART])
        rev_end = rev_start + len(rseq)

        gap = forw_end < rev_start

        if forw_end > rev_end:#{
            return 2
        #}

        # ===  Handle alignment results ====

        # Overlap region is long enough
        if not gap and forw_end - rev_start > _MIN_OVERLAP:#{
            if V3V4 and forw_start < rev_start:#{
                # Try to merge them and search for constant region in merged sequence
                loffset = rev_start - forw_start
                overl = forw_end - rev_start
                reply = _search_for_constant_region(loffset, overl, fseq, fqual, rseq, rqual)
                return 0 if reply else 1
            #}
            else:#{
                # can't merge reads
                return 1
            #}
        #}


        # Length of the overlapping region is short,
        # but reverse read alignes as they should do it, so let it be.
        elif not gap and forw_end - rev_start <= _MIN_OVERLAP:#{
            
            loffset = rev_start - forw_start
            overl = forw_end - rev_start 

            merged_seq, merged_qual = _merge_by_overlap(loffset, overl, fseq, fqual, rseq, rqual)
            globals()["_curr_merged_seq"] = merged_seq
            globals()["_curr_merged_qual"] = merged_qual
            return 0
        #}
            
        # Here we have a gap. 
        elif gap:#{
            
            gap_len = rev_start - forw_end
            # Very long gap isn't probable enough.
            if gap_len > _MAX_CRED_GAP_LEN:#{
                return 1
            #}
            else:#{
                # Fill in the gap with 'N'.
                merged_seq = fseq + 'N' * (gap_len - 1) + rseq
                merged_qual = fqual + chr(33) * (gap_len - 1) + rqual   # Illumina uses Phred33

                globals()["_curr_merged_seq"] = merged_seq
                globals()["_curr_merged_qual"] = merged_qual
                return 0
            #}
        #}

        # === Unforseen case. This code should not be ran. But if it does-- report about it. ===
        _handle_unforseen_case(f_id, fseq, r_id, rseq)
        return 3
    #}

    # Too short sequence, even if is merged. We do not need it.
    elif too_short:#{
        # --FFFFFFF
        # RRRRRRR--
        return 2
    #}

    # === Unforseen case. This code should not be ran. But if it'll do -- report about it. ===
    _handle_unforseen_case(f_id, fseq, r_id, rseq)
    return 3
#}


def _handle_merge_pair_result(merging_result, fastq_recs, result_files, accurate=False):#{
    """
    This function handles the result of read merging and writes sequences in corresponding files.

    :param merging_result: a tuple of strings (if reads are merged), or an integer, if they cannot be merged
    :param fastq_reqs: a dictionary of two fastq-records stored as dictionary of it's fields
    :type fastq_reads: dict<str: dict<str, str>>
    :type result_files: dict<str: _io.TextIOWrapper>
    :type read_files: dict<str: _io.TextIOWrapper>
    """

    # if reads are merged
    if merging_result == 0:#{

        merged_rec = {
            "seq_id": fastq_recs["R1"]["seq_id"],
            "seq": _curr_merged_seq,
            "optional_id": '+',
            "quality_str": _curr_merged_qual
        }
        _write_fastq_record(result_files["merg"], merged_rec)
        _merging_stats[merging_result] += 1
        if accurate:#{
            _merging_stats[1] -= 1
        #}
        return 0
    #}

    # can't merge reads
    elif merging_result == 1:#{
        _write_fastq_record(result_files["umR1"], fastq_recs["R1"])
        _write_fastq_record(result_files["umR2"], fastq_recs["R2"])
        if not accurate:#{
            _merging_stats[merging_result] += 1
        #}
        return 1
    #}

    # if resulting sequence is too short to distinguish taxa
    elif merging_result == 2:#{
        _write_fastq_record(result_files["shrtR1"], fastq_recs["R1"])
        _write_fastq_record(result_files["shrtR2"], fastq_recs["R2"])
        _merging_stats[merging_result] += 1
        if accurate:#{
            _merging_stats[1] -= 1
        #}
        return 2
    #}

    # if unforseen situation occured in 'read_merging_16S'
    elif merging_result == 3:#{
        _close_files(result_files)
        input("Press ENTER to exit:")
        exit(1)
    #}

    # if 'read_merging_16S' returnes something unexpected and undesigned
    else:#{
        print_red("ERROR!!!\n\tModule that was merging reads returned an unexpected and undesigned value.")
        try:#{
            print("\tHere it is: {}".format(str(merging_result)))
        #}
        except:#{
            print("\tUnfortunately, it cannot be printed. It is of {} type".format(type(merging_result)))
        #}
        print("It is my fault. Report to me -- I will fix it.")
        _close_files(read_files, result_files)
        input("Press ENTER to exit:")
        exit(1)
    #}
#}


def _check_file_existance(path):#{

    if os.path.exists(path):#{
        return
    #}
    else:#{
        print_red("\nFile '{}' does not exist!".format(path))
        exit(1)
    #}
#}


def _open_files(read_paths, result_paths, hto=open, wmode='w'):#{

    file_type = int(match(r".*\.gz$", read_paths["R1"]) is not None)
    how_to_open = _OPEN_FUNCS[file_type]

    read_files = dict()
    result_files = dict()
    try:#{
        for key in read_paths.keys():#{
            read_files[key] = how_to_open(read_paths[key])
        #}
        for key in result_paths.keys():#{
            result_files[key] = open(result_paths[key], wmode)
        #}
    #}
    except OSError as oserror:#{
        print_red("Error while opening file", str(oserror))
        _close_files(read_files, result_files)
        input("Press enter to exit:")
        exit(1)
    #}
    return read_files, result_files
#}


# ===============================  "Public" stuff  ===============================


class AccessStatsBeforeMergingError(Exception):#{
    "Exception raising on the attempt of accessing merging statisttics before merging itself"
    pass
#}


def get_merging_stats():#{
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
    if _merging_stats is not None:#{
        return _merging_stats
    #}
    else:#{
        raise AccessStatsBeforeMergingError(RED + "You cannot access merging statistics before merging is completed!\a" + ENDCLR)
    #}
#}


def merge_reads(R1_path, R2_path,
    outdir_path="read_merging_result_{}".format(_now).replace(' ', '_'), 
    V3V4=False):#{
    """
    This is the function that you should actually call from the outer scope in order to merge reads
    (and 'get_merging_stats' after it, if you want).

    :param R1_path: path to a file with forward reads
    :type R1_path: str
    :param R2_path: path to a file with reverse reads
    :type R2_path: str
    :param outdir_path: path to a directory, in which result files will be stored
    :type outdir_path: str
    :param V3V4: if True, more accurate merging can be performed,
        but only if target sequences contain V3 and V4 regions and a constant region between them.
        This parameter is False by default.
    :type V3V4: bool

    :return: function returns a dict<str: str> of the following format:
    {   
        "merg": path to a file with successfully merged reads,
        "umR1": path to a file with forward unmerged reads,
        "umR2": path to a file with reverse unmerged reads,
        "shrtR1": path to a file with forward reads considered as too short,
        "shrtR2": path to a file with reverse reads considered as too short
    }
    """

    # Create output directory
    if not os.path.exists(outdir_path):#{
        try:#{
            os.mkdir(outdir_path)
        #}
        except OSError as oserror:#{
            print_red("Error while creating result directory\n", str(oserror))
            input("Press ENTER to exit:")
            exit(1)
        #}
    #}

    # Create a directory for putative artifacts
    artif_dir = "{}{}putative_artifacts".format(outdir_path, os.sep)
    if not os.path.exists(artif_dir):#{
        try:#{
            os.mkdir(artif_dir)
        #}
        except OSError as oserror:#{
            print_red("Error while creating result directory\n", str(oserror))
            input("Press ENTER to exit:")
            exit(1)
        #}
    #}

    # The followig dictionary represents some statistics of merging. 
    # Keys in this dictionary are in consistency with return codes of the function "_merge_pair"
    globals()["_merging_stats"] = {
        0: 0,           # number of merged reads
        1: 0,           # number of unmerged reads
        2: 0            # number of too short sequences
    }

    read_paths = {
        "R1": R1_path,
        "R2": R2_path
    }

    names = {
        "R1": os.path.basename(read_paths["R1"])[: search(r"f(ast)?q", os.path.basename(read_paths["R1"])).span()[0]-1],
        "R2": os.path.basename(read_paths["R2"])[: search(r"f(ast)?q", os.path.basename(read_paths["R2"])).span()[0]-1]
    }
    # Name without "__R1__" and "__R2__":
    more_common_name = os.path.basename(read_paths["R1"])[: os.path.basename(read_paths["R1"]).find("R1")]

    merg_path = "{}{}{}.16S.merged.fastq".format(outdir_path, os.sep, more_common_name)
    shrtR1_path = "{}{}{}.too_short.fastq".format(artif_dir, os.sep, names["R1"])
    shrtR2_path = "{}{}{}.too_short.fastq".format(artif_dir, os.sep, names["R2"])

    result_paths = {
        # File for merged sequences:
        "merg": merg_path,
        # Files for unmerged sequences:
        "umR1": "{}{}{}.roughly_unmerged.fastq".format(artif_dir, os.sep, names["R1"]),
        "umR2": "{}{}{}.roughly_unmerged.fastq".format(artif_dir, os.sep, names["R2"]),
        # Files for those reads, that form too short sequence after merging:
        "shrtR1": shrtR1_path,
        "shrtR2": shrtR2_path
    }

    # Collect some info
    file_type = int(match(r".*\.gz$", read_paths["R1"]) is not None)
    how_to_open = _OPEN_FUNCS[file_type]
    actual_format_func = _FORMATTING_FUNCS[file_type]
    readfile_length = sum(1 for line in how_to_open(read_paths["R1"], 'r'))  # lengths of read files are equal
    read_pairs_num = int(readfile_length / 4)            # divizion by 4, because there are 4 lines per one fastq-record

    # Open files
    read_files, result_files = _open_files(read_paths, result_paths, wmode='w')

    print_yellow("\nRead merging started")

    # Proceed
    inc_percentage = 0.01
    reads_processed = 0
    spaces = 50
    next_done_percentage = inc_percentage
    print("\nProceeding...\n\n")
    print("[" + " "*50 + "]" + "  0%", end="")
    while reads_processed < read_pairs_num:#{

        fastq_recs = _read_fastq_record(read_files)

        for file_key in fastq_recs.keys():#{
            for field_key in fastq_recs[file_key].keys():
                fastq_recs[file_key][field_key] = actual_format_func(fastq_recs[file_key][field_key])
        #}

        merging_result = _naive_merging(fastq_recs)
        _handle_merge_pair_result(merging_result, fastq_recs, result_files, accurate=False)

        reads_processed += 1
        if reads_processed / read_pairs_num >= next_done_percentage:#{
            count = round(next_done_percentage * 100)
            spaces = 50 - int(count/2)
            next_done_percentage += inc_percentage
            print("\r[" + "="*int(count/2) + ">" + " "*spaces + "]" + "  {}% ({}/{} read pairs are processed)"
                .format(count, reads_processed, read_pairs_num), end="")
        #}
    #}

    print("\r[" + "="*50 + "]" + "  100% ({}/{} read pairs are processed)\n\n"
            .format(reads_processed, read_pairs_num))
    print_green("\nRead merging is completed")
    print("""\n{} read pairs have been merged together
{} read pairs haven't been merged together
{} read pairs have been considered as too short"""
.format(_merging_stats[0], _merging_stats[1], _merging_stats[2]))
    print('\n' + '~' * 50 + '\n')
    _close_files(read_files, result_files)
    _del_temp_files(None)


    # |==== Accurate merging ===|

    read_paths = {
        "R1": result_paths["umR1"],
        "R2": result_paths["umR2"]
    }

    result_paths["umR1"] = result_paths["umR1"].replace("roughly_", "")
    result_paths["umR2"] = result_paths["umR2"].replace("roughly_", "")

    # Collect some info
    file_type = int(match(r".*\.gz$", read_paths["R1"]) is not None)
    how_to_open = _OPEN_FUNCS[file_type]
    actual_format_func = _FORMATTING_FUNCS[file_type]
    readfile_length = sum(1 for line in how_to_open(read_paths["R1"], 'r'))  # lengths of read files are equal
    read_pairs_num = int(readfile_length / 4)            # divizion by 4, because there are 4 lines per one fastq-record

    # Open files
    read_files, result_files = _open_files(read_paths, result_paths, wmode='a')

    print_yellow("""\nNow the program will merge the rest of reads 
    with accurate-but-slow strategy""")
    print("\tIt will take a while")

    # Proceed
    reads_processed = 0
    spaces = 50
    next_done_percentage = inc_percentage
    count = 0
    print("\nProceeding...\n\n")
    print("[" + " "*50 + "]" + "  0%", end="")
    while reads_processed < read_pairs_num:#{

        fastq_recs = _read_fastq_record(read_files)

        for file_key in fastq_recs.keys():#{
            for field_key in fastq_recs[file_key].keys():#{
                fastq_recs[file_key][field_key] = actual_format_func(fastq_recs[file_key][field_key])
            #}
        #}

        merging_result = _accurate_merging(fastq_recs, V3V4=V3V4)
        _handle_merge_pair_result(merging_result, fastq_recs, result_files, accurate=True)

        reads_processed += 1
        if reads_processed / read_pairs_num >= next_done_percentage:#{
            count = round(next_done_percentage * 100)
            spaces = 50 - int(count/2)
            next_done_percentage += inc_percentage
        #}
        # always print read pairs number
        print("\r[" + "="*int(count/2) + ">" + " "*spaces + "]" + "  {}% ({}/{} read pairs are processed)"
            .format(count, reads_processed, read_pairs_num), end="")
    #}


    print("\r[" + "="*50 + "]" + "  100% ({}/{} read pairs are processed)\n\n"
            .format(reads_processed, read_pairs_num))

    print_green("\nRead merging is completed")
    print("\nFinally,")
    print("""\t{} read pairs have been merged together
\t{} read pairs haven't been merged together
\t{} read pairs have been considered as too short"""
.format(_merging_stats[0], _merging_stats[1], _merging_stats[2]))
    print('\n' + '~' * 50 + '\n')
    _close_files(read_files, result_files)
    _del_temp_files(artif_dir)


    return result_paths
#}



# |==== 'read_merging_16S.py' can be used as script ====|


if __name__ == "__main__":#{

    from sys import version_info
    if (version_info.major + 0.1 * version_info.minor) < (3.0 - 1e-6):#{  # just-in-case 1e-6 substraction
        print_red("\t\nATTENTION!\nThis script cannot be executed by python interpreter version < 3.0!\a")
        print_red("\tYour python version: {}.{}".format(version_info.major, version_info.minor))
        print_red("Try '$ python3 preprocess16S.py' or update your interpreter.")
        exit(0)
    #}

    import os
    import getopt
    from sys import argv

    usage_msg = "Usage:\n\tpython read_merging_16S.py -1 forward_R1_reads.fastq.gz -2 reverse_R2_reads.fastq.gz [-o ourdir]"

    try:#{
        opts, args = getopt.getopt(argv[1:], "h1:2:o:", ["V3-V4", "R1=", "R2=", "outdir="])
    #}
    except getopt.GetoptError as opt_err:#{
        print_red(opt_err + '\a')
        print(usage_msg)
        exit(2)
    #}

    outdir_path = "{}{}read_merging_16S_result_{}".format(os.getcwd(), os.sep, _now).replace(" ", "_") # default path
    read_paths = dict()
    V3V4 = False

    for opt, arg in opts:#{
        if opt in ("-h", "--help"):#{
            pritn(usage_msg)
            exit(0)
        #}

        elif opt in ("-o", "--outdir"):#{
            outdir_path = os.path.abspath(arg)
        #}

        elif opt == "--V3-V4":#{
            V3V4 = True
        #}

        elif opt in ("-1", "--R1"):#{
            if not "-2" in argv and not "--R2" in argv:#{
                print_red("\nATTENTION!\n\tYou should specify both forward and reverse reads!\a")
                exit(1)
            #}
            _check_file_existance(arg)
            read_paths["R1"] = arg
        #}

        elif opt in ("-2", "--R2"):#{
            if not "-1" in argv and not "--R1" in argv:#{
                print_red("\nATTENTION!\n\tYou should specify both forward and reverse reads!\a")
                exit(1)
            #}
            _check_file_existance(arg)
            read_paths["R2"] = arg
        #}
    #}

    # Check utilities for read merging
    pathdirs = os.environ["PATH"].split(os.pathsep)
    for utility in ("fasta36", "blastn", "blastdbcmd"):#{
        utility_found = False
        for directory in pathdirs:#{
            if utility in os.listdir(directory):#{
                utility_found = True
                break
            #}
        #}
        if not utility_found:#{
            print_red("\tAttention!\n{} is not found in your system.\a".format(utility))
            print("If you want to use this feature, please install {}".format(utility))
            print("""If this error still occure although you have installed everything 
    -- make sure that this program is added to PATH)""")
        #}
    #}

    print_yellow("\nResult files will be placed in the following directory:\n\t'{}'".format(outdir_path))

    # Proceed
    result_files = merge_reads(read_paths["R1"], read_paths["R2"], outdir_path=outdir_path, V3V4=V3V4)

    # Remove empty files
    for file in result_files.values():#{
        if os.stat(file).st_size == 0:#{
            os.remove(file)
            print("'{}' is removed since it is empty".format(file))
        #}
    #}

    # Gzip result files
    print_yellow("\nGzipping result files...")
    for file in result_files.values():#{
        if os.path.exists(file):#{
            os.system("gzip -f {}".format(file))
            print("\'{}\' is gzipped".format(file))
        #}
    #}
    print_green("Gzipping is completed\n")
    print("Result files are placed in the following directory:\n\t{}".format(outdir_path))

    # Write log file
    with open("{}{}read_merging_16S_{}.log".format(outdir_path, os.sep, _now).replace(" ", "_"), 'w') as logfile:#{
        logfile.write("Script 'read_merging_16S.py' was ran on {}\n".format(_now.replace('.', ':')))
        logfile.write("{} read pairs have been merged.\n".format(_merging_stats[0]))
        logfile.write("{} read pairs haven't been merged together.\n".format(_merging_stats[1]))
        logfile.write("{} read pairs have been considered as too short for merging.\n".format(_merging_stats[2]))
    #}
#}