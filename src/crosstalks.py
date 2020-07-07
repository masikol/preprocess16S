# -*- coding: utf-8 -*-
# Module for detecting cross-talks depending on presence of proper primer sequences in reads.

from src.fastq import *
import re


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
    "T": "TUNYWKBDH",
    "U": "TUNYWKBDH",
    "N": ""
}
# This dictionary is necessary beacause primers used in 16S amplicon researches are often redundant.
# Keys in this dictionary are nucleotides that read sequence may contain.
# Values are lists (I do not mean build-in python type, obviously) of symbols a redundant primer
#   sequence may contain and match corresponding key-nucleotide.
# E.g., if we have got A from read sequence and one of ANRWMDHV symbols is located at
#    corresponding position in primer sequence, this A-nucleotide will be considered as matching.
# 'N'-nucleotide in read sequence can not match anything.


def get_primers(argv, primer_path):
    """
    Function returns primers sequeces and IDs. It returns default Illumina 16S amplicon primers
      if no primers were specified with '-p' option. Otherwise it parses primers from '-p' file.

    :param argv: command-line arguments [1:];
    :type argv: same type as sys.argv;
    :param primer_path: path to multi-fasta file with primers;
    :type primer_path: str;

    Returns tuple of two lists:
    ([<primer_1_seq>, <primer_2_seq>], [<primer_1_ID>, <primer_2_ID>])
    """

    # Use Illumina V3-V4 primers if they are not specified
    if not "-p" in argv and not "--primers" in argv:
        primers = ["CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"]     # list of required primer sequences
        primer_ids = [">16S Amplicon PCR Forward Primer", ">16S Amplicon PCR Reverse Primer"]  # list of their names
    else:
        # Variable named 'file_type' below will be 0(zero) if primers are in plain FASTA file
        #   and 1(one) if they are in fasta.gz file
        # If primers are in .gz file, it will be opened as .gz file and as standard text file otherwise.
        file_type = get_archv_fmt_indx(primer_path)
        how_to_open = OPEN_FUNCS[file_type]
        actual_format_func = FORMATTING_FUNCS[file_type]

        primers = list()
        primer_ids = list()
        # Assuming that each primer sequence is written in one line
        primer_file = None
        try:
            primer_file = how_to_open(primer_path, 'r')
            for i, line in enumerate(primer_file):
                if i % 2 == 0:                                            # if line is sequense id
                    primer_ids.append(actual_format_func(line))
                
                else:                                                     # if line is a sequence
                    line = actual_format_func(line).upper()
                    err_set = set(fre.indall(r"[^ATGCRYSWKMBDHVN]", line))     # primer validation
                    if len(err_set) != 0:
                        print_error("There are some inappropriate symbols in your primers!")
                        print("Here they are:")
                        for err_symb in err_set:
                            print(" " + err_symb)
                        # end for
                        
                        print("See you later with correct data")
                        input("Press enter to exit:")
                        exit(1)
                    # end if
                    
                    primers.append(line)
                # end if
            # end for
        except OSError as oserror:
            print("Error while reading primer file.\n", str(oserror))
            input("Press enter to exit:")
            exit(1)
        finally:
            if primer_file is not None:
                primer_file.close()
            # end if
        # end try
    # end if


    # The last line in fasta file can be a blank line.
    # We do not need it.
    try:
        primer_ids.remove('')
    except ValueError:
        pass
    # end try

    return ( primers, primer_ids )

# end def get_primers


def find_primer(primer, fastq_rec, keep_primers):
    """
    This function figures out, whether a primer sequence is at 5'-end of a read passed to it.
    Moreover it trims primer sequences if there are any.
    Function modifies fastq_rec dictionary if only keep_primers is not True.

    :param primer: required primer sequence
    :type primer: str
    :param fastq_rec: fastq-record containing read, which this function searches for required primer in
    :type fastq_rec: dict<str: str>
    :param keep_primers: logical value indicating whether we need to keep primer sequences in reads. False -- trim them;
    :type keep_primers: bool;


    Returns logical value:
    True if primer sequence is found in read, otherwise returns False.
    """

    read = fastq_rec["seq"]
    primer_len = len(primer)

    for shift in range(0, MAX_SHIFT + 1):
        score_1, score_2 = 0, 0

        for pos in range(0, primer_len - shift):
            # if match, 1(one) will be added to score, 0(zero) otherwise
            score_1 += primer[pos+shift] in MATCH_DICT[read[pos]]
        # end for
        
        if score_1 / (primer_len-shift) >= RECOGN_PERCENTAGE:
            if not keep_primers:
                cutlen = len(primer) - shift
                fastq_rec["seq"] = read[cutlen : ]
                fastq_rec["qual_str"] = fastq_rec["qual_str"][cutlen : ]
                return True
            
            else:
                return True
            # end if
        # end if

        # If there is no primer in 0-position, loop below will do unnecessary work.
        # Let it go a bit further.
        shift += 1

        for pos in range(0, primer_len - shift):
            score_2 += primer[pos] in MATCH_DICT[read[pos + shift]]
        # end for
        
        if score_2 / (primer_len) >= RECOGN_PERCENTAGE:
            if not keep_primers:
                cutlen = len(primer) - shift
                fastq_rec["seq"] = read[cutlen : ]
                fastq_rec["qual_str"] = fastq_rec["qual_str"][cutlen : ]
                return True
            
            else:
                return True
            # end if
        # end if
    # end for
    return False
# end def find_primer


def find_primer_organizer(fastq_recs, result_files, **kwargs):
    """
    Function to be passed to progress_counter for "cross-talk" removing.
    """

    # "Parse" kwargs
    stats = kwargs["stats"]
    primers = kwargs["primers"]
    keep_primers = kwargs["keep_primers"]

    primer_in_R1 = find_primer(primers[0], fastq_recs["R1"], keep_primers)
    if primer_in_R1:

        primer_in_R2 = find_primer(primers[1], fastq_recs["R2"], keep_primers)

        if primer_in_R2:
            write_fastq_record(result_files["mR1"], fastq_recs["R1"])
            write_fastq_record(result_files["mR2"], fastq_recs["R2"])
            stats["match"] += 1
            return
        # end if
    # end if
    write_fastq_record(result_files["trR1"], fastq_recs["R1"])
    write_fastq_record(result_files["trR2"], fastq_recs["R2"])
    stats["trash"] += 1
# end def find_primer_organizer