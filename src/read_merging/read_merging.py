# -*- coding: utf-8 -*-

from src.read_merging.ngmerge_quality_profile import match_prof, mismatch_prof


def _merge_by_overlap(loffset, fseq, fqual, rseq, rqual, phred_offset):
    # Function merges two reads according to their overlapping redion.
    # Function leaves nucleotide with higher quality.
    # :param loffset: number of nucleotides in forward read from it's beginning to start of the overlapping region;
    # :type loffset: int;
    # :param fseq: forward read itself;
    # :type fseq: str;
    # :param fqual: quality string of the forward read;
    # :type fqual: str;
    # :param rseq: reverse read itself;
    # :type rseq: str;
    # :param rqual: quality string of the reverse read;
    # :type rqual: str;
    # Returns tuple of merged sequence itself and it's quality line;
    # Return type: tuple<str>;

    # 5'-end comes from forward read
    merged_seq = fseq[: loffset]
    merged_qual = fqual[: loffset]

    flen = len(fseq)
    overl = flen - loffset

    # "Recover" overlapping region according to quality profile
    #   adopted from NGmerge.
    for i, j in zip(range(loffset, flen), range(overl)):

        fqual_i = ord(fqual[i]) - phred_offset
        rqual_j = ord(rqual[j]) - phred_offset
        fbase = fseq[i]
        rbase = rseq[j]

        if fbase == rbase:
            new_base = fbase
            new_qual = match_prof[fqual_i][rqual_j]
        else:
            if rqual_j > fqual_i:
                new_base = rbase
            else:
                new_base = fbase
            # end if

            new_qual = mismatch_prof[fqual_i][rqual_j]
        # end if

        merged_seq += new_base
        merged_qual += chr(new_qual + phred_offset)
    # end for

    # 3'-end comes from reverse read
    merged_seq += rseq[overl :]
    merged_qual += rqual[overl :]

    return (merged_seq, merged_qual)
# end def _merge_by_overlap
