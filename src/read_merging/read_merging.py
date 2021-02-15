# -*- coding: utf-8 -*-

from src.read_merging.ngmerge_quality_profile import match_prof, mismatch_prof


def merge_read_pair(faref_report, raref_report, fseq, rseq, fqual, rqual, min_overlap, num_N, phred_offset):

    # Calculate some "features".
    # All these "features" are in coordinates of the reference sequence
    forw_start = int(faref_report.s_start) - int(faref_report.q_start)
    forw_end = forw_start + faref_report.qlen
    rev_start = raref_report.s_start - raref_report.q_start
    rev_end = rev_start + raref_report.qlen

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
        merged_seq, merged_qual = _merge_by_overlap(loffset, fseq, fqual, rseq, rqual, phred_offset)

        return 0, {"seq": merged_seq, "qual_str": merged_qual}

    # Here we have a gap.
    elif gap:
        
        gap_len = rev_start - forw_end
        # If gap is too long -- discard this pair.
        if gap_len > num_N:
            return 1, None
        
        else:
            # Fill in the gap with 'N'.
            merged_seq = fseq + 'N' * gap_len + rseq
            merged_qual = fqual + chr(phred_offset+3) * gap_len + rqual

            return 0, {"seq": merged_seq, "qual_str": merged_qual}
        # end if
    # end if

    # === Unforseen case. This code should not be ran. But if it does-- report about it. ===
    _handle_unforseen_case(f_id, fseq, r_id, rseq)
    return 3, None
# end def _try_merge


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


def _fill_gap(fseq, fqual, rseq, rqual, gap_len, phred_offset):
    merged_seq = fseq + 'N' * gap_len + rseq
    merged_qual = fqual + chr(phred_offset) * gap_len + rqual

    return (merged_seq, merged_qual)
# end def _fill_gap

