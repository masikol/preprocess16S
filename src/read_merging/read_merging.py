# -*- coding: utf-8 -*-

from src.read_merging.ngmerge_quality_profile import match_prof, mismatch_prof
import src.read_merging.blastoutfmt as bof


def merge_read_pair(faref_report, raref_report, fseq, rseq, fqual, rqual, max_gap, phred_offset):

    # Calculate some "features".
    # All these "features" are in coordinates of the reference sequence
    forw_start = faref_report[bof.SSTART] - faref_report[bof.QSTART]
    forw_end = forw_start + faref_report[bof.QLEN]
    revr_start = raref_report[bof.SSTART] - raref_report[bof.QSTART]
    revr_end = revr_start + raref_report[bof.QLEN]

    if forw_end > revr_end:
        # Dovetailed read pair
        return 1, None
    # end if

    # ===  Handle alignment results ====

    if not forw_end < revr_start:
        # No gap -- reads do overlap
        forw_offset = revr_start - forw_start
        merged_seq, merged_qual = _merge_by_overlap(forw_offset, fseq, fqual, rseq, rqual, phred_offset)

        return 0, (merged_seq, merged_qual)
    else:
        # Here we have a gap.
        gap_len = revr_start - forw_end
        # If gap is too long -- discard this pair.
        if gap_len > max_gap:
            return 1, None
        else:
            # Fill in the gap with 'N's.
            merged_seq, merged_qual = _fill_gap(fseq, fqual, rseq, rqual, gap_len, phred_offset)
            return 0, (merged_seq, merged_qual)
        # end if
    # end if

    # === Unforseen case. This code should not be ran. But if it does-- report about it. ===
    # _handle_unforseen_case(f_id, fseq, r_id, rseq)
    return 3, None
# end def merge_read_pair


def _merge_by_overlap(forw_offset, fseq, fqual, rseq, rqual, phred_offset):
    # Function merges two reads according to their overlapping redion.
    # Function leaves nucleotide with higher quality.
    # :param forw_offset: number of nucleotides in forward read from it's beginning to start of the overlapping region;
    # :type forw_offset: int;
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
    merged_seq = fseq[: forw_offset]
    merged_qual = fqual[: forw_offset]

    flen = len(fseq)
    overl = flen - forw_offset

    # "Recover" overlapping region according to quality profile
    #   adopted from NGmerge.
    for i, j in zip(range(forw_offset, flen), range(overl)):

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

