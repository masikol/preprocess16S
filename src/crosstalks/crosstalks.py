# -*- coding: utf-8 -*-

import itertools


def detect_crosstalk(primers, fastq_records, threshold, max_offset):

    response_generator = (
        _find_primer(p, r.seq, threshold, max_offset) \
        for p, r in zip(primers, fastq_records)
    )

    valid_response_inerator = itertools.takewhile(
        _is_not_none,
        response_generator
    )

    valid_responses = tuple(valid_response_inerator)

    if len(valid_responses) < len(primers):
        return None
    else:
        return valid_responses
    # end if
# end def detect_crosstalk


def cut_off_primers(fastq_records, responses):
    for fastq_record, trim_rght_idx in zip(fastq_records, responses):
        fastq_record.seq = fastq_record.seq[trim_rght_idx :]
        fastq_record.qual_str = fastq_record.qual_str[trim_rght_idx :]
    # end for
# end def cut_off_primers


_MATCH_DICT = {
    'A': 'ARWMDHVN',
    'G': 'GRSKBDVN',
    'C': 'CYSMBHVN',
    'T': 'TUYWKBDHN',
    'U': 'TUYWKBDHN',
    'R': 'RDVN',
    'Y': 'YBHN',
    'S': 'SBVN',
    'W': 'WDHN',
    'K': 'KBDN',
    'M': 'MHVN',
    'B': 'BN',
    'D': 'DN',
    'H': 'HN',
    'V': 'VN',
    'N': ''
}


def _is_not_none(x):
    return not x is None
# end def _is_not_none


def _bases_are_equal(primer_base, read_base):
    return primer_base in _MATCH_DICT[read_base]
# end def _bases_are_equal


def _offset_generator(max_offset):
    yield 0
    for i in range(1, max_offset + 1):
        yield  i
        yield -i
    # end for
# end def _offset_generator


def _find_primer(primer, read, threshold=0.52, max_offset=2):

    for offset in _offset_generator(max_offset):

        primer_subseq = primer[max(0, -offset):]

        read_offset = max(0, offset)
        read_subseq_rindex = len(primer_subseq) + read_offset
        read_subseq = read[read_offset : read_subseq_rindex]

        score = sum(
            map(
                _bases_are_equal,
                primer_subseq, read_subseq
            )
        ) / len(primer_subseq)

        if score >= threshold:
            return read_subseq_rindex
        # end if
    # end for

    return None
# end def _find_primer
