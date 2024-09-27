# -*- coding: utf-8 -*-
# Module contains functions for detecting crosstalks.

import itertools


def detect_crosstalk(primers, fastq_records, threshold, max_offset):
    # Function determines if read pair (or a single read) is a cross-talk.
    #
    # :param primers: primser sequence(s);
    # :type primers: list<str>;
    # :param fastq_records: list of read(s);
    # :type fastq_records: list<FastqRecord>;
    # :param threshold: score threshold;
    # :type threshold: float;
    # :param max_offset: maximum offset;
    # :type max_offset: int;
    #
    # Returns None if `fastq_records` is as crosstalk.
    # Returns tuple of rightmost primer coordinates (0-based) in a read.
    # -CCTACGGGNGGCWGCAG
    # CCTACGGGAGGCTGCAGTGGGGAATATTGCACAA...
    #                  ^
    #                  this coordinate

    # Create generator that yields responses of '_find_primer' function
    #   first for forward read and then -- for reverse one.
    response_generator = (
        _find_primer(p, r.seq, threshold, max_offset) \
        for p, r in zip(primers, fastq_records)
    )

    # If primer is detected in forward read -- check reverse one.
    # If no -- exit immediately and don't even check reverse one.
    valid_response_iterator = itertools.takewhile(
        _is_not_none,
        response_generator
    )

    # Run generator and get responses
    valid_responses = tuple(valid_response_iterator)

    if len(valid_responses) < len(primers):
        return None # two valid responses (or one for single read)
    else:
        return valid_responses # < 2 valid responses
    # end if
# end def detect_crosstalk


def cut_off_primers(fastq_records, responses):
    # Function cuts off primers sequence from FASTQ record(s)
    # :param fastq_records: list of fastq records to process;
    # :type fastq_records: list<FastqRecord>;
    # :param responses: list of rightmost primer positions in a read;
    # :type responses: list<int>, tuple<int>;
    for fastq_record, trim_rght_idx in zip(fastq_records, responses):
        # Cut off from sequence
        fastq_record.seq = fastq_record.seq[trim_rght_idx :]
        # Cut off from quality string
        fastq_record.qual_str = fastq_record.qual_str[trim_rght_idx :]
    # end for
# end def cut_off_primers


# Matching dictionary for defining 'equals' operator
#   for bases, considering degenerate bases.
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
#    ^   ^
# read   primer
# base   base


def _is_not_none(x):
    return not x is None
# end def _is_not_none


def _bases_are_equal(primer_base, read_base):
    # Function is an 'equals' operator for bases, considering degenerate bases.
    # :type primer_base: str;
    # :type read_base: str;
    return primer_base in _MATCH_DICT[read_base]
# end def _bases_are_equal


def _offset_generator(max_offset):
    # Function yields offsets in following way:
    #   0,1,-1,2,-2...max_offset,-max_offset
    # :param max_offset: maximum offset;
    # :type max_offset: int;
    yield 0
    for i in range(1, max_offset + 1):
        yield  i
        yield -i
    # end for
# end def _offset_generator


def _find_primer(primer, read, threshold, max_offset):
    # Function seeks for primer in a read sequence.
    #
    # :param primer: primer sequence;
    # :type primer: str;
    # :param read: read sequence;
    # :type read: str;
    # :param threshold: score threshold;
    # :type threshold: float;
    # :param max_offset: maximum offset;
    # :type max_offset: int;

    for offset in _offset_generator(max_offset):

        # Get primer sequence to compare
        #   (it is a shorter substring of primer sequence
        #   if offset is negative)
        primer_subseq = primer[max(0, -offset):]
        primer_subseq_len = len(primer_subseq)

        # Get read sequence to compare
        read_offset = max(0, offset)
        read_subseq_rindex = primer_subseq_len + read_offset
        read_subseq = read[read_offset : read_subseq_rindex]

        # Calculate score
        score = sum(
            map(
                _bases_are_equal,
                primer_subseq, read_subseq
            )
        ) / primer_subseq_len

        # Primer is found -- return rightmost coordinate of primer sequence in a read
        if score >= threshold:
            return read_subseq_rindex
        # end if
    # end for

    return None
# end def _find_primer
