# -*- coding: utf-8 -*-
# Module contains a function that returns function called "crosstalk pipe".
# It's purpose is to not to check 'cut_off_primers_flag' each time for every damned
#   fastq record.

import src.crosstalks.crosstalks as crt

def provide_crosstalk_pipe(cut_off_primers_flag):
    # Purposeof this function is to not to check 'cut_off_primers_flag'
    #   each time for every damned fastq record.
    #
    # :param cut_off_primers_flag: flag indicating whether to cut off primers;
    # :type cut_off_primers_flag: bool;
    #
    # It returns a function for processing single pair of fastq records
    #   (or single fastq record for single-end reads).

    def cut_off_primers_pipe(primers, fastq_records, threshold, max_offset):
        # "Crosstalk pipe" which cuts primers off after crosstalk detection.
        #
        # :param primers: primers sequences;
        # :type primers: list<str>;
        # :param fastq_records: collection of fastq records to process;
        # :type fastq_records: list<FastqRecord>;
        # :param threshold: score threshold;
        # :type threshold: float;
        # :param max_offset: maxinum offset;
        # :type max_offset: int;
        #
        # Returns two values:
        # True/False: False if `fastq_records` is a crosstalks. Else True.
        # Pair of processed (maybe changed) `fastq_records`.
        cutoff_rght_idxs = crt.detect_crosstalk(primers, fastq_records, threshold, max_offset)
        if not cutoff_rght_idxs is None:
            crt.cut_off_primers(fastq_records, cutoff_rght_idxs)
            return (True, fastq_records)
        else:
            return (False, fastq_records)
        # end if
    # end def cut_off_primers_pipe

    def keep_primers_pipe(primers, fastq_records, threshold, max_offset):
        # "Crosstalk pipe" which keeps primers after crosstalk detection.
        #
        # :param primers: primers sequences;
        # :type primers: list<str>;
        # :param fastq_records: collection of fastq records to process;
        # :type fastq_records: list<FastqRecord>;
        # :param threshold: score threshold;
        # :type threshold: float;
        # :param max_offset: maxinum offset;
        # :type max_offset: int;
        #
        # Returns two values:
        # True/False: False if `fastq_records` is a crosstalks. Else True.
        # Pair of processed `fastq_records` (for campatibility with `cut_off_primers_pipe`).
        response = crt.detect_crosstalk(primers, fastq_records, threshold, max_offset)
        return (not response is None, fastq_records)
    # end def keep_primers_pipe

    return cut_off_primers_pipe if cut_off_primers_flag else keep_primers_pipe
# end def
