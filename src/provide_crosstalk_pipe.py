# -*- coding: utf-8 -*-

import src.crosstalks as crt

def provide_crosstalk_pipe(cut_off_primers_flag):

    def cut_off_primers_pipe(primers, fastq_records, threshold, max_offset):
        cutoff_rght_idxs = crt.detect_crosstalk(primers, fastq_records, threshold, max_offset)
        if not cutoff_rght_idxs is None:
            crt.cut_off_primers(fastq_records, cutoff_rght_idxs)
            return (True, fastq_records)
        else:
            return (False, fastq_records)
        # end if
    # end def cut_off_primers_pipe

    def keep_primers_pipe(primers, fastq_records, threshold, max_offset):
        response = crt.detect_crosstalk(primers, fastq_records, threshold, max_offset)
        return (not response is None, fastq_records)
    # end def keep_primers_pipe

    return cut_off_primers_pipe if cut_off_primers_flag else keep_primers_pipe
# end def
