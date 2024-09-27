# -*- coding: utf-8 -*-
# Mudule contains container classes for storing arguments for preprocess16S tasks.

class CrosstalksArguments:
    # Container class for storing arguments for crosstalks detection task.
    # Fields:
    # infpath: collection of paths to input files (list<str>).
    # primers: collection of primer sequences (list<str>).
    # threshold: score threshold (float).
    # max_offset: maximum offset (int).
    # cut_off_primers: flag for cutting off primers (bool).
    # outdir: path to output directory (str).

    def __init__(self, infpaths, primers, threshold, max_offset, cut_off_primers, outdir):
        self.infpaths = infpaths
        self.primers = primers
        self.threshold = threshold
        self.max_offset = max_offset
        self.cut_off_primers = cut_off_primers
        self.outdir = outdir
    # end def __init__
# end class CrosstalksArguments


class NGmergeArguments:
    # Container class for storing arguments for NGmerge task.
    # Fields:
    # infpath: collection of paths to input files (list<str>).
    # ngmerge: path to NGmerge executable (str).
    # n_thr: number of threads (int).
    # min_overlap: minimum overlap (int).
    # mismatch_frac: maximum mismatch fraction in overlappint region (float).
    # phred_offset: Phred offset (int).
    # outdir: path to output directory (str).

    def __init__(self, infpaths, ngmerge, n_thr, min_overlap, mismatch_frac, phred_offset, outdir):
        self.infpaths = infpaths
        self.ngmerge = ngmerge
        self.n_thr = n_thr
        self.min_overlap = min_overlap
        self.mismatch_frac = mismatch_frac
        self.phred_offset = phred_offset
        self.outdir = outdir
    # end def __init__
# end class NGMergeArguments
