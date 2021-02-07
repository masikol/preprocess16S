# -*- coding: utf-8 -*-

class CrosstalksArguments:

    def __init__(self, infpaths, primers, threshold, max_offset, cut_off_primers, outdir):
        self.infpaths = infpaths
        self.primers = primers
        self.threshold = threshold
        self.max_offset = max_offset
        self.cut_off_primers = cut_off_primers
        self.outdir = outdir
    # end def __init__
# end class CrosstalksArguments
