# -*- coding: utf-8 -*-

import os
from math import log # for log scale in plot


try:
    import numpy as np
except ImportError as imperr:
    print_error("'numpy' package is not installed!")
    print( str(imperr) )
    print("Please install numpy (e.g. pip3 install numpy)")
    exit(1)
# end try

try:
    import matplotlib.pyplot as plt
except ImportError as imperr:
    print_error("'matplotlib' package is not installed!")
    print( str(imperr) )
    print("Please install matplotlib (e.g. pip3 install matplotlib)")
    exit(1)
# end try


# List of propabilities corresponding to indices (index is Q, value is the propability):
q2p_map = [10 ** (-q/10) for q in range(128)] # 127 -- max value of a signed byte
# Function for accessing propabilities by Q:
qual2prop = lambda q: q2p_map[q]
# Function for accessing Q by propability:
prop2qual = lambda p: round(-10 * log(p, 10), 2)

top_x_scale, step = 40.0, 0.5
# average read quality
X = np.arange(0, top_x_scale + step, step)
# amount of reads with sertain average quality
Y = np.zeros(int(top_x_scale / step), dtype=int)


# This function will be used as "organizer"
def calc_qual_disrib(fastq_recs, **kwargs):

    substr_phred_offs = kwargs["substr_phred_offs"]

    for rec in fastq_recs.values():

        qual_str = rec["qual_str"]

        qual_array = map(substr_phred_offs, qual_str)
        qual_array = tuple( map(qual2prop, qual_array) )

        avg_qual = prop2qual( np.mean(qual_array) )
        min_indx = ( np.abs(X - avg_qual) ).argmin()
        
        Y[min_indx] += 1
    # end for
# end def calc_qual_disrib


def create_plot(outdir_path, phred_offset):

    image_path = os.sep.join((outdir_path, "quality_plot.png"))

    fig, ax = plt.subplots()

    Y_buff = np.array(list(map( lambda y: log(y) if y > 0 else 0 , Y)))

    ax.plot(X[:len(Y_buff)], Y_buff, 'r')
    ax.set(title="Quality per read distribution",
        xlabel="avg quality (Phred{})".format(phred_offset),
        ylabel="ln(amount of reads)")
    ax.grid()
    fig.savefig(image_path)

    # Open image as image via xdg-open in order not to pause executing main script while
    #   you are watching at pretty plot.

    img_open_util = "xdg-open"
    util_found = False
    for directory in os.environ["PATH"].split(os.pathsep):
        if os.path.isdir(directory) and img_open_util in os.listdir(directory):
            util_found = True
            break
        # end if
    # end for

    if util_found:
        os.system("{} {}".format(img_open_util, image_path))
    # end if
    print("Plot image is here:\n  '{}'\n".format(image_path))

    return image_path
# end def create_plot