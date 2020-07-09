# -*- coding: utf-8 -*-
# __last_update_date__ = "2020-08-07"

import os

import multiprocessing as mp
from functools import reduce
from math import log # for log scale in plot

from src.filesystem import *
from src.fastq import *
from src.printing import *


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


# In 2020 sequencators do not read better that Q40.
top_x_scale, step = 40.0, 0.5
# average read quality
X = np.arange(0, top_x_scale + step, step)


def fastq_read_packets(read_paths, reads_at_all, n_thr):
    """
    Function-generator for retrieving FASTQ records from PE files
        and distributing them evenly between 'n_thr' processes
        for further parallel processing.
        # end for

    :param read_paths: dictionary (dict<str: str> of the following structure:
    {
        "R1": path_to_file_with_forward_reads,
        "R2": path_to_file_with_reverse_reads
    }
    :param reads_at_all: number of read pairs in these files;
    :type reads_at_all: int;

    Yields lists of FASTQ-records (structure of these records is described in 'write_fastq_record' function).
    Returns None when end of file(s) is reached.
    """
    how_to_open = OPEN_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]
    fmt_func = FORMATTING_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]
    read_files = dict() # dictionary for file objects

    # Open files that contain reads meant to be processed.
    for key, path in read_paths.items():
        read_files[key] = how_to_open(path)
    # end for

    # Compute packet size (one packet -- one thread).
    pack_size = reads_at_all // n_thr
    if reads_at_all % n_thr != 0:
        pack_size += 1
    # end if

    for i in range(n_thr):
        packet = list()
        for j in range(pack_size):
            tmp = read_fastq_pair(read_files, fmt_func)
            # if the end of the line is reached
            if tmp is None:
                # Yield partial packet and return 'None' on next call
                yield packet
                return
            else:
                packet.append(tmp)
            # end if
        # end for
        yield packet # yield full packet
    # end for
# end def fastq_read_packets


def single_qual_calcer(data, reads_at_all, substr_phred_offs):
    """
    Function that performs task meant to be done by one process while parallel quality calculation.

    :param data: list of FASTQ-records
        (structure of these records is described in 'write_fastq_record' function);
    :type data: list< dict<str: str> >;
    :param reads_at_all: total number of read pairs in input files;
    :type reads_at_all: int;

    Returns numpy.ndarray<int> performing quality distribution of reads.
    """

    # amount of reads with sertain average quality
    Y = np.zeros(int(top_x_scale / step), dtype=int)

    # Processes will print number of processed reads every 'delay' reads.
    delay, i = 1000, 0

    for fastq_recs in data:

        for rec in fastq_recs.values():

            qual_str = rec["qual_str"]

            qual_array = map(substr_phred_offs, qual_str)
            qual_array = tuple( map(qual2prop, qual_array) )

            avg_qual = prop2qual( np.mean(qual_array) )
            min_indx = ( np.abs(X - avg_qual) ).argmin()

            Y[min_indx] += 1
        # end for

        # if next 'delay' reads are processed
        if i == delay:

            # synchronized incrementing
            with count_lock:
                counter.value += delay
            # end with

            bar_len = 50 # length of status bar
            eqs = int( bar_len*(counter.value / reads_at_all) ) # number of '=' characters
            spaces = bar_len - eqs # number of ' ' characters

            with print_lock:
                printn("\r[" + "="*eqs + '>' + ' '*spaces +"] {}% ({}/{})".format(int(counter.value/reads_at_all*100),
                    counter.value, reads_at_all))
            # end with

            i = 0 # reset i
        # end if

        i += 1
    # end for
    return Y
# end def single_qual_calcer


def proc_init(print_lock_buff, counter_buff, count_lock_buff):
    """
    Function that initializes global variables that all processes shoud have access to.
    This function is meant to be passed as 'initializer' argument to 'multiprocessing.Pool' function.

    :param print_lock_buff: lock that synchronizes printing to the console;
    :type print_lock_buff: multiprocessing.Lock;
    :param counter_buff: integer number representing number of processed reads;
    :type counter_buff: multiprocessing.Value;
    :param count_lock_buff: lock that synchronizes incrementing 'counter' variable;
    :type count_lock_buff: multiprocessing.Lock;
    """

    global print_lock
    print_lock = print_lock_buff

    global counter
    counter = counter_buff

    global count_lock
    count_lock = count_lock_buff
# end def proc_init



def parallel_qual(read_paths, n_thr, substr_phred_offs):
    """
    Function launches parallel quality calculations.

    :param read_paths: dict of paths to read files. it's structure is described in 'fastq_read_packets' function;
    :type read_paths: dict<str: str>;
    :param n_thr: int;
    :type n_thr: int:
    """
    
    print_lock = mp.Lock()     # lock that synchronizes printing to console;
    count_lock = mp.Lock()     # lock that synchronizes incrementing 'counter';
    counter = mp.Value('i', 0) # integer number representing number of processed reads;

    how_to_open = OPEN_FUNCS[ get_archv_fmt_indx(read_paths["R1"]) ]

    # Count number of read pairs
    reads_at_all = int(sum(1 for line in how_to_open(read_paths["R1"])) / 4)

    # Create pool of processes
    pool = mp.Pool(n_thr, initializer=proc_init, initargs=(print_lock, counter, count_lock))
    # Run parallel calculations
    Y = pool.starmap(single_qual_calcer, [(data, reads_at_all, substr_phred_offs) for data in fastq_read_packets(read_paths, reads_at_all, n_thr)])
    print("\r["+"="*50+"] 100% ({}/{})\n".format(reads_at_all, reads_at_all))

    # Reaping zombies
    pool.close()
    pool.join()

    def arr_sum(Y1, Y2): # function to perform 'functools.reduce' sum
        return Y1 + Y2
    # end def arr_sum

    return reduce(arr_sum, Y) # element-wise sum of all processes' results
# end def parallel_qual


def create_plot(Y, outdir_path, phred_offset):

    image_path = os.sep.join((outdir_path, "quality_plot.png"))

    fig, ax = plt.subplots()

    Y = np.array(list(map( lambda y: log(y) if y > 0 else 0 , Y)))

    ax.plot(X[:len(Y)], Y, 'r')
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