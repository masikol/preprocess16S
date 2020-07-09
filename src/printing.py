# -*- coding: utf-8 -*-
# Module defines some functions for printing,
# __last_update_date__ = "2020-08-07"

import sys

# |===== Stuff for dealing with time =====|

from time import time, strftime, localtime, gmtime
start_time = time()


def get_work_time():
    return strftime("%H:%M:%S", gmtime( time() - start_time))
# end def get_work_time


def printn(text=""):
    """
    Function prints text to the console without adding '\n' in the end of the line.
    Why not just to use 'print(text, end="")'?
    In order to display informative error message if Python 2.X is launched
        instead if awful error traceback.
    """
    sys.stdout.write(text)
# end def printn


def print_error(text):
    """Function for printing pretty error messages"""
    print("\n   \a!! - ERROR: " + text + '\n')
# end def print_error