# -*- coding: utf-8 -*-
# Module for dealing with filesystem: close/open files, detect and handlt files in different formats.
# __last_update_date__ = "2020-08-07"

import sys
from gzip import open as open_as_gzip
from gzip import GzipFile
from bz2 import open as open_as_bz2
from _io import TextIOWrapper

def get_archv_fmt_indx(fpath):
    if fpath.endswith(".gz"):
        return 1
    elif fpath.endswith(".bz2"):
        return 2
    else:
        return 0
    # end if
# end def get_archv_fmt_indx


# Following tuples of functions are here in order to process '.fastq', '.fastq.gz' and '.fastq.bz2'
#    files using the same interface.

OPEN_FUNCS = (open, open_as_gzip, open_as_bz2)

FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format .fastq line
    lambda line: line.decode("utf-8").strip(),   # format .fastq.gz line
    lambda line: line.decode("utf-8").strip(),   # format .fastq.bz2 line (same as previous, but for consistency let it be so)
)


# |========================= Functions =========================|

def close_files(*files):
    """
    Function closes all files passed to it, no matter whether they are single file objects or collenctions.

    :type of the argumets:
        'dict<str: _io.TextIOWrapper' or 'gzip.GzipFile>'' or
        'list<_io.TextIOWrapper or gzip.GzipFile>' or
        'tuple<_io.TextIOWrapper or gzip.GzipFile>' or
        '_io.TextIOWrapper' or
        'gzip.GzipFile'
    """
    for obj in files:

        if isinstance(obj, dict) or isinstance(obj, list) or isinstance(obj, tuple):
            for indx_or_key in obj:
                obj[indx_or_key].close()
            # end for

        elif isinstance(obj, TextIOWrapper) or isinstance(obj, GzipFile):
            obj.close()

        else:
            print("""If you use function 'close_files', please, store file objects in 
    lists, tuples or dictionaries or pass file objects itself to the function.""")
            print("You have passed '{}'".format( type(obj) ))
            try:
                obj.close()
            
            except:
                print_error("Object passed to 'close_files' function can't be closed")
                sys.exit(1)
            # end try
        # end if
    # end for
# end def close_files


def open_files(fpaths, how_to_open, mode=None):
    """
    Function opens files. Path to files should be passed to it in a dictionary.

    :param fpaths: dictionary of path to files;
    :type fpaths: dict<str: str>;
    :param how_to_open: file-opening function;
    :type how_to_open: function;
    :param mode: file-opening mode (same as in 'open' function);
    :type mode: str;

    Returns dictionary <key: file_object>
    """

    files = dict()

    try:
        for key in fpaths.keys():
            if mode is None:
                files[key] = how_to_open(fpaths[key]) # for uniform opening plain and compressed files for reading
            else:
                files[key] = how_to_open(fpaths[key], mode)
            # end if
        # end for
    except OSError as oserror:
        print("Error while opening file", str(oserror))
        close_files(files)
        sys.exit(1)
    # end try

    return files
# end def open_files