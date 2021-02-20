# -*- coding: utf-8 -*-
# Module contains functions for interacting with (maybe) compressed files.

import os
import re
import bz2
import glob
import gzip
import shutil
import functools

from src.platform import platf_depend_exit
from src.printlog import printlog_error, printlog_info, printlog_info_time


def provide_open_funcs(fpaths):
    # Function, which returns opening function(s) for input file(s).
    #
    # :param fpaths: collection of paths to input files;
    # :type fpaths: list<str>;

    open_funcs = list()

    try:
        for fpath in fpaths:
            # Check if input file is gzipped
            if _is_gzipped(fpath):
                open_funcs.append(
                    functools.partial(gzip.open, mode='rt', encoding='utf-8')
                )
            # Check if input file is bzipped2
            elif _is_bzipped(fpath):
                open_funcs.append(
                    functools.partial(bz2.open, mode='rt', encoding='utf-8')
                )
            # Check if input file is plain text file
            elif _is_plain_text(fpath):
                open_funcs.append(
                    functools.partial(open, mode='r', encoding='utf-8')
                )
            else:
                # Raise a super terrifying exception
                raise _InvalidFileError('Error: cannot read file `{}`: \
it is neither plain text file, nor gzipped, nor bzipped2.'.format(fpath))
            # end if
        # end for
    except _InvalidFileError as err:
        printlog_error(str(err))
        platf_depend_exit(1)
    # end try

    return open_funcs
# end def get_compression_tools


def gzip_outfiles(outdir):
    # Function gzips all fastq files in directory `outdir`.
    #
    # :param outdir: path to outdir;
    # :type outdir: str;

    # Get gzipping function
    gzip_func = _get_gzip_func()
    print()
    printlog_info_time('Gzipping output files...')

    # Get fastq files
    is_fastq = lambda x: not re.match(r'.+\.f(ast)?q$', x) is None
    fq_fpaths = filter(is_fastq, glob.iglob(os.path.join(outdir, '*')))

    # Gzip them!
    for fpath in fq_fpaths:
        try:
            gzip_func(fpath)
        except OSError as err:
            printlog_info('Error: cannot gzip file `{}`: {}.'.format(fpath, err))
            platf_depend_exit(1)
        # end try
    # end for

    printlog_info_time('Output files are gzipped.')
# end def gzip_outfiles


def _get_gzip_func():
    # Function returns gzipping function depending on whether GNU gzip utility is available.
    # It not, the program will compress files ot it's own. Oh yeah, it will be slow.

    # Search for GNU gzip utility in PATH
    gzip_util = 'gzip'
    util_found = False
    for directory in os.environ['PATH'].split(os.pathsep):
        if os.path.isdir(directory) and gzip_util in os.listdir(directory):
            util_found = True
            break
        # end if
    # end for

    def gzip_with_gnu_gzip(fpath):
        # Funtion for gzipping with GNU gzip.
        printlog_info('Gzipping `{}`'.format(fpath))
        os.system('{} {}'.format(gzip_util, fpath))
    # end def gzip_with_gnu_gzip

    def gzip_with_shutil(fpath):
        # Funtion for gzipping using Python funtionality.
        printlog_info('Gzipping `{}`'.format(fpath))
        with open(fpath, 'rb') as plain_file, gzip.open(fpath+'.gz', 'wb') as gz_file:
            shutil.copyfileobj(plain_file, gz_file)
        # end with
        os.unlink(fpath)
    # end def gzip_with_shutil

    # Return appropriate function
    if util_found:
        return gzip_with_gnu_gzip
    else:
        return gzip_with_shutil
    # end if
# end def _get_gzip_func


def _is_gzipped(fpath):
    # Function, which returs True if passed file is a file compressed with gzip.
    # :param fpath: path to file ot check;
    # :type fpath: str;

    if not fpath.endswith('.gz'):
        return False
    # end if

    # Test this file
    try:
        with gzip.open(fpath, 'rt', encoding='utf-8') as gzfile:
            gzfile.readline() # this will fail if file is not valid
        # end with
    except (gzip.BadGzipFile, OSError, UnicodeDecodeError) as err:
        raise _InvalidFileError('Error: file `{}` is not a valid gzip file: {}'\
            .format(fpath, err))
    else:
        return True
    # end try
# end def _is_gzipped


def _is_bzipped(fpath):
    # Function, which returs True if passed file is a file compressed with bzip2.
    # :param fpath: path to file ot check;
    # :type fpath: str;

    if not fpath.endswith('.bz2'):
        return False
    # end if

    # Test this file
    try:
        with bz2.open(fpath, 'rt', encoding='utf-8') as bz2file:
            bz2file.readline() # this will fail if file is not valid
        # end with
    except (OSError, UnicodeDecodeError) as err:
        raise _InvalidFileError('Error: file `{}` is not a valid bzip2 file: {}'\
            .format(fpath, err))
    else:
        return True
    # end try
# end def _is_bzipped


def _is_plain_text(fpath):
    # Function, which returs True if passed file is a plain text file.
    # :param fpath: path to file ot check;
    # :type fpath: str;

    # Test this file
    try:
        with open(fpath, 'r', encoding='utf-8') as file:
            file.readline() # this will fail if file is not valid
        # end with
    except (OSError, UnicodeDecodeError) as err:
        raise _InvalidFileError('Error: cannot read file `{}`: {}.'\
            .format(fpath, err))
    else:
        return True
    # end try
# end def _is_plain_text


class _InvalidFileError(OSError):
    # Custom exception for indicating that file compression type is not supported.
    pass
# end class _InvalidFileError
