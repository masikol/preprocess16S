# -*- coding: utf-8 -*-

import os
import sys
import bz2
import gzip
import shutil
import functools

from src.printlog import printlog_error, printlog_info
from src.platform import platf_depend_exit


def provide_open_funcs(fpaths):

    open_funcs = list()

    try:
        for fpath in fpaths:
            if _is_gzipped(fpath):
                open_funcs.append(
                    functools.partial(gzip.open, mode='rt', encoding='utf-8')
                )
            elif _is_bzipped(fpath):
                open_funcs.append(
                    functools.partial(bz2.open, mode='rt', encoding='utf-8')
                )
            elif _is_plain_text(fpath):
                open_funcs.append(
                    functools.partial(open, mode='r', encoding='utf-8')
                )
            else:
                raise _InvalidFileError('Error: cannot read file `{}`: {}.'\
                    .format(fpath, err))
            # end if
        # end for
    except _InvalidFileError as err:
        printlog_error(str(err))
        platf_depend_exit(1)
    # end try

    return open_funcs
# end def get_compression_tools


def get_gzip_func():
    # Gzip result files
    gzip_util = 'gzip'
    util_found = False
    for directory in os.environ['PATH'].split(os.pathsep):
        if os.path.isdir(directory) and gzip_util in os.listdir(directory):
            util_found = True
            break
        # end if
    # end for

    def gzip_with_gnu_gzip(fpath):
        printlog_info('Gzipping `{}`'.format(fpath))
        os.system('{} {}'.format(gzip_util, fpath))
    # end def gzip_with_gnu_gzip

    def gzip_with_shutil(fpath):
        printlog_info('Gzipping `{}`'.format(fpath))
        with open(fpath, 'rb') as plain_file, gzip.open(fpath+'.gz', 'wb') as gz_file:
            shutil.copyfileobj(plain_file, gz_file)
        # end with
        os.unlink(fpath)
    # end def gzip_with_shutil

    if util_found:
        return gzip_with_gnu_gzip
    else:
        return gzip_with_shutil
    # end if
# end def get_gzip_func


def _is_gzipped(fpath):

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
    pass
# end class _InvalidFileError
