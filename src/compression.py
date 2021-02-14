# -*- coding: utf-8 -*-

import os
import re
import sys
import bz2
import glob
import gzip
import shutil
import functools

from src.printlog import printlog_error, printlog_info, printlog_info_time
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


def _get_gzip_func():
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
# end def _get_gzip_func


def gzip_outfiles(outdir):
    gzip_func = _get_gzip_func()
    print()
    printlog_info_time('Gzipping output files...')

    is_fastq = lambda x: not re.match(r'.+\.f(ast)?q$', x) is None
    fq_fpaths = filter(is_fastq, glob.iglob(os.path.join(outdir, '*')))

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
