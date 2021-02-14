# -*- coding: utf-8 -*-

import os
import sys
import logging
import traceback as tb
from time import time, strftime, localtime, gmtime


def config_logging(outdir, version, last_update_date):
    logfile_path = os.path.join(outdir, 'preprocess16S.log')
    logging.basicConfig(filename=logfile_path,
        format='%(levelname)s: %(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
        level=logging.INFO, filemode='w')
    log_info(sys.platform)
    log_info(sys.implementation)
    log_info(sys.version)
    log_info('preprocess16S. Version {}. {} Edition.'.format(version, last_update_date))
# end def config_logging


start_time = time() # consider time of importing as start time

def getwt():
    # Function (get work time) returns time HH:MM:SS that has passed from start_time.
    return strftime('%H:%M:%S', gmtime( time() - start_time))
# end def getwt


def get_full_time():
    # Function returns current time HH:MM:SS (YYYY-mm-dd HH:MM:SS).
    # The first HH:MM:SS is from 'getwt'.
    return '{} ({})'\
        .format(getwt(), strftime('%Y-%m-%d %H:%M:%S', localtime(time())))
# end def get_full_time


def printlog_info_time(msg):
    sys.stdout.write('{} - {}\n'.format(getwt(), msg))
    logging.info('(%s from start)\t%s', getwt(), msg)
# end def printlog_info_time


def printlog_info(msg):
    sys.stdout.write('{}\n'.format(msg))
    logging.info('(%s from start)\t%s', getwt(), msg)
# end def printlog_info


def log_info(msg):
    logging.info('(%s from start)\t%s', getwt(), msg)
# end def log_info


def printlog_error(msg):
    sys.stderr.write('{}\n'.format(msg))
    logging.error('(%s from start)\t%s', getwt(), msg)
    logging.error('Stack trace:\n%s', ''.join(tb.format_stack()))
# end def printlog_error


def printlog_warning(msg):
    sys.stdout.write('{}\n'.format(msg))
    logging.warning('(%s from start)\t%s', getwt(), msg)
# end def printlog_warning
