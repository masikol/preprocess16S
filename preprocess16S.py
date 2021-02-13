#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = '5.0.a'
# Year, month, day
__last_update_date__ = '2021-XX-XX'
__author__ = 'Maxim Sikolenko'


import sys

if sys.version_info.major < 3:
    print( 'Your python interpreter version is ' + '%d.%d' % (sys.version_info.major,
        sys.version_info.minor) )
    print('   Please, use Python 3.\a')
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    if sys.platform.startswith('win'):
        raw_input('Press ENTER to exit:')
    # end if
    sys.exit(1)
# end if


import os
import re
import glob
import getopt
import itertools

import src.runners
import src.compression
from src.platform import platf_depend_exit
import src.arguments as srcargs
from src.crosstalks.get_primers_seqs import get_primers_seqs
from src.printlog import config_logging, printlog_info_time, printlog_info, getwt


if '-v' in sys.argv[1:] or '--version' in sys.argv[1:]:
    print(__version__)
    platf_depend_exit(0)
# end if

if '-h' in sys.argv[1:] or '--help' in sys.argv[1:]:
    print(' <<< preprocess16S >>>')
    print(' Version {}. {} edition.'.format(__version__, __last_update_date__))
    platf_depend_exit(0)
# end if


def handle_args():

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:],
            'hv1:2:o:t:r:x:s:c:m:p:',
            ['help', 'version',
            'tasks=', 'R1=', 'R2=', 'outdir=', 'n-thr=',
            'primers=', 'threshold=', 'max-offset=', 'cut-off-primers=',
            'ngmerge-path=', 'min-overlap=', 'mismatch-frac='])
    except getopt.GetoptError as err:
        print( str(err) )
        platf_depend_exit(2)
    # end try

    args = {
        # general opts
        'tasks': ['rm-crosstalks'],
        '1': None,                                              # forward reads
        '2': None,                                              # reverse reads
        'o': os.path.join(os.getcwd(), 'preprocess16S-outdir'), # outdir
        't': 1,                                                 # threads number
        # crosstalks opts
        'r': None,                                              # primers file
        'x': 0.52,                                              # threshold
        's': 2,                                                 # max-offset
        'c': True,                                              # cut-off-primers
        # read merging opts
        'ngmerge-path': os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'binaries', 'NGmerge'
        ),                                                      # ngmerge path
        'm': 20,                                                # min-overlap
        'p': 0.1,                                               # mismatch-frac
    }

    for opt, arg in opts:

        # General opts
        if opt == '--tasks':
            args['tasks'] = set(arg.split(','))
            permitted_tasks = ('rm-crosstalks', 'ngmerge')
            isnot_permitted_task = lambda x: not x in permitted_tasks
            if any(map(isnot_permitted_task, args['tasks'])):
                print('Error: task(s) not recognized: {}.'\
                    .format(', '.join(filter(isnot_permitted_task, args['tasks']))
                    )
                )
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end if

        elif opt in ('-1', '--R1'):
            if os.path.exists(arg):
                args['1'] = arg
            else:
                print('\aError: file `{}` does not exist.'.format(arg))
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end if

        elif opt in ('-2', '--R2'):
            if os.path.exists(arg):
                args['2'] = arg
            else:
                print('\aError: file `{}` does not exist.'.format(arg))
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end if

        elif opt in ('-o', '--outdir'):
            if os.path.isfile(arg):
                print('Error: `{}` cannot be output directory. \
Regular file with the same name already exists.'.format(arg))
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            else:
                args['p'] = arg
            # end if

        elif opt in ('-t', '--n-thr'):
            try:
                args['t'] = int(arg)
                if args['t'] < 0:
                    raise ValueError
                # end if
            except ValueError:
                print('Invalid threads number value: `{}`.'.format(arg))
                print('It must be integer > 0.')
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end try
            max_threads = len(os.sched_getaffinity(0))
            if args['t'] > max_threads:
                print('Warning! You have specified {} threads to use \
although {} are available.'.format(args['t'], max_threads))
                print('Switching threads number to {}.\n'.format(max_threads))
                args['t'] = max_threads
            # end if

        # Crosstalks opts
        elif opt in ('-r', '--primers'):
            if os.path.exists(arg):
                args['p'] = arg
            else:
                print('\aError: file `{}` does not exist.'.format(arg))
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end if


        elif opt in ('-x', '--threshold'):
            try:
                args['x'] = float(arg)
                if args['x'] < 0 or args['x'] > 1:
                    raise ValueError
                # end if
            except ValueError:
                print('Invalid threshold: `{}`.'.format(arg))
                print('It must be a number in following interval: [0, 1].')
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end try

        elif opt in ('-s', '--max-offset'):
            try:
                args['s'] = int(arg)
                if args['s'] < 0:
                    raise ValueError
                # end if
            except ValueError:
                print('Invalid max offset value: `{}`.'.format(arg))
                print('It must be integer > 0.')
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end try

        elif opt in ('-c', '--cut-off-primers'):
            if arg == '0':
                args['c'] = False
            elif arg == '1':
                args['c'] = True
            else:
                print('Error: invalid value passed with option {}: `{}`.'\
                    .format(opt, arg))
                print('It must be 0 or 1.')
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end if

        # Read merging opts
        elif opt == '--ngmerge-path':
            if not os.path.exists(arg):
                print('\aError: file `{}` does not exist.'.format(arg))
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end if
            if not os.access(arg, os.X_OK):
                print('\aError: NGmerge file is not executable: `{}`'.format(arg))
                print('Please, make it executable. You can do it in this way:')
                print(' chmod +x {}'.format(ngmerge))
            # end if
            args['ngmerge-path'] = arg

        elif opt in ('-m', '--min-overlap'):
            try:
                args['m'] = int(arg)
                if args['m'] < 0:
                    raise ValueError
                # end if
            except ValueError:
                print('Invalid min overlap value: `{}`.'.format(arg))
                print('It must be integer > 0.')
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end try

        elif opt in ('-p', '--mismatch-frac'):
            try:
                args['p'] = float(arg)
                if args['p'] < 0 or args['p'] > 1:
                    raise ValueError
                # end if
            except ValueError:
                print('Invalid mismatch fraction: `{}`.'.format(arg))
                print('It must be a number in following interval: [0, 1].')
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end try
        # end if
    # end for

    mandatory_opts = ['1']
    if 'ngmerge' in args['tasks']:
        mandatory_opts.append('2')
    # end if

    for opt in mandatory_opts:
        if args[opt] is None:
            print('Error: option `-{}` is mandatory.'.format(opt))
            print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
            platf_depend_exit(0)
        # end if
    # end for

    # end if
    return args
# end def handle_args


def make_outdir(outdir):

    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except OSError as err:
            print('Cannot create output directory: {}'.format(err))
            platf_depend_exit(1)
        # end try

    elif len(os.listdir(outdir)) != 0:

        print('\nOutput directory `{}` is not empty.'.format(outdir))

        error = True
        while error:
            reply = input("""Press ENTER to remove all files in it or
 enter `q` to exit\n >> """)

            if reply == '':
                error = False
                for fpath in glob.iglob(os.path.join(outdir, '*')):
                    print('Removing `{}`'.format(fpath))
                    try:
                        os.unlink(fpath)
                    except OSError as err:
                        print('Error. Cannot remove file `{}`: {}'\
                            .format(fpath, err))
                    # end try
                # end for
                print()
            elif reply.lower() == 'q':
                sys.exit(0)
            else:
                print('Invalid reply: `{}`'.format(reply))
            # end if
        # end while
# end def make_outdir


def gzip_outfiles(outdir):
    gzip_func = src.compression.get_gzip_func()
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


def run_task_chain(args):

    curr_infpaths = tuple(itertools.takewhile(
                        lambda x: not x is None,
                        [args['1'], args['2']])
    )
    tasks_sorted = sorted(args['tasks'], reverse=True)

    for task in tasks_sorted:

        if task == 'rm-crosstalks':
            arguments = srcargs.CrosstalksArguments(
                curr_infpaths,                # infpaths
                get_primers_seqs(args['r']),  # primers
                args['x'],                    # threshold
                args['s'],                    # max-offset
                args['c'],                    # cut-off-primers
                args['o']                     # outdir
            )
            curr_infpaths = src.runners.crosstalks_runner(arguments)

        elif task == 'ngmerge':
            arguments = srcargs.NGmergeArguments(
                curr_infpaths,                # infpaths
                args['ngmerge-path'],         # ngmerge
                args['t'],                    # n_thr
                args['m'],                    # min_overlap
                args['p'],                    # mismatch_frac
                args['o']                     # outdir
            )
            curr_infpaths = src.runners.ngmerge_runner(arguments)
    # end if

# end def run_task_chain


def main():
    args = handle_args()
    make_outdir(args['o'])
    config_logging(args['o'])
    run_task_chain(args)
    gzip_outfiles(args['o'])
    printlog_info_time('Work is completed.')
# end def main


if __name__ == '__main__':
    main()
# end if
