#!/usr/bin/env python3

__version__ = '5.0.h'
# Year, month, day
__last_update_date__ = '2024-09-27'
__author__ = 'Maksim Sikolenko'


import sys

# Check Python version.
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
import getopt
import itertools

import src.runners
import src.print_help
import src.filesystem
import src.compression
import src.arguments as srcargs
from src.platform import platf_depend_exit
from src.crosstalks.get_primers_seqs import get_primers_seqs
from src.printlog import config_logging, printlog_info_time, printlog_info


# Print-and-exit
if '-v' in sys.argv[1:] or '--version' in sys.argv[1:]:
    print(__version__)
    platf_depend_exit(0)
# end if

if '-h' in sys.argv[1:] or '--help' in sys.argv[1:]:
    src.print_help.print_help(__version__, __last_update_date__)
    platf_depend_exit(0)
# end if

# preprocess16S tasks
TASKS = ('rm-crosstalks', 'ngmerge')


def handle_args():
    # Function handles command-line arguments.
    # Returns dictionary of parameters.

    # Parse arguments
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:],
            'hv1:2:o:t:z:r:x:s:c:m:p:q:',
            ['help', 'version',
            'tasks=', 'R1=', 'R2=', 'outdir=', 'threads=', 'gzip-output=',
            'primers=', 'threshold=', 'max-offset=', 'cut-off-primers=',
            'ngmerge-path=', 'min-overlap=', 'mismatch-frac=', 'phred-offset='])
    except getopt.GetoptError as err:
        print( str(err) )
        platf_depend_exit(2)
    # end try

    # Arguments dictionary with default values.
    args = {
        # general opts
        'tasks': ['rm-crosstalks'],
        '1': None,                                              # forward reads
        '2': None,                                              # reverse reads
        'o': os.path.join(os.getcwd(), 'preprocess16S-outdir'), # outdir
        't': 1,                                                 # threads number
        'z': True,                                              # gzip output
        # crosstalks opts
        'r': None,                                              # primers file
        'x': 0.60,                                              # threshold
        's': 3,                                                 # max-offset
        'c': True,                                              # cut-off-primers
        # read merging opts
        'ngmerge-path': os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'binaries', 'NGmerge'
        ),                                                      # ngmerge path
        'm': 20,                                                # min-overlap
        'p': 0.1,                                               # mismatch-frac
        'q': 33,                                                # phred-offset
    }

    # Validate arguments and fill in `args` dictionary.
    for opt, arg in opts:

        # General opts
        if opt == '--tasks':
            args['tasks'] = arg.split(',')
            isnot_permitted_task = lambda x: not x in TASKS
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
                args['1'] = os.path.abspath(arg)
            else:
                print('\aError: file `{}` does not exist.'.format(arg))
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end if

        elif opt in ('-2', '--R2'):
            if os.path.exists(arg):
                args['2'] = os.path.abspath(arg)
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
                args['o'] = os.path.abspath(arg)
            # end if

        elif opt in ('-t', '--threads'):
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

        elif opt in ('-z', '--gzip-output'):
            if arg == '0':
                args['z'] = False
            elif arg == '1':
                args['z'] = True
            else:
                print('Error: invalid value passed with option {}: `{}`.'\
                    .format(opt, arg))
                print('It must be 0 or 1.')
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end if

        # Crosstalks opts
        elif opt in ('-r', '--primers'):
            if os.path.exists(arg):
                args['r'] = os.path.abspath(arg)
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

        elif opt in ('-q', '--phred-offset'):
            try:
                args['q'] = int(arg)
                if not args['q'] in (33, 64):
                    raise ValueError
                # end if
            except ValueError:
                print('Invalid Phred offset: `{}`.'.format(arg))
                print('It must be 33 or 64.')
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end try

        elif opt == '--ngmerge-path':
            if not os.path.exists(arg):
                print('\aError: file `{}` does not exist.'.format(arg))
                print('Type `python3 {} -h` for help.'.format(sys.argv[0]))
                platf_depend_exit(1)
            # end if
            if not os.access(arg, os.X_OK):
                print('\aError: NGmerge file is not executable: `{}`'.format(arg))
                print('Please, make it executable. You can do it in this way:')
                print(' chmod +x {}'.format(arg))
            # end if
            args['ngmerge-path'] = os.path.abspath(arg)
        # end if
    # end for

    # Check if all mandatory options are passed
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

    # Check if "forward" file is the same as "reverse" file.
    if not args['2'] is None:
        if args['1'] == args['2']:
            print('Error: "forward" file and "reverse" file are the same file!')
            platf_depend_exit(1)
        # end if
    # end if

    return args
# end def handle_args


def report_run_params(args):
    # Funtion prints run parameters.
    #
    # :param args: argument dictionary returned by handle_args;
    # :type args: dict;

    print()
    printlog_info('  -- Run parameters --')
    printlog_info('  General:')
    printlog_info('- Tasks: {}.'.format(', '.join(args['tasks'])))
    printlog_info('- Forward reads: `{}`.'.format(args['1']))
    if not args['2'] is None:
        printlog_info('- Reverse reads: `{}`.'.format(args['2']))
    else:
        printlog_info('- Reverse reads: none.')
    # end if
    printlog_info('- Output directory: `{}`.'.format(args['o']))
    printlog_info('- Threads: {}.'.format(args['t']))
    printlog_info('- Gzip output files afterwards: {}.'.format(args['z']))
    printlog_info('  Crosstalks detection:')
    if args['r'] is None:
        printlog_info('- Primers: standard Illumina 16S rRNA V3-V4 primers.')
    else:
        printlog_info('- Primers file: `{}`.'.format(args['r']))
    # end if
    printlog_info('- Threshold: {}.'.format(args['x']))
    printlog_info('- Max offset: {}.'.format(args['s']))
    printlog_info('- Cut off primers: {}.'.format(args['c']))
    printlog_info('  Read merging:')
    printlog_info('- Minimum overlap: {}.'.format(args['m']))
    printlog_info('- Mismatch fraction: {}.'.format(args['p']))
    printlog_info('- Phred offset: {}.'.format(args['q']))
    printlog_info('- NGmerge path: {}.'.format(args['ngmerge-path']))
    print('-'*10+'\n')
# end def report_run_params


def run_task_chain(args):
    # Function for running tasks in proper order and with proper input data.

    # Deduplicate tasks.
    nr_tasks = set(args['tasks'])

    # Configure input files for first task.
    # First task is crosstalk detection.
    curr_valid_fpaths = tuple(itertools.takewhile(
                        lambda x: not x is None,
                        [args['1'], args['2']])
    )
    curr_trash_fpaths = tuple()

    # Run crosstalks detection.
    if 'rm-crosstalks' in nr_tasks:
        arguments = srcargs.CrosstalksArguments(
            curr_valid_fpaths,                # infpaths
            get_primers_seqs(args['r']),      # primers
            args['x'],                        # threshold
            args['s'],                        # max-offset
            args['c'],                        # cut-off-primers
            args['o']                         # outdir
        )
        curr_valid_fpaths, curr_trash_fpaths = src.runners.crosstalks_runner(arguments)
    # end if

    # Run NGmerge.
    if 'ngmerge' in nr_tasks:
        arguments = srcargs.NGmergeArguments(
            curr_valid_fpaths,                # infpaths
            args['ngmerge-path'],             # ngmerge
            args['t'],                        # threads
            args['m'],                        # min_overlap
            args['p'],                        # mismatch_frac
            args['q'],                        # phred offset
            args['o']                         # outdir
        )
        curr_valid_fpaths, curr_trash_fpaths = src.runners.ngmerge_runner(arguments)
    # end if

# end def run_task_chain


def main():
    args = handle_args()
    src.filesystem.make_outdir(args['o'])
    config_logging(args['o'], __version__, __last_update_date__)
    report_run_params(args)

    run_task_chain(args)

    if args['z']:
        src.compression.gzip_outfiles(args['o'])
    # end if

    printlog_info_time('Work is completed.')
# end def main


if __name__ == '__main__':
    main()
# end if
