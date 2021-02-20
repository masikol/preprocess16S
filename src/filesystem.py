# -*- coding: utf-8 -*-
# Module contains functions for basic interation with file system objects.

import os
import sys
import glob

from src.platform import platf_depend_exit


def open_files_for_appending(fpaths):
    # Function opens files for appending and returns them.
    #
    # :param fpaths: collection of paths to files to open;
    # :type fpaths: list<str>;
    #
    # Returns collection of files opened for apending (list<_io.TextIOWrapper>).

    outfiles = list()

    for fpath in fpaths:
        outfiles.append(open(fpath, 'a', encoding='utf-8'))
    # end for

    return outfiles
# end def open_files_for_appending


def close_files(files):
    # Function just closes files.
    #
    # :param files: collection of file objects;
    # :type files: list<_io.TextIOWrapper>, list<gzip.GzipFile>, list<bz2.BZ2File>;
    for file in files:
        file.close()
    # end for
# end def close_files


def make_outdir(outdir):
    # Funciton makes output directory.
    # Function warns a user if outdir is not empty.
    #
    # :aram outdir: path to outdir;
    # :type outdir: str;

    if not os.path.exists(outdir):
        # Create outdir if it doesn't exist.
        try:
            os.makedirs(outdir)
        except OSError as err:
            print('Cannot create output directory: {}'.format(err))
            platf_depend_exit(1)
        # end try

    elif len(os.listdir(outdir)) != 0:
        # It outdir is not empty -- warn user and ask if he/she wants to empty it now.

        print('\nOutput directory `{}` is not empty.'.format(outdir))

        error = True
        while error:
            reply = input("""Press ENTER to remove all files in it and proceed
 or enter `q` to exit\n >> """)

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
            elif reply.lower() == 'q':
                # Just exit
                sys.exit(0)
            else:
                print('Invalid reply: `{}`'.format(reply))
            # end if
        # end while
# end def make_outdir
