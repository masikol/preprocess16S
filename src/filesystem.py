# -*- coding: utf-8 -*-

import os
import sys
import glob

from src.platform import platf_depend_exit

def open_files_for_appending(fpaths):

    outfiles = list()

    for fpath in fpaths:
        outfiles.append(open(fpath, 'a', encoding='utf-8'))
    # end for

    return outfiles
# end def open_files_for_appending


def close_files(files):
    for file in files:
        file.close()
    # end for
# end def close_files


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
                sys.exit(0)
            else:
                print('Invalid reply: `{}`'.format(reply))
            # end if
        # end while
# end def make_outdir
