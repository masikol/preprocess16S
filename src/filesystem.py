# -*- coding: utf-8 -*-

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