# -*- coding: utf-8 -*-

import os
import re


def _get_extention(fpath):
    re_obj = re.search(r'(\.f(ast)?q(\.[bg]z(2)?)?)', fpath)
    return re_obj.group(1)
# end def _get_extention


def get_crosstalk_outfpaths(outdir, infpaths):
    valid_fpaths = list()
    trash_fpaths = list()

    for infpath in infpaths:

        infpath_bname = os.path.basename(infpath)
        extention = _get_extention(infpath_bname)

        valid_fpaths.append(
            os.path.join(
                outdir,
                infpath_bname.replace(extention, '.16S.fastq')
            )
        )

        trash_fpaths.append(
            os.path.join(
                outdir,
                infpath_bname.replace(extention, '.trash.fastq')
            )
        )
    # end for

    return valid_fpaths, trash_fpaths
# end def get_crosstalk_outfpaths