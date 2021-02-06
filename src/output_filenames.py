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

        extention = _get_extention(infpath)

        valid_fpaths.append(
            os.path.join(
                outdir,
                infpath.replace(extention, '.16S.fastq')
            )
        )

        trash_fpaths.append(
            os.path.join(
                outdir,
                infpath.replace(extention, '.trash.fastq')
            )
        )
    # end for

    return valid_fpaths, trash_fpaths
# end def get_crosstalk_outfpaths