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


def get_ngmerge_outprefixes(forw_fpath):
    forw_bname = os.path.basename(forw_fpath)
    more_common_name = forw_bname.partition('R1')[0].strip('_')

    merged_basename = '{}.merged.fastq'.format(more_common_name)
    unmerged_prefix = '{}.unmerged'.format(more_common_name)

    return merged_basename, unmerged_prefix
# end def get_ngmerge_outprefixes
