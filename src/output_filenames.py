# -*- coding: utf-8 -*-
# Module contains function, which decide how output files should be named.

import os
import re


def _get_extention(fpath):
    # Function detects and returns extention of fastq file.
    # If file name is, e.g. `file.fastq.gz`, extention is `.fastq.gz`, not just `.gz`.
    #
    # :param fpath: path to fastq file;
    # :type fpath: str;
    re_obj = re.search(r'(\.f(ast)?q(\.[bg]z(2)?)?)', fpath)
    return re_obj.group(1)
# end def _get_extention


def get_crosstalk_outfpaths(outdir, infpaths):
    # Function confugures paths to output files for crosstalks detection task.
    #
    # :param outdir: path to outdir;
    # :type outdir: str;
    # :param infpaths: collection of paths to input files;
    # :type infpaths: list<str>;
    #
    # Returns two collections:
    # 1. A collection of valid (non-crosstalk) paths.
    # 2. A collection of trash (crosstalk) paths.

    valid_fpaths = list()
    trash_fpaths = list()

    for infpath in infpaths:

        # Get basename
        infpath_bname = os.path.basename(infpath)
        # Get extention
        extention = _get_extention(infpath_bname)

        # Configure outfpaths
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
    # Function confugures paths and prefixes to output files for NGmerge task.
    #
    # :param forw_fpath: path to forward file;
    # :type forw_fpath: str;
    #
    # Returns two strings:
    # 1. Basename of "merged" file.
    # 2. Prefix for "unmerged" files.

    # Get basename
    forw_bname = os.path.basename(forw_fpath)
    # Remove _R1 stuff
    more_common_name = forw_bname.partition('R1')[0].strip('_')

    # Make name for "merged" file
    merged_basename = '{}.merged.fastq'.format(more_common_name)
    # Make prefix for "unmerged" files
    unmerged_prefix = '{}.unmerged'.format(more_common_name)

    return merged_basename, unmerged_prefix
# end def get_ngmerge_outprefixes
