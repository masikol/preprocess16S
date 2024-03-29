# -*- coding: utf-8 -*-
# Module contains function for obtaining primer sequence(s).

import re

from src.printlog import printlog_error
from src.platform import platf_depend_exit


def get_primers_seqs(primers_fpath):
    # Function for for obtaining primer sequence(s).
    # If primer_fpath is None, it returns refault primers:
    #   Illumina 16S V3-V4 primers.
    # Otherwise it parses primers from provided fasta file.

    # Use Illumina V3-V4 primers by dafault
    if primers_fpath is None:
        primers = ('CCTACGGGNGGCWGCAG', 'GACTACHVGGGTATCTAATCC')
    else:
        primers = list()

        # Get lines
        try:
            with open(primers_fpath, 'r') as primers_file:
                lines = primers_file.readlines()
            # end with
        except OSError as oserror:
            printlog_error('Error while reading file of primers: {}'\
                .format(oserror))
            platf_depend_exit(1)
        # end try

        # Remove blank lines
        lines = list(filter(lambda x: x != '', lines))

        # There must be 1 or 2 primers in primers file.
        if len(lines) not in (2, 4):
            printlog_error('Error: invalid format of primers file.\
It should be single (2 lines at all) or "double" (4 lines at all) fasta file.\
Bu there are {} lines in your file.'.format(len(lines)))
            platf_depend_exit(1)
        # end if

        bases = 'AGCTUTYSWKMBDHVN'

        # Validate sequence(s).
        for i in range(1, len(lines), 2):
            seq = lines[i].strip().upper()
            if re.match(r'[{}]+'.format(bases), seq) is None:
                printlog_error('Error: invalid character in primer sequence.\
Here is invalid primer sequence: `{}`. Permitted characters: `{}`'\
                    .format(seq, bases))
                platf_depend_exit(1)
            # end if
            primers.append(seq)
        # end for
    # end if

    return primers
# end def get_primer_seqs
