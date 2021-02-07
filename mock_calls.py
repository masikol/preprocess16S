#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import functools

import src.fastq
import src.provide_crosstalk_pipe
import src.compression
import src.output_filenames
import src.filesystem
import src.status_bar

cut_off_primers = True
threshold = 0.52
max_offset = 2

# False:
# primer_R1 = 'GACTACHVGGGTATCTAATCC'
# primer_R2 = 'CCTACGGGNGGCWGCAG'

primer_R1 = 'CCTACGGGNGGCWGCAG'
primer_R2 = 'GACTACHVGGGTATCTAATCC'


fastq_record_R1 = src.fastq.FastqRecord(
    '@M02531:11:000000000-AM5DG:1:1101:13628:1558 1:N:0:8',
    'CCTACGGGAGGCTGCAGTGGGGAATATTGCACAATTGTTGAAACCCTTTTGCTTCCTCTCCTCTTTTTTTTTTTTTTTCTTCTTTTTTTATTTCTCTTTCTTTTTTTTCTTTTTTTACTTTACCTATTTTTTATTCCTCTTCTATCTCCTTTCCTTCATCCTCGTTCATTCGTTTTTTTCTTTCTTTTTCCGTTTTTACTTTTTTTAATTTTCTCTTTTTCTTTTTCTCATTTCTGTTGTTCAAACTTTTTTCTTTACCCTTTTACTTCATTTTATTCTTTTTTTCTTTTTTTCTTTTT',
    '+',
    '<BCCCEC7;B@FGFGG@F8<<7,,;,CF,C,6,,;,,,,,,,,;66C,,;,6,,6,,6,66,,,;,,,7++8+++++74?B?,,,<9+6,,,,5:5@,9<A?<,,+++8,,:,,8+,3,,7,339,,,,,,+,,,,7,8,,37,,5@,3,7,33,,3,,,3,2**4,,5,2,5,*****5,,,1,5:+3*0*+**2:**3+++2)2*)+++**11*0***0****00)*))).4))).).())))..).))(.66))))().)))-.)-).66:)))).-6)),,(.5)).13()46))'
)

fastq_record_R2 = src.fastq.FastqRecord(
    '@M02531:11:000000000-AM5DG:1:1101:13628:1558 2:N:0:8',
    'CCCTCCTTTTTTCTCTCTTCCTCTTTTCTCCCCCCTCTTTCTTGCCTCTTCTTCTCTTCTCCCCCTTTCCTCCTCCTTCTCCTCTTCTCTTCCTCCTCCTTTCTACTCTTTTCCCCTCTACTCTCTTCTTTCCGCTTTCCTCTCCTTCCCTCTACTCCCCCCTTTTCTCTTTCCTTTTCTCTTTTATCCCCCTCGTTTTCACTTCTCTCTTTCCTTCCCCCCTCCTCTCCCTTTCCCCCTCTTACTTCCGTCCACTCCCCCCTCCCTCCTTCTTCCCTCTCCTCTCTCCCCCTCTTCTCC',
    '+',
    '--88,8,,,,,6,;6;,,<6<;,<;C,;6666,+8+8;;;@,,,66;;,,6,<6,,5;,,,,;,8+,9,,,54,99=?A,94,49,,9,999>>;4>,,9,994,9,9,9939,9,,60,,,26,,,,,66<6**33*33041*3139)390)***)10))))000;***2****120*0***10***0)1)1))10007*0*003***13*3****)1)0)0)/)0*00-:6.*1*-/(((0**0*93(.(.(.)(((((((((((((()).))).).((..()()).(,(-((.))--'
)


primers = (primer_R1, primer_R2)
fastq_records = (fastq_record_R1, fastq_record_R2)


# primers = [primer_R1]
# fastq_records = [fastq_record_R1]


# print(crosstalk_pipe(primers, fastq_records, threshold, max_offset))

outdir = '/home/deynonih/cager/Metagenomics/test'
# fpath_R1 = '/home/deynonih/cager/Metagenomics/data/test_R1.fastq'
# fpath_R2 = '/home/deynonih/cager/Metagenomics/data/test_R2.fastq'
# fpath_R1 = '/home/deynonih/cager/Metagenomics/data/08_MG-fecal-A2_S8_L001_R1_001.fastq.bz2'
fpath_R1 = '/home/deynonih/cager/Metagenomics/data/08_MG-fecal-A2_S8_L001_R1_001.fastq.gz'
fpath_R2 = '/home/deynonih/cager/Metagenomics/data/08_MG-fecal-A2_S8_L001_R2_001.fastq.gz'
# fpath_R2 = '/home/deynonih/cager/Metagenomics/data/08_MG-fecal-A2_S8_L001_R2_001.fastq'

# infpaths = [fpath_R1]
infpaths = (fpath_R1, fpath_R2)



def crosstalk_iteration(fastq_records, primers, threshold, max_offset, valid_files, trash_files):
    not_crosstalk, fastq_records = crosstalk_pipe(primers, fastq_records, threshold, max_offset)

    if not_crosstalk:
        src.fastq.write_fastq_records(fastq_records, valid_files)
    else:
        src.fastq.write_fastq_records(fastq_records, trash_files)
    # end if
# end def



crosstalk_pipe = src.provide_crosstalk_pipe.provide_crosstalk_pipe(cut_off_primers)
valid_fpaths, trash_fpaths = src.output_filenames.get_crosstalk_outfpaths(outdir, infpaths)
valid_files = src.filesystem.open_files_for_appending(valid_fpaths)
trash_files = src.filesystem.open_files_for_appending(trash_fpaths)

loaded_crosstalk_iteration = functools.partial(
    crosstalk_iteration,
    primers=primers, threshold=threshold, max_offset=max_offset,
    valid_files=valid_files, trash_files=trash_files
)

src.status_bar.status_bar(loaded_crosstalk_iteration, infpaths)

for file_collection in (valid_files, trash_files):
    src.filesystem.close_files(file_collection)
# end for
