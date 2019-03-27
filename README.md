# preprocess16S

Attention! This script cannot be executed by python interpreter version < 3.0!

This script processes reads from 16S regions of rDNA. It works with Illumina pair-end reads.
Excactly: it removes those reads, that came from other loci, relying on the information,
    whether there are primer sequences in these reads or not. If required primer sequence is
    found in a read, therefore, it is a read from 16S rDNA and we need it.
Moreover, it can cut these primers off, if you specify [-c|--cutoff] option.
Reads should be stored in files, both of which are of fastq format or both be gzipped (i.e. fastq.gz).
Sequences of required primers are retrieved from .fa or fa.gz file.

## Usage:
    python preprocess16S.py [--cutoff] -p <primer_file> -1 <forward_reads> -2 <reverse_reads> [-o output_dir]

## Options:
    1) -c or --cutoff   ---   script cuts primer sequences off;
    2) -p or --primers   ---   file, in which primer sequences are stored;
    3) -1 or --R1   ---   file, in which forward reads are stored;
    4) -2 or --R2   ---   file, in which reverse reads are stored;
    5) -o or --outdir   ---   directory, in which result files will be plased.

If you do not specify some of input files or outpur directory via CL arguments, follow suggestions below:
1) It is recommended to place primer file and read files in the same directory with this .py file.
2) It is recommended to name primer file as 'primers.fasta' or something like it.
    Excactly: the program will find primer file automatically, if it's name
    contains word "primers" and has .fa, .fasta or .mfa extention.
3) It is recommended to keep names of read files in standard Illumina format.
3) If these files will be found automatically, you should confirm utilizing of these files
    by pressing ENTER is appropriate moments.
4) Result files named '...16S.fastq.gz' and '...trash.fastq.gz' will be
    placed in the directory nested in directory, where this .py file is located.
5) This output directory will be named preprocess16S_result... and so on according to time it was ran.

Uses:
- gzip;
