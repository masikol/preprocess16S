# preprocess16S


This script preprocesses reads from 16S regions of rDNA. It works with Illumina pair-end reads.
Excactly: it detects and removes reads, that came from other sample, relying on the information,
    whether there are PCR primer sequences in these reads or not. If required primer sequence is
    found in a read, therefore, it is a read from 16S rDNA and we need it.

Moreover, it can cut these primers off.
Reads should be stored in files, both of which are of fastq format or both are gzipped (i.e. .fastq.gz).
Sequences of required primers are retrieved from .fa or .fa.gz file.

Attention! This script cannot be executed by python interpreter version < 3.0!

Usage:

    python preprocess16S.py [--cutoff] [-p primer_file] [-1 forward_reads -2 reverse_reads] [-o output_dir]

Options:

    -c or --cutoff 
        script cuts primer sequences off;
    -p or --primers
        file, in which primer sequences are stored;
    -1 or --R1
        file, in which forward reads are stored;
    -2 or --R2
        file, in which reverse reads are stored;
    -o or --outdir
        directory, in which result files will be plased.

If you do not specify primer file or files containing reads, script will run in interactive mode.
If you run it in interactive mode, follow suggestions below:
1) It is recommended to place primer file and read files in your directory, because the program will automatically detect them
	if they are in the current directory.
2) It is recommended to name primer file as 'primers.fasta' or something like it.
    Excactly: the program will find primer file automatically, if it's name
    contains word "primers" and has .fa, .fasta or .mfa extention.
3) It is recommended to keep names of read files in standard Illumina format.
3) If these files will be found automatically, you should confirm utilizing them
    by pressing ENTER in appropriate moments.
4) Result files named '...16S.fastq.gz' and '...trash.fastq.gz' will be
    placed in the directory nested in the current directory, if -o option is not specified.
5) This output directory will be named preprocess16S_result... and so on according to time it was ran.

This script uses:
- gzip;
