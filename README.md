## Preprocessing of reads from 16S rDNA


Script "preprocess_16S" preprocesses reads from 16S regions of rDNA. It works with Illumina pair-end reads.

Excactly: it detects and removes reads, that came from other sample, relying on the information,
    whether there are PCR primer sequences in these reads or not. If required primer sequence is
    found in a read, therefore, it is a read from 16S rDNA and we need it.

Moreover, it can cut these primers off and merge reads by using 'read_merging_16S' module.

Reads should be stored in files, both of which are of fastq format or both are gzipped (i.e. .fastq.gz).
Sequences of required primers are retrieved from .fa or .fa.gz file.

Attention! This script cannot be executed by python interpreter version < 3.0!


### Using Silva:

Silva SSU Ref_Nr99 release 132 is used to merge reads in 'read_merging_16S' module.

Link to Silva: https://www.arb-silva.de/

Silva license information: https://www.arb-silva.de/silva-license-information


### Installation:
    
    bash install.sh [-o Silva_db_dir]

You need this installation, if you want to merge reads.

During the installation the 16S SSU Silva database will be downloaded and reconfigured into blast-like database.

Attention! For installation you need to place files 'constant_region_V3-V4.fasta' and 'read_merging_16S.py' in the same directory, where 'install.sh' is located.

You can specify a directory for Silva database installation by using -o option, for example:

    bash install.sh -o directory_for_Silva_db

If you do not specify installation directory, Silva database will be stored in a directory nested in your currrent directory.

### Usage:

    python preprocess16S.py [--cutoff] [-p primer_file] [-1 forward_reads -2 reverse_reads] [-o output_dir]

Interactive mode:

    python preprocess16S.py

Silent mode:

    python preprocess16S.py -p primer_file.fasta -1 reads_R1_fastq -2 reads_R1_fastq -o output_dir

### Options:

    -c or --cutoff 
        script cuts primer sequences off;
    -m or --merge-reads
        all of a sudden, if this option is specified, script will merge reads together
    -q or --quality-plot
        plot a graph of number of reads as a function of average read quality
    -p or --primers
        file, in which primer sequences are stored;
    -1 or --R1
        file, in which forward reads are stored;
    -2 or --R2
        file, in which reverse reads are stored;
    -o or --outdir
        directory, in which result files will be plased.

### Read merging:

Module 'read_merging_16S' firstly tries to merge reads by overlapping region, keeping nucleotide with higher quality. Norlmal overlap is considered as a situation, when reads align ageinst one another in the followig way:

    FFFFFFFFF-------
    ------RRRRRRRRRR,

where F means a nuleotide from forward read, and R, correspondingly, -- a nucleotide from reverse read.

If reads are considered as non-overlapping, algorithm searches for a reference in Silva database by blasting forward read (because forward read usually preforms higher sequencing quality). Then algorithm alignes the reverse read against the reference that have been found above. If forward and reverse reads align with a gap, this gap is filled with N-s, so merged read will look like:

    FFFFFFFFNNNNNRRRRRRR

Not a very convenient sequence to use, but at least we will know the length of this gap.

Reads that: 

1) do not overlap properly;
2) do not contain a constant region after merging (I mean region whish is located between V3 and V4 variable regions of 16S rRNA gene);
3) do not align against reference with credible gap

are considered as chimeras and are places in corresponding files so you can deal with them further.

Reads that align against one another in the following way (it cat be a result of incorrect primer annealing):

    --FFFFFFFFFF
    RRRRRRRRR---

are considered as too short to distinguish taxa and are places in corresponding files so you can deal with them further.


It is strongly recommended to trim your reads before merging.
If you do not do it -- be ready to see many chimeras in result files, because quality of Illumina reads decreadses towards the end of read.
Therefore, low-quality reads can be considered as non-overlappting due to sequencimg errors.

Do not forget, that it is expedient to remove short reads (< 200 b. p.) while trimming.

### Suggestions for running script in interactive mode:

If you do not specify primer file or files containing reads, script will run in interactive mode.
If you run it in interactive mode, follow suggestions below:
1) It is recommended to place primer file and read files in your current directory, because the program will automatically detect them
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


Script "preprocess_16S.py":
- gzip;
- gnuplot, if you want to plot a graph (-q option);

Module "merging_reads_16S" uses:
- fasta36;
- blastn, blastdbcmd;

For correct work of the 'read_merging_16S' module you need no have fasta36 and blastn installed on your computer.
Moreover, they should be added to the path.

It is quite awkward to use both these aligners in one module, so I'll work in this direction.

Script "install.sh" uses:
- makeblastdb;

