# 16S rDNA reads preprocessing

- [Description](#description)
- [Pre-requirements](#pre-requirements)
- [Installation](#installation)
- [Usage and examples](#usage)
- [Read merging](#read-merging)
- [Silva](#silva)
- [Plotting]($plotting)

## Description

1. Script **"preprocess16S.py"** is designed to preprocess reads from 16S regions of rDNA. It works with Illumina pair-end reads.

Version `3.3.0`; 2020.01.10 edition;

It's main dedication is to detect and remove reads, that came from other samples (aka **"cross-talks"**), relying on the information,
whether there are PCR primer sequences in these reads or not. More precisely: if required primer sequence is
found at least in one of two corresponding PE reads, therefore, it is a read from 16S rDNA and we need it.
Otherwise this pair of reads is considered as cross-talk.

Moreover, script can:
1) trim these primer sequences (`-c` option);
2) merge reads by using `read_merging_16S` module (`-m` option);
3) plot a graph representing distribution of average read quality (`-q` option).

Reads should be stored in files, both of which are of `fastq` format or both are gzipped (i.e. `.fastq.gz`).
Sequences of required primers are retrieved from multi-fasta (`.mfa` or `.fasta` or `.fa`) file.
Gzipped primer files are allowed.


2. Module **"read_merging_16S.py"** is designed to merge Illumina (MiSeq) pair-end reads from 16S rDNA.

Version `3.3.0`; 2020.01.10 edition;

It can be used as script as well as be imported like Python module from outer Python program and then used via calling `merge_reads()` function.

**Attention!** These scripts cannot be executed by Python interpreter version **< 3.0**!

## Pre-requirements

- `BLAST+` tookit is reqiured for read merging (`-m` option of "preprocess16S.py").

`BLAST+`toolkit can be downloaded [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download);

Download tar archive, extract it and add `bin` directory to PATH. Then go to ["Installation"](#installation) section.

- Script `preprocess16S.py` numpy and matplotlib Python packages to plot a graph (`-q` option); See ["Plotting"]($plotting) section for details.

## Installation

### ~~~
### Script `preprocess16S.py` itself can be used without any installation. Installation is necessary only for read merging (`-m` option).
### ~~~

To **install `read_merging_16S`** module you need to go to `preprocess16S/read_merging_16S` directory and run:

`bash install_read_merging_16S.sh [-o Silva_db_dir]`

During the installation the **16S SSU Silva database** will be downloaded and reconfigured into blast-like database.

You can specify a directory for Silva database installation by using `-o` option, for example:

`bash install_read_merging_16S.sh -o directory_for_Silva_db`

If Silva database is already downloaded to the directory specified with `-o` option,
downloading and configuring of this database will be omitted.

If you do not specify installation directory, Silva database will be stored in a directory nested in your current directory.

## Usage:

**Attention!**

You might have problems with paths completion while calling Python interpreter explicitly.

Therefore I recommend to make files executable (`chmod +x preprocess16S.py`) and run them in the following way:

`./preprocess16S.py -p primers.mfa -1 forward_R1_reads.fastq.gz -2 reverse_R2_reads.fastq.gz -o outdir/`

### 1. preprocess16S.py:

#### Options:

```
-h (--help) --- show help message.
    '-h' -- brief, '--help' -- full;

-v (--version) --- show version;

-c (--trim-primers) --- Flag option. If specified, primer sequences will be trimmed;

-m (--merge-reads) --- Flag option. If specified, reads will be merged together
    with 'read_merging_16S' module;

-q (--quality-plot) --- Flag option. If specified, a graph of read quality distribution
    will be plotted. Requires 'numpy' and 'matplotlib' Python packages;

-p (--primers) --- FASTA file, in which primer sequences are stored;
    Illumina V3-V4 primer sequences are used by default:
    https://support.illumina.com/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf
    See "Amplicon primers" section.

-1 (--R1) --- FASTQ file, in which forward reads are stored;

-2 (--R2) --- FASTQ file, in which reverse reads are stored;

-o (--outdir) --- directory, in which result files will be placed.

-t (--threads) <int> --- number of threads to launch;
    Default value is 1.
    Attention: whole files meant to be processed will be loaded to memory if number of theads if grated than 1.
    So, if you are limited to memory, be careful while using this option.
    Anyway, an option that performs more memory-efficient behaviour wil be disigned soon.
```

#### Examples:

1) Filter cross-talks from files 'forw_R1_reads.fastq.gz' and 'rev_R2_reads.fastq.gz' depending on
     default Illumina V3-V4 primer sequenes. Trim primer sequences and put result files in
     the directory named 'outdir':

`./preprocess16S.py -1 forw_R1_reads.fastq.gz -2 rev_R2_reads.fastq.gz -o outdir -c`

2) Filter cross-talks from files 'forw_R1_reads.fastq.gz' and 'rev_R2_reads.fastq.gz' depending on
     primer sequenes in file 'V3V4_primers.fasta'. Trim primer sequences, merge reads,
     plot a quality graph and put result files in the directory named 'outdir'. Launch 4 threads to calculations:

```
./preprocess16S.py -p V3V4_primers.fasta -1 forw_R1_reads.fastq.gz -2 rev_R2_reads.fastq.gz \
-o outdir -cmq -t 4
```

### 2. read_merging_16S.py:

`./read_merging_16S.py -1 forward_reads -2 reverse_reads [-o output_dir]`

For example:

`./read_merging_16S.py -1 forward_R1_reads.fastq.gz -2 reverse_R2_reads.fastq.gz -o outdir/`

#### Options:
```
    -h (--help) --- show help message;

    -v (--version) --- show version;

    -1 (--R1) --- FASTQ file, in which forward reads are stored;

    -2 (--R2) --- FASTQ file, in which reverse reads are stored;

    -o (--outdir) --- directory, in which result files will be placed;

    -t (--threads) <int> --- number of threads to launch;
        Default value is 1.
        Attention: whole files meant to be processed will be loaded to memory if number of theads if grated than 1.
        So, if you are limited to memory, be careful while using this option.
        Anyway, an option that performs more memory-efficient behaviour will be disigned soon.
```

See ["Read merging"](#read-merging) section below for more presice information about read merging.


#### Using `read_merging_16S` as Python module

If you want to use read_merging_16S as imported Python module, follow steps below:

1. Configure paths to both read files and (optionally) path to output directory and store these paths in `str` objects.

2. Call `read_merging_16S.merge_reads()` function, passing to it paths mentioned above as follows:

```
result_files = read_merging_16S.merge_reads(forward_R1_reads, reverse_R2_reads, outdir)
```

This function returns a `dict<str: str>` of the following format:

```
    {   
        "merg": path to a file with successfully merged reads,
        "umR1": path to a file with forward unmerged reads,
        "umR2": path to a file with reverse unmerged reads,
        "shrtR1": path to a file with forward reads considered as too short,
        "shrtR2": path to a file with reverse reads considered as too short
    }
```

3. After read merging you can get some statistics via calling `get_merging_stats()` function, as follows:

`merging_stats = read_merging_16S.get_merging_stats()`

Function `get_merging_stats()` does not take any argument and returns a `dict<int: int>` of the following format:

```
    {
        0: merged_read_pairs_number,
        1: unmerged_read_pairs_number,
        2: too_short_reads_number
    }
```

Function rises an `AccessStatsBeforeMergingError` on the attempt of 
accessing merging statistics before merging (e.i. before calling `merge_reads()` function).


## Read merging:

Module `read_merging_16S` firstly tries to merge reads by finding overlapping region using naive algorithm:
    it takes the tail of forward read and slides it through the reverse read until significant similarity.
If no significant similarity is found, function tries to merge reads by overlapping region (using Smith-Waterman algorithm),
    keeping nucleotide with higher quality. 
Normal overlap is considered as a situation, when reads align against one another in the following way:

```
    FFFFFFFFF-------
    ------RRRRRRRRRR,
```

where F means a nucleotide from forward read, and R, correspondingly -- a nucleotide from reverse read.

If reads are considered as non-overlapping, algorithm searches for a reference in Silva 
database by blasting forward read (because forward read usually performs higher sequencing quality).
Then algorithm aligns the reverse read against the reference that have been found above. 
If forward and reverse reads align with a gap, this gap is filled with N-s, so merged read will look like:

```
    FFFFFFFFNNNNNRRRRRRR
```

Not a very convenient sequence to use, but at least we will know the length of this gap.

Reads are considered as unmergable according to following criterions:

1) reads do not overlap properly;
2) reads do not align against reference with credible gap

Reads that align against one another in the following way (it can be a result of incorrect primer annealing):

```
    --FFFFFFFFFF
    RRRRRRRRR---
```

are considered as too short to distinguish taxa and are placed
in corresponding files so you can deal with them further.


It is strongly recommended to **trim** your reads before merging.
If you do not do it -- be ready to wait a while, because reads with high error rate
can't be properly merged without blasting them against Silva database,
and this procedure in turn requires long time to be done.

## Silva:

Silva SSU Ref_Nr99 (trunc) release 138 is used to merge reads in `read_merging_16S` module.

Link to Silva: https://www.arb-silva.de/

Silva license information: https://www.arb-silva.de/silva-license-information

## Plotting

You will need to install **matplotlib** and **numpy** if you want to use this feature.
These packages can be installed via `pip` or `conda` (e.g. `pip3 install numpy`).

As it has been mentioned above, if `-q` (`--quality-plot`) option is specified, a quality graph will be plotted.

It is a graph representing distribution of reads by their average quality.

Distribution is shown for positive result reads only (e.i. for those from 16S rDNA or for successfully merged reads if you merge them).
