# Preprocessing 16S rDNA Illumina reads with preprocess16S

- [Description](#description)
- [Pre-requirements](#pre-requirements)
- [Installation](#installation)
- [Usage and examples](#usage)
- [Read merging](#read-merging)
- [Silva](#silva)
- [Plotting]($plotting)

## Description

1. Script **"preprocess16S.py"** is designed for preprocessing Illumina paired-end reads of 16S rDNA amplicons.

It's main purpose is to detect and remove reads, that originate from other samples (aka **"cross-talks"**) depending on presence/absence of PCR primer sequences in these reads. If PCR primer sequences are found at 5'-ends of both PE reads, therefore, it is a read from 16S rDNA and we need it. Otherwise this pair of reads is considered as cross-talk and discarded. Primer sequences are then trimmed.

Moreover, script can:
1) merge paired-end reads together with `read_merging_16S` module (`-m` option);
2) make a plot representing distribution of average read quality (`-q` option).

Input files should be in fastq format (or `fastq.gz` or `.fastq.bz2`).

Script is written in Python 3 and does not support Python 2.

2. Module **"read_merging_16S.py"** is designed for merging Illumina paired-end reads of 16S rDNA.

It can be used as separate script (see ["Examples"](#examples) section below).

## Pre-requirements

- [NGmerge](https://github.com/harvardinformatics/NGmerge) is required for read merging. It is bundled (version 0.3) with preprocess16S, and there is no need to install it separately.

- Script `preprocess16S.py` uses numpy and matplotlib Python packages to create a plot (`-q` option); See ["Plotting"]($plotting) section for details.

- `BLAST+` tookit is reqiured for [gap-filling](#read-merging) read merging.

    `BLAST+` toolkit can be downloaded [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download);

    Download archive, extract it and add `bin` directory to PATH. Then go to ["Installation"](#installation) section.


## Installation

Script `preprocess16S.py` itself can be used in-place without any installation.

Installation is necessary only for [gap-filling](#read-merging) read merging. It needs Silva SSU database to be downloaded and configured as BLAST database.

To configure Silva database for gap-filling merging, you need to go to downloaded `preprocess16S` directory and run:

`bash configure_Silva_db.sh`

You can specify installation directory for Silva database with `-o` option, for example:

`bash configure_Silva_db.sh -o directory_for_Silva_db`

If configured Silva database is already in directory specified with `-o` option, downloading and configuring of this database will be omitted.

If you do not specify installation directory, Silva database will be stored in a directory nested in your current directory.

## Usage:

You might have problems with paths completion while calling Python interpreter explicitly.

Therefore I recommend to make files executable (`chmod +x preprocess16S.py`) and run them in the following way:

`./preprocess16S.py -p primers.mfa -1 forward_R1_reads.fastq.gz -2 reverse_R2_reads.fastq.gz -o outdir/`

### 1. preprocess16S.py:

#### Options:

```
-h (--help) --- show help message.
    '-h' -- brief, '--help' -- full;

-v (--version) --- show version;

-k (--keep-primers) --- Flag option. If specified, primer sequences will not be trimmed;

-q (--quality-plot) --- Flag option. If specified, a graph of read quality distribution
    will be plotted. Requires 'numpy' and 'matplotlib' Python packages;

-p (--primers) --- multi-FASTA file, in which primer sequences are stored (one line per sequence);
    Illumina V3-V4 primer sequences are used by default.

-1 (--R1) --- FASTQ file, in which forward reads are stored;

-2 (--R2) --- FASTQ file, in which reverse reads are stored;

-o (--outdir) --- directory, in which result files will be placed.

-t (--threads) <int> --- number of threads to launch;
    Default value is 1.

-f (--phred-offset) [33, 64] --- Phred quality offset (default -- 33);

Read merging options

-m (--merge-reads) --- Flag option. If specified, reads will be merged together;

--ngmerge-path -- path to NGmerge executable.
    You can specify it if bundled NGmerge 0.3 is not siutable for you.

-N (--num-N) --- maximum length of a gap that preprocess16S can fill with Ns.
    Default value: 35 -- length of conservative region between V3 and V4
    variable regions (DOI 10.1371/journal.pone.0007401);

-m (--min-overlap) --- minimum overlap of the paired-end reads to be merged with NGmerge.
    Default value: 20 nt.

-p (--mismatch-frac) -- fraction mismatches to allow in the overlapped region
    (a fraction of the overlap length).
    Default value: 0.1.

--fill-gaps -- Flag option. If specified, the second, gap-filling step
    of read merging will be applied after NGmerge.
    Disabled by default.
```

#### Examples:

1) Filter cross-talks from files 'forw_R1_reads.fastq.gz' and 'rev_R2_reads.fastq.gz' depending on
     default Illumina V3-V4 primer sequenes:

`./preprocess16S.py -1 forw_R1_reads.fastq.gz -2 rev_R2_reads.fastq.gz -o outdir`

2) Filter cross-talks from files 'forw_R1_reads.fastq.gz' and 'rev_R2_reads.fastq.gz' depending on
     primer sequenes in file 'V3V4_primers.fasta'. Keep primer sequences, merge reads,
     create a quality plot and put result files in the directory named 'outdir'. Use 4 threads:

```
./preprocess16S.py -p V3V4_primers.fasta -1 forw_R1_reads.fastq.gz -2 rev_R2_reads.fastq.gz \
-o outdir -kmq -t 4
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

-f (--phred-offset) [33, 64] --- Phred quality offset (default -- 33);

--ngmerge-path -- path to NGmerge executable.
    You can specify it if bundled NGmerge 0.3 is not siutable for you.

-N (--num-N) --- maximum length of a gap that preprocess16S can fill with Ns.
    Default value: 35 -- length of conservative region between V3 and V4
    variable regions (DOI 10.1371/journal.pone.0007401);

-m (--min-overlap) --- minimum overlap of the paired-end reads to be merged with NGmerge.
    Default value: 20 nt.

-p (--mismatch-frac) -- fraction mismatches to allow in the overlapped region
    (a fraction of the overlap length).
    Default value: 0.1.

--fill-gaps -- Flag option. If specified, the second, gap-filling step
    of read merging will be applied after NGmerge.
    Disabled by default.
```

See ["Read merging"](#read-merging) section below for details.


## Read merging

#### 1. First merging step.

Module `read_merging_16S` firstly merges reads with [NGmerge](https://github.com/harvardinformatics/NGmerge), discarding [dovetailed](https://github.com/harvardinformatics/NGmerge) alignments.

#### 2. Second merging step (optional).

This step can be enabled with `--fill-gaps` flag.

For reads, discarded by NGmerge, reference-guided algorithm is applied. It can merge reads with short overlaps and even without overlap, filling the gap with N's.

This step is beneficial if reads do not overlap or length of overlap is small. Also this step can merge reads, whose 3'-ends were previously trimmed by QC software due to low quality. However, this step is very (!) slow, since it aligns each forward read against Silva SSU Ref_Nr99 database with BLAST.

At this step script uses Silva SSU database to determine location of forward and reverse reads in 16S rRNA gene and merge them according to this location.

If forward and reverse reads align against reference sequence with a gap, this gap is filled with Ns, so merged read will look like:

```
    FFFFFFFFNNNNNRRRRRRR
```

## Silva:

Silva SSU Ref_Nr99 (trunc) release 138 is used to merge reads in `read_merging_16S` module.

Silva: [https://www.arb-silva.de/](https://www.arb-silva.de/)

Silva license information: [https://www.arb-silva.de/silva-license-information](https://www.arb-silva.de/silva-license-information)

## Plotting

You will need to install **matplotlib** and **numpy** Python packages if you want to create a plot.
These packages can be installed via `pip` or `conda` (e.g. `pip3 install numpy`).

Plotting is final step of script's work. Whilst plotting, script ignores cross-talks and unmerged reads.
