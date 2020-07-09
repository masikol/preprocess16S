# 16S rDNA reads preprocessing

- [Description](#description)
- [Pre-requirements](#pre-requirements)
- [Installation](#installation)
- [Usage and examples](#usage)
- [Read merging](#read-merging)
- [Silva](#silva)
- [Plotting]($plotting)

## Description

1. Script **"preprocess16S.py"** is designed for preprocessing Illumina paired-end reads of 16S rDNA amplicons.

It's main dedication is to detect and remove reads, that originate from other samples (aka **"cross-talks"**) depending on presence/absence of PCR primer sequences in these reads. If PCR primer sequences are found at 5'-ends of both PE reads, therefore, it is a read from 16S rDNA and we need it. Otherwise this pair of reads is considered as cross-talk and discarded. Primer sequences are then trimmed.

Moreover, script can:
1) merge paired-end reads together with `read_merging_16S` module (`-m` option);
2) make a plot representing distribution of average read quality (`-q` option).

Input files should be in fastq format (or `fastq.gz` or `.fastq.bz2`).

It is written in Python 3 and does npt support Python 2.

2. Module **"read_merging_16S.py"** is designed for merging Illumina paired-end reads of 16S rDNA.

It can be used as separate script (see ["Examples"](#examples) section below).

## Pre-requirements

- `BLAST+` tookit is reqiured for read merging.

`BLAST+` toolkit can be downloaded [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download);

Download archive, extract it and add `bin` directory to PATH. Then go to ["Installation"](#installation) section.

- [NGmerge](https://github.com/harvardinformatics/NGmerge) is required for read merging. It is bundled (version 0.3) with preprocess16S, and there is no need to install it separately.

- Script `preprocess16S.py` numpy and matplotlib Python packages to create a plot (`-q` option); See ["Plotting"]($plotting) section for details.

## Installation

### ~~~
### Script `preprocess16S.py` itself can be used without any installation. Installation is necessary only for read merging (`-m` option).
### ~~~

To **install `read_merging_16S`** module you need to go to downloaded `preprocess16S` directory and run:

`bash install_read_merging_16S.sh [-o Silva_db_dir]`

During the installation the **16S SSU Silva database** will be downloaded and reconfigured into blast-like database.

You can specify a directory for Silva database installation by using `-o` option, for example:

`bash install_read_merging_16S.sh -o directory_for_Silva_db`

If Silva database is already downloaded to the directory specified with `-o` option, downloading and configuring of this database will be omitted.

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
    variable regions (DOI [10.1371/journal.pone.0007401](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0007401));

-m (--min-overlap) --- minimum overlap of the paired-end reads to be merged with NGmerge.
    Default value: 20 nt.

-p (--mismatch-frac) -- fraction mismatches to allow in the overlapped region
    (a fraction of the overlap length).
    Default value: 0.1.
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
        variable regions (DOI [10.1371/journal.pone.0007401](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0007401));

    -m (--min-overlap) --- minimum overlap of the paired-end reads to be merged with NGmerge.
        Default value: 20 nt.

    -p (--mismatch-frac) -- fraction mismatches to allow in the overlapped region
        (a fraction of the overlap length).
        Default value: 0.1.
```

See ["Read merging"](#read-merging) section below for details.


## Read merging:

Module `read_merging_16S` firstly merges reads with [NGmerge](https://github.com/harvardinformatics/NGmerge), discarding [dovetailed](https://github.com/harvardinformatics/NGmerge) alignments.

For reads, discarded by NGmerge, reference-guided algorithm is applied. It can merge reads with short overlaps and even without overlap, filling the gap with N's.

Firstly, script aligns forward read against reference one with Smith-Waterman algorithm.

If reads overlap normally, script merges them. Normal overlap is alignment like this:

```
    FFFFFFFFF-------
    ------RRRRRRRRRR,
```

where F means a nucleotide from forward read, and R, correspondingly -- a nucleotide from reverse read.

If reads are considered as non-overlapping, algorithm searches for a best-hit reference in Silva database by "blasting" forward read (because forward read usually shows higher sequencing quality). Then algorithm aligns the reverse read against bets-hit reference and merges reads according to their alignments against the reference. If forward and reverse reads align with a gap, this gap is filled with Ns, so merged read will look like:

```
    FFFFFFFFNNNNNRRRRRRR
```

## Silva:

Silva SSU Ref_Nr99 (trunc) release 138 is used to merge reads in `read_merging_16S` module.

Link to Silva: https://www.arb-silva.de/

Silva license information: https://www.arb-silva.de/silva-license-information

## Plotting

You will need to install **matplotlib** and **numpy** Python packages if you want to create a plot.
These packages can be installed via `pip` or `conda` (e.g. `pip3 install numpy`).

It is a plot representing distribution of average read quality.

Plotting is final step of script's work. Whilst plotting, script ignores cross-talks and unmerged reads.
