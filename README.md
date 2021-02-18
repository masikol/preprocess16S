# Preprocessing 16S rDNA Illumina reads with preprocess16S

- [Description](#description)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Options](#options)
- [Examples](#examples)
- [Algorithm details](#algorithm-details)

## Description

Script preprocess16S is designed for preprocessing Illumina paired-end reads of 16S rDNA amplicons.

It's main purpose is to detect and remove reads, which originate from other samples (aka "**crosstalks**") depending on presence/absence of PCR primer sequences in these reads. If PCR primer sequences are found at 5'-ends of both PE reads, therefore, it is a read from 16S rDNA and we need it. Otherwise preprocess16S considers this pair of reads as crosstalks and discards it. After that preprocess16S removes primer sequences from reads.

Also, it can do the same thing with single-end reads. In this case the program seeks for primer only at one end (5' one) of a read, therefore it is not suitable for merged PE reads.

Limitation: preprocess16S can only detect crosstalks originating from non-16S samples.

Moreover, script can merge paired-end reads together with [NGmerge](https://github.com/jsh58/NGmerge) program after removing crosstalks.

Script is written in Python 3 and does not support Python 2.

## Dependencies

- Python 3 ([https://www.python.org/](https://www.python.org/)).

- [NGmerge](https://github.com/harvardinformatics/NGmerge) is required for read merging. It is bundled (version 0.3) with preprocess16S, and there is no need to install it separately.

## Usage:

Basic usage is:

```
./preprocess16S.py --tasks task1,task2 -1 forward_R1_reads.fastq.gz -2 reverse_R2_reads.fastq.gz -o outdir/
```

## Options:

```
Print-and-exit options:

  -h (--help) -- show help message.

  -v (--version) -- show version.

General:

  --tasks -- comma-separated list of tasks to run.
      Permitted values: rm-crosstalks, ngmerge.
      Default: rm-crosstalks.

  -1 (--R1) -- FASTQ file of forward reads.

  -2 (--R2) -- FASTQ file of reverse reads.

  -o (--outdir) -- output directory.

  -z (--gzip-output) -- [0, 1]. 1 -- gzip output files after work is done.
      0 -- keep output files uncompressed.
      Default: 1.

Crosstalks detection:

* -r (--primers) -- FASTA file, where primers sequences are stored
      (one line per sequence).
      Illumina V3-V4 primer sequences are used by default.

  -x (--threshold) -- threshold value used in crosstalks detection;
      Real number from 0 to 1. Default: 0.52.
      See "Algorithm details" section below for details.

  -s (--max-offset) -- maximum offset used in crosstalks detection;
      Integer > 0. Default: 2.
      See "Algorithm details" section below for details.

  -с (--cut-off-primers) [0, 1]. 0 -- keep primers, 1 -- cut them off.
      Default -- 1.

Read merging:

  -m (--min-overlap) -- minimum overlap of the paired-end reads to be merged with NGmerge.
      Default: 20 nt.

  -p (--mismatch-frac) -- fraction of mismatches to allow in the overlapped region
      (a fraction of the overlap length).
      Default: 0.1.

* -t (--threads) -- number of threads to launch.
      Default: 1.

   -q (--phred-offset) [33, 64] -- Phred quality offset.
      Default: 33.

  --ngmerge-path -- path to NGmerge executable.
      You can specify it if bundled NGmerge 0.3 is not suitable for you.
```

### Note

`*` Illumina V3-V4 primer sequences used by preprocess16S by default can be found [here](https://support.illumina.com/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf).

`*` Removing cross-talks in parallel makes no profit, so preprocess16S removes cross-talks in single thread anyway. Only read merging with NGmerge can be executed in parallel.


## Examples

1) Remove non-16S crosstalks from files `forw_R1_reads.fastq.gz` and `rev_R2_reads.fastq.gz` with default Illumina V3-V4 primer sequenes:

```
./preprocess16S.py -1 forw_R1_reads.fastq.gz -2 rev_R2_reads.fastq.gz -o outdir
```

2) Do the same thing only with forward reads:

```
./preprocess16S.py -1 forw_R1_reads.fastq.gz -o outdir
```

3) Remove cross-talks from files `forw_R1_reads.fastq.gz` and `rev_R2_reads.fastq.gz` with primer sequenes from your own file `my_V3V4_primers.fasta`. Then merge reads with NGmerge. Use 4 threads:

```
./preprocess16S.py --tasks rm-crosstalks,ngmerge \
-1 forw_R1_reads.fastq.gz -2 rev_R2_reads.fastq.gz \
-r my_V3V4_primers.fasta \
-o outdir -t 4
```


## Algorithm details

The idea of the crosstalk detection algorithm is to find sequence of PCR primer at the start of the read.

Given a read:

`CCTACGGGAGCCTGCAGTGGGGAATATTGCACAATTGTTGAAACCCTTTTGCTTCCTCT...`

Given a primer:

`CCTACGGGNGGCWGCAG`

Now we need to introduce some values. With examples below, they should be clear for you, so go ahead ;)

1. Length of strings compared (denoted as `L` below).

2. Score (denoted as `score` below): number of matching characters in compared strings divided by `L`. In other words, it is identity ratio.

3. Offset (denoted as `offset` below): an offset, with which primer is "applied" to read sequence. preprocess16S starts with zero offset, and then increases it's absolute value in following way: `1,-1,2,-2...` untill maximum value passed with `-s` option.

4. Score threshold. If `score` is above this threshold, the algorithm decides that primer is detected. This is the threshold from `-x` option.

**Here are some illustrations of algorithm work with different parameters.**

In all examples, let threshold be 0.7.

Example 1:

```
CCTACGGGNGGCWGCAG
|||||||||| ||||||
CCTACGGGAGCCTGCAGTGGGGAATATTGCACAATTGTTGAAACCCTTTTGCTTCCTCT...

offset = 0
L = 17 (is equal to length of primer sequence)
score = 16 / L = 16 / 17 = 0.94

score > 0.7, threrefore primer is detected.
```

Example 2:
```
-CCTACGGGNGGCWGCAG
 |    || |        
CCTACGGGAGCCTGCAGTGGGGAATATTGCACAATTGTTGAAACCCTTTTGCTTCCTCT...

offset = 1
L = 17 (is equal to length of primer sequence)
score = 4 / L = 4 / 17 = 0.24

score < 0.7, threrefore there is no primer at this position.
```

Example 3:
```
CCTACGGGNGGCWGCAG
 |    ||| ||     
-CCTACGGGAGCCTGCAGTGGGGAATATTGCACAATTGTTGAAACCCTTTTGCTTCCTCT...

offset = -1
L = 16 (is equal to length of primer sequence minus 1)
score = 6 / L = 6 / 16 = 0.38

score < 0.7, threrefore there is no primer at this position.
```


