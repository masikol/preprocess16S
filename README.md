# Preprocessing 16S rDNA Illumina reads with preprocess16S

- [Description](#description)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Options](#options)
- [Examples](#examples)

## Description

Script preprocess16S.py is designed for preprocessing Illumina paired-end reads of 16S rDNA amplicons.

It's main purpose is to detect and remove reads, which originate from other samples (aka "**crosstalks**") depending on presence/absence of PCR primer sequences in these reads. If PCR primer sequences are found at 5'-ends of both PE reads, therefore, it is a read from 16S rDNA and we need it. Otherwise preprocess16S considers this pair of reads crosstalk and discards it. After that preprocess16S removes primer sequences from reads.

Limitation: preprocess16S can only detect crosstalks originating from non-16S samples.

Moreover, script can merge paired-end reads together with [NGmerge](https://github.com/jsh58/NGmerge) program;

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

Crosstalks detection:

  -r (--primers) -- FASTA file, where primers sequences are stored
      (one line per sequence).
      Illumina V3-V4 primer sequences are used by default.

  -x (--threshold) -- threshold value used in crosstalks detection;
      Real number from 0 to 1. Default: 0.52.

  -s (--max-offset) -- maximum offset used in crosstalks detection;
      Integer > 0. Default: 2.

  -с (--cut-off-primers) [0, 1]. 0 -- keep primers, 1 -- cut them off.
      Default -- 1.

Read merging:

  -m (--min-overlap) -- minimum overlap of the paired-end reads to be merged with NGmerge.
      Default: 20 nt.

  -p (--mismatch-frac) -- fraction of mismatches to allow in the overlapped region
      (a fraction of the overlap length).
      Default: 0.1.

* -t (--threads) <int> -- number of threads to launch.
      Default: 1.

  -q (--phred-offset) [33, 64] -- Phred quality offset.
      Default: 33.

  --ngmerge-path -- path to NGmerge executable.
      You can specify it if bundled NGmerge 0.3 is not suitable for you.
```

### Note

`*` Removing cross-talks in parallel makes no profit, so preprocess16S removes cross-talks in single thread anyway. Only read merging with NGmerge can be executed in parallel.


## Examples

1) Remove non-16S crosstalks from files `forw_R1_reads.fastq.gz` and `rev_R2_reads.fastq.gz` with default Illumina V3-V4 primer sequenes:

```
./preprocess16S.py -1 forw_R1_reads.fastq.gz -2 rev_R2_reads.fastq.gz -o outdir
```

2) Remove cross-talks from files `forw_R1_reads.fastq.gz` and `rev_R2_reads.fastq.gz` with primer sequenes from your own file `my_V3V4_primers.fasta`. Then merge reads with NGmerge. Use 4 threads:

```
./preprocess16S.py --tasks rm-crosstalks,ngmerge \
-1 forw_R1_reads.fastq.gz -2 rev_R2_reads.fastq.gz \
-r my_V3V4_primers.fasta \
-o outdir -t 4
```
