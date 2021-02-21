## preprocess16S changelog

## 2021-02-21 edition.

- fixed bug that would cause preprocess16S not to print error message if `NGmerge` binary file is not executable.

### Version changes:

`5.0.c -> 5.0.d`

## 2021-02-20 evening edition.

- fixed bug that would cause preprocess16S to ignore option `-o`.

### Version changes:

`5.0.b -> 5.0.c`

## 2021-02-20 edition.

- preprocess16S now checks it "forward" and "reverse" files are actually different files.
- Small bug fix.

### Version changes:

`5.0.a -> 5.0.b`

## 2021-02-18 edition.

A major revision was done.

- Maximum offset and threshold values in crosstalk detection can now be defined by user.
- Removed reference-guided merging of non-overlapping reads until better time.
- Removed creation of quality plot (it is by no means necessary here).
- Changed names of some options.

### Version changes:

`4.0.d -> 5.0.a`

## 2020-12-16 edition.

Flag `--fill-gaps` is replaced with `--no-ovlp-merge` with the same functionality. The reason is that users used to have hard time understanding it's description in `README.md`, and inappropriate name of the flag just exacerbated the frustration.

### Version changes:

- preprocess16S: `4.0.c --> 4.0.d`
- read_mering_16S: `4.0.d --> 4.0.c`

## 2020-10-23 edition. Version 4.0.c

1. Fixed bug that would cause `preprocess16S` to write "primers were not trimmed" to log if primers actually WERE trimmed.

## 2020-07-14 edition. Version 4.0.b

1. Fixed bug that would cause scripts to terminate execution if BLAST+ is not installed, but `--fill-gaps` is not specified.

## 2020-07-13 edition. Version 4.0.a

1. First smerging step is replaced with NGmerge program.
2. Second merging step is improved: it now uses information from alignment of reverse read if alignment of forward read is ambiguous.
3. NGmerge's quality profile is adopted.

Unfortunately, changes were not tracked in this way before 2020-07-13.