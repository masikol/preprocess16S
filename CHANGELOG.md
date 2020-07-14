## preprocess16S changelog

## 2020-07-14 edition. Version 4.0.b

1. Fixed bug that would cause scripts to terminate execution if BLAST+ is not installed, but `--fill-gaps` is not specified.

## 2020-07-13 edition. Version 4.0.a

1. First smerging step is replaced with NGmerge program.
2. Second merging step is improved: it now uses information from alignment of reverse read if alignment of forward read is ambiguous.
3. NGmerge's quality profile is adopted.