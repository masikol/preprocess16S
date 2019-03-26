# preprocess16S

Attention! This script cannot be executed by python interpreter version < 3.6!
If your interpreter throws an unhandled SyntaxError -- make sure you are using python 3.6+.

This script processes reads from 16S regions of rDNA.
Excactly: it removes those reads, that came from other loci, relying on the information,
    whether there are primer sequences in these reads or not. If required primer sequence is
    found in a read, therefore, it is a read from 16S rDNA and we need it.
Moreover, it cuts these primers off.
For now you can utilize this script only by interaction with it (by pressing ENTER in appropriate moments).

Reads should be stored in files, both of which are of fastq format or both be gzipped (i.e. fastq.gz).
Sequences of required primers are retrieved from .fa or fa.gz file.
It is recommended to place file with primers and files with reads
    in the same directory with this .py file.
It is recommended to name file with primers as 'primers.fasta' or something like it.
Excactly: the program will find file with primers automatically, if it's name
    contains word "primers" and has .fa, .fasta or .mfa extention.
Result files named '...16S.fastq.gz' and '...trash.fastq.gz' will be
    placed in the directory nested in directory, where this .py file is located.
This result directory will be named 'preprocess16S_result...' and so on according to time script was ran.

TODO:
- get rid of unhandled SyntaxError in case of utilizing old python interpreter versions;
- add ability to specify all required files and options via CL arguments;
