from distutils.core import setup
from shutil import rmtree

setup(
name="read_merging_16S",
version="2.0",
description="Merging together Illumina PE reads from 16S rDNA",
author="Maxim Sikolenko",
author_email="maximdeynonih@gmail.com",
py_modules=["read_merging_16S"] 
)

rmtree("build") # seems this directory isn't necessary
