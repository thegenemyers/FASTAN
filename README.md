# FasTAN: A Fast Tandem Repeat Finder  

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://bioconda.github.io/recipes/fastan/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/fastan/badges/version.svg)](https://anaconda.org/bioconda/fastan)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/fastan/badges/downloads.svg)](https://anaconda.org/bioconda/fastan)

<font size ="4">**_Author:  Gene Myers_**<br>
**_First:   Sept 30, 2025_**<br>
**_Release v0.8: April 25, 2026_**<br>

- [FasTAN](#FasTAN) 

## Overview

**FasTAN** is a tandem array finder under development.
The capabilities of such a utlity can be partitioned into 3 parts:
(a) detection of intervals of the subject sequence that contain a tandem array with
an estimate of the unit size,
(b) partitioning of a detected interval into individual subunits, and
(c) reporting a model of the units composing a specific array.
The current version of FasTAN can perform (a) and (b) but not yet (c).

FasTAN is part of the "Fast" series of tools developed primarily by EWM.
As such its outputs are ONE files describing alignments (files with extension .1aln)
or annotations (files with extension .1ano) and a number of useful utilities for
manipulating these are in the [FASTGA](https://www.github.com/thegenemyers/FASTGA)
repo which we strongly suggest you also install
for this reason.  In addition there is a QT-based viewer
[ALNVIEW](https://www.github.com/thegenemyers/ALNVIEW) that allows you to
visualize and peruse the alignments in a .1aln file.

A .1ano file is our version of a BED-file and there are utilites in the FASTGA repo
(ANOtoBED and BEDtoANO) that allow you to convert between the two.  However, .1ano
files can also contain additional information not contemplated by the BED format,
such as the partitioning of a TA-interval into units and the display of a model of
a unit.  For those of you that are not C-programmers and/or do not want to directly
work with ONE files, you can convert a binary ONE file into an easy to parse ASCII
version with ONEview (which is also in the FASTGA repo).

<a name="FasTAN"></a>

## FasTAN Reference

```
Usage: FasTAN [-vamp] [-T(8)] [-o<target>] <source:path>[<fa_extn>|<1_extn>|.1gdb]

           <fa_extn> = (.fa|.fna|.fasta)[.gz]
           <1_extn>  = any valid 1-code sequence file type

      -v: Verbose mode, output statistics as proceed.
      -a: Make a .1aln of all the hits found.
      -m: Make a .1ano mask of the hits found.
      -p: Parse hits into units.

      -o: Root path of .1aln/.1ano file (default root path of input).
      -T: Number of threads to use.
```

FasTAN takes a FASTA, a ONE sequence, or a 1GDB file as input.
With the -a option set it produces a .1aln file of the "first-wave" alignments
of the tandem arrays it finds.  The first wave alignment is that sequence alignment
that aligns consequtive units to each other and as such its level of similarity is
an indicator of the level of similarity/conservation between the units of the array.
With the -m option set it produces a .1ano file of the intervals containining
the detected tandem arrays along with an estimate of the average unit size and
the average identity of the alignment.
If the -p option is set, then the -m option must also be set, and the .1ano file
also contains a P-line that describes the partitioning of the TA interval into
units.

By default 8 threads are used and the output file names are the root name of the
input with a .1aln or .1ano suffix as per the options.  You can direct the output
to a particular root path with the -o option, and control the number of threads
with the -T option.