# FasTAN: A Fast Tandem Repeat Finder  
<font size ="4">**_Author:  Gene Myers_**<br>
**_First:   Sept 30, 2025_**<br>
**_Release v0.5: Oct. 18, 2025_**<br>

- [FasTAN](#FasTAN) 

## Overview

Once a production version is developed a proper description of **FasTAN**
will be found here.

<a name="FasTAN"></a>

## FasTAN Reference

```
FasTAN [-vM] [-T(8)] <source:path>[<fa_extn>|<1_extn>|.1gdb] <target>[.1aln]
                  
           <fa_extn> = (.fa|.fna|.fasta)[.gz]
           <1_extn>  = any valid 1-code sequence file type

      -v: Verbose mode, output statistics as proceed.
      -T: Number of threads to use.
      -M: Make a .1ano mask of the hits found.
```

FasTAN takes a FASTA or 1GDB as input and outputs the off-diagonal alignment delimiting
tandem repeats and their period as a .1aln file.
If the -M option is set then FasTAN further produces a .1ano file with the same root as
the .1aln file, that contains the intervals spanned by the alignments.
For details on a .1aln and .1ano files and how to convert the former to PAF or PSL format, see
the repo [FASTGA](https://www.github.com/thegenemyers/FASTGA).
