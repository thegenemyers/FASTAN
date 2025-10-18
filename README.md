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
FasTAN <source:path>[<precursor] <target:path>[.1aln]
                  
    <precursor> = .1gdb | <fa_extn> | <1_extn>
    
    <fa_extn> = (.fa|.fna|.fasta)[.gz]
    <1_extn>  = any valid 1-code sequence file type
```

FasTAN takes a FASTA or 1GDB as input and outputs the off-diagonal alignment delimiting
tandem repeats and their period as a .1aln file.
For details on a .1aln file and how to convert is to PAF or PSL format, see
the repo [FASTGA](www.github.com/thegenemyers/FASTGA).  For a program **tanbed** for
converting the .1aln to a BED-file with the period see
Richard Durbin's [alntools](www.github.com/richarddurbin/alntools).