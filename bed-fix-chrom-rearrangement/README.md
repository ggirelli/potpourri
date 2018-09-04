bed-fix-chrom-rearrangement
===

```
usage: bedFix_chrRear.py [-h] -1 suffix -2 suffix [-S suffix] [--version]
                         bedfile [bedfile ...]

Corrects chromosomal rearrangements in bed files. Corrected chromosome sites
format checked with assertions. Only first three columns of the bed files are
expected and manipulated (i.e., BED3). If preasent, header track line is
preserved and translocation site information are added as "transSite1" and
"transSite2". Corrected bed files are exported with suffix ".transCorrected".

positional arguments:
  bedfile               Bed file(s).

optional arguments:
  -h, --help            show this help message and exit
  -1 suffix, --first-site suffix
                        First translocation site, format chrA:NNNNNN.
  -2 suffix, --second-site suffix
                        Second translocation site, format chrB:NNNNNN.
  -S suffix, --suffix suffix
                        Suffix for output bed files. Default:
                        '.transCorrected'
  --version             show program's version number and exit
```
