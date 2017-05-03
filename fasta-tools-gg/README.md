fasta-tools-gg
===

A few scripts to manage fasta files, more details on the format are available [here](http://zhanglab.ccmb.med.umich.edu/FASTA/).

Briefly:
> a sequence in FASTA format begins with a single-line description, followed by lines of sequence data. The description line is distinguished from the sequence data by a greater-than (">") symbol in the first column. It is recommended that all lines of text be shorter than 80 characters in length. An example sequence in

Cheers!

└[∵┌]└[ ∵ ]┘[┐∵]┘

## `./fasta_filter.py`

```
usage: fasta_filter.py [-h] FASTA regexp

Select sequences from the input FASTA file that match a given regular
expression in their header.

positional arguments:
  FASTA       Path to a fasta file.
  regexp      Regular expression to match to the header.

optional arguments:
  -h, --help  show this help message and exit
```
