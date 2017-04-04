bed-tools-gg
===

A few scripts to manage bed files, more details on the format are available [here](https://genome.ucsc.edu/FAQ/FAQformat.html#format1).

Cheers!

└[∵┌]└[ ∵ ]┘[┐∵]┘

## `bed2matrix.sh`

```
 usage: ./bed2matrix.sh [-hn] [BEDFILEs]...

 Description:
  Merge bedfiles into a matrix. The score column is merged based on the positon
  given by the chr+start+end columns (default) or by the name column (-n option)

 Mandatory arguments:
  BEDFILEs  Bed file(s). Expected to be ordered per condition.

 Optional arguments:
  -h        Show this help page.
  -n        Merge bedfiles based on name instead of location.
```

## `bed_shuffle.R`

Shuffle the reads (i.e., the score values) of the score column of a bed file. The score values are not merely shuffled but considered as read counts, then the reads are shuffled.

Saves the current random number generator seed status at the end of the script in `OUTDIR/.Random.seed.RData`. Subsequent runs of the script, with the same `OUTDIR`, will re-load the seed status unless `-k F` is used.

```
usage: bed_shuffle.R [-h][-n NITER][-p PERC][-o OUTDIR]
                     [-k KEEPSEED][-t THREADS] seed bedfile

Shuffle bed file reads.

positional arguments:
  seed              Seed for random number generation.
  bedfile           Path to bedfile.

flags:
  -h, --help            show this help message and exit

optional arguments:
  -n NITER       Number of iterations.
                 [default: 100]
  -p PERC        Percentage of reads to shuffle.
                 [default: 10]
  -o OUTDIR      Output directory.
                 [default: ./shuffled/]
  -k KEEPSEED    Reload previous seed state. Use on subsequent runs.
                 [default: TRUE]
  -t THREADS     Number of threads for parallelization.
                 [default: 1]

```

## `beds_shuffle.sh`

Runs `bed_shuffle.R` on the given bedfile(s). The initial seed is consistently kept through the R script runs as it is saved and re-loaded.

```
 usage: ./beds_shuffle.sh [-h][-n nIter][-p perc][-o outDir]
                          -s seed [BEDFILE]...

 Description:
  Shuffle a certain percentage of reads in the given bed files.

 Mandatory arguments:
  -s seed   Seed for random number generation.
  BEDFILE   Bed file(s). In any order.

 Optional arguments:
  -h    Show this help page.
  -n nIter  Number of iterations. Default: 100
  -p perc   Percentage of reads to shuffle. Default: 10
  -o outDir Output directory. Default: ./shuffled/
  -t threads    Number of threads for parallelization. Default: 1
```