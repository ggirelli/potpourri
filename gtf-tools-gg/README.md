gtf-tools-gg
===

A few scripts to manage gtf files, more details on the format are available [here](http://www.ensembl.org/info/website/upload/gff.html).

Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.':

1. seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
2. source - name of the program that generated this feature, or the data source (database or project name)
3. feature - feature type name, e.g. Gene, Variation, Similarity
4. start - Start position of the feature, with sequence numbering starting at 1.
5. end - End position of the feature, with sequence numbering starting at 1.
6. score - A floating point value.
7. strand - defined as + (forward) or - (reverse).
8. frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
9. attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

Cheers!

└[∵┌]└[ ∵ ]┘[┐∵]┘

## `./gtf_convert_ids.sh`

```
 usage: ./gtf_convert_id.sh [-hf] -g gtf -c conv

 Description:
  Convert IDs in the specified field of a GTF attribute column. Requires a GTF
  file and a conversion file obtained from bioDBnet.
  bioDBnet: bioDBnet: https://biodbnet-abcc.ncifcrf.gov/db/dbOrg.php

 Mandatory arguments:
  -g gtf    GTF file.
  -c conv   Conversion file, obtained from bioDBnet.

 Optional arguments:
  -h        Show this help page.
  -v        Skip rows without the specified attribute field.
  -r        Skip rows with no conversion ID.
  -f field  Label of the attribute column field containing the ID to convert.
            Default: ref_gene_id
```

## `./gtf_copy_attr.sh`

```
 usage: ./gtf_copy_attr.sh [-h] -g gtf -1 from -2 to

 Description:
  Overwrites the value of attribute (-2) with the value of attribute (-1)

 Mandatory arguments:
  -g gtf    GTF file.
  -1 from   Attribute value to copy from.
  -2 to     Attribute value to copy to.

 Optional arguments:
  -h        Show this help page.
```

