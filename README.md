[![Build Status](https://travis-ci.org/tleonardi/bedparse.svg?branch=master)](https://travis-ci.org/tleonardi/bedparse)

# Bedparse
Bedparse is a simple python module and CLI tool to perform common operations on BED files.

*This program is under active development and is likely to contain bugs.*

## Installation

```
git clone git@github.com:tleonardi/bedparse.git
cd bedparse
```
and then:
```
pip3 install .
```

## Usage

```
usage: bedparse [-h] [--version]
                {3pUTR,5pUTR,cds,promoter,introns,filter,join,gtf2bed,bed12tobed6,convertChr}
                ...

Perform various simple operations on BED files.

positional arguments:
  {3pUTR,5pUTR,cds,promoter,introns,filter,join,gtf2bed,bed12tobed6,convertChr}
                        sub-command help
    3pUTR               Prints the 3' of coding genes.
    5pUTR               Prints the 5' of coding genes.
    cds                 Prints the CDS of coding genes.
    promoter            Prints the promoters of coding genes.
    introns             Prints BED records corresponding to the introns of
                        each transcript in the original file.
    filter              Filters a BED file based on an annotation.
    join                Joins a BED file with an annotation file using the BED
                        name (col4) as the joining key.
    gtf2bed             Converts a GTF file to BED12 format.
    bed12tobed6         Converts a BED12 file to BED6 format.
    convertChr          Convert chromosome names between UCSC and Ensembl
                        formats

optional arguments:
  -h, --help            show this help message and exit
  --version, -v         show program's version number and exit

```
