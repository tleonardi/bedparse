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
                {3pUTR,5pUTR,cds,promoter,introns,filter,join,gtf2bed,bed12tobed6}

Perform various simple operations on BED files.

positional arguments:
  {3pUTR,5pUTR,cds,promoter,introns,filter,join,gtf2bed,bed12tobed6}
                        sub-command help
    3pUTR               Prints the 3' of coding genes.
    5pUTR               Prints the 5' of coding genes.
    cds                 Prints the CDS of coding genes.
    promoter            Prints the promoters of coding genes.
    introns             Prints BED records corresponding to the introns of
                        each transcript in the original file.
    filter              Filters a BED file based on an annotation. BED entries
                        with a name (i.e. col4) that appears in the specified
                        column of the annotation are printed to stdout. For
                        efficiency reasons this command doesn't perform BED
                        validation.
    join                Joins a BED file with an annotation file using the BED
                        name (col4) as the joining key.
    gtf2bed             Converts a GTF file to BED12 format. This tool
                        supports the Ensembl GTF format. The GTF file must
                        contain 'transcript' and 'exon' features in field 3.
                        If the GTF file also annotates 'CDS' 'start_codon' or
                        'stop_codon' these are used to annotate the thickStart
                        and thickEnd in the BED file.
    bed12tobed6         Converts a BED12 file to BED6 format.

optional arguments:
  -h, --help            show this help message and exit
  --version, -v         show program's version number and exit
```
