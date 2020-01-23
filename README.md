[![Build Status](https://travis-ci.org/tleonardi/bedparse.svg?branch=master)](https://travis-ci.org/tleonardi/bedparse)
[![Docs Status](https://readthedocs.org/projects/bedparse/badge/?version=latest&style=flat)](https://bedparse.readthedocs.io/en/latest/)
[![JOSS Status](http://joss.theoj.org/papers/22763a3b37fde13e548e884edd3221fa/status.svg)](http://joss.theoj.org/papers/22763a3b37fde13e548e884edd3221fa)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2578820.svg)](https://doi.org/10.5281/zenodo.2578820)
[![License: MIT](https://img.shields.io/badge/License-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
# Bedparse

![](docs/bedparse.svg)

Bedparse is a simple python module and CLI tool to perform common operations on BED files.

It offers 11 sub-commands that implement the following functionality:
* `filter`: Filtering of transcripts based on annotations
* `join`: Joining of annotation files based on transcript names
* `gtf2bed`: Conversion from GTF to BED format
* `convertChr`: Conversion from UCSC to Ensembl chromosome names (and viceversa)
* `bed12tobed6`: Conversion from bed12 to bed6
* `promoter`: Promoter reporting
* `introns`: Intron reporting
* `cds`: CDS reporting
* `3pUTR` and `5pUTR`: UTR reporting 
* `validateFormat`: Check that the file conforms with the BED format

## Installation

Installing is as simple as:

```
pip install bedparse
```

## Basic usage

The basic syntax in the form: `bedparse subcommand [parameters]`.

For a list of all subcommands and a brief explanation of what they do, use: `bedparse --help`.

For a detailed explanation of each subcommand and a list of its parameters, use the `--help` option after the subcommand's name, e.g.: `bedparse promoter --help`

## Documentation

Our documentation is hosted on [Read the Docs](https://bedparse.readthedocs.io/en/latest/).

We also have a short [tutorial](https://bedparse.readthedocs.io/en/latest/Tutorial.html) to guide you through the basic functions.

