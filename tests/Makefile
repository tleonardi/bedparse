SHELL:=/bin/bash

# Extended funcionality tests for bedparse
#
# In order to run these tests the following files are needed:
# ensembl.bed (e.g. from ucsc)
# ensembl_5putr: (e.g. from ucsc)
# ensembl_3putr: (e.g. from ucsc)
# ensembl_introns: (e.g. from ucsc)
# gencode.gtf: from the gencode website
# gencode.bed: converted from the gtf with gtfToGenePred and genePredToBed
# bad_bed.bed: example of bed file with bad formatting
# good_bed: same as bad_bed.bed but with fixed formatting
# This Makefile automatically downloads the file above as 'test_data.tgz'
#

all: fivep threep introns gtf2bed operations

.PHONY: data

test_data.tgz:
	wget -O test_data.tgz "https://raw.githubusercontent.com/tleonardi/bedparse/master/tests/test_data.tgz"

data: ensembl.bed 
ensembl.bed: test_data.tgz
	tar -zxvf test_data.tgz

.PHONY: fivep
fivep: data
	bedparse 5pUTR ensembl.bed | bedparse bed12tobed6 > extracted_5p
	awk 'BEGIN{FS=OFS="\t"}{gsub("_utr.*", "", $$4); print $$1,$$2,$$3,$$4,$$5,$$6}' ensembl_5putr > ensembl_5putr_renamed
	bedparse filter -a <( cut -f4 extracted_5p) ensembl_5putr_renamed > ensembl_5putr_renamed_filtered
	diff <(sort -k1,1 -k2,2n -k3,3n ensembl_5putr_renamed_filtered) <(sort -k1,1 -k2,2n -k3,3n extracted_5p)>/dev/null

.PHONY: threep
threep: data
	bedparse 3pUTR ensembl.bed | bedparse bed12tobed6 > extracted_3p
	awk 'BEGIN{FS=OFS="\t"}{gsub("_utr.*", "", $$4); print $$1,$$2,$$3,$$4,$$5,$$6}' ensembl_3putr > ensembl_3putr_renamed
	bedparse filter -a <( cut -f4 extracted_3p) ensembl_3putr_renamed > ensembl_3putr_renamed_filtered
	diff <(sort -k1,1 -k2,2n -k3,3n ensembl_3putr_renamed_filtered) <(sort -k1,1 -k2,2n -k3,3n extracted_3p)>/dev/null

.PHONY: introns
introns: data 
	bedparse introns ensembl.bed | bedparse bed12tobed6 > extracted_introns
	awk 'BEGIN{FS=OFS="\t"}{gsub("_intron.*", "", $$4); print $$1,$$2,$$3,$$4,$$5,$$6}' ensembl_introns > ensembl_introns_renamed
	bedparse filter -a <( cut -f4 extracted_introns) ensembl_introns_renamed > ensembl_introns_renamed_filtered
	diff <(sort -k1,1 -k2,2n -k3,3n ensembl_introns_renamed_filtered) <(sort -k1,1 -k2,2n -k3,3n extracted_introns)>/dev/null

.PHONY: gtf2bed
gtf2bed: data
	diff <(bedparse gtf2bed gencode.gtf | sort -k1,1 -k2,2n -k3,3n) <(sort -k1,1 -k2,2n -k3,3n gencode.bed) >/dev/null

.PHONY: validate
validate: data
	diff <(bedparse validateFormat -f bad_bed.bed) good_bed.bed >/dev/null

.PHONY: operations
operations:
	python3 tests.py

success: fivep threep gtf2bed operations

clean:
	rm -f \
	ensembl.bed \
	ensembl_5putr \
	ensembl_3putr \
	ensembl_introns \
	gencode.gtf \
	gencode.bed \
	extracted_5p \
	ensembl_5putr_renamed \
	ensembl_5putr_renamed_filtered \
	extracted_3p \
	ensembl_3putr_renamed \
	ensembl_3putr_renamed_filtered \
	ensembl_introns_renamed \
	extracted_introns \
	ensembl_introns_renamed \
	ensembl_introns_renamed_filtered \
	good_bed.bed \
	bad_bed.bed
