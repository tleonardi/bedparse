SHELL:=/bin/bash
# In order to run these tests the following files are needed:
# ensembl.bed (e.g. from ucsc)
# ensembl_5putr: (e.g. from ucsc)
# ensembl_3putr: (e.g. from ucsc)
# gencode.gtf: from the gencode website
# gencode.bed: converted from the gtf with gtfToGenePred and genePredToBed
# bedtools binary
all: fivep threep gtf2bed operations

.PHONY: fivep
fivep:
	bedparse 5pUTR ensembl.bed | bedtools bed12tobed6 > extracted_5p
	awk 'BEGIN{FS=OFS="\t"}{gsub("_utr.*", "", $$4); print $$1,$$2,$$3,$$4,$$5,$$6}' ensembl_5putr > ensembl_5putr_renamed
	bedparse filter -a <( cut -f4 extracted_5p) ensembl_5putr_renamed > ensembl_5putr_renamed_filtered
	diff <(sort -k1,1 -k2,2n -k3,3n ensembl_5putr_renamed_filtered) <(sort -k1,1 -k2,2n -k3,3n extracted_5p)>/dev/null

.PHONY: threep
threep:
	bedparse 3pUTR ensembl.bed | bedtools bed12tobed6 > extracted_3p
	awk 'BEGIN{FS=OFS="\t"}{gsub("_utr.*", "", $$4); print $$1,$$2,$$3,$$4,$$5,$$6}' ensembl_3putr > ensembl_3putr_renamed
	bedparse filter -a <( cut -f4 extracted_3p) ensembl_3putr_renamed > ensembl_3putr_renamed_filtered
	diff <(sort -k1,1 -k2,2n -k3,3n ensembl_3putr_renamed_filtered) <(sort -k1,1 -k2,2n -k3,3n extracted_3p)>/dev/null

.PHONY: gtf2bed
gtf2bed:
	diff <(bedparse gtf2bed gencode.gtf | sort -k1,1 -k2,2n -k3,3n) <(sort -k1,1 -k2,2n -k3,3n gencode.bed) >/dev/null

.PHONY: operations
operations:
	python3 tests.py

success: fivep threep gtf2bed operations

clean:
	rm -f extracted_5p ensembl_5putr_renamed ensembl_5putr_renamed_filtered extracted_3p ensembl_3putr_renamed ensembl_3putr_renamed_filtered
