Motivation
===============

The BED (Browser Extensible Data) format is a plain text file format commonly used in bioinformatics to represent genomic features (e.g. genes, transcripts, peaks, regulatory regions, etc.). Each line in the file represents a genomic feature and consists of up to 12 tab-separated fields:

1. chromosome name
2. start coordinate in the chromosome
3. end coordinate in the chromosome
4. feature name
5. feature score
6. strand
7. thick start (conventionally the start codon for protein coding transcripts)
8. thick end (conventionally the stop codon for protein coding transcripts)
9. rgb color for visualisation in genome browsers
10. number of connected blocks (conventionally the number of exons)
11. comma separated list of blocks size
12. comma separated list of block starts relative to field 2 (i.e. genomic start of the feature)

One of the major advantages of the BED format over many of its alternatives is that each line includes all the information required to define an individual transcript. This characteristic allows to perform numerous operations on BED a file as part of unix pipes, for example using GNU awk.

For example, the following is a common approach to extract gene promoters (here defined as 500bp around the gene start)::

   awk 'BEGIN{OFS=FS="\t"}{print $1,$2-500,$3+500,$4,$5}' transcritpome.bed > promoters.bed

However, these one-liners can quickly get long and hard to read. For example, if we wanted to do the same as before but keeping the strand into considerations::

   awk 'BEGIN{OFS=FS="\t"}{if($6=="+"){print $1,$2-500,$2+500,$4,$5}else{print $1,$3-500,$3+500,$4,$5}}' transcritpome.bed > promoters_stranded.bed

These and other more complex operations quicly get long to type and prone to errors and typos. Bedparse greatly simplifies the process::

   bedparse promoter transcritpome.bed > promoters_stranded.bed

or::
   
   bedparse promoter --unstranded transcritpome.bed > promoters.bed

Despite the simplicity of most of its operations, all functions in bedparse are thouroughly and rigourously tested through an automated test suit to ensure the accuracy and correctness of the results. Additionally, bedparse performs syntax validation checks on the input BED files and warns the user in case of malformed or unsupported formats.

Additionally, bedparse also provides two format conversion operations:

* gtf2bed allows converting Ensembl/Gencode Gene transfer format (GTF) files into bed format
* convertChr implements an internal dictionary that allows conversion of human and mouse chromosome names between the two most widely used formats, i.e. the Ensembl and the UCSC naming schemes.



