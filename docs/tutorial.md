# Bedparse tutorial

Hi, thanks for your interest in `bedparse1 `!

The following is a short tutorial that will guide you through the functionality of `bedparse`. You can find the `example.bed` file in this repo under docs/example.bed.
This files contains 6 human transcript models from Gencode. The first three are non-coding transcripts (i.e. field 7 and 8 are the same coordinate), wheread the last three are coding (i.e. fields 7 and 8 indicate the start and end of the [CDS](https://en.wikipedia.org/wiki/Coding_region)).

## Extracting the promoters

The `bedparse promoter` command reports the promoter of each transcript, defined as user specified interval around the [TSS](https://en.wikipedia.org/wiki/Transcription_start_site).
For example, we can extract promoters consisting of 1000bp upstream and 500bp downstream of the CDS:

```
$ bedparse promoter --up 1000 --down 500 example.bed 
chr1    10868   12368   ENST00000456328.2
chr1    11009   12509   ENST00000450305.2
chr1    29070   30570   ENST00000488147.1
chr1    922927  924427  ENST00000420190.6
chr1    924149  925649  ENST00000437963.5
chr1    924737  926237  ENST00000342066.7
```

Note how the TSS (and as a consequence the promoter) depends on the strand: for transcripts on the negative strand the TSS is the end coordinate, i.e. column 3. The `--unstranded` allows you to override this behaviour and report promoters as an interval around column 2, thus disregarding the strand.


## Extracting the CDS

To extract the CDS of the coding transcripts in the BED file use the `bedparse cds` command:

```
bedparse cds example.bed 
chr1    924431  939291  ENST00000420190.6       0       +       924431  939291  0       7       517,92,182,51,125,90,17,        0,1490,5723,6607,11340,14608,14843,
chr1    925941  935793  ENST00000437963.5       0       +       925941  935793  0       4       72,182,51,22,   0,4213,5097,9830,
chr1    925941  944153  ENST00000342066.7       0       +       925941  944153  0       13      72,182,51,125,90,186,163,116,79,500,125,111,246,        0,4213,5097,9830,13098,13333,15202,16194,16468,16617,17311,17756,17966,

```

Note how non-coding transcripts are not reported (because by definition they don't have a CDS). Also, note how the number of exons (column 10) and exon lengths and starts (columns 11 and 12) have been readjusted to reflect the fact that the transcripts have "lost" the UTRs. To visualise this operation you can save the output of the command above to a new text file and upload it as a custom track in the ucsc genome broser: you'll see that the new transcripts only correspond to the thick portion of the original Gencode transcripts.


## Extracting introns

In a BED file introns are implicitly defined as the genomic regions inbetween exons. The `bedparse introns` command created new "artificial" transcripts that correspond to the introns of the original transcripts:

```
$ bedparse introns example.bed 
chr1    12227   13220   ENST00000456328.2       0       +       12227   12227   0       2       385,499,        0,494,
chr1    12057   13452   ENST00000450305.2       0       +       12057   12057   0       5       121,385,277,168,78,     0,170,640,995,1317,
chr1    14501   29533   ENST00000488147.1       0       -       14501   14501   0       10      503,757,659,92,177,237,172,206,6371,4642,       0,537,1446,2264,2554,2867,3241,3560,3865,10390,
chr1    924948  939274  ENST00000420190.6       0       +       924948  924948  0       6       973,4141,702,4682,3143,145,     0,1065,5388,6141,10948,14181,
chr1    925189  935771  ENST00000437963.5       0       +       925189  925189  0       4       732,4141,702,4682,      0,824,5147,5900,
chr1    925800  943907  ENST00000342066.7       0       +       925800  925800  0       13      121,4141,702,4682,3143,145,1683,829,158,70,194,320,99,  0,213,4536,5289,10096,13329,13660,15506,16451,16688,17258,17577,18008,
```



