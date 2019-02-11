# Bedparse tutorial


![](bedparse.svg)



Hi, thanks for your interest in `bedparse`!

The following is a short tutorial that will guide you through the functionality of `bedparse`. You can find the `example.bed` file in this repo under docs/example.bed.
This file contains 6 human transcript models from Gencode. The first three are non-coding transcripts (i.e. field 7 and 8 contain the same coordinate), whereas the last three are coding (i.e. fields 7 and 8 indicate the thickStart and thickEnd, i.e. start and end of the [CDS](https://en.wikipedia.org/wiki/Coding_region)).

## Extracting the promoters

The `bedparse promoter` command reports the promoter of each transcript, defined as a user specified interval around the [TSS](https://en.wikipedia.org/wiki/Transcription_start_site).
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

Note how the TSS (and as a consequence the promoter) depends on the strand: for transcripts on the negative strand the TSS is the end coordinate, i.e. column 3. The `--unstranded` option allows you to override this behaviour and report promoters as an interval around column 2, thus disregarding the strand.

## Extracting 5' or 3' UTRs

The [UTRs](https://en.wikipedia.org/wiki/Untranslated_region) are defined in a BED file as the region between the start (column 2) and the thickStart (column 7) for the 5' and between the thickEnd (column 8) and the end (column 3) for the 3'. These rules are reversed for transcripts on the - strand, and `bedparse` automatically takes care of this. Additionally, `bedparse` also handles correctly UTRs that span multiple exons: in these cases `bedparse` recomputes all exon starts and exon lengths as sets coulmns 11 and 12 accordingly.

```
$ bedparse 5pUTR example.bed 
chr1    923927  924431  ENST00000420190.6       0       +       923927  923927  0       1       504,    0,
chr1    925149  925941  ENST00000437963.5       0       +       925149  925149  0       2       40,20,  0,772,
chr1    925737  925941  ENST00000342066.7       0       +       925737  925737  0       2       63,20,  0,184,
```

Clearly, as you can see from the output above, UTRs are only reported for coding transcripts.

## Extracting the CDS

To extract the CDS of the coding transcripts in the BED file use the `bedparse cds` command:

```
$ bedparse cds example.bed 
chr1    924431  939291  ENST00000420190.6       0       +       924431  939291  0       7       517,92,182,51,125,90,17,        0,1490,5723,6607,11340,14608,14843,
chr1    925941  935793  ENST00000437963.5       0       +       925941  935793  0       4       72,182,51,22,   0,4213,5097,9830,
chr1    925941  944153  ENST00000342066.7       0       +       925941  944153  0       13      72,182,51,125,90,186,163,116,79,500,125,111,246,        0,4213,5097,9830,13098,13333,15202,16194,16468,16617,17311,17756,17966,

```

Note how non-coding transcripts are not reported (because by definition they don't have a CDS). Also, note how the number of exons (column 10) and exon lengths and starts (columns 11 and 12) have been readjusted to reflect the fact that the transcripts have "lost" the UTRs. To visualise this operation you can save the output of the command above to a new text file and upload it as a custom track in the UCSC genome browser: you'll see that the new transcripts only correspond to the thick portion of the original Gencode transcripts.


## Extracting introns

In a BED file introns are implicitly defined as the genomic regions between exons. The `bedparse introns` command creates new "artificial" transcripts that correspond to the introns of the original transcripts:

```
$ bedparse introns example.bed 
chr1    12227   13220   ENST00000456328.2       0       +       12227   12227   0       2       385,499,        0,494,
chr1    12057   13452   ENST00000450305.2       0       +       12057   12057   0       5       121,385,277,168,78,     0,170,640,995,1317,
chr1    14501   29533   ENST00000488147.1       0       -       14501   14501   0       10      503,757,659,92,177,237,172,206,6371,4642,       0,537,1446,2264,2554,2867,3241,3560,3865,10390,
chr1    924948  939274  ENST00000420190.6       0       +       924948  924948  0       6       973,4141,702,4682,3143,145,     0,1065,5388,6141,10948,14181,
chr1    925189  935771  ENST00000437963.5       0       +       925189  925189  0       4       732,4141,702,4682,      0,824,5147,5900,
chr1    925800  943907  ENST00000342066.7       0       +       925800  925800  0       13      121,4141,702,4682,3143,145,1683,829,158,70,194,320,99,  0,213,4536,5289,10096,13329,13660,15506,16451,16688,17258,17577,18008,
```

## Convert BED12 to BED6

It's often convenient to convert a BED12 file into BED6, where each exon appears on its own line. This is easily done with `bedparse bed12tobed6`:

```
$ bedparse bed12tobed6 --appendExN example.bed 
chr1    11868   12227   ENST00000456328.2_Exon001       0       +
chr1    12612   12721   ENST00000456328.2_Exon002       0       +
chr1    13220   14409   ENST00000456328.2_Exon003       0       +
chr1    12009   12057   ENST00000450305.2_Exon001       0       +
chr1    12178   12227   ENST00000450305.2_Exon002       0       +
chr1    12612   12697   ENST00000450305.2_Exon003       0       +
chr1    12974   13052   ENST00000450305.2_Exon004       0       +
chr1    13220   13374   ENST00000450305.2_Exon005       0       +
chr1    13452   13670   ENST00000450305.2_Exon006       0       +
chr1    14403   14501   ENST00000488147.1_Exon001       0       -
chr1    15004   15038   ENST00000488147.1_Exon002       0       -
chr1    15795   15947   ENST00000488147.1_Exon003       0       -
chr1    16606   16765   ENST00000488147.1_Exon004       0       -
chr1    16857   17055   ENST00000488147.1_Exon005       0       -
chr1    17232   17368   ENST00000488147.1_Exon006       0       -
chr1    17605   17742   ENST00000488147.1_Exon007       0       -
chr1    17914   18061   ENST00000488147.1_Exon008       0       -
chr1    18267   18366   ENST00000488147.1_Exon009       0       -
chr1    24737   24891   ENST00000488147.1_Exon010       0       -
chr1    29533   29570   ENST00000488147.1_Exon011       0       -
chr1    923927  924948  ENST00000420190.6_Exon001       0       +
chr1    925921  926013  ENST00000420190.6_Exon002       0       +
chr1    930154  930336  ENST00000420190.6_Exon003       0       +
chr1    931038  931089  ENST00000420190.6_Exon004       0       +
chr1    935771  935896  ENST00000420190.6_Exon005       0       +
chr1    939039  939129  ENST00000420190.6_Exon006       0       +
chr1    939274  939291  ENST00000420190.6_Exon007       0       +
chr1    925149  925189  ENST00000437963.5_Exon001       0       +
chr1    925921  926013  ENST00000437963.5_Exon002       0       +
chr1    930154  930336  ENST00000437963.5_Exon003       0       +
chr1    931038  931089  ENST00000437963.5_Exon004       0       +
chr1    935771  935793  ENST00000437963.5_Exon005       0       +
chr1    925737  925800  ENST00000342066.7_Exon001       0       +
chr1    925921  926013  ENST00000342066.7_Exon002       0       +
chr1    930154  930336  ENST00000342066.7_Exon003       0       +
chr1    931038  931089  ENST00000342066.7_Exon004       0       +
chr1    935771  935896  ENST00000342066.7_Exon005       0       +
chr1    939039  939129  ENST00000342066.7_Exon006       0       +
chr1    939274  939460  ENST00000342066.7_Exon007       0       +
chr1    941143  941306  ENST00000342066.7_Exon008       0       +
chr1    942135  942251  ENST00000342066.7_Exon009       0       +
chr1    942409  942488  ENST00000342066.7_Exon010       0       +
chr1    942558  943058  ENST00000342066.7_Exon011       0       +
chr1    943252  943377  ENST00000342066.7_Exon012       0       +
chr1    943697  943808  ENST00000342066.7_Exon013       0       +
chr1    943907  944575  ENST00000342066.7_Exon014       0       +
```

The optional flag --appendExN adds ExonNNN to the end of each transcript name.
