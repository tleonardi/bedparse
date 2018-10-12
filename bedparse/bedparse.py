#!/usr/bin/python3
import signal
import argparse
import sys
import csv
import re
from pkg_resources import get_distribution
from bedparse import bedline
from bedparse import gtf2bed
from bedparse import BEDexception
# This allows using the program in a pipe
# The program is killed when it receives a sigpipe
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
__version__ = get_distribution('bedparse').version

def introns(args):
    with args.bedfile as tsvfile:
        for line in tsvfile:
            introns=bedline(line.split('\t')).introns()
            if(introns): introns.print()
    tsvfile.close()

def threeP(args):
    with args.bedfile as tsvfile:
        for line in tsvfile:
            utr=bedline(line.split('\t')).utr(which=3)
            if(utr): utr.print()
    tsvfile.close()

def fiveP(args):
    with args.bedfile as tsvfile:
        for line in tsvfile:
            utr=bedline(line.split('\t')).utr(which=5)
            if(utr): utr.print()
    tsvfile.close()

def cds(args):
    with args.bedfile as tsvfile:
        for line in tsvfile:
            utr=bedline(line.split('\t')).cds(ignoreCDSonly=args.ignoreCDSonly)
            if(utr): utr.print()
    tsvfile.close()

def prom(args):
    with args.bedfile as tsvfile:
        for line in tsvfile:
            bedline(line.split('\t')).promoter(up=args.up, down=args.down, strand=(not args.unstranded)).print()
    tsvfile.close()

def bed12tobed6(args):
    if args.whichExon is not "all" and args.keepIntrons:
        raise BEDexception("--keepIntrons is only allowed with --whichExon all")
    with args.bedfile as tsvfile:
        for line in tsvfile:
            tx = bedline(line.split('\t'))
            exon_list = tx.bed12tobed6(appendExN=args.appendExN, whichExon=args.whichExon)
            for el in exon_list:
                el.print()
            if(args.keepIntrons):
                nameSub=re.compile("_Exon([0-9]+)")
                for el in tx.introns().bed12tobed6(appendExN=args.appendExN):
                    el.name=nameSub.sub(r"_Intron\1", el.name)
                    el.print()
    tsvfile.close()

def filter(args):
    col=args.column-1
    filterset=set()
    try:
        annotation=open(args.annotation)
    except:
        raise BEDexception("Annotation file not valid")
    annotationReader = csv.reader(annotation, delimiter="\t")
    for line in annotationReader:
        filterset.add(line[col])
    annotation.close()
    with args.bedfile as tsvfile:
        for line in tsvfile:
            if(line.split('\t')[3] in filterset):
                print(line.rstrip())
    tsvfile.close()

def join(args):
    col=args.column-1
    annot=dict()
    try:
        annotation=open(args.annotation)
    except:
        raise BEDexception("Annotation file not valid")
    annotationReader = csv.reader(annotation, delimiter=args.separator)
    for line in annotationReader:
        if(len(line)<=col):
            raise BEDexception("Some lines don't contain the annotation column")
        annot.setdefault(line[col], []).append(line[0:col]+line[col+1:])
    annotation.close()
    with args.bedfile as tsvfile:
        for line in tsvfile:
            line=line.split('\t')
            if(args.noUnmatched==False or line[3] in annot.keys()):
                record=bedline(line)
                if(record):
                        nrec=len(annot.setdefault(record.name, []))
                        if(nrec==0):
                            if(args.empty==''):
                                record.print()
                            else:
                                record.print(end='')
                                print('',args.empty,sep="\t")
                        else:
                            for i in range(0,nrec):
                                record.print(end='')
                                print('',*annot[record.name][i], sep='\t')
    tsvfile.close()

def convertChr(args):
    with args.bedfile as tsvfile:
        for line in tsvfile:
            translatedLine=bedline(line.split('\t')).translateChr(assembly=args.assembly, target=args.target, suppress=args.suppressMissing, all=args.allowMissing, patches=args.patches)
            if(translatedLine):
               translatedLine.print()
    tsvfile.close()

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
            description="""Perform various simple operations on BED files.""")
    
    parser.add_argument('--version', '-v', action='version', version='v'+__version__)
    subparsers = parser.add_subparsers(help='sub-command help', dest='sub-command')
    subparsers.required = True
    
    parser_3pUTR = subparsers.add_parser('3pUTR', help="Prints the 3' of coding genes.")
    parser_3pUTR.add_argument("bedfile", type=argparse.FileType('r'), nargs='?', default=sys.stdin, help="Path to the BED file.")
    parser_3pUTR.set_defaults(func=threeP)
    
    parser_5pUTR = subparsers.add_parser('5pUTR', help="Prints the 5' of coding genes.")
    parser_5pUTR.add_argument("bedfile", type=argparse.FileType('r'), nargs='?', default=sys.stdin, help="Path to the BED file.")
    parser_5pUTR.set_defaults(func=fiveP)
    
    parser_cds = subparsers.add_parser('cds', help="Prints the CDS of coding genes.")
    parser_cds.add_argument("--ignoreCDSonly",action="store_true", help="Ignore transcripts that only consist of CDS.")
    parser_cds.add_argument("bedfile", type=argparse.FileType('r'), nargs='?', default=sys.stdin, help="Path to the BED file.")
    parser_cds.set_defaults(func=cds)
    
    parser_prom = subparsers.add_parser('promoter', help="Prints the promoters of coding genes.")
    parser_prom.add_argument("--up",type=int, default=500, help="Get this many nt upstream of each feature.")
    parser_prom.add_argument("--down",type=int, default=500, help="Get this many nt downstream of each feature.")
    parser_prom.add_argument("--unstranded",action="store_true", help="Do not consider strands.")
    parser_prom.add_argument("bedfile", type=argparse.FileType('r'), nargs='?', default=sys.stdin, help="Path to the BED file.")
    parser_prom.set_defaults(func=prom)
    
    parser_introns = subparsers.add_parser('introns', help="Prints BED records corresponding to the introns of each transcript in the original file.")
    parser_introns.add_argument("bedfile", type=argparse.FileType('r'), nargs='?', default=sys.stdin, help="Path to the BED file.")
    parser_introns.set_defaults(func=introns)
    
    parser_filter = subparsers.add_parser('filter', 
            help="Filters a BED file based on an annotation.",
            description="""Filters a BED file based on an annotation.
            BED entries with a name (i.e. col4) that appears
            in the specified column of the annotation are
            printed to stdout. For efficiency reasons this
            command doesn't perform BED validation.""")
    parser_filter.add_argument("--annotation", "-a", type=str, help="Path to the annotation file.", required=True)
    parser_filter.add_argument("--column","-c",type=int, default=1, help="Column of the annotation file (1-based, default=1).")
    parser_filter.set_defaults(func=filter)
    parser_filter.add_argument("bedfile", type=argparse.FileType('r'), nargs='?', default=sys.stdin,
    help="Path to the BED file.")
    
    parser_join = subparsers.add_parser('join', 
            help="""Joins a BED file with an annotation file using
            the BED name (col4) as the joining key.""")
    parser_join.add_argument("--annotation", "-a", type=str, help="Path to the annotation file.", required=True)
    parser_join.add_argument("--column","-c",type=int, default=1, help="Column of the annotation file (1-based, default=1).")
    parser_join.add_argument("--separator","-s",type=str, default='\t', help="Field separator for the annotation file (default tab)")
    parser_join.add_argument("--empty","-e",type=str, default='.', help="String to append to empty records (default '.').")
    parser_join.add_argument("--noUnmatched", "-n" ,action="store_true", help="Do not print unmatched lines.")
    parser_join.set_defaults(func=join)
    parser_join.add_argument("bedfile", type=argparse.FileType('r'), nargs='?', default=sys.stdin,
    help="Path to the BED file.")
 
 
    parser_gtf2bed = subparsers.add_parser('gtf2bed', 
            help="Converts a GTF file to BED12 format.",
            description="""Converts a GTF file to BED12 format.
            This tool supports the Ensembl GTF format.
            The GTF file must contain 'transcript' and 'exon' 
            features in field 3. If the GTF file also annotates
            'CDS' 'start_codon' or 'stop_codon' these are used
            to annotate the thickStart and thickEnd in the BED
            file.""")
    parser_gtf2bed.add_argument("gtf", type=argparse.FileType('r'), nargs='?', default=sys.stdin, help="Path to the GTF file.")
    parser_gtf2bed.add_argument("--extraFields",type=str, default='', help="Comma separated list of extra GTF fields to be added after col 12 (e.g. gene_id,gene_name).")
    parser_gtf2bed.add_argument("--filterType",type=str, default='', help="Comma separated list of 'transcript_type' types to retain.")
    parser_gtf2bed.set_defaults(func=lambda args: gtf2bed(args.gtf, args.extraFields.split(','), args.filterType.split(',')))
 
    parser_bed12tobed6 = subparsers.add_parser('bed12tobed6', 
            help="""Converts a BED12 file to BED6 format.
            """)
    parser_bed12tobed6.add_argument("bedfile", type=argparse.FileType('r'), nargs='?', default=sys.stdin, help="Path to the GTF file.")
    parser_bed12tobed6.add_argument("--appendExN", action="store_true", help="Appends the exon number to the transcript name.")
    parser_bed12tobed6.add_argument("--whichExon",type=str, default='all', choices=["all", "first", "last"], help="Which exon to return. First and last respectively report the first or last exon relative to the TSS (i.e. taking strand into account).")
    parser_bed12tobed6.add_argument("--keepIntrons", action="store_true", help="Add records for introns as well. Only allowed if --whichExon all")
    parser_bed12tobed6.set_defaults(func=bed12tobed6)
    
    parser_convertChr = subparsers.add_parser('convertChr', help="Convert chromosome names between UCSC and Ensembl formats",
                                                   description="""Convert chromosome names between UCSC and Ensembl formats.
                                                                  The conversion supports the hg38 assembly up to patch 11 and
                                                                  the mm10 assembly up to patch 4. By default patches
                                                                  are not converted (because the UCSC genome browser does not
                                                                  support them), but can be enabled using the -p flag.
                                                                  When the BED file contains a chromsome that is not recognised,
                                                                  by default the program stops and throws an error. Alternatively,
                                                                  urnecognised chrosomes can be suppressed (-s) or artificially set
                                                                  to 'NA' (-a).""")
    parser_convertChr.add_argument("bedfile", type=argparse.FileType('r'), nargs='?', default=sys.stdin, help="Path to the BED file.")
    parser_convertChr.add_argument("--assembly", type=str, help="Assembly of the BED file (either hg38 or mm10).", required=True)
    parser_convertChr.add_argument("--target", type=str, help="Desidered chromosome name convention (ucsc or ens).", required=True)
    parser_convertChr.add_argument("--allowMissing", "-a" ,action="store_true", help="""When a chromsome name can't be matched between USCS and Ensembl
                                                                                        set it to 'NA' (by default thrown as error).""")
    parser_convertChr.add_argument("--suppressMissing", "-s" ,action="store_true", help="""When a chromsome name can't be matched between USCS and Ensembl
                                                                                           do not report it in the output (by default throws an error).""")
    parser_convertChr.add_argument("--patches", "-p" ,action="store_true", help="""Allows conversion of all patches up to p11 for hg38 and p4 for mm10.
                                                                                   Without this option, if the BED file contains contigs added by a patch
                                                                                   the conversion terminates with an error (unless the -a or -s flags are
                                                                                   present).""")

    parser_convertChr.set_defaults(func=convertChr)
    
 
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
