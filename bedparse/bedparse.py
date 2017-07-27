#!/usr/bin/python3
import signal
import argparse
import sys
import csv
from pkg_resources import get_distribution
from bedparse import bedline
from bedparse import gtf2bed
# This allows using the program in a pipe
# The program is killed when it receives a sigpipe
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
__version__ = get_distribution('bedparse').version

bedline([1,1,1])
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
    parser_cds.add_argument("--ignoreCDSonly",action="store_true", help="Ignore transcripts that only consist of CDS")
    parser_cds.add_argument("bedfile", type=argparse.FileType('r'), nargs='?', default=sys.stdin, help="Path to the BED file.")
    parser_cds.set_defaults(func=cds)
    
    parser_prom = subparsers.add_parser('promoter', help="Prints the promoters of coding genes.")
    parser_prom.add_argument("--up",type=int, default=500, help="Get this many nt upstream of each feature.")
    parser_prom.add_argument("--down",type=int, default=500, help="Get this many nt downstream of each feature.")
    parser_prom.add_argument("--unstranded",action="store_true", help="Do not consider strands.")
    parser_prom.add_argument("bedfile", type=argparse.FileType('r'), nargs='?', default=sys.stdin, help="Path to the BED file.")
    parser_prom.set_defaults(func=prom)
    
    parser_filter = subparsers.add_parser('filter', 
            help="""Filters a BED file based on an annotation.
            BED entries with a name (i.e. col4) that appears
            in the specified column of the annotation are
            printed to stdout. For efficiency reasons this
            command doesn't perform BED validation.""")
    parser_filter.add_argument("--annotation", "-a", type=str, help="Path to the annotation file", required=True)
    parser_filter.add_argument("--column","-c",type=int, default=1, help="Column of the annotation file (1-based, default=1)")
    parser_filter.set_defaults(func=filter)
    parser_filter.add_argument("bedfile", type=argparse.FileType('r'), nargs='?', default=sys.stdin,
    help="Path to the BED file.")
 
    parser_gtf2bed = subparsers.add_parser('gtf2bed', 
            help="""Converts a GTF file to BED12 format.
            This tool supports the Ensembl GTF format.
            The GTF file must contain 'transcript' and 'exon' 
            features in field 3. If the GTF file also annotates
            'CDS' 'start_codon' or 'stop_codon' these are used
            to annotate the thickStart and thickEnd in the BED
            file.""")
    parser_gtf2bed.add_argument("gtf", type=argparse.FileType('r'), nargs='?', default=sys.stdin, help="Path to the GTF file")
    parser_gtf2bed.add_argument("--extraFields",type=str, default='', help="Comma separated list of extra GTF fields to be added after col 12 (e.g. gene_id,gene_name).")
    parser_gtf2bed.set_defaults(func=lambda args: gtf2bed(args.gtf, args.extraFields.split(',')))
    
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
