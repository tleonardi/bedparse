#!/usr/bin/python3
import signal
import argparse
import sys
import csv
from bedClass import * 
# This allows using the program in a pipe
# The program is killed when it receives a sigpipe
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

parser = argparse.ArgumentParser(description="Filters a BED file based on an annotation. BED entries with a name (i.e. col4) that appears in the specified column of the annotation are printed to stdout.")
parser.add_argument("bedfile", type=argparse.FileType('r'), nargs='?', default=sys.stdin,
                    help="Path to the BED file.")

parser.add_argument("-p","--promoter", action="store_true",
        help="Return the promoter. Modifiers: --up, --down, --stranded")
parser.add_argument("--fivePutr",action="store_true",
        help="Return the 5' UTR of coding genes.")
parser.add_argument("--threePutr",action="store_true",
        help="Return the 3' UTR of coding genes.")
parser.add_argument("--cds",action="store_true",
        help="Return the CDS of coding genes.")
parser.add_argument("--filter",action="store_true",
        help="""Filters a BED file based on an annotation.
        BED entries with a name (i.e. col4) that appears
        in the specified column of the annotation are
        printed to stdout. For efficiency reasons this
        command doesn't perform BED validation.""")
parser.add_argument("--unstranded",action="store_true",
        help="Do not consider strands.")
parser.add_argument("--up",type=int, default=500,
        help="Get this many nt upstream of each feature.")
parser.add_argument("--down",type=int, default=500,
        help="Get this many nt downstream of each feature.")
parser.add_argument("--annotation", "-a", type=str,
        help="Path to the annotation file")
parser.add_argument("--column","-c",type=int, default=1,
        help="Column of the annotation file (1-based, default=1)")

args = parser.parse_args()

if(args.promoter+args.fivePutr+args.threePutr+args.filter> 1):
    raise BEDexception("Incompatible options")

if(args.filter):
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
        if(args.promoter):
            bedLine(line.split('\t')).promoter(up=args.up, down=args.down, strand=(not args.unstranded)).print()

        elif(args.fivePutr):
            utr=bedLine(line.split('\t')).utr(which=5)
            if(utr): utr.print()

        elif(args.threePutr):
            utr=bedLine(line.split('\t')).utr(which=3)
            if(utr): utr.print()

        elif(args.cds):
            utr=bedLine(line.split('\t')).cds()
            if(utr): utr.print()

        elif(args.filter):
            #bedline=bedLine(line.split('\t'))
            #if(bedline.name in filterset):
            #    bedline.print()
            if(line.split('\t')[3] in filterset):
                print(line.rstrip())


tsvfile.close()
