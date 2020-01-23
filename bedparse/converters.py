#!/usr/bin/python3
import csv
import re
from bedparse import bedline

def gtf2bed(gtf, extra=[''], filterKey="transcript_biotype", filterType=[''], transcript_feature_name= "transcript"):
    gtfRecords={'exon':list(), 'transcript': list(), 'cds':list()}
    transcripts=dict()
    exons=dict()
    cds=dict()
    extrainfo=dict()
    gtfReader = csv.reader((row for row in gtf if not row.startswith('#')), delimiter="\t")
    for line in gtfReader:
        # Store all transcript lines
        if(line[2]== transcript_feature_name):
            txName=re.sub('.*transcript_id "([^"]+)";.*', "\\1", line[8])
            if(line[6]!="+" and line[6]!="-"):
                raise BEDexception("Transcript with unrecognized strand: "+txName)
            # Start-1 converts from 1-based to 0-based
            transcripts[txName] = [line[0], int(line[3])-1, int(line[4]), txName, 0, line[6]]
            # Parse the extra fields
            if(extra!=['']):
                for field in extra:
                    extrainfo.setdefault(txName, dict())[field]=re.sub('.*'+field+' "?([^"]+)"?;.*', "\\1", line[8])
                    # If no substitution occured, set the extra field to '.'
                    if(extrainfo[txName][field] == line[8]): extrainfo[txName][field] = "."
            if filterType!=[''] and filterKey not in extra:
                extrainfo.setdefault(txName, dict())[filterKey]=re.sub('.*'+filterKey+' "?([^"]+)"?;.*', "\\1", line[8])
                # If no substitution occured, set the extra field to '.'
                if(extrainfo[txName][filterKey] == line[8]): extrainfo[txName][field] = "."
        # Parse exon lines
        if(line[2]=='exon'):
            txName=re.sub('.*transcript_id "([^"]+)";.*', "\\1", line[8])
            if(line[6]!=transcripts[txName][5]):
                raise BEDexception("Exon has different strand from parent transcript: "+txName)
            start=int(line[3])-1
            length=int(line[4])-int(line[3])
            exons.setdefault(txName, []).append([start,length])

        # Start CDS, start and stop codons
        # Any of these features extends the CDS
        if(line[2] in ['CDS', 'start_codon', 'stop_codon']):
            txName=re.sub('.*transcript_id "([^"]+)";.*', "\\1", line[8])
            if(line[6]!=transcripts[txName][5]):
                raise BEDexception(("%s has different strand from parent transcript: %s" % line[2], txName))
            start=int(line[3])
            stop=int(line[4])
            if(txName not in cds.keys()):
                cds[txName] = [start,stop]
            else:
                # Is the current start/stop are up/down compared to
                # the previous one, update cds.
                if(start < cds[txName][0]): cds[txName][0] = start
                if(stop > cds[txName][1]):cds[txName][1] = stop
    gtf.close()

    for transcript in transcripts.keys():
        if(filterType!=['']):
            if extrainfo[transcript][filterKey] not in filterType:
                continue
        if(transcript in cds.keys()):
            cdsStart=int(cds[transcript][0])-1
            cdsEnd = int(cds[transcript][1])
        else:
            cdsStart=transcripts[transcript][2]
            cdsEnd=transcripts[transcript][2]

        transcripts[transcript].append(cdsStart)
        transcripts[transcript].append(cdsEnd)
        transcripts[transcript].append('0')
        # Add the number of exons in field 10
        transcripts[transcript].append(len(exons[transcript]))
        # Sort the [start, length] pairs by start position
        exons[transcript].sort(key=lambda x: x[0])
        starts=""
        lens=""
        for exon in exons[transcript]:
            starts=starts+str(exon[0]-int(transcripts[transcript][1]))+","
            lens=lens+str(exon[1]+1)+","
        transcripts[transcript].append(lens)
        transcripts[transcript].append(starts)
	# Convert to bedline for format validation
        out=list()
        bed = bedline(transcripts[transcript])
        for key in bed._bedline__fields[:bed.bedType]:
            if key=="exLengths": bed.exLengths+=","
            if key=="exStarts": bed.exStarts+=","
            out.append(bed.__dict__[key])
        if(extra!=['']):
            for field in extra:
                out.append(extrainfo[bed.name][field])
        print(*out, sep="\t")
