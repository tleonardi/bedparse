import re

class BEDexception(Exception):
    pass




class bedLine(object):
    """An object to represent a BED 12 line"""
    fields = ("chr", "start", "end", "name", "score", "strand", "cdsStart", "cdsEnd", "color", "nEx", "exLengths", "exStarts")
    def __init__(self, line):
        self.bedType=len(line)
        for n in range(self.bedType):
           self.__dict__[self.fields[n]] = line[n]
      
        # Check bed type
        if(not self.bedType in (3,4,6,12)):
            raise BEDexception("Only BED3,4,6,12 are supported. "+self.name+" is neither.")

        # Validate Start and End
        try:
            self.start=int(self.start)
            self.end=int(self.end)
        except:
            raise BEDexception("Start or End are not an int for transcript "+self.name)
        if(self.start>self.end):
            raise BEDexception("Start is greater than End for transcript "+self.name)
      
        #Validate the strand and set stranded property
        if(self.bedType>=6):
            if(self.strand=="+" or self.strand=="-"):
                self.stranded=1
            elif(self.strand=="" or self.strand == "."):
                self.strand=0
            else:
                raise BEDexception("The strand is not any of '+', '-', '.' or '' for transcript: "+self.name)

        if(self.bedType==12):
            # Validate nEx, and CDS fields
            try:
                self.nEx=int(self.nEx)
            except:    
                raise BEDexception("Number of exons is not an int for transcript "+self.name)
            try:
                self.cdsStart=int(self.cdsStart) 
                self.cdsEnd=int(self.cdsEnd)
            except:
                raise BEDexception("CDSstart or CDSend are not int for transcript "+self.name)
            if(self.cdsStart>self.cdsEnd):
                raise BEDexception("CDSstart is greater than CDSend for transcript "+self.name)
            if(self.cdsStart<self.start or self.cdsEnd>self.end):
                raise BEDexception("The CDS range is bigger than the transcript for transcript "+self.name)
            if(re.search(',$', self.exLengths) == None or re.search(',$', self.exStarts) == None ):
                raise BEDexception("Exon lengths or starts do not end with ',' for transcript "+self.name)
            # Check that number of blocks corresponds to the conent of fields 11 and 12
            if(len(self.exLengths.split(","))-1!=self.nEx):
                raise BEDexception("Exon lengths and number of exons mismatch for transcript "+self.name)
            if(len(self.exStarts.split(","))-1!=self.nEx):
                raise BEDexception("Exon starts and number of exons mismatch for transcript "+self.name)
            # Check that every element of exLengths and exStarts can be coerced to int
            for ex in self.exStarts.split(",")[1:-1]:
                try: int(ex)
                except ValueError:
                    raise BEDexception("Exon starts are not int for transcript "+self.name)
            for ex in self.exLengths.split(",")[1:-1]:
                try: int(ex)
                except ValueError:
                    raise BEDexception("Exon lengths are not int for transcript "+self.name)

            # Remove trailing new line
            self.exStarts=self.exStarts.rstrip()

            # If cds start and end are the same set hasORF to 0
            if(self.cdsStart==self.cdsEnd):
                self.hasORF=0
            else:
                self.hasORF=1

    def __str__(self):
        out=[]
        for key in self.fields[:self.bedType]:
            out.append(self.__dict__[key])
        return str(out)

    def print(self):
        out=[]
        for key in self.fields[:self.bedType]:
            out.append(self.__dict__[key])
        return print(*out, sep="\t")

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def promoter(self, up=500, down=500, strand=1):
        """ Returns a bedLine of the promoters"""
        if(not strand or self.strand=="+"):
            start = self.start-up if self.start-up>0 else 0
            end = self.start+down
        elif(strand and self.strand=="-"):
            start= self.end-down if self.end-down>0 else 0
            end=self.end+up
        else:
            raise BEDexception("Strand not recognised for transcript "+self.name)
        return bedLine([self.chr, start, end, self.name])

    def utr(self, which=None):
        """ Returns the 5p UTR of coding transcripts (i.e. those with a CDS) """
        if(not self.stranded):
            raise BEDexception("UTRs for an unstranded transcript make little sense: "+self.name)
        if(which!=5 and which!=3):
            raise BEDexception("'which' needs to be 3 or 5")

        if(self.hasORF==0 or (self.cdsStart == self.start and self.cdsEnd == self.end)): 
            return None 
        # This block return the first UTR, i.e. the 5'UTR of + transcripts
        # or the 3' UTR of - transcripts
        if((self.strand=="+" and which==5) or (self.strand=="-" and which==3)):
            # The UTR starts and the beginning of the transcript
            start=self.start
            # and it ends at the beginning of the CDS
            end=self.cdsStart
            # This is the UTR end in transcripts coordinates
            relEnd=self.cdsStart-start
            exStarts=[]
            exLens=[]
            oldStarts=self.exStarts.split(",")
            oldLengths=self.exLengths.split(",")
            # Add exons one by one until we pass relEnd
            for i in range(0,self.nEx):
                exStarts.append(oldStarts[i])
                # If the current exons end after relEnd, stop the loop
                if(relEnd < int(oldStarts[i])+int(oldLengths[i])):
                    exLens.append(relEnd-int(oldStarts[i]))
                    nEx=i+1
                    break
                else:
                    exLens.append(int(oldLengths[i]))
        # This block returns the second UTR, i.e the 3'UTR of + transcripts
        # or the 5'UTR of - transcripts
        elif((self.strand=="+" and which==3) or (self.strand=="-" and which==5)):
            # The UTR starts at the end of the CDS
            start=self.cdsEnd
            # and it ends at the end of the transcript
            end=self.end
            # This is the UTR start in transcript coordinates
            relStart=self.cdsEnd-self.start
            exStarts=[]
            exLens=[]
            oldStarts=self.exStarts.split(",")
            oldLengths=self.exLengths.split(",")
            nEx=self.nEx
            # Loop through all exons and skip them until we reach relStart
            for i in range(0,self.nEx):
                # if the current exon ends before relStart, skip it
                if(relStart > int(oldStarts[i])+int(oldLengths[i])):
                    nEx=nEx-1
                    next
                # otherwise if the current exon starts before relStart
                # (i.e. the current exon contain relStart), add it to the
                # list of exons starts and lengths
                elif(int(oldStarts[i])<relStart):
                    exStarts.append("0")
                    exLens.append(int(oldStarts[i])+int(oldLengths[i])-relStart)
                # otherwise (i.e. the current exon is past relStart)
                # add its start and length to the lists
                else:
                    exStarts.append(int(oldStarts[i])-relStart)
                    exLens.append(int(oldLengths[i]))

        # The list of starts and lengths has to end with a comma
        exStarts.append("")
        exLens.append("")
        result=bedLine([self.chr, start, end, self.name, self.score, self.strand, start, start, self.color, nEx, ','.join(str(x) for x in exLens), ','.join(str(x) for x in exStarts)])
        return result


    def cds(self, ignoreCDSonly=False):
        """Return the CDS of a coding transcript. Transcripts that are only CDS are NOT reported."""
        if(not self.stranded):
            raise BEDexception("CDS for an unstranded transcript makes little sense: "+self.name)

        if(self.hasORF==0): 
            return None 
        if(ignoreCDSonly == True and (self.cdsStart == self.start and self.cdsEnd == self.end)):
            return None

        start=self.cdsStart
        end=self.cdsEnd
        exStarts=[]
        exLens=[]
        oldStarts=self.exStarts.split(",")
        oldLengths=self.exLengths.split(",")
        nEx=0
            # Loop through all exons and skip them until we reach relStart
        # This is the relative cds Start
        relStart=self.cdsStart-self.start
        relEnd=self.cdsEnd-self.start
        # Loop through exons and skip them till we reach relStart
        for i in range(0,self.nEx):
            # If the current exons ends before the start of the CDS, skip it
            if(relStart>int(oldStarts[i])+int(oldLengths[i])):
                next
            # Else, if the current exon starts before the CDS
            # add it to the list
            elif(int(oldStarts[i])<relStart):
                exStarts.append('0')
                exLens.append(int(oldStarts[i])+int(oldLengths[i])-relStart)
                nEx=nEx+1
            # otherwise (i.e. the current exon is past relStart)
            # add its start and length to the lists
            else:
                exStarts.append(int(oldStarts[i])-relStart)
                exLens.append(int(oldLengths[i]))
                nEx=nEx+1

            # If the current exon ends after relEnd, stop the loop,
            # remove the last length, and add the correct one
            if(relEnd < int(oldStarts[i])+int(oldLengths[i])):
                curLen=exLens.pop()
                # The final length is the current length (i.e. without the 1st UTR in case the
                # transcript is mono-exonic) minus the difference between the old length and the
                # end of the CDS
                #100           200            300              500
                #---------------|--------------|-----------------
                #---------------@@@@@@@@@@@@@@@@-----------------
                # Rel:   100(relStart)     200(relEnd)
                # length=400(curLen) - (500(old size) - 200) = 100
                exLens.append(curLen-(int(oldStarts[i])+int(oldLengths[i])-relEnd))
                break
        # The list of starts and lengths has to end with a comma
        exStarts.append("")
        exLens.append("")
        result=bedLine([self.chr, start, end, self.name, self.score, self.strand, start, end, self.color, nEx, ','.join(str(x) for x in exLens), ','.join(str(x) for x in exStarts)])
        return result

