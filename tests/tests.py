import unittest
import bedParse.bedLine

class KnownValues(unittest.TestCase):
    known_promoters = ( 
            (["chr1", 1000, 2000, "Name", 0, "+"], ["chr1", 500, 1500, "Name"]),
            (["chr1", 1000, 2000, "Name", 0, "-"], ["chr1", 1500, 2500, "Name"]),
            (["chr1", 100, 200, "Name", 0, "+"], ["chr1", 0, 600, "Name"]),
            (["chr1", 100, 200, "Name", 0, "-"], ["chr1", 0, 700, "Name"])
            )
    known_promoters100 = ( 
            (["chr1", 1000, 2000, "Name", 0, "+"], ["chr1", 900, 1100, "Name"]),
            (["chr1", 1000, 2000, "Name", 0, "-"], ["chr1", 1900, 2100, "Name"]),
            (["chr1", 50, 100, "Name", 0, "+"], ["chr1", 0, 150, "Name"]),
            (["chr1", 10, 80, "Name", 0, "-"], ["chr1", 0, 180, "Name"])
            )
    known_promotersUnstranded= ( 
            (["chr1", 1000, 2000, "Name", 0, "+"], ["chr1", 500, 1500, "Name"]),
            (["chr1", 1000, 2000, "Name", 0, "-"], ["chr1", 500, 1500, "Name"]),
            (["chr1", 1000, 2000, "Name", 0, "."], ["chr1", 500, 1500, "Name"])
            )

    badBed= (
            ["chr1", 1000, 200, "Name", 0, "+"],
            ["chr1", 1000, "a", "Name", 0, "+"],
            ["chr1", "a", 2000, "Name", 0, "+"],
            ["chr1", "1000", 2000, "Name", 0, "a"],
            ["chr1", 1000, 2000, "Name", 0, "+", 1000, 1000],
            ["chr1", 1000, 2000, "Name", 0, "+", 1000, 1000, ".", "a", "10,10,10,", "0,100,200,"],
            ["chr1", 1000, 2000, "Name", 0, "+", "a", 1000, ".", 3, "10,10,10,", "0,100,200,"],
            ["chr1", 1000, 2000, "Name", 0, "+", 1000, "a", ".", 3, "10,10,10,", "0,100,200,"],
            ["chr1", 1000, 2000, "Name", 0, "+", 1000, 900, ".", 3, "10,10,10,", "0,100,200,"],
            ["chr1", 1000, 2000, "Name", 0, "+", 100, 900, ".", 3, "10,10,10,", "0,100,200,"],
            ["chr1", 1000, 2000, "Name", 0, "+", 1500, 2001, ".", 3, "10,10,10,", "0,100,200,"],
            ["chr1", 1000, 2000, "Name", 0, "+", 1500, 2000, ".", 3, "10,10,10", "0,100,200,"],
            ["chr1", 1000, 2000, "Name", 0, "+", 1500, 2000, ".", 3, "10,10,10,", "0,100,200"],
            ["chr1", 1000, 2000, "Name", 0, "+", 1500, 2000, ".", 3, "10,10,10,10,", "0,100,200,"],
            ["chr1", 1000, 2000, "Name", 0, "+", 1500, 2000, ".", 3, "10,10,", "0,100,200,"],
            ["chr1", 1000, 2000, "Name", 0, "+", 1500, 2000, ".", 3, "10,10,10,", "0,100,"],
            ["chr1", 1000, 2000, "Name", 0, "+", 1500, 2000, ".", 3, "10,10,10,", "0,100,200,300,"],
            ["chr1", 1000, 2000, "Name", 0, "+", 1500, 2000, ".", 3, "10,10,a,", "0,100,200,"],
            ["chr1", 1000, 2000, "Name", 0, "+", 1500, 2000, ".", 3, "10,10,10,", "0,100,b,"]
            )

    known_5pUTRs =(
            (["chr1", 100, 420, "Name", 0, "+", 210, 310, ".", 4, "20,20,20,20,", "0,100,200,300,"], ["chr1", 100, 210, "Name", 0, "+", 100,100, ".", 2, "20,10,", "0,100,"]),
            (["chr1", 100, 500, "Name", 0, "+", 200, 300, ".", 1, "400,", "0,"], ["chr1", 100, 200, "Name", 0, "+", 100,100, ".", 1, "100,", "0,"]),
            (["chr1", 100, 500, "Name", 0, "-", 200, 300, ".", 1, "400,", "0,"], ["chr1", 300, 500, "Name", 0, "-", 300,300, ".", 1, "200,", "0,"]),
            (["chr1", 100, 420, "Name", 0, "-", 210, 310, ".", 4, "20,20,20,20,", "0,100,200,300,"], ["chr1", 310, 420, "Name", 0, "-", 310,310, ".", 2, "10,20,", "0,90,"])
            )

    known_3pUTRs =(
            (["chr1", 100, 420, "Name", 0, "+", 210, 310, ".", 4, "20,20,20,20,", "0,100,200,300,"], ["chr1", 310, 420, "Name", 0, "+", 310,310, ".", 2, "10,20,", "0,90,"]),
            (["chr1", 100, 500, "Name", 0, "-", 200, 300, ".", 1, "400,", "0,"], ["chr1", 100, 200, "Name", 0, "-", 100,100, ".", 1, "100,", "0,"]),
            (["chr1", 100, 500, "Name", 0, "+", 200, 300, ".", 1, "400,", "0,"], ["chr1", 300, 500, "Name", 0, "+", 300,300, ".", 1, "200,", "0,"]),
            (["chr1", 100, 420, "Name", 0, "-", 210, 310, ".", 4, "20,20,20,20,", "0,100,200,300,"], ["chr1", 100, 210, "Name", 0, "-", 100,100, ".", 2, "20,10,", "0,100,"]),
            )
    
    known_CDSs =(
            (["chr1", 100, 420, "Name", 0, "+", 210, 310, ".", 4, "20,20,20,20,", "0,100,200,300,"], ["chr1", 210, 310, "Name", 0, "+", 210, 310, ".", 2, "10,10,", "0,90,"]),
            (["chr1", 100, 420, "Name", 0, "-", 210, 310, ".", 4, "20,20,20,20,", "0,100,200,300,"], ["chr1", 210, 310, "Name", 0, "-", 210, 310, ".", 2, "10,10,", "0,90,"]),
            (["chr1", 100, 500, "Name", 0, "-", 200, 300, ".", 1, "400,", "0,"], ["chr1", 200, 300, "Name", 0, "-", 200,300, ".", 1, "100,", "0,"]),
            (["chr1", 100, 500, "Name", 0, "+", 200, 300, ".", 1, "400,", "0,"], ["chr1", 200, 300, "Name", 0, "+", 200,300, ".", 1, "100,", "0,"])
            )
 
    def test_promoter(self):
        '''promoters() should return correct promoters with known input'''
        for ((bed), (prom)) in self.known_promoters:
            result = bedParse.bedLine(bed).promoter()
            self.assertEqual(result, bedParse.bedLine(prom))

    def test_promoter100(self):
        '''promoters() should return correct promoters100 with known input'''
        for ((bed), (prom)) in self.known_promoters100:
            result = bedParse.bedLine(bed).promoter(up=100, down=100)
            self.assertEqual(result, bedParse.bedLine(prom))

    def test_promoterUnstranded(self):
        '''promoters() should return correct promotersUnstranded with known input'''
        for ((bed), (prom)) in self.known_promotersUnstranded:
            result = bedParse.bedLine(bed).promoter(strand=0)
            self.assertEqual(result, bedParse.bedLine(prom))

    def test_badBed(self):
        '''Bad BED lines should throw and exception'''
        for bed in self.badBed:
            self.assertRaises(bedParse.BEDexception, bedParse.bedLine, bed)

    def test_5pUTRs(self):
        '''fivePutr should return correct UTR for know cases'''
        for ((bed), (utr)) in self.known_5pUTRs:
            result = bedParse.bedLine(bed).utr(which=5)
            self.assertEqual(result, bedParse.bedLine(utr))

    def test_3pUTRs(self):
        '''threePutr should return correct UTR for know cases'''
        for ((bed), (utr)) in self.known_3pUTRs:
            result = bedParse.bedLine(bed).utr(which=3)
            self.assertEqual(result, bedParse.bedLine(utr))

    def test_CDSs(self):
        '''cds() should return correct CDSs for know cases'''
        for ((bed), (cds)) in self.known_CDSs:
            result = bedParse.bedLine(bed).cds()
            self.assertEqual(result, bedParse.bedLine(cds))


if __name__ == '__main__':
    unittest.main()
