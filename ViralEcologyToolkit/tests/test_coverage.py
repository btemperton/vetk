from unittest import TestCase

import ViralEcologyToolkit


class TestCoverage(TestCase):
    def testCalculatePctIdentity(self):
        x = ViralEcologyToolkit.calculate_pct_identity("NM:i:11", "MD:Z:0T0G2A13A38G6C0T22A4G4C1")
        self.assertAlmostEqual(x, 89.1, 1)

    def testIsGoodEnoughRead(self):
        line = "H3:D1BKNACXX:1:1213:12925:62082	73	HTVC008M	5052	0	4M2I3M2D91M	=	5052	0	ATTTACCTGTCTGGCTGAACAACGATTAAGAGGTTGTAAGGCATGGGTAGAGAGGGTTATGGCCGAATGGCTGAAGACACTCTATTCAGTTGTTAGTAAG	@@@FDFFFGHHFHFIJIJJJJGIEIEGIIIHIJ<DE?FEEGEHGDEH?BDBCDHEG;@EG>ECEFDFD>>A@C?>C;ACDCDD<CDDCCFECDBCDAC>:	AS:i:-60	XN:i:0	XM:i:8	XO:i:2	XG:i:4	NM:i:12	MD:Z:3G3^CT0A13A38G6C0T22A4G1	YT:Z:UP"
        self.assertFalse(ViralEcologyToolkit.is_good_enough_read(line, 90))
        line = "H3:D1BKNACXX:1:2214:9047:2888	163	uvMED-CGR-U-MedDCM-OCT-S46-C34	37279	8	101M	=	37330	152	ATTAAAACAGTTATGGATCAAGTATCAAGTACACCTATTAAAGAGCTACCTGGAGATGATGGGAATAATAAAGTCTCAAGCTGAAAGAACGGAAAAGAGAA	CCCFFFFFHHHHHJJJJJJJIJHJJJJJIHIJJJJJJJJIJJJJJJJIJJIJJIJIJJJJJJJIJJJIJJJJJGHIJIIHHHHHFFFFFEDCDDDDDDCDC	AS:i:-24	XN:i:0	XM:i:4	XO:i:0	XG:i:0	NM:i:4	MD:Z:58T0C3G13G23	YS:i:-43	YT:Z:CP"
        self.assertTrue(ViralEcologyToolkit.is_good_enough_read(line, 90))

    def testFilterBAM(self):
        ViralEcologyToolkit.extract_good_reads(
            '/Users/bt273/Dropbox/TempertonLab/ViralEcologyToolkit/ViralEcologyToolkit/tests/testfiles/clean_read_mapping_sorted.bam',
            output='/Users/bt273/Dropbox/TempertonLab/ViralEcologyToolkit/ViralEcologyToolkit/tests/testfiles/filtered')
