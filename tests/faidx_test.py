import pysam
import unittest
import os

DATADIR = "pysam_data"


class TestFastaFile(unittest.TestCase):

    sequences = {
        'chr1':
        b"CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCGTCCATGGCCCAGCATTAGGGAGCTGTGGACCCTGCAGCCTGGCTGTGGGGGCCGCAGTGGCTGAGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCACAGGGGCTTAACCTCTGGTGACTGCCAGAGCTGCTGGCAAGCTAGAGTCCCATTTGGAGCCCCTCTAAGCCGTTCTATTTGTAATGAAAACTATATTTATGCTATTCAGTTCTAAATATAGAAATTGAAACAGCTGTGTTTAGTGCCTTTGTTCAACCCCCTTGCAACAACCTTGAGAACCCCAGGGAATTTGTCAATGTCAGGGAAGGAGCATTTTGTCAGTTACCAAATGTGTTTATTACCAGAGGGATGGAGGGAAGAGGGACGCTGAAGAACTTTGATGCCCTCTTCTTCCAAAGATGAAACGCGTAACTGCGCTCTCATTCACTCCAGCTCCCTGTCACCCAATGGACCTGTGATATCTGGATTCTGGGAAATTCTTCATCCTGGACCCTGAGAGATTCTGCAGCCCAGCTCCAGATTGCTTGTGGTCTGACAGGCTGCAACTGTGAGCCATCACAATGAACAACAGGAAGAAAAGGTCTTTCAAAAGGTGATGTGTGTTCTCATCAACCTCATACACACACATGGTTTAGGGGTATAATACCTCTACATGGCTGATTATGAAAACAATGTTCCCCAGATACCATCCCTGTCTTACTTCCAGCTCCCCAGAGGGAAAGCTTTCAACGCTTCTAGCCATTTCTTTTGGCATTTGCCTTCAGACCCTACACGAATGCGTCTCTACCACAGGGGGCTGCGCGGTTTCCCATCATGAAGCACTGAACTTCCACGTCTCATCTAGGGGAACAGGGAGGTGCACTAATGCGCTCCACGCCCAAGCCCTTCTCACAGTTTCTGCCCCCAGCATGGTTGTACTGGGCAATACATGAGATTATTAGGAAATGCTTTACTGTCATAACTATGAAGAGACTATTGCCAGATGAACCACACATTAATACTATGTTTCTTATCTGCACATTACTACCCTGCAATTAATATAATTGTGTCCATGTACACACGCTGTCCTATGTACTTATCATGACTCTATCCCAAATTCCCAATTACGTCCTATCTTCTTCTTAGGGAAGAACAGCTTAGGTATCAATTTGGTGTTCTGTGTAAAGTCTCAGGGAGCCGTCCGTGTCCTCCCATCTGGCCTCGTCCACACTGGTTCTCTTGAAAGCTTGGGCTGTAATGATGCCCCTTGGCCATCACCCAGTCCCTGCCCCATCTCTTGTAATCTCTCTCCTTTTTGCTGCATCCCTGTCTTCCTCTGTCTTGATTTACTTGTTGTTGGTTTTCTGTTTCTTTGTTTGATTTGGTGGAAGACATAATCCCACGCTTCCTATGGAAAGGTTGTTGGGAGATTTTTAATGATTCCTCAATGTTAAAATGTCTATTTTTGTCTTGACACCCAACTAATATTTGTCTGAGCAAAACAGTCTAGATGAGAGAGAACTTCCCTGGAGGTCTGATGGCGTTTCTCCCTCGTCTTCTTA",
        'chr2':
        b"TTCAAATGAACTTCTGTAATTGAAAAATTCATTTAAGAAATTACAAAATATAGTTGAAAGCTCTAACAATAGACTAAACCAAGCAGAAGAAAGAGGTTCAGAACTTGAAGACAAGTCTCTTATGAATTAACCCAGTCAGACAAAAATAAAGAAAAAAATTTTAAAAATGAACAGAGCTTTCAAGAAGTATGAGATTATGTAAAGTAACTGAACCTATGAGTCACAGGTATTCCTGAGGAAAAAGAAAAAGTGAGAAGTTTGGAAAAACTATTTGAGGAAGTAATTGGGGAAAACCTCTTTAGTCTTGCTAGAGATTTAGACATCTAAATGAAAGAGGCTCAAAGAATGCCAGGAAGATACATTGCAAGACAGACTTCATCAAGATATGTAGTCATCAGACTATCTAAAGTCAACATGAAGGAAAAAAATTCTAAAATCAGCAAGAGAAAAGCATACAGTCATCTATAAAGGAAATCCCATCAGAATAACAATGGGCTTCTCAGCAGAAACCTTACAAGCCAGAAGAGATTGGATCTAATTTTTGGACTTCTTAAAGAAAAAAAAACCTGTCAAACACGAATGTTATGCCCTGCTAAACTAAGCATCATAAATGAAGGGGAAATAAAGTCAAGTCTTTCCTGACAAGCAAATGCTAAGATAATTCATCATCACTAAACCAGTCCTATAAGAAATGCTCAAAAGAATTGTAAAAGTCAAAATTAAAGTTCAATACTCACCATCATAAATACACACAAAAGTACAAAACTCACAGGTTTTATAAAACAATTGAGACTACAGAGCAACTAGGTAAAAAATTAACATTACAACAGGAACAAAACCTCATATATCAATATTAACTTTGAATAAAAAGGGATTAAATTCCCCCACTTAAGAGATATAGATTGGCAGAACAGATTTAAAAACATGAACTAACTATATGCTGTTTACAAGAAACTCATTAATAAAGACATGAGTTCAGGTAAAGGGGTGGAAAAAGATGTTCTACGCAAACAGAAACCAAATGAGAGAAGGAGTAGCTATACTTATATCAGATAAAGCACACTTTAAATCAACAACAGTAAAATAAAACAAAGGAGGTCATCATACAATGATAAAAAGATCAATTCAGCAAGAAGATATAACCATCCTACTAAATACATATGCACCTAACACAAGACTACCCAGATTCATAAAACAAATACTACTAGACCTAAGAGGGATGAGAAATTACCTAATTGGTACAATGTACAATATTCTGATGATGGTTACACTAAAAGCCCATACTTTACTGCTACTCAATATATCCATGTAACAAATCTGCGCTTGTACTTCTAAATCTATAAAAAAATTAAAATTTAACAAAAGTAAATAAAACACATAGCTAAAACTAAAAAAGCAAAAACAAAAACTATGCTAAGTATTGGTAAAGATGTGGGGAAAAAAGTAAACTCTCAAATATTGCTAGTGGGAGTATAAATTGTTTTCCACTTTGGAAAACAATTTGGTAATTTCGTTTTTTTTTTTTTCTTTTCTCTTTTTTTTTTTTTTTTTTTTGCATGCCAGAAAAAAATATTTACAGTAACT",
    }

    def setUp(self):
        self.file = pysam.FastaFile(os.path.join(DATADIR, "ex1.fa"))

    def testFetch(self):
        for id, seq in list(self.sequences.items()):
            self.assertEqual(seq, self.file.fetch(id))
            for x in range(0, len(seq), 10):
                self.assertEqual(seq[x:x + 10], self.file.fetch(id, x, x + 10))
                # test x:end
                self.assertEqual(seq[x:], self.file.fetch(id, x))
                # test 0:x
                self.assertEqual(seq[:x], self.file.fetch(id, None, x))

        # unknown sequence returns ""
        # change: should be an IndexError
        self.assertEqual(b"", self.file.fetch("chr12"))

    def testOutOfRangeAccess(self):
        '''test out of range access.'''
        # out of range access returns an empty string
        for contig, s in self.sequences.items():
            self.assertEqual(self.file.fetch(contig, len(s), len(s) + 1), b"")

        self.assertEqual(self.file.fetch("chr3", 0, 100), b"")

    def testFetchErrors(self):
        self.assertRaises(ValueError, self.file.fetch)
        self.assertRaises(IndexError, self.file.fetch, "chr1", -1, 10)
        self.assertRaises(ValueError, self.file.fetch, "chr1", 20, 10)

        # does not work yet
        # self.assertRaises( KeyError, self.file.fetch, "chrX" )

    def testLength(self):
        self.assertEqual(len(self.file), 2)

    def testSequenceLengths(self):
        self.assertEqual(1575, self.file.get_reference_length("chr1"))
        self.assertEqual(1584, self.file.get_reference_length("chr2"))

    def tearDown(self):
        self.file.close()


class TestFastxFileFastq(unittest.TestCase):

    filetype = pysam.FastxFile
    filename = "faidx_ex1.fq"
    persist = True

    def setUp(self):
        self.file = self.filetype(os.path.join(DATADIR, self.filename),
                                  persist=self.persist)
        self.has_quality = self.filename.endswith('.fq')

    def checkFirst(self, s):
        # test first entry
        self.assertEqual(s.sequence, b"GGGAACAGGGGGGTGCACTAATGCGCTCCACGCCC")
        self.assertEqual(s.name, b"B7_589:1:101:825:28")
        if self.has_quality:
            self.assertEqual(s.quality, b"<<86<<;<78<<<)<;4<67<;<;<74-7;,;8,;")
            self.assertEqual(list(s.get_quality_array()),
                             [ord(x) - 33 for x in s.quality])
        else:
            self.assertEqual(s.quality, None)
            self.assertEqual(s.get_quality_array(), None)
        
    def checkLast(self, s):
        self.assertEqual(s.sequence, b"TAATTGAAAAATTCATTTAAGAAATTACAAAATAT")
        self.assertEqual(s.name, b"EAS56_65:8:64:507:478")
        if self.has_quality:
            self.assertEqual(s.quality, b"<<<<<;<<<<<<<<<<<<<<<;;;<<<;<<8;<;<")
            self.assertEqual(list(s.get_quality_array()),
                             [ord(x) - 33 for x in s.quality])
        else:
            self.assertEqual(s.quality, None)
            self.assertEqual(s.get_quality_array(), None)

    def testCounts(self):
        self.assertEqual(len([x for x in self.file]), 3270)

    def testMissingFile(self):
        self.assertRaises(IOError, self.filetype, "nothere.fq")

    def testSequence(self):
        first = self.file.__next__()
        self.checkFirst(first)
        for last in self.file:
            pass
        self.checkLast(last)

        # test for persistence
        if self.persist:
            self.checkFirst(first)
        else:
            self.checkLast(first)

# Test for backwards compatibility
class TestFastqFileFastq(TestFastxFileFastq):
    filetype = pysam.FastqFile

# Test for backwards compatibility
class TestFastxFileFasta(TestFastxFileFastq):
    filetype = pysam.FastqFile
    filename = "faidx_ex1.fa"

class TestFastxFileFastqStream(TestFastxFileFastq):
    persist = False

if __name__ == "__main__":
    unittest.main()
