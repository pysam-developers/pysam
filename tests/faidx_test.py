import pysam
import unittest
import os
import gzip
import shutil

from TestUtils import checkURL

DATADIR = "pysam_data"


class TestFastaFile(unittest.TestCase):

    sequences = {
        'chr1':
        "CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCGTCCATGGCCCAGCATTAGGGAGCTGTGGACCCTGCAGCCTGGCTGTGGGGGCCGCAGTGGCTGAGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCACAGGGGCTTAACCTCTGGTGACTGCCAGAGCTGCTGGCAAGCTAGAGTCCCATTTGGAGCCCCTCTAAGCCGTTCTATTTGTAATGAAAACTATATTTATGCTATTCAGTTCTAAATATAGAAATTGAAACAGCTGTGTTTAGTGCCTTTGTTCAACCCCCTTGCAACAACCTTGAGAACCCCAGGGAATTTGTCAATGTCAGGGAAGGAGCATTTTGTCAGTTACCAAATGTGTTTATTACCAGAGGGATGGAGGGAAGAGGGACGCTGAAGAACTTTGATGCCCTCTTCTTCCAAAGATGAAACGCGTAACTGCGCTCTCATTCACTCCAGCTCCCTGTCACCCAATGGACCTGTGATATCTGGATTCTGGGAAATTCTTCATCCTGGACCCTGAGAGATTCTGCAGCCCAGCTCCAGATTGCTTGTGGTCTGACAGGCTGCAACTGTGAGCCATCACAATGAACAACAGGAAGAAAAGGTCTTTCAAAAGGTGATGTGTGTTCTCATCAACCTCATACACACACATGGTTTAGGGGTATAATACCTCTACATGGCTGATTATGAAAACAATGTTCCCCAGATACCATCCCTGTCTTACTTCCAGCTCCCCAGAGGGAAAGCTTTCAACGCTTCTAGCCATTTCTTTTGGCATTTGCCTTCAGACCCTACACGAATGCGTCTCTACCACAGGGGGCTGCGCGGTTTCCCATCATGAAGCACTGAACTTCCACGTCTCATCTAGGGGAACAGGGAGGTGCACTAATGCGCTCCACGCCCAAGCCCTTCTCACAGTTTCTGCCCCCAGCATGGTTGTACTGGGCAATACATGAGATTATTAGGAAATGCTTTACTGTCATAACTATGAAGAGACTATTGCCAGATGAACCACACATTAATACTATGTTTCTTATCTGCACATTACTACCCTGCAATTAATATAATTGTGTCCATGTACACACGCTGTCCTATGTACTTATCATGACTCTATCCCAAATTCCCAATTACGTCCTATCTTCTTCTTAGGGAAGAACAGCTTAGGTATCAATTTGGTGTTCTGTGTAAAGTCTCAGGGAGCCGTCCGTGTCCTCCCATCTGGCCTCGTCCACACTGGTTCTCTTGAAAGCTTGGGCTGTAATGATGCCCCTTGGCCATCACCCAGTCCCTGCCCCATCTCTTGTAATCTCTCTCCTTTTTGCTGCATCCCTGTCTTCCTCTGTCTTGATTTACTTGTTGTTGGTTTTCTGTTTCTTTGTTTGATTTGGTGGAAGACATAATCCCACGCTTCCTATGGAAAGGTTGTTGGGAGATTTTTAATGATTCCTCAATGTTAAAATGTCTATTTTTGTCTTGACACCCAACTAATATTTGTCTGAGCAAAACAGTCTAGATGAGAGAGAACTTCCCTGGAGGTCTGATGGCGTTTCTCCCTCGTCTTCTTA",
        'chr2':
        "TTCAAATGAACTTCTGTAATTGAAAAATTCATTTAAGAAATTACAAAATATAGTTGAAAGCTCTAACAATAGACTAAACCAAGCAGAAGAAAGAGGTTCAGAACTTGAAGACAAGTCTCTTATGAATTAACCCAGTCAGACAAAAATAAAGAAAAAAATTTTAAAAATGAACAGAGCTTTCAAGAAGTATGAGATTATGTAAAGTAACTGAACCTATGAGTCACAGGTATTCCTGAGGAAAAAGAAAAAGTGAGAAGTTTGGAAAAACTATTTGAGGAAGTAATTGGGGAAAACCTCTTTAGTCTTGCTAGAGATTTAGACATCTAAATGAAAGAGGCTCAAAGAATGCCAGGAAGATACATTGCAAGACAGACTTCATCAAGATATGTAGTCATCAGACTATCTAAAGTCAACATGAAGGAAAAAAATTCTAAAATCAGCAAGAGAAAAGCATACAGTCATCTATAAAGGAAATCCCATCAGAATAACAATGGGCTTCTCAGCAGAAACCTTACAAGCCAGAAGAGATTGGATCTAATTTTTGGACTTCTTAAAGAAAAAAAAACCTGTCAAACACGAATGTTATGCCCTGCTAAACTAAGCATCATAAATGAAGGGGAAATAAAGTCAAGTCTTTCCTGACAAGCAAATGCTAAGATAATTCATCATCACTAAACCAGTCCTATAAGAAATGCTCAAAAGAATTGTAAAAGTCAAAATTAAAGTTCAATACTCACCATCATAAATACACACAAAAGTACAAAACTCACAGGTTTTATAAAACAATTGAGACTACAGAGCAACTAGGTAAAAAATTAACATTACAACAGGAACAAAACCTCATATATCAATATTAACTTTGAATAAAAAGGGATTAAATTCCCCCACTTAAGAGATATAGATTGGCAGAACAGATTTAAAAACATGAACTAACTATATGCTGTTTACAAGAAACTCATTAATAAAGACATGAGTTCAGGTAAAGGGGTGGAAAAAGATGTTCTACGCAAACAGAAACCAAATGAGAGAAGGAGTAGCTATACTTATATCAGATAAAGCACACTTTAAATCAACAACAGTAAAATAAAACAAAGGAGGTCATCATACAATGATAAAAAGATCAATTCAGCAAGAAGATATAACCATCCTACTAAATACATATGCACCTAACACAAGACTACCCAGATTCATAAAACAAATACTACTAGACCTAAGAGGGATGAGAAATTACCTAATTGGTACAATGTACAATATTCTGATGATGGTTACACTAAAAGCCCATACTTTACTGCTACTCAATATATCCATGTAACAAATCTGCGCTTGTACTTCTAAATCTATAAAAAAATTAAAATTTAACAAAAGTAAATAAAACACATAGCTAAAACTAAAAAAGCAAAAACAAAAACTATGCTAAGTATTGGTAAAGATGTGGGGAAAAAAGTAAACTCTCAAATATTGCTAGTGGGAGTATAAATTGTTTTCCACTTTGGAAAACAATTTGGTAATTTCGTTTTTTTTTTTTTCTTTTCTCTTTTTTTTTTTTTTTTTTTTGCATGCCAGAAAAAAATATTTACAGTAACT",
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

        # unknown sequence raises IndexError
        self.assertRaises(KeyError, self.file.fetch, "chr12")

    def testOutOfRangeAccess(self):
        '''test out of range access.'''
        # out of range access returns an empty string
        for contig, s in self.sequences.items():
            self.assertEqual(self.file.fetch(contig, len(s), len(s) + 1), "")

    def testFetchErrors(self):
        self.assertRaises(ValueError, self.file.fetch)
        self.assertRaises(ValueError, self.file.fetch, "chr1", -1, 10)
        self.assertRaises(ValueError, self.file.fetch, "chr1", 20, 10)
        self.assertRaises(KeyError, self.file.fetch, "chr3", 0, 100)

    def testLength(self):
        self.assertEqual(len(self.file), 2)

    def testSequenceLengths(self):
        self.assertEqual(1575, self.file.get_reference_length("chr1"))
        self.assertEqual(1584, self.file.get_reference_length("chr2"))

    def tearDown(self):
        self.file.close()


class TestFastaFilePathIndex(unittest.TestCase):

    filename = os.path.join(DATADIR, "ex1.fa")

    def testGarbageIndex(self):
        self.assertRaises(NotImplementedError,
                          pysam.FastaFile,
                          self.filename,
                          filepath_index="garbage.fa.fai")
        return

        self.assertRaises(ValueError,
                          pysam.FastaFile,
                          self.filename,
                          filepath_index="garbage.fa.fai")

    def testOpenWithoutIndex(self):
        faidx = pysam.FastaFile(self.filename)
        faidx.close()

    def testOpenWithStandardIndex(self):
        self.assertRaises(NotImplementedError,
                          pysam.FastaFile,
                          self.filename,
                          filepath_index=self.filename + ".fai")
        return

        faidx = pysam.FastaFile(self.filename,
                                filepath_index=self.filename + ".fai")
        faidx.close()

    def testOpenWithOtherIndex(self):
        return
        tmpfilename = "tmp_" + os.path.basename(self.filename)
        shutil.copyfile(self.filename, tmpfilename)
        faidx = pysam.FastaFile(tmpfilename,
                                filepath_index=self.filename + ".fai")
        faidx.close()
        # index should not be auto-generated
        self.assertFalse(os.path.exists(tmpfilename + ".fai"))
        os.unlink(tmpfilename)

class TestFastaFilePathIndexCompressed(TestFastaFilePathIndex):
    
    filename = os.path.join(DATADIR, "ex1.fa.gz")


class TestFastxFileFastq(unittest.TestCase):

    filetype = pysam.FastxFile
    filename = "faidx_ex1.fq"
    persist = True

    def setUp(self):
        self.file = self.filetype(os.path.join(DATADIR, self.filename),
                                  persist=self.persist)
        self.has_quality = self.filename.endswith('.fq')

    def tearDown(self):
        self.file.close()

    def checkFirst(self, s):
        # test first entry
        self.assertEqual(s.sequence, "GGGAACAGGGGGGTGCACTAATGCGCTCCACGCCC")
        self.assertEqual(s.name, "B7_589:1:101:825:28")
        if self.has_quality:
            self.assertEqual(s.quality, "<<86<<;<78<<<)<;4<67<;<;<74-7;,;8,;")
            self.assertEqual(list(s.get_quality_array()),
                             [ord(x) - 33 for x in s.quality])
            self.assertEqual(str(s),
                             "@B7_589:1:101:825:28\n"
                             "GGGAACAGGGGGGTGCACTAATGCGCTCCACGCCC\n"
                             "+\n"
                             "<<86<<;<78<<<)<;4<67<;<;<74-7;,;8,;")

        else:
            self.assertEqual(s.quality, None)
            self.assertEqual(s.get_quality_array(), None)
            self.assertEqual(str(s),
                             ">B7_589:1:101:825:28\n"
                             "GGGAACAGGGGGGTGCACTAATGCGCTCCACGCCC")

    def checkLast(self, s):
        self.assertEqual(s.sequence, "TAATTGAAAAATTCATTTAAGAAATTACAAAATAT")
        self.assertEqual(s.name, "EAS56_65:8:64:507:478")
        if self.has_quality:
            self.assertEqual(s.quality, "<<<<<;<<<<<<<<<<<<<<<;;;<<<;<<8;<;<")
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

    def testManager(self):
        with self.filetype(os.path.join(DATADIR, self.filename),
                           persist=self.persist) as inf:
            first = inf.__next__()
            self.checkFirst(first)
            for last in inf:
                pass
            self.checkLast(last)

        self.assertEqual(inf.closed, True)


# Test for backwards compatibility
class TestFastqFileFastq(TestFastxFileFastq):
    filetype = pysam.FastqFile


# Test for backwards compatibility
class TestFastxFileFasta(TestFastxFileFastq):
    filetype = pysam.FastqFile
    filename = "faidx_ex1.fa"


class TestFastxFileFastqStream(TestFastxFileFastq):
    persist = False


class TestFastxFileWithEmptySequence(unittest.TestCase):
    """see issue 204:

    iteration over fastq file with empty sequence stops prematurely
    """

    filetype = pysam.FastxFile
    filename = "faidx_empty_seq.fq.gz"

    def testIteration(self):
        fn = os.path.join(DATADIR, self.filename)

        with gzip.open(fn) as inf:
            ref_num = len(list(inf)) / 4

        with self.filetype(fn) as f:
            l = len(list(f))
        self.assertEqual(ref_num, l)


class TestRemoteFileFTP(unittest.TestCase):
    '''test remote access.
    '''

    url = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"


    def testFTPView(self):
        if not checkURL(self.url):
            return
        with pysam.Fastafile(self.url) as f:
            self.assertEqual(
                len(f.fetch("chr1", 0, 1000)),
                1000)


if __name__ == "__main__":
    unittest.main()
