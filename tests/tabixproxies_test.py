import unittest
import pysam
import os
import sys
import re
import copy
import gzip
from TestUtils import load_and_convert

DATADIR = 'tabix_data'


class TestParser(unittest.TestCase):

    filename = os.path.join(DATADIR, "example.gtf.gz")

    def setUp(self):

        self.tabix = pysam.TabixFile(self.filename)
        self.compare = load_and_convert(self.filename)

    def tearDown(self):
        self.tabix.close()

    def testRead(self):

        for x, r in enumerate(self.tabix.fetch(parser=pysam.asTuple())):
            c = self.compare[x]
            self.assertEqual(c, list(r))
            self.assertEqual(len(c), len(r))

            # test indexing
            for y in range(0, len(r)):
                self.assertEqual(c[y], r[y])

            # test slicing access
            for y in range(0, len(r) - 1):
                for cc in range(y + 1, len(r)):
                    self.assertEqual(c[y:cc],
                                     r[y:cc])
            self.assertEqual("\t".join(map(str, c)),
                             str(r))

    def testWrite(self):

        for x, r in enumerate(self.tabix.fetch(parser=pysam.asTuple())):
            self.assertEqual(self.compare[x], list(r))
            c = list(r)
            for y in range(len(r)):
                r[y] = "test_%05i" % y
                c[y] = "test_%05i" % y
            self.assertEqual([x for x in c], list(r))
            self.assertEqual("\t".join(c), str(r))
            # check second assignment
            for y in range(len(r)):
                r[y] = "test_%05i" % y
            self.assertEqual([x for x in c], list(r))
            self.assertEqual("\t".join(c), str(r))

    def testUnset(self):
        for x, r in enumerate(self.tabix.fetch(parser=pysam.asTuple())):
            self.assertEqual(self.compare[x], list(r))
            c = list(r)
            e = list(r)
            for y in range(len(r)):
                r[y] = None
                c[y] = None
                e[y] = ""
                self.assertEqual(c, list(r))
                self.assertEqual("\t".join(e), str(r))

    def testIteratorCompressed(self):
        '''test iteration from compressed file.'''
        with gzip.open(self.filename) as infile:
            for x, r in enumerate(pysam.tabix_iterator(
                    infile, pysam.asTuple())):
                self.assertEqual(self.compare[x], list(r))
                self.assertEqual(len(self.compare[x]), len(r))

                # test indexing
                for c in range(0, len(r)):
                    self.assertEqual(self.compare[x][c], r[c])

                # test slicing access
                for c in range(0, len(r) - 1):
                    for cc in range(c + 1, len(r)):
                        self.assertEqual(self.compare[x][c:cc],
                                         r[c:cc])

    def testIteratorUncompressed(self):
        '''test iteration from uncompressed file.'''
        tmpfilename = 'tmp_testIteratorUncompressed'
        with gzip.open(self.filename, "rb") as infile, \
             open(tmpfilename, "wb") as outfile:
            outfile.write(infile.read())

        with open(tmpfilename) as infile:
            for x, r in enumerate(pysam.tabix_iterator(
                    infile, pysam.asTuple())):
                self.assertEqual(self.compare[x], list(r))
                self.assertEqual(len(self.compare[x]), len(r))

                # test indexing
                for c in range(0, len(r)):
                    self.assertEqual(self.compare[x][c], r[c])

                # test slicing access
                for c in range(0, len(r) - 1):
                    for cc in range(c + 1, len(r)):
                        self.assertEqual(self.compare[x][c:cc],
                                         r[c:cc])

        os.unlink(tmpfilename)

    def testCopy(self):
        a = self.tabix.fetch(parser=pysam.asTuple()).next()
        b = copy.copy(a)
        self.assertEqual(a, b)

        a = self.tabix.fetch(parser=pysam.asGTF()).next()
        b = copy.copy(a)
        self.assertEqual(a, b)


class TestGTF(TestParser):

    parser = pysam.asGTF

    def testRead(self):

        for x, r in enumerate(self.tabix.fetch(parser=self.parser())):
            c = self.compare[x]
            self.assertEqual(len(c), len(r))
            self.assertEqual(list(c), list(r))
            self.assertEqual(c, str(r).split("\t"))
            self.assertTrue(r.gene_id.startswith("ENSG"))
            if r.feature != 'gene':
                self.assertTrue(r.transcript_id.startswith("ENST"))
            self.assertEqual(c[0], r.contig)
            self.assertEqual("\t".join(map(str, c)),
                             str(r))

    def testSetting(self):

        r = self.tabix.fetch(parser=self.parser()).next()

        r.contig = r.contig + "_test_contig"          
        r.source = r.source + "_test_source"
        r.feature = r.feature + "_test_feature"
        r.start += 10
        r.end += 10
        r.score = 20
        r.strand = "+"
        r.frame = 0
        r.attributes = 'gene_id "0001";'
        r.transcript_id = "0002"
        sr = str(r)
        self.assertTrue("_test_contig" in sr)
        self.assertTrue("_test_source" in sr)
        self.assertTrue("_test_feature" in sr)
        self.assertTrue("gene_id \"0001\"" in sr)
        self.assertTrue("transcript_id \"0002\"" in sr)

    def test_added_attribute_is_output(self):
        r = self.tabix.fetch(parser=self.parser()).next()

        r.new_int_attribute = 12
        self.assertTrue("new_int_attribute 12" in str(r).split("\t")[8])

        r.new_float_attribute = 12.0
        self.assertTrue("new_float_attribute 12.0" in str(r).split("\t")[8])

        r.new_text_attribute = "abc"
        self.assertTrue("new_text_attribute \"abc\"" in str(r).split("\t")[8])

    def test_setting_start_is_one_based(self):
        
        r = self.tabix.fetch(parser=self.parser()).next()
        r.start = 1800
        self.assertEqual(r.start, 1800)
        self.assertEqual(str(r).split("\t")[3], "1801")

    def test_setting_end_is_one_based(self):
        
        r = self.tabix.fetch(parser=self.parser()).next()
        r.end = 2100
        self.assertEqual(r.end, 2100)
        self.assertEqual(str(r).split("\t")[4], "2100")

    def test_setting_frame_to_none_produces_dot(self):

        r = self.tabix.fetch(parser=self.parser()).next()
        r.frame = None
        self.assertEqual(str(r).split("\t")[7], ".")

        r.frame = 2
        self.assertEqual(str(r).split("\t")[7], "2")

        r = self.tabix.fetch(parser=self.parser()).next()
        r.frame = "."
        self.assertEqual(r.frame, None)
        self.assertEqual(str(r).split("\t")[7], ".")

    def test_setting_source_to_none_produces_dot(self):

        r = self.tabix.fetch(parser=self.parser()).next()
        r.source = None
        self.assertEqual(str(r).split("\t")[1], ".")

        r.source = "source"
        self.assertEqual(str(r).split("\t")[1], "source")

        r = self.tabix.fetch(parser=self.parser()).next()
        r.source = "."
        self.assertEqual(r.source, None)
        self.assertEqual(str(r).split("\t")[1], ".")

    def test_setting_feature_to_none_produces_dot(self):

        r = self.tabix.fetch(parser=self.parser()).next()
        r.feature = None
        self.assertEqual(str(r).split("\t")[2], ".")

        r.feature = "feature"
        self.assertEqual(str(r).split("\t")[2], "feature")

        r = self.tabix.fetch(parser=self.parser()).next()
        r.feature = "."
        self.assertEqual(r.feature, None)
        self.assertEqual(str(r).split("\t")[2], ".")

    def test_setting_strand_to_none_produces_dot(self):

        r = self.tabix.fetch(parser=self.parser()).next()
        r.strand = None
        self.assertEqual(str(r).split("\t")[6], ".")

        r.strand = "-"
        self.assertEqual(str(r).split("\t")[6], "-")

        r = self.tabix.fetch(parser=self.parser()).next()
        r.strand = "."
        self.assertEqual(r.strand, None)
        self.assertEqual(str(r).split("\t")[6], ".")

    def test_setting_score_to_none_produces_dot(self):

        r = self.tabix.fetch(parser=self.parser()).next()
        r.score = None
        self.assertEqual(str(r).split("\t")[5], ".")

        r.score = 12.0
        self.assertEqual(str(r).split("\t")[5], "12.0")

        r.score = -12.0
        self.assertEqual(str(r).split("\t")[5], "-12.0")

        r = self.tabix.fetch(parser=self.parser()).next()
        r.score = "."
        self.assertEqual(r.score, None)
        self.assertEqual(str(r).split("\t")[5], ".")

        r.score = 12
        self.assertEqual(str(r).split("\t")[5], "12")

        r.score = -12
        self.assertEqual(str(r).split("\t")[5], "-12")


class TestGFF3(TestGTF):

    parser = pysam.asGFF3
    filename = os.path.join(DATADIR, "example.gff3.gz")

    def testRead(self):
        for x, r in enumerate(self.tabix.fetch(parser=self.parser())):
            c = self.compare[x]
            self.assertEqual(len(c), len(r))
            self.assertEqual(list(c), list(r))
            self.assertEqual(c, str(r).split("\t"))
            self.assertEqual(c[0], r.contig)
            self.assertEqual("\t".join(map(str, c)),
                             str(r))
            self.assertTrue(r.ID.startswith("MI00"))

    def testSetting(self):

        for r in self.tabix.fetch(parser=self.parser()):
            r.contig = r.contig + "_test_contig"          
            r.source = "test_source"
            r.feature = "test_feature"
            r.start += 10
            r.end += 10
            r.score = 20
            r.strand = "+"
            r.frame = 0
            r.ID="test"
            sr = str(r)
            self.assertTrue("test_contig" in sr)
            self.assertTrue("test_source" in sr)
            self.assertTrue("test_feature" in sr)
            self.assertTrue("ID=test" in sr)
            
    def test_added_attribute_is_output(self):
        r = self.tabix.fetch(parser=self.parser()).next()

        r.new_int_attribute = 12
        self.assertTrue("new_int_attribute=12" in str(r).split("\t")[8])

        r.new_float_attribute = 12.0
        self.assertTrue("new_float_attribute=12.0" in str(r).split("\t")[8])

        r.new_text_attribute = "abc"
        self.assertTrue("new_text_attribute=abc" in str(r).split("\t")[8])


if __name__ == "__main__":
    unittest.main()
