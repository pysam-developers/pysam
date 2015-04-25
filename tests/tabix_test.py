#!/usr/bin/env python
'''unit testing code for pysam.

Execute in the :file:`tests` directory as it requires the Makefile
and data files located there.
'''

import sys
import os
import shutil
import gzip
import pysam
import unittest
import glob
import re

DATADIR = 'tabix_data'

IS_PYTHON3 = sys.version_info[0] >= 3


def myzip_open(infile, mode="r"):
    '''open compressed file and decode.'''

    def _convert(f):
        for l in f:
            yield l.decode("ascii")

    if IS_PYTHON3:
        if mode == "r":
            return _convert(gzip.open(infile, "r"))
    else:
        return gzip.open(mode)


def loadAndConvert(infile, encode=True):
    '''load data from infile and convert all fields to string.

    Infile can be either plain or compressed (ending in .gz).
    '''
    data = []
    if infile.endswith(".gz"):
        for line in gzip.open(infile):
            line = line.decode("ascii")
            if line.startswith("#"):
                continue
            d = line.strip().split("\t")
            data.append(d)
    else:
        with open(infile) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                d = line.strip().split("\t")
                data.append(d)

    return data


def splitToBytes(s):
    '''split string and return list of bytes.'''
    return [x.encode("ascii") for x in s.split("\t")]


def checkBinaryEqual(filename1, filename2):
    '''return true if the two files are binary equal.'''
    if os.path.getsize(filename1) != os.path.getsize(filename2):
        return False

    infile1 = open(filename1, "rb")
    infile2 = open(filename2, "rb")

    d1, d2 = infile1.read(), infile2.read()
    found = False
    for c1, c2 in zip(d1, d2):
        if c1 != c2:
            break
    else:
        found = True

    infile1.close()
    infile2.close()
    return found


class TestIndexing(unittest.TestCase):
    filename = os.path.join(DATADIR, "example.gtf.gz")
    filename_idx = os.path.join(DATADIR, "example.gtf.gz.tbi")

    def setUp(self):

        self.tmpfilename = "tmp_%i.gtf.gz" % id(self)
        shutil.copyfile(self.filename, self.tmpfilename)

    def testIndexPreset(self):
        '''test indexing via preset.'''

        pysam.tabix_index(self.tmpfilename, preset="gff")
        checkBinaryEqual(self.tmpfilename + ".tbi", self.filename_idx)

    def tearDown(self):
        os.unlink(self.tmpfilename)
        os.unlink(self.tmpfilename + ".tbi")


class TestCompression(unittest.TestCase):
    filename = os.path.join(DATADIR, "example.gtf.gz")
    filename_idx = os.path.join(DATADIR, "example.gtf.gz.tbi")
    preset = "gff"

    def setUp(self):

        self.tmpfilename = "tmp_%i" % id(self)
        infile = gzip.open(self.filename, "rb")
        outfile = open(self.tmpfilename, "wb")
        outfile.write(infile.read())
        outfile.close()
        infile.close()

    def testCompression(self):
        '''see also issue 106'''
        pysam.tabix_compress(self.tmpfilename, self.tmpfilename + ".gz")
        checkBinaryEqual(self.tmpfilename, self.tmpfilename + ".gz")

    def testIndexPresetUncompressed(self):
        '''test indexing via preset.'''

        pysam.tabix_index(self.tmpfilename, preset=self.preset)
        # check if uncompressed file has been removed
        self.assertEqual(os.path.exists(self.tmpfilename), False)
        checkBinaryEqual(self.tmpfilename + ".gz", self.filename)
        checkBinaryEqual(self.tmpfilename + ".gz.tbi", self.filename_idx)

    def testIndexPresetCompressed(self):
        '''test indexing via preset.'''

        pysam.tabix_compress(self.tmpfilename, self.tmpfilename + ".gz")
        pysam.tabix_index(self.tmpfilename + ".gz", preset=self.preset)
        checkBinaryEqual(self.tmpfilename + ".gz", self.filename)
        checkBinaryEqual(self.tmpfilename + ".gz.tbi", self.filename_idx)

    def tearDown(self):

        try:
            os.unlink(self.tmpfilename)
            os.unlink(self.tmpfilename + ".gz")
            os.unlink(self.tmpfilename + ".gz.tbi")
        except OSError:
            pass


class TestCompressionSam(TestCompression):
    filename = os.path.join(DATADIR, "example.sam.gz")
    filename_index = os.path.join(DATADIR, "example.sam.gz.tbi")
    preset = "sam"


class TestCompressionBed(TestCompression):
    filename = os.path.join(DATADIR, "example.bed.gz")
    filename_index = os.path.join(DATADIR, "example.bed.gz.tbi")
    preset = "bed"


class TestCompressionVCF(TestCompression):
    filename = os.path.join(DATADIR, "example.vcf.gz")
    filename_index = os.path.join(DATADIR, "example.vcf.gz.tbi")
    preset = "vcf"


class IterationTest(unittest.TestCase):

    with_comments = False

    def setUp(self):

        lines = []
        inf = gzip.open(self.filename, "rb")
        for line in inf:
            line = line.decode('ascii')
            if line.startswith("#"):
                if not self.with_comments:
                    continue
            lines.append(line)
        inf.close()

        # creates index of contig, start, end, adds content without newline.
        self.compare = [
            (x[0][0], int(x[0][3]), int(x[0][4]), x[1])
            for x in [(y.split("\t"), y[:-1]) for y in lines
                      if not y.startswith("#")]]

        self.comments = [x[:-1] for x in lines if x.startswith("#")]

    def getSubset(self, contig=None, start=None, end=None):

        if contig is None:
            # all lines
            subset = [x[3] for x in self.compare]
        else:
            if start is not None and end is None:
                # until end of contig
                subset = [x[3]
                          for x in self.compare if x[0] == contig
                          and x[2] > start]
            elif start is None and end is not None:
                # from start of contig
                subset = [x[3]
                          for x in self.compare if x[0] == contig
                          and x[1] <= end]
            elif start is None and end is None:
                subset = [x[3] for x in self.compare if x[0] == contig]
            else:
                # all within interval
                subset = [x[3] for x in self.compare if x[0] == contig
                          and min(x[2], end) - max(x[1], start) > 0]

        if self.with_comments:
            subset.extend(self.comments)

        return subset

    def checkPairwise(self, result, ref):
        '''check pairwise results.
        '''
        result.sort()
        ref.sort()
        a = set(result)
        b = set(ref)

        self.assertEqual(
            len(result), len(ref),
            "unexpected number of results: "
            "result=%i, expected ref=%i, differences are %s: %s"
            % (len(result), len(ref),
               a.difference(b),
               b.difference(a)))

        for x, d in enumerate(list(zip(result, ref))):
            self.assertEqual(
                d[0], d[1],
                "unexpected results in pair %i:\n'%s', expected\n'%s'" %
                (x, d[0], d[1]))


class TestGZFile(IterationTest):

    filename = os.path.join(DATADIR, "example.gtf.gz")
    with_comments = True

    def setUp(self):

        IterationTest.setUp(self)

        self.iter = pysam.GZIterator(self.filename)

    def testAll(self):
        result = list(self.iter)
        ref = self.getSubset()
        self.checkPairwise(result, ref)


class TestIterationWithoutComments(IterationTest):

    '''test iterating with TabixFile.fetch() when
    there are no comments in the file.'''

    filename = os.path.join(DATADIR,
                            "example.gtf.gz")

    def setUp(self):
        IterationTest.setUp(self)
        self.tabix = pysam.TabixFile(self.filename)

    def testRegionStrings(self):
        """test if access with various region strings
        works"""

        self.assertEqual(218, len(list(
            self.tabix.fetch("chr1"))))
        self.assertEqual(218, len(list(
            self.tabix.fetch("chr1", 1000))))
        self.assertEqual(218, len(list(
            self.tabix.fetch("chr1", end=1000000))))
        self.assertEqual(218, len(list(
            self.tabix.fetch("chr1", 1000, 1000000))))

    def testAll(self):
        result = list(self.tabix.fetch())
        ref = self.getSubset()
        self.checkPairwise(result, ref)

    def testPerContig(self):
        for contig in ("chr1", "chr2", "chr1", "chr2"):
            result = list(self.tabix.fetch(contig))
            ref = self.getSubset(contig)
            self.checkPairwise(result, ref)

    def testPerContigToEnd(self):

        end = None
        for contig in ("chr1", "chr2", "chr1", "chr2"):
            for start in range(0, 200000, 1000):
                result = list(self.tabix.fetch(contig, start, end))
                ref = self.getSubset(contig, start, end)
                self.checkPairwise(result, ref)

    def testPerContigFromStart(self):

        start = None
        for contig in ("chr1", "chr2", "chr1", "chr2"):
            for end in range(0, 200000, 1000):
                result = list(self.tabix.fetch(contig, start, end))
                ref = self.getSubset(contig, start, end)
                self.checkPairwise(result, ref)

    def testPerContig2(self):

        start, end = None, None
        for contig in ("chr1", "chr2", "chr1", "chr2"):
            result = list(self.tabix.fetch(contig, start, end))
            ref = self.getSubset(contig, start, end)
            self.checkPairwise(result, ref)

    def testPerInterval(self):

        start, end = None, None
        for contig in ("chr1", "chr2", "chr1", "chr2"):
            for start in range(0, 200000, 2000):
                for end in range(start, start + 2000, 500):
                    result = list(self.tabix.fetch(contig, start, end))
                    ref = self.getSubset(contig, start, end)
                    self.checkPairwise(result, ref)

    def testInvalidIntervals(self):

        # invalid intervals (start > end)
        self.assertRaises(ValueError, self.tabix.fetch, "chr1", 0, -10)
        self.assertRaises(ValueError, self.tabix.fetch, "chr1", 200, 0)

        # out of range intervals
        self.assertRaises(ValueError, self.tabix.fetch, "chr1", -10, 200)
        self.assertRaises(ValueError, self.tabix.fetch, "chr1", -10, -20)

        # unknown chromosome
        self.assertRaises(ValueError, self.tabix.fetch, "chrUn")

        # out of range access
        # to be implemented
        # self.assertRaises(IndexError, self.tabix.fetch, "chr1", 1000000, 2000000)

        # raise no error for invalid intervals
        self.tabix.fetch("chr1", 100, 100)

    def testGetContigs(self):
        self.assertEqual(sorted(self.tabix.contigs), [b"chr1", b"chr2"])
        # check that contigs is read-only
        self.assertRaises(
            AttributeError, setattr, self.tabix, "contigs", ["chr1", "chr2"])

    def testHeader(self):
        ref = []
        inf = gzip.open(self.filename)
        for x in inf:
            x = x.decode("ascii")
            if not x.startswith("#"):
                break
            ref.append(x[:-1].encode('ascii'))
        inf.close()

        header = list(self.tabix.header)
        self.assertEqual(ref, header)

    def testReopening(self):
        '''test repeated opening of the same file.'''
        def func1():
            # opens any tabix file
            inf = pysam.TabixFile(self.filename)
            return

        for i in range(10000):
            func1()


class TestIterationWithComments(TestIterationWithoutComments):

    '''test iterating with TabixFile.fetch() when
    there are comments in the file.

    Tests will create plenty of warnings on stderr.
    '''

    filename = os.path.join(DATADIR, "example_comments.gtf.gz")

    def setUp(self):
        TestIterationWithoutComments.setUp(self)


class TestParser(unittest.TestCase):

    filename = os.path.join(DATADIR, "example.gtf.gz")

    def setUp(self):

        self.tabix = pysam.TabixFile(self.filename)
        self.compare = loadAndConvert(self.filename)

    def testRead(self):

        for x, r in enumerate(self.tabix.fetch(parser=pysam.asTuple())):
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
        infile = gzip.open(self.filename, "rb")
        outfile = open(tmpfilename, "wb")
        outfile.write(infile.read())
        outfile.close()
        infile.close()

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


class TestIterators(unittest.TestCase):

    filename = os.path.join(DATADIR, "example.gtf.gz")

    iterator = pysam.tabix_generic_iterator
    parser = pysam.asTuple
    is_compressed = False

    def setUp(self):

        self.tabix = pysam.TabixFile(self.filename)
        self.compare = loadAndConvert(self.filename)
        self.tmpfilename_uncompressed = 'tmp_TestIterators'
        infile = gzip.open(self.filename, "rb")
        outfile = open(self.tmpfilename_uncompressed, "wb")
        outfile.write(infile.read())
        outfile.close()
        infile.close()

    def open(self):

        if self.is_compressed:
            infile = gzip.open(self.filename)
        else:
            infile = open(self.tmpfilename_uncompressed)
        return infile

    def testIteration(self):

        infile = self.open()

        for x, r in enumerate(self.iterator(infile, self.parser())):
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

    def testClosedFile(self):
        '''test for error when iterating from closed file.'''
        infile = self.open()
        infile.close()

        # iterating from a closed file should raise a value error
        self.assertRaises(ValueError, self.iterator, infile, self.parser())

    def testClosedFileIteration(self):
        '''test for error when iterating from file that has been closed'''

        infile = self.open()

        i = self.iterator(infile, self.parser())
        x = i.next()
        infile.close()
        # Not implemented
        # self.assertRaises(ValueError, i.next)

    def tearUp(self):
        os.unlink(self.tmpfilename_uncompressed)


class TestIteratorsGenericCompressed(TestIterators):
    is_compressed = True


class TestIteratorsFileCompressed(TestIterators):
    iterator = pysam.tabix_file_iterator
    is_compressed = True


class TestIteratorsFileUncompressed(TestIterators):
    iterator = pysam.tabix_file_iterator
    is_compressed = False


class TestGTF(TestParser):

    def testRead(self):

        for x, r in enumerate(self.tabix.fetch(parser=pysam.asGTF())):
            c = self.compare[x]
            self.assertEqual(len(c), len(r))
            self.assertEqual(list(c), list(r))
            self.assertEqual(c, str(r).split("\t"))
            self.assertTrue(r.gene_id.startswith("ENSG"))
            if r.feature != 'gene':
                self.assertTrue(r.transcript_id.startswith("ENST"))
            self.assertEqual(c[0], r.contig)


class TestIterationMalformattedGTFFiles(unittest.TestCase):

    '''test reading from malformatted gtf files.'''

    parser = pysam.asGTF
    iterator = pysam.tabix_generic_iterator
    parser = pysam.asGTF

    def testGTFTooManyFields(self):

        iterator = self.iterator(
            gzip.open(os.path.join(DATADIR,
                                   "gtf_toomany_fields.gtf.gz")),
            parser=self.parser())
        self.assertRaises(ValueError, iterator.next)

    def testGTFTooFewFields(self):

        iterator = self.iterator(
            gzip.open(os.path.join(DATADIR,
                                   "gtf_toofew_fields.gtf.gz")),
            parser=self.parser())
        self.assertRaises(ValueError, iterator.next)


class TestBed(unittest.TestCase):
    filename = os.path.join(DATADIR, "example.bed.gz")

    def setUp(self):

        self.tabix = pysam.TabixFile(self.filename)
        self.compare = loadAndConvert(self.filename)

    def testRead(self):

        for x, r in enumerate(self.tabix.fetch(parser=pysam.asBed())):
            c = self.compare[x]
            self.assertEqual(len(c), len(r))
            self.assertEqual(c, str(r).split("\t"))
            self.assertEqual(c[0], r.contig)
            self.assertEqual(int(c[1]), r.start)
            self.assertEqual(int(c[2]), r.end)
            self.assertEqual(list(c), list(r))

    def testWrite(self):

        for x, r in enumerate(self.tabix.fetch(parser=pysam.asBed())):
            c = self.compare[x]
            self.assertEqual(c, str(r).split("\t"))
            self.assertEqual(list(c), list(r))

            r.contig = "test"
            self.assertEqual("test", r.contig)
            self.assertEqual("test", r[0])

            r.start += 1
            self.assertEqual(int(c[1]) + 1, r.start)
            self.assertEqual(str(int(c[1]) + 1), r[1])

            r.end += 1
            self.assertEqual(int(c[2]) + 1, r.end)
            self.assertEqual(str(int(c[2]) + 1), r[2])


class TestVCF(unittest.TestCase):

    filename = os.path.join(DATADIR, "example.vcf40")

    def setUp(self):
        self.tmpfilename = "tmp_%s.vcf" % id(self)
        shutil.copyfile(self.filename, self.tmpfilename)
        pysam.tabix_index(self.tmpfilename, preset="vcf")

    def tearDown(self):
        os.unlink(self.tmpfilename + ".gz")
        if os.path.exists(self.tmpfilename + ".gz.tbi"):
            os.unlink(self.tmpfilename + ".gz.tbi")


if IS_PYTHON3:
    class TestUnicode(unittest.TestCase):

        '''test reading from a file with non-ascii characters.'''

        filename = os.path.join(DATADIR, "example_unicode.vcf")

        def setUp(self):
            self.tmpfilename = "tmp_%s.vcf" % id(self)
            shutil.copyfile(self.filename, self.tmpfilename)
            pysam.tabix_index(self.tmpfilename, preset="vcf")

        def testFromTabix(self):

            # use ascii encoding - should raise error
            t = pysam.TabixFile(self.tmpfilename + ".gz", encoding="ascii")
            results = list(t.fetch(parser=pysam.asVCF()))
            self.assertRaises(UnicodeDecodeError, getattr, results[1], "id")

            t = pysam.TabixFile(self.tmpfilename + ".gz", encoding="utf-8")
            results = list(t.fetch(parser=pysam.asVCF()))
            self.assertEqual(getattr(results[1], "id"), u"Rene\xe9")

        def testFromVCF(self):
            self.vcf = pysam.VCF()
            self.assertRaises(
                UnicodeDecodeError,
                self.vcf.connect, self.tmpfilename + ".gz", "ascii")
            self.vcf.connect(self.tmpfilename + ".gz", encoding="utf-8")
            v = self.vcf.getsamples()[0]


class TestVCFFromTabix(TestVCF):

    columns = ("contig", "pos", "id",
               "ref", "alt", "qual",
               "filter", "info", "format")

    def setUp(self):

        TestVCF.setUp(self)

        self.tabix = pysam.TabixFile(self.tmpfilename + ".gz")
        self.compare = loadAndConvert(self.filename)

    def testRead(self):

        ncolumns = len(self.columns)

        for x, r in enumerate(self.tabix.fetch(parser=pysam.asVCF())):
            c = self.compare[x]
            for y, field in enumerate(self.columns):
                # it is ok to have a missing format column
                if y == 8 and y == len(c):
                    continue
                if field == "pos":
                    self.assertEqual(int(c[y]) - 1, getattr(r, field))
                    self.assertEqual(int(c[y]) - 1, r.pos)
                else:
                    self.assertEqual(c[y], getattr(r, field),
                                     "mismatch in field %s: %s != %s" %
                                     (field, c[y], getattr(r, field)))
            if len(c) == 8:
                self.assertEqual(0, len(r))
            else:
                self.assertEqual(len(c), len(r) + ncolumns)

            for y in range(len(c) - ncolumns):
                self.assertEqual(c[ncolumns + y], r[y])

    def testWrite(self):

        ncolumns = len(self.columns)

        for x, r in enumerate(self.tabix.fetch(parser=pysam.asVCF())):
            c = self.compare[x]
            # check unmodified string
            cmp_string = str(r)
            ref_string = "\t".join([x for x in c])

            self.assertEqual(ref_string, cmp_string)

            # set fields and compare field-wise
            for y, field in enumerate(self.columns):
                # it is ok to have a missing format column
                if y == 8 and y == len(c):
                    continue
                if field == "pos":
                    rpos = getattr(r, field)
                    self.assertEqual(int(c[y]) - 1, rpos)
                    self.assertEqual(int(c[y]) - 1, r.pos)
                    # increment pos by 1
                    setattr(r, field, rpos + 1)
                    self.assertEqual(getattr(r, field), rpos + 1)
                    c[y] = str(int(c[y]) + 1)
                else:
                    setattr(r, field, "test_%i" % y)
                    c[y] = "test_%i" % y
                    self.assertEqual(c[y], getattr(r, field),
                                     "mismatch in field %s: %s != %s" %
                                     (field, c[y], getattr(r, field)))

            if len(c) == 8:
                self.assertEqual(0, len(r))
            else:
                self.assertEqual(len(c), len(r) + ncolumns)

            for y in range(len(c) - ncolumns):
                c[ncolumns + y] = "test_%i" % y
                r[y] = "test_%i" % y
                self.assertEqual(c[ncolumns + y], r[y])


class TestVCFFromVCF(TestVCF):

    columns = ("chrom", "pos", "id",
               "ref", "alt", "qual",
               "filter", "info", "format")

    # tests failing while parsing
    fail_on_parsing = (
        (5, "Flag fields should not have a value"),
        (9, "aouao"),
        (13, "aoeu"),
        (18, "Error BAD_NUMBER_OF_PARAMETERS"),
        (24, "Error HEADING_NOT_SEPARATED_BY_TABS"))

    # tests failing on opening
    fail_on_opening = ((24, "Error HEADING_NOT_SEPARATED_BY_TABS"),
                       )

    def setUp(self):

        TestVCF.setUp(self)

        self.vcf = pysam.VCF()
        self.compare = loadAndConvert(self.filename, encode=False)

    def testConnecting(self):

        fn = os.path.basename(self.filename)
        for x, msg in self.fail_on_opening:
            if "%i.vcf" % x == fn:
                self.assertRaises(ValueError,
                                  self.vcf.connect,
                                  self.tmpfilename + ".gz")
            else:
                self.vcf.connect(self.tmpfilename + ".gz")

    def testParsing(self):

        fn = os.path.basename(self.filename)
        with open(self.filename) as f:

            for x, msg in self.fail_on_opening:
                if "%i.vcf" % x == fn:
                    self.assertRaises(ValueError, self.vcf.parse, f)
                    return
            else:
                iter = self.vcf.parse(f)

            for x, msg in self.fail_on_parsing:
                if "%i.vcf" % x == fn:
                    self.assertRaises(ValueError, list, iter)
                    break
                    # python 2.7
                    # self.assertRaisesRegexp(
                    # ValueError, re.compile(msg), self.vcf.parse, f)
            else:
                # do the actual parsing
                for x, r in enumerate(iter):
                    c = self.compare[x]
                    for y, field in enumerate(self.columns):
                        # it is ok to have a missing format column
                        if y == 8 and y == len(c):
                            continue

                        val = r[field]
                        if field == "pos":
                            self.assertEqual(int(c[y]) - 1, val)
                        elif field == "alt":
                            if c[y] == ".":
                                # convert . to empty list
                                self.assertEqual(
                                    [], val,
                                    "mismatch in field %s: expected %s, got %s" %
                                    (field, c[y], val))
                            else:
                                # convert to list
                                self.assertEqual(
                                    c[y].split(","), val,
                                    "mismatch in field %s: expected %s, got %s" %
                                    (field, c[y], val))

                        elif field == "filter":
                            if c[y] == "PASS" or c[y] == ".":
                                # convert PASS to empty list
                                self.assertEqual(
                                    [], val,
                                    "mismatch in field %s: expected %s, got %s" %
                                    (field, c[y], val))
                            else:
                                # convert to list
                                self.assertEqual(
                                    c[y].split(";"), val,
                                    "mismatch in field %s: expected %s, got %s" %
                                    (field, c[y], val))

                        elif field == "info":
                            # tests for info field not implemented
                            pass
                        elif field == "qual" and c[y] == ".":
                            self.assertEqual(
                                -1, val,
                                "mismatch in field %s: expected %s, got %s" %
                                (field, c[y], val))

                        elif field == "format":
                            # format field converted to list
                            self.assertEqual(
                                c[y].split(":"), val,
                                "mismatch in field %s: expected %s, got %s" %
                                (field, c[y], val))

                        elif type(val) in (int, float):
                            if c[y] == ".":
                                self.assertEqual(
                                    None, val,
                                    "mismatch in field %s: expected %s, got %s" %
                                    (field, c[y], val))
                            else:
                                self.assertEqual(
                                    float(c[y]), float(val),
                                    "mismatch in field %s: expected %s, got %s" %
                                    (field, c[y], val))
                        else:
                            self.assertEqual(
                                c[y], val,
                                "mismatch in field %s: expected %s, got %s" %
                                (field, c[y], val))

############################################################################
# create a test class for each example vcf file.
# Two samples are created -
# 1. Testing pysam/tabix access
# 2. Testing the VCF class
vcf_files = glob.glob(os.path.join(DATADIR, "vcf", "*.vcf"))

for vcf_file in vcf_files:
    n = "VCFFromTabixTest_%s" % os.path.basename(vcf_file[:-4])
    globals()[n] = type(n, (TestVCFFromTabix,), dict(filename=vcf_file,))
    n = "VCFFromVCFTest_%s" % os.path.basename(vcf_file[:-4])
    globals()[n] = type(n, (TestVCFFromVCF,), dict(filename=vcf_file,))

############################################################################


class TestRemoteFileHTTP(unittest.TestCase):

    url = "http://genserv.anat.ox.ac.uk/downloads/pysam/test/example_htslib.gtf.gz"
    region = "chr1:1-1000"
    local = os.path.join(DATADIR, "example.gtf.gz")

    def testFetchAll(self):
        remote_file = pysam.TabixFile(self.url, "r")
        remote_result = list(remote_file.fetch())
        local_file = pysam.TabixFile(self.local, "r")
        local_result = list(local_file.fetch())

        self.assertEqual(len(remote_result), len(local_result))
        for x, y in zip(remote_result, local_result):
            self.assertEqual(x, y)


class TestIndexArgument(unittest.TestCase):

    filename_src = os.path.join(DATADIR, "example.vcf.gz")
    filename_dst = "tmp_example.vcf.gz"
    index_src = os.path.join(DATADIR, "example.vcf.gz.tbi")
    index_dst = "tmp_index_example.vcf.gz.tbi"
    preset = "vcf"

    def testFetchAll(self):
        shutil.copyfile(self.filename_src, self.filename_dst)
        shutil.copyfile(self.index_src, self.index_dst)
        same_basename_file = pysam.TabixFile(
            self.filename_src, "r", index=self.index_src)
        same_basename_results = list(same_basename_file.fetch())
        diff_index_file = pysam.TabixFile(
            self.filename_dst, "r", index=self.index_dst)
        diff_index_result = list(diff_index_file.fetch())

        self.assertEqual(len(same_basename_results), len(diff_index_result))
        for x, y in zip(same_basename_results, diff_index_result):
            self.assertEqual(x, y)


def _TestMultipleIteratorsHelper(filename, multiple_iterators):
    '''open file within scope, return iterator.'''

    tabix = pysam.TabixFile(filename)
    iterator = tabix.fetch(parser=pysam.asGTF(),
                           multiple_iterators=multiple_iterators)
    tabix.close()
    return iterator


class TestMultipleIterators(unittest.TestCase):

    filename = os.path.join(DATADIR, "example.gtf.gz")

    def testJoinedIterators(self):

        # two iterators working on the same file
        tabix = pysam.TabixFile(self.filename)
        a = tabix.fetch(parser=pysam.asGTF()).next()
        b = tabix.fetch(parser=pysam.asGTF()).next()
        # the first two lines differ only by the feature field
        self.assertEqual(a.feature, "UTR")
        self.assertEqual(b.feature, "exon")
        self.assertEqual(re.sub("UTR", "", str(a)),
                         re.sub("exon", "", str(b)))

    def testDisjointIterators(self):
        # two iterators working on the same file
        tabix = pysam.TabixFile(self.filename)
        a = tabix.fetch(parser=pysam.asGTF(), multiple_iterators=True).next()
        b = tabix.fetch(parser=pysam.asGTF(), multiple_iterators=True).next()
        # both iterators are at top of file
        self.assertEqual(str(a), str(b))

    def testScope(self):
        # technically it does not really test if the scope is correct
        i = _TestMultipleIteratorsHelper(self.filename,
                                         multiple_iterators=True)
        self.assertTrue(i.next())
        i = _TestMultipleIteratorsHelper(self.filename,
                                         multiple_iterators=False)
        self.assertRaises(IOError, i.next)

    def testDoubleFetch(self):

        f = pysam.TabixFile(self.filename)

        for a, b in zip(f.fetch(multiple_iterators=True),
                        f.fetch(multiple_iterators=True)):
            self.assertEqual(str(a), str(b))


if __name__ == "__main__":
    unittest.main()
