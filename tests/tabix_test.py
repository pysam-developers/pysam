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
import copy
from TestUtils import checkURL

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


def loadAndConvert(filename, encode=True):
    '''load data from filename and convert all fields to string.

    Filename can be either plain or compressed (ending in .gz).
    '''
    data = []
    if filename.endswith(".gz"):
        with gzip.open(filename) as inf:
            for line in inf:
                line = line.decode("ascii")
                if line.startswith("#"):
                    continue
                d = line.strip().split("\t")
                data.append(d)
    else:
        with open(filename) as f:
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

    with open(filename1, "rb") as infile:
        d1 = infile.read()
       
    with open(filename2, "rb") as infile:
        d2 = infile.read()
 
    found = False
    for c1, c2 in zip(d1, d2):
        if c1 != c2:
            break
    else:
        found = True

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

        self.tmpfilename = "tmp_TestCompression_%i" % id(self)
        with gzip.open(self.filename, "rb") as infile, \
             open(self.tmpfilename, "wb") as outfile:
            outfile.write(infile.read())

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
        with gzip.open(self.filename, "rb") as inf:
            for line in inf:
                line = line.decode('ascii')
                if line.startswith("#"):
                    if not self.with_comments:
                        continue
                lines.append(line)

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
        self.gzfile = pysam.GZIterator(self.filename)

    def testAll(self):
        result = list(self.gzfile)
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

    def tearDown(self):
        self.tabix.close()

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

        # raise no error for empty intervals
        self.tabix.fetch("chr1", 100, 100)

    def testGetContigs(self):
        self.assertEqual(sorted(self.tabix.contigs), ["chr1", "chr2"])
        # check that contigs is read-only
        self.assertRaises(
            AttributeError, setattr, self.tabix, "contigs", ["chr1", "chr2"])

    def testHeader(self):
        ref = []
        with gzip.open(self.filename) as inf:
            for x in inf:
                x = x.decode("ascii")
                if not x.startswith("#"):
                    break
                ref.append(x[:-1].encode('ascii'))

        header = list(self.tabix.header)
        self.assertEqual(ref, header)

    def testReopening(self):
        '''test repeated opening of the same file.'''
        def func1():
            # opens any tabix file
            with pysam.TabixFile(self.filename) as inf:
                pass

        for i in range(1000):
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
            self.assertEqual("\t".join(map(str, c)),
                             str(r))

    def testSetting(self):

        for r in self.tabix.fetch(parser=pysam.asGTF()):
            r.contig = r.contig + "_test"          
            r.source = r.source + "_test"
            r.feature = r.feature + "_test"
            r.start += 10
            r.end += 10
            r.score = 20
            r.strand = "+"
            r.frame = 0
            r.attributes = 'gene_id "0001";'


class TestIterators(unittest.TestCase):

    filename = os.path.join(DATADIR, "example.gtf.gz")

    iterator = pysam.tabix_generic_iterator
    parser = pysam.asTuple
    is_compressed = False

    def setUp(self):

        self.tabix = pysam.TabixFile(self.filename)
        self.compare = loadAndConvert(self.filename)
        self.tmpfilename_uncompressed = 'tmp_TestIterators'
        with gzip.open(self.filename, "rb") as infile, \
             open(self.tmpfilename_uncompressed, "wb") as outfile:
            outfile.write(infile.read())

    def tearDown(self):
        self.tabix.close()
        os.unlink(self.tmpfilename_uncompressed)

    def open(self):

        if self.is_compressed:
            infile = gzip.open(self.filename)
        else:
            infile = open(self.tmpfilename_uncompressed)
        return infile

    def testIteration(self):

        with self.open() as infile:
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


class TestIteratorsGenericCompressed(TestIterators):
    is_compressed = True


class TestIteratorsFileCompressed(TestIterators):
    iterator = pysam.tabix_file_iterator
    is_compressed = True


class TestIteratorsFileUncompressed(TestIterators):
    iterator = pysam.tabix_file_iterator
    is_compressed = False


class TestIterationMalformattedGTFFiles(unittest.TestCase):

    '''test reading from malformatted gtf files.'''

    parser = pysam.asGTF
    iterator = pysam.tabix_generic_iterator
    parser = pysam.asGTF

    def testGTFTooManyFields(self):

        with gzip.open(os.path.join(
                DATADIR,
                "gtf_toomany_fields.gtf.gz")) as infile:
            iterator = self.iterator(
                infile,
                parser=self.parser())
            self.assertRaises(ValueError, iterator.next)

    def testGTFTooFewFields(self):

        with gzip.open(os.path.join(
                DATADIR,
                "gtf_toofew_fields.gtf.gz")) as infile:
            iterator = self.iterator(
                infile,
                parser=self.parser())
            self.assertRaises(ValueError, iterator.next)


class TestBed(unittest.TestCase):
    filename = os.path.join(DATADIR, "example.bed.gz")

    def setUp(self):

        self.tabix = pysam.TabixFile(self.filename)
        self.compare = loadAndConvert(self.filename)

    def tearDown(self):
        self.tabix.close()

    def testRead(self):

        for x, r in enumerate(self.tabix.fetch(parser=pysam.asBed())):
            c = self.compare[x]
            self.assertEqual(len(c), len(r))
            self.assertEqual(c, str(r).split("\t"))
            self.assertEqual(c[0], r.contig)
            self.assertEqual(int(c[1]), r.start)
            self.assertEqual(int(c[2]), r.end)
            self.assertEqual(list(c), list(r))
            self.assertEqual("\t".join(map(str, c)),
                             str(r))

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
            with pysam.TabixFile(
                    self.tmpfilename + ".gz", encoding="ascii") as t:
                results = list(t.fetch(parser=pysam.asVCF()))
                self.assertRaises(UnicodeDecodeError, getattr, results[1], "id")

            with pysam.TabixFile(
                    self.tmpfilename + ".gz", encoding="utf-8") as t:
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

    def tearDown(self):
        self.tabix.close()

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
            self.assertEqual("\t".join(map(str, c)),
                             str(r))

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

    fail_on_samples = []

    check_samples = False
    coordinate_offset = 1

    # value returned for missing values
    missing_value = "."
    missing_quality = -1

    def setUp(self):

        TestVCF.setUp(self)

        self.vcf = pysam.VCF()
        self.compare = loadAndConvert(self.filename, encode=False)

    def tearDown(self):
        self.vcf.close()

    def testConnecting(self):

        fn = os.path.basename(self.filename)
        for x, msg in self.fail_on_opening:
            if "%i.vcf" % x == fn:
                self.assertRaises(ValueError,
                                  self.vcf.connect,
                                  self.tmpfilename + ".gz")
            else:
                self.vcf.connect(self.tmpfilename + ".gz")

    def get_iterator(self):

        with open(self.filename) as f:
            fn = os.path.basename(self.filename)

            for x, msg in self.fail_on_opening:
                if "%i.vcf" % x == fn:
                    self.assertRaises(ValueError, self.vcf.parse, f)
                    return

            for vcf_code, msg in self.fail_on_parsing:
                if "%i.vcf" % vcf_code == fn:
                    self.assertRaises((ValueError,
                                       AssertionError),
                                      list, self.vcf.parse(f))
                    return
                # python 2.7
                # self.assertRaisesRegexp(
                # ValueError, re.compile(msg), self.vcf.parse, f)

            return list(self.vcf.parse(f))

    def get_field_value(self, record, field):
        return record[field]

    def sample2value(self, r, v):
        return r, v

    def alt2value(self, r, v):
        if r == ".":
            return [], v
        else:
            return r.split(","), list(v)

    def filter2value(self, r, v):
        if r == "PASS":
            return [], v
        elif r == ".":
            return [], v
        else:
            return r.split(";"), v

    def testParsing(self):

        itr = self.get_iterator()
        if itr is None:
            return

        fn = os.path.basename(self.filename)

        for vcf_code, msg in self.fail_on_parsing:
            if "%i.vcf" % vcf_code == fn:
                self.assertRaises((ValueError,
                                   AssertionError),
                                  list, itr)
                return
                # python 2.7
                # self.assertRaisesRegexp(
                # ValueError, re.compile(msg), self.vcf.parse, f)

        check_samples = self.check_samples
        for vcf_code, msg in self.fail_on_samples:
            if "%i.vcf" % vcf_code == fn:
                check_samples = False

        for x, r in enumerate(itr):
            c = self.compare[x]
            for y, field in enumerate(self.columns):
                # it is ok to have a missing format column
                if y == 8 and y == len(c):
                    continue

                val = self.get_field_value(r, field)
                if field == "pos":
                    self.assertEqual(int(c[y]) - self.coordinate_offset,
                                     val)
                elif field == "alt" or field == "alts":
                    cc, vv = self.alt2value(c[y], val)
                    if cc != vv:
                        # import pdb; pdb.set_trace()
                        pass
                    self.assertEqual(
                        cc, vv,
                        "mismatch in field %s: expected %s, got %s" %
                        (field, cc, vv))

                elif field == "filter":
                    cc, vv = self.filter2value(c[y], val)
                    self.assertEqual(
                        cc, vv,
                        "mismatch in field %s: expected %s, got %s" %
                        (field, cc, vv))

                elif field == "info":
                    # tests for info field not implemented
                    pass

                elif field == "qual" and c[y] == ".":
                    self.assertEqual(
                        self.missing_quality, val,
                        "mismatch in field %s: expected %s, got %s" %
                        (field, c[y], val))

                elif field == "format":
                    # format field converted to list
                    self.assertEqual(
                        c[y].split(":"), list(val),
                        "mismatch in field %s: expected %s, got %s" %
                        (field, c[y], val))

                elif type(val) in (int, float):
                    if c[y] == ".":
                        self.assertEqual(
                            None, val,
                            "mismatch in field %s: expected %s, got %s" %
                            (field, c[y], val))
                    else:
                        self.assertAlmostEqual(
                            float(c[y]), float(val), 2,
                            "mismatch in field %s: expected %s, got %s" %
                            (field, c[y], val))
                else:
                    if c[y] == ".":
                        ref_val = self.missing_value
                    else:
                        ref_val = c[y]
                    self.assertEqual(
                        ref_val, val,
                        "mismatch in field %s: expected %s(%s), got %s(%s)" %
                        (field, ref_val, type(ref_val), val, type(val)))
            # parse samples
            if check_samples:
                if len(c) == 8:
                    for x, s in enumerate(r.samples):
                        self.assertEqual(
                            [], r.samples[s].values(),
                            "mismatch in sample {}: "
                            "expected [], got {}, src={}, line={}".format(
                                s, r.samples[s].values(),
                                r.samples[s].items(), r))
                else:
                    for x, s in enumerate(r.samples):
                        ref, comp = self.sample2value(
                            c[9 + x],
                            r.samples[s])
                        self.compare_samples(ref, comp, s, r)

    def compare_samples(self, ref, comp, s, r):

        if ref != comp:

            # check if GT not at start, not VCF conform and
            # not supported by cbcf.pyx
            k = r.format.keys()
            if "GT" in k and k[0] != "GT":
                return

            # perform an element-wise checto work around rounding differences
            for a, b in zip(re.split("[:,;]", ref),
                            re.split("[:,;]", comp)):
                is_float = True
                try:
                    a = float(a)
                    b = float(b)
                except ValueError:
                    is_float = False

                if is_float:
                    self.assertAlmostEqual(
                        a, b, 2,
                        "mismatch in sample {}: "
                        "expected {}, got {}, src={}, line={}"
                        .format(
                            s, ref, comp,
                            r.samples[s].items(), r))
                else:
                    self.assertEqual(
                        a, b,
                        "mismatch in sample {}: "
                        "expected {}, got {}, src={}, line={}"
                        .format(
                            s, ref, comp,
                            r.samples[s].items(), r))


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


class TestVCFFromVariantFile(TestVCFFromVCF):

    columns = ("chrom", "pos", "id",
               "ref", "alts", "qual",
               "filter", "info", "format")

    fail_on_parsing = []
    fail_on_opening = []
    coordinate_offset = 0
    check_samples = True
    fail_on_samples = [
        (9, "PL field not defined. Expected to be scalar, but is array"),
        (12, "PL field not defined. Expected to be scalar, but is array"),
        (18, "PL field not defined. Expected to be scalar, but is array"),
    ]

    # value returned for missing values
    missing_value = None
    missing_quality = None

    vcf = None

    def filter2value(self, r, v):
        if r == "PASS":
            return ["PASS"], list(v)
        elif r == ".":
            return [], list(v)
        else:
            return r.split(";"), list(v)

    def alt2value(self, r, v):
        if r == ".":
            return None, v
        else:
            return r.split(","), list(v)

    def sample2value(self, r, smp):

        def convert_field(f):
            if f is None:
                return "."
            elif isinstance(f, tuple):
                return ",".join(map(convert_field, f))
            else:
                return str(f)

        v = smp.values()

        if 'GT' in smp:
            alleles = [str(a) if a is not None else '.' for a in smp.allele_indices]
            v[0] = '/|'[smp.phased].join(alleles)

        comp = ":".join(map(convert_field, v))

        if comp.endswith(":."):
            comp = comp[:-2]

        return r, comp

    def setUp(self):
        TestVCF.setUp(self)
        self.compare = loadAndConvert(self.filename, encode=False)

    def tearDown(self):
        if self.vcf:
            self.vcf.close()
        self.vcf = None

    def get_iterator(self):
        self.vcf = pysam.VariantFile(self.filename)
        return self.vcf.fetch()

    def get_field_value(self, record, field):
        return getattr(record, field)


for vcf_file in vcf_files:
    n = "TestVCFFromVariantFile_%s" % os.path.basename(vcf_file[:-4])
    globals()[n] = type(n, (TestVCFFromVariantFile,), dict(filename=vcf_file,))


class TestRemoteFileHTTP(unittest.TestCase):

    url = "http://genserv.anat.ox.ac.uk/downloads/pysam/test/example_htslib.gtf.gz"
    region = "chr1:1-1000"
    local = os.path.join(DATADIR, "example.gtf.gz")

    def setUp(self):
        if not checkURL(self.url):
            self.remote_file = None
            return

        self.remote_file = pysam.TabixFile(self.url, "r")
        self.local_file = pysam.TabixFile(self.local, "r")

    def tearDown(self):
        if self.remote_file is None:
            return

        self.remote_file.close()
        self.local_file.close()

    def testFetchAll(self):
        if self.remote_file is None:
            return

        remote_result = list(self.remote_file.fetch())
        local_result = list(self.local_file.fetch())

        self.assertEqual(len(remote_result), len(local_result))
        for x, y in zip(remote_result, local_result):
            self.assertEqual(x, y)

    def testHeader(self):
        if self.remote_file is None:
            return

        self.assertEqual(list(self.local_file.header), [])
        self.assertRaises(AttributeError,
                          getattr,
                          self.remote_file,
                          "header")


class TestIndexArgument(unittest.TestCase):

    filename_src = os.path.join(DATADIR, "example.vcf.gz")
    filename_dst = "tmp_example.vcf.gz"
    index_src = os.path.join(DATADIR, "example.vcf.gz.tbi")
    index_dst = "tmp_index_example.vcf.gz.tbi"
    preset = "vcf"

    def testFetchAll(self):
        shutil.copyfile(self.filename_src, self.filename_dst)
        shutil.copyfile(self.index_src, self.index_dst)

        with pysam.TabixFile(
                self.filename_src, "r", index=self.index_src) as same_basename_file:
            same_basename_results = list(same_basename_file.fetch())

        with pysam.TabixFile(
                self.filename_dst, "r", index=self.index_dst) as diff_index_file:
            diff_index_result = list(diff_index_file.fetch())

        self.assertEqual(len(same_basename_results), len(diff_index_result))
        for x, y in zip(same_basename_results, diff_index_result):
            self.assertEqual(x, y)

        os.unlink(self.filename_dst)
        os.unlink(self.index_dst)


def _TestMultipleIteratorsHelper(filename, multiple_iterators):
    '''open file within scope, return iterator.'''

    tabix = pysam.TabixFile(filename)
    iterator = tabix.fetch(parser=pysam.asGTF(),
                           multiple_iterators=multiple_iterators)
    tabix.close()
    return iterator


class TestBackwardsCompatibility(unittest.TestCase):
    """check if error is raised if a tabix file from an
    old version is accessed from pysam"""

    def check(self, filename, raises=None):
        with pysam.TabixFile(filename) as tf:
            ref = loadAndConvert(filename)
            if raises is None:
                self.assertEqual(len(list(tf.fetch())), len(ref))
            else:
                self.assertRaises(raises, tf.fetch)

    def testVCF0v23(self):
        self.check(os.path.join(DATADIR, "example_0v23.vcf.gz"),
                   ValueError)

    def testBED0v23(self):
        self.check(os.path.join(DATADIR, "example_0v23.bed.gz"),
                   ValueError)

    def testVCF0v26(self):
        self.check(os.path.join(DATADIR, "example_0v26.vcf.gz"),
                   ValueError)

    def testBED0v26(self):
        self.check(os.path.join(DATADIR, "example_0v26.bed.gz"),
                   ValueError)

    def testVCF(self):
        self.check(os.path.join(DATADIR, "example.vcf.gz"))

    def testBED(self):
        self.check(os.path.join(DATADIR, "example.bed.gz"))

    def testEmpty(self):
        self.check(os.path.join(DATADIR, "empty.bed.gz"))


class TestMultipleIterators(unittest.TestCase):

    filename = os.path.join(DATADIR, "example.gtf.gz")

    def testJoinedIterators(self):

        # two iterators working on the same file
        with pysam.TabixFile(self.filename) as tabix:
            a = tabix.fetch(parser=pysam.asGTF()).next()
            b = tabix.fetch(parser=pysam.asGTF()).next()
            # the first two lines differ only by the feature field
            self.assertEqual(a.feature, "UTR")
            self.assertEqual(b.feature, "exon")
            self.assertEqual(re.sub("UTR", "", str(a)),
                             re.sub("exon", "", str(b)))

    def testDisjointIterators(self):
        # two iterators working on the same file
        with pysam.TabixFile(self.filename) as tabix:
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

        with pysam.TabixFile(self.filename) as f:

            for a, b in zip(f.fetch(multiple_iterators=True),
                            f.fetch(multiple_iterators=True)):
                self.assertEqual(str(a), str(b))


class TestContextManager(unittest.TestCase):

    filename = os.path.join(DATADIR, "example.gtf.gz")

    def testManager(self):

        with pysam.TabixFile(self.filename) as tabixfile:
            tabixfile.fetch()
        self.assertEqual(tabixfile.closed, True)


if __name__ == "__main__":
    unittest.main()
