"""Benchmarking module for AlignmentFile functionality"""
import os
import pysam
import unittest
from TestUtils import BAM_DATADIR, IS_PYTHON3, force_str, flatten_nested_list
import PileupTestUtils


class TestPileupReadSelection(unittest.TestCase):
    '''test pileup functionality.'''

    samfilename = os.path.join(BAM_DATADIR, "ex1.bam")
    fastafilename = os.path.join(BAM_DATADIR, "ex1.fa")

    def setUp(self):

        self.samfile = pysam.AlignmentFile(self.samfilename)
        self.fastafile = pysam.FastaFile(self.fastafilename)

    def tearDown(self):
        self.samfile.close()
        self.fastafile.close()

    def check_equal(self, references, iterator):

        for x, column in enumerate(iterator):
            v = references[x][:-1].split("\t")
            self.assertEqual(
                len(v), 6,
                "expected 6 values, got {}".format(v))
            (contig, pos, reference_base,
             read_bases, read_qualities, alignment_mapping_qualities) \
                = v
            self.assertEqual(int(pos) - 1, column.reference_pos)

    def test_samtools_stepper(self):
        refs = force_str(
            pysam.samtools.mpileup(
                "-f", self.fastafilename,
                self.samfilename)).splitlines(True)
        iterator = self.samfile.pileup(
            stepper="samtools",
            fastafile=self.fastafile)
        self.check_equal(refs, iterator)

    def test_all_stepper(self):
        refs = force_str(
            pysam.samtools.mpileup(
                "-f", self.fastafilename,
                "-A", "-B",
                self.samfilename)).splitlines(True)

        iterator = self.samfile.pileup(
            stepper="all",
            fastafile=self.fastafile)
        self.check_equal(refs, iterator)

    def test_ignore_overlaps(self):
        refs = force_str(
            pysam.samtools.mpileup(
                "-f", self.fastafilename,
                "-A", "-B", "-x",
                self.samfilename)).splitlines(True)

        iterator = self.samfile.pileup(
            stepper="all",
            fastafile=self.fastafile,
            ignore_overlaps=False)
        self.check_equal(refs, iterator)

    def test_samtools_stepper_mapping_quality_threshold(self):
        refs = force_str(
            pysam.samtools.mpileup(
                "-f", self.fastafilename,
                "--min-MQ", "15",
                self.samfilename)).splitlines(True)
        iterator = self.samfile.pileup(
            stepper="samtools",
            fastafile=self.fastafile,
            min_mapping_quality=15)
        self.check_equal(refs, iterator)

    def test_samtools_stepper_base_quality_threshold(self):
        refs = force_str(
            pysam.samtools.mpileup(
                "-f", self.fastafilename,
                "--min-BQ", "20",
                self.samfilename)).splitlines(True)
        iterator = self.samfile.pileup(
            stepper="samtools",
            fastafile=self.fastafile,
            min_base_quality=20)
        self.check_equal(refs, iterator)

    def test_samtools_stepper_ignore_orphans(self):
        refs = force_str(
            pysam.samtools.mpileup(
                "-f", self.fastafilename,
                "--count-orphans",
                self.samfilename)).splitlines(True)
        iterator = self.samfile.pileup(
            stepper="samtools",
            fastafile=self.fastafile,
            ignore_orphans=False)
        self.check_equal(refs, iterator)

    def test_samtools_stepper_redo_baq(self):
        refs = force_str(
            pysam.samtools.mpileup(
                "-f", self.fastafilename,
                "--redo-BAQ",
                self.samfilename)).splitlines(True)
        iterator = self.samfile.pileup(
            stepper="samtools",
            fastafile=self.fastafile,
            redo_baq=True)
        self.check_equal(refs, iterator)
        

class TestPileupReadSelectionFastafile(TestPileupReadSelection):
    '''test pileup functionality - backwards compatibility'''

    samfilename = os.path.join(BAM_DATADIR, "ex1.bam")
    fastafilename = os.path.join(BAM_DATADIR, "ex1.fa")

    def setUp(self):

        self.samfile = pysam.AlignmentFile(self.samfilename)
        self.fastafile = pysam.Fastafile(self.fastafilename)


class TestPileupQueryPosition(unittest.TestCase):

    filename = "test_query_position.bam"

    def testPileup(self):
        last = {}
        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, self.filename)) as inf:
            for col in inf.pileup():
                for r in col.pileups:
                    # print r.alignment.query_name
                    # print r.query_position, r.query_position_or_next, r.is_del
                    if r.is_del:
                        self.assertEqual(r.query_position, None)
                        self.assertEqual(r.query_position_or_next,
                                         last[r.alignment.query_name] + 1)
                    else:
                        self.assertNotEqual(r.query_position, None)
                        last[r.alignment.query_name] = r.query_position


class TestPileupObjects(unittest.TestCase):

    def setUp(self):
        self.samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex1.bam"),
                                           "rb")

    def testPileupColumn(self):
        for pcolumn1 in self.samfile.pileup(region="chr1:105-106"):
            if pcolumn1.reference_pos == 104:
                self.assertEqual(
                    pcolumn1.reference_id, 0,
                    "chromosome/target id mismatch in position 1: %s != %s" %
                    (pcolumn1.reference_id, 0))
                self.assertEqual(
                    pcolumn1.reference_name, "chr1",
                    "chromosome mismatch in position 1: %s != %s" %
                    (pcolumn1.reference_name, "chr1"))
                self.assertEqual(
                    pcolumn1.reference_pos, 105 - 1,
                    "position mismatch in position 1: %s != %s" %
                    (pcolumn1.reference_pos, 105 - 1))
                self.assertEqual(
                    pcolumn1.nsegments, 1,
                    "# reads mismatch in position 1: %s != %s" %
                    (pcolumn1.nsegments, 1))
                self.assertEqual(
                    len(pcolumn1.pileups), 1,
                    "# reads aligned to column mismatch in position 1"
                    ": %s != %s" %
                    (len(pcolumn1.pileups), 1))

        for pcolumn2 in self.samfile.pileup(region="chr2:1480-1481"):
            if pcolumn2.reference_pos == 1479:
                self.assertEqual(
                    pcolumn2.reference_id, 1,
                    "chromosome/target id mismatch in position 1: %s != %s" %
                    (pcolumn2.reference_id, 1))
                self.assertEqual(
                    pcolumn2.reference_name, "chr2",
                    "chromosome mismatch in position 1: %s != %s" %
                    (pcolumn2.reference_name, "chr2"))
                self.assertEqual(
                    pcolumn2.reference_pos, 1480 - 1,
                    "position mismatch in position 1: %s != %s" %
                    (pcolumn2.reference_pos, 1480 - 1))
                self.assertEqual(
                    pcolumn2.nsegments, 12,
                    "# reads mismatch in position 1: %s != %s" %
                    (pcolumn2.nsegments, 12))

    def tearDown(self):
        self.samfile.close()

    def testIteratorOutOfScope(self):
        '''test if exception is raised if pileup col is accessed after
        iterator is exhausted.'''

        for pileupcol in self.samfile.pileup():
            pass

        self.assertRaises(ValueError, getattr, pileupcol, "pileups")


class TestIteratorColumnBAM(unittest.TestCase):

    '''test iterator column against contents of ex4.bam.'''

    # note that samfile contains 1-based coordinates
    # 1D means deletion with respect to reference sequence
    #
    mCoverages = {'chr1': [0] * 20 + [1] * 36 + [0] * (100 - 20 - 35),
                  'chr2': [0] * 20 + [1] * 35 + [0] * (100 - 20 - 35),
                  }

    def setUp(self):
        self.samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex4.bam"),
                                           "rb")

    def checkRange(self, contig, start=None, end=None, truncate=False):
        '''compare results from iterator with those from samtools.'''
        # check if the same reads are returned and in the same order
        for column in self.samfile.pileup(
                contig, start, end, truncate=truncate, min_base_quality=0):
            if truncate:
                self.assertGreaterEqual(column.reference_pos, start)
                self.assertLess(column.reference_pos, end)
            thiscov = len(column.pileups)
            refcov = self.mCoverages[
                self.samfile.getrname(column.reference_id)][column.reference_pos]
            self.assertEqual(thiscov, refcov,
                             "wrong coverage at pos %s:%i %i should be %i" % (
                                 self.samfile.getrname(column.reference_id),
                                 column.reference_pos, thiscov, refcov))

    def testIterateAll(self):
        '''check random access per contig'''
        self.checkRange(None)

    def testIteratePerContig(self):
        '''check random access per contig'''
        for contig in self.samfile.references:
            self.checkRange(contig)

    def testIterateRanges(self):
        '''check random access per range'''
        for contig, length in zip(
                self.samfile.references, self.samfile.lengths):
            for start in range(1, length, 90):
                # this includes empty ranges
                self.checkRange(contig, start, start + 90)

    def testInverse(self):
        '''test the inverse, is point-wise pileup accurate.'''
        for contig, refseq in list(self.mCoverages.items()):
            refcolumns = sum(refseq)
            for pos, refcov in enumerate(refseq):
                columns = list(self.samfile.pileup(contig, pos, pos + 1))
                if refcov == 0:
                    # if no read, no coverage
                    self.assertEqual(
                        len(columns),
                        refcov,
                        "wrong number of pileup columns returned for position %s:%i, %i should be %i" % (
                            contig, pos,
                            len(columns), refcov))
                elif refcov == 1:
                    # one read, all columns of the read are returned
                    self.assertEqual(
                        len(columns),
                        refcolumns,
                        "pileup incomplete at position %i: got %i, expected %i " %
                        (pos, len(columns), refcolumns))

    def testIterateTruncate(self):
        '''check random access per range'''
        for contig, length in zip(self.samfile.references,
                                  self.samfile.lengths):
            for start in range(1, length, 90):
                # this includes empty ranges
                self.checkRange(contig, start, start + 90, truncate=True)

    def tearDown(self):
        self.samfile.close()


class TestIteratorColumn2(unittest.TestCase):

    '''test iterator column against contents of ex1.bam.'''

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex1.bam"),
            "rb")

    def testStart(self):
        # print self.samfile.fetch().next().reference_start
        # print self.samfile.pileup().next().reference_start
        pass

    def testTruncate(self):
        '''see issue 107.'''
        # note that ranges in regions start from 1
        p = self.samfile.pileup(region='chr1:170-172', truncate=True)
        columns = [x.reference_pos for x in p]
        self.assertEqual(len(columns), 3)
        self.assertEqual(columns, [169, 170, 171])

        p = self.samfile.pileup('chr1', 169, 172, truncate=True)
        columns = [x.reference_pos for x in p]

        self.assertEqual(len(columns), 3)
        self.assertEqual(columns, [169, 170, 171])

    def testAccessOnClosedIterator(self):
        '''see issue 131

        Accessing pileup data after iterator has closed.
        '''
        pcolumn = self.samfile.pileup('chr1', 170, 180).__next__()
        self.assertRaises(ValueError, getattr, pcolumn, "pileups")

    def testStr(self):
        '''test if PileupRead can be printed.'''
        iter = self.samfile.pileup('chr1', 170, 180)
        pcolumn = iter.__next__()
        s = str(pcolumn)
        self.assertEqual(len(s.split("\n")), 2)


@unittest.skipIf(not IS_PYTHON3,
                 "tests requires at least python3 for subprocess context manager")
class PileUpColumnTests(unittest.TestCase):

    fn = os.path.join(BAM_DATADIR, "ex2.bam")
    fn_fasta = os.path.join(BAM_DATADIR, "ex1.fa")
    
    def test_pileup_depths_are_equal(self):
        samtools_result = PileupTestUtils.build_depth_with_samtoolspipe(self.fn)
        pysam_result = PileupTestUtils.build_depth_with_filter_with_pysam(self.fn)
        self.assertEqual(pysam_result, samtools_result)

    def test_pileup_query_bases_without_reference_are_equal(self):
        samtools_result = PileupTestUtils.build_query_bases_with_samtoolspipe(self.fn)
        pysam_result = PileupTestUtils.build_query_bases_with_pysam(self.fn)
        self.assertEqual(["".join(x) for x in pysam_result], samtools_result)

    def test_pileup_query_bases_with_reference_are_equal(self):
        samtools_result = PileupTestUtils.build_query_bases_with_samtoolspipe(self.fn, "-f", self.fn_fasta)
        with pysam.FastaFile(self.fn_fasta) as fasta:
            pysam_result = PileupTestUtils.build_query_bases_with_pysam(self.fn, fastafile=fasta, stepper="samtools")
        self.assertEqual(["".join(x) for x in pysam_result], samtools_result)
        
    def test_pileup_query_qualities_are_equal(self):
        samtools_result = PileupTestUtils.build_query_qualities_with_samtoolspipe(self.fn)
        pysam_result = PileupTestUtils.build_query_qualities_with_pysam(self.fn)
        pysam_result = [
            [chr(min(126, x + 33)) for x in l] for l in pysam_result]
        self.assertEqual("".join(flatten_nested_list(pysam_result)),
                         "".join(flatten_nested_list(samtools_result)))

    def test_pileup_mapping_qualities_are_equal(self):
        samtools_result = PileupTestUtils.build_mapping_qualities_with_samtoolspipe(self.fn)
        pysam_result = PileupTestUtils.build_mapping_qualities_with_pysam(self.fn)
        # convert to chars
        pysam_result = [
            [chr(min(126, x + 33)) for x in l] for l in pysam_result]

        self.assertEqual("".join(flatten_nested_list(pysam_result)),
                         "".join(flatten_nested_list(samtools_result)))

    def test_pileup_query_qualities_from_pileups_are_equal(self):
        samtools_result = PileupTestUtils.build_query_qualities_with_samtoolspipe(self.fn)
        pysam_result = PileupTestUtils.build_query_qualities_with_pysam_pileups(self.fn)
        pysam_result = [
            "".join([chr(min(126, x + 33)) for x in l]) for l in pysam_result]

        self.assertEqual(pysam_result, samtools_result)


if __name__ == "__main__":
    unittest.main()
