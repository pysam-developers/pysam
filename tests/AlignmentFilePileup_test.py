"""Benchmarking module for AlignmentFile functionality"""
import os
import pysam
import pytest

from TestUtils import make_data_files, BAM_DATADIR, force_str, flatten_nested_list
import PileupTestUtils


def setUpModule():
    make_data_files(BAM_DATADIR)


class TestPileupReadSelection:
    '''test pileup functionality.'''

    samfilename = os.path.join(BAM_DATADIR, "ex1.bam")
    fastafilename = os.path.join(BAM_DATADIR, "ex1.fa")

    def setup_method(self):
        self.samfile = pysam.AlignmentFile(self.samfilename)
        self.fastafile = pysam.FastaFile(self.fastafilename)

    def teardown_method(self):
        self.samfile.close()
        self.fastafile.close()

    def check_equal(self, references, iterator):

        for x, column in enumerate(iterator):
            v = references[x][:-1].split("\t")
            assert len(v) == 6
            (contig, pos, reference_base,
             read_bases, read_qualities, alignment_mapping_qualities) \
                = v
            assert int(pos) - 1 == column.reference_pos

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

    def setup_method(self):
        self.samfile = pysam.AlignmentFile(self.samfilename)
        self.fastafile = pysam.Fastafile(self.fastafilename)


class TestPileupQueryPosition:

    filename = "test_query_position.bam"

    def testPileup(self):
        last = {}
        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, self.filename)) as inf:
            for col in inf.pileup():
                for r in col.pileups:
                    # print r.alignment.query_name
                    # print r.query_position, r.query_position_or_next, r.is_del
                    if r.is_del:
                        assert r.query_position is None
                        assert r.query_position_or_next == last[r.alignment.query_name] + 1
                    else:
                        assert r.query_position is not None
                        last[r.alignment.query_name] = r.query_position


class TestPileupObjects:
    def setup_method(self):
        self.samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex1.bam"),
                                           "rb")

    def testPileupColumn(self):
        for pcolumn1 in self.samfile.pileup(region="chr1:105-106"):
            if pcolumn1.reference_pos == 104:
                assert pcolumn1.reference_id == 0, "chromosome/target id mismatch in position 1"
                assert pcolumn1.reference_name == "chr1", "chromosome mismatch in position 1"
                assert pcolumn1.reference_pos == 105 - 1, "position mismatch in position 1"
                assert pcolumn1.nsegments == 1, "# reads mismatch in position 1"
                assert len(pcolumn1.pileups) == 1, "# reads aligned to column mismatch in position 1"

        for pcolumn2 in self.samfile.pileup(region="chr2:1480-1481"):
            if pcolumn2.reference_pos == 1479:
                assert pcolumn2.reference_id == 1, "chromosome/target id mismatch in position 1"
                assert pcolumn2.reference_name == "chr2", "chromosome mismatch in position 1"
                assert pcolumn2.reference_pos == 1480 - 1, "position mismatch in position 1"
                assert pcolumn2.nsegments == 12, "# reads mismatch in position 1"

    def teardown_method(self):
        self.samfile.close()

    def testIteratorOutOfScope(self):
        '''test if exception is raised if pileup col is accessed after
        iterator is exhausted.

        see issue 1151: used to crash instead of raising, depending on
        platform/Python version.
        '''

        max_n = 0
        for pileupcol in self.samfile.pileup():
            if max_n < pileupcol.n:
                max_col = pileupcol
                max_n = pileupcol.n

        with pytest.raises(ValueError): max_col.pileups
        with pytest.raises(ValueError): max_col.get_query_sequences()
        with pytest.raises(ValueError): max_col.get_num_aligned()
        with pytest.raises(ValueError): max_col.get_query_qualities()
        with pytest.raises(ValueError): max_col.get_mapping_qualities()
        with pytest.raises(ValueError): max_col.get_query_positions()
        with pytest.raises(ValueError): max_col.get_query_names()

    def testColumnFromExhaustedIteratorAfterAnotherIteratorRuns(self):
        '''see issue 1151: columns from one exhausted pileup() call,
        accessed only after a second, independent pileup() call has
        also run to completion.
        '''
        cols_a = list(self.samfile.pileup())
        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex1.bam"), "rb") as samfile2:
            cols_b = list(samfile2.pileup())

        for col in (cols_a[-1], cols_b[-1]):
            with pytest.raises(ValueError):
                col.get_num_aligned()

    def testColumnFromExhaustedUnindexedIterator(self, tmp_path):
        '''see issue 1151: an unindexed file uses IteratorColumnAll
        rather than IteratorColumnAllRefs, with the same gap.
        '''
        unindexed = tmp_path / "unindexed.bam"
        with pysam.AlignmentFile(self.samfile.filename) as src, \
                pysam.AlignmentFile(str(unindexed), "wb", template=src) as dst:
            for read in src:
                dst.write(read)

        with pysam.AlignmentFile(str(unindexed)) as inf:
            assert not inf.has_index()
            columns = list(inf.pileup())

        with pytest.raises(ValueError):
            columns[-1].get_num_aligned()


class TestIteratorColumnBAM:

    '''test iterator column against contents of ex4.bam.'''

    # note that samfile contains 1-based coordinates
    # 1D means deletion with respect to reference sequence
    #
    mCoverages = {'chr1': [0] * 20 + [1] * 36 + [0] * (100 - 20 - 35),
                  'chr2': [0] * 20 + [1] * 35 + [0] * (100 - 20 - 35),
                  }

    def setup_method(self):
        self.samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex4.bam"),
                                           "rb")

    def checkRange(self, contig, start=None, end=None, truncate=False):
        '''compare results from iterator with those from samtools.'''
        # check if the same reads are returned and in the same order
        for column in self.samfile.pileup(
                contig, start, end, truncate=truncate, min_base_quality=0):
            if truncate:
                assert column.reference_pos >= start
                assert column.reference_pos < end
            thiscov = len(column.pileups)
            refcov = self.mCoverages[
                self.samfile.getrname(column.reference_id)][column.reference_pos]
            assert thiscov == refcov, \
                   f"wrong coverage at {self.samfile.getrname(column.reference_id)}:{column.reference_pos}"

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
                    assert len(columns) == refcov, \
                           f"wrong number of pileup columns returned for position {contig}:{pos}"
                elif refcov == 1:
                    # one read, all columns of the read are returned
                    assert len(columns) == refcolumns, f"pileup incomplete at position {pos}"

    def testIterateTruncate(self):
        '''check random access per range'''
        for contig, length in zip(self.samfile.references,
                                  self.samfile.lengths):
            for start in range(1, length, 90):
                # this includes empty ranges
                self.checkRange(contig, start, start + 90, truncate=True)

    def teardown_method(self):
        self.samfile.close()


class TestIteratorColumn2:

    '''test iterator column against contents of ex1.bam.'''

    def setup_method(self):
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
        assert len(columns) == 3
        assert columns == [169, 170, 171]

        p = self.samfile.pileup('chr1', 169, 172, truncate=True)
        columns = [x.reference_pos for x in p]

        assert len(columns) == 3
        assert columns == [169, 170, 171]

    def testAccessOnClosedIterator(self):
        '''see issue 131 and issue 1151

        Accessing pileup data after iterator has closed.
        '''
        pcolumn = self.samfile.pileup('chr1', 170, 180).__next__()
        n = pcolumn.n
        assert len(pcolumn.pileups) == n
        assert pcolumn.get_num_aligned() == n
        assert len(pcolumn.get_query_sequences()) == n
        assert len(pcolumn.get_query_qualities()) == n
        assert len(pcolumn.get_mapping_qualities()) == n
        assert len(pcolumn.get_query_positions()) == n
        assert len(pcolumn.get_query_names()) == n

    def testStr(self):
        '''test if PileupRead can be printed.'''
        iter = self.samfile.pileup('chr1', 170, 180)
        pcolumn = iter.__next__()
        s = str(pcolumn)
        assert len(s.split("\n")) == 2


class TestPileUpColumns:

    fn = os.path.join(BAM_DATADIR, "ex2.bam")
    fn_fasta = os.path.join(BAM_DATADIR, "ex1.fa")

    def test_pileup_depths_are_equal(self):
        samtools_result = PileupTestUtils.build_depth_with_samtoolspipe(self.fn)
        pysam_result = PileupTestUtils.build_depth_with_filter_with_pysam(self.fn)
        assert pysam_result == samtools_result

    def test_pileup_query_bases_without_reference_are_equal(self):
        samtools_result = PileupTestUtils.build_query_bases_with_samtoolspipe(self.fn)
        pysam_result = PileupTestUtils.build_query_bases_with_pysam(self.fn)
        assert ["".join(x) for x in pysam_result] == samtools_result

    def test_pileup_query_bases_with_reference_are_equal(self):
        samtools_result = PileupTestUtils.build_query_bases_with_samtoolspipe(self.fn, "-f", self.fn_fasta)
        with pysam.FastaFile(self.fn_fasta) as fasta:
            pysam_result = PileupTestUtils.build_query_bases_with_pysam(self.fn, fastafile=fasta, stepper="samtools")
        assert ["".join(x) for x in pysam_result] == samtools_result

    def test_pileup_query_qualities_are_equal(self):
        samtools_result = PileupTestUtils.build_query_qualities_with_samtoolspipe(self.fn)
        pysam_result = PileupTestUtils.build_query_qualities_with_pysam(self.fn)
        pysam_result = [
            [chr(min(126, x + 33)) for x in l] for l in pysam_result]
        assert "".join(flatten_nested_list(pysam_result)) == "".join(flatten_nested_list(samtools_result))

    def test_pileup_mapping_qualities_are_equal(self):
        samtools_result = PileupTestUtils.build_mapping_qualities_with_samtoolspipe(self.fn)
        pysam_result = PileupTestUtils.build_mapping_qualities_with_pysam(self.fn)
        # convert to chars
        pysam_result = [
            [chr(min(126, x + 33)) for x in l] for l in pysam_result]

        assert "".join(flatten_nested_list(pysam_result)) == "".join(flatten_nested_list(samtools_result))

    def test_pileup_query_qualities_from_pileups_are_equal(self):
        samtools_result = PileupTestUtils.build_query_qualities_with_samtoolspipe(self.fn)
        pysam_result = PileupTestUtils.build_query_qualities_with_pysam_pileups(self.fn)
        pysam_result = [
            "".join([chr(min(126, x + 33)) for x in l]) for l in pysam_result]

        assert pysam_result == samtools_result
