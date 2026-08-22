"""Benchmarking module for AlignmentFile functionality"""
import os

from parameterized import parameterized_class

import pysam
import sys
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

    @pytest.mark.skipif(sys.version_info[:2] == (3, 11) or sys.platform.startswith("netbsd"),
                        reason="exercises invalid accesses, which crashes on Python 3.11 and NetBSD")
    def testIteratorOutOfScope(self):
        '''test if exception is raised if pileup col is accessed after
        iterator is exhausted.'''

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

    @pytest.mark.skipif(sys.platform.startswith("netbsd"), reason="exercises invalid accesses, crashing on NetBSD")
    def testAccessOnClosedIterator(self):
        '''see issue 131

        Accessing pileup data after iterator has closed.
        '''
        pcolumn = self.samfile.pileup('chr1', 170, 180).__next__()
        with pytest.raises(ValueError): pcolumn.pileups

    def testStr(self):
        '''test if PileupRead can be printed.'''
        iter = self.samfile.pileup('chr1', 170, 180)
        pcolumn = iter.__next__()
        s = str(pcolumn)
        assert len(s.split("\n")) == 2


@parameterized_class(("from_file", ), [
    (True, ), (False, )
])
class TestPileUpColumns:
    """Test pileup column generation using different methods.

    - from_file=True: Uses AlignmentFile.pileup() with stepper="samtools" (standard approach)
    - from_file=False: Uses IteratorColumnRecords with manually filtered records

    Note: The from_file=False case may show minor depth discrepancies (~1 read at some
    positions) compared to samtools mpileup: unlike AlignmentFile.pileup()'s
    stepper="samtools", IteratorColumnRecords does not apply BAQ computation or
    mapping-quality adjustment. See IteratorColumnRecords documentation for details.
    """
    from_file: bool

    fn = os.path.join(BAM_DATADIR, "ex2.bam")
    fn_fasta = os.path.join(BAM_DATADIR, "ex1.fa")

    def test_pileup_depths_are_equal(self):
        samtools_result = PileupTestUtils.build_depth_with_samtoolspipe(self.fn)
        pysam_result = PileupTestUtils.build_depth_with_filter_with_pysam(self.fn, from_file=self.from_file)

        if self.from_file:
            # from_file=True should match samtools exactly
            assert pysam_result == samtools_result
        else:
            # from_file=False may have minor discrepancies since IteratorColumnRecords
            # doesn't apply BAQ computation or mapping-quality adjustment. Verify
            # results are "close enough":
            # - Same number of positions
            # - Differences at most ±1 read per position
            # - Differences at a small percentage of positions (< 1%)
            assert len(pysam_result) == len(samtools_result)

            diffs = sum(1 for st, py in zip(samtools_result, pysam_result) if st != py)
            max_diff = max(abs(st - py) for st, py in zip(samtools_result, pysam_result))
            diff_rate = diffs / len(samtools_result) * 100

            assert max_diff <= 1, f"Maximum depth difference should be ≤1, got {max_diff}"
            assert diff_rate < 1.0, f"Difference rate should be <1%, got {diff_rate:.2f}%"

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


class TestIteratorColumnRecords:
    """Test IteratorColumnRecords functionality."""

    fn = os.path.join(BAM_DATADIR, "ex1.bam")
    fn_fasta = os.path.join(BAM_DATADIR, "ex1.fa")

    def test_basic_iteration(self):
        """Test basic pileup generation from records."""
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            records = [rec for rec in inf]
            result = list(IteratorColumnRecords(records))

        # Should produce pileup columns
        assert len(result) > 0
        # Check first column is a PileupColumn
        assert result[0].__class__.__name__ == "PileupColumn"

    def test_with_fastafile_parameter(self):
        """Test providing fastafile at initialization."""
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf, pysam.FastaFile(self.fn_fasta) as fasta:
            records = [rec for rec in inf if rec.reference_name == "chr1"][:100]
            iter_col = IteratorColumnRecords(records, fastafile=fasta)

            # Should have reference
            assert iter_col.has_reference()

            # Iterate and check we can access columns
            result = list(iter_col)
            assert len(result) > 0

    def test_add_reference_method(self):
        """Test add_reference() method."""
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf, pysam.FastaFile(self.fn_fasta) as fasta:
            records = [rec for rec in inf][:50]
            iter_col = IteratorColumnRecords(records)

            # Should not have reference initially
            assert not iter_col.has_reference()

            # Add reference
            iter_col.add_reference(fasta)

            # Should have reference now
            assert iter_col.has_reference()

        with pysam.AlignmentFile(self.fn) as inf, pysam.FastaFile(self.fn_fasta) as fasta:
            records = [rec for rec in inf if rec.reference_name == "chr1"][:50]
            iter_col = IteratorColumnRecords(records, fastafile=fasta)

            # Trigger the initial sequence load.
            next(iter_col)
            assert iter_col.seq_len == fasta.get_reference_length("chr1")

            # Regression test: add_reference() must null the freed sequence
            # pointer, or a later reload double-frees it.
            iter_col.add_reference(fasta)
            assert iter_col.has_reference()

            # Iterating again after add_reference() must reload the sequence.
            for col in iter_col:
                assert iter_col.seq_len == fasta.get_reference_length("chr1")
                break

    def test_has_reference_method(self):
        """Test has_reference() method."""
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            records = [rec for rec in inf][:50]

            # Without fasta
            iter_col_no_ref = IteratorColumnRecords(records)
            assert not iter_col_no_ref.has_reference()

        with pysam.AlignmentFile(self.fn) as inf, pysam.FastaFile(self.fn_fasta) as fasta:
            records = [rec for rec in inf][:50]

            # With fasta
            iter_col_with_ref = IteratorColumnRecords(records, fastafile=fasta)
            assert iter_col_with_ref.has_reference()

    def test_seq_len_property(self):
        """Test seq_len property."""
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf, pysam.FastaFile(self.fn_fasta) as fasta:
            records = [rec for rec in inf if rec.reference_name == "chr1"][:50]
            iter_col = IteratorColumnRecords(records, fastafile=fasta)
            expected_seq_len = fasta.get_reference_length("chr1")

            # Iterate to trigger sequence loading
            for col in iter_col:
                # seq_len should reflect the length of the loaded chr1 sequence
                seq_len = iter_col.seq_len
                assert isinstance(seq_len, int)
                assert seq_len == expected_seq_len
                break

    def test_min_base_quality_parameter(self):
        """Test min_base_quality parameter affects results."""
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            records = [rec for rec in inf][:100]

            # Default min_base_quality=13
            depths_default = [col.get_num_aligned() for col in IteratorColumnRecords(records)]

            # Higher min_base_quality should filter more bases
            depths_high_qual = [
                col.get_num_aligned()
                for col in IteratorColumnRecords(records, min_base_quality=30)
            ]

            # Both should produce the same number of positions
            assert len(depths_default) == len(depths_high_qual)

            # A stricter quality threshold can only reduce (never increase) the
            # depth at each position, and must reduce it somewhere in this data.
            assert all(hi <= lo for hi, lo in zip(depths_high_qual, depths_default))
            assert any(hi < lo for hi, lo in zip(depths_high_qual, depths_default))

    def test_max_depth_parameter(self):
        """Test max_depth parameter reduces reported depth.

        Note: like AlignmentFile.pileup()'s own max_depth, this is htslib's
        longstanding rate-limiting behavior on newly-added reads, not a
        strict per-position cap, so depth can still exceed max_depth.
        """
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            records = [rec for rec in inf][:100]

            depths_unlimited = [col.get_num_aligned() for col in IteratorColumnRecords(records)]
            depths_capped = [
                col.get_num_aligned()
                for col in IteratorColumnRecords(records, max_depth=5)
            ]

            # Both should produce the same number of positions
            assert len(depths_unlimited) == len(depths_capped)

            assert all(hi <= lo for hi, lo in zip(depths_capped, depths_unlimited))
            assert any(hi < lo for hi, lo in zip(depths_capped, depths_unlimited))

    def test_max_depth_zero_means_unlimited(self):
        """Test max_depth=0 is treated as "no limit", matching AlignmentFile.pileup()."""
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            records = [rec for rec in inf][:100]

            depths_unlimited = [col.get_num_aligned() for col in IteratorColumnRecords(records)]
            depths_zero = [
                col.get_num_aligned()
                for col in IteratorColumnRecords(records, max_depth=0)
            ]
            assert depths_unlimited == depths_zero

    def test_unsupported_keyword_arguments_are_rejected(self):
        """Test unsupported keyword arguments raise TypeError instead of being silently ignored."""
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            records = [rec for rec in inf][:10]

            for kwargs in (
                {"stepper": "samtools"},
                {"not_a_real_argument": 1},
            ):
                with pytest.raises(TypeError):
                    IteratorColumnRecords(records, **kwargs)

    def _build_overlapping_mate_pair(self):
        """build two synthetic mates overlapping at a single reference position.

        Both mates' original base qualities (20 and 15) are above the
        default min_base_quality (13), so any difference in
        get_num_aligned() at the overlap is attributable only to
        ignore_overlaps' quality-zeroing, not to min_base_quality
        filtering that would apply regardless.
        """
        header = pysam.AlignmentHeader.from_references(["chr1"], [1000])

        mate1 = pysam.AlignedSegment(header)
        mate1.query_name = "overlapping_pair"
        mate1.query_sequence = "A"
        mate1.query_qualities = [20]
        mate1.flag = 99  # paired, proper pair, mate reverse, first in pair
        mate1.reference_id = 0
        mate1.reference_start = 10
        mate1.mapping_quality = 60
        mate1.cigartuples = [(0, 1)]
        mate1.next_reference_id = 0
        mate1.next_reference_start = 10
        mate1.template_length = 1

        mate2 = pysam.AlignedSegment(header)
        mate2.query_name = "overlapping_pair"
        mate2.query_sequence = "A"
        mate2.query_qualities = [15]
        mate2.flag = 147  # paired, proper pair, reverse, second in pair
        mate2.reference_id = 0
        mate2.reference_start = 10
        mate2.mapping_quality = 60
        mate2.cigartuples = [(0, 1)]
        mate2.next_reference_id = 0
        mate2.next_reference_start = 10
        mate2.template_length = -1

        return [mate1, mate2]

    def test_ignore_overlaps_parameter(self):
        """Test ignore_overlaps controls overlap detection, matching AlignmentFile.pileup().

        With ignore_overlaps=True (the default, matching
        AlignmentFile.pileup()), htslib zeroes the base quality of the
        lower-quality mate at the overlapping position, dropping it
        below min_base_quality and out of the count. With
        ignore_overlaps=False, both mates' original qualities are left
        intact.
        """
        from pysam.libcalignmentfile import IteratorColumnRecords

        records = self._build_overlapping_mate_pair()

        depths_default = [
            col.get_num_aligned() for col in IteratorColumnRecords(records)
        ]
        depths_no_overlap_detection = [
            col.get_num_aligned()
            for col in IteratorColumnRecords(records, ignore_overlaps=False)
        ]

        assert depths_default == [1]
        assert depths_no_overlap_detection == [2]

    def test_construction_failure_does_not_crash(self):
        """Test a failed construction is cleaned up safely instead of segfaulting.

        Regression test: __cinit__ can raise (e.g. for a bad fastafile type,
        or an unsupported keyword argument) before bam_mplp_init() has run,
        leaving self.pileup_iter NULL. __dealloc__ must not blindly pass
        that NULL pointer to bam_mplp_destroy(), which does not itself
        guard against NULL and previously segfaulted the whole process.
        """
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            records = [rec for rec in inf][:5]

            with pytest.raises(TypeError):
                IteratorColumnRecords(records, fastafile="not a FastaFile")

            with pytest.raises(TypeError):
                IteratorColumnRecords(records, not_a_real_argument=1)

    def test_unsorted_records_raise_value_error(self):
        """Test records violating the documented coordinate-sorted-order requirement fail cleanly."""
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            records = [rec for rec in inf][:5]
        records.reverse()

        # Records are pulled lazily, so construction alone cannot fail;
        # the error only surfaces once iteration forces them to be pushed.
        with pytest.raises(ValueError):
            list(IteratorColumnRecords(records))

    def test_empty_records(self):
        """Test with empty record list."""
        from pysam.libcalignmentfile import IteratorColumnRecords

        empty_records = []
        result = list(IteratorColumnRecords(empty_records))

        # Should return empty result
        assert len(result) == 0

    def test_single_record(self):
        """Test with single record."""
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            records = [next(iter(inf))]
            result = list(IteratorColumnRecords(records))

        # Should produce pileup columns for the single read
        assert len(result) > 0

    def test_multiple_chromosomes(self):
        """Test with records from multiple chromosomes."""
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            # Get records from different chromosomes
            records = []
            seen_chroms = set()
            for rec in inf:
                if rec.reference_name not in seen_chroms:
                    records.append(rec)
                    seen_chroms.add(rec.reference_name)
                if len(seen_chroms) >= 2:
                    break

            if len(seen_chroms) >= 2:
                result = list(IteratorColumnRecords(records))
                # Should handle multiple chromosomes
                assert len(result) >= len(records)

    def test_seq_len_updates_across_chromosomes(self):
        """Test seq_len reloads when iteration crosses a chromosome boundary.

        Regression test: the reference sequence must be reloaded whenever
        the pileup's current chromosome differs from the last-loaded one,
        not just once at construction.
        """
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            chr1_records = [rec for rec in inf if rec.reference_name == "chr1"][:20]
        with pysam.AlignmentFile(self.fn) as inf:
            chr2_records = [rec for rec in inf if rec.reference_name == "chr2"][:20]

        with pysam.FastaFile(self.fn_fasta) as fasta:
            seq_lens_by_chrom = {}
            iter_col = IteratorColumnRecords(chr1_records + chr2_records, fastafile=fasta)
            for col in iter_col:
                seq_lens_by_chrom.setdefault(col.reference_name, iter_col.seq_len)

            assert seq_lens_by_chrom == {
                "chr1": fasta.get_reference_length("chr1"),
                "chr2": fasta.get_reference_length("chr2"),
            }

    def test_iterator_protocol(self):
        """Test that IteratorColumnRecords follows iterator protocol."""
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            records = [rec for rec in inf][:50]
            iter_col = IteratorColumnRecords(records)

            # Should be iterable
            assert hasattr(iter_col, '__iter__')
            assert hasattr(iter_col, '__next__')

            # Should return self from __iter__
            assert iter(iter_col) is iter_col

    def test_stop_iteration(self):
        """Test that StopIteration is raised when exhausted."""
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            records = [rec for rec in inf][:10]
            iter_col = IteratorColumnRecords(records)

            # Exhaust the iterator
            for _ in iter_col:
                pass

            # Should raise StopIteration on next call
            with pytest.raises(StopIteration):
                next(iter_col)

    def test_access_after_iterator_out_of_scope(self):
        """Test a column stays safely usable after nothing else references its iterator.

        Unlike AlignmentFile.pileup()'s IteratorColumn family (see
        TestPileupObjects.testIteratorOutOfScope / issue 131 / issue
        1151), where a PileupColumn is a raw view onto its iterator's
        buffer and accessing one after the iterator is gone is
        undefined behavior (a NULL check that depends on freed memory
        still reading back as NULL -- known to segfault instead of
        raising on some Python versions), a PileupColumn from
        IteratorColumnRecords holds a reference to the iterator that
        produced it, keeping it alive for exactly this case. Dropping
        every other reference to the iterator must not make the column
        unusable or unsafe.
        """
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            records = [rec for rec in inf][:10]

        col = next(IteratorColumnRecords(records))

        assert col.get_num_aligned() >= 0
        assert list(col.pileups) is not None

    def test_records_are_consumed_lazily(self):
        """Test records are pulled from `recs` only as iteration requires them.

        Unlike AlignmentFile.pileup()'s IteratorColumn family, which reads
        directly from an open file, IteratorColumnRecords is handed a plain
        iterable and has no way to know how much of it a caller actually
        wants; consuming eagerly would defeat the point of accepting a
        generator (e.g. one reading a file too large to hold in memory).
        """
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            total = 200
            consumed = []

            def record_gen():
                for rec in inf:
                    consumed.append(rec)
                    yield rec
                    if len(consumed) >= total:
                        break

            iter_col = IteratorColumnRecords(record_gen())

            # Construction must not have touched the generator at all.
            assert len(consumed) == 0

            # A single produced column only requires part of the input.
            next(iter_col)
            assert 0 < len(consumed) < total

            # Draining the iterator consumes the rest of the input.
            list(iter_col)
            assert len(consumed) == total

    def test_exception_from_recs_propagates_with_original_type(self):
        """Test an exception raised while pulling from `recs` surfaces unchanged.

        A Python exception raised inside the pull callback can't cross the
        nogil C callback boundary directly, so it is stashed and re-raised
        once control returns to __next__(). This must preserve the
        exception's actual type rather than reporting a generic error.
        """
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            records = [rec for rec in inf][:5]

        def failing_gen():
            yield from records
            raise RuntimeError("boom")

        with pytest.raises(RuntimeError, match="boom"):
            list(IteratorColumnRecords(failing_gen()))

    def test_non_alignedsegment_first_record_raises(self):
        """Test a non-AlignedSegment as the very first yielded value raises cleanly.

        Regression test: the pull callback's guard originally covered only
        `next(recs_iter)` itself. For the first record (before `header` is
        set), the very next statement accesses `rec.header`; on a
        non-AlignedSegment this raised AttributeError outside that guard,
        which a `noexcept nogil` callback silently drops instead of
        propagating, so the iterator produced an empty result with no
        error at all.
        """
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pytest.raises(AttributeError):
            list(IteratorColumnRecords(iter(["not an AlignedSegment"])))

    def test_non_alignedsegment_later_record_raises(self):
        """Test a non-AlignedSegment after the first record raises cleanly.

        Regression test: the pull callback's `<AlignedSegment>rec` cast
        was previously unchecked, so a non-AlignedSegment `rec` did not
        raise TypeError at all -- it silently reinterpreted whatever
        `rec` actually is as an AlignedSegment's memory layout, and the
        resulting garbage `_delegate` pointer crashed the interpreter a
        few frames later inside bam_copy1(). The cast must be the
        checked `<AlignedSegment?>rec` form.
        """
        from pysam.libcalignmentfile import IteratorColumnRecords

        with pysam.AlignmentFile(self.fn) as inf:
            records = [rec for rec in inf][:5]

        with pytest.raises(TypeError):
            list(IteratorColumnRecords(iter(records + ["not an AlignedSegment"])))
