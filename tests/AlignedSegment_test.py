import os
import pysam
import unittest
import json
import collections
import string
import copy
import array

from TestUtils import checkFieldEqual, BAM_DATADIR, get_temp_filename, get_temp_context, IS_PYTHON3


if IS_PYTHON3:
    maketrans = str.maketrans
else:
    maketrans = string.maketrans

class ReadTest(unittest.TestCase):

    def build_read(self):
        '''build an example read.'''

        header = pysam.AlignmentHeader.from_references(
            ["chr1", "chr2"],
            [10000000, 10000000])
        
        a = pysam.AlignedSegment(header)
        a.query_name = "read_12345"
        a.query_sequence = "ATGC" * 10
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 20
        a.mapping_quality = 20
        a.cigartuples = ((0, 10), (2, 1), (0, 9), (1, 1), (0, 20))
        a.next_reference_id = 0
        a.next_reference_start = 200
        a.template_length = 167
        a.query_qualities = pysam.qualitystring_to_array("1234") * 10
        return a


class TestAlignedSegment(ReadTest):

    '''tests to check if aligned read can be constructed
    and manipulated.
    '''

    def testEmpty(self):

        a = pysam.AlignedSegment()
        self.assertEqual(a.query_name, None)
        self.assertEqual(a.query_sequence, None)
        self.assertEqual(pysam.qualities_to_qualitystring(
            a.query_qualities), None)
        self.assertEqual(a.flag, 0)
        self.assertEqual(a.reference_id, -1)
        self.assertEqual(a.mapping_quality, 0)
        self.assertEqual(a.cigartuples, None)
        self.assertEqual(a.tags, [])
        self.assertEqual(a.next_reference_id, -1)
        self.assertEqual(a.next_reference_start, -1)
        self.assertEqual(a.template_length, 0)
        
    def testStrOfEmptyRead(self):
        a = pysam.AlignedSegment()
        s = str(a)
        self.assertEqual(
            "None\t0\t-1\t-1\t0\tNone\t-1\t-1\t0\tNone\tNone\t[]",
            s)

    def testSettingTagInEmptyRead(self):
        '''see issue 62'''
        a = pysam.AlignedSegment()
        a.tags = (("NM", 1),)
        a.query_qualities = None
        self.assertEqual(a.tags, [("NM", 1), ])

    def testCompare(self):
        '''check comparison functions.'''
        a = self.build_read()
        b = self.build_read()

        self.assertEqual(0, a.compare(b))
        self.assertEqual(0, b.compare(a))
        self.assertTrue(a == b)
        self.assertTrue(b == a)
        self.assertFalse(a != b)
        self.assertFalse(b != a)

        b.tid = 1
        self.assertFalse(a == b)
        self.assertFalse(b == a)
        self.assertTrue(a != b)
        self.assertTrue(b != a)

    def testHashing(self):
        a = self.build_read()
        b = self.build_read()
        self.assertEqual(hash(a), hash(b))
        b.tid = 1
        self.assertNotEqual(hash(a), hash(b))

    def testUpdate(self):
        '''check if updating fields affects other variable length data
        '''
        a = self.build_read()
        b = self.build_read()

        # check qname
        b.query_name = "read_123"
        checkFieldEqual(self, a, b, "query_name")
        b.query_name = "read_12345678"
        checkFieldEqual(self, a, b, "query_name")
        b.query_name = "read_12345"
        checkFieldEqual(self, a, b)

        # check cigar
        b.cigartuples = ((0, 10), )
        checkFieldEqual(self, a, b, "cigartuples")
        b.cigartuples = ((0, 10), (2, 1), (0, 10))
        checkFieldEqual(self, a, b, "cigartuples")
        b.cigartuples = ((0, 10), (2, 1), (0, 9), (1, 1), (0, 20))
        checkFieldEqual(self, a, b)

        # check seq
        b.query_sequence = "ATGC"
        checkFieldEqual(self,
                        a, b,
                        ("query_sequence", "query_qualities", "query_length"))
        b.query_sequence = "ATGC" * 3
        checkFieldEqual(self,
                        a, b,
                        ("query_sequence", "query_qualities", "query_length"))
        b.query_sequence = "ATGC" * 10
        checkFieldEqual(self, a, b, ("query_qualities",))

        # reset qual
        b = self.build_read()

        # check flags:
        for x in (
                "is_paired", "is_proper_pair",
                "is_unmapped", "mate_is_unmapped",
                "is_reverse", "mate_is_reverse",
                "is_read1", "is_read2",
                "is_secondary", "is_qcfail",
                "is_duplicate", "is_supplementary"):
            setattr(b, x, True)
            self.assertEqual(getattr(b, x), True)
            checkFieldEqual(self, a, b, ("flag", x,))
            setattr(b, x, False)
            self.assertEqual(getattr(b, x), False)
            checkFieldEqual(self, a, b)

    def testUpdate2(self):
        '''issue 135: inplace update of sequence and quality score.

        This does not work as setting the sequence will erase
        the quality scores.
        '''
        a = self.build_read()
        a.query_sequence = a.query_sequence[5:10]
        self.assertEqual(pysam.qualities_to_qualitystring(
            a.query_qualities), None)

        a = self.build_read()
        s = pysam.qualities_to_qualitystring(a.query_qualities)
        a.query_sequence = a.query_sequence[5:10]
        a.query_qualities = pysam.qualitystring_to_array(s[5:10])

        self.assertEqual(pysam.qualities_to_qualitystring(
            a.query_qualities), s[5:10])

    def testLargeRead(self):
        '''build an example read.'''

        a = pysam.AlignedSegment()
        a.query_name = "read_12345"
        a.query_sequence = "ATGC" * 200
        a.flag = 0
        a.reference_id = -1
        a.reference_start = 20
        a.mapping_quality = 20
        a.cigartuples = ((0, 4 * 200), )
        a.next_reference_id = 0
        a.next_reference_start = 200
        a.template_length = 167
        a.query_qualities = pysam.qualitystring_to_array("1234") * 200

        return a

    def testUpdateTlen(self):
        '''check if updating tlen works'''
        a = self.build_read()
        oldlen = a.template_length
        oldlen *= 2
        a.template_length = oldlen
        self.assertEqual(a.template_length, oldlen)

    def testPositions(self):
        a = self.build_read()
        self.assertEqual(a.get_reference_positions(),
                         [20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                          31, 32, 33, 34, 35, 36, 37, 38, 39,
                          40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                          50, 51, 52, 53, 54, 55, 56, 57, 58, 59])

        self.assertEqual(a.get_aligned_pairs(),
                         [(0, 20), (1, 21), (2, 22), (3, 23), (4, 24),
                          (5, 25), (6, 26), (7, 27), (8, 28), (9, 29),
                          (None, 30),
                          (10, 31), (11, 32), (12, 33), (13, 34), (14, 35),
                          (15, 36), (16, 37), (17, 38), (18, 39), (19, None),
                          (20, 40), (21, 41), (22, 42), (23, 43), (24, 44),
                          (25, 45), (26, 46), (27, 47), (28, 48), (29, 49),
                          (30, 50), (31, 51), (32, 52), (33, 53), (34, 54),
                          (35, 55), (36, 56), (37, 57), (38, 58), (39, 59)])

        self.assertEqual(
            a.get_reference_positions(),
            [x[1] for x in a.get_aligned_pairs()
             if x[0] is not None and x[1] is not None])
        # alen is the length of the aligned read in genome
        self.assertEqual(a.reference_length,
                         a.get_aligned_pairs()[-1][0] + 1)
        # aend points to one beyond last aligned base in ref
        self.assertEqual(a.get_reference_positions()[-1],
                         a.reference_end - 1)

    def testFullReferencePositions(self):
        '''see issue 26'''
        a = self.build_read()
        a.cigar = [(4, 30), (0, 20), (1, 3), (0, 47)]

        self.assertEqual(100,
                         len(a.get_reference_positions(full_length=True)))

    def testBlocks(self):
        a = self.build_read()
        self.assertEqual(a.get_blocks(),
                         [(20, 30), (31, 40), (40, 60)])

    def test_infer_query_length(self):
        '''Test infer_query_length on M|=|X|I|D|H|S cigar ops'''
        a = self.build_read()
        a.cigarstring = '40M'
        self.assertEqual(a.infer_query_length(), 40)
        a.cigarstring = '40='
        self.assertEqual(a.infer_query_length(), 40)
        a.cigarstring = '40X'
        self.assertEqual(a.infer_query_length(), 40)
        a.cigarstring = '20M5I20M'
        self.assertEqual(a.infer_query_length(), 45)
        a.cigarstring = '20M5D20M'
        self.assertEqual(a.infer_query_length(), 40)
        a.cigarstring = '5H35M'
        self.assertEqual(a.infer_query_length(), 35)
        a.cigarstring = '5S35M'
        self.assertEqual(a.infer_query_length(), 40)
        a.cigarstring = '35M5H'
        self.assertEqual(a.infer_query_length(), 35)
        a.cigarstring = '35M5S'
        self.assertEqual(a.infer_query_length(), 40)
        a.cigarstring = None
        self.assertEqual(a.infer_query_length(), None)

    def test_infer_read_length(self):
        '''Test infer_read_length on M|=|X|I|D|H|S cigar ops'''
        a = self.build_read()
        a.cigarstring = '40M'
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = '40='
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = '40X'
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = '20M5I20M'
        self.assertEqual(a.infer_read_length(), 45)
        a.cigarstring = '20M5D20M'
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = '5H35M'
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = '5S35M'
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = '35M5H'
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = '35M5S'
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = None
        self.assertEqual(a.infer_read_length(), None)

    def test_get_aligned_pairs_soft_clipping(self):
        a = self.build_read()
        a.cigartuples = ((4, 2), (0, 35), (4, 3))
        self.assertEqual(a.get_aligned_pairs(),
                         [(0, None), (1, None)] +
                         [(qpos, refpos) for (qpos, refpos) in zip(
                             range(2, 2 + 35), range(20, 20 + 35))] +
                         [(37, None), (38, None), (39, None)]
                         )
        self.assertEqual(a.get_aligned_pairs(True),
                         # [(0, None), (1, None)] +
                         [(qpos, refpos) for (qpos, refpos) in zip(
                             range(2, 2 + 35), range(20, 20 + 35))]
                         # [(37, None), (38, None), (39, None)]
                         )

    def test_get_aligned_pairs_hard_clipping(self):
        a = self.build_read()
        a.cigartuples = ((5, 2), (0, 35), (5, 3))
        self.assertEqual(a.get_aligned_pairs(),
                         # No seq, no seq pos
                         [(qpos, refpos) for (qpos, refpos) in zip(
                             range(0, 0 + 35), range(20, 20 + 35))])
        self.assertEqual(a.get_aligned_pairs(True),
                         [(qpos, refpos) for (qpos, refpos) in zip(
                             range(0, 0 + 35), range(20, 20 + 35))])

    def test_get_aligned_pairs_skip(self):
        a = self.build_read()
        a.cigarstring = "2M100D38M"
        self.assertEqual(a.get_aligned_pairs(),
                         [(0, 20), (1, 21)] +
                         [(None, refpos) for refpos in range(22, 22 + 100)] +
                         [(qpos, refpos) for (qpos, refpos) in zip(
                             range(2, 2 + 38),
                             range(20 + 2 + 100, 20 + 2 + 100 + 38))])
        self.assertEqual(a.get_aligned_pairs(True),
                         [(0, 20), (1, 21)] +
                         # [(None, refpos) for refpos in range(21, 21+100)] +
                         [(qpos, refpos) for (qpos, refpos) in zip(
                             range(2, 2 + 38),
                             range(20 + 2 + 100, 20 + 2 + 100 + 38))])

    def test_get_aligned_pairs_match_mismatch(self):
        a = self.build_read()
        a.cigartuples = ((7, 20), (8, 20))
        self.assertEqual(a.get_aligned_pairs(),
                         [(qpos, refpos) for (qpos, refpos) in zip(
                             range(0, 0 + 40), range(20, 20 + 40))])
        self.assertEqual(a.get_aligned_pairs(True),
                         [(qpos, refpos) for (qpos, refpos) in zip(
                             range(0, 0 + 40), range(20, 20 + 40))])

    def test_get_aligned_pairs_padding(self):
        a = self.build_read()
        a.cigartuples = ((7, 20), (6, 1), (8, 19))

        def inner():
            a.get_aligned_pairs()
        # padding is not bein handled right now
        self.assertRaises(NotImplementedError, inner)

    def test_get_aligned_pairs(self):
        a = self.build_read()
        a.query_sequence = "A" * 9
        a.cigarstring = "9M"
        a.set_tag("MD", "9")
        self.assertEqual(
            a.get_aligned_pairs(with_seq=True),
            [(0, 20, 'A'), (1, 21, 'A'), (2, 22, 'A'),
             (3, 23, 'A'), (4, 24, 'A'), (5, 25, 'A'),
             (6, 26, 'A'), (7, 27, 'A'), (8, 28, 'A')])

        a.set_tag("MD", "4C4")
        self.assertEqual(
            a.get_aligned_pairs(with_seq=True),
            [(0, 20, 'A'), (1, 21, 'A'), (2, 22, 'A'),
             (3, 23, 'A'), (4, 24, 'c'), (5, 25, 'A'),
             (6, 26, 'A'), (7, 27, 'A'), (8, 28, 'A')])

        a.cigarstring = "5M2D4M"
        a.set_tag("MD", "4C^TT4")
        self.assertEqual(
            a.get_aligned_pairs(with_seq=True),
            [(0, 20, 'A'), (1, 21, 'A'), (2, 22, 'A'),
             (3, 23, 'A'), (4, 24, 'c'),
             (None, 25, 'T'), (None, 26, 'T'),
             (5, 27, 'A'), (6, 28, 'A'), (7, 29, 'A'), (8, 30, 'A')]
        )

        a.cigarstring = "5M2D2I2M"
        a.set_tag("MD", "4C^TT2")
        self.assertEqual(
            a.get_aligned_pairs(with_seq=True),
            [(0, 20, 'A'), (1, 21, 'A'), (2, 22, 'A'),
             (3, 23, 'A'), (4, 24, 'c'),
             (None, 25, 'T'), (None, 26, 'T'),
             (5, None, None), (6, None, None),
             (7, 27, 'A'), (8, 28, 'A')]
        )

    def test_get_aligned_pairs_with_malformed_MD_tag(self):

        a = self.build_read()
        a.query_sequence = "A" * 9

        # out of range issue, see issue #560
        a.cigarstring = "64M2D85M2S"
        a.set_tag("MD", "64^TG86A0")
        self.assertRaises(
            AssertionError,
            a.get_aligned_pairs,
            with_seq=True)

    def test_get_aligned_pairs_skip_reference(self):
        a = self.build_read()
        a.query_sequence = "A" * 10
        a.cigarstring = "5M1N5M"
        a.set_tag("MD", "10")

        self.assertEqual(
            a.get_aligned_pairs(with_seq=True),
            [(0, 20, 'A'), (1, 21, 'A'), (2, 22, 'A'),
             (3, 23, 'A'), (4, 24, 'A'), (None, 25, None),
             (5, 26, 'A'), (6, 27, 'A'), (7, 28, 'A'),
             (8, 29, 'A'), (9, 30, 'A')])

        self.assertEqual(
            a.get_aligned_pairs(with_seq=False),
            [(0, 20), (1, 21), (2, 22),
             (3, 23), (4, 24), (None, 25),
             (5, 26), (6, 27), (7, 28),
             (8, 29), (9, 30)])

        self.assertEqual(
            a.get_aligned_pairs(matches_only=True, with_seq=False),
            [(0, 20), (1, 21),
             (2, 22), (3, 23),
             (4, 24), (5, 26),
             (6, 27), (7, 28),
             (8, 29), (9, 30)])

    def testNoSequence(self):
        '''issue 176: retrieving length without query sequence
        with soft-clipping.
        '''
        a = self.build_read()
        a.query_sequence = None
        a.cigarstring = "20M"
        self.assertEqual(a.query_alignment_length, 20)
        a.cigarstring = "20M1S"
        self.assertEqual(a.query_alignment_length, 20)
        a.cigarstring = "20M1H"
        self.assertEqual(a.query_alignment_length, 20)
        a.cigarstring = "1S20M"
        self.assertEqual(a.query_alignment_length, 20)
        a.cigarstring = "1H20M"
        self.assertEqual(a.query_alignment_length, 20)
        a.cigarstring = "1S20M1S"
        self.assertEqual(a.query_alignment_length, 20)
        a.cigarstring = "1H20M1H"
        self.assertEqual(a.query_alignment_length, 20)

    def test_query_length_is_limited(self):
        a = self.build_read()
        a.query_name = "A" * 1
        a.query_name = "A" * 251
        self.assertRaises(
            ValueError,
            setattr,
            a,
            "query_name",
            "A" * 252)

    def test_header_accessible(self):
        a = self.build_read()
        self.assertTrue(isinstance(a.header, pysam.AlignmentHeader))

    def test_bin_values_for_unmapped_reads_ignore_length(self):
        a = self.build_read()
        # use a long read
        a.cigarstring="2000000M"
        self.assertEqual(a.bin, 9)
        # changing unmapped flag changes bin because length is 0
        a.is_unmapped = True
        self.assertTrue(a.is_unmapped)
        self.assertEqual(a.bin, 4681)
        
        # unmapped read without chromosomal location
        a.reference_start = -1
        self.assertEqual(a.reference_start, -1)
        self.assertEqual(a.bin, 4680)

    def test_bin_values_for_mapped_reads_are_updated(self):
        a = self.build_read()
        a.pos = 20000
        self.assertFalse(a.is_unmapped)
        self.assertEqual(a.bin, 4682)

        # updating length updates bin
        a.cigarstring="2000000M"
        self.assertEqual(a.bin, 9)

        # updating length updates bin
        a.cigarstring="20M"
        self.assertEqual(a.bin, 4682)

        # updating length updates bin
        a.reference_start = 2000000
        self.assertEqual(a.bin, 4803)
        

class TestTidMapping(ReadTest):

    def test_reference_name_can_be_set_to_none(self):
        a = self.build_read()
        a.reference_name = None
        self.assertEqual(a.reference_name, None)
        self.assertEqual(a.reference_id, -1)

    def test_reference_name_can_be_set_to_asterisk(self):
        a = self.build_read()
        a.reference_name = "*"
        self.assertEqual(a.reference_name, None)
        self.assertEqual(a.reference_id, -1)

    def test_reference_name_can_be_set_to_chromosome(self):
        a = self.build_read()
        a.reference_name = "chr1"
        self.assertEqual(a.reference_name, "chr1")
        self.assertEqual(a.reference_id, 0)

    def test_reference_name_can_not_be_set_to_unknown_chromosome(self):
        a = self.build_read()
        self.assertRaises(ValueError,
                          setattr,
                          a,
                          "reference_name",
                          "chrX")
        
    def test_tid_can_be_set_to_missing(self):
        a = self.build_read()
        a.reference_id = -1
        self.assertEqual(a.reference_id, -1)
        self.assertEqual(a.reference_name, None)

    def test_tid_can_be_set_to_missing_without_header(self):
        a = pysam.AlignedSegment()
        a.reference_id = -1
        self.assertEqual(a.reference_id, -1)
        self.assertEqual(a.reference_name, None)

    def test_tid_can_be_set_without_header(self):
        a = pysam.AlignedSegment()
        a.reference_id = 1
        self.assertRaises(ValueError, getattr, a, "reference_name")
        
    def test_tid_can_be_set_to_chromosome(self):
        a = self.build_read()
        a.reference_id = 0
        self.assertEqual(a.reference_id, 0)
        self.assertEqual(a.reference_name, "chr1")

    def test_tid_can_not_be_set_to_unknown_chromosome(self):
        a = self.build_read()
        self.assertRaises(ValueError,
                          setattr,
                          a,
                          "reference_id",
                          2)
        
    def test_unmapped_tid_is_asterisk_in_output(self):
        a = self.build_read()
        a.reference_id = -1
        self.assertEqual(a.to_string().split("\t")[2], "*")


class TestNextTidMapping(ReadTest):

    def test_next_reference_name_can_be_set_to_none(self):
        a = self.build_read()
        a.next_reference_name = None
        self.assertEqual(a.next_reference_name, None)
        self.assertEqual(a.next_reference_id, -1)

    def test_next_reference_name_can_be_set_to_asterisk(self):
        a = self.build_read()
        a.next_reference_name = "*"
        self.assertEqual(a.next_reference_name, None)
        self.assertEqual(a.next_reference_id, -1)

    def test_next_reference_name_can_be_set_to_chromosome(self):
        a = self.build_read()
        a.next_reference_name = "chr1"
        self.assertEqual(a.next_reference_name, "chr1")
        self.assertEqual(a.next_reference_id, 0)

    def test_next_reference_name_can_not_be_set_to_unknown_chromosome(self):
        a = self.build_read()
        self.assertRaises(ValueError,
                          setattr,
                          a,
                          "next_reference_name",
                          "chrX")
        
    def test_next_tid_can_be_set_to_missing(self):
        a = self.build_read()
        a.next_reference_id = -1
        self.assertEqual(a.next_reference_id, -1)
        self.assertEqual(a.next_reference_name, None)

    def test_next_tid_can_be_set_to_equal(self):
        a = self.build_read()
        a.reference_name = "chr1"
        a.next_reference_name = "="
        self.assertEqual(a.next_reference_id, a.reference_id)
        self.assertEqual(a.next_reference_name, a.reference_name)
        self.assertEqual(a.to_string().split("\t")[6], "=")
        
    def test_next_tid_can_be_set_to_missing_without_header(self):
        a = pysam.AlignedSegment()
        a.next_reference_id = -1
        self.assertEqual(a.next_reference_id, -1)
        self.assertEqual(a.next_reference_name, None)

    def test_next_tid_can_be_set_without_header(self):
        a = pysam.AlignedSegment()
        a.next_reference_id = 1
        self.assertRaises(ValueError, getattr, a, "next_reference_name")
        
    def test_next_tid_can_be_set_to_chromosome(self):
        a = self.build_read()
        a.next_reference_id = 0
        self.assertEqual(a.next_reference_id, 0)
        self.assertEqual(a.next_reference_name, "chr1")

    def test_next_tid_can_not_be_set_to_unknown_chromosome(self):
        a = self.build_read()
        self.assertRaises(ValueError,
                          setattr,
                          a,
                          "next_reference_id",
                          2)
        
    def test_next_unmapped_tid_is_asterisk_in_output(self):
        a = self.build_read()
        a.next_reference_id = -1
        self.assertEqual(a.to_string().split("\t")[6], "*")
        
        
class TestCigar(ReadTest):

    def testCigarString(self):
        r = self.build_read()
        self.assertEqual(r.cigarstring, "10M1D9M1I20M")
        r.cigarstring = "20M10D20M"
        self.assertEqual(r.cigartuples, [(0, 20), (2, 10), (0, 20)])
        # unsetting cigar string
        r.cigarstring = None
        self.assertEqual(r.cigarstring, None)

    def testCigar(self):
        r = self.build_read()
        self.assertEqual(
            r.cigartuples, [(0, 10), (2, 1), (0, 9), (1, 1), (0, 20)])
        # unsetting cigar string
        r.cigartuples = None
        self.assertEqual(r.cigartuples, None)


class TestCigarStats(ReadTest):

    def testStats(self):
        a = self.build_read()

        a.cigarstring = None
        self.assertEqual(
            [list(x) for x in a.get_cigar_stats()],
            [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

        a.cigarstring = "10M"
        self.assertEqual(
            [list(x) for x in a.get_cigar_stats()],
            [[10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

        a.cigarstring = "10M2I2M"
        self.assertEqual(
            [list(x) for x in a.get_cigar_stats()],
            [[12, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

        for i, x in enumerate("MIDNSHP=X"):
            a.cigarstring = "2{}".format(x)
            expected = [[0] * 11, [0] * 11]
            expected[0][i] = 2
            expected[1][i] = 1
            self.assertEqual(
                [list(x) for x in a.get_cigar_stats()],
                expected)

        a.cigarstring = "10M"
        a.set_tag("NM", 5)
        self.assertEqual(
            [list(x) for x in a.get_cigar_stats()],
            [[10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5],
             [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

        a.cigarstring = None
        self.assertEqual(
            [list(x) for x in a.get_cigar_stats()],
            [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])


class TestAlignedPairs(unittest.TestCase):
    filename = os.path.join(BAM_DATADIR, "example_aligned_pairs.bam")

    def testReferenceBases(self):
        """reference bases should always be the same nucleotide
        """
        reference_bases = collections.defaultdict(list)
        with pysam.AlignmentFile(self.filename) as inf:
            for c in inf.pileup():
                for r in c.pileups:
                    for read, ref, base in r.alignment.get_aligned_pairs(
                            with_seq=True):
                        if ref is None:
                            continue
                        reference_bases[ref].append(base.upper())

        for x, y in reference_bases.items():
            self.assertEqual(len(set(y)), 1)


class TestTags(ReadTest):

    def testMissingTag(self):
        a = self.build_read()
        self.assertRaises(KeyError, a.get_tag, "XP")

    def testEmptyTag(self):
        a = self.build_read()
        self.assertRaises(KeyError, a.get_tag, "XT")

    def testSetTag(self):
        a = self.build_read()
        self.assertEqual(False, a.has_tag("NM"))
        a.set_tag("NM", 2)
        self.assertEqual(True, a.has_tag("NM"))
        self.assertEqual(a.get_tag("NM"), 2)
        a.set_tag("NM", 3)
        self.assertEqual(a.get_tag("NM"), 3)
        a.set_tag("NM", None)
        self.assertEqual(False, a.has_tag("NM"))
        # check if deleting a non-existing tag is fine
        a.set_tag("NM", None)
        a.set_tag("NM", None)

    def testArrayTags(self):
        read = self.build_read()
        supported_dtypes = "bhBHf"
        unsupported_dtypes = "lLd"

        for dtype in supported_dtypes:
            key = "F" + dtype
            read.set_tag(key, array.array(dtype, range(10)))
            ary = read.get_tag(key)

        for dtype in unsupported_dtypes:
            key = "F" + dtype
            self.assertRaises(ValueError,
                              read.set_tag,
                              key,
                              array.array(dtype, range(10)))

    def testAddTagsType(self):
        a = self.build_read()
        a.tags = None
        self.assertEqual(a.tags, [])

        a.setTag('X1', 5.0)
        a.setTag('X2', "5.0")
        a.setTag('X3', 5)

        self.assertEqual(sorted(a.tags),
                         sorted([('X1', 5.0),
                                 ('X2', "5.0"),
                                 ('X3', 5)]))

        # test setting float for int value
        a.setTag('X4', 5, value_type='d')
        self.assertEqual(sorted(a.tags),
                         sorted([('X1', 5.0),
                                 ('X2', "5.0"),
                                 ('X3', 5),
                                 ('X4', 5.0)]))

        # test setting int for float value - the
        # value will be rounded.
        a.setTag('X5', 5.2, value_type='i')
        self.assertEqual(sorted(a.tags),
                         sorted([('X1', 5.0),
                                 ('X2', "5.0"),
                                 ('X3', 5),
                                 ('X4', 5.0),
                                 ('X5', 5)]))

        # test setting invalid type code
        self.assertRaises(ValueError, a.set_tag, 'X6', 5.2, 'g')

    def testTagsUpdatingFloat(self):
        a = self.build_read()
        a.tags = [('NM', 1), ('RG', 'L1'),
                  ('PG', 'P1'), ('XT', 'U')]

        self.assertEqual(a.tags,
                         [('NM', 1), ('RG', 'L1'),
                          ('PG', 'P1'), ('XT', 'U')])
        a.tags += [('XC', 5.0)]
        self.assertEqual(a.tags,
                         [('NM', 1), ('RG', 'L1'),
                          ('PG', 'P1'), ('XT', 'U'), ('XC', 5.0)])

    def testAddTags(self):
        a = self.build_read()
        a.tags = [('NM', 1), ('RG', 'L1'),
                  ('PG', 'P1'), ('XT', 'U')]

        self.assertEqual(sorted(a.tags),
                         sorted([('NM', 1), ('RG', 'L1'),
                                 ('PG', 'P1'), ('XT', 'U')]))

        a.setTag('X1', 'C')
        self.assertEqual(sorted(a.tags),
                         sorted([('X1', 'C'), ('NM', 1), ('RG', 'L1'),
                                 ('PG', 'P1'), ('XT', 'U'), ]))
        a.setTag('X2', 5)
        self.assertEqual(sorted(a.tags),
                         sorted([('X2', 5), ('X1', 'C'),
                                 ('NM', 1), ('RG', 'L1'),
                                 ('PG', 'P1'), ('XT', 'U'), ]))
        # add with replacement
        a.setTag('X2', 10)
        self.assertEqual(sorted(a.tags),
                         sorted([('X2', 10), ('X1', 'C'),
                                 ('NM', 1), ('RG', 'L1'),
                                 ('PG', 'P1'), ('XT', 'U'), ]))

        # add without replacement
        a.setTag('X2', 5, replace=False)
        self.assertEqual(sorted(a.tags),
                         sorted([('X2', 10), ('X1', 'C'),
                                 ('X2', 5),
                                 ('NM', 1), ('RG', 'L1'),
                                 ('PG', 'P1'), ('XT', 'U'), ]))

    def testTagParsing(self):
        '''test for tag parsing

        see http://groups.google.com/group/pysam-user-group/browse_thread/thread/67ca204059ea465a
        '''
        samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex8.bam"),
            "rb")

        for entry in samfile:
            before = entry.get_tags()
            entry.set_tags(before)
            after = entry.get_tags()
            self.assertEqual(after, before)

    def testMDTagMatchOnly(self):
        a = self.build_read()

        # Substitutions only
        a.cigarstring = "21M"
        a.query_sequence = "A" * 21
        a.set_tag('MD', "5C0T0G05C0G0T5")
        self.assertEqual(
            "AAAAActgAAAAAcgtAAAAA",
            a.get_reference_sequence())

        a.cigarstring = "21M"
        a.query_sequence = "A" * 21
        a.set_tag('MD', "5CTG5CGT5")
        self.assertEqual(
            "AAAAActgAAAAAcgtAAAAA",
            a.get_reference_sequence())

        a.cigarstring = "11M"
        a.query_sequence = "A" * 11
        a.set_tag('MD', "CTG5CGT")
        self.assertEqual(
            "ctgAAAAAcgt",
            a.get_reference_sequence())

    def testMDTagInsertions(self):
        a = self.build_read()

        # insertions are silent in the reference sequence
        a.cigarstring = "5M1I5M"
        a.query_sequence = "A" * 5 + "C" + "A" * 5
        a.set_tag('MD', "10")
        self.assertEqual(
            a.get_reference_sequence(),
            "A" * 10)

        a.cigarstring = "1I10M"
        a.query_sequence = "C" * 1 + "A" * 10
        self.assertEqual(
            a.get_reference_sequence(),
            "A" * 10)

        a.cigarstring = "10M1I"
        a.query_sequence = "A" * 10 + "C" * 1
        self.assertEqual(
            a.get_reference_sequence(),
            "A" * 10)

    def testMDTagDeletions(self):
        a = self.build_read()

        a.cigarstring = "5M1D5M"
        a.query_sequence = "A" * 10
        a.set_tag('MD', "5^C5")
        self.assertEqual(
            "A" * 5 + "C" + "A" * 5,
            a.get_reference_sequence())

        a.cigarstring = "5M3D5M"
        a.query_sequence = "A" * 10
        a.set_tag('MD', "5^CCC5")
        self.assertEqual(
            "A" * 5 + "C" * 3 + "A" * 5,
            a.get_reference_sequence())

    def testMDTagRefSkipping(self):
        a = self.build_read()

        a.cigarstring = "5M1N5M"
        a.query_sequence = "A" * 10
        a.set_tag('MD', "10")
        self.assertEqual(
            "A" * 10,
            a.get_reference_sequence())

        a.cigarstring = "5M3N5M"
        a.query_sequence = "A" * 10
        a.set_tag('MD', "10")
        self.assertEqual(
            "A" * 10,
            a.get_reference_sequence())

    def testMDTagSoftClipping(self):
        a = self.build_read()

        # softclipping
        a.cigarstring = "5S5M1D5M5S"
        a.query_sequence = "G" * 5 + "A" * 10 + "G" * 5
        a.set_tag('MD', "5^C5")
        self.assertEqual(
            "A" * 5 + "C" + "A" * 5,
            a.get_reference_sequence())

        # all together
        a.cigarstring = "5S5M1D5M1I5M5S"
        a.query_sequence = "G" * 5 + "A" * 16 + "G" * 5
        a.set_tag('MD', "2C2^T10")
        self.assertEqual(
            "AAcAATAAAAAAAAAA",
            a.get_reference_sequence())

    def testMDTagComplex(self):
        a = self.build_read()

        a.cigarstring = "5S5M1I2D5M5S"
        a.query_sequence = "G" * 5 + "A" * 11 + "G" * 5
        a.set_tag('MD', "2C2^TC5")
        self.assertEqual(
            "AAcAATCAAAAA",
            a.get_reference_sequence())

        a.cigarstring = "5S5M2D1I5M5S"
        a.query_sequence = "G" * 5 + "A" * 11 + "G" * 5
        a.set_tag('MD', "2C2^TC5")
        self.assertEqual(
            "AAcAATCAAAAA",
            a.get_reference_sequence())

        # insertion in reference overlapping deletion in reference
        # read: AACCCCA---AAA
        # ref:  AA----AGGGAAA
        a.cigarstring = "2M4I1M3D3M"
        a.set_tag("MD", "3^GGG3")
        a.query_sequence = "AACCCCAAAA"
        self.assertEqual(
            "AAAGGGAAA",
            a.get_reference_sequence())

        a.cigarstring = "5M2D2I2M"
        a.set_tag("MD", "4C^TT2")
        a.query_sequence = "A" * 9
        self.assertEqual(
            "AAAAcTTAA",
            a.get_reference_sequence())

    def testArrayTags(self):

        r = self.build_read()

        def c(r, l):
            r.tags = [('ZM', l)]
            self.assertEqual(list(r.opt("ZM")), list(l))

        # signed integers
        c(r, (-1, 1))
        c(r, (-1, 100))
        c(r, (-1, 200))
        c(r, (-1, 1000))
        c(r, (-1, 30000))
        c(r, (-1, 50000))
        c(r, (1, -1))
        c(r, (1, -100))
        c(r, (1, -200))
        c(r, (1, -1000))
        c(r, (1, -30000))
        c(r, (1, -50000))

        # unsigned integers
        c(r, (1, 100))
        c(r, (1, 1000))
        c(r, (1, 10000))
        c(r, (1, 100000))

        # floats
        c(r, (1.0, 100.0))

    def testLongTags(self):
        '''see issue 115'''

        r = self.build_read()
        rg = 'HS2000-899_199.L3'
        tags = [('XC', 85), ('XT', 'M'), ('NM', 5),
                ('SM', 29), ('AM', 29), ('XM', 1),
                ('XO', 1), ('XG', 4), ('MD', '37^ACCC29T18'),
                ('XA', '5,+11707,36M1I48M,2;21,-48119779,46M1I38M,2;hs37d5,-10060835,40M1D45M,3;5,+11508,36M1I48M,3;hs37d5,+6743812,36M1I48M,3;19,-59118894,46M1I38M,3;4,-191044002,6M1I78M,3;')]  # noqa

        r.tags = tags
        r.tags += [("RG", rg)] * 100
        tags += [("RG", rg)] * 100

        self.assertEqual(tags, r.tags)

    def testNegativeIntegers(self):
        x = -2
        aligned_read = self.build_read()
        aligned_read.tags = [("XD", int(x))]
        self.assertEqual(aligned_read.opt('XD'), x)
        # print (aligned_read.tags)

    def testNegativeIntegersWrittenToFile(self):
        r = self.build_read()
        x = -2
        r.tags = [("XD", x)]
        with get_temp_context("negative_integers.bam") as fn:
            with pysam.AlignmentFile(fn,
                                     "wb",
                                     referencenames=("chr1",),
                                     referencelengths=(1000,)) as outf:
                outf.write(r)
            with pysam.AlignmentFile(fn) as inf:
                r = next(inf)
            self.assertEqual(r.tags, [("XD", x)])


class TestCopy(ReadTest):

    def testCopy(self):
        a = self.build_read()
        b = copy.copy(a)
        # check if a and be are the same
        self.assertEqual(a, b)

        # check if they map to different objects
        a.query_name = 'ReadA'
        b.query_name = 'ReadB'
        self.assertEqual(a.query_name, 'ReadA')
        self.assertEqual(b.query_name, 'ReadB')

    def testDeepCopy(self):
        a = self.build_read()
        b = copy.deepcopy(a)
        # check if a and be are the same
        self.assertEqual(a, b)

        # check if they map to different objects
        a.query_name = 'ReadA'
        b.query_name = 'ReadB'
        self.assertEqual(a.query_name, 'ReadA')
        self.assertEqual(b.query_name, 'ReadB')


class TestSetTagGetTag(ReadTest):

    def check_tag(self, tag, value, value_type, alt_value_type=None):
        a = self.build_read()
        a.set_tag(tag, value, value_type=value_type)
        v, t = a.get_tag(tag, with_value_type=True)
        self.assertEqual(v, value)

        if alt_value_type:
            self.assertEqual(t, alt_value_type)
        else:
            self.assertEqual(t, value_type)

    def test_set_tag_with_A(self):
        self.check_tag('TT', "x", value_type="A")

    def test_set_tag_with_a(self):
        self.check_tag('TT', "x", value_type="a", alt_value_type="A")

    def test_set_tag_with_C(self):
        self.check_tag('TT', 12, value_type="C")

    def test_set_tag_with_c(self):
        self.check_tag('TT', 12, value_type="c")

    def test_set_tag_with_S(self):
        self.check_tag('TT', 12, value_type="S")

    def test_set_tag_with_s(self):
        self.check_tag('TT', 12, value_type="s")

    def test_set_tag_with_I(self):
        self.check_tag('TT', 12, value_type="I")

    def test_set_tag_with_i(self):
        self.check_tag('TT', 12, value_type="i")

    def test_set_tag_with_f(self):
        self.check_tag('TT', 2.5, value_type="f")

    def test_set_tag_with_d(self):
        self.check_tag('TT', 2.5, value_type="d")

    def test_set_tag_with_H(self):
        self.check_tag('TT', "AE12", value_type="H")

    def test_set_tag_with_automated_type_detection(self):
        self.check_tag('TT', -(1 << 7), value_type=None, alt_value_type="c")
        self.check_tag('TT', -(1 << 7) - 1,
                       value_type=None, alt_value_type="s")
        self.check_tag('TT', -(1 << 15), value_type=None, alt_value_type="s")
        self.check_tag('TT', -(1 << 15) - 1,
                       value_type=None, alt_value_type="i")
        self.check_tag('TT', -(1 << 31), value_type=None, alt_value_type="i")
        self.assertRaises(
            ValueError,
            self.check_tag,
            'TT',
            -(1 << 31) - 1,
            value_type=None,
            alt_value_type="i")

        self.check_tag('TT', (1 << 8) - 1, value_type=None, alt_value_type="C")
        self.check_tag('TT', (1 << 8), value_type=None, alt_value_type="S")
        self.check_tag('TT', (1 << 16) - 1,
                       value_type=None, alt_value_type="S")
        self.check_tag('TT', (1 << 16), value_type=None, alt_value_type="I")
        self.check_tag('TT', (1 << 32) - 1,
                       value_type=None, alt_value_type="I")
        self.assertRaises(
            ValueError,
            self.check_tag,
            'TT',
            (1 << 32),
            value_type=None,
            alt_value_type="I")


class TestSetTagsGetTag(TestSetTagGetTag):

    def check_tag(self, tag, value, value_type, alt_value_type=None):
        a = self.build_read()
        a.set_tags([(tag, value, value_type)])
        v, t = a.get_tag(tag, with_value_type=True)
        if alt_value_type:
            self.assertEqual(t, alt_value_type)
        else:
            self.assertEqual(t, value_type)
        self.assertEqual(v, value)


class TestEnums(unittest.TestCase):

    def test_cigar_enums_are_defined(self):
        self.assertEqual(pysam.CMATCH, 0)
        self.assertEqual(pysam.CINS, 1)
        self.assertEqual(pysam.CDEL, 2)
        self.assertEqual(pysam.CREF_SKIP, 3)
        self.assertEqual(pysam.CSOFT_CLIP, 4)
        self.assertEqual(pysam.CHARD_CLIP, 5)
        self.assertEqual(pysam.CPAD, 6)
        self.assertEqual(pysam.CEQUAL, 7)
        self.assertEqual(pysam.CDIFF, 8)
        self.assertEqual(pysam.CBACK, 9)

    def test_sam_flags_are_defined(self):
        self.assertEqual(pysam.FPAIRED, 1)
        self.assertEqual(pysam.FPROPER_PAIR, 2)
        self.assertEqual(pysam.FUNMAP, 4)
        self.assertEqual(pysam.FMUNMAP, 8)
        self.assertEqual(pysam.FREVERSE, 16)
        self.assertEqual(pysam.FMREVERSE, 32)
        self.assertEqual(pysam.FREAD1, 64)
        self.assertEqual(pysam.FREAD2, 128)
        self.assertEqual(pysam.FSECONDARY, 256)
        self.assertEqual(pysam.FQCFAIL, 512)
        self.assertEqual(pysam.FDUP, 1024)
        self.assertEqual(pysam.FSUPPLEMENTARY, 2048)


class TestBuildingReadsWithoutHeader(unittest.TestCase):
    
    def build_read(self):
        '''build an example read, but without header information.'''

        a = pysam.AlignedSegment()
        a.query_name = "read_12345"
        a.query_sequence = "ATGC" * 10
        a.flag = 0
        a.reference_id = -1
        a.reference_start = 20
        a.mapping_quality = 20
        a.cigartuples = ((0, 10), (2, 1), (0, 9), (1, 1), (0, 20))
        a.next_reference_id = 0
        a.next_reference_start = 200
        a.template_length = 167
        a.query_qualities = pysam.qualitystring_to_array("1234") * 10
        # todo: create tags
        return a

    def test_read_can_be_constructed_without_header(self):
        read = self.build_read()
        self.assertEqual(read.query_name, "read_12345")
            
    def test_reference_id_can_be_set(self):
        read = self.build_read()
        read.reference_id = 2
        self.assertEqual(read.reference_id, 2)

    def test_reference_name_is_not_available(self):
        read = self.build_read()
        self.assertRaises(
            ValueError,
            setattr,
            read,
            "reference_name",
            "chr2")
        
    def test_read_can_be_written_to_file(self):
        tmpfilename = get_temp_filename(".bam")
        with pysam.AlignmentFile(tmpfilename, "wb",
                                 reference_names=["chr1", "chr2", "chr3"],
                                 reference_lengths=[1000, 2000, 3000]) as outf:
            read = self.build_read()
            read.reference_id = 2
            outf.write(read)

        stdout = pysam.samtools.view(tmpfilename)
        chromosome = stdout.split("\t")[2]
        self.assertEqual(chromosome, "chr3")
        os.unlink(tmpfilename)


class TestForwardStrandValues(ReadTest):

    def test_sequence_is_complemented(self):
        a = self.build_read()
        a.is_reverse = False
        fwd_seq = a.query_sequence
        
        rev_seq = fwd_seq.translate(maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))[::-1]
        self.assertEqual(fwd_seq, a.get_forward_sequence())
        a.is_reverse = True        
        self.assertEqual(fwd_seq, a.query_sequence)
        self.assertEqual(rev_seq, a.get_forward_sequence())

    def test_qualities_are_complemented(self):
        a = self.build_read()
        a.is_reverse = False
        fwd_qual = a.query_qualities
        rev_qual = fwd_qual[::-1]
        self.assertEqual(fwd_qual, a.get_forward_qualities())
        a.is_reverse = True        
        self.assertEqual(fwd_qual, a.query_qualities)
        self.assertEqual(rev_qual, a.get_forward_qualities())


class TestExportImport(ReadTest):

    def test_string_export(self):
        a = self.build_read()
        self.assertEqual(a.to_string(),
                         "read_12345\t0\tchr1\t21\t20\t10M1D9M1I20M\t=\t201\t167\t"
                         "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\t1234123412341234123412341234123412341234")

    def test_string_export_import_without_tags(self):
        a = self.build_read()
        a.tags = []
        b = pysam.AlignedSegment.fromstring(a.to_string(), a.header)
        self.assertEqual(a, b)

    def test_string_export_import_with_tags(self):
        a = self.build_read()
        a.tags = [("XD", 12), ("RF", "abc")]
        b = pysam.AlignedSegment.fromstring(a.to_string(), a.header)
        self.assertEqual(a, b)
        
    def test_to_string_without_alignment_file(self):
        with open(os.path.join(BAM_DATADIR, "ex2.sam")) as samf:
            reference = [x[:-1] for x in samf if not x.startswith("@")]

        with pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex2.bam"), "r") as pysamf:
            for s, p in zip(reference, pysamf):
                self.assertEqual(s, p.to_string())
                
    def test_dict_export(self):
        a = self.build_read()
        a.tags = [("XD", 12), ("RF", "abc")]
        
        self.assertEqual(
            a.to_dict(),
            json.loads(
                '{"name": "read_12345", "flag": "0", "ref_name": "chr1", "ref_pos": "21", '
                '"map_quality": "20", "cigar": "10M1D9M1I20M", "next_ref_name": "=", '
                '"next_ref_pos": "201", "length": "167", '
                '"seq": "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC", '
                '"qual": "1234123412341234123412341234123412341234", "tags": ["XD:i:12", "RF:Z:abc"]}'))

    def test_string_export_import_without_tags(self):
        a = self.build_read()
        a.tags = []
        b = pysam.AlignedSegment.from_dict(a.to_dict(), a.header)
        self.assertEqual(a, b)

    def test_string_export_import_with_tags(self):
        a = self.build_read()
        a.tags = [("XD", 12), ("RF", "abc")]
        b = pysam.AlignedSegment.from_dict(a.to_dict(), a.header)
        self.assertEqual(a, b)


class TestFillMD(unittest.TestCase):
    def test_fillmd(self):
        """
        Ensure MD tags come out identical to those generated with calmd
        """

        bam_no_md_path = os.path.join(BAM_DATADIR, "ex2.bam")
        bam_with_md_path = os.path.join(BAM_DATADIR, "with_md.bam")
        fasta_path = os.path.join(BAM_DATADIR, "ex1.fa")

        with pysam.AlignmentFile(bam_no_md_path, "rb") as pysamf_no_md, \
            pysam.AlignmentFile(bam_with_md_path, "rb") as pysamf_with_md, \
            pysam.FastaFile(fasta_path) as fa:

            # bam file has multiple chromosomes, so load them all upfront for easy
            # look up
            reference_fa = {
                'chr1': fa.fetch('chr1', return_bytes=True),
                'chr2': fa.fetch('chr2', return_bytes=True),
            }
 
            for read in pysamf_no_md:
                expected_read = next(pysamf_with_md)


                self.assertFalse(read.has_tag('MD'))
                read.fill_md(reference_fa[read.reference_name])

                # unmapped reads don't get an MD tag
                if read.is_unmapped:
                    self.assertFalse(read.has_tag('MD'))
                else:
                    self.assertTrue(read.has_tag('MD'))
                    self.assertEqual(read.get_tag('MD'), expected_read.get_tag('MD'))

        
if __name__ == "__main__":
    unittest.main()
