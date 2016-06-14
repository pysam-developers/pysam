import os
import pysam
import unittest
import collections
import copy
import array

from TestUtils import checkFieldEqual


SAMTOOLS = "samtools"
WORKDIR = "pysam_test_work"
DATADIR = "pysam_data"


class ReadTest(unittest.TestCase):

    def buildRead(self):
        '''build an example read.'''

        a = pysam.AlignedSegment()
        a.query_name = "read_12345"
        a.query_sequence = "ACGT" * 10
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 20
        a.mapping_quality = 20
        a.cigartuples = ((0, 10), (2, 1), (0, 9), (1, 1), (0, 20))
        a.next_reference_id = 0
        a.next_reference_start = 200
        a.template_length = 167
        a.query_qualities = pysam.qualitystring_to_array("1234") * 10
        # todo: create tags
        return a


class TestAlignedSegment(ReadTest):

    '''tests to check if aligned read can be constructed
    and manipulated.
    '''

    def testEmpty(self):
        a = pysam.AlignedSegment()
        self.assertEqual(a.query_name, None)
        self.assertEqual(a.query_sequence, None)
        self.assertEqual(pysam.qualities_to_qualitystring(a.query_qualities), None)
        self.assertEqual(a.flag, 0)
        self.assertEqual(a.reference_id, 0)
        self.assertEqual(a.mapping_quality, 0)
        self.assertEqual(a.cigartuples, None)
        self.assertEqual(a.tags, [])
        self.assertEqual(a.next_reference_id, 0)
        self.assertEqual(a.next_reference_start, 0)
        self.assertEqual(a.template_length, 0)

    def testStrOfEmptyRead(self):
        a = pysam.AlignedSegment()
        s = str(a)
        self.assertEqual(
            "None\t0\t0\t0\t0\tNone\t0\t0\t0\tNone\tNone\t[]",
            s)

    def testSettingTagInEmptyRead(self):
        '''see issue 62'''
        a = pysam.AlignedSegment()
        a.tags = (("NM", 1),)
        a.query_qualities = None
        self.assertEqual(a.tags, [("NM", 1), ])

    def testCompare(self):
        '''check comparison functions.'''
        a = self.buildRead()
        b = self.buildRead()

        self.assertEqual(0, a.compare(b))
        self.assertEqual(0, b.compare(a))
        self.assertTrue(a == b)
        self.assertTrue(b == a)
        self.assertFalse(a != b)
        self.assertFalse(b != a)

        b.tid = 2
        self.assertFalse(a == b)
        self.assertFalse(b == a)
        self.assertTrue(a != b)
        self.assertTrue(b != a)

    def testHashing(self):
        a = self.buildRead()
        b = self.buildRead()
        self.assertEqual(hash(a), hash(b))
        b.tid = 2
        self.assertNotEqual(hash(a), hash(b))

    def testUpdate(self):
        '''check if updating fields affects other variable length data
        '''
        a = self.buildRead()
        b = self.buildRead()

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
        b.query_sequence = "ACGT"
        checkFieldEqual(self,
                        a, b,
                        ("query_sequence", "query_qualities", "query_length"))
        b.query_sequence = "ACGT" * 3
        checkFieldEqual(self,
                        a, b,
                        ("query_sequence", "query_qualities", "query_length"))
        b.query_sequence = "ACGT" * 10
        checkFieldEqual(self, a, b, ("query_qualities",))

        # reset qual
        b = self.buildRead()

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
        a = self.buildRead()
        a.query_sequence = a.query_sequence[5:10]
        self.assertEqual(pysam.qualities_to_qualitystring(a.query_qualities), None)

        a = self.buildRead()
        s = pysam.qualities_to_qualitystring(a.query_qualities)
        a.query_sequence = a.query_sequence[5:10]
        a.query_qualities = pysam.qualitystring_to_array(s[5:10])

        self.assertEqual(pysam.qualities_to_qualitystring(a.query_qualities), s[5:10])

    def testLargeRead(self):
        '''build an example read.'''

        a = pysam.AlignedSegment()
        a.query_name = "read_12345"
        a.query_sequence = "ACGT" * 200
        a.flag = 0
        a.reference_id = 0
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
        a = self.buildRead()
        oldlen = a.template_length
        oldlen *= 2
        a.template_length = oldlen
        self.assertEqual(a.template_length, oldlen)

    def testPositions(self):
        a = self.buildRead()
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
        a = self.buildRead()
        a.cigar = [(4, 30), (0, 20), (1, 3), (0, 47)]

        self.assertEqual(100,
                         len(a.get_reference_positions(full_length=True)))

    def testBlocks(self):
        a = self.buildRead()
        self.assertEqual(a.get_blocks(),
                         [(20, 30), (31, 40), (40, 60)])

    def test_get_aligned_pairs_soft_clipping(self):
        a = self.buildRead()
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
        a = self.buildRead()
        a.cigartuples = ((5, 2), (0, 35), (5, 3))
        self.assertEqual(a.get_aligned_pairs(),
                         # No seq, no seq pos
                         [(qpos, refpos) for (qpos, refpos) in zip(
                             range(0, 0 + 35), range(20, 20 + 35))])
        self.assertEqual(a.get_aligned_pairs(True),
                         [(qpos, refpos) for (qpos, refpos) in zip(
                             range(0, 0 + 35), range(20, 20 + 35))])

    def test_get_aligned_pairs_skip(self):
        a = self.buildRead()
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
        a = self.buildRead()
        a.cigartuples = ((7, 20), (8, 20))
        self.assertEqual(a.get_aligned_pairs(),
                         [(qpos, refpos) for (qpos, refpos) in zip(
                             range(0, 0 + 40), range(20, 20 + 40))])
        self.assertEqual(a.get_aligned_pairs(True),
                         [(qpos, refpos) for (qpos, refpos) in zip(
                             range(0, 0 + 40), range(20, 20 + 40))])

    def test_get_aligned_pairs_padding(self):
        a = self.buildRead()
        a.cigartuples = ((7, 20), (6, 1), (8, 19))

        def inner():
            a.get_aligned_pairs()
        # padding is not bein handled right now
        self.assertRaises(NotImplementedError, inner)

    def test_get_aligned_pairs(self):
        a = self.buildRead()
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

    def test_get_aligned_pairs_skip_reference(self):
        a = self.buildRead()
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
        a = self.buildRead()
        a.query_sequence = None
        a.cigarstring = "20M"
        self.assertEqual(a.query_alignment_length, 20)
        a.cigarstring = "20M1S"
        self.assertEqual(a.query_alignment_length, 20)
        a.cigarstring = "1S20M"
        self.assertEqual(a.query_alignment_length, 20)
        a.cigarstring = "1S20M1S"
        self.assertEqual(a.query_alignment_length, 20)


class TestCigarStats(ReadTest):
    
    def testStats(self):
        
        a = self.buildRead()

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
    filename = os.path.join(DATADIR, "example_aligned_pairs.bam")

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
        a = self.buildRead()
        self.assertRaises(KeyError, a.get_tag, "XP")

    def testEmptyTag(self):
        a = self.buildRead()
        self.assertRaises(KeyError, a.get_tag, "XT")

    def testSetTag(self):
        a = self.buildRead()
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
        read = self.buildRead()
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
        a = self.buildRead()
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
        self.assertRaises(ValueError, a.setTag, 'X6', 5.2, 'g')

    def testTagsUpdatingFloat(self):
        a = self.buildRead()
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
        a = self.buildRead()
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
            os.path.join(DATADIR, "ex8.bam"),
            "rb")

        for entry in samfile:
            before = entry.get_tags()
            entry.set_tags(before)
            after = entry.get_tags()
            self.assertEqual(after, before)

    def testMDTagMatchOnly(self):
        a = self.buildRead()

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
        a = self.buildRead()

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
        a = self.buildRead()

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
        a = self.buildRead()

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
        a = self.buildRead()

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
        a = self.buildRead()

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


class TestCopy(ReadTest):

    def testCopy(self):
        a = self.buildRead()
        b = copy.copy(a)
        # check if a and be are the same
        self.assertEqual(a, b)

        # check if they map to different objects
        a.query_name = 'ReadA'
        b.query_name = 'ReadB'
        self.assertEqual(a.query_name, 'ReadA')
        self.assertEqual(b.query_name, 'ReadB')

    def testDeepCopy(self):
        a = self.buildRead()
        b = copy.deepcopy(a)
        # check if a and be are the same
        self.assertEqual(a, b)

        # check if they map to different objects
        a.query_name = 'ReadA'
        b.query_name = 'ReadB'
        self.assertEqual(a.query_name, 'ReadA')
        self.assertEqual(b.query_name, 'ReadB')


class TestAsString(unittest.TestCase):

    def testAsString(self):
        with open(os.path.join(DATADIR, "ex2.sam")) as samf:
            reference = [x[:-1] for x in samf if not x.startswith("@")]

        with pysam.AlignmentFile(
            os.path.join(DATADIR, "ex2.bam"), "r") as pysamf:
            for s, p in zip(reference, pysamf):
                self.assertEqual(s, p.tostring(pysamf))

if __name__ == "__main__":
    unittest.main()
