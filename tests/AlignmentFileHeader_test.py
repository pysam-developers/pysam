#!/usr/bin/env python
'''unit testing code for pysam.

Execute in the :file:`tests` directory as it requires the Makefile
and data files located there.
'''

import unittest
import os
import shutil
import sys
import re
import copy
import collections
from collections import OrderedDict as odict
import subprocess
import logging
import array
if sys.version_info.major >= 3:
    from io import StringIO
else:
    from StringIO import StringIO

from functools import partial

import pysam
import pysam.samtools
from TestUtils import checkBinaryEqual, checkURL, \
    check_samtools_view_equal, checkFieldEqual, force_str, \
    get_temp_filename, BAM_DATADIR


class TestHeaderConstruction(unittest.TestCase):
    """testing header construction."""

    header_dict = odict(
        [('SQ', [odict([('LN', 1575), ('SN', 'chr1'), ('AH', 'chr1:5000000-5010000')]),
                 odict([('LN', 1584), ('SN', 'chr2'), ('AH', '*')])]),
         ('RG', [odict([('LB', 'SC_1'), ('ID', 'L1'), ('SM', 'NA12891'),
                        ('PU', 'SC_1_10'), ("CN", "name:with:colon")]),
                 odict([('LB', 'SC_2'), ('ID', 'L2'), ('SM', 'NA12891'),
                        ('PU', 'SC_2_12'), ("CN", "name:with:colon")])]),
         ('PG', [odict([('ID', 'P1'), ('VN', '1.0')]),
                 odict([('ID', 'P2'), ('VN', '1.1')])]),
         ('HD', odict([('VN', '1.0')])),
         ('CO', ['this is a comment', 'this is another comment']),
        ])

    header_text = ("@HD\tVN:1.0\n"
                   "@SQ\tSN:chr1\tLN:1575\tAH:chr1:5000000-5010000\n"
                   "@SQ\tSN:chr2\tLN:1584\tAH:*\n"
                   "@RG\tID:L1\tPU:SC_1_10\tLB:SC_1\tSM:NA12891\tCN:name:with:colon\n"
                   "@RG\tID:L2\tPU:SC_2_12\tLB:SC_2\tSM:NA12891\tCN:name:with:colon\n"
                   "@PG\tID:P1\tVN:1.0\n"
                   "@PG\tID:P2\tVN:1.1\n"
                   "@CO\tthis is a comment\n"
                   "@CO\tthis is another comment\n")


    header_from_references = odict(
        [('SQ', [odict([('LN', 1575), ('SN', 'chr1')]),
                 odict([('LN', 1584), ('SN', 'chr2')])]),
         ('RG', [odict([('LB', 'SC_1'), ('ID', 'L1'), ('SM', 'NA12891'),
                        ('PU', 'SC_1_10'), ("CN", "name:with:colon")]),
                 odict([('LB', 'SC_2'), ('ID', 'L2'), ('SM', 'NA12891'),
                        ('PU', 'SC_2_12'), ("CN", "name:with:colon")])]),
         ('PG', [odict([('ID', 'P1'), ('VN', '1.0')]),
                 odict([('ID', 'P2'), ('VN', '1.1')])]),
         ('HD', odict([('VN', '1.0')])),
         ('CO', ['this is a comment', 'this is another comment']),
        ])

    header_without_text = odict(
        [('SQ', [odict([('LN', 1575), ('SN', 'chr1')]),
                 odict([('LN', 1584), ('SN', 'chr2')])]),
        ])
    
    def compare_headers(self, test_header, ref_header=None):
        '''compare two headers a and b.'''
        test_header_dict = test_header.as_dict()
        if ref_header is None:
            ref_header = self.header_dict
            
        for ak, av in test_header_dict.items():
            self.assertTrue(ak in self.header_dict, "key '%s' not in '%s' " % (ak, ref_header))
            self.assertEqual(av, ref_header[ak])
        for ak, av in ref_header.items():
            self.assertTrue(ak in test_header_dict, "key '%s' not in '%s' " % (ak, test_header_dict))
            self.assertEqual(av, test_header_dict[ak])

    def check_name_mapping(self, test_header):
        for x, y in enumerate(("chr1", "chr2")):
            tid = test_header.get_tid(y)
            ref = test_header.get_reference_name(x)
            self.assertEqual(tid, x)
            self.assertEqual(ref, y)

        self.assertEqual(test_header.get_tid("chr?"), -1)
        self.assertRaises(ValueError, test_header.get_reference_name, 2)
            
    def test_header_constructed_from_dict(self):
        header = pysam.AlignmentHeader.from_dict(self.header_dict)
        self.compare_headers(header)
        self.check_name_mapping(header)
        
    def test_header_constructed_from_text(self):
        header = pysam.AlignmentHeader.from_text(self.header_text)
        self.compare_headers(header)
        self.check_name_mapping(header)
        
    def test_header_constructed_from_header(self):
        header = pysam.AlignmentHeader.from_text(self.header_text)
        self.compare_headers(header.copy())
        self.check_name_mapping(header)

    def test_header_constructed_from_references(self):
        text = re.sub("@SQ[^\n]+\n", "", self.header_text)
        assert "@SQ" not in text
        header = pysam.AlignmentHeader.from_references(
            reference_names=["chr1", "chr2"],
            reference_lengths=[1575, 1584],
            text=text)
        self.compare_headers(header, self.header_from_references)
        self.check_name_mapping(header)

    def test_header_constructed_from_references_without_text(self):
        header = pysam.AlignmentHeader.from_references(
            reference_names=["chr1", "chr2"],
            reference_lengths=[1575, 1584])
        self.compare_headers(header, self.header_without_text)
        self.check_name_mapping(header)
        
        
class TestHeaderSAM(unittest.TestCase):
    """testing header manipulation"""

    header = {'SQ': [{'LN': 1575, 'SN': 'chr1', 'AH': 'chr1:5000000-5010000'},
                     {'LN': 1584, 'SN': 'chr2', 'AH': '*'}],
              'RG': [{'LB': 'SC_1', 'ID': 'L1', 'SM': 'NA12891',
                      'PU': 'SC_1_10', "CN": "name:with:colon"},
                     {'LB': 'SC_2', 'ID': 'L2', 'SM': 'NA12891',
                      'PU': 'SC_2_12', "CN": "name:with:colon"}],
              'PG': [{'ID': 'P1', 'VN': '1.0'}, {'ID': 'P2', 'VN': '1.1'}],
              'HD': {'VN': '1.0'},
              'CO': ['this is a comment', 'this is another comment'],
              }

    def compare_headers(self, a, b):
        '''compare two headers a and b.'''
        for ak, av in a.items():
            self.assertTrue(ak in b, "key '%s' not in '%s' " % (ak, b))
            self.assertEqual(av, b[ak])

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.sam"),
            "r")

    def test_header_content_is_as_expected(self):
        self.compare_headers(self.header, self.samfile.header.to_dict())
        self.compare_headers(self.samfile.header.to_dict(), self.header)

    def test_text_access_works(self):
        self.assertEqual(self.samfile.text, self.samfile.header.__str__())
        
    def test_name_mapping(self):
        for x, y in enumerate(("chr1", "chr2")):
            tid = self.samfile.gettid(y)
            ref = self.samfile.getrname(x)
            self.assertEqual(tid, x)
            self.assertEqual(ref, y)

        self.assertEqual(self.samfile.gettid("chr?"), -1)
        self.assertRaises(ValueError, self.samfile.getrname, 2)

    def tearDown(self):
        self.samfile.close()


class TestHeaderBAM(TestHeaderSAM):

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.bam"),
            "rb")


class TestHeaderCRAM(TestHeaderSAM):

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.cram"),
            "rc")

    def compare_headers(self, a, b):
        '''compare two headers a and b.'''
        def _strip(dd):
            for x in dd:
                for y in ("M5", "UR"):
                    if y in x:
                        del x[y]
        for ak, av in a.items():
            _strip(av)
            self.assertTrue(ak in b, "key '%s' not in '%s' " % (ak, b))
            _strip(b[ak])

            self.assertEqual(av, b[ak])


class TestHeaderFromRefs(unittest.TestCase):
    '''see issue 144

    reference names need to be converted to string for python 3
    '''

    # def testHeader( self ):
    #     refs = ['chr1', 'chr2']
    #     tmpfile = "tmp_%i" % id(self)
    #     s = pysam.AlignmentFile(tmpfile, 'wb',
    #                       referencenames=refs,
    #                       referencelengths=[100]*len(refs))
    #     s.close()

    #     self.assertTrue( checkBinaryEqual( 'issue144.bam', tmpfile ),
    #                      'bam files differ')
    #     os.unlink( tmpfile )



class TestHeaderWriteRead(unittest.TestCase):
    header = {'SQ': [{'LN': 1575, 'SN': 'chr1'},
                     {'LN': 1584, 'SN': 'chr2'}],
              'RG': [{'LB': 'SC_1', 'ID': 'L1', 'SM': 'NA12891',
                      'PU': 'SC_1_10', "CN": "name:with:colon"},
                     {'LB': 'SC_2', 'ID': 'L2', 'SM': 'NA12891',
                      'PU': 'SC_2_12', "CN": "name:with:colon"}],
              'PG': [{'ID': 'P1', 'VN': '1.0', 'CL': 'tool'},
                     {'ID': 'P2', 'VN': '1.1', 'CL': 'tool with in option -R a\tb',
                      'PP': 'P1'}],
              'HD': {'VN': '1.0'},
              'CO': ['this is a comment', 'this is another comment'],
              }

    def compare_headers(self, a, header_b):
        '''compare two headers a and b.

        Ignore M5 and UR field as they are set application specific.
        '''
        b = header_b.to_dict()
        for ak, av in a.items():
            self.assertTrue(ak in b, "key '%s' not in '%s' " % (ak, b))
            self.assertEqual(
                len(av), len(b[ak]),
                "unequal number of entries for key {}: {} vs {}"
                .format(ak, av, b[ak]))

            for row_a, row_b in zip(av, b[ak]):
                if isinstance(row_b, dict):
                    for x in ["M5", "UR"]:
                        try:
                            del row_b[x]
                        except KeyError:
                            pass
                self.assertEqual(row_a, row_b)

    def check_read_write(self, flag_write, header):

        fn = get_temp_filename()
        with pysam.AlignmentFile(
                fn,
                flag_write,
                header=header,
                reference_filename=os.path.join(BAM_DATADIR, "ex1.fa")) as outf:
            a = pysam.AlignedSegment()
            a.query_name = "abc"
            outf.write(a)

        with pysam.AlignmentFile(fn) as inf:
            read_header = inf.header

        os.unlink(fn)
        self.compare_headers(header, read_header)

    def test_SAM(self):
        self.check_read_write("wh", self.header)

    def test_BAM(self):
        self.check_read_write("wb", self.header)

    def test_CRAM(self):
        header = copy.copy(self.header)
        # for CRAM, \t needs to be quoted:
        header['PG'][1]['CL'] = re.sub(r"\t", r"\\\\t", header['PG'][1]['CL'])
        self.check_read_write("wc", header)
