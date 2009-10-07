#!/usr/bin/python

import unittest
import pysam
from commands import getoutput

class TestIteratorFunction(unittest.TestCase):

    # build all test files to be used
    def setUp(self):
        print getoutput("samtools view -bT ex1.fa ex1.sam.gz > ex1.bam")
        print getoutput("samtools index ex1.bam ex1.bam.bai")
        self.samfile=pysam.Samfile()
        self.samfile.open("ex1.bam","rb")

    # test that the IteratorRow returns as many alignments as samtools on the first contig
    def testIter1(self):
        iter=pysam.IteratorRow(self.samfile,"seq1")
        mCounts=len(list(iter))
        nCounts=int(getoutput("samtools view ex1.bam seq1 | wc -l").strip())
        self.assertEqual(mCounts, nCounts)

    # test that the iterator returns as many alignments as samtools on subsequent contigs
    def testIter2(self):
        iter=pysam.IteratorRow(self.samfile,"seq2")
        mCounts=len(list(iter))
        nCounts=int(getoutput("samtools view ex1.bam seq2 | wc -l").strip())
        self.assertEqual(mCounts, nCounts)

    # test that the iterator returns as many alignments as samtools on all contigs
    def testIterAll(self):
        iter=pysam.IteratorRowAll(self.samfile)
        mCounts=len(list(iter))
        nCounts=int(getoutput("samtools view ex1.bam | wc -l").strip())
        self.assertEqual(mCounts, nCounts)

    # tests that pysam doesn't die when it tries to grab alignments from an empty region
    def testGrabNoAligns(self):
        iter=pysam.IteratorRow(self.samfile,"seq1:1-10")
        mCounts=len(list(iter))
        nCounts=int(getoutput("samtools view ex1.bam seq1:1-10 | wc -l").strip())
        self.assertEqual(mCounts, nCounts)

    # test a bunch of attributes for a single read
    

if __name__ == '__main__':
    unittest.main()

