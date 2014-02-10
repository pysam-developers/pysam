#!/usr/bin/env python
'''unit testing code for pysam.

Execute in the :file:`tests` directory as it requires the Makefile
and data files located there.
'''

import sys, os, shutil, gzip
import pysam
import unittest
import itertools
import subprocess

class TestVCFIterator( unittest.TestCase ):

    filename = "example.vcf40.gz"
    columns = ("contig", "pos", "id", 
               "ref", "alt", "qual", 
               "filter", "info", "format" )

    def testRead( self ):
        self.vcf = pysam.VCF()
        self.vcf.connect( self.filename )
        
        for x in self.vcf.fetch():
            print str(x)
            print x.pos
            print x.alt
            print x.id
            print x.qual
            print x.filter
            print x.info
            print x.format

            for s in x.samples:
                print s, x[s]
        
if __name__ == "__main__":

    unittest.main()


def Test():

    vcf33 = """##fileformat=VCFv3.3
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=NS,1,Integer,"Number of Samples With Data"
##INFO=DP,1,Integer,"Total Depth"
##INFO=AF,-1,Float,"Allele Frequency"
##INFO=AA,1,String,"Ancestral Allele"
##INFO=DB,0,Flag,"dbSNP membership, build 129"
##INFO=H2,0,Flag,"HapMap2 membership"
##FILTER=q10,"Quality below 10"
##FILTER=s50,"Less than 50% of samples have data"
##FORMAT=GT,1,String,"Genotype"
##FORMAT=GQ,1,Integer,"Genotype Quality"
##FORMAT=DP,1,Integer,"Read Depth"
##FORMAT=HQ,2,Integer,"Haplotype Quality"
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
20\t14370\trs6054257\tG\tA\t29\t0\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:-1,-1
17\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3:-1,-1
20\t1110696\trs6040355\tA\tG,T\t67\t0\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4:-1,-1
17\t1230237\t.\tT\t.\t47\t0\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2:-1,-1
20\t1234567\tmicrosat1\tG\tD4,IGA\t50\t0\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3"""

    vcf40 = """##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
M\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t1234567\tmicrosat1\tGTCT\tG,GTACT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
17\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4"""

    if False:
        print "Parsing v3.3 file:"
        print vcf33
        vcf = VCFFile()
        lines = [data for data in vcf.parse( (line+"\n" for line in vcf33.split('\n') ) )]
        print "Writing v3.3 file:"
        vcf.write( sys.stdout, lines )

    if False:
        print "Parsing v4.0 file:"
        print vcf40
        vcf = VCFFile()
        lines = [data for data in vcf.parse( (line+"\n" for line in vcf40.split('\n') ) )]
        print "Writing v4.0 file:"
        vcf.write( sys.stdout, lines )

    if True:
        print "Parsing v3.3 file:"
        print vcf33
        vcf = sortedVCFFile()
        lines = [data for data in vcf.parse( (line+"\n" for line in vcf33.split('\n') ) )]
        print "Writing v3.3 file:"
        vcf.write( sys.stdout, lines )

    if True:
        print "Parsing v4.0 file:"
        print vcf40
        vcf = sortedVCFFile()
        lines = [data for data in vcf.parse( (line+"\n" for line in vcf40.split('\n') ) )]
        print "Writing v4.0 file:"
        vcf.write( sys.stdout, lines )
