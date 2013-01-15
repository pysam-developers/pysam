from pysam.csamtools cimport Samfile
from pysam.ctabix cimport Tabixfile
from pysam.cvcf cimport VCF
# from pysam.TabProxies cimport TupleProxy

cdef Samfile samfile
cdef Tabixfile tabixfile
cdef VCF vcf
cdef TupleProxy proxy

