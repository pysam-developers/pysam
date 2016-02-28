.. _Usage: 

=========================================
Working with BAM/CRAM/SAM-formatted files
=========================================

Opening a file
==============

To begin with, import the pysam module and open a
:class:`pysam.AlignmentFile`::

   import pysam
   samfile = pysam.AlignmentFile("ex1.bam", "rb")

The above command opens the file :file:`ex1.bam` for reading.
The ``b`` qualifier indicates that this is a :term:`BAM` file. 
To open a :term:`SAM` file, type::

   import pysam
   samfile = pysam.AlignmentFile("ex1.sam", "r")

:term:`CRAM` files are identified by a ``c`` qualifier::

   import pysam
   samfile = pysam.AlignmentFile("ex1.cram", "rc")

Fetching reads mapped to a :term:`region`
=========================================

Reads are obtained through a call to the
:meth:`pysam.AlignmentFile.fetch` method which returns an iterator.
Each call to the iterator will returns a :class:`pysam.AlignedSegment`
object::

   iter = samfile.fetch("seq1", 10, 20)
   for x in iter:
       print (str(x))

:meth:`pysam.AlignmentFile.fetch` returns all reads overlapping a
region sorted by the first aligned base in the :term:`reference`
sequence.  Note that it will also return reads that are only partially
overlapping with the :term:`region`. Thus the reads returned might
span a region that is larger than the one queried.

Using the pileup-engine
=======================

In contrast to :term:`fetching`, the :term:`pileup` engine returns for
each base in the :term:`reference` sequence the reads that map to that
particular position. In the typical view of reads stacking vertically
on top of the reference sequence similar to a multiple alignment,
:term:`fetching` iterates over the rows of this implied multiple
alignment while a :term:`pileup` iterates over the :term:`columns`.

Calling :meth:`~pysam.AlignmentFile.pileup` will return an iterator
over each :term:`column` (reference base) of a specified
:term:`region`. Each call to the iterator returns an object of the
type :class:`pysam.PileupColumn` that provides access to all the
reads aligned to that particular reference position as well as
some additional information::

   iter = samfile.pileup('seq1', 10, 20)
   for x in iter:
      print (str(x))
 

Creating BAM/CRAM/SAM files from scratch
========================================

The following example shows how a new :term:`BAM` file is constructed
from scratch.  The important part here is that the
:class:`pysam.AlignmentFile` class needs to receive the sequence
identifiers. These can be given either as a dictionary in a header
structure, as lists of names and sizes, or from a template file.
Here, we use a header dictionary::

   header = { 'HD': {'VN': '1.0'},
               'SQ': [{'LN': 1575, 'SN': 'chr1'}, 
                      {'LN': 1584, 'SN': 'chr2'}] }

   with pysam.AlignmentFile(tmpfilename, "wb", header=header) as outf:
       a = pysam.AlignedSegment()
       a.query_name = "read_28833_29006_6945"
       a.query_sequence="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
       a.flag = 99
       a.reference_id = 0
       a.reference_start = 32
       a.mapping_quality = 20
       a.cigar = ((0,10), (2,1), (0,25))
       a.next_reference_id = 0
       a.next_reference_start=199
       a.template_length=167
       a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
       a.tags = (("NM", 1),
		 ("RG", "L1"))
       outf.write(a)

Using streams
=============

Pysam does not support reading and writing from true python file
objects, but it does support reading and writing from stdin and
stdout. The following example reads from stdin and writes to stdout::

   infile = pysam.AlignmentFile("-", "r")
   outfile = pysam.AlignmentFile("-", "w", template=infile)
   for s in infile:
       outfile.write(s)

It will also work with :term:`BAM` files. The following script
converts a :term:`BAM` formatted file on stdin to a :term:`SAM`
formatted file on stdout::

   infile = pysam.AlignmentFile("-", "rb")
   outfile = pysam.AlignmentFile("-", "w", template=infile)
   for s in infile:
       outfile.write(s)

Note that the file open mode needs to changed from ``r`` to ``rb``.

=====================================
Using samtools commands within python
=====================================

Commands available in :term:`csamtools` are available
as simple function calls. For example::

   pysam.sort("ex1.bam", "output")

corresponds to the command line::

   samtools sort ex1.bam output

Command line options can be provided as arguments::
   
   pysam.sort("-n", "ex1.bam", "output")

or::

   pysam.sort("-m", "1000000", "ex1.bam", "output")

In order to get usage information, try::

   print pysam.sort.usage()

Argument errors raise a :class:`pysam.SamtoolsError`::

   pysam.sort()

   Traceback (most recent call last):
   File "x.py", line 12, in <module>
     pysam.sort()
   File "/build/lib.linux-x86_64-2.6/pysam/__init__.py", line 37, in __call__
     if retval: raise SamtoolsError( "\n".join( stderr ) )
   pysam.SamtoolsError: 'Usage: samtools sort [-n] [-m <maxMem>] <in.bam> <out.prefix>\n'

Messages from :term:`csamtools` on stderr are captured and are
available using the :meth:`getMessages` method::

   pysam.sort.getMessage()

Note that only the output from the last invocation of a command is
stored.

In order for pysam to make the output of samtools commands accessible
the stdout stream needs to be redirected. This is the default
behaviour, but can cause problems in environments such as the ipython
notebook. A solution is to pass the ``catch_stdout`` keyword
argument::

   pysam.sort(catch_stdout=False)

Note that this means that output from commands which produce output on
stdout will not be available. The only solution is to run samtools
commands through subprocess.

================================
Working with tabix-indexed files
================================

To open a tabular file that has been indexed with tabix_, use
:class:`~pysam.TabixFile`::

    import pysam
    tbx = pysam.TabixFile("example.bed.gz")

Similar to :class:`~pysam.AlignmentFile.fetch`, intervals within a
region can be retrieved by calling :meth:`~pysam.TabixFile.fetch()`::

    for row in tbx.fetch("chr1", 1000, 2000):
         print (str(row))

This will return a tuple-like data structure in which columns can
be retrieved by numeric index:

    for row in tbx.fetch("chr1", 1000, 2000):
         print ("chromosome is", row[0])

By providing a parser argument to :class:`~pysam.AlignmentFile.fetch`
or :class:`~pysam.TabixFile`, the data will we presented in parsed
form:

    for row in tbx.fetch("chr1", 1000, 2000, parser=pysam.asTuple()):
         print ("chromosome is", row.contig)

.. Currently inactivated as pileup deprecated
.. Using the samtools SNP caller
.. -----------------------------

.. There are two ways to access the samtools SNP caller. The :class:`pysam.IteratorSNPCalls`
.. is appropriate when calling many consecutive SNPs, while :class:`pysam.SNPCaller` is
.. best when calling SNPs at non-consecutive genomic positions. Each snp caller returns objects of
.. type :class:`pysam.SNPCall`.

.. To use :class:`pysam.IteratorSNPCalls`, associate it with a :class:`pysam.IteratorColumn`::

..     samfile = pysam.AlignmentFile( "ex1.bam", "rb")  
..     fastafile = pysam.Fastafile( "ex1.fa" )
..     pileup_iter = samfile.pileup( stepper = "samtools", fastafile = fastafile )
..     sncpall_iter = pysam.IteratorSNPCalls(pileup_iter)
..     for call in snpcall_iter:
..         print str(call)

.. Usage of :class:`pysam.SNPCaller` is similar::

..     samfile = pysam.AlignmentFile( "ex1.bam", "rb")  
..     fastafile = pysam.Fastafile( "ex1.fa" )
..     pileup_iter = samfile.pileup( stepper = "samtools", fastafile = fastafile )
..     snpcaller = pysam.SNPCaller.call(pileup_iter)
..     print snpcaller( "chr1", 100 )

.. Note the use of the option *stepper* to control which reads are included in the 
.. in the :term:`pileup`. The ``samtools`` stepper implements the same read selection
.. and processing as in the samtools pileup command.

.. Calling indels works along the same lines, using the :class:`pysam.IteratorIndelCalls`
.. and :class:`pysam.IteratorIndelCaller`.


====================================
Working with VCF/BCF formatted files
====================================

To iterate through a VCF/BCF formatted file use
:class:`~pysam.VariantFile`::

   from pysam import VariantFile

   bcf_in = VariantFile("test.bcf")  # auto-detect input format
   bcf_out = VariantFile('-', 'w', header=bcf_in.header)
   
   for rec in bcf_in.fetch('chr1', 100000, 200000):
       bcf_out.write(rec)

:meth:`_pysam.VariantFile.fetch()` iterates over
:class:`~pysam.VariantRecord` objects which provides access to
simple variant attributes such as :class:`~pysam.VariantRecord.contig`,
:class:`~pysam.VariantRecord.pos`, :class:`~pysam.VariantRecord.ref`::

   for rec in bcf_in.fetch():
       print (rec.pos)

but also to complex attributes such as the contents to the
:term:`info`, :term:`format` and :term:`genotype` columns. These
complex attributes are views on the underlying htslib data structures
and provide dictionary-like access to the data::

   for rec in bcf_in.fetch():
       print (rec.info)
       print (rec.info.keys())
       print (rec.info["DP"])

The :py:attr:`~pysam.VariantFile.header` attribute
(:class:`~pysam.VariantHeader`) provides access information
stored in the :term:`vcf` header. The complete header can be printed::

   >>> print (bcf_in.header)
   ##fileformat=VCFv4.2
   ##FILTER=<ID=PASS,Description="All filters passed">
   ##fileDate=20090805
   ##source=myImputationProgramV3.1
   ##reference=1000GenomesPilot-NCBI36
   ##phasing=partial
   ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples
   With Data">
   ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
   ##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
   ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
   ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build
   129">
   ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
   ##FILTER=<ID=q10,Description="Quality below 10">
   ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
   ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
   ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
   ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
   ##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
   ##contig=<ID=M>
   ##contig=<ID=17>
   ##contig=<ID=20>
   ##bcftools_viewVersion=1.3+htslib-1.3
   ##bcftools_viewCommand=view -O b -o example_vcf42.bcf
   example_vcf42.vcf.gz
   #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO   FORMAT    NA00001 NA00002 NA0000
  
Individual contents such as contigs, info fields, samples, formats can
be retrieved as attributes from :py:attr:`~pysam.VariantFile.header`::

   >>> print (bcf_in.header.contigs)
   <pysam.cbcf.VariantHeaderContigs object at 0xf250f8>

To convert these views to native python types, iterate through the views::

   >>> print list((bcf_in.header.contigs))
   ['M', '17', '20']
   >>> print list((bcf_in.header.filters))
   ['PASS', 'q10', 's50']
   >>> print list((bcf_in.header.info))
   ['NS', 'DP', 'AF', 'AA', 'DB', 'H2']
   >>> print list((bcf_in.header.samples))
   ['NA00001', 'NA00002', 'NA00003']

Alternatively, it is possible to iterate through all records in the
header returning objects of type :py:class:`~pysam.VariantHeaderRecord`:: ::

   >>> for x in bcf_in.header.records:
   >>>    print (x)
   >>>    print (x.type, x.key)
   GENERIC fileformat
   FILTER FILTER
   GENERIC fileDate
   GENERIC source
   GENERIC reference
   GENERIC phasing
   INFO INFO
   INFO INFO
   INFO INFO
   INFO INFO
   INFO INFO
   INFO INFO
   FILTER FILTER
   FILTER FILTER
   FORMAT FORMAT
   FORMAT FORMAT
   FORMAT FORMAT
   FORMAT FORMAT
   CONTIG contig
   CONTIG contig
   CONTIG contig
   GENERIC bcftools_viewVersion
   GENERIC bcftools_viewCommand

===============
Extending pysam
===============

Using pyximport_, it is (relatively) straight-forward to access pysam
internals and the underlying samtools library. An example is provided
in the :file:`tests` directory. The example emulates the samtools
flagstat command and consists of three files:

1. The main script :file:`pysam_flagstat.py`. The important lines in
   this script are::

      import pyximport
      pyximport.install()
      import _pysam_flagstat

      ...
   
      flag_counts = _pysam_flagstat.count(pysam_in)

   The first part imports, sets up pyximport_ and imports the cython
   module :file:`_pysam_flagstat`.  The second part calls the
   ``count`` method in :file:`_pysam_flagstat`.
 
2. The cython implementation :file:`_pysam_flagstat.pyx`. This script
   imports the pysam API via::

      from pysam.calignmentfile cimport AlignmentFile, AlignedSegment

   This statement imports, amongst others, :class:`AlignedSegment`
   into the namespace. Speed can be gained from declaring
   variables. For example, to efficiently iterate over a file, an
   :class:`AlignedSegment` object is declared::

      # loop over samfile
      cdef AlignedSegment read
      for read in samfile:
          ...

3. A :file:`pyxbld` providing pyximport_ with build information.
   Required are the locations of the samtools and pysam header
   libraries of a source installation of pysam plus the
   :file:`csamtools.so` shared library. For example::

     def make_ext(modname, pyxfilename):
	 from distutils.extension import Extension
	 import pysam
	 return Extension(name=modname,
               sources=[pyxfilename],
               extra_link_args=pysam.get_libraries(),
	       include_dirs=pysam.get_include(),
	       define_macros=pysam.get_defines())

If the script :file:`pysam_flagstat.py` is called the first time,
pyximport_ will compile the cython_ extension
:file:`_pysam_flagstat.pyx` and make it available to the
script. Compilation requires a working compiler and cython_
installation.  Each time :file:`_pysam_flagstat.pyx` is modified, a
new compilation will take place.

pyximport_ comes with cython_.

