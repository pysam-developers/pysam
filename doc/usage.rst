.. _Usage: 

====================================
Working with BAM/SAM-formatted files
====================================

In order to follow the examples below, execute
   
   make

in the :file:`test` directory.

Opening a samfile
=================

The first operation is to import the pysam module and open a 
:class:`pysam.Samfile`::

   import pysam
   samfile = pysam.Samfile( "ex1.bam", "rb" )

The above command opens the file :file:`ex1.bam` for reading.
The ``b`` qualifier indicates that this is a :term:`BAM` file. 
To open a :term:`SAM` file, type::

   import pysam
   samfile = pysam.Samfile( "ex1.bam", "r" )

Fetching reads mapped to a :term:`region`
=========================================

There are two ways to obtain the reads mapped to a genomic region. The
first method follows the :term:`csamtools` API and  works 
via a callback function. The callback will be executed for each 
alignment in a :term:`region`::

   def my_fetch_callback( alignment ):
       print str(alignment)

   samfile.fetch( 'seq1', 10, 20, callback = my_fetch_callback )

Using a function object, work can be done on the alignments. The
code below simply counts aligned reads::

   class Counter:
       def __init__(self):
           self.counts = 0
       def __call__(self, alignment):
           self.counts += 1
   
   c = Counter()
   samfile.fetch( 'seq1', 10, 20, callback = c )
   print "counts=", c.counts

The second method uses python iterators. If you call :meth:`pysam.Samfile.fetch`
without a callback, an iterator of the type :class:`pysam.IteratorRow` is returned.
It will iterate through mapped reads
and return a :class:`pysam.AlignedRead` object for each::

   iter = samfile.fetch( 'seq1', 10, 20)
   for x in iter: print str(x)

Note that both methods iterate through a :term:`BAM` file
on a read-by-read basis. 

:meth:`pysam.Samfile.fetch` returns all reads overlapping a region sorted
by the first aligned base in the :term:`reference` sequence.
Note that it will also return reads that are only partially
overlapping with the :term:`region`. Thus the reads returned
might span a region that is larger than the one queried.

Using the pileup-engine
=======================

The :term:`pileup` engine of :term:`csamtools` iterates
over all reads that are aligned to a :term:`region`. In
contrast to :term:`fetching`, the :term:`pileup` engine 
returns for each base in the :term:`reference` sequence the reads that
map to that particular position.

Again, there are two principal methods to iterate.
The first works via a callback function::

   def my_pileup_callback( pileups ):
       print str(pileups)
   samfile.pileup( 'seq1', 10, 20, callback = my_pileup_callback )

The second method uses python iterators. The iterator
:class:`pysam.IteratorColumn` will iterate through each :term:`column`
(reference bases) and return a list of aligned reads::

   iter = samfile.pileup( 'seq1', 10, 20 )
   for x in iter: print str(x)

Aligned reads are returned as a :class:`pysam.PileupColumn`.

Using samtools commands within python
=====================================

Commands available in :term:`csamtools` are available
as simple function calls. For example::

   pysam.sort( "ex1.bam", "output" )

corresponds to the command line::

   samtools sort ex1.bam output

Command line options can be provided as arguments::
   
   pysam.sort( "-n", "ex1.bam", "output" )

or::

   pysam.sort( "-m", "1000000", "ex1.bam", "output" )

In order to get usage information, try::

   print pysam.sort.usage()

Argument errors raise a :class:`pysam.SamtoolsError`::

   pysam.sort()

   Traceback (most recent call last):
   File "x.py", line 12, in <module>
     pysam.sort()
   File "/home/andreas/pysam/build/lib.linux-x86_64-2.6/pysam/__init__.py", line 37, in __call__
     if retval: raise SamtoolsError( "\n".join( stderr ) )
   pysam.SamtoolsError: 'Usage: samtools sort [-n] [-m <maxMem>] <in.bam> <out.prefix>\n'

Messages from :term:`csamtools` on stderr are captured and are
available using the :meth:`getMessages` method::

   pysam.sort.getMessage()

Note that only the output from the last invocation of a command
is stored.

Creating SAM/BAM files from scratch
===================================

The following example shows how a new BAM file is constructed from scratch.
The important part here is that the :class:`pysam.Samfile` class needs to receive
the sequence identifiers. These can be given either as a dictionary in a
header structure, as lists of names and sizes, or from a template file.
Here, we use a header dictionary::

   header = { 'HD': {'VN': '1.0'},
               'SQ': [{'LN': 1575, 'SN': 'chr1'}, 
                      {'LN': 1584, 'SN': 'chr2'}] }

   outfile = pysam.Samfile( tmpfilename, "wh", header = header )
   a = pysam.AlignedRead()
   a.qname = "read_28833_29006_6945"
   a.seq="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
   a.flag = 99
   a.rname = 0
   a.pos = 32
   a.mapq = 20
   a.cigar = ( (0,10), (2,1), (0,25) )
   a.mrnm = 0
   a.mpos=199
   a.isize=167
   a.qual="<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<"
   a.tags = ( ("NM", 1),
	      ("RG", "L1") )
   outfile.write(a)
   outfile.close()

Using streams
=============

Pysam does not support reading and writing from true python file objects, but
it does support reading and writing from stdin and stdout. The following example reads 
from stdin and writes to stdout::

   infile = pysam.Samfile( "-", "r" )
   outfile = pysam.Samfile( "-", "w", template = infile )
   for s in infile: outfile.write(s)

It will also work with :term:`BAM` files. The following script converts a :term:`BAM` formatted file
on stdin to a :term:`SAM` formatted file on stdout::

   infile = pysam.Samfile( "-", "rb" )
   outfile = pysam.Samfile( "-", "w", template = infile )
   for s in infile: outfile.write(s)

Note, only the file open mode needs to changed from ``r`` to ``rb``.

.. Currently inactivated as pileup deprecated
.. Using the samtools SNP caller
.. -----------------------------

.. There are two ways to access the samtools SNP caller. The :class:`pysam.IteratorSNPCalls`
.. is appropriate when calling many consecutive SNPs, while :class:`pysam.SNPCaller` is
.. best when calling SNPs at non-consecutive genomic positions. Each snp caller returns objects of
.. type :class:`pysam.SNPCall`.

.. To use :class:`pysam.IteratorSNPCalls`, associate it with a :class:`pysam.IteratorColumn`::

..     samfile = pysam.Samfile( "ex1.bam", "rb")  
..     fastafile = pysam.Fastafile( "ex1.fa" )
..     pileup_iter = samfile.pileup( stepper = "samtools", fastafile = fastafile )
..     sncpall_iter = pysam.IteratorSNPCalls(pileup_iter)
..     for call in snpcall_iter:
..         print str(call)

.. Usage of :class:`pysam.SNPCaller` is similar::

..     samfile = pysam.Samfile( "ex1.bam", "rb")  
..     fastafile = pysam.Fastafile( "ex1.fa" )
..     pileup_iter = samfile.pileup( stepper = "samtools", fastafile = fastafile )
..     snpcaller = pysam.SNPCaller.call(pileup_iter)
..     print snpcaller( "chr1", 100 )

.. Note the use of the option *stepper* to control which reads are included in the 
.. in the :term:`pileup`. The ``samtools`` stepper implements the same read selection
.. and processing as in the samtools pileup command.

.. Calling indels works along the same lines, using the :class:`pysam.IteratorIndelCalls`
.. and :class:`pysam.IteratorIndelCaller`.

Extending pysam
===============

Using pyximport_, it is (relatively) straight-forward to access pysam
internals and the underlying samtools library. An example is provided
in the :file:`test` directory. The example emulates the samtools flagstat command
and consists of three files:

1. The main script :file:`pysam_flagstat.py`. The important lines in this script
   are::

      import pyximport
      pyximport.install()
      import _pysam_flagstat

      ...
   
      flag_counts = _pysam_flagstat.count( pysam_in )

   The first part imports, sets up pyximport_ and imports the cython module :file:`_pysam_flagstat`. 
   The second part calls the ``count`` method in :file:`_pysam_flagstat`. 
 
2. The cython implementation :file:`_pysam_flagstat.pyx`. This script imports the pysam API via::

      from csamtools cimport *

   This statement imports, amongst others, :class:`AlignedRead` into the namespace. Speed can be
   gained from declaring variables. For example, to efficiently iterate
   over a file, an :class:`AlignedRead` object is declared::

      # loop over samfile
      cdef AlignedRead read
      for read in samfile:
          ...

3. A :file:`pyxbld` providing pyximport_ with build information. 
   Required are the locations of the samtools and pysam header libraries 
   of a source installation of pysam plus the :file:`csamtools.so` 
   shared library. For example::

     def make_ext(modname, pyxfilename):
	 from distutils.extension import Extension
	 import pysam, os
	 dirname = os.path.dirname( pysam.__file__ )[:-len("pysam")]
	 return Extension(name = modname,
			  sources=[pyxfilename],
			  extra_link_args=[ os.path.join( dirname, "csamtools.so")],
			  include_dirs =  pysam.get_include(),
			  define_macros = pysam.get_defines() )

If the script :file:`pysam_flagstat.py` is called the first time, pyximport_ will 
compile the cython_ extension :file:`_pysam_flagstat.pyx` and make it available 
to the script. Compilation requires a working compiler and cython_ installation.
Each time :file:`_pysam_flagstat.pyx` is modified, a new compilation will take place.

pyximport_ comes with cython_.

.. _cython: http://cython.org/

.. _pyximport: http://www.prescod.net/pyximport/
