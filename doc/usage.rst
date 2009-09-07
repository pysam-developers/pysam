*****
Usage
*****

In order to follow the examples, execute::
   
   make

in the :file:`test` directory.

Opening a sampfile
------------------

The first operation is to open a :class:`pysam.Samfile`::

   import pysam

   samfile = pysam.Samfile()
   samfile.open( "ex1.bam", "rb" )

Fetching mapped reads in a :term:`region`
-----------------------------------------

There are two ways to iterate through reads in a region. The
first method follows the :term:`csamtools` API and  works 
via a callback function. The callback will be executed for each 
alignment in a :term:`region`::

   def my_fetch_callback( alignment ):
       print str(alignment)

   samfile.fetch( "seq1:10-20", my_fetch_callback )

Using a function object, work can be done on the alignments. The
code below simply counts aligned reads::

   class Counter:
       mCounts = 0
       def __call__(self, alignment):
           self.mCounts += 1
   
   c = Counter()
   samfile.fetch( "seq1:10-20", c )
   print "counts=", c.mCounts

The second method uses python iterators. The iterator
:class:`pysam.IteratorRow` will iterate through mapped reads
and return a :class:`pysam.AlignedRead` object for each::

   iter = pysam.IteratorRow( samfile, "seq1:10-20")
   for x in iter: print str(x)

Note that both methods iterate through a :term:`bam file`
on a read-by-read basis. They need not load all data into
memory.

Fetching returns all reads overlapping a region sorted
by the lowest aligned base in the :term:`target` sequence.
Note that it will also return reads that are only partially
overlapping with the :term:`region`. Thus the reads returned
might span a region that is larger than the one queried.

Using the pileup-engine
-----------------------

The :term:`pileup` engine of :term:`csamtools` iterates
over all reads that are aligned to a :term:`region`. In
contrast to :term:`fetching`, the :term:`pileup` engine 
returns for each base in the target sequence the reads that
map to that particular position.

Again, there are two principal methods to iterate.
The first works via a callback function::

   def my_pileup_callback( pileups ):
       print str(pileups)
   samfile.pileup( "seq1:10-20", my_pileup_callback )

The second method uses python iterators. The iterator
:class:`pysam.IteratorColumn` will iterate through each :term:`column`
(target bases) and return a list of aligned reads::

   iter = pysam.IteratorRow( samfile, "seq1:10-20")
   for x in iter: print str(x)

Aligned reads are returned as a :class:`pysam.PileupColumn`.

Using samtools within python
----------------------------

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

