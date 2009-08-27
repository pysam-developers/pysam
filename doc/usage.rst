*****
Usage
*****

In order to follow the examples, execute in the :file:`test` directory::
   
   make


Opening a sampfile
------------------

The basic operation is to open a :class:`Samfile`::

   import pysam

   samfile = pysam.Samfile()
   samfile.open( "ex1.bam", "rb" )

Fetching mapped reads in a :term:`region`
-----------------------------------------

There are two ways to process reads in a region. The
first method follows the :term:`csamtools` API and  works 
via a callback function. The callback will be executed for each 
alignment in a :term:`region`::

   def my_fetch_callback( alignment ):
       print str(alignment)

   samfile.fetch( "seq1:10-20", my_fetch_callback )

The second method uses python iterators. The iterator
:class:`IteratorRow` will iterate through mapped reads
and return a :class:`Alignment` object for each::

   iter = pysam.IteratorRow( samfile, "seq1:10-20")
   for x in iter: print str(x)

Note that both cases can iterate through a :term:`bam file`
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
return for each base in the target sequence the reads that
map to that particular position.

Again, there are two principal methods to use this iteration.
The first works via a callback function::

   def my_pileup_callback( pileups ):
       print str(pileups)
   samfile.pileup( "seq1:10-20", my_pileup_callback )

The second method uses python iterators. The iterator
:class:`IteratorColumn` will iterate through each :term:`column`
(target bases) and return list of aligned reads::

   iter = pysam.IteratorRow( samfile, "seq1:10-20")
   for x in iter: print str(x)


