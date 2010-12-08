.. _Usage: 

*****
Usage
*****

In order to follow the examples below, execute::
   
   make

in the :file:`test` directory.

Opening a sampfile
------------------

The first operation is to open a :class:`pysam.Samfile`::

   import pysam

   samfile = pysam.Samfile( "ex1.bam", "rb" )

Fetching mapped reads in a :term:`region`
-----------------------------------------

There are two ways to iterate through reads in a region. The
first method follows the :term:`csamtools` API and  works 
via a callback function. The callback will be executed for each 
alignment in a :term:`region`::

   def my_fetch_callback( alignment ):
       print str(alignment)

   samfile.fetch( 'seq1', 10, 20, callback = my_fetch_callback )

Using a function object, work can be done on the alignments. The
code below simply counts aligned reads::

   class Counter:
       mCounts = 0
       def __call__(self, alignment):
           self.mCounts += 1
   
   c = Counter()
   samfile.fetch( 'seq1', 10, 20, callback = c )
   print "counts=", c.mCounts

The second method uses python iterators. If you call :meth:`samtools.samfile.fetch`
without a callback, an iterator of the type :class:`pysam.IteratorRow` is returned.
It will iterate through mapped reads
and return a :class:`pysam.AlignedRead` object for each::

   iter = samfile.fetch( 'seq1', 10, 20)
   for x in iter: print str(x)

Note that both methods iterate through a :term:`BAM` file
on a read-by-read basis. They need not load all data into
memory.

Fetching returns all reads overlapping a region sorted
by the lowest aligned base in the :term:`reference` sequence.
Note that it will also return reads that are only partially
overlapping with the :term:`region`. Thus the reads returned
might span a region that is larger than the one queried.

Using the pileup-engine
-----------------------

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

Some samtools commands that create output on stdout are 
associated with parsers. For example, the pysam equivalent of
calling "samtools pileup -c"::

   for p in pysam.pileup( "-c", "ex1.bam" ):
      print str(p)
      
returns an iterator over SNP calls. The iterator return objects of
type :class:`pysam.PileupEntry`. The output of the two lines of code above 
is::

   PileupEntry(chromosome='seq1', position=2, reference_base='N', consensus_base='A', consensus_quality=27, snp_quality=0, rms_mapping_quality=60, coverage=1, read_bases='A', base_qualities='<')
   PileupEntry(chromosome='seq1', position=3, reference_base='N', consensus_base='C', consensus_quality=33, snp_quality=0, rms_mapping_quality=60, coverage=2, read_bases='C^~C', base_qualities='<<')
   PileupEntry(chromosome='seq1', position=4, reference_base='N', consensus_base='T', consensus_quality=33, snp_quality=0, rms_mapping_quality=60, coverage=2, read_bases='TT', base_qualities='<<')
   PileupEntry(chromosome='seq1', position=5, reference_base='N', consensus_base='A', consensus_quality=36, snp_quality=0, rms_mapping_quality=60, coverage=3, read_bases='AA^~A', base_qualities='<<<')
   PileupEntry(chromosome='seq1', position=6, reference_base='N', consensus_base='G', consensus_quality=39, snp_quality=0, rms_mapping_quality=60, coverage=4, read_bases='GGG^`G', base_qualities='<<<(')
   ...

Messages from :term:`csamtools` on stderr are captured and are
available using the :meth:`getMessages` method::

   pysam.pileup.getMessage()

Note that only the output from the last invocation of a command
is stored.

In order to get the unparsed output, use the *raw* argument::

   for p in pysam.pileup( "-c", "ex1.bam", raw=True ):
      print str(p),

Creating SAM/BAM files from scratch
-----------------------------------

The following example shows how a new BAM file is constructed from scratch.
The important part here is that the :class:`Samfile` class needs to receive
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
-------------

Pysam does not support reading and writing from true python file objects, but
it does support reading and writing from stdin and stdout. The following example reads 
from stdin and writes to stdout::

   infile = pysam.Samfile( "-", "r" )
   outfile = pysam.Samfile( "-", "w", template = infile )
   for s in infile: outfile.write(s)

It will also work with BAM files. The following script converts a BAM formatted file
on stdin to a SAM formatted file on stdout::

   infile = pysam.Samfile( "-", "rb" )
   outfile = pysam.Samfile( "-", "w", template = infile )
   for s in infile: outfile.write(s)


Note that only the file open mode needs to changed from ``r`` to ``rb``.
