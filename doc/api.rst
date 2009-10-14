pysam - An interface for reading and writing SAM files
=====================================================



Pysam is a python module that makes it easy to read and manipulate mapped short read sequence data stored in SAM files. It's a lightweight wrapper of the samtools_ C-API.

To use the module to read a file in SAM format, just open the file with :class:`~pysam.Samfile`::

   import pysam
   samfile = pysam.Samfile( "sample1.bam", "rb" )
   

Now the file is open you can iterate over all of the read mapping to a specified region using :meth:`~pysam.Samfile.fetch`.
Each iteration returns a :class:`~pysam.AlignedRead` object which represents a single read along with its fields and optional tags::

   for alignedreadread in samfile.fetch('chr1', 100, 120):
	print alignedread
   samfile.close()


    EAS56_57:6:190:289:82	0	99	<<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<;	69	CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA	0	192	1
    EAS56_57:6:190:289:82	0	99	<<<<<<;<<<<<<<<<<;<<;<<<<;8<6;9;;2;	137	AGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC	73	64	1
    EAS51_64:3:190:727:308	0	102	<<<<<<<<<<<<<<<<<<<<<<<<<<<::<<<844	99	GGTGCAGAGCCGAGTCACGGGGTTGCCAGCACAGG	99	18	1
    ...



You can also write to a :class:`~pysam.Samfile`::

   import pysam
   samfile = pysam.Samfile("sample1.bam", "rb")
   pairedreads = pysam.Samfile("allpaired.bam", "wb")
   for read in samfile.fetch():
	if read.is_paired:
		pairedreads.write(read)
   pairedreads.close()
   samfile.close()


An alternative way of accessing the data in a SAM file is by iterating over each base of a specified region using the :meth:`~pysam.Samfile.pileup` method. Each iteration returns a :class:`~pysam.PileupColumn` which represents all the reads in the SAM file that map to a single base in the reference sequence. The list of reads are represented as :class:`~pysam.PileupRead` objects in the :attr:`PileupColumn.pileups <pysam.PileupColumn.pileups>` property::

    samfile = pysam.Samfile("/net/cpp-group/martin/projects/pysam/pysam/tests/ex1.bam", "rb" )
    for pileupcolumn in samfile.pileup( 'chr1', 100, 120):
	print
	print 'coverage at base %s = %s' % (pileupcolumn.pos , pileupcolumn.n)
	for pileupread in pileupcolumn.pileups:
	    print '\tbase in read %s = %s' % (pileupread.alignment.qname, pileupread.alignment.seq[pileupread.qpos])


    coverage at base 99 = 1
        base in read EAS56_57:6:190:289:82 = A

    coverage at base 100 = 1
        base in read EAS56_57:6:190:289:82 = G

    coverage at base 101 = 1
        base in read EAS56_57:6:190:289:82 = G

    coverage at base 102 = 2
        base in read EAS56_57:6:190:289:82 = G
        base in read EAS51_64:3:190:727:308 = G
    ...


.. _samtools: http://samtools.sourceforge.net/ 

C-API
*****

.. automodule:: pysam
   :members:
   :undoc-members:

Pileup
******

.. automodule:: pysam.Pileup
   :members:
   :undoc-members:
