=============
Release notes
=============

Release 0.8.1
=============

* Pysam now wraps htslib and samtools versions 1.1

* Bugfixes
  * issue #43: uncompressed BAM output
  * issue #42: skip tests requiring network if none available
    
* API changes for consistency. The old API is still present,
  but deprecated.
  In particular:

  * Tabixfile -> TabixFile
  * Fastafile -> FastaFile
  * Fastqfile -> FastqFile
  * Samfile -> AlignmentFile
  * AlignedRead -> AlignedSegment
     * qname -> query_name
     * tid -> reference_id
     * pos -> reference_start
     * mapq -> mapping_quality
     * rnext -> next_reference_id
     * pnext -> next_reference_start
     * cigar = alignment (now returns CigarAlignment object)
     * cigarstring = cigarstring
     * tlen -> query_length
     * seq -> query_sequence
     * qual -> query_qualities, now returns array
     * tags = tags (now returns Tags object)
     * alen -> reference_length, reference is always "alignment", so removed
     * aend -> reference_end
     * rlen -> query_length 
     * query -> query_alignment_sequence
     * qqual -> query_alignment_qualities, now returns array
     * qstart -> query_alignment_start
     * qend -> query_alignment_end
     * qlen -> query_alignment_length
     * mrnm -> next_reference_id   
     * mpos -> next_reference_start
     * rname -> reference_id
     * isize -> query_length
     * blocks -> getBlocks()
     * aligned_pairs -> getAlignedPairs()
     * inferred_length -> inferQueryLength()
     * positions -> getReferencePositions()
     * overlap() -> getOverlap()

  * improved Cigar handling
  * improved Tags handling
  * Text in SAM/BAM/Tabix files is now a string - no more bytes.

Backwards incompatible changes

* Empty cigarstring now returns None (intstead of '')
* Empty cigar now returns None (instead of [])
* When using the extension classes in cython modules, AlignedRead
  needs to be substituted with AlignedSegment. 
* fancy_str() has been removed


Release 0.8.0
=============

* Disabled features
  * IteratorColumn.setMask() disabled as htslib does not implement
    this functionality?

* Not implemented yet:
  * reading SAM files without header

Tabix files between version 0.7.8 and 0.8.0 are
not compatible and need to be re-indexed.

While version 0.7.8 and 0.8.0 should be mostly
compatible, there are some notable exceptions:

* tabix iterators will fail if there are comments
  in the middle or the end of a file.

* tabix raises always ValueError for invalid intervals.
  Previously, different types of errors were raised
  (KeyError, IndexError, ValueError) depending on
  the type of invalid intervals (missing chromosome,
  out-of-range, malformatted interval).


Release 0.7.8
=============

* added AlignedRead.setTag method
* added AlignedRead.blocks
* unsetting CIGAR strings is now possible
* empty CIGAR string returns empty list
* added reopen flag to Samfile.fetch()
* various bugfixes

Release 0.7.7
=============

* added Fastafile.references, .nreferences and .lengths
* tabix_iterator now uses kseq.h for python 2.7

Release 0.7.6
=============

* added inferred_length property
* issue 122: MACOSX getline missing, now it works?
* seq and qual can be set None
* added Fastqfile

Release 0.7.5
=============

* switch to samtools 0.1.19
* issue 122: MACOSX getline missing
* issue 130: clean up tempfiles
* various other bugfixes

Release 0.7.4
=============
	
* further bugfixes to setup.py and package layout

Release 0.7.3
=============
	
* further bugfixes to setup.py
* upgraded distribute_setup.py to 0.6.34

Release 0.7.2
=============
  
* bugfix in installer - failed when cython not present
* changed installation locations of shared libraries

Release 0.7.1
=============

* bugfix: missing PP tag PG records in header
* added pre-built .c files to distribution

Release 0.7
===========

* switch to tabix 0.2.6
* added cigarstring field
* python3 compatibility
* added B tag handling
* added check_sq and check_header options to Samfile.__init__
* added lazy GTF parsing to tabix
* reworked support for VCF format parsing
* bugfixes

Release 0.6
===========

* switch to samtools 0.1.18
* various bugfixes
* removed references to deprecated 'samtools pileup' functionality
* AlignedRead.tags now returns an empty list if there are no tags.
* added pnext, rnext and tlen

Release 0.5
===========

* switch to samtools 0.1.16 and tabix 0.2.5
* improved tabix parsing, added vcf support
* re-organized code to permit linking against pysam
* various bugfixes
* added Samfile.positions and Samfile.overlap

Release 0.4
===========

* switch to samtools 0.1.12a and tabix 0.2.3
* added snp and indel calling.
* switch from pyrex to cython
* changed handling of samtools stderr
* various bugfixes
* added Samfile.count and Samfile.mate
* deprecated AlignedRead.rname, added AlignedRead.tid

Release 0.3
===========

* switch to samtools 0.1.8
* added support for tabix files
* numerous bugfixes including
* permit simultaneous iterators on the same file
* working access to remote files
