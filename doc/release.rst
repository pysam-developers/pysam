=============
Release notes
=============

Release 0.8.4
=============

This release contains numerous bugfixes and a first implementation of
a pythonic interface to VCF/BCF files. Note that this code is still
incomplete and preliminary, but does offer a nearly complete immutable
Pythonic interface to VCF/BCF metadata and data with reading and
writing capability.

Potential isses when upgrading from v0.8.3:

* binary tags are now returned as python arrays
* renamed several methods for pep8 compatibility, old names still retained for	
  backwards compatibility, but should be considered deprecated.
   * gettid() is now get_tid()
   * getrname() is now get_reference_name()
   * parseRegion() is now parse_region()
* faidx now returns strings and not binary strings in py3.
* The cython components have been broken up into smaller files with
  more specific content. This will affect users using the cython
  interfaces.

Edited list of commit log changes:

*    fixes AlignmentFile.check_index to return True
*    add RG/PM header tag - closes #179
*    add with_seq option to get_aligned_pairs
*    use char * inside reconsituteReferenceSequence
*    add soft clipping for get_reference_sequence
*    add get_reference_sequence
*    queryEnd now computes length from cigar string if no sequence present, closes #176
*    tolerate missing space at end of gtf files, closes #162
*    do not raise Error when receiving output on stderr
*    add docu about fetching without index, closes #170
*    FastaFile and FastxFile now return strings in python3, closes #173
*    py3 compat: relative -> absolute imports.
*    add reference_name and next_reference_name attributes to AlignedSegment
*    add function signatures to cvcf cython.  Added note about other VCF code.
*    add context manager functions to FastaFile
*    add reference_name and next_reference_name attributes to AlignedSegment
*    PileupColumn also gets a reference_name attribute.
*    add context manager functions to FastaFile
*    TabixFile.header for remote files raises AttributeError, fixes #157
*    add context manager interface to TabixFile, closes #165
*    change ctypedef enum to typedef enum for cython 0.23
*    add function signatures to cvcf cython, also added note about other VCF code
*    remove exception for custom upper-case header record tags.
*    rename VALID_HEADER_FIELDS to KNOWN_HEADER_FIELDS
*    fix header record tag parsing for custom tags.
*    use cython.str in count_coverage, fixes #141
*    avoid maketrans (issues with python3)
*    refactoring: AlignedSegment now in separate module
*    do not execute remote tests if URL not available
*    fix the unmapped count, incl reads with no SQ group
*    add raw output to tags
*    added write access for binary tags
*    bugfix in call to resize
*    implemented writing of binary tags from arrays
*    implemented convert_binary_tag to use arrays
*    add special cases for reads that are unmapped or whose mates are unmapped.
*    rename TabProxies to ctabixproxies
*    remove underscores from utility functions
*    move utility methods into cutils
*    remove callback argument to fetch - closes #128
*    avoid calling close in dealloc
*    add unit tests for File object opening
*    change AlignmentFile.open to filepath_or_object
*    implement copy.copy, close #65
*    add chaching of array attributes in AlignedSegment, closes #121
*    add export of Fastafile
*    remove superfluous pysam_dispatch
*    use persist option in FastqFile
*    get_tag: expose tag type if requested with `with_value_type`
*    fix to allow reading vcf record info via tabix-based vcf reader
*    add pFastqProxy and pFastqFile objects to make it possible to work with multiple fastq records per file handle, unlike FastqProxy/FastqFile.
*    release GIL around htslib IO operations
*    More work on read/write support, API improvements
*    add `phased` property on `VariantRecordSample`
*    add mutable properties to VariantRecord
*    BCF fixes and start of read/write support
*    VariantHeaderRecord objects now act like mappings for attributes.
*    add VariantHeader.alts dict from alt ID->Record.
*    Bug fix to strong representation of structured header records.
*    VariantHeader is now mutable


Release 0.8.3
=============

* samtools command now accept the "catch_stdout" option.

* get_aligned_pairs now works for soft-clipped reads.

* query_position is now None when a PileupRead is not aligned
  to a particular position.

* AlignedSegments are now comparable and hashable.

Release 0.8.2.1
===============

* Installation bugfix release.

Release 0.8.2
=============

* Pysam now wraps htslib 1.2.1 and samtools version 1.2.

* Added CRAM file support to pysam.

* New alignment info interface.
   * opt() and setTag are deprecated, use get_tag() and set_tag()
     instead.
   * added has_tag()
   * tags is deprecated, use get_tags() and set_tags() instead.

* FastqFile is now FastxFile to reflect that the latter permits
  iteration over both fastq- and fasta-formatted files.

* A Cython wrapper for htslib VCF/BCF reader/writer. The wrapper
  provides a nearly complete Pythonic interface to VCF/BCF metadata
  with reading and writing capability. However, the interface is still
  incomplete and preliminary and lacks capability to mutate the
  resulting data.
  
Release 0.8.1
=============

* Pysam now wraps htslib and samtools versions 1.1.

* Bugfixes, most notable:
  * issue #43: uncompressed BAM output
  * issue #42: skip tests requiring network if none available
  * issue #19: multiple iterators can now be made to work on the same tabix file
  * issue #24: All strings returned from/passed to the pysam API are now unicode in python 3
  * issue #5:  type guessing for lists of integers fixed    
    
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
     * cigar -> cigartuples
     * cigarstring -> cigarstring
     * tlen -> template_length
     * seq -> query_sequence
     * qual -> query_qualities, now returns array
     * qqual -> query_alignment_qualities, now returns array
     * tags -> tags
     * alen -> reference_length, reference is always "alignment", so removed
     * aend -> reference_end
     * rlen -> query_length
     * query -> query_alignment_sequence
     * qstart -> query_alignment_start
     * qend -> query_alignment_end
     * qlen -> query_alignment_length
     * mrnm -> next_reference_id   
     * mpos -> next_reference_start
     * rname -> reference_id
     * isize -> template_length
     * blocks -> get_blocks()
     * aligned_pairs -> get_aligned_pairs()
     * inferred_length -> infer_query_length()
     * positions -> get_reference_positions()
     * overlap() -> get_overlap()

  * All strings are now passed to or received from the pysam API
    as strings, no more bytes.

Other changes:
   * AlignmentFile.fetch(reopen) option is now multiple_iterators. The
     default changed to not reopen a file unless requested by the user.
   * FastaFile.getReferenceLength is now FastaFile.get_reference_length

Backwards incompatible changes

* Empty cigarstring now returns None (intstead of '')
* Empty cigar now returns None (instead of [])
* When using the extension classes in cython modules, AlignedRead
  needs to be substituted with AlignedSegment. 
* fancy_str() has been removed
* qual, qqual now return arrays




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
