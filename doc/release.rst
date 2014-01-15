=============
Release notes
=============

Release 0.7.8
=============

   * added AlignedRead.setTag method

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
   * re-organized code to permit linking againts pysam
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

