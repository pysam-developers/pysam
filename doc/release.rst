=============
Release notes
=============

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

