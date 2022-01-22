========
Glossary
========

.. glossary::
   :sorted:

   cigar
      Stands for Compact Idiosyncratic Gapped Alignment Report and
      represents a compressed (run-length encoded) pairwise alignment
      format.  It was first defined by the Exonerate Aligner, but was alter
      adapted and adopted as part of the :term:`SAM` standard and many other
      aligners.  In the Python API, the cigar alignment is presented as a
      list of tuples ``(operation,length)``.  For example, the tuple ``[
      (0,3), (1,5), (0,2) ]`` refers to an alignment with 3 matches, 5
      insertions and another 2 matches.

   region
      A genomic region, stated relative to a :term:`reference` sequence. A
      region consists of reference name ('chr1'), start (15000), and
      end (20000). Start and end can be omitted for regions spanning
      a whole chromosome. If ``end`` is missing, the region will span from
      ``start`` to the end of the chromosome. Within pysam, coordinates
      are 0-based half-open intervals, i.e., the first base of the
      reference sequence is numbered zero; and the base at position
      ``start`` is part of the interval, but the base at ``end`` is not.

      When a region is written as a single string using
      `samtools`_-compatible notation, e.g., 'chr1:15001-20000',
      the string's coordinates instead represent a 1-based closed interval,
      i.e., both (1-based) positions 15,001 and 20,000 are part of the
      interval. (This example denotes the same 5,000-base region as the
      example in the previous paragraph.)

   genotype
      An individual's collection of genes. It can also refer to the two alleles
	  inherited for a particular gene.

   column
      The portion of reads aligned to a single base in the 
	  :term:`reference` sequence.

   tid
      The :term:`target` id. The target id is 0 or a positive integer mapping to
      entries within the sequence dictionary in the header section of
      a :term:`TAM` file or :term:`BAM` file.

   contig
      The sequence that a :term:`tid` refers to. For example ``chr1``, ``contig123``.

   reference
      Synonym for contig.

   BED
      Browser Extensible Data format. A text file format used to store genomic
      :term:`regions<region>` as coordinates and associated notations.

   GTF
      The Gene Transfer Format is a file format used to hold information
	  about gene structure.

   SAM
       A textual format for storing genomic alignment information.

   BAM
       Binary SAM format. BAM files are binary formatted, indexed and
       allow random access.

   CRAM
       CRAM is a binary format representing the same sequence alignment
       information as SAM and BAM, but offering significantly better
       lossless compression than BAM.

   TAM
       Text SAM file. TAM files are human readable files of
       tab-separated fields. TAM files do not allow random access.

   sam file
       A file containing aligned reads. The :term:`sam file` can either
       be a :term:`BAM` file or a :term:`TAM` file.

   pileup
      Pileup

   fetching
      Retrieving all mapped reads mapped to a :term:`region`.

   target
      The sequence that a read has been aligned to. Target
      sequences have bot a numerical identifier (:term:`tid`)
      and an alphanumeric name (:term:`Reference`).

   tabix file
      A sorted, compressed and indexed tab-separated file created
      by the command line tool :file:`tabix` or the commands
      :meth:`tabix_compress` and :meth:`tabix_index`. The file
      is indexed by chromosomal coordinates.

   tabix row
      A row in a :term:`tabix file`. Fields within a row are
      tab-separated.

   soft clipping
   soft clipped

      In alignments with soft clipping part of the query sequence
      are not aligned. The unaligned query sequence is still part
      of the alignment record. This is in difference to
      :term:`hard clipped` reads.

   hard clipping
   hard clipped

      In hard clipped reads, part of the sequence has been removed
      prior to alignment. That only a subsequence is aligend might be
      recorded in the :term:`cigar` alignment, but the removed
      sequence will not be part of the alignment record, in contrast
      to :term:`soft clipped` reads.

   VCF
      Variant Call Format.

   BCF
      Binary :term:`VCF`.

   FASTA
      Simple text format containing sequence data, with only the bare
      minimum of metadata. Typically used for reference sequence data.

   FASTQ
      Simple text format containing sequence data and associated base
      qualities.

   tabix
      Utility in the htslib package to index :term:`bgzip` compressed
      files.

   faidx
      Utility in the `samtools`_ package to index :term:`fasta` formatted
      files.

   bgzip
      Utility in the htslib package to block compress genomic data
      files.
