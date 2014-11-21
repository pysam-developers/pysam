========
Glossary
========

.. glossary::
   :sorted:

   cigar
      An alignment format string. In the python API, the cigar alignment is 
      presented as a list of tuples ``(operation,length)``. For example, the tuple
      ``[ (0,3), (1,5), (0,2) ]`` refers to an alignment with 3 matches, 5 insertions
      and another 2 matches.

   region
      A genomic region, stated relative to a reference sequence. A region consists of 
      reference name ('chr1'), start (100000), and end (200000). 
      0-based coordinates. Can be expressed as a string ('chr1:10000:20000')

   column
      Reads that are aligned to a base in the :term:`reference` sequence.
     
   tid
      The :term:`target` id. The target id is 0 or a positive integer mapping to
      entries within the sequence dictionary in the header section of 
      a :term:`TAM` file or :term:`BAM` file.

   Reference
      The sequence that a :term:`tid` refers to. For example ``chr1``, ``contig123``.

   SAM
       A textual format for storing genomic alignment information.

   BAM
       Binary SAM format. BAM files are binary formatted, indexed and 
       allow random access.

   TAM
       Text SAM file. TAM files are human readable files of 
       tab-separated fields. TAM files do not allow random access.

   sam file
       A file containing aligned reads. The :term:`sam file` can either
       be a :term:`BAM` file or a :term:`TAM` file.

   pileup
      Pileup     

   samtools
      The samtools_ package.

   csamtools
      The samtools_ C-API.

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
      of the alignment record. This is in difference to hard clipped reads.

   hard clipping
   hard clipped

      In hard clipped reads, part of the sequence has been removed
      prior to alignment. That only a subsequence is aligend might be
      recorded in the :term:`cigar` alignment, but the removed
      sequence will not be part of the alignment record, in contrast
      to :term:`soft clipped` reads.
     
