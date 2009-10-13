*********
Glossary
*********

.. glossary::
   :sorted:

   cigar
      An alignment format string. In the python API, the cigar alignment is 
      presented as a list of tuples ``(operation,length)``. For example, the tuple
      ``[ (0,3), (1,5), (0,2) ]`` refers to an alignment with 3 matches, 5 insertions
      and another 2 matches.

   region
      A region string according to samtools convention, for example
      ``chr1:10-200`` and ``chr1``. The region can not be empty.

   column
      Reads that are aligned to a base in the :term:`target` sequence.
     
   target
      The sequence that reads have been mapped onto. Synonymous to
      :term:`reference`.

   bam file
       A file in bam format. Bam files are binary formatted, indexed and 
       allow random access.

   tam file
       A file in tam format. Tam files are human readable files of 
       tab-separated fields. Tam files do not allow random access.

   sam file
       A file containing aligned reads. The :term:`sam file` can either
       be a :term:`bam file` or a :term:`tam file`.

   pileup
      Pileup     

   samtools
      The samtools_ package.

   csamtools
      The samtools_ C-API.

   fetching
      Retrieving all mapped reads mapped to a :term:`region`.

   tid
      A :term:`target` id. This numerical identifier is used internally
      by samtools and can be translated into a :term:`reference` using
      the method :meth:`pysam.Samfile.getTarget`.

   reference
      A reference sequence id, for example ``chr1``, ``contig123``.

.. _samtools: http://samtools.sourceforge.net
