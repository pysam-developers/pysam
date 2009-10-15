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
      A genomic region, stated relative to a reference sequence. Consists of reference name ('chr1'), start (100000), and end (200000). 0-based coordinates. Can be expressed as a string ('chr1:10000:20000')

   column
      Reads that are aligned to a base in the :term:`target` sequence.
     
   Reference
      The sequence that reads have been mapped onto. For example ``chr1``, ``contig123``.

   BAM
       Binary SAM format. BAM files are binary formatted, indexed and 
       allow random access.

   TAM
       Text SAM file. TAM files are human readable files of 
       tab-separated fields. TAM files do not allow random access.

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


.. _samtools: http://samtools.sourceforge.net
