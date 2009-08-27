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
      A region string according to samtools convention. 

   column
      Reads that are aligned to a base in the :term:`target` sequence.
     
   target
      The sequence that reads have been mapped onto.

   bam file
       A file in bam format.

   pileup
      Pileup     

   samtools
      The samtools_ package.

   csamtools
      The samtools_ C-API.

   fetching
      Retrieving all mapped reads mapped to a :term:`region`.

.. _samtools: http://samtools.sourceforge.net
