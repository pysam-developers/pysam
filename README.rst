=====
Pysam
=====

Pysam is a python module for reading and manipulating files in the
SAM/BAM format. The SAM/BAM format is a way to store efficiently large
numbers of alignments (`Li 2009`_), such as those routinely created by
next-generation sequencing methods.

Pysam is a lightweight wrapper of the samtools_ C-API. Pysam also includes an
interface for tabix_.

The latest version is available through 
`pypi <https://pypi.python.org/pypi/pysam>`_. To install, simply
type::
  
   pip install pysam
                                                      .
Pysam documentation is available through https://readthedocs.org/ from
`here <http://pysam.readthedocs.org/en/latest/>`_

Questions and comments are very welcome and should be sent to the
`pysam user group <http://groups.google.com/group/pysam-user-group>`_

.. _samtools: http://samtools.sourceforge.net/
.. _tabix: http://samtools.sourceforge.net/tabix.shtml
.. _Li 2009: http://www.ncbi.nlm.nih.gov/pubmed/19505943
