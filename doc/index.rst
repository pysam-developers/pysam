pysam: htslib interface for python
==================================

:Author: Andreas Heger, Kevin Jacobs and contributors
:Date: |today|
:Version: |version|

Pysam is a python module for reading, manipulating and writing
genomic data sets.

Pysam is a wrapper of the htslib_ C-API and provides facilities to
read and write SAM/BAM/VCF/BCF/BED/GFF/GTF/FASTA/FASTQ files as well
as access to the command line functionality of the samtools_ and
bcftools_ packages. The module supports compression and random access
through indexing.

This module provides a low-level wrapper around the htslib_ C-API as
using cython and a high-level, pythonic API for convenient access to
the data within genomic file formats. 

The current version wraps *htslib-1.16*, *samtools-1.16.1*, and *bcftools-1.16*.

To install the latest release, type::

    pip install pysam

See the :ref:`Installation notes <installation>` for details.

Contents
--------

.. toctree::
   :maxdepth: 2

   api.rst
   usage.rst
   installation.rst
   faq.rst
   developer.rst
   release.rst
   benchmarking.rst
   glossary.rst

Indices and tables
------------------

Contents:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

References
----------

.. [Li.2009] *The Sequence Alignment/Map format and SAMtools.*
   Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup.
   Bioinformatics. 2009 Aug 15;25(16):2078-9.
   Epub 2009 Jun 8 `btp352 <https://doi.org/10.1093/bioinformatics/btp352>`_.
   PMID: `19505943 <https://pubmed.ncbi.nlm.nih.gov/19505943>`_.

.. [Bonfield.2021] *HTSlib: C library for reading/writing high-throughput sequencing data.*
   Bonfield JK, Marshall J, Danecek P, Li H, Ohan V, Whitwham A, Keane T, Davies RM.
   GigaScience (2021) 10(2) `giab007 <https://doi.org/10.1093/gigascience/giab007>`_.
   PMID: `33594436 <https://pubmed.ncbi.nlm.nih.gov/33594436>`_.

.. [Danecek.2021] *Twelve years of SAMtools and BCFtools.*
   Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H.
   GigaScience (2021) 10(2) `giab008 <https://doi.org/10.1093/gigascience/giab008>`_.
   PMID: `33590861 <https://pubmed.ncbi.nlm.nih.gov/33590861>`_.

.. seealso::
 
   Information about htslib
      http://www.htslib.org

   The samtools homepage
      http://samtools.sourceforge.net
   
   The cython C-extensions for python
      https://cython.org/

   The python language
      https://www.python.org
