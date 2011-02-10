.. samtools documentation master file, created by
   sphinx-quickstart on Wed Aug 12 14:43:42 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pysam: samtools interface for python
====================================

:Author: Andreas Heger, Tildon Grant Belgrad, Martin Goodson, Leo Goodstad, Kevin Jacobs
:Date: |today|
:Version: |version|

The *SAM/BAM* format is a way to store efficiently large numbers of alignments
[Li2009]_, such as those routinely are created by next-generation sequencing 
methods.

This module provides a low-level wrapper around the samtools_ C-API using `cython`_
and a high-level API for convenient access to the data in *SAM/BAM* formatted files.
Also included is an interface to the tabix_ C-API for reading compressed and
indexed tabular data.

The current version wraps *samtools-0.1.10a* and *tabix-0.2.3*.

Contents
--------

.. toctree::
   :maxdepth: 2

   api.rst
   usage.rst
   glossary.rst
   developer.rst
   release.rst

Indices and tables
------------------

Contents:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

References
----------

.. [Li2009] The Sequence Alignment/Map format and SAMtools. Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup.
   	    Bioinformatics. 2009 Aug 15;25(16):2078-9. Epub 2009 Jun 8.
	    `PMID: 19505943 <http://www.ncbi.nlm.nih.gov/pubmed/19505943?dopt=Abstract>`_

.. seealso::

   The samtools homepage
      http://samtools.sourceforge.net
   
   The cython C-extensions for python
      http://cython.org/

   The python language
      http://www.python.org

.. _samtools: http://samtools.sourceforge.net

.. _tabix: http://samtools.sourceforge.net/tabix.shtml

.. _cython: http://cython.org/

.. _python: http://www.python.org


