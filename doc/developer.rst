=================
Developer's guide
=================

Code organization
=================

The top level directory is organized in the following 
directories:

:file:`pysam`
   Code specific to pysam

:file:`doc`
   The documentation. To build the latest documention type::

       make -C doc html

:file:`tests`
   Code and data for testing

:file:`htslib`
   Source code from htslib_ shipped with pysam. See
   :file:`setup.py` about importing.

:file:`samtools`
   Source code from :term:`csamtools` shipped with pysam. See
   :file:`setup.py` about importing.


Importing new versions of htslib and samtools
=============================================

See instructions in :file:`setup.py` to import the latest
version of htslib_ and samtools_.

Unit testing
============

Unit tests are in the :file:`tests` directory. To run all unit tests,
run::

   nosetests -s -v tests

Note to use the ``-s/--nocapture`` option to prevent nosetests from
captpuring standard output.

Contributors
============

Please see github for a list of all contributors:

https://github.com/pysam-developers/pysam/graphs/contributors

Many thanks to all contributors for helping in making pysam
useful.






