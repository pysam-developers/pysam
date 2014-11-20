=================
Developer's guide
=================

Code organization
=================

The top level directory is organized in the following 
directories:

pysam
   Code specific to pysam

samtools
   Original and unmodified source code from :term:`csamtools`. Use 
   :file:`setup.py` to obtain the latest code.

tests
   Examples and data for testing

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

The following people have contributed to pysam:

* Andreas Heger
* Tildon Grant Belgrad
* Kevin Jacobs
* Florian Finkernagel
* Ben Schiller
* Marcel Martin
* Gerton Lunter
* Martin Goodson
* Leo Goodstadt








