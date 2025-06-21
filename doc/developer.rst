=================
Developer's guide
=================

Code organization
=================

The top level directory is organized in the following 
directories:

:file:`pysam`
   Code specific to pysam.

:file:`doc`
   The documentation. To build the latest documentation, first install 
   `Sphinx`_ and then type::

       make -C doc html

:file:`tests`
   Code and data for testing and benchmarking.

:file:`htslib`
   Source code from `htslib`_ shipped with pysam. See
   :file:`import.py` about importing.

:file:`samtools`
   Source code from `samtools`_ shipped with pysam. See
   :file:`import.py` about importing.

:file:`bcftools`
   Source code from `bcftools`_ shipped with pysam. See
   :file:`import.py` about importing.


Python language level
=====================

Pysam currently requires Python 3.8 as a minimum language level.
For example, this means that the following comparatively recent
language features and library functions are available for use:

* f-strings
* ``raise ... from None``
* :meth:`str.startswith`, :meth:`str.endswith`, :meth:`str.rstrip`, etc
* walrus ``:=`` operator in Python code

However in particular the following should not be used in
pysam source code or infrastructure scripts:

* :meth:`str.removeprefix`, :meth:`str.removesuffix` (new in 3.9)
* walrus ``:=`` operator in Cython code (requires Cython 3)
* ``Optional[type]`` type hints written as ``type | None`` etc (new in 3.10)


Importing new versions of htslib and samtools
=============================================

See instructions in :file:`import.py` to import the latest
versions of `htslib`_, `samtools`_ and `bcftools`_.

Unit testing
============

Unit tests are in the :file:`tests` directory. To run all unit tests,
run::

   pytest tests

Most tests use test data from the :file:`tests/*_data` directories.
Some of these test data files are generated from other files in these
directories, which is done by running ``make`` in each directory::

   make -C tests/pysam_data
   # etc

Alternatively if any :file:`tests/*_data/all.stamp` file is not already
present, running the unit tests should generate that directory's data
files automatically.

Benchmarking
============

To run the benchmarking suite, make sure that `pytest-benchmark
<https://github.com/ionelmc/pytest-benchmark>`_ is installed. To run
all benchmarks, type::

   pytest tests/*_bench.py

See :ref:`Benchmarking` for more on this topic.

Contributors
============

Please see Github for a list of all contributors:

https://github.com/pysam-developers/pysam/graphs/contributors

Many thanks to all contributors for helping in making pysam
useful.
