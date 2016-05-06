.. _installation:

================
Installing pysam
================

Pysam provides a python interface to the functionality contained
within the htslib_ C library. There are two ways that these two
can be combined, ``builtin`` and ``external``.

Builtin
=======

The typical installation will be through pypi_::

   pip install pysam

This will compile the ``builtin`` htslib source code within pysam.

htslib_ can be configured at compilation to turn on additional
features such support using encrypted configurations, enable plugins,
and more. See the htslib_ project for more information on these.

Pysam will attempt to configure htslib_ to turn on some advanced
features. If these fail, for example due to missing library
dependencies (`libcurl`, `libcrypto`), it will fall back to
conservative defaults.

Options can be passed to the configure script explicitely by
setting the environment variable `HTSLIB_CONFIGURE_OPTIONS`.
For example::

  export HTSLIB_CONFIGURE_OPTIONS=--enable-plugins
  pip install pysam

External
========

pysam can be combined with an externally installed htslib_
library. This is a good way to avoid duplication of libraries. To link
against an externally installed library, set the environment variables
`HTSLIB_LIBRARY_DIR` and `HTSLIB_INCLUDE_DIR` before installing::

   export HTSLIB_LIBRARY_DIR=/usr/local/lib
   export HTSLIB_INCLUDE_DIR=/usr/local/include
   pip install pysam

Note that the location of the file :file:`libhts.so` needs to be known
to the linker once you run pysam, for example by setting the
environment-varirable `LD_LIBRARY_PATH`.

cython
======

pysam depends on cython_ to provide the connectivity to the htslib_ C
library. The installation of the source tarball (:file:`.tar.gz`)
python 2.7 contains pre-built C-files and cython needs not be present
during installation. However, when installing the source tarball on
python 3 or building from the repository, these pre-built C-files are
not present and cython needs to be installed beforehand.
