.. _installation:

================
Installing pysam
================

Pysam can be installed through conda_, PyPI_ and from the repository.
The recommended way to install pysam is through conda/bioconda.

Conda installation
==================

To install pysam in your current conda_ environment, type::

   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda install pysam

This will install pysam from the bioconda_ channel and automatically
makes sure that dependencies are installed. Also, compilation flags
will be set automatically, which will potentially save a lot of
trouble on OS X.

PyPI installation
=================

Pysam provides a python interface to the functionality contained
within the htslib_ C library. There are two ways that these two
can be combined, ``builtin`` and ``external``.

Builtin
-------

The typical installation will be through PyPI_::

   pip install pysam

Generally you will have the ``wheel`` package installed and
this command will speedily install pysam from a pre-built wheel.
Otherwise, or if you use pip's ``--no-binary`` option, this will
compile the ``builtin`` htslib source code within pysam and allow
the configuration facilities described below to be used.

htslib_ can be configured at compilation to turn on additional
features such support using encrypted configurations, enable plugins,
and more. See the htslib_ project for more information on these.

Pysam will attempt to configure htslib_ to turn on some advanced
features. If these fail, for example due to missing library
dependencies (`libcurl`, `libcrypto`), it will fall back to
conservative defaults.

Options can be passed to the configure script explicitly by
setting the environment variable `HTSLIB_CONFIGURE_OPTIONS`.
For example::

  export HTSLIB_CONFIGURE_OPTIONS=--enable-plugins
  pip install pysam

External
--------

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

Note that generally the pysam and htslib version need to be
compatible. See the release notes for more information.

Installation from repository
============================

pysam depends on cython_ to provide the connectivity to the htslib_ C
library. The installation of the source tarball (:file:`.tar.gz`)
contains pre-built C-files and cython needs not be present
during installation. However, when installing from the repository,
cython needs to be installed beforehand.

To install from repository, type::

    python setup.py install

For compilation options, see the section on PyPI installation above.

Windows
=======

``pip install pysam`` installs a pre-built wheel on 64-bit x86 Windows
for Python 3.9 and later. The wheel needs no other software: htslib and
the libraries it uses are linked into pysam's modules. Remote files
(``https://``, ``s3://``, ``gs://``) are supported, and HTTPS
certificates are checked against the Windows certificate store.

Building from source on Windows
-------------------------------

pysam is built with MinGW-w64 from MSYS2_; Microsoft Visual C++ is not
supported. The build uses MSYS2's UCRT64 environment and a regular
Windows CPython, for example from python.org.

1. Install MSYS2_ and, in an MSYS2 shell, install the build
   prerequisites::

      pacman -S --needed make mingw-w64-ucrt-x86_64-gcc \
          mingw-w64-ucrt-x86_64-bzip2 mingw-w64-ucrt-x86_64-libsystre \
          mingw-w64-ucrt-x86_64-xz mingw-w64-ucrt-x86_64-zlib

   For remote file support, also install
   ``mingw-w64-ucrt-x86_64-curl-winssl`` and ``mingw-w64-ucrt-x86_64-pkgconf``.

2. Open a UCRT64 shell that keeps the Windows ``PATH``, so that your
   Windows Python can be found, by setting ``MSYS2_PATH_TYPE=inherit``
   before starting it (or with ``msys2_shell.cmd -ucrt64 -use-full-path``).

3. Build and install pysam as usual, using the Windows Python::

      pip install .

   or ``python setup.py build`` from a checkout. pysam selects the
   MinGW-w64 compiler on Windows automatically.

The resulting modules are statically linked, so they do not need MSYS2's
DLLs at run time.

When pysam is imported on Windows it ignores ``SIGTERM``, if the
application has not set a handler for it. htslib raises ``SIGTERM`` in
place of ``SIGPIPE`` when writing to a pipe whose reader has exited, and
ignoring it makes such writes fail with an error, as they do on other
platforms, rather than ending the process.

.. _MSYS2: https://www.msys2.org/

Requirements
============

Depending on the installation method, requirements for building pysam differ.

When installing through conda_, dependencies will be resolved by the
package manager. The pip_ installation and installation from source
require a C compiler and its standard libraries as well as all
requirements for building htslib. Htslib requirements are listed in
the htslib/INSTALL file.

Installing from the repository will require cython_ to be installed.

.. _pysam-profile:

Profiling
=========

Pysam's Cython code is no longer built with Python profiling enabled by default.
Profiling can be enabled when building it yourself from source by setting the
environment variable `PYSAM_PROFILE`, for example::

    export PYSAM_PROFILE=1
    pip install --no-binary pysam pysam

Note that when installing via pip you will need to use ``--no-binary pysam`` to
ensure pip does not simply install a pre-built wheel.
