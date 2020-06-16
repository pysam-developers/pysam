#!/bin/bash
#
# Build manylinux1 wheels for pysam. Based on the example at
# <https://github.com/pypa/python-manylinux-demo>
#
# It is best to run this in a fresh clone of the repository!
#
# Before running, make sure to update image:
#
#   docker pull quay.io/pypa/manylinux1_x86_64
#
# Run this within the repository root:
#   docker run --rm -v $(pwd):/io quay.io/pypa/manylinux1_x86_64 /io/buildwheels.sh
#
# The wheels will be put into the wheelhouse/ subdirectory.
#
# For interactive tests:
#   docker run -it -v $(pwd):/io quay.io/pypa/manylinux1_x86_64 /bin/bash

set -xeuo pipefail

# For convenience, if this script is called from outside of a docker container,
# it starts a container and runs itself inside of it.
if ! grep -q docker /proc/1/cgroup; then
  # We are not inside a container
  exec docker run --rm -v $(pwd):/io quay.io/pypa/manylinux1_x86_64 /io/$0
fi

yum install -y zlib-devel bzip2-devel xz-devel

# Without libcurl support, htslib can open files from HTTP and FTP URLs.
# With libcurl support, it also supports HTTPS and S3 URLs, but libcurl needs a
# current version of OpenSSL, and we do not want to be responsible for
# updating the wheels as soon as there are any security issues. So disable
# libcurl for now.
# See also <https://github.com/pypa/manylinux/issues/74>.
#
export HTSLIB_CONFIGURE_OPTIONS="--disable-libcurl"

PYBINS="/opt/python/*/bin"
echo $PYBINS

for PYBIN in ${PYBINS}; do
    ${PYBIN}/pip install -r /io/requirements.txt
    ${PYBIN}/pip wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
#
# The '-L ""' option is a workaround. By default, auditwheel puts all external
# libraries (.so files) into a .libs directory and sets the RUNPATH to $ORIGIN/.libs.
# When HTSLIB_MODE is 'shared' (now the default), then all so libraries part of
# pysam require that RUNPATH is set to $ORIGIN (without the .libs). It seems
# auditwheel overwrites $ORIGIN with $ORIGIN/.libs. This workaround makes
# auditwheel keeps the RUNPATH at "$ORIGIN" and put all the external libraries into
# the pysam directory.
for whl in wheelhouse/*.whl; do
    auditwheel repair -L "" $whl -w /io/wheelhouse/
done

# mkdir -p /io/wheelhouse
# cp wheelhouse/*.whl /io/wheelhouse

# Created files are owned by root, so fix permissions.
chown -R --reference=/io/setup.py /io/wheelhouse/

# TODO Install packages and test them
for PYBIN in ${PYBINS}; do
   ${PYBIN}/pip install pysam --no-index -f /io/wheelhouse
   # smoketest
   (cd $HOME; ${PYBIN}/python -c 'import pysam')
   # todo: add more tests
done
