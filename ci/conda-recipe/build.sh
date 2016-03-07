#!/bin/bash

# Use internal htslib
chmod a+x ./htslib/configure
export CFLAGS="-I${PREFIX}/include -L${PREFIX}/lib"
export HTSLIB_CONFIGURE_OPTIONS="--disable-libcurl"

# Use --old-and-unmanageable to disable egg-style folder naming, which is a minor hassle for RPATH stuff.
$PYTHON setup.py install --old-and-unmanageable

# Hack to find SO filenames, which are different on py35 and greater.
# We only manually adjust rpath on OSX, so we don't need to deal with linux sufix.

if [ $PY_VER == "3.5" ] && [ "$(uname)" == "Darwin" ]; then
  SO_SUFFIX=".cpython-35m-darwin.so"
else
  SO_SUFFIX=".so"
fi

echo SO_SUFFIX

# Hacky workaround to fix rpath pathing issues
if [ "$(uname)" == "Darwin" ]; then
  otool -L ${SP_DIR}/pysam/*.so
  install_name_tool -change @rpath/pysam/libchtslib.so @rpath/python${PY_VER}/site-packages/pysam/libchtslib${SO_SUFFIX} ${SP_DIR}/pysam/calignedsegment${SO_SUFFIX}
  install_name_tool -change @rpath/pysam/libchtslib.so @rpath/python${PY_VER}/site-packages/pysam/libchtslib${SO_SUFFIX} ${SP_DIR}/pysam/calignmentfile${SO_SUFFIX}
  install_name_tool -change @rpath/pysam/libchtslib.so @rpath/python${PY_VER}/site-packages/pysam/libchtslib${SO_SUFFIX} ${SP_DIR}/pysam/cbcf${SO_SUFFIX}
  install_name_tool -change @rpath/pysam/libchtslib.so @rpath/python${PY_VER}/site-packages/pysam/libchtslib${SO_SUFFIX} ${SP_DIR}/pysam/cfaidx${SO_SUFFIX}
  install_name_tool -change @rpath/pysam/libchtslib.so @rpath/python${PY_VER}/site-packages/pysam/libchtslib${SO_SUFFIX} ${SP_DIR}/pysam/csamfile${SO_SUFFIX}
  install_name_tool -change @rpath/pysam/libchtslib.so @rpath/python${PY_VER}/site-packages/pysam/libchtslib${SO_SUFFIX} ${SP_DIR}/pysam/ctabix${SO_SUFFIX}
  install_name_tool -change @rpath/pysam/libchtslib.so @rpath/python${PY_VER}/site-packages/pysam/libchtslib${SO_SUFFIX} ${SP_DIR}/pysam/ctabixproxies${SO_SUFFIX}
  install_name_tool -change @rpath/pysam/libchtslib.so @rpath/python${PY_VER}/site-packages/pysam/libchtslib${SO_SUFFIX} ${SP_DIR}/pysam/cutils${SO_SUFFIX}
  install_name_tool -change @rpath/pysam/libchtslib.so @rpath/python${PY_VER}/site-packages/pysam/libchtslib${SO_SUFFIX} ${SP_DIR}/pysam/cvcf${SO_SUFFIX}
  otool -L ${SP_DIR}/pysam/*.so
else
	echo "Skipping rpath workaround on linux"
fi
