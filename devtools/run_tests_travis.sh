#!/usr/bin/env bash

# test script for pysam.
# The script performs the following tasks:
# 1. Setup a conda environment and install dependencies via conda
# 2. Build pysam via the conda recipe
# 3. Build pysam via setup.py from repository
# 4. Run tests on the setup.py version
# 5. Additional build tests
# 5.1 pip install with cython
# 5.2 pip install without cython
# 5.3 pip install without cython and without configure options

pushd .

WORKDIR=`pwd`

#Install miniconda python
if [ $TRAVIS_OS_NAME == "osx" ]; then
	wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O Miniconda3.sh
else
	wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3.sh --no-check-certificate  # Default OS versions are old and have SSL / CERT issues
fi

bash Miniconda3.sh -b

# Create a new conda environment with the target python version
~/miniconda3/bin/conda install conda-build -y
~/miniconda3/bin/conda create -q -y --name testenv python=$CONDA_PY cython numpy pytest psutil pip

# activate testenv environment
source ~/miniconda3/bin/activate testenv

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# pin versions, so that tests do not fail when pysam/htslib out of step
# add htslib dependencies
# NB: force conda-forge:blas due to conda/conda#7548
conda install -y "samtools>=1.11" "bcftools>=1.11" "htslib>=1.11" xz curl bzip2 "conda-forge::blas=*=openblas"

# As HTSLIB_MODE is (defaulted to) 'shared', ensure we don't pick up
# the external headers from the Conda-installed htslib package.
mv $CONDA_PREFIX/include/htslib $CONDA_PREFIX/include/htslib-disable

export HTSLIB_CONFIGURE_OPTIONS="--disable-libcurl"

echo "show samtools, htslib, and bcftools versions"
samtools --version
htsfile --version
bcftools --version

# Try building conda recipe first
~/miniconda3/bin/conda-build devtools/conda-recipe/ --python=$CONDA_PY

# install code from the repository via setup.py
echo
echo "============ installing via setup.py from repository ============"
echo
python setup.py install || exit

# create auxiliary data
echo
echo 'building test data'
echo
make -C tests/pysam_data
make -C tests/cbcf_data

# echo any limits that are in place
ulimit -a

# run tests
pytest

if [ $? != 0 ]; then
    exit 1
fi

# build source tar-ball. Make sure to run 'build' target so that .pyx
# files are cythonized.
python setup.py build sdist

if [ $? != 0 ]; then
    exit 1
fi

# check for presence of config.h files
echo "checking for presence of config.h files in tar-ball"
tar -tvzf dist/pysam-*.tar.gz | grep "config.h$"

if [ $? != 1 ]; then
    echo "ERROR: found config.h in tar-ball"
    tar -tvzf dist/pysam-*.tar.gz | grep "config.h%"
    exit 1
fi

# test pip installation from tar-ball with cython
echo "pip installing with cython"
pip install --verbose --no-deps --no-binary=:all: dist/pysam-*.tar.gz

if [ $? != 0 ]; then
    exit 1
fi

# attempt pip installation without cython
echo "pip installing without cython"
~/miniconda3/bin/conda remove -y cython
~/miniconda3/bin/conda list
echo "python is" `which python`
pip install --verbose --no-deps --no-binary=:all: --force-reinstall --upgrade dist/pysam-*.tar.gz

if [ $? != 0 ]; then
    exit 1
fi

# attempt pip installation without cython and without
# command line options
echo "pip installing without cython and no configure options"
export HTSLIB_CONFIGURE_OPTIONS=""
pip install --verbose --no-deps --no-binary=:all: --force-reinstall --upgrade dist/pysam-*.tar.gz

if [ $? != 0 ]; then
    exit 1
fi
