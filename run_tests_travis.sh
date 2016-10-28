#!/usr/bin/env bash

pushd .

WORKDIR=`pwd`

#Install miniconda python
if [ $TRAVIS_OS_NAME == "osx" ]; then
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O Miniconda3.sh
else
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3.sh
fi

bash Miniconda3.sh -b

# Create a new conda environment with the target python version
~/miniconda3/bin/conda install conda-build -y
~/miniconda3/bin/conda create -q -y --name testenv python=$CONDA_PY cython numpy nose psutil pip

# activate testenv environment
source ~/miniconda3/bin/activate testenv

# Need to make C compiler and linker use the anaconda includes and libraries:
export PREFIX=~/miniconda3/
export CFLAGS="-I${PREFIX}/include -L${PREFIX}/lib"
export HTSLIB_CONFIGURE_OPTIONS="--disable-libcurl"

# create a new folder to store external tools
mkdir -p $WORKDIR/external-tools

# install htslib
cd $WORKDIR/external-tools
curl --retry-delay 60 --retry 5 -L https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2 > htslib-1.3.1.tar.bz2
tar xjvf htslib-1.3.1.tar.bz2
cd htslib-1.3.1
make
PATH=$PATH:$WORKDIR/external-tools/htslib-1.3.1
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WORKDIR/external-tools/htslib-1.3.1

# install samtools, compile against htslib
cd $WORKDIR/external-tools
curl --retry-delay 60 --retry 5 -L http://downloads.sourceforge.net/project/samtools/samtools/1.3.1/samtools-1.3.1.tar.bz2 > samtools-1.3.1.tar.bz2

tar xjvf samtools-1.3.1.tar.bz2
cd samtools-1.3.1
./configure --with-htslib=../htslib-1.3.1
make
PATH=$PATH:$WORKDIR/external-tools/samtools-1.3.1

echo "installed samtools"
samtools --version

if [ $? != 0 ]; then
    exit 1
fi

# install bcftools
cd $WORKDIR/external-tools
curl --retry-delay 60 --retry 5 -L https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 > bcftools-1.3.1.tar.bz2

tar xjf bcftools-1.3.1.tar.bz2
cd bcftools-1.3.1
./configure --with-htslib=../htslib-1.3.1
make
PATH=$PATH:$WORKDIR/external-tools/bcftools-1.3.1

echo "installed bcftools"
bcftools --version

if [ $? != 0 ]; then
    exit 1
fi

popd

# Try building conda recipe first
~/miniconda3/bin/conda-build ci/conda-recipe/ --python=$CONDA_PY

# install code from the repository
python setup.py install

# find build/

# change into tests directory. Otherwise,
# 'import pysam' will import the repository,
# not the installed version. This causes
# problems in the compilation test.
cd tests

# create auxilliary data
echo
echo 'building test data'
echo
make -C pysam_data
make -C cbcf_data

# run nosetests
# -s: do not capture stdout, conflicts with pysam.dispatch
# -v: verbose output
nosetests -s -v

if [ $? != 0 ]; then
    exit 1
fi

# build source tar-ball
cd ..
python setup.py sdist

if [ $? != 0 ]; then
    exit 1
fi

# check for presence of config.h files
echo "checking for presence of config.h files in tar-ball"
tar -tvzf dist/pysam-*.tar.gz | grep "config.h$"

if [ $? != 1 ]; then
    exit 1
fi

# test pip installation from tar-ball with cython
echo "pip installing with cython"
pip install --verbose --no-deps --no-use-wheel dist/pysam-*.tar.gz

if [ $? != 0 ]; then
    exit 1
fi

# attempt pip installation without cython
echo "pip installing without cython"
~/miniconda3/bin/conda remove cython
~/miniconda3/bin/conda list
echo "python is" `which python`
pip install --verbose --no-deps --no-use-wheel --force-reinstall --upgrade dist/pysam-*.tar.gz

if [ $? != 0 ]; then
    exit 1
fi

# attempt pip installation without cython and without
# command line options
echo "pip installing without cython and no configure options"
export HTSLIB_CONFIGURE_OPTIONS=""
pip install --verbose --no-deps --no-use-wheel --force-reinstall --upgrade dist/pysam-*.tar.gz

if [ $? != 0 ]; then
    exit 1
fi
