#!/usr/bin/env bash

pushd .

WORKDIR=`pwd`

#Install miniconda python
if [ $TRAVIS_OS_NAME == "osx" ]; then
	curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
	bash Miniconda3-latest-MacOSX-x86_64.sh -b
else
	curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b
fi

# Create a new conda environment with the target python version
~/miniconda3/bin/conda install conda-build -y
~/miniconda3/bin/conda create -q -y --name testenv python=$CONDA_PY cython numpy nose

# Add new conda environment to PATH
export PATH=~/miniconda3/envs/testenv/bin/:$PATH

# Hack to force linking to anaconda libraries rather than system libraries
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/miniconda3/envs/testenv/lib/
#export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:~/miniconda3/envs/testenv/lib/

# Need to make C compiler and linker use the anaconda includes and libraries:
export PREFIX=~/miniconda3/
export CFLAGS="-I${PREFIX}/include -L${PREFIX}/lib"
export HTSLIB_CONFIGURE_OPTIONS="--disable-libcurl"

# create a new folder to store external tools
mkdir -p $WORKDIR/external-tools

# install htslib
cd $WORKDIR/external-tools
curl -L https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2 > htslib-1.3.1.tar.bz2
tar xjvf htslib-1.3.1.tar.bz2
cd htslib-1.3.1
make
PATH=$PATH:$WORKDIR/external-tools/htslib-1.3.1
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WORKDIR/external-tools/htslib-1.3.1

# install samtools, compile against htslib
cd $WORKDIR/external-tools
curl -L http://downloads.sourceforge.net/project/samtools/samtools/1.3.1/samtools-1.3.1.tar.bz2 > samtools-1.3.1.tar.bz2
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
curl -L https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 > bcftools-1.3.1.tar.bz2
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

# build source tar-ball and test installation from tar-ball
cd ..
python setup.py sdist
tar -xvzf dist/pysam-*.tar.gz
cd pysam-*
python setup.py install
