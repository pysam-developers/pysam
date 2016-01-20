#!/usr/bin/env bash

pushd .

# create a new folder to store external tools
mkdir -p $HOME/external-tools
cd $HOME/external-tools

# install samtools
curl -L http://downloads.sourceforge.net/project/samtools/samtools/1.2/samtools-1.2.tar.bz2 > samtools-1.2.tar.bz2
tar xjvf samtools-1.2.tar.bz2 
cd samtools-1.2
make
PATH=$PATH:$HOME/external-tools/samtools-1.2

popd

# install code
python setup.py install

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

