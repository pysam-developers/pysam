#!/usr/bin/env bash

# function to detect the Operating System
detect_os(){

if [ -f /etc/os-release ]; then

   OS=$(cat /etc/os-release | awk '/VERSION_ID/ {sub("="," "); print $2;}' | sed 's/\"//g' | awk '{sub("\\."," "); print $1;}')
   if [ "$OS" != "12" ] ; then

      echo       
      echo " Ubuntu version not supported "
      echo
      echo " Only Ubuntu 12.x has been tested so far "
      echo 
      exit 1;

   fi

   OS="ubuntu"

elif [ -f /etc/system-release ]; then

   OS=$(cat /etc/system-release | awk ' {print $4;}' | awk '{sub("\\."," "); print $1;}')
   if [ "$OS" != "6" ] ; then
      echo
      echo " Scientific Linux version not supported "
      echo
      echo " Only 6.x Scientific Linux has been tested so far "
      echo
      exit 1;
   fi

   OS="sl"

else

   echo
   echo " Operating system not supported "
   echo
   echo " Exiting installation "
   echo
   exit 1;

fi
} # detect_os

# message to display when the OS is not correct
sanity_check_os() {
   echo
   echo " Unsupported operating system: $OS "
   echo " Installation aborted "
   echo
   exit 1;
} # sanity_check_os

# function to install operating system dependencies
install_os_packages() {

if [ "$OS" == "ubuntu" -o "$OS" == "travis" ] ; then

   echo
   echo " Installing packages for Ubuntu "
   echo

   apt-get install -y gcc g++

elif [ "$OS" == "sl" ] ; then

   echo 
   echo " Installing packages for Scientific Linux "
   echo

   yum -y install gcc zlib-devel gcc-c++

else

   sanity_check_os

fi # if-OS
} # install_os_packages

# function to install Python dependencies
install_python_deps() {

if [ "$OS" == "ubuntu" -o "$OS" == "sl" ] ; then

   echo
   echo " Installing Python dependencies for $1 "
   echo

   # Create virtual environment
   cd
   mkdir CGAT
   cd CGAT
   wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.10.1.tar.gz
   tar xvfz virtualenv-1.10.1.tar.gz
   rm virtualenv-1.10.1.tar.gz
   cd virtualenv-1.10.1
   python virtualenv.py cgat-venv
   source cgat-venv/bin/activate

   # Install Python prerequisites
   pip install cython

elif [ "$OS" == "travis" ] ; then
   # Travis-CI provides a virtualenv with Python 2.7
   echo 
   echo " Installing Python dependencies in travis "
   echo

   # Install Python prerequisites
   pip install cython
   pip install nose

else

   sanity_check_os

fi # if-OS
} # install_python_deps

# common set of tasks to prepare external dependencies
nosetests_external_deps() {
echo
echo " Running nosetests for $1 "
echo

pushd .

# create a new folder to store external tools
mkdir -p $HOME/CGAT/external-tools

# install samtools
cd $HOME/CGAT/external-tools
curl -L http://downloads.sourceforge.net/project/samtools/samtools/1.3/samtools-1.3.tar.bz2 > samtools-1.3.tar.bz2
tar xjf samtools-1.3.tar.bz2 
cd samtools-1.3
make
PATH=$PATH:$HOME/CGAT/external-tools/samtools-1.3

echo "installed samtools"
samtools --version

if [ $? != 0 ]; then
    exit 1
fi

# install bcftools
cd $HOME/CGAT/external-tools
curl -L https://github.com/samtools/bcftools/releases/download/1.3/bcftools-1.3.tar.bz2 > bcftools-1.3.tar.bz2
tar xjf bcftools-1.3.tar.bz2
cd bcftools-1.3
make
PATH=$PATH:$HOME/CGAT/external-tools/bcftools-1.3

echo "installed bcftools"
bcftools --version

if [ $? != 0 ]; then
    exit 1
fi

popd

} # nosetests_external_deps


# function to run nosetests
run_nosetests() {

echo
echo " Running nosetests for $1 "
echo

# prepare external dependencies
nosetests_external_deps $OS

# install code
python setup.py install

# change into tests directory. Otherwise,
# 'import pysam' will import the repository,
# not the installed version. This causes
# problems in the compilation test.
cd tests

# create auxiliary data
echo
echo 'building test data'
echo 
make -C pysam_data all
make -C cbcf_data all
make -C tabix_data all

# run nosetests
# -s: do not capture stdout, conflicts with pysam.dispatch
# -v: verbose output
nosetests -s -v 

} # run_nosetests

# function to display help message
help_message() {
echo
echo " Use this script as follows: "
echo
echo " 1) Become root and install the operating system* packages: "
echo " ./install-CGAT-tools.sh --install-os-packages"
echo
echo " 2) Now, as a normal user (non root), install the Python dependencies**: "
echo " ./install-CGAT-tools.sh --install-python-deps"
echo
echo " At this stage the CGAT Code Collection is ready to go and you do not need further steps. Please type the following for more information:"
echo " source $HOME/CGAT/virtualenv-1.10.1/cgat-venv/bin/activate"
echo " cgat --help "
echo
echo " The CGAT Code Collection tests the software with nosetests. If you are interested in running those, please continue with the following steps:"
echo
echo " 3) Become root to install external tools and set up the environment: "
echo " ./install-CGAT-tools.sh --install-nosetests-deps"
echo
echo " 4) Then, back as a normal user (non root), run nosetests as follows:"
echo " ./install-CGAT-tools.sh --run-nosetests"
echo 
echo " NOTES: "
echo " * Supported operating systems: Ubuntu 12.x and Scientific Linux 6.x "
echo " ** An isolated virtual environment will be created to install Python dependencies "
echo
exit 1;
}


# the main script starts here

if [ $# -eq 0 -o $# -gt 1 ] ; then

   help_message

else

   if [ "$1" == "--help" ] ; then

      help_message

   elif [ "$1" == "--travis" ] ; then

      OS="travis"
      install_os_packages
      install_python_deps
      run_nosetests

   elif [ "$1" == "--install-os-packages" ] ; then

      detect_os
      install_os_packages

   elif [ "$1" == "--install-python-deps" ] ; then

      detect_os
      install_python_deps

   elif [ "$1" == "--install-nosetests-deps" ] ; then

      detect_os
      install_nosetests_deps

   elif [ "$1" == "--run-nosetests" ] ; then

      detect_os
      run_nosetests

   else 

      echo 
      echo " Incorrect input parameter: $1 "
      help_message

   fi # if argument 1

fi # if number of input parameters

