*****************
Developer's guide
*****************

Code organization
*****************

The top level directory is organized in the following 
directories:

pysam
   Code specific to pysam

samtools
   Original and unmodified source code from :term:`csamtools`. Use 
   :file:`setup.py` to obtain the latest code.

tests
   Examples and data for testing

Importing :term:`csamtools`
***************************

Running :file:`setup.py` will import the csamtools source code. 
The command::

   python setup.py import PATH

where ``PATH`` points to a :term:`csamtools` source directory. For example::

   python setup.py import ~/samtools-0.1.6

Note that files will not be overwritten. To import all anew, 
delete all :file:`*.c` and :file:`*.h` files in the :file:`samtools`
directory first. 

Unit testing
************

Unit tests are currently in the script :file:`pysam_test.py`. 

Unit tests for the python methods
---------------------------------

Few implemented yet. This will require also writing functionality,
hence postponed

Unit tests for the command line interface
-----------------------------------------

The current suite of tests compare the binary files of selected
samtools commands against running the same commands from within
the pysam module. The general expectation is that the files
are binary identical given that most of the code is
from :term:`csamtools` anyway. However, differences might be
found if the installed :term:`csamtools` version is different
from the one wrapped with pysam.

The tests create files in the current test directory. They
are modeled on the example given within the :term:`csamtools`
distribution. Two files are required in the working directory
of the test script:

1. :file:`ex1.fa`: a fasta file
2. :file:`ex1.sam.gz`: a sam file










