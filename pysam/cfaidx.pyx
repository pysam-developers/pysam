# cython: embedsignature=True
# cython: profile=True
# adds doc-strings for sphinx
import sys
import os

cdef class FastqProxy
cdef makeFastqProxy(kseq_t * src):
    '''enter src into AlignedRead.'''
    cdef FastqProxy dest = FastqProxy.__new__(FastqProxy)
    dest._delegate = src
    return dest


from cpython cimport PyErr_SetString, \
    PyBytes_Check, \
    PyUnicode_Check, \
    PyBytes_FromStringAndSize

from cpython.version cimport PY_MAJOR_VERSION

from chtslib cimport \
    faidx_nseq, fai_load, fai_destroy, fai_fetch, \
    faidx_fetch_seq, gzopen, gzclose, \
    kseq_init, kseq_destroy, kseq_read


########################################################################
########################################################################
########################################################################
## Python 3 compatibility functions
########################################################################
IS_PYTHON3 = PY_MAJOR_VERSION >= 3

# filename encoding (copied from lxml.etree.pyx)
cdef str _FILENAME_ENCODING
_FILENAME_ENCODING = sys.getfilesystemencoding()
if _FILENAME_ENCODING is None:
    _FILENAME_ENCODING = sys.getdefaultencoding()
if _FILENAME_ENCODING is None:
    _FILENAME_ENCODING = 'ascii'

#cdef char* _C_FILENAME_ENCODING
#_C_FILENAME_ENCODING = <char*>_FILENAME_ENCODING

cdef bytes _encodeFilename(object filename):
    """Make sure a filename is 8-bit encoded (or None)."""
    if filename is None:
        return None
    elif PyBytes_Check(filename):
        return filename
    elif PyUnicode_Check(filename):
        return filename.encode(_FILENAME_ENCODING)
    else:
        raise TypeError, u"Argument must be string or unicode."



#####################################################################
# hard-coded constants
cdef int max_pos = 2 << 29

## TODO:
##        add automatic indexing.
##        add function to get sequence names.
cdef class FastaFile:
    '''*(filename)*

    A *FASTA* file. The file is automatically opened.

    This class expects an indexed fasta file and permits
    random access to fasta sequences.
    '''

    def __cinit__(self, *args, **kwargs ):
        self.fastafile = NULL
        self._filename = None
        self._references = None
        self._lengths = None
        self.reference2length = None
        self._open( *args, **kwargs )

    def _isOpen( self ):
        '''return true if samfile has been opened.'''
        return self.fastafile != NULL

    def __len__(self):
        if self.fastafile == NULL:
            raise ValueError("calling len() on closed file")

        return faidx_nseq(self.fastafile)

    def _open(self, filename):
        '''open an indexed fasta file.

        This method expects an indexed fasta file.
        '''

        # close a previously opened file
        if self.fastafile != NULL: self.close()
        self._filename = _encodeFilename(filename)
        self.fastafile = fai_load(self._filename)

        if self.fastafile == NULL:
            raise IOError("could not open file `%s`" % filename)

        # read index
        if not os.path.exists( self._filename + b".fai" ):
            raise ValueError("could not locate index file")

        with open( self._filename + b".fai" ) as inf:
            data = [ x.split("\t") for x in inf ]
            self._references = tuple(x[0] for x in data)
            self._lengths = tuple(int(x[1]) for x in data)
            self.reference2length = dict(zip(self._references, self._lengths))

    def close( self ):
        if self.fastafile != NULL:
            fai_destroy( self.fastafile )
            self.fastafile = NULL

    def __dealloc__(self):
        self.close()

    property filename:
        '''filename associated with this object.'''
        def __get__(self):
            return self._filename

    property references:
        '''tuple with the names of :term:`reference` sequences.'''
        def __get__(self):
            return self._references

    property nreferences:
        '''number of :term:`reference` sequences in the file.'''
        def __get__(self):
            return len(self._references) if self.references else None

    property lengths:
        '''tuple with the lengths of :term:`reference` sequences.'''
        def __get__(self):
            return self._lengths

    def fetch(self,
              reference=None,
              start=None,
              end=None,
              region=None):

        '''*(reference = None, start = None, end = None, region = None)*

        fetch sequences in a :term:`region` using 0-based indexing.

        The region is specified by :term:`reference`, *start* and *end*.

        fetch returns an empty string if the region is out of range or
        addresses an unknown *reference*.

        If *reference* is given and *start* is None, the sequence from the
        first base is returned. Similarly, if *end* is None, the sequence
        until the last base is returned.

        Alternatively, a samtools :term:`region` string can be supplied.
        '''

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        cdef int length
        cdef char * seq

        if not region:
            if reference is None:
                raise ValueError('no sequence/region supplied.')
            if start is None:
                start = 0
            if end is None:
                end = max_pos - 1

            if start > end:
                raise ValueError(
                    'invalid region: start (%i) > end (%i)' % (start, end))
            if start == end:
                return b""
            # valid ranges are from 0 to 2^29-1
            if not 0 <= start < max_pos:
                raise IndexError('start out of range (%i)' % start)
            if not 0 <= end < max_pos:
                raise IndexError('end out of range (%i)' % end)
            # note: faidx_fetch_seq has a bug such that out-of-range access
            # always returns the last residue. Hence do not use faidx_fetch_seq,
            # but use fai_fetch instead
            # seq = faidx_fetch_seq(self.fastafile,
            #                       reference,
            #                       start,
            #                       end-1,
            #                       &length)
            region = "%s:%i-%i" % (reference, start+1, end)
            if PY_MAJOR_VERSION >= 3:
                region = region.encode('ascii')
            seq = fai_fetch( self.fastafile,
                             region,
                             &length )
        else:
            # samtools adds a '\0' at the end
            seq = fai_fetch( self.fastafile, region, &length )

        # copy to python
        if seq == NULL:
            return b""
        else:
            try:
                py_seq = seq[:length]
            finally:
                free(seq)

        return py_seq

    cdef char * _fetch( self, char * reference, int start, int end, int * length ):
        '''fetch sequence for reference, start and end'''

        return faidx_fetch_seq(self.fastafile,
                               reference,
                               start,
                               end-1,
                               length )

    def get_reference_length(self, reference):
        '''return the length of reference.'''
        return self.reference2length[reference]

    def __getitem__(self, reference):
        return self.fetch(reference)

    def __contains__(self, reference):
        '''return true if reference in fasta file.'''
        return reference in self.reference2length


cdef class FastqProxy:
    def __init__(self): pass

    property name:
        def __get__(self):
            return self._delegate.name.s

    property sequence:
        def __get__(self):
            return self._delegate.seq.s

    property comment:
        def __get__(self):
            if self._delegate.comment.l:
                return self._delegate.comment.s
            else: return None

    property quality:
        def __get__(self):
            if self._delegate.qual.l:
                return self._delegate.qual.s
            else: return None


cdef class FastxFile:
    '''*(filename)*

    A :term:`fastq` or :term:`fasta` formatted file. The file
    is automatically opened.

    Entries in the file can be both fastq or fasta formatted
    or even a mixture of the two.

    This file object permits iterating over all entries in
    the file. Random access is not implemented. The iteration
    returns objects of type :class:`FastqProxy`

    '''
    def __cinit__(self, *args, **kwargs):
        # self.fastqfile = <gzFile*>NULL
        self._filename = None
        self.entry = NULL
        self._open(*args, **kwargs)

    def _isOpen( self ):
        '''return true if samfile has been opened.'''
        return self.entry != NULL

    def _open(self, filename):
        '''open a fastq/fasta file.
        '''
        self.close()

        if not os.path.exists(filename):
            raise IOError("no such file or directory: %s" % filename)

        filename = _encodeFilename(filename)
        self.fastqfile = gzopen(filename, "r")
        self.entry = kseq_init(self.fastqfile)
        self._filename = filename

    def close( self ):
        '''close file.'''
        if self.entry != NULL:
            gzclose(self.fastqfile)
            if self.entry:
                kseq_destroy(self.entry)
                self.entry = NULL

    def __dealloc__(self):
        self.close()

    property filename:
        '''filename associated with this object.'''
        def __get__(self):
            return self._filename

    def __iter__(self):
        if not self._isOpen():
            raise ValueError("I/O operation on closed file")
        return self

    cdef kseq_t * getCurrent(self):
        return self.entry

    cdef int cnext(self):
        '''C version of iterator
        '''
        return kseq_read(self.entry)

    def __next__(self):
        """
        python version of next().
        """
        cdef int l
        l = kseq_read(self.entry)
        if (l > 0):
            return makeFastqProxy(self.entry)
        else:
            raise StopIteration

# Compatibility Layer for pysam 0.8.1
cdef class FastqFile(FastxFile):
    pass

# Compatibility Layer for pysam < 0.8
cdef class Fastafile(FastaFile):
    pass

cdef class Fastqfile(FastxFile):
    pass

__all__ = ["FastaFile",
           "FastqFile",
           "Fastafile",
           "Fastqfile"]
