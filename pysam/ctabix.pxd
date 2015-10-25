from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup
from libc.stdio cimport FILE, printf

# Note: this replaces python "open"!
cdef extern from "fcntl.h":
    int open(char *pathname, int flags)

cdef extern from "unistd.h" nogil:
    ctypedef int ssize_t
    ssize_t read(int fd, void *buf, size_t count)
    int close(int fd)

from pysam.chtslib cimport hts_idx_t, hts_itr_t, htsFile, \
    gzFile, tbx_t, kstring_t

# These functions are put here and not in chtslib.pxd in order
# to avoid warnings for unused functions.
cdef extern from "pysam_stream.h" nogil:

    ctypedef struct kstream_t:
        pass

    ctypedef struct kseq_t:
        kstring_t name
        kstring_t comment
        kstring_t seq
        kstring_t qual

    gzFile gzopen(char *, char *)
    kseq_t *kseq_init(gzFile)
    int kseq_read(kseq_t *)
    void kseq_destroy(kseq_t *)
    int gzclose(gzFile)

    kstream_t *ks_init(gzFile)
    void ks_destroy(kstream_t *)

    # Retrieve characters from stream until delimiter
    # is reached placing results in str.
    int ks_getuntil(kstream_t *,
                    int delimiter,
                    kstring_t * str,
                    int * dret)


cdef class tabix_file_iterator:
    cdef gzFile fh
    cdef kstream_t * kstream
    cdef kstring_t buffer
    cdef size_t size
    cdef Parser parser
    cdef int fd
    cdef int duplicated_fd
    cdef infile

    cdef __cnext__(self)

cdef class TabixFile:

    # pointer to tabixfile
    cdef htsFile * tabixfile
    # pointer to index structure
    cdef tbx_t * index

    # flag indicating whether file is remote
    cdef int isremote

    cdef _filename
    cdef _filename_index

    cdef Parser parser

    cdef encoding    

cdef class Parser:
    cdef encoding

    cdef parse(self, char * buffer, int len)

cdef class asTuple(Parser):
    cdef parse(self, char * buffer, int len)

cdef class asGTF(Parser):
    pass

cdef class asBed(Parser):
    pass

cdef class asVCF(Parser):
    pass

cdef class TabixIterator:
    cdef hts_itr_t * iterator
    cdef TabixFile tabixfile
    cdef kstring_t buffer
    cdef encoding
    cdef int __cnext__(self)

cdef class TabixIteratorParsed(TabixIterator):
    cdef Parser parser

cdef class GZIterator:
    cdef object _filename
    cdef gzFile gzipfile
    cdef kstream_t * kstream
    cdef kstring_t buffer
    cdef int __cnext__(self)
    cdef encoding

cdef class GZIteratorHead(GZIterator):
    pass

cdef class GZIteratorParsed(GZIterator):
    cdef Parser parser

# Compatibility Layer for pysam < 0.8
cdef class Tabixfile(TabixFile):
    pass

