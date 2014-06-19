from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup
from libc.stdio cimport FILE, printf
                     
cdef extern from "Python.h":
   long _Py_HashPointer(void*)
   FILE* PyFile_AsFile(object)

cdef extern from "zlib.h":
  ctypedef void * gzFile
  ctypedef int64_t z_off_t

  int gzclose(gzFile fp)
  int gzread(gzFile fp, void *buf, unsigned int n)
  char *gzerror(gzFile fp, int *errnum)

  gzFile gzopen( char *path, char *mode)
  gzFile gzdopen (int fd, char *mode)
  char * gzgets(gzFile file, char *buf, int len)
  int gzeof( gzFile file )

cdef extern from "htslib/kstring.h" nogil:
    ctypedef struct kstring_t:
        size_t l, m
        char *s

cdef extern from "pysam_stream.h" nogil:

    ctypedef struct kstream_t:
        pass

    ctypedef struct kseq_t:
        kstring_t name
        kstring_t comment
        kstring_t seq
        kstring_t qual

    gzFile gzopen(char *, char *)
    kseq_t * kseq_init(gzFile)
    int kseq_read(kseq_t *)
    void kseq_destroy(kseq_t *)
    int gzclose(gzFile)

    kstream_t * ks_init(gzFile)

    # Retrieve characters from stream until delimiter
    # is reached placing results in str.
    int ks_getuntil(kstream_t *, int delimiter,
                    kstring_t str, int * dret)

cdef extern from "htslib/faidx.h":

   ctypedef struct faidx_t:
      pass

   int fai_build(char *fn)

   void fai_destroy(faidx_t *fai)

   faidx_t *fai_load(char *fn)

   char *fai_fetch(faidx_t *fai,
                   char *reg,
                   int *len)

   int faidx_fetch_nseq(faidx_t *fai)

   char *faidx_fetch_seq(faidx_t *fai,
                         char *c_name,
                         int p_beg_i,
                         int p_end_i,
                         int *len)


# Exposing pysam extension classes
#
# Note: need to declare all C fields and methods here
cdef class Fastafile:
    cdef object _filename, _references, _lengths, reference2length
    cdef faidx_t* fastafile
    cdef char* _fetch(self, char* reference, int start, int end, int* length)


cdef class FastqProxy:
    cdef kseq_t * _delegate


cdef class Fastqfile:
    cdef object _filename
    cdef gzFile fastqfile
    cdef kseq_t * entry

    cdef kseq_t * getCurrent( self )
    cdef int cnext(self)

