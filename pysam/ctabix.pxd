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

cdef extern from "zlib.h" nogil:
  ctypedef void * gzFile
  ctypedef int64_t z_off_t

  int gzclose(gzFile fp)
  int gzread(gzFile fp, void *buf, unsigned int n)
  char *gzerror(gzFile fp, int *errnum)

  gzFile gzopen( char *path, char *mode)
  gzFile gzdopen (int fd, char *mode)
  char * gzgets(gzFile file, char *buf, int len)
  int gzeof( gzFile file )

cdef extern from "pysam_stream.h" nogil:

    ctypedef struct kstream_t:
        pass

    ctypedef struct kstring_t:
        size_t l
        size_t m
        char * s

    kstream_t * ks_init(gzFile)

    int ks_read(kstream_t *)
    void ks_destroy(kstream_t *)
    int ks_getuntil(kstream_t *, int, kstring_t *, int *) 

# tabix support
cdef extern from "htslib/tbx.h" nogil:
    
    # redefinitions from chtslib.pxd 
    ctypedef struct hts_idx_t
    ctypedef struct hts_itr_t
    ctypedef struct htsFile

    # tbx.h definitions
    int8_t TBX_MAX_SHIFT
    int8_t TBX_GENERIC
    int8_t TBX_SAM
    int8_t TBX_VCF
    int8_t TBX_UCSC

    ctypedef struct tbx_conf_t:
        int32_t preset
        int32_t sc, bc, ec   # seq col., beg col. and end col.
        int32_t meta_char, line_skip

    ctypedef struct tbx_t:
        tbx_conf_t conf
        hts_idx_t *idx
        void * dict

    tbx_conf_t tbx_conf_gff
    tbx_conf_t tbx_conf_bed
    tbx_conf_t tbx_conf_psltbl
    tbx_conf_t tbx_conf_sam
    tbx_conf_t tbx_conf_vcf
    
    void tbx_itr_destroy(hts_itr_t * iter)
    hts_itr_t * tbx_itr_queryi(tbx_t * t, int tid, int bed, int end)
    hts_itr_t * tbx_itr_querys(tbx_t * t, char * s)
    int tbx_itr_next(htsFile * fp, tbx_t * t, hts_itr_t * iter, void * data)

    int tbx_name2id(tbx_t *tbx, char *ss)

    int tbx_index_build(char *fn,
                        int min_shift,
                        tbx_conf_t *conf)
    
    tbx_t * tbx_index_load(char *fn)

    # free the array but not the values
    char **tbx_seqnames(tbx_t *tbx, int *n)

    void tbx_destroy(tbx_t *tbx)


cdef class tabix_file_iterator:
    cdef gzFile fh
    cdef kstream_t * kstream
    cdef kstring_t buffer
    cdef size_t size
    cdef Parser parser
    cdef int fd
    cdef infile

    cdef __cnext__(self)

cdef class Tabixfile:

    # pointer to tabixfile
    cdef htsFile * tabixfile
    # pointer to index structure
    cdef tbx_t * index

    # flag indicating whether file is remote
    cdef int isremote

    cdef char * _filename

    cdef Parser parser
    
cdef class Parser:
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
    cdef Tabixfile tabixfile
    cdef kstring_t buffer
    cdef int __cnext__(self)

cdef class TabixIteratorParsed(TabixIterator):
    cdef Parser parser

cdef class GZIterator:
    cdef object _filename
    cdef gzFile gzipfile
    cdef kstream_t * kstream
    cdef kstring_t buffer
    cdef int __cnext__(self)

cdef class GZIteratorHead(GZIterator):
    pass

cdef class GZIteratorParsed(GZIterator):
    cdef Parser parser

cdef _force_str(object s)
