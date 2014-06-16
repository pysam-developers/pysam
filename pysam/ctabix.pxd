
cdef extern from "string.h":
  ctypedef int size_t
  void *memcpy(void *dst,void *src,size_t len)
  void *memmove(void *dst,void *src,size_t len)
  void *memset(void *b,int c,size_t len)
  char *strtok_r(char *str, char *delim, char **saveptr)
  char *strncpy(char *dest, char *src, size_t n)
  void *memchr(void *s, int c, size_t n)

cdef extern from "stdlib.h":
  void free(void *)
  void *malloc(size_t)
  void *calloc(size_t,size_t)
  void *realloc(void *,size_t)
  void qsort(void *base, size_t nmemb, size_t size,
             int (*compar)(void *,void *))
  int c_abs "abs" (int)
  int atoi( char *nptr)
  long atol( char *nptr)
  double atof( char *nptr)

cdef extern from "stdio.h":
  ctypedef struct FILE:
    pass
  FILE *fopen(char *,char *)
  FILE *freopen(char *path, char *mode, FILE *stream)
  int fileno(FILE *stream)
  int dup2(int oldfd, int newfd)
  int fflush(FILE *stream)

  FILE * stderr
  FILE * stdout
  int fclose(FILE *)
  int sscanf(char *str,char *fmt,...)
  int printf(char *str,char *fmt,...)
  int sprintf(char *str,char *fmt,...)
  int fprintf(FILE *ifile,char *fmt,...)
  char *fgets(char *str,int size,FILE *ifile)

cdef extern from "ctype.h":
  int toupper(int c)
  int tolower(int c)

cdef extern from "sys/types.h":
  pass

cdef extern from "sys/stat.h":
  pass

cdef extern from "fcntl.h":
  int open(char *pathname, int flags)
  
cdef extern from "unistd.h":
  ctypedef int ssize_t
  char *ttyname(int fd)
  int isatty(int fd)  
  ssize_t read(int fd, void *buf, size_t count)

cdef extern from "string.h":
  int strcmp(char *s1, char *s2)
  int strncmp(char *s1,char *s2,size_t len)
  char *strcpy(char *dest,char *src)
  char *strncpy(char *dest,char *src, size_t len)
  char *strdup(char *)
  char *strcat(char *,char *)
  size_t strlen(char *s)
  int memcmp( void * s1, void *s2, size_t len )

cdef extern from "stdint.h":
  ctypedef int int64_t
  ctypedef int int32_t
  ctypedef int uint32_t
  ctypedef int uint8_t
  ctypedef int uint64_t

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

cdef extern from "Python.h":
    ctypedef struct FILE
    char *fgets(char *str, int size, FILE *ifile)
    int feof(FILE *stream)
    size_t strlen(char *s)
    size_t getline(char **lineptr, size_t *n, FILE *stream)
    char *strstr(char *, char *)
    char *strchr(char *string, int c)
    int fileno(FILE *stream)
    FILE *fdopen(int fd, char *mode)

from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t

# tabix support
cdef extern from "htslib/tbx.h":
    
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

cdef extern from "pysam_stream.h":

    ctypedef struct kstring_t:
        size_t l
        size_t m
        char * s

    ctypedef struct kstream_t:
        pass

    kstream_t * ks_init( gzFile )

    int ks_read( kstream_t * )
    void ks_destroy( kstream_t * )
    int ks_getuntil( kstream_t *, int, kstring_t *, int * ) 

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
