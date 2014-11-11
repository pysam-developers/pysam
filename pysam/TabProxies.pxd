cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)       
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)
    int c_abs "abs" (int)  
    int c_abs "abs" (int)
    int atoi( char *nptr)
    long atol( char *nptr)
    double atof( char *nptr)

cdef extern from "Python.h":
    ctypedef struct FILE
    char *fgets(char *str, int size, FILE *ifile)
    int feof(FILE *stream)
    size_t strlen(char *s)
    size_t getline(char **lineptr, size_t *n, FILE *stream)
    char *strstr(char *, char *)
    char *strchr(char *string, int c)
    int fileno(FILE *stream)

cdef extern from "string.h":
  int strcmp(char *s1, char *s2)
  int strncmp(char *s1,char *s2,size_t len)
  char *strcpy(char *dest,char *src)
  char *strncpy(char *dest,char *src, size_t len)
  char *strdup(char *)
  char *strcat(char *,char *)
  size_t strlen(char *s)
  int memcmp( void * s1, void *s2, size_t len )
  void *memcpy(void *dest, void *src, size_t n)
  void *memchr(void *s, int c, size_t n)

cdef extern from "stdint.h":
  ctypedef int int64_t
  ctypedef int int32_t
  ctypedef int uint32_t
  ctypedef int uint8_t
  ctypedef int uint64_t

cdef class TupleProxy:

    cdef:
        char * data
        char ** fields
        int nfields
        int index
        int nbytes
        int offset
        bint is_modified

    cdef encoding

    cdef int getMaxFields(self)
    cdef int getMinFields(self)
#    cdef char * _getindex(self, int idx)

    cdef take(self, char * buffer, size_t nbytes)
    cdef present(self, char * buffer, size_t nbytes)
    cdef copy(self, char * buffer, size_t nbytes)
    cdef update(self, char * buffer, size_t nbytes)

cdef class GTFProxy(TupleProxy) :

    cdef:
        char * _attributes
        cdef bint hasOwnAttributes

    cdef int getMaxFields(self)
    cdef int getMinFields(self)
    cdef char * getAttributes( self )

cdef class NamedTupleProxy(TupleProxy) :
    pass

cdef class BedProxy(NamedTupleProxy) :

    cdef:
        char * contig
        uint32_t start
        uint32_t end
        int bedfields

    cdef int getMaxFields(self)
    cdef int getMinFields(self)
    cdef update(self, char * buffer, size_t nbytes)

cdef class VCFProxy(NamedTupleProxy) :

    cdef:
        char * contig
        uint32_t pos

    cdef update( self, char * buffer, size_t nbytes )
