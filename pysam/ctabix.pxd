
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

cdef extern from "Python.h":
    ctypedef struct FILE
    FILE* PyFile_AsFile(object)
    char *fgets(char *str, int size, FILE *ifile)
    int feof(FILE *stream)
    size_t strlen(char *s)
    size_t getline(char **lineptr, size_t *n, FILE *stream)
    char *strstr(char *, char *)
    char *strchr(char *string, int c)
    int fileno(FILE *stream)

cdef extern from "bgzf.h":

  ctypedef struct BGZF:
    pass

  int64_t bgzf_seek(BGZF* fp, int64_t pos, int where)

  BGZF * bgzf_open(char * path, char * mode)

  int bgzf_write(BGZF * fp, void* data, int length)

  int bgzf_close(BGZF* fp)

# tabix support
cdef extern from "tabix.h":

  ctypedef struct ti_conf_t:
    int32_t preset
    int32_t sc, bc, ec
    int32_t meta_char, line_skip

  ctypedef struct ti_index_t:
     pass
      
  ctypedef struct tabix_t: 
    BGZF *fp
    ti_index_t *idx
    char *fn
    char *fnidx

  ctypedef struct ti_iter_t:
    pass

  tabix_t *ti_open(char *fn, char *fnidx)

  int ti_lazy_index_load(tabix_t *t)

  void ti_close(tabix_t *t)

  ti_iter_t ti_query(tabix_t *t, char *name, int beg, int end)
  ti_iter_t ti_queryi(tabix_t *t, int tid, int beg, int end)
  ti_iter_t ti_querys(tabix_t *t, char *reg)
  char * ti_read(tabix_t *t, ti_iter_t iter, int *len)

  # Get the list of sequence names. Each "char*" pointer points to a
  #	internal member of the index, so DO NOT modify the returned
  #	 pointer; otherwise the index will be corrupted. The returned
  #	pointer should be freed by a single free() call by the routine
  #	calling this function. The number of sequences is returned at *n
  char **ti_seqname(ti_index_t *idx, int *n)
  
  # Destroy the iterator
  void ti_iter_destroy(ti_iter_t iter)

  # Build the index for file <fn>. File <fn>.tbi will be generated
  # and overwrite the file of the same name. Return -1 on failure. */
  int ti_index_build(char *fn, ti_conf_t *conf)

  #/* Load the index from file <fn>.tbi. If <fn> is a URL and the index
  #   * file is not in the working directory, <fn>.tbi will be
  #   * downloaded. Return NULL on failure. */
  ti_index_t *ti_index_load( char *fn)

  ti_index_t *ti_index_load_local(char *fnidx)

  #/* Destroy the index */
  void ti_index_destroy(ti_index_t *idx)

  #/* Parse a region like: chr2, chr2:100, chr2:100-200. Return -1 on failure. */
  int ti_parse_region( ti_index_t *idx,  char *str, int *tid, int *begin, int *end)

  int ti_get_tid( ti_index_t *idx,  char *name)

  #  /* Get the iterator pointing to the first record at the current file
  #   * position. If the file is just openned, the iterator points to the
  #   * first record in the file. */
  ti_iter_t ti_iter_first()

  #  /* Get the iterator pointing to the first record in region tid:beg-end */
  ti_iter_t ti_iter_query( ti_index_t *idx, int tid, int beg, int end)

  #  /* Get the data line pointed by the iterator and iterate to the next record. */
  # char *ti_iter_read(BGZF *fp, ti_iter_t iter, int *len)

cdef class Tabixfile:
    cdef char * _filename

    # pointer to tabixfile
    cdef tabix_t * tabixfile
     
cdef class Parser:
     pass
