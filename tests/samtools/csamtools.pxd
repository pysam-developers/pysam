
cdef extern from "string.h":
  ctypedef int size_t
  void *memcpy(void *dst,void *src,size_t len)
  void *memmove(void *dst,void *src,size_t len)
  void *memset(void *b,int c,size_t len)

cdef extern from "stdlib.h":
  void free(void *)
  void *malloc(size_t)
  void *calloc(size_t,size_t)
  void *realloc(void *,size_t)
  int c_abs "abs" (int)
  void qsort(void *base, size_t nmemb, size_t size,
             int (*compar)(void *,void *))

cdef extern from "stdio.h":
  ctypedef struct FILE:
    pass
  FILE *fopen(char *,char *)
  int fclose(FILE *)
  int sscanf(char *str,char *fmt,...)
  int sprintf(char *str,char *fmt,...)
  int fprintf(FILE *ifile,char *fmt,...)
  char *fgets(char *str,int size,FILE *ifile)

cdef extern from "string.h":
  int strcmp(char *s1, char *s2)
  int strncmp(char *s1,char *s2,size_t len)
  char *strcpy(char *dest,char *src)
  char *strdup(char *)
  char *strcat(char *,char *)

cdef extern from "bam.h":

  # IF _IOLIB=2, bamFile = BGZF, see bgzf.h
  # samtools uses KNETFILE, check how this works

  ctypedef int int64_t
  ctypedef int int32_t
  ctypedef int uint32_t

  ctypedef struct tamFile:
    pass

  ctypedef struct bamFile:
    pass

  ctypedef struct bam_pileup1_t:
    pass

  ctypedef struct bam1_t:
    pass

  ctypedef int (*bam_pileup_f)(uint32_t tid, uint32_t pos, int n, bam_pileup1_t *pl, void *data)

  ctypedef int (*bam_fetch_f)(bam1_t *b, void *data)

  ctypedef struct bam_header_t:
     int32_t n_targets
     char **target_name
     uint32_t *target_len
     void *hash
     void *rg2lib
     int l_text
     char *text

  ctypedef struct bam_index_t:
    pass

  ctypedef struct bam_plbuf_t:
    pass

  bamFile razf_dopen(int data_fd, char *mode)

  bam_index_t *bam_index_load(char *f )

  void bam_index_destroy(bam_index_t *idx)

  int bam_parse_region(bam_header_t *header, char *str, int *ref_id, int *begin, int *end)

  bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data)

  int bam_fetch(bamFile fp, bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func)

  int bam_plbuf_push(bam1_t *b, bam_plbuf_t *buf)

  void bam_plbuf_destroy(bam_plbuf_t *buf)


cdef extern from "sam.h":

  cdef struct samfile_t_un:
    tamFile tamr
    bamFile bam
    FILE *tamw
    
  ctypedef struct samfile_t:
     int type
     samfile_t_un x
     bam_header_t *header

  samfile_t *samopen( char *fn, char * mode, void *aux)

  int sampileup( samfile_t *fp, int mask, bam_pileup_f func, void *data)

  void samclose(samfile_t *fp)




