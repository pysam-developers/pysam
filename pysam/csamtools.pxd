
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
  FILE *freopen(char *path, char *mode, FILE *stream)
  int fileno(FILE *stream)
  int dup2(int oldfd, int newfd)
  int fflush(FILE *stream)

  FILE * stderr
  FILE * stdout
  int fclose(FILE *)
  int sscanf(char *str,char *fmt,...)
  int sprintf(char *str,char *fmt,...)
  int fprintf(FILE *ifile,char *fmt,...)
  char *fgets(char *str,int size,FILE *ifile)

cdef extern from "unistd.h":
  char *ttyname(int fd)
  int isatty(int fd)  

cdef extern from "string.h":
  int strcmp(char *s1, char *s2)
  int strncmp(char *s1,char *s2,size_t len)
  char *strcpy(char *dest,char *src)
  char *strncpy(char *dest,char *src, size_t len)
  char *strdup(char *)
  char *strcat(char *,char *)
  size_t strlen(char *s)

cdef extern from "razf.h":
  pass

cdef extern from "bam.h":

  # IF _IOLIB=2, bamFile = BGZF, see bgzf.h
  # samtools uses KNETFILE, check how this works
  ctypedef int int64_t
  ctypedef int int32_t
  ctypedef int uint32_t
  ctypedef int uint8_t
  ctypedef int uint64_t

  ctypedef struct tamFile:
      pass

  ctypedef struct bamFile:
      pass

  ctypedef struct bam1_core_t:
      int32_t tid 
      int32_t pos
      uint32_t bin
      uint32_t qual
      uint32_t l_qname
      uint32_t flag
      uint32_t n_cigar
      int32_t l_qseq
      int32_t mtid 
      int32_t mpos 
      int32_t isize

  ctypedef struct bam1_t:
    bam1_core_t core
    int l_aux
    int data_len
    int m_data
    uint8_t *data

  ctypedef struct bam_pileup1_t:
      bam1_t *b 
      int32_t qpos 
      int indel
      int level
      uint32_t is_del
      uint32_t is_head
      uint32_t is_tail

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

  int64_t bam_seek( bamFile fp, uint64_t voffset, int where)

  int64_t bam_tell( bamFile fp )

  bam1_t * bam_dup1( bam1_t *src ) 
  
  bam1_t * bam_copy1(bam1_t *bdst, bam1_t *bsrc)
  bam_index_t *bam_index_load(char *f )

  void bam_index_destroy(bam_index_t *idx)

  void bam_destroy1( bam1_t * b) 

  int bam_parse_region(bam_header_t *header, char *str, int *ref_id, int *begin, int *end)

  bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data)

  ctypedef struct bam_fetch_iterator_t:
      pass
  
  bam_fetch_iterator_t* bam_init_fetch_iterator(bamFile fp, bam_index_t *idx, int tid, int beg, int end)
  
  bam1_t * bam_fetch_iterate(bam_fetch_iterator_t *iter)
  
  void bam_cleanup_fetch_iterator(bam_fetch_iterator_t *iter)

  bam_fetch_iterator_t* bam_init_fetchall_iterator(bamFile fp, bam_index_t *idx)

  bam1_t * bam_fetchall_iterate(bam_fetch_iterator_t *iter)

  int bam_fetch(bamFile fp, bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func)

  int bam_plbuf_push(bam1_t *b, bam_plbuf_t *buf)

  void bam_plbuf_destroy(bam_plbuf_t *buf)

  int bam_read1(bamFile fp, bam1_t *b)

  int bam_write1( bamFile fp, bam1_t *b)

  bam_header_t *bam_header_init()

  int bam_header_write( bamFile fp, bam_header_t *header)

  bam_header_t *bam_header_read( bamFile fp )

  void bam_header_destroy(bam_header_t *header)

  bam1_t * bam_dup1( bam1_t *src ) 
  
  bam1_t * bam_copy1(bam1_t *bdst, bam1_t *bsrc)

  uint8_t *bam_aux_get(bam1_t *b,  char tag[2])
  int bam_aux2i(uint8_t *s)
  float bam_aux2f(uint8_t *s)
  double bam_aux2d(uint8_t *s)
  char bam_aux2A( uint8_t *s)
  char *bam_aux2Z( uint8_t *s)
  
  void bam_init_header_hash(bam_header_t *header)

cdef extern from "sam.h":

  ctypedef struct samfile_t_un:
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

  int samread(samfile_t *fp, bam1_t *b)

  int samwrite(samfile_t *fp, bam1_t *b)

cdef extern from "faidx.h":

   ctypedef struct faidx_t:
      pass

   int fai_build(char *fn)

   void fai_destroy(faidx_t *fai)

   faidx_t *fai_load(char *fn)

   char *fai_fetch(faidx_t *fai, char *reg, int *len)

cdef extern from "pysam_util.h":

    ctypedef struct pair64_t:
        uint64_t u
        uint64_t v

    int pysam_bam_fetch_init(bamFile fp, bam_index_t *idx, int tid, int beg, int end, pair64_t ** offp )

    int pysam_bam_fetch_is_overlap(uint32_t beg, uint32_t end, bam1_t *b)

    int pysam_bam_plbuf_push(bam1_t *b, bam_plbuf_t *buf, int cont)

    int pysam_get_pos( bam_plbuf_t *buf)

    int pysam_get_tid( bam_plbuf_t *buf)

    bam_pileup1_t * pysam_get_pileup( bam_plbuf_t *buf)

    int pysam_dispatch(int argc, char *argv[] )


