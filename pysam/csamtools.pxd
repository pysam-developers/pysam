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

cdef extern from "math.h":
   double sqrt(double x)

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
  int printf(char *fmt,...)
  int sprintf(char *str,char *fmt,...)
  int fprintf(FILE *ifile,char *fmt,...)
  char *fgets(char *str,int size,FILE *ifile)

cdef extern from "ctype.h":
  int toupper(int c)
  int tolower(int c)
  
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
  int memcmp( void * s1, void *s2, size_t len )

cdef extern from "Python.h":
   long _Py_HashPointer(void*)
   FILE* PyFile_AsFile(object)

cdef extern from "fileobject.h":
   ctypedef class __builtin__.file [object PyFileObject]:
        pass

cdef extern from "razf.h":
  pass

cdef extern from "stdint.h":
  ctypedef int int64_t
  ctypedef int int32_t
  ctypedef int uint32_t
  ctypedef int uint8_t
  ctypedef int uint64_t

cdef extern from "bam.h":

  # constants
  int BAM_DEF_MASK
  # IF _IOLIB=2, bamFile = BGZF, see bgzf.h
  # samtools uses KNETFILE, check how this works

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
      int32_t n
      uint64_t n_no_coor

  ctypedef struct bam_plbuf_t:
      pass

  ctypedef struct pair64_t:
      uint64_t u, v
      
  ctypedef struct bam_iter_t:
      int from_first
      int tid, beg, end, n_off, i, finished
      uint64_t curr_off
      pair64_t *off

  # ctypedef __bam_iter_t * bam_iter_t

  bam1_t * bam_init1()
  void bam_destroy1(bam1_t *)

  bamFile razf_dopen(int data_fd, char *mode)

  int64_t bam_seek( bamFile fp, uint64_t voffset, int where)
  int64_t bam_tell( bamFile fp )

  # void bam_init_header_hash(bam_header_t *header)

  ###############################################
  # stand-ins for samtools macros
  uint32_t * bam1_cigar( bam1_t * b)
  char * bam1_qname( bam1_t * b)
  uint8_t * bam1_seq( bam1_t * b)
  uint8_t * bam1_qual( bam1_t * b)
  uint8_t * bam1_aux( bam1_t * b)

  ###############################################
  # bam iterator interface
  bam_iter_t bam_iter_query( bam_index_t *idx, int tid, int beg, int end)

  int bam_iter_read(bamFile fp, bam_iter_t iter, bam1_t *b)

  void bam_iter_destroy(bam_iter_t iter)

  ###############################################

  bam1_t * bam_dup1( bam1_t *src ) 
  
  bam1_t * bam_copy1(bam1_t *bdst, bam1_t *bsrc)
  bam_index_t *bam_index_load(char *f )

  void bam_index_destroy(bam_index_t *idx)

  int bam_parse_region(bam_header_t *header, char *str, int *ref_id, int *begin, int *end)

  ###############################################
  bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data)

  int bam_fetch(bamFile fp, bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func)

  int bam_plbuf_push(bam1_t *b, bam_plbuf_t *buf)

  void bam_plbuf_destroy(bam_plbuf_t *buf)
  ########################################
  # pileup iterator interface
  ctypedef struct bam_plp_t:
      pass

  ctypedef bam_pileup1_t * const_bam_pileup1_t_ptr "const bam_pileup1_t *"

  ctypedef int (*bam_plp_auto_f)(void *data, bam1_t *b)

  bam_plp_t bam_plp_init( bam_plp_auto_f func, void *data)
  int bam_plp_push( bam_plp_t iter,  bam1_t *b)
  bam_pileup1_t * bam_plp_next( bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
  bam_pileup1_t * bam_plp_auto( bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
  void bam_plp_set_mask(bam_plp_t iter, int mask)
  void bam_plp_reset(bam_plp_t iter)
  void bam_plp_destroy(bam_plp_t iter)
  void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt)

  ##################################################

  int bam_read1( bamFile fp, bam1_t *b)
  int bam_validate1( bam_header_t *header, bam1_t *b)
  int bam_write1( bamFile fp, bam1_t *b)

  bam_header_t *bam_header_init()

  int bam_header_write( bamFile fp, bam_header_t *header)

  bam_header_t *bam_header_read( bamFile fp )

  void bam_header_destroy(bam_header_t *header)

  bam1_t * bam_dup1( bam1_t *src ) 
  
  bam1_t * bam_copy1(bam1_t *bdst, bam1_t *bsrc)

  uint8_t *bam_aux_get(bam1_t *b,  char tag[2])

  int32_t bam_aux2i(uint8_t *s)
  float bam_aux2f(uint8_t *s)
  double bam_aux2d(uint8_t *s)
  char bam_aux2A( uint8_t *s)
  char *bam_aux2Z( uint8_t *s)
  
  int bam_reg2bin(uint32_t beg, uint32_t end)

  uint32_t bam_calend(bam1_core_t *c, uint32_t *cigar)

cdef extern from *:
    ctypedef char* const_char_ptr "const char*"

cdef extern from "sam.h":

  ctypedef struct samfile_t_un:
    tamFile tamr
    bamFile bam
    FILE *tamw
    
  ctypedef struct samfile_t:
     int type
     samfile_t_un x
     bam_header_t *header

  samfile_t *samopen( const_char_ptr fn, char * mode, void *aux)

  int sampileup( samfile_t *fp, int mask, bam_pileup_f func, void *data)

  void samclose(samfile_t *fp)

  int samread(samfile_t *fp, bam1_t *b)

  int samwrite(samfile_t *fp, bam1_t *b)

  int bam_prob_realn(bam1_t *b, char *ref)
  int bam_cap_mapQ(bam1_t *b, char *ref, int thres)


#cdef extern from "glf.h":
#   ctypedef struct glf1_t:
#      pass

#cdef extern from "bam_maqcns.h":
#
#  ctypedef struct bam_maqcns_t:
#     float het_rate, theta
#     int n_hap, cap_mapQ, errmod, min_baseQ
#     float eta, q_r
#     double *fk, *coef
#     double *lhet
#     void *aux

#  glf1_t *bam_maqcns_glfgen(int n, 
#                            bam_pileup1_t *pl, 
#                            uint8_t ref_base, 
#                            bam_maqcns_t *bm)

#  ctypedef struct bam_maqindel_opt_t:
#      int q_indel
#      float r_indel
#      float r_snp
#      int mm_penalty, indel_err, ambi_thres
     
#  uint32_t bam_maqcns_call(int n, bam_pileup1_t *pl, bam_maqcns_t *bm)
#  bam_maqcns_t * bam_maqcns_init()
#  void bam_maqcns_destroy(bam_maqcns_t *bm)
#  void bam_maqcns_prepare(bam_maqcns_t *bm)
  
#  uint32_t glf2cns(glf1_t *g, int q_r)

#  int BAM_ERRMOD_MAQ2
#  int BAM_ERRMOD_MAQ
#  int BAM_ERRMOD_SOAP

#  ctypedef struct bam_maqindel_ret_t: 
#    int indel1
#    int indel2        
#    int cnt1
#    int cnt2
#    int cnt_anti
#    int cnt_ref
#    int cnt_ambi
#    char *s[2]
#    int gt
#    int gl[2]
#    int q_cns
#    int q_ref
    
#  void bam_maqindel_ret_destroy( bam_maqindel_ret_t * )

#  bam_maqindel_opt_t *bam_maqindel_opt_init()

#  bam_maqindel_ret_t * bam_maqindel(int n, 
#  		     int pos, 
#  		     bam_maqindel_opt_t * mi, 
#  		     bam_pileup1_t * pl, 
#		     char *ref,
#		     int _n_types, 
#		     int * _types )
                                                               

cdef extern from "faidx.h":

   ctypedef struct faidx_t:
      pass

   int fai_build(char *fn)

   void fai_destroy(faidx_t *fai)

   faidx_t *fai_load(char *fn)

   char *fai_fetch(faidx_t *fai, char *reg, int *len)

   int faidx_fetch_nseq(faidx_t *fai)

   char *faidx_fetch_seq(faidx_t *fai, char *c_name, 
                         int p_beg_i, int p_end_i, int *len)


cdef extern from "pysam_util.h":

    int pysam_pileup_next(bam1_t *b, 
                          bam_plbuf_t *buf, 
                          bam_pileup1_t ** plp,
                          int * tid,
                          int * pos,
                          int * n_plp )


    int pysam_dispatch(int argc, char *argv[] )

    # stand-in functions for samtools macros
    void pysam_bam_destroy1( bam1_t * b) 

    # add *nbytes* into the variable length data of *src* at *pos*
    bam1_t * pysam_bam_update( bam1_t * b, 
                               size_t nbytes_old,
                               size_t nbytes_new,
                               uint8_t * pos )

    # translate char to unsigned char
    unsigned char pysam_translate_sequence( char s )

    unsigned char * bam_nt16_table

    int pysam_reference2tid( bam_header_t *header, char * s )

    void pysam_set_stderr( FILE * file )

    # return mapped/unmapped reads on tid
    uint32_t pysam_get_mapped( bam_index_t *idx, int tid )
    uint32_t pysam_get_unmapped( bam_index_t *idx, int tid )

#    uint32_t pysam_glf_depth( glf1_t * g )

#    void pysam_dump_glf( glf1_t * g, bam_maqcns_t * c )

# need to declare all C fields and methods here
cdef class AlignedRead:

    # object that this AlignedRead represents
    cdef bam1_t * _delegate

cdef class Samfile:
    cdef char * _filename
    # pointer to samfile
    cdef samfile_t * samfile
    # pointer to index
    cdef bam_index_t *index
    # true if file is a bam file
    cdef int isbam
    # true if not a file but a stream
    cdef int isstream
    # true if file is not on the local filesystem
    cdef int isremote
    # current read within iteration
    cdef bam1_t * b
    # file opening mode
    cdef char * mode

    # beginning of read section
    cdef int64_t start_offset 

    cdef bam_header_t * _buildHeader( self, new_header )
    cdef bam1_t * getCurrent( self )
    cdef int cnext(self)

    # write an aligned read
    cpdef int write( self, AlignedRead read )

    cdef char * _getrname( self, int tid )

cdef class IteratorRow:
    pass

cdef class IteratorRowAll(IteratorRow):
    cdef bam1_t * b
    cdef samfile_t * fp
    # true if samfile belongs to this object
    cdef int owns_samfile

    cdef bam1_t * getCurrent( self )

    cdef int cnext(self)


