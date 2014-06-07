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

cdef extern from "stdint.h":
  ctypedef int int8_t
  ctypedef int int16_t
  ctypedef int int32_t
  ctypedef int int64_t
  ctypedef int uint8_t
  ctypedef int uint16_t
  ctypedef int uint32_t
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


cdef extern from "pysam_stream.h":

    ctypedef struct kstring_t:
      size_t l
      size_t m
      char *s

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


cdef extern from "hts.h":

  ctypedef struct BGZF:
      pass
  ctypedef struct cram_fd:
      pass
  ctypedef struct hFILE:
      pass

  unsigned char * seq_nt16_table

  ctypedef int hts_readrec_func(BGZF *fp,
                                void *data,
                                void *r,
                                int *tid,
                                int *beg,
                                int *end)

  ctypedef struct hts_idx_t:
      pass

  ctypedef union FilePointerUnion:
      BGZF * bgzf
      cram_fd * cram
      hFILE * hfile
      void * voidp
      
  ctypedef struct htsFile:
      uint32_t is_bin
      int64_t lineno
      kstring_t line
      char * fn
      char * fn_aux
      FilePointerUnion fp

# /*!
#   @abstract  Get the htslib version number
#   @return    For released versions, a string like "N.N[.N]"; or git describe
#   output if using a library built within a Git repository.
# */
  char *hts_version()

# /*!
#   @abstract       Open a SAM/BAM/CRAM/VCF/BCF/etc file
#   @param fn       The file name or "-" for stdin/stdout
#   @param mode     Mode matching /[rwa][bcuz0-9]+/
#   @discussion
#       With 'r' opens for reading; any further format mode letters are ignored
#       as the format is detected by checking the first few bytes or BGZF blocks
#       of the file.  With 'w' or 'a' opens for writing or appending, with format
#       specifier letters:
#         b  binary format (BAM, BCF, etc) rather than text (SAM, VCF, etc)
#         c  CRAM format
#         u  uncompressed
#         z  compressed
#         [0-9]  zlib compression level
#       Note that there is a distinction between 'u' and '0': the first yields
#       plain uncompressed output whereas the latter outputs uncompressed data
#       wrapped in the zlib format.
#   @example
#       [rw]b .. compressed BCF, BAM, FAI
#       [rw]u .. uncompressed BCF
#       [rw]z .. compressed VCF
#       [rw]  .. uncompressed VCF
# */
  htsFile *hts_open(char *fn, char *mode)

# /*!
#   @abstract  Close a file handle, flushing buffered data for output streams
#   @param fp  The file handle to be closed
#   @return    0 for success, or negative if an error occurred.
# */
  int hts_close(htsFile *fp)

# int hts_getline(htsFile *fp, int delimiter, kstring_t *str);
# char **hts_readlines(const char *fn, int *_n);
# /*!
#     @abstract       Parse comma-separated list or read list from a file
#     @param list     File name or comma-separated list
#     @param is_file
#     @param _n       Size of the output array (number of items read)
#     @return         NULL on failure or pointer to newly allocated array of
#                     strings
# */
  char **hts_readlist(char *fn, int is_file, int *_n)

# /*!
#   @abstract  Create extra threads to aid compress/decompression for this file
#   @param fp  The file handle
#   @param n   The number of worker threads to create
#   @return    0 for success, or negative if an error occurred.
#   @notes     THIS THREADING API IS LIKELY TO CHANGE IN FUTURE.
# */
  int hts_set_threads(htsFile *fp, int n)

# /*!
#   @abstract  Set .fai filename for a file opened for reading
#   @return    0 for success, negative on failure
#   @discussion
#       Called before *_hdr_read(), this provides the name of a .fai file
#       used to provide a reference list if the htsFile contains no @SQ headers.
# */
  int hts_set_fai_filename(htsFile *fp, char *fn_aux)

  ctypedef struct hts_pair64_t:
      uint64_t u
      uint64_t v

  ctypedef struct _binstruct:
      int n
      int m
      int *a

  ctypedef struct hts_itr_t:
      #uint32_t read_rest:1, finished:1, dummy:29
      uint32_t flags
      int tid
      int beg
      int end
      int n_off
      int i
      uint64_t curr_off
      hts_pair64_t *off
      hts_readrec_func *readrec
      _binstruct bins

  hts_idx_t *hts_idx_init(int n,
                          int fmt,
                          uint64_t offset0,
                          int min_shift,
                          int n_lvls)

  void hts_idx_destroy(hts_idx_t *idx)

  int hts_idx_push(hts_idx_t *idx,
                   int tid,
                   int beg,
                   int end,
                   uint64_t offset,
                   int is_mapped)

  void hts_idx_finish(hts_idx_t *idx,
                      uint64_t final_offset)
  
  void hts_idx_save(hts_idx_t *idx,
                    char *fn,
                    int fmt)

  hts_idx_t *hts_idx_load(char *fn,
                          int fmt)
        
  uint8_t *hts_idx_get_meta(hts_idx_t *idx,
                            int *l_meta)

  void hts_idx_set_meta(hts_idx_t *idx,
                        int l_meta,
                        uint8_t *meta,
                        int is_copy)

  char *hts_parse_reg(char *s,
                      int *beg,
                      int *end)

  hts_itr_t *hts_itr_query(hts_idx_t *idx,
                           int tid,
                           int beg,
                           int end,
                           hts_readrec_func *readrec)

  void hts_itr_destroy(hts_itr_t *iter)

  ctypedef int (*hts_name2id_f)(void*, char*)
  ctypedef const char *(*hts_id2name_f)(void*, int)
  ctypedef hts_itr_t *hts_itr_query_func(hts_idx_t *idx,
                                         int tid,
                                         int beg,
                                         int end,
                                         hts_readrec_func *readrec)
  
  hts_itr_t *hts_itr_querys(hts_idx_t *idx,
                            char *reg,
                            hts_name2id_f getid,
                            void *hdr,
                            hts_itr_query_func *itr_query,
                            hts_readrec_func *readrec)

  int hts_itr_next(BGZF *fp,
                   hts_itr_t *iter,
                   void *r,
                   void *data)

  char **hts_idx_seqnames(hts_idx_t *idx,
                          int *n,
                          hts_id2name_f getid,
                          void *hdr) 
  
  int hts_file_type(const char *fname)

  int hts_reg2bin(int64_t beg,
                  int64_t end,
                  int min_shift,
                  int n_lvls)

cdef extern from "sam.h":

  # constants
  int BAM_DEF_MASK
  # IF _IOLIB=2, bamFile = BGZF, see bgzf.h
  # samtools uses KNETFILE, check how this works

  ctypedef struct bamFile:
      pass

 # @abstract Structure for core alignment information.
 # @field  tid     chromosome ID, defined by bam_hdr_t
 # @field  pos     0-based leftmost coordinate
 # @field  bin     bin calculated by bam_reg2bin()
 # @field  qual    mapping quality
 # @field  l_qname length of the query name
 # @field  flag    bitwise flag
 # @field  n_cigar number of CIGAR operations
 # @field  l_qseq  length of the query sequence (read)
 # @field  mtid    chromosome ID of next read in template, defined by bam_hdr_t
 # @field  mpos    0-based leftmost coordinate of next read in template
 # Note: 	uint32_t bin:16, qual:8, l_qname:8
 #	        uint32_t flag:16, n_cigar:16
  ctypedef struct bam1_core_t:
      int32_t tid 
      int32_t pos
      uint32_t bin
      uint32_t flag
      int32_t l_qseq
      int32_t mtid 
      int32_t mpos 
      int32_t isize

 # @abstract Structure for one alignment.
 # @field  core       core information about the alignment
 # @field  l_data     current length of bam1_t::data
 # @field  m_data     maximum length of bam1_t::data
 # @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
 
 # @discussion Notes:
 
 # 1. qname is zero tailing and core.l_qname includes the tailing '\0'.
 # 2. l_qseq is calculated from the total length of an alignment block
 # on reading or from CIGAR.
 # 3. cigar data is encoded 4 bytes per CIGAR operation.
 # 4. seq is nybble-encoded according to bam_nt16_table.
  ctypedef struct bam1_t:
    bam1_core_t core
    int l_data
    int m_data
    uint8_t *data
    uint64_t id

  # @abstract Structure for the alignment header.
  # @field n_targets   number of reference sequences
  # @field l_text      length of the plain text in the header
  # @field target_len  lengths of the reference sequences
  # @field target_name names of the reference sequences
  # @field text        plain text
  # @field sdict       header dictionary
  ctypedef struct bam_hdr_t:
      int32_t n_targets
      int32_t ignore_sam_err
      uint32_t l_text
      uint32_t *target_len
      uint8_t * cigar_tab
      char **target_name
      char *text
      void *sdict

  ctypedef struct bam_plbuf_t:
      pass

  ctypedef struct pair64_t:
      uint64_t u, v
      
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

  # BGZF => bamFile

  # HTSLIB definitions
  bam_hdr_t *bam_hdr_init()
  bam_hdr_t *bam_hdr_read(BGZF *fp)
  int bam_hdr_write(BGZF *fp, bam_hdr_t *h)
  void bam_hdr_destroy(bam_hdr_t *h)
  int bam_name2id(bam_hdr_t *h, const char *ref)
  bam_hdr_t* bam_hdr_dup(bam_hdr_t *h0)

  bam1_t *bam_init1()
  void bam_destroy1(bam1_t *b)
  int bam_read1(BGZF *fp, bam1_t *b)
  int bam_write1(BGZF *fp, bam1_t *b)
  bam1_t *bam_copy1(bam1_t *bdst, bam1_t *bsrc)
  bam1_t *bam_dup1(const bam1_t *bsrc)
  
  int bam_cigar2qlen(int n_cigar, uint32_t *cigar)
  int bam_cigar2rlen(int n_cigar, uint32_t *cigar)


  # Iterator interface


  # unmapped functions
  # int64_t bam_seek( bamFile fp, uint64_t voffset, int where)
  # int64_t bam_tell( bamFile fp )

  ###############################################
  # bam iterator interface
  hts_itr_t * sam_itr_queryi(hts_idx_t *idx,
                             int tid,
                             int beg,
                             int end)

  #int bam_iter_read(bamFile fp,
  #                  hts_iter_t iter,
  #                  bam1_t *b)

  ###############################################
  # stand-ins for samtools macros
  uint32_t * bam1_cigar( bam1_t * b)
  char * bam1_qname( bam1_t * b)
  uint8_t * bam1_seq( bam1_t * b)
  uint8_t * bam1_qual( bam1_t * b)
  uint8_t * bam1_aux( bam1_t * b)

  # ###############################################
  # bam_plbuf_t *bam_plbuf_init(bam_pileup_f func,
  #                             void *data)

  # int bam_fetch(bamFile fp,
  #               bam_index_t *idx,
  #               int tid,
  #               int beg,
  #               int end,
  #               void *data,
  #               bam_fetch_f func)

  # int bam_plbuf_push(bam1_t *b,
  #                    bam_plbuf_t *buf)

  # void bam_plbuf_destroy(bam_plbuf_t *buf)
  ########################################
  # pileup iterator interface
  # ctypedef struct bam_plp_t:
  #     pass

  # ctypedef bam_pileup1_t * const_bam_pileup1_t_ptr "const bam_pileup1_t *"

  # ctypedef int (*bam_plp_auto_f)(void *data, bam1_t *b)

  # bam_plp_t bam_plp_init( bam_plp_auto_f func, void *data)
  # int bam_plp_push( bam_plp_t iter,  bam1_t *b)
  # bam_pileup1_t * bam_plp_next( bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
  # bam_pileup1_t * bam_plp_auto( bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
  # void bam_plp_set_mask(bam_plp_t iter, int mask)
  # void bam_plp_reset(bam_plp_t iter)
  # void bam_plp_destroy(bam_plp_t iter)
  # void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt)

  ##################################################
  ## SAM file definitions

  bam_hdr_t *sam_hdr_parse(int l_text,
                           char *text)

  bam_hdr_t *sam_hdr_read(htsFile *fp)

  int sam_hdr_write(htsFile *fp,
                    bam_hdr_t *h)

  # int sam_parse1(kstring_t *s, bam_hdr_t *h, bam1_t *b)

  # int sam_format1(bam_hdr_t *h,
  #                 bam1_t *b,
  #                 kstring_t *str)
  # Read from BAM, CRAM or SAM files
  int sam_read1(htsFile *fp,
                bam_hdr_t *h,
                bam1_t *b)

  # Read from BAM, CRAM or SAM files
  int sam_write1(htsFile *fp,
                 bam_hdr_t *h,
                 bam1_t *b)

  # functions for dealing with the auxillary string
  uint8_t *bam_aux_get(bam1_t *b, char tag[2])
  int32_t bam_aux2i(uint8_t *s)
  float bam_aux2f(uint8_t *s)
  char bam_aux2A( uint8_t *s)
  char *bam_aux2Z( uint8_t *s)

  void bam_aux_append(bam1_t *b,
                      char tag[2],
                      char type,
                      int len,
                      uint8_t *data)

  int bam_aux_del(bam1_t *b, uint8_t *s)

  int32_t bam_endpos(bam1_t *b)

cdef extern from *:
    ctypedef char* const_char_ptr "const char*"




#   ctypedef struct samfile_t_un:
#     tamFile tamr
#     bamFile bam
#     FILE *tamw
    
#   ctypedef struct samfile_t:
#      int type
#      samfile_t_un x
#      bam_hdr_t *header

#   samfile_t *samopen( const_char_ptr fn, char * mode, void *aux)

#   int sampileup( samfile_t *fp, int mask, bam_pileup_f func, void *data)

#   void samclose(samfile_t *fp)

#   int samread(samfile_t *fp, bam1_t *b)

#   int samwrite(samfile_t *fp, bam1_t *b)

#   # functions not declared in sam.h but available as extern
#   int bam_prob_realn(bam1_t *b, char *ref)
#   int bam_cap_mapQ(bam1_t *b, char *ref, int thres)

cdef extern from "faidx.h":

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


cdef extern from "htslib_util.h":

    int pysam_pileup_next(bam1_t *b, 
                          bam_plbuf_t *buf, 
                          bam_pileup1_t ** plp,
                          int * tid,
                          int * pos,
                          int * n_plp )


    # stand-in functions for samtools macros
    void pysam_bam_destroy1( bam1_t * b) 

    # add *nbytes* into the variable length data of *src* at *pos*
    bam1_t * pysam_bam_update( bam1_t * b, 
                               size_t nbytes_old,
                               size_t nbytes_new,
                               uint8_t * pos )

    # FILE * pysam_set_stderr(int fd)
    # void pysam_unset_stderr()

    # return mapped/unmapped reads on tid
    uint32_t pysam_get_mapped(hts_idx_t *idx, int tid)
    uint32_t pysam_get_unmapped(hts_idx_t *idx, int tid)

    # now: static
    int aux_type2size(int)

    char * pysam_bam_get_qname(bam1_t * b)
    uint32_t * pysam_bam_get_cigar(bam1_t * b)
    uint8_t * pysam_bam_get_seq(bam1_t * b)
    uint8_t * pysam_bam_get_qual(bam1_t * b)
    uint8_t * pysam_bam_get_aux(bam1_t * b)
    int pysam_bam_get_l_aux(bam1_t * b)
    char pysam_bam_seqi(uint8_t * s, int i)

    uint16_t pysam_get_bin(bam1_t * b)
    uint8_t pysam_get_qual(bam1_t * b)
    uint8_t pysam_get_l_qname(bam1_t * b)
    uint16_t pysam_get_flag(bam1_t * b)
    uint16_t pysam_get_n_cigar(bam1_t * b)
    void pysam_set_bin(bam1_t * b, uint16_t v)
    void pysam_set_qual(bam1_t * b, uint8_t v)
    void pysam_set_l_qname(bam1_t * b, uint8_t v)
    void pysam_set_flag(bam1_t * b, uint8_t v)
    void pysam_set_n_cigar(bam1_t * b, uint16_t v)
    void pysam_update_flag(bam1_t * b, uint16_t v, uint16_t flag)

####################################################################
# Utility types

# ctypedef struct __iterdata:
#     samfile_t * samfile
#     hts_itr_t iter
#     faidx_t * fastafile
#     int tid
#     char * seq
#     int seq_len



####################################################################
#
# Exposing pysam extension classes
#
# Note: need to declare all C fields and methods here
#
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

cdef class AlignedRead:

    # object that this AlignedRead represents
    cdef bam1_t * _delegate

    # add an alignment tag with value to the AlignedRead 
    # an existing tag of the same name will be replaced.
    cpdef setTag( self, tag, value, value_type = ?, replace = ? )

cdef class Samfile:

    cdef object _filename
    # pointer to htsFile structure
    cdef htsFile * samfile

    # pointer to compressed file
    cdef BGZF * fp

    # pointer to index
    cdef hts_idx_t *index
    # header structure
    cdef bam_hdr_t * header
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

    cdef bam_hdr_t * _buildHeader(self, new_header)
    cdef bam1_t * getCurrent(self)
    cdef int cnext(self)

    # write an aligned read
    cpdef int write(self, AlignedRead read)

    cdef char * _getrname(self, int tid)

cdef class PileupProxy:
    cdef bam_pileup1_t ** plp
    cdef int tid
    cdef int pos
    cdef int n_pu

cdef class PileupRead:
    cdef AlignedRead _alignment
    cdef int32_t  _qpos
    cdef int _indel
    cdef int _level
    cdef uint32_t _is_del
    cdef uint32_t _is_head
    cdef uint32_t _is_tail

cdef class IteratorRow:
    pass

cdef class IteratorRowRegion(IteratorRow):
    cdef hts_itr_t * iter
    cdef bam1_t * b
    cdef int retval
    cdef Samfile samfile
    cdef htsFile * fp
    # true if samfile belongs to this object
    cdef int owns_samfile

    cdef bam1_t * getCurrent( self )

    cdef int cnext(self)

cdef class IteratorRowHead(IteratorRow):
    cdef bam1_t *               b
    cdef int                    retval
    cdef Samfile samfile
    cdef htsFile              * fp
    # true if samfile belongs to this object
    cdef int owns_samfile
    cdef int max_rows
    cdef int current_row

    cdef bam1_t * getCurrent(self)
    cdef int cnext(self)

cdef class IteratorRowAll(IteratorRow):
    cdef bam1_t * b
    cdef htsFile * fp
    cdef Samfile samfile
    cdef int owns_samfile
    cdef bam1_t * getCurrent( self )
    cdef int cnext(self)

cdef class IteratorRowAllRefs(IteratorRow):
    cdef Samfile samfile
    cdef int         tid
    cdef IteratorRowRegion rowiter

cdef class IteratorRowSelection(IteratorRow):
    cdef Samfile samfile
    cdef bam1_t * b
    cdef int current_pos
    cdef htsFile * fp
    cdef positions
    # true if samfile belongs to this object
    cdef int owns_samfile

    cdef bam1_t * getCurrent( self )

    cdef int cnext(self)

# cdef class IteratorColumn:

#     # result of the last plbuf_push
#     cdef IteratorRowRegion iter
#     cdef int tid
#     cdef int pos
#     cdef int n_plp
#     cdef int mask
#     cdef const_bam_pileup1_t_ptr plp
#     cdef bam_plp_t pileup_iter
#     cdef __iterdata iterdata
#     cdef Samfile samfile
#     cdef Fastafile fastafile
#     cdef stepper
#     cdef int max_depth

#     cdef int cnext(self)
#     cdef char * getSequence( self )
#     cdef setMask( self, mask )
#     cdef setupIteratorData( self,
#                             int tid,
#                             int start,
#                             int end,
#                             int reopen = ? )

#     cdef reset( self, tid, start, end )

# cdef class IteratorColumnRegion(IteratorColumn):
#     cdef int start
#     cdef int end
#     cdef int truncate

# cdef class IteratorColumnAllRefs(IteratorColumn):
#     pass

cdef class IndexedReads:
    cdef Samfile samfile
    cdef htsFile * fp
    cdef index
    # true if samfile belongs to this object
    cdef int owns_samfile


