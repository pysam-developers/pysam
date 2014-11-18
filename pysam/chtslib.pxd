from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup
from libc.stdio cimport FILE, printf

cdef extern from "Python.h":
   long _Py_HashPointer(void*)
   FILE* PyFile_AsFile(object)


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

cdef extern from "htslib/kstring.h" nogil:
    ctypedef struct kstring_t:
        size_t l, m
        char *s

cdef extern from "htslib/hfile.h" nogil:
    ctypedef struct hFILE

cdef extern from "htslib/bgzf.h" nogil:
    ctypedef struct bgzf_mtaux_t
    ctypedef struct bgzidx_t
    ctypedef struct z_stream

    ctypedef struct BGZF:
        int           errcode
        int           is_write
        int           is_be
        int           compress_level
        int           is_compressed
        int           is_gzip
        int           cache_size
        int64_t       block_address
        int64_t       uncompressed_address
        void         *uncompressed_block
        void         *compressed_block
        void         *cache
        hFILE        *fp
        bgzf_mtaux_t *mt
        bgzidx_t     *idx
        int           idx_build_otf
        z_stream     *gz_stream

    #*****************
    #  Basic routines *
    # *****************/

    #  Open an existing file descriptor for reading or writing.
    #
    #  @param fd    file descriptor
    #  @param mode  mode matching /[rwa][u0-9]+/: 'r' for reading, 'w' for
    #               writing, or 'a' for appending, while a digit specifies
    #               the zlib compression level.
    #               Note that there is a distinction between 'u' and '0': the
    #               first yields plain uncompressed output whereas the latter
    #               outputs uncompressed data wrapped in the zlib format.
    #  @return      BGZF file handler; 0 on error

    BGZF* bgzf_dopen(int fd, const char *mode)
    BGZF* bgzf_fdopen(int fd, const char *mode) # for backward compatibility

    #  Open the specified file for reading or writing.
    BGZF* bgzf_open(const char* path, const char *mode)

    #  Open an existing hFILE stream for reading or writing.
    BGZF* bgzf_hopen(hFILE *fp, const char *mode)

    #  Close the BGZF and free all associated resources.
    #
    #  @param fp    BGZF file handler
    #  @return      0 on success and -1 on error
    int bgzf_close(BGZF *fp)

    #  Read up to _length_ bytes from the file storing into _data_.
    #
    #  @param fp     BGZF file handler
    #  @param data   data array to read into
    #  @param length size of data to read
    #  @return       number of bytes actually read; 0 on end-of-file and -1 on error
    ssize_t bgzf_read(BGZF *fp, void *data, size_t length)

    #  Write _length_ bytes from _data_ to the file.  If no I/O errors occur,
    #  the complete _length_ bytes will be written (or queued for writing).
    #
    #  @param fp     BGZF file handler
    #  @param data   data array to write
    #  @param length size of data to write
    #  @return       number of bytes written (i.e., _length_); negative on error
    ssize_t bgzf_write(BGZF *fp, const void *data, size_t length)

    #  Read up to _length_ bytes directly from the underlying stream without
    #  decompressing.  Bypasses BGZF blocking, so must be used with care in
    #  specialised circumstances only.
    #
    #  @param fp     BGZF file handler
    #  @param data   data array to read into
    #  @param length number of raw bytes to read
    #  @return       number of bytes actually read; 0 on end-of-file and -1 on error
    ssize_t bgzf_raw_read(BGZF *fp, void *data, size_t length)

    #  Write _length_ bytes directly to the underlying stream without
    #  compressing.  Bypasses BGZF blocking, so must be used with care
    #  in specialised circumstances only.
    #
    #  @param fp     BGZF file handler
    #  @param data   data array to write
    #  @param length number of raw bytes to write
    #  @return       number of bytes actually written; -1 on error
    ssize_t bgzf_raw_write(BGZF *fp, const void *data, size_t length)

    #  Write the data in the buffer to the file.
    int bgzf_flush(BGZF *fp)

    #  Return a virtual file pointer to the current location in the file.
    #  No interpetation of the value should be made, other than a subsequent
    #  call to bgzf_seek can be used to position the file at the same point.
    #  Return value is non-negative on success.
    #define bgzf_tell(fp) (((fp)->block_address << 16) | ((fp)->block_offset & 0xFFFF))
    int64_t bgzf_tell(BGZF * fp)

    #  Set the file to read from the location specified by _pos_.
    #
    #  @param fp     BGZF file handler
    #  @param pos    virtual file offset returned by bgzf_tell()
    #  @param whence must be SEEK_SET
    #  @return       0 on success and -1 on error
    # /
    int64_t bgzf_seek(BGZF *fp, int64_t pos, int whence)

    #  Check if the BGZF end-of-file (EOF) marker is present
    #
    #  @param fp    BGZF file handler opened for reading
    #  @return      1 if the EOF marker is present and correct
    #               2 if it can't be checked, e.g., because fp isn't seekable
    #               0 if the EOF marker is absent
    #               -1 (with errno set) on error
    int bgzf_check_EOF(BGZF *fp)

    #  Check if a file is in the BGZF format
    #
    #  @param fn    file name
    #  @return      1 if _fn_ is BGZF; 0 if not or on I/O error
    int bgzf_is_bgzf(const char *fn)

    #*********************
    #  Advanced routines *
    #*********************

    #  Set the cache size. Only effective when compiled with -DBGZF_CACHE.
    #
    #  @param fp    BGZF file handler
    #  @param size  size of cache in bytes; 0 to disable caching (default)
    void bgzf_set_cache_size(BGZF *fp, int size)

    #  Flush the file if the remaining buffer size is smaller than _size_
    #  @return      0 if flushing succeeded or was not needed; negative on error
    int bgzf_flush_try(BGZF *fp, ssize_t size)

    #  Read one byte from a BGZF file. It is faster than bgzf_read()
    #  @param fp     BGZF file handler
    #  @return       byte read; -1 on end-of-file or error
    int bgzf_getc(BGZF *fp)

    #  Read one line from a BGZF file. It is faster than bgzf_getc()
    #
    #  @param fp     BGZF file handler
    #  @param delim  delimitor
    #  @param str    string to write to; must be initialized
    #  @return       length of the string; 0 on end-of-file; negative on error
    int bgzf_getline(BGZF *fp, int delim, kstring_t *str)

    #  Read the next BGZF block.
    int bgzf_read_block(BGZF *fp)

    #  Enable multi-threading (only effective on writing and when the
    #  library was compiled with -DBGZF_MT)
    #
    #  @param fp          BGZF file handler; must be opened for writing
    #  @param n_threads   #threads used for writing
    #  @param n_sub_blks  #blocks processed by each thread; a value 64-256 is recommended
    int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks)


    #*******************
    #  bgzidx routines *
    #   BGZF at the uncompressed offset
    #
    #   @param fp           BGZF file handler; must be opened for reading
    #   @param uoffset      file offset in the uncompressed data
    #   @param where        SEEK_SET supported atm
    #
    #   Returns 0 on success and -1 on error.
    int bgzf_useek(BGZF *fp, long uoffset, int where)

    #   Position in uncompressed BGZF
    #
    #   @param fp           BGZF file handler; must be opened for reading
    #
    #   Returns the current offset on success and -1 on error.
    long bgzf_utell(BGZF *fp)

    #  Tell BGZF to build index while compressing.
    #
    #  @param fp          BGZF file handler; can be opened for reading or writing.
    #
    #  Returns 0 on success and -1 on error.
    int bgzf_index_build_init(BGZF *fp)

    #  Load BGZF index
    #
    #  @param fp          BGZF file handler
    #  @param bname       base name
    #  @param suffix      suffix to add to bname (can be NULL)
    #
    #  Returns 0 on success and -1 on error.
    int bgzf_index_load(BGZF *fp, const char *bname, const char *suffix)

    #  Save BGZF index
    #
    #  @param fp          BGZF file handler
    #  @param bname       base name
    #  @param suffix      suffix to add to bname (can be NULL)
    #
    #  Returns 0 on success and -1 on error.
    int bgzf_index_dump(BGZF *fp, const char *bname, const char *suffix)


cdef extern from "htslib/hts.h" nogil:
    uint32_t kroundup32(uint32_t x)

    ctypedef struct cram_fd

    ctypedef union FilePointerUnion:
        BGZF * bgzf
        cram_fd * cram
        hFILE * hfile
        void * voidp

    ctypedef struct htsFile:
	# uint32_t is_bin:1, is_write:1, is_be:1, is_cram:1, is_compressed:2, is_kstream:1, dummy:25;
        uint32_t  is_bin
        int64_t lineno
        kstring_t line
        char * fn
        char * fn_aux
        FilePointerUnion fp

    int hts_verbose

    # @abstract Table for converting a nucleotide character to the 4-bit encoding.
    const unsigned char *seq_nt16_table

    # @abstract Table for converting a 4-bit encoded nucleotide to a letter.
    const char *seq_nt16_str

    # @abstract  Get the htslib version number
    # @return    For released versions, a string like "N.N[.N]"; or git describe
    # output if using a library built within a Git repository.
    const char *hts_version()

    # @abstract       Open a SAM/BAM/CRAM/VCF/BCF/etc file
    # @param fn       The file name or "-" for stdin/stdout
    # @param mode     Mode matching /[rwa][bcuz0-9]+/
    # @discussion
    #     With 'r' opens for reading; any further format mode letters are ignored
    #     as the format is detected by checking the first few bytes or BGZF blocks
    #     of the file.  With 'w' or 'a' opens for writing or appending, with format
    #     specifier letters:
    #       b  binary format (BAM, BCF, etc) rather than text (SAM, VCF, etc)
    #       c  CRAM format
    #       u  uncompressed
    #       z  compressed
    #       [0-9]  zlib compression level
    #     Note that there is a distinction between 'u' and '0': the first yields
    #     plain uncompressed output whereas the latter outputs uncompressed data
    #     wrapped in the zlib format.
    # @example
    #     [rw]b .. compressed BCF, BAM, FAI
    #     [rw]u .. uncompressed BCF
    #     [rw]z .. compressed VCF
    #     [rw]  .. uncompressed VCF
    htsFile *hts_open(const char *fn, const char *mode)

    # @abstract  Close a file handle, flushing buffered data for output streams
    # @param fp  The file handle to be closed
    # @return    0 for success, or negative if an error occurred.
    int hts_close(htsFile *fp)

    int hts_getline(htsFile *fp, int delimiter, kstring_t *str)
    char **hts_readlines(const char *fn, int *_n)

    #   @abstract       Parse comma-separated list or read list from a file
    #   @param list     File name or comma-separated list
    #   @param is_file
    #   @param _n       Size of the output array (number of items read)
    #   @return         NULL on failure or pointer to newly allocated array of
    #                   strings
    char **hts_readlist(const char *fn, int is_file, int *_n)

    # @abstract  Create extra threads to aid compress/decompression for this file
    # @param fp  The file handle
    # @param n   The number of worker threads to create
    # @return    0 for success, or negative if an error occurred.
    # @notes     THIS THREADING API IS LIKELY TO CHANGE IN FUTURE.
    int hts_set_threads(htsFile *fp, int n)

    # @abstract  Set .fai filename for a file opened for reading
    # @return    0 for success, negative on failure
    # @discussion
    #     Called before *_hdr_read(), this provides the name of a .fai file
    #     used to provide a reference list if the htsFile contains no @SQ headers.
    int hts_set_fai_filename(htsFile *fp, const char *fn_aux)

    int8_t HTS_IDX_NOCOOR
    int8_t HTS_IDX_START
    int8_t HTS_IDX_REST
    int8_t HTS_IDX_NONE

    int8_t HTS_FMT_CSI
    int8_t HTS_FMT_BAI
    int8_t HTS_FMT_TBI
    int8_t HTS_FMT_CRAI

    ctypedef struct hts_idx_t

    ctypedef struct hts_pair64_t:
        uint64_t u, v

    ctypedef int hts_readrec_func(BGZF *fp, void *data, void *r, int *tid, int *beg, int *end)

    ctypedef struct hts_bins_t:
        int n, m
        int *a

    ctypedef struct hts_itr_t:
        uint32_t read_rest
        uint32_t finished
        int tid, bed, end, n_off, i
        uint64_t curr_off
        hts_pair64_t *off
        hts_readrec_func *readfunc
        hts_bins_t bins

    hts_idx_t *hts_idx_init(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls)
    void hts_idx_destroy(hts_idx_t *idx)
    int hts_idx_push(hts_idx_t *idx, int tid, int beg, int end, uint64_t offset, int is_mapped)
    void hts_idx_finish(hts_idx_t *idx, uint64_t final_offset)

    void hts_idx_save(const hts_idx_t *idx, const char *fn, int fmt)
    hts_idx_t *hts_idx_load(const char *fn, int fmt)

    uint8_t *hts_idx_get_meta(hts_idx_t *idx, int *l_meta)
    void hts_idx_set_meta(hts_idx_t *idx, int l_meta, uint8_t *meta, int is_copy)

    int hts_idx_get_stat(const hts_idx_t* idx, int tid,
                         uint64_t* mapped, uint64_t* unmapped)

    uint64_t hts_idx_get_n_no_coor(const hts_idx_t* idx)

    const char *hts_parse_reg(const char *s, int *beg, int *end)
    hts_itr_t *hts_itr_query(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec)
    void hts_itr_destroy(hts_itr_t *iter)

    ctypedef int (*hts_name2id_f)(void*, const char*)
    ctypedef const char *(*hts_id2name_f)(void*, int)
    ctypedef hts_itr_t *hts_itr_query_func(
        const hts_idx_t *idx,
        int tid,
        int beg,
        int end,
        hts_readrec_func *readrec)

    hts_itr_t *hts_itr_querys(
        const hts_idx_t *idx,
        const char *reg,
        hts_name2id_f getid,
        void *hdr,
        hts_itr_query_func *itr_query,
        hts_readrec_func *readrec)

    int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data)
    const char **hts_idx_seqnames(const hts_idx_t *idx, int *n, hts_id2name_f getid, void *hdr)  # free only the array, not the values

    # hts_file_type() - Convenience function to determine file type
    # @fname: the file name
    #
    # Returns one of the FT_* defines.
    #
    # This function was added in order to avoid the need for excessive command
    # line switches.
    int FT_UNKN
    int FT_GZ
    int FT_VCF
    int FT_VCF_GZ
    int FT_BCF
    int FT_BCF_GZ
    int FT_STDIN

    int hts_file_type(const char *fname)

    inline int hts_reg2bin(int64_t beg, int64_t end, int min_shift, int n_lvls)
    inline int hts_bin_bot(int bin, int n_lvls)

    # * Endianness *
    inline int ed_is_big()
    inline uint16_t ed_swap_2(uint16_t v)
    inline void *ed_swap_2p(void *x)
    inline uint32_t ed_swap_4(uint32_t v)
    inline void *ed_swap_4p(void *x)
    inline uint64_t ed_swap_8(uint64_t v)
    inline void *ed_swap_8p(void *x)


cdef extern from "htslib/sam.h" nogil:
    #**********************
    #*** SAM/BAM header ***
    #**********************

    # @abstract Structure for the alignment header.
    # @field n_targets   number of reference sequences
    # @field l_text      length of the plain text in the header
    # @field target_len  lengths of the reference sequences
    # @field target_name names of the reference sequences
    # @field text        plain text
    # @field sdict       header dictionary

    ctypedef struct bam_hdr_t:
         int32_t n_targets, ignore_sam_err
         uint32_t l_text
         uint32_t *target_len
         uint8_t *cigar_tab
         char **target_name
         char *text
         void *sdict

    #****************************
    #*** CIGAR related macros ***
    #****************************

    int BAM_CMATCH
    int BAM_CINS
    int BAM_CDEL
    int BAM_CREF_SKIP
    int BAM_CSOFT_CLIP
    int BAM_CHARD_CLIP
    int BAM_CPAD
    int BAM_CEQUAL
    int BAM_CDIFF
    int BAM_CBACK

    char    *BAM_CIGAR_STR
    int      BAM_CIGAR_SHIFT
    uint32_t BAM_CIGAR_MASK
    uint32_t BAM_CIGAR_TYPE

    char bam_cigar_op(uint32_t c)
    uint32_t bam_cigar_oplen(uint32_t c)
    char bam_cigar_opchr(uint32_t)
    uint32_t bam_cigar_gen(char, uint32_t)
    int bam_cigar_type(char o)

    # @abstract the read is paired in sequencing, no matter whether it is mapped in a pair
    int BAM_FPAIRED
    # @abstract the read is mapped in a proper pair
    int BAM_FPROPER_PAIR
    # @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR
    int BAM_FUNMAP
    # @abstract the mate is unmapped
    int BAM_FMUNMAP
    # @abstract the read is mapped to the reverse strand
    int BAM_FREVERSE
    # @abstract the mate is mapped to the reverse strand
    int BAM_FMREVERSE
    # @abstract this is read1
    int BAM_FREAD1
    # @abstract this is read2
    int BAM_FREAD2
    # @abstract not primary alignment
    int BAM_FSECONDARY
    # @abstract QC failure
    int BAM_FQCFAIL
    # @abstract optical or PCR duplicate
    int BAM_FDUP
    # @abstract supplementary alignment
    int BAM_FSUPPLEMENTARY

    #*************************
    #*** Alignment records ***
    #*************************

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

    ctypedef struct bam1_core_t:
        int32_t tid
        int32_t pos
        uint16_t bin
        uint8_t qual
        uint8_t l_qname
        uint16_t flag
        uint16_t n_cigar
        int32_t l_qseq
        int32_t mtid
        int32_t mpos
        int32_t isize

    # @abstract Structure for one alignment.
    # @field  core       core information about the alignment
    # @field  l_data     current length of bam1_t::data
    # @field  m_data     maximum length of bam1_t::data
    # @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
    #
    # @discussion Notes:
    #
    # 1. qname is zero tailing and core.l_qname includes the tailing '\0'.
    # 2. l_qseq is calculated from the total length of an alignment block
    # on reading or from CIGAR.
    # 3. cigar data is encoded 4 bytes per CIGAR operation.
    # 4. seq is nybble-encoded according to seq_nt16_table.
    ctypedef struct bam1_t:
        bam1_core_t core
        int l_data, m_data
        uint8_t *data
        uint64_t id

    # @abstract  Get whether the query is on the reverse strand
    # @param  b  pointer to an alignment
    # @return    boolean true if query is on the reverse strand
    int bam_is_rev(bam1_t *b)

    # @abstract  Get whether the query's mate is on the reverse strand
    # @param  b  pointer to an alignment
    # @return    boolean true if query's mate on the reverse strand
    int bam_is_mrev(bam1_t *b)

    # @abstract  Get the name of the query
    # @param  b  pointer to an alignment
    # @return    pointer to the name string, null terminated
    char *bam_get_qname(bam1_t *b)

    # @abstract  Get the CIGAR array
    # @param  b  pointer to an alignment
    # @return    pointer to the CIGAR array
    #
    # @discussion In the CIGAR array, each element is a 32-bit integer. The
    # lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
    # length of a CIGAR.
    uint32_t *bam_get_cigar(bam1_t *b)

    # @abstract  Get query sequence
    # @param  b  pointer to an alignment
    # @return    pointer to sequence
    #
    # @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
    # 8 for T and 15 for N. Two bases are packed in one byte with the base
    # at the higher 4 bits having smaller coordinate on the read. It is
    # recommended to use bam_seqi() macro to get the base.
    char *bam_get_seq(bam1_t *b)

    # @abstract  Get query quality
    # @param  b  pointer to an alignment
    # @return    pointer to quality string
    uint8_t *bam_get_qual(bam1_t *b)

    # @abstract  Get auxiliary data
    # @param  b  pointer to an alignment
    # @return    pointer to the concatenated auxiliary data
    uint8_t *bam_get_aux(bam1_t *b)

    # @abstract  Get length of auxiliary data
    # @param  b  pointer to an alignment
    # @return    length of the concatenated auxiliary data
    int bam_get_l_aux(bam1_t *b)

    # @abstract  Get a base on read
    # @param  s  Query sequence returned by bam1_seq()
    # @param  i  The i-th position, 0-based
    # @return    4-bit integer representing the base.
    char bam_seqi(char *s, int i)

    #**************************
    #*** Exported functions ***
    #**************************

    #***************
    #*** BAM I/O ***
    #***************

    bam_hdr_t *bam_hdr_init()
    bam_hdr_t *bam_hdr_read(BGZF *fp)
    int bam_hdr_write(BGZF *fp, const bam_hdr_t *h)
    void bam_hdr_destroy(bam_hdr_t *h)
    int bam_name2id(bam_hdr_t *h, const char *ref)
    bam_hdr_t* bam_hdr_dup(const bam_hdr_t *h0)

    bam1_t *bam_init1()
    void bam_destroy1(bam1_t *b)
    int bam_read1(BGZF *fp, bam1_t *b)
    int bam_write1(BGZF *fp, const bam1_t *b)
    bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc)
    bam1_t *bam_dup1(const bam1_t *bsrc)

    int bam_cigar2qlen(int n_cigar, const uint32_t *cigar)
    int bam_cigar2rlen(int n_cigar, const uint32_t *cigar)

    # @abstract Calculate the rightmost base position of an alignment on the
    # reference genome.

    # @param  b  pointer to an alignment
    # @return    the coordinate of the first base after the alignment, 0-based

    # @discussion For a mapped read, this is just b->core.pos + bam_cigar2rlen.
    # For an unmapped read (either according to its flags or if it has no cigar
    # string), we return b->core.pos + 1 by convention.
    int32_t bam_endpos(const bam1_t *b)

    int   bam_str2flag(const char *str)  # returns negative value on error
    char *bam_flag2str(int flag)         # The string must be freed by the user

    #*************************
    #*** BAM/CRAM indexing ***
    #*************************

    # These BAM iterator functions work only on BAM files.  To work with either
    # BAM or CRAM files use the sam_index_load() & sam_itr_*() functions.
    void bam_itr_destroy(hts_itr_t *iter)
    hts_itr_t *bam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end)
    hts_itr_t *bam_itr_querys(const hts_idx_t *idx, bam_hdr_t *hdr, const char *region)
    int bam_itr_next(htsFile *htsfp, hts_itr_t *itr, void *r)

    # Load .csi or .bai BAM index file.
    hts_idx_t *bam_idx_load(const char *fn)
    int bam_index_build(const char *fn, int min_shift)

    # Load BAM (.csi or .bai) or CRAM (.crai) index file.
    hts_idx_t *sam_index_load(htsFile *fp, const char *fn)

    void sam_itr_destroy(hts_itr_t *iter)
    hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end)
    hts_itr_t *sam_itr_querys(const hts_idx_t *idx, bam_hdr_t *hdr, const char *region)
    int sam_itr_next(htsFile *htsfp, hts_itr_t *itr, void *r)

    #***************
    #*** SAM I/O ***
    #***************

    htsFile *sam_open(const char *fn, const char *mode)
    int sam_close(htsFile *fp)

    bam_hdr_t *sam_hdr_parse(int l_text, const char *text)
    bam_hdr_t *sam_hdr_read(htsFile *fp)
    int sam_hdr_write(htsFile *fp, const bam_hdr_t *h)

    int sam_parse1(kstring_t *s, bam_hdr_t *h, bam1_t *b)
    int sam_format1(const bam_hdr_t *h, const bam1_t *b, kstring_t *str)
    int sam_read1(htsFile *fp, bam_hdr_t *h, bam1_t *b)
    int sam_write1(htsFile *fp, const bam_hdr_t *h, const bam1_t *b)

    #*************************************
    #*** Manipulating auxiliary fields ***
    #*************************************

    uint8_t *bam_aux_get(const bam1_t *b, const char *tag)
    int32_t  bam_aux2i(const uint8_t *s)
    double   bam_aux2f(const uint8_t *s)
    char     bam_aux2A(const uint8_t *s)
    char    *bam_aux2Z(const uint8_t *s)

    void bam_aux_append(bam1_t *b, const char *tag, char type, int len, uint8_t *data)
    int bam_aux_del(bam1_t *b, uint8_t *s)

    #**************************
    #*** Pileup and Mpileup ***
    #**************************

    # @abstract Structure for one alignment covering the pileup position.
    # @field  b          pointer to the alignment
    # @field  qpos       position of the read base at the pileup site, 0-based
    # @field  indel      indel length; 0 for no indel, positive for ins and negative for del
    # @field  level      the level of the read in the "viewer" mode
    # @field  is_del     1 iff the base on the padded read is a deletion
    # @field  is_head    ???
    # @field  is_tail    ???
    # @field  is_refskip ???
    # @field  aux        ???
    #
    # @discussion See also bam_plbuf_push() and bam_lplbuf_push(). The
    # difference between the two functions is that the former does not
    # set bam_pileup1_t::level, while the later does. Level helps the
    # implementation of alignment viewers, but calculating this has some
    # overhead.
    # 
    # is_del, is_head, etc are a bit field, declaring as below should
    # work as expected, see
    # https://groups.google.com/forum/#!msg/cython-users/24tD1kwRY7A/pmoPuSmanM0J

    ctypedef struct bam_pileup1_t:
        bam1_t *b
        int32_t qpos
        int indel, level
        uint32_t is_del
        uint32_t is_head
        uint32_t is_tail
        uint32_t is_refskip
        uint32_t aux

    ctypedef int (*bam_plp_auto_f)(void *data, bam1_t *b)
    ctypedef int (*bam_test_f)()

    ctypedef struct __bam_plp_t
    ctypedef __bam_plp_t *bam_plp_t

    ctypedef struct __bam_mplp_t
    ctypedef __bam_mplp_t *bam_mplp_t

    # bam_plp_init() - sets an iterator over multiple
    # @func:      see mplp_func in bam_plcmd.c in samtools for an example. Expected return
    #             status: 0 on success, -1 on end, < -1 on non-recoverable errors
    # @data:      user data to pass to @func
    bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data)
    void bam_plp_destroy(bam_plp_t iter)
    int bam_plp_push(bam_plp_t iter, const bam1_t *b)
    const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
    const bam_pileup1_t *bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
    void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt)
    void bam_plp_reset(bam_plp_t iter)

    bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data)

    # bam_mplp_init_overlaps() - if called, mpileup will detect overlapping
    # read pairs and for each base pair set the base quality of the
    # lower-quality base to zero, thus effectively discarding it from
    # calling. If the two bases are identical, the quality of the other base
    # is increased to the sum of their qualities (capped at 200), otherwise
    # it is multiplied by 0.8.
    void bam_mplp_init_overlaps(bam_mplp_t iter)
    void bam_mplp_destroy(bam_mplp_t iter)
    void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt)
    int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp)

    # Added by AH
    # ctypedef bam_pileup1_t * const_bam_pileup1_t_ptr "const bam_pileup1_t *"

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
    void ks_destroy(kstream_t *)

    # Retrieve characters from stream until delimiter
    # is reached placing results in str.
    int ks_getuntil(kstream_t *,
                    int delimiter,
                    kstring_t * str,
                    int * dret)

cdef extern from "htslib/faidx.h":

    ctypedef struct faidx_t:
       pass

    int fai_build(char *fn)

    void fai_destroy(faidx_t *fai)

    faidx_t *fai_load(char *fn)

    char *fai_fetch(faidx_t *fai,
                    char *reg,
                    int *len)

    int faidx_nseq(faidx_t *fai)

    int faidx_has_seq(faidx_t *fai, const char *seq)

    char *faidx_fetch_seq(faidx_t *fai,
                         char *c_name,
                         int p_beg_i,
                         int p_end_i,
                         int *len)

    int faidx_seq_len(faidx_t *fai, const char *seq)

# tabix support
cdef extern from "htslib/tbx.h" nogil:
    
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

