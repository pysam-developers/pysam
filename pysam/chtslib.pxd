from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup
from libc.stdio cimport FILE, printf
from posix.types cimport off_t

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
  int gzeof(gzFile file)

cdef extern from "htslib/kstring.h" nogil:
    ctypedef struct kstring_t:
        size_t l, m
        char *s

cdef extern from "htslib_util.h" nogil:
    ctypedef uint32_t khint32_t
    ctypedef uint32_t khint_t
    ctypedef khint_t  khiter_t

    # Used to manage BCF Header info
    ctypedef struct vdict_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        khint32_t *flags
        const char *keys
        bcf_idinfo_t *vals

    # Used to manage indexed contigs in Tabix
    ctypedef struct s2i_t:
        khint_t n_buckets, size, n_occupied, upper_bound
        khint32_t *flags
        const char *keys
        int64_t *vals

    # Generic khash methods
    khint_t kh_size(void *d)
    khint_t kh_begin(void *d)
    khint_t kh_end(void *d)
    int kh_exist(void *d, khiter_t i)

    # Specialized khash methods for vdict
    khint_t kh_get_vdict(vdict_t *d, const char *key)
    const char *kh_key_vdict "kh_key" (vdict_t *d, khint_t i)
    bcf_idinfo_t kh_val_vdict "kh_val" (vdict_t *d, khint_t i)


cdef extern from "htslib/hfile.h" nogil:
    ctypedef struct hFILE

    # @abstract  Open the named file or URL as a stream
    # @return    An hFILE pointer, or NULL (with errno set) if an error occurred.
    hFILE *hopen(const char *filename, const char *mode)

    # @abstract  Associate a stream with an existing open file descriptor
    # @return    An hFILE pointer, or NULL (with errno set) if an error occurred.
    # @notes     For socket descriptors (on Windows), mode should contain 's'.
    hFILE *hdopen(int fd, const char *mode)

    # @abstract  Report whether the file name or URL denotes remote storage
    # @return    0 if local, 1 if remote.
    # @notes     "Remote" means involving e.g. explicit network access, with the
    #   implication that callers may wish to cache such files' contents locally.
    int hisremote(const char *filename)

    # @abstract  Flush (for output streams) and close the stream
    # @return    0 if successful, or EOF (with errno set) if an error occurred.
    int hclose(hFILE *fp)

    # @abstract  Close the stream, without flushing or propagating errors
    # @notes     For use while cleaning up after an error only.  Preserves errno.
    void hclose_abruptly(hFILE *fp)

    # @abstract  Return the stream's error indicator
    # @return    Non-zero (in fact, an errno value) if an error has occurred.
    # @notes     This would be called herror() and return true/false to parallel
    #   ferror(3), but a networking-related herror(3) function already exists.  */
    int herrno(hFILE *fp)

    # @abstract  Clear the stream's error indicator
    void hclearerr(hFILE *fp)

    # @abstract  Reposition the read/write stream offset
    # @return    The resulting offset within the stream (as per lseek(2)),
    #   or negative if an error occurred.
    off_t hseek(hFILE *fp, off_t offset, int whence)

    # @abstract  Report the current stream offset
    # @return    The offset within the stream, starting from zero.
    off_t htell(hFILE *fp)

    # @abstract  Read one character from the stream
    # @return    The character read, or EOF on end-of-file or error
    int hgetc(hFILE *fp)

    # @abstract  Peek at characters to be read without removing them from buffers
    # @param fp      The file stream
    # @param buffer  The buffer to which the peeked bytes will be written
    # @param nbytes  The number of bytes to peek at; limited by the size of the
    #   internal buffer, which could be as small as 4K.
    # @return    The number of bytes peeked, which may be less than nbytes if EOF
    #   is encountered; or negative, if there was an I/O error.
    # @notes  The characters peeked at remain in the stream's internal buffer,
    #   and will be returned by later hread() etc calls.
    ssize_t hpeek(hFILE *fp, void *buffer, size_t nbytes)

    # @abstract  Read a block of characters from the file
    # @return    The number of bytes read, or negative if an error occurred.
    # @notes     The full nbytes requested will be returned, except as limited
    #   by EOF or I/O errors.
    ssize_t hread(hFILE *fp, void *buffer, size_t nbytes)

    # @abstract  Write a character to the stream
    # @return    The character written, or EOF if an error occurred.
    int hputc(int c, hFILE *fp)

    # @abstract  Write a string to the stream
    # @return    0 if successful, or EOF if an error occurred.
    int hputs(const char *text, hFILE *fp)

    # @abstract  Write a block of characters to the file
    # @return    Either nbytes, or negative if an error occurred.
    # @notes     In the absence of I/O errors, the full nbytes will be written.
    ssize_t hwrite(hFILE *fp, const void *buffer, size_t nbytes)

    # @abstract  For writing streams, flush buffered output to the underlying stream
    # @return    0 if successful, or EOF if an error occurred.
    int hflush(hFILE *fp)


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
    #  @param mode  mode matching /[rwag][u0-9]+/: 'r' for reading, 'w' for
    #               writing, 'a' for appending, 'g' for gzip rather than BGZF
    #               compression (with 'w' only), and digit specifies the zlib
    #               compression level.
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

    int SEEK_SET

    #  Return a virtual file pointer to the current location in the file.
    #  No interpetation of the value should be made, other than a subsequent
    #  call to bgzf_seek can be used to position the file at the same point.
    #  Return value is non-negative on success.
    int64_t bgzf_tell(BGZF *fp)

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
        BGZF    *bgzf
        cram_fd *cram
        hFILE   *hfile
        void    *voidp

    ctypedef enum htsFormatCategory:
        unknown_category
        sequence_data    # Sequence data -- SAM, BAM, CRAM, etc
        variant_data     # Variant calling data -- VCF, BCF, etc
        index_file       # Index file associated with some data file
        region_list      # Coordinate intervals or regions -- BED, etc
        category_maximum

    ctypedef enum htsExactFormat:
        unknown_format
        binary_format
        text_format
        sam, bam, bai, cram, crai, vcf, bcf, csi, gzi, tbi, bed
        format_maximum

    ctypedef enum htsCompression:
        no_compression, gzip, bgzf, custom
        compression_maximum

    cdef struct htsVersion:
        short major, minor

    ctypedef struct htsFormat:
        htsFormatCategory category
        htsExactFormat    format
        htsVersion        version
        htsCompression    compression

    ctypedef struct htsFile:
        uint8_t  is_bin
        uint8_t  is_write
        uint8_t  is_be
        uint8_t  is_cram
        int64_t lineno
        kstring_t line
        char *fn
        char *fn_aux
        FilePointerUnion fp
        htsFormat format

    int hts_verbose

    # @abstract Table for converting a nucleotide character to 4-bit encoding.
    # The input character may be either an IUPAC ambiguity code, '=' for 0, or
    # '0'/'1'/'2'/'3' for a result of 1/2/4/8.  The result is encoded as 1/2/4/8
    # for A/C/G/T or combinations of these bits for ambiguous bases.
    const unsigned char *seq_nt16_table

    # @abstract Table for converting a 4-bit encoded nucleotide to an IUPAC
    # ambiguity code letter (or '=' when given 0).
    const char *seq_nt16_str

    # @abstract Table for converting a 4-bit encoded nucleotide to about 2 bits.
    # Returns 0/1/2/3 for 1/2/4/8 (i.e., A/C/G/T), or 4 otherwise (0 or ambiguous).
    const int *seq_nt16_int

    # @abstract  Get the htslib version number
    # @return    For released versions, a string like "N.N[.N]"; or git describe
    # output if using a library built within a Git repository.
    const char *hts_version()

    # @abstract    Determine format by peeking at the start of a file
    # @param fp    File opened for reading, positioned at the beginning
    # @param fmt   Format structure that will be filled out on return
    # @return      0 for success, or negative if an error occurred.
    int hts_detect_format(hFILE *fp, htsFormat *fmt)

    # @abstract    Get a human-readable description of the file format
    # @return      Description string, to be freed by the caller after use.
    char *hts_format_description(const htsFormat *format)

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
    #       g  gzip compressed
    #       u  uncompressed
    #       z  bgzf compressed
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

    # @abstract       Open an existing stream as a SAM/BAM/CRAM/VCF/BCF/etc file
    # @param fp       The already-open file handle
    # @param fn       The file name or "-" for stdin/stdout
    # @param mode     Open mode, as per hts_open()
    htsFile *hts_hopen(hFILE *fp, const char *fn, const char *mode)

    # @abstract  Close a file handle, flushing buffered data for output streams
    # @param fp  The file handle to be closed
    # @return    0 for success, or negative if an error occurred.
    int hts_close(htsFile *fp)

    # @abstract  Returns the file's format information
    # @param fp  The file handle
    # @return    Read-only pointer to the file's htsFormat.
    const htsFormat *hts_get_format(htsFile *fp)

    # @abstract  Sets a specified CRAM option on the open file handle.
    # @param fp  The file handle open the open file.
    # @param opt The CRAM_OPT_* option.
    # @param ... Optional arguments, dependent on the option used.
    # @return    0 for success, or negative if an error occurred.
    #int hts_set_opt(htsFile *fp, enum cram_option opt, ...)

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

    BGZF *hts_get_bgzfp(htsFile *fp)
    int hts_useek(htsFile *fp, long uoffset, int where)
    long hts_utell(htsFile *fp)

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
        int curr_tid, curr_beg, curr_end
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
    # DEPRECATED:  This function has been replaced by hts_detect_format().
    # It and these FT_* macros will be removed in a future HTSlib release.
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
    kseq_t *kseq_init(gzFile)
    int kseq_read(kseq_t *)
    void kseq_destroy(kseq_t *)
    int gzclose(gzFile)

    kstream_t *ks_init(gzFile)
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


# VCF/BCF API
cdef extern from "htslib/vcf.h" nogil:

    # Header struct

    uint8_t BCF_HL_FLT   # header line
    uint8_t BCF_HL_INFO
    uint8_t BCF_HL_FMT
    uint8_t BCF_HL_CTG
    uint8_t BCF_HL_STR   # structured header line TAG=<A=..,B=..>
    uint8_t BCF_HL_GEN   # generic header line

    uint8_t BCF_HT_FLAG  # header type
    uint8_t BCF_HT_INT
    uint8_t BCF_HT_REAL
    uint8_t BCF_HT_STR

    uint8_t BCF_VL_FIXED # variable length
    uint8_t BCF_VL_VAR
    uint8_t BCF_VL_A
    uint8_t BCF_VL_G
    uint8_t BCF_VL_R

    # === Dictionary ===
    #
    # The header keeps three dictonaries. The first keeps IDs in the
    # "FILTER/INFO/FORMAT" lines, the second keeps the sequence names and lengths
    # in the "contig" lines and the last keeps the sample names. bcf_hdr_t::dict[]
    # is the actual hash table, which is opaque to the end users. In the hash
    # table, the key is the ID or sample name as a C string and the value is a
    # bcf_idinfo_t struct. bcf_hdr_t::id[] points to key-value pairs in the hash
    # table in the order that they appear in the VCF header. bcf_hdr_t::n[] is the
    # size of the hash table or, equivalently, the length of the id[] arrays.

    uint8_t BCF_DT_ID       # dictionary type
    uint8_t BCF_DT_CTG
    uint8_t BCF_DT_SAMPLE

    # Complete textual representation of a header line
    ctypedef struct bcf_hrec_t:
        int type            # One of the BCF_HL_* type
        char *key           # The part before '=', i.e. FILTER/INFO/FORMAT/contig/fileformat etc.
        char *value         # Set only for generic lines, NULL for FILTER/INFO, etc.
        int nkeys           # Number of structured fields
        char **keys         # The key=value pairs
        char **vals

    ctypedef struct bcf_idinfo_t:
        uint32_t info[3]     # stores Number:20, var:4, Type:4, ColType:4 in info[0..2]
        bcf_hrec_t *hrec[3]  # for BCF_HL_FLT,INFO,FMT and contig length in info[0] for BCF_HL_CTG
        int id

    ctypedef struct bcf_idpair_t:
        const char *key
        const bcf_idinfo_t *val

    ctypedef struct bcf_hdr_t:
        int32_t n[3]
        bcf_idpair_t *id[3]
        void *dict[3]               # ID dictionary, contig dict and sample dict
        char **samples
        bcf_hrec_t **hrec
        int nhrec, dirty
        int ntransl
        int *transl[2]              # for bcf_translate()
        int nsamples_ori            # for bcf_hdr_set_samples()
        uint8_t *keep_samples
        kstring_t mem

    uint8_t bcf_type_shift[]

    # * VCF record *

    uint8_t BCF_BT_NULL
    uint8_t BCF_BT_INT8
    uint8_t BCF_BT_INT16
    uint8_t BCF_BT_INT32
    uint8_t BCF_BT_FLOAT
    uint8_t BCF_BT_CHAR

    uint8_t VCF_REF
    uint8_t VCF_SNP
    uint8_t VCF_MNP
    uint8_t VCF_INDEL
    uint8_t VCF_OTHER

    ctypedef struct variant_t:
        int type, n     # variant type and the number of bases affected, negative for deletions

    ctypedef struct bcf_fmt_t:
        int id             # id: numeric tag id, the corresponding string is bcf_hdr_t::id[BCF_DT_ID][$id].key
        int n, size, type  # n: number of values per-sample; size: number of bytes per-sample; type: one of BCF_BT_* types
        uint8_t *p         # same as vptr and vptr_* in bcf_info_t below
        uint32_t p_len
        uint32_t p_off
        uint8_t p_free

    ctypedef union bcf_info_union_t:
        int32_t i      # integer value
        float f        # float value

    ctypedef struct bcf_info_t:
        int key        # key: numeric tag id, the corresponding string is bcf_hdr_t::id[BCF_DT_ID][$key].key
        int type, len  # type: one of BCF_BT_* types; len: vector length, 1 for scalars

        # v1 union only set if $len==1; for easier access
        bcf_info_union_t v1
        uint8_t *vptr           # pointer to data array in bcf1_t->shared.s, excluding the size+type and tag id bytes
        uint32_t vptr_len       # length of the vptr block or, when set, of the vptr_mod block, excluding offset
        uint32_t vptr_off       # vptr offset, i.e., the size of the INFO key plus size+type bytes
        uint8_t  vptr_free      # indicates that vptr-vptr_off must be freed; set only when modified and the new
                                # data block is bigger than the original

    uint8_t BCF1_DIRTY_ID
    uint8_t BCF1_DIRTY_ALS
    uint8_t BCF1_DIRTY_FLT
    uint8_t BCF1_DIRTY_INF

    ctypedef struct bcf_dec_t:
        int m_fmt, m_info, m_id, m_als, m_allele, m_flt  # allocated size (high-water mark); do not change
        int n_flt           # Number of FILTER fields
        int *flt            # FILTER keys in the dictionary
        char *id            # ID
        char *als           # REF+ALT block (\0-seperated)
        char **allele       # allele[0] is the REF (allele[] pointers to the als block); all null terminated
        bcf_info_t *info    # INFO
        bcf_fmt_t *fmt      # FORMAT and individual sample
        variant_t *var      # $var and $var_type set only when set_variant_types called
        int n_var, var_type
        int shared_dirty    # if set, shared.s must be recreated on BCF output
        int indiv_dirty     # if set, indiv.s must be recreated on BCF output

    uint8_t BCF_ERR_CTG_UNDEF
    uint8_t BCF_ERR_TAG_UNDEF
    uint8_t BCF_ERR_NCOLS

    # The bcf1_t structure corresponds to one VCF/BCF line. Reading from VCF file
    # is slower because the string is first to be parsed, packed into BCF line
    # (done in vcf_parse), then unpacked into internal bcf1_t structure. If it
    # is known in advance that some of the fields will not be required (notably
    # the sample columns), parsing of these can be skipped by setting max_unpack
    # appropriately.
    # Similarly, it is fast to output a BCF line because the columns (kept in
    # shared.s, indiv.s, etc.) are written directly by bcf_write, whereas a VCF
    # line must be formatted in vcf_format.

    ctypedef struct bcf1_t:
        int32_t rid               # CHROM
        int32_t pos               # POS
        int32_t rlen              # length of REF
        float qual                # QUAL
        uint32_t n_info, n_allele
        uint32_t n_fmt, n_sample
        kstring_t shared, indiv
        bcf_dec_t d               # lazy evaluation: $d is not generated by bcf_read(), but by explicitly calling bcf_unpack()
        int max_unpack            # Set to BCF_UN_STR, BCF_UN_FLT, or BCF_UN_INFO to boost performance of vcf_parse when some of the fields won't be needed
        int unpacked              # remember what has been unpacked to allow calling bcf_unpack() repeatedly without redoing the work
        int unpack_size[3]        # the original block size of ID, REF+ALT and FILTER
        int errcode               # one of BCF_ERR_* codes

    ####### API #######

    # BCF and VCF I/O
    #
    # A note about naming conventions: htslib internally represents VCF
    # records as bcf1_t data structures, therefore most functions are
    # prefixed with bcf_. There are a few exceptions where the functions must
    # be aware of both BCF and VCF worlds, such as bcf_parse vs vcf_parse. In
    # these cases, functions prefixed with bcf_ are more general and work
    # with both BCF and VCF.

    # bcf_hdr_init() - create an empty BCF header.
    # @param mode    "r" or "w"
    #
    # When opened for writing, the mandatory fileFormat and
    # FILTER=PASS lines are added automatically.
    bcf_hdr_t *bcf_hdr_init(const char *mode)

    # Destroy a BCF header struct
    void bcf_hdr_destroy(bcf_hdr_t *h)

    # Initialize a bcf1_t object; equivalent to calloc(1, sizeof(bcf1_t))
    bcf1_t *bcf_init()

    # Deallocate a bcf1_t object
    void bcf_destroy(bcf1_t *v)

    # Same as bcf_destroy() but frees only the memory allocated by bcf1_t,
    # not the bcf1_t object itself.
    void bcf_empty(bcf1_t *v)

    # Make the bcf1_t object ready for next read. Intended mostly for
    # internal use, the user should rarely need to call this function
    # directly.
    void bcf_clear(bcf1_t *v)

    # Reads VCF or BCF header
    bcf_hdr_t *bcf_hdr_read(htsFile *fp)

    # bcf_hdr_set_samples() - for more efficient VCF parsing when only one/few samples are needed
    # @samples: samples to include or exclude from file or as a comma-separated string.
    #             LIST|FILE   .. select samples in list/file
    #             ^LIST|FILE  .. exclude samples from list/file
    #             -           .. include all samples
    #             NULL        .. exclude all samples
    # @is_file: @samples is a file (1) or a comma-separated list (1)
    #
    # The bottleneck of VCF reading is parsing of genotype fields. If the
    # reader knows in advance that only subset of samples is needed (possibly
    # no samples at all), the performance of bcf_read() can be significantly
    # improved by calling bcf_hdr_set_samples after bcf_hdr_read().
    # The function bcf_read() will subset the VCF/BCF records automatically
    # with the notable exception when reading records via bcf_itr_next().
    # In this case, bcf_subset_format() must be called explicitly, because
    # bcf_readrec() does not see the header.
    #
    # Returns 0 on success, -1 on error or a positive integer if the list
    # contains samples not present in the VCF header. In such a case, the
    # return value is the index of the offending sample.
    #
    int bcf_hdr_set_samples(bcf_hdr_t *hdr, const char *samples, int is_file)
    int bcf_subset_format(const bcf_hdr_t *hdr, bcf1_t *rec)

    # Writes VCF or BCF header
    int bcf_hdr_write(htsFile *fp, bcf_hdr_t *h)

    # Parse VCF line contained in kstring and populate the bcf1_t struct
    int vcf_parse(kstring_t *s, const bcf_hdr_t *h, bcf1_t *v)

    # The opposite of vcf_parse. It should rarely be called directly, see vcf_write
    int vcf_format(const bcf_hdr_t *h, const bcf1_t *v, kstring_t *s)

    # bcf_read() - read next VCF or BCF record
    #
    # Returns -1 on critical errors, 0 otherwise. On errors which are not
    # critical for reading, such as missing header definitions, v->errcode is
    # set to one of BCF_ERR* code and must be checked before calling
    # vcf_write().
    int bcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)

    # bcf_unpack() - unpack/decode a BCF record (fills the bcf1_t::d field)
    #
    # Note that bcf_unpack() must be called even when reading VCF. It is safe
    # to call the function repeatedly, it will not unpack the same field
    # twice.
    uint8_t BCF_UN_STR        # up to ALT inclusive
    uint8_t BCF_UN_FLT        # up to FILTER
    uint8_t BCF_UN_INFO       # up to INFO
    uint8_t BCF_UN_SHR        # all shared information
    uint8_t BCF_UN_FMT        # unpack format and each sample
    uint8_t BCF_UN_IND        # a synonymo of BCF_UN_FMT
    uint8_t BCF_UN_ALL        # everything

    int bcf_unpack(bcf1_t *b, int which)

    # bcf_dup() - create a copy of BCF record.
    #
    # Note that bcf_unpack() must be called on the returned copy as if it was
    # obtained from bcf_read(). Also note that bcf_dup() calls bcf_sync1(src)
    # internally to reflect any changes made by bcf_update_* functions.
    bcf1_t *bcf_dup(bcf1_t *src)
    bcf1_t *bcf_copy(bcf1_t *dst, bcf1_t *src)

    # bcf_write() - write one VCF or BCF record. The type is determined at the open() call.
    int bcf_write(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)

    # The following functions work only with VCFs and should rarely be called
    # directly. Usually one wants to use their bcf_* alternatives, which work
    # transparently with both VCFs and BCFs.
    bcf_hdr_t *vcf_hdr_read(htsFile *fp)
    int vcf_hdr_write(htsFile *fp, const bcf_hdr_t *h)
    int vcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
    int vcf_write(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)

    #************************************************************************
    # Header querying and manipulation routines
    #************************************************************************

    # Create a new header using the supplied template
    bcf_hdr_t *bcf_hdr_dup(const bcf_hdr_t *hdr)

    # Copy header lines from src to dst if not already present in dst. See also bcf_translate().
    # Returns 0 on success or sets a bit on error:
    #     1 .. conflicting definitions of tag length
    #     # todo
    int bcf_hdr_combine(bcf_hdr_t *dst, const bcf_hdr_t *src)

    # bcf_hdr_add_sample() - add a new sample.
    # @param sample:  sample name to be added
    int bcf_hdr_add_sample(bcf_hdr_t *hdr, const char *sample)

    # Read VCF header from a file and update the header
    int bcf_hdr_set(bcf_hdr_t *hdr, const char *fname)

    # Returns formatted header (newly allocated string) and its length,
    # excluding the terminating \0. If is_bcf parameter is unset, IDX
    # fields are discarded.
    char *bcf_hdr_fmt_text(const bcf_hdr_t *hdr, int is_bcf, int *len)

    # Append new VCF header line, returns 0 on success
    int bcf_hdr_append(bcf_hdr_t *h, const char *line)
    int bcf_hdr_printf(bcf_hdr_t *h, const char *format, ...)

    # VCF version, e.g. VCFv4.2
    const char *bcf_hdr_get_version(const bcf_hdr_t *hdr)
    void bcf_hdr_set_version(bcf_hdr_t *hdr, const char *version)

    # bcf_hdr_remove() - remove VCF header tag
    # @param type:      one of BCF_HL_*
    # @param key:       tag name
    void bcf_hdr_remove(bcf_hdr_t *h, int type, const char *key)

    # bcf_hdr_subset() - creates a new copy of the header removing unwanted samples
    # @param n:        number of samples to keep
    # @param samples:  names of the samples to keep
    # @param imap:     mapping from index in @samples to the sample index in the original file
    #
    # Sample names not present in h0 are ignored. The number of unmatched samples can be checked
    # by comparing n and bcf_hdr_nsamples(out_hdr).
    # This function can be used to reorder samples.
    # See also bcf_subset() which subsets individual records.
    #
    bcf_hdr_t *bcf_hdr_subset(const bcf_hdr_t *h0, int n, char *const* samples, int *imap)

    # Creates a list of sequence names. It is up to the caller to free the list (but not the sequence names)
    const char **bcf_hdr_seqnames(const bcf_hdr_t *h, int *nseqs)

    # Get number of samples
    int32_t bcf_hdr_nsamples(const bcf_hdr_t *h)

    # The following functions are for internal use and should rarely be called directly
    int bcf_hdr_parse(bcf_hdr_t *hdr, char *htxt)
    int bcf_hdr_sync(bcf_hdr_t *h)
    bcf_hrec_t *bcf_hdr_parse_line(const bcf_hdr_t *h, const char *line, int *len)
    void bcf_hrec_format(const bcf_hrec_t *hrec, kstring_t *str)
    int bcf_hdr_add_hrec(bcf_hdr_t *hdr, bcf_hrec_t *hrec)

    # bcf_hdr_get_hrec() - get header line info
    # @param type:  one of the BCF_HL_* types: FLT,INFO,FMT,CTG,STR,GEN
    # @param key:   the header key for generic lines (e.g. "fileformat"), any field
    #                 for structured lines, typically "ID".
    # @param value: the value which pairs with key. Can be be NULL for BCF_HL_GEN
    # @param str_class: the class of BCF_HL_STR line (e.g. "ALT" or "SAMPLE"), otherwise NULL
    #
    bcf_hrec_t *bcf_hdr_get_hrec(const bcf_hdr_t *hdr, int type, const char *key, const char *value, const char *str_class)
    bcf_hrec_t *bcf_hrec_dup(bcf_hrec_t *hrec)
    void bcf_hrec_add_key(bcf_hrec_t *hrec, const char *str, int len)
    void bcf_hrec_set_val(bcf_hrec_t *hrec, int i, const char *str, int len, int is_quoted)
    int bcf_hrec_find_key(bcf_hrec_t *hrec, const char *key)
    void hrec_add_idx(bcf_hrec_t *hrec, int idx)
    void bcf_hrec_destroy(bcf_hrec_t *hrec)

    #************************************************************************
    # Individual record querying and manipulation routines
    #************************************************************************

    # See the description of bcf_hdr_subset()
    int bcf_subset(const bcf_hdr_t *h, bcf1_t *v, int n, int *imap)

    # bcf_translate() - translate tags ids to be consistent with different header. This function
    #                   is useful when lines from multiple VCF need to be combined.
    # @dst_hdr:   the destination header, to be used in bcf_write(), see also bcf_hdr_combine()
    # @src_hdr:   the source header, used in bcf_read()
    # @src_line:  line obtained by bcf_read()
    int bcf_translate(const bcf_hdr_t *dst_hdr, bcf_hdr_t *src_hdr, bcf1_t *src_line)

    # bcf_get_variant_type[s]()  - returns one of VCF_REF, VCF_SNP, etc
    int bcf_get_variant_types(bcf1_t *rec)
    int bcf_get_variant_type(bcf1_t *rec, int ith_allele)
    int bcf_is_snp(bcf1_t *v)

    # bcf_update_filter() - sets the FILTER column
    # @flt_ids:  The filter IDs to set, numeric IDs returned by bcf_id2int(hdr, BCF_DT_ID, "PASS")
    # @n:        Number of filters. If n==0, all filters are removed
    int bcf_update_filter(const bcf_hdr_t *hdr, bcf1_t *line, int *flt_ids, int n)

    # bcf_add_filter() - adds to the FILTER column
    # @flt_id:   filter ID to add, numeric ID returned by bcf_id2int(hdr, BCF_DT_ID, "PASS")
    #
    # If flt_id is PASS, all existing filters are removed first. If other than PASS, existing PASS is removed.
    int bcf_add_filter(const bcf_hdr_t *hdr, bcf1_t *line, int flt_id)

    # bcf_remove_filter() - removes from the FILTER column
    # @flt_id:   filter ID to remove, numeric ID returned by bcf_id2int(hdr, BCF_DT_ID, "PASS")
    # @pass:     when set to 1 and no filters are present, set to PASS
    int bcf_remove_filter(const bcf_hdr_t *hdr, bcf1_t *line, int flt_id, int set_pass)

    # Returns 1 if present, 0 if absent, or -1 if filter does not exist. "PASS" and "." can be used interchangeably.
    int bcf_has_filter(const bcf_hdr_t *hdr, bcf1_t *line, char *filter)

    # bcf_update_alleles() and bcf_update_alleles_str() - update REF and ALT column
    # @alleles:           Array of alleles
    # @nals:              Number of alleles
    # @alleles_string:    Comma-separated alleles, starting with the REF allele
    #
    # Not that in order for indexing to work correctly in presence of INFO/END tag,
    # the length of reference allele (line->rlen) must be set explicitly by the caller,
    # or otherwise, if rlen is zero, strlen(line->d.allele[0]) is used to set the length
    # on bcf_write().
    #
    int bcf_update_alleles(const bcf_hdr_t *hdr, bcf1_t *line, const char **alleles, int nals)
    int bcf_update_alleles_str(const bcf_hdr_t *hdr, bcf1_t *line, const char *alleles_string)
    int bcf_update_id(const bcf_hdr_t *hdr, bcf1_t *line, const char *id)

    # bcf_update_info_*() - functions for updating INFO fields
    # @hdr:       the BCF header
    # @line:      VCF line to be edited
    # @key:       the INFO tag to be updated
    # @values:    pointer to the array of values. Pass NULL to remove the tag.
    # @n:         number of values in the array. When set to 0, the INFO tag is removed
    #
    # The @string in bcf_update_info_flag() is optional, @n indicates whether
    # the flag is set or removed.
    #
    # Returns 0 on success or negative value on error.
    #
    int bcf_update_info_int32(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const int32_t *values, int n)
    int bcf_update_info_float(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const float *values, int n)
    int bcf_update_info_flag(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const char *values, int n)
    int bcf_update_info_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const char *values, int n)
    int bcf_update_info(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type)

    # bcf_update_format_*() - functions for updating FORMAT fields
    # @values:    pointer to the array of values, the same number of elements
    #             is expected for each sample. Missing values must be padded
    #             with bcf_*_missing or bcf_*_vector_end values.
    # @n:         number of values in the array. If n==0, existing tag is removed.
    #
    # The function bcf_update_format_string() is a higher-level (slower) variant of
    # bcf_update_format_char(). The former accepts array of \0-terminated strings
    # whereas the latter requires that the strings are collapsed into a single array
    # of fixed-length strings. In case of strings with variable length, shorter strings
    # can be \0-padded. Note that the collapsed strings passed to bcf_update_format_char()
    # are not \0-terminated.
    #
    # Returns 0 on success or negative value on error.
    #
    int bcf_update_format_int32(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const int32_t *values, int n)
    int bcf_update_format_float(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const float *values, int n)
    int bcf_update_format_char(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const char *values, int n)
    int bcf_update_genotypes(const bcf_hdr_t *hdr, bcf1_t *line, const int32_t *values, int n)
    int bcf_update_format_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const char **values, int n)
    int bcf_update_format(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type)

    # Macros for setting genotypes correctly, for use with bcf_update_genotypes only; idx corresponds
    # to VCF's GT (1-based index to ALT or 0 for the reference allele) and val is the opposite, obtained
    # from bcf_get_genotypes() below.
    uint32_t bcf_gt_phased(uint32_t idx)
    uint32_t bcf_gt_unphased(uint32_t idx)
    uint32_t bcf_gt_missing
    uint32_t bcf_gt_is_missing(uint32_t val)
    uint32_t bcf_gt_is_phased(uint32_t idx)
    uint32_t bcf_gt_allele(uint32_t val)

    # Conversion between alleles indexes to Number=G genotype index (assuming diploid, all 0-based)
    uint32_t bcf_alleles2gt(uint32_t a, uint32_t b)
    void bcf_gt2alleles(int igt, int *a, int *b)

    # bcf_get_fmt() - returns pointer to FORMAT's field data
    # @header: for access to BCF_DT_ID dictionary
    # @line:   VCF line obtained from vcf_parse1
    # @fmt:    one of GT,PL,...
    #
    # Returns bcf_fmt_t* if the call succeeded, or returns NULL when the field
    # is not available.
    #
    bcf_fmt_t *bcf_get_fmt(const bcf_hdr_t *hdr, bcf1_t *line, const char *key)
    bcf_info_t *bcf_get_info(const bcf_hdr_t *hdr, bcf1_t *line, const char *key)

    # bcf_get_*_id() - returns pointer to FORMAT/INFO field data given the header index instead of the string ID
    # @line: VCF line obtained from vcf_parse1
    # @id:  The header index for the tag, obtained from bcf_hdr_id2int()
    #
    # Returns bcf_fmt_t* / bcf_info_t*. These functions do not check if the index is valid
    # as their goal is to avoid the header lookup.
    #
    bcf_fmt_t *bcf_get_fmt_id(bcf1_t *line, const int id)
    bcf_info_t *bcf_get_info_id(bcf1_t *line, const int id)

    # bcf_get_info_*() - get INFO values, integers or floats
    # @hdr:       BCF header
    # @line:      BCF record
    # @tag:       INFO tag to retrieve
    # @dst:       *dst is pointer to a memory location, can point to NULL
    # @ndst:      pointer to the size of allocated memory
    #
    # Returns negative value on error or the number of written values on
    # success. bcf_get_info_string() returns on success the number of
    # characters written excluding the null-terminating byte. bcf_get_info_flag()
    # returns 1 when flag is set or 0 if not.
    #
    # List of return codes:
    #     -1 .. no such INFO tag defined in the header
    #     -2 .. clash between types defined in the header and encountered in the VCF record
    #     -3 .. tag is not present in the VCF record
    #
    int bcf_get_info_int32(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, int32_t **dst, int *ndst)
    int bcf_get_info_float(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, float **dst, int *ndst)
    int bcf_get_info_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, char **dst, int *ndst)
    int bcf_get_info_flag(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, int **dst, int *ndst)
    int bcf_get_info_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type)

    # bcf_get_format_*() - same as bcf_get_info*() above
    #
    # The function bcf_get_format_string() is a higher-level (slower) variant of bcf_get_format_char().
    # see the description of bcf_update_format_string() and bcf_update_format_char() above.
    # Unlike other bcf_get_format__*() functions, bcf_get_format_string() allocates two arrays:
    # a single block of \0-terminated strings collapsed into a single array and an array of pointers
    # to these strings. Both arrays must be cleaned by the user.
    #
    # Returns negative value on error or the number of written values on success.
    #
    # Example:
    #     int ndst = 0; char **dst = NULL
    #     if ( bcf_get_format_string(hdr, line, "XX", &dst, &ndst) > 0 )
    #         for (i=0; i<bcf_hdr_nsamples(hdr); i++) printf("%s\n", dst[i])
    #     free(dst[0]); free(dst)
    #
    # Example:
    #     int ngt, *gt_arr = NULL, ngt_arr = 0
    #     ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr)
    #
    int bcf_get_format_int32(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, int32_t **dst, int *ndst)
    int bcf_get_format_float(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, float **dst, int *ndst)
    int bcf_get_format_char(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, char **dst, int *ndst)
    int bcf_get_genotypes(const bcf_hdr_t *hdr, bcf1_t *line, int **dst, int *ndst)
    int bcf_get_format_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, char ***dst, int *ndst)
    int bcf_get_format_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type)

    #************************************************************************
    # Helper functions
    #************************************************************************

    #
    # bcf_hdr_id2int() - Translates string into numeric ID
    # bcf_hdr_int2id() - Translates numeric ID into string
    # @type:     one of BCF_DT_ID, BCF_DT_CTG, BCF_DT_SAMPLE
    # @id:       tag name, such as: PL, DP, GT, etc.
    #
    # Returns -1 if string is not in dictionary, otherwise numeric ID which identifies
    # fields in BCF records.
    #
    int bcf_hdr_id2int(const bcf_hdr_t *hdr, int type, const char *id)
    const char *bcf_hdr_int2id(const bcf_hdr_t *hdr, int type, int int_id)

    # bcf_hdr_name2id() - Translates sequence names (chromosomes) into numeric ID
    # bcf_hdr_id2name() - Translates numeric ID to sequence name
    #
    int bcf_hdr_name2id(const bcf_hdr_t *hdr, const char *id)
    const char *bcf_hdr_id2name(const bcf_hdr_t *hdr, int rid)
    const char *bcf_seqname(const bcf_hdr_t *hdr, bcf1_t *rec)

    #
    # bcf_hdr_id2*() - Macros for accessing bcf_idinfo_t
    # @type:      one of BCF_HL_FLT, BCF_HL_INFO, BCF_HL_FMT
    # @int_id:    return value of bcf_id2int, must be >=0
    #
    # The returned values are:
    #    bcf_hdr_id2length   ..  whether the number of values is fixed or variable, one of BCF_VL_*
    #    bcf_hdr_id2number   ..  the number of values, 0xfffff for variable length fields
    #    bcf_hdr_id2type     ..  the field type, one of BCF_HT_*
    #    bcf_hdr_id2coltype  ..  the column type, one of BCF_HL_*
    #
    # Notes: Prior to using the macros, the presence of the info should be
    # tested with bcf_hdr_idinfo_exists().
    #
    int bcf_hdr_id2length(const bcf_hdr_t *hdr, int type, int int_id)
    int bcf_hdr_id2number(const bcf_hdr_t *hdr, int type, int int_id)
    int bcf_hdr_id2type(const bcf_hdr_t *hdr, int type, int int_id)
    int bcf_hdr_id2coltype(const bcf_hdr_t *hdr, int type, int int_id)
    int bcf_hdr_idinfo_exists(const bcf_hdr_t *hdr, int type, int int_id)
    bcf_hrec_t *bcf_hdr_id2hrec(const bcf_hdr_t *hdr, int type, int col_type, int int_id)

    void bcf_fmt_array(kstring_t *s, int n, int type, void *data)
    uint8_t *bcf_fmt_sized_array(kstring_t *s, uint8_t *ptr)

    void bcf_enc_vchar(kstring_t *s, int l, const char *a)
    void bcf_enc_vint(kstring_t *s, int n, int32_t *a, int wsize)
    void bcf_enc_vfloat(kstring_t *s, int n, float *a)

    #************************************************************************
    # BCF index
    #
    # Note that these functions work with BCFs only. See synced_bcf_reader.h
    # which provides (amongst other things) an API to work transparently with
    # both indexed BCFs and VCFs.
    #************************************************************************

    int bcf_index_build(const char *fn, int min_shift)

    #*******************
    # Typed value I/O *
    #******************

    # Note that in contrast with BCFv2.1 specification, HTSlib implementation
    # allows missing values in vectors. For integer types, the values 0x80,
    # 0x8000, 0x80000000 are interpreted as missing values and 0x81, 0x8001,
    # 0x80000001 as end-of-vector indicators.  Similarly for floats, the value of
    # 0x7F800001 is interpreted as a missing value and 0x7F800002 as an
    # end-of-vector indicator.
    # Note that the end-of-vector byte is not part of the vector.

    # This trial BCF version (v2.2) is compatible with the VCF specification and
    # enables to handle correctly vectors with different ploidy in presence of
    # missing values.

    int32_t bcf_int8_vector_end
    int32_t bcf_int16_vector_end
    int32_t bcf_int32_vector_end
    int32_t bcf_str_vector_end
    int32_t bcf_int8_missing
    int32_t bcf_int16_missing
    int32_t bcf_int32_missing
    int32_t bcf_str_missing

    uint32_t bcf_float_vector_end
    uint32_t bcf_float_missing

    void bcf_float_set(float *ptr, uint32_t value)
    void bcf_float_set_vector_end(float *x)
    void bcf_float_set_missing(float *x)

    int bcf_float_is_missing(float f)
    int bcf_float_is_vector_end(float f)
    void bcf_format_gt(bcf_fmt_t *fmt, int isample, kstring_t *str)
    void bcf_enc_size(kstring_t *s, int size, int type)
    int bcf_enc_inttype(long x)
    void bcf_enc_int1(kstring_t *s, int32_t x)
    int32_t bcf_dec_int1(const uint8_t *p, int type, uint8_t **q)
    int32_t bcf_dec_typed_int1(const uint8_t *p, uint8_t **q)
    int32_t bcf_dec_size(const uint8_t *p, uint8_t **q, int *type)

    # These trivial wrappers are defined only for consistency with other parts of htslib
    bcf1_t *bcf_init1()
    int bcf_read1(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
    int vcf_read1(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
    int bcf_write1(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
    int vcf_write1(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
    void bcf_destroy1(bcf1_t *v)
    int vcf_parse1(kstring_t *s, const bcf_hdr_t *h, bcf1_t *v)
    void bcf_clear1(bcf1_t *v)
    int vcf_format1(const bcf_hdr_t *h, const bcf1_t *v, kstring_t *s)

    # Other nice wrappers
    void bcf_itr_destroy(hts_itr_t *iter)
    hts_itr_t *bcf_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end)
    hts_itr_t *bcf_itr_querys(const hts_idx_t *idx, const bcf_hdr_t *hdr, char *s)
    int bcf_itr_next(htsFile *fp, hts_itr_t *iter, void *r)
    hts_idx_t *bcf_index_load(const char *fn)
    const char **bcf_index_seqnames(const hts_idx_t *idx, const bcf_hdr_t *hdr, int *nptr)
