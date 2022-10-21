# cython: language_level=3
from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup
from libc.stdio cimport FILE, printf
from posix.types cimport off_t

cdef extern from "Python.h":
   FILE* PyFile_AsFile(object)


# cython does not wrap stdarg
cdef extern from "stdarg.h":
    ctypedef struct va_list:
        pass

   
cdef extern from "htslib/kstring.h" nogil:
    ctypedef struct kstring_t:
        size_t l, m
        char *s

    int kputc(int c, kstring_t *s)
    int kputw(int c, kstring_t *s)
    int kputl(long c, kstring_t *s)
    int ksprintf(kstring_t *s, const char *fmt, ...)


cdef extern from "htslib_util.h" nogil:
    int hts_set_verbosity(int verbosity)
    int hts_get_verbosity()

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
    hFILE *hopen(const char *filename, const char *mode, ...)

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

    # Read from the stream until the delimiter, up to a maximum length
    #    @param buffer  The buffer into which bytes will be written
    #    @param size    The size of the buffer
    #    @param delim   The delimiter (interpreted as an `unsigned char`)
    #    @param fp      The file stream
    #    @return  The number of bytes read, or negative on error.
    #    @since   1.4
    #
    # Bytes will be read into the buffer up to and including a delimiter, until
    # EOF is reached, or _size-1_ bytes have been written, whichever comes first.
    # The string will then be terminated with a NUL byte (`\0`).
    ssize_t hgetdelim(char *buffer, size_t size, int delim, hFILE *fp)

    # Read a line from the stream, up to a maximum length
    #    @param buffer  The buffer into which bytes will be written
    #    @param size    The size of the buffer
    #    @param fp      The file stream
    #    @return  The number of bytes read, or negative on error.
    #    @since   1.4
    #
    # Specialization of hgetdelim() for a `\n` delimiter.
    ssize_t hgetln(char *buffer, size_t size, hFILE *fp)

    # Read a line from the stream, up to a maximum length
    #    @param buffer  The buffer into which bytes will be written
    #    @param size    The size of the buffer (must be > 1 to be useful)
    #    @param fp      The file stream
    #    @return  _buffer_ on success, or `NULL` if an error occurred.
    #    @since   1.4
    #
    # This function can be used as a replacement for `fgets(3)`, or together with
    # kstring's `kgetline()` to read arbitrarily-long lines into a _kstring_t_.
    char *hgets(char *buffer, int size, hFILE *fp)

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
        unsigned           errcode
        unsigned           is_write
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
    #  No interpretation of the value should be made, other than a subsequent
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
    #  @param delim  delimiter
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


    # Compress a single BGZF block.
    #
    # @param dst    output buffer (must have size >= BGZF_MAX_BLOCK_SIZE)
    # @param dlen   size of output buffer; updated on return to the number
    #               of bytes actually written to dst
    # @param src    buffer to be compressed
    # @param slen   size of data to compress (must be <= BGZF_BLOCK_SIZE)
    # @param level  compression level
    # @return       0 on success and negative on error
    #
    int bgzf_compress(void *dst, size_t *dlen, const void *src, size_t slen, int level)

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

    union FilePointerUnion:
        BGZF    *bgzf
        cram_fd *cram
        hFILE   *hfile
        void    *voidp

    enum htsFormatCategory:
        unknown_category
        sequence_data    # Sequence data -- SAM, BAM, CRAM, etc
        variant_data     # Variant calling data -- VCF, BCF, etc
        index_file       # Index file associated with some data file
        region_list      # Coordinate intervals or regions -- BED, etc
        category_maximum

    enum htsExactFormat:
        unknown_format
        binary_format
        text_format
        sam, bam, bai, cram, crai, vcf, bcf, csi, gzi, tbi, bed
        format_maximum

    enum htsCompression:
        no_compression, gzip, bgzf, custom
        compression_maximum

    cdef enum hts_fmt_option:
        CRAM_OPT_DECODE_MD,
        CRAM_OPT_PREFIX,
        CRAM_OPT_VERBOSITY,
        CRAM_OPT_SEQS_PER_SLICE,
        CRAM_OPT_SLICES_PER_CONTAINER,
        CRAM_OPT_RANGE,
        CRAM_OPT_VERSION,
        CRAM_OPT_EMBED_REF,
        CRAM_OPT_IGNORE_MD5,
        CRAM_OPT_REFERENCE,
        CRAM_OPT_MULTI_SEQ_PER_SLICE,
        CRAM_OPT_NO_REF,
        CRAM_OPT_USE_BZIP2,
        CRAM_OPT_SHARED_REF,
        CRAM_OPT_NTHREADS,
        CRAM_OPT_THREAD_POOL,
        CRAM_OPT_USE_LZMA,
        CRAM_OPT_USE_RANS,
        CRAM_OPT_REQUIRED_FIELDS,
        HTS_OPT_COMPRESSION_LEVEL,
        HTS_OPT_NTHREADS,

    ctypedef struct htsVersion:
        short major, minor

    ctypedef struct htsFormat:
        htsFormatCategory category
        htsExactFormat    format
        htsVersion        version
        htsCompression    compression
        short             compression_level
        void              *specific  

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

    cdef union hts_opt_val_union:
        int i
        char *s

    ctypedef struct hts_opt:
        char *arg
        hts_fmt_option opt
        hts_opt_val_union val
        void *next

    # @abstract Parses arg and appends it to the option list.
    # @return   0 on success and -1 on failure
    int hts_opt_add(hts_opt **opts, const char *c_arg)

    # @abstract Applies an hts_opt option list to a given htsFile.
    # @return   0 on success and -1 on failure
    int hts_opt_apply(htsFile *fp, hts_opt *opts)

    # @abstract Frees an hts_opt list.
    void hts_opt_free(hts_opt *opts)

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
    # @param mode     Mode matching / [rwa][bceguxz0-9]* /
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
    #     and with non-format option letters (for any of 'r'/'w'/'a'):
    #       e  close the file on exec(2) (opens with O_CLOEXEC, where supported)
    #       x  create the file exclusively (opens with O_EXCL, where supported)
    #     Note that there is a distinction between 'u' and '0': the first yields
    #     plain uncompressed output whereas the latter outputs uncompressed data
    #     wrapped in the zlib format.
    # @example
    #     [rw]b  .. compressed BCF, BAM, FAI
    #     [rw]bu .. uncompressed BCF
    #     [rw]z  .. compressed VCF
    #     [rw]   .. uncompressed VCF
    htsFile *hts_open(const char *fn, const char *mode)

    # @abstract       Open a SAM/BAM/CRAM/VCF/BCF/etc file
    # @param fn       The file name or "-" for stdin/stdout
    # @param mode     Open mode, as per hts_open()
    # @param fmt      Optional format specific parameters
    # @discussion
    #     See hts_open() for description of fn and mode.
    #     // TODO Update documentation for s/opts/fmt/
    #     Opts contains a format string (sam, bam, cram, vcf, bcf) which will,
    #     if defined, override mode.  Opts also contains a linked list of hts_opt
    #     structures to apply to the open file handle.  These can contain things
    #     like pointers to the reference or information on compression levels,
    #     block sizes, etc.
    htsFile *hts_open_format(const char *fn, const char *mode, const htsFormat *fmt)

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

    # @ abstract      Returns a string containing the file format extension.
    # @ param format  Format structure containing the file type.
    # @ return        A string ("sam", "bam", etc) or "?" for unknown formats.
    const char *hts_format_file_extension(const htsFormat *format)

    # @abstract  Sets a specified CRAM option on the open file handle.
    # @param fp  The file handle open the open file.
    # @param opt The CRAM_OPT_* option.
    # @param ... Optional arguments, dependent on the option used.
    # @return    0 for success, or negative if an error occurred.
    int hts_set_opt(htsFile *fp, hts_fmt_option opt, ...)

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

    #### Save an index to a file
    #    @param idx  Index to be written
    #    @param fn   Input BAM/BCF/etc filename, to which .bai/.csi/etc will be added
    #    @param fmt  One of the HTS_FMT_* index formats
    #    @return  0 if successful, or negative if an error occurred.
    int hts_idx_save(const hts_idx_t *idx, const char *fn, int fmt)

    #### Save an index to a specific file
    #    @param idx    Index to be written
    #    @param fn     Input BAM/BCF/etc filename
    #    @param fnidx  Output filename, or NULL to add .bai/.csi/etc to @a fn
    #    @param fmt    One of the HTS_FMT_* index formats
    #    @return  0 if successful, or negative if an error occurred.
    int hts_idx_save_as(const hts_idx_t *idx, const char *fn, const char *fnidx, int fmt)

    #### Load an index file
    #    @param fn   BAM/BCF/etc filename, to which .bai/.csi/etc will be added or
    #                the extension substituted, to search for an existing index file
    #    @param fmt  One of the HTS_FMT_* index formats
    #    @return  The index, or NULL if an error occurred.
    hts_idx_t *hts_idx_load(const char *fn, int fmt)

    #### Load a specific index file
    #    @param fn     Input BAM/BCF/etc filename
    #    @param fnidx  The input index filename
    #    @return  The index, or NULL if an error occurred.
    hts_idx_t *hts_idx_load2(const char *fn, const char *fnidx)

    uint8_t *hts_idx_get_meta(hts_idx_t *idx, uint32_t *l_meta)
    void hts_idx_set_meta(hts_idx_t *idx, int l_meta, uint8_t *meta, int is_copy)

    int hts_idx_get_stat(const hts_idx_t* idx, int tid,
                         uint64_t* mapped, uint64_t* unmapped)

    uint64_t hts_idx_get_n_no_coor(const hts_idx_t* idx)

    int HTS_PARSE_THOUSANDS_SEP  # Ignore ',' separators within numbers

    # Parse a numeric string
    #    The number may be expressed in scientific notation, and optionally may
    #    contain commas in the integer part (before any decimal point or E notation).
    #    @param str     String to be parsed
    #    @param strend  If non-NULL, set on return to point to the first character
    #                   in @a str after those forming the parsed number
    #    @param flags   Or'ed-together combination of HTS_PARSE_* flags
    #    @return  Converted value of the parsed number.
    #
    #    When @a strend is NULL, a warning will be printed (if hts_verbose is 2
    #    or more) if there are any trailing characters after the number.
    long long hts_parse_decimal(const char *str, char **strend, int flags)

    # Parse a "CHR:START-END"-style region string
    #    @param str  String to be parsed
    #    @param beg  Set on return to the 0-based start of the region
    #    @param end  Set on return to the 1-based end of the region
    #    @return  Pointer to the colon or '\0' after the reference sequence name,
    #             or NULL if @a str could not be parsed.
    const char *hts_parse_reg(const char *str, int *beg, int *end)

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

    # /***************************
    #  * Revised MAQ error model *
    #  ***************************/

    ctypedef struct errmod_t

    errmod_t *errmod_init(double depcorr)
    void errmod_destroy(errmod_t *em)

    # /*
    #     n: number of bases
    #     m: maximum base
    #     bases[i]: qual:6, strand:1, base:4
    #     q[i*m+j]: phred-scaled likelihood of (i,j)
    #  */
    int errmod_cal(const errmod_t *em, int n, int m, uint16_t *bases, float *Probabilistic)

    # /*****************************************
    #  * q banded glocal alignment *
    #  *****************************************/

    ctypedef struct probaln_par_t:
        float d, e
        int bw

    int probaln_glocal(const uint8_t *ref,
                       int l_ref,
                       const uint8_t *query,
                       int l_query, const uint8_t *iqual,
                       const probaln_par_t *c,
                       int *state, uint8_t *q)

    # /**********************
    #  * MD5 implementation *
    #  **********************/

    ctypedef struct hts_md5_context

    # /*! @abstract   Initialises an MD5 context.
    #  *  @discussion
    #  *    The expected use is to allocate an hts_md5_context using
    #  *    hts_md5_init().  This pointer is then passed into one or more calls
    #  *    of hts_md5_update() to compute successive internal portions of the
    #  *    MD5 sum, which can then be externalised as a full 16-byte MD5sum
    #  *    calculation by calling hts_md5_final().  This can then be turned
    #  *    into ASCII via hts_md5_hex().
    #  *
    #  *    To dealloate any resources created by hts_md5_init() call the
    #  *    hts_md5_destroy() function.
    #  *
    #  *  @return     hts_md5_context pointer on success, NULL otherwise.
    #  */
    hts_md5_context *hts_md5_init()

    # /*! @abstract Updates the context with the MD5 of the data. */
    void hts_md5_update(hts_md5_context *ctx, const void *data, unsigned long size)

    # /*! @abstract Computes the final 128-bit MD5 hash from the given context */
    void hts_md5_final(unsigned char *digest, hts_md5_context *ctx)

    # /*! @abstract Resets an md5_context to the initial state, as returned
    #  *            by hts_md5_init().
    #  */
    void hts_md5_reset(hts_md5_context *ctx)

    # /*! @abstract Converts a 128-bit MD5 hash into a 33-byte nul-termninated
    #  *            hex string.
    #  */
    void hts_md5_hex(char *hex, const unsigned char *digest)

    # /*! @abstract Deallocates any memory allocated by hts_md5_init. */
    void hts_md5_destroy(hts_md5_context *ctx)

    int hts_reg2bin(int64_t beg, int64_t end, int min_shift, int n_lvls)
    int hts_bin_bot(int bin, int n_lvls)

    # * Endianness *
    int ed_is_big()
    uint16_t ed_swap_2(uint16_t v)
    void *ed_swap_2p(void *x)
    uint32_t ed_swap_4(uint32_t v)
    void *ed_swap_4p(void *x)
    uint64_t ed_swap_8(uint64_t v)
    void *ed_swap_8p(void *x)


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
        uint8_t unused1
        uint8_t l_extranul
        uint32_t n_cigar
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
        int l_data
        uint32_t m_data
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

    # Load/build .csi or .bai BAM index file.  Does not work with CRAM.
    # It is recommended to use the sam_index_* functions below instead.
    hts_idx_t *bam_index_load(const char *fn)
    int bam_index_build(const char *fn, int min_shift)

    # Load a BAM (.csi or .bai) or CRAM (.crai) index file
    # @param fp  File handle of the data file whose index is being opened
    # @param fn  BAM/CRAM/etc filename to search alongside for the index file
    # @return  The index, or NULL if an error occurred.
    hts_idx_t *sam_index_load(htsFile *fp, const char *fn)

    # Load a specific BAM (.csi or .bai) or CRAM (.crai) index file
    # @param fp     File handle of the data file whose index is being opened
    # @param fn     BAM/CRAM/etc data file filename
    # @param fnidx  Index filename, or NULL to search alongside @a fn
    # @return  The index, or NULL if an error occurred.
    hts_idx_t *sam_index_load2(htsFile *fp, const char *fn, const char *fnidx)

    # Generate and save an index file
    # @param fn        Input BAM/etc filename, to which .csi/etc will be added
    # @param min_shift Positive to generate CSI, or 0 to generate BAI
    # @return  0 if successful, or negative if an error occurred (usually -1; or
    #         -2: opening fn failed; -3: format not indexable)
    int sam_index_build(const char *fn, int min_shift)

    # Generate and save an index to a specific file
    # @param fn        Input BAM/CRAM/etc filename
    # @param fnidx     Output filename, or NULL to add .bai/.csi/etc to @a fn
    # @param min_shift Positive to generate CSI, or 0 to generate BAI
    # @return  0 if successful, or negative if an error occurred.
    int sam_index_build2(const char *fn, const char *fnidx, int min_shift)

    void sam_itr_destroy(hts_itr_t *iter)
    hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end)
    hts_itr_t *sam_itr_querys(const hts_idx_t *idx, bam_hdr_t *hdr, const char *region)
    int sam_itr_next(htsFile *htsfp, hts_itr_t *itr, void *r)

    #***************
    #*** SAM I/O ***
    #***************

    htsFile *sam_open(const char *fn, const char *mode)
    htsFile *sam_open_format(const char *fn, const char *mode, const htsFormat *fmt)
    int sam_close(htsFile *fp)

    int sam_open_mode(char *mode, const char *fn, const char *format)

    # A version of sam_open_mode that can handle ,key=value options.
    # The format string is allocated and returned, to be freed by the caller.
    # Prefix should be "r" or "w",
    char *sam_open_mode_opts(const char *fn, const char *mode, const char *format)

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
    int64_t  bam_aux2i(const uint8_t *s)
    double   bam_aux2f(const uint8_t *s)
    char     bam_aux2A(const uint8_t *s)
    char    *bam_aux2Z(const uint8_t *s)

    void bam_aux_append(bam1_t *b, const char *tag, char type, int len, uint8_t *data)
    int bam_aux_del(bam1_t *b, uint8_t *s)

    #**************************
    #*** Pileup and Mpileup ***
    #**************************

    #  @abstract Generic pileup 'client data'.
    #  @discussion The pileup iterator allows setting a constructor and
    #  destructor function, which will be called every time a sequence is
    #  fetched and discarded.  This permits caching of per-sequence data in
    #  a tidy manner during the pileup process.  This union is the cached
    #  data to be manipulated by the "client" (the caller of pileup).
    # 
    union bam_pileup_cd:
        void *p
        int64_t i
        double f

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
        bam_pileup_cd cd

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
    void bam_mplp_reset(bam_mplp_t iter)
    void bam_mplp_constructor(bam_mplp_t iter,
          		      int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd))
    void bam_mplp_destructor(bam_mplp_t iter,
			     int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd))
    
    # Added by AH
    # ctypedef bam_pileup1_t * const_bam_pileup1_t_ptr "const bam_pileup1_t *"




    # // ---------------------------
    # // Base modification retrieval

    # /*! @typedef
    #  @abstract Holds a single base modification.
    #  @field modified_base     The short base code (m, h, etc) or -ChEBI (negative)
    #  @field canonical_base    The canonical base referred to in the MM tag.
    #                           One of A, C, G, T or N.  Note this may not be the
    #                           explicit base recorded in the SEQ column (esp. if N).
    #  @field strand            0 or 1, indicating + or - strand from MM tag.
    #  @field qual              Quality code (256*probability), or -1 if unknown

    #  @discussion
    #  Note this doesn't hold any location data or information on which other
    #  modifications may be possible at this site.
    ctypedef struct hts_base_mod:
        int modified_base
        int canonical_base
        int strand
        int qual

    # /// Allocates an hts_base_mode_state.
    # /**
    # * @return An hts_base_mode_state pointer on success,
    # *         NULL on failure.
    # *
    # * This just allocates the memory.  The initialisation of the contents is
    # * done using bam_parse_basemod.  Successive calls may be made to that
    # * without the need to free and allocate a new state.
    # *
    # * The state be destroyed using the hts_base_mode_state_free function.
    # */
    ctypedef struct hts_base_mod_state 
    hts_base_mod_state *hts_base_mod_state_alloc()


    # /// Destroys an  hts_base_mode_state.
    # /**
    # * @param state    The base modification state pointer.
    # *
    # * The should have previously been created by hts_base_mode_state_alloc.
    # */
    void hts_base_mod_state_free(hts_base_mod_state *state)

    # /// Parses the Mm and Ml tags out of a bam record.
    # /**
    # * @param b        BAM alignment record
    # * @param state    The base modification state pointer.
    # * @return 0 on success,
    # *         -1 on failure.
    # *
    # * This fills out the contents of the modification state, resetting the
    # * iterator location to the first sequence base.
    # */
    int bam_parse_basemod(const bam1_t *b, hts_base_mod_state *state)

    # /// Finds the next location containing base modifications and returns them
    # /**
    # * @param b        BAM alignment record
    # * @param state    The base modification state pointer.
    # * @param mods     A supplied array for returning base modifications
    # * @param n_mods   The size of the mods array
    # * @return The number of modifications found on success,
    # *         0 if no more modifications are present,
    # *         -1 on failure.
    # *
    # * Unlike bam_mods_at_next_pos this skips ahead to the next site
    # * with modifications.
    # *
    # * If more than n_mods modifications are found, the total found is returned.
    # * Note this means the caller needs to check whether this is higher than
    # * n_mods.
    # */

    int bam_next_basemod(const bam1_t *b, hts_base_mod_state *state,hts_base_mod *mods, int n_mods, int *pos)

    # ***********************************
    # * BAQ calculation and realignment *
    # ***********************************/
    int sam_cap_mapq(bam1_t *b, const char *ref, int ref_len, int thres)
    int sam_prob_realn(bam1_t *b, const char *ref, int ref_len, int flag)


cdef extern from "htslib/faidx.h" nogil:

    ctypedef struct faidx_t:
       pass

    # /// Build index for a FASTA or bgzip-compressed FASTA file.
    # /**  @param  fn  FASTA file name
    # @param  fnfai Name of .fai file to build.
    # @param  fngzi Name of .gzi file to build (if fn is bgzip-compressed).
    # @return     0 on success; or -1 on failure

    # If fnfai is NULL, ".fai" will be appended to fn to make the FAI file name.
    # If fngzi is NULL, ".gzi" will be appended to fn for the GZI file.  The GZI
    # file will only be built if fn is bgzip-compressed.
    # */
    int fai_build3(const char *fn,
                   const char *fnfai,
                   const char *fngzi)

    # /// Build index for a FASTA or bgzip-compressed FASTA file.
    # /** @param  fn  FASTA file name
    # @return     0 on success; or -1 on failure
    #
    # File "fn.fai" will be generated.  This function is equivalent to
    # fai_build3(fn, NULL, NULL);
    # */
    int fai_build(char *fn)

    # /// Destroy a faidx_t struct
    void fai_destroy(faidx_t *fai)

    # /// Load FASTA indexes.
    # /** @param  fn  File name of the FASTA file (can be compressed with bgzip).
    #     @param  fnfai File name of the FASTA index.
    #     @param  fngzi File name of the bgzip index.
    #     @param  flags Option flags to control index file caching and creation.
    #     @return Pointer to a faidx_t struct on success, NULL on failure.
    
    # If fnfai is NULL, ".fai" will be appended to fn to make the FAI file name.
    # If fngzi is NULL, ".gzi" will be appended to fn for the bgzip index name.
    # The bgzip index is only needed if fn is compressed.
    
    # If (flags & FAI_CREATE) is true, the index files will be built using
    # fai_build3() if they are not already present.
    # */
    faidx_t *fai_load3(const char *fn,
                       const char *fnfai,
                       const char *fngzi,
                       int flags)

    # /// Load index from "fn.fai".
    # /** @param  fn  File name of the FASTA file
    #     @return Pointer to a faidx_t struct on success, NULL on failure.
    # This function is equivalent to fai_load3(fn, NULL, NULL, FAI_CREATE|FAI_CACHE);
    # */
    faidx_t *fai_load(char *fn)

    # /// Fetch the sequence in a region
    # /** @param  fai  Pointer to the faidx_t struct
    #     @param  reg  Region in the format "chr2:20,000-30,000"
    #     @param  len  Length of the region; -2 if seq not present, -1 general error
    #     @return      Pointer to the sequence; `NULL` on failure
    # The returned sequence is allocated by `malloc()` family and should be destroyed
    # by end users by calling `free()` on it.
    # */
    char *fai_fetch(faidx_t *fai,
                    char *reg,
                    int *len)

    # /// Fetch the sequence in a region
    # /** @param  fai  Pointer to the faidx_t struct
    #     @param  c_name Region name
    #     @param  p_beg_i  Beginning position number (zero-based)
    #     @param  p_end_i  End position number (zero-based)
    #     @param  len  Length of the region; -2 if c_name not present, -1 general error
    #     @return      Pointer to the sequence; null on failure
    # The returned sequence is allocated by `malloc()` family and should be destroyed
    # by end users by calling `free()` on it.
    # */
    char *faidx_fetch_seq(faidx_t *fai,
                         char *c_name,
                         int p_beg_i,
                         int p_end_i,
                         int *len)

    # /// Query if sequence is present
    # /**   @param  fai  Pointer to the faidx_t struct
    #   @param  seq  Sequence name
    #   @return      1 if present or 0 if absent
    #   */
    int faidx_has_seq(faidx_t *fai, const char *seq)

    # /// Fetch the number of sequences
    # /** @param  fai  Pointer to the faidx_t struct
    # @return      The number of sequences
    # */
    int faidx_nseq(const faidx_t *fai)

    # /// Return name of i-th sequence
    const char *faidx_iseq(const faidx_t *fai, int i)

    # /// Return sequence length, -1 if not present
    int faidx_seq_len(faidx_t *fai, const char *seq)

# tabix support
cdef extern from "htslib/tbx.h" nogil:

    # tbx.h definitions
    int8_t TBX_MAX_SHIFT
    int32_t TBX_GENERIC
    int32_t TBX_SAM
    int32_t TBX_VCF
    int32_t TBX_UCSC

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

    int tbx_index_build(char *fn, int min_shift, tbx_conf_t *conf)
    int tbx_index_build2(const char *fn, const char *fnidx, int min_shift, const tbx_conf_t *conf)

    tbx_t * tbx_index_load(char *fn)
    tbx_t *tbx_index_load2(const char *fn, const char *fnidx)

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
    # The header keeps three dictionaries. The first keeps IDs in the
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
        int32_t n[3]                # n:the size of the dictionary block in use, (allocated size, m, is below to preserve ABI)
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
        int32_t m[3]                # m: allocated size of the dictionary block in use (see n above)

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
    uint8_t VCF_BND
    uint8_t VCF_OVERLAP


    ctypedef struct variant_t:
        int type, n     # variant type and the number of bases affected, negative for deletions

    ctypedef struct bcf_fmt_t:
        int id             # id: numeric tag id, the corresponding string is bcf_hdr_t::id[BCF_DT_ID][$id].key
        int n, size, type  # n: number of values per-sample; size: number of bytes per-sample; type: one of BCF_BT_* types
        uint8_t *p         # same as vptr and vptr_* in bcf_info_t below
        uint32_t p_len
        uint32_t p_off
        uint8_t p_free

    union bcf_info_union_t:
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
    uint8_t BCF_ERR_LIMITS
    uint8_t BCF_ERR_CHAR
    uint8_t BCF_ERR_CTG_INVALID
    uint8_t BCF_ERR_TAG_INVALID

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
    # @is_file: @samples is a file (1) or a comma-separated list (0)
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
    int bcf_write(htsFile *fp, bcf_hdr_t *h, bcf1_t *v)

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

    # bcf_hdr_merge() - copy header lines from src to dst, see also bcf_translate()
    # @param dst: the destination header to be merged into, NULL on the first pass
    # @param src: the source header
    #
    # Notes:
    #     - use as:
    #         bcf_hdr_t *dst = NULL;
    #         for (i=0; i<nsrc; i++) dst = bcf_hdr_merge(dst,src[i]);
    #
    #     - bcf_hdr_merge() replaces bcf_hdr_combine() which had a problem when
    #     combining multiple BCF headers. The current bcf_hdr_combine()
    #     does not have this problem, but became slow when used for many files.
    bcf_hdr_t *bcf_hdr_merge(bcf_hdr_t *dst, const bcf_hdr_t *src)

    # bcf_hdr_add_sample() - add a new sample.
    # @param sample:  sample name to be added
    int bcf_hdr_add_sample(bcf_hdr_t *hdr, const char *sample)

    # Read VCF header from a file and update the header
    int bcf_hdr_set(bcf_hdr_t *hdr, const char *fname)

    # Appends formatted header text to _str_.
    # If _is_bcf_ is zero, `IDX` fields are discarded.
    #  @return 0 if successful, or negative if an error occurred
    #  @since 1.4
    int bcf_hdr_format(const bcf_hdr_t *hdr, int is_bcf, kstring_t *str);
    
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
    # @param key:       tag name or NULL to remove all tags of the given type
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
    # @flt_ids:  The filter IDs to set, numeric IDs returned by bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS")
    # @n:        Number of filters. If n==0, all filters are removed
    int bcf_update_filter(const bcf_hdr_t *hdr, bcf1_t *line, int *flt_ids, int n)

    # bcf_add_filter() - adds to the FILTER column
    # @flt_id:   The filter IDs to add, numeric IDs returned by bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS")
    #
    # If flt_id is PASS, all existing filters are removed first. If other than PASS, existing PASS is removed.
    int bcf_add_filter(const bcf_hdr_t *hdr, bcf1_t *line, int flt_id)

    # bcf_remove_filter() - removes from the FILTER column
    # @flt_id:   filter ID to remove, numeric ID returned by bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS")
    # @pass:     when set to 1 and no filters are present, set to PASS
    int bcf_remove_filter(const bcf_hdr_t *hdr, bcf1_t *line, int flt_id, int set_pass)

    # Returns 1 if present, 0 if absent, or -1 if filter does not exist. "PASS" and "." can be used interchangeably.
    int bcf_has_filter(const bcf_hdr_t *hdr, bcf1_t *line, char *filter)

    # bcf_update_alleles() and bcf_update_alleles_str() - update REF and ALT column
    # @alleles:           Array of alleles
    # @nals:              Number of alleles
    # @alleles_string:    Comma-separated alleles, starting with the REF allele
    int bcf_update_alleles(const bcf_hdr_t *hdr, bcf1_t *line, const char **alleles, int nals)
    int bcf_update_alleles_str(const bcf_hdr_t *hdr, bcf1_t *line, const char *alleles_string)

    # bcf_update_id() - sets new ID string
    # bcf_add_id() - adds to the ID string checking for duplicates
    int bcf_update_id(const bcf_hdr_t *hdr, bcf1_t *line, const char *id)
    int bcf_add_id(const bcf_hdr_t *hdr, bcf1_t *line, const char *id)

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
    int bcf_get_genotypes(const bcf_hdr_t *hdr, bcf1_t *line, int32_t **dst, int *ndst)
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
    # @int_id:    return value of bcf_hdr_id2int, must be >=0
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

    hts_idx_t *bcf_index_load2(const char *fn, const char *fnidx)
    int bcf_index_build(const char *fn, int min_shift)
    int bcf_index_build2(const char *fn, const char *fnidx, int min_shift)

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
    void bcf_empty1(bcf1_t *v)
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


# VCF/BCF utility functions
cdef extern from "htslib/vcfutils.h" nogil:
    struct kbitset_t

    # bcf_trim_alleles() - remove ALT alleles unused in genotype fields
    # @header:  for access to BCF_DT_ID dictionary
    # @line:    VCF line obtain from vcf_parse1
    #
    # Returns the number of removed alleles on success or negative
    # on error:
    #     -1 .. some allele index is out of bounds
    int bcf_trim_alleles(const bcf_hdr_t *header, bcf1_t *line)

    # bcf_remove_alleles() - remove ALT alleles according to bitmask @mask
    # @header:  for access to BCF_DT_ID dictionary
    # @line:    VCF line obtained from vcf_parse1
    # @mask:    alleles to remove
    #
    # If you have more than 31 alleles, then the integer bit mask will
    # overflow, so use bcf_remove_allele_set instead
    void bcf_remove_alleles(const bcf_hdr_t *header, bcf1_t *line, int mask)

    # bcf_remove_allele_set() - remove ALT alleles according to bitset @rm_set
    # @header:  for access to BCF_DT_ID dictionary
    # @line:    VCF line obtained from vcf_parse1
    # @rm_set:  pointer to kbitset_t object with bits set for allele
    #           indexes to remove
    #
    # Number=A,R,G INFO and FORMAT fields will be updated accordingly.
    void bcf_remove_allele_set(const bcf_hdr_t *header, bcf1_t *line, kbitset_t *rm_set)

    # bcf_calc_ac() - calculate the number of REF and ALT alleles
    # @header:  for access to BCF_DT_ID dictionary
    # @line:    VCF line obtained from vcf_parse1
    # @ac:      array of length line->n_allele
    # @which:   determine if INFO/AN,AC and indv fields be used
    #
    # Returns 1 if the call succeeded, or 0 if the value could not
    # be determined.
    #
    # The value of @which determines if existing INFO/AC,AN can be
    # used (BCF_UN_INFO) and and if indv fields can be split (BCF_UN_FMT).
    int bcf_calc_ac(const bcf_hdr_t *header, bcf1_t *line, int *ac, int which)

    # bcf_gt_type() - determines type of the genotype
    # @fmt_ptr:  the GT format field as set for example by set_fmt_ptr
    # @isample:  sample index (starting from 0)
    # @ial:      index of the 1st non-reference allele (starting from 1)
    # @jal:      index of the 2nd non-reference allele (starting from 1)
    #
    # Returns the type of the genotype (one of GT_HOM_RR, GT_HET_RA,
    # GT_HOM_AA, GT_HET_AA, GT_HAPL_R, GT_HAPL_A or GT_UNKN). If $ial
    # is not NULL and the genotype has one or more non-reference
    # alleles, $ial will be set. In case of GT_HET_AA, $ial is the
    # position of the allele which appeared first in ALT. If $jal is
    # not null and the genotype is GT_HET_AA, $jal will be set and is
    # the position of the second allele in ALT.
    uint8_t GT_HOM_RR    # note: the actual value of GT_* matters, used in dosage r2 calculation
    uint8_t GT_HOM_AA
    uint8_t GT_HET_RA
    uint8_t GT_HET_AA
    uint8_t GT_HAPL_R
    uint8_t GT_HAPL_A
    uint8_t GT_UNKN
    int bcf_gt_type(bcf_fmt_t *fmt_ptr, int isample, int *ial, int *jal)

    int bcf_acgt2int(char c)
    char bcf_int2acgt(int i)

    # bcf_ij2G() - common task: allele indexes to Number=G index (diploid)
    # @i,j:  allele indexes, 0-based, i<=j
    # Returns index to the Number=G diploid array
    uint32_t bcf_ij2G(uint32_t i, uint32_t j)


cdef extern from "htslib/cram.h" nogil:

    enum cram_block_method:
        ERROR
        RAW
        GZIP
        BZIP2
        LZMA
        RANS
        RANS0
        RANS1
        GZIP_RLE

    enum cram_content_type:
        CT_ERROR
        FILE_HEADER
        COMPRESSION_HEADER
        MAPPED_SLICE
        UNMAPPED_SLICE
        EXTERNAL
        CORE

    # Opaque data types, see cram_structs for the fully fledged versions.
    ctypedef struct SAM_hdr
    ctypedef struct cram_file_def
    ctypedef struct cram_fd
    ctypedef struct cram_container
    ctypedef struct cram_block
    ctypedef struct cram_slice
    ctypedef struct cram_metrics
    ctypedef struct cram_block_slice_hdr
    ctypedef struct cram_block_compression_hdr
    ctypedef struct refs_t

    # Accessor functions

    #
    #-----------------------------------------------------------------------------
    # cram_fd
    #
    SAM_hdr *cram_fd_get_header(cram_fd *fd)
    void cram_fd_set_header(cram_fd *fd, SAM_hdr *hdr)

    int cram_fd_get_version(cram_fd *fd)
    void cram_fd_set_version(cram_fd *fd, int vers)

    int cram_major_vers(cram_fd *fd)
    int cram_minor_vers(cram_fd *fd)

    hFILE *cram_fd_get_fp(cram_fd *fd)
    void cram_fd_set_fp(cram_fd *fd, hFILE *fp)

    #
    #-----------------------------------------------------------------------------
    # cram_container
    #
    int32_t cram_container_get_length(cram_container *c)
    void cram_container_set_length(cram_container *c, int32_t length)
    int32_t cram_container_get_num_blocks(cram_container *c)
    void cram_container_set_num_blocks(cram_container *c, int32_t num_blocks)
    int32_t *cram_container_get_landmarks(cram_container *c, int32_t *num_landmarks)
    void cram_container_set_landmarks(cram_container *c, int32_t num_landmarks,
				      int32_t *landmarks)

    # Returns true if the container is empty (EOF marker) */
    int cram_container_is_empty(cram_fd *fd)


    #
    #-----------------------------------------------------------------------------
    # cram_block
    #
    int32_t cram_block_get_content_id(cram_block *b)
    int32_t cram_block_get_comp_size(cram_block *b)
    int32_t cram_block_get_uncomp_size(cram_block *b)
    int32_t cram_block_get_crc32(cram_block *b)
    void *  cram_block_get_data(cram_block *b)

    cram_content_type cram_block_get_content_type(cram_block *b)

    void cram_block_set_content_id(cram_block *b, int32_t id)
    void cram_block_set_comp_size(cram_block *b, int32_t size)
    void cram_block_set_uncomp_size(cram_block *b, int32_t size)
    void cram_block_set_crc32(cram_block *b, int32_t crc)
    void cram_block_set_data(cram_block *b, void *data)

    int cram_block_append(cram_block *b, void *data, int size)
    void cram_block_update_size(cram_block *b)

    # Offset is known as "size" internally, but it can be confusing.
    size_t cram_block_get_offset(cram_block *b)
    void cram_block_set_offset(cram_block *b, size_t offset)

    #
    # Computes the size of a cram block, including the block
    # header itself.
    #
    uint32_t cram_block_size(cram_block *b)

    #
    # Renumbers RG numbers in a cram compression header.
    #
    # CRAM stores RG as the Nth number in the header, rather than a
    # string holding the ID: tag.  This is smaller in space, but means
    # "samtools cat" to join files together that contain single but
    # different RG lines needs a way of renumbering them.
    #
    # The file descriptor is expected to be immediately after the
    # cram_container structure (ie before the cram compression header).
    # Due to the nature of the CRAM format, this needs to read and write
    # the blocks itself.  Note that there may be multiple slices within
    # the container, meaning multiple compression headers to manipulate.
    # Changing RG may change the size of the compression header and
    # therefore the length field in the container.  Hence we rewrite all
    # blocks just in case and also emit the adjusted container.
    #
    # The current implementation can only cope with renumbering a single
    # RG (and only then if it is using HUFFMAN or BETA codecs).  In
    # theory it *may* be possible to renumber multiple RGs if they use
    # HUFFMAN to the CORE block or use an external block unshared by any
    # other data series.  So we have an API that can be upgraded to
    # support this, but do not implement it for now.  An example
    # implementation of RG as an EXTERNAL block would be to find that
    # block and rewrite it, returning the number of blocks consumed.
    #
    # Returns 0 on success;
    #        -1 if unable to edit;
    #        -2 on other errors (eg I/O).
    #
    int cram_transcode_rg(cram_fd *input, cram_fd *output,
    			  cram_container *c,
			  int nrg, int *in_rg, int *out_rg)

    #
    # Copies the blocks representing the next num_slice slices from a
    # container from 'in' to 'out'.  It is expected that the file pointer
    # is just after the read of the cram_container and cram compression
    # header.
    #
    # Returns 0 on success
    #        -1 on failure
    #
    int cram_copy_slice(cram_fd *input, cram_fd *output, int32_t num_slice)

    #
    #-----------------------------------------------------------------------------
    # SAM_hdr
    #

    # Tokenises a SAM header into a hash table.
    #
    # Also extracts a few bits on specific data types, such as @RG lines.
    #
    # @return
    # Returns a SAM_hdr struct on success (free with sam_hdr_free())
    #         NULL on failure
    #
    SAM_hdr *sam_hdr_parse_(const char *hdr, int len)


    #
    #-----------------------------------------------------------------------------
    # cram_io basics
    #

    # CRAM blocks - the dynamically growable data block. We have code to
    # create, update, (un)compress and read/write.
    #
    # These are derived from the deflate_interlaced.c blocks, but with the
    # CRAM extension of content types and IDs.
    #

    # Allocates a new cram_block structure with a specified content_type and
    # id.
    #
    # @return
    # Returns block pointer on success;
    #         NULL on failure
    #
    cram_block *cram_new_block(cram_content_type content_type,
			       int content_id)

    # Reads a block from a cram file.
    #
    # @return
    # Returns cram_block pointer on success;
    #         NULL on failure
    #
    cram_block *cram_read_block(cram_fd *fd)

    # Writes a CRAM block.
    #
    # @return
    # Returns 0 on success;
    #        -1 on failure
    #
    int cram_write_block(cram_fd *fd, cram_block *b)

    # Frees a CRAM block, deallocating internal data too.
    #
    void cram_free_block(cram_block *b)

    # Uncompresses a CRAM block, if compressed.
    #
    # @return
    # Returns 0 on success;
    #        -1 on failure
    #
    int cram_uncompress_block(cram_block *b)

    # Compresses a block.
    #
    # Compresses a block using one of two different zlib strategies. If we only
    # want one choice set strat2 to be -1.
    #
    # The logic here is that sometimes Z_RLE does a better job than Z_FILTERED
    # or Z_DEFAULT_STRATEGY on quality data. If so, we'd rather use it as it is
    # significantly faster.
    #
    # @return
    # Returns 0 on success;
    #        -1 on failure
    #
    int cram_compress_block(cram_fd *fd, cram_block *b, cram_metrics *metrics,
			    int method, int level)

    # Containers
    #

    # Creates a new container, specifying the maximum number of slices
    # and records permitted.
    #
    # @return
    # Returns cram_container ptr on success;
    #         NULL on failure
    #
    cram_container *cram_new_container(int nrec, int nslice)
    void cram_free_container(cram_container *c)

    # Reads a container header.
    #
    # @return
    # Returns cram_container on success;
    #         NULL on failure or no container left (fd->err == 0).
    #
    cram_container *cram_read_container(cram_fd *fd)

    # Writes a container structure.
    #
    # @return
    # Returns 0 on success;
    #        -1 on failure
    #
    int cram_write_container(cram_fd *fd, cram_container *h)

    #
    # Stores the container structure in dat and returns *size as the
    # number of bytes written to dat[].  The input size of dat is also
    # held in *size and should be initialised to cram_container_size(c).
    #
    # Returns 0 on success;
    #        -1 on failure
    #
    int cram_store_container(cram_fd *fd, cram_container *c, char *dat, int *size)

    int cram_container_size(cram_container *c)

    # The top-level cram opening, closing and option handling
    #

    # Opens a CRAM file for read (mode "rb") or write ("wb").
    #
    # The filename may be "-" to indicate stdin or stdout.
    #
    # @return
    # Returns file handle on success;
    #         NULL on failure.
    #
    cram_fd *cram_open(const char *filename, const char *mode)

    # Opens an existing stream for reading or writing.
    #
    # @return
    # Returns file handle on success;
    #         NULL on failure.
    #
    cram_fd *cram_dopen(hFILE *fp, const char *filename, const char *mode)

    # Closes a CRAM file.
    #
    # @return
    # Returns 0 on success;
    #        -1 on failure
    #
    int cram_close(cram_fd *fd)

    #
    # Seek within a CRAM file.
    #
    # Returns 0 on success
    #        -1 on failure
    #
    int cram_seek(cram_fd *fd, off_t offset, int whence)

    #
    # Flushes a CRAM file.
    # Useful for when writing to stdout without wishing to close the stream.
    #
    # Returns 0 on success
    #        -1 on failure
    #
    int cram_flush(cram_fd *fd)

    # Checks for end of file on a cram_fd stream.
    #
    # @return
    # Returns 0 if not at end of file
    #         1 if we hit an expected EOF (end of range or EOF block)
    #         2 for other EOF (end of stream without EOF block)
    #
    int cram_eof(cram_fd *fd)

    # Sets options on the cram_fd.
    #
    # See CRAM_OPT_* definitions in hts.h.
    # Use this immediately after opening.
    #
    # @return
    # Returns 0 on success;
    #        -1 on failure
    #
    int cram_set_option(cram_fd *fd, hts_fmt_option opt, ...)

    # Sets options on the cram_fd.
    #
    # See CRAM_OPT_* definitions in hts.h.
    # Use this immediately after opening.
    #
    # @return
    # Returns 0 on success;
    #        -1 on failure
    #
    int cram_set_voption(cram_fd *fd, hts_fmt_option opt, va_list args)

    #
    # Attaches a header to a cram_fd.
    #
    # This should be used when creating a new cram_fd for writing where
    # we have an SAM_hdr already constructed (eg from a file we've read
    # in).
    #
    # @return
    # Returns 0 on success;
    #        -1 on failure
    #
    int cram_set_header(cram_fd *fd, SAM_hdr *hdr)

    # Check if this file has a proper EOF block
    #
    # @return
    # Returns 3 if the file is a version of CRAM that does not contain EOF blocks
    #         2 if the file is a stream and thus unseekable
    #         1 if the file contains an EOF block
    #         0 if the file does not contain an EOF block
    #        -1 if an error occurred whilst reading the file or we could not seek back to where we were
    #
    #
    int cram_check_EOF(cram_fd *fd)

    # As int32_decoded/encode, but from/to blocks instead of cram_fd */
    int int32_put_blk(cram_block *b, int32_t val)

    # Deallocates all storage used by a SAM_hdr struct.
    #
    # This also decrements the header reference count. If after decrementing
    # it is still non-zero then the header is assumed to be in use by another
    # caller and the free is not done.
    #
    # This is a synonym for sam_hdr_dec_ref().
    #
    void sam_hdr_free(SAM_hdr *hdr)

    # Returns the current length of the SAM_hdr in text form.
    #
    # Call sam_hdr_rebuild() first if editing has taken place.
    #
    int sam_hdr_length(SAM_hdr *hdr)

    # Returns the string form of the SAM_hdr.
    #
    # Call sam_hdr_rebuild() first if editing has taken place.
    #
    char *sam_hdr_str(SAM_hdr *hdr)

    # Appends a formatted line to an existing SAM header.
    #
    # Line is a full SAM header record, eg "@SQ\tSN:foo\tLN:100", with
    # optional new-line. If it contains more than 1 line then multiple lines
    # will be added in order.
    #
    # Len is the length of the text data, or 0 if unknown (in which case
    # it should be null terminated).
    #
    # @return
    # Returns 0 on success;
    #        -1 on failure
    #

    # Add an @PG line.
    #
    # If we wish complete control over this use sam_hdr_add() directly. This
    # function uses that, but attempts to do a lot of tedious house work for
    # you too.
    #
    # - It will generate a suitable ID if the supplied one clashes.
    # - It will generate multiple @PG records if we have multiple PG chains.
    #
    # Call it as per sam_hdr_add() with a series of key,value pairs ending
    # in NULL.
    #
    # @return
    # Returns 0 on success;
    #        -1 on failure
    #
    int sam_hdr_add_PG(SAM_hdr *sh, const char *name, ...)

    #
    # A function to help with construction of CL tags in @PG records.
    # Takes an argc, argv pair and returns a single space-separated string.
    # This string should be deallocated by the calling function.
    #
    # @return
    # Returns malloced char * on success;
    #         NULL on failure
    #
    char *stringify_argv(int argc, char *argv[])

    #
    # Returns the refs_t structure used by a cram file handle.
    #
    # This may be used in conjunction with option CRAM_OPT_SHARED_REF to
    # share reference memory between multiple file handles.
    #
    # @return
    # Returns NULL if none exists or the file handle is not a CRAM file.
    #
    refs_t *cram_get_refs(htsFile *fd)


cdef class HTSFile(object):
    cdef          htsFile *htsfile       # pointer to htsFile structure
    cdef          int64_t start_offset   # BGZF offset of first record

    cdef readonly object  filename       # filename as supplied by user
    cdef readonly object  mode           # file opening mode
    cdef readonly object  threads        # number of threads to use
    cdef readonly object  index_filename # filename of index, if supplied by user

    cdef readonly bint    is_stream      # Is htsfile a non-seekable stream
    cdef readonly bint    is_remote      # Is htsfile a remote stream
    cdef readonly bint	  duplicate_filehandle   # Duplicate filehandle when opening via fh

    cdef htsFile *_open_htsfile(self) except? NULL
