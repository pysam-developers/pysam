# cython: language_level=3
from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup
from posix.types cimport off_t


# cython does not wrap stdarg
cdef extern from "<stdarg.h>":

    ctypedef struct va_list:
        pass


# See the corresponding htslib/*.h headers for documentation of these
# HTSlib API functions.


cdef extern from "htslib/kstring.h" nogil:

    ctypedef struct kstring_t:
        size_t l, m
        char *s

    int ksprintf(kstring_t *s, const char *fmt, ...)
    int kputc(int c, kstring_t *s)
    int kputw(int c, kstring_t *s)
    int kputl(long c, kstring_t *s)


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

    hFILE *hopen(const char *filename, const char *mode, ...)
    hFILE *hdopen(int fd, const char *mode)
    int hisremote(const char *filename)
    int hclose(hFILE *fp)
    void hclose_abruptly(hFILE *fp)

    int herrno(hFILE *fp)
    void hclearerr(hFILE *fp)

    off_t hseek(hFILE *fp, off_t offset, int whence)
    off_t htell(hFILE *fp)
    int hgetc(hFILE *fp)
    ssize_t hgetdelim(char *buffer, size_t size, int delim, hFILE *fp)
    ssize_t hgetln(char *buffer, size_t size, hFILE *fp)
    char *hgets(char *buffer, int size, hFILE *fp)
    ssize_t hpeek(hFILE *fp, void *buffer, size_t nbytes)
    ssize_t hread(hFILE *fp, void *buffer, size_t nbytes)
    int hputc(int c, hFILE *fp)
    int hputs(const char *text, hFILE *fp)
    ssize_t hwrite(hFILE *fp, const void *buffer, size_t nbytes)
    int hflush(hFILE *fp)


cdef extern from "htslib/bgzf.h" nogil:

    ctypedef struct bgzf_mtaux_t
    ctypedef struct bgzidx_t
    ctypedef struct z_stream

    ctypedef struct BGZF:
        unsigned      errcode
        unsigned      is_write
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

    BGZF *bgzf_dopen(int fd, const char *mode)
    BGZF *bgzf_fdopen(int fd, const char *mode) # for backward compatibility
    BGZF *bgzf_open(const char *path, const char *mode)
    BGZF *bgzf_hopen(hFILE *fp, const char *mode)
    int bgzf_close(BGZF *fp)

    ssize_t bgzf_read(BGZF *fp, void *data, size_t length)
    ssize_t bgzf_write(BGZF *fp, const void *data, size_t length)
    ssize_t bgzf_raw_read(BGZF *fp, void *data, size_t length)
    ssize_t bgzf_raw_write(BGZF *fp, const void *data, size_t length)
    int bgzf_flush(BGZF *fp)
    int64_t bgzf_tell(const BGZF *fp)
    int64_t bgzf_seek(BGZF *fp, int64_t pos, int whence)

    int bgzf_check_EOF(BGZF *fp)
    int bgzf_is_bgzf(const char *fn)

    void bgzf_set_cache_size(BGZF *fp, int size)
    int bgzf_flush_try(BGZF *fp, ssize_t size)
    int bgzf_getc(BGZF *fp)
    int bgzf_getline(BGZF *fp, int delim, kstring_t *str)
    int bgzf_read_block(BGZF *fp)
    int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks)
    int bgzf_compress(void *dst, size_t *dlen, const void *src, size_t slen, int level)

    int bgzf_useek(BGZF *fp, long uoffset, int where)
    long bgzf_utell(BGZF *fp)

    int bgzf_index_build_init(BGZF *fp)
    int bgzf_index_load(BGZF *fp, const char *bname, const char *suffix)
    int bgzf_index_dump(BGZF *fp, const char *bname, const char *suffix)


cdef extern from "htslib/hts.h" nogil:

    int hts_verbose

    uint32_t kroundup32(uint32_t x)

    ctypedef struct cram_fd

    cdef enum htsFormatCategory:
        unknown_category
        sequence_data
        variant_data
        index_file
        region_list
        category_maximum

    cdef enum htsExactFormat:
        unknown_format
        binary_format, text_format
        sam, bam, bai, cram, crai, vcf, bcf, csi, gzi, tbi, bed
        format_maximum

    cdef enum htsCompression:
        no_compression, gzip, bgzf, custom
        compression_maximum

    ctypedef struct htsVersion:
        short major, minor

    ctypedef struct htsFormat:
        htsFormatCategory category
        htsExactFormat    format
        htsVersion        version
        htsCompression    compression
        short             compression_level
        void              *specific

    ctypedef struct hts_idx_t

    cdef union FilePointerUnion:
        BGZF    *bgzf
        cram_fd *cram
        hFILE   *hfile
        void    *voidp

    ctypedef struct htsFile:
        uint8_t is_bin, is_write, is_be, is_cram
        int64_t lineno
        kstring_t line
        char *fn
        char *fn_aux
        FilePointerUnion fp
        htsFormat format

    cdef enum hts_fmt_option:
        CRAM_OPT_DECODE_MD
        CRAM_OPT_PREFIX
        CRAM_OPT_VERBOSITY
        CRAM_OPT_SEQS_PER_SLICE
        CRAM_OPT_SLICES_PER_CONTAINER
        CRAM_OPT_RANGE
        CRAM_OPT_VERSION
        CRAM_OPT_EMBED_REF
        CRAM_OPT_IGNORE_MD5
        CRAM_OPT_REFERENCE
        CRAM_OPT_MULTI_SEQ_PER_SLICE
        CRAM_OPT_NO_REF
        CRAM_OPT_USE_BZIP2
        CRAM_OPT_SHARED_REF
        CRAM_OPT_NTHREADS
        CRAM_OPT_THREAD_POOL
        CRAM_OPT_USE_LZMA
        CRAM_OPT_USE_RANS
        CRAM_OPT_REQUIRED_FIELDS
        HTS_OPT_COMPRESSION_LEVEL
        HTS_OPT_NTHREADS

    cdef union hts_opt_val_union:
        int i
        char *s

    ctypedef struct hts_opt:
        char *arg
        hts_fmt_option opt
        hts_opt_val_union val
        void *next

    int hts_opt_add(hts_opt **opts, const char *c_arg)
    int hts_opt_apply(htsFile *fp, hts_opt *opts)
    void hts_opt_free(hts_opt *opts)

    const unsigned char seq_nt16_table[]
    const char seq_nt16_str[]
    const int seq_nt16_int[]

    const char *hts_version()
    int hts_detect_format(hFILE *fp, htsFormat *fmt)
    char *hts_format_description(const htsFormat *format)

    htsFile *hts_open(const char *fn, const char *mode)
    htsFile *hts_open_format(const char *fn, const char *mode, const htsFormat *fmt)
    htsFile *hts_hopen(hFILE *fp, const char *fn, const char *mode)
    int hts_flush(htsFile *fp)
    int hts_close(htsFile *fp)

    const htsFormat *hts_get_format(htsFile *fp)
    const char *hts_format_file_extension(const htsFormat *format)
    int hts_set_opt(htsFile *fp, hts_fmt_option opt, ...)

    int hts_getline(htsFile *fp, int delimiter, kstring_t *str)
    char **hts_readlines(const char *fn, int *_n)
    char **hts_readlist(const char *fn, int is_file, int *_n)

    int hts_set_threads(htsFile *fp, int n)
    int hts_set_fai_filename(htsFile *fp, const char *fn_aux)

    int8_t HTS_IDX_NOCOOR
    int8_t HTS_IDX_START
    int8_t HTS_IDX_REST
    int8_t HTS_IDX_NONE

    int8_t HTS_FMT_CSI
    int8_t HTS_FMT_BAI
    int8_t HTS_FMT_TBI
    int8_t HTS_FMT_CRAI

    ctypedef struct hts_pair64_t:
        uint64_t u, v



    ctypedef int hts_readrec_func(BGZF *fp, void *data, void *r, int *tid, int *beg, int *end)

    ctypedef struct hts_bins_t:
        int n, m
        int *a

    ctypedef struct hts_itr_t:
        uint32_t read_rest, finished
        int tid, beg, end, n_off, i
        int curr_tid, curr_beg, curr_end
        uint64_t curr_off
        hts_pair64_t *off
        hts_readrec_func *readrec
        hts_bins_t bins

    hts_idx_t *hts_idx_init(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls)
    void hts_idx_destroy(hts_idx_t *idx)
    int hts_idx_push(hts_idx_t *idx, int tid, int beg, int end, uint64_t offset, int is_mapped)
    void hts_idx_finish(hts_idx_t *idx, uint64_t final_offset)
    int hts_idx_save(const hts_idx_t *idx, const char *fn, int fmt)
    int hts_idx_save_as(const hts_idx_t *idx, const char *fn, const char *fnidx, int fmt)
    hts_idx_t *hts_idx_load(const char *fn, int fmt)
    hts_idx_t *hts_idx_load2(const char *fn, const char *fnidx)
    hts_idx_t *hts_idx_load3(const char *fn, const char *fnidx, int fmt, int flags)

    int HTS_IDX_SAVE_REMOTE
    int HTS_IDX_SILENT_FAIL

    ctypedef const char *(*hts_id2name_f)(void *, int)

    uint8_t *hts_idx_get_meta(hts_idx_t *idx, uint32_t *l_meta)
    void hts_idx_set_meta(hts_idx_t *idx, int l_meta, uint8_t *meta, int is_copy)

    int hts_idx_get_stat(const hts_idx_t *idx, int tid, uint64_t *mapped, uint64_t *unmapped)
    uint64_t hts_idx_get_n_no_coor(const hts_idx_t *idx)
    const char **hts_idx_seqnames(const hts_idx_t *idx, int *n, hts_id2name_f getid, void *hdr)

    int HTS_PARSE_THOUSANDS_SEP

    long long hts_parse_decimal(const char *str, char **strend, int flags)

    ctypedef int (*hts_name2id_f)(void *, const char *)

    const char *hts_parse_reg(const char *str, int *beg, int *end)

    hts_itr_t *hts_itr_query(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec)
    void hts_itr_destroy(hts_itr_t *iter)

    ctypedef hts_itr_t *hts_itr_query_func(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec)

    hts_itr_t *hts_itr_querys(const hts_idx_t *idx, const char *reg, hts_name2id_f getid, void *hdr, hts_itr_query_func *itr_query, hts_readrec_func *readrec)

    int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data)

    int FT_UNKN
    int FT_GZ
    int FT_VCF
    int FT_VCF_GZ
    int FT_BCF
    int FT_BCF_GZ
    int FT_STDIN

    int hts_file_type(const char *fname)

    ctypedef struct errmod_t

    errmod_t *errmod_init(double depcorr)
    void errmod_destroy(errmod_t *em)
    int errmod_cal(const errmod_t *em, int n, int m, uint16_t *bases, float *q)

    ctypedef struct probaln_par_t:
        float d, e
        int bw

    int probaln_glocal(const uint8_t *ref,
                       int l_ref,
                       const uint8_t *query,
                       int l_query, const uint8_t *iqual,
                       const probaln_par_t *c,
                       int *state, uint8_t *q)

    ctypedef struct hts_md5_context

    hts_md5_context *hts_md5_init()
    void hts_md5_update(hts_md5_context *ctx, const void *data, unsigned long size)
    void hts_md5_final(unsigned char *digest, hts_md5_context *ctx)
    void hts_md5_reset(hts_md5_context *ctx)
    void hts_md5_hex(char *hex, const unsigned char *digest)
    void hts_md5_destroy(hts_md5_context *ctx)

    int hts_reg2bin(int64_t beg, int64_t end, int min_shift, int n_lvls)
    int hts_bin_bot(int bin, int n_lvls)

    int ed_is_big()
    uint16_t ed_swap_2(uint16_t v)
    void *ed_swap_2p(void *x)
    uint32_t ed_swap_4(uint32_t v)
    void *ed_swap_4p(void *x)
    uint64_t ed_swap_8(uint64_t v)
    void *ed_swap_8p(void *x)


cdef extern from "htslib/sam.h" nogil:

    ctypedef struct sam_hdr_t:
        int32_t n_targets, ignore_sam_err
        uint32_t l_text
        uint32_t *target_len
        uint8_t *cigar_tab
        char **target_name
        char *text
        void *sdict

    ctypedef sam_hdr_t bam_hdr_t

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

    const char BAM_CIGAR_STR[]
    int      BAM_CIGAR_SHIFT
    uint32_t BAM_CIGAR_MASK
    uint32_t BAM_CIGAR_TYPE

    char bam_cigar_op(uint32_t c)
    uint32_t bam_cigar_oplen(uint32_t c)
    char bam_cigar_opchr(uint32_t c)
    uint32_t bam_cigar_gen(char l, uint32_t o)
    int bam_cigar_type(char o)

    int BAM_FPAIRED
    int BAM_FPROPER_PAIR
    int BAM_FUNMAP
    int BAM_FMUNMAP
    int BAM_FREVERSE
    int BAM_FMREVERSE
    int BAM_FREAD1
    int BAM_FREAD2
    int BAM_FSECONDARY
    int BAM_FQCFAIL
    int BAM_FDUP
    int BAM_FSUPPLEMENTARY

    ctypedef struct bam1_core_t:
        int32_t pos
        int32_t tid
        uint16_t bin
        uint8_t qual
        uint8_t l_extranul
        uint16_t flag
        uint8_t l_qname
        uint8_t unused1
        uint32_t n_cigar
        int32_t l_qseq
        int32_t mtid
        int32_t mpos
        int32_t isize

    ctypedef struct bam1_t:
        bam1_core_t core
        uint64_t id
        uint8_t *data
        int l_data
        uint32_t m_data

    int bam_is_rev(const bam1_t *b)
    int bam_is_mrev(const bam1_t *b)

    char *bam_get_qname(bam1_t *b)
    uint32_t *bam_get_cigar(bam1_t *b)
    char *bam_get_seq(bam1_t *b)
    uint8_t *bam_get_qual(bam1_t *b)
    uint8_t *bam_get_aux(bam1_t *b)
    int bam_get_l_aux(const bam1_t *b)
    char bam_seqi(char *s, int i)

    sam_hdr_t *bam_hdr_init()
    sam_hdr_t *bam_hdr_read(BGZF *fp)
    int bam_hdr_write(BGZF *fp, const sam_hdr_t *h)
    void bam_hdr_destroy(sam_hdr_t *h)
    sam_hdr_t *bam_hdr_dup(const sam_hdr_t *h0)

    ctypedef htsFile samFile

    sam_hdr_t *sam_hdr_parse(int l_text, const char *text)
    sam_hdr_t *sam_hdr_read(samFile *fp)
    int sam_hdr_write(samFile *fp, const sam_hdr_t *h)
    int sam_hdr_length(sam_hdr_t *h)
    const char *sam_hdr_str(sam_hdr_t *h)

    int bam_name2id(sam_hdr_t *h, const char *ref)
    char *stringify_argv(int argc, char *argv[])

    bam1_t *bam_init1()
    void bam_destroy1(bam1_t *b)
    int bam_read1(BGZF *fp, bam1_t *b)
    int bam_write1(BGZF *fp, const bam1_t *b)
    bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc)
    bam1_t *bam_dup1(const bam1_t *bsrc)

    int bam_cigar2qlen(int n_cigar, const uint32_t *cigar)
    int bam_cigar2rlen(int n_cigar, const uint32_t *cigar)
    int32_t bam_endpos(const bam1_t *b)

    int   bam_str2flag(const char *str)
    char *bam_flag2str(int flag)


    void bam_itr_destroy(hts_itr_t *iter)
    hts_itr_t *bam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end)
    hts_itr_t *bam_itr_querys(const hts_idx_t *idx, sam_hdr_t *hdr, const char *region)
    int bam_itr_next(htsFile *htsfp, hts_itr_t *itr, void *r)

    hts_idx_t *bam_index_load(const char *fn)
    int bam_index_build(const char *fn, int min_shift)
    hts_idx_t *sam_index_load(htsFile *fp, const char *fn)
    hts_idx_t *sam_index_load2(htsFile *fp, const char *fn, const char *fnidx)
    hts_idx_t *sam_index_load3(htsFile *fp, const char *fn, const char *fnidx, int flags)
    int sam_index_build(const char *fn, int min_shift)
    int sam_index_build2(const char *fn, const char *fnidx, int min_shift)

    void sam_itr_destroy(hts_itr_t *iter)
    hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end)
    hts_itr_t *sam_itr_querys(const hts_idx_t *idx, sam_hdr_t *hdr, const char *region)
    int sam_itr_next(htsFile *htsfp, hts_itr_t *itr, void *r)

    samFile *sam_open(const char *fn, const char *mode)
    samFile *sam_open_format(const char *fn, const char *mode, const htsFormat *fmt)
    int sam_close(samFile *fp)

    int sam_open_mode(char *mode, const char *fn, const char *format)
    char *sam_open_mode_opts(const char *fn, const char *mode, const char *format)

    int sam_parse1(kstring_t *s, sam_hdr_t *h, bam1_t *b)
    int sam_format1(const sam_hdr_t *h, const bam1_t *b, kstring_t *str)
    int sam_read1(samFile *fp, sam_hdr_t *h, bam1_t *b)
    int sam_write1(samFile *fp, const sam_hdr_t *h, const bam1_t *b)

    uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])
    int64_t  bam_aux2i(const uint8_t *s)
    double   bam_aux2f(const uint8_t *s)
    char     bam_aux2A(const uint8_t *s)
    char    *bam_aux2Z(const uint8_t *s)

    void bam_aux_append(bam1_t *b, const char tag[2], char type, int len, uint8_t *data)
    int bam_aux_del(bam1_t *b, uint8_t *s)

    ctypedef union bam_pileup_cd:
        void *p
        int64_t i
        double f

    ctypedef struct bam_pileup1_t:
        bam1_t *b
        int32_t qpos
        int indel, level
        uint32_t is_del, is_head, is_tail, is_refskip
        uint32_t aux
        bam_pileup_cd cd

    ctypedef int (*bam_plp_auto_f)(void *data, bam1_t *b)
    ctypedef int (*bam_test_f)()

    ctypedef struct __bam_plp_t
    ctypedef __bam_plp_t *bam_plp_t

    ctypedef struct __bam_mplp_t
    ctypedef __bam_mplp_t *bam_mplp_t

    bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data)
    void bam_plp_destroy(bam_plp_t iter)
    int bam_plp_push(bam_plp_t iter, const bam1_t *b)
    const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
    const bam_pileup1_t *bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
    void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt)
    void bam_plp_reset(bam_plp_t iter)

    ctypedef struct hts_base_mod_state

    bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data)
    void bam_mplp_init_overlaps(bam_mplp_t iter)
    void bam_mplp_destroy(bam_mplp_t iter)
    void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt)
    int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp)
    void bam_mplp_reset(bam_mplp_t iter)
    void bam_mplp_constructor(bam_mplp_t iter, int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd))
    void bam_mplp_destructor (bam_mplp_t iter, int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd))

    int sam_cap_mapq(bam1_t *b, const char *ref, int ref_len, int thres)
    int sam_prob_realn(bam1_t *b, const char *ref, int ref_len, int flag)

    ctypedef struct hts_base_mod:
        int modified_base
        int canonical_base
        int strand
        int qual

    hts_base_mod_state *hts_base_mod_state_alloc()
    void hts_base_mod_state_free(hts_base_mod_state *state)
    int bam_parse_basemod(const bam1_t *b, hts_base_mod_state *state)
    int bam_next_basemod(const bam1_t *b, hts_base_mod_state *state, hts_base_mod *mods, int n_mods, int *pos)


cdef extern from "htslib/cram.h" nogil:

    cdef enum cram_block_method:
        ERROR
        RAW
        GZIP
        BZIP2
        LZMA
        RANS
        RANS0
        RANS1
        GZIP_RLE

    cdef enum cram_content_type:
        CT_ERROR
        FILE_HEADER
        COMPRESSION_HEADER
        MAPPED_SLICE
        UNMAPPED_SLICE
        EXTERNAL
        CORE

    ctypedef struct cram_file_def
    ctypedef struct cram_fd
    ctypedef struct cram_container
    ctypedef struct cram_block
    ctypedef struct cram_slice
    ctypedef struct cram_metrics
    ctypedef struct cram_block_slice_hdr
    ctypedef struct cram_block_compression_hdr
    ctypedef struct refs_t

    sam_hdr_t *cram_fd_get_header(cram_fd *fd)
    void cram_fd_set_header(cram_fd *fd, sam_hdr_t *hdr)

    int cram_fd_get_version(cram_fd *fd)
    void cram_fd_set_version(cram_fd *fd, int vers)

    int cram_major_vers(cram_fd *fd)
    int cram_minor_vers(cram_fd *fd)

    hFILE *cram_fd_get_fp(cram_fd *fd)
    void cram_fd_set_fp(cram_fd *fd, hFILE *fp)

    int32_t cram_container_get_length(cram_container *c)
    void cram_container_set_length(cram_container *c, int32_t length)
    int32_t cram_container_get_num_blocks(cram_container *c)
    void cram_container_set_num_blocks(cram_container *c, int32_t num_blocks)
    int32_t *cram_container_get_landmarks(cram_container *c, int32_t *num_landmarks)
    void cram_container_set_landmarks(cram_container *c, int32_t num_landmarks, int32_t *landmarks)
    int cram_container_is_empty(cram_fd *fd)

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

    int cram_block_append(cram_block *b, const void *data, int size)
    void cram_block_update_size(cram_block *b)

    size_t cram_block_get_offset(cram_block *b)
    void cram_block_set_offset(cram_block *b, size_t offset)
    uint32_t cram_block_size(cram_block *b)

    int cram_transcode_rg(cram_fd *input, cram_fd *output, cram_container *c, int nrg, int *in_rg, int *out_rg)

    int cram_copy_slice(cram_fd *input, cram_fd *output, int32_t num_slice)

    cram_block *cram_new_block(cram_content_type content_type, int content_id)
    cram_block *cram_read_block(cram_fd *fd)
    int cram_write_block(cram_fd *fd, cram_block *b)
    void cram_free_block(cram_block *b)
    int cram_uncompress_block(cram_block *b)
    int cram_compress_block(cram_fd *fd, cram_block *b, cram_metrics *metrics, int method, int level)

    cram_container *cram_new_container(int nrec, int nslice)
    void cram_free_container(cram_container *c)
    cram_container *cram_read_container(cram_fd *fd)
    int cram_write_container(cram_fd *fd, cram_container *h)
    int cram_store_container(cram_fd *fd, cram_container *c, char *dat, int *size)
    int cram_container_size(cram_container *c)

    cram_fd *cram_open(const char *filename, const char *mode)
    cram_fd *cram_dopen(hFILE *fp, const char *filename, const char *mode)
    int cram_close(cram_fd *fd)
    int cram_seek(cram_fd *fd, off_t offset, int whence)
    int cram_flush(cram_fd *fd)
    int cram_eof(cram_fd *fd)
    int cram_set_option(cram_fd *fd, hts_fmt_option opt, ...)
    int cram_set_voption(cram_fd *fd, hts_fmt_option opt, va_list args)
    int cram_set_header(cram_fd *fd, sam_hdr_t *hdr)
    int cram_check_EOF(cram_fd *fd)

    int int32_put_blk(cram_block *b, int32_t val)

    ctypedef sam_hdr_t SAM_hdr
    SAM_hdr *sam_hdr_parse_(const char *hdr, int len)
    void sam_hdr_free(SAM_hdr *hdr)
    int sam_hdr_add_PG(SAM_hdr *sh, const char *name, ...)

    refs_t *cram_get_refs(htsFile *fd)


cdef extern from "htslib/faidx.h" nogil:

    ctypedef struct faidx_t:
       pass

    int fai_build3(const char *fn, const char *fnfai, const char *fngzi)
    int fai_build(const char *fn)
    void fai_destroy(faidx_t *fai)
    faidx_t *fai_load3(const char *fn, const char *fnfai, const char *fngzi, int flags)
    faidx_t *fai_load(const char *fn)
    char *fai_fetch(const faidx_t *fai, const char *reg, int *len)
    char *faidx_fetch_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len)
    int faidx_has_seq(const faidx_t *fai, const char *seq)
    int faidx_nseq(const faidx_t *fai)
    const char *faidx_iseq(const faidx_t *fai, int i)
    int faidx_seq_len(const faidx_t *fai, const char *seq)


cdef extern from "htslib/tbx.h" nogil:

    int8_t TBX_MAX_SHIFT
    int32_t TBX_GENERIC
    int32_t TBX_SAM
    int32_t TBX_VCF
    int32_t TBX_UCSC

    ctypedef struct tbx_conf_t:
        int32_t preset
        int32_t sc, bc, ec
        int32_t meta_char, line_skip

    ctypedef struct tbx_t:
        tbx_conf_t conf
        hts_idx_t *idx
        void *dict

    tbx_conf_t tbx_conf_gff
    tbx_conf_t tbx_conf_bed
    tbx_conf_t tbx_conf_psltbl
    tbx_conf_t tbx_conf_sam
    tbx_conf_t tbx_conf_vcf

    void tbx_itr_destroy(hts_itr_t *iter)
    hts_itr_t *tbx_itr_queryi(tbx_t *tbx, int tid, int beg, int end)
    hts_itr_t *tbx_itr_querys(tbx_t *tbx, const char *s)
    int tbx_itr_next(htsFile *fp, tbx_t *tbx, hts_itr_t *iter, void *r)

    int tbx_name2id(tbx_t *tbx, const char *ss)
    BGZF *hts_get_bgzfp(htsFile *fp)

    int tbx_index_build(const char *fn, int min_shift, const tbx_conf_t *conf)
    int tbx_index_build2(const char *fn, const char *fnidx, int min_shift, const tbx_conf_t *conf)

    tbx_t *tbx_index_load(const char *fn)
    tbx_t *tbx_index_load2(const char *fn, const char *fnidx)
    tbx_t *tbx_index_load3(const char *fn, const char *fnidx, int flags)

    const char **tbx_seqnames(tbx_t *tbx, int *n)

    void tbx_destroy(tbx_t *tbx)


cdef extern from "htslib/vcf.h" nogil:

    uint8_t BCF_HL_FLT
    uint8_t BCF_HL_INFO
    uint8_t BCF_HL_FMT
    uint8_t BCF_HL_CTG
    uint8_t BCF_HL_STR
    uint8_t BCF_HL_GEN

    uint8_t BCF_HT_FLAG
    uint8_t BCF_HT_INT
    uint8_t BCF_HT_REAL
    uint8_t BCF_HT_STR

    uint8_t BCF_VL_FIXED
    uint8_t BCF_VL_VAR
    uint8_t BCF_VL_A
    uint8_t BCF_VL_G
    uint8_t BCF_VL_R

    uint8_t BCF_DT_ID
    uint8_t BCF_DT_CTG
    uint8_t BCF_DT_SAMPLE

    ctypedef struct bcf_hrec_t:
        int type
        char *key
        char *value
        int nkeys
        char **keys
        char **vals

    ctypedef struct bcf_idinfo_t:
        uint32_t info[3]
        bcf_hrec_t *hrec[3]
        int id

    ctypedef struct bcf_idpair_t:
        const char *key
        const bcf_idinfo_t *val

    ctypedef struct bcf_hdr_t:
        int32_t n[3]
        bcf_idpair_t *id[3]
        void *dict[3]
        char **samples
        bcf_hrec_t **hrec
        int nhrec, dirty
        int ntransl
        int *transl[2]
        int nsamples_ori
        uint8_t *keep_samples
        kstring_t mem
        int32_t m[3]

    uint8_t bcf_type_shift[]

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
        int type, n

    ctypedef struct bcf_fmt_t:
        int id
        int n, size, type
        uint8_t *p
        uint32_t p_len
        uint32_t p_off
        uint8_t p_free

    cdef union bcf_info_union_t:
        int32_t i
        float f

    ctypedef struct bcf_info_t:
        int key
        int type
        bcf_info_union_t v1
        uint8_t *vptr
        uint32_t vptr_len
        uint32_t vptr_off
        uint8_t  vptr_free
        int len

    uint8_t BCF1_DIRTY_ID
    uint8_t BCF1_DIRTY_ALS
    uint8_t BCF1_DIRTY_FLT
    uint8_t BCF1_DIRTY_INF

    ctypedef struct bcf_dec_t:
        int m_fmt, m_info, m_id, m_als, m_allele, m_flt
        int n_flt
        int *flt
        char *id
        char *als
        char **allele
        bcf_info_t *info
        bcf_fmt_t *fmt
        variant_t *var
        int n_var, var_type
        int shared_dirty
        int indiv_dirty

    uint8_t BCF_ERR_CTG_UNDEF
    uint8_t BCF_ERR_TAG_UNDEF
    uint8_t BCF_ERR_NCOLS
    uint8_t BCF_ERR_LIMITS
    uint8_t BCF_ERR_CHAR
    uint8_t BCF_ERR_CTG_INVALID
    uint8_t BCF_ERR_TAG_INVALID

    ctypedef struct bcf1_t:
        int32_t pos
        int32_t rlen
        int32_t rid
        float qual
        uint32_t n_info, n_allele
        uint32_t n_fmt, n_sample
        kstring_t shared, indiv
        bcf_dec_t d
        int max_unpack
        int unpacked
        int unpack_size[3]
        int errcode

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

    bcf_hdr_t *bcf_hdr_init(const char *mode)
    void bcf_hdr_destroy(bcf_hdr_t *h)
    bcf1_t *bcf_init()
    void bcf_destroy(bcf1_t *v)
    void bcf_empty(bcf1_t *v)
    void bcf_clear(bcf1_t *v)

    bcf_hdr_t *bcf_hdr_read(htsFile *fp)
    int bcf_hdr_set_samples(bcf_hdr_t *hdr, const char *samples, int is_file)
    int bcf_subset_format(const bcf_hdr_t *hdr, bcf1_t *rec)
    int bcf_hdr_write(htsFile *fp, bcf_hdr_t *h)

    int vcf_parse(kstring_t *s, const bcf_hdr_t *h, bcf1_t *v)
    int vcf_format(const bcf_hdr_t *h, const bcf1_t *v, kstring_t *s)
    int bcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)

    uint8_t BCF_UN_STR
    uint8_t BCF_UN_FLT
    uint8_t BCF_UN_INFO
    uint8_t BCF_UN_SHR
    uint8_t BCF_UN_FMT
    uint8_t BCF_UN_IND
    uint8_t BCF_UN_ALL

    int bcf_unpack(bcf1_t *b, int which)

    bcf1_t *bcf_dup(bcf1_t *src)
    bcf1_t *bcf_copy(bcf1_t *dst, bcf1_t *src)
    int bcf_write(htsFile *fp, bcf_hdr_t *h, bcf1_t *v)

    bcf_hdr_t *vcf_hdr_read(htsFile *fp)
    int vcf_hdr_write(htsFile *fp, const bcf_hdr_t *h)
    int vcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
    int vcf_write(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)

    bcf_hdr_t *bcf_hdr_dup(const bcf_hdr_t *hdr)
    int bcf_hdr_combine(bcf_hdr_t *dst, const bcf_hdr_t *src)
    bcf_hdr_t *bcf_hdr_merge(bcf_hdr_t *dst, const bcf_hdr_t *src)
    int bcf_hdr_add_sample(bcf_hdr_t *hdr, const char *sample)
    int bcf_hdr_set(bcf_hdr_t *hdr, const char *fname)
    int bcf_hdr_format(const bcf_hdr_t *hdr, int is_bcf, kstring_t *str)
    char *bcf_hdr_fmt_text(const bcf_hdr_t *hdr, int is_bcf, int *len)
    int bcf_hdr_append(bcf_hdr_t *h, const char *line)
    int bcf_hdr_printf(bcf_hdr_t *h, const char *format, ...)
    const char *bcf_hdr_get_version(const bcf_hdr_t *hdr)
    void bcf_hdr_set_version(bcf_hdr_t *hdr, const char *version)
    void bcf_hdr_remove(bcf_hdr_t *h, int type, const char *key)
    bcf_hdr_t *bcf_hdr_subset(const bcf_hdr_t *h0, int n, char *const* samples, int *imap)
    const char **bcf_hdr_seqnames(const bcf_hdr_t *h, int *nseqs)
    int32_t bcf_hdr_nsamples(const bcf_hdr_t *hdr)
    int bcf_hdr_parse(bcf_hdr_t *hdr, char *htxt)
    int bcf_hdr_sync(bcf_hdr_t *h)
    bcf_hrec_t *bcf_hdr_parse_line(const bcf_hdr_t *h, const char *line, int *len)
    void bcf_hrec_format(const bcf_hrec_t *hrec, kstring_t *str)
    int bcf_hdr_add_hrec(bcf_hdr_t *hdr, bcf_hrec_t *hrec)

    bcf_hrec_t *bcf_hdr_get_hrec(const bcf_hdr_t *hdr, int type, const char *key, const char *value, const char *str_class)
    bcf_hrec_t *bcf_hrec_dup(bcf_hrec_t *hrec)
    int bcf_hrec_add_key(bcf_hrec_t *hrec, const char *str, size_t len)
    int bcf_hrec_set_val(bcf_hrec_t *hrec, int i, const char *str, size_t len, int is_quoted)
    int bcf_hrec_find_key(bcf_hrec_t *hrec, const char *key)
    int hrec_add_idx(bcf_hrec_t *hrec, int idx)
    void bcf_hrec_destroy(bcf_hrec_t *hrec)

    int bcf_subset(const bcf_hdr_t *h, bcf1_t *v, int n, int *imap)
    int bcf_translate(const bcf_hdr_t *dst_hdr, bcf_hdr_t *src_hdr, bcf1_t *src_line)
    int bcf_get_variant_types(bcf1_t *rec)
    int bcf_get_variant_type(bcf1_t *rec, int ith_allele)
    int bcf_is_snp(bcf1_t *v)

    int bcf_update_filter(const bcf_hdr_t *hdr, bcf1_t *line, int *flt_ids, int n)
    int bcf_add_filter(const bcf_hdr_t *hdr, bcf1_t *line, int flt_id)
    int bcf_remove_filter(const bcf_hdr_t *hdr, bcf1_t *line, int flt_id, int pass_)
    int bcf_has_filter(const bcf_hdr_t *hdr, bcf1_t *line, char *filter)

    int bcf_update_alleles(const bcf_hdr_t *hdr, bcf1_t *line, const char **alleles, int nals)
    int bcf_update_alleles_str(const bcf_hdr_t *hdr, bcf1_t *line, const char *alleles_string)

    int bcf_update_id(const bcf_hdr_t *hdr, bcf1_t *line, const char *id)
    int bcf_add_id(const bcf_hdr_t *hdr, bcf1_t *line, const char *id)

    int bcf_update_info_int32(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const int32_t *values, int n)
    int bcf_update_info_float(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const float *values, int n)
    int bcf_update_info_flag(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const char *string, int n)
    int bcf_update_info_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const char *string, int n)
    int bcf_update_info(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type)

    int bcf_update_format_int32(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const int32_t *values, int n)
    int bcf_update_format_float(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const float *values, int n)
    int bcf_update_format_char(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const char *values, int n)
    int bcf_update_genotypes(const bcf_hdr_t *hdr, bcf1_t *line, const int32_t *gts, int n)
    int bcf_update_format_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const char **values, int n)
    int bcf_update_format(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type)

    uint32_t bcf_gt_phased(uint32_t idx)
    uint32_t bcf_gt_unphased(uint32_t idx)
    uint32_t bcf_gt_missing
    uint32_t bcf_gt_is_missing(uint32_t val)
    uint32_t bcf_gt_is_phased(uint32_t idx)
    uint32_t bcf_gt_allele(uint32_t val)

    uint32_t bcf_alleles2gt(uint32_t a, uint32_t b)
    void bcf_gt2alleles(int igt, int *a, int *b)

    bcf_fmt_t *bcf_get_fmt(const bcf_hdr_t *hdr, bcf1_t *line, const char *key)
    bcf_info_t *bcf_get_info(const bcf_hdr_t *hdr, bcf1_t *line, const char *key)
    bcf_fmt_t *bcf_get_fmt_id(bcf1_t *line, const int id)
    bcf_info_t *bcf_get_info_id(bcf1_t *line, const int id)

    int bcf_get_info_int32(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, int32_t **dst, int *ndst)
    int bcf_get_info_float(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, float **dst, int *ndst)
    int bcf_get_info_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, char **dst, int *ndst)
    int bcf_get_info_flag(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, int **dst, int *ndst)
    int bcf_get_info_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type)

    int bcf_get_format_int32(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, int32_t **dst, int *ndst)
    int bcf_get_format_float(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, float **dst, int *ndst)
    int bcf_get_format_char(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, char **dst, int *ndst)
    int bcf_get_genotypes(const bcf_hdr_t *hdr, bcf1_t *line, int32_t **dst, int *ndst)
    int bcf_get_format_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, char ***dst, int *ndst)
    int bcf_get_format_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type)

    int bcf_hdr_id2int(const bcf_hdr_t *hdr, int type, const char *id)
    const char *bcf_hdr_int2id(const bcf_hdr_t *hdr, int type, int int_id)
    int bcf_hdr_name2id(const bcf_hdr_t *hdr, const char *id)
    const char *bcf_hdr_id2name(const bcf_hdr_t *hdr, int rid)
    const char *bcf_seqname(const bcf_hdr_t *hdr, const bcf1_t *rec)

    int bcf_hdr_id2length(const bcf_hdr_t *hdr, int type, int int_id)
    int bcf_hdr_id2number(const bcf_hdr_t *hdr, int type, int int_id)
    int bcf_hdr_id2type(const bcf_hdr_t *hdr, int type, int int_id)
    int bcf_hdr_id2coltype(const bcf_hdr_t *hdr, int type, int int_id)
    int bcf_hdr_idinfo_exists(const bcf_hdr_t *hdr, int type, int int_id)
    bcf_hrec_t *bcf_hdr_id2hrec(const bcf_hdr_t *hdr, int dict_type, int col_type, int int_id)

    void bcf_fmt_array(kstring_t *s, int n, int type, void *data)
    uint8_t *bcf_fmt_sized_array(kstring_t *s, uint8_t *ptr)

    void bcf_enc_vchar(kstring_t *s, int l, const char *a)
    void bcf_enc_vint(kstring_t *s, int n, int32_t *a, int wsize)
    void bcf_enc_vfloat(kstring_t *s, int n, float *a)

    void bcf_itr_destroy(hts_itr_t *iter)
    hts_itr_t *bcf_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end)
    hts_itr_t *bcf_itr_querys(const hts_idx_t *idx, const bcf_hdr_t *hdr, const char *s)
    int bcf_itr_next(htsFile *htsfp, hts_itr_t *itr, void *r)

    hts_idx_t *bcf_index_load(const char *fn)
    const char **bcf_index_seqnames(const hts_idx_t *idx, const bcf_hdr_t *hdr, int *nptr)
    hts_idx_t *bcf_index_load2(const char *fn, const char *fnidx)
    hts_idx_t *bcf_index_load3(const char *fn, const char *fnidx, int flags)
    int bcf_index_build(const char *fn, int min_shift)
    int bcf_index_build2(const char *fn, const char *fnidx, int min_shift)

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


cdef extern from "htslib/vcfutils.h" nogil:

    cdef struct kbitset_t

    int bcf_trim_alleles(const bcf_hdr_t *header, bcf1_t *line)
    void bcf_remove_alleles(const bcf_hdr_t *header, bcf1_t *line, int mask)
    void bcf_remove_allele_set(const bcf_hdr_t *header, bcf1_t *line, kbitset_t *rm_set)
    int bcf_calc_ac(const bcf_hdr_t *header, bcf1_t *line, int *ac, int which)

    uint8_t GT_HOM_RR
    uint8_t GT_HOM_AA
    uint8_t GT_HET_RA
    uint8_t GT_HET_AA
    uint8_t GT_HAPL_R
    uint8_t GT_HAPL_A
    uint8_t GT_UNKN
    int bcf_gt_type(bcf_fmt_t *fmt_ptr, int isample, int *ial, int *jal)

    int bcf_acgt2int(char c)
    char bcf_int2acgt(int i)

    uint32_t bcf_ij2G(uint32_t i, uint32_t j)


cdef class HTSFile(object):
    cdef          htsFile *htsfile       # pointer to htsFile structure
    cdef          int64_t start_offset   # BGZF offset of first record

    cdef readonly object  filename       # filename as supplied by user
    cdef readonly object  mode           # file opening mode
    cdef readonly object  threads        # number of threads to use
    cdef readonly object  index_filename # filename of index, if supplied by user

    cdef readonly bint    is_stream      # Is htsfile a non-seekable stream
    cdef readonly bint    is_remote      # Is htsfile a remote stream
    cdef readonly bint    duplicate_filehandle   # Duplicate filehandle when opening via fh

    cdef htsFile *_open_htsfile(self) except? NULL
