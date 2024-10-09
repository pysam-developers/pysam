/*  vcfsort.c -- sort subcommand

   Copyright (C) 2017-2024 Genome Research Ltd.

   Author: Petr Danecek <pd3@sanger.ac.uk>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <math.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/hts_os.h>
#include <htslib/hts_defs.h>
#include <htslib/bgzf.h>
#include "kheap.h"
#include "bcftools.h"

#define MAX_TMP_FILES_PER_LAYER 32
#define MERGE_LAYERS 12
#define MAX_TMP_FILES (MAX_TMP_FILES_PER_LAYER * MERGE_LAYERS)

typedef struct
{
    char *fname;
    htsFile *fh;
    BGZF *bgz;
    size_t idx;
    bcf1_t *rec;
    int is_merged;
}
blk_t;

typedef struct
{
    size_t len;
    hts_pos_t pos;
    int rid;
    float qual;
    uint8_t data[];
}
packed_bcf_t;

typedef struct _args_t
{
    bcf_hdr_t *hdr;
    char **argv, *fname, *output_fname, *tmp_dir;
    int argc, output_type, clevel;
    size_t max_mem, mem;
    packed_bcf_t **buf;
    uint8_t *mem_block;

    size_t nbuf, mbuf, nblk, tmp_count;
    blk_t blk[MAX_TMP_FILES];
    uint32_t tmp_layers[MERGE_LAYERS];
    int write_index;
}
args_t;

void clean_files(args_t *args)
{
    int i;
    fprintf(stderr,"Cleaning\n");
    for (i=0; i<MAX_TMP_FILES; i++)
    {
        blk_t *blk = &args->blk[i];
        if ( blk->fname )
        {
            unlink(blk->fname);
            free(blk->fname);
        }
        if ( blk->rec )
            bcf_destroy(blk->rec);
    }
    rmdir(args->tmp_dir);
}
void HTS_FORMAT(HTS_PRINTF_FMT, 2, 3) HTS_NORETURN
clean_files_and_throw(args_t *args, const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    clean_files(args);
    exit(-1);
}

int cmp_bcf_pos(const void *aptr, const void *bptr)
{
    bcf1_t *a = *((bcf1_t**)aptr);
    bcf1_t *b = *((bcf1_t**)bptr);
    if ( a->rid < b->rid ) return -1;
    if ( a->rid > b->rid ) return 1;
    if ( a->pos < b->pos ) return -1;
    if ( a->pos > b->pos ) return 1;
    return 0;
}
int cmp_bcf_pos_ref_alt(const void *aptr, const void *bptr)
{
    bcf1_t *a = *((bcf1_t**)aptr);
    bcf1_t *b = *((bcf1_t**)bptr);
    if ( a->rid < b->rid ) return -1;
    if ( a->rid > b->rid ) return 1;
    if ( a->pos < b->pos ) return -1;
    if ( a->pos > b->pos ) return 1;

    // Sort the same chr:pos records lexicographically by ref,alt.
    // This will be called rarely so should not slow the sorting down
    // noticeably.

    int i;
    for (i=0; i<a->n_allele; i++)
    {
        if ( i >= b->n_allele ) return 1;
        int ret = strcasecmp(a->d.allele[i],b->d.allele[i]);
        if ( ret ) return ret;
    }
    if ( a->n_allele < b->n_allele ) return -1;
    return 0;
}

static int cmp_packed_bcf_pos_ref_alt(const void *aptr, const void *bptr)
{
    packed_bcf_t *a = *(packed_bcf_t **) aptr;
    packed_bcf_t *b = *(packed_bcf_t **) bptr;

    if ( a->rid < b->rid ) return -1;
    if ( a->rid > b->rid ) return 1;
    if ( a->pos < b->pos ) return -1;
    if ( a->pos > b->pos ) return 1;
    
    // Sort lexicographically by ref,alt.  These are stored tab-separated
    // as the first item in packed_bcf_t::data
    return strcmp((char *) a->data, (char *) b->data);
}

static int cmp_packed_bcf_pos_ref_alt_stable(const void *aptr, const void *bptr)
{
    // cmp_bcf_pos_ref_alt() with tie-breaker to make qsort stable
    int res = cmp_packed_bcf_pos_ref_alt(aptr, bptr);
    if (res != 0) return res;

    // Got a tie - use the position in the original input to break it.
    // As everything is read into a big memory buffer, for most records
    // we can just compare the pointers directly.  The exception is
    // any record that didn't quite fit in the memory buffer, causing it to be
    // flushed.  Those are flagged by setting packed_bcf_t::len = SIZE_MAX, and
    // as they were the last record in the segment, they should always sort
    // after unflagged records.

    packed_bcf_t *a = *(packed_bcf_t **) aptr;
    packed_bcf_t *b = *(packed_bcf_t **) bptr;

    if (a->len == SIZE_MAX) return 1;
    if (b->len == SIZE_MAX) return -1;

    return a < b ? -1 : 1;
}

static uint8_t *pack_unsigned(uint8_t *data, uint64_t val)
{
    do {
        *data++ = (val & 0x7f) | ((val > 0x7f) ? 0x80 : 0);
        val >>= 7;
    } while (val > 0);
    return data;    
}

static uint8_t *pack_hts_pos(uint8_t *data, hts_pos_t val)
{
    uint64_t sign = val < 0;
    uint64_t v = val < 0 ? -(val + 1) : val;
    v = v << 1 | sign;
    return pack_unsigned(data, v);
}

static uint8_t *pack_bcf_data(packed_bcf_t *dest, const bcf1_t *src,
                              int outside_buffer)
{
    uint32_t i;
    uint8_t *data = dest->data;
    uint8_t *start = dest->data;
    dest->pos = src->pos;
    dest->rid = src->rid;
    dest->qual = src->qual;

    // Copy in alleles, for the comparison function
    for (i = 0; i < src->n_allele; i++)
    {
        size_t l = strlen(src->d.allele[i]);
        if (i > 0) *data++ = '\t';
        memcpy(data, src->d.allele[i], l);
        data += l;
    }
    *data++ = '\0';

    if (outside_buffer)
    {
        dest->len = SIZE_MAX;
        memcpy(data, &src, sizeof(src));
        data += sizeof(src);
        return data;
    }

    data = pack_hts_pos(data, src->rlen);
    data = pack_unsigned(data, src->n_info);
    data = pack_unsigned(data, src->n_allele);
    data = pack_unsigned(data, src->n_fmt);
    data = pack_unsigned(data, src->n_sample);
    data = pack_unsigned(data, src->shared.l);
    data = pack_unsigned(data, src->indiv.l);
    if (src->shared.l)
        memcpy(data, src->shared.s, src->shared.l);
    data += src->shared.l;
    if (src->indiv.l)
        memcpy(data, src->indiv.s, src->indiv.l);
    data += src->indiv.l;
    dest->len = data - start;
    return data;
}

static int write_packed_bcf(BGZF *fp, packed_bcf_t *src)
{
    // Write pos, rid, qual
    size_t len = src->data - (uint8_t *) &src->pos;
    if (bgzf_write_small(fp, &src->pos, len) < len)
        return -1;

    // Skip the copy of the alleles
    size_t skip = strlen((char *) src->data) + 1;
 
    // Write everything else
    if (src->len < SIZE_MAX)
    {
        // In main memory block
        len = src->len - skip;
        if (bgzf_write_small(fp, src->data + skip, len) < len)
            return -1;
    }
    else
    {
        // Record didn't fit in the main block.  To minimize the
        // overflow, its packed_bcf_t data will be imcomplete.  A pointer to
        // its bcf1_t struct will have been placed after the allele data
        // so we can finish the packing job and write it in the same format
        // the rest of the data
        bcf1_t *rec;
        uint8_t tmp[100], *data = tmp;
        memcpy(&rec, src->data + skip, sizeof(rec));

        data = pack_hts_pos(data, rec->rlen);
        data = pack_unsigned(data, rec->n_info);
        data = pack_unsigned(data, rec->n_allele);
        data = pack_unsigned(data, rec->n_fmt);
        data = pack_unsigned(data, rec->n_sample);
        data = pack_unsigned(data, rec->shared.l);
        data = pack_unsigned(data, rec->indiv.l);
        if (bgzf_write_small(fp, tmp, data - tmp) < data - tmp)
            return -1;
        if (rec->shared.l > 0 &&
            bgzf_write_small(fp, rec->shared.s, rec->shared.l) < rec->shared.l)
            return -1;
        if (rec->indiv.l > 0 &&
            bgzf_write_small(fp, rec->indiv.s, rec->indiv.l) < rec->indiv.l)
            return -1;
    }

    return 0;
}

static uint64_t unpack_unsigned(BGZF *fp, int *err)
{
    uint8_t data;
    uint64_t val = 0;
    uint32_t i = 0;

    if (bgzf_read_small(fp, &data, sizeof(data)) <= 0)
        goto short_read;

    while (data & 0x80)
    {
        val |= (uint64_t)(data & 0x7f) << i;
        i += 7;
        if (bgzf_read_small(fp, &data, sizeof(data)) <= 0)
            goto short_read;
    }
    val |= (uint64_t)data << i;
    return val;

 short_read:
    *err = 1;
    return 0;
}

static hts_pos_t unpack_hts_pos(BGZF *fp, int *err)
{
    uint64_t v = unpack_unsigned(fp, err);

    if ((v & 1) == 0)
        return (hts_pos_t)(v >> 1);
    else
        return -(hts_pos_t)(v >> 1) - 1;
}

static int read_packed_bcf(BGZF *fp, bcf1_t *dest)
{
    int err = 0;
    packed_bcf_t tmp;
    size_t len = tmp.data - (uint8_t *) &tmp.pos;
    
    bcf_clear(dest);
    ssize_t got = bgzf_read_small(fp, &tmp.pos, len);
    if (got == 0)
        return -1;  // EOF
    if (got < len)
        return -2;  // Error or short read
    dest->pos = tmp.pos;
    dest->rid = tmp.rid;
    dest->qual = tmp.qual;
    dest->rlen = unpack_hts_pos(fp, &err);
    dest->n_info = unpack_unsigned(fp, &err);
    dest->n_allele = unpack_unsigned(fp, &err);
    dest->n_fmt = unpack_unsigned(fp, &err);
    dest->n_sample = unpack_unsigned(fp, &err);
    len = unpack_unsigned(fp, &err);
    if (ks_resize(&dest->shared, len) != 0)
        return -2;
    dest->shared.l = len;
    len = unpack_unsigned(fp, &err);
    if (ks_resize(&dest->indiv, len) != 0)
        return -2;
    dest->indiv.l = len;
    err |= bgzf_read_small(fp, dest->shared.s, dest->shared.l) < dest->shared.l;
    err |= bgzf_read_small(fp, dest->indiv.s, dest->indiv.l) < dest->indiv.l;
    return err == 0 ? 0 : -2;
}

void open_tmp_file(args_t *args, blk_t *blk, int is_merged)
{
    kstring_t str = {0,0,0};
    int tries = 1000;

    blk->fh = NULL;
    blk->bgz = NULL;

    do {
        if (ksprintf(ks_clear(&str), "%s/%05zd%s",
                     args->tmp_dir, args->tmp_count++,
                     is_merged ? ".bcf" : "") < 0) {
            clean_files_and_throw(args, "%s", strerror(errno));
        }

        if (is_merged)
            blk->fh = hts_open(str.s, "wbx1");
        else
            blk->bgz = bgzf_open(str.s, "wx1");
        if ( blk->fh == NULL && blk->bgz == NULL && (errno != EEXIST || --tries <= 0)) {
            clean_files_and_throw(args, "Cannot write %s: %s\n",
                                  str.s, strerror(errno));
        }
    } while (blk->fh == NULL && blk->bgz == NULL);

    blk->fname = ks_release(&str);
    blk->idx = args->tmp_count - 1;
}

void do_partial_merge(args_t *args);

void buf_flush(args_t *args, bcf1_t *last_rec)
{
    if ( !args->nbuf ) return;

    qsort(args->buf, args->nbuf, sizeof(*args->buf), cmp_packed_bcf_pos_ref_alt_stable);

    if (args->tmp_layers[0] >= MAX_TMP_FILES_PER_LAYER)
        do_partial_merge(args);

    assert(args->nblk < MAX_TMP_FILES);
    blk_t *blk = &args->blk[args->nblk];
    blk->is_merged = 0;
    args->nblk++;
    args->tmp_layers[0]++;

    assert(blk->fname == NULL && blk->fh == NULL && blk->bgz == NULL);

    open_tmp_file(args, blk, 0);
    int i;
    for (i=0; i<args->nbuf; i++)
    {
        if ( write_packed_bcf(blk->bgz, args->buf[i])!=0 ) clean_files_and_throw(args, "[%s] Error: cannot write to %s\n", __func__,blk->fname);
    }

    if ( bgzf_close(blk->bgz)!=0 ) clean_files_and_throw(args, "[%s] Error: close failed .. %s\n", __func__,blk->fname);
    blk->bgz = NULL;

    args->nbuf = 0;
    args->mem  = 0;
}


static inline uint8_t *_align_up(uint8_t *ptr)
{
    return (uint8_t*)(((size_t)ptr + 8 - 1) & ~((size_t)(8 - 1)));
}

#define varint_size(X) ((sizeof(X) * 8 + 7) / 7) // worst case

void buf_push(args_t *args, bcf1_t *rec)
{
    size_t delta = sizeof(rec->pos)
        + sizeof(rec->rid)
        + sizeof(rec->qual)
        + varint_size(rec->rlen)
        + varint_size(2) // n_info
        + varint_size(2) // n_allele
        + varint_size(1) // n_fmt
        + varint_size(3) // n_sample
        + varint_size(rec->shared.l)
        + varint_size(rec->indiv.l)
        + rec->shared.l + rec->indiv.l
        + rec->unpack_size[1]   // Alleles
        + 8;                    // the number of _align_up() calls

    if ( delta > args->max_mem - args->mem )
    {
        packed_bcf_t *tmp = malloc(sizeof(*tmp) + rec->unpack_size[1] * sizeof(bcf1_t *));
        if (!tmp)
            clean_files_and_throw(args, "[%s] Out of memory\n", __func__);
        pack_bcf_data(tmp, rec, 1);

        args->nbuf++;
        hts_expand(bcf1_t*, args->nbuf, args->mbuf, args->buf);
        args->buf[args->nbuf-1] = tmp;

        buf_flush(args, rec);

        free(tmp);
        bcf_destroy(rec);
        return;
    }

    // make sure nothing has changed in htslib
    assert( rec->unpacked==BCF_UN_STR && !rec->d.flt && !rec->d.info && !rec->d.fmt && !rec->d.var );

    uint8_t *ptr_beg = args->mem_block + args->mem;
    uint8_t *ptr = _align_up(ptr_beg);
    packed_bcf_t *packed_rec = (packed_bcf_t *) ptr;

    ptr = pack_bcf_data(packed_rec, rec, 0);

    args->nbuf++;
    hts_expand(bcf1_t*, args->nbuf, args->mbuf, args->buf);
    args->buf[args->nbuf-1] = packed_rec;

    delta = ptr - ptr_beg;
    args->mem += delta;

    assert( args->mem <= args->max_mem );

    bcf_destroy(rec);
}

void sort_blocks(args_t *args)
{
    htsFile *in = hts_open(args->fname, "r");
    if ( !in ) clean_files_and_throw(args, "Could not read %s\n", args->fname);
    args->hdr = bcf_hdr_read(in);
    if ( !args->hdr) clean_files_and_throw(args, "Could not read VCF/BCF headers from %s\n", args->fname);

    while ( 1 )
    {
        bcf1_t *rec = bcf_init();
        int ret = bcf_read1(in, args->hdr, rec);
        if ( ret < -1 ) clean_files_and_throw(args,"Error encountered while parsing the input\n");
        if ( ret == -1 )
        {
            bcf_destroy(rec);
            break;
        }
        if ( rec->errcode ) clean_files_and_throw(args,"Error encountered while parsing the input at %s:%"PRIhts_pos"\n",bcf_seqname(args->hdr,rec),rec->pos+1);
        bcf_unpack(rec, BCF_UN_STR);
        buf_push(args, rec);
    }
    buf_flush(args, NULL);
    free(args->buf);

    if ( hts_close(in)!=0 ) clean_files_and_throw(args,"Close failed: %s\n", args->fname);
}

static inline int blk_is_smaller(blk_t **aptr, blk_t **bptr)
{
    blk_t *a = *aptr;
    blk_t *b = *bptr;
    int ret = cmp_bcf_pos_ref_alt(&a->rec, &b->rec);
    if ( ret < 0 ) return 1;
    if (ret == 0 && a->idx < b->idx) return 1;
    return 0;
}
KHEAP_INIT(blk, blk_t*, blk_is_smaller)

void blk_read(args_t *args, khp_blk_t *bhp, bcf_hdr_t *hdr, blk_t *blk)
{
    int ret;
    if (blk->is_merged)
    {
        if ( !blk->fh ) return;
        ret = bcf_read(blk->fh, hdr, blk->rec);
    }
    else
    {
        if ( !blk->bgz ) return;
        ret = read_packed_bcf(blk->bgz, blk->rec);
    }
    if ( ret < -1 ) clean_files_and_throw(args, "Error reading %s\n", blk->fname);
    if ( ret == -1 )
    {
        if (blk->is_merged)
        {
            if ( hts_close(blk->fh)!=0 ) clean_files_and_throw(args, "Close failed: %s\n", blk->fname);
            blk->fh = NULL;
            return;
        }
        else
        {
            if ( bgzf_close(blk->bgz) != 0)
                clean_files_and_throw(args, "Close failed: %s\n", blk->fname);
            blk->bgz = NULL;
            return;
        }
    }
    bcf_unpack(blk->rec, BCF_UN_STR);
    khp_insert(blk, bhp, &blk);
}

void merge_blocks(args_t *args, htsFile *out, const char *output_fname,
                  int idx_fmt, size_t from)
{
    khp_blk_t *bhp = khp_init(blk);
    char *index_fn = NULL;
    size_t i;

    for (i=from; i<args->nblk; i++)
    {
        blk_t *blk = &args->blk[i];
        if (blk->is_merged)
        {
            blk->fh = hts_open(blk->fname, "r");
            if ( !blk->fh ) clean_files_and_throw(args, "Could not read %s: %s\n", blk->fname, strerror(errno));
            bcf_hdr_t *hdr = bcf_hdr_read(blk->fh);
            bcf_hdr_destroy(hdr);
        }
        else
        {
            blk->bgz = bgzf_open(blk->fname, "r");
            if (!blk->bgz)
                clean_files_and_throw(args, "Could not read %s: %s\n", blk->fname, strerror(errno));
        }
        blk_read(args, bhp, args->hdr, blk);
    }

    if ( bcf_hdr_write(out, args->hdr)!=0 ) clean_files_and_throw(args, "[%s] Error: cannot write to %s\n", __func__, output_fname);

    if (idx_fmt) {
        if ( init_index2(out,args->hdr,output_fname,&index_fn,idx_fmt)<0 )
            error("Error: failed to initialise index for %s\n",output_fname);
    }

    while ( bhp->ndat )
    {
        blk_t *blk = bhp->dat[0];
        if ( bcf_write(out, args->hdr, blk->rec)!=0 ) clean_files_and_throw(args, "[%s] Error: cannot write to %s\n", __func__,args->output_fname);
        khp_delete(blk, bhp);
        blk_read(args, bhp, args->hdr, blk);
    }
    if ( idx_fmt )
    {
        if ( bcf_idx_save(out)<0 )
        {
            if ( hts_close(out)!=0 ) error("Error: close failed .. %s\n", output_fname);
            error("Error: cannot write to index %s\n", index_fn);
        }
        free(index_fn);
    }

    for (i = from; i < args->nblk; i++)
    {
        blk_t *blk = &args->blk[i];
        if (unlink(blk->fname) != 0)
            clean_files_and_throw(args, "Couldn't remove temporary file %s\n", blk->fname);
        free(blk->fname);
        blk->fname = NULL;
    }

    khp_destroy(blk, bhp);
}

void do_partial_merge(args_t *args)
{
    uint32_t to_layer = 0;
    size_t to_merge = 0;

    // Temp. files are arranged in layers of at most MAX_TMP_FILES_PER_LAYER.
    // When a layer is full, it is merged into the next layer up.  Each
    // layer will therefore contain files with exponentially more records
    // then the previous one, but will be merged exponentially less frequently.
    // The result is that the overall complexity will remain O(n*log(n))
    // even if we need to do lots of partial merges.

    while (to_layer < MERGE_LAYERS
           && args->tmp_layers[to_layer] >= MAX_TMP_FILES_PER_LAYER)
    {
        to_merge += args->tmp_layers[to_layer];
        args->tmp_layers[to_layer] = 0;
        to_layer++;
    }

    assert(to_merge > 0 && to_merge <= args->nblk);

    if (to_layer == MERGE_LAYERS) {
        // Edge case - if we've got here, we've completely used the
        // temp file allocation, so merge absolutely everything and
        // leave one file at the highest level.  Strictly this breaks
        // the O(n*log(n)) complexity, but unless MERGE_LAYERS and
        // MAX_TMP_FILES_PER_LAYER are too small it would take so long
        // to get here it should never actually happen...
        assert(to_merge == MAX_TMP_FILES_PER_LAYER * MERGE_LAYERS);
        to_layer = MERGE_LAYERS - 1;
    }

    blk_t tmp = { NULL };
    open_tmp_file(args, &tmp, 1);
    merge_blocks(args, tmp.fh, tmp.fname, 0, args->nblk - to_merge);
    if (hts_close(tmp.fh) != 0)
        clean_files_and_throw(args, "Close failed: %s\n", tmp.fname);

    args->nblk -= to_merge;
    assert(args->blk[args->nblk].fh == NULL);
    assert(args->blk[args->nblk].fname == NULL);
    args->blk[args->nblk].is_merged = 1;
    args->blk[args->nblk].idx = tmp.idx;
    args->blk[args->nblk++].fname = tmp.fname;
    args->tmp_layers[to_layer]++;
}

void merge_to_output(args_t *args)
{
    char wmode[8] = { 0 };
    set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
    const char *output_fname = args->output_fname ? args->output_fname : "-";

    htsFile *out = hts_open(output_fname, wmode);
    if (!out) clean_files_and_throw(args, "[%s] Error: cannot open %s\n", __func__, output_fname);

    fprintf(stderr,"Merging %zd temporary files\n", args->nblk);
    merge_blocks(args, out, output_fname, args->write_index, 0);
    fprintf(stderr,"Done\n");

    if ( hts_close(out)!=0 )
        clean_files_and_throw(args, "Close failed: %s\n", output_fname);

    clean_files(args);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Sort VCF/BCF file.\n");
    fprintf(stderr, "Usage:   bcftools sort [OPTIONS] <FILE.vcf>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -m, --max-mem FLOAT[kMG]       maximum memory to use [768M]\n");    // using metric units, 1M=1e6
    fprintf(stderr, "    -o, --output FILE              output file name [stdout]\n");
    fprintf(stderr, "    -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n");

#ifdef _WIN32
    fprintf(stderr, "    -T, --temp-dir DIR             temporary files [/bcftools.XXXXXX]\n");
#else
    fprintf(stderr, "    -T, --temp-dir DIR             temporary files [/tmp/bcftools.XXXXXX]\n");
#endif
    fprintf(stderr, "    -W, --write-index[=FMT]        Automatically index the output files [off]\n");
    fprintf(stderr, "\n");
    exit(1);
}

size_t parse_mem_string(const char *str)
{
    char *tmp;
    double mem = strtod(str, &tmp);
    if ( tmp==str ) error("Could not parse the memory string: \"%s\"\n", str);
    if ( !strcasecmp("k",tmp) ) mem *= 1000;
    else if ( !strcasecmp("m",tmp) ) mem *= 1000*1000;
    else if ( !strcasecmp("g",tmp) ) mem *= 1000*1000*1000;
    return mem;
}

static void init(args_t *args)
{
    size_t i;
    args->max_mem *= 0.9;
    args->mem_block = malloc(args->max_mem);
    if ( !args->mem_block ) error("Error: could not allocate %zu bytes of memory, try reducing --max-mem\n",args->max_mem);
    args->mem = 0;

    for (i = 0; i < MAX_TMP_FILES; i++)
    {
        args->blk[i].fname = NULL;
        args->blk[i].rec = bcf_init();
        if (!args->blk[i].rec)
            clean_files_and_throw(args,"Couldn't allocate bcf record\n");
    }


    args->tmp_dir = init_tmp_prefix(args->tmp_dir);

#ifdef _WIN32
        int ret = mkdir(mktemp(args->tmp_dir), 0700);
        if ( ret ) error("mkdir(%s) failed: %s\n", args->tmp_dir,strerror(errno));
#else
        char *tmp = mkdtemp(args->tmp_dir);
        if ( !tmp ) error("mkdtemp(%s) failed: %s\n",  args->tmp_dir,strerror(errno));
        int ret = chmod(tmp, S_IRUSR|S_IWUSR|S_IXUSR);
        if ( ret ) error("chmod(%s,S_IRUSR|S_IWUSR|S_IXUSR) failed: %s\n", args->tmp_dir,strerror(errno));
#endif

    fprintf(stderr,"Writing to %s\n", args->tmp_dir);
}
static void destroy(args_t *args)
{
    bcf_hdr_destroy(args->hdr);
    free(args->mem_block);
    free(args->tmp_dir);
    free(args);
}

int main_sort(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->max_mem = 768*1000*1000;
    args->output_fname = "-";
    args->clevel = -1;

    static struct option loptions[] =
    {
        {"max-mem",required_argument,NULL,'m'},
        {"temp-dir",required_argument,NULL,'T'},
        {"output-type",required_argument,NULL,'O'},
        {"output-file",required_argument,NULL,'o'},
        {"output",required_argument,NULL,'o'},
        {"help",no_argument,NULL,'h'},
        {"write-index",optional_argument,NULL,'W'},
        {0,0,0,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "m:T:O:o:W::h?",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'm': args->max_mem = parse_mem_string(optarg); break;
            case 'T': args->tmp_dir = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': args->output_type = FT_BCF_GZ; break;
                          case 'u': args->output_type = FT_BCF; break;
                          case 'z': args->output_type = FT_VCF_GZ; break;
                          case 'v': args->output_type = FT_VCF; break;
                          default:
                          {
                              args->clevel = strtol(optarg,&tmp,10);
                              if ( *tmp || args->clevel<0 || args->clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                          }
                      };
                      if ( optarg[1] )
                      {
                          args->clevel = strtol(optarg+1,&tmp,10);
                          if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                      }
                      break;
            case 'W':
                if (!(args->write_index = write_index_parse(optarg)))
                    error("Unsupported index format '%s'\n", optarg);
                break;
            case 'h':
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else usage(args);
    }
    else args->fname = argv[optind];

    init(args);
    sort_blocks(args);
    merge_to_output(args);
    destroy(args);

    return 0;
}
