/* The MIT License

   Copyright (c) 2016-2024 Genome Research Ltd.

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

#include <assert.h>
#include <strings.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/hts_os.h>
#include <htslib/kbitset.h>
#include "bcftools.h"
#include "vcfbuf.h"
#include "rbuf.h"

typedef struct
{
    double max[VCFBUF_LD_N];
    int rand_missing, filter1;
}
ld_t;

typedef struct
{
    bcf1_t *rec;
    double af;
    unsigned int af_set:1, filter:1, idx:30;
}
vcfrec_t;

#define PRUNE_MODE_MAX_AF 1
#define PRUNE_MODE_1ST    2
#define PRUNE_MODE_RAND   3
typedef struct
{
    int max_sites, mvrec, mac, mfarr, mode;
    int *ac, *idx;
    float *farr;
    char *af_tag;
    vcfrec_t **vrec;
}
prune_t;


#define MARK_OVERLAP 1
#define MARK_DUP     2
#define MARK_EXPR    3

#define MARK_MISSING_SCALAR 0   // actual value to use
#define MARK_MISSING_MAX_DP 1   // max overlap_t.value scaled by INFO/DP

// temporary internal structure for iterative overlap removal by mark_t.expr
typedef struct
{
    double value;       // the sort value
    int rmme, idx;      // mark for removal, index in vcfbuf_t.rbuf
    int dp;             // with MARK_MISSING_MAX_DP, INFO/DP is used extrapolate missing QUAL
    kbitset_t *bset;    // mark which records it overlaps with, given as 0-based indexes to vcfbuf_t.rbuf
    bcf1_t *rec;
}
overlap_t;
typedef struct
{
    // modes
    int mode;
    char *expr;

    // sites marked according to expr, returned to the caller via vcfbuf_get()
    rbuf_t rbuf;
    uint8_t *mark;
    int last;

    // MARK_OVERLAP
    int overlap_rid, overlap_end;

    // MARK_EXPR
    int nbuf;
    overlap_t *buf, **buf_ptr;
    int missing_expr;       // the value to use when min(QUAL) encounters a missing value
    float missing_value;    // the default missing value
    float max_qual;         // with MARK_MISSING_MAX_DP
    int max_qual_dp;        //
    int ntmpi;              // temporary int array and the allocated memory
    int32_t *tmpi;
}
mark_t;

struct _vcfbuf_t
{
    int win,            // maximum number of sites in the buffer, either number of sites (<0) or bp (<0)
        dummy;          // the caller maintains the buffer via push/peek/flush
    bcf_hdr_t *hdr;
    vcfrec_t *vcf;
    rbuf_t rbuf;
    ld_t ld;
    prune_t prune;
    mark_t mark;
    enum { clean, dirty } status;
};

vcfbuf_t *vcfbuf_init(bcf_hdr_t *hdr, int win)
{
    vcfbuf_t *buf = (vcfbuf_t*) calloc(1,sizeof(vcfbuf_t));
    buf->hdr = hdr;
    buf->win = win;
    buf->status = clean;
    buf->mark.overlap_rid = -1;
    int i;
    for (i=0; i<VCFBUF_LD_N; i++) buf->ld.max[i] = HUGE_VAL;
    rbuf_init(&buf->rbuf, 0);
    return buf;
}

void vcfbuf_destroy(vcfbuf_t *buf)
{
    int i;
    for (i=0; i<buf->rbuf.m; i++)
        if ( buf->vcf[i].rec ) bcf_destroy(buf->vcf[i].rec);
    free(buf->vcf);
    free(buf->prune.farr);
    free(buf->prune.vrec);
    free(buf->prune.ac);
    free(buf->prune.af_tag);
    free(buf->prune.idx);
    free(buf->mark.mark);
    free(buf->mark.expr);
    for (i=0; i<buf->mark.nbuf; i++) kbs_destroy(buf->mark.buf[i].bset);
    free(buf->mark.buf);
    free(buf->mark.buf_ptr);
    free(buf->mark.tmpi);
    free(buf);
}

int vcfbuf_set(vcfbuf_t *buf, vcfbuf_opt_t key, ...)
{
    va_list args;
    switch (key)
    {
        case LD_FILTER1:
            va_start(args, key);
            buf->ld.filter1 = va_arg(args,int);
            va_end(args);
            return 0;

        case LD_RAND_MISSING:
            va_start(args, key);
            buf->ld.rand_missing = va_arg(args,int);
            va_end(args);
            return 0;

        case LD_MAX_R2:
            va_start(args, key);
            buf->ld.max[VCFBUF_LD_IDX_R2] = va_arg(args,double);
            va_end(args);
            return 0;

        case LD_MAX_LD:
            va_start(args, key);
            buf->ld.max[VCFBUF_LD_IDX_LD] = va_arg(args,double);
            va_end(args);
            return 0;

        case LD_MAX_HD:
            va_start(args, key);
            buf->ld.max[VCFBUF_LD_IDX_HD] = va_arg(args,double);
            va_end(args);
            return 0;

        case VCFBUF_DUMMY:
            va_start(args, key);
            buf->dummy = va_arg(args,int);
            va_end(args);
            return 0;

        case PRUNE_NSITES:
            va_start(args, key);
            buf->prune.max_sites = va_arg(args,int);
            if ( !buf->prune.mode ) buf->prune.mode = PRUNE_MODE_MAX_AF;
            va_end(args);
            return 0;

        case PRUNE_NSITES_MODE:
            va_start(args, key);
            char *mode = va_arg(args,char*);
            va_end(args);
            if ( !strcasecmp(mode,"maxAF") ) buf->prune.mode = PRUNE_MODE_MAX_AF;
            else if ( !strcasecmp(mode,"1st") ) buf->prune.mode = PRUNE_MODE_1ST;
            else if ( !strcasecmp(mode,"rand") ) buf->prune.mode = PRUNE_MODE_RAND;
            else error("The mode \"%s\" is not recognised\n",mode);
            return 0;

        case PRUNE_AF_TAG:
            va_start(args, key);
            buf->prune.af_tag = strdup(va_arg(args,char*));
            va_end(args);
            return 0;

        case MARK:
            va_start(args, key);
            buf->mark.expr = strdup(va_arg(args,char*));
            if ( !strcasecmp(buf->mark.expr,"overlap") ) buf->mark.mode = MARK_OVERLAP;
            else if ( !strcasecmp(buf->mark.expr,"dup") ) buf->mark.mode = MARK_DUP;
            else buf->mark.mode = MARK_EXPR;
            va_end(args);
            return 0;

        case MARK_MISSING_EXPR:
            va_start(args, key);
            char *expr = va_arg(args,char*);
            if ( !strcasecmp(expr,"0") )
            {
                buf->mark.missing_expr = MARK_MISSING_SCALAR;
                buf->mark.missing_value = 0;
            }
            else if ( !strcasecmp(expr,"DP") )
            {
                if ( buf->mark.mode!=MARK_EXPR ) error("Only the combination of --mark 'min(QUAL)' with --missing DP is currently supported\n");
                buf->mark.missing_expr = MARK_MISSING_MAX_DP;
            }
            else
                error("todo: MARK_MISSING_EXPR=%s\n",expr);
            va_end(args);
            return 0;
    }
    return 0;
}

void *vcfbuf_get(vcfbuf_t *buf, vcfbuf_opt_t key, ...)
{
    va_list args;
    va_start(args, key);
    if ( key==MARK )
        return &buf->mark.last;
    va_end(args);
    return NULL;
}

int vcfbuf_nsites(vcfbuf_t *buf)
{
    return buf->rbuf.n;
}

bcf1_t *vcfbuf_push(vcfbuf_t *buf, bcf1_t *rec)
{
    // make sure the caller is using the buffer correctly and calls vcfbuf_flush()
    // before placing next vcfbuf_push() call
    assert(buf->status!=dirty);
    if ( !buf->dummy ) buf->status = dirty;

    rbuf_expand0(&buf->rbuf, vcfrec_t, buf->rbuf.n+1, buf->vcf);
    int i = rbuf_append(&buf->rbuf);
    if ( !buf->vcf[i].rec ) buf->vcf[i].rec = bcf_init1();

    bcf1_t *ret = buf->vcf[i].rec;
    buf->vcf[i].rec = rec;
    buf->vcf[i].af_set = 0;
    buf->vcf[i].filter = buf->ld.filter1;
    buf->ld.filter1 = 0;

    return ret;
}

bcf1_t *vcfbuf_peek(vcfbuf_t *buf, int idx)
{
    buf->status = clean;
    int i = rbuf_kth(&buf->rbuf, idx);
    return i<0 ? NULL : buf->vcf[i].rec;
}

bcf1_t *vcfbuf_remove(vcfbuf_t *buf, int idx)
{
    int i = rbuf_kth(&buf->rbuf, idx);
    if ( i<0 ) return NULL;
    bcf1_t *rec = buf->vcf[i].rec;
	rbuf_remove_kth(&buf->rbuf, vcfrec_t, idx, buf->vcf);
    return rec;
}

static int cmpvrec(const void *_a, const void *_b)
{
    vcfrec_t *a = *((vcfrec_t**) _a);
    vcfrec_t *b = *((vcfrec_t**) _b);
    if ( a->af < b->af ) return -1;
    if ( a->af == b->af ) return 0;
    return 1;
}
static int cmpint_desc(const void *_a, const void *_b)
{
    int a = *((int*)_a);
    int b = *((int*)_b);
    if ( a < b ) return 1;
    if ( a == b ) return 0;
    return -1;
}

static void _prune_sites(vcfbuf_t *buf, int flush_all)
{

    int nbuf = flush_all ? buf->rbuf.n : buf->rbuf.n - 1;

    int nprune = nbuf - buf->prune.max_sites;
    int i,k,irec = 0;
    if ( buf->prune.mode==PRUNE_MODE_1ST )
    {
        int eoff = flush_all ? 1 : 2;
        for (i=0; i<nprune; i++)
            rbuf_remove_kth(&buf->rbuf, vcfrec_t, buf->rbuf.n - eoff, buf->vcf);
        return;
    }
    if ( buf->prune.mode==PRUNE_MODE_RAND )
    {
        int eoff = flush_all ? 0 : 1;
        for (i=0; i<nprune; i++)
        {
            int j = (buf->rbuf.n - eoff) * hts_drand48();
            rbuf_remove_kth(&buf->rbuf, vcfrec_t, j, buf->vcf);
        }
        return;
    }

    if ( nbuf > buf->prune.mvrec )
    {
        buf->prune.idx   = (int*) realloc(buf->prune.idx, nbuf*sizeof(int));
        buf->prune.vrec  = (vcfrec_t**) realloc(buf->prune.vrec, nbuf*sizeof(vcfrec_t*));
        buf->prune.mvrec = nbuf;
    }

    // set allele frequency and prepare buffer for sorting
    for (i=-1; rbuf_next(&buf->rbuf,&i) && irec<nbuf; )
    {
        bcf1_t *line = buf->vcf[i].rec;
        if ( line->n_allele > buf->prune.mac )
        {
            buf->prune.ac = (int*) realloc(buf->prune.ac, line->n_allele*sizeof(*buf->prune.ac));
            buf->prune.mac = line->n_allele;
        }
        if ( !buf->vcf[i].af_set )
        {
            buf->vcf[i].af = 0;
            if ( buf->prune.af_tag )
            {
                if ( bcf_get_info_float(buf->hdr,line,buf->prune.af_tag,&buf->prune.farr, &buf->prune.mfarr) > 0 ) buf->vcf[i].af = buf->prune.farr[0];
            }
            else if ( bcf_calc_ac(buf->hdr, line, buf->prune.ac, BCF_UN_INFO|BCF_UN_FMT) )
            {
                int ntot = buf->prune.ac[0], nalt = 0;
                for (k=1; k<line->n_allele; k++) nalt += buf->prune.ac[k];
                buf->vcf[i].af = ntot ? (float)nalt/ntot : 0;
            }
            buf->vcf[i].af_set = 1;
        }
        buf->vcf[i].idx = irec;
        buf->prune.vrec[irec++] = &buf->vcf[i];
    }

    // sort by allele frequency, low AF will be removed preferentially
    qsort(buf->prune.vrec, nbuf, sizeof(*buf->prune.vrec), cmpvrec);

    // sort the rbuf indexes to be pruned descendently so that j-th rbuf index
    // is removed before i-th index if i<j
    for (i=0; i<nprune; i++)
        buf->prune.idx[i] = buf->prune.vrec[i]->idx;

    qsort(buf->prune.idx, nprune, sizeof(int), cmpint_desc);

    for (i=0; i<nprune; i++)
        rbuf_remove_kth(&buf->rbuf, vcfrec_t, buf->prune.idx[i], buf->vcf);
}

static int mark_dup_can_flush_(vcfbuf_t *buf, int flush_all)
{
    int flush = flush_all;
    mark_t *mark = &buf->mark;
    if  ( buf->status==dirty )
    {
        // a new site was just added by vcfbuf_push()
        rbuf_expand0(&mark->rbuf, uint8_t, buf->rbuf.n, mark->mark);
        int i = rbuf_append(&mark->rbuf);
        mark->mark[i] = 0;

        if ( buf->rbuf.n==1 ) goto flush;

        // there is at least one previous site, check if it's a duplicate
        int k1 = rbuf_kth(&buf->rbuf, -1);
        int k2 = rbuf_kth(&buf->rbuf, -2);
        vcfrec_t *rec1 = &buf->vcf[k1];
        vcfrec_t *rec2 = &buf->vcf[k2];

        int is_dup = 1;
        if ( rec1->rec->rid!=rec2->rec->rid ) is_dup = 0;
        else if ( rec1->rec->pos!=rec2->rec->pos ) is_dup = 0;

        if ( is_dup )
        {
            // it is, mark the last two sites as duplicates
            int k1 = rbuf_kth(&mark->rbuf, -1);
            int k2 = rbuf_kth(&mark->rbuf, -2);
            mark->mark[k1] = 1;
            mark->mark[k2] = 1;
            goto flush;
        }

        // the last site is not a duplicate with the previous, all sites but the last one can be flushed
        flush = 1;
    }
    else if ( buf->rbuf.n > 1 ) flush = 1;

flush:
    if ( !flush ) return 0;

    int i = rbuf_shift(&mark->rbuf);
    mark->last = mark->mark[i];
    return 1;
}

static int mark_overlap_helper_(vcfbuf_t *buf, int flush_all)
{
    if ( buf->status!=dirty ) return flush_all;

    int flush = flush_all;
    mark_t *mark = &buf->mark;

    // a new site was just added by vcfbuf_push()
    buf->status = clean;

    rbuf_expand0(&mark->rbuf, uint8_t, buf->rbuf.n, mark->mark);
    int i = rbuf_append(&mark->rbuf);
    mark->mark[i] = 0;

    // determine beg and end of the last record that was just added
    i = rbuf_last(&buf->rbuf);
    vcfrec_t *last = &buf->vcf[i];
    if ( mark->overlap_rid != last->rec->rid ) mark->overlap_end = 0;
    int beg_pos = last->rec->pos;
    int end_pos = last->rec->pos + last->rec->rlen - 1;

    // Assuming left-aligned indels. In case it is a deletion, the real variant
    // starts one base after. If an insertion, the overlap with previous is zero
    int imin = last->rec->rlen;
    for (i=0; i<last->rec->n_allele; i++)
    {
        char *ref = last->rec->d.allele[0];
        char *alt = last->rec->d.allele[i];
        if ( *alt == '<' ) continue;    // ignore symbolic alleles
        while ( *ref && *alt && nt_to_upper(*ref)==nt_to_upper(*alt) ) { ref++; alt++; }
        if ( imin > ref - last->rec->d.allele[0] ) imin = ref - last->rec->d.allele[0];
    }
    if ( beg_pos <= mark->overlap_end )
    {
        // the new site overlaps with the previous
        beg_pos += imin;
        if ( beg_pos > end_pos ) end_pos = beg_pos;
    }
    if ( buf->rbuf.n==1 )
    {
        mark->overlap_rid = last->rec->rid;
        mark->overlap_end = end_pos;
        return flush;
    }
    if ( beg_pos <= mark->overlap_end )
    {
        if ( mark->overlap_end < end_pos ) mark->overlap_end = end_pos;
        int k1 = rbuf_kth(&mark->rbuf, -1);
        int k2 = rbuf_kth(&mark->rbuf, -2);
        mark->mark[k1] = 1;
        mark->mark[k2] = 1;
    }
    else
    {
        if ( mark->overlap_end < end_pos ) mark->overlap_end = end_pos;
        flush = 1;
    }
    return flush;
}


static int mark_overlap_can_flush_(vcfbuf_t *buf, int flush_all)
{
    int flush = flush_all;
    if  ( buf->status==dirty ) flush = mark_overlap_helper_(buf,flush_all);
    else if ( buf->rbuf.n > 1 ) flush = 1;
    if ( !flush ) return 0;

    mark_t *mark = &buf->mark;
    int i = rbuf_shift(&mark->rbuf);
    mark->last = mark->mark[i];
    return 1;
}


static int records_overlap(bcf1_t *a, bcf1_t *b)
{
    if ( a->rid != b->rid ) return 0;
    if ( a->pos + a->rlen - 1 < b->pos ) return 0;
    return 1;
}

static int cmp_overlap_ptr_asc(const void *aptr, const void *bptr)
{
    overlap_t *a = *((overlap_t**)aptr);
    overlap_t *b = *((overlap_t**)bptr);
    if ( a->value < b->value ) return -1;
    if ( a->value > b->value ) return 1;
    return 0;
}
static void mark_expr_missing_reset_(vcfbuf_t *buf)
{
    buf->mark.max_qual = 0;
    buf->mark.max_qual_dp = 0;
}
static void mark_expr_missing_prep_(vcfbuf_t *buf, overlap_t *olap)
{
    int nval = bcf_get_info_int32(buf->hdr,olap->rec,"DP",&buf->mark.tmpi,&buf->mark.ntmpi);
    if ( nval!=1 ) return;

    olap->dp = buf->mark.tmpi[0];
    if ( bcf_float_is_missing(olap->rec->qual) ) return;
    if ( buf->mark.max_qual < olap->rec->qual )
    {
        buf->mark.max_qual = olap->rec->qual;
        buf->mark.max_qual_dp = olap->dp;
    }
}
static void mark_expr_missing_set_(vcfbuf_t *buf, overlap_t *olap)
{
    if ( !bcf_float_is_missing(olap->rec->qual) ) return;
    if ( !buf->mark.max_qual_dp ) return;

    // scale QUAL of the most confident variant in the overlap proportionally to the coverage
    // and use that to prioritize the records
    olap->value = buf->mark.max_qual * olap->dp / buf->mark.max_qual_dp;
}
static int mark_expr_can_flush_(vcfbuf_t *buf, int flush_all)
{
    mark_t *mark = &buf->mark;
    if ( strcasecmp("min(QUAL)",mark->expr) ) error("Todo; at this time only min(QUAL) is supported\n");

    int flush = flush_all;
    if  ( buf->status==dirty )
    {
        flush = mark_overlap_helper_(buf,flush_all);
        if ( !flush ) return 0;

        if ( mark->missing_expr==MARK_MISSING_MAX_DP ) mark_expr_missing_reset_(buf);

        // init overlaps, each overlap_t structure keeps a list of overlapping records, symmetrical
        size_t nori = mark->nbuf;
        hts_resize(overlap_t,  buf->rbuf.n, &mark->nbuf, &mark->buf, HTS_RESIZE_CLEAR);
        hts_resize(overlap_t*, buf->rbuf.n, &nori, &mark->buf_ptr, HTS_RESIZE_CLEAR);
        int i;
        for (i=0; i<buf->rbuf.n; i++)
        {
            overlap_t *oi = &mark->buf[i];
            int j = rbuf_kth(&buf->rbuf, i);
            assert(j>=0);
            bcf1_t *rec = buf->vcf[j].rec;
            assert(rec);
            oi->rec = rec;

            // todo: other than QUAL values
            oi->value = bcf_float_is_missing(rec->qual) ? mark->missing_value : rec->qual;
            if ( mark->missing_expr==MARK_MISSING_MAX_DP ) mark_expr_missing_prep_(buf,oi);
            if ( oi->bset )
            {
                kbs_resize(&oi->bset,buf->rbuf.n);
                kbs_clear(oi->bset);
            }
            else
                oi->bset = kbs_init(buf->rbuf.n);
            oi->idx  = i;
            mark->buf_ptr[i] = oi;
            mark->mark[oi->idx] = 0;
        }
        int nolap = 0;
        for (i=0; i<buf->rbuf.n; i++)
        {
            overlap_t *oi = &mark->buf[i];
            if ( mark->missing_expr==MARK_MISSING_MAX_DP ) mark_expr_missing_set_(buf,oi);
            int j;
            for (j=i+1; j<buf->rbuf.n; j++)
            {
                overlap_t *oj = &mark->buf[j];
                if ( !records_overlap(oi->rec,oj->rec) ) continue;
                kbs_insert(oi->bset,j);
                kbs_insert(oj->bset,i);
                nolap++;
            }
        }

        // sort according to the requested criteria, currently only min(QUAL)
        qsort(mark->buf_ptr,buf->rbuf.n,sizeof(*mark->buf_ptr),cmp_overlap_ptr_asc);   // todo: other than min()

        // go through the list sorted by overlap_t.value, eg QUAL
        for (i=0; nolap && i<buf->rbuf.n; i++)
        {
            kbitset_iter_t itr;
            overlap_t *oi = mark->buf_ptr[i];
            kbs_start(&itr);
            int j;
            while ((j=kbs_next(oi->bset, &itr)) >= 0)
            {
                kbs_delete(oi->bset,j);
                assert(nolap);
                assert(kbs_exists(mark->buf[j].bset,oi->idx));
                kbs_delete(mark->buf[j].bset,oi->idx);
                nolap--;
            }
            j = rbuf_kth(&mark->rbuf,oi->idx);
            mark->mark[j] = 1;
        }
    }
    else if ( buf->rbuf.n > 1 ) flush = 1;
    if ( !flush ) return 0;

    int i = rbuf_shift(&mark->rbuf);
    mark->last = mark->mark[i];
    return 1;
}

bcf1_t *vcfbuf_flush(vcfbuf_t *buf, int flush_all)
{
    int i,j;

    // nothing to do, no lines in the buffer
    if ( buf->rbuf.n==0 ) return NULL;

    // dummy mode, always flushing
    if ( buf->dummy ) goto ret;

    // pruning mode
    if ( buf->win )
    {
        int can_flush = flush_all;
        i = rbuf_kth(&buf->rbuf, 0);    // first
        j = rbuf_last(&buf->rbuf);      // last
        if ( buf->vcf[i].rec->rid != buf->vcf[j].rec->rid ) can_flush = 1;
        else if ( buf->win > 0 )
        {
            if ( buf->rbuf.n > buf->win ) can_flush = 1;
        }
        else if ( buf->win < 0 )
        {
            if ( !(buf->vcf[i].rec->pos - buf->vcf[j].rec->pos > buf->win) ) can_flush = 1;
        }
        buf->status = clean;
        if ( !can_flush ) return NULL;
        if ( buf->prune.max_sites && buf->prune.max_sites < buf->rbuf.n ) _prune_sites(buf, flush_all);
        goto ret;
    }

    // overlaps and duplicates
    if ( buf->mark.mode )
    {
        int can_flush = 0;
        if ( buf->mark.mode==MARK_OVERLAP )
        {
            if ( mark_overlap_can_flush_(buf,flush_all) ) can_flush = 1;
        }
        else if ( buf->mark.mode==MARK_DUP )
        {
            if ( mark_dup_can_flush_(buf,flush_all) ) can_flush = 1;
        }
        if ( buf->mark.mode==MARK_EXPR )
        {
            if ( mark_expr_can_flush_(buf,flush_all) ) can_flush = 1;
        }
        buf->status = clean;
        if ( !can_flush ) return NULL;
        goto ret;
    }

ret:
    buf->status = clean;
    i = rbuf_shift(&buf->rbuf);
    return buf->vcf[i].rec;
}

static double _estimate_af(int8_t *ptr, int size, int nvals, int nsamples)
{
    int i,j, nref = 0, nalt = 0;
    for (i=0; i<nsamples; i++)
    {
        for (j=0; j<nvals; j++)
        {
            if ( ptr[j]==bcf_gt_missing ) break;
            if ( ptr[j]==bcf_int8_vector_end ) break;
            if ( bcf_gt_allele(ptr[j]) ) nalt++;
            else nref++;
        }
        ptr += size;
    }
    if ( nref+nalt == 0 ) return 0;
    return (double)nalt/(nref+nalt);
}

/*
    The `ld` is set to D approximated as suggested in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2710162/
        D =~ (GT correlation) * sqrt(Pa*(1-Pa)*Pb*(1-Pb))

    and `hd` as proposed in Ragsdale, A. P., & Gravel, S. (2019). Unbiased estimation of linkage
    disequilibrium from unphased data.  Molecular Biology and Evolution. doi:10.1093/molbev/msz265

        \hat{D} = 1/[n*(n+1)]*[
                             (n1 + n2/2 + n4/2 + n5/4)*(n5/4 + n6/2 + n8/2 + n9)
                            -(n2/2 + n3 + n5/4 + n6/2)*(n4/2 + n5/4 + n7 + n8/2)
                        ]
    where n1,n2,..n9 are counts of RR/RR,RR/RA,..,AA/AA genotypes.

    Returns 0 on success, -1 if the values could not be determined (missing genotypes)
*/
static int _calc_r2_ld(vcfbuf_t *buf, bcf1_t *arec, bcf1_t *brec, vcfbuf_ld_t *ld)
{
    if ( arec->n_sample!=brec->n_sample ) error("Different number of samples: %d vs %d\n",arec->n_sample,brec->n_sample);
    assert( arec->n_sample );

    int i,j,igt = bcf_hdr_id2int(buf->hdr, BCF_DT_ID, "GT");
    bcf_unpack(arec, BCF_UN_FMT);
    bcf_unpack(brec, BCF_UN_FMT);
    bcf_fmt_t *afmt = NULL, *bfmt = NULL;
    for (i=0; i<arec->n_fmt; i++)
        if ( arec->d.fmt[i].id==igt ) { afmt = &arec->d.fmt[i]; break; }
    if ( !afmt ) return -1;  // no GT tag
    for (i=0; i<brec->n_fmt; i++)
        if ( brec->d.fmt[i].id==igt ) { bfmt = &brec->d.fmt[i]; break; }
    if ( !bfmt ) return -1;  // no GT tag

    if ( afmt->n==0 ) return -1;   // empty?!
    if ( bfmt->n==0 ) return -1;   // empty?!
    if ( afmt->type!=BCF_BT_INT8 ) error("TODO: the GT fmt_type is not int8!\n");
    if ( bfmt->type!=BCF_BT_INT8 ) error("TODO: the GT fmt_type is not int8!\n");

    // Determine allele frequencies, this is to sample randomly missing genotypes
    double aaf = 0, baf = 0;
    if ( buf->ld.rand_missing )
    {
        aaf = _estimate_af((int8_t*)afmt->p, afmt->size, afmt->n, arec->n_sample);
        baf = _estimate_af((int8_t*)bfmt->p, bfmt->size, bfmt->n, brec->n_sample);
    }

    // Calculate r2, lf, hd
    double nhd[] = {0,0,0,0,0,0,0,0,0};
    double ab = 0, aa = 0, bb = 0, a = 0, b = 0;
    int nab = 0, ndiff = 0;
    int an_tot = 0, bn_tot = 0;
    for (i=0; i<arec->n_sample; i++)
    {
        int8_t *aptr = (int8_t*) (afmt->p + i*afmt->size);
        int8_t *bptr = (int8_t*) (bfmt->p + i*bfmt->size);
        int adsg = 0, bdsg = 0;     // dosages (0,1,2) at sites (a,b)
        int an = 0, bn = 0;         // number of alleles at sites (a,b)
        for (j=0; j<afmt->n; j++)
        {
            if ( aptr[j]==bcf_int8_vector_end ) break;
            if ( aptr[j]==bcf_gt_missing )
            {
                if ( !buf->ld.rand_missing ) break;
                if ( hts_drand48() >= aaf ) adsg += 1;
            }
            else if ( bcf_gt_allele(aptr[j]) ) adsg += 1;
            an++;
        }
        for (j=0; j<bfmt->n; j++)
        {
            if ( bptr[j]==bcf_int8_vector_end ) break;
            if ( bptr[j]==bcf_gt_missing )
            {
                if ( !buf->ld.rand_missing ) break;
                if ( hts_drand48() >= baf ) bdsg += 1;
            }
            else if ( bcf_gt_allele(bptr[j]) ) bdsg += 1;
            bn++;
        }
        if ( an && bn )
        {
            an_tot += an;
            aa += adsg*adsg;
            a  += adsg;

            bn_tot += bn;
            bb += bdsg*bdsg;
            b  += bdsg;

            if ( adsg!=bdsg ) ndiff++;
            ab += adsg*bdsg;
            nab++;
        }
        if ( an==2 && bn==2 )   // for now only diploid genotypes
        {
            assert( adsg<=2 && bdsg<=2 );
            nhd[ bdsg*3 + adsg ]++;
        }
    }
    if ( !nab ) return -1;  // no data in common for the two sites

    double pa = a/an_tot;
    double pb = b/bn_tot;
    double cor;
    if ( !ndiff ) cor = 1;
    else
    {
        if ( aa == a*a/nab || bb == b*b/nab )     // zero variance, add small noise
        {
            aa += 1e-4;
            bb += 1e-4;
            ab += 1e-4;
            a  += 1e-2;
            b  += 1e-2;
            nab++;
        }
        cor = (ab - a*b/nab) / sqrt(aa - a*a/nab) / sqrt(bb - b*b/nab);
    }

    ld->val[VCFBUF_LD_IDX_R2] = cor * cor;

    // Lewontin's normalization of D. Also we cap at 1 as the calculation
    // can result in values bigger than 1 for high AFs.
    ld->val[VCFBUF_LD_IDX_LD] = cor * sqrt(pa*(1-pa)*pb*(1-pb));
    double norm;
    if ( ld->val[VCFBUF_LD_IDX_LD] < 0 )
        norm = -pa*pb > -(1-pa)*(1-pb) ? -pa*pb : -(1-pa)*(1-pb);
    else
        norm = pa*(1-pb) > (1-pa)*pb ? pa*(1-pb) : (1-pa)*pb;
    if ( norm )
        ld->val[VCFBUF_LD_IDX_LD] = fabs(norm) > fabs(ld->val[VCFBUF_LD_IDX_LD]) ? ld->val[VCFBUF_LD_IDX_LD]/norm : 1;
    if ( !ld->val[VCFBUF_LD_IDX_LD] )
        ld->val[VCFBUF_LD_IDX_LD] = fabs(ld->val[VCFBUF_LD_IDX_LD]);    // avoid "-0" on output

    ld->val[VCFBUF_LD_IDX_HD] =
        (nhd[0] + nhd[1]/2. + nhd[3]/2. + nhd[4]/4.)*(nhd[4]/4. + nhd[5]/2. + nhd[7]/2. + nhd[8])
        - (nhd[1]/2. + nhd[2] + nhd[4]/4. + nhd[5]/2.)*(nhd[3]/2. + nhd[4]/4. + nhd[6] + nhd[7]/2.);
    ld->val[VCFBUF_LD_IDX_HD] /= nab;
    ld->val[VCFBUF_LD_IDX_HD] /= nab+1;

    return 0;
}

int vcfbuf_ld(vcfbuf_t *buf, bcf1_t *rec, vcfbuf_ld_t *ld)
{
    int ret = -1;
    if ( !buf->rbuf.n ) return ret;

    int j, i = buf->rbuf.f;

    // Relying on vcfbuf being properly flushed - all sites in the buffer
    // must come from the same chromosome
    if ( buf->vcf[i].rec->rid != rec->rid ) return ret;

    vcfbuf_ld_t tmp;
    for (j=0; j<VCFBUF_LD_N; j++)
    {
        ld->val[j] = -HUGE_VAL;
        ld->rec[j] = NULL;
    }

    for (i=-1; rbuf_next(&buf->rbuf,&i); )
    {
        if ( buf->vcf[i].filter ) continue;
        if ( _calc_r2_ld(buf, buf->vcf[i].rec, rec, &tmp) < 0 ) continue;   // missing genotypes

        int done = 0;
        for (j=0; j<VCFBUF_LD_N; j++)
        {
            if ( ld->val[j] < tmp.val[j] )
            {
                ld->val[j] = tmp.val[j];
                ld->rec[j] = buf->vcf[i].rec;
            }
            if ( buf->ld.max[j] < tmp.val[j] ) done = 1;
            ret = 0;
        }
        if ( done ) return ret;
    }
    return ret;
}


