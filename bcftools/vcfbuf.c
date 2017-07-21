/* The MIT License

   Copyright (c) 2016 Genome Research Ltd.

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

#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "vcfbuf.h"
#include "rbuf.h"

typedef struct
{
    double max;
    int rand_missing, skip_filter;
}
ld_t;

typedef struct
{
    bcf1_t *rec;
    double af;
    int af_set:1, idx:31;
}
vcfrec_t;

typedef struct
{
    int max_sites, mvrec, mac, mfarr;
    int *ac, *idx;
    float *farr;
    char *af_tag;
    vcfrec_t **vrec;
}
prune_t;

typedef struct
{
    int active, rid, end;
}
overlap_t;

struct _vcfbuf_t
{
    int win;
    bcf_hdr_t *hdr;
    vcfrec_t *vcf;
    rbuf_t rbuf;
    ld_t ld;
    prune_t prune;
    overlap_t overlap;
};

vcfbuf_t *vcfbuf_init(bcf_hdr_t *hdr, int win)
{
    vcfbuf_t *buf = (vcfbuf_t*) calloc(1,sizeof(vcfbuf_t));
    buf->hdr = hdr;
    buf->win = win;
    buf->overlap.rid = -1;
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
    free(buf->prune.idx);
    free(buf);
}

void vcfbuf_set(vcfbuf_t *buf, vcfbuf_opt_t key, void *value)
{
    if ( key==VCFBUF_LD_MAX ) { buf->ld.max = *((double*)value); return; }
    if ( key==VCFBUF_RAND_MISSING ) { buf->ld.rand_missing = *((int*)value); return; }
    if ( key==VCFBUF_SKIP_FILTER ) { buf->ld.skip_filter = *((int*)value); return; }
    if ( key==VCFBUF_NSITES ) { buf->prune.max_sites = *((int*)value); return; }
    if ( key==VCFBUF_AF_TAG ) { buf->prune.af_tag = *((char**)value); return; }
    if ( key==VCFBUF_OVERLAP_WIN ) { buf->overlap.active = *((int*)value); return; }
}

int vcfbuf_nsites(vcfbuf_t *buf)
{
    return buf->rbuf.n;
}

bcf1_t *vcfbuf_push(vcfbuf_t *buf, bcf1_t *rec, int swap)
{
    if ( !swap ) error("todo: swap=%d\n", swap);

    rbuf_expand0(&buf->rbuf, vcfrec_t, buf->rbuf.n+1, buf->vcf);

    int i = rbuf_append(&buf->rbuf);
    if ( !buf->vcf[i].rec ) buf->vcf[i].rec = bcf_init1();
    
    bcf1_t *ret = buf->vcf[i].rec;
    buf->vcf[i].rec = rec;
    buf->vcf[i].af_set = 0;

    return ret;
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

    if ( nbuf > buf->prune.mvrec )
    {
        buf->prune.idx   = (int*) realloc(buf->prune.idx, nbuf*sizeof(int));
        buf->prune.vrec  = (vcfrec_t**) realloc(buf->prune.vrec, nbuf*sizeof(vcfrec_t*));
        buf->prune.mvrec = nbuf;
    }

    // set allele frequency and prepare buffer for sorting
    int i,k,irec = 0;
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
    int nprune = nbuf - buf->prune.max_sites;
    for (i=0; i<nprune; i++)
        buf->prune.idx[i] = buf->prune.vrec[i]->idx;

    qsort(buf->prune.idx, nprune, sizeof(int), cmpint_desc);

    for (i=0; i<nprune; i++)
        rbuf_remove_kth(&buf->rbuf, vcfrec_t, buf->prune.idx[i], buf->vcf);
}

static int _overlap_can_flush(vcfbuf_t *buf, int flush_all)
{
    if ( flush_all ) { buf->overlap.rid = -1; return 1; }

    int i = rbuf_last(&buf->rbuf);
    vcfrec_t *last = &buf->vcf[i];
    if ( buf->overlap.rid != last->rec->rid ) buf->overlap.end = 0;

    int beg_pos = last->rec->pos;
    int end_pos = last->rec->pos + last->rec->rlen - 1;

    // Assuming left-aligned indels. In case it is a deletion, the real variant
    // starts one base after. If an insertion, the overlap with previous zero length.
    int imin = last->rec->rlen;
    for (i=0; i<last->rec->n_allele; i++)
    {
        char *ref = last->rec->d.allele[0];
        char *alt = last->rec->d.allele[i];
        if ( *alt == '<' ) continue;    // ignore symbolic alleles
        while ( *ref && *alt && nt_to_upper(*ref)==nt_to_upper(*alt) ) { ref++; alt++; }
        if ( imin > ref - last->rec->d.allele[0] ) imin = ref - last->rec->d.allele[0];
    }

    if ( beg_pos <= buf->overlap.end )
    {
        beg_pos += imin;
        if ( beg_pos > end_pos ) end_pos = beg_pos;
    }

    if ( buf->rbuf.n==1 )
    {
        buf->overlap.rid = last->rec->rid;
        buf->overlap.end = end_pos;
        return 0; 
    }
    if ( beg_pos <= buf->overlap.end )
    {
        if ( buf->overlap.end < end_pos ) buf->overlap.end = end_pos;
        return 0;
    }
    return 1;
}

bcf1_t *vcfbuf_flush(vcfbuf_t *buf, int flush_all)
{
    int i,j;

    if ( buf->rbuf.n==0 ) return NULL;
    if ( flush_all ) goto ret;

    i = rbuf_kth(&buf->rbuf, 0);    // first
    j = rbuf_last(&buf->rbuf);      // last

    if ( buf->vcf[i].rec->rid != buf->vcf[j].rec->rid ) goto ret;
    if ( buf->overlap.active )
    {
        int ret = _overlap_can_flush(buf, flush_all);
        //printf("can_flush: %d  %d - %d\n", ret, buf->vcf[i].rec->pos+1, buf->vcf[j].rec->pos+1);
        if ( ret ) goto ret;
    }
    //if ( buf->overlap.active && _overlap_can_flush(buf, flush_all) ) goto ret;

    if ( buf->win > 0 )
    {
        if ( buf->rbuf.n <= buf->win ) return NULL;
        goto ret;
    }
    else if ( buf->win < 0 )
    {
        if ( buf->vcf[i].rec->pos - buf->vcf[j].rec->pos > buf->win ) return NULL;
    }
    else return NULL;
    
ret:
    if ( buf->prune.max_sites && buf->prune.max_sites < buf->rbuf.n ) _prune_sites(buf, flush_all);

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
    For unphased genotypes D is approximated as suggested in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2710162/
        D =~ (GT correlation) * sqrt(Pa*(1-Pa)*Pb*(1-Pb))
*/
static double _calc_ld(vcfbuf_t *buf, bcf1_t *arec, bcf1_t *brec)
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

    // Calculate correlation 
    double ab = 0, aa = 0, bb = 0, a = 0, b = 0;
    int nab = 0, na = 0, nb = 0, ndiff = 0;
    for (i=0; i<arec->n_sample; i++)
    {
        int8_t *aptr = (int8_t*) (afmt->p + i*afmt->size);
        int8_t *bptr = (int8_t*) (bfmt->p + i*bfmt->size);
        int adsg = 0, bdsg = 0, an = 0, bn = 0;
        for (j=0; j<afmt->n; j++)
        {
            if ( aptr[j]==bcf_int8_vector_end ) break;
            if ( aptr[j]==bcf_gt_missing )
            {
                if ( !buf->ld.rand_missing ) break;
                if ( rand()/RAND_MAX >= aaf ) adsg += 1;
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
                if ( rand()/RAND_MAX >= baf ) bdsg += 1;
            }
            else if ( bcf_gt_allele(bptr[j]) ) bdsg += 1;
            bn++;
        }
        if ( an )
        {
            aa += adsg*adsg;
            a  += adsg;
            na++;
        }
        if ( bn )
        {
            bb += bdsg*bdsg;
            b  += bdsg;
            nb++;
        }
        if ( an && bn )
        {
            if ( adsg!=bdsg ) ndiff++;
            ab += adsg*bdsg;
            nab++;
        }
    }
    if ( !nab ) return -1;

    double cor;
    if ( !ndiff ) cor = 1;
    else
    {
        // Don't know how to deal with zero variance. Since this the purpose is filtering,
        // it is not enough to say the value is undefined. Therefore an artificial noise is
        // added to make the denominator non-zero.
        if ( aa == a*a/na || bb == b*b/nb )
        {
            aa += 3*3;
            bb += 3*3;
            ab += 3*3;
            a  += 3;
            b  += 3;
            na++;
            nb++;
            nab++;
        }
        cor = (ab/nab - a/na*b/nb) / sqrt(aa/na - a/na*a/na) / sqrt(bb/nb - b/nb*b/nb);
    }
    return cor*cor;
}

bcf1_t *vcfbuf_max_ld(vcfbuf_t *buf, bcf1_t *rec, double *ld)
{
    *ld = -1;
    if ( !buf->rbuf.n ) return NULL;

    int i = buf->rbuf.f;

    // Relying on vcfbuf being properly flushed - all sites in the buffer
    // must come from the same chromosome
    if ( buf->vcf[i].rec->rid != rec->rid ) return NULL;

    int imax = 0;
    double max = 0;
    for (i=-1; rbuf_next(&buf->rbuf,&i); )
    {   
        if ( buf->ld.skip_filter )
        {
            if ( buf->vcf[i].rec->d.n_flt > 1 ) continue;   // multiple filters are set
            if ( buf->vcf[i].rec->d.n_flt==1 && buf->vcf[i].rec->d.flt[0]!=0 ) continue;    // not PASS
        }
        double val = _calc_ld(buf, buf->vcf[i].rec, rec);
        if ( buf->ld.max && buf->ld.max < val ) 
        {
            *ld = val;
            return buf->vcf[i].rec;
        }
        if ( val > max )
        {
            max  = val;
            imax = i;
        }
    }
    *ld = max;
    return buf->vcf[imax].rec;
}


