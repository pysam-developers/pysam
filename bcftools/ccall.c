/*  ccall.c -- consensus variant calling.

    Copyright (C) 2013-2014 Genome Research Ltd.
    Portions copyright (C) 2010 Broad Institute.

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
THE SOFTWARE.  */

#include <math.h>
#include <htslib/kfunc.h>
#include "call.h"
#include "kmin.h"
#include "prob1.h"

// Most of the original -c calling was moved to bcftools as it was
// and its data structures were wrapped into the ccal_t to make it
// functional quickly. This is not the desired state.
struct _ccall_t
{
    bcf_p1aux_t *p1;
};

void ccall_init(call_t *call)
{
    call->cdat = (ccall_t*) calloc(1,sizeof(ccall_t));
    call_init_pl2p(call);
    call->cdat->p1 = bcf_p1_init(bcf_hdr_nsamples(call->hdr), call->ploidy);
    call->gts = (int*) calloc(bcf_hdr_nsamples(call->hdr)*2,sizeof(int));   // assuming at most diploid everywhere
    call->nals_map = 5;
    call->als_map  = (int*) malloc(sizeof(int)*call->nals_map);

    bcf_hdr_append(call->hdr,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    if ( call->output_tags & CALL_FMT_GQ )
    {
        bcf_hdr_append(call->hdr,"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
        call->GQs = (int32_t*) malloc(sizeof(int32_t)*bcf_hdr_nsamples(call->hdr));
    }
    if ( call->output_tags & CALL_FMT_GP )
        error("Sorry, -f GP is not supported with -c\n");
    bcf_hdr_append(call->hdr,"##INFO=<ID=AF1,Number=1,Type=Float,Description=\"Max-likelihood estimate of the first ALT allele frequency (assuming HWE)\">");
    // Todo: groups not migrated to 'bcftools call' yet
    bcf_hdr_append(call->hdr,"##INFO=<ID=AF2,Number=1,Type=Float,Description=\"Max-likelihood estimate of the first and second group ALT allele frequency (assuming HWE)\">");
    bcf_hdr_append(call->hdr,"##INFO=<ID=AC1,Number=1,Type=Float,Description=\"Max-likelihood estimate of the first ALT allele count (no HWE assumption)\">");
    bcf_hdr_append(call->hdr,"##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Root-mean-square mapping quality of covering reads\">\n");
    bcf_hdr_append(call->hdr,"##INFO=<ID=FQ,Number=1,Type=Float,Description=\"Phred probability of all samples being the same\">\n");
    bcf_hdr_append(call->hdr,"##INFO=<ID=PV4,Number=4,Type=Float,Description=\"P-values for strand bias, baseQ bias, mapQ bias and tail distance bias\">\n");
    bcf_hdr_append(call->hdr,"##INFO=<ID=G3,Number=3,Type=Float,Description=\"ML estimate of genotype frequencies\">\n");
    bcf_hdr_append(call->hdr,"##INFO=<ID=HWE,Number=1,Type=Float,Description=\"Chi^2 based HWE test P-value based on G3\">\n");
    // bcf_hdr_append(call->hdr,);
    // bcf_hdr_append(call->hdr,);
    bcf_hdr_append(call->hdr,"##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases\">");

    return;
}
void ccall_destroy(call_t *call)
{
    free(call->itmp);
    free(call->als_map);
    free(call->gts);
    free(call->anno16);
    free(call->PLs);
    free(call->GQs);
    free(call->pdg);
    bcf_p1_destroy(call->cdat->p1);
    free(call->cdat);
    return;
}

// Inits P(D|G): convert PLs from log space, only two alleles (three PLs) are used.
// NB: The original samtools calling code uses pdgs in reverse order (AA comes
// first, RR last), while the -m calling model uses the canonical order.
static void set_pdg3(double *pl2p, int *PLs, double *pdg, int n_smpl, int n_gt)
{
    int i;
    for (i=0; i<n_smpl; i++)
    {
        pdg[2] = pl2p[ PLs[0] ];
        pdg[1] = pl2p[ PLs[1] ];
        pdg[0] = pl2p[ PLs[2] ];
        PLs += n_gt;
        pdg += 3;
    }
}

static double ttest(int n1, int n2, float a[4])
{
    extern double kf_betai(double a, double b, double x);
    double t, v, u1, u2;
    if (n1 == 0 || n2 == 0 || n1 + n2 < 3) return 1.0;
    u1 = (double)a[0] / n1; u2 = (double)a[2] / n2;
    if (u1 <= u2) return 1.;
    t = (u1 - u2) / sqrt(((a[1] - n1 * u1 * u1) + (a[3] - n2 * u2 * u2)) / (n1 + n2 - 2) * (1./n1 + 1./n2));
    v = n1 + n2 - 2;
    return t < 0.? 1. : .5 * kf_betai(.5*v, .5, v/(v+t*t));
}

static int test16_core(float anno[16], anno16_t *a)
{
    extern double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two);
    double left, right;
    int i;
    a->p[0] = a->p[1] = a->p[2] = a->p[3] = 1.;
    for (i=0; i<4; i++) a->d[i] = anno[i];
    a->depth = anno[0] + anno[1] + anno[2] + anno[3];
    a->is_tested = (anno[0] + anno[1] > 0 && anno[2] + anno[3] > 0);
    if (a->depth == 0) return -1;
    a->mq = (int)(sqrt((anno[9] + anno[11]) / a->depth) + .499);
    kt_fisher_exact(anno[0], anno[1], anno[2], anno[3], &left, &right, &a->p[0]);
    for (i = 1; i < 4; ++i)
        a->p[i] = ttest(anno[0] + anno[1], anno[2] + anno[3], anno+4*i);
    return 0;
}

int test16(float *anno16, anno16_t *a)
{
    a->p[0] = a->p[1] = a->p[2] = a->p[3] = 1.;
    a->d[0] = a->d[1] = a->d[2] = a->d[3] = 0.;
    a->mq = a->depth = a->is_tested = 0;
    return test16_core(anno16, a);
}
static int update_bcf1(call_t *call, bcf1_t *rec, const bcf_p1rst_t *pr, double em[10])
{
    int has_I16, is_var;
    float fq, r;
    anno16_t a;
    float tmpf[4], tmpi;

    bcf_get_info_float(call->hdr, rec, "I16", &call->anno16, &call->n16);

    has_I16 = test16(call->anno16, &a) >= 0? 1 : 0;

    // print EM
    if (em[0] >= 0)
    {
        tmpf[0] = 1 - em[0];
        bcf_update_info_float(call->hdr, rec, "AF1", tmpf, 1);
    }
    if (em[4] >= 0 && em[4] <= 0.05)
    {
        tmpf[0] = em[3]; tmpf[1] = em[2]; tmpf[2] = em[1]; tmpf[3] = em[4];
        bcf_update_info_float(call->hdr, rec, "G3", tmpf, 3);
        bcf_update_info_float(call->hdr, rec, "HWE", &tmpf[3], 1);
    }
    if (em[5] >= 0 && em[6] >= 0)
    {
        tmpf[0] = 1 - em[5]; tmpf[1] = 1 - em[6];
        bcf_update_info_float(call->hdr, rec, "AF2", tmpf, 2);
    }
    if (em[7] >= 0)
    {
        tmpf[0] = em[7];
        bcf_update_info_float(call->hdr, rec, "LRT", tmpf, 1);
    }
    if (em[8] >= 0)
    {
        tmpf[0] = em[8];
        bcf_update_info_float(call->hdr, rec, "LRT2", tmpf, 1);
    }

    bcf_p1aux_t *p1 = call->cdat->p1;
    if (p1->cons_llr > 0)
    {
        tmpi = p1->cons_llr;
        bcf_update_info_int32(call->hdr, rec, "CLR", &tmpi, 1);
        // todo: trio calling with -c
        if (p1->cons_gt > 0)
        {
            char tmp[4];
            tmp[0] = p1->cons_gt&0xff; tmp[1] = p1->cons_gt>>8&0xff; tmp[2] = p1->cons_gt>>16&0xff; tmp[3] = 0;
            bcf_update_info_string(call->hdr, rec, "UGT", tmp);
            tmp[0] = p1->cons_gt>>32&0xff; tmp[1] = p1->cons_gt>>40&0xff; tmp[2] = p1->cons_gt>>48&0xff;
            bcf_update_info_string(call->hdr, rec, "CGT", tmp);
        }
    }
    is_var = (pr->p_ref < call->pref);
    r = is_var? pr->p_ref : pr->p_var;

    bcf_update_info_int32(call->hdr, rec, "AC1", &pr->ac, 1);
    int32_t dp[4]; dp[0] = call->anno16[0]; dp[1] = call->anno16[1]; dp[2] = call->anno16[2]; dp[3] = call->anno16[3];
    bcf_update_info_int32(call->hdr, rec, "DP4", dp, 4);
    bcf_update_info_int32(call->hdr, rec, "MQ", &a.mq, 1);

    fq = pr->p_ref_folded < 0.5? -4.343 * log(pr->p_ref_folded) : 4.343 * log(pr->p_var_folded);
    if (fq < -999) fq = -999;
    if (fq > 999) fq = 999;
    bcf_update_info_float(call->hdr, rec, "FQ", &fq, 1);

    assert( pr->cmp[0]<0 );
    // todo
    //  if (pr->cmp[0] >= 0.) { // two sample groups
    //      int i, q[3];
    //      for (i = 1; i < 3; ++i) {
    //          double x = pr->cmp[i] + pr->cmp[0]/2.;
    //          q[i] = x == 0? 255 : (int)(-4.343 * log(x) + .499);
    //          if (q[i] > 255) q[i] = 255;
    //      }
    //      if (pr->perm_rank >= 0) ksprintf(&s, "PR=%d;", pr->perm_rank);
    //
    //      ksprintf(&s, "PCHI2=%.3g;PC2=%d,%d;", q[1], q[2], pr->p_chi2);
    //  }

    if (has_I16 && a.is_tested)
    {
        int i;
        for (i=0; i<4; i++) tmpf[i] = a.p[i];
        bcf_update_info_float(call->hdr, rec, "PV4", tmpf, 4);
    }
    bcf_update_info_int32(call->hdr, rec, "I16", NULL, 0);     // remove I16 tag
    bcf_update_info_int32(call->hdr, rec, "QS", NULL, 0);      // remove QS tag

    rec->qual = r < 1e-100? 999 : -4.343 * log(r);
    if (rec->qual > 999) rec->qual = 999;

    // Remove unused alleles
    int nals_ori = rec->n_allele, nals = !is_var && !(call->flag & CALL_KEEPALT) ? 1 : pr->rank0 < 2? 2 : pr->rank0+1;
    if ( call->flag & CALL_KEEPALT && call->unseen==nals-1 ) nals--;
    
    if ( nals<rec->n_allele )
    {
        bcf_update_alleles(call->hdr, rec, (const char**)rec->d.allele, nals);

        // Update PLs
        int npls_src = call->nPLs / rec->n_sample, npls_dst = nals*(nals+1)/2;
        int *pls_src = call->PLs - npls_src, *pls_dst = call->PLs - npls_dst;
        int isample, i;
        for (isample = 0; isample < rec->n_sample; isample++)
        {
            pls_src += npls_src;
            pls_dst += npls_dst;
            if ( !call->ploidy || call->ploidy[isample]==2 )
            {
                for (i=0; i<npls_dst; i++)
                    pls_dst[i] =  pls_src[i];
            }
            else
            {
                for (i=0; i<nals; i++)
                {
                    int isrc = (i+1)*(i+2)/2-1;
                    pls_dst[i] = pls_src[isrc];
                }
                if (i<npls_dst) pls_dst[i] = bcf_int32_vector_end;
            }
        }
        bcf_update_format_int32(call->hdr, rec, "PL", call->PLs, npls_dst*rec->n_sample);
    }

    // Call genotypes
    int i;
    for (i=0; i<rec->n_sample; i++)
    {
        int x = ( is_var || call->output_tags & CALL_FMT_GQ ) ? bcf_p1_call_gt(p1, pr->f_exp, i, is_var) : 2;
        int gt = x&3;
        if ( !call->ploidy || call->ploidy[i]==2 )
        {
            if ( gt==1 )
            {
                call->gts[2*i]   = bcf_gt_unphased(0);
                call->gts[2*i+1] = bcf_gt_unphased(1);
            }
            else if ( gt==0 )
            {
                call->gts[2*i]   = bcf_gt_unphased(1);
                call->gts[2*i+1] = bcf_gt_unphased(1);
            }
            else
            {
                call->gts[2*i]   = bcf_gt_unphased(0);
                call->gts[2*i+1] = bcf_gt_unphased(0);
            }
            if ( call->output_tags & CALL_FMT_GQ ) call->GQs[i] = x>>2;
        }
        else
        {
            if ( gt==0 ) call->gts[2*i] = bcf_gt_unphased(1);
            else call->gts[2*i] = bcf_gt_unphased(0);
            call->gts[2*i+1] = bcf_int32_vector_end;
            if ( call->output_tags & CALL_FMT_GQ ) call->GQs[i] = bcf_int32_missing;
        }
    }
    bcf_update_genotypes(call->hdr, rec, call->gts, rec->n_sample*2);
    if ( call->output_tags & CALL_FMT_GQ )
        bcf_update_format_int32(call->hdr, rec, "GQ", call->GQs, rec->n_sample);

    // trim Number=R tags
    int out_als = 0;
    for (i=0; i<nals; i++) out_als |= 1<<i;
    init_allele_trimming_maps(call, out_als, nals_ori);
    mcall_trim_numberR(call, rec, nals_ori, nals, out_als);

    return is_var;
}


int ccall(call_t *call, bcf1_t *rec)
{
    int nsmpl = bcf_hdr_nsamples(call->hdr);

    // Get the genotype likelihoods
    int nals = rec->n_allele;
    call->nPLs = bcf_get_format_int32(call->hdr, rec, "PL", &call->PLs, &call->mPLs);
    if ( call->nPLs!=nsmpl*nals*(nals+1)/2 && call->nPLs!=nsmpl*nals )  // diploid+haploid or haploid only
        error("Wrong number of PL fields? nals=%d npl=%d\n", nals,call->nPLs);

    // Convert PLs to probabilities, only first two alleles are considered
    int ngts = nals*(nals+1)/2;
    hts_expand(double, 3*nsmpl, call->npdg, call->pdg);
    set_pdg3(call->pl2p, call->PLs, call->pdg, nsmpl, ngts);

    double em[10] = {-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.};
    int ret = bcf_em1(call, rec, call->ngrp1_samples, 0x1ff, em);

    bcf_p1rst_t pr;
    int do_contrast = (em[7] >= 0 && em[7] < call->min_lrt) ? 1 : 0;
    ret = bcf_p1_cal(call, rec, do_contrast, call->cdat->p1, &pr);
    if (pr.p_ref >= call->pref && (call->flag & CALL_VARONLY)) return 0;
    if (ret >= 0) ret = update_bcf1(call, rec, &pr, em);
    return ret;
}

