/*  prob1.h -- mathematical utility functions.

    Copyright (C) 2010, 2011 Broad Institute.
    Copyright (C) 2012, 2013 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#ifndef BCF_PROB1_H
#define BCF_PROB1_H

#include <htslib/vcf.h>
#include "call.h"

typedef struct {
    int n; // Number of samples
    int M; // Total number of chromosomes across all samples (n*2 if all samples are diploid)
    int n1;
    int is_indel;
    uint8_t *ploidy; // haploid or diploid ONLY
    double *q2p, *pdg; // q2p maps from phread scaled to real likelihood, pdg -> P(D|g)
    double *phi; // Probability of seeing k reference alleles
    double *phi_indel;
    double *z, *zswap; // aux for afs
    double *z1, *z2, *phi1, *phi2; // only calculated when n1 is set
    double **hg; // hypergeometric distribution
    double *lf; // log factorial
    double t, t1, t2;
    double *afs, *afs1; // afs: accumulative allele frequency spectrum (AFS); afs1: site posterior distribution
    const int *PL; // point to PL
    int PL_len;
    int cons_llr;       // pair and trio calling
    int64_t cons_gt;
} bcf_p1aux_t;

typedef struct {
    int rank0, perm_rank; // NB: perm_rank is always set to -1 by bcf_p1_cal()
    int ac; // ML alternative allele count
    double f_exp, f_flat, p_ref_folded, p_ref, p_var_folded, p_var;
    double cil, cih;
    double cmp[3], p_chi2, lrt; // used by contrast2()
} bcf_p1rst_t;

typedef struct {
    double p[4];
    double edb, mqb, bqb;   // end distance bias, mapQ bias, baseQ bias
    int mq, depth, is_tested, d[4];
} anno16_t;

#define MC_PTYPE_FULL  1
#define MC_PTYPE_COND2 2
#define MC_PTYPE_FLAT  3

#ifdef __cplusplus
extern "C" {
#endif

    bcf_p1aux_t *bcf_p1_init(int n_smpl, uint8_t *ploidy);
    void bcf_p1_init_prior(bcf_p1aux_t *ma, int type, double theta);
    void bcf_p1_init_subprior(bcf_p1aux_t *ma, int type, double theta);
    void bcf_p1_destroy(bcf_p1aux_t *ma);
    void bcf_p1_set_ploidy(bcf1_t *b, bcf_p1aux_t *ma);
    int bcf_p1_cal(call_t *call, bcf1_t *b, int do_contrast, bcf_p1aux_t *ma, bcf_p1rst_t *rst);
    int bcf_p1_call_gt(const bcf_p1aux_t *ma, double f0, int k, int is_var);
    void bcf_p1_dump_afs(bcf_p1aux_t *ma);
    int bcf_p1_read_prior(bcf_p1aux_t *ma, const char *fn);
    int bcf_p1_set_n1(bcf_p1aux_t *b, int n1);
    void bcf_p1_set_folded(bcf_p1aux_t *p1a); // only effective when set_n1() is not called

    int bcf_em1(call_t *call, const bcf1_t *b, int n1, int flag, double x[10]);

#ifdef __cplusplus
}
#endif

#endif
