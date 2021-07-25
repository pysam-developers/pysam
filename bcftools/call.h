/*  call.h -- variant calling declarations.

    Copyright (C) 2013-2015, 2019-2020 Genome Research Ltd.

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

#ifndef __CALL_H__
#define __CALL_H__

#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include "vcmp.h"

#define CALL_KEEPALT        1
#define CALL_VARONLY        (1<<1)
#define CALL_CONSTR_TRIO    (1<<2)
#define CALL_CONSTR_ALLELES (1<<3)
//
#define CALL_FMT_PV4        (1<<5)
#define CALL_FMT_GQ         (1<<6)
#define CALL_FMT_GP         (1<<7)

#define FATHER 0
#define MOTHER 1
#define CHILD  2
typedef struct
{
    char *name;
    int sample[3];  // father, mother, child
    int type;       // see FTYPE_* definitions in mcall.c
}
family_t;

// For the single-sample and grouped -G calling
typedef struct
{
    double ref_lk, max_lk, lk_sum;
    float *qsum;    // QS(quality sum) values
    int nqsum;
    uint32_t *smpl, nsmpl;
    uint32_t nals, als;
}
smpl_grp_t;

// For the `-C alleles -i` constrained calling
typedef struct
{
    uint32_t n:31, used:1;
    char **allele;
}
tgt_als_t;

typedef struct _ccall_t ccall_t;
typedef struct
{
    // mcall only
    int npdg;
    int *als_map, nals_map; // mapping from full set of alleles to trimmed set of alleles (old -> new)
    int *pl_map, npl_map;   // same as above for PLs, but reverse (new -> old)
    char **als;             // array to hold the trimmed set of alleles to appear on output
    int nals;               // size of the als array
    int als_new, nals_new;  // bitmask with final alleles and their number
    family_t *fams;         // list of families and samples for trio calling
    int nfams, mfams;
    int ntrio[5][5];        // possible trio genotype combinations and their counts; first idx:
    uint16_t *trio[5][5];   //  family type, second index: allele count (2-4, first two are unused)
    double *GLs;
    float *GPs;             // FORMAT/GP: posterior probabilities
    int32_t *GQs, *ADs;     // FORMAT/GQ: genotype qualities; AD: allelic depth for -G
    int32_t *itmp;          // temporary int array, used for new PLs with CALL_CONSTR_ALLELES
    int n_itmp, nGPs, nADs;
    vcmp_t *vcmp;
    double trio_Pm_SNPs, trio_Pm_del, trio_Pm_ins;      // P(mendelian) for trio calling, see mcall_call_trio_genotypes()
    int32_t *ugts, *cgts;   // unconstraind and constrained GTs
    uint32_t output_tags;
    char *prior_AN, *prior_AC;  // reference panel AF tags (AF=AC/AN)
    tgt_als_t *tgt_als;         // for CALL_CONSTR_ALLELES
    char *sample_groups;        // for single-sample or grouped calling with -G
    char *sample_groups_tag;    // for -G [AD|QS:]
    smpl_grp_t *smpl_grp;
    int nsmpl_grp;

    // ccall only
    double indel_frac, min_perm_p, min_lrt;
    double prior_type, pref;
    int ngrp1_samples, n_perm;
    char *prior_file;
    ccall_t *cdat;

    // shared
    bcf_srs_t *srs;         // BCF synced readers holding target alleles for CALL_CONSTR_ALLELES
    bcf1_t *rec;
    bcf_hdr_t *hdr;
    uint32_t flag;          // One or more of the CALL_* flags defined above
    uint8_t *ploidy, all_diploid, unseen;

    double pl2p[256];       // PL to 10^(-PL/10) table
    int32_t *PLs;           // VCF PL likelihoods (rw)
    int nPLs, mPLs, nac;
    int32_t *gts, *ac;      // GTs and AC (w)
    double *pdg;            // PLs converted to P(D|G)
    float *anno16; int n16; // see anno[16] in bam2bcf.h
    double theta;           // prior
}
call_t;

void error(const char *format, ...);

/*
 *  call() - return -1 value on critical error; -2 to skip the site; or the number of non-reference
 *            alleles on success.
 */
int mcall(call_t *call, bcf1_t *rec);    // multiallic and rare-variant calling model
int ccall(call_t *call, bcf1_t *rec);    // the default consensus calling model
int qcall(call_t *call, bcf1_t *rec);    // QCall output

void mcall_init(call_t *call);
void ccall_init(call_t *call);
void qcall_init(call_t *call);

void mcall_destroy(call_t *call);
void ccall_destroy(call_t *call);
void qcall_destroy(call_t *call);

void call_init_pl2p(call_t *call);
uint32_t *call_trio_prep(int is_x, int is_son);

void init_allele_trimming_maps(call_t *call, int nals_ori, int als_out);
void mcall_trim_and_update_numberR(call_t *call, bcf1_t *rec, int nals_ori, int nals_new);

#endif
