/*  mcall.c -- multiallelic and rare variant calling.

    Copyright (C) 2012-2022 Genome Research Ltd.

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

#include <assert.h>
#include <math.h>
#include <inttypes.h>
#include <ctype.h>
#include <htslib/kfunc.h>
#include <htslib/khash_str2int.h>
#include "call.h"
#include "prob1.h"

// Using priors for GTs does not seem to be mathematically justified. Although
// it seems effective in removing false calls, it also flips a significant
// proportion of HET genotypes. Better is to filter by FORMAT/GQ using
// `bcftools filter`.
#define USE_PRIOR_FOR_GTS 0

// Go with uniform PLs for samples with no coverage. If unset, missing
// genotypes is reported instead.
#define FLAT_PDG_FOR_MISSING 0

int test16(float *anno16, anno16_t *a);

void qcall_init(call_t *call) { return; }
void qcall_destroy(call_t *call) { return; }
int qcall(call_t *call, bcf1_t *rec)
{
    // QCall format:
    //  chromosome, position, reference allele, depth, mapping quality, 0, ..
    error("TODO: qcall output\n");
    return 0;
}

void call_init_pl2p(call_t *call)
{
    int i;
    for (i=0; i<256; i++)
        call->pl2p[i] = pow(10., -i/10.);
}

// Macros for accessing call->trio and call->ntrio
#define FTYPE_222 0     // family type: all diploid
#define FTYPE_121 1     // chrX, the child is a boy
#define FTYPE_122 2     // chrX, a girl
#define FTYPE_101 3     // chrY, boy
#define FTYPE_100 4     // chrY, girl

#define GT_SKIP 0xf     // empty genotype (chrY in females)

#define IS_POW2(x) (!((x) & ((x) - 1)))    // zero is permitted
#define IS_HOM(x)  IS_POW2(x)

// Pkij = P(k|i,j) tells how likely it is to be a het if the parents
// are homs etc. The consistency of i,j,k has been already checked.
// Parameters are alleles and ploidy of father, mother, kid
// Returns 2/Pkij.
int calc_Pkij(int fals, int mals, int kals, int fpl, int mpl, int kpl)
{
    int als = fals|mals|kals;
    if ( IS_HOM(als) ) return 2;    // all are the same: child must be a HOM, P=1

    if ( fpl==1 )
    {
        if ( kpl==1 )   // chr X, the child is a boy, the copy is inherited from the mother
        {
            if ( IS_HOM(mals) ) return 2;   // 0 11 -> P(1) = 1
            return 4;                       // 0 01 -> P(0) = P(1) = 1/2
        }
        // chr X, the child is a girl
        if ( IS_HOM(mals) ) return 2;       // 0 11 -> P(01) = 1
        return 4;                           // 0 01 -> P(00) = P(01) = 1/2
    }

    if ( IS_HOM(fals) && IS_HOM(mals) ) return 2;   // 00 11 01, the child must be a HET, P=1
    if ( !IS_HOM(fals) && !IS_HOM(mals) )
    {
        if ( IS_HOM(kals) ) return 8;   // 01 01 00 or 01 01 11, P(k=HOM) = 1/4
        return 4;                       // 01 01 01, P(k=HET) = 1/2
    }
    return 4;   // 00 01, P(k=HET) = P(k=HOM) = 1/2
}

// Initialize ntrio and trio: ntrio lists the number of possible
// genotypes given combination of haploid/diploid genomes and the
// number of alleles. trio lists allowed genotype combinations:
//      4bit: 2/Pkij, 4: father, 4: mother, 4: child
// See also mcall_call_trio_genotypes()
//
static void mcall_init_trios(call_t *call)
{
    if ( call->prior_AN )
    {
        int id = bcf_hdr_id2int(call->hdr,BCF_DT_ID,call->prior_AN);
        if ( id==-1 ) error("No such tag \"%s\"\n", call->prior_AN);
        if ( !bcf_hdr_idinfo_exists(call->hdr,BCF_HL_FMT,id) )  error("No such FORMAT tag \"%s\"\n", call->prior_AN);
        id = bcf_hdr_id2int(call->hdr,BCF_DT_ID,call->prior_AC);
        if ( id==-1 ) error("No such tag \"%s\"\n", call->prior_AC);
        if ( !bcf_hdr_idinfo_exists(call->hdr,BCF_HL_FMT,id) )  error("No such FORMAT tag \"%s\"\n", call->prior_AC);
    }

    // 23, 138, 478 possible diploid trio genotypes with 2, 3, 4 alleles
    call->ntrio[FTYPE_222][2] = 15; call->ntrio[FTYPE_222][3] = 78;  call->ntrio[FTYPE_222][4] = 250;
    call->ntrio[FTYPE_121][2] = 8;  call->ntrio[FTYPE_121][3] = 27;  call->ntrio[FTYPE_121][4] = 64;
    call->ntrio[FTYPE_122][2] = 8;  call->ntrio[FTYPE_122][3] = 27;  call->ntrio[FTYPE_122][4] = 64;
    call->ntrio[FTYPE_101][2] = 2;  call->ntrio[FTYPE_101][3] = 3;   call->ntrio[FTYPE_101][4] = 4;
    call->ntrio[FTYPE_100][2] = 2;  call->ntrio[FTYPE_100][3] = 3;   call->ntrio[FTYPE_100][4] = 4;

    int nals, itype;
    for (itype=0; itype<=4; itype++)
    {
        for (nals=2; nals<=4; nals++)
            call->trio[itype][nals] = (uint16_t*) malloc(sizeof(uint16_t)*call->ntrio[itype][nals]);
    }

    // max 10 possible diploid genotypes
    int gts[10];
    for (nals=2; nals<=4; nals++)
    {
        int i,j,k, n = 0, ngts = 0;
        for (i=0; i<nals; i++)
            for (j=0; j<=i; j++)
                gts[ngts++] = 1<<i | 1<<j;

        // 222: all diploid
        // i,j,k: father, mother, child
        for (i=0; i<ngts; i++)
            for (j=0; j<ngts; j++)
                for (k=0; k<ngts; k++)
                {
                    if ( ((gts[i]|gts[j])&gts[k]) != gts[k] ) continue;             // k not present in neither i nor j
                    if ( !(gts[i] & gts[k]) || !(gts[j] & gts[k]) ) continue;       // one copy from father, one from mother
                    int Pkij = calc_Pkij(gts[i],gts[j],gts[k], 2,2,2);
                    call->trio[FTYPE_222][nals][n++] = Pkij<<12 | i<<8 | j<<4 | k;  // father, mother, child
                }
        assert( n==call->ntrio[FTYPE_222][nals] );

        // 121: chrX, boy
        n = 0;
        for (i=0; i<ngts; i++)
            for (j=0; j<ngts; j++)
                for (k=0; k<ngts; k++)
                {
                    if ( !IS_HOM(gts[i]) || !IS_HOM(gts[k]) ) continue;   // father nor boy can be diploid
                    if ( ((gts[i]|gts[j])&gts[k]) != gts[k] ) continue;
                    if ( !(gts[j] & gts[k]) ) continue;     // boy must inherit the copy from mother
                    int Pkij = calc_Pkij(gts[i],gts[j],gts[k], 1,2,1);
                    call->trio[FTYPE_121][nals][n++] = Pkij<<12 | i<<8 | j<<4 | k;
                }
        assert( n==call->ntrio[FTYPE_121][nals] );

        // 122: chrX, girl
        n = 0;
        for (i=0; i<ngts; i++)
            for (j=0; j<ngts; j++)
                for (k=0; k<ngts; k++)
                {
                    if ( !IS_HOM(gts[i]) ) continue;
                    if ( ((gts[i]|gts[j])&gts[k]) != gts[k] ) continue;
                    if ( !(gts[i] & gts[k]) ) continue;     // girl must inherit one copy from the father and one from the mother
                    if ( !(gts[j] & gts[k]) ) continue;
                    int Pkij = calc_Pkij(gts[i],gts[j],gts[k], 1,2,2);
                    call->trio[FTYPE_122][nals][n++] = Pkij<<12 | i<<8 | j<<4 | k;
                }
        assert( n==call->ntrio[FTYPE_122][nals] );

        // 101: chrY, boy
        n = 0;
        for (i=0; i<ngts; i++)
            for (k=0; k<ngts; k++)
            {
                if ( !IS_HOM(gts[i]) || !IS_HOM(gts[k]) ) continue;
                if ( (gts[i]&gts[k]) != gts[k] ) continue;
                call->trio[FTYPE_101][nals][n++] = 1<<12 | i<<8 | GT_SKIP<<4 | k;
            }
        assert( n==call->ntrio[FTYPE_101][nals] );

        // 100: chrY, girl
        n = 0;
        for (i=0; i<ngts; i++)
        {
            if ( !IS_POW2(gts[i]) ) continue;
            call->trio[FTYPE_100][nals][n++] = 1<<12 | i<<8 | GT_SKIP<<4 | GT_SKIP;
        }
        assert( n==call->ntrio[FTYPE_100][nals] );

    }
    call->GLs = (double*) calloc(bcf_hdr_nsamples(call->hdr)*10,sizeof(double));

    int i, j;
    for (i=0; i<call->nfams; i++)
    {
        family_t *fam = &call->fams[i];
        int ploidy[3];
        for (j=0; j<3; j++)
            ploidy[j] = call->ploidy[fam->sample[j]];

        if ( ploidy[FATHER]==2 )    // not X, not Y
        {
            if ( ploidy[MOTHER]!=2 || ploidy[CHILD]!=2 )
                error("Incorrect ploidy: %d %d %d\n", ploidy[FATHER],ploidy[MOTHER],ploidy[CHILD]);
            fam->type = FTYPE_222;
            continue;
        }
        if ( ploidy[FATHER]!=1 || ploidy[MOTHER]==1 )
                error("Incorrect ploidy: %d %d %d\n", ploidy[FATHER],ploidy[MOTHER],ploidy[CHILD]);
        if ( ploidy[MOTHER]==2 )    // X
        {
            if ( ploidy[CHILD]==0 )
                error("Incorrect ploidy: %d %d %d\n", ploidy[FATHER],ploidy[MOTHER],ploidy[CHILD]);
            fam->type = ploidy[CHILD]==2 ? FTYPE_122 : FTYPE_121;   // a girl or a boy
        }
        else    // Y
        {
            if ( ploidy[CHILD]==2 )
                error("Incorrect ploidy: %d %d %d\n", ploidy[FATHER],ploidy[MOTHER],ploidy[CHILD]);
            fam->type = ploidy[CHILD]==0 ? FTYPE_100 : FTYPE_101;   // a girl or a boy
        }
    }
}
static void mcall_destroy_trios(call_t *call)
{
    int i, j;
    for (i=2; i<=4; i++)
        for (j=0; j<=4; j++)
            free(call->trio[j][i]);
}

static void init_sample_groups(call_t *call)
{
    int i, nsmpl = bcf_hdr_nsamples(call->hdr);
    if ( !call->sample_groups )
    {
        // standard pooled calling, all samples in the same group
        call->nsmpl_grp = 1;
        call->smpl_grp  = (smpl_grp_t*)calloc(1,sizeof(*call->smpl_grp));
        call->smpl_grp[0].nsmpl = nsmpl;
        call->smpl_grp[0].smpl  = (uint32_t*)calloc(call->smpl_grp[0].nsmpl,sizeof(uint32_t));
        for (i=0; i<nsmpl; i++)
            call->smpl_grp[0].smpl[i] = i;
        return;
    }

    if ( call->sample_groups_tag )
    {
        // Is the tag defined in the header?
        int tag_id = bcf_hdr_id2int(call->hdr,BCF_DT_ID,call->sample_groups_tag);
        if ( tag_id==-1 ) error("No such tag \"%s\"\n",call->sample_groups_tag);
        if ( !bcf_hdr_idinfo_exists(call->hdr,BCF_HL_FMT,tag_id) )  error("No such FORMAT tag \"%s\"\n", call->sample_groups_tag);
    }
    else
    {
        int tag_id = bcf_hdr_id2int(call->hdr,BCF_DT_ID,"QS");
        if ( tag_id >= 0 && bcf_hdr_idinfo_exists(call->hdr,BCF_HL_FMT,tag_id) ) call->sample_groups_tag = "QS";
        else
        {
            tag_id = bcf_hdr_id2int(call->hdr,BCF_DT_ID,"AD");
            if ( tag_id >= 0 && bcf_hdr_idinfo_exists(call->hdr,BCF_HL_FMT,tag_id) ) call->sample_groups_tag = "AD";
            else error("Error: neither \"AD\" nor \"QS\" FORMAT tag exists and no alternative given with -G\n");
        }
    }

    // Read samples/groups
    if ( !strcmp("-",call->sample_groups) )
    {
        // single-sample calling, each sample creates its own group
        call->nsmpl_grp = nsmpl;
        call->smpl_grp  = (smpl_grp_t*)calloc(nsmpl,sizeof(*call->smpl_grp));
        for (i=0; i<nsmpl; i++)
        {
            call->smpl_grp[i].nsmpl = 1;
            call->smpl_grp[i].smpl  = (uint32_t*)calloc(call->smpl_grp[i].nsmpl,sizeof(uint32_t));
            call->smpl_grp[i].smpl[0] = i;
        }
    }
    else
    {
        int nlines;
        char **lines = hts_readlist(call->sample_groups, 1, &nlines);
        if ( !lines ) error("Could not read the file: %s\n", call->sample_groups);

        uint32_t *smpl2grp = (uint32_t*)calloc(nsmpl,sizeof(uint32_t));
        uint32_t *grp2n = (uint32_t*)calloc(nsmpl,sizeof(uint32_t));
        void *grp2idx = khash_str2int_init();

        call->nsmpl_grp = 0;
        for (i=0; i<nlines; i++)
        {
            char *ptr = lines[i];
            while ( *ptr && !isspace(*ptr) ) ptr++;
            if ( !*ptr ) error("Could not parse the line in %s, expected a sample name followed by tab and a population name: %s\n",call->sample_groups,lines[i]);
            char *tmp = ptr;
            while ( *ptr && isspace(*ptr) ) ptr++;
            if ( !*ptr ) error("Could not parse the line in %s, expected a sample name followed by tab and a population name: %s\n",call->sample_groups,lines[i]);
            *tmp = 0;
            int ismpl = bcf_hdr_id2int(call->hdr, BCF_DT_SAMPLE, lines[i]);
            if ( ismpl<0 ) continue;
            if ( smpl2grp[ismpl] ) error("Error: the sample \"%s\" is listed twice in %s\n", lines[i],call->sample_groups);
            if ( !khash_str2int_has_key(grp2idx,ptr+1) )
            {
                khash_str2int_set(grp2idx, ptr+1, call->nsmpl_grp);
                call->nsmpl_grp++;
            }
            int igrp = -1;
            if ( khash_str2int_get(grp2idx, ptr+1, &igrp)!=0 )
                error("This should not happen, fixme: %s\n",ptr+1);
            grp2n[igrp]++;
            smpl2grp[ismpl] = igrp+1;   // +1 to distinguish unlisted samples
        }
        khash_str2int_destroy(grp2idx);
        if ( !call->nsmpl_grp ) error("Could not parse the file, no matching samples found: %s\n", call->sample_groups);

        call->smpl_grp = (smpl_grp_t*)calloc(call->nsmpl_grp,sizeof(*call->smpl_grp));
        for (i=0; i<nsmpl; i++)
        {
            if ( !smpl2grp[i] ) error("Error: The sample \"%s\" is not listed in %s\n",call->hdr->samples[i],call->sample_groups);
            int igrp = smpl2grp[i] - 1;
            if ( !call->smpl_grp[igrp].nsmpl )
                call->smpl_grp[igrp].smpl = (uint32_t*)calloc(grp2n[igrp],sizeof(uint32_t));
            call->smpl_grp[igrp].smpl[call->smpl_grp[igrp].nsmpl] = i;
            call->smpl_grp[igrp].nsmpl++;
        }
        free(smpl2grp);
        free(grp2n);
        for (i=0; i<nlines; i++) free(lines[i]);
        free(lines);
    }
}
static void destroy_sample_groups(call_t *call)
{
    int i;
    for (i=0; i<call->nsmpl_grp; i++)
    {
        free(call->smpl_grp[i].qsum);
        free(call->smpl_grp[i].smpl);
    }
    free(call->smpl_grp);
}

void mcall_init(call_t *call)
{
    init_sample_groups(call);
    call_init_pl2p(call);

    call->nals_map = 5;
    call->als_map  = (int*) malloc(sizeof(int)*call->nals_map);
    call->npl_map  = 5*(5+1)/2;     // will be expanded later if necessary
    call->pl_map   = (int*) malloc(sizeof(int)*call->npl_map);
    call->gts  = (int32_t*) calloc(bcf_hdr_nsamples(call->hdr)*2,sizeof(int32_t));   // assuming at most diploid everywhere

    if ( call->flag & CALL_CONSTR_TRIO )
    {
        call->cgts = (int32_t*) calloc(bcf_hdr_nsamples(call->hdr),sizeof(int32_t));
        call->ugts = (int32_t*) calloc(bcf_hdr_nsamples(call->hdr),sizeof(int32_t));
        mcall_init_trios(call);
        bcf_hdr_append(call->hdr,"##FORMAT=<ID=CGT,Number=1,Type=Integer,Description=\"Constrained Genotype (0-based index to Number=G ordering).\">");
        bcf_hdr_append(call->hdr,"##FORMAT=<ID=UGT,Number=1,Type=Integer,Description=\"Unconstrained Genotype (0-based index to Number=G ordering).\">");
    }
    if ( call->flag & CALL_CONSTR_ALLELES ) call->vcmp = vcmp_init();

    bcf_hdr_append(call->hdr,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    if ( call->output_tags & CALL_FMT_GQ )
        bcf_hdr_append(call->hdr,"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Phred-scaled Genotype Quality\">");
    if ( call->output_tags & CALL_FMT_GP )
        bcf_hdr_append(call->hdr,"##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype posterior probabilities in the range 0 to 1\">");
    if ( call->output_tags & (CALL_FMT_GQ|CALL_FMT_GP) )
        call->GQs = (int32_t*) malloc(sizeof(int32_t)*bcf_hdr_nsamples(call->hdr));
    bcf_hdr_append(call->hdr,"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes for each ALT allele, in the same order as listed\">");
    bcf_hdr_append(call->hdr,"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    bcf_hdr_append(call->hdr,"##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases\">");
    bcf_hdr_append(call->hdr,"##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Average mapping quality\">");
    if ( call->output_tags & CALL_FMT_PV4 )
        bcf_hdr_append(call->hdr,"##INFO=<ID=PV4,Number=4,Type=Float,Description=\"P-values for strand bias, baseQ bias, mapQ bias and tail distance bias\">\n");

    // init the prior
    if ( call->theta>0 )
    {
        int i, n = 0;
        if ( !call->ploidy ) n = 2*bcf_hdr_nsamples(call->hdr); // all are diploid
        else
        {
            for (i=0; i<bcf_hdr_nsamples(call->hdr); i++)
                n += call->ploidy[i];
        }
        // Watterson factor, here aM_1 = aM_2 = 1
        double aM = 1;
        for (i=2; i<n; i++) aM += 1./i;
        call->theta *= aM;
        if ( call->theta >= 1 )
        {
            fprintf(stderr,"The prior is too big (theta*aM=%.2f), going with 0.99\n", call->theta);
            call->theta = 0.99;
        }
        call->theta = log(call->theta);
    }
}

void mcall_destroy(call_t *call)
{
    destroy_sample_groups(call);
    if (call->vcmp) vcmp_destroy(call->vcmp);
    free(call->itmp);
    mcall_destroy_trios(call);
    free(call->GPs);
    free(call->ADs);
    free(call->GLs);
    free(call->GQs);
    free(call->anno16);
    free(call->PLs);
    free(call->als_map);
    free(call->pl_map);
    free(call->gts); free(call->cgts); free(call->ugts);
    free(call->pdg);
    free(call->als);
    free(call->ac);
    return;
}


// Inits P(D|G): convert PLs from log space and normalize. In case of zero
// depth, missing PLs are all zero. In this case, pdg's are set to 0
// so that the corresponding genotypes can be set as missing and the
// qual calculation is not affected.
// Missing values are replaced by generic likelihoods when X (unseen allele) is
// present.
// NB: While the -m callig model uses the pdgs in canonical order,
// the original samtools -c calling code uses pdgs in reverse order (AA comes
// first, RR last).
// NB: Ploidy is not taken into account here, which is incorrect.
void set_pdg(double *pl2p, int *PLs, double *pdg, int n_smpl, int n_gt, int unseen)
{
    int i, j, nals;

    // find out the number of alleles, expecting diploid genotype likelihoods
    bcf_gt2alleles(n_gt-1, &i, &nals);
    assert( i==nals );
    nals++;

    for (i=0; i<n_smpl; i++)
    {
        double sum = 0;
        for (j=0; j<n_gt; j++)
        {
            if ( PLs[j]==bcf_int32_vector_end )
            {
                // We expect diploid genotype likelihoods. If not diploid, treat as missing
                j = 0;
                break;
            }
            if ( PLs[j]==bcf_int32_missing ) break;
            pdg[j] = PLs[j] < 256 ? pl2p[PLs[j]] : pow(10., -PLs[j]/10.);
            sum += pdg[j];
        }

        if ( j==0 )
        {
            // First value is missing (LK of RR), this indicates that
            // all values are missing.
            j = sum = n_gt;
        }
        else if ( j<n_gt && unseen<0 )
        {
            // Some of the values are missing and the unseen allele LK is not
            // available. In such a case, we set LK to a very small value.
            sum = 0;
            for (j=0; j<n_gt; j++)
            {
                assert( PLs[j]!=bcf_int32_vector_end );
                if ( PLs[j]==bcf_int32_missing ) PLs[j] = 255;
                pdg[j] = PLs[j] < 256 ? pl2p[PLs[j]] : pow(10., -PLs[j]/10.);
                sum += pdg[j];
            }
        }
        if ( j<n_gt )
        {
            // Missing values present, fill with unseen allele LK. This can be only
            // as good as the merge was.
            int ia,ib, k;
            j = 0;
            sum = 0;
            for (ia=0; ia<nals; ia++)
            {
                for (ib=0; ib<=ia; ib++)
                {
                    if ( PLs[j]==bcf_int32_missing )
                    {
                        k = bcf_alleles2gt(ia,unseen);
                        if ( PLs[k]==bcf_int32_missing ) k = bcf_alleles2gt(ib,unseen);
                        if ( PLs[k]==bcf_int32_missing ) k = bcf_alleles2gt(unseen,unseen);
                        if ( PLs[k]==bcf_int32_missing )
                        {
                            // The PLs for unseen allele X are not present as well as for ia, ib.
                            // This can happen with incremental calling, when one of the merged
                            // files had all alleles A,C,G,T, in such a case, X was not present.
                            // Use a very small value instead.
                            PLs[j] = 255;
                        }
                        else
                            PLs[j] = PLs[k];
                    }
                    pdg[j] = pl2p[ PLs[j] ];
                    sum += pdg[j];
                    j++;
                }
            }
        }
        // Normalize: sum_i pdg_i = 1
        if ( sum==n_gt )
        {
            // all missing
            #if FLAT_PDG_FOR_MISSING
                for (j=0; j<n_gt; j++) pdg[j] = 1./n_gt;
            #else
                for (j=0; j<n_gt; j++) pdg[j] = 0;
            #endif
        }
        else
            for (j=0; j<n_gt; j++) pdg[j] /= sum;

        PLs += n_gt;
        pdg += n_gt;
    }
}

// Create mapping between old and new (trimmed) alleles
void init_allele_trimming_maps(call_t *call, int nals_ori, int als_out)
{
    int i, j, nout = 0;

    // als_map: old(i) -> new(j)
    for (i=0; i<nals_ori; i++)
    {
        if ( als_out & (1<<i) ) call->als_map[i] = nout++;
        else call->als_map[i] = -1;
    }

    if ( !call->pl_map ) return;

    // pl_map: new(k) -> old(l)
    int k = 0, l = 0;
    for (i=0; i<nals_ori; i++)
    {
        for (j=0; j<=i; j++)
        {
            if ( (als_out & (1<<i)) && (als_out & (1<<j)) ) call->pl_map[k++] = l;
            l++;
        }
    }
}

/** log(exp(a)+exp(b)) */
static inline double logsumexp2(double a, double b)
{
    if ( a>b )
        return log(1 + exp(b-a)) + a;
    else
        return log(1 + exp(a-b)) + b;
}

// Macro to set the most likely alleles
#define UPDATE_MAX_LKs(als,sum) { \
     if ( max_lk<lk_tot && lk_tot_set ) { max_lk = lk_tot; max_als = (als); } \
     if ( sum ) lk_sum = logsumexp2(lk_tot,lk_sum); \
}

#define SWAP(type_t,x,y) {type_t tmp; tmp = x; x = y; y = tmp; }

// Determine the most likely combination of alleles. In this implementation,
// at most tri-allelic sites are considered. Returns the number of alleles.
static int mcall_find_best_alleles(call_t *call, int nals, smpl_grp_t *grp)
{
    int ia,ib,ic;   // iterators over up to three alleles
    int max_als=0;  // most likely combination of alleles
    double ref_lk = -HUGE_VAL, max_lk = -HUGE_VAL; // likelihood of the reference and of most likely combination of alleles
    double lk_sum = -HUGE_VAL;    // for normalizing the likelihoods
    int nsmpl = grp->nsmpl;
    int ngts  = nals*(nals+1)/2;

    // Single allele
    for (ia=0; ia<nals; ia++)
    {
        double lk_tot  = 0;
        int lk_tot_set = 0;
        int iaa = (ia+1)*(ia+2)/2-1;    // index in PL which corresponds to the homozygous "ia/ia" genotype
        int ismpl;
        for (ismpl=0; ismpl<nsmpl; ismpl++)
        {
            double *pdg = call->pdg + grp->smpl[ismpl]*ngts + iaa;
            if ( *pdg ) { lk_tot += log(*pdg); lk_tot_set = 1; }
        }
        if ( ia==0 ) ref_lk = lk_tot;   // likelihood of 0/0 for all samples
        else lk_tot += call->theta; // the prior
        UPDATE_MAX_LKs(1<<ia, ia>0 && lk_tot_set);
    }

    // Two alleles
    if ( nals>1 )
    {
        for (ia=0; ia<nals; ia++)
        {
            if ( grp->qsum[ia]==0 ) continue;
            int iaa = (ia+1)*(ia+2)/2-1;
            for (ib=0; ib<ia; ib++)
            {
                if ( grp->qsum[ib]==0 ) continue;
                double lk_tot  = 0;
                int lk_tot_set = 0;
                double fa  = grp->qsum[ia]/(grp->qsum[ia] + grp->qsum[ib]);
                double fb  = grp->qsum[ib]/(grp->qsum[ia] + grp->qsum[ib]);
                double fa2 = fa*fa;
                double fb2 = fb*fb;
                double fab = 2*fa*fb;
                int is, ibb = (ib+1)*(ib+2)/2-1, iab = iaa - ia + ib;
                for (is=0; is<nsmpl; is++)
                {
                    int ismpl = grp->smpl[is];
                    double *pdg = call->pdg + ismpl*ngts;
                    double val = 0;
                    if ( !call->ploidy || call->ploidy[ismpl]==2 )
                        val = fa2*pdg[iaa] + fb2*pdg[ibb] + fab*pdg[iab];
                    else if ( call->ploidy && call->ploidy[ismpl]==1 )
                        val = fa*pdg[iaa] + fb*pdg[ibb];
                    if ( val ) { lk_tot += log(val); lk_tot_set = 1; }
                }
                if ( ia!=0 ) lk_tot += call->theta;    // the prior
                if ( ib!=0 ) lk_tot += call->theta;
                UPDATE_MAX_LKs(1<<ia|1<<ib, lk_tot_set);
            }
        }
    }

    // Three alleles
    if ( nals>2 )
    {
        for (ia=0; ia<nals; ia++)
        {
            if ( grp->qsum[ia]==0 ) continue;
            int iaa = (ia+1)*(ia+2)/2-1;
            for (ib=0; ib<ia; ib++)
            {
                if ( grp->qsum[ib]==0 ) continue;
                int ibb = (ib+1)*(ib+2)/2-1;
                int iab = iaa - ia + ib;
                for (ic=0; ic<ib; ic++)
                {
                    if ( grp->qsum[ic]==0 ) continue;
                    double lk_tot  = 0;
                    int lk_tot_set = 0;

                    double fa  = grp->qsum[ia]/(grp->qsum[ia] + grp->qsum[ib] + grp->qsum[ic]);
                    double fb  = grp->qsum[ib]/(grp->qsum[ia] + grp->qsum[ib] + grp->qsum[ic]);
                    double fc  = grp->qsum[ic]/(grp->qsum[ia] + grp->qsum[ib] + grp->qsum[ic]);
                    double fa2 = fa*fa;
                    double fb2 = fb*fb;
                    double fc2 = fc*fc;
                    double fab = 2*fa*fb, fac = 2*fa*fc, fbc = 2*fb*fc;
                    int is, icc = (ic+1)*(ic+2)/2-1;
                    int iac = iaa - ia + ic, ibc = ibb - ib + ic;
                    for (is=0; is<nsmpl; is++)
                    {
                        int ismpl = grp->smpl[is];
                        double *pdg = call->pdg + ismpl*ngts;
                        double val = 0;
                        if ( !call->ploidy || call->ploidy[ismpl]==2 )
                            val = fa2*pdg[iaa] + fb2*pdg[ibb] + fc2*pdg[icc] + fab*pdg[iab] + fac*pdg[iac] + fbc*pdg[ibc];
                        else if ( call->ploidy && call->ploidy[ismpl]==1 )
                            val = fa*pdg[iaa] + fb*pdg[ibb] + fc*pdg[icc];
                        if ( val ) { lk_tot += log(val); lk_tot_set = 1; }
                    }
                    if ( ia!=0 ) lk_tot += call->theta;    // the prior
                    if ( ib!=0 ) lk_tot += call->theta;    // the prior
                    if ( ic!=0 ) lk_tot += call->theta;    // the prior
                    UPDATE_MAX_LKs(1<<ia|1<<ib|1<<ic, lk_tot_set);
                }
            }
        }
    }

    int i, n = 0;
    for (i=0; i<nals; i++) if ( max_als & 1<<i) n++;

    grp->max_lk = max_lk;
    grp->ref_lk = ref_lk;
    grp->lk_sum = lk_sum;
    grp->als  = max_als;
    grp->nals = n;

    return n;
}

// Sets GT=0/0 or GT=. if PL=0,0,0
static void mcall_set_ref_genotypes(call_t *call, int nals_ori)
{
    int i;
    int ngts  = nals_ori*(nals_ori+1)/2;            // need this to distinguish between GT=0/0 vs GT=.
    int nsmpl = bcf_hdr_nsamples(call->hdr);

    for (i=0; i<nals_ori; i++) call->ac[i] = 0;     // nals_new<=nals_ori, never mind setting extra 0's

    // Set all genotypes to 0/0 or 0
    int *gts    = call->gts;
    double *pdg = call->pdg;
    int isample;
    for (isample = 0; isample < nsmpl; isample++)
    {
        int ploidy = call->ploidy ? call->ploidy[isample] : 2;
        for (i=0; i<ngts; i++) if ( pdg[i]!=0.0 ) break;
        if ( i==ngts || !ploidy )
        {
            gts[0] = bcf_gt_missing;
            gts[1] = ploidy==2 ? bcf_gt_missing : bcf_int32_vector_end;
        }
        else
        {
            gts[0] = bcf_gt_unphased(0);
            gts[1] = ploidy==2 ? bcf_gt_unphased(0) : bcf_int32_vector_end;
            call->ac[0] += ploidy;
        }
        gts += 2;
        pdg += ngts;
    }
}

static void mcall_call_genotypes(call_t *call, int nals_ori, smpl_grp_t *grp)
{
    int ia, ib, i;
    int ngts_ori = nals_ori*(nals_ori+1)/2;
    int ngts_new = call->nals_new*(call->nals_new+1)/2;
    int nsmpl = grp->nsmpl;

    #if USE_PRIOR_FOR_GTS
        float prior = exp(call->theta);
    #endif

    int is;
    for (is = 0; is < nsmpl; is++)
    {
        int ismpl   = grp->smpl[is];
        double *pdg = call->pdg + ismpl*ngts_ori;
        float *gps  = call->GPs + ismpl*ngts_new;
        int *gts    = call->gts + ismpl*2;

        int ploidy = call->ploidy ? call->ploidy[ismpl] : 2;
        assert( ploidy>=0 && ploidy<=2 );

        if ( !ploidy )
        {
            gts[0] = bcf_gt_missing;
            gts[1] = bcf_int32_vector_end;
            gps[0] = -1;
            continue;
        }

        #if !FLAT_PDG_FOR_MISSING
            // Skip samples with zero depth, they have all pdg's equal to 0
            for (i=0; i<ngts_ori; i++) if ( pdg[i]!=0.0 ) break;
            if ( i==ngts_ori )
            {
                gts[0] = bcf_gt_missing;
                gts[1] = ploidy==2 ? bcf_gt_missing : bcf_int32_vector_end;
                gps[0] = -1;
                continue;
            }
        #endif

        // Default fallback for the case all LKs are the same
        gts[0] = bcf_gt_unphased(0);
        gts[1] = ploidy==2 ? bcf_gt_unphased(0) : bcf_int32_vector_end;

        // Non-zero depth, determine the most likely genotype
        double best_lk = 0;
        for (ia=0; ia<nals_ori; ia++)
        {
            if ( !(grp->als & 1<<ia) ) continue;    // ia-th allele not in the final selection, skip
            int iaa = (ia+1)*(ia+2)/2-1;                // PL index of the ia/ia genotype
            double lk = ploidy==2 ? pdg[iaa]*grp->qsum[ia]*grp->qsum[ia] : pdg[iaa]*grp->qsum[ia];
            #if USE_PRIOR_FOR_GTS
                if ( ia!=0 ) lk *= prior;
            #endif
            int igt  = ploidy==2 ? bcf_alleles2gt(call->als_map[ia],call->als_map[ia]) : call->als_map[ia];
            gps[igt] = lk;
            if ( best_lk < lk )
            {
                best_lk = lk;
                gts[0] = bcf_gt_unphased(call->als_map[ia]);
            }
        }
        if ( ploidy==2 )
        {
            gts[1] = gts[0];
            for (ia=0; ia<nals_ori; ia++)
            {
                if ( !(grp->als & 1<<ia) ) continue;
                int iaa = (ia+1)*(ia+2)/2-1;
                for (ib=0; ib<ia; ib++)
                {
                    if ( !(grp->als & 1<<ib) ) continue;
                    int iab = iaa - ia + ib;
                    double lk = 2*pdg[iab]*grp->qsum[ia]*grp->qsum[ib];
                    #if USE_PRIOR_FOR_GTS
                        if ( ia!=0 ) lk *= prior;
                        if ( ib!=0 ) lk *= prior;
                    #endif
                    int igt  = bcf_alleles2gt(call->als_map[ia],call->als_map[ib]);
                    gps[igt] = lk;
                    if ( best_lk < lk )
                    {
                        best_lk = lk;
                        gts[0] = bcf_gt_unphased(call->als_map[ib]);
                        gts[1] = bcf_gt_unphased(call->als_map[ia]);
                    }
                }
            }
        }
        else
            gts[1] = bcf_int32_vector_end;

        call->ac[ bcf_gt_allele(gts[0]) ]++;
        if ( gts[1]!=bcf_int32_vector_end ) call->ac[ bcf_gt_allele(gts[1]) ]++;
    }
    if ( !(call->output_tags & (CALL_FMT_GQ|CALL_FMT_GP)) ) return;
    double max, sum;
    for (is=0; is<nsmpl; is++)
    {
        int ismpl  = grp->smpl[is];
        float *gps = call->GPs + ismpl*ngts_new;

        int nmax;
        if ( call->ploidy )
        {
            if ( call->ploidy[ismpl]==2 ) nmax = ngts_new;
            else if ( call->ploidy[ismpl]==1 ) nmax = grp->nals;
            else nmax = 0;
        }
        else nmax = ngts_new;

        max = gps[0];
        if ( max<0 || nmax==0 )
        {
            // no call
            if ( call->output_tags & CALL_FMT_GP )
            {
                for (i=0; i<nmax; i++) gps[i] = 0;
                if ( nmax==0 ) { bcf_float_set_missing(gps[i]); nmax++; }
                if ( nmax < ngts_new ) bcf_float_set_vector_end(gps[nmax]);
            }
            call->GQs[ismpl] = 0;
            continue;
        }
        sum = gps[0];
        for (i=1; i<nmax; i++)
        {
            if ( max < gps[i] ) max = gps[i];
            sum += gps[i];
        }
        max = -4.34294*log(1 - max/sum);
        call->GQs[ismpl] = max<=INT8_MAX ? max : INT8_MAX;
        if ( call->output_tags & CALL_FMT_GP )
        {
            assert( max );
            for (i=0; i<nmax; i++) gps[i] = gps[i]/sum;
            for (; i<ngts_new; i++) bcf_float_set_vector_end(gps[i]);
        }
    }
}


/**
    Pm = P(mendelian) .. parameter to vary, 1-Pm is the probability of novel mutation.
                         When trio_Pm_ins is negative, Pm is calculated dynamically
                         according to indel length. For simplicity, only the
                         first ALT is considered.
    Pkij = P(k|i,j)   .. probability that the genotype combination i,j,k is consistent
                         with mendelian inheritance (the likelihood that offspring
                         of two HETs is a HOM is smaller than it being a HET)

    P_uc(F=i,M=j,K=k) = P(F=i) . P(M=j) . P(K=k)  .. unconstrained P
    P_c(F=i,M=j,K=k) = P_uc . Pkij                .. constrained P
    P(F=i,M=j,K=k) = P_uc . (1 - Pm) + P_c . Pm
                   = P_uc . [1 - Pm + Pkij . Pm]

    We choose genotype combination i,j,k which maximizes P(F=i,M=j,K=k). This
    probability gives the quality GQ(Trio).
    Individual qualities are calculated as
        GQ(F=i,M=j,K=k) = P(F=i,M=j,K=k) / \sum_{x,y} P(F=i,M=x,K=y)
 */
#if 0
static void mcall_call_trio_genotypes(call_t *call, bcf1_t *rec, int nals, int nals_new, int als_new)
{
    int ia, ib, i;
    int nsmpl    = bcf_hdr_nsamples(call->hdr);
    int ngts     = nals*(nals+1)/2;
    int nout_gts = nals_new*(nals_new+1)/2;
    double *gls  = call->GLs - nout_gts;
    double *pdg  = call->pdg - ngts;

    // Calculate individuals' genotype likelihoods P(X=i)
    int isample;
    for (isample = 0; isample < nsmpl; isample++)
    {
        int ploidy = call->ploidy ? call->ploidy[isample] : 2;
        int32_t *gts = call->ugts + isample;

        gls += nout_gts;
        pdg += ngts;

        // Skip samples with all pdg's equal to 1. These have zero depth.
        for (i=0; i<ngts; i++) if ( pdg[i]!=0.0 ) break;
        if ( i==ngts || !ploidy )
        {
            gts[0] = -1;
            gls[0] = 1;
            continue;
        }

        for (i=0; i<nout_gts; i++) gls[i] = -HUGE_VAL;

        grp1_t *grp = &call->smpl_grp.grp[call->smpl_grp.smpl2grp[isample]];
        double sum_lk  = 0;
        double best_lk = 0;
        for (ia=0; ia<nals; ia++)
        {
            if ( !(als_new & 1<<ia) ) continue;     // ia-th allele not in the final selection, skip
            int iaa   = bcf_alleles2gt(ia,ia);      // PL index of the ia/ia genotype
            int idx   = bcf_alleles2gt(call->als_map[ia],call->als_map[ia]);
            double lk = ploidy==2 ? pdg[iaa]*grp->qsum[ia]*grp->qsum[ia] : pdg[iaa]*grp->qsum[ia];
            sum_lk   += lk;
            gls[idx]  = lk;
            if ( best_lk < lk )
            {
                best_lk = lk;
                gts[0] = bcf_alleles2gt(call->als_map[ia],call->als_map[ia]);
            }
        }
        if ( ploidy==2 )
        {
            for (ia=0; ia<nals; ia++)
            {
                if ( !(als_new & 1<<ia) ) continue;
                for (ib=0; ib<ia; ib++)
                {
                    if ( !(als_new & 1<<ib) ) continue;
                    int iab   = bcf_alleles2gt(ia,ib);
                    int idx   = bcf_alleles2gt(call->als_map[ia],call->als_map[ib]);
                    double lk = 2*pdg[iab]*grp->qsum[ia]*grp->qsum[ib];
                    sum_lk   += lk;
                    gls[idx]  = lk;
                    if ( best_lk < lk )
                    {
                        best_lk = lk;
                        gts[0] = bcf_alleles2gt(call->als_map[ib],call->als_map[ia]);
                    }
                }
            }
        }
        for (i=0; i<nout_gts; i++)
            if ( gls[i]!=-HUGE_VAL ) gls[i] = log(gls[i]/sum_lk);
    }

    // Set novel mutation rate for this site: using first ALT allele for simplicity.
    double trio_Pm;
    if ( call->trio_Pm_ins<0 && call->trio_Pm_del<0 ) trio_Pm = call->trio_Pm_SNPs;     // the same Pm for indels and SNPs requested
    else
    {
        int ret = bcf_get_variant_types(rec);
        if ( !(ret & VCF_INDEL) ) trio_Pm = call->trio_Pm_SNPs;
        else
        {
            if ( call->trio_Pm_ins<0 )  // dynamic calculation, trio_Pm_del holds the scaling factor
            {
                trio_Pm = rec->d.var[1].n<0 ? -21.9313 - 0.2856*rec->d.var[1].n : -22.8689 + 0.2994*rec->d.var[1].n;
                trio_Pm = 1 - call->trio_Pm_del * exp(trio_Pm);
            }
            else                        // snps and indels set explicitly
            {
                trio_Pm = rec->d.var[1].n<0 ? call->trio_Pm_del : call->trio_Pm_ins;
            }
        }
    }

    // Calculate constrained likelihoods and determine genotypes
    int ifm;
    for (ifm=0; ifm<call->nfams; ifm++)
    {
        family_t *fam = &call->fams[ifm];
        int ntrio = call->ntrio[fam->type][nals_new];
        uint16_t *trio = call->trio[fam->type][nals_new];

        // Unconstrained likelihood
        int uc_itr = 0;
        double uc_lk = 0;
        for (i=0; i<3; i++)     // for father, mother, child
        {
            int ismpl = fam->sample[i];
            double *gl = call->GLs + nout_gts*ismpl;
            if ( gl[0]==1 ) continue;
            int j, jmax = 0;
            double max  = gl[0];
            for (j=1; j<nout_gts; j++)
                if ( max < gl[j] ) { max = gl[j]; jmax = j; }
            uc_lk += max;
            uc_itr |= jmax << ((2-i)*4);
        }

        // Best constrained likelihood
        int c_itr = -1, itr, uc_is_mendelian = 0;
        double c_lk = -HUGE_VAL;
        for (itr=0; itr<ntrio; itr++)   // for each trio genotype combination
        {
            double lk = 0;
            int npresent = 0;
            for (i=0; i<3; i++)     // for father, mother, child
            {
                int ismpl = fam->sample[i];
                double *gl = call->GLs + nout_gts*ismpl;
                if ( gl[0]==1 ) continue;
                int igt = trio[itr]>>((2-i)*4) & 0xf;
                assert( !call->ploidy || call->ploidy[ismpl]>0 );
                if ( igt==GT_SKIP ) continue;
                lk += gl[igt];
                npresent++;
                // fprintf(stderr," %e", gl[igt]);
            }
            // fprintf(stderr,"\t\t");
            double Pkij = npresent==3 ? (double)2/(trio[itr]>>12) : 1;  // with missing genotypes Pkij's are different
            lk += log(1 - trio_Pm * (1 - Pkij));
            // fprintf(stderr,"%d%d%d\t%e\t%.2f\n", trio[itr]>>8&0xf,trio[itr]>>4&0xf,trio[itr]&0xf, lk, Pkij);
            if ( c_lk < lk ) { c_lk = lk; c_itr = trio[itr]; }
            if ( uc_itr==trio[itr] ) uc_is_mendelian = 1;
        }

        if ( !uc_is_mendelian )
        {
            uc_lk += log(1 - trio_Pm);
            // fprintf(stderr,"c_lk=%e uc_lk=%e c_itr=%d%d%d uc_itr=%d%d%d\n", c_lk,uc_lk,c_itr>>8&0xf,c_itr>>4&0xf,c_itr&0xf,uc_itr>>8&0xf,uc_itr>>4&0xf,uc_itr&0xf);
            if ( c_lk < uc_lk ) { c_lk = uc_lk; c_itr = uc_itr; }
        }
        // fprintf(stderr,"best_lk=%e best_itr=%d%d%d uc_itr=%d%d%d\n", c_lk,c_itr>>8&0xf,c_itr>>4&0xf,c_itr&0xf,uc_itr>>8&0xf,uc_itr>>4&0xf,uc_itr&0xf);

        // Set genotypes for father, mother, child and calculate genotype qualities
        for (i=0; i<3; i++)
        {
            // GT
            int ismpl    = fam->sample[i];
            int igt      = c_itr>>((2-i)*4) & 0xf;
            double *gl   = call->GLs + nout_gts*ismpl;
            int32_t *gts = call->cgts + ismpl;
            if ( gl[0]==1 || igt==GT_SKIP )    // zero depth, set missing genotypes
            {
                gts[0] = -1;
                // bcf_float_set_missing(call->GQs[ismpl]);
                continue;
            }
            gts[0] = igt;

            #if 0
                // todo: Genotype Qualities
                //
                // GQ: for each family member i sum over all genotypes j,k keeping igt fixed
                double lk_sum = 0;
                for (itr=0; itr<ntrio; itr++)
                {
                    if ( igt != (trio[itr]>>((2-i)*4) & 0xf) ) continue;
                    double lk = 0;
                    int j;
                    for (j=0; j<3; j++)
                    {
                        int jsmpl = fam->sample[j];
                        double *gl = call->GLs + ngts*jsmpl;
                        if ( gl[0]==1 ) continue;
                        int jgt = trio[itr]>>((2-j)*4) & 0xf;
                        if ( jgt==GT_SKIP ) continue;
                        lk += gl[jgt];
                    }
                    double Pkij = (double)2/(trio[itr]>>12);
                    lk += log(1 - trio_Pm * (1 - Pkij));
                    lk_sum = logsumexp2(lk_sum, lk);
                }
                if ( !uc_is_mendelian && (best_itr>>((2-i)*4)&0xf)==(uc_itr>>((2-i)*4)&0xf) ) lk_sum = logsumexp2(lk_sum,uc_lk);
                call->GQs[ismpl] = -4.3429*(best_lk - lk_sum);
            #endif
        }
    }

    for (i=0; i<4; i++) call->ac[i] = 0;
    call->nhets = 0;
    call->ndiploid = 0;

    // Test if CGT,UGT are needed
    int ucgts_needed = 0;
    int32_t *cgts = call->cgts - 1;
    int32_t *ugts = call->ugts - 1;
    int32_t *gts  = call->gts - 2;
    for (isample = 0; isample < nsmpl; isample++)
    {
        int ploidy = call->ploidy ? call->ploidy[isample] : 2;
        cgts++;
        ugts++;
        gts += 2;
        if ( ugts[0]==-1 )
        {
            gts[0] = bcf_gt_missing;
            gts[1] = ploidy==2 ? bcf_gt_missing : bcf_int32_vector_end;
            continue;
        }
        int a,b;
        if ( cgts[0]!=ugts[0] )
        {
            bcf_gt2alleles(cgts[0], &a, &b);
            gts[0] = bcf_gt_unphased(a);
            gts[1] = ploidy==1 ? bcf_int32_vector_end : bcf_gt_unphased(b);
        }
        else
        {
            bcf_gt2alleles(ugts[0], &a, &b);
            gts[0] = bcf_gt_unphased(a);
            gts[1] = ploidy==1 ? bcf_int32_vector_end : bcf_gt_unphased(b);
        }
        if ( cgts[0]!=ugts[0] ) ucgts_needed = 1;
        call->ac[a]++;
        if ( ploidy==2 )
        {
            call->ac[b]++;
            call->ndiploid++;
            if ( a!=b ) call->nhets++;
        }
    }
    if ( ucgts_needed )
    {
        // Some GTs are different
        bcf_update_format_int32(call->hdr,rec,"UGT",call->ugts,nsmpl);
        bcf_update_format_int32(call->hdr,rec,"CGT",call->cgts,nsmpl);
    }
}
#endif

static void mcall_trim_and_update_PLs(call_t *call, bcf1_t *rec, int nals_ori, int nals_new)
{
    int npls_src = nals_ori*(nals_ori+1)/2;
    int npls_dst = nals_new*(nals_new+1)/2;     // number of PL values in diploid samples, ori and new
    if ( call->all_diploid && npls_src == npls_dst ) return;

    int *pls_src = call->PLs, *pls_dst = call->PLs;

    int nsmpl = bcf_hdr_nsamples(call->hdr);
    int isample, ia;
    for (isample = 0; isample < nsmpl; isample++)
    {
        int ploidy = call->ploidy ? call->ploidy[isample] : 2;
        if ( ploidy==2 )
        {
            for (ia=0; ia<npls_dst; ia++)
                pls_dst[ia] =  pls_src[ call->pl_map[ia] ];
        }
        else if ( ploidy==1 )
        {
            for (ia=0; ia<nals_new; ia++)
            {
                int isrc = (ia+1)*(ia+2)/2-1;
                pls_dst[ia] = pls_src[ call->pl_map[isrc] ];
            }
            if ( ia<npls_dst ) pls_dst[ia] = bcf_int32_vector_end;
        }
        else
        {
            pls_dst[0] = bcf_int32_missing;
            pls_dst[1] = bcf_int32_vector_end;  // relying on nals_new>1 in mcall()
        }
        pls_src += npls_src;
        pls_dst += npls_dst;
    }
    bcf_update_format_int32(call->hdr, rec, "PL", call->PLs, npls_dst*nsmpl);
}

void mcall_trim_and_update_numberR(call_t *call, bcf1_t *rec, int nals_ori, int nals_new)
{
    if ( nals_ori==nals_new ) return;

    int i,j, nret, size = sizeof(float);

    void *tmp_ori = call->itmp, *tmp_new = call->PLs;  // reusing PLs storage which is not used at this point
    int ntmp_ori = call->n_itmp, ntmp_new = call->mPLs;

    // INFO fields
    for (i=0; i<rec->n_info; i++)
    {
        bcf_info_t *info = &rec->d.info[i];
        int vlen = bcf_hdr_id2length(call->hdr,BCF_HL_INFO,info->key);
        if ( vlen!=BCF_VL_R ) continue; // not a Number=R tag

        int type  = bcf_hdr_id2type(call->hdr,BCF_HL_INFO,info->key);
        const char *key = bcf_hdr_int2id(call->hdr,BCF_DT_ID,info->key);
        nret = bcf_get_info_values(call->hdr, rec, key, &tmp_ori, &ntmp_ori, type);
        if ( nret<=0 ) continue;

        if ( nals_new==1 )
            bcf_update_info_int32(call->hdr, rec, key, tmp_ori, 1);     // has to be the REF, the order could not change
        else
        {
            for (j=0; j<nals_ori; j++)
            {
                int k = call->als_map[j];
                if ( k==-1 ) continue;   // to be dropped
                memcpy((char *)tmp_new+size*k, (char *)tmp_ori+size*j, size);
            }
            bcf_update_info_int32(call->hdr, rec, key, tmp_new, nals_new);
        }
    }

    // FORMAT fields
    for (i=0; i<rec->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &rec->d.fmt[i];
        int vlen = bcf_hdr_id2length(call->hdr,BCF_HL_FMT,fmt->id);
        if ( vlen!=BCF_VL_R ) continue; // not a Number=R tag

        int type = bcf_hdr_id2type(call->hdr,BCF_HL_FMT,fmt->id);
        const char *key = bcf_hdr_int2id(call->hdr,BCF_DT_ID,fmt->id);
        nret = bcf_get_format_values(call->hdr, rec, key, &tmp_ori, &ntmp_ori, type);
        if (nret<=0) continue;
        int nsmpl = bcf_hdr_nsamples(call->hdr);

        assert( nret==nals_ori*nsmpl );

        for (j=0; j<nsmpl; j++)
        {
            char *ptr_src = (char *)tmp_ori + j*nals_ori*size;
            char *ptr_dst = (char *)tmp_new + j*nals_new*size;
            int k;
            for (k=0; k<nals_ori; k++)
            {
                int l = call->als_map[k];
                if ( l==-1 ) continue;   // to be dropped
                memcpy(ptr_dst+size*l, ptr_src+size*k, size);
            }
        }
        bcf_update_format_int32(call->hdr, rec, key, tmp_new, nals_new*nsmpl);
    }

    call->PLs    = (int32_t*) tmp_new;
    call->mPLs   = ntmp_new;
    call->itmp   = (int32_t*) tmp_ori;
    call->n_itmp = ntmp_ori;
}


// NB: in this function we temporarily use calls->als_map for a different
// purpose to store mapping from new (target) alleles to original alleles.
//
static int mcall_constrain_alleles(call_t *call, bcf1_t *rec, int *unseen)
{
    assert( call->tgt_als->n );
    hts_expand(char*,call->tgt_als->n+1,call->nals,call->als);
    hts_expand(int,call->tgt_als->n+1,call->nals_map,call->als_map);
    hts_expand(int,(call->tgt_als->n+1)*(call->tgt_als->n+2)/2,call->npl_map,call->pl_map);

    int has_new = 0;

    int i, j, nals = 1;
    for (i=1; i<call->nals_map; i++) call->als_map[i] = -1;

    if ( vcmp_set_ref(call->vcmp, rec->d.allele[0], call->tgt_als->allele[0]) < 0 )
        error("The reference alleles are not compatible at %s:%d .. %s vs %s\n", call->hdr->id[BCF_DT_CTG][rec->rid].key,rec->pos+1,call->tgt_als->allele[0],rec->d.allele[0]);

    // create mapping from new to old alleles
    call->als[0] = call->tgt_als->allele[0];
    call->als_map[0] = 0;

    for (i=1; i<call->tgt_als->n; i++)
    {
        call->als[nals] = call->tgt_als->allele[i];
        j = vcmp_find_allele(call->vcmp, rec->d.allele+1, rec->n_allele - 1, call->tgt_als->allele[i]);
        if ( j>=0 )
        {
            // existing allele
            call->als_map[nals] = j+1;
        }
        else
        {
            // There is a new allele in targets which is not present in VCF.
            // We use the X allele to estimate PLs. Note that X may not be
            // present at multiallelic indels sites. In that case we use the
            // last allele anyway, because the least likely allele comes last
            // in mpileup's ALT output.
            call->als_map[nals] = (*unseen)>=0 ? *unseen : rec->n_allele - 1;
            has_new = 1;
        }
        nals++;
    }
    if ( *unseen )
    {
        call->als_map[nals] = *unseen;
        call->als[nals] = rec->d.allele[*unseen];
        nals++;
    }

    if ( !has_new && nals==rec->n_allele ) return 0;
    bcf_update_alleles(call->hdr, rec, (const char**)call->als, nals);

    // create mapping from new PL to old PL
    int k = 0;
    for (i=0; i<nals; i++)
    {
        for (j=0; j<=i; j++)
        {
            int a = call->als_map[i], b = call->als_map[j];
            call->pl_map[k++] = a>b ? a*(a+1)/2 + b : b*(b+1)/2 + a;
        }
    }

    // update PL
    call->nPLs = bcf_get_format_int32(call->hdr, rec, "PL", &call->PLs, &call->mPLs);
    int nsmpl  = bcf_hdr_nsamples(call->hdr);
    int npls_ori = call->nPLs / nsmpl;
    int npls_new = k;
    hts_expand(int32_t,npls_new*nsmpl,call->n_itmp,call->itmp);
    int *ori_pl = call->PLs, *new_pl = call->itmp;
    for (i=0; i<nsmpl; i++)
    {
        for (k=0; k<npls_new; k++)
        {
            new_pl[k] = ori_pl[call->pl_map[k]];
            if ( new_pl[k]==bcf_int32_missing && *unseen>=0 )
            {
                // missing value, and there is an unseen allele: identify the
                // alleles and use the lk of either AX or XX
                int k_ori = call->pl_map[k], ia, ib;
                bcf_gt2alleles(k_ori, &ia, &ib);
                k_ori = bcf_alleles2gt(ia,*unseen);
                if ( ori_pl[k_ori]==bcf_int32_missing ) k_ori = bcf_alleles2gt(ib,*unseen);
                if ( ori_pl[k_ori]==bcf_int32_missing ) k_ori = bcf_alleles2gt(*unseen,*unseen);
                new_pl[k] = ori_pl[k_ori];
            }
            if ( !k && new_pl[k]==bcf_int32_vector_end ) new_pl[k]=bcf_int32_missing;
        }
        ori_pl += npls_ori;
        new_pl += npls_new;
    }
    bcf_update_format_int32(call->hdr, rec, "PL", call->itmp, npls_new*nsmpl);

    // update QS, use temporarily call->GPs to store the values
    int nqs = bcf_get_info_float(call->hdr, rec, "QS", &call->smpl_grp[0].qsum, &call->smpl_grp[0].nqsum);
    hts_expand(float,nals,call->nGPs,call->GPs);
    for (i=0; i<nals; i++)
        call->GPs[i] = call->als_map[i]<nqs ? call->smpl_grp[0].qsum[call->als_map[i]] : 0;
    bcf_update_info_float(call->hdr, rec, "QS", call->GPs, nals);

    // update any Number=R tags
    void *tmp_ori = call->itmp, *tmp_new = call->PLs;  // reusing PLs storage which is not used at this point
    int ntmp_ori = call->n_itmp, ntmp_new = call->mPLs;
    for (i=0; i<rec->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &rec->d.fmt[i];
        int vlen = bcf_hdr_id2length(call->hdr,BCF_HL_FMT,fmt->id);
        if ( vlen!=BCF_VL_R ) continue; // not a Number=R tag

        // NB:works only for BCF_HT_INT and BCF_HT_REAL
        int type = bcf_hdr_id2type(call->hdr,BCF_HL_FMT,fmt->id);
        assert( type==BCF_HT_INT || type==BCF_HT_REAL );
        assert( sizeof(float)==sizeof(int32_t) );

        const char *key = bcf_hdr_int2id(call->hdr,BCF_DT_ID,fmt->id);
        int nret = bcf_get_format_values(call->hdr, rec, key, &tmp_ori, &ntmp_ori, type);
        if (nret<=0) continue;
        int nsmpl = bcf_hdr_nsamples(call->hdr);
        int size1 = sizeof(float);
        hts_expand(float, nsmpl * nals, ntmp_new, tmp_new);
        for (j=0; j<nsmpl; j++)
        {
            uint8_t *ptr_ori = (uint8_t *) tmp_ori + j*size1*fmt->n;
            uint8_t *ptr_new = (uint8_t *) tmp_new + j*nals*size1;
            for (k=0; k<nals; k++)
            {
                uint8_t *dst = ptr_new + size1*k;
                uint8_t *src = ptr_ori + size1*call->als_map[k];
                memcpy(dst,src,size1);
            }
        }
        nret = bcf_update_format(call->hdr, rec, key, tmp_new, nsmpl*nals, type);
        assert( nret==0 );
    }
    call->PLs    = (int32_t*) tmp_new;
    call->mPLs   = ntmp_new;
    call->itmp   = (int32_t*) tmp_ori;
    call->n_itmp = ntmp_ori;

    if ( *unseen ) *unseen = nals-1;
    return 0;
}


/**
  *  This function implements the multiallelic calling model. It has two major parts:
  *   1) determine the most likely set of alleles and calculate the quality of ref/non-ref site
  *   2) determine and set the genotypes
  *  In various places in between, the BCF record gets updated.
  */
int mcall(call_t *call, bcf1_t *rec)
{
    int i,j, unseen = call->unseen;

    // Force alleles when calling genotypes given alleles was requested
    if ( call->flag & CALL_CONSTR_ALLELES && mcall_constrain_alleles(call, rec, &unseen)!=0 ) return -2;

    int nsmpl    = bcf_hdr_nsamples(call->hdr);
    int nals_ori = rec->n_allele;
    hts_expand(int,nals_ori,call->nac,call->ac);
    hts_expand(int,nals_ori,call->nals_map,call->als_map);
    hts_expand(int,nals_ori*(nals_ori+1)/2,call->npl_map,call->pl_map);

    // Get the genotype likelihoods
    call->nPLs = bcf_get_format_int32(call->hdr, rec, "PL", &call->PLs, &call->mPLs);
    if ( call->nPLs!=nsmpl*nals_ori*(nals_ori+1)/2 && call->nPLs!=nsmpl*nals_ori )  // a mixture of diploid and haploid or haploid only
        error("Wrong number of PL fields? nals=%d npl=%d\n", nals_ori,call->nPLs);

    // Convert PLs to probabilities
    int ngts_ori = nals_ori*(nals_ori+1)/2;
    hts_expand(double, call->nPLs, call->npdg, call->pdg);
    set_pdg(call->pl2p, call->PLs, call->pdg, nsmpl, ngts_ori, unseen);

    // Get sum of qualities, serves as an AF estimate, f_x = QS/N in Eq. 1 in call-m math notes.
    if ( call->nsmpl_grp == 1  )
    {
        int nqs = bcf_get_info_float(call->hdr, rec, "QS", &call->smpl_grp[0].qsum, &call->smpl_grp[0].nqsum);
        if ( nqs<=0 ) error("The QS annotation not present at %s:%d\n", bcf_seqname(call->hdr,rec),rec->pos+1);
        if ( nqs < nals_ori )
        {
            // Some of the listed alleles do not have the corresponding QS field. This is
            // typically ref-only site with <*> in ALT.
            hts_expand(float,nals_ori,call->smpl_grp[0].nqsum,call->smpl_grp[0].qsum);
            for (i=nqs; i<nals_ori; i++) call->smpl_grp[0].qsum[i] = 0;
        }
    }
    else
    {
        for (j=0; j<call->nsmpl_grp; j++)
        {
            hts_expand(float,nals_ori,call->smpl_grp[j].nqsum,call->smpl_grp[j].qsum);
            memset(call->smpl_grp[j].qsum, 0, sizeof(float)*nals_ori);
        }

        // Use FORMAT/AD or FORMAT/QS
        int nad = bcf_get_format_int32(call->hdr, rec, call->sample_groups_tag, &call->ADs, &call->nADs);
        if ( nad<1 ) error("Error: FORMAT/%s is required with the -G option, mpileup must be run with \"-a AD\" or \"-a QS\"\n",call->sample_groups_tag);
        nad /= bcf_hdr_nsamples(call->hdr);
        for (i=0; i<call->nsmpl_grp; i++)
        {
            int is;
            smpl_grp_t *grp = &call->smpl_grp[i];
            hts_expand(float,nals_ori,grp->nqsum,grp->qsum);
            for (j=0; j<nals_ori; j++) grp->qsum[j] = 0;
            for (is=0; is<grp->nsmpl; is++)
            {
                int ismpl = grp->smpl[is];
                int32_t *ptr = call->ADs + ismpl*nad;
                float sum = 0;
                for (j=0; j<nad; j++)
                {
                    if ( ptr[j]==bcf_int32_vector_end ) break;
                    if ( ptr[j]!=bcf_int32_missing ) sum += ptr[j];
                }
                if ( sum )
                {
                    for (j=0; j<nad; j++)
                    {
                        if ( ptr[j]==bcf_int32_vector_end ) break;
                        if ( ptr[j]!=bcf_int32_missing ) grp->qsum[j] += ptr[j]/sum;
                    }
                }
            }
        }
    }

    // If available, take into account reference panel AFs
    if ( call->prior_AN && bcf_get_info_int32(call->hdr, rec, call->prior_AN ,&call->ac, &call->nac)==1 )
    {
        int an = call->ac[0];   // number of alleles total, procede only if not zero; reuse call->ac
        if ( an > 0 && bcf_get_info_int32(call->hdr, rec, call->prior_AC ,&call->ac, &call->nac)==nals_ori-1 )    // number of ALT alleles
        {
            int ac0 = an;       // this will become the number of REFs
            for (i=0; i<nals_ori-1; i++)
            {
                if ( call->ac[i]==bcf_int32_vector_end ) break;
                if ( call->ac[i]==bcf_int32_missing ) continue;
                ac0 -= call->ac[i];

                // here an*0.5 is the number of samples in the populatio and ac*0.5 is the AF weighted by the number of samples
                for (j=0; j<call->nsmpl_grp; j++)
                    call->smpl_grp[j].qsum[i+1] = (call->smpl_grp[j].qsum[i+1] + 0.5*call->ac[i]) / (call->smpl_grp[j].nsmpl + 0.5*an);
            }
            if ( ac0<0 ) error("Incorrect %s,%s values at %s:%d\n", call->prior_AN,call->prior_AC,bcf_seqname(call->hdr,rec),rec->pos+1);
            for (j=0; j<call->nsmpl_grp; j++)
                call->smpl_grp[j].qsum[0] = (call->smpl_grp[j].qsum[0] + 0.5*ac0) / (call->smpl_grp[j].nsmpl + 0.5*an);
        }
    }

    // normalize so that QS sums to 1 for each group
    for (j=0; j<call->nsmpl_grp; j++)
    {
        float sum = 0;
        for (i=0; i<nals_ori; i++) sum += call->smpl_grp[j].qsum[i];
        if ( sum ) for (i=0; i<nals_ori; i++) call->smpl_grp[j].qsum[i] /= sum;
    }

    bcf_update_info_int32(call->hdr, rec, "QS", NULL, 0);      // remove QS tag

    if ( nals_ori > 8*sizeof(call->als_new) )
    {
        fprintf(stderr,"Too many alleles at %s:%"PRId64", skipping.\n", bcf_seqname(call->hdr,rec),(int64_t) rec->pos+1);
        return 0;
    }

    // For each group find the best combination of alleles
    call->als_new = 0;
    double ref_lk = -HUGE_VAL, lk_sum = -HUGE_VAL, max_qual = -HUGE_VAL;
    for (j=0; j<call->nsmpl_grp; j++)
    {
        smpl_grp_t *grp = &call->smpl_grp[j];
        mcall_find_best_alleles(call, nals_ori, grp);
        call->als_new |= grp->als;
        if ( grp->max_lk==-HUGE_VAL ) continue;
        double qual = -4.343*(grp->ref_lk - logsumexp2(grp->lk_sum,grp->ref_lk));
        if ( max_qual < qual )
        {
            max_qual = qual;
            lk_sum = grp->lk_sum;
            ref_lk = grp->ref_lk;
        }
    }

    // Make sure the REF allele is always present
    if ( !(call->als_new&1) ) call->als_new |= 1;

    int is_variant = call->als_new==1 ? 0 : 1;
    if ( call->flag & CALL_VARONLY && !is_variant ) return 0;

    call->nals_new = 0;
    for (i=0; i<nals_ori; i++)
    {
        if ( i>0 && i==unseen ) continue;
        if ( call->flag & CALL_KEEPALT ) call->als_new |= 1<<i;
        if ( call->als_new & (1<<i) ) call->nals_new++;
    }

    init_allele_trimming_maps(call,nals_ori,call->als_new);

    int nAC = 0;
    if ( call->als_new==1 )   // only REF allele on output
    {
        mcall_set_ref_genotypes(call,nals_ori);
        bcf_update_format_int32(call->hdr, rec, "PL", NULL, 0);    // remove PL, useless now
    }
    else if ( !is_variant )
    {
        mcall_set_ref_genotypes(call,nals_ori);     // running with -A, prevent mcall_call_genotypes from putting some ALT back
        mcall_trim_and_update_PLs(call, rec, nals_ori, call->nals_new);
    }
    else
    {
        // The most likely set of alleles includes non-reference allele (or was enforced), call genotypes.
        // Note that it is a valid outcome if the called genotypes exclude some of the ALTs.
        int ngts_new = call->nals_new*(call->nals_new+1)/2;
        hts_expand(float,ngts_new*nsmpl,call->nGPs,call->GPs);
        for (i=0; i<call->nals_new; i++) call->ac[i] = 0;

        if ( call->flag & CALL_CONSTR_TRIO && call->nals_new>4 )
        {
            fprintf(stderr,"Too many alleles at %s:%"PRId64", skipping.\n", bcf_seqname(call->hdr,rec),(int64_t) rec->pos+1);
            return 0;
        }
        if ( call->output_tags & (CALL_FMT_GQ|CALL_FMT_GP) )
        {
            memset(call->GPs,0,nsmpl*ngts_new*sizeof(*call->GPs));
            memset(call->GQs,0,nsmpl*sizeof(*call->GQs));
        }
        for (i=0; i<call->nsmpl_grp; i++)
        {
            if ( call->flag & CALL_CONSTR_TRIO )
                error("todo: constrained trio calling temporarily disabled\n");   //mcall_call_trio_genotypes(call,rec,nals,&call->smpl_grp[i]);
            else
                mcall_call_genotypes(call,nals_ori,&call->smpl_grp[i]);
        }

        // Skip the site if all samples are 0/0. This can happen occasionally.
        for (i=1; i<call->nals_new; i++) nAC += call->ac[i];
        if ( !nAC && call->flag & CALL_VARONLY ) return 0;

        if ( call->output_tags & CALL_FMT_GP )
            bcf_update_format_float(call->hdr, rec, "GP", call->GPs, nsmpl*ngts_new);
        if ( call->output_tags & CALL_FMT_GQ )
            bcf_update_format_int32(call->hdr, rec, "GQ", call->GQs, nsmpl);

        mcall_trim_and_update_PLs(call,rec,nals_ori,call->nals_new);
    }
    if ( nals_ori!=call->nals_new )
        mcall_trim_and_update_numberR(call,rec,nals_ori,call->nals_new);

    // Set QUAL
    if ( nAC )
    {
        // Quality of a variant site. fabs() to avoid negative zeros in VCF output when CALL_KEEPALT is set
        rec->qual = max_qual;
    }
    else
    {
        // Set the quality of a REF site
        if ( lk_sum!=-HUGE_VAL )  // no support from (high quality) reads, so QUAL=1-prior
            rec->qual = -4.343*(lk_sum - logsumexp2(lk_sum,ref_lk));
        else if ( call->ac[0] )
            rec->qual = call->theta ? -4.343*call->theta : 0;
        else
            bcf_float_set_missing(rec->qual);
    }

    // AC, AN
    if ( call->nals_new>1 ) bcf_update_info_int32(call->hdr, rec, "AC", call->ac+1, call->nals_new-1);
    nAC += call->ac[0];
    bcf_update_info_int32(call->hdr, rec, "AN", &nAC, 1);

    // Remove unused alleles
    hts_expand(char*,call->nals_new,call->nals,call->als);
    for (i=0; i<nals_ori; i++)
        if ( call->als_map[i]>=0 ) call->als[call->als_map[i]] = rec->d.allele[i];
    bcf_update_alleles(call->hdr, rec, (const char**)call->als, call->nals_new);
    bcf_update_genotypes(call->hdr, rec, call->gts, nsmpl*2);

    // DP4 and PV4 tags
    if ( bcf_get_info_float(call->hdr, rec, "I16", &call->anno16, &call->n16)==16 )
    {
        int32_t dp[4]; dp[0] = call->anno16[0]; dp[1] = call->anno16[1]; dp[2] = call->anno16[2]; dp[3] = call->anno16[3];
        bcf_update_info_int32(call->hdr, rec, "DP4", dp, 4);

        int32_t mq = (call->anno16[8]+call->anno16[10])/(call->anno16[0]+call->anno16[1]+call->anno16[2]+call->anno16[3]);
        bcf_update_info_int32(call->hdr, rec, "MQ", &mq, 1);

        if ( call->output_tags & CALL_FMT_PV4 )
        {
            anno16_t a;
            float tmpf[4];
            int is_tested = test16(call->anno16, &a) >= 0 && a.is_tested ? 1 : 0;
            if ( is_tested )
            {
                for (i=0; i<4; i++) tmpf[i] = a.p[i];
                bcf_update_info_float(call->hdr, rec, "PV4", tmpf, 4);
            }
        }
    }

    bcf_update_info_int32(call->hdr, rec, "I16", NULL, 0);     // remove I16 tag

    return call->nals_new;
}

