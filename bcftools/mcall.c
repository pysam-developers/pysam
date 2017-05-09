/*  mcall.c -- multiallelic and rare variant calling.

    Copyright (C) 2012-2016 Genome Research Ltd.

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

// Using priors for GTs does not seem to be mathematically justified. Although
// it seems effective in removing false calls, it also flips a significant
// proportion of HET genotypes. Better is to filter by FORMAT/GQ using
// `bcftools filter`.
#define USE_PRIOR_FOR_GTS 0

// Go with uniform PLs for samples with no coverage. If unset, missing
// genotypes is reported instead.
#define FLAT_PDG_FOR_MISSING 0

// Estimate QS (combined quality and allele frequencies) from PLs
#define QS_FROM_PDG 0


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

void mcall_init(call_t *call)
{
    call_init_pl2p(call);

    call->nqsum = 5;
    call->qsum  = (float*) malloc(sizeof(float)*call->nqsum); // will be expanded later if ncessary
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
        bcf_hdr_append(call->hdr,"##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Phred-scaled genotype posterior probabilities\">");
    if ( call->output_tags & (CALL_FMT_GQ|CALL_FMT_GP) )
        call->GQs = (int32_t*) malloc(sizeof(int32_t)*bcf_hdr_nsamples(call->hdr));
    bcf_hdr_append(call->hdr,"##INFO=<ID=ICB,Number=1,Type=Float,Description=\"Inbreeding Coefficient Binomial test (bigger is better)\">");
    bcf_hdr_append(call->hdr,"##INFO=<ID=HOB,Number=1,Type=Float,Description=\"Bias in the number of HOMs number (smaller is better)\">");
    bcf_hdr_append(call->hdr,"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes for each ALT allele, in the same order as listed\">");
    bcf_hdr_append(call->hdr,"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    bcf_hdr_append(call->hdr,"##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases\">");
    bcf_hdr_append(call->hdr,"##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Average mapping quality\">");

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

    return;
}

void mcall_destroy(call_t *call)
{
    if (call->vcmp) vcmp_destroy(call->vcmp);
    free(call->itmp);
    mcall_destroy_trios(call);
    free(call->GPs);
    free(call->GLs);
    free(call->GQs);
    free(call->anno16);
    free(call->PLs);
    free(call->qsum);
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

/*
    Allele frequency estimated as:
        #A  = \sum_i (2*P_AA + P_AB)
        F_A = #A / ( #A + #B )
    where i runs across all samples
*/
void estimate_qsum(call_t *call, bcf1_t *rec)
{
    double *pdg  = call->pdg;
    int ngts = rec->n_allele*(rec->n_allele+1)/2;
    int i,nsmpl = bcf_hdr_nsamples(call->hdr);

    hts_expand(float,rec->n_allele,call->nqsum,call->qsum);
    for (i=0; i<rec->n_allele; i++) call->qsum[i] = 0;

    for (i=0; i<nsmpl; i++)
    {
        int a, b, k = 0;
        for (a=0; a<rec->n_allele; a++)
        {
            for (b=0; b<=a; b++)
            {
                call->qsum[a] += pdg[k];
                call->qsum[b] += pdg[k];
                k++;
            }
        }
        pdg += ngts;
    }
    float sum = 0;
    for (i=0; i<rec->n_allele; i++) sum += call->qsum[i];
    if ( sum ) for (i=0; i<rec->n_allele; i++) call->qsum[i] /= sum;
}

// Create mapping between old and new (trimmed) alleles
void init_allele_trimming_maps(call_t *call, int als, int nals)
{
    int i, j;

    // als_map: old(i) -> new(j)
    for (i=0, j=0; i<nals; i++)
    {
        if ( als & 1<<i ) call->als_map[i] = j++;
        else call->als_map[i] = -1;
    }

    if ( !call->pl_map ) return;

    // pl_map: new(k) -> old(l)
    int k = 0, l = 0;
    for (i=0; i<nals; i++)
    {
        for (j=0; j<=i; j++)
        {
            if ( (als & 1<<i) && (als & 1<<j) ) call->pl_map[k++] = l;
            l++;
        }
    }
}

double binom_dist(int N, double p, int k)
{
    int mean = (int) (N*p);
    if ( mean==k ) return 1.0;

    double log_p = (k-mean)*log(p) + (mean-k)*log(1.0-p);
    if ( k > N - k ) k = N - k;
    if ( mean > N - mean ) mean = N - mean;

    if ( k < mean ) { int tmp = k; k = mean; mean = tmp; }
    double diff = k - mean;

    double val = 1.0;
    int i;
    for (i=0; i<diff; i++)
        val = val * (N-mean-i) / (k-i);

    return exp(log_p)/val;
}


// Inbreeding Coefficient, binomial test
float calc_ICB(int nref, int nalt, int nhets, int ndiploid)
{
    if ( !nref || !nalt || !ndiploid ) return HUGE_VAL;

    double fref = (double)nref/(nref+nalt); // fraction of reference allelels
    double falt = (double)nalt/(nref+nalt); // non-ref als
    double q = 2*fref*falt;                 // probability of a het, assuming HWE
    double mean = q*ndiploid;

    //fprintf(stderr,"\np=%e N=%d k=%d  .. nref=%d nalt=%d nhets=%d ndiploid=%d\n", q,ndiploid,nhets, nref,nalt,nhets,ndiploid);

    // Can we use normal approximation? The second condition is for performance only
    // and is not well justified.
    if ( (mean>10 && (1-q)*ndiploid>10 ) || ndiploid>200 )
    {
        //fprintf(stderr,"out: mean=%e  p=%e\n", mean,exp(-0.5*(nhets-mean)*(nhets-mean)/(mean*(1-q))));
        return exp(-0.5*(nhets-mean)*(nhets-mean)/(mean*(1-q)));
    }

    return binom_dist(ndiploid, q, nhets);
}

float calc_HOB(int nref, int nalt, int nhets, int ndiploid)
{
    if ( !nref || !nalt || !ndiploid ) return HUGE_VAL;

    double fref = (double)nref/(nref+nalt); // fraction of reference allelels
    double falt = (double)nalt/(nref+nalt); // non-ref als
    return fabs((double)nhets/ndiploid - 2*fref*falt);
}

/**
  *  log(sum_i exp(a_i))
  */
// static inline double logsumexp(double *vals, int nvals)
// {
//     int i;
//     double max_exp = vals[0];
//     for (i=1; i<nvals; i++)
//         if ( max_exp < vals[i] ) max_exp = vals[i];

//     double sum = 0;
//     for (i=0; i<nvals; i++)
//         sum += exp(vals[i] - max_exp);

//     return log(sum) + max_exp;
// }
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
     if ( max_lk<lk_tot ) { max_lk = lk_tot; max_als = (als); } \
     if ( sum ) lk_sum = logsumexp2(lk_tot,lk_sum); \
}

#define SWAP(type_t,x,y) {type_t tmp; tmp = x; x = y; y = tmp; }

// Determine the most likely combination of alleles. In this implementation,
// at most tri-allelic sites are considered. Returns the number of alleles.
static int mcall_find_best_alleles(call_t *call, int nals, int *out_als)
{
    int ia,ib,ic;   // iterators over up to three alleles
    int max_als=0;  // most likely combination of alleles
    double ref_lk = 0, max_lk = -HUGE_VAL; // likelihood of the reference and of most likely combination of alleles
    double lk_sum = -HUGE_VAL;    // for normalizing the likelihoods
    int nsmpl = bcf_hdr_nsamples(call->hdr);
    int ngts  = nals*(nals+1)/2;

    // Single allele
    for (ia=0; ia<nals; ia++)
    {
        double lk_tot  = 0;
        int lk_tot_set = 0;
        int iaa = (ia+1)*(ia+2)/2-1;    // index in PL which corresponds to the homozygous "ia/ia" genotype
        int isample;
        double *pdg = call->pdg + iaa;
        for (isample=0; isample<nsmpl; isample++)
        {
            if ( *pdg ) { lk_tot += log(*pdg); lk_tot_set = 1; }
            pdg += ngts;
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
            if ( call->qsum[ia]==0 ) continue;
            int iaa = (ia+1)*(ia+2)/2-1;
            for (ib=0; ib<ia; ib++)
            {
                if ( call->qsum[ib]==0 ) continue;
                double lk_tot  = 0;
                int lk_tot_set = 0;
                double fa  = call->qsum[ia]/(call->qsum[ia]+call->qsum[ib]);
                double fb  = call->qsum[ib]/(call->qsum[ia]+call->qsum[ib]);
                double fa2 = fa*fa;
                double fb2 = fb*fb;
                double fab = 2*fa*fb;
                int isample, ibb = (ib+1)*(ib+2)/2-1, iab = iaa - ia + ib;
                double *pdg  = call->pdg;
                for (isample=0; isample<nsmpl; isample++)
                {
                    double val = 0;
                    if ( !call->ploidy || call->ploidy[isample]==2 )
                        val = fa2*pdg[iaa] + fb2*pdg[ibb] + fab*pdg[iab];
                    else if ( call->ploidy && call->ploidy[isample]==1 )
                        val = fa*pdg[iaa] + fb*pdg[ibb];
                    if ( val ) { lk_tot += log(val); lk_tot_set = 1; }
                    pdg += ngts;
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
            if ( call->qsum[ia]==0 ) continue;
            int iaa = (ia+1)*(ia+2)/2-1;
            for (ib=0; ib<ia; ib++)
            {
                if ( call->qsum[ib]==0 ) continue;
                int ibb = (ib+1)*(ib+2)/2-1;
                int iab = iaa - ia + ib;
                for (ic=0; ic<ib; ic++)
                {
                    if ( call->qsum[ic]==0 ) continue;
                    double lk_tot  = 0;
                    int lk_tot_set = 1;
                    double fa  = call->qsum[ia]/(call->qsum[ia]+call->qsum[ib]+call->qsum[ic]);
                    double fb  = call->qsum[ib]/(call->qsum[ia]+call->qsum[ib]+call->qsum[ic]);
                    double fc  = call->qsum[ic]/(call->qsum[ia]+call->qsum[ib]+call->qsum[ic]);
                    double fa2 = fa*fa;
                    double fb2 = fb*fb;
                    double fc2 = fc*fc;
                    double fab = 2*fa*fb, fac = 2*fa*fc, fbc = 2*fb*fc;
                    int isample, icc = (ic+1)*(ic+2)/2-1;
                    int iac = iaa - ia + ic, ibc = ibb - ib + ic;
                    double *pdg = call->pdg;
                    for (isample=0; isample<nsmpl; isample++)
                    {
                        double val = 0;
                        if ( !call->ploidy || call->ploidy[isample]==2 )
                            val = fa2*pdg[iaa] + fb2*pdg[ibb] + fc2*pdg[icc] + fab*pdg[iab] + fac*pdg[iac] + fbc*pdg[ibc];
                        else if ( call->ploidy && call->ploidy[isample]==1 )
                            val = fa*pdg[iaa] + fb*pdg[ibb] + fc*pdg[icc];
                        if ( val ) { lk_tot += log(val); lk_tot_set = 1; }
                        pdg += ngts;
                    }
                    if ( ia!=0 ) lk_tot += call->theta;    // the prior
                    if ( ib!=0 ) lk_tot += call->theta;    // the prior
                    if ( ic!=0 ) lk_tot += call->theta;    // the prior
                    UPDATE_MAX_LKs(1<<ia|1<<ib|1<<ic, lk_tot_set);
                }
            }
        }
    }

    call->ref_lk = ref_lk;
    call->lk_sum = lk_sum;
    *out_als = max_als;

    int i, n = 0;
    for (i=0; i<nals; i++) if ( max_als & 1<<i) n++;

    return n;
}

static void mcall_set_ref_genotypes(call_t *call, int nals)
{
    int i;
    int ngts  = nals*(nals+1)/2;
    int nsmpl = bcf_hdr_nsamples(call->hdr);

    for (i=0; i<nals; i++) call->ac[i] = 0;
    call->nhets = 0;
    call->ndiploid = 0;

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

static void mcall_call_genotypes(call_t *call, bcf1_t *rec, int nals, int nout_als, int out_als)
{
    int ia, ib, i;
    int ngts  = nals*(nals+1)/2;
    int nsmpl = bcf_hdr_nsamples(call->hdr);
    int nout_gts = nout_als*(nout_als+1)/2;
    hts_expand(float,nout_gts*nsmpl,call->nGPs,call->GPs);

    for (i=0; i<nout_als; i++) call->ac[i] = 0;
    call->nhets = 0;
    call->ndiploid = 0;

    #if USE_PRIOR_FOR_GTS
        float prior = exp(call->theta);
    #endif
    float *gps  = call->GPs - nout_gts;
    double *pdg = call->pdg - ngts;
    int *gts  = call->gts - 2;

    int isample;
    for (isample = 0; isample < nsmpl; isample++)
    {
        int ploidy = call->ploidy ? call->ploidy[isample] : 2;
        assert( ploidy>=0 && ploidy<=2 );

        pdg += ngts;
        gts += 2;
        gps += nout_gts;

        if ( !ploidy )
        {
            gts[0] = bcf_gt_missing;
            gts[1] = bcf_int32_vector_end;
            gps[0] = -1;
            continue;
        }

        #if !FLAT_PDG_FOR_MISSING
            // Skip samples with zero depth, they have all pdg's equal to 0
            for (i=0; i<ngts; i++) if ( pdg[i]!=0.0 ) break;
            if ( i==ngts )
            {
                gts[0] = bcf_gt_missing;
                gts[1] = ploidy==2 ? bcf_gt_missing : bcf_int32_vector_end;
                gps[0] = -1;
                continue;
            }
        #endif

        if ( ploidy==2 ) call->ndiploid++;

        // Default fallback for the case all LKs are the same
        gts[0] = bcf_gt_unphased(0);
        gts[1] = ploidy==2 ? bcf_gt_unphased(0) : bcf_int32_vector_end;

        // Non-zero depth, determine the most likely genotype
        double best_lk = 0;
        for (ia=0; ia<nals; ia++)
        {
            if ( !(out_als & 1<<ia) ) continue;     // ia-th allele not in the final selection, skip
            int iaa = (ia+1)*(ia+2)/2-1;            // PL index of the ia/ia genotype
            double lk = ploidy==2 ? pdg[iaa]*call->qsum[ia]*call->qsum[ia] : pdg[iaa]*call->qsum[ia];
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
            for (ia=0; ia<nals; ia++)
            {
                if ( !(out_als & 1<<ia) ) continue;
                int iaa = (ia+1)*(ia+2)/2-1;
                for (ib=0; ib<ia; ib++)
                {
                    if ( !(out_als & 1<<ib) ) continue;
                    int iab = iaa - ia + ib;
                    double lk = 2*pdg[iab]*call->qsum[ia]*call->qsum[ib];
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
            if ( gts[0] != gts[1] ) call->nhets++;
        }
        else
            gts[1] = bcf_int32_vector_end;

        call->ac[ bcf_gt_allele(gts[0]) ]++;
        if ( gts[1]!=bcf_int32_vector_end ) call->ac[ bcf_gt_allele(gts[1]) ]++;
    }
    if ( call->output_tags & (CALL_FMT_GQ|CALL_FMT_GP) )
    {
        double max, sum;
        for (isample=0; isample<nsmpl; isample++)
        {
            gps = call->GPs + isample*nout_gts;

            int nmax;
            if ( call->ploidy )
            {
                if ( call->ploidy[isample]==2 ) nmax = nout_gts;
                else if ( call->ploidy[isample]==1 ) nmax = nout_als;
                else nmax = 0;
            }
            else nmax = nout_gts;

            max = gps[0];
            if ( max<0 || nmax==0 )
            {
                // no call
                if ( call->output_tags & CALL_FMT_GP )
                {
                    for (i=0; i<nmax; i++) gps[i] = 0;
                    if ( nmax==0 ) { bcf_float_set_missing(gps[i]); nmax++; }
                    if ( nmax < nout_gts ) bcf_float_set_vector_end(gps[nmax]);
                }
                call->GQs[isample] = 0;
                continue;
            }
            sum = gps[0];
            for (i=1; i<nmax; i++)
            {
                if ( max < gps[i] ) max = gps[i];
                sum += gps[i];
            }
            max = -4.34294*log(1 - max/sum);
            call->GQs[isample] = max<=INT8_MAX ? max : INT8_MAX;
            if ( call->output_tags & CALL_FMT_GP )
            {
                assert( max );
                for (i=0; i<nmax; i++) gps[i] = (int)(-4.34294*log(gps[i]/sum));
                if ( nmax < nout_gts ) bcf_float_set_vector_end(gps[nmax]);
            }
        }
    }
    if ( call->output_tags & CALL_FMT_GP )
        bcf_update_format_float(call->hdr, rec, "GP", call->GPs, nsmpl*nout_gts);
    if ( call->output_tags & CALL_FMT_GQ )
        bcf_update_format_int32(call->hdr, rec, "GQ", call->GQs, nsmpl);
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
static void mcall_call_trio_genotypes(call_t *call, bcf1_t *rec, int nals, int nout_als, int out_als)
{
    int ia, ib, i;
    int nsmpl    = bcf_hdr_nsamples(call->hdr);
    int ngts     = nals*(nals+1)/2;
    int nout_gts = nout_als*(nout_als+1)/2;
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

        double sum_lk  = 0;
        double best_lk = 0;
        for (ia=0; ia<nals; ia++)
        {
            if ( !(out_als & 1<<ia) ) continue;     // ia-th allele not in the final selection, skip
            int iaa   = bcf_alleles2gt(ia,ia);      // PL index of the ia/ia genotype
            int idx   = bcf_alleles2gt(call->als_map[ia],call->als_map[ia]);
            double lk = ploidy==2 ? pdg[iaa]*call->qsum[ia]*call->qsum[ia] : pdg[iaa]*call->qsum[ia];
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
                if ( !(out_als & 1<<ia) ) continue;
                for (ib=0; ib<ia; ib++)
                {
                    if ( !(out_als & 1<<ib) ) continue;
                    int iab   = bcf_alleles2gt(ia,ib);
                    int idx   = bcf_alleles2gt(call->als_map[ia],call->als_map[ib]);
                    double lk = 2*pdg[iab]*call->qsum[ia]*call->qsum[ib];
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
        int ntrio = call->ntrio[fam->type][nout_als];
        uint16_t *trio = call->trio[fam->type][nout_als];

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

static void mcall_trim_PLs(call_t *call, bcf1_t *rec, int nals, int nout_als, int out_als)
{
    int ngts  = nals*(nals+1)/2;
    int npls_src = ngts, npls_dst = nout_als*(nout_als+1)/2;     // number of PL values in diploid samples, ori and new
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
            for (ia=0; ia<nout_als; ia++)
            {
                int isrc = (ia+1)*(ia+2)/2-1;
                pls_dst[ia] = pls_src[ call->pl_map[isrc] ];
            }
            if ( ia<npls_dst ) pls_dst[ia] = bcf_int32_vector_end;
        }
        else
        {
            pls_dst[0] = bcf_int32_missing;
            pls_dst[1] = bcf_int32_vector_end;  // relying on nout_als>1 in mcall()
        }
        pls_src += npls_src;
        pls_dst += npls_dst;
    }
    bcf_update_format_int32(call->hdr, rec, "PL", call->PLs, npls_dst*nsmpl);
}

void mcall_trim_numberR(call_t *call, bcf1_t *rec, int nals, int nout_als, int out_als)
{
    if ( nals==nout_als ) return;

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

        if ( nout_als==1 )
            bcf_update_info_int32(call->hdr, rec, key, tmp_ori, 1);     // has to be the REF, the order could not change
        else
        {
            for (j=0; j<nals; j++)
            {
                int k = call->als_map[j];
                if ( k==-1 ) continue;   // to be dropped
                memcpy((char *)tmp_new+size*k, (char *)tmp_ori+size*j, size);
            }
            bcf_update_info_int32(call->hdr, rec, key, tmp_new, nout_als);
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

        assert( nret==nals*nsmpl );

        for (j=0; j<nsmpl; j++)
        {
            char *ptr_src = (char *)tmp_ori + j*nals*size;
            char *ptr_dst = (char *)tmp_new + j*nout_als*size;
            int k;
            for (k=0; k<nals; k++)
            {
                int l = call->als_map[k];
                if ( l==-1 ) continue;   // to be dropped
                memcpy(ptr_dst+size*l, ptr_src+size*k, size);
            }
        }
        bcf_update_format_int32(call->hdr, rec, key, tmp_new, nout_als*nsmpl);
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
    bcf_sr_regions_t *tgt = call->srs->targets;
    if ( tgt->nals>5 ) error("Maximum accepted number of alleles is 5, got %d\n", tgt->nals);
    hts_expand(char*,tgt->nals+1,call->nals,call->als);

    int has_new = 0;

    int i, j, nals = 1;
    for (i=1; i<call->nals_map; i++) call->als_map[i] = -1;

    if ( vcmp_set_ref(call->vcmp, rec->d.allele[0], tgt->als[0]) < 0 )
        error("The reference alleles are not compatible at %s:%d .. %s vs %s\n", call->hdr->id[BCF_DT_CTG][rec->rid].key,rec->pos+1,tgt->als[0],rec->d.allele[0]);

    // create mapping from new to old alleles
    call->als[0] = tgt->als[0];
    call->als_map[0] = 0;

    for (i=1; i<tgt->nals; i++)
    {
        call->als[nals] = tgt->als[i];
        j = vcmp_find_allele(call->vcmp, rec->d.allele+1, rec->n_allele - 1, tgt->als[i]);

        if ( j+1==*unseen ) { fprintf(stderr,"fixme? Cannot constrain to %s\n",tgt->als[i]); return -1; }
        
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

    // update QS
    float qsum[5];
    int nqs = bcf_get_info_float(call->hdr, rec, "QS", &call->qsum, &call->nqsum);
    for (i=0; i<nals; i++)
        qsum[i] = call->als_map[i]<nqs ? call->qsum[call->als_map[i]] : 0;
    bcf_update_info_float(call->hdr, rec, "QS", qsum, nals);

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
    int i, unseen = call->unseen;

    // Force alleles when calling genotypes given alleles was requested
    if ( call->flag & CALL_CONSTR_ALLELES && mcall_constrain_alleles(call, rec, &unseen)!=0 ) return -2;

    int nsmpl = bcf_hdr_nsamples(call->hdr);
    int nals  = rec->n_allele;
    hts_expand(int,nals,call->nac,call->ac);
    hts_expand(int,nals,call->nals_map,call->als_map);
    hts_expand(int,nals*(nals+1)/2,call->npl_map,call->pl_map);

    // Get the genotype likelihoods
    call->nPLs = bcf_get_format_int32(call->hdr, rec, "PL", &call->PLs, &call->mPLs);
    if ( call->nPLs!=nsmpl*nals*(nals+1)/2 && call->nPLs!=nsmpl*nals )  // a mixture of diploid and haploid or haploid only
        error("Wrong number of PL fields? nals=%d npl=%d\n", nals,call->nPLs);

    // Convert PLs to probabilities
    int ngts = nals*(nals+1)/2;
    hts_expand(double, call->nPLs, call->npdg, call->pdg);
    set_pdg(call->pl2p, call->PLs, call->pdg, nsmpl, ngts, unseen);

    #if QS_FROM_PDG
        estimate_qsum(call, rec);
    #else
        // Get sum of qualities, serves as an AF estimate, f_x = QS/N in Eq. 1 in call-m math notes.
        int nqs = bcf_get_info_float(call->hdr, rec, "QS", &call->qsum, &call->nqsum);
        if ( nqs<=0 ) error("The QS annotation not present at %s:%d\n", bcf_seqname(call->hdr,rec),rec->pos+1);
        if ( nqs < nals )
        {
            // Some of the listed alleles do not have the corresponding QS field. This is
            // typically ref-only site with X in ALT.

            hts_expand(float,nals,call->nqsum,call->qsum);
            for (i=nqs; i<nals; i++) call->qsum[i] = 0;
        }

        // If available, take into account reference panel AFs
        if ( call->prior_AN && bcf_get_info_int32(call->hdr, rec, call->prior_AN ,&call->ac, &call->nac)==1 )
        {
            int an = call->ac[0];
            if ( bcf_get_info_int32(call->hdr, rec, call->prior_AC ,&call->ac, &call->nac)==nals-1 )
            {
                int ac0 = an;   // number of alleles in the reference population
                for (i=0; i<nals-1; i++)
                {
                    if ( call->ac[i]==bcf_int32_vector_end ) break;
                    if ( call->ac[i]==bcf_int32_missing ) continue;
                    ac0 -= call->ac[i];
                    call->qsum[i+1] += call->ac[i]*0.5;
                }
                if ( ac0<0 ) error("Incorrect %s,%s values at %s:%d\n", call->prior_AN,call->prior_AC,bcf_seqname(call->hdr,rec),rec->pos+1);
                call->qsum[0] += ac0*0.5;
                for (i=0; i<nals; i++) call->qsum[i] /= nsmpl + 0.5*an;
            }
        }

        float qsum_tot = 0;
        for (i=0; i<nals; i++) qsum_tot += call->qsum[i];

        // Is this still necessary??
        //
        //  if (0&& !call->qsum[0] )
        //  {
        //      // As P(RR)!=0 even for QS(ref)=0, we set QS(ref) to a small value,
        //      // an equivalent of a single reference read.
        //      if ( bcf_get_info_int32(call->hdr, rec, "DP", &call->itmp, &call->n_itmp)!=1 )
        //          error("Could not read DP at %s:%d\n", call->hdr->id[BCF_DT_CTG][rec->rid].key,rec->pos+1);
        //      if ( call->itmp[0] )
        //      {
        //          call->qsum[0] = 1.0 / call->itmp[0] / nsmpl;
        //          qsum_tot += call->qsum[0];
        //      }
        //  }

        if ( qsum_tot ) for (i=0; i<nals; i++) call->qsum[i] /= qsum_tot;
    #endif

    bcf_update_info_int32(call->hdr, rec, "QS", NULL, 0);      // remove QS tag

    // Find the best combination of alleles
    int out_als, nout;
    if ( nals > 8*sizeof(out_als) )
    { 
        fprintf(stderr,"Too many alleles at %s:%d, skipping.\n", bcf_seqname(call->hdr,rec),rec->pos+1); 
        return 0; 
    }
    nout = mcall_find_best_alleles(call, nals, &out_als);

    // Make sure the REF allele is always present
    if ( !(out_als&1) )
    {
        out_als |= 1;
        nout++;
    }
    int is_variant = out_als==1 ? 0 : 1;
    if ( call->flag & CALL_VARONLY && !is_variant ) return 0;

    // With -A, keep all ALTs except X
    if ( call->flag & CALL_KEEPALT )
    {
        nout = 0;
        for (i=0; i<nals; i++)
        {
            if ( i>0 && i==unseen ) continue;
            out_als |= 1<<i;
            nout++;
        }
    }

    int nAC = 0;
    if ( out_als==1 )   // only REF allele on output
    {
        init_allele_trimming_maps(call, 1, nals);
        mcall_set_ref_genotypes(call,nals);
        bcf_update_format_int32(call->hdr, rec, "PL", NULL, 0);    // remove PL, useless now
    }
    else
    {
        // The most likely set of alleles includes non-reference allele (or was enforced), call genotypes.
        // Note that it is a valid outcome if the called genotypes exclude some of the ALTs.
        init_allele_trimming_maps(call, out_als, nals);
        if ( !is_variant )
            mcall_set_ref_genotypes(call,nals);     // running with -A, prevent mcall_call_genotypes from putting some ALT back
        else if ( call->flag & CALL_CONSTR_TRIO )
        {
            if ( nout>4 ) 
            { 
                fprintf(stderr,"Too many alleles at %s:%d, skipping.\n", bcf_seqname(call->hdr,rec),rec->pos+1); 
                return 0; 
            }
            mcall_call_trio_genotypes(call, rec, nals,nout,out_als);
        }
        else
            mcall_call_genotypes(call,rec,nals,nout,out_als);

        // Skip the site if all samples are 0/0. This can happen occasionally.
        nAC = 0;
        for (i=1; i<nout; i++) nAC += call->ac[i];
        if ( !nAC && call->flag & CALL_VARONLY ) return 0;
        mcall_trim_PLs(call, rec, nals, nout, out_als);
    }
    if ( nals!=nout ) mcall_trim_numberR(call, rec, nals, nout, out_als);

    // Set QUAL and calculate HWE-related annotations
    if ( nAC )
    {
        float icb = calc_ICB(call->ac[0],nAC, call->nhets, call->ndiploid);
        if ( icb != HUGE_VAL ) bcf_update_info_float(call->hdr, rec, "ICB", &icb, 1);

        float hob = calc_HOB(call->ac[0],nAC, call->nhets, call->ndiploid);
        if ( hob != HUGE_VAL ) bcf_update_info_float(call->hdr, rec, "HOB", &hob, 1);

        // Quality of a variant site. fabs() to avoid negative zeros in VCF output when CALL_KEEPALT is set
        rec->qual = -4.343*(call->ref_lk - logsumexp2(call->lk_sum,call->ref_lk));
    }
    else
    {
        // Set the quality of a REF site
        if ( call->lk_sum==-HUGE_VAL )  // no support from (high quality) reads, so QUAL=1-prior
            rec->qual = call->theta ? -4.343*call->theta : 0;
        else
            rec->qual = -4.343*(call->lk_sum - logsumexp2(call->lk_sum,call->ref_lk));
    }

    if ( rec->qual>999 ) rec->qual = 999;
    if ( rec->qual>50 ) rec->qual = rint(rec->qual);

    // AC, AN
    if ( nout>1 ) bcf_update_info_int32(call->hdr, rec, "AC", call->ac+1, nout-1);
    nAC += call->ac[0];
    bcf_update_info_int32(call->hdr, rec, "AN", &nAC, 1);

    // Remove unused alleles
    hts_expand(char*,nout,call->nals,call->als);
    for (i=0; i<nals; i++)
        if ( call->als_map[i]>=0 ) call->als[call->als_map[i]] = rec->d.allele[i];
    bcf_update_alleles(call->hdr, rec, (const char**)call->als, nout);
    bcf_update_genotypes(call->hdr, rec, call->gts, nsmpl*2);

    // DP4 tag
    if ( bcf_get_info_float(call->hdr, rec, "I16", &call->anno16, &call->n16)==16 )
    {
        int32_t dp[4]; dp[0] = call->anno16[0]; dp[1] = call->anno16[1]; dp[2] = call->anno16[2]; dp[3] = call->anno16[3];
        bcf_update_info_int32(call->hdr, rec, "DP4", dp, 4);

        int32_t mq = (call->anno16[8]+call->anno16[10])/(call->anno16[0]+call->anno16[1]+call->anno16[2]+call->anno16[3]);
        bcf_update_info_int32(call->hdr, rec, "MQ", &mq, 1);
    }

    bcf_update_info_int32(call->hdr, rec, "I16", NULL, 0);     // remove I16 tag

    return nout;
}

