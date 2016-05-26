/*  vcfroh.c -- HMM model for detecting runs of autozygosity.

    Copyright (C) 2013-2015 Genome Research Ltd.

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

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include "bcftools.h"
#include "HMM.h"

#define STATE_HW 0        // normal state, follows Hardy-Weinberg allele frequencies
#define STATE_AZ 1        // autozygous state

/** Genetic map */
typedef struct
{
    int pos;
    double rate;
}
genmap_t;

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr;
    double t2AZ, t2HW;      // P(AZ|HW) and P(HW|AZ) parameters
    double unseen_PL, dflt_AF;

    char *genmap_fname;
    genmap_t *genmap;
    int ngenmap, mgenmap, igenmap;
    double rec_rate;        // constant recombination rate if > 0

    hmm_t *hmm;
    double *eprob;          // emission probs [2*nsites,msites]
    uint32_t *sites;        // positions [nsites,msites]
    int nsites, msites;
    int nrids, *rids, *rid_offs;    // multiple chroms with vi_training

    int32_t *itmp;
    int nitmp, mitmp;
    float *AFs;
    int mAFs;

    double pl2p[256], *pdg;
    int32_t skip_rid, prev_rid, prev_pos;

    int ntot, nused;            // some stats to detect if things didn't go awfully wrong
    int ismpl, nsmpl;           // index of query sample
    char *estimate_AF, *sample; // list of samples for AF estimate and query sample
    char **argv, *targets_list, *regions_list, *af_fname, *af_tag;
    int argc, fake_PLs, snps_only, vi_training;
}
args_t;

void set_tprob_genmap(hmm_t *hmm, uint32_t prev_pos, uint32_t pos, void *data, double *tprob);
void set_tprob_recrate(hmm_t *hmm, uint32_t prev_pos, uint32_t pos, void *data, double *tprob);

void *smalloc(size_t size)
{
    void *mem = malloc(size);
    if ( !mem ) error("malloc: Could not allocate %d bytes\n", (int)size);
    return mem;
}

static void init_data(args_t *args)
{
    args->prev_rid = args->skip_rid = -1;
    args->hdr = args->files->readers[0].header;

    if ( !args->sample )
    {
        if ( bcf_hdr_nsamples(args->hdr)>1 ) error("Missing the option -s, --sample\n");
        args->sample = strdup(args->hdr->samples[0]);
    }
    if ( !bcf_hdr_nsamples(args->hdr) ) error("No samples in the VCF?\n");

    // Set samples
    kstring_t str = {0,0,0};
    if ( args->estimate_AF && strcmp("-",args->estimate_AF) )
    {
        int i, n;
        char **smpls = hts_readlist(args->estimate_AF, 1, &n);

        // Make sure the query sample is included
        for (i=0; i<n; i++)
            if ( !strcmp(args->sample,smpls[i]) ) break;

        // Add the query sample if not present
        if ( i!=n ) kputs(args->sample, &str);

        for (i=0; i<n; i++)
        {
            if ( str.l ) kputc(',', &str);
            kputs(smpls[i], &str);
            free(smpls[i]);
        }
        free(smpls);
    }
    else if ( !args->estimate_AF )
        kputs(args->sample, &str);

    if ( str.l )
    {
        int ret = bcf_hdr_set_samples(args->hdr, str.s, 0);
        if ( ret<0 ) error("Error parsing the list of samples: %s\n", str.s);
        else if ( ret>0 ) error("The %d-th sample not found in the VCF\n", ret);
    }

    if ( args->af_tag )
        if ( !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_INFO,bcf_hdr_id2int(args->hdr,BCF_DT_ID,args->af_tag)) )
            error("No such INFO tag in the VCF: %s\n", args->af_tag);

    args->nsmpl = bcf_hdr_nsamples(args->hdr);
    args->ismpl = bcf_hdr_id2int(args->hdr, BCF_DT_SAMPLE, args->sample);
    free(str.s);

    int i;
    for (i=0; i<256; i++) args->pl2p[i] = pow(10., -i/10.);

    // Init transition matrix and HMM
    double tprob[4];
    MAT(tprob,2,STATE_HW,STATE_HW) = 1 - args->t2AZ;
    MAT(tprob,2,STATE_HW,STATE_AZ) = args->t2HW;
    MAT(tprob,2,STATE_AZ,STATE_HW) = args->t2AZ;
    MAT(tprob,2,STATE_AZ,STATE_AZ) = 1 - args->t2HW; 

    if ( args->genmap_fname ) 
    {
        args->hmm = hmm_init(2, tprob, 0);
        hmm_set_tprob_func(args->hmm, set_tprob_genmap, args);
    }
    else if ( args->rec_rate > 0 )
    {
        args->hmm = hmm_init(2, tprob, 0);
        hmm_set_tprob_func(args->hmm, set_tprob_recrate, args);

    }
    else
        args->hmm = hmm_init(2, tprob, 10000);

    // print header
    printf("# This file was produced by: bcftools roh(%s+htslib-%s)\n", bcftools_version(),hts_version());
    printf("# The command line was:\tbcftools %s", args->argv[0]);
    for (i=1; i<args->argc; i++)
        printf(" %s",args->argv[i]);
    printf("\n#\n");
    printf("# [1]Chromosome\t[2]Position\t[3]State (0:HW, 1:AZ)\t[4]Quality\n");
}

static void destroy_data(args_t *args)
{
    free(args->sites);
    free(args->eprob);
    free(args->sample);
    free(args->rids);
    free(args->rid_offs);
    hmm_destroy(args->hmm);
    bcf_sr_destroy(args->files);
    free(args->itmp); free(args->AFs); free(args->pdg);
    free(args->genmap);
}

static int load_genmap(args_t *args, bcf1_t *line)
{
    if ( !args->genmap_fname ) { args->ngenmap = 0; return 0; }

    kstring_t str = {0,0,0};
    char *fname = strstr(args->genmap_fname,"{CHROM}");
    if ( fname )
    {
        kputsn(args->genmap_fname, fname - args->genmap_fname, &str);
        kputs(bcf_seqname(args->hdr,line), &str);
        kputs(fname+7,&str);
        fname = str.s;
    }
    else
        fname = args->genmap_fname;

    htsFile *fp = hts_open(fname, "rb");
    if ( !fp )
    {
        args->ngenmap = 0;
        return -1;
    }

    hts_getline(fp, KS_SEP_LINE, &str);
    if ( strcmp(str.s,"position COMBINED_rate(cM/Mb) Genetic_Map(cM)") )
        error("Unexpected header, found:\n\t[%s], but expected:\n\t[position COMBINED_rate(cM/Mb) Genetic_Map(cM)]\n", fname, str.s);

    args->ngenmap = args->igenmap = 0;
    while ( hts_getline(fp, KS_SEP_LINE, &str) > 0 )
    {
        args->ngenmap++;
        hts_expand(genmap_t,args->ngenmap,args->mgenmap,args->genmap);
        genmap_t *gm = &args->genmap[args->ngenmap-1];

        char *tmp, *end;
        gm->pos = strtol(str.s, &tmp, 10);
        if ( str.s==tmp ) error("Could not parse %s: %s\n", fname, str.s);

        // skip second column
        tmp++;
        while ( *tmp && !isspace(*tmp) ) tmp++;

        // read the genetic map in cM
        gm->rate = strtod(tmp+1, &end);
        if ( tmp+1==end ) error("Could not parse %s: %s\n", fname, str.s);
    }
    if ( !args->ngenmap ) error("Genetic map empty?\n");
    int i;
    for (i=0; i<args->ngenmap; i++) args->genmap[i].rate /= args->genmap[args->ngenmap-1].rate; // scale to 1
    if ( hts_close(fp) ) error("Close failed\n");
    free(str.s);
    return 0;
}

static double get_genmap_rate(args_t *args, int start, int end)
{
    // position i to be equal to or smaller than start
    int i = args->igenmap;
    if ( args->genmap[i].pos > start )
    {
        while ( i>0 && args->genmap[i].pos > start ) i--;
    }
    else
    {
        while ( i+1<args->ngenmap && args->genmap[i+1].pos < start ) i++;
    }
    // position j to be equal or larger than end
    int j = i;
    while ( j+1<args->ngenmap && args->genmap[j].pos < end ) j++;

    if ( i==j )
    {
        args->igenmap = i;
        return 0;
    }

    if ( start <  args->genmap[i].pos ) start = args->genmap[i].pos;
    if ( end >  args->genmap[j].pos ) end = args->genmap[j].pos;
    double rate = (args->genmap[j].rate - args->genmap[i].rate)/(args->genmap[j].pos - args->genmap[i].pos) * (end-start);
    args->igenmap = j;
    return rate;
}

void set_tprob_genmap(hmm_t *hmm, uint32_t prev_pos, uint32_t pos, void *data, double *tprob)
{
    args_t *args = (args_t*) data;
    double ci = get_genmap_rate(args, pos - prev_pos, pos);
    MAT(tprob,2,STATE_HW,STATE_AZ) *= ci;
    MAT(tprob,2,STATE_AZ,STATE_HW) *= ci;
    MAT(tprob,2,STATE_AZ,STATE_AZ)  = 1 - MAT(tprob,2,STATE_HW,STATE_AZ);
    MAT(tprob,2,STATE_HW,STATE_HW)  = 1 - MAT(tprob,2,STATE_AZ,STATE_HW);
}

void set_tprob_recrate(hmm_t *hmm, uint32_t prev_pos, uint32_t pos, void *data, double *tprob)
{
    args_t *args = (args_t*) data;
    double ci = (pos - prev_pos) * args->rec_rate;
    MAT(tprob,2,STATE_HW,STATE_AZ) *= ci;
    MAT(tprob,2,STATE_AZ,STATE_HW) *= ci;
    MAT(tprob,2,STATE_AZ,STATE_AZ)  = 1 - MAT(tprob,2,STATE_HW,STATE_AZ);
    MAT(tprob,2,STATE_HW,STATE_HW)  = 1 - MAT(tprob,2,STATE_AZ,STATE_HW);
}


/**
 *  This function implements the HMM model:
 *    D = Data, AZ = autozygosity, HW = Hardy-Weinberg (non-autozygosity),
 *    f = non-ref allele frequency
 *
 *  Emission probabilities:
 *    oAZ = P_i(D|AZ) = (1-f)*P(D|RR) + f*P(D|AA)
 *    oHW = P_i(D|HW) = (1-f)^2 * P(D|RR) + f^2 * P(D|AA) + 2*f*(1-f)*P(D|RA)
 *
 *  Transition probabilities:
 *    tAZ = P(AZ|HW)  .. parameter
 *    tHW = P(HW|AZ)  .. parameter
 *
 *    ci  = P_i(C)    .. probability of cross-over at site i, from genetic map
 *
 *    AZi = P_i(AZ)   .. probability of site i being AZ/non-AZ, scaled so that AZi+HWi = 1
 *    HWi = P_i(HW)
 *
 *    P_i(AZ|HW) = P(AZ|HW) * ci * HW{i-1}     = tAZ * ci * (1 - AZ{i-1})
 *    P_i(HW|AZ) = P(HW|AZ) * ci * AZ{i-1}     = tHW * ci * AZ{i-1}
 *    P_i(AZ|AZ) = 1 - P_i(HW|AZ)
 *    P_i(HW|HW) = 1 - P_i(AZ|HW)
 *
 */

static void flush_viterbi(args_t *args)
{
    int i,j;

    if ( !args->nsites ) return; 

    if ( !args->vi_training )
    {
        // single viterbi pass, one chromsome
        hmm_run_viterbi(args->hmm, args->nsites, args->eprob, args->sites);
        hmm_run_fwd_bwd(args->hmm, args->nsites, args->eprob, args->sites);
        double *fwd = hmm_get_fwd_bwd_prob(args->hmm);

        const char *chr = bcf_hdr_id2name(args->hdr,args->prev_rid);
        uint8_t *vpath = hmm_get_viterbi_path(args->hmm);
        for (i=0; i<args->nsites; i++)
        {
            int state = vpath[i*2]==STATE_AZ ? 1 : 0;
            double *pval = fwd + i*2;
            printf("%s\t%d\t%d\t%.1f\n", chr,args->sites[i]+1, state, phred_score(1.0-pval[state]));
        }
        return;
    }

    // viterbi training, multiple chromosomes
    double t2az_prev, t2hw_prev;
    double deltaz, delthw;
    int niter = 0;
    do
    {
        double *tprob_arr = hmm_get_tprob(args->hmm);
        t2az_prev = MAT(tprob_arr,2,1,0); //args->t2AZ;
        t2hw_prev = MAT(tprob_arr,2,0,1); //args->t2HW;
        double tcounts[] = { 0,0,0,0 };
        for (i=0; i<args->nrids; i++)
        {
            // run viterbi for each chromosomes. eprob and sites contain
            // multiple chromosomes, rid_offs mark the boundaries
            int ioff = args->rid_offs[i];
            int nsites = (i+1==args->nrids ? args->nsites : args->rid_offs[i+1]) - ioff;
            hmm_run_viterbi(args->hmm, nsites, args->eprob+ioff*2, args->sites+ioff);

            // what transitions were observed: add to the total counts
            uint8_t *vpath = hmm_get_viterbi_path(args->hmm);
            for (j=1; j<nsites; j++)
            {
                // count the number of transitions
                int prev_state = vpath[2*(j-1)];
                int curr_state = vpath[2*j];
                MAT(tcounts,2,curr_state,prev_state) += 1;
            }
        }

        // update the transition matrix
        int n = 1;
        for (i=0; i<2; i++)
        {
            for (j=0; j<2; j++) n += MAT(tcounts,2,i,j);
        }
        for (i=0; i<2; i++)
        {
            for (j=0; j<2; j++)
            {
                // no transition to i-th state was observed, set to a small number
                if ( !MAT(tcounts,2,i,j) ) MAT(tcounts,2,i,j) = 0.1/n;
                else MAT(tcounts,2,i,j) /= n;
            }
        }

        // normalize
        for (i=0; i<2; i++)
        {
            double norm = 0;
            for (j=0; j<2; j++) norm += MAT(tcounts,2,j,i);
            assert( norm!=0 );
            for (j=0; j<2; j++) MAT(tcounts,2,j,i) /= norm;
        }

        if ( args->genmap_fname || args->rec_rate > 0 )
            hmm_set_tprob(args->hmm, tcounts, 0);
        else
            hmm_set_tprob(args->hmm, tcounts, 10000);

        tprob_arr = hmm_get_tprob(args->hmm);
        deltaz = fabs(MAT(tprob_arr,2,1,0)-t2az_prev);
        delthw = fabs(MAT(tprob_arr,2,0,1)-t2hw_prev);
        niter++;
        fprintf(stderr,"Viterbi training, iteration %d: dAZ=%e dHW=%e\tP(HW|HW)=%e  P(AZ|HW)=%e  P(AZ|AZ)=%e  P(HW|AZ)=%e\n", 
            niter,deltaz,delthw,
            MAT(tprob_arr,2,STATE_HW,STATE_HW),MAT(tprob_arr,2,STATE_AZ,STATE_HW),
            MAT(tprob_arr,2,STATE_AZ,STATE_AZ),MAT(tprob_arr,2,STATE_HW,STATE_AZ));
    }
    while ( deltaz > 0.0 || delthw > 0.0 );
    double *tprob_arr = hmm_get_tprob(args->hmm);
    fprintf(stderr, "Viterbi training converged in %d iterations to P(HW|HW)=%e  P(AZ|HW)=%e  P(AZ|AZ)=%e  P(HW|AZ)=%e\n", niter,
            MAT(tprob_arr,2,STATE_HW,STATE_HW),MAT(tprob_arr,2,STATE_AZ,STATE_HW),
            MAT(tprob_arr,2,STATE_AZ,STATE_AZ),MAT(tprob_arr,2,STATE_HW,STATE_AZ));
    
    // output the results
    for (i=0; i<args->nrids; i++)
    {
        int ioff = args->rid_offs[i];
        int nsites = (i+1==args->nrids ? args->nsites : args->rid_offs[i+1]) - ioff;
        hmm_run_viterbi(args->hmm, nsites, args->eprob+ioff*2, args->sites+ioff);
        hmm_run_fwd_bwd(args->hmm, nsites, args->eprob+ioff*2, args->sites+ioff);
        uint8_t *vpath = hmm_get_viterbi_path(args->hmm);
        double  *fwd   = hmm_get_fwd_bwd_prob(args->hmm);

        const char *chr = bcf_hdr_id2name(args->hdr,args->rids[i]);
        for (j=0; j<nsites; j++)
        {
            int state = vpath[j*2];
            double pval = fwd[j*2 + state];
            printf("%s\t%d\t%d\t%e\n", chr,args->sites[ioff+j]+1,state==STATE_AZ ? 1 : 0, pval);
        }
    }
}

static void push_rid(args_t *args, int rid)
{
    args->nrids++;
    args->rids = (int*) realloc(args->rids, args->nrids*sizeof(int));
    args->rid_offs = (int*) realloc(args->rid_offs, args->nrids*sizeof(int));
    args->rids[ args->nrids-1 ] = rid;
    args->rid_offs[ args->nrids-1 ] = args->nsites;
}

int read_AF(bcf_sr_regions_t *tgt, bcf1_t *line, double *alt_freq)
{
    if ( tgt->nals != line->n_allele ) return -1;    // number of alleles does not match

    int i;
    for (i=0; i<tgt->nals; i++)
        if ( strcmp(line->d.allele[i],tgt->als[i]) ) break; // we could be smarter, see vcmp
    if ( i<tgt->nals ) return -1;

    char *tmp, *str = tgt->line.s;
    i = 0;
    while ( *str && i<3 ) 
    {
        if ( *str=='\t' ) i++;
        str++;
    }
    *alt_freq = strtod(str, &tmp);
    if ( *tmp && !isspace(*tmp) )
    {
        if ( str[0]=='.' && (!str[1] || isspace(str[1])) ) return -1; // missing value
        error("Could not parse: [%s]\n", tgt->line.s);
    }
    if ( *alt_freq<0 || *alt_freq>1 ) error("Could not parse AF: [%s]\n", tgt->line.s);
    return 0;
}

int estimate_AF(args_t *args, bcf1_t *line, double *alt_freq)
{
    if ( !args->nitmp )
    {
        args->nitmp = bcf_get_genotypes(args->hdr, line, &args->itmp, &args->mitmp);
        if ( args->nitmp != 2*args->nsmpl ) return -1;     // not diploid?
        args->nitmp /= args->nsmpl;
    }

    int i, nalt = 0, nref = 0;
    for (i=0; i<args->nsmpl; i++)
    {
        int32_t *gt = &args->itmp[i*args->nitmp];

        if ( bcf_gt_is_missing(gt[0]) || bcf_gt_is_missing(gt[1]) ) continue;

        if ( bcf_gt_allele(gt[0]) ) nalt++;
        else nref++;

        if ( bcf_gt_allele(gt[1]) ) nalt++;
        else nref++;
    }
    if ( !nalt && !nref ) return -1;

    *alt_freq = (double)nalt / (nalt + nref);
    return 0;
}


int parse_line(args_t *args, bcf1_t *line, double *alt_freq, double *pdg)
{
    args->nitmp = 0;

    // Set allele frequency
    int ret;
    if ( args->af_tag )
    {
        // Use an INFO tag provided by the user
        ret = bcf_get_info_float(args->hdr, line, args->af_tag, &args->AFs, &args->mAFs);
        if ( ret==1 )
            *alt_freq = args->AFs[0];
        if ( ret==-2 )
            error("Type mismatch for INFO/%s tag at %s:%d\n", args->af_tag, bcf_seqname(args->hdr,line), line->pos+1);
    }
    else if ( args->af_fname ) 
    {
        // Read AF from a file
        ret = read_AF(args->files->targets, line, alt_freq);
    }
    else
    {
        // Use GTs or AC/AN: GTs when AC/AN not present or when GTs explicitly requested by --estimate-AF
        ret = -1;
        if ( !args->estimate_AF )
        {
            int AC = -1, AN = 0;
            ret = bcf_get_info_int32(args->hdr, line, "AN", &args->itmp, &args->mitmp);
            if ( ret==1 )
            {
                AN = args->itmp[0];
                ret = bcf_get_info_int32(args->hdr, line, "AC", &args->itmp, &args->mitmp);
                if ( ret>0 )
                    AC = args->itmp[0];
            }
            if ( AN<=0 || AC<0 ) 
                ret = -1;
            else 
                *alt_freq = (double) AC/AN;
        }
        if ( ret==-1 )
            ret = estimate_AF(args, line, alt_freq);    // reads GTs into args->itmp
    }

    if ( ret<0 ) return ret;
    if ( *alt_freq==0.0 )
    {
        if ( args->dflt_AF==0 ) return -1;       // we skip sites with AF=0
        *alt_freq = args->dflt_AF;
    }

    // Set P(D|G)
    if ( args->fake_PLs )
    {
        if ( !args->nitmp )
        {
            args->nitmp = bcf_get_genotypes(args->hdr, line, &args->itmp, &args->mitmp);
            if ( args->nitmp != 2*args->nsmpl ) return -1;     // not diploid?
            args->nitmp /= args->nsmpl;
        }

        int32_t *gt = &args->itmp[args->ismpl*args->nitmp];
        if ( bcf_gt_is_missing(gt[0]) || bcf_gt_is_missing(gt[1]) ) return -1;

        int a = bcf_gt_allele(gt[0]);
        int b = bcf_gt_allele(gt[1]);
        if ( a!=b )
        {
            pdg[0] = pdg[2] = args->unseen_PL;
            pdg[1] = 1 - 2*args->unseen_PL;
        }
        else if ( a==0 )
        {
            pdg[0] = 1 - 2*args->unseen_PL;
            pdg[1] = pdg[2] = args->unseen_PL;
        }
        else
        {
            pdg[0] = pdg[1] = args->unseen_PL;
            pdg[2] = 1 - 2*args->unseen_PL;
        }
    }
    else
    {
        args->nitmp = bcf_get_format_int32(args->hdr, line, "PL", &args->itmp, &args->mitmp);
        if ( args->nitmp != args->nsmpl*line->n_allele*(line->n_allele+1)/2. ) return -1;     // not diploid?
        args->nitmp /= args->nsmpl;

        int32_t *pl = &args->itmp[args->ismpl*args->nitmp];
        pdg[0] = pl[0] < 256 ? args->pl2p[ pl[0] ] : 1.0;
        pdg[1] = pl[1] < 256 ? args->pl2p[ pl[1] ] : 1.0;
        pdg[2] = pl[2] < 256 ? args->pl2p[ pl[2] ] : 1.0;

        double sum = pdg[0] + pdg[1] + pdg[2];
        if ( !sum ) return -1;
        pdg[0] /= sum;
        pdg[1] /= sum;
        pdg[2] /= sum;
    }

    return 0;
}

static void vcfroh(args_t *args, bcf1_t *line)
{
    // Are we done?
    if ( !line )
    { 
        flush_viterbi(args);
        return; 
    }
    args->ntot++;

    // Skip unwanted lines
    if ( line->rid == args->skip_rid ) return;
    if ( line->n_allele==1 ) return;    // no ALT allele
    if ( line->n_allele!=2 ) return;    // only biallelic sites
    if ( args->snps_only && !bcf_is_snp(line) ) return;

    // Initialize genetic map
    int skip_rid = 0;
    if ( args->prev_rid<0 )
    {
        args->prev_rid = line->rid;
        args->prev_pos = line->pos;
        skip_rid = load_genmap(args, line);
        if ( !skip_rid && args->vi_training ) push_rid(args, line->rid);
    }

    // New chromosome?
    if ( args->prev_rid!=line->rid )
    {
        skip_rid = load_genmap(args, line);
        if ( args->vi_training )
        {
            if ( !skip_rid ) push_rid(args, line->rid);
        }
        else
        {
            flush_viterbi(args);
            args->nsites = 0;
        }
        args->prev_rid = line->rid;
        args->prev_pos = line->pos;
    }

    if ( skip_rid )
    {
        fprintf(stderr,"Skipping the sequence, no genmap for %s\n", bcf_seqname(args->hdr,line));
        args->skip_rid = line->rid;
        return;
    }
    if ( args->prev_pos > line->pos ) error("The file is not sorted?!\n");

    args->prev_rid = line->rid;
    args->prev_pos = line->pos;


    // Ready for the new site
    int m = args->msites;
    hts_expand(uint32_t,args->nsites+1,args->msites,args->sites);
    if ( args->msites!=m )
        args->eprob = (double*) realloc(args->eprob,sizeof(double)*args->msites*2);

    // Set likelihoods and alternate allele frequencies
    double alt_freq, pdg[3];
    if ( parse_line(args, line, &alt_freq, pdg)<0 ) return; // something went wrong

    args->nused++;

    // Calculate emission probabilities P(D|AZ) and P(D|HW)
    double *eprob = &args->eprob[2*args->nsites];
    eprob[STATE_AZ] = pdg[0]*(1-alt_freq) + pdg[2]*alt_freq;
    eprob[STATE_HW] = pdg[0]*(1-alt_freq)*(1-alt_freq) + 2*pdg[1]*(1-alt_freq)*alt_freq + pdg[2]*alt_freq*alt_freq;

    args->sites[args->nsites] = line->pos;
    args->nsites++;
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   HMM model for detecting runs of autozygosity.\n");
    fprintf(stderr, "Usage:   bcftools roh [options] <in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "General Options:\n");
    fprintf(stderr, "        --AF-dflt <float>              if AF is not known, use this allele frequency [skip]\n");
    fprintf(stderr, "        --AF-tag <TAG>                 use TAG for allele frequency\n");
    fprintf(stderr, "        --AF-file <file>               read allele frequencies from file (CHR\\tPOS\\tREF,ALT\\tAF)\n");
    fprintf(stderr, "    -e, --estimate-AF <file>           calculate AC,AN counts on the fly, using either all samples (\"-\") or samples listed in <file>\n");
    fprintf(stderr, "    -G, --GTs-only <float>             use GTs, ignore PLs, use <float> for PL of unseen genotypes. Safe value to use is 30 to account for GT errors.\n");
    fprintf(stderr, "    -I, --skip-indels                  skip indels as their genotypes are enriched for errors\n");
    fprintf(stderr, "    -m, --genetic-map <file>           genetic map in IMPUTE2 format, single file or mask, where string \"{CHROM}\" is replaced with chromosome name\n");
    fprintf(stderr, "    -M, --rec-rate <float>             constant recombination rate per bp\n");
    fprintf(stderr, "    -r, --regions <region>             restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>          restrict to regions listed in a file\n");
    fprintf(stderr, "    -s, --sample <sample>              sample to analyze\n");
    fprintf(stderr, "    -t, --targets <region>             similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file <file>          similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "HMM Options:\n");
    fprintf(stderr, "    -a, --hw-to-az <float>             P(AZ|HW) transition probability from HW (Hardy-Weinberg) to AZ (autozygous) state [6.7e-8]\n");
    fprintf(stderr, "    -H, --az-to-hw <float>             P(HW|AZ) transition probability from AZ to HW state [5e-9]\n");
    fprintf(stderr, "    -V, --viterbi-training             perform Viterbi training to estimate transition probabilities\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfroh(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->t2AZ    = 6.7e-8;
    args->t2HW    = 5e-9;
    args->rec_rate = 0;
    int regions_is_file = 0, targets_is_file = 0;

    static struct option loptions[] =
    {
        {"AF-tag",1,0,0},
        {"AF-file",1,0,1},
        {"AF-dflt",1,0,2},
        {"estimate-AF",1,0,'e'},
        {"GTs-only",1,0,'G'},
        {"sample",1,0,'s'},
        {"hw-to-az",1,0,'a'},
        {"az-to-hw",1,0,'H'},
        {"viterbi-training",0,0,'V'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"genetic-map",1,0,'m'},
        {"rec-rate",1,0,'M'},
        {"skip-indels",0,0,'I'},
        {0,0,0,0}
    };

    int naf_opts = 0;
    char *tmp;
    while ((c = getopt_long(argc, argv, "h?r:R:t:T:H:a:s:m:M:G:Ia:e:V",loptions,NULL)) >= 0) {
        switch (c) {
            case 0: args->af_tag = optarg; naf_opts++; break;
            case 1: args->af_fname = optarg; naf_opts++; break;
            case 2: 
                args->dflt_AF = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: --AF-dflt %s\n", optarg);
                break;
            case 'e': args->estimate_AF = optarg; naf_opts++; break;
            case 'I': args->snps_only = 1; break;
            case 'G':
                args->fake_PLs = 1; 
                args->unseen_PL = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -G %s\n", optarg);
                args->unseen_PL = pow(10,-args->unseen_PL/10.); 
                break;
            case 'm': args->genmap_fname = optarg; break;
            case 'M':
                args->rec_rate = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -M %s\n", optarg);
                break;
            case 's': args->sample = strdup(optarg); break;
            case 'a':
                args->t2AZ = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -a %s\n", optarg);
                break;
            case 'H':
                args->t2HW = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -H %s\n", optarg);
                break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 'V': args->vi_training = 1; break;
            case 'h': 
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( argc<optind+1 ) usage(args);
    if ( args->t2AZ<0 || args->t2AZ>1 ) error("Error: The parameter --hw-to-az is not in [0,1]\n", args->t2AZ);
    if ( args->t2HW<0 || args->t2HW>1 ) error("Error: The parameter --az-to-hw is not in [0,1]\n", args->t2HW);
    if ( naf_opts>1 ) error("Error: The options --AF-tag, --AF-file and -e are mutually exclusive\n");
    if ( args->af_fname && args->targets_list ) error("Error: The options --AF-file and -t are mutually exclusive\n");
    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->targets_list )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }
    if ( args->af_fname )
    {
        if ( bcf_sr_set_targets(args->files, args->af_fname, 1, 3)<0 )
            error("Failed to read the targets: %s\n", args->af_fname);
    }
    if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open %s: %s\n", argv[optind],bcf_sr_strerror(args->files->errnum));

    init_data(args);
    while ( bcf_sr_next_line(args->files) )
    {
        vcfroh(args, args->files->readers[0].buffer[0]);
    }
    vcfroh(args, NULL);
    fprintf(stderr,"Number of lines: total/processed: %d/%d\n", args->ntot,args->nused);
    destroy_data(args);
    free(args);
    return 0;
}


