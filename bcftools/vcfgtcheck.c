/*  vcfgtcheck.c -- Check sample identity.

    Copyright (C) 2013-2014 Genome Research Ltd.

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
#include <stdarg.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"

typedef struct
{
    bcf_srs_t *files;           // first reader is the query VCF - single sample normally or multi-sample for cross-check
    bcf_hdr_t *gt_hdr, *sm_hdr; // VCF with genotypes to compare against and the query VCF
    int ntmp_arr, npl_arr;
    int32_t *tmp_arr, *pl_arr;
    double *lks, *sites;
    int *cnts, *dps, hom_only, cross_check, all_sites;
    char *cwd, **argv, *gt_fname, *plot, *query_sample, *target_sample;
    int argc, no_PLs;
}
args_t;

FILE *open_file(char **fname, const char *mode, const char *fmt, ...);
char *msprintf(const char *fmt, ...);
void mkdir_p(const char *fmt, ...);

void py_plot(char *script)
{
    mkdir_p(script);
    int len = strlen(script);
    char *cmd = !strcmp(".py",script+len-3) ? msprintf("python %s", script) : msprintf("python %s.py", script);
    int ret = system(cmd);
    if ( ret ) fprintf(stderr, "The command returned non-zero status %d: %s\n", ret, cmd);
    free(cmd);
}

static void plot_check(args_t *args, char *target_sample, char *query_sample)
{
    char *fname;
    FILE *fp = open_file(&fname, "w", "%s.py", args->plot);
    fprintf(fp,
            "import matplotlib as mpl\n"
            "mpl.use('Agg')\n"
            "import matplotlib.pyplot as plt\n"
            "import matplotlib.gridspec as gridspec\n"
            "import csv\n"
            "csv.register_dialect('tab', delimiter='\\t', quoting=csv.QUOTE_NONE)\n"
            "\n"
            "sample_ids = False\n"
            "\n"
            "dat = []\n"
            "with open('%s.tab', 'rb') as f:\n"
            "    reader = csv.reader(f, 'tab')\n"
            "    for row in reader:\n"
            "        if row[0][0]=='#': continue\n"
            "        if row[0]!='CN': continue\n"
            "        tgt = 0\n"
            "        if row[4]=='%s': tgt = 1\n"
            "        dat.append([float(row[1]), float(row[2]), float(row[3]), tgt, row[4]])\n"
            "\n"
            "dat = sorted(dat)\n"
            "\n"
            "iq = -1; dp = 0\n"
            "for i in range(len(dat)):\n"
            "    if iq==-1 and dat[i][3]==1: iq = i\n"
            "    dp += dat[i][2]\n"
            "dp /= len(dat)\n"
            "\n"
            "fig,ax1 = plt.subplots(figsize=(8,5))\n"
            "ax2 = ax1.twinx()\n"
            "plots  = ax1.plot([x[0] for x in dat],'o-', ms=3, color='g', mec='g', label='Discordance (total)')\n"
            "plots += ax1.plot([x[1] for x in dat], '^', ms=3, color='r', mec='r', label='Discordance (avg per site)')\n"
            "plots += ax2.plot([x[2] for x in dat],'v', ms=3, color='k', label='Number of sites')\n"
            "if iq!=-1:\n"
            "   ax1.plot([iq],[dat[iq][0]],'o',color='orange', ms=9)\n"
            "   ax1.annotate('%s',xy=(iq,dat[iq][0]), xytext=(5,5), textcoords='offset points',fontsize='xx-small',rotation=45,va='bottom',ha='left')\n"
            "   ax1.plot([iq],[dat[iq][1]],'^',color='red', ms=5)\n"
            "for tl in ax1.get_yticklabels(): tl.set_color('g')\n"
            "for tl in ax2.get_yticklabels(): tl.set_color('k'); tl.set_fontsize(9)\n"
            "min_dp = min([x[2] for x in dat])\n"
            "max_dp = max([x[2] for x in dat])\n"
            "ax2.set_ylim(min_dp-1,max_dp+1)\n"
            "ax1.set_title('Discordance with %s')\n"
            "ax1.set_xlim(-0.05*len(dat),1.05*(len(dat)-1))\n"
            "ax1.set_xlabel('Sample ID')\n"
            "plt.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9)\n"
            "if sample_ids:\n"
            "   ax1.set_xticks(range(len(dat)))\n"
            "   ax1.set_xticklabels([x[4] for x in dat],**{'rotation':45, 'ha':'right', 'fontsize':8})\n"
            "   plt.subplots_adjust(bottom=0.2)\n"
            "ax1.set_ylabel('Discordance',color='g')\n"
            "ax2.set_ylabel('Number of sites',color='k')\n"
            "ax2.ticklabel_format(style='sci', scilimits=(-3,2), axis='y')\n"
            "ax1.ticklabel_format(style='sci', scilimits=(-3,2), axis='y')\n"
            "labels = [l.get_label() for l in plots]\n"
            "plt.legend(plots,labels,numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)\n"
            "plt.savefig('%s.png')\n"
            "plt.close()\n"
            "\n", args->plot, target_sample, target_sample, query_sample, args->plot
           );
    fclose(fp);
    py_plot(fname);
    free(fname);
}

static void plot_cross_check(args_t *args)
{
    char *fname;
    FILE *fp = open_file(&fname, "w", "%s.py", args->plot);
    fprintf(fp,
            "import matplotlib as mpl\n"
            "mpl.use('Agg')\n"
            "import matplotlib.pyplot as plt\n"
            "import matplotlib.gridspec as gridspec\n"
            "import csv\n"
            "csv.register_dialect('tab', delimiter='\\t', quoting=csv.QUOTE_NONE)\n"
            "avg   = []\n"
            "dp    = []\n"
            "sm2id = {}\n"
            "dat   = None\n"
            "min   = None\n"
            "max   = None\n"
            "with open('%s.tab', 'rb') as f:\n"
            "   reader = csv.reader(f, 'tab')\n"
            "   i = 0\n"
            "   for row in reader:\n"
            "       if row[0]=='SM':\n"
            "           sm2id[row[4]] = i\n"
            "           avg.append([i,float(row[1])])\n"
            "           dp.append([i,float(row[2])])\n"
            "           i += 1\n"
            "       elif row[0]=='CN':\n"
            "           val = 0\n"
            "           if int(row[2])!=0: val = float(row[1])/int(row[2])\n"
            "           if not dat:\n"
            "               dat = [[0]*len(sm2id) for x in xrange(len(sm2id))]\n"
            "               min = val\n"
            "               max = val\n"
            "           id_i = sm2id[row[4]]\n"
            "           id_j = sm2id[row[5]]\n"
            "           dat[id_i][id_j] = val\n"
            "           dat[id_j][id_i] = val\n"
            "           if min > val: min = val\n"
            "           if max < val: max = val\n"
            "\n"
            "if len(sm2id)<=1: exit(1)\n"
            "if min==max: exit(1)\n"
            "\n"
            "fig = plt.figure(figsize=(6,7))\n"
            "gs  = gridspec.GridSpec(2, 1, height_ratios=[1, 1.5])\n"
            "ax1 = plt.subplot(gs[0])\n"
            "ax2 = plt.subplot(gs[1])\n"
            "\n"
            "ax1.plot([x[0] for x in avg],[x[1] for x in avg],'^-', ms=3, color='k')\n"
            "ax3 = ax1.twinx()\n"
            "ax3.plot([x[0] for x in dp],[x[1] for x in dp],'^-', ms=3, color='r',mec='r')\n"
            "for tl in ax3.get_yticklabels():\n"
            "   tl.set_color('r')\n"
            "   tl.set_fontsize(9)\n"
            "\n"
            "im = ax2.imshow(dat,clim=(min),interpolation='nearest',origin='lower')\n"
            "cb1  = plt.colorbar(im,ax=ax2)\n"
            "cb1.set_label('Pairwise discordance')\n"
            "for t in cb1.ax.get_yticklabels(): t.set_fontsize(9)\n"
            "\n"
            "ax1.tick_params(axis='both', which='major', labelsize=9)\n"
            "ax1.tick_params(axis='both', which='minor', labelsize=9)\n"
            "ax2.tick_params(axis='both', which='major', labelsize=9)\n"
            "ax2.tick_params(axis='both', which='minor', labelsize=9)\n"
            "\n"
            "ax1.set_title('Sample Discordance Score')\n"
            "ax2.set_ylabel('Sample ID')\n"
            "ax2.set_xlabel('Sample ID')\n"
            "ax3.set_ylabel('Average Depth',color='r')\n"
            "ax1.set_xlabel('Sample ID')\n"
            "ax1.set_ylabel('Average discordance')\n"
            "\n"
            "plt.subplots_adjust(left=0.15,right=0.87,bottom=0.08,top=0.93,hspace=0.25)\n"
            "plt.savefig('%s.png')\n"
            "plt.close()\n"
            "\n", args->plot,args->plot
           );
    fclose(fp);
    py_plot(fname);
    free(fname);
}

static void init_data(args_t *args)
{
    args->sm_hdr = args->files->readers[0].header;
    if ( !bcf_hdr_nsamples(args->sm_hdr) ) error("No samples in %s?\n", args->files->readers[0].fname);

    if ( !args->cross_check )
    {
        args->gt_hdr = args->files->readers[1].header;
        int nsamples = bcf_hdr_nsamples(args->gt_hdr);
        if ( !nsamples ) error("No samples in %s?\n", args->files->readers[1].fname);
        args->lks   = (double*) calloc(nsamples,sizeof(double));
        args->cnts  = (int*) calloc(nsamples,sizeof(int));
        args->sites = (double*) calloc(nsamples,sizeof(double));
        args->dps   = (int*) calloc(nsamples,sizeof(int));
    }
    else
    {
        int nsamples = bcf_hdr_nsamples(args->sm_hdr);
        int narr = (nsamples-1)*nsamples/2;
        args->lks  = (double*) calloc(narr,sizeof(double));
        args->cnts = (int*) calloc(narr,sizeof(int));
        args->dps  = (int*) calloc(narr,sizeof(int));
    }
}

static void destroy_data(args_t *args)
{
    free(args->lks); free(args->cnts); free(args->dps); free(args->cwd); free(args->sites);
}

static int allele_to_int(bcf1_t *line, char *allele)
{
    int i;
    for (i=0; i<line->n_allele; i++)
        if ( !strcmp(allele,line->d.allele[i]) ) return i;
    if ( strcmp(line->d.allele[i-1],"X") ) return -1;
    return i-1;
}

static int init_gt2ipl(args_t *args, bcf1_t *gt_line, bcf1_t *sm_line, int *gt2ipl, int n_gt2ipl)
{
    int i, j;
    for (i=0; i<n_gt2ipl; i++) gt2ipl[i] = -1;
    for (i=0; i<gt_line->n_allele; i++)
    {
        // find which of the sm_alleles (k) corresponds to the gt_allele (i)
        int k = allele_to_int(sm_line, gt_line->d.allele[i]);
        if ( k<0 ) return 0;
        for (j=0; j<=i; j++)
        {
            int l = allele_to_int(sm_line, gt_line->d.allele[j]);
            if ( l<0 ) return 0;
            gt2ipl[ bcf_ij2G(j,i) ] = k<=l ? bcf_ij2G(k,l) : bcf_ij2G(l,k);
        }
    }
    //for (i=0; i<n_gt2ipl; i++) printf("%d .. %d\n", i,gt2ipl[i]);
    return 1;
}

static void set_cwd(args_t *args)
{
    int i;
    char *buf;
    size_t nbuf = 500;
    args->cwd = (char*) malloc(sizeof(char)*nbuf);
    for (i=0; i<5; i++)
    {
        if ( (buf = getcwd(args->cwd, nbuf)) ) break;
        nbuf *= 2;
        args->cwd = (char*) realloc(args->cwd, sizeof(char)*nbuf);
    }
    assert(buf);
}

static void print_header(args_t *args, FILE *fp)
{
    fprintf(fp, "# This file was produced by bcftools (%s+htslib-%s), the command line was:\n", bcftools_version(), hts_version());
    fprintf(fp, "# \t bcftools %s ", args->argv[0]);
    int i;
    for (i=1; i<args->argc; i++)
        fprintf(fp, " %s",args->argv[i]);
    fprintf(fp, "\n# and the working directory was:\n");
    fprintf(fp, "# \t %s\n#\n", args->cwd);
}

static int fake_PLs(args_t *args, bcf_hdr_t *hdr, bcf1_t *line)
{
    // PLs not present, use GTs instead.
    int fake_PL = args->no_PLs ? args->no_PLs : 99;    // with 1, discordance is the number of non-matching GTs
    int nsm_gt, i;
    if ( (nsm_gt=bcf_get_genotypes(hdr, line, &args->tmp_arr, &args->ntmp_arr)) <= 0 )
        error("GT not present at %s:%d?\n", hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1);
    nsm_gt /= bcf_hdr_nsamples(hdr);
    int npl = line->n_allele*(line->n_allele+1)/2;
    hts_expand(int,npl*bcf_hdr_nsamples(hdr),args->npl_arr,args->pl_arr);
    for (i=0; i<bcf_hdr_nsamples(hdr); i++)
    {
        int *gt_ptr = args->tmp_arr + i*nsm_gt;
        int j, *pl_ptr = args->pl_arr + i*npl;
        if ( bcf_gt_is_missing(gt_ptr[0]) || bcf_gt_is_missing(gt_ptr[1]) ) // missing genotype
        {
            for (j=0; j<npl; j++) pl_ptr[j] = -1;
        }
        else
        {
            int a = bcf_gt_allele(gt_ptr[0]);
            int b = bcf_gt_allele(gt_ptr[1]);
            for (j=0; j<npl; j++) pl_ptr[j] = fake_PL;
            int idx = bcf_alleles2gt(a,b);
            pl_ptr[idx] = 0;
        }
    }
    return npl;
}

static int cmp_doubleptr(const void *_a, const void *_b)
{
    double *a = *((double**)_a);
    double *b = *((double**)_b);
    if ( *a < *b ) return -1;
    else if ( *a == *b ) return 0;
    return 1;
}

static void check_gt(args_t *args)
{
    int i,ret, *gt2ipl = NULL, m_gt2ipl = 0, *gt_arr = NULL, ngt_arr = 0;
    int fake_pls = args->no_PLs;

    // Initialize things: check which tags are defined in the header, sample names etc.
    if ( bcf_hdr_id2int(args->gt_hdr, BCF_DT_ID, "GT")<0 ) error("[E::%s] GT not present in the header of %s?\n", __func__, args->files->readers[1].fname);
    if ( bcf_hdr_id2int(args->sm_hdr, BCF_DT_ID, "PL")<0 )
    {
        if ( bcf_hdr_id2int(args->sm_hdr, BCF_DT_ID, "GT")<0 )
            error("[E::%s] Neither PL nor GT present in the header of %s\n", __func__, args->files->readers[0].fname);
        if ( !args->no_PLs )
            fprintf(stderr,"Warning: PL not present in the header of %s, using GT instead\n", args->files->readers[0].fname);
        fake_pls = 1;
    }

    FILE *fp = args->plot ? open_file(NULL, "w", "%s.tab", args->plot) : stdout;
    print_header(args, fp);

    int tgt_isample = -1, query_isample = 0;
    if ( args->target_sample )
    {
        tgt_isample = bcf_hdr_id2int(args->gt_hdr, BCF_DT_SAMPLE, args->target_sample);
        if ( tgt_isample<0 ) error("No such sample in %s: [%s]\n", args->files->readers[1].fname, args->target_sample);
    }
    if ( args->all_sites )
    {
        if ( tgt_isample==-1 )
        {
            fprintf(stderr,"No target sample selected for comparison, using the first sample in %s: %s\n", args->gt_fname,args->gt_hdr->samples[0]);
            tgt_isample = 0;
        }
    }
    if ( args->query_sample )
    {
        query_isample = bcf_hdr_id2int(args->sm_hdr, BCF_DT_SAMPLE, args->query_sample);
        if ( query_isample<0 ) error("No such sample in %s: [%s]\n", args->files->readers[0].fname, args->query_sample);
    }
    if ( args->all_sites )
        fprintf(fp, "# [1]SC, Site by Site Comparison\t[2]Chromosome\t[3]Position\t[4]-g alleles\t[5]-g GT (%s)\t[6]match log LK\t[7]Query alleles\t[8-]Query PLs (%s)\n",
                args->gt_hdr->samples[tgt_isample],args->sm_hdr->samples[query_isample]);

    // Main loop
    float prev_lk = 0;
    while ( (ret=bcf_sr_next_line(args->files)) )
    {
        if ( ret!=2 ) continue;
        bcf1_t *sm_line = args->files->readers[0].buffer[0];    // the query file
        bcf1_t *gt_line = args->files->readers[1].buffer[0];    // the -g target file
        bcf_unpack(sm_line, BCF_UN_FMT);
        bcf_unpack(gt_line, BCF_UN_FMT);

        // Init mapping from target genotype index to the sample's PL fields
        int n_gt2ipl = gt_line->n_allele*(gt_line->n_allele + 1)/2;
        if ( n_gt2ipl > m_gt2ipl )
        {
            m_gt2ipl = n_gt2ipl;
            gt2ipl   = (int*) realloc(gt2ipl, sizeof(int)*m_gt2ipl);
        }
        if ( !init_gt2ipl(args, gt_line, sm_line, gt2ipl, n_gt2ipl) ) continue;

        // Target genotypes
        int ngt, npl;
        if ( (ngt=bcf_get_genotypes(args->gt_hdr, gt_line, &gt_arr, &ngt_arr)) <= 0 )
            error("GT not present at %s:%d?", args->gt_hdr->id[BCF_DT_CTG][gt_line->rid].key, gt_line->pos+1);
        ngt /= bcf_hdr_nsamples(args->gt_hdr);
        if ( ngt!=2 ) continue; // checking only diploid genotypes

        // Sample PLs
        if ( !fake_pls )
        {
            if ( (npl=bcf_get_format_int32(args->sm_hdr, sm_line, "PL", &args->pl_arr, &args->npl_arr)) <= 0 )
            {
                if ( sm_line->n_allele==1 )
                {
                    // PL values may not be present when ALT=. (mpileup/bcftools output), in that case 
                    // switch automatically to GT at these sites
                    npl = fake_PLs(args, args->sm_hdr, sm_line);
                }
                else
                    error("PL not present at %s:%d?\n", args->sm_hdr->id[BCF_DT_CTG][sm_line->rid].key, sm_line->pos+1);
            }
            else
                npl /= bcf_hdr_nsamples(args->sm_hdr);
        }
        else
            npl = fake_PLs(args, args->sm_hdr, sm_line);

        // Calculate likelihoods for all samples, assuming diploid genotypes

        // For faster access to genotype likelihoods (PLs) of the query sample
        int max_ipl, *pl_ptr = args->pl_arr + query_isample*npl;
        double sum_pl = 0; // for converting PLs to probs
        for (max_ipl=0; max_ipl<npl; max_ipl++)
        {
            if ( pl_ptr[max_ipl]==bcf_int32_vector_end ) break;
            if ( pl_ptr[max_ipl]==bcf_int32_missing ) continue;
            sum_pl += pow(10, -0.1*pl_ptr[max_ipl]);
        }
        if ( sum_pl==0 ) continue; // no PLs present
        if ( fake_pls && args->no_PLs==1 ) sum_pl = -1;

        // The main stats: concordance of the query sample with the target -g samples
        for (i=0; i<bcf_hdr_nsamples(args->gt_hdr); i++)
        {
            int *gt_ptr = gt_arr + i*ngt;
            if ( gt_ptr[1]==bcf_int32_vector_end ) continue;    // skip haploid genotypes
            if ( bcf_gt_is_missing(gt_ptr[0]) || bcf_gt_is_missing(gt_ptr[1]) ) continue;
            int a = bcf_gt_allele(gt_ptr[0]);
            int b = bcf_gt_allele(gt_ptr[1]);
            if ( args->hom_only && a!=b ) continue; // heterozygous genotype
            int igt_tgt = igt_tgt = bcf_alleles2gt(a,b); // genotype index in the target file
            int igt_qry = gt2ipl[igt_tgt];  // corresponding genotype in query file
            if ( igt_qry>=max_ipl || pl_ptr[igt_qry]<0 ) continue;   // genotype not present in query sample: haploid or missing
            args->lks[i] += sum_pl<0 ? -pl_ptr[igt_qry] : log(pow(10, -0.1*pl_ptr[igt_qry])/sum_pl);
            args->sites[i]++;
        }
        if ( args->all_sites )
        {
            // Print LKs at all sites for debugging
            int *gt_ptr = gt_arr + tgt_isample*ngt;
            if ( gt_ptr[1]==bcf_int32_vector_end ) continue;    // skip haploid genotypes
            int a = bcf_gt_allele(gt_ptr[0]);
            int b = bcf_gt_allele(gt_ptr[1]);
            if ( args->hom_only && a!=b ) continue; // heterozygous genotype
            fprintf(fp, "SC\t%s\t%d", args->gt_hdr->id[BCF_DT_CTG][gt_line->rid].key, gt_line->pos+1);
            for (i=0; i<gt_line->n_allele; i++) fprintf(fp, "%c%s", i==0?'\t':',', gt_line->d.allele[i]);
            fprintf(fp, "\t%s/%s", a>=0 ? gt_line->d.allele[a] : ".", b>=0 ? gt_line->d.allele[b] : ".");
            fprintf(fp, "\t%f", args->lks[query_isample]-prev_lk);
            prev_lk = args->lks[query_isample];

            int igt, *pl_ptr = args->pl_arr + query_isample*npl; // PLs of the query sample
            for (i=0; i<sm_line->n_allele; i++) fprintf(fp, "%c%s", i==0?'\t':',', sm_line->d.allele[i]);
            for (igt=0; igt<npl; igt++)
                if ( pl_ptr[igt]==bcf_int32_vector_end ) break;
                else if ( pl_ptr[igt]==bcf_int32_missing ) fprintf(fp, ".");
                else fprintf(fp, "\t%d", pl_ptr[igt]);
            fprintf(fp, "\n");
        }
    }
    free(gt2ipl);
    free(gt_arr);
    free(args->pl_arr);
    free(args->tmp_arr);

    // To be able to plot total discordance (=number of mismatching GTs with -G1) in the same
    // plot as discordance per site, the latter must be scaled to the same range
    int nsamples = bcf_hdr_nsamples(args->gt_hdr);
    double extreme_lk = 0, extreme_lk_per_site = 0;
    for (i=0; i<nsamples; i++)
    {
        if ( args->lks[i] < extreme_lk ) extreme_lk = args->lks[i];
        if ( args->sites[i] && args->lks[i]/args->sites[i] < extreme_lk_per_site ) extreme_lk_per_site = args->lks[i]/args->sites[i];
    }

    // Sorted output
    double **p = (double**) malloc(sizeof(double*)*nsamples);
    for (i=0; i<nsamples; i++) p[i] = &args->lks[i];
    qsort(p, nsamples, sizeof(int*), cmp_doubleptr);

    fprintf(fp, "# [1]CN\t[2]Discordance with %s (total)\t[3]Discordance (avg score per site)\t[4]Number of sites compared\t[5]Sample\t[6]Sample ID\n", args->sm_hdr->samples[query_isample]);
    for (i=0; i<nsamples; i++)
    {
        int idx = p[i] - args->lks;
        double per_site = 0;
        if ( args->sites[idx] )
        {
            if ( args->sites[idx] && extreme_lk_per_site )
            {
                per_site = args->lks[idx]/args->sites[idx];
                per_site *= extreme_lk / extreme_lk_per_site;
            }
            else
                per_site = 0;
        }
        fprintf(fp, "CN\t%e\t%e\t%.0f\t%s\t%d\n", fabs(args->lks[idx]), fabs(per_site), args->sites[idx], args->gt_hdr->samples[idx], i);
    }

    if ( args->plot )
    {
        fclose(fp);
        plot_check(args, args->target_sample ? args->target_sample : "", args->sm_hdr->samples[query_isample]);
    }
}

static inline int is_hom_most_likely(int nals, int *pls)
{
    int ia, ib, idx = 1, min_is_hom = 1, min_pl = pls[0];
    for (ia=1; ia<nals; ia++)
    {
        for (ib=0; ib<ia; ib++)
        {
            if ( pls[idx] < min_pl ) { min_pl = pls[idx]; min_is_hom = 0; }
            idx++;
        }
        if ( pls[idx] < min_pl ) { min_pl = pls[idx]; min_is_hom = 1; }
        idx++;
    }
    return min_is_hom;
}

static void cross_check_gts(args_t *args)
{
    int nsamples = bcf_hdr_nsamples(args->sm_hdr), ndp_arr = 0;
    unsigned int *dp = (unsigned int*) calloc(nsamples,sizeof(unsigned int)), *ndp = (unsigned int*) calloc(nsamples,sizeof(unsigned int)); // this will overflow one day...
    int fake_pls = args->no_PLs, ignore_dp = 0;

    int i,j,k,idx, pl_warned = 0, dp_warned = 0;
    int32_t *dp_arr = NULL;
    int *is_hom = args->hom_only ? (int*) malloc(sizeof(int)*nsamples) : NULL;
    if ( bcf_hdr_id2int(args->sm_hdr, BCF_DT_ID, "PL")<0 )
    {
        if ( bcf_hdr_id2int(args->sm_hdr, BCF_DT_ID, "GT")<0 )
            error("[E::%s] Neither PL nor GT present in the header of %s\n", __func__, args->files->readers[0].fname);
        if ( !args->no_PLs )
            fprintf(stderr,"Warning: PL not present in the header of %s, using GT instead\n", args->files->readers[0].fname);
        fake_pls = 1;
    }
    if ( bcf_hdr_id2int(args->sm_hdr, BCF_DT_ID, "DP")<0 ) ignore_dp = 1;

    FILE *fp = args->plot ? open_file(NULL, "w", "%s.tab", args->plot) : stdout;
    print_header(args, fp);
    if ( args->all_sites ) fprintf(fp,"# [1]SD, Average Site Discordance\t[2]Chromosome\t[3]Position\t[4]Number of available pairs\t[5]Average discordance\n");

    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = args->files->readers[0].buffer[0];
        bcf_unpack(line, BCF_UN_FMT);

        int npl;
        if ( !fake_pls )
        {
            npl = bcf_get_format_int32(args->sm_hdr, line, "PL", &args->pl_arr, &args->npl_arr);
            if ( npl<=0 ) { pl_warned++; continue; }
            npl /= nsamples;
        }
        else
            npl = fake_PLs(args, args->sm_hdr, line);
        int mdp = 0;
        if ( !ignore_dp && (mdp=bcf_get_format_int32(args->sm_hdr, line, "DP", &dp_arr, &ndp_arr)) <= 0 ) dp_warned++;

        if ( args->hom_only )
        {
            for (i=0; i<nsamples; i++)
                is_hom[i] = is_hom_most_likely(line->n_allele, args->pl_arr+i*npl);
        }

        double sum = 0; int nsum = 0;
        idx = 0;
        for (i=0; i<nsamples; i++)
        {
            int *ipl = &args->pl_arr[i*npl];
            if ( *ipl==-1 ) { idx += i; continue; } // missing genotype
            if ( mdp>0 && (dp_arr[i]==bcf_int32_missing || !dp_arr[i]) ) { idx += i; continue; }
            if ( args->hom_only && !is_hom[i] ) { idx += i; continue; }

            for (j=0; j<i; j++)
            {
                int *jpl = &args->pl_arr[j*npl];
                if ( *jpl==-1 ) { idx++; continue; } // missing genotype
                if ( mdp>0 && (dp_arr[j]==bcf_int32_missing || !dp_arr[j]) ) { idx++; continue; }
                if ( args->hom_only && !is_hom[j] ) { idx++; continue; }

                int min_pl = INT_MAX;
                for (k=0; k<npl; k++)
                {
                    if ( ipl[k]==bcf_int32_missing || jpl[k]==bcf_int32_missing ) break;
                    if ( ipl[k]==bcf_int32_vector_end || jpl[k]==bcf_int32_vector_end ) { k = npl; break; }
                    if ( min_pl > ipl[k]+jpl[k] ) min_pl = ipl[k]+jpl[k];
                }
                if ( k!=npl ) { idx++; continue; }

                if ( args->all_sites ) { sum += min_pl; nsum++; }
                args->lks[idx] += min_pl;
                args->cnts[idx]++;

                if ( mdp>0 )
                {
                    args->dps[idx] += dp_arr[i] < dp_arr[j] ? dp_arr[i] : dp_arr[j];
                    dp[i] += dp_arr[i]; ndp[i]++;
                    dp[j] += dp_arr[j]; ndp[j]++;
                }
                else
                {
                    args->dps[idx]++;
                    dp[i]++; ndp[i]++;
                    dp[j]++; ndp[j]++;
                }
                idx++;
            }
        }
        if ( args->all_sites )
            fprintf(fp,"SD\t%s\t%d\t%d\t%.0f\n", args->sm_hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1, nsum, nsum?sum/nsum:0);
    }
    if ( dp_arr ) free(dp_arr);
    if ( args->pl_arr ) free(args->pl_arr);
    if ( args->tmp_arr ) free(args->tmp_arr);
    if ( is_hom ) free(is_hom);

    if ( pl_warned ) fprintf(stderr, "[W::%s] PL was not found at %d site(s)\n", __func__, pl_warned);
    if ( dp_warned ) fprintf(stderr, "[W::%s] DP was not found at %d site(s)\n", __func__, dp_warned);

    // Output samples sorted by average discordance
    double *score  = (double*) calloc(nsamples,sizeof(double));
    args->sites = (double*) calloc(nsamples,sizeof(double));
    idx = 0;
    for (i=0; i<nsamples; i++)
    {
        for (j=0; j<i; j++)
        {
            score[i] += args->lks[idx];
            score[j] += args->lks[idx];
            args->sites[i] += args->cnts[idx];
            args->sites[j] += args->cnts[idx];
            idx++;
        }
    }
    for (i=0; i<nsamples; i++)
        if ( args->sites[i] ) score[i] /= args->sites[i];
    double **p = (double**) malloc(sizeof(double*)*nsamples), avg_score = 0;
    for (i=0; i<nsamples; i++) p[i] = &score[i];
    qsort(p, nsamples, sizeof(int*), cmp_doubleptr);
    // The average discordance gives the number of differing sites in % with -G1
    fprintf(fp, "# [1]SM\t[2]Average Discordance\t[3]Average depth\t[4]Average number of sites\t[5]Sample\t[6]Sample ID\n");
    for (i=0; i<nsamples; i++)
    {
        idx = p[i] - score;
        double adp = ndp[idx] ? (double)dp[idx]/ndp[idx] : 0;
        double nsites = args->sites[idx]/(nsamples-1);
        avg_score += score[idx];
        fprintf(fp, "SM\t%f\t%.2lf\t%.0lf\t%s\t%d\n", score[idx]*100., adp, nsites, args->sm_hdr->samples[idx],i);
    }

    //  // Overall score: maximum absolute deviation from the average score
    //  fprintf(fp, "# [1] MD\t[2]Maximum deviation\t[3]The culprit\n");
    //  fprintf(fp, "MD\t%f\t%s\n", (score[idx] - avg_score/nsamples)*100., args->sm_hdr->samples[idx]);    // idx still set
    free(p);
    free(score);
    free(dp);
    free(ndp);

    // Pairwise discordances
    fprintf(fp, "# [1]CN\t[2]Discordance\t[3]Number of sites\t[4]Average minimum depth\t[5]Sample i\t[6]Sample j\n");
    idx = 0;
    for (i=0; i<nsamples; i++)
    {
        for (j=0; j<i; j++)
        {
            fprintf(fp, "CN\t%.0f\t%d\t%.2f\t%s\t%s\n", args->lks[idx], args->cnts[idx], args->cnts[idx]?(double)args->dps[idx]/args->cnts[idx]:0.0,
                    args->sm_hdr->samples[i],args->sm_hdr->samples[j]);
            idx++;
        }
    }
    fclose(fp);
    if ( args->plot )
        plot_cross_check(args);
}

static char *init_prefix(char *prefix)
{
    int len = strlen(prefix);
    if ( prefix[len-1] == '/' || prefix[len-1] == '\\' )
        return msprintf("%sgtcheck", prefix);
    return strdup(prefix);
}

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Check sample identity. With no -g BCF given, multi-sample cross-check is performed.\n");
    fprintf(stderr, "Usage:   bcftools gtcheck [options] [-g <genotypes.vcf.gz>] <query.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -a, --all-sites                 output comparison for all sites\n");
    fprintf(stderr, "    -g, --genotypes <file>          genotypes to compare against\n");
    fprintf(stderr, "    -G, --GTs-only <int>            use GTs, ignore PLs, using <int> for unseen genotypes [99]\n");
    fprintf(stderr, "    -H, --homs-only                 homozygous genotypes only (useful for low coverage data)\n");
    fprintf(stderr, "    -p, --plot <prefix>             plot\n");
    fprintf(stderr, "    -r, --regions <region>          restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>       restrict to regions listed in a file\n");
    fprintf(stderr, "    -s, --query-sample <string>     query sample (by default the first sample is checked)\n");
    fprintf(stderr, "    -S, --target-sample <string>    target sample in the -g file (used only for plotting)\n");
    fprintf(stderr, "    -t, --targets <region>          similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file <file>       similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfgtcheck(int argc, char *argv[])
{
    int c;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->files  = bcf_sr_init();
    args->argc   = argc; args->argv = argv; set_cwd(args);
    char *regions = NULL, *targets = NULL;
    int regions_is_file = 0, targets_is_file = 0;

    static struct option loptions[] =
    {
        {"GTs-only",1,0,'G'},
        {"all-sites",0,0,'a'},
        {"homs-only",0,0,'H'},
        {"help",0,0,'h'},
        {"genotypes",1,0,'g'},
        {"plot",1,0,'p'},
        {"target-sample",1,0,'S'},
        {"query-sample",1,0,'s'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {0,0,0,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "hg:p:s:S:Hr:R:at:T:G:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'G':
                args->no_PLs = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --GTs-only %s\n", optarg);
                break;
            case 'a': args->all_sites = 1; break;
            case 'H': args->hom_only = 1; break;
            case 'g': args->gt_fname = optarg; break;
            case 'p': args->plot = optarg; break;
            case 'S': args->target_sample = optarg; break;
            case 's': args->query_sample = optarg; break;
            case 'r': regions = optarg; break;
            case 'R': regions = optarg; regions_is_file = 1; break;
            case 't': targets = optarg; break;
            case 'T': targets = optarg; targets_is_file = 1; break;
            case 'h':
            case '?': usage();
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    char *fname = NULL;
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else usage();   // no files given
    }
    else fname = argv[optind];
    if ( argc>optind+1 )  usage();  // too many files given
    if ( !args->gt_fname ) args->cross_check = 1;   // no genotype file, run in cross-check mode
    else args->files->require_index = 1;
    if ( regions && bcf_sr_set_regions(args->files, regions, regions_is_file)<0 ) error("Failed to read the regions: %s\n", regions);
    if ( targets && bcf_sr_set_targets(args->files, targets, targets_is_file, 0)<0 ) error("Failed to read the targets: %s\n", targets);
    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->files->errnum));
    if ( args->gt_fname && !bcf_sr_add_reader(args->files, args->gt_fname) ) error("Failed to open %s: %s\n", args->gt_fname,bcf_sr_strerror(args->files->errnum));
    args->files->collapse = COLLAPSE_SNPS|COLLAPSE_INDELS;
    if ( args->plot ) args->plot = init_prefix(args->plot);
    init_data(args);
    if ( args->cross_check )
        cross_check_gts(args);
    else
        check_gt(args);
    destroy_data(args);
    bcf_sr_destroy(args->files);
    if (args->plot) free(args->plot);
    free(args);
    return 0;
}

