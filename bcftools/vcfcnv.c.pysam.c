#include "bcftools.pysam.h"

/* The MIT License

   Copyright (c) 2014-2022 Genome Research Ltd.

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

/*
    Known issues:
    - The --AF-file option behaves like --targets-file, sites not listed in the AFs
      are skipped.
*/

#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <inttypes.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kstring.h>
#include <htslib/kfunc.h>
#include <htslib/khash_str2int.h>
#include "bcftools.h"
#include "HMM.h"
#include "rbuf.h"

#define DBG0 0

#define N_STATES 4
#define CN0 0
#define CN1 1
#define CN2 2
#define CN3 3

typedef struct
{
    float mean, dev2, norm;
}
gauss_param_t;

typedef struct
{
    char *name;
    int idx;    // VCF sample index
    float *lrr,*baf, baf_dev2, baf_dev2_dflt, lrr_dev2;
    float cell_frac, cell_frac_dflt;
    gauss_param_t gauss_param[18];
    double pobs[N_STATES];
    FILE *dat_fh, *cn_fh, *summary_fh;
    char *dat_fname, *cn_fname, *summary_fname;
}
sample_t;

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr;
    int prev_rid, ntot, nused;
    sample_t query_sample, control_sample;

    int nstates;    // number of states: N_STATES for one sample, N_STATES^2 for two samples
    double lrr_bias, baf_bias;              // LRR/BAF weights
    double same_prob, ij_prob;              // prior of both samples being the same and the transition probability P(i|j)
    double err_prob;                        // constant probability of erroneous measurement
    float *nonref_afs, nonref_af, nonref_af_dflt, fRR, fRA, fAA, *tmpf;
    unsigned long int nRR, nRA, nAA;
    int mtmpf;

    double *tprob, *tprob_arr;  // array of transition matrices, precalculated up to ntprob_arr positions
    double *iprobs;             // states' initial probabilities
    int ntprob_arr;

    hmm_t *hmm;
    double *eprob;          // emission probs [nstates*nsites,msites]
    uint32_t *sites;        // positions [nsites,msites]
    int nsites, msites;

    double baum_welch_th, optimize_frac; 
    float plot_th;
    FILE *summary_fh;
    char **argv, *regions_list, *summary_fname, *output_dir;
    char *targets_list, *af_fname;
    int argc, verbose, lrr_smooth_win;
}
args_t;

FILE *open_file(char **fname, const char *mode, const char *fmt, ...);

static inline void hmm2cn_state(int nstates, int i, int *a, int *b)
{
    *a = i / N_STATES;
    *b = i - (*a)*N_STATES;
}
static double *init_tprob_matrix(int ndim, double ij_prob, double same_prob)
{
    int i,j;
    double *mat = (double*) malloc(sizeof(double)*ndim*ndim);

    assert( ndim==N_STATES || ndim==N_STATES*N_STATES);

    if ( ndim==N_STATES )   // one sample
    {
        double pii = 1 - ij_prob*(N_STATES-1);
        if ( pii < ij_prob ) error("Error: -x set a bit too high, P(x|x) < P(x|y): %e vs %e\n", pii,ij_prob);
        for (j=0; j<ndim; j++)
        {
            double sum = 0;
            for (i=0; i<ndim; i++)
            {
                // transition from j-th to i-th state
                if ( i==j )
                    MAT(mat,ndim,i,j) = pii;
                else
                    MAT(mat,ndim,i,j) = ij_prob;

                sum += MAT(mat,ndim,i,j);
            }
            assert( fabs(sum - 1.0)<1e-15 );
        }
    }
    else    // two samples
    {
        // interpret ij_prob differently, as ii_prob in fact, so that for two
        // samples the behaviour is somewhat closer to single sample calling
        // with s=0. 
        double pii = 1 - ij_prob*(N_STATES-1);
        ij_prob = (1 - pii) / (ndim - 1);
        for (j=0; j<ndim; j++)
        {
            int ja,jb;
            hmm2cn_state(ndim, j, &ja, &jb);

            double sum = 0;
            for (i=0; i<ndim; i++)
            {
                int ia,ib;
                hmm2cn_state(ndim, i, &ia, &ib);

                // transition from (ja,jb)-th to (ia,ib)-th state
                double pa = ja==ia ? pii : ij_prob;
                double pb = jb==ib ? pii : ij_prob;

                if ( ia==ib && ja==jb )
                    MAT(mat,ndim,i,j) = pa*pb - pa*pb*same_prob + sqrt(pa*pb)*same_prob;
                else if ( ia==ib )
                    MAT(mat,ndim,i,j) = pa*pb;
                else
                    MAT(mat,ndim,i,j) = pa*pb*(1-same_prob);

                sum += MAT(mat,ndim,i,j);
            }
            for (i=0; i<ndim; i++) MAT(mat,ndim,i,j) /= sum;
        }
    }
    return mat;
}

static double *init_iprobs(int ndim, double same_prob)
{
    int i;
    double *probs = (double*) malloc(sizeof(double)*ndim);

    assert( ndim==N_STATES || ndim==N_STATES*N_STATES);

    if ( ndim==N_STATES )   
    {
        // one sample: prior on CN2
        for (i=0; i<ndim; i++) 
            probs[i] = i==CN2 ? 0.5 : 0.5/3;
    }
    else
    {
        // two samples
        double norm = 0;
        for (i=0; i<ndim; i++) 
        {
            int ia,ib;
            hmm2cn_state(ndim, i, &ia, &ib);

            double pa = ia==CN2 ? 0.5 : 0.5/3;
            double pb = ib==CN2 ? 0.5 : 0.5/3;

            probs[i] = pa*pb;
            if ( ia!=ib ) probs[i] *= 1-same_prob;

            norm += probs[i];
        }
        for (i=0; i<ndim; i++) probs[i] /= norm;
    }
    return probs;
}

static void init_sample_files(sample_t *smpl, char *dir)
{
    smpl->dat_fh = open_file(&smpl->dat_fname,"w","%s/dat.%s.tab",dir,smpl->name);
    if ( !smpl->dat_fh ) error("Error opening file: %s/dat.%s.tab\n",dir,smpl->name);

    smpl->cn_fh  = open_file(&smpl->cn_fname,"w","%s/cn.%s.tab",dir,smpl->name);
    if ( !smpl->cn_fh ) error("Error opening file: %s/cn.%s.tab\n",dir,smpl->name);

    smpl->summary_fh = open_file(&smpl->summary_fname,"w","%s/summary.%s.tab",dir,smpl->name);
    if ( !smpl->summary_fh ) error("Error opening file: %s/summary.%s.tab\n",dir,smpl->name);

    fprintf(smpl->dat_fh,"# [1]Chromosome\t[2]Position\t[3]BAF\t[4]LRR\n");
    fprintf(smpl->cn_fh,"# [1]Chromosome\t[2]Position\t[3]CN\t[4]P(CN0)\t[5]P(CN1)\t[6]P(CN2)\t[7]P(CN3)\n");
    fprintf(smpl->summary_fh,"# RG, Regions [2]Chromosome\t[3]Start\t[4]End\t[5]Copy Number state\t[6]Quality\t[7]nSites\t[8]nHETs\n");
}
static void close_sample_files(sample_t *smpl)
{
    if ( fclose(smpl->dat_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,smpl->dat_fname);
    if ( fclose(smpl->cn_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,smpl->cn_fname);
    if ( fclose(smpl->summary_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,smpl->summary_fname);
}

static double norm_cdf(double mean, double dev);
static void init_data(args_t *args)
{
    args->prev_rid = -1;
    args->hdr = args->files->readers[0].header;

    if ( !args->query_sample.name )
    {
        if ( bcf_hdr_nsamples(args->hdr)>1 ) error("Multi-sample VCF, missing the -s option\n");
        args->query_sample.name = strdup(args->hdr->samples[0]);
    }
    else 
        if ( bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,args->query_sample.name)<0 ) error("The sample \"%s\" not found\n", args->query_sample.name);
    if ( !args->files->readers[0].file->is_bin )
    {
        int ret;
        kstring_t tmp = {0,0,0};
        if ( args->control_sample.name )
        {
            ksprintf(&tmp, "%s,%s", args->query_sample.name,args->control_sample.name);
            ret = bcf_hdr_set_samples(args->hdr, tmp.s, 0);
        }
        else
        {
            ret = bcf_hdr_set_samples(args->hdr, args->query_sample.name, 0);
            tmp.s = args->query_sample.name;
        }
        if ( ret<0 ) error("Error parsing the list of samples: %s\n", tmp.s);
        else if ( ret>0 ) error("The sample not found in the VCF: %s\n", ret==1 ? args->query_sample.name : args->control_sample.name);

        if ( args->control_sample.name ) free(tmp.s);
    }
    args->query_sample.idx = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,args->query_sample.name);
    args->control_sample.idx = args->control_sample.name ? bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,args->control_sample.name) : -1;
    args->nstates = args->control_sample.name ? N_STATES*N_STATES : N_STATES;
    args->tprob  = init_tprob_matrix(args->nstates, args->ij_prob, args->same_prob);
    args->iprobs = init_iprobs(args->nstates, args->same_prob);
    args->hmm = hmm_init(args->nstates, args->tprob, 10000);
    hmm_init_states(args->hmm, args->iprobs);

    args->summary_fh = bcftools_stdout;
    init_sample_files(&args->query_sample, args->output_dir);
    if ( args->control_sample.name )
    {
        init_sample_files(&args->control_sample, args->output_dir);
        args->summary_fh = open_file(&args->summary_fname,"w","%s/summary.tab",args->output_dir);
    }
    else
        args->summary_fh = NULL;    // one sample only, no two-file summary
        

    int i;
    FILE *fh = args->summary_fh ? args->summary_fh : args->query_sample.summary_fh;

    fprintf(fh, "# This file was produced by: bcftools cnv(%s+htslib-%s)\n", bcftools_version(),hts_version());
    fprintf(fh, "# The command line was:\tbcftools %s", args->argv[0]);
    for (i=1; i<args->argc; i++) fprintf(fh, " %s",args->argv[i]);
    if ( args->control_sample.name )
        fprintf(fh, "\n#\n"
                "# RG, Regions\t[2]Chromosome\t[3]Start\t[4]End\t[5]Copy number:%s\t[6]Copy number:%s\t[7]Quality"
                "\t[8]nSites in (5)\t[9]nHETs in (5)\t[10]nSites in (6)\t[11]nHETs in(6)\n",
                args->query_sample.name,args->control_sample.name
               );
    else
        fprintf(fh, "\n#\n"
                "# RG, Regions\t[2]Chromosome\t[3]Start\t[4]End\t[5]Copy number:%s\t[6]Quality\t[7]nSites\t[8]nHETs\n",
                args->query_sample.name
               );
    if ( args->optimize_frac )
    {
        fprintf(args->query_sample.summary_fh, "# CF, cell fraction estimate\t[2]Chromosome\t[3]Start\t[4]End\t[5]Cell fraction\t[6]BAF deviation\n");
        if ( args->control_sample.name )
        {
            fprintf(args->control_sample.summary_fh, "# CF, cell fraction estimate\t[2]Chromosome\t[3]Start\t[4]End\t[5]Cell fraction\t[6]BAF deviation\n");
            fprintf(args->summary_fh, "# CF, cell fraction estimate\t[2]Chromosome\t[3]Start\t[4]End\t"
                "[5]Cell fraction:%s\t[6]Cell fraction:%s\t[7]BAF deviation:%s\t[8]BAF deviation:%s\n",
                args->query_sample.name,args->control_sample.name,
                args->query_sample.name,args->control_sample.name
                );
        }
    }
}

char *msprintf(const char *fmt, ...);
static void py_plot_cnv(char *script, float th)
{
    if ( th>100 ) return;   // create no plots

    char *cmd = msprintf("python %s -p %f", script, th);
    int ret = system(cmd);
    if ( ret) fprintf(bcftools_stderr, "The command returned non-zero status %d: %s\n", ret, cmd);
    free(cmd);
}

static void plot_sample(args_t *args, sample_t *smpl)
{
    char *fname;
    FILE *fp = open_file(&fname,"w","%s/plot.%s.py",args->output_dir,smpl->name);
    fprintf(fp,
            "import matplotlib as mpl\n"
            "mpl.use('Agg')\n"
            "import matplotlib.pyplot as plt\n"
            "import csv\n"
            "import numpy as np\n"
            "csv.register_dialect('tab', delimiter='\\t', quoting=csv.QUOTE_NONE)\n"
            "\n"
            "dat = {}\n"
            "with open('%s', 'r') as f:\n"
            "    reader = csv.reader(f, 'tab')\n"
            "    for row in reader:\n"
            "        chr = row[0]\n"
            "        if chr[0]=='#': continue\n"
            "        if chr not in dat: dat[chr] = []\n"
            "        dat[chr].append([int(row[1]), float(row[2]), float(row[3])])\n"
            "\n"
            "cnv = {}\n"
            "with open('%s', 'r') as f:\n"
            "    reader = csv.reader(f, 'tab')\n"
            "    for row in reader:\n"
            "        chr = row[0]\n"
            "        if chr[0]=='#': continue\n"
            "        if chr not in cnv: cnv[chr] = []\n"
            "        row[2] = int(row[2]) + 0.5\n"
            "        row[1] = int(row[1])\n"
            "        cnv[chr].append(row[1:])\n"
            "\n"
            "for chr in dat:\n"
            "    fig,(ax1, ax2, ax3) = plt.subplots(3,1,figsize=(10,8),sharex=True)\n"
            "    ax1.plot([x[0] for x in dat[chr]],[x[2] for x in dat[chr]], '.', ms=3)\n"
            "    ax2.plot([x[0] for x in dat[chr]],[x[1] for x in dat[chr]], '.', ms=3)\n"
            "    cn_dat = cnv[chr]\n"
            "    xgrid = [float(x[0]) for x in cn_dat]\n"
            "    ygrid = np.linspace(0,5,6)\n"
            "    xgrid, ygrid = np.meshgrid(xgrid, ygrid)\n"
            "    heat = np.zeros_like(xgrid)\n"
            "    for x in range(len(heat[0])-1):\n"
            "       heat[0][x] = cn_dat[x][2]\n"
            "       heat[1][x] = cn_dat[x][3]\n"
            "       heat[2][x] = cn_dat[x][4]\n"
            "       heat[3][x] = cn_dat[x][5]\n"
            "    mesh = ax3.pcolormesh(xgrid, ygrid, heat, cmap='bwr_r', shading='auto', alpha=0)\n"
            "    mesh.set_clim(vmin=-1,vmax=1)\n"
            "    ax3.plot([x[0] for x in cn_dat],[x[1] for x in cn_dat],'.-',ms=3,color='black')\n"
            "    fig.suptitle('%s (chr '+chr+')')\n"
            "    ax1.set_title('Log-R intensities Ratio',fontsize=10)\n"
            "    ax2.set_title('B-Allele Frequency',fontsize=10)\n"
            "    ax3.set_title('Copy Number Variation',fontsize=10)\n"
            "    ax1.set_ylabel('LRR')\n"
            "    ax2.set_ylabel('BAF')\n"
            "    ax3.set_ylabel('CN')\n"
            "    ax3.set_xlabel('Coordinate (chrom '+chr+')',fontsize=10)\n"
            "    ax3.set_ylim(-0.1,4.1)\n"
            "    ax3.set_yticks([0.5,1.5,2.5,3.5])\n"
            "    ax3.set_yticklabels(['CN0','CN1','CN2','CN3'])\n"
            "    plt.subplots_adjust(left=0.08,right=0.95,bottom=0.08,top=0.92)\n"
            "    plt.savefig('%s/plot.%s.chr'+chr+'.png')\n"
            "    plt.close()\n"
            "\n", 
            smpl->dat_fname,smpl->cn_fname,smpl->name,args->output_dir,smpl->name
    );
    fclose(fp);

    py_plot_cnv(fname, args->plot_th);
    free(fname);
}

static void create_plots(args_t *args)
{
    close_sample_files(&args->query_sample);
    if ( args->control_sample.name ) close_sample_files(&args->control_sample);
    if ( args->summary_fh ) fclose(args->summary_fh);

    if ( !args->control_sample.name )
    {
        plot_sample(args, &args->query_sample);
        return;
    }

    char *fname;
    FILE *fp = open_file(&fname,"w","%s/plot.%s.%s.py",args->output_dir,args->control_sample.name,args->query_sample.name);
    fprintf(fp,
            "import matplotlib as mpl\n"
            "mpl.use('Agg')\n"
            "import matplotlib.pyplot as plt\n"
            "import csv,argparse\n"
            "import numpy as np\n"
            "csv.register_dialect('tab', delimiter='\\t', quoting=csv.QUOTE_NONE)\n"
            "\n"
            "control_sample = '%s'\n"
            "query_sample   = '%s'\n"
            "\n"
            "parser = argparse.ArgumentParser()\n"
            "parser.add_argument('-p', '--plot-threshold', type=float)\n"
            "parser.add_argument('-c', '--chromosome')\n"
            "args = parser.parse_args()\n"
            "if args.plot_threshold==None: args.plot_threshold = 0\n"
            "\n"
            "def chroms_to_plot(th):\n"
            "   dat = {}\n"
            "   with open('%s/summary.tab', 'r') as f:\n"
            "       reader = csv.reader(f, 'tab')\n"
            "       for row in reader:\n"
            "           if row[0]!='RG': continue\n"
            "           chr   = row[1]\n"
            "           start = row[2]\n"
            "           end   = row[3]\n"
            "           qual  = float(row[6])\n"
            "           if row[4]==row[5] and args.plot_threshold!=0: continue\n"
            "           if chr not in dat: dat[chr] = 0.0\n"
            "           if qual > dat[chr]: dat[chr] = qual\n"
            "   out = {}\n"
            "   for chr in dat:\n"
            "       if (chr not in dat) or dat[chr]<th: continue\n"
            "       out[chr] = 1\n"
            "   return out\n"
            "if args.chromosome!=None:\n"
            "   plot_chroms = { args.chromosome:1 }\n"
            "else:\n"
            "   plot_chroms = chroms_to_plot(args.plot_threshold)\n"
            "\n"
            "def read_dat(file,dat,plot_chr):\n"
            "   with open(file, 'r') as f:\n"
            "       reader = csv.reader(f, 'tab')\n"
            "       for row in reader:\n"
            "           chr = row[0]\n"
            "           if chr != plot_chr: continue\n"
            "           dat.append([int(row[1]), float(row[2]), float(row[3])])\n"
            "def read_cnv(file,cnv,plot_chr):\n"
            "   with open(file, 'r') as f:\n"
            "       reader = csv.reader(f, 'tab')\n"
            "       for row in reader:\n"
            "           chr = row[0]\n"
            "           if chr != plot_chr: continue\n"
            "           row[2] = int(row[2]) + 0.5\n"
            "           row[1] = int(row[1])\n"
            "           cnv.append(row[1:])\n"
            "def find_diffs(a,b):\n"
            "    out = []\n"
            "    diff = []\n"
            "    for i in range(len(a)):\n"
            "        if a[i][1]!=b[i][1]:\n"
            "            if i>0: diff.append([b[i-1][0],b[i-1][1],a[i-1][1]])\n"
            "            diff.append([b[i][0],b[i][1],a[i][1]])\n"
            "        elif len(diff):\n"
            "            diff.append([b[i][0],b[i][1],a[i][1]])\n"
            "            out.append(diff)\n"
            "            diff = []\n"
            "    if len(diff): out.append(diff)\n"
            "    return out\n"
            "\n"
            "for chr in sorted(plot_chroms.keys()):\n"
            "    control_dat = []\n"
            "    control_cnv = []\n"
            "    query_dat   = []\n"
            "    query_cnv   = []\n"
            "    read_dat('%s',control_dat,chr)\n"
            "    read_dat('%s',query_dat,chr)\n"
            "    read_cnv('%s',control_cnv,chr)\n"
            "    read_cnv('%s',query_cnv,chr)\n"
            "\n"
            "    fig,(ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(6,1,figsize=(10,8),sharex=True)\n"
            "    ax1.plot([x[0] for x in control_dat],[x[2] for x in control_dat], '.', ms=3,color='red')\n"
            "    ax2.plot([x[0] for x in control_dat],[x[1] for x in control_dat], '.', ms=3,color='red')\n"
            "    cn_dat = control_cnv\n"
            "    xgrid = [float(x[0]) for x in cn_dat]\n"
            "    ygrid = np.linspace(0,5,6)\n"
            "    xgrid, ygrid = np.meshgrid(xgrid, ygrid)\n"
            "    heat = np.zeros_like(xgrid)\n"
            "    for x in range(len(heat[0])-1):\n"
            "       heat[0][x] = cn_dat[x][2]\n"
            "       heat[1][x] = cn_dat[x][3]\n"
            "       heat[2][x] = cn_dat[x][4]\n"
            "       heat[3][x] = cn_dat[x][5]\n"
            "    mesh = ax3.pcolormesh(xgrid, ygrid, heat, cmap='bwr', shading='auto', alpha=0)\n"
            "    mesh.set_clim(vmin=-1,vmax=1)\n"
            "    ax3.plot([x[0] for x in cn_dat],[x[1] for x in cn_dat],'-',ms=3,color='black',lw=1.7)\n"
            "\n"
            "    ax6.plot([x[0] for x in query_dat],[x[2] for x in query_dat], '.', ms=3)\n"
            "    ax5.plot([x[0] for x in query_dat],[x[1] for x in query_dat], '.', ms=3)\n"
            "    cn_dat = query_cnv\n"
            "    xgrid = [float(x[0]) for x in cn_dat]\n"
            "    ygrid = np.linspace(0,5,6)\n"
            "    xgrid, ygrid = np.meshgrid(xgrid, ygrid)\n"
            "    heat = np.zeros_like(xgrid)\n"
            "    for x in range(len(heat[0])-1):\n"
            "       heat[0][x] = cn_dat[x][2]\n"
            "       heat[1][x] = cn_dat[x][3]\n"
            "       heat[2][x] = cn_dat[x][4]\n"
            "       heat[3][x] = cn_dat[x][5]\n"
            "    mesh = ax4.pcolormesh(xgrid, ygrid, heat, cmap='bwr_r')\n"
            "    mesh.set_clim(vmin=-1,vmax=1)\n"
            "    ax4.plot([x[0] for x in cn_dat],[x[1] for x in cn_dat],'-',ms=3,color='black',lw=1.7)\n"
            "    ax3.annotate(control_sample, xy=(0.02,0.1), xycoords='axes fraction', color='red',fontsize=12, va='bottom',ha='left')\n"
            "    ax4.annotate(query_sample, xy=(0.02,0.9), xycoords='axes fraction', color='blue',fontsize=12, va='top',ha='left')\n"
            "\n"
            "    diffs = find_diffs(control_cnv,query_cnv)\n"
            "    for diff in diffs:\n"
            "        ax3.plot([x[0] for x in diff],[x[1] for x in diff],'-',ms=3,color='blue',lw=1.7)\n"
            "        ax4.plot([x[0] for x in diff],[x[2] for x in diff],'-',ms=3,color='red',lw=1.7)\n"
            "\n"
            "    fig.suptitle('chr '+chr+', '+control_sample+' vs '+query_sample)\n"
            "    ax1.tick_params(axis='both', labelsize=8)\n"
            "    ax2.tick_params(axis='both', labelsize=8)\n"
            "    ax3.tick_params(axis='both', labelsize=8)\n"
            "    ax4.tick_params(axis='both', labelsize=8)\n"
            "    ax5.tick_params(axis='both', labelsize=8)\n"
            "    ax6.tick_params(axis='both', labelsize=8)\n"
            "    ax6.set_xlabel('Coordinate (chrom '+chr+')',fontsize=8)\n"
            "    ax1.set_ylabel('LRR')\n"
            "    ax2.set_ylabel('BAF')\n"
            "    ax3.set_ylabel('CN')\n"
            "    ax6.set_ylabel('LRR')\n"
            "    ax5.set_ylabel('BAF')\n"
            "    ax4.set_ylabel('CN')\n"
            "    ax3.set_ylim(-0.1,4.1)\n"
            "    ax3.set_yticks([0.5,1.5,2.5,3.5])\n"
            "    ax3.set_yticklabels(['CN0','CN1','CN2','CN3'])\n"
            "    ax4.set_ylim(-0.1,4.1)\n"
            "    ax4.set_yticks([0.5,1.5,2.5,3.5])\n"
            "    ax4.set_yticklabels(['CN0','CN1','CN2','CN3'])\n"
            "    plt.subplots_adjust(left=0.08,right=0.95,bottom=0.08,top=0.92,hspace=0)\n"
            "    plt.savefig('%s/plot.%s.%s.chr'+chr+'.png')\n"
            "    plt.close()\n"
            "\n", 
            args->control_sample.name,args->query_sample.name,
            args->output_dir,
            args->control_sample.dat_fname,args->query_sample.dat_fname,
            args->control_sample.cn_fname,args->query_sample.cn_fname,
            args->output_dir,args->control_sample.name,args->query_sample.name
        );
        fclose(fp);

    py_plot_cnv(fname,args->plot_th);
    free(fname);
}

static void destroy_data(args_t *args)
{
    bcf_sr_destroy(args->files);
    hmm_destroy(args->hmm);
    free(args->tmpf);
    free(args->sites);
    free(args->eprob);
    free(args->tprob);
    free(args->iprobs);
    free(args->summary_fname);
    free(args->nonref_afs);
    free(args->query_sample.baf);
    free(args->query_sample.lrr);
    free(args->control_sample.baf);
    free(args->control_sample.lrr);
    free(args->query_sample.name);
    free(args->query_sample.dat_fname);
    free(args->query_sample.cn_fname);
    free(args->query_sample.summary_fname);
    free(args->control_sample.dat_fname);
    free(args->control_sample.cn_fname);
    free(args->control_sample.summary_fname);
}

static inline char copy_number_state(args_t *args, int istate, int ismpl)
{
    char code[] = "01234";
    if ( !args->control_sample.name ) return code[istate];
    int idx = ismpl ? istate - (istate/N_STATES)*N_STATES : istate/N_STATES;
    return code[idx];
}

static double avg_ii_prob(int n, double *mat)
{
    int i;
    double avg = 0;
    for (i=0; i<n; i++) avg += MAT(mat,n,i,i);
    return avg/n;
}

#define GAUSS_CN1_PK_R(smpl)    (&((smpl)->gauss_param[0]))
#define GAUSS_CN1_PK_A(smpl)    (&((smpl)->gauss_param[1]))
#define GAUSS_CN2_PK_RR(smpl)   (&((smpl)->gauss_param[2]))
#define GAUSS_CN2_PK_RA(smpl)   (&((smpl)->gauss_param[3]))
#define GAUSS_CN2_PK_AA(smpl)   (&((smpl)->gauss_param[4]))
#define GAUSS_CN3_PK_RRR(smpl)  (&((smpl)->gauss_param[5]))
#define GAUSS_CN3_PK_RRA(smpl)  (&((smpl)->gauss_param[6]))
#define GAUSS_CN3_PK_RAA(smpl)  (&((smpl)->gauss_param[7]))
#define GAUSS_CN3_PK_AAA(smpl)  (&((smpl)->gauss_param[8]))

static inline double norm_prob(double baf, gauss_param_t *param)
{
    return exp(-(baf-param->mean)*(baf-param->mean)*0.5/param->dev2) / param->norm / sqrt(2*M_PI*param->dev2);
}

static int set_observed_prob(args_t *args, sample_t *smpl, int isite)
{
    float baf = smpl->baf[isite];
    float lrr = args->lrr_bias>0 ? smpl->lrr[isite] : 0;

    float fRR = args->fRR;
    float fRA = args->fRA;
    float fAA = args->fAA;

    if ( baf<0 )
    {
        // no call: either some technical issue or the call could not be made because it is CN0
        int i;
        smpl->pobs[CN0] = 0.5;
        for (i=1; i<N_STATES; i++) smpl->pobs[i] = (1.0-smpl->pobs[CN0])/(N_STATES-1);
        return 0;
    }

    double cn1_baf = 
        norm_prob(baf,GAUSS_CN1_PK_R(smpl)) * (fRR + fRA*0.5) +
        norm_prob(baf,GAUSS_CN1_PK_A(smpl)) * (fAA + fRA*0.5) ;
    double cn2_baf = 
        norm_prob(baf,GAUSS_CN2_PK_RR(smpl)) * fRR + 
        norm_prob(baf,GAUSS_CN2_PK_RA(smpl)) * fRA + 
        norm_prob(baf,GAUSS_CN2_PK_AA(smpl)) * fAA;
    double cn3_baf = 
        norm_prob(baf,GAUSS_CN3_PK_RRR(smpl)) * fRR + 
        norm_prob(baf,GAUSS_CN3_PK_RRA(smpl)) * fRA*0.5 + 
        norm_prob(baf,GAUSS_CN3_PK_RAA(smpl)) * fRA*0.5 + 
        norm_prob(baf,GAUSS_CN3_PK_AAA(smpl)) * fAA;

    double norm = cn1_baf + cn2_baf + cn3_baf;
    cn1_baf /= norm;
    cn2_baf /= norm;
    cn3_baf /= norm;

    #if DBG0
    if ( args->verbose ) fprintf(bcftools_stderr,"%f\t%f %f %f\n", baf,cn1_baf,cn2_baf,cn3_baf);
    #endif

    double cn1_lrr = exp(-(lrr + 0.45)*(lrr + 0.45)/smpl->lrr_dev2);
    double cn2_lrr = exp(-(lrr - 0.00)*(lrr - 0.00)/smpl->lrr_dev2);
    double cn3_lrr = exp(-(lrr - 0.30)*(lrr - 0.30)/smpl->lrr_dev2);

    smpl->pobs[CN0] = 0;
    smpl->pobs[CN1] = args->err_prob + (1 - args->baf_bias + args->baf_bias*cn1_baf)*(1 - args->lrr_bias + args->lrr_bias*cn1_lrr);
    smpl->pobs[CN2] = args->err_prob + (1 - args->baf_bias + args->baf_bias*cn2_baf)*(1 - args->lrr_bias + args->lrr_bias*cn2_lrr);
    smpl->pobs[CN3] = args->err_prob + (1 - args->baf_bias + args->baf_bias*cn3_baf)*(1 - args->lrr_bias + args->lrr_bias*cn3_lrr);

    return 0;
}

static void set_emission_prob(args_t *args, int isite)
{
    double *eprob = &args->eprob[args->nstates*isite];
    int i;
    for (i=0; i<N_STATES; i++)
        eprob[i] = args->query_sample.pobs[i];
}

static void set_emission_prob2(args_t *args, int isite)
{
    double *eprob = &args->eprob[args->nstates*isite];
    int i, j;
    for (i=0; i<N_STATES; i++)
    {
        for (j=0; j<N_STATES; j++)
        {
            eprob[i*N_STATES+j] = args->query_sample.pobs[i]*args->control_sample.pobs[j];
        }
    }
}

static void set_gauss_params(args_t *args, sample_t *smpl);
static double norm_cdf(double mean, double dev)
{
    double bot = 0, top = 1;
    top = 1 - 0.5*erfc((top-mean)/(dev*sqrt(2)));
    bot = 1 - 0.5*erfc((bot-mean)/(dev*sqrt(2)));
    return top-bot;
}

static void set_emission_probs(args_t *args)
{
    if ( !args->af_fname )
    {
        args->fRR = 0.76;
        args->fRA = 0.14;
        args->fAA = 0.098;
    }

    set_gauss_params(args, &args->query_sample);
    if ( args->control_sample.name ) set_gauss_params(args, &args->control_sample);

    #if DBG0
    args->verbose = 1;
    args->query_sample.baf[0] = 0; set_observed_prob(args,&args->query_sample,0);
    args->query_sample.baf[0] = 1/3.; set_observed_prob(args,&args->query_sample,0);
    args->query_sample.baf[0] = 1/2.; set_observed_prob(args,&args->query_sample,0);
    args->query_sample.baf[0] = 2/3.; set_observed_prob(args,&args->query_sample,0);
    args->query_sample.baf[0] = 1; set_observed_prob(args,&args->query_sample,0);
    args->verbose = 0;
    #endif

    int i;
    for (i=0; i<args->nsites; i++)
    {
        if ( args->af_fname )
        {
            args->fRR = (1-args->nonref_afs[i])*(1-args->nonref_afs[i]);
            args->fRA = 2*args->nonref_afs[i]*(1-args->nonref_afs[i]);
            args->fAA = args->nonref_afs[i]*args->nonref_afs[i];
        }
        set_observed_prob(args,&args->query_sample,i);
        if ( args->control_sample.name )
        {
            set_observed_prob(args,&args->control_sample,i);
            set_emission_prob2(args,i);
        }
        else
            set_emission_prob(args,i);
    }
}

static void smooth_data(float *dat, int ndat, int win)
{
    if ( win<=1 ) return;

    int i,j, k1 = win/2, k2 = win-k1;
    rbuf_t rbuf;
    rbuf_init(&rbuf,win);
    float sum = 0, *buf = (float*)malloc(sizeof(float)*win);
    for (i=0; i<k2; i++)
    {
        sum += dat[i];
        int j = rbuf_append(&rbuf);
        buf[j] = dat[i];
    }
    for (i=0; i<ndat; i++)
    {
        dat[i] = sum/rbuf.n;
        if ( i>=k1 )
        {
            j = rbuf_shift(&rbuf);
            sum -= buf[j];
        }
        if ( i+k2<ndat )
        {
            sum += dat[i+k2];
            j = rbuf_append(&rbuf);
            buf[j] = dat[i+k2];
        }
    }
    free(buf);
}

static void set_gauss_params(args_t *args, sample_t *smpl)
{
    int i;
    for (i=0; i<18; i++) smpl->gauss_param[i].dev2 = smpl->baf_dev2;

    double dev = sqrt(smpl->baf_dev2);

    GAUSS_CN1_PK_R(smpl)->mean = 0;
    GAUSS_CN1_PK_A(smpl)->mean = 1;
    GAUSS_CN1_PK_R(smpl)->norm = norm_cdf(GAUSS_CN1_PK_R(smpl)->mean,dev);
    GAUSS_CN1_PK_A(smpl)->norm = norm_cdf(GAUSS_CN1_PK_A(smpl)->mean,dev);

    GAUSS_CN2_PK_RR(smpl)->mean = 0;
    GAUSS_CN2_PK_RA(smpl)->mean = 0.5;
    GAUSS_CN2_PK_AA(smpl)->mean = 1;
    GAUSS_CN2_PK_RR(smpl)->norm = norm_cdf(GAUSS_CN2_PK_RR(smpl)->mean,dev);
    GAUSS_CN2_PK_RA(smpl)->norm = norm_cdf(GAUSS_CN2_PK_RA(smpl)->mean,dev);
    GAUSS_CN2_PK_AA(smpl)->norm = norm_cdf(GAUSS_CN2_PK_AA(smpl)->mean,dev);

    GAUSS_CN3_PK_RRR(smpl)->mean = 0;
    GAUSS_CN3_PK_RRA(smpl)->mean = 1.0/(2+smpl->cell_frac);
    GAUSS_CN3_PK_RAA(smpl)->mean = (1.0+smpl->cell_frac)/(2+smpl->cell_frac);
    GAUSS_CN3_PK_AAA(smpl)->mean = 1;
    GAUSS_CN3_PK_RRR(smpl)->norm = norm_cdf(GAUSS_CN3_PK_RRR(smpl)->mean,dev);
    GAUSS_CN3_PK_RRA(smpl)->norm = norm_cdf(GAUSS_CN3_PK_RRA(smpl)->mean,dev);
    GAUSS_CN3_PK_RAA(smpl)->norm = norm_cdf(GAUSS_CN3_PK_RAA(smpl)->mean,dev);
    GAUSS_CN3_PK_AAA(smpl)->norm = norm_cdf(GAUSS_CN3_PK_AAA(smpl)->mean,dev);
}

static int update_sample_args(args_t *args, sample_t *smpl, int ismpl)
{
    hmm_t *hmm = args->hmm;
    double *fwd = hmm_get_fwd_bwd_prob(hmm);
    int nstates = hmm_get_nstates(hmm);

    // estimate the BAF mean and deviation for CN3
    double mean_cn3 = 0, norm_cn3 = 0;
    double baf_dev2 = 0, baf_AA_dev2 = 0, norm_baf_AA_dev2 = 0;

    // experimental: smooth CN3 probs to bias toward bigger events, this lowers
    // the FP rate when the data is noisy
    hts_expand(float,args->nsites,args->mtmpf,args->tmpf);
    int i, j, k = 0;
    for (i=0; i<args->nsites; i++)
    {
        float baf = smpl->baf[i];
        if ( baf>4/5.) continue;        // skip AA genotypes
        if ( baf>0.5 ) baf = 1 - baf;   // the bands should be symmetric
        if ( baf<1/5.) continue;        // skip RR genotypes

        double prob_cn3 = 0, *probs = fwd + i*nstates;
        if ( !args->control_sample.name )
        {
            prob_cn3 = probs[CN3];
        }
        else if ( ismpl==0 )
        {
            // query sample: CN3 probability must be recovered from all states of the control sample
            for (j=0; j<N_STATES; j++) prob_cn3 += probs[CN3*N_STATES+j];
        }
        else
        {
            // same as above but for control sample
            for (j=0; j<N_STATES; j++) prob_cn3 += probs[CN3+j*N_STATES];
        }
        args->tmpf[k++] = prob_cn3;
    }
    smooth_data(args->tmpf, k, 50);
    k = 0;
    for (i=0; i<args->nsites; i++)
    {
        float baf = smpl->baf[i];
        if ( baf>4/5.) { baf_AA_dev2 += (1.0-baf)*(1.0-baf); norm_baf_AA_dev2++; continue; }       // skip AA genotypes
        if ( baf>0.5 ) baf = 1 - baf;   // the bands should be symmetric
        if ( baf<1/5.) continue;        // skip RR genotypes

        double prob_cn3 = args->tmpf[k++];
        mean_cn3 += prob_cn3 * baf;
        norm_cn3 += prob_cn3;
    }
    if ( !norm_cn3 )
    {
        smpl->cell_frac = 1.0;
        return 1;
    }
    mean_cn3 /= norm_cn3;
    k = 0;
    for (i=0; i<args->nsites; i++)
    {
        float baf = smpl->baf[i];
        if ( baf>0.5 ) baf = 1 - baf;   // the bands should be symmetric
        if ( baf<1/5.) continue;        // skip RR,AA genotypes

        double prob_cn3 = args->tmpf[k++];
        baf_dev2 += prob_cn3 * (baf - mean_cn3)*(baf - mean_cn3);
    }

    /*
        A noisy CN2 band is hard to distinguish from two CN3 bands which are
        close to each other. Set a treshold on the minimum separation based
        on the BAF deviation at p=0.95
    */
    baf_dev2 /= norm_cn3;
    baf_AA_dev2 /= norm_baf_AA_dev2;
    if ( baf_dev2 < baf_AA_dev2 )  baf_dev2 = baf_AA_dev2;
    double max_mean_cn3 = 0.5 - sqrt(baf_dev2)*1.644854;    // R: qnorm(0.95)=1.644854
    //fprintf(bcftools_stderr,"dev=%f  AA_dev=%f  max_mean_cn3=%f  mean_cn3=%f\n", baf_dev2,baf_AA_dev2,max_mean_cn3,mean_cn3);
    assert( max_mean_cn3>0 );

    double new_frac = 1./mean_cn3 - 2;
    if ( mean_cn3 > max_mean_cn3 || new_frac < args->optimize_frac )
    {
        // out of bounds, beyond our detection limits. Give up and say it converged
        smpl->cell_frac = 1.0;
        return 1;
    }
    if ( new_frac>1 ) new_frac = 1;
    int converged = fabs(new_frac - smpl->cell_frac) < 1e-1 ? 1 : 0;

    // Update dev2, but stay within safe limits
    if ( baf_dev2 > 3*smpl->baf_dev2_dflt ) baf_dev2 = 3*smpl->baf_dev2_dflt;
    else if ( baf_dev2 < 0.5*smpl->baf_dev2_dflt ) baf_dev2 = 0.5*smpl->baf_dev2_dflt;

    smpl->cell_frac = new_frac;
    smpl->baf_dev2  = baf_dev2;

    return converged;
}

// Update parameters which depend on the estimated fraction of aberrant cells
// in CN3.  Returns 0 if the current estimate did not need to be updated or 1
// if there was a change.
static int update_args(args_t *args)
{
    int converged = update_sample_args(args, &args->query_sample, 0);
    if ( args->control_sample.name )
    {
        converged += update_sample_args(args, &args->control_sample, 1);
        return converged==2 ? 0 : 1;
    }
    return converged ? 0 : 1;
}

// for an approximate estimate of the number of het genotypes in a region
#define BAF_LIKELY_HET(val)   (val)>0.25 && (val)<0.75

static void cnv_flush_viterbi(args_t *args)
{
    if ( !args->nsites ) return;

    // Set HMM transition matrix for the new chromsome again. This is for case
    // Baum-Welch was used, which is experimental, largerly unsupported and not
    // done by default.
    hmm_t *hmm = args->hmm;
    hmm_set_tprob(args->hmm, args->tprob, 10000);

    // Smooth LRR values to reduce noise
    if ( args->lrr_bias > 0 )
    {
        smooth_data(args->query_sample.lrr,args->nsites, args->lrr_smooth_win);
        if ( args->control_sample.name ) smooth_data(args->control_sample.lrr,args->nsites, args->lrr_smooth_win);
    }

    // Set the BAF peak likelihoods, such as P(RRR|CN3), taking account the
    // estimated fraction of aberrant cells in the mixture. With the new chromosome,
    // reset the fraction to the default value.
    args->query_sample.cell_frac   = args->query_sample.cell_frac_dflt;
    args->control_sample.cell_frac = args->control_sample.cell_frac_dflt;
    args->query_sample.baf_dev2    = args->query_sample.baf_dev2_dflt;
    args->control_sample.baf_dev2  = args->control_sample.baf_dev2_dflt;
    set_gauss_params(args, &args->query_sample);
    if ( args->control_sample.name ) set_gauss_params(args, &args->control_sample);

    if ( args->optimize_frac )
    {
        int niter = 0;
        fprintf(bcftools_stderr,"Attempting to estimate the fraction of aberrant cells (chr %s):\n", bcf_hdr_id2name(args->hdr,args->prev_rid));
        do
        {
            fprintf(bcftools_stderr,"\t.. %f %f", args->query_sample.cell_frac,args->query_sample.baf_dev2);
            if ( args->control_sample.name )
                fprintf(bcftools_stderr,"\t.. %f %f", args->control_sample.cell_frac,args->control_sample.baf_dev2);
            fprintf(bcftools_stderr,"\n");
            set_emission_probs(args);
            hmm_run_fwd_bwd(hmm, args->nsites, args->eprob, args->sites);
        }
        while ( update_args(args) && ++niter<20 );
        if ( niter>=20 )
        {
            // no convergence
            args->query_sample.cell_frac   = args->query_sample.cell_frac_dflt;
            args->control_sample.cell_frac = args->control_sample.cell_frac_dflt;
            args->query_sample.baf_dev2    = args->query_sample.baf_dev2_dflt;
            args->control_sample.baf_dev2  = args->control_sample.baf_dev2_dflt;
            set_gauss_params(args, &args->query_sample);
            if ( args->control_sample.name ) set_gauss_params(args, &args->control_sample);
        }

        fprintf(bcftools_stderr,"\t.. %f %f", args->query_sample.cell_frac,args->query_sample.baf_dev2);
        if ( args->control_sample.name )
            fprintf(bcftools_stderr,"\t.. %f %f", args->control_sample.cell_frac,args->control_sample.baf_dev2);
        fprintf(bcftools_stderr,"\n");

        fprintf(args->query_sample.summary_fh,"CF\t%s\t%d\t%d\t%.2f\t%f\n",
            bcf_hdr_id2name(args->hdr,args->prev_rid),args->sites[0]+1,args->sites[args->nsites-1]+1,
            args->query_sample.cell_frac,sqrt(args->query_sample.baf_dev2));
        if ( args->control_sample.name )
        {
            fprintf(args->control_sample.summary_fh,"CF\t%s\t%d\t%d\t%.2f\t%f\n",
                    bcf_hdr_id2name(args->hdr,args->prev_rid),args->sites[0]+1,args->sites[args->nsites-1]+1,
                    args->control_sample.cell_frac,sqrt(args->control_sample.baf_dev2));
            fprintf(args->summary_fh,"CF\t%s\t%d\t%d\t%.2f\t%.2f\t%f\t%f\n",
                    bcf_hdr_id2name(args->hdr,args->prev_rid),args->sites[0]+1,args->sites[args->nsites-1]+1,
                    args->query_sample.cell_frac, args->control_sample.cell_frac,
                    sqrt(args->query_sample.baf_dev2), sqrt(args->control_sample.baf_dev2));
        }
    }
    set_emission_probs(args);

    while ( args->baum_welch_th!=0 )
    {
        int nstates = hmm_get_nstates(hmm);
        double ori_ii = avg_ii_prob(nstates,hmm_get_tprob(hmm));
        hmm_run_baum_welch(hmm, args->nsites, args->eprob, args->sites);
        double new_ii = avg_ii_prob(nstates,hmm_get_tprob(hmm));
        fprintf(bcftools_stderr,"%e\t%e\t%e\n", ori_ii,new_ii,new_ii-ori_ii);
        double *tprob = init_tprob_matrix(nstates, 1-new_ii, args->same_prob);
        hmm_set_tprob(args->hmm, tprob, 10000);
        double *tprob_arr = hmm_get_tprob(hmm);
        free(tprob);
        if ( fabs(new_ii - ori_ii) < args->baum_welch_th )
        {
            int i,j;
            for (i=0; i<nstates; i++)
            {
                for (j=0; j<nstates; j++)
                {
                    fprintf(bcftools_stdout, " %.15f", MAT(tprob_arr,nstates,j,i));
                }
                fprintf(bcftools_stdout, "\n");
            }
            break;
        }
    }
    hmm_run_viterbi(hmm, args->nsites, args->eprob, args->sites);
    hmm_run_fwd_bwd(hmm, args->nsites, args->eprob, args->sites);


    // Output the results
    uint8_t *vpath = hmm_get_viterbi_path(hmm);
    double qual = 0, *fwd = hmm_get_fwd_bwd_prob(hmm);
    int i,j, isite, start_cn = vpath[0], start_pos = args->sites[0], istart_pos = 0;
    int ctrl_ntot = 0, smpl_ntot = 0, ctrl_nhet = 0, smpl_nhet = 0;
    for (isite=0; isite<args->nsites; isite++)
    {
        int state = vpath[args->nstates*isite];
        double *pval = fwd + isite*args->nstates;

        qual += pval[start_cn];

        // output CN and fwd-bwd likelihood for each site
        if ( args->query_sample.cn_fh )
        {
            fprintf(args->query_sample.cn_fh, "%s\t%d\t%c", bcf_hdr_id2name(args->hdr,args->prev_rid), args->sites[isite]+1, copy_number_state(args,state,0));
            if ( !args->control_sample.cn_fh )
                for (i=0; i<args->nstates; i++) fprintf(args->query_sample.cn_fh, "\t%f", pval[i]);
            else
                for (i=0; i<N_STATES; i++)
                {
                    double sum = 0;
                    for (j=0; j<N_STATES; j++) sum += pval[i*N_STATES+j];
                    fprintf(args->query_sample.cn_fh, "\t%f", sum);
                }
            fprintf(args->query_sample.cn_fh, "\n");
            if ( args->query_sample.baf[isite]>=0 )     // if non-missing
            {
                if ( BAF_LIKELY_HET(args->query_sample.baf[isite]) ) smpl_nhet++;
                smpl_ntot++;
            }
        }
        if ( args->control_sample.cn_fh )
        {
            fprintf(args->control_sample.cn_fh, "%s\t%d\t%c", bcf_hdr_id2name(args->hdr,args->prev_rid), args->sites[isite]+1, copy_number_state(args,state,1));
            for (i=0; i<N_STATES; i++)
            {
                double sum = 0;
                for (j=0; j<N_STATES; j++) sum += pval[i+N_STATES*j];
                fprintf(args->control_sample.cn_fh, "\t%f", sum);
            }
            fprintf(args->control_sample.cn_fh, "\n");
            if ( args->control_sample.baf[isite]>=0 )     // if non-missing
            {
                if ( BAF_LIKELY_HET(args->control_sample.baf[isite]) ) ctrl_nhet++;
                ctrl_ntot++;
            }
        }

        if ( start_cn != state )
        {
            char start_cn_query = copy_number_state(args,start_cn,0);
            qual = phred_score(1 - qual/(isite - istart_pos));
            fprintf(args->query_sample.summary_fh,"RG\t%s\t%d\t%d\t%c\t%.1f\t%d\t%d\n",
                bcf_hdr_id2name(args->hdr,args->prev_rid), start_pos+1, args->sites[isite],start_cn_query,qual,smpl_ntot,smpl_nhet);

            if ( args->control_sample.name )
            {
                // regions 0-based, half-open
                char start_cn_ctrl = copy_number_state(args,start_cn,1);
                fprintf(args->control_sample.summary_fh,"RG\t%s\t%d\t%d\t%c\t%.1f\t%d\t%d\n",
                    bcf_hdr_id2name(args->hdr,args->prev_rid), start_pos+1, args->sites[isite],start_cn_ctrl,qual,ctrl_ntot,ctrl_nhet);
                fprintf(args->summary_fh,"RG\t%s\t%d\t%d\t%c\t%c\t%.1f\t%d\t%d\t%d\t%d\n",
                    bcf_hdr_id2name(args->hdr,args->prev_rid), start_pos+1, args->sites[isite],start_cn_query,start_cn_ctrl,qual,smpl_ntot,smpl_nhet,ctrl_ntot,ctrl_nhet);
            }

            istart_pos = isite;
            start_pos = args->sites[isite];
            start_cn = state;
            qual = 0;
            smpl_ntot = smpl_nhet = ctrl_ntot = ctrl_nhet = 0;
        }
    }
    qual = phred_score(1 - qual/(isite - istart_pos));
    char start_cn_query = copy_number_state(args,start_cn,0);
    fprintf(args->query_sample.summary_fh,"RG\t%s\t%d\t%d\t%c\t%.1f\t%d\t%d\n",
        bcf_hdr_id2name(args->hdr,args->prev_rid), start_pos+1, args->sites[isite-1]+1,start_cn_query,qual,smpl_ntot,smpl_nhet);
    if ( args->control_sample.name )
    {
        char start_cn_ctrl = copy_number_state(args,start_cn,1);
        fprintf(args->control_sample.summary_fh,"RG\t%s\t%d\t%d\t%c\t%.1f\t%d\t%d\n",
            bcf_hdr_id2name(args->hdr,args->prev_rid), start_pos+1, args->sites[isite-1]+1,start_cn_ctrl,qual,ctrl_ntot,ctrl_nhet);
        fprintf(args->summary_fh,"RG\t%s\t%d\t%d\t%c\t%c\t%.1f\t%d\t%d\t%d\t%d\n",
            bcf_hdr_id2name(args->hdr,args->prev_rid), start_pos+1, args->sites[isite-1]+1,start_cn_query,start_cn_ctrl,qual,smpl_ntot,smpl_nhet,ctrl_ntot,ctrl_nhet);
    }
}

static int parse_lrr_baf(sample_t *smpl, bcf_fmt_t *baf_fmt, bcf_fmt_t *lrr_fmt, float *baf, float *lrr)
{
    *baf = ((float*)(baf_fmt->p + baf_fmt->size*smpl->idx))[0];
    if ( bcf_float_is_missing(*baf) || isnan(*baf) ) *baf = -0.1;    // arbitrary negative value == missing value

    if ( lrr_fmt )
    {
        *lrr = ((float*)(lrr_fmt->p + lrr_fmt->size*smpl->idx))[0];
        if ( bcf_float_is_missing(*lrr) || isnan(*lrr) ) { *lrr = 0; *baf = -0.1; }
    }
    else
        *lrr = 0;

    return *baf<0 ? 0 : 1;
}

static void cnv_next_line(args_t *args, bcf1_t *line)
{
    if ( !line ) 
    {
        // Done, flush viterbi
        cnv_flush_viterbi(args);
        return;
    }

    if ( line->rid!=args->prev_rid )
    {
        // New chromosome
        cnv_flush_viterbi(args);
        args->prev_rid = line->rid;
        args->nsites = 0;
        args->nRR = args->nAA = args->nRA = 0;
    }

    // Process line
    args->ntot++;

    bcf_fmt_t *baf_fmt, *lrr_fmt = NULL;
    if ( !(baf_fmt = bcf_get_fmt(args->hdr, line, "BAF")) ) return; 
    if ( args->lrr_bias>0 && !(lrr_fmt = bcf_get_fmt(args->hdr, line, "LRR")) ) return;

    float baf1,lrr1,baf2,lrr2;
    int ret = 0;
    ret += parse_lrr_baf(&args->query_sample,  baf_fmt,lrr_fmt,&baf1,&lrr1);
    ret += parse_lrr_baf(&args->control_sample,baf_fmt,lrr_fmt,&baf2,&lrr2);
    if ( !ret ) return;

    // Realloc buffers needed to store observed data and used by viterbi and fwd-bwd
    args->nsites++;
    int m = args->msites;
    hts_expand(uint32_t,args->nsites,args->msites,args->sites);
    if ( args->msites!=m )
    {
        args->eprob = (double*) realloc(args->eprob,sizeof(double)*args->msites*args->nstates);
        if ( args->control_sample.name )
        {
            args->control_sample.lrr = (float*) realloc(args->control_sample.lrr,sizeof(float)*args->msites);
            args->control_sample.baf = (float*) realloc(args->control_sample.baf,sizeof(float)*args->msites);
        }
        args->query_sample.lrr = (float*) realloc(args->query_sample.lrr,sizeof(float)*args->msites);
        args->query_sample.baf = (float*) realloc(args->query_sample.baf,sizeof(float)*args->msites);
        if ( args->af_fname )
            args->nonref_afs = (float*) realloc(args->nonref_afs,sizeof(float)*args->msites);
    }
    args->sites[args->nsites-1] = line->pos;
    args->query_sample.lrr[args->nsites-1] = lrr1;
    args->query_sample.baf[args->nsites-1] = baf1;
    if ( args->af_fname )
    {
        double alt_freq;
        args->nonref_afs[args->nsites-1] = read_AF(args->files->targets,line,&alt_freq)<0 ? args->nonref_af_dflt : alt_freq;
    }
    if ( args->control_sample.name )
    {
        args->control_sample.lrr[args->nsites-1] = lrr2;
        args->control_sample.baf[args->nsites-1] = baf2;
        if ( baf2>=0 )  // skip missing values
            fprintf(args->control_sample.dat_fh,"%s\t%"PRId64"\t%.3f\t%.3f\n",bcf_hdr_id2name(args->hdr,args->prev_rid), (int64_t) line->pos+1,baf2,lrr2);
    }
    if ( baf1>=0 )  // skip missing values
        fprintf(args->query_sample.dat_fh,"%s\t%"PRId64"\t%.3f\t%.3f\n",bcf_hdr_id2name(args->hdr,args->prev_rid), (int64_t) line->pos+1,baf1,lrr1);

    if ( baf1>=0 )
    {
        if ( baf1<1/5. ) args->nRR++;
        else if ( baf1>4/5. ) args->nAA++;
        else args->nRA++;
    }
    args->nused++;
}

static void usage(args_t *args)
{
    fprintf(bcftools_stderr, "\n");
    fprintf(bcftools_stderr, "About:   Copy number variation caller, requires Illumina's B-allele frequency (BAF) and Log R\n");
    fprintf(bcftools_stderr, "         Ratio intensity (LRR). The HMM considers the following copy number states: CN 2\n");
    fprintf(bcftools_stderr, "         (normal), 1 (single-copy loss), 0 (complete loss), 3 (single-copy gain)\n");
    fprintf(bcftools_stderr, "Usage:   bcftools cnv [OPTIONS] FILE.vcf\n");
    fprintf(bcftools_stderr, "General Options:\n");
    fprintf(bcftools_stderr, "    -c, --control-sample STRING      Optional control sample name to highlight differences\n");
    fprintf(bcftools_stderr, "    -f, --AF-file FILE               Read allele frequencies from file (CHR\\tPOS\\tREF,ALT\\tAF)\n");
    fprintf(bcftools_stderr, "    -o, --output-dir PATH            \n");
    fprintf(bcftools_stderr, "    -p, --plot-threshold FLOAT       Plot aberrant chromosomes with quality at least FLOAT\n");
    fprintf(bcftools_stderr, "    -r, --regions REGION             Restrict to comma-separated list of regions\n");
    fprintf(bcftools_stderr, "    -R, --regions-file FILE          Restrict to regions listed in a file\n");
    fprintf(bcftools_stderr, "        --regions-overlap 0|1|2      Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n");
    fprintf(bcftools_stderr, "    -s, --query-sample STRING        Query samply name\n");
    fprintf(bcftools_stderr, "    -t, --targets REGION             Similar to -r but streams rather than index-jumps\n");
    fprintf(bcftools_stderr, "    -T, --targets-file FILE          Similar to -R but streams rather than index-jumps\n");
    fprintf(bcftools_stderr, "        --targets-overlap 0|1|2      Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n");
    fprintf(bcftools_stderr, "HMM Options:\n");
    fprintf(bcftools_stderr, "    -a, --aberrant FLOAT[,FLOAT]     Fraction of aberrant cells in query and control [1.0,1.0]\n");
    fprintf(bcftools_stderr, "    -b, --BAF-weight FLOAT           Relative contribution from BAF [1]\n");
    fprintf(bcftools_stderr, "    -d, --BAF-dev FLOAT[,FLOAT]      Expected BAF deviation in query and control [0.04,0.04]\n"); // experimental
    fprintf(bcftools_stderr, "    -e, --err-prob FLOAT             Uniform error probability [1e-4]\n");
    fprintf(bcftools_stderr, "    -k, --LRR-dev FLOAT[,FLOAT]      Expected LRR deviation [0.2,0.2]\n"); // experimental
    fprintf(bcftools_stderr, "    -l, --LRR-weight FLOAT           Relative contribution from LRR [0.2]\n");
    fprintf(bcftools_stderr, "    -L, --LRR-smooth-win INT         Window of LRR moving average smoothing [10]\n");
    fprintf(bcftools_stderr, "    -O, --optimize FLOAT             Estimate fraction of aberrant cells down to FLOAT [1.0]\n");
    fprintf(bcftools_stderr, "    -P, --same-prob FLOA>            Prior probability of -s/-c being the same [0.5]\n");
    fprintf(bcftools_stderr, "    -x, --xy-prob FLOAT              P(x|y) transition probability [1e-9]\n");
    fprintf(bcftools_stderr, "\n");
    bcftools_exit(1);
}

int main_vcfcnv(int argc, char *argv[])
{
    int c;
    args_t *args    = (args_t*) calloc(1,sizeof(args_t));
    args->argc      = argc; args->argv = argv;
    args->files     = bcf_sr_init();
    args->plot_th   = 1e9;   // by default plot none
    args->nonref_af_dflt = 0.1;
    args->lrr_smooth_win = 10;

    args->query_sample.cell_frac_dflt = 1;
    args->control_sample.cell_frac_dflt = 1;

    // How much FORMAT/LRR and FORMAT/BAF matter
    args->lrr_bias  = 0.2;
    args->baf_bias  = 1.0;
    args->err_prob  = 1e-4;

    // Transition probability to a different state and the prior of both samples being the same
    args->ij_prob   = 1e-9;
    args->same_prob = 0.5;

    // Squared std dev of BAF and LRR values (gaussian noise), estimated from real data (hets, one sample, one chr)
    args->query_sample.baf_dev2_dflt = args->control_sample.baf_dev2_dflt = 0.04*0.04; // illumina: 0.03
    args->query_sample.lrr_dev2 = args->control_sample.lrr_dev2 = 0.2*0.2; //0.20*0.20;   // illumina: 0.18

    int regions_is_file = 0, targets_is_file = 0;
    int regions_overlap = 1;
    int targets_overlap = 0;

    static struct option loptions[] = 
    {
        {"BAF-dev",1,0,'d'},
        {"LRR-dev",1,0,'k'},
        {"LRR-smooth-win",1,0,'L'},
        {"AF-file",1,0,'f'},
        {"baum-welch",1,0,'W'}, // hidden
        {"optimize",1,0,'O'},
        {"aberrant",1,0,'a'},
        {"err-prob",1,0,'e'},
        {"BAF-weight",1,0,'b'},
        {"LRR-weight",1,0,'l'},
        {"same-prob",1,0,'P'},
        {"xy-prob",1,0,'x'},
        {"query-sample",1,0,'s'},
        {"control-sample",1,0,'c'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"targets-overlap",required_argument,NULL,4},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"regions-overlap",required_argument,NULL,3},
        {"plot-threshold",1,0,'p'},
        {"output-dir",1,0,'o'},
        {0,0,0,0}
    };
    char *tmp = NULL;
    while ((c = getopt_long(argc, argv, "h?r:R:t:T:s:o:p:l:T:c:b:P:x:e:O:W:f:a:L:d:k:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'L': 
                args->lrr_smooth_win = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse: --LRR-smooth-win %s\n", optarg);
                break;
            case 'f': args->af_fname = optarg; break;
            case 'O': 
                args->optimize_frac = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -O %s\n", optarg);
                break;
            case 'd':
                args->query_sample.baf_dev2_dflt = strtod(optarg,&tmp);
                if ( *tmp )
                {
                    if ( *tmp!=',') error("Could not parse: -d %s\n", optarg);
                    args->control_sample.baf_dev2_dflt = strtod(tmp+1,&tmp);
                    if ( *tmp ) error("Could not parse: -d %s\n", optarg);
                }
                else
                    args->control_sample.baf_dev2_dflt = args->query_sample.baf_dev2_dflt;
                args->query_sample.baf_dev2_dflt   *= args->query_sample.baf_dev2_dflt;
                args->control_sample.baf_dev2_dflt *= args->control_sample.baf_dev2_dflt;
                break;
            case 'k':
                args->query_sample.lrr_dev2 = strtod(optarg,&tmp);
                if ( *tmp )
                {
                    if ( *tmp!=',') error("Could not parse: -k %s\n", optarg);
                    args->control_sample.lrr_dev2 = strtod(tmp+1,&tmp);
                    if ( *tmp ) error("Could not parse: -d %s\n", optarg);
                }
                else
                    args->control_sample.lrr_dev2 = args->query_sample.lrr_dev2;
                args->query_sample.lrr_dev2   *= args->query_sample.lrr_dev2;
                args->control_sample.lrr_dev2 *= args->control_sample.lrr_dev2;
                break;
            case 'a':
                args->query_sample.cell_frac_dflt = strtod(optarg,&tmp);
                if ( *tmp )
                {
                    if ( *tmp!=',') error("Could not parse: -a %s\n", optarg);
                    args->control_sample.cell_frac_dflt = strtod(tmp+1,&tmp);
                    if ( *tmp ) error("Could not parse: -a %s\n", optarg);
                }
                break;
            case 'W':
                args->baum_welch_th = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -W %s\n", optarg);
                break;
            case 'e': 
                args->err_prob = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -e %s\n", optarg);
                break;
            case 'b': 
                args->baf_bias = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -b %s\n", optarg);
                break;
            case 'x': 
                args->ij_prob = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -x %s\n", optarg);
                break;
            case 'P': 
                args->same_prob = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -P %s\n", optarg);
                break;
            case 'l': 
                args->lrr_bias = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -l %s\n", optarg);
                break;
            case 'p': 
                args->plot_th = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -p %s\n", optarg);
                break;
            case 'o': args->output_dir = optarg; break;
            case 's': args->query_sample.name = strdup(optarg); break;
            case 'c': args->control_sample.name = optarg; break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case  3 :
                regions_overlap = parse_overlap_option(optarg);
                if ( regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case  4 :
                targets_overlap = parse_overlap_option(optarg);
                if ( targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                break;
            case 'h': 
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";
    }
    else fname = argv[optind];
    if ( !fname ) usage(args);

    if ( !args->output_dir ) error("Expected -o option\n");
    if ( args->regions_list )
    {
        bcf_sr_set_opt(args->files,BCF_SR_REGIONS_OVERLAP,regions_overlap);
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->targets_list )
    {
        bcf_sr_set_opt(args->files,BCF_SR_TARGETS_OVERLAP,targets_overlap);
        if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }
    if ( args->af_fname )
    {
        if ( bcf_sr_set_targets(args->files, args->af_fname, 1, 3)<0 )
            error("Failed to read the targets: %s\n", args->af_fname);
    }
    if ( !bcf_sr_add_reader(args->files, fname) )
        error("Failed to read from %s: %s\n", !strcmp("-",fname)?"standard input":fname,bcf_sr_strerror(args->files->errnum));
    
    init_data(args);
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        cnv_next_line(args, line);
    }
    cnv_next_line(args, NULL);
    create_plots(args);
    fprintf(bcftools_stderr,"Number of lines: total/processed: %d/%d\n", args->ntot,args->nused);
    destroy_data(args);
    free(args);
    return 0;
}


