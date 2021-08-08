#include "bcftools.pysam.h"

/*  vcfsom.c -- SOM (Self-Organizing Map) filtering.

    Copyright (C) 2013-2014, 2020 Genome Research Ltd.

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
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <time.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <htslib/hts_os.h>
#include <inttypes.h>
#include "bcftools.h"

#define SOM_TRAIN    1
#define SOM_CLASSIFY 2

typedef struct
{
    int ndim;       // dimension of the map (2D, 3D, ...)
    int nbin;       // number of bins in th map
    int size;       // pow(nbin,ndim)
    int kdim;       // dimension of the input vectors
    int nt, t;      // total number of learning cycles and the current cycle
    double *w, *c;  // weights and counts (sum of learning influence)
    double learn;   // learning rate
    double bmu_th;  // best-matching unit threshold
    int *a_idx, *b_idx; // temp arrays for traversing variable number of nested loops
    double *div;        // dtto
}
som_t;

typedef struct
{
    // SOM parameters
    double bmu_th, learn;
    int ndim, nbin, ntrain, t;
    int nfold;                  // n-fold cross validation = the number of SOMs
    som_t **som;

    // annots reader's data
    htsFile *file;              // reader
    kstring_t str;              // temporary string for the reader
    int dclass, mvals;
    double *vals;

    // training data
    double *train_dat;
    int *train_class, mtrain_class, mtrain_dat;

    int rand_seed, good_class, bad_class;
    char **argv, *fname, *prefix;
    int argc, action, train_bad, merge;
}
args_t;

static void usage(void);
FILE *open_file(char **fname, const char *mode, const char *fmt, ...);
void mkdir_p(const char *fmt, ...);

char *msprintf(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(NULL, 0, fmt, ap) + 2;
    va_end(ap);

    char *str = (char*)malloc(n);
    va_start(ap, fmt);
    vsnprintf(str, n, fmt, ap);
    va_end(ap);

    return str;
}

/*
 *  char *t, *p = str;
 *  t = column_next(p, '\t');
 *  if ( strlen("<something>")==t-p && !strncmp(p,"<something>",t-p) ) fprintf(bcftools_stdout, "found!\n");
 *
 *  char *t;
 *  t = column_next(str, '\t'); if ( !*t ) error("expected field\n", str);
 *  t = column_next(t+1, '\t'); if ( !*t ) error("expected field\n", str);
 */
static inline char *column_next(char *start, char delim)
{
    char *end = start;
    while (*end && *end!=delim) end++;
    return end;
}
/**
 *  annots_reader_next() - reads next line from annots.tab.gz and sets: class, vals
 *   Returns 1 on successful read or 0 if no further record could be read.
 */
int annots_reader_next(args_t *args)
{
    args->str.l = 0;
    if ( hts_getline(args->file,'\n',&args->str)<=0 ) return 0;

    char *t, *line = args->str.s;

    if ( !args->mvals )
    {
        t = line;
        while ( *t )
        {
            if ( *t=='\t' ) args->mvals++;
            t++;
        }
        args->vals = (double*) malloc(args->mvals*sizeof(double));
    }

    // class
    args->dclass = atoi(line);
    t = column_next(line, '\t');

    // values
    int i;
    for (i=0; i<args->mvals; i++)
    {
        if ( !*t ) error("Could not parse %d-th data field: is the line truncated?\nThe line was: [%s]\n",i+2,line);
        args->vals[i] = atof(++t);
        t = column_next(t,'\t');
    }
    return 1;
}
void annots_reader_reset(args_t *args)
{
    if ( args->file ) hts_close(args->file);
    if ( !args->fname ) error("annots_reader_reset: no fname\n");
    args->file = hts_open(args->fname, "r");
}
void annots_reader_close(args_t *args)
{
    hts_close(args->file);
}

static void som_write_map(char *prefix, som_t **som, int nsom)
{
    FILE *fp = open_file(NULL,"w","%s.som",prefix);
    size_t nw;
    if ( (nw=fwrite("SOMv1",5,1,fp))!=5 ) error("Failed to write 5 bytes\n");
    if ( (nw=fwrite(&nsom,sizeof(int),1,fp))!=sizeof(int) ) error("Failed to write %zu bytes\n",sizeof(int));
    int i;
    for (i=0; i<nsom; i++)
    {
        if ( (nw=fwrite(&som[i]->size,sizeof(int),1,fp))!=sizeof(int) ) error("Failed to write %zu bytes\n",sizeof(int));
        if ( (nw=fwrite(&som[i]->kdim,sizeof(int),1,fp))!=sizeof(int) ) error("Failed to write %zu bytes\n",sizeof(int));
        if ( (nw=fwrite(som[i]->w,sizeof(double),som[i]->size*som[i]->kdim,fp))!=sizeof(double)*som[i]->size*som[i]->kdim ) error("Failed to write %zu bytes\n",sizeof(double)*som[i]->size*som[i]->kdim);
        if ( (nw=fwrite(som[i]->c,sizeof(double),som[i]->size,fp))!=sizeof(double)*som[i]->size ) error("Failed to write %zu bytes\n",sizeof(double)*som[i]->size);
    }
    if ( fclose(fp) ) error("%s.som: fclose failed\n",prefix);
}
static som_t** som_load_map(char *prefix, int *nsom)
{
    FILE *fp = open_file(NULL,"r","%s.som",prefix);
    char buf[5];
    if ( fread(buf,5,1,fp)!=1 || strncmp(buf,"SOMv1",5) ) error("Could not parse %s.som\n", prefix);

    if ( fread(nsom,sizeof(int),1,fp)!=1 ) error("Could not read %s.som\n", prefix);
    som_t **som = (som_t**)malloc(*nsom*sizeof(som_t*));

    int i;
    for (i=0; i<*nsom; i++)
    {
        som[i] = (som_t*) calloc(1,sizeof(som_t));
        if ( fread(&som[i]->size,sizeof(int),1,fp) != 1 ) error("Could not read %s.som\n", prefix);
        if ( fread(&som[i]->kdim,sizeof(int),1,fp) != 1 ) error("Could not read %s.som\n", prefix);
        som[i]->w = (double*) malloc(sizeof(double)*som[i]->size*som[i]->kdim);
        som[i]->c = (double*) malloc(sizeof(double)*som[i]->size);
        if ( fread(som[i]->w,sizeof(double),som[i]->size*som[i]->kdim,fp) != som[i]->size*som[i]->kdim ) error("Could not read from %s.som\n", prefix);
        if ( fread(som[i]->c,sizeof(double),som[i]->size,fp) != som[i]->size ) error("Could not read from %s.som\n", prefix);
    }
    if ( fclose(fp) ) error("%s.som: fclose failed\n",prefix);
    return som;
}
static void som_create_plot(som_t *som, char *prefix)
{
    if ( som->ndim!=2 ) return;

    char *fname;
    FILE *fp = open_file(&fname,"w","%s.py",prefix);
    fprintf(fp,
            "import matplotlib as mpl\n"
            "mpl.use('Agg')\n"
            "import matplotlib.pyplot as plt\n"
            "\n"
            "dat = [\n"
           );
    int i,j;
    double *val = som->c;
    for (i=0; i<som->nbin; i++)
    {
        fprintf(fp,"[");
        for (j=0; j<som->nbin; j++)
        {
            if ( j>0 ) fprintf(fp,",");
            fprintf(fp,"%e", *val);
            val++;
        }
        fprintf(fp,"],\n");
    }
    fprintf(fp,
            "]\n"
            "fig = plt.figure()\n"
            "ax1 = plt.subplot(111)\n"
            "im1 = ax1.imshow(dat)\n"
            "fig.colorbar(im1)\n"
            "plt.savefig('%s.png')\n"
            "plt.close()\n"
            "\n", prefix
           );
    fclose(fp);
    free(fname);
}
// Find the best matching unit: the node with minimum distance from the input vector
static inline int som_find_bmu(som_t *som, double *vec, double *dist)
{
    double *ptr = som->w;
    double min_dist = HUGE_VAL;
    int min_idx = 0;

    int i, k;
    for (i=0; i<som->size; i++)
    {
        double dist = 0;
        for (k=0; k<som->kdim; k++)
            dist += (vec[k] - ptr[k]) * (vec[k] - ptr[k]);
        if ( dist < min_dist )
        {
            min_dist = dist;
            min_idx  = i;
        }
        ptr += som->kdim;
    }

    if ( dist ) *dist = min_dist;
    return min_idx;
}
static inline double som_get_score(som_t *som, double *vec, double bmu_th)
{
    double *ptr = som->w;
    double min_dist = HUGE_VAL;

    int i, k;
    for (i=0; i<som->size; i++)
    {
        if ( som->c[i] >= bmu_th )
        {
            double dist = 0;
            for (k=0; k<som->kdim; k++)
                dist += (vec[k] - ptr[k]) * (vec[k] - ptr[k]);
            if ( dist < min_dist ) min_dist = dist;
        }
        ptr += som->kdim;
    }
    return sqrt(min_dist);
}
// Convert flat index to that of a k-dimensional cube
static inline void som_idx_to_ndim(som_t *som, int idx, int *ndim)
{
    int i;
    double sub = 0;

    ndim[0] = idx/som->div[0];
    for (i=1; i<som->ndim; i++)
    {
        sub += ndim[i-1] * som->div[i-1];
        ndim[i] = (idx - sub)/som->div[i];
    }
}
static void som_train_site(som_t *som, double *vec, int update_counts)
{
    // update learning rate and learning radius
    som->t++;
    double dt = exp(-som->t/som->nt);
    double learning_rate = som->learn * dt;
    double radius = som->nbin * dt; radius *= radius;

    // find the best matching unit and its indexes
    int min_idx = som_find_bmu(som, vec, NULL);
    som_idx_to_ndim(som, min_idx, som->a_idx);

    // update the weights: traverse the map and make all nodes within the
    // radius more similar to the input vector
    double *ptr = som->w;
    double *cnt = som->c;
    int i, j, k;
    for (i=0; i<som->size; i++)
    {
        som_idx_to_ndim(som, i, som->b_idx);
        double dist = 0;
        for (j=0; j<som->ndim; j++)
            dist += (som->a_idx[j] - som->b_idx[j]) * (som->a_idx[j] - som->b_idx[j]);
        if ( dist <= radius )
        {
            double influence = exp(-dist*dist*0.5/radius) * learning_rate;
            for (k=0; k<som->kdim; k++)
                ptr[k] += influence * (vec[k] - ptr[k]);

            // Bad sites may help to shape the map, but only nodes with big enough
            // influence will be used for classification.
            if ( update_counts ) *cnt += influence;
        }
        ptr += som->kdim;
        cnt++;
    }
}
static void som_norm_counts(som_t *som)
{
    int i;
    double max = 0;
    for (i=0; i<som->size; i++)
        if ( max < som->c[i] ) max = som->c[i];
    for (i=0; i<som->size; i++)
        som->c[i] /= max;
}
static som_t *som_init(args_t *args)
{
    som_t *som  = (som_t*) calloc(1,sizeof(som_t));
    som->ndim   = args->ndim;
    som->nbin   = args->nbin;
    som->kdim   = args->mvals;
    som->nt     = args->ntrain;
    som->learn  = args->learn;
    som->bmu_th = args->bmu_th;
    som->size   = pow(som->nbin,som->ndim);
    som->w = (double*) malloc(sizeof(double)*som->size*som->kdim);
    if ( !som->w ) error("Could not alloc %"PRIu64" bytes [nbin=%d ndim=%d kdim=%d]\n", (uint64_t)(sizeof(double)*som->size*som->kdim),som->nbin,som->ndim,som->kdim);
    som->c = (double*) calloc(som->size,sizeof(double));
    if ( !som->w ) error("Could not alloc %"PRIu64" bytes [nbin=%d ndim=%d]\n", (uint64_t)(sizeof(double)*som->size),som->nbin,som->ndim);
    int i;
    for (i=0; i<som->size*som->kdim; i++)
        som->w[i] = random();
    som->a_idx = (int*) malloc(sizeof(int)*som->ndim);
    som->b_idx = (int*) malloc(sizeof(int)*som->ndim);
    som->div   = (double*) malloc(sizeof(double)*som->ndim);
    for (i=0; i<som->ndim; i++)
        som->div[i] = pow(som->nbin,som->ndim-i-1);
    return som;
}
static void som_destroy(som_t *som)
{
    free(som->a_idx); free(som->b_idx); free(som->div);
    free(som->w); free(som->c);
    free(som);
}

static void init_data(args_t *args)
{
    // Get first line to learn the vector size
    annots_reader_reset(args);
    annots_reader_next(args);

    if ( args->action==SOM_CLASSIFY )
        args->som = som_load_map(args->prefix,&args->nfold);
}
static void destroy_data(args_t *args)
{
    int i;
    if ( args->som )
    {
        for (i=0; i<args->nfold; i++) som_destroy(args->som[i]);
    }
    free(args->train_dat);
    free(args->train_class);
    free(args->som);
    free(args->vals);
    free(args->str.s);
}

#define MERGE_MIN 0
#define MERGE_MAX 1
#define MERGE_AVG 2
static double get_min_score(args_t *args, int iskip)
{
    int i;
    double score, min_score = HUGE_VAL;
    for (i=0; i<args->nfold; i++)
    {
        if ( i==iskip ) continue;
        score = som_get_score(args->som[i], args->vals, args->bmu_th);
        if ( i==0 || score < min_score ) min_score = score;
    }
    return min_score;
}
static double get_max_score(args_t *args, int iskip)
{
    int i;
    double score, max_score = -HUGE_VAL;
    for (i=0; i<args->nfold; i++)
    {
        if ( i==iskip ) continue;
        score = som_get_score(args->som[i], args->vals, args->bmu_th);
        if ( i==0 || max_score < score ) max_score = score;
    }
    return max_score;
}
static double get_avg_score(args_t *args, int iskip)
{
    int i, n = 0;
    double score = 0;
    for (i=0; i<args->nfold; i++)
    {
        if ( i==iskip ) continue;
        score += som_get_score(args->som[i], args->vals, args->bmu_th);
        n++;
    }
    return score/n;
}
static int cmpfloat_desc(const void *a, const void *b)
{
    float fa = *((float*)a);
    float fb = *((float*)b);
    if ( fa<fb ) return 1;
    if ( fa>fb ) return -1;
    return 0;
}

static void create_eval_plot(args_t *args)
{
    FILE *fp = open_file(NULL,"w","%s.eval.py", args->prefix);
    fprintf(fp,
            "import matplotlib as mpl\n"
            "mpl.use('Agg')\n"
            "import matplotlib.pyplot as plt\n"
            "\n"
            "import csv\n"
            "csv.register_dialect('tab', delimiter='\\t', quoting=csv.QUOTE_NONE)\n"
            "dat = []\n"
            "with open('%s.eval', 'r') as f:\n"
            "\treader = csv.reader(f, 'tab')\n"
            "\tfor row in reader:\n"
            "\t\tif row[0][0]!='#': dat.append(row)\n"
            "\n"
            "fig = plt.figure()\n"
            "ax1 = plt.subplot(111)\n"
            "ax1.plot([x[0] for x in dat],[x[1] for x in dat],'g',label='Good')\n"
            "ax1.plot([x[0] for x in dat],[x[2] for x in dat],'r',label='Bad')\n"
            "ax1.set_xlabel('SOM score')\n"
            "ax1.set_ylabel('Number of training sites')\n"
            "ax1.legend(loc='best',prop={'size':8},frameon=False)\n"
            "plt.savefig('%s.eval.png')\n"
            "plt.close()\n"
            "\n", args->prefix,args->prefix
           );
    fclose(fp);
}

static void do_train(args_t *args)
{
    // read training sites
    int i, igood = 0, ibad = 0, ngood = 0, nbad = 0, ntrain = 0;
    annots_reader_reset(args);
    while ( annots_reader_next(args) )
    {
        // determine which of the nfold's SOMs to train
        int isom = 0;
        if ( args->dclass == args->good_class )
        {
            if ( ++igood >= args->nfold ) igood = 0;
            isom = igood;
            ngood++;
        }
        else if ( args->dclass == args->bad_class )
        {
            if ( ++ibad >= args->nfold ) ibad = 0;
            isom = ibad;
            nbad++;
        }
        else
            error("Could not determine the class: %d (vs %d and %d)\n", args->dclass,args->good_class,args->bad_class);

        // save the values for evaluation
        ntrain++;
        hts_expand(double, ntrain*args->mvals, args->mtrain_dat, args->train_dat);
        hts_expand(int, ntrain, args->mtrain_class, args->train_class);
        memcpy(args->train_dat+(ntrain-1)*args->mvals, args->vals, args->mvals*sizeof(double));
        args->train_class[ntrain-1] = (args->dclass==args->good_class ? 1 : 0) | isom<<1;  // store class + chunk used for training
    }
    annots_reader_close(args);

    // init maps
    if ( !args->ntrain ) args->ntrain = ngood/args->nfold;
    srandom(args->rand_seed);
    args->som = (som_t**) malloc(sizeof(som_t*)*args->nfold);
    for (i=0; i<args->nfold; i++) args->som[i] = som_init(args);

    // train
    for (i=0; i<ntrain; i++)
    {
        int is_good = args->train_class[i] & 1;
        int isom    = args->train_class[i] >> 1;
        if ( is_good || args->train_bad )
            som_train_site(args->som[isom], args->train_dat+i*args->mvals, is_good);
    }

    // norm and create plots
    for (i=0; i<args->nfold; i++)
    {
        som_norm_counts(args->som[i]);
        if ( args->prefix )
        {
            char *bname = msprintf("%s.som.%d", args->prefix,i);
            som_create_plot(args->som[i], bname);
            free(bname);
        }
    }

    // evaluate
    float *good = (float*) malloc(sizeof(float)*ngood); assert(good);
    float *bad  = (float*) malloc(sizeof(float)*nbad); assert(bad);
    igood = ibad = 0;
    double max_score = sqrt(args->som[0]->kdim);
    for (i=0; i<ntrain; i++)
    {
        double score = 0;
        int is_good = args->train_class[i] & 1;
        int isom    = args->train_class[i] >> 1;    // this vector was used for training isom-th SOM, skip
        if ( args->nfold==1 ) isom = -1;
        memcpy(args->vals, args->train_dat+i*args->mvals, args->mvals*sizeof(double));
        switch (args->merge)
        {
            case MERGE_MIN: score = get_min_score(args, isom); break;
            case MERGE_MAX: score = get_max_score(args, isom); break;
            case MERGE_AVG: score = get_avg_score(args, isom); break;
        }
        score = 1.0 - score/max_score;
        if ( is_good )
            good[igood++] = score;
        else
            bad[ibad++] = score;
    }
    qsort(good, ngood, sizeof(float), cmpfloat_desc);
    qsort(bad, nbad, sizeof(float), cmpfloat_desc);
    FILE *fp = NULL;
    if ( args->prefix ) fp = open_file(NULL,"w","%s.eval", args->prefix);
    igood = 0;
    ibad  = 0;
    float prev_score = good[0]>bad[0] ? good[0] : bad[0];
    int printed = 0;
    while ( igood<ngood || ibad<nbad )
    {
        if ( igood<ngood && good[igood]==prev_score ) { igood++; continue; }
        if ( ibad<nbad && bad[ibad]==prev_score ) { ibad++; continue; }
        if ( fp )
            fprintf(fp,"%e\t%f\t%f\n", prev_score, (float)igood/ngood, (float)ibad/nbad);
        if ( !printed && (float)igood/ngood > 0.9 )
        {
            fprintf(bcftools_stdout, "%.2f\t%.2f\t%e\t# %% of bad [1] and good [2] sites at a cutoff [3]\n", 100.*ibad/nbad,100.*igood/ngood,prev_score);
            printed = 1;
        }

        if ( igood<ngood && ibad<nbad ) prev_score = good[igood]>bad[ibad] ? good[igood] : bad[ibad];
        else if ( igood<ngood ) prev_score = good[igood];
        else prev_score = bad[ibad];
    }
    if ( !printed ) fprintf(bcftools_stdout, "%.2f\t%.2f\t%e\t# %% of bad [1] and good [2] sites at a cutoff [3]\n", 100.*ibad/nbad,100.*igood/ngood,prev_score);
    if ( fp )
    {
        if ( fclose(fp) ) error("%s.eval: fclose failed: %s\n",args->prefix,strerror(errno));
        create_eval_plot(args);
        som_write_map(args->prefix, args->som, args->nfold);
    }

    free(good);
    free(bad);
}

static void do_classify(args_t *args)
{
    annots_reader_reset(args);
    double max_score = sqrt(args->som[0]->kdim);
    while ( annots_reader_next(args) )
    {
        double score = 0;
        switch (args->merge)
        {
            case MERGE_MIN: score = get_min_score(args, -1); break;
            case MERGE_MAX: score = get_max_score(args, -1); break;
            case MERGE_AVG: score = get_avg_score(args, -1); break;
        }
        fprintf(bcftools_stdout, "%e\n", 1.0 - score/max_score);
    }
    annots_reader_close(args);
}

static void usage(void)
{
    fprintf(bcftools_stderr, "\n");
    fprintf(bcftools_stderr, "About:   SOM (Self-Organizing Map) filtering.\n");
    fprintf(bcftools_stderr, "Usage:   bcftools som --train    [options] <annots.tab.gz>\n");
    fprintf(bcftools_stderr, "         bcftools som --classify [options]\n");
    fprintf(bcftools_stderr, "\n");
    fprintf(bcftools_stderr, "Model training options:\n");
    fprintf(bcftools_stderr, "    -f, --nfold <int>                  n-fold cross-validation (number of maps) [5]\n");
    fprintf(bcftools_stderr, "    -p, --prefix <string>              prefix of output files\n");
    fprintf(bcftools_stderr, "    -s, --size <int>                   map size [20]\n");
    fprintf(bcftools_stderr, "    -t, --train                        \n");
    fprintf(bcftools_stderr, "\n");
    fprintf(bcftools_stderr, "Classifying options:\n");
    fprintf(bcftools_stderr, "    -c, --classify                     \n");
    fprintf(bcftools_stderr, "\n");
    fprintf(bcftools_stderr, "Experimental training options (no reason to change):\n");
    fprintf(bcftools_stderr, "    -b, --bmu-threshold <float>        threshold for selection of best-matching unit [0.9]\n");
    fprintf(bcftools_stderr, "    -d, --som-dimension <int>          SOM dimension [2]\n");
    fprintf(bcftools_stderr, "    -e, --exclude-bad                  exclude bad sites from training, use for evaluation only\n");
    fprintf(bcftools_stderr, "    -l, --learning-rate <float>        learning rate [1.0]\n");
    fprintf(bcftools_stderr, "    -m, --merge <min|max|avg>          -f merge algorithm [avg]\n");
    fprintf(bcftools_stderr, "    -n, --ntrain-sites <int>           effective number of training sites [number of good sites]\n");
    fprintf(bcftools_stderr, "    -r, --random-seed <int>            random seed, 0 for time() [1]\n");
    fprintf(bcftools_stderr, "\n");
    bcftools_exit(1);
}

int main_vcfsom(int argc, char *argv[])
{
    int c;
    args_t *args     = (args_t*) calloc(1,sizeof(args_t));
    args->argc       = argc; args->argv = argv;
    args->nbin       = 20;
    args->learn      = 1.0;
    args->bmu_th     = 0.9;
    args->nfold      = 5;
    args->rand_seed  = 1;
    args->ndim       = 2;
    args->bad_class  = 1;
    args->good_class = 2;
    args->merge      = MERGE_AVG;
    args->train_bad  = 1;

    static struct option loptions[] =
    {
        {"help",0,0,'h'},
        {"prefix",1,0,'p'},
        {"ntrain-sites",1,0,'n'},
        {"random-seed",1,0,'r'},
        {"bmu-threshold",1,0,'b'},
        {"exclude-bad",0,0,'e'},
        {"learning-rate",1,0,'l'},
        {"size",1,0,'s'},
        {"som-dimension",1,0,'d'},
        {"nfold",1,0,'f'},
        {"merge",1,0,'m'},
        {"train",0,0,'t'},
        {"classify",0,0,'c'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "htcp:n:r:b:l:s:f:d:m:e",loptions,NULL)) >= 0) {
        switch (c) {
            case 'e': args->train_bad = 0; break;
            case 'm':
                if ( !strcmp(optarg,"min") ) args->merge = MERGE_MIN;
                else if ( !strcmp(optarg,"max") ) args->merge = MERGE_MAX;
                else if ( !strcmp(optarg,"avg") ) args->merge = MERGE_AVG;
                else error("The -m method not recognised: %s\n", optarg);
                break;
            case 'p': args->prefix = optarg; break;
            case 'n': args->ntrain = atoi(optarg); break;
            case 'r': args->rand_seed = atoi(optarg); break;
            case 'b': args->bmu_th = atof(optarg); break;
            case 'l': args->learn = atof(optarg); break;
            case 's': args->nbin = atoi(optarg); break;
            case 'f': args->nfold = atoi(optarg); break;
            case 'd':
                args->ndim = atoi(optarg);
                if ( args->ndim<2 ) error("Expected -d >=2, got %d\n", args->ndim);
                if ( args->ndim>3 ) fprintf(bcftools_stderr,"Warning: This will take a long time and is not going to make the results better: -d %d\n", args->ndim);
                break;
            case 't': args->action = SOM_TRAIN; break;
            case 'c': args->action = SOM_CLASSIFY; break;
            case 'h':
            case '?': usage(); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( !args->rand_seed ) args->rand_seed = time(NULL);
    if ( argc!=optind+1 ) usage();
    args->fname = argv[optind];
    init_data(args);

    if ( args->action == SOM_TRAIN ) do_train(args);
    else if ( args->action == SOM_CLASSIFY ) do_classify(args);

    destroy_data(args);
    free(args);
    return 0;
}

