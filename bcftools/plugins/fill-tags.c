/* The MIT License

   Copyright (c) 2015 Genome Research Ltd.

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

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <getopt.h>
#include <math.h>
#include <htslib/hts.h>
#include <htslib/kseq.h>
#include <htslib/vcf.h>
#include <htslib/khash_str2int.h>
#include "bcftools.h"

#define SET_AN      (1<<0)
#define SET_AC      (1<<1)
#define SET_AC_Hom  (1<<2)
#define SET_AC_Het  (1<<3)
#define SET_AC_Hemi (1<<4)
#define SET_AF      (1<<5)
#define SET_NS      (1<<6)
#define SET_MAF     (1<<7)
#define SET_HWE     (1<<8)

typedef struct
{
    int nhom, nhet, nhemi, nac;
}
counts_t;

typedef struct
{
    int ns;
    int ncounts, mcounts;
    counts_t *counts;
    char *name, *suffix;
    int nsmpl, *smpl;
}
pop_t;

typedef struct
{
    bcf_hdr_t *in_hdr, *out_hdr;
    int npop, tags, drop_missing, gt_id;
    pop_t *pop, **smpl2pop;
    float *farr;
    int32_t *iarr, niarr, miarr, nfarr, mfarr;
    double *hwe_probs;
    int mhwe_probs;
    kstring_t str;
}
args_t;

static args_t *args;

const char *about(void)
{
    return "Set INFO tags AF, AC, AC_Hemi, AC_Hom, AC_Het, AN, HWE, MAF, NS.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Set INFO tags AF, AC, AC_Hemi, AC_Hom, AC_Het, AN, HWE, MAF, NS.\n"
        "Usage: bcftools +fill-tags [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -d, --drop-missing          do not count half-missing genotypes \"./1\" as hemizygous\n"
        "   -t, --tags LIST             list of output tags. By default, all tags are filled.\n"
        "   -S, --samples-file FILE     list of samples (first column) and comma-separated list of populations (second column)\n"
        "\n"
        "Example:\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -t AN,AC\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -d\n"
        "   bcftools +fill-tags in.bcf -Ob -o out.bcf -- -S sample-group.txt -t HWE\n"
        "\n";
}

void parse_samples(args_t *args, char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);

    void *pop2i = khash_str2int_init();
    void *smpli = khash_str2int_init();
    kstring_t str = {0,0,0};

    int moff = 0, *off = NULL, nsmpl = 0;
    while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 )
    {
        // NA12400 GRP1
        // NA18507 GRP1,GRP2
        char *pop_names = str.s + str.l - 1;
        while ( pop_names >= str.s && isspace(*pop_names) ) pop_names--;
        if ( pop_names <= str.s ) error("Could not parse the file: %s\n", str.s);
        pop_names[1] = 0;   // trailing spaces
        while ( pop_names >= str.s && !isspace(*pop_names) ) pop_names--;
        if ( pop_names <= str.s ) error("Could not parse the file: %s\n", str.s);

        char *smpl = pop_names++;
        while ( smpl >= str.s && isspace(*smpl) ) smpl--;
        if ( smpl <= str.s+1 ) error("Could not parse the file: %s\n", str.s);
        smpl[1] = 0;
        smpl = str.s;

        int ismpl = bcf_hdr_id2int(args->in_hdr,BCF_DT_SAMPLE,smpl);
        if ( ismpl<0 ) 
        {
            fprintf(stderr,"Warning: The sample not present in the VCF: %s\n",smpl);
            continue;
        }
        if ( khash_str2int_has_key(smpli,smpl) )
        {
            fprintf(stderr,"Warning: The sample is listed twice in %s: %s\n",fname,smpl);
            continue;
        }
        khash_str2int_inc(smpli,strdup(smpl));

        int i,npops = ksplit_core(pop_names,',',&moff,&off);
        for (i=0; i<npops; i++)
        {
            char *pop_name = &pop_names[off[i]];
            if ( !khash_str2int_has_key(pop2i,pop_name) )
            {
                pop_name = strdup(pop_name);
                khash_str2int_set(pop2i,pop_name,args->npop);
                args->npop++;
                args->pop = (pop_t*) realloc(args->pop,args->npop*sizeof(*args->pop));
                memset(args->pop+args->npop-1,0,sizeof(*args->pop));
                args->pop[args->npop-1].name = pop_name;
                args->pop[args->npop-1].suffix = (char*)malloc(strlen(pop_name)+2);
                memcpy(args->pop[args->npop-1].suffix+1,pop_name,strlen(pop_name)+1);
                args->pop[args->npop-1].suffix[0] = '_';
            }
            int ipop = 0;
            khash_str2int_get(pop2i,pop_name,&ipop);
            pop_t *pop = &args->pop[ipop];
            pop->nsmpl++;
            pop->smpl = (int*) realloc(pop->smpl,pop->nsmpl*sizeof(*pop->smpl));
            pop->smpl[pop->nsmpl-1] = ismpl;
        }
        nsmpl++;
    }

    if ( nsmpl != bcf_hdr_nsamples(args->in_hdr) )
        fprintf(stderr,"Warning: %d samples in the list, %d samples in the VCF.\n", nsmpl,bcf_hdr_nsamples(args->in_hdr));

    if ( !args->npop ) error("No populations given?\n");

    khash_str2int_destroy(pop2i);
    khash_str2int_destroy_free(smpli);
    free(str.s);
    free(off);
    hts_close(fp);
}

void init_pops(args_t *args)
{
    int i,j, nsmpl;

    // add the population "ALL", which is a summary population for all samples
    args->npop++;
    args->pop = (pop_t*) realloc(args->pop,args->npop*sizeof(*args->pop));
    memset(args->pop+args->npop-1,0,sizeof(*args->pop));
    args->pop[args->npop-1].name   = strdup("");
    args->pop[args->npop-1].suffix = strdup("");

    nsmpl = bcf_hdr_nsamples(args->in_hdr);
    args->smpl2pop = (pop_t**) calloc(nsmpl*(args->npop+1),sizeof(pop_t*));
    for (i=0; i<nsmpl; i++)
        args->smpl2pop[i*(args->npop+1)] = &args->pop[args->npop-1];

    for (i=0; i<args->npop; i++)
    {
        for (j=0; j<args->pop[i].nsmpl; j++)
        {
            int ismpl = args->pop[i].smpl[j];
            pop_t **smpl2pop = &args->smpl2pop[ismpl*(args->npop+1)];
            while (*smpl2pop) smpl2pop++;
            *smpl2pop = &args->pop[i];
        }
    }
}

int parse_tags(args_t *args, const char *str)
{
    int i, flag = 0, n_tags;
    char **tags = hts_readlist(str, 0, &n_tags);
    for(i=0; i<n_tags; i++)
    {
        if ( !strcasecmp(tags[i],"AN") ) flag |= SET_AN;
        else if ( !strcasecmp(tags[i],"AC") ) flag |= SET_AC;
        else if ( !strcasecmp(tags[i],"NS") ) flag |= SET_NS;
        else if ( !strcasecmp(tags[i],"AC_Hom") ) flag |= SET_AC_Hom;
        else if ( !strcasecmp(tags[i],"AC_Het") ) flag |= SET_AC_Het;
        else if ( !strcasecmp(tags[i],"AC_Hemi") ) flag |= SET_AC_Hemi;
        else if ( !strcasecmp(tags[i],"AF") ) flag |= SET_AF;
        else if ( !strcasecmp(tags[i],"MAF") ) flag |= SET_MAF;
        else if ( !strcasecmp(tags[i],"HWE") ) flag |= SET_HWE;
        else
        {
            fprintf(stderr,"Error parsing \"--tags %s\": the tag \"%s\" is not supported\n", str,tags[i]);
            exit(1);
        }
        free(tags[i]);
    }
    if (n_tags) free(tags);
    return flag;
}

void hdr_append(args_t *args, char *fmt)
{
    int i;
    for (i=0; i<args->npop; i++)
        bcf_hdr_printf(args->out_hdr, fmt, args->pop[i].suffix,*args->pop[i].name ? " in " : "",args->pop[i].name);
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    args = (args_t*) calloc(1,sizeof(args_t));
    args->in_hdr  = in;
    args->out_hdr = out;
    char *samples_fname = NULL;
    static struct option loptions[] =
    {
        {"drop-missing",0,0,'d'},
        {"tags",1,0,'t'},
        {"samples-file",1,0,'S'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "?ht:dS:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'd': args->drop_missing = 1; break;
            case 't': args->tags |= parse_tags(args,optarg); break;
            case 'S': samples_fname = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }

    if ( optind != argc ) error(usage());

    args->gt_id = bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,"GT");
    if ( args->gt_id<0 ) error("Error: GT field is not present\n");

    if ( !args->tags )
        for (c=0; c<=8; c++) args->tags |= 1<<c;    // by default all tags will be filled

    if ( samples_fname ) parse_samples(args, samples_fname);
    init_pops(args);

    if ( args->tags & SET_AN ) hdr_append(args, "##INFO=<ID=AN%s,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes%s%s\">");
    if ( args->tags & SET_AC ) hdr_append(args, "##INFO=<ID=AC%s,Number=A,Type=Integer,Description=\"Allele count in genotypes%s%s\">");
    if ( args->tags & SET_NS ) hdr_append(args, "##INFO=<ID=NS%s,Number=1,Type=Integer,Description=\"Number of samples with data%s%s\">");
    if ( args->tags & SET_AC_Hom ) hdr_append(args, "##INFO=<ID=AC_Hom%s,Number=A,Type=Integer,Description=\"Allele counts in homozygous genotypes%s%s\">");
    if ( args->tags & SET_AC_Het ) hdr_append(args, "##INFO=<ID=AC_Het%s,Number=A,Type=Integer,Description=\"Allele counts in heterozygous genotypes%s%s\">");
    if ( args->tags & SET_AC_Hemi ) hdr_append(args, "##INFO=<ID=AC_Hemi%s,Number=A,Type=Integer,Description=\"Allele counts in hemizygous genotypes%s%s\">");
    if ( args->tags & SET_AF ) hdr_append(args, "##INFO=<ID=AF%s,Number=A,Type=Float,Description=\"Allele frequency%s%s\">");
    if ( args->tags & SET_MAF ) hdr_append(args, "##INFO=<ID=MAF%s,Number=A,Type=Float,Description=\"Minor Allele frequency%s%s\">");
    if ( args->tags & SET_HWE ) hdr_append(args, "##INFO=<ID=HWE%s,Number=A,Type=Float,Description=\"HWE test%s%s (PMID:15789306)\">");

    return 0;
}

/* 
    Wigginton 2005, PMID: 15789306 

    nref .. number of reference alleles
    nalt .. number of alt alleles
    nhet .. number of het genotypes, assuming number of genotypes = (nref+nalt)*2

*/
float calc_hwe(args_t *args, int nref, int nalt, int nhet)
{
    int ngt   = (nref+nalt) / 2;
    int nrare = nref < nalt ? nref : nalt;

    // sanity check: there is odd/even number of rare alleles iff there is odd/even number of hets
    if ( (nrare & 1) ^ (nhet & 1) ) error("nrare/nhet should be both odd or even: nrare=%d nref=%d nalt=%d nhet=%d\n",nrare,nref,nalt,nhet);
    if ( nrare < nhet ) error("Fewer rare alleles than hets? nrare=%d nref=%d nalt=%d nhet=%d\n",nrare,nref,nalt,nhet);
    if ( (nref+nalt) & 1 ) error("Expected diploid genotypes: nref=%d nalt=%d\n",nref,nalt);

    // initialize het probs
    hts_expand(double,nrare+1,args->mhwe_probs,args->hwe_probs);
    memset(args->hwe_probs, 0, sizeof(*args->hwe_probs)*(nrare+1));
    double *probs = args->hwe_probs;

    // start at midpoint
    int mid = nrare * (nref + nalt - nrare) / (nref + nalt);

    // check to ensure that midpoint and rare alleles have same parity
    if ( (nrare & 1) ^ (mid & 1) ) mid++;

    int het = mid;
    int hom_r  = (nrare - mid) / 2;
    int hom_c  = ngt - het - hom_r;
    double sum = probs[mid] = 1.0;

    for (het = mid; het > 1; het -= 2)
    {
        probs[het - 2] = probs[het] * het * (het - 1.0) / (4.0 * (hom_r + 1.0) * (hom_c + 1.0));
        sum += probs[het - 2];

        // 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
        hom_r++;
        hom_c++;
    }

    het = mid;
    hom_r = (nrare - mid) / 2;
    hom_c = ngt - het - hom_r;
    for (het = mid; het <= nrare - 2; het += 2)
    {
        probs[het + 2] = probs[het] * 4.0 * hom_r * hom_c / ((het + 2.0) * (het + 1.0));
        sum += probs[het + 2];

        // add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
        hom_r--;
        hom_c--;
    }

    for (het=0; het<nrare+1; het++) probs[het] /= sum;

    double p_rank = 0.0;
    for (het=0; het <= nrare; het++)
    {
        if ( probs[het] > probs[nhet]) continue;
        p_rank += probs[het];
    }

    return p_rank > 1 ? 1.0 : p_rank;
}

static inline void set_counts(pop_t *pop, int is_half, int is_hom, int is_hemi, int als)
{
    int ial;
    for (ial=0; als; ial++)
    {
        if ( als&1 )
        { 
            if ( is_half ) pop->counts[ial].nac++;
            else if ( !is_hom ) pop->counts[ial].nhet++;
            else if ( !is_hemi ) pop->counts[ial].nhom += 2;
            else pop->counts[ial].nhemi++;
        }
        als >>= 1;
    }
    pop->ns++;
}
static void clean_counts(pop_t *pop, int nals)
{
    pop->ns = 0;
    memset(pop->counts,0,sizeof(counts_t)*nals);
}

bcf1_t *process(bcf1_t *rec)
{
    int i,j, nsmpl = bcf_hdr_nsamples(args->in_hdr);;

    bcf_unpack(rec, BCF_UN_FMT);
    bcf_fmt_t *fmt_gt = NULL;
    for (i=0; i<rec->n_fmt; i++)
        if ( rec->d.fmt[i].id==args->gt_id ) { fmt_gt = &rec->d.fmt[i]; break; }
    if ( !fmt_gt ) return rec;    // no GT tag

    hts_expand(int32_t,rec->n_allele, args->miarr, args->iarr);
    hts_expand(float,rec->n_allele, args->mfarr, args->farr);
    for (i=0; i<args->npop; i++)
        hts_expand(counts_t,rec->n_allele,args->pop[i].mcounts, args->pop[i].counts);

    for (i=0; i<args->npop; i++)
        clean_counts(&args->pop[i], rec->n_allele);

    assert( rec->n_allele < 8*sizeof(int) );

    #define BRANCH_INT(type_t,vector_end) \
    { \
        for (i=0; i<nsmpl; i++) \
        { \
            type_t *p = (type_t*) (fmt_gt->p + i*fmt_gt->size); \
            int ial, als = 0, nals = 0, is_half, is_hom, is_hemi; \
            for (ial=0; ial<fmt_gt->n; ial++) \
            { \
                if ( p[ial]==vector_end ) break; /* smaller ploidy */ \
                if ( bcf_gt_is_missing(p[ial]) ) continue; /* missing allele */ \
                int idx = bcf_gt_allele(p[ial]); \
                nals++; \
                \
                if ( idx >= rec->n_allele ) \
                    error("Incorrect allele (\"%d\") in %s at %s:%d\n",idx,args->in_hdr->samples[i],bcf_seqname(args->in_hdr,rec),rec->pos+1); \
                als |= (1<<idx);  /* this breaks with too many alleles */ \
            } \
            if ( nals==0 ) continue; /* missing genotype */ \
            is_hom = als && !(als & (als-1)); /* only one bit is set */ \
            if ( nals!=ial ) \
            { \
                if ( args->drop_missing ) is_hemi = 0, is_half = 1; \
                else is_hemi = 1, is_half = 0; \
            } \
            else if ( nals==1 ) is_hemi = 1, is_half = 0; \
            else is_hemi = 0, is_half = 0; \
            pop_t **pop = &args->smpl2pop[i*(args->npop+1)]; \
            while ( *pop ) { set_counts(*pop,is_half,is_hom,is_hemi,als); pop++; }\
        } \
    }
    switch (fmt_gt->type) {
        case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_vector_end); break;
        default: error("The GT type is not recognised: %d at %s:%d\n",fmt_gt->type, bcf_seqname(args->in_hdr,rec),rec->pos+1); break;
    }
    #undef BRANCH_INT

    if ( args->tags & SET_NS )
    {
        for (i=0; i<args->npop; i++)
        {
            args->str.l = 0;
            ksprintf(&args->str, "NS%s", args->pop[i].suffix);
            if ( bcf_update_info_int32(args->out_hdr,rec,args->str.s,&args->pop[i].ns,1)!=0 )
                error("Error occurred while updating %s at %s:%d\n", args->str.s,bcf_seqname(args->in_hdr,rec),rec->pos+1);
        }
    }
    if ( args->tags & SET_AN )
    {
        for (i=0; i<args->npop; i++)
        {
            pop_t *pop = &args->pop[i];
            int32_t an = 0;
            for (j=0; j<rec->n_allele; j++) 
                an += pop->counts[j].nhet + pop->counts[j].nhom + pop->counts[j].nhemi + pop->counts[j].nac;

            args->str.l = 0;
            ksprintf(&args->str, "AN%s", args->pop[i].suffix);
            if ( bcf_update_info_int32(args->out_hdr,rec,args->str.s,&an,1)!=0 )
                error("Error occurred while updating %s at %s:%d\n", args->str.s,bcf_seqname(args->in_hdr,rec),rec->pos+1);
        }
    }
    if ( args->tags & (SET_AF | SET_MAF) )
    {
        for (i=0; i<args->npop; i++)
        {
            int32_t an = 0;
            if ( rec->n_allele > 1 )
            {
                pop_t *pop = &args->pop[i];
                memset(args->farr, 0, sizeof(*args->farr)*(rec->n_allele-1));
                for (j=1; j<rec->n_allele; j++) 
                    args->farr[j-1] += pop->counts[j].nhet + pop->counts[j].nhom + pop->counts[j].nhemi + pop->counts[j].nac;
                an = pop->counts[0].nhet + pop->counts[0].nhom + pop->counts[0].nhemi + pop->counts[0].nac;
                for (j=1; j<rec->n_allele; j++) an += args->farr[j-1];
                if ( !an ) continue;
                for (j=1; j<rec->n_allele; j++) args->farr[j-1] /= an;
            }
            if ( args->tags & SET_AF )
            {
                args->str.l = 0;
                ksprintf(&args->str, "AF%s", args->pop[i].suffix);
                if ( bcf_update_info_float(args->out_hdr,rec,args->str.s,args->farr,rec->n_allele-1)!=0 )
                    error("Error occurred while updating %s at %s:%d\n", args->str.s,bcf_seqname(args->in_hdr,rec),rec->pos+1);
            }
            if ( args->tags & SET_MAF )
            {
                if ( !an ) continue;
                for (j=1; j<rec->n_allele; j++)
                    if ( args->farr[j-1] > 0.5 ) args->farr[j-1] = 1 - args->farr[j-1];     // todo: this is incorrect for multiallelic sites
                args->str.l = 0;
                ksprintf(&args->str, "MAF%s", args->pop[i].suffix);
                if ( bcf_update_info_float(args->out_hdr,rec,args->str.s,args->farr,rec->n_allele-1)!=0 )
                    error("Error occurred while updating %s at %s:%d\n", args->str.s,bcf_seqname(args->in_hdr,rec),rec->pos+1);
            }
        }
    }
    if ( args->tags & SET_AC )
    {
        for (i=0; i<args->npop; i++)
        {
            if ( rec->n_allele > 1 )
            {
                pop_t *pop = &args->pop[i];
                memset(args->iarr, 0, sizeof(*args->iarr)*(rec->n_allele-1));
                for (j=1; j<rec->n_allele; j++) 
                    args->iarr[j-1] += pop->counts[j].nhet + pop->counts[j].nhom + pop->counts[j].nhemi + pop->counts[j].nac;
            }
            args->str.l = 0;
            ksprintf(&args->str, "AC%s", args->pop[i].suffix);
            if ( bcf_update_info_int32(args->out_hdr,rec,args->str.s,args->iarr,rec->n_allele-1)!=0 )
                error("Error occurred while updating %s at %s:%d\n", args->str.s,bcf_seqname(args->in_hdr,rec),rec->pos+1);
        }
    }
    if ( args->tags & SET_AC_Het )
    {
        for (i=0; i<args->npop; i++)
        {
            if ( rec->n_allele > 1 )
            {
                pop_t *pop = &args->pop[i];
                memset(args->iarr, 0, sizeof(*args->iarr)*(rec->n_allele-1));
                for (j=1; j<rec->n_allele; j++) 
                    args->iarr[j-1] += pop->counts[j].nhet;
            }
            args->str.l = 0;
            ksprintf(&args->str, "AC_Het%s", args->pop[i].suffix);
            if ( bcf_update_info_int32(args->out_hdr,rec,args->str.s,args->iarr,rec->n_allele-1)!=0 )
                error("Error occurred while updating %s at %s:%d\n", args->str.s,bcf_seqname(args->in_hdr,rec),rec->pos+1);
        }
    }
    if ( args->tags & SET_AC_Hom )
    {
        for (i=0; i<args->npop; i++)
        {
            if ( rec->n_allele > 1 )
            {
                pop_t *pop = &args->pop[i];
                memset(args->iarr, 0, sizeof(*args->iarr)*(rec->n_allele-1));
                for (j=1; j<rec->n_allele; j++) 
                    args->iarr[j-1] += pop->counts[j].nhom;
            }
            args->str.l = 0;
            ksprintf(&args->str, "AC_Hom%s", args->pop[i].suffix);
            if ( bcf_update_info_int32(args->out_hdr,rec,args->str.s,args->iarr,rec->n_allele-1)!=0 )
                error("Error occurred while updating %s at %s:%d\n", args->str.s,bcf_seqname(args->in_hdr,rec),rec->pos+1);
        }
    }
    if ( args->tags & SET_AC_Hemi && rec->n_allele > 1 )
    {
        for (i=0; i<args->npop; i++)
        {
            if ( rec->n_allele > 1 )
            {
                pop_t *pop = &args->pop[i];
                memset(args->iarr, 0, sizeof(*args->iarr)*(rec->n_allele-1));
                for (j=1; j<rec->n_allele; j++) 
                    args->iarr[j-1] += pop->counts[j].nhemi;
            }
            args->str.l = 0;
            ksprintf(&args->str, "AC_Hemi%s", args->pop[i].suffix);
            if ( bcf_update_info_int32(args->out_hdr,rec,args->str.s,args->iarr,rec->n_allele-1)!=0 )
                error("Error occurred while updating %s at %s:%d\n", args->str.s,bcf_seqname(args->in_hdr,rec),rec->pos+1);
        }
    }
    if ( args->tags & SET_HWE )
    {
        for (i=0; i<args->npop; i++)
        {
            if ( rec->n_allele > 1 )
            {
                pop_t *pop = &args->pop[i];
                memset(args->farr, 0, sizeof(*args->farr)*(rec->n_allele-1));
                int nref_tot = pop->counts[0].nhom;
                for (j=0; j<rec->n_allele; j++) nref_tot += pop->counts[j].nhet;   // NB this neglects multiallelic genotypes
                for (j=1; j<rec->n_allele; j++) 
                {
                    int nref = nref_tot - pop->counts[j].nhet;
                    int nalt = pop->counts[j].nhet + pop->counts[j].nhom;
                    int nhet = pop->counts[j].nhet;
                    args->farr[j-1] = (nref>0 && nalt>0) ? calc_hwe(args, nref, nalt, nhet) : 1;
                }
            }
            args->str.l = 0;
            ksprintf(&args->str, "HWE%s", args->pop[i].suffix);
            if ( bcf_update_info_float(args->out_hdr,rec,args->str.s,args->farr,rec->n_allele-1)!=0 )
                error("Error occurred while updating %s at %s:%d\n", args->str.s,bcf_seqname(args->in_hdr,rec),rec->pos+1);
        }
    }

    return rec;
}

void destroy(void)
{
    int i; 
    for (i=0; i<args->npop; i++)
    {
        free(args->pop[i].name);
        free(args->pop[i].suffix);
        free(args->pop[i].smpl);
        free(args->pop[i].counts);
    }
    free(args->str.s);
    free(args->pop);
    free(args->smpl2pop);
    free(args->iarr);
    free(args->farr);
    free(args->hwe_probs);
    free(args);
}



