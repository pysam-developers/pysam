/*  vcfmerge.c -- Merge multiple VCF/BCF files to create one multi-sample file.

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

#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <unistd.h>
#include <getopt.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <htslib/faidx.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include "bcftools.h"
#include "regidx.h"
#include "vcmp.h"

#define DBG 0

#include <htslib/khash.h>
KHASH_MAP_INIT_STR(strdict, int)
typedef khash_t(strdict) strdict_t;

#define FLT_LOGIC_ADD    0
#define FLT_LOGIC_REMOVE 1

#define SKIP_DONE 1     // the record was processed
#define SKIP_DIFF 2     // not compatible, merge later

#define IS_VL_G(hdr,id) (bcf_hdr_id2length(hdr,BCF_HL_FMT,id) == BCF_VL_G)
#define IS_VL_A(hdr,id) (bcf_hdr_id2length(hdr,BCF_HL_FMT,id) == BCF_VL_A)
#define IS_VL_R(hdr,id) (bcf_hdr_id2length(hdr,BCF_HL_FMT,id) == BCF_VL_R)

#define SWAP(type_t,a,b) { type_t tmp = (a); (a) = (b); (b) = tmp; }

// For merging INFO Number=A,G,R tags
typedef struct
{
    const char *hdr_tag;
    int type, nvals;
    int nbuf, mbuf;
    uint8_t *buf;
}
AGR_info_t;

// Rules for merging arbitrary INFO tags
typedef struct _info_rule_t
{
    char *hdr_tag;
    void (*merger)(bcf_hdr_t *hdr, bcf1_t *line, struct _info_rule_t *rule);
    int type;           // one of BCF_HT_*
    int block_size;     // number of values in a block
    int type_size;      // size of the corresponding BCF_HT_* type
    int nblocks;        // number of blocks in nvals (the number of merged files)
    int nvals, mvals;   // used and total size of vals array
    void *vals;         // the info tag values
}
info_rule_t;

typedef struct
{
    bcf1_t *line;
    int end, active;
}
gvcf_aux_t;

// Auxiliary merge data for selecting the right combination
//  of buffered records across multiple readers. maux1_t
//  corresponds to one buffered line.
typedef struct
{
    int skip;
    int *map;   // mapping from input alleles to the array of output alleles (set by merge_alleles)
    int mmap;   // size of map array (only buffer[i].n_allele is actually used)
    int als_differ;
}
maux1_t;
typedef struct
{
    int rid;        // current rid
    int beg,end;    // valid ranges in reader's buffer [beg,end). Maintained by maux_reset and gvcf_flush.
    int cur;        // current line or -1 if none
    int mrec;       // allocated size of buf
    maux1_t *rec;   // buffer to keep reader's lines
    bcf1_t **lines; // source buffer: either gvcf or readers' buffer
}
buffer_t;
typedef struct
{
    int n, pos, var_types;  // number of readers, current position, currently available variant types
    char *chr;              // current chromosome
    char **als, **out_als;  // merged alleles (temp, may contain empty records) and merged alleles ready for output
    int nals, mals, nout_als, mout_als; // size of the output array
    int *cnt, ncnt; // number of records that refer to the alleles
    int *smpl_ploidy, *smpl_nGsize; // ploidy and derived number of values in Number=G tags, updated for each line (todo: cache for missing cases)
    bcf_fmt_t **fmt_map; // i-th output FORMAT field corresponds in j-th reader to i*nreader+j, first row is reserved for GT
    int nfmt_map;        // number of rows in the fmt_map array
    int *agr_map, nagr_map, magr_map;   // mapping between Number=AGR element indexes
    void *tmp_arr;
    int ntmp_arr;
    buffer_t *buf;
    AGR_info_t *AGR_info;
    int nAGR_info, mAGR_info;
    bcf_srs_t *files;
    int gvcf_min, gvcf_break;   // min buffered gvcf END position (NB: gvcf_min is 1-based) or 0 if no active lines are present
    gvcf_aux_t *gvcf;           // buffer of gVCF lines
}
maux_t;

typedef struct
{
    vcmp_t *vcmp;
    maux_t *maux;
    regidx_t *regs;    // apply regions only after the blocks are expanded
    regitr_t *regs_itr;
    int header_only, collapse, output_type, force_samples, merge_by_id, do_gvcf, filter_logic, missing_to_ref;
    char *header_fname, *output_fname, *regions_list, *info_rules, *file_list;
    faidx_t *gvcf_fai;
    info_rule_t *rules;
    int nrules;
    strdict_t *tmph;
    kstring_t tmps;
    bcf_srs_t *files;
    bcf1_t *out_line;
    htsFile *out_fh;
    bcf_hdr_t *out_hdr;
    char **argv;
    int argc, n_threads, record_cmd_line;
}
args_t;

static bcf1_t *maux_get_line(args_t *args, int i)
{
    maux_t *ma = args->maux;
    int ibuf = ma->buf[i].cur;
    if ( ibuf >= 0 ) return ma->buf[i].lines[ibuf];
    return NULL;
}

static void info_rules_merge_sum(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    int i, j, ndim = rule->block_size;
    #define BRANCH(type_t,is_missing) { \
        type_t *ptr = (type_t*) rule->vals; \
        for (i=0; i<rule->nvals; i++) if ( is_missing ) ptr[i] = 0; \
        for (i=1; i<rule->nblocks; i++) \
        { \
            for (j=0; j<ndim; j++) ptr[j] += ptr[j+i*ndim]; \
        } \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int32_t, ptr[i]==bcf_int32_missing); break;
        case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(ptr[i])); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,ndim,rule->type);
}
static void info_rules_merge_avg(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    int i, j, ndim = rule->block_size;
    #define BRANCH(type_t,is_missing) { \
        type_t *ptr = (type_t*) rule->vals; \
        for (i=0; i<rule->nvals; i++) if ( is_missing ) ptr[i] = 0; \
        for (j=0; j<ndim; j++) \
        { \
            double sum = 0; \
            for (i=0; i<rule->nblocks; i++) sum += ptr[j+i*ndim]; \
            ptr[j] = sum / rule->nblocks; \
        } \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int32_t, ptr[i]==bcf_int32_missing); break;
        case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(ptr[i])); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,ndim,rule->type);
}
static void info_rules_merge_min(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    int i, j, ndim = rule->block_size;
    #define BRANCH(type_t,is_missing,set_missing,huge_val) { \
        type_t *ptr = (type_t*) rule->vals; \
        for (i=0; i<rule->nvals; i++) if ( is_missing ) ptr[i] = huge_val; \
        for (i=1; i<rule->nblocks; i++) \
        { \
            for (j=0; j<ndim; j++) if ( ptr[j] > ptr[j+i*ndim] ) ptr[j] = ptr[j+i*ndim]; \
        } \
        for (i=0; i<rule->nvals; i++) if ( ptr[i]==huge_val ) set_missing; \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int32_t, ptr[i]==bcf_int32_missing, ptr[i]=bcf_int32_missing, INT32_MAX); break;
        case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(ptr[i]), bcf_float_set_missing(ptr[i]), HUGE_VAL); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,ndim,rule->type);
}
static void info_rules_merge_max(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    int i, j, ndim = rule->block_size;
    #define BRANCH(type_t,is_missing,set_missing,huge_val) { \
        type_t *ptr = (type_t*) rule->vals; \
        for (i=0; i<rule->nvals; i++) if ( is_missing ) ptr[i] = huge_val; \
        for (i=1; i<rule->nblocks; i++) \
        { \
            for (j=0; j<ndim; j++) if ( ptr[j] < ptr[j+i*ndim] ) ptr[j] = ptr[j+i*ndim]; \
        } \
        for (i=0; i<rule->nvals; i++) if ( ptr[i]==huge_val ) set_missing; \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int32_t, ptr[i]==bcf_int32_missing, ptr[i]=bcf_int32_missing, INT32_MIN); break;
        case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(ptr[i]), bcf_float_set_missing(ptr[i]), -HUGE_VAL); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,ndim,rule->type);
}
static void info_rules_merge_join(bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule)
{
    if ( !rule->nvals ) return;
    if ( rule->type==BCF_HT_STR )
    {
        ((char*)rule->vals)[rule->nvals] = 0;
        bcf_update_info_string(hdr,line,rule->hdr_tag,rule->vals);
    }
    else
        bcf_update_info(hdr,line,rule->hdr_tag,rule->vals,rule->nvals,rule->type);
}

static int info_rules_comp_key2(const void *a, const void *b)
{
    info_rule_t *rule1 = (info_rule_t*) a;
    info_rule_t *rule2 = (info_rule_t*) b;
    return strcmp(rule1->hdr_tag, rule2->hdr_tag);
}
static int info_rules_comp_key(const void *a, const void *b)
{
    char *key = (char*) a;
    info_rule_t *rule = (info_rule_t*) b;
    return strcmp(key, rule->hdr_tag);
}
static void info_rules_init(args_t *args)
{
    if ( args->info_rules && !strcmp("-",args->info_rules) ) return;

    kstring_t str = {0,0,0};
    if ( !args->info_rules )
    {
        if ( bcf_hdr_idinfo_exists(args->out_hdr,BCF_HL_INFO,bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, "DP")) ) kputs("DP:sum",&str);
        if ( bcf_hdr_idinfo_exists(args->out_hdr,BCF_HL_INFO,bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, "DP4")) )
        {
            if ( str.l ) kputc(',',&str);
            kputs("DP4:sum",&str);
        }
        if ( args->do_gvcf && bcf_hdr_idinfo_exists(args->out_hdr,BCF_HL_INFO,bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, "QS")) )
        {
            if ( str.l ) kputc(',',&str);
            kputs("QS:sum",&str);
        }
        if ( args->do_gvcf && bcf_hdr_idinfo_exists(args->out_hdr,BCF_HL_INFO,bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, "MinDP")) )
        {
            if ( str.l ) kputc(',',&str);
            kputs("MinDP:min",&str);
        }
        if ( args->do_gvcf && bcf_hdr_idinfo_exists(args->out_hdr,BCF_HL_INFO,bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, "I16")) )
        {
            if ( str.l ) kputc(',',&str);
            kputs("I16:sum",&str);
        }
        if ( args->do_gvcf && bcf_hdr_idinfo_exists(args->out_hdr,BCF_HL_INFO,bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, "IDV")) )
        {
            if ( str.l ) kputc(',',&str);
            kputs("IDV:max",&str);
        }
        if ( args->do_gvcf && bcf_hdr_idinfo_exists(args->out_hdr,BCF_HL_INFO,bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, "IMF")) )
        {
            if ( str.l ) kputc(',',&str);
            kputs("IMF:max",&str);
        }

        if ( !str.l ) return;
        args->info_rules = str.s;
    }

    args->nrules = 1;
    char *ss = strdup(args->info_rules), *tmp = ss;
    int n = 0;
    while ( *ss )
    {
        if ( *ss==':' ) { *ss = 0; n++; if ( n%2==0 ) error("Could not parse INFO rules: \"%s\"\n", args->info_rules); }
        else if ( *ss==',' ) { *ss = 0; args->nrules++; n++; if ( n%2==1 ) error("Could not parse INFO rules: \"%s\"\n", args->info_rules); }
        ss++;
    }
    if ( n%2==0 ) error("Could not parse INFO rules: \"%s\"\n", args->info_rules);
    args->rules = (info_rule_t*) calloc(args->nrules,sizeof(info_rule_t));

    n = 0;
    ss = tmp;
    while ( n < args->nrules )
    {
        info_rule_t *rule = &args->rules[n];
        rule->hdr_tag = strdup(ss);
        int id = bcf_hdr_id2int(args->out_hdr, BCF_DT_ID, rule->hdr_tag);
        if ( !bcf_hdr_idinfo_exists(args->out_hdr,BCF_HL_INFO,id) ) error("The tag is not defined in the header: \"%s\"\n", rule->hdr_tag);
        rule->type = bcf_hdr_id2type(args->out_hdr,BCF_HL_INFO,id);
        if ( rule->type==BCF_HT_INT ) rule->type_size = sizeof(int32_t);
        else if ( rule->type==BCF_HT_REAL ) rule->type_size = sizeof(float);
        else if ( rule->type==BCF_HT_STR ) rule->type_size = sizeof(char); 
        else error("The type is not supported: \"%s\"\n", rule->hdr_tag);

        ss = strchr(ss, '\0'); ss++;
        if ( !*ss ) error("Could not parse INFO rules, missing logic of \"%s\"\n", rule->hdr_tag);

        int is_join = 0;
        if ( !strcasecmp(ss,"sum") ) rule->merger = info_rules_merge_sum;
        else if ( !strcasecmp(ss,"avg") ) rule->merger = info_rules_merge_avg;
        else if ( !strcasecmp(ss,"min") ) rule->merger = info_rules_merge_min;
        else if ( !strcasecmp(ss,"max") ) rule->merger = info_rules_merge_max;
        else if ( !strcasecmp(ss,"join") ) { rule->merger = info_rules_merge_join; is_join = 1; }
        else error("The rule logic \"%s\" not recognised\n", ss);

        if ( !is_join && rule->type==BCF_HT_STR )
            error("Numeric operation \"%s\" requested on non-numeric field: %s\n", ss, rule->hdr_tag);
        if ( bcf_hdr_id2number(args->out_hdr,BCF_HL_INFO,id)==0xfffff )
        {
            int is_agr = (
                    bcf_hdr_id2length(args->out_hdr,BCF_HL_INFO,id)==BCF_VL_A ||
                    bcf_hdr_id2length(args->out_hdr,BCF_HL_INFO,id)==BCF_VL_G ||
                    bcf_hdr_id2length(args->out_hdr,BCF_HL_INFO,id)==BCF_VL_R
                    ) ? 1 : 0;
            if ( is_join && is_agr )
                error("Cannot -i %s:join on Number=[AGR] tags is not supported.\n", rule->hdr_tag);
            if ( !is_join && !is_agr )
                error("Only fixed-length vectors are supported with -i %s:%s\n", ss, rule->hdr_tag);
        }

        ss = strchr(ss, '\0'); ss++;
        n++;
    }
    free(str.s);
    free(tmp);

    qsort(args->rules, args->nrules, sizeof(*args->rules), info_rules_comp_key2);
}
static void info_rules_destroy(args_t *args)
{
    int i;
    for (i=0; i<args->nrules; i++)
    {
        info_rule_t *rule = &args->rules[i];
        free(rule->hdr_tag);
        free(rule->vals);
    }
    free(args->rules);
}
static void info_rules_reset(args_t *args)
{
    int i;
    for (i=0; i<args->nrules; i++)
        args->rules[i].nblocks = args->rules[i].nvals = args->rules[i].block_size = 0;
}
static int info_rules_add_values(args_t *args, bcf_hdr_t *hdr, bcf1_t *line, info_rule_t *rule, maux1_t *als, int var_len)
{
    int msize = args->maux->ntmp_arr / rule->type_size;
    int ret = bcf_get_info_values(hdr, line, rule->hdr_tag, &args->maux->tmp_arr, &msize, rule->type);
    if ( ret<=0 ) error("FIXME: error parsing %s at %s:%d .. %d\n", rule->hdr_tag,bcf_seqname(hdr,line),line->pos+1,ret);
    args->maux->ntmp_arr = msize * rule->type_size;

    rule->nblocks++;

    if ( rule->type==BCF_HT_STR )
    {
        int need_comma = rule->nblocks==1 ? 0 : 1;
        hts_expand(char,rule->nvals+ret+need_comma+1,rule->mvals,rule->vals);            // 1 for null-termination
        char *tmp = (char*) rule->vals + rule->nvals;
        if ( rule->nvals>0 ) { *tmp = ','; tmp++; }
        strncpy(tmp,(char*)args->maux->tmp_arr,ret);
        rule->nvals += ret + need_comma;
        return 1;
    }

    int i, j;
    if ( var_len==BCF_VL_A )
    {
        if ( ret!=line->n_allele-1 ) error("Wrong number of %s fields at %s:%d\n",rule->hdr_tag,bcf_seqname(hdr,line),line->pos+1);
        args->maux->nagr_map = ret;
        hts_expand(int,args->maux->nagr_map,args->maux->magr_map,args->maux->agr_map);
        // create mapping from source file ALT indexes to dst file indexes
        for (i=0; i<ret; i++) args->maux->agr_map[i] = als->map[i+1] - 1;
        rule->block_size = args->maux->nout_als - 1;
    }
    else if ( var_len==BCF_VL_R )
    {
        if ( ret!=line->n_allele ) error("Wrong number of %s fields at %s:%d\n",rule->hdr_tag,bcf_seqname(hdr,line),line->pos+1);
        args->maux->nagr_map = ret;
        hts_expand(int,args->maux->nagr_map,args->maux->magr_map,args->maux->agr_map);
        for (i=0; i<ret; i++) args->maux->agr_map[i] = als->map[i];
        rule->block_size = args->maux->nout_als;
    }
    else if ( var_len==BCF_VL_G )
    {
        args->maux->nagr_map = bcf_alleles2gt(line->n_allele-1,line->n_allele-1)+1;
        assert( ret==line->n_allele || ret==args->maux->nagr_map );
        if ( ret==line->n_allele ) // haploid
        {
            args->maux->nagr_map = line->n_allele;
            hts_expand(int,args->maux->nagr_map,args->maux->magr_map,args->maux->agr_map);
            for (i=0; i<ret; i++) args->maux->agr_map[i] = als->map[i];
            rule->block_size = args->maux->nout_als;
        }
        else
        {
            hts_expand(int,args->maux->nagr_map,args->maux->magr_map,args->maux->agr_map);
            int k_src = 0;
            for (i=0; i<line->n_allele; i++)
            {
                for (j=0; j<=i; j++)
                {
                    args->maux->agr_map[k_src] = bcf_alleles2gt(als->map[i],als->map[j]);
                    k_src++;
                }
            }
            rule->block_size = bcf_alleles2gt(args->maux->nout_als-1,args->maux->nout_als-1)+1;
        }
    }
    else
    {
        if ( rule->nblocks>1 && ret!=rule->block_size )
            error("Mismatch in number of values for INFO/%s at %s:%d\n", rule->hdr_tag,bcf_seqname(hdr,line),line->pos+1);
        rule->block_size = ret;
        args->maux->nagr_map = 0;
    }

    #define BRANCH(src_type_t,dst_type_t,set_missing) { \
        src_type_t *src = (src_type_t *) args->maux->tmp_arr; \
        hts_expand0(dst_type_t,(rule->nvals+rule->block_size),rule->mvals,rule->vals); \
        dst_type_t *dst = (dst_type_t *) rule->vals + rule->nvals; \
        rule->nvals += rule->block_size; \
        if ( !args->maux->nagr_map ) \
        { \
            for (i=0; i<ret; i++) dst[i] = src[i]; \
        } \
        else \
        { \
            for (i=0; i<rule->block_size; i++) set_missing; \
            for (i=0; i<ret; i++) dst[args->maux->agr_map[i]] = src[i]; \
        } \
    }
    switch (rule->type) {
        case BCF_HT_INT:  BRANCH(int, int32_t, dst[i] = bcf_int32_missing); break;
        case BCF_HT_REAL: BRANCH(float, float, bcf_float_set_missing(dst[i])); break;
        default: error("TODO: %s:%d .. type=%d\n", __FILE__,__LINE__, rule->type);
    }
    #undef BRANCH

    return 1;
}

int bcf_hdr_sync(bcf_hdr_t *h);

void merge_headers(bcf_hdr_t *hw, const bcf_hdr_t *hr, const char *clash_prefix, int force_samples)
{
    // header lines
    hw = bcf_hdr_merge(hw, hr);

    // samples
    int i;
    for (i=0; i<bcf_hdr_nsamples(hr); i++)
    {
        char *name = hr->samples[i];
        if ( bcf_hdr_id2int(hw, BCF_DT_SAMPLE, name)!=-1 )
        {
            // there is a sample with the same name
            if ( !force_samples ) error("Error: Duplicate sample names (%s), use --force-samples to proceed anyway.\n", name);

            int len = strlen(hr->samples[i]) + strlen(clash_prefix) + 1;
            name = (char*) malloc(sizeof(char)*(len+1));
            sprintf(name,"%s:%s",clash_prefix,hr->samples[i]);
            bcf_hdr_add_sample(hw,name);
            free(name);
        }
        else
            bcf_hdr_add_sample(hw,name);
    }
}

void debug_als(char **als, int nals)
{
    int k; for (k=0; k<nals; k++) fprintf(stderr,"%s ", als[k]);
    fprintf(stderr,"\n");
}

/**
 * normalize_alleles() - create smallest possible representation of the alleles
 * @als:    alleles to be merged, first is REF (rw)
 * @nals:   number of $a alleles
 *
 * Best explained on an example:
 *      In:  REF=GTTT  ALT=GTT
 *      Out: REF=GT    ALT=G
 *
 * Note: the als array will be modified
 */
void normalize_alleles(char **als, int nals)
{
    if ( !als[0][1] ) return;   // ref is 1base long, we're done

    int j, i = 1, done = 0;
    int *lens = (int*) malloc(sizeof(int)*nals);
    for (j=0; j<nals; j++) lens[j] = strlen(als[j]);

    while ( i<lens[0] )
    {
        for (j=1; j<nals; j++)
        {
            if ( i>=lens[j] ) done = 1;
            if ( als[j][lens[j]-i] != als[0][lens[0]-i] ) { done = 1; break; }
        }
        if ( done ) break;
        i++;
    }
    if ( i>1 )
    {
        i--;
        als[0][lens[0]-i] = 0;
        for (j=1; j<nals; j++) als[j][lens[j]-i] = 0;
    }
    free(lens);
}

 /**
 * merge_alleles() - merge two REF,ALT records, $a and $b into $b.
 * @a:      alleles to be merged, first is REF
 * @na:     number of $a alleles
 * @map:    map from the original $a indexes to new $b indexes (0-based)
 * @b:      alleles to be merged, the array will be expanded as required
 * @nb:     number of $b alleles
 * @mb:     size of $b
 *
 * Returns NULL on error or $b expanded to incorporate $a alleles and sets
 * $map. Best explained on an example:
 *      In:     REF   ALT
 *           a: ACG,  AC,A    (1bp and 2bp deletion)
 *           b: ACGT, A       (3bp deletion)
 *      Out:
 *           b: ACGT, A,ACT,AT (3bp, 1bp and 2bp deletion)
 *           map: 0,2,3
 * Here the mapping from the original $a alleles to the new $b alleles is 0->0,
 * 1->2, and 2->3.
 */
char **merge_alleles(char **a, int na, int *map, char **b, int *nb, int *mb)
{
    // reference allele never changes
    map[0] = 0;

    int i,j;
    int rla = !a[0][1] ? 1 : strlen(a[0]);
    int rlb = !b[0][1] ? 1 : strlen(b[0]);

    // the most common case: same SNPs
    if ( na==2 && *nb==2 && rla==1 && rlb==1 && a[1][0]==b[1][0] && !a[1][1] && !b[1][1] )
    {
        map[1] = 1;
        return b;
    }

    // Sanity check: reference prefixes must be identical
    if ( strncmp(a[0],b[0],rla<rlb?rla:rlb) )
    {
        if ( strncasecmp(a[0],b[0],rla<rlb?rla:rlb) )
        {
            fprintf(stderr, "The REF prefixes differ: %s vs %s (%d,%d)\n", a[0],b[0],rla,rlb);
            return NULL;
        }
        // Different case, change to uppercase
        for (i=0; i<na; i++)
        {
            int len = strlen(a[i]);
            for (j=0; j<len; j++) a[i][j] = toupper(a[i][j]);
        }
        for (i=0; i<*nb; i++)
        {
            int len = strlen(b[i]);
            for (j=0; j<len; j++) b[i][j] = toupper(b[i][j]);
        }
    }

    int n = *nb + na;
    hts_expand0(char*,n,*mb,b);

    // $b alleles need expanding
    if ( rla>rlb )
    {
        for (i=0; i<*nb; i++)
        {
            if ( b[i][0]=='<' ) continue;   // symbolic allele, do not modify
            if ( b[i][0]=='*' ) continue;   // overlapping deletion (*), do not modify
            int l = strlen(b[i]);
            b[i] = (char*) realloc(b[i],l+rla-rlb+1);
            memcpy(b[i]+l,a[0]+rlb,rla-rlb+1);
        }
    }

    // now check if the $a alleles are present and if not add them
    for (i=1; i<na; i++)
    {
        int const_ai = 1;
        char *ai;
        if ( rlb>rla && a[i][0]!='<' && a[i][0]!='*' )  // $a alleles need expanding and not a symbolic allele or *
        {
            int l = strlen(a[i]);
            ai = (char*) malloc(l+rlb-rla+1);
            memcpy(ai,a[i],l);
            memcpy(ai+l,b[0]+rla,rlb-rla+1);
            const_ai = 0;
        }
        else
            ai = a[i];

        for (j=1; j<*nb; j++)
            if ( !strcasecmp(ai,b[j]) ) break;

        if ( j<*nb ) // $b already has the same allele
        {
            map[i] = j;
            if ( !const_ai ) free(ai);
            continue;
        }
        // new allele
        map[i] = *nb;
        b[*nb] = const_ai ? strdup(ai) : ai;
        (*nb)++;
    }
    return b;
}

maux_t *maux_init(args_t *args)
{
    bcf_srs_t *files = args->files;
    maux_t *ma = (maux_t*) calloc(1,sizeof(maux_t));
    ma->n      = files->nreaders;
    ma->files  = files;
    int i, n_smpl = 0;
    for (i=0; i<ma->n; i++)
        n_smpl += bcf_hdr_nsamples(files->readers[i].header);
    if ( args->do_gvcf )
    {
        ma->gvcf = (gvcf_aux_t*) calloc(ma->n,sizeof(gvcf_aux_t));
        for (i=0; i<ma->n; i++)
            ma->gvcf[i].line = bcf_init1();
    }
    ma->smpl_ploidy = (int*) calloc(n_smpl,sizeof(int));
    ma->smpl_nGsize = (int*) malloc(n_smpl*sizeof(int));
    ma->buf = (buffer_t*) calloc(ma->n,sizeof(buffer_t));
    for (i=0; i<ma->n; i++)
        ma->buf[i].rid = -1;
    return ma;
}
void maux_destroy(maux_t *ma)
{
    int i,j;
    for (i=0; i<ma->mals; i++)
    {
        free(ma->als[i]);
        ma->als[i] = NULL;
    }
    for (i=0; i<ma->n; i++) // for each reader
    {
        for (j=0; j<ma->buf[i].mrec; j++)  // for each buffered line
            free(ma->buf[i].rec[j].map);
        free(ma->buf[i].rec);
    }
    free(ma->buf);
    if ( ma->gvcf )
    {
        for (i=0; i<ma->n; i++) bcf_destroy(ma->gvcf[i].line);
        free(ma->gvcf);
    }
    for (i=0; i<ma->mAGR_info; i++)
        free(ma->AGR_info[i].buf);
    free(ma->agr_map);
    free(ma->AGR_info);
    if (ma->ntmp_arr) free(ma->tmp_arr);
    if (ma->nfmt_map) free(ma->fmt_map);
    // ma->inf freed in bcf_destroy1
    for (i=0; i<ma->mals; i++) free(ma->als[i]);
    if (ma->mout_als) free(ma->out_als);
    free(ma->als);
    free(ma->cnt);
    free(ma->smpl_ploidy);
    free(ma->smpl_nGsize);
    free(ma->chr);
    free(ma);
}
void maux_expand1(buffer_t *buf, int size)
{
    if ( buf->mrec < size )
    {
        hts_expand0(maux1_t,size,buf->mrec,buf->rec);
        buf->mrec = size;
    }
}
void maux_reset(maux_t *ma)
{
    int i,j;
    for (i=0; i<ma->n; i++) maux_expand1(&ma->buf[i],ma->files->readers[i].nbuffer+1);
    for (i=0; i<ma->ncnt; i++) ma->cnt[i] = 0;
    for (i=0; i<ma->mals; i++)
    {
        free(ma->als[i]);
        ma->als[i] = NULL;
    }
    const char *chr = NULL;
    ma->nals  = 0;
    ma->pos   = -1;
    for (i=0; i<ma->n; i++)
    {
        if ( !bcf_sr_has_line(ma->files,i) ) continue;
        bcf1_t *line = bcf_sr_get_line(ma->files,i);
        bcf_hdr_t *hdr = bcf_sr_get_header(ma->files,i);
        chr = bcf_seqname(hdr,line);
        ma->pos = line->pos;
        break;
    }
    int new_chr = 0;
    if ( chr && (!ma->chr || strcmp(ma->chr,chr)) )
    {
        free(ma->chr);
        ma->chr = strdup(chr);
        new_chr = 1;
    }
    for (i=0; i<ma->n; i++)
    {
        bcf_hdr_t *hdr = bcf_sr_get_header(ma->files,i);
        ma->buf[i].rid = bcf_hdr_name2id(hdr,chr);
        ma->buf[i].beg = bcf_sr_has_line(ma->files,i) ? 0 : 1;
        for (j=ma->buf[i].beg; j<=ma->files->readers[i].nbuffer; j++)
        {
            ma->buf[i].rec[j].skip = 0;
            bcf1_t *line = ma->files->readers[i].buffer[j];
            if ( line->rid!=ma->buf[i].rid || line->pos!=ma->pos ) break;
        }
        ma->buf[i].end = j;
        ma->buf[i].cur = -1;
        if ( ma->buf[i].beg < ma->buf[i].end ) 
        {
            ma->buf[i].lines = ma->files->readers[i].buffer;
            if ( ma->gvcf ) ma->gvcf[i].active = 0;     // gvcf block cannot overlap with the next record
        }
        if ( new_chr && ma->gvcf ) ma->gvcf[i].active = 0;  // make sure to close active gvcf block on new chr
    }
}
void maux_debug(maux_t *ma, int ir, int ib)
{
    printf("[%d,%d]\t", ir,ib);
    int i;
    for (i=0; i<ma->nals; i++)
    {
        printf(" %s [%d]", ma->als[i], ma->cnt[i]);
    }
    printf("\n");
}

void merge_chrom2qual(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;

    int i, ret;
    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);
    kstring_t *tmps = &args->tmps;
    tmps->l = 0;

    maux_t *ma = args->maux;
    int *al_idxs = (int*) calloc(ma->nals,sizeof(int));
    bcf_float_set_missing(out->qual);

    // CHROM, POS, ID, QUAL
    out->pos = -1;
    for (i=0; i<files->nreaders; i++)
    {
        bcf1_t *line = maux_get_line(args, i);
        if ( !line ) continue;
        bcf_unpack(line, BCF_UN_ALL);

        bcf_sr_t *reader = &files->readers[i];
        bcf_hdr_t *hdr = reader->header;

        // not all maux alleles are always used, mark the ones we'll need
        int j;
        for (j=1; j<line->n_allele; j++)
        {
            int irec = ma->buf[i].cur;
            al_idxs[ ma->buf[i].rec[irec].map[j] ] = 1;
        }

        // position
        if ( out->pos==-1 )
        {
            const char *chr = hdr->id[BCF_DT_CTG][line->rid].key;
            out->rid = bcf_hdr_name2id(out_hdr, chr);
            if ( strcmp(chr,out_hdr->id[BCF_DT_CTG][out->rid].key) ) error("Uh\n");
            out->pos = line->pos;
        }

        // ID
        if ( line->d.id[0]!='.' || line->d.id[1] )
        {
            kitr = kh_get(strdict, tmph, line->d.id);
            if ( kitr == kh_end(tmph) )
            {
                if ( tmps->l ) kputc(';', tmps);
                kputs(line->d.id, tmps);
                kh_put(strdict, tmph, line->d.id, &ret);
            }
        }

        // set QUAL to the max qual value. Not exactly correct, but good enough for now
        if ( !bcf_float_is_missing(line->qual) )
        {
            if ( bcf_float_is_missing(out->qual) || out->qual < line->qual ) out->qual = line->qual;
        }
    }

    // set ID
    if ( !tmps->l ) kputs(".", tmps);
    bcf_update_id(out_hdr, out, tmps->s);

    // set alleles
    ma->nout_als = 0;
    for (i=1; i<ma->nals; i++)
    {
        if ( !al_idxs[i] ) continue;
        ma->nout_als++;

        // Adjust the indexes, the allele map could be created for multiple collapsed records,
        //  some of which might be unused for this output line
        int ir, j;
        for (ir=0; ir<files->nreaders; ir++)
        {
            bcf1_t *line = maux_get_line(args,ir);
            if ( !line ) continue;
            for (j=1; j<line->n_allele; j++)
            {
                int irec = ma->buf[ir].cur;
                if ( ma->buf[ir].rec[irec].map[j]==i ) ma->buf[ir].rec[irec].map[j] = ma->nout_als;
            }
        }
    }
    // Expand the arrays and realloc the alleles string. Note that all alleles are in a single allocated block.
    ma->nout_als++;
    hts_expand0(char*, ma->nout_als, ma->mout_als, ma->out_als);
    int k = 0;
    for (i=0; i<ma->nals; i++)
        if ( i==0 || al_idxs[i] ) ma->out_als[k++] = strdup(ma->als[i]);
    assert( k==ma->nout_als );
    normalize_alleles(ma->out_als, ma->nout_als);
    bcf_update_alleles(out_hdr, out, (const char**) ma->out_als, ma->nout_als);
    free(al_idxs);
    for (i=0; i<ma->nout_als; i++) free(ma->out_als[i]);
}

void merge_filter(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;

    int i, ret;
    if ( args->filter_logic == FLT_LOGIC_REMOVE )
    {
        for (i=0; i<files->nreaders; i++)
        {
            bcf1_t *line = maux_get_line(args, i);
            if ( !line ) continue;
            bcf_sr_t *reader = &files->readers[i];
            bcf_hdr_t *hdr = reader->header;
            if ( bcf_has_filter(hdr, line, "PASS") ) break;
        }
        if ( i<files->nreaders )
        {
            int flt_id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, "PASS");
            bcf_add_filter(out_hdr, out, flt_id);
            return;
        }
    }

    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);

    out->d.n_flt = 0;
    for (i=0; i<files->nreaders; i++)
    {
        bcf1_t *line = maux_get_line(args, i);
        if ( !line ) continue;

        bcf_sr_t *reader = &files->readers[i];
        bcf_hdr_t *hdr = reader->header;

        int k;
        for (k=0; k<line->d.n_flt; k++)
        {
            const char *flt = hdr->id[BCF_DT_ID][line->d.flt[k]].key;
            kitr = kh_get(strdict, tmph, flt);
            if ( kitr == kh_end(tmph) )
            {
                int id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, flt);
                if ( id==-1 ) error("Error: The filter is not defined in the header: %s\n", flt);
                hts_expand(int,out->d.n_flt+1,out->d.m_flt,out->d.flt);
                out->d.flt[out->d.n_flt] = id;
                out->d.n_flt++;
                kh_put(strdict, tmph, flt, &ret);
            }
        }
    }
    // Check if PASS is not mixed with other filters
    if ( out->d.n_flt>1 )
    {
        int id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, "PASS");
        for (i=0; i<out->d.n_flt; i++)
            if ( out->d.flt[i]==id ) break;
        if ( i<out->d.n_flt )
        {
            out->d.n_flt--;
            for (; i<out->d.n_flt; i++) out->d.flt[i] = out->d.flt[i+1];
        }
    }
}

static void bcf_info_set_id(bcf1_t *line, bcf_info_t *info, int id, kstring_t *tmp_str)
{
    uint8_t *ptr = info->vptr - info->vptr_off;
    bcf_dec_typed_int1(ptr, &ptr);

    tmp_str->l = 0;
    bcf_enc_int1(tmp_str, id);

    if ( tmp_str->l == ptr - info->vptr + info->vptr_off )
    {
        // the new id is represented with the same number of bytes
        memcpy(info->vptr - info->vptr_off, tmp_str->s, tmp_str->l);
        return;
    }

    kputsn_(ptr, info->vptr - ptr, tmp_str);
    info->vptr_off = tmp_str->l;
    kputsn_(info->vptr, info->len << bcf_type_shift[info->type], tmp_str);

    info->vptr = (uint8_t*) tmp_str->s + info->vptr_off;
    tmp_str->s = NULL;
    tmp_str->m = 0;
    tmp_str->l = 0;
}

/*
 *  copy_string_field() - copy a comma-separated field
 *  @param src:     source string
 *  @param isrc:    index of the field to copy 
 *  @param src_len: length of source string (excluding the terminating \0) 
 *  @param dst:     destination kstring (must be initialized)
 *  @param idst:    index of the destination field
 */
int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst)
{
    int ith_src = 0, start_src = 0;    // i-th field in src string
    while ( ith_src<isrc && start_src<src_len )
    {
        if ( src[start_src]==',' ) { ith_src++; }
        start_src++;
    }
    if ( ith_src!=isrc ) return -1; // requested field not found
    int end_src = start_src;
    while ( end_src<src_len && src[end_src] && src[end_src]!=',' ) end_src++;

    int nsrc_cpy = end_src - start_src;
    if ( nsrc_cpy==1 && src[start_src]=='.' ) return 0;   // don't write missing values, dst is already initialized

    int ith_dst = 0, start_dst = 0;
    while ( ith_dst<idst && start_dst<dst->l )
    {
        if ( dst->s[start_dst]==',' ) { ith_dst++; }
        start_dst++;
    }
    if ( ith_dst!=idst ) return -2;
    int end_dst = start_dst;
    while ( end_dst<dst->l && dst->s[end_dst]!=',' ) end_dst++;

    if ( end_dst - start_dst>1 || dst->s[start_dst]!='.' ) return 0;   // do not overwrite non-empty values

    // Now start_dst and end_dst are indexes to the destination memory area
    // which needs to be replaced with nsrc_cpy
    // source bytes, end_dst points just after.
    int ndst_shift = nsrc_cpy - (end_dst - start_dst);
    int ndst_move  = dst->l - end_dst + 1;  // how many bytes must be moved (including \0)
    if ( ndst_shift )
    {
        ks_resize(dst, dst->l + ndst_shift + 1);    // plus \0
        memmove(dst->s+end_dst+ndst_shift, dst->s+end_dst, ndst_move);
    }
    memcpy(dst->s+start_dst, src+start_src, nsrc_cpy);
    dst->l += ndst_shift;
    return 0;
}

static void merge_AGR_info_tag(bcf_hdr_t *hdr, bcf1_t *line, bcf_info_t *info, int len, maux1_t *als, AGR_info_t *agr)
{
    int i;
    if ( !agr->nbuf )
    {
        if ( info->type==BCF_BT_INT8 || info->type==BCF_BT_INT16 || info->type==BCF_BT_INT32 || info->type==BCF_BT_FLOAT )
        {
            agr->nbuf = 4 * agr->nvals;
            hts_expand(uint8_t,agr->nbuf,agr->mbuf,agr->buf);
            if ( info->type!=BCF_BT_FLOAT )
            {
                int32_t *tmp = (int32_t*) agr->buf;
                for (i=0; i<agr->nvals; i++) tmp[i] = bcf_int32_missing;
            }
            else
            {
                float *tmp = (float*) agr->buf;
                for (i=0; i<agr->nvals; i++) bcf_float_set_missing(tmp[i]);
            }
        }
        else if ( info->type==BCF_BT_CHAR )
        {
            kstring_t tmp; tmp.l = 0; tmp.m = agr->mbuf; tmp.s = (char*)agr->buf;
            kputc('.',&tmp);
            for (i=1; i<agr->nvals; i++) kputs(",.",&tmp);
            agr->mbuf = tmp.m; agr->nbuf = tmp.l; agr->buf = (uint8_t*)tmp.s;
        }
        else
            error("Not ready for type [%d]: %s at %d\n", info->type,agr->hdr_tag,line->pos+1);
    }

    if ( info->type==BCF_BT_INT8 || info->type==BCF_BT_INT16 || info->type==BCF_BT_INT32 || info->type==BCF_BT_FLOAT )
    {
        if ( len==BCF_VL_A || len==BCF_VL_R )
        {
            int ifrom = len==BCF_VL_A ? 1 : 0;
            #define BRANCH(type_t, is_missing, is_vector_end, out_type_t) { \
                type_t *src = (type_t *) info->vptr; \
                out_type_t *tgt = (out_type_t *) agr->buf; \
                int iori, inew; \
                for (iori=ifrom; iori<line->n_allele; iori++) \
                { \
                    if ( is_vector_end ) break; \
                    if ( is_missing ) continue; \
                    inew = als->map[iori] - ifrom; \
                    tgt[inew] = *src; \
                    src++; \
                } \
            }
            switch (info->type) {
                case BCF_BT_INT8:  BRANCH(int8_t,  *src==bcf_int8_missing,  *src==bcf_int8_vector_end,  int); break;
                case BCF_BT_INT16: BRANCH(int16_t, *src==bcf_int16_missing, *src==bcf_int16_vector_end, int); break;
                case BCF_BT_INT32: BRANCH(int32_t, *src==bcf_int32_missing, *src==bcf_int32_vector_end, int); break;
                case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(*src), bcf_float_is_vector_end(*src), float); break;
                default: fprintf(stderr,"TODO: %s:%d .. info->type=%d\n", __FILE__,__LINE__, info->type); exit(1);
            }
            #undef BRANCH
        }
        else
        {
            #define BRANCH(type_t, is_missing, is_vector_end, out_type_t) { \
                type_t *src = (type_t *) info->vptr; \
                out_type_t *tgt = (out_type_t *) agr->buf; \
                int iori,jori, inew,jnew; \
                for (iori=0; iori<line->n_allele; iori++) \
                { \
                    inew = als->map[iori]; \
                    for (jori=0; jori<=iori; jori++) \
                    { \
                        jnew = als->map[jori]; \
                        int kori = iori*(iori+1)/2 + jori; \
                        if ( is_vector_end ) break; \
                        if ( is_missing ) continue; \
                        int knew = inew>jnew ? inew*(inew+1)/2 + jnew : jnew*(jnew+1)/2 + inew; \
                        tgt[knew] = src[kori]; \
                    } \
                    if ( jori<=iori ) break; \
                } \
            }
            switch (info->type) {
                case BCF_BT_INT8:  BRANCH(int8_t,  src[kori]==bcf_int8_missing,  src[kori]==bcf_int8_vector_end,  int); break;
                case BCF_BT_INT16: BRANCH(int16_t, src[kori]==bcf_int16_missing, src[kori]==bcf_int16_vector_end, int); break;
                case BCF_BT_INT32: BRANCH(int32_t, src[kori]==bcf_int32_missing, src[kori]==bcf_int32_vector_end, int); break;
                case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(src[kori]), bcf_float_is_vector_end(src[kori]), float); break;
                default: fprintf(stderr,"TODO: %s:%d .. info->type=%d\n", __FILE__,__LINE__, info->type); exit(1);
            }
            #undef BRANCH
        }
    }
    else
    {
        kstring_t tmp; tmp.l = agr->nbuf; tmp.m = agr->mbuf; tmp.s = (char*)agr->buf;
        if ( len==BCF_VL_A || len==BCF_VL_R )
        {
            int iori, ifrom = len==BCF_VL_A ? 1 : 0;
            for (iori=ifrom; iori<line->n_allele; iori++)
            {
                int ret = copy_string_field((char*)info->vptr, iori-ifrom, info->len, &tmp, als->map[iori]-ifrom);
                if ( ret )
                    error("Error at %s:%d: wrong number of fields in %s?\n", bcf_seqname(hdr,line),line->pos+1,agr->hdr_tag);
            }
        }
        else
        {
            int iori,jori, inew,jnew;
            for (iori=0; iori<line->n_allele; iori++)
            {
                inew = als->map[iori];
                for (jori=0; jori<=iori; jori++)
                {
                    jnew = als->map[jori];
                    int kori = iori*(iori+1)/2 + jori;
                    int knew = bcf_alleles2gt(inew,jnew);
                    int ret  = copy_string_field((char*)info->vptr, kori, info->len, &tmp, knew);
                    if ( ret )
                        error("Error at %s:%d: wrong number of fields in %s?\n", bcf_seqname(hdr,line),line->pos+1,agr->hdr_tag);
                }
            }
        }
        agr->mbuf = tmp.m; agr->nbuf = tmp.l; agr->buf = (uint8_t*)tmp.s;
    }
}

void merge_info(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;

    int i, j, ret;
    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);

    maux_t *ma = args->maux;
    ma->nAGR_info = 0;
    out->n_info   = 0;
    info_rules_reset(args);
    for (i=0; i<files->nreaders; i++)
    {
        bcf1_t *line = maux_get_line(args,i);
        if ( !line ) continue;
        int irec = ma->buf[i].cur;
        bcf_sr_t *reader = &files->readers[i];
        bcf_hdr_t *hdr = reader->header;
        for (j=0; j<line->n_info; j++)
        {
            bcf_info_t *inf = &line->d.info[j];

            const char *key = hdr->id[BCF_DT_ID][inf->key].key;
            if ( !strcmp("AC",key) || !strcmp("AN",key) ) continue;  // AC and AN are done in merge_format() after genotypes are done

            int id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, key);
            if ( id==-1 ) error("Error: The INFO field is not defined in the header: %s\n", key);

            kitr = kh_get(strdict, tmph, key);  // have we seen the tag in one of the readers?
            int len = bcf_hdr_id2length(hdr,BCF_HL_INFO,inf->key);
            if ( args->nrules )
            {
                info_rule_t *rule = (info_rule_t*) bsearch(key, args->rules, args->nrules, sizeof(*args->rules), info_rules_comp_key);
                if ( rule )
                {
                    maux1_t *als = ( len==BCF_VL_A || len==BCF_VL_G || len==BCF_VL_R ) ? &ma->buf[i].rec[irec] : NULL;
                    if ( info_rules_add_values(args, hdr, line, rule, als, len) ) continue;
                }
            }

            // Todo: Number=AGR tags should use the newer info_rules_* functions (info_rules_merge_first to be added)
            // and merge_AGR_info_tag to be made obsolete.
            if ( len==BCF_VL_A || len==BCF_VL_G || len==BCF_VL_R  ) // Number=R,G,A requires special treatment
            {
                if ( kitr == kh_end(tmph) )
                {
                    // seeing this key for the first time
                    ma->nAGR_info++;
                    hts_expand0(AGR_info_t,ma->nAGR_info,ma->mAGR_info,ma->AGR_info);
                    kitr = kh_put(strdict, tmph, key, &ret);
                    kh_val(tmph,kitr) = ma->nAGR_info - 1;
                    ma->AGR_info[ma->nAGR_info-1].hdr_tag = key;
                    ma->AGR_info[ma->nAGR_info-1].type  = bcf_hdr_id2type(hdr,BCF_HL_INFO,inf->key);
                    ma->AGR_info[ma->nAGR_info-1].nbuf  = 0;    // size of the buffer
                    switch (len)
                    {
                        case BCF_VL_A: ma->AGR_info[ma->nAGR_info-1].nvals = ma->nout_als - 1; break;
                        case BCF_VL_G: ma->AGR_info[ma->nAGR_info-1].nvals = bcf_alleles2gt(ma->nout_als-1,ma->nout_als-1)+1; break;
                        case BCF_VL_R: ma->AGR_info[ma->nAGR_info-1].nvals = ma->nout_als; break;
                    }
                }
                kitr = kh_get(strdict, tmph, key);
                int idx = kh_val(tmph, kitr);
                if ( idx<0 ) error("Error occurred while processing INFO tag \"%s\" at %s:%d\n", key,bcf_seqname(hdr,line),line->pos+1);
                merge_AGR_info_tag(hdr, line,inf,len,&ma->buf[i].rec[irec],&ma->AGR_info[idx]);
                continue;
            }

            if ( kitr == kh_end(tmph) )
            {
                // Seeing this key for the first time.  Although quite hacky,
                // this is faster than anything else given the data structures..

                hts_expand0(bcf_info_t,out->n_info+1,out->d.m_info,out->d.info);
                out->d.info[out->n_info].key  = id;
                out->d.info[out->n_info].type = inf->type;
                out->d.info[out->n_info].len  = inf->len;
                out->d.info[out->n_info].v1.i = inf->v1.i;
                out->d.info[out->n_info].v1.f = inf->v1.f;
                out->d.info[out->n_info].vptr_off  = inf->vptr_off;
                out->d.info[out->n_info].vptr_len  = inf->vptr_len;
                out->d.info[out->n_info].vptr_free = 1;
                out->d.info[out->n_info].vptr = (uint8_t*) malloc(inf->vptr_len+inf->vptr_off); 
                memcpy(out->d.info[out->n_info].vptr,inf->vptr-inf->vptr_off, inf->vptr_len+inf->vptr_off);
                out->d.info[out->n_info].vptr += inf->vptr_off;
                if ( (args->output_type & FT_BCF) && id!=bcf_hdr_id2int(hdr, BCF_DT_ID, key) )
                    bcf_info_set_id(out, &out->d.info[out->n_info], id, &args->tmps);
                out->d.shared_dirty |= BCF1_DIRTY_INF;
                out->n_info++;
                kitr = kh_put(strdict, tmph, key, &ret);
                kh_val(tmph,kitr) = -(out->n_info-1);   // arbitrary negative value
            }
        }
    }
    for (i=0; i<args->nrules; i++)
        args->rules[i].merger(args->out_hdr, out, &args->rules[i]);
    for (i=0; i<ma->nAGR_info; i++)
    {
        AGR_info_t *agr = &ma->AGR_info[i];
        bcf_update_info(out_hdr,out,agr->hdr_tag,agr->buf,agr->nvals,agr->type);
    }
}

void update_AN_AC(bcf_hdr_t *hdr, bcf1_t *line)
{
    int32_t an = 0, *tmp = (int32_t*) malloc(sizeof(int)*line->n_allele);
    int ret = bcf_calc_ac(hdr, line, tmp, BCF_UN_FMT);
    if ( ret>0 )
    {
        int i;
        for (i=0; i<line->n_allele; i++) an += tmp[i];
        bcf_update_info_int32(hdr, line, "AN", &an, 1);
        bcf_update_info_int32(hdr, line, "AC", tmp+1, line->n_allele-1);
    }
    free(tmp);
}

void merge_GT(args_t *args, bcf_fmt_t **fmt_map, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    int i, ismpl = 0, nsamples = bcf_hdr_nsamples(out_hdr);

    int nsize = 0, msize = sizeof(int32_t);
    for (i=0; i<files->nreaders; i++)
    {
        if ( !fmt_map[i] ) continue;
        if ( fmt_map[i]->n > nsize ) nsize = fmt_map[i]->n;
    }

    if ( ma->ntmp_arr < nsamples*nsize*msize )
    {
        ma->ntmp_arr = nsamples*nsize*msize;
        ma->tmp_arr  = realloc(ma->tmp_arr, ma->ntmp_arr);
    }
    memset(ma->smpl_ploidy,0,nsamples*sizeof(int));

    int default_gt = args->missing_to_ref ? bcf_gt_unphased(0) : bcf_gt_missing;
    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        bcf_hdr_t *hdr = reader->header;
        bcf_fmt_t *fmt_ori = fmt_map[i];
        int32_t *tmp  = (int32_t *) ma->tmp_arr + ismpl*nsize;
        int irec = ma->buf[i].cur;

        int j, k;
        if ( !fmt_ori )
        {
            // missing values: assume maximum ploidy
            for (j=0; j<bcf_hdr_nsamples(hdr); j++)
            {
                for (k=0; k<nsize; k++) { tmp[k] = default_gt; ma->smpl_ploidy[ismpl+j]++; }
                tmp += nsize;
            }
            ismpl += bcf_hdr_nsamples(hdr);
            continue;
        }

        #define BRANCH(type_t, vector_end) { \
            type_t *p_ori  = (type_t*) fmt_ori->p; \
            if ( !ma->buf[i].rec[irec].als_differ ) \
            { \
                /* the allele numbering is unchanged */ \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    for (k=0; k<fmt_ori->n; k++) \
                    { \
                        if ( p_ori[k]==vector_end ) break; /* smaller ploidy */ \
                        ma->smpl_ploidy[ismpl+j]++; \
                        if ( bcf_gt_is_missing(p_ori[k]) ) tmp[k] = 0; /* missing allele */ \
                        else tmp[k] = p_ori[k]; \
                    } \
                    for (; k<nsize; k++) tmp[k] = bcf_int32_vector_end; \
                    tmp += nsize; \
                    p_ori += fmt_ori->n; \
                } \
                ismpl += bcf_hdr_nsamples(hdr); \
                continue; \
            } \
            /* allele numbering needs to be changed */ \
            for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
            { \
                for (k=0; k<fmt_ori->n; k++) \
                { \
                    if ( p_ori[k]==vector_end ) break; /* smaller ploidy */ \
                    ma->smpl_ploidy[ismpl+j]++; \
                    if ( bcf_gt_is_missing(p_ori[k]) ) tmp[k] = 0; /* missing allele */ \
                    else \
                    { \
                        int al = (p_ori[k]>>1) - 1; \
                        al = al<=0 ? al + 1 : ma->buf[i].rec[irec].map[al] + 1; \
                        tmp[k] = (al << 1) | ((p_ori[k])&1); \
                    } \
                } \
                for (; k<nsize; k++) tmp[k] = bcf_int32_vector_end; \
                tmp += nsize; \
                p_ori += fmt_ori->n; \
            } \
            ismpl += bcf_hdr_nsamples(hdr); \
        }
        switch (fmt_ori->type)
        {
            case BCF_BT_INT8: BRANCH(int8_t,   bcf_int8_vector_end); break;
            case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_vector_end); break;
            case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_vector_end); break;
            default: error("Unexpected case: %d\n", fmt_ori->type);
        }
        #undef BRANCH
    }
    bcf_update_format_int32(out_hdr, out, "GT", (int32_t*)ma->tmp_arr, nsamples*nsize);
}

void merge_format_field(args_t *args, bcf_fmt_t **fmt_map, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    int i, ismpl = 0, nsamples = bcf_hdr_nsamples(out_hdr);

    const char *key = NULL;
    int nsize = 0, length = BCF_VL_FIXED, type = -1;
    for (i=0; i<files->nreaders; i++)
    {
        if ( !maux_get_line(args,i) ) continue;
        if ( !fmt_map[i] ) continue;
        if ( !key ) key = files->readers[i].header->id[BCF_DT_ID][fmt_map[i]->id].key;
        type = fmt_map[i]->type;
        if ( IS_VL_G(files->readers[i].header, fmt_map[i]->id) )
        {
            length = BCF_VL_G;
            nsize = out->n_allele*(out->n_allele + 1)/2;
            break;
        }
        if ( IS_VL_A(files->readers[i].header, fmt_map[i]->id) )
        {
            length = BCF_VL_A;
            nsize = out->n_allele - 1;
            break;
        }
        if ( IS_VL_R(files->readers[i].header, fmt_map[i]->id) )
        {
            length = BCF_VL_R;
            nsize = out->n_allele;
            break;
        }
        if ( fmt_map[i]->n > nsize ) nsize = fmt_map[i]->n;
    }

    int msize = sizeof(float)>sizeof(int32_t) ? sizeof(float) : sizeof(int32_t);
    if ( ma->ntmp_arr < nsamples*nsize*msize )
    {
        ma->ntmp_arr = nsamples*nsize*msize;
        ma->tmp_arr  = realloc(ma->tmp_arr, ma->ntmp_arr);
    }

    // Fill the temp array for all samples by collecting values from all files
    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        bcf_hdr_t *hdr = reader->header;
        bcf_fmt_t *fmt_ori = fmt_map[i];
        bcf1_t *line = maux_get_line(args, i);
        int irec = ma->buf[i].cur;
        if ( fmt_ori )
        {
            type = fmt_ori->type;
            int nals_ori = line->n_allele;
            if ( length==BCF_VL_G )
            {
                // if all fields are missing then n==1 is valid
                if ( fmt_ori->n!=1 && fmt_ori->n != nals_ori*(nals_ori+1)/2 && fmt_map[i]->n != nals_ori )
                    error("Incorrect number of %s fields (%d) at %s:%d, cannot merge.\n", key,fmt_ori->n,bcf_seqname(args->out_hdr,out),out->pos+1);
            }
            else if ( length==BCF_VL_A )
            {
                if ( fmt_ori->n!=1 && fmt_ori->n != nals_ori-1 )
                    error("Incorrect number of %s fields (%d) at %s:%d, cannot merge.\n", key,fmt_ori->n,bcf_seqname(args->out_hdr,out),out->pos+1);
            }
            else if ( length==BCF_VL_R )
            {
                if ( fmt_ori->n!=1 && fmt_ori->n != nals_ori )
                    error("Incorrect number of %s fields (%d) at %s:%d, cannot merge.\n", key,fmt_ori->n,bcf_seqname(args->out_hdr,out),out->pos+1);
            }
        }

        // set the values
        #define BRANCH(tgt_type_t, src_type_t, src_is_missing, src_is_vector_end, tgt_set_missing, tgt_set_vector_end) { \
            int j, l, k; \
            tgt_type_t *tgt = (tgt_type_t *) ma->tmp_arr + ismpl*nsize; \
            if ( !fmt_ori ) \
            { \
                /* the field is not present in this file, set missing values */ \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    tgt_set_missing; tgt++; for (l=1; l<nsize; l++) { tgt_set_vector_end; tgt++; } \
                } \
                ismpl += bcf_hdr_nsamples(hdr); \
                continue; \
            } \
            src_type_t *src = (src_type_t*) fmt_ori->p; \
            if ( (length!=BCF_VL_G && length!=BCF_VL_A && length!=BCF_VL_R) || (line->n_allele==out->n_allele && !ma->buf[i].rec[irec].als_differ) ) \
            { \
                /* alleles unchanged, copy over */ \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    for (l=0; l<fmt_ori->n; l++) \
                    { \
                        if ( src_is_vector_end ) break; \
                        else if ( src_is_missing ) tgt_set_missing; \
                        else *tgt = *src; \
                        tgt++; src++; \
                    } \
                    for (k=l; k<nsize; k++) { tgt_set_vector_end; tgt++; } \
                    src += fmt_ori->n - l; \
                } \
                ismpl += bcf_hdr_nsamples(hdr); \
                continue; \
            } \
            /* allele numbering needs to be changed */ \
            if ( length==BCF_VL_G ) \
            { \
                /* Number=G tags */ \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize; \
                    src = (src_type_t*) fmt_ori->p + j*fmt_ori->n; \
                    if ( (src_is_missing && fmt_ori->n==1) || (++src && src_is_vector_end) ) \
                    { \
                        /* tag with missing value "." */ \
                        tgt_set_missing; \
                        for (l=1; l<nsize; l++) { tgt++; tgt_set_vector_end; } \
                        continue; \
                    } \
                    int ngsize = ma->smpl_ploidy[ismpl+j]==1 ? out->n_allele : out->n_allele*(out->n_allele + 1)/2; \
                    for (l=0; l<ngsize; l++) { tgt_set_missing; tgt++; } \
                    for (; l<nsize; l++) { tgt_set_vector_end; tgt++; } \
                    if ( ma->smpl_ploidy[ismpl+j]==1 ) \
                    { \
                        /* Haploid */ \
                        int iori, inew; \
                        for (iori=0; iori<line->n_allele; iori++) \
                        { \
                            inew = ma->buf[i].rec[irec].map[iori]; \
                            src = (src_type_t*) fmt_ori->p + j*fmt_ori->n + iori; \
                            tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize + inew; \
                            if ( src_is_vector_end ) break; \
                            if ( src_is_missing ) tgt_set_missing; \
                            else *tgt = *src; \
                        } \
                    } \
                    else \
                    { \
                        /* Diploid */ \
                        int iori,jori, inew,jnew; \
                        for (iori=0; iori<line->n_allele; iori++) \
                        { \
                            inew = ma->buf[i].rec[irec].map[iori]; \
                            for (jori=0; jori<=iori; jori++) \
                            { \
                                jnew = ma->buf[i].rec[irec].map[jori]; \
                                int kori = iori*(iori+1)/2 + jori; \
                                int knew = inew>jnew ? inew*(inew+1)/2 + jnew : jnew*(jnew+1)/2 + inew; \
                                src = (src_type_t*) fmt_ori->p + j*fmt_ori->n + kori; \
                                tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize + knew; \
                                if ( src_is_vector_end ) \
                                { \
                                    iori = line->n_allele; \
                                    break; \
                                } \
                                if ( src_is_missing ) tgt_set_missing; \
                                else *tgt = *src; \
                            } \
                        } \
                    } \
                } \
            } \
            else \
            { \
                /* Number=A or Number=R tags */ \
                int ifrom = length==BCF_VL_A ? 1 : 0; \
                for (j=0; j<bcf_hdr_nsamples(hdr); j++) \
                { \
                    tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize; \
                    src = (src_type_t*) (fmt_ori->p + j*fmt_ori->size); \
                    if ( (src_is_missing && fmt_ori->n==1) || (++src && src_is_vector_end) ) \
                    { \
                        /* tag with missing value "." */ \
                        tgt_set_missing; \
                        for (l=1; l<nsize; l++) { tgt++; tgt_set_vector_end; } \
                        continue; \
                    } \
                    src = (src_type_t*) (fmt_ori->p + j*fmt_ori->size); \
                    for (l=0; l<nsize; l++) { tgt_set_missing; tgt++; } \
                    int iori,inew; \
                    for (iori=ifrom; iori<line->n_allele; iori++) \
                    { \
                        inew = ma->buf[i].rec[irec].map[iori] - ifrom; \
                        tgt = (tgt_type_t *) ma->tmp_arr + (ismpl+j)*nsize + inew; \
                        if ( src_is_vector_end ) break; \
                        if ( src_is_missing ) tgt_set_missing; \
                        else *tgt = *src; \
                        src++; \
                    } \
                } \
            } \
            ismpl += bcf_hdr_nsamples(hdr); \
        }
        switch (type)
        {
            case BCF_BT_INT8:  BRANCH(int32_t,  int8_t, *src==bcf_int8_missing,  *src==bcf_int8_vector_end,  *tgt=bcf_int32_missing, *tgt=bcf_int32_vector_end); break;
            case BCF_BT_INT16: BRANCH(int32_t, int16_t, *src==bcf_int16_missing, *src==bcf_int16_vector_end, *tgt=bcf_int32_missing, *tgt=bcf_int32_vector_end); break;
            case BCF_BT_INT32: BRANCH(int32_t, int32_t, *src==bcf_int32_missing, *src==bcf_int32_vector_end, *tgt=bcf_int32_missing, *tgt=bcf_int32_vector_end); break;
            case BCF_BT_FLOAT: BRANCH(float, float, bcf_float_is_missing(*src), bcf_float_is_vector_end(*src), bcf_float_set_missing(*tgt), bcf_float_set_vector_end(*tgt)); break;
            case BCF_BT_CHAR:  BRANCH(uint8_t, uint8_t, *src==bcf_str_missing, *src==bcf_str_vector_end, *tgt=bcf_str_missing, *tgt=bcf_str_vector_end); break;
            default: error("Unexpected case: %d, %s\n", type, key);
        }
        #undef BRANCH
    }
    if ( type==BCF_BT_FLOAT )
        bcf_update_format_float(out_hdr, out, key, (float*)ma->tmp_arr, nsamples*nsize);
    else if ( type==BCF_BT_CHAR )
        bcf_update_format_char(out_hdr, out, key, (float*)ma->tmp_arr, nsamples*nsize);
    else
        bcf_update_format_int32(out_hdr, out, key, (int32_t*)ma->tmp_arr, nsamples*nsize);
}

void merge_format(args_t *args, bcf1_t *out)
{
    bcf_srs_t *files = args->files;
    bcf_hdr_t *out_hdr = args->out_hdr;
    maux_t *ma = args->maux;
    if ( !ma->nfmt_map )
    {
        ma->nfmt_map = 2;
        ma->fmt_map  = (bcf_fmt_t**) calloc(ma->nfmt_map*files->nreaders, sizeof(bcf_fmt_t*));
    }
    else
        memset(ma->fmt_map, 0, ma->nfmt_map*files->nreaders*sizeof(bcf_fmt_t**));

    khiter_t kitr;
    strdict_t *tmph = args->tmph;
    kh_clear(strdict, tmph);
    int i, j, ret, has_GT = 0, max_ifmt = 0; // max fmt index
    for (i=0; i<files->nreaders; i++)
    {
        bcf1_t *line = maux_get_line(args,i);
        if ( !line ) continue;
        bcf_sr_t *reader = &files->readers[i];
        bcf_hdr_t *hdr = reader->header;
        for (j=0; j<line->n_fmt; j++)
        {
            // Wat this tag already seen?
            bcf_fmt_t *fmt = &line->d.fmt[j];
            const char *key = hdr->id[BCF_DT_ID][fmt->id].key;
            kitr = kh_get(strdict, tmph, key);

            int ifmt;
            if ( kitr != kh_end(tmph) )
                ifmt = kh_value(tmph, kitr);    // seen
            else
            {
                // new FORMAT tag
                if ( key[0]=='G' && key[1]=='T' && key[2]==0 ) { has_GT = 1; ifmt = 0; }
                else
                {
                    ifmt = ++max_ifmt;
                    if ( max_ifmt >= ma->nfmt_map )
                    {
                        ma->fmt_map = (bcf_fmt_t**) realloc(ma->fmt_map, sizeof(bcf_fmt_t*)*(max_ifmt+1)*files->nreaders);
                        memset(ma->fmt_map+ma->nfmt_map*files->nreaders, 0, (max_ifmt-ma->nfmt_map+1)*files->nreaders*sizeof(bcf_fmt_t*));
                        ma->nfmt_map = max_ifmt+1;
                    }
                }
                kitr = kh_put(strdict, tmph, key, &ret);
                kh_value(tmph, kitr) = ifmt;
            }
            ma->fmt_map[ifmt*files->nreaders+i] = fmt;
        }
        // Check if the allele numbering must be changed
        int irec = ma->buf[i].cur;
        for (j=1; j<line->n_allele; j++)
            if ( ma->buf[i].rec[irec].map[j]!=j ) break;
        ma->buf[i].rec[irec].als_differ = j==line->n_allele ? 0 : 1;
    }

    out->n_sample = bcf_hdr_nsamples(out_hdr);
    if ( has_GT )
        merge_GT(args, ma->fmt_map, out);
    update_AN_AC(out_hdr, out);

    for (i=1; i<=max_ifmt; i++)
        merge_format_field(args, &ma->fmt_map[i*files->nreaders], out);
    out->d.indiv_dirty = 1;
}

void gvcf_set_alleles(args_t *args)
{
    int i,k;
    bcf_srs_t *files = args->files;
    maux_t *maux = args->maux;
    gvcf_aux_t *gaux = maux->gvcf;
    for (i=0; i<maux->nals; i++) 
    {
        free(maux->als[i]);
        maux->als[i] = NULL;
    }
    maux->nals = 0;

    for (i=0; i<files->nreaders; i++)
    {
        if ( !gaux[i].active ) continue;
        bcf1_t *line = maux_get_line(args, i);
        int irec = maux->buf[i].cur;

        hts_expand(int, line->n_allele, maux->buf[i].rec[irec].mmap, maux->buf[i].rec[irec].map);
        if ( !maux->nals )    // first record, copy the alleles to the output
        {
            maux->nals = line->n_allele;
            hts_expand0(char*, maux->nals, maux->mals, maux->als);
            hts_expand0(int, maux->nals, maux->ncnt, maux->cnt);
            for (k=0; k<maux->nals; k++)
            {
                if ( maux->als[k] ) free(maux->als[k]);
                maux->als[k] = strdup(line->d.allele[k]);
                maux->buf[i].rec[irec].map[k] = k;
            }
        }
        else
        {
            maux->als = merge_alleles(line->d.allele, line->n_allele, maux->buf[i].rec[irec].map, maux->als, &maux->nals, &maux->mals);
            if ( !maux->als )
            {
                bcf_hdr_t *hdr = bcf_sr_get_header(args->files,i);
                error("Failed to merge alleles at %s:%d\n",bcf_seqname(hdr,line),line->pos+1);
            }
        }
    }
}

/*
    Output staged gVCF blocks, end is the last position of the block. Assuming
    gaux[i].active flags are set and maux_get_line returns correct lines.
*/
void gvcf_write_block(args_t *args, int start, int end)
{
    int i;
    maux_t *maux = args->maux;
    gvcf_aux_t *gaux = maux->gvcf;
    assert(gaux);

    // Update POS
    int min = INT_MAX;
    char ref = 'N';
    for (i=0; i<args->files->nreaders; i++)
    {
        if ( !gaux[i].active ) continue;
        if ( ref=='N' && gaux[i].line->pos==start ) ref = gaux[i].line->d.allele[0][0];
        gaux[i].line->pos = start;
    }
    for (i=0; i<args->files->nreaders; i++)
    {
        if ( !gaux[i].active ) continue;
        if ( gaux[i].end < start ) 
        { 
            gaux[i].active = 0; 
            maux->buf[i].cur = -1;
            continue; 
        }
        gaux[i].line->d.allele[0][0] = ref;
        if ( min > gaux[i].end ) min = gaux[i].end;
    }
    // Check for valid gVCF blocks in this region
    if ( min==INT_MAX )
    {
    assert(0);
        maux->gvcf_min = 0;
        return;
    }

    bcf1_t *out = args->out_line;

    gvcf_set_alleles(args);

    // Merge the staged lines
    merge_chrom2qual(args, out);
    merge_filter(args, out);
    merge_info(args, out);
    merge_format(args, out);

    if ( args->gvcf_fai && out->d.allele[0][0]=='N' )
    {
        int slen  = 0;
        char *seq = faidx_fetch_seq(args->gvcf_fai,maux->chr,out->pos,out->pos,&slen);
        if (slen)
        {
            out->d.allele[0][0] = seq[0];
            free(seq);
        }
    }

    // Update END boundary
    if ( end > start )
    {
        end++;
        bcf_update_info_int32(args->out_hdr, out, "END", &end, 1);
    }
    else
        bcf_update_info_int32(args->out_hdr, out, "END", NULL, 0);
    bcf_write1(args->out_fh, args->out_hdr, out);
    bcf_clear1(out);


    // Inactivate blocks which do not extend beyond END and find new gvcf_min
    min = INT_MAX;
    for (i=0; i<args->files->nreaders; i++)
    {
        if ( !gaux[i].active ) continue;
        if ( gaux[i].end < end )
        {
            gaux[i].active = 0; 
            maux->buf[i].cur = -1;
            continue; 
        }
        // next min END position bigger than the current one
        if ( maux->gvcf_min < gaux[i].end+1 && min > gaux[i].end+1 ) min = gaux[i].end + 1;
    }
    maux->gvcf_min = min==INT_MAX ? 0 : min;
}

/*
    Flush staged gVCF blocks. Flush everything if there are no more lines
    (done=1) or if there is a new chromosome. If still on the same chromosome,
    all hanging blocks must be ended by creating new records:
        A
            1 END=10
        B
            3 END=7
        C
            3 END=5
        out
            1 END=2  A . .
            3 END=5  A B C
            6 END=7  A B .
            8 END=10 A . .
    
*/
void gvcf_flush(args_t *args, int done)
{
    int i;
    maux_t *maux = args->maux;

    if ( !maux->chr ) return;   // first time here, nothing to flush

    int flush_until = INT_MAX;
    if ( !done )
    {
        // Get current position and chromosome
        for (i=0; i<maux->n; i++)
            if ( bcf_sr_has_line(maux->files,i) ) break;
        bcf1_t *line = bcf_sr_get_line(maux->files,i);
        bcf_hdr_t *hdr = bcf_sr_get_header(maux->files,i);

        if ( !strcmp(maux->chr,bcf_seqname(hdr,line)) ) flush_until = line->pos;    // still on the same chr
    }

    // When called on a region, trim the blocks accordingly
    int start = maux->gvcf_break>=0 ? maux->gvcf_break + 1 : maux->pos;
    if ( args->regs )
    {
        int rstart = -1, rend = -1;
        if ( regidx_overlap(args->regs,maux->chr,start,flush_until,args->regs_itr) )
        {
            // In case there are multiple regions, we treat them as one
            rstart = args->regs_itr->beg;
            while ( regitr_overlap(args->regs_itr) ) rend = args->regs_itr->end;
        }
        if ( rstart > start ) start = rstart;
        if ( rend < flush_until ) flush_until = rend+1;
    }

    // output all finished blocks
    while ( maux->gvcf_min && start < flush_until )
    {
        // does the block end before the new line or is it interrupted?
        int tmp = maux->gvcf_min < flush_until ? maux->gvcf_min : flush_until;
        if ( start > tmp-1 ) break;
        gvcf_write_block(args,start,tmp-1); // gvcf_min is 1-based
        start = tmp;
    }
}

/*
    Check incoming lines for new gVCF blocks, set pointer to the current source
    buffer (gvcf or readers).  In contrast to gvcf_flush, this function can be
    called only after maux_reset as it relies on updated maux buffers.
*/
void gvcf_stage(args_t *args, int pos)
{
    maux_t *maux = args->maux;
    gvcf_aux_t *gaux = maux->gvcf;
    bcf_srs_t *files = args->files;
    int32_t *end = (int32_t*) maux->tmp_arr;
    int i, nend = maux->ntmp_arr / sizeof(int32_t);

    maux->gvcf_break = -1;
    maux->gvcf_min = INT_MAX;
    for (i=0; i<files->nreaders; i++)
    {
        if ( gaux[i].active )
        {
            // gvcf block should not overlap with another record
            if ( maux->gvcf_min > gaux[i].end+1 ) maux->gvcf_min = gaux[i].end + 1;
            maux->buf[i].beg = 0;
            maux->buf[i].end = 1;
            maux->buf[i].cur = 0;
            continue;
        }

        // Does any of the lines have END set? It is enough to check only the
        // first line, there should be no duplicate records with END in gVCF

        if ( maux->buf[i].beg==maux->buf[i].end ) continue; // no new record

        int irec = maux->buf[i].beg;
        bcf_hdr_t *hdr = bcf_sr_get_header(files, i);
        bcf1_t *line = args->files->readers[i].buffer[irec];
        int ret = bcf_get_info_int32(hdr,line,"END",&end,&nend);
        if ( ret==1 )
        {
            // END is set, this is a new gVCF block. Cache this line in gaux[i] and swap with
            // an empty record: the gaux line must be kept until we reach its END.
            gaux[i].active = 1;
            gaux[i].end = end[0] - 1;
            SWAP(bcf1_t*,args->files->readers[i].buffer[irec],gaux[i].line);
            gaux[i].line->pos = pos;

            maux->buf[i].lines = &gaux[i].line;
            maux->buf[i].beg = 0;
            maux->buf[i].end = 1;
            maux->buf[i].cur = 0;

            // Set the rid,pos of the swapped line in the buffer or else the
            // synced reader will have a problem with the next line
            //
            args->files->readers[i].buffer[irec]->rid = maux->buf[i].rid;
            args->files->readers[i].buffer[irec]->pos = maux->pos;

            // Update block offsets
            if ( maux->gvcf_min > gaux[i].end+1 ) maux->gvcf_min = gaux[i].end + 1;
        }
        else
            maux->gvcf_break = line->pos;   // must break the gvcf block 
    }
    maux->ntmp_arr = nend * sizeof(int32_t);
    maux->tmp_arr  = end;
    if ( maux->gvcf_min==INT_MAX ) maux->gvcf_min = 0;
}


void debug_buffers(FILE *fp, bcf_srs_t *files);
void debug_buffer(FILE *fp, bcf_srs_t *files, int reader);

/*
    Flush all buffered and processed records with the same coordinate.
    Note that synced reader discards buffer[0], so that needs to stay
    untouched.
*/
void clean_buffer(args_t *args)
{
    maux_t *ma = args->maux;

    int ir;
    for (ir=0; ir<ma->n; ir++)
    {
        // Invalidate pointer to reader's buffer or else gvcf_flush will attempt
        // to use the old lines via maux_get_line()
        if ( ma->gvcf && !ma->gvcf[ir].active ) ma->buf[ir].cur = -1;

        bcf_sr_t *reader = bcf_sr_get_reader(args->files,ir);
        if ( !reader->nbuffer ) continue;   // nothing to clean

        bcf1_t **buf = reader->buffer;
        if ( buf[1]->rid!=ma->buf[ir].rid || buf[1]->pos!=ma->pos ) continue;    // nothing to flush

        int a = 1, b = 2;
        while ( b<=reader->nbuffer && buf[b]->rid==ma->buf[ir].rid && buf[b]->pos==ma->pos ) b++;
        // b now points to the first line we want to preserve
        while ( b<=reader->nbuffer )
        {
            SWAP(bcf1_t*, buf[a], buf[b]);
            a++; b++;
        }
        reader->nbuffer -= b-a;
    }
}

void debug_maux(args_t *args)
{
    bcf_srs_t *files = args->files;
    maux_t *maux = args->maux;
    int j,k,l;

    fprintf(stderr,"Alleles to merge at %d, nals=%d\n", maux->pos+1,maux->nals);
    for (j=0; j<files->nreaders; j++)
    {
        bcf_sr_t *reader = &files->readers[j];
        buffer_t *buf = &maux->buf[j];
        fprintf(stderr," reader %d: ", j);
        for (k=buf->beg; k<buf->end; k++)
        {
            if ( buf->rec[k].skip & SKIP_DONE ) continue;
            bcf1_t *line = reader->buffer[k];
            fprintf(stderr,"\t");
            if ( buf->rec[k].skip ) fprintf(stderr,"[");  // this record will not be merged in this round
            for (l=0; l<line->n_allele; l++)
                fprintf(stderr,"%s%s", l==0?"":",", line->d.allele[l]);
            if ( buf->rec[k].skip ) fprintf(stderr,"]");
        }
        fprintf(stderr,"\n");
    }
    fprintf(stderr," counts: ");
    for (j=0; j<maux->nals; j++) fprintf(stderr,"%s   %dx %s", j==0?"":",",maux->cnt[j], maux->als[j]);
    fprintf(stderr,"\n\n");
}

void debug_state(args_t *args)
{
    maux_t *maux = args->maux;
    int i,j;
    for (i=0; i<args->files->nreaders; i++)
    {
        fprintf(stderr,"reader %d:\tcur,beg,end=% d,%d,%d", i,maux->buf[i].cur,maux->buf[i].beg,maux->buf[i].end);
        if ( maux->buf[i].cur >=0 )
        {
            bcf_hdr_t *hdr = bcf_sr_get_header(args->files,i);
            const char *chr = bcf_hdr_id2name(hdr, maux->buf[i].rid);
            fprintf(stderr,"\t");
            for (j=maux->buf[i].beg; j<maux->buf[i].end; j++) fprintf(stderr," %s:%d",chr,maux->buf[i].lines[j]->pos+1);
        }
        fprintf(stderr,"\n");
    }
    for (i=0; i<args->files->nreaders; i++)
    {
        fprintf(stderr,"reader %d:\tgvcf_active=%d", i,maux->gvcf[i].active);
        if ( maux->gvcf[i].active ) fprintf(stderr,"\tpos,end=%d,%d", maux->gvcf[i].line->pos+1,maux->gvcf[i].end+1);
        fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
}


/*
   Determine which line should be merged from which reader: go through all
   readers and all buffered lines, expand REF,ALT and try to match lines with
   the same ALTs.
 */
int can_merge(args_t *args)
{
    bcf_srs_t *files = args->files;
    int snp_mask = (VCF_SNP<<1)|(VCF_MNP<<1), indel_mask = VCF_INDEL<<1, ref_mask = 1;
    maux_t *maux = args->maux;
    gvcf_aux_t *gaux = maux->gvcf;
    char *id = NULL, ref = 'N';
    int i,j,k, ntodo = 0;

    for (i=0; i<maux->nals; i++) 
    {
        free(maux->als[i]);
        maux->als[i] = NULL;
    }
    maux->var_types = maux->nals = 0;

    for (i=0; i<files->nreaders; i++)
    {
        buffer_t *buf = &maux->buf[i];

        if ( gaux && gaux[i].active )
        {
            // skip readers with active gvcf blocks
            buf->rec[buf->beg].skip = SKIP_DIFF;
            continue;
        }
        for (j=buf->beg; j<buf->end; j++)
        {
            if ( buf->rec[j].skip & SKIP_DONE ) continue;

            buf->rec[j].skip = SKIP_DIFF;
            ntodo++;

            if ( args->merge_by_id )
                id = buf->lines[j]->d.id;
            else
            {
                int var_type = bcf_get_variant_types(buf->lines[j]);
                maux->var_types |= var_type ? var_type<<1 : 1;
            }
        }

        // for gvcf: find out REF at this position
        if ( buf->beg < buf->end && ref=='N' )
            ref = buf->lines[buf->beg]->d.allele[0][0];
    }
    if ( !ntodo ) return 0;

    // In this loop we select from each reader compatible candidate lines.
    // (i.e. SNPs or indels). Go through all files and all lines at this
    // position and normalize relevant alleles.
    // REF-only sites may be associated with both SNPs and indels.
    for (i=0; i<files->nreaders; i++)
    {
        bcf_sr_t *reader = &files->readers[i];
        buffer_t *buf = &maux->buf[i];

        if ( gaux && gaux[i].active )
        {
            gaux[i].line->d.allele[0][0] = ref;
            gaux[i].line->pos = maux->pos;
        }

        for (j=buf->beg; j<buf->end; j++)
        {
            if ( buf->rec[j].skip & SKIP_DONE ) continue;

            bcf1_t *line = buf->lines[j]; // ptr to reader's buffer or gvcf buffer

            int line_type = bcf_get_variant_types(line);
            line_type = line_type ? line_type<<1 : 1;

            // select relevant lines
            if ( args->merge_by_id )
            {
                if ( strcmp(id,line->d.id) ) continue;
            }
            else
            {
                if ( args->collapse==COLLAPSE_NONE && maux->nals )
                {
                    // All alleles of the tested record must be present in the
                    // selected maux record plus variant types must be the same
                    if ( (maux->var_types & line_type) != line_type ) continue;
                    if ( vcmp_set_ref(args->vcmp,maux->als[0],line->d.allele[0]) < 0 ) continue;   // refs not compatible
                    for (k=1; k<line->n_allele; k++)
                    {
                        if ( vcmp_find_allele(args->vcmp,maux->als+1,maux->nals-1,line->d.allele[k])>=0 ) break;
                    }
                    if ( !(line_type&ref_mask) && k==line->n_allele ) continue;  // not a REF-only site and there is no matching allele
                }
                if ( !(args->collapse&COLLAPSE_ANY) )
                {
                    // Merge:
                    //  - SNPs+SNPs+MNPs+REF if -m both,snps
                    //  - indels+indels+REF  if -m both,indels, REF only if SNPs are not present
                    //  - SNPs come first
                    if ( line_type & indel_mask )
                    {
                        if ( !(line_type&snp_mask) && maux->var_types&snp_mask ) continue;  // SNPs come first
                        if ( args->do_gvcf && maux->var_types&ref_mask ) continue;  // never merge indels with gVCF blocks
                    }
                }
            }
            buf->rec[j].skip = 0;

            hts_expand(int, line->n_allele, buf->rec[j].mmap, buf->rec[j].map);
            if ( !maux->nals )    // first record, copy the alleles to the output
            {
                maux->nals = line->n_allele;
                hts_expand0(char*, maux->nals, maux->mals, maux->als);
                hts_expand0(int, maux->nals, maux->ncnt, maux->cnt);
                for (k=0; k<maux->nals; k++)
                {
                    free(maux->als[k]);
                    maux->als[k] = strdup(line->d.allele[k]);
                    buf->rec[j].map[k] = k;
                    maux->cnt[k] = 1;
                }
                continue;
            }
            // normalize alleles
            maux->als = merge_alleles(line->d.allele, line->n_allele, buf->rec[j].map, maux->als, &maux->nals, &maux->mals);
            if ( !maux->als ) error("Failed to merge alleles at %s:%d in %s\n",bcf_seqname(args->out_hdr,line),line->pos+1,reader->fname);
            hts_expand0(int, maux->nals, maux->ncnt, maux->cnt);
            for (k=1; k<line->n_allele; k++)
                maux->cnt[ buf->rec[j].map[k] ]++;    // how many times an allele appears in the files
            maux->cnt[0]++;
        }
    }
    return 1;
}

/*
   Select records that have the same alleles; the input ordering of indels
   must not matter. Multiple VCF lines can be emitted from this loop.
   We expect only very few alleles and not many records with the same
   position in the buffers, therefore the nested loops should not slow us
   much.
*/
void stage_line(args_t *args)
{
    int snp_mask = (VCF_SNP<<1)|(VCF_MNP<<1), indel_mask = VCF_INDEL<<1, ref_mask = 1;
    bcf_srs_t *files = args->files;
    maux_t *maux = args->maux;

    // debug_maux(args);

    // take the most frequent allele present in multiple files, REF is skipped
    int i,j,k,icnt = 1;
    for (i=2; i<maux->nals; i++)
        if ( maux->cnt[i] > maux->cnt[icnt] ) icnt = i;

    int nout = 0;
    for (i=0; i<files->nreaders; i++)
    {
        buffer_t *buf = &maux->buf[i];
        buf->cur = -1;
        if ( buf->beg >= buf->end ) continue;   // no lines in the buffer

        // find lines with the same allele
        for (j=buf->beg; j<buf->end; j++)
        {
            if ( buf->rec[j].skip ) continue;   // done or not compatible
            if ( args->merge_by_id ) break;
            if ( maux->nals==1 && buf->lines[j]->n_allele==1 ) break;   // REF-only record

            for (k=0; k<buf->lines[j]->n_allele; k++)
                if ( icnt==buf->rec[j].map[k] ) break;

            if ( k<buf->lines[j]->n_allele ) break;
        }
        if ( j>=buf->end )
        {
            // no matching allele found in this file
            if ( args->collapse==COLLAPSE_NONE ) continue;

            for (j=buf->beg; j<buf->end; j++)
            {
                if ( buf->rec[j].skip ) continue;   // done or not compatible
                if ( args->collapse&COLLAPSE_ANY ) break;   // anything can be merged
                int line_type = bcf_get_variant_types(buf->lines[j]);
                if ( maux->var_types&snp_mask && line_type&VCF_SNP && (args->collapse&COLLAPSE_SNPS) ) break;
                if ( maux->var_types&indel_mask && line_type&VCF_INDEL && (args->collapse&COLLAPSE_INDELS) ) break;
                if ( line_type==VCF_REF )
                {
                    if ( maux->var_types&snp_mask && (args->collapse&COLLAPSE_SNPS) ) break;
                    if ( maux->var_types&indel_mask && (args->collapse&COLLAPSE_INDELS) ) break;
                    if ( maux->var_types&ref_mask ) break;
                }
                else if ( maux->var_types&ref_mask )
                {
                    if ( line_type&snp_mask && (args->collapse&COLLAPSE_SNPS) ) break;
                    if ( line_type&indel_mask && (args->collapse&COLLAPSE_INDELS) ) break;
                }
            }
        }
        if ( j<buf->end )
        {
            // found a suitable line for merging
            buf->cur = j;

            // mark as finished so that it's ignored next time
            buf->rec[j].skip  = SKIP_DONE;
            nout++;
        }
    }
    assert( nout );
}

void merge_line(args_t *args)
{
    if ( args->regs )
    {
        if ( !regidx_overlap(args->regs,args->maux->chr,args->maux->pos,args->maux->pos,NULL) ) return;
    }

    bcf1_t *out = args->out_line;
    merge_chrom2qual(args, out);
    merge_filter(args, out);
    merge_info(args, out);
    if ( args->do_gvcf )
        bcf_update_info_int32(args->out_hdr, out, "END", NULL, 0);
    merge_format(args, out);
    bcf_write1(args->out_fh, args->out_hdr, out);
    bcf_clear1(out);
}

void bcf_hdr_append_version(bcf_hdr_t *hdr, int argc, char **argv, const char *cmd)
{
    kstring_t str = {0,0,0};
    ksprintf(&str,"##%sVersion=%s+htslib-%s\n", cmd, bcftools_version(), hts_version());
    bcf_hdr_append(hdr,str.s);

    str.l = 0;
    ksprintf(&str,"##%sCommand=%s", cmd, argv[0]);
    int i;
    for (i=1; i<argc; i++)
    {
        if ( strchr(argv[i],' ') )
            ksprintf(&str, " '%s'", argv[i]);
        else
            ksprintf(&str, " %s", argv[i]);
    }
    kputs("; Date=", &str);
    time_t tm; time(&tm); kputs(ctime(&tm), &str);
    kputc('\n', &str);
    bcf_hdr_append(hdr,str.s);
    free(str.s);

    bcf_hdr_sync(hdr);
}

void merge_vcf(args_t *args)
{
    args->out_fh  = hts_open(args->output_fname, hts_bcf_wmode(args->output_type));
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    if ( args->n_threads ) hts_set_opt(args->out_fh, HTS_OPT_THREAD_POOL, args->files->p); //hts_set_threads(args->out_fh, args->n_threads);
    args->out_hdr = bcf_hdr_init("w");

    if ( args->header_fname )
    {
        if ( bcf_hdr_set(args->out_hdr,args->header_fname) ) error("Could not read/parse the header: %s\n", args->header_fname);
    }
    else
    {
        int i;
        for (i=0; i<args->files->nreaders; i++)
        {
            char buf[10]; snprintf(buf,10,"%d",i+1);
            merge_headers(args->out_hdr, args->files->readers[i].header,buf,args->force_samples);
        }
        if (args->record_cmd_line) bcf_hdr_append_version(args->out_hdr, args->argc, args->argv, "bcftools_merge");
        bcf_hdr_sync(args->out_hdr);
    }
    info_rules_init(args);

    bcf_hdr_set_version(args->out_hdr, bcf_hdr_get_version(args->files->readers[0].header));
    bcf_hdr_write(args->out_fh, args->out_hdr);
    if ( args->header_only )
    {
        bcf_hdr_destroy(args->out_hdr);
        hts_close(args->out_fh);
        return;
    }

    if ( args->collapse==COLLAPSE_NONE ) args->vcmp = vcmp_init();
    args->maux = maux_init(args);
    args->out_line = bcf_init1();
    args->tmph = kh_init(strdict);

    while ( bcf_sr_next_line(args->files) )
    {
        // output cached gVCF blocks which end before the new record
        if ( args->do_gvcf )
            gvcf_flush(args,0);

        maux_reset(args->maux);

        // determine which of the new records are gvcf blocks
        if ( args->do_gvcf )
            gvcf_stage(args, args->maux->pos);

        while ( can_merge(args) )
        {
            stage_line(args);
            merge_line(args);
        }
        clean_buffer(args);
        // debug_state(args);
    }
    if ( args->do_gvcf )
        gvcf_flush(args,1);

    info_rules_destroy(args);
    maux_destroy(args->maux);
    bcf_hdr_destroy(args->out_hdr);
    hts_close(args->out_fh);
    bcf_destroy1(args->out_line);
    kh_destroy(strdict, args->tmph);
    if ( args->tmps.m ) free(args->tmps.s);
    if ( args->vcmp ) vcmp_destroy(args->vcmp);
}

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file.\n");
    fprintf(stderr, "         Note that only records from different files can be merged, never from the same file. For\n");
    fprintf(stderr, "         \"vertical\" merge take a look at \"bcftools norm\" instead.\n");
    fprintf(stderr, "Usage:   bcftools merge [options] <A.vcf.gz> <B.vcf.gz> [...]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "        --force-samples                resolve duplicate sample names\n");
    fprintf(stderr, "        --print-header                 print only the merged header and exit\n");
    fprintf(stderr, "        --use-header <file>            use the provided header\n");
    fprintf(stderr, "    -0  --missing-to-ref               assume genotypes at missing sites are 0/0\n");
    fprintf(stderr, "    -f, --apply-filters <list>         require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n");
    fprintf(stderr, "    -F, --filter-logic <x|+>           remove filters if some input is PASS (\"x\"), or apply all filters (\"+\") [+]\n");
    fprintf(stderr, "    -g, --gvcf <-|ref.fa>              merge gVCF blocks, INFO/END tag is expected. Implies -i QS:sum,MinDP:min,I16:sum,IDV:max,IMF:max\n");
    fprintf(stderr, "    -i, --info-rules <tag:method,..>   rules for merging INFO fields (method is one of sum,avg,min,max,join) or \"-\" to turn off the default [DP:sum,DP4:sum]\n");
    fprintf(stderr, "    -l, --file-list <file>             read file names from the file\n");
    fprintf(stderr, "    -m, --merge <string>               allow multiallelic records for <snps|indels|both|all|none|id>, see man page for details [both]\n");
    fprintf(stderr, "        --no-version                   do not append version and command line to the header\n");
    fprintf(stderr, "    -o, --output <file>                write output to a file [standard output]\n");
    fprintf(stderr, "    -O, --output-type <b|u|z|v>        'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]\n");
    fprintf(stderr, "    -r, --regions <region>             restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>          restrict to regions listed in a file\n");
    fprintf(stderr, "        --threads <int>                number of extra output compression threads [0]\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfmerge(int argc, char *argv[])
{
    int c;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->files  = bcf_sr_init();
    args->argc   = argc; args->argv = argv;
    args->output_fname = "-";
    args->output_type = FT_VCF;
    args->n_threads = 0;
    args->record_cmd_line = 1;
    args->collapse = COLLAPSE_BOTH;
    int regions_is_file = 0;

    static struct option loptions[] =
    {
        {"help",no_argument,NULL,'h'},
        {"merge",required_argument,NULL,'m'},
        {"gvcf",required_argument,NULL,'g'},
        {"file-list",required_argument,NULL,'l'},
        {"missing-to-ref",no_argument,NULL,'0'},
        {"apply-filters",required_argument,NULL,'f'},
        {"use-header",required_argument,NULL,1},
        {"print-header",no_argument,NULL,2},
        {"force-samples",no_argument,NULL,3},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"threads",required_argument,NULL,9},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"info-rules",required_argument,NULL,'i'},
        {"no-version",no_argument,NULL,8},
        {"filter-logic",required_argument,NULL,'F'},
        {NULL,0,NULL,0}
    };
    while ((c = getopt_long(argc, argv, "hm:f:r:R:o:O:i:l:g:F:0",loptions,NULL)) >= 0) {
        switch (c) {
            case 'F': 
                if ( !strcmp(optarg,"+") ) args->filter_logic = FLT_LOGIC_ADD;
                else if ( !strcmp(optarg,"x") ) args->filter_logic = FLT_LOGIC_REMOVE;
                else error("Filter logic not recognised: %s\n", optarg);
                break;
            case '0': args->missing_to_ref = 1; break;
            case 'g':
                args->do_gvcf = 1;
                if ( strcmp("-",optarg) )
                {
                    args->gvcf_fai = fai_load(optarg);
                    if ( !args->gvcf_fai ) error("Failed to load the fai index: %s\n", optarg);
                }
                break;
            case 'l': args->file_list = optarg; break;
            case 'i': args->info_rules = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                }
                break;
            case 'm':
                args->collapse = COLLAPSE_NONE;
                if ( !strcmp(optarg,"snps") ) args->collapse |= COLLAPSE_SNPS;
                else if ( !strcmp(optarg,"indels") ) args->collapse |= COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"both") ) args->collapse |= COLLAPSE_BOTH;
                else if ( !strcmp(optarg,"any") ) args->collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"all") ) args->collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"none") ) args->collapse = COLLAPSE_NONE;
                else if ( !strcmp(optarg,"id") ) { args->collapse = COLLAPSE_NONE; args->merge_by_id = 1; }
                else error("The -m type \"%s\" is not recognised.\n", optarg);
                break;
            case 'f': args->files->apply_filters = optarg; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case  1 : args->header_fname = optarg; break;
            case  2 : args->header_only = 1; break;
            case  3 : args->force_samples = 1; break;
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case  8 : args->record_cmd_line = 0; break;
            case 'h':
            case '?': usage();
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( argc==optind && !args->file_list ) usage();
    if ( argc-optind<2 && !args->file_list ) usage();

    args->files->require_index = 1;
    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
        if ( regions_is_file )
            args->regs = regidx_init(args->regions_list,NULL,NULL,sizeof(char*),NULL);
        else
        {
            args->regs = regidx_init(NULL,regidx_parse_reg,NULL,sizeof(char*),NULL);
            if ( regidx_insert_list(args->regs,args->regions_list,',') !=0 ) error("Could not parse the regions: %s\n", args->regions_list);
            regidx_insert(args->regs,NULL);
        }
        if ( !args->regs ) error("Could not parse the regions: %s\n", args->regions_list);
        args->regs_itr = regitr_init(args->regs);
    }

    if ( bcf_sr_set_threads(args->files, args->n_threads)<0 ) error("Failed to create threads\n");
    while (optind<argc)
    {
        if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open %s: %s\n", argv[optind],bcf_sr_strerror(args->files->errnum));
        optind++;
    }
    if ( args->file_list )
    {
        int nfiles, i;
        char **files = hts_readlines(args->file_list, &nfiles);
        if ( !files ) error("Failed to read from %s\n", args->file_list);
        for (i=0;i<nfiles; i++)
            if ( !bcf_sr_add_reader(args->files, files[i]) ) error("Failed to open %s: %s\n", files[i],bcf_sr_strerror(args->files->errnum));
        for (i=0; i<nfiles; i++) free(files[i]);
        free(files);
    }
    merge_vcf(args);
    bcf_sr_destroy(args->files);
    if ( args->regs ) regidx_destroy(args->regs);
    if ( args->regs_itr ) regitr_destroy(args->regs_itr);
    if ( args->gvcf_fai ) fai_destroy(args->gvcf_fai);
    free(args);
    return 0;
}

