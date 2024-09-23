/*  filter.c -- filter expressions.

    Copyright (C) 2013-2024 Genome Research Ltd.

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

#include <ctype.h>
#include <stdlib.h>
#include <strings.h>
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <sys/types.h>
#include <stdint.h>
#include <inttypes.h>
#ifndef _WIN32
#include <pwd.h>
#endif
#include <regex.h>
#include <htslib/khash_str2int.h>
#include <htslib/hts_defs.h>
#include <htslib/vcfutils.h>
#include <htslib/kfunc.h>
#include <htslib/hts_endian.h>
#include "config.h"
#include "filter.h"
#include "bcftools.h"

#if ENABLE_PERL_FILTERS
// Work around clang warning problems
#  if defined(__clang__)
#    define PERL_GCC_BRACE_GROUPS_FORBIDDEN
#  endif

#  define filter_t perl_filter_t
#  include <EXTERN.h>
#  include <perl.h>
#  undef filter_t
#  define my_perl perl

static int filter_ninit = 0;
#endif


#ifndef __FUNCTION__
#  define __FUNCTION__ __func__
#endif

typedef struct _token_t
{
    // read-only values, same for all VCF lines
    int tok_type;       // one of the TOK_* keys below
    int nargs;          // with TOK_PERLSUB the first argument is the name of the subroutine
    char *key;          // set only for string constants, otherwise NULL
    char *tag;          // for debugging and printout only, VCF tag name
    double threshold;   // filtering threshold
    int is_constant;    // the threshold is set
    int hdr_id, hl_type, ht_type;   // BCF header lookup ID and one of BCF_HL_* types and BCF_HT_* types
    int idx;            // 0-based index to VCF vectors,
                        //  -2: list (e.g. [0,1,2] or [1..3] or [1..] or any field[*], which is equivalent to [0..])
                        //  -3: select indices on the fly based on values in GT
    int *idxs;          // set indexes to 0 to exclude, to 1 to include, and last element negative if unlimited; used by VCF retrievers only
    int nidxs, nuidxs;  // size of idxs array and the number of elements set to 1
    uint8_t *usmpl;     // bitmask of used samples as set by idx, set for FORMAT fields, NULL otherwise
    int nsamples;       // number of samples for format fields, 0 for info and other fields
    void (*setter)(filter_t *, bcf1_t *, struct _token_t *);
    int (*func)(filter_t *, bcf1_t *, struct _token_t *rtok, struct _token_t **stack, int nstack);
    void (*comparator)(struct _token_t *, struct _token_t *, struct _token_t *rtok, bcf1_t *);
    void *hash;         // test presence of str value in the hash via comparator
    regex_t *regex;     // precompiled regex for string comparison
    int iext;           // for the use with filter_test_ext(), 1-based index to external values, 0=don't use

    // modified on filter evaluation at each VCF line
    double *values;
    kstring_t str_value;
    int is_str, is_missing; // is_missing is set only for constants, variables are controlled via nvalues
    int pass_site;          // -1 not applicable, 0 fails, >0 pass
    uint8_t *pass_samples;  // status of individual samples
    int nvalues, mvalues;   // number of used values: n=0 for missing values, n=1 for scalars, for strings n=str_value.l
    int nval1;              // number of per-sample fields or string length
}
token_t;

struct _filter_t
{
    bcf_hdr_t *hdr;
    char *str;
    int nfilters;
    token_t *filters, **flt_stack;  // filtering input tokens (in RPN) and evaluation stack
    int32_t *tmpi;
    float   *tmpf;
    kstring_t tmps;
    int max_unpack, mtmpi, mtmpf, nsamples;
    struct {
        bcf1_t *line;
        int32_t *buf, nbuf, mbuf;   // GTs as obtained by bcf_get_genotypes()
        uint64_t *mask;             // GTs as mask, e.g 0/0 is 1; 0/1 is 3, max 63 unique alleles
    } cached_GT;
#if ENABLE_PERL_FILTERS
    PerlInterpreter *perl;
#endif
    char **undef_tag, **used_tag;
    int nundef_tag, nused_tag;
    int status, exit_on_error;
    int n_ext;      // number of external values to fill via filter_test_ext()
    int *ext;       // types of external values to fill via filter_test_ext()
};


#define TOK_VAL     0
#define TOK_LFT     1       // (
#define TOK_RGT     2       // )
#define TOK_LE      3       // less or equal
#define TOK_LT      4       // less than
#define TOK_EQ      5       // equal
#define TOK_BT      6       // bigger than
#define TOK_BE      7       // bigger or equal
#define TOK_NE      8       // not equal
#define TOK_OR      9       // |
#define TOK_AND     10      // &
#define TOK_ADD     11      // +
#define TOK_SUB     12      // -
#define TOK_MULT    13      // *
#define TOK_DIV     14      // /
#define TOK_MAX     15
#define TOK_MIN     16
#define TOK_AVG     17
#define TOK_AND_VEC 18      // &&   (operator applied in samples)
#define TOK_OR_VEC  19      // ||   (operator applied in samples)
#define TOK_LIKE    20      //  ~ regular expression
#define TOK_NLIKE   21      // !~ regular expression
#define TOK_SUM     22
#define TOK_ABS     23
#define TOK_LEN     24
#define TOK_FUNC    25
#define TOK_CNT     26
#define TOK_PERLSUB 27
#define TOK_BINOM   28
#define TOK_PHRED   29
#define TOK_MEDIAN  30
#define TOK_STDEV   31
#define TOK_sMAX    32
#define TOK_sMIN    33
#define TOK_sAVG    34
#define TOK_sMEDIAN 35
#define TOK_sSTDEV  36
#define TOK_sSUM    37
#define TOK_IN      38      // contains, e.g. FILTER~"A"
#define TOK_NOT_IN  39      // does not contain, e.g. FILTER!~"A"
#define TOK_MODULO  40      // %
#define TOK_EXT     41      // external values set before each filter_test_ext() call, can be one of {},{str},{int},{float}

//                      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41
//                        ( ) [ < = > ] ! | &  +  -  *  /  M  m  a  A  O  ~  ^  S  .  l  f  c  p  b  P  i  s                          %
static int op_prec[] = {0,1,1,5,5,5,5,5,5,2,3, 6, 6, 7, 7, 8, 8, 8, 3, 2, 5, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 0};
#define TOKEN_STRING "x()[<=>]!|&+-*/MmaAO~^S.lfcpis"       // this is only for debugging, not maintained diligently

static void cmp_vector_strings(token_t *atok, token_t *btok, token_t *rtok);
inline static void tok_init_samples(token_t *atok, token_t *btok, token_t *rtok);


// Return negative values if it is a function with variable number of arguments
static int filters_next_token(char **str, int *len)
{
    char *tmp = *str;
    while ( *tmp && isspace(*tmp) ) tmp++;
    *str = tmp;
    *len = 0;

    // test for doubles: d.ddde[+-]dd
    if ( isdigit(*str[0]) || *str[0]=='.' )   // strtod would eat +/-
    {
        double HTS_UNUSED v = strtod(*str, &tmp);
        if ( *str!=tmp && (!tmp[0] || !isalnum(tmp[0])) )
        {
            *len = tmp - (*str);
            return TOK_VAL;
        }
        tmp = *str;
    }

    if ( !strncasecmp(tmp,"SMPL_MAX(",9) ) { (*str) += 8; return TOK_sMAX; }
    if ( !strncasecmp(tmp,"SMPL_MIN(",9) ) { (*str) += 8; return TOK_sMIN; }
    if ( !strncasecmp(tmp,"SMPL_MEAN(",10) ) { (*str) += 9; return TOK_sAVG; }
    if ( !strncasecmp(tmp,"SMPL_MEDIAN(",12) ) { (*str) += 11; return TOK_sMEDIAN; }
    if ( !strncasecmp(tmp,"SMPL_AVG(",9) ) { (*str) += 8; return TOK_sAVG; }
    if ( !strncasecmp(tmp,"SMPL_STDEV(",11) ) { (*str) += 10; return TOK_sSTDEV; }
    if ( !strncasecmp(tmp,"SMPL_SUM(",9) ) { (*str) += 8; return TOK_sSUM; }
    if ( !strncasecmp(tmp,"sMAX(",5) ) { (*str) += 4; return TOK_sMAX; }
    if ( !strncasecmp(tmp,"sMIN(",5) ) { (*str) += 4; return TOK_sMIN; }
    if ( !strncasecmp(tmp,"sMEAN(",6) ) { (*str) += 5; return TOK_sAVG; }
    if ( !strncasecmp(tmp,"sMEDIAN(",8) ) { (*str) += 7; return TOK_sMEDIAN; }
    if ( !strncasecmp(tmp,"sAVG(",5) ) { (*str) += 4; return TOK_sAVG; }
    if ( !strncasecmp(tmp,"sSTDEV(",7) ) { (*str) += 6; return TOK_sSTDEV; }
    if ( !strncasecmp(tmp,"sSUM(",5) ) { (*str) += 4; return TOK_sSUM; }
    if ( !strncasecmp(tmp,"MAX(",4) ) { (*str) += 3; return TOK_MAX; }
    if ( !strncasecmp(tmp,"MIN(",4) ) { (*str) += 3; return TOK_MIN; }
    if ( !strncasecmp(tmp,"MEAN(",5) ) { (*str) += 4; return TOK_AVG; }
    if ( !strncasecmp(tmp,"MEDIAN(",7) ) { (*str) += 6; return TOK_MEDIAN; }
    if ( !strncasecmp(tmp,"AVG(",4) ) { (*str) += 3; return TOK_AVG; }
    if ( !strncasecmp(tmp,"STDEV(",6) ) { (*str) += 5; return TOK_STDEV; }
    if ( !strncasecmp(tmp,"SUM(",4) ) { (*str) += 3; return TOK_SUM; }
    if ( !strncasecmp(tmp,"ABS(",4) ) { (*str) += 3; return TOK_ABS; }
    if ( !strncasecmp(tmp,"COUNT(",4) ) { (*str) += 5; return TOK_CNT; }
    if ( !strncasecmp(tmp,"STRLEN(",7) ) { (*str) += 6; return TOK_LEN; }
    if ( !strncasecmp(tmp,"BINOM(",6) ) { (*str) += 5; return -TOK_BINOM; }
    if ( !strncasecmp(tmp,"PHRED(",6) ) { (*str) += 5; return TOK_PHRED; }
    if ( !strncasecmp(tmp,"%MAX(",5) ) { (*str) += 4; return TOK_MAX; } // for backward compatibility
    if ( !strncasecmp(tmp,"%MIN(",5) ) { (*str) += 4; return TOK_MIN; } // for backward compatibility
    if ( !strncasecmp(tmp,"%AVG(",5) ) { (*str) += 4; return TOK_AVG; } // for backward compatibility
    if ( !strncasecmp(tmp,"%SUM(",5) ) { (*str) += 4; return TOK_SUM; } // for backward compatibility
    if ( !strncasecmp(tmp,"INFO/",5) ) tmp += 5;
    if ( !strncasecmp(tmp,"FORMAT/",7) ) tmp += 7;
    if ( !strncasecmp(tmp,"FMT/",4) ) tmp += 4;
    if ( !strncasecmp(tmp,"PERL.",5) ) { (*str) += 5; return -TOK_PERLSUB; }
    if ( !strncasecmp(tmp,"N_PASS(",7) ) { *len = 6; (*str) += 6; return -TOK_FUNC; }
    if ( !strncasecmp(tmp,"F_PASS(",7) ) { *len = 6; (*str) += 6; return -TOK_FUNC; }
    if ( !strncasecmp(tmp,"%ILEN",5) ) { *len = 5; return TOK_VAL; }    // to be able to distinguish between INFO/ILEN and on-the-fly ILEN
    if ( !strncasecmp(tmp,"{}",2) ) { *len = 2; return TOK_EXT; }
    if ( !strncasecmp(tmp,"{STR}",5) ) { *len = 5; return TOK_EXT; }
    if ( !strncasecmp(tmp,"{INT}",5) ) { *len = 5; return TOK_EXT; }
    if ( !strncasecmp(tmp,"{FLOAT}",7) ) { *len = 7; return TOK_EXT; }

    if ( tmp[0]=='@' )  // file name
    {
        while ( *tmp && !isspace(*tmp) && *tmp!='=' && *tmp!='!' ) tmp++;
        *len = tmp - (*str);
        return TOK_VAL;
    }

    int square_brackets = 0;
    while ( tmp[0] )
    {
        if ( !square_brackets )
        {
            if ( tmp[0]=='"' ) break;
            if ( tmp[0]=='\'' ) break;
            if ( isspace(tmp[0]) ) break;
            if ( tmp[0]=='<' ) break;
            if ( tmp[0]=='>' ) break;
            if ( tmp[0]=='=' ) break;
            if ( tmp[0]=='!' ) break;
            if ( tmp[0]=='&' ) break;
            if ( tmp[0]=='|' ) break;
            if ( tmp[0]=='(' ) break;
            if ( tmp[0]==')' ) break;
            if ( tmp[0]=='+' ) break;
            if ( tmp[0]=='*' ) break;
            if ( tmp[0]=='-' ) break;
            if ( tmp[0]=='/' ) break;
            if ( tmp[0]=='~' ) break;
            if ( tmp[0]=='%' ) break;
        }
        if ( tmp[0]==']' ) { if (square_brackets) tmp++; break; }
        if ( tmp[0]=='[' ) square_brackets++;
        tmp++;
    }
    if ( tmp > *str )
    {
        *len = tmp - (*str);
        return TOK_VAL;
    }
    if ( tmp[0]=='"' || tmp[0]=='\'' )
    {
        int quote = tmp[0];
        tmp++;
        while ( *tmp && tmp[0]!=quote ) tmp++;
        if ( !*tmp ) return -1;     // missing quotes
        *len = tmp - (*str) + 1;
        return TOK_VAL;
    }
    if ( tmp[0]=='!' )
    {
        if ( tmp[1]=='=' ) { (*str) += 2; return TOK_NE; }
        if ( tmp[1]=='~' ) { (*str) += 2; return TOK_NLIKE; }
    }
    if ( tmp[0]=='<' )
    {
        if ( tmp[1]=='=' ) { (*str) += 2; return TOK_LE; }
        (*str) += 1; return TOK_LT;
    }
    if ( tmp[0]=='>' )
    {
        if ( tmp[1]=='=' ) { (*str) += 2; return TOK_BE; }
        (*str) += 1; return TOK_BT;
    }
    if ( tmp[0]=='=' )
    {
        if ( tmp[1]=='=' ) { (*str) += 2; return TOK_EQ; }
        (*str) += 1; return TOK_EQ;
    }
    if ( tmp[0]=='(' ) { (*str) += 1; return TOK_LFT; }
    if ( tmp[0]==')' ) { (*str) += 1; return TOK_RGT; }
    if ( tmp[0]=='&' && tmp[1]=='&' ) { (*str) += 2; return TOK_AND_VEC; }
    if ( tmp[0]=='|' && tmp[1]=='|' ) { (*str) += 2; return TOK_OR_VEC; }
    if ( tmp[0]=='&' ) { (*str) += 1; return TOK_AND; }
    if ( tmp[0]=='|' ) { (*str) += 1; return TOK_OR; }
    if ( tmp[0]=='+' ) { (*str) += 1; return TOK_ADD; }
    if ( tmp[0]=='-' ) { (*str) += 1; return TOK_SUB; }
    if ( tmp[0]=='*' ) { (*str) += 1; return TOK_MULT; }
    if ( tmp[0]=='/' ) { (*str) += 1; return TOK_DIV; }
    if ( tmp[0]=='~' ) { (*str) += 1; return TOK_LIKE; }
    if ( tmp[0]=='%' ) { (*str) += 1; return TOK_MODULO; }

    *len = tmp - (*str);
    return TOK_VAL;
}

#define FILTER_OK 0
#define FILTER_ERR_UNKN_TAGS 1
#define FILTER_ERR_OTHER 2

static void filter_add_undef_tag(filter_t *filter, char *str)
{
    int i;
    for (i=0; i<filter->nundef_tag; i++)
        if ( !strcmp(str,filter->undef_tag[i]) ) break;
    if ( i<filter->nundef_tag ) return;
    filter->nundef_tag++;
    filter->undef_tag = (char**)realloc(filter->undef_tag,sizeof(*filter->undef_tag)*filter->nundef_tag);
    if ( !filter->undef_tag ) error("Could not allocate memory\n");
    filter->undef_tag[filter->nundef_tag-1] = strdup(str);
    if ( !filter->undef_tag[filter->nundef_tag-1] ) error("Could not allocate memory\n");
}
const char **filter_list_undef_tags(filter_t *filter, int *ntags)
{
    *ntags = filter->nundef_tag;
    return (const char**)filter->undef_tag;
}
static void filter_add_used_tag(filter_t *filter, const char *prefix, char *str)
{
    int i;
    kstring_t tmp = {0,0,0};
    if ( prefix ) kputs(prefix,&tmp);
    kputs(str,&tmp);
    for (i=0; i<filter->nused_tag; i++)
        if ( !strcmp(tmp.s,filter->used_tag[i]) ) break;
    if ( i<filter->nused_tag )
    {
        free(tmp.s);
        return;
    }

    filter->nused_tag++;
    filter->used_tag = (char**)realloc(filter->used_tag,sizeof(*filter->used_tag)*filter->nused_tag);
    if ( !filter->used_tag ) error("Could not allocate memory\n");
    filter->used_tag[filter->nused_tag-1] = tmp.s;
    if ( !filter->used_tag[filter->nused_tag-1] ) error("Could not allocate memory\n");
}
const char **filter_list_used_tags(filter_t *filter, int *ntags)
{
    *ntags = filter->nused_tag;
    return (const char**)filter->used_tag;
}



/*
    Simple path expansion, expands ~/, ~user, $var. The result must be freed by the caller.

    Based on jkb's staden code with some adjustments.
    https://sourceforge.net/p/staden/code/HEAD/tree/staden/trunk/src/Misc/getfile.c#l123
*/
char *expand_path(char *path)
{
    kstring_t str = {0,0,0};

    if ( path[0] == '~' )
    {
        if ( !path[1] || path[1] == '/' )
        {
#ifdef _WIN32
            kputs(getenv("HOMEDRIVE"), &str);
            kputs(getenv("HOMEPATH"), &str);
#else
            // ~ or ~/path
            kputs(getenv("HOME"), &str);
            if ( path[1] ) kputs(path+1, &str);
#endif
        }
#ifndef _WIN32
        else
        {
            // user name: ~pd3/path
            char *end = path;
            while ( *end && *end!='/' ) end++;
            kputsn(path+1, end-path-1, &str);
            struct passwd *pwentry = getpwnam(str.s);
            str.l = 0;

            if ( !pwentry ) kputsn(path, end-path, &str);
            else kputs(pwentry->pw_dir, &str);
            kputs(end, &str);
        }
#endif
        return ks_release(&str);
    }
    if ( path[0] == '$' )
    {
        char *var = getenv(path+1);
        if ( var ) {
            kputs(var, &str);
            return ks_release(&str);
        }
    }

    return strdup(path);
}
static int filters_cache_genotypes(filter_t *flt, bcf1_t *line)
{
    if ( flt->cached_GT.line==line ) return flt->cached_GT.nbuf > 0 ? 0 : -1;
    flt->cached_GT.line = line;
    flt->cached_GT.nbuf = bcf_get_genotypes(flt->hdr, line, &flt->cached_GT.buf, &flt->cached_GT.mbuf);
    if ( flt->cached_GT.nbuf<=0 ) return -1;
    if ( !flt->cached_GT.mask )
    {
        flt->cached_GT.mask = (uint64_t*) malloc(sizeof(*flt->cached_GT.mask)*flt->nsamples);
        if ( !flt->cached_GT.mask ) error("Could not alloc %zu bytes\n",sizeof(*flt->cached_GT.mask)*flt->nsamples);
    }
    int i,j, ngt1 = flt->cached_GT.nbuf / line->n_sample;
    for (i=0; i<line->n_sample; i++)
    {
        int32_t *ptr = flt->cached_GT.buf + i*ngt1;
        flt->cached_GT.mask[i] = 0;
        for (j=0; j<ngt1; j++)
        {
            if ( bcf_gt_is_missing(ptr[j]) ) continue;
            if ( ptr[j]==bcf_int32_vector_end ) break;
            int allele = bcf_gt_allele(ptr[j]);
            if ( allele > 63 )
            {
                static int warned = 0;
                if ( !warned )
                {
                    fprintf(stderr,"Too many alleles, skipping GT filtering at this site %s:%"PRId64". "
                                   "(This warning is printed only once.)\n", bcf_seqname(flt->hdr,line),line->pos+1);
                    warned = 1;
                }
                flt->cached_GT.nbuf = 0;
                return -1;
            }
            flt->cached_GT.mask[i] |= 1<<allele;
        }
    }
    return 0;
}

static void filters_set_qual(filter_t *flt, bcf1_t *line, token_t *tok)
{
    float *ptr = &line->qual;
    if ( bcf_float_is_missing(*ptr) )
        bcf_double_set_missing(tok->values[0]);
    else
        tok->values[0] = (double)line->qual;
    tok->nvalues  = 1;
}
static void filters_set_type(filter_t *flt, bcf1_t *line, token_t *tok)
{
    tok->values[0] = bcf_get_variant_types(line);
    if ( !tok->values[0] ) tok->values[0] = 1;      // mistake in htslib: VCF_* should start with 1
    else tok->values[0] = ((int)tok->values[0]) << 1;
    tok->nvalues = 1;
}
static void filters_set_info(filter_t *flt, bcf1_t *line, token_t *tok)
{
    assert( tok->hdr_id >=0  );
    int i;
    for (i=0; i<line->n_info; i++)
        if ( line->d.info[i].key == tok->hdr_id ) break;

    if ( i==line->n_info )
        tok->nvalues = tok->str_value.l = 0;
    else if ( line->d.info[i].type==BCF_BT_CHAR )
    {
        int n = line->d.info[i].len;
        if ( n >= tok->str_value.m )
        {
            tok->str_value.m = n + 1;
            tok->str_value.s = (char*) realloc(tok->str_value.s, tok->str_value.m);
            if ( !tok->str_value.s ) error("Failed to alloc %d bytes\n", (int)tok->str_value.m);
        }
        memcpy(tok->str_value.s, line->d.info[i].vptr, n);
        tok->str_value.s[n] = 0;
        tok->nvalues = tok->str_value.l = n;
    }
    else if ( line->d.info[i].type==BCF_BT_FLOAT )
    {
        if ( bcf_float_is_missing(line->d.info[i].v1.f) ) tok->nvalues = 0;
        else
        {
            tok->values[0] = line->d.info[i].v1.f;
            tok->nvalues   = 1;
        }
        tok->str_value.l = 0;
    }
    else
    {
        tok->str_value.l = 0;
        if ( line->d.info[i].type==BCF_BT_INT8 && line->d.info[i].v1.i==bcf_int8_missing ) tok->nvalues = 0;
        else if ( line->d.info[i].type==BCF_BT_INT16 && line->d.info[i].v1.i==bcf_int16_missing ) tok->nvalues = 0;
        else if ( line->d.info[i].type==BCF_BT_INT32 && line->d.info[i].v1.i==bcf_int32_missing ) tok->nvalues = 0;
        else
        {
            tok->values[0] = line->d.info[i].v1.i;
            tok->nvalues   = 1;
        }
    }
}
static void filters_cmp_bit_and(token_t *atok, token_t *btok, token_t *rtok, bcf1_t *line)
{
    int a = (int)(atok->nvalues?atok->values[0]:atok->threshold);
    int b = (int)(btok->nvalues?btok->values[0]:btok->threshold);
    if ( rtok->tok_type==TOK_LIKE )
        rtok->pass_site = a&b ? 1 : 0;
    else
        rtok->pass_site = a&b ? 0 : 1;
}
static void filters_cmp_filter(token_t *atok, token_t *btok, token_t *rtok, bcf1_t *line)
{
    // the btok values contain FILTER ids obtained by parsing the user expression
    int i,j;
    if ( rtok->tok_type==TOK_NOT_IN )   // fail if the query expression is a subset of the VCF FILTER
    {
        if ( !btok->nvalues )   // the query expression is ".", pass everything unless the VCF is also "."
        {
            if ( line->d.n_flt ) rtok->pass_site = 1;
            return;
        }
        if ( !line->d.n_flt )   // no filters at this VCF line and the query expression has a value
        {
            rtok->pass_site = 1;
            return;
        }
        for (j=0; j<btok->nvalues; j++) // some query expression value must be absent from VCF in order to pass
        {
            for (i=0; i<line->d.n_flt; i++)
                if ( btok->values[j]==line->d.flt[i] ) break;
            if ( i==line->d.n_flt ) break;  // the query is not in the VCF
        }
        if ( j!=btok->nvalues ) rtok->pass_site = 1;
        return;
    }
    else if ( rtok->tok_type==TOK_IN )
    {
        if ( !btok->nvalues )   // the query expression is ".", fail everything unless the VCF is also "."
        {
            if ( !line->d.n_flt ) rtok->pass_site = 1;
            return;
        }
        if ( !line->d.n_flt ) return;   // no filters at this VCF line and the query expression has a value
        for (j=0; j<btok->nvalues; j++) // all of the query values must be present in the VCF in order to pass
        {
            for (i=0; i<line->d.n_flt; i++)
                if ( btok->values[j]==line->d.flt[i] ) break;
            if ( i==line->d.n_flt ) break;  // the query is not in the VCF
        }
        if ( j==btok->nvalues ) rtok->pass_site = 1;
        return;
    }
    else if ( rtok->tok_type==TOK_NE )  // require anything but exact match
    {
        if ( btok->nvalues != line->d.n_flt )
        {
            rtok->pass_site = 1;
            return;
        }
        if ( !btok->nvalues ) return;
        for (j=0; j<btok->nvalues; j++) // some of the query values must be absent from the VCF in order to pass
        {
            for (i=0; i<line->d.n_flt; i++)
                if ( btok->values[j]==line->d.flt[i] ) break;
            if ( i==line->d.n_flt ) break;  // the query is not in the VCF
        }
        if ( j!=btok->nvalues ) rtok->pass_site = 1;
        return;
    }
    else if ( rtok->tok_type==TOK_EQ )  // require exact match
    {
        if ( btok->nvalues != line->d.n_flt ) return;
        if ( !btok->nvalues )
        {
            rtok->pass_site = 1;
            return;
        }
        for (j=0; j<btok->nvalues; j++) // all of the query values must be present in the VCF in order to pass
        {
            for (i=0; i<line->d.n_flt; i++)
                if ( btok->values[j]==line->d.flt[i] ) break;
            if ( i==line->d.n_flt ) break;  // the query is not in the VCF
        }
        if ( j==btok->nvalues ) rtok->pass_site = 1;
        return;
    }
    else
        error("Only ==, !=, ~, and !~ operators are supported for FILTER\n");
    return;
}
static void filters_cmp_string_hash(token_t *atok, token_t *btok, token_t *rtok, bcf1_t *line)
{
    if ( btok->hash )
    {
        token_t *tmp = atok; atok = btok; btok = tmp;
    }
    if ( rtok->tok_type!=TOK_EQ && rtok->tok_type!=TOK_NE )
        error("Only == and != operators are supported for strings read from a file\n");

    // INFO
    if ( !btok->nsamples )
    {
        // there is only one string value, e.g. STR[1]=@list.txt
        if ( btok->idx >= 0 )
        {
            int ret = btok->str_value.s ? khash_str2int_has_key(atok->hash, btok->str_value.s) : 0;
            if ( rtok->tok_type==TOK_NE ) ret = ret ? 0 : 1;
            rtok->pass_site = ret;
            return;
        }

        // there can be multiple comma-separated string values, e.g. STR=@list.txt or STR[*]=@list.txt
        int ret = 0;
        char *ptr = btok->str_value.s;
        while ( *ptr )
        {
            char *eptr = ptr + 1;
            while ( *eptr && *eptr!=',' ) eptr++;
            char keep = *eptr;
            *eptr = 0;
            ret |= khash_str2int_has_key(atok->hash, ptr);
            *eptr = keep;
            if ( !keep ) break;
            ptr = eptr + 1;
        }
        if ( rtok->tok_type==TOK_NE ) ret = ret ? 0 : 1;
        rtok->pass_site = ret;
        return;
    }


    // FORMAT
    tok_init_samples(atok, btok, rtok);
    rtok->pass_site = 0;
    int i;

    // there is only one string value, e.g. FMT/STR[*:1]=@list.txt
    if ( btok->idx >= 0 )
    {
        for (i=0; i<btok->nsamples; i++)
        {
            if ( !rtok->usmpl[i] ) continue;
            char *str = btok->str_value.s + i*btok->nval1;
            char keep = str[btok->nval1];
            str[btok->nval1] = 0;
            int ret = khash_str2int_has_key(atok->hash, str);
            str[btok->nval1] = keep;
            if ( rtok->tok_type==TOK_NE ) ret = ret ? 0 : 1;
            rtok->pass_samples[i] = ret;
            rtok->pass_site |= ret;
        }
        return;
    }

    // there can be multiple comma-separated string values, e.g. FMT/STR=@list.txt
    for (i=0; i<btok->nsamples; i++)
    {
        if ( !rtok->usmpl[i] ) continue;
        char *str = btok->str_value.s + i*btok->nval1;
        char keep = str[btok->nval1];
        str[btok->nval1] = 0;

        // now str contains the block of per-sample comma-separated strings to loop over
        int ret = 0;
        char *ptr = str;
        while ( *ptr )
        {
            char *eptr = ptr + 1;
            while ( *eptr && *eptr!=',' ) eptr++;
            char keep0 = *eptr;
            *eptr = 0;
            ret |= khash_str2int_has_key(atok->hash, ptr);
            *eptr = keep0;
            if ( !keep0 ) break;
            ptr = eptr + 1;
        }
        str[btok->nval1] = keep;
        if ( rtok->tok_type==TOK_NE ) ret = ret ? 0 : 1;
        rtok->pass_samples[i] = ret;
        rtok->pass_site |= ret;
    }
}
static void filters_cmp_id(token_t *atok, token_t *btok, token_t *rtok, bcf1_t *line)
{
    if ( btok->hash )
    {
        token_t *tmp = atok; atok = btok; btok = tmp;
    }

    char *id = line->d.id;
    int pass = 0;

    while ( id )
    {
        char *ep = strchr(id,';');
        if ( ep ) *ep = 0;

        if ( atok->hash )
        {
            if ( rtok->tok_type!=TOK_EQ && rtok->tok_type!=TOK_NE )
                error("Only == and != operators are supported for strings read from a file\n");

            pass = khash_str2int_has_key(atok->hash, id);
        }
        else
        {
            if ( !btok->str_value.l ) error("Error occurred while evaluating the expression\n");

            if ( rtok->tok_type==TOK_EQ || rtok->tok_type==TOK_NE )
                pass = strcmp(btok->str_value.s,id) ? 0 : 1;
            else
            {
                if ( rtok->tok_type!=TOK_LIKE && rtok->tok_type!=TOK_NLIKE )
                    error("Only the following operators are supported for querying ID: ==, !=, ~, !~; the operator type %d is not supported (%p %p)\n",
                            rtok->tok_type,atok->regex,btok->regex);

                regex_t *regex = atok->regex ? atok->regex : (btok->regex ? btok->regex : NULL);
                if ( !regex ) error("fixme: regex initialization failed\n");
                pass = regexec(regex,id, 0,NULL,0) ? 0 : 1;
            }
        }
        if ( ep )
        {
            *ep = ';';
            id = ep + 1;
        }
        if ( pass || !ep ) break;
    }
    if ( rtok->tok_type==TOK_NE || rtok->tok_type==TOK_NE) pass = pass ? 0 : 1;
    rtok->pass_site = pass;
}

/**
 *  bcf_get_info_value() - get single INFO value, int64_t or double
 *  @line:      BCF line
 *  @info_id:   tag ID, as returned by bcf_hdr_id2int
 *  @ivec:      0-based index to retrieve, -1 when single value is expected
 *  @vptr:      pointer to memory location of sufficient size to accommodate
 *              info_id's type
 *
 *  The returned value is -1 if tag is not present, 0 if present but
 *  values is missing or ivec is out of range, and 1 on success.
 */
static int bcf_get_info_value(bcf1_t *line, int info_id, int ivec, void *value)
{
    int j;
    for (j=0; j<line->n_info; j++)
        if ( line->d.info[j].key == info_id ) break;
    if ( j==line->n_info ) return -1;

    bcf_info_t *info = &line->d.info[j];
    if ( info->len == 1 )
    {
        if ( info->type==BCF_BT_FLOAT ) *((double*)value) = info->v1.f;
        else if ( info->type==BCF_BT_INT8 || info->type==BCF_BT_INT16 || info->type==BCF_BT_INT32 ) *((int64_t*)value) = info->v1.i;
        return 1;
    }

    if ( ivec<0 ) ivec = 0;
    if ( ivec>=info->len) return 0;

    #define BRANCH(type_t, convert, is_missing, is_vector_end, out_type_t) { \
        uint8_t *p = info->vptr; \
        for (j=0; j<ivec; j++) \
        { \
            type_t val = convert(&p[j * sizeof(type_t)]); \
            if ( is_vector_end ) return 0; \
        } \
        type_t val = convert(&p[ivec * sizeof(type_t)]); \
        if ( is_missing ) return 0; \
        *((out_type_t*)value) = val; \
        return 1; \
    }
    switch (info->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  le_to_i8,  val==bcf_int8_missing,  val==bcf_int8_vector_end,  int64_t); break;
        case BCF_BT_INT16: BRANCH(int16_t, le_to_i16, val==bcf_int16_missing, val==bcf_int16_vector_end, int64_t); break;
        case BCF_BT_INT32: BRANCH(int32_t, le_to_i32, val==bcf_int32_missing, val==bcf_int32_vector_end, int64_t); break;
        case BCF_BT_FLOAT: BRANCH(float,   le_to_float, bcf_float_is_missing(val), bcf_float_is_vector_end(val), double); break;
        default: fprintf(stderr,"todo: type %d\n", info->type); exit(1); break;
    }
    #undef BRANCH
    return -1;  // this shouldn't happen
}

static void filters_set_chrom(filter_t *flt, bcf1_t *line, token_t *tok)
{
    tok->str_value.l = 0;
    kputs(bcf_seqname(flt->hdr,line), &tok->str_value);
    tok->nvalues = tok->str_value.l;
    tok->is_str  = 1;
}

static void filters_set_pos(filter_t *flt, bcf1_t *line, token_t *tok)
{
    tok->values[0] = line->pos+1;
    tok->nvalues = 1;
}

static void filters_set_info_int(filter_t *flt, bcf1_t *line, token_t *tok)
{
    if ( tok->idx==-2 )
    {
        tok->nvalues = bcf_get_info_int32(flt->hdr,line,tok->tag,&flt->tmpi,&flt->mtmpi);
        if ( tok->nvalues<=0 ) tok->nvalues = 0;
        else
        {
            hts_expand(double,tok->nvalues,tok->mvalues,tok->values);
            int i, j = 0, end = tok->idxs[tok->nidxs-1] < 0 ? tok->nvalues - 1 : tok->nidxs - 1;
            if ( end >= tok->nvalues ) end = tok->nvalues - 1;
            for (i=0; i<=end; i++)
                if ( i>=tok->nidxs || tok->idxs[i] ) tok->values[j++] = flt->tmpi[i];
            tok->nvalues = j;
        }
    }
    else
    {
        int64_t value = 0;
        if ( bcf_get_info_value(line,tok->hdr_id,tok->idx,&value) <= 0 )
            tok->nvalues = 0;
        else
        {
            tok->values[0] = value;
            tok->nvalues = 1;
        }
    }
}

static void filters_set_info_float(filter_t *flt, bcf1_t *line, token_t *tok)
{
    if ( tok->idx==-2 )
    {
        tok->nvalues = bcf_get_info_float(flt->hdr,line,tok->tag,&flt->tmpf,&flt->mtmpf);
        if ( tok->nvalues<=0 ) tok->nvalues = 0;
        else
        {
            hts_expand(double,tok->nvalues,tok->mvalues,tok->values);
            int i, j = 0, end = tok->idxs[tok->nidxs-1] < 0 ? tok->nvalues - 1 : tok->nidxs - 1;
            if ( end >= tok->nvalues ) end = tok->nvalues - 1;
            for (i=0; i<=end; i++)
                if ( i>=tok->nidxs || tok->idxs[i] )
                {
                    if ( bcf_float_is_missing(flt->tmpf[i]) ) bcf_double_set_missing(tok->values[j]);
                    else tok->values[j] = flt->tmpf[i];
                    j++;
                }
            tok->nvalues = j;
        }
    }
    else
    {
        double value;
        if ( bcf_get_info_value(line,tok->hdr_id,tok->idx,&value) <= 0 )
            tok->nvalues = 0;
        else
        {
            tok->values[0] = value;
            tok->nvalues = 1;
        }
    }
}

static void filters_set_info_string(filter_t *flt, bcf1_t *line, token_t *tok)
{
    int32_t m = tok->str_value.m;
    int n = bcf_get_info_string(flt->hdr,line,tok->tag,&tok->str_value.s,&m);
    tok->str_value.m = m;
    if ( n<0 ) { tok->nvalues = tok->str_value.l = 0; return; }

    if ( tok->idx>=0 )
    {
        // get ith field (i=tok->idx)
        int i = 0;
        char *ss = tok->str_value.s, *se = tok->str_value.s + n;
        while ( ss<se && i<tok->idx )
        {
            if ( *ss==',' ) i++;
            ss++;
        }
        if ( ss==se || i!=tok->idx ) { tok->nvalues = tok->str_value.l = 0; return; }
        se = ss;
        while ( se - tok->str_value.s < n && *se!=',' ) se++;
        if ( ss==tok->str_value.s ) *se = 0;
        else
        {
            memmove(tok->str_value.s, ss, se-ss);
            tok->str_value.s[se-ss] = 0;
        }
        tok->str_value.l = se - ss;
    }
    else if ( tok->idx==-2 && tok->idxs[0]==-1 )    // keep all values, TAG[*]
        tok->str_value.l = n;
    else if ( tok->idx==-2 )
    {
        flt->tmps.l = 0;
        ks_resize(&flt->tmps, n);
        int i, iend = tok->idxs[tok->nidxs-1] < 0 ? n - 1 : tok->nidxs - 1;
        if ( iend >= n ) iend = n - 1;
        char *beg = tok->str_value.s, *dst = flt->tmps.s;
        for (i=0; i<=iend; i++)
        {
            char *end = beg;
            while ( *end && *end!=',' ) end++;
            if ( i>=tok->nidxs || tok->idxs[i] )
            {
                memcpy(dst, beg, end - beg);
                dst += end - beg;
                dst[0] = ',';
                dst++;
            }
            beg = end+1;
        }
        dst[0] = 0;
        tok->str_value.l = dst - flt->tmps.s;

        #define SWAP(type_t, a, b) { type_t t = a; a = b; b = t; }
        SWAP(char *, flt->tmps.s, tok->str_value.s);
        SWAP(size_t, flt->tmps.m, tok->str_value.m);
    }
    tok->nvalues = tok->str_value.l;
}

static void filters_set_info_flag(filter_t *flt, bcf1_t *line, token_t *tok)
{
    int j;
    for (j=0; j<line->n_info; j++)
        if ( line->d.info[j].key == tok->hdr_id ) break;
    tok->values[0] = j==line->n_info ? 0 : 1;
    tok->nvalues = 1;
}

static void filters_set_format_int(filter_t *flt, bcf1_t *line, token_t *tok)
{
    if ( line->n_sample != tok->nsamples )
        error("Incorrect number of FORMAT fields at %s:%"PRId64" .. %s, %d vs %d\n", bcf_seqname(flt->hdr,line),(int64_t) line->pos+1,tok->tag,line->n_sample,tok->nsamples);

    int nvals;
    if ( (nvals=bcf_get_format_int32(flt->hdr,line,tok->tag,&flt->tmpi,&flt->mtmpi))<0 )
    {
        tok->nvalues = 0;
        return;
    }
    int i, nsrc1 = nvals / tok->nsamples;
    tok->nval1 = tok->idx >= 0 ? 1 : (tok->nuidxs ? tok->nuidxs : nsrc1);
    tok->nvalues = tok->nval1*tok->nsamples;
    hts_expand(double, tok->nvalues, tok->mvalues, tok->values);

    if ( tok->idx >= 0 )    // scalar or vector index
    {
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) continue;
            int32_t *ptr = flt->tmpi + i*nsrc1;
            if ( tok->idx>=nsrc1 || ptr[tok->idx]==bcf_int32_missing )
                bcf_double_set_missing(tok->values[i]);
            else if ( ptr[tok->idx]==bcf_int32_vector_end )
                bcf_double_set_vector_end(tok->values[i]);
            else
                tok->values[i] = ptr[tok->idx];
        }
    }
    else if ( tok->idx==-3 )
    {
        if ( filters_cache_genotypes(flt,line)!=0 )
        {
            tok->nvalues = 0;
            return;
        }
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) continue;
            int32_t *src = flt->tmpi + i*nsrc1;
            double  *dst = tok->values + i*tok->nval1;
            int k, j = 0;
            for (k=0; k<nsrc1; k++) // source values are AD[0..nsrc1]
            {
                if ( !(flt->cached_GT.mask[i] & (1<<k)) ) continue;
                dst[j++] = src[k];
            }
            if ( !j ) { bcf_double_set_missing(dst[j]); j++; }
            for (; j<tok->nval1; j++) bcf_double_set_vector_end(dst[j]);
        }
    }
    else
    {
        int kend = tok->idxs[tok->nidxs-1] < 0 ? tok->nval1 : tok->nidxs;
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) continue;
            int32_t *src = flt->tmpi + i*nsrc1;
            double  *dst = tok->values + i*tok->nval1;
            int k, j = 0;
            for (k=0; k<kend; k++)
            {
                if ( k<tok->nidxs && !tok->idxs[k] ) continue;
                if ( src[k]==bcf_int32_missing )
                    bcf_double_set_missing(dst[j]);
                else if ( src[k]==bcf_int32_vector_end )
                    bcf_double_set_vector_end(dst[j]);
                else
                    dst[j] = src[k];
                j++;
            }
            if ( j==0 )
            {
                bcf_double_set_missing(dst[j]);
                j++;
            }
            while (j < tok->nval1)
            {
                bcf_double_set_vector_end(dst[j]);
                j++;
            }
        }
    }
}
static void filters_set_format_float(filter_t *flt, bcf1_t *line, token_t *tok)
{
    if ( line->n_sample != tok->nsamples )
        error("Incorrect number of FORMAT fields at %s:%"PRId64" .. %s, %d vs %d\n", bcf_seqname(flt->hdr,line),(int64_t) line->pos+1,tok->tag,line->n_sample,tok->nsamples);

    int nvals;
    if ( (nvals=bcf_get_format_float(flt->hdr,line,tok->tag,&flt->tmpf,&flt->mtmpf))<0 )
    {
        tok->nvalues = 0;
        return;
    }
    int i, nsrc1 = nvals / tok->nsamples;
    tok->nval1 = tok->idx >= 0 ? 1 : (tok->nuidxs ? tok->nuidxs : nsrc1);
    tok->nvalues = tok->nval1*tok->nsamples;
    hts_expand(double, tok->nvalues, tok->mvalues, tok->values);

    if ( tok->idx >= 0 )    // scalar or vector index
    {
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) continue;
            float *ptr = flt->tmpf + i*nsrc1;
            if ( tok->idx>=nsrc1 || bcf_float_is_missing(ptr[tok->idx]) )
                bcf_double_set_missing(tok->values[i]);
            else if ( bcf_float_is_vector_end(ptr[tok->idx]) )
                bcf_double_set_vector_end(tok->values[i]);
            else
                tok->values[i] = ptr[tok->idx];
        }
    }
    else if ( tok->idx==-3 )
    {
        if ( filters_cache_genotypes(flt,line)!=0 )
        {
            tok->nvalues = 0;
            return;
        }
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) continue;
            float  *src = flt->tmpf + i*nsrc1;
            double *dst = tok->values + i*tok->nval1;
            int k, j = 0;
            for (k=0; k<nsrc1; k++) // source values are AF[0..nsrc1]
            {
                if ( !(flt->cached_GT.mask[i] & (1<<k)) ) continue;
                if ( bcf_float_is_missing(src[k]) )
                    bcf_double_set_missing(dst[j]);
                else if ( bcf_float_is_vector_end(src[k]) )
                    bcf_double_set_vector_end(dst[j]);
                else
                    dst[j] = src[k];
                j++;
            }
            if ( !j ) { bcf_double_set_missing(dst[j]); j++; }
            for (; j<tok->nval1; j++) bcf_double_set_vector_end(dst[j]);
        }
    }
    else
    {
        int kend = tok->idxs[tok->nidxs-1] < 0 ? tok->nval1 : tok->nidxs;
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) continue;
            float *src = flt->tmpf + i*nsrc1;
            double  *dst = tok->values + i*tok->nval1;
            int k, j = 0;
            for (k=0; k<kend; k++)
            {
                if ( k<tok->nidxs && !tok->idxs[k] ) continue;
                if ( bcf_float_is_missing(src[k]) )
                    bcf_double_set_missing(dst[j]);
                else if ( bcf_float_is_vector_end(src[k]) )
                    bcf_double_set_vector_end(dst[j]);
                else
                    dst[j] = src[k];
                j++;
            }
            if ( j==0 )
            {
                bcf_double_set_missing(dst[j]);
                j++;
            }
            while (j < tok->nval1)
            {
                bcf_double_set_vector_end(dst[j]);
                j++;
            }
        }
    }
}
static void filters_set_format_string(filter_t *flt, bcf1_t *line, token_t *tok)
{
    if ( line->n_sample != tok->nsamples )
        error("Incorrect number of FORMAT fields at %s:%"PRId64" .. %s, %d vs %d\n", bcf_seqname(flt->hdr,line),(int64_t) line->pos+1,tok->tag,line->n_sample,tok->nsamples);

    int i, ndim = tok->str_value.m;
    int nstr = bcf_get_format_char(flt->hdr, line, tok->tag, &tok->str_value.s, &ndim);
    tok->str_value.m = tok->str_value.l = ndim;
    kputc(0,&tok->str_value); // append the nul byte
    tok->str_value.l = tok->nvalues = 0;

    if ( nstr<0 ) return;

    if ( tok->idx==-3 && filters_cache_genotypes(flt,line)!=0 ) return;

    tok->nvalues = tok->str_value.l = nstr;
    tok->nval1   = nstr / tok->nsamples;
    for (i=0; i<tok->nsamples; i++)
    {
        if ( !tok->usmpl[i] ) continue;
        char *src = tok->str_value.s + i*tok->nval1, *dst = src;
        int ibeg = 0;
        int idx = 0;        // idx-th field
        while ( ibeg < tok->nval1 )
        {
            int iend = ibeg;
            while ( iend < tok->nval1 && src[iend] && src[iend]!=',' ) iend++;

            int keep = 0;
            if ( tok->idx >= 0 )    // single index, given explicitly, e.g. AD[:1]
            {
                if ( tok->idx==idx ) keep = 1;
            }
            else if ( tok->idx == -3 )  // given by GT index, e.g. AD[:GT]
            {
                if ( flt->cached_GT.mask[i] & (1<<idx) ) keep = 1;
            }
            else    // given as a list, e.g. AD[:0,3]
            {
                if ( idx < tok->nidxs )
                {
                    if ( tok->idxs[idx] != 0 ) keep = 1;
                }
                else if ( tok->idxs[tok->nidxs-1] < 0 )
                    keep = 1;
            }
            if ( keep )
            {
                if ( ibeg!=0 ) memmove(dst, src+ibeg, iend-ibeg+1);
                dst += iend - ibeg + 1;
                if ( tok->idx>=0 ) break;
            }
            if ( !src[iend] ) break;
            ibeg = iend + 1;
            idx++;
        }
        if ( dst==src ) { dst[0] = '.'; dst+=2; }
        if ( dst - src < tok->nval1 ) memset(dst-1, 0, tok->nval1 - (dst - src));
    }
}
static void _filters_set_genotype(filter_t *flt, bcf1_t *line, token_t *tok, int type)
{
    bcf_fmt_t *fmt = bcf_get_fmt(flt->hdr, line, "GT");
    if ( !fmt )
    {
        tok->nvalues = tok->str_value.l = 0;
        return;
    }

    int i,j, nsmpl = bcf_hdr_nsamples(flt->hdr), nvals1 = type==2 ? 3 : 4;
    if ( tok->str_value.m <= nvals1*nsmpl )
    {
        tok->str_value.m = nvals1*nsmpl + 1;
        tok->str_value.s = (char*)realloc(tok->str_value.s, tok->str_value.m);
    }

#define BRANCH_INT(type_t, convert, vector_end) \
    { \
        for (i=0; i<line->n_sample; i++) \
        { \
            uint8_t *ptr = fmt->p + i*fmt->size; \
            int is_het = 0, has_ref = 0, missing = 0, jal = 0; \
            for (j=0; j<fmt->n; j++) \
            { \
                type_t ial = convert(&ptr[j * sizeof(type_t)]); \
                if ( ial==vector_end ) break; /* smaller ploidy */ \
                if ( bcf_gt_is_missing(ial) ) { missing=1; break; } /* missing allele */ \
                if ( bcf_gt_allele(ial)==0 ) has_ref = 1; \
                if ( j>0 ) \
                { \
                    if ( bcf_gt_allele(ial)!=bcf_gt_allele(jal) ) is_het = 1; \
                } \
                jal = ial; \
            } \
            char *dst = &tok->str_value.s[nvals1*i]; \
            if ( type==4 ) \
            { \
                if ( !j || missing ) dst[0]='m', dst[1]='i', dst[2]='s', dst[3] = 0; /* mis, missing genotype */ \
                else if ( !has_ref ) dst[0]='a', dst[1]='l', dst[2]='t', dst[3] = 0; /* alt, no ref, must have alt allele */ \
                else if ( !is_het ) dst[0]='r', dst[1]='e', dst[2]='f', dst[3] = 0; /* ref, must be ref-only, no alt alelle */ \
                else dst[0]='a', dst[1]='l', dst[2]='t', dst[3] = 0; /* alt, is het, has alt allele */ \
            } \
            else if ( !j || missing ) dst[0]='.', dst[1]=0; /* ., missing genotype */ \
            else if ( type==3 ) \
            { \
                if ( j==1 ) dst[0]='h', dst[1]='a', dst[2]='p', dst[3] = 0; /* hap, haploid */ \
                else if ( !is_het ) dst[0]='h', dst[1]='o', dst[2]='m', dst[3] = 0; /* hom */ \
                else dst[0]='h', dst[1]='e', dst[2]='t', dst[3] = 0; /* het */ \
            } \
            else \
            { \
                if ( j==1 ) \
                { \
                    if ( has_ref ) dst[0]='r', dst[1]=0; /* r, haploid */ \
                    else dst[0]='a', dst[1]=0; /* a, haploid */ \
                } \
                else if ( !is_het ) \
                { \
                    if ( has_ref ) dst[0]='r', dst[1]='r', dst[2] = 0; /* rr */ \
                    else dst[0]='a', dst[1]='a', dst[2] = 0; /* aa */ \
                } \
                else \
                { \
                    if ( has_ref ) dst[0]='r', dst[1]='a', dst[2] = 0; /* ra */ \
                    else dst[0]='a', dst[1]='A', dst[2] = 0; /* aA */ \
                } \
            } \
        } \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH_INT(int8_t,  le_to_i8,  bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH_INT(int16_t, le_to_i16, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH_INT(int32_t, le_to_i32, bcf_int32_vector_end); break;
        default: error("The GT type is not lineognised: %d at %s:%"PRId64"\n",fmt->type, bcf_seqname(flt->hdr,line),(int64_t) line->pos+1); break;
    }
#undef BRANCH_INT
    assert( tok->nsamples == nsmpl );
    tok->nvalues = tok->str_value.l = nvals1*nsmpl;
    tok->str_value.s[tok->str_value.l] = 0;
    tok->nval1 = nvals1;
}
static void filters_set_genotype2(filter_t *flt, bcf1_t *line, token_t *tok) { _filters_set_genotype(flt, line, tok, 2); }  // rr, ra, aa, aA etc
static void filters_set_genotype3(filter_t *flt, bcf1_t *line, token_t *tok) { _filters_set_genotype(flt, line, tok, 3); }  // hap, hom, het
static void filters_set_genotype4(filter_t *flt, bcf1_t *line, token_t *tok) { _filters_set_genotype(flt, line, tok, 4); }  // mis, alt, ref

static void filters_set_genotype_string(filter_t *flt, bcf1_t *line, token_t *tok)
{
    bcf_fmt_t *fmt = bcf_get_fmt(flt->hdr, line, "GT");
    if ( !fmt )
    {
        tok->nvalues = 0;
        return;
    }
    int i, blen = 4, nsmpl = line->n_sample;

gt_length_too_big:
    tok->str_value.l = 0;
    for (i=0; i<nsmpl; i++)
    {
        size_t plen = tok->str_value.l;
        bcf_format_gt(fmt, i, &tok->str_value);
        kputc_(0, &tok->str_value);
        if ( tok->str_value.l - plen > blen )
        {
            // too many alternate alleles or ploidy is too large, the genotype does not fit
            // three characters ("0/0" vs "10/10").
            blen *= 2;
            goto gt_length_too_big;
        }

        plen = tok->str_value.l - plen;
        while ( plen < blen )
        {
            kputc_(0, &tok->str_value);
            plen++;
        }
    }
    assert( tok->nsamples == nsmpl );
    tok->nvalues = tok->str_value.l;
    tok->nval1 = blen;
}
static void filters_set_ilen(filter_t *flt, bcf1_t *line, token_t *tok)
{
    tok->nvalues = line->n_allele - 1;
    hts_expand(double,tok->nvalues,tok->mvalues,tok->values);

    int i, rlen = strlen(line->d.allele[0]);
    for (i=1; i<line->n_allele; i++)
    {
        if ( line->d.allele[i][0]=='<' )
        {
            bcf_double_set_missing(tok->values[i-1]);
            continue;
        }
        int alen = strlen(line->d.allele[i]);
        tok->values[i-1] = alen - rlen;
    }
}
static void filters_set_ref_string(filter_t *flt, bcf1_t *line, token_t *tok)
{
    tok->str_value.l = 0;
    kputs(line->d.allele[0], &tok->str_value);
    tok->nvalues = tok->str_value.l;
}
static void filters_set_alt_string(filter_t *flt, bcf1_t *line, token_t *tok)
{
    tok->str_value.l = 0;
    if ( tok->idx>=0 )
    {
        if ( line->n_allele > tok->idx + 1 )
            kputs(line->d.allele[tok->idx + 1], &tok->str_value);
        else
            kputc('.', &tok->str_value);
    }
    else if ( tok->idx==-2 )
    {
        int i, end = tok->nuidxs ? tok->nuidxs : line->n_allele - 1;
        if ( end >= line->n_allele - 1 ) end = line->n_allele - 2;
        for (i=0; i<=end; i++)
            if ( i>=tok->nidxs || tok->idxs[i] )
            {
                if ( tok->str_value.l ) kputc(',', &tok->str_value);
                kputs(line->d.allele[i+1], &tok->str_value);
            }
    }
    else if ( line->n_allele>1 )
    {
        kputs(line->d.allele[1], &tok->str_value);
        int i;
        for (i=2; i<line->n_allele; i++)
        {
            kputc(',', &tok->str_value);
            kputs(line->d.allele[i], &tok->str_value);
        }
    }
    else if ( line->n_allele==1 )
        kputc('.', &tok->str_value);
    tok->nvalues = tok->str_value.l;
}
static void filters_set_nmissing(filter_t *flt, bcf1_t *line, token_t *tok)
{
    bcf_unpack(line, BCF_UN_FMT);
    if ( !line->n_sample )
    {
        tok->nvalues = 1;
        tok->values[0] = 0;
        return;
    }

    int i,igt = bcf_hdr_id2int(flt->hdr, BCF_DT_ID, "GT");
    bcf_fmt_t *fmt = NULL;
    for (i=0; i<line->n_fmt; i++)
        if ( line->d.fmt[i].id==igt ) { fmt = &line->d.fmt[i]; break; }
    if ( !fmt )
    {
        tok->nvalues = 0;
        return;
    }

    int j,nmissing = 0;
    #define BRANCH(type_t, convert, is_vector_end) { \
        for (i=0; i<line->n_sample; i++) \
        { \
            uint8_t *ptr = fmt->p + i*fmt->size; \
            for (j=0; j<fmt->n; j++) \
            { \
                type_t val = convert(&ptr[j * sizeof(type_t)]); \
                if ( val==is_vector_end ) break; \
                if ( val==bcf_gt_missing ) { nmissing++; break; } \
            } \
        } \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  le_to_i8,  bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH(int16_t, le_to_i16, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH(int32_t, le_to_i32, bcf_int32_vector_end); break;
        default: fprintf(stderr,"todo: type %d\n", fmt->type); exit(1); break;
    }
    #undef BRANCH
    tok->nvalues = 1;
    tok->values[0] = tok->tag[0]=='N' ? nmissing : (double)nmissing / line->n_sample;
}
static int func_npass(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    if ( nstack==0 ) error("Error parsing the expression\n");
    token_t *tok = stack[nstack - 1];
    if ( !tok->nsamples ) error("The function %s works with FORMAT fields\n", rtok->tag);
    assert(tok->usmpl);

    int i, npass = 0;
    for (i=0; i<tok->nsamples; i++)
    {
        if ( !tok->usmpl[i] ) continue;
        if ( tok->pass_samples[i] ) npass++;
    }
    hts_expand(double,1,rtok->mvalues,rtok->values);
    rtok->nsamples = 0;
    rtok->nvalues = 1;
    rtok->values[0] = rtok->tag[0]=='N' ? npass : (line->n_sample ? 1.0*npass/line->n_sample : 0);

    return 1;
}
static void filters_set_nalt(filter_t *flt, bcf1_t *line, token_t *tok)
{
    tok->nvalues = 1;
    tok->values[0] = line->n_allele - 1;
}
static void filters_set_ac(filter_t *flt, bcf1_t *line, token_t *tok)
{
    hts_expand(int32_t, line->n_allele, flt->mtmpi, flt->tmpi);
    if ( !bcf_calc_ac(flt->hdr, line, flt->tmpi, BCF_UN_INFO|BCF_UN_FMT) )
    {
        tok->nvalues = 0;
        return;
    }
    int i, an = flt->tmpi[0];
    for (i=1; i<line->n_allele; i++) an += flt->tmpi[i];
    if ( !an )
    {
        tok->nvalues = 0;
        return;
    }
    flt->tmpi[0] = an;  // for filters_set_[mac|af|maf]
    if ( tok->idx>=0 )
    {
        tok->nvalues = 1;
        tok->values[0] = tok->idx+1<line->n_allele ? flt->tmpi[tok->idx+1] : 0;
    }
    else if ( line->n_allele==1 )   // no ALT
    {
        tok->nvalues = 1;
        tok->values[0] = 0;
    }
    else
    {
        hts_expand(double,line->n_allele,tok->mvalues,tok->values);
        for (i=1; i<line->n_allele; i++)
            tok->values[i-1] = flt->tmpi[i];
        tok->nvalues = line->n_allele - 1;
    }
}
static void filters_set_an(filter_t *flt, bcf1_t *line, token_t *tok)
{
    filters_set_ac(flt,line,tok);
    tok->values[0] = tok->nvalues ? flt->tmpi[0] : 0;
    tok->nvalues = 1;
}
static void filters_set_mac(filter_t *flt, bcf1_t *line, token_t *tok)
{
    filters_set_ac(flt,line,tok);
    if ( !tok->nvalues ) return;
    int i, an = flt->tmpi[0];
    for (i=0; i<tok->nvalues; i++)
        if ( tok->values[i] > an*0.5 ) tok->values[i] = an - tok->values[i];
}
static void filters_set_af(filter_t *flt, bcf1_t *line, token_t *tok)
{
    filters_set_ac(flt,line,tok);
    if ( !tok->nvalues ) return;
    int i, an = flt->tmpi[0];
    for (i=0; i<tok->nvalues; i++)
        tok->values[i] /= (double)an;
}
static void filters_set_maf(filter_t *flt, bcf1_t *line, token_t *tok)
{
    filters_set_ac(flt,line,tok);
    if ( !tok->nvalues ) return;
    int i, an = flt->tmpi[0];
    for (i=0; i<tok->nvalues; i++)
    {
        tok->values[i] /= (double)an;
        if ( tok->values[i] > 0.5 ) tok->values[i] = 1 - tok->values[i];
    }
}

static int func_max(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    rtok->nvalues = 0;
    if ( !tok->nvalues ) return 1;
    double *ptr, val = -HUGE_VAL;
    int i,j, has_value = 0;
    if ( tok->nsamples )
    {
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) continue;
            ptr = tok->values + i*tok->nval1;
            for (j=0; j<tok->nval1; j++)
            {
                if ( bcf_double_is_missing_or_vector_end(ptr[j]) ) continue;
                has_value = 1;
                if ( val < ptr[j] ) val = ptr[j];
            }
        }
    }
    else
    {
        for (i=0; i<tok->nvalues; i++)
        {
            if ( bcf_double_is_missing_or_vector_end(tok->values[i]) ) continue;
            has_value = 1;
            if ( val < tok->values[i] ) val = tok->values[i];
        }
    }
    if ( has_value )
    {
        rtok->values[0] = val;
        rtok->nvalues = has_value;
    }
    return 1;
}
static int func_smpl_max(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    if ( !tok->nsamples ) return func_max(flt,line,rtok,stack,nstack);
    rtok->nsamples = tok->nsamples;
    rtok->nvalues  = tok->nsamples;
    rtok->nval1 = 1;
    hts_expand(double,rtok->nvalues,rtok->mvalues,rtok->values);
    assert(tok->usmpl);
    if ( !rtok->usmpl ) rtok->usmpl = (uint8_t*) malloc(tok->nsamples);
    memcpy(rtok->usmpl, tok->usmpl, tok->nsamples);
    int i, j, has_value;
    double val, *ptr;
    for (i=0; i<tok->nsamples; i++)
    {
        if ( !rtok->usmpl[i] ) continue;
        val = -HUGE_VAL;
        has_value = 0;
        ptr = tok->values + i*tok->nval1;
        for (j=0; j<tok->nval1; j++)
        {
            if ( bcf_double_is_missing_or_vector_end(ptr[j]) ) continue;
            has_value = 1;
            if ( val < ptr[j] ) val = ptr[j];
        }
        if ( has_value ) rtok->values[i] = val;
        else bcf_double_set_missing(rtok->values[i]);
    }
    return 1;
}
static int func_min(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    rtok->nvalues = 0;
    if ( !tok->nvalues ) return 1;
    double *ptr, val = HUGE_VAL;
    int i,j, has_value = 0;
    if ( tok->nsamples )
    {
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) continue;
            ptr = tok->values + i*tok->nval1;
            for (j=0; j<tok->nval1; j++)
            {
                if ( bcf_double_is_missing_or_vector_end(ptr[j]) ) continue;
                has_value = 1;
                if ( val > ptr[j] ) val = ptr[j];
            }
        }
    }
    else
    {
        for (i=0; i<tok->nvalues; i++)
        {
            if ( bcf_double_is_missing_or_vector_end(tok->values[i]) ) continue;
            has_value = 1;
            if ( val > tok->values[i] ) val = tok->values[i];
        }
    }
    if ( has_value )
    {
        rtok->values[0] = val;
        rtok->nvalues = has_value;
    }
    return 1;
}
static int func_smpl_min(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    if ( !tok->nsamples ) return func_min(flt,line,rtok,stack,nstack);
    rtok->nsamples = tok->nsamples;
    rtok->nvalues  = tok->nsamples;
    rtok->nval1 = 1;
    hts_expand(double,rtok->nvalues,rtok->mvalues,rtok->values);
    assert(tok->usmpl);
    if ( !rtok->usmpl ) rtok->usmpl = (uint8_t*) malloc(tok->nsamples);
    memcpy(rtok->usmpl, tok->usmpl, tok->nsamples);
    int i, j, has_value;
    double val, *ptr;
    for (i=0; i<tok->nsamples; i++)
    {
        if ( !rtok->usmpl[i] ) continue;
        val = HUGE_VAL;
        has_value = 0;
        ptr = tok->values + i*tok->nval1;
        for (j=0; j<tok->nval1; j++)
        {
            if ( bcf_double_is_missing_or_vector_end(ptr[j]) ) continue;
            has_value = 1;
            if ( val > ptr[j] ) val = ptr[j];
        }
        if ( has_value ) rtok->values[i] = val;
        else bcf_double_set_missing(rtok->values[i]);
    }
    return 1;
}
static int func_avg(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    rtok->nvalues = 0;
    if ( !tok->nvalues ) return 1;
    double *ptr, val = 0;
    int i,j, n = 0;
    if ( tok->nsamples )
    {
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) continue;
            ptr = tok->values + i*tok->nval1;
            for (j=0; j<tok->nval1; j++)
            {
                if ( bcf_double_is_missing_or_vector_end(ptr[j]) ) continue;
                val += ptr[j];
                n++;
            }
        }
    }
    else
    {
        for (i=0; i<tok->nvalues; i++)
            if ( !bcf_double_is_missing_or_vector_end(tok->values[i]) ) { val += tok->values[i]; n++; }
    }
    if ( n )
    {
        rtok->values[0] = val / n;
        rtok->nvalues   = 1;
    }
    return 1;
}
static int func_smpl_avg(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    if ( !tok->nsamples ) return func_avg(flt,line,rtok,stack,nstack);
    rtok->nsamples = tok->nsamples;
    rtok->nvalues  = tok->nsamples;
    rtok->nval1 = 1;
    hts_expand(double,rtok->nvalues,rtok->mvalues,rtok->values);
    assert(tok->usmpl);
    if ( !rtok->usmpl ) rtok->usmpl = (uint8_t*) malloc(tok->nsamples);
    memcpy(rtok->usmpl, tok->usmpl, tok->nsamples);
    int i, j, n;
    double val, *ptr;
    for (i=0; i<tok->nsamples; i++)
    {
        if ( !rtok->usmpl[i] ) continue;
        val = 0;
        n = 0;
        ptr = tok->values + i*tok->nval1;
        for (j=0; j<tok->nval1; j++)
        {
            if ( !bcf_double_is_missing_or_vector_end(ptr[j]) ) { val += ptr[j]; n++; }
        }
        if ( n ) rtok->values[i] = val / n;
        else bcf_double_set_missing(rtok->values[i]);
    }
    return 1;
}
static int compare_doubles(const void *lhs, const void *rhs)
{
    double arg1 = *(const double*) lhs;
    double arg2 = *(const double*) rhs;
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}
static int func_median(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    rtok->nvalues = 0;
    if ( !tok->nvalues ) return 1;
    // sweep through all tok->values and while excluding all missing values reuse the very same array
    int i,j,k = 0, n = 0;
    if ( tok->nsamples )
    {
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) { k += tok->nval1; continue; }
            for (j=0; j<tok->nval1; k++,j++)
            {
                if ( bcf_double_is_missing_or_vector_end(tok->values[k]) ) continue;
                if ( n < k ) tok->values[n] = tok->values[k];
                n++;
            }
        }
    }
    else
    {
        for (i=0; i<tok->nvalues; i++)
        {
            if ( bcf_double_is_missing_or_vector_end(tok->values[i]) ) continue;
            if ( n < i ) tok->values[n] = tok->values[i];
            n++;
        }
    }
    if ( !n ) return 1;
    if ( n==1 ) rtok->values[0] = tok->values[0];
    else
    {
        qsort(tok->values, n, sizeof(double), compare_doubles);
        rtok->values[0] = n % 2 ? tok->values[n/2] : (tok->values[n/2-1] + tok->values[n/2]) * 0.5;
    }
    rtok->nvalues = 1;
    return 1;
}
static int func_smpl_median(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    if ( !tok->nsamples ) return func_avg(flt,line,rtok,stack,nstack);
    rtok->nsamples = tok->nsamples;
    rtok->nvalues  = tok->nsamples;
    rtok->nval1 = 1;
    hts_expand(double,rtok->nvalues,rtok->mvalues,rtok->values);
    assert(tok->usmpl);
    if ( !rtok->usmpl ) rtok->usmpl = (uint8_t*) malloc(tok->nsamples);
    memcpy(rtok->usmpl, tok->usmpl, tok->nsamples);
    int i, j, n;
    double *ptr;
    for (i=0; i<tok->nsamples; i++)
    {
        if ( !rtok->usmpl[i] ) continue;
        n = 0;
        ptr = tok->values + i*tok->nval1;
        for (j=0; j<tok->nval1; j++)
        {
            if ( bcf_double_is_missing_or_vector_end(ptr[j]) ) continue;
            if ( n < j ) ptr[n] = ptr[j];
            n++;
        }
        if ( n==0 )
            bcf_double_set_missing(rtok->values[i]);
        else if ( n==1 )
            rtok->values[i] = ptr[0];
        else
        {
            qsort(ptr, n, sizeof(double), compare_doubles);
            rtok->values[i] = n % 2 ? ptr[n/2] : (ptr[n/2-1] + ptr[n/2]) * 0.5;
        }
    }
    return 1;
}
static int func_stddev(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    rtok->nvalues = 0;
    if ( !tok->nvalues ) return 1;
    // sweep through all tok->values and while excluding all missing values reuse the very same array
    int i,j,k = 0, n = 0;
    if ( tok->nsamples )
    {
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) { k += tok->nval1; continue; }
            for (j=0; j<tok->nval1; k++,j++)
            {
                if ( bcf_double_is_missing_or_vector_end(tok->values[k]) ) continue;
                if ( n < k ) tok->values[n] = tok->values[k];
                n++;
            }
        }
    }
    else
    {
        for (i=0; i<tok->nvalues; i++)
        {
            if ( bcf_double_is_missing_or_vector_end(tok->values[i]) ) continue;
            if ( n < i ) tok->values[n] = tok->values[i];
            n++;
        }
    }
    if ( !n ) return 1;
    if ( n==1 ) rtok->values[0] = 0;
    else
    {
        double sdev = 0, avg = 0;
        for (i=0; i<n; i++) avg += tok->values[i];
        avg /= n;
        for (i=0; i<n; i++) sdev += (tok->values[i] - avg) * (tok->values[i] - avg);
        rtok->values[0] = sqrt(sdev/n);
    }
    rtok->nvalues = 1;
    return 1;
}
static int func_smpl_stddev(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    if ( !tok->nsamples ) return func_avg(flt,line,rtok,stack,nstack);
    rtok->nsamples = tok->nsamples;
    rtok->nvalues  = tok->nsamples;
    rtok->nval1 = 1;
    hts_expand(double,rtok->nvalues,rtok->mvalues,rtok->values);
    assert(tok->usmpl);
    if ( !rtok->usmpl ) rtok->usmpl = (uint8_t*) malloc(tok->nsamples);
    memcpy(rtok->usmpl, tok->usmpl, tok->nsamples);
    int i, j, n;
    double *ptr;
    for (i=0; i<tok->nsamples; i++)
    {
        if ( !rtok->usmpl[i] ) continue;
        n = 0;
        ptr = tok->values + i*tok->nval1;
        for (j=0; j<tok->nval1; j++)
        {
            if ( bcf_double_is_missing_or_vector_end(ptr[j]) ) continue;
            if ( n < j ) ptr[n] = ptr[j];
            n++;
        }
        if ( n==0 )
            bcf_double_set_missing(rtok->values[i]);
        else if ( n==1 )
            rtok->values[i] = 0;
        else
        {
            double sdev = 0, avg = 0;
            for (j=0; j<n; j++) avg += ptr[j];
            avg /= n;
            for (j=0; j<n; j++) sdev += (ptr[j] - avg) * (ptr[j] - avg);
            rtok->values[i] = sqrt(sdev/n);
        }
    }
    return 1;
}
static int func_sum(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    rtok->nvalues = 0;
    token_t *tok = stack[nstack - 1];
    if ( !tok->nvalues ) return 1;
    double *ptr, val = 0;
    int i,j, n = 0;
    if ( tok->nsamples )
    {
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) continue;
            ptr = tok->values + i*tok->nval1;
            for (j=0; j<tok->nval1; j++)
            {
                if ( bcf_double_is_missing_or_vector_end(ptr[j]) ) continue;
                val += ptr[j];
                n++;
            }
        }
    }
    else
    {
        for (i=0; i<tok->nvalues; i++)
            if ( !bcf_double_is_missing_or_vector_end(tok->values[i]) ) { val += tok->values[i]; n++; }
    }
    if ( n )
    {
        rtok->values[0] = val;
        rtok->nvalues   = 1;
    }
    return 1;
}
static int func_smpl_sum(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    if ( !tok->nsamples ) return func_avg(flt,line,rtok,stack,nstack);
    rtok->nsamples = tok->nsamples;
    rtok->nvalues  = tok->nsamples;
    rtok->nval1 = 1;
    hts_expand(double,rtok->nvalues,rtok->mvalues,rtok->values);
    assert(tok->usmpl);
    if ( !rtok->usmpl ) rtok->usmpl = (uint8_t*) malloc(tok->nsamples);
    memcpy(rtok->usmpl, tok->usmpl, tok->nsamples);
    int i, j, has_value;
    double val, *ptr;
    for (i=0; i<tok->nsamples; i++)
    {
        if ( !rtok->usmpl[i] ) continue;
        val = 0;
        has_value = 0;
        ptr = tok->values + i*tok->nval1;
        for (j=0; j<tok->nval1; j++)
        {
            if ( bcf_double_is_missing_or_vector_end(ptr[j]) ) continue;
            has_value = 1;
            val += ptr[j];
        }
        if ( has_value ) rtok->values[i] = val;
        else bcf_double_set_missing(rtok->values[i]);
    }
    return 1;
}
static int func_abs(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    if ( tok->is_str ) error("ABS() can be applied only on numeric values\n");
    rtok->nsamples = tok->nsamples;
    rtok->nvalues = tok->nvalues;
    rtok->nval1 = tok->nval1;
    hts_expand(double,rtok->nvalues,rtok->mvalues,rtok->values);
    if ( tok->usmpl )
    {
        if ( !rtok->usmpl ) rtok->usmpl = (uint8_t*) malloc(tok->nsamples);
        memcpy(rtok->usmpl, tok->usmpl, tok->nsamples);
    }
    if ( !tok->nvalues ) return 1;
    hts_expand(double, rtok->nvalues, rtok->mvalues, rtok->values);
    int i,j,k = 0;
    if ( tok->usmpl )
    {
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) { k+= tok->nval1; continue; }
            for (j=0; j<tok->nval1; k++,j++)
            {
                if ( bcf_double_is_missing_or_vector_end(tok->values[k]) ) bcf_double_set_missing(rtok->values[k]);
                else rtok->values[k] = fabs(tok->values[k]);
            }
        }
    }
    else
    {
        for (i=0; i<tok->nvalues; i++)
        {
            if ( tok->usmpl && !tok->usmpl[i] ) continue;
            if ( bcf_double_is_missing(tok->values[i]) ) bcf_double_set_missing(rtok->values[i]);
            else if ( !bcf_double_is_vector_end(tok->values[i]) ) rtok->values[i] = fabs(tok->values[i]);
        }
    }
    return 1;
}
static int func_count(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    int i,j, cnt = 0;
    if ( tok->tag && tok->nsamples )
    {
        // raw number of values in a FMT tag, e.g. COUNT(FMT/TAG)
        if ( tok->is_str ) error("todo: Type=String for COUNT on FORMAT fields?\n");
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) continue;
            double *ptr = tok->values + i*tok->nval1;
            for (j=0; j<tok->nval1; j++)
                if ( !bcf_double_is_missing_or_vector_end(ptr[j]) ) cnt++;
        }
    }
    else if ( tok->nsamples )
    {
        // number of samples that pass a processed FMT tag
        for (i=0; i<tok->nsamples; i++)
            if ( tok->pass_samples[i] ) cnt++;
    }
    else if ( tok->is_str )
    {
        if ( tok->str_value.l ) cnt = 1;
        for (i=0; i<tok->str_value.l; i++) if ( tok->str_value.s[i]==',' ) cnt++;
    }
    else
        cnt = tok->nvalues;

    rtok->nvalues = 1;
    rtok->values[0] = cnt;
    return 1;
}
static int func_strlen(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    rtok->nvalues = rtok->str_value.l = 0;
    if ( !tok->str_value.l ) return 1;

    if ( tok->idx==-2 )
    {
        int i = 0;
        char *ss = tok->str_value.s;
        while ( *ss )
        {
            char *se = ss;
            while ( *se && *se!=',' ) se++;
            hts_expand(double, i+1, rtok->mvalues, rtok->values);
            if ( !*se ) rtok->values[i] = strlen(ss);
            else
            {
                *se = 0;
                rtok->values[i] = strlen(ss);
                *se = ',';
            }
            ss = *se ? se + 1 : se;
            i++;
        }
        rtok->nvalues = i;
    }
    else
    {
        if ( !tok->str_value.s[1] && tok->str_value.s[0]=='.' )
            rtok->values[0] = 0;
        else
            rtok->values[0] = strlen(tok->str_value.s);
        rtok->nvalues = 1;
    }
    return 1;
}
static int func_binom(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    int i, istack = nstack - rtok->nargs;

    if ( rtok->nargs!=2 && rtok->nargs!=1 ) error("Error: binom() takes one or two arguments\n");
    assert( istack>=0 );

    // The expected mean is 0.5. Should we support also prob!=0.5?
    //
    //  double prob = 0.5;
    //  if ( nstack==3 )
    //  {
    //      // three parameters, the first one must be a scalar:  binom(0.25,AD[0],AD[1])
    //      if ( !stack[istack]->is_constant )
    //          error("The first argument of binom() must be a numeric constant if three parameters are given\n");
    //      prob = stack[istack]->threshold;
    //      istack++;
    //  }
    //  else if ( nstack==2 && stack[istack]->is_constant )
    //  {
    //      // two parameters, the first can be a scalar:  binom(0.25,AD) or binom(AD[0],AD[1])
    //      prob = stack[istack]->threshold;
    //      istack++;
    //  }

    token_t *tok = stack[istack];
    if ( tok->nsamples )
    {
        // working with a FORMAT tag
        rtok->nval1    = 1;
        rtok->nvalues  = tok->nsamples;
        rtok->nsamples = tok->nsamples;
        hts_expand(double, rtok->nvalues, rtok->mvalues, rtok->values);
        assert(tok->usmpl);
        if ( !rtok->usmpl ) rtok->usmpl = (uint8_t*) malloc(tok->nsamples);
        memcpy(rtok->usmpl, tok->usmpl, tok->nsamples);

        if ( istack+1==nstack )
        {
            // determine the index from the GT field: binom(AD)
            int ngt = bcf_get_genotypes(flt->hdr, line, &flt->tmpi, &flt->mtmpi);
            int max_ploidy = ngt/line->n_sample;
            if ( ngt <= 0 || max_ploidy < 2 ) // GT not present or not diploid, cannot set
            {
                for (i=0; i<rtok->nsamples; i++)
                    if ( rtok->usmpl[i] ) bcf_double_set_missing(rtok->values[i]);
                return rtok->nargs;
            }
            for (i=0; i<rtok->nsamples; i++)
            {
                if ( !rtok->usmpl[i] ) continue;
                int32_t *ptr = flt->tmpi + i*max_ploidy;
                if ( bcf_gt_is_missing(ptr[0]) || bcf_gt_is_missing(ptr[1]) || ptr[1]==bcf_int32_vector_end )
                {
                    bcf_double_set_missing(rtok->values[i]);
                    continue;
                }
                int idx1 = bcf_gt_allele(ptr[0]);
                int idx2 = bcf_gt_allele(ptr[1]);
                if ( idx1>=line->n_allele ) error("Incorrect allele index at %s:%"PRId64", sample %s\n", bcf_seqname(flt->hdr,line),(int64_t) line->pos+1,flt->hdr->samples[i]);
                if ( idx2>=line->n_allele ) error("Incorrect allele index at %s:%"PRId64", sample %s\n", bcf_seqname(flt->hdr,line),(int64_t) line->pos+1,flt->hdr->samples[i]);
                double *vals = tok->values + tok->nval1*i;
                if ( bcf_double_is_missing_or_vector_end(vals[idx1]) || bcf_double_is_missing_or_vector_end(vals[idx2]) )
                {
                    bcf_double_set_missing(rtok->values[i]);
                    continue;
                }
                rtok->values[i] = calc_binom_two_sided(vals[idx1],vals[idx2],0.5);
                if ( rtok->values[i] < 0 )
                {
                    bcf_double_set_missing(rtok->values[i]);
                    continue;
                }
            }
        }
        else
        {
            // the fields given explicitly: binom(AD[:0],AD[:1])
            token_t *tok2 = stack[istack+1];
            if ( tok->nval1!=1 || tok2->nval1!=1 )
                error("Expected one value per binom() argument, found %d and %d at %s:%"PRId64"\n",tok->nval1,tok2->nval1, bcf_seqname(flt->hdr,line),(int64_t) line->pos+1);
            for (i=0; i<rtok->nsamples; i++)
            {
                if ( !rtok->usmpl[i] ) continue;
                double *ptr1 = tok->values + tok->nval1*i;
                double *ptr2 = tok2->values + tok2->nval1*i;
                if ( bcf_double_is_missing_or_vector_end(ptr1[0]) || bcf_double_is_missing_or_vector_end(ptr2[0]) )
                {
                    bcf_double_set_missing(rtok->values[i]);
                    continue;
                }
                rtok->values[i] = calc_binom_two_sided(ptr1[0],ptr2[0],0.5);
                if ( rtok->values[i] < 0 )
                {
                    bcf_double_set_missing(rtok->values[i]);
                    continue;
                }
            }
        }
    }
    else
    {
        // working with an INFO tag
        rtok->nvalues  = 1;
        hts_expand(double, rtok->nvalues, rtok->mvalues, rtok->values);

        double *ptr1 = NULL, *ptr2 = NULL;
        if ( istack+1==nstack )
        {
            // only one tag, expecting two values: binom(INFO/AD)
            if ( tok->nvalues==2 )
            {
                ptr1 = &tok->values[0];
                ptr2 = &tok->values[1];
            }
        }
        else
        {
            // two tags, expecting one value in each: binom(INFO/AD[0],INFO/AD[1])
            token_t *tok2 = stack[istack+1];
            if ( tok->nvalues==1 && tok2->nvalues==1 )
            {
                ptr1 = &tok->values[0];
                ptr2 = &tok2->values[0];
            }
        }
        if ( !ptr1 || !ptr2 || bcf_double_is_missing_or_vector_end(ptr1[0]) || bcf_double_is_missing_or_vector_end(ptr2[0]) )
            bcf_double_set_missing(rtok->values[0]);
        else
        {
            rtok->values[0] = calc_binom_two_sided(ptr1[0],ptr2[0],0.5);
            if ( rtok->values[0] < 0 )
                bcf_double_set_missing(rtok->values[0]);
        }
    }
    return rtok->nargs;
}
static int func_phred(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    token_t *tok = stack[nstack - 1];
    if ( tok->is_str ) error("PHRED() can be applied only on numeric values\n");

    rtok->nsamples = tok->nsamples;
    rtok->nval1 = tok->nval1;
    memcpy(rtok->pass_samples, tok->pass_samples, rtok->nsamples*sizeof(*rtok->pass_samples));
    assert(tok->usmpl);
    if ( !rtok->usmpl )
    {
        rtok->usmpl = (uint8_t*) malloc(tok->nsamples*sizeof(*rtok->usmpl));
        memcpy(rtok->usmpl, tok->usmpl, tok->nsamples*sizeof(*rtok->usmpl));
    }
    rtok->nvalues = tok->nvalues;
    if ( !tok->nvalues ) return 1;

    hts_expand(double, rtok->nvalues, rtok->mvalues, rtok->values);
    int i,j,k = 0;
    if ( tok->usmpl )
    {
        for (i=0; i<tok->nsamples; i++)
        {
            if ( !tok->usmpl[i] ) { k+= tok->nval1; continue; }
            for (j=0; j<tok->nval1; k++,j++)
            {
                if ( bcf_double_is_missing_or_vector_end(tok->values[k]) ) bcf_double_set_missing(rtok->values[k]);
                else rtok->values[k] = -4.34294481903*log(tok->values[k]);
            }
        }
    }
    else
    {
        for (i=0; i<tok->nvalues; i++)
        {
            if ( bcf_double_is_missing_or_vector_end(tok->values[i]) ) bcf_double_set_missing(rtok->values[i]);
            else rtok->values[i] = -4.34294481903*log(tok->values[i]);
        }
    }
    return 1;
}
inline static void tok_init_values(token_t *atok, token_t *btok, token_t *rtok)
{
    token_t *tok;
    if ( (atok->nsamples || btok->nsamples) && (!atok->nsamples || !btok->nsamples) )
        tok = atok->nsamples ? atok : btok;
    else
        tok = atok->nvalues > btok->nvalues ? atok : btok;
    rtok->nvalues = tok->nvalues;
    rtok->nval1   = tok->nval1;
    hts_expand(double, rtok->nvalues, rtok->mvalues, rtok->values);
}
inline static void tok_init_samples(token_t *atok, token_t *btok, token_t *rtok)
{
    if ( (atok->nsamples || btok->nsamples) && !rtok->nsamples )
    {
        rtok->nsamples = atok->nsamples ? atok->nsamples : btok->nsamples;
        rtok->usmpl = (uint8_t*) calloc(rtok->nsamples,1);
        int i;
        for (i=0; i<atok->nsamples; i++) rtok->usmpl[i] |= atok->usmpl[i];
        for (i=0; i<btok->nsamples; i++) rtok->usmpl[i] |= btok->usmpl[i];
    }
    if (rtok->nsamples)
        memset(rtok->pass_samples, 0, rtok->nsamples);
}

#define VECTOR_ARITHMETICS(atok,btok,_rtok,AOP,TYPE) \
{ \
    token_t *rtok = _rtok; \
    int i, has_values = 0; \
    if ( atok->nvalues && btok->nvalues ) \
    { \
        tok_init_values(atok, btok, rtok); \
        tok_init_samples(atok, btok, rtok); \
        if ( !atok->nsamples && !btok->nsamples ) \
        { \
            if ( atok->nvalues!=btok->nvalues && atok->nvalues!=1 && btok->nvalues!=1 ) \
                error("Cannot run numeric operator in -i/-e filtering on vectors of different lengths: %d vs %d\n",atok->nvalues,btok->nvalues); \
            int ir,ia = 0, ib = 0; \
            for (ir=0; ir<rtok->nvalues; ir++) \
            { \
                if ( atok->nvalues > 1 ) ia = ir; \
                if ( btok->nvalues > 1 ) ib = ir; \
                if ( bcf_double_is_missing_or_vector_end(atok->values[ia]) || bcf_double_is_missing_or_vector_end(btok->values[ib]) ) \
                { \
                    bcf_double_set_missing(rtok->values[ir]); \
                    continue; \
                } \
                has_values = 1; \
                rtok->values[ir] = TYPE atok->values[ia] AOP TYPE btok->values[ib]; \
            } \
        } \
        else if ( atok->nsamples && btok->nsamples ) \
        { \
            assert( atok->nsamples==btok->nsamples ); \
            if ( atok->nval1!=btok->nval1 && atok->nval1!=1 && btok->nval1!=1 ) \
                error("Cannot run numeric operator in -i/-e filtering on vectors of different lengths: %d vs %d\n",atok->nval1,btok->nval1); \
            for (i=0; i<rtok->nsamples; i++) \
            { \
                double *rval = rtok->values + i*rtok->nval1; \
                double *aval = atok->values + i*atok->nval1; \
                double *bval = btok->values + i*btok->nval1; \
                int ir,ia = 0, ib = 0; \
                for (ir=0; ir<rtok->nval1; ir++) \
                { \
                    if ( atok->nval1 > 1 ) ia = ir; \
                    if ( btok->nval1 > 1 ) ib = ir; \
                    if ( bcf_double_is_missing_or_vector_end(aval[ia]) || bcf_double_is_missing_or_vector_end(bval[ib]) ) \
                    { \
                        bcf_double_set_missing(rval[ir]); \
                        continue; \
                    } \
                    has_values = 1; \
                    rval[ir] = TYPE aval[ia] AOP TYPE bval[ib]; \
                } \
            } \
        } \
        else if ( atok->nsamples ) \
        { \
            assert( btok->nvalues==1 ); \
            if ( !bcf_double_is_missing_or_vector_end(btok->values[0]) ) \
            { \
                for (i=0; i<atok->nvalues; i++) \
                { \
                    if ( bcf_double_is_missing_or_vector_end(atok->values[i]) ) \
                    { \
                        bcf_double_set_missing(rtok->values[i]); \
                        continue; \
                    } \
                    has_values = 1; \
                    rtok->values[i] = TYPE atok->values[i] AOP TYPE btok->values[0]; \
                } \
            } \
        } \
        else \
        { \
            assert( atok->nvalues==1 ); \
            if ( !bcf_double_is_missing_or_vector_end(atok->values[0]) ) \
            { \
                for (i=0; i<btok->nvalues; i++) \
                { \
                    if ( bcf_double_is_missing_or_vector_end(btok->values[i]) ) \
                    { \
                        bcf_double_set_missing(rtok->values[i]); \
                        continue; \
                    } \
                    has_values = 1; \
                    rtok->values[i] = TYPE atok->values[0] AOP TYPE btok->values[i]; \
                } \
            } \
        } \
    } \
    if ( !has_values ) rtok->nvalues = 0; \
}

static int vector_logic_or(filter_t *filter, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    if ( nstack < 2 ) error("Error occurred while processing the filter \"%s\"\n", filter->str);

    token_t *atok = stack[nstack-2];
    token_t *btok = stack[nstack-1];
    tok_init_samples(atok, btok, rtok);

    if ( !atok->pass_site && !btok->pass_site ) return 2;

    rtok->pass_site = 1;
    if ( !atok->nsamples && !btok->nsamples ) return 2;

    int i;
    if ( rtok->tok_type==TOK_OR_VEC )   // ||, select all samples if one is true
    {
        if ( (!atok->nsamples && !atok->pass_site) || (!btok->nsamples && !btok->pass_site) )
        {
            // These two conditions are to ensure the following does not set all samples
            // at sites with QUAL<=30:
            //      QUAL>30 || FMT/GQ>30

            token_t *tok = atok->nsamples ? atok : btok;
            for (i=0; i<rtok->nsamples; i++)
            {
                if ( !rtok->usmpl[i] ) continue;
                rtok->pass_samples[i] = tok->pass_samples[i];
            }
        }
        else
        {
            for (i=0; i<rtok->nsamples; i++)
            {
                if ( !rtok->usmpl[i] ) continue;
                rtok->pass_samples[i] = 1;
            }
        }
        return 2;
    }

    //  |, only select samples which are actually true

    if ( !atok->nsamples || !btok->nsamples )
    {
        token_t *tok = atok->nsamples ? atok : btok;
        for (i=0; i<rtok->nsamples; i++)
        {
            if ( !rtok->usmpl[i] ) continue;
            rtok->pass_samples[i] = tok->pass_samples[i];
        }
        return 2;
    }

    assert( atok->nsamples==btok->nsamples );

    for (i=0; i<rtok->nsamples; i++)
    {
        if ( !rtok->usmpl[i] ) continue;
        rtok->pass_samples[i] = atok->pass_samples[i] | btok->pass_samples[i];
    }
    return 2;
}
static int vector_logic_and(filter_t *filter, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
    if ( nstack < 2 ) error("Error occurred while processing the filter \"%s\". (nstack=%d)\n", filter->str,nstack);

    token_t *atok = stack[nstack-2];
    token_t *btok = stack[nstack-1];
    tok_init_samples(atok, btok, rtok);

    if ( !atok->pass_site || !btok->pass_site ) return 2;
    if ( !atok->nsamples && !btok->nsamples ) { rtok->pass_site = 1; return 2; }

    int i;
    if ( !atok->nsamples || !btok->nsamples )
    {
        token_t *tok = atok->nsamples ? atok : btok;
        for (i=0; i<rtok->nsamples; i++)
        {
            if ( !rtok->usmpl[i] ) continue;
            rtok->pass_samples[i] = tok->pass_samples[i];
        }
        rtok->pass_site = 1;
        return 2;
    }

    assert( atok->nsamples==btok->nsamples );
    if ( rtok->tok_type==TOK_AND_VEC )  // &&, can be true in different samples
    {
        for (i=0; i<rtok->nsamples; i++)
        {
            if ( !rtok->usmpl[i] ) continue;
            rtok->pass_samples[i] = atok->pass_samples[i] | btok->pass_samples[i];
        }
        rtok->pass_site = 1;
    }
    else    // &, must be true within one sample
    {
        for (i=0; i<rtok->nsamples; i++)
        {
            if ( !rtok->usmpl[i] ) continue;
            rtok->pass_samples[i] = atok->pass_samples[i] & btok->pass_samples[i];
            if ( rtok->pass_samples[i] ) rtok->pass_site = 1;
        }
    }
    return 2;
}

// A note about comparisons:
// When setting value by determining index from the genotype, we face the problem
// of how to interpret truncating arrays. Say we have TAG defined as Number=. and
//      GT:TAG   1/1:0,1,2  0/0:0
// Then when querying we expect the following expression to evaluate for the second
// sample as
//      -i 'TAG[1:1]="."'  .. true
//      -i 'TAG[1:GT]="."' .. false
// The problem is that the implementation truncates the number of fields, filling
// usually fewer than the original number of per-sample values. This is fixed by
// adding an exception that makes the code aware of this: the GT indexing can be
// recognised by having tok->idx==-3
#define CMP_VECTORS(atok,btok,_rtok,CMP_OP,missing_logic) \
{ \
    token_t *rtok = _rtok; \
    int i, j, k; \
    tok_init_samples(atok, btok, rtok); \
    if ( !atok->nsamples && !btok->nsamples ) \
    { \
        if ( !atok->nvalues && !btok->nvalues ) { rtok->pass_site = missing_logic[2]; } \
        else if ( !atok->nvalues || !btok->nvalues ) \
        { \
            token_t *tok = atok->nvalues ? atok : btok; \
            for (j=0; j<tok->nvalues; j++) \
            { \
                if ( bcf_double_is_missing_or_vector_end(tok->values[j]) ) \
                { \
                    if ( missing_logic[2] ) { rtok->pass_site = 1; break; } \
                } \
                else if ( missing_logic[1] ) {  rtok->pass_site = 1; break; } \
            } \
        } \
        else \
        { \
            for (i=0; i<atok->nvalues; i++) \
            { \
                int amiss = bcf_double_is_missing_or_vector_end(atok->values[i]) ? 1 : 0; \
                for (j=0; j<btok->nvalues; j++) \
                { \
                    int nmiss = amiss + (bcf_double_is_missing_or_vector_end(btok->values[j]) ? 1 : 0); \
                    if ( nmiss ) \
                    { \
                        if ( missing_logic[nmiss] ) { rtok->pass_site = 1; i = atok->nvalues; break; } \
                    } \
                    else if ( atok->values[i] > 16777216 || btok->values[j] > 16777216 ) /* Ugly, see #871 */ \
                    { \
                        if ( atok->values[i] CMP_OP btok->values[j] ) { rtok->pass_site = 1; i = atok->nvalues; break; } \
                    } \
                    else if ( (float)atok->values[i] CMP_OP (float)btok->values[j] ) { rtok->pass_site = 1; i = atok->nvalues; break; } \
                } \
            } \
        } \
    } \
    else \
    { \
        if ( !atok->nvalues && !btok->nvalues ) \
        { \
            if ( missing_logic[2] ) \
            { \
                for (i=0; i<rtok->nsamples; i++) \
                    if ( rtok->usmpl[i] ) { rtok->pass_samples[i] = missing_logic[2]; rtok->pass_site = 1; } \
            } \
        } \
        else if ( !atok->nvalues || !btok->nvalues ) \
        { \
            token_t *tok = atok->nvalues ? atok : btok; \
            if ( !tok->nsamples ) \
            { \
                int miss = 0; \
                for (j=0; j<tok->nvalues; j++) \
                    miss |= bcf_double_is_missing_or_vector_end(tok->values[j]) ? 1 : 0; \
                if ( missing_logic[++miss] ) \
                { \
                    for (i=0; i<rtok->nsamples; i++) \
                        if ( rtok->usmpl[i] ) { rtok->pass_samples[i] = missing_logic[miss]; rtok->pass_site = 1; } \
                } \
            } \
            else \
                for (i=0; i<tok->nsamples; i++) \
                { \
                    if ( !rtok->usmpl[i] ) continue; \
                    double *ptr = tok->values + i*tok->nval1; \
                    int miss = 0; \
                    for (j=0; j<tok->nval1; j++) \
                        miss |= bcf_double_is_missing_or_vector_end(ptr[j]) ? 1 : 0; \
                    if ( missing_logic[++miss] ) { rtok->pass_samples[i] = missing_logic[miss]; rtok->pass_site = 1; } \
                } \
        } \
        else if ( atok->nsamples && btok->nsamples ) \
        { \
            if ( atok->nval1!=btok->nval1 ) error("Incompatible number of per-sample values in comparison: %d vs %d\n",atok->nval1,btok->nval1); \
            if ( atok->nsamples!=btok->nsamples ) error("Incompatible number samples in comparison: %d vs %d\n",atok->nsamples,btok->nsamples); \
            for (i=0; i<atok->nsamples; i++) \
            { \
                if ( !atok->usmpl[i] || !btok->usmpl[i] ) { rtok->usmpl[i] = 0; continue; } \
                double *aptr = atok->values + i*atok->nval1; \
                double *bptr = btok->values + i*btok->nval1; \
                for (j=0; j<atok->nval1; j++) \
                { \
                    if ( atok->idx==-3 && bcf_double_is_vector_end(aptr[j]) ) break; /* explained above */ \
                    if ( btok->idx==-3 && bcf_double_is_vector_end(bptr[j]) ) break; /* explained above */ \
                    int nmiss = bcf_double_is_missing_or_vector_end(aptr[j]) ? 1 : 0; \
                    if ( nmiss && !missing_logic[0] ) continue; /* any is missing => result is false */ \
                    nmiss += (bcf_double_is_missing_or_vector_end(bptr[j]) ? 1 : 0); \
                    if ( nmiss ) \
                    { \
                        if ( missing_logic[nmiss] ) { rtok->pass_samples[i] = 1; rtok->pass_site = 1; break; } \
                    } \
                    else if ( aptr[j] > 16777216 || bptr[j] > 16777216 ) /* Ugly, see #871 */ \
                    { \
                        if ( aptr[j] CMP_OP bptr[j] ) { rtok->pass_samples[i] = 1; rtok->pass_site = 1; break; } \
                    } \
                    else if ( (float)aptr[j] CMP_OP (float)bptr[j] ) { rtok->pass_samples[i] = 1; rtok->pass_site = 1; break; } \
                } \
            } \
        } \
        else if ( atok->nsamples )\
        { \
            for (i=0; i<atok->nsamples; i++) \
            { \
                if ( !rtok->usmpl[i] ) continue; \
                double *aptr = atok->values + i*atok->nval1; \
                double *bptr = btok->values; \
                for (j=0; j<atok->nval1; j++) \
                { \
                    if ( atok->idx==-3 && bcf_double_is_vector_end(aptr[j]) ) break; /* explained above */ \
                    int miss = bcf_double_is_missing_or_vector_end(aptr[j]) ? 1 : 0; \
                    if ( miss && !missing_logic[0] ) continue; /* any is missing => result is false */ \
                    for (k=0; k<btok->nvalues; k++) \
                    { \
                        int nmiss = miss + (bcf_double_is_missing_or_vector_end(bptr[k]) ? 1 : 0); \
                        if ( nmiss ) \
                        { \
                            if ( missing_logic[nmiss] ) { rtok->pass_samples[i] = 1; rtok->pass_site = 1; j = atok->nval1; break; } \
                        } \
                        else if ( aptr[j] > 16777216 || bptr[k] > 16777216 ) /* Ugly, see #871 */ \
                        { \
                            if ( aptr[j] CMP_OP bptr[k] ) { rtok->pass_samples[i] = 1; rtok->pass_site = 1; j = atok->nval1; break; } \
                        } \
                        else if ( (float)aptr[j] CMP_OP (float)bptr[k] ) { rtok->pass_samples[i] = 1; rtok->pass_site = 1; j = atok->nval1; break; } \
                    } \
                } \
            } \
        } \
        else /* btok->nsamples */ \
        { \
            for (i=0; i<btok->nsamples; i++) \
            { \
                if ( !rtok->usmpl[i] ) continue; \
                double *aptr = atok->values; \
                double *bptr = btok->values + i*btok->nval1; \
                for (j=0; j<btok->nval1; j++) \
                { \
                    if ( atok->idx==-3 && bcf_double_is_vector_end(bptr[j]) ) break; /* explained above */ \
                    int miss = bcf_double_is_missing_or_vector_end(bptr[j]) ? 1 : 0; \
                    if ( miss && !missing_logic[0] ) continue; /* any is missing => result is false */ \
                    for (k=0; k<atok->nvalues; k++) \
                    { \
                        int nmiss = miss + (bcf_double_is_missing_or_vector_end(aptr[k]) ? 1 : 0); \
                        if ( nmiss ) \
                        { \
                            if ( missing_logic[nmiss] ) { rtok->pass_samples[i] = 1; rtok->pass_site = 1; j = btok->nval1; break; } \
                        } \
                        else if ( bptr[j] > 16777216 || aptr[k] > 16777216 ) /* Ugly, see #871 */ \
                        { \
                            if ( aptr[k] CMP_OP bptr[j] ) { rtok->pass_samples[i] = 1; rtok->pass_site = 1; j = btok->nval1; break; } \
                        } \
                        else if ( (float)aptr[k] CMP_OP (float)bptr[j] ) { rtok->pass_samples[i] = 1; rtok->pass_site = 1; j = btok->nval1; break; } \
                    } \
                } \
            } \
        } \
    } \
}
static int _regex_vector_strings(regex_t *regex, char *str, size_t len, int logic, int *missing_logic)
{
    char *end = str + len;
    while ( str < end && *str )
    {
        char *mid = str;
        while ( mid < end && *mid && *mid!=',' ) mid++;
        int miss = mid - str == 1 && str[0]=='.' ? 1 : 0;
        if ( miss && missing_logic[miss] ) return 1;
        char tmp = *mid; *mid = 0;
        int match = regexec(regex, str, 0,NULL,0) ? 0 : 1;
        *mid = tmp;
        if ( logic==TOK_NLIKE ) match = match ? 0 : 1;
        if ( match ) return 1;
        if ( !*mid ) break;
        str = mid + 1;
    }
    return 0;
}
static inline int _has_missing_string(char *beg)
{
    while ( *beg )
    {
        char *end = beg;
        while ( *end && *end!=',' ) end++;
        if ( end-beg==1 && beg[0]=='.' ) return 1;
        if ( !*end ) break;
        beg = end + 1;
    }
    return 0;
}

// Compare two strings with multiple fields, for example "A,B,.,C"=="X,Y,A".
// Returns 1 if any field matches, otherwise returns 0
static inline int _match_vector_strings(char *abeg, size_t alen, char *bstr, size_t blen, int logic, int *missing_logic)
{
    char *aend = abeg + alen;
    char *bend = bstr + blen;
    while ( abeg < aend && *abeg )
    {
        char *amid = abeg;
        while ( amid < aend && *amid && *amid!=',' ) amid++;
        int miss = amid - abeg == 1 && abeg[0]=='.' ? 1 : 0;
        char *bbeg = bstr;
        while ( bbeg < bend && *bbeg )
        {
            char *bmid = bbeg;
            while ( bmid < bend && *bmid && *bmid!=',' ) bmid++;
            int nmiss = miss + (bmid - bbeg == 1 && bbeg[0]=='.' ? 1 : 0);
            if ( nmiss )
            {
                if ( missing_logic[nmiss] ) return 1;
            }
            else
            {
                int match = amid-abeg==bmid-bbeg && !strncmp(abeg,bbeg,amid-abeg) ? 1 : 0;
                if ( logic==TOK_NE ) match = match==1 ? 0 : 1;
                if ( match ) return 1;
            }
            if ( !*bmid ) break;
            bbeg = bmid + 1;
        }
        if ( !*amid ) break;
        abeg = amid + 1;
    }
    return 0;
}
static void cmp_vector_strings(token_t *atok, token_t *btok, token_t *rtok)
{
    tok_init_samples(atok, btok, rtok);

    int i, logic = rtok->tok_type;     // TOK_EQ, TOK_NE, TOK_LIKE, TOK_NLIKE
    regex_t *regex = atok->regex ? atok->regex : (btok->regex ? btok->regex : NULL);

    assert( atok->nvalues==atok->str_value.l && btok->nvalues==btok->str_value.l );
    assert( !atok->nsamples || !btok->nsamples );
    assert( (!regex && (logic==TOK_EQ || logic==TOK_NE)) || (regex && (logic==TOK_LIKE || logic==TOK_NLIKE)) );

    int missing_logic[] = {0,0,0};
    if ( logic==TOK_EQ || logic==TOK_LIKE ) missing_logic[0] = missing_logic[2] = 1;
    else if ( logic==TOK_NE || logic==TOK_NLIKE ) missing_logic[0] = missing_logic[1] = 1;

    if ( !atok->nsamples && !btok->nsamples )
    {
        if ( !atok->nvalues && !btok->nvalues ) { rtok->pass_site = missing_logic[2]; return; }
        if ( !atok->nvalues || !btok->nvalues )
        {
            int miss = _has_missing_string(atok->nvalues ? atok->str_value.s : btok->str_value.s);
            if ( missing_logic[miss+1] ) rtok->pass_site = 1;
            return;
        }
        if ( !regex )
            rtok->pass_site = _match_vector_strings(atok->str_value.s, atok->str_value.l, btok->str_value.s, btok->str_value.l, logic, missing_logic);
        else
        {
            token_t *tok = atok->regex ? btok : atok;
            rtok->pass_site = _regex_vector_strings(regex, tok->str_value.s, tok->str_value.l, logic, missing_logic);
        }
        return;
    }

    // The case of (!atok->nsamples || !btok->nsamples)

    if ( !atok->nvalues && !btok->nvalues )
    {
        if ( missing_logic[2] )
        {
            for (i=0; i<rtok->nsamples; i++)
                if ( rtok->usmpl[i] ) { rtok->pass_samples[i] = missing_logic[2]; rtok->pass_site = 1; }
        }
        return;
    }
    if ( !atok->nvalues || !btok->nvalues )
    {
        token_t *tok = atok->nvalues ? atok : btok;
        if ( !tok->nsamples )
        {
            int miss = _has_missing_string(tok->str_value.s);
            if ( !missing_logic[miss+1] ) return;
            for (i=0; i<rtok->nsamples; i++)
                if ( rtok->usmpl[i] ) { rtok->pass_samples[i] = 1; rtok->pass_site = 1; }
        }
        else
            for (i=0; i<tok->nsamples; i++)
            {
                if ( !rtok->usmpl[i] ) continue;
                int miss = _has_missing_string(tok->str_value.s + i*tok->nval1);
                if ( missing_logic[miss+1] ) { rtok->pass_samples[i] = 1; rtok->pass_site = 1; }
            }
        return;
    }

    // The case of (!atok->nsamples || !btok->nsamples) && (atok->nvalues && btok->nvalues)
    token_t *xtok = atok->nsamples ? atok : btok;
    token_t *ytok = atok->nsamples ? btok : atok;
    assert( regex==ytok->regex );
    for (i=0; i<xtok->nsamples; i++)
    {
        if ( !rtok->usmpl[i] ) continue;
        int match;
        if ( regex )
            match = _regex_vector_strings(regex, xtok->str_value.s + i*xtok->nval1, xtok->nval1, logic, missing_logic);
        else
            match = _match_vector_strings(xtok->str_value.s + i*xtok->nval1, xtok->nval1, ytok->str_value.s, ytok->str_value.l, logic, missing_logic);
        if ( match ) { rtok->pass_samples[i] = 1; rtok->pass_site = 1; }
    }
}

static int parse_idxs(char *tag_idx, int **idxs, int *nidxs, int *idx)
{
    // TAG[], TAG[*] .. any field; sets idx=-2, idxs[0]=-1
    if ( *tag_idx==0 || !strcmp("*", tag_idx) )
    {
        *idxs = (int*) malloc(sizeof(int));
        (*idxs)[0] = -1;
        *nidxs     = 1;
        *idx       = -2;
        return 0;
    }
    if ( !strcmp("GT", tag_idx) )
    {
        *idxs = (int*) malloc(sizeof(int));
        (*idxs)[0] = -1;
        *nidxs     = 1;
        *idx = -3;
        return 0;
    }

    // TAG[integer] .. one field; idx positive
    char *end, *beg = tag_idx;
    *idx = strtol(tag_idx, &end, 10);
    if ( *idx >= 0 && *end==0 ) return 0;

    // TAG[0,1] or TAG[0-2] or [1-] etc; idx=-2, idxs[...]=0,0,1,1,..
    int i, ibeg = -1;
    while ( *beg )
    {
        int num = strtol(beg, &end, 10);
        if ( end[0]==',' ) beg = end + 1;
        else if ( end[0]==0 ) beg = end;
        else if ( end[0]=='-' ) { beg = end + 1; ibeg = num; continue; }
        else return -1;
        if ( num >= *nidxs )
        {
            *idxs = (int*) realloc(*idxs, sizeof(int)*(num+1));
            memset(*idxs + *nidxs, 0, sizeof(int)*(num - *nidxs + 1));
            *nidxs = num + 1;
        }
        if ( ibeg>=0 )
        {
            for (i=ibeg; i<=num; i++) (*idxs)[i] = 1;
            ibeg = -1;
        }
        (*idxs)[num] = 1;
    }
    if ( ibeg >=0 )
    {
        if ( ibeg >= *nidxs )
        {
            *idxs = (int*) realloc(*idxs, sizeof(int)*(ibeg+1));
            memset(*idxs + *nidxs, 0, sizeof(int)*(ibeg - *nidxs + 1));
            *nidxs = ibeg + 1;
        }
        (*idxs)[ibeg] = -1;
    }
    *idx = -2;
    return 0;
}

static void parse_tag_idx(bcf_hdr_t *hdr, int is_fmt, char *tag, char *tag_idx, token_t *tok)   // tag_idx points just after "TAG["
{
    int i, len = strlen(tag_idx);
    if ( tag_idx[len-1] == ']' ) tag_idx[len-1] = 0;
    char *ori = strdup(tag_idx);

    assert( !tok->idxs && !tok->usmpl );
    int *idxs1 = NULL, nidxs1 = 0, idx1 = 0;
    int *idxs2 = NULL, nidxs2 = 0, idx2 = 0;

    int set_samples = 0;
    char *colon = strrchr(tag_idx, ':');
    if ( tag_idx[0]=='@' )     // file list with sample names
    {
        if ( !is_fmt ) error("Could not parse \"%s\". (Not a FORMAT tag yet a sample list provided.)\n", ori);
        char *fname = expand_path(tag_idx+1);
#ifdef _WIN32
        if (fname && strlen(fname) > 2 && fname[1] == ':') // Deal with Windows paths, such as 'C:\..'
            colon = strrchr(fname+2, ':');
#endif
        int nsmpl;
        char **list = hts_readlist(fname, 1, &nsmpl);
        if ( !list && colon )
        {
            if ( parse_idxs(colon+1, &idxs2, &nidxs2, &idx2) != 0 ) error("Could not parse the index: %s\n", ori);
            tok->idxs  = idxs2;
            tok->nidxs = nidxs2;
            tok->idx   = idx2;
            colon = strrchr(fname, ':');
            *colon = 0;
            list = hts_readlist(fname, 1, &nsmpl);
        }
        if ( !list ) error("Could not read: %s\n", fname);
        free(fname);
        tok->nsamples = bcf_hdr_nsamples(hdr);
        tok->usmpl = (uint8_t*) calloc(tok->nsamples,1);
        for (i=0; i<nsmpl; i++)
        {
            int ismpl = bcf_hdr_id2int(hdr,BCF_DT_SAMPLE,list[i]);
            if ( ismpl<0 ) error("No such sample in the VCF: \"%s\"\n", list[i]);
            free(list[i]);
            tok->usmpl[ismpl] = 1;
        }
        free(list);
        if ( !colon )
        {
            tok->idxs = (int*) malloc(sizeof(int));
            tok->idxs[0] = -1;
            tok->nidxs   = 1;
            tok->idx     = -2;
        }
    }
    else if ( colon )
    {
        if ( !is_fmt ) error("Could not parse the index \"%s\". (Not a FORMAT tag yet sample index implied.)\n", ori);
        *colon = 0;
        if ( parse_idxs(tag_idx, &idxs1, &nidxs1, &idx1) != 0 ) error("Could not parse the index: %s\n", ori);
        if ( parse_idxs(colon+1, &idxs2, &nidxs2, &idx2) != 0 ) error("Could not parse the index: %s\n", ori);
        tok->idxs  = idxs2;
        tok->nidxs = nidxs2;
        tok->idx   = idx2;
        set_samples = 1;
    }
    else
    {
        if ( parse_idxs(tag_idx, &idxs1, &nidxs1, &idx1) != 0 ) error("Could not parse the index: %s\n", ori);
        if ( is_fmt )
        {
            if ( nidxs1==1 && idxs1[0]==-1 )
            {
                tok->idxs = (int*) malloc(sizeof(int));
                tok->idxs[0] = -1;
                tok->nidxs   = 1;
                tok->idx     = idx1;
            }
            else if ( bcf_hdr_id2number(hdr,BCF_HL_FMT,tok->hdr_id)!=1 )
                error("The FORMAT tag %s can have multiple subfields, run as %s[sample:subfield]\n", tag,tag);
            else
                tok->idx = 0;
            set_samples = 1;
        }
        else
        {
            tok->idxs  = idxs1;
            tok->nidxs = nidxs1;
            tok->idx   = idx1;
        }
    }

    if ( set_samples )
    {
        tok->nsamples = bcf_hdr_nsamples(hdr);
        tok->usmpl = (uint8_t*) calloc(tok->nsamples,1);
        if ( idx1>=0 )
        {
            if ( idx1 >= bcf_hdr_nsamples(hdr) ) error("The sample index is too large: %s\n", ori);
            tok->usmpl[idx1] = 1;
        }
        else if ( idx1==-2 || idx1==-3 )
        {
            for (i=0; i<nidxs1; i++)
            {
                if ( idxs1[i]==0 ) continue;
                if ( idxs1[i]==-1 ) break;
                if ( i >= bcf_hdr_nsamples(hdr) ) error("The sample index is too large: %s\n", ori);
                tok->usmpl[i] = 1;
            }
            if ( nidxs1 && idxs1[nidxs1-1]==-1 )    // open range, such as "7-"
            {
                for (; i<tok->nsamples; i++) tok->usmpl[i] = 1;
            }
        }
        else error("todo: %s:%d .. %d\n", __FILE__,__LINE__, idx2);
        free(idxs1);
    }
    free(ori);

    if ( tok->nidxs && tok->idxs[tok->nidxs-1]!=-1 )
    {
        for (i=0; i<tok->nidxs; i++) if ( tok->idxs[i] ) tok->nuidxs++;
    }
}
static int max_ac_an_unpack(bcf_hdr_t *hdr)
{
    int hdr_id = bcf_hdr_id2int(hdr,BCF_DT_ID,"AC");
    if ( hdr_id<0 ) return BCF_UN_FMT;
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,hdr_id) ) return BCF_UN_FMT;

    hdr_id = bcf_hdr_id2int(hdr,BCF_DT_ID,"AN");
    if ( hdr_id<0 ) return BCF_UN_FMT;
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,hdr_id) ) return BCF_UN_FMT;

    return BCF_UN_INFO;
}
static int filters_init1_ext(filter_t *filter, char *str, int len, token_t *tok)
{
    tok->hl_type  = -1;
    tok->ht_type  = -1;
    tok->tok_type  = TOK_VAL;
    tok->hdr_id    = -1;
    tok->pass_site = -1;
    tok->idx       = 0;
    tok->iext = ++filter->n_ext;
    filter->ext = realloc(filter->ext,sizeof(*filter->ext)*filter->n_ext);
    if ( !strncasecmp(str,"{str}",len) ) { tok->ht_type = BCF_HT_STR; tok->is_str = 1; }
    else if ( !strncasecmp(str,"{int}",len) ) tok->ht_type = BCF_HT_INT;
    else if ( !strncasecmp(str,"{float}",len) ) tok->ht_type = BCF_HT_REAL;
    filter->ext[filter->n_ext-1] = tok->ht_type;
    return 0;
}
static int filters_init1(filter_t *filter, char *str, int len, token_t *tok)
{
    tok->ht_type  = -1;
    tok->hl_type  = -1;
    tok->tok_type  = TOK_VAL;
    tok->hdr_id    = -1;
    tok->pass_site = -1;
    tok->idx       = 0;

    // is this a string constant?
    if ( str[0]=='"' || str[0]=='\'' )
    {
        int quote = str[0];
        if ( str[len-1] != quote ) error("TODO: [%s]\n", filter->str);
        tok->key = (char*) calloc(len-1,sizeof(char));
        hts_expand(double,1,tok->mvalues,tok->values);
        tok->values[0] = len-2;
        memcpy(tok->key,str+1,len-2);
        tok->key[len-2] = 0;
        tok->is_str = 1;
        tok->ht_type = BCF_HT_STR;
        tok->nvalues = len-2;
        if ( !strcmp(".",tok->key) ) tok->is_missing = 1;
        return 0;
    }

    // is it a file?
    if ( str[0]=='@' )
    {
        tok->tag = (char*) calloc(len+1,sizeof(char));
        memcpy(tok->tag,str,len);
        tok->tag[len] = 0;
        char *fname = expand_path(tok->tag+1);
        int i, n;
        char **list = hts_readlist(fname, 1, &n);
        if ( !list ) error("Could not read: %s\n", fname);
        free(fname);
        tok->hash = khash_str2int_init();
        for (i=0; i<n; i++)
        {
            char *se = list[i];
            while ( *se && !isspace(*se) ) se++;
            *se = 0;
            if ( !khash_str2int_has_key(tok->hash,list[i]) )
                khash_str2int_inc(tok->hash,list[i]);
            else
                free(list[i]);
        }
        free(list);
        return 0;
    }

    int is_fmt = -1;
    if ( !strncasecmp(str,"FMT/",4) ) { str += 4; len -= 4; is_fmt = 1; }
    else if ( !strncasecmp(str,"FORMAT/",7) ) { str += 7; len -= 7; is_fmt = 1; }
    else
    {
        if ( !strncasecmp(str,"INFO/",5) ) { is_fmt = 0; str += 5; len -= 5; }
        else if ( !strncasecmp(str,"QUAL",len) || !strncmp(str,"%QUAL",len) /* for backward compatibility */ )
        {
            tok->setter = filters_set_qual;
            tok->tag = strdup("QUAL");
            tok->ht_type = BCF_HT_REAL;
            filter_add_used_tag(filter,NULL,tok->tag);
            return 0;
        }
        else if ( !strncasecmp(str,"TYPE",len) || !strncmp(str,"%TYPE",len) /* for backward compatibility */ )
        {
            tok->setter = filters_set_type;
            tok->tag = strdup("TYPE");
            tok->ht_type = BCF_HT_STR;
            return 0;
        }
        else if ( !strncasecmp(str,"FILTER",len) || !strncmp(str,"%FILTER",len) /* for backward compatibility */ )
        {
            tok->comparator = filters_cmp_filter;
            tok->tag = strdup("FILTER");
            filter->max_unpack |= BCF_UN_FLT;
            tok->hl_type = BCF_HL_FLT;
            tok->ht_type = BCF_HT_STR;
            filter_add_used_tag(filter,NULL,tok->tag);
            return 0;
        }
        else if ( !strncasecmp(str,"ID",len) || !strncasecmp(str,"%ID",len) /* for backward compatibility */ )
        {
            tok->comparator = filters_cmp_id;
            tok->tag = strdup("ID");
            tok->ht_type = BCF_HT_STR;
            filter_add_used_tag(filter,NULL,tok->tag);
            return 0;
        }
        else if ( !strncasecmp(str,"CHROM",len) )
        {
            tok->setter = &filters_set_chrom;
            tok->tag = strdup("CHROM");
            tok->ht_type = BCF_HT_STR;
            filter_add_used_tag(filter,NULL,tok->tag);
            return 0;
        }
        else if ( !strncasecmp(str,"POS",len) )
        {
            tok->setter = &filters_set_pos;
            tok->tag = strdup("POS");
            tok->ht_type = BCF_HT_INT;
            filter_add_used_tag(filter,NULL,tok->tag);
            return 0;
        }
        else if ( !strncasecmp(str,"REF",len) )
        {
            tok->setter = &filters_set_ref_string;
            tok->is_str = 1;
            tok->tag = strdup("REF");
            tok->ht_type = BCF_HT_STR;
            filter_add_used_tag(filter,NULL,tok->tag);
            return 0;
        }
        else if ( !strncasecmp(str,"ALT",len) )
        {
            tok->setter = &filters_set_alt_string;
            tok->is_str = 1;
            tok->tag = strdup("ALT");
            tok->ht_type = BCF_HT_STR;
            tok->idxs = (int*) malloc(sizeof(int));
            tok->idxs[0] = -1;
            tok->nidxs   = 1;
            tok->idx     = -2;
            filter_add_used_tag(filter,NULL,tok->tag);
            return 0;
        }
        else if ( !strncasecmp(str,"N_ALT",len) )
        {
            tok->setter = &filters_set_nalt;
            tok->tag = strdup("N_ALT");
            tok->ht_type = BCF_HT_INT;
            return 0;
        }
        else if ( !strncasecmp(str,"N_SAMPLES",len) )
        {
            tok->tok_type = TOK_VAL;
            tok->threshold = bcf_hdr_nsamples(filter->hdr);
            tok->is_constant = 1;
            tok->ht_type = BCF_HT_INT;
            return 0;
        }
        else if ( !strncasecmp(str,"N_MISSING",len) )
        {
            filter->max_unpack |= BCF_UN_FMT;
            tok->setter = &filters_set_nmissing;
            tok->tag = strdup("N_MISSING");
            tok->ht_type = BCF_HT_INT;
            return 0;
        }
        else if ( !strncasecmp(str,"F_MISSING",len) )
        {
            filter->max_unpack |= BCF_UN_FMT;
            tok->setter = &filters_set_nmissing;
            tok->tag = strdup("F_MISSING");
            tok->ht_type = BCF_HT_REAL;
            return 0;
        }
    }

    // does it have array subscript?
    int is_array = 0;
    kstring_t tmp = {0,0,0};
    kputsn(str, len, &tmp);
    if ( tmp.s[tmp.l-1] == ']' )
    {
        int i;
        for (i=0; i<tmp.l; i++)
            if ( tmp.s[i]=='[' ) { tmp.s[i] = 0; is_array = i+1; break; }
    }
    tok->hdr_id = bcf_hdr_id2int(filter->hdr,BCF_DT_ID,tmp.s);
    if ( is_fmt==-1 )
    {
        if ( tok->hdr_id >=0 )
        {
            int is_info = bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_INFO,tok->hdr_id) ? 1 : 0;
            is_fmt = bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_FMT,tok->hdr_id) ? 1 : 0;
            if ( is_info && is_fmt )
                error("Error: ambiguous filtering expression, both INFO/%s and FORMAT/%s are defined in the VCF header.\n" , tmp.s,tmp.s);
        }
        if ( is_fmt==-1 ) is_fmt = 0;
    }
    if ( is_array )
    {
        parse_tag_idx(filter->hdr, is_fmt, tmp.s, tmp.s+is_array, tok);
        if ( tok->idx==-3 && bcf_hdr_id2length(filter->hdr,BCF_HL_FMT,tok->hdr_id)!=BCF_VL_R )
            error("Error: GT subscripts can be used only with Number=R tags\n");
    }
    else if ( is_fmt && !tok->nsamples )
    {
        int i;
        tok->nsamples = bcf_hdr_nsamples(filter->hdr);
        tok->usmpl = (uint8_t*) malloc(tok->nsamples);
        for (i=0; i<tok->nsamples; i++) tok->usmpl[i] = 1;
    }

    tok->hl_type = is_fmt ? BCF_HL_FMT : BCF_HL_INFO;
    if ( is_fmt ) filter->max_unpack |= BCF_UN_FMT;
    if ( tok->hdr_id>=0 )
    {
        if ( is_fmt && !strcmp("GT",tmp.s) )
        {
            tok->setter = &filters_set_genotype_string; tok->is_str = 1;
            tok->ht_type = BCF_HT_STR;
        }
        else if ( is_fmt )
        {
            if ( !bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_FMT,tok->hdr_id) )
                error("No such FORMAT field: %s\n", tmp.s);
            if ( bcf_hdr_id2number(filter->hdr,BCF_HL_FMT,tok->hdr_id)!=1 && !is_array )
            {
                tok->idxs = (int*) malloc(sizeof(int));
                tok->idxs[0] = -1;
                tok->nidxs   = 1;
                tok->idx     = -2;
            }
            switch ( bcf_hdr_id2type(filter->hdr,BCF_HL_FMT,tok->hdr_id) )
            {
                case BCF_HT_INT:  tok->setter = &filters_set_format_int; tok->ht_type = BCF_HT_INT; break;
                case BCF_HT_REAL: tok->setter = &filters_set_format_float; tok->ht_type = BCF_HT_REAL; break;
                case BCF_HT_STR:  tok->setter = &filters_set_format_string; tok->ht_type = BCF_HT_STR; tok->is_str = 1; break;
                default: error("[%s:%d %s] FIXME\n", __FILE__,__LINE__,__FUNCTION__);
            }
        }
        else if ( !bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_INFO,tok->hdr_id) )
            error("No such INFO field: %s\n", tmp.s);
        else
        {
            if ( bcf_hdr_id2type(filter->hdr,BCF_HL_INFO,tok->hdr_id) == BCF_HT_FLAG )
            {
                tok->setter = filters_set_info_flag;
                tok->ht_type = BCF_HT_INT;
            }
            else
            {
                tok->ht_type = bcf_hdr_id2type(filter->hdr,BCF_HL_INFO,tok->hdr_id);
                if ( tok->ht_type == BCF_HT_STR ) tok->is_str = 1;
                if ( bcf_hdr_id2number(filter->hdr,BCF_HL_INFO,tok->hdr_id)==1 )
                    tok->setter = filters_set_info;
                else
                {
                    switch ( bcf_hdr_id2type(filter->hdr,BCF_HL_INFO,tok->hdr_id) )
                    {
                        case BCF_HT_INT:  tok->setter = &filters_set_info_int; break;
                        case BCF_HT_REAL: tok->setter = &filters_set_info_float; break;
                        case BCF_HT_STR:  tok->setter = &filters_set_info_string; tok->is_str = 1; break;
                        default: error("[%s:%d %s] FIXME\n", __FILE__,__LINE__,__FUNCTION__);
                    }
                    if (!is_array)
                    {
                        tok->idx = -2;
                        tok->idxs = (int*) malloc(sizeof(int));
                        tok->idxs[0] = -1;
                        tok->nidxs   = 1;
                    }
                }
            }
            filter->max_unpack |= BCF_UN_INFO;
        }
        tok->tag = strdup(tmp.s);
        if ( tmp.s ) free(tmp.s);
        filter_add_used_tag(filter,is_fmt ? "FORMAT/" : "INFO/",tok->tag);
        return 0;
    }
    else if ( !strcasecmp(tmp.s,"ALT") )
    {
        tok->setter = &filters_set_alt_string;
        tok->is_str = 1;
        tok->ht_type = BCF_HT_STR;
        tok->tag = strdup(tmp.s);
        free(tmp.s);
        filter_add_used_tag(filter,NULL,tok->tag);
        return 0;
    }
    else if ( !strcasecmp(tmp.s,"AN") )
    {
        filter->max_unpack |= BCF_UN_FMT;
        tok->setter = &filters_set_an;
        tok->tag = strdup("AN");
        tok->ht_type = BCF_HT_INT;
        free(tmp.s);
        return 0;
    }
    else if ( !strcasecmp(tmp.s,"AC") )
    {
        filter->max_unpack |= BCF_UN_FMT;
        tok->setter = &filters_set_ac;
        tok->tag = strdup("AC");
        tok->ht_type = BCF_HT_INT;
        free(tmp.s);
        return 0;
    }
    else if ( !strcasecmp(tmp.s,"MAC") )
    {
        filter->max_unpack |= max_ac_an_unpack(filter->hdr);
        tok->setter = &filters_set_mac;
        tok->tag = strdup("MAC");
        tok->ht_type = BCF_HT_INT;
        free(tmp.s);
        return 0;
    }
    else if ( !strcasecmp(tmp.s,"AF") )
    {
        filter->max_unpack |= max_ac_an_unpack(filter->hdr);
        tok->setter = &filters_set_af;
        tok->tag = strdup("AF");
        tok->ht_type = BCF_HT_REAL;
        free(tmp.s);
        return 0;
    }
    else if ( !strcasecmp(tmp.s,"MAF") )
    {
        filter->max_unpack |= max_ac_an_unpack(filter->hdr);
        tok->setter = &filters_set_maf;
        tok->tag = strdup("MAF");
        tok->ht_type = BCF_HT_REAL;
        free(tmp.s);
        return 0;
    }
    else if ( !strcasecmp(tmp.s,"ILEN") || !strcasecmp(tmp.s,"%ILEN") )
    {
        filter->max_unpack |= BCF_UN_STR;
        tok->setter = &filters_set_ilen;
        tok->tag = strdup("ILEN");
        tok->ht_type = BCF_HT_INT;
        free(tmp.s);
        return 0;
    }

    // is it a value? Here we parse as integer/float separately and use strtof
    // rather than strtod, because the more accurate double representation
    // would invalidate floating point comparisons like QUAL=59.2, obtained via
    // htslib/vcf parser.
    // Update: use strtod() and force floats only in comparisons
    char *end;
    tok->threshold = strtol(tmp.s, &end, 10);   // integer?
    if ( end - tmp.s != strlen(tmp.s) )
    {
        errno = 0;
        tok->threshold = strtod(tmp.s, &end);   // float?
        if ( errno!=0 || end!=tmp.s+len )
        {
            if ( filter->exit_on_error )
                error("[%s:%d %s] Error: the tag \"%s\" is not defined in the VCF header\n", __FILE__,__LINE__,__FUNCTION__,tmp.s);
            filter->status |= FILTER_ERR_UNKN_TAGS;
            filter_add_undef_tag(filter,tmp.s);
        }
        tok->ht_type = BCF_HT_REAL;
    }
    else
        tok->ht_type = BCF_HT_INT;
    tok->is_constant = 1;

    if ( tmp.s ) free(tmp.s);
    return 0;
}

static void filter_debug_print(token_t *toks, token_t **tok_ptrs, int ntoks)
{
    int i;
    for (i=0; i<ntoks; i++)
    {
        token_t *tok = toks ? &toks[i] : tok_ptrs[i];
        if ( tok->tok_type==TOK_VAL )
        {
            if ( tok->key )
                fprintf(stderr,"%s", tok->key);
            else if ( tok->tag )
                fprintf(stderr,"%s", tok->tag);
            else
                fprintf(stderr,"%e", tok->threshold);
        }
        else
            fprintf(stderr,"%c", TOKEN_STRING[tok->tok_type]);
        if ( tok->setter ) fprintf(stderr,"\t[setter %p]", tok->setter);
        fprintf(stderr,"\n");
    }
}

static void str_to_lower(char *str)
{
    while ( *str ) { *str = tolower(*str); str++; }
}
static int perl_exec(filter_t *flt, bcf1_t *line, token_t *rtok, token_t **stack, int nstack)
{
#if ENABLE_PERL_FILTERS

    PerlInterpreter *perl = flt->perl;
    if ( !perl ) error("Error: perl expression without a perl script name\n");

    dSP;
    ENTER;
    SAVETMPS;

    PUSHMARK(SP);
    int i,j, istack = nstack - rtok->nargs;
    for (i=istack+1; i<nstack; i++)
    {
        token_t *tok = stack[i];
        if ( tok->is_str )
            XPUSHs(sv_2mortal(newSVpvn(tok->str_value.s,tok->str_value.l)));
        else if ( tok->nvalues==1 )
            XPUSHs(sv_2mortal(newSVnv(tok->values[0])));
        else if ( tok->nvalues>1 )
        {
            AV *av = newAV();
            for (j=0; j<tok->nvalues; j++) av_push(av, newSVnv(tok->values[j]));
            SV *rv = newRV_inc((SV*)av);
            XPUSHs(rv);
        }
        else
        {
            bcf_double_set_missing(tok->values[0]);
            XPUSHs(sv_2mortal(newSVnv(tok->values[0])));
        }
    }
    PUTBACK;

    // A possible future todo: provide a means to select samples and indexes,
    // expressions like this don't work yet
    //          perl.filter(FMT/AD)[1:0]

    int nret = call_pv(stack[istack]->str_value.s, G_ARRAY);

    SPAGAIN;

    rtok->nvalues = nret;
    hts_expand(double, rtok->nvalues, rtok->mvalues, rtok->values);
    for (i=nret; i>0; i--)
    {
        rtok->values[i-1] = (double) POPn;
        if ( isnan(rtok->values[i-1]) ) bcf_double_set_missing(rtok->values[i-1]);
    }

    PUTBACK;
    FREETMPS;
    LEAVE;

#else
    error("\nPerl filtering requires running `configure --enable-perl-filters` at compile time.\n\n");
#endif
    return rtok->nargs;
}
static void perl_init(filter_t *filter, char **str)
{
    char *beg = *str;
    while ( *beg && isspace(*beg) ) beg++;
    if ( !*beg ) return;
    if ( strncasecmp("perl:", beg, 5) ) return;
#if ENABLE_PERL_FILTERS
    beg += 5;
    char *end = beg;
    while ( *end && *end!=';' ) end++;  // for now not escaping semicolons
    *str = end+1;

    if ( ++filter_ninit == 1 )
    {
        // must be executed only once, even for multiple filters; first time here
        int argc = 0;
        char **argv = NULL;
        char **env  = NULL;
        PERL_SYS_INIT3(&argc, &argv, &env);
    }

    filter->perl = perl_alloc();
    PerlInterpreter *perl = filter->perl;

    if ( !perl ) error("perl_alloc failed\n");
    perl_construct(perl);

    // name of the perl script to run
    char *rmme = (char*) calloc(end - beg + 1,1);
    memcpy(rmme, beg, end - beg);
    char *argv[] = { "", "" };
    argv[1] = expand_path(rmme);
    free(rmme);

    PL_origalen = 1;    // don't allow $0 change
    int ret = perl_parse(filter->perl, NULL, 2, argv, NULL);
    PL_exit_flags |= PERL_EXIT_DESTRUCT_END;
    if ( ret ) error("Failed to parse: %s\n", argv[1]);
    free(argv[1]);

    perl_run(perl);
#else
    error("\nPerl filtering requires running `configure --enable-perl-filters` at compile time.\n\n");
#endif
}
static void perl_destroy(filter_t *filter)
{
#if ENABLE_PERL_FILTERS
    if ( !filter->perl ) return;

    PerlInterpreter *perl = filter->perl;
    perl_destruct(perl);
    perl_free(perl);
    if ( --filter_ninit <= 0  )
    {
        // similarly to PERL_SYS_INIT3, can must be executed only once? todo: test
        PERL_SYS_TERM();
    }
#endif
}

// A very rudimentary heuristics to determine type, e.g. STR_TAG={} implies {str}.
// Throws an error on anything more complex and asks for an explicit type.
static void determine_ext_types(filter_t *filter, int ntok, token_t *tok)
{
    int i;
    for (i=0; i<ntok; i++)
    {
        if ( !tok[i].iext || tok[i].ht_type!=-1 ) continue;
        if ( !i || i+1==ntok ) break;       // first or last in the RPN
        if ( tok[i-1].ht_type==-1 ) break;  // previous type not set
        // todo: check if the next is an operator
        // todo: case when the order is reversed, {}=TAG
        tok[i].ht_type = tok[i-1].ht_type;
    }
    if ( i!=ntok )
        error("[%s:%d %s] Error: unable to determine the type, use explicit notation: %s\n",__FILE__,__LINE__,__FUNCTION__,filter->str);
    for (i=0; i<ntok; i++)
    {
        int j = tok[i].iext - 1;
        if ( j<0 ) continue;
        if ( filter->ext[j]!=-1 && filter->ext[j]!=tok[i].ht_type  )
            error("[%s:%d %s] FIXME: this should not happen %d vs %d, iext=%d\n",__FILE__,__LINE__,__FUNCTION__,filter->ext[j],tok[i].ht_type,j);
        filter->ext[j] = tok[i].ht_type;
        if ( tok[i].ht_type==BCF_HT_STR ) tok[i].is_str = 1;
    }
}


// Parse filter expression and convert to reverse polish notation. Dijkstra's shunting-yard algorithm
static filter_t *filter_init_(bcf_hdr_t *hdr, const char *str, int exit_on_error)
{
    filter_t *filter = (filter_t *) calloc(1,sizeof(filter_t));
    filter->str = strdup(str);
    filter->hdr = hdr;
    filter->max_unpack |= BCF_UN_STR;
    filter->exit_on_error = exit_on_error;

    int nops = 0, mops = 0;    // operators stack
    int nout = 0, mout = 0;    // filter tokens, RPN
    token_t *out = NULL;
    token_t *ops = NULL;
    char *tmp = filter->str;
    perl_init(filter, &tmp);

    int last_op = -1;
    while ( *tmp )
    {
        int len, ret;
        ret = filters_next_token(&tmp, &len);
        if ( ret==-1 ) error("Missing quotes in: %s\n", str);

        // fprintf(stderr,"token=[%c] .. [%s] %d\n", TOKEN_STRING[ret], tmp, len);
        // int i; for (i=0; i<nops; i++) fprintf(stderr," .%c", TOKEN_STRING[ops[i].tok_type]); fprintf(stderr,"\n");

        if ( ret==TOK_LFT )         // left bracket
        {
            nops++;
            hts_expand0(token_t, nops, mops, ops);
            ops[nops-1].tok_type = ret;
        }
        else if ( ret==TOK_RGT )    // right bracket
        {
            while ( nops>0 && ops[nops-1].tok_type!=TOK_LFT )
            {
                nout++;
                hts_expand0(token_t, nout, mout, out);
                out[nout-1] = ops[nops-1];
                memset(&ops[nops-1],0,sizeof(token_t));
                nops--;
            }
            if ( nops<=0 ) error("Could not parse: %s\n", str);
            memset(&ops[nops-1],0,sizeof(token_t));
            nops--;
        }
        else if ( ret==TOK_EXT )    // external value
        {
            nout++;
            hts_expand0(token_t, nout, mout, out);
            filters_init1_ext(filter, tmp, len, &out[nout-1]);
            tmp += len;
        }
        else if ( ret!=TOK_VAL )    // one of the operators
        {
            // detect unary minus: replace -value with -1*(value)
            if ( ret==TOK_SUB && last_op!=TOK_VAL && last_op!=TOK_RGT )
            {
                nout++;
                hts_expand0(token_t, nout, mout, out);
                token_t *tok = &out[nout-1];
                tok->tok_type  = TOK_VAL;
                tok->hdr_id    = -1;
                tok->pass_site = -1;
                tok->threshold = -1.0;
                tok->is_constant = 1;
                ret = TOK_MULT;
            }
            else if ( ret == -TOK_FUNC )
            {
                // this is different from TOK_PERLSUB,TOK_BINOM in that the expression inside the
                // brackets gets evaluated as normal expression
                nops++;
                hts_expand0(token_t, nops, mops, ops);
                token_t *tok = &ops[nops-1];
                tok->tok_type  = -ret;
                tok->hdr_id    = -1;
                tok->pass_site = -1;
                tok->threshold = -1.0;
                if ( !strncasecmp(tmp-len,"N_PASS",6) )
                {
                    filter->max_unpack |= BCF_UN_FMT;
                    tok->func = func_npass;
                    tok->tag = strdup("N_PASS");
                }
                else if ( !strncasecmp(tmp-len,"F_PASS",6) )
                {
                    filter->max_unpack |= BCF_UN_FMT;
                    tok->func = func_npass;
                    tok->tag = strdup("F_PASS");
                }
                else error("The function \"%s\" is not supported\n", tmp-len);
                continue;
            }
            else if ( ret < 0 )     // variable number of arguments: TOK_PERLSUB,TOK_BINOM
            {
                ret = -ret;

                tmp += len;
                char *beg = tmp;
                kstring_t rmme = {0,0,0};
                int i, margs, nargs = 0;

                if ( ret == TOK_PERLSUB )
                {
                    while ( *beg && ((isalnum(*beg) && !ispunct(*beg)) || *beg=='_') ) beg++;
                    if ( *beg!='(' ) error("Could not parse the expression: %s\n", str);

                    // the subroutine name
                    kputc('"', &rmme);
                    kputsn(tmp, beg-tmp, &rmme);
                    kputc('"', &rmme);
                    nout++;
                    hts_expand0(token_t, nout, mout, out);
                    filters_init1(filter, rmme.s, rmme.l, &out[nout-1]);
                    nargs++;
                }
                char *end = beg;
                while ( *end && *end!=')' ) end++;
                if ( !*end ) error("Could not parse the expression: %s\n", str);

                // subroutine arguments
                rmme.l = 0;
                kputsn(beg+1, end-beg-1, &rmme);
                char **rmme_list = hts_readlist(rmme.s, 0, &margs);
                for (i=0; i<margs; i++)
                {
                    nargs++;
                    nout++;
                    hts_expand0(token_t, nout, mout, out);
                    filters_init1(filter, rmme_list[i], strlen(rmme_list[i]), &out[nout-1]);
                    free(rmme_list[i]);
                }
                free(rmme_list);
                free(rmme.s);

                nout++;
                hts_expand0(token_t, nout, mout, out);
                token_t *tok = &out[nout-1];
                tok->tok_type  = ret;
                tok->nargs     = nargs;
                tok->hdr_id    = -1;
                tok->pass_site = -1;
                tok->threshold = -1.0;

                tmp = end + 1;
                continue;
            }
            else
            {
                while ( nops>0 && op_prec[ret] < op_prec[ops[nops-1].tok_type] )
                {
                    nout++;
                    hts_expand0(token_t, nout, mout, out);
                    out[nout-1] = ops[nops-1];
                    memset(&ops[nops-1],0,sizeof(token_t));
                    nops--;
                }
            }
            nops++;
            hts_expand0(token_t, nops, mops, ops);
            ops[nops-1].tok_type = ret;
        }
        else if ( !len )    // all tokes read or an error
        {
            if ( *tmp && !isspace(*tmp) ) error("Could not parse the expression: [%s]\n", str);
            break;     // all tokens read
        }
        else           // TOK_VAL: annotation name or value
        {
            nout++;
            hts_expand0(token_t, nout, mout, out);
            if ( tmp[len-1]==',' )
                filters_init1(filter, tmp, len-1, &out[nout-1]);
            else
                filters_init1(filter, tmp, len, &out[nout-1]);
            tmp += len;
        }
        last_op = ret;
    }
    while ( nops>0 )
    {
        if ( ops[nops-1].tok_type==TOK_LFT || ops[nops-1].tok_type==TOK_RGT ) error("Could not parse the expression: [%s]\n", filter->str);
        nout++;
        hts_expand0(token_t, nout, mout, out);
        out[nout-1] = ops[nops-1];
        memset(&ops[nops-1],0,sizeof(token_t));
        nops--;
    }

    if ( filter->status != FILTER_OK )
    {
        if ( mops ) free(ops);
        filter->filters   = out;
        filter->nfilters  = nout;
        return filter;
    }

    // Determine types of external variables from the context
    determine_ext_types(filter,nout,out);

    // In the special cases of TYPE and FILTER the BCF header IDs are yet unknown. Walk through the
    // list of operators and convert the strings (e.g. "PASS") to BCF ids. The string value token must be
    // just before or after the FILTER token and they must be followed with a comparison operator.
    // At this point we also initialize regex expressions which, in RPN, must precede the LIKE/NLIKE operator.
    // Additionally, treat "." as missing value rather than a string in numeric equalities; that
    // @file is only used with ID; etc.
    // This code is fragile: improve me.
    static int comma_separator_warned = 0;
    int i;
    for (i=0; i<nout; i++)
    {
        if ( i+1<nout && (out[i].tok_type==TOK_LT || out[i].tok_type==TOK_BT) && out[i+1].tok_type==TOK_EQ )
            error("Error parsing the expression: \"%s\"\n", filter->str);

        if ( out[i].hash )
        {
            int j = out[i+1].tok_type==TOK_VAL ? i+1 : i-1;
            if ( out[j].comparator!=filters_cmp_id )
            {
                if ( out[j].comparator ) error("Error: could not parse the expression with \"@file_name\" syntax (possible todo)\n");
                out[j].comparator = filters_cmp_string_hash;
            }
        }
        if ( out[i].tok_type==TOK_OR || out[i].tok_type==TOK_OR_VEC )
            out[i].func = vector_logic_or;
        if ( out[i].tok_type==TOK_AND || out[i].tok_type==TOK_AND_VEC )
            out[i].func = vector_logic_and;
        if ( out[i].tok_type==TOK_EQ || out[i].tok_type==TOK_NE )
        {
            // Look for j="." and k numeric type
            int j = i-1, k = i-2;
            if ( !out[j].is_str ) { k = i-1, j = i-2; }
            if ( out[j].is_str && out[j].key && !strcmp(".",out[j].key) )
            {
                int set_missing = 0;
                if ( out[k].hdr_id>0 )
                {
                    int type = bcf_hdr_id2type(filter->hdr,out[k].hl_type,out[k].hdr_id);
                    if ( type==BCF_HT_INT ) set_missing = 1;
                    else if ( type==BCF_HT_REAL ) set_missing = 1;
                }
                else if ( !strcmp("QUAL",out[k].tag) ) set_missing = 1;
                if ( set_missing ) { out[j].is_str = 0; out[j].is_missing = 1; bcf_double_set_missing(out[j].values[0]); }
            }
        }
        if ( out[i].tok_type==TOK_LIKE || out[i].tok_type==TOK_NLIKE )
        {
            int j = i-1;
            if ( !out[j].key )
                error("Could not parse the expression, wrong value for regex operator: %s\n", filter->str);
            out[j].regex = (regex_t *) malloc(sizeof(regex_t));
            int cflags = REG_NOSUB;
            int len = strlen(out[j].key);
            if ( len>2 && out[j].key[len-1]=='i' && out[j].key[len-2]=='/' && out[j].key[len-3]!='\\'  )
            {
                out[j].key[len-2] = 0;
                cflags |= REG_ICASE;
            }
            if ( regcomp(out[j].regex, out[j].key, cflags) )
                error("Could not compile the regex expression \"%s\": %s\n", out[j].key,filter->str);
        }
        if ( out[i].is_str && out[i].tok_type==TOK_VAL && out[i].key && strchr(out[i].key,',') )
        {
            int print_note = 0;
            if ( out[i+1].tok_type==TOK_EQ || (out[i+1].is_str && out[i+2].tok_type==TOK_EQ) ) print_note = 1;
            else if ( out[i+1].tok_type==TOK_NE || (out[i+1].is_str && out[i+2].tok_type==TOK_NE) ) print_note = 1;
            if ( print_note && !comma_separator_warned )
            {
                comma_separator_warned = 1;
                fprintf(stderr,
                    "Warning: comma is interpreted as a separator and OR logic is used in string comparisons.\n"
                    "         (Search the manual for \"Comma in strings\" to learn more.)\n");
            }
        }
        if ( out[i].tok_type!=TOK_VAL ) continue;
        if ( !out[i].tag ) continue;
        if ( out[i].setter==filters_set_type )
        {
            if ( i+1==nout || !out[i+1].key )
                error("Could not parse the expression: %s\n", filter->str);
            int itok, ival;
            if ( out[i+1].tok_type==TOK_EQ || out[i+1].tok_type==TOK_NE ) ival = i - 1, itok = i + 1;
            else if ( out[i+1].tok_type==TOK_LIKE || out[i+1].tok_type==TOK_NLIKE ) ival = i - 1, itok = i + 1;
            else if ( out[i+2].tok_type==TOK_EQ || out[i+2].tok_type==TOK_NE ) itok = i + 2, ival = i + 1;
            else if ( out[i+2].tok_type==TOK_LIKE || out[i+2].tok_type==TOK_NLIKE ) itok = i + 2, ival = i + 1;
            else error("[%s:%d %s] Could not parse the expression: %s\n",  __FILE__,__LINE__,__FUNCTION__, filter->str);
            if ( !strcasecmp(out[ival].key,"snp") || !strcasecmp(out[ival].key,"snps") ) { out[ival].threshold = VCF_SNP<<1; out[ival].is_str = 0; }
            else if ( !strcasecmp(out[ival].key,"indel") || !strcasecmp(out[ival].key,"indels") ) { out[ival].threshold = VCF_INDEL<<1; out[ival].is_str = 0; }
            else if ( !strcasecmp(out[ival].key,"mnp") || !strcasecmp(out[ival].key,"mnps") ) { out[ival].threshold = VCF_MNP<<1; out[ival].is_str = 0; }
            else if ( !strcasecmp(out[ival].key,"other") ) { out[ival].threshold = VCF_OTHER<<1; out[ival].is_str = 0; }
            else if ( !strcasecmp(out[ival].key,"bnd") ) { out[ival].threshold = VCF_BND<<1; out[ival].is_str = 0; }
            else if ( !strcasecmp(out[ival].key,"overlap") ) { out[ival].threshold = VCF_OVERLAP<<1; out[ival].is_str = 0; }
            else if ( !strcasecmp(out[ival].key,"ref") ) { out[ival].threshold = 1; out[ival].is_str = 0; }
            else error("The type \"%s\" not recognised: %s\n", out[ival].key, filter->str);
            out[ival].is_constant = 1;
            if ( out[itok].tok_type==TOK_LIKE || out[itok].tok_type==TOK_NLIKE ) out[itok].comparator = filters_cmp_bit_and;
            out[ival].tag = out[ival].key; out[ival].key = NULL;
            i = itok;
            continue;
        }
        if ( !strcmp(out[i].tag,"GT") )
        {
            if ( i+1==nout ) error("Could not parse the expression: %s\n", filter->str);
            int ival;
            if ( out[i+1].tok_type==TOK_EQ || out[i+1].tok_type==TOK_NE ) ival = i - 1;
            else if ( out[i+1].tok_type==TOK_LIKE || out[i+1].tok_type==TOK_NLIKE ) ival = i - 1;
            else if ( out[i+2].tok_type==TOK_EQ || out[i+2].tok_type==TOK_NE ) ival = i + 1;
            else if ( out[i+2].tok_type==TOK_LIKE || out[i+2].tok_type==TOK_NLIKE ) ival = i + 1;
            else error("[%s:%d %s] Could not parse the expression: %s\n",  __FILE__,__LINE__,__FUNCTION__, filter->str);

            if ( !out[ival].key ) error("Comparison between samples is not supported, sorry!\n");

            // assign correct setters and unify expressions, eg ar->ra, HOM->hom, etc
            if ( !strcasecmp(out[ival].key,"hom") ) { out[i].setter = filters_set_genotype3;  str_to_lower(out[ival].key); }
            else if ( !strcasecmp(out[ival].key,"het") ) { out[i].setter = filters_set_genotype3;  str_to_lower(out[ival].key); }
            else if ( !strcasecmp(out[ival].key,"hap") ) { out[i].setter = filters_set_genotype3;  str_to_lower(out[ival].key); }
            else if ( !strcasecmp(out[ival].key,"mis") ) { out[i].setter = filters_set_genotype4;  str_to_lower(out[ival].key); }
            else if ( !strcasecmp(out[ival].key,"ref") ) { out[i].setter = filters_set_genotype4;  str_to_lower(out[ival].key); }
            else if ( !strcasecmp(out[ival].key,"alt") ) { out[i].setter = filters_set_genotype4;  str_to_lower(out[ival].key); }
            else if ( !strcasecmp(out[ival].key,"rr") ) { out[i].setter = filters_set_genotype2;  str_to_lower(out[ival].key); }
            else if ( !strcasecmp(out[ival].key,"ra") || !strcasecmp(out[ival].key,"ar") ) { out[i].setter = filters_set_genotype2; out[ival].key[0]='r'; out[ival].key[1]='a'; }   // ra
            else if ( !strcmp(out[ival].key,"aA") || !strcmp(out[ival].key,"Aa") ) { out[i].setter = filters_set_genotype2; out[ival].key[0]='a'; out[ival].key[1]='A'; }   // aA
            else if ( !strcasecmp(out[ival].key,"aa") ) { out[i].setter = filters_set_genotype2; out[ival].key[0]='a'; out[ival].key[1]='a'; }  // aa
            else if ( !strcasecmp(out[ival].key,"a") ) { out[i].setter = filters_set_genotype2; out[ival].key[0]='a'; out[ival].key[1]=0; }  // a
            else if ( !strcasecmp(out[ival].key,"r") ) { out[i].setter = filters_set_genotype2; out[ival].key[0]='r'; out[ival].key[1]=0; }  // r
            continue;
        }
        if ( out[i].hl_type==BCF_HL_FLT )
        {
            if ( i+1==nout ) error("Could not parse the expression: %s\n", filter->str);
            int itok = i, ival;
            if ( out[i+1].tok_type==TOK_EQ || out[i+1].tok_type==TOK_NE ) ival = i - 1;
            else if ( out[i+1].tok_type==TOK_LIKE ) out[i+1].tok_type = TOK_IN, ival = i - 1;
            else if ( out[i+1].tok_type==TOK_NLIKE ) out[i+1].tok_type = TOK_NOT_IN, ival = i - 1;
            else if ( out[i+2].tok_type==TOK_EQ || out[i+2].tok_type==TOK_NE ) ival = ++i;
            else if ( out[i+2].tok_type==TOK_LIKE ) out[i+2].tok_type = TOK_IN, ival = ++i;
            else if ( out[i+2].tok_type==TOK_NLIKE ) out[i+2].tok_type = TOK_NOT_IN, ival = ++i;
            else error("[%s:%d %s] Could not parse the expression: %s\n",  __FILE__,__LINE__,__FUNCTION__, filter->str);
            if ( out[ival].tok_type!=TOK_VAL || !out[ival].key )
                error("[%s:%d %s] Could not parse the expression, an unquoted string value perhaps? %s\n", __FILE__,__LINE__,__FUNCTION__, filter->str);
            token_t *tok = &out[ival];
            char *bp = tok->key;
            tok->nvalues = 0;
            int has_missing = 0;
            while ( *bp )
            {
                char tmp, *ep = bp;
                while ( *ep && *ep!=';' ) ep++;
                tmp = *ep;
                *ep = 0;
                if ( !strcmp(".",bp) ) has_missing = 1;
                else
                {
                    tok->nvalues++;
                    hts_expand(double,tok->nvalues,tok->mvalues,tok->values);
                    int id = bcf_hdr_id2int(filter->hdr, BCF_DT_ID, bp);
                    if ( !bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_FLT,id) )
                        error("The filter \"%s\" not present in the VCF header\n", bp);
                    tok->values[tok->nvalues-1] = id;
                }
                *ep = tmp;
                if ( !tmp ) break;
                bp = ep + 1;
            }
            if ( has_missing && tok->nvalues ) error("The FILTER expression cannot contain missing value AND filters: \"%s\" (%d)\n",tok->key,tok->nvalues);
            out[ival].tag = tok->key;
            tok->key = NULL;
            out[itok].hdr_id = tok->hdr_id;
            continue;
        }
    }
    filter->nsamples = filter->max_unpack&BCF_UN_FMT ? bcf_hdr_nsamples(filter->hdr) : 0;
    for (i=0; i<nout; i++)
    {
        if ( out[i].tok_type==TOK_MAX )      { out[i].func = func_max; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_MIN ) { out[i].func = func_min; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_AVG ) { out[i].func = func_avg; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_MEDIAN ) { out[i].func = func_median; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_AVG ) { out[i].func = func_avg; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_STDEV ) { out[i].func = func_stddev; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_SUM ) { out[i].func = func_sum; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_ABS ) { out[i].func = func_abs; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_CNT ) { out[i].func = func_count; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_LEN ) { out[i].func = func_strlen; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_PHRED ) { out[i].func = func_phred; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_BINOM ) { out[i].func = func_binom; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_PERLSUB ) { out[i].func = perl_exec; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_sMAX ) { out[i].func = func_smpl_max; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_sMIN ) { out[i].func = func_smpl_min; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_sAVG ) { out[i].func = func_smpl_avg; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_sMEDIAN ) { out[i].func = func_smpl_median; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_sSTDEV ) { out[i].func = func_smpl_stddev; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_sSUM ) { out[i].func = func_smpl_sum; out[i].tok_type = TOK_FUNC; }
        hts_expand0(double,1,out[i].mvalues,out[i].values);
        if ( filter->nsamples )
        {
            out[i].pass_samples = (uint8_t*)malloc(filter->nsamples);
            int j;
            for (j=0; j<filter->nsamples; j++) out[i].pass_samples[j] = 0;
        }
    }

    if (0) filter_debug_print(out, NULL, nout);

    if ( mops ) free(ops);
    filter->filters   = out;
    filter->nfilters  = nout;
    filter->flt_stack = (token_t **)malloc(sizeof(token_t*)*nout);
    return filter;
}
filter_t *filter_parse(bcf_hdr_t *hdr, const char *str)
{
    return filter_init_(hdr, str, 0);
}
filter_t *filter_init(bcf_hdr_t *hdr, const char *str)
{
    return filter_init_(hdr, str, 1);
}

void filter_destroy(filter_t *filter)
{
    perl_destroy(filter);
    int i;
    for (i=0; i<filter->nfilters; i++)
    {
        if ( filter->filters[i].key ) free(filter->filters[i].key);
        free(filter->filters[i].str_value.s);
        free(filter->filters[i].tag);
        free(filter->filters[i].idxs);
        free(filter->filters[i].usmpl);
        free(filter->filters[i].values);
        free(filter->filters[i].pass_samples);
        if (filter->filters[i].hash) khash_str2int_destroy_free(filter->filters[i].hash);
        if (filter->filters[i].regex)
        {
            regfree(filter->filters[i].regex);
            free(filter->filters[i].regex);
        }
    }
    for (i=0; i<filter->nundef_tag; i++) free(filter->undef_tag[i]);
    for (i=0; i<filter->nused_tag; i++) free(filter->used_tag[i]);
    free(filter->ext);
    free(filter->undef_tag);
    free(filter->used_tag);
    free(filter->cached_GT.buf);
    free(filter->cached_GT.mask);
    free(filter->filters);
    free(filter->flt_stack);
    free(filter->str);
    free(filter->tmpi);
    free(filter->tmpf);
    free(filter->tmps.s);
    free(filter);
}

int filter_test_ext(filter_t *filter, bcf1_t *rec, const uint8_t **samples, const void **ext)
{
    if ( !filter->n_ext )
        return filter_test(filter,rec,samples);

    int i;
    for (i=0; i<filter->nfilters; i++)
    {
        token_t *tok = &filter->filters[i];
        if ( !tok->iext ) continue;
        if ( !ext[tok->iext-1] )
        {
            tok->is_missing = 1;
            tok->nvalues = 0;
            if ( filter->ext[tok->iext-1]==BCF_HT_STR ) tok->str_value.l = 0;
            continue;
        }
        tok->is_missing = 0;
        tok->nvalues = 1;
        if ( filter->ext[tok->iext-1]==BCF_HT_STR )
        {
            tok->str_value.l = 0;
            kputs((const char*)ext[tok->iext-1],&tok->str_value);
            tok->nvalues = tok->str_value.l;
        }
        else if ( filter->ext[tok->iext-1]==BCF_HT_INT ) tok->values[0] = *((const int*)ext[tok->iext-1]);
        else if ( filter->ext[tok->iext-1]==BCF_HT_REAL ) tok->values[0] = *((const float*)ext[tok->iext-1]);
    }
    return filter_test(filter,rec,samples);
}

int filter_test(filter_t *filter, bcf1_t *line, const uint8_t **samples)
{
    if ( filter->status != FILTER_OK ) error("Error: the caller did not check the filter status\n");
    bcf_unpack(line, filter->max_unpack);

    int i, nstack = 0;
    for (i=0; i<filter->nfilters; i++)
    {
        filter->filters[i].pass_site = 0;
        if ( filter->filters[i].tok_type == TOK_VAL )
        {
            if ( filter->filters[i].setter )    // variable, query the VCF line
                filter->filters[i].setter(filter, line, &filter->filters[i]);
            else if ( filter->filters[i].key )  // string constant
            {
                filter->filters[i].str_value.l = 0;
                kputs(filter->filters[i].key, &filter->filters[i].str_value);
                filter->filters[i].nvalues = filter->filters[i].str_value.l;
            }
            else if ( filter->filters[i].is_constant )   // numeric constant
            {
                filter->filters[i].values[0] = filter->filters[i].threshold;
                filter->filters[i].nvalues   = 1;
            }
            filter->flt_stack[nstack++] = &filter->filters[i];
            continue;
        }
        else if ( filter->filters[i].func )
        {
            int nargs = filter->filters[i].func(filter, line, &filter->filters[i], filter->flt_stack, nstack);
            filter->flt_stack[nstack-nargs] = &filter->filters[i];
            if ( --nargs > 0 ) nstack -= nargs;
            continue;
        }
        if ( nstack<2 )
            error("Error occurred while processing the filter \"%s\" (1:%d)\n", filter->str,nstack);  // too few values left on the stack

        int is_str  = filter->flt_stack[nstack-1]->is_str + filter->flt_stack[nstack-2]->is_str;

        if ( filter->filters[i].tok_type == TOK_ADD )
        {
            VECTOR_ARITHMETICS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],&filter->filters[i],+,(double));
            filter->flt_stack[nstack-2] = &filter->filters[i];
            nstack--;
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_SUB )
        {
            VECTOR_ARITHMETICS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],&filter->filters[i],-,(double));
            filter->flt_stack[nstack-2] = &filter->filters[i];
            nstack--;
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_MULT )
        {
            VECTOR_ARITHMETICS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],&filter->filters[i],*,(double));
            filter->flt_stack[nstack-2] = &filter->filters[i];
            nstack--;
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_DIV )
        {
            VECTOR_ARITHMETICS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],&filter->filters[i],/,(double));
            filter->flt_stack[nstack-2] = &filter->filters[i];
            nstack--;
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_MODULO )
        {
            VECTOR_ARITHMETICS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],&filter->filters[i],%,(int));
            filter->flt_stack[nstack-2] = &filter->filters[i];
            nstack--;
            continue;
        }

        // ideally, these comparators would become func, but this would require more work in init1()
        if ( filter->filters[i].comparator )
        {
            filter->filters[i].comparator(filter->flt_stack[nstack-1],filter->flt_stack[nstack-2],&filter->filters[i],line);
        }
        else if ( filter->flt_stack[nstack-1]->comparator )
        {
            filter->flt_stack[nstack-1]->comparator(filter->flt_stack[nstack-1],filter->flt_stack[nstack-2],&filter->filters[i],line);
        }
        else if ( filter->flt_stack[nstack-2]->comparator )
        {
            filter->flt_stack[nstack-2]->comparator(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],&filter->filters[i],line);
        }
        else if ( is_str==2 )
        {
            cmp_vector_strings(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],&filter->filters[i]);
        }
        else
        {
            if ( is_str==1 ) error("Error: cannot use arithmetic operators to compare strings and numbers\n");

            // Determine what to do with one [1] or both [2] sides missing. The first field [0] gives [1]|[2]
            int missing_logic[] = {0,0,0};
            if ( filter->filters[i].tok_type == TOK_EQ ) { missing_logic[0] = missing_logic[2] = 1; }
            if ( filter->filters[i].tok_type == TOK_NE ) { missing_logic[0] = missing_logic[1] = 1; }

            if ( filter->filters[i].tok_type == TOK_EQ )
                CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],&filter->filters[i],==,missing_logic)
            else if ( filter->filters[i].tok_type == TOK_NE )
                CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],&filter->filters[i],!=,missing_logic)
            else if ( filter->filters[i].tok_type == TOK_LE )
                CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],&filter->filters[i],<=,missing_logic)
            else if ( filter->filters[i].tok_type == TOK_LT )
                CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],&filter->filters[i],<,missing_logic)
            else if ( filter->filters[i].tok_type == TOK_BT )
                CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],&filter->filters[i],>,missing_logic)
            else if ( filter->filters[i].tok_type == TOK_BE )
                CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],&filter->filters[i],>=,missing_logic)
            else
                error("todo: %s:%d .. type=%d\n", __FILE__,__LINE__,filter->filters[i].tok_type);
        }
        filter->flt_stack[nstack-2] = &filter->filters[i];
        nstack--;
    }
    if ( nstack>1 ) error("Error occurred while processing the filter \"%s\" (2:%d)\n", filter->str,nstack);    // too few values left on the stack
    if ( samples )
    {
        *samples = filter->max_unpack&BCF_UN_FMT ? filter->flt_stack[0]->pass_samples : NULL;
        if ( *samples && !filter->flt_stack[0]->nsamples )
        {
            for (i=0; i<filter->nsamples; i++)
                filter->flt_stack[0]->pass_samples[i] = filter->flt_stack[0]->pass_site;
        }
    }
    return filter->flt_stack[0]->pass_site;
}

int filter_max_unpack(filter_t *flt)
{
    return flt->max_unpack;
}
const int *filter_ext_types(filter_t *filter, int *n_ext)
{
    *n_ext = filter->n_ext;
    return filter->ext;
}
const double *filter_get_doubles(filter_t *filter, int *nval, int *nval1)
{
    token_t *tok = filter->flt_stack[0];
    if ( tok->nvalues )
    {
        *nval  = tok->nvalues;
        *nval1 = tok->nval1;
    }
    else
    {
        if ( !tok->values ) error("fixme in filter_get_doubles(): %s\n", filter->str);
        *nval  = 1;
        *nval1 = 1;
        tok->values[0] = filter->flt_stack[0]->pass_site;
    }
    return tok->values;
}

void filter_set_samples(filter_t *filter, const uint8_t *samples)
{
    int i,j;
    for (i=0; i<filter->nfilters; i++)
    {
        if ( !filter->filters[i].nsamples ) continue;
        for (j=0; j<filter->filters[i].nsamples; j++) filter->filters[i].usmpl[j] = samples[j];
    }
}

int filter_status(filter_t *filter)
{
    return filter->status;
}

