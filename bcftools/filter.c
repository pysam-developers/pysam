/*  filter.c -- filter expressions.

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

#include <ctype.h>
#include <stdlib.h>
#include <strings.h>
#include <errno.h>
#include <math.h>
#include <wordexp.h>
#include <regex.h>
#include <htslib/khash_str2int.h>
#include "filter.h"
#include "bcftools.h"
#include <htslib/hts_defs.h>
#include <htslib/vcfutils.h>

#ifndef __FUNCTION__
#  define __FUNCTION__ __func__
#endif

uint64_t bcf_double_missing    = 0x7ff0000000000001;
uint64_t bcf_double_vector_end = 0x7ff0000000000002;
static inline void bcf_double_set(double *ptr, uint64_t value)
{
    union { uint64_t i; double d; } u;
    u.i = value;
    *ptr = u.d;
}
static inline int bcf_double_test(double d, uint64_t value)
{
    union { uint64_t i; double d; } u;
    u.d = d;
    return u.i==value ? 1 : 0;
}
#define bcf_double_set_vector_end(x) bcf_double_set(&(x),bcf_double_vector_end)
#define bcf_double_set_missing(x)    bcf_double_set(&(x),bcf_double_missing)
#define bcf_double_is_vector_end(x)  bcf_double_test((x),bcf_double_vector_end)
#define bcf_double_is_missing(x)     bcf_double_test((x),bcf_double_missing)


typedef struct _token_t
{
    // read-only values, same for all VCF lines
    int tok_type;       // one of the TOK_* keys below
    char *key;          // set only for string constants, otherwise NULL
    char *tag;          // for debugging and printout only, VCF tag name
    double threshold;   // filtering threshold
    int hdr_id, type;   // BCF header lookup ID and one of BCF_HT_* types
    int idx;            // 0-based index to VCF vectors, -1: not a vector, 
                        //  -2: list (e.g. [0,1,2] or [1..3] or [1..] or any field[*], which is equivalent to [0..])
    int *idxs, nidxs;   // set indexes to 0 to exclude, to 1 to include, and last element negative if unlimited
    void (*setter)(filter_t *, bcf1_t *, struct _token_t *);
    int (*comparator)(struct _token_t *, struct _token_t *, int op_type, bcf1_t *);
    void *hash;         // test presence of str value in the hash via comparator
    regex_t *regex;     // precompiled regex for string comparison

    // modified on filter evaluation at each VCF line
    double *values;
    kstring_t str_value;
    int is_str, is_missing; // is_missing is set only for constants, variables are controled via nvalues
    int pass_site;          // -1 not applicable, 0 fails, >0 pass
    uint8_t *pass_samples;  // status of individual samples
    int nsamples;           // number of samples
    int nvalues, mvalues;   // number of used values: n=0 for missing values, n=1 for scalars, for strings n=str_value.l
    int nstr1;              // per-sample string length, set only with str_value.l>0 && nsamples>1
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

//                      0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
//                        ( ) [ < = > ] ! | &  +  -  *  /  M  m  a  A  O  ~  ^  S  .  l
static int op_prec[] = {0,1,1,5,5,5,5,5,5,2,3, 6, 6, 7, 7, 8, 8, 8, 3, 2, 5, 5, 8, 8, 8};
#define TOKEN_STRING "x()[<=>]!|&+-*/MmaAO~^f"

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

    if ( !strncasecmp(tmp,"MAX(",4) ) { (*str) += 3; return TOK_MAX; }
    if ( !strncasecmp(tmp,"MIN(",4) ) { (*str) += 3; return TOK_MIN; }
    if ( !strncasecmp(tmp,"AVG(",4) ) { (*str) += 3; return TOK_AVG; }
    if ( !strncasecmp(tmp,"SUM(",4) ) { (*str) += 3; return TOK_SUM; }
    if ( !strncasecmp(tmp,"ABS(",4) ) { (*str) += 3; return TOK_ABS; }
    if ( !strncasecmp(tmp,"STRLEN(",7) ) { (*str) += 6; return TOK_LEN; }
    if ( !strncasecmp(tmp,"%MAX(",5) ) { (*str) += 4; return TOK_MAX; } // for backward compatibility
    if ( !strncasecmp(tmp,"%MIN(",5) ) { (*str) += 4; return TOK_MIN; } // for backward compatibility
    if ( !strncasecmp(tmp,"%AVG(",5) ) { (*str) += 4; return TOK_AVG; } // for backward compatibility
    if ( !strncasecmp(tmp,"%SUM(",5) ) { (*str) += 4; return TOK_SUM; } // for backward compatibility
    if ( !strncasecmp(tmp,"INFO/",5) ) tmp += 5;
    if ( !strncasecmp(tmp,"FORMAT/",7) ) tmp += 7;
    if ( !strncasecmp(tmp,"FMT/",4) ) tmp += 4;

    if ( tmp[0]=='@' )  // file name
    {
        while ( *tmp && !isspace(*tmp) && *tmp!='=' && *tmp!='!' ) tmp++;
        *len = tmp - (*str);
        return TOK_VAL;
    }

    int square_brackets = 0;
    while ( tmp[0] )
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
        if ( tmp[0]=='*' && !square_brackets ) break;
        if ( tmp[0]=='-' && !square_brackets ) break;
        if ( tmp[0]=='/' ) break;
        if ( tmp[0]=='~' ) break;
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

    *len = tmp - (*str);
    return TOK_VAL;
}

static void filters_set_qual(filter_t *flt, bcf1_t *line, token_t *tok)
{
    float *ptr = &line->qual;
    if ( bcf_float_is_missing(*ptr) )
        tok->nvalues = 0;
    else
    {
        tok->values[0] = (double)line->qual;
        tok->nvalues = 1;
    }
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
        tok->nvalues = 0;
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
static int filters_cmp_bit_and(token_t *atok, token_t *btok, int op_type, bcf1_t *line)
{
    int a = (int)(atok->nvalues?atok->values[0]:atok->threshold);
    int b = (int)(btok->nvalues?btok->values[0]:btok->threshold);
    if ( op_type==TOK_LIKE ) return a&b ? 1 : 0;
    return a&b ? 0 : 1;
}
static int filters_cmp_filter(token_t *atok, token_t *btok, int op_type, bcf1_t *line)
{
    int i;
    if ( op_type==TOK_NE )  // AND logic: none of the filters can match
    {
        if ( !line->d.n_flt )
        {
            if ( atok->hdr_id==-1 ) return 0;   // missing value
            return 1; // no filter present, eval to true
        }
        for (i=0; i<line->d.n_flt; i++)
            if ( atok->hdr_id==line->d.flt[i] ) return 0;
        return 1;
    }
    // TOK_EQ with OR logic: at least one of the filters must match
    if ( !line->d.n_flt )
    {
        if ( atok->hdr_id==-1 ) return 1;
        return 0; // no filter present, eval to false
    }
    for (i=0; i<line->d.n_flt; i++)
        if ( atok->hdr_id==line->d.flt[i] ) return 1;
    return 0;
}
static int filters_cmp_id(token_t *atok, token_t *btok, int op_type, bcf1_t *line)
{
    // multiple IDs not supported yet (easy to add though)

    if ( btok->hash )
    {
        token_t *tmp = atok; atok = btok; btok = tmp;
    }
    if ( atok->hash )
    {
        int ret = khash_str2int_has_key(atok->hash, line->d.id);
        if ( op_type==TOK_EQ ) return ret;
        return ret ? 0 : 1;
    }

    if ( op_type==TOK_EQ ) return strcmp(btok->str_value.s,line->d.id) ? 0 : 1;
    return strcmp(btok->str_value.s,line->d.id) ? 1 : 0;
}

/**
 *  bcf_get_info_value() - get single INFO value, int64_t or double
 *  @line:      BCF line
 *  @info_id:   tag ID, as returned by bcf_hdr_id2int
 *  @ivec:      0-based index to retrieve, -1 when single value is expected
 *  @vptr:      pointer to memory location of sufficient size to accomodate
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

    #define BRANCH(type_t, is_missing, is_vector_end, out_type_t) { \
        type_t *p = (type_t *) info->vptr; \
        for (j=0; j<ivec && j<info->len; j++) \
        { \
            if ( is_vector_end ) return 0; \
        } \
        if ( is_missing ) return 0; \
        *((out_type_t*)value) = p[j]; \
        return 1; \
    }
    switch (info->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  p[j]==bcf_int8_missing,  p[j]==bcf_int8_vector_end,  int64_t); break;
        case BCF_BT_INT16: BRANCH(int16_t, p[j]==bcf_int16_missing, p[j]==bcf_int16_vector_end, int64_t); break;
        case BCF_BT_INT32: BRANCH(int32_t, p[j]==bcf_int32_missing, p[j]==bcf_int32_vector_end, int64_t); break;
        case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(p[j]), bcf_float_is_vector_end(p[j]), double); break;
        default: fprintf(stderr,"todo: type %d\n", info->type); exit(1); break;
    }
    #undef BRANCH
    return -1;  // this shouldn't happen
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
        int64_t value;
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
        int i, end = tok->idxs[tok->nidxs-1] < 0 ? n - 1 : tok->nidxs - 1;
        if ( end >= n ) end = n - 1;
        char *beg = tok->str_value.s, *dst = flt->tmps.s;
        for (i=0; i<=end; i++)
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
    int i;
    if ( (tok->nvalues=bcf_get_format_int32(flt->hdr,line,tok->tag,&flt->tmpi,&flt->mtmpi))<0 )
    {
        tok->nvalues = tok->nsamples = 0;
        return;
    }
    if ( tok->idx >= -1 )    // scalar or vector index
    {
        hts_expand(double,flt->nsamples,tok->mvalues,tok->values);
        int nvals = tok->nvalues / line->n_sample;
        int idx = tok->idx >= 0 ? tok->idx : 0;
        int is_missing = 1;
        int32_t *ptr = flt->tmpi;
        for (i=0; i<line->n_sample; i++)
        {
            if ( ptr[idx]==bcf_int32_missing || ptr[idx]==bcf_int32_vector_end ) 
                bcf_double_set_missing(tok->values[i]);
            else
            {
                tok->values[i] = ptr[idx];
                is_missing = 0;
            }
            ptr += nvals;
        }
        if ( is_missing ) tok->nvalues = 0;
        else tok->nvalues = line->n_sample;
        tok->nsamples = tok->nvalues;
        return;
    }
    if ( tok->idx == -2 )
    {
        hts_expand(double,tok->nvalues,tok->mvalues,tok->values);
        int nvals = tok->nvalues / line->n_sample;
        int idx = tok->idx >= 0 ? tok->idx : 0;
        int is_missing = 1;
        int k, j = 0, end = tok->idxs[tok->nidxs-1] < 0 ? nvals - 1 : tok->nidxs - 1;
        if ( end >= nvals ) end = nvals - 1;
        int32_t *ptr = flt->tmpi;
        for (i=0; i<line->n_sample; i++)
        {
            for (k=0; k<=end; k++)
                if ( k>=tok->nidxs || tok->idxs[k] )
                {
                    if ( ptr[k]==bcf_int32_missing || ptr[k]==bcf_int32_vector_end )
                        bcf_double_set_missing(tok->values[j]);
                    else
                    {
                        tok->values[j] = ptr[k];
                        is_missing = 0;
                    }
                    j++;
                }
            ptr += nvals;
        }
        if ( is_missing ) tok->nvalues = tok->nsamples = 0;
        else
        {
            tok->nsamples = line->n_sample;
            tok->nvalues = j;
        }
        return;
    }
}
static void filters_set_format_float(filter_t *flt, bcf1_t *line, token_t *tok)
{
    int i;
    if ( (tok->nvalues=bcf_get_format_float(flt->hdr,line,tok->tag,&flt->tmpf,&flt->mtmpf))<0 )
    {
        tok->nvalues = tok->nsamples = 0;
        return;
    }
    if ( tok->idx >= -1 )    // scalar or vector index
    {
        hts_expand(double,flt->nsamples,tok->mvalues,tok->values);
        int nvals = tok->nvalues / line->n_sample;
        int idx = tok->idx >= 0 ? tok->idx : 0;
        int is_missing = 1;
        float *ptr = flt->tmpf;
        for (i=0; i<line->n_sample; i++)
        {
            if ( bcf_float_is_missing(ptr[idx]) || bcf_float_is_vector_end(ptr[idx]) )
                bcf_double_set_missing(tok->values[i]);
            else
            {
                tok->values[i] = ptr[idx];
                is_missing = 0;
            }
            ptr += nvals;
        }
        if ( is_missing ) tok->nvalues = 0;
        else tok->nvalues = line->n_sample;
        tok->nsamples = tok->nvalues;
        return;
    }
    if ( tok->idx == -2 )
    {
        hts_expand(double,tok->nvalues,tok->mvalues,tok->values);
        int nvals = tok->nvalues / line->n_sample;
        int idx = tok->idx >= 0 ? tok->idx : 0;
        int is_missing = 1;
        int k, j = 0, end = tok->idxs[tok->nidxs-1] < 0 ? nvals - 1 : tok->nidxs - 1;
        if ( end >= nvals ) end = nvals - 1;
        float *ptr = flt->tmpf;
        for (i=0; i<line->n_sample; i++)
        {
            for (k=0; k<=end; k++)
                if ( k>=tok->nidxs || tok->idxs[k] )
                {
                    if ( bcf_float_is_missing(ptr[k]) || bcf_float_is_vector_end(ptr[k]) )
                        bcf_double_set_missing(tok->values[j]);
                    else
                    {
                        tok->values[j] = ptr[k];
                        is_missing = 0;
                    }
                    j++;
                }
            ptr += nvals;
        }
        if ( is_missing ) tok->nvalues = tok->nsamples = 0;
        else
        {
            tok->nsamples = line->n_sample;
            tok->nvalues = j;
        }
        return;
    }
}
static void filters_set_format_string(filter_t *flt, bcf1_t *line, token_t *tok)
{
    tok->str_value.l = tok->nvalues = 0;
    if ( !line->n_sample ) return;

    int ndim = tok->str_value.m;
    int nstr = bcf_get_format_char(flt->hdr, line, tok->tag, &tok->str_value.s, &ndim);
    tok->str_value.m = ndim;

    if ( nstr<=0 ) return;

    if ( tok->idx == -1 || (tok->idx==-2 && tok->idxs[0]==-1) ) // scalar or keep all values of a vector: TAG[*]
    {
        tok->nsamples = line->n_sample;
        tok->nstr1 = ndim / line->n_sample;
        tok->nvalues = tok->str_value.l = nstr;
        return;
    }

    int nstr1 = nstr / line->n_sample;

    // vector, one or multiple indices
    int i;
    for (i=0; i<line->n_sample; i++)
    {
        char *dst = tok->str_value.s + i*nstr1, *str = dst;
        int nval = 0, ibeg = 0;
        while ( ibeg < nstr1 )
        {
            int iend = ibeg + 1;
            while ( iend < nstr1 && str[iend] && str[iend]!=',' ) iend++;

            int keep = 0;
            if ( tok->idx >=0 )
                keep = tok->idx==nval ? 1 : 0;
            else if ( nval < tok->nidxs )
                keep = tok->idxs[nval] ? 1 : 0;
            else if ( tok->idxs[tok->nidxs-1] < 0 )
                keep = 1;

            if ( keep )
            {
                if ( ibeg>0 ) memmove(dst, str+ibeg, iend-ibeg+1);
                dst += iend - ibeg + 1;
                if ( tok->idx>=0 ) break;
            }
            if ( !str[iend] ) break;
            ibeg = iend + 1;
            nval++;
        }
        if ( dst==str ) { dst[0] = '.'; dst+=2; }
        if ( dst - str < nstr1 ) memset(dst-1, 0, nstr1 - (dst - str));
    }
    tok->nvalues = tok->str_value.l = nstr;
    tok->nstr1 = nstr1;
    tok->nsamples = line->n_sample;
}
static void _filters_set_genotype(filter_t *flt, bcf1_t *line, token_t *tok, int type)
{
    bcf_fmt_t *fmt = bcf_get_fmt(flt->hdr, line, "GT");
    if ( !fmt )
    {
        tok->nvalues = tok->str_value.l = 0;
        return;
    }
    
    int i,j, nsmpl = bcf_hdr_nsamples(flt->hdr), nvals = type==2 ? 3 : 4;
    if ( tok->str_value.m <= nvals*nsmpl )
    {
        tok->str_value.m = nvals*nsmpl + 1;
        tok->str_value.s = (char*)realloc(tok->str_value.s, tok->str_value.m);
    }

#define BRANCH_INT(type_t,vector_end) \
    { \
        for (i=0; i<line->n_sample; i++) \
        { \
            type_t *ptr = (type_t*) (fmt->p + i*fmt->size); \
            int is_het = 0, has_ref = 0, missing = 0; \
            for (j=0; j<fmt->n; j++) \
            { \
                if ( ptr[j]==vector_end ) break; /* smaller ploidy */ \
                if ( bcf_gt_is_missing(ptr[j]) ) { missing=1; break; } /* missing allele */ \
                int ial = ptr[j]; \
                if ( bcf_gt_allele(ial)==0 ) has_ref = 1; \
                if ( j>0 ) \
                { \
                    int jal = ptr[j-1]; \
                    if ( bcf_gt_allele(ial)!=bcf_gt_allele(jal) ) is_het = 1; \
                } \
            } \
            char *dst = &tok->str_value.s[nvals*i]; \
            if ( !j || missing ) dst[0]='.', dst[1]=0; /* ., missing genotype */ \
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
        case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_vector_end); break;
        default: error("The GT type is not lineognised: %d at %s:%d\n",fmt->type, bcf_seqname(flt->hdr,line),line->pos+1); break;
    }
#undef BRANCH_INT
    tok->nsamples = nsmpl;
    tok->nvalues = tok->str_value.l = nvals*nsmpl;
    tok->str_value.s[tok->str_value.l] = 0;
    tok->nstr1 = nvals;
}
static void filters_set_genotype2(filter_t *flt, bcf1_t *line, token_t *tok) { _filters_set_genotype(flt, line, tok, 2); }
static void filters_set_genotype3(filter_t *flt, bcf1_t *line, token_t *tok) { _filters_set_genotype(flt, line, tok, 3); }

static void filters_set_genotype_string(filter_t *flt, bcf1_t *line, token_t *tok)
{
    bcf_fmt_t *fmt = bcf_get_fmt(flt->hdr, line, "GT");
    if ( !fmt )
    {
        tok->nvalues = tok->nsamples = 0;
        return;
    }
    int i, blen = 4, nsmpl = bcf_hdr_nsamples(flt->hdr);
    kstring_t str;

gt_length_too_big:
    tok->str_value.l = 0;
    for (i=0; i<nsmpl; i++)
    {
        int plen = tok->str_value.l;

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
    tok->nsamples = nsmpl;
    tok->nvalues = tok->str_value.l;
    tok->nstr1 = blen;
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
        tok->idx = 0;
    }
    else if ( tok->idx==-2 )
    {
        int i, j = 0, end = tok->idxs[tok->nidxs-1] < 0 ? line->n_allele - 1 : tok->nidxs - 1;
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
    if ( fmt->type!=BCF_BT_INT8 ) error("TODO: the GT fmt_type is not int8\n");

    int j,nmissing = 0;
    for (i=0; i<line->n_sample; i++)
    {
        int8_t *ptr = (int8_t*) (fmt->p + i*fmt->size);
        for (j=0; j<fmt->n; j++)
        {
            if ( ptr[j]==bcf_int8_vector_end ) break;
            if ( ptr[j]==bcf_gt_missing ) { nmissing++; break; }
        }
    }
    tok->nvalues = 1;
    tok->values[0] = tok->tag[0]=='N' ? nmissing : (double)nmissing / line->n_sample;
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

static void set_max(filter_t *flt, bcf1_t *line, token_t *tok)
{
    double val = -HUGE_VAL;
    int i;
    for (i=0; i<tok->nvalues; i++)
    {
        if ( !bcf_double_is_missing(tok->values[i]) && val < tok->values[i] ) val = tok->values[i];
    }
    tok->values[0] = val;
    tok->nvalues   = 1;
    tok->nsamples  = 0;
}
static void set_min(filter_t *flt, bcf1_t *line, token_t *tok)
{
    double val = HUGE_VAL;
    int i;
    for (i=0; i<tok->nvalues; i++)
        if ( !bcf_double_is_missing(tok->values[i]) && val > tok->values[i] ) val = tok->values[i];
    tok->values[0] = val;
    tok->nvalues   = 1;
    tok->nsamples  = 0;
}
static void set_avg(filter_t *flt, bcf1_t *line, token_t *tok)
{
    double val = 0;
    int i, n = 0;
    for (i=0; i<tok->nvalues; i++)
        if ( !bcf_double_is_missing(tok->values[i]) ) { val += tok->values[i]; n++; }
    tok->values[0] = n ? val / n : 0;
    tok->nvalues   = 1;
    tok->nsamples  = 0;
}
static void set_sum(filter_t *flt, bcf1_t *line, token_t *tok)
{
    double val = 0;
    int i, n = 0;
    for (i=0; i<tok->nvalues; i++)
        if ( !bcf_double_is_missing(tok->values[i]) ) { val += tok->values[i]; n++; }
    tok->values[0] = val;
    tok->nvalues   = 1;
    tok->nsamples  = 0;
}
static void set_abs(filter_t *flt, bcf1_t *line, token_t *tok)
{
    if ( tok->is_str ) error("ABS() can be applied only on numeric values\n");
    int i;
    for (i=0; i<tok->nvalues; i++)
        tok->values[i] = fabs(tok->values[i]);
}
static void set_strlen(filter_t *flt, bcf1_t *line, token_t *tok)
{
    tok->is_str = 0;
    if ( !tok->str_value.l ) return;
    if ( tok->idx==-2 )
    {
        int i = 0;
        char *ss = tok->str_value.s;
        while ( *ss )
        {
            char *se = ss;
            while ( *se && *se!=',' ) se++;
            hts_expand(double, i+1, tok->mvalues, tok->values);
            if ( !*se ) tok->values[i] = strlen(ss);
            else
            {
                *se = 0;
                tok->values[i] = strlen(ss);
                *se = ',';
            }
            ss = *se ? se + 1 : se;
            i++;
        }
        tok->nvalues = i;
    }
    else
    {
        tok->values[0] = strlen(tok->str_value.s);
        tok->nvalues = 1;
    }
    tok->str_value.l = 0;
}
#define VECTOR_ARITHMETICS(atok,btok,AOP) \
{ \
    int i, has_values = 0; \
    if ( !(atok)->nvalues || !(btok)->nvalues ) /* missing values */ \
    { \
        (atok)->nvalues = 0; (atok)->nsamples = 0; \
    } \
    else \
    { \
        if ( ((atok)->nsamples && (btok)->nsamples) || (!(atok)->nsamples && !(btok)->nsamples)) \
        { \
            for (i=0; i<(atok)->nvalues; i++) \
            { \
                if ( bcf_double_is_missing((atok)->values[i]) ) continue; \
                if ( bcf_double_is_missing((btok)->values[i]) ) { bcf_double_set_missing((atok)->values[i]); continue; } \
                has_values = 1; \
                (atok)->values[i] = (atok)->values[i] AOP (btok)->values[i]; \
            } \
        } \
        else if ( (btok)->nsamples ) \
        { \
            hts_expand(double,(btok)->nvalues,(atok)->mvalues,(atok)->values); \
            for (i=0; i<(btok)->nvalues; i++) \
            { \
                if ( bcf_double_is_missing((atok)->values[0]) || bcf_double_is_missing((btok)->values[i]) ) \
                { \
                    bcf_double_set_missing((atok)->values[i]); \
                    continue; \
                } \
                has_values = 1; \
                (atok)->values[i] = (atok)->values[0] AOP (btok)->values[i]; \
            } \
            (atok)->nvalues  = (btok)->nvalues; \
            (atok)->nsamples = (btok)->nsamples; \
        } \
        else if ( (atok)->nsamples ) \
        { \
            for (i=0; i<(atok)->nvalues; i++) \
            { \
                if ( bcf_double_is_missing((atok)->values[i]) || bcf_double_is_missing((btok)->values[0]) ) \
                { \
                    bcf_double_set_missing((atok)->values[i]); \
                    continue; \
                } \
                has_values = 1; \
                (atok)->values[i] = (atok)->values[i] AOP (btok)->values[0]; \
            } \
        } \
    } \
    if ( !has_values ) { (atok)->nvalues = 0; (atok)->nsamples = 0; } \
}

static int vector_logic_and(token_t *atok, token_t *btok, int and_type)
{
    // We are comparing either two scalars (result of INFO tag vs a threshold), two vectors (two FORMAT fields),
    // or a vector and a scalar (FORMAT field vs threshold)
    int i, pass_site = 0;
    if ( !atok->nvalues || !btok->nvalues )
    {
        atok->nvalues = atok->nsamples = 0;
        return 0;
    }
    if ( !atok->nsamples && !btok->nsamples ) return atok->pass_site && btok->pass_site;
    if ( atok->nsamples && btok->nsamples )
    {
        if ( and_type==TOK_AND )
        {
            // perform AND within a sample
            for (i=0; i<atok->nsamples; i++)
            {
                atok->pass_samples[i] = atok->pass_samples[i] && btok->pass_samples[i];
                if ( !pass_site && atok->pass_samples[i] ) pass_site = 1;
            }
        }
        else
        {
            // perform AND across samples
            int pass_a = 0, pass_b = 0;
            for (i=0; i<atok->nsamples; i++)
            {
                if ( atok->pass_samples[i] ) pass_a = 1;
                atok->pass_samples[i] = atok->pass_samples[i] && btok->pass_samples[i];
            }
            for (i=0; i<btok->nsamples; i++)
            {
                if ( btok->pass_samples[i] ) { pass_b = 1; break; }
            }
            pass_site = pass_a && pass_b;
        }
        return pass_site;
    }
    if ( btok->nsamples )
    {
        for (i=0; i<btok->nsamples; i++)
        {
            atok->pass_samples[i] = atok->pass_site && btok->pass_samples[i];
            if ( !pass_site && atok->pass_samples[i] ) pass_site = 1;
        }
        atok->nsamples = btok->nsamples;
        return pass_site;
    }
    /* atok->nsamples!=0 */
    for (i=0; i<atok->nsamples; i++)
    {
        atok->pass_samples[i] = atok->pass_samples[i] && btok->pass_site;
        if ( !pass_site && atok->pass_samples[i] ) pass_site = 1;
    }
    return pass_site;
}
static int vector_logic_or(token_t *atok, token_t *btok, int or_type)
{
    int i, pass_site = 0;
    if ( !atok->nvalues && !btok->nvalues )   // missing sites in both
    {
        atok->nvalues = atok->nsamples = 0;
        return 0;
    }
    if ( !atok->nvalues ) // missing value in a
    {
        for (i=0; i<btok->nsamples; i++)
            atok->pass_samples[i] = btok->pass_samples[i];
        atok->nsamples = btok->nsamples;
        atok->nvalues  = 1;
        return btok->pass_site;
    }
    if ( !btok->nvalues ) // missing value in b
    {
        btok->nvalues = 1;
        return atok->pass_site;
    }

    if ( !atok->nsamples && !btok->nsamples ) return atok->pass_site || btok->pass_site;
    if ( !atok->nsamples )
    {
        if ( or_type==TOK_OR )
        {
            for (i=0; i<btok->nsamples; i++)
            {
                atok->pass_samples[i] = btok->pass_samples[i];
                if ( atok->pass_site || atok->pass_samples[i] ) pass_site = 1;
            }
        }
        else
        {
            for (i=0; i<btok->nsamples; i++)
            {
                atok->pass_samples[i] = atok->pass_site || btok->pass_samples[i];
                if ( atok->pass_samples[i] ) pass_site = 1;
            }
        }
        atok->nsamples = btok->nsamples;
        return pass_site;
    }
    if ( !btok->nsamples )  // vector vs site
    {
        if ( or_type==TOK_OR )
        {
            for (i=0; i<atok->nsamples; i++)
                if ( btok->pass_site || atok->pass_samples[i] ) pass_site = 1;
        }
        else
        {
            for (i=0; i<atok->nsamples; i++)
            {
                atok->pass_samples[i] = atok->pass_samples[i] || btok->pass_site;
                if ( atok->pass_samples[i] ) pass_site = 1;
            }
        }
        return pass_site;
    }
    for (i=0; i<atok->nsamples; i++)
    {
        atok->pass_samples[i] = atok->pass_samples[i] || btok->pass_samples[i];
        if ( !pass_site && atok->pass_samples[i] ) pass_site = 1;
    }
    return pass_site;
}

#define CMP_MISSING(atok,btok,CMP_OP,ret) \
{ \
    if ( (atok)->nsamples || (btok)->nsamples ) error("todo: Querying of missing values in FORMAT\n"); \
    token_t *tok = (atok)->is_missing ? (btok) : (atok); \
    (ret) = ( tok->nvalues CMP_OP 1 ) ? 0 : 1; \
    tok->nvalues = 1; \
}

#define CMP_VECTORS(atok,btok,CMP_OP,ret) \
{ \
    int i, j, has_values = 0, pass_site = 0; \
    if ( !(atok)->nvalues || !(btok)->nvalues ) { (atok)->nvalues = 0; (atok)->nsamples = 0; (ret) = 0; } \
    else \
    { \
        if ( (atok)->idx<=-2 || (btok)->idx<=-2 ) \
        { \
            /* any field can match: [*] */ \
            for (i=0; i<(atok)->nvalues; i++) \
            { \
                for (j=0; j<(btok)->nvalues; j++) \
                    if ( (atok)->values[i] CMP_OP (btok)->values[j] ) { pass_site = 1; i = (atok)->nvalues; break; } \
            } \
        } \
        else if ( (atok)->nsamples && (btok)->nsamples ) \
        { \
            for (i=0; i<(atok)->nsamples; i++) \
            { \
                has_values = 1; \
                if ( (atok)->values[i] CMP_OP (btok)->values[i] ) { (atok)->pass_samples[i] = 1; pass_site = 1; } \
                else (atok)->pass_samples[i] = 0; \
            } \
            if ( !has_values ) (atok)->nvalues = 0; \
        } \
        else if ( (atok)->nsamples ) \
        { \
            for (i=0; i<(atok)->nsamples; i++) \
            { \
                /*if ( bcf_double_is_missing((atok)->values[i]) ) { (atok)->pass_samples[i] = 0; continue; }*/ \
                has_values = 1; \
                if ( (atok)->values[i] CMP_OP (btok)->values[0] ) { (atok)->pass_samples[i] = 1; pass_site = 1; } \
                else (atok)->pass_samples[i] = 0; \
            } \
            if ( !has_values ) (atok)->nvalues = 0; \
        } \
        else if ( (btok)->nsamples ) \
        { \
            for (i=0; i<(btok)->nsamples; i++) \
            { \
                if ( bcf_double_is_missing((btok)->values[i]) ) { (atok)->pass_samples[i] = 0; continue; } \
                has_values = 1; \
                if ( (atok)->values[0] CMP_OP (btok)->values[i] ) { (atok)->pass_samples[i] = 1; pass_site = 1; } \
                else (atok)->pass_samples[i] = 0; \
            } \
            (atok)->nvalues  = (btok)->nvalues; \
            (atok)->nsamples = (btok)->nsamples; \
            if ( !has_values ) (atok)->nvalues = 0; \
        } \
        else \
        { \
            if ( (atok)->values[0] CMP_OP (btok)->values[0] ) { pass_site = 1; } \
        } \
        /*fprintf(stderr,"pass=%d\n", pass_site);*/ \
        (ret) = pass_site; \
    } \
}
static int cmp_vector_strings(token_t *atok, token_t *btok, int logic)    // logic: TOK_EQ or TOK_NE
{
    if ( !atok->str_value.l ) { return 0; }
    if ( !btok->str_value.l ) { atok->str_value.l = 0; return 0; }
    int i, pass_site = 0;
    if ( atok->nsamples && atok->nsamples==btok->nsamples )
    {
        for (i=0; i<atok->nsamples; i++)
        {
            char *astr = atok->str_value.s + i*atok->nstr1;
            char *bstr = btok->str_value.s + i*btok->nstr1;
            char *aend = astr + atok->str_value.l, *a = astr;
            while ( a<aend && *a ) a++;
            char *bend = bstr + btok->str_value.l, *b = bstr;
            while ( b<bend && *b ) b++;
            if ( a-astr != b-bstr ) atok->pass_samples[i] = 0;
            else atok->pass_samples[i] = strncmp(astr,bstr,a-astr)==0 ? 1 : 0;
            if ( logic!=TOK_EQ )
                atok->pass_samples[i] = atok->pass_samples[i] ? 0 : 1;
            pass_site |= atok->pass_samples[i];
        }
        if ( !atok->nsamples ) atok->nsamples = btok->nsamples;
    }
    else if ( !atok->nsamples && !btok->nsamples )
    {
        if ( atok->idx==-2 || btok->idx==-2 )
        {
            // any field can match: [*]
            if ( atok->idx==-2 && btok->idx==-2 )
                error("fixme: Expected at least one scalar value [%s %s %s]\n", atok->tag ? atok->tag : btok->tag, atok->str_value,btok->str_value);
            token_t *xtok, *ytok;   // xtok is scalar, ytok array
            if ( btok->idx==-2 ) { xtok = atok; ytok = btok; }
            else { xtok = btok; ytok = atok; }
            char *xstr = xtok->str_value.s, *xend = xstr + xtok->str_value.l;
            char *ystr = ytok->str_value.s, *yend = ystr + ytok->str_value.l, *y = ystr;
            while ( y<=yend )
            {
                if ( y==yend || *y==',' )
                {
                    if ( y-ystr==xend-xstr && !strncmp(xstr,ystr,xend-xstr) )
                    {
                        pass_site = 1;
                        break;
                    }
                    ystr = y+1;
                }
                y++;
            }
        }
        else
            pass_site = strcmp(atok->str_value.s,btok->str_value.s) ? 0 : 1;
        if ( logic!=TOK_EQ ) pass_site = pass_site ? 0 : 1;
    }
    else
    {
        token_t *xtok, *ytok;
        if ( !atok->nsamples ) { xtok = atok; ytok = btok; }
        else { xtok = btok; ytok = atok; }
        char *xstr = xtok->str_value.s;
        char *xend = xstr + xtok->str_value.l, *x = xstr;
        while ( x<xend && *x ) x++;
        for (i=0; i<ytok->nsamples; i++)
        {
            char *ystr = ytok->str_value.s + i*ytok->nstr1;
            char *ybeg = ystr, *yend = ystr + ytok->nstr1;
            int pass = 0;
            while ( ybeg < yend )
            {
                char *y = ybeg;
                while ( y<yend && *y && *y!=',' ) y++;
                if ( y - ybeg != x - xstr ) pass = 0;
                else pass = strncmp(xstr,ybeg,x-xstr)==0 ? 1 : 0;
                if ( logic!=TOK_EQ ) pass = pass ? 0 : 1;
                if ( pass || !*y ) break;
                ybeg = y+1;
            }
            atok->pass_samples[i] = pass;
            pass_site |= pass;
        }
        if ( !atok->nsamples )
            atok->nvalues = atok->nsamples = btok->nsamples; // is it a bug? not sure if atok->nvalues should be set
    }
    return pass_site;
}
static int regex_vector_strings(token_t *atok, token_t *btok, int negate)
{
    int i, pass_site = 0;
    if ( atok->nsamples )
    {
        for (i=0; i<atok->nsamples; i++)
        {
            char *ptr = atok->str_value.s + i*atok->nstr1;
            atok->pass_samples[i] = regexec(btok->regex, ptr, 0,NULL,0) ? 0 : 1;
            if ( negate ) atok->pass_samples[i] = atok->pass_samples[i] ? 0 : 1;
            pass_site |= atok->pass_samples[i];
        }
        return pass_site;
    }
    pass_site = regexec(btok->regex, atok->str_value.s, 0,NULL,0) ? 0 : 1;
    if ( negate ) pass_site = pass_site ? 0 : 1;
    return pass_site;
}

static void parse_tag_idx(char *tag, char *tag_idx, token_t *tok)   // tag_idx points just after "TAG["
{
    // TAG[*] .. any field
    if ( !strncmp("*]", tag_idx, 3) )
    {
        tok->idxs = (int*) malloc(sizeof(int));
        tok->idxs[0] = -1;
        tok->nidxs   = 1;
        tok->idx     = -2;
        return;
    }

    // TAG[integer] .. one field
    char *end, *beg = tag_idx;
    tok->idx = strtol(tag_idx, &end, 10);
    if ( tok->idx >= 0 && *end==']' ) return;


    // TAG[0,1] or TAG[0-2] or [1-] etc
    int i, ibeg = -1;
    while ( *beg && *beg!=']' )
    {
        int idx = strtol(beg, &end, 10);
        if ( end[0]==',' ) beg = end + 1;
        else if ( end[0]==']' ) beg = end;
        else if ( end[0]=='-' ) { beg = end + 1; ibeg = idx; continue; }
        else error("Could not parse the index: %s[%s\n", tag, tag_idx+1);
        if ( idx >= tok->nidxs )
        {
            tok->idxs = (int*) realloc(tok->idxs, sizeof(int)*(idx+1));
            memset(tok->idxs + tok->nidxs, 0, sizeof(int)*(idx - tok->nidxs + 1));
            tok->nidxs = idx + 1;
        }
        if ( ibeg>=0 )
        {
            for (i=ibeg; i<=idx; i++) tok->idxs[i] = 1;
            ibeg = -1;
        }
        tok->idxs[idx] = 1;
    }
    if ( ibeg >=0 )
    {
        if ( ibeg >= tok->nidxs )
        {
            tok->idxs = (int*) realloc(tok->idxs, sizeof(int)*(ibeg+1));
            memset(tok->idxs + tok->nidxs, 0, sizeof(int)*(ibeg - tok->nidxs + 1));
            tok->nidxs = ibeg + 1;
        }
        tok->idxs[ibeg] = -1;
    }
    tok->idx = -2;
}
static int filters_init1(filter_t *filter, char *str, int len, token_t *tok)
{
    tok->tok_type  = TOK_VAL;
    tok->hdr_id    = -1;
    tok->pass_site = -1;
    tok->idx       = -1;

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
        wordexp_t wexp;
        wordexp(tok->tag+1, &wexp, 0);
        if ( !wexp.we_wordc ) error("No such file: %s\n", tok->tag+1);
        int i, n;
        char **list = hts_readlist(wexp.we_wordv[0], 1, &n);
        if ( !list ) error("Could not read: %s\n", wexp.we_wordv[0]);
        wordfree(&wexp);
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
            return 0;
        }
        else if ( !strncasecmp(str,"TYPE",len) || !strncmp(str,"%TYPE",len) /* for backward compatibility */ )
        {
            tok->setter = filters_set_type;
            tok->tag = strdup("TYPE");
            return 0;
        }
        else if ( !strncasecmp(str,"FILTER",len) || !strncmp(str,"%FILTER",len) /* for backward compatibility */ )
        {
            tok->comparator = filters_cmp_filter;
            tok->tag = strdup("FILTER");
            filter->max_unpack |= BCF_UN_FLT;
            return 0;
        }
        else if ( !strncasecmp(str,"ID",len) || !strncasecmp(str,"%ID",len) /* for backward compatibility */ )
        {
            tok->comparator = filters_cmp_id;
            tok->tag = strdup("ID");
            return 0;
        }
        else if ( !strncasecmp(str,"POS",len) )
        {
            tok->setter = &filters_set_pos;
            tok->tag = strdup("POS");
            return 0;
        }
        else if ( !strncasecmp(str,"REF",len) )
        {
            tok->setter = &filters_set_ref_string;
            tok->is_str = 1;
            tok->tag = strdup("REF");
            return 0;
        }
        else if ( !strncasecmp(str,"ALT",len) )
        {
            tok->setter = &filters_set_alt_string;
            tok->is_str = 1;
            tok->tag = strdup("ALT");
            return 0;
        }
        else if ( !strncasecmp(str,"N_ALT",len) )
        {
            tok->setter = &filters_set_nalt;
            tok->tag = strdup("N_ALT");
            return 0;
        }
        else if ( !strncasecmp(str,"N_SAMPLES",len) )
        {
            tok->tok_type = TOK_VAL;
            tok->threshold = bcf_hdr_nsamples(filter->hdr);
            return 0;
        }
        else if ( !strncasecmp(str,"N_MISSING",len) )
        {
            tok->setter = &filters_set_nmissing;
            tok->tag = strdup("N_MISSING");
            return 0;
        }
        else if ( !strncasecmp(str,"F_MISSING",len) )
        {
            tok->setter = &filters_set_nmissing;
            tok->tag = strdup("F_MISSING");
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
        if ( is_array ) 
            parse_tag_idx(tmp.s, tmp.s+is_array, tok);
    }
    tok->hdr_id = bcf_hdr_id2int(filter->hdr,BCF_DT_ID,tmp.s);
    if ( is_fmt==-1 )
    {
        if ( tok->hdr_id >=0 )
        {
            if ( bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_INFO,tok->hdr_id) ) is_fmt = 0;
            else if ( bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_FMT,tok->hdr_id) ) is_fmt = 1;
        }
        if ( is_fmt==-1 ) is_fmt = 0;
    }
    tok->type = is_fmt ? BCF_HL_FMT : BCF_HL_INFO;
    if ( is_fmt ) filter->max_unpack |= BCF_UN_FMT;
    if ( tok->hdr_id>=0 )
    {
        if ( is_fmt && !strcmp("GT",tmp.s) )
        {
            tok->setter = &filters_set_genotype_string; tok->is_str = 1;
        }
        else if ( is_fmt )
        {
            if ( !bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_FMT,tok->hdr_id) )
                error("No such FORMAT field: %s\n", tmp.s);
            if ( bcf_hdr_id2number(filter->hdr,BCF_HL_FMT,tok->hdr_id)!=1 && !is_array )
                error("Error: FORMAT vectors must be subscripted, e.g. %s[0] or %s[*]\n", tmp.s, tmp.s);
            switch ( bcf_hdr_id2type(filter->hdr,BCF_HL_FMT,tok->hdr_id) )
            {
                case BCF_HT_INT:  tok->setter = &filters_set_format_int; break;
                case BCF_HT_REAL: tok->setter = &filters_set_format_float; break;
                case BCF_HT_STR:  tok->setter = &filters_set_format_string; tok->is_str = 1; break;
                default: error("[%s:%d %s] FIXME\n", __FILE__,__LINE__,__FUNCTION__);
            }
        }
        else if ( !bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_INFO,tok->hdr_id) )
            error("No such INFO field: %s\n", tmp.s);
        else
        {
            if ( bcf_hdr_id2type(filter->hdr,BCF_HL_INFO,tok->hdr_id) == BCF_HT_FLAG )
                tok->setter = filters_set_info_flag;
            else
            {
                if ( bcf_hdr_id2type(filter->hdr,BCF_HL_INFO,tok->hdr_id) == BCF_HT_STR ) tok->is_str = 1;
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
        return 0;
    }
    else if ( !strcasecmp(tmp.s,"ALT") )
    {
        tok->setter = &filters_set_alt_string;
        tok->is_str = 1;
        tok->tag = strdup(tmp.s);
        free(tmp.s);
        return 0;
    }
    else if ( !strcasecmp(tmp.s,"AN") )
    {
        tok->setter = &filters_set_an;
        tok->tag = strdup("AN");
        free(tmp.s);
        return 0;
    }
    else if ( !strcasecmp(tmp.s,"AC") )
    {
        tok->setter = &filters_set_ac;
        tok->tag = strdup("AC");
        free(tmp.s);
        return 0;
    }
    else if ( !strcasecmp(tmp.s,"MAC") )
    {
        tok->setter = &filters_set_mac;
        tok->tag = strdup("MAC");
        free(tmp.s);
        return 0;
    }
    else if ( !strcasecmp(tmp.s,"AF") )
    {
        tok->setter = &filters_set_af;
        tok->tag = strdup("AF");
        free(tmp.s);
        return 0;
    }
    else if ( !strcasecmp(tmp.s,"MAF") )
    {
        tok->setter = &filters_set_maf;
        tok->tag = strdup("MAF");
        free(tmp.s);
        return 0;
    }

    // is it a value? Here we parse as integer/float separately and use strtof
    // rather than strtod, because the more accurate double representation
    // would invalidate floating point comparisons like QUAL=59.2, obtained via
    // htslib/vcf parser
    char *end;
    tok->threshold = strtol(tmp.s, &end, 10);   // integer?
    if ( end - tmp.s != strlen(tmp.s) )
    {
        errno = 0;
        tok->threshold = strtof(tmp.s, &end);   // float?
        if ( errno!=0 || end!=tmp.s+len ) error("[%s:%d %s] Error: the tag \"INFO/%s\" is not defined in the VCF header\n", __FILE__,__LINE__,__FUNCTION__,tmp.s);
    }

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


// Parse filter expression and convert to reverse polish notation. Dijkstra's shunting-yard algorithm
filter_t *filter_init(bcf_hdr_t *hdr, const char *str)
{
    filter_t *filter = (filter_t *) calloc(1,sizeof(filter_t));
    filter->str = strdup(str);
    filter->hdr = hdr;
    filter->max_unpack |= BCF_UN_STR;

    int nops = 0, mops = 0, *ops = NULL;    // operators stack
    int nout = 0, mout = 0;                 // filter tokens, RPN
    token_t *out = NULL;
    char *tmp = filter->str;
    int last_op = -1;
    while ( *tmp )
    {
        int len, ret;
        ret = filters_next_token(&tmp, &len);
        if ( ret==-1 ) error("Missing quotes in: %s\n", str);

        // fprintf(stderr,"token=[%c] .. [%s] %d\n", TOKEN_STRING[ret], tmp, len);
        // int i; for (i=0; i<nops; i++) fprintf(stderr," .%c.", TOKEN_STRING[ops[i]]); fprintf(stderr,"\n");

        if ( ret==TOK_LFT )         // left bracket
        {
            nops++;
            hts_expand(int, nops, mops, ops);
            ops[nops-1] = ret;
        }
        else if ( ret==TOK_RGT )    // right bracket
        {
            while ( nops>0 && ops[nops-1]!=TOK_LFT )
            {
                nout++;
                hts_expand0(token_t, nout, mout, out);
                out[nout-1].tok_type = ops[nops-1];
                nops--;
            }
            if ( nops<=0 ) error("Could not parse: %s\n", str);
            nops--;
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
                ret = TOK_MULT;
            }
            else
            {
                while ( nops>0 && op_prec[ret] < op_prec[ops[nops-1]] )
                {
                    nout++;
                    hts_expand0(token_t, nout, mout, out);
                    out[nout-1].tok_type = ops[nops-1];
                    nops--;
                }
            }
            nops++;
            hts_expand(int, nops, mops, ops);
            ops[nops-1] = ret;
        }
        else if ( !len )
        {
            if ( *tmp && !isspace(*tmp) ) error("Could not parse the expression: [%s]\n", str);
            break;     // all tokens read
        }
        else           // annotation name or filtering value
        {
            nout++;
            hts_expand0(token_t, nout, mout, out);
            filters_init1(filter, tmp, len, &out[nout-1]);
            tmp += len;
        }
        last_op = ret;
    }
    while ( nops>0 )
    {
        if ( ops[nops-1]==TOK_LFT || ops[nops-1]==TOK_RGT ) error("Could not parse the expression: [%s]\n", filter->str);
        nout++;
        hts_expand0(token_t, nout, mout, out);
        out[nout-1].tok_type = ops[nops-1];
        nops--;
    }

    // In the special cases of TYPE and FILTER the BCF header IDs are yet unknown. Walk through the
    // list of operators and convert the strings (e.g. "PASS") to BCF ids. The string value token must be
    // just before or after the FILTER token and they must be followed with a comparison operator.
    // At this point we also initialize regex expressions which, in RPN, must preceed the LIKE/NLIKE operator.
    // Additionally, treat "." as missing value rather than a string in numeric equalities.
    // This code is fragile: improve me.
    int i;
    for (i=0; i<nout; i++)
    {
        if ( out[i].tok_type==TOK_EQ || out[i].tok_type==TOK_NE )
        {
            // Look for j="." and k numeric type
            int j = i-1, k = i-2;
            if ( !out[j].is_str ) { k = i-1, j = i-2; }
            if ( out[k].hdr_id>0 && out[j].is_str && out[j].key && !strcmp(".",out[j].key) )
            {
                int type = bcf_hdr_id2type(filter->hdr,out[k].type,out[k].hdr_id);
                if ( type==BCF_HT_INT ) { out[j].is_str = 0; out[j].is_missing = 1; bcf_double_set_missing(out[j].values[0]); }
                if ( type==BCF_HT_REAL ) { out[j].is_str = 0; out[j].is_missing = 1; bcf_double_set_missing(out[j].values[0]); }
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
        if ( out[i].tok_type!=TOK_VAL ) continue;
        if ( !out[i].tag ) continue;
        if ( !strcmp(out[i].tag,"TYPE") )
        {
            if ( i+1==nout ) error("Could not parse the expression: %s\n", filter->str);
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
            else if ( !strcasecmp(out[ival].key,"ref") ) { out[ival].threshold = 1; out[ival].is_str = 0; }
            else error("The type \"%s\" not recognised: %s\n", out[ival].key, filter->str);
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

            // assign correct setters and unify expressions, eg ar->ra, HOM->hom, etc
            if ( !strcasecmp(out[ival].key,"hom") ) { out[i].setter = filters_set_genotype3;  str_to_lower(out[ival].key); }
            else if ( !strcasecmp(out[ival].key,"het") ) { out[i].setter = filters_set_genotype3;  str_to_lower(out[ival].key); }
            else if ( !strcasecmp(out[ival].key,"hap") ) { out[i].setter = filters_set_genotype3;  str_to_lower(out[ival].key); }
            else if ( !strcasecmp(out[ival].key,"rr") ) { out[i].setter = filters_set_genotype2;  str_to_lower(out[ival].key); }
            else if ( !strcasecmp(out[ival].key,"ra") || !strcasecmp(out[ival].key,"ar") ) { out[i].setter = filters_set_genotype2; out[ival].key[0]='r'; out[ival].key[1]='a'; }   // ra
            else if ( !strcmp(out[ival].key,"aA") || !strcmp(out[ival].key,"Aa") ) { out[i].setter = filters_set_genotype2; out[ival].key[0]='a'; out[ival].key[1]='A'; }   // aA
            else if ( !strcasecmp(out[ival].key,"aa") ) { out[i].setter = filters_set_genotype2; out[ival].key[0]='a'; out[ival].key[1]='a'; }  // aa
            else if ( !strcasecmp(out[ival].key,"a") ) { out[i].setter = filters_set_genotype2; out[ival].key[0]='a'; out[ival].key[1]=0; }  // a
            else if ( !strcasecmp(out[ival].key,"r") ) { out[i].setter = filters_set_genotype2; out[ival].key[0]='r'; out[ival].key[1]=0; }  // r
            continue;
        }
        if ( !strcmp(out[i].tag,"FILTER") )
        {
            if ( i+1==nout ) error("Could not parse the expression: %s\n", filter->str);
            int itok = i, ival;
            if ( out[i+1].tok_type==TOK_EQ || out[i+1].tok_type==TOK_NE ) ival = i - 1;
            else if ( out[i+1].tok_type==TOK_LIKE ) out[i+1].tok_type = TOK_EQ, ival = i - 1;
            else if ( out[i+1].tok_type==TOK_NLIKE ) out[i+1].tok_type = TOK_NE, ival = i - 1;
            else if ( out[i+2].tok_type==TOK_EQ || out[i+2].tok_type==TOK_NE ) ival = ++i;
            else if ( out[i+2].tok_type==TOK_LIKE ) out[i+2].tok_type = TOK_EQ, ival = ++i;
            else if ( out[i+2].tok_type==TOK_NLIKE ) out[i+2].tok_type = TOK_NE, ival = ++i;
            else error("[%s:%d %s] Could not parse the expression: %s\n",  __FILE__,__LINE__,__FUNCTION__, filter->str);
            if ( out[ival].tok_type!=TOK_VAL || !out[ival].key )
                error("[%s:%d %s] Could not parse the expression, an unquoted string value perhaps? %s\n", __FILE__,__LINE__,__FUNCTION__, filter->str);
            if ( strcmp(".",out[ival].key) )
            {
                out[ival].hdr_id = bcf_hdr_id2int(filter->hdr, BCF_DT_ID, out[ival].key);
                if ( !bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_FLT,out[ival].hdr_id) )
                    error("The filter \"%s\" not present in the VCF header\n", out[ival].key);
            }
            else
                out[ival].hdr_id = -1;
            out[ival].tag = out[ival].key; out[ival].key = NULL;
            out[itok].hdr_id = out[ival].hdr_id;
            continue;
        }
    }
    filter->nsamples = filter->max_unpack&BCF_UN_FMT ? bcf_hdr_nsamples(filter->hdr) : 0;
    for (i=0; i<nout; i++)
    {
        if ( out[i].tok_type==TOK_MAX )      { out[i].setter = set_max; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_MIN ) { out[i].setter = set_min; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_AVG ) { out[i].setter = set_avg; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_SUM ) { out[i].setter = set_sum; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_ABS ) { out[i].setter = set_abs; out[i].tok_type = TOK_FUNC; }
        else if ( out[i].tok_type==TOK_LEN ) { out[i].setter = set_strlen; out[i].tok_type = TOK_FUNC; }
        hts_expand0(double,1,out[i].mvalues,out[i].values);
        if ( filter->nsamples )
        {
            out[i].pass_samples = (uint8_t*)malloc(filter->nsamples);
            int j;
            for (j=0; j<filter->nsamples; j++) out[i].pass_samples[j] = 1;
        }
    }

    if (0) filter_debug_print(out, NULL, nout);

    if ( mops ) free(ops);
    filter->filters   = out;
    filter->nfilters  = nout;
    filter->flt_stack = (token_t **)malloc(sizeof(token_t*)*nout);
    return filter;
}

void filter_destroy(filter_t *filter)
{
    int i;
    for (i=0; i<filter->nfilters; i++)
    {
        if ( filter->filters[i].key ) free(filter->filters[i].key);
        free(filter->filters[i].str_value.s);
        free(filter->filters[i].tag);
        free(filter->filters[i].idxs);
        free(filter->filters[i].values);
        free(filter->filters[i].pass_samples);
        if (filter->filters[i].hash) khash_str2int_destroy_free(filter->filters[i].hash);
        if (filter->filters[i].regex)
        {
            regfree(filter->filters[i].regex);
            free(filter->filters[i].regex);
        }
    }
    free(filter->filters);
    free(filter->flt_stack);
    free(filter->str);
    free(filter->tmpi);
    free(filter->tmpf);
    free(filter->tmps.s);
    free(filter);
}

int filter_test(filter_t *filter, bcf1_t *line, const uint8_t **samples)
{
    bcf_unpack(line, filter->max_unpack);

    int i, nstack = 0;
    for (i=0; i<filter->nfilters; i++)
    {
        filter->filters[i].nsamples  = 0;
        filter->filters[i].nvalues   = 0;
        filter->filters[i].pass_site = -1;

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
            else    // numeric constant
            {
                filter->filters[i].values[0] = filter->filters[i].threshold;
                filter->filters[i].nvalues   = 1;
            }
            filter->flt_stack[nstack++] = &filter->filters[i];
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_FUNC ) // all functions take only one argument
        {
            filter->filters[i].setter(filter, line, filter->flt_stack[nstack-1]);
            continue;
        }
        if ( nstack<2 )
            error("Error occurred while processing the filter \"%s\" (1:%d)\n", filter->str,nstack);  // too few values left on the stack

        int is_str  = filter->flt_stack[nstack-1]->is_str + filter->flt_stack[nstack-2]->is_str;

        if ( filter->filters[i].tok_type == TOK_OR || filter->filters[i].tok_type == TOK_OR_VEC )
        {
            if ( filter->flt_stack[nstack-1]->pass_site<0 || filter->flt_stack[nstack-2]->pass_site<0 )
                error("Error occurred while processing the filter \"%s\" (%d %d OR)\n", filter->str,filter->flt_stack[nstack-2]->pass_site,filter->flt_stack[nstack-1]->pass_site);
            filter->flt_stack[nstack-2]->pass_site = vector_logic_or(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1], filter->filters[i].tok_type);
            nstack--;
            continue;
        }
        if ( filter->filters[i].tok_type == TOK_AND || filter->filters[i].tok_type == TOK_AND_VEC )
        {
            if ( filter->flt_stack[nstack-1]->pass_site<0 || filter->flt_stack[nstack-2]->pass_site<0 )
                error("Error occurred while processing the filter \"%s\" (%d %d AND)\n", filter->str,filter->flt_stack[nstack-2]->pass_site,filter->flt_stack[nstack-1]->pass_site);
            filter->flt_stack[nstack-2]->pass_site = vector_logic_and(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1], filter->filters[i].tok_type);
            nstack--;
            continue;
        }

        if ( filter->filters[i].tok_type == TOK_ADD )
        {
            VECTOR_ARITHMETICS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],+);
            nstack--;
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_SUB )
        {
            VECTOR_ARITHMETICS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],-);
            nstack--;
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_MULT )
        {
            VECTOR_ARITHMETICS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],*);
            nstack--;
            continue;
        }
        else if ( filter->filters[i].tok_type == TOK_DIV )
        {
            VECTOR_ARITHMETICS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],/);
            nstack--;
            continue;
        }

        int is_true = 0;
        if ( filter->filters[i].comparator )
            is_true = filter->filters[i].comparator(filter->flt_stack[nstack-1],filter->flt_stack[nstack-2],filter->filters[i].tok_type,line);
        else if ( !filter->flt_stack[nstack-1]->nvalues || !filter->flt_stack[nstack-2]->nvalues )
        {
            int skip = 0;
            if ( !filter->flt_stack[nstack-2]->is_missing && !filter->flt_stack[nstack-1]->is_missing ) skip = 1;
            if ( filter->filters[i].tok_type != TOK_EQ  && filter->filters[i].tok_type != TOK_NE ) skip = 1;

            if ( skip ) 
                filter->flt_stack[nstack-2]->nvalues = filter->flt_stack[nstack-2]->nsamples = 0;
            else if ( filter->filters[i].tok_type == TOK_EQ )
                CMP_MISSING(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],==,is_true)
            else if ( filter->filters[i].tok_type == TOK_NE )
                CMP_MISSING(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],!=,is_true)
        }
        else if ( filter->filters[i].tok_type == TOK_EQ )
        {
            if ( filter->flt_stack[nstack-1]->comparator )
                is_true = filter->flt_stack[nstack-1]->comparator(filter->flt_stack[nstack-1],filter->flt_stack[nstack-2],TOK_EQ,line);
            else if ( filter->flt_stack[nstack-2]->comparator )
                is_true = filter->flt_stack[nstack-2]->comparator(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],TOK_EQ,line);
            else if ( is_str==2 )   // both are strings
                is_true = cmp_vector_strings(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],TOK_EQ);
            else if ( is_str==1 )
                error("Comparing string to numeric value: %s\n", filter->str);
            else
                CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],==,is_true);
        }
        else if ( filter->filters[i].tok_type == TOK_NE )
        {
            if ( filter->flt_stack[nstack-1]->comparator )
                is_true = filter->flt_stack[nstack-1]->comparator(filter->flt_stack[nstack-1],filter->flt_stack[nstack-2],TOK_NE,line);
            else if ( filter->flt_stack[nstack-2]->comparator )
                is_true = filter->flt_stack[nstack-2]->comparator(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],TOK_NE,line);
            else if ( is_str==2 )
                is_true = cmp_vector_strings(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],TOK_NE);
            else if ( is_str==1 )
                error("Comparing string to numeric value: %s\n", filter->str);
            else
                CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],!=,is_true);
        }
        else if ( filter->filters[i].tok_type == TOK_LIKE || filter->filters[i].tok_type == TOK_NLIKE )
        {
            if ( is_str==2 )
                is_true = regex_vector_strings(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1], filter->filters[i].tok_type == TOK_LIKE ? 0 : 1);
            else
                error("The regex operator can be used on strings only: %s\n", filter->str);
        }
        else if ( is_str>0 )
            error("Wrong operator in string comparison: %s [%s,%s]\n", filter->str, filter->flt_stack[nstack-1]->str_value, filter->flt_stack[nstack-2]->str_value);
        else if ( filter->filters[i].tok_type == TOK_LE )
            CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],<=,is_true)
        else if ( filter->filters[i].tok_type == TOK_LT )
            CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],<,is_true)
        else if ( filter->filters[i].tok_type == TOK_BT )
            CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],>,is_true)
        else if ( filter->filters[i].tok_type == TOK_BE )
            CMP_VECTORS(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1],>=,is_true)
        else
            error("FIXME: did not expect this .. tok_type %d = %d\n", i, filter->filters[i].tok_type);

        filter->flt_stack[nstack-2]->pass_site = is_true;
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
