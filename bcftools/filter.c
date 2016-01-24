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
#include <errno.h>
#include <math.h>
#include <wordexp.h>
#include <regex.h>
#include <htslib/khash_str2int.h>
#include "filter.h"
#include "bcftools.h"
#include <htslib/hts_defs.h>
#include <htslib/vcfutils.h>

typedef struct _token_t
{
    // read-only values, same for all VCF lines
    int tok_type;       // one of the TOK_* keys below
    char *key;          // set only for string constants, otherwise NULL
    char *tag;          // for debugging and printout only, VCF tag name
    float threshold;    // filtering threshold
    int hdr_id, type;   // BCF header lookup ID and one of BCF_HT_* types
    int idx;            // 0-based index to VCF vectors, -1: not a vector, -2: any field ([*])
    void (*setter)(filter_t *, bcf1_t *, struct _token_t *);
    int (*comparator)(struct _token_t *, struct _token_t *, int op_type, bcf1_t *);
    void *hash;         // test presence of str value in the hash via comparator
    regex_t *regex;     // precompiled regex for string comparison

    // modified on filter evaluation at each VCF line
    float *values;      // In case str_value is set, values[0] is one sample's string length
    char *str_value;    //  and values[0]*nsamples gives the total length;
    int is_str, is_missing; // is_missing is set only for constants, variables are controled via nvalues
    int pass_site;          // -1 not applicable, 0 fails, >0 pass
    uint8_t *pass_samples;  // status of individual samples
    int nsamples;           // number of samples
    int nvalues, mvalues;   // number of used values, n=0 for missing values, n=1 for scalars
                            // for strings, total length of str_value
}
token_t;

struct _filter_t
{
    bcf_hdr_t *hdr;
    char *str;
    int nfilters;
    token_t *filters, **flt_stack;  // filtering input tokens (in RPN) and evaluation stack
    int32_t *tmpi;
    int max_unpack, mtmpi, nsamples;
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
        // hacky: so that [*] is not split, the tokenizer does not recognise square brackets []
        if ( tmp[0]=='*' && (tmp==*str || tmp[-1]!='[') ) break;
        if ( tmp[0]=='-' ) break;
        if ( tmp[0]=='/' ) break;
        if ( tmp[0]=='~' ) break;
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
        tok->values[0] = line->qual;
        tok->nvalues = 1;
    }
}
static void filters_set_type(filter_t *flt, bcf1_t *line, token_t *tok)
{
    tok->values[0] = bcf_get_variant_types(line);
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
        int m = (int)tok->values[0];
        hts_expand(char,n+1,m,tok->str_value);
        memcpy(tok->str_value,line->d.info[i].vptr,n);
        tok->str_value[n] = 0;
        tok->values[0] = m;
        tok->nvalues   = n;
    }
    else if ( line->d.info[i].type==BCF_BT_FLOAT )
    {
        if ( bcf_float_is_missing(line->d.info[i].v1.f) ) tok->nvalues = 0;
        else
        {
            tok->values[0] = line->d.info[i].v1.f;
            tok->nvalues   = 1;
        }
        tok->str_value = NULL;
    }
    else
    {
        if ( line->d.info[i].type==BCF_BT_INT8 && line->d.info[i].v1.i==bcf_int8_missing ) tok->nvalues = 0;
        else if ( line->d.info[i].type==BCF_BT_INT16 && line->d.info[i].v1.i==bcf_int16_missing ) tok->nvalues = 0;
        else if ( line->d.info[i].type==BCF_BT_INT32 && line->d.info[i].v1.i==bcf_int32_missing ) tok->nvalues = 0;
        else
        {
            tok->values[0] = line->d.info[i].v1.i;
            tok->nvalues   = 1;
        }
        tok->str_value = NULL;
    }
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

    if ( op_type==TOK_EQ ) return strcmp(btok->str_value,line->d.id) ? 0 : 1;
    return strcmp(btok->str_value,line->d.id) ? 1 : 0;
}

/**
 *  bcf_get_info_value() - get single INFO value, int or float
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
        if ( info->type==BCF_BT_FLOAT ) *((float*)value) = info->v1.f;
        else if ( info->type==BCF_BT_INT8 || info->type==BCF_BT_INT16 || info->type==BCF_BT_INT32 ) *((int*)value) = info->v1.i;
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
        case BCF_BT_INT8:  BRANCH(int8_t,  p[j]==bcf_int8_missing,  p[j]==bcf_int8_vector_end,  int); break;
        case BCF_BT_INT16: BRANCH(int16_t, p[j]==bcf_int16_missing, p[j]==bcf_int16_vector_end, int); break;
        case BCF_BT_INT32: BRANCH(int32_t, p[j]==bcf_int32_missing, p[j]==bcf_int32_vector_end, int); break;
        case BCF_BT_FLOAT: BRANCH(float,   bcf_float_is_missing(p[j]), bcf_float_is_vector_end(p[j]), float); break;
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
        int i, n = bcf_get_info_int32(flt->hdr,line,tok->tag,&flt->tmpi,&flt->mtmpi);
        tok->nvalues = n;
        hts_expand(float,n,tok->mvalues,tok->values);
        for (i=0; i<n; i++) tok->values[i] = flt->tmpi[i];
    }
    else
    {
        int32_t value;
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
        tok->nvalues = bcf_get_info_float(flt->hdr,line,tok->tag,&tok->values,&tok->mvalues);
        if ( tok->nvalues<0 ) tok->nvalues = 0;
    }
    else
    {
        float value;
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
    int m = (int)tok->values[0];
    int n = bcf_get_info_string(flt->hdr,line,tok->tag,&tok->str_value,&m);
    if ( n<0 ) { tok->nvalues = 0; return; }
    tok->values[0] = m;     // allocated length

    if ( tok->idx>=0 )
    {
        // get ith field (i=tok->idx)
        int i = 0;
        char *ss = tok->str_value, *se = tok->str_value + n;
        while ( ss<se && i<tok->idx )
        {
            if ( *ss==',' ) i++;
            ss++;
        }
        if ( ss==se || i!=tok->idx ) { tok->nvalues = 0; return; }
        se = ss;
        while ( se-tok->str_value<n && *se!=',' ) se++;
        if ( ss==tok->str_value ) *se = 0;
        else
        {
            memmove(tok->str_value,ss,se-ss);
            tok->str_value[se-ss] = 0;
        }
        tok->nvalues = se-ss;
    }
    else if ( tok->idx==-2 ) tok->nvalues = n;
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
        tok->nvalues = 0;
    else
    {
        int is_missing = 1;
        hts_expand(float,tok->nvalues,tok->mvalues,tok->values);
        for (i=0; i<tok->nvalues; i++)
        {
            if ( flt->tmpi[i]==bcf_int32_missing || flt->tmpi[i]==bcf_int32_vector_end )
                bcf_float_set_missing(tok->values[i]);
            else
            {
                tok->values[i] = flt->tmpi[i];
                is_missing = 0;
            }
        }
        if ( is_missing ) tok->nvalues = 0;
        else if ( tok->idx >= 0 )
        {
            int nsmpl = bcf_hdr_nsamples(flt->hdr);
            int nvals = tok->nvalues / nsmpl;
            if ( tok->idx >= nvals )
                tok->nvalues = 0;  // the index is too big
            else
            {
                for (i=0; i<nsmpl; i++)
                    tok->values[i] = tok->values[i*nvals+tok->idx];
                tok->nvalues = nsmpl;
            }
        }
    }
    tok->nsamples = tok->nvalues;
}
static void filters_set_format_float(filter_t *flt, bcf1_t *line, token_t *tok)
{
    if ( (tok->nvalues=bcf_get_format_float(flt->hdr,line,tok->tag,&tok->values,&tok->mvalues))<=0 )
        tok->nvalues = tok->nsamples = 0;   // missing values
    else if ( tok->idx >= 0 )
    {
        int i, nsmpl, nvals;
        nsmpl = bcf_hdr_nsamples(flt->hdr);
        nvals = tok->nvalues / nsmpl;
        if ( tok->idx >= nvals )
            tok->nsamples = tok->nvalues = 0;  // the index is too big
        else
        {
            for (i=0; i<nsmpl; i++)
                tok->values[i] = tok->values[i*nvals+tok->idx];
            tok->nsamples = tok->nvalues = nsmpl;
        }
    }
    tok->nsamples = tok->nvalues;
}
static void filters_set_format_string(filter_t *flt, bcf1_t *line, token_t *tok)
{
    int ndim = tok->nsamples * (int)tok->values[0];
    int ret = bcf_get_format_char(flt->hdr,line,tok->tag,&tok->str_value,&ndim);

    int nsmpl = bcf_hdr_nsamples(flt->hdr);
    ndim /= nsmpl;
    tok->values[0] = ndim;

    if ( ret<=0 )
    {
        tok->nvalues = 0;
        return;
    }

    if ( tok->idx < 0 ) // scalar
    {
        tok->nvalues = tok->nsamples = nsmpl;
        return;
    }

    // vector
    int i;
    for (i=0; i<nsmpl; i++)
    {
        char *ss = tok->str_value + i*ndim;
        int is = 0, ivec = 0;
        while ( ivec<tok->idx && is<ndim && ss[is] )
        {
            if ( ss[is]==',' ) ivec++;
            is++;
        }
        if ( ivec!=tok->idx || is==ndim || !ss[is] )
        {
            ss[0] = '.';
            ss[1] = 0;
            continue;
        }
        int ie = is;
        while ( ie<ndim && ss[ie] && ss[ie]!=',' ) ie++;
        if ( is ) memmove(ss,&ss[is],ie-is);
        if ( ndim-(ie-is) ) memset(ss+ie-is,0,ndim-(ie-is));
    }
    if ( !ndim )
    {
        tok->nvalues = 0;
        return;
    }
    tok->nvalues  = ret;
    tok->nsamples = nsmpl;
}
static void filters_set_genotype_string(filter_t *flt, bcf1_t *line, token_t *tok)
{
    bcf_fmt_t *fmt = bcf_get_fmt(flt->hdr, line, "GT");
    if ( !fmt )
    {
        tok->nvalues = tok->nsamples = 0;
        return;
    }
    int i, blen = 3, nsmpl = bcf_hdr_nsamples(flt->hdr);
    kstring_t str;

gt_length_too_big:
    str.s = tok->str_value; str.m = tok->values[0] * nsmpl; str.l = 0;
    for (i=0; i<nsmpl; i++)
    {
        int plen = str.l;

        #define BRANCH(type_t) { \
            type_t *ptr = (type_t*) (fmt->p + i*fmt->size); \
            if ( !(ptr[0]>>1) ) kputc('.',&str); \
        }
        switch (fmt->type) {
            case BCF_BT_INT8:  BRANCH(int8_t); break;
            case BCF_BT_INT16: BRANCH(int16_t); break;
            case BCF_BT_INT32: BRANCH(int32_t); break;
            default: fprintf(stderr,"FIXME: type %d in bcf_format_gt?\n", fmt->type); abort(); break;
        }
        #undef BRANCH

        if ( plen==str.l )
        {
            bcf_format_gt(fmt, i, &str);
            if ( str.l - plen > blen )
            {
                // too many alternate alleles or ploidy is too large, the genotype does not fit
                // three characters ("0/0" vs "10/10").
                tok->str_value = str.s;
                blen *= 2;
                goto gt_length_too_big;
            }
        }

        plen = str.l - plen;
        while ( plen<blen )
        {
            kputc_(0, &str);
            plen++;
        }
    }
    tok->nvalues = str.l;
    tok->nsamples = nsmpl;
    tok->values[0] = blen;
    tok->str_value = str.s;
}
static void filters_set_ref_string(filter_t *flt, bcf1_t *line, token_t *tok)
{
    kstring_t str; str.s = tok->str_value; str.m = tok->values[0]; str.l = 0;
    kputs(line->d.allele[0], &str);
    tok->nvalues = str.l;
    tok->values[0] = str.m;
    tok->str_value = str.s;
}
static void filters_set_alt_string(filter_t *flt, bcf1_t *line, token_t *tok)
{
    kstring_t str; str.s = tok->str_value; str.m = tok->values[0]; str.l = 0;
    if ( tok->idx>=0 )
    {
        if ( line->n_allele >= tok->idx )
            kputs(line->d.allele[tok->idx], &str);
        else
            kputc('.', &str);
    }
    else if ( line->n_allele>1 )
    {
        kputs(line->d.allele[1], &str);
        int i;
        for (i=2; i<line->n_allele; i++)
        {
            kputc(',', &str);
            kputs(line->d.allele[i], &str);
        }
    }
    else if ( line->n_allele==1 )
        kputc('.', &str);
    tok->nvalues = str.l;
    tok->values[0] = str.m;
    tok->str_value = str.s;
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
        hts_expand(float,line->n_allele,tok->mvalues,tok->values);
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
        tok->values[i] /= (float)an;
}
static void filters_set_maf(filter_t *flt, bcf1_t *line, token_t *tok)
{
    filters_set_ac(flt,line,tok);
    if ( !tok->nvalues ) return;
    int i, an = flt->tmpi[0];
    for (i=0; i<tok->nvalues; i++)
    {
        tok->values[i] /= (float)an;
        if ( tok->values[i] > 0.5 ) tok->values[i] = 1 - tok->values[i];
    }
}

static void set_max(filter_t *flt, bcf1_t *line, token_t *tok)
{
    float val = -HUGE_VAL;
    int i;
    for (i=0; i<tok->nvalues; i++)
    {
        if ( !bcf_float_is_missing(tok->values[i]) && val < tok->values[i] ) val = tok->values[i];
    }
    tok->values[0] = val;
    tok->nvalues   = 1;
    tok->nsamples  = 0;
}
static void set_min(filter_t *flt, bcf1_t *line, token_t *tok)
{
    float val = HUGE_VAL;
    int i;
    for (i=0; i<tok->nvalues; i++)
        if ( !bcf_float_is_missing(tok->values[i]) && val > tok->values[i] ) val = tok->values[i];
    tok->values[0] = val;
    tok->nvalues   = 1;
    tok->nsamples  = 0;
}
static void set_avg(filter_t *flt, bcf1_t *line, token_t *tok)
{
    float val = 0;
    int i, n = 0;
    for (i=0; i<tok->nvalues; i++)
        if ( !bcf_float_is_missing(tok->values[i]) ) { val += tok->values[i]; n++; }
    tok->values[0] = n ? val / n : 0;
    tok->nvalues   = 1;
    tok->nsamples  = 0;
}
static void set_sum(filter_t *flt, bcf1_t *line, token_t *tok)
{
    float val = 0;
    int i, n = 0;
    for (i=0; i<tok->nvalues; i++)
        if ( !bcf_float_is_missing(tok->values[i]) ) { val += tok->values[i]; n++; }
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
    if ( !tok->nvalues ) return;
    if ( tok->idx==-2 )
    {
        int i = 0;
        char *ss = tok->str_value;
        while ( *ss )
        {
            char *se = ss;
            while ( *se && *se!=',' ) se++;
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
        tok->values[0] = strlen(tok->str_value);
        tok->nvalues = 1;
    }
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
                if ( bcf_float_is_missing((atok)->values[i]) ) continue; \
                if ( bcf_float_is_missing((btok)->values[i]) ) { bcf_float_set_missing((atok)->values[i]); continue; } \
                has_values = 1; \
                (atok)->values[i] = (atok)->values[i] AOP (btok)->values[i]; \
            } \
        } \
        else if ( (btok)->nsamples ) \
        { \
            hts_expand(float,(btok)->nvalues,(atok)->mvalues,(atok)->values); \
            for (i=0; i<(btok)->nvalues; i++) \
            { \
                if ( bcf_float_is_missing((atok)->values[0]) || bcf_float_is_missing((btok)->values[i]) ) \
                { \
                    bcf_float_set_missing((atok)->values[i]); \
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
                if ( bcf_float_is_missing((atok)->values[i]) || bcf_float_is_missing((btok)->values[0]) ) \
                { \
                    bcf_float_set_missing((atok)->values[i]); \
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
        return btok->pass_site;
    }
    if ( !btok->nvalues ) // missing value in b
        return atok->pass_site;

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
}

#define CMP_VECTORS(atok,btok,CMP_OP,ret) \
{ \
    int i, j, has_values = 0, pass_site = 0; \
    if ( !(atok)->nvalues || !(btok)->nvalues ) { (atok)->nvalues = 0; (atok)->nsamples = 0; (ret) = 0; } \
    else \
    { \
        if ( (atok)->nsamples && (btok)->nsamples ) \
        { \
            for (i=0; i<(atok)->nsamples; i++) \
            { \
                if ( bcf_float_is_missing((atok)->values[i]) ) { (atok)->pass_samples[i] = 0; continue; } \
                if ( bcf_float_is_missing((btok)->values[i]) ) { (atok)->pass_samples[i] = 0; continue; } \
                has_values = 1; \
                if ( (atok)->values[i] CMP_OP (btok)->values[i] ) { (atok)->pass_samples[i] = 1; pass_site = 1; } \
                else (atok)->pass_samples[i] = 0; \
            } \
            if ( !has_values ) (atok)->nvalues = 0; \
        } \
        else if ( (atok)->nsamples ) \
        { \
            if ( bcf_float_is_missing((btok)->values[0]) ) { (atok)->nvalues = 0; (atok)->nsamples = 0; (ret) = 0; } \
            else \
            { \
                for (i=0; i<(atok)->nsamples; i++) \
                { \
                    if ( bcf_float_is_missing((atok)->values[i]) ) { (atok)->pass_samples[i] = 0; continue; } \
                    has_values = 1; \
                    if ( (atok)->values[i] CMP_OP (btok)->values[0] ) { (atok)->pass_samples[i] = 1; pass_site = 1; } \
                    else (atok)->pass_samples[i] = 0; \
                } \
            } \
            if ( !has_values ) (atok)->nvalues = 0; \
        } \
        else if ( (btok)->nsamples ) \
        { \
            if ( bcf_float_is_missing((atok)->values[0]) ) { (atok)->nvalues = 0; (atok)->nsamples = 0; (ret) = 0; } \
            else \
            { \
                for (i=0; i<(btok)->nsamples; i++) \
                { \
                    if ( bcf_float_is_missing((btok)->values[i]) ) { (atok)->pass_samples[i] = 0; continue; } \
                    has_values = 1; \
                    if ( (atok)->values[0] CMP_OP (btok)->values[i] ) { (atok)->pass_samples[i] = 1; pass_site = 1; } \
                    else (atok)->pass_samples[i] = 0; \
                } \
                (atok)->nvalues  = (btok)->nvalues; \
                (atok)->nsamples = (btok)->nsamples; \
            } \
            if ( !has_values ) (atok)->nvalues = 0; \
        } \
        else if ( (atok)->idx==-2 || (btok)->idx==-2 ) \
        { \
            /* any field can match: [*] */ \
            for (i=0; i<(atok)->nvalues; i++) \
            { \
                for (j=0; j<(btok)->nvalues; j++) \
                    if ( (atok)->values[i] CMP_OP (btok)->values[j] ) { pass_site = 1; i = (atok)->nvalues; break; } \
            } \
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
    if ( !atok->nvalues ) { return 0; }
    if ( !btok->nvalues ) { atok->nvalues = 0; return 0; }
    int i, pass_site = 0;
    if ( atok->nsamples && atok->nsamples==btok->nsamples )
    {
        for (i=0; i<atok->nsamples; i++)
        {
            char *astr = atok->str_value + i*(int)atok->values[0];
            char *bstr = btok->str_value + i*(int)btok->values[0];
            char *aend = astr + (int)atok->values[0], *a = astr;
            while ( a<aend && *a ) a++;
            char *bend = bstr + (int)btok->values[0], *b = bstr;
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
            char *xstr = xtok->str_value, *xend = xstr + xtok->nvalues;
            char *ystr = ytok->str_value, *yend = ystr + ytok->nvalues, *y = ystr;
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
            pass_site = strcmp(atok->str_value,btok->str_value) ? 0 : 1;
        if ( logic!=TOK_EQ ) pass_site = pass_site ? 0 : 1;
    }
    else
    {
        token_t *xtok, *ytok;
        if ( !atok->nsamples ) { xtok = atok; ytok = btok; }
        else { xtok = btok; ytok = atok; }
        char *xstr = xtok->str_value;
        char *xend = xstr + (int)xtok->values[0], *x = xstr;
        while ( x<xend && *x ) x++;
        for (i=0; i<ytok->nsamples; i++)
        {
            char *ystr = ytok->str_value + i*(int)ytok->values[0];
            char *yend = ystr + (int)ytok->values[0], *y = ystr;
            while ( y<yend && *y ) y++;
            if ( x-xstr != y-ystr ) atok->pass_samples[i] = 0;
            else atok->pass_samples[i] = strncmp(xstr,ystr,x-xstr)==0 ? 1 : 0;
            if ( logic!=TOK_EQ )
                atok->pass_samples[i] = atok->pass_samples[i] ? 0 : 1;
            pass_site |= atok->pass_samples[i];
        }
        if ( !atok->nsamples )
            atok->nvalues = atok->nsamples = btok->nsamples; // is it a bug? not sure if atok->nvalues should be set
    }
    return pass_site;
}
static int regex_vector_strings(token_t *atok, token_t *btok)
{
    int ret = regexec(btok->regex, atok->str_value, 0,NULL,0);
    return ret==0 ? 1 : 0;
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
        hts_expand(float,1,tok->mvalues,tok->values);
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
        {
            if ( tmp.s[is_array]=='*' )
                tok->idx = -2;      // tag[*] .. any field
            else
            {
                char *end;
                tok->idx = strtol(tmp.s+is_array, &end, 10);
                if ( *end!=']' ) error("Could not parse the index: %s[%s\n", tmp.s,tmp.s+is_array);
            }
        }
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
                    if(!is_array) tok->idx = -2;
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

    // is it a value?
    char *end;
    errno = 0;
    tok->threshold = strtod(tmp.s, &end);
    if ( errno!=0 || end!=tmp.s+len ) error("[%s:%d %s] Error: the tag \"INFO/%s\" is not defined in the VCF header\n", __FILE__,__LINE__,__FUNCTION__,tmp.s);

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

        //fprintf(stderr,"token=[%c] .. [%s] %d\n", TOKEN_STRING[ret], tmp, len);
        //int i; for (i=0; i<nops; i++) fprintf(stderr," .%c.", TOKEN_STRING[ops[i]]); fprintf(stderr,"\n");

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
            if ( out[k].hdr_id>0 && out[j].is_str && !strcmp(".",out[j].key) )
            {
                int type = bcf_hdr_id2type(filter->hdr,out[k].type,out[k].hdr_id);
                if ( type==BCF_HT_INT ) { out[j].is_str = 0; out[j].is_missing = 1; bcf_float_set_missing(out[j].values[0]); }
                if ( type==BCF_HT_REAL ) { out[j].is_str = 0; out[j].is_missing = 1; bcf_float_set_missing(out[j].values[0]); }
            }
        }
        if ( out[i].tok_type==TOK_LIKE || out[i].tok_type==TOK_NLIKE )
        {
            int j = i-1;
            if ( !out[j].key )
                error("Could not parse the expression, wrong value for regex operator: %s\n", filter->str);
            out[j].regex = (regex_t *) malloc(sizeof(regex_t));
            if ( regcomp(out[j].regex, out[j].key, REG_ICASE|REG_NOSUB) )
                error("Could not compile the regex expression \"%s\": %s\n", out[j].key,filter->str);
        }
        if ( out[i].tok_type!=TOK_VAL ) continue;
        if ( !out[i].tag ) continue;
        if ( !strcmp(out[i].tag,"TYPE") )
        {
            if ( i+1==nout ) error("Could not parse the expression: %s\n", filter->str);
            int j = i+1;
            if ( out[j].tok_type==TOK_EQ || out[j].tok_type==TOK_NE ) j = i - 1;
            if ( out[j].tok_type!=TOK_VAL || !out[j].key ) error("[%s:%d %s] Could not parse the expression: %s\n",  __FILE__,__LINE__,__FUNCTION__, filter->str);
            if ( !strcasecmp(out[j].key,"snp") || !strcasecmp(out[j].key,"snps") ) { out[j].threshold = VCF_SNP; out[j].is_str = 0; }
            else if ( !strcasecmp(out[j].key,"indel") || !strcasecmp(out[j].key,"indels") ) { out[j].threshold = VCF_INDEL; out[j].is_str = 0; }
            else if ( !strcasecmp(out[j].key,"mnp") || !strcasecmp(out[j].key,"mnps") ) { out[j].threshold = VCF_MNP; out[j].is_str = 0; }
            else if ( !strcasecmp(out[j].key,"other") ) { out[j].threshold = VCF_OTHER; out[j].is_str = 0; }
            else if ( !strcasecmp(out[j].key,"ref") ) { out[j].threshold = VCF_REF; out[j].is_str = 0; }
            else error("The type \"%s\" not recognised: %s\n", out[j].key, filter->str);
            out[j].tag = out[j].key; out[j].key = NULL;
            i = j;
            continue;
        }
        if ( !strcmp(out[i].tag,"FILTER") )
        {
            if ( i+1==nout ) error("Could not parse the expression: %s\n", filter->str);
            int j = i+1;
            if ( out[j].tok_type==TOK_EQ || out[j].tok_type==TOK_NE ) j = i - 1;    // the expression has "value"=FILTER rather than FILTER="value"
            if ( out[j].tok_type==TOK_LIKE ) out[j].tok_type = TOK_EQ;              // for FILTER, ~ and !~ work the same way as = and !=
            if ( out[j].tok_type==TOK_NLIKE ) out[j].tok_type = TOK_NE;
            if ( out[j+1].tok_type==TOK_LIKE ) out[j+1].tok_type = TOK_EQ;
            if ( out[j+1].tok_type==TOK_NLIKE ) out[j+1].tok_type = TOK_NE;
            if ( out[j].tok_type!=TOK_VAL || !out[j].key )
                error("[%s:%d %s] Could not parse the expression, an unquoted string value perhaps? %s\n", __FILE__,__LINE__,__FUNCTION__, filter->str);
            if ( strcmp(".",out[j].key) )
            {
                out[j].hdr_id = bcf_hdr_id2int(filter->hdr, BCF_DT_ID, out[j].key);
                if ( !bcf_hdr_idinfo_exists(filter->hdr,BCF_HL_FLT,out[j].hdr_id) )
                    error("The filter \"%s\" not present in the VCF header\n", out[j].key);
            }
            else
                out[j].hdr_id = -1;
            out[j].tag = out[j].key; out[j].key = NULL;
            out[i].hdr_id = out[j].hdr_id;
            i = j;
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
        hts_expand0(float,1,out[i].mvalues,out[i].values);
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
        //if ( filter->filters[i].key ) free(filter->filters[i].key);
        free(filter->filters[i].str_value);
        free(filter->filters[i].tag);
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
                filter->filters[i].str_value = filter->filters[i].key;
                filter->filters[i].values[0] = filter->filters[i].values[0];
                filter->filters[i].nvalues   = strlen(filter->filters[i].key);
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
        if ( !filter->flt_stack[nstack-1]->nvalues || !filter->flt_stack[nstack-2]->nvalues )
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
            {
                is_true = regex_vector_strings(filter->flt_stack[nstack-2],filter->flt_stack[nstack-1]);
                if ( filter->filters[i].tok_type == TOK_NLIKE ) is_true = is_true ? 0 : 1;
            }
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
