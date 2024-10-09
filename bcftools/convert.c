/*  convert.c -- functions for converting between VCF/BCF and related formats.

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

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <htslib/kfunc.h>
#include <htslib/khash_str2int.h>
#include <htslib/hts_endian.h>
#include "bcftools.h"
#include "variantkey.h"
#include "convert.h"
#include "filter.h"

#define T_CHROM   1
#define T_POS     2
#define T_ID      3
#define T_REF     4
#define T_ALT     5
#define T_QUAL    6
#define T_FILTER  7
#define T_INFO    8
#define T_FORMAT  9
#define T_SAMPLE  10
#define T_SEP     11
#define T_IS_TS   12
#define T_TYPE    13
#define T_MASK    14
#define T_GT      15
#define T_TGT     16
#define T_LINE    17
#define T_CHROM_POS_ID 18   // not publicly advertised
#define T_GT_TO_PROB3  19   // not publicly advertised
#define T_PL_TO_PROB3  20   // not publicly advertised
#define T_GP_TO_PROB3  21   // not publicly advertised
#define T_FIRST_ALT    22   // not publicly advertised
#define T_IUPAC_GT     23
#define T_GT_TO_HAP    24   // not publicly advertised
#define T_GT_TO_HAP2   25   // not publicly advertised
#define T_TBCSQ        26
#define T_END          27
#define T_POS0         28
#define T_END0         29
#define T_RSX          30   // RSID HEX
#define T_VKX          31   // VARIANTKEY HEX
#define T_PBINOM       32
#define T_NPASS        33

typedef struct _fmt_t
{
    int type, id, is_gt_field, ready, subscript;
    char *key;
    bcf_fmt_t *fmt;
    void *usr;                  // user data (optional)
    void (*handler)(convert_t *, bcf1_t *, struct _fmt_t *, int, kstring_t *);
    void (*destroy)(void*);     // clean user data (optional)
}
fmt_t;

struct _convert_t
{
    fmt_t *fmt;
    int nfmt, mfmt;
    int nsamples, *samples;
    bcf_hdr_t *header;
    int max_unpack;
    char *format_str;
    bcf_srs_t *readers; // required only for %MASK
    int nreaders;
    void *dat;
    int ndat;
    char *undef_info_tag;
    void *used_tags_hash;
    char **used_tags_list;
    char *print_filtered;
    int nused_tags;
    int allow_undef_tags;
    int force_newline;
    int header_samples;
    int no_hdr_indices;
    uint8_t **subset_samples;
};

typedef struct
{
    kstring_t hap1,hap2;
    char **str;
    int n, m;
}
bcsq_t;

static void process_chrom(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputs(convert->header->id[BCF_DT_CTG][line->rid].key, str); }
static void process_pos(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputw(line->pos+1, str); }
static void process_pos0(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputw(line->pos, str); }
static void process_end(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputw(line->pos+line->rlen, str); }
static void process_end0(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputw(line->pos+line->rlen-1, str); }
static void process_id(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputs(line->d.id, str); }
static void process_ref(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { kputs(line->d.allele[0], str); }
static void process_alt(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int i;
    if ( line->n_allele==1 )
    {
        kputc('.', str);
        return;
    }
    if ( fmt->subscript>=0 )
    {
        if ( line->n_allele > fmt->subscript+1 )
            kputs(line->d.allele[fmt->subscript+1], str);
        else
            kputc('.', str);
        return;
    }
    for (i=1; i<line->n_allele; i++)
    {
        if ( i>1 ) kputc(',', str);
        kputs(line->d.allele[i], str);
    }
}
static void process_first_alt(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( line->n_allele==1 )
        kputc('.', str);
    else
        kputs(line->d.allele[1], str);
}
static void process_qual(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( bcf_float_is_missing(line->qual) ) kputc('.', str);
    else kputd(line->qual, str);
}
static void process_filter(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int i;
    if ( line->d.n_flt )
    {
        for (i=0; i<line->d.n_flt; i++)
        {
            if (i) kputc(';', str);
            kputs(convert->header->id[BCF_DT_ID][line->d.flt[i]].key, str);
        }
    }
    else kputc('.', str);
}
static inline int32_t bcf_array_ivalue(uint8_t *bcf_array, int type, int idx)
{
    if ( type==BCF_BT_INT8 )
    {
        int8_t val = le_to_i8(&bcf_array[idx * sizeof(val)]);
        if ( val==bcf_int8_missing ) return bcf_int32_missing;
        if ( val==bcf_int8_vector_end ) return bcf_int32_vector_end;
        return val;
    }
    if ( type==BCF_BT_INT16 )
    {
        int16_t val = le_to_i16(&bcf_array[idx * sizeof(val)]);
        if ( val==bcf_int16_missing ) return bcf_int32_missing;
        if ( val==bcf_int16_vector_end ) return bcf_int32_vector_end;
        return val;
    }
    return le_to_i32(&bcf_array[idx * sizeof(int32_t)]);
}
static inline void _copy_field(char *src, uint32_t len, int idx, kstring_t *str)
{
    int n = 0, ibeg = 0;
    while ( src[ibeg] && ibeg<len && n < idx )
    {
        if ( src[ibeg]==',' ) n++;
        ibeg++;
    }
    if ( ibeg==len ) { kputc('.', str); return; }

    int iend = ibeg;
    while ( src[iend] && src[iend]!=',' && iend<len ) iend++;

    if ( iend>ibeg )
        kputsn(src+ibeg, iend-ibeg, str);
    else
        kputc('.', str);
}
static void process_info(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int i;
    if ( !fmt->key )    // the whole INFO column
    {
        int first = 1;
        for (i=0; i<line->n_info; i++)
        {
            bcf_info_t *inf = &line->d.info[i];
            if ( !inf->vptr ) continue;
            if ( !first ) kputc(';', str);
            first = 0;
            if ( inf->key >= convert->header->n[BCF_DT_ID] ) continue;
            kputs(convert->header->id[BCF_DT_ID][inf->key].key, str);
            if ( inf->len <= 0 ) continue;
            kputc('=', str);
            if ( inf->len == 1 )
            {
                switch (inf->type)
                {
                    case BCF_BT_INT8:  if ( inf->v1.i==bcf_int8_missing ) kputc('.', str); else kputw(inf->v1.i, str); break;
                    case BCF_BT_INT16: if ( inf->v1.i==bcf_int16_missing ) kputc('.', str); else kputw(inf->v1.i, str); break;
                    case BCF_BT_INT32: if ( inf->v1.i==bcf_int32_missing ) kputc('.', str); else kputw(inf->v1.i, str); break;
                    case BCF_BT_FLOAT: if ( bcf_float_is_missing(inf->v1.f) ) kputc('.', str); else kputd(inf->v1.f, str); break;
                    case BCF_BT_CHAR:  kputc(inf->v1.i, str); break;
                    default: error("Unexpected type %d", inf->type); break;
                }
            }
            else bcf_fmt_array(str, inf->len, inf->type, inf->vptr);
        }
        if ( first ) kputc('.', str);
        return;
    }

    if ( fmt->id<0 )
    {
        kputc('.', str);
        return;
    }

    for (i=0; i<line->n_info; i++)
        if ( line->d.info[i].key == fmt->id ) break;

    // output "." if the tag is not present
    if ( i==line->n_info )
    {
        kputc('.', str);
        return;
    }

    bcf_info_t *info = &line->d.info[i];

    // if this is a flag, output 1
    if ( info->len <=0 )
    {
        kputc('1', str);
        return;
    }

    if ( info->len == 1 )
    {
        switch (info->type)
        {
            case BCF_BT_INT8:  if ( info->v1.i==bcf_int8_missing ) kputc('.', str); else kputw(info->v1.i, str); break;
            case BCF_BT_INT16: if ( info->v1.i==bcf_int16_missing ) kputc('.', str); else kputw(info->v1.i, str); break;
            case BCF_BT_INT32: if ( info->v1.i==bcf_int32_missing ) kputc('.', str); else kputw(info->v1.i, str); break;
            case BCF_BT_FLOAT: if ( bcf_float_is_missing(info->v1.f) ) kputc('.', str); else kputd(info->v1.f, str); break;
            case BCF_BT_CHAR:  kputc(info->v1.i, str); break;
            default: fprintf(stderr,"todo: type %d\n", info->type); exit(1); break;
        }
    }
    else if ( fmt->subscript >=0 )
    {
        if ( info->len <= fmt->subscript )
        {
            kputc('.', str);
            return;
        }
        #define BRANCH(type_t, convert, is_missing, is_vector_end, kprint) { \
            type_t val = convert(&info->vptr[fmt->subscript * sizeof(type_t)]); \
            if ( is_missing || is_vector_end ) kputc('.',str); \
            else kprint; \
        }
        switch (info->type)
        {
            case BCF_BT_INT8:  BRANCH(int8_t,  le_to_i8,  val==bcf_int8_missing,  val==bcf_int8_vector_end,  kputw(val, str)); break;
            case BCF_BT_INT16: BRANCH(int16_t, le_to_i16, val==bcf_int16_missing, val==bcf_int16_vector_end, kputw(val, str)); break;
            case BCF_BT_INT32: BRANCH(int32_t, le_to_i32, val==bcf_int32_missing, val==bcf_int32_vector_end, kputw(val, str)); break;
            case BCF_BT_FLOAT: BRANCH(float,   le_to_float, bcf_float_is_missing(val), bcf_float_is_vector_end(val), kputd(val, str)); break;
            case BCF_BT_CHAR:  _copy_field((char*)info->vptr, info->vptr_len, fmt->subscript, str); break;
            default: fprintf(stderr,"todo: type %d\n", info->type); exit(1); break;
        }
        #undef BRANCH
    }
    else
        bcf_fmt_array(str, info->len, info->type, info->vptr);
}
static void init_format(convert_t *convert, bcf1_t *line, fmt_t *fmt)
{
    fmt->id = bcf_hdr_id2int(convert->header, BCF_DT_ID, fmt->key);
    if ( !bcf_hdr_idinfo_exists(convert->header,BCF_HL_FMT,fmt->id) ) fmt->id = -1;
    fmt->fmt = NULL;
    if ( fmt->id >= 0 )
    {
        int i;
        for (i=0; i<(int)line->n_fmt; i++)
            if ( line->d.fmt[i].id==fmt->id ) { fmt->fmt = &line->d.fmt[i]; break; }
    }
    else if ( !convert->allow_undef_tags )
        error("Error: no such tag defined in the VCF header: FORMAT/%s\n", fmt->key);

    fmt->ready = 1;
}
static void process_complete_format(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( convert->nsamples )
    {
        int i,j;
        if ( line->n_fmt)
        {
            int gt_i = -1;
            bcf_fmt_t *fmt = line->d.fmt;
            int first = 1;
            for (i=0; i<(int)line->n_fmt; i++)
            {
                if ( !fmt[i].p || fmt[i].id<0 ) continue;
                if ( !first ) kputc(':', str);
                first = 0;
                kputs(convert->header->id[BCF_DT_ID][fmt[i].id].key, str);
                if ( strcmp(convert->header->id[BCF_DT_ID][fmt[i].id].key, "GT") == 0) gt_i = i;
            }
            if ( first ) kputc('.', str);
            for (j=0; j<convert->nsamples; j++)
            {
                kputc('\t', str);
                first = 1;
                for (i=0; i<(int)line->n_fmt; i++)
                {
                    bcf_fmt_t *f = &fmt[i];
                    if ( !f->p ) continue;
                    if ( !first ) kputc(':', str);
                    first = 0;
                    if (gt_i == i)
                        bcf_format_gt(f,convert->samples[j],str);
                    else
                        bcf_fmt_array(str, f->n, f->type, f->p + convert->samples[j] * f->size);
                }
                if ( first ) kputc('.', str);
            }
        }
        else
            for (j=0; j<=line->n_sample; j++)
                kputs("\t.", str);
    }
    else
        kputc('.',str);
}
static void process_format(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( !fmt->ready )
        init_format(convert, line, fmt);

    if ( fmt->fmt==NULL )
    {
        kputc('.', str);
        return;
    }
    else if ( fmt->subscript >=0 )
    {
        if ( fmt->fmt->n <= fmt->subscript )
        {
            kputc('.', str);
            return;
        }
        if ( fmt->fmt->type == BCF_BT_FLOAT )
        {
            uint8_t *ptr = fmt->fmt->p + isample*fmt->fmt->size;
            float val = le_to_float(&ptr[fmt->subscript * sizeof(float)]);
            if ( bcf_float_is_missing(val) || bcf_float_is_vector_end(val) )
                kputc('.', str);
            else
                kputd(val, str);
        }
        else if ( fmt->fmt->type != BCF_BT_CHAR )
        {
            int32_t ival = bcf_array_ivalue(fmt->fmt->p+isample*fmt->fmt->size,fmt->fmt->type,fmt->subscript);
            if ( ival==bcf_int32_missing || ival==bcf_int32_vector_end )
                kputc('.', str);
            else
                kputw(ival, str);
        }
        else if ( fmt->fmt->type == BCF_BT_CHAR )
            _copy_field((char*)(fmt->fmt->p + isample*fmt->fmt->size), fmt->fmt->size, fmt->subscript, str);
        else error("TODO: %s:%d .. fmt->type=%d\n", __FILE__,__LINE__, fmt->fmt->type);
    }
    else
        bcf_fmt_array(str, fmt->fmt->n, fmt->fmt->type, fmt->fmt->p + isample*fmt->fmt->size);
}
static void process_gt(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( !fmt->ready )
        init_format(convert, line, fmt);

    if ( fmt->fmt==NULL )
    {
        kputc('.', str);
        return;
    }
    bcf_format_gt(fmt->fmt, isample, str);
}
static void process_tgt(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( !fmt->ready )
        init_format(convert, line, fmt);

    if ( fmt->fmt==NULL )
    {
        kputc('.', str);
        return;
    }

    assert( fmt->fmt->type==BCF_BT_INT8 );

    int l;
    int8_t *x = (int8_t*)(fmt->fmt->p + isample*fmt->fmt->size); // FIXME: does not work with n_alt >= 64
    for (l = 0; l < fmt->fmt->n && x[l] != bcf_int8_vector_end; ++l)
    {
        if (l) kputc("/|"[x[l]&1], str);
        if (x[l]>>1)
        {
            int ial = (x[l]>>1) - 1;
            kputs(line->d.allele[ial], str);
        }
        else
            kputc('.', str);
    }
    if (l == 0) kputc('.', str);
}
static void destroy_tbcsq(void *usr)
{
    if ( !usr ) return;
    bcsq_t *csq = (bcsq_t*) usr;
    free(csq->hap1.s);
    free(csq->hap2.s);
    if ( csq->n )
        free(csq->str[0]);
    free(csq->str);
    free(csq);
}
static void process_tbcsq(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( !fmt->ready )
    {
        init_format(convert, line, fmt);

        bcsq_t *csq;
        if ( fmt->usr )
        {
            csq = (bcsq_t*) fmt->usr;
            if ( csq->n )
                free(csq->str[0]);
            csq->n = 0;
        }
        else
            csq = (bcsq_t*) calloc(1,sizeof(bcsq_t));
        fmt->usr = csq;

        int i=0, len = 0;
        char *tmp = NULL;
        if ( bcf_get_info_string(convert->header,line,fmt->key,&tmp,&len)<0 )
        {
            csq->n = 0;
            return;
        }
        do
        {
            csq->n++;
            hts_expand(char*, csq->n, csq->m, csq->str);
            csq->str[ csq->n-1 ] = tmp + i;
            while ( i<len && tmp[i]!=',' ) i++;
            if ( i<len && tmp[i]==',' ) tmp[i++] = 0;
        }
        while ( i<len );
    }

    bcsq_t *csq = (bcsq_t*)fmt->usr;

    if ( fmt->fmt==NULL || !csq->n ) return;

    csq->hap1.l = 0;
    csq->hap2.l = 0;

    int mask = fmt->subscript==0 ? 3 : 1;   // merge both haplotypes if subscript==0

    #define BRANCH(type_t, convert, nbits) { \
        uint8_t *x = fmt->fmt->p + isample*fmt->fmt->size; \
        int i,j; \
        if ( fmt->subscript<=0 || fmt->subscript==1 ) \
        { \
            for (j=0; j < fmt->fmt->n; j++) \
            { \
                type_t val = convert(&x[j * sizeof(type_t)]); \
                if ( !val ) continue; \
                for (i=0; i<nbits; i+=2) \
                    if ( val & (mask<<i) ) { kputs(csq->str[(j*30+i)/2], &csq->hap1); kputc_(',', &csq->hap1); } \
            } \
        } \
        if ( fmt->subscript<0 || fmt->subscript==2 ) \
        { \
            for (j=0; j < fmt->fmt->n; j++) \
            { \
                type_t val = convert(&x[j * sizeof(type_t)]); \
                if ( !val ) continue; \
                for (i=1; i<nbits; i+=2) \
                    if ( val & (1<<i) ) { kputs(csq->str[(j*30+i)/2], &csq->hap2); kputc_(',', &csq->hap2); } \
            } \
        } \
    }
    switch (fmt->fmt->type)
    {
        case BCF_BT_INT8:  BRANCH(uint8_t,  le_to_u8,   8); break;
        case BCF_BT_INT16: BRANCH(uint16_t, le_to_u16, 16); break;
        case BCF_BT_INT32: BRANCH(uint32_t, le_to_u32, 30); break;  // 2 bits unused to account for the reserved BCF values
        default: error("Unexpected type: %d\n", fmt->fmt->type); exit(1); break;
    }
    #undef BRANCH

    if ( !csq->hap1.l && !csq->hap2.l ) return;

    if ( csq->hap1.l ) csq->hap1.s[--csq->hap1.l] = 0;
    if ( csq->hap2.l ) csq->hap2.s[--csq->hap2.l] = 0;

    if ( fmt->subscript<0 )
    {
        kputs(csq->hap1.l?csq->hap1.s:".", str);
        kputc_('\t', str);
        kputs(csq->hap2.l?csq->hap2.s:".", str);
    }
    else if ( fmt->subscript<2 )
        kputs(csq->hap1.l?csq->hap1.s:".", str);
    else
        kputs(csq->hap2.l?csq->hap2.s:".", str);
}
static void init_format_iupac(convert_t *convert, bcf1_t *line, fmt_t *fmt)
{
    init_format(convert, line, fmt);
    if ( fmt->fmt==NULL ) return;

    // Init mapping between alleles and IUPAC table
    hts_expand(uint8_t, line->n_allele, convert->ndat, convert->dat);
    int8_t *dat = (int8_t*)convert->dat;
    int i;
    for (i=0; i<line->n_allele; i++)
    {
        if ( line->d.allele[i][1] ) dat[i] = -1;
        else
        {
            switch (line->d.allele[i][0])
            {
                case 'A': dat[i] = 0; break;
                case 'C': dat[i] = 1; break;
                case 'G': dat[i] = 2; break;
                case 'T': dat[i] = 3; break;
                case 'a': dat[i] = 0; break;
                case 'c': dat[i] = 1; break;
                case 'g': dat[i] = 2; break;
                case 't': dat[i] = 3; break;
                default: dat[i] = -1;
            }
        }
    }
}
static void process_iupac_gt(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( !fmt->ready )
        init_format_iupac(convert, line, fmt);

    if ( fmt->fmt==NULL )
    {
        kputc('.', str);
        return;
    }

    assert( fmt->fmt->type==BCF_BT_INT8 );

    static const char iupac[4][4] = { {'A','M','R','W'},{'M','C','S','Y'},{'R','S','G','K'},{'W','Y','K','T'} };
    int8_t *dat = (int8_t*)convert->dat;

    int8_t *x = (int8_t*)(fmt->fmt->p + isample*fmt->fmt->size); // FIXME: does not work with n_alt >= 64
    int l = 0;
    while ( l<fmt->fmt->n && x[l]!=bcf_int8_vector_end && x[l]!=bcf_int8_missing ) l++;

    if ( l==2 )
    {
        // diploid
        int ia = (x[0]>>1) - 1, ib = (x[1]>>1) - 1;
        if ( ia>=0 && ia<line->n_allele && ib>=0 && ib<line->n_allele && dat[ia]>=0 && dat[ib]>=0 )
        {
            kputc(iupac[dat[ia]][dat[ib]], str);
            return;
        }
    }
    for (l = 0; l < fmt->fmt->n && x[l] != bcf_int8_vector_end; ++l)
    {
        if (l) kputc("/|"[x[l]&1], str);
        if (x[l]>>1)
        {
            int ial = (x[l]>>1) - 1;
            kputs(line->d.allele[ial], str);
        }
        else
            kputc('.', str);
    }
    if (l == 0) kputc('.', str);
}
static void process_sample(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    kputs(convert->header->samples[isample], str);
}
static void process_sep(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str) { if (fmt->key) kputs(fmt->key, str); }
static void process_is_ts(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int is_ts = 0;
    if ( bcf_get_variant_types(line) & (VCF_SNP|VCF_MNP) )
        is_ts = abs(bcf_acgt2int(*line->d.allele[0])-bcf_acgt2int(*line->d.allele[1])) == 2 ? 1 : 0;
    kputc(is_ts ? '1' : '0', str);
}
static void process_type(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int line_type = bcf_get_variant_types(line);
    int i = 0;
    if ( line_type == VCF_REF ) { kputs("REF", str); i++; }
    if ( line_type & VCF_SNP ) { if (i) kputc(',',str); kputs("SNP", str); i++; }
    if ( line_type & VCF_MNP ) { if (i) kputc(',',str); kputs("MNP", str); i++; }
    if ( line_type & VCF_INDEL ) { if (i) kputc(',',str); kputs("INDEL", str); i++; }
    if ( line_type & VCF_OTHER ) { if (i) kputc(',',str); kputs("OTHER", str); i++; }
    if ( line_type & VCF_BND ) { if (i) kputc(',',str); kputs("BND", str); i++; }
    if ( line_type & VCF_OVERLAP ) { if (i) kputc(',',str); kputs("OVERLAP", str); i++; }
}
static void process_line(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    vcf_format1(convert->header, line, str);
    if ( str->s[str->l-1]=='\n' ) str->l--;
}
static void process_chrom_pos_id(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    if ( line->d.id[0]!='.' || line->d.id[1] )
    {
        // ID is present
        kputs(line->d.id, str);
    }
    else
    {
        // use CHROM:POS instead of ID
        kputs(convert->header->id[BCF_DT_CTG][line->rid].key, str);
        kputc(':', str);
        kputw(line->pos+1, str);
    }
}
static void process_gt_to_prob3(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int m,n,i;

    m = convert->ndat / sizeof(int32_t);
    n = bcf_get_genotypes(convert->header,line,&convert->dat,&m);
    convert->ndat = m * sizeof(int32_t);

    if ( n<=0 )
    {
        // Throw an error or silently proceed?
        //
        // for (i=0; i<convert->nsamples; i++) kputs(" 0.33 0.33 0.33", str);
        // return;

        error("Error parsing GT tag at %s:%"PRId64"\n", bcf_seqname(convert->header,line),(int64_t) line->pos+1);
    }

    n /= convert->nsamples;
    for (i=0; i<convert->nsamples; i++)
    {
        int32_t *ptr = (int32_t*)convert->dat + i*n;
        int j;
        for (j=0; j<n; j++)
            if ( ptr[j]==bcf_int32_vector_end ) break;

        if ( j==2 )
        {
            // diploid
            if ( bcf_gt_is_missing(ptr[0]) )
                kputs(" 0.33 0.33 0.33", str);
            else if ( bcf_gt_allele(ptr[0])!=bcf_gt_allele(ptr[1]) )
                kputs(" 0 1 0", str);       // HET
            else if ( bcf_gt_allele(ptr[0])==1 )
                kputs(" 0 0 1", str);       // ALT HOM, first ALT allele
            else
                kputs(" 1 0 0", str);       // REF HOM or something else than first ALT
        }
        else if ( j==1 )
        {
            // haploid
            if ( bcf_gt_is_missing(ptr[0]) )
                kputs(" 0.5 0.0 0.5", str);
            else if ( bcf_gt_allele(ptr[0])==1 )
                kputs(" 0 0 1", str);       // first ALT allele
            else
                kputs(" 1 0 0", str);       // REF or something else than first ALT
        }
        else error("FIXME: not ready for ploidy %d\n", j);
    }
}
static void process_pl_to_prob3(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int m,n,i;

    m = convert->ndat / sizeof(int32_t);
    n = bcf_get_format_int32(convert->header,line,"PL",&convert->dat,&m);
    convert->ndat = m * sizeof(int32_t);

    if ( n<=0 )
    {
        // Throw an error or silently proceed?
        //
        // for (i=0; i<convert->nsamples; i++) kputs(" 0.33 0.33 0.33", str);
        // return;

        error("Error parsing PL tag at %s:%"PRId64"\n", bcf_seqname(convert->header,line),(int64_t) line->pos+1);
    }

    n /= convert->nsamples;
    for (i=0; i<convert->nsamples; i++)
    {
        int32_t *ptr = (int32_t*)convert->dat + i*n;
        int j;
        float sum = 0;
        for (j=0; j<n; j++)
        {
            if ( ptr[j]==bcf_int32_vector_end ) break;
            sum += pow(10,-0.1*ptr[j]);
        }
        if ( j==line->n_allele )
        {
            // haploid
            kputc(' ',str);
            ksprintf(str,"%f",pow(10,-0.1*ptr[0])/sum);
            kputs(" 0 ", str);
            ksprintf(str,"%f",pow(10,-0.1*ptr[1])/sum);
        }
        else
        {
            // diploid
            kputc(' ',str);
            ksprintf(str,"%f",pow(10,-0.1*ptr[0])/sum);
            kputc(' ',str);
            ksprintf(str,"%f",pow(10,-0.1*ptr[1])/sum);
            kputc(' ',str);
            ksprintf(str,"%f",pow(10,-0.1*ptr[2])/sum);
        }
    }
}
static void process_gp_to_prob3(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int m,n,i;

    m = convert->ndat / sizeof(float);
    n = bcf_get_format_float(convert->header,line,"GP",&convert->dat,&m);
    convert->ndat = m * sizeof(float);

    if ( n<=0 )
    {
        // Throw an error or silently proceed?
        //
        // for (i=0; i<convert->nsamples; i++) kputs(" 0.33 0.33 0.33", str);
        // return;

        error("Error parsing GP tag at %s:%"PRId64"\n", bcf_seqname(convert->header,line),(int64_t) line->pos+1);
    }

    n /= convert->nsamples;
    for (i=0; i<convert->nsamples; i++)
    {
        float *ptr = (float*)convert->dat + i*n;
        int j;
        for (j=0; j<n; j++)
        {
            if ( bcf_float_is_vector_end(ptr[j]) ) break;
            if ( bcf_float_is_missing(ptr[j]) ) { ptr[j]=0; continue; }
            if ( ptr[j]<0 || ptr[j]>1 ) error("[%s:%"PRId64":%f] GP value outside range [0,1]; bcftools convert expects the VCF4.3+ spec for the GP field encoding genotype posterior probabilities", bcf_seqname(convert->header,line),(int64_t) line->pos+1,ptr[j]);
        }
        if ( j==line->n_allele )
            ksprintf(str," %f %f %f",ptr[0],0.,ptr[1]); // haploid
        else
            ksprintf(str," %f %f %f",ptr[0],ptr[1],ptr[2]); // diploid
    }
}

static void process_gt_to_hap(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    // https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-known_haps_g

    // File containing known haplotypes for the study cohort. The format
    // is the same as the output format from IMPUTE2's -phase option:
    // five header columns (as in the -g file) followed by two columns
    // (haplotypes) per individual. Allowed values in the haplotype
    // columns are 0, 1, and ?.

    // If your study dataset is fully phased, you can replace the -g file
    // with a -known_haps_g file. This will cause IMPUTE2 to perform
    // haploid imputation, although it will still report diploid imputation
    // probabilities in the main output file. If any genotypes are missing,
    // they can be marked as '? ?' (two question marks separated by one
    // space) in the input file. (The program does not allow just one
    // allele from a diploid genotype to be missing.) If the reference
    // panels are also phased, IMPUTE2 will perform a single, fast
    // imputation step rather than its standard MCMC module this is how
    // the program imputes into pre-phased GWAS haplotypes.

    // The -known_haps_g file can also be used to specify study
    // genotypes that are "partially" phased, in the sense that some
    // genotypes are phased relative to a fixed reference point while
    // others are not. We anticipate that this will be most useful when
    // trying to phase resequencing data onto a scaffold of known
    // haplotypes. To mark a known genotype as unphased, place an
    // asterisk immediately after each allele, with no space between
    // the allele (0/1) and the asterisk (*); e.g., "0* 1*" for a
    // heterozygous genotype of unknown phase.

    int i, gt_id = bcf_hdr_id2int(convert->header, BCF_DT_ID, "GT");
    if ( !bcf_hdr_idinfo_exists(convert->header,BCF_HL_FMT,gt_id) )
        error("FORMAT/GT tag not present at %s:%"PRId64"\n", bcf_seqname(convert->header, line),(int64_t) line->pos+1);
    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);
    bcf_fmt_t *fmt_gt = NULL;
    for (i=0; i<line->n_fmt; i++)
        if ( line->d.fmt[i].id==gt_id ) { fmt_gt = &line->d.fmt[i]; break; }
    if ( !fmt_gt )
        error("FORMAT/GT tag not present at %s:%"PRId64"\n", bcf_seqname(convert->header, line),(int64_t) line->pos+1);

    // Alloc all memory in advance to avoid kput routines. The biggest allowed allele index is 99
    if ( line->n_allele > 100 )
        error("Too many alleles (%d) at %s:%"PRId64"\n", line->n_allele, bcf_seqname(convert->header, line),(int64_t) line->pos+1);
    if ( ks_resize(str, str->l+convert->nsamples*8) != 0 )
        error("Could not alloc %" PRIu64 " bytes\n", (uint64_t)(str->l + convert->nsamples*8));

    if ( fmt_gt->type!=BCF_BT_INT8 )    // todo: use BRANCH_INT if the VCF is valid
        error("Uh, too many alleles (%d) or redundant BCF representation at %s:%"PRId64"\n", line->n_allele, bcf_seqname(convert->header, line),(int64_t) line->pos+1);
    if ( fmt_gt->n!=1 && fmt_gt->n!=2 )
        error("Uh, ploidy of %d not supported, see %s:%"PRId64"\n", fmt_gt->n, bcf_seqname(convert->header, line),(int64_t) line->pos+1);

    int8_t *ptr = ((int8_t*) fmt_gt->p) - fmt_gt->n;
    for (i=0; i<convert->nsamples; i++)
    {
        ptr += fmt_gt->n;
        if ( fmt_gt->n==1 ) // haploid genotypes
        {
            if ( ptr[0]==2 ) /* 0 */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = ' '; str->s[str->l++] = '-'; str->s[str->l++] = ' ';
            }
            else if ( ptr[0]==bcf_int8_missing )   /* . */
            {
                str->s[str->l++] = '?'; str->s[str->l++] = ' '; str->s[str->l++] = '?'; str->s[str->l++] = ' ';
            }
            else if ( ptr[0]==4 ) /* 1 */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = ' '; str->s[str->l++] = '-'; str->s[str->l++] = ' ';
            }
            else
            {
                kputw(bcf_gt_allele(ptr[0]),str); str->s[str->l++] = ' '; str->s[str->l++] = '-'; str->s[str->l++] = ' ';
            }
        }
        else if ( ptr[0]==2 )
        {
            if ( ptr[1]==3 ) /* 0|0 */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = ' '; str->s[str->l++] = '0'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==5 ) /* 0|1 */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = ' '; str->s[str->l++] = '1'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==bcf_int8_vector_end ) /* 0 */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = ' '; str->s[str->l++] = '-'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==2 ) /* 0/0 */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = '*'; str->s[str->l++] = ' '; str->s[str->l++] = '0'; str->s[str->l++] = '*'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==4 ) /* 0/1 */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = '*'; str->s[str->l++] = ' '; str->s[str->l++] = '1'; str->s[str->l++] = '*'; str->s[str->l++] = ' ';
            }
            else if ( bcf_gt_is_missing(ptr[1]) ) /* 0/. */
            {
                str->s[str->l++] = '?'; str->s[str->l++] = ' '; str->s[str->l++] = '?'; str->s[str->l++] = ' ';
            }
            else if ( bcf_gt_is_phased(ptr[1]) ) /* 0|x */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = ' ';
                kputw(bcf_gt_allele(ptr[1]),str);
                str->s[str->l++] = ' ';
            }
            else /* 0/x */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = '*'; str->s[str->l++] = ' ';
                kputw(bcf_gt_allele(ptr[1]),str);
                str->s[str->l++] = '*'; str->s[str->l++] = ' ';
            }
        }
        else if ( ptr[0]==4 )
        {
            if ( ptr[1]==3 ) /* 1|0 */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = ' '; str->s[str->l++] = '0'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==5 ) /* 1|1 */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = ' '; str->s[str->l++] = '1'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==bcf_int8_vector_end ) /* 1 */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = ' '; str->s[str->l++] = '-'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==2 ) /* 1/0 */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = '*'; str->s[str->l++] = ' '; str->s[str->l++] = '0'; str->s[str->l++] = '*'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==4 ) /* 1/1 */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = '*'; str->s[str->l++] = ' '; str->s[str->l++] = '1'; str->s[str->l++] = '*'; str->s[str->l++] = ' ';
            }
            else if ( bcf_gt_is_missing(ptr[1]) ) /* 1/. */
            {
                str->s[str->l++] = '?'; str->s[str->l++] = ' '; str->s[str->l++] = '?'; str->s[str->l++] = ' ';
            }
            else if ( bcf_gt_is_phased(ptr[1]) )    /* 1|x */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = ' ';
                kputw(bcf_gt_allele(ptr[1]),str);
                str->s[str->l++] = ' ';
            }
            else /* 1/x */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = '*'; str->s[str->l++] = ' ';
                kputw(bcf_gt_allele(ptr[1]),str);
                str->s[str->l++] = '*'; str->s[str->l++] = ' ';
            }
        }
        else if ( bcf_gt_is_missing(ptr[0]) )
        {
            if ( ptr[1]==bcf_int8_vector_end )
            {
                str->s[str->l++] = '?'; str->s[str->l++] = ' '; str->s[str->l++] = '-'; str->s[str->l++] = ' ';
            }
            else
            {
                str->s[str->l++] = '?'; str->s[str->l++] = ' '; str->s[str->l++] = '?'; str->s[str->l++] = ' ';
            }
        }
        else if ( ptr[1]==bcf_int8_vector_end )
        {
            /* use REF for something else than first ALT */
            str->s[str->l++] = '0'; str->s[str->l++] = ' '; str->s[str->l++] = '-'; str->s[str->l++] = ' ';
        }
        else
        {
            kputw(bcf_gt_allele(ptr[0]),str);
            if ( bcf_gt_is_phased(ptr[1]) ) str->s[str->l++] = '*';
            str->s[str->l++] = ' ';
            kputw(bcf_gt_allele(ptr[1]),str);
            if ( bcf_gt_is_phased(ptr[1]) ) str->s[str->l++] = '*';
            str->s[str->l++] = ' ';
        }
    }
    str->s[--str->l] = 0;     // delete the last space
}
static void process_gt_to_hap2(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    // same as process_gt_to_hap but converts haploid genotypes into diploid

    int i, gt_id = bcf_hdr_id2int(convert->header, BCF_DT_ID, "GT");
    if ( !bcf_hdr_idinfo_exists(convert->header,BCF_HL_FMT,gt_id) )
        error("FORMAT/GT tag not present at %s:%"PRId64"\n", bcf_seqname(convert->header, line),(int64_t) line->pos+1);
    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);
    bcf_fmt_t *fmt_gt = NULL;
    for (i=0; i<line->n_fmt; i++)
        if ( line->d.fmt[i].id==gt_id ) { fmt_gt = &line->d.fmt[i]; break; }
    if ( !fmt_gt )
        error("FORMAT/GT tag not present at %s:%"PRId64"\n", bcf_seqname(convert->header, line),(int64_t)  line->pos+1);

    // Alloc all memory in advance to avoid kput routines. The biggest allowed allele index is 99
    if ( line->n_allele > 100 )
        error("Too many alleles (%d) at %s:%"PRId64"\n", line->n_allele, bcf_seqname(convert->header, line),(int64_t) line->pos+1);
    if ( ks_resize(str, str->l+convert->nsamples*8) != 0 )
        error("Could not alloc %" PRIu64 " bytes\n", (uint64_t)(str->l + convert->nsamples*8));

    if ( fmt_gt->type!=BCF_BT_INT8 )    // todo: use BRANCH_INT if the VCF is valid
        error("Uh, too many alleles (%d) or redundant BCF representation at %s:%"PRId64"\n", line->n_allele, bcf_seqname(convert->header, line),(int64_t) line->pos+1);

    int8_t *ptr = ((int8_t*) fmt_gt->p) - fmt_gt->n;
    for (i=0; i<convert->nsamples; i++)
    {
        ptr += fmt_gt->n;
        if ( ptr[0]==2 )
        {
            if ( ptr[1]==3 ) /* 0|0 */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = ' '; str->s[str->l++] = '0'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==5 ) /* 0|1 */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = ' '; str->s[str->l++] = '1'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==bcf_int8_vector_end ) /* 0 -> 0|0 */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = ' '; str->s[str->l++] = '0'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==2 ) /* 0/0 */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = '*'; str->s[str->l++] = ' '; str->s[str->l++] = '0'; str->s[str->l++] = '*'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==4 ) /* 0/1 */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = '*'; str->s[str->l++] = ' '; str->s[str->l++] = '1'; str->s[str->l++] = '*'; str->s[str->l++] = ' ';
            }
            else if ( bcf_gt_is_missing(ptr[1]) ) /* 0/. */
            {
                str->s[str->l++] = '?'; str->s[str->l++] = ' '; str->s[str->l++] = '?'; str->s[str->l++] = ' ';
            }
            else if ( bcf_gt_is_phased(ptr[1]) ) /* 0|x */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = ' ';
                kputw(bcf_gt_allele(ptr[1]),str);
                str->s[str->l++] = ' ';
            }
            else /* 0/x */
            {
                str->s[str->l++] = '0'; str->s[str->l++] = '*'; str->s[str->l++] = ' ';
                kputw(bcf_gt_allele(ptr[1]),str);
                str->s[str->l++] = '*'; str->s[str->l++] = ' ';
            }
        }
        else if ( ptr[0]==4 )
        {
            if ( ptr[1]==3 ) /* 1|0 */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = ' '; str->s[str->l++] = '0'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==5 ) /* 1|1 */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = ' '; str->s[str->l++] = '1'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==bcf_int8_vector_end ) /* 1 -> 1|1 */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = ' '; str->s[str->l++] = '1'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==2 ) /* 1/0 */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = '*'; str->s[str->l++] = ' '; str->s[str->l++] = '0'; str->s[str->l++] = '*'; str->s[str->l++] = ' ';
            }
            else if ( ptr[1]==4 ) /* 1/1 */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = '*'; str->s[str->l++] = ' '; str->s[str->l++] = '1'; str->s[str->l++] = '*'; str->s[str->l++] = ' ';
            }
            else if ( bcf_gt_is_missing(ptr[1]) ) /* 1/. */
            {
                str->s[str->l++] = '?'; str->s[str->l++] = ' '; str->s[str->l++] = '?'; str->s[str->l++] = ' ';
            }
            else if ( bcf_gt_is_phased(ptr[1]) )    /* 1|x */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = ' ';
                kputw(bcf_gt_allele(ptr[1]),str);
                str->s[str->l++] = ' ';
            }
            else /* 1/x */
            {
                str->s[str->l++] = '1'; str->s[str->l++] = '*'; str->s[str->l++] = ' ';
                kputw(bcf_gt_allele(ptr[1]),str);
                str->s[str->l++] = '*'; str->s[str->l++] = ' ';
            }
        }
        else if ( bcf_gt_is_missing(ptr[0]) )
        {
            str->s[str->l++] = '?'; str->s[str->l++] = ' '; str->s[str->l++] = '?'; str->s[str->l++] = ' ';
        }
        else if ( ptr[1]==bcf_int8_vector_end )
        {
            /* use REF for something else than first ALT */
            str->s[str->l++] = '0'; str->s[str->l++] = ' '; str->s[str->l++] = '0'; str->s[str->l++] = ' ';
        }
        else
        {
            kputw(bcf_gt_allele(ptr[0]),str);
            if ( bcf_gt_is_phased(ptr[1]) ) str->s[str->l++] = '*';
            str->s[str->l++] = ' ';
            kputw(bcf_gt_allele(ptr[1]),str);
            if ( bcf_gt_is_phased(ptr[1]) ) str->s[str->l++] = '*';
            str->s[str->l++] = ' ';
        }
    }
    str->s[--str->l] = 0;     // delete the last space
}

static void process_rsid_hex(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    char *ptr = line->d.id;
    ptr += 2; // remove 'rs'
    ksprintf(str, "%08" PRIx32 "", (uint32_t)strtoul(ptr, NULL, 10));
}

static void process_variantkey_hex(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    const char *alt = NULL;
    size_t sizealt = 0;
    if ( line->n_allele>1 )
    {
        alt = line->d.allele[1];
        sizealt = strlen(line->d.allele[1]);
    }
    uint64_t vk = variantkey(
        convert->header->id[BCF_DT_CTG][line->rid].key,
        strlen(convert->header->id[BCF_DT_CTG][line->rid].key),
        line->pos,
        line->d.allele[0],
        strlen(line->d.allele[0]),
        alt,
        sizealt);
    ksprintf(str, "%016" PRIx64 "", vk);
}

static void process_npass(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int i, nsmpl = 0;
    filter_t *flt = (filter_t*) fmt->usr;
    const uint8_t *smpl;
    filter_test(flt,line,&smpl);
    for (i=0; i<convert->nsamples; i++)
        if ( smpl[i] ) nsmpl++;
    kputd(nsmpl, str);
}
static void destroy_npass(void *usr)
{
    filter_destroy((filter_t*)usr);
}

static void process_pbinom(convert_t *convert, bcf1_t *line, fmt_t *fmt, int isample, kstring_t *str)
{
    int i;
    if ( !fmt->ready )
    {
        fmt->fmt = NULL;    // AD
        fmt->usr = NULL;    // GT

        for (i=0; i<(int)line->n_fmt; i++)
            if ( line->d.fmt[i].id==fmt->id ) { fmt->fmt = &line->d.fmt[i]; break; }

        // Check that the first field is GT
        int gt_id = bcf_hdr_id2int(convert->header, BCF_DT_ID, "GT");
        if ( !bcf_hdr_idinfo_exists(convert->header, BCF_HL_FMT, fmt->id)  ) error("Error: FORMAT/GT is not defined in the header\n");
        for (i=0; i<(int)line->n_fmt; i++)
            if ( line->d.fmt[i].id==gt_id ) { fmt->usr = &line->d.fmt[i]; break; }  // it should always be first according to VCF spec, but...

        if ( fmt->usr && line->d.fmt[i].type!=BCF_BT_INT8 )   // skip sites with many alleles
            fmt->usr = NULL;

        fmt->ready = 1;
    }
    bcf_fmt_t *gt_fmt = (bcf_fmt_t*) fmt->usr;
    if ( !fmt->fmt || !gt_fmt || gt_fmt->n!=2 ) goto invalid;

    int n[2] = {0,0};
    int8_t *gt = (int8_t*)(gt_fmt->p + isample*gt_fmt->size);
    for (i=0; i<2; i++)
    {
        if ( bcf_gt_is_missing(gt[i]) || gt[i] == bcf_int8_vector_end ) goto invalid;
        int al = bcf_gt_allele(gt[i]);
        if ( al > line->n_allele || al >= fmt->fmt->n ) goto invalid;

        #define BRANCH(type_t, convert, missing, vector_end) { \
            type_t val = convert(&fmt->fmt->p[(al + isample*fmt->fmt->n)*sizeof(type_t)]); \
            if ( val==missing || val==vector_end ) goto invalid; \
            else n[i] = val; \
        }
        switch (fmt->fmt->type)
        {
            case BCF_BT_INT8:  BRANCH(int8_t,  le_to_i8,  bcf_int8_missing,  bcf_int8_vector_end); break;
            case BCF_BT_INT16: BRANCH(int16_t, le_to_i16, bcf_int16_missing, bcf_int16_vector_end); break;
            case BCF_BT_INT32: BRANCH(int32_t, le_to_i32, bcf_int32_missing, bcf_int32_vector_end); break;
            default: goto invalid; break;
        }
        #undef BRANCH
    }

    if ( n[0]==n[1] ) kputc(n[0]==0 ? '.':'0', str);
    else
    {
        double pval = calc_binom_two_sided(n[0],n[1],0.5);

        // convrt to phred
        if ( pval>=1 ) pval = 0;
        else pval = -4.34294481903*log(pval);
        kputd(pval, str);
    }
    return;

invalid:
    kputc('.', str);
}

static void _used_tags_add(convert_t *convert, int type, char *key)
{
    kstring_t str = {0,0,0};
    ksprintf(&str,"%s/%s",type==T_INFO?"INFO":"FORMAT",key);
    khash_str2int_inc(convert->used_tags_hash,str.s);
    convert->nused_tags++;
    convert->used_tags_list = (char**)realloc(convert->used_tags_list,sizeof(*convert->used_tags_list)*convert->nused_tags);
    convert->used_tags_list[convert->nused_tags-1] = str.s;
}


#define _SET_NON_FORMAT_TAGS(function,key,...) \
    if ( !strcmp("CHROM",key) ) { function(__VA_ARGS__, T_CHROM); } \
    else if ( !strcmp("POS",key) ) { function(__VA_ARGS__, T_POS); } \
    else if ( !strcmp("POS0",key) ) { function(__VA_ARGS__, T_POS0); } \
    else if ( !strcmp("END",key) ) { function(__VA_ARGS__, T_END); } \
    else if ( !strcmp("END0",key) ) { function(__VA_ARGS__, T_END0); } \
    else if ( !strcmp("ID",key) ) { function(__VA_ARGS__, T_ID); } \
    else if ( !strcmp("REF",key) ) { function(__VA_ARGS__, T_REF); } \
    else if ( !strcmp("FIRST_ALT",key) ) { function(__VA_ARGS__, T_FIRST_ALT); } \
    else if ( !strcmp("QUAL",key) ) { function(__VA_ARGS__, T_QUAL); } \
    else if ( !strcmp("TYPE",key) ) { function(__VA_ARGS__, T_TYPE); } \
    else if ( !strcmp("FILTER",key) ) { function(__VA_ARGS__, T_FILTER); } \
    else if ( !strcmp("IS_TS",key) ) { function(__VA_ARGS__, T_IS_TS); } \
    else if ( !strcmp("MASK",key) ) { function(__VA_ARGS__, T_MASK); } \
    else if ( !strcmp("LINE",key) ) { function(__VA_ARGS__, T_LINE); }

static void set_type(fmt_t *fmt, int type) { fmt->type = type; }
static fmt_t *register_tag(convert_t *convert, char *key, int is_gtf, int type)
{
    convert->nfmt++;
    if ( convert->nfmt > convert->mfmt )
    {
        convert->mfmt += 10;
        convert->fmt   = (fmt_t*) realloc(convert->fmt, convert->mfmt*sizeof(fmt_t));
    }
    fmt_t *fmt = &convert->fmt[ convert->nfmt-1 ];
    fmt->type  = type;
    fmt->key   = key ? strdup(key) : NULL;
    fmt->is_gt_field = is_gtf;
    fmt->subscript = -1;
    fmt->usr     = NULL;
    fmt->destroy = NULL;

    // Allow non-format tags, such as CHROM, INFO, etc., to appear amongst the format tags.
    if ( key )
    {
        int id = bcf_hdr_id2int(convert->header, BCF_DT_ID, key);
        if ( fmt->type==T_FORMAT && !bcf_hdr_idinfo_exists(convert->header,BCF_HL_FMT,id) )
        {
            _SET_NON_FORMAT_TAGS(set_type,key,fmt)
           else if ( !strcmp("ALT",key) ) { fmt->type = T_ALT; }
           else if ( !strcmp("_CHROM_POS_ID",key) ) { fmt->type = T_CHROM_POS_ID; }
            else if ( !strcmp("RSX",key) ) { fmt->type = T_RSX; }
            else if ( !strcmp("VKX",key) ) { fmt->type = T_VKX; }
            else if ( id>=0 && bcf_hdr_idinfo_exists(convert->header,BCF_HL_INFO,id) )
            {
                fmt->type = T_INFO;
                _used_tags_add(convert,T_INFO,key);
            }
        }
        else if ( fmt->type==T_PBINOM )
        {
            fmt->id = bcf_hdr_id2int(convert->header, BCF_DT_ID, fmt->key);
            if ( !bcf_hdr_idinfo_exists(convert->header,BCF_HL_FMT, fmt->id)  ) error("No such FORMAT tag defined in the header: %s\n", fmt->key);
            _used_tags_add(convert,T_FORMAT,key);
        }
        else if ( fmt->type==T_NPASS )
        {
            filter_t *flt = filter_init(convert->header,key);
            convert->max_unpack |= filter_max_unpack(flt);
            fmt->usr = (void*) flt;
        }
    }

    switch (fmt->type)
    {
        case T_FIRST_ALT: fmt->handler = &process_first_alt; break;
        case T_CHROM_POS_ID: fmt->handler = &process_chrom_pos_id; break;
        case T_GT_TO_PROB3: fmt->handler = &process_gt_to_prob3; break;
        case T_PL_TO_PROB3: fmt->handler = &process_pl_to_prob3; break;
        case T_GP_TO_PROB3: fmt->handler = &process_gp_to_prob3; break;
        case T_CHROM: fmt->handler = &process_chrom; break;
        case T_POS: fmt->handler = &process_pos; break;
        case T_POS0: fmt->handler = &process_pos0; break;
        case T_END: fmt->handler = &process_end; convert->max_unpack |= BCF_UN_INFO; break;
        case T_END0: fmt->handler = &process_end0; convert->max_unpack |= BCF_UN_INFO; break;
        case T_ID: fmt->handler = &process_id; break;
        case T_REF: fmt->handler = &process_ref; break;
        case T_ALT: fmt->handler = &process_alt; break;
        case T_QUAL: fmt->handler = &process_qual; break;
        case T_FILTER: fmt->handler = &process_filter; convert->max_unpack |= BCF_UN_FLT; break;
        case T_INFO: fmt->handler = &process_info; convert->max_unpack |= BCF_UN_INFO; break;
        case T_FORMAT: fmt->handler = fmt->key ? &process_format : &process_complete_format; convert->max_unpack |= BCF_UN_FMT; break;
        case T_SAMPLE: fmt->handler = &process_sample; break;
        case T_SEP: fmt->handler = &process_sep; break;
        case T_IS_TS: fmt->handler = &process_is_ts; break;
        case T_TYPE: fmt->handler = &process_type; break;
        case T_MASK: fmt->handler = NULL; break;
        case T_GT: fmt->handler = &process_gt; convert->max_unpack |= BCF_UN_FMT; break;
        case T_TGT: fmt->handler = &process_tgt; convert->max_unpack |= BCF_UN_FMT; break;
        case T_IUPAC_GT: fmt->handler = &process_iupac_gt; convert->max_unpack |= BCF_UN_FMT; break;
        case T_GT_TO_HAP: fmt->handler = &process_gt_to_hap; convert->max_unpack |= BCF_UN_FMT; break;
        case T_GT_TO_HAP2: fmt->handler = &process_gt_to_hap2; convert->max_unpack |= BCF_UN_FMT; break;
        case T_TBCSQ: fmt->handler = &process_tbcsq; fmt->destroy = &destroy_tbcsq; convert->max_unpack |= BCF_UN_FMT; break;
        case T_LINE: fmt->handler = &process_line; convert->max_unpack |= BCF_UN_FMT; break;
        case T_RSX: fmt->handler = &process_rsid_hex; break;
        case T_VKX: fmt->handler = &process_variantkey_hex; break;
        case T_PBINOM: fmt->handler = &process_pbinom; convert->max_unpack |= BCF_UN_FMT; break;
        case T_NPASS: fmt->handler = &process_npass; fmt->destroy = &destroy_npass; break;
        default: error("TODO: handler for type %d\n", fmt->type);
    }
    if ( key && fmt->type==T_INFO )
    {
        fmt->id = bcf_hdr_id2int(convert->header, BCF_DT_ID, key);
        if ( !bcf_hdr_idinfo_exists(convert->header,BCF_HL_INFO,fmt->id) )
        {
            fmt->id = -1;
            convert->undef_info_tag = strdup(key);
        }
    }
    return fmt;
}

static int parse_subscript(char **p)
{
    char *q = *p;
    if ( *q!='{' ) return -1;
    q++;
    while ( *q && *q!='}' && isdigit(*q) ) q++;
    if ( *q!='}' ) return -1;
    int idx = atoi((*p)+1);
    *p = q+1;
    return idx;
}

static char *parse_tag(convert_t *convert, char *p, int is_gtf)
{
    char *q = ++p;
    while ( *q && (isalnum(*q) || *q=='_' || *q=='.') ) q++;
    kstring_t str = {0,0,0};
    if ( q-p==0 ) error("Could not parse format string: %s\n", convert->format_str);
    kputsn(p, q-p, &str);
    if ( is_gtf )
    {
        if ( !strcmp(str.s, "SAMPLE") ) register_tag(convert, "SAMPLE", is_gtf, T_SAMPLE);
        else if ( !strcmp(str.s, "GT") ) register_tag(convert, "GT", is_gtf, T_GT);
        else if ( !strcmp(str.s, "TGT") ) register_tag(convert, "GT", is_gtf, T_TGT);
        else if ( !strcmp(str.s, "TBCSQ") )
        {
            fmt_t *fmt = register_tag(convert, "BCSQ", is_gtf, T_TBCSQ);
            fmt->subscript = parse_subscript(&q);
            if ( fmt->subscript==-1 )
            {
                if ( !strncmp(q,"{*}",3) ) { fmt->subscript = 0; q += 3; }
            }
            else fmt->subscript++;
        }
        else if ( !strcmp(str.s, "IUPACGT") ) register_tag(convert, "GT", is_gtf, T_IUPAC_GT);
        else if ( !strcmp(str.s, "INFO") )
        {
            if ( *q!='/' )
            {
                int id = bcf_hdr_id2int(convert->header, BCF_DT_ID, str.s);
                if ( bcf_hdr_idinfo_exists(convert->header,BCF_HL_INFO,id) )
                    error("Could not parse format string \"%s\". Did you mean %%INFO/%s?\n", convert->format_str,str.s);
                else
                    error("Could not parse format string: %s\n", convert->format_str);
            }
            p = ++q;
            str.l = 0;
            while ( *q && (isalnum(*q) || *q=='_' || *q=='.') ) q++;
            if ( q-p==0 ) error("Could not parse format string: %s\n", convert->format_str);
            kputsn(p, q-p, &str);
            fmt_t *fmt = register_tag(convert, str.s, is_gtf, T_INFO);
            fmt->subscript = parse_subscript(&q);
            _used_tags_add(convert,T_INFO,str.s);
        }
        else if ( !strcmp(str.s,"PBINOM") )
        {
            if ( *q!='(' ) error("Could not parse the expression: %s\n", convert->format_str);
            p = ++q;
            str.l = 0;
            while ( *q && *q!=')' ) q++;
            if ( q-p==0 ) error("Could not parse format string: %s\n", convert->format_str);
            kputsn(p, q-p, &str);
            register_tag(convert, str.s, is_gtf, T_PBINOM);
            q++;
        }
        else if ( !strcmp(str.s,"N_PASS") )
            error("N_PASS() must be placed outside the square brackets\n");
        else
        {
            fmt_t *fmt = register_tag(convert, str.s, is_gtf, T_FORMAT);
            fmt->subscript = parse_subscript(&q);
        }
    }
    else
    {
        _SET_NON_FORMAT_TAGS(register_tag, str.s, convert, str.s, is_gtf)
        else if ( !strcmp(str.s, "ALT") )
        {
            fmt_t *fmt = register_tag(convert, str.s, is_gtf, T_ALT);
            fmt->subscript = parse_subscript(&q);
        }
        else if ( !strcmp(str.s, "_CHROM_POS_ID") ) register_tag(convert, str.s, is_gtf, T_CHROM_POS_ID);
        else if ( !strcmp(str.s, "_GT_TO_PROB3") ) register_tag(convert, str.s, is_gtf, T_GT_TO_PROB3);
        else if ( !strcmp(str.s, "_PL_TO_PROB3") ) register_tag(convert, str.s, is_gtf, T_PL_TO_PROB3);
        else if ( !strcmp(str.s, "_GP_TO_PROB3") ) register_tag(convert, str.s, is_gtf, T_GP_TO_PROB3);
        else if ( !strcmp(str.s, "_GT_TO_HAP") ) register_tag(convert, str.s, is_gtf, T_GT_TO_HAP);
        else if ( !strcmp(str.s, "_GT_TO_HAP2") ) register_tag(convert, str.s, is_gtf, T_GT_TO_HAP2);
        else if ( !strcmp(str.s, "RSX") ) register_tag(convert, str.s, is_gtf, T_RSX);
        else if ( !strcmp(str.s, "VKX") ) register_tag(convert, str.s, is_gtf, T_VKX);
        else if ( !strcmp(str.s,"PBINOM") ) error("Error: PBINOM() is currently supported only with FORMAT tags. (todo)\n");
        else if ( !strcmp(str.s, "INFO") )
        {
            if ( *q=='/' )
            {
                p = ++q;
                str.l = 0;
                while ( *q && (isalnum(*q) || *q=='_' || *q=='.') ) q++;
                if ( q-p==0 ) error("Could not parse format string: %s\n", convert->format_str);
                kputsn(p, q-p, &str);
                fmt_t *fmt = register_tag(convert, str.s, is_gtf, T_INFO);
                fmt->subscript = parse_subscript(&q);
                _used_tags_add(convert,T_INFO,str.s);
            }
            else
                register_tag(convert, NULL, is_gtf, T_INFO);    // the whole INFO
        }
        else if ( !strcmp(str.s, "FORMAT") )
             register_tag(convert, NULL, 0, T_FORMAT);
        else if ( !strcmp(str.s,"N_PASS") )
        {
            if ( *q!='(' ) error("Could not parse the expression: %s\n", convert->format_str);
            p = ++q;
            str.l = 0;
            int nopen = 1;
            while ( *q && nopen )
            {
                if ( *q=='(' ) nopen++;
                else if ( *q==')' ) nopen--;
                q++;
            }
            if ( q-p==0 || nopen ) error("Could not parse format string: %s\n", convert->format_str);
            kputsn(p, q-p-1, &str);
            register_tag(convert, str.s, is_gtf, T_NPASS);
        }
        else
        {
            fmt_t *fmt = register_tag(convert, str.s, is_gtf, T_INFO);
            fmt->subscript = parse_subscript(&q);
            _used_tags_add(convert,T_INFO,str.s);
        }
    }
    free(str.s);
    return q;
}

static char *parse_sep(convert_t *convert, char *p, int is_gtf)
{
    char *q = p;
    kstring_t str = {0,0,0};
    while ( *q && *q!='[' && *q!=']' && *q!='%' )
    {
        if ( *q=='\\' )
        {
            q++;
            if ( *q=='n' ) kputc('\n', &str);
            else if ( *q=='t' ) kputc('\t', &str);
            else kputc(*q, &str);
        }
        else kputc(*q, &str);
        q++;
    }
    if ( !str.l ) error("Could not parse format string: %s\n", convert->format_str);
    register_tag(convert, str.s, is_gtf, T_SEP);
    free(str.s);
    return q;
}

convert_t *convert_init(bcf_hdr_t *hdr, int *samples, int nsamples, const char *format_str)
{
    convert_t *convert = (convert_t*) calloc(1,sizeof(convert_t));
    convert->header = hdr;
    convert->format_str = strdup(format_str);
    convert->max_unpack = BCF_UN_STR;
    convert->used_tags_hash = khash_str2int_init();

    int i, is_gtf = 0;
    char *p = convert->format_str;
    while ( *p )
    {
        //fprintf(stderr,"<%s>\n", p);
        switch (*p)
        {
            case '[': is_gtf = 1; p++; break;
            case ']': is_gtf = 0; register_tag(convert, NULL, 0, T_SEP); p++; break;
            case '%': p = parse_tag(convert, p, is_gtf); break;
            default:  p = parse_sep(convert, p, is_gtf); break;
        }
    }
    if ( is_gtf )
        error("Could not parse the format string, missing the square bracket \"]\": %s\n", convert->format_str);

    if ( nsamples )
    {
        convert->nsamples = nsamples;
        convert->samples = (int*) malloc(sizeof(int)*nsamples);
        for (i=0; i<convert->nsamples; i++) convert->samples[i] = samples[i];
    }
    else
    {
        convert->nsamples = bcf_hdr_nsamples(convert->header);
        convert->samples = (int*) malloc(sizeof(int)*convert->nsamples);
        for (i=0; i<convert->nsamples; i++) convert->samples[i] = i;
    }
    return convert;
}

void convert_destroy(convert_t *convert)
{
    int i;
    for (i=0; i<convert->nfmt; i++)
    {
        if ( convert->fmt[i].destroy ) convert->fmt[i].destroy(convert->fmt[i].usr);
        free(convert->fmt[i].key);
    }
    if ( convert->nused_tags )
    {
        for (i=0; i<convert->nused_tags; i++) free(convert->used_tags_list[i]);
        free(convert->used_tags_list);
    }
    khash_str2int_destroy(convert->used_tags_hash);
    free(convert->print_filtered);
    free(convert->fmt);
    free(convert->undef_info_tag);
    free(convert->dat);
    free(convert->samples);
    free(convert->format_str);
    free(convert);
}


int convert_header(convert_t *convert, kstring_t *str)
{
    int i, icol = 0, l_ori = str->l;
    bcf_hdr_t *hdr = convert->header;

    // Suppress the header output if LINE is present
    for (i=0; i<convert->nfmt; i++)
        if ( convert->fmt[i].type == T_LINE ) break;
    if ( i!=convert->nfmt )
        return str->l - l_ori;

    // Header formatting becomes problematic when the formatting expression contains a newline.
    // Simple cases like
    //      -f'[%CHROM %POS %SAMPLE\n]'
    // can be handled quite easily with has_fmt_newline. Note this will not work if multiple newlines
    // are present.
    int has_fmt_newline = 0;
    kputc('#', str);
    for (i=0; i<convert->nfmt; i++)
    {
        // Genotype fields
        if ( convert->fmt[i].is_gt_field )
        {
            int j = i, js, k;
            while ( convert->fmt[j].is_gt_field ) j++;
            for (js=0; js<convert->nsamples; js++)
            {
                int ks = convert->samples[js];
                for (k=i; k<j; k++)
                {
                    if ( convert->fmt[k].type == T_SEP )
                    {
                        if ( convert->fmt[k].key )
                        {
                            char *tmp = convert->fmt[k].key;
                            while ( *tmp )
                            {
                                if ( *tmp=='\n' ) has_fmt_newline = 1;
                                else kputc(*tmp,str);
                                tmp++;
                            }
                        }
                    }
                    else if ( convert->header_samples )
                    {
                        icol++;
                        if ( !convert->no_hdr_indices ) ksprintf(str,"[%d]",icol);
                        ksprintf(str,"%s:%s", hdr->samples[ks], convert->fmt[k].key);
                    }
                    else
                    {
                        icol++;
                        if ( !convert->no_hdr_indices ) ksprintf(str,"[%d]",icol);
                        ksprintf(str,"%s", convert->fmt[k].key);
                    }
                }
                if ( has_fmt_newline )
                {
                    if ( !convert->header_samples ) break;

                    // this is unfortunate: the formatting expression breaks the per-sample output into separate lines,
                    // therefore including a sample name in the header makes no sense anymore
                    convert->header_samples = 0;
                    str->l = l_ori;
                    return convert_header(convert, str);
                }
            }
            i = j-1;
            continue;
        }
        // Fixed fields
        if ( convert->fmt[i].type == T_SEP )
        {
            if ( convert->fmt[i].key ) kputs(convert->fmt[i].key, str);
            continue;
        }
        icol++;
        if ( !convert->no_hdr_indices ) ksprintf(str,"[%d]",icol);
        ksprintf(str,"%s", convert->fmt[i].key);
    }
    if ( has_fmt_newline ) kputc('\n',str);
    return str->l - l_ori;
}

int convert_line(convert_t *convert, bcf1_t *line, kstring_t *str)
{
    if ( !convert->allow_undef_tags && convert->undef_info_tag )
    {
        kstring_t msg = {0,0,0};
        ksprintf(&msg,"Error: no such tag defined in the VCF header: INFO/%s", convert->undef_info_tag);

        int hdr_id = bcf_hdr_id2int(convert->header,BCF_DT_ID,convert->undef_info_tag);
        if ( hdr_id>=0 && bcf_hdr_idinfo_exists(convert->header,BCF_HL_FMT,hdr_id) )
            ksprintf(&msg,". FORMAT fields must be enclosed in square brackets, e.g. \"[ %%%s]\"", convert->undef_info_tag);
        error("%s\n", msg.s);
    }

    int l_ori = str->l;
    bcf_unpack(line, convert->max_unpack);

    int i, ir;
    str->l = 0;
    for (i=0; i<convert->nfmt; i++)
    {
        // Genotype fields.
        if ( convert->fmt[i].is_gt_field )
        {
            int j = i, js, k;
            while ( j<convert->nfmt && convert->fmt[j].is_gt_field )
            {
                convert->fmt[j].ready = 0;
                j++;
            }
            for (js=0; js<convert->nsamples; js++)
            {
                // Skip samples when filtering was requested
                int ks = convert->samples[js];
                if ( convert->subset_samples && *convert->subset_samples && !(*convert->subset_samples)[ks] )
                {
                    if ( !convert->print_filtered ) continue;

                    for (k=i; k<j; k++)
                        if ( convert->fmt[k].type==T_SEP )
                            convert->fmt[k].handler(convert, line, &convert->fmt[k], ks, str);
                        else
                            kputs(convert->print_filtered, str);
                    continue;
                }

                // Here comes a hack designed for TBCSQ. When running on large files,
                // such as 1000GP, there are too many empty fields in the output and
                // it's very very slow. Therefore in case the handler does not add
                // anything to the string, we trim all genotype fields enclosed in square
                // brackets here. This may be changed in future, time will show...
                size_t l_start = str->l;

                for (k=i; k<j; k++)
                {
                    if ( convert->fmt[k].type == T_MASK )
                    {
                        for (ir=0; ir<convert->nreaders; ir++)
                            kputc(bcf_sr_has_line(convert->readers,ir)?'1':'0', str);
                    }
                    else if ( convert->fmt[k].handler )
                    {
                        size_t l = str->l;
                        convert->fmt[k].handler(convert, line, &convert->fmt[k], ks, str);
                        if ( l==str->l ) { str->l = l_start; break; }  // only TBCSQ does this
                    }
                }
            }
            i = j-1;
            continue;
        }
        // Fixed fields
        if ( convert->fmt[i].type == T_MASK )
        {
            for (ir=0; ir<convert->nreaders; ir++)
                kputc(bcf_sr_has_line(convert->readers,ir)?'1':'0', str);
        }
        else if ( convert->fmt[i].handler )
            convert->fmt[i].handler(convert, line, &convert->fmt[i], -1, str);

    }
    return str->l - l_ori;
}

static void force_newline_(convert_t *convert)
{
    int i, has_newline = 0;
    for (i=0; i<convert->nfmt; i++)
    {
        if ( !convert->fmt[i].key ) continue;
        char *tmp = convert->fmt[i].key;
        while (*tmp)
        {
            if ( *tmp=='\n' ) { has_newline = 1; break; }
            tmp++;
        }
        if ( has_newline ) break;
    }
    if ( has_newline ) return;

    // A newline is not present, force it. But where to add it? Always at the end.
    //
    // Briefly, in 1.18, we considered the following automatic behavior, which for
    // per-site output it would add it at the end of the expression and for per-sample
    // output it would add it inside the square brackets:
    //           -f'%CHROM[ %SAMPLE]\n'
    //           -f'[%CHROM %SAMPLE\n]'
    //
    // However, this is an annoyance for users, as it is not entirely clear what
    // will happen unless one understands the internals well (#1969)

    register_tag(convert, "\n", 0, T_SEP);
}

int convert_set_option(convert_t *convert, enum convert_option opt, ...)
{
    int ret = 0;
    va_list args;

    va_start(args, opt);
    switch (opt)
    {
        case allow_undef_tags:
            convert->allow_undef_tags = va_arg(args, int);
            break;
        case subset_samples:
            convert->subset_samples = va_arg(args, uint8_t**);
            break;
        case header_samples:
            convert->header_samples = va_arg(args, int);
            break;
        case print_filtered:
            convert->print_filtered = strdup(va_arg(args, char*));
            break;
        case force_newline:
            convert->force_newline = va_arg(args, int);
            if ( convert->force_newline ) force_newline_(convert);
            break;
        case no_hdr_indices:
            convert->no_hdr_indices = va_arg(args, int);
            break;
        default:
            ret = -1;
    }
    va_end(args);
    return ret;
}

int convert_max_unpack(convert_t *convert)
{
    return convert->max_unpack;
}

int convert_is_tag_used(convert_t *convert, char *tag)
{
    return khash_str2int_has_key(convert->used_tags_hash, tag);
}
const char **convert_list_used_tags(convert_t *convert, int *ntags)
{
    *ntags = convert->nused_tags;
    return (const char **)convert->used_tags_list;
}

