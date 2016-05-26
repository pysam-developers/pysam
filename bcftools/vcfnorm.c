/*  vcfnorm.c -- Left-align and normalize indels.

    Copyright (C) 2013-2016 Genome Research Ltd.

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
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/faidx.h>
#include "bcftools.h"
#include "rbuf.h"

#define CHECK_REF_EXIT 0
#define CHECK_REF_WARN 1
#define CHECK_REF_SKIP 2
#define CHECK_REF_FIX  4

#define MROWS_SPLIT 1
#define MROWS_MERGE  2

// for -m+, mapping from allele indexes of a single input record
// to allele indexes of output record
typedef struct
{
    int nals, mals, *map;
}
map_t;

typedef struct
{
    char *tseq, *seq;
    int mseq;
    bcf1_t **lines, **tmp_lines, **alines, **blines, *mrow_out;
    int ntmp_lines, mtmp_lines, nalines, malines, nblines, mblines;
    map_t *maps;     // mrow map for each buffered record
    char **als;
    int mmaps, nals, mals;
    uint8_t *tmp_arr1, *tmp_arr2, *diploid;
    int ntmp_arr1, ntmp_arr2;
    kstring_t *tmp_str;
    kstring_t *tmp_als, tmp_als_str;
    int ntmp_als;
    rbuf_t rbuf;
    int buf_win;            // maximum distance between two records to consider
    int aln_win;            // the realignment window size (maximum repeat size)
    bcf_srs_t *files;       // using the synced reader only for -r option
    bcf_hdr_t *hdr;
    faidx_t *fai;
    struct { int tot, set, swap; } nref;
    char **argv, *output_fname, *ref_fname, *vcf_fname, *region, *targets;
    int argc, rmdup, output_type, n_threads, check_ref, strict_filter, do_indels;
    int nchanged, nskipped, nsplit, ntotal, mrows_op, mrows_collapse, parsimonious;
    int record_cmd_line;
}
args_t;

static inline int replace_iupac_codes(char *seq, int nseq)
{
    // Replace ambiguity codes with N for now, it awaits to be seen what the VCF spec codifies in the end
    int i, n = 0;
    for (i=0; i<nseq; i++)
    {
        char c = toupper(seq[i]);
        if ( c!='A' && c!='C' && c!='G' && c!='T' ) { seq[i] = 'N'; n++; }
    }
    return n;
}

static void fix_ref(args_t *args, bcf1_t *line)
{
    int reflen = strlen(line->d.allele[0]);
    int i, maxlen = reflen, len;
    for (i=1; i<line->n_allele; i++)
    {
        int len = strlen(line->d.allele[i]);
        if ( maxlen < len ) maxlen = len;
    }

    char *ref = faidx_fetch_seq(args->fai, (char*)bcf_seqname(args->hdr,line), line->pos, line->pos+maxlen-1, &len);
    if ( !ref ) error("faidx_fetch_seq failed at %s:%d\n", bcf_seqname(args->hdr,line),line->pos+1);
    replace_iupac_codes(ref,len);

    args->nref.tot++;

    // is the REF different?
    if ( !strncasecmp(line->d.allele[0],ref,reflen) ) { free(ref); return; }

    // is the REF allele missing or N?
    if ( reflen==1 && (line->d.allele[0][0]=='.' || line->d.allele[0][0]=='N' || line->d.allele[0][0]=='n') ) 
    { 
        line->d.allele[0][0] = ref[0]; 
        args->nref.set++; 
        free(ref);
        bcf_update_alleles(args->hdr,line,(const char**)line->d.allele,line->n_allele);
        return;
    }

    // does REF contain non-standard bases?
    if ( replace_iupac_codes(line->d.allele[0],strlen(line->d.allele[0])) )
    {
        args->nref.set++;
        bcf_update_alleles(args->hdr,line,(const char**)line->d.allele,line->n_allele);
        if ( !strncasecmp(line->d.allele[0],ref,reflen) ) { free(ref); return; }
    }

    // is it swapped?
    for (i=1; i<line->n_allele; i++)
    {
        int len = strlen(line->d.allele[i]);
        if ( !strncasecmp(line->d.allele[i],ref,len) ) break;
    }

    kstring_t str = {0,0,0};
    if ( i==line->n_allele )
    {
        // none of the alternate alleles matches the reference
        if ( line->n_allele>1 )
            args->nref.set++;
        else
            args->nref.swap++;

        kputs(line->d.allele[0],&str);
        kputc(',',&str);
        for (i=1; i<line->n_allele; i++)
        {
            kputs(line->d.allele[i],&str);
            kputc(',',&str);
        }
        kputc(ref[0],&str);
        bcf_update_alleles_str(args->hdr,line,str.s);
        str.l = 0;
    }
    else
        args->nref.swap++;
    free(ref);

    // swap the alleles
    int j;
    kputs(line->d.allele[i],&str);
    for (j=1; j<i; j++)
    {
        kputc(',',&str);
        kputs(line->d.allele[j],&str);
    }
    kputc(',',&str);
    kputs(line->d.allele[0],&str);
    for (j=i+1; j<line->n_allele; j++)
    {
        kputc(',',&str);
        kputs(line->d.allele[j],&str);
    }
    bcf_update_alleles_str(args->hdr,line,str.s);

    // swap genotypes
    int ntmp = args->ntmp_arr1 / sizeof(int32_t); // reuse tmp_arr declared as uint8_t
    int ngts = bcf_get_genotypes(args->hdr, line, &args->tmp_arr1, &ntmp);
    args->ntmp_arr1 = ntmp * sizeof(int32_t);
    int32_t *gts = (int32_t*) args->tmp_arr1;
    int ni = 0;
    for (j=0; j<ngts; j++)
    {
        if ( gts[j]==bcf_gt_unphased(0) ) { gts[j] = bcf_gt_unphased(i); ni++; }
        else if ( gts[j]==bcf_gt_phased(0) ) { gts[j] = bcf_gt_phased(i); ni++; }
        else if ( gts[j]==bcf_gt_unphased(i) ) gts[j] = bcf_gt_unphased(0);
        else if ( gts[j]==bcf_gt_phased(i) ) gts[j] = bcf_gt_phased(0);
    }
    bcf_update_genotypes(args->hdr,line,gts,ngts);

    // update AC
    int nac = bcf_get_info_int32(args->hdr, line, "AC", &args->tmp_arr1, &ntmp);
    args->ntmp_arr1 = ntmp * sizeof(int32_t);
    if ( i <= nac )
    {
        int32_t *ac = (int32_t*)args->tmp_arr1;
        ac[i-1] = ni;
        bcf_update_info_int32(args->hdr, line, "AC", ac, nac);
    }
    
    free(str.s);
}

static void fix_dup_alt(args_t *args, bcf1_t *line)
{
    // update alleles, create a mapping between old and new indexes
    hts_expand(uint8_t,line->n_allele,args->ntmp_arr1,args->tmp_arr1);
    args->tmp_arr1[0] = 0;  // ref always unchanged

    int i, j, nals = line->n_allele, nals_ori = line->n_allele;
    for (i=1, j=1; i<line->n_allele; i++)
    {
        if ( strcmp(line->d.allele[0],line->d.allele[i]) )
        {
            args->tmp_arr1[i] = j++;
            continue;
        }
        args->tmp_arr1[i] = 0;
        nals--;
    }
    for (i=1, j=1; i<line->n_allele; i++)
    {
        if ( !args->tmp_arr1[i] ) continue;
        line->d.allele[j++] = line->d.allele[i];
    }
    bcf_update_alleles(args->hdr, line, (const char**)line->d.allele, nals);


    // update genotypes
    int ntmp = args->ntmp_arr2 / sizeof(int32_t); // reuse tmp_arr declared as uint8_t
    int ngts = bcf_get_genotypes(args->hdr, line, &args->tmp_arr2, &ntmp);
    args->ntmp_arr2 = ntmp * sizeof(int32_t);
    int32_t *gts = (int32_t*) args->tmp_arr2;
    int changed = 0;
    for (i=0; i<ngts; i++)
    {
        if ( bcf_gt_is_missing(gts[i]) || gts[i]==bcf_int32_vector_end ) continue;
        int ial = bcf_gt_allele(gts[i]);
        if ( ial<nals_ori && ial==args->tmp_arr1[ial] ) continue;
        int ial_new = ial<nals_ori ? args->tmp_arr1[ial] : 0;
        gts[i] = bcf_gt_is_phased(gts[i]) ? bcf_gt_phased(ial_new) : bcf_gt_unphased(ial_new);
        changed = 1;
    }
    if ( changed ) bcf_update_genotypes(args->hdr,line,gts,ngts);
}

#define ERR_DUP_ALLELE      -2
#define ERR_REF_MISMATCH    -1
#define ERR_OK              0
#define ERR_SYMBOLIC        1

static int realign(args_t *args, bcf1_t *line)
{
    bcf_unpack(line, BCF_UN_STR);

    // Sanity check REF
    int i, nref, reflen = strlen(line->d.allele[0]);
    char *ref = faidx_fetch_seq(args->fai, (char*)args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos, line->pos+reflen-1, &nref);
    if ( !ref ) error("faidx_fetch_seq failed at %s:%d\n", args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos+1);
    replace_iupac_codes(ref,nref);

    // does REF contain non-standard bases?
    if ( replace_iupac_codes(line->d.allele[0],reflen) )
    {
        args->nchanged++;
        bcf_update_alleles(args->hdr,line,(const char**)line->d.allele,line->n_allele);
    }
    if ( strcasecmp(ref,line->d.allele[0]) )
    {
        if ( args->check_ref==CHECK_REF_EXIT )
            error("Reference allele mismatch at %s:%d .. REF_SEQ:'%s' vs VCF:'%s'\n", bcf_seqname(args->hdr,line),line->pos+1,ref,line->d.allele[0]);
        if ( args->check_ref & CHECK_REF_WARN )
            fprintf(stderr,"REF_MISMATCH\t%s\t%d\t%s\n", bcf_seqname(args->hdr,line),line->pos+1,line->d.allele[0]);
        free(ref);
        return ERR_REF_MISMATCH;
    }
    free(ref);
    ref = NULL;

    if ( line->n_allele == 1 ) return ERR_OK;    // a REF

    // make a copy of each allele for trimming
    hts_expand0(kstring_t,line->n_allele,args->ntmp_als,args->tmp_als);
    kstring_t *als = args->tmp_als;
    for (i=0; i<line->n_allele; i++)
    {
        if ( line->d.allele[i][0]=='<' ) return ERR_SYMBOLIC;  // symbolic allele

        als[i].l = 0;
        kputs(line->d.allele[i], &als[i]);

        if ( i>0 && als[i].l==als[0].l && !strcasecmp(als[0].s,als[i].s) ) return ERR_DUP_ALLELE;
    }

    // trim from right
    int ori_pos = line->pos;
    while (1)
    {
        // is the rightmost base identical in all alleles?
        int min_len = als[0].l;
        for (i=1; i<line->n_allele; i++)
        {
            if ( als[0].s[ als[0].l-1 ]!=als[i].s[ als[i].l-1 ] ) break;
            if ( als[i].l < min_len ) min_len = als[i].l;
        }
        if ( i!=line->n_allele ) break; // there are differences, cannot be trimmed
        if ( min_len<=1 && line->pos==0 ) break;

        int pad_from_left = 0;
        for (i=0; i<line->n_allele; i++) // trim all alleles
        {
            als[i].l--;
            if ( !als[i].l ) pad_from_left = 1;
        }
        if ( pad_from_left )
        {
            int npad = line->pos >= args->aln_win ? args->aln_win : line->pos;
            free(ref);
            ref = faidx_fetch_seq(args->fai, (char*)args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos-npad, line->pos-1, &nref);
            if ( !ref ) error("faidx_fetch_seq failed at %s:%d\n", args->hdr->id[BCF_DT_CTG][line->rid].key, line->pos-npad+1);
            replace_iupac_codes(ref,nref);
            for (i=0; i<line->n_allele; i++)
            {
                ks_resize(&als[i], als[i].l + npad);
                if ( als[i].l ) memmove(als[i].s+npad,als[i].s,als[i].l);
                memcpy(als[i].s,ref,npad);
                als[i].l += npad;
            }
            line->pos -= npad;
        }
    }
    free(ref);

    // trim from left
    int ntrim_left = 0;
    while (1)
    {
        // is the first base identical in all alleles?
        int min_len = als[0].l - ntrim_left;
        for (i=1; i<line->n_allele; i++)
        {
            if ( als[0].s[ntrim_left]!=als[i].s[ntrim_left] ) break;
            if ( min_len > als[i].l - ntrim_left ) min_len = als[i].l - ntrim_left;
        }
        if ( i!=line->n_allele || min_len<=1 ) break; // there are differences, cannot be trimmed
        ntrim_left++;
    }
    if ( ntrim_left )
    {
        for (i=0; i<line->n_allele; i++)
        {
            memmove(als[i].s,als[i].s+ntrim_left,als[i].l-ntrim_left);
            als[i].l -= ntrim_left;
        }
        line->pos += ntrim_left;
    }

    // Have the alleles changed?
    als[0].s[ als[0].l ] = 0;  // in order for strcmp to work
    if ( ori_pos==line->pos && !strcasecmp(line->d.allele[0],als[0].s) ) return ERR_OK;

    // Create new block of alleles and update
    args->tmp_als_str.l = 0;
    for (i=0; i<line->n_allele; i++)
    {
        if (i>0) kputc(',',&args->tmp_als_str);
        kputsn(als[i].s,als[i].l,&args->tmp_als_str);
    }
    args->tmp_als_str.s[ args->tmp_als_str.l ] = 0;
    bcf_update_alleles_str(args->hdr,line,args->tmp_als_str.s);
    args->nchanged++;

    return ERR_OK;
}

static void split_info_numeric(args_t *args, bcf1_t *src, bcf_info_t *info, int ialt, bcf1_t *dst)
{
    #define BRANCH_NUMERIC(type,type_t) \
    { \
        const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key); \
        int ntmp = args->ntmp_arr1 / sizeof(type_t); \
        int ret = bcf_get_info_##type(args->hdr,src,tag,&args->tmp_arr1,&ntmp); \
        args->ntmp_arr1 = ntmp * sizeof(type_t); \
        assert( ret>0 ); \
        type_t *vals = (type_t*) args->tmp_arr1; \
        int len = bcf_hdr_id2length(args->hdr,BCF_HL_INFO,info->key); \
        if ( len==BCF_VL_A ) \
        { \
            assert( ret==src->n_allele-1); \
            bcf_update_info_##type(args->hdr,dst,tag,vals+ialt,1); \
        } \
        else if ( len==BCF_VL_R ) \
        { \
            assert( ret==src->n_allele); \
            if ( ialt!=0 ) vals[1] = vals[ialt+1]; \
            bcf_update_info_##type(args->hdr,dst,tag,vals,2); \
        } \
        else if ( len==BCF_VL_G ) \
        { \
            assert( ret==src->n_allele*(src->n_allele+1)/2 ); \
            if ( ialt!=0 ) \
            { \
                vals[1] = vals[bcf_alleles2gt(0,ialt+1)]; \
                vals[2] = vals[bcf_alleles2gt(ialt+1,ialt+1)]; \
            } \
            bcf_update_info_##type(args->hdr,dst,tag,vals,3); \
        } \
        else \
            bcf_update_info_##type(args->hdr,dst,tag,vals,ret); \
    }
    switch (bcf_hdr_id2type(args->hdr,BCF_HL_INFO,info->key))
    {
        case BCF_HT_INT:  BRANCH_NUMERIC(int32, int32_t); break;
        case BCF_HT_REAL: BRANCH_NUMERIC(float, float); break;
    }
    #undef BRANCH_NUMERIC
}
// Find n-th field in a comma-separated list and move it to dst.
// The memory areas may overlap.
#define STR_MOVE_NTH(dst,src,end,nth,len) \
{ \
    char *ss = src, *se = src; \
    int j = 0; \
    while ( *se && se<(end) ) \
    { \
        if ( *se==',' ) \
        { \
            if ( j==nth ) break; \
            j++; \
            ss = se+1; \
        } \
        se++; \
    } \
    if ( j==nth ) \
    { \
        int n = se - ss; \
        memmove((dst),ss,n); \
        src = se; \
        len += n; \
    } \
    else len = -1; \
}
static void split_info_string(args_t *args, bcf1_t *src, bcf_info_t *info, int ialt, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key);
    int ret = bcf_get_info_string(args->hdr,src,tag,&args->tmp_arr1,&args->ntmp_arr1);
    assert( ret>0 );

    kstring_t str;
    str.m = args->ntmp_arr1;
    str.l = ret;
    str.s = (char*) args->tmp_arr1;

    int len = bcf_hdr_id2length(args->hdr,BCF_HL_INFO,info->key);
    if ( len==BCF_VL_A )
    {
        char *tmp = str.s;
        int len = 0;
        STR_MOVE_NTH(str.s,tmp,str.s+str.l,ialt,len);
        if ( len<0 ) return;   // wrong number of fields: skip
        str.s[len] = 0;
        bcf_update_info_string(args->hdr,dst,tag,str.s);
    }
    else if ( len==BCF_VL_R )
    {
        char *tmp = str.s;
        int len = 0;
        STR_MOVE_NTH(str.s,tmp,str.s+str.l,0,len);
        str.s[len]=','; tmp++; len++;
        STR_MOVE_NTH(&str.s[len],tmp,str.s+str.l,ialt,len);
        if ( len<0 ) return;   // wrong number of fields: skip
        str.s[len] = 0;
        bcf_update_info_string(args->hdr,dst,tag,str.s);
    }
    else if ( len==BCF_VL_G )
    {
        int i0a = bcf_alleles2gt(0,ialt+1), iaa = bcf_alleles2gt(ialt+1,ialt+1);
        char *tmp = str.s;
        int len = 0;
        STR_MOVE_NTH(str.s,tmp,str.s+str.l,0,len);
        str.s[len]=','; tmp++; len++;
        STR_MOVE_NTH(&str.s[len],tmp,str.s+str.l,i0a-1,len);
        if ( len<0 ) return;   // wrong number of fields: skip
        str.s[len]=','; tmp++; len++;
        STR_MOVE_NTH(&str.s[len],tmp,str.s+str.l,iaa-i0a-1,len);
        if ( len<0 ) return;   // wrong number of fields: skip
        str.s[len] = 0;
        bcf_update_info_string(args->hdr,dst,tag,str.s);
    }
    else
        bcf_update_info_string(args->hdr,dst,tag,str.s);
}
static void split_info_flag(args_t *args, bcf1_t *src, bcf_info_t *info, int ialt, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key);
    int ret = bcf_get_info_flag(args->hdr,src,tag,&args->tmp_arr1,&args->ntmp_arr1);
    bcf_update_info_flag(args->hdr,dst,tag,NULL,ret);
}

static void split_format_genotype(args_t *args, bcf1_t *src, bcf_fmt_t *fmt, int ialt, bcf1_t *dst)
{
    int ntmp = args->ntmp_arr1 / 4;
    int ngts = bcf_get_genotypes(args->hdr,src,&args->tmp_arr1,&ntmp);
    args->ntmp_arr1 = ntmp * 4;
    assert( ngts >0 );

    int32_t *gt = (int32_t*) args->tmp_arr1;
    int i, j, nsmpl = bcf_hdr_nsamples(args->hdr);
    ngts /= nsmpl;
    for (i=0; i<nsmpl; i++)
    {
        for (j=0; j<ngts; j++)
        {
            if ( gt[j]==bcf_int32_vector_end ) break;
            if ( bcf_gt_is_missing(gt[j]) || bcf_gt_allele(gt[j])==0 ) continue; // missing allele or ref: leave as is
            if ( bcf_gt_allele(gt[j])==ialt+1 )
                gt[j] = bcf_gt_unphased(1) | bcf_gt_is_phased(gt[j]); // set to first ALT
            else
                gt[j] = bcf_gt_unphased(0) | bcf_gt_is_phased(gt[j]); // set to REF
        }
        gt += ngts;
    }
    bcf_update_genotypes(args->hdr,dst,args->tmp_arr1,ngts*nsmpl);
}
static void split_format_numeric(args_t *args, bcf1_t *src, bcf_fmt_t *fmt, int ialt, bcf1_t *dst)
{
    #define BRANCH_NUMERIC(type,type_t,is_vector_end,set_vector_end) \
    { \
        const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id); \
        int ntmp = args->ntmp_arr1 / sizeof(type_t); \
        int nvals = bcf_get_format_##type(args->hdr,src,tag,&args->tmp_arr1,&ntmp); \
        args->ntmp_arr1 = ntmp * sizeof(type_t); \
        assert( nvals>0 ); \
        type_t *vals = (type_t *) args->tmp_arr1; \
        int len = bcf_hdr_id2length(args->hdr,BCF_HL_FMT,fmt->id); \
        int i, nsmpl = bcf_hdr_nsamples(args->hdr); \
        if ( nvals==nsmpl ) /* all values are missing */ \
        { \
            bcf_update_format_##type(args->hdr,dst,tag,vals,nsmpl); \
            return; \
        } \
        if ( len==BCF_VL_A ) \
        { \
            assert( nvals==(src->n_allele-1)*nsmpl); \
            nvals /= nsmpl; \
            type_t *src_vals = vals, *dst_vals = vals; \
            for (i=0; i<nsmpl; i++) \
            { \
                dst_vals[0] = src_vals[ialt]; \
                dst_vals += 1; \
                src_vals += nvals; \
            } \
            bcf_update_format_##type(args->hdr,dst,tag,vals,nsmpl); \
        } \
        else if ( len==BCF_VL_R ) \
        { \
            assert( nvals==src->n_allele*nsmpl); \
            nvals /= nsmpl; \
            type_t *src_vals = vals, *dst_vals = vals; \
            for (i=0; i<nsmpl; i++) \
            { \
                dst_vals[0] = src_vals[0]; \
                dst_vals[1] = src_vals[ialt+1]; \
                dst_vals += 2; \
                src_vals += nvals; \
            } \
            bcf_update_format_##type(args->hdr,dst,tag,vals,nsmpl*2); \
        } \
        else if ( len==BCF_VL_G ) \
        { \
            if ( nvals!=src->n_allele*(src->n_allele+1)/2*nsmpl && nvals!=src->n_allele*nsmpl ) \
                error("Error at %s:%d, the tag %s has wrong number of fields\n", bcf_seqname(args->hdr,src),src->pos+1,bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id)); \
            nvals /= nsmpl; \
            int all_haploid = nvals==src->n_allele ? 1 : 0; \
            type_t *src_vals = vals, *dst_vals = vals; \
            for (i=0; i<nsmpl; i++) \
            { \
                int haploid = all_haploid; \
                if ( !haploid ) \
                { \
                    int j; \
                    for (j=0; j<nvals; j++) if ( is_vector_end ) break; \
                    if ( j!=nvals ) haploid = 1; \
                } \
                dst_vals[0] = src_vals[0]; \
                if ( haploid ) \
                { \
                    dst_vals[1] = src_vals[ialt+1]; \
                    if ( !all_haploid ) set_vector_end; \
                } \
                else \
                { \
                    dst_vals[1] = src_vals[bcf_alleles2gt(0,ialt+1)]; \
                    dst_vals[2] = src_vals[bcf_alleles2gt(ialt+1,ialt+1)]; \
                } \
                dst_vals += all_haploid ? 2 : 3; \
                src_vals += nvals; \
            } \
            bcf_update_format_##type(args->hdr,dst,tag,vals,all_haploid ? nsmpl*2 : nsmpl*3); \
        } \
        else \
            bcf_update_format_##type(args->hdr,dst,tag,vals,nvals); \
    }
    switch (bcf_hdr_id2type(args->hdr,BCF_HL_FMT,fmt->id))
    {
        case BCF_HT_INT:  BRANCH_NUMERIC(int32, int32_t, src_vals[j]==bcf_int32_vector_end, dst_vals[2]=bcf_int32_vector_end); break;
        case BCF_HT_REAL: BRANCH_NUMERIC(float, float, bcf_float_is_vector_end(src_vals[j]), bcf_float_set_vector_end(dst_vals[2])); break;
    }
    #undef BRANCH_NUMERIC
}
static void squeeze_format_char(char *str, int src_blen, int dst_blen, int n)
{
    int i, isrc = 0, idst = 0;
    for (i=0; i<n; i++)
    {
        memmove(str+idst,str+isrc,dst_blen);
        idst += dst_blen;
        isrc += src_blen;
    }
}
static void split_format_string(args_t *args, bcf1_t *src, bcf_fmt_t *fmt, int ialt, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id);
    int ret = bcf_get_format_char(args->hdr,src,tag,&args->tmp_arr1,&args->ntmp_arr1);
    assert( ret>0 );

    kstring_t str;
    str.m = args->ntmp_arr1;
    str.l = ret;
    str.s = (char*) args->tmp_arr1;

    int nsmpl = bcf_hdr_nsamples(args->hdr);
    int len = bcf_hdr_id2length(args->hdr,BCF_HL_FMT,fmt->id);
    if ( len==BCF_VL_A )
    {
        int i, blen = ret/nsmpl, maxlen = 0;
        char *ptr = str.s;
        for (i=0; i<nsmpl; i++)
        {
            char *tmp = ptr;
            int len = 0;
            STR_MOVE_NTH(tmp,tmp,ptr+blen,ialt,len);
            if ( len<0 ) return;   // wrong number of fields: skip
            if ( maxlen < len ) maxlen = len;
            ptr += blen;
        }
        if ( maxlen<blen ) squeeze_format_char(str.s,blen,maxlen,nsmpl);
        bcf_update_format_char(args->hdr,dst,tag,str.s,nsmpl*maxlen);
    }
    else if ( len==BCF_VL_R )
    {
        int i, blen = ret/nsmpl, maxlen = 0;
        char *ptr = str.s;
        for (i=0; i<nsmpl; i++)
        {
            char *tmp = ptr;
            int len = 0;
            STR_MOVE_NTH(ptr,tmp,ptr+blen,0,len);
            ptr[len]=','; tmp++; len++;
            STR_MOVE_NTH(&ptr[len],tmp,ptr+blen,ialt,len);
            if ( len<0 ) return;   // wrong number of fields: skip
            if ( maxlen < len ) maxlen = len;
            ptr += blen;
        }
        if ( maxlen<blen ) squeeze_format_char(str.s,blen,maxlen,nsmpl);
        bcf_update_format_char(args->hdr,dst,tag,str.s,nsmpl*maxlen);
    }
    else if ( len==BCF_VL_G )
    {
        int i, blen = ret/nsmpl, maxlen = 0, i0a = bcf_alleles2gt(0,ialt+1), iaa = bcf_alleles2gt(ialt+1,ialt+1);
        char *ptr = str.s;
        for (i=0; i<nsmpl; i++)
        {
            char *se = ptr, *sx = ptr+blen;
            int nfields = 1;
            while ( *se && se<sx )
            {
                if ( *se==',' ) nfields++;
                se++;
            }
            assert( nfields==src->n_allele*(src->n_allele+1)/2 || nfields==src->n_allele );
            int len = 0;
            if ( nfields==src->n_allele )   // haploid
            {
                char *tmp = ptr;
                STR_MOVE_NTH(&ptr[len],tmp,ptr+blen,0,len);
                ptr[len]=','; tmp++; len++;
                STR_MOVE_NTH(&ptr[len],tmp,ptr+blen,ialt,len);
                if ( len<0 ) return;   // wrong number of fields: skip
            }
            else    // diploid
            {
                char *tmp = ptr;
                STR_MOVE_NTH(&ptr[len],tmp,ptr+blen,0,len);
                ptr[len]=','; tmp++; len++;
                STR_MOVE_NTH(&ptr[len],tmp,ptr+blen,i0a-1,len);
                if ( len<0 ) return;   // wrong number of fields: skip
                ptr[len]=','; tmp++; len++;
                STR_MOVE_NTH(&ptr[len],tmp,ptr+blen,iaa-i0a-1,len);
                if ( len<0 ) return;   // wrong number of fields: skip
            }
            if ( maxlen < len ) maxlen = len;
            ptr += blen;
        }
        if ( maxlen<blen ) squeeze_format_char(str.s,blen,maxlen,nsmpl);
        bcf_update_format_char(args->hdr,dst,tag,str.s,nsmpl*maxlen);
    }
    else
        bcf_update_format_char(args->hdr,dst,tag,str.s,str.l);
}


static void split_multiallelic_to_biallelics(args_t *args, bcf1_t *line)
{
    int i;

    bcf_unpack(line, BCF_UN_ALL);

    // Init the target biallelic lines
    args->ntmp_lines = line->n_allele-1;
    if ( args->mtmp_lines < args->ntmp_lines )
    {
        args->tmp_lines = (bcf1_t **)realloc(args->tmp_lines,sizeof(bcf1_t*)*args->ntmp_lines);
        for (i=args->mtmp_lines; i<args->ntmp_lines; i++)
            args->tmp_lines[i] = NULL;
        args->mtmp_lines = args->ntmp_lines;
    }
    kstring_t tmp = {0,0,0};
    kputs(line->d.allele[0], &tmp);
    kputc(',', &tmp);
    int rlen  = tmp.l;
    int gt_id = bcf_hdr_id2int(args->hdr,BCF_DT_ID,"GT");
    for (i=0; i<args->ntmp_lines; i++)  // for each ALT allele
    {
        if ( !args->tmp_lines[i] ) args->tmp_lines[i] = bcf_init1();
        bcf1_t *dst = args->tmp_lines[i];
        bcf_clear(dst);

        dst->rid  = line->rid;
        dst->pos  = line->pos;
        dst->qual = line->qual;

        // Not quite sure how to handle IDs, they can be assigned to a specific
        // ALT.  For now we leave the ID unchanged for all.
        bcf_update_id(args->hdr, dst, line->d.id ? line->d.id : ".");

        tmp.l = rlen;
        kputs(line->d.allele[i+1],&tmp);
        bcf_update_alleles_str(args->hdr,dst,tmp.s);

        if ( line->d.n_flt ) bcf_update_filter(args->hdr, dst, line->d.flt, line->d.n_flt);

        int j;
        for (j=0; j<line->n_info; j++)
        {
            bcf_info_t *info = &line->d.info[j];
            int type = bcf_hdr_id2type(args->hdr,BCF_HL_INFO,info->key);
            if ( type==BCF_HT_INT || type==BCF_HT_REAL ) split_info_numeric(args, line, info, i, dst);
            else if ( type==BCF_HT_FLAG ) split_info_flag(args, line, info, i, dst);
            else split_info_string(args, line, info, i, dst);
        }

        dst->n_sample = line->n_sample;
        for (j=0; j<line->n_fmt; j++)
        {
            bcf_fmt_t *fmt = &line->d.fmt[j];
            int type = bcf_hdr_id2type(args->hdr,BCF_HL_FMT,fmt->id);
            if ( fmt->id==gt_id ) split_format_genotype(args, line, fmt, i, dst);
            else if ( type==BCF_HT_INT || type==BCF_HT_REAL ) split_format_numeric(args, line, fmt, i, dst);
            else split_format_string(args, line, fmt, i, dst);
        }
    }
    free(tmp.s);
}

// Enlarge FORMAT array containing nsmpl samples each with nals_ori values
// to accommodate nvals values for each sample, filling the gaps with missing
// values. Works also for INFO arrays, with nsmpl set to 1.
#define ENLARGE_ARRAY(type_t,set_missing,arr,narr_bytes,nsmpl,nvals_ori,nvals) \
{ \
    int nbytes_new = (nsmpl)*(nvals)*sizeof(type_t); \
    hts_expand(uint8_t,nbytes_new,narr_bytes,arr); \
    int ismpl, k; \
    for (ismpl=nsmpl-1; ismpl>=0; ismpl--) \
    { \
        type_t *dst_ptr = ((type_t*)arr) + ismpl*(nvals); \
        type_t *src_ptr = ((type_t*)arr) + ismpl*nvals_ori; \
        memmove(dst_ptr,src_ptr,sizeof(type_t)*nvals_ori); \
        for (k=nvals_ori; k<nvals; k++) set_missing; \
    } \
}
static void merge_info_numeric(args_t *args, bcf1_t **lines, int nlines, bcf_info_t *info, bcf1_t *dst)
{
    #define BRANCH_NUMERIC(type,type_t,set_missing,is_vector_end) \
    { \
        const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key); \
        int ntmp = args->ntmp_arr1 / sizeof(type_t); \
        int nvals_ori = bcf_get_info_##type(args->hdr,lines[0],tag,&args->tmp_arr1,&ntmp); \
        args->ntmp_arr1 = ntmp * sizeof(type_t); \
        assert( nvals_ori>0 ); \
        type_t *vals = (type_t*) args->tmp_arr1, *vals2; \
        int i,k,len = bcf_hdr_id2length(args->hdr,BCF_HL_INFO,info->key);  \
        if ( len==BCF_VL_A ) \
        { \
            if (nvals_ori!=lines[0]->n_allele - 1) \
                error("vcfnorm: number of fields in first record at position %s:%d for INFO tag %s not as expected [found: %d vs expected:%d]\n", bcf_seqname(args->hdr,lines[0]),lines[0]->pos+1, tag, nvals_ori, lines[0]->n_allele-1); \
            int nvals = dst->n_allele - 1; \
            ENLARGE_ARRAY(type_t,set_missing,args->tmp_arr1,args->ntmp_arr1,1,nvals_ori,nvals); \
            vals = (type_t*) args->tmp_arr1; \
            for (i=1; i<nlines; i++) \
            { \
                int ntmp2 = args->ntmp_arr2 / sizeof(type_t); \
                int nvals2 = bcf_get_info_##type(args->hdr,lines[i],tag,&args->tmp_arr2,&ntmp2); \
                if (nvals2<0) continue; /* info tag does not exist in this record, skip */ \
                args->ntmp_arr2 = ntmp2 * sizeof(type_t); \
                if (nvals2!=lines[i]->n_allele-1) \
                    error("vcfnorm: could not merge INFO tag %s at position %s:%d\n", tag, bcf_seqname(args->hdr,lines[i]),lines[i]->pos+1); \
                vals2 = (type_t*) args->tmp_arr2; \
                for (k=0; k<nvals2; k++) \
                { \
                    if ( is_vector_end ) break; \
                    vals[ args->maps[i].map[k+1] - 1 ] = vals2[k]; \
                } \
            } \
            bcf_update_info_##type(args->hdr,dst,tag,args->tmp_arr1,nvals); \
        } \
        else if ( len==BCF_VL_R ) \
        { \
            if (nvals_ori!=lines[0]->n_allele) \
                error("vcfnorm: number of fields in first record at position %s:%d for INFO tag %s not as expected [found: %d vs expected:%d]\n", bcf_seqname(args->hdr,lines[0]),lines[0]->pos+1, tag, nvals_ori, lines[0]->n_allele); \
            int nvals = dst->n_allele; \
            ENLARGE_ARRAY(type_t,set_missing,args->tmp_arr1,args->ntmp_arr1,1,nvals_ori,nvals); \
            vals = (type_t*) args->tmp_arr1; \
            for (i=1; i<nlines; i++) \
            { \
                int ntmp2 = args->ntmp_arr2 / sizeof(type_t); \
                int nvals2 = bcf_get_info_##type(args->hdr,lines[i],tag,&args->tmp_arr2,&ntmp2); \
                if (nvals2<0) continue; /* info tag does not exist in this record, skip */ \
                args->ntmp_arr2 = ntmp2 * sizeof(type_t); \
                if (nvals2!=lines[i]->n_allele) \
                    error("vcfnorm: could not merge INFO tag %s at position %s:%d\n", tag, bcf_seqname(args->hdr,lines[i]),lines[i]->pos+1); \
                vals2 = (type_t*) args->tmp_arr2; \
                for (k=0; k<nvals2; k++) \
                { \
                    if ( is_vector_end ) break; \
                    vals[ args->maps[i].map[k] ] = vals2[k]; \
                } \
            } \
            bcf_update_info_##type(args->hdr,dst,tag,args->tmp_arr1,nvals); \
        } \
        else if ( len==BCF_VL_G ) \
        { \
            /* expecting diploid gt in INFO */ \
            if (nvals_ori!=lines[0]->n_allele*(lines[0]->n_allele+1)/2) { \
                fprintf(stderr, "todo: merge Number=G INFO fields for haploid sites\n"); \
                error("vcfnorm: number of fields in first record at position %s:%d for INFO tag %s not as expected [found: %d vs expected:%d]\n", bcf_seqname(args->hdr,lines[0]),lines[0]->pos+1, tag, nvals_ori, lines[0]->n_allele*(lines[0]->n_allele+1)/2); \
            } \
            int nvals = dst->n_allele*(dst->n_allele+1)/2; \
            ENLARGE_ARRAY(type_t,set_missing,args->tmp_arr1,args->ntmp_arr1,1,nvals_ori,nvals); \
            vals = (type_t*) args->tmp_arr1; \
            for (i=1; i<nlines; i++) \
            { \
                int ntmp2 = args->ntmp_arr2 / sizeof(type_t); \
                int nvals2 = bcf_get_info_##type(args->hdr,lines[i],tag,&args->tmp_arr2,&ntmp2); \
                if (nvals2<0) continue; /* info tag does not exist in this record, skip */ \
                args->ntmp_arr2 = ntmp2 * sizeof(type_t); \
                if (nvals2!=lines[i]->n_allele*(lines[i]->n_allele+1)/2) \
                    error("vcfnorm: could not merge INFO tag %s at position %s:%d\n", tag, bcf_seqname(args->hdr,lines[i]),lines[i]->pos+1); \
                vals2 = (type_t*) args->tmp_arr2; \
                int ia,ib; \
                k = 0; \
                for (ia=0; ia<lines[i]->n_allele; ia++) \
                { \
                    for (ib=0; ib<=ia; ib++) \
                    { \
                        if ( is_vector_end ) break; \
                        int l = bcf_alleles2gt(args->maps[i].map[ia],args->maps[i].map[ib]); \
                        vals[l] = vals2[k];  \
                        k++; \
                    } \
                } \
            } \
            bcf_update_info_##type(args->hdr,dst,tag,args->tmp_arr1,nvals); \
        } \
        else \
            bcf_update_info_##type(args->hdr,dst,tag,vals,nvals_ori); \
    }
    switch (bcf_hdr_id2type(args->hdr,BCF_HL_INFO,info->key))
    {
        case BCF_HT_INT:  BRANCH_NUMERIC(int32, int32_t, dst_ptr[k]=bcf_int32_missing, vals2[k]==bcf_int32_vector_end); break;
        case BCF_HT_REAL: BRANCH_NUMERIC(float, float, bcf_float_set_missing(dst_ptr[k]), bcf_float_is_vector_end(vals2[k])); break;
    }
    #undef BRANCH_NUMERIC
}
static void merge_info_flag(args_t *args, bcf1_t **lines, int nlines, bcf_info_t *info, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key);
    int ret = bcf_get_info_flag(args->hdr,lines[0],tag,&args->tmp_arr1,&args->ntmp_arr1);
    bcf_update_info_flag(args->hdr,dst,tag,NULL,ret);
}
int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst); // see vcfmerge.c
static void merge_info_string(args_t *args, bcf1_t **lines, int nlines, bcf_info_t *info, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,info->key);

    kstring_t str;
    str.m = args->ntmp_arr1;
    str.l = 0;
    str.s = (char*) args->tmp_arr1;

    int i, j, len = bcf_hdr_id2length(args->hdr,BCF_HL_INFO,info->key);
    if ( len==BCF_VL_A || len==BCF_VL_R )
    {
        int jfrom = len==BCF_VL_A ? 1 : 0;
        kputc('.',&str);
        for (i=jfrom+1; i<dst->n_allele; i++) kputs(",.",&str);
        for (i=0; i<nlines; i++)
        {
            bcf_info_t *src = bcf_get_info(args->hdr,lines[i],tag);
            if (!src) continue;
            for (j=jfrom; j<lines[i]->n_allele; j++)
                copy_string_field((char*)src->vptr, j-jfrom, src->len, &str, args->maps[i].map[j]-jfrom);
        }
        str.s[str.l] = 0;
        args->tmp_arr1  = (uint8_t*) str.s;
        args->ntmp_arr1 = str.m;
        bcf_update_info_string(args->hdr,dst,tag,str.s);
    }
    else if ( len==BCF_VL_G )
    {
        int ngts = dst->n_allele*(dst->n_allele+1)/2;
        kputc('.',&str);
        for (i=1; i<ngts; i++) kputs(",.",&str);
        for (i=0; i<nlines; i++)
        {
            bcf_info_t *src = bcf_get_info(args->hdr,lines[i],tag);
            if (!src) continue;
            int iori, jori, kori = 0;
            for (iori=0; iori<lines[i]->n_allele; iori++)
            {
                int inew = args->maps[i].map[iori];
                for (jori=0; jori<=iori; jori++)
                {
                    int jnew = args->maps[i].map[jori];
                    int knew = bcf_alleles2gt(inew,jnew);
                    copy_string_field((char*)src->vptr,kori,src->len,&str,knew);
                    kori++;
                }
            }
        }
        str.s[str.l] = 0;
        args->tmp_arr1  = (uint8_t*) str.s;
        args->ntmp_arr1 = str.m;
        bcf_update_info_string(args->hdr,dst,tag,str.s);
    }
    else
    {
        bcf_get_info_string(args->hdr,lines[0],tag,&args->tmp_arr1,&args->ntmp_arr1);
        bcf_update_info_string(args->hdr,dst,tag,args->tmp_arr1);
    }
}
static void merge_format_genotype(args_t *args, bcf1_t **lines, int nlines, bcf_fmt_t *fmt, bcf1_t *dst)
{
    int ntmp = args->ntmp_arr1 / 4;
    int ngts = bcf_get_genotypes(args->hdr,lines[0],&args->tmp_arr1,&ntmp);
    args->ntmp_arr1 = ntmp * 4;
    assert( ngts >0 );

    int nsmpl = bcf_hdr_nsamples(args->hdr);
    ngts /= nsmpl;

    int i, j, k;
    for (i=1; i<nlines; i++)
    {
        int ntmp2 = args->ntmp_arr2 / 4;
        int ngts2 = bcf_get_genotypes(args->hdr,lines[i],&args->tmp_arr2,&ntmp2);
        args->ntmp_arr2 = ntmp2 * 4;
        ngts2 /= nsmpl;
        if ( ngts!=ngts2 ) error("Error at %s:%d: cannot combine diploid with haploid genotype\n", bcf_seqname(args->hdr,lines[i]),lines[i]->pos+1);

        int32_t *gt  = (int32_t*) args->tmp_arr1;
        int32_t *gt2 = (int32_t*) args->tmp_arr2;
        for (j=0; j<nsmpl; j++)
        {
            for (k=0; k<ngts; k++)
            {
                if ( gt2[k]==bcf_int32_vector_end ) break;
                if ( bcf_gt_is_missing(gt2[k]) || bcf_gt_allele(gt2[k])==0 ) continue;
                if ( gt2[k]==0 ) gt[k] = 0; // missing genotype
                else
                {
                    int ial = bcf_gt_allele(gt2[k]);
                    assert( ial<args->maps[i].nals );
                    gt[k] = bcf_gt_unphased( args->maps[i].map[ial] ) | bcf_gt_is_phased(gt[k]);
                }
            }
            gt  += ngts;
            gt2 += ngts;
        }
    }
    bcf_update_genotypes(args->hdr,dst,args->tmp_arr1,ngts*nsmpl);
}
static int diploid_to_haploid(int size, int nsmpl, int nals, uint8_t *vals)
{
    int i, dsrc = size*nals*(nals+1)/2, ddst = size*nals;
    uint8_t *src_ptr = vals + dsrc, *dst_ptr = vals + ddst;
    for (i=1; i<nsmpl; i++)
    {
        memmove(dst_ptr,src_ptr,ddst);
        dst_ptr += ddst;
        src_ptr += dsrc;
    }
    return nals;
}
static void merge_format_numeric(args_t *args, bcf1_t **lines, int nlines, bcf_fmt_t *fmt, bcf1_t *dst)
{
    #define BRANCH_NUMERIC(type,type_t,set_missing,is_vector_end,set_vector_end) \
    { \
        const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id); \
        int ntmp = args->ntmp_arr1 / sizeof(type_t); \
        int nvals_ori = bcf_get_format_##type(args->hdr,lines[0],tag,&args->tmp_arr1,&ntmp); \
        args->ntmp_arr1 = ntmp * sizeof(type_t); \
        assert( nvals_ori>0 ); \
        type_t *vals2, *vals = (type_t *) args->tmp_arr1; \
        int len = bcf_hdr_id2length(args->hdr,BCF_HL_FMT,fmt->id); \
        int i, j, k, nsmpl = bcf_hdr_nsamples(args->hdr); \
        nvals_ori /= nsmpl; \
        if ( len==BCF_VL_A ) \
        { \
            int nvals = dst->n_allele - 1; \
            ENLARGE_ARRAY(type_t,set_missing,args->tmp_arr1,args->ntmp_arr1,nsmpl,nvals_ori,nvals); \
            for (i=1; i<nlines; i++) \
            { \
                int ntmp2 = args->ntmp_arr2 / sizeof(type_t); \
                int nvals2 = bcf_get_format_##type(args->hdr,lines[i],tag,&args->tmp_arr2,&ntmp2); \
                if (nvals2<0) continue; /* format tag does not exist in this record, skip */ \
                args->ntmp_arr2 = ntmp2 * sizeof(type_t); \
                nvals2 /= nsmpl; \
                if (nvals2!=lines[i]->n_allele-1) \
                    error("vcfnorm: could not merge FORMAT tag %s at position %s:%d\n", tag, bcf_seqname(args->hdr,lines[i]),lines[i]->pos+1); \
                vals  = (type_t*) args->tmp_arr1; \
                vals2 = (type_t*) args->tmp_arr2; \
                for (j=0; j<nsmpl; j++) \
                { \
                    for (k=0; k<nvals2; k++) \
                    { \
                        if ( is_vector_end ) break; \
                        vals[ args->maps[i].map[k+1] - 1 ] = vals2[k]; \
                    } \
                    vals  += nvals; \
                    vals2 += nvals2; \
                } \
            } \
            bcf_update_format_##type(args->hdr,dst,tag,args->tmp_arr1,nvals*nsmpl); \
        } \
        else if ( len==BCF_VL_R ) \
        { \
            int nvals = dst->n_allele; \
            ENLARGE_ARRAY(type_t,set_missing,args->tmp_arr1,args->ntmp_arr1,nsmpl,nvals_ori,nvals); \
            for (i=1; i<nlines; i++) \
            { \
                int ntmp2 = args->ntmp_arr2 / sizeof(type_t); \
                int nvals2 = bcf_get_format_##type(args->hdr,lines[i],tag,&args->tmp_arr2,&ntmp2); \
                if (nvals2<0) continue; /* format tag does not exist in this record, skip */ \
                args->ntmp_arr2 = ntmp2 * sizeof(type_t); \
                nvals2 /= nsmpl; \
                if (nvals2!=lines[i]->n_allele) \
                    error("vcfnorm: could not merge FORMAT tag %s at position %s:%d\n", tag, bcf_seqname(args->hdr,lines[i]),lines[i]->pos+1); \
                vals  = (type_t*) args->tmp_arr1; \
                vals2 = (type_t*) args->tmp_arr2; \
                for (j=0; j<nsmpl; j++) \
                { \
                    for (k=0; k<nvals2; k++) \
                    { \
                        if ( is_vector_end ) break; \
                        vals[ args->maps[i].map[k] ] = vals2[k]; \
                    } \
                    vals  += nvals; \
                    vals2 += nvals2; \
                } \
            } \
            bcf_update_format_##type(args->hdr,dst,tag,args->tmp_arr1,nvals*nsmpl); \
        } \
        else if ( len==BCF_VL_G ) \
        { \
            /* which samples are diploid */ \
            memset(args->diploid,0,nsmpl); \
            int all_haploid = 1; \
            if ( nvals_ori > lines[0]->n_allele ) /* line possibly diploid */ \
            { \
                vals2 = (type_t*) args->tmp_arr1; \
                int ndiploid = lines[0]->n_allele*(lines[0]->n_allele+1)/2; \
                for (i=0; i<nsmpl; i++) \
                { \
                    if ( !args->diploid[i] ) \
                    { \
                        for (k=0; k<nvals_ori; k++) if ( is_vector_end ) break; \
                        if ( k==ndiploid ) { args->diploid[i] = 1; all_haploid = 0; }\
                    } \
                    vals2 += nvals_ori; \
                } \
            } \
            int nvals = dst->n_allele*(dst->n_allele+1)/2; \
            ENLARGE_ARRAY(type_t,set_missing,args->tmp_arr1,args->ntmp_arr1,nsmpl,nvals_ori,nvals); \
            for (i=1; i<nlines; i++) \
            { \
                int ntmp2 = args->ntmp_arr2 / sizeof(type_t); \
                int nvals2 = bcf_get_format_##type(args->hdr,lines[i],tag,&args->tmp_arr2,&ntmp2); \
                if (nvals2<0) continue; /* format tag does not exist in this record, skip */ \
                args->ntmp_arr2 = ntmp2 * sizeof(type_t); \
                nvals2 /= nsmpl; \
                int ndiploid = lines[i]->n_allele*(lines[i]->n_allele+1)/2; \
                int line_diploid = nvals2==ndiploid ? 1 : 0; \
                if (!(nvals2==1 || nvals2==lines[i]->n_allele || nvals2==lines[i]->n_allele*(lines[i]->n_allele+1)/2)) \
                    error("vcfnorm: could not merge FORMAT tag %s at position %s:%d\n", tag, bcf_seqname(args->hdr,lines[i]),lines[i]->pos+1); \
                vals  = (type_t*) args->tmp_arr1; \
                vals2 = (type_t*) args->tmp_arr2; \
                for (j=0; j<nsmpl; j++) \
                { \
                    int smpl_diploid = line_diploid; \
                    if ( smpl_diploid ) \
                    { \
                        for (k=0; k<nvals2; k++) if ( is_vector_end ) break; \
                        if ( k!=ndiploid ) smpl_diploid = 0; \
                    } \
                    if ( smpl_diploid && !args->diploid[j] ) { args->diploid[j] = 1; all_haploid = 0; } \
                    if ( !smpl_diploid ) \
                    { \
                        for (k=0; k<lines[i]->n_allele; k++) vals[args->maps[i].map[k]] = vals2[k]; \
                    } \
                    else \
                    { \
                        k = 0; \
                        int ia,ib; \
                        for (ia=0; ia<lines[i]->n_allele; ia++) \
                        { \
                            for (ib=0; ib<=ia; ib++) \
                            { \
                                int l = bcf_alleles2gt(args->maps[i].map[ia],args->maps[i].map[ib]); \
                                vals[l] = vals2[k]; \
                                k++; \
                            } \
                        } \
                    } \
                    vals  += nvals; \
                    vals2 += nvals2; \
                } \
            } \
            if ( all_haploid ) \
                nvals = diploid_to_haploid(sizeof(type_t),nsmpl,dst->n_allele,args->tmp_arr1); \
            else \
            {\
                k = dst->n_allele;\
                vals2 = (type_t*) args->tmp_arr1;\
                for (i=0; i<nsmpl; i++)\
                {\
                    if ( !args->diploid[i] ) set_vector_end;\
                    vals2 += nvals;\
                }\
            }\
            bcf_update_format_##type(args->hdr,dst,tag,args->tmp_arr1,nvals*nsmpl); \
        } \
        else \
            bcf_update_format_##type(args->hdr,dst,tag,args->tmp_arr1,nvals_ori*nsmpl); \
    }
    switch (bcf_hdr_id2type(args->hdr,BCF_HL_FMT,fmt->id))
    {
        case BCF_HT_INT:  BRANCH_NUMERIC(int32, int32_t, dst_ptr[k]=bcf_int32_missing, vals2[k]==bcf_int32_vector_end, vals2[k]=bcf_int32_vector_end); break;
        case BCF_HT_REAL: BRANCH_NUMERIC(float, float, bcf_float_set_missing(dst_ptr[k]), bcf_float_is_vector_end(vals2[k]), bcf_float_set_vector_end(vals2[k])); break;
    }
    #undef BRANCH_NUMERIC
}
static void merge_format_string(args_t *args, bcf1_t **lines, int nlines, bcf_fmt_t *fmt, bcf1_t *dst)
{
    const char *tag = bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id);

    int i, j, k, len = bcf_hdr_id2length(args->hdr,BCF_HL_FMT,fmt->id);
    if ( len!=BCF_VL_A && len!=BCF_VL_R && len!=BCF_VL_G )
    {
        int nret = bcf_get_format_char(args->hdr,lines[0],tag,&args->tmp_arr1,&args->ntmp_arr1);
        bcf_update_format_char(args->hdr,dst,tag,args->tmp_arr1,nret);
        return;
    }

    int nsmpl = bcf_hdr_nsamples(args->hdr);
    for (i=0; i<nsmpl; i++) args->tmp_str[i].l = 0;

    if ( len==BCF_VL_A || len==BCF_VL_R )
    {
        int jfrom = len==BCF_VL_A ? 1 : 0;
        for (i=0; i<nsmpl; i++)
        {
            kstring_t *tmp = &args->tmp_str[i];
            kputc('.',tmp);
            for (k=jfrom+1; k<dst->n_allele; k++) kputs(",.",tmp);
        }
        for (i=0; i<nlines; i++)
        {
            int nret = bcf_get_format_char(args->hdr,lines[i],tag,&args->tmp_arr1,&args->ntmp_arr1);
            if (nret<0) continue; /* format tag does not exist in this record, skip */ \
            nret /= nsmpl;
            for (k=0; k<nsmpl; k++)
            {
                kstring_t *tmp = &args->tmp_str[k];
                char *src = (char*)args->tmp_arr1 + k*nret;
                for (j=jfrom; j<lines[i]->n_allele; j++)
                    copy_string_field(src, j-jfrom, nret, tmp, args->maps[i].map[j]-jfrom);
            }
        }
    }
    else if ( len==BCF_VL_G )
    {
        hts_expand(uint8_t,nsmpl,args->ntmp_arr2,args->tmp_arr2);
        uint8_t *haploid = args->tmp_arr2;
        int nret = bcf_get_format_char(args->hdr,lines[0],tag,&args->tmp_arr1,&args->ntmp_arr1);
        nret /= nsmpl;
        for (i=0; i<nsmpl; i++)
        {
            char *ss = (char*)args->tmp_arr1 + i*nret, *se = ss+nret;
            int nfields = 1;
            while ( *ss && ss<se )
            {
                if ( *ss==',' ) nfields++;
                ss++;
            }
            if ( nfields==lines[0]->n_allele )
            {
                haploid[i] = 1;
                nfields = dst->n_allele;
            }
            else if ( nfields==lines[0]->n_allele*(lines[0]->n_allele+1)/2 )
            {
                haploid[i] = 0;
                nfields = dst->n_allele*(dst->n_allele+1)/2;
            }
            else error("The field %s at %s:%d neither diploid nor haploid?\n", tag,bcf_seqname(args->hdr,dst),dst->pos+1);

            kstring_t *tmp = &args->tmp_str[i];
            kputc('.',tmp);
            for (j=1; j<nfields; j++) kputs(",.",tmp);
        }
        for (i=0; i<nlines; i++)
        {
            if ( i ) // we already have a copy
            {
                nret = bcf_get_format_char(args->hdr,lines[i],tag,&args->tmp_arr1,&args->ntmp_arr1);
                if (nret<0) continue; /* format tag does not exist in this record, skip */ \
                nret /= nsmpl;
            }
            for (k=0; k<nsmpl; k++)
            {
                kstring_t *tmp = &args->tmp_str[k];
                char *src = (char*)args->tmp_arr1 + k*nret;
                if ( haploid[k] )
                {
                    for (j=0; j<lines[i]->n_allele; j++)
                        copy_string_field(src,j,nret, tmp, args->maps[i].map[j]);
                }
                else
                {
                    int iori, jori, kori = 0;
                    for (iori=0; iori<lines[i]->n_allele; iori++)
                    {
                        int inew = args->maps[i].map[iori];
                        for (jori=0; jori<=iori; jori++)
                        {
                            int jnew = args->maps[i].map[jori];
                            int knew = bcf_alleles2gt(inew,jnew);
                            copy_string_field(src,kori,nret,tmp,knew);
                            kori++;
                        }
                    }
                }
            }
        }
    }
    kstring_t str;
    str.m = args->ntmp_arr2;
    str.l = 0;
    str.s = (char*) args->tmp_arr2;

    int max_len = 0;
    for (i=0; i<nsmpl; i++)
        if ( max_len < args->tmp_str[i].l ) max_len = args->tmp_str[i].l;
    for (i=0; i<nsmpl; i++)
    {
        kstring_t *tmp = &args->tmp_str[i];
        kputsn(tmp->s,tmp->l,&str);
        for (j=tmp->l; j<max_len; j++) kputc('\0',&str);
    }
    args->ntmp_arr2 = str.m;
    args->tmp_arr2  = (uint8_t*)str.s;
    bcf_update_format_char(args->hdr,dst,tag,str.s,str.l);
}

char **merge_alleles(char **a, int na, int *map, char **b, int *nb, int *mb);   // see vcfmerge.c
static void merge_biallelics_to_multiallelic(args_t *args, bcf1_t *dst, bcf1_t **lines, int nlines)
{
    int i;
    for (i=0; i<nlines; i++)
        bcf_unpack(lines[i], BCF_UN_ALL);

    dst->rid  = lines[0]->rid;
    dst->pos  = lines[0]->pos;

    // take max for QUAL
    bcf_float_set_missing(dst->qual);
    for (i=0; i<nlines; i++) {
        if (bcf_float_is_missing(lines[i]->qual)) continue;
        if (bcf_float_is_missing(dst->qual) || dst->qual<lines[i]->qual)
            dst->qual = lines[i]->qual;
    }

    bcf_update_id(args->hdr, dst, lines[0]->d.id);

    // Merge and set the alleles, create a mapping from source allele indexes to dst idxs
    hts_expand0(map_t,nlines,args->mmaps,args->maps);   // a mapping for each line
    args->nals = args->maps[0].nals = lines[0]->n_allele;
    hts_expand(int,args->maps[0].nals,args->maps[0].mals,args->maps[0].map);
    hts_expand(char*,args->nals,args->mals,args->als);
    for (i=0; i<args->maps[0].nals; i++)
    {
        args->maps[0].map[i] = i;
        args->als[i] = strdup(lines[0]->d.allele[i]);
    }
    for (i=1; i<nlines; i++)
    {
        if (lines[i]->d.id[0]!='.' || lines[i]->d.id[1]) bcf_add_id(args->hdr, dst, lines[i]->d.id);
        args->maps[i].nals = lines[i]->n_allele;
        hts_expand(int,args->maps[i].nals,args->maps[i].mals,args->maps[i].map);
        args->als = merge_alleles(lines[i]->d.allele, lines[i]->n_allele, args->maps[i].map, args->als, &args->nals, &args->mals);
        if ( !args->als ) error("Failed to merge alleles at %s:%d\n", bcf_seqname(args->hdr,dst),dst->pos+1);
    }
    bcf_update_alleles(args->hdr, dst, (const char**)args->als, args->nals);
    for (i=0; i<args->nals; i++)
    {
        free(args->als[i]);
        args->als[i] = NULL;
    }

    if ( lines[0]->d.n_flt ) bcf_update_filter(args->hdr, dst, lines[0]->d.flt, lines[0]->d.n_flt);
    for (i=1; i<nlines; i++) {
        int j;
        for (j=0; j<lines[i]->d.n_flt; j++) {
            // if strict_filter, set FILTER to PASS if any site PASS
            // otherwise accumulate FILTERs
            if (lines[i]->d.flt[j] == bcf_hdr_id2int(args->hdr, BCF_DT_ID, "PASS")) {
                if (args->strict_filter) {
                    bcf_update_filter(args->hdr, dst, lines[i]->d.flt, lines[i]->d.n_flt);
                    break;
                }
                else
                    continue;
            }
            bcf_add_filter(args->hdr, dst, lines[i]->d.flt[j]);
        }
    }

    // merge info
    for (i=0; i<lines[0]->n_info; i++)
    {
        bcf_info_t *info = &lines[0]->d.info[i];
        int type = bcf_hdr_id2type(args->hdr,BCF_HL_INFO,info->key);
        if ( type==BCF_HT_INT || type==BCF_HT_REAL ) merge_info_numeric(args, lines, nlines, info, dst);
        else if ( type==BCF_HT_FLAG ) merge_info_flag(args, lines, nlines, info, dst);
        else merge_info_string(args, lines, nlines, info, dst);
    }

    // merge format
    int gt_id = bcf_hdr_id2int(args->hdr,BCF_DT_ID,"GT");
    dst->n_sample = lines[0]->n_sample;
    for (i=0; i<lines[0]->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &lines[0]->d.fmt[i];
        int type = bcf_hdr_id2type(args->hdr,BCF_HL_FMT,fmt->id);
        if ( fmt->id==gt_id ) merge_format_genotype(args, lines, nlines, fmt, dst);
        else if ( type==BCF_HT_INT || type==BCF_HT_REAL ) merge_format_numeric(args, lines, nlines, fmt, dst);
        else merge_format_string(args, lines, nlines, fmt, dst);
    }
}

#define SWAP(type_t, a, b) { type_t t = a; a = b; b = t; }
static void mrows_schedule(args_t *args, bcf1_t **line)
{
    int i,m;
    if ( args->mrows_collapse==COLLAPSE_ANY         // merge all record types together
        || bcf_get_variant_types(*line)&VCF_SNP     // SNP, put into alines
        || bcf_get_variant_types(*line)==VCF_REF )  // ref
    {
        args->nalines++;
        m = args->malines;
        hts_expand(bcf1_t*,args->nalines,args->malines,args->alines);
        for (i=m; i<args->malines; i++) args->alines[i] = bcf_init1();
        SWAP(bcf1_t*, args->alines[args->nalines-1], *line);
    }
    else
    {
        args->nblines++;
        m = args->mblines;
        hts_expand(bcf1_t*,args->nblines,args->mblines,args->blines);
        for (i=m; i<args->mblines; i++) args->blines[i] = bcf_init1();
        SWAP(bcf1_t*, args->blines[args->nblines-1], *line);
    }
}
static int mrows_ready_to_flush(args_t *args, bcf1_t *line)
{
    if ( args->nalines && (args->alines[0]->rid!=line->rid || args->alines[0]->pos!=line->pos) ) return 1;
    if ( args->nblines && (args->blines[0]->rid!=line->rid || args->blines[0]->pos!=line->pos) ) return 1;
    return 0;
}
static bcf1_t *mrows_flush(args_t *args)
{
    if ( args->nblines && args->nalines==1 && bcf_get_variant_types(args->alines[0])==VCF_REF )
    {
        // By default, REF lines are merged with SNPs if SNPs and indels are to be kept separately.
        // However, if there are indels only and a single REF line, merge it with indels.
        args->nblines++;
        int i,m = args->mblines;
        hts_expand(bcf1_t*,args->nblines,args->mblines,args->blines);
        for (i=m; i<args->mblines; i++) args->blines[i] = bcf_init1();
        SWAP(bcf1_t*, args->blines[args->nblines-1], args->alines[0]);
        args->nalines--;
    }
    if ( args->nalines )
    {
        if ( args->nalines==1 )
        {
            args->nalines = 0;
            return args->alines[0];
        }
        bcf_clear(args->mrow_out);
        merge_biallelics_to_multiallelic(args, args->mrow_out, args->alines, args->nalines);
        args->nalines = 0;
        return args->mrow_out;
    }
    else if ( args->nblines )
    {
        if ( args->nblines==1 )
        {
            args->nblines = 0;
            return args->blines[0];
        }
        bcf_clear(args->mrow_out);
        merge_biallelics_to_multiallelic(args, args->mrow_out, args->blines, args->nblines);
        args->nblines = 0;
        return args->mrow_out;
    }
    return NULL;
}
static void flush_buffer(args_t *args, htsFile *file, int n)
{
    bcf1_t *line;
    int i, k;
    for (i=0; i<n; i++)
    {
        k = rbuf_shift(&args->rbuf);
        if ( args->mrows_op==MROWS_MERGE )
        {
            if ( mrows_ready_to_flush(args, args->lines[k]) )
            {
                while ( (line=mrows_flush(args)) ) bcf_write1(file, args->hdr, line);
            }
            int merge = 1;
            if ( args->mrows_collapse!=COLLAPSE_BOTH && args->mrows_collapse!=COLLAPSE_ANY )
            {
                if ( !(bcf_get_variant_types(args->lines[k]) & args->mrows_collapse) ) merge = 0;
            }
            if ( merge )
            {
                mrows_schedule(args, &args->lines[k]);
                continue;
            }
        }
        bcf_write1(file, args->hdr, args->lines[k]);
    }
    if ( args->mrows_op==MROWS_MERGE && !args->rbuf.n )
    {
        while ( (line=mrows_flush(args)) ) bcf_write1(file, args->hdr, line);
    }
}

static void init_data(args_t *args)
{
    args->hdr = args->files->readers[0].header;
    rbuf_init(&args->rbuf, 100);
    args->lines = (bcf1_t**) calloc(args->rbuf.m, sizeof(bcf1_t*));
    if ( args->ref_fname )
    {
        args->fai = fai_load(args->ref_fname);
        if ( !args->fai ) error("Failed to load the fai index: %s\n", args->ref_fname);
    }
    if ( args->mrows_op==MROWS_MERGE )
    {
        args->mrow_out = bcf_init1();
        args->tmp_str = (kstring_t*) calloc(bcf_hdr_nsamples(args->hdr),sizeof(kstring_t));
        args->diploid = (uint8_t*) malloc(bcf_hdr_nsamples(args->hdr));
    }
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->rbuf.m; i++)
        if ( args->lines[i] ) bcf_destroy1(args->lines[i]);
    free(args->lines);
    for (i=0; i<args->mtmp_lines; i++)
        if ( args->tmp_lines[i] ) bcf_destroy1(args->tmp_lines[i]);
    free(args->tmp_lines);
    for (i=0; i<args->malines; i++)
        bcf_destroy1(args->alines[i]);
    free(args->alines);
    for (i=0; i<args->mblines; i++)
        bcf_destroy1(args->blines[i]);
    free(args->blines);
    for (i=0; i<args->mmaps; i++)
        free(args->maps[i].map);
    for (i=0; i<args->ntmp_als; i++)
        free(args->tmp_als[i].s);
    free(args->tmp_als);
    free(args->tmp_als_str.s);
    if ( args->tmp_str )
    {
        for (i=0; i<bcf_hdr_nsamples(args->hdr); i++) free(args->tmp_str[i].s);
        free(args->tmp_str);
    }
    free(args->maps);
    free(args->als);
    free(args->tmp_arr1);
    free(args->tmp_arr2);
    free(args->diploid);
    if ( args->mrow_out ) bcf_destroy1(args->mrow_out);
    if ( args->fai ) fai_destroy(args->fai);
    if ( args->mseq ) free(args->seq);
}


static void normalize_line(args_t *args, bcf1_t **line_ptr)
{
    bcf1_t *line = *line_ptr;
    if ( args->fai )
    {
        if ( args->check_ref & CHECK_REF_FIX ) fix_ref(args, line);
        if ( args->do_indels )
        {
            int ret = realign(args, line);

            // exclude broken VCF lines
            if ( ret==ERR_REF_MISMATCH && args->check_ref & CHECK_REF_SKIP )
            {
                args->nskipped++;
                return;
            }
            if ( ret==ERR_DUP_ALLELE )
            {
                if ( args->check_ref & CHECK_REF_FIX )
                    fix_dup_alt(args, line);
                else if ( args->check_ref==CHECK_REF_EXIT )
                    error("Duplicate alleles at %s:%d; run with -cw to turn the error into warning or with -cs to fix.\n", bcf_seqname(args->hdr,line),line->pos+1);
                else if ( args->check_ref & CHECK_REF_WARN )
                    fprintf(stderr,"ALT_DUP\t%s\t%d\n", bcf_seqname(args->hdr,line),line->pos+1);
            }
        }
    }

    // insert into sorted buffer
    rbuf_expand0(&args->rbuf,bcf1_t*,args->rbuf.n+1,args->lines);
    int i,j;
    i = j = rbuf_append(&args->rbuf);
    if ( !args->lines[i] ) args->lines[i] = bcf_init1();
    SWAP(bcf1_t*, (*line_ptr), args->lines[i]);
    while ( rbuf_prev(&args->rbuf,&i) )
    {
        if ( args->lines[i]->pos > args->lines[j]->pos ) SWAP(bcf1_t*, args->lines[i], args->lines[j]);
        j = i;
    }
}

static void normalize_vcf(args_t *args)
{
    htsFile *out = hts_open(args->output_fname, hts_bcf_wmode(args->output_type));
    if ( out == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    if ( args->n_threads ) hts_set_threads(out, args->n_threads);
    if (args->record_cmd_line) bcf_hdr_append_version(args->hdr, args->argc, args->argv, "bcftools_norm");
    bcf_hdr_write(out, args->hdr);

    int prev_rid = -1, prev_pos = -1, prev_type = 0;
    while ( bcf_sr_next_line(args->files) )
    {
        args->ntotal++;

        bcf1_t *line = args->files->readers[0].buffer[0];
        if ( args->rmdup )
        {
            int line_type = bcf_get_variant_types(line);
            if ( prev_rid>=0 && prev_rid==line->rid && prev_pos==line->pos )
            {
                if ( (args->rmdup>>1)&COLLAPSE_ANY ) continue;
                if ( (args->rmdup>>1)&COLLAPSE_SNPS && line_type&(VCF_SNP|VCF_MNP) && prev_type&(VCF_SNP|VCF_MNP) ) continue;
                if ( (args->rmdup>>1)&COLLAPSE_INDELS && line_type&(VCF_INDEL) && prev_type&(VCF_INDEL) ) continue;
            }
            else
            {
                prev_rid  = line->rid;
                prev_pos  = line->pos;
                prev_type = 0;
            }
            prev_type |= line_type;
        }

        // still on the same chromosome?
        int i,j,ilast = rbuf_last(&args->rbuf);
        if ( ilast>=0 && line->rid != args->lines[ilast]->rid ) flush_buffer(args, out, args->rbuf.n); // new chromosome

        int split = 0;
        if ( args->mrows_op==MROWS_SPLIT )
        {
            split = 1;
            if ( args->mrows_collapse!=COLLAPSE_BOTH && args->mrows_collapse!=COLLAPSE_ANY )
            {
                if ( !(bcf_get_variant_types(line) & args->mrows_collapse) ) split = 0;
            }
            if ( split && line->n_allele>2 )
            {
                args->nsplit++;
                split_multiallelic_to_biallelics(args, line);
                for (j=0; j<args->ntmp_lines; j++)
                    normalize_line(args, &args->tmp_lines[j]);
            }
            else
                split = 0;
        }
        if ( !split )
            normalize_line(args, &args->files->readers[0].buffer[0]);

        // find out how many sites to flush
        ilast = rbuf_last(&args->rbuf);
        j = 0;
        for (i=-1; rbuf_next(&args->rbuf,&i); )
        {
            if ( args->lines[ilast]->pos - args->lines[i]->pos < args->buf_win ) break;
            j++;
        }
        if ( j>0 ) flush_buffer(args, out, j);
    }
    flush_buffer(args, out, args->rbuf.n);
    hts_close(out);

    fprintf(stderr,"Lines   total/split/realigned/skipped:\t%d/%d/%d/%d\n", args->ntotal,args->nsplit,args->nchanged,args->nskipped);
    if ( args->check_ref & CHECK_REF_FIX )
        fprintf(stderr,"REF/ALT total/modified/added:  \t%d/%d/%d\n", args->nref.tot,args->nref.swap,args->nref.set);
}

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Left-align and normalize indels; check if REF alleles match the reference;\n");
    fprintf(stderr, "         split multiallelic sites into multiple rows; recover multiallelics from\n");
    fprintf(stderr, "         multiple rows.\n");
    fprintf(stderr, "Usage:   bcftools norm [options] <in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -c, --check-ref <e|w|x|s>         check REF alleles and exit (e), warn (w), exclude (x), or set (s) bad sites [e]\n");
    fprintf(stderr, "    -D, --remove-duplicates           remove duplicate lines of the same type.\n");
    fprintf(stderr, "    -d, --rm-dup <type>               remove duplicate snps|indels|both|any\n");
    fprintf(stderr, "    -f, --fasta-ref <file>            reference sequence\n");
    fprintf(stderr, "    -m, --multiallelics <-|+>[type]   split multiallelics (-) or join biallelics (+), type: snps|indels|both|any [both]\n");
    fprintf(stderr, "        --no-version                  do not append version and command line to the header\n");
    fprintf(stderr, "    -N, --do-not-normalize            do not normalize indels (with -m or -c s)\n");
    fprintf(stderr, "    -o, --output <file>               write output to a file [standard output]\n");
    fprintf(stderr, "    -O, --output-type <type>          'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]\n");
    fprintf(stderr, "    -r, --regions <region>            restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>         restrict to regions listed in a file\n");
    fprintf(stderr, "    -s, --strict-filter               when merging (-m+), merged site is PASS only if all sites being merged PASS\n");
    fprintf(stderr, "    -t, --targets <region>            similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file <file>         similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "        --threads <int>               number of extra output compression threads [0]\n");
    fprintf(stderr, "    -w, --site-win <int>              buffer for sorting lines which changed position during realignment [1000]\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfnorm(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->output_fname = "-";
    args->output_type = FT_VCF;
    args->n_threads = 0;
    args->record_cmd_line = 1;
    args->aln_win = 100;
    args->buf_win = 1000;
    args->mrows_collapse = COLLAPSE_BOTH;
    args->do_indels = 1;
    int region_is_file  = 0;
    int targets_is_file = 0;

    static struct option loptions[] =
    {
        {"help",no_argument,NULL,'h'},
        {"fasta-ref",required_argument,NULL,'f'},
        {"do-not-normalize",no_argument,NULL,'N'},
        {"multiallelics",required_argument,NULL,'m'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"targets",required_argument,NULL,'t'},
        {"targets-file",required_argument,NULL,'T'},
        {"site-win",required_argument,NULL,'w'},
        {"remove-duplicates",no_argument,NULL,'D'},
        {"rm-dup",required_argument,NULL,'d'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"threads",required_argument,NULL,9},
        {"check-ref",required_argument,NULL,'c'},
        {"strict-filter",no_argument,NULL,'s'},
        {"no-version",no_argument,NULL,8},
        {NULL,0,NULL,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "hr:R:f:w:Dd:o:O:c:m:t:T:sN",loptions,NULL)) >= 0) {
        switch (c) {
            case 'N': args->do_indels = 0; break;
            case 'd':
                if ( !strcmp("snps",optarg) ) args->rmdup = COLLAPSE_SNPS<<1;
                else if ( !strcmp("indels",optarg) ) args->rmdup = COLLAPSE_INDELS<<1;
                else if ( !strcmp("both",optarg) ) args->rmdup = COLLAPSE_BOTH<<1;
                else if ( !strcmp("any",optarg) ) args->rmdup = COLLAPSE_ANY<<1;
                else error("The argument to -d not recognised: %s\n", optarg);
                break;
            case 'm':
                if ( optarg[0]=='-' ) args->mrows_op = MROWS_SPLIT;
                else if ( optarg[0]=='+' ) args->mrows_op = MROWS_MERGE;
                else error("Expected '+' or '-' with -m\n");
                if ( optarg[1]!=0 )
                {
                    if ( !strcmp("snps",optarg+1) ) args->mrows_collapse = COLLAPSE_SNPS;
                    else if ( !strcmp("indels",optarg+1) ) args->mrows_collapse = COLLAPSE_INDELS;
                    else if ( !strcmp("both",optarg+1) ) args->mrows_collapse = COLLAPSE_BOTH;
                    else if ( !strcmp("any",optarg+1) ) args->mrows_collapse = COLLAPSE_ANY;
                    else error("The argument to -m not recognised: %s\n", optarg);
                }
                break;
            case 'c':
                if ( strchr(optarg,'w') ) args->check_ref |= CHECK_REF_WARN;
                if ( strchr(optarg,'x') ) args->check_ref |= CHECK_REF_SKIP;
                if ( strchr(optarg,'s') ) args->check_ref |= CHECK_REF_FIX;
                if ( strchr(optarg,'e') ) args->check_ref = CHECK_REF_EXIT; // overrides the above
                break;
            case 'O':
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                }
                break;
            case 'o': args->output_fname = optarg; break;
            case 'D':
                fprintf(stderr,"Warning: `-D` is functional but deprecated, replaced by `-d both`.\n"); 
                args->rmdup = COLLAPSE_NONE<<1;
                break;
            case 's': args->strict_filter = 1; break;
            case 'f': args->ref_fname = optarg; break;
            case 'r': args->region = optarg; break;
            case 'R': args->region = optarg; region_is_file = 1; break;
            case 't': args->targets = optarg; break;
            case 'T': args->targets = optarg; targets_is_file = 1; break;
            case 'w':
                args->buf_win = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --site-win %s\n", optarg);
                break;
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case  8 : args->record_cmd_line = 0; break;
            case 'h':
            case '?': usage();
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( argc>optind+1 ) usage();
    if ( !args->ref_fname && !args->mrows_op && !args->rmdup ) usage();
    if ( !args->ref_fname && args->check_ref&CHECK_REF_FIX ) error("Expected --fasta-ref with --check-ref s\n");
    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else usage();
    }
    else fname = argv[optind];

    if ( args->region )
    {
        if ( bcf_sr_set_regions(args->files, args->region,region_is_file)<0 )
            error("Failed to read the regions: %s\n", args->region);
    }
    if ( args->targets )
    {
        if ( bcf_sr_set_targets(args->files, args->targets,targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets);
    }

    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->files->errnum));
    if ( args->mrows_op&MROWS_SPLIT && args->rmdup ) error("Cannot combine -D and -m-\n");
    init_data(args);
    normalize_vcf(args);
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}

