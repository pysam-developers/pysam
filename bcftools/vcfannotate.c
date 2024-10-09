/*  vcfannotate.c -- Annotate and edit VCF/BCF files.

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
#include <strings.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <inttypes.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kseq.h>
#include <htslib/khash_str2int.h>
#include "bcftools.h"
#include "vcmp.h"
#include "filter.h"
#include "convert.h"
#include "smpl_ilist.h"
#include "regidx.h"
#include "dbuf.h"

struct _args_t;

typedef struct _rm_tag_t
{
    char *key;
    int hdr_id;
    void (*handler)(struct _args_t *, bcf1_t *, struct _rm_tag_t *);
}
rm_tag_t;

typedef struct
{
    char **cols;
    int ncols, mcols;
    char **als;
    int nals, mals;
    kstring_t line;
    int rid, start, end;
}
annot_line_t;

#define REPLACE_MISSING     (1<<0)   // -c +TAG  .. replace only missing values
#define REPLACE_ALL         (1<<1)   // -c TAG   .. replace both missing and existing values
#define REPLACE_NON_MISSING (1<<2)   // -c -TAG  .. replace only if tgt is not missing
#define SET_OR_APPEND       (1<<3)   // -c =TAG  .. set new value if missing or non-existent, append otherwise
#define MATCH_VALUE         (1<<4)   // -c ~ID   .. do not set, just match the value
#define CARRY_OVER_MISSING  (1<<5)   // -c .TAG  .. carry over source missing values as well
#define MM_FIRST   0    // if multiple annotation lines overlap a VCF record, use the first, discarding the rest
#define MM_APPEND  1    // append, possibly multiple times
#define MM_UNIQUE  2    // append, only unique values
#define MM_SUM     3
#define MM_AVG     4
#define MM_MIN     5
#define MM_MAX     6
#define MM_APPEND_MISSING 7     // missing values will be transferred as well
typedef struct _annot_col_t
{
    int icol, replace, number;  // number: one of BCF_VL_* types
    char *hdr_key_src, *hdr_key_dst;
    // The setters return 0 on successful update of the bcf record, negative value (bcf_update_* return status) on errors,
    // or 1 on (repeated partial updates) concluded with a src=NULL call
    int (*setter)(struct _args_t *, bcf1_t *dst, struct _annot_col_t *, void *src); // the last is the annotation line, either src bcf1_t or annot_line_t
    int (*getter)(struct _args_t *, bcf1_t *src, struct _annot_col_t *, void **ptr, int *mptr);
    int merge_method;               // one of the MM_* defines
    khash_t(str2int) *mm_str_hash;  // lookup table to ensure uniqueness of added string values
    kstring_t mm_kstr;
    size_t
        mm_dbl_nalloc,  // the allocated size --merge-logic values array
        mm_dbl_nused,   // the number of used elements in the mm_dbl array
        mm_dbl_ndat;    // the number of merged rows (for calculating the average)
    double
        *mm_dbl;
    void *ptr;
    int mptr, done;
}
annot_col_t;

typedef struct
{
    char *name;     // column name
    int ht_type;    // type, one of BCF_HT_STR,BCF_HT_INT,BCF_HT_REAL
    int icol;       // index of the annotation column to use
    union {         // memory area with the current annotation value to pass to filter_test_ext
        int i;
        float f;
        char *s;
    };
}
ext_t;

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define MARK_LISTED   1
#define MARK_UNLISTED 2

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr, *hdr_out, *tgts_hdr;
    htsFile *out_fh;
    int output_type, n_threads, clevel;
    bcf_sr_regions_t *tgts;
    char *index_fn;
    int write_index;

    regidx_t *tgt_idx;  // keep everything in memory only with .tab annotation file and -c BEG,END columns
    regitr_t *tgt_itr;
    int tgt_is_bed;

    filter_t *filter, *filter_ext;  // only one is initialized, the latter contains external values to set dynamically on the fly
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE
    int keep_sites;

    rm_tag_t *rm;           // tags scheduled for removal
    int nrm;
    int flt_keep_pass;      // when all filters removed, reset to PASS

    vcmp_t *vcmp;           // for matching annotation and VCF lines by allele
    annot_line_t *alines;   // buffered annotation lines
    annot_line_t *aline_missing;
    uint32_t *srt_alines;   // sorted indexes (iALT<<16 || iAline)
    int nalines, malines, nsrt_alines, msrt_alines;
    int ref_idx, alt_idx, chr_idx, beg_idx, end_idx;   // -1 if not present
    annot_col_t *cols;      // column indexes and setters
    int ncols;
    int match_id;           // set iff `-c ~ID` given, -1 otherwise
    int match_end;          // set iff `-c ~INFO/END` is given, -1 otherwise

    char *set_ids_fmt;
    convert_t *set_ids;
    int set_ids_replace;

    // external values for dynamic -i/-e expressions
    int n_ext;
    ext_t *ext;
    void **ext_ptr;

    int nsmpl_annot;
    int *sample_map, nsample_map, sample_is_file;   // map[idst] -> isrc
    uint8_t *src_smpl_pld, *dst_smpl_pld;   // for Number=G format fields
    int mtmpi, mtmpf, mtmps;
    int mtmpi2, mtmpf2, mtmps2;
    int mtmpi3, mtmpf3, mtmps3;
    int32_t *tmpi, *tmpi2, *tmpi3;
    float *tmpf, *tmpf2, *tmpf3;
    char *tmps, *tmps2, **tmpp, **tmpp2;
    kstring_t tmpks;

    char **argv, *output_fname, *targets_fname, *regions_list, *header_fname;
    char *remove_annots, *columns, *rename_chrs, *rename_annots, *sample_names, *mark_sites;
    char **rename_annots_map;
    char *min_overlap_str;
    float min_overlap_ann, min_overlap_vcf;
    int rename_annots_nmap;
    kstring_t merge_method_str;
    int argc, drop_header, record_cmd_line, tgts_is_vcf, mark_sites_logic, force, single_overlaps;
    int columns_is_file, has_append_mode, pair_logic;
    dbuf_t *header_lines;
    bcf1_t *current_rec;    // current record for local setters
}
args_t;

char *msprintf(const char *fmt, ...);

int parse_with_payload(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, void *payload, void *usr)
{
    args_t *args = (args_t*) usr;
    int ret = args->tgt_is_bed ? regidx_parse_bed(line, chr_beg, chr_end, beg, end, NULL, NULL) : regidx_parse_tab(line, chr_beg, chr_end, beg, end, NULL, NULL);
    if ( ret<0 ) return ret;
    *((char **)payload) = strdup(line);
    return 0;
}
void free_payload(void *payload)
{
    char *str = *((char**)payload);
    free(str);
}

void remove_id(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_id(args->hdr,line,NULL);
}
void remove_filter(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    if ( tag->key && tag->hdr_id<0 )
    {
        error("Error: Cannot proceed, not even with the --force option, bad things could happen.\n"
              "       Note that \"bcftools annotate -x FILTER\" can be used to remove ALL filters.\n"
              "       Even better, use \"bcftools view -h\" and \"bcftools reheader\" to fix the header!\n"
              );
    }
    if ( !tag->key ) bcf_update_filter(args->hdr, line, NULL, args->flt_keep_pass);
    else bcf_remove_filter(args->hdr, line, tag->hdr_id, args->flt_keep_pass);
}
void remove_qual(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_float_set_missing(line->qual);
}
void remove_info(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    // remove all INFO fields
    if ( !(line->unpacked & BCF_UN_INFO) ) bcf_unpack(line, BCF_UN_INFO);

    int i;
    for (i=0; i<line->n_info; i++)
    {
        bcf_info_t *inf = &line->d.info[i];
        if (  !strcmp("END",bcf_hdr_int2id(args->hdr,BCF_DT_ID,inf->key)) )
            line->rlen = line->n_allele ? strlen(line->d.allele[0]) : 0;
        if ( inf->vptr_free )
        {
            free(inf->vptr - inf->vptr_off);
            inf->vptr_free = 0;
        }
        line->d.shared_dirty |= BCF1_DIRTY_INF;
        inf->vptr = NULL;
        inf->vptr_off = inf->vptr_len = 0;
    }
}
void remove_info_tag(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_info(args->hdr, line, tag->key, NULL, 0, BCF_HT_INT);  // the type does not matter with n=0
}
void remove_format_tag(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_format(args->hdr, line, tag->key, NULL, 0, BCF_HT_INT);  // the type does not matter with n=0
}
void remove_format(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    // remove all FORMAT fields except GT
    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);

    int i;
    for (i=0; i<line->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &line->d.fmt[i];
        const char *key = bcf_hdr_int2id(args->hdr,BCF_DT_ID,fmt->id);
        if ( key[0]=='G' && key[1]=='T' && !key[2] ) continue;

        if ( fmt->p_free )
        {
            free(fmt->p - fmt->p_off);
            fmt->p_free = 0;
        }
        line->d.indiv_dirty = 1;
        fmt->p = NULL;
    }
}

#include "htslib/khash.h"
KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

static void remove_hdr_lines(bcf_hdr_t *hdr, int type)
{
    int i = 0, nrm = 0;
    while ( i<hdr->nhrec )
    {
        if ( hdr->hrec[i]->type!=type ) { i++; continue; }
        bcf_hrec_t *hrec = hdr->hrec[i];
        if ( type==BCF_HL_FMT || type==BCF_HL_INFO || type==BCF_HL_FMT || type== BCF_HL_CTG )
        {
            // everything except FORMAT/GT
            int id = bcf_hrec_find_key(hrec, "ID");
            if ( id>=0 )
            {
                if ( type==BCF_HL_FMT && !strcmp(hrec->vals[id],"GT") ) { i++; continue; }
                vdict_t *d = type==BCF_HL_CTG ? (vdict_t*)hdr->dict[BCF_DT_CTG] : (vdict_t*)hdr->dict[BCF_DT_ID];
                khint_t k = kh_get(vdict, d, hdr->hrec[i]->vals[id]);
                kh_val(d, k).hrec[type==BCF_HL_CTG?0:type] = NULL;
                kh_val(d, k).info[type] |= 0xf;
            }
        }
        nrm++;
        hdr->nhrec--;
        if ( i < hdr->nhrec )
            memmove(&hdr->hrec[i],&hdr->hrec[i+1],(hdr->nhrec-i)*sizeof(bcf_hrec_t*));
        bcf_hrec_destroy(hrec);
    }
    if ( nrm ) {
        if (bcf_hdr_sync(hdr) < 0)
            error_errno("[%s] Failed to update header", __func__);
    }
}

static void init_remove_annots(args_t *args)
{
    int keep_info = 0, keep_fmt = 0, keep_flt = 0;
    void *keep = khash_str2int_init();
    kstring_t str = {0,0,0};
    char *ss = args->remove_annots;

    int i, ntags, needs_info = 0;
    if ( args->set_ids )
    {
        const char **tags = convert_list_used_tags(args->set_ids,&ntags);
        for (i=0; i<ntags; i++)
            if ( !strncmp("INFO/",tags[i],4) ) needs_info = 1;
    }

    while ( *ss )
    {
        args->nrm++;
        args->rm = (rm_tag_t*) realloc(args->rm,sizeof(rm_tag_t)*args->nrm);
        rm_tag_t *tag = &args->rm[args->nrm-1];
        tag->key = NULL;

        int type = BCF_HL_GEN;
        if ( !strncasecmp("INFO/",ss,5) ) { type = BCF_HL_INFO; ss += 5; }
        else if ( !strncasecmp("INF/",ss,4) ) { type = BCF_HL_INFO; ss += 4; }
        else if ( !strncasecmp("FORMAT/",ss,7) ) { type = BCF_HL_FMT; ss += 7; }
        else if ( !strncasecmp("FMT/",ss,4) ) { type = BCF_HL_FMT; ss += 4; }
        else if ( !strncasecmp("FILTER/",ss,7) ) { type = BCF_HL_FLT; ss += 7; }
        else if ( !strncasecmp("^INFO/",ss,6) ) { type = BCF_HL_INFO; ss += 6; keep_info = 1; }
        else if ( !strncasecmp("^INF/",ss,5) ) { type = BCF_HL_INFO; ss += 5; keep_info = 1; }
        else if ( !strncasecmp("^FORMAT/",ss,8) ) { type = BCF_HL_FMT; ss += 8; keep_fmt = 1; }
        else if ( !strncasecmp("^FMT/",ss,5) ) { type = BCF_HL_FMT; ss += 5; keep_fmt = 1; }
        else if ( !strncasecmp("^FILTER/",ss,8) ) { type = BCF_HL_FLT; ss += 8; keep_flt = 1; }

        char *se = ss;
        while ( *se && *se!=',' ) se++;
        str.l = 0;
        kputsn(ss, se-ss, &str);

        if ( type==BCF_HL_FLT )
        {
            if ( !keep_flt )
            {
                args->flt_keep_pass = 1;
                tag->handler = remove_filter;
                tag->key = strdup(str.s);
                tag->hdr_id = bcf_hdr_id2int(args->hdr, BCF_DT_ID, tag->key);
                if ( !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FLT,tag->hdr_id) )
                {
                    if ( args->keep_sites )
                        error("Error: The filter \"%s\" is not defined in the header, cannot use the -k option\n", str.s);
                    else
                        fprintf(stderr,"Warning: The filter \"%s\" is not defined in the header\n", str.s);
                }
                else if ( !args->keep_sites ) bcf_hdr_remove(args->hdr_out,BCF_HL_FLT,tag->key);
            }
            else
            {
                int value, ret = khash_str2int_get(keep, str.s, &value);
                if ( ret==-1 ) khash_str2int_set(keep, strdup(str.s),1<<BCF_HL_FLT);
                else khash_str2int_set(keep, str.s, value | 1<<BCF_HL_FLT);
                args->nrm--;
            }
        }
        else if ( type!=BCF_HL_GEN )
        {
            int id = bcf_hdr_id2int(args->hdr,BCF_DT_ID,str.s);
            if ( !bcf_hdr_idinfo_exists(args->hdr,type,id) )
            {
                if ( args->keep_sites )
                    error("Error: The tag \"%s\" is not defined in the header, cannot use the -k option\n", str.s);
                else
                    fprintf(stderr,"Warning: The tag \"%s\" not defined in the header\n", str.s);

                tag->key = strdup(str.s);
                if ( type==BCF_HL_INFO )
                {
                    tag->handler = remove_info_tag;
                    if ( needs_info ) error("Error: `--remove INFO/%s` is executed first, cannot combine with `--set-id %s`\n",tag->key,args->set_ids_fmt);
                }
                else if ( type==BCF_HL_FMT ) tag->handler = remove_format_tag;
            }
            else if ( (type==BCF_HL_FMT && keep_fmt) || (type==BCF_HL_INFO && keep_info) )
            {
                int value, ret = khash_str2int_get(keep, str.s, &value);
                if ( ret==-1 ) khash_str2int_set(keep, strdup(str.s),1<<type);
                else khash_str2int_set(keep, str.s, value | 1<<type);
                args->nrm--;
            }
            else
            {
                tag->key = strdup(str.s);
                if ( type==BCF_HL_INFO )
                {
                    tag->handler = remove_info_tag;
                    if ( needs_info ) error("Error: `--remove INFO/%s` is executed first, cannot combine with `--set-id %s`\n",tag->key,args->set_ids_fmt);
                }
                else if ( type==BCF_HL_FMT ) tag->handler = remove_format_tag;
                if ( !args->keep_sites ) bcf_hdr_remove(args->hdr_out,type,tag->key);
            }
        }
        else if ( !strcasecmp("ID",str.s) ) tag->handler = remove_id;
        else if ( !strcasecmp("FILTER",str.s) )
        {
            tag->handler = remove_filter;
            if ( !args->keep_sites ) remove_hdr_lines(args->hdr_out,BCF_HL_FLT);
        }
        else if ( !strcasecmp("QUAL",str.s) ) tag->handler = remove_qual;
        else if ( !strcasecmp("INFO",str.s) )
        {
            if ( needs_info ) error("Error: `--remove INFO` is executed first, cannot combine with `--set-id %s`\n",args->set_ids_fmt);
            tag->handler = remove_info;
            if ( !args->keep_sites ) remove_hdr_lines(args->hdr_out,BCF_HL_INFO);
        }
        else if ( !strcasecmp("FMT",str.s) || !strcasecmp("FORMAT",str.s) )
        {
            tag->handler = remove_format;
            if ( !args->keep_sites ) remove_hdr_lines(args->hdr_out,BCF_HL_FMT);
        }
        else if ( str.l )
        {
            int id = bcf_hdr_id2int(args->hdr, BCF_DT_ID, str.s);
            if ( bcf_hdr_idinfo_exists(args->hdr,BCF_HL_INFO,id) ) error("Error: did you mean INFO/%s?\n",str.s);
            if ( bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FMT,id) ) error("Error: did you mean FORMAT/%s?\n",str.s);

            if ( !args->keep_sites )
            {
                if ( str.s[0]=='#' && str.s[1]=='#' )
                    bcf_hdr_remove(args->hdr_out,BCF_HL_GEN,str.s+2);
                else
                    bcf_hdr_remove(args->hdr_out,BCF_HL_STR,str.s);
            }
            args->nrm--;
        }

        ss = *se ? se+1 : se;
    }
    free(str.s);
    if ( keep_flt || keep_info || keep_fmt )
    {
        int j;
        for (j=0; j<args->hdr->nhrec; j++)
        {
            bcf_hrec_t *hrec = args->hdr->hrec[j];
            if ( hrec->type!=BCF_HL_FLT && hrec->type!=BCF_HL_INFO && hrec->type!=BCF_HL_FMT ) continue;
            if ( !keep_flt && hrec->type==BCF_HL_FLT ) continue;
            if ( !keep_info && hrec->type==BCF_HL_INFO ) continue;
            if ( !keep_fmt && hrec->type==BCF_HL_FMT ) continue;
            int k = bcf_hrec_find_key(hrec,"ID");
            assert( k>=0 ); // this should always be true for valid VCFs
            int value, ret = khash_str2int_get(keep,hrec->vals[k],&value);
            if ( ret==0 && value>>hrec->type ) // keep
            {
                if ( hrec->type==BCF_HL_FLT && !strcmp("PASS",hrec->vals[k]) ) args->flt_keep_pass = 1;
                continue;
            }
            args->nrm++;
            args->rm = (rm_tag_t*) realloc(args->rm,sizeof(rm_tag_t)*args->nrm);
            rm_tag_t *tag = &args->rm[args->nrm-1];
            if ( hrec->type==BCF_HL_INFO ) tag->handler = remove_info_tag;
            else if ( hrec->type==BCF_HL_FMT ) tag->handler = remove_format_tag;
            else
            {
                tag->handler = remove_filter;
                tag->hdr_id = bcf_hdr_id2int(args->hdr, BCF_DT_ID, hrec->vals[k]);
            }
            tag->key = strdup(hrec->vals[k]);
            if ( !args->keep_sites ) bcf_hdr_remove(args->hdr_out,hrec->type,tag->key);
        }
    }
    khash_str2int_destroy_free(keep);
    if ( !args->nrm ) error("No matching tag in -x %s\n", args->remove_annots);
    if (bcf_hdr_sync(args->hdr_out) < 0)
        error_errno("[%s] Failed to update header", __func__);
}
static void init_header_lines(args_t *args)
{
    if ( args->header_fname )
    {
        htsFile *file = hts_open(args->header_fname, "rb");
        if ( !file ) error("Error reading %s\n", args->header_fname);
        kstring_t str = {0,0,0};
        while ( hts_getline(file, KS_SEP_LINE, &str) > 0 )
        {
            if ( bcf_hdr_append(args->hdr_out,str.s) ) error("Could not parse %s: %s\n", args->header_fname, str.s);
            bcf_hdr_append(args->hdr,str.s);    // the input file may not have the header line if run with -h (and nothing else)
        }
        if ( hts_close(file)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->header_fname);
        free(str.s);
    }
    if ( args->header_lines )
    {
        int i, n = dbuf_n(args->header_lines);
        for (i=0; i<n; i++)
        {
            char *line = dbuf_ith(args->header_lines,i);
            if ( bcf_hdr_append(args->hdr_out,line) ) error("Could not parse the header line: %s\n", line);
            bcf_hdr_append(args->hdr,line);    // the input file may not have the header line if run with -H (and nothing else)
        }
        dbuf_destroy_free(args->header_lines);
        args->header_lines = NULL;
    }
    if (bcf_hdr_sync(args->hdr_out) < 0)
        error_errno("[%s] Failed to update output header", __func__);
    if (bcf_hdr_sync(args->hdr) < 0)
        error_errno("[%s] Failed to update input header", __func__);
}
static int vcf_getter_info_str2str(args_t *args, bcf1_t *rec, annot_col_t *col, void **ptr, int *mptr)
{
    return bcf_get_info_string(args->tgts_hdr,rec,col->hdr_key_src,ptr,mptr);
}
static int vcf_getter_id2str(args_t *args, bcf1_t *rec, annot_col_t *col, void **ptr, int *mptr)
{
    char *str = *((char**)ptr);
    int i, len = strlen(rec->d.id);
    if ( len >= *mptr ) str = realloc(str, len+1);
    for (i=0; i<len; i++)
        str[i] = rec->d.id[i]==';' ? ',' : rec->d.id[i];
    str[len] = 0;
    *((char**)ptr) = str;
    *mptr = len+1;
    return len;
}
inline static int vcf_getter_filter2str_core(bcf_hdr_t *hdr, bcf1_t *rec, char **ptr, int *mptr)
{
    if ( !(rec->unpacked & BCF_UN_FLT) ) bcf_unpack(rec, BCF_UN_FLT);

    kstring_t str;
    str.s = *ptr;
    str.m = *mptr;
    str.l = 0;

    int i;
    if ( rec->d.n_flt )
    {
        for (i=0; i<rec->d.n_flt; i++)
        {
            if (i) kputc(',', &str);
            kputs(bcf_hdr_int2id(hdr,BCF_DT_ID,rec->d.flt[i]), &str);
        }
    }
    else kputc('.', &str);

    *ptr  = str.s;
    *mptr = str.m;
    return str.l;
}
static int vcf_getter_filter2str_local(args_t *args, bcf1_t *rec, annot_col_t *col, void **ptr, int *mptr)
{
    return vcf_getter_filter2str_core(args->hdr_out, args->current_rec, (char**)ptr, mptr);
}
static int vcf_getter_filter2str(args_t *args, bcf1_t *rec, annot_col_t *col, void **ptr, int *mptr)
{
    return vcf_getter_filter2str_core(args->tgts_hdr, rec, (char**)ptr, mptr);
}
static int setter_filter(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    if ( !data ) error("Error: the --merge-logic option cannot be used with FILTER (yet?)\n");

    // note: so far this works only with one filter, not a list of filters
    annot_line_t *tab = (annot_line_t*) data;
    if ( tab->cols[col->icol][0]=='.' && !tab->cols[col->icol][1] ) // don't overwrite with a missing value unless asked
    {
        if ( (col->replace & CARRY_OVER_MISSING) && (col->replace & (REPLACE_ALL|REPLACE_NON_MISSING)) ) bcf_update_filter(args->hdr_out,line,NULL,0);
        return 0;
    }
    hts_expand(int,1,args->mtmpi,args->tmpi);
    args->tmpi[0] = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, tab->cols[col->icol]);
    if ( args->tmpi[0]<0 ) error("The FILTER \"%s\" is not defined in the header, was the -h option provided?\n", tab->cols[col->icol]);
    if ( col->replace & SET_OR_APPEND ) return bcf_add_filter(args->hdr_out,line,args->tmpi[0]);
    if ( !(col->replace & REPLACE_MISSING) )
    {
        bcf_update_filter(args->hdr_out,line,NULL,0);
        return bcf_update_filter(args->hdr_out,line,args->tmpi,1);
    }

    // only update missing FILTER
    if ( !(line->unpacked & BCF_UN_FLT) ) bcf_unpack(line, BCF_UN_FLT);
    if ( !line->d.n_flt )
        return bcf_update_filter(args->hdr_out,line,args->tmpi,1);

    return 0;
}
static int vcf_setter_filter(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    int i, ret = 0;
    bcf1_t *rec = (bcf1_t*) data;
    if ( !(rec->unpacked & BCF_UN_FLT) ) bcf_unpack(rec, BCF_UN_FLT);
    if ( !(line->unpacked & BCF_UN_FLT) ) bcf_unpack(line, BCF_UN_FLT);
    if ( !rec->d.n_flt ) // don't overwrite with a missing value unless asked
    {
        if ( (col->replace & CARRY_OVER_MISSING) && (col->replace & (REPLACE_ALL|REPLACE_NON_MISSING)) ) bcf_update_filter(args->hdr_out,line,NULL,0);
        return 0;
    }
    if ( col->replace & (SET_OR_APPEND|REPLACE_MISSING) )
    {
        if ( (col->replace & REPLACE_MISSING) && line->d.n_flt ) return 0; // only update missing FILTER
        for (i=0; i<rec->d.n_flt; i++)
        {
            const char *flt = bcf_hdr_int2id(args->files->readers[1].header, BCF_DT_ID, rec->d.flt[i]);
            if ( bcf_add_filter(args->hdr_out,line,bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, flt)) < 0 ) ret = -1;
        }
        return ret;
    }
    hts_expand(int,rec->d.n_flt,args->mtmpi,args->tmpi);
    for (i=0; i<rec->d.n_flt; i++)
    {
        const char *flt = bcf_hdr_int2id(args->files->readers[1].header, BCF_DT_ID, rec->d.flt[i]);
        args->tmpi[i] = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, flt);
    }
    bcf_update_filter(args->hdr_out,line,NULL,0);
    return bcf_update_filter(args->hdr_out,line,args->tmpi,rec->d.n_flt);
}
static int setter_pos(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    if ( tab->cols[col->icol] && tab->cols[col->icol][0]=='.' && !tab->cols[col->icol][1] ) return 0; // don't replace with "."
    char *tmp;
    int pos = strtol(tab->cols[col->icol], &tmp, 10);
    if ( tmp==tab->cols[col->icol] )
        error("Could not parse -POS at %s:%"PRId64" .. [%s]\n",bcf_seqname(args->hdr,line),(int64_t)line->pos+1,tab->cols[col->icol]);
    line->pos = pos - 1;
    return 0;
}
static int setter_id(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    if ( !data ) error("Error: the --merge-logic option cannot be used with ID (yet?)\n");
    if ( col->replace & MATCH_VALUE ) return 0;

    // possible cases:
    //      IN  ANNOT   OUT     ACHIEVED_BY
    //      x   y       x        -c +ID
    //      x   y       y        -c ID
    //      x   y       x,y      -c =ID
    //      x   .       x        -c +ID, ID
    //      x   .       .        -x ID
    //      .   y       y        -c +ID, -c ID
    //
    annot_line_t *tab = (annot_line_t*) data;
    if ( tab->cols[col->icol] && tab->cols[col->icol][0]=='.' && !tab->cols[col->icol][1] ) return 0; // don't replace with "."
    if ( col->replace & SET_OR_APPEND ) return bcf_add_id(args->hdr_out,line,tab->cols[col->icol]);
    if ( !(col->replace & REPLACE_MISSING) ) return bcf_update_id(args->hdr_out,line,tab->cols[col->icol]);

    // running with +ID, only update missing ids
    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) )
        return bcf_update_id(args->hdr_out,line,tab->cols[col->icol]);
    return 0;
}
static int vcf_setter_id(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    if ( col->replace & MATCH_VALUE ) return 0;

    bcf1_t *rec = (bcf1_t*) data;

    char *id;
    if ( col->getter )
    {
        int nret = col->getter(args,rec,col,&col->ptr,&col->mptr);
        id = (char*) col->ptr;
        if ( nret<=0 || (nret==1 && *id=='.') ) return 0;   // don't replace with "."
    }
    else
    {
        if ( rec->d.id && rec->d.id[0]=='.' && !rec->d.id[1] ) return 0;    // don't replace with "."
        id = rec->d.id;
    }
    if ( col->replace & SET_OR_APPEND ) return bcf_add_id(args->hdr_out,line,id);
    if ( !(col->replace & REPLACE_MISSING) ) return bcf_update_id(args->hdr_out,line,id);

    // running with +ID, only update missing ids
    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) )
        return bcf_update_id(args->hdr_out,line,id);
    return 0;
}
static int vcf_setter_ref(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    if ( !strcmp(rec->d.allele[0],line->d.allele[0]) ) return 0;    // no update necessary
    const char **als = (const char**) malloc(sizeof(char*)*line->n_allele);
    als[0] = rec->d.allele[0];
    int i;
    for (i=1; i<line->n_allele; i++) als[i] = line->d.allele[i];
    int ret = bcf_update_alleles(args->hdr_out, line, als, line->n_allele);
    free(als);
    return ret;
}
static int vcf_setter_alt(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int i;
    if ( line->n_allele>1 && (col->replace & REPLACE_MISSING) ) return 0;
    if ( rec->n_allele==line->n_allele )
    {
        for (i=1; i<rec->n_allele; i++) if ( strcmp(rec->d.allele[i],line->d.allele[i]) ) break;
        if ( i==rec->n_allele ) return 0;   // no update necessary
    }
    const char **als = (const char**) malloc(sizeof(char*)*rec->n_allele);
    als[0] = line->d.allele[0];
    for (i=1; i<rec->n_allele; i++) als[i] = rec->d.allele[i];
    int ret = bcf_update_alleles(args->hdr_out, line, als, rec->n_allele);
    free(als);
    return ret;
}
static int setter_qual(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    if ( !data ) error("Error: the --merge-logic option cannot be used with QUAL (yet?)\n");

    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol];
    if ( str[0]=='.' && str[1]==0 ) // don't overwrite with a missing value unless asked
    {
        if ( (col->replace & CARRY_OVER_MISSING) && (col->replace & (REPLACE_ALL|REPLACE_NON_MISSING)) ) bcf_float_set_missing(line->qual);
        return 0;
    }
    if ( (col->replace & REPLACE_MISSING) && !bcf_float_is_missing(line->qual) ) return 0;

    line->qual = strtod(str, &str);
    if ( str == tab->cols[col->icol] )
        error("Could not parse %s at %s:%"PRId64" .. [%s]\n", col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1,tab->cols[col->icol]);
    return 0;
}
static int vcf_setter_qual(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    if ( bcf_float_is_missing(rec->qual) )  // don't overwrite with a missing value unless asked
    {
        if ( (col->replace & CARRY_OVER_MISSING) && (col->replace & (REPLACE_ALL|REPLACE_NON_MISSING)) ) bcf_float_set_missing(line->qual);
        return 0;
    }
    if ( (col->replace & REPLACE_MISSING) && !bcf_float_is_missing(line->qual) ) return 0;
    line->qual = rec->qual;
    return 0;
}
static int setter_info_flag(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    if ( !data ) error("Error: the --merge-logic option cannot be used with INFO type=Flag (yet?)\n");

    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol];
    if ( str[0]=='.' && str[1]==0 ) // don't overwrite with a missing value unless asked
    {
        if ( (col->replace & CARRY_OVER_MISSING) && (col->replace & (REPLACE_ALL|REPLACE_NON_MISSING)) ) bcf_update_info_flag(args->hdr_out,line,col->hdr_key_dst,NULL,0);
        return 0;
    }

    if ( str[0]=='1' && str[1]==0 ) return bcf_update_info_flag(args->hdr_out,line,col->hdr_key_dst,NULL,1);
    if ( str[0]=='0' && str[1]==0 ) return bcf_update_info_flag(args->hdr_out,line,col->hdr_key_dst,NULL,0);
    error("Could not parse %s at %s:%"PRId64" .. [%s]\n", col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1,tab->cols[col->icol]);
    return -1;
}
static int vcf_setter_info_flag(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int flag = bcf_get_info_flag(args->files->readers[1].header,rec,col->hdr_key_src,NULL,NULL);
    bcf_update_info_flag(args->hdr_out,line,col->hdr_key_dst,NULL,flag);
    return 0;
}
static int setter_ARinfo_int32(args_t *args, bcf1_t *line, annot_col_t *col, int nals, char **als, int ntmpi)
{
    if ( col->number==BCF_VL_A && ntmpi!=nals-1 && (ntmpi!=1 || args->tmpi[0]!=bcf_int32_missing || args->tmpi[1]!=bcf_int32_vector_end) )
        error("Incorrect number of values (%d) for the %s tag at %s:%"PRId64"\n", ntmpi,col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
    else if ( col->number==BCF_VL_R && ntmpi!=nals && (ntmpi!=1 || args->tmpi[0]!=bcf_int32_missing || args->tmpi[1]!=bcf_int32_vector_end) )
        error("Incorrect number of values (%d) for the %s tag at %s:%"PRId64"\n", ntmpi,col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(args->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%"PRId64"\n", bcf_seqname(args->hdr,line),(int64_t) line->pos+1);

    // fill in any missing values in the target VCF (or all, if not present)
    int ntmpi2 = bcf_get_info_float(args->hdr, line, col->hdr_key_dst, &args->tmpi2, &args->mtmpi2);
    if ( ntmpi2 < ndst ) hts_expand(int32_t,ndst,args->mtmpi2,args->tmpi2);

    int i;
    for (i=0; i<ndst; i++)
    {
        if ( map[i]<0 )
        {
            if ( ntmpi2 < ndst ) args->tmpi2[i] = bcf_int32_missing;
            continue;
        }
        if ( ntmpi2==ndst && (col->replace & REPLACE_MISSING)
                && args->tmpi2[i]!=bcf_int32_missing
                && args->tmpi2[i]!=bcf_int32_vector_end ) continue;

        args->tmpi2[i] = args->tmpi[ map[i] ];
    }
    return bcf_update_info_int32(args->hdr_out,line,col->hdr_key_dst,args->tmpi2,ndst);
}
static int setter_info_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;

    // This is a bit hacky, only to reuse existing code with minimal changes:
    //      -c =TAG will now behave as -l TAG:APPEND for integers
    if ( col->replace & SET_OR_APPEND ) col->merge_method=MM_APPEND;

    if ( !tab )
    {
        if ( col->merge_method!=MM_SUM && col->merge_method!=MM_AVG &&
             col->merge_method!=MM_MIN && col->merge_method!=MM_MAX &&
             col->merge_method!=MM_APPEND &&
             col->merge_method!=MM_APPEND_MISSING )
            error("Error: at the moment only the sum,avg,min,max,append,append-missing options are supported with --merge-logic for INFO type=Integer\n");
    }

    int i,ntmpi = 0;
    if ( (col->replace & SET_OR_APPEND) && !col->mm_dbl_nused )
    {
        ntmpi = bcf_get_info_int32(args->hdr, line, col->hdr_key_dst, &args->tmpi, &args->mtmpi);
        if ( ntmpi>0 && (args->tmpi[0]!=bcf_int32_missing || (col->replace & CARRY_OVER_MISSING)) )
        {
            col->mm_dbl_nused = col->mm_dbl_ndat = ntmpi;
            hts_expand(double,col->mm_dbl_nused,col->mm_dbl_nalloc,col->mm_dbl);
            for (i=0; i<ntmpi; i++)
                col->mm_dbl[i] = args->tmpi[i];
            col->mm_dbl_ndat = 1;
        }
        ntmpi = 0;
    }
    if ( tab )  // has data, not flushing yet
    {
        char *str = tab->cols[col->icol], *end = str;
        if ( str[0]=='.' && str[1]==0 && col->merge_method!=MM_APPEND_MISSING && !(col->replace & CARRY_OVER_MISSING) ) return 1;

        while ( *end )
        {
            ntmpi++;
            hts_expand(int32_t,ntmpi,args->mtmpi,args->tmpi);
            if ( str[0]=='.' && (str[1]==0 || str[1]==',') )
            {
                if ( col->merge_method==MM_APPEND_MISSING || (col->replace & CARRY_OVER_MISSING) )
                    args->tmpi[ntmpi-1] = bcf_int32_missing;
                else
                    ntmpi--;
                if ( str[1]==0 ) end = str+1;
                str += 2;
            }
            else
            {
                args->tmpi[ntmpi-1] = strtol(str, &end, 10);
                if ( end==str )
                    error("Could not parse %s at %s:%"PRId64" .. [%s]\n", col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1,tab->cols[col->icol]);
                str = end+1;
            }
        }
        if ( col->merge_method!=MM_FIRST )
        {
            if ( !col->mm_dbl_nused )
            {
                col->mm_dbl_nused = ntmpi;
                hts_expand(double,col->mm_dbl_nused,col->mm_dbl_nalloc,col->mm_dbl);
                for (i=0; i<ntmpi; i++)
                    col->mm_dbl[i] = args->tmpi[i];
            }
            else
            {
                if ( col->merge_method==MM_APPEND || col->merge_method==MM_APPEND_MISSING )
                {
                    int nori = col->mm_dbl_nused;
                    col->mm_dbl_nused += ntmpi;
                    hts_expand(double,col->mm_dbl_nused,col->mm_dbl_nalloc,col->mm_dbl);
                    for (i=0; i<ntmpi; i++)
                        col->mm_dbl[i+nori] = args->tmpi[i];
                }
                else
                {
                    if ( ntmpi!=col->mm_dbl_nused ) error("Error: cannot merge fields of unequal length\n");
                    if ( col->merge_method==MM_SUM || col->merge_method==MM_AVG )
                        for (i=0; i<ntmpi; i++) col->mm_dbl[i] += args->tmpi[i];
                    else if ( col->merge_method==MM_MIN )
                        for (i=0; i<ntmpi; i++) { if ( col->mm_dbl[i] > args->tmpi[i] ) col->mm_dbl[i] = args->tmpi[i]; }
                    else if ( col->merge_method==MM_MAX )
                        for (i=0; i<ntmpi; i++) { if ( col->mm_dbl[i] < args->tmpi[i] ) col->mm_dbl[i] = args->tmpi[i]; }
                }
            }
            col->mm_dbl_ndat++;
            return 1;
        }
    }
    else if ( col->merge_method==MM_SUM || col->merge_method==MM_MIN || col->merge_method==MM_MAX || col->merge_method==MM_APPEND || col->merge_method==MM_APPEND_MISSING )
    {
        ntmpi = col->mm_dbl_nused;
        hts_expand(int32_t,ntmpi,args->mtmpi,args->tmpi);
        for (i=0; i<ntmpi; i++) args->tmpi[i] = col->mm_dbl[i];
        col->mm_dbl_nused = col->mm_dbl_ndat = 0;
    }
    else if ( col->merge_method==MM_AVG )
    {
        ntmpi = col->mm_dbl_nused;
        hts_expand(int32_t,ntmpi,args->mtmpi,args->tmpi);
        for (i=0; i<ntmpi; i++) args->tmpi[i] = col->mm_dbl[i]/col->mm_dbl_ndat;
        col->mm_dbl_nused = col->mm_dbl_ndat = 0;
    }

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R )
        return setter_ARinfo_int32(args,line,col,tab->nals,tab->als,ntmpi);

    if ( col->replace & REPLACE_MISSING )
    {
        int ret = bcf_get_info_int32(args->hdr, line, col->hdr_key_dst, &args->tmpi2, &args->mtmpi2);
        if ( ret>0 && args->tmpi2[0]!=bcf_int32_missing ) return 0;
    }
    return bcf_update_info_int32(args->hdr_out,line,col->hdr_key_dst,args->tmpi,ntmpi);
}
static int vcf_setter_info_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int ntmpi = bcf_get_info_int32(args->files->readers[1].header,rec,col->hdr_key_src,&args->tmpi,&args->mtmpi);
    if ( ntmpi < 0 ) return 0;    // nothing to add

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R )
        return setter_ARinfo_int32(args,line,col,rec->n_allele,rec->d.allele,ntmpi);

    if ( col->replace & REPLACE_MISSING )
    {
        int ret = bcf_get_info_int32(args->hdr, line, col->hdr_key_dst, &args->tmpi2, &args->mtmpi2);
        if ( ret>0 && args->tmpi2[0]!=bcf_int32_missing ) return 0;
    }

    return bcf_update_info_int32(args->hdr_out,line,col->hdr_key_dst,args->tmpi,ntmpi);
}
static int setter_ARinfo_real(args_t *args, bcf1_t *line, annot_col_t *col, int nals, char **als, int ntmpf)
{
    if ( col->number==BCF_VL_A && ntmpf!=nals-1 && (ntmpf!=1 || !bcf_float_is_missing(args->tmpf[0]) || !bcf_float_is_vector_end(args->tmpf[0])) )
        error("Incorrect number of values (%d) for the %s tag at %s:%"PRId64"\n", ntmpf,col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
    else if ( col->number==BCF_VL_R && ntmpf!=nals && (ntmpf!=1 || !bcf_float_is_missing(args->tmpf[0]) || !bcf_float_is_vector_end(args->tmpf[0])) )
        error("Incorrect number of values (%d) for the %s tag at %s:%"PRId64"\n", ntmpf,col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(args->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%"PRId64"\n", bcf_seqname(args->hdr,line),(int64_t) line->pos+1);

    // fill in any missing values in the target VCF (or all, if not present)
    int ntmpf2 = bcf_get_info_float(args->hdr, line, col->hdr_key_dst, &args->tmpf2, &args->mtmpf2);
    if ( ntmpf2 < ndst ) hts_expand(float,ndst,args->mtmpf2,args->tmpf2);

    int i;
    for (i=0; i<ndst; i++)
    {
        if ( map[i]<0 )
        {
            if ( ntmpf2 < ndst ) bcf_float_set_missing(args->tmpf2[i]);
            continue;
        }
        if ( ntmpf2==ndst && (col->replace & REPLACE_MISSING)
                && !bcf_float_is_missing(args->tmpf2[i])
                && !bcf_float_is_vector_end(args->tmpf2[i]) ) continue;

        args->tmpf2[i] = args->tmpf[ map[i] ];
    }
    return bcf_update_info_float(args->hdr_out,line,col->hdr_key_dst,args->tmpf2,ndst);
}
static int setter_info_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;

    // This is a bit hacky, only to reuse existing code with minimal changes:
    //      -c =TAG will now behave as -l TAG:APPEND for floats
    if ( col->replace & SET_OR_APPEND ) col->merge_method=MM_APPEND;

    if ( !tab )
    {
        if ( col->merge_method!=MM_SUM && col->merge_method!=MM_AVG &&
             col->merge_method!=MM_MIN && col->merge_method!=MM_MAX &&
             col->merge_method!=MM_APPEND &&
             col->merge_method!=MM_APPEND_MISSING )
            error("Error: at the moment only the sum,avg,min,max,append,append-missing options are supported with --merge-logic for INFO type=Float\n");
    }

    int i,ntmpf = 0;
    if ( (col->replace & SET_OR_APPEND) && !col->mm_dbl_nused )
    {
        ntmpf = bcf_get_info_float(args->hdr, line, col->hdr_key_dst, &args->tmpf, &args->mtmpf);
        if ( ntmpf>0 && (!bcf_float_is_missing(args->tmpf[0]) || (col->replace & CARRY_OVER_MISSING)) )
        {
            col->mm_dbl_nused = ntmpf;
            hts_expand(double,col->mm_dbl_nused,col->mm_dbl_nalloc,col->mm_dbl);
            for (i=0; i<ntmpf; i++)
                if ( bcf_float_is_missing(args->tmpf[i]) )
                    bcf_double_set_missing(col->mm_dbl[i]);
                else
                    col->mm_dbl[i] = args->tmpf[i];
            col->mm_dbl_ndat = 1;
        }
        ntmpf = 0;
    }
    if ( tab )  // data row, not just flushing
    {
        char *str = tab->cols[col->icol], *end = str;
        if ( str[0]=='.' && str[1]==0 && col->merge_method!=MM_APPEND_MISSING && !(col->replace & CARRY_OVER_MISSING) ) return 1;

        while ( *end )
        {
            ntmpf++;
            hts_expand(float,ntmpf,args->mtmpf,args->tmpf);
            if ( str[0]=='.' && (str[1]==0 || str[1]==',') )
            {
                if ( col->merge_method==MM_APPEND_MISSING || (col->replace & CARRY_OVER_MISSING) )
                    bcf_float_set_missing(args->tmpf[ntmpf-1]);
                else
                    ntmpf--;
                if ( str[1]==0 ) end = str+1;
                str += 2;
            }
            else
            {
                args->tmpf[ntmpf-1] = strtod(str, &end);
                if ( end==str )
                    error("Could not parse %s at %s:%"PRId64" .. [%s]\n", col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1,tab->cols[col->icol]);
                str = end+1;
            }
        }
        if ( col->merge_method!=MM_FIRST )
        {
            if ( !col->mm_dbl_nused )
            {
                col->mm_dbl_nused = ntmpf;
                hts_expand(double,col->mm_dbl_nused,col->mm_dbl_nalloc,col->mm_dbl);
                for (i=0; i<ntmpf; i++)
                {
                    if ( bcf_float_is_missing(args->tmpf[i]) )
                        bcf_double_set_missing(col->mm_dbl[i]);
                    else
                        col->mm_dbl[i] = args->tmpf[i];
                }
            }
            else
            {
                if ( col->merge_method==MM_APPEND || col->merge_method==MM_APPEND_MISSING )
                {
                    int nori = col->mm_dbl_nused;
                    col->mm_dbl_nused += ntmpf;
                    hts_expand(double,col->mm_dbl_nused,col->mm_dbl_nalloc,col->mm_dbl);
                    for (i=0; i<ntmpf; i++)
                    {
                        if ( bcf_float_is_missing(args->tmpf[i]) )
                            bcf_double_set_missing(col->mm_dbl[i+nori]);
                        else
                            col->mm_dbl[i+nori] = args->tmpf[i];
                    }
                }
                else
                {
                    if ( ntmpf!=col->mm_dbl_nused ) error("Error: cannot merge fields of unequal length\n");
                    if ( col->merge_method==MM_SUM || col->merge_method==MM_AVG )
                        for (i=0; i<ntmpf; i++) col->mm_dbl[i] += args->tmpf[i];
                    else if ( col->merge_method==MM_MIN )
                        for (i=0; i<ntmpf; i++) { if ( col->mm_dbl[i] > args->tmpf[i] ) col->mm_dbl[i] = args->tmpf[i]; }
                    else if ( col->merge_method==MM_MAX )
                        for (i=0; i<ntmpf; i++) { if ( col->mm_dbl[i] < args->tmpf[i] ) col->mm_dbl[i] = args->tmpf[i]; }
                }
            }
            col->mm_dbl_ndat++;
            return 1;
        }
    }
    else if ( col->merge_method==MM_SUM || col->merge_method==MM_MIN || col->merge_method==MM_MAX || col->merge_method==MM_APPEND || col->merge_method==MM_APPEND_MISSING )
    {
        ntmpf = col->mm_dbl_nused;
        hts_expand(int32_t,ntmpf,args->mtmpf,args->tmpf);
        for (i=0; i<ntmpf; i++)
        {
            if ( bcf_double_is_missing(col->mm_dbl[i]) )
                bcf_float_set_missing(args->tmpf[i]);
            else
                args->tmpf[i] = col->mm_dbl[i];
        }
        col->mm_dbl_nused = col->mm_dbl_ndat = 0;
    }
    else if ( col->merge_method==MM_AVG )
    {
        ntmpf = col->mm_dbl_nused;
        hts_expand(int32_t,ntmpf,args->mtmpf,args->tmpf);
        for (i=0; i<ntmpf; i++) args->tmpf[i] = col->mm_dbl[i]/col->mm_dbl_ndat;
        col->mm_dbl_nused = col->mm_dbl_ndat = 0;
    }

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R )
        return setter_ARinfo_real(args,line,col,tab->nals,tab->als,ntmpf);

    if ( col->replace & REPLACE_MISSING )
    {
        int ret = bcf_get_info_float(args->hdr, line, col->hdr_key_dst, &args->tmpf2, &args->mtmpf2);
        if ( ret>0 && !bcf_float_is_missing(args->tmpf2[0]) ) return 0;
    }

    return bcf_update_info_float(args->hdr_out,line,col->hdr_key_dst,args->tmpf,ntmpf);
}
static int vcf_setter_info_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int ntmpf = bcf_get_info_float(args->files->readers[1].header,rec,col->hdr_key_src,&args->tmpf,&args->mtmpf);
    if ( ntmpf < 0 ) return 0;    // nothing to add

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R )
        return setter_ARinfo_real(args,line,col,rec->n_allele,rec->d.allele,ntmpf);

    if ( col->replace & REPLACE_MISSING )
    {
        int ret = bcf_get_info_float(args->hdr, line, col->hdr_key_dst, &args->tmpf2, &args->mtmpf2);
        if ( ret>0 && !bcf_float_is_missing(args->tmpf2[0]) ) return 0;
    }

    return bcf_update_info_float(args->hdr_out,line,col->hdr_key_dst,args->tmpf,ntmpf);
}
int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst); // see vcfmerge.c
static int setter_ARinfo_string(args_t *args, bcf1_t *line, annot_col_t *col, int nals, char **als)
{
    assert( col->merge_method==MM_FIRST );

    int nsrc = 1, lsrc = 0;
    while ( args->tmps[lsrc] )
    {
        if ( args->tmps[lsrc]==',' ) nsrc++;
        lsrc++;
    }
    if ( col->number==BCF_VL_A && nsrc!=nals-1 && (nsrc!=1 || args->tmps[0]!='.' || args->tmps[1]!=0 ) )
        error("Incorrect number of values (%d) for the %s tag at %s:%"PRId64"\n", nsrc,col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
    else if ( col->number==BCF_VL_R && nsrc!=nals && (nsrc!=1 || args->tmps[0]!='.' || args->tmps[1]!=0 ) )
        error("Incorrect number of values (%d) for the %s tag at %s:%"PRId64"\n", nsrc,col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(args->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%"PRId64"\n", bcf_seqname(args->hdr,line),(int64_t) line->pos+1);

    // fill in any missing values in the target VCF (or all, if not present)
    int i, empty = 0, nstr, mstr = args->tmpks.m;
    nstr = bcf_get_info_string(args->hdr, line, col->hdr_key_dst, &args->tmpks.s, &mstr);
    args->tmpks.m = mstr;
    if ( nstr<0 || (nstr==1 && args->tmpks.s[0]=='.' && args->tmpks.s[1]==0) )
    {
        empty = 0;
        args->tmpks.l = 0;
        kputc('.',&args->tmpks);
        for (i=1; i<ndst; i++) kputs(",.",&args->tmpks);
    }
    else args->tmpks.l = nstr;
    for (i=0; i<ndst; i++)
    {
        if ( map[i]<0 )
        {
            if ( empty ) copy_string_field(".",0,1,&args->tmpks,i);
            continue;
        }
        if ( col->replace & REPLACE_MISSING )
        {
            // Do not replace filled values. The field must be looked up again because
            // of realloc in copy_string_field
            int n = 0;
            char *str = args->tmpks.s;
            while ( *str && n<i )
            {
                if ( *str==',' ) n++;
                str++;
            }
            if ( str[0]!='.' || (str[1]!=',' && str[1]!=0) ) continue;  // value already set
        }
        int ret = copy_string_field(args->tmps,map[i],lsrc,&args->tmpks,i);
        if ( ret!=0 ) error("[%s:%d %s] Failed to copy a string field\n",  __FILE__,__LINE__,__func__);
    }
    return bcf_update_info_string(args->hdr_out,line,col->hdr_key_dst,args->tmpks.s);
}
void khash_str2int_clear_free(void *_hash)
{
    khash_t(str2int) *hash = (khash_t(str2int)*)_hash;
    khint_t k;
    if (hash == 0) return;
    for (k = 0; k < kh_end(hash); ++k)
        if (kh_exist(hash, k)) free((char*)kh_key(hash, k));
    kh_clear(str2int, hash);
}
static const char *escape_string(const char *str, char needle[], char **rmme, size_t *len)
{
    kstring_t tmp = {0,0,0};
    const char *bp = str, *ep = str;
    while ( *ep )
    {
        int i = 0;
        while ( needle[i] && needle[i]!=*ep ) i++;
        if ( !needle[i] ) { ep++; continue; }
        kputsn(bp,ep-bp,&tmp);
        ksprintf(&tmp,"%%%X",*ep);
        bp = ++ep;
    }
    if ( !tmp.l )
    {
        *len = strlen(str);
        return str;
    }
    kputs(bp,&tmp);
    *len  = tmp.l;
    *rmme = tmp.s;
    return tmp.s;
}
static int setter_info_str(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    if ( (col->replace & REPLACE_MISSING) && col->number!=BCF_VL_A && col->number!=BCF_VL_R )
    {
        int ret = bcf_get_info_string(args->hdr, line, col->hdr_key_dst, &args->tmps2, &args->mtmps2);
        if ( ret>0 && (args->tmps2[0]!='.' || args->tmps2[1]!=0) ) return 0;
    }

    // This is a bit hacky, only to reuse existing code with minimal changes:
    //      -c =TAG will now behave as -l TAG:unique for strings
    if ( col->replace & SET_OR_APPEND ) col->merge_method=MM_UNIQUE;

    annot_line_t *tab = (annot_line_t*) data;
    const char *escaped = NULL;
    char *rmme = NULL;

    size_t len = 0;
    if ( tab )
    {
        char *str = tab->cols[col->icol];
        if ( !str || !*str ) return 0;
        if ( !str[1] && str[0]=='.' && col->merge_method!=MM_APPEND_MISSING && !(col->replace & CARRY_OVER_MISSING) ) return 1;
        char needle[] = {';','=',0};
        escaped = escape_string(tab->cols[col->icol],needle,&rmme,&len);
    }

    if ( col->merge_method!=MM_FIRST )
    {
        if ( col->number==BCF_VL_A || col->number==BCF_VL_R )
            error("Error: the --merge-logic option cannot be used with INFO tags Type=String,Number={A,R,G}\n");

        if ( data )
        {
            assert( col->merge_method==MM_APPEND || col->merge_method==MM_APPEND_MISSING || col->merge_method==MM_UNIQUE );
            if ( col->merge_method==MM_UNIQUE )
            {
                if ( !col->mm_str_hash ) col->mm_str_hash = (khash_t(str2int)*)khash_str2int_init();
                if ( khash_str2int_has_key(col->mm_str_hash, escaped) )
                {
                    free(rmme);
                    return 1;
                }
                khash_str2int_inc(col->mm_str_hash, strdup(escaped));
            }

            if ( (col->replace & SET_OR_APPEND) && !col->mm_kstr.l )
            {
                int m = col->mm_kstr.m;
                int n = bcf_get_info_string(args->hdr, line, col->hdr_key_dst, &col->mm_kstr.s, &m);
                col->mm_kstr.m = m;
                if ( n>0 && ((col->replace & CARRY_OVER_MISSING) || col->mm_kstr.s[0]!='.' || col->mm_kstr.s[1]) ) col->mm_kstr.l = n;
            }

            if ( col->mm_kstr.l ) kputc(',',&col->mm_kstr);
            kputs(escaped, &col->mm_kstr);
            free(rmme);
            return 1;
        }
        if ( col->mm_kstr.l )
        {
            hts_expand(char,col->mm_kstr.l+1,args->mtmps,args->tmps);
            memcpy(args->tmps,col->mm_kstr.s,col->mm_kstr.l+1);
        }
        else
        {
            free(rmme);
            return 0;
        }

        // flush the line
        if ( col->merge_method==MM_UNIQUE )
            khash_str2int_clear_free(col->mm_str_hash);
        col->mm_kstr.l = 0;
    }
    else
    {
        assert(tab);
        hts_expand(char,len+1,args->mtmps,args->tmps);
        memcpy(args->tmps,escaped,len+1);
        if ( col->number==BCF_VL_A || col->number==BCF_VL_R )
            return setter_ARinfo_string(args,line,col,tab->nals,tab->als);
    }
    int ret = bcf_update_info_string(args->hdr_out,line,col->hdr_key_dst,args->tmps);
    free(rmme);
    return ret;
}
static int vcf_setter_info_str(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;

    if ( col->getter )
        col->getter(args,rec,col,(void**)&args->tmps, &args->mtmps);
    else
    {
        int ntmps = bcf_get_info_string(args->files->readers[1].header,rec,col->hdr_key_src,&args->tmps,&args->mtmps);
        if ( ntmps < 0 ) return 0;    // nothing to add
    }

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R )
        return setter_ARinfo_string(args,line,col,rec->n_allele,rec->d.allele);

    if ( col->replace & REPLACE_MISSING )
    {
        int ret = bcf_get_info_string(args->hdr, line, col->hdr_key_dst, &args->tmps2, &args->mtmps2);
        if ( ret>0 && (args->tmps2[0]!='.' || args->tmps2[1]!=0) ) return 0;
    }

    return bcf_update_info_string(args->hdr_out,line,col->hdr_key_dst,args->tmps);
}
static int genotypes_to_string(args_t *args, int nsrc1, int32_t *src, int nsmpl_dst, kstring_t *str)
{
    int i, isrc, idst;
    int blen = nsrc1 > 1 ? nsrc1 + 1 : 1;   // typically the genotypes take three bytes 0/1, no 0-termination is needed

gt_length_too_big:
    str->l = 0;
    for (idst=0; idst<nsmpl_dst; idst++)
    {
        isrc = args->sample_map ? args->sample_map[idst] : idst;
        if ( isrc==-1 )
        {
            kputc_('.', str);
            for (i=1; i < blen; i++) kputc_(0, str);
            continue;
        }

        size_t plen = str->l;
        int32_t *ptr = src + isrc*nsrc1;
        for (i=0; i<nsrc1 && ptr[i]!=bcf_int32_vector_end; i++)
        {
            if ( i ) kputc("/|"[bcf_gt_is_phased(ptr[i])], str);
            if ( bcf_gt_is_missing(ptr[i]) ) kputc('.', str);
            else kputw(bcf_gt_allele(ptr[i]), str);
        }
        if ( i==0 ) kputc('.', str);
        if ( str->l - plen > blen )
        {
            // too many alternate alleles or ploidy is too large, the genotype does not fit
            // three characters ("0/0" vs "10/10").
            blen *= 2;
            goto gt_length_too_big;
        }
        plen = str->l - plen;
        while ( plen < blen )
        {
            kputc_(0, str);
            plen++;
        }
    }
    return 0;
}
static int vcf_setter_format_gt(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int nsrc = bcf_get_genotypes(args->files->readers[1].header,rec,&args->tmpi,&args->mtmpi);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    // Genotypes are internally represented as integers. This is a complication when
    // adding as a different Type=String field, such as FMT/newGT:=GT
    if ( strcmp(col->hdr_key_src,col->hdr_key_dst) )
    {
        int nsmpl_dst = bcf_hdr_nsamples(args->hdr_out);
        int nsmpl_src = bcf_hdr_nsamples(args->files->readers[1].header);
        genotypes_to_string(args,nsrc/nsmpl_src,args->tmpi,nsmpl_dst,&args->tmpks);
        return bcf_update_format_char(args->hdr_out,line,col->hdr_key_dst,args->tmpks.s,args->tmpks.l);
    }

    if ( !args->sample_map )
        return bcf_update_genotypes(args->hdr_out,line,args->tmpi,nsrc);

    int i, j, ndst = bcf_get_genotypes(args->hdr,line,&args->tmpi2,&args->mtmpi2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(args->hdr_out);
    nsrc /= bcf_hdr_nsamples(args->files->readers[1].header);
    if ( ndst<=0 )  // field not present in dst file
    {
        if ( col->replace & REPLACE_NON_MISSING ) return 0;
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(args->hdr_out), args->mtmpi2, args->tmpi2);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            int32_t *dst = args->tmpi2 + nsrc*i;
            if ( args->sample_map[i]==-1 )
            {
                dst[0] = bcf_gt_missing;
                for (j=1; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = args->tmpi + nsrc*args->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_genotypes(args->hdr_out,line,args->tmpi2,nsrc*bcf_hdr_nsamples(args->hdr_out));
    }
    else if ( ndst >= nsrc )
    {
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            if ( args->sample_map[i]==-1 ) continue;
            int32_t *src = args->tmpi  + nsrc*args->sample_map[i];
            int32_t *dst = args->tmpi2 + ndst*i;
            if ( (col->replace & REPLACE_NON_MISSING) && bcf_gt_is_missing(dst[0]) ) continue;
            if ( (col->replace & REPLACE_MISSING)  && !bcf_gt_is_missing(dst[0]) ) continue;
            for (j=0; j<nsrc; j++) dst[j] = src[j];
            for (; j<ndst; j++) dst[j] = bcf_int32_vector_end;
        }
        return bcf_update_genotypes(args->hdr_out,line,args->tmpi2,ndst*bcf_hdr_nsamples(args->hdr_out));
    }
    else    // ndst < nsrc
    {
        hts_expand(int32_t, nsrc*bcf_hdr_nsamples(args->hdr_out), args->mtmpi3, args->tmpi3);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            int32_t *ori = args->tmpi2 + ndst*i;
            int32_t *dst = args->tmpi3 + nsrc*i;
            int keep_ori = 0;
            if ( args->sample_map[i]==-1 ) keep_ori = 1;
            else if ( (col->replace & REPLACE_NON_MISSING) && bcf_gt_is_missing(ori[0]) ) keep_ori = 1;
            else if ( (col->replace & REPLACE_MISSING)  && !bcf_gt_is_missing(ori[0]) ) keep_ori = 1;
            if ( keep_ori )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nsrc; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = args->tmpi + nsrc*args->sample_map[i];
                for (j=0; j<nsrc; j++) dst[j] = src[j];
            }
        }
        return bcf_update_genotypes(args->hdr_out,line,args->tmpi3,nsrc*bcf_hdr_nsamples(args->hdr_out));
    }
}
static int count_vals(annot_line_t *tab, int icol_beg, int icol_end)
{
    int i, nmax = 1;
    for (i=icol_beg; i<icol_end; i++)
    {
        char *str = tab->cols[i], *end = str;
        if ( str[0]=='.' && !str[1] )
        {
            // missing value
            if ( !nmax ) nmax = 1;
            continue;
        }
        int n = 1;
        while ( *end )
        {
            if ( *end==',' ) n++;
            end++;
        }
        if ( nmax<n ) nmax = n;
    }
    return nmax;
}
static int core_setter_format_int(args_t *args, bcf1_t *line, annot_col_t *col, int32_t *vals, int nvals)
{
    if ( !args->sample_map )
        return bcf_update_format_int32(args->hdr_out,line,col->hdr_key_dst,vals,nvals*args->nsmpl_annot);

    int i, j, ndst = bcf_get_format_int32(args->hdr,line,col->hdr_key_dst,&args->tmpi2,&args->mtmpi2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(args->hdr_out);
    if ( ndst<=0 )
    {
        if ( col->replace & REPLACE_NON_MISSING ) return 0;    // overwrite only if present
        hts_expand(int32_t, nvals*bcf_hdr_nsamples(args->hdr_out), args->mtmpi2, args->tmpi2);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            int32_t *dst = args->tmpi2 + nvals*i;
            if ( args->sample_map[i]==-1 )
            {
                dst[0] = bcf_int32_missing;
                for (j=1; j<nvals; j++) dst[j] = bcf_int32_vector_end;
            }
            else
            {
                int32_t *src = vals + nvals*args->sample_map[i];
                for (j=0; j<nvals; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_int32(args->hdr_out,line,col->hdr_key_dst,args->tmpi2,nvals*bcf_hdr_nsamples(args->hdr_out));
    }
    else if ( ndst >= nvals )
    {
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            if ( args->sample_map[i]==-1 ) continue;
            int32_t *src = vals  + nvals*args->sample_map[i];
            int32_t *dst = args->tmpi2 + ndst*i;
            // possible cases:
            //      in annot out
            //       x  y     x     TAG,-TAG,=TAG    .. REPLACE_ALL, REPLACE_NON_MISSING, SET_OR_APPEND
            //       x  y     y    +TAG              .. REPLACE_MISSING
            //       .  y     .    =TAG              .. SET_OR_APPEND
            //       .  y     y     TAG,+TAG,-TAG    .. REPLACE_ALL, REPLACE_MISSING, REPLACE_NON_MISSING
            //       x  .     x     TAG,+TAG         .. REPLACE_ALL, REPLACE_MISSING
            //       x  .     .    -TAG              .. REPLACE_NON_MISSING
            if ( col->replace & REPLACE_NON_MISSING ) { if ( dst[0]==bcf_int32_missing ) continue; }
            else if ( col->replace & REPLACE_MISSING ) { if ( dst[0]!=bcf_int32_missing ) continue; }
            else if ( col->replace & REPLACE_ALL ) { if ( src[0]==bcf_int32_missing ) continue; }
            for (j=0; j<nvals; j++) dst[j] = src[j];
            for (; j<ndst; j++) dst[j] = bcf_int32_vector_end;
        }
        return bcf_update_format_int32(args->hdr_out,line,col->hdr_key_dst,args->tmpi2,ndst*bcf_hdr_nsamples(args->hdr_out));
    }
    else    // ndst < nvals
    {
        hts_expand(int32_t, nvals*bcf_hdr_nsamples(args->hdr_out), args->mtmpi3, args->tmpi3);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            int32_t *ann = vals + nvals*args->sample_map[i];
            int32_t *ori = args->tmpi2 + ndst*i;                // ori vcf line
            int32_t *dst = args->tmpi3 + nvals*i;               // expanded buffer
            int use_new_ann = 1;
            if ( args->sample_map[i]==-1 ) use_new_ann = 0;
            else if ( col->replace & REPLACE_NON_MISSING ) { if ( ori[0]==bcf_int32_missing ) use_new_ann = 0; }
            else if ( col->replace & REPLACE_MISSING ) { if ( ori[0]!=bcf_int32_missing ) use_new_ann = 0; }
            else if ( col->replace & REPLACE_ALL ) { if ( ann[0]==bcf_int32_missing ) use_new_ann = 0; }
            if ( !use_new_ann )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nvals; j++) dst[j] = bcf_int32_vector_end;
            }
            else
                for (j=0; j<nvals; j++) dst[j] = ann[j];
        }
        return bcf_update_format_int32(args->hdr_out,line,col->hdr_key_dst,args->tmpi3,nvals*bcf_hdr_nsamples(args->hdr_out));
    }
}
static int core_setter_format_real(args_t *args, bcf1_t *line, annot_col_t *col, float *vals, int nvals)
{
    if ( !args->sample_map )
        return bcf_update_format_float(args->hdr_out,line,col->hdr_key_dst,vals,nvals*args->nsmpl_annot);

    int i, j, ndst = bcf_get_format_float(args->hdr,line,col->hdr_key_dst,&args->tmpf2,&args->mtmpf2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(args->hdr_out);
    if ( ndst<=0 )
    {
        if ( col->replace & REPLACE_NON_MISSING ) return 0;    // overwrite only if present
        hts_expand(float, nvals*bcf_hdr_nsamples(args->hdr_out), args->mtmpf2, args->tmpf2);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            float *dst = args->tmpf2 + nvals*i;
            if ( args->sample_map[i]==-1 )
            {
                bcf_float_set_missing(dst[0]);
                for (j=1; j<nvals; j++) bcf_float_set_vector_end(dst[j]);
            }
            else
            {
                float *src = vals + nvals*args->sample_map[i];
                for (j=0; j<nvals; j++) dst[j] = src[j];
            }
        }
        return bcf_update_format_float(args->hdr_out,line,col->hdr_key_dst,args->tmpf2,nvals*bcf_hdr_nsamples(args->hdr_out));
    }
    else if ( ndst >= nvals )
    {
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            if ( args->sample_map[i]==-1 ) continue;
            float *src = vals  + nvals*args->sample_map[i];
            float *dst = args->tmpf2 + ndst*i;
            if ( col->replace & REPLACE_NON_MISSING ) { if ( bcf_float_is_missing(dst[0]) ) continue; }
            else if ( col->replace & REPLACE_MISSING ) { if ( !bcf_float_is_missing(dst[0]) ) continue; }
            else if ( col->replace & REPLACE_ALL ) { if ( bcf_float_is_missing(src[0]) ) continue; }
            for (j=0; j<nvals; j++) dst[j] = src[j];
            for (; j<ndst; j++) bcf_float_set_vector_end(dst[j]);
        }
        return bcf_update_format_float(args->hdr_out,line,col->hdr_key_dst,args->tmpf2,ndst*bcf_hdr_nsamples(args->hdr_out));
    }
    else    // ndst < nvals
    {
        hts_expand(float, nvals*bcf_hdr_nsamples(args->hdr_out), args->mtmpf3, args->tmpf3);
        for (i=0; i<bcf_hdr_nsamples(args->hdr_out); i++)
        {
            float *ann = vals + nvals*args->sample_map[i];
            float *ori = args->tmpf2 + ndst*i;                // ori vcf line
            float *dst = args->tmpf3 + nvals*i;               // expanded buffer
            int use_new_ann = 1;
            if ( args->sample_map[i]==-1 ) use_new_ann = 0;
            else if ( col->replace & REPLACE_NON_MISSING ) { if ( bcf_float_is_missing(ori[0]) ) use_new_ann = 0; }
            else if ( col->replace & REPLACE_MISSING ) { if ( !bcf_float_is_missing(ori[0]) ) use_new_ann = 0; }
            else if ( col->replace & REPLACE_ALL ) { if ( bcf_float_is_missing(ann[0]) ) use_new_ann = 0; }
            if ( !use_new_ann )
            {
                for (j=0; j<ndst; j++) dst[j] = ori[j];
                for (; j<nvals; j++) bcf_float_set_vector_end(dst[j]);
            }
            else
                for (j=0; j<nvals; j++) dst[j] = ann[j];
        }
        return bcf_update_format_float(args->hdr_out,line,col->hdr_key_dst,args->tmpf3,nvals*bcf_hdr_nsamples(args->hdr_out));
    }
}
static int core_setter_format_str(args_t *args, bcf1_t *line, annot_col_t *col, char **vals)
{
    if ( !args->sample_map )
        return bcf_update_format_string(args->hdr_out,line,col->hdr_key_dst,(const char**)vals,args->nsmpl_annot);

    int i;
    args->tmpp2[0] = args->tmps2;
    int ret = bcf_get_format_string(args->hdr,line,col->hdr_key_dst,&args->tmpp2,&args->mtmps2);
    args->tmps2 = args->tmpp2[0];   // tmps2 might be realloced

    int nsmpl = bcf_hdr_nsamples(args->hdr_out);
    if ( ret<=0 )   // not present in dst
    {
        hts_expand(char,bcf_hdr_nsamples(args->hdr_out)*2,args->mtmps2,args->tmps2);
        char *tmp = args->tmps2;
        for (i=0; i<nsmpl; i++)
        {
            tmp[0] = '.';
            tmp[1] = 0;
            args->tmpp2[i] = tmp;
            tmp += 2;
        }
    }
    for (i=0; i<nsmpl; i++)
    {
        if ( args->sample_map[i]==-1 ) continue;
        char **src = vals + args->sample_map[i];
        char **dst = args->tmpp2 + i;

        if ( col->replace & REPLACE_NON_MISSING ) { if ( (*dst)[0]=='.' && (*dst)[1]==0 ) continue; }
        else if ( col->replace & REPLACE_MISSING ) { if ( (*dst)[0]!='.' || (*dst)[1]!=0 ) continue; }
        else if ( col->replace & REPLACE_ALL ) { if ( (*src)[0]=='.' && (*src)[1]==0 ) continue; }
        *dst = *src;
    }
    return bcf_update_format_string(args->hdr_out,line,col->hdr_key_dst,(const char**)args->tmpp2,nsmpl);
}
static int setter_format_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    if ( !data ) error("Error: the --merge-logic option cannot be used with FORMAT tags (yet?)\n");

    annot_line_t *tab = (annot_line_t*) data;
    if ( col->icol+args->nsmpl_annot > tab->ncols )
        error("Incorrect number of values for %s at %s:%"PRId64"\n",col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
    int nvals = count_vals(tab,col->icol,col->icol+args->nsmpl_annot);
    hts_expand(int32_t,nvals*args->nsmpl_annot,args->mtmpi,args->tmpi);

    int icol = col->icol, ismpl;
    for (ismpl=0; ismpl<args->nsmpl_annot; ismpl++)
    {
        int32_t *ptr = args->tmpi + ismpl*nvals;
        int ival = 0;

        char *str = tab->cols[icol];
        while ( *str )
        {
            if ( str[0]=='.' && (!str[1] || str[1]==',') )  // missing value
            {
                ptr[ival++] = bcf_int32_missing;
                str += str[1] ? 2 : 1;
                continue;
            }

            char *end = str;
            ptr[ival] = strtol(str, &end, 10);
            if ( end==str )
                error("Could not parse %s at %s:%"PRId64" .. [%s]\n", col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1,tab->cols[col->icol]);

            ival++;
            str = *end ? end+1 : end;
        }
        while ( ival<nvals ) ptr[ival++] = bcf_int32_vector_end;
        icol++;
    }
    return core_setter_format_int(args,line,col,args->tmpi,nvals);
}
static int setter_format_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    if ( !data ) error("Error: the --merge-logic option cannot be used with FORMAT tags (yet?)\n");

    annot_line_t *tab = (annot_line_t*) data;
    if ( col->icol+args->nsmpl_annot > tab->ncols )
        error("Incorrect number of values for %s at %s:%"PRId64"\n",col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
    int nvals = count_vals(tab,col->icol,col->icol+args->nsmpl_annot);
    hts_expand(float,nvals*args->nsmpl_annot,args->mtmpf,args->tmpf);

    int icol = col->icol, ismpl;
    for (ismpl=0; ismpl<args->nsmpl_annot; ismpl++)
    {
        float *ptr = args->tmpf + ismpl*nvals;
        int ival = 0;

        char *str = tab->cols[icol];
        while ( *str )
        {
            if ( str[0]=='.' && (!str[1] || str[1]==',') )  // missing value
            {
                bcf_float_set_missing(ptr[ival]);
                ival++;
                str += str[1] ? 2 : 1;
                continue;
            }

            char *end = str;
            ptr[ival] = strtod(str, &end);
            if ( end==str )
                error("Could not parse %s at %s:%"PRId64" .. [%s]\n", col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1,tab->cols[col->icol]);

            ival++;
            str = *end ? end+1 : end;
        }
        while ( ival<nvals ) { bcf_float_set_vector_end(ptr[ival]); ival++; }
        icol++;
    }
    return core_setter_format_real(args,line,col,args->tmpf,nvals);
}
static int setter_format_str(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    if ( !data ) error("Error: the --merge-logic option cannot be used with FORMAT tags (yet?)\n");

    annot_line_t *tab = (annot_line_t*) data;
    if ( col->icol+args->nsmpl_annot > tab->ncols )
        error("Incorrect number of values for %s at %s:%"PRId64"\n",col->hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);

    char needle[] = {':',0};
    int ismpl;
    for (ismpl=0; ismpl<args->nsmpl_annot; ismpl++)
    {
        size_t len;
        char *rmme = NULL;
        const char *str = escape_string(tab->cols[col->icol + ismpl],needle,&rmme,&len);
        args->tmpp[ismpl] = rmme ? rmme : strdup(str);
    }
    int ret = core_setter_format_str(args,line,col,args->tmpp);
    for (ismpl=0; ismpl<args->nsmpl_annot; ismpl++) free(args->tmpp[ismpl]);
    return ret;
}
static int determine_ploidy(int nals, int *vals, int nvals1, uint8_t *smpl, int nsmpl)
{
    int i, j, ndip = nals*(nals+1)/2, max_ploidy = 0;
    for (i=0; i<nsmpl; i++)
    {
        int *ptr = vals + i*nvals1;
        int has_value = 0;
        for (j=0; j<nvals1; j++)
        {
            if ( ptr[j]==bcf_int32_vector_end ) break;
            if ( ptr[j]!=bcf_int32_missing ) has_value = 1;
        }
        if ( has_value )
        {
            if ( j==ndip )
            {
                smpl[i] = 2;
                max_ploidy = 2;
            }
            else if ( j==nals )
            {
                smpl[i] = 1;
                if ( !max_ploidy ) max_ploidy = 1;
            }
            else return -j;
        }
        else smpl[i] = 0;
    }
    return max_ploidy;
}
static int vcf_setter_format_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int nsrc = bcf_get_format_int32(args->files->readers[1].header,rec,col->hdr_key_src,&args->tmpi,&args->mtmpi);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error
    int nsmpl_src = bcf_hdr_nsamples(args->files->readers[1].header);
    int nsrc1 = nsrc / nsmpl_src;
    if ( col->number!=BCF_VL_G && col->number!=BCF_VL_R && col->number!=BCF_VL_A )
        return core_setter_format_int(args,line,col,args->tmpi,nsrc1);

    // create mapping from src to dst genotypes, haploid and diploid version
    int nmap_hap = col->number==BCF_VL_G || col->number==BCF_VL_R ? rec->n_allele : rec->n_allele - 1;
    int *map_hap = vcmp_map_ARvalues(args->vcmp,nmap_hap,line->n_allele,line->d.allele,rec->n_allele,rec->d.allele);
    if ( !map_hap ) error("REF alleles not compatible at %s:%"PRId64"\n", bcf_seqname(args->hdr,line),(int64_t) line->pos+1);

    int i, j;
    if ( rec->n_allele==line->n_allele )
    {
        // alleles unchanged?
        for (i=0; i<rec->n_allele; i++) if ( map_hap[i]!=i ) break;
        if ( i==rec->n_allele )
            return core_setter_format_int(args,line,col,args->tmpi,nsrc1);
    }

    int nsmpl_dst = rec->n_sample;
    int ndst  = bcf_get_format_int32(args->hdr,line,col->hdr_key_dst,&args->tmpi2,&args->mtmpi2);
    int ndst1 = ndst / nsmpl_dst;
    if ( ndst <= 0 )
    {
        if ( col->replace & REPLACE_NON_MISSING ) return 0;  // overwrite only if present
        if ( col->number==BCF_VL_G )
            ndst1 = line->n_allele*(line->n_allele+1)/2;
        else
            ndst1 = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
        hts_expand(int, ndst1*nsmpl_dst, args->mtmpi2, args->tmpi2);
        for (i=0; i<nsmpl_dst; i++)
        {
            int32_t *dst = args->tmpi2 + i*ndst1;
            for (j=0; j<ndst1; j++) dst[j] = bcf_int32_missing;
        }
    }

    int nmap_dip = 0, *map_dip = NULL;
    if ( col->number==BCF_VL_G )
    {
        map_dip = vcmp_map_dipGvalues(args->vcmp, &nmap_dip);
        if ( !args->src_smpl_pld )
        {
            args->src_smpl_pld = (uint8_t*) malloc(nsmpl_src);
            args->dst_smpl_pld = (uint8_t*) malloc(nsmpl_dst);
        }
        int pld_src = determine_ploidy(rec->n_allele, args->tmpi, nsrc1, args->src_smpl_pld, nsmpl_src);
        if ( pld_src<0 )
            error("Unexpected number of %s values (%d) for %d alleles at %s:%"PRId64"\n", col->hdr_key_src,-pld_src, rec->n_allele, bcf_seqname(bcf_sr_get_header(args->files,1),rec),(int64_t) rec->pos+1);
        int pld_dst = determine_ploidy(line->n_allele, args->tmpi2, ndst1, args->dst_smpl_pld, nsmpl_dst);
        if ( pld_dst<0 )
            error("Unexpected number of %s values (%d) for %d alleles at %s:%"PRId64"\n", col->hdr_key_src,-pld_dst, line->n_allele, bcf_seqname(args->hdr,line),(int64_t) line->pos+1);

        int ndst1_new = pld_dst==1 ? line->n_allele : line->n_allele*(line->n_allele+1)/2;
        if ( ndst1_new != ndst1 )
        {
            if ( ndst1 ) error("todo: %s ndst1!=ndst .. %d %d  at %s:%"PRId64"\n",col->hdr_key_src,ndst1_new,ndst1,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
            ndst1 = ndst1_new;
            hts_expand(int32_t, ndst1*nsmpl_dst, args->mtmpi2, args->tmpi2);
        }
    }
    else if ( !ndst1 )
    {
        ndst1 = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
        hts_expand(int32_t, ndst1*nsmpl_dst, args->mtmpi2, args->tmpi2);
    }

    for (i=0; i<nsmpl_dst; i++)
    {
        int ii = args->sample_map ? args->sample_map[i] : i;
        int32_t *ptr_src = args->tmpi + i*nsrc1;
        int32_t *ptr_dst = args->tmpi2 + ii*ndst1;

        if ( col->number==BCF_VL_G )
        {
            if ( args->src_smpl_pld[ii] > 0 && args->dst_smpl_pld[i] > 0 && args->src_smpl_pld[ii]!=args->dst_smpl_pld[i] )
                error("Sample ploidy differs at %s:%"PRId64"\n", bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
            if ( !args->dst_smpl_pld[i] )
                for (j=0; j<ndst1; j++) ptr_dst[j] = bcf_int32_missing;
        }
        if ( col->number!=BCF_VL_G || args->src_smpl_pld[i]==1 )
        {
            for (j=0; j<nmap_hap; j++)
            {
                int k = map_hap[j];
                if ( k>=0 ) ptr_dst[k] = ptr_src[j];
            }
            if ( col->number==BCF_VL_G )
                for (j=line->n_allele; j<ndst1; j++) ptr_dst[j++] = bcf_int32_vector_end;
        }
        else
        {
            for (j=0; j<nmap_dip; j++)
            {
                int k = map_dip[j];
                if ( k>=0 ) ptr_dst[k] = ptr_src[j];
            }
        }
    }
    return bcf_update_format_int32(args->hdr_out,line,col->hdr_key_dst,args->tmpi2,nsmpl_dst*ndst1);
}
static int vcf_setter_format_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int nsrc = bcf_get_format_float(args->files->readers[1].header,rec,col->hdr_key_src,&args->tmpf,&args->mtmpf);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error
    int nsmpl_src = bcf_hdr_nsamples(args->files->readers[1].header);
    int nsrc1 = nsrc / nsmpl_src;
    if ( col->number!=BCF_VL_G && col->number!=BCF_VL_R && col->number!=BCF_VL_A )
        return core_setter_format_real(args,line,col,args->tmpf,nsrc1);

    // create mapping from src to dst genotypes, haploid and diploid version
    int nmap_hap = col->number==BCF_VL_G || col->number==BCF_VL_R ? rec->n_allele : rec->n_allele - 1;
    int *map_hap = vcmp_map_ARvalues(args->vcmp,nmap_hap,line->n_allele,line->d.allele,rec->n_allele,rec->d.allele);
    if ( !map_hap ) error("REF alleles not compatible at %s:%"PRId64"\n", bcf_seqname(args->hdr,line),(int64_t) line->pos+1);

    int i, j;
    if ( rec->n_allele==line->n_allele )
    {
        // alleles unchanged?
        for (i=0; i<rec->n_allele; i++) if ( map_hap[i]!=i ) break;
        if ( i==rec->n_allele )
            return core_setter_format_real(args,line,col,args->tmpf,nsrc1);
    }

    int nsmpl_dst = rec->n_sample;
    int ndst  = bcf_get_format_float(args->hdr,line,col->hdr_key_dst,&args->tmpf2,&args->mtmpf2);
    int ndst1 = ndst / nsmpl_dst;
    if ( ndst <= 0 )
    {
        if ( col->replace & REPLACE_NON_MISSING ) return 0;  // overwrite only if present
        if ( col->number==BCF_VL_G )
            ndst1 = line->n_allele*(line->n_allele+1)/2;
        else
            ndst1 = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
        hts_expand(float, ndst1*nsmpl_dst, args->mtmpf2, args->tmpf2);
        for (i=0; i<nsmpl_dst; i++)
        {
            float *dst = args->tmpf2 + i*ndst1;
            for (j=0; j<ndst1; j++) bcf_float_set_missing(dst[j]);
        }
    }

    int nmap_dip = 0, *map_dip = NULL;
    if ( col->number==BCF_VL_G )
    {
        map_dip = vcmp_map_dipGvalues(args->vcmp, &nmap_dip);
        if ( !args->src_smpl_pld )
        {
            args->src_smpl_pld = (uint8_t*) malloc(nsmpl_src);
            args->dst_smpl_pld = (uint8_t*) malloc(nsmpl_dst);
        }
        int pld_src = determine_ploidy(rec->n_allele, args->tmpi, nsrc1, args->src_smpl_pld, nsmpl_src);
        if ( pld_src<0 )
            error("Unexpected number of %s values (%d) for %d alleles at %s:%"PRId64"\n", col->hdr_key_src,-pld_src, rec->n_allele, bcf_seqname(bcf_sr_get_header(args->files,1),rec),(int64_t) rec->pos+1);
        int pld_dst = determine_ploidy(line->n_allele, args->tmpi2, ndst1, args->dst_smpl_pld, nsmpl_dst);
        if ( pld_dst<0 )
            error("Unexpected number of %s values (%d) for %d alleles at %s:%"PRId64"\n", col->hdr_key_src,-pld_dst, line->n_allele, bcf_seqname(args->hdr,line),(int64_t) line->pos+1);

        int ndst1_new = pld_dst==1 ? line->n_allele : line->n_allele*(line->n_allele+1)/2;
        if ( ndst1_new != ndst1 )
        {
            if ( ndst1 ) error("todo: %s ndst1!=ndst .. %d %d  at %s:%"PRId64"\n",col->hdr_key_src,ndst1_new,ndst1,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
            ndst1 = ndst1_new;
            hts_expand(float, ndst1*nsmpl_dst, args->mtmpf2, args->tmpf2);
        }
    }
    else if ( !ndst1 )
    {
        ndst1 = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
        hts_expand(float, ndst1*nsmpl_dst, args->mtmpf2, args->tmpf2);
    }

    for (i=0; i<nsmpl_dst; i++)
    {
        int ii = args->sample_map ? args->sample_map[i] : i;
        float *ptr_src = args->tmpf + i*nsrc1;
        float *ptr_dst = args->tmpf2 + ii*ndst1;

        if ( col->number==BCF_VL_G )
        {
            if ( args->src_smpl_pld[ii] > 0 && args->dst_smpl_pld[i] > 0 && args->src_smpl_pld[ii]!=args->dst_smpl_pld[i] )
                error("Sample ploidy differs at %s:%"PRId64"\n", bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
            if ( !args->dst_smpl_pld[i] )
                for (j=0; j<ndst1; j++) bcf_float_set_missing(ptr_dst[j]);
        }
        if ( col->number!=BCF_VL_G || args->src_smpl_pld[i]==1 )
        {
            for (j=0; j<nmap_hap; j++)
            {
                int k = map_hap[j];
                if ( k>=0 )
                {
                    if ( bcf_float_is_missing(ptr_src[j]) ) bcf_float_set_missing(ptr_dst[k]);
                    else if ( bcf_float_is_vector_end(ptr_src[j]) ) bcf_float_set_vector_end(ptr_dst[k]);
                    else ptr_dst[k] = ptr_src[j];
                }
            }
            if ( col->number==BCF_VL_G )
                for (j=line->n_allele; j<ndst1; j++) bcf_float_set_vector_end(ptr_dst[j]);
        }
        else
        {
            for (j=0; j<nmap_dip; j++)
            {
                int k = map_dip[j];
                if ( k>=0 )
                {
                    if ( bcf_float_is_missing(ptr_src[j]) ) bcf_float_set_missing(ptr_dst[k]);
                    else if ( bcf_float_is_vector_end(ptr_src[j]) ) bcf_float_set_vector_end(ptr_dst[k]);
                    else ptr_dst[k] = ptr_src[j];
                }
            }
        }
    }
    return bcf_update_format_float(args->hdr_out,line,col->hdr_key_dst,args->tmpf2,nsmpl_dst*ndst1);
}

static int vcf_setter_format_str(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    args->tmpp[0] = args->tmps;
    int ret = bcf_get_format_string(args->files->readers[1].header,rec,col->hdr_key_src,&args->tmpp,&args->mtmps);
    args->tmps = args->tmpp[0]; // tmps might be realloced
    if ( ret==-3 ) return 0;    // the tag is not present
    if ( ret<=0 ) return 1;     // error
    if ( strcmp("GT",col->hdr_key_dst) )
        return core_setter_format_str(args,line,col,args->tmpp);

    // Genotypes are internally represented as integers. This is a complication for FMT/GT:=oldGT
    // First determine the maximum number of alleles per-sample ndst1
    int nsmpl_src = bcf_hdr_nsamples(args->files->readers[1].header);
    int nsmpl_dst = bcf_hdr_nsamples(args->hdr_out);
    int isrc,idst, ndst1 = 0, nsrc1 = ret / nsmpl_src;
    char *ptr = args->tmps, *ptr_end = ptr + ret;
    while ( ptr < ptr_end )
    {
        char *smpl_end = ptr + nsrc1;
        int n = 1;
        while ( ptr < smpl_end )
        {
            if ( *ptr=='/' || *ptr=='|' ) n++;
            ptr++;
        }
        if ( ndst1 < n ) ndst1 = n;
    }
    assert( ndst1 );

    int ndst = ndst1*nsmpl_dst;
    hts_expand(int32_t,ndst,args->mtmpi,args->tmpi);
    hts_expand(char,ret+1,args->mtmps,args->tmps); args->tmps[ret] = 0; // the FORMAT string may not be 0-terminated
    for (idst=0; idst<nsmpl_dst; idst++)
    {
        int i = 0, is_phased = 0;
        int32_t *dst = args->tmpi + idst*ndst1;
        isrc = args->sample_map ? args->sample_map[idst] : idst;
        if ( isrc==-1 )
        {
            dst[0] = bcf_gt_missing;
            for (i=1; i<ndst1; i++) dst[i] = bcf_int32_vector_end;
            continue;
        }
        char *beg = args->tmps + isrc*nsrc1, *tmp;
        char *keep_ptr = beg+nsrc1, keep = *keep_ptr; *keep_ptr = 0;
        while ( *beg )
        {
            char *end = beg;
            while ( *end && *end!='/' && *end!='|' ) end++;
            if ( *beg=='.' && end-beg==1 ) dst[i] = bcf_gt_missing;
            else
            {
                if ( *end=='|' ) is_phased = 1;
                dst[i] = strtol(beg, &tmp, 10);
                if ( tmp!=end )
                    error("Could not parse the %s field at %s:%"PRId64" in %s\n", col->hdr_key_src,bcf_seqname(args->files->readers[1].header,rec),(int64_t) rec->pos+1,args->targets_fname);
                if ( dst[i] >= line->n_allele )
                    error("The source allele index is bigger than the number of destination alleles at %s:%"PRId64"\n", bcf_seqname(args->files->readers[1].header,rec),(int64_t) rec->pos+1);
                dst[i] = is_phased ? bcf_gt_phased(dst[i]) : bcf_gt_unphased(dst[i]);
            }
            beg = *end ? end+1 : end;
            i++;
        }
        *keep_ptr = keep;
        for (; i<ndst1; i++) dst[i] = bcf_int32_vector_end;
    }
    return bcf_update_genotypes(args->hdr_out,line,args->tmpi,ndst);
}
static int init_sample_map(args_t *args, bcf_hdr_t *src, bcf_hdr_t *dst)
{
    int i;
    if ( !args->sample_names )
    {
        args->nsmpl_annot = bcf_hdr_nsamples(dst);

        // tab annotation file, expecting that all samples are present: sample map not needed
        if ( !src ) return 0;

        int nmatch = 0;
        for (i=0; i<bcf_hdr_nsamples(src); i++)
        {
            int id = bcf_hdr_id2int(dst, BCF_DT_SAMPLE, src->samples[i]);
            if ( id!=-1 ) nmatch++;
        }
        if ( !nmatch ) return -1;   // No matching samples found in the source and the destination file

        args->nsample_map = bcf_hdr_nsamples(dst);
        args->sample_map  = (int*) malloc(sizeof(int)*args->nsample_map);
        for (i=0; i<args->nsample_map; i++)
        {
            int id = bcf_hdr_id2int(src, BCF_DT_SAMPLE, dst->samples[i]);
            args->sample_map[i] = id;   // idst -> isrc, -1 if not present
        }
        return 1;
    }

    args->nsample_map = bcf_hdr_nsamples(dst);
    args->sample_map  = (int*) malloc(sizeof(int)*args->nsample_map);
    for (i=0; i<args->nsample_map; i++) args->sample_map[i] = -1;

    int flags = !src ? SMPL_STRICT|SMPL_SINGLE|SMPL_REORDER : SMPL_STRICT|SMPL_SINGLE|SMPL_PAIR2;   // is tab vs vcf annotation file
    smpl_ilist_t *ilist = smpl_ilist_init(dst, args->sample_names, args->sample_is_file, flags);    // gives mapping dst->src
    if ( !ilist || !ilist->n ) error("Could not parse the samples: %s\n", args->sample_names);
    args->nsmpl_annot = ilist->n;
    int need_sample_map = args->nsmpl_annot==bcf_hdr_nsamples(dst) ? 0 : 1;
    for (i=0; i<args->nsmpl_annot; i++)
    {
        int idst = ilist->idx[i];
        const char *src_name = ilist->pair && ilist->pair[i] ? ilist->pair[i] : bcf_hdr_int2id(dst, BCF_DT_SAMPLE, idst);
        int isrc = i;
        if ( src )     // the annotation file is a VCF, not a tab-delimited file
        {
            isrc = bcf_hdr_id2int(src, BCF_DT_SAMPLE, src_name);
            if ( isrc==-1 ) error("Sample \"%s\" not found in the annotation file\n", src_name);
        }
        if ( isrc!=idst ) need_sample_map = 1;
        args->sample_map[idst] = isrc;
    }
    smpl_ilist_destroy(ilist);
    return need_sample_map;
}
static char *columns_complement(char *columns, void **skip_info, void **skip_fmt)
{
    kstring_t str = {0,0,0};
    char *ss = columns, *se = ss;
    while ( *ss )
    {
        if ( *se && *se!=',' ) { se++; continue; }
        if ( *ss!='^' )
        {
            if ( str.l ) kputc(',',&str);
            kputsn(ss, se-ss, &str);
            if ( !*se ) break;
            ss = ++se;
            continue;
        }

        if ( !strncasecmp("^INFO/",ss,6) )
        {
            if ( !*skip_info )
            {
                *skip_info = khash_str2int_init();
                if ( str.l ) kputc(',',&str);
                kputs("INFO",&str);
            }
            char tmp = *se; *se = 0;
            khash_str2int_inc(*skip_info, strdup(ss+6));
            *se = tmp;
        }
        else if ( !strncasecmp("^FORMAT/",ss,8) || !strncasecmp("^FMT/",ss,5) )
        {
            int n = !strncasecmp("^FMT/",ss,5) ? 5 : 8;
            if ( !*skip_fmt )
            {
                *skip_fmt = khash_str2int_init();
                if ( str.l ) kputc(',',&str);
                kputs("FORMAT",&str);
            }
            char tmp = *se; *se = 0;
            khash_str2int_inc(*skip_fmt, strdup(ss+n));
            *se = tmp;
        }
        else
        {
            if ( !*skip_info )
            {
                *skip_info = khash_str2int_init();
                if ( str.l ) kputc(',',&str);
                kputs("INFO",&str);
            }
            char tmp = *se; *se = 0;
            khash_str2int_inc(*skip_info, strdup(ss+1));
            *se = tmp;
        }

        if ( !*se ) break;
        ss = ++se;
    }
    free(columns);
    return str.s;
}
static void bcf_hrec_format_rename(bcf_hrec_t *hrec, char *tag, kstring_t *str)
{
    int j, nout = 0;
    ksprintf(str, "##%s=<", hrec->key);
    for (j=0; j<hrec->nkeys; j++)
    {
        if ( !strcmp("IDX",hrec->keys[j]) ) continue;
        if ( nout ) kputc(',',str);
        if ( !strcmp("ID", hrec->keys[j]) )
            ksprintf(str,"%s=%s", hrec->keys[j], tag);
        else
            ksprintf(str,"%s=%s", hrec->keys[j], hrec->vals[j]);
        nout++;
    }
    ksprintf(str,">\n");
}
static char *set_replace_mode(char *ss, int *replace)
{
    int mode = 0;
    while (*ss)
    {
        if ( *ss=='+' ) mode |= REPLACE_MISSING;
        else if ( *ss=='-' ) mode |= REPLACE_NON_MISSING;
        else if ( *ss=='=' ) mode |= SET_OR_APPEND;
        else if ( *ss=='.' ) mode |= CARRY_OVER_MISSING;
        else break;
        ss++;
    }
    if ( !mode ) mode = REPLACE_ALL;
// is exactly one bit set?
//    if ( mode && !(mode && ((mode & mode-1) == 0)) )
    *replace = mode;
    return ss;
}
static void rename_annots_push(args_t *args, char *src, char *dst);
static void init_columns(args_t *args)
{
    int need_sample_map = 0;
    int sample_map_ok = init_sample_map(args, args->tgts_is_vcf?args->files->readers[1].header:NULL, args->hdr);

    kstring_t tmp = {0,0,0};
    if ( args->columns_is_file )
    {
        int i,n;
        char **str = hts_readlist(args->columns, args->columns_is_file, &n);
        if ( !str ) error("Could not parse %s\n", args->columns);
        for (i=0; i<n; i++)
        {
            char *ptr = str[i];
            while ( *ptr && !isspace(*ptr) ) ptr++;
            if ( *ptr )
            {
                *ptr = 0;
                ptr++;
                while ( *ptr && isspace(*ptr) ) ptr++;
                if ( *ptr )
                {
                    if ( args->merge_method_str.l ) kputc(',',&args->merge_method_str);
                    kputs(str[i],&args->merge_method_str);
                    kputc(':',&args->merge_method_str);
                    kputs(ptr,&args->merge_method_str);
                }
            }
            if ( tmp.l ) kputc(',',&tmp);
            kputs(str[i],&tmp);
            free(str[i]);
        }
        free(str);
        free(args->columns);
        args->columns = tmp.s;
        tmp.l = tmp.m = 0;
        tmp.s = NULL;
    }

    void *skip_fmt = NULL, *skip_info = NULL;
    if ( args->tgts_is_vcf )
        args->columns = columns_complement(args->columns, &skip_info, &skip_fmt);

    kstring_t str = {0,0,0};
    char *ss = args->columns, *se = ss;
    args->ncols = 0;
    int icol = -1, has_fmt_str = 0;
    while ( *ss )
    {
        char *ptr;
        if ( *se && *se!=',' ) { se++; continue; }
        int replace;
        ss = set_replace_mode(ss, &replace);
        icol++;
        str.l = 0;
        kputsn(ss, se-ss, &str);
        if ( !str.s[0] || !strcasecmp("-",str.s) ) ;
        else if ( !strcasecmp("CHROM",str.s) ) args->chr_idx = icol;
        else if ( !strcasecmp("POS",str.s) )
        {
            if ( replace==REPLACE_NON_MISSING && !args->tgts_is_vcf )
            {
                args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
                annot_col_t *col = &args->cols[args->ncols-1];
                memset(col,0,sizeof(*col));
                col->icol = icol;
                col->replace = replace;
                col->setter  = setter_pos;
                col->hdr_key_src = strdup(str.s);
                col->hdr_key_dst = strdup(str.s);
                args->match_end = icol;
            }
            else
                args->beg_idx = icol;
        }
        else if ( !strcasecmp("FROM",str.s) || !strcasecmp("BEG",str.s) ) args->beg_idx = icol;
        else if ( !strcasecmp("TO",str.s) || !strcasecmp("END",str.s) ) args->end_idx = icol;
        else if ( !strcasecmp("REF",str.s) )
        {
            if ( args->tgts_is_vcf )
            {
                args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
                annot_col_t *col = &args->cols[args->ncols-1];
                memset(col,0,sizeof(*col));
                col->setter = vcf_setter_ref;
                col->hdr_key_src = strdup(str.s);
                col->hdr_key_dst = strdup(str.s);
            }
            else args->ref_idx = icol;
        }
        else if ( !strcasecmp("ALT",str.s) )
        {
            if ( args->tgts_is_vcf )
            {
                args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
                annot_col_t *col = &args->cols[args->ncols-1];
                memset(col,0,sizeof(*col));
                col->setter = vcf_setter_alt;
                col->hdr_key_src = strdup(str.s);
                col->hdr_key_dst = strdup(str.s);
                col->replace = replace;
                if ( args->pair_logic==-1 ) bcf_sr_set_opt(args->files,BCF_SR_PAIR_LOGIC,BCF_SR_PAIR_BOTH_REF);
            }
            else args->alt_idx = icol;
        }
        else if ( !strcasecmp("ID",str.s) || !strcasecmp("~ID",str.s) )
        {
            if ( replace & REPLACE_NON_MISSING ) error("Apologies, the -ID feature has not been implemented yet.\n");
            if ( str.s[0]=='~' ) replace = MATCH_VALUE;
            if ( args->tgts_is_vcf && (replace & MATCH_VALUE) ) error("todo: -c ~ID with -a VCF?\n");
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            memset(col,0,sizeof(*col));
            col->icol = icol;
            col->replace = replace;
            col->setter = args->tgts_is_vcf ? vcf_setter_id : setter_id;
            col->hdr_key_src = strdup(str.s);
            col->hdr_key_dst = strdup(str.s);
            if ( replace & MATCH_VALUE ) args->match_id = icol;
        }
        else if ( !strcasecmp("~INFO/END",str.s) && !args->tgts_is_vcf )
        {
            replace = MATCH_VALUE;
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            memset(col,0,sizeof(*col));
            col->icol = icol;
            col->replace = replace;
            col->setter  = NULL;
            col->hdr_key_src = strdup(str.s);
            col->hdr_key_dst = strdup(str.s);
            args->match_end = icol;
        }
        else if ( !strcasecmp("~POS",str.s) )
        {
            error("Error: the use of ~POS has been deprecated, use -POS to transfer the column POS.\n");
        }
        else if ( str.s[0]=='~' )
        {
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            memset(col,0,sizeof(*col));
            col->icol = icol;
            col->replace = MATCH_VALUE;
            col->setter  = NULL;
            col->hdr_key_src = strdup(str.s+1);
        }
        else if ( !strcasecmp("-POS",str.s) && !args->tgts_is_vcf )
        {
            if ( args->tgts_is_vcf ) error("Error: cannot use -POS, position can be replaced only from a tab-delimited file\n");
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            memset(col,0,sizeof(*col));
            col->icol = icol;
            col->replace = replace;
            col->setter  = setter_pos;
            col->hdr_key_src = strdup(str.s);
            col->hdr_key_dst = strdup(str.s);
            args->match_end = icol;
        }
        else if ( !strncasecmp("ID:=",str.s,4) )    // transfer a tag from INFO to ID column
        {
            if ( !args->tgts_is_vcf ) error("The annotation source must be a VCF for \"%s\"\n",str.s);
            if ( replace & REPLACE_NON_MISSING ) error("Apologies, the -ID feature has not been implemented yet.\n");
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            memset(col,0,sizeof(*col));
            col->icol = icol;
            col->replace = replace;
            col->setter = vcf_setter_id;
            col->getter = vcf_getter_info_str2str;
            str.s[2] = 0;
            col->hdr_key_dst = strdup(str.s);
            col->hdr_key_src = strncasecmp("INFO/",str.s+4,5) ? strdup(str.s+4) : strdup(str.s+4+5);
            int hdr_id = bcf_hdr_id2int(args->tgts_hdr, BCF_DT_ID,col->hdr_key_src);
            if ( !bcf_hdr_idinfo_exists(args->tgts_hdr,BCF_HL_INFO,hdr_id) )
                error("The INFO tag \"%s\" is not defined in %s\n", col->hdr_key_src, args->targets_fname);
            if ( bcf_hdr_id2type(args->tgts_hdr,BCF_HL_INFO,hdr_id)!=BCF_HT_STR )
                error("Only Type=String tags can be used to annotate the ID column\n");
        }
        else if ( (ptr=strstr(str.s,":=")) && (!args->targets_fname || !strncasecmp(ptr+2,"./",2)) )
        {
            *ptr = 0;
            if ( !strncasecmp(str.s,"INFO/",5) && (!strcasecmp(ptr+2,"FILTER") || !strcasecmp(ptr+2,"./FILTER")) )
            {
                // -a not present and transferring filter, needs to be a local transfer
                args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
                annot_col_t *col = &args->cols[args->ncols-1];
                memset(col,0,sizeof(*col));
                col->icol = icol;
                col->replace = replace;
                col->setter = vcf_setter_info_str;
                col->getter = vcf_getter_filter2str_local;
                col->hdr_key_src = strdup(ptr+2);
                col->hdr_key_dst = strdup(str.s+5);
                tmp.l = 0;
                ksprintf(&tmp,"##INFO=<ID=%s,Number=1,Type=String,Description=\"Transferred FILTER column\">",col->hdr_key_dst);
                bcf_hdr_append(args->hdr_out, tmp.s);
                if (bcf_hdr_sync(args->hdr_out) < 0) error_errno("[%s] Failed to update header", __func__);
                int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, col->hdr_key_dst);
                col->number = bcf_hdr_id2length(args->hdr_out,BCF_HL_INFO,hdr_id);
            }
            else
                rename_annots_push(args,ptr+2,str.s);
            *ptr = ':';
        }
        else if ( !strcasecmp("FILTER",str.s) )
        {
            if ( replace & REPLACE_NON_MISSING ) error("Apologies, the -FILTER feature has not been implemented yet.\n");
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            memset(col,0,sizeof(*col));
            col->icol = icol;
            col->replace = replace;
            col->setter = args->tgts_is_vcf ? vcf_setter_filter : setter_filter;
            col->hdr_key_src = strdup(str.s);
            col->hdr_key_dst = strdup(str.s);
            if ( args->tgts_is_vcf )
            {
                bcf_hdr_t *tgts_hdr = args->files->readers[1].header;
                int j;
                for (j=0; j<tgts_hdr->nhrec; j++)
                {
                    bcf_hrec_t *hrec = tgts_hdr->hrec[j];
                    if ( hrec->type!=BCF_HL_FLT ) continue;
                    int k = bcf_hrec_find_key(hrec,"ID");
                    if ( k<0 ) error("[%s] Failed to parse the header, the ID attribute not found", __func__);
                    tmp.l = 0;
                    bcf_hrec_format(hrec, &tmp);
                    bcf_hdr_append(args->hdr_out, tmp.s);
                }
                if (bcf_hdr_sync(args->hdr_out) < 0)
                    error_errno("[%s] Failed to update header", __func__);
            }
        }
        else if ( !strcasecmp("QUAL",str.s) )
        {
            if ( replace & REPLACE_NON_MISSING ) error("Apologies, the -QUAL feature has not been implemented yet.\n");
            if ( replace & SET_OR_APPEND ) error("Apologies, the =QUAL feature has not been implemented yet.\n");
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            memset(col,0,sizeof(*col));
            col->icol = icol;
            col->replace = replace;
            col->setter = args->tgts_is_vcf ? vcf_setter_qual : setter_qual;
            col->hdr_key_src = strdup(str.s);
            col->hdr_key_dst = strdup(str.s);
        }
        else if ( args->tgts_is_vcf && !strcasecmp("INFO",str.s) ) // All INFO fields
        {
            if ( replace & REPLACE_NON_MISSING ) error("Apologies, the -INFO/TAG feature has not been implemented yet.\n");
            if ( replace & SET_OR_APPEND ) error("Apologies, the =INFO feature has not been implemented yet.\n");
            bcf_hdr_t *tgts_hdr = args->files->readers[1].header;
            int j;
            for (j=0; j<tgts_hdr->nhrec; j++)
            {
                bcf_hrec_t *hrec = tgts_hdr->hrec[j];
                if ( hrec->type!=BCF_HL_INFO ) continue;
                int k = bcf_hrec_find_key(hrec,"ID");
                assert( k>=0 ); // this should always be true for valid VCFs
                if ( skip_info && khash_str2int_has_key(skip_info,hrec->vals[k]) ) continue;
                tmp.l = 0;
                bcf_hrec_format(hrec, &tmp);
                bcf_hdr_append(args->hdr_out, tmp.s);
                if (bcf_hdr_sync(args->hdr_out) < 0)
                    error_errno("[%s] Failed to update header", __func__);
                int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, hrec->vals[k]);
                args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
                annot_col_t *col = &args->cols[args->ncols-1];
                memset(col,0,sizeof(*col));
                col->icol = -1;
                col->replace = replace;
                col->hdr_key_src = strdup(hrec->vals[k]);
                col->hdr_key_dst = strdup(hrec->vals[k]);
                col->number  = bcf_hdr_id2length(args->hdr_out,BCF_HL_INFO,hdr_id);
                switch ( bcf_hdr_id2type(args->hdr_out,BCF_HL_INFO,hdr_id) )
                {
                    case BCF_HT_FLAG:   col->setter = vcf_setter_info_flag; break;
                    case BCF_HT_INT:    col->setter = vcf_setter_info_int; break;
                    case BCF_HT_REAL:   col->setter = vcf_setter_info_real; break;
                    case BCF_HT_STR:    col->setter = vcf_setter_info_str; break;
                    default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(args->hdr_out,BCF_HL_INFO,hdr_id));
                }
            }
        }
        else if ( args->tgts_is_vcf && (!strcasecmp("FORMAT",str.s) || !strcasecmp("FMT",str.s)) ) // All FORMAT fields
        {
            bcf_hdr_t *tgts_hdr = args->files->readers[1].header;
            need_sample_map = 1;
            int j;
            for (j=0; j<tgts_hdr->nhrec; j++)
            {
                bcf_hrec_t *hrec = tgts_hdr->hrec[j];
                if ( hrec->type!=BCF_HL_FMT) continue;
                int k = bcf_hrec_find_key(hrec,"ID");
                assert( k>=0 ); // this should always be true for valid VCFs
                if ( skip_fmt && khash_str2int_has_key(skip_fmt,hrec->vals[k]) ) continue;
                tmp.l = 0;
                bcf_hrec_format(hrec, &tmp);
                bcf_hdr_append(args->hdr_out, tmp.s);
                if (bcf_hdr_sync(args->hdr_out) < 0)
                    error_errno("[%s] Failed to update header", __func__);
                int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, hrec->vals[k]);
                args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
                annot_col_t *col = &args->cols[args->ncols-1];
                memset(col,0,sizeof(*col));
                col->icol = -1;
                col->replace = replace;
                col->hdr_key_src = strdup(hrec->vals[k]);
                col->hdr_key_dst = strdup(hrec->vals[k]);
                if ( !strcasecmp("GT",col->hdr_key_src) )
                {
                    if ( !args->tgts_is_vcf ) error("The FORMAT/GT field can be currently populated only from a VCF\n");
                    col->setter = vcf_setter_format_gt;
                }
                else
                    switch ( bcf_hdr_id2type(args->hdr_out,BCF_HL_FMT,hdr_id) )
                    {
                        case BCF_HT_INT:    col->setter = vcf_setter_format_int; break;
                        case BCF_HT_REAL:   col->setter = vcf_setter_format_real; break;
                        case BCF_HT_STR:    col->setter = vcf_setter_format_str; has_fmt_str = 1; break;
                        default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(args->hdr_out,BCF_HL_FMT,hdr_id));
                    }
                hdr_id = bcf_hdr_id2int(tgts_hdr, BCF_DT_ID, hrec->vals[k]);
                col->number = bcf_hdr_id2length(tgts_hdr,BCF_HL_FMT,hdr_id);
            }
        }
        else if ( !strncasecmp("FORMAT/",str.s, 7) || !strncasecmp("FMT/",str.s,4) )
        {
            char *key_dst = str.s + (!strncasecmp("FMT/",str.s,4) ? 4 : 7);
            char *key_src = strstr(key_dst,":=");
            if ( key_src )
            {
                *key_src = 0;
                key_src += 2;
                if ( !strncasecmp("FORMAT/",key_src,7) ) key_src += 7;
                else if ( !strncasecmp("FMT/",key_src,4) ) key_src += 4;
            }
            else
                key_src = key_dst;
            need_sample_map = 1;
            if ( args->tgts_is_vcf )
            {
                bcf_hrec_t *hrec = bcf_hdr_get_hrec(args->files->readers[1].header, BCF_HL_FMT, "ID", key_src, NULL);
                if ( !hrec ) error("No such annotation \"%s\" in %s\n", key_src,args->targets_fname);
                tmp.l = 0;
                bcf_hrec_format_rename(hrec, key_dst, &tmp);
                bcf_hdr_append(args->hdr_out, tmp.s);
                if (bcf_hdr_sync(args->hdr_out) < 0)
                    error_errno("[%s] Failed to update header", __func__);
            }
            int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, key_dst);
            if ( !bcf_hdr_idinfo_exists(args->hdr_out,BCF_HL_FMT,hdr_id) )
                error("The FORMAT tag \"%s\" is not defined in %s, was the -h option provided?\n", str.s, args->targets_fname);
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            memset(col,0,sizeof(*col));
            if ( !args->tgts_is_vcf )
            {
                col->icol = icol;
                icol += args->nsmpl_annot - 1;
            }
            else
                col->icol = -1;
            col->replace = replace;
            col->hdr_key_src = strdup(key_src);
            col->hdr_key_dst = strdup(key_dst);
            if ( !strcasecmp("GT",key_src) )
            {
                if ( !args->tgts_is_vcf ) error("The FORMAT/GT field can be currently populated only from a VCF\n");
                col->setter = vcf_setter_format_gt;
            }
            else
                switch ( bcf_hdr_id2type(args->hdr_out,BCF_HL_FMT,hdr_id) )
                {
                    case BCF_HT_INT:    col->setter = args->tgts_is_vcf ? vcf_setter_format_int  : setter_format_int; break;
                    case BCF_HT_REAL:   col->setter = args->tgts_is_vcf ? vcf_setter_format_real : setter_format_real; break;
                    case BCF_HT_STR:    col->setter = args->tgts_is_vcf ? vcf_setter_format_str  : setter_format_str; has_fmt_str = 1; break;
                    default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(args->hdr_out,BCF_HL_FMT,hdr_id));
                }
            if ( args->tgts_is_vcf )
            {
                bcf_hdr_t *tgts_hdr = args->files->readers[1].header;
                hdr_id = bcf_hdr_id2int(tgts_hdr, BCF_DT_ID, col->hdr_key_src);
                col->number = bcf_hdr_id2length(tgts_hdr,BCF_HL_FMT,hdr_id);
            }
        }
        else
        {
            if ( replace & REPLACE_NON_MISSING ) error("Apologies, the -INFO/TAG feature has not been implemented yet.\n");
            if ( replace & SET_OR_APPEND )
            {
                if ( args->tgts_is_vcf )
                    error("Error: the =INFO/TAG feature is currently supported only with TAB annotation files and has limitations\n"
                          "       (the annotation type is modified to \"Number=.\" and allele ordering is disregarded)\n");
                fprintf(stderr,"Warning: the =INFO/TAG feature modifies the annotation to \"Number=.\" and disregards allele ordering\n");
            }

            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            memset(col,0,sizeof(*col));
            col->icol = icol;
            col->replace = replace;

            int explicit_src_info = 0;
            int explicit_dst_info = 0;
            char *key_dst;
            if ( !strncasecmp("INFO/",str.s,5) )
            {
                key_dst = str.s + 5;
                explicit_dst_info = 1;
            }
            else if ( !strcasecmp("~INFO/END",str.s) )
            {
                key_dst = str.s + 6;
                explicit_dst_info = 1;
            }
            else
                key_dst = str.s;
            char *key_src = strstr(key_dst,":=");
            if ( key_src )
            {
                *key_src = 0;
                key_src += 2;
                if ( !strncasecmp("INFO/",key_src,5) )
                {
                    key_src += 5;
                    explicit_src_info = 1;
                }
                else if ( !strncasecmp("FMT/",key_src,4) || !strncasecmp("FORMAT/",key_src,5) )
                {
                    key_src[-2] = ':';
                    error("Did you mean \"FMT/%s\" rather than \"%s\"?\n",str.s,str.s);
                }
                else if ( !strcasecmp("FILTER",key_src) && args->tgts_is_vcf )
                {
                    col->getter = vcf_getter_filter2str;
                }
            }
            else
                key_src = key_dst;

            col->hdr_key_src = strdup(key_src);
            col->hdr_key_dst = strdup(key_dst);

            int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, key_dst);
            if ( !bcf_hdr_idinfo_exists(args->hdr_out,BCF_HL_INFO,hdr_id) )
            {
                if ( args->tgts_is_vcf ) // reading annotations from a VCF, add a new header line
                {
                    if ( !strcasecmp("ID",key_src) && !explicit_src_info )
                    {
                        // transferring ID column into a new INFO tag
                        tmp.l = 0;
                        ksprintf(&tmp,"##INFO=<ID=%s,Number=1,Type=String,Description=\"Transferred ID column\">",key_dst);
                    }
                    else if ( !strcasecmp("FILTER",key_src) && !explicit_src_info )
                    {
                        // transferring FILTER column into a new INFO tag
                        tmp.l = 0;
                        ksprintf(&tmp,"##INFO=<ID=%s,Number=1,Type=String,Description=\"Transferred FILTER column\">",key_dst);
                    }
                    else
                    {
                        bcf_hrec_t *hrec = bcf_hdr_get_hrec(args->files->readers[1].header, BCF_HL_INFO, "ID", key_src, NULL);
                        if ( !hrec )
                        {
                            if ( explicit_dst_info+explicit_src_info==0 && bcf_hdr_get_hrec(args->files->readers[1].header, BCF_HL_FMT, "ID", key_src, NULL) )
                                error("Did you mean \"FMT/%s\" rather than \"%s\"?\n",str.s,str.s);
                            char *ptr = strchr(key_src,'=');
                            if ( ptr )
                            {
                                *ptr = 0; tmp.l = 0; ksprintf(&tmp,"%s:=%s",key_src,ptr+1); *ptr = '=';
                                error("The INFO tag \"%s\" is not defined, is this what you want \"%s\" ?\n",key_src,tmp.s);
                            }
                            error("The INFO tag \"%s\" is not defined in %s, was the -h option provided?\n", key_src,args->files->readers[1].fname);
                        }
                        tmp.l = 0;
                        bcf_hrec_format_rename(hrec, key_dst, &tmp);
                    }
                    bcf_hdr_append(args->hdr_out, tmp.s);
                    if (bcf_hdr_sync(args->hdr_out) < 0)
                        error_errno("[%s] Failed to update header", __func__);
                    hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, key_dst);
                }
                else
                    error("The INFO tag \"%s\" is not defined in %s, was the -h option provided?\n", key_dst, args->targets_fname);
                assert( bcf_hdr_idinfo_exists(args->hdr_out,BCF_HL_INFO,hdr_id) );
            }
            if  ( args->tgts_is_vcf )
            {
                if ( !strcasecmp("ID",key_src) && !explicit_src_info ) col->getter = vcf_getter_id2str;
                else if ( !strcasecmp("FILTER",key_src) && !explicit_src_info ) col->getter = vcf_getter_filter2str;
            }
            col->number = bcf_hdr_id2length(args->hdr_out,BCF_HL_INFO,hdr_id);
            switch ( bcf_hdr_id2type(args->hdr_out,BCF_HL_INFO,hdr_id) )
            {
                case BCF_HT_FLAG:   col->setter = args->tgts_is_vcf ? vcf_setter_info_flag : setter_info_flag; break;
                case BCF_HT_INT:    col->setter = args->tgts_is_vcf ? vcf_setter_info_int  : setter_info_int; break;
                case BCF_HT_REAL:   col->setter = args->tgts_is_vcf ? vcf_setter_info_real : setter_info_real; break;
                case BCF_HT_STR:    col->setter = args->tgts_is_vcf ? vcf_setter_info_str  : setter_info_str; break;
                default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(args->hdr_out,BCF_HL_INFO,hdr_id));
            }
            if ( replace & SET_OR_APPEND )   // change to Number=.
            {
                bcf_hrec_t *hrec = bcf_hdr_get_hrec(args->hdr_out, BCF_HL_INFO, "ID", key_dst, NULL);
                if ( !hrec ) error("Uh, could not find the new tag \"%s\" in the header\n", key_dst);
                hrec = bcf_hrec_dup(hrec);
                int j = bcf_hrec_find_key(hrec, "Number");
                if ( j<0 ) error("Uh, could not find the entry Number in the header record of %s\n",key_dst);
                free(hrec->vals[j]);
                hrec->vals[j] = strdup(".");
                bcf_hdr_remove(args->hdr_out,BCF_HL_INFO, key_dst);
                bcf_hdr_add_hrec(args->hdr_out, hrec);
            }
        }
        if ( !*se ) break;
        ss = ++se;
    }
    free(str.s);
    free(tmp.s);
    free(args->columns);
    if ( skip_info ) khash_str2int_destroy_free(skip_info);
    if ( skip_fmt ) khash_str2int_destroy_free(skip_fmt);
    if ( has_fmt_str )
    {
        int n = bcf_hdr_nsamples(args->hdr_out);
        if ( args->tgts_is_vcf && n<bcf_hdr_nsamples(args->files->readers[1].header) ) n = bcf_hdr_nsamples(args->files->readers[1].header);
        args->tmpp  = (char**)malloc(sizeof(char*)*n);
        args->tmpp2 = (char**)malloc(sizeof(char*)*n);
    }
    if ( !need_sample_map )
    {
        free(args->sample_map);
        args->sample_map = NULL;
    }
    else if ( sample_map_ok<0 )
        error("No matching samples in source and destination file?\n");
}
static void init_merge_method(args_t *args)
{
    int i;
    for (i=0; i<args->ncols; i++)
    {
        args->cols[i].merge_method = MM_FIRST;
        args->cols[i].mm_str_hash = NULL;
        args->cols[i].mm_dbl = NULL;
        args->cols[i].mm_dbl_nalloc = args->cols[i].mm_dbl_nused = args->cols[i].mm_dbl_ndat = 0;
        memset(&args->cols[i].mm_kstr, 0, sizeof(args->cols[i].mm_kstr));
    }
    if ( !args->merge_method_str.l ) return;
    if ( args->tgts_is_vcf ) error("Error: the --merge-logic is intended for use with BED or TAB-delimited files only.\n");
    if ( !args->tgt_idx && !args->tgts ) error("Error: BEG,END (or FROM,TO) columns or REF,ALT columns are expected with the --merge-logic option.\n");
    char *sb = args->merge_method_str.s;
    while ( *sb )
    {
        char *se = sb;
        while ( *se && *se!=',' ) se++;
        args->tmpks.l = 0;
        kputsn(sb, se-sb, &args->tmpks);
        kputc(0, &args->tmpks);
        char *mm_type_str = args->tmpks.s + args->tmpks.l;
        while ( *mm_type_str!=':' && mm_type_str > args->tmpks.s ) mm_type_str--;
        if ( *mm_type_str!=':' )
            error("Error: could not parse the argument to --merge-logic: %s\n", args->merge_method_str.s);
        *mm_type_str = 0;
        mm_type_str++;
        int mm_type = MM_FIRST;
        if ( !strcasecmp("unique",mm_type_str) ) mm_type = MM_UNIQUE;
        else if ( !strcasecmp("first",mm_type_str) ) mm_type = MM_FIRST;
        else if ( !strcasecmp("append",mm_type_str) ) mm_type = MM_APPEND;
        else if ( !strcasecmp("append-missing",mm_type_str) )
        {
            mm_type = MM_APPEND_MISSING;
            if ( args->ref_idx!=-1 ) args->has_append_mode = 1;
        }
        else if ( !strcasecmp("sum",mm_type_str) ) mm_type = MM_SUM;
        else if ( !strcasecmp("avg",mm_type_str) ) mm_type = MM_AVG;
        else if ( !strcasecmp("min",mm_type_str) ) mm_type = MM_MIN;
        else if ( !strcasecmp("max",mm_type_str) ) mm_type = MM_MAX;
        else error("Error: could not parse --merge-logic %s, the logic \"%s\" is not recognised\n", args->merge_method_str.s,mm_type_str);
        for (i=0; i<args->ncols; i++)
        {
            if ( strcmp(args->cols[i].hdr_key_dst,args->tmpks.s) ) continue;
            if ( (mm_type==MM_APPEND || mm_type==MM_APPEND_MISSING) && args->cols[i].number!=BCF_VL_VAR )
                error("Error: --merge-logic append can be requested only for tags of variable length (Number=.)\n");
            args->cols[i].merge_method = mm_type;
            break;
        }
        if ( i==args->ncols ) error("No such tag in the destination file: %s\n", args->tmpks.s);
        sb = *se ? se + 1 : se;
    }
    if ( args->has_append_mode )
    {
        // create a missing line to insert missing values when VCF ALT finds no match in the annotation file
        args->aline_missing = (annot_line_t*)calloc(1,sizeof(*args->aline_missing));
        int ncol = 0;
        for (i=0; i<args->ncols; i++)
            if ( ncol < args->cols[i].icol + 1 ) ncol = args->cols[i].icol + 1;
        if ( ncol < args->ref_idx + 1 ) ncol = args->ref_idx + 1;
        args->aline_missing->mcols = ncol;
        args->aline_missing->ncols = ncol;
        args->aline_missing->cols = (char**) malloc(ncol*sizeof(char*));
        for (i=0; i<ncol; i++)
            args->aline_missing->cols[i] = strdup(".");
    }
}

static void rename_chrs(args_t *args, char *fname)
{
    int n, i;
    char **map = hts_readlist(fname, 1, &n);
    if ( !map ) error("Could not read: %s\n", fname);
    for (i=0; i<n; i++)
    {
        char *ss = map[i];
        while ( *ss && !isspace(*ss) ) ss++;
        if ( !*ss ) error("Could not parse: %s\n", fname);
        *ss = 0;
        int rid = bcf_hdr_name2id(args->hdr_out, map[i]);
        bcf_hrec_t *hrec = bcf_hdr_get_hrec(args->hdr_out, BCF_HL_CTG, "ID", map[i], NULL);
        if ( !hrec ) continue;  // the sequence not present
        int j = bcf_hrec_find_key(hrec, "ID");
        assert( j>=0 );
        free(hrec->vals[j]);
        ss++;
        while ( *ss && isspace(*ss) ) ss++;
        char *se = ss;
        while ( *se && !isspace(*se) ) se++;
        *se = 0;
        hrec->vals[j] = strdup(ss);
        args->hdr_out->id[BCF_DT_CTG][rid].key = hrec->vals[j];
    }
    for (i=0; i<n; i++) free(map[i]);
    free(map);
}
// Dirty: this relies on bcf_hdr_sync NOT being called
static int rename_annots_core(args_t *args, char *ori_tag, char *new_tag)
{
    int type;
    if ( !strncasecmp("info/",ori_tag,5) ) type = BCF_HL_INFO, ori_tag += 5;
    else if ( !strncasecmp("format/",ori_tag,7) ) type = BCF_HL_FMT, ori_tag += 7;
    else if ( !strncasecmp("fmt/",ori_tag,4) ) type = BCF_HL_FMT, ori_tag += 4;
    else if ( !strncasecmp("filter/",ori_tag,7) ) type = BCF_HL_FLT, ori_tag += 7;
    else return -1;
    if ( !strncasecmp("info/",new_tag,5) )
    {
        if ( type != BCF_HL_INFO ) error("Cannot transfer %s to INFO\n", ori_tag);
        new_tag += 5;
    }
    else if ( !strncasecmp("format/",new_tag,7) )
    {
        if ( type != BCF_HL_FMT ) error("Cannot transfer %s to FORMAT\n", ori_tag);
        new_tag += 7;
    }
    else if ( !strncasecmp("fmt/",new_tag,4) )
    {
        if ( type != BCF_HL_FMT ) error("Cannot transfer %s to FORMAT\n", ori_tag);
        new_tag += 4;
    }
    else if ( !strncasecmp("filter/",new_tag,7) )
    {
        if ( type != BCF_HL_FLT ) error("Cannot transfer %s to FILTER\n", ori_tag);
        new_tag += 7;
    }
    int id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, ori_tag);
    if ( id<0 ) return 1;
    bcf_hrec_t *hrec = bcf_hdr_get_hrec(args->hdr_out, type, "ID", ori_tag, NULL);
    if ( !hrec ) return 1;  // the ID attribute not present
    int j = bcf_hrec_find_key(hrec, "ID");
    assert( j>=0 );
    free(hrec->vals[j]);
    char *ptr = new_tag;
    while ( *ptr && !isspace(*ptr) ) ptr++;
    *ptr = 0;
    hrec->vals[j] = strdup(new_tag);
    args->hdr_out->id[BCF_DT_ID][id].key = hrec->vals[j];
    return 0;
}
static void rename_annots(args_t *args)
{
    int i;
    if ( args->rename_annots )
    {
        args->rename_annots_map = hts_readlist(args->rename_annots, 1, &args->rename_annots_nmap);
        if ( !args->rename_annots_map ) error("Could not read: %s\n", args->rename_annots);
    }
    for (i=0; i<args->rename_annots_nmap; i++)
    {
        char *ptr = args->rename_annots_map[i];
        while ( *ptr && !isspace(*ptr) ) ptr++;
        if ( !*ptr ) error("Could not parse: %s\n", args->rename_annots_map[i]);
        char *rmme = ptr;
        *ptr = 0;
        ptr++;
        while ( *ptr && isspace(*ptr) ) ptr++;
        if ( !*ptr ) { *rmme = ' '; error("Could not parse: %s\n", args->rename_annots_map[i]); }
        if ( rename_annots_core(args, args->rename_annots_map[i], ptr) < 0 )
            error("Cannot rename \"%s\" to \"%s\"\n",args->rename_annots_map[i],ptr);
    }
}
static void rename_annots_push(args_t *args, char *src, char *dst)
{
    args->rename_annots_nmap++;
    args->rename_annots_map = (char**)realloc(args->rename_annots_map,sizeof(*args->rename_annots_map)*args->rename_annots_nmap);
    kstring_t str = {0,0,0};
    ksprintf(&str,"%s %s",src,dst);
    args->rename_annots_map[ args->rename_annots_nmap - 1 ] = str.s;
}
static void init_filters(args_t *args)
{
    // Check if the -i/-e expressions contain external values that should be determined
    // on the fly from the annotation file. The expressions can be given as
    //      TAG={NAME}
    //      TAG={str:NAME}
    //      TAG={int:NAME}
    //      TAG={float:NAME}
    kstring_t str = {0,0,0};
    char *src = strdup(args->filter_str);
    int len = 0;
    while (1)
    {
        char *beg = strchr(src+len,'{');
        if ( !beg ) break;

        // check if "{" appears inside quotes, in such case do not modify
        char skip = 0;
        char *tmp = src;
        while ( tmp<beg )
        {
            if ( tmp[0]!='"' && tmp[0]!='\'' ) { tmp++; continue; }

            // quote character found
            int quote = tmp[0];
            tmp++;
            while ( *tmp && tmp[0]!=quote ) tmp++;
            if ( !*tmp ) error("Could not parse the expression: %s\n",args->filter_str);    // unbalanced quotation; todo: check for escape char
            len = tmp - src + 1;
            skip = 1;
        }
        if ( skip ) continue;

        char *end = ++beg;
        while ( *end && *end!='}' ) end++;
        if ( !*end ) error("Could not parse the expression: %s\n",args->filter_str);
        *end = 0;

        // explicit typing?
        int type = -1;
        tmp = beg;
        while ( *tmp && *tmp!=':' ) tmp++;
        if ( *tmp )
        {
            *tmp = 0;
            if ( !strcasecmp(beg,"str") ) type = BCF_HT_STR;
            else if ( !strcasecmp(beg,"int") ) type = BCF_HT_INT;
            else if ( !strcasecmp(beg,"float") ) type = BCF_HT_REAL;
        }
        args->n_ext++;
        args->ext = (ext_t*)realloc(args->ext,sizeof(*args->ext)*args->n_ext);
        ext_t *ext = &args->ext[args->n_ext-1];
        ext->ht_type = type;
        ext->name = strdup(beg);
        if ( beg-1 > src ) kputsn(src,beg-1-src,&str);
        if ( type==-1 ) kputs("{}",&str);
        else if ( type==BCF_HT_STR ) kputs("{str}",&str);
        else if ( type==BCF_HT_INT ) kputs("{int}",&str);
        else if ( type==BCF_HT_REAL ) kputs("{float}",&str);
        len = str.l;
        kputs(end+1,&str);
        free(src);
        src = strdup(str.s);
        str.l = 0;
    }
    args->filter = filter_init(args->hdr, src);
    free(src);
    free(str.s);

    int i,j,n_ext;
    const int *ext_type = filter_ext_types(args->filter, &n_ext);
    if ( n_ext != args->n_ext )
        error("Failed to parse the expression, unexpected number of dynamic variables (%d vs %d): %s\n",n_ext,args->n_ext,args->filter_str);

    if ( !args->n_ext ) return;

    if ( !args->tgts )
        error("Error: dynamic variables in -i/-e expressions can be currently used only with tab-delimited file, not with VCF (todo)\n");

    // contains external values
    args->ext_ptr = malloc(sizeof(*args->ext_ptr)*args->n_ext);
    for (i=0; i<args->n_ext; i++) args->ext[i].ht_type = ext_type[i];
    args->filter_ext = args->filter;
    args->filter = NULL;

    // set the column idx
    if ( args->ncols )
    {
        for (i=0; i<args->n_ext; i++)
        {
            for (j=0; j<args->ncols; j++)
            {
                if ( strcmp(args->ext[i].name,args->cols[j].hdr_key_src) ) continue;
                args->ext[i].icol = args->cols[j].icol;
                break;
            }
            if ( j==args->ncols ) error("No such column: %s\n",args->ext[i].name);
        }
    }
}

static void init_data(args_t *args)
{
    args->hdr = args->files->readers[0].header;
    args->hdr_out = bcf_hdr_dup(args->hdr);

    if ( args->set_ids_fmt )
    {
        if ( args->set_ids_fmt[0]=='+' ) { args->set_ids_replace = 0; args->set_ids_fmt++; }
        args->set_ids = convert_init(args->hdr_out, NULL, 0, args->set_ids_fmt);
    }
    if ( args->remove_annots ) init_remove_annots(args);
    if ( args->header_fname || args->header_lines ) init_header_lines(args);
    if ( args->targets_fname && args->tgts_is_vcf )
    {
        // reading annots from a VCF
        if ( !bcf_sr_add_reader(args->files, args->targets_fname) )
            error("Failed to open %s: %s\n", args->targets_fname,bcf_sr_strerror(args->files->errnum));
        args->tgts_hdr = args->files->readers[1].header;
    }
    if ( args->columns ) init_columns(args);
    if ( args->targets_fname && !args->tgts_is_vcf )
    {
        if ( !args->columns ) error("The -c option not given\n");
        if ( args->chr_idx==-1 ) error("The -c CHROM option not given\n");
        if ( args->beg_idx==-1 ) error("The -c POS option not given\n");
        if ( args->single_overlaps && args->merge_method_str.l ) error("The options --merge-logic and --single-overlaps cannot be combined\n");
        if ( args->end_idx==-1 || (args->single_overlaps && !args->merge_method_str.l) )
        {
            args->end_idx = -args->beg_idx - 1;
            args->tgts = bcf_sr_regions_init(args->targets_fname,1,args->chr_idx,args->beg_idx,args->end_idx);
            if ( !args->tgts ) error("Could not initialize the annotation file: %s\n", args->targets_fname);
            if ( !args->tgts->tbx ) error("Expected tabix-indexed annotation file: %s\n", args->targets_fname);
        }
        else
        {
            if ( args->ref_idx!=-1 ) error("Error: the REF columns will be ignored when BEG,END (or FROM,TO) is present. Replace END (or TO) with \"-\".\n");
            int len = strlen(args->targets_fname);
            if ( len>=7 && !strcasecmp(".bed.gz",args->targets_fname+len-7) ) args->tgt_is_bed = 1;
            else if ( len>=8 && !strcasecmp(".bed.bgz",args->targets_fname+len-8) ) args->tgt_is_bed = 1;
            else if ( len>=4 && !strcasecmp(".bed",args->targets_fname+len-4) ) args->tgt_is_bed = 1;
            args->tgt_idx = regidx_init(args->targets_fname,parse_with_payload,free_payload,sizeof(char*),args);
            if ( !args->tgt_idx ) error("Failed to parse: %s\n", args->targets_fname);
            args->tgt_itr = regitr_init(args->tgt_idx);
            args->nalines++;
            hts_expand0(annot_line_t,args->nalines,args->malines,args->alines);
        }
        if ( args->min_overlap_str )
        {
            char *tmp = args->min_overlap_str;
            if ( args->min_overlap_str[0] != ':' )
            {
                args->min_overlap_ann = strtod(args->min_overlap_str,&tmp);
                if ( args->min_overlap_ann < 0 || args->min_overlap_ann > 1 || (*tmp && *tmp!=':') )
                    error("Could not parse \"--min-overlap %s\", expected value(s) between 0-1\n", args->min_overlap_str);
            }
            if ( *tmp && *tmp==':' )
            {
                args->min_overlap_vcf = strtod(tmp+1,&tmp);
                if ( args->min_overlap_vcf < 0 || args->min_overlap_vcf > 1 || *tmp )
                    error("Could not parse \"--min-overlap %s\", expected value(s) between 0-1\n", args->min_overlap_str);
            }
        }
    }
    init_merge_method(args);
    args->vcmp = vcmp_init();

    if ( args->filter_str )
        init_filters(args);

    if ( args->mark_sites )
    {
        if ( !args->targets_fname )
        {
            if ( args->mark_sites_logic!=MARK_LISTED ) error("The -a option not given but -%s logic was requested\n",args->mark_sites);
            fprintf(stderr,"Note: The -a option not given, all sites will be annotated with INFO/%s\n",args->mark_sites);
            bcf_hdr_printf(args->hdr_out,"##INFO=<ID=%s,Number=0,Type=Flag,Description=\"Sites marked with `bcftools annotate -m %s`\">",
                    args->mark_sites,args->mark_sites);
        }
        else
            bcf_hdr_printf(args->hdr_out,"##INFO=<ID=%s,Number=0,Type=Flag,Description=\"Sites %slisted in %s\">",
                args->mark_sites,args->mark_sites_logic==MARK_LISTED?"":"not ",args->mark_sites);
    }

    if (args->record_cmd_line) bcf_hdr_append_version(args->hdr_out, args->argc, args->argv, "bcftools_annotate");
    if ( !args->drop_header )
    {
        if ( args->rename_chrs ) rename_chrs(args, args->rename_chrs);
        if ( args->rename_annots || args->rename_annots_map ) rename_annots(args);

        char wmode[8];
        set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
        args->out_fh = hts_open(args->output_fname ? args->output_fname : "-", wmode);
        if ( args->out_fh == NULL ) error("[%s] Error: cannot write to \"%s\": %s\n", __func__,args->output_fname, strerror(errno));
        if ( args->n_threads )
            hts_set_opt(args->out_fh, HTS_OPT_THREAD_POOL, args->files->p);
        if ( bcf_hdr_write(args->out_fh, args->hdr_out)!=0 ) error("[%s] Error: failed to write the header to %s\n", __func__,args->output_fname);
        if ( init_index2(args->out_fh,args->hdr,args->output_fname,
                         &args->index_fn, args->write_index) < 0 )
            error("Error: failed to initialise index for %s\n",args->output_fname);
    }
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->n_ext; i++)
    {
        free(args->ext[i].name);
        if ( args->ext[i].ht_type!=BCF_HT_STR ) continue;
    }
    free(args->ext_ptr);
    free(args->ext);
    for (i=0; i<args->nrm; i++) free(args->rm[i].key);
    free(args->rm);
    if ( args->hdr_out ) bcf_hdr_destroy(args->hdr_out);
    if (args->vcmp) vcmp_destroy(args->vcmp);
    for (i=0; i<args->ncols; i++)
    {
        free(args->cols[i].hdr_key_src);
        free(args->cols[i].hdr_key_dst);
        free(args->cols[i].mm_kstr.s);
        if ( args->cols[i].mm_str_hash ) khash_str2int_destroy_free(args->cols[i].mm_str_hash);
        free(args->cols[i].mm_dbl);
        free(args->cols[i].ptr);
    }
    free(args->cols);
    if ( args->aline_missing )
    {
        for (i=0; i<args->aline_missing->ncols; i++) free(args->aline_missing->cols[i]);
        free(args->aline_missing->cols);
        free(args->aline_missing);
    }
    for (i=0; i<args->malines; i++)
    {
        free(args->alines[i].cols);
        free(args->alines[i].als);
        free(args->alines[i].line.s);
    }
    free(args->alines);
    free(args->srt_alines);
    if ( args->tgt_idx )
    {
        regidx_destroy(args->tgt_idx);
        regitr_destroy(args->tgt_itr);
    }
    if ( args->rename_annots_map )
    {
        for (i=0; i<args->rename_annots_nmap; i++) free(args->rename_annots_map[i]);
        free(args->rename_annots_map);
    }
    if ( args->tgts ) bcf_sr_regions_destroy(args->tgts);
    free(args->tmpks.s);
    free(args->tmpi);
    free(args->tmpf);
    free(args->tmps);
    free(args->tmpp);
    free(args->tmpi2);
    free(args->tmpf2);
    free(args->tmps2);
    free(args->tmpp2);
    free(args->tmpi3);
    free(args->tmpf3);
    free(args->src_smpl_pld);
    free(args->dst_smpl_pld);
    if ( args->set_ids )
        convert_destroy(args->set_ids);
    if ( args->filter ) filter_destroy(args->filter);
    if ( args->filter_ext ) filter_destroy(args->filter_ext);
    if (args->out_fh)
    {
        if ( args->write_index )
        {
            if ( bcf_idx_save(args->out_fh)<0 )
            {
                if ( hts_close(args->out_fh)!=0 ) error("Error: close failed .. %s\n", args->output_fname?args->output_fname:"stdout");
                error("Error: cannot write to index %s\n", args->index_fn);
            }
            free(args->index_fn);
        }
        if ( hts_close(args->out_fh)!=0 ) error("Error: close failed .. %s\n", args->output_fname?args->output_fname:"stdout");
    }
    free(args->sample_map);
    free(args->merge_method_str.s);
}

static void parse_annot_line(args_t *args, char *str, annot_line_t *tmp)
{
    tmp->line.l = 0;
    kputs(str, &tmp->line);
    char *s = tmp->line.s;
    tmp->ncols = 1;
    hts_expand(char*,tmp->ncols,tmp->mcols,tmp->cols);
    tmp->cols[0] = s;
    while ( *s )
    {
        if ( *s=='\t' )
        {
            tmp->ncols++;
            hts_expand(char*,tmp->ncols,tmp->mcols,tmp->cols);
            tmp->cols[tmp->ncols-1] = s+1;
            *s = 0;
        }
        s++;
    }
    if ( args->ref_idx != -1 )
    {
        if ( args->ref_idx >= tmp->ncols )
            error("Could not parse the line, expected %d+ columns, found %d:\n\t%s\n",args->ref_idx+1,tmp->ncols,str);
        if ( args->alt_idx >= tmp->ncols )
            error("Could not parse the line, expected %d+ columns, found %d:\n\t%s\n",args->alt_idx+1,tmp->ncols,str);
        tmp->nals = 2;
        hts_expand(char*,tmp->nals,tmp->mals,tmp->als);
        tmp->als[0] = tmp->cols[args->ref_idx];
        tmp->als[1] = s = tmp->cols[args->alt_idx];
        while ( *s )
        {
            if ( *s==',' )
            {
                tmp->nals++;
                hts_expand(char*,tmp->nals,tmp->mals,tmp->als);
                tmp->als[tmp->nals-1] = s+1;
                *s = 0;
            }
            s++;
        }
    }
}
static void buffer_annot_lines(args_t *args, bcf1_t *line, int start_pos, int end_pos)
{
    if ( args->nalines && args->alines[0].rid != line->rid ) args->nalines = 0;

    int i = 0;
    while ( i<args->nalines )
    {
        if ( line->pos > args->alines[i].end )
        {
            args->nalines--;
            if ( args->nalines && i<args->nalines )
            {
                annot_line_t tmp = args->alines[i];
                memmove(&args->alines[i],&args->alines[i+1],(args->nalines-i)*sizeof(annot_line_t));
                args->alines[args->nalines] = tmp;
            }
        }
        else i++;
    }
    if ( !args->filter_ext && args->ref_idx==-1 && args->nalines ) return;

    while ( !bcf_sr_regions_overlap(args->tgts, bcf_seqname(args->hdr,line), start_pos,end_pos) )
    {
        if ( args->nalines + 1 == 0xffff ) break;   // likely a symbolic allele, don't let the buffer overflow
        args->nalines++;
        hts_expand0(annot_line_t,args->nalines,args->malines,args->alines);
        annot_line_t *tmp = &args->alines[args->nalines-1];
        tmp->rid   = line->rid;
        tmp->start = args->tgts->start;
        tmp->end   = args->tgts->end;
        parse_annot_line(args, args->tgts->line.s, tmp);
        if ( args->filter_ext || args->ref_idx != -1 )
        {
            int iseq = args->tgts->iseq;
            if ( bcf_sr_regions_next(args->tgts)<0 || args->tgts->iseq!=iseq ) break;
        }
        else break;
    }
}

// search string in semicolon separated strings (xx vs aa;bb)
static int str_match(char *needle, char *haystack)
{
    int len = strlen(needle);
    char *ptr = haystack;
    while ( *ptr && (ptr=strstr(ptr,needle)) )
    {
        if ( ptr[len]!=0 && ptr[len]!=';' ) ptr++;          // a prefix, not a match
        else if ( ptr==haystack || ptr[-1]==';' ) return 1; // a match
        ptr++;  // a suffix, not a match
    }
    return 0;
}
// search common string in semicolon separated strings (xx;yy;zz vs aa;bb)
static int strstr_match(char *a, char *b)
{
    char *beg = a;
    while ( *beg )
    {
        char *end = beg;
        while ( *end && *end!=';' ) end++;
        char tmp = *end;
        if ( *end==';' ) *end = 0;
        int ret = str_match(beg,b);
        *end = tmp;
        if ( ret || !*end ) return ret;
        beg = end + 1;
    }
    return 0;
}
static int annotate_from_regidx(args_t *args, bcf1_t *line)
{
    int j;
    int has_overlap = 0;

    for (j=0; j<args->ncols; j++) args->cols[j].done = 0;
    if ( regidx_overlap(args->tgt_idx, bcf_seqname(args->hdr,line),line->pos,line->pos+line->rlen-1, args->tgt_itr) )
    {
        hts_pos_t vcf_end = line->pos + line->rlen - 1;
        while ( regitr_overlap(args->tgt_itr) )
        {
            annot_line_t *tmp = &args->alines[0];
            tmp->rid   = line->rid;
            tmp->start = args->tgt_itr->beg;
            tmp->end   = args->tgt_itr->end;

            // Check min overlap
            int len_ann = tmp->end - tmp->start + 1;
            int len_vcf = line->rlen;
            int isec = (tmp->end < vcf_end ? tmp->end : vcf_end) - (tmp->start > line->pos ? tmp->start : line->pos) + 1;
            assert( isec > 0 );
            if ( args->min_overlap_ann && args->min_overlap_ann > (float)isec/len_ann ) continue;
            if ( args->min_overlap_vcf && args->min_overlap_vcf > (float)isec/len_vcf ) continue;

            parse_annot_line(args, regitr_payload(args->tgt_itr,char*), tmp);

            // If a plain BED file is provided and we are asked to just mark overlapping sites, there are
            // no additional columns. Not sure if there can be any side effects for ill-formatted BED files
            // with variable number of columns
            if ( !args->ncols && args->mark_sites ) has_overlap = 1;

            for (j=0; j<args->ncols; j++)
            {
                if ( args->cols[j].done==1 ) continue;
                int ret = args->cols[j].setter(args,line,&args->cols[j],tmp);
                if ( ret < 0 )
                    error("fixme: Could not set %s at %s:%"PRId64"\n", args->cols[j].hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
                if ( ret==0 )
                    args->cols[j].done = 1;
                has_overlap = 1;
            }
        }
    }
    for (j=0; j<args->ncols; j++)
    {
        if ( args->cols[j].done==1 || args->cols[j].merge_method == MM_FIRST ) continue;
        if ( !args->cols[j].setter ) continue;
        if ( args->cols[j].setter(args,line,&args->cols[j],NULL) < 0 )
            error("fixme: Could not set %s at %s:%"PRId64"\n", args->cols[j].hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
    }
    return has_overlap;
}
static int pass_filter_test_ext(args_t *args, bcf1_t *line, annot_line_t *ann)
{
    char *tmp;
    int i;
    for (i=0; i<args->n_ext; i++)
    {
        int j = args->ext[i].icol;
        if ( args->ext[i].ht_type==BCF_HT_STR ) args->ext_ptr[i] = args->ext[i].s = ann->cols[j];
        else if ( args->ext[i].ht_type==BCF_HT_INT )
        {
            args->ext[i].i = strtol(ann->cols[j],&tmp,10);
            if ( *tmp )
            {
                if ( strcmp(".",ann->cols[j]) ) error("Error: could not parse the annotation file, expected an integer, found \"%s\"\n",ann->cols[j]);
                args->ext_ptr[i] = NULL;
            }
            else
                args->ext_ptr[i] = &args->ext[i].i;
        }
        else if ( args->ext[i].ht_type==BCF_HT_REAL )
        {
            args->ext[i].f = strtod(ann->cols[j],&tmp);
            if ( *tmp )
            {
                if ( strcmp(".",ann->cols[j]) ) error("Error: could not parse the annotation file, expected a float, found \"%s\"\n",ann->cols[j]);
                args->ext_ptr[i] = NULL;
            }
            else
                args->ext_ptr[i] = &args->ext[i].f;
        }
    }
    int pass = filter_test_ext(args->filter_ext,line,NULL,(const void**)args->ext_ptr);
    if ( args->filter_logic==FLT_EXCLUDE ) pass = pass ? 0 : 1;
    return pass;
}
static int annotate_from_tab(args_t *args, bcf1_t *line)
{
    int i,j;
    int has_overlap = 0;

    // Buffer annotation lines. When multiple ALT alleles are present in the annotation file, at least one
    // must match some of the VCF alleles. If the append-missing mode is set (and REF+ALT is requested), the
    // buffered lines will annotate the VCF respecting the order in ALT and when no matching line is found
    // for an ALT, missing value is appended instead.
    int end_pos = line->pos + line->rlen - 1;
    buffer_annot_lines(args, line, line->pos, end_pos);

    args->nsrt_alines = 0;
    hts_expand(uint32_t,args->nalines,args->msrt_alines,args->srt_alines);
    if ( args->nalines >= 0xffff || line->n_allele >= 0xffff )
        error("Error: too many alleles or annotation lines in the buffer at %s:%"PRId64" (todo:skip?)\n",bcf_seqname(args->hdr,line),(int64_t) line->pos+1);

    kstring_t match_end = {0,0,0};
    if ( args->match_end>=0 && bcf_get_info_int32(args->hdr,line,"END",&args->tmpi,&args->mtmpi)==1 )
        kputw(args->tmpi[0],&match_end);

    // Find matching lines
    for (i=0; i<args->nalines; i++)
    {
        if ( line->pos > args->alines[i].end || end_pos < args->alines[i].start ) continue;
        if ( args->ref_idx != -1 )  // REF+ALT matching requested
        {
            if ( line->pos!=args->alines[i].start || vcmp_set_ref(args->vcmp, line->d.allele[0], args->alines[i].als[0]) < 0 ) continue;   // refs are not compatible
            for (j=1; j<args->alines[i].nals; j++)
            {
                int ialt;
                if ( line->n_allele==1 && args->alines[i].als[j][0]=='.' && args->alines[i].als[j][1]==0 )  // match: no ALT allele in VCF and annot file has "."
                    ialt = 0;
                else
                {
                    ialt = vcmp_find_allele(args->vcmp, line->d.allele+1, line->n_allele - 1, args->alines[i].als[j]);
                    if ( ialt < 0 ) continue;
                    ialt++;
                }
                if ( args->match_id>=0 && !strstr_match(line->d.id,args->alines[i].cols[args->match_id]) ) continue;
                if ( args->match_end>=0 && match_end.l && strcmp(match_end.s,args->alines[i].cols[args->match_end]) ) continue;
                if ( args->filter_ext && !pass_filter_test_ext(args,line,&args->alines[i]) ) continue;
                args->srt_alines[args->nsrt_alines++] = (ialt<<16) | i;
                has_overlap = 1;
                break;
            }
        }
        else if ( args->filter_ext )
        {
            if ( pass_filter_test_ext(args,line,&args->alines[i]) )
            {
                args->srt_alines[args->nsrt_alines++] = (0xffff<<16) | i;
                has_overlap = 1;
            }
        }
        else    // overlap, REF+ALT matching not requested
        {
            args->srt_alines[args->nsrt_alines++] = (0xffff<<16) | i;
            has_overlap = 1;
        }
    }

    free(match_end.s);
    if ( !has_overlap && args->filter_ext && !args->keep_sites ) return has_overlap;

    // Sort lines if needed
    if ( args->has_append_mode )
    {
        // insertion sort by VCF ALT index (top bits) and alines index (low bits)
        uint32_t tmp;
        for (i=1; i<args->nsrt_alines; i++)
            for (j=i; j>0 && args->srt_alines[j] < args->srt_alines[j-1]; j--)
                tmp = args->srt_alines[j], args->srt_alines[j] = args->srt_alines[j-1], args->srt_alines[j-1] = tmp;
    }
    // Annotate
    for (j=0; j<args->ncols; j++) args->cols[j].done = 0;
    int ialt_exp = 1;
    for (i=0; i<args->nsrt_alines; i++)
    {
        int ialt = args->srt_alines[i] >> 16;
        int ilin = args->srt_alines[i] & 0xffff;
        if ( args->has_append_mode )
        {
            if ( ialt_exp > ialt ) continue;    // multiple annotation lines for the same position
            if ( ialt_exp < ialt )
            {
                // REF+ALT matching requested, append-missing mode: insert "." if no annotation line was found for the ALT
                while ( ialt_exp++ < ialt )
                {
                    for (j=0; j<args->ncols; j++)
                    {
                        if ( args->cols[j].merge_method != MM_APPEND_MISSING ) continue;
                        if ( args->cols[j].done==1 ) continue;
                        if ( !args->cols[j].setter ) continue;
                        int ret = args->cols[j].setter(args,line,&args->cols[j],args->aline_missing);
                        if ( ret < 0 )
                            error("fixme: Could not set missing %s at %s:%"PRId64"\n", args->cols[j].hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
                        if ( ret==0 )
                            args->cols[j].done = 1;
                    }
                }
            }
        }
        for (j=0; j<args->ncols; j++)
        {
            if ( args->cols[j].done==1 ) continue;
            if ( !args->cols[j].setter ) continue;
            int ret = args->cols[j].setter(args,line,&args->cols[j],&args->alines[ilin]);
            if ( ret < 0 )
                error("fixme: Could not set %s at %s:%"PRId64"\n", args->cols[j].hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
            if ( ret==0 )
                args->cols[j].done = 1;
        }
        ialt_exp = ialt + 1;
    }
    if ( args->nsrt_alines )
    {
        // In the append-missing mode fill missing values to all trailing ALTs, but only if at least one
        // record was found. Otherwise leave the row will be left without annotation.
        if ( args->has_append_mode && ialt_exp < line->n_allele )
        {
            while ( ialt_exp++ < line->n_allele )
            {
                for (j=0; j<args->ncols; j++)
                {
                    if ( args->cols[j].merge_method != MM_APPEND_MISSING ) continue;
                    if ( args->cols[j].done==1 ) continue;
                    if ( !args->cols[j].setter ) continue;
                    int ret = args->cols[j].setter(args,line,&args->cols[j],args->aline_missing);
                    if ( ret < 0 )
                        error("fixme: Could not set missing %s at %s:%"PRId64"\n", args->cols[j].hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
                    if ( ret==0 )
                        args->cols[j].done = 1;
                }
            }
        }
        // Flush
        for (j=0; j<args->ncols; j++)
        {
            if ( args->cols[j].done==1 || args->cols[j].merge_method == MM_FIRST ) continue;
            if ( !args->cols[j].setter ) continue;
            int ret = args->cols[j].setter(args,line,&args->cols[j],NULL);
            if ( ret < 0 )
                error("fixme: Could not set %s at %s:%"PRId64"\n", args->cols[j].hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
        }
    }
    return has_overlap;
}
static int annotate_from_vcf(args_t *args, bcf1_t *line)
{
    if ( !bcf_sr_has_line(args->files,1) ) return 0;
    int j;
    bcf1_t *aline = bcf_sr_get_line(args->files,1);
    for (j=0; j<args->ncols; j++)
    {
        if ( !args->cols[j].setter ) continue;
        if ( args->cols[j].setter(args,line,&args->cols[j],aline) )
            error("fixme: Could not set %s at %s:%"PRId64"\n", args->cols[j].hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
    }
    return 1;
}
static int annotate_from_self(args_t *args, bcf1_t *line)
{
    int j;
    for (j=0; j<args->ncols; j++)
    {
        if ( !args->cols[j].setter ) continue;
        if ( args->cols[j].setter(args,line,&args->cols[j],NULL) )
            error("fixme: Could not set %s at %s:%"PRId64"\n", args->cols[j].hdr_key_src,bcf_seqname(args->hdr,line),(int64_t) line->pos+1);
    }
    return 0;
}
static int annotate_line(args_t *args, bcf1_t *line)
{
    args->current_rec = line;

    int i;
    for (i=0; i<args->nrm; i++)
        args->rm[i].handler(args, line, &args->rm[i]);

    int has_overlap = 0;
    if ( args->tgt_idx )
        has_overlap = annotate_from_regidx(args,line);

    else if ( args->tgts )
        has_overlap = annotate_from_tab(args,line);

    else if ( args->files->nreaders == 2 )
        has_overlap = annotate_from_vcf(args,line);

    else if ( args->ncols )
        has_overlap = annotate_from_self(args,line);

    if ( args->set_ids )
    {
        args->tmpks.l = 0;
        convert_line(args->set_ids, line, &args->tmpks);
        if ( args->tmpks.l )
        {
            int replace = 0;
            if ( args->set_ids_replace ) replace = 1;
            else if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) ) replace = 1;
            if ( replace )
                bcf_update_id(args->hdr_out,line,args->tmpks.s);
        }
    }

    if ( args->mark_sites )
    {
        if ( !args->targets_fname ) has_overlap = 1;

        // ideally, we'd like to be far more general than this in future, see https://github.com/samtools/bcftools/issues/87
        if ( args->mark_sites_logic==MARK_LISTED )
            bcf_update_info_flag(args->hdr_out,line,args->mark_sites,NULL,has_overlap?1:0);
        else
            bcf_update_info_flag(args->hdr_out,line,args->mark_sites,NULL,has_overlap?0:1);
    }

    return has_overlap;
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Annotate and edit VCF/BCF files.\n");
    fprintf(stderr, "Usage:   bcftools annotate [options] VCF\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -a, --annotations FILE          VCF file or tabix-indexed FILE with annotations: CHR\\tPOS[\\tVALUE]+\n");
    fprintf(stderr, "   -c, --columns LIST              List of columns in the annotation file, e.g. CHROM,POS,REF,ALT,-,INFO/TAG. See man page for details\n");
    fprintf(stderr, "   -C, --columns-file FILE         Read -c columns from FILE, one name per row, with optional --merge-logic TYPE: NAME[ TYPE]\n");
    fprintf(stderr, "   -e, --exclude EXPR              Exclude sites for which the expression is true (see man page for details)\n");
    fprintf(stderr, "       --force                     Continue despite parsing error (at your own risk!)\n");
    fprintf(stderr, "   -H, --header-line STR           Header line which should be appended to the VCF header, can be given multiple times\n");
    fprintf(stderr, "   -h, --header-lines FILE         Lines which should be appended to the VCF header\n");
    fprintf(stderr, "   -I, --set-id [+]FORMAT          Set ID column using a `bcftools query`-like expression, see man page for details\n");
    fprintf(stderr, "   -i, --include EXPR              Select sites for which the expression is true (see man page for details)\n");
    fprintf(stderr, "   -k, --keep-sites                Leave -i/-e sites unchanged instead of discarding them\n");
    fprintf(stderr, "   -l, --merge-logic TAG:TYPE      Merge logic for multiple overlapping regions (see man page for details), EXPERIMENTAL\n");
    fprintf(stderr, "   -m, --mark-sites [+-]TAG        Add INFO/TAG flag to sites which are (\"+\") or are not (\"-\") listed in the -a file\n");
    fprintf(stderr, "       --min-overlap ANN:VCF       Required overlap as a fraction of variant in the -a file (ANN), the VCF (:VCF), or reciprocal (ANN:VCF)\n");
    fprintf(stderr, "       --no-version                Do not append version and command line to the header\n");
    fprintf(stderr, "   -o, --output FILE               Write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n");
    fprintf(stderr, "       --pair-logic STR            Matching records by <snps|indels|both|all|some|exact>, see man page for details [some]\n");
    fprintf(stderr, "   -r, --regions REGION            Restrict to comma-separated list of regions\n");
    fprintf(stderr, "   -R, --regions-file FILE         Restrict to regions listed in FILE\n");
    fprintf(stderr, "       --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n");
    fprintf(stderr, "       --rename-annots FILE        Rename annotations: TYPE/old\\tnew, where TYPE is one of FILTER,INFO,FORMAT\n");
    fprintf(stderr, "       --rename-chrs FILE          Rename sequences according to the mapping: old\\tnew\n");
    fprintf(stderr, "   -s, --samples [^]LIST           Comma separated list of samples to annotate (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "   -S, --samples-file [^]FILE      File of samples to annotate (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "       --single-overlaps           Keep memory low by avoiding complexities arising from handling multiple overlapping intervals\n");
    fprintf(stderr, "   -x, --remove LIST               List of annotations (e.g. ID,INFO/DP,FORMAT/DP,FILTER) to remove (or keep with \"^\" prefix). See man page for details\n");
    fprintf(stderr, "       --threads INT               Number of extra output compression threads [0]\n");
    fprintf(stderr, "   -W, --write-index[=FMT]         Automatically index the output files [off]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "   http://samtools.github.io/bcftools/howtos/annotate.html\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfannotate(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->output_fname = "-";
    args->output_type = FT_VCF;
    args->n_threads = 0;
    args->record_cmd_line = 1;
    args->ref_idx = args->alt_idx = args->chr_idx = args->beg_idx = args->end_idx = -1;
    args->set_ids_replace = 1;
    args->match_id = -1;
    args->match_end = -1;
    args->clevel = -1;
    args->pair_logic = -1;
    int regions_is_file = 0;
    int regions_overlap = 1;

    static struct option loptions[] =
    {
        {"keep-sites",no_argument,NULL,'k'},
        {"mark-sites",required_argument,NULL,'m'},
        {"set-id",required_argument,NULL,'I'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"threads",required_argument,NULL,9},
        {"annotations",required_argument,NULL,'a'},
        {"merge-logic",required_argument,NULL,'l'},
        {"collapse",required_argument,NULL,2},
        {"pair-logic",required_argument,NULL,2},
        {"include",required_argument,NULL,'i'},
        {"exclude",required_argument,NULL,'e'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"regions-overlap",required_argument,NULL,3},
        {"remove",required_argument,NULL,'x'},
        {"columns-file",required_argument,NULL,'C'},
        {"columns",required_argument,NULL,'c'},
        {"rename-annots",required_argument,NULL,11},
        {"rename-chrs",required_argument,NULL,1},
        {"header-lines",required_argument,NULL,'h'},
        {"header-line",required_argument,NULL,'H'},
        {"samples",required_argument,NULL,'s'},
        {"samples-file",required_argument,NULL,'S'},
        {"single-overlaps",no_argument,NULL,10},
        {"min-overlap",required_argument,NULL,12},
        {"no-version",no_argument,NULL,8},
        {"force",no_argument,NULL,'f'},
        {"write-index",optional_argument,NULL,'W'},
        {NULL,0,NULL,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "h:H:?o:O:r:R:a:x:c:C:i:e:S:s:I:m:kl:fW::",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'f': args->force = 1; break;
            case 'k': args->keep_sites = 1; break;
            case 'm':
                args->mark_sites_logic = MARK_LISTED;
                if ( optarg[0]=='+' ) args->mark_sites = optarg+1;
                else if ( optarg[0]=='-' ) { args->mark_sites = optarg+1; args->mark_sites_logic = MARK_UNLISTED; }
                else args->mark_sites = optarg;
                break;
            case 'l':
                if ( args->merge_method_str.l ) kputc(',',&args->merge_method_str);
                kputs(optarg,&args->merge_method_str);
                break;
            case 'I': args->set_ids_fmt = optarg; break;
            case 's': args->sample_names = optarg; break;
            case 'S': args->sample_names = optarg; args->sample_is_file = 1; break;
            case 'c': args->columns = strdup(optarg); break;
            case 'C': args->columns = strdup(optarg); args->columns_is_file = 1; break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default:
                    {
                        args->clevel = strtol(optarg,&tmp,10);
                        if ( *tmp || args->clevel<0 || args->clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                    }
                };
                if ( optarg[1] )
                {
                    args->clevel = strtol(optarg+1,&tmp,10);
                    if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                }
                break;
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'x': args->remove_annots = optarg; break;
            case 'a': args->targets_fname = optarg; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 'h': args->header_fname = optarg; break;
            case 'H': args->header_lines = dbuf_push(args->header_lines,strdup(optarg)); break;
            case  1 : args->rename_chrs = optarg; break;
            case  2 :
                if ( args->pair_logic==-1 ) args->pair_logic = 0;
                if ( !strcmp(optarg,"snps") ) args->pair_logic |= BCF_SR_PAIR_SNP_REF;
                else if ( !strcmp(optarg,"indels") ) args->pair_logic |= BCF_SR_PAIR_INDEL_REF;
                else if ( !strcmp(optarg,"both") ) args->pair_logic |= BCF_SR_PAIR_BOTH_REF;
                else if ( !strcmp(optarg,"any") ) args->pair_logic |= BCF_SR_PAIR_ANY;
                else if ( !strcmp(optarg,"all") ) args->pair_logic |= BCF_SR_PAIR_ANY;
                else if ( !strcmp(optarg,"some") ) args->pair_logic |= BCF_SR_PAIR_SOME;
                else if ( !strcmp(optarg,"none") ) args->pair_logic = BCF_SR_PAIR_EXACT;
                else if ( !strcmp(optarg,"exact") ) args->pair_logic = BCF_SR_PAIR_EXACT;
                else error("The --pair-logic string \"%s\" not recognised.\n", optarg);
                break;
            case  3 :
                regions_overlap = parse_overlap_option(optarg);
                if ( regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case  8 : args->record_cmd_line = 0; break;
            case 10 : args->single_overlaps = 1; break;
            case 11 : args->rename_annots = optarg; break;
            case 12 : args->min_overlap_str = optarg; break;
            case 'W':
                if (!(args->write_index = write_index_parse(optarg)))
                    error("Unsupported index format '%s'\n", optarg);
                break;
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else usage(args);
    }
    else fname = argv[optind];

    if ( args->regions_list )
    {
        bcf_sr_set_opt(args->files,BCF_SR_REGIONS_OVERLAP,regions_overlap);
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->targets_fname )
    {
        htsFile *fp = hts_open(args->targets_fname,"r");
        if ( !fp ) error("Failed to open %s\n", args->targets_fname);
        htsFormat type = *hts_get_format(fp);
        hts_close(fp);

        if ( type.format==vcf || type.format==bcf )
        {
            args->tgts_is_vcf = 1;
            args->files->require_index = 1;
            bcf_sr_set_opt(args->files,BCF_SR_PAIR_LOGIC,args->pair_logic>=0 ? args->pair_logic : BCF_SR_PAIR_SOME);
            if ( args->min_overlap_str ) error("The --min-overlap option cannot be used when annotating from a VCF\n");
        }
    }
    if ( args->min_overlap_str && args->single_overlaps ) error("The options --single-overlaps and --min-overlap cannot be combined\n");
    if ( bcf_sr_set_threads(args->files, args->n_threads)<0 ) error("Failed to create threads\n");
    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to read from %s: %s\n", !strcmp("-",fname)?"standard input":fname,bcf_sr_strerror(args->files->errnum));

    static int line_errcode_warned = 0;
    init_data(args);
    while ( bcf_sr_next_line(args->files) )
    {
        if ( !bcf_sr_has_line(args->files,0) ) continue;
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        if ( line->errcode )
        {
            if ( !args->force )
                error("Encountered an error, cannot proceed. Please check the error output above.\n"
                      "If feeling adventurous, use the --force option. (At your own risk!)\n");
            else if ( !line_errcode_warned )
            {
                fprintf(stderr,
                    "Warning: Encountered an error, proceeding only because --force was given.\n"
                    "         Note that this can result in a segfault or a silent corruption of the output file!\n");
                line_errcode_warned = 1;
                line->errcode = 0;
            }
        }
        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass )
            {
                if ( args->keep_sites && bcf_write1(args->out_fh, args->hdr_out, line)!=0 ) error("[%s] Error: failed to write to %s\n", __func__,args->output_fname);
                continue;
            }
        }
        int keep = annotate_line(args, line);
        if ( args->filter_ext && !args->keep_sites && !keep ) continue;
        if ( bcf_write1(args->out_fh, args->hdr_out, line)!=0 ) error("[%s] Error: failed to write to %s\n", __func__,args->output_fname);
    }
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}
