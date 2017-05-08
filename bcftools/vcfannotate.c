/*  vcfannotate.c -- Annotate and edit VCF/BCF files.

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
#include <strings.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kseq.h>
#include <htslib/khash_str2int.h>
#include <dlfcn.h>
#include "bcftools.h"
#include "vcmp.h"
#include "filter.h"
#include "convert.h"
#include "smpl_ilist.h"

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

#define REPLACE_MISSING  0  // replace only missing values
#define REPLACE_ALL      1  // replace both missing and existing values
#define REPLACE_NON_MISSING 2  // replace only if tgt is not missing
#define SET_OR_APPEND    3  // set new value if missing or non-existent, append otherwise
typedef struct _annot_col_t
{
    int icol, replace, number;  // number: one of BCF_VL_* types
    char *hdr_key_src, *hdr_key_dst;
    int (*setter)(struct _args_t *, bcf1_t *, struct _annot_col_t *, void*);
}
annot_col_t;

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define MARK_LISTED   1
#define MARK_UNLISTED 2

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr, *hdr_out;
    htsFile *out_fh;
    int output_type, n_threads;
    bcf_sr_regions_t *tgts;

    filter_t *filter;
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE

    rm_tag_t *rm;           // tags scheduled for removal
    int nrm;
    int flt_keep_pass;      // when all filters removed, reset to PASS

    vcmp_t *vcmp;           // for matching annotation and VCF lines by allele
    annot_line_t *alines;   // buffered annotation lines
    int nalines, malines;
    int ref_idx, alt_idx, chr_idx, from_idx, to_idx;   // -1 if not present
    annot_col_t *cols;      // column indexes and setters
    int ncols;

    char *set_ids_fmt;
    convert_t *set_ids;
    int set_ids_replace;

    int nsmpl_annot;
    int *sample_map, nsample_map, sample_is_file;   // map[idst] -> isrc
    int mtmpi, mtmpf, mtmps;
    int mtmpi2, mtmpf2, mtmps2;
    int mtmpi3, mtmpf3, mtmps3;
    int32_t *tmpi, *tmpi2, *tmpi3;
    float *tmpf, *tmpf2, *tmpf3;
    char *tmps, *tmps2, **tmpp, **tmpp2;
    kstring_t tmpks;

    char **argv, *output_fname, *targets_fname, *regions_list, *header_fname;
    char *remove_annots, *columns, *rename_chrs, *sample_names, *mark_sites;
    int argc, drop_header, record_cmd_line, tgts_is_vcf, mark_sites_logic;
}
args_t;

char *msprintf(const char *fmt, ...);

void remove_id(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
    bcf_update_id(args->hdr,line,NULL);
}
void remove_filter(args_t *args, bcf1_t *line, rm_tag_t *tag)
{
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
    if ( nrm ) bcf_hdr_sync(hdr);
}

static void init_remove_annots(args_t *args)
{
    int keep_info = 0, keep_fmt = 0, keep_flt = 0;
    void *keep = khash_str2int_init();
    kstring_t str = {0,0,0};
    char *ss = args->remove_annots;
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
                if ( !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FLT,tag->hdr_id) ) error("Cannot remove %s, not defined in the header.\n", str.s);
                bcf_hdr_remove(args->hdr_out,BCF_HL_FLT,tag->key);
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
                fprintf(stderr,"Warning: The tag \"%s\" not defined in the header\n", str.s);
                args->nrm--;
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
                if ( type==BCF_HL_INFO ) tag->handler = remove_info_tag;
                else if ( type==BCF_HL_FMT ) tag->handler = remove_format_tag;
                bcf_hdr_remove(args->hdr_out,type,tag->key);
            }
        }
        else if ( !strcasecmp("ID",str.s) ) tag->handler = remove_id;
        else if ( !strcasecmp("FILTER",str.s) )
        {
            tag->handler = remove_filter;
            remove_hdr_lines(args->hdr_out,BCF_HL_FLT);
        }
        else if ( !strcasecmp("QUAL",str.s) ) tag->handler = remove_qual;
        else if ( !strcasecmp("INFO",str.s) ) 
        {
            tag->handler = remove_info;
            remove_hdr_lines(args->hdr_out,BCF_HL_INFO);
        }
        else if ( !strcasecmp("FMT",str.s) || !strcasecmp("FORMAT",str.s) )
        {
            tag->handler = remove_format;
            remove_hdr_lines(args->hdr_out,BCF_HL_FMT);
        }
        else if ( str.l )
        {
            if ( str.s[0]=='#' && str.s[1]=='#' )
                bcf_hdr_remove(args->hdr_out,BCF_HL_GEN,str.s+2);
            else
                bcf_hdr_remove(args->hdr_out,BCF_HL_STR,str.s);
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
            bcf_hdr_remove(args->hdr_out,hrec->type,tag->key);
        }
    }
    khash_str2int_destroy_free(keep);
    if ( !args->nrm ) error("No matching tag in -x %s\n", args->remove_annots);
    bcf_hdr_sync(args->hdr_out);
}
static void init_header_lines(args_t *args)
{
    htsFile *file = hts_open(args->header_fname, "rb");
    if ( !file ) error("Error reading %s\n", args->header_fname);
    kstring_t str = {0,0,0};
    while ( hts_getline(file, KS_SEP_LINE, &str) > 0 )
    {
        if ( bcf_hdr_append(args->hdr_out,str.s) ) error("Could not parse %s: %s\n", args->header_fname, str.s);
        bcf_hdr_append(args->hdr,str.s);    // the input file may not have the header line if run with -h (and nothing else)
    }
    hts_close(file);
    free(str.s);
    bcf_hdr_sync(args->hdr_out);
    bcf_hdr_sync(args->hdr);
}
static int setter_filter(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    // note: so far this works only with one filter, not a list of filters
    annot_line_t *tab = (annot_line_t*) data;
    if ( tab->cols[col->icol] && tab->cols[col->icol][0]=='.' && !tab->cols[col->icol][1] ) return 0; // don't replace with "."
    hts_expand(int,1,args->mtmpi,args->tmpi);
    args->tmpi[0] = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, tab->cols[col->icol]);
    if ( args->tmpi[0]<0 ) error("The FILTER is not defined in the header: %s\n", tab->cols[col->icol]);
    if ( col->replace==SET_OR_APPEND ) { bcf_add_filter(args->hdr_out,line,args->tmpi[0]); return 0; }
    if ( col->replace!=REPLACE_MISSING )
    {
        bcf_update_filter(args->hdr_out,line,NULL,0);
        bcf_update_filter(args->hdr_out,line,args->tmpi,1); 
        return 0; 
    }
    
    // only update missing FILTER
    if ( !(line->unpacked & BCF_UN_FLT) ) bcf_unpack(line, BCF_UN_FLT);
    if ( !line->d.n_flt )
        bcf_update_filter(args->hdr_out,line,args->tmpi,1);
    return 0;
}
static int vcf_setter_filter(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    int i;
    bcf1_t *rec = (bcf1_t*) data;
    if ( !(rec->unpacked & BCF_UN_FLT) ) bcf_unpack(rec, BCF_UN_FLT);
    if ( !(line->unpacked & BCF_UN_FLT) ) bcf_unpack(line, BCF_UN_FLT);
    if ( !rec->d.n_flt ) return 0;  // don't overwrite with a missing value
    if ( col->replace==SET_OR_APPEND || col->replace==REPLACE_MISSING )
    {
        if ( col->replace==REPLACE_MISSING && line->d.n_flt ) return 0; // only update missing FILTER
        for (i=0; i<rec->d.n_flt; i++)
        {
            const char *flt = bcf_hdr_int2id(args->files->readers[1].header, BCF_DT_ID, rec->d.flt[i]);
            bcf_add_filter(args->hdr_out,line,bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, flt));
        }
        return 0;
    }
    hts_expand(int,rec->d.n_flt,args->mtmpi,args->tmpi);
    for (i=0; i<rec->d.n_flt; i++)
    {
        const char *flt = bcf_hdr_int2id(args->files->readers[1].header, BCF_DT_ID, rec->d.flt[i]);
        args->tmpi[i] = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, flt);
    }
    bcf_update_filter(args->hdr_out,line,NULL,0);
    bcf_update_filter(args->hdr_out,line,args->tmpi,rec->d.n_flt);
    return 0;
}
static int setter_id(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
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
    if ( col->replace==SET_OR_APPEND ) return bcf_add_id(args->hdr_out,line,tab->cols[col->icol]);
    if ( col->replace!=REPLACE_MISSING ) return bcf_update_id(args->hdr_out,line,tab->cols[col->icol]);

    // running with +ID, only update missing ids
    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) )
        return bcf_update_id(args->hdr_out,line,tab->cols[col->icol]);
    return 0;
}
static int vcf_setter_id(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    if ( rec->d.id && rec->d.id[0]=='.' && !rec->d.id[1] ) return 0;    // don't replace with "."
    if ( col->replace==SET_OR_APPEND ) return bcf_add_id(args->hdr_out,line,rec->d.id);
    if ( col->replace!=REPLACE_MISSING ) return bcf_update_id(args->hdr_out,line,rec->d.id);

    // running with +ID, only update missing ids
    if ( !line->d.id || (line->d.id[0]=='.' && !line->d.id[1]) )
        return bcf_update_id(args->hdr_out,line,rec->d.id);
    return 0;
}
static int setter_qual(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol];
    if ( str[0]=='.' && str[1]==0 ) return 0;   // empty

    if ( col->replace==REPLACE_MISSING && !bcf_float_is_missing(line->qual) ) return 0;

    line->qual = strtod(str, &str);
    if ( str == tab->cols[col->icol] )
        error("Could not parse %s at %s:%d .. [%s]\n", col->hdr_key_src,bcf_seqname(args->hdr,line),line->pos+1,tab->cols[col->icol]);
    return 0;
}
static int vcf_setter_qual(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    if ( bcf_float_is_missing(rec->qual) ) return 0;
    if ( col->replace==REPLACE_MISSING && !bcf_float_is_missing(line->qual) ) return 0;
    line->qual = rec->qual;
    return 0;
}
static int setter_info_flag(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol];
    if ( str[0]=='.' && str[1]==0 ) return 0;

    if ( str[0]=='1' && str[1]==0 ) return bcf_update_info_flag(args->hdr_out,line,col->hdr_key_dst,NULL,1);
    if ( str[0]=='0' && str[1]==0 ) return bcf_update_info_flag(args->hdr_out,line,col->hdr_key_dst,NULL,0);
    error("Could not parse %s at %s:%d .. [%s]\n", bcf_seqname(args->hdr,line),line->pos+1,tab->cols[col->icol]);
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
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpi,col->hdr_key_src,bcf_seqname(args->hdr,line),line->pos+1);
    else if ( col->number==BCF_VL_R && ntmpi!=nals && (ntmpi!=1 || args->tmpi[0]!=bcf_int32_missing || args->tmpi[1]!=bcf_int32_vector_end) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpi,col->hdr_key_src,bcf_seqname(args->hdr,line),line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(args->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%d\n");

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
        if ( ntmpi2==ndst && col->replace==REPLACE_MISSING
                && args->tmpi2[i]!=bcf_int32_missing
                && args->tmpi2[i]!=bcf_int32_vector_end ) continue;

        args->tmpi2[i] = args->tmpi[ map[i] ];
    }
    bcf_update_info_int32(args->hdr_out,line,col->hdr_key_dst,args->tmpi2,ndst);
    return 0;
}
static int setter_info_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol], *end = str;
    if ( str[0]=='.' && str[1]==0 ) return 0;

    int ntmpi = 0;
    while ( *end )
    {
        int val = strtol(str, &end, 10); 
        if ( end==str )
            error("Could not parse %s at %s:%d .. [%s]\n", bcf_seqname(args->hdr,line),line->pos+1,tab->cols[col->icol]);
        ntmpi++;
        hts_expand(int32_t,ntmpi,args->mtmpi,args->tmpi);
        args->tmpi[ntmpi-1] = val;
        str = end+1;
    }

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_int32(args,line,col,tab->nals,tab->als,ntmpi);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_int32(args->hdr, line, col->hdr_key_dst, &args->tmpi2, &args->mtmpi2);
        if ( ret>0 && args->tmpi2[0]!=bcf_int32_missing ) return 0;
    }

    bcf_update_info_int32(args->hdr_out,line,col->hdr_key_dst,args->tmpi,ntmpi);
    return 0;
}
static int vcf_setter_info_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int ntmpi = bcf_get_info_int32(args->files->readers[1].header,rec,col->hdr_key_src,&args->tmpi,&args->mtmpi);
    if ( ntmpi < 0 ) return 0;    // nothing to add

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_int32(args,line,col,rec->n_allele,rec->d.allele,ntmpi);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_int32(args->hdr, line, col->hdr_key_dst, &args->tmpi2, &args->mtmpi2);
        if ( ret>0 && args->tmpi2[0]!=bcf_int32_missing ) return 0;
    }

    bcf_update_info_int32(args->hdr_out,line,col->hdr_key_dst,args->tmpi,ntmpi);
    return 0;
}
static int setter_ARinfo_real(args_t *args, bcf1_t *line, annot_col_t *col, int nals, char **als, int ntmpf)
{
    if ( col->number==BCF_VL_A && ntmpf!=nals-1 && (ntmpf!=1 || !bcf_float_is_missing(args->tmpf[0]) || !bcf_float_is_vector_end(args->tmpf[0])) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpf,col->hdr_key_src,bcf_seqname(args->hdr,line),line->pos+1);
    else if ( col->number==BCF_VL_R && ntmpf!=nals && (ntmpf!=1 || !bcf_float_is_missing(args->tmpf[0]) || !bcf_float_is_vector_end(args->tmpf[0])) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", ntmpf,col->hdr_key_src,bcf_seqname(args->hdr,line),line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(args->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%d\n");

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
        if ( ntmpf2==ndst && col->replace==REPLACE_MISSING
                && !bcf_float_is_missing(args->tmpf2[i])
                && !bcf_float_is_vector_end(args->tmpf2[i]) ) continue;

        args->tmpf2[i] = args->tmpf[ map[i] ];
    }
    bcf_update_info_float(args->hdr_out,line,col->hdr_key_dst,args->tmpf2,ndst);
    return 0;
}
static int setter_info_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    char *str = tab->cols[col->icol], *end = str;
    if ( str[0]=='.' && str[1]==0 ) return 0;

    int ntmpf = 0;
    while ( *end )
    {
        double val = strtod(str, &end);
        if ( end==str )
            error("Could not parse %s at %s:%d .. [%s]\n", bcf_seqname(args->hdr,line),line->pos+1,tab->cols[col->icol]);
        ntmpf++;
        hts_expand(float,ntmpf,args->mtmpf,args->tmpf);
        args->tmpf[ntmpf-1] = val;
        str = end+1;
    }

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_real(args,line,col,tab->nals,tab->als,ntmpf);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_float(args->hdr, line, col->hdr_key_dst, &args->tmpf2, &args->mtmpf2);
        if ( ret>0 && !bcf_float_is_missing(args->tmpf2[0]) ) return 0;
    }

    bcf_update_info_float(args->hdr_out,line,col->hdr_key_dst,args->tmpf,ntmpf);
    return 0;
}
static int vcf_setter_info_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int ntmpf = bcf_get_info_float(args->files->readers[1].header,rec,col->hdr_key_src,&args->tmpf,&args->mtmpf);
    if ( ntmpf < 0 ) return 0;    // nothing to add

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_real(args,line,col,rec->n_allele,rec->d.allele,ntmpf);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_float(args->hdr, line, col->hdr_key_dst, &args->tmpf2, &args->mtmpf2);
        if ( ret>0 && !bcf_float_is_missing(args->tmpf2[0]) ) return 0;
    }

    bcf_update_info_float(args->hdr_out,line,col->hdr_key_dst,args->tmpf,ntmpf);
    return 0;
}
int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst); // see vcfmerge.c
static int setter_ARinfo_string(args_t *args, bcf1_t *line, annot_col_t *col, int nals, char **als)
{
    int nsrc = 1, lsrc = 0;
    while ( args->tmps[lsrc] )
    {
        if ( args->tmps[lsrc]==',' ) nsrc++;
        lsrc++;
    }
    if ( col->number==BCF_VL_A && nsrc!=nals-1 && (nsrc!=1 || args->tmps[0]!='.' || args->tmps[1]!=0 ) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", nsrc,col->hdr_key_src,bcf_seqname(args->hdr,line),line->pos+1);
    else if ( col->number==BCF_VL_R && nsrc!=nals && (nsrc!=1 || args->tmps[0]!='.' || args->tmps[1]!=0 ) )
        error("Incorrect number of values (%d) for the %s tag at %s:%d\n", nsrc,col->hdr_key_src,bcf_seqname(args->hdr,line),line->pos+1);

    int ndst = col->number==BCF_VL_A ? line->n_allele - 1 : line->n_allele;
    int *map = vcmp_map_ARvalues(args->vcmp,ndst,nals,als,line->n_allele,line->d.allele);
    if ( !map ) error("REF alleles not compatible at %s:%d\n");

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
        if ( col->replace==REPLACE_MISSING )
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
        assert( ret==0 );
    }
    bcf_update_info_string(args->hdr_out,line,col->hdr_key_dst,args->tmpks.s);
    return 0;
}
static int setter_info_str(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    int len = strlen(tab->cols[col->icol]);
    if ( !len ) return 0;
    hts_expand(char,len+1,args->mtmps,args->tmps);
    memcpy(args->tmps,tab->cols[col->icol],len+1);
    if ( args->tmps[0]=='.' && args->tmps[1]==0 ) return 0;

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_string(args,line,col,tab->nals,tab->als);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_string(args->hdr, line, col->hdr_key_dst, &args->tmps2, &args->mtmps2);
        if ( ret>0 && (args->tmps2[0]!='.' || args->tmps2[1]!=0) ) return 0;
    }

    bcf_update_info_string(args->hdr_out,line,col->hdr_key_dst,args->tmps);
    return 0;
}
static int vcf_setter_info_str(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int ntmps = bcf_get_info_string(args->files->readers[1].header,rec,col->hdr_key_src,&args->tmps,&args->mtmps);
    if ( ntmps < 0 ) return 0;    // nothing to add

    if ( col->number==BCF_VL_A || col->number==BCF_VL_R ) 
        return setter_ARinfo_string(args,line,col,rec->n_allele,rec->d.allele);

    if ( col->replace==REPLACE_MISSING )
    {
        int ret = bcf_get_info_string(args->hdr, line, col->hdr_key_dst, &args->tmps2, &args->mtmps2);
        if ( ret>0 && (args->tmps2[0]!='.' || args->tmps2[1]!=0) ) return 0;
    }

    bcf_update_info_string(args->hdr_out,line,col->hdr_key_dst,args->tmps);
    return 0;
}
static int vcf_setter_format_gt(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int nsrc = bcf_get_genotypes(args->files->readers[1].header,rec,&args->tmpi,&args->mtmpi);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error

    if ( !args->sample_map )
        return bcf_update_genotypes(args->hdr_out,line,args->tmpi,nsrc);

    int i, j, ndst = bcf_get_genotypes(args->hdr,line,&args->tmpi2,&args->mtmpi2);
    if ( ndst > 0 ) ndst /= bcf_hdr_nsamples(args->hdr_out);
    nsrc /= bcf_hdr_nsamples(args->files->readers[1].header);
    if ( ndst<=0 )  // field not present in dst file
    {
        if ( col->replace==REPLACE_NON_MISSING ) return 0;
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
            if ( col->replace==REPLACE_NON_MISSING && bcf_gt_is_missing(dst[0]) ) continue;
            if ( col->replace==REPLACE_MISSING  && !bcf_gt_is_missing(dst[0]) ) continue;
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
            else if ( col->replace==REPLACE_NON_MISSING && bcf_gt_is_missing(ori[0]) ) keep_ori = 1;
            else if ( col->replace==REPLACE_MISSING  && !bcf_gt_is_missing(ori[0]) ) keep_ori = 1;
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
        if ( col->replace==REPLACE_NON_MISSING ) return 0;    // overwrite only if present
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
            if ( col->replace==REPLACE_NON_MISSING ) { if ( dst[0]==bcf_int32_missing ) continue; } 
            else if ( col->replace==REPLACE_MISSING ) { if ( dst[0]!=bcf_int32_missing ) continue; }
            else if ( col->replace==REPLACE_ALL ) { if ( src[0]==bcf_int32_missing ) continue; }
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
            else if ( col->replace==REPLACE_NON_MISSING ) { if ( ori[0]==bcf_int32_missing ) use_new_ann = 0; }
            else if ( col->replace==REPLACE_MISSING ) { if ( ori[0]!=bcf_int32_missing ) use_new_ann = 0; }
            else if ( col->replace==REPLACE_ALL ) { if ( ann[0]==bcf_int32_missing ) use_new_ann = 0; }
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
        if ( col->replace==REPLACE_NON_MISSING ) return 0;    // overwrite only if present
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
            if ( col->replace==REPLACE_NON_MISSING ) { if ( bcf_float_is_missing(dst[0]) ) continue; } 
            else if ( col->replace==REPLACE_MISSING ) { if ( !bcf_float_is_missing(dst[0]) ) continue; }
            else if ( col->replace==REPLACE_ALL ) { if ( bcf_float_is_missing(src[0]) ) continue; }
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
            else if ( col->replace==REPLACE_NON_MISSING ) { if ( bcf_float_is_missing(ori[0]) ) use_new_ann = 0; }
            else if ( col->replace==REPLACE_MISSING ) { if ( !bcf_float_is_missing(ori[0]) ) use_new_ann = 0; }
            else if ( col->replace==REPLACE_ALL ) { if ( bcf_float_is_missing(ann[0]) ) use_new_ann = 0; }
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

        if ( col->replace==REPLACE_NON_MISSING ) { if ( (*dst)[0]=='.' && (*dst)[1]==0 ) continue; } 
        else if ( col->replace==REPLACE_MISSING ) { if ( (*dst)[0]!='.' || (*dst)[1]!=0 ) continue; }
        else if ( col->replace==REPLACE_ALL ) { if ( (*src)[0]=='.' && (*src)[1]==0 ) continue; }
        *dst = *src;
    }
    return bcf_update_format_string(args->hdr_out,line,col->hdr_key_dst,(const char**)args->tmpp2,nsmpl);
}
static int setter_format_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    annot_line_t *tab = (annot_line_t*) data;
    if ( col->icol+args->nsmpl_annot > tab->ncols ) 
        error("Incorrect number of values for %s at %s:%d\n",col->hdr_key_src,bcf_seqname(args->hdr,line),line->pos+1);
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
                error("Could not parse %s at %s:%d .. [%s]\n", col->hdr_key_src,bcf_seqname(args->hdr,line),line->pos+1,tab->cols[col->icol]);

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
    annot_line_t *tab = (annot_line_t*) data;
    if ( col->icol+args->nsmpl_annot > tab->ncols ) 
        error("Incorrect number of values for %s at %s:%d\n",col->hdr_key_src,bcf_seqname(args->hdr,line),line->pos+1);
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
                error("Could not parse %s at %s:%d .. [%s]\n", col->hdr_key_src,bcf_seqname(args->hdr,line),line->pos+1,tab->cols[col->icol]);

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
    annot_line_t *tab = (annot_line_t*) data;
    if ( col->icol+args->nsmpl_annot > tab->ncols ) 
        error("Incorrect number of values for %s at %s:%d\n",col->hdr_key_src,bcf_seqname(args->hdr,line),line->pos+1);

    int ismpl;
    for (ismpl=0; ismpl<args->nsmpl_annot; ismpl++)
        args->tmpp[ismpl] = tab->cols[col->icol + ismpl];

    return core_setter_format_str(args,line,col,args->tmpp);
}
static int vcf_setter_format_int(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int nsrc = bcf_get_format_int32(args->files->readers[1].header,rec,col->hdr_key_src,&args->tmpi,&args->mtmpi);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error
    return core_setter_format_int(args,line,col,args->tmpi,nsrc/bcf_hdr_nsamples(args->files->readers[1].header));
}
static int vcf_setter_format_real(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    int nsrc = bcf_get_format_float(args->files->readers[1].header,rec,col->hdr_key_src,&args->tmpf,&args->mtmpf);
    if ( nsrc==-3 ) return 0;    // the tag is not present
    if ( nsrc<=0 ) return 1;     // error
    return core_setter_format_real(args,line,col,args->tmpf,nsrc/bcf_hdr_nsamples(args->files->readers[1].header));
}

static int vcf_setter_format_str(args_t *args, bcf1_t *line, annot_col_t *col, void *data)
{
    bcf1_t *rec = (bcf1_t*) data;
    args->tmpp[0] = args->tmps;
    int ret = bcf_get_format_string(args->files->readers[1].header,rec,col->hdr_key_src,&args->tmpp,&args->mtmps);
    args->tmps = args->tmpp[0]; // tmps might be realloced
    if ( ret==-3 ) return 0;    // the tag is not present
    if ( ret<=0 ) return 1;     // error
    return core_setter_format_str(args,line,col,args->tmpp);
}
static int init_sample_map(args_t *args, bcf_hdr_t *src, bcf_hdr_t *dst)
{
    int i;
    if ( !args->sample_names )
    {
        args->nsmpl_annot = bcf_hdr_nsamples(dst);

        // tab annotation file, expecting that all samples are present: sample map not needed
        if ( !src ) return 0;

        int nmatch = 0, order_ok = 1;
        for (i=0; i<bcf_hdr_nsamples(src); i++)
        {
            int id = bcf_hdr_id2int(dst, BCF_DT_SAMPLE, src->samples[i]);
            if ( id!=-1 ) 
            {
                nmatch++;
                if ( i!=id ) order_ok = 0;
            }
        }
        if ( bcf_hdr_nsamples(src)==bcf_hdr_nsamples(dst) && nmatch==bcf_hdr_nsamples(src) && order_ok ) return 0;  // not needed
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

    // possible todo: could do with smpl_ilist only
    smpl_ilist_t *ilist = smpl_ilist_init(dst, args->sample_names, args->sample_is_file, SMPL_STRICT);
    if ( !ilist || !ilist->n ) error("Could not parse: %s\n", args->sample_names);
    char **samples = (char**) malloc(sizeof(char*)*ilist->n);
    for (i=0; i<ilist->n; i++) samples[i] = strdup(dst->samples[i]);
    args->nsmpl_annot = ilist->n;
    smpl_ilist_destroy(ilist);
    int need_sample_map = args->nsmpl_annot==bcf_hdr_nsamples(dst) ? 0 : 1;
    if ( !src )
    {
        // tab annotation file
        for (i=0; i<args->nsmpl_annot; i++)
        {
            int idst = bcf_hdr_id2int(dst, BCF_DT_SAMPLE, samples[i]);
            if ( idst==-1 ) error("Sample \"%s\" not found in the destination file\n", samples[i]);
            args->sample_map[idst] = i;
            if ( idst!=i ) need_sample_map = 1;
        }
    }
    else
    {
        // vcf annotation file
        for (i=0; i<args->nsmpl_annot; i++)
        {
            int isrc, idst;
            char *ss = samples[i], *se = samples[i];
            while ( *se && !isspace(*se) ) se++;
            if ( !*se ) 
            {
                // only one sample name
                isrc = bcf_hdr_id2int(src, BCF_DT_SAMPLE,ss);
                if ( isrc==-1 ) error("Sample \"%s\" not found in the source file\n", ss);
                idst = bcf_hdr_id2int(dst, BCF_DT_SAMPLE,ss);
                if ( idst==-1 ) error("Sample \"%s\" not found in the destination file\n", ss);
                args->sample_map[idst] = isrc;
                if ( idst!=isrc ) need_sample_map = 1;
                continue;
            }
            *se = 0;
            isrc = bcf_hdr_id2int(src, BCF_DT_SAMPLE,ss);
            if ( isrc==-1 ) error("Sample \"%s\" not found in the source file\n", ss);

            ss = se+1;
            while ( isspace(*ss) ) ss++;
            se = ss;
            while ( *se && !isspace(*se) ) se++;

            idst = bcf_hdr_id2int(dst, BCF_DT_SAMPLE,ss);
            if ( idst==-1 ) error("Sample \"%s\" not found in the destination file\n", ss);

            args->sample_map[idst] = isrc;
            if ( idst!=isrc ) need_sample_map = 1;
        }
    }
    for (i=0; i<args->nsmpl_annot; i++) free(samples[i]);
    free(samples);
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
static void init_columns(args_t *args)
{
    int need_sample_map = 0;
    int sample_map_ok = init_sample_map(args, args->tgts_is_vcf?args->files->readers[1].header:NULL, args->hdr);

    void *skip_fmt = NULL, *skip_info = NULL;
    if ( args->tgts_is_vcf )
        args->columns = columns_complement(args->columns, &skip_info, &skip_fmt);

    kstring_t str = {0,0,0}, tmp = {0,0,0};
    char *ss = args->columns, *se = ss;
    args->ncols = 0;
    int icol = -1, has_fmt_str = 0;
    while ( *ss )
    {
        if ( *se && *se!=',' ) { se++; continue; }
        int replace = REPLACE_ALL;
        if ( *ss=='+' ) { replace = REPLACE_MISSING; ss++; }
        else if ( *ss=='-' ) { replace = REPLACE_NON_MISSING; ss++; }
        else if ( *ss=='=' ) { replace = SET_OR_APPEND; ss++; }
        icol++;
        str.l = 0;
        kputsn(ss, se-ss, &str);
        if ( !str.s[0] || !strcasecmp("-",str.s) ) ;
        else if ( !strcasecmp("CHROM",str.s) ) args->chr_idx = icol;
        else if ( !strcasecmp("POS",str.s) ) args->from_idx = icol;
        else if ( !strcasecmp("FROM",str.s) ) args->from_idx = icol;
        else if ( !strcasecmp("TO",str.s) ) args->to_idx = icol;
        else if ( !strcasecmp("REF",str.s) ) args->ref_idx = icol;
        else if ( !strcasecmp("ALT",str.s) ) args->alt_idx = icol;
        else if ( !strcasecmp("ID",str.s) )
        {
            if ( replace==REPLACE_NON_MISSING ) error("Apologies, the -ID feature has not been implemented yet.\n");
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            col->icol = icol;
            col->replace = replace;
            col->setter = args->tgts_is_vcf ? vcf_setter_id : setter_id;
            col->hdr_key_src = strdup(str.s);
            col->hdr_key_dst = strdup(str.s);
        }
        else if ( !strcasecmp("FILTER",str.s) )
        {
            if ( replace==REPLACE_NON_MISSING ) error("Apologies, the -FILTER feature has not been implemented yet.\n");
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
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
                    assert( k>=0 ); // this should always be true for valid VCFs
                    tmp.l = 0;
                    bcf_hrec_format(hrec, &tmp);
                    bcf_hdr_append(args->hdr_out, tmp.s);
                }
                bcf_hdr_sync(args->hdr_out);
            }
        }
        else if ( !strcasecmp("QUAL",str.s) )
        {
            if ( replace==REPLACE_NON_MISSING ) error("Apologies, the -QUAL feature has not been implemented yet.\n");
            if ( replace==SET_OR_APPEND ) error("Apologies, the =QUAL feature has not been implemented yet.\n");
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            col->icol = icol;
            col->replace = replace;
            col->setter = args->tgts_is_vcf ? vcf_setter_qual : setter_qual;
            col->hdr_key_src = strdup(str.s);
            col->hdr_key_dst = strdup(str.s);
        }
        else if ( args->tgts_is_vcf && !strcasecmp("INFO",str.s) ) // All INFO fields
        {
            if ( replace==REPLACE_NON_MISSING ) error("Apologies, the -INFO/TAG feature has not been implemented yet.\n");
            if ( replace==SET_OR_APPEND ) error("Apologies, the =INFO/TAG feature has not been implemented yet.\n");
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
                bcf_hdr_sync(args->hdr_out);
                int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, hrec->vals[k]);
                args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
                annot_col_t *col = &args->cols[args->ncols-1];
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
                bcf_hdr_sync(args->hdr_out);
                int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, hrec->vals[k]);
                args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
                annot_col_t *col = &args->cols[args->ncols-1];
                col->icol = -1;
                col->replace = replace;
                col->hdr_key_src = strdup(hrec->vals[k]);
                col->hdr_key_dst = strdup(hrec->vals[k]);
                if ( !strcasecmp("GT",col->hdr_key_src) ) col->setter = vcf_setter_format_gt;
                else
                    switch ( bcf_hdr_id2type(args->hdr_out,BCF_HL_FMT,hdr_id) )
                    {
                        case BCF_HT_INT:    col->setter = vcf_setter_format_int; break;
                        case BCF_HT_REAL:   col->setter = vcf_setter_format_real; break;
                        case BCF_HT_STR:    col->setter = vcf_setter_format_str; has_fmt_str = 1; break;
                        default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(args->hdr_out,BCF_HL_FMT,hdr_id));
                    }
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
                tmp.l = 0;
                bcf_hrec_format_rename(hrec, key_dst, &tmp);
                bcf_hdr_append(args->hdr_out, tmp.s);
                bcf_hdr_sync(args->hdr_out);
            }
            int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, key_dst);
            if ( !bcf_hdr_idinfo_exists(args->hdr_out,BCF_HL_FMT,hdr_id) )
                error("The tag \"%s\" is not defined in %s\n", str.s, args->targets_fname);
            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
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
            if ( !strcasecmp("GT",key_src) ) col->setter = vcf_setter_format_gt;
            else
                switch ( bcf_hdr_id2type(args->hdr_out,BCF_HL_FMT,hdr_id) )
                {
                    case BCF_HT_INT:    col->setter = args->tgts_is_vcf ? vcf_setter_format_int  : setter_format_int; break;
                    case BCF_HT_REAL:   col->setter = args->tgts_is_vcf ? vcf_setter_format_real : setter_format_real; break;
                    case BCF_HT_STR:    col->setter = args->tgts_is_vcf ? vcf_setter_format_str  : setter_format_str; has_fmt_str = 1; break;
                    default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(args->hdr_out,BCF_HL_FMT,hdr_id));
                }
        }
        else
        {
            if ( replace==REPLACE_NON_MISSING ) error("Apologies, the -INFO/TAG feature has not been implemented yet.\n");
            if ( replace==SET_OR_APPEND ) error("Apologies, the =INFO/TAG feature has not been implemented yet.\n");
            char *key_dst = !strncasecmp("INFO/",str.s,5) ? str.s + 5 : str.s;
            char *key_src = strstr(key_dst,":=");
            if ( key_src )
            {
                *key_src = 0;
                key_src += 2;
                if ( !strncasecmp("INFO/",key_src,5) ) key_src += 5;
            }
            else
                key_src = key_dst;
            int hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, key_dst);
            if ( !bcf_hdr_idinfo_exists(args->hdr_out,BCF_HL_INFO,hdr_id) )
            {
                if ( args->tgts_is_vcf ) // reading annotations from a VCF, add a new header line
                {
                    bcf_hrec_t *hrec = bcf_hdr_get_hrec(args->files->readers[1].header, BCF_HL_INFO, "ID", key_src, NULL);
                    if ( !hrec ) error("The tag \"%s\" is not defined in %s\n", str.s,args->files->readers[1].fname);
                    tmp.l = 0;
                    bcf_hrec_format_rename(hrec, key_dst, &tmp);
                    bcf_hdr_append(args->hdr_out, tmp.s);
                    bcf_hdr_sync(args->hdr_out);
                    hdr_id = bcf_hdr_id2int(args->hdr_out, BCF_DT_ID, key_dst);
                }
                else
                    error("The tag \"%s\" is not defined in %s\n", key_src, args->targets_fname);
                assert( bcf_hdr_idinfo_exists(args->hdr_out,BCF_HL_INFO,hdr_id) );
            }

            args->ncols++; args->cols = (annot_col_t*) realloc(args->cols,sizeof(annot_col_t)*args->ncols);
            annot_col_t *col = &args->cols[args->ncols-1];
            col->icol = icol;
            col->replace = replace;
            col->hdr_key_src = strdup(key_src);
            col->hdr_key_dst = strdup(key_dst);
            col->number  = bcf_hdr_id2length(args->hdr_out,BCF_HL_INFO,hdr_id);
            switch ( bcf_hdr_id2type(args->hdr_out,BCF_HL_INFO,hdr_id) )
            {
                case BCF_HT_FLAG:   col->setter = args->tgts_is_vcf ? vcf_setter_info_flag : setter_info_flag; break;
                case BCF_HT_INT:    col->setter = args->tgts_is_vcf ? vcf_setter_info_int  : setter_info_int; break;
                case BCF_HT_REAL:   col->setter = args->tgts_is_vcf ? vcf_setter_info_real : setter_info_real; break;
                case BCF_HT_STR:    col->setter = args->tgts_is_vcf ? vcf_setter_info_str  : setter_info_str; break;
                default: error("The type of %s not recognised (%d)\n", str.s,bcf_hdr_id2type(args->hdr_out,BCF_HL_INFO,hdr_id));
            }
        }
        if ( !*se ) break;
        ss = ++se;
    }
    free(str.s);
    free(tmp.s);
    if ( args->to_idx==-1 ) args->to_idx = args->from_idx;
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

static void init_data(args_t *args)
{
    args->hdr = args->files->readers[0].header;
    args->hdr_out = bcf_hdr_dup(args->hdr);

    if ( args->remove_annots ) init_remove_annots(args);
    if ( args->header_fname ) init_header_lines(args);
    if ( args->targets_fname && args->tgts_is_vcf )
    {
        // reading annots from a VCF
        if ( !bcf_sr_add_reader(args->files, args->targets_fname) )
            error("Failed to open %s: %s\n", args->targets_fname,bcf_sr_strerror(args->files->errnum));
    }
    if ( args->columns ) init_columns(args);
    if ( args->targets_fname && !args->tgts_is_vcf )
    {
        if ( !args->columns ) error("The -c option not given\n");
        if ( args->chr_idx==-1 ) error("The -c CHROM option not given\n");
        if ( args->from_idx==-1 ) error("The -c POS option not given\n");
        if ( args->to_idx==-1 ) args->to_idx = -args->from_idx - 1;

        args->tgts = bcf_sr_regions_init(args->targets_fname,1,args->chr_idx,args->from_idx,args->to_idx);
        if ( !args->tgts ) error("Could not initialize the annotation file: %s\n", args->targets_fname);
        if ( !args->tgts->tbx ) error("Expected tabix-indexed annotation file: %s\n", args->targets_fname);
    }
    args->vcmp = vcmp_init();

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);

    if ( args->set_ids_fmt )
    {
        if ( args->set_ids_fmt[0]=='+' ) { args->set_ids_replace = 0; args->set_ids_fmt++; }
        args->set_ids = convert_init(args->hdr_out, NULL, 0, args->set_ids_fmt);
    }

    if ( args->mark_sites )
    {
        if ( !args->targets_fname ) error("The -a option not given\n");
        bcf_hdr_printf(args->hdr_out,"##INFO=<ID=%s,Number=0,Type=Flag,Description=\"Sites %slisted in %s\">",
            args->mark_sites,args->mark_sites_logic==MARK_LISTED?"":"not ",args->mark_sites);
    }

     if (args->record_cmd_line) bcf_hdr_append_version(args->hdr_out, args->argc, args->argv, "bcftools_annotate");
    if ( !args->drop_header )
    {
        if ( args->rename_chrs ) rename_chrs(args, args->rename_chrs);

        args->out_fh = hts_open(args->output_fname,hts_bcf_wmode(args->output_type));
        if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
        if ( args->n_threads )
            hts_set_opt(args->out_fh, HTS_OPT_THREAD_POOL, args->files->p);
        bcf_hdr_write(args->out_fh, args->hdr_out);
    }
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->nrm; i++) free(args->rm[i].key);
    free(args->rm);
    if ( args->hdr_out ) bcf_hdr_destroy(args->hdr_out);
    if (args->vcmp) vcmp_destroy(args->vcmp);
    for (i=0; i<args->ncols; i++)
    {
        free(args->cols[i].hdr_key_src);
        free(args->cols[i].hdr_key_dst);
    }
    free(args->cols);
    for (i=0; i<args->malines; i++)
    {
        free(args->alines[i].cols);
        free(args->alines[i].als);
        free(args->alines[i].line.s);
    }
    free(args->alines);
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
    if ( args->set_ids )
        convert_destroy(args->set_ids);
    if ( args->filter )
        filter_destroy(args->filter);
    if (args->out_fh) hts_close(args->out_fh);
    free(args->sample_map);
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

    if ( args->ref_idx==-1 && args->nalines ) return;

    while ( !bcf_sr_regions_overlap(args->tgts, bcf_seqname(args->hdr,line), start_pos,end_pos) )
    {
        args->nalines++;
        hts_expand0(annot_line_t,args->nalines,args->malines,args->alines);
        annot_line_t *tmp = &args->alines[args->nalines-1];
        tmp->rid   = line->rid;
        tmp->start = args->tgts->start;
        tmp->end   = args->tgts->end;
        tmp->line.l = 0;
        kputs(args->tgts->line.s, &tmp->line);
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
                error("Could not parse the line, expected %d+ columns, found %d:\n\t%s\n",args->ref_idx+1,tmp->ncols,args->tgts->line.s);
            if ( args->alt_idx >= tmp->ncols )
                error("Could not parse the line, expected %d+ columns, found %d:\n\t%s\n",args->alt_idx+1,tmp->ncols,args->tgts->line.s);
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
            int iseq = args->tgts->iseq;
            if ( bcf_sr_regions_next(args->tgts)<0 || args->tgts->iseq!=iseq ) break;
        }
        else break;
    }
}

static void annotate(args_t *args, bcf1_t *line)
{
    int i, j;
    for (i=0; i<args->nrm; i++)
        args->rm[i].handler(args, line, &args->rm[i]);

    if ( args->tgts )
    {
        // Buffer annotation lines. When multiple ALT alleles are present in the
        // annotation file, at least one must match one of the VCF alleles.
        int len = 0;
        bcf_get_variant_types(line);
        for (i=1; i<line->n_allele; i++)
            if ( len > line->d.var[i].n ) len = line->d.var[i].n;
        int end_pos = len<0 ? line->pos - len : line->pos;
        buffer_annot_lines(args, line, line->pos, end_pos);
        for (i=0; i<args->nalines; i++)
        {
            if ( line->pos > args->alines[i].end || end_pos < args->alines[i].start ) continue;
            if ( args->ref_idx != -1 )
            {
                if ( vcmp_set_ref(args->vcmp, line->d.allele[0], args->alines[i].als[0]) < 0 ) continue;   // refs not compatible
                for (j=1; j<args->alines[i].nals; j++)
                {
                    if ( line->n_allele==1 && args->alines[i].als[j][0]=='.' && args->alines[i].als[j][1]==0 ) break;   // no ALT allele in VCF and annot file has "."
                    if ( vcmp_find_allele(args->vcmp, line->d.allele+1, line->n_allele - 1, args->alines[i].als[j]) >= 0 ) break;
                }
                if ( j==args->alines[i].nals ) continue;    // none of the annot alleles present in VCF's ALT
            }
            break;
        }

        if ( i<args->nalines )
        {
            // there is a matching line
            for (j=0; j<args->ncols; j++)
                if ( args->cols[j].setter(args,line,&args->cols[j],&args->alines[i]) )
                    error("fixme: Could not set %s at %s:%d\n", args->cols[j].hdr_key_src,bcf_seqname(args->hdr,line),line->pos+1);

        }

        if ( args->mark_sites )
        {
            // ideally, we'd like to be far more general than this in future, see https://github.com/samtools/bcftools/issues/87
            if ( args->mark_sites_logic==MARK_LISTED )
                bcf_update_info_flag(args->hdr_out,line,args->mark_sites,NULL,i<args->nalines?1:0);
            else
                bcf_update_info_flag(args->hdr_out,line,args->mark_sites,NULL,i<args->nalines?0:1);
        }
    }
    else if ( args->files->nreaders == 2 )
    {
        if ( bcf_sr_has_line(args->files,1) )
        {
            bcf1_t *aline = bcf_sr_get_line(args->files,1);
            for (j=0; j<args->ncols; j++)
                if ( args->cols[j].setter(args,line,&args->cols[j],aline) )
                    error("fixme: Could not set %s at %s:%d\n", args->cols[j].hdr_key_src,bcf_seqname(args->hdr,line),line->pos+1);

            if ( args->mark_sites )
                bcf_update_info_flag(args->hdr_out,line,args->mark_sites,NULL,args->mark_sites_logic==MARK_LISTED ? 1 : 0);
        }
        else if ( args->mark_sites )
            bcf_update_info_flag(args->hdr_out,line,args->mark_sites,NULL, args->mark_sites_logic==MARK_UNLISTED ? 1 : 0);
    }
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
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Annotate and edit VCF/BCF files.\n");
    fprintf(stderr, "Usage:   bcftools annotate [options] <in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -a, --annotations <file>       VCF file or tabix-indexed file with annotations: CHR\\tPOS[\\tVALUE]+\n");
    fprintf(stderr, "       --collapse <string>        matching records by <snps|indels|both|all|some|none>, see man page for details [some]\n");
    fprintf(stderr, "   -c, --columns <list>           list of columns in the annotation file, e.g. CHROM,POS,REF,ALT,-,INFO/TAG. See man page for details\n");
    fprintf(stderr, "   -e, --exclude <expr>           exclude sites for which the expression is true (see man page for details)\n");
    fprintf(stderr, "   -h, --header-lines <file>      lines which should be appended to the VCF header\n");
    fprintf(stderr, "   -I, --set-id [+]<format>       set ID column, see man page for details\n");
    fprintf(stderr, "   -i, --include <expr>           select sites for which the expression is true (see man page for details)\n");
    fprintf(stderr, "   -m, --mark-sites [+-]<tag>     add INFO/tag flag to sites which are (\"+\") or are not (\"-\") listed in the -a file\n");
    fprintf(stderr, "       --no-version               do not append version and command line to the header\n");
    fprintf(stderr, "   -o, --output <file>            write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "   -r, --regions <region>         restrict to comma-separated list of regions\n");
    fprintf(stderr, "   -R, --regions-file <file>      restrict to regions listed in a file\n");
    fprintf(stderr, "       --rename-chrs <file>       rename sequences according to map file: from\\tto\n");
    fprintf(stderr, "   -s, --samples [^]<list>        comma separated list of samples to annotate (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "   -S, --samples-file [^]<file>   file of samples to annotate (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "   -x, --remove <list>            list of annotations to remove (e.g. ID,INFO/DP,FORMAT/DP,FILTER). See man page for details\n");
    fprintf(stderr, "       --threads <int>            number of extra output compression threads [0]\n");
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
    args->ref_idx = args->alt_idx = args->chr_idx = args->from_idx = args->to_idx = -1;
    args->set_ids_replace = 1;
    int regions_is_file = 0, collapse = 0;

    static struct option loptions[] =
    {
        {"mark-sites",required_argument,NULL,'m'},
        {"set-id",required_argument,NULL,'I'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"threads",required_argument,NULL,9},
        {"annotations",required_argument,NULL,'a'},
        {"collapse",required_argument,NULL,2},
        {"include",required_argument,NULL,'i'},
        {"exclude",required_argument,NULL,'e'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"remove",required_argument,NULL,'x'},
        {"columns",required_argument,NULL,'c'},
        {"rename-chrs",required_argument,NULL,1},
        {"header-lines",required_argument,NULL,'h'},
        {"samples",required_argument,NULL,'s'},
        {"samples-file",required_argument,NULL,'S'},
        {"no-version",no_argument,NULL,8},
        {NULL,0,NULL,0}
    };
    while ((c = getopt_long(argc, argv, "h:?o:O:r:R:a:x:c:i:e:S:s:I:m:",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'm': 
                args->mark_sites_logic = MARK_LISTED;
                if ( optarg[0]=='+' ) args->mark_sites = optarg+1;
                else if ( optarg[0]=='-' ) { args->mark_sites = optarg+1; args->mark_sites_logic = MARK_UNLISTED; }
                else args->mark_sites = optarg; 
                break;
            case 'I': args->set_ids_fmt = optarg; break;
            case 's': args->sample_names = optarg; break;
            case 'S': args->sample_names = optarg; args->sample_is_file = 1; break;
            case 'c': args->columns = strdup(optarg); break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                };
                break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'x': args->remove_annots = optarg; break;
            case 'a': args->targets_fname = optarg; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 'h': args->header_fname = optarg; break;
            case  1 : args->rename_chrs = optarg; break;
            case  2 :
                if ( !strcmp(optarg,"snps") ) collapse |= COLLAPSE_SNPS;
                else if ( !strcmp(optarg,"indels") ) collapse |= COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"both") ) collapse |= COLLAPSE_SNPS | COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"any") ) collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"all") ) collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"some") ) collapse |= COLLAPSE_SOME;
                else if ( !strcmp(optarg,"none") ) collapse = COLLAPSE_NONE;
                else error("The --collapse string \"%s\" not recognised.\n", optarg);
                break;
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case  8 : args->record_cmd_line = 0; break;
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
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->targets_fname )
    {
        htsFile *fp = hts_open(args->targets_fname,"r"); 
        htsFormat type = *hts_get_format(fp);
        hts_close(fp);

        if ( type.format==vcf || type.format==bcf )
        {
            args->tgts_is_vcf = 1;
            args->files->require_index = 1;
            args->files->collapse = collapse ? collapse : COLLAPSE_SOME;
        }
    }
    if ( bcf_sr_set_threads(args->files, args->n_threads)<0 ) error("Failed to create threads\n");
    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->files->errnum));

    init_data(args);
    while ( bcf_sr_next_line(args->files) )
    {
        if ( !bcf_sr_has_line(args->files,0) ) continue;
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        if ( line->errcode ) error("Encountered error, cannot proceed. Please check the error output above.\n");
        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) continue;
        }
        annotate(args, line);
        bcf_write1(args->out_fh, args->hdr_out, line);
    }
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}
