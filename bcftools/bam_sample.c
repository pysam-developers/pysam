/*  bam_sample.c -- group data by sample.

    Copyright (C) 2010, 2011 Broad Institute.
    Copyright (C) 2013, 2016 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>, Petr Danecek <pd3@sanger.ac.uk>

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/khash_str2int.h>
#include <khash_str2str.h>
#include "bam_sample.h"
#include "bcftools.h"


typedef struct
{
    char *fname;
    void *rg2idx;       // hash: read group name to BCF output sample index. Maintained by bsmpl_add_readgroup
    int default_idx;    // default BCF output sample index, set only when all readgroups are treated as one sample
}
file_t;

struct _bam_smpl_t
{
    kstring_t tmp;
    file_t *files;
    int ignore_rg, nsmpl, nfiles;
    char **smpl;        // list of BCF output sample names. Maintained by bsmpl_add_readgroup
    void *sample_list;  // hash: BAM input sample name to BCF output sample name. This is the -s/-S list
    int sample_logic;   // the -s/-S logic, 1: include, 0: exclude
    void *rg_list;      // hash: BAM/rg_id to sample name or */rg_id for global ids. This is the -G list
    int rg_logic;       // the -G logic, 1: include, 0: exclude
    void *name2idx;     // hash: BCF output sample name to BCF output sample index. Maintained by bsmpl_add_readgroup
};

bam_smpl_t *bam_smpl_init(void)
{
    bam_smpl_t *bsmpl;
    bsmpl = (bam_smpl_t*) calloc(1, sizeof(bam_smpl_t));
    bsmpl->name2idx = khash_str2int_init();
    return bsmpl;
}

void bam_smpl_destroy(bam_smpl_t *bsmpl)
{
    if ( !bsmpl ) return;
    if ( bsmpl->name2idx ) khash_str2int_destroy_free(bsmpl->name2idx);
    if ( bsmpl->sample_list ) khash_str2str_destroy_free_all(bsmpl->sample_list);
    if ( bsmpl->rg_list ) khash_str2str_destroy_free_all(bsmpl->rg_list);
    int i;
    for (i=0; i<bsmpl->nfiles; i++)
    {
        file_t *file = &bsmpl->files[i];
        if ( file->rg2idx ) khash_str2int_destroy_free(file->rg2idx);
        free(file->fname);
    }
    free(bsmpl->smpl);
    free(bsmpl->files);
    free(bsmpl->tmp.s);
    free(bsmpl);
}

void bam_smpl_ignore_readgroups(bam_smpl_t* bsmpl)
{
    bsmpl->ignore_rg = 1;
}

static void bsmpl_add_readgroup(bam_smpl_t *bsmpl, file_t *file, const char *rg_id, const char *smpl_name)
{
    int ismpl = -1;
    if ( smpl_name )
    {
        if ( khash_str2int_get(bsmpl->name2idx,smpl_name,&ismpl) < 0 )
        {
            // new sample
            bsmpl->nsmpl++;
            bsmpl->smpl = (char**) realloc(bsmpl->smpl,sizeof(char*)*bsmpl->nsmpl);
            bsmpl->smpl[bsmpl->nsmpl-1] = strdup(smpl_name);
            ismpl = khash_str2int_inc(bsmpl->name2idx,bsmpl->smpl[bsmpl->nsmpl-1]);
        }
    }
    if ( !strcmp("*",rg_id) )
    {
        // all read groups in the bam treated as the same sample
        file->default_idx = ismpl;
        return;
    }
    if ( !file->rg2idx ) file->rg2idx = khash_str2int_init();
    if ( khash_str2int_has_key(file->rg2idx,rg_id) ) return;    // duplicate @RG:ID
    khash_str2int_set(file->rg2idx, strdup(rg_id), ismpl);
}
static int bsmpl_keep_readgroup(bam_smpl_t *bsmpl, file_t *file, const char *rg_id, const char **smpl_name)
{
    char *rg_smpl = khash_str2str_get(bsmpl->rg_list,rg_id);    // unique read group present in one bam only
    if ( !rg_smpl )
    {
        // read group specific to this bam
        bsmpl->tmp.l = 0;
        ksprintf(&bsmpl->tmp,"%s\t%s",rg_id,file->fname);
        rg_smpl = khash_str2str_get(bsmpl->rg_list,bsmpl->tmp.s);
    }
    if ( !rg_smpl )
    {
        // any read group in this file?
        bsmpl->tmp.l = 0;
        ksprintf(&bsmpl->tmp,"*\t%s",file->fname);
        rg_smpl = khash_str2str_get(bsmpl->rg_list,bsmpl->tmp.s);
    }
    if ( !rg_smpl && bsmpl->rg_logic ) return 0;
    if ( rg_smpl && !bsmpl->rg_logic ) return 0;

    if ( rg_smpl && rg_smpl[0]!='\t' ) *smpl_name = rg_smpl;    // rename the sample
    return 1;
}

/*
    The logic of this function is a bit complicated because we want to work
    also with broken bams containing read groups that are not listed in the
    header. The desired behavior is as follows:
        - when -G is given, read groups which are not listed in the header must
          be given explicitly using the "?" symbol in -G.
          Otherwise:
        - if the bam has no header, all reads in the file are assigned to a
          single sample named after the file
        - if there is at least one sample defined in the header, reads with no
          read group id or with a read group id not listed in the header are
          assigned to the first sample encountered in the header
*/
int bam_smpl_add_bam(bam_smpl_t *bsmpl, char *bam_hdr, const char *fname)
{
    bsmpl->nfiles++;
    bsmpl->files = (file_t*) realloc(bsmpl->files,bsmpl->nfiles*sizeof(file_t));
    file_t *file = &bsmpl->files[bsmpl->nfiles-1];
    memset(file,0,sizeof(file_t));
    file->fname  = strdup(fname);
    file->default_idx = -1;

    if ( bsmpl->ignore_rg || !bam_hdr )
    {
        // The option --ignore-RG is set or there is no BAM header: use the file name as the sample name
        bsmpl_add_readgroup(bsmpl,file,"*",file->fname);
        return bsmpl->nfiles-1;
    }

    void *bam_smpls = khash_str2int_init();
    int first_smpl = -1, nskipped = 0;
    const char *p = bam_hdr, *q, *r;
    while ((q = strstr(p, "@RG")) != 0) 
    {
        p = q + 3;
        r = q = 0;
        if ((q = strstr(p, "\tID:")) != 0) q += 4;
        if ((r = strstr(p, "\tSM:")) != 0) r += 4;
        if (r && q)
        {
            char *u, *v;
            int ioq, ior;
            for (u = (char*)q; *u && *u != '\t' && *u != '\n'; ++u);
            for (v = (char*)r; *v && *v != '\t' && *v != '\n'; ++v);
            ioq = *u; ior = *v; *u = *v = '\0';

            // q now points to a null terminated read group id
            // r points to a null terminated sample name
            if ( !strcmp("*",q) || !strcmp("?",q) )
                error("Error: the read group IDs \"*\" and \"?\" have a special meaning in the mpileup code. Please fix the code or the bam: %s\n", fname);

            int accept_rg = 1;
            if ( bsmpl->sample_list )
            {
                // restrict samples based on the -s/-S options
                char *name = khash_str2str_get(bsmpl->sample_list,r);
                if ( bsmpl->sample_logic==0 )
                    accept_rg = name ? 0 : 1;
                else if ( !name )
                    accept_rg = 0;
                else
                    r = name;
            }
            if ( accept_rg && bsmpl->rg_list )
            {
                // restrict readgroups based on the -G option, possibly renaming the sample
                accept_rg = bsmpl_keep_readgroup(bsmpl,file,q,&r);
            }
            if ( accept_rg )
                bsmpl_add_readgroup(bsmpl,file,q,r);
            else
            {
                bsmpl_add_readgroup(bsmpl,file,q,NULL); // ignore this RG but note that it was seen in the header
                nskipped++;
            }

            if ( first_smpl<0 )
                khash_str2int_get(bsmpl->name2idx,r,&first_smpl);
            if ( !khash_str2int_has_key(bam_smpls,r) )
                khash_str2int_inc(bam_smpls,strdup(r));

            *u = ioq; *v = ior;
        }
        else
            break;
        p = q > r ? q : r;
    }
    int nsmpls = khash_str2int_size(bam_smpls);
    khash_str2int_destroy_free(bam_smpls);

    const char *smpl_name = NULL;
    int accept_null_rg = 1;
    if ( bsmpl->rg_list && !bsmpl_keep_readgroup(bsmpl,file,"?",&smpl_name) ) accept_null_rg = 0;
    if ( bsmpl->sample_list && first_smpl==-1 ) accept_null_rg = 0;

    if ( !accept_null_rg && first_smpl==-1 )
    {
        // no suitable read group is available in this bam: ignore the whole file.
        free(file->fname);
        bsmpl->nfiles--;
        return -1;
    }
    if ( !accept_null_rg ) return bsmpl->nfiles-1;
    if ( nsmpls==1 && !nskipped )
    {
        file->default_idx = first_smpl;
        return bsmpl->nfiles-1;
    }
    if ( !smpl_name ) smpl_name = first_smpl==-1 ? file->fname : bsmpl->smpl[first_smpl];

    bsmpl_add_readgroup(bsmpl,file,"?",smpl_name);
    return bsmpl->nfiles-1;
}

const char **bam_smpl_get_samples(bam_smpl_t *bsmpl, int *nsmpl)
{
    *nsmpl = bsmpl->nsmpl;
    return (const char**)bsmpl->smpl;
}

int bam_smpl_get_sample_id(bam_smpl_t *bsmpl, int bam_id, bam1_t *bam_rec)
{
    file_t *file = &bsmpl->files[bam_id];
    if ( file->default_idx >= 0 ) return file->default_idx;

    char *aux_rg = (char*) bam_aux_get(bam_rec, "RG");
    aux_rg = aux_rg ? aux_rg+1 : "?";

    int rg_id;
    if ( khash_str2int_get(file->rg2idx, aux_rg, &rg_id)==0 ) return rg_id;
    if ( khash_str2int_get(file->rg2idx, "?", &rg_id)==0 ) return rg_id;
    return -1;
}

int bam_smpl_add_samples(bam_smpl_t *bsmpl, char *list, int is_file)
{
    if ( list[0]!='^' ) bsmpl->sample_logic = 1;
    else list++;

    int i, nsamples = 0;
    char **samples = hts_readlist(list, is_file, &nsamples);
    if ( !nsamples ) return 0;

    kstring_t ori = {0,0,0};
    kstring_t ren = {0,0,0};

    bsmpl->sample_list = khash_str2str_init();
    for (i=0; i<nsamples; i++)
    {
        char *ptr = samples[i];
        ori.l = ren.l = 0;
        int escaped = 0;
        while ( *ptr )
        {
            if ( *ptr=='\\' && !escaped ) { escaped = 1; ptr++; continue; }
            if ( isspace(*ptr) && !escaped ) break;
            kputc(*ptr, &ori);
            escaped = 0;
            ptr++;
        }
        if ( *ptr )
        {
            while ( *ptr && isspace(*ptr) ) ptr++;
            while ( *ptr )
            {
                if ( *ptr=='\\' && !escaped ) { escaped = 1; ptr++; continue; }
                if ( isspace(*ptr) && !escaped ) break;
                kputc(*ptr, &ren);
                escaped = 0;
                ptr++;
            }
        }
        khash_str2str_set(bsmpl->sample_list,strdup(ori.s),strdup(ren.l?ren.s:ori.s));
        free(samples[i]);
    }
    free(samples);
    free(ori.s);
    free(ren.s);
    return nsamples;
}

int bam_smpl_add_readgroups(bam_smpl_t *bsmpl, char *list, int is_file)
{
    if ( list[0]!='^' ) bsmpl->rg_logic = 1;
    else list++;

    int i, nrows  = 0;
    char **rows = hts_readlist(list, is_file, &nrows);
    if ( !nrows ) return 0;

    kstring_t fld1 = {0,0,0};
    kstring_t fld2 = {0,0,0};
    kstring_t fld3 = {0,0,0};

    bsmpl->rg_list = khash_str2str_init();
    for (i=0; i<nrows; i++)
    {
        char *ptr = rows[i];
        fld1.l = fld2.l = fld3.l = 0;
        int escaped = 0;
        while ( *ptr )
        {
            if ( *ptr=='\\' && !escaped ) { escaped = 1; ptr++; continue; }
            if ( isspace(*ptr) && !escaped ) break;
            kputc(*ptr, &fld1);
            escaped = 0;
            ptr++;
        }
        if ( *ptr )
        {
            while ( *ptr && isspace(*ptr) ) ptr++;
            while ( *ptr )
            {
                if ( *ptr=='\\' && !escaped ) { escaped = 1; ptr++; continue; }
                if ( isspace(*ptr) && !escaped ) break;
                kputc(*ptr, &fld2);
                escaped = 0;
                ptr++;
            }
        }
        if ( *ptr )
        {
            while ( *ptr && isspace(*ptr) ) ptr++;
            while ( *ptr )
            {
                if ( *ptr=='\\' && !escaped ) { escaped = 1; ptr++; continue; }
                if ( isspace(*ptr) && !escaped ) break;
                kputc(*ptr, &fld3);
                escaped = 0;
                ptr++;
            }
        }
        if ( fld3.l )
        {
            // ID FILE SAMPLE
            kputc('\t',&fld1);
            kputs(fld2.s,&fld1);
            fld2.l = 0;
            kputs(fld3.s,&fld2);
        }
        // fld2.s now contains a new sample name. If NULL, use \t to keep the bam header name
        char *value = khash_str2str_get(bsmpl->rg_list,fld1.s);
        if ( !value )
            khash_str2str_set(bsmpl->rg_list,strdup(fld1.s),strdup(fld2.l?fld2.s:"\t"));
        else if ( strcmp(value,fld2.l?fld2.s:"\t") )
            error("Error: The read group \"%s\" was assigned to two different samples: \"%s\" and \"%s\"\n", fld1.s,value,fld2.l?fld2.s:"\t");
        free(rows[i]);
    }
    free(rows);
    free(fld1.s);
    free(fld2.s);
    free(fld3.s);
    return nrows;
}


