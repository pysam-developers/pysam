/*  vcfconcat.c -- Concatenate or combine VCF/BCF files.

    Copyright (C) 2013-2021 Genome Research Ltd.

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
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <inttypes.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kseq.h>
#include <htslib/bgzf.h>
#include <htslib/tbx.h> // for hts_get_bgzfp()
#include <htslib/thread_pool.h>
#include <sys/time.h>
#include "bcftools.h"

typedef struct _args_t
{
    bcf_srs_t *files;
    htsFile *out_fh;
    int output_type, n_threads, record_cmd_line, clevel;
    bcf_hdr_t *out_hdr;
    int *seen_seq;

    // phasing
    int *start_pos, start_tid, ifname;
    int *swap_phase, nswap, *nmatch, *nmism;
    bcf1_t **buf;
    uint8_t *buf_mask;
    int nbuf, mbuf, prev_chr, min_PQ, prev_pos_check;
    int32_t *GTa, *GTb, mGTa, mGTb, *phase_qual, *phase_set;

    char **argv, *output_fname, *file_list, **fnames, *remove_dups, *regions_list;
    int argc, nfnames, allow_overlaps, phased_concat, regions_is_file, regions_overlap;
    int compact_PS, phase_set_changed, naive_concat, naive_concat_trust_headers;
    int verbose, explicit_output_type, ligate_force, ligate_warn;
    htsThreadPool *tpool;
}
args_t;

static void init_data(args_t *args)
{
    bcf1_t *line = NULL;

    // With phased concat, the chunks overlap and come in the right order.  To
    // avoid opening all files at once, store start positions to recognise need
    // for the next one. This way we can keep only two open chunks at once.
    if ( args->phased_concat )
    {
        args->start_pos = (int*) malloc(sizeof(int)*args->nfnames);
        line = bcf_init();
    }

    if ( args->verbose ) fprintf(stderr,"Checking the headers and starting positions of %d files\n", args->nfnames);
    kstring_t str = {0,0,0};
    int i, prev_chrid = -1;
    for (i=0; i<args->nfnames; i++)
    {
        htsFile *fp = hts_open(args->fnames[i], "r"); if ( !fp ) error("Failed to open: %s\n", args->fnames[i]);
        bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) error("Failed to parse header: %s\n", args->fnames[i]);
        args->out_hdr = bcf_hdr_merge(args->out_hdr,hdr);
        if ( bcf_hdr_nsamples(hdr) != bcf_hdr_nsamples(args->out_hdr) )
            error("Different number of samples in %s. Perhaps \"bcftools merge\" is what you are looking for?\n", args->fnames[i]);

        int j;
        for (j=0; j<bcf_hdr_nsamples(hdr); j++)
            if ( strcmp(args->out_hdr->samples[j],hdr->samples[j]) )
                error("Different sample names in %s. Perhaps \"bcftools merge\" is what you are looking for?\n", args->fnames[i]);

        if ( args->phased_concat )
        {
            int ret = bcf_read(fp, hdr, line);
            if ( ret!=0 ) args->start_pos[i] = -2;  // empty file
            else
            {
                int chrid = bcf_hdr_id2int(args->out_hdr,BCF_DT_CTG,bcf_seqname(hdr,line));
                args->start_pos[i] = chrid==prev_chrid ? line->pos : -1;
                prev_chrid = chrid;
            }
        }
        bcf_hdr_destroy(hdr);
        if ( hts_close(fp)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->fnames[i]);
    }
    free(str.s);
    if ( line ) bcf_destroy(line);

    args->seen_seq = (int*) calloc(args->out_hdr->n[BCF_DT_CTG],sizeof(int));

    if ( args->phased_concat )
    {
        bcf_hdr_append(args->out_hdr,"##FORMAT=<ID=PQ,Number=1,Type=Integer,Description=\"Phasing Quality (bigger is better)\">");
        bcf_hdr_append(args->out_hdr,"##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">");
    }
    if (args->record_cmd_line) bcf_hdr_append_version(args->out_hdr, args->argc, args->argv, "bcftools_concat");
    char wmode[8];
    set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
    args->out_fh = hts_open(args->output_fname ? args->output_fname : "-", wmode);
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    if ( args->allow_overlaps || args->phased_concat )
    {
        args->files = bcf_sr_init();
        args->files->require_index = 1;
    }
    if ( args->n_threads )
    {
        if ( args->files )
        {
            if ( bcf_sr_set_threads(args->files, args->n_threads)<0 ) error("Failed to create threads\n");
            args->tpool = args->files->p;
        }
        else
        {
            args->tpool = (htsThreadPool*) calloc(1, sizeof(htsThreadPool));
            if ( !args->tpool ) error("Failed to allocate memory\n");
            if ( !(args->tpool->pool = hts_tpool_init(args->n_threads)) ) error("Failed to initialize %d threads\n",args->n_threads);
        }
        hts_set_opt(args->out_fh, HTS_OPT_THREAD_POOL, args->tpool);
    }
    if ( bcf_hdr_write(args->out_fh, args->out_hdr)!=0 ) error("[%s] Error: cannot write the header to %s\n", __func__,args->output_fname);

    if ( args->allow_overlaps )
    {
        if ( args->regions_list )
        {
            bcf_sr_set_opt(args->files,BCF_SR_REGIONS_OVERLAP,args->regions_overlap);
            if ( bcf_sr_set_regions(args->files, args->regions_list, args->regions_is_file)<0 )
                error("Failed to read the regions: %s\n", args->regions_list);
        }
        if ( args->remove_dups )
        {
            if ( !strcmp(args->remove_dups,"snps") ) args->files->collapse |= COLLAPSE_SNPS;
            else if ( !strcmp(args->remove_dups,"indels") ) args->files->collapse |= COLLAPSE_INDELS;
            else if ( !strcmp(args->remove_dups,"both") ) args->files->collapse |= COLLAPSE_SNPS | COLLAPSE_INDELS;
            else if ( !strcmp(args->remove_dups,"any") ) args->files->collapse |= COLLAPSE_ANY;
            else if ( !strcmp(args->remove_dups,"all") ) args->files->collapse |= COLLAPSE_ANY;
            else if ( !strcmp(args->remove_dups,"none") ) args->files->collapse = COLLAPSE_NONE;
            else if ( !strcmp(args->remove_dups,"exact") ) args->files->collapse = COLLAPSE_NONE;
            else error("The -D string \"%s\" not recognised.\n", args->remove_dups);
        }
        for (i=0; i<args->nfnames; i++)
            if ( !bcf_sr_add_reader(args->files,args->fnames[i]) ) error("Failed to open %s: %s\n", args->fnames[i],bcf_sr_strerror(args->files->errnum));
    }
    else if ( args->phased_concat )
    {
        // Remove empty files from the list
        int nok = 0;
        while (1)
        {
            while ( nok<args->nfnames && args->start_pos[nok]!=-2 ) nok++;
            if ( nok==args->nfnames ) break;

            i = nok;
            while ( i<args->nfnames && args->start_pos[i]==-2 ) i++;
            if ( i==args->nfnames ) break;

            int tmp = args->start_pos[nok]; args->start_pos[nok] = args->start_pos[i]; args->start_pos[i] = tmp;
            char *str = args->fnames[nok]; args->fnames[nok] = args->fnames[i]; args->fnames[i] = str;
        }
        for (i=nok; i<args->nfnames; i++) free(args->fnames[i]);
        args->nfnames = nok;

        for (i=1; i<args->nfnames; i++)
            if ( args->start_pos[i-1]!=-1 && args->start_pos[i]!=-1 && args->start_pos[i]<args->start_pos[i-1] )
                error("The files not in ascending order: %d in %s, %d in %s\n", args->start_pos[i-1]+1,args->fnames[i-1],args->start_pos[i]+1,args->fnames[i]);

        args->prev_chr = -1;
        args->swap_phase = (int*) calloc(bcf_hdr_nsamples(args->out_hdr),sizeof(int));
        args->nmatch = (int*) calloc(bcf_hdr_nsamples(args->out_hdr),sizeof(int));
        args->nmism  = (int*) calloc(bcf_hdr_nsamples(args->out_hdr),sizeof(int));
        args->phase_qual = (int32_t*) malloc(bcf_hdr_nsamples(args->out_hdr)*sizeof(int32_t));
        args->phase_set  = (int32_t*) malloc(bcf_hdr_nsamples(args->out_hdr)*sizeof(int32_t));
        args->ifname = 0;
    }
}

static void destroy_data(args_t *args)
{
    int i;
    if ( args->out_fh )
    {
        if ( hts_close(args->out_fh)!=0 ) error("hts_close error\n");
    }
    if ( args->tpool && !args->files )
    {
        hts_tpool_destroy(args->tpool->pool);
        free(args->tpool);
    }
    if ( args->files ) bcf_sr_destroy(args->files);
    if ( args->out_hdr ) bcf_hdr_destroy(args->out_hdr);
    free(args->seen_seq);
    free(args->start_pos);
    free(args->swap_phase);
    for (i=0; i<args->mbuf; i++) bcf_destroy(args->buf[i]);
    free(args->buf);
    free(args->buf_mask);
    free(args->GTa);
    free(args->GTb);
    free(args->nmatch);
    free(args->nmism);
    free(args->phase_qual);
    free(args->phase_set);
    for (i=0; i<args->nfnames; i++) free(args->fnames[i]);
    free(args->fnames);
}

int vcf_write_line(htsFile *fp, kstring_t *line);

#define SWAP(type_t, a, b) { type_t t = a; a = b; b = t; }
static void phase_update(args_t *args, bcf_hdr_t *hdr, bcf1_t *rec)
{
    int i, nGTs = bcf_get_genotypes(hdr, rec, &args->GTa, &args->mGTa);
    if ( nGTs <= 0 ) return;    // GT field is not present
    for (i=0; i<bcf_hdr_nsamples(hdr); i++)
    {
        if ( !args->swap_phase[i] ) continue;
        int *gt = &args->GTa[i*2];
        if ( bcf_gt_is_missing(gt[0]) || gt[1]==bcf_int32_vector_end ) continue;
        if ( !bcf_gt_is_phased(gt[1]) ) continue;
        SWAP(int, gt[0], gt[1]);
        gt[1] |= 1;
    }
    bcf_update_genotypes(hdr,rec,args->GTa,nGTs);
}

static void phased_flush(args_t *args)
{
    if ( !args->nbuf ) return;

    bcf_hdr_t *ahdr = args->files->readers[0].header;
    bcf_hdr_t *bhdr = args->files->readers[1].header;

    int i, j, nsmpl = bcf_hdr_nsamples(args->out_hdr);
    static int gt_absent_warned = 0;
    for (i=0; i<args->nbuf; i+=2)
    {
        if ( args->buf_mask[i/2]!=3 ) continue;

        bcf1_t *arec = args->buf[i];
        bcf1_t *brec = args->buf[i+1];

        int nGTs = bcf_get_genotypes(ahdr, arec, &args->GTa, &args->mGTa);
        if ( nGTs < 0 ) 
        {
            if ( !gt_absent_warned )
            {
                fprintf(stderr,"GT is not present at %s:%"PRId64". (This warning is printed only once.)\n", bcf_seqname(ahdr,arec), (int64_t) arec->pos+1);
                gt_absent_warned = 1;
            }
            continue;
        }
        if ( nGTs != 2*nsmpl ) continue;    // not diploid
        nGTs = bcf_get_genotypes(bhdr, brec, &args->GTb, &args->mGTb);
        if ( nGTs < 0 )
        {
            if ( !gt_absent_warned )
            {
                fprintf(stderr,"GT is not present at %s:%"PRId64". (This warning is printed only once.)\n", bcf_seqname(bhdr,brec), (int64_t) brec->pos+1);
                gt_absent_warned = 1;
            }
            continue;
        }
        if ( nGTs != 2*nsmpl ) continue;    // not diploid

        for (j=0; j<nsmpl; j++)
        {
            int *gta = &args->GTa[j*2];
            int *gtb = &args->GTb[j*2];
            if ( gta[1]==bcf_int32_vector_end || gtb[1]==bcf_int32_vector_end ) continue;
            if ( bcf_gt_is_missing(gta[0]) || bcf_gt_is_missing(gta[1]) || bcf_gt_is_missing(gtb[0]) || bcf_gt_is_missing(gtb[1]) ) continue;
            if ( !bcf_gt_is_phased(gta[1]) || !bcf_gt_is_phased(gtb[1]) ) continue;
            if ( bcf_gt_allele(gta[0])==bcf_gt_allele(gta[1]) || bcf_gt_allele(gtb[0])==bcf_gt_allele(gtb[1]) ) continue;
            if ( bcf_gt_allele(gta[0])==bcf_gt_allele(gtb[0]) && bcf_gt_allele(gta[1])==bcf_gt_allele(gtb[1]) )
            {
                if ( args->swap_phase[j] ) args->nmism[j]++; else args->nmatch[j]++;
            }
            if ( bcf_gt_allele(gta[0])==bcf_gt_allele(gtb[1]) && bcf_gt_allele(gta[1])==bcf_gt_allele(gtb[0]) )
            {
                if ( args->swap_phase[j] ) args->nmatch[j]++; else args->nmism[j]++;
            }
        }
    }
    for (i=0; i<args->nbuf/2; i+=2)
    {
        bcf1_t *rec;
        bcf_hdr_t *hdr;
        int mask = args->buf_mask[i/2];
        if ( mask & 1 ) { rec = args->buf[i]; hdr = args->files->readers[0].header; }
        else { rec = args->buf[i+1]; hdr = args->files->readers[1].header; }
        bcf_translate(args->out_hdr, hdr, rec);
        if ( args->nswap && (mask&1) )
            phase_update(args, args->out_hdr, rec);
        if ( !args->compact_PS || args->phase_set_changed )
        {
            bcf_update_format_int32(args->out_hdr,rec,"PS",args->phase_set,nsmpl);
            args->phase_set_changed = 0;
        }
        if ( bcf_write(args->out_fh, args->out_hdr, rec)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);

        if ( rec->pos < args->prev_pos_check ) error("FIXME, disorder: %s:%"PRId64" vs %d  [1]\n", bcf_seqname(hdr,rec),(int64_t)rec->pos+1,args->prev_pos_check+1);
        args->prev_pos_check = rec->pos;
    }
    args->nswap = 0;
    for (j=0; j<nsmpl; j++)
    {
        if ( args->nmatch[j] >= args->nmism[j] )
            args->swap_phase[j] = 0;
        else
        {
            args->swap_phase[j] = 1;
            args->nswap++;
        }
        if ( args->nmatch[j] && args->nmism[j] )
        {
            // Entropy-inspired quality. The factor 0.7 shifts and scales to (0,1)
            double f = (double)args->nmatch[j]/(args->nmatch[j]+args->nmism[j]);
            args->phase_qual[j] = 99*(0.7 + f*log(f) + (1-f)*log(1-f))/0.7;
        }
        else
            args->phase_qual[j] = 99;
        args->nmatch[j] = 0;
        args->nmism[j]  = 0;
    }
    int PQ_printed = 0;
    for (; i<args->nbuf; i+=2)
    {
        bcf1_t *rec;
        bcf_hdr_t *hdr;
        int mask = args->buf_mask[i/2];
        if ( mask & 2 ) { rec = args->buf[i+1]; hdr = args->files->readers[1].header; }
        else { rec = args->buf[i]; hdr = args->files->readers[0].header; }
        bcf_translate(args->out_hdr, hdr, rec);
        if ( !PQ_printed && mask==3 )
        {
            bcf_update_format_int32(args->out_hdr,rec,"PQ",args->phase_qual,nsmpl);
            PQ_printed = 1;
            for (j=0; j<nsmpl; j++)
                if ( args->phase_qual[j] < args->min_PQ ) 
                {
                    args->phase_set[j] = rec->pos+1;
                    args->phase_set_changed = 1;
                }
                else if ( args->compact_PS ) args->phase_set[j] = bcf_int32_missing;
        }
        if ( args->nswap )
            phase_update(args, args->out_hdr, rec);
        if ( !args->compact_PS || args->phase_set_changed )
        {
            bcf_update_format_int32(args->out_hdr,rec,"PS",args->phase_set,nsmpl);
            args->phase_set_changed = 0;
        }
        if ( bcf_write(args->out_fh, args->out_hdr, rec)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);

        if ( rec->pos < args->prev_pos_check ) error("FIXME, disorder: %s:%"PRId64" vs %d  [2]\n", bcf_seqname(hdr,rec),(int64_t)rec->pos+1,args->prev_pos_check+1);
        args->prev_pos_check = rec->pos;
    }
    args->nbuf = 0;
}

static void phased_push(args_t *args, bcf1_t *arec, bcf1_t *brec, int is_overlap)
{
    bcf_hdr_t *ahdr = arec ? bcf_sr_get_header(args->files,0) : NULL;
    bcf_hdr_t *bhdr = brec ? bcf_sr_get_header(args->files,1) : NULL;

    if ( arec && arec->errcode )
        error("Parse error at %s:%"PRId64", cannot proceed: %s\n", bcf_seqname(ahdr,arec),(int64_t) arec->pos+1, args->files->readers[0].fname);
    if ( brec && brec->errcode )
        error("Parse error at %s:%"PRId64", cannot proceed: %s\n", bcf_seqname(bhdr,brec),(int64_t) brec->pos+1, args->files->readers[1].fname);

    int i, nsmpl = bcf_hdr_nsamples(args->out_hdr);
    int chr_id = arec ? bcf_hdr_name2id(args->out_hdr,bcf_seqname(ahdr,arec)) : bcf_hdr_name2id(args->out_hdr,bcf_seqname(bhdr,brec));
    if ( args->prev_chr<0 || args->prev_chr!=chr_id )
    {
        if ( args->prev_chr>=0 ) phased_flush(args);

        for (i=0; i<nsmpl; i++)
            args->phase_set[i] = arec ? arec->pos+1 : brec->pos+1;
        args->phase_set_changed = 1;

        if ( args->seen_seq[chr_id] ) error("The chromosome block %s is not contiguous\n", arec ? bcf_seqname(ahdr,arec) : bcf_seqname(bhdr,brec));
        args->seen_seq[chr_id] = 1;
        args->prev_chr = chr_id;
        args->prev_pos_check = -1;
    }

    if ( !is_overlap )
    {
        assert(arec);

        bcf_translate(args->out_hdr, ahdr, arec);
        if ( args->nswap )
            phase_update(args, args->out_hdr, arec);
        if ( !args->compact_PS || args->phase_set_changed )
        {
            bcf_update_format_int32(args->out_hdr,arec,"PS",args->phase_set,nsmpl);
            args->phase_set_changed = 0;
        }
        if ( bcf_write(args->out_fh, args->out_hdr, arec)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);

        if ( arec->pos < args->prev_pos_check )
            error("FIXME, disorder: %s:%"PRId64" in %s vs %d written  [3]\n", bcf_seqname(ahdr,arec), (int64_t) arec->pos+1,args->files->readers[0].fname, args->prev_pos_check+1);
        args->prev_pos_check = arec->pos;
        return;
    }

    int m = args->mbuf;
    args->nbuf += 2;
    hts_expand(bcf1_t*,args->nbuf,args->mbuf,args->buf);
    if ( m < args->mbuf ) args->buf_mask = (uint8_t*)realloc(args->buf_mask,sizeof(*args->buf_mask)*args->mbuf);
    for (i=m; i<args->mbuf; i++)
        args->buf[i] = bcf_init1();

    if ( arec ) SWAP(bcf1_t*, args->files->readers[0].buffer[0], args->buf[args->nbuf-2]);
    if ( brec ) SWAP(bcf1_t*, args->files->readers[1].buffer[0], args->buf[args->nbuf-1]);
    args->buf_mask[args->nbuf/2-1] = (arec?1:0) | (brec?2:0);
}

static int _get_active_index(bcf_srs_t *sr)
{
    int i;
    for (i=0; i<sr->nreaders; i++)
        if ( bcf_sr_has_line(sr,i) ) return i;
    return -1;
}

static void concat(args_t *args)
{
    static int site_drop_warned = 0;
    int i;
    if ( args->phased_concat )  // phased concat
    {
        // keep only two open files at a time
        while ( args->ifname < args->nfnames )
        {
            int new_file = 0;
            while ( args->files->nreaders < 2 && args->ifname < args->nfnames )
            {
                if ( !bcf_sr_add_reader(args->files,args->fnames[args->ifname]) ) error("Failed to open %s: %s\n", args->fnames[args->ifname],bcf_sr_strerror(args->files->errnum));
                new_file = 1;

                args->ifname++;
                if ( args->start_pos[args->ifname-1]==-1 ) break;   // new chromosome, start with only one file open
                if ( args->ifname < args->nfnames && args->start_pos[args->ifname]==-1 ) break; // next file starts on a different chromosome
            }

            // is there a line from the previous run? Seek the newly opened reader to that position
            int seek_pos = -1;
            int seek_chr = -1;
            if ( bcf_sr_has_line(args->files,0) )
            {
                bcf1_t *line = bcf_sr_get_line(args->files,0);
                bcf_sr_seek(args->files, bcf_seqname(args->files->readers[0].header,line), line->pos);
                seek_pos = line->pos;
                seek_chr = bcf_hdr_name2id(args->out_hdr, bcf_seqname(args->files->readers[0].header,line));
            }
            else if ( new_file )
                bcf_sr_seek(args->files,NULL,0);  // set to start

            int nret, ir;
            while ( (nret = bcf_sr_next_line(args->files)) )
            {
                int is_overlap = args->files->nreaders==1 ? 0 : 1;
                if ( !bcf_sr_has_line(args->files,0) )  // no input from the first reader
                {
                    // We are assuming that there is a perfect overlap, sites which are not present in both files are dropped
                    if ( bcf_sr_region_done(args->files,0) )
                    {
                        phased_flush(args);
                        bcf_sr_remove_reader(args->files, 0);
                        is_overlap = 0;
                    }
                    else if ( args->ligate_warn )
                    {
                        if ( !site_drop_warned )
                        {
                            ir = _get_active_index(args->files);
                            fprintf(stderr,
                                "Warning: Dropping the site %s:%"PRId64". The --ligate option is intended for VCFs with perfect\n"
                                "         overlap, sites in overlapping regions present in one but missing in other are dropped.\n"
                                "         This warning is printed only once.\n",
                                bcf_seqname(bcf_sr_get_header(args->files,ir),bcf_sr_get_line(args->files,ir)), (int64_t) bcf_sr_get_line(args->files,ir)->pos+1);
                            site_drop_warned = 1;
                        }
                        continue;
                    }
                    else if ( !args->ligate_force )
                    {
                        ir = _get_active_index(args->files);
                        error("Error: The --ligate option is intended for VCFs with perfect overlap, the site %s:%"PRId64" breaks the assumption\n",
                            bcf_seqname(bcf_sr_get_header(args->files,ir),bcf_sr_get_line(args->files,ir)), (int64_t) bcf_sr_get_line(args->files,ir)->pos+1);
                    }
                }

                // Get a line to learn about current position
                ir = _get_active_index(args->files);
                bcf1_t *line = bcf_sr_get_line(args->files,ir);

                // This can happen after bcf_sr_seek: indel may start before the coordinate which we seek to.
                if ( seek_chr>=0 && seek_pos>line->pos && seek_chr==bcf_hdr_name2id(args->out_hdr, bcf_seqname(args->files->readers[ir].header,line)) ) continue;
                seek_pos = seek_chr = -1;

                //  Check if the position overlaps with the next, yet unopened, reader
                int must_seek = 0;
                while ( args->ifname < args->nfnames && args->start_pos[args->ifname]!=-1 && line->pos >= args->start_pos[args->ifname] )
                {
                    must_seek = 1;
                    if ( !bcf_sr_add_reader(args->files,args->fnames[args->ifname]) ) error("Failed to open %s: %s\n", args->fnames[args->ifname],bcf_sr_strerror(args->files->errnum));
                    args->ifname++;
                }
                if ( must_seek )
                {
                    bcf_sr_seek(args->files, bcf_seqname(args->files->readers[ir].header,line), line->pos);
                    seek_pos = line->pos;
                    seek_chr = bcf_hdr_name2id(args->out_hdr, bcf_seqname(args->files->readers[ir].header,line));
                    continue;
                }

                // We are assuming that there is a perfect overlap, sites which are not present in both files are dropped
                if ( args->files->nreaders>1 && !bcf_sr_has_line(args->files,1) && !bcf_sr_region_done(args->files,1) && !args->ligate_force )
                {
                    if ( args->ligate_warn && !site_drop_warned )
                    {
                        ir = _get_active_index(args->files);
                        fprintf(stderr,
                                "Warning: Dropping the site %s:%"PRId64". The --ligate option is intended for VCFs with perfect\n"
                                "         overlap, sites in overlapping regions present in one but missing in other are dropped.\n"
                                "         This warning is printed only once.\n",
                                bcf_seqname(bcf_sr_get_header(args->files,ir),line), (int64_t) line->pos+1);
                        site_drop_warned = 1;
                    }
                    else if ( !args->ligate_warn )
                    {
                        ir = _get_active_index(args->files);
                        error("Error: The --ligate option is intended for VCFs with perfect overlap, the site %s:%"PRId64" breaks the assumption\n",
                            bcf_seqname(bcf_sr_get_header(args->files,ir),bcf_sr_get_line(args->files,ir)), (int64_t) bcf_sr_get_line(args->files,ir)->pos+1);
                    }
                    continue;
                }

                bcf1_t *line0 = bcf_sr_get_line(args->files,0);
                bcf1_t *line1 = args->files->nreaders > 1 ? bcf_sr_get_line(args->files,1) : NULL;
                phased_push(args, line0, line1, is_overlap);
            }

            if ( args->files->nreaders )
            {
                phased_flush(args);
                while ( args->files->nreaders )
                    bcf_sr_remove_reader(args->files, 0);
            }
        }
    }
    else if ( args->files )  // combining overlapping files, using synced reader
    {
        while ( bcf_sr_next_line(args->files) )
        {
            for (i=0; i<args->files->nreaders; i++)
            {
                bcf1_t *line = bcf_sr_get_line(args->files,i);
                if ( !line ) continue;
                bcf_translate(args->out_hdr, args->files->readers[i].header, line);
                if ( bcf_write1(args->out_fh, args->out_hdr, line)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
                if ( args->remove_dups ) break;
            }
        }
    }
    else    // concatenating
    {
        struct timeval t0, t1;
        kstring_t tmp = {0,0,0};
        int prev_chr_id = -1, prev_pos;
        bcf1_t *line = bcf_init();
        for (i=0; i<args->nfnames; i++)
        {
            if ( args->verbose )
            {
                fprintf(stderr,"Concatenating %s", args->fnames[i]);
                gettimeofday(&t0, NULL);
            }
            htsFile *fp = hts_open(args->fnames[i], "r"); if ( !fp ) error("\nFailed to open: %s\n", args->fnames[i]);
            if ( args->n_threads ) hts_set_opt(fp, HTS_OPT_THREAD_POOL, args->tpool);
            bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) error("\nFailed to parse header: %s\n", args->fnames[i]);
            if ( !fp->is_bin && args->output_type&FT_VCF )
            {
                line->max_unpack = BCF_UN_STR;
                // if VCF is on both input and output, avoid VCF to BCF conversion
                while ( hts_getline(fp, KS_SEP_LINE, &fp->line) >=0 )
                {
                    char *str = fp->line.s;
                    while ( *str && *str!='\t' ) str++;
                    tmp.l = 0;
                    kputsn(fp->line.s,str-fp->line.s,&tmp);
                    int chr_id = bcf_hdr_name2id(args->out_hdr, tmp.s);
                    if ( chr_id<0 ) error("\nThe sequence \"%s\" not defined in the header: %s\n(Quick workaround: index the file.)\n", tmp.s, args->fnames[i]);
                    if ( prev_chr_id!=chr_id )
                    {
                        prev_pos = -1;
                        if ( args->seen_seq[chr_id] )
                            error("\nThe chromosome block %s is not contiguous, consider running with -a.\n", tmp.s);
                    }
                    char *end;
                    int pos = strtol(str+1,&end,10) - 1;
                    if ( end==str+1 ) error("Could not parse line: %s\n", fp->line.s);
                    if ( prev_pos > pos )
                        error("\nThe chromosome block %s is not sorted, consider running with -a.\n", tmp.s);
                    args->seen_seq[chr_id] = 1;
                    prev_chr_id = chr_id;

                    if ( vcf_write_line(args->out_fh, &fp->line)!=0 ) error("\nFailed to write %"PRIu64" bytes\n", (uint64_t)fp->line.l);
                }
            }
            else
            {
                // BCF conversion is required
                line->max_unpack = 0;
                while ( bcf_read(fp, hdr, line)==0 )
                {
                    bcf_translate(args->out_hdr, hdr, line);

                    if ( prev_chr_id!=line->rid )
                    {
                        prev_pos = -1;
                        if ( args->seen_seq[line->rid] )
                            error("\nThe chromosome block %s is not contiguous, consider running with -a.\n", bcf_seqname(args->out_hdr, line));
                    }
                    if ( prev_pos > line->pos )
                        error("\nThe chromosome block %s is not sorted, consider running with -a.\n", bcf_seqname(args->out_hdr, line));
                    args->seen_seq[line->rid] = 1;
                    prev_chr_id = line->rid;

                    if ( bcf_write(args->out_fh, args->out_hdr, line)!=0 ) error("\nFailed to write\n");
                }
            }
            bcf_hdr_destroy(hdr);
            hts_close(fp);
            if ( args->verbose )
            {
                gettimeofday(&t1, NULL);
                double delta = (t1.tv_sec - t0.tv_sec) * 1e6 + (t1.tv_usec - t0.tv_usec);
                fprintf(stderr,"\t%f seconds\n",delta/1e6);
            }
        }
        bcf_destroy(line);
        free(tmp.s);
    }
}

int print_vcf_gz_header(BGZF *fp, BGZF *bgzf_out, int print_header, kstring_t *tmp)
{
    char *buffer = (char*) fp->uncompressed_block;

    // Read the header and find the position of the data block
    if ( buffer[0]!='#' ) error("Could not parse the header, expected '#', found '%c'\n", buffer[0]);

    int nskip = 1;     // end of the header in the current uncompressed block
    while (1)
    {
        if ( buffer[nskip]=='\n' )
        {
            nskip++;
            if ( nskip>=fp->block_length )
            {
                kputsn(buffer,nskip,tmp);
                if ( bgzf_read_block(fp) != 0 ) return -1;
                if ( !fp->block_length ) break;
                nskip = 0;
            }
            // The header has finished
            if ( buffer[nskip]!='#' )
            {
                kputsn(buffer,nskip,tmp);
                break;
            }
        }
        nskip++;
        if ( nskip>=fp->block_length )
        {
            kputsn(buffer,fp->block_length,tmp);
            if ( bgzf_read_block(fp) != 0 ) return -1;
            if ( !fp->block_length ) break;
            nskip = 0;
        }
    }
    if ( print_header )
    {
        if ( bgzf_write(bgzf_out,tmp->s,tmp->l) != tmp->l ) error("Failed to write %"PRIu64" bytes\n", (uint64_t)tmp->l);
        tmp->l = 0;
    }
    return nskip;
}

static inline int unpackInt16(const uint8_t *buffer)
{
    return buffer[0] | buffer[1] << 8;
}
static int check_header(const uint8_t *header)
{
    if ( header[0] != 31 || header[1] != 139 || header[2] != 8 ) return -2;
    return ((header[3] & 4) != 0
            && unpackInt16((uint8_t*)&header[10]) == 6
            && header[12] == 'B' && header[13] == 'C'
            && unpackInt16((uint8_t*)&header[14]) == 2) ? 0 : -1;
}
static void _check_hrecs(const bcf_hdr_t *hdr0, const bcf_hdr_t *hdr, char *fname0, char *fname)
{
    int j;
    for (j=0; j<hdr0->nhrec; j++)
    {
        bcf_hrec_t *hrec0 = hdr0->hrec[j];
        if ( hrec0->type!=BCF_HL_FLT && hrec0->type!=BCF_HL_INFO && hrec0->type!=BCF_HL_FMT && hrec0->type!=BCF_HL_CTG ) continue;    // skip fiels w/o IDX
        int itag = bcf_hrec_find_key(hrec0, "ID");
        bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, hrec0->type, "ID", hrec0->vals[itag], NULL);

        char *type = NULL;
        if ( hrec0->type==BCF_HL_FLT ) type = "FILTER";
        if ( hrec0->type==BCF_HL_INFO ) type = "INFO";
        if ( hrec0->type==BCF_HL_FMT ) type = "FORMAT";
        if ( hrec0->type==BCF_HL_CTG ) type = "contig";

        if ( !hrec )
            error("Cannot use --naive, incompatible headers, the tag %s/%s not present in %s\n",type,hrec0->vals[itag],fname);

        int idx0 = bcf_hrec_find_key(hrec0, "IDX");
        int idx  = bcf_hrec_find_key(hrec,  "IDX");
        if ( idx0<0 || idx<0 )
            error("fixme: unexpected IDX<0 for %s/%s in %s or %s\n",type,hrec0->vals[itag],fname0,fname);
        if ( strcmp(hrec0->vals[idx0],hrec->vals[idx]) )
            error("Cannot use --naive, use --naive-force instead: different order the tag %s/%s in %s vs %s\n",type,hrec0->vals[itag],fname0,fname);
    }
}
static void naive_concat_check_headers(args_t *args)
{
    fprintf(stderr,"Checking the headers of %d files.\n",args->nfnames);
    bcf_hdr_t *hdr0 = NULL;
    int i,j;
    for (i=0; i<args->nfnames; i++)
    {
        htsFile *fp = hts_open(args->fnames[i], "r"); if ( !fp ) error("Failed to open: %s\n", args->fnames[i]);
        bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) error("Failed to parse header: %s\n", args->fnames[i]);
        htsFormat type = *hts_get_format(fp);
        hts_close(fp);

        if ( i==0 )
        {
            hdr0 = hdr;
            continue;
        }

        // check the samples
        if ( bcf_hdr_nsamples(hdr0)!=bcf_hdr_nsamples(hdr) )
            error("Cannot concatenate, different number of samples: %d vs %d in %s vs %s\n",bcf_hdr_nsamples(hdr0),bcf_hdr_nsamples(hdr),args->fnames[0],args->fnames[i]);
        for (j=0; j<bcf_hdr_nsamples(hdr0); j++)
            if ( strcmp(hdr0->samples[j],hdr->samples[j]) )
                error("Cannot concatenate, different samples in %s vs %s\n",args->fnames[0],args->fnames[i]);

        // if BCF, check if tag IDs are consistent in the dictionary of strings
        if ( type.compression!=bgzf )
            error("The --naive option works only for compressed BCFs or VCFs, sorry :-/\n");

        _check_hrecs(hdr0,hdr,args->fnames[0],args->fnames[i]);
        _check_hrecs(hdr,hdr0,args->fnames[i],args->fnames[0]);

        bcf_hdr_destroy(hdr);
    }
    if ( hdr0 ) bcf_hdr_destroy(hdr0);
    fprintf(stderr,"Done, the headers are compatible.\n");
}
static void naive_concat(args_t *args)
{
    if ( !args->naive_concat_trust_headers )
        naive_concat_check_headers(args);

    // only compressed BCF atm
    BGZF *bgzf_out = bgzf_open(args->output_fname,"w");;

    htsFormat output_type;
    output_type.format = (args->output_type & FT_VCF) ? vcf : bcf;
    output_type.compression = (args->output_type & FT_GZ) ? bgzf : no_compression;

    struct timeval t0, t1;
    const size_t page_size = BGZF_MAX_BLOCK_SIZE;
    uint8_t *buf = (uint8_t*) malloc(page_size);
    kstring_t tmp = {0,0,0};
    int i, file_types = 0;
    for (i=0; i<args->nfnames; i++)
    {
        if ( args->verbose )
        {
            fprintf(stderr,"Concatenating %s", args->fnames[i]);
            gettimeofday(&t0, NULL);
        }
        htsFile *hts_fp = hts_open(args->fnames[i],"r");
        if ( !hts_fp ) error("\nFailed to open: %s\n", args->fnames[i]);
        htsFormat type = *hts_get_format(hts_fp);

        if ( type.compression!=bgzf )
            error("\nThe --naive option works only for compressed BCFs or VCFs\n");
        file_types |= type.format==vcf ? 1 : 2;
        if ( file_types==3 )
            error("\nThe --naive option works only for compressed files of the same type, all BCFs or all VCFs\n");
        if ( args->explicit_output_type )
        {
            if ( output_type.format!=type.format )
                error("\nThe --naive option works only for the output of the same type, all BCFs or all VCFs\n");
            if ( output_type.compression!=type.compression )
                error("\nThe --naive option works only for the output of the same compression type\n");
        }

        BGZF *fp = hts_get_bgzfp(hts_fp);
        if ( !fp || bgzf_read_block(fp) != 0 || !fp->block_length )
            error("\nFailed to read %s: %s\n", args->fnames[i], strerror(errno));

        int nskip;
        if ( type.format==bcf )
        {
            uint8_t magic[5];
            if ( bgzf_read(fp, magic, 5) != 5 ) error("\nFailed to read the BCF header in %s\n", args->fnames[i]);
            if (strncmp((char*)magic, "BCF\2\2", 5) != 0) error("\nInvalid BCF magic string in %s\n", args->fnames[i]);

            if ( bgzf_read(fp, &tmp.l, 4) != 4 ) error("\nFailed to read the BCF header in %s\n", args->fnames[i]);
            hts_expand(char,tmp.l,tmp.m,tmp.s);
            if ( bgzf_read(fp, tmp.s, tmp.l) != tmp.l ) error("\nFailed to read the BCF header in %s\n", args->fnames[i]);

            // write only the first header
            if ( i==0 )
            {
                if ( bgzf_write(bgzf_out, "BCF\2\2", 5) !=5 ) error("\nFailed to write %d bytes to %s\n", 5,args->output_fname);
                if ( bgzf_write(bgzf_out, &tmp.l, 4) !=4 ) error("\nFailed to write %d bytes to %s\n", 4,args->output_fname);
                if ( bgzf_write(bgzf_out, tmp.s, tmp.l) != tmp.l) error("\nFailed to write %"PRId64" bytes to %s\n", (uint64_t)tmp.l,args->output_fname);
            }
            nskip = fp->block_offset;
        }
        else
        {
            nskip = print_vcf_gz_header(fp, bgzf_out, i==0?1:0, &tmp);
            if ( nskip==-1 ) error("\nError reading %s\n", args->fnames[i]);
        }

        // Output all non-header data that were read together with the header block
        if ( fp->block_length - nskip > 0 )
        {
            if ( bgzf_write(bgzf_out, (char *)fp->uncompressed_block+nskip, fp->block_length-nskip)<0 ) error("\nError: %d\n",fp->errcode);
        }
        if ( bgzf_flush(bgzf_out)<0 ) error("\nError: %d\n",bgzf_out->errcode);


        // Stream the rest of the file as it is, without recompressing, but remove BGZF EOF blocks
        // The final bgzf eof block will be added by bgzf_close.
        ssize_t nread, nblock, nwr;
        const int nheader = 18, neof = 28;
        const uint8_t *eof = (uint8_t*) "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0";
        while (1)
        {
            nread = bgzf_raw_read(fp, buf, nheader);
            if ( !nread ) break;
            if ( nread != nheader || check_header(buf)!=0 ) error("\nCould not parse the header of a bgzf block: %s\n",args->fnames[i]);
            nblock = unpackInt16(buf+16) + 1;
            assert( nblock <= page_size && nblock >= nheader );
            nread += bgzf_raw_read(fp, buf+nheader, nblock - nheader);
            if ( nread!=nblock ) error("\nCould not read %"PRId64" bytes: %s\n",(uint64_t)nblock,args->fnames[i]);
            if ( nread==neof && !memcmp(buf,eof,neof) ) continue;
            nwr = bgzf_raw_write(bgzf_out, buf, nread);
            if ( nwr != nread ) error("\nWrite failed, wrote %"PRId64" instead of %d bytes.\n", (uint64_t)nwr,(int)nread);
        }
        if (hts_close(hts_fp)) error("\nClose failed: %s\n",args->fnames[i]);
        if ( args->verbose )
        {
            gettimeofday(&t1, NULL);
            double delta = (t1.tv_sec - t0.tv_sec) * 1e6 + (t1.tv_usec - t0.tv_usec);
            fprintf(stderr,"\t%f seconds\n",delta/1e6);
        }
    }
    free(buf);
    free(tmp.s);
    if (bgzf_close(bgzf_out) < 0) error("Error: %d\n",bgzf_out->errcode);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Concatenate or combine VCF/BCF files. All source files must have the same sample\n");
    fprintf(stderr, "         columns appearing in the same order. The program can be used, for example, to\n");
    fprintf(stderr, "         concatenate chromosome VCFs into one VCF, or combine a SNP VCF and an indel\n");
    fprintf(stderr, "         VCF into one. The input files must be sorted by chr and position. The files\n");
    fprintf(stderr, "         must be given in the correct order to produce sorted VCF on output unless\n");
    fprintf(stderr, "         the -a, --allow-overlaps option is specified. With the --naive option, the files\n");
    fprintf(stderr, "         are concatenated without being recompressed, which is very fast.\n");
    fprintf(stderr, "Usage:   bcftools concat [options] <A.vcf.gz> [<B.vcf.gz> [...]]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -a, --allow-overlaps           First coordinate of the next file can precede last record of the current file.\n");
    fprintf(stderr, "   -c, --compact-PS               Do not output PS tag at each site, only at the start of a new phase set block.\n");
    fprintf(stderr, "   -d, --rm-dups STRING           Output duplicate records present in multiple files only once: <snps|indels|both|all|exact>\n");
    fprintf(stderr, "   -D, --remove-duplicates        Alias for -d exact\n");
    fprintf(stderr, "   -f, --file-list FILE           Read the list of files from a file.\n");
    fprintf(stderr, "   -l, --ligate                   Ligate phased VCFs by matching phase at overlapping haplotypes\n");
    fprintf(stderr, "       --ligate-force             Ligate even non-overlapping chunks, keep all sites\n");
    fprintf(stderr, "       --ligate-warn              Drop sites in imperfect overlaps\n");
    fprintf(stderr, "       --no-version               Do not append version and command line to the header\n");
    fprintf(stderr, "   -n, --naive                    Concatenate files without recompression, a header check compatibility is performed\n");
    fprintf(stderr, "       --naive-force              Same as --naive, but header compatibility is not checked. Dangerous, use with caution.\n");
    fprintf(stderr, "   -o, --output FILE              Write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n");
    fprintf(stderr, "   -q, --min-PQ INT               Break phase set if phasing quality is lower than <int> [30]\n");
    fprintf(stderr, "   -r, --regions REGION           Restrict to comma-separated list of regions\n");
    fprintf(stderr, "   -R, --regions-file FILE        Restrict to regions listed in a file\n");
    fprintf(stderr, "       --regions-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n");
    fprintf(stderr, "       --threads INT              Use multithreading with <int> worker threads [0]\n");
    fprintf(stderr, "   -v, --verbose 0|1              Set verbosity level [1]\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfconcat(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->output_fname = "-";
    args->output_type = FT_VCF;
    args->n_threads = 0;
    args->record_cmd_line = 1;
    args->min_PQ  = 30;
    args->verbose = 1;
    args->clevel = -1;

    static struct option loptions[] =
    {
        {"verbose",required_argument,NULL,'v'},
        {"naive",no_argument,NULL,'n'},
        {"naive-force",no_argument,NULL,7},
        {"compact-PS",no_argument,NULL,'c'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"regions-overlap",required_argument,NULL,12},
        {"remove-duplicates",no_argument,NULL,'D'},
        {"rm-dups",required_argument,NULL,'d'},
        {"allow-overlaps",no_argument,NULL,'a'},
        {"ligate",no_argument,NULL,'l'},
        {"ligate-force",no_argument,NULL,10},
        {"ligate-warn",no_argument,NULL,11},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"threads",required_argument,NULL,9},
        {"file-list",required_argument,NULL,'f'},
        {"min-PQ",required_argument,NULL,'q'},
        {"no-version",no_argument,NULL,8},
        {NULL,0,NULL,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "h:?o:O:f:alq:Dd:r:R:cnv:",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'c': args->compact_PS = 1; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; args->regions_is_file = 1; break;
            case 'd': args->remove_dups = optarg; break;
            case 'D': args->remove_dups = "exact"; break;
            case 'q': 
                args->min_PQ = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --min-PQ %s\n", optarg);
                break;
            case 'n': args->naive_concat = 1; break;
            case 'a': args->allow_overlaps = 1; break;
            case 'l': args->phased_concat = 1; break;
            case 'f': args->file_list = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                args->explicit_output_type = 1;
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
            case 10 : args->ligate_force = 1; break;
            case 11 : args->ligate_warn  = 1; break;
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case  8 : args->record_cmd_line = 0; break;
            case  7 : args->naive_concat = 1; args->naive_concat_trust_headers = 1; break;
            case 12 :
                args->regions_overlap = parse_overlap_option(optarg);
                if ( args->regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case 'v':
                      args->verbose = strtol(optarg, &tmp, 0);
                      if ( *tmp || args->verbose<0 || args->verbose>1 ) error("Error: currently only --verbose 0 or --verbose 1 is supported\n");
                      break;
            case 'h':
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    while ( optind<argc )
    {
        args->nfnames++;
        args->fnames = (char **)realloc(args->fnames,sizeof(char*)*args->nfnames);
        args->fnames[args->nfnames-1] = strdup(argv[optind]);
        optind++;
    }
    if ( args->ligate_force && args->ligate_warn ) error("The options cannot be combined: --ligate-force and --ligate-warn\n");
    if ( args->allow_overlaps && args->phased_concat ) error("The options -a and -l should not be combined. Please run with -l only.\n");
    if ( args->compact_PS && !args->phased_concat ) error("The -c option is intended only with -l\n");
    if ( args->file_list )
    {
        if ( args->nfnames ) error("Cannot combine -l with file names on command line.\n");
        args->fnames = hts_readlines(args->file_list, &args->nfnames);
        if ( !args->fnames ) error("Could not read the file: %s\n", args->file_list);
    }
    if ( !args->nfnames ) usage(args);
    if ( args->remove_dups && !args->allow_overlaps ) error("The -D option is supported only with -a\n");
    if ( args->regions_list && !args->allow_overlaps ) error("The -r/-R option is supported only with -a\n");
    if ( args->naive_concat )
    {
        if ( args->allow_overlaps ) error("The option --naive cannot be combined with --allow-overlaps\n");
        if ( args->phased_concat ) error("The option --naive cannot be combined with --ligate\n");
        naive_concat(args);
        destroy_data(args);
        free(args);
        return 0;
    }
    init_data(args);
    concat(args);
    destroy_data(args);
    free(args);
    return 0;
}
