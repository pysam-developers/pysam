/*  vcfconcat.c -- Concatenate or combine VCF/BCF files.

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

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kseq.h>
#include <htslib/bgzf.h>
#include <htslib/tbx.h> // for hts_get_bgzfp()
#include "bcftools.h"

typedef struct _args_t
{
    bcf_srs_t *files;
    htsFile *out_fh;
    int output_type, n_threads, record_cmd_line;
    bcf_hdr_t *out_hdr;
    int *seen_seq;

    // phasing
    int *start_pos, start_tid, ifname;
    int *swap_phase, nswap, *nmatch, *nmism;
    bcf1_t **buf;
    int nbuf, mbuf, prev_chr, min_PQ, prev_pos_check;
    int32_t *GTa, *GTb, mGTa, mGTb, *phase_qual, *phase_set;

    char **argv, *output_fname, *file_list, **fnames, *remove_dups, *regions_list;
    int argc, nfnames, allow_overlaps, phased_concat, regions_is_file;
    int compact_PS, phase_set_changed, naive_concat;
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
        hts_close(fp);
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
    args->out_fh = hts_open(args->output_fname,hts_bcf_wmode(args->output_type));
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    if ( args->n_threads ) hts_set_threads(args->out_fh, args->n_threads);

    bcf_hdr_write(args->out_fh, args->out_hdr);

    if ( args->allow_overlaps )
    {
        args->files = bcf_sr_init();
        args->files->require_index = 1;
        if ( args->regions_list )
        {
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
        args->files = bcf_sr_init();
        args->files->require_index = 1;
        args->ifname = 0;
    }
}

static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->nfnames; i++) free(args->fnames[i]);
    free(args->fnames);
    if ( args->files ) bcf_sr_destroy(args->files);
    if ( args->out_fh )
    {
        if ( hts_close(args->out_fh)!=0 ) error("hts_close error\n");
    }
    if ( args->out_hdr ) bcf_hdr_destroy(args->out_hdr);
    free(args->seen_seq);
    free(args->start_pos);
    free(args->swap_phase);
    for (i=0; i<args->mbuf; i++) bcf_destroy(args->buf[i]);
    free(args->buf);
    free(args->GTa);
    free(args->GTb);
    free(args->nmatch);
    free(args->nmism);
    free(args->phase_qual);
    free(args->phase_set);
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
        bcf1_t *arec = args->buf[i];
        bcf1_t *brec = args->buf[i+1];

        int nGTs = bcf_get_genotypes(ahdr, arec, &args->GTa, &args->mGTa);
        if ( nGTs < 0 ) 
        {
            if ( !gt_absent_warned )
            {
                fprintf(stderr,"GT is not present at %s:%d. (This warning is printed only once.)\n", bcf_seqname(ahdr,arec), arec->pos+1);
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
                fprintf(stderr,"GT is not present at %s:%d. (This warning is printed only once.)\n", bcf_seqname(bhdr,brec), brec->pos+1);
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
        bcf1_t *arec = args->buf[i];
        bcf_translate(args->out_hdr, args->files->readers[0].header, arec);
        if ( args->nswap )
            phase_update(args, args->out_hdr, arec);
        if ( !args->compact_PS || args->phase_set_changed )
        {
            bcf_update_format_int32(args->out_hdr,arec,"PS",args->phase_set,nsmpl);
            args->phase_set_changed = 0;
        }
        bcf_write(args->out_fh, args->out_hdr, arec);

        if ( arec->pos < args->prev_pos_check ) error("FIXME, disorder: %s:%d vs %d  [1]\n", bcf_seqname(args->files->readers[0].header,arec),arec->pos+1,args->prev_pos_check+1);
        args->prev_pos_check = arec->pos;
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
        bcf1_t *brec = args->buf[i+1];
        bcf_translate(args->out_hdr, args->files->readers[1].header, brec);
        if ( !PQ_printed )
        {
            bcf_update_format_int32(args->out_hdr,brec,"PQ",args->phase_qual,nsmpl);
            PQ_printed = 1;
            for (j=0; j<nsmpl; j++)
                if ( args->phase_qual[j] < args->min_PQ ) 
                {
                    args->phase_set[j] = brec->pos+1;
                    args->phase_set_changed = 1;
                }
                else if ( args->compact_PS ) args->phase_set[j] = bcf_int32_missing;
        }
        if ( args->nswap )
            phase_update(args, args->out_hdr, brec);
        if ( !args->compact_PS || args->phase_set_changed )
        {
            bcf_update_format_int32(args->out_hdr,brec,"PS",args->phase_set,nsmpl);
            args->phase_set_changed = 0;
        }
        bcf_write(args->out_fh, args->out_hdr, brec);

        if ( brec->pos < args->prev_pos_check ) error("FIXME, disorder: %s:%d vs %d  [2]\n", bcf_seqname(args->files->readers[1].header,brec),brec->pos+1,args->prev_pos_check+1);
        args->prev_pos_check = brec->pos;
    }
    args->nbuf = 0;
}

static void phased_push(args_t *args, bcf1_t *arec, bcf1_t *brec)
{
    if ( arec && arec->errcode )
        error("Parse error at %s:%d, cannot proceed: %s\n", bcf_seqname(args->files->readers[0].header,arec),arec->pos+1, args->files->readers[0].fname);
    if ( brec && brec->errcode )
        error("Parse error at %s:%d, cannot proceed: %s\n", bcf_seqname(args->files->readers[1].header,brec),brec->pos+1, args->files->readers[1].fname);

    int i, nsmpl = bcf_hdr_nsamples(args->out_hdr);
    int chr_id = bcf_hdr_name2id(args->out_hdr, bcf_seqname(args->files->readers[0].header,arec));
    if ( args->prev_chr<0 || args->prev_chr!=chr_id )
    {
        if ( args->prev_chr>=0 ) phased_flush(args);

        for (i=0; i<nsmpl; i++)
            args->phase_set[i] = arec->pos+1;
        args->phase_set_changed = 1;

        if ( args->seen_seq[chr_id] ) error("The chromosome block %s is not contiguous\n", bcf_seqname(args->files->readers[0].header,arec));
        args->seen_seq[chr_id] = 1;
        args->prev_chr = chr_id;
        args->prev_pos_check = -1;
    }

    if ( !brec )
    {
        bcf_translate(args->out_hdr, args->files->readers[0].header, arec);
        if ( args->nswap )
            phase_update(args, args->out_hdr, arec);
        if ( !args->compact_PS || args->phase_set_changed )
        {
            bcf_update_format_int32(args->out_hdr,arec,"PS",args->phase_set,nsmpl);
            args->phase_set_changed = 0;
        }
        bcf_write(args->out_fh, args->out_hdr, arec);

        if ( arec->pos < args->prev_pos_check )
            error("FIXME, disorder: %s:%d in %s vs %d written  [3]\n", bcf_seqname(args->files->readers[0].header,arec), arec->pos+1,args->files->readers[0].fname, args->prev_pos_check+1);
        args->prev_pos_check = arec->pos;
        return;
    }

    int m = args->mbuf;
    args->nbuf += 2;
    hts_expand(bcf1_t*,args->nbuf,args->mbuf,args->buf);
    for (i=m; i<args->mbuf; i++)
        args->buf[i] = bcf_init1();

    SWAP(bcf1_t*, args->files->readers[0].buffer[0], args->buf[args->nbuf-2]);
    SWAP(bcf1_t*, args->files->readers[1].buffer[0], args->buf[args->nbuf-1]);
}

static void concat(args_t *args)
{
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

            int nret;
            while ( (nret = bcf_sr_next_line(args->files)) )
            {
                if ( !bcf_sr_has_line(args->files,0) )  // no input from the first reader
                {
                    // We are assuming that there is a perfect overlap, sites which are not present in both files are dropped
                    if ( ! bcf_sr_region_done(args->files,0) ) continue;

                    phased_flush(args);
                    bcf_sr_remove_reader(args->files, 0);
                }

                // Get a line to learn about current position
                for (i=0; i<args->files->nreaders; i++)
                    if ( bcf_sr_has_line(args->files,i) ) break;
                bcf1_t *line = bcf_sr_get_line(args->files,i);

                // This can happen after bcf_sr_seek: indel may start before the coordinate which we seek to.
                if ( seek_chr>=0 && seek_pos>line->pos && seek_chr==bcf_hdr_name2id(args->out_hdr, bcf_seqname(args->files->readers[i].header,line)) ) continue;
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
                    bcf_sr_seek(args->files, bcf_seqname(args->files->readers[i].header,line), line->pos);
                    seek_pos = line->pos;
                    seek_chr = bcf_hdr_name2id(args->out_hdr, bcf_seqname(args->files->readers[i].header,line));
                    continue;
                }

                // We are assuming that there is a perfect overlap, sites which are not present in both files are dropped
                if ( args->files->nreaders>1 && !bcf_sr_has_line(args->files,1) && !bcf_sr_region_done(args->files,1) ) continue;

                phased_push(args, bcf_sr_get_line(args->files,0), args->files->nreaders>1 ? bcf_sr_get_line(args->files,1) : NULL);
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
                bcf_write1(args->out_fh, args->out_hdr, line);
                if ( args->remove_dups ) break;
            }
        }
    }
    else    // concatenating
    {
        kstring_t tmp = {0,0,0};
        int prev_chr_id = -1, prev_pos;
        bcf1_t *line = bcf_init();
        for (i=0; i<args->nfnames; i++)
        {
            htsFile *fp = hts_open(args->fnames[i], "r"); if ( !fp ) error("Failed to open: %s\n", args->fnames[i]);
            bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) error("Failed to parse header: %s\n", args->fnames[i]);
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
                    if ( chr_id<0 ) error("The sequence \"%s\" not defined in the header: %s\n(Quick workaround: index the file.)\n", tmp.s, args->fnames[i]);
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
                        error("The chromosome block %s is not sorted, consider running with -a.\n", tmp.s);
                    args->seen_seq[chr_id] = 1;
                    prev_chr_id = chr_id;

                    if ( vcf_write_line(args->out_fh, &fp->line)!=0 ) error("Failed to write %d bytes\n", fp->line.l);
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
                        error("The chromosome block %s is not sorted, consider running with -a.\n", bcf_seqname(args->out_hdr, line));
                    args->seen_seq[line->rid] = 1;
                    prev_chr_id = line->rid;

                    if ( bcf_write(args->out_fh, args->out_hdr, line)!=0 ) error("Failed to write\n");
                }
            }
            bcf_hdr_destroy(hdr);
            hts_close(fp);
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
        if ( bgzf_write(bgzf_out,tmp->s,tmp->l) != tmp->l ) error("Failed to write %d bytes\n", tmp->l);
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
static void naive_concat(args_t *args)
{
    // only compressed BCF atm
    BGZF *bgzf_out = bgzf_open(args->output_fname,"w");;

    const size_t page_size = BGZF_MAX_BLOCK_SIZE;
    uint8_t *buf = (uint8_t*) malloc(page_size);
    kstring_t tmp = {0,0,0};
    int i, file_types = 0;
    for (i=0; i<args->nfnames; i++)
    {
        htsFile *hts_fp = hts_open(args->fnames[i],"r");
        if ( !hts_fp ) error("Failed to open: %s\n", args->fnames[i]);
        htsFormat type = *hts_get_format(hts_fp);

        if ( type.compression!=bgzf )
            error("The --naive option works only for compressed BCFs or VCFs, sorry :-/\n");
        file_types |= type.format==vcf ? 1 : 2;
        if ( file_types==3 )
            error("The --naive option works only for compressed files of the same type, all BCFs or all VCFs :-/\n");

        BGZF *fp = hts_get_bgzfp(hts_fp);
        if ( !fp || bgzf_read_block(fp) != 0 || !fp->block_length )
            error("Failed to read %s: %s\n", args->fnames[i], strerror(errno));

        int nskip;
        if ( type.format==bcf )
        {
            uint8_t magic[5];
            if ( bgzf_read(fp, magic, 5) != 5 ) error("Failed to read the BCF header in %s\n", args->fnames[i]);
            if (strncmp((char*)magic, "BCF\2\2", 5) != 0) error("Invalid BCF magic string in %s\n", args->fnames[i]);

            if ( bgzf_read(fp, &tmp.l, 4) != 4 ) error("Failed to read the BCF header in %s\n", args->fnames[i]);
            hts_expand(char,tmp.l,tmp.m,tmp.s);
            if ( bgzf_read(fp, tmp.s, tmp.l) != tmp.l ) error("Failed to read the BCF header in %s\n", args->fnames[i]);

            // write only the first header
            if ( i==0 )
            {
                if ( bgzf_write(bgzf_out, "BCF\2\2", 5) !=5 ) error("Failed to write %d bytes to %s\n", 5,args->output_fname);
                if ( bgzf_write(bgzf_out, &tmp.l, 4) !=4 ) error("Failed to write %d bytes to %s\n", 4,args->output_fname);
                if ( bgzf_write(bgzf_out, tmp.s, tmp.l) != tmp.l) error("Failed to write %d bytes to %s\n", tmp.l,args->output_fname);
            }
            nskip = fp->block_offset;
        }
        else
        {
            nskip = print_vcf_gz_header(fp, bgzf_out, i==0?1:0, &tmp);
            if ( nskip==-1 ) error("Error reading %s\n", args->fnames[i]);
        }

        // Output all non-header data that were read together with the header block
        if ( fp->block_length - nskip > 0 )
        {
            if ( bgzf_write(bgzf_out, (char *)fp->uncompressed_block+nskip, fp->block_length-nskip)<0 ) error("Error: %d\n",fp->errcode);
        }
        if ( bgzf_flush(bgzf_out)<0 ) error("Error: %d\n",bgzf_out->errcode);


        // Stream the rest of the file as it is, without recompressing, but remove BGZF EOF blocks
        // The final bgzf eof block will be added by bgzf_close.
        ssize_t nread, nblock, nwr;
        const int nheader = 18, neof = 28;
        const uint8_t *eof = (uint8_t*) "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0";
        while (1)
        {
            nread = bgzf_raw_read(fp, buf, nheader);
            if ( !nread ) break;
            if ( nread != nheader || check_header(buf)!=0 ) error("Could not parse the header of a bgzf block: %s\n",args->fnames[i]);
            nblock = unpackInt16(buf+16) + 1;
            assert( nblock <= page_size && nblock >= nheader );
            nread += bgzf_raw_read(fp, buf+nheader, nblock - nheader);
            if ( nread!=nblock ) error("Could not read %d bytes: %s\n",nblock,args->fnames[i]);
            if ( nread==neof && !memcmp(buf,eof,neof) ) continue;
            nwr = bgzf_raw_write(bgzf_out, buf, nread);
            if ( nwr != nread ) error("Write failed, wrote %d instead of %d bytes.\n", nwr,(int)nread);
        }
        if (hts_close(hts_fp)) error("Close failed: %s\n",args->fnames[i]);
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
    fprintf(stderr, "         are concatenated without being recompressed, which is very fast but dangerous\n");
    fprintf(stderr, "         if the BCF headers differ.\n");
    fprintf(stderr, "Usage:   bcftools concat [options] <A.vcf.gz> [<B.vcf.gz> [...]]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -a, --allow-overlaps           First coordinate of the next file can precede last record of the current file.\n");
    fprintf(stderr, "   -c, --compact-PS               Do not output PS tag at each site, only at the start of a new phase set block.\n");
    fprintf(stderr, "   -d, --rm-dups <string>         Output duplicate records present in multiple files only once: <snps|indels|both|all|none>\n");
    fprintf(stderr, "   -D, --remove-duplicates        Alias for -d none\n");
    fprintf(stderr, "   -f, --file-list <file>         Read the list of files from a file.\n");
    fprintf(stderr, "   -l, --ligate                   Ligate phased VCFs by matching phase at overlapping haplotypes\n");
    fprintf(stderr, "       --no-version               Do not append version and command line to the header\n");
    fprintf(stderr, "   -n, --naive                    Concatenate files without recompression (dangerous, use with caution)\n");
    fprintf(stderr, "   -o, --output <file>            Write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "   -q, --min-PQ <int>             Break phase set if phasing quality is lower than <int> [30]\n");
    fprintf(stderr, "   -r, --regions <region>         Restrict to comma-separated list of regions\n");
    fprintf(stderr, "   -R, --regions-file <file>      Restrict to regions listed in a file\n");
    fprintf(stderr, "       --threads <int>            Number of extra output compression threads [0]\n");
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

    static struct option loptions[] =
    {
        {"naive",no_argument,NULL,'n'},
        {"compact-PS",no_argument,NULL,'c'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"remove-duplicates",no_argument,NULL,'D'},
        {"rm-dups",required_argument,NULL,'d'},
        {"allow-overlaps",no_argument,NULL,'a'},
        {"ligate",no_argument,NULL,'l'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"threads",required_argument,NULL,9},
        {"file-list",required_argument,NULL,'f'},
        {"min-PQ",required_argument,NULL,'q'},
        {"no-version",no_argument,NULL,8},
        {NULL,0,NULL,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "h:?o:O:f:alq:Dd:r:R:cn",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'c': args->compact_PS = 1; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; args->regions_is_file = 1; break;
            case 'd': args->remove_dups = optarg; break;
            case 'D': args->remove_dups = "none"; break;
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
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                };
                break;
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case  8 : args->record_cmd_line = 0; break;
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
    if ( args->allow_overlaps && args->phased_concat ) args->allow_overlaps = 0;
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
