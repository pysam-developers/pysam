/*  reheader.c -- reheader subcommand.

    Copyright (C) 2014-2022 Genome Research Ltd.

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
#include <inttypes.h>
#include <fcntl.h>
#include <math.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/tbx.h> // for hts_get_bgzfp()
#include <htslib/kseq.h>
#include <htslib/thread_pool.h>
#include <htslib/faidx.h>
#include <htslib/khash_str2int.h>
#include "bcftools.h"
#include "khash_str2str.h"

typedef struct _args_t
{
    char **argv, *fname, *samples_fname, *header_fname, *output_fname;
    char *fai_fname, *rm_tmpfile, *tmp_prefix;
    htsFile *fp;
    htsFormat type;
    htsThreadPool *threads;
    int argc, n_threads;
}
args_t;

static inline int is_escaped(const char *min, const char *str)
{
    int n = 0;
    while ( --str>=min && *str=='\\' ) n++;
    return n%2;
}
static char *copy_and_update_contig_line(faidx_t *fai, char *line, void *chr_seen, kstring_t *dst)
{
    kstring_t key = {0,0,0}, val = {0,0,0}, tmp = {0,0,0};
    char *chr_name = NULL, *p, *q = line + 9;   // skip ##contig=
    char *end = q;
    int nopen = 1, chr_len = 0;
    while ( *end && *end!='\n' ) end++;
    while ( *q && *q!='\n' && nopen>0 )
    {
        p = ++q;
        while ( *q && (*q==' ' || *q=='\t') ) { p++; q++; }
        // ^[A-Za-z_][0-9A-Za-z_.]*$
        if (p==q && *q && (isalpha(*q) || *q=='_'))
        {
            q++;
            while ( *q && (isalnum(*q) || *q=='_' || *q=='.') ) q++;
        }
        int n = q-p;
        int m = 0;
        while ( *q && (*q==' ' || *q=='\t') ) { q++; m++; }
        if ( *q!='=' || !n )
        {
            char *x = q;
            while ( *x && *x!='\n' ) x++;
            *x = '\0';
            error("Could not parse the line: %s [%s][%s]\n", line,p,q);
        }
        key.l = 0;
        kputsn(p,q-p-m,&key);
        p = ++q;
        while ( *q && (*q==' ' || *q=='\t') ) { p++; q++; }
        int quoted = *p=='"' ? 1 : 0;
        if ( quoted ) p++, q++;
        while ( *q && *q != '\n' )
        {
            if ( quoted ) { if ( *q=='"' && !is_escaped(p,q) ) break; }
            else
            {
                if ( *q=='<' ) nopen++;
                if ( *q=='>' ) nopen--;
                if ( !nopen ) break;
                if ( *q==',' && nopen==1 ) break;
            }
            q++;
        }
        char *r = q;
        while ( r > p && r[-1] == ' ' ) r--;
        val.l = 0;
        kputsn(p,r-p,&val);
        if ( quoted && *q=='"' ) q++;
        if ( *q=='>' ) { nopen--; q++; }
        if ( !strcmp("length",key.s) ) continue;
        if ( !strcmp("ID",key.s) )
        {
            if ( khash_str2int_has_key(chr_seen,val.s) ) continue;
            chr_len = faidx_seq_len(fai, val.s);
            if ( chr_len==-1 )
            {
                free(val.s); free(key.s); free(tmp.s);
                return end;   // the sequence is not in fai, remove
            }
            chr_name = strdup(val.s);
            khash_str2int_inc(chr_seen, chr_name);
            continue;
        }
        kputc(',',&tmp);
        kputs(key.s,&tmp);
        kputc('=',&tmp);
        if ( quoted ) kputc('"',&tmp);
        kputs(val.s,&tmp);
        if ( quoted ) kputc('"',&tmp);
    }
    if ( !chr_name ) return end;
    ksprintf(dst,"##contig=<ID=%s,length=%d%s>",chr_name,chr_len,tmp.l ? tmp.s : "");
    free(key.s); free(val.s); free(tmp.s);
    return q;
}
char *init_tmp_prefix(const char *tmp_prefix)
{
    kstring_t prefix = {0,0,0};
    if ( tmp_prefix )
    {
        ksprintf(&prefix,"%sXXXXXX",tmp_prefix);
        return prefix.s;
    }

    char *tmpdir = getenv("TMPDIR");
    if ( tmpdir )
        kputs(tmpdir, &prefix);
    else
    {
        #ifdef _WIN32
            char tmp_path[MAX_PATH];
            int ret = GetTempPath(MAX_PATH, tmp_path);
            if (!ret || ret > MAX_PATH)
                error("Could not get the path to the temporary folder\n");
            kputs(tmp_path, &prefix);
        #else
            kputs("/tmp", &prefix);
        #endif
    }
    kputs("/bcftools.XXXXXX", &prefix);
    return prefix.s;
}
static void update_from_fai(args_t *args)
{
    if ( !strcmp("-",args->fname) )
        error("Cannot use the --fai option when reading from standard input.\n");

    faidx_t *fai = fai_load3(args->fai_fname,args->fai_fname,NULL,FAI_FASTA);
    if ( !fai ) error("Could not parse %s\n", args->fai_fname);
    args->rm_tmpfile = init_tmp_prefix(args->tmp_prefix);
    int fd = mkstemp(args->rm_tmpfile);
    if ( fd<0 ) error("Could not open a temporary file for writing: %s\n", args->rm_tmpfile);

    // get a template header: either from the original VCF or from --header
    char *ori_hdr_fname = args->header_fname ? args->header_fname : args->fname;
    htsFile *fp = hts_open(ori_hdr_fname,"r");
    if ( !fp ) error("Failed to open: %s\n", ori_hdr_fname);
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if ( !hdr ) error("Failed to read the header: %s\n", ori_hdr_fname);
    hts_close(fp);  // no need to check the return status here

    // put the header in a text buffer
    kstring_t hdr_txt_ori = {0,0,0}, hdr_txt_new = {0,0,0};
    bcf_hdr_format(hdr, 0, &hdr_txt_ori);
    bcf_hdr_destroy(hdr);

    // update the existing contig lines and remove lines not present in the fai file
    void *chr_seen = khash_str2int_init();
    char *tmp, *beg = hdr_txt_ori.s;
    while ( beg && *beg )
    {
        tmp = strstr(beg, "\n##contig=<");
        if ( !tmp ) break;
        kputsn(beg, tmp-beg+1, &hdr_txt_new);
        size_t l_prev = hdr_txt_new.l;
        beg = copy_and_update_contig_line(fai,tmp+1,chr_seen, &hdr_txt_new);
        if ( l_prev==hdr_txt_new.l ) hdr_txt_new.l--;   // nothing was added, remove the newline
    }
    if ( !beg || !(tmp=strstr(beg,"\n#CHROM")) ) error("Failed to parse the header, #CHROM not found\n");
    kputsn(beg, tmp-beg+1, &hdr_txt_new);

    // add any new contig lines
    int i, n = faidx_nseq(fai);
    for (i=0; i<n; i++)
    {
        if ( khash_str2int_has_key(chr_seen,faidx_iseq(fai,i)) ) continue;
        ksprintf(&hdr_txt_new,"##contig=<ID=%s,length=%d>\n",faidx_iseq(fai,i),faidx_seq_len(fai,faidx_iseq(fai,i)));
    }
    kputs(tmp+1,&hdr_txt_new);

    if ( write(fd, hdr_txt_new.s, hdr_txt_new.l)!=hdr_txt_new.l ) error("Failed to write %zu bytes to %s\n", hdr_txt_new.l,args->rm_tmpfile);
    if ( close(fd)!=0 ) error("Failed to close %s\n", args->rm_tmpfile);
    args->header_fname = args->rm_tmpfile;

    free(hdr_txt_ori.s);
    free(hdr_txt_new.s);
    fai_destroy(fai);
    khash_str2int_destroy_free(chr_seen);
}

static void read_header_file(char *fname, kstring_t *hdr)
{
    kstring_t tmp = {0,0,0};
    hdr->l = 0;

    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);
    while ( hts_getline(fp, KS_SEP_LINE, &tmp) > 0 )
    {
        kputsn(tmp.s,tmp.l,hdr);
        kputc('\n',hdr);
    }
    if ( hts_close(fp) ) error("Close failed: %s\n", fname);
    free(tmp.s);

    while ( hdr->l>0 && isspace(hdr->s[hdr->l-1]) ) hdr->l--;  // remove trailing newlines
    kputc('\n',hdr);
}

static int set_sample_pairs(char **samples, int nsamples, kstring_t *hdr, int idx)
{
    int i, j, n;
    kstring_t key = {0,0,0};
    kstring_t val = {0,0,0};

    // Are these samples "old-name new-name" pairs?
    void *hash = khash_str2str_init();
    for (i=0; i<nsamples; i++)
    {
        char *ptr = samples[i];
        key.l = val.l = 0;
        int escaped = 0;
        while ( *ptr )
        {
            if ( *ptr=='\\' && !escaped ) { escaped = 1; ptr++; continue; }
            if ( isspace(*ptr) && !escaped ) break;
            kputc(*ptr, &key);
            escaped = 0;
            ptr++;
        }
        if ( !*ptr ) break;
        while ( *ptr && isspace(*ptr) ) ptr++;
        while ( *ptr )
        {
            if ( *ptr=='\\' && !escaped ) { escaped = 1; ptr++; continue; }
            if ( isspace(*ptr) && !escaped ) break;
            kputc(*ptr, &val);
            escaped = 0;
            ptr++;
        }
        khash_str2str_set(hash,strdup(key.s),strdup(val.s));
    }
    free(key.s);
    free(val.s);
    if ( i!=nsamples )  // not "old-name new-name" pairs
    {
        khash_str2str_destroy_free_all(hash);
        return 0;
    }

    while ( hdr->l>0 && isspace(hdr->s[hdr->l-1]) ) hdr->l--;  // remove trailing newlines
    hdr->s[hdr->l] = 0;

    kstring_t tmp = {0,0,0};
    i = j = n = 0;  // i:traverse the #CHROM line 1 by 1; j:points to the last column
    while ( hdr->s[idx+i] )
    {
        if ( hdr->s[idx+i]=='\t' )
        {
            hdr->s[idx+i] = 0;

            if ( ++n>9 )
            {
                char *new_name = khash_str2str_get(hash,hdr->s+idx+j);
                kputs(new_name ? new_name : hdr->s+idx+j, &tmp);
            }
            else
                kputs(hdr->s+idx+j, &tmp);

            kputc('\t',&tmp);

            j = ++i;
            continue;
        }
        i++;
    }
    char *new_name = khash_str2str_get(hash,hdr->s+idx+j);
    kputs(new_name ? new_name : hdr->s+idx+j, &tmp);

    khash_str2str_destroy_free_all(hash);

    hdr->l = idx;
    kputs(tmp.s, hdr);
    kputc('\n', hdr);
    free(tmp.s);

    return 1;
}

static void set_samples(char **samples, int nsamples, kstring_t *hdr)
{
    // Find the beginning of the #CHROM line
    int i = hdr->l - 2, ncols = 0;
    while ( i>=0 && hdr->s[i]!='\n' )
    {
        if ( hdr->s[i]=='\t' ) ncols++;
        i--;
    }
    if ( i<0 || strncmp(hdr->s+i+1,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",45) )
    {
        if ( i>0 && !strncmp(hdr->s+i+1,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",38) )
            error("Error: missing FORMAT fields, cowardly refusing to add samples\n");

        error("Could not parse the header: %s\n", hdr->s);
    }

    // Are the samples "old-sample new-sample" pairs?
    if ( set_sample_pairs(samples,nsamples,hdr, i+1) ) return;

    // Replace all samples
    if ( ncols!=nsamples+8 )
        fprintf(stderr, "Warning: different number of samples: %d vs %d\n", nsamples,ncols-8);

    ncols = 0;
    while ( ncols!=9 )
    {
        i++;
        if ( hdr->s[i]=='\t' ) ncols++;
    }
    hdr->l = i;

    for (i=0; i<nsamples; i++)
    {
        kputc('\t', hdr);
        kputs(samples[i], hdr);
    }
    kputc('\n', hdr);
}

BGZF *hts_get_bgzfp(htsFile *fp);
static void reheader_vcf_gz(args_t *args)
{
    BGZF *fp = hts_get_bgzfp(args->fp);
    if ( !fp || bgzf_read_block(fp) != 0 || !fp->block_length )
        error("Failed to read %s: %s\n", args->fname, strerror(errno));

    kstring_t hdr = {0,0,0};
    char *buffer = (char*) fp->uncompressed_block;

    // Read the header and find the position of the data block
    if ( buffer[0]!='#' ) error("Could not parse the header, expected '#', found '%c'\n", buffer[0]);

    int skip_until = 1;     // end of the header in the current uncompressed block
    while (1)
    {
        if ( buffer[skip_until]=='\n' )
        {
            skip_until++;
            if ( skip_until>=fp->block_length )
            {
                kputsn(buffer,skip_until,&hdr);
                if ( bgzf_read_block(fp) != 0 ) error("Error reading %s\n", args->fname);
                if ( !fp->block_length ) break;
                skip_until = 0;
            }
            // The header has finished
            if ( buffer[skip_until]!='#' )
            {
                kputsn(buffer,skip_until,&hdr);
                break;
            }
        }
        skip_until++;
        if ( skip_until>=fp->block_length )
        {
            kputsn(buffer,fp->block_length,&hdr);
            if ( bgzf_read_block(fp) != 0 ) error("Error reading %s\n", args->fname);
            if ( !fp->block_length ) break;
            skip_until = 0;
        }
    }

    int nsamples = 0;
    char **samples = NULL;
    if ( args->samples_fname )
    {
        samples = hts_readlines(args->samples_fname, &nsamples);
        if ( !samples || !nsamples ) error("Error reading the --samples file \"%s\"\n", args->samples_fname);
    }
    if ( args->header_fname )
    {
        free(hdr.s); hdr.s = NULL; hdr.l = hdr.m = 0;
        read_header_file(args->header_fname, &hdr);
    }
    if ( samples )
    {
        set_samples(samples, nsamples, &hdr);
        int i;
        for (i=0; i<nsamples; i++) free(samples[i]);
        free(samples);
    }

    // Output the modified header
    BGZF *bgzf_out = bgzf_open(args->output_fname ? args->output_fname : "-","w");;
    if ( bgzf_write(bgzf_out, hdr.s, hdr.l) < 0 ) error("Can't write BGZF header (code %d)\n", bgzf_out->errcode);
    free(hdr.s);

    // Output all remainig data read with the header block
    if ( fp->block_length - skip_until > 0 )
    {
        if ( bgzf_write(bgzf_out, buffer+skip_until, fp->block_length-skip_until)<0 ) error("Error: %d\n",fp->errcode);
    }
    if ( bgzf_flush(bgzf_out)<0 ) error("Error: %d\n",bgzf_out->errcode);

    // Stream the rest of the file without as it is, without decompressing
    ssize_t nread;
    const size_t page_size = 32768;
    char *buf = (char*) malloc(page_size);
    while (1)
    {
        nread = bgzf_raw_read(fp, buf, page_size);
        if ( nread<=0 ) break;

        int count = bgzf_raw_write(bgzf_out, buf, nread);
        if (count != nread) error("Write failed, wrote %d instead of %d bytes.\n", count,(int)nread);
    }
    if (bgzf_close(bgzf_out) < 0) error("Error closing %s: %d\n",args->output_fname ? args->output_fname : "-",bgzf_out->errcode);
    if (hts_close(args->fp)) error("Error closing %s: %d\n",args->fname,fp->errcode);
    free(buf);
}
static void reheader_vcf(args_t *args)
{
    kstring_t hdr = {0,0,0};
    htsFile *fp = args->fp;
    while ( hts_getline(fp, KS_SEP_LINE, &fp->line) >=0 )
    {
        kputc('\n',&fp->line);  // hts_getline eats the newline character
        if ( fp->line.s[0]!='#' ) break;
        kputsn(fp->line.s,fp->line.l,&hdr);
    }

    int nsamples = 0;
    char **samples = NULL;
    if ( args->samples_fname )
    {
        samples = hts_readlines(args->samples_fname, &nsamples);
        if ( !samples || !nsamples ) error("Error reading the --samples file \"%s\"\n", args->samples_fname);
    }
    if ( args->header_fname )
    {
        free(hdr.s); hdr.s = NULL; hdr.l = hdr.m = 0;
        read_header_file(args->header_fname, &hdr);
    }
    if ( samples )
    {
        set_samples(samples, nsamples, &hdr);
        int i;
        for (i=0; i<nsamples; i++) free(samples[i]);
        free(samples);
    }

    int out = args->output_fname ? open(args->output_fname, O_WRONLY|O_CREAT|O_TRUNC, 0666) : STDOUT_FILENO;
    if ( out==-1 ) error("%s: %s\n", args->output_fname,strerror(errno));
    if ( write(out, hdr.s, hdr.l)!=hdr.l ) error("Failed to write %"PRIu64" bytes\n", (uint64_t)hdr.l);
    free(hdr.s);
    if ( fp->line.l )
    {
        if ( write(out, fp->line.s, fp->line.l)!=fp->line.l ) error("Failed to write %"PRIu64" bytes\n", (uint64_t)fp->line.l);
    }
    while ( hts_getline(fp, KS_SEP_LINE, &fp->line) >=0 )   // uncompressed file implies small size, we don't worry about speed
    {
        kputc('\n',&fp->line);
        if ( write(out, fp->line.s, fp->line.l)!=fp->line.l ) error("Failed to write %"PRIu64" bytes\n", (uint64_t)fp->line.l);
    }
    if ( hts_close(fp)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->fname);
    if ( close(out)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->output_fname);
}

static bcf_hdr_t *strip_header(bcf_hdr_t *src, bcf_hdr_t *dst)
{
    bcf_hrec_t *src_hrec, *dst_hrec, *tmp;
    bcf_hdr_t *out = bcf_hdr_init("r");
    int i;
    for (i=0; i<dst->nhrec; i++)
    {
        // first insert lines which do not code BCF ids, their order does not matter
        dst_hrec = dst->hrec[i];
        if ( dst_hrec->type==BCF_HL_FLT || dst_hrec->type==BCF_HL_INFO || dst_hrec->type==BCF_HL_FMT || dst_hrec->type== BCF_HL_CTG ) continue;
        bcf_hdr_add_hrec(out, bcf_hrec_dup(dst_hrec));
    }
    for (i=0; i<src->nhrec; i++)
    {
        // now transfer header lines which define BCF ids
        src_hrec = src->hrec[i];

        if ( src_hrec->type==BCF_HL_FLT || src_hrec->type==BCF_HL_INFO || src_hrec->type==BCF_HL_FMT || src_hrec->type== BCF_HL_CTG )
        {
            int j = bcf_hrec_find_key(src_hrec, "ID");
            dst_hrec = bcf_hdr_get_hrec(dst, src_hrec->type, "ID", src_hrec->vals[j], NULL);
            if ( !dst_hrec ) continue;

            tmp = bcf_hrec_dup(dst_hrec);

            j = bcf_hrec_find_key(src_hrec, "IDX");
            if ( j>=0 )
            {
                j = atoi(src_hrec->vals[j]);
                if (hrec_add_idx(tmp, j) < 0)
                    error_errno("[%s] Failed to add IDX header", __func__);
            }
            bcf_hdr_add_hrec(out, tmp);
        }
    }
    if (bcf_hdr_sync(out) < 0)
        error_errno("[%s] Failed to update header", __func__);
    for (i=0; i<dst->nhrec; i++)
    {
        // finally add new structured fields
        dst_hrec = dst->hrec[i];
        if ( dst_hrec->type==BCF_HL_FLT || dst_hrec->type==BCF_HL_INFO || dst_hrec->type==BCF_HL_FMT || dst_hrec->type== BCF_HL_CTG )
        {
            int j = bcf_hrec_find_key(dst_hrec, "ID");
            tmp = bcf_hdr_get_hrec(out, dst_hrec->type, "ID", dst_hrec->vals[j], NULL);
            if ( !tmp )
                bcf_hdr_add_hrec(out, bcf_hrec_dup(dst_hrec));
        }
    }
    for (i=0; i<dst->n[BCF_DT_SAMPLE]; i++) bcf_hdr_add_sample(out, dst->samples[i]);
    bcf_hdr_destroy(dst);
    return out;
}

static void reheader_bcf(args_t *args, int is_compressed)
{
    htsFile *fp = args->fp;

    if ( args->n_threads > 0 )
    {
        args->threads = (htsThreadPool *) calloc(1, sizeof(htsThreadPool));
        if ( !args->threads ) error("Could not allocate memory\n");
        if ( !(args->threads->pool = hts_tpool_init(args->n_threads)) ) error("Could not initialize threading\n");
        hts_set_thread_pool(fp, args->threads);
    }

    bcf_hdr_t *hdr = bcf_hdr_read(fp); if ( !hdr ) error("Failed to read the header: %s\n", args->fname);
    kstring_t htxt = {0,0,0};
    bcf_hdr_format(hdr, 1, &htxt);

    int i, nsamples = 0;
    char **samples = NULL;
    if ( args->samples_fname )
    {
        samples = hts_readlines(args->samples_fname, &nsamples);
        if ( !samples || !nsamples ) error("Error reading the --samples file \"%s\"\n", args->samples_fname);
    }
    if ( args->header_fname )
    {
        free(htxt.s); htxt.s = NULL; htxt.l = htxt.m = 0;
        read_header_file(args->header_fname, &htxt);
    }
    if ( samples )
    {
        set_samples(samples, nsamples, &htxt);
        for (i=0; i<nsamples; i++) free(samples[i]);
        free(samples);
    }

    bcf_hdr_t *hdr_out = bcf_hdr_init("r");
    if ( bcf_hdr_parse(hdr_out, htxt.s) < 0 ) error("An error occurred while parsing the header\n");
    if ( args->header_fname ) hdr_out = strip_header(hdr, hdr_out);

    // write the header and the body
    htsFile *fp_out = hts_open(args->output_fname ? args->output_fname : "-",is_compressed ? "wb" : "wbu");
    if ( !fp_out ) error("%s: %s\n", args->output_fname ? args->output_fname : "-", strerror(errno));
    if ( args->threads )
        hts_set_thread_pool(fp_out, args->threads);
    if ( bcf_hdr_write(fp_out, hdr_out)!=0 ) error("[%s] Error: cannot write the header to %s\n", __func__,args->output_fname ? args->output_fname : "standard output");

    bcf1_t *rec = bcf_init();
    while ( bcf_read(fp, hdr, rec)==0 )
    {
        // sanity checking, this slows things down. Make it optional?
        bcf_unpack(rec, BCF_UN_ALL);
        if ( rec->rid >= hdr_out->n[BCF_DT_CTG] || strcmp(bcf_hdr_int2id(hdr,BCF_DT_CTG,rec->rid),bcf_hdr_int2id(hdr_out,BCF_DT_CTG,rec->rid)) )
            error("The CHROM is not defined: \"%s\"\n", bcf_hdr_int2id(hdr,BCF_DT_CTG,rec->rid));

        for (i=0; i<rec->d.n_flt; i++)
        {
            int id = rec->d.flt[i];
            if ( id >= hdr_out->n[BCF_DT_ID] ) break;
            if ( !bcf_hdr_idinfo_exists(hdr_out,BCF_HL_FLT,id) ) break;
            if ( strcmp(hdr->id[BCF_DT_ID][id].key,hdr_out->id[BCF_DT_ID][id].key) )
                error("FIXME: Broken FILTER ids: %s vs %s\n", hdr->id[BCF_DT_ID][id].key,hdr_out->id[BCF_DT_ID][id].key);
        }
        if ( i!=rec->d.n_flt )
            error("The FILTER is not defined: \"%s\"\n", bcf_hdr_int2id(hdr,BCF_DT_ID,rec->d.flt[i]));

        for (i=0; i<rec->n_info; i++)
        {
            int id = rec->d.info[i].key;
            if ( id >= hdr_out->n[BCF_DT_ID] ) break;
            if ( !hdr_out->id[BCF_DT_ID][id].key ) break;
            if ( !bcf_hdr_idinfo_exists(hdr_out,BCF_HL_INFO,id) ) break;
            if ( strcmp(hdr->id[BCF_DT_ID][id].key,hdr_out->id[BCF_DT_ID][id].key) )
                error("FIXME: Broken INFO ids: %s vs %s\n", hdr->id[BCF_DT_ID][id].key,hdr_out->id[BCF_DT_ID][id].key);
        }
        if ( i!=rec->n_info )
            error("The INFO tag is not defined: \"%s\"\n", bcf_hdr_int2id(hdr,BCF_DT_ID,rec->d.info[i].key));

        for (i=0; i<rec->n_fmt; i++)
        {
            int id = rec->d.fmt[i].id;
            if ( id >= hdr_out->n[BCF_DT_ID] ) break;
            if ( !hdr_out->id[BCF_DT_ID][id].key ) break;
            if ( !bcf_hdr_idinfo_exists(hdr_out,BCF_HL_FMT,id) ) break;
            if ( strcmp(hdr->id[BCF_DT_ID][id].key,hdr_out->id[BCF_DT_ID][id].key) )
                error("FIXME: Broken FORMAT ids: %s vs %s\n", hdr->id[BCF_DT_ID][id].key,hdr_out->id[BCF_DT_ID][id].key);
        }
        if ( i!=rec->n_fmt )
            error("The FORMAT tag is not defined: \"%s\"\n", bcf_hdr_int2id(hdr,BCF_DT_ID,rec->d.fmt[i].id));

        if ( bcf_write(fp_out,hdr_out,rec)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname ? args->output_fname : "standard output");
    }
    bcf_destroy(rec);

    free(htxt.s);
    if ( hts_close(fp_out)!=0 ) error("[%s] Error: failed to close the file %s\n",__func__,args->output_fname ? args->output_fname : "standard output");
    if ( hts_close(fp)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->fname);
    bcf_hdr_destroy(hdr_out);
    bcf_hdr_destroy(hdr);
    if ( args->threads )
    {
        hts_tpool_destroy(args->threads->pool);
        free(args->threads);
    }
}


static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Modify header of VCF/BCF files, change sample names.\n");
    fprintf(stderr, "Usage:   bcftools reheader [OPTIONS] <in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -f, --fai FILE             update sequences and their lengths from the .fai file\n");
    fprintf(stderr, "    -h, --header FILE          new header\n");
    fprintf(stderr, "    -o, --output FILE          write output to a file [standard output]\n");
    fprintf(stderr, "    -s, --samples FILE         new sample names\n");
#ifdef _WIN32
    fprintf(stderr, "    -T, --temp-prefix PATH     template for temporary file name [/bcftools.XXXXXX]\n");
#else
    fprintf(stderr, "    -T, --temp-prefix PATH     template for temporary file name [/tmp/bcftools.XXXXXX]\n");
#endif
    fprintf(stderr, "        --threads INT          use multithreading with <int> worker threads (BCF only) [0]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Example:\n");
    fprintf(stderr, "   # Write out the header to be modified\n");
    fprintf(stderr, "   bcftools view -h old.bcf > header.txt\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   # Edit the header using your favorite text editor\n");
    fprintf(stderr, "   vi header.txt\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   # Reheader the file\n");
    fprintf(stderr, "   bcftools reheader -h header.txt -o new.bcf old.bcf\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_reheader(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    
    static struct option loptions[] =
    {
        {"temp-prefix",1,0,'T'},
        {"fai",1,0,'f'},
        {"output",1,0,'o'},
        {"header",1,0,'h'},
        {"samples",1,0,'s'},
        {"threads",1,NULL,1},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "s:h:o:f:T:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case  1 : args->n_threads = strtol(optarg, 0, 0); break;
            case 'T': args->tmp_prefix = optarg; break;
            case 'f': args->fai_fname = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 's': args->samples_fname = optarg; break;
            case 'h': args->header_fname = optarg; break;
            case '?': usage(args); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else usage(args);
    }
    else args->fname = argv[optind];

    if ( args->fai_fname ) update_from_fai(args);
    if ( !args->samples_fname && !args->header_fname ) usage(args);
    if ( !args->fname ) usage(args);

    args->fp = hts_open(args->fname,"r");
    if ( !args->fp ) error("Failed to read from %s\n", !strcmp("-",args->fname)?"standard input":args->fname);
    args->type = *hts_get_format(args->fp);

    if ( args->type.format==vcf )
    {
        if ( args->type.compression==bgzf )
            reheader_vcf_gz(args);
        else if ( args->type.compression==no_compression )
            reheader_vcf(args);
        else if ( args->type.compression==gzip )
            error("Error: cannot reheader gzip-compressed files, first convert with `bcftools view --output-type` to a supported format\n");
        else
            error("Error: the compression type of \"%s\" is not recognised/supported\n", args->fname);
    }
    else
        reheader_bcf(args, args->type.compression==bgzf || args->type.compression==gzip);

    if ( args->rm_tmpfile )
    {
        unlink(args->rm_tmpfile);
        free(args->rm_tmpfile);
    }
    free(args);
    return 0;
}
