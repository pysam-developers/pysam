/*  vcfsort.c -- sort subcommand

   Copyright (C) 2017-2022 Genome Research Ltd.

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
   THE SOFTWARE.

 */

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <math.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/hts_os.h>
#include "kheap.h"
#include "bcftools.h"

typedef struct
{
    char *fname;
    htsFile *fh;
    bcf1_t *rec;
}
blk_t;

typedef struct _args_t
{
    bcf_hdr_t *hdr;
    char **argv, *fname, *output_fname, *tmp_dir;
    int argc, output_type, clevel;
    size_t max_mem, mem;
    bcf1_t **buf;
    uint8_t *mem_block;
    size_t nbuf, mbuf, nblk;
    blk_t *blk;
}
args_t;

void clean_files(args_t *args)
{
    int i;
    fprintf(stderr,"Cleaning\n");
    for (i=0; i<args->nblk; i++)
    {
        blk_t *blk = args->blk + i;
        if ( blk->fname )
        {
            unlink(blk->fname);
            free(blk->fname);
        }
        if ( blk->rec )
            bcf_destroy(blk->rec);
    }
    rmdir(args->tmp_dir);
}
void clean_files_and_throw(args_t *args, const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    clean_files(args);
    exit(-1);
}

int cmp_bcf_pos(const void *aptr, const void *bptr)
{
    bcf1_t *a = *((bcf1_t**)aptr);
    bcf1_t *b = *((bcf1_t**)bptr);
    if ( a->rid < b->rid ) return -1;
    if ( a->rid > b->rid ) return 1;
    if ( a->pos < b->pos ) return -1;
    if ( a->pos > b->pos ) return 1;

    // Sort the same chr:pos records lexicographically by ref,alt.
    // This will be called rarely so should not slow the sorting down
    // noticeably.

    int i;
    for (i=0; i<a->n_allele; i++)
    {
        if ( i >= b->n_allele ) return 1;
        int ret = strcasecmp(a->d.allele[i],b->d.allele[i]);
        if ( ret ) return ret;
    }
    if ( a->n_allele < b->n_allele ) return -1;
    return 0;
}

void buf_flush(args_t *args)
{
    if ( !args->nbuf ) return;

    qsort(args->buf, args->nbuf, sizeof(*args->buf), cmp_bcf_pos);

    args->nblk++;
    args->blk = (blk_t*) realloc(args->blk, sizeof(blk_t)*args->nblk);
    if ( !args->blk ) error("Error: could not allocate %zu bytes of memory, try reducing --max-mem\n",sizeof(blk_t)*args->nblk);
    blk_t *blk = args->blk + args->nblk - 1;

    kstring_t str = {0,0,0};
    ksprintf(&str, "%s/%05d.bcf", args->tmp_dir, (int)args->nblk);
    blk->fname = str.s;
    blk->rec   = NULL;
    blk->fh    = NULL;

    htsFile *fh = hts_open(blk->fname, "wbu");
    if ( fh == NULL ) clean_files_and_throw(args, "Cannot write %s: %s\n", blk->fname, strerror(errno));
    if ( bcf_hdr_write(fh, args->hdr)!=0 ) clean_files_and_throw(args, "[%s] Error: cannot write to %s\n", __func__,blk->fname);

    int i;
    for (i=0; i<args->nbuf; i++)
    {
        if ( bcf_write(fh, args->hdr, args->buf[i])!=0 ) clean_files_and_throw(args, "[%s] Error: cannot write to %s\n", __func__,blk->fname);
    }
    if ( hts_close(fh)!=0 ) clean_files_and_throw(args, "[%s] Error: close failed .. %s\n", __func__,blk->fname);

    args->nbuf = 0;
    args->mem  = 0;
}


static inline uint8_t *_align_up(uint8_t *ptr)
{
    return (uint8_t*)(((size_t)ptr + 8 - 1) & ~((size_t)(8 - 1)));
}

void buf_push(args_t *args, bcf1_t *rec)
{
    size_t delta = sizeof(bcf1_t) + rec->shared.l + rec->indiv.l + rec->unpack_size[0] + rec->unpack_size[1]
        + sizeof(*rec->d.allele)*rec->d.m_allele
        + sizeof(bcf1_t*)       // args->buf
        + 8;                    // the number of _align_up() calls

    if ( delta > args->max_mem - args->mem )
    {
        args->nbuf++;
        hts_expand(bcf1_t*, args->nbuf, args->mbuf, args->buf);
        args->buf[args->nbuf-1] = rec;
        buf_flush(args);
        bcf_destroy(rec);
        return;
    }

    // make sure nothing has changed in htslib
    assert( rec->unpacked==BCF_UN_STR && !rec->d.flt && !rec->d.info && !rec->d.fmt && !rec->d.var );

    uint8_t *ptr_beg = args->mem_block + args->mem;
    uint8_t *ptr = _align_up(ptr_beg);
    bcf1_t *new_rec = (bcf1_t*)ptr;
    memcpy(new_rec,rec,sizeof(*rec));
    ptr += sizeof(*rec);

    // The array of allele pointers does not need alignment as bcf1_t is already padded to the biggest
    // data type in the structure
    char **allele = (char**)ptr;
    ptr += rec->n_allele*sizeof(*allele);

    // This is just to prevent valgrind from complaining about memcpy, unpack_size is a high-water mark
    // and the end may be uninitialized
    delta = rec->d.allele[rec->n_allele-1] - rec->d.allele[0];
    while ( delta < rec->unpack_size[1] ) if ( !rec->d.als[delta++] ) break;
    memcpy(ptr,rec->d.als,delta);
    new_rec->d.als = (char*)ptr;
    ptr = ptr + delta;

    int i;
    for (i=0; i<rec->n_allele; i++) allele[i] = new_rec->d.als + (ptrdiff_t)(rec->d.allele[i] - rec->d.allele[0]);
    new_rec->d.allele = allele;

    memcpy(ptr,rec->shared.s,rec->shared.l);
    new_rec->shared.s = (char*)ptr;
    new_rec->shared.m = rec->shared.l;
    ptr += rec->shared.l;

    memcpy(ptr,rec->indiv.s,rec->indiv.l);
    new_rec->indiv.s = (char*)ptr;
    new_rec->indiv.m = rec->indiv.l;
    ptr += rec->indiv.l;

    // This is just to prevent valgrind from complaining about memcpy, unpack_size is a high-water mark
    // and the end may be uninitialized
    i = 0;
    while ( i < rec->unpack_size[0] ) if ( !rec->d.id[i++] ) break;
    memcpy(ptr,rec->d.id,i);
    new_rec->d.id = (char*)ptr;
    ptr += i;

    args->nbuf++;
    hts_expand(bcf1_t*, args->nbuf, args->mbuf, args->buf);
    args->buf[args->nbuf-1] = new_rec;

    delta = ptr - ptr_beg;
    args->mem += delta;

    assert( args->mem <= args->max_mem );

    bcf_destroy(rec);
}

void sort_blocks(args_t *args)
{
    htsFile *in = hts_open(args->fname, "r");
    if ( !in ) clean_files_and_throw(args, "Could not read %s\n", args->fname);
    args->hdr = bcf_hdr_read(in);
    if ( !args->hdr) clean_files_and_throw(args, "Could not read VCF/BCF headers from %s\n", args->fname);

    while ( 1 )
    {
        bcf1_t *rec = bcf_init();
        int ret = bcf_read1(in, args->hdr, rec);
        if ( ret < -1 ) clean_files_and_throw(args,"Error encountered while parsing the input\n");
        if ( ret == -1 )
        {
            bcf_destroy(rec);
            break;
        }
        if ( rec->errcode ) clean_files_and_throw(args,"Error encountered while parsing the input at %s:%d\n",bcf_seqname(args->hdr,rec),rec->pos+1);
        bcf_unpack(rec, BCF_UN_STR);
        buf_push(args, rec);
    }
    buf_flush(args);
    free(args->buf);

    if ( hts_close(in)!=0 ) clean_files_and_throw(args,"Close failed: %s\n", args->fname);
}

static inline int blk_is_smaller(blk_t **aptr, blk_t **bptr)
{
    blk_t *a = *aptr;
    blk_t *b = *bptr;
    int ret = cmp_bcf_pos(&a->rec, &b->rec);
    if ( ret < 0 ) return 1;
    return 0;
}
KHEAP_INIT(blk, blk_t*, blk_is_smaller)

void blk_read(args_t *args, khp_blk_t *bhp, bcf_hdr_t *hdr, blk_t *blk)
{
    if ( !blk->fh ) return;
    int ret = bcf_read(blk->fh, hdr, blk->rec);
    if ( ret < -1 ) clean_files_and_throw(args, "Error reading %s\n", blk->fname);
    if ( ret == -1 )
    {
        if ( hts_close(blk->fh)!=0 ) clean_files_and_throw(args, "Close failed: %s\n", blk->fname);
        blk->fh = 0;
        return;
    }
    bcf_unpack(blk->rec, BCF_UN_STR);
    khp_insert(blk, bhp, &blk);
}

void merge_blocks(args_t *args)
{
    fprintf(stderr,"Merging %d temporary files\n", (int)args->nblk);
    khp_blk_t *bhp = khp_init(blk);

    int i;
    for (i=0; i<args->nblk; i++)
    {
        blk_t *blk = args->blk + i;
        blk->fh = hts_open(blk->fname, "r");
        if ( !blk->fh ) clean_files_and_throw(args, "Could not read %s: %s\n", blk->fname, strerror(errno));
        bcf_hdr_t *hdr = bcf_hdr_read(blk->fh);
        bcf_hdr_destroy(hdr);
        blk->rec = bcf_init();
        blk_read(args, bhp, args->hdr, blk);
    }

    char wmode[8];
    set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
    htsFile *out = hts_open(args->output_fname ? args->output_fname : "-", wmode);
    if ( bcf_hdr_write(out, args->hdr)!=0 ) clean_files_and_throw(args, "[%s] Error: cannot write to %s\n", __func__,args->output_fname);
    while ( bhp->ndat )
    {
        blk_t *blk = bhp->dat[0];
        if ( bcf_write(out, args->hdr, blk->rec)!=0 ) clean_files_and_throw(args, "[%s] Error: cannot write to %s\n", __func__,args->output_fname);
        khp_delete(blk, bhp);
        blk_read(args, bhp, args->hdr, blk);
    }
    if ( hts_close(out)!=0 ) clean_files_and_throw(args, "Close failed: %s\n", args->output_fname);

    clean_files(args);

    free(args->blk);
    khp_destroy(blk, bhp);
    fprintf(stderr,"Done\n");
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Sort VCF/BCF file.\n");
    fprintf(stderr, "Usage:   bcftools sort [OPTIONS] <FILE.vcf>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -m, --max-mem FLOAT[kMG]       maximum memory to use [768M]\n");    // using metric units, 1M=1e6
    fprintf(stderr, "    -o, --output FILE              output file name [stdout]\n");
    fprintf(stderr, "    -O, --output-type b|u|z|v      b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "    -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n");

#ifdef _WIN32
    fprintf(stderr, "    -T, --temp-dir DIR             temporary files [/bcftools.XXXXXX]\n");
#else
    fprintf(stderr, "    -T, --temp-dir DIR             temporary files [/tmp/bcftools.XXXXXX]\n");
#endif
    fprintf(stderr, "\n");
    exit(1);
}

size_t parse_mem_string(const char *str)
{
    char *tmp;
    double mem = strtod(str, &tmp);
    if ( tmp==str ) error("Could not parse the memory string: \"%s\"\n", str);
    if ( !strcasecmp("k",tmp) ) mem *= 1000;
    else if ( !strcasecmp("m",tmp) ) mem *= 1000*1000;
    else if ( !strcasecmp("g",tmp) ) mem *= 1000*1000*1000;
    return mem;
}

void mkdir_p(const char *fmt, ...);
static void init(args_t *args)
{
    args->max_mem *= 0.9;
    args->mem_block = malloc(args->max_mem);
    if ( !args->mem_block ) error("Error: could not allocate %zu bytes of memory, try reducing --max-mem\n",args->max_mem);
    args->mem = 0;

    args->tmp_dir = init_tmp_prefix(args->tmp_dir);

#ifdef _WIN32
        int ret = mkdir(mktemp(args->tmp_dir), 0700);
        if ( ret ) error("mkdir(%s) failed: %s\n", args->tmp_dir,strerror(errno));
#else
        char *tmp = mkdtemp(args->tmp_dir);
        if ( !tmp ) error("mkdtemp(%s) failed: %s\n",  args->tmp_dir,strerror(errno));
        int ret = chmod(tmp, S_IRUSR|S_IWUSR|S_IXUSR);
        if ( ret ) error("chmod(%s,S_IRUSR|S_IWUSR|S_IXUSR) failed: %s\n", args->tmp_dir,strerror(errno));
#endif

    fprintf(stderr,"Writing to %s\n", args->tmp_dir);
}
static void destroy(args_t *args)
{
    bcf_hdr_destroy(args->hdr);
    free(args->mem_block);
    free(args->tmp_dir);
    free(args);
}

int main_sort(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->max_mem = 768*1000*1000;
    args->output_fname = "-";
    args->clevel = -1;

    static struct option loptions[] =
    {
        {"max-mem",required_argument,NULL,'m'},
        {"temp-dir",required_argument,NULL,'T'},
        {"output-type",required_argument,NULL,'O'},
        {"output-file",required_argument,NULL,'o'},
        {"output",required_argument,NULL,'o'},
        {"help",no_argument,NULL,'h'},
        {0,0,0,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "m:T:O:o:h?",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'm': args->max_mem = parse_mem_string(optarg); break;
            case 'T': args->tmp_dir = optarg; break;
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
            case 'h':
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

    init(args);
    sort_blocks(args);
    merge_blocks(args);
    destroy(args);

    return 0;
}
