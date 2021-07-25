/*  ext-sort.h -- sort on disk

   Copyright (C) 2020-2021 Genome Research Ltd.

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
#include <unistd.h>     // for unlink()
#include <sys/stat.h>   // for chmod()
#include <assert.h>
#include <fcntl.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include "bcftools.h"
#include "extsort.h"
#include "kheap.h"

typedef struct
{
    extsort_t *es;  // this is to get access to extsort_cmp_f from kheap
    int fd;
    char *fname;
    void *dat;
}
blk_t;

static inline int blk_is_smaller(blk_t **aptr, blk_t **bptr);
KHEAP_INIT(blk, blk_t*, blk_is_smaller)     /* defines khp_blk_t */

struct _extsort_t
{
    size_t dat_size, mem, max_mem;
    char *tmp_prefix;
    extsort_cmp_f cmp;

    size_t nbuf, mbuf, nblk;
    blk_t **blk;
    void **buf, *tmp_dat;
    khp_blk_t *bhp;
};

static inline int blk_is_smaller(blk_t **aptr, blk_t **bptr)
{
    blk_t *a = *aptr;
    blk_t *b = *bptr;
    int ret = a->es->cmp(&a->dat,&b->dat);
    if ( ret < 0 ) return 1;
    return 0;
}

size_t parse_mem_string(const char *str);

void extsort_set(extsort_t *es, extsort_opt_t key, void *value)
{
    if ( key==DAT_SIZE ) { es->dat_size = *((size_t*)value); return; }
    if ( key==MAX_MEM )
    {
        es->max_mem = parse_mem_string(*((const char**)value));
        if ( es->max_mem <=0 ) error("Could not parse the memory string, expected positive number: %s\n",*((const char**)value));
        return;
    }
    if ( key==TMP_PREFIX ) { es->tmp_prefix = init_tmp_prefix(*((const char**)value)); return; }
    if ( key==FUNC_CMP ) { es->cmp = *((extsort_cmp_f*)value); return; }
}

extsort_t *extsort_alloc(void)
{
    extsort_t *es = (extsort_t*) calloc(1,sizeof(*es));
    es->max_mem = 100e6;
    return es;
}
void extsort_init(extsort_t *es)
{
    assert( es->cmp );
    assert( es->dat_size );
    if ( !es->tmp_prefix ) es->tmp_prefix = init_tmp_prefix(NULL);
    es->tmp_dat = malloc(es->dat_size);
}

void extsort_destroy(extsort_t *es)
{
    int i;
    for (i=0; i<es->nblk; i++)
    {
        blk_t *blk = es->blk[i];
        if ( blk->fd!=-1 )
#ifdef _WIN32
            _close(blk->fd);
#else
            close(blk->fd);
#endif
        free(blk->fname);
        free(blk->dat);
        free(blk);
    }
    free(es->tmp_dat);
    free(es->tmp_prefix);
    free(es->blk);
    khp_destroy(blk, es->bhp);
    free(es);
}

static void _buf_flush(extsort_t *es)
{
    int i;
    if ( !es->nbuf ) return;

    qsort(es->buf, es->nbuf, sizeof(void*), es->cmp);

    es->nblk++;
    es->blk = (blk_t**) realloc(es->blk, sizeof(blk_t*)*es->nblk);
    es->blk[es->nblk-1] = (blk_t*) calloc(1,sizeof(blk_t));
    blk_t *blk = es->blk[es->nblk-1];
    blk->es    = es;
    blk->dat   = malloc(es->dat_size);
    blk->fname = strdup(es->tmp_prefix);
    #ifdef _WIN32
        for (i=0; i<100000; i++)
        {
            memcpy(blk->fname,es->tmp_prefix,strlen(es->tmp_prefix));
            mktemp(blk->fname);
            blk->fd = _open(blk->fname, O_RDWR|O_CREAT|O_EXCL|O_BINARY|O_TEMPORARY, 0600);
            if ( blk->fd==-1 )
            {
                if ( errno==EEXIST ) continue; 
                error("Error: failed to open a temporary file %s\n",blk->fname);
            }
            break;
        }
        if ( !blk->fd ) error("Error: failed to create a unique temporary file name from %s\n",es->tmp_prefix);
        if ( _chmod(blk->fname, S_IRUSR|S_IWUSR)!=0 ) error("Error: failed to set permissions of the temporary file %s\n",blk->fname);
    #else
        if ( (blk->fd = mkstemp(blk->fname))==-1 )
            error("Error: failed to open a temporary file %s\n",blk->fname);
        if ( fchmod(blk->fd,S_IRUSR|S_IWUSR)!=0 ) error("Error: failed to set permissions of the temporary file %s\n",blk->fname);
        unlink(blk->fname); // should auto delete when closed on linux, the descriptor remains open
    #endif

    for (i=0; i<es->nbuf; i++)
    {
        #ifdef _WIN32
            if ( _write(blk->fd, es->buf[i], es->dat_size)!=es->dat_size ) error("Error: failed to write %zu bytes to the temporary file %s\n",es->dat_size,blk->fname);
        #else
            if ( write(blk->fd, es->buf[i], es->dat_size)!=es->dat_size ) error("Error: failed to write %zu bytes to the temporary file %s\n",es->dat_size,blk->fname);
        #endif
        free(es->buf[i]);
    }
#ifdef _WIN32
    if ( _lseek(blk->fd,0,SEEK_SET)!=0 ) error("Error: failed to lseek() to the start of the temporary file %s\n", blk->fname);
#else
    if ( lseek(blk->fd,0,SEEK_SET)!=0 ) error("Error: failed to lseek() to the start of the temporary file %s\n", blk->fname);
#endif

    es->nbuf = 0;
    es->mem  = 0;
}

void extsort_push(extsort_t *es, void *dat)
{
    int delta = sizeof(void*) + es->dat_size;
    if ( es->nbuf && es->mem + delta > es->max_mem ) _buf_flush(es);
    es->nbuf++;
    es->mem += delta;
    hts_expand(void*, es->nbuf, es->mbuf, es->buf);
    es->buf[es->nbuf-1] = dat;
}

// return number of elements read
static ssize_t _blk_read(extsort_t *es, blk_t *blk)
{
    ssize_t ret = 0;
    if ( blk->fd==-1 ) return ret;
#ifdef _WIN32
    ret = _read(blk->fd, blk->dat, es->dat_size);
#else
    ret = read(blk->fd, blk->dat, es->dat_size);
#endif
    if ( ret < 0 ) error("Error: failed to read from the temporary file %s\n", blk->fname);
    if ( ret == 0 )
    {
#ifdef _WIN32
        if ( _close(blk->fd)!=0 ) error("Error: failed to close the temporary file %s\n", blk->fname);
#else
        if ( close(blk->fd)!=0 ) error("Error: failed to close the temporary file %s\n", blk->fname);
#endif
        blk->fd = -1;
        return ret;
    }
    if ( ret < es->dat_size ) error("Error: failed to read %zu bytes from the temporary file %s\n",es->dat_size,blk->fname);
    return ret;
}

void extsort_sort(extsort_t *es)
{
    _buf_flush(es);
    free(es->buf);
    es->buf = NULL;
    es->bhp = khp_init(blk);

    // open all blocks, read one record from each, create a heap
    int i;
    for (i=0; i<es->nblk; i++)
    {
        blk_t *blk = es->blk[i];
#ifdef _WIN32
        if ( _lseek(blk->fd,0,SEEK_SET)!=0 ) error("Error: failed to lseek() to the start of the temporary file %s\n", blk->fname);
#else
        if ( lseek(blk->fd,0,SEEK_SET)!=0 ) error("Error: failed to lseek() to the start of the temporary file %s\n", blk->fname);
#endif
        int ret = _blk_read(es, blk);
        if ( ret ) khp_insert(blk, es->bhp, &blk);
    }
}

void *extsort_shift(extsort_t *es)
{
    if ( !es->bhp->ndat ) return NULL;
    blk_t *blk = es->bhp->dat[0];

    // swap the pointer which keeps the location of user data so that it is not overwritten by the next read
    void *tmp = es->tmp_dat; es->tmp_dat = blk->dat; blk->dat = tmp;
    khp_delete(blk, es->bhp);

    int ret = _blk_read(es, blk);
    if ( ret ) khp_insert(blk, es->bhp, &blk);

    return es->tmp_dat;
}

