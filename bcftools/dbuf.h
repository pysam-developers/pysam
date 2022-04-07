/* The MIT License

   Copyright (c) 2022 Genome Research Ltd.

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

/*
    Simple data buffer
*/

#ifndef __DBUF_H__
#define __DBUF_H__

#include <htslib/hts.h>

typedef struct
{
    size_t n,m;
    void **dat;
}
dbuf_t;

static inline dbuf_t *dbuf_push(dbuf_t *buf, void *ptr)
{
    if ( !buf ) buf = calloc(1,sizeof(dbuf_t));
    buf->n++;
    hts_expand(void*,buf->n,buf->m,buf->dat);
    buf->dat[buf->n-1] = ptr;
    return buf;
}

static inline void *dbuf_ith(dbuf_t *buf, int i)
{
    return buf->dat[i];
}

static inline size_t dbuf_n(dbuf_t *buf)
{
    return buf->n;
}

static inline void dbuf_destroy_free(dbuf_t *buf)
{
    int i;
    for (i=0; i<buf->n; i++) free(buf->dat[i]);
    free(buf->dat);
    free(buf);
}

#endif

