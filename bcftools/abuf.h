/* The MIT License

   Copyright (c) 2021 Genome Research Ltd.

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
    Atomize/deatomize complex variants
*/

#ifndef __ABUF_H__
#define __ABUF_H__

#include <htslib/vcf.h>

typedef struct _abuf_t abuf_t;

// Modes of operation
typedef enum
{
    NONE,

    // mode of operation, to be passed to abuf_init
    SPLIT,
    JOIN,

    BCF_HDR,        // should the records be annotated, a writable bcf header is required
    INFO_TAG,       // set BCF_HDR first
    STAR_ALLELE     // 1: use STAR allele (the default), 0: set overlaps to missing
}
abuf_opt_t;

#define abuf_set_opt(buf,type,key,value) { type tmp = value; abuf_set(buf, key, (void*)&tmp); }
void abuf_set(abuf_t *buf, abuf_opt_t key, void *value);

/*
 *  abuf_init() - init buffer
 *  @win:   number of sites (>0) or bp (<0)
 */
abuf_t *abuf_init(const bcf_hdr_t *hdr, abuf_opt_t mode);
void abuf_destroy(abuf_t *buf);

/*
 *  abuf_push() - Push a new site for analysis
 */
void abuf_push(abuf_t *buf, bcf1_t *rec);

/*
 *  abuf_flush() - Return next buffered record
 *  @flush_all: Set to 1 if no more overlapping records are coming (e.g. end of chromosome or end of file),
 *              the buffer can be emptied.
 *  return:     The next atomized/deatomized VCF record or NULL if no record is ready. The returned
 *              structure will be cleaned by abuf.
 */
bcf1_t *abuf_flush(abuf_t *buf, int flush_all);

#endif

