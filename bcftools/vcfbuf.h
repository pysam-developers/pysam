/* The MIT License

   Copyright (c) 2017 Genome Research Ltd.

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
    Buffer VCF records and perform operations on the buffer
*/

#ifndef __VCFBUF_H__
#define __VCFBUF_H__

#include <htslib/vcf.h>

typedef struct _vcfbuf_t vcfbuf_t;

// Modes of operation
typedef enum
{
    VCFBUF_LD_MAX,          // vcfbuf_max_ld() stops at the first record that exceeds the threshold
    VCFBUF_RAND_MISSING,    // randomize rather than ignore missing genotypes
    VCFBUF_SKIP_FILTER,     // skip sites with FILTER diferent from "PASS" or "."
    VCFBUF_NSITES,          // leave at max this many sites in the window
    VCFBUF_AF_TAG,          // use this INFO tag with LD_NSITES
    VCFBUF_OVERLAP_WIN,     // keep only overlapping variants in the window
}
vcfbuf_opt_t;

#define vcfbuf_set_opt(buf,type,key,value) { type tmp = value; vcfbuf_set(buf, key, (void*)&tmp); }
void vcfbuf_set(vcfbuf_t *buf, vcfbuf_opt_t key, void *value);


/*
 *  vcfbuf_init() - init buffer
 *  @win:   number of sites (>0) or bp (<0)
 */
vcfbuf_t *vcfbuf_init(bcf_hdr_t *hdr, int win);
void vcfbuf_destroy(vcfbuf_t *buf);

/*
 *  vcfbuf_push() - push a new site for analysis
 *  @swap:  if set, do not create a copy, but return a substitute
 */
bcf1_t *vcfbuf_push(vcfbuf_t *buf, bcf1_t *rec, int swap);

bcf1_t *vcfbuf_flush(vcfbuf_t *buf, int flush_all);

/*
 *  vcfbuf_nsites() - return the number of sites in the buffer
 */
int vcfbuf_nsites(vcfbuf_t *buf);

/*
 *  vcfbuf_max_ld() - return a record that has maximum D or first record exceeding the threshold
 *  @ld:        will be filled with the maximum D found
 */
bcf1_t *vcfbuf_max_ld(vcfbuf_t *buf, bcf1_t *rec, double *ld);

#endif

