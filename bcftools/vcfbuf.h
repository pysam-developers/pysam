/* The MIT License

   Copyright (c) 2017-2021 Genome Research Ltd.

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
    VCFBUF_OVERLAP_WIN,     // keep only overlapping variants in the window
    VCFBUF_RMDUP,           // remove duplicate sites (completely)
    VCFBUF_NSITES,          // leave at max this many sites in the window
    VCFBUF_NSITES_MODE,     // one of: maxAF (keep sites with max AF), 1st (sites that come first), rand (pick randomly)
    VCFBUF_AF_TAG,          // use this INFO tag with VCFBUF_NSITES

    // LD related options
    LD_RAND_MISSING,        // randomize rather than ignore missing genotypes
    LD_FILTER1,             // exclude the next record inserted by vcfbuf_push() from LD analysis
    LD_MAX_R2,              // If set, vcfbuf_ld() will stop at the first record that exceeds the R2,
    LD_MAX_LD,              //      LD, or HD threshold. When multiple are set, the OR logic is applied
    LD_MAX_HD,              //      
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
 */
bcf1_t *vcfbuf_push(vcfbuf_t *buf, bcf1_t *rec);

/*
 *  vcfbuf_peek() - return pointer to i-th record in the buffer but do not remove it from the buffer
 *  @idx:  0-based index to buffered lines
 */
bcf1_t *vcfbuf_peek(vcfbuf_t *buf, int idx);

/*
 *  vcfbuf_remove() - return pointer to i-th record in the buffer and remove it from the buffer
 *  @idx:  0-based index to buffered lines
 */
bcf1_t *vcfbuf_remove(vcfbuf_t *buf, int idx);

bcf1_t *vcfbuf_flush(vcfbuf_t *buf, int flush_all);

/*
 *  vcfbuf_nsites() - return the number of sites in the buffer
 */
int vcfbuf_nsites(vcfbuf_t *buf);

/*
 *  vcfbuf_ld() - find records with maximum LD values or the values in first record that exceeds thresholds
 *                set by vcfbuf_set_opt(..,LD_MAX*,..)
 *
 *  Returns 0 on success or -1 if no values were filled.
 *
 *  @val:  will be filled with the values
 *          .. correlation coefficient r-squared
 *          .. Lewontin's D' (PMID: 19433632)
 *          .. Ragsdale's \hat{D} (doi:10.1093/molbev/msz265)
 *  @rec: corresponding positions or NULL if the value(s) has not been set
 */
#define VCFBUF_LD_N 3
#define VCFBUF_LD_IDX_R2 0
#define VCFBUF_LD_IDX_LD 1
#define VCFBUF_LD_IDX_HD 2
typedef struct
{
    double val[VCFBUF_LD_N];    // r2, ld, hd
    bcf1_t *rec[VCFBUF_LD_N];   // record with max r2, ld, hd
}
vcfbuf_ld_t;
int vcfbuf_ld(vcfbuf_t *buf, bcf1_t *rec, vcfbuf_ld_t *ld);

#endif

