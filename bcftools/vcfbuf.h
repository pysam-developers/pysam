/* The MIT License

   Copyright (c) 2017-2024 Genome Research Ltd.

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
    VCFBUF_DUMMY,           // int {0,1}, the caller maintains the buffer via push/peek/flush, nothing is removed by vcfbuf

    // pruning
    PRUNE_NSITES,           // int, leave max this many sites in the window
    PRUNE_NSITES_MODE,      // char *, maxAF (keep sites with max AF), 1st (sites that come first), rand (pick randomly)
    PRUNE_AF_TAG,           // char *, use this INFO/AF tag with VCFBUF_NSITES

    // duplicates and overlaps
    MARK,                   // w: char *, resolve overlaps by preferentially removing sites according to EXPR:
                            //      min(QUAL) .. remove sites with lowest QUAL until overlaps are resolved
                            //      overlap   .. select all overlapping sites
                            //      dup       .. select duplicate sites
                            // r: use as
                            //      while ( (rec=vcfbuf_flush(buf,flush_all)) )
                            //      {
                            //          int is_marked = vcfbuf_get_val(buf,int,MARK);
                            //          if ( is_marked ) do_something(rec);
                            //      }
    MARK_MISSING_EXPR,      // char *, what to do when missing value are encountered with min(QUAL)
                            //      0   .. set to 0 (the default)
                            //      DP  .. scale max quality in the window proportionally to INFO/DP

    // LD related options
    LD_RAND_MISSING,        // randomize rather than ignore missing genotypes
    LD_FILTER1,             // exclude the next record inserted by vcfbuf_push() from LD analysis
    LD_MAX_R2,              // If set, vcfbuf_ld() will stop at the first record that exceeds the R2,
    LD_MAX_LD,              //      LD, or HD threshold. When multiple are set, the OR logic is applied
    LD_MAX_HD,              //
}
vcfbuf_opt_t;


/**
 *  vcfbuf_set() - set various options, see the vcfbuf_opt_t keys for the complete list
 *
 *  Returns 0 if the call succeeded, or negative number on error.
 */
int vcfbuf_set(vcfbuf_t *buf, vcfbuf_opt_t key, ...);   // returns 0 on success

/**
 *  vcfbuf_get()     - get various options, see the vcfbuf_opt_t keys
 *  vcfbuf_get_val() - wrapper for `vcfbuf_get()` to return typed value
 *
 *  The former returns pointer to the memory area populated by the requested setting,
 *  its type can be inferred from the vcfbuf_opt_t documentation.
 */
void *vcfbuf_get(vcfbuf_t *buf, vcfbuf_opt_t key, ...);
#define vcfbuf_get_val(buf,type,key) (*(type*)vcfbuf_get(buf, key))


/*
 *  vcfbuf_init() - init buffer
 *  @win:   number of sites (>0), bp (<0)
 */
vcfbuf_t *vcfbuf_init(bcf_hdr_t *hdr, int win);
void vcfbuf_destroy(vcfbuf_t *buf);

/*
 *  vcfbuf_push() - push a new site for analysis.
 *
 *  Note that vcfbuf_flush() or vcfbuf_peek() must be called before next site is pushed.
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

/*
 *  vcfbuf_flush() - returns the next record or NULL, depending on the mode of operation and
 *      the content of the buffer
 *
 *  @flush_all: 1 if no more vcfbuf_push() calls will follow, 0 otherwise
 */
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

