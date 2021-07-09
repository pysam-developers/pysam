/* The MIT License

   Copyright (c) 2016-2020 Genome Research Ltd.

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
    Logarithmic binning

    Example of usage:

        // Initialize, make the binning exact up to 10^4, then add a log-step
        dist_t *dist = dist_init(4);

        // Insert values
        int i;
        for (i=0; i<1e6; i++)
            dist_insert(dist, i);

        // Number of bins used
        int n = dist_n(dist);

        // Now print the distribution
        uint32_t beg, end;
        for (i=0; i<n; i++)
        {
            // Raw count in the bin. The boundaries beg,end are optional, 
            // and can be used to plot correctly the density
            uint64_t cnt = dist_get(dist, i, &beg, &end);
            if ( !cnt ) continue;

            // Print the interval, count and density
            printf("%u\t%u\t%"PRIu64"\t%f\n", beg, end, cnt, (double)cnt/(end-beg));
        }

        // Clean up
        dist_destroy(dist);
 */

#ifndef __DIST_H__
#define __DIST_H__

#include <stdio.h>
#include <inttypes.h>

typedef struct _dist_t dist_t;

/*
 *  dist_init() - init bins
 */
dist_t *dist_init(int npow);
void dist_destroy(dist_t *dist);

/*
    dist_nbins() - get the number of bins
 */
int dist_nbins(dist_t *dist);

/*
    dist_nvalues() - get the total number of values inserted
 */
int dist_nvalues(dist_t *dist);

/*
    dist_insert()   - insert new value
    dist_insert_n() - insert new value n times
 */
uint32_t dist_insert(dist_t *dist, uint32_t value);
uint32_t dist_insert_n(dist_t *dist, uint32_t value, uint32_t cnt);

/*
   dist_get() 
   @idx:        from the interval [0,dist_n-1]
   @beg,end:    [beg,end)
 */
uint64_t dist_get(dist_t *dist, uint32_t idx, uint32_t *beg, uint32_t *end);

#endif

