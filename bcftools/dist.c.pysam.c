#include "bcftools.pysam.h"

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "dist.h"

extern void error(const char *format, ...);

struct _dist_t
{
    uint64_t *bins, nvalues;
    int nbins;
    int npow;   // the number of orders of magnitude to represent exactly
    int nexact; // pow(10,npow)
    int nlevel;
};

dist_t *dist_init(int npow)
{
    dist_t *dist = (dist_t*) calloc(1,sizeof(dist_t));
    dist->npow   = npow;
    dist->nexact = pow(10,npow);
    dist->nlevel = dist->nexact - pow(10,npow-1);
    return dist;
}

void dist_destroy(dist_t *dist)
{
    if ( !dist ) return;
    free(dist->bins);
    free(dist);
}

int dist_nbins(dist_t *dist)
{
    return dist->nbins;
}

int dist_nvalues(dist_t *dist)
{
    return dist->nvalues;
}

uint32_t dist_insert(dist_t *dist, uint32_t value)
{
    int ibin;

    if ( value <= dist->nexact ) 
        ibin = value;
    else
    {
        int npow  = (int) log10(value);
        int level = npow - dist->npow + 1;
        uint32_t step = pow(10, level);
        ibin = dist->nexact + dist->nlevel*(level-1) + (value - pow(10,npow)) / step;
    }

    if ( ibin >= dist->nbins )
    {
        dist->bins = (uint64_t*) realloc(dist->bins, sizeof(*dist->bins)*(ibin+1));
        memset(dist->bins + dist->nbins, 0, (ibin+1 - dist->nbins)*sizeof(*dist->bins));
        dist->nbins = ibin+1;
    }
    dist->bins[ibin]++;
    dist->nvalues++;
    return ibin;
}
uint32_t dist_insert_n(dist_t *dist, uint32_t value, uint32_t cnt)
{
    if ( !cnt ) return 0;
    int ibin = dist_insert(dist, value);
    dist->bins[ibin] += cnt - 1;
    dist->nvalues += cnt;
    return ibin;
}

uint64_t dist_get(dist_t *dist, uint32_t idx, uint32_t *beg, uint32_t *end)
{
    if ( idx < dist->nexact )
    {
        if ( beg ) *beg = idx;
        if ( end ) *end = idx + 1;
    }
    else
    {
        int level = (idx - dist->nexact) / dist->nlevel + 1;
        int bin   = idx - dist->nexact - dist->nlevel*(level-1);

        uint32_t step  = pow(10, level);
        uint32_t value = pow(10, level + dist->npow - 1) + step*bin;

        if ( beg ) *beg = value;
        if ( end ) *end = value + step;
    }
    return dist->bins[idx];
}

