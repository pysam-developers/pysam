/* The MIT License

   Copyright (c) 2016 Genome Research Ltd.

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
#include "bcftools.h"
#include "bin.h"

struct _bin_t
{
    float *bins;
    int nbins;
};

bin_t *bin_init(const char *list_def, float min, float max)
{
    bin_t *bin = (bin_t*) calloc(1,sizeof(bin_t));

    // a comma indicates a list, otherwise a file
    int is_file = strchr(list_def,',') ? 0 : 1;
    int i, nlist;
    char **list = hts_readlist(list_def, is_file, &nlist);
    bin->nbins = nlist;
    bin->bins  = (float*) malloc(sizeof(float)*nlist);
    for (i=0; i<nlist; i++)
    {
        char *tmp;
        bin->bins[i] = strtod(list[i],&tmp);
        if ( !tmp ) error("Could not parse %s: %s\n", list_def, list[i]);
        if ( min!=max && (bin->bins[i]<min || bin->bins[i]>max) )
            error("Expected values from the interval [%f,%f], found %s\n", list[i]);
        free(list[i]); 
    }
    free(list);

    if ( min!=max )
    {
        // make sure we've got both boundaries: min,max.
        assert( nlist>1 );
        float max_err = (bin->bins[1] - bin->bins[0])*1e-6;
        if ( fabs(bin->bins[0] - min) > max_err )
        {
            bin->bins = (float*) realloc(bin->bins, (++bin->nbins)*sizeof(float));
            memmove(bin->bins+1, bin->bins, sizeof(float)*(bin->nbins-1));
            bin->bins[0] = min;
        }
        if ( fabs(bin->bins[bin->nbins-1] - max) > max_err ) 
        {
            bin->bins = (float*) realloc(bin->bins, (++bin->nbins)*sizeof(float));
            bin->bins[bin->nbins-1] = max;
        }
    }
    return bin;
}

void bin_destroy(bin_t *bin)
{
    free(bin->bins);
    free(bin);
}

int bin_get_size(bin_t *bin) { return bin->nbins; }

float bin_get_value(bin_t *bin, int idx) { return bin->bins[idx]; }

int bin_get_idx(bin_t *bin, float value)
{
    if ( bin->bins[bin->nbins-1] < value ) return bin->nbins-1;

    // Binary search in half-closed,half-open intervals [)
    int imin = 0, imax = bin->nbins - 2;
    while ( imin<imax )
    {
        int i = (imin+imax)/2;
        if ( value < bin->bins[i] ) imax = i - 1;
        else if ( value > bin->bins[i] ) imin = i + 1;
        else return i;
    }
    if ( bin->bins[imax] <= value ) return imax;
    return imin - 1;
}

