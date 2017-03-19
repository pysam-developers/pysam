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

/*
    Simple binning of float values into predefined bins
*/

#ifndef __BIN_H__
#define __BIN_H__

#include <stdio.h>

typedef struct _bin_t bin_t;

/*
 *  bin_init() - init bins
 *  @list: list of half-open intervals [). If the list does not contain commas,
 *      it is interpreted as a file name.
 *  @min,max:  extreme values. This is for user convenience so that well-known
 *      extremes can be left out from the list. Ignored if min=max
 */
bin_t *bin_init(const char *list, float min, float max);
void bin_destroy(bin_t *bin);

/*
 *  bin_get_size() - number of boundaries, subtract 1 to get the number of bins
 */
int bin_get_size(bin_t *bin);

/*
   bin_get_idx() - find the bin index which corresponds to the value (binary search)
   Returns the bin index 0 <= idx <= size-2 or -1,size-1 for out of range values.
 */
int bin_get_idx(bin_t *bin, float value);

/*
   bin_get_value() - get the i-th boundary value, i=0,..,size-1
 */
float bin_get_value(bin_t *bin, int ith);

#endif

