/* 
    Copyright (C) 2016 Genome Research Ltd.

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
    Parse --samples and --samples-file
*/

#ifndef __SMPL_ILIST_H__
#define __SMPL_ILIST_H__

#include <htslib/vcf.h>

#define SMPL_NONE   0   // flexible error recovery
#define SMPL_STRICT 1   // samples must exist

typedef struct
{
    int *idx;  // index to bcf_hdr_t.samples 
    int n;
}
smpl_ilist_t;

smpl_ilist_t *smpl_ilist_init(bcf_hdr_t *hdr, char *sample_list, int is_file, int flags);
smpl_ilist_t *smpl_ilist_map(bcf_hdr_t *hdr_a, bcf_hdr_t *hdr_b, int flags);
void smpl_ilist_destroy(smpl_ilist_t *smpl);

#endif
