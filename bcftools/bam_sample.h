/*  bam_sample.h -- group data by sample.

    Copyright (C) 2010 Broad Institute.
    Copyright (C) 2016 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>, Petr Danecek <pd3@sanger.ac.uk>

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef BAM_SAMPLE_H
#define BAM_SAMPLE_H

#include <htslib/sam.h>

typedef struct _bam_smpl_t bam_smpl_t;

bam_smpl_t *bam_smpl_init(void);

int bam_smpl_add_samples(bam_smpl_t *bsmpl, char *list, int is_file);
int bam_smpl_add_readgroups(bam_smpl_t *bsmpl, char *list, int is_file);
void bam_smpl_ignore_readgroups(bam_smpl_t* bsmpl);

// The above should be called only before bams are added. Returns the BAM id
// to be passed to bam_smpl_get_sample_id() later. It is safe to assume
// sequential numbering, starting from 0.
//
int bam_smpl_add_bam(bam_smpl_t *bsmpl, char *bam_hdr, const char *fname);

const char **bam_smpl_get_samples(bam_smpl_t *bsmpl, int *nsmpl);
int bam_smpl_get_sample_id(bam_smpl_t *bsmpl, int bam_id, bam1_t *bam_rec);

void bam_smpl_destroy(bam_smpl_t *bsmpl);

#endif
