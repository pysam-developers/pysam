/*  read_consensus.h -- create and maintain consensus of reads

    Copyright (C) 2022 Genome Research Ltd.

    Author: pd3@sanger

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

#ifndef READ_CONSENSUS_H
#define READ_CONSENSUS_H

#include <stdint.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#ifndef DEBUG_RCNS
#define DEBUG_RCNS 0
#endif

typedef struct
{
    char *seq;          // nt5 sequence: "ACGTN"[(int)seq[i]]
    int nseq, ipos;     // the sequence length and the `pos` index relative to seq
}
cns_seq_t;

typedef struct _read_cns_t read_cns_t;

// Init and destroy read consensus
read_cns_t *rcns_init(hts_pos_t pos, hts_pos_t beg, hts_pos_t end);
void rcns_destroy(read_cns_t *rcns);

// Reset the structures for new sample and/or position
int rcns_reset(read_cns_t *rcns, hts_pos_t pos, hts_pos_t beg, hts_pos_t end);

// Add reads to consensus. The provided structures must continue to exist
// until rcns_get_consensus() is called.
//
// Todo (easy): allow it to be called once or multiple times, eg for
// creating a shared consensus for multiple samples
int rcns_set_reads(read_cns_t *rcns, bam_pileup1_t *plp, int nplp);

// Generate up to two consensus sequences, cns_seq[1].nseq is 0 when only
// the first is set
cns_seq_t *rcns_get_consensus(read_cns_t *rcns, const char *ref);

#endif
