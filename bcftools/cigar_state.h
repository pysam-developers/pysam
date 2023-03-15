/*  cigar_state.h -- API for efficient parsing of CIGAR strings

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

#ifndef CIGAR_STATE_H
#define CIGAR_STATE_H

#include <stdint.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

typedef struct
{
    bam1_t *bam;
    uint32_t *cigar;
    uint8_t *seq;
    int ncig;
    int icig;           // position in the cigar string
    int iseq;           // the cigar[icigar] operation refers to &seq[iseq]
    hts_pos_t ref_pos;  // reference coordinate, corresponds to iseq; points to
                        // the first base after the read when consumed
}
cigar_state_t;

static inline void cstate_init(cigar_state_t *cs, bam1_t *bam)
{
    cs->bam   = bam;
    cs->cigar = bam_get_cigar(bam);
    cs->seq   = bam_get_seq(bam);
    cs->ncig  = bam->core.n_cigar;
    cs->icig  = 0;
    cs->iseq  = 0;
    cs->ref_pos = bam->core.pos;
}

/**
 *  cstate_seek_fwd() - Move in the cigar forward to find query index that
 *  matches the reference position.
 *
 *  When the position is not contained within the sequence, either because there
 *  is a deletion or there is no overlap, the behavior is controlled by the value
 *  of trim_left:
 *      - read starts after: qry_beg > pos && trim_left=1 .. returns 0 and sets pos to qry_beg
 *      - read starts after: qry_beg > pos && trim_left=0 .. returns -1
 *      - read ends before: qry_end < pos && trim_left=1 .. returns -2
 *      - read ends before: qry_end < pos && trim_left=0 .. returns qry_len-1 and sets pos to qry_end
 *      - pos inside a deletion && trim_left=1 .. returns position after the deletion
 *      - pos inside a deletion && trim_left=0 .. returns position before the deletion
 */
static inline int cstate_seek_fwd(cigar_state_t *cs, hts_pos_t *pos_ptr, int trim_left)
{
    hts_pos_t pos = *pos_ptr;
    while ( cs->ref_pos <= pos )
    {
        if ( cs->icig >= cs->ncig )     // the read ends before pos
        {
            if ( trim_left ) return -2;
            *pos_ptr = cs->ref_pos - 1;
            return cs->iseq - 1;
        }

        int op  = cs->cigar[cs->icig] &  BAM_CIGAR_MASK;
        int len = cs->cigar[cs->icig] >> BAM_CIGAR_SHIFT;
        if ( op==BAM_CMATCH || op==BAM_CEQUAL || op==BAM_CDIFF )
        {
            if ( cs->ref_pos + len > pos ) return pos - cs->ref_pos + cs->iseq;  // the cigar op overlaps pos
            cs->ref_pos += len;
            cs->iseq += len;
            cs->icig++;
            continue;
        }
        if ( op==BAM_CINS || op==BAM_CSOFT_CLIP )
        {
            cs->iseq += len;
            cs->icig++;
            continue;
        }
        if ( op==BAM_CDEL || op==BAM_CREF_SKIP )
        {
            if ( cs->ref_pos + len > pos )
            {
                // The deletion overlaps the position. NB: assuming del is never the first or last op
                *pos_ptr = trim_left ? cs->ref_pos + len : cs->ref_pos - 1;
                return trim_left ? cs->iseq : cs->iseq - 1;
            }
            cs->ref_pos += len;
            cs->icig++;
            continue;
        }
    }
    // the read starts after pos
    if ( trim_left )
    {
        *pos_ptr = cs->bam->core.pos;
        return 0;
    }
    return -1;
}


/**
 *  cstate_seek_op_fwd() - Move in the cigar forward to find query index that
 *  matches the seek operator and the reference position.
 *
 *  In order to match a deletion, pass the position of the first deleted base.
 *  In order to match an insertion, pass the reference coordinate of the base
 *  after the inserted sequence.
 *
 *  Returns the index to the query sequence cs->seq
 *  on success; -1 when there is no such matching position but the cigar
 *  is still not entirely consumed (e.g. a deletion or a soft-clip); -2
 *  when there is no overlap (i.e. the read ends before the position).
 */
static inline int cstate_seek_op_fwd(cigar_state_t *cs, hts_pos_t pos, int seek_op, int *oplen)
{
    while ( cs->ref_pos <= pos )
    {
        if ( cs->icig >= cs->ncig ) return -2;

        int op  = cs->cigar[cs->icig] &  BAM_CIGAR_MASK;
        int len = cs->cigar[cs->icig] >> BAM_CIGAR_SHIFT;
        if ( op==BAM_CMATCH || op==BAM_CEQUAL || op==BAM_CDIFF )
        {
            if ( cs->ref_pos + len <= pos )
            {
                cs->ref_pos += len;
                cs->iseq += len;
                cs->icig++;
                continue;
            }
            if ( seek_op==BAM_CMATCH ) return pos - cs->ref_pos + cs->iseq;
            return -1;
        }
        if ( op==BAM_CINS || op==BAM_CSOFT_CLIP )
        {
            if ( cs->ref_pos == pos && seek_op==op )
            {
                if ( oplen ) *oplen = len;
                return cs->iseq;
            }
            if ( cs->ref_pos >= pos ) return -1;
            cs->iseq += len;
            cs->icig++;
            continue;
        }
        if ( op==BAM_CDEL || op==BAM_CREF_SKIP )
        {
            if ( cs->ref_pos == pos && seek_op==op )
            {
                if ( oplen ) *oplen = len;
                return cs->iseq;
            }
            if ( cs->ref_pos >= pos ) return -1;
            cs->ref_pos += len;
            cs->icig++;
            continue;
        }
    }
    return cs->icig < cs->ncig ? -1 : -2;
}

#endif
