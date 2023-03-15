/*
 * Copyright (c) 2019-2022 Genome Research Ltd.
 * Author(s): James Bonfield
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *       Institute nor the names of its contributors may be used to endorse
 *       or promote products derived from this software without specific
 *       prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// As per standard rANS_static but using optional RLE or bit-packing
// techniques prior to entropy encoding.  This is a significant
// reduction in some data sets.

// top bits in order byte
#define X_PACK   0x80    // Pack 2,4,8 or infinite symbols into a byte.
#define X_RLE    0x40    // Run length encoding with runs & lits encoded separately
#define X_CAT    0x20    // Nop; for tiny segments where rANS overhead is too big
#define X_NOSZ   0x10    // Don't store the original size; used by STRIPE mode
#define X_STRIPE 0x08    // For 4-byte integer data; rotate & encode 4 streams.
#define X_EXT    0x04    // External compression codec via magic num (gz, xz, bz2)
#define X_ORDER  0x03    // Mask to obtain order

#include "config.h"

#ifdef HAVE_LIBBZ2
#include <bzlib.h>
#endif

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <limits.h>

#include "arith_dynamic.h"
#include "varint.h"
#include "pack.h"
#include "utils.h"

#define MIN(a,b) ((a)<(b)?(a):(b))

/*-----------------------------------------------------------------------------
 * Memory to memory compression functions.
 *
 * These are original versions without any manual loop unrolling. They
 * are easier to understand, but can be up to 2x slower.
 */
#define MAGIC 8

unsigned int arith_compress_bound(unsigned int size, int order) {
    int N = (order>>8) & 0xff;
    if (!N) N=4;
    return (order == 0
        ? 1.05*size + 257*3 + 4
        : 1.05*size + 257*257*3 + 4 + 257*3+4) + 5 +
        ((order & X_PACK) ? 1 : 0) +
        ((order & X_RLE) ? 1 + 257*3+4: 0) +
        ((order & X_STRIPE) ? 7 + 5*N: 0);
}

#ifndef MODEL_256 // see fqzcomp_qual_fuzz.c
#define NSYM 256
#include "c_simple_model.h"
#endif

// Compresses in_size bytes from 'in' to *out_size bytes in 'out'.
//
// NB: The output buffer does not hold the original size, so it is up to
// the caller to store this.
static
unsigned char *arith_compress_O0(unsigned char *in, unsigned int in_size,
                                 unsigned char *out, unsigned int *out_size) {
    int i, bound = arith_compress_bound(in_size,0)-5; // -5 for order/size

    if (!out) {
        *out_size = bound;
        out = malloc(*out_size);
    }
    if (!out || bound > *out_size)
        return NULL;

    unsigned int m = 0;
    for (i = 0; i < in_size; i++)
        if (m < in[i])
            m = in[i];
    m++;
    *out = m;

    SIMPLE_MODEL(256,_) byte_model;
    SIMPLE_MODEL(256,_init)(&byte_model, m);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out+1);
    RC_StartEncode(&rc);

    for (i = 0; i < in_size; i++)
        SIMPLE_MODEL(256, _encodeSymbol)(&byte_model, &rc, in[i]);

    RC_FinishEncode(&rc);

    // Finalise block size and return it
    *out_size = RC_OutSize(&rc)+1;

    return out;
}

static
unsigned char *arith_uncompress_O0(unsigned char *in, unsigned int in_size,
                                   unsigned char *out, unsigned int out_sz) {
    RangeCoder rc;
    int i;
    unsigned int m = in[0] ? in[0] : 256;

    SIMPLE_MODEL(256,_) byte_model;
    SIMPLE_MODEL(256,_init)(&byte_model, m);

    if (!out)
        out = malloc(out_sz);
    if (!out)
        return NULL;

    RC_SetInput(&rc, (char *)in+1, (char *)in+in_size);
    RC_StartDecode(&rc);

    for (i = 0; i < out_sz; i++)
        out[i] = SIMPLE_MODEL(256, _decodeSymbol)(&byte_model, &rc);

    RC_FinishDecode(&rc);
    
    return out;
}


//-----------------------------------------------------------------------------
static
unsigned char *arith_compress_O1(unsigned char *in, unsigned int in_size,
                                 unsigned char *out, unsigned int *out_size) {
    int i, bound = arith_compress_bound(in_size,0)-5; // -5 for order/size
    unsigned char *out_free = NULL;

    if (!out) {
        *out_size = bound;
        out_free = out = malloc(*out_size);
    }
    if (!out || bound > *out_size)
        return NULL;

    SIMPLE_MODEL(256,_) *byte_model =
        htscodecs_tls_alloc(256 * sizeof(*byte_model));
    if (!byte_model) {
        free(out_free);
        return NULL;
    }
    unsigned int m = 0;
    if (1 || in_size > 1000) {
        for (i = 0; i < in_size; i++)
            if (m < in[i])
                m = in[i];
        //fprintf(stderr, "%d max %d\n", in_size, m);
        m++;
    }
    *out = m;
    for (i = 0; i < 256; i++)
        SIMPLE_MODEL(256,_init)(&byte_model[i], m);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out+1);
    RC_StartEncode(&rc);

    uint8_t last = 0;
    for (i = 0; i < in_size; i++) {
        SIMPLE_MODEL(256, _encodeSymbol)(&byte_model[last], &rc, in[i]);
        last = in[i];
    }

    RC_FinishEncode(&rc);

    // Finalise block size and return it
    *out_size = RC_OutSize(&rc)+1;

    htscodecs_tls_free(byte_model);
    return out;
}

static
unsigned char *arith_uncompress_O1(unsigned char *in, unsigned int in_size,
                                   unsigned char *out, unsigned int out_sz) {
    RangeCoder rc;
    unsigned char *out_free = NULL;

    if (!out)
        out_free = out = malloc(out_sz);
    if (!out)
        return NULL;


    SIMPLE_MODEL(256,_) *byte_model =
        htscodecs_tls_alloc(256 * sizeof(*byte_model));
    if (!byte_model) {
        free(out_free);
        return NULL;
    }

    unsigned int m = in[0] ? in[0] : 256, i;
    for (i = 0; i < 256; i++)
        SIMPLE_MODEL(256,_init)(&byte_model[i], m);

    RC_SetInput(&rc, (char *)in+1, (char *)in+in_size);
    RC_StartDecode(&rc);

    unsigned char last = 0;
    for (i = 0; i < out_sz; i++) {
        out[i] = SIMPLE_MODEL(256, _decodeSymbol)(&byte_model[last], &rc);
        last = out[i];
    }

    RC_FinishDecode(&rc);
    
    htscodecs_tls_free(byte_model);
    return out;
}

//-----------------------------------------------------------------------------

// Disable O2 for now
#if 0

#if 0
unsigned char *arith_compress_O2(unsigned char *in, unsigned int in_size,
                                 unsigned char *out, unsigned int *out_size) {
    fprintf(stderr, "WARNING: using undocumented O2 arith\n");

    int i, j;
    int bound = arith_compress_bound(in_size,0)-5; // -5 for order/size

    if (!out) {
        *out_size = bound;
        out = malloc(*out_size);
    }
    if (!out || bound > *out_size)
        return NULL;

    unsigned int m = 0;
    if (1 || in_size > 1000) {
        for (i = 0; i < in_size; i++)
            if (m < in[i])
                m = in[i];
        //fprintf(stderr, "%d max %d\n", in_size, m);
        m++;
    }
    *out = m;

    SIMPLE_MODEL(256,_) *byte_model;
    byte_model = malloc(256*256*sizeof(*byte_model));
    for (i = 0; i < 256; i++)
        for (j = 0; j < 256; j++)
            SIMPLE_MODEL(256,_init)(&byte_model[i*256+j], m);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out+1);
    RC_StartEncode(&rc);

    unsigned char last1 = 0, last2 = 0;
    for (i = 0; i < in_size; i++) {
        SIMPLE_MODEL(256, _encodeSymbol)(&byte_model[last1*256 + last2], &rc, in[i]);
        last2 = last1;
        last1 = in[i];
    }

    free(byte_model);
    RC_FinishEncode(&rc);

    // Finalise block size and return it
    *out_size = RC_OutSize(&rc)+1;

    return out;
}
#else
unsigned char *arith_compress_O2(unsigned char *in, unsigned int in_size,
                                 unsigned char *out, unsigned int *out_size) {
    fprintf(stderr, "WARNING: using undocumented O2 arith\n");

    int i, j;
    int bound = arith_compress_bound(in_size,0)-5; // -5 for order/size

    if (!out) {
        *out_size = bound;
        out = malloc(*out_size);
    }
    if (!out || bound > *out_size)
        return NULL;

    unsigned int m = 0;
    if (1 || in_size > 1000) {
        for (i = 0; i < in_size; i++)
            if (m < in[i])
                m = in[i];
        //fprintf(stderr, "%d max %d\n", in_size, m);
        m++;
    }
    *out = m;

    SIMPLE_MODEL(256,_) *byte_model;
    byte_model = malloc(256*256*sizeof(*byte_model));
    for (i = 0; i < 256; i++)
        for (j = 0; j < 256; j++)
            SIMPLE_MODEL(256,_init)(&byte_model[i*256+j], m);
    SIMPLE_MODEL(256,_) byte_model1[256];
    for (i = 0; i < 256; i++)
        SIMPLE_MODEL(256,_init)(&byte_model1[i], m);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out+1);
    RC_StartEncode(&rc);

    unsigned char last1 = 0, last2 = 0;
    for (i = 0; i < in_size; i++) {
        // Use Order-1 is order-2 isn't sufficiently advanced yet (75+ symbols)
        if (byte_model[last1*256+last2].TotFreq <= m+75*16) {
            SIMPLE_MODEL(256, _encodeSymbol)(&byte_model1[last1], &rc, in[i]);
            SIMPLE_MODEL(256, _updateSymbol)(&byte_model[last1*256 + last2], &rc, in[i]);
        } else {
            SIMPLE_MODEL(256, _encodeSymbol)(&byte_model[last1*256 + last2], &rc, in[i]);
            //SIMPLE_MODEL(256, _updateSymbol)(&byte_model1[last1], &rc, in[i]);
        }
        last2 = last1;
        last1 = in[i];
    }

    free(byte_model);
    RC_FinishEncode(&rc);

    // Finalise block size and return it
    *out_size = RC_OutSize(&rc)+1;

    return out;
}
#endif

unsigned char *arith_uncompress_O2(unsigned char *in, unsigned int in_size,
                                   unsigned char *out, unsigned int out_sz) {
    RangeCoder rc;

    SIMPLE_MODEL(256,_) *byte_model;
    byte_model = malloc(256*256*sizeof(*byte_model));
    unsigned int m = in[0] ? in[0] : 256, i, j;
    for (i = 0; i < 256; i++)
        for (j = 0; j < 256; j++)
            SIMPLE_MODEL(256,_init)(&byte_model[i*256+j], m);
    
    if (!out)
        out = malloc(out_sz);
    if (!out)
        return NULL;

    RC_SetInput(&rc, (char *)in+1, (char *)in+in_size);
    RC_StartDecode(&rc);

    unsigned char last1 = 0, last2 = 0;
    for (i = 0; i < out_sz; i++) {
        out[i] = SIMPLE_MODEL(256, _decodeSymbol)(&byte_model[last1*256 + last2], &rc);
        last2 = last1;
        last1 = out[i];
    }

    free(byte_model);
    RC_FinishDecode(&rc);
    
    return out;
}

#endif // Disable O2
/*-----------------------------------------------------------------------------
 */

#undef NSYM
#define NSYM 258
#include "c_simple_model.h"
#define MAX_RUN 4

static
unsigned char *arith_compress_O0_RLE(unsigned char *in, unsigned int in_size,
                                     unsigned char *out, unsigned int *out_size) {
    int i, bound = arith_compress_bound(in_size,0)-5; // -5 for order/size
    unsigned char *out_free = NULL;

    if (!out) {
        *out_size = bound;
        out_free = out = malloc(*out_size);
    }
    if (!out || bound > *out_size)
        return NULL;

    unsigned int m = 0;
    for (i = 0; i < in_size; i++)
        if (m < in[i])
            m = in[i];
    m++;
    *out = m;

    SIMPLE_MODEL(256,_) byte_model;
    SIMPLE_MODEL(256,_init)(&byte_model, m);

    SIMPLE_MODEL(NSYM,_) *run_model =
        htscodecs_tls_alloc(NSYM * sizeof(*run_model));
    if (!run_model) {
        free(out_free);
        return NULL;
    }

    for (i = 0; i < NSYM; i++)
        SIMPLE_MODEL(NSYM,_init)(&run_model[i], MAX_RUN);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out+1);
    RC_StartEncode(&rc);

    unsigned char last = 0;
    for (i = 0; i < in_size;) {
        //SIMPLE_MODEL(256, _encodeSymbol)(&byte_model, &rc, in[i]);
        SIMPLE_MODEL(256, _encodeSymbol)(&byte_model, &rc, in[i]);
        //fprintf(stderr, "lit %c (ctx %c)\n", in[i], last);
        int run = 0;
        last = in[i++];
        while (i < in_size && in[i] == last/* && run < MAX_RUN-1*/)
            run++, i++;
        int rctx = last;
        do {
            int c = run < MAX_RUN ? run : MAX_RUN-1;
            SIMPLE_MODEL(NSYM, _encodeSymbol)(&run_model[rctx], &rc, c);
            run -= c;

            if (rctx == last)
                rctx = 256;
            else
                rctx += (rctx < NSYM-1);
            if (c == MAX_RUN-1 && run == 0)
                SIMPLE_MODEL(NSYM, _encodeSymbol)(&run_model[rctx], &rc, 0);
        } while (run);
    }

    RC_FinishEncode(&rc);

    // Finalise block size and return it
    *out_size = RC_OutSize(&rc)+1;

    //fprintf(stderr, "RLE %d to %d\n", in_size, *out_size);

    htscodecs_tls_free(run_model);
    return out;
}

static
unsigned char *arith_uncompress_O0_RLE(unsigned char *in, unsigned int in_size,
                                       unsigned char *out, unsigned int out_sz) {
    RangeCoder rc;
    int i;
    unsigned int m = in[0] ? in[0] : 256;
    unsigned char *out_free = NULL;

    if (!out)
        out_free = out = malloc(out_sz);
    if (!out)
        return NULL;

    SIMPLE_MODEL(256,_) byte_model;
    SIMPLE_MODEL(256,_init)(&byte_model, m);

    SIMPLE_MODEL(NSYM,_) *run_model =
        htscodecs_tls_alloc(NSYM * sizeof(*run_model));
    if (!run_model) {
        free(out_free);
        return NULL;
    }

    for (i = 0; i < NSYM; i++)
        SIMPLE_MODEL(NSYM,_init)(&run_model[i], MAX_RUN);

    RC_SetInput(&rc, (char *)in+1, (char *)in+in_size);
    RC_StartDecode(&rc);

    for (i = 0; i < out_sz; i++) {
        unsigned char last;
        last = out[i] = SIMPLE_MODEL(256, _decodeSymbol)(&byte_model, &rc);
        //fprintf(stderr, "lit %c\n", last);
        int run = 0, r = 0, rctx = out[i];
        do {
            r = SIMPLE_MODEL(NSYM, _decodeSymbol)(&run_model[rctx], &rc);
            if (rctx == last)
                rctx = 256;
            else
                rctx += (rctx < NSYM-1);
            //fprintf(stderr, "run %d (ctx %d, %d)\n", r, last, l);
            run += r;
        } while (r == MAX_RUN-1 && run < out_sz);
        while (run-- && i+1 < out_sz)
            out[++i] = last;
    }

    RC_FinishDecode(&rc);

    htscodecs_tls_free(run_model);
    return out;
}

static
unsigned char *arith_compress_O1_RLE(unsigned char *in, unsigned int in_size,
                                     unsigned char *out, unsigned int *out_size) {
    int i, bound = arith_compress_bound(in_size,0)-5; // -5 for order/size
    unsigned char *out_free = NULL;

    if (!out) {
        *out_size = bound;
        out_free = out = malloc(*out_size);
    }
    if (!out || bound > *out_size)
        return NULL;

    unsigned int m = 0;
    for (i = 0; i < in_size; i++)
        if (m < in[i])
            m = in[i];
    m++;
    *out = m;

    SIMPLE_MODEL(256,_) *byte_model =
        htscodecs_tls_alloc(256 * sizeof(*byte_model));
    if (!byte_model) {
        free(out_free);
        return NULL;
    }
    for (i = 0; i < 256; i++)
        SIMPLE_MODEL(256,_init)(&byte_model[i], m);

    SIMPLE_MODEL(NSYM,_) *run_model =
        htscodecs_tls_alloc(NSYM * sizeof(*run_model));
    if (!run_model) {
        htscodecs_tls_free(byte_model);
        free(out_free);
        return NULL;
    }
    for (i = 0; i < NSYM; i++)
        SIMPLE_MODEL(NSYM,_init)(&run_model[i], MAX_RUN);

    RangeCoder rc;
    RC_SetOutput(&rc, (char *)out+1);
    RC_StartEncode(&rc);

    unsigned char last = 0;
    for (i = 0; i < in_size;) {
        //SIMPLE_MODEL(256, _encodeSymbol)(&byte_model, &rc, in[i]);
        SIMPLE_MODEL(256, _encodeSymbol)(&byte_model[last], &rc, in[i]);
        //fprintf(stderr, "lit %c (ctx %c)\n", in[i], last);
        int run = 0;
        last = in[i++];
        while (i < in_size && in[i] == last/* && run < MAX_RUN-1*/)
            run++, i++;
        int rctx = last;
        do {
            int c = run < MAX_RUN ? run : MAX_RUN-1;
            SIMPLE_MODEL(NSYM, _encodeSymbol)(&run_model[rctx], &rc, c);
            run -= c;

            if (rctx == last)
                rctx = 256;
            else
                rctx += (rctx < NSYM-1);
            if (c == MAX_RUN-1 && run == 0)
                SIMPLE_MODEL(NSYM, _encodeSymbol)(&run_model[rctx], &rc, 0);
        } while (run);
    }

    RC_FinishEncode(&rc);

    // Finalise block size and return it
    *out_size = RC_OutSize(&rc)+1;

    //fprintf(stderr, "RLE %d to %d\n", in_size, *out_size);

    htscodecs_tls_free(byte_model);
    htscodecs_tls_free(run_model);
    return out;
}

static
unsigned char *arith_uncompress_O1_RLE(unsigned char *in, unsigned int in_size,
                                       unsigned char *out, unsigned int out_sz) {
    RangeCoder rc;
    int i;
    unsigned int m = in[0] ? in[0] : 256;
    unsigned char *out_free = NULL;

    if (!out)
        out_free = out = malloc(out_sz);
    if (!out)
        return NULL;

    SIMPLE_MODEL(256,_) *byte_model =
        htscodecs_tls_alloc(256 * sizeof(*byte_model));
    if (!byte_model) {
        free(out_free);
        return NULL;
    }
    for (i = 0; i < 256; i++)
        SIMPLE_MODEL(256,_init)(&byte_model[i], m);

    SIMPLE_MODEL(NSYM,_) *run_model =
        htscodecs_tls_alloc(NSYM * sizeof(*run_model));
    if (!run_model) {
        htscodecs_tls_free(byte_model);
        free(out_free);
        return NULL;
    }
    for (i = 0; i < NSYM; i++)
        SIMPLE_MODEL(NSYM,_init)(&run_model[i], MAX_RUN);

    RC_SetInput(&rc, (char *)in+1, (char *)in+in_size);
    RC_StartDecode(&rc);

    unsigned char last = 0;
    for (i = 0; i < out_sz; i++) {
        out[i] = SIMPLE_MODEL(256, _decodeSymbol)(&byte_model[last], &rc);
        //fprintf(stderr, "lit %c (ctx %c)\n", out[i], last);
        last = out[i];
        int run = 0, r = 0, rctx = last;

        do {
            r = SIMPLE_MODEL(NSYM, _decodeSymbol)(&run_model[rctx], &rc);
            if (rctx == last)
                rctx = 256;
            else
                rctx += (rctx < NSYM-1);
            run += r;
        } while (r == MAX_RUN-1 && run < out_sz);
        while (run-- && i+1 < out_sz)
            out[++i] = last;
    }

    RC_FinishDecode(&rc);

    htscodecs_tls_free(byte_model);
    htscodecs_tls_free(run_model);
    return out;
}

/*-----------------------------------------------------------------------------
 * Simple interface to the order-0 vs order-1 encoders and decoders.
 *
 * Smallest is method, <in_size> <input>, so worst case 2 bytes longer.
 */
unsigned char *arith_compress_to(unsigned char *in,  unsigned int in_size,
                                 unsigned char *out, unsigned int *out_size,
                                 int order) {
    unsigned int c_meta_len;
    uint8_t *rle = NULL, *packed = NULL;

    if (in_size > INT_MAX) {
        *out_size = 0;
        return NULL;
    }

    if (!out) {
        *out_size = arith_compress_bound(in_size, order);
        if (!(out = malloc(*out_size)))
            return NULL;
    }
    unsigned char *out_end = out + *out_size;

    if (in_size <= 20)
        order &= ~X_STRIPE;

    if (order & X_CAT) {
        out[0] = X_CAT;
        c_meta_len = 1 + var_put_u32(&out[1], out_end, in_size);
        memcpy(out+c_meta_len, in, in_size);
        *out_size = in_size+c_meta_len;
    }

    if (order & X_STRIPE) {
        int N = (order>>8);
        if (N == 0) N = 4; // default for compatibility with old tests

        if (N > 255)
            return NULL;

        unsigned char *transposed = malloc(in_size);
        unsigned int part_len[256];
        unsigned int idx[256];
        if (!transposed)
            return NULL;
        int i, j, x;

        for (i = 0; i < N; i++) {
            part_len[i] = in_size / N + ((in_size % N) > i);
            idx[i] = i ? idx[i-1] + part_len[i-1] : 0; // cumulative index
        }

        for (i = x = 0; i < in_size-N; i += N, x++) {
            for (j = 0; j < N; j++)
                transposed[idx[j]+x] = in[i+j];
        }
        for (; i < in_size; i += N, x++) {
            for (j = 0; i+j < in_size; j++)
                transposed[idx[j]+x] = in[i+j];
        }

        unsigned int olen2;
        unsigned char *out2, *out2_start;
        c_meta_len = 1;
        *out = order & ~X_NOSZ;
        c_meta_len += var_put_u32(out+c_meta_len, out_end, in_size);
        out[c_meta_len++] = N;

        out2_start = out2 = out+7+5*N; // shares a buffer with c_meta
        for (i = 0; i < N; i++) {
            // Brute force try all methods.
            // FIXME: optimise this bit.  Maybe learn over time?
            int j, best_j = 0, best_sz = INT_MAX;

            // Works OK with read names. The first byte is the most important,
            // as it has most variability (little-endian).  After that it's
            // often quite predictable.
            //
            // Do we gain in any other context in CRAM? Aux tags maybe?
            int m[][4] = {{3, 1,64,0},
                          {2, 1,0},
                          {2, 1,128},
                          {2, 1,128}};

//          int m[][6] = {{4, 1,64,2,0},  //test of adding in an order-2 codec
//                        {3, 1,2,0},
//                        {3, 1,2,128},
//                        {3, 1,2,128}};

// Other possibilities for methods to try.
//          int m[][10] = {{8, 1,128,129,64,65,192,193,4,0},
//                         {8, 1,128,129,64,65,192,193,4,0},
//                         {8, 1,128,129,64,65,192,193,4,0},
//                         {8, 1,128,129,64,65,192,193,4,0}};

//          int m[][9] = {{5, 1,128,64,65,0},
//                        {5, 1,128,64,65,0},
//                        {5, 1,128,64,65,0},
//                        {5, 1,128,64,65,0}};

//          int m[][6] = {{4, 0,1,128,64},
//                        {5, 0,1,128,65,193},
//                        {3, 0,1,128},
//                        {3, 0,1,128}};

//          int m[][6] = {{4, 1,128,64,0},
//                        {4, 1,128,65,0},
//                        {2, 128,0},
//                        {2, 128,0}};

//          int m[][6] = {{2, 64,0},
//                        {1, 0},
//                        {1, 128},
//                        {1, 128}};

//          int m[][6] = {{1, 0},
//                        {2, 128,0},
//                        {1, 128},
//                        {1, 128}};

            for (j = 1; j <= m[MIN(i,3)][0]; j++) {
                olen2 = *out_size - (out2 - out);
                //fprintf(stderr, "order=%d m=%d\n", order&3, m[MIN(i,4)][j]);
                if ((order&3) == 0 && (m[MIN(i,3)][j]&1))
                    continue;

                arith_compress_to(transposed+idx[i], part_len[i],
                                  out2, &olen2, m[MIN(i,3)][j] | X_NOSZ);
                if (best_sz > olen2) {
                    best_sz = olen2;
                    best_j = j;
                }
            }
//          if (best_j == 0) // none desireable
//              return NULL;
            if (best_j != j-1) {
                olen2 = *out_size - (out2 - out);
                arith_compress_to(transposed+idx[i], part_len[i],
                                  out2, &olen2, m[MIN(i,3)][best_j] | X_NOSZ);
            }
            out2 += olen2;
            c_meta_len += var_put_u32(out+c_meta_len, out_end, olen2);
        }
        memmove(out+c_meta_len, out2_start, out2-out2_start);
        free(transposed);
        *out_size = c_meta_len + out2-out2_start;
        return out;
    }

    int do_pack = order & X_PACK;
    int do_rle  = order & X_RLE;
    int no_size = order & X_NOSZ;
    int do_ext  = order & X_EXT;

    out[0] = order;
    c_meta_len = 1;

    if (!no_size)
        c_meta_len += var_put_u32(&out[1], out_end, in_size);

    order &= 0x3;

    // Format is compressed meta-data, compressed data.
    // Meta-data can be empty, pack, rle lengths, or pack + rle lengths.
    // Data is either the original data, bit-packed packed, rle literals or
    // packed + rle literals.

    if (do_pack && in_size) {
        // PACK 2, 4 or 8 symbols into one byte.
        int pmeta_len;
        uint64_t packed_len;
        packed = hts_pack(in, in_size, out+c_meta_len, &pmeta_len, &packed_len);
        if (!packed) {
            out[0] &= ~X_PACK;
            do_pack = 0;
            free(packed);
            packed = NULL;
        } else {
            in = packed;
            in_size = packed_len;
            c_meta_len += pmeta_len;

            // Could derive this rather than storing verbatim.
            // Orig size * 8/nbits (+1 if not multiple of 8/n)
            int sz = var_put_u32(out+c_meta_len, out_end, in_size);
            c_meta_len += sz;
            *out_size -= sz;
        }
    } else if (do_pack) {
        out[0] &= ~X_PACK;
    }

    if (do_rle && !in_size) {
        out[0] &= ~X_RLE;
    }

    *out_size -= c_meta_len;
    if (order && in_size < 8) {
        out[0] &= ~3;
        order  &= ~3;
    }

    if (do_ext) {
        // Use an external compression library instead.
        // For now, bzip2
#ifdef HAVE_LIBBZ2
        if (BZ_OK != BZ2_bzBuffToBuffCompress((char *)out+c_meta_len, out_size,
                                              (char *)in, in_size, 9, 0, 30))
            *out_size = in_size; // Didn't fit with bz2; force X_CAT below instead
#else
        fprintf(stderr, "Htscodecs has been compiled without libbz2 support\n");
        free(out);
        return NULL;
#endif

//      // lzma doesn't help generally, at least not for the name tokeniser
//      size_t lzma_size = 0;
//      lzma_easy_buffer_encode(9, LZMA_CHECK_CRC32, NULL,
//                              in, in_size, out+c_meta_len, &lzma_size,
//                              *out_size);
//      *out_size = lzma_size;

    } else {
        if (do_rle) {
            if (order == 0)
                arith_compress_O0_RLE(in, in_size, out+c_meta_len, out_size);
            else
                arith_compress_O1_RLE(in, in_size, out+c_meta_len, out_size);
        } else {
            //if (order == 2)
            //  arith_compress_O2(in, in_size, out+c_meta_len, out_size);
            //else
            if (order == 1)
                arith_compress_O1(in, in_size, out+c_meta_len, out_size);
            else
                arith_compress_O0(in, in_size, out+c_meta_len, out_size);
        }
    }

    if (*out_size >= in_size) {
        out[0] &= ~(3|X_EXT); // no entropy encoding, but keep e.g. PACK
        out[0] |= X_CAT | no_size;
        memcpy(out+c_meta_len, in, in_size);
        *out_size = in_size;
    }

    free(rle);
    free(packed);

    *out_size += c_meta_len;

    return out;
}

unsigned char *arith_compress(unsigned char *in, unsigned int in_size,
                              unsigned int *out_size, int order) {
    return arith_compress_to(in, in_size, NULL, out_size, order);
}

unsigned char *arith_uncompress_to(unsigned char *in,  unsigned int in_size,
                                   unsigned char *out, unsigned int *out_size) {
    unsigned char *in_end = in + in_size;
    unsigned char *out_free = NULL;
    unsigned char *tmp_free = NULL;

    if (in_size == 0)
        return NULL;

    if (*in & X_STRIPE) {
        unsigned int ulen, olen, c_meta_len = 1;
        int i;
        uint64_t clen_tot = 0;

        // Decode lengths
        c_meta_len += var_get_u32(in+c_meta_len, in_end, &ulen);
        if (c_meta_len >= in_size)
            return NULL;
        unsigned int N = in[c_meta_len++];
        if (N < 1)  // Must be at least one stripe
            return NULL;
        unsigned int clenN[256], ulenN[256], idxN[256];
        if (!out) {
            if (ulen >= INT_MAX)
                return NULL;
            if (!(out_free = out = malloc(ulen))) {
                return NULL;
            }
            *out_size = ulen;
        }
        if (ulen != *out_size) {
            free(out_free);
            return NULL;
        }

        for (i = 0; i < N; i++) {
            ulenN[i] = ulen / N + ((ulen % N) > i);
            idxN[i] = i ? idxN[i-1] + ulenN[i-1] : 0;
            c_meta_len += var_get_u32(in+c_meta_len, in_end, &clenN[i]);
            clen_tot += clenN[i];
            if (c_meta_len > in_size || clenN[i] > in_size || clenN[i] < 1) {
                free(out_free);
                return NULL;
            }
        }

        // We can call this with a larger buffer, but once we've determined
        // how much we really use we limit it so the recursion becomes easier
        // to limit.
        if (c_meta_len + clen_tot > in_size) {
            free(out_free);
            return NULL;
        }
        in_size = c_meta_len + clen_tot;

        //fprintf(stderr, "    stripe meta %d\n", c_meta_len); //c-size

        // Uncompress the N streams
        unsigned char *outN = malloc(ulen);
        if (!outN) {
            free(out_free);
            return NULL;
        }
        for (i = 0; i < N; i++) {
            olen = ulenN[i];
            if (in_size < c_meta_len) {
                free(out_free);
                free(outN);
                return NULL;
            }
            if (!arith_uncompress_to(in+c_meta_len, in_size-c_meta_len, outN + idxN[i], &olen)
                || olen != ulenN[i]) {
                free(out_free);
                free(outN);
                return NULL;
            }
            c_meta_len += clenN[i];
        }

        unstripe(out, outN, ulen, N, idxN);

        free(outN);
        *out_size = ulen;
        return out;
    }

    int order = *in++;  in_size--;
    int do_pack = order & X_PACK;
    int do_rle  = order & X_RLE;
    int do_cat  = order & X_CAT;
    int no_size = order & X_NOSZ;
    int do_ext  = order & X_EXT;
    order &= 3;

    int sz = 0;
    unsigned int osz;
    if (!no_size)
        sz = var_get_u32(in, in_end, &osz);
    else
        sz = 0, osz = *out_size;
    in += sz;
    in_size -= sz;

    if (osz >= INT_MAX)
        return NULL;

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
        // Limit maximum size to get fast turnaround on fuzzing test cases
        if (osz > 100000)
            goto err;
#endif

    if (no_size && !out)
        return NULL; // Need one or the other

    if (!out) {
        *out_size = osz;
        if (!(out_free = out = malloc(*out_size)))
            return NULL;
    } else {
        if (*out_size < osz)
            return NULL;
        *out_size = osz;
    }

    uint32_t c_meta_size = 0;
    unsigned int tmp1_size = *out_size;
    unsigned int tmp2_size = *out_size;
    unsigned char *tmp1 = NULL, *tmp2 = NULL, *tmp = NULL;

    // Need In, Out and Tmp buffers with temporary buffer of the same size
    // as output.  Our entropy decode is either arithmetic (with/without RLE)
    // or external (bz2, gzip, lzma) but with an optional unPACK transform
    // at the end.
    //
    // To avoid pointless memcpy when unpacking we switch around which
    // buffers we're writing to accordingly.

    // Format is pack meta data if present, followed by compressed data.
    if (do_pack) {
        if (!(tmp_free = tmp = malloc(*out_size)))
            goto err;
        tmp1 = tmp;  // uncompress
        tmp2 = out;  // unpack
    } else {
        // no pack
        tmp  = NULL;
        tmp1 = out;  // uncompress
        tmp2 = out;  // NOP
    }

    
    // Decode the bit-packing map.
    uint8_t map[16] = {0};
    int npacked_sym = 0;
    uint64_t unpacked_sz = 0; // FIXME: rename to packed_per_byte
    if (do_pack) {
        c_meta_size = hts_unpack_meta(in, in_size, *out_size, map, &npacked_sym);
        if (c_meta_size == 0)
            goto err;

        unpacked_sz = osz;
        in      += c_meta_size;
        in_size -= c_meta_size;

        // New unpacked size.  We could derive this bit from *out_size
        // and npacked_sym.
        unsigned int osz;
        sz = var_get_u32(in, in_end, &osz);
        in += sz;
        in_size -= sz;
        if (osz > tmp1_size)
            goto err;
        tmp1_size = osz;
    }

    //fprintf(stderr, "    meta_size %d bytes\n", (int)(in - orig_in)); //c-size

    // uncompress RLE data.  in -> tmp1
    if (in_size) {
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
        // Limit maximum size to get fast turnaround on fuzzing test cases
        if (tmp1_size > 100000)
            goto err;
#endif
        if (do_cat) {
            //fprintf(stderr, "    CAT %d\n", tmp1_size); //c-size
            if (tmp1_size > in_size)
                goto err;
            if (tmp1_size > *out_size)
                goto err;
            memcpy(tmp1, in, tmp1_size);
        } else if (do_ext) {
#ifdef HAVE_LIBBZ2
            if (BZ_OK != BZ2_bzBuffToBuffDecompress((char *)tmp1, &tmp1_size,
                                                    (char *)in, in_size, 0, 0))
                goto err;
#else
            fprintf(stderr, "Htscodecs has been compiled without libbz2 support\n");
            goto err;
#endif
          } else {
            // in -> tmp1
            if (do_rle) {
                tmp1 = order == 1
                    ? arith_uncompress_O1_RLE(in, in_size, tmp1, tmp1_size)
                    : arith_uncompress_O0_RLE(in, in_size, tmp1, tmp1_size);
            } else {
                //if (order == 2)
                //    tmp1 = arith_uncompress_O2(in, in_size, tmp1, tmp1_size)
                //else
                tmp1 = order == 1
                    ? arith_uncompress_O1(in, in_size, tmp1, tmp1_size)
                    : arith_uncompress_O0(in, in_size, tmp1, tmp1_size);
            }
            if (!tmp1)
                goto err;
        }
    } else {
        tmp1_size = 0;
    }

    if (do_pack) {
        // Unpack bits via pack-map.  tmp1 -> tmp2
        if (npacked_sym == 1)
            unpacked_sz = tmp1_size;
        //uint8_t *porig = unpack(tmp2, tmp2_size, unpacked_sz, npacked_sym, map);
        //memcpy(tmp3, porig, unpacked_sz);
        if (!hts_unpack(tmp1, tmp1_size, tmp2, unpacked_sz, npacked_sym, map))
            goto err;
        tmp2_size = unpacked_sz;
    } else {
        tmp2_size = tmp1_size;
    }

    if (tmp)
        free(tmp);

    *out_size = tmp2_size;
    return tmp2;

 err:
    free(tmp_free);
    free(out_free);
    return NULL;
}

unsigned char *arith_uncompress(unsigned char *in, unsigned int in_size,
                                unsigned int *out_size) {
    return arith_uncompress_to(in, in_size, NULL, out_size);
}
