/*
 * Copyright (c) 2017-2023 Genome Research Ltd.
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

// FIXME Can we get decoder to return the compressed sized read, avoiding
// us needing to store it?  Yes we can.  See c-size comments.  If we added all these
// together we could get rans_uncompress_to_4x16 to return the number of bytes
// consumed, avoiding the calling code from needed to explicitly stored the size.
// However the effect on name tokeniser is to save 0.1 to 0.2% so not worth it.

/*-------------------------------------------------------------------------- */
/*
 * Example wrapper to use the rans_byte.h functions included above.
 *
 * This demonstrates how to use, and unroll, an order-0 and order-1 frequency
 * model.
 */

#include "config.h"

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <limits.h>
#include <math.h>

#ifndef NO_THREADS
#include <pthread.h>
#endif

#include "rANS_word.h"
#include "rANS_static4x16.h"
#include "rANS_static16_int.h"
#include "pack.h"
#include "rle.h"
#include "utils.h"

#define TF_SHIFT 12
#define TOTFREQ (1<<TF_SHIFT)

// 9-11 is considerably faster in the O1 variant due to reduced table size.
// We auto-tune between 10 and 12 though.  Anywhere from 9 to 14 are viable.
#ifndef TF_SHIFT_O1
#define TF_SHIFT_O1 12
#endif
#ifndef TF_SHIFT_O1_FAST
#define TF_SHIFT_O1_FAST 10
#endif
#define TOTFREQ_O1 (1<<TF_SHIFT_O1)
#define TOTFREQ_O1_FAST (1<<TF_SHIFT_O1_FAST)


/*-----------------------------------------------------------------------------
 * Memory to memory compression functions.
 *
 * These are original versions without any manual loop unrolling. They
 * are easier to understand, but can be up to 2x slower.
 */

unsigned int rans_compress_bound_4x16(unsigned int size, int order) {
    int N = (order>>8) & 0xff;
    if (!N) N=4;

    order &= 0xff;
    unsigned int sz = (order == 0
        ? 1.05*size + 257*3 + 4
        : 1.05*size + 257*257*3 + 4 + 257*3+4) +
        ((order & RANS_ORDER_PACK) ? 1 : 0) +
        ((order & RANS_ORDER_RLE) ? 1 + 257*3+4: 0) + 20 +
        ((order & RANS_ORDER_X32) ? (32-4)*4 : 0) +
        ((order & RANS_ORDER_STRIPE) ? 7 + 5*N: 0);
    return sz + (sz&1) + 2; // make this even so buffers are word aligned
}

// Compresses in_size bytes from 'in' to *out_size bytes in 'out'.
//
// NB: The output buffer does not hold the original size, so it is up to
// the caller to store this.
unsigned char *rans_compress_O0_4x16(unsigned char *in, unsigned int in_size,
                                     unsigned char *out, unsigned int *out_size) {
    unsigned char *cp, *out_end;
    RansEncSymbol syms[256];
    RansState rans0;
    RansState rans2;
    RansState rans1;
    RansState rans3;
    uint8_t* ptr;
    uint32_t F[256+MAGIC] = {0};
    int i, j, tab_size = 0, rle, x;
    // -20 for order/size/meta
    uint32_t bound = rans_compress_bound_4x16(in_size,0)-20;

    if (!out) {
        *out_size = bound;
        out = malloc(*out_size);
    }
    if (!out || bound > *out_size)
        return NULL;

    // If "out" isn't word aligned, tweak out_end/ptr to ensure it is.
    // We already added more round in bound to allow for this.
    if (((size_t)out)&1)
        bound--;
    ptr = out_end = out + bound;

    if (in_size == 0)
        goto empty;

    // Compute statistics
    if (hist8(in, in_size, F) < 0)
        return NULL;

    // Normalise so frequences sum to power of 2
    uint32_t fsum = in_size;
    uint32_t max_val = round2(fsum);
    if (max_val > TOTFREQ)
        max_val = TOTFREQ;

    if (normalise_freq(F, fsum, max_val) < 0)
        return NULL;
    fsum=max_val;

    cp = out;
    cp += encode_freq(cp, F);
    tab_size = cp-out;
    //write(2, out+4, cp-(out+4));

    if (normalise_freq(F, fsum, TOTFREQ) < 0)
        return NULL;

    // Encode statistics.
    for (x = rle = j = 0; j < 256; j++) {
        if (F[j]) {
            RansEncSymbolInit(&syms[j], x, F[j], TF_SHIFT);
            x += F[j];
        }
    }

    RansEncInit(&rans0);
    RansEncInit(&rans1);
    RansEncInit(&rans2);
    RansEncInit(&rans3);

    switch (i=(in_size&3)) {
    case 3: RansEncPutSymbol(&rans2, &ptr, &syms[in[in_size-(i-2)]]);
        // fall-through
    case 2: RansEncPutSymbol(&rans1, &ptr, &syms[in[in_size-(i-1)]]);
        // fall-through
    case 1: RansEncPutSymbol(&rans0, &ptr, &syms[in[in_size-(i-0)]]);
        // fall-through
    case 0:
        break;
    }
    for (i=(in_size &~3); i>0; i-=4) {
        RansEncSymbol *s3 = &syms[in[i-1]];
        RansEncSymbol *s2 = &syms[in[i-2]];
        RansEncSymbol *s1 = &syms[in[i-3]];
        RansEncSymbol *s0 = &syms[in[i-4]];

#if 1
        RansEncPutSymbol(&rans3, &ptr, s3);
        RansEncPutSymbol(&rans2, &ptr, s2);
        RansEncPutSymbol(&rans1, &ptr, s1);
        RansEncPutSymbol(&rans0, &ptr, s0);
#else
        // Slightly beter on gcc, much better on clang
        uint16_t *ptr16 = (uint16_t *)ptr;

        if (rans3 >= s3->x_max) *--ptr16 = (uint16_t)rans3, rans3 >>= 16;
        if (rans2 >= s2->x_max) *--ptr16 = (uint16_t)rans2, rans2 >>= 16;
        uint32_t q3 = (uint32_t) (((uint64_t)rans3 * s3->rcp_freq) >> s3->rcp_shift);
        uint32_t q2 = (uint32_t) (((uint64_t)rans2 * s2->rcp_freq) >> s2->rcp_shift);
        rans3 += s3->bias + q3 * s3->cmpl_freq;
        rans2 += s2->bias + q2 * s2->cmpl_freq;

        if (rans1 >= s1->x_max) *--ptr16 = (uint16_t)rans1, rans1 >>= 16;
        if (rans0 >= s0->x_max) *--ptr16 = (uint16_t)rans0, rans0 >>= 16;
        uint32_t q1 = (uint32_t) (((uint64_t)rans1 * s1->rcp_freq) >> s1->rcp_shift);
        uint32_t q0 = (uint32_t) (((uint64_t)rans0 * s0->rcp_freq) >> s0->rcp_shift);
        rans1 += s1->bias + q1 * s1->cmpl_freq;
        rans0 += s0->bias + q0 * s0->cmpl_freq;

        ptr = (uint8_t *)ptr16;
#endif
    }

    RansEncFlush(&rans3, &ptr);
    RansEncFlush(&rans2, &ptr);
    RansEncFlush(&rans1, &ptr);
    RansEncFlush(&rans0, &ptr);

 empty:
    // Finalise block size and return it
    *out_size = (out_end - ptr) + tab_size;

    memmove(out + tab_size, ptr, out_end-ptr);

    return out;
}

unsigned char *rans_uncompress_O0_4x16(unsigned char *in, unsigned int in_size,
                                       unsigned char *out, unsigned int out_sz) {
    if (in_size < 16) // 4-states at least
        return NULL;

    if (out_sz >= INT_MAX)
        return NULL; // protect against some overflow cases

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (out_sz > 100000)
        return NULL;
#endif

    /* Load in the static tables */
    unsigned char *cp = in, *out_free = NULL;
    unsigned char *cp_end = in + in_size - 8; // within 8 => be extra safe
    int i, j;
    unsigned int x, y;
    uint16_t sfreq[TOTFREQ+32];
    uint16_t sbase[TOTFREQ+32]; // faster to use 32-bit on clang
    uint8_t  ssym [TOTFREQ+64]; // faster to use 16-bit on clang

    if (!out)
        out_free = out = malloc(out_sz);
    if (!out)
        return NULL;

    // Precompute reverse lookup of frequency.
    uint32_t F[256] = {0}, fsum;
    int fsz = decode_freq(cp, cp_end, F, &fsum);
    if (!fsz)
        goto err;
    cp += fsz;

    normalise_freq_shift(F, fsum, TOTFREQ);

    // Build symbols; fixme, do as part of decode, see the _d variant
    for (j = x = 0; j < 256; j++) {
        if (F[j]) {
            if (F[j] > TOTFREQ - x)
                goto err;
            for (y = 0; y < F[j]; y++) {
                ssym [y + x] = j;
                sfreq[y + x] = F[j];
                sbase[y + x] = y;
            }
            x += F[j];
        }
    }

    if (x != TOTFREQ)
        goto err;

    if (cp+16 > cp_end+8)
        goto err;

    RansState R[4];
    RansDecInit(&R[0], &cp); if (R[0] < RANS_BYTE_L) goto err;
    RansDecInit(&R[1], &cp); if (R[1] < RANS_BYTE_L) goto err;
    RansDecInit(&R[2], &cp); if (R[2] < RANS_BYTE_L) goto err;
    RansDecInit(&R[3], &cp); if (R[3] < RANS_BYTE_L) goto err;

// Simple version is comparable to below, but only with -O3
//
//    for (i = 0; cp < cp_end-8 && i < (out_sz&~7); i+=8) {
//        for(j=0; j<8;j++) {
//          RansState m = RansDecGet(&R[j%4], TF_SHIFT);
//          R[j%4] = sfreq[m] * (R[j%4] >> TF_SHIFT) + sbase[m];
//          out[i+j] = ssym[m];
//          RansDecRenorm(&R[j%4], &cp);
//        }
//    }

    for (i = 0; cp < cp_end-8 && i < (out_sz&~7); i+=8) {
        for (j = 0; j < 8; j+=4) {
            RansState m0 = RansDecGet(&R[0], TF_SHIFT);
            RansState m1 = RansDecGet(&R[1], TF_SHIFT);
            out[i+j+0] = ssym[m0];
            out[i+j+1] = ssym[m1];

            R[0] = sfreq[m0] * (R[0] >> TF_SHIFT) + sbase[m0];
            R[1] = sfreq[m1] * (R[1] >> TF_SHIFT) + sbase[m1];

            RansState m2 = RansDecGet(&R[2], TF_SHIFT);
            RansState m3 = RansDecGet(&R[3], TF_SHIFT);

            RansDecRenorm(&R[0], &cp);
            RansDecRenorm(&R[1], &cp);

            R[2] = sfreq[m2] * (R[2] >> TF_SHIFT) + sbase[m2];
            R[3] = sfreq[m3] * (R[3] >> TF_SHIFT) + sbase[m3];

            RansDecRenorm(&R[2], &cp);
            RansDecRenorm(&R[3], &cp);

            out[i+j+2] = ssym[m2];
            out[i+j+3] = ssym[m3];
        }
    }

    // remainder
    for (; i < out_sz; i++) {
        RansState m = RansDecGet(&R[i%4], TF_SHIFT);
        R[i%4] = sfreq[m] * (R[i%4] >> TF_SHIFT) + sbase[m];
        out[i] = ssym[m];
        RansDecRenormSafe(&R[i%4], &cp, cp_end+8);
    }

    //fprintf(stderr, "    0 Decoded %d bytes\n", (int)(cp-in)); //c-size

    return out;

 err:
    free(out_free);
    return NULL;
}

//-----------------------------------------------------------------------------

// Compute the entropy of 12-bit vs 10-bit frequency tables.
// 10 bit means smaller memory footprint when decoding and
// more speed due to cache hits, but it *may* be a poor
// compression fit.
int rans_compute_shift(uint32_t *F0, uint32_t (*F)[256], uint32_t *T,
                       uint32_t *S) {
    int i, j;

    double e10 = 0, e12 = 0;
    int max_tot = 0;
    for (i = 0; i < 256; i++) {
        if (F0[i] == 0)
            continue;
        unsigned int max_val = round2(T[i]);
        int ns = 0;
#define MAX(a,b) ((a)>(b)?(a):(b))

        // Number of samples that get their freq bumped to 1
        int sm10 = 0, sm12 = 0;
        for (j = 0; j < 256; j++) {
            if (F[i][j] && max_val / F[i][j] > TOTFREQ_O1_FAST)
                sm10++;
            if (F[i][j] && max_val / F[i][j] > TOTFREQ_O1)
                sm12++;
        }

        double l10 = log(TOTFREQ_O1_FAST + sm10);
        double l12 = log(TOTFREQ_O1      + sm12);
        double T_slow = (double)TOTFREQ_O1/T[i];
        double T_fast = (double)TOTFREQ_O1_FAST/T[i];

        for (j = 0; j < 256; j++) {
            if (F[i][j]) {
                ns++;

                e10 -= F[i][j] * (fast_log(MAX(F[i][j]*T_fast,1)) - l10);
                e12 -= F[i][j] * (fast_log(MAX(F[i][j]*T_slow,1)) - l12);

                // Estimation of compressed symbol freq table too.
                e10 += 1.3;
                e12 += 4.7;
            }
        }

        // Order-1 frequencies often end up totalling under TOTFREQ.
        // In this case it's smaller to output the real frequencies
        // prior to normalisation and normalise after (with an extra
        // normalisation step needed in the decoder too).
        //
        // Thus we normalise to a power of 2 only, store those,
        // and renormalise later here (and in decoder) by bit-shift
        // to get to the fixed size.
        if (ns < 64 && max_val > 128) max_val /= 2;
        if (max_val > 1024)           max_val /= 2;
        if (max_val > TOTFREQ_O1)     max_val = TOTFREQ_O1;
        S[i] = max_val; // scale to max this
        if (max_tot < max_val)
            max_tot = max_val;
    }
    int shift = e10/e12 < 1.01 || max_tot <= TOTFREQ_O1_FAST
        ? TF_SHIFT_O1_FAST
        : TF_SHIFT_O1;

//    fprintf(stderr, "e10/12 = %f %f %f, shift %d\n",
//          e10/log(256), e12/log(256), e10/e12, shift);

    return shift;
}

static
unsigned char *rans_compress_O1_4x16(unsigned char *in, unsigned int in_size,
                                     unsigned char *out, unsigned int *out_size) {
    unsigned char *cp, *out_end, *out_free = NULL;
    unsigned int tab_size;
    
    // -20 for order/size/meta
    uint32_t bound = rans_compress_bound_4x16(in_size,1)-20;

    if (!out) {
        *out_size = bound;
        out_free = out = malloc(*out_size);
    }
    if (!out || bound > *out_size)
        return NULL;

    if (((size_t)out)&1)
        bound--;
    out_end = out + bound;

    RansEncSymbol (*syms)[256] = htscodecs_tls_alloc(256 * (sizeof(*syms)));
    if (!syms) {
        free(out_free);
        return NULL;
    }

    cp = out;
    int shift = encode_freq1(in, in_size, 4, syms, &cp); 
    if (shift < 0) {
        htscodecs_tls_free(syms);
        return NULL;
    }
    tab_size = cp - out;

    RansState rans0, rans1, rans2, rans3;
    RansEncInit(&rans0);
    RansEncInit(&rans1);
    RansEncInit(&rans2);
    RansEncInit(&rans3);

    uint8_t* ptr = out_end;

    int isz4 = in_size>>2;
    int i0 = 1*isz4-2;
    int i1 = 2*isz4-2;
    int i2 = 3*isz4-2;
    int i3 = 4*isz4-2;

    unsigned char l0 = in[i0+1];
    unsigned char l1 = in[i1+1];
    unsigned char l2 = in[i2+1];
    unsigned char l3 = in[i3+1];

    // Deal with the remainder
    l3 = in[in_size-1];
    for (i3 = in_size-2; i3 > 4*isz4-2; i3--) {
        unsigned char c3 = in[i3];
        RansEncPutSymbol(&rans3, &ptr, &syms[c3][l3]);
        l3 = c3;
    }

    for (; i0 >= 0; i0--, i1--, i2--, i3--) {
        unsigned char c0, c1, c2, c3;
        RansEncSymbol *s3 = &syms[c3 = in[i3]][l3];
        RansEncSymbol *s2 = &syms[c2 = in[i2]][l2];
        RansEncSymbol *s1 = &syms[c1 = in[i1]][l1];
        RansEncSymbol *s0 = &syms[c0 = in[i0]][l0];

        RansEncPutSymbol(&rans3, &ptr, s3);
        RansEncPutSymbol(&rans2, &ptr, s2);
        RansEncPutSymbol(&rans1, &ptr, s1);
        RansEncPutSymbol(&rans0, &ptr, s0);

        l0 = c0;
        l1 = c1;
        l2 = c2;
        l3 = c3;
    }

    RansEncPutSymbol(&rans3, &ptr, &syms[0][l3]);
    RansEncPutSymbol(&rans2, &ptr, &syms[0][l2]);
    RansEncPutSymbol(&rans1, &ptr, &syms[0][l1]);
    RansEncPutSymbol(&rans0, &ptr, &syms[0][l0]);

    RansEncFlush(&rans3, &ptr);
    RansEncFlush(&rans2, &ptr);
    RansEncFlush(&rans1, &ptr);
    RansEncFlush(&rans0, &ptr);

    *out_size = (out_end - ptr) + tab_size;

    cp = out;
    memmove(out + tab_size, ptr, out_end-ptr);

    htscodecs_tls_free(syms);
    return out;
}

//#define MAGIC2 111
#define MAGIC2 179
//#define MAGIC2 0

static
unsigned char *rans_uncompress_O1_4x16(unsigned char *in, unsigned int in_size,
                                       unsigned char *out, unsigned int out_sz) {
    if (in_size < 16) // 4-states at least
        return NULL;

    if (out_sz >= INT_MAX)
        return NULL; // protect against some overflow cases

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (out_sz > 100000)
        return NULL;
#endif

    /* Load in the static tables */
    unsigned char *cp = in, *cp_end = in+in_size, *out_free = NULL;
    unsigned char *c_freq = NULL;
    int i, j = -999;
    unsigned int x;

    uint8_t *sfb_ = htscodecs_tls_alloc(256*(TOTFREQ_O1+MAGIC2)*sizeof(*sfb_));
    uint32_t (*s3)[TOTFREQ_O1_FAST] = (uint32_t (*)[TOTFREQ_O1_FAST])sfb_;
    // reuse the same memory for the fast mode lookup, but this only works
    // if we're on e.g. 12-bit freqs vs 10-bit freqs as needs 4x larger array.
    //uint32_t s3[256][TOTFREQ_O1_FAST];

    if (!sfb_)
        return NULL;
    fb_t (*fb)[256] = htscodecs_tls_alloc(256 * sizeof(*fb));
    if (!fb)
        goto err;
    uint8_t *sfb[256];
    if ((*cp >> 4) == TF_SHIFT_O1) {
        for (i = 0; i < 256; i++)
            sfb[i]=  sfb_ + i*(TOTFREQ_O1+MAGIC2);
    } else {
        for (i = 0; i < 256; i++)
            sfb[i]=  sfb_ + i*(TOTFREQ_O1_FAST+MAGIC2);
    }

    if (!out)
        out_free = out = malloc(out_sz);

    if (!out)
        goto err;

    //fprintf(stderr, "out_sz=%d\n", out_sz);

    // compressed header? If so uncompress it
    unsigned char *tab_end = NULL;
    unsigned char *c_freq_end = cp_end;
    unsigned int shift = *cp >> 4;
    if (*cp++ & 1) {
        uint32_t u_freq_sz, c_freq_sz;
        cp += var_get_u32(cp, cp_end, &u_freq_sz);
        cp += var_get_u32(cp, cp_end, &c_freq_sz);
        if (c_freq_sz > cp_end - cp)
            goto err;
        tab_end = cp + c_freq_sz;
        if (!(c_freq = rans_uncompress_O0_4x16(cp, c_freq_sz, NULL, u_freq_sz)))
            goto err;
        cp = c_freq;
        c_freq_end = c_freq + u_freq_sz;
    }

    // Decode order-0 symbol list; avoids needing in order-1 tables
    uint32_t F0[256] = {0};
    int fsz = decode_alphabet(cp, c_freq_end, F0);
    if (!fsz)
        goto err;
    cp += fsz;

    if (cp >= c_freq_end)
        goto err;

    const int s3_fast_on = in_size >= 100000;

    for (i = 0; i < 256; i++) {
        if (F0[i] == 0)
            continue;

        uint32_t F[256] = {0}, T = 0;
        fsz = decode_freq_d(cp, c_freq_end, F0, F, &T);
        if (!fsz)
            goto err;
        cp += fsz;

        if (!T) {
            //fprintf(stderr, "No freq for F_%d\n", i);
            continue;
        }

        normalise_freq_shift(F, T, 1<<shift);

        // Build symbols; fixme, do as part of decode, see the _d variant
        for (j = x = 0; j < 256; j++) {
            if (F[j]) {
                if (F[j] > (1<<shift) - x)
                    goto err;

                if (shift == TF_SHIFT_O1_FAST && s3_fast_on) {
                    int y;
                    for (y = 0; y < F[j]; y++)
                        s3[i][y+x] = (((uint32_t)F[j])<<(shift+8)) |(y<<8) |j;
                } else {
                    memset(&sfb[i][x], j, F[j]);
                    fb[i][j].f = F[j];
                    fb[i][j].b = x;
                }
                x += F[j];
            }
        }
        if (x != (1<<shift))
            goto err;
    }

    if (tab_end)
        cp = tab_end;
    free(c_freq);
    c_freq = NULL;

    if (cp+16 > cp_end)
        goto err;

    RansState rans0, rans1, rans2, rans3;
    uint8_t *ptr = cp, *ptr_end = in + in_size - 8;
    RansDecInit(&rans0, &ptr); if (rans0 < RANS_BYTE_L) goto err;
    RansDecInit(&rans1, &ptr); if (rans1 < RANS_BYTE_L) goto err;
    RansDecInit(&rans2, &ptr); if (rans2 < RANS_BYTE_L) goto err;
    RansDecInit(&rans3, &ptr); if (rans3 < RANS_BYTE_L) goto err;

    unsigned int isz4 = out_sz>>2;
    int l0 = 0, l1 = 0, l2 = 0, l3 = 0;
    unsigned int i4[] = {0*isz4, 1*isz4, 2*isz4, 3*isz4};

    RansState R[4];
    R[0] = rans0;
    R[1] = rans1;
    R[2] = rans2;
    R[3] = rans3;

    // Around 15% faster to specialise for 10/12 than to have one
    // loop with shift as a variable.
    if (shift == TF_SHIFT_O1) {
        // TF_SHIFT_O1 = 12

        const uint32_t mask = ((1u << TF_SHIFT_O1)-1);
        for (; i4[0] < isz4; i4[0]++, i4[1]++, i4[2]++, i4[3]++) {
            uint16_t m, c;
            c = sfb[l0][m = R[0] & mask];
            R[0] = fb[l0][c].f * (R[0]>>TF_SHIFT_O1) + m - fb[l0][c].b;
            out[i4[0]] = l0 = c;

            c = sfb[l1][m = R[1] & mask];
            R[1] = fb[l1][c].f * (R[1]>>TF_SHIFT_O1) + m - fb[l1][c].b;
            out[i4[1]] = l1 = c;

            c = sfb[l2][m = R[2] & mask];
            R[2] = fb[l2][c].f * (R[2]>>TF_SHIFT_O1) + m - fb[l2][c].b;
            out[i4[2]] = l2 = c;

            c = sfb[l3][m = R[3] & mask];
            R[3] = fb[l3][c].f * (R[3]>>TF_SHIFT_O1) + m - fb[l3][c].b;
            out[i4[3]] = l3 = c;

            if (ptr < ptr_end) {
                RansDecRenorm(&R[0], &ptr);
                RansDecRenorm(&R[1], &ptr);
                RansDecRenorm(&R[2], &ptr);
                RansDecRenorm(&R[3], &ptr);
            } else {
                RansDecRenormSafe(&R[0], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[1], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[2], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[3], &ptr, ptr_end+8);
            }
        }

        // Remainder
        for (; i4[3] < out_sz; i4[3]++) {
            uint32_t m3 = R[3] & ((1u<<TF_SHIFT_O1)-1);
            unsigned char c3 = sfb[l3][m3];
            out[i4[3]] = c3;
            R[3] = fb[l3][c3].f * (R[3]>>TF_SHIFT_O1) + m3 - fb[l3][c3].b;
            RansDecRenormSafe(&R[3], &ptr, ptr_end + 8);
            l3 = c3;
        }
    } else if (!s3_fast_on) {
        // TF_SHIFT_O1 = 10 with sfb[256][1024] & fb[256]256] array lookup
        // Slightly faster for -o193 on q4 (high comp), but also less
        // initialisation cost for smaller data
        const uint32_t mask = ((1u << TF_SHIFT_O1_FAST)-1);
        for (; i4[0] < isz4; i4[0]++, i4[1]++, i4[2]++, i4[3]++) {
            uint16_t m, c;
            c = sfb[l0][m = R[0] & mask];
            R[0] = fb[l0][c].f * (R[0]>>TF_SHIFT_O1_FAST) + m - fb[l0][c].b;
            out[i4[0]] = l0 = c;

            c = sfb[l1][m = R[1] & mask];
            R[1] = fb[l1][c].f * (R[1]>>TF_SHIFT_O1_FAST) + m - fb[l1][c].b;
            out[i4[1]] = l1 = c;

            c = sfb[l2][m = R[2] & mask];
            R[2] = fb[l2][c].f * (R[2]>>TF_SHIFT_O1_FAST) + m - fb[l2][c].b;
            out[i4[2]] = l2 = c;

            c = sfb[l3][m = R[3] & mask];
            R[3] = fb[l3][c].f * (R[3]>>TF_SHIFT_O1_FAST) + m - fb[l3][c].b;
            out[i4[3]] = l3 = c;

            if (ptr < ptr_end) {
                RansDecRenorm(&R[0], &ptr);
                RansDecRenorm(&R[1], &ptr);
                RansDecRenorm(&R[2], &ptr);
                RansDecRenorm(&R[3], &ptr);
            } else {
                RansDecRenormSafe(&R[0], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[1], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[2], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[3], &ptr, ptr_end+8);
            }
        }

        // Remainder
        for (; i4[3] < out_sz; i4[3]++) {
            uint32_t m3 = R[3] & ((1u<<TF_SHIFT_O1_FAST)-1);
            unsigned char c3 = sfb[l3][m3];
            out[i4[3]] = c3;
            R[3] = fb[l3][c3].f * (R[3]>>TF_SHIFT_O1_FAST) + m3 - fb[l3][c3].b;
            RansDecRenormSafe(&R[3], &ptr, ptr_end + 8);
            l3 = c3;
        }
    } else {
        // TF_SHIFT_O1_FAST.
        // Significantly faster for -o1 on q40 (low comp).
        // Higher initialisation cost, so only use if big blocks.
        const uint32_t mask = ((1u << TF_SHIFT_O1_FAST)-1);
        for (; i4[0] < isz4; i4[0]++, i4[1]++, i4[2]++, i4[3]++) {
            uint32_t S0 = s3[l0][R[0] & mask];
            uint32_t S1 = s3[l1][R[1] & mask];
            l0 = out[i4[0]] = S0;
            l1 = out[i4[1]] = S1;
            uint16_t F0 = S0>>(TF_SHIFT_O1_FAST+8);
            uint16_t F1 = S1>>(TF_SHIFT_O1_FAST+8);
            uint16_t B0 = (S0>>8) & mask;
            uint16_t B1 = (S1>>8) & mask;

            R[0] = F0 * (R[0]>>TF_SHIFT_O1_FAST) + B0;
            R[1] = F1 * (R[1]>>TF_SHIFT_O1_FAST) + B1;

            uint32_t S2 = s3[l2][R[2] & mask];
            uint32_t S3 = s3[l3][R[3] & mask];
            l2 = out[i4[2]] = S2;
            l3 = out[i4[3]] = S3;
            uint16_t F2 = S2>>(TF_SHIFT_O1_FAST+8);
            uint16_t F3 = S3>>(TF_SHIFT_O1_FAST+8);
            uint16_t B2 = (S2>>8) & mask;
            uint16_t B3 = (S3>>8) & mask;

            R[2] = F2 * (R[2]>>TF_SHIFT_O1_FAST) + B2;
            R[3] = F3 * (R[3]>>TF_SHIFT_O1_FAST) + B3;

            if (ptr < ptr_end) {
                RansDecRenorm(&R[0], &ptr);
                RansDecRenorm(&R[1], &ptr);
                RansDecRenorm(&R[2], &ptr);
                RansDecRenorm(&R[3], &ptr);
            } else {
                RansDecRenormSafe(&R[0], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[1], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[2], &ptr, ptr_end+8);
                RansDecRenormSafe(&R[3], &ptr, ptr_end+8);
            }
        }

        // Remainder
        for (; i4[3] < out_sz; i4[3]++) {
            uint32_t S = s3[l3][R[3] & ((1u<<TF_SHIFT_O1_FAST)-1)];
            l3 = out[i4[3]] = S;
            R[3] = (S>>(TF_SHIFT_O1_FAST+8)) * (R[3]>>TF_SHIFT_O1_FAST)
                + ((S>>8) & ((1u<<TF_SHIFT_O1_FAST)-1));
            RansDecRenormSafe(&R[3], &ptr, ptr_end + 8);
        }
    }
    //fprintf(stderr, "    1 Decoded %d bytes\n", (int)(ptr-in)); //c-size

    htscodecs_tls_free(fb);
    htscodecs_tls_free(sfb_);
    return out;

 err:
    htscodecs_tls_free(fb);
    htscodecs_tls_free(sfb_);
    free(out_free);
    free(c_freq);

    return NULL;
}

/*-----------------------------------------------------------------------------
 * r32x16 implementation, included here for now for simplicity
 */
#include "rANS_static32x16pr.h"

static int rans_cpu = 0xFFFF; // all

#if defined(__x86_64__) && \
    defined(HAVE_DECL___CPUID_COUNT)   && HAVE_DECL___CPUID_COUNT && \
    defined(HAVE_DECL___GET_CPUID_MAX) && HAVE_DECL___GET_CPUID_MAX
#include <cpuid.h>

#if defined(__clang__) && defined(__has_attribute)
#  if __has_attribute(unused)
#    define UNUSED __attribute__((unused))
#  else
#    define UNUSED
#  endif
#elif defined(__GNUC__) && __GNUC__ >= 3
#  define UNUSED __attribute__((unused))
#else
#  define UNUSED
#endif

// CPU detection is performed once.  NB this has an assumption that we're
// not migrating between processes with different instruction stes, but
// to date the only systems I know of that support this don't have different
// capabilities (that we use) per core.
#ifndef NO_THREADS
static pthread_once_t rans_cpu_once = PTHREAD_ONCE_INIT;
#endif

static int have_ssse3   UNUSED = 0;
static int have_sse4_1  UNUSED = 0;
static int have_popcnt  UNUSED = 0;
static int have_avx2    UNUSED = 0;
static int have_avx512f UNUSED = 0;
static int is_amd       UNUSED = 0;

#define HAVE_HTSCODECS_TLS_CPU_INIT
static void htscodecs_tls_cpu_init(void) {
    unsigned int eax = 0, ebx = 0, ecx = 0, edx = 0;
    // These may be unused, depending on HAVE_* config.h macros

    int level = __get_cpuid_max(0, NULL);
    __cpuid_count(0, 0, eax, ebx, ecx, edx);
    is_amd = (ecx == 0x444d4163);
    if (level >= 1) {
        __cpuid_count(1, 0, eax, ebx, ecx, edx);
#if defined(bit_SSSE3)
        have_ssse3 = ecx & bit_SSSE3;
#endif
#if defined(bit_POPCNT)
        have_popcnt = ecx & bit_POPCNT;
#endif
#if defined(bit_SSE4_1)
        have_sse4_1 = ecx & bit_SSE4_1;
#endif
    }
    if (level >= 7) {
        __cpuid_count(7, 0, eax, ebx, ecx, edx);
#if defined(bit_AVX2)
        have_avx2 = ebx & bit_AVX2;
#endif
#if defined(bit_AVX512F)
        have_avx512f = ebx & bit_AVX512F;
#endif
    }

    if (!have_popcnt) have_avx512f = have_avx2 = have_sse4_1 = 0;
    if (!have_ssse3)  have_sse4_1 = 0;
}

static inline
unsigned char *(*rans_enc_func(int do_simd, int order))
    (unsigned char *in,
     unsigned int in_size,
     unsigned char *out,
     unsigned int *out_size) {

    int have_e_sse4_1  = have_sse4_1;
    int have_e_avx2    = have_avx2;
    int have_e_avx512f = have_avx512f;

    if (!(rans_cpu & RANS_CPU_ENC_AVX512)) have_e_avx512f = 0;
    if (!(rans_cpu & RANS_CPU_ENC_AVX2))   have_e_avx2    = 0;
    if (!(rans_cpu & RANS_CPU_ENC_SSE4))   have_e_sse4_1  = 0;

    if (!do_simd) { // SIMD disabled
        return order & 1
            ? rans_compress_O1_4x16
            : rans_compress_O0_4x16;
    }

#ifdef NO_THREADS
    htscodecs_tls_cpu_init();
#else
    int err = pthread_once(&rans_cpu_once, htscodecs_tls_cpu_init);
    if (err != 0) {
        fprintf(stderr, "Initialising TLS data failed: pthread_once: %s\n",
                strerror(err));
        fprintf(stderr, "Using scalar code only\n");
    }
#endif

    if (order & 1) {
        // With simulated gathers, the AVX512 is now slower than AVX2, so
        // we avoid using it unless asking for the real avx512 gather.
        // Note for testing we do -c 0x0404 to enable AVX512 and disable AVX2.
        // We then need to call the avx512 func regardless.
        int use_gather;
#ifdef USE_GATHER
        use_gather = 1;
#else
        use_gather = !have_e_avx2;
#endif

#if defined(HAVE_AVX512)
        if (have_e_avx512f && (!is_amd || !have_e_avx2) && use_gather)
            return rans_compress_O1_32x16_avx512;
#endif
#if defined(HAVE_AVX2)
        if (have_e_avx2)
            return rans_compress_O1_32x16_avx2;
#endif
#if defined(HAVE_SSE4_1) && defined(HAVE_SSSE3) && defined(HAVE_POPCNT)
        if (have_e_sse4_1)
            return rans_compress_O1_32x16;
#endif
        return rans_compress_O1_32x16;
    } else {
#if defined(HAVE_AVX512)
        if (have_e_avx512f && (!is_amd || !have_e_avx2))
            return rans_compress_O0_32x16_avx512;
#endif
#if defined(HAVE_AVX2)
        if (have_e_avx2)
            return rans_compress_O0_32x16_avx2;
#endif
#if defined(HAVE_SSE4_1) && defined(HAVE_SSSE3) && defined(HAVE_POPCNT)
        if (have_e_sse4_1)
            return rans_compress_O0_32x16;
#endif
        return rans_compress_O0_32x16;
    }
}

static inline
unsigned char *(*rans_dec_func(int do_simd, int order))
    (unsigned char *in,
     unsigned int in_size,
     unsigned char *out,
     unsigned int out_size) {

    int have_d_sse4_1  = have_sse4_1;
    int have_d_avx2    = have_avx2;
    int have_d_avx512f = have_avx512f;

    if (!(rans_cpu & RANS_CPU_DEC_AVX512)) have_d_avx512f = 0;
    if (!(rans_cpu & RANS_CPU_DEC_AVX2))   have_d_avx2    = 0;
    if (!(rans_cpu & RANS_CPU_DEC_SSE4))   have_d_sse4_1  = 0;

    if (!do_simd) { // SIMD disabled
        return order & 1
            ? rans_uncompress_O1_4x16
            : rans_uncompress_O0_4x16;
    }

#ifdef NO_THREADS
    htscodecs_tls_cpu_init();
#else
    int err = pthread_once(&rans_cpu_once, htscodecs_tls_cpu_init);
    if (err != 0) {
        fprintf(stderr, "Initialising TLS data failed: pthread_once: %s\n",
                strerror(err));
        fprintf(stderr, "Using scalar code only\n");
    }
#endif

    if (order & 1) {
#if defined(HAVE_AVX512)
        if (have_d_avx512f)
            return rans_uncompress_O1_32x16_avx512;
#endif
#if defined(HAVE_AVX2)
        if (have_d_avx2)
            return rans_uncompress_O1_32x16_avx2;
#endif
#if defined(HAVE_SSE4_1) && defined(HAVE_SSSE3) && defined(HAVE_POPCNT)
        if (have_d_sse4_1)
            return rans_uncompress_O1_32x16_sse4;
#endif
        return rans_uncompress_O1_32x16;
    } else {
#if defined(HAVE_AVX512)
        if (have_d_avx512f)
            return rans_uncompress_O0_32x16_avx512;
#endif
#if defined(HAVE_AVX2)
        if (have_d_avx2)
            return rans_uncompress_O0_32x16_avx2;
#endif
#if defined(HAVE_SSE4_1) && defined(HAVE_SSSE3) && defined(HAVE_POPCNT)
        if (have_d_sse4_1)
            return rans_uncompress_O0_32x16_sse4;
#endif
        return rans_uncompress_O0_32x16;
    }
}

#elif defined(__ARM_NEON) && defined(__aarch64__)

#if defined(__linux__) || defined(__FreeBSD__)
#include <sys/auxv.h>
#elif defined(_WIN32)
#include <processthreadsapi.h>
#endif

static inline int have_neon(void) {
#if defined(__linux__) && defined(__arm__)
    return (getauxval(AT_HWCAP) & HWCAP_NEON) != 0;
#elif defined(__linux__) && defined(__aarch64__) && defined(HWCAP_ASIMD)
    return (getauxval(AT_HWCAP) & HWCAP_ASIMD) != 0;
#elif defined(__APPLE__)
    return 1;
#elif defined(__FreeBSD__) && defined(__arm__)
    unsigned long cap;
    if (elf_aux_info(AT_HWCAP, &cap, sizeof cap) != 0) return 0;
    return (cap & HWCAP_NEON) != 0;
#elif defined(__FreeBSD__) && defined(__aarch64__) && defined(HWCAP_ASIMD)
    unsigned long cap;
    if (elf_aux_info(AT_HWCAP, &cap, sizeof cap) != 0) return 0;
    return (cap & HWCAP_ASIMD) != 0;
#elif defined(_WIN32)
    return IsProcessorFeaturePresent(PF_ARM_V8_INSTRUCTIONS_AVAILABLE) != 0;
#else
    return 0;
#endif
}

static inline
unsigned char *(*rans_enc_func(int do_simd, int order))
    (unsigned char *in,
     unsigned int in_size,
     unsigned char *out,
     unsigned int *out_size) {

    if (do_simd) {
        if ((rans_cpu & RANS_CPU_ENC_NEON) && have_neon())
            return order & 1
                ? rans_compress_O1_32x16_neon
                : rans_compress_O0_32x16_neon;
        else
            return order & 1
                ? rans_compress_O1_32x16
                : rans_compress_O0_32x16;
    } else {
        return order & 1
            ? rans_compress_O1_4x16
            : rans_compress_O0_4x16;
    }
}

static inline
unsigned char *(*rans_dec_func(int do_simd, int order))
    (unsigned char *in,
     unsigned int in_size,
     unsigned char *out,
     unsigned int out_size) {

    if (do_simd) {
        if ((rans_cpu & RANS_CPU_DEC_NEON) && have_neon())
            return order & 1
                ? rans_uncompress_O1_32x16_neon
                : rans_uncompress_O0_32x16_neon;
        else
            return order & 1
                ? rans_uncompress_O1_32x16
                : rans_uncompress_O0_32x16;
    } else {
        return order & 1
            ? rans_uncompress_O1_4x16
            : rans_uncompress_O0_4x16;
    }
}

#else // !(defined(__GNUC__) && defined(__x86_64__)) && !defined(__ARM_NEON)

static inline
unsigned char *(*rans_enc_func(int do_simd, int order))
    (unsigned char *in,
     unsigned int in_size,
     unsigned char *out,
     unsigned int *out_size) {

    if (do_simd) {
        return order & 1
            ? rans_compress_O1_32x16
            : rans_compress_O0_32x16;
    } else {
        return order & 1
            ? rans_compress_O1_4x16
            : rans_compress_O0_4x16;
    }
}

static inline
unsigned char *(*rans_dec_func(int do_simd, int order))
    (unsigned char *in,
     unsigned int in_size,
     unsigned char *out,
     unsigned int out_size) {

    if (do_simd) {
        return order & 1
            ? rans_uncompress_O1_32x16
            : rans_uncompress_O0_32x16;
    } else {
        return order & 1
            ? rans_uncompress_O1_4x16
            : rans_uncompress_O0_4x16;
    }
}

#endif

// Test interface for restricting the auto-detection methods so we
// can forcibly compare different implementations on the same machine.
// See RANS_CPU_ defines in rANS_static4x16.h
void rans_set_cpu(int opts) {
    rans_cpu = opts;
#ifdef HAVE_HTSCODECS_TLS_CPU_INIT
    htscodecs_tls_cpu_init();
#endif
}

/*-----------------------------------------------------------------------------
 * Simple interface to the order-0 vs order-1 encoders and decoders.
 *
 * Smallest is method, <in_size> <input>, so worst case 2 bytes longer.
 */
unsigned char *rans_compress_to_4x16(unsigned char *in, unsigned int in_size,
                                     unsigned char *out,unsigned int *out_size,
                                     int order) {
    if (in_size > INT_MAX) {
        *out_size = 0;
        return NULL;
    }

    unsigned int c_meta_len;
    uint8_t *meta = NULL, *rle = NULL, *packed = NULL;
    uint8_t *out_free = NULL;

    if (!out) {
        *out_size = rans_compress_bound_4x16(in_size, order);
        if (*out_size == 0)
            return NULL;
        if (!(out_free = out = malloc(*out_size)))
            return NULL;
    }

    unsigned char *out_end = out + *out_size;

    // Permit 32-way unrolling for large blocks, paving the way for
    // AVX2 and AVX512 SIMD variants.
    if ((order & RANS_ORDER_SIMD_AUTO) && in_size >= 50000
        && !(order & RANS_ORDER_STRIPE))
        order |= X_32;

    if (in_size <= 20)
        order &= ~RANS_ORDER_STRIPE;
#ifndef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (in_size <= 1000)
        order &= ~RANS_ORDER_X32;
#endif
    if (order & RANS_ORDER_STRIPE) {
        int N = (order>>8) & 0xff;
        if (N == 0) N = 4; // default for compatibility with old tests

        unsigned char *transposed = malloc(in_size);
        unsigned int part_len[256];
        unsigned int idx[256];
        if (!transposed) {
            free(out_free);
            return NULL;
        }
        int i, j, x;

        for (i = 0; i < N; i++) {
            part_len[i] = in_size / N + ((in_size % N) > i);
            idx[i] = i ? idx[i-1] + part_len[i-1] : 0; // cumulative index
        }

#define KN 8
        i = x = 0;
        if (in_size >= N*KN) {
            for (; i < in_size-N*KN;) {
                int k;
                unsigned char *ink = in+i;
                for (j = 0; j < N; j++)
                    for (k = 0; k < KN; k++)
                        transposed[idx[j]+x+k] = ink[j+N*k];
                x += KN; i+=N*KN;
            }
        }
#undef KN
        for (; i < in_size-N; i += N, x++) {
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
        *out = order & ~RANS_ORDER_NOSZ;
        c_meta_len += var_put_u32(out+c_meta_len, out_end, in_size);
        out[c_meta_len++] = N;
        
        unsigned char *out_best = NULL;
        unsigned int out_best_len = 0;

        out2_start = out2 = out+7+5*N; // shares a buffer with c_meta
        for (i = 0; i < N; i++) {
            // Brute force try all methods.
            int j, m[] = {1,64,128,0}, best_j = 0, best_sz = in_size+10;
            for (j = 0; j < sizeof(m)/sizeof(*m); j++) {
                if ((order & m[j]) != m[j])
                    continue;

                // order-1 *only*; bit check above cannot elide order-0
                if ((order & RANS_ORDER_STRIPE_NO0) && (m[j]&1) == 0)
                    continue;
                olen2 = *out_size - (out2 - out);
                rans_compress_to_4x16(transposed+idx[i], part_len[i],
                                      out2, &olen2,
                                      m[j] | RANS_ORDER_NOSZ
                                      | (order&RANS_ORDER_X32));
                if (best_sz > olen2) {
                    best_sz = olen2;
                    best_j = j;
                    if (j < sizeof(m)/sizeof(*m) && olen2 > out_best_len) {
                        unsigned char *tmp = realloc(out_best, olen2);
                        if (!tmp) {
                            free(out_free);
                            return NULL;
                        }
                        out_best = tmp;
                        out_best_len = olen2;
                    }

                    // Cache a copy of the best so far
                    memcpy(out_best, out2, olen2);
                }
            }
            if (best_j < sizeof(m)/sizeof(*m)) {
                // Copy the best compression to output buffer if not current
                memcpy(out2, out_best, best_sz);
                olen2 = best_sz;
            }

            out2 += olen2;
            c_meta_len += var_put_u32(out+c_meta_len, out_end, olen2);
        }
        if (out_best)
            free(out_best);

        memmove(out+c_meta_len, out2_start, out2-out2_start);
        free(transposed);
        *out_size = c_meta_len + out2-out2_start;
        return out;
    }

    if (order & RANS_ORDER_CAT) {
        out[0] = RANS_ORDER_CAT;
        c_meta_len = 1;
        c_meta_len += var_put_u32(&out[1], out_end, in_size);
        if (in_size)
            memcpy(out+c_meta_len, in, in_size);
        *out_size = c_meta_len + in_size;
        return out;
    }

    int do_pack = order & RANS_ORDER_PACK;
    int do_rle  = order & RANS_ORDER_RLE;
    int no_size = order & RANS_ORDER_NOSZ;
    int do_simd = order & RANS_ORDER_X32;

    out[0] = order;
    c_meta_len = 1;

    if (!no_size)
        c_meta_len += var_put_u32(&out[1], out_end, in_size);

    order &= 3;

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
            out[0] &= ~RANS_ORDER_PACK;
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
        out[0] &= ~RANS_ORDER_PACK;
    }

    if (do_rle && in_size) {
        // RLE 'in' -> rle_length + rle_literals arrays
        unsigned int rmeta_len, c_rmeta_len;
        uint64_t rle_len;
        c_rmeta_len = in_size+257;
        if (!(meta = malloc(c_rmeta_len))) {
            free(out_free);
            return NULL;
        }

        uint8_t rle_syms[256];
        int rle_nsyms = 0;
        uint64_t rmeta_len64;
        rle = hts_rle_encode(in, in_size, meta, &rmeta_len64,
                             rle_syms, &rle_nsyms, NULL, &rle_len);
        memmove(meta+1+rle_nsyms, meta, rmeta_len64);
        meta[0] = rle_nsyms;
        memcpy(meta+1, rle_syms, rle_nsyms);
        rmeta_len = rmeta_len64 + rle_nsyms+1;

        if (!rle || rle_len + rmeta_len >= .99*in_size) {
            // Not worth the speed hit.
            out[0] &= ~RANS_ORDER_RLE;
            do_rle = 0;
            free(rle);
            rle = NULL;
        } else {
            // Compress lengths with O0 and literals with O0/O1 ("order" param)
            int sz = var_put_u32(out+c_meta_len, out_end, rmeta_len*2), sz2;
            sz += var_put_u32(out+c_meta_len+sz, out_end, rle_len);
            c_rmeta_len = *out_size - (c_meta_len+sz+5);
            rans_enc_func(do_simd, 0)(meta, rmeta_len, out+c_meta_len+sz+5, &c_rmeta_len);
            if (c_rmeta_len < rmeta_len) {
                sz2 = var_put_u32(out+c_meta_len+sz, out_end, c_rmeta_len);
                memmove(out+c_meta_len+sz+sz2, out+c_meta_len+sz+5, c_rmeta_len);
            } else {
                // Uncompressed RLE meta-data as too small
                sz = var_put_u32(out+c_meta_len, out_end, rmeta_len*2+1);
                sz2 = var_put_u32(out+c_meta_len+sz, out_end, rle_len);
                memcpy(out+c_meta_len+sz+sz2, meta, rmeta_len);
                c_rmeta_len = rmeta_len;
            }

            c_meta_len += sz + sz2 + c_rmeta_len;

            in = rle;
            in_size = rle_len;
        }

        free(meta);
    } else if (do_rle) {
        out[0] &= ~RANS_ORDER_RLE;
    }

    *out_size -= c_meta_len;
    if (order && in_size < 8) {
        out[0] &= ~1;
        order  &= ~1;
    }

    rans_enc_func(do_simd, order)(in, in_size, out+c_meta_len, out_size);

    if (*out_size >= in_size) {
        out[0] &= ~3;
        out[0] |= RANS_ORDER_CAT | no_size;
        if (in_size)
            memcpy(out+c_meta_len, in, in_size);
        *out_size = in_size;
    }

    free(rle);
    free(packed);

    *out_size += c_meta_len;

    return out;
}

unsigned char *rans_compress_4x16(unsigned char *in, unsigned int in_size,
                                  unsigned int *out_size, int order) {
    return rans_compress_to_4x16(in, in_size, NULL, out_size, order);
}

unsigned char *rans_uncompress_to_4x16(unsigned char *in,  unsigned int in_size,
                                       unsigned char *out, unsigned int *out_size) {
    unsigned char *in_end = in + in_size;
    unsigned char *out_free = NULL, *tmp_free = NULL, *meta_free = NULL;

    if (in_size == 0)
        return NULL;

    if (*in & RANS_ORDER_STRIPE) {
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
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
            if (ulen > 100000)
                return NULL;
#endif
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
            if (!rans_uncompress_to_4x16(in+c_meta_len, in_size-c_meta_len, outN + idxN[i], &olen)
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
    int do_pack = order & RANS_ORDER_PACK;
    int do_rle  = order & RANS_ORDER_RLE;
    int do_cat  = order & RANS_ORDER_CAT;
    int no_size = order & RANS_ORDER_NOSZ;
    int do_simd = order & RANS_ORDER_X32;
    order &= 1;

    int sz = 0;
    unsigned int osz;
    if (!no_size) {
        sz = var_get_u32(in, in_end, &osz);
    } else
        sz = 0, osz = *out_size;
    in += sz;
    in_size -= sz;

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (osz > 100000)
        return NULL;
#endif

    if (no_size && !out)
        goto err; // Need one or the other

    if (!out) {
        *out_size = osz;
        if (!(out = out_free = malloc(*out_size)))
            return NULL;
    } else {
        if (*out_size < osz)
        goto err;
        *out_size = osz;
    }

//    if (do_pack || do_rle) {
//      in += sz; // size field not needed when pure rANS
//      in_size -= sz;
//    }

    uint32_t c_meta_size = 0;
    unsigned int tmp1_size = *out_size;
    unsigned int tmp2_size = *out_size;
    unsigned int tmp3_size = *out_size;
    unsigned char *tmp1 = NULL, *tmp2 = NULL, *tmp3 = NULL, *tmp = NULL;

    // Need In, Out and Tmp buffers with temporary buffer of the same size
    // as output.  All use rANS, but with optional transforms (none, RLE,
    // Pack, or both).
    //
    //                    rans   unrle  unpack
    // If none:     in -> out
    // If RLE:      in -> tmp -> out
    // If Pack:     in -> tmp        -> out
    // If RLE+Pack: in -> out -> tmp -> out
    //                    tmp1   tmp2   tmp3
    //
    // So rans is in   -> tmp1
    // RLE     is tmp1 -> tmp2
    // Unpack  is tmp2 -> tmp3

    // Format is meta data (Pack and RLE in that order if present),
    // followed by rANS compressed data.

    if (do_pack || do_rle) {
        if (!(tmp = tmp_free = malloc(*out_size)))
            goto err;
        if (do_pack && do_rle) {
            tmp1 = out;
            tmp2 = tmp;
            tmp3 = out;
        } else if (do_pack) {
            tmp1 = tmp;
            tmp2 = tmp1;
            tmp3 = out;
        } else if (do_rle) {
            tmp1 = tmp;
            tmp2 = out;
            tmp3 = out;
        }
    } else {
        // neither
        tmp  = NULL;
        tmp1 = out;
        tmp2 = out;
        tmp3 = out;
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

    uint8_t *meta = NULL;
    uint32_t u_meta_size = 0;
    if (do_rle) {
        // Uncompress meta data
        uint32_t c_meta_size, rle_len, sz;
        sz  = var_get_u32(in,    in_end, &u_meta_size);
        sz += var_get_u32(in+sz, in_end, &rle_len);
        if (rle_len > tmp1_size) // should never grow
            goto err;
        if (u_meta_size & 1) {
            meta = in + sz;
            u_meta_size = u_meta_size/2 > (in_end-meta) ? (in_end-meta) : u_meta_size/2;
            c_meta_size = u_meta_size;
        } else {
            sz += var_get_u32(in+sz, in_end, &c_meta_size);
            u_meta_size /= 2;

            meta_free = meta = rans_dec_func(do_simd, 0)(in+sz, in_size-sz, NULL, u_meta_size);
            if (!meta)
                goto err;
        }
        if (c_meta_size+sz > in_size)
            goto err;
        in      += c_meta_size+sz;
        in_size -= c_meta_size+sz;
        tmp1_size = rle_len;
    }
    //fprintf(stderr, "    meta_size %d bytes\n", (int)(in - orig_in)); //c-size

    // uncompress RLE data.  in -> tmp1
    if (in_size) {
        if (do_cat) {
            //fprintf(stderr, "    CAT %d\n", tmp1_size); //c-size
            if (tmp1_size > in_size)
                goto err;
            if (tmp1_size > *out_size)
                goto err;
            memcpy(tmp1, in, tmp1_size);
        } else {
            tmp1 = rans_dec_func(do_simd, order)(in, in_size, tmp1, tmp1_size);
            if (!tmp1)
                goto err;
        }
    } else {
        tmp1_size = 0;
    }
    tmp2_size = tmp3_size = tmp1_size;

    if (do_rle) {
        // Unpack RLE.  tmp1 -> tmp2.
        if (u_meta_size == 0)
            goto err;
        uint64_t unrle_size = *out_size;
        int rle_nsyms = *meta ? *meta : 256;
        if (u_meta_size < 1+rle_nsyms)
            goto err;
        if (!hts_rle_decode(tmp1, tmp1_size,
                            meta+1+rle_nsyms, u_meta_size-(1+rle_nsyms),
                            meta+1, rle_nsyms, tmp2, &unrle_size))
            goto err;
        tmp3_size = tmp2_size = unrle_size;
        free(meta_free);
        meta_free = NULL;
    }
    if (do_pack) {
        // Unpack bits via pack-map.  tmp2 -> tmp3
        if (npacked_sym == 1)
            unpacked_sz = tmp2_size;
        //uint8_t *porig = unpack(tmp2, tmp2_size, unpacked_sz, npacked_sym, map);
        //memcpy(tmp3, porig, unpacked_sz);
        if (!hts_unpack(tmp2, tmp2_size, tmp3, unpacked_sz, npacked_sym, map))
            goto err;
        tmp3_size = unpacked_sz;
    }

    if (tmp)
        free(tmp);

    *out_size = tmp3_size;
    return tmp3;

 err:
    free(meta_free);
    free(out_free);
    free(tmp_free);
    return NULL;
}

unsigned char *rans_uncompress_4x16(unsigned char *in, unsigned int in_size,
                                    unsigned int *out_size) {
    return rans_uncompress_to_4x16(in, in_size, NULL, out_size);
}
