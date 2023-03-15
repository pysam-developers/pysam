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

#include "config.h"

#if defined(HAVE_SSE4_1) && defined(HAVE_SSSE3) && defined(HAVE_POPCNT)

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <x86intrin.h>

#include "rANS_word.h"
#include "rANS_static4x16.h"
#include "rANS_static16_int.h"
#include "varint.h"
#include "utils.h"

/* Uses: SSE, SSE2, SSSE3, SSE4.1 and POPCNT
SSE:
_mm_movemask_ps

SSE2:
    _mm_load_si128 _mm_store_si128
    _mm_set_epi32  _mm_set1_epi32
    _mm_and_si128  _mm_or_si128
    _mm_srli_epi32 _mm_slli_epi32 _mm_srli_epi64 _mm_set1_epi64x
    _mm_add_epi32
    _mm_packus_epi32
    _mm_andnot_si128
    _mm_cmpeq_epi32
    _mm_mul_epu32
    _mm_shuffle_epi32

SSSE3:
    _mm_shuffle_epi8

SSE4.1:
    _mm_mullo_epi32
    _mm_packus_epi32
    _mm_max_epu32
    _mm_cvtepu16_epi32
    _mm_blendv_epi8

POPCNT:
    _mm_popcnt_u32
 */


#define NX 32

#define LOAD128(a,b)                                    \
    __m128i a##1 = _mm_load_si128((__m128i *)&b[0]);    \
    __m128i a##2 = _mm_load_si128((__m128i *)&b[4]);    \
    __m128i a##3 = _mm_load_si128((__m128i *)&b[8]);    \
    __m128i a##4 = _mm_load_si128((__m128i *)&b[12]);   \
    __m128i a##5 = _mm_load_si128((__m128i *)&b[16]);   \
    __m128i a##6 = _mm_load_si128((__m128i *)&b[20]);   \
    __m128i a##7 = _mm_load_si128((__m128i *)&b[24]);   \
    __m128i a##8 = _mm_load_si128((__m128i *)&b[28]);

#define STORE128(a,b)                                   \
    _mm_store_si128((__m128i *)&b[ 0], a##1);           \
    _mm_store_si128((__m128i *)&b[ 4], a##2);           \
    _mm_store_si128((__m128i *)&b[ 8], a##3);           \
    _mm_store_si128((__m128i *)&b[12], a##4);           \
    _mm_store_si128((__m128i *)&b[16], a##5);           \
    _mm_store_si128((__m128i *)&b[20], a##6);           \
    _mm_store_si128((__m128i *)&b[24], a##7);           \
    _mm_store_si128((__m128i *)&b[28], a##8);

static inline __m128i _mm_i32gather_epi32x(int *b, __m128i idx, int size) {
    int c[4] __attribute__((aligned(32)));
    _mm_store_si128((__m128i *)c, idx);
    return _mm_set_epi32(b[c[3]], b[c[2]], b[c[1]], b[c[0]]);
}

// SSE4 implementation of the Order-0 encoder is poorly performing.
// Disabled for now.
#if 0
#define LOAD128v(a,b)                                   \
    __m128i a[8];                                       \
    a[0] = _mm_load_si128((__m128i *)&b[0]);            \
    a[1] = _mm_load_si128((__m128i *)&b[4]);            \
    a[2] = _mm_load_si128((__m128i *)&b[8]);            \
    a[3] = _mm_load_si128((__m128i *)&b[12]);           \
    a[4] = _mm_load_si128((__m128i *)&b[16]);           \
    a[5] = _mm_load_si128((__m128i *)&b[20]);           \
    a[6] = _mm_load_si128((__m128i *)&b[24]);           \
    a[7] = _mm_load_si128((__m128i *)&b[28]);

#define STORE128v(a,b)                                  \
    _mm_store_si128((__m128i *)&b[ 0], a[0]);           \
    _mm_store_si128((__m128i *)&b[ 4], a[1]);           \
    _mm_store_si128((__m128i *)&b[ 8], a[2]);           \
    _mm_store_si128((__m128i *)&b[12], a[3]);           \
    _mm_store_si128((__m128i *)&b[16], a[4]);           \
    _mm_store_si128((__m128i *)&b[20], a[5]);           \
    _mm_store_si128((__m128i *)&b[24], a[6]);           \
    _mm_store_si128((__m128i *)&b[28], a[7]);

static inline __m128i _mm_mulhi_epu32(__m128i a, __m128i b) {
    // Multiply bottom 4 items and top 4 items together.
    __m128i ab_hm = _mm_mul_epu32(_mm_srli_epi64(a, 32),_mm_srli_epi64(b, 32));
    __m128i ab_lm = _mm_srli_epi64(_mm_mul_epu32(a, b), 32);

    // Blend or and/or seems to make no difference.
    return _mm_blend_epi16(ab_lm, ab_hm, 0xcc);

//    // Shift to get hi 32-bit of each 64-bit product
//    ab_hm = _mm_and_si128(ab_hm,_mm_set1_epi64x((uint64_t)0xffffffff00000000));
//
//    return _mm_or_si128(ab_lm, ab_hm);
}

// Shift A>>B for non-constant B exists in AVX2, but not SSE world.
// We simulate this for now by store, shift, and load.  Ugly!
static inline __m128i _mm_srlv_epi32x(__m128i a, __m128i b) {
// Extract and inline shift.  Slowest clang, joint fastest gcc
//    return _mm_set_epi32(_mm_extract_epi32(a,3)>>_mm_extract_epi32(b,3),
//                       _mm_extract_epi32(a,2)>>_mm_extract_epi32(b,2),
//                       _mm_extract_epi32(a,1)>>_mm_extract_epi32(b,1),
//                       _mm_extract_epi32(a,0)>>_mm_extract_epi32(b,0));

// Half store and inline shift; Fastest gcc, comparable to others below clang
//    uint32_t A[4];
//    _mm_storeu_si128((__m128i *)&A, a);
//
//    return _mm_set_epi32(A[3]>>_mm_extract_epi32(b,3),
//                       A[2]>>_mm_extract_epi32(b,2),
//                       A[1]>>_mm_extract_epi32(b,1),
//                       A[0]>>_mm_extract_epi32(b,0));

// Other half
    uint32_t B[4];
    _mm_storeu_si128((__m128i *)&B, b);
    return _mm_set_epi32(_mm_extract_epi32(a,3)>>B[3],
                         _mm_extract_epi32(a,2)>>B[2],
                         _mm_extract_epi32(a,1)>>B[1],
                         _mm_extract_epi32(a,0)>>B[0]);

// Check if all b[] match, and constant shift if so.
// Too costly, even on q4 where it's common for all shift to be identical.
//    __m128i cmp = _mm_cmpeq_epi32(b, _mm_shuffle_epi32(b, 0x39));
//    if (_mm_movemask_ps((__m128)cmp) == 15) {
//      return _mm_srl_epi32(a, _mm_set1_epi64x(_mm_extract_epi32(b,0)));
//      //_mm_storeu_si128((__m128i *)&B,_mm_set1_epi32(_mm_extract_epi32(b,0)));
//    } else {
//      uint32_t B[4];
//      _mm_storeu_si128((__m128i *)&B, b);
//      return _mm_set_epi32(_mm_extract_epi32(a,3)>>B[3],
//                           _mm_extract_epi32(a,2)>>B[2],
//                           _mm_extract_epi32(a,1)>>B[1],
//                           _mm_extract_epi32(a,0)>>B[0]);
//    }


// Full store and inline shift
//   uint32_t A[4], B[4] __attribute__((aligned(16)));
//   _mm_storeu_si128((__m128i *)&A, a);
//   _mm_storeu_si128((__m128i *)&B, b);
//
//   return _mm_set_epi32(A[3]>>B[3], A[2]>>B[2], A[1]>>B[1], A[0]>>B[0]);

// Full store, shift and load
//    uint32_t A[4], B[4] __attribute__((aligned(16)));
//    _mm_storeu_si128((__m128i *)&A, a);
//    _mm_storeu_si128((__m128i *)&B, b);
//    A[0]>>=B[0];
//    A[1]>>=B[1];
//    A[2]>>=B[2];
//    A[3]>>=B[3];
//    return _mm_loadu_si128((__m128i *)A);
}

unsigned char *rans_compress_O0_32x16_sse4(unsigned char *in,
                                           unsigned int in_size,
                                           unsigned char *out,
                                           unsigned int *out_size) {
    unsigned char *cp, *out_end;
    RansEncSymbol syms[256];
    RansState ransN[NX];
    uint8_t* ptr;
    uint32_t F[256+MAGIC] = {0};
    int i, j, tab_size = 0, x, z;
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
    for (x = j = 0; j < 256; j++) {
        if (F[j]) {
            RansEncSymbolInit(&syms[j], x, F[j], TF_SHIFT);
            x += F[j];
        }
    }

    for (z = 0; z < NX; z++)
      RansEncInit(&ransN[z]);

    z = i = in_size&(NX-1);
    while (z-- > 0)
      RansEncPutSymbol(&ransN[z], &ptr, &syms[in[in_size-(i-z)]]);

    uint32_t SB[256], SA[256], SD[256], SC[256];

    // Build lookup tables for SIMD encoding
    uint16_t *ptr16 = (uint16_t *)ptr;
    for (i = 0; i < 256; i++) {
        SB[i] = syms[i].x_max;
        SA[i] = syms[i].rcp_freq;
        SD[i] = (syms[i].cmpl_freq<<0) | (syms[i].rcp_shift<<16);
        SC[i] = syms[i].bias;
    }

    LOAD128v(Rv, ransN);

    const __m128i shuf = _mm_set_epi8(0x80, 0x80, 0x80, 0x80,
                                      0x80, 0x80, 0x80, 0x80,
                                      0x0d, 0x0c, 0x09, 0x08,
                                      0x05, 0x04, 0x01, 0x00);

    // FIXME: slower!
    // q40:  340 (scalar) vs 300 (this)
    // q4:   430 (scalar) vs 304 (this)
    for (i=(in_size &~(NX-1)); i>0; i-=NX) {
      uint8_t *c = &in[i-NX];

      int h;
      // Slightly better in 4x8 instead of 2x16 cycles with clang,
      // but the reverse with gcc.
      //for (h=24; h >= 0; h -= 8) {
      for (h=16; h >= 0; h -= 16) {
        int H = h/4; // rans index
        uint8_t *C = &in[i-NX+h];

#define SET(i,a) _mm_set_epi32(a[C[i+3]],a[C[i+2]],a[C[i+1]],a[C[i+0]])
        __m128i xmax8 = SET(12, SB);
        __m128i xmax7 = SET( 8, SB);
        __m128i xmax6 = SET( 4, SB);
        __m128i xmax5 = SET( 0, SB);

        __m128i cv8 = _mm_cmpgt_epi32(Rv[H+3], xmax8);
        __m128i cv7 = _mm_cmpgt_epi32(Rv[H+2], xmax7);
        __m128i cv6 = _mm_cmpgt_epi32(Rv[H+1], xmax6);
        __m128i cv5 = _mm_cmpgt_epi32(Rv[H+0], xmax5);

        // Store bottom 16-bits at ptr16
        unsigned int imask8 = _mm_movemask_ps((__m128)cv8);
        unsigned int imask7 = _mm_movemask_ps((__m128)cv7);
        unsigned int imask6 = _mm_movemask_ps((__m128)cv6);
        unsigned int imask5 = _mm_movemask_ps((__m128)cv5);

#define X(A) 4*A,4*A+1,0x80,0x80
#define _ 0x80,0x80,0x80,0x80
        uint8_t permutec[16][16] __attribute__((aligned(16))) = {
            {  _ ,  _ ,  _ ,  _ },
            {  _ ,  _ ,  _ ,X(0)},
            {  _ ,  _ ,  _ ,X(1)},
            {  _ ,  _ ,X(0),X(1)},

            {  _ ,  _ ,  _ ,X(2)},
            {  _ ,  _ ,X(0),X(2)},
            {  _ ,  _ ,X(1),X(2)},
            {  _ ,X(0),X(1),X(2)},

            {  _ ,  _ ,  _ ,X(3)},
            {  _ ,  _ ,X(0),X(3)},
            {  _ ,  _ ,X(1),X(3)},
            {  _ ,X(0),X(1),X(3)},

            {  _ ,  _ ,X(2),X(3)},
            {  _ ,X(0),X(2),X(3)},
            {  _ ,X(1),X(2),X(3)},
            {X(0),X(1),X(2),X(3)},
        };
#undef X
#undef _

        __m128i idx8 = _mm_load_si128((__m128i *)permutec[imask8]);
        __m128i idx7 = _mm_load_si128((__m128i *)permutec[imask7]);
        __m128i idx6 = _mm_load_si128((__m128i *)permutec[imask6]);
        __m128i idx5 = _mm_load_si128((__m128i *)permutec[imask5]);

        // Permute; to gather together the rans states that need flushing
        __m128i V1, V2, V3, V4, V5, V6, V7, V8;
        V8 = _mm_shuffle_epi8(_mm_and_si128(Rv[H+3], cv8), idx8);
        V7 = _mm_shuffle_epi8(_mm_and_si128(Rv[H+2], cv7), idx7);
        V6 = _mm_shuffle_epi8(_mm_and_si128(Rv[H+1], cv6), idx6);
        V5 = _mm_shuffle_epi8(_mm_and_si128(Rv[H+0], cv5), idx5);

        // Shuffle alternating shorts together to collect low 16-bit
        // elements together.  ... 9 8 5 4 1 0.
        // Or as with avx2 code use packus instead.
        V8 = _mm_shuffle_epi8(V8, shuf);
        V7 = _mm_shuffle_epi8(V7, shuf);
        V6 = _mm_shuffle_epi8(V6, shuf);
        V5 = _mm_shuffle_epi8(V5, shuf);

        _mm_storeu_si64(ptr16-4, V8); ptr16 -= _mm_popcnt_u32(imask8);
        _mm_storeu_si64(ptr16-4, V7); ptr16 -= _mm_popcnt_u32(imask7);
        _mm_storeu_si64(ptr16-4, V6); ptr16 -= _mm_popcnt_u32(imask6);
        _mm_storeu_si64(ptr16-4, V5); ptr16 -= _mm_popcnt_u32(imask5);

        Rv[H+3] = _mm_blendv_epi8(Rv[H+3], _mm_srli_epi32(Rv[H+3], 16), cv8);
        Rv[H+2] = _mm_blendv_epi8(Rv[H+2], _mm_srli_epi32(Rv[H+2], 16), cv7);
        Rv[H+1] = _mm_blendv_epi8(Rv[H+1], _mm_srli_epi32(Rv[H+1], 16), cv6);
        Rv[H+0] = _mm_blendv_epi8(Rv[H+0], _mm_srli_epi32(Rv[H+0], 16), cv5);

        // Cannot trivially replace the multiply as mulhi_epu32 doesn't
        // exist (only mullo).
        // However we can use _mm_mul_epu32 twice to get 64bit results
        // (h our lanes) and shift/or to get the answer.
        //
        // (AVX512 allows us to hold it all in 64-bit lanes and use mullo_epi64
        // plus a shift.  KNC has mulhi_epi32, but not sure if this is
        // available.)
        __m128i rfv8 = SET(12, SA);
        __m128i rfv7 = SET( 8, SA);
        __m128i rfv6 = SET( 4, SA);
        __m128i rfv5 = SET( 0, SA);

        rfv8 = _mm_mulhi_epu32(Rv[H+3], rfv8);
        rfv7 = _mm_mulhi_epu32(Rv[H+2], rfv7);
        rfv6 = _mm_mulhi_epu32(Rv[H+1], rfv6);
        rfv5 = _mm_mulhi_epu32(Rv[H+0], rfv5);

        __m128i SDv8 = SET(12, SD);
        __m128i SDv7 = SET( 8, SD);
        __m128i SDv6 = SET( 4, SD);
        __m128i SDv5 = SET( 0, SD);

        __m128i shiftv8 = _mm_srli_epi32(SDv8, 16);
        __m128i shiftv7 = _mm_srli_epi32(SDv7, 16);
        __m128i shiftv6 = _mm_srli_epi32(SDv6, 16);
        __m128i shiftv5 = _mm_srli_epi32(SDv5, 16);

        __m128i freqv8 = _mm_and_si128(SDv8, _mm_set1_epi32(0xffff));
        __m128i freqv7 = _mm_and_si128(SDv7, _mm_set1_epi32(0xffff));
        __m128i freqv6 = _mm_and_si128(SDv6, _mm_set1_epi32(0xffff));
        __m128i freqv5 = _mm_and_si128(SDv5, _mm_set1_epi32(0xffff));

        // Bake this into the tabel to start with?
        shiftv8 = _mm_sub_epi32(shiftv8, _mm_set1_epi32(32));
        shiftv7 = _mm_sub_epi32(shiftv7, _mm_set1_epi32(32));
        shiftv6 = _mm_sub_epi32(shiftv6, _mm_set1_epi32(32));
        shiftv5 = _mm_sub_epi32(shiftv5, _mm_set1_epi32(32));

        // No way to shift by varying amounts.  Store, shift, load? Simulated
        __m128i qv8 = _mm_srlv_epi32x(rfv8, shiftv8);
        __m128i qv7 = _mm_srlv_epi32x(rfv7, shiftv7);
        __m128i qv6 = _mm_srlv_epi32x(rfv6, shiftv6);
        __m128i qv5 = _mm_srlv_epi32x(rfv5, shiftv5);

        qv8 = _mm_mullo_epi32(qv8, freqv8);
        qv7 = _mm_mullo_epi32(qv7, freqv7);
        qv6 = _mm_mullo_epi32(qv6, freqv6);
        qv5 = _mm_mullo_epi32(qv5, freqv5);

        qv8 = _mm_add_epi32(qv8, SET(12, SC));
        qv7 = _mm_add_epi32(qv7, SET( 8, SC));
        qv6 = _mm_add_epi32(qv6, SET( 4, SC));
        qv5 = _mm_add_epi32(qv5, SET( 0, SC));

        Rv[H+3] = _mm_add_epi32(Rv[H+3], qv8);
        Rv[H+2] = _mm_add_epi32(Rv[H+2], qv7);
        Rv[H+1] = _mm_add_epi32(Rv[H+1], qv6);
        Rv[H+0] = _mm_add_epi32(Rv[H+0], qv5);
      }
    }

    STORE128v(Rv, ransN);

    ptr = (uint8_t *)ptr16;

    for (z = NX-1; z >= 0; z--)
      RansEncFlush(&ransN[z], &ptr);

 empty:
    // Finalise block size and return it
    *out_size = (out_end - ptr) + tab_size;

//    cp = out;
//    *cp++ = (in_size>> 0) & 0xff;
//    *cp++ = (in_size>> 8) & 0xff;
//    *cp++ = (in_size>>16) & 0xff;
//    *cp++ = (in_size>>24) & 0xff;

    memmove(out + tab_size, ptr, out_end-ptr);

    return out;
}
#endif // disable SSE4 encoder

unsigned char *rans_uncompress_O0_32x16_sse4(unsigned char *in,
                                             unsigned int in_size,
                                             unsigned char *out,
                                             unsigned int out_sz) {
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
    unsigned char *cp_end = in + in_size;
    int i;
    uint32_t s3[TOTFREQ] __attribute__((aligned(32))); // For TF_SHIFT <= 12

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
    if (rans_F_to_s3(F, TF_SHIFT, s3))
        goto err;

    if (cp_end - cp < NX * 4)
        goto err;

    int z;
    RansState R[NX] __attribute__((aligned(32)));
    for (z = 0; z < NX; z++) {
        RansDecInit(&R[z], &cp);
        if (R[z] < RANS_BYTE_L)
            goto err;
    }

    uint16_t *sp = (uint16_t *)cp;

    int out_end = (out_sz&~(NX-1));
    const uint32_t mask = (1u << TF_SHIFT)-1;

    __m128i maskv  = _mm_set1_epi32(mask); // set mask in all lanes
    LOAD128(Rv, R);

    uint8_t overflow[72+64] = {0};
    for (i=0; i < out_end; i+=NX) {
        //for (z = 0; z < NX; z++)
        //  m[z] = R[z] & mask;
        __m128i masked1 = _mm_and_si128(Rv1, maskv);
        __m128i masked2 = _mm_and_si128(Rv2, maskv);
        __m128i masked3 = _mm_and_si128(Rv3, maskv);
        __m128i masked4 = _mm_and_si128(Rv4, maskv);

        //  S[z] = s3[m[z]];
        __m128i Sv1 = _mm_i32gather_epi32x((int *)s3, masked1, sizeof(*s3));
        __m128i Sv2 = _mm_i32gather_epi32x((int *)s3, masked2, sizeof(*s3));
        __m128i Sv3 = _mm_i32gather_epi32x((int *)s3, masked3, sizeof(*s3));
        __m128i Sv4 = _mm_i32gather_epi32x((int *)s3, masked4, sizeof(*s3));

        //  f[z] = S[z]>>(TF_SHIFT+8);
        __m128i fv1 = _mm_srli_epi32(Sv1, TF_SHIFT+8);
        __m128i fv2 = _mm_srli_epi32(Sv2, TF_SHIFT+8);
        __m128i fv3 = _mm_srli_epi32(Sv3, TF_SHIFT+8);
        __m128i fv4 = _mm_srli_epi32(Sv4, TF_SHIFT+8);

        //  b[z] = (S[z]>>8) & mask;
        __m128i bv1 = _mm_and_si128(_mm_srli_epi32(Sv1, 8), maskv);
        __m128i bv2 = _mm_and_si128(_mm_srli_epi32(Sv2, 8), maskv);
        __m128i bv3 = _mm_and_si128(_mm_srli_epi32(Sv3, 8), maskv);
        __m128i bv4 = _mm_and_si128(_mm_srli_epi32(Sv4, 8), maskv);

        //  s[z] = S[z] & 0xff;
        __m128i sv1 = _mm_and_si128(Sv1, _mm_set1_epi32(0xff));
        __m128i sv2 = _mm_and_si128(Sv2, _mm_set1_epi32(0xff));
        __m128i sv3 = _mm_and_si128(Sv3, _mm_set1_epi32(0xff));
        __m128i sv4 = _mm_and_si128(Sv4, _mm_set1_epi32(0xff));

        //  R[z] = f[z] * (R[z] >> TF_SHIFT) + b[z];
        Rv1 = _mm_add_epi32(
                  _mm_mullo_epi32(
                      _mm_srli_epi32(Rv1,TF_SHIFT), fv1), bv1);
        Rv2 = _mm_add_epi32(
                  _mm_mullo_epi32(
                      _mm_srli_epi32(Rv2,TF_SHIFT), fv2), bv2);
        Rv3 = _mm_add_epi32(
                  _mm_mullo_epi32(
                      _mm_srli_epi32(Rv3,TF_SHIFT), fv3), bv3);
        Rv4 = _mm_add_epi32(
                  _mm_mullo_epi32(
                      _mm_srli_epi32(Rv4,TF_SHIFT), fv4), bv4);

        // Tricky one:  out[i+z] = s[z];
        //             ---d---c ---b---a  sv1
        //             ---h---g ---f---e  sv2
        // packs_epi32 -h-g-f-e -d-c-b-a  sv1(2)
        // packs_epi16 ponmlkji hgfedcba  sv1(2) / sv3(4)
        sv1 = _mm_packus_epi32(sv1, sv2);
        sv3 = _mm_packus_epi32(sv3, sv4);
        sv1 = _mm_packus_epi16(sv1, sv3);

        // c =  R[z] < RANS_BYTE_L;
        // A little tricky as we only have signed comparisons.
        // See https://stackoverflow.com/questions/32945410/sse2-intrinsics-comparing-unsigned-integers

#define _mm_cmplt_epu32_imm(a,b) _mm_andnot_si128(_mm_cmpeq_epi32(_mm_max_epu32((a),_mm_set1_epi32(b)), (a)), _mm_set1_epi32(-1));

//#define _mm_cmplt_epu32_imm(a,b) _mm_cmpgt_epi32(_mm_set1_epi32((b)-0x80000000), _mm_xor_si128((a), _mm_set1_epi32(0x80000000)))

        __m128i renorm_mask1, renorm_mask2, renorm_mask3, renorm_mask4;
        renorm_mask1 = _mm_cmplt_epu32_imm(Rv1, RANS_BYTE_L);
        renorm_mask2 = _mm_cmplt_epu32_imm(Rv2, RANS_BYTE_L);
        renorm_mask3 = _mm_cmplt_epu32_imm(Rv3, RANS_BYTE_L);
        renorm_mask4 = _mm_cmplt_epu32_imm(Rv4, RANS_BYTE_L);

//#define P(A,B,C,D) ((A)+((B)<<2) + ((C)<<4) + ((D)<<6))
#define P(A,B,C,D)                              \
        { A+0,A+1,A+2,A+3,                      \
          B+0,B+1,B+2,B+3,                      \
          C+0,C+1,C+2,C+3,                      \
          D+0,D+1,D+2,D+3}
#ifdef _
#undef _
#endif
#define _ 0x80
        uint8_t pidx[16][16] = {
            P(_,_,_,_),
            P(0,_,_,_),
            P(_,0,_,_),
            P(0,4,_,_),

            P(_,_,0,_),
            P(0,_,4,_),
            P(_,0,4,_),
            P(0,4,8,_),

            P(_,_,_,0),
            P(0,_,_,4),
            P(_,0,_,4),
            P(0,4,_,8),

            P(_,_,0,4),
            P(0,_,4,8),
            P(_,0,4,8),
            P(0,4,8,12),
        };
#undef _

        // Protect against running off the end of in buffer.
        // We copy it to a worst-case local buffer when near the end.
        // 72 = 7*8(imask1..7) + 16;  worse case for 8th _mm_loadu_si128 call.
        // An extra 64 bytes is to avoid triggering this multiple times
        // after we swap sp/cp_end over.
        if ((uint8_t *)sp+72 > cp_end) {
            memmove(overflow, sp, cp_end - (uint8_t *)sp);
            sp = (uint16_t *)overflow;
            cp_end = (uint8_t *)overflow + sizeof(overflow);
        }

        // Shuffle the renorm values to correct lanes and incr sp pointer
        __m128i Vv1 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
        unsigned int imask1 = _mm_movemask_ps((__m128)renorm_mask1);
        Vv1 = _mm_shuffle_epi8(Vv1, _mm_load_si128((__m128i*)pidx[imask1]));
        sp += _mm_popcnt_u32(imask1);

        __m128i Vv2 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
        unsigned int imask2 = _mm_movemask_ps((__m128)renorm_mask2);
        sp += _mm_popcnt_u32(imask2);
        Vv2 = _mm_shuffle_epi8(Vv2, _mm_load_si128((__m128i*)pidx[imask2]));

        __m128i Vv3 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
        unsigned int imask3 = _mm_movemask_ps((__m128)renorm_mask3);
        Vv3 = _mm_shuffle_epi8(Vv3, _mm_load_si128((__m128i*)pidx[imask3]));
        sp += _mm_popcnt_u32(imask3);

        __m128i Vv4 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
        unsigned int imask4 = _mm_movemask_ps((__m128)renorm_mask4);
        sp += _mm_popcnt_u32(imask4);
        Vv4 = _mm_shuffle_epi8(Vv4, _mm_load_si128((__m128i*)pidx[imask4]));

        __m128i Yv1 = _mm_slli_epi32(Rv1, 16);
        __m128i Yv2 = _mm_slli_epi32(Rv2, 16);
        __m128i Yv3 = _mm_slli_epi32(Rv3, 16);
        __m128i Yv4 = _mm_slli_epi32(Rv4, 16);

        // y = (R[z] << 16) | V[z];
        Yv1 = _mm_or_si128(Yv1, Vv1);
        Yv2 = _mm_or_si128(Yv2, Vv2);
        Yv3 = _mm_or_si128(Yv3, Vv3);
        Yv4 = _mm_or_si128(Yv4, Vv4);

        // R[z] = c ? Y[z] : R[z];
        Rv1 = _mm_blendv_epi8(Rv1, Yv1, renorm_mask1);
        Rv2 = _mm_blendv_epi8(Rv2, Yv2, renorm_mask2);
        Rv3 = _mm_blendv_epi8(Rv3, Yv3, renorm_mask3);
        Rv4 = _mm_blendv_epi8(Rv4, Yv4, renorm_mask4);

        // ------------------------------------------------------------

        //  m[z] = R[z] & mask;
        __m128i masked5 = _mm_and_si128(Rv5, maskv);
        __m128i masked6 = _mm_and_si128(Rv6, maskv);
        __m128i masked7 = _mm_and_si128(Rv7, maskv);
        __m128i masked8 = _mm_and_si128(Rv8, maskv);

        //  S[z] = s3[m[z]];
        __m128i Sv5 = _mm_i32gather_epi32x((int *)s3, masked5, sizeof(*s3));
        __m128i Sv6 = _mm_i32gather_epi32x((int *)s3, masked6, sizeof(*s3));
        __m128i Sv7 = _mm_i32gather_epi32x((int *)s3, masked7, sizeof(*s3));
        __m128i Sv8 = _mm_i32gather_epi32x((int *)s3, masked8, sizeof(*s3));

        //  f[z] = S[z]>>(TF_SHIFT+8);
        __m128i fv5 = _mm_srli_epi32(Sv5, TF_SHIFT+8);
        __m128i fv6 = _mm_srli_epi32(Sv6, TF_SHIFT+8);
        __m128i fv7 = _mm_srli_epi32(Sv7, TF_SHIFT+8);
        __m128i fv8 = _mm_srli_epi32(Sv8, TF_SHIFT+8);

        //  b[z] = (S[z]>>8) & mask;
        __m128i bv5 = _mm_and_si128(_mm_srli_epi32(Sv5, 8), maskv);
        __m128i bv6 = _mm_and_si128(_mm_srli_epi32(Sv6, 8), maskv);
        __m128i bv7 = _mm_and_si128(_mm_srli_epi32(Sv7, 8), maskv);
        __m128i bv8 = _mm_and_si128(_mm_srli_epi32(Sv8, 8), maskv);

        //  s[z] = S[z] & 0xff;
        __m128i sv5 = _mm_and_si128(Sv5, _mm_set1_epi32(0xff));
        __m128i sv6 = _mm_and_si128(Sv6, _mm_set1_epi32(0xff));
        __m128i sv7 = _mm_and_si128(Sv7, _mm_set1_epi32(0xff));
        __m128i sv8 = _mm_and_si128(Sv8, _mm_set1_epi32(0xff));

        //  R[z] = f[z] * (R[z] >> TF_SHIFT) + b[z];
        Rv5 = _mm_add_epi32(
                  _mm_mullo_epi32(
                      _mm_srli_epi32(Rv5,TF_SHIFT), fv5), bv5);
        Rv6 = _mm_add_epi32(
                  _mm_mullo_epi32(
                      _mm_srli_epi32(Rv6,TF_SHIFT), fv6), bv6);
        Rv7 = _mm_add_epi32(
                  _mm_mullo_epi32(
                      _mm_srli_epi32(Rv7,TF_SHIFT), fv7), bv7);
        Rv8 = _mm_add_epi32(
                  _mm_mullo_epi32(
                      _mm_srli_epi32(Rv8,TF_SHIFT), fv8), bv8);

        // Tricky one:  out[i+z] = s[z];
        //             ---d---c ---b---a  sv1
        //             ---h---g ---f---e  sv2
        // packs_epi32 -h-g-f-e -d-c-b-a  sv1(2)
        // packs_epi16 ponmlkji hgfedcba  sv1(2) / sv3(4)
        sv5 = _mm_packus_epi32(sv5, sv6);
        sv7 = _mm_packus_epi32(sv7, sv8);
        sv5 = _mm_packus_epi16(sv5, sv7);

        // c =  R[z] < RANS_BYTE_L;
        __m128i renorm_mask5, renorm_mask6, renorm_mask7, renorm_mask8;
        renorm_mask5 = _mm_cmplt_epu32_imm(Rv5, RANS_BYTE_L);
        renorm_mask6 = _mm_cmplt_epu32_imm(Rv6, RANS_BYTE_L);
        renorm_mask7 = _mm_cmplt_epu32_imm(Rv7, RANS_BYTE_L);
        renorm_mask8 = _mm_cmplt_epu32_imm(Rv8, RANS_BYTE_L);
        
        // Shuffle the renorm values to correct lanes and incr sp pointer
        __m128i Vv5 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
        unsigned int imask5 = _mm_movemask_ps((__m128)renorm_mask5);
        Vv5 = _mm_shuffle_epi8(Vv5, _mm_load_si128((__m128i*)pidx[imask5]));
        sp += _mm_popcnt_u32(imask5);

        __m128i Vv6 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
        unsigned int imask6 = _mm_movemask_ps((__m128)renorm_mask6);
        sp += _mm_popcnt_u32(imask6);
        Vv6 = _mm_shuffle_epi8(Vv6, _mm_load_si128((__m128i*)pidx[imask6]));

        __m128i Vv7 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
        unsigned int imask7 = _mm_movemask_ps((__m128)renorm_mask7);
        Vv7 = _mm_shuffle_epi8(Vv7, _mm_load_si128((__m128i*)pidx[imask7]));
        sp += _mm_popcnt_u32(imask7);

        __m128i Vv8 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
        unsigned int imask8 = _mm_movemask_ps((__m128)renorm_mask8);
        sp += _mm_popcnt_u32(imask8);
        Vv8 = _mm_shuffle_epi8(Vv8, _mm_load_si128((__m128i*)pidx[imask8]));

        __m128i Yv5 = _mm_slli_epi32(Rv5, 16);
        __m128i Yv6 = _mm_slli_epi32(Rv6, 16);
        __m128i Yv7 = _mm_slli_epi32(Rv7, 16);
        __m128i Yv8 = _mm_slli_epi32(Rv8, 16);

        // y = (R[z] << 16) | V[z];
        Yv5 = _mm_or_si128(Yv5, Vv5);
        Yv6 = _mm_or_si128(Yv6, Vv6);
        Yv7 = _mm_or_si128(Yv7, Vv7);
        Yv8 = _mm_or_si128(Yv8, Vv8);

        // R[z] = c ? Y[z] : R[z];
        Rv5 = _mm_blendv_epi8(Rv5, Yv5, renorm_mask5);
        Rv6 = _mm_blendv_epi8(Rv6, Yv6, renorm_mask6);
        Rv7 = _mm_blendv_epi8(Rv7, Yv7, renorm_mask7);
        Rv8 = _mm_blendv_epi8(Rv8, Yv8, renorm_mask8);

        // Maybe just a store128 instead?
        _mm_storeu_si128((__m128i *)&out[i+ 0], sv1);
        _mm_storeu_si128((__m128i *)&out[i+16], sv5);
//      *(uint64_t *)&out[i+ 0] = _mm_extract_epi64(sv1, 0);
//      *(uint64_t *)&out[i+ 8] = _mm_extract_epi64(sv1, 1);
//      *(uint64_t *)&out[i+16] = _mm_extract_epi64(sv5, 0);
//      *(uint64_t *)&out[i+24] = _mm_extract_epi64(sv5, 1);
    }

    STORE128(Rv, R);

    for (z = out_sz & (NX-1); z-- > 0; )
      out[out_end + z] = s3[R[z] & mask];

    //fprintf(stderr, "    0 Decoded %d bytes\n", (int)(cp-in)); //c-size

    return out;

 err:
    free(out_free);
    return NULL;
}

//#define MAGIC2 111
#define MAGIC2 179
//#define MAGIC2 0

/*
 * A 32 x 32 matrix transpose and serialise from t[][] to out.
 * Storing in the other orientation speeds up the decoder, and we
 * can then flush to out in 1KB blocks.
 */
static inline void transpose_and_copy(uint8_t *out, int iN[32],
                                      uint8_t t[32][32]) {
    int z;
#ifdef UBSAN
    // Simplified version to avoid undefined behaviour sanitiser warnings.
    for (z = 0; z < NX; z++) {
        int k;
        for (k = 0; k < 32; k++)
            out[iN[z]+k] = t[k][z];
        iN[z] += 32;
    }
#else
    // Unaligned access.  We know we can get away with this as this
    // code is only ever executed on x86 platforms which permit this.
    for (z = 0; z < NX; z+=4) {
        *(uint64_t *)&out[iN[z]] =
            ((uint64_t)(t[0][z])<< 0) +
            ((uint64_t)(t[1][z])<< 8) +
            ((uint64_t)(t[2][z])<<16) +
            ((uint64_t)(t[3][z])<<24) +
            ((uint64_t)(t[4][z])<<32) +
            ((uint64_t)(t[5][z])<<40) +
            ((uint64_t)(t[6][z])<<48) +
            ((uint64_t)(t[7][z])<<56);
        *(uint64_t *)&out[iN[z+1]] =
            ((uint64_t)(t[0][z+1])<< 0) +
            ((uint64_t)(t[1][z+1])<< 8) +
            ((uint64_t)(t[2][z+1])<<16) +
            ((uint64_t)(t[3][z+1])<<24) +
            ((uint64_t)(t[4][z+1])<<32) +
            ((uint64_t)(t[5][z+1])<<40) +
            ((uint64_t)(t[6][z+1])<<48) +
            ((uint64_t)(t[7][z+1])<<56);
        *(uint64_t *)&out[iN[z+2]] =
            ((uint64_t)(t[0][z+2])<< 0) +
            ((uint64_t)(t[1][z+2])<< 8) +
            ((uint64_t)(t[2][z+2])<<16) +
            ((uint64_t)(t[3][z+2])<<24) +
            ((uint64_t)(t[4][z+2])<<32) +
            ((uint64_t)(t[5][z+2])<<40) +
            ((uint64_t)(t[6][z+2])<<48) +
            ((uint64_t)(t[7][z+2])<<56);
        *(uint64_t *)&out[iN[z+3]] =
            ((uint64_t)(t[0][z+3])<< 0) +
            ((uint64_t)(t[1][z+3])<< 8) +
            ((uint64_t)(t[2][z+3])<<16) +
            ((uint64_t)(t[3][z+3])<<24) +
            ((uint64_t)(t[4][z+3])<<32) +
            ((uint64_t)(t[5][z+3])<<40) +
            ((uint64_t)(t[6][z+3])<<48) +
            ((uint64_t)(t[7][z+3])<<56);

        *(uint64_t *)&out[iN[z]+8] =
            ((uint64_t)(t[8+0][z])<< 0) +
            ((uint64_t)(t[8+1][z])<< 8) +
            ((uint64_t)(t[8+2][z])<<16) +
            ((uint64_t)(t[8+3][z])<<24) +
            ((uint64_t)(t[8+4][z])<<32) +
            ((uint64_t)(t[8+5][z])<<40) +
            ((uint64_t)(t[8+6][z])<<48) +
            ((uint64_t)(t[8+7][z])<<56);
        *(uint64_t *)&out[iN[z+1]+8] =
            ((uint64_t)(t[8+0][z+1])<< 0) +
            ((uint64_t)(t[8+1][z+1])<< 8) +
            ((uint64_t)(t[8+2][z+1])<<16) +
            ((uint64_t)(t[8+3][z+1])<<24) +
            ((uint64_t)(t[8+4][z+1])<<32) +
            ((uint64_t)(t[8+5][z+1])<<40) +
            ((uint64_t)(t[8+6][z+1])<<48) +
            ((uint64_t)(t[8+7][z+1])<<56);
        *(uint64_t *)&out[iN[z+2]+8] =
            ((uint64_t)(t[8+0][z+2])<< 0) +
            ((uint64_t)(t[8+1][z+2])<< 8) +
            ((uint64_t)(t[8+2][z+2])<<16) +
            ((uint64_t)(t[8+3][z+2])<<24) +
            ((uint64_t)(t[8+4][z+2])<<32) +
            ((uint64_t)(t[8+5][z+2])<<40) +
            ((uint64_t)(t[8+6][z+2])<<48) +
            ((uint64_t)(t[8+7][z+2])<<56);
        *(uint64_t *)&out[iN[z+3]+8] =
            ((uint64_t)(t[8+0][z+3])<< 0) +
            ((uint64_t)(t[8+1][z+3])<< 8) +
            ((uint64_t)(t[8+2][z+3])<<16) +
            ((uint64_t)(t[8+3][z+3])<<24) +
            ((uint64_t)(t[8+4][z+3])<<32) +
            ((uint64_t)(t[8+5][z+3])<<40) +
            ((uint64_t)(t[8+6][z+3])<<48) +
            ((uint64_t)(t[8+7][z+3])<<56);

        *(uint64_t *)&out[iN[z]+16] =
            ((uint64_t)(t[16+0][z])<< 0) +
            ((uint64_t)(t[16+1][z])<< 8) +
            ((uint64_t)(t[16+2][z])<<16) +
            ((uint64_t)(t[16+3][z])<<24) +
            ((uint64_t)(t[16+4][z])<<32) +
            ((uint64_t)(t[16+5][z])<<40) +
            ((uint64_t)(t[16+6][z])<<48) +
            ((uint64_t)(t[16+7][z])<<56);
        *(uint64_t *)&out[iN[z+1]+16] =
            ((uint64_t)(t[16+0][z+1])<< 0) +
            ((uint64_t)(t[16+1][z+1])<< 8) +
            ((uint64_t)(t[16+2][z+1])<<16) +
            ((uint64_t)(t[16+3][z+1])<<24) +
            ((uint64_t)(t[16+4][z+1])<<32) +
            ((uint64_t)(t[16+5][z+1])<<40) +
            ((uint64_t)(t[16+6][z+1])<<48) +
            ((uint64_t)(t[16+7][z+1])<<56);
        *(uint64_t *)&out[iN[z+2]+16] =
            ((uint64_t)(t[16+0][z+2])<< 0) +
            ((uint64_t)(t[16+1][z+2])<< 8) +
            ((uint64_t)(t[16+2][z+2])<<16) +
            ((uint64_t)(t[16+3][z+2])<<24) +
            ((uint64_t)(t[16+4][z+2])<<32) +
            ((uint64_t)(t[16+5][z+2])<<40) +
            ((uint64_t)(t[16+6][z+2])<<48) +
            ((uint64_t)(t[16+7][z+2])<<56);
        *(uint64_t *)&out[iN[z+3]+16] =
            ((uint64_t)(t[16+0][z+3])<< 0) +
            ((uint64_t)(t[16+1][z+3])<< 8) +
            ((uint64_t)(t[16+2][z+3])<<16) +
            ((uint64_t)(t[16+3][z+3])<<24) +
            ((uint64_t)(t[16+4][z+3])<<32) +
            ((uint64_t)(t[16+5][z+3])<<40) +
            ((uint64_t)(t[16+6][z+3])<<48) +
            ((uint64_t)(t[16+7][z+3])<<56);

        *(uint64_t *)&out[iN[z]+24] =
            ((uint64_t)(t[24+0][z])<< 0) +
            ((uint64_t)(t[24+1][z])<< 8) +
            ((uint64_t)(t[24+2][z])<<16) +
            ((uint64_t)(t[24+3][z])<<24) +
            ((uint64_t)(t[24+4][z])<<32) +
            ((uint64_t)(t[24+5][z])<<40) +
            ((uint64_t)(t[24+6][z])<<48) +
            ((uint64_t)(t[24+7][z])<<56);
        *(uint64_t *)&out[iN[z+1]+24] =
            ((uint64_t)(t[24+0][z+1])<< 0) +
            ((uint64_t)(t[24+1][z+1])<< 8) +
            ((uint64_t)(t[24+2][z+1])<<16) +
            ((uint64_t)(t[24+3][z+1])<<24) +
            ((uint64_t)(t[24+4][z+1])<<32) +
            ((uint64_t)(t[24+5][z+1])<<40) +
            ((uint64_t)(t[24+6][z+1])<<48) +
            ((uint64_t)(t[24+7][z+1])<<56);
        *(uint64_t *)&out[iN[z+2]+24] =
            ((uint64_t)(t[24+0][z+2])<< 0) +
            ((uint64_t)(t[24+1][z+2])<< 8) +
            ((uint64_t)(t[24+2][z+2])<<16) +
            ((uint64_t)(t[24+3][z+2])<<24) +
            ((uint64_t)(t[24+4][z+2])<<32) +
            ((uint64_t)(t[24+5][z+2])<<40) +
            ((uint64_t)(t[24+6][z+2])<<48) +
            ((uint64_t)(t[24+7][z+2])<<56);
        *(uint64_t *)&out[iN[z+3]+24] =
            ((uint64_t)(t[24+0][z+3])<< 0) +
            ((uint64_t)(t[24+1][z+3])<< 8) +
            ((uint64_t)(t[24+2][z+3])<<16) +
            ((uint64_t)(t[24+3][z+3])<<24) +
            ((uint64_t)(t[24+4][z+3])<<32) +
            ((uint64_t)(t[24+5][z+3])<<40) +
            ((uint64_t)(t[24+6][z+3])<<48) +
            ((uint64_t)(t[24+7][z+3])<<56);

        iN[z+0] += 32;
        iN[z+1] += 32;
        iN[z+2] += 32;
        iN[z+3] += 32;
    }
#endif
}

unsigned char *rans_uncompress_O1_32x16_sse4(unsigned char *in,
                                             unsigned int in_size,
                                             unsigned char *out,
                                             unsigned int out_sz) {
    if (in_size < NX*4) // 4-states at least
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

    uint32_t (*s3)[TOTFREQ_O1] = htscodecs_tls_alloc(256*TOTFREQ_O1*4);
    if (!s3)
        return NULL;
    uint32_t (*s3F)[TOTFREQ_O1_FAST] = (uint32_t (*)[TOTFREQ_O1_FAST])s3;

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
        if (!(c_freq = rans_uncompress_O0_4x16(cp, c_freq_sz, NULL,u_freq_sz)))
            goto err;
        cp = c_freq;
        c_freq_end = c_freq + u_freq_sz;
    }

    // Decode order-0 symbol list; avoids needing in order-1 tables
    cp += decode_freq1(cp, c_freq_end, shift, s3, s3F, NULL, NULL);

    if (tab_end)
        cp = tab_end;
    free(c_freq);
    c_freq = NULL;

    if (cp_end - cp < NX * 4)
        goto err;

    RansState R[NX];
    uint8_t *ptr = cp, *ptr_end = in + in_size;
    int z;
    for (z = 0; z < NX; z++) {
        RansDecInit(&R[z], &ptr);
        if (R[z] < RANS_BYTE_L)
            goto err;
    }

    int isz4 = out_sz/NX;
    int i4[NX], l[NX] = {0};
    for (z = 0; z < NX; z++)
        i4[z] = z*isz4;

    // Around 15% faster to specialise for 10/12 than to have one
    // loop with shift as a variable.
    if (shift == TF_SHIFT_O1) {
        // TF_SHIFT_O1 = 12
        uint16_t *sp = (uint16_t *)ptr;
        const uint32_t mask = ((1u << TF_SHIFT_O1)-1);
        __m128i maskv  = _mm_set1_epi32(mask); // set mask in all lanes
        uint8_t tbuf[32][32];
        int tidx = 0;
        LOAD128(Rv, R);
        LOAD128(Lv, l);

        isz4 -= 64;
        for (; i4[0] < isz4 && (uint8_t *)sp+72 < ptr_end; ) {
            //for (z = 0; z < NX; z++)
            //  m[z] = R[z] & mask;
            __m128i masked1 = _mm_and_si128(Rv1, maskv);
            __m128i masked2 = _mm_and_si128(Rv2, maskv);
            __m128i masked3 = _mm_and_si128(Rv3, maskv);
            __m128i masked4 = _mm_and_si128(Rv4, maskv);

            Lv1 = _mm_slli_epi32(Lv1, TF_SHIFT_O1);
            Lv2 = _mm_slli_epi32(Lv2, TF_SHIFT_O1);
            Lv3 = _mm_slli_epi32(Lv3, TF_SHIFT_O1);
            Lv4 = _mm_slli_epi32(Lv4, TF_SHIFT_O1);
            masked1 = _mm_add_epi32(masked1, Lv1);
            masked2 = _mm_add_epi32(masked2, Lv2);
            masked3 = _mm_add_epi32(masked3, Lv3);
            masked4 = _mm_add_epi32(masked4, Lv4);

            //  S[z] = s3[l[z]][m[z]];
            __m128i Sv1 = _mm_i32gather_epi32x((int *)s3, masked1, sizeof(*s3));
            __m128i Sv2 = _mm_i32gather_epi32x((int *)s3, masked2, sizeof(*s3));
            __m128i Sv3 = _mm_i32gather_epi32x((int *)s3, masked3, sizeof(*s3));
            __m128i Sv4 = _mm_i32gather_epi32x((int *)s3, masked4, sizeof(*s3));

            //  f[z] = S[z]>>(TF_SHIFT+8);
            __m128i fv1 = _mm_srli_epi32(Sv1, TF_SHIFT_O1+8);
            __m128i fv2 = _mm_srli_epi32(Sv2, TF_SHIFT_O1+8);
            __m128i fv3 = _mm_srli_epi32(Sv3, TF_SHIFT_O1+8);
            __m128i fv4 = _mm_srli_epi32(Sv4, TF_SHIFT_O1+8);

            //  b[z] = (S[z]>>8) & mask;
            __m128i bv1 = _mm_and_si128(_mm_srli_epi32(Sv1, 8), maskv);
            __m128i bv2 = _mm_and_si128(_mm_srli_epi32(Sv2, 8), maskv);
            __m128i bv3 = _mm_and_si128(_mm_srli_epi32(Sv3, 8), maskv);
            __m128i bv4 = _mm_and_si128(_mm_srli_epi32(Sv4, 8), maskv);

            //  s[z] = S[z] & 0xff;
            __m128i sv1 = _mm_and_si128(Sv1, _mm_set1_epi32(0xff));
            __m128i sv2 = _mm_and_si128(Sv2, _mm_set1_epi32(0xff));
            __m128i sv3 = _mm_and_si128(Sv3, _mm_set1_epi32(0xff));
            __m128i sv4 = _mm_and_si128(Sv4, _mm_set1_epi32(0xff));

            // A maximum frequency of 4096 doesn't fit in our s3 array.
            // as it's 12 bit + 12 bit + 8 bit.  It wraps around to zero.
            // (We don't have this issue for TOTFREQ_O1_FAST.)
            //
            // Solution 1 is to change to spec to forbid freq of 4096.
            // Easy hack is to add an extra symbol so it sums correctly.
            //
            // Solution 2 implemented here is to look for the wrap around
            // and fix it.
            __m128i max_freq = _mm_set1_epi32(TOTFREQ_O1);
            __m128i zero = _mm_set1_epi32(0);
            __m128i cmp1 = _mm_cmpeq_epi32(fv1, zero);
            fv1 = _mm_blendv_epi8(fv1, max_freq, cmp1);
            __m128i cmp2 = _mm_cmpeq_epi32(fv2, zero);
            fv2 = _mm_blendv_epi8(fv2, max_freq, cmp2);
            __m128i cmp3 = _mm_cmpeq_epi32(fv3, zero);
            fv3 = _mm_blendv_epi8(fv3, max_freq, cmp3);
            __m128i cmp4 = _mm_cmpeq_epi32(fv4, zero);
            fv4 = _mm_blendv_epi8(fv4, max_freq, cmp4);

            //  R[z] = f[z] * (R[z] >> TF_SHIFT_O1) + b[z];
            Rv1 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv1,TF_SHIFT_O1), fv1), bv1);
            Rv2 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv2,TF_SHIFT_O1), fv2), bv2);
            Rv3 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv3,TF_SHIFT_O1), fv3), bv3);
            Rv4 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv4,TF_SHIFT_O1), fv4), bv4);

            Lv1 = sv1;
            Lv2 = sv2;
            Lv3 = sv3;
            Lv4 = sv4;

            // Tricky one:  out[i+z] = s[z];
            //             ---d---c ---b---a  sv1
            //             ---h---g ---f---e  sv2
            // packs_epi32 -h-g-f-e -d-c-b-a  sv1(2)
            // packs_epi16 ponmlkji hgfedcba  sv1(2) / sv3(4)
            sv1 = _mm_packus_epi32(sv1, sv2);
            sv3 = _mm_packus_epi32(sv3, sv4);
            sv1 = _mm_packus_epi16(sv1, sv3);

            // c =  R[z] < RANS_BYTE_L;
            // A little tricky as we only have signed comparisons.
            // See https://stackoverflow.com/questions/32945410/sse2-intrinsics-comparing-unsigned-integers

//#define _mm_cmplt_epu32_imm(a,b) _mm_andnot_si128(_mm_cmpeq_epi32(_mm_max_epu32((a),_mm_set1_epi32(b)), (a)), _mm_set1_epi32(-1));

            //#define _mm_cmplt_epu32_imm(a,b) _mm_cmpgt_epi32(_mm_set1_epi32((b)-0x80000000), _mm_xor_si128((a), _mm_set1_epi32(0x80000000)))

            __m128i renorm_mask1, renorm_mask2, renorm_mask3, renorm_mask4;
            renorm_mask1 = _mm_cmplt_epu32_imm(Rv1, RANS_BYTE_L);
            renorm_mask2 = _mm_cmplt_epu32_imm(Rv2, RANS_BYTE_L);
            renorm_mask3 = _mm_cmplt_epu32_imm(Rv3, RANS_BYTE_L);
            renorm_mask4 = _mm_cmplt_epu32_imm(Rv4, RANS_BYTE_L);

            //#define P(A,B,C,D) ((A)+((B)<<2) + ((C)<<4) + ((D)<<6))
#define P(A,B,C,D)                              \
            { A+0,A+1,A+2,A+3,                  \
              B+0,B+1,B+2,B+3,                  \
              C+0,C+1,C+2,C+3,                  \
              D+0,D+1,D+2,D+3}
#ifdef _
#undef _
#endif
#define _ 0x80
            uint8_t pidx[16][16] = {
                P(_,_,_,_),
                P(0,_,_,_),
                P(_,0,_,_),
                P(0,4,_,_),

                P(_,_,0,_),
                P(0,_,4,_),
                P(_,0,4,_),
                P(0,4,8,_),

                P(_,_,_,0),
                P(0,_,_,4),
                P(_,0,_,4),
                P(0,4,_,8),

                P(_,_,0,4),
                P(0,_,4,8),
                P(_,0,4,8),
                P(0,4,8,12),
            };
#undef _

            // Shuffle the renorm values to correct lanes and incr sp pointer
            __m128i Vv1 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask1 = _mm_movemask_ps((__m128)renorm_mask1);
            Vv1 = _mm_shuffle_epi8(Vv1, _mm_load_si128((__m128i*)pidx[imask1]));
            sp += _mm_popcnt_u32(imask1);

            __m128i Vv2 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask2 = _mm_movemask_ps((__m128)renorm_mask2);
            sp += _mm_popcnt_u32(imask2);
            Vv2 = _mm_shuffle_epi8(Vv2, _mm_load_si128((__m128i*)pidx[imask2]));

            __m128i Vv3 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask3 = _mm_movemask_ps((__m128)renorm_mask3);
            Vv3 = _mm_shuffle_epi8(Vv3, _mm_load_si128((__m128i*)pidx[imask3]));
            sp += _mm_popcnt_u32(imask3);

            __m128i Vv4 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask4 = _mm_movemask_ps((__m128)renorm_mask4);
            sp += _mm_popcnt_u32(imask4);
            Vv4 = _mm_shuffle_epi8(Vv4, _mm_load_si128((__m128i*)pidx[imask4]));

            __m128i Yv1 = _mm_slli_epi32(Rv1, 16);
            __m128i Yv2 = _mm_slli_epi32(Rv2, 16);
            __m128i Yv3 = _mm_slli_epi32(Rv3, 16);
            __m128i Yv4 = _mm_slli_epi32(Rv4, 16);

            // y = (R[z] << 16) | V[z];
            Yv1 = _mm_or_si128(Yv1, Vv1);
            Yv2 = _mm_or_si128(Yv2, Vv2);
            Yv3 = _mm_or_si128(Yv3, Vv3);
            Yv4 = _mm_or_si128(Yv4, Vv4);

            // R[z] = c ? Y[z] : R[z];
            Rv1 = _mm_blendv_epi8(Rv1, Yv1, renorm_mask1);
            Rv2 = _mm_blendv_epi8(Rv2, Yv2, renorm_mask2);
            Rv3 = _mm_blendv_epi8(Rv3, Yv3, renorm_mask3);
            Rv4 = _mm_blendv_epi8(Rv4, Yv4, renorm_mask4);

            // ------------------------------------------------------------

            //  m[z] = R[z] & mask;
            __m128i masked5 = _mm_and_si128(Rv5, maskv);
            __m128i masked6 = _mm_and_si128(Rv6, maskv);
            __m128i masked7 = _mm_and_si128(Rv7, maskv);
            __m128i masked8 = _mm_and_si128(Rv8, maskv);


            Lv5 = _mm_slli_epi32(Lv5, TF_SHIFT_O1);
            Lv6 = _mm_slli_epi32(Lv6, TF_SHIFT_O1);
            Lv7 = _mm_slli_epi32(Lv7, TF_SHIFT_O1);
            Lv8 = _mm_slli_epi32(Lv8, TF_SHIFT_O1);
            masked5 = _mm_add_epi32(masked5, Lv5);
            masked6 = _mm_add_epi32(masked6, Lv6);
            masked7 = _mm_add_epi32(masked7, Lv7);
            masked8 = _mm_add_epi32(masked8, Lv8);

            //  S[z] = s3[m[z]];
            __m128i Sv5 = _mm_i32gather_epi32x((int *)s3, masked5, sizeof(*s3));
            __m128i Sv6 = _mm_i32gather_epi32x((int *)s3, masked6, sizeof(*s3));
            __m128i Sv7 = _mm_i32gather_epi32x((int *)s3, masked7, sizeof(*s3));
            __m128i Sv8 = _mm_i32gather_epi32x((int *)s3, masked8, sizeof(*s3));

            //  f[z] = S[z]>>(TF_SHIFT_O1+8);
            __m128i fv5 = _mm_srli_epi32(Sv5, TF_SHIFT_O1+8);
            __m128i fv6 = _mm_srli_epi32(Sv6, TF_SHIFT_O1+8);
            __m128i fv7 = _mm_srli_epi32(Sv7, TF_SHIFT_O1+8);
            __m128i fv8 = _mm_srli_epi32(Sv8, TF_SHIFT_O1+8);

            //  b[z] = (S[z]>>8) & mask;
            __m128i bv5 = _mm_and_si128(_mm_srli_epi32(Sv5, 8), maskv);
            __m128i bv6 = _mm_and_si128(_mm_srli_epi32(Sv6, 8), maskv);
            __m128i bv7 = _mm_and_si128(_mm_srli_epi32(Sv7, 8), maskv);
            __m128i bv8 = _mm_and_si128(_mm_srli_epi32(Sv8, 8), maskv);

            //  s[z] = S[z] & 0xff;
            __m128i sv5 = _mm_and_si128(Sv5, _mm_set1_epi32(0xff));
            __m128i sv6 = _mm_and_si128(Sv6, _mm_set1_epi32(0xff));
            __m128i sv7 = _mm_and_si128(Sv7, _mm_set1_epi32(0xff));
            __m128i sv8 = _mm_and_si128(Sv8, _mm_set1_epi32(0xff));

            // A maximum frequency of 4096 doesn't fit in our s3 array. Fix
            __m128i cmp5 = _mm_cmpeq_epi32(fv5, zero);
            fv5 = _mm_blendv_epi8(fv5, max_freq, cmp5);
            __m128i cmp6 = _mm_cmpeq_epi32(fv6, zero);
            fv6 = _mm_blendv_epi8(fv6, max_freq, cmp6);
            __m128i cmp7 = _mm_cmpeq_epi32(fv7, zero);
            fv7 = _mm_blendv_epi8(fv7, max_freq, cmp7);
            __m128i cmp8 = _mm_cmpeq_epi32(fv8, zero);
            fv8 = _mm_blendv_epi8(fv8, max_freq, cmp8);

            //  R[z] = f[z] * (R[z] >> TF_SHIFT_O1) + b[z];
            Rv5 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv5,TF_SHIFT_O1), fv5), bv5);
            Rv6 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv6,TF_SHIFT_O1), fv6), bv6);
            Rv7 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv7,TF_SHIFT_O1), fv7), bv7);
            Rv8 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv8,TF_SHIFT_O1), fv8), bv8);

            Lv5 = sv5;
            Lv6 = sv6;
            Lv7 = sv7;
            Lv8 = sv8;

            // Tricky one:  out[i+z] = s[z];
            //             ---d---c ---b---a  sv1
            //             ---h---g ---f---e  sv2
            // packs_epi32 -h-g-f-e -d-c-b-a  sv1(2)
            // packs_epi16 ponmlkji hgfedcba  sv1(2) / sv3(4)
            sv5 = _mm_packus_epi32(sv5, sv6);
            sv7 = _mm_packus_epi32(sv7, sv8);
            sv5 = _mm_packus_epi16(sv5, sv7);

            // c =  R[z] < RANS_BYTE_L;
            __m128i renorm_mask5, renorm_mask6, renorm_mask7, renorm_mask8;
            renorm_mask5 = _mm_cmplt_epu32_imm(Rv5, RANS_BYTE_L);
            renorm_mask6 = _mm_cmplt_epu32_imm(Rv6, RANS_BYTE_L);
            renorm_mask7 = _mm_cmplt_epu32_imm(Rv7, RANS_BYTE_L);
            renorm_mask8 = _mm_cmplt_epu32_imm(Rv8, RANS_BYTE_L);
        
            // Shuffle the renorm values to correct lanes and incr sp pointer
            __m128i Vv5 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask5 = _mm_movemask_ps((__m128)renorm_mask5);
            Vv5 = _mm_shuffle_epi8(Vv5, _mm_load_si128((__m128i*)pidx[imask5]));
            sp += _mm_popcnt_u32(imask5);

            __m128i Vv6 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask6 = _mm_movemask_ps((__m128)renorm_mask6);
            sp += _mm_popcnt_u32(imask6);
            Vv6 = _mm_shuffle_epi8(Vv6, _mm_load_si128((__m128i*)pidx[imask6]));

            __m128i Vv7 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask7 = _mm_movemask_ps((__m128)renorm_mask7);
            Vv7 = _mm_shuffle_epi8(Vv7, _mm_load_si128((__m128i*)pidx[imask7]));
            sp += _mm_popcnt_u32(imask7);

            __m128i Vv8 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask8 = _mm_movemask_ps((__m128)renorm_mask8);
            sp += _mm_popcnt_u32(imask8);
            Vv8 = _mm_shuffle_epi8(Vv8, _mm_load_si128((__m128i*)pidx[imask8]));

            __m128i Yv5 = _mm_slli_epi32(Rv5, 16);
            __m128i Yv6 = _mm_slli_epi32(Rv6, 16);
            __m128i Yv7 = _mm_slli_epi32(Rv7, 16);
            __m128i Yv8 = _mm_slli_epi32(Rv8, 16);

            // y = (R[z] << 16) | V[z];
            Yv5 = _mm_or_si128(Yv5, Vv5);
            Yv6 = _mm_or_si128(Yv6, Vv6);
            Yv7 = _mm_or_si128(Yv7, Vv7);
            Yv8 = _mm_or_si128(Yv8, Vv8);

            // R[z] = c ? Y[z] : R[z];
            Rv5 = _mm_blendv_epi8(Rv5, Yv5, renorm_mask5);
            Rv6 = _mm_blendv_epi8(Rv6, Yv6, renorm_mask6);
            Rv7 = _mm_blendv_epi8(Rv7, Yv7, renorm_mask7);
            Rv8 = _mm_blendv_epi8(Rv8, Yv8, renorm_mask8);

            // Maybe just a store128 instead?
            _mm_store_si128((__m128i *)&tbuf[tidx][ 0], sv1);
            _mm_store_si128((__m128i *)&tbuf[tidx][16], sv5);
            //  *(uint64_t *)&out[i+ 0] = _mm_extract_epi64(sv1, 0);
            //  *(uint64_t *)&out[i+ 8] = _mm_extract_epi64(sv1, 1);
            //  *(uint64_t *)&out[i+16] = _mm_extract_epi64(sv5, 0);
            //  *(uint64_t *)&out[i+24] = _mm_extract_epi64(sv5, 1);

            // WRONG - need to reorder these periodically.

            i4[0]++;
            if (++tidx == 32) {
                i4[0]-=32;
                transpose_and_copy(out, i4, tbuf);
                tidx = 0;
            }

        }
        isz4 += 64;

        STORE128(Rv, R);
        STORE128(Lv, l);
        ptr = (uint8_t *)sp;

        i4[0]-=tidx;
        int T;
        for (z = 0; z < NX; z++)
            for (T = 0; T < tidx; T++)
                out[i4[z]++] = tbuf[T][z];

        // Scalar version for close to the end of in[] array so we don't
        // do SIMD loads beyond the end of the buffer
        for (; i4[0] < isz4;) {
            for (z = 0; z < NX; z++) {
                uint32_t m = R[z] & ((1u<<TF_SHIFT_O1)-1);
                uint32_t S = s3[l[z]][m];
                unsigned char c = S & 0xff;
                out[i4[z]++] = c;
                uint32_t F = S>>(TF_SHIFT_O1+8);
                R[z] = (F?F:4096) * (R[z]>>TF_SHIFT_O1) +
                    ((S>>8) & ((1u<<TF_SHIFT_O1)-1));
                RansDecRenormSafe(&R[z], &ptr, ptr_end);
                l[z] = c;
            }
        }

        // Remainder
        z = NX-1;
        for (; i4[z] < out_sz; ) {
            uint32_t m = R[z] & ((1u<<TF_SHIFT_O1)-1);
            uint32_t S = s3[l[z]][m];
            unsigned char c = S & 0xff;
            out[i4[z]++] = c;
            R[z] = (S>>(TF_SHIFT_O1+8)) * (R[z]>>TF_SHIFT_O1) +
                ((S>>8) & ((1u<<TF_SHIFT_O1)-1));
            RansDecRenormSafe(&R[z], &ptr, ptr_end);
            l[z] = c;
        }
    } else {
        // TF_SHIFT_O1 = 10
        uint16_t *sp = (uint16_t *)ptr;
        const uint32_t mask = ((1u << TF_SHIFT_O1_FAST)-1);
        __m128i maskv  = _mm_set1_epi32(mask); // set mask in all lanes
        uint8_t tbuf[32][32] __attribute__((aligned(32)));
        int tidx = 0;
        LOAD128(Rv, R);
        LOAD128(Lv, l);

        isz4 -= 64;
        for (; i4[0] < isz4 && (uint8_t *)sp+72 < ptr_end; ) {
            //for (z = 0; z < NX; z++)
            //  m[z] = R[z] & mask;
            __m128i masked1 = _mm_and_si128(Rv1, maskv);
            __m128i masked2 = _mm_and_si128(Rv2, maskv);
            __m128i masked3 = _mm_and_si128(Rv3, maskv);
            __m128i masked4 = _mm_and_si128(Rv4, maskv);

            Lv1 = _mm_slli_epi32(Lv1, TF_SHIFT_O1_FAST);
            Lv2 = _mm_slli_epi32(Lv2, TF_SHIFT_O1_FAST);
            Lv3 = _mm_slli_epi32(Lv3, TF_SHIFT_O1_FAST);
            Lv4 = _mm_slli_epi32(Lv4, TF_SHIFT_O1_FAST);
            masked1 = _mm_add_epi32(masked1, Lv1);
            masked2 = _mm_add_epi32(masked2, Lv2);
            masked3 = _mm_add_epi32(masked3, Lv3);
            masked4 = _mm_add_epi32(masked4, Lv4);

            //  S[z] = s3[l[z]][m[z]];
            __m128i Sv1 = _mm_i32gather_epi32x((int *)s3F, masked1, sizeof(*s3));
            __m128i Sv2 = _mm_i32gather_epi32x((int *)s3F, masked2, sizeof(*s3));
            __m128i Sv3 = _mm_i32gather_epi32x((int *)s3F, masked3, sizeof(*s3));
            __m128i Sv4 = _mm_i32gather_epi32x((int *)s3F, masked4, sizeof(*s3));

            //  f[z] = S[z]>>(TF_SHIFT+8);
            __m128i fv1 = _mm_srli_epi32(Sv1, TF_SHIFT_O1_FAST+8);
            __m128i fv2 = _mm_srli_epi32(Sv2, TF_SHIFT_O1_FAST+8);
            __m128i fv3 = _mm_srli_epi32(Sv3, TF_SHIFT_O1_FAST+8);
            __m128i fv4 = _mm_srli_epi32(Sv4, TF_SHIFT_O1_FAST+8);

            //  b[z] = (S[z]>>8) & mask;
            __m128i bv1 = _mm_and_si128(_mm_srli_epi32(Sv1, 8), maskv);
            __m128i bv2 = _mm_and_si128(_mm_srli_epi32(Sv2, 8), maskv);
            __m128i bv3 = _mm_and_si128(_mm_srli_epi32(Sv3, 8), maskv);
            __m128i bv4 = _mm_and_si128(_mm_srli_epi32(Sv4, 8), maskv);

            //  s[z] = S[z] & 0xff;
            __m128i sv1 = _mm_and_si128(Sv1, _mm_set1_epi32(0xff));
            __m128i sv2 = _mm_and_si128(Sv2, _mm_set1_epi32(0xff));
            __m128i sv3 = _mm_and_si128(Sv3, _mm_set1_epi32(0xff));
            __m128i sv4 = _mm_and_si128(Sv4, _mm_set1_epi32(0xff));

            //  R[z] = f[z] * (R[z] >> TF_SHIFT_O1_FAST) + b[z];
            Rv1 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv1,TF_SHIFT_O1_FAST), fv1), bv1);
            Rv2 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv2,TF_SHIFT_O1_FAST), fv2), bv2);
            Rv3 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv3,TF_SHIFT_O1_FAST), fv3), bv3);
            Rv4 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv4,TF_SHIFT_O1_FAST), fv4), bv4);

            Lv1 = sv1;
            Lv2 = sv2;
            Lv3 = sv3;
            Lv4 = sv4;

            // Tricky one:  out[i+z] = s[z];
            //             ---d---c ---b---a  sv1
            //             ---h---g ---f---e  sv2
            // packs_epi32 -h-g-f-e -d-c-b-a  sv1(2)
            // packs_epi16 ponmlkji hgfedcba  sv1(2) / sv3(4)
            sv1 = _mm_packus_epi32(sv1, sv2);
            sv3 = _mm_packus_epi32(sv3, sv4);
            sv1 = _mm_packus_epi16(sv1, sv3);

            // c =  R[z] < RANS_BYTE_L;
            // A little tricky as we only have signed comparisons.
            // See https://stackoverflow.com/questions/32945410/sse2-intrinsics-comparing-unsigned-integers

//#define _mm_cmplt_epu32_imm(a,b) _mm_andnot_si128(_mm_cmpeq_epi32(_mm_max_epu32((a),_mm_set1_epi32(b)), (a)), _mm_set1_epi32(-1));

            //#define _mm_cmplt_epu32_imm(a,b) _mm_cmpgt_epi32(_mm_set1_epi32((b)-0x80000000), _mm_xor_si128((a), _mm_set1_epi32(0x80000000)))

            __m128i renorm_mask1, renorm_mask2, renorm_mask3, renorm_mask4;
            renorm_mask1 = _mm_cmplt_epu32_imm(Rv1, RANS_BYTE_L);
            renorm_mask2 = _mm_cmplt_epu32_imm(Rv2, RANS_BYTE_L);
            renorm_mask3 = _mm_cmplt_epu32_imm(Rv3, RANS_BYTE_L);
            renorm_mask4 = _mm_cmplt_epu32_imm(Rv4, RANS_BYTE_L);

            //#define P(A,B,C,D) ((A)+((B)<<2) + ((C)<<4) + ((D)<<6))
#define P(A,B,C,D)                              \
            { A+0,A+1,A+2,A+3,                  \
              B+0,B+1,B+2,B+3,                  \
              C+0,C+1,C+2,C+3,                  \
              D+0,D+1,D+2,D+3}
#ifdef _
#undef _
#endif
#define _ 0x80
            uint8_t pidx[16][16] = {
                P(_,_,_,_),
                P(0,_,_,_),
                P(_,0,_,_),
                P(0,4,_,_),

                P(_,_,0,_),
                P(0,_,4,_),
                P(_,0,4,_),
                P(0,4,8,_),

                P(_,_,_,0),
                P(0,_,_,4),
                P(_,0,_,4),
                P(0,4,_,8),

                P(_,_,0,4),
                P(0,_,4,8),
                P(_,0,4,8),
                P(0,4,8,12),
            };
#undef _

            // Shuffle the renorm values to correct lanes and incr sp pointer
            __m128i Vv1 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask1 = _mm_movemask_ps((__m128)renorm_mask1);
            Vv1 = _mm_shuffle_epi8(Vv1, _mm_load_si128((__m128i*)pidx[imask1]));
            sp += _mm_popcnt_u32(imask1);

            __m128i Vv2 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask2 = _mm_movemask_ps((__m128)renorm_mask2);
            sp += _mm_popcnt_u32(imask2);
            Vv2 = _mm_shuffle_epi8(Vv2, _mm_load_si128((__m128i*)pidx[imask2]));

            __m128i Vv3 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask3 = _mm_movemask_ps((__m128)renorm_mask3);
            Vv3 = _mm_shuffle_epi8(Vv3, _mm_load_si128((__m128i*)pidx[imask3]));
            sp += _mm_popcnt_u32(imask3);

            __m128i Vv4 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask4 = _mm_movemask_ps((__m128)renorm_mask4);
            sp += _mm_popcnt_u32(imask4);
            Vv4 = _mm_shuffle_epi8(Vv4, _mm_load_si128((__m128i*)pidx[imask4]));

            __m128i Yv1 = _mm_slli_epi32(Rv1, 16);
            __m128i Yv2 = _mm_slli_epi32(Rv2, 16);
            __m128i Yv3 = _mm_slli_epi32(Rv3, 16);
            __m128i Yv4 = _mm_slli_epi32(Rv4, 16);

            // y = (R[z] << 16) | V[z];
            Yv1 = _mm_or_si128(Yv1, Vv1);
            Yv2 = _mm_or_si128(Yv2, Vv2);
            Yv3 = _mm_or_si128(Yv3, Vv3);
            Yv4 = _mm_or_si128(Yv4, Vv4);

            // R[z] = c ? Y[z] : R[z];
            Rv1 = _mm_blendv_epi8(Rv1, Yv1, renorm_mask1);
            Rv2 = _mm_blendv_epi8(Rv2, Yv2, renorm_mask2);
            Rv3 = _mm_blendv_epi8(Rv3, Yv3, renorm_mask3);
            Rv4 = _mm_blendv_epi8(Rv4, Yv4, renorm_mask4);

            // ------------------------------------------------------------

            //  m[z] = R[z] & mask;
            __m128i masked5 = _mm_and_si128(Rv5, maskv);
            __m128i masked6 = _mm_and_si128(Rv6, maskv);
            __m128i masked7 = _mm_and_si128(Rv7, maskv);
            __m128i masked8 = _mm_and_si128(Rv8, maskv);


            Lv5 = _mm_slli_epi32(Lv5, TF_SHIFT_O1_FAST);
            Lv6 = _mm_slli_epi32(Lv6, TF_SHIFT_O1_FAST);
            Lv7 = _mm_slli_epi32(Lv7, TF_SHIFT_O1_FAST);
            Lv8 = _mm_slli_epi32(Lv8, TF_SHIFT_O1_FAST);
            masked5 = _mm_add_epi32(masked5, Lv5);
            masked6 = _mm_add_epi32(masked6, Lv6);
            masked7 = _mm_add_epi32(masked7, Lv7);
            masked8 = _mm_add_epi32(masked8, Lv8);

            //  S[z] = s3[m[z]];
            __m128i Sv5 = _mm_i32gather_epi32x((int *)s3F, masked5, sizeof(*s3));
            __m128i Sv6 = _mm_i32gather_epi32x((int *)s3F, masked6, sizeof(*s3));
            __m128i Sv7 = _mm_i32gather_epi32x((int *)s3F, masked7, sizeof(*s3));
            __m128i Sv8 = _mm_i32gather_epi32x((int *)s3F, masked8, sizeof(*s3));

            //  f[z] = S[z]>>(TF_SHIFT_O1_FAST+8);
            __m128i fv5 = _mm_srli_epi32(Sv5, TF_SHIFT_O1_FAST+8);
            __m128i fv6 = _mm_srli_epi32(Sv6, TF_SHIFT_O1_FAST+8);
            __m128i fv7 = _mm_srli_epi32(Sv7, TF_SHIFT_O1_FAST+8);
            __m128i fv8 = _mm_srli_epi32(Sv8, TF_SHIFT_O1_FAST+8);

            //  b[z] = (S[z]>>8) & mask;
            __m128i bv5 = _mm_and_si128(_mm_srli_epi32(Sv5, 8), maskv);
            __m128i bv6 = _mm_and_si128(_mm_srli_epi32(Sv6, 8), maskv);
            __m128i bv7 = _mm_and_si128(_mm_srli_epi32(Sv7, 8), maskv);
            __m128i bv8 = _mm_and_si128(_mm_srli_epi32(Sv8, 8), maskv);

            //  s[z] = S[z] & 0xff;
            __m128i sv5 = _mm_and_si128(Sv5, _mm_set1_epi32(0xff));
            __m128i sv6 = _mm_and_si128(Sv6, _mm_set1_epi32(0xff));
            __m128i sv7 = _mm_and_si128(Sv7, _mm_set1_epi32(0xff));
            __m128i sv8 = _mm_and_si128(Sv8, _mm_set1_epi32(0xff));

            //  R[z] = f[z] * (R[z] >> TF_SHIFT_O1_FAST) + b[z];
            Rv5 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv5,TF_SHIFT_O1_FAST), fv5), bv5);
            Rv6 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv6,TF_SHIFT_O1_FAST), fv6), bv6);
            Rv7 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv7,TF_SHIFT_O1_FAST), fv7), bv7);
            Rv8 = _mm_add_epi32(
                      _mm_mullo_epi32(
                          _mm_srli_epi32(Rv8,TF_SHIFT_O1_FAST), fv8), bv8);

            Lv5 = sv5;
            Lv6 = sv6;
            Lv7 = sv7;
            Lv8 = sv8;

            // Tricky one:  out[i+z] = s[z];
            //             ---d---c ---b---a  sv1
            //             ---h---g ---f---e  sv2
            // packs_epi32 -h-g-f-e -d-c-b-a  sv1(2)
            // packs_epi16 ponmlkji hgfedcba  sv1(2) / sv3(4)
            sv5 = _mm_packus_epi32(sv5, sv6);
            sv7 = _mm_packus_epi32(sv7, sv8);
            sv5 = _mm_packus_epi16(sv5, sv7);

            // c =  R[z] < RANS_BYTE_L;
            __m128i renorm_mask5, renorm_mask6, renorm_mask7, renorm_mask8;
            renorm_mask5 = _mm_cmplt_epu32_imm(Rv5, RANS_BYTE_L);
            renorm_mask6 = _mm_cmplt_epu32_imm(Rv6, RANS_BYTE_L);
            renorm_mask7 = _mm_cmplt_epu32_imm(Rv7, RANS_BYTE_L);
            renorm_mask8 = _mm_cmplt_epu32_imm(Rv8, RANS_BYTE_L);
        
            // Shuffle the renorm values to correct lanes and incr sp pointer
            __m128i Vv5 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask5 = _mm_movemask_ps((__m128)renorm_mask5);
            Vv5 = _mm_shuffle_epi8(Vv5, _mm_load_si128((__m128i*)pidx[imask5]));
            sp += _mm_popcnt_u32(imask5);

            __m128i Vv6 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask6 = _mm_movemask_ps((__m128)renorm_mask6);
            sp += _mm_popcnt_u32(imask6);
            Vv6 = _mm_shuffle_epi8(Vv6, _mm_load_si128((__m128i*)pidx[imask6]));

            __m128i Vv7 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask7 = _mm_movemask_ps((__m128)renorm_mask7);
            Vv7 = _mm_shuffle_epi8(Vv7, _mm_load_si128((__m128i*)pidx[imask7]));
            sp += _mm_popcnt_u32(imask7);

            __m128i Vv8 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)sp));
            unsigned int imask8 = _mm_movemask_ps((__m128)renorm_mask8);
            sp += _mm_popcnt_u32(imask8);
            Vv8 = _mm_shuffle_epi8(Vv8, _mm_load_si128((__m128i*)pidx[imask8]));

            __m128i Yv5 = _mm_slli_epi32(Rv5, 16);
            __m128i Yv6 = _mm_slli_epi32(Rv6, 16);
            __m128i Yv7 = _mm_slli_epi32(Rv7, 16);
            __m128i Yv8 = _mm_slli_epi32(Rv8, 16);

            // y = (R[z] << 16) | V[z];
            Yv5 = _mm_or_si128(Yv5, Vv5);
            Yv6 = _mm_or_si128(Yv6, Vv6);
            Yv7 = _mm_or_si128(Yv7, Vv7);
            Yv8 = _mm_or_si128(Yv8, Vv8);

            // R[z] = c ? Y[z] : R[z];
            Rv5 = _mm_blendv_epi8(Rv5, Yv5, renorm_mask5);
            Rv6 = _mm_blendv_epi8(Rv6, Yv6, renorm_mask6);
            Rv7 = _mm_blendv_epi8(Rv7, Yv7, renorm_mask7);
            Rv8 = _mm_blendv_epi8(Rv8, Yv8, renorm_mask8);

            // Maybe just a store128 instead?
            _mm_store_si128((__m128i *)&tbuf[tidx][ 0], sv1);
            _mm_store_si128((__m128i *)&tbuf[tidx][16], sv5);
            //  *(uint64_t *)&out[i+ 0] = _mm_extract_epi64(sv1, 0);
            //  *(uint64_t *)&out[i+ 8] = _mm_extract_epi64(sv1, 1);
            //  *(uint64_t *)&out[i+16] = _mm_extract_epi64(sv5, 0);
            //  *(uint64_t *)&out[i+24] = _mm_extract_epi64(sv5, 1);

            // WRONG - need to reorder these periodically.

            i4[0]++;
            if (++tidx == 32) {
                i4[0]-=32;
                transpose_and_copy(out, i4, tbuf);
                tidx = 0;
            }
        }
        isz4 += 64;

        STORE128(Rv, R);
        STORE128(Lv, l);
        ptr = (uint8_t *)sp;

        i4[0]-=tidx;
        int T;
        for (z = 0; z < NX; z++)
            for (T = 0; T < tidx; T++)
                out[i4[z]++] = tbuf[T][z];

        // Scalar version for close to the end of in[] array so we don't
        // do SIMD loads beyond the end of the buffer
        for (; i4[0] < isz4;) {
            for (z = 0; z < NX; z++) {
                uint32_t m = R[z] & ((1u<<TF_SHIFT_O1_FAST)-1);
                uint32_t S = s3F[l[z]][m];
                unsigned char c = S & 0xff;
                out[i4[z]++] = c;
                R[z] = (S>>(TF_SHIFT_O1_FAST+8)) * (R[z]>>TF_SHIFT_O1_FAST) +
                    ((S>>8) & ((1u<<TF_SHIFT_O1_FAST)-1));
                RansDecRenormSafe(&R[z], &ptr, ptr_end);
                l[z] = c;
            }
        }

        // Remainder
        z = NX-1;
        for (; i4[z] < out_sz; ) {
            uint32_t m = R[z] & ((1u<<TF_SHIFT_O1_FAST)-1);
            uint32_t S = s3F[l[z]][m];
            unsigned char c = S & 0xff;
            out[i4[z]++] = c;
            R[z] = (S>>(TF_SHIFT_O1_FAST+8)) * (R[z]>>TF_SHIFT_O1_FAST) +
                ((S>>8) & ((1u<<TF_SHIFT_O1_FAST)-1));
            RansDecRenormSafe(&R[z], &ptr, ptr_end);
            l[z] = c;
        }
    }
    //fprintf(stderr, "    1 Decoded %d bytes\n", (int)(ptr-in)); //c-size

    htscodecs_tls_free(s3);
    return out;

 err:
    htscodecs_tls_free(s3);
    free(out_free);
    free(c_freq);

    return NULL;
}

#endif // HAVE_SSE4_1 and HAVE_SSSE3
