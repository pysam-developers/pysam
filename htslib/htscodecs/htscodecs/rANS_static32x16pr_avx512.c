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

/*
 * This is an AVX512 implementation of the 32-way interleaved 16-bit rANS.
 * For now it only contains an order-0 implementation.  The AVX2 code may
 * be used for order-1.
 */

#include "config.h"

#ifdef HAVE_AVX512

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <x86intrin.h>

#include "rANS_word.h"
#include "rANS_static4x16.h"
#define ROT32_SIMD
#include "rANS_static16_int.h"
#include "varint.h"
#include "utils.h"

unsigned char *rans_compress_O0_32x16_avx512(unsigned char *in,
                                             unsigned int in_size,
                                             unsigned char *out,
                                             unsigned int *out_size) {
    unsigned char *cp, *out_end;
    RansEncSymbol syms[256];
    RansState ransN[32] __attribute__((aligned(64)));
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

    // Encode statistics and build lookup tables for SIMD encoding.
    uint32_t SB[256], SA[256], SD[256], SC[256];
    for (x = j = 0; j < 256; j++) {
        if (F[j]) {
            RansEncSymbolInit(&syms[j], x, F[j], TF_SHIFT);
            SB[j] = syms[j].x_max;
            SA[j] = syms[j].rcp_freq;
            SD[j] = (syms[j].cmpl_freq<<0) | ((syms[j].rcp_shift-32)<<16);
            SC[j] = syms[j].bias;
            x += F[j];
        }
    }

    for (z = 0; z < 32; z++)
      RansEncInit(&ransN[z]);

    z = i = in_size&(32-1);
    while (z-- > 0)
      RansEncPutSymbol(&ransN[z], &ptr, &syms[in[in_size-(i-z)]]);

#define LOAD512(a,b)                                     \
    __m512i a##1 = _mm512_load_si512((__m512i *)&b[0]); \
    __m512i a##2 = _mm512_load_si512((__m512i *)&b[16]);

#define STORE512(a,b)                           \
    _mm512_store_si512((__m256i *)&b[0], a##1); \
    _mm512_store_si512((__m256i *)&b[16], a##2);

    LOAD512(Rv, ransN);

    uint16_t *ptr16 = (uint16_t *)ptr;
    for (i=(in_size &~(32-1)); i>0; i-=32) {
        uint8_t *c = &in[i-32];

        // GATHER versions
        // Much faster now we have an efficient loadu mechanism in place,
        // BUT...
        // Try this for avx2 variant too?  Better way to populate the mm256
        // regs for mix of avx2 and avx512 opcodes.
        __m256i c12 = _mm256_loadu_si256((__m256i const *)c);
        __m512i c1 = _mm512_cvtepu8_epi32(_mm256_extracti128_si256(c12,0));
        __m512i c2 = _mm512_cvtepu8_epi32(_mm256_extracti128_si256(c12,1));
#define SET512(a,b) \
        __m512i a##1 = _mm512_i32gather_epi32(c1, b, 4); \
        __m512i a##2 = _mm512_i32gather_epi32(c2, b, 4)

        SET512(xmax, SB);

        uint16_t gt_mask1 = _mm512_cmpgt_epi32_mask(Rv1, xmax1);
        int pc1 = _mm_popcnt_u32(gt_mask1);
        __m512i Rp1 = _mm512_and_si512(Rv1, _mm512_set1_epi32(0xffff));
        __m512i Rp2 = _mm512_and_si512(Rv2, _mm512_set1_epi32(0xffff));
        uint16_t gt_mask2 = _mm512_cmpgt_epi32_mask(Rv2, xmax2);
        SET512(SDv,  SD);
        int pc2 = _mm_popcnt_u32(gt_mask2);

        //Rp1 = _mm512_maskz_compress_epi32(gt_mask1, Rp1);
        Rp1 = _mm512_maskz_compress_epi32(gt_mask1, Rp1);
        Rp2 = _mm512_maskz_compress_epi32(gt_mask2, Rp2);

        _mm512_mask_cvtepi32_storeu_epi16(ptr16-pc2, (1<<pc2)-1, Rp2);
        ptr16 -= pc2;
        _mm512_mask_cvtepi32_storeu_epi16(ptr16-pc1, (1<<pc1)-1, Rp1);
        ptr16 -= pc1;

        SET512(rfv,  SA);
        Rv1 = _mm512_mask_srli_epi32(Rv1, gt_mask1, Rv1, 16);
        Rv2 = _mm512_mask_srli_epi32(Rv2, gt_mask2, Rv2, 16);

        // interleaved form of this, helps on icc a bit
        //rfv1 = _mm512_mulhi_epu32(Rv1, rfv1);
        //rfv2 = _mm512_mulhi_epu32(Rv2, rfv2);

        // Alternatives here:
        // SHIFT right/left instead of AND: (very marginally slower)
        //   rf1_hm = _mm512_and_epi32(
        //           rf1_hm, _mm512_set1_epi64((uint64_t)0xffffffff00000000));
        // vs
        //   rf1_hm = _mm512_srli_epi64(rf1_hm, 32);
        //   rf1_hm = _mm512_slli_epi64(rf1_hm, 32);
        __m512i rf1_hm = _mm512_mul_epu32(_mm512_srli_epi64(Rv1,  32),
                                          _mm512_srli_epi64(rfv1, 32));
        __m512i rf2_hm = _mm512_mul_epu32(_mm512_srli_epi64(Rv2,  32),
                                          _mm512_srli_epi64(rfv2, 32));
        __m512i rf1_lm = _mm512_srli_epi64(_mm512_mul_epu32(Rv1, rfv1), 32);
        __m512i rf2_lm = _mm512_srli_epi64(_mm512_mul_epu32(Rv2, rfv2), 32);
        rf1_hm = _mm512_and_epi32(rf1_hm,
                                  _mm512_set1_epi64((uint64_t)0xffffffff<<32));
        rf2_hm = _mm512_and_epi32(rf2_hm,
                                  _mm512_set1_epi64((uint64_t)0xffffffff<<32));
        rfv1 = _mm512_or_epi32(rf1_lm, rf1_hm);
        rfv2 = _mm512_or_epi32(rf2_lm, rf2_hm);

        // Or a pure masked blend approach, sadly slower.
        // rfv1 = _mm512_mask_blend_epi32(0x5555,
        //            _mm512_mul_epu32(
        //                _mm512_srli_epi64(Rv1, 32),
        //                _mm512_srli_epi64(rfv1, 32)),
        //            _mm512_srli_epi64(
        //                _mm512_mul_epu32(Rv1, rfv1),
        //                32));
        // rfv2 = _mm512_mask_blend_epi32(0xaaaa,
        //            _mm512_srli_epi64(
        //                _mm512_mul_epu32(Rv2, rfv2),
        //                32),
        //            _mm512_mul_epu32(
        //                _mm512_srli_epi64(Rv2, 32),
        //                _mm512_srli_epi64(rfv2, 32)));

        SET512(biasv, SC);
        __m512i shiftv1 = _mm512_srli_epi32(SDv1, 16);
        __m512i shiftv2 = _mm512_srli_epi32(SDv2, 16);

        __m512i qv1 = _mm512_srlv_epi32(rfv1, shiftv1);
        __m512i qv2 = _mm512_srlv_epi32(rfv2, shiftv2);

        qv1 = _mm512_mullo_epi32(qv1,
                  _mm512_and_si512(SDv1,_mm512_set1_epi32(0xffff)));
        qv1 = _mm512_add_epi32(qv1, biasv1);
        Rv1 = _mm512_add_epi32(Rv1, qv1);

        qv2 = _mm512_mullo_epi32(qv2,
                  _mm512_and_si512(SDv2, _mm512_set1_epi32(0xffff)));
        qv2 = _mm512_add_epi32(qv2, biasv2);
        Rv2 = _mm512_add_epi32(Rv2, qv2);
    }
    ptr = (uint8_t *)ptr16;
    STORE512(Rv, ransN);

    for (z = 32-1; z >= 0; z--)
        RansEncFlush(&ransN[z], &ptr);
    
 empty:
    // Finalise block size and return it
    *out_size = (out_end - ptr) + tab_size;

    memmove(out + tab_size, ptr, out_end-ptr);

    return out;
}

unsigned char *rans_uncompress_O0_32x16_avx512(unsigned char *in,
                                               unsigned int in_size,
                                               unsigned char *out,
                                               unsigned int out_sz) {
    if (in_size < 32*4) // 32-states at least
        return NULL;

    if (out_sz >= INT_MAX)
        return NULL; // protect against some overflow cases

    /* Load in the static tables */
    unsigned char *cp = in, *out_free = NULL;
    unsigned char *cp_end = in + in_size;
    int i;
    uint32_t s3[TOTFREQ]  __attribute__((aligned(64))); // For TF_SHIFT <= 12

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

    if (cp_end - cp < 32 * 4)
        goto err;

    int z;
    RansState Rv[32] __attribute__((aligned(64)));
    for (z = 0; z < 32; z++) {
        RansDecInit(&Rv[z], &cp);
        if (Rv[z] < RANS_BYTE_L)
            goto err;
    }

    uint16_t *sp = (uint16_t *)cp;

    int out_end = (out_sz&~(32-1));
    const uint32_t mask = (1u << TF_SHIFT)-1;

    __m512i maskv = _mm512_set1_epi32(mask); // set mask in all lanes
    __m512i R1 = _mm512_load_epi32(&Rv[0]);
    __m512i R2 = _mm512_load_epi32(&Rv[16]);

    // Start of the first loop iteration, which we do move to the end of the
    // loop for the next cycle so we can remove some of the instr. latency.
    __m512i masked1 = _mm512_and_epi32(R1, maskv);
    __m512i masked2 = _mm512_and_epi32(R2, maskv);
    __m512i S1 = _mm512_i32gather_epi32(masked1, (int *)s3, sizeof(*s3));
    __m512i S2 = _mm512_i32gather_epi32(masked2, (int *)s3, sizeof(*s3));

    uint8_t overflow[64+64] = {0};
    for (i=0; i < out_end; i+=32) {
      //for (z = 0; z < 16; z++) {

      // Protect against running off the end of in buffer.
      // We copy it to a worst-case local buffer when near the end.
      if ((uint8_t *)sp+64 > cp_end) {
        memmove(overflow, sp, cp_end - (uint8_t *)sp);
        sp = (uint16_t *)overflow;
        cp_end = overflow + sizeof(overflow);
      }

      //uint32_t S = s3[R[z] & mask];
      __m512i renorm_words1 = _mm512_cvtepu16_epi32(_mm256_loadu_si256((const __m256i *)sp)); // next 16 words

      //uint16_t f = S>>(TF_SHIFT+8), b = (S>>8) & mask;
      __m512i f1 = _mm512_srli_epi32(S1, TF_SHIFT+8);
      __m512i f2 = _mm512_srli_epi32(S2, TF_SHIFT+8);
      __m512i b1 = _mm512_and_epi32(_mm512_srli_epi32(S1, 8), maskv);
      __m512i b2 = _mm512_and_epi32(_mm512_srli_epi32(S2, 8), maskv);

      //R[z] = f * (R[z] >> TF_SHIFT) + b;
      // approx 10 cycle latency on mullo.
      R1 = _mm512_add_epi32(
               _mm512_mullo_epi32(
                   _mm512_srli_epi32(R1, TF_SHIFT), f1), b1);
      R2 = _mm512_add_epi32(
               _mm512_mullo_epi32(
                   _mm512_srli_epi32(R2, TF_SHIFT), f2), b2);

      // renorm. this is the interesting part:
      __mmask16 renorm_mask1, renorm_mask2;
      renorm_mask1=_mm512_cmplt_epu32_mask(R1, _mm512_set1_epi32(RANS_BYTE_L));
      renorm_mask2=_mm512_cmplt_epu32_mask(R2, _mm512_set1_epi32(RANS_BYTE_L));
      // advance by however many words we actually read
      sp += _mm_popcnt_u32(renorm_mask1);
      __m512i renorm_words2 = _mm512_cvtepu16_epi32(_mm256_loadu_si256(
                                      (const __m256i *)sp));

      // select masked only
      __m512i renorm_vals1, renorm_vals2;
      renorm_vals1 = _mm512_maskz_expand_epi32(renorm_mask1, renorm_words1);
      renorm_vals2 = _mm512_maskz_expand_epi32(renorm_mask2, renorm_words2);
      // shift & add selected words
      R1 = _mm512_mask_slli_epi32(R1, renorm_mask1, R1, 16);
      R2 = _mm512_mask_slli_epi32(R2, renorm_mask2, R2, 16);
      R1 = _mm512_add_epi32(R1, renorm_vals1);
      R2 = _mm512_add_epi32(R2, renorm_vals2);

      // For start of next loop iteration.  This has been moved here
      // (and duplicated to before the loop starts) so we can do something
      // with the latency period of gather, such as finishing up the
      // renorm offset and writing the results. 
      __m512i S1_ = S1; // temporary copy for use in out[]=S later
      __m512i S2_ = S2;

      masked1 = _mm512_and_epi32(R1, maskv);
      masked2 = _mm512_and_epi32(R2, maskv);
      // Gather is slow bit (half total time) - 30 cycle latency.
      S1 = _mm512_i32gather_epi32(masked1, (int *)s3, sizeof(*s3));
      S2 = _mm512_i32gather_epi32(masked2, (int *)s3, sizeof(*s3));

      // advance by however many words we actually read
      sp += _mm_popcnt_u32(renorm_mask2);

      //out[i+z] = S;
      _mm_storeu_si128((__m128i *)(out+i),    _mm512_cvtepi32_epi8(S1_));
      _mm_storeu_si128((__m128i *)(out+i+16), _mm512_cvtepi32_epi8(S2_));
    }      

    _mm512_store_epi32(&Rv[ 0], R1);
    _mm512_store_epi32(&Rv[16], R2);

    for (z = out_sz & (32-1); z-- > 0; )
      out[out_end + z] = s3[Rv[z] & mask];

    return out;

 err:
    free(out_free);
    return NULL;
}

#define TBUF8
#ifdef TBUF8
// 15% quicker overall O1 decode now due to rot32_simd below.

// NB: This uses AVX2 though and we could rewrite using AVX512 for
// further speed gains.
static inline void transpose_and_copy(uint8_t *out, int iN[32],
                                      uint8_t t[32][32]) {
//  int z;
//  for (z = 0; z < 32; z++) {
//      int k;
//      for (k = 0; k < 32; k++)
//          out[iN[z]+k] = t[k][z];
//      iN[z] += 32;
//  }

    rot32_simd(t, out, iN);
}

#else
// Implemented using AVX512 gathers.
// This is faster than a naive scalar implementation, but doesn't beat the
// AVX2 vectorised 32x32 transpose function.
static inline void transpose_and_copy_avx512(uint8_t *out, int iN[32],
                                             uint32_t t32[32][32]) {
    int z;
//  for (z = 0; z < 32; z++) {
//      int k;
//      for (k = 0; k < 32; k++)
//          out[iN[z]+k] = t32[k][z];
//      iN[z] += 32;
//  }

    
    __m512i v1 = _mm512_set_epi32(15, 14, 13, 12, 11, 10,  9,  8,
                                   7,  6,  5,  4,  3,  2,  1,  0);
    v1 = _mm512_slli_epi32(v1, 5);
    
    for (z = 0; z < 32; z++) {
        __m512i t1 = _mm512_i32gather_epi32(v1, &t32[ 0][z], 4);
        __m512i t2 = _mm512_i32gather_epi32(v1, &t32[16][z], 4);
        _mm_storeu_si128((__m128i*)(&out[iN[z]   ]), _mm512_cvtepi32_epi8(t1));
        _mm_storeu_si128((__m128i*)(&out[iN[z]+16]), _mm512_cvtepi32_epi8(t2));
        iN[z] += 32;
    }
}
#endif // TBUF

unsigned char *rans_compress_O1_32x16_avx512(unsigned char *in,
                                             unsigned int in_size,
                                             unsigned char *out,
                                             unsigned int *out_size) {
    unsigned char *cp, *out_end, *out_free = NULL;
    unsigned int tab_size;
    uint32_t bound = rans_compress_bound_4x16(in_size,1)-20;
    int z;
    RansState ransN[32] __attribute__((aligned(64)));

    if (in_size < 32) // force O0 instead
        return NULL;

    if (!out) {
        *out_size = bound;
        out = malloc(*out_size);
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
    int shift = encode_freq1(in, in_size, 32, syms, &cp); 
    if (shift < 0) {
        free(out_free);
        htscodecs_tls_free(syms);
        return NULL;
    }
    tab_size = cp - out;

    for (z = 0; z < 32; z++)
      RansEncInit(&ransN[z]);

    uint8_t* ptr = out_end;

    int iN[32], isz4 = in_size/32;
    for (z = 0; z < 32; z++)
        iN[z] = (z+1)*isz4-2;

    uint32_t lN[32] __attribute__((aligned(64)));
    for (z = 0; z < 32; z++)
        lN[z] = in[iN[z]+1];

    // Deal with the remainder
    z = 32-1;
    lN[z] = in[in_size-1];
    for (iN[z] = in_size-2; iN[z] > 32*isz4-2; iN[z]--) {
        unsigned char c = in[iN[z]];
        RansEncPutSymbol(&ransN[z], &ptr, &syms[c][lN[z]]);
        lN[z] = c;
    }

    LOAD512(Rv, ransN);

    uint16_t *ptr16 = (uint16_t *)ptr;
    __m512i last2 = _mm512_set_epi32(lN[31], lN[30], lN[29], lN[28],
                                     lN[27], lN[26], lN[25], lN[24],
                                     lN[23], lN[22], lN[21], lN[20],
                                     lN[19], lN[18], lN[17], lN[16]);
    __m512i last1 = _mm512_set_epi32(lN[15], lN[14], lN[13], lN[12],
                                     lN[11], lN[10], lN[ 9], lN[ 8],
                                     lN[ 7], lN[ 6], lN[ 5], lN[ 4],
                                     lN[ 3], lN[ 2], lN[ 1], lN[ 0]);
    
    __m512i iN2 = _mm512_set_epi32(iN[31], iN[30], iN[29], iN[28],
                                   iN[27], iN[26], iN[25], iN[24],
                                   iN[23], iN[22], iN[21], iN[20],
                                   iN[19], iN[18], iN[17], iN[16]);
    __m512i iN1 = _mm512_set_epi32(iN[15], iN[14], iN[13], iN[12],
                                   iN[11], iN[10], iN[ 9], iN[ 8],
                                   iN[ 7], iN[ 6], iN[ 5], iN[ 4],
                                   iN[ 3], iN[ 2], iN[ 1], iN[ 0]);

    __m512i c1 = _mm512_i32gather_epi32(iN1, in, 1);
    __m512i c2 = _mm512_i32gather_epi32(iN2, in, 1);

    for (; iN[0] >= 0; iN[0]--) {
        // Note, consider doing the same approach for the AVX2 encoder.
        // Maybe we can also get gather working well there?
        // Gather here is still a major latency bottleneck, consuming
        // around 40% of CPU cycles overall.

        // FIXME: maybe we need to cope with in[31] read over-flow
        // on loop cycles 0, 1, 2 where gather reads 32-bits instead of
        // 8 bits.  Use set instead there on c2?
        c1 = _mm512_and_si512(c1, _mm512_set1_epi32(0xff));
        c2 = _mm512_and_si512(c2, _mm512_set1_epi32(0xff));

        // index into syms[0][0] array, used for x_max, rcp_freq, and bias
        __m512i vidx1 = _mm512_slli_epi32(c1, 8);
        __m512i vidx2 = _mm512_slli_epi32(c2, 8);
        vidx1 = _mm512_add_epi32(vidx1, last1);
        vidx2 = _mm512_add_epi32(vidx2, last2);
        vidx1 = _mm512_slli_epi32(vidx1, 2);
        vidx2 = _mm512_slli_epi32(vidx2, 2);

        // ------------------------------------------------------------
        //      for (z = NX-1; z >= 0; z--) {
        //          if (ransN[z] >= x_max[z]) {
        //              *--ptr16 = ransN[z] & 0xffff;
        //              ransN[z] >>= 16;
        //          }
        //      }

#define SET512x(a,x) \
        __m512i a##1 = _mm512_i32gather_epi32(vidx1, &syms[0][0].x, 4); \
        __m512i a##2 = _mm512_i32gather_epi32(vidx2, &syms[0][0].x, 4)

        // Start of next loop, moved here to remove latency.
        // last[z] = c[z]
        // iN[z]--
        // c[z] = in[iN[z]]
        last1 = c1;
        last2 = c2;
        iN1 = _mm512_sub_epi32(iN1, _mm512_set1_epi32(1));
        iN2 = _mm512_sub_epi32(iN2, _mm512_set1_epi32(1));
        c1 = _mm512_i32gather_epi32(iN1, in, 1);
        c2 = _mm512_i32gather_epi32(iN2, in, 1);

        SET512x(xmax, x_max); // high latency

        uint16_t gt_mask1 = _mm512_cmpgt_epi32_mask(Rv1, xmax1);
        int pc1 = _mm_popcnt_u32(gt_mask1);
        __m512i Rp1 = _mm512_and_si512(Rv1, _mm512_set1_epi32(0xffff));
        __m512i Rp2 = _mm512_and_si512(Rv2, _mm512_set1_epi32(0xffff));
        uint16_t gt_mask2 = _mm512_cmpgt_epi32_mask(Rv2, xmax2);
        SET512x(SDv, cmpl_freq); // good
        int pc2 = _mm_popcnt_u32(gt_mask2);

        Rp1 = _mm512_maskz_compress_epi32(gt_mask1, Rp1);
        Rp2 = _mm512_maskz_compress_epi32(gt_mask2, Rp2);

        _mm512_mask_cvtepi32_storeu_epi16(ptr16-pc2, (1<<pc2)-1, Rp2);
        ptr16 -= pc2;
        _mm512_mask_cvtepi32_storeu_epi16(ptr16-pc1, (1<<pc1)-1, Rp1);
        ptr16 -= pc1;

        Rv1 = _mm512_mask_srli_epi32(Rv1, gt_mask1, Rv1, 16);
        Rv2 = _mm512_mask_srli_epi32(Rv2, gt_mask2, Rv2, 16);

        // ------------------------------------------------------------
        // uint32_t q = (uint32_t) (((uint64_t)ransN[z] * rcp_freq[z])
        //                          >> rcp_shift[z]);
        // ransN[z] = ransN[z] + bias[z] + q * cmpl_freq[z];
        SET512x(rfv, rcp_freq); // good-ish

        __m512i rf1_hm = _mm512_mul_epu32(_mm512_srli_epi64(Rv1, 32),
                                          _mm512_srli_epi64(rfv1, 32));
        __m512i rf2_hm = _mm512_mul_epu32(_mm512_srli_epi64(Rv2, 32),
                                          _mm512_srli_epi64(rfv2, 32));
        __m512i rf1_lm = _mm512_srli_epi64(_mm512_mul_epu32(Rv1, rfv1), 32);
        __m512i rf2_lm = _mm512_srli_epi64(_mm512_mul_epu32(Rv2, rfv2), 32);

        const __m512i top32 = _mm512_set1_epi64((uint64_t)0xffffffff00000000);
        rf1_hm = _mm512_and_epi32(rf1_hm, top32);
        rf2_hm = _mm512_and_epi32(rf2_hm, top32);
        rfv1 = _mm512_or_epi32(rf1_lm, rf1_hm);
        rfv2 = _mm512_or_epi32(rf2_lm, rf2_hm);

        SET512x(biasv, bias); // good
        __m512i shiftv1 = _mm512_srli_epi32(SDv1, 16);
        __m512i shiftv2 = _mm512_srli_epi32(SDv2, 16);

        shiftv1 = _mm512_sub_epi32(shiftv1, _mm512_set1_epi32(32));
        shiftv2 = _mm512_sub_epi32(shiftv2, _mm512_set1_epi32(32));

        __m512i qv1 = _mm512_srlv_epi32(rfv1, shiftv1);
        __m512i qv2 = _mm512_srlv_epi32(rfv2, shiftv2);

        const __m512i bot16 = _mm512_set1_epi32(0xffff);
        qv1 = _mm512_mullo_epi32(qv1, _mm512_and_si512(SDv1, bot16));
        qv2 = _mm512_mullo_epi32(qv2, _mm512_and_si512(SDv2, bot16));

        qv1 = _mm512_add_epi32(qv1, biasv1);
        Rv1 = _mm512_add_epi32(Rv1, qv1);

        qv2 = _mm512_add_epi32(qv2, biasv2);
        Rv2 = _mm512_add_epi32(Rv2, qv2);
    }

    STORE512(Rv, ransN);
    STORE512(last, lN);

    ptr = (uint8_t *)ptr16;

    for (z = 32-1; z>=0; z--)
        RansEncPutSymbol(&ransN[z], &ptr, &syms[0][lN[z]]);

    for (z = 32-1; z >= 0; z--)
        RansEncFlush(&ransN[z], &ptr);

    // Finalise block size and return it
    *out_size = (out_end - ptr) + tab_size;

//    cp = out;
//    *cp++ = (in_size>> 0) & 0xff;
//    *cp++ = (in_size>> 8) & 0xff;
//    *cp++ = (in_size>>16) & 0xff;
//    *cp++ = (in_size>>24) & 0xff;

    memmove(out + tab_size, ptr, out_end-ptr);

    htscodecs_tls_free(syms);
    return out;
}

#define NX 32
unsigned char *rans_uncompress_O1_32x16_avx512(unsigned char *in,
                                               unsigned int in_size,
                                               unsigned char *out,
                                               unsigned int out_sz) {
    if (in_size < NX*4) // 4-states at least
        return NULL;

    if (out_sz >= INT_MAX)
        return NULL; // protect against some overflow cases

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
        if (!(c_freq = rans_uncompress_O0_4x16(cp, c_freq_sz, NULL,
                                               u_freq_sz)))
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

    RansState R[NX] __attribute__((aligned(64)));
    uint8_t *ptr = cp, *ptr_end = in + in_size;
    int z;
    for (z = 0; z < NX; z++) {
        RansDecInit(&R[z], &ptr);
        if (R[z] < RANS_BYTE_L)
            goto err;
    }

    int isz4 = out_sz/NX;
    int iN[NX], lN[NX] __attribute__((aligned(64))) = {0};
    for (z = 0; z < NX; z++)
        iN[z] = z*isz4;

    uint16_t *sp = (uint16_t *)ptr;
    const uint32_t mask = (1u << shift)-1;

    __m512i _maskv  = _mm512_set1_epi32(mask);
    LOAD512(_Rv, R);
    LOAD512(_Lv, lN);

#ifdef TBUF8
    union {
        unsigned char tbuf[32][32];
        uint64_t tbuf64[32][4];
    } u __attribute__((aligned(32)));
#else
    uint32_t tbuf[32][32];
#endif

    unsigned int tidx = 0;

    if (shift == TF_SHIFT_O1) {
        isz4 -= 64;
        for (; iN[0] < isz4 && (uint8_t *)sp+64 < ptr_end; ) {
            // m[z] = R[z] & mask;
            __m512i _masked1 = _mm512_and_si512(_Rv1, _maskv);
            __m512i _masked2 = _mm512_and_si512(_Rv2, _maskv);

            //  S[z] = s3[lN[z]][m[z]];
            _Lv1 = _mm512_slli_epi32(_Lv1, TF_SHIFT_O1);
            _Lv2 = _mm512_slli_epi32(_Lv2, TF_SHIFT_O1);

            _masked1 = _mm512_add_epi32(_masked1, _Lv1);
            _masked2 = _mm512_add_epi32(_masked2, _Lv2);

            // This is the biggest bottleneck
            __m512i _Sv1 = _mm512_i32gather_epi32(_masked1, (int *)&s3F[0][0],
                                                  sizeof(s3F[0][0]));
            __m512i _Sv2 = _mm512_i32gather_epi32(_masked2, (int *)&s3F[0][0],
                                                  sizeof(s3F[0][0]));

            //  f[z] = S[z]>>(TF_SHIFT_O1+8);
            __m512i _fv1 = _mm512_srli_epi32(_Sv1, TF_SHIFT_O1+8);
            __m512i _fv2 = _mm512_srli_epi32(_Sv2, TF_SHIFT_O1+8);

            //  b[z] = (S[z]>>8) & mask;
            __m512i _bv1 = _mm512_and_si512(_mm512_srli_epi32(_Sv1,8), _maskv);
            __m512i _bv2 = _mm512_and_si512(_mm512_srli_epi32(_Sv2,8), _maskv);

            //  s[z] = S[z] & 0xff;
            __m512i _sv1 = _mm512_and_si512(_Sv1, _mm512_set1_epi32(0xff));
            __m512i _sv2 = _mm512_and_si512(_Sv2, _mm512_set1_epi32(0xff));

            // A maximum frequency of 4096 doesn't fit in our s3 array.
            // as it's 12 bit + 12 bit + 8 bit.  It wraps around to zero.
            // (We don't have this issue for TOTFREQ_O1_FAST.)
            //
            // Solution 1 is to change to spec to forbid freq of 4096.
            // Easy hack is to add an extra symbol so it sums correctly.
            // => 572 MB/s on q40 (deskpro).
            //
            // Solution 2 implemented here is to look for the wrap around
            // and fix it.
            // => 556 MB/s on q40
            // cope with max freq of 4096.  Only 3% hit
            __m512i max_freq = _mm512_set1_epi32(TOTFREQ_O1);
            __m512i zero = _mm512_setzero_si512();
            __mmask16 cmp1 = _mm512_cmpeq_epi32_mask(_fv1, zero);
            __mmask16 cmp2 = _mm512_cmpeq_epi32_mask(_fv2, zero);
            _fv1 = _mm512_mask_blend_epi32(cmp1, _fv1, max_freq);
            _fv2 = _mm512_mask_blend_epi32(cmp2, _fv2, max_freq);

            //  R[z] = f[z] * (R[z] >> TF_SHIFT_O1) + b[z];
            _Rv1 = _mm512_add_epi32(
                       _mm512_mullo_epi32(
                           _mm512_srli_epi32(_Rv1,TF_SHIFT_O1), _fv1), _bv1);
            _Rv2 = _mm512_add_epi32(
                       _mm512_mullo_epi32(
                           _mm512_srli_epi32(_Rv2,TF_SHIFT_O1), _fv2), _bv2);

            //for (z = 0; z < NX; z++) lN[z] = c[z];
            _Lv1 = _sv1;
            _Lv2 = _sv2;

            // RansDecRenorm(&R[z], &ptr);
            __m512i _renorm_mask1 = _mm512_xor_si512(_Rv1,
                                        _mm512_set1_epi32(0x80000000));
            __m512i _renorm_mask2 = _mm512_xor_si512(_Rv2,
                                        _mm512_set1_epi32(0x80000000));

            int _imask1 =_mm512_cmpgt_epi32_mask
                (_mm512_set1_epi32(RANS_BYTE_L-0x80000000), _renorm_mask1);
            int _imask2 = _mm512_cmpgt_epi32_mask
                (_mm512_set1_epi32(RANS_BYTE_L-0x80000000), _renorm_mask2);

            __m512i renorm_words1 = _mm512_cvtepu16_epi32
                (_mm256_loadu_si256((const __m256i *)sp));
            sp += _mm_popcnt_u32(_imask1);

            __m512i renorm_words2 = _mm512_cvtepu16_epi32
                (_mm256_loadu_si256((const __m256i *)sp));
            sp += _mm_popcnt_u32(_imask2);

            __m512i _renorm_vals1 =
                _mm512_maskz_expand_epi32(_imask1, renorm_words1);
            __m512i _renorm_vals2 =
                _mm512_maskz_expand_epi32(_imask2, renorm_words2);

            _Rv1 = _mm512_mask_slli_epi32(_Rv1, _imask1, _Rv1, 16);
            _Rv2 = _mm512_mask_slli_epi32(_Rv2, _imask2, _Rv2, 16);

            _Rv1 = _mm512_add_epi32(_Rv1, _renorm_vals1);
            _Rv2 = _mm512_add_epi32(_Rv2, _renorm_vals2);

#ifdef TBUF8
            _mm_storeu_si128((__m128i *)(&u.tbuf64[tidx][0]),
                             _mm512_cvtepi32_epi8(_Sv1)); // or _sv1?
            _mm_storeu_si128((__m128i *)(&u.tbuf64[tidx][2]),
                             _mm512_cvtepi32_epi8(_Sv2));
#else
            _mm512_storeu_si512((__m512i *)(&tbuf[tidx][ 0]), _sv1);
            _mm512_storeu_si512((__m512i *)(&tbuf[tidx][16]), _sv2);
#endif

            iN[0]++;
            if (++tidx == 32) {
                iN[0]-=32;

                // We have tidx[x][y] which we want to store in
                // memory in out[y][z] instead.  This is an unrolled
                // transposition.
#ifdef TBUF8
                transpose_and_copy(out, iN, u.tbuf);
#else
                transpose_and_copy_avx512(out, iN, tbuf);
#endif
                tidx = 0;
            }
        }
        isz4 += 64;

        STORE512(_Rv, R);
        STORE512(_Lv, lN);
        ptr = (uint8_t *)sp;

        if (1) {
            iN[0]-=tidx;
            int T;
            for (z = 0; z < NX; z++)
                for (T = 0; T < tidx; T++)
#ifdef TBUF8
                    out[iN[z]++] = u.tbuf[T][z];
#else
                    out[iN[z]++] = tbuf[T][z];
#endif
        }

        // Scalar version for close to the end of in[] array so we don't
        // do SIMD loads beyond the end of the buffer
        for (; iN[0] < isz4;) {
            for (z = 0; z < NX; z++) {
                uint32_t m = R[z] & ((1u<<TF_SHIFT_O1)-1);
                uint32_t S = s3[lN[z]][m];
                unsigned char c = S & 0xff;
                out[iN[z]++] = c;
                uint32_t F = S>>(TF_SHIFT_O1+8);
                R[z] = (F?F:4096) * (R[z]>>TF_SHIFT_O1) +
                    ((S>>8) & ((1u<<TF_SHIFT_O1)-1));
                RansDecRenormSafe(&R[z], &ptr, ptr_end);
                lN[z] = c;
            }
        }

        // Remainder
        z = NX-1;
        for (; iN[z] < out_sz; ) {
            uint32_t m = R[z] & ((1u<<TF_SHIFT_O1)-1);
            uint32_t S = s3[lN[z]][m];
            unsigned char c = S & 0xff;
            out[iN[z]++] = c;
            uint32_t F = S>>(TF_SHIFT_O1+8);
            R[z] = (F?F:4096) * (R[z]>>TF_SHIFT_O1) +
                ((S>>8) & ((1u<<TF_SHIFT_O1)-1));
            RansDecRenormSafe(&R[z], &ptr, ptr_end);
            lN[z] = c;
        }
    } else {
        // TF_SHIFT_O1_FAST.  This is the most commonly used variant.

        // SIMD version ends decoding early as it reads at most 64 bytes
        // from input via 4 vectorised loads.
        isz4 -= 64;
        for (; iN[0] < isz4 && (uint8_t *)sp+64 < ptr_end; ) {
            // m[z] = R[z] & mask;
            __m512i _masked1 = _mm512_and_si512(_Rv1, _maskv);
            __m512i _masked2 = _mm512_and_si512(_Rv2, _maskv);

            //  S[z] = s3[lN[z]][m[z]];
            _Lv1 = _mm512_slli_epi32(_Lv1, TF_SHIFT_O1_FAST);
            _Lv2 = _mm512_slli_epi32(_Lv2, TF_SHIFT_O1_FAST);

            _masked1 = _mm512_add_epi32(_masked1, _Lv1);
            _masked2 = _mm512_add_epi32(_masked2, _Lv2);

            // This is the biggest bottleneck
            __m512i _Sv1 = _mm512_i32gather_epi32(_masked1, (int *)&s3F[0][0],
                                                  sizeof(s3F[0][0]));
            __m512i _Sv2 = _mm512_i32gather_epi32(_masked2, (int *)&s3F[0][0],
                                                  sizeof(s3F[0][0]));

            //  f[z] = S[z]>>(TF_SHIFT_O1+8);
            __m512i _fv1 = _mm512_srli_epi32(_Sv1, TF_SHIFT_O1_FAST+8);
            __m512i _fv2 = _mm512_srli_epi32(_Sv2, TF_SHIFT_O1_FAST+8);

            //  b[z] = (S[z]>>8) & mask;
            __m512i _bv1 = _mm512_and_si512(_mm512_srli_epi32(_Sv1,8), _maskv);
            __m512i _bv2 = _mm512_and_si512(_mm512_srli_epi32(_Sv2,8), _maskv);

            //  s[z] = S[z] & 0xff;
            __m512i _sv1 = _mm512_and_si512(_Sv1, _mm512_set1_epi32(0xff));
            __m512i _sv2 = _mm512_and_si512(_Sv2, _mm512_set1_epi32(0xff));

            //  R[z] = f[z] * (R[z] >> TF_SHIFT_O1) + b[z];
            _Rv1 = _mm512_add_epi32(
                       _mm512_mullo_epi32(
                           _mm512_srli_epi32(_Rv1,TF_SHIFT_O1_FAST),
                           _fv1), _bv1);
            _Rv2 = _mm512_add_epi32(
                       _mm512_mullo_epi32(
                           _mm512_srli_epi32(_Rv2,TF_SHIFT_O1_FAST),
                           _fv2), _bv2);

            //for (z = 0; z < NX; z++) lN[z] = c[z];
            _Lv1 = _sv1;
            _Lv2 = _sv2;

            // RansDecRenorm(&R[z], &ptr);
            __m512i _renorm_mask1 = _mm512_xor_si512(_Rv1,
                                        _mm512_set1_epi32(0x80000000));
            __m512i _renorm_mask2 = _mm512_xor_si512(_Rv2,
                                        _mm512_set1_epi32(0x80000000));

            int _imask1 =_mm512_cmpgt_epi32_mask
                (_mm512_set1_epi32(RANS_BYTE_L-0x80000000), _renorm_mask1);
            int _imask2 = _mm512_cmpgt_epi32_mask
                (_mm512_set1_epi32(RANS_BYTE_L-0x80000000), _renorm_mask2);

            __m512i renorm_words1 = _mm512_cvtepu16_epi32
                (_mm256_loadu_si256((const __m256i *)sp));
            sp += _mm_popcnt_u32(_imask1);

            __m512i renorm_words2 = _mm512_cvtepu16_epi32
                (_mm256_loadu_si256((const __m256i *)sp));
            sp += _mm_popcnt_u32(_imask2);

            __m512i _renorm_vals1 =
                _mm512_maskz_expand_epi32(_imask1, renorm_words1);
            __m512i _renorm_vals2 =
                _mm512_maskz_expand_epi32(_imask2, renorm_words2);

            _Rv1 = _mm512_mask_slli_epi32(_Rv1, _imask1, _Rv1, 16);
            _Rv2 = _mm512_mask_slli_epi32(_Rv2, _imask2, _Rv2, 16);

            _Rv1 = _mm512_add_epi32(_Rv1, _renorm_vals1);
            _Rv2 = _mm512_add_epi32(_Rv2, _renorm_vals2);

#ifdef TBUF8
            _mm_storeu_si128((__m128i *)(&u.tbuf64[tidx][0]),
                             _mm512_cvtepi32_epi8(_Sv1)); // or _sv1?
            _mm_storeu_si128((__m128i *)(&u.tbuf64[tidx][2]),
                             _mm512_cvtepi32_epi8(_Sv2));
#else
            _mm512_storeu_si512((__m512i *)(&tbuf[tidx][ 0]), _sv1);
            _mm512_storeu_si512((__m512i *)(&tbuf[tidx][16]), _sv2);
#endif

            iN[0]++;
            if (++tidx == 32) {
                iN[0]-=32;
#ifdef TBUF8
                transpose_and_copy(out, iN, u.tbuf);
#else
                transpose_and_copy_avx512(out, iN, tbuf);
#endif
                tidx = 0;
            }
        }
        isz4 += 64;

        STORE512(_Rv, R);
        STORE512(_Lv, lN);
        ptr = (uint8_t *)sp;

        if (1) {
            iN[0]-=tidx;
            int T;
            for (z = 0; z < NX; z++)
                for (T = 0; T < tidx; T++)
#ifdef TBUF8
                    out[iN[z]++] = u.tbuf[T][z];
#else
                    out[iN[z]++] = tbuf[T][z];
#endif
        }

        // Scalar version for close to the end of in[] array so we don't
        // do SIMD loads beyond the end of the buffer
        for (; iN[0] < isz4;) {
            for (z = 0; z < NX; z++) {
                uint32_t m = R[z] & ((1u<<TF_SHIFT_O1_FAST)-1);
                uint32_t S = s3F[lN[z]][m];
                unsigned char c = S & 0xff;
                out[iN[z]++] = c;
                R[z] = (S>>(TF_SHIFT_O1_FAST+8)) * (R[z]>>TF_SHIFT_O1_FAST) +
                    ((S>>8) & ((1u<<TF_SHIFT_O1_FAST)-1));
                RansDecRenormSafe(&R[z], &ptr, ptr_end);
                lN[z] = c;
            }
        }

        // Remainder
        z = NX-1;
        for (; iN[z] < out_sz; ) {
            uint32_t m = R[z] & ((1u<<TF_SHIFT_O1_FAST)-1);
            uint32_t S = s3F[lN[z]][m];
            unsigned char c = S & 0xff;
            out[iN[z]++] = c;
            R[z] = (S>>(TF_SHIFT_O1_FAST+8)) * (R[z]>>TF_SHIFT_O1_FAST) +
                ((S>>8) & ((1u<<TF_SHIFT_O1_FAST)-1));
            RansDecRenormSafe(&R[z], &ptr, ptr_end);
            lN[z] = c;
        }
    }

    htscodecs_tls_free(s3);
    return out;

 err:
    htscodecs_tls_free(s3);
    free(out_free);
    free(c_freq);

    return NULL;
}

#endif // HAVE_AVX512
