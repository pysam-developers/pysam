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
#if defined(__ARM_NEON) && defined(__aarch64__)
#include <arm_neon.h>

#include <limits.h>
#include <assert.h>

#include "rANS_word.h"
#include "rANS_static4x16.h"
#include "rANS_static16_int.h"
#include "varint.h"
#include "utils.h"

#define NX 32

// TODO: get access to MVE architecture so we can tune for the newer
// SIMD instructions.
//
// #if __ARM_FEATURE_MVE & 1
// #include <arm_mve.h> // Helium, eg for use of vcreateq_u32
// #endif

#define _ 99
static uint8x8_t vtab[16] = {
    {_,_,  _,_,  _,_,  _,_ },
    {_,_,  _,_,  _,_, 12,13},
    {_,_,  _,_,  _,_,  8,9 },
    {_,_,  _,_,  8,9, 12,13},
    {_,_,  _,_,  _,_,  4,5 },
    {_,_,  _,_,  4,5, 12,13},
    {_,_,  _,_,  4,5,  8,9 },
    {_,_,  4,5,  8,9, 12,13},
    {_,_,  _,_,  _,_,  0,1 },
    {_,_,  _,_,  0,1, 12,13},
    {_,_,  _,_,  0,1,  8,9 },
    {_,_,  0,1 , 8,9, 12,13},
    {_,_,  _,_,  0,1,  4,5 },
    {_,_,  0,1,  4,5, 12,13},
    {_,_,  0,1,  4,5,  8,9 },
    {0,1,  4,5,  8,9, 12,13},
};
#undef _

unsigned char *rans_compress_O0_32x16_neon(unsigned char *in,
                                           unsigned int in_size,
                                           unsigned char *out,
                                           unsigned int *out_size) {
    unsigned char *cp, *out_end;
    RansEncSymbol syms[256];
    RansState R[NX];
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
      RansEncInit(&R[z]);

    z = i = in_size&(NX-1);
    while (z-- > 0)
      RansEncPutSymbol(&R[z], &ptr, &syms[in[in_size-(i-z)]]);

    for (i=(in_size &~(NX-1)); i>0; i-=NX) {
//      // Scalar equivalent
//      for (z = NX-1; z >= 0; z-=4) {
//          // 327 / 272
//          RansEncSymbol *s0 = &syms[in[i-(NX-z+0)]];
//          RansEncSymbol *s1 = &syms[in[i-(NX-z+1)]];
//          RansEncSymbol *s2 = &syms[in[i-(NX-z+2)]];
//          RansEncSymbol *s3 = &syms[in[i-(NX-z+3)]];
//
//          RansEncPutSymbol(&R[z-0], &ptr, s0);
//          RansEncPutSymbol(&R[z-1], &ptr, s1);
//          RansEncPutSymbol(&R[z-2], &ptr, s2);
//          RansEncPutSymbol(&R[z-3], &ptr, s3);
//      }

        // SIMD with 16-way unrolling
        for (z = NX-1; z >= 0; z-=8) {
            RansEncSymbol *s0 = &syms[in[i-(NX-z+0)]];
            RansEncSymbol *s1 = &syms[in[i-(NX-z+1)]];
            RansEncSymbol *s2 = &syms[in[i-(NX-z+2)]];
            RansEncSymbol *s3 = &syms[in[i-(NX-z+3)]];

            RansEncSymbol *s4 = &syms[in[i-(NX-z+4)]];
            RansEncSymbol *s5 = &syms[in[i-(NX-z+5)]];
            RansEncSymbol *s6 = &syms[in[i-(NX-z+6)]];
            RansEncSymbol *s7 = &syms[in[i-(NX-z+7)]];

            uint32x4_t Rv1 = vld1q_u32(&R[z-3]);
            uint32x4_t Rv2 = vld1q_u32(&R[z-7]);

            // Sym bit sizes = 128bits
            // 32: x_max
            // 32: rcp_freq
            // 32: bias
            // 16: cmpl_freq
            // 16: rcp_shift

            // Load and shuffle around
            //   A <---Xmax---><---RFreq--><---Bias---><-cf-><-rs->
            //   B <---Xmax---><---RFreq--><---Bias---><-cf-><-rs->
            //   C <---Xmax---><---RFreq--><---Bias---><-cf-><-rs->
            //   D <---Xmax---><---RFreq--><---Bias---><-cf-><-rs->
            // vtrn1q_u32 vtrn2q_u32  (A1 = A+B)
            //   A1  <---Xmax---><---Xmax---><---Bias---><---Bias--->
            //   C1  <---Xmax---><---Xmax---><---Bias---><---Bias--->
            //   A2  <---RFreq--><---RFreq--><-cf-><-rs-><-cf-><-rs->
            //   C2  <---RFreq--><---RFreq--><-cf-><-rs-><-cf-><-rs->
            // vtrn1q_u64 vtrn2q_u64  (A11 = A1+C1)
            //   A11 <---Xmax---><---Xmax---><---Xmax---><---Xmax--->
            //   A12 <---Bias---><---Bias---><---Bias---><---Bias--->
            //   A21 <---RFreq--><---RFreq--><---RFreq--><---RFreq-->
            //   A22 <-cf-><-rs-><-cf-><-rs-><-cf-><-rs-><-cf-><-rs->
            uint32x4_t A_1 = vld1q_u32((void *)s3);
            uint32x4_t B_1 = vld1q_u32((void *)s2);
            uint32x4_t C_1 = vld1q_u32((void *)s1);
            uint32x4_t D_1 = vld1q_u32((void *)s0);

            uint32x4_t A1_1 = vtrn1q_u32(A_1, B_1);
            uint32x4_t C1_1 = vtrn1q_u32(C_1, D_1);
            uint32x4_t A2_1 = vtrn2q_u32(A_1, B_1);
            uint32x4_t C2_1 = vtrn2q_u32(C_1, D_1);

#define u32_u64(x) vreinterpretq_u32_u64((x))
#define u64_u32(x) vreinterpretq_u64_u32((x))
            uint32x4_t Xmaxv1=u32_u64(vtrn1q_u64(u64_u32(A1_1),u64_u32(C1_1)));
            uint32x4_t Biasv1=u32_u64(vtrn2q_u64(u64_u32(A1_1),u64_u32(C1_1)));
            uint32x4_t RFv1  =u32_u64(vtrn1q_u64(u64_u32(A2_1),u64_u32(C2_1)));
            uint32x4_t FSv1  =u32_u64(vtrn2q_u64(u64_u32(A2_1),u64_u32(C2_1)));

            uint32x4_t A_2 = vld1q_u32((void *)s7);
            uint32x4_t B_2 = vld1q_u32((void *)s6);
            uint32x4_t C_2 = vld1q_u32((void *)s5);
            uint32x4_t D_2 = vld1q_u32((void *)s4);

            uint32x4_t A1_2 = vtrn1q_u32(A_2, B_2);
            uint32x4_t C1_2 = vtrn1q_u32(C_2, D_2);
            uint32x4_t A2_2 = vtrn2q_u32(A_2, B_2);
            uint32x4_t C2_2 = vtrn2q_u32(C_2, D_2);

            uint32x4_t Xmaxv2=u32_u64(vtrn1q_u64(u64_u32(A1_2),u64_u32(C1_2)));
            uint32x4_t Biasv2=u32_u64(vtrn2q_u64(u64_u32(A1_2),u64_u32(C1_2)));
            uint32x4_t RFv2  =u32_u64(vtrn1q_u64(u64_u32(A2_2),u64_u32(C2_2)));
            uint32x4_t FSv2  =u32_u64(vtrn2q_u64(u64_u32(A2_2),u64_u32(C2_2)));
            
            // Turn multi R<xmax checks into a bit-field (imask)
            uint32x4_t Cv1 = vcgtq_u32(Rv1, Xmaxv1);
            uint32x4_t Cv2 = vcgtq_u32(Rv2, Xmaxv2);

            uint32x4_t bit = {8,4,2,1};
            uint32_t imask1 = vaddvq_u32(vandq_u32(Cv1, bit));
            uint32_t imask2 = vaddvq_u32(vandq_u32(Cv2, bit));

            // Select low 16-bits from Rv based on imask, using tbl
            uint8x8_t norm1, norm2;
            norm1 = vqtbl1_u8(vreinterpretq_u8_u32(Rv1),vtab[imask1]);
            norm2 = vqtbl1_u8(vreinterpretq_u8_u32(Rv2),vtab[imask2]);

            static int nbits[16] = { 0,2,2,4, 2,4,4,6, 2,4,4,6, 4,6,6,8 };
            vst1_u8(ptr-8, norm1);   ptr -= nbits[imask1];
            vst1_u8(ptr-8, norm2);   ptr -= nbits[imask2];

            // R' = R>>16
            uint32x4_t Rv1_r = vshrq_n_u32(Rv1, 16);
            uint32x4_t Rv2_r = vshrq_n_u32(Rv2, 16);

            // Blend R and R' based on Cv.
            Rv1 = vbslq_u32(Cv1, Rv1_r, Rv1);
            Rv2 = vbslq_u32(Cv2, Rv2_r, Rv2);

            // R -> R' update
            //   q = (uint32_t) (((uint64_t)x * rcp_freq) >> rcp_shift);
            //   R' = R + sym->bias + q * sym->cmpl_freq;

            // Mix SIMD (mul) & scalar (shift). 365MB/s

            // We do 32 x 32 mul to get 64-bit, but then extract this
            // a 64-bit quantity and shift as scalar, before
            // recreating the 32x4 result.  Despite SIMD-scalar-SIMD reg
            // it's slightly quicker.

            uint64x2_t qvl1 = vmull_u32(vget_low_u32(Rv1), vget_low_u32(RFv1));
            uint64x2_t qvh1 = vmull_high_u32(Rv1, RFv1);

            uint64x2_t qvl2 = vmull_u32(vget_low_u32(Rv2), vget_low_u32(RFv2));
            uint64x2_t qvh2 = vmull_high_u32(Rv2, RFv2);

            uint32x2_t qv1a =
                 vcreate_u32(vgetq_lane_u64(qvl1, 1) >> s2->rcp_shift << 32 |
                             vgetq_lane_u64(qvl1, 0) >> s3->rcp_shift);
            uint32x2_t qv1b =
                 vcreate_u32(vgetq_lane_u64(qvh1, 1) >> s0->rcp_shift << 32 |
                             vgetq_lane_u64(qvh1, 0) >> s1->rcp_shift);

            uint32x2_t qv2a =
                 vcreate_u32(vgetq_lane_u64(qvl2, 1) >> s6->rcp_shift << 32 |
                             vgetq_lane_u64(qvl2, 0) >> s7->rcp_shift);
            uint32x2_t qv2b =
                 vcreate_u32(vgetq_lane_u64(qvh2, 1) >> s4->rcp_shift << 32 |
                             vgetq_lane_u64(qvh2, 0) >> s5->rcp_shift);

            uint32x4_t qv1 = vcombine_u32(qv1a, qv1b);
            uint32x4_t qv2 = vcombine_u32(qv2a, qv2b);
                
            FSv1 = vandq_u32(FSv1, vdupq_n_u32(0xffff)); // cmpl_freq
            FSv2 = vandq_u32(FSv2, vdupq_n_u32(0xffff));
        
            qv1 = vmlaq_u32(Biasv1, qv1, FSv1);
            qv2 = vmlaq_u32(Biasv2, qv2, FSv2);

            Rv1 = vaddq_u32(Rv1, qv1);
            Rv2 = vaddq_u32(Rv2, qv2);

            vst1q_u32(&R[z-3], Rv1);
            vst1q_u32(&R[z-7], Rv2);
        }
        if (z < -1) abort();
    }
    for (z = NX-1; z >= 0; z--)
      RansEncFlush(&R[z], &ptr);

 empty:
    // Finalise block size and return it
    *out_size = (out_end - ptr) + tab_size;

    memmove(out + tab_size, ptr, out_end-ptr);

    return out;
}

#define _ 99
static uint8x8_t idx[16] = {
    { _,_,_,_,_,_,_,_ }, // 0000
    { _,_,_,_,_,_,0,1 }, // 0001
    { _,_,_,_,0,1,_,_ }, // 0010
    { _,_,_,_,0,1,2,3 }, // 0011

    { _,_,0,1,_,_,_,_ }, // 0100
    { _,_,0,1,_,_,2,3 }, // 0101
    { _,_,0,1,2,3,_,_ }, // 0110
    { _,_,0,1,2,3,4,5 }, // 0111

    { 0,1,_,_,_,_,_,_ }, // 1000
    { 0,1,_,_,_,_,2,3 }, // 1001
    { 0,1,_,_,2,3,_,_ }, // 1010
    { 0,1,_,_,2,3,4,5 }, // 1011

    { 0,1,2,3,_,_,_,_ }, // 1100
    { 0,1,2,3,_,_,4,5 }, // 1101
    { 0,1,2,3,4,5,_,_ }, // 1110
    { 0,1,2,3,4,5,6,7 }, // 1111
};

// norm2 with norm1 in top 4 bits already consumed
static uint8x8_t idx2[256] = {
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 0, 1, },
    { _, _, _, _, 0, 1, _, _, },
    { _, _, _, _, 0, 1, 2, 3, },
    { _, _, 0, 1, _, _, _, _, },
    { _, _, 0, 1, _, _, 2, 3, },
    { _, _, 0, 1, 2, 3, _, _, },
    { _, _, 0, 1, 2, 3, 4, 5, },
    { 0, 1, _, _, _, _, _, _, },
    { 0, 1, _, _, _, _, 2, 3, },
    { 0, 1, _, _, 2, 3, _, _, },
    { 0, 1, _, _, 2, 3, 4, 5, },
    { 0, 1, 2, 3, _, _, _, _, },
    { 0, 1, 2, 3, _, _, 4, 5, },
    { 0, 1, 2, 3, 4, 5, _, _, },
    { 0, 1, 2, 3, 4, 5, 6, 7, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 2, 3, },
    { _, _, _, _, 2, 3, _, _, },
    { _, _, _, _, 2, 3, 4, 5, },
    { _, _, 2, 3, _, _, _, _, },
    { _, _, 2, 3, _, _, 4, 5, },
    { _, _, 2, 3, 4, 5, _, _, },
    { _, _, 2, 3, 4, 5, 6, 7, },
    { 2, 3, _, _, _, _, _, _, },
    { 2, 3, _, _, _, _, 4, 5, },
    { 2, 3, _, _, 4, 5, _, _, },
    { 2, 3, _, _, 4, 5, 6, 7, },
    { 2, 3, 4, 5, _, _, _, _, },
    { 2, 3, 4, 5, _, _, 6, 7, },
    { 2, 3, 4, 5, 6, 7, _, _, },
    { 2, 3, 4, 5, 6, 7, 8, 9, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 2, 3, },
    { _, _, _, _, 2, 3, _, _, },
    { _, _, _, _, 2, 3, 4, 5, },
    { _, _, 2, 3, _, _, _, _, },
    { _, _, 2, 3, _, _, 4, 5, },
    { _, _, 2, 3, 4, 5, _, _, },
    { _, _, 2, 3, 4, 5, 6, 7, },
    { 2, 3, _, _, _, _, _, _, },
    { 2, 3, _, _, _, _, 4, 5, },
    { 2, 3, _, _, 4, 5, _, _, },
    { 2, 3, _, _, 4, 5, 6, 7, },
    { 2, 3, 4, 5, _, _, _, _, },
    { 2, 3, 4, 5, _, _, 6, 7, },
    { 2, 3, 4, 5, 6, 7, _, _, },
    { 2, 3, 4, 5, 6, 7, 8, 9, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 4, 5, },
    { _, _, _, _, 4, 5, _, _, },
    { _, _, _, _, 4, 5, 6, 7, },
    { _, _, 4, 5, _, _, _, _, },
    { _, _, 4, 5, _, _, 6, 7, },
    { _, _, 4, 5, 6, 7, _, _, },
    { _, _, 4, 5, 6, 7, 8, 9, },
    { 4, 5, _, _, _, _, _, _, },
    { 4, 5, _, _, _, _, 6, 7, },
    { 4, 5, _, _, 6, 7, _, _, },
    { 4, 5, _, _, 6, 7, 8, 9, },
    { 4, 5, 6, 7, _, _, _, _, },
    { 4, 5, 6, 7, _, _, 8, 9, },
    { 4, 5, 6, 7, 8, 9, _, _, },
    { 4, 5, 6, 7, 8, 9,10,11, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 2, 3, },
    { _, _, _, _, 2, 3, _, _, },
    { _, _, _, _, 2, 3, 4, 5, },
    { _, _, 2, 3, _, _, _, _, },
    { _, _, 2, 3, _, _, 4, 5, },
    { _, _, 2, 3, 4, 5, _, _, },
    { _, _, 2, 3, 4, 5, 6, 7, },
    { 2, 3, _, _, _, _, _, _, },
    { 2, 3, _, _, _, _, 4, 5, },
    { 2, 3, _, _, 4, 5, _, _, },
    { 2, 3, _, _, 4, 5, 6, 7, },
    { 2, 3, 4, 5, _, _, _, _, },
    { 2, 3, 4, 5, _, _, 6, 7, },
    { 2, 3, 4, 5, 6, 7, _, _, },
    { 2, 3, 4, 5, 6, 7, 8, 9, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 4, 5, },
    { _, _, _, _, 4, 5, _, _, },
    { _, _, _, _, 4, 5, 6, 7, },
    { _, _, 4, 5, _, _, _, _, },
    { _, _, 4, 5, _, _, 6, 7, },
    { _, _, 4, 5, 6, 7, _, _, },
    { _, _, 4, 5, 6, 7, 8, 9, },
    { 4, 5, _, _, _, _, _, _, },
    { 4, 5, _, _, _, _, 6, 7, },
    { 4, 5, _, _, 6, 7, _, _, },
    { 4, 5, _, _, 6, 7, 8, 9, },
    { 4, 5, 6, 7, _, _, _, _, },
    { 4, 5, 6, 7, _, _, 8, 9, },
    { 4, 5, 6, 7, 8, 9, _, _, },
    { 4, 5, 6, 7, 8, 9,10,11, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 4, 5, },
    { _, _, _, _, 4, 5, _, _, },
    { _, _, _, _, 4, 5, 6, 7, },
    { _, _, 4, 5, _, _, _, _, },
    { _, _, 4, 5, _, _, 6, 7, },
    { _, _, 4, 5, 6, 7, _, _, },
    { _, _, 4, 5, 6, 7, 8, 9, },
    { 4, 5, _, _, _, _, _, _, },
    { 4, 5, _, _, _, _, 6, 7, },
    { 4, 5, _, _, 6, 7, _, _, },
    { 4, 5, _, _, 6, 7, 8, 9, },
    { 4, 5, 6, 7, _, _, _, _, },
    { 4, 5, 6, 7, _, _, 8, 9, },
    { 4, 5, 6, 7, 8, 9, _, _, },
    { 4, 5, 6, 7, 8, 9,10,11, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 6, 7, },
    { _, _, _, _, 6, 7, _, _, },
    { _, _, _, _, 6, 7, 8, 9, },
    { _, _, 6, 7, _, _, _, _, },
    { _, _, 6, 7, _, _, 8, 9, },
    { _, _, 6, 7, 8, 9, _, _, },
    { _, _, 6, 7, 8, 9,10,11, },
    { 6, 7, _, _, _, _, _, _, },
    { 6, 7, _, _, _, _, 8, 9, },
    { 6, 7, _, _, 8, 9, _, _, },
    { 6, 7, _, _, 8, 9,10,11, },
    { 6, 7, 8, 9, _, _, _, _, },
    { 6, 7, 8, 9, _, _,10,11, },
    { 6, 7, 8, 9,10,11, _, _, },
    { 6, 7, 8, 9,10,11,12,13, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 2, 3, },
    { _, _, _, _, 2, 3, _, _, },
    { _, _, _, _, 2, 3, 4, 5, },
    { _, _, 2, 3, _, _, _, _, },
    { _, _, 2, 3, _, _, 4, 5, },
    { _, _, 2, 3, 4, 5, _, _, },
    { _, _, 2, 3, 4, 5, 6, 7, },
    { 2, 3, _, _, _, _, _, _, },
    { 2, 3, _, _, _, _, 4, 5, },
    { 2, 3, _, _, 4, 5, _, _, },
    { 2, 3, _, _, 4, 5, 6, 7, },
    { 2, 3, 4, 5, _, _, _, _, },
    { 2, 3, 4, 5, _, _, 6, 7, },
    { 2, 3, 4, 5, 6, 7, _, _, },
    { 2, 3, 4, 5, 6, 7, 8, 9, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 4, 5, },
    { _, _, _, _, 4, 5, _, _, },
    { _, _, _, _, 4, 5, 6, 7, },
    { _, _, 4, 5, _, _, _, _, },
    { _, _, 4, 5, _, _, 6, 7, },
    { _, _, 4, 5, 6, 7, _, _, },
    { _, _, 4, 5, 6, 7, 8, 9, },
    { 4, 5, _, _, _, _, _, _, },
    { 4, 5, _, _, _, _, 6, 7, },
    { 4, 5, _, _, 6, 7, _, _, },
    { 4, 5, _, _, 6, 7, 8, 9, },
    { 4, 5, 6, 7, _, _, _, _, },
    { 4, 5, 6, 7, _, _, 8, 9, },
    { 4, 5, 6, 7, 8, 9, _, _, },
    { 4, 5, 6, 7, 8, 9,10,11, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 4, 5, },
    { _, _, _, _, 4, 5, _, _, },
    { _, _, _, _, 4, 5, 6, 7, },
    { _, _, 4, 5, _, _, _, _, },
    { _, _, 4, 5, _, _, 6, 7, },
    { _, _, 4, 5, 6, 7, _, _, },
    { _, _, 4, 5, 6, 7, 8, 9, },
    { 4, 5, _, _, _, _, _, _, },
    { 4, 5, _, _, _, _, 6, 7, },
    { 4, 5, _, _, 6, 7, _, _, },
    { 4, 5, _, _, 6, 7, 8, 9, },
    { 4, 5, 6, 7, _, _, _, _, },
    { 4, 5, 6, 7, _, _, 8, 9, },
    { 4, 5, 6, 7, 8, 9, _, _, },
    { 4, 5, 6, 7, 8, 9,10,11, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 6, 7, },
    { _, _, _, _, 6, 7, _, _, },
    { _, _, _, _, 6, 7, 8, 9, },
    { _, _, 6, 7, _, _, _, _, },
    { _, _, 6, 7, _, _, 8, 9, },
    { _, _, 6, 7, 8, 9, _, _, },
    { _, _, 6, 7, 8, 9,10,11, },
    { 6, 7, _, _, _, _, _, _, },
    { 6, 7, _, _, _, _, 8, 9, },
    { 6, 7, _, _, 8, 9, _, _, },
    { 6, 7, _, _, 8, 9,10,11, },
    { 6, 7, 8, 9, _, _, _, _, },
    { 6, 7, 8, 9, _, _,10,11, },
    { 6, 7, 8, 9,10,11, _, _, },
    { 6, 7, 8, 9,10,11,12,13, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 4, 5, },
    { _, _, _, _, 4, 5, _, _, },
    { _, _, _, _, 4, 5, 6, 7, },
    { _, _, 4, 5, _, _, _, _, },
    { _, _, 4, 5, _, _, 6, 7, },
    { _, _, 4, 5, 6, 7, _, _, },
    { _, _, 4, 5, 6, 7, 8, 9, },
    { 4, 5, _, _, _, _, _, _, },
    { 4, 5, _, _, _, _, 6, 7, },
    { 4, 5, _, _, 6, 7, _, _, },
    { 4, 5, _, _, 6, 7, 8, 9, },
    { 4, 5, 6, 7, _, _, _, _, },
    { 4, 5, 6, 7, _, _, 8, 9, },
    { 4, 5, 6, 7, 8, 9, _, _, },
    { 4, 5, 6, 7, 8, 9,10,11, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 6, 7, },
    { _, _, _, _, 6, 7, _, _, },
    { _, _, _, _, 6, 7, 8, 9, },
    { _, _, 6, 7, _, _, _, _, },
    { _, _, 6, 7, _, _, 8, 9, },
    { _, _, 6, 7, 8, 9, _, _, },
    { _, _, 6, 7, 8, 9,10,11, },
    { 6, 7, _, _, _, _, _, _, },
    { 6, 7, _, _, _, _, 8, 9, },
    { 6, 7, _, _, 8, 9, _, _, },
    { 6, 7, _, _, 8, 9,10,11, },
    { 6, 7, 8, 9, _, _, _, _, },
    { 6, 7, 8, 9, _, _,10,11, },
    { 6, 7, 8, 9,10,11, _, _, },
    { 6, 7, 8, 9,10,11,12,13, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 6, 7, },
    { _, _, _, _, 6, 7, _, _, },
    { _, _, _, _, 6, 7, 8, 9, },
    { _, _, 6, 7, _, _, _, _, },
    { _, _, 6, 7, _, _, 8, 9, },
    { _, _, 6, 7, 8, 9, _, _, },
    { _, _, 6, 7, 8, 9,10,11, },
    { 6, 7, _, _, _, _, _, _, },
    { 6, 7, _, _, _, _, 8, 9, },
    { 6, 7, _, _, 8, 9, _, _, },
    { 6, 7, _, _, 8, 9,10,11, },
    { 6, 7, 8, 9, _, _, _, _, },
    { 6, 7, 8, 9, _, _,10,11, },
    { 6, 7, 8, 9,10,11, _, _, },
    { 6, 7, 8, 9,10,11,12,13, },
    { _, _, _, _, _, _, _, _, },
    { _, _, _, _, _, _, 8, 9, },
    { _, _, _, _, 8, 9, _, _, },
    { _, _, _, _, 8, 9,10,11, },
    { _, _, 8, 9, _, _, _, _, },
    { _, _, 8, 9, _, _,10,11, },
    { _, _, 8, 9,10,11, _, _, },
    { _, _, 8, 9,10,11,12,13, },
    { 8, 9, _, _, _, _, _, _, },
    { 8, 9, _, _, _, _,10,11, },
    { 8, 9, _, _,10,11, _, _, },
    { 8, 9, _, _,10,11,12,13, },
    { 8, 9,10,11, _, _, _, _, },
    { 8, 9,10,11, _, _,12,13, },
    { 8, 9,10,11,12,13, _, _, },
    { 8, 9,10,11,12,13,14,15, },
};

// SIMD: 650MB/s
unsigned char *rans_uncompress_O0_32x16_neon(unsigned char *in,
                                             unsigned int in_size,
                                             unsigned char *out,
                                             unsigned int out_sz) {
    if (in_size < 16) // 4-states at least
        return NULL;

    if (out_sz >= INT_MAX)
        return NULL; // protect against some overflow cases

    /* Load in the static tables */
    unsigned char *cp = in, *out_free = NULL;
    unsigned char *cp_end = in + in_size;
    int i;
    uint32_t s3[TOTFREQ]; // For TF_SHIFT <= 12

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
    RansState R[NX];
    for (z = 0; z < NX; z++) {
        RansDecInit(&R[z], &cp);
        if (R[z] < RANS_BYTE_L)
            goto err;
    }

    int out_end = (out_sz&~(NX-1));
    const uint32_t mask = (1u << TF_SHIFT)-1;
    uint32x4_t maskv = vdupq_n_u32((1u << TF_SHIFT)-1);

    // assume NX is divisible by 4
    assert(NX%4==0);

    uint32x4_t Rv1 = vld1q_u32(&R[0]);
    uint32x4_t Rv2 = vld1q_u32(&R[4]);
    uint32x4_t Rv3 = vld1q_u32(&R[8]);
    uint32x4_t Rv4 = vld1q_u32(&R[12]);
    uint32x4_t Rv5 = vld1q_u32(&R[16]);
    uint32x4_t Rv6 = vld1q_u32(&R[20]);
    uint32x4_t Rv7 = vld1q_u32(&R[24]);
    uint32x4_t Rv8 = vld1q_u32(&R[28]);

    // Note this has a considerable amount of manual instruction reordering
    // to avoid latency.  We have 8 lanes of 4 rans states, but process 4
    // lanes at a time with the two sets of 4-lane steps interleaved to
    // ensure best use of processor pipeline units and latency removal.
    //
    // Unfortunately this is still poor with Clang-10.  Gcc more or less
    // honours this order and on my test set operates at ~675MB/s decode.
    // Clang without the manual reordering was at 440MB/s and with it at
    // 500MB/s.  Clang does a lot of reordering of this code, removing some
    // of the manual tuning benefits.  Short of dropping to assembly, for now
    // I would recommend using gcc to compile this file.
    uint16_t *sp = (uint16_t *)cp;
    uint8_t overflow[64+64] = {0};
    for (i=0; i < out_end; i+=NX) {
        // Decode freq, bias and symbol from s3 lookups
        uint32x4_t Sv1, Sv2, Sv3, Sv4, Sv5, Sv6, Sv7, Sv8;
        uint32x4_t Fv1, Fv2, Fv3, Fv4, Fv5, Fv6, Fv7, Fv8;
        uint32x4_t Bv1, Bv2, Bv3, Bv4, Bv5, Bv6, Bv7, Bv8;

        // Note we could check __ARM_FEATURE_MVE & 1 and use
        // vcreateq_u32 here, but I don't have a system to test with
        // so cannot validate if the code works.
        uint32x2_t s1a, s1b, s2a, s2b, s3a, s3b, s4a, s4b;
        s1a = vcreate_u32((uint64_t)(s3[R[ 1]&mask])<<32 | (s3[R[ 0]&mask]));
        s2a = vcreate_u32((uint64_t)(s3[R[ 5]&mask])<<32 | (s3[R[ 4]&mask]));
        s1b = vcreate_u32((uint64_t)(s3[R[ 3]&mask])<<32 | (s3[R[ 2]&mask]));
        s2b = vcreate_u32((uint64_t)(s3[R[ 7]&mask])<<32 | (s3[R[ 6]&mask]));
        s3a = vcreate_u32((uint64_t)(s3[R[ 9]&mask])<<32 | (s3[R[ 8]&mask]));
        s3b = vcreate_u32((uint64_t)(s3[R[11]&mask])<<32 | (s3[R[10]&mask]));
        s4a = vcreate_u32((uint64_t)(s3[R[13]&mask])<<32 | (s3[R[12]&mask]));
        s4b = vcreate_u32((uint64_t)(s3[R[15]&mask])<<32 | (s3[R[14]&mask]));

        Sv1 = vcombine_u32(s1a, s1b);
        Sv2 = vcombine_u32(s2a, s2b);
        Sv3 = vcombine_u32(s3a, s3b);
        Sv4 = vcombine_u32(s4a, s4b);

        Fv1 = vshrq_n_u32(Sv1, TF_SHIFT+8); // Freq = S >> TF_SHIFT+8
        Fv2 = vshrq_n_u32(Sv2, TF_SHIFT+8);
        Fv3 = vshrq_n_u32(Sv3, TF_SHIFT+8);
        Fv4 = vshrq_n_u32(Sv4, TF_SHIFT+8);

        uint32x2_t s5a, s5b, s6a, s6b, s7a, s7b, s8a, s8b;
        s5a = vcreate_u32((uint64_t)(s3[R[17]&mask])<<32 | (s3[R[16]&mask]));
        s5b = vcreate_u32((uint64_t)(s3[R[19]&mask])<<32 | (s3[R[18]&mask]));
        s6a = vcreate_u32((uint64_t)(s3[R[21]&mask])<<32 | (s3[R[20]&mask]));
        s6b = vcreate_u32((uint64_t)(s3[R[23]&mask])<<32 | (s3[R[22]&mask]));
        s7a = vcreate_u32((uint64_t)(s3[R[25]&mask])<<32 | (s3[R[24]&mask]));
        s7b = vcreate_u32((uint64_t)(s3[R[27]&mask])<<32 | (s3[R[26]&mask]));
        s8a = vcreate_u32((uint64_t)(s3[R[29]&mask])<<32 | (s3[R[28]&mask]));
        s8b = vcreate_u32((uint64_t)(s3[R[31]&mask])<<32 | (s3[R[30]&mask]));

        Bv1 = vshrq_n_u32(Sv1, 8);          // Bias = (S >> 8)
        Bv2 = vshrq_n_u32(Sv2, 8);
        Bv3 = vshrq_n_u32(Sv3, 8);
        Bv4 = vshrq_n_u32(Sv4, 8);

        // R[0] = (freq * (R[0] >> TF_SHIFT) + bias;
        Rv1 = vshrq_n_u32(Rv1, TF_SHIFT);   // R >> TF_SHIFT
        Rv2 = vshrq_n_u32(Rv2, TF_SHIFT);
        Rv3 = vshrq_n_u32(Rv3, TF_SHIFT);
        Rv4 = vshrq_n_u32(Rv4, TF_SHIFT);

        Sv5 = vcombine_u32(s5a, s5b);
        Sv6 = vcombine_u32(s6a, s6b);
        Sv7 = vcombine_u32(s7a, s7b);
        Sv8 = vcombine_u32(s8a, s8b);

        Bv1 = vandq_u32(Bv1, maskv);        //      & mask
        Bv2 = vandq_u32(Bv2, maskv);
        Bv3 = vandq_u32(Bv3, maskv);
        Bv4 = vandq_u32(Bv4, maskv);

        Fv5 = vshrq_n_u32(Sv5, TF_SHIFT+8);
        Fv6 = vshrq_n_u32(Sv6, TF_SHIFT+8);
        Fv7 = vshrq_n_u32(Sv7, TF_SHIFT+8);
        Fv8 = vshrq_n_u32(Sv8, TF_SHIFT+8);

        // A mix of mul+add and mla instructions seems to win.
        //Rv1 = vmulq_u32(Fv1, Rv1); Rv1 = vaddq_u32(Rv1, Bv1);
        Rv1 = vmlaq_u32(Bv1, Fv1, Rv1);     // R = R*Freq + Bias
        Rv2 = vmulq_u32(Fv2, Rv2); Rv2 = vaddq_u32(Rv2, Bv2);
        //Rv2 = vmlaq_u32(Bv2, Fv2, Rv2);
        //Rv3 = vmulq_u32(Fv3, Rv3); Rv3 = vaddq_u32(Rv3, Bv3);
        Rv3 = vmlaq_u32(Bv3, Fv3, Rv3);
        Rv4 = vmulq_u32(Fv4, Rv4); Rv4 = vaddq_u32(Rv4, Bv4);
        //Rv4 = vmlaq_u32(Bv4, Fv4, Rv4);

        Bv5 = vshrq_n_u32(Sv5, 8);
        Bv6 = vshrq_n_u32(Sv6, 8);
        Bv7 = vshrq_n_u32(Sv7, 8);
        Bv8 = vshrq_n_u32(Sv8, 8);

        // Renorm
        uint32x4_t Rlt1 = vcltq_u32(Rv1, vdupq_n_u32(RANS_BYTE_L)); // R<L
        uint32x4_t Rlt2 = vcltq_u32(Rv2, vdupq_n_u32(RANS_BYTE_L));
        uint32x4_t Rlt3 = vcltq_u32(Rv3, vdupq_n_u32(RANS_BYTE_L));
        uint32x4_t Rlt4 = vcltq_u32(Rv4, vdupq_n_u32(RANS_BYTE_L));

        Rv5 = vshrq_n_u32(Rv5, TF_SHIFT);
        Rv6 = vshrq_n_u32(Rv6, TF_SHIFT);
        Rv7 = vshrq_n_u32(Rv7, TF_SHIFT);
        Rv8 = vshrq_n_u32(Rv8, TF_SHIFT);

        static int nbits[16] = { 0,1,1,2, 1,2,2,3, 1,2,2,3, 2,3,3,4 };

        // Combine with Rlt & {8,4,2,1} to get single bits.
        // Sum lanes to git a bit-field indicating which R < L
        uint32x4_t bit = {8,4,2,1};
        uint32_t imask1 = vaddvq_u32(vandq_u32(Rlt1, bit));
        uint32_t imask2 = vaddvq_u32(vandq_u32(Rlt2, bit));
        uint32_t imask3 = vaddvq_u32(vandq_u32(Rlt3, bit));
        uint32_t imask4 = vaddvq_u32(vandq_u32(Rlt4, bit));

        // Protect against running off the end of in buffer.
        // We copy it to a worst-case local buffer when near the end.
        if ((uint8_t *)sp+64 > cp_end) {
            memmove(overflow, sp, cp_end - (uint8_t *)sp);
            sp = (uint16_t *)overflow;
            cp_end = overflow + sizeof(overflow);
        }

        uint16x8_t norm12 =  vld1q_u16(sp);
        sp += nbits[imask1] + nbits[imask2];
        uint16x8_t norm34 =  vld1q_u16(sp);
        sp += nbits[imask3] + nbits[imask4];

        Bv5 = vandq_u32(Bv5, maskv);
        Bv6 = vandq_u32(Bv6, maskv);
        Bv7 = vandq_u32(Bv7, maskv);
        Bv8 = vandq_u32(Bv8, maskv);

        // Shuffle norm to the corresponding R lanes, via imask
        //Rv5 = vmulq_u32(Fv5, Rv5); Rv5 = vaddq_u32(Rv5, Bv5);
        Rv5 = vmlaq_u32(Bv5, Fv5, Rv5);
        Rv6 = vmulq_u32(Fv6, Rv6); Rv6 = vaddq_u32(Rv6, Bv6);
        //Rv6 = vmlaq_u32(Bv6, Fv6, Rv6);
        //Rv7 = vmulq_u32(Fv7, Rv7); Rv7 = vaddq_u32(Rv7, Bv7);
        Rv7 = vmlaq_u32(Bv7, Fv7, Rv7);
        Rv8 = vmulq_u32(Fv8, Rv8); Rv8 = vaddq_u32(Rv8, Bv8);
        //Rv8 = vmlaq_u32(Bv8, Fv8, Rv8);
        
        uint32_t imask12 = (imask1<<4)|imask2;
        uint32_t imask34 = (imask3<<4)|imask4;

        // #define for brevity and formatting
#define cast_u16_u8 vreinterpret_u16_u8
#define cast_u8_u16 vreinterpretq_u8_u16
        uint16x4_t norm1, norm2, norm3, norm4, norm5, norm6, norm7, norm8;
        norm1 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm12),idx [imask1]));
        norm2 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm12),idx2[imask12]));
        norm3 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm34),idx [imask3]));
        norm4 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm34),idx2[imask34]));

        uint32x4_t Rlt5 = vcltq_u32(Rv5, vdupq_n_u32(RANS_BYTE_L));
        uint32x4_t Rlt6 = vcltq_u32(Rv6, vdupq_n_u32(RANS_BYTE_L));
        uint32x4_t Rlt7 = vcltq_u32(Rv7, vdupq_n_u32(RANS_BYTE_L));
        uint32x4_t Rlt8 = vcltq_u32(Rv8, vdupq_n_u32(RANS_BYTE_L));

        // Add norm to R<<16 (Rsl) and blend back in with R
        uint32x4_t Rsl1 = vshlq_n_u32(Rv1, 16); // Rsl = R << 16
        uint32x4_t Rsl2 = vshlq_n_u32(Rv2, 16);
        uint32x4_t Rsl3 = vshlq_n_u32(Rv3, 16);
        uint32x4_t Rsl4 = vshlq_n_u32(Rv4, 16);

        uint16x8_t norm56 =  vld1q_u16(sp);
        uint32_t imask5 = vaddvq_u32(vandq_u32(Rlt5, bit));
        uint32_t imask6 = vaddvq_u32(vandq_u32(Rlt6, bit));
        uint32_t imask7 = vaddvq_u32(vandq_u32(Rlt7, bit));
        uint32_t imask8 = vaddvq_u32(vandq_u32(Rlt8, bit));

        sp += nbits[imask5] + nbits[imask6];
        uint16x8_t norm78 =  vld1q_u16(sp);
        sp += nbits[imask7] + nbits[imask8];

        Rsl1 = vaddw_u16(Rsl1, norm1);          // Rsl += norm
        Rsl2 = vaddw_u16(Rsl2, norm2);
        Rsl3 = vaddw_u16(Rsl3, norm3);
        Rsl4 = vaddw_u16(Rsl4, norm4);

        uint32_t imask56 = (imask5<<4)|imask6;
        uint32_t imask78 = (imask7<<4)|imask8;

        Rv1 = vbslq_u32(Rlt1, Rsl1, Rv1);       // R = R<L ? Rsl : R
        Rv2 = vbslq_u32(Rlt2, Rsl2, Rv2);
        Rv3 = vbslq_u32(Rlt3, Rsl3, Rv3);
        Rv4 = vbslq_u32(Rlt4, Rsl4, Rv4);

        norm5 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm56),idx [imask5]));
        norm6 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm56),idx2[imask56]));
        norm7 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm78),idx [imask7]));
        norm8 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm78),idx2[imask78]));

        uint32x4_t Rsl5 = vshlq_n_u32(Rv5, 16);
        uint32x4_t Rsl6 = vshlq_n_u32(Rv6, 16);
        uint32x4_t Rsl7 = vshlq_n_u32(Rv7, 16);
        uint32x4_t Rsl8 = vshlq_n_u32(Rv8, 16);

        Rsl5 = vaddw_u16(Rsl5, norm5);
        Rsl6 = vaddw_u16(Rsl6, norm6);
        Rsl7 = vaddw_u16(Rsl7, norm7);
        Rsl8 = vaddw_u16(Rsl8, norm8);

        // Can we avoid storing this in memory?
        // It's tricky due to the need to do s3[R[z] & mask] and no
        // vector gathers.  So we need to have something in memory.
        vst1q_u32(&R[ 0], Rv1); // store R for Sv creation
        vst1q_u32(&R[ 4], Rv2);
        vst1q_u32(&R[ 8], Rv3);
        vst1q_u32(&R[12], Rv4);

        Rv5 = vbslq_u32(Rlt5, Rsl5, Rv5);
        Rv6 = vbslq_u32(Rlt6, Rsl6, Rv6);
        Rv7 = vbslq_u32(Rlt7, Rsl7, Rv7);
        Rv8 = vbslq_u32(Rlt8, Rsl8, Rv8);

        vst1q_u32(&R[16], Rv5);
        vst1q_u32(&R[20], Rv6);
        vst1q_u32(&R[24], Rv7);
        vst1q_u32(&R[28], Rv8);

        // out[i+z+?] = S[?]
        uint16x4_t p16_1 = vmovn_u32(Sv1);
        uint16x4_t p16_2 = vmovn_u32(Sv2);
        uint16x4_t p16_3 = vmovn_u32(Sv3);
        uint16x4_t p16_4 = vmovn_u32(Sv4);
        uint16x4_t p16_5 = vmovn_u32(Sv5);
        uint16x4_t p16_6 = vmovn_u32(Sv6);
        uint16x4_t p16_7 = vmovn_u32(Sv7);
        uint16x4_t p16_8 = vmovn_u32(Sv8);
        uint8x8_t  p8_12  = vmovn_u16(vcombine_u16(p16_1,p16_2));
        uint8x8_t  p8_34  = vmovn_u16(vcombine_u16(p16_3,p16_4));
        uint8x8_t  p8_56  = vmovn_u16(vcombine_u16(p16_5,p16_6));
        uint8x8_t  p8_78  = vmovn_u16(vcombine_u16(p16_7,p16_8));
        uint8x16_t p8_a   = vcombine_u8(p8_12, p8_34);
        uint8x16_t p8_b   = vcombine_u8(p8_56, p8_78);
        vst1q_u8(&out[i   ], p8_a);
        vst1q_u8(&out[i+16], p8_b);
    }

    for (z = out_sz & (NX-1); z-- > 0; )
      out[out_end + z] = s3[R[z] & mask];

    //fprintf(stderr, "    0 Decoded %d bytes\n", (int)(cp-in)); //c-size

    return out;

 err:
    free(out_free);
    return NULL;
}

//-----------------------------------------------------------------------------

unsigned char *rans_compress_O1_32x16_neon(unsigned char *in,
                                           unsigned int in_size,
                                           unsigned char *out,
                                           unsigned int *out_size) {
    unsigned char *cp, *out_end, *out_free = NULL;
    unsigned int tab_size;
    uint32_t bound = rans_compress_bound_4x16(in_size,1)-20;
    int z;
    RansState ransN[NX];

    if (in_size < NX) // force O0 instead
        return NULL;

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
    int shift = encode_freq1(in, in_size, 32, syms, &cp); 
    if (shift < 0) {
        free(out_free);
        htscodecs_tls_free(syms);
        return NULL;
    }
    tab_size = cp - out;

    for (z = 0; z < NX; z++)
      RansEncInit(&ransN[z]);

    uint8_t* ptr = out_end;

    int iN[NX], isz4 = in_size/NX;
    for (z = 0; z < NX; z++)
        iN[z] = (z+1)*isz4-2;

    unsigned char lN[NX];
    for (z = 0; z < NX; z++)
        lN[z] = in[iN[z]+1];

    // Deal with the remainder
    z = NX-1;
    lN[z] = in[in_size-1];
    for (iN[z] = in_size-2; iN[z] > NX*isz4-2; iN[z]--) {
        unsigned char c = in[iN[z]];
        RansEncPutSymbol(&ransN[z], &ptr, &syms[c][lN[z]]);
        lN[z] = c;
    }

#if 0
    // Scalar code equivalent
    for (; iN[0] >= 0; ) {
        for (z = NX-1; z >= 0; z-=4) {
            unsigned char c0;
            unsigned char c1;
            unsigned char c2;
            unsigned char c3;

            RansEncSymbol *s0 = &syms[c0=in[iN[z-0]--]][lN[z-0]]; lN[z-0] = c0;
            RansEncSymbol *s1 = &syms[c1=in[iN[z-1]--]][lN[z-1]]; lN[z-1] = c1;
            RansEncSymbol *s2 = &syms[c2=in[iN[z-2]--]][lN[z-2]]; lN[z-2] = c2;
            RansEncSymbol *s3 = &syms[c3=in[iN[z-3]--]][lN[z-3]]; lN[z-3] = c3;

            RansEncPutSymbol(&ransN[z-0], &ptr, s0);
            RansEncPutSymbol(&ransN[z-1], &ptr, s1);
            RansEncPutSymbol(&ransN[z-2], &ptr, s2);
            RansEncPutSymbol(&ransN[z-3], &ptr, s3);
        }
    }
#else
    // SIMD code
    for (; iN[0] >= 0; ) {
        for (z = NX-1; z >= 0; z-=16) {
            unsigned char c;

            RansEncSymbol *s0 = &syms[c=in[iN[z- 0]--]][lN[z- 0]]; lN[z- 0]=c;
            RansEncSymbol *s1 = &syms[c=in[iN[z- 1]--]][lN[z- 1]]; lN[z- 1]=c;
            RansEncSymbol *s2 = &syms[c=in[iN[z- 2]--]][lN[z- 2]]; lN[z- 2]=c;
            RansEncSymbol *s3 = &syms[c=in[iN[z- 3]--]][lN[z- 3]]; lN[z- 3]=c;

            uint32x4_t Rv1 = vld1q_u32(&ransN[z-3]);
            uint32x4_t Rv2 = vld1q_u32(&ransN[z-7]);
            uint32x4_t Rv3 = vld1q_u32(&ransN[z-11]);
            uint32x4_t Rv4 = vld1q_u32(&ransN[z-15]);
        
            RansEncSymbol *s4 = &syms[c=in[iN[z- 4]--]][lN[z- 4]]; lN[z- 4]=c;
            RansEncSymbol *s5 = &syms[c=in[iN[z- 5]--]][lN[z- 5]]; lN[z- 5]=c;
            RansEncSymbol *s6 = &syms[c=in[iN[z- 6]--]][lN[z- 6]]; lN[z- 6]=c;
            RansEncSymbol *s7 = &syms[c=in[iN[z- 7]--]][lN[z- 7]]; lN[z- 7]=c;

            uint32x4_t A_1 = vld1q_u32((void *)s3);
            uint32x4_t B_1 = vld1q_u32((void *)s2);
            uint32x4_t C_1 = vld1q_u32((void *)s1);
            uint32x4_t D_1 = vld1q_u32((void *)s0);

            uint32x4_t A1_1 = vtrn1q_u32(A_1, B_1);
            uint32x4_t C1_1 = vtrn1q_u32(C_1, D_1);
            uint32x4_t A2_1 = vtrn2q_u32(A_1, B_1);
            uint32x4_t C2_1 = vtrn2q_u32(C_1, D_1);

            uint32x4_t Xmaxv1=u32_u64(vtrn1q_u64(u64_u32(A1_1),u64_u32(C1_1)));
            uint32x4_t Biasv1=u32_u64(vtrn2q_u64(u64_u32(A1_1),u64_u32(C1_1)));
            uint32x4_t RFv1  =u32_u64(vtrn1q_u64(u64_u32(A2_1),u64_u32(C2_1)));
            uint32x4_t FSv1  =u32_u64(vtrn2q_u64(u64_u32(A2_1),u64_u32(C2_1)));

            uint32x4_t A_2 = vld1q_u32((void *)s7);
            uint32x4_t B_2 = vld1q_u32((void *)s6);
            uint32x4_t C_2 = vld1q_u32((void *)s5);
            uint32x4_t D_2 = vld1q_u32((void *)s4);
            uint32x4_t A1_2 = vtrn1q_u32(A_2, B_2);
            uint32x4_t C1_2 = vtrn1q_u32(C_2, D_2);
            uint32x4_t A2_2 = vtrn2q_u32(A_2, B_2);
            uint32x4_t C2_2 = vtrn2q_u32(C_2, D_2);

            uint32x4_t Xmaxv2=u32_u64(vtrn1q_u64(u64_u32(A1_2),u64_u32(C1_2)));
            uint32x4_t Biasv2=u32_u64(vtrn2q_u64(u64_u32(A1_2),u64_u32(C1_2)));
            uint32x4_t RFv2  =u32_u64(vtrn1q_u64(u64_u32(A2_2),u64_u32(C2_2)));
            uint32x4_t FSv2  =u32_u64(vtrn2q_u64(u64_u32(A2_2),u64_u32(C2_2)));

            uint32x4_t Cv1 = vcgtq_u32(Rv1, Xmaxv1);
            uint32x4_t Cv2 = vcgtq_u32(Rv2, Xmaxv2);
            uint32x4_t bit = {8,4,2,1};
            uint32_t imask1 = vaddvq_u32(vandq_u32(Cv1, bit));
            uint32_t imask2 = vaddvq_u32(vandq_u32(Cv2, bit));

            RansEncSymbol *s8 = &syms[c=in[iN[z- 8]--]][lN[z- 8]]; lN[z- 8]=c;
            RansEncSymbol *s9 = &syms[c=in[iN[z- 9]--]][lN[z- 9]]; lN[z- 9]=c;
            RansEncSymbol *s10= &syms[c=in[iN[z-10]--]][lN[z-10]]; lN[z-10]=c;
            RansEncSymbol *s11= &syms[c=in[iN[z-11]--]][lN[z-11]]; lN[z-11]=c;

            RansEncSymbol *s12= &syms[c=in[iN[z-12]--]][lN[z-12]]; lN[z-12]=c;
            RansEncSymbol *s13= &syms[c=in[iN[z-13]--]][lN[z-13]]; lN[z-13]=c;
            RansEncSymbol *s14= &syms[c=in[iN[z-14]--]][lN[z-14]]; lN[z-14]=c;
            RansEncSymbol *s15= &syms[c=in[iN[z-15]--]][lN[z-15]]; lN[z-15]=c;

            uint32x4_t A_3 = vld1q_u32((void *)s11);
            uint32x4_t B_3 = vld1q_u32((void *)s10);
            uint32x4_t C_3 = vld1q_u32((void *)s9);
            uint32x4_t D_3 = vld1q_u32((void *)s8);


            uint32x4_t A1_3 = vtrn1q_u32(A_3, B_3);
            uint32x4_t C1_3 = vtrn1q_u32(C_3, D_3);
            uint32x4_t A2_3 = vtrn2q_u32(A_3, B_3);
            uint32x4_t C2_3 = vtrn2q_u32(C_3, D_3);

            uint32x4_t Xmaxv3=u32_u64(vtrn1q_u64(u64_u32(A1_3),u64_u32(C1_3)));
            uint32x4_t Biasv3=u32_u64(vtrn2q_u64(u64_u32(A1_3),u64_u32(C1_3)));
            uint32x4_t RFv3  =u32_u64(vtrn1q_u64(u64_u32(A2_3),u64_u32(C2_3)));
            uint32x4_t FSv3  =u32_u64(vtrn2q_u64(u64_u32(A2_3),u64_u32(C2_3)));

            uint32x4_t A_4 = vld1q_u32((void *)s15);
            uint32x4_t B_4 = vld1q_u32((void *)s14);
            uint32x4_t C_4 = vld1q_u32((void *)s13);
            uint32x4_t D_4 = vld1q_u32((void *)s12);

            uint32x4_t A1_4 = vtrn1q_u32(A_4, B_4);
            uint32x4_t C1_4 = vtrn1q_u32(C_4, D_4);
            uint32x4_t A2_4 = vtrn2q_u32(A_4, B_4);
            uint32x4_t C2_4 = vtrn2q_u32(C_4, D_4);

            uint32x4_t Xmaxv4=u32_u64(vtrn1q_u64(u64_u32(A1_4),u64_u32(C1_4)));
            uint32x4_t Biasv4=u32_u64(vtrn2q_u64(u64_u32(A1_4),u64_u32(C1_4)));
            uint32x4_t RFv4  =u32_u64(vtrn1q_u64(u64_u32(A2_4),u64_u32(C2_4)));
            uint32x4_t FSv4  =u32_u64(vtrn2q_u64(u64_u32(A2_4),u64_u32(C2_4)));

            uint32x4_t Cv3 = vcgtq_u32(Rv3, Xmaxv3);
            uint32x4_t Cv4 = vcgtq_u32(Rv4, Xmaxv4);
            uint32_t imask3 = vaddvq_u32(vandq_u32(Cv3, bit));
            uint32_t imask4 = vaddvq_u32(vandq_u32(Cv4, bit));

            // Select low 16-bits from Rv based on imask, using tbl
            uint8x8_t norm1, norm2, norm3, norm4;
            static int nbits[16] = { 0,2,2,4, 2,4,4,6, 2,4,4,6, 4,6,6,8 };
            norm1 = vqtbl1_u8(vreinterpretq_u8_u32(Rv1),vtab[imask1]);
            norm2 = vqtbl1_u8(vreinterpretq_u8_u32(Rv2),vtab[imask2]);

            vst1_u8(ptr-8, norm1);   ptr -= nbits[imask1];
            vst1_u8(ptr-8, norm2);   ptr -= nbits[imask2];

            norm3 = vqtbl1_u8(vreinterpretq_u8_u32(Rv3),vtab[imask3]);
            norm4 = vqtbl1_u8(vreinterpretq_u8_u32(Rv4),vtab[imask4]);

            vst1_u8(ptr-8, norm3);   ptr -= nbits[imask3];
            vst1_u8(ptr-8, norm4);   ptr -= nbits[imask4];

            // R' = R>>16
            uint32x4_t Rv1_r = vshrq_n_u32(Rv1, 16);
            uint32x4_t Rv2_r = vshrq_n_u32(Rv2, 16);
            uint32x4_t Rv3_r = vshrq_n_u32(Rv3, 16);
            uint32x4_t Rv4_r = vshrq_n_u32(Rv4, 16);

            // Blend R and R' based on Cv.
            Rv1 = vbslq_u32(Cv1, Rv1_r, Rv1);
            Rv2 = vbslq_u32(Cv2, Rv2_r, Rv2);
            Rv3 = vbslq_u32(Cv3, Rv3_r, Rv3);
            Rv4 = vbslq_u32(Cv4, Rv4_r, Rv4);

            uint64x2_t qvl1 = vmull_u32(vget_low_u32(Rv1), vget_low_u32(RFv1));
            uint64x2_t qvh1 = vmull_high_u32(Rv1, RFv1);
            uint64x2_t qvl2 = vmull_u32(vget_low_u32(Rv2), vget_low_u32(RFv2));
            uint64x2_t qvh2 = vmull_high_u32(Rv2, RFv2);

            int32x4_t RSv1 = vnegq_s32(vreinterpretq_s32_u32(
                                           vshrq_n_u32(FSv1, 16)));
            int32x4_t RSv2 = vnegq_s32(vreinterpretq_s32_u32(
                                           vshrq_n_u32(FSv2, 16)));

            uint64x2_t qvl3 = vmull_u32(vget_low_u32(Rv3), vget_low_u32(RFv3));
            uint64x2_t qvh3 = vmull_high_u32(Rv3, RFv3);
            uint64x2_t qvl4 = vmull_u32(vget_low_u32(Rv4), vget_low_u32(RFv4));
            uint64x2_t qvh4 = vmull_high_u32(Rv4, RFv4);

            int32x4_t RSv3 = vnegq_s32(vreinterpretq_s32_u32(
                                           vshrq_n_u32(FSv3, 16)));
            int32x4_t RSv4 = vnegq_s32(vreinterpretq_s32_u32(
                                           vshrq_n_u32(FSv4, 16)));

            qvl1 = vreinterpretq_u64_s64(
                       vshlq_s64(vreinterpretq_s64_u64(qvl1),
                                 vmovl_s32(vget_low_s32(RSv1))));
            qvh1 = vreinterpretq_u64_s64(
                       vshlq_s64(vreinterpretq_s64_u64(qvh1),
                                 vmovl_s32(vget_high_s32(RSv1))));

            qvl2 = vreinterpretq_u64_s64(
                       vshlq_s64(vreinterpretq_s64_u64(qvl2),
                                 vmovl_s32(vget_low_s32(RSv2))));
            qvh2 = vreinterpretq_u64_s64(
                       vshlq_s64(vreinterpretq_s64_u64(qvh2),
                                 vmovl_s32(vget_high_s32(RSv2))));

            uint32x4_t qv1 = vcombine_u32(vmovn_u64(qvl1),
                                          vmovn_u64(qvh1));
            uint32x4_t qv2 = vcombine_u32(vmovn_u64(qvl2),
                                          vmovn_u64(qvh2));

            qvl3 = vreinterpretq_u64_s64(
                       vshlq_s64(vreinterpretq_s64_u64(qvl3),
                                 vmovl_s32(vget_low_s32(RSv3))));
            qvh3 = vreinterpretq_u64_s64(
                       vshlq_s64(vreinterpretq_s64_u64(qvh3),
                                 vmovl_s32(vget_high_s32(RSv3))));

            qvl4 = vreinterpretq_u64_s64(
                       vshlq_s64(vreinterpretq_s64_u64(qvl4),
                                 vmovl_s32(vget_low_s32(RSv4))));
            qvh4 = vreinterpretq_u64_s64(
                       vshlq_s64(vreinterpretq_s64_u64(qvh4),
                                 vmovl_s32(vget_high_s32(RSv4))));

            uint32x4_t qv3 = vcombine_u32(vmovn_u64(qvl3),
                                          vmovn_u64(qvh3));
            uint32x4_t qv4 = vcombine_u32(vmovn_u64(qvl4),
                                          vmovn_u64(qvh4));

            FSv1 = vandq_u32(FSv1, vdupq_n_u32(0xffff)); // cmpl_freq
            FSv2 = vandq_u32(FSv2, vdupq_n_u32(0xffff));
            FSv3 = vandq_u32(FSv3, vdupq_n_u32(0xffff));
            FSv4 = vandq_u32(FSv4, vdupq_n_u32(0xffff));
        
            qv1 = vmlaq_u32(Biasv1, qv1, FSv1);
            qv2 = vmlaq_u32(Biasv2, qv2, FSv2);
            qv3 = vmlaq_u32(Biasv3, qv3, FSv3);
            qv4 = vmlaq_u32(Biasv4, qv4, FSv4);

            Rv1 = vaddq_u32(Rv1, qv1);
            Rv2 = vaddq_u32(Rv2, qv2);
            Rv3 = vaddq_u32(Rv3, qv3);
            Rv4 = vaddq_u32(Rv4, qv4);

            vst1q_u32(&ransN[z-3], Rv1);
            vst1q_u32(&ransN[z-7], Rv2);
            vst1q_u32(&ransN[z-11],Rv3);
            vst1q_u32(&ransN[z-15],Rv4);
        }
    }
#endif

    for (z = NX-1; z>=0; z--)
        RansEncPutSymbol(&ransN[z], &ptr, &syms[0][lN[z]]);

    for (z = NX-1; z>=0; z--)
        RansEncFlush(&ransN[z], &ptr);

    *out_size = (out_end - ptr) + tab_size;

    cp = out;
    memmove(out + tab_size, ptr, out_end-ptr);

    htscodecs_tls_free(syms);
    return out;
}

//#define MAGIC2 111
#define MAGIC2 179
//#define MAGIC2 0
typedef struct {
  union {
    struct {
      uint16_t f;
      uint16_t b;
    } s;
    uint32_t fb;
  } u;
} bf_t;

static inline void transpose_and_copy(uint8_t *out, int iN[32],
                                      uint8_t t[32][32]) {
    int z;
//  for (z = 0; z < NX; z++) {
//      int k;
//      for (k = 0; k < 32; k++)
//      out[iN[z]+k] = t[k][z];
//      iN[z] += 32;
//  }

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
}

unsigned char *rans_uncompress_O1_32x16_neon(unsigned char *in,
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
    int i, j = -999;
    unsigned int x;

    uint8_t *sfb_ = htscodecs_tls_alloc(256*(TOTFREQ_O1+MAGIC2)*sizeof(*sfb_));
    uint32_t s3[256][TOTFREQ_O1_FAST];

    if (!sfb_)
        return NULL;
    bf_t fb[256][256];
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
#if 0
    // Disable inline for now as this is ~10% slower under gcc.  Why?
    cp += decode_freq1(cp, c_freq_end, shift, NULL, s3, sfb, fb);
#else
    uint32_t F0[256] = {0};
    int fsz = decode_alphabet(cp, c_freq_end, F0);
    if (!fsz)
        goto err;
    cp += fsz;

    if (cp >= c_freq_end)
        goto err;

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

                if (shift == TF_SHIFT_O1_FAST) {
                    int y;
                    for (y = 0; y < F[j]; y++)
                        s3[i][y+x] = (((uint32_t)F[j])<<(shift+8)) |(y<<8) |j;
                }
                memset(&sfb[i][x], j, F[j]);
                fb[i][j].u.s.f = F[j];
                fb[i][j].u.s.b = x;

                x += F[j];
            }
        }
        if (x != (1<<shift))
            goto err;
    }
#endif

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
    int i4[NX];
    uint8_t l[NX] = {0};
    for (z = 0; z < NX; z++)
        i4[z] = z*isz4;

    // Around 15% faster to specialise for 10/12 than to have one
    // loop with shift as a variable.
    if (shift == TF_SHIFT_O1) {
        // TF_SHIFT_O1 = 12
        const uint32_t mask = ((1u << TF_SHIFT_O1)-1);
        uint32x4_t maskv = vdupq_n_u32((1u << TF_SHIFT_O1)-1);

        // FIXME: plus room for "safe" renorm.
        // Follow with 2nd copy doing scalar code instead?
        unsigned char tbuf[32][32] = {0};
        int tidx = 0;
        for (; i4[0] < isz4 && ptr+64 < ptr_end;) {
            for (z = 0; z < NX; z+=16) {
                uint32x4_t Rv1 = vld1q_u32(&R[z+0]);
                uint32x4_t Rv2 = vld1q_u32(&R[z+4]);
                uint32x4_t Rv3 = vld1q_u32(&R[z+8]);
                uint32x4_t Rv4 = vld1q_u32(&R[z+12]);

                uint32x4_t Sv1, Sv2, Sv3, Sv4;
                uint32x4_t Fv1, Fv2, Fv3, Fv4;
                uint32x4_t Bv1, Bv2, Bv3, Bv4;

                // Use sfb[256][256] + fb[256][256] instead of s3[256][4096].
                // It's lower memory and thus less cache contention and also
                // avoids the 12-bit issue of FREQ==0 vs 4096.
                // However we need an extra step to store m = R&mask.
                // Overall it's a big win for gcc (no difference to clang, but
                // that's lagging behind gcc on this code).
                //
                // TODO: consider this for AVX2?
                //
                // sfb layout is char (sym)
                // fb layout is 16-bit freq, 16-bit bias.
                uint8_t *c = &tbuf[tidx][z];
                uint32_t f32[16];
                int i;
                for (i = 0; i < 16; i++) {
                  c[i] = sfb[l[z+i]][R[z+i] & mask];
                  f32[i] = *(uint32_t *)&fb[l[z+i]][c[i]];
                }
                
                // vcreate faster than vld1q_u32(&f32[0])
                uint32x2_t s1a, s1b, s2a, s2b, s3a, s3b, s4a, s4b;
                s1a = vcreate_u32((uint64_t)(f32[ 3])<<32 |
                                  (uint64_t)(f32[ 2]));
                s1b = vcreate_u32((uint64_t)(f32[ 1])<<32 |
                                  (uint64_t)(f32[ 0]));
                s2a = vcreate_u32((uint64_t)(f32[ 7])<<32 |
                                  (uint64_t)(f32[ 6]));
                s2b = vcreate_u32((uint64_t)(f32[ 5])<<32 |
                                  (uint64_t)(f32[ 4]));
                s3a = vcreate_u32((uint64_t)(f32[11])<<32 |
                                  (uint64_t)(f32[10]));
                s3b = vcreate_u32((uint64_t)(f32[ 9])<<32 |
                                  (uint64_t)(f32[ 8]));
                s4a = vcreate_u32((uint64_t)(f32[15])<<32 |
                                  (uint64_t)(f32[14]));
                s4b = vcreate_u32((uint64_t)(f32[13])<<32 |
                                  (uint64_t)(f32[12]));

                memcpy(l+z, c, 16);

                // vcreate   = INS   = throughput 2, latency 2, pipe V
                // vcombine  = DUP+INS = thr    2+2, lat   2+2, pipe V+V
                // vandq     = AND   = throughput 2, latency 1, pipe V
                // vshrq     = USHR  = throughput 2, latency 1, pipe V
                Sv1 = vcombine_u32(s1b, s1a);
                Sv2 = vcombine_u32(s2b, s2a);
                Sv3 = vcombine_u32(s3b, s3a);
                Sv4 = vcombine_u32(s4b, s4a);

                Bv1 = vshrq_n_u32(Sv1, 16);
                Bv2 = vshrq_n_u32(Sv2, 16);

                Bv3 = vshrq_n_u32(Sv3, 16);
                Bv4 = vshrq_n_u32(Sv4, 16);

                Fv1 = vandq_u32(Sv1, vdupq_n_u32(0xffff));
                Fv2 = vandq_u32(Sv2, vdupq_n_u32(0xffff));
                Fv3 = vandq_u32(Sv3, vdupq_n_u32(0xffff));
                Fv4 = vandq_u32(Sv4, vdupq_n_u32(0xffff));

                Bv1 = vsubq_u32(vandq_u32(Rv1, maskv), Bv1);
                Bv2 = vsubq_u32(vandq_u32(Rv2, maskv), Bv2);
                Bv3 = vsubq_u32(vandq_u32(Rv3, maskv), Bv3);
                Bv4 = vsubq_u32(vandq_u32(Rv4, maskv), Bv4);

                Rv1 = vshrq_n_u32(Rv1, TF_SHIFT_O1);
                Rv2 = vshrq_n_u32(Rv2, TF_SHIFT_O1);
                Rv3 = vshrq_n_u32(Rv3, TF_SHIFT_O1);
                Rv4 = vshrq_n_u32(Rv4, TF_SHIFT_O1);

                Rv1 = vmlaq_u32(Bv1, Fv1, Rv1);
                Rv2 = vmlaq_u32(Bv2, Fv2, Rv2);
                Rv3 = vmlaq_u32(Bv3, Fv3, Rv3);
                Rv4 = vmlaq_u32(Bv4, Fv4, Rv4);

                // Renorm
                uint32x4_t Rlt1 = vcltq_u32(Rv1, vdupq_n_u32(RANS_BYTE_L)); // R<L
                uint32x4_t Rlt2 = vcltq_u32(Rv2, vdupq_n_u32(RANS_BYTE_L));
                uint32x4_t Rlt3 = vcltq_u32(Rv3, vdupq_n_u32(RANS_BYTE_L));
                uint32x4_t Rlt4 = vcltq_u32(Rv4, vdupq_n_u32(RANS_BYTE_L));
                uint32x4_t all2 = {2,2,2,2};
                // load 8 lanes of renorm data
                uint16x8_t norm12 =  vld1q_u16((uint16_t *)ptr);
                // move ptr by no. renorm lanes used
                ptr += vaddvq_u32(vandq_u32(Rlt1, all2))
                    +  vaddvq_u32(vandq_u32(Rlt2, all2));
                uint16x8_t norm34 =  vld1q_u16((uint16_t *)ptr);
                ptr += vaddvq_u32(vandq_u32(Rlt3, all2))
                    +  vaddvq_u32(vandq_u32(Rlt4, all2));

                // Compute lookup table index
                uint32x4_t bit = {8,4,2,1};
                uint32_t imask1 = vaddvq_u32(vandq_u32(Rlt1, bit));
                uint32_t imask2 = vaddvq_u32(vandq_u32(Rlt2, bit));
                uint32_t imask3 = vaddvq_u32(vandq_u32(Rlt3, bit));
                uint32_t imask4 = vaddvq_u32(vandq_u32(Rlt4, bit));

                uint32_t imask12 = (imask1<<4)|imask2;
                uint32_t imask34 = (imask3<<4)|imask4;

                // Shuffle norm to the corresponding R lanes, via imask
                // #define for brevity and formatting
                uint16x4_t norm1, norm2, norm3, norm4;
                norm1 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm12),
                                              idx [imask1]));
                norm2 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm12),
                                              idx2[imask12]));
                norm3 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm34),
                                              idx [imask3]));
                norm4 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm34),
                                              idx2[imask34]));

                // Add norm to R<<16 and blend back in with R
                uint32x4_t Rsl1 = vshlq_n_u32(Rv1, 16); // Rsl = R << 16
                uint32x4_t Rsl2 = vshlq_n_u32(Rv2, 16);
                uint32x4_t Rsl3 = vshlq_n_u32(Rv3, 16);
                uint32x4_t Rsl4 = vshlq_n_u32(Rv4, 16);

                Rsl1 = vaddw_u16(Rsl1, norm1);          // Rsl += norm
                Rsl2 = vaddw_u16(Rsl2, norm2);
                Rsl3 = vaddw_u16(Rsl3, norm3);
                Rsl4 = vaddw_u16(Rsl4, norm4);
//              Rsl1 = vaddq_u32(Rsl1, vmovl_u16(norm1));
//              Rsl2 = vaddq_u32(Rsl2, vmovl_u16(norm2));
//              Rsl3 = vaddq_u32(Rsl3, vmovl_u16(norm3));
//              Rsl4 = vaddq_u32(Rsl4, vmovl_u16(norm4));

                Rv1 = vbslq_u32(Rlt1, Rsl1, Rv1);       // R = R<L ? Rsl : R
                Rv2 = vbslq_u32(Rlt2, Rsl2, Rv2);
                Rv3 = vbslq_u32(Rlt3, Rsl3, Rv3);
                Rv4 = vbslq_u32(Rlt4, Rsl4, Rv4);

                vst1q_u32(&R[z+ 0], Rv1);
                vst1q_u32(&R[z+ 4], Rv2);
                vst1q_u32(&R[z+ 8], Rv3);
                vst1q_u32(&R[z+12], Rv4);
            }

            i4[0]++;
            if (++tidx == 32) {
                i4[0] -= 32;

                transpose_and_copy(out, i4, tbuf);
                tidx = 0;
            }
        }

        i4[0]-=tidx;
        int T;
        for (z = 0; z < NX; z++)
            for (T = 0; T < tidx; T++)
                out[i4[z]++] = tbuf[T][z];

        // Scalar version for close to end of in[] array so we don't do
        // SIMD loads beyond the end of the buffer
        for (; i4[0] < isz4; ) {
            for (z = 0; z < NX; z++) {
                uint32_t m = R[z] & ((1u<<TF_SHIFT_O1)-1);
                unsigned char c = sfb[l[z]][m];
                out[i4[z]++] = c;
                R[z] = fb[l[z]][c].u.s.f * (R[z]>>TF_SHIFT_O1) + m - fb[l[z]][c].u.s.b;
                RansDecRenormSafe(&R[z], &ptr, ptr_end);
                l[z] = c;
            }
        }


        // Remainder
        for (; i4[NX-1] < out_sz; i4[NX-1]++) {
            uint32_t m = R[NX-1] & ((1u<<TF_SHIFT_O1)-1);
            unsigned char c = sfb[l[NX-1]][m];
            out[i4[NX-1]] = c;
            R[NX-1] = fb[l[NX-1]][c].u.s.f * (R[NX-1]>>TF_SHIFT_O1) + m - fb[l[NX-1]][c].u.s.b;
            RansDecRenormSafe(&R[NX-1], &ptr, ptr_end);
            l[NX-1] = c;
        }
    } else {
        // TF_SHIFT_O1 = 10
        const uint32_t mask = ((1u << TF_SHIFT_O1_FAST)-1);
        uint32x4_t maskv = vdupq_n_u32((1u << TF_SHIFT_O1_FAST)-1);

        // FIXME: plus room for "safe" renorm.
        // Follow with 2nd copy doing scalar code instead?
        unsigned char tbuf[32][32];
        int tidx = 0;

        uint32x4_t RV[8] = {
            vld1q_u32(&R[0]),
            vld1q_u32(&R[4]),
            vld1q_u32(&R[8]),
            vld1q_u32(&R[12]),
            vld1q_u32(&R[16]),
            vld1q_u32(&R[20]),
            vld1q_u32(&R[24]),
            vld1q_u32(&R[28]),
        };

//      uint32x4_t MV[8] = {
//            vandq_u32(RV[0], maskv),
//            vandq_u32(RV[1], maskv),
//            vandq_u32(RV[2], maskv),
//            vandq_u32(RV[3], maskv),
//            vandq_u32(RV[4], maskv),
//            vandq_u32(RV[5], maskv),
//            vandq_u32(RV[6], maskv),
//            vandq_u32(RV[7], maskv),
//      };

        uint32_t m[NX];
        for (z = 0; z < NX; z++)
            m[z] = l[z]*TOTFREQ_O1_FAST + (R[z] & mask);

        uint32_t *S3 = (uint32_t *)s3;
        
        for (; i4[0] < isz4 && ptr+64 < ptr_end;) {
            int Z = 0;
            for (z = 0; z < NX; z+=16, Z+=4) {
                // streamline these.  Could swap between two banks and pre-load
                uint32x4_t Sv1, Sv2, Sv3, Sv4;
                uint32x4_t Fv1, Fv2, Fv3, Fv4;
                uint32x4_t Bv1, Bv2, Bv3, Bv4;
                uint32x2_t s1a, s1b, s2a, s2b, s3a, s3b, s4a, s4b;

                s1a = vcreate_u32((uint64_t)(S3[m[z+1]])<<32  | (S3[m[z+0]]));
                s1b = vcreate_u32((uint64_t)(S3[m[z+3]])<<32  | (S3[m[z+2]]));
                s2a = vcreate_u32((uint64_t)(S3[m[z+5]])<<32  | (S3[m[z+4]]));
                s2b = vcreate_u32((uint64_t)(S3[m[z+7]])<<32  | (S3[m[z+6]]));
                s3a = vcreate_u32((uint64_t)(S3[m[z+9]])<<32  | (S3[m[z+8]]));
                s3b = vcreate_u32((uint64_t)(S3[m[z+11]])<<32 | (S3[m[z+10]]));
                s4a = vcreate_u32((uint64_t)(S3[m[z+13]])<<32 | (S3[m[z+12]]));
                s4b = vcreate_u32((uint64_t)(S3[m[z+15]])<<32 | (S3[m[z+14]]));

                Sv1 = vcombine_u32(s1a, s1b);
                Sv2 = vcombine_u32(s2a, s2b);
                Sv3 = vcombine_u32(s3a, s3b);
                Sv4 = vcombine_u32(s4a, s4b);
                
                uint16x4_t p16_1 = vmovn_u32(Sv1);
                uint16x4_t p16_2 = vmovn_u32(Sv2);
                uint16x4_t p16_3 = vmovn_u32(Sv3);
                uint16x4_t p16_4 = vmovn_u32(Sv4);

                uint8x8_t  p8_12  = vmovn_u16(vcombine_u16(p16_1,p16_2));
                uint8x8_t  p8_34  = vmovn_u16(vcombine_u16(p16_3,p16_4));
                uint8x16_t p8_a   = vcombine_u8(p8_12, p8_34);
                vst1q_u8(l+z, p8_a);

                Fv1 = vshrq_n_u32(Sv1, TF_SHIFT_O1_FAST+8);
                Fv2 = vshrq_n_u32(Sv2, TF_SHIFT_O1_FAST+8);
                Fv3 = vshrq_n_u32(Sv3, TF_SHIFT_O1_FAST+8);
                Fv4 = vshrq_n_u32(Sv4, TF_SHIFT_O1_FAST+8);

                Bv1 = vandq_u32(vshrq_n_u32(Sv1, 8), maskv);
                Bv2 = vandq_u32(vshrq_n_u32(Sv2, 8), maskv);
                Bv3 = vandq_u32(vshrq_n_u32(Sv3, 8), maskv);
                Bv4 = vandq_u32(vshrq_n_u32(Sv4, 8), maskv);

                // Add in transpose here.
                memcpy(&tbuf[tidx][z], &l[z], 16);
                
                RV[Z+0] = vshrq_n_u32(RV[Z+0], TF_SHIFT_O1_FAST);
                RV[Z+1] = vshrq_n_u32(RV[Z+1], TF_SHIFT_O1_FAST);
                RV[Z+2] = vshrq_n_u32(RV[Z+2], TF_SHIFT_O1_FAST);
                RV[Z+3] = vshrq_n_u32(RV[Z+3], TF_SHIFT_O1_FAST);

                // Ready for use in S3[] offset
                Sv1 = vshlq_n_u32(vandq_u32(Sv1, vdupq_n_u32(0xff)), TF_SHIFT_O1_FAST);
                Sv2 = vshlq_n_u32(vandq_u32(Sv2, vdupq_n_u32(0xff)), TF_SHIFT_O1_FAST);
                Sv3 = vshlq_n_u32(vandq_u32(Sv3, vdupq_n_u32(0xff)), TF_SHIFT_O1_FAST);
                Sv4 = vshlq_n_u32(vandq_u32(Sv4, vdupq_n_u32(0xff)), TF_SHIFT_O1_FAST);

                RV[Z+0] = vmlaq_u32(Bv1, Fv1, RV[Z+0]);
                RV[Z+1] = vmlaq_u32(Bv2, Fv2, RV[Z+1]);
                RV[Z+2] = vmlaq_u32(Bv3, Fv3, RV[Z+2]);
                RV[Z+3] = vmlaq_u32(Bv4, Fv4, RV[Z+3]);

                // Renorm
                uint32x4_t Rlt1 = vcltq_u32(RV[Z+0], vdupq_n_u32(RANS_BYTE_L));
                uint32x4_t Rlt2 = vcltq_u32(RV[Z+1], vdupq_n_u32(RANS_BYTE_L));
                uint32x4_t Rlt3 = vcltq_u32(RV[Z+2], vdupq_n_u32(RANS_BYTE_L));
                uint32x4_t Rlt4 = vcltq_u32(RV[Z+3], vdupq_n_u32(RANS_BYTE_L));

                // Compute lookup table index
                static int nbits[16] = { 0,2,2,4, 2,4,4,6, 2,4,4,6, 4,6,6,8 };
                uint32x4_t bit = {8,4,2,1};
                uint32_t imask1 = vaddvq_u32(vandq_u32(Rlt1, bit));
                uint32_t imask2 = vaddvq_u32(vandq_u32(Rlt2, bit));
                uint32_t imask3 = vaddvq_u32(vandq_u32(Rlt3, bit));
                uint32_t imask4 = vaddvq_u32(vandq_u32(Rlt4, bit));

                // load 8 lanes of renorm data
                uint16x8_t norm12 =  vld1q_u16((uint16_t *)ptr);
                // move ptr by no. renorm lanes used
                ptr += nbits[imask1] + nbits[imask2];
                uint16x8_t norm34 =  vld1q_u16((uint16_t *)ptr);
                ptr += nbits[imask3] + nbits[imask4];

                uint32_t imask12 = (imask1<<4)|imask2;
                uint32_t imask34 = (imask3<<4)|imask4;

                // Shuffle norm to the corresponding R lanes, via imask
                // #define for brevity and formatting
                uint16x4_t norm1, norm2, norm3, norm4;
                norm1 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm12),idx [imask1]));
                norm2 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm12),idx2[imask12]));
                norm3 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm34),idx [imask3]));
                norm4 = cast_u16_u8(vqtbl1_u8(cast_u8_u16(norm34),idx2[imask34]));

                // Add norm to R<<16 and blend back in with R
                uint32x4_t Rsl1 = vshlq_n_u32(RV[Z+0], 16); // Rsl = R << 16
                uint32x4_t Rsl2 = vshlq_n_u32(RV[Z+1], 16);
                uint32x4_t Rsl3 = vshlq_n_u32(RV[Z+2], 16);
                uint32x4_t Rsl4 = vshlq_n_u32(RV[Z+3], 16);

                Rsl1 = vaddw_u16(Rsl1, norm1);          // Rsl += norm
                Rsl2 = vaddw_u16(Rsl2, norm2);
                Rsl3 = vaddw_u16(Rsl3, norm3);
                Rsl4 = vaddw_u16(Rsl4, norm4);

                RV[Z+0] = vbslq_u32(Rlt1, Rsl1, RV[Z+0]);       // R = R<L ? Rsl : R
                RV[Z+1] = vbslq_u32(Rlt2, Rsl2, RV[Z+1]);
                RV[Z+2] = vbslq_u32(Rlt3, Rsl3, RV[Z+2]);
                RV[Z+3] = vbslq_u32(Rlt4, Rsl4, RV[Z+3]);

                // Offset into s3[l][c] => s3 + l*TOTFREQ_O1_FAST + c.
                uint32x4_t off1 = vandq_u32(RV[Z+0], maskv);
                uint32x4_t off2 = vandq_u32(RV[Z+1], maskv);
                uint32x4_t off3 = vandq_u32(RV[Z+2], maskv);
                uint32x4_t off4 = vandq_u32(RV[Z+3], maskv);

                off1 = vaddq_u32(off1, Sv1);
                off2 = vaddq_u32(off2, Sv2);
                off3 = vaddq_u32(off3, Sv3);
                off4 = vaddq_u32(off4, Sv4);
                
                vst1q_u32(&m[z+ 0], off1);
                vst1q_u32(&m[z+ 4], off2);
                vst1q_u32(&m[z+ 8], off3);
                vst1q_u32(&m[z+12], off4);
            }

            i4[0]++;
            if (++tidx == 32) {
                i4[0] -= 32;

                transpose_and_copy(out, i4, tbuf);
                tidx = 0;
            }
        }

        vst1q_u32(&R[ 0], RV[0]);
        vst1q_u32(&R[ 4], RV[1]);
        vst1q_u32(&R[ 8], RV[2]);
        vst1q_u32(&R[12], RV[3]);
        vst1q_u32(&R[16], RV[4]);
        vst1q_u32(&R[20], RV[5]);
        vst1q_u32(&R[24], RV[6]);
        vst1q_u32(&R[28], RV[7]);

        i4[0]-=tidx;
        int T;
        for (z = 0; z < NX; z++)
            for (T = 0; T < tidx; T++)
                out[i4[z]++] = tbuf[T][z];

        // Scalar version for close to end of in[] array so we don't do
        // SIMD loads beyond the end of the buffer
        for (; i4[0] < isz4; ) {
            for (z = 0; z < NX; z++) {
                uint32_t m = R[z] & ((1u<<TF_SHIFT_O1_FAST)-1);
                uint32_t S = s3[l[z]][m];
                unsigned char c = S & 0xff;
                out[i4[z]++] = c;
                R[z] = (S>>(TF_SHIFT_O1_FAST+8)) * (R[z]>>TF_SHIFT_O1_FAST) +
                    ((S>>8) & ((1u<<TF_SHIFT_O1_FAST)-1));
                RansDecRenormSafe(&R[z], &ptr, ptr_end);
                l[z] = c;
            }
        }

        // Remainder
        for (; i4[NX-1] < out_sz; i4[NX-1]++) {
            uint32_t S = s3[l[NX-1]][R[NX-1] & ((1u<<TF_SHIFT_O1_FAST)-1)];
            out[i4[NX-1]] = l[NX-1] = S&0xff;
            R[NX-1] = (S>>(TF_SHIFT_O1_FAST+8)) * (R[NX-1]>>TF_SHIFT_O1_FAST)
                + ((S>>8) & ((1u<<TF_SHIFT_O1_FAST)-1));
            RansDecRenormSafe(&R[NX-1], &ptr, ptr_end);
        }
    }
    //fprintf(stderr, "    1 Decoded %d bytes\n", (int)(ptr-in)); //c-size

    htscodecs_tls_free(sfb_);
    return out;

 err:
    htscodecs_tls_free(sfb_);
    free(out_free);
    free(c_freq);

    return NULL;
}

#undef MAGIC2
#else  /* __ARM_NEON */
// Prevent "empty translation unit" errors when building without NEON
const char *rANS_static32x16pr_neon_disabled = "No NEON";
#endif /* __ARM_NEON */
