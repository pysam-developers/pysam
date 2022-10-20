/*
 * Copyright (c) 2017-2019 Genome Research Ltd.
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

#ifndef RANS_STATIC4x16_H
#define RANS_STATIC4x16_H

#ifdef __cplusplus
extern "C" {
#endif

unsigned int rans_compress_bound_4x16(unsigned int size, int order);
unsigned char *rans_compress_to_4x16(unsigned char *in,  unsigned int in_size,
                                     unsigned char *out, unsigned int *out_size,
                                     int order);
unsigned char *rans_compress_4x16(unsigned char *in, unsigned int in_size,
                                  unsigned int *out_size, int order);
unsigned char *rans_uncompress_to_4x16(unsigned char *in,  unsigned int in_size,
                                       unsigned char *out, unsigned int *out_size);
unsigned char *rans_uncompress_4x16(unsigned char *in, unsigned int in_size,
                                    unsigned int *out_size);

// CPU detection control.  Used for testing and benchmarking.
// These bitfields control what methods are permitted to be used.
#define RANS_CPU_ENC_SSE4     (1<<0)
#define RANS_CPU_ENC_AVX2     (2<<0)
#define RANS_CPU_ENC_AVX512   (4<<0)
#define RANS_CPU_ENC_NEON     (8<<0)

#define RANS_CPU_DEC_SSE4     (1<<8)
#define RANS_CPU_DEC_AVX2     (2<<8)
#define RANS_CPU_DEC_AVX512   (4<<8)
#define RANS_CPU_DEC_NEON     (8<<8)

void rans_set_cpu(int opts);

// "Order" byte options. ORed into the order byte.
// The bottom bits are the order itself, currently
// supporting order-0 and order-1 but with expansion room
// up to order-3 (unlikely).

//--
// The values below are stored in the file format

// Pack 2,4,8 or infinite symbols into a byte.
#define RANS_ORDER_PACK   0x80

// Run length encoding with runs & lits encoded separately
#define RANS_ORDER_RLE    0x40

// Nop; for tiny segments where rANS overhead is too big
#define RANS_ORDER_CAT    0x20

// Don't store the original size; used by STRIPE mode
#define RANS_ORDER_NOSZ   0x10

// For N-byte integer data; rotate & encode N streams.
#define RANS_ORDER_STRIPE 0x08

// 32-way unrolling instead of 4-way
#define RANS_ORDER_X32    0x04

//--
// order values below are not directly part of the file format, but control
// the behaviour of the encoder.

// Bit 8-15 of order hold the stripe size (N).
// Note: N is stored separately after the order byte

// Used to disable order-0 in the STRIPE sub-methods.
#define RANS_ORDER_STRIPE_NO0 (1<<16)

// Used to request automatic selection between 4-way and 32-way
#define RANS_ORDER_SIMD_AUTO  (1<<17)

#ifdef __cplusplus
}
#endif

#endif /* RANS_STATIC4x16_H */
