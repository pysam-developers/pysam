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

#ifndef RANS_STATIC32x16PR_H
#define RANS_STATIC32x16PR_H

/*
 * This header contains standard scalar implementations of the 32-way
 * unrolled rANS codec as well as declarations for the custom SIMD
 * implementations of x86_64 and Arm Aarch64 CPUs.
 *
 * The AVX2 and AVX512 source files need to be compiled separately as
 * we have per-file -march= compiler options and we don't wish to
 * accidentally get AVX instructions in the scalar variant.  The x86-64
 * binary then contains all 3 variants at the same time and selected
 * automatically at run time.
 *
 * The ARM Neon version is currently different, as we don't have any test
 * machines without this capability. I think it's default on all 64-bit
 * CPUs, so it's not something we're concerned with.  For simplicity
 * therefore the ARM code is simply #included into the scalar file.
 */

#ifdef __cplusplus
extern "C" {
#endif

//----------------------------------------------------------------------
// Standard scalar versions
unsigned char *rans_compress_O0_32x16(unsigned char *in,
                                      unsigned int in_size,
                                      unsigned char *out,
                                      unsigned int *out_size);

unsigned char *rans_uncompress_O0_32x16(unsigned char *in,
                                        unsigned int in_size,
                                        unsigned char *out,
                                        unsigned int out_sz);

unsigned char *rans_compress_O1_32x16(unsigned char *in,
                                      unsigned int in_size,
                                      unsigned char *out,
                                      unsigned int *out_size);

unsigned char *rans_uncompress_O1_32x16(unsigned char *in,
                                        unsigned int in_size,
                                        unsigned char *out,
                                        unsigned int out_sz);

//----------------------------------------------------------------------
// Intel SSE4 implementation.  Only the O0 decoder for now
#if defined(HAVE_SSE4_1) && defined(HAVE_SSSE3) && defined(HAVE_POPCNT)
unsigned char *rans_compress_O0_32x16_sse4(unsigned char *in,
                                           unsigned int in_size,
                                           unsigned char *out,
                                           unsigned int *out_size);

unsigned char *rans_uncompress_O0_32x16_sse4(unsigned char *in,
                                             unsigned int in_size,
                                             unsigned char *out,
                                             unsigned int out_sz);

unsigned char *rans_uncompress_O1_32x16_sse4(unsigned char *in,
                                             unsigned int in_size,
                                             unsigned char *out,
                                             unsigned int out_sz);
#endif

//----------------------------------------------------------------------
// Intel AVX2 implementation
#ifdef HAVE_AVX2
unsigned char *rans_compress_O0_32x16_avx2(unsigned char *in,
                                           unsigned int in_size,
                                           unsigned char *out,
                                           unsigned int *out_size);

unsigned char *rans_uncompress_O0_32x16_avx2(unsigned char *in,
                                             unsigned int in_size,
                                             unsigned char *out,
                                             unsigned int out_sz);

unsigned char *rans_compress_O1_32x16_avx2(unsigned char *in,
                                           unsigned int in_size,
                                           unsigned char *out,
                                           unsigned int *out_size);

unsigned char *rans_uncompress_O1_32x16_avx2(unsigned char *in,
                                             unsigned int in_size,
                                             unsigned char *out,
                                             unsigned int out_sz);
#endif // HAVE_AVX2

//----------------------------------------------------------------------
// Intel AVX512 implementation
#ifdef HAVE_AVX512
unsigned char *rans_compress_O0_32x16_avx512(unsigned char *in,
                                             unsigned int in_size,
                                             unsigned char *out,
                                             unsigned int *out_size);

unsigned char *rans_uncompress_O0_32x16_avx512(unsigned char *in,
                                               unsigned int in_size,
                                               unsigned char *out,
                                               unsigned int out_sz);

unsigned char *rans_compress_O1_32x16_avx512(unsigned char *in,
                                             unsigned int in_size,
                                             unsigned char *out,
                                             unsigned int *out_size);

unsigned char *rans_uncompress_O1_32x16_avx512(unsigned char *in,
                                               unsigned int in_size,
                                               unsigned char *out,
                                               unsigned int out_sz);
#endif // HAVE_AVX512

//----------------------------------------------------------------------
// Arm Neon implementation
#ifdef __ARM_NEON
unsigned char *rans_compress_O0_32x16_neon(unsigned char *in,
                                           unsigned int in_size,
                                           unsigned char *out,
                                           unsigned int *out_size);

unsigned char *rans_uncompress_O0_32x16_neon(unsigned char *in,
                                             unsigned int in_size,
                                             unsigned char *out,
                                             unsigned int out_sz);

unsigned char *rans_compress_O1_32x16_neon(unsigned char *in,
                                           unsigned int in_size,
                                           unsigned char *out,
                                           unsigned int *out_size);

unsigned char *rans_uncompress_O1_32x16_neon(unsigned char *in,
                                             unsigned int in_size,
                                             unsigned char *out,
                                             unsigned int out_sz);
#endif // ARM_NEON

#ifdef __cplusplus
}
#endif

#endif /* RANS_STATIC32x16PR_H */
