#ifndef RANS_INTERNAL_H
#define RANS_INTERNAL_H

#include "config.h"
#include "varint.h"
#include "utils.h"

/*
 * Copyright (c) 2017-2022 Genome Research Ltd.
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

// Internal: common parts to all the rANSNx16pr implementations.

// As per standard rANS_static but using optional RLE or bit-packing
// techniques prior to entropy encoding.  This is a significant
// reduction in some data sets.

// top bits in order byte
#define X_PACK   0x80    // Pack 2,4,8 or infinite symbols into a byte.
#define X_RLE    0x40    // Run length encoding with runs & lits encoded separately
#define X_CAT    0x20    // Nop; for tiny segments where rANS overhead is too big
#define X_NOSZ   0x10    // Don't store the original size; used by STRIPE mode
#define X_STRIPE 0x08    // For N-byte integer data; rotate & encode N streams.
#define X_32     0x04    // 32-way unrolling instead of 4-way

// Not part of the file format, but used to direct the encoder
#define X_SIMD_AUTO 0x100 // automatically enable X_32 if we deem it worthy
#define X_SW32_ENC  0x200 // forcibly use the software version of X_32
#define X_SW32_DEC  0x400 // forcibly use the software version of X_32
#define X_NO_AVX512 0x800 // turn off avx512, but permits AVX2

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

unsigned char *rans_compress_O0_4x16(unsigned char *in, unsigned int in_size,
                                     unsigned char *out, unsigned int *out_size);
unsigned char *rans_uncompress_O0_4x16(unsigned char *in, unsigned int in_size,
                                       unsigned char *out, unsigned int out_sz);

int rans_compute_shift(uint32_t *F0, uint32_t (*F)[256], uint32_t *T,
                       uint32_t *S);

// Rounds to next power of 2.
// credit to http://graphics.stanford.edu/~seander/bithacks.html
static inline uint32_t round2(uint32_t v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}

static inline int normalise_freq(uint32_t *F, int size, uint32_t tot) {
    int m, M, j, loop = 0;
    uint64_t tr;
    if (!size)
        return 0;

 again:
    tr = ((uint64_t)tot<<31)/size + (1<<30)/size;

    for (size = m = M = j = 0; j < 256; j++) {
        if (!F[j])
            continue;

        if (m < F[j])
            m = F[j], M = j;

        if ((F[j] = (F[j]*tr)>>31) == 0)
            F[j] = 1;
        size += F[j];
//      if (F[j] == tot)
//          F[j]--;
    }

    int adjust = tot - size;
    if (adjust > 0) {
        F[M] += adjust;
    } else if (adjust < 0) {
        if (F[M] > -adjust && (loop == 1 || F[M]/2 >= -adjust)) {
            F[M] += adjust;
        } else {
            if (loop < 1) {
                loop++;
                goto again;
            }
            adjust += F[M]-1;
            F[M] = 1;
            for (j = 0; adjust && j < 256; j++) {
                if (F[j] < 2) continue;

                int d = F[j] > -adjust;
                int m = d ? adjust : 1-F[j];
                F[j]   += m;
                adjust -= m;
            }
        }
    }

    //printf("F[%d]=%d\n", M, F[M]);
    return F[M]>0 ? 0 : -1;
}

// A specialised version of normalise_freq_shift where the input size
// is already normalised to a power of 2, meaning we can just perform
// shifts instead of hard to define multiplications and adjustments.
static inline void normalise_freq_shift(uint32_t *F, uint32_t size,
                                        uint32_t max_tot) {
    if (size == 0 || size == max_tot)
        return;

    int shift = 0, i;
    while (size < max_tot)
        size*=2, shift++;

    for (i = 0; i < 256; i++)
        F[i] <<= shift;
}

// symbols only
static inline int encode_alphabet(uint8_t *cp, uint32_t *F) {
    uint8_t *op = cp;
    int rle, j;

    for (rle = j = 0; j < 256; j++) {
        if (F[j]) {
            // j
            if (rle) {
                rle--;
            } else {
                *cp++ = j;
                if (!rle && j && F[j-1])  {
                    for(rle=j+1; rle<256 && F[rle]; rle++)
                        ;
                    rle -= j+1;
                    *cp++ = rle;
                }
                //fprintf(stderr, "%d: %d %d\n", j, rle, N[j]);
            }
        }
    }
    *cp++ = 0;
    
    return cp - op;
}

static inline int decode_alphabet(uint8_t *cp, uint8_t *cp_end, uint32_t *F) {
    if (cp == cp_end)
        return 0;

    uint8_t *op = cp;
    int rle = 0;
    int j = *cp++;
    if (cp+2 >= cp_end)
        goto carefully;

    do {
        F[j] = 1;
        if (!rle && j+1 == *cp) {
            j = *cp++;
            rle = *cp++;
        } else if (rle) {
            rle--;
            j++;
            if (j > 255)
                return 0;
        } else {
            j = *cp++;
        }
    } while(j && cp+2 < cp_end);

 carefully:
    if (j) {
        do {
            F[j] = 1;
            if(cp >= cp_end) return 0;
            if (!rle && j+1 == *cp) {
                if (cp+1 >= cp_end) return 0;
                j = *cp++;
                rle = *cp++;
            } else if (rle) {
                rle--;
                j++;
                if (j > 255)
                    return 0;
            } else {
                if (cp >= cp_end) return 0;
                j = *cp++;
            }
        } while(j && cp < cp_end);
    }

    return cp - op;
}

static inline int encode_freq(uint8_t *cp, uint32_t *F) {
    uint8_t *op = cp;
    int j;

    cp += encode_alphabet(cp, F);

    for (j = 0; j < 256; j++) {
        if (F[j])
            cp += var_put_u32(cp, NULL, F[j]);
    }

    return cp - op;
}

static inline int decode_freq(uint8_t *cp, uint8_t *cp_end, uint32_t *F,
                              uint32_t *fsum) {
    if (cp == cp_end)
        return 0;

    uint8_t *op = cp;
    cp += decode_alphabet(cp, cp_end, F);

    int j, tot = 0;
    for (j = 0; j < 256; j++) {
        if (F[j]) {
            cp += var_get_u32(cp, cp_end, (unsigned int *)&F[j]);
            tot += F[j];
        }
    }

    *fsum = tot;
    return cp - op;
}


// Use the order-0 freqs in F0 to encode the order-1 stats in F.
// All symbols present in F are present in F0, but some in F0 will
// be empty in F.  Thus we run-length encode the 0 frequencies.
static inline int encode_freq_d(uint8_t *cp, uint32_t *F0, uint32_t *F) {
    uint8_t *op = cp;
    int j, dz;

    for (dz = j = 0; j < 256; j++) {
        if (F0[j]) {
            if (F[j] != 0) {
                if (dz) {
                    // Replace dz zeros with zero + dz-1 run length
                    cp -= dz-1;
                    *cp++ = dz-1;
                }
                dz = 0;
                cp += var_put_u32(cp, NULL, F[j]);
            } else {
                //fprintf(stderr, "2: j=%d F0[j]=%d, F[j]=%d, dz=%d\n", j, F0[j], F[j], dz);
                dz++;
                *cp++ = 0;
            }
        }
    }
    
    if (dz) {
        cp -= dz-1;
        *cp++ = dz-1;
    }

    return cp - op;
}

// Normalise frequency total T[i] to match TOTFREQ_O1 and encode.
// Also initialises the RansEncSymbol structs.
//
// Returns the desired TF_SHIFT; 10 or 12 bit, or -1 on error.
static inline int encode_freq1(uint8_t *in, uint32_t in_size, int Nway,
                               RansEncSymbol syms[256][256], uint8_t **cp_p) {
    int i, j, z;
    uint8_t *out = *cp_p, *cp = out;

    // Compute O1 frequency statistics
    uint32_t (*F)[256] = htscodecs_tls_calloc(256, (sizeof(*F)));
    if (!F)
        return -1;
    uint32_t T[256+MAGIC] = {0};
    int isz4 = in_size/Nway;
    if (hist1_4(in, in_size, F, T) < 0)
        goto err;
    for (z = 1; z < Nway; z++)
        F[0][in[z*isz4]]++;
    T[0]+=Nway-1;

    // Potential fix for the wrap-around bug in AVX2 O1 encoder with shift=12.
    // This occurs when we have one single symbol, giving freq=4096.
    // We fix it elsewhere for now by looking for the wrap-around.
    // See "if (1)" statements in the AVX2 code, which is an alternative
    // to the "if (0)" here.
//    if (0) {
//      int x = -1, y = -1;
//      int n1, n2;
//      for (x = 0; x < 256; x++) {
//          n1 = n2 = -1;
//          for (y = 0; y < 256; y++) {
//              if (F[x][y])
//                  n2 = n1, n1 = y;
//          }
//          if (n2!=-1 || n1 == -1)
//              continue;
//
//          for (y = 0; y < 256; y++)
//              if (!F[x][y])
//                  break;
//          assert(y<256);
//          F[x][y]++;
//          F[0][y]++; T[y]++; F0[y]=1;
//          F[0][x]++; T[x]++; F0[x]=1;
//      }
//    }

    // Encode the order-0 stats
    int tmp_T0 = T[0];
    T[0] = 1;
    *cp++ = 0; // marker for uncompressed (may change)
    cp += encode_alphabet(cp, T);
    T[0] = tmp_T0;

    // Decide between 10-bit and 12-bit freqs.
    // Fills out S[] to hold the new scaled maximum value.
    uint32_t S[256] = {0};
    int shift = rans_compute_shift(T, F, T, S);

    // Normalise so T[i] == TOTFREQ_O1
    for (i = 0; i < 256; i++) {
        unsigned int x;

        if (T[i] == 0)
            continue;

        uint32_t max_val = S[i];
        if (shift == TF_SHIFT_O1_FAST && max_val > TOTFREQ_O1_FAST)
            max_val = TOTFREQ_O1_FAST;

        if (normalise_freq(F[i], T[i], max_val) < 0)
            goto err;
        T[i]=max_val;

        // Encode our frequency array
        cp += encode_freq_d(cp, T, F[i]);

        normalise_freq_shift(F[i], T[i], 1<<shift); T[i]=1<<shift;

        // Initialise Rans Symbol struct too.
        uint32_t *F_i_ = F[i];
        for (x = j = 0; j < 256; j++) {
            RansEncSymbolInit(&syms[i][j], x, F_i_[j], shift);
            x += F_i_[j];
        }
    }

    *out = shift<<4;
    if (cp - out > 1000) {
        uint8_t *op = out;
        // try rans0 compression of header
        unsigned int u_freq_sz = cp-(op+1);
        unsigned int c_freq_sz;
        unsigned char *c_freq = rans_compress_O0_4x16(op+1, u_freq_sz, NULL,
                                                      &c_freq_sz);
        if (c_freq && c_freq_sz + 6 < cp-op) {
            *op++ |= 1; // compressed
            op += var_put_u32(op, NULL, u_freq_sz);
            op += var_put_u32(op, NULL, c_freq_sz);
            memcpy(op, c_freq, c_freq_sz);
            cp = op+c_freq_sz;
        }
        free(c_freq);
    }

    *cp_p = cp;
    htscodecs_tls_free(F);
    return shift;

 err:
    htscodecs_tls_free(F);
    return -1;
}

// Part of decode_freq1 below.  This decodes an order-1 frequency table
// using an order-0 table to determine which stats may be stored.
static inline int decode_freq_d(uint8_t *cp, uint8_t *cp_end, uint32_t *F0,
                                uint32_t *F, uint32_t *total) {
    if (cp == cp_end)
        return 0;

    uint8_t *op = cp;
    int j, dz, T = 0;

    for (j = dz = 0; j < 256 && cp < cp_end; j++) {
        //if (F0[j]) fprintf(stderr, "F0[%d]=%d\n", j, F0[j]);
        if (!F0[j])
            continue;

        uint32_t f;
        if (dz) {
            f = 0;
            dz--;
        } else {
            if (cp >= cp_end) return 0;
            cp += var_get_u32(cp, cp_end, &f);
            if (f == 0) {
                if (cp >= cp_end) return 0;
                dz = *cp++;
            }
        }
        F[j] = f;
        T += f;
    }

    if (total) *total = T;
    return cp - op;
}

typedef struct {
    uint16_t f;
    uint16_t b;
} fb_t;

// Decode order-1 frequency table, filling out various lookup tables
// in the process. (Which will depend on shift and which values have
// been passed in.)
//
// Returns the number of bytes decoded.
static inline int decode_freq1(uint8_t *cp, uint8_t *cp_end, int shift,
                               uint32_t s3 [256][TOTFREQ_O1],
                               uint32_t s3F[256][TOTFREQ_O1_FAST],
                               uint8_t *sfb[256], fb_t fb[256][256]) {
    uint8_t *cp_start = cp;
    int i, j, x;
    uint32_t F0[256] = {0};
    int fsz = decode_alphabet(cp, cp_end, F0);
    if (!fsz)
        goto err;
    cp += fsz;

    if (cp >= cp_end)
        goto err;

    // silence false gcc warnings
    if (fb) {fb [0][0].b= 0;}
    if (s3) {s3 [0][0]  = 0;}
    if (s3F){s3F[0][0]  = 0;}

    for (i = 0; i < 256; i++) {
        if (F0[i] == 0)
            continue;

        uint32_t F[256] = {0}, T = 0;
        fsz = decode_freq_d(cp, cp_end, F0, F, &T);
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

                if (sfb && shift == TF_SHIFT_O1) {
                    memset(&sfb[i][x], j, F[j]);
                    fb[i][j].f = F[j];
                    fb[i][j].b = x;
                } else if (s3 && shift == TF_SHIFT_O1) {
                    int y;
                    for (y = 0; y < F[j]; y++)
                        s3[i][y+x] = (((uint32_t)F[j])<<(shift+8)) |(y<<8) |j;
                } else if (s3F && shift == TF_SHIFT_O1_FAST) {
                    int y;
                    for (y = 0; y < F[j]; y++)
                        s3F[i][y+x] = (((uint32_t)F[j])<<(shift+8)) |(y<<8) |j;
                }

                x += F[j];
            }
        }
        if (x != (1<<shift))
            goto err;
    }

    return cp - cp_start;

 err:
    return 0;
}

// Build s3 symbol lookup table.
// This is 12 bit freq, 12 bit bias and 8 bit symbol.
static inline int rans_F_to_s3(const uint32_t *F, int shift, uint32_t *s3) {
    int j, x;
    for (j = x = 0; j < 256; j++) {
        if (F[j] && F[j] <= (1<<shift) - x) {
            uint32_t base = (((uint32_t)F[j])<<(shift+8))|j, y;
            for (y = 0; y < F[j]; y++, x++)
                s3[x] = base + (y<<8);
        }
    }

    return x == (1<<shift) ? 0 : 1;
}

#ifdef ROT32_SIMD
#include <x86intrin.h>

// Our own implementation of _mm256_set_m128i as it's not there on older
// gcc implementations.  This is basically the same thing.
static inline __m256i _mm256_set_m128ix(__m128i H, __m128i L) {
    return _mm256_insertf128_si256(_mm256_castsi128_si256(L), H, 1);
}

static inline void rot32_simd(uint8_t t[32][32], uint8_t *out, int iN[32]) {
    int z;

    __m256i lh8[32];
    for (z = 0; z < 32/2; z+=2) {
        __m256i a, b, c, d;
        a = _mm256_loadu_si256((__m256i *)&t[z*2+0]);
        b = _mm256_loadu_si256((__m256i *)&t[z*2+1]);
        c = _mm256_loadu_si256((__m256i *)&t[z*2+2]);
        d = _mm256_loadu_si256((__m256i *)&t[z*2+3]);

        lh8[z+0]  = _mm256_unpacklo_epi8(a, b);
        lh8[z+16] = _mm256_unpackhi_epi8(a, b);
        lh8[z+1]  = _mm256_unpacklo_epi8(c, d);
        lh8[z+17] = _mm256_unpackhi_epi8(c, d);
    }

    __m256i lh32[32];
    for (z = 0; z < 32/4; z+=2) {
        __m256i a, b, c, d;
        a = _mm256_unpacklo_epi16(lh8[z*4+0], lh8[z*4+1]);
        b = _mm256_unpacklo_epi16(lh8[z*4+2], lh8[z*4+3]);
        c = _mm256_unpackhi_epi16(lh8[z*4+0], lh8[z*4+1]);
        d = _mm256_unpackhi_epi16(lh8[z*4+2], lh8[z*4+3]);

        __m256i e, f, g, h;
        e = _mm256_unpacklo_epi16(lh8[(z+1)*4+0], lh8[(z+1)*4+1]);
        f = _mm256_unpacklo_epi16(lh8[(z+1)*4+2], lh8[(z+1)*4+3]);
        g = _mm256_unpackhi_epi16(lh8[(z+1)*4+0], lh8[(z+1)*4+1]);
        h = _mm256_unpackhi_epi16(lh8[(z+1)*4+2], lh8[(z+1)*4+3]);

        lh32[z+0]  = _mm256_unpacklo_epi32(a,b);
        lh32[z+8]  = _mm256_unpacklo_epi32(c,d);
        lh32[z+16] = _mm256_unpackhi_epi32(a,b);
        lh32[z+24] = _mm256_unpackhi_epi32(c,d);

        lh32[z+1+0]  = _mm256_unpacklo_epi32(e,f);
        lh32[z+1+8]  = _mm256_unpacklo_epi32(g,h);
        lh32[z+1+16] = _mm256_unpackhi_epi32(e,f);
        lh32[z+1+24] = _mm256_unpackhi_epi32(g,h);
    }

    // Final unpack 64 and store
    int idx[] = {0, 8, 4, 12, 2, 10, 6, 14};
    for (z = 0; z < 8; z++) {
        int i = idx[z];

        // Putting this here doesn't soeed things up
        __m256i a = _mm256_unpacklo_epi64(lh32[i*2+0], lh32[i*2+1]);
        __m256i b = _mm256_unpacklo_epi64(lh32[i*2+2], lh32[i*2+3]);
        __m256i c = _mm256_unpackhi_epi64(lh32[i*2+0], lh32[i*2+1]);
        __m256i d = _mm256_unpackhi_epi64(lh32[i*2+2], lh32[i*2+3]);

        __m256i p = _mm256_set_m128ix(_mm256_extracti128_si256(b,0),
                                      _mm256_extracti128_si256(a,0));
        __m256i q = _mm256_set_m128ix(_mm256_extracti128_si256(d,0),
                                      _mm256_extracti128_si256(c,0));
        __m256i r = _mm256_set_m128ix(_mm256_extracti128_si256(b,1),
                                      _mm256_extracti128_si256(a,1));
        __m256i s = _mm256_set_m128ix(_mm256_extracti128_si256(d,1),
                                      _mm256_extracti128_si256(c,1));

        _mm256_storeu_si256((__m256i *)(&out[iN[z*2+0]]),  p);
        _mm256_storeu_si256((__m256i *)(&out[iN[z*2+1]]),  q);
        _mm256_storeu_si256((__m256i *)(&out[iN[z*2+16]]), r);
        _mm256_storeu_si256((__m256i *)(&out[iN[z*2+17]]), s);
    }

    // Store
    for (z = 0; z < 32; z++)
        iN[z] += 32;
}
#endif

#endif // RANS_INTERNAL_H
