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

#ifndef RANS_UTILS_H
#define RANS_UTILS_H

#include <string.h>
#include <stdlib.h>
#include <math.h>

#if defined(__GNUC__) || defined(__clang__)
#  if !defined(__clang__) && __GNUC__ >= 100
     // better still on gcc10 for O1 decode of old rans 4x8
     // gcc 10=246/205 11=243/205 12=230/197
#    define likely(x)      __builtin_expect_with_probability((x), 1, 0.99)
#  else
     // gcc 10=193/168 11=195/161 12=199/176
#    define likely(x)      __builtin_expect((x), 1)
#  endif
#  define unlikely(x)     __builtin_expect((x), 0)
#else
#  define likely(x)   (x)
#  define unlikely(x) (x)
#endif

/*
 * Allocates size bytes from the global Thread Local Storage pool.
 * This is shared by all subsequent calls within this thread.
 *
 * Note this is NOT a general purpose allocator and usage outside of this
 * library is not advised due to assumptions and limitations in the design.
 */
void *htscodecs_tls_alloc(size_t size);
void *htscodecs_tls_calloc(size_t nmemb, size_t size);
void  htscodecs_tls_free(void *ptr);


/* Fast approximate log base 2 */
static inline double fast_log(double a) {
  union { double d; long long x; } u = { a };
  return (u.x - 4606921278410026770) * 1.539095918623324e-16; /* 1 / 6497320848556798.0; */
}

/*
 * Data transpose by N.  Common to rANS4x16 and arith_dynamic decoders.
 *
 * Tuned for specific common cases of N.
 */
static inline void unstripe(unsigned char *out, unsigned char *outN,
                            unsigned int ulen, unsigned int N,
                            unsigned int idxN[256]) {
    int j = 0, k;

    if (ulen >= N) {
        switch (N) {
        case 4:
#define LLN 16
            if (ulen >= 4*LLN) {
                while (j < ulen-4*LLN) {
                    int l;
                    for (l = 0; l < LLN; l++) {
                        for (k = 0; k < 4; k++)
                            out[j+k+l*4] = outN[idxN[k]+l];
                    }
                    for (k = 0; k < 4; k++)
                        idxN[k] += LLN;
                    j += 4*LLN;
                }
            }
            while (j < ulen-4) {
                for (k = 0; k < 4; k++)
                    out[j++] = outN[idxN[k]++];
            }
#undef LLN
            break;

        case 2:
#define LLN 4
            if (ulen >= 2*LLN) {
                while (j < ulen-2*LLN) {
                    int l;
                    for (l = 0; l < LLN; l++) {
                        for (k = 0; k < 2; k++)
                            out[j++] = outN[idxN[k]+l];
                    }
                    for (k = 0; k < 2; k++)
                        idxN[k] += l;
                }
            }
            while (j < ulen-2) {
                for (k = 0; k < 2; k++)
                    out[j++] = outN[idxN[k]++];
            }
#undef LLN
            break;

        default:
            // General case, around 25% slower overall decode
            while (j < ulen-N) {
                for (k = 0; k < N; k++)
                    out[j++] = outN[idxN[k]++];
            }
            break;
        }
    }
    for (k = 0; j < ulen; k++)
        out[j++] = outN[idxN[k]++];
}

#define MAGIC 8

/*
 * Order 0 histogram construction.  8-way unrolled to avoid cache collisions.
 */
static inline
int hist8(unsigned char *in, unsigned int in_size, uint32_t F0[256]) {
    if (in_size > 500000) {
        uint32_t *f0 = htscodecs_tls_calloc((65536+37)*3, sizeof(*f0));
        if (f0 == NULL)
            return -1;
        uint32_t *f1 = f0 + 65536+37;
        uint32_t *f2 = f1 + 65536+37;

        uint32_t i, i8 = in_size & ~15;

        for (i = 0; i < i8; i+=16) {
            uint16_t i16a[4], i16b[4];
            memcpy(i16a, in+i, 8);
            f0[i16a[0]]++;
            f1[i16a[1]]++;
            f2[i16a[2]]++;
            f0[i16a[3]]++;

            memcpy(i16b, in+i+8, 8);
            f1[i16b[0]]++;
            f0[i16b[1]]++;
            f1[i16b[2]]++;
            f2[i16b[3]]++;
        }

        while (i < in_size)
            F0[in[i++]]++;

        for (i = 0; i < 65536; i++) {
            F0[i & 0xff] += f0[i] + f1[i] + f2[i];
            F0[i >> 8  ] += f0[i] + f1[i] + f2[i];
        }
        htscodecs_tls_free(f0);
    } else {
        uint32_t F1[256+MAGIC] = {0}, F2[256+MAGIC] = {0}, F3[256+MAGIC] = {0};
        uint32_t i, i8 = in_size & ~7;

        for (i = 0; i < i8; i+=8) {
            F0[in[i+0]]++;
            F1[in[i+1]]++;
            F2[in[i+2]]++;
            F3[in[i+3]]++;
            F0[in[i+4]]++;
            F1[in[i+5]]++;
            F2[in[i+6]]++;
            F3[in[i+7]]++;
        }

        while (i < in_size)
            F0[in[i++]]++;

        for (i = 0; i < 256; i++)
            F0[i] += F1[i] + F2[i] + F3[i];
    }

    return 0;
}

// Hist8 with a crude entropy (bits / byte) estimator.
static inline
double hist8e(unsigned char *in, unsigned int in_size, uint32_t F0[256]) {
    uint32_t F1[256+MAGIC] = {0}, F2[256+MAGIC] = {0}, F3[256+MAGIC] = {0};
    uint32_t F4[256+MAGIC] = {0}, F5[256+MAGIC] = {0}, F6[256+MAGIC] = {0};
    uint32_t F7[256+MAGIC] = {0};

#ifdef __GNUC__
    double e = 0, in_size_r2 = log(1.0/in_size)/log(2);
#else
    double e = 0, in_size_r2 = log(1.0/in_size);
#endif

    unsigned int i, i8 = in_size & ~7;
    for (i = 0; i < i8; i+=8) {
        F0[in[i+0]]++;
        F1[in[i+1]]++;
        F2[in[i+2]]++;
        F3[in[i+3]]++;
        F4[in[i+4]]++;
        F5[in[i+5]]++;
        F6[in[i+6]]++;
        F7[in[i+7]]++;
    }
    while (i < in_size)
        F0[in[i++]]++;

    for (i = 0; i < 256; i++) {
        F0[i] += F1[i] + F2[i] + F3[i] + F4[i] + F5[i] + F6[i] + F7[i];
#ifdef __GNUC__
        e -= F0[i] * (32 - __builtin_clz(F0[i]|1) + in_size_r2);
#else
        e -= F0[i] * (fast_log(F0[i]) + in_size_r2);
#endif
    }

#ifndef __GNUC__
    e /= log(2);
#endif
    return e/in_size;
}

/*
 * A variant of hist8 that simply marks the presence of a symbol rather
 * than its frequency.
 */
static inline
void present8(unsigned char *in, unsigned int in_size,
              uint32_t F0[256]) {
    uint32_t F1[256+MAGIC] = {0}, F2[256+MAGIC] = {0}, F3[256+MAGIC] = {0};
    uint32_t F4[256+MAGIC] = {0}, F5[256+MAGIC] = {0}, F6[256+MAGIC] = {0};
    uint32_t F7[256+MAGIC] = {0};

    unsigned int i, i8 = in_size & ~7;
    for (i = 0; i < i8; i+=8) {
        F0[in[i+0]]=1;
        F1[in[i+1]]=1;
        F2[in[i+2]]=1;
        F3[in[i+3]]=1;
        F4[in[i+4]]=1;
        F5[in[i+5]]=1;
        F6[in[i+6]]=1;
        F7[in[i+7]]=1;
    }
    while (i < in_size)
        F0[in[i++]]=1;

    for (i = 0; i < 256; i++)
        F0[i] += F1[i] + F2[i] + F3[i] + F4[i] + F5[i] + F6[i] + F7[i];
}

/*
 * Order 1 histogram construction.  4-way unrolled to avoid cache collisions.
 */
#if 1
static inline
int hist1_4(unsigned char *in, unsigned int in_size,
            uint32_t F0[256][256], uint32_t *T0) {
    unsigned char l = 0, c;
    unsigned char *in_end = in + in_size;

    unsigned char cc[5] = {0};
    if (in_size > 500000) {
        uint32_t (*F1)[259] = htscodecs_tls_calloc(256, sizeof(*F1));
        if (F1 == NULL)
            return -1;
        while (in < in_end-8) {
            memcpy(cc, in, 4); in += 4;
            F0[cc[4]][cc[0]]++;
            F1[cc[0]][cc[1]]++;
            F0[cc[1]][cc[2]]++;
            F1[cc[2]][cc[3]]++;
            cc[4] = cc[3];

            memcpy(cc, in, 4); in += 4;
            F0[cc[4]][cc[0]]++;
            F1[cc[0]][cc[1]]++;
            F0[cc[1]][cc[2]]++;
            F1[cc[2]][cc[3]]++;
            cc[4] = cc[3];
        }
        l = cc[3];

        while (in < in_end) {
            F0[l][c = *in++]++;
            l = c;
        }
        T0[l]++;

        int i, j;
        for (i = 0; i < 256; i++) {
            int tt = 0;
            for (j = 0; j < 256; j++) {
                F0[i][j] += F1[i][j];
                tt += F0[i][j];
            }
            T0[i]+=tt;
        }
        htscodecs_tls_free(F1);
    } else {
        while (in < in_end-8) {
            memcpy(cc, in, 4); in += 4;
            F0[cc[4]][cc[0]]++;
            F0[cc[0]][cc[1]]++;
            F0[cc[1]][cc[2]]++;
            F0[cc[2]][cc[3]]++;
            cc[4] = cc[3];

            memcpy(cc, in, 4); in += 4;
            F0[cc[4]][cc[0]]++;
            F0[cc[0]][cc[1]]++;
            F0[cc[1]][cc[2]]++;
            F0[cc[2]][cc[3]]++;
            cc[4] = cc[3];
        }
        l = cc[3];

        while (in < in_end) {
            F0[l][c = *in++]++;
            l = c;
        }
        T0[l]++;

        int i, j;
        for (i = 0; i < 256; i++) {
            int tt = 0;
            for (j = 0; j < 256; j++)
                tt += F0[i][j];
            T0[i]+=tt;
        }
    }

    return 0;
}

#else
// 16 bit mode, similar to O0 freq.
// This is better on some low entropy data, but generally we prefer to do
// bit-packing and/or RLE to turn it into higher-entropy data first.
//
// Kept here for posterity incase we need it again, as it's quick tricky.
static inline
int hist1_4(unsigned char *in, unsigned int in_size,
            uint32_t F0[256][256], uint32_t *T0) {
    uint32_t f0[65536+MAGIC] = {0};
    uint32_t f1[65536+MAGIC] = {0};

    uint32_t i, i8 = (in_size-1) & ~15;

    T0[0]++; f0[in[0]<<8]++;
    for (i = 0; i < i8; i+=16) {
        uint16_t i16a[16];
        memcpy(i16a,   in+i,   16); // faster in 2 as gcc recognises this
        memcpy(i16a+8, in+i+1, 16); // faster in 2 as gcc recognises this

        f0[i16a[0]]++;
        f1[i16a[1]]++;
        f0[i16a[2]]++;
        f1[i16a[3]]++;
        f0[i16a[4]]++;
        f1[i16a[5]]++;
        f0[i16a[6]]++;
        f1[i16a[7]]++;
        f0[i16a[8]]++;
        f1[i16a[9]]++;
        f0[i16a[10]]++;
        f1[i16a[11]]++;
        f0[i16a[12]]++;
        f1[i16a[13]]++;
        f0[i16a[14]]++;
        f1[i16a[15]]++;
    }

    while (i < in_size-1) {
        F0[in[i]][in[i+1]]++;
        T0[in[i+1]]++;
        i++;
    }

    for (i = 0; i < 65536; i++) {
        F0[i&0xff][i>>8] += f0[i] + f1[i];
        T0[i>>8]         += f0[i] + f1[i];
    }

    return 0;
}
#endif

#endif /* RANS_UTILS_H */
