/*
 * Copyright (c) 2019,2021 Genome Research Ltd.
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

#include <string.h>

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
	    while (j < ulen-4) {
		for (k = 0; k < 4; k++)
		    out[j++] = outN[idxN[k]++];
	    }
	    break;

	case 2:
	    while (j < ulen-2) {
		for (k = 0; k < 2; k++)
		    out[j++] = outN[idxN[k]++];
	    }
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
void hist8(unsigned char *in, unsigned int in_size, uint32_t F0[256]) {
    uint32_t F1[256+MAGIC] = {0}, F2[256+MAGIC] = {0}, F3[256+MAGIC] = {0};
    uint32_t F4[256+MAGIC] = {0}, F5[256+MAGIC] = {0}, F6[256+MAGIC] = {0};
    uint32_t F7[256+MAGIC] = {0};

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

    for (i = 0; i < 256; i++)
	F0[i] += F1[i] + F2[i] + F3[i] + F4[i] + F5[i] + F6[i] + F7[i];
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
static inline
void hist1_4(unsigned char *in, unsigned int in_size,
	     uint32_t F0[256][256], uint32_t *T0) {
    uint32_t T1[256+MAGIC] = {0}, T2[256+MAGIC] = {0}, T3[256+MAGIC] = {0};

    unsigned char l = 0, c;
    unsigned char *in_end = in + in_size;

    unsigned char cc[5] = {0};
    if (in_size > 500000) {
	uint32_t F1[256][259] = {{0}};
	while (in < in_end-8) {
	    memcpy(cc, in, 4); in += 4;
	    T0[cc[4]]++; F0[cc[4]][cc[0]]++;
	    T1[cc[0]]++; F1[cc[0]][cc[1]]++;
	    T2[cc[1]]++; F0[cc[1]][cc[2]]++;
	    T3[cc[2]]++; F1[cc[2]][cc[3]]++;
	    cc[4] = cc[3];

	    memcpy(cc, in, 4); in += 4;
	    T0[cc[4]]++; F0[cc[4]][cc[0]]++;
	    T1[cc[0]]++; F1[cc[0]][cc[1]]++;
	    T2[cc[1]]++; F0[cc[1]][cc[2]]++;
	    T3[cc[2]]++; F1[cc[2]][cc[3]]++;
	    cc[4] = cc[3];
	}
	l = cc[3];

	while (in < in_end) {
	    F0[l][c = *in++]++;
	    T0[l]++;
	    l = c;
	}

	int i, j;
	for (i = 0; i < 256; i++)
	    for (j = 0; j < 256; j++)
		F0[i][j] += F1[i][j];
    } else {
	while (in < in_end-8) {
	    memcpy(cc, in, 4); in += 4;
	    T0[cc[4]]++; F0[cc[4]][cc[0]]++;
	    T1[cc[0]]++; F0[cc[0]][cc[1]]++;
	    T2[cc[1]]++; F0[cc[1]][cc[2]]++;
	    T3[cc[2]]++; F0[cc[2]][cc[3]]++;
	    cc[4] = cc[3];

	    memcpy(cc, in, 4); in += 4;
	    T0[cc[4]]++; F0[cc[4]][cc[0]]++;
	    T1[cc[0]]++; F0[cc[0]][cc[1]]++;
	    T2[cc[1]]++; F0[cc[1]][cc[2]]++;
	    T3[cc[2]]++; F0[cc[2]][cc[3]]++;
	    cc[4] = cc[3];
	}
	l = cc[3];

	while (in < in_end) {
	    F0[l][c = *in++]++;
	    T0[l]++;
	    l = c;
	}
    }

    int i;
    for (i = 0; i < 256; i++)
	T0[i]+=T1[i]+T2[i]+T3[i];
}
