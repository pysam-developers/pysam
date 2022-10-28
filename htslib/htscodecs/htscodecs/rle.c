/*
 * Copyright (c) 2019-2021 Genome Research Ltd.
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

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "varint.h"
#include "rle.h"

#define MAGIC 8

//-----------------------------------------------------------------------------
// Auto compute rle_syms / rle_nsyms
static void rle_find_syms(uint8_t *data, uint64_t data_len,
                          int64_t *saved, // dim >= 256 
                          uint8_t *rle_syms, int *rle_nsyms) {
    int last = -1, n;
    uint64_t i;

    if (data_len > 256) {
        // 186/450
        // Interleaved buffers to avoid cache collisions
        int64_t saved2[256+MAGIC] = {0};
        int64_t saved3[256+MAGIC] = {0};
        int64_t saved4[256+MAGIC] = {0};
        int64_t len4 = data_len&~3;
        for (i = 0; i < len4; i+=4) {
            int d1 = (data[i+0] == last)     <<1;
            int d2 = (data[i+1] == data[i+0])<<1;
            int d3 = (data[i+2] == data[i+1])<<1;
            int d4 = (data[i+3] == data[i+2])<<1;
            last = data[i+3];
            saved [data[i+0]] += d1-1;
            saved2[data[i+1]] += d2-1;
            saved3[data[i+2]] += d3-1;
            saved4[data[i+3]] += d4-1;
        }
        while (i < data_len) {
            int d = (data[i] == last)<<1;
            saved[data[i]] += d - 1;
            last = data[i];
            i++;
        }
        for (i = 0; i < 256; i++)
            saved[i] += saved2[i] + saved3[i] + saved4[i];
    } else {
        // 163/391
        for (i = 0; i < data_len; i++) {
            if (data[i] == last) {
                saved[data[i]]++;
            } else {
                saved[data[i]]--;
                last = data[i];
            }
        }
    }

    // Map back to a list
    for (i = n = 0; i < 256; i++) {
        if (saved[i] > 0)
            rle_syms[n++] = i;
    }
    *rle_nsyms = n;
}

uint8_t *hts_rle_encode(uint8_t *data, uint64_t data_len,
                        uint8_t *run,  uint64_t *run_len,
                        uint8_t *rle_syms, int *rle_nsyms,
                        uint8_t *out, uint64_t *out_len) {
    uint64_t i, j, k;
    if (!out)
        if (!(out = malloc(data_len*2)))
            return NULL;

    // Two pass:  Firstly compute which symbols are worth using RLE on.
    int64_t saved[256+MAGIC] = {0};

    if (*rle_nsyms) {
        for (i = 0; i < *rle_nsyms; i++)
            saved[rle_syms[i]] = 1;
    } else {
        // Writes back to rle_syms and rle_nsyms
        rle_find_syms(data, data_len, saved, rle_syms, rle_nsyms);
    }

    // 2nd pass: perform RLE itself to out[] and run[] arrays.
    for (i = j = k = 0; i < data_len; i++) {
        out[k++] = data[i];
        if (saved[data[i]] > 0) {
            int rlen = i;
            int last = data[i];
            while (i < data_len && data[i] == last)
                i++;
            i--;
            rlen = i-rlen;

            j += var_put_u32(&run[j], NULL, rlen);
        }
    }
    
    *run_len = j;
    *out_len = k;
    return out;
}

// On input *out_len holds the allocated size of out[].
// On output it holds the used size of out[].
uint8_t *hts_rle_decode(uint8_t *lit, uint64_t lit_len,
                        uint8_t *run, uint64_t run_len,
                        uint8_t *rle_syms, int rle_nsyms,
                        uint8_t *out, uint64_t *out_len) {
    uint64_t j;
    uint8_t *run_end = run + run_len;

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (*out_len > 100000)
        return NULL;
#endif

    int saved[256] = {0};
    for (j = 0; j < rle_nsyms; j++)
        saved[rle_syms[j]] = 1;

    uint8_t *lit_end = lit + lit_len;
    uint8_t *out_end = out + *out_len;
    uint8_t *outp = out;

    while (lit < lit_end) {
        if (outp >= out_end)
            goto err;

        uint8_t b = *lit;
        if (saved[b]) {
            uint32_t rlen;
            run += var_get_u32(run, run_end, &rlen);
            if (rlen) {
                if (outp + rlen >= out_end)
                    goto err;
                memset(outp, b, rlen+1);
                outp += rlen+1;
            } else {
                *outp++ = b;
            }
        } else {
            *outp++ = b;
        }
        lit++;
    }

    *out_len = outp-out;
    return out;

 err:
    return NULL;
}

// Deprecated interface; to remove when we next to an ABI breakage
uint8_t *rle_encode(uint8_t *data, uint64_t data_len,
                    uint8_t *run,  uint64_t *run_len,
                    uint8_t *rle_syms, int *rle_nsyms,
                    uint8_t *out, uint64_t *out_len) {
    return hts_rle_encode(data, data_len, run, run_len,
                          rle_syms, rle_nsyms, out, out_len);
}

// Deprecated interface; to remove when we next to an ABI breakage
uint8_t *rle_decode(uint8_t *lit, uint64_t lit_len,
                    uint8_t *run, uint64_t run_len,
                    uint8_t *rle_syms, int rle_nsyms,
                    uint8_t *out, uint64_t *out_len) {
    return hts_rle_decode(lit, lit_len, run, run_len,
                          rle_syms, rle_nsyms, out, out_len);
}
