/*
 * Copyright (c) 2019-2020, 2022 Genome Research Ltd.
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
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "pack.h"

//-----------------------------------------------------------------------------

/*
 * Packs multiple symbols into a single byte if the total alphabet of symbols
 * used is <= 16.  Each new symbol takes up 1, 2, 4 or 8 bits, or 0 if the
 * alphabet used is 1 (constant).
 *
 * If successful, out_meta/out_meta_len are set to hold the mapping table
 * to be used during decompression.
 *
 * Returns the packed buffer on success with new length in out_len,
 *         NULL of failure
 */
uint8_t *hts_pack(uint8_t *data, int64_t len,
                  uint8_t *out_meta, int *out_meta_len, uint64_t *out_len) {
    int p[256] = {0}, n;
    uint64_t i, j;

    // count syms
    for (i = 0; i < len; i++)
        p[data[i]]=1;
    
    for (i = n = 0; i < 256; i++) {
        if (p[i]) {
            p[i] = n++; // p[i] is now the code number
            out_meta[n] = i;
        }
    }
    out_meta[0] = n; // 256 wraps to 0
    j = n+1;

    // 1 value per byte
    if (n > 16)
        return NULL;

    uint8_t *out = malloc(len+1);
    if (!out)
        return NULL;

    // Work out how many values per byte to encode.
    int val_per_byte;
    if (n > 4)
        val_per_byte = 2;
    else if (n > 2)
        val_per_byte = 4;
    else if (n > 1)
        val_per_byte = 8;
    else
        val_per_byte = 0; // infinite

    *out_meta_len = j;
    j = 0;

    switch (val_per_byte) {
    case 2:
        for (i = 0; i < (len & ~1); i+=2)
            out[j++] = (p[data[i]]<<0) | (p[data[i+1]]<<4);
        switch (len-i) {
        case 1: out[j++] = p[data[i]];
        }
        *out_len = j;
        return out;

    case 4: {
        for (i = 0; i < (len & ~3); i+=4)
            out[j++] = (p[data[i]]<<0) | (p[data[i+1]]<<2) | (p[data[i+2]]<<4) | (p[data[i+3]]<<6);
        out[j] = 0;
        int s = len-i, x = 0;
        switch (s) {
        case 3: out[j] |= p[data[i++]] << x; x+=2; // fall-through
        case 2: out[j] |= p[data[i++]] << x; x+=2; // fall-through
        case 1: out[j] |= p[data[i++]] << x; x+=2;
            j++;
        }
        *out_len = j;
        return out;
    }

    case 8: {
        for (i = 0; i < (len & ~7); i+=8)
            out[j++] = (p[data[i+0]]<<0) | (p[data[i+1]]<<1) | (p[data[i+2]]<<2) | (p[data[i+3]]<<3)
                     | (p[data[i+4]]<<4) | (p[data[i+5]]<<5) | (p[data[i+6]]<<6) | (p[data[i+7]]<<7);
        out[j] = 0;
        int s = len-i, x = 0;
        switch (s) {
        case 7: out[j] |= p[data[i++]] << x++; // fall-through
        case 6: out[j] |= p[data[i++]] << x++; // fall-through
        case 5: out[j] |= p[data[i++]] << x++; // fall-through
        case 4: out[j] |= p[data[i++]] << x++; // fall-through
        case 3: out[j] |= p[data[i++]] << x++; // fall-through
        case 2: out[j] |= p[data[i++]] << x++; // fall-through
        case 1: out[j] |= p[data[i++]] << x++;
            j++;
        }
        *out_len = j;
        return out;
    }

    case 0:
        *out_len = j;
        return out;
    }

    return NULL;
}


/*
 * Unpacks the meta-data portions of the hts_pack algorithm.
 * This consists of the count of symbols and their values.
 *
 * The "map" array is filled out with the used symbols.
 * "nsym" is set to contain the number of symbols per byte;
 * 0, 1, 2, 4 or 8.
 *
 * Returns number of bytes of data[] consumed on success,
 *         zero on failure.
 */
uint8_t hts_unpack_meta(uint8_t *data, uint32_t data_len,
                        uint64_t udata_len, uint8_t *map, int *nsym) {
    if (data_len == 0)
        return 0;

    // Number of symbols used
    unsigned int n = data[0];
    if (n == 0)
        n = 256;

    // Symbols per byte
    if (n <= 1)
        *nsym = 0;
    else if (n <= 2)
        *nsym = 8;
    else if (n <= 4)
        *nsym = 4;
    else if (n <= 16)
        *nsym = 2;
    else {
        *nsym = 1; // no packing
        return 1;
    }

    if (data_len <= 1)
        return 0;

    int j = 1, c = 0;
    do {
        map[c++] = data[j++];
    } while (c < n && j < data_len);

    return c < n ? 0 : j;
}

/*
 * Unpacks a packed data steam (given the unpacked meta-data).
 *
 * "map" is the pack map, mapping 0->n to the expanded symbols.
 * The "out" buffer must be preallocated by the caller to be the correct
 * size.  For error checking purposes, out_len is set to the size of
 * this buffer.
 *
 * Returns uncompressed data (out) on success,
 *         NULL on failure.
 */
uint8_t *hts_unpack(uint8_t *data, int64_t len, uint8_t *out, uint64_t out_len, int nsym, uint8_t *p) {
    //uint8_t *out;
    uint8_t c = 0;
    int64_t i, j = 0, olen;

    if (nsym == 1) {
        // raw data; FIXME: shortcut the need for malloc & memcpy here
        memcpy(out, data, len);
        return out;
    }

    switch(nsym) {
    case 8: {
        union {
            uint64_t w;
            uint8_t c[8];
        } map[256];
        int x;
        for (x = 0; x < 256; x++) {
            map[x].c[0] = p[x>>0&1];
            map[x].c[1] = p[x>>1&1];
            map[x].c[2] = p[x>>2&1];
            map[x].c[3] = p[x>>3&1];
            map[x].c[4] = p[x>>4&1];
            map[x].c[5] = p[x>>5&1];
            map[x].c[6] = p[x>>6&1];
            map[x].c[7] = p[x>>7&1];
        }
        if ((out_len+7)/8 > len)
            return NULL;
        olen = out_len & ~7;

        for (i = 0; i < olen; i+=8)
            memcpy(&out[i], &map[data[j++]].w, 8);

        if (out_len != olen) {
            c = data[j++];
            while (i < out_len) {
                out[i++] = p[c & 1];
                c >>= 1;
            }
        }
        break;
    }

    case 4: {
        union {
            uint32_t w;
            uint8_t c[4];
        } map[256];

        int x, y, z, _, P=0;
        for (x = 0; x < 4; x++)
            for (y = 0; y < 4; y++)
                for (z = 0; z < 4; z++)
                    for (_ = 0; _ < 4; _++, P++) {
                        map[P].c[0] = p[_];
                        map[P].c[1] = p[z];
                        map[P].c[2] = p[y];
                        map[P].c[3] = p[x];
                    }

        if ((out_len+3)/4 > len)
            return NULL;
        olen = out_len & ~3;

        for (i = 0; i < olen-12; i+=16) {
            uint32_t w[] = {
                map[data[j+0]].w,
                map[data[j+1]].w,
                map[data[j+2]].w,
                map[data[j+3]].w
            };
            j += 4;
            memcpy(&out[i], &w, 16);
        }

        for (; i < olen; i+=4)
            memcpy(&out[i], &map[data[j++]].w, 4);

        if (out_len != olen) {
            c = data[j++];
            while (i < out_len) {
                out[i++] = p[c & 3];
                c >>= 2;
            }
        }
        break;
    }

    case 2: {
        union {
            uint16_t w;
            uint8_t c[2];
        } map[256];

        int x, y;
        for (x = 0; x < 16; x++) {
            for (y = 0; y < 16; y++) {
                map[x*16+y].c[0] = p[y];
                map[x*16+y].c[1] = p[x];
            }
        }

        if ((out_len+1)/2 > len)
            return NULL;
        olen = out_len & ~1;

        for (i = j = 0; i+2 < olen; i+=4) {
            uint16_t w[] = {
                map[data[j+0]].w,
                map[data[j+1]].w
            };
            memcpy(&out[i], &w, 4);

            j += 2;
        }

        for (; i < olen; i+=2)
            memcpy(&out[i], &map[data[j++]].w, 2);

        if (out_len != olen) {
            c = data[j++];
            out[i+0] = p[c&15];
        }
        break;
    }

    case 0:
        memset(out, p[0], out_len);
        break;

    default:
        return NULL;
    }

    return out;
}


uint8_t *hts_unpack_(uint8_t *data, int64_t len, uint8_t *out, uint64_t out_len, int nsym, uint8_t *p) {
    //uint8_t *out;
    uint8_t c = 0;
    int64_t i, j = 0, olen;

    if (nsym == 1) {
        // raw data; FIXME: shortcut the need for malloc & memcpy here
        memcpy(out, data, len);
        return out;
    }

    switch(nsym) {
    case 2: {
        uint16_t map[256], x, y;
        for (x = 0; x < 16; x++)
            for (y = 0; y < 16; y++)
                map[x*16+y] = p[x]*256+p[y];

        if ((out_len+1)/2 > len)
            return NULL;
        olen = out_len & ~1;

        uint16_t *o16 = (uint16_t *)out;
        for (i = 0; i+4 < olen/2; i+=4) {
            int k;
            for (k = 0; k < 4; k++)
                o16[i+k] = map[data[i+k]];
        }
        j = i; i *= 2;

        for (; i < olen; i+=2) {
            uint16_t w1 = map[data[j++]];
            *(uint16_t *)&out[i] = w1;
        }

        if (out_len != olen) {
            c = data[j++];
            out[i+0] = p[c&15];
        }
        break;
    }

    default:
        return NULL;
    }

    return out;
}
