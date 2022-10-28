//#include <stdio.h>

// FIXME: make get functions const uint8_t *

/*
 * Copyright (c) 2019 Genome Research Ltd.
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

#ifndef VARINT2_H
#define VARINT2_H

#include <stdint.h>

// General API scheme is var_{get,put}_{s,u}{32,64}
// s/u for signed/unsigned;  32/64 for integer size.

// The ideas here are taken from the vbenc code in TurboPFor
// (https://github.com/powturbo/TurboPFor) with analysis at
// https://github.com/stoklund/varint.

// Unlike the ITF8 and standard 7-bit at a time encodings, this
// tries to ensure a larger portion of small numbers still fit in 1 byte.
// This trades more space for long integers with less space for short ones,
// which seems like a good tradeoff given the typical distribution curves.
//
// Like ITF8 and LTF8, the first byte also indicates the total number of
// bytes we need to decode, but unlike those it uses the same format for
// both meaning changing data type doesn't change encoding.
//
// Size comparison examples.
//
//              Max value
// Bytes        ITF8/7bit               This
// 1                  127                176
// 2               16,383             16,560
// 3            2,097,151            540,848
// 4          268,435,455         16,777,215
// 5       34,359,738,368      4,294,967,296
// 6    4,398,046,511,104  1,099,511,627,776
// ...
//
// The format is as follows:
// 0-176                     1 byte:  0 + 8 bit
// 177-16560 (14 bit range)  2 bytes: 177 + 6bit, 0 + 8bit, for x-177
// 16561-540848 (19 bits)    3 bytes: 241 + 3bit, 0+8, 0+8, for x-16561
// 540849-16777215 (~24 bit) 4 bytes: 249, 0+8, 0+8, 0+8, for x
// 2^24 - 2^32-1             5 bytes: 250, 0+8 x4
// 2^32 - 2^40-1             6 bytes: 251, 0+8 x5
// 2^40 - 2^48-1             7 bytes: 252, 0+8 x6
// 2^48 - 2^56-1             8 bytes: 253, 0+8 x7
// 2^56 - 2^64-1             9 bytes: 254, 0+8 x8
//
// Hence first byte value 255 is not possible and permits future
// escape code.


// FIXME: consider returning the value and having nbytes passed in by
// reference instead of vice-versa.
//
// ie uint64_t var_get_u64(uint8_t *cp, int *nbytes)
// vs int      var_get_u64(uint8_t *cp, uint64_t *val)
//
// The return value can then be assigned to 32-bit or 64-bit type
// without need of a new function name.  The cost is we can't then
// do "cp += var_get_u32(cp, endp, &u_freq_sz);".  Maybe we can't do
// overflow detection with former? (Want 32-bit but got, say, 40 bit)


// static inline char *var_dump(const uint8_t *cp, int n) {
//     static char buf[1000];
//     int i, o = 0;
//     for (i = 0; i < n; i++)
//      o += sprintf(&buf[o], " %d", cp[i]);
//     return buf;
// }

static inline int var_put_u64(uint8_t *cp, const uint8_t *endp, uint64_t x) {
    uint8_t *op = cp;

    if (x < 177) {
        if (endp && endp - cp < 1) return 0;
        // 0 to 176 in single byte as-is
        *cp++ = x;
    } else if (x < 16561) {
        if (endp && endp - cp < 2) return 0;
        *cp++ = ((x-177)>>8)+177;
        *cp++ = x-177;
    } else if (x < 540849) {
        if (endp && endp - cp < 3) return 0;
        *cp++ = ((x-16561)>>16)+241;
        *cp++ = (x-16561)>>8;
        *cp++ = x-16561;
    } else if (x < (1<<24)) {
        if (endp && endp - cp < 4) return 0;
        *cp++ = 249;
        *cp++ = x>>16;
        *cp++ = x>>8;
        *cp++ = x;
    } else if (x < (1LL<<32)) {
        if (endp && endp - cp < 5) return 0;
        *cp++ = 250;
        *cp++ = x>>24;
        *cp++ = x>>16;
        *cp++ = x>>8;
        *cp++ = x;
    } else if (x < (1LL<<40)) {
        if (endp && endp - cp < 6) return 0;
        *cp++ = 251;
        *cp++ = x>>32;
        *cp++ = x>>24;
        *cp++ = x>>16;
        *cp++ = x>>8;
        *cp++ = x;
    } else if (x < (1LL<<48)) {
        if (endp && endp - cp < 7) return 0;
        *cp++ = 252;
        *cp++ = x>>40;
        *cp++ = x>>32;
        *cp++ = x>>24;
        *cp++ = x>>16;
        *cp++ = x>>8;
        *cp++ = x;
    } else if (x < (1LL<<56)) {
        if (endp && endp - cp < 8) return 0;
        *cp++ = 253;
        *cp++ = x>>48;
        *cp++ = x>>40;
        *cp++ = x>>32;
        *cp++ = x>>24;
        *cp++ = x>>16;
        *cp++ = x>>8;
        *cp++ = x;
    } else {
        if (endp && endp - cp < 9) return 0;
        *cp++ = 254;
        *cp++ = x>>56;
        *cp++ = x>>48;
        *cp++ = x>>40;
        *cp++ = x>>32;
        *cp++ = x>>24;
        *cp++ = x>>16;
        *cp++ = x>>8;
        *cp++ = x;
    }

//    fprintf(stderr, "Put64 %d (%s)\n", x, var_dump(op, cp-op));

    return cp-op;
}

static inline int var_put_u32(uint8_t *cp, const uint8_t *endp, uint32_t x) {
    uint8_t *op = cp;

    if (x < 177) {
        if (endp && endp - cp < 1) abort();//return 0;
        // 0 to 176 in single byte as-is
        *cp++ = x;
    } else if (x < 16561) {
        if (endp && endp - cp < 2) abort();//return 0;
        *cp++ = ((x-177)>>8)+177;
        *cp++ = x-177;
    } else if (x < 540849) {
        if (endp && endp - cp < 3) abort();//return 0;
        *cp++ = ((x-16561)>>16)+241;
        *cp++ = (x-16561)>>8;
        *cp++ = x-16561;
    } else if (x < (1<<24)) {
        if (endp && endp - cp < 4) abort();//return 0;
        *cp++ = 249;
        *cp++ = x>>16;
        *cp++ = x>>8;
        *cp++ = x;
    } else {
        if (endp && endp - cp < 5) abort();//return 0;
        *cp++ = 250;
        *cp++ = x>>24;
        *cp++ = x>>16;
        *cp++ = x>>8;
        *cp++ = x;
    }

//    fprintf(stderr, "Put32 %d (%s)\n", x, var_dump(op, cp-op));

    return cp-op;
}

static inline int var_get_u64(uint8_t *cp, const uint8_t *endp, uint64_t *i) {
    uint8_t *op = cp;
    uint64_t j = 0;

    if (endp && cp >= endp) {
        *i = 0;
        return 0;
    }
    if (*cp < 177) {
        j = *cp++;
    } else if (*cp < 241) {
        j = ((cp[0] - 177)<<8) + cp[1] + 177;
        cp += 2;
    } else if (*cp < 249) {
        j = ((cp[0] - 241)<<16) + (cp[1]<<8) + cp[2] + 16561;
        cp += 3;
    } else {
        int n = *cp++ - 249 + 3;
        while (n--)
            j = (j<<8) + *cp++;
    }

//    fprintf(stderr, "Get64 %ld (%s)\n", j, var_dump(op, cp-op));

    *i = j;
    return cp-op;
}

static inline int var_get_u32(uint8_t *cp, const uint8_t *endp, uint32_t *i) {
    uint8_t *op = cp;
    uint32_t j = 0;

    if (endp && cp >= endp) {
        *i = 0;
        return 0;
    }
    if (*cp < 177) {
        j = *cp++;
    } else if (*cp < 241) {
        j = ((cp[0] - 177)<<8) + cp[1] + 177;
        cp += 2;
    } else if (*cp < 249) {
        j = ((cp[0] - 241)<<16) + (cp[1]<<8) + cp[2] + 16561;
        cp += 3;
    } else {
        int n = *cp++ - 249 + 3;
        while (n--)
            j = (j<<8) + *cp++;
    }

//    fprintf(stderr, "Get32 %d (%s)\n", j, var_dump(op, cp-op));

    *i = j;
    return cp-op;
}

// Signed versions of the above using zig-zag integer encoding.
// This folds the sign bit into the bottom bit so we iterate
// 0, -1, +1, -2, +2, etc.
static inline int var_put_s32(uint8_t *cp, const uint8_t *endp, int32_t i) {
    return var_put_u32(cp, endp, (i << 1) ^ (i >> 31));
}
static inline int var_put_s64(uint8_t *cp, const uint8_t *endp, int64_t i) {
    return var_put_u64(cp, endp, (i << 1) ^ (i >> 63));
}

static inline int var_get_s32(uint8_t *cp, const uint8_t *endp, int32_t *i) {
    int b = var_get_u32(cp, endp, (uint32_t *)i);
    *i = (*i >> 1) ^ -(*i & 1);
    return b;
}
static inline int var_get_s64(uint8_t *cp, const uint8_t *endp, int64_t *i) {
    int b = var_get_u64(cp, endp, (uint64_t *)i);
    *i = (*i >> 1) ^ -(*i & 1);
    return b;
}

static inline int var_size_u64(uint64_t v) {
    if (v < 177)
        return 1;
    else if (v < 16561)
        return 2;
    else if (v < 540849)
        return 3;

    int i = 0;
    do {
        v >>= 8;
        i++;
    } while (v);

//    fprintf(stderr, "Size %ld (%d)\n", v, i+1);

    return i+1;
}
#define var_size_u32 var_size_u64

static inline int var_size_s64(int64_t v) {
    return var_size_u64((v >> 63) ^ (v << 1));
}
#define var_size_s32 var_size_s64

#endif /* VARINT2_H */
