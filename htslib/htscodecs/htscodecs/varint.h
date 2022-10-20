// FIXME: make get functions const uint8_t *

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

#ifndef VARINT_H
#define VARINT_H

#include <stdint.h>

#ifdef VARINT2
#include "varint2.h"
#else

// General API scheme is var_{get,put}_{s,u}{32,64}
// s/u for signed/unsigned;  32/64 for integer size.

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


// Big endian.
// Harder for encoding, but a simpler and faster decoder.
#define BIG_END
#ifdef BIG_END

static inline
int var_put_u64_safe(uint8_t *cp, const uint8_t *endp, uint64_t i) {
    uint8_t *op = cp;
    int s = 0;
    uint64_t X = i;

    // safe method when we're near end of buffer
    do {
        s += 7;
        X >>= 7;
    } while (X);

    if (endp && (endp-cp)*7 < s)
        return 0;

    int n;
    for (n = 0; n < 10; n++) {
        s -= 7;
        *cp++ = ((i>>s) & 0x7f) + (s?128:0);
        if (!s)
            break;
    }

    return cp-op;
}

// This can be optimised further with __builtin_clzl(i) and goto various
// bits of the if/else-if structure, but it's not a vast improvement and
// we are dominated by small values.  Simplicity wins for now
static inline
int var_put_u64(uint8_t *cp, const uint8_t *endp, uint64_t i) {
    if (endp && (endp-cp) < 10)
        return var_put_u64_safe(cp, endp, i);

    // maximum of 10 bytes written
    if (i < (1<<7)) {
        *cp = i;
        return 1;
    } else if (i < (1<<14)) {
        *cp++ = ((i>> 7) & 0x7f) | 128;
        *cp++ =   i      & 0x7f;
        return 2;
    } else if (i < (1<<21)) {
        *cp++ = ((i>>14) & 0x7f) | 128;
        *cp++ = ((i>> 7) & 0x7f) | 128;
        *cp++ =   i      & 0x7f;
        return 3;
    } else if (i < (1<<28)) {
        *cp++ = ((i>>21) & 0x7f) | 128;
        *cp++ = ((i>>14) & 0x7f) | 128;
        *cp++ = ((i>> 7) & 0x7f) | 128;
        *cp++ =   i      & 0x7f;
        return 4;
    } else if (i < (1LL<<35)) {
        *cp++ = ((i>>28) & 0x7f) | 128;
        *cp++ = ((i>>21) & 0x7f) | 128;
        *cp++ = ((i>>14) & 0x7f) | 128;
        *cp++ = ((i>> 7) & 0x7f) | 128;
        *cp++ =   i      & 0x7f;
        return 5;
    } else if (i < (1LL<<42)) {
        *cp++ = ((i>>35) & 0x7f) | 128;
        *cp++ = ((i>>28) & 0x7f) | 128;
        *cp++ = ((i>>21) & 0x7f) | 128;
        *cp++ = ((i>>14) & 0x7f) | 128;
        *cp++ = ((i>> 7) & 0x7f) | 128;
        *cp++ =   i      & 0x7f;
        return 6;
    } else if (i < (1LL<<49)) {
        *cp++ = ((i>>42) & 0x7f) | 128;
        *cp++ = ((i>>35) & 0x7f) | 128;
        *cp++ = ((i>>28) & 0x7f) | 128;
        *cp++ = ((i>>21) & 0x7f) | 128;
        *cp++ = ((i>>14) & 0x7f) | 128;
        *cp++ = ((i>> 7) & 0x7f) | 128;
        *cp++ =   i      & 0x7f;
        return 7;
    } else if (i < (1LL<<56)) {
        *cp++ = ((i>>49) & 0x7f) | 128;
        *cp++ = ((i>>42) & 0x7f) | 128;
        *cp++ = ((i>>35) & 0x7f) | 128;
        *cp++ = ((i>>28) & 0x7f) | 128;
        *cp++ = ((i>>21) & 0x7f) | 128;
        *cp++ = ((i>>14) & 0x7f) | 128;
        *cp++ = ((i>> 7) & 0x7f) | 128;
        *cp++ =   i      & 0x7f;
        return 8;
    } else if (i < (1LL<<63)) {
        *cp++ = ((i>>56) & 0x7f) | 128;
        *cp++ = ((i>>49) & 0x7f) | 128;
        *cp++ = ((i>>42) & 0x7f) | 128;
        *cp++ = ((i>>35) & 0x7f) | 128;
        *cp++ = ((i>>28) & 0x7f) | 128;
        *cp++ = ((i>>21) & 0x7f) | 128;
        *cp++ = ((i>>14) & 0x7f) | 128;
        *cp++ = ((i>> 7) & 0x7f) | 128;
        *cp++ =   i      & 0x7f;
        return 9;
    } else {
        *cp++ = ((i>>63) & 0x7f) | 128;
        *cp++ = ((i>>56) & 0x7f) | 128;
        *cp++ = ((i>>49) & 0x7f) | 128;
        *cp++ = ((i>>42) & 0x7f) | 128;
        *cp++ = ((i>>35) & 0x7f) | 128;
        *cp++ = ((i>>28) & 0x7f) | 128;
        *cp++ = ((i>>21) & 0x7f) | 128;
        *cp++ = ((i>>14) & 0x7f) | 128;
        *cp++ = ((i>> 7) & 0x7f) | 128;
        *cp++ =   i      & 0x7f;
    }

    return 10;
}

static inline
int var_put_u32_safe(uint8_t *cp, const uint8_t *endp, uint32_t i) {
    uint8_t *op = cp;
    int s = 0;
    uint32_t X = i;

    // safe method when we're near end of buffer
    do {
        s += 7;
        X >>= 7;
    } while (X);

    if (endp && (endp-cp)*7 < s)
        return 0;

    int n;
    for (n = 0; n < 5; n++) {
        s -= 7;
        *cp++ = ((i>>s) & 0x7f) + (s?128:0);
        if (!s)
            break;
    }

    return cp-op;
}

static inline
int var_put_u32(uint8_t *cp, const uint8_t *endp, uint32_t i) {
    if (endp && (endp-cp) < 5)
        return var_put_u32_safe(cp, endp, i);

    if (i < (1<<7)) {
        *cp = i;
        return 1;
    } else if (i < (1<<14)) {
        *cp++ = ((i>> 7) & 0x7f) | 128;
        *cp++ =   i      & 0x7f;
        return 2;
    } else if (i < (1<<21)) {
        *cp++ = ((i>>14) & 0x7f) | 128;
        *cp++ = ((i>> 7) & 0x7f) | 128;
        *cp++ =   i      & 0x7f;
        return 3;
    } else if (i < (1<<28)) {
        *cp++ = ((i>>21) & 0x7f) | 128;
        *cp++ = ((i>>14) & 0x7f) | 128;
        *cp++ = ((i>> 7) & 0x7f) | 128;
        *cp++ =   i      & 0x7f;
        return 4;
    } else {
        *cp++ = ((i>>28) & 0x7f) | 128;
        *cp++ = ((i>>21) & 0x7f) | 128;
        *cp++ = ((i>>14) & 0x7f) | 128;
        *cp++ = ((i>> 7) & 0x7f) | 128;
        *cp++ =   i      & 0x7f;
    }

    return 5;
}

static inline
int var_get_u64(uint8_t *cp, const uint8_t *endp, uint64_t *i) {
    uint8_t *op = cp, c;
    uint64_t j = 0;

    if (!endp || endp - cp >= 10) {
        int n = 10;
        do {
            c = *cp++;
            j = (j<<7) | (c & 0x7f);
        } while ((c & 0x80) && n-- > 0);
    } else {
        if (cp >= endp) {
            *i = 0;
            return 0;
        }

        do {
            c = *cp++;
            j = (j<<7) | (c & 0x7f);
        } while ((c & 0x80) && cp < endp);
    }

    *i = j;
    return cp-op;
}

static inline
int var_get_u32(uint8_t *cp, const uint8_t *endp, uint32_t *i) {
    uint8_t *op = cp, c;
    uint32_t j = 0;

    if (!endp || endp - cp >= 6) {
        // Known maximum loop count helps optimiser.
        // NB: this helps considerably at -O3 level, but may harm -O2.
        // (However we optimise for those that want optimal code.)
        int n = 5;
        do {
            c = *cp++;
            j = (j<<7) | (c & 0x7f);
        } while ((c & 0x80) && n-- > 0);
    } else {
        if (cp >= endp) {
            *i = 0;
            return 0;
        }

        if (*cp < 128) {
            *i = *cp;
            return 1;
        }

        do {
            c = *cp++;
            j = (j<<7) | (c & 0x7f);
        } while ((c & 0x80) && cp < endp);
    }

    *i = j;
    return cp-op;
}

//-----------------------------------------------------------------------------
#else // BIG_END

// Little endian 7-bit variable sized integer encoding.
// The unsigned value is equivalent to LEB128 encoding.
// For signed, see below.
// This is also the Google Protocol Buffer and WebAssembly format.
static inline int var_put_u64(uint8_t *cp, const uint8_t *endp, uint64_t i) {
    uint8_t *op = cp;

    if (!endp || (endp-cp)*7 >= 10) {
        // Unsafe or big-enough anyway
        do {
            *cp++ = (i&0x7f) + ((i>=0x80)<<7);
            i >>= 7;
        } while (i);
    } else if (cp < endp) {
        // End checked variant
        do {
            *cp++ = (i&0x7f) + ((i>=0x80)<<7);
            i >>= 7;
        } while (i && cp < endp);
    }

    return cp-op;
}

static inline int var_put_u32(uint8_t *cp, const uint8_t *endp, uint32_t i) {
    uint8_t *op = cp;

    if (!endp || (endp-cp)*7 >= 5) {
        // Unsafe or big-enough anyway
        do {
            *cp++ = (i&0x7f) + ((i>=0x80)<<7);
            i >>= 7;
        } while (i);
    } else if (cp < endp) {
        // End checked variant
        do {
            *cp++ = (i&0x7f) + ((i>=0x80)<<7);
            i >>= 7;
        } while (i && cp < endp);
    }

    return cp-op;
}

static inline int var_get_u64(uint8_t *cp, const uint8_t *endp, uint64_t *i) {
    uint8_t *op = cp, c;
    uint64_t j = 0, s = 0;

    if (endp) {
        // Safe variant
        if (cp >= endp) {
            *i = 0;
            return 0;
        }

        do {
            c = *cp++;
            j |= (c & 0x7f) << s;
            s += 7;
        } while ((c & 0x80) && cp < endp);
    } else {
        // Unsafe variant
        do {
            c = *cp++;
            j |= (c & 0x7f) << s;
            s += 7;
        } while ((c & 0x80));
    }

    *i = j;
    return cp-op;
}

static inline int var_get_u32(uint8_t *cp, const uint8_t *endp, uint32_t *i) {
    uint8_t *op = cp, c;
    uint32_t j = 0, s = 0;

    if (endp) {
        // Safe variant
        if (cp >= endp) {
            *i = 0;
            return 0;
        }

        do {
            c = *cp++;
            j |= (c & 0x7f) << s;
            s += 7;
        } while ((c & 0x80) && cp < endp);
    } else {
        // Unsafe variant
        do {
            c = *cp++;
            j |= (c & 0x7f) << s;
            s += 7;
        } while ((c & 0x80));
    }

    *i = j;
    return cp-op;
}
#endif // BIG_END

//-----------------------------------------------------------------------------
// Signed versions of the above using zig-zag integer encoding.
// This folds the sign bit into the bottom bit so we iterate
// 0, -1, +1, -2, +2, etc.
static inline int var_put_s32(uint8_t *cp, const uint8_t *endp, int32_t i) {
    return var_put_u32(cp, endp, ((uint32_t)i << 1) ^ (i >> 31));
}
static inline int var_put_s64(uint8_t *cp, const uint8_t *endp, int64_t i) {
    return var_put_u64(cp, endp, ((uint64_t)i << 1) ^ (i >> 63));
}

static inline int var_get_s32(uint8_t *cp, const uint8_t *endp, int32_t *i) {
    int b = var_get_u32(cp, endp, (uint32_t *)i);
    *i = ((uint32_t)*i >> 1) ^ -(int32_t)(*i & 1);
    return b;
}
static inline int var_get_s64(uint8_t *cp, const uint8_t *endp, int64_t *i) {
    int b = var_get_u64(cp, endp, (uint64_t *)i);
    *i = ((uint64_t)*i >> 1) ^ -(int64_t)(*i & 1);
    return b;
}

static inline int var_size_u64(uint64_t v) {
    int i = 0;
    do {
        i++;
        v >>= 7;
    } while (v);
    return i;
}
#define var_size_u32 var_size_u64

static inline int var_size_s64(int64_t v) {
    return var_size_u64(((uint64_t)v << 1) ^ (v >> 63));
}
#define var_size_s32 var_size_s64

#endif /* VARINT2 */

#endif /* VARINT_H */
