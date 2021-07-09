// FIXME: make get functions const uint8_t *

/*
 * Copyright (c) 2019,2020 Genome Research Ltd.
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
static inline int var_put_u64(uint8_t *cp, const uint8_t *endp, uint64_t i) {
    uint8_t *op = cp;
    int s = 0;
    uint64_t X = i;

    do {
	s += 7;
	X >>= 7;
    } while (X);

    if (endp && (endp-cp)*7 < s)
	return 0;

    do {
	s -= 7;
	*cp++ = ((i>>s) & 0x7f) + (s?128:0);
    } while (s);

    return cp-op;
}

static inline int var_put_u32(uint8_t *cp, const uint8_t *endp, uint32_t i) {
    uint8_t *op = cp;
    int s = 0;
    uint32_t X = i;

    do {
	s += 7;
	X >>= 7;
    } while (X);

    if (endp && (endp-cp)*7 < s)
	return 0;

    do {
	s -= 7;
	*cp++ = ((i>>s) & 0x7f) + (s?128:0);
    } while (s);

    return cp-op;
}

static inline int var_get_u64(uint8_t *cp, const uint8_t *endp, uint64_t *i) {
    uint8_t *op = cp, c;
    uint64_t j = 0;

    if (endp) {
	if (cp >= endp) {
	    *i = 0;
	    return 0;
	}

	do {
	    c = *cp++;
	    j = (j<<7) | (c & 0x7f);
	} while ((c & 0x80) && cp < endp);
    } else {
	// unsafe variant
	do {
	    c = *cp++;
	    j = (j<<7) | (c & 0x7f);
	} while ((c & 0x80));
    }
    *i = j;
    return cp-op;
}

static inline int var_get_u32(uint8_t *cp, const uint8_t *endp, uint32_t *i) {
    uint8_t *op = cp, c;
    uint32_t j = 0;

    if (endp) {
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
    } else {
	// unsafe variant
	do {
	    c = *cp++;
	    j = (j<<7) | (c & 0x7f);
	} while ((c & 0x80));
    }

    *i = j;
    return cp-op;
}
#else

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
#endif

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
