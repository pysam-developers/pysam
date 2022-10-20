/* rans_byte.h originally from https://github.com/rygorous/ryg_rans
 *
 * This is a public-domain implementation of several rANS variants. rANS is an
 * entropy coder from the ANS family, as described in Jarek Duda's paper
 * "Asymmetric numeral systems" (http://arxiv.org/abs/1311.2540).
 */

/*-------------------------------------------------------------------------- */
/* rans_byte.h from https://github.com/rygorous/ryg_rans */

// Simple byte-aligned rANS encoder/decoder - public domain - Fabian 'ryg' Giesen 2014
//
// Not intended to be "industrial strength"; just meant to illustrate the general
// idea.

#ifndef RANS_BYTE_HEADER
#define RANS_BYTE_HEADER

#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include "utils.h"

#ifdef assert
#define RansAssert assert
#else
#define RansAssert(x)
#endif

// READ ME FIRST:
//
// This is designed like a typical arithmetic coder API, but there's three
// twists you absolutely should be aware of before you start hacking:
//
// 1. You need to encode data in *reverse* - last symbol first. rANS works
//    like a stack: last in, first out.
// 2. Likewise, the encoder outputs bytes *in reverse* - that is, you give
//    it a pointer to the *end* of your buffer (exclusive), and it will
//    slowly move towards the beginning as more bytes are emitted.
// 3. Unlike basically any other entropy coder implementation you might
//    have used, you can interleave data from multiple independent rANS
//    encoders into the same bytestream without any extra signaling;
//    you can also just write some bytes by yourself in the middle if
//    you want to. This is in addition to the usual arithmetic encoder
//    property of being able to switch models on the fly. Writing raw
//    bytes can be useful when you have some data that you know is
//    incompressible, and is cheaper than going through the rANS encode
//    function. Using multiple rANS coders on the same byte stream wastes
//    a few bytes compared to using just one, but execution of two
//    independent encoders can happen in parallel on superscalar and
//    Out-of-Order CPUs, so this can be *much* faster in tight decoding
//    loops.
//
//    This is why all the rANS functions take the write pointer as an
//    argument instead of just storing it in some context struct.

// --------------------------------------------------------------------------

// L ('l' in the paper) is the lower bound of our normalization interval.
// Between this and our byte-aligned emission, we use 31 (not 32!) bits.
// This is done intentionally because exact reciprocals for 31-bit uints
// fit in 32-bit uints: this permits some optimizations during encoding.
#define RANS_BYTE_L (1u << 23)  // lower bound of our normalization interval

// State for a rANS encoder. Yep, that's all there is to it.
typedef uint32_t RansState;

// Initialize a rANS encoder.
static inline void RansEncInit(RansState* r)
{
    *r = RANS_BYTE_L;
}

#if 0 /* Curently unused */
// Renormalize the encoder. Internal function.
static inline RansState RansEncRenorm(RansState x, uint8_t** pptr, uint32_t freq, uint32_t scale_bits)
{
    uint32_t x_max = ((RANS_BYTE_L >> scale_bits) << 8) * freq; // this turns into a shift.
    if (x >= x_max) {
        uint8_t* ptr = *pptr;
        do {
            *--ptr = (uint8_t) (x & 0xff);
            x >>= 8;
        } while (x >= x_max);
        *pptr = ptr;
    }
    return x;
}

// Encodes a single symbol with range start "start" and frequency "freq".
// All frequencies are assumed to sum to "1 << scale_bits", and the
// resulting bytes get written to ptr (which is updated).
//
// NOTE: With rANS, you need to encode symbols in *reverse order*, i.e. from
// beginning to end! Likewise, the output bytestream is written *backwards*:
// ptr starts pointing at the end of the output buffer and keeps decrementing.
static inline void RansEncPut(RansState* r, uint8_t** pptr, uint32_t start, uint32_t freq, uint32_t scale_bits)
{
    // renormalize
    RansState x = RansEncRenorm(*r, pptr, freq, scale_bits);

    // x = C(s,x)
    *r = ((x / freq) << scale_bits) + (x % freq) + start;
}
#endif /* Curently unused */

// Flushes the rANS encoder.
static inline void RansEncFlush(RansState* r, uint8_t** pptr)
{
    uint32_t x = *r;
    uint8_t* ptr = *pptr;

    ptr -= 4;
    ptr[0] = (uint8_t) (x >> 0);
    ptr[1] = (uint8_t) (x >> 8);
    ptr[2] = (uint8_t) (x >> 16);
    ptr[3] = (uint8_t) (x >> 24);

    *pptr = ptr;
}

// Initializes a rANS decoder.
// Unlike the encoder, the decoder works forwards as you'd expect.
static inline void RansDecInit(RansState* r, uint8_t** pptr)
{
    uint32_t x;
    uint8_t* ptr = *pptr;

    x  = ptr[0] << 0;
    x |= ptr[1] << 8;
    x |= ptr[2] << 16;
    x |= ((uint32_t)ptr[3]) << 24;
    ptr += 4;

    *pptr = ptr;
    *r = x;
}

// Returns the current cumulative frequency (map it to a symbol yourself!)
static inline uint32_t RansDecGet(RansState* r, uint32_t scale_bits)
{
    return *r & ((1u << scale_bits) - 1);
}

// Advances in the bit stream by "popping" a single symbol with range start
// "start" and frequency "freq". All frequencies are assumed to sum to "1 << scale_bits",
// and the resulting bytes get written to ptr (which is updated).
static inline void RansDecAdvance(RansState* r, uint8_t** pptr, uint32_t start, uint32_t freq, uint32_t scale_bits)
{
    uint32_t mask = (1u << scale_bits) - 1;

    // s, x = D(x)
    uint32_t x = *r;
    x = freq * (x >> scale_bits) + (x & mask) - start;

    // renormalize
    if (x < RANS_BYTE_L) {
        uint8_t* ptr = *pptr;
        do x = (x << 8) | *ptr++; while (x < RANS_BYTE_L);
        *pptr = ptr;
    }

    *r = x;
}

// --------------------------------------------------------------------------

// That's all you need for a full encoder; below here are some utility
// functions with extra convenience or optimizations.

// Encoder symbol description
// This (admittedly odd) selection of parameters was chosen to make
// RansEncPutSymbol as cheap as possible.
typedef struct {
    uint32_t x_max;     // (Exclusive) upper bound of pre-normalization interval
    uint32_t rcp_freq;  // Fixed-point reciprocal frequency
    uint32_t bias;      // Bias
    uint16_t cmpl_freq; // Complement of frequency: (1 << scale_bits) - freq
    uint16_t rcp_shift; // Reciprocal shift
} RansEncSymbol;

// Decoder symbols are straightforward.
// 32-bit means more memory, but oddly faster on old gcc? Why?
// 322MB/s vs 309MB/s for order-1.
typedef struct {
    uint16_t freq;      // Symbol frequency.
    uint16_t start;     // Start of range.
} RansDecSymbol;

typedef struct {
    uint32_t freq;      // Symbol frequency.
    uint32_t start;     // Start of range.
} RansDecSymbol32;

// Initializes an encoder symbol to start "start" and frequency "freq"
static inline void RansEncSymbolInit(RansEncSymbol* s, uint32_t start, uint32_t freq, uint32_t scale_bits)
{
    RansAssert(scale_bits <= 16);
    RansAssert(start <= (1u << scale_bits));
    RansAssert(freq <= (1u << scale_bits) - start);

    // Say M := 1 << scale_bits.
    //
    // The original encoder does:
    //   x_new = (x/freq)*M + start + (x%freq)
    //
    // The fast encoder does (schematically):
    //   q     = mul_hi(x, rcp_freq) >> rcp_shift   (division)
    //   r     = x - q*freq                         (remainder)
    //   x_new = q*M + bias + r                     (new x)
    // plugging in r into x_new yields:
    //   x_new = bias + x + q*(M - freq)
    //        =: bias + x + q*cmpl_freq             (*)
    //
    // and we can just precompute cmpl_freq. Now we just need to
    // set up our parameters such that the original encoder and
    // the fast encoder agree.

    s->x_max = ((RANS_BYTE_L >> scale_bits) << 8) * freq;
    s->cmpl_freq = (uint16_t) ((1 << scale_bits) - freq);
    if (freq < 2) {
        // freq=0 symbols are never valid to encode, so it doesn't matter what
        // we set our values to.
        //
        // freq=1 is tricky, since the reciprocal of 1 is 1; unfortunately,
        // our fixed-point reciprocal approximation can only multiply by values
        // smaller than 1.
        //
        // So we use the "next best thing": rcp_freq=0xffffffff, rcp_shift=0.
        // This gives:
        //   q = mul_hi(x, rcp_freq) >> rcp_shift
        //     = mul_hi(x, (1<<32) - 1)) >> 0
        //     = floor(x - x/(2^32))
        //     = x - 1 if 1 <= x < 2^32
        // and we know that x>0 (x=0 is never in a valid normalization interval).
        //
        // So we now need to choose the other parameters such that
        //   x_new = x*M + start
        // plug it in:
        //     x*M + start                   (desired result)
        //   = bias + x + q*cmpl_freq        (*)
        //   = bias + x + (x - 1)*(M - 1)    (plug in q=x-1, cmpl_freq)
        //   = bias + 1 + (x - 1)*M
        //   = x*M + (bias + 1 - M)
        //
        // so we have start = bias + 1 - M, or equivalently
        //   bias = start + M - 1.
        s->rcp_freq = ~0u;
        s->rcp_shift = 0;
        s->bias = start + (1 << scale_bits) - 1;
    } else {
        // Alverson, "Integer Division using reciprocals"
        // shift=ceil(log2(freq))
        uint32_t shift = 0;
        while (freq > (1u << shift))
            shift++;

        s->rcp_freq = (uint32_t) (((1ull << (shift + 31)) + freq-1) / freq);
        s->rcp_shift = shift - 1;

        // With these values, 'q' is the correct quotient, so we
        // have bias=start.
        s->bias = start;
    }

    s->rcp_shift += 32; // Avoid the extra >>32 in RansEncPutSymbol
}

// Initialize a decoder symbol to start "start" and frequency "freq"
static inline void RansDecSymbolInit(RansDecSymbol* s, uint32_t start, uint32_t freq)
{
    RansAssert(start <= (1 << 16));
    RansAssert(freq <= (1 << 16) - start);
    s->start = (uint16_t) start;
    s->freq = (uint16_t) freq;
}

// Encodes a given symbol. This is faster than straight RansEnc since we can do
// multiplications instead of a divide.
//
// See RansEncSymbolInit for a description of how this works.
static inline void RansEncPutSymbol(RansState* r, uint8_t** pptr, RansEncSymbol const* sym)
{
    RansAssert(sym->x_max != 0); // can't encode symbol with freq=0

    // renormalize
    uint32_t x = *r;
    uint32_t x_max = sym->x_max;

    // This is better for 40-qual illumina (3.7% quicker overall CRAM).
    // The old method was better for low complexity data such as NovaSeq
    // quals (2.6% quicker overall CRAM).
    int o = x >= x_max;
    uint8_t* ptr = *pptr;
    ptr[-1] = x & 0xff;
    ptr -= o;
    x >>= o*8;

    if (unlikely(x >= x_max)) {
        *--ptr = (uint8_t) (x & 0xff);
        x >>= 8;
    }
    *pptr = ptr;

    //uint32_t q = (uint32_t) (((uint64_t)x * sym->rcp_freq) >> sym->rcp_shift);
    //*r = q * sym->cmpl_freq + x + sym->bias;

    // x = C(s,x)
    // NOTE: written this way so we get a 32-bit "multiply high" when
    // available. If you're on a 64-bit platform with cheap multiplies
    // (e.g. x64), just bake the +32 into rcp_shift.
    //uint32_t q = (uint32_t) (((uint64_t)x * sym->rcp_freq) >> 32) >> sym->rcp_shift;

    // The extra >>32 has already been added to RansEncSymbolInit
    uint32_t q = (uint32_t) (((uint64_t)x * sym->rcp_freq) >> sym->rcp_shift);
    *r = q * sym->cmpl_freq + x + sym->bias;
}

// A 4-way version of RansEncPutSymbol, renormalising 4 states
// simulatenously with their results written to the same ptr buffer.
// (This is perhaps a failing as it makes optmisation tricky.)
static inline void RansEncPutSymbol4(RansState *r0,
                                     RansState *r1,
                                     RansState *r2,
                                     RansState *r3,
                                     uint8_t** pptr,
                                     RansEncSymbol const *sym0,
                                     RansEncSymbol const *sym1,
                                     RansEncSymbol const *sym2,
                                     RansEncSymbol const *sym3)
{
    RansAssert(sym0->x_max != 0); // can't encode symbol with freq=0
    RansAssert(sym1->x_max != 0); // can't encode symbol with freq=0
    RansAssert(sym2->x_max != 0); // can't encode symbol with freq=0
    RansAssert(sym3->x_max != 0); // can't encode symbol with freq=0

    // renormalize
    uint32_t x0, x1, x2, x3;
    uint8_t* ptr = *pptr;

    int o;
    uint32_t m[4] = {
        sym0->x_max,
        sym1->x_max,
        sym2->x_max,
        sym3->x_max
    };

    x0 = *r0;
    o = x0 >= m[0];
    ptr[-1] = x0;
    ptr -= o;
    x0 >>= o*8;
    if (x0 >= m[0]) {
        *--ptr = x0;
        x0 >>= 8;
    }

    x1 = *r1;
    o = x1 >= m[1];
    ptr[-1] = x1;
    ptr -= o;
    x1 >>= o*8;
    if (x1 >= m[1]) {
        *--ptr = x1;
        x1 >>= 8;
    }

    x2 = *r2;
    o = x2 >= m[2];
    ptr[-1] = x2;
    ptr -= o;
    x2 >>= o*8;
    if (x2 >= m[2]) {
        *--ptr = x2;
        x2 >>= 8;
    }

    x3 = *r3;
    o = x3 >= m[3];
    ptr[-1] = x3;
    ptr -= o;
    x3 >>= o*8;
    if (x3 >= m[3]) {
        *--ptr = x3;
        x3 >>= 8;
    }

    *pptr = ptr;

    // x = C(s,x)
    uint32_t qa, qb;
    qa = (uint32_t) (((uint64_t)x0 * sym0->rcp_freq) >> sym0->rcp_shift);
    uint32_t X0 = qa * sym0->cmpl_freq;
    qb = (uint32_t) (((uint64_t)x1 * sym1->rcp_freq) >> sym1->rcp_shift);
    uint32_t X1 = qb * sym1->cmpl_freq;

    *r0 = X0 + x0 + sym0->bias;
    *r1 = X1 + x1 + sym1->bias;

    qa = (uint32_t) (((uint64_t)x2 * sym2->rcp_freq) >> sym2->rcp_shift);
    uint32_t X2 = qa * sym2->cmpl_freq;
    qb = (uint32_t) (((uint64_t)x3 * sym3->rcp_freq) >> sym3->rcp_shift);
    uint32_t X3 = qb * sym3->cmpl_freq;

    *r2 = X2 + x2 + sym2->bias;
    *r3 = X3 + x3 + sym3->bias;
}

// Equivalent to RansDecAdvance that takes a symbol.
static inline void RansDecAdvanceSymbol(RansState* r, uint8_t** pptr, RansDecSymbol const* sym, uint32_t scale_bits)
{
    RansDecAdvance(r, pptr, sym->start, sym->freq, scale_bits);
}

// Advances in the bit stream by "popping" a single symbol with range start
// "start" and frequency "freq". All frequencies are assumed to sum to "1 << scale_bits".
// No renormalization or output happens.
static inline void RansDecAdvanceStep(RansState* r, uint32_t start, uint32_t freq, uint32_t scale_bits)
{
    uint32_t mask = (1u << scale_bits) - 1;

    // s, x = D(x)
    uint32_t x = *r;
    *r = freq * (x >> scale_bits) + (x & mask) - start;
}

// Equivalent to RansDecAdvanceStep that takes a symbol.
static inline void RansDecAdvanceSymbolStep(RansState* r, RansDecSymbol const* sym, uint32_t scale_bits)
{
    RansDecAdvanceStep(r, sym->start, sym->freq, scale_bits);
}

// Renormalize.
#if defined(__x86_64) && !defined(__ILP32__)
/*
 * Assembly variants of the RansDecRenorm code.
 * These are based on joint ideas from Rob Davies and from looking at
 * the clang assembly output.
 */
static inline void RansDecRenorm(RansState* r, uint8_t** pptr) {
    uint32_t  x   = *r;
    uint8_t  *ptr = *pptr;

    __asm__ ("movzbl (%0), %%eax\n\t"
             "mov    %1, %%edx\n\t"
             "shl    $0x8,%%edx\n\t"
             "or     %%eax,%%edx\n\t"
             "cmp    $0x800000,%1\n\t"
             "cmovb  %%edx,%1\n\t"
             "adc    $0x0,%0\n\t"
             : "=r" (ptr), "=r" (x)
             : "0" (ptr), "1" (x)
             : "eax", "edx"
             );
    if (x < 0x800000) x = (x << 8) | *ptr++;
    *pptr = ptr;
    *r = x;
}

/*
 * A variant that normalises two rans states.
 * The only minor tweak here is to adjust the reorder a few opcodes
 * to reduce dependency delays.
 */
static inline void RansDecRenorm2(RansState* r1, RansState* r2, uint8_t** pptr) {
    uint32_t  x1   = *r1;
    uint32_t  x2   = *r2;
    uint8_t  *ptr = *pptr;

    __asm__ ("movzbl (%0), %%eax\n\t"
             "mov    %1, %%edx\n\t"
             "shl    $0x8, %%edx\n\t"
             "or     %%eax, %%edx\n\t"
             "cmp    $0x800000, %1\n\t"
             "cmovb  %%edx, %1\n\t"
             "adc    $0x0, %0\n\t"
             "mov    %2, %%edx\n\t"
             "shl    $0x8, %%edx\n\t"
             "cmp    $0x800000, %1\n\t"
             "jae    1f\n\t"
             "movzbl (%0), %%eax\n\t"
             "shl    $0x8, %1\n\t"
             "or     %%eax, %1\n\t"
             "add    $0x1, %0\n\t"
             "1:\n\t"
             "movzbl (%0), %%eax\n\t"
             "or     %%eax, %%edx\n\t"
             "cmp    $0x800000, %2\n\t"
             "cmovb  %%edx, %2\n\t"
             "adc    $0x0, %0\n\t"
             "cmp    $0x800000, %2\n\t"
             "jae    2f\n\t"
             "movzbl (%0), %%eax\n\t"
             "shl    $0x8, %2\n\t"
             "or     %%eax, %2\n\t"
             "add    $0x1, %0\n\t"
             "2:\n\t"
             : "=r" (ptr), "=r" (x1), "=r" (x2)
             : "0" (ptr), "1" (x1), "2" (x2)
             : "eax", "edx"
             );

    *pptr = ptr;
    *r1 = x1;
    *r2 = x2;
}

#else /* __x86_64 */

static inline void RansDecRenorm(RansState* r, uint8_t** pptr)
{
    // renormalize
    uint32_t x = *r;

#ifdef __clang__
    // Generates cmov instructions on clang, but alas not gcc
    uint8_t* ptr = *pptr;
    uint32_t y = (x << 8) | *ptr;
    uint32_t cond = x < RANS_BYTE_L;
    x    = cond ? y : x;
    ptr += cond ? 1 : 0;
    if (x < RANS_BYTE_L) x = (x<<8) | *ptr++;
    *pptr = ptr;
#else
    if (x >= RANS_BYTE_L) return;
    uint8_t* ptr = *pptr;
    x = (x << 8) | *ptr++;
    if (x < RANS_BYTE_L) x = (x << 8) | *ptr++;
    *pptr = ptr;
#endif /* __clang__ */

    *r = x;
}

static inline void RansDecRenorm2(RansState* r1, RansState* r2, uint8_t** pptr) {
    RansDecRenorm(r1, pptr);
    RansDecRenorm(r2, pptr);
}

#endif /* __x86_64 */

static inline void RansDecRenormSafe(RansState* r, uint8_t** pptr, uint8_t *ptr_end)
{
    uint32_t x = *r;
    uint8_t* ptr = *pptr;
    if (x >= RANS_BYTE_L || ptr >= ptr_end) return;
    x = (x << 8) | *ptr++;
    if (x < RANS_BYTE_L && ptr < ptr_end)
        x = (x << 8) | *ptr++;
    *pptr = ptr;
    *r = x;
}

static inline void RansDecSymbolInit32(RansDecSymbol32* s, uint32_t start, uint32_t freq)
{
    RansAssert(start <= (1 << 16));
    RansAssert(freq <= (1 << 16) - start);
    s->start = (uint16_t) start;
    s->freq = (uint16_t) freq;
}

static inline void RansDecAdvanceSymbol32(RansState* r, uint8_t** pptr, RansDecSymbol32 const* sym, uint32_t scale_bits)
{
    RansDecAdvance(r, pptr, sym->start, sym->freq, scale_bits);
}

#endif // RANS_BYTE_HEADER
