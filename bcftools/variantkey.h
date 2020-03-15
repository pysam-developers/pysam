// VariantKey
//
// variantkey.h
//
// @category   Libraries
// @author     Nicola Asuni <nicola.asuni@genomicsplc.com>
// @copyright  2017-2018 GENOMICS plc
// @license    MIT (see LICENSE)
// @link       https://github.com/genomicsplc/variantkey
//
// LICENSE
//
// Copyright (c) 2017-2018 GENOMICS plc
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

/**
 * @file variantkey.h
 * @brief VariantKey main functions.
 *
 * The functions provided here allows to generate and process a 64 bit Unsigned Integer Keys for Human Genetic Variants.
 * The VariantKey is sortable for chromosome and position,
 * and it is also fully reversible for variants with up to 11 bases between Reference and Alternate alleles.
 * It can be used to sort, search and match variant-based data easily and very quickly.
 */

#ifndef VARIANTKEY_H
#define VARIANTKEY_H

#include <inttypes.h>
#include <stddef.h>
#include <stdio.h>
#include "hex.h"

#define VKMASK_CHROM    0xF800000000000000  //!< VariantKey binary mask for CHROM     [ 11111000 00000000 00000000 00000000 00000000 00000000 00000000 00000000 ]
#define VKMASK_POS      0x07FFFFFF80000000  //!< VariantKey binary mask for POS       [ 00000111 11111111 11111111 11111111 10000000 00000000 00000000 00000000 ]
#define VKMASK_CHROMPOS 0xFFFFFFFF80000000  //!< VariantKey binary mask for CHROM+POS [ 11111111 11111111 11111111 11111111 10000000 00000000 00000000 00000000 ]
#define VKMASK_REFALT   0x000000007FFFFFFF  //!< VariantKey binary mask for REF+ALT   [ 00000000 00000000 00000000 00000000 01111111 11111111 11111111 11111111 ]
#define VKSHIFT_CHROM   59 //!< CHROM LSB position from the VariantKey LSB
#define VKSHIFT_POS     31 //!< POS LSB position from the VariantKey LSB

/**
 * VariantKey struct.
 * Contains the numerically encoded VariantKey components (CHROM, POS, REF+ALT).
 */
typedef struct variantkey_t
{
    uint8_t chrom;   //!< Chromosome encoded number (only the LSB 5 bit are used)
    uint32_t pos;    //!< Reference position, with the first base having position 0 (only the LSB 28 bit are used)
    uint32_t refalt; //!< Code for Reference and Alternate allele (only the LSB 31 bits are used)
} variantkey_t;

/**
 * Struct containing the minimum and maximum VariantKey values for range searches.
 */
typedef struct vkrange_t
{
    uint64_t min; //!< Minimum VariantKey value for any given REF+ALT encoding
    uint64_t max; //!< Maximum VariantKey value for any given REF+ALT encoding
} vkrange_t;

/** @brief Returns chromosome numerical encoding.
 *
 * @param chrom  Chromosome. An identifier from the reference genome, no white-space permitted.
 * @param size   Length of the chrom string, excluding the terminating null byte.
 *
 * @return CHROM code
 */
static inline uint8_t encode_chrom(const char *chrom, size_t size)
{
    // X > 23 ; Y > 24 ; M > 25
    static const uint8_t onecharmap[] =
    {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        /*                                    M                                X  Y                  */
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,23,24, 0, 0, 0, 0, 0, 0,
        /*                                    m                                x  y                  */
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,23,24, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };
    // remove "chr" prefix
    if ((size > 3)
            && ((chrom[0] == 'c') || (chrom[0] == 'C'))
            && ((chrom[1] == 'h') || (chrom[1] == 'H'))
            && ((chrom[2] == 'r') || (chrom[2] == 'R')))
    {
        chrom += 3;
        size -= 3;
    }
    if (size == 0)
    {
        return 0;
    }
    if ((chrom[0] <= '9') && (chrom[0] >= '0')) // Number
    {
        size_t i;
        uint8_t v = (chrom[0] - '0');
        for (i = 1; i < size; i++)
        {
            if ((chrom[i] > '9') || (chrom[i] < '0'))
            {
                return 0; // NA
            }
            v = ((v * 10) + (chrom[i] - '0'));
        }
        return v;
    }
    if ((size == 1) || ((size == 2) && ((chrom[1] == 'T') || (chrom[1] == 't'))))
    {
        return onecharmap[((uint8_t)chrom[0])];
    }
    return 0; // NA
}

/** @brief Decode the chromosome numerical code.
 *
 * @param code   CHROM code.
 * @param chrom  CHROM string buffer to be returned. Its size should be enough to contain the results (max 4 bytes).
 *
 * @return If successful, the total number of characters written is returned,
 *         excluding the null-character appended at the end of the string,
 *         otherwise a negative number is returned in case of failure.
 */
static inline size_t decode_chrom(uint8_t code, char *chrom)
{
    if ((code < 1) || (code > 25))
    {
        return sprintf(chrom, "NA");
    }
    if (code < 23)
    {
        return sprintf(chrom, "%" PRIu8, code);
    }
    static const char *map[] = {"X", "Y", "MT"};
    return sprintf(chrom, "%s", map[(code - 23)]);
}

static inline uint32_t encode_base(const uint8_t c)
{
    /*
      Encode base:
      A > 0
      C > 1
      G > 2
      T > 3
    */
    static const uint32_t map[] =
    {
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        /*A   C       G                         T*/
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
        /*a   c       g                         t*/
        4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
        4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    };
    return map[c];
}

static inline int encode_allele(uint32_t *h, uint8_t *bitpos, const char *str, size_t size)
{
    uint32_t v;
    while (size--)
    {
        v = encode_base(*str++);
        if (v > 3)
        {
            return -1;
        }
        *bitpos -= 2;
        *h |= (v << *bitpos);
    }
    return 0;
}

static inline uint32_t encode_refalt_rev(const char *ref, size_t sizeref, const char *alt, size_t sizealt)
{
    //[******** ******** ******** ******** *RRRRAAA A1122334 45566778 8990011*]
    uint32_t h = 0;
    h |= ((uint32_t)(sizeref) << 27); // RRRR: length of (REF - 1)
    h |= ((uint32_t)(sizealt) << 23); // AAAA: length of (ALT - 1)
    uint8_t bitpos = 23;
    if ((encode_allele(&h, &bitpos, ref, sizeref) < 0) || (encode_allele(&h, &bitpos, alt, sizealt) < 0))
    {
        return 0; // error code
    }
    return h;
}

// Mix two 32 bit hash numbers using a MurmurHash3-like algorithm
static inline uint32_t muxhash(uint32_t k, uint32_t h)
{
    k *= 0xcc9e2d51;
    k = (k >> 17) | (k << 15);
    k *= 0x1b873593;
    h ^= k;
    h = (h >> 19) | (h << 13);
    return ((h * 5) + 0xe6546b64);
}

static inline uint32_t encode_packchar(int c)
{
    if (c < 'A')
    {
        return 27;
    }
    if (c >= 'a')
    {
        return (uint32_t)(c - 'a' + 1);
    }
    return (uint32_t)(c - 'A' + 1);
}

// pack blocks of 6 characters in 32 bit (6 x 5 bit + 2 spare bit) [ 01111122 22233333 44444555 55666660 ]
static inline uint32_t pack_chars_tail(const char *str, size_t size)
{
    uint32_t h = 0;
    const char *pos = (str + size - 1);
    switch (size)
    {
    case 5:
        h ^= encode_packchar(*pos--) << (1 + (5 * 1));
    // fall through
    case 4:
        h ^= encode_packchar(*pos--) << (1 + (5 * 2));
    // fall through
    case 3:
        h ^= encode_packchar(*pos--) << (1 + (5 * 3));
    // fall through
    case 2:
        h ^= encode_packchar(*pos--) << (1 + (5 * 4));
    // fall through
    case 1:
        h ^= encode_packchar(*pos) << (1 + (5 * 5));
    }
    return h;
}

static inline uint32_t pack_chars(const char *str)
{
    const char *pos = (str + 5);
    return ((encode_packchar(*pos) << 1)
            ^ (encode_packchar(*(pos-1)) << (1 + (5 * 1)))
            ^ (encode_packchar(*(pos-2)) << (1 + (5 * 2)))
            ^ (encode_packchar(*(pos-3)) << (1 + (5 * 3)))
            ^ (encode_packchar(*(pos-4)) << (1 + (5 * 4)))
            ^ (encode_packchar(*(pos-5)) << (1 + (5 * 5))));
}

// Return a 32 bit hash of a nucleotide string
static inline uint32_t hash32(const char *str, size_t size)
{
    uint32_t h = 0;
    size_t len = 6;
    while (size >= len)
    {
        h = muxhash(pack_chars(str), h);
        str += len;
        size -= len;
    }
    if (size > 0)
    {
        h = muxhash(pack_chars_tail(str, size), h);
    }
    return h;
}

static inline uint32_t encode_refalt_hash(const char *ref, size_t sizeref, const char *alt, size_t sizealt)
{
    // 0x3 is the separator character between REF and ALT [00000000 00000000 00000000 00000011]
    uint32_t h = muxhash(hash32(alt, sizealt), muxhash(0x3, hash32(ref, sizeref)));
    // MurmurHash3 finalization mix - force all bits of a hash block to avalanche
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return ((h >> 1) | 0x1); // 0x1 is the set bit to indicate HASH mode [00000000 00000000 00000000 00000001]
}

/** @brief Returns reference+alternate numerical encoding.
 *
 * @param ref      Reference allele. String containing a sequence of nucleotide letters.
 *                 The value in the pos field refers to the position of the first nucleotide in the String.
 *                 Characters must be A-Z, a-z or *
 * @param sizeref  Length of the ref string, excluding the terminating null byte.
 * @param alt      Alternate non-reference allele string.
 *                 Characters must be A-Z, a-z or *
 * @param sizealt  Length of the alt string, excluding the terminating null byte.
 *
 * @return REF+ALT code
 */
static inline uint32_t encode_refalt(const char *ref, size_t sizeref, const char *alt, size_t sizealt)
{
    if ((sizeref + sizealt) <= 11)
    {
        uint32_t h = encode_refalt_rev(ref, sizeref, alt, sizealt);
        if (h != 0)
        {
            return h;
        }
    }
    return encode_refalt_hash(ref, sizeref, alt, sizealt);
}

static inline char decode_base(uint32_t code, int bitpos)
{
    static const char base[4] = {'A', 'C', 'G', 'T'};
    return base[((code >> bitpos) & 0x3)]; // 0x3 is the 2 bit mask [00000011]
}

static inline size_t decode_refalt_rev(uint32_t code, char *ref, size_t *sizeref, char *alt, size_t *sizealt)
{
    *sizeref = (size_t)((code & 0x78000000) >> 27); // [01111000 00000000 00000000 00000000]
    *sizealt = (size_t)((code & 0x07800000) >> 23); // [00000111 10000000 00000000 00000000]
    switch (*sizeref)
    {
    case 10:
        ref[9] = decode_base(code, (3 + (2 * 0)));
    // fall through
    case 9:
        ref[8] = decode_base(code, (3 + (2 * 1)));
    // fall through
    case 8:
        ref[7] = decode_base(code, (3 + (2 * 2)));
    // fall through
    case 7:
        ref[6] = decode_base(code, (3 + (2 * 3)));
    // fall through
    case 6:
        ref[5] = decode_base(code, (3 + (2 * 4)));
    // fall through
    case 5:
        ref[4] = decode_base(code, (3 + (2 * 5)));
    // fall through
    case 4:
        ref[3] = decode_base(code, (3 + (2 * 6)));
    // fall through
    case 3:
        ref[2] = decode_base(code, (3 + (2 * 7)));
    // fall through
    case 2:
        ref[1] = decode_base(code, (3 + (2 * 8)));
    // fall through
    case 1:
        ref[0] = decode_base(code, (3 + (2 * 9)));
    }
    ref[*sizeref] = 0;
    uint8_t bitpos = (23 - ((*sizeref) << 1));
    switch (*sizealt)
    {
    case 10:
        alt[9] = decode_base(code, bitpos - (2 * 10));
    // fall through
    case 9:
        alt[8] = decode_base(code, bitpos - (2 * 9));
    // fall through
    case 8:
        alt[7] = decode_base(code, bitpos - (2 * 8));
    // fall through
    case 7:
        alt[6] = decode_base(code, bitpos - (2 * 7));
    // fall through
    case 6:
        alt[5] = decode_base(code, bitpos - (2 * 6));
    // fall through
    case 5:
        alt[4] = decode_base(code, bitpos - (2 * 5));
    // fall through
    case 4:
        alt[3] = decode_base(code, bitpos - (2 * 4));
    // fall through
    case 3:
        alt[2] = decode_base(code, bitpos - (2 * 3));
    // fall through
    case 2:
        alt[1] = decode_base(code, bitpos - (2 * 2));
    // fall through
    case 1:
        alt[0] = decode_base(code, bitpos - (2 * 1));
    }
    alt[*sizealt] = 0;
    return (*sizeref + *sizealt);
}

/** @brief Decode the 32 bit REF+ALT code if reversible (if it has 11 or less bases in total and only contains ACGT letters).
 *
 * @param code     REF+ALT code
 * @param ref      REF string buffer to be returned.
 * @param sizeref  Pointer to the size of the ref buffer, excluding the terminating null byte.
 *                 This will contain the final ref size.
 * @param alt      ALT string buffer to be returned.
 * @param sizealt  Pointer to the size of the alt buffer, excluding the terminating null byte.
 *                 This will contain the final alt size.
 *
 * @return      If the code is reversible, then the total number of characters of REF+ALT is returned.
 *              Otherwise 0 is returned.
 */
static inline size_t decode_refalt(uint32_t code, char *ref, size_t *sizeref, char *alt, size_t *sizealt)
{
    if (code & 0x1) // check last bit
    {
        return 0; // non-reversible encoding
    }
    return decode_refalt_rev(code, ref, sizeref, alt, sizealt);
}

/** @brief Returns a 64 bit variant key based on the pre-encoded CHROM, POS (0-based) and REF+ALT.
 *
 * @param chrom      Encoded Chromosome (see encode_chrom).
 * @param pos        Position. The reference position, with the first base having position 0.
 * @param refalt     Encoded Reference + Alternate (see encode_refalt).
 *
 * @return      VariantKey 64 bit code.
 */
static inline uint64_t encode_variantkey(uint8_t chrom, uint32_t pos, uint32_t refalt)
{
    return (((uint64_t)chrom << VKSHIFT_CHROM) | ((uint64_t)pos << VKSHIFT_POS) | (uint64_t)refalt);
}

/** @brief Extract the CHROM code from VariantKey.
 *
 * @param vk VariantKey code.
 *
 * @return CHROM code.
 */
static inline uint8_t extract_variantkey_chrom(uint64_t vk)
{
    return (uint8_t)((vk & VKMASK_CHROM) >> VKSHIFT_CHROM);
}

/** @brief Extract the POS code from VariantKey.
 *
 * @param vk VariantKey code.
 *
 * @return POS.
 */
static inline uint32_t extract_variantkey_pos(uint64_t vk)
{
    return (uint32_t)((vk & VKMASK_POS) >> VKSHIFT_POS);
}

/** @brief Extract the REF+ALT code from VariantKey.
 *
 * @param vk VariantKey code.
 *
 * @return REF+ALT code.
 */
static inline uint32_t extract_variantkey_refalt(uint64_t vk)
{
    return (uint32_t)(vk & VKMASK_REFALT);
}

/** @brief Decode a VariantKey code and returns the components as variantkey_t structure.
 *
 * @param code VariantKey code.
 * @param vk   Decoded variantkey structure.
 */
static inline void decode_variantkey(uint64_t code, variantkey_t *vk)
{
    vk->chrom = extract_variantkey_chrom(code);
    vk->pos = extract_variantkey_pos(code);
    vk->refalt = extract_variantkey_refalt(code);
}

/** @brief Returns a 64 bit variant key based on CHROM, POS (0-based), REF, ALT.
 *
 * @param chrom      Chromosome. An identifier from the reference genome, no white-space or leading zeros permitted.
 * @param sizechrom  Length of the chrom string, excluding the terminating null byte.
 * @param pos        Position. The reference position, with the first base having position 0.
 * @param ref        Reference allele. String containing a sequence of nucleotide letters.
 *                   The value in the pos field refers to the position of the first nucleotide in the String.
 *                   Characters must be A-Z, a-z or *
 * @param sizeref    Length of the ref string, excluding the terminating null byte.
 * @param alt        Alternate non-reference allele string.
 *                   Characters must be A-Z, a-z or *
 * @param sizealt    Length of the alt string, excluding the terminating null byte.
 *
 * @return      VariantKey 64 bit code.
 */
static inline uint64_t variantkey(const char *chrom, size_t sizechrom, uint32_t pos, const char *ref, size_t sizeref, const char *alt, size_t sizealt)
{
    return encode_variantkey(encode_chrom(chrom, sizechrom), pos, encode_refalt(ref, sizeref, alt, sizealt));
}

/** @brief Returns minimum and maximum VariantKeys for range searches.
 *
 * @param chrom     Chromosome encoded number.
 * @param pos_min   Start reference position, with the first base having position 0.
 * @param pos_max   End reference position, with the first base having position 0.
 * @param range     VariantKey range values.
 */
static inline void variantkey_range(uint8_t chrom, uint32_t pos_min, uint32_t pos_max, vkrange_t *range)
{
    uint64_t c = ((uint64_t)chrom << VKSHIFT_CHROM);
    range->min = (c | ((uint64_t)pos_min << VKSHIFT_POS));
    range->max = (c | ((uint64_t)pos_max << VKSHIFT_POS) | VKMASK_REFALT);
}

static inline int8_t compare_uint64_t(uint64_t a, uint64_t b)
{
    return (a < b) ? -1 : (a > b);
}

/** @brief Compares two VariantKeys by chromosome only.
 *
 * @param vka    The first VariantKey to be compared.
 * @param vkb    The second VariantKey to be compared.
 *
 * @return -1 if the first chromosome is smaller than the second, 0 if they are equal and 1 if the first is greater than the second.
 */
static inline int8_t compare_variantkey_chrom(uint64_t vka, uint64_t vkb)
{
    return compare_uint64_t((vka >> VKSHIFT_CHROM), (vkb >> VKSHIFT_CHROM));
}

/** @brief Compares two VariantKeys by chromosome and position.
 *
 * @param vka    The first VariantKey to be compared.
 * @param vkb    The second VariantKey to be compared.
 *
 * @return -1 if the first CHROM+POS is smaller than the second, 0 if they are equal and 1 if the first is greater than the second.
 */
static inline int8_t compare_variantkey_chrom_pos(uint64_t vka, uint64_t vkb)
{
    return compare_uint64_t((vka >> VKSHIFT_POS), (vkb >> VKSHIFT_POS));
}

/** @brief Returns VariantKey hexadecimal string (16 characters).
 *
 * The string represent a 64 bit number or:
 *   -  5 bit for CHROM
 *   - 28 bit for POS
 *   - 31 bit for REF+ALT
 *
 * @param vk    VariantKey code.
 * @param str   String buffer to be returned (it must be sized 17 bytes at least).
 *
 * @return      Upon successful return, these function returns the number of characters processed
 *              (excluding the null byte used to end output to strings).
 *              If the buffer size is not sufficient, then the return value is the number of characters required for
 *              buffer string, including the terminating null byte.
 */
static inline size_t variantkey_hex(uint64_t vk, char *str)
{
    return hex_uint64_t(vk, str);
}

/** @brief Parses a VariantKey hexadecimal string and returns the code.
 *
 * @param vs    VariantKey hexadecimal string (it must contain 16 hexadecimal characters).
 *
 * @return A VariantKey code.
 */
static inline uint64_t parse_variantkey_hex(const char *vs)
{
    return parse_hex_uint64_t(vs);
}

#endif  // VARIANTKEY_H
