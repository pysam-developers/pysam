// VariantKey
//
// hex.h
//
// @category   Libraries
// @author     Nicola Asuni <info@tecnick.com>
// @link       https://github.com/tecnickcom/variantkey
// @license    MIT [LICENSE](https://raw.githubusercontent.com/tecnickcom/variantkey/main/LICENSE)
// @copyright  2017-2018 GENOMICS plc
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
 * @file hex.h
 * @brief Utility functions to manipulate strings.
 *
 * Collection of utility functions to manipulate strings.
 */

#ifndef ASTRING_H
#define ASTRING_H

#include <inttypes.h>
#include <stdio.h>

/** @brief Returns uint64_t hexadecimal string (16 characters).
 *
 * @param n     Number to parse
 * @param str   String buffer to be returned (it must be sized 17 bytes at least).
 *
 * @return      Upon successful return, these function returns the number of characters processed
 *              (excluding the null byte used to end output to strings).
 *              If the buffer size is not sufficient, then the return value is the number of characters required for
 *              buffer string, including the terminating null byte.
 */
static inline size_t hex_uint64_t(uint64_t n, char *str)
{
    return sprintf(str, "%016" PRIx64, n);
}

/** @brief Parses a 16 chars hexadecimal string and returns the code.
 *
 * @param s    Hexadecimal string to parse (it must contain 16 hexadecimal characters).
 *
 * @return uint64_t unsigned integer number.
 */
static inline uint64_t parse_hex_uint64_t(const char *s)
{
    uint64_t v = 0;
    uint8_t b;
    size_t i;
    for (i = 0; i < 16; i++)
    {
        b = s[i];
        if (b >= 'a')
        {
            b -= ('a' - 10); // a-f
        }
        else
        {
            if (b >= 'A')
            {
                b -= ('A' - 10); // A-F
            }
            else
            {
                b -= '0'; // 0-9
            }
        }
        v = ((v << 4) | b);
    }
    return v;
}

#endif  // ASTRING_H
