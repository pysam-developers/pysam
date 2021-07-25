/*  str_finder.c -- Short Tandem Repeat finder.
    Originally from Crumble (https://github.com/jkbonfield/crumble)

    Copyright (C) 2015-2016, 2021 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef _STR_FINDER_H_
#define _STR_FINDER_H_

#include "utlist.h"

typedef struct rep_ele {
    int start, end, rep_len;
    struct rep_ele *prev;
    struct rep_ele *next;
} rep_ele;

/*
 * Finds repeated homopolymers up to 8-mers.
 *
 * If lower_only is true then it only adds STRs for regions that
 * contain at least one lower-case base. This can be used as a marker
 * for looking for specific types of repeats.
 * (One use for this is to only mark STRs that overlap a heterozygous
 * indel region.)
 *
 * Returns a list of rep_ele structs holding the start,end tuples of repeats;
 *         NULL on failure.
 */
rep_ele *find_STR(char *cons, int len, int lower_only);

/*
 * Returns an array of STR vs no-STR values.
 *         0  => non repetitive.
 *         1+ => repeat with consecutive bit-number for repeat size.
 *
 * Eg:  AGGGGAGGAGAAGAC
 *       1111  1111
 *         2222222
 *              444444
 * =>   011331137754440
 */
char *cons_mark_STR(char *cons, int len, int lower_only);

#endif /* _STR_FINDER_H_ */
