/* 
    Copyright (C) 2014 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

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
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/

/*
    Lookup from region and sex to ploidy.

    Example of usage:

        int default_ploidy = 2;
        ploidy_t *pld = ploidy_init(fname, default_ploidy);

        int nsex = ploidy_nsex(pld);
        int *sex2ploidy = malloc(sizeof(int)*nsex);

        ploidy_query(pld, "X",60000, sex2ploidy, NULL, NULL);
        for (i=0; i<nsex; i++)
            printf("ploidy of %s is %d\n", ploidy_id2sex(pld,i), sex2ploidy[i]);

        ploidy_destroy(pld);

    An example of ploidy file format follows. The coordinates are 1-based and
    inclusive. The "*" records define the default ploidy for each sex. If not
    present, the default_ploidy passed to ploidy_init is used instead:
        X 1 60000 M 1
        X 2699521 154931043 M 1
        Y 1 59373566 M 1
        Y 1 59373566 F 0
        MT 1 16569 M 1
        MT 1 16569 F 1
        *  * *     M 2
        *  * *     F 2
*/

#ifndef __PLOIDY_H__
#define __PLOIDY_H__

#include "regidx.h"

typedef struct _ploidy_t ploidy_t;

/*
 *  ploidy_init()
 *  @param fname:   input file name or NULL if default ploidy from example above should be used
 *  @param dflt:    default ploidy to use for unlisted regions (the '* * *' records have precedence).
 *
 *  Returns new structure on success or NULL on error.
 */
ploidy_t *ploidy_init(const char *fname, int dflt);

/* Same as ploidy_init() but the whole file is passed as a single string */
ploidy_t *ploidy_init_string(const char *str, int dflt);

/*
 *  ploidy_destroy() - free memory allocated by ploidy_init
 */
void ploidy_destroy(ploidy_t *ploidy);

/*
 *  ploidy_query() - query ploidy at a position for all genders at once
 *  @param seq: chromosome name
 *  @param pos: 0-based position
 *  @param sex2ploidy:  if not NULL, array will be filled with mapping from sex id to ploidy
 *  @param min: if not NULL, minimum encountered encountered will be set
 *  @param max: if not NULL, maximum encountered encountered will be set
 *
 *  Returns 1 if the position is listed in the regions or 0 otherwise.
 */
int ploidy_query(ploidy_t *ploidy, char *seq, int pos, int *sex2ploidy, int *min, int *max);

/*
 *  ploidy_nsex() - return number of recognised genders
 */
int ploidy_nsex(ploidy_t *ploidy);

/*
 *  ploidy_id2sex() - mapping between numeric gender id and name
 *
 *  Returns gender name (e.g. "M" or "F" in the example above)
 *  or NULL if there is no such id.
 */
char *ploidy_id2sex(ploidy_t *ploidy, int id);

/*
 *  ploidy_sex2id() - mapping between gender name and its numeric id
 *
 *  Returns numeric id or -1 if not present.
 */
int ploidy_sex2id(ploidy_t *ploidy, char *sex);

/*
 *  ploidy_add_sex() - register new gender name. This function is
 *      useful when gender has the default ploidy for all regions
 *      and is not listed in the file passed to ploidy_init()
 *
 *  Returns numeric id of the added sex, regardless of whether the string was
 *  newly added or was already present in the dictionary.
 */
int ploidy_add_sex(ploidy_t *ploidy, const char *sex);

/** Returns region index for raw access */
regidx_t *ploidy_regions(ploidy_t *ploidy);

/** Return the minimum / maximum recognised ploidy  */
int ploidy_max(ploidy_t *ploidy);
int ploidy_min(ploidy_t *ploidy);

#endif

