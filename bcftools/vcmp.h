/*  vcmp.h -- reference allele utility functions.

    Copyright (C) 2013-2014 Genome Research Ltd.

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
THE SOFTWARE.  */

#ifndef __VCMP_H__
#define __VCMP_H__

typedef struct _vcmp_t vcmp_t;

vcmp_t *vcmp_init(void);
void vcmp_destroy(vcmp_t *vcmp);

/*
 *  vcmp_set_ref() - sets and compares reference alleles
 *  Returns 0 on success or -1 if alleles not compatible
 */
int vcmp_set_ref(vcmp_t *vcmp, char *ref1, char *ref2);

/*
 *  vcmp_find_allele()
 *  @param als1:  alternate alleles to ref1 above
 *  @param al2:   alternate allele to ref2 above
 *  Returns -1 if not found or 0-based index to als1 of matching allele
 */
int vcmp_find_allele(vcmp_t *vcmp, char **als1, int nals1, char *al2);

/*
 *  vcmp_map_ARvalues() - Create mapping for Number=A,R tag values
 *  @param number:  nals2 for Number=R, nals2-1 for Number=A
 *  @param nals1:   number of alleles
 *  @param als1:    alleles
 *
 *  Returns pointer to an array of size nals2 with mapping from als2
 *  to als1 or NULL if REFs (als1[0] and als2[0]) are not compatible.
 *  If i is the index of an allele in als2, ret[i] is the index of matching
 *  allele in als1 or -1 if als2 does not have a matching allele.
 *  The caller must not free the array.
 */
int *vcmp_map_ARvalues(vcmp_t *vcmp, int number, int nals1, char **als1, int nals2, char **als2);


#endif
