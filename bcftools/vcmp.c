/*  vcmp.c -- reference allele utility functions.

    Copyright (C) 2013-2015, 2018 Genome Research Ltd.

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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <ctype.h>
#include "vcmp.h"

struct _vcmp_t
{
    char *dref;
    int ndref, mdref;   // ndref: positive when ref1 longer, negative when ref2 is longer
    int nmatch;
    int *map, mmap, nmap;
    int *map_dip, mmap_dip, nmap_dip;
};

vcmp_t *vcmp_init()
{
    return (vcmp_t*)calloc(1,sizeof(vcmp_t));
}

void vcmp_destroy(vcmp_t *vcmp)
{
    free(vcmp->map_dip);
    free(vcmp->map);
    free(vcmp->dref);
    free(vcmp);
}

int vcmp_set_ref(vcmp_t *vcmp, char *ref1, char *ref2)
{
    vcmp->ndref = 0;

    char *a = ref1, *b = ref2;
    while ( *a && *b && toupper(*a)==toupper(*b) ) { a++; b++; }
    if ( !*a && !*b ) return 0;
    if ( *a && *b ) return -1;  // refs not compatible

    int i;
    if ( *a )   // ref1 is longer
    {
        vcmp->nmatch = b-ref2;
        while ( *a ) a++;
        vcmp->ndref = (a-ref1) - vcmp->nmatch;
        hts_expand(char,vcmp->ndref+1,vcmp->mdref,vcmp->dref);
        for (i=0; i<vcmp->ndref; i++) vcmp->dref[i] = toupper(ref1[vcmp->nmatch+i]);
        vcmp->dref[vcmp->ndref] = 0;
        return 0;
    }

    // ref2 is longer
    vcmp->nmatch = a-ref1;
    while ( *b ) b++;
    vcmp->ndref = (b-ref2) - vcmp->nmatch;
    hts_expand(char,vcmp->ndref+1,vcmp->mdref,vcmp->dref);
    for (i=0; i<vcmp->ndref; i++) vcmp->dref[i] = toupper(ref2[vcmp->nmatch+i]);
    vcmp->dref[vcmp->ndref] = 0;
    vcmp->ndref *= -1;
    return 0;
}

int vcmp_find_allele(vcmp_t *vcmp, char **als1, int nals1, char *al2)
{
    int i, j;
    for (i=0; i<nals1; i++)
    {
        char *a = als1[i], *b = al2;
        while ( *a && *b && toupper(*a)==toupper(*b) ) { a++; b++; }
        if ( *a && *b ) continue;   // mismatch
        if ( !vcmp->ndref )
        {
            if ( !*a && !*b ) break;    // found
            continue;
        }

        // the prefixes match
        if ( *a )
        {
            if ( vcmp->ndref<0 ) continue;
            for (j=0; j<vcmp->ndref; j++)
                if ( !a[j] || toupper(a[j])!=vcmp->dref[j] ) break;
            if ( j!=vcmp->ndref || a[j] ) continue;
            break;  // found
        }

        if ( vcmp->ndref>0 ) continue;
        for (j=0; j<-vcmp->ndref; j++)
            if ( !b[j] || toupper(b[j])!=vcmp->dref[j] ) break;
        if ( j!=-vcmp->ndref || b[j] ) continue;
        break;  // found
    }
    if (i==nals1) return -1;
    return i;
}


int *vcmp_map_ARvalues(vcmp_t *vcmp, int n, int nals1, char **als1, int nals2, char **als2)
{
    if ( vcmp_set_ref(vcmp,als1[0],als2[0]) < 0 ) return NULL;

    vcmp->nmap = n;
    hts_expand(int, vcmp->nmap, vcmp->mmap, vcmp->map);

    int i, ifrom = n==nals2 ? 0 : 1;
    for (i=ifrom; i<nals2; i++)
    {
        vcmp->map[i-ifrom] = vcmp_find_allele(vcmp, als1+ifrom, nals1-ifrom, als2[i]);
    }
    return vcmp->map;
}

int *vcmp_map_dipGvalues(vcmp_t *vcmp, int *nmap)
{
    vcmp->nmap_dip = vcmp->nmap*(vcmp->nmap+1)/2;
    hts_expand(int, vcmp->nmap_dip, vcmp->mmap_dip, vcmp->map_dip);

    int i, j, k = 0;
    for (i=0; i<vcmp->nmap; i++)
    {
        for (j=0; j<=i; j++)
        {
            vcmp->map_dip[k] = vcmp->map[i]>=0 && vcmp->map[j]>=0 ? bcf_alleles2gt(vcmp->map[i],vcmp->map[j]) : -1;
            k++;
        }
    }
    *nmap = k;
    return vcmp->map_dip;
}


