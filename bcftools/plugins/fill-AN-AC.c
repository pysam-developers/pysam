/*  plugins/fill-AN-AC.c -- fills AN and AC INFO fields.

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

bcf_hdr_t *in_hdr, *out_hdr;
int *arr = NULL, marr = 0;

const char *about(void)
{
    return "Fill INFO fields AN and AC.\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    in_hdr  = in;
    out_hdr = out;
    bcf_hdr_append(out_hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
    bcf_hdr_append(out_hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    return 0;
}

bcf1_t *process(bcf1_t *rec)
{
    hts_expand(int,rec->n_allele,marr,arr);
    int ret = bcf_calc_ac(in_hdr,rec,arr,BCF_UN_FMT);
    if ( ret>0 )
    {
        int i, an = 0;
        for (i=0; i<rec->n_allele; i++) an += arr[i];
        bcf_update_info_int32(out_hdr, rec, "AN", &an, 1);
        bcf_update_info_int32(out_hdr, rec, "AC", arr+1, rec->n_allele-1);
    }
    return rec;
}

void destroy(void)
{
    free(arr);
}


