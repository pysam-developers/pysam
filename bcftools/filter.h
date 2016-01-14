/*  filter.h -- filter expressions.

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

#ifndef __FILTER_H__
#define __FILTER_H__

#include <htslib/vcf.h>

typedef struct _filter_t filter_t;

/**
  *  @hdr:  BCF header file
  *  @str:  see the bcftools filter command help for description
  */
filter_t *filter_init(bcf_hdr_t *hdr, const char *str);

void filter_destroy(filter_t *filter);

/**
  *  filter_test() - test whether the BCF record passes the test
  *  @samples:  if not NULL, a pointer to an array with samples statuses is
  *             stored in the location referenced by @samples. The pointer
  *             will be set to NULL if the FORMAT fields were not queried.
  *  Returns 1 if the expression is true and 0 if false.
  */
int filter_test(filter_t *filter, bcf1_t *rec, const uint8_t **samples);

void filter_expression_info(FILE *fp);
int filter_max_unpack(filter_t *filter);

#endif
