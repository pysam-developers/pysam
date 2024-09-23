/*  filter.h -- filter expressions.

    Copyright (C) 2013-2024 Genome Research Ltd.

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
  *  @str:  see the bcftools filter command help for description.
  *         See also the extended usage described in filter_test_ext(),
  *         intended for programmatic access
  *  Same as filter_parse() but exits on errors
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

/**
  *  filter_test_ext() - same as filter_test(), but sets some of the terms
  *        on the fly. An expression initialized with, say,
  *        "STR_TAG={} | INT_TAG={} | FLT_TAG={}" takes three
  *        additional pointer arguments which are expected to point to memory
  *        area occupied by the appropriate type, see also filter_ext_types().
  *        The type determination is not fool-proof, in such case the type can
  *        be given explicitly as eg "TAG={str}".
  *  @ext: array of size 'n_ext' occupied with pointers to the data types
  *        inferred from the expression given at the time of initialization.
  *        The pointers set to NULL will be treated as if missing value "."
  *        was given.
  *  @n_ext: the size of 'ext' array
  */
int filter_test_ext(filter_t *filter, bcf1_t *rec, const uint8_t **samples, const void **ext);

/**
  *  filter_set_samples() - restrict filtering expression to samples.
  *             Call after filter_init().
  *  @samples:  use samples set to 1, ignore samples set 0
  */
void filter_set_samples(filter_t *filter, const uint8_t *samples);

/**
  *  filter_get_doubles() - return a pointer to values from the last filter_test() evaluation
  */
const double *filter_get_doubles(filter_t *filter, int *nval, int *nval1);

int filter_max_unpack(filter_t *filter);

/**
  *  filter_ext_types() - returns the number and BCF_HT_* types of external values
  *         found in the filtering expression
  */
const int *filter_ext_types(filter_t *filter, int *n_ext);

/**
  *  Same as filter_init() but may not exit on some type of errors. The caller
  *  must check if the returned value is not NULL and if the consequent call
  *  of filter_status() returns FILTER_OK before the filter_pass() can be called.
  */
filter_t *filter_parse(bcf_hdr_t *hdr, const char *str);

#define FILTER_OK 0
#define FILTER_ERR_UNKN_TAGS 1
#define FILTER_ERR_OTHER 2

/**
  *  Check if filter_parse() was successful
  */
int filter_status(filter_t *filter);
const char **filter_list_undef_tags(filter_t *filter, int *nundef);
const char **filter_list_used_tags(filter_t *filter, int *nused);

#endif
