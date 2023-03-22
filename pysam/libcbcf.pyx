# cython: language_level=3
# cython: embedsignature=True
# cython: profile=True
###############################################################################
###############################################################################
## Cython wrapper for htslib VCF/BCF reader/writer
###############################################################################
#
# NOTICE: This code is incomplete and preliminary.  It offers a nearly
#         complete Pythonic interface to VCF/BCF metadata and data with
#         reading and writing capability.  Documentation and a unit test suite
#         are in the works.  The code is best tested under Python 2, but
#         should also work with Python 3.  Please report any remaining
#         str/bytes issues on the github site when using Python 3 and I'll
#         fix them promptly.
#
# Here is a minimal example of how to use the API:
#
#     $ cat bcfview.py
#     import sys
#     from pysam import VariantFile
#
#     bcf_in = VariantFile(sys.argv[1]) # auto-detect input format
#     bcf_out = VariantFile('-', 'w', header=bcf_in.header)
#
#     for rec in bcf_in:
#         bcf_out.write(rec)
#
# Performance is fairly close to that of bcftools view.  Here is an example
# using some 1k Genomes data:
#
#     $ time python bcfview.py ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf |wc -l
#      1103799
#
#     real	0m56.114s
#     user	1m4.489s
#     sys	0m3.102s
#
#     $ time bcftools view ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf |wc -l
#      1103800  # bcftools adds an extra header
#
#     real	0m55.126s
#     user	1m3.502s
#     sys	0m3.459s
#
###############################################################################
#
# TODO list:
#
#   * more genotype methods
#   * unit test suite (perhaps py.test based)
#   * documentation
#   * pickle support
#   * left/right locus normalization
#   * fix reopen to re-use fd
#
###############################################################################
#
# The MIT License
#
# Copyright (c) 2015,2016 Kevin Jacobs (jacobs@bioinformed.com)
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
#
###############################################################################

from __future__ import division, print_function

import os
import sys

from libc.errno  cimport errno, EPIPE
from libc.string cimport strcmp, strpbrk, strerror
from libc.stdint cimport INT8_MAX, INT16_MAX, INT32_MAX

cimport cython

from cpython.object  cimport PyObject
from cpython.ref     cimport Py_INCREF
from cpython.dict    cimport PyDict_GetItemString, PyDict_SetItemString
from cpython.tuple   cimport PyTuple_New, PyTuple_SET_ITEM
from cpython.bytes   cimport PyBytes_FromStringAndSize
from cpython.unicode cimport PyUnicode_DecodeUTF8

from pysam.libchtslib cimport HTSFile, hisremote

from pysam.utils import unquoted_str


__all__ = ['VariantFile',
           'VariantHeader',
           'VariantHeaderRecord',
           'VariantHeaderRecords',
           'VariantMetadata',
           'VariantHeaderMetadata',
           'VariantContig',
           'VariantHeaderContigs',
           'VariantHeaderSamples',
           'VariantRecordFilter',
           'VariantRecordFormat',
           'VariantRecordInfo',
           'VariantRecordSamples',
           'VariantRecord',
           'VariantRecordSample',
           'BaseIndex',
           'BCFIndex',
           'TabixIndex',
           'BaseIterator',
           'BCFIterator',
           'TabixIterator',
           'VariantRecord']

########################################################################
########################################################################
## Constants
########################################################################

cdef int MAX_POS = (1 << 31) - 1
cdef tuple VALUE_TYPES = ('Flag', 'Integer', 'Float', 'String')
cdef tuple METADATA_TYPES = ('FILTER', 'INFO', 'FORMAT', 'CONTIG', 'STRUCTURED', 'GENERIC')
cdef tuple METADATA_LENGTHS = ('FIXED', 'VARIABLE', 'A', 'G', 'R')


########################################################################
########################################################################
## Python 3 compatibility functions
########################################################################

from pysam.libcutils cimport force_bytes, force_str, charptr_to_str, charptr_to_str_w_len
from pysam.libcutils cimport encode_filename, from_string_and_size, decode_bytes


########################################################################
########################################################################
## Sentinel object
########################################################################

cdef object _nothing = object()

########################################################################
########################################################################
## VCF/BCF string intern system
########################################################################

cdef dict bcf_str_cache = {}

cdef inline bcf_str_cache_get_charptr(const char* s):
    if s == NULL:
        return None

    cdef PyObject *pystr = PyDict_GetItemString(bcf_str_cache, s)
    if pystr:
        return <object>pystr

    val = PyUnicode_DecodeUTF8(s, strlen(s), NULL)

    PyDict_SetItemString(bcf_str_cache, s, val)

    return val


########################################################################
########################################################################
## Genotype math
########################################################################

cdef int comb(int n, int k) except -1:
    """Return binomial coefficient: n choose k

    >>> comb(5, 1)
    5
    >>> comb(5, 2)
    10
    >>> comb(2, 2)
    1
    >>> comb(100, 2)
    4950
    """
    if k > n:
        return 0
    elif k == n:
        return 1
    elif k > n // 2:
        k = n - k

    cdef d, result

    d = result = n - k + 1
    for i in range(2, k + 1):
        d += 1
        result  *= d
        result //= i
    return result


cdef inline int bcf_geno_combinations(int ploidy, int alleles) except -1:
    """Return the count of genotypes expected for the given ploidy and number of alleles.

    >>> bcf_geno_combinations(1, 2)
    2
    >>> bcf_geno_combinations(2, 2)
    3
    >>> bcf_geno_combinations(2, 3)
    6
    >>> bcf_geno_combinations(3, 2)
    4
    """
    return comb(alleles + ploidy - 1, ploidy)


########################################################################
########################################################################
## Low level type conversion helpers
########################################################################


cdef inline bint check_header_id(bcf_hdr_t *hdr, int hl_type, int id):
    return id >= 0 and id < hdr.n[BCF_DT_ID] and bcf_hdr_idinfo_exists(hdr, hl_type, id)


cdef inline int is_gt_fmt(bcf_hdr_t *hdr, int fmt_id):
    return strcmp(bcf_hdr_int2id(hdr, BCF_DT_ID, fmt_id), 'GT') == 0


cdef inline int bcf_genotype_count(bcf_hdr_t *hdr, bcf1_t *rec, int sample) except -1:

    if sample < 0:
        raise ValueError('genotype is only valid as a format field')

    cdef int32_t *gt_arr = NULL
    cdef int ngt = 0
    ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt)

    if ngt <= 0 or not gt_arr:
        return 0

    assert ngt % rec.n_sample == 0
    cdef int max_ploidy = ngt // rec.n_sample
    cdef int32_t *gt = gt_arr + sample * max_ploidy
    cdef int ploidy = 0

    while ploidy < max_ploidy and gt[0] != bcf_int32_vector_end:
        gt += 1
        ploidy += 1

    free(<void*>gt_arr)

    return bcf_geno_combinations(ploidy, rec.n_allele)


cdef tuple char_array_to_tuple(const char **a, ssize_t n, int free_after=0):
    if not a:
        return None
    try:
         return tuple(charptr_to_str(a[i]) for i in range(n))
    finally:
        if free_after and a:
            free(a)


cdef bcf_array_to_object(void *data, int type, ssize_t n, ssize_t count, int scalar):
    cdef char    *datac
    cdef int8_t  *data8
    cdef int16_t *data16
    cdef int32_t *data32
    cdef float   *dataf
    cdef int      i
    cdef bytes    b

    if not data or n <= 0:
        return None

    if type == BCF_BT_CHAR:
        datac = <char *>data

        if not n:
            value = ()
        else:
            # Check if at least one null terminator is present
            if datac[n-1] == bcf_str_vector_end:
                # If so, create a string up to the first null terminator
                b = datac
            else:
                # Otherwise, copy the entire block
                b = datac[:n]
            value = tuple(decode_bytes(v, 'utf-8') if v and v != bcf_str_missing else None for v in b.split(b','))
    else:
        value = []
        if type == BCF_BT_INT8:
            data8 = <int8_t *>data
            for i in range(n):
                if data8[i] == bcf_int8_vector_end:
                    break
                value.append(data8[i] if data8[i] != bcf_int8_missing else None)
        elif type == BCF_BT_INT16:
            data16 = <int16_t *>data
            for i in range(n):
                if data16[i] == bcf_int16_vector_end:
                    break
                value.append(data16[i] if data16[i] != bcf_int16_missing else None)
        elif type == BCF_BT_INT32:
            data32 = <int32_t *>data
            for i in range(n):
                if data32[i] == bcf_int32_vector_end:
                    break
                value.append(data32[i] if data32[i] != bcf_int32_missing else None)
        elif type == BCF_BT_FLOAT:
            dataf = <float *>data
            for i in range(n):
                if bcf_float_is_vector_end(dataf[i]):
                    break
                value.append(dataf[i] if not bcf_float_is_missing(dataf[i]) else None)
        else:
            raise TypeError('unsupported info type code')

    # FIXME: Need to know length?  Report errors?  Pad with missing values?  Not clear what to do.
    if not value:
        if scalar:
            value = None
        elif count <= 0:
            value = ()
        else:
            value = (None,)*count
    elif scalar and len(value) == 1:
        value = value[0]
    else:
        value = tuple(value)

    return value


cdef bcf_object_to_array(values, void *data, int bt_type, ssize_t n, int vlen):
    cdef char    *datac
    cdef int8_t  *data8
    cdef int16_t *data16
    cdef int32_t *data32
    cdef float   *dataf
    cdef ssize_t  i, value_count = len(values)

    assert value_count <= n

    if bt_type == BCF_BT_CHAR:
        if not isinstance(values, (str, bytes)):
            values = b','.join(force_bytes(v) if v else bcf_str_missing for v in values)
            value_count = len(values)
        assert value_count <= n
        datac = <char *>data
        memcpy(datac, <char *>values, value_count)
        for i in range(value_count, n):
            datac[i] = 0
    elif bt_type == BCF_BT_INT8:
        datai8 = <int8_t *>data
        for i in range(value_count):
            val = values[i]
            datai8[i] = val if val is not None else bcf_int8_missing
        for i in range(value_count, n):
            datai8[i] = bcf_int8_vector_end
    elif bt_type == BCF_BT_INT16:
        datai16 = <int16_t *>data
        for i in range(value_count):
            val = values[i]
            datai16[i] = val if val is not None else bcf_int16_missing
        for i in range(value_count, n):
            datai16[i] = bcf_int16_vector_end
    elif bt_type == BCF_BT_INT32:
        datai32 = <int32_t *>data
        for i in range(value_count):
            val = values[i]
            datai32[i] = val if val is not None else bcf_int32_missing
        for i in range(value_count, n):
            datai32[i] = bcf_int32_vector_end
    elif bt_type == BCF_BT_FLOAT:
        dataf = <float *>data
        for i in range(value_count):
            val = values[i]
            if val is None:
                bcf_float_set(dataf + i, bcf_float_missing)
            else:
                dataf[i] = val
        for i in range(value_count, n):
            bcf_float_set(dataf + i, bcf_float_vector_end)
    else:
        raise TypeError('unsupported type')


cdef bcf_empty_array(int type, ssize_t n, int vlen):
    cdef char    *datac
    cdef int32_t *data32
    cdef float   *dataf
    cdef int      i

    if n <= 0:
        raise ValueError('Cannot create empty array')

    if type == BCF_HT_STR:
        value = PyBytes_FromStringAndSize(NULL, sizeof(char)*n)
        datac = <char *>value
        for i in range(n):
            datac[i] = bcf_str_missing if not vlen else bcf_str_vector_end
    elif type == BCF_HT_INT:
        value = PyBytes_FromStringAndSize(NULL, sizeof(int32_t)*n)
        data32 = <int32_t *><char *>value
        for i in range(n):
            data32[i] = bcf_int32_missing if not vlen else bcf_int32_vector_end
    elif type == BCF_HT_REAL:
        value = PyBytes_FromStringAndSize(NULL, sizeof(float)*n)
        dataf = <float *><char *>value
        for i in range(n):
            bcf_float_set(dataf + i, bcf_float_missing if not vlen else bcf_float_vector_end)
    else:
        raise TypeError('unsupported header type code')

    return value


cdef bcf_copy_expand_array(void *src_data, int src_type, size_t src_values,
                           void *dst_data, int dst_type, size_t dst_values,
                           int vlen):
    """copy data from src to dest where the size of the elements (src_type/dst_type) differ
    as well as the number of elements (src_values/dst_values).
    """

    cdef char    *src_datac
    cdef char    *dst_datac
    cdef int8_t  *src_datai8
    cdef int16_t *src_datai16
    cdef int32_t *src_datai32
    cdef int32_t *dst_datai
    cdef float   *src_dataf
    cdef float   *dst_dataf
    cdef ssize_t src_size, dst_size, i, j
    cdef int val

    if src_values > dst_values:
        raise ValueError('Cannot copy arrays with src_values={} > dst_values={}'.format(src_values, dst_values))

    if src_type == dst_type == BCF_BT_CHAR:
        src_datac = <char *>src_data
        dst_datac = <char *>dst_data
        memcpy(dst_datac, src_datac, src_values)
        for i in range(src_values, dst_values):
            dst_datac[i] = 0
    elif src_type == BCF_BT_INT8 and dst_type == BCF_BT_INT32:
        src_datai8 = <int8_t *>src_data
        dst_datai  = <int32_t *>dst_data
        for i in range(src_values):
            val = src_datai8[i]
            if val == bcf_int8_missing:
                val = bcf_int32_missing
            elif val == bcf_int8_vector_end:
                val = bcf_int32_vector_end
            dst_datai[i] = val
        for i in range(src_values, dst_values):
            dst_datai[i] = bcf_int32_missing if not vlen else bcf_int32_vector_end
    elif src_type == BCF_BT_INT16 and dst_type == BCF_BT_INT32:
        src_datai16 = <int16_t *>src_data
        dst_datai   = <int32_t *>dst_data
        for i in range(src_values):
            val = src_datai16[i]
            if val == bcf_int16_missing:
                val = bcf_int32_missing
            elif val == bcf_int16_vector_end:
                val = bcf_int32_vector_end
            dst_datai[i] = val
        for i in range(src_values, dst_values):
            dst_datai[i] = bcf_int32_missing if not vlen else bcf_int32_vector_end
    elif src_type == BCF_BT_INT32 and dst_type == BCF_BT_INT32:
        src_datai32 = <int32_t *>src_data
        dst_datai   = <int32_t *>dst_data
        for i in range(src_values):
            dst_datai[i] = src_datai32[i]
        for i in range(src_values, dst_values):
            dst_datai[i] = bcf_int32_missing if not vlen else bcf_int32_vector_end
    elif src_type == BCF_BT_FLOAT and dst_type == BCF_BT_FLOAT:
        src_dataf = <float *>src_data
        dst_dataf = <float *>dst_data
        for i in range(src_values):
            dst_dataf[i] = src_dataf[i]
        for i in range(src_values, dst_values):
            bcf_float_set(dst_dataf + i, bcf_float_missing if not vlen else bcf_float_vector_end)
    else:
        raise TypeError('unsupported types')


cdef bcf_get_value_count(VariantRecord record, int hl_type, int id, ssize_t *count, int *scalar, int sample):
    if record is None:
        raise ValueError('record must not be None')

    cdef bcf_hdr_t *hdr = record.header.ptr
    cdef bcf1_t *r = record.ptr

    if not check_header_id(hdr, hl_type, id):
        raise ValueError('Invalid header')

    cdef int length = bcf_hdr_id2length(hdr, hl_type, id)
    cdef int number = bcf_hdr_id2number(hdr, hl_type, id)

    scalar[0] = 0

    if hl_type == BCF_HL_FMT and is_gt_fmt(hdr, id):
        count[0] = number
    elif length == BCF_VL_FIXED:
        if number == 1:
            scalar[0] = 1
        count[0] = number
    elif length == BCF_VL_R:
        count[0] = r.n_allele
    elif length == BCF_VL_A:
        count[0] = r.n_allele - 1
    elif length == BCF_VL_G:
        count[0] = bcf_genotype_count(hdr, r, sample)
    elif length == BCF_VL_VAR:
        count[0] = -1
    else:
        raise ValueError('Unknown format length')


cdef object bcf_info_get_value(VariantRecord record, const bcf_info_t *z):
    if record is None:
        raise ValueError('record must not be None')

    cdef bcf_hdr_t *hdr = record.header.ptr

    cdef char *s
    cdef ssize_t count
    cdef int scalar

    bcf_get_value_count(record, BCF_HL_INFO, z.key, &count, &scalar, -1)

    if z.len == 0:
        if  bcf_hdr_id2type(hdr, BCF_HL_INFO, z.key) == BCF_HT_FLAG:
            value = True
        elif scalar:
            value = None
        else:
            value = ()
    elif z.len == 1:
        if z.type == BCF_BT_INT8:
            value = z.v1.i if z.v1.i != bcf_int8_missing else None
        elif z.type == BCF_BT_INT16:
            value = z.v1.i if z.v1.i != bcf_int16_missing else None
        elif z.type == BCF_BT_INT32:
            value = z.v1.i if z.v1.i != bcf_int32_missing else None
        elif z.type == BCF_BT_FLOAT:
            value = z.v1.f if not bcf_float_is_missing(z.v1.f) else None
        elif z.type == BCF_BT_CHAR:
            value = force_str(chr(z.v1.i))
        else:
            raise TypeError('unsupported info type code')

        if not scalar and value != ():
            value = (value,)
    else:
        value = bcf_array_to_object(z.vptr, z.type, z.len, count, scalar)

    return value


cdef object bcf_check_values(VariantRecord record, value, int sample,
                             int hl_type, int ht_type,
                             int id, int bt_type, ssize_t bt_len,
                             ssize_t *value_count, int *scalar, int *realloc):

    if record is None:
        raise ValueError('record must not be None')

    bcf_get_value_count(record, hl_type, id, value_count, scalar, sample)

    # Validate values now that we know the type and size
    values = (value,) if not isinstance(value, (list, tuple)) else value

    # Validate values now that we know the type and size
    if ht_type == BCF_HT_FLAG:
        value_count[0] = 1
    elif hl_type == BCF_HL_FMT and is_gt_fmt(record.header.ptr, id):
        # KBJ: htslib lies about the cardinality of GT fields-- they're really VLEN (-1)
        value_count[0] = -1

    cdef int given = len(values)
    if value_count[0] != -1 and value_count[0] != given:
        if scalar[0]:
            raise TypeError('value expected to be scalar, given len={}'.format(given))
        else:
            raise TypeError('values expected to be {}-tuple, given len={}'.format(value_count[0], given))

    if ht_type == BCF_HT_REAL:
        for v in values:
            if not(v is None or isinstance(v, (float, int))):
                raise TypeError('invalid value for Float format')
    elif ht_type == BCF_HT_INT:
        for v in values:
            if not(v is None or (isinstance(v, (float, int)) and int(v) == v)):
                raise TypeError('invalid value for Integer format')
        for v in values:
            if not(v is None or bcf_int32_missing < v <= INT32_MAX):
                raise ValueError('Integer value too small/large to store in VCF/BCF')
    elif ht_type == BCF_HT_STR:
        values = b','.join(force_bytes(v) if v is not None else b'' for v in values)
    elif ht_type == BCF_HT_FLAG:
        if values[0] not in (True, False, None, 1, 0):
            raise ValueError('Flag values must be: True, False, None, 1, 0')
    else:
        raise TypeError('unsupported type')

    realloc[0] = 0
    if len(values) <= 1 and hl_type == BCF_HL_INFO:
        realloc[0] = 0
    elif len(values) > bt_len:
        realloc[0] = 1
    elif bt_type == BCF_BT_INT8:
        for v in values:
            if v is not None and not(bcf_int8_missing < v <= INT8_MAX):
                realloc[0] = 1
                break
    elif bt_type == BCF_BT_INT16:
        for v in values:
            if v is not None and not(bcf_int16_missing < v <= INT16_MAX):
                realloc[0] = 1
                break

    return values


cdef bcf_encode_alleles(VariantRecord record, values):
    if record is None:
        raise ValueError('record must not be None')

    cdef bcf1_t *r = record.ptr
    cdef int32_t nalleles = r.n_allele
    cdef list gt_values = []
    cdef char *s
    cdef int i

    if values is None:
        return ()

    if not isinstance(values, (list, tuple)):
        values = (values,)

    for value in values:
        if value is None:
            gt_values.append(bcf_gt_missing)
        elif isinstance(value, (str, bytes)):
            bvalue = force_bytes(value)
            s = bvalue
            for i in range(r.n_allele):
                if strcmp(r.d.allele[i], s) != 0:
                    gt_values.append(bcf_gt_unphased(i))
                    break
            else:
                raise ValueError('Unknown allele')
        else:
            i = value
            if not (0 <= i < nalleles):
                raise ValueError('Invalid allele index')
            gt_values.append(bcf_gt_unphased(i))

    return gt_values


cdef bcf_info_set_value(VariantRecord record, key, value):
    if record is None:
        raise ValueError('record must not be None')

    cdef bcf_hdr_t *hdr = record.header.ptr
    cdef bcf1_t *r = record.ptr
    cdef int info_id, info_type, scalar, dst_type, realloc, vlen = 0
    cdef ssize_t i, value_count, alloc_len, alloc_size, dst_size

    if bcf_unpack(r, BCF_UN_INFO) < 0:
        raise ValueError('Error unpacking VariantRecord')

    cdef bytes bkey = force_bytes(key)
    cdef bcf_info_t *info = bcf_get_info(hdr, r, bkey)

    if info:
        info_id = info.key
    else:
        info_id = bcf_header_get_info_id(hdr, bkey)

    if info_id < 0:
        raise KeyError('unknown INFO: {}'.format(key))

    if not check_header_id(hdr, BCF_HL_INFO, info_id):
        raise ValueError('Invalid header')

    info_type = bcf_hdr_id2type(hdr, BCF_HL_INFO, info_id)
    values = bcf_check_values(record, value, -1,
                              BCF_HL_INFO, info_type, info_id,
                              info.type if info else -1,
                              info.len  if info else -1,
                              &value_count, &scalar, &realloc)

    if info_type == BCF_HT_FLAG:
        if bcf_update_info(hdr, r, bkey, NULL, bool(values[0]), info_type) < 0:
            raise ValueError('Unable to update INFO values')
        return

    vlen = value_count < 0
    value_count = len(values)

    # DISABLED DUE TO ISSUES WITH THE CRAZY POINTERS
    # If we can, write updated values to existing allocated storage
    if 0 and info and not realloc:
        r.d.shared_dirty |= BCF1_DIRTY_INF

        if value_count == 0:
            info.len = 0
            if not info.vptr:
                info.vptr = <uint8_t *>&info.v1.i

        elif value_count == 1:
            # FIXME: Check if need to free vptr if info.len > 0?
            if info.type == BCF_BT_INT8 or info.type == BCF_BT_INT16 or info.type == BCF_BT_INT32:
                bcf_object_to_array(values, &info.v1.i, BCF_BT_INT32, 1, vlen)
            elif info.type == BCF_BT_FLOAT:
                bcf_object_to_array(values, &info.v1.f, BCF_BT_FLOAT, 1, vlen)
            else:
                raise TypeError('unsupported info type code')

            info.len = 1
            if not info.vptr:
                info.vptr = <uint8_t *>&info.v1.i
        else:
            bcf_object_to_array(values, info.vptr, info.type, info.len, vlen)

        return

    alloc_len = max(1, value_count)
    if info and info.len > alloc_len:
        alloc_len = info.len

    new_values = bcf_empty_array(info_type, alloc_len, vlen)
    cdef char *valp = <char *>new_values

    if info_type == BCF_HT_INT:
        dst_type = BCF_BT_INT32
    elif info_type == BCF_HT_REAL:
        dst_type = BCF_BT_FLOAT
    elif info_type == BCF_HT_STR:
        dst_type = BCF_BT_CHAR
    else:
        raise ValueError('Unsupported INFO type')

    bcf_object_to_array(values, valp, dst_type, alloc_len, vlen)

    if bcf_update_info(hdr, r, bkey, valp, <int>alloc_len, info_type) < 0:
        raise ValueError('Unable to update INFO values')


cdef bcf_info_del_value(VariantRecord record, key):
    if record is None:
        raise ValueError('record must not be None')

    cdef bcf_hdr_t *hdr = record.header.ptr
    cdef bcf1_t *r = record.ptr
    cdef ssize_t value_count
    cdef int scalar

    if bcf_unpack(r, BCF_UN_INFO) < 0:
        raise ValueError('Error unpacking VariantRecord')

    cdef bytes bkey = force_bytes(key)
    cdef bcf_info_t *info = bcf_get_info(hdr, r, bkey)

    if not info:
        raise KeyError(key)

    bcf_get_value_count(record, BCF_HL_INFO, info.key, &value_count, &scalar, -1)

    if value_count <= 0:
        null_value = ()
    elif scalar:
        null_value = None
    else:
        null_value = (None,)*value_count

    bcf_info_set_value(record, bkey, null_value)


cdef bcf_format_get_value(VariantRecordSample sample, key):
    if sample is None:
        raise ValueError('sample must not be None')

    cdef bcf_hdr_t *hdr = sample.record.header.ptr
    cdef bcf1_t *r = sample.record.ptr
    cdef ssize_t count
    cdef int scalar

    if bcf_unpack(r, BCF_UN_ALL) < 0:
        raise ValueError('Error unpacking VariantRecord')

    cdef bytes bkey = force_bytes(key)
    cdef bcf_fmt_t *fmt = bcf_get_fmt(hdr, r, bkey)

    if not fmt or not fmt.p:
        raise KeyError('invalid FORMAT: {}'.format(key))

    if is_gt_fmt(hdr, fmt.id):
        return bcf_format_get_allele_indices(sample)

    bcf_get_value_count(sample.record, BCF_HL_FMT, fmt.id, &count, &scalar, sample.index)

    if fmt.p and fmt.n and fmt.size:
        return bcf_array_to_object(fmt.p + sample.index * fmt.size, fmt.type, fmt.n, count, scalar)
    elif scalar:
        return None
    elif count <= 0:
        return ()
    else:
        return (None,)*count


cdef bcf_format_set_value(VariantRecordSample sample, key, value):
    if sample is None:
        raise ValueError('sample must not be None')

    if key == 'phased':
        sample.phased = bool(value)
        return

    cdef bcf_hdr_t *hdr = sample.record.header.ptr
    cdef bcf1_t *r = sample.record.ptr
    cdef int fmt_id
    cdef vdict_t *d
    cdef khiter_t k
    cdef int fmt_type, scalar, realloc, dst_type, vlen = 0
    cdef ssize_t i, nsamples, value_count, alloc_size, alloc_len, dst_size

    if bcf_unpack(r, BCF_UN_ALL) < 0:
        raise ValueError('Error unpacking VariantRecord')

    cdef bytes bkey = force_bytes(key)
    cdef bcf_fmt_t *fmt = bcf_get_fmt(hdr, r, bkey)

    if fmt:
        fmt_id = fmt.id
    else:
        d = <vdict_t *>hdr.dict[BCF_DT_ID]
        k = kh_get_vdict(d, bkey)

        if k == kh_end(d) or kh_val_vdict(d, k).info[BCF_HL_FMT] & 0xF == 0xF:
            raise KeyError('unknown format: {}'.format(key))

        fmt_id = kh_val_vdict(d, k).id

    if not check_header_id(hdr, BCF_HL_FMT, fmt_id):
        raise ValueError('Invalid header')

    fmt_type = bcf_hdr_id2type(hdr, BCF_HL_FMT, fmt_id)

    if fmt_type == BCF_HT_FLAG:
        raise ValueError('Flag types are not allowed on FORMATs')

    if is_gt_fmt(hdr, fmt_id):
        value = bcf_encode_alleles(sample.record, value)
        # KBJ: GT field is considered to be a string by the VCF header but BCF represents it as INT.
        fmt_type = BCF_HT_INT

    values = bcf_check_values(sample.record, value, sample.index,
                              BCF_HL_FMT, fmt_type, fmt_id,
                              fmt.type if fmt else -1,
                              fmt.n    if fmt else -1,
                              &value_count, &scalar, &realloc)
    vlen = value_count < 0
    value_count = len(values)

    # If we can, write updated values to existing allocated storage.
    if fmt and not realloc:
        r.d.indiv_dirty = 1
        bcf_object_to_array(values, fmt.p + sample.index * fmt.size, fmt.type, fmt.n, vlen)
        return

    alloc_len = max(1, value_count)
    if fmt and fmt.n > alloc_len:
        alloc_len = fmt.n

    nsamples = r.n_sample
    new_values = bcf_empty_array(fmt_type, nsamples * alloc_len, vlen)
    cdef char *new_values_p = <char *>new_values

    if fmt_type == BCF_HT_INT:
        dst_type = BCF_BT_INT32
        dst_size = sizeof(int32_t) * alloc_len
    elif fmt_type == BCF_HT_REAL:
        dst_type = BCF_BT_FLOAT
        dst_size = sizeof(float) * alloc_len
    elif fmt_type == BCF_HT_STR:
        dst_type = BCF_BT_CHAR
        dst_size = sizeof(char) * alloc_len
    else:
        raise ValueError('Unsupported FORMAT type')

    if fmt and nsamples > 1:
        for i in range(nsamples):
            bcf_copy_expand_array(fmt.p + i * fmt.size, fmt.type, fmt.n,
                                  new_values_p  + i * dst_size, dst_type, alloc_len,
                                  vlen)

    bcf_object_to_array(values, new_values_p + sample.index * dst_size, dst_type, alloc_len, vlen)

    if bcf_update_format(hdr, r, bkey, new_values_p, <int>(nsamples * alloc_len), fmt_type) < 0:
        raise ValueError('Unable to update format values')


cdef bcf_format_del_value(VariantRecordSample sample, key):
    if sample is None:
        raise ValueError('sample must not be None')

    cdef bcf_hdr_t *hdr = sample.record.header.ptr
    cdef bcf1_t *r = sample.record.ptr
    cdef ssize_t value_count
    cdef int scalar

    if bcf_unpack(r, BCF_UN_ALL) < 0:
        raise ValueError('Error unpacking VariantRecord')

    cdef bytes bkey = force_bytes(key)
    cdef bcf_fmt_t *fmt = bcf_get_fmt(hdr, r, bkey)

    if not fmt or not fmt.p:
        raise KeyError(key)

    bcf_get_value_count(sample.record, BCF_HL_FMT, fmt.id, &value_count, &scalar, sample.index)

    if value_count <= 0:
        null_value = ()
    elif scalar:
        null_value = None
    else:
        null_value = (None,)*value_count

    bcf_format_set_value(sample, bkey, null_value)


cdef bcf_format_get_allele_indices(VariantRecordSample sample):
    if sample is None:
        raise ValueError('sample must not be None')

    cdef bcf_hdr_t *hdr = sample.record.header.ptr
    cdef bcf1_t *r = sample.record.ptr
    cdef int32_t n = r.n_sample

    if bcf_unpack(r, BCF_UN_ALL) < 0:
        raise ValueError('Error unpacking VariantRecord')

    if sample.index < 0 or sample.index >= n or not r.n_fmt:
        return ()

    cdef bcf_fmt_t *fmt0 = r.d.fmt
    cdef int gt0 = is_gt_fmt(hdr, fmt0.id)

    if not gt0 or not fmt0.n:
        return ()

    cdef int8_t  *data8
    cdef int16_t *data16
    cdef int32_t *data32
    cdef int32_t a, nalleles = r.n_allele
    cdef list alleles = []

    if fmt0.type == BCF_BT_INT8:
        data8 = <int8_t *>(fmt0.p + sample.index * fmt0.size)
        for i in range(fmt0.n):
            if data8[i] == bcf_int8_vector_end:
                break
            elif data8[i] == bcf_gt_missing:
                a = -1
            else:
                a = bcf_gt_allele(data8[i])
            alleles.append(a if 0 <= a < nalleles else None)
    elif fmt0.type == BCF_BT_INT16:
        data16 = <int16_t *>(fmt0.p + sample.index * fmt0.size)
        for i in range(fmt0.n):
            if data16[i] == bcf_int16_vector_end:
                break
            elif data16[i] == bcf_gt_missing:
                a = -1
            else:
                a = bcf_gt_allele(data16[i])
            alleles.append(a if 0 <= a < nalleles else None)
    elif fmt0.type == BCF_BT_INT32:
        data32 = <int32_t *>(fmt0.p + sample.index * fmt0.size)
        for i in range(fmt0.n):
            if data32[i] == bcf_int32_vector_end:
                break
            elif data32[i] == bcf_gt_missing:
                a = -1
            else:
                a = bcf_gt_allele(data32[i])
            alleles.append(a if 0 <= a < nalleles else None)

    return tuple(alleles)


cdef bcf_format_get_alleles(VariantRecordSample sample):
    if sample is None:
        raise ValueError('sample must not be None')

    cdef bcf_hdr_t *hdr = sample.record.header.ptr
    cdef bcf1_t *r = sample.record.ptr
    cdef int32_t nsamples = r.n_sample

    if bcf_unpack(r, BCF_UN_ALL) < 0:
        raise ValueError('Error unpacking VariantRecord')

    cdef int32_t nalleles = r.n_allele

    if sample.index < 0 or sample.index >= nsamples or not r.n_fmt:
        return ()

    cdef bcf_fmt_t *fmt0 = r.d.fmt
    cdef int gt0 = is_gt_fmt(hdr, fmt0.id)

    if not gt0 or not fmt0.n:
        return ()

    cdef int32_t  a
    cdef int8_t  *data8
    cdef int16_t *data16
    cdef int32_t *data32
    alleles = []
    if fmt0.type == BCF_BT_INT8:
        data8 = <int8_t *>(fmt0.p + sample.index * fmt0.size)
        for i in range(fmt0.n):
            if data8[i] == bcf_int8_vector_end:
                break
            a = bcf_gt_allele(data8[i])
            alleles.append(charptr_to_str(r.d.allele[a]) if 0 <= a < nalleles else None)
    elif fmt0.type == BCF_BT_INT16:
        data16 = <int16_t *>(fmt0.p + sample.index * fmt0.size)
        for i in range(fmt0.n):
            if data16[i] == bcf_int16_vector_end:
                break
            a = bcf_gt_allele(data16[i])
            alleles.append(charptr_to_str(r.d.allele[a]) if 0 <= a < nalleles else None)
    elif fmt0.type == BCF_BT_INT32:
        data32 = <int32_t *>(fmt0.p + sample.index * fmt0.size)
        for i in range(fmt0.n):
            if data32[i] == bcf_int32_vector_end:
                break
            a = bcf_gt_allele(data32[i])
            alleles.append(charptr_to_str(r.d.allele[a]) if 0 <= a < nalleles else None)
    return tuple(alleles)


cdef bint bcf_sample_get_phased(VariantRecordSample sample):
    if sample is None:
        raise ValueError('sample must not be None')

    cdef bcf_hdr_t *hdr = sample.record.header.ptr
    cdef bcf1_t *r = sample.record.ptr
    cdef int32_t n = r.n_sample

    if bcf_unpack(r, BCF_UN_ALL) < 0:
        raise ValueError('Error unpacking VariantRecord')

    if sample.index < 0 or sample.index >= n or not r.n_fmt:
        return False

    cdef bcf_fmt_t *fmt0 = r.d.fmt
    cdef int gt0 = is_gt_fmt(hdr, fmt0.id)

    if not gt0 or not fmt0.n:
        return False

    cdef int8_t  *data8
    cdef int16_t *data16
    cdef int32_t *data32

    cdef bint phased = False

    if fmt0.type == BCF_BT_INT8:
        data8 = <int8_t *>(fmt0.p + sample.index * fmt0.size)
        for i in range(fmt0.n):
            if data8[i] == bcf_int8_vector_end:
                break
            elif data8[i] == bcf_int8_missing:
                continue
            elif i and not bcf_gt_is_phased(data8[i]):
                return False
            else:
                phased = True
    elif fmt0.type == BCF_BT_INT16:
        data16 = <int16_t *>(fmt0.p + sample.index * fmt0.size)
        for i in range(fmt0.n):
            if data16[i] == bcf_int16_vector_end:
                break
            elif data16[i] == bcf_int16_missing:
                continue
            elif i and not bcf_gt_is_phased(data16[i]):
                return False
            else:
                phased = True
    elif fmt0.type == BCF_BT_INT32:
        data32 = <int32_t *>(fmt0.p + sample.index * fmt0.size)
        for i in range(fmt0.n):
            if data32[i] == bcf_int32_vector_end:
                break
            elif data32[i] == bcf_int32_missing:
                continue
            elif i and not bcf_gt_is_phased(data32[i]):
                return False
            else:
                phased = True

    return phased


cdef bcf_sample_set_phased(VariantRecordSample sample, bint phased):
    if sample is None:
        raise ValueError('sample must not be None')

    cdef bcf_hdr_t *hdr = sample.record.header.ptr
    cdef bcf1_t *r = sample.record.ptr
    cdef int32_t n = r.n_sample

    if bcf_unpack(r, BCF_UN_ALL) < 0:
        raise ValueError('Error unpacking VariantRecord')

    if sample.index < 0 or sample.index >= n or not r.n_fmt:
        return

    cdef bcf_fmt_t *fmt0 = r.d.fmt
    cdef int gt0 = is_gt_fmt(hdr, fmt0.id)

    if not gt0 or not fmt0.n:
        raise ValueError('Cannot set phased before genotype is set')

    cdef int8_t  *data8
    cdef int16_t *data16
    cdef int32_t *data32

    if fmt0.type == BCF_BT_INT8:
        data8 = <int8_t *>(fmt0.p + sample.index * fmt0.size)
        for i in range(fmt0.n):
            if data8[i] == bcf_int8_vector_end:
                break
            elif data8[i] == bcf_int8_missing:
                continue
            elif i:
                data8[i] = (data8[i] & 0xFE) | phased
    elif fmt0.type == BCF_BT_INT16:
        data16 = <int16_t *>(fmt0.p + sample.index * fmt0.size)
        for i in range(fmt0.n):
            if data16[i] == bcf_int16_vector_end:
                break
            elif data16[i] == bcf_int16_missing:
                continue
            elif i:
                data16[i] = (data16[i] & 0xFFFE) | phased
    elif fmt0.type == BCF_BT_INT32:
        data32 = <int32_t *>(fmt0.p + sample.index * fmt0.size)
        for i in range(fmt0.n):
            if data32[i] == bcf_int32_vector_end:
                break
            elif data32[i] == bcf_int32_missing:
                continue
            elif i:
                data32[i] = (data32[i] & 0xFFFFFFFE) | phased


cdef inline bcf_sync_end(VariantRecord record):
    cdef bcf_hdr_t *hdr = record.header.ptr
    cdef bcf_info_t *info
    cdef int end_id = bcf_header_get_info_id(record.header.ptr, b'END')
    cdef int ref_len

    # allow missing ref when instantiating a new record
    if record.ref is not None:
        ref_len = len(record.ref)
    else:
        ref_len = 0

    # Delete INFO/END if no alleles are present or if rlen is equal to len(ref)
    # Always keep END for symbolic alleles
    if not has_symbolic_allele(record) and (not record.ptr.n_allele or record.ptr.rlen == ref_len):
        # If INFO/END is not defined in the header, it doesn't exist in the record
        if end_id >= 0:
            info = bcf_get_info(hdr, record.ptr, b'END')
            if info and info.vptr:
                if bcf_update_info(hdr, record.ptr, b'END', NULL, 0, info.type) < 0:
                    raise ValueError('Unable to delete END')
    else:
        # Create END header, if not present
        if end_id < 0:
            record.header.info.add('END', number=1, type='Integer', description='Stop position of the interval')

        # Update to reflect stop position
        bcf_info_set_value(record, b'END', record.ptr.pos + record.ptr.rlen)


cdef inline int has_symbolic_allele(VariantRecord record):
    """Return index of first symbolic allele. 0 if no symbolic alleles."""

    for i in range(1, record.ptr.n_allele):
        alt = record.ptr.d.allele[i]
        if alt[0] == b'<' and alt[len(alt) - 1] == b'>':
            return i

    return 0


########################################################################
########################################################################
## Variant Header objects
########################################################################


cdef bcf_header_remove_hrec(VariantHeader header, int i):
    if header is None:
        raise ValueError('header must not be None')

    cdef bcf_hdr_t *hdr = header.ptr

    if i < 0 or i >= hdr.nhrec:
        raise ValueError('Invalid header record index')

    cdef bcf_hrec_t *hrec = hdr.hrec[i]
    hdr.nhrec -= 1

    if i < hdr.nhrec:
        memmove(&hdr.hrec[i], &hdr.hrec[i+1], (hdr.nhrec-i)*sizeof(bcf_hrec_t*))

    bcf_hrec_destroy(hrec)
    hdr.hrec[hdr.nhrec] = NULL
    hdr.dirty = 1


#FIXME: implement a full mapping interface
#FIXME: passing bcf_hrec_t* is not safe, since we cannot control the
#       object lifetime.
cdef class VariantHeaderRecord(object):
    """header record from a :class:`VariantHeader` object"""
    def __init__(self, *args, **kwargs):
        raise TypeError('this class cannot be instantiated from Python')

    @property
    def type(self):
        """header type: FILTER, INFO, FORMAT, CONTIG, STRUCTURED, or GENERIC"""
        cdef bcf_hrec_t *r = self.ptr
        if not r:
            return None
        return METADATA_TYPES[r.type]

    @property
    def key(self):
        """header key (the part before '=', in FILTER/INFO/FORMAT/contig/fileformat etc.)"""
        cdef bcf_hrec_t *r = self.ptr
        return bcf_str_cache_get_charptr(r.key) if r and r.key else None

    @property
    def value(self):
        """header value.  Set only for generic lines, None for FILTER/INFO, etc."""
        cdef bcf_hrec_t *r = self.ptr
        return charptr_to_str(r.value) if r and r.value else None

    @property
    def attrs(self):
        """sequence of additional header attributes"""
        cdef bcf_hrec_t *r = self.ptr
        if not r:
            return ()
        cdef int i
        return tuple((bcf_str_cache_get_charptr(r.keys[i]) if r.keys[i] else None,
                      charptr_to_str(r.vals[i]) if r.vals[i] else None)
                     for i in range(r.nkeys))

    def __len__(self):
        cdef bcf_hrec_t *r = self.ptr
        return r.nkeys if r else 0

    def __bool__(self):
        cdef bcf_hrec_t *r = self.ptr
        return r != NULL and r.nkeys != 0

    def __getitem__(self, key):
        """get attribute value"""
        cdef bcf_hrec_t *r = self.ptr
        cdef int i
        if r:
            bkey = force_bytes(key)
            for i in range(r.nkeys):
                if r.keys[i] and r.keys[i] == bkey:
                    return charptr_to_str(r.vals[i]) if r.vals[i] else None
        raise KeyError('cannot find metadata key')

    def __iter__(self):
        cdef bcf_hrec_t *r = self.ptr
        if not r:
            return
        cdef int i
        for i in range(r.nkeys):
            if r.keys[i]:
                yield bcf_str_cache_get_charptr(r.keys[i])

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        try:
            self[key]
        except KeyError:
            return False
        else:
            return True

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

    def itervalues(self):
        """D.itervalues() -> an iterator over the values of D"""
        cdef bcf_hrec_t *r = self.ptr
        if not r:
            return
        cdef int i
        for i in range(r.nkeys):
            if r.keys[i]:
                yield charptr_to_str(r.vals[i]) if r.vals[i] else None

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        cdef bcf_hrec_t *r = self.ptr
        if not r:
            return
        cdef int i
        for i in range(r.nkeys):
            if r.keys[i]:
                yield (bcf_str_cache_get_charptr(r.keys[i]), charptr_to_str(r.vals[i]) if r.vals[i] else None)

    def keys(self):
        """D.keys() -> list of D's keys"""
        return list(self)

    def items(self):
        """D.items() -> list of D's (key, value) pairs, as 2-tuples"""
        return list(self.iteritems())

    def values(self):
        """D.values() -> list of D's values"""
        return list(self.itervalues())

    def update(self, items=None, **kwargs):
        """D.update([E, ]**F) -> None.

        Update D from dict/iterable E and F.
        """
        for k, v in items.items():
            self[k] = v

        if kwargs:
            for k, v in kwargs.items():
                self[k] = v

    def pop(self, key, default=_nothing):
        try:
            value = self[key]
            del self[key]
            return value
        except KeyError:
            if default is not _nothing:
                return default
            raise

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None

    #TODO: implement __richcmp__

    def __str__(self):
        cdef bcf_hrec_t *r = self.ptr

        if not r:
            raise ValueError('cannot convert deleted record to str')

        cdef kstring_t hrec_str
        hrec_str.l = hrec_str.m = 0
        hrec_str.s = NULL

        bcf_hrec_format(r, &hrec_str)

        ret = charptr_to_str_w_len(hrec_str.s, hrec_str.l)

        if hrec_str.m:
            free(hrec_str.s)

        return ret

    # FIXME: Not safe -- causes trivial segfaults at the moment
    def remove(self):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef bcf_hrec_t *r = self.ptr
        if not r:
            return
        assert r.key
        cdef char *key = r.key if r.type == BCF_HL_GEN else r.value
        bcf_hdr_remove(hdr, r.type, key)
        self.ptr = NULL


cdef VariantHeaderRecord makeVariantHeaderRecord(VariantHeader header, bcf_hrec_t *hdr):
    if not header:
        raise ValueError('invalid VariantHeader')

    if not hdr:
        return None

    cdef VariantHeaderRecord record = VariantHeaderRecord.__new__(VariantHeaderRecord)
    record.header = header
    record.ptr = hdr

    return record


cdef class VariantHeaderRecords(object):
    """sequence of :class:`VariantHeaderRecord` object from a :class:`VariantHeader` object"""
    def __init__(self, *args, **kwargs):
        raise TypeError('this class cannot be instantiated from Python')

    def __len__(self):
        return self.header.ptr.nhrec

    def __bool__(self):
        return self.header.ptr.nhrec != 0

    def __getitem__(self, index):
        cdef int32_t i = index
        if i < 0 or i >= self.header.ptr.nhrec:
            raise IndexError('invalid header record index')
        return makeVariantHeaderRecord(self.header, self.header.ptr.hrec[i])

    def __iter__(self):
        cdef int32_t i
        for i in range(self.header.ptr.nhrec):
            yield makeVariantHeaderRecord(self.header, self.header.ptr.hrec[i])

    __hash__ = None


cdef VariantHeaderRecords makeVariantHeaderRecords(VariantHeader header):
    if not header:
        raise ValueError('invalid VariantHeader')

    cdef VariantHeaderRecords records = VariantHeaderRecords.__new__(VariantHeaderRecords)
    records.header = header
    return records


cdef class VariantMetadata(object):
    """filter, info or format metadata record from a :class:`VariantHeader` object"""
    def __init__(self, *args, **kwargs):
        raise TypeError('this class cannot be instantiated from Python')

    @property
    def name(self):
        """metadata name"""
        cdef bcf_hdr_t *hdr = self.header.ptr
        return bcf_str_cache_get_charptr(hdr.id[BCF_DT_ID][self.id].key)

    # Q: Should this be exposed?
    @property
    def id(self):
        """metadata internal header id number"""
        return self.id

    @property
    def number(self):
        """metadata number (i.e. cardinality)"""
        cdef bcf_hdr_t *hdr = self.header.ptr

        if not check_header_id(hdr, self.type, self.id):
            raise ValueError('Invalid header id')

        if self.type == BCF_HL_FLT:
            return None

        cdef int l = bcf_hdr_id2length(hdr, self.type, self.id)
        if l == BCF_VL_FIXED:
            return bcf_hdr_id2number(hdr, self.type, self.id)
        elif l == BCF_VL_VAR:
            return '.'
        else:
            return METADATA_LENGTHS[l]

    @property
    def type(self):
        """metadata value type"""
        cdef bcf_hdr_t *hdr = self.header.ptr
        if not check_header_id(hdr, self.type, self.id):
            raise ValueError('Invalid header id')

        if self.type == BCF_HL_FLT:
            return None
        return VALUE_TYPES[bcf_hdr_id2type(hdr, self.type, self.id)]

    @property
    def description(self):
        """metadata description (or None if not set)"""
        descr = self.record.get('Description')
        if descr:
            descr = descr.strip('"')
        return force_str(descr)

    @property
    def record(self):
        """:class:`VariantHeaderRecord` associated with this :class:`VariantMetadata` object"""
        cdef bcf_hdr_t *hdr = self.header.ptr
        if not check_header_id(hdr, self.type, self.id):
            raise ValueError('Invalid header id')
        cdef bcf_hrec_t *hrec = hdr.id[BCF_DT_ID][self.id].val.hrec[self.type]
        if not hrec:
            return None
        return makeVariantHeaderRecord(self.header, hrec)

    def remove_header(self):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef const char *key = hdr.id[BCF_DT_ID][self.id].key
        bcf_hdr_remove(hdr, self.type, key)


cdef VariantMetadata makeVariantMetadata(VariantHeader header, int type, int id):
    if not header:
        raise ValueError('invalid VariantHeader')

    if type != BCF_HL_FLT and type != BCF_HL_INFO and type != BCF_HL_FMT:
        raise ValueError('invalid metadata type')

    if id < 0 or id >= header.ptr.n[BCF_DT_ID]:
        raise ValueError('invalid metadata id')

    cdef VariantMetadata meta = VariantMetadata.__new__(VariantMetadata)
    meta.header = header
    meta.type = type
    meta.id = id

    return meta


cdef class VariantHeaderMetadata(object):
    """mapping from filter, info or format name to :class:`VariantMetadata` object"""
    def __init__(self, *args, **kwargs):
        raise TypeError('this class cannot be instantiated from Python')

    def add(self, id, number, type, description, **kwargs):
        """Add a new filter, info or format record"""
        if id in self:
            raise ValueError('Header already exists for id={}'.format(id))

        if self.type == BCF_HL_FLT:
            if number is not None:
                raise ValueError('Number must be None when adding a filter')
            if type is not None:
                raise ValueError('Type must be None when adding a filter')

            items = [('ID', unquoted_str(id)), ('Description', description)]
        else:
            if type not in VALUE_TYPES:
                raise ValueError('unknown type specified: {}'.format(type))
            if number is None:
                number = '.'

            items = [('ID', unquoted_str(id)),
                     ('Number', unquoted_str(number)),
                     ('Type', unquoted_str(type)),
                     ('Description', description)]

        items += kwargs.items()
        self.header.add_meta(METADATA_TYPES[self.type], items=items)

    def __len__(self):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef bcf_idpair_t *idpair
        cdef int32_t i, n = 0

        for i in range(hdr.n[BCF_DT_ID]):
            idpair = hdr.id[BCF_DT_ID] + i
            if idpair.key and idpair.val and idpair.val.info[self.type] & 0xF != 0xF:
                n += 1
        return n

    def __bool__(self):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef bcf_idpair_t *idpair
        cdef int32_t i

        for i in range(hdr.n[BCF_DT_ID]):
            idpair = hdr.id[BCF_DT_ID] + i
            if idpair.key and idpair.val and idpair.val.info[self.type] & 0xF != 0xF:
                return True
        return False

    def __getitem__(self, key):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef vdict_t *d = <vdict_t *>hdr.dict[BCF_DT_ID]

        cdef bytes bkey = force_bytes(key)
        cdef khiter_t k = kh_get_vdict(d, bkey)

        if k == kh_end(d) or kh_val_vdict(d, k).info[self.type] & 0xF == 0xF:
            raise KeyError('invalid key: {}'.format(key))

        return makeVariantMetadata(self.header, self.type, kh_val_vdict(d, k).id)

    def remove_header(self, key):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef vdict_t *d = <vdict_t *>hdr.dict[BCF_DT_ID]

        cdef bytes bkey = force_bytes(key)
        cdef khiter_t k = kh_get_vdict(d, bkey)

        if k == kh_end(d) or kh_val_vdict(d, k).info[self.type] & 0xF == 0xF:
            raise KeyError('invalid key: {}'.format(key))

        bcf_hdr_remove(hdr, self.type, bkey)
        #bcf_hdr_sync(hdr)

    def clear_header(self):
        cdef bcf_hdr_t *hdr = self.header.ptr
        bcf_hdr_remove(hdr, self.type, NULL)
        #bcf_hdr_sync(hdr)

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef bcf_idpair_t *idpair
        cdef int32_t i

        for i in range(hdr.n[BCF_DT_ID]):
            idpair = hdr.id[BCF_DT_ID] + i
            if idpair.key and idpair.val and idpair.val.info[self.type] & 0xF != 0xF:
                yield bcf_str_cache_get_charptr(idpair.key)

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        try:
            self[key]
        except KeyError:
            return False
        else:
            return True

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

    def itervalues(self):
        """D.itervalues() -> an iterator over the values of D"""
        for key in self:
            yield self[key]

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        for key in self:
            yield (key, self[key])

    def keys(self):
        """D.keys() -> list of D's keys"""
        return list(self)

    def items(self):
        """D.items() -> list of D's (key, value) pairs, as 2-tuples"""
        return list(self.iteritems())

    def values(self):
        """D.values() -> list of D's values"""
        return list(self.itervalues())

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None

    #TODO: implement __richcmp__


cdef VariantHeaderMetadata makeVariantHeaderMetadata(VariantHeader header, int32_t type):
    if not header:
        raise ValueError('invalid VariantHeader')

    cdef VariantHeaderMetadata meta = VariantHeaderMetadata.__new__(VariantHeaderMetadata)
    meta.header = header
    meta.type = type

    return meta


cdef class VariantContig(object):
    """contig metadata from a :class:`VariantHeader`"""
    def __init__(self, *args, **kwargs):
        raise TypeError('this class cannot be instantiated from Python')

    @property
    def name(self):
        """contig name"""
        cdef bcf_hdr_t *hdr = self.header.ptr
        return bcf_str_cache_get_charptr(hdr.id[BCF_DT_CTG][self.id].key)

    @property
    def id(self):
        """contig internal id number"""
        return self.id

    @property
    def length(self):
        """contig length or None if not available"""
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef uint32_t length = hdr.id[BCF_DT_CTG][self.id].val.info[0]
        return length if length else None

    @property
    def header_record(self):
        """:class:`VariantHeaderRecord` associated with this :class:`VariantContig` object"""
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef bcf_hrec_t *hrec = hdr.id[BCF_DT_CTG][self.id].val.hrec[0]
        return makeVariantHeaderRecord(self.header, hrec)

    def remove_header(self):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef const char *key = hdr.id[BCF_DT_CTG][self.id].key
        bcf_hdr_remove(hdr, BCF_HL_CTG, key)


cdef VariantContig makeVariantContig(VariantHeader header, int id):
    if not header:
        raise ValueError('invalid VariantHeader')

    if id < 0 or id >= header.ptr.n[BCF_DT_CTG]:
        raise ValueError('invalid contig id')

    cdef VariantContig contig = VariantContig.__new__(VariantContig)
    contig.header = header
    contig.id = id

    return contig


cdef class VariantHeaderContigs(object):
    """mapping from contig name or index to :class:`VariantContig` object."""
    def __init__(self, *args, **kwargs):
        raise TypeError('this class cannot be instantiated from Python')

    def __len__(self):
        cdef bcf_hdr_t *hdr = self.header.ptr
        assert kh_size(<vdict_t *>hdr.dict[BCF_DT_CTG]) == hdr.n[BCF_DT_CTG]
        return hdr.n[BCF_DT_CTG]

    def __bool__(self):
        cdef bcf_hdr_t *hdr = self.header.ptr
        assert kh_size(<vdict_t *>hdr.dict[BCF_DT_CTG]) == hdr.n[BCF_DT_CTG]
        return hdr.n[BCF_DT_CTG] != 0

    def __getitem__(self, key):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef int index

        if isinstance(key, int):
            index = key
            if index < 0 or index >= hdr.n[BCF_DT_CTG]:
                raise IndexError('invalid contig index')
            return makeVariantContig(self.header, index)

        cdef vdict_t *d = <vdict_t *>hdr.dict[BCF_DT_CTG]
        cdef bytes bkey = force_bytes(key)
        cdef khiter_t k = kh_get_vdict(d, bkey)

        if k == kh_end(d):
            raise KeyError('invalid contig: {}'.format(key))

        cdef int id = kh_val_vdict(d, k).id

        return makeVariantContig(self.header, id)

    def remove_header(self, key):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef int index
        cdef const char *ckey
        cdef vdict_t *d
        cdef khiter_t k

        if isinstance(key, int):
            index = key
            if index < 0 or index >= hdr.n[BCF_DT_CTG]:
                raise IndexError('invalid contig index')
            ckey = hdr.id[BCF_DT_CTG][self.id].key
        else:
            d = <vdict_t *>hdr.dict[BCF_DT_CTG]
            key = force_bytes(key)
            if kh_get_vdict(d, key) == kh_end(d):
                raise KeyError('invalid contig: {}'.format(key))
            ckey = key

        bcf_hdr_remove(hdr, BCF_HL_CTG, ckey)

    def clear_header(self):
        cdef bcf_hdr_t *hdr = self.header.ptr
        bcf_hdr_remove(hdr, BCF_HL_CTG, NULL)
        #bcf_hdr_sync(hdr)

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef vdict_t *d = <vdict_t *>hdr.dict[BCF_DT_CTG]
        cdef uint32_t n = kh_size(d)

        assert n == hdr.n[BCF_DT_CTG]

        for i in range(n):
            yield bcf_str_cache_get_charptr(bcf_hdr_id2name(hdr, i))

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        try:
            self[key]
        except KeyError:
            return False
        else:
            return True

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

    def itervalues(self):
        """D.itervalues() -> an iterator over the values of D"""
        for key in self:
            yield self[key]

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        for key in self:
            yield (key, self[key])

    def keys(self):
        """D.keys() -> list of D's keys"""
        return list(self)

    def items(self):
        """D.items() -> list of D's (key, value) pairs, as 2-tuples"""
        return list(self.iteritems())

    def values(self):
        """D.values() -> list of D's values"""
        return list(self.itervalues())

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None

    #TODO: implement __richcmp__

    def add(self, id, length=None, **kwargs):
        """Add a new contig record"""
        if id in self:
            raise ValueError('Header already exists for contig {}'.format(id))

        items = [('ID', unquoted_str(id))]
        if length is not None:
            items.append(("length", unquoted_str(length)))
        items += kwargs.items()
        self.header.add_meta('contig', items=items)


cdef VariantHeaderContigs makeVariantHeaderContigs(VariantHeader header):
    if not header:
        raise ValueError('invalid VariantHeader')

    cdef VariantHeaderContigs contigs = VariantHeaderContigs.__new__(VariantHeaderContigs)
    contigs.header = header

    return contigs


cdef class VariantHeaderSamples(object):
    """sequence of sample names from a :class:`VariantHeader` object"""
    def __init__(self, *args, **kwargs):
        raise TypeError('this class cannot be instantiated from Python')

    def __len__(self):
        return bcf_hdr_nsamples(self.header.ptr)

    def __bool__(self):
        return bcf_hdr_nsamples(self.header.ptr) != 0

    def __getitem__(self, index):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef int32_t n = bcf_hdr_nsamples(hdr)
        cdef int32_t i = index

        if i < 0 or i >= n:
            raise IndexError('invalid sample index')

        return charptr_to_str(hdr.samples[i])

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef int32_t i, n = bcf_hdr_nsamples(hdr)

        for i in range(n):
            yield charptr_to_str(hdr.samples[i])

    def __contains__(self, key):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef vdict_t *d = <vdict_t *>hdr.dict[BCF_DT_SAMPLE]
        cdef bytes bkey = force_bytes(key)
        cdef khiter_t k = kh_get_vdict(d, bkey)

        return k != kh_end(d)

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None

    #TODO: implement __richcmp__

    def add(self, name):
        """Add a new sample"""
        self.header.add_sample(name)


cdef VariantHeaderSamples makeVariantHeaderSamples(VariantHeader header):
    if not header:
        raise ValueError('invalid VariantHeader')

    cdef VariantHeaderSamples samples = VariantHeaderSamples.__new__(VariantHeaderSamples)
    samples.header = header

    return samples


cdef class VariantHeader(object):
    """header information for a :class:`VariantFile` object"""
    #FIXME: Add structured proxy
    #FIXME: Add generic proxy
    #FIXME: Add mutable methods

    # See makeVariantHeader for C constructor
    def __cinit__(self):
        self.ptr = NULL

    # Python constructor
    def __init__(self):
        self.ptr = bcf_hdr_init(b'w')
        if not self.ptr:
            raise ValueError('cannot create VariantHeader')

    def __dealloc__(self):
        if self.ptr:
            bcf_hdr_destroy(self.ptr)
            self.ptr = NULL

    def __bool__(self):
        return self.ptr != NULL

    def copy(self):
        return makeVariantHeader(bcf_hdr_dup(self.ptr))

    def merge(self, VariantHeader header):
        if header is None:
            raise ValueError('header must not be None')
        bcf_hdr_merge(self.ptr, header.ptr)

    @property
    def version(self):
        """VCF version"""
        return force_str(bcf_hdr_get_version(self.ptr))

    @property
    def samples(self):
        """samples (:class:`VariantHeaderSamples`)"""
        return makeVariantHeaderSamples(self)

    @property
    def records(self):
        """header records (:class:`VariantHeaderRecords`)"""
        return makeVariantHeaderRecords(self)

    @property
    def contigs(self):
        """contig information (:class:`VariantHeaderContigs`)"""
        return makeVariantHeaderContigs(self)

    @property
    def filters(self):
        """filter metadata (:class:`VariantHeaderMetadata`)"""
        return makeVariantHeaderMetadata(self, BCF_HL_FLT)

    @property
    def info(self):
        """info metadata (:class:`VariantHeaderMetadata`)"""
        return makeVariantHeaderMetadata(self, BCF_HL_INFO)

    @property
    def formats(self):
        """format metadata (:class:`VariantHeaderMetadata`)"""
        return makeVariantHeaderMetadata(self, BCF_HL_FMT)

    @property
    def alts(self):
        """alt metadata (:class:`dict` ID->record).

        The data returned just a snapshot of alt records, is created
        every time the property is requested, and modifications will
        not be reflected in the header metadata and vice versa.

        i.e. it is just a dict that reflects the state of alt records
        at the time it is created.
        """
        return {record['ID']:record for record in self.records
                if record.key.upper() == 'ALT' }

    # only safe to do when opening an htsfile
    cdef _subset_samples(self, include_samples):
        keep_samples    = set(self.samples)
        include_samples = set(include_samples)
        missing_samples = include_samples - keep_samples
        keep_samples   &= include_samples

        if missing_samples:
            # FIXME: add specialized exception with payload
            raise ValueError(
                'missing {:d} requested samples'.format(
                    len(missing_samples)))

        keep_samples = force_bytes(','.join(keep_samples))
        cdef char *keep = <char *>keep_samples if keep_samples else NULL
        cdef ret = bcf_hdr_set_samples(self.ptr, keep, 0)

        if ret != 0:
            raise ValueError(
                'bcf_hdr_set_samples failed: ret = {}'.format(ret))

    def __str__(self):
        cdef int hlen
        cdef kstring_t line
        line.l = line.m = 0
        line.s = NULL

        if bcf_hdr_format(self.ptr, 0, &line) < 0:
            if line.m:
                free(line.s)
            raise ValueError('bcf_hdr_format failed')

        ret = charptr_to_str_w_len(line.s, line.l)

        if line.m:
            free(line.s)
        return ret

    def new_record(self, contig=None, start=0, stop=0, alleles=None,
                         id=None, qual=None, filter=None, info=None, samples=None,
                         **kwargs):
        """Create a new empty VariantRecord.

        Arguments are currently experimental.  Use with caution and expect
        changes in upcoming releases.

        """
        rec = makeVariantRecord(self, bcf_init())

        if not rec:
            raise MemoryError('unable to allocate BCF record')

        rec.ptr.n_sample = bcf_hdr_nsamples(self.ptr)

        if contig is not None:
            rec.contig  = contig

        rec.start = start
        rec.stop  = stop
        rec.id    = id
        rec.qual  = qual

        if alleles is not None:
            rec.alleles = alleles

        if filter is not None:
            if isinstance(filter, (list, tuple, VariantRecordFilter)):
                for f in filter:
                    rec.filter.add(f)
            else:
                rec.filter.add(filter)

        if info:
            rec.info.update(info)

        if kwargs:
            if 'GT' in kwargs:
                rec.samples[0]['GT'] = kwargs.pop('GT')
            rec.samples[0].update(kwargs)

        if samples:
            for i, sample in enumerate(samples):
                if 'GT' in sample:
                    rec.samples[i]['GT'] = sample.pop('GT')
                rec.samples[i].update(sample)

        return rec

    def add_record(self, VariantHeaderRecord record):
        """Add an existing :class:`VariantHeaderRecord` to this header"""
        if record is None:
            raise ValueError('record must not be None')

        cdef bcf_hrec_t *hrec = bcf_hrec_dup(record.ptr)

        bcf_hdr_add_hrec(self.ptr, hrec)

        self._hdr_sync()

    def add_line(self, line):
        """Add a metadata line to this header"""
        bline = force_bytes(line)
        if bcf_hdr_append(self.ptr, bline) < 0:
            raise ValueError('invalid header line')

        self._hdr_sync()


    def add_meta(self, key, value=None, items=None):
        """Add metadata to this header"""
        if not ((value is not None) ^ (items is not None)):
            raise ValueError('either value or items must be specified')

        cdef bcf_hrec_t *hrec = <bcf_hrec_t*>calloc(1, sizeof(bcf_hrec_t))
        cdef int quoted

        try:
            key = force_bytes(key)
            hrec.key = strdup(key)

            if value is not None:
                hrec.value = strdup(force_bytes(value))
            else:
                for key, value in items:
                    quoted = not isinstance(value, unquoted_str) and key not in ("ID", "Number", "Type")

                    key = force_bytes(key)
                    bcf_hrec_add_key(hrec, key, <int>len(key))

                    value = force_bytes(str(value))
                    bcf_hrec_set_val(hrec, hrec.nkeys-1, value, <int>len(value), quoted)
        except:
            bcf_hrec_destroy(hrec)
            raise

        bcf_hdr_add_hrec(self.ptr, hrec)

        self._hdr_sync()

    cdef _add_sample(self, name):
        bname = force_bytes(name)
        if bcf_hdr_add_sample(self.ptr, bname) < 0:
            raise ValueError('Duplicated sample name: {}'.format(name))

    cdef _hdr_sync(self):
        cdef bcf_hdr_t *hdr = self.ptr
        if hdr.dirty:
            if bcf_hdr_sync(hdr) < 0:
                raise MemoryError('unable to reallocate VariantHeader')

    def add_sample(self, name):
        """Add a new sample to this header"""
        self._add_sample(name)
        self._hdr_sync()

    def add_samples(self, *args):
        """Add several new samples to this header.
        This function takes multiple arguments, each of which may
        be either a sample name or an iterable returning sample names
        (e.g., a list of sample names).
        """
        for arg in args:
            if isinstance(arg, str):
                self._add_sample(arg)
            else:
                for name in arg:
                    self._add_sample(name)
        self._hdr_sync()


cdef VariantHeader makeVariantHeader(bcf_hdr_t *hdr):
    if not hdr:
        raise ValueError('cannot create VariantHeader')

    cdef VariantHeader header = VariantHeader.__new__(VariantHeader)
    header.ptr = hdr

    return header


cdef inline int bcf_header_get_info_id(bcf_hdr_t *hdr, key) except? -2:
    cdef vdict_t *d
    cdef khiter_t k
    cdef int info_id

    if isinstance(key, str):
        key = force_bytes(key)

    d = <vdict_t *>hdr.dict[BCF_DT_ID]
    k = kh_get_vdict(d, key)

    if k == kh_end(d) or kh_val_vdict(d, k).info[BCF_HL_INFO] & 0xF == 0xF:
        return -1

    return kh_val_vdict(d, k).id


########################################################################
########################################################################
## Variant Record objects
########################################################################

cdef class VariantRecordFilter(object):
    """Filters set on a :class:`VariantRecord` object, presented as a mapping from
       filter index or name to :class:`VariantMetadata` object"""
    def __init__(self, *args, **kwargs):
        raise TypeError('this class cannot be instantiated from Python')

    def __len__(self):
        return self.record.ptr.d.n_flt

    def __bool__(self):
        return self.record.ptr.d.n_flt != 0

    def __getitem__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int index, id
        cdef int n = r.d.n_flt

        if isinstance(key, int):
            index = key

            if index < 0 or index >= n:
                raise IndexError('invalid filter index')

            id = r.d.flt[index]
        else:
            if key == '.':
                key = 'PASS'

            bkey = force_bytes(key)
            id = bcf_hdr_id2int(hdr, BCF_DT_ID, bkey)

            if not check_header_id(hdr, BCF_HL_FLT, id) or not bcf_has_filter(hdr, r, bkey):
                raise KeyError('Invalid filter: {}'.format(key))

        return makeVariantMetadata(self.record.header, BCF_HL_FLT, id)

    def add(self, key):
        """Add a new filter"""
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int id

        if key == '.':
            key = 'PASS'

        cdef bytes bkey = force_bytes(key)
        id = bcf_hdr_id2int(hdr, BCF_DT_ID, bkey)

        if not check_header_id(hdr, BCF_HL_FLT, id):
            raise KeyError('Invalid filter: {}'.format(key))

        bcf_add_filter(hdr, r, id)

    def __delitem__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int index, id
        cdef int n = r.d.n_flt

        if isinstance(key, int):
            index = key

            if index < 0 or index >= n:
                raise IndexError('invalid filter index')

            id = r.d.flt[index]
        else:
            if key == '.':
                key = 'PASS'

            bkey = force_bytes(key)
            id = bcf_hdr_id2int(hdr, BCF_DT_ID, bkey)

            if not check_header_id(hdr, BCF_HL_FLT, id) or not bcf_has_filter(hdr, r, bkey):
                raise KeyError('Invalid filter: {}'.format(key))

        bcf_remove_filter(hdr, r, id, 0)

    def clear(self):
        """Clear all filters"""
        cdef bcf1_t *r = self.record.ptr
        r.d.shared_dirty |= BCF1_DIRTY_FLT
        r.d.n_flt = 0

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i

        for i in range(r.d.n_flt):
            yield bcf_str_cache_get_charptr(bcf_hdr_int2id(hdr, BCF_DT_ID, r.d.flt[i]))

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bytes bkey = force_bytes(key)
        return bcf_has_filter(hdr, r, bkey) == 1

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

    def itervalues(self):
        """D.itervalues() -> an iterator over the values of D"""
        for key in self:
            yield self[key]

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        for key in self:
            yield (key, self[key])

    def keys(self):
        """D.keys() -> list of D's keys"""
        return list(self)

    def items(self):
        """D.items() -> list of D's (key, value) pairs, as 2-tuples"""
        return list(self.iteritems())

    def values(self):
        """D.values() -> list of D's values"""
        return list(self.itervalues())

    def __richcmp__(VariantRecordFilter self not None, VariantRecordFilter other not None, int op):
        if op != 2 and op != 3:
            return NotImplemented

        cdef bcf1_t *s = self.record.ptr
        cdef bcf1_t *o = other.record.ptr

        cdef bint cmp = (s.d.n_flt == o.d.n_flt and list(self) == list(other))

        if op == 3:
            cmp = not cmp

        return cmp

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None

    #TODO: implement __richcmp__


cdef VariantRecordFilter makeVariantRecordFilter(VariantRecord record):
    if not record:
        raise ValueError('invalid VariantRecord')

    cdef VariantRecordFilter filter = VariantRecordFilter.__new__(VariantRecordFilter)
    filter.record = record

    return filter


cdef class VariantRecordFormat(object):
    """Format data present for each sample in a :class:`VariantRecord` object,
       presented as mapping from format name to :class:`VariantMetadata` object."""
    def __init__(self, *args, **kwargs):
        raise TypeError('this class cannot be instantiated from Python')

    def __len__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i, n = 0

        for i in range(r.n_fmt):
            if r.d.fmt[i].p:
                n += 1
        return n

    def __bool__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i

        for i in range(r.n_fmt):
            if r.d.fmt[i].p:
                return True
        return False

    def __getitem__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr

        cdef bytes bkey = force_bytes(key)
        cdef bcf_fmt_t *fmt = bcf_get_fmt(hdr, r, bkey)

        if not fmt or not fmt.p:
            raise KeyError('unknown format: {}'.format(key))

        return makeVariantMetadata(self.record.header, BCF_HL_FMT, fmt.id)

    def __delitem__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr

        cdef bytes bkey = force_bytes(key)
        cdef bcf_fmt_t *fmt = bcf_get_fmt(hdr, r, bkey)

        if not fmt or not fmt.p:
            raise KeyError('unknown format: {}'.format(key))

        if bcf_update_format(hdr, r, bkey, fmt.p, 0, fmt.type) < 0:
            raise ValueError('Unable to delete FORMAT')

    def clear(self):
        """Clear all formats for all samples within the associated
           :class:`VariantRecord` instance"""
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_fmt_t *fmt
        cdef const char *key
        cdef int i

        for i in reversed(range(r.n_fmt)):
            fmt = &r.d.fmt[i]
            if fmt.p:
                key = bcf_hdr_int2id(hdr, BCF_DT_ID, fmt.id)
                if bcf_update_format(hdr, r, key, fmt.p, 0, fmt.type) < 0:
                    raise ValueError('Unable to delete FORMAT')

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_fmt_t *fmt
        cdef int i

        for i in range(r.n_fmt):
            fmt = &r.d.fmt[i]
            if fmt.p:
                yield bcf_str_cache_get_charptr(bcf_hdr_int2id(hdr, BCF_DT_ID, fmt.id))

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bytes bkey = force_bytes(key)
        cdef bcf_fmt_t *fmt = bcf_get_fmt(hdr, r, bkey)
        return fmt != NULL and fmt.p != NULL

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

    def itervalues(self):
        """D.itervalues() -> an iterator over the values of D"""
        for key in self:
            yield self[key]

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        for key in self:
            yield (key, self[key])

    def keys(self):
        """D.keys() -> list of D's keys"""
        return list(self)

    def items(self):
        """D.items() -> list of D's (key, value) pairs, as 2-tuples"""
        return list(self.iteritems())

    def values(self):
        """D.values() -> list of D's values"""
        return list(self.itervalues())

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None

    #TODO: implement __richcmp__


cdef VariantRecordFormat makeVariantRecordFormat(VariantRecord record):
    if not record:
        raise ValueError('invalid VariantRecord')

    cdef VariantRecordFormat format = VariantRecordFormat.__new__(VariantRecordFormat)
    format.record = record

    return format


#TODO: Add a getmeta method to return the corresponding VariantMetadata?
cdef class VariantRecordInfo(object):
    """Info data stored in a :class:`VariantRecord` object, presented as a
       mapping from info metadata name to value."""

    def __init__(self, *args, **kwargs):
        raise TypeError('this class cannot be instantiated from Python')

    def __len__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info
        cdef const char *key
        cdef int i, count = 0

        if bcf_unpack(r, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        for i in range(r.n_info):
            info = &r.d.info[i]
            key = bcf_hdr_int2id(hdr, BCF_DT_ID, info.key)
            if info != NULL and info.vptr != NULL and strcmp(key, b'END') != 0:
                count += 1

        return count

    def __bool__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info
        cdef const char *key
        cdef int i

        if bcf_unpack(r, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        for i in range(r.n_info):
            info = &r.d.info[i]
            key = bcf_hdr_int2id(hdr, BCF_DT_ID, info.key)
            if info != NULL and info.vptr != NULL and strcmp(key, b'END') != 0:
                return True

        return False

    def __getitem__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr

        if bcf_unpack(r, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        cdef bytes bkey = force_bytes(key)

        if strcmp(bkey, b'END') == 0:
            raise KeyError('END is a reserved attribute; access is via record.stop')

        cdef bcf_info_t *info = bcf_get_info(hdr, r, bkey)

        # Cannot stop here if info == NULL, since flags must return False
        cdef int info_id = bcf_header_get_info_id(hdr, bkey) if not info else info.key

        if info_id < 0:
            raise KeyError('Unknown INFO field: {}'.format(key))

        if not check_header_id(hdr, BCF_HL_INFO, info_id):
            raise ValueError('Invalid header')

        # Handle type=Flag values
        if bcf_hdr_id2type(hdr, BCF_HL_INFO, info_id) == BCF_HT_FLAG:
            return info != NULL and info.vptr != NULL

        if not info or not info.vptr:
            raise KeyError('Invalid INFO field: {}'.format(key))

        return bcf_info_get_value(self.record, info)

    def __setitem__(self, key, value):
        cdef bytes bkey = force_bytes(key)

        if strcmp(bkey, b'END') == 0:
            raise KeyError('END is a reserved attribute; access is via record.stop')

        if bcf_unpack(self.record.ptr, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        bcf_info_set_value(self.record, key, value)

    def __delitem__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr

        cdef bytes bkey = force_bytes(key)
        if strcmp(bkey, b'END') == 0:
            raise KeyError('END is a reserved attribute; access is via record.stop')

        if bcf_unpack(r, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        cdef bcf_info_t *info = bcf_get_info(hdr, r, bkey)

        # Cannot stop here if info == NULL, since flags must return False
        cdef int info_id = bcf_header_get_info_id(hdr, bkey) if not info else info.key

        if info_id < 0:
            raise KeyError('Unknown INFO field: {}'.format(key))

        if not check_header_id(hdr, BCF_HL_INFO, info_id):
            raise ValueError('Invalid header')

        # Handle flags
        if bcf_hdr_id2type(hdr, BCF_HL_INFO, info_id) == BCF_HT_FLAG and (not info or not info.vptr):
            return

        if not info or not info.vptr:
            raise KeyError('Unknown INFO field: {}'.format(key))

        if bcf_update_info(hdr, r, bkey, NULL, 0, info.type) < 0:
            raise ValueError('Unable to delete INFO')

    def clear(self):
        """Clear all info data"""
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info
        cdef const char *key
        cdef int i

        if bcf_unpack(r, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        for i in range(r.n_info):
            info = &r.d.info[i]
            if info and info.vptr:
                key = bcf_hdr_int2id(hdr, BCF_DT_ID, info.key)
                if strcmp(key, b'END') == 0:
                    continue
                if bcf_update_info(hdr, r, key, NULL, 0, info.type) < 0:
                    raise ValueError('Unable to delete INFO')

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info
        cdef const char *key
        cdef int i

        if bcf_unpack(r, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        for i in range(r.n_info):
            info = &r.d.info[i]
            if info and info.vptr:
                key = bcf_hdr_int2id(hdr, BCF_DT_ID, info.key)
                if strcmp(key, b'END') != 0:
                    yield bcf_str_cache_get_charptr(key)

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr

        if bcf_unpack(r, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        cdef bytes bkey = force_bytes(key)

        if strcmp(bkey, b'END') == 0:
            return default

        cdef bcf_info_t *info = bcf_get_info(hdr, r, bkey)

        # Cannot stop here if info == NULL, since flags must return False
        cdef int info_id = bcf_header_get_info_id(hdr, bkey) if not info else info.key

        if not check_header_id(hdr, BCF_HL_INFO, info_id):
            raise ValueError('Invalid header')

        # Handle flags
        if bcf_hdr_id2type(hdr, BCF_HL_INFO, info_id) == BCF_HT_FLAG:
            return info != NULL and info.vptr != NULL

        if not info or not info.vptr:
            return default

        return bcf_info_get_value(self.record, info)

    def __contains__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr

        if bcf_unpack(r, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        cdef bytes bkey = force_bytes(key)

        if strcmp(bkey, b'END') == 0:
            return False

        cdef bcf_info_t *info = bcf_get_info(hdr, r, bkey)

        return info != NULL and info.vptr != NULL

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

    def itervalues(self):
        """D.itervalues() -> an iterator over the values of D"""
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info
        cdef const char *key
        cdef int i

        if bcf_unpack(r, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        for i in range(r.n_info):
            info = &r.d.info[i]
            if info and info.vptr:
                key = bcf_hdr_int2id(hdr, BCF_DT_ID, info.key)
                if strcmp(key, b'END') != 0:
                    yield bcf_info_get_value(self.record, info)

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info
        cdef const char *key
        cdef int i

        if bcf_unpack(r, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        for i in range(r.n_info):
            info = &r.d.info[i]
            if info and info.vptr:
                key = bcf_hdr_int2id(hdr, BCF_DT_ID, info.key)
                if strcmp(key, b'END') != 0:
                    value = bcf_info_get_value(self.record, info)
                    yield bcf_str_cache_get_charptr(key), value

    def keys(self):
        """D.keys() -> list of D's keys"""
        return list(self)

    def items(self):
        """D.items() -> list of D's (key, value) pairs, as 2-tuples"""
        return list(self.iteritems())

    def values(self):
        """D.values() -> list of D's values"""
        return list(self.itervalues())

    def update(self, items=None, **kwargs):
        """D.update([E, ]**F) -> None.

        Update D from dict/iterable E and F.
        """
        for k, v in items.items():
            if k != 'END':
                self[k] = v

        if kwargs:
            kwargs.pop('END', None)
            for k, v in kwargs.items():
                self[k] = v

    def pop(self, key, default=_nothing):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr

        if bcf_unpack(r, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        cdef bytes bkey = force_bytes(key)
        cdef bcf_info_t *info = bcf_get_info(hdr, r, bkey)

        # Cannot stop here if info == NULL, since flags must return False
        cdef int info_id = bcf_header_get_info_id(hdr, bkey) if not info else info.key

        if info_id < 0:
            if default is _nothing:
                raise KeyError('Unknown INFO field: {}'.format(key))
            return default

        if not check_header_id(hdr, BCF_HL_INFO, info_id):
            raise ValueError('Invalid header')

        # Handle flags
        if bcf_hdr_id2type(hdr, BCF_HL_INFO, info_id) == BCF_HT_FLAG and (not info or not info.vptr):
            return

        if not info or not info.vptr:
            if default is _nothing:
                raise KeyError('Unknown INFO field: {}'.format(key))
            return default

        value = bcf_info_get_value(self.record, info)

        if bcf_update_info(hdr, r, bkey, NULL, 0, info.type) < 0:
            raise ValueError('Unable to delete INFO')

        return value

    def __richcmp__(VariantRecordInfo self not None, VariantRecordInfo other not None, int op):
        if op != 2 and op != 3:
            return NotImplemented

        cdef bcf1_t *s = self.record.ptr
        cdef bcf1_t *o = other.record.ptr

        # Cannot use n_info as shortcut logic, since null values may remain
        cdef bint cmp = dict(self) == dict(other)

        if op == 3:
            cmp = not cmp

        return cmp

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None


cdef VariantRecordInfo makeVariantRecordInfo(VariantRecord record):
    if not record:
        raise ValueError('invalid VariantRecord')

    cdef VariantRecordInfo info = VariantRecordInfo.__new__(VariantRecordInfo)
    info.record = record

    return info


cdef class VariantRecordSamples(object):
    """mapping from sample index or name to :class:`VariantRecordSample` object."""
    def __init__(self, *args, **kwargs):
        raise TypeError('this class cannot be instantiated from Python')

    def __len__(self):
        return self.record.ptr.n_sample  # bcf_hdr_nsamples(self.record.header.ptr)

    def __bool__(self):
        return self.record.ptr.n_sample != 0  # bcf_hdr_nsamples(self.record.header.ptr) != 0

    def __getitem__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int n = self.record.ptr.n_sample
        cdef int sample_index
        cdef vdict_t *d
        cdef khiter_t k

        if isinstance(key, int):
            sample_index = key
        else:
            bkey = force_bytes(key)
            sample_index = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, bkey)
            if sample_index < 0:
                raise KeyError('invalid sample name: {}'.format(key))

        if sample_index < 0 or sample_index >= n:
            raise IndexError('invalid sample index')

        return makeVariantRecordSample(self.record, sample_index)

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int32_t i, n = self.record.ptr.n_sample

        for i in range(n):
            yield charptr_to_str(hdr.samples[i])

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int n = self.record.ptr.n_sample
        cdef int sample_index
        cdef vdict_t *d
        cdef khiter_t k

        if isinstance(key, int):
            sample_index = key
        else:
            bkey = force_bytes(key)
            sample_index = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, bkey)
            if sample_index < 0:
                raise KeyError('invalid sample name: {}'.format(key))

        return 0 <= sample_index < n

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

    def itervalues(self):
        """D.itervalues() -> an iterator over the values of D"""
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int32_t i, n = self.record.ptr.n_sample

        for i in range(n):
            yield makeVariantRecordSample(self.record, i)

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int32_t i, n = self.record.ptr.n_sample

        for i in range(n):
            yield (charptr_to_str(hdr.samples[i]), makeVariantRecordSample(self.record, i))

    def keys(self):
        """D.keys() -> list of D's keys"""
        return list(self)

    def items(self):
        """D.items() -> list of D's (key, value) pairs, as 2-tuples"""
        return list(self.iteritems())

    def values(self):
        """D.values() -> list of D's values"""
        return list(self.itervalues())

    def update(self, items=None, **kwargs):
        """D.update([E, ]**F) -> None.

        Update D from dict/iterable E and F.
        """
        for k, v in items.items():
            self[k] = v

        if kwargs:
            for k, v in kwargs.items():
                self[k] = v

    def pop(self, key, default=_nothing):
        try:
            value = self[key]
            del self[key]
            return value
        except KeyError:
            if default is not _nothing:
                return default
            raise

    def __richcmp__(VariantRecordSamples self not None, VariantRecordSamples other not None, int op):
        if op != 2 and op != 3:
            return NotImplemented

        cdef bcf1_t *s = self.record.ptr
        cdef bcf1_t *o = other.record.ptr

        cdef bint cmp = (s.n_sample == o.n_sample and self.values() == other.values())

        if op == 3:
            cmp = not cmp

        return cmp

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None


cdef VariantRecordSamples makeVariantRecordSamples(VariantRecord record):
    if not record:
        raise ValueError('invalid VariantRecord')

    cdef VariantRecordSamples samples = VariantRecordSamples.__new__(
        VariantRecordSamples)
    samples.record = record

    return samples


cdef class VariantRecord(object):
    """Variant record"""
    def __init__(self, *args, **kwargs):
        raise TypeError('this class cannot be instantiated from Python')

    def __dealloc__(self):
        if self.ptr:
            bcf_destroy1(self.ptr)
            self.ptr = NULL

    def copy(self):
        """return a copy of this VariantRecord object"""
        return makeVariantRecord(self.header, bcf_dup(self.ptr))

    def translate(self, VariantHeader dst_header):
        if dst_header is None:
            raise ValueError('dst_header must not be None')

        cdef bcf_hdr_t *src_hdr = self.header.ptr
        cdef bcf_hdr_t *dst_hdr = dst_header.ptr

        if src_hdr != dst_hdr:
            if self.ptr.n_sample != bcf_hdr_nsamples(dst_hdr):
                msg = 'Cannot translate record.  Number of samples does not match header ({} vs {})'
                raise ValueError(msg.format(self.ptr.n_sample, bcf_hdr_nsamples(dst_hdr)))

            bcf_translate(dst_hdr, src_hdr, self.ptr)
            self.header = dst_header

    @property
    def rid(self):
        """internal reference id number"""
        return self.ptr.rid

    @rid.setter
    def rid(self, value):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef int r = value
        if r < 0 or r >= hdr.n[BCF_DT_CTG] or not hdr.id[BCF_DT_CTG][r].val:
            raise ValueError('invalid reference id')
        self.ptr.rid = r

    @property
    def chrom(self):
        """chromosome/contig name"""
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef int rid = self.ptr.rid
        if rid < 0 or rid >= hdr.n[BCF_DT_CTG]:
            raise ValueError('Invalid header')
        return bcf_str_cache_get_charptr(bcf_hdr_id2name(hdr, rid))

    @chrom.setter
    def chrom(self, value):
        cdef vdict_t *d = <vdict_t*>self.header.ptr.dict[BCF_DT_CTG]
        bchrom = force_bytes(value)
        cdef khint_t k = kh_get_vdict(d, bchrom)
        if k == kh_end(d):
            raise ValueError('Invalid chromosome/contig')
        self.ptr.rid = kh_val_vdict(d, k).id

    @property
    def contig(self):
        """chromosome/contig name"""
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef int rid = self.ptr.rid
        if rid < 0 or rid >= hdr.n[BCF_DT_CTG]:
            raise ValueError('Invalid header')
        return bcf_str_cache_get_charptr(bcf_hdr_id2name(hdr, rid))

    @contig.setter
    def contig(self, value):
        cdef vdict_t *d = <vdict_t*>self.header.ptr.dict[BCF_DT_CTG]
        bchrom = force_bytes(value)
        cdef khint_t k = kh_get_vdict(d, bchrom)
        if k == kh_end(d):
            raise ValueError('Invalid chromosome/contig')
        self.ptr.rid = kh_val_vdict(d, k).id

    @property
    def pos(self):
        """record start position on chrom/contig (1-based inclusive)"""
        return self.ptr.pos + 1

    @pos.setter
    def pos(self, value):
        cdef int p = value
        if p < 1:
            raise ValueError('Position must be positive')
        self.ptr.pos = p - 1
        bcf_sync_end(self)

    @property
    def start(self):
        """record start position on chrom/contig (0-based inclusive)"""
        return self.ptr.pos

    @start.setter
    def start(self, value):
        cdef int s = value
        if s < 0:
            raise ValueError('Start coordinate must be non-negative')
        self.ptr.pos = s
        bcf_sync_end(self)

    @property
    def stop(self):
        """record stop position on chrom/contig (0-based exclusive)"""
        return self.ptr.pos + self.ptr.rlen

    @stop.setter
    def stop(self, value):
        cdef int s = value
        if s < 0:
            raise ValueError('Stop coordinate must be non-negative')
        self.ptr.rlen = s - self.ptr.pos
        bcf_sync_end(self)

    @property
    def rlen(self):
        """record length on chrom/contig (aka rec.stop - rec.start)"""
        return self.ptr.rlen

    @rlen.setter
    def rlen(self, value):
        cdef int r = value
        self.ptr.rlen = r
        bcf_sync_end(self)

    @property
    def qual(self):
        """phred scaled quality score or None if not available"""
        return self.ptr.qual if not bcf_float_is_missing(self.ptr.qual) else None

    @qual.setter
    def qual(self, value):
        if value is not None:
            self.ptr.qual = value
        else:
            bcf_float_set(&self.ptr.qual, bcf_float_missing)


#   @property
#   def n_allele(self):
#       return self.ptr.n_allele

#   @property
#   def n_sample(self):
#       return self.ptr.n_sample

    @property
    def id(self):
        """record identifier or None if not available"""
        cdef bcf1_t *r = self.ptr
        if bcf_unpack(r, BCF_UN_STR) < 0:
            raise ValueError('Error unpacking VariantRecord')
        # causes a memory leak https://github.com/pysam-developers/pysam/issues/773
        # return bcf_str_cache_get_charptr(r.d.id) if r.d.id != b'.' else None
        if (r.d.m_id == 0):
            raise ValueError('Error extracting ID')
        return charptr_to_str(r.d.id) if r.d.id != b'.' else None

    @id.setter
    def id(self, value):
        cdef bcf1_t *r = self.ptr
        if bcf_unpack(r, BCF_UN_STR) < 0:
            raise ValueError('Error unpacking VariantRecord')
        cdef char *idstr = NULL
        if value is not None:
            bid = force_bytes(value)
            idstr = bid
        if bcf_update_id(self.header.ptr, self.ptr, idstr) < 0:
            raise ValueError('Error updating id')

    @property
    def ref(self):
        """reference allele"""
        cdef bcf1_t *r = self.ptr
        if bcf_unpack(r, BCF_UN_STR) < 0:
            raise ValueError('Error unpacking VariantRecord')
        return charptr_to_str(r.d.allele[0]) if r.d.allele else None

    @ref.setter
    def ref(self, value):
        cdef bcf1_t *r = self.ptr
        if bcf_unpack(r, BCF_UN_STR) < 0:
            raise ValueError('Error unpacking VariantRecord')
        #FIXME: Set alleles directly -- this is stupid
        if not value:
            raise ValueError('ref allele must not be null')
        value = force_bytes(value)
        if r.d.allele and r.n_allele:
            alleles = [r.d.allele[i] for i in range(r.n_allele)]
            alleles[0] = value
        else:
            alleles = [value, '<NON_REF>']
        self.alleles = alleles
        bcf_sync_end(self)

    @property
    def alleles(self):
        """tuple of reference allele followed by alt alleles"""
        cdef bcf1_t *r = self.ptr
        if bcf_unpack(r, BCF_UN_STR) < 0:
            raise ValueError('Error unpacking VariantRecord')
        if not r.d.allele:
            return None
        cdef tuple res = PyTuple_New(r.n_allele)
        for i in range(r.n_allele):
            a = charptr_to_str(r.d.allele[i])
            PyTuple_SET_ITEM(res, i, a)
            Py_INCREF(a)
        return res

    @alleles.setter
    def alleles(self, values):
        cdef bcf1_t *r = self.ptr

        # Cache rlen of symbolic alleles before call to bcf_update_alleles_str
        cdef int rlen = r.rlen

        if bcf_unpack(r, BCF_UN_STR) < 0:
            raise ValueError('Error unpacking VariantRecord')

        values = [force_bytes(v) for v in values]

        if len(values) < 2:
            raise ValueError('must set at least 2 alleles')

        if b'' in values:
            raise ValueError('cannot set null allele')

        value = b','.join(values)

        if bcf_update_alleles_str(self.header.ptr, r, value) < 0:
            raise ValueError('Error updating alleles')

        # Reset rlen if alternate allele isn't symbolic, otherwise used cached
        if has_symbolic_allele(self):
            self.ptr.rlen = rlen
        else:
            self.ptr.rlen = len(values[0])
        r.d.var_type = -1
        bcf_sync_end(self)

    @property
    def alts(self):
        """tuple of alt alleles"""
        cdef bcf1_t *r = self.ptr
        if bcf_unpack(r, BCF_UN_STR) < 0:
            raise ValueError('Error unpacking VariantRecord')
        if r.n_allele < 2 or not r.d.allele:
            return None
        cdef tuple res = PyTuple_New(r.n_allele - 1)
        for i in range(1, r.n_allele):
            a = charptr_to_str(r.d.allele[i])
            PyTuple_SET_ITEM(res, i - 1, a)
            Py_INCREF(a)
        return res

    @alts.setter
    def alts(self, value):
        #FIXME: Set alleles directly -- this is stupid
        cdef bcf1_t *r = self.ptr
        if bcf_unpack(r, BCF_UN_STR) < 0:
            raise ValueError('Error unpacking VariantRecord')
        value = [force_bytes(v) for v in value]
        if b'' in value:
            raise ValueError('cannot set null alt allele')
        ref  = [r.d.allele[0] if r.d.allele and r.n_allele else b'.']
        self.alleles = ref + value
        r.d.var_type = -1

    @property
    def filter(self):
        """filter information (see :class:`VariantRecordFilter`)"""
        if bcf_unpack(self.ptr, BCF_UN_FLT) < 0:
            raise ValueError('Error unpacking VariantRecord')
        return makeVariantRecordFilter(self)

    @property
    def info(self):
        """info data (see :class:`VariantRecordInfo`)"""
        if bcf_unpack(self.ptr, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')
        return makeVariantRecordInfo(self)

    @property
    def format(self):
        """sample format metadata (see :class:`VariantRecordFormat`)"""
        if bcf_unpack(self.ptr, BCF_UN_FMT) < 0:
            raise ValueError('Error unpacking VariantRecord')
        return makeVariantRecordFormat(self)

    @property
    def samples(self):
        """sample data (see :class:`VariantRecordSamples`)"""
        if bcf_unpack(self.ptr, BCF_UN_ALL) < 0:
            raise ValueError('Error unpacking VariantRecord')
        return makeVariantRecordSamples(self)

    property alleles_variant_types:
        def __get__(self):
            cdef bcf1_t *r = self.ptr
            cdef tuple result = PyTuple_New(r.n_allele)

            for i in range(r.n_allele):
                tp = bcf_get_variant_type(r, i)

                if tp == VCF_REF:
                    v_type = "REF"
                elif tp == VCF_SNP:
                    v_type = "SNP"
                elif tp == VCF_MNP:
                    v_type = "MNP"
                elif tp == VCF_INDEL:
                    v_type = "INDEL"
                elif tp == VCF_BND:
                    v_type = "BND"
                elif tp == VCF_OVERLAP:
                    v_type = "OVERLAP"
                else:
                    v_type = "OTHER"

                PyTuple_SET_ITEM(result, i, v_type)
                Py_INCREF(v_type)

            return result

    def __richcmp__(VariantRecord self not None, VariantRecord other not None, int op):
        if op != 2 and op != 3:
            return NotImplemented

        cdef bcf1_t *s = self.ptr
        cdef bcf1_t *o = other.ptr

        cdef bint cmp = self is other or (
                             s.pos        == o.pos
                        and  s.rlen       == o.rlen
                        and ((bcf_float_is_missing(s.qual) and bcf_float_is_missing(o.qual))
                          or s.qual       == o.qual)
                        and  s.n_sample   == o.n_sample
                        and  s.n_allele   == o.n_allele
                        and  self.contig  == other.contig
                        and  self.alleles == other.alleles
                        and  self.id      == other.id
                        and  self.info    == other.info
                        and  self.filter  == other.filter
                        and  self.samples == other.samples)

        if op == 3:
            cmp = not cmp

        return cmp

    def __str__(self):
        cdef kstring_t line
        cdef char c

        line.l = line.m = 0
        line.s = NULL

        if vcf_format(self.header.ptr, self.ptr, &line) < 0:
            if line.m:
                free(line.s)
            raise ValueError('vcf_format failed')

        # Strip CR/LF?
        #while line.l:
        #    c = line.s[line.l - 1]
        #    if c != b'\n' and c != b'\r':
        #        break
        #    line.l -= 1

        ret = charptr_to_str_w_len(line.s, line.l)

        if line.m:
            free(line.s)

        return ret


cdef VariantRecord makeVariantRecord(VariantHeader header, bcf1_t *r):
    if not header:
        raise ValueError('invalid VariantHeader')

    if not r:
        raise ValueError('cannot create VariantRecord')

    if r.errcode:
        msg = []
        #if r.errcode & BCF_ERR_CTG_UNDEF:
        #    msg.append('undefined contig')
        #if r.errcode & BCF_ERR_TAG_UNDEF:
        #    msg.append('undefined tag')
        if r.errcode & BCF_ERR_NCOLS:
            msg.append('invalid number of columns')
        if r.errcode & BCF_ERR_LIMITS:
            msg.append('limits violated')
        if r.errcode & BCF_ERR_CHAR:
            msg.append('invalid character found')
        if r.errcode & BCF_ERR_CTG_INVALID:
            msg.append('invalid contig')
        if r.errcode & BCF_ERR_TAG_INVALID:
            msg.append('invalid tag')

        if msg:
            msg = ', '.join(msg)
            raise ValueError('Error(s) reading record: {}'.format(msg))

    cdef VariantRecord record = VariantRecord.__new__(VariantRecord)
    record.header = header
    record.ptr = r

    return record


########################################################################
########################################################################
## Variant Sampletype object
########################################################################


cdef class VariantRecordSample(object):
    """Data for a single sample from a :class:`VariantRecord` object.
       Provides data accessors for genotypes and a mapping interface
       from format name to values.
    """
    def __init__(self, *args, **kwargs):
        raise TypeError('this class cannot be instantiated from Python')

    @property
    def name(self):
        """sample name"""
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int32_t n = r.n_sample

        if self.index < 0 or self.index >= n:
            raise ValueError('invalid sample index')

        return charptr_to_str(hdr.samples[self.index])

    @property
    def allele_indices(self):
        """allele indices for called genotype, if present.  Otherwise None"""
        return bcf_format_get_allele_indices(self)

    @allele_indices.setter
    def allele_indices(self, value):
        self['GT'] = value

    @allele_indices.deleter
    def allele_indices(self):
        self['GT'] = ()

    @property
    def alleles(self):
        """alleles for called genotype, if present.  Otherwise None"""
        return bcf_format_get_alleles(self)

    @alleles.setter
    def alleles(self, value: tuple):
        # Sets the genotype, supply a tuple of alleles to set.
        # The supplied alleles need to be defined in the correspoding pysam.libcbcf.VariantRecord
        # The genotype is reset when an empty tuple, None or (None,) is supplied

        if value==(None,) or value==tuple() or value is None:
            self['GT'] = ()
            return

        if any((type(x) == int for x in value)):
            raise ValueError('Use .allele_indices to set integer allele indices')

        # determine and set allele indices:    
        try:
            self['GT'] = tuple( (self.record.alleles.index(allele) for allele in value) )
        except ValueError:
            raise ValueError("One or more of the supplied sample alleles are not defined as alleles of the corresponding pysam.libcbcf.VariantRecord."
                             "First set the .alleles of this record to define the alleles")

    @alleles.deleter
    def alleles(self):
        self['GT'] = ()

    @property
    def phased(self):
        """False if genotype is missing or any allele is unphased.  Otherwise True."""
        return bcf_sample_get_phased(self)

    @phased.setter
    def phased(self, value):
        bcf_sample_set_phased(self, value)

    def __len__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i, n = 0

        if bcf_unpack(r, BCF_UN_FMT) < 0:
            raise ValueError('Error unpacking VariantRecord')

        for i in range(r.n_fmt):
            if r.d.fmt[i].p:
                n += 1
        return n

    def __bool__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i

        if bcf_unpack(r, BCF_UN_FMT) < 0:
            raise ValueError('Error unpacking VariantRecord')

        for i in range(r.n_fmt):
            if r.d.fmt[i].p:
                return True
        return False

    def __getitem__(self, key):
        return bcf_format_get_value(self, key)

    def __setitem__(self, key, value):
        bcf_format_set_value(self, key, value)

    def __delitem__(self, key):
        bcf_format_del_value(self, key)

    def clear(self):
        """Clear all format data (including genotype) for this sample"""
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_fmt_t *fmt
        cdef int i

        for i in range(r.n_fmt):
            fmt = &r.d.fmt[i]
            if fmt.p:
                bcf_format_del_value(self, bcf_hdr_int2id(hdr, BCF_DT_ID, fmt.id))

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_fmt_t *fmt
        cdef int i

        for i in range(r.n_fmt):
            fmt = &r.d.fmt[i]
            if r.d.fmt[i].p:
                yield bcf_str_cache_get_charptr(bcf_hdr_int2id(hdr, BCF_DT_ID, fmt.id))

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bytes bkey = force_bytes(key)
        cdef bcf_fmt_t *fmt = bcf_get_fmt(hdr, r, bkey)
        return fmt != NULL and fmt.p != NULL

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

    def itervalues(self):
        """D.itervalues() -> an iterator over the values of D"""
        for key in self:
            yield self[key]

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        for key in self:
            yield (key, self[key])

    def keys(self):
        """D.keys() -> list of D's keys"""
        return list(self)

    def items(self):
        """D.items() -> list of D's (key, value) pairs, as 2-tuples"""
        return list(self.iteritems())

    def values(self):
        """D.values() -> list of D's values"""
        return list(self.itervalues())

    def update(self, items=None, **kwargs):
        """D.update([E, ]**F) -> None.

        Update D from dict/iterable E and F.
        """
        for k, v in items.items():
            self[k] = v

        if kwargs:
            for k, v in kwargs.items():
                self[k] = v

    def pop(self, key, default=_nothing):
        try:
            value = self[key]
            del self[key]
            return value
        except KeyError:
            if default is not _nothing:
                return default
            raise

    def __richcmp__(VariantRecordSample self not None, VariantRecordSample other not None, int op):
        if op != 2 and op != 3:
            return NotImplemented

        cdef bint cmp = dict(self) == dict(other)

        if op == 3:
            cmp = not cmp

        return cmp

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None


cdef VariantRecordSample makeVariantRecordSample(VariantRecord record, int32_t sample_index):
    if not record or sample_index < 0:
        raise ValueError('cannot create VariantRecordSample')

    cdef VariantRecordSample sample = VariantRecordSample.__new__(VariantRecordSample)
    sample.record = record
    sample.index = sample_index

    return sample


########################################################################
########################################################################
## Index objects
########################################################################


cdef class BaseIndex(object):
    def __init__(self):
        self.refs = ()
        self.remap = {}

    def __len__(self):
        return len(self.refs)

    def __bool__(self):
        return len(self.refs) != 0

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.refs[key]
        else:
            return self.refmap[key]

    def __iter__(self):
        return iter(self.refs)

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        try:
            self[key]
        except KeyError:
            return False
        else:
            return True

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

    def itervalues(self):
        """D.itervalues() -> an iterator over the values of D"""
        for key in self:
            yield self[key]

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        for key in self:
            yield (key, self[key])

    def keys(self):
        """D.keys() -> list of D's keys"""
        return list(self)

    def items(self):
        """D.items() -> list of D's (key, value) pairs, as 2-tuples"""
        return list(self.iteritems())

    def values(self):
        """D.values() -> list of D's values"""
        return list(self.itervalues())

    def update(self, items=None, **kwargs):
        """D.update([E, ]**F) -> None.

        Update D from dict/iterable E and F.
        """
        for k, v in items.items():
            self[k] = v

        if kwargs:
            for k, v in kwargs.items():
                self[k] = v

    def pop(self, key, default=_nothing):
        try:
            value = self[key]
            del self[key]
            return value
        except KeyError:
            if default is not _nothing:
                return default
            raise

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None

    #TODO: implement __richcmp__


cdef class BCFIndex(object):
    """CSI index data structure for BCF files"""
    def __init__(self):
        self.refs = ()
        self.refmap = {}

        if not self.ptr:
            raise ValueError('Invalid index object')

        cdef int n
        cdef const char **refs = bcf_index_seqnames(self.ptr, self.header.ptr, &n)

        self.refs = char_array_to_tuple(refs, n, free_after=1) if refs else ()
        self.refmap = { r:i for i,r in enumerate(self.refs) }

    def __dealloc__(self):
        if self.ptr:
            hts_idx_destroy(self.ptr)
            self.ptr = NULL

    def fetch(self, bcf, contig, start, stop, reopen):
        return BCFIterator(bcf, contig, start, stop, reopen)


cdef BCFIndex makeBCFIndex(VariantHeader header, hts_idx_t *idx):
    if not idx:
        return None

    if not header:
        raise ValueError('invalid VariantHeader')

    cdef BCFIndex index = BCFIndex.__new__(BCFIndex)
    index.header = header
    index.ptr = idx
    index.__init__()

    return index


cdef class TabixIndex(BaseIndex):
    """Tabix index data structure for VCF files"""
    def __init__(self):
        self.refs = ()
        self.refmap = {}

        if not self.ptr:
            raise ValueError('Invalid index object')

        cdef int n
        cdef const char **refs = tbx_seqnames(self.ptr, &n)

        self.refs = char_array_to_tuple(refs, n, free_after=1) if refs else ()
        self.refmap = { r:i for i,r in enumerate(self.refs) }

    def __dealloc__(self):
        if self.ptr:
            tbx_destroy(self.ptr)
            self.ptr = NULL

    def fetch(self, bcf, contig, start, stop, reopen):
        return TabixIterator(bcf, contig, start, stop, reopen)


cdef TabixIndex makeTabixIndex(tbx_t *idx):
    if not idx:
        return None

    cdef TabixIndex index = TabixIndex.__new__(TabixIndex)
    index.ptr = idx
    index.__init__()

    return index


########################################################################
########################################################################
## Iterators
########################################################################


cdef class BaseIterator(object):
    pass


# Internal function to clean up after iteration stop or failure.
# This would be a nested function if it weren't a cdef function.
cdef void _stop_BCFIterator(BCFIterator self, bcf1_t *record):
    bcf_destroy1(record)

    # destroy iter so future calls to __next__ raise StopIteration
    bcf_itr_destroy(self.iter)
    self.iter = NULL


cdef class BCFIterator(BaseIterator):
    def __init__(self, VariantFile bcf, contig=None, start=None, stop=None, reopen=True):
        if bcf is None:
            raise ValueError('bcf must not be None')

        if contig is None:
            raise ValueError('contig must be specified')

        if not isinstance(bcf.index, BCFIndex):
            raise ValueError('bcf index required')

        cdef BCFIndex index = bcf.index

        self.bcf = bcf
        self.index = index

        cdef int rid, cstart, cstop

        try:
            rid = index.refmap[contig]
        except KeyError:
            # A query for a non-existent contig yields an empty iterator, does not raise an error
            self.iter = NULL
            return

        if reopen:
            self.bcf = self.bcf.copy()

        cstart = start if start is not None else 0
        cstop  = stop  if stop  is not None else MAX_POS

        with nogil:
            self.iter = bcf_itr_queryi(index.ptr, rid, cstart, cstop)

        if not self.iter:
            if errno:
                raise IOError(errno, strerror(errno))
            else:
                raise IOError('unable to fetch {}:{}-{}'.format(contig, start+1, stop))

    def __dealloc__(self):
        if self.iter:
            bcf_itr_destroy(self.iter)
            self.iter = NULL

    def __iter__(self):
        return self

    def __next__(self):
        if not self.iter:
            raise StopIteration

        cdef bcf1_t *record = bcf_init1()

        if not record:
            raise MemoryError('unable to allocate BCF record')

        record.pos = -1
        if self.bcf.drop_samples:
            record.max_unpack = BCF_UN_SHR

        cdef int ret

        with nogil:
            ret = bcf_itr_next(self.bcf.htsfile, self.iter, record)

        if ret < 0:
            _stop_BCFIterator(self, record)
            if ret == -1:
                raise StopIteration
            elif ret == -2:
                raise IOError('truncated file')
            elif errno:
                raise IOError(errno, strerror(errno))
            else:
                raise IOError('unable to fetch next record')

        ret = bcf_subset_format(self.bcf.header.ptr, record)

        if ret < 0:
            _stop_BCFIterator(self, record)
            raise ValueError('error in bcf_subset_format')

        return makeVariantRecord(self.bcf.header, record)


cdef class TabixIterator(BaseIterator):
    def __cinit__(self, *args, **kwargs):
        self.line_buffer.l = 0
        self.line_buffer.m = 0
        self.line_buffer.s = NULL

    def __init__(self, VariantFile bcf, contig=None, start=None, stop=None, reopen=True):
        if bcf is None:
            raise ValueError('bcf must not be None')

        if not isinstance(bcf.index, TabixIndex):
            raise ValueError('tabix index required')

        cdef TabixIndex index = bcf.index

        self.bcf = bcf
        self.index = index

        cdef int rid, cstart, cstop

        try:
            rid = index.refmap[contig]
        except KeyError:
            # A query for a non-existent contig yields an empty iterator, does not raise an error
            self.iter = NULL
            return

        if reopen:
            self.bcf = self.bcf.copy()

        cstart = start if start is not None else 0
        cstop  = stop  if stop  is not None else MAX_POS

        self.iter = tbx_itr_queryi(index.ptr, rid, start, stop)

        if not self.iter:
            if errno:
                raise IOError(errno, strerror(errno))
            else:
                raise IOError('unable to fetch {}:{}-{}'.format(contig, start+1, stop))

    def __dealloc__(self):
        if self.iter:
            tbx_itr_destroy(self.iter)
            self.iter = NULL

        if self.line_buffer.m:
            free(self.line_buffer.s)

        self.line_buffer.l = 0
        self.line_buffer.m = 0
        self.line_buffer.s = NULL

    def __iter__(self):
        return self

    def __next__(self):
        if not self.iter:
            raise StopIteration

        cdef int ret

        with nogil:
            ret = tbx_itr_next(self.bcf.htsfile, self.index.ptr, self.iter, &self.line_buffer)

        if ret < 0:
            tbx_itr_destroy(self.iter)
            self.iter = NULL
            if ret == -1:
                raise StopIteration
            elif ret == -2:
                raise IOError('truncated file')
            elif errno:
                raise IOError(errno, strerror(errno))
            else:
                raise IOError('unable to fetch next record')

        cdef bcf1_t *record = bcf_init1()

        if not record:
            raise MemoryError('unable to allocate BCF record')

        record.pos = -1
        if self.bcf.drop_samples:
            record.max_unpack = BCF_UN_SHR

        ret = vcf_parse1(&self.line_buffer, self.bcf.header.ptr, record)

        # FIXME: stop iteration on parse failure?
        if ret < 0:
            bcf_destroy1(record)
            raise ValueError('error in vcf_parse')

        return makeVariantRecord(self.bcf.header, record)


########################################################################
########################################################################
## Variant File
########################################################################


cdef class VariantFile(HTSFile):
    """*(filename, mode=None, index_filename=None, header=None, drop_samples=False,
    duplicate_filehandle=True, ignore_truncation=False, threads=1)*

    A :term:`VCF`/:term:`BCF` formatted file. The file is automatically
    opened.

    If an index for a variant file exists (.csi or .tbi), it will be
    opened automatically.  Without an index random access to records
    via :meth:`fetch` is disabled.

    For writing, a :class:`VariantHeader` object must be provided,
    typically obtained from another :term:`VCF` file/:term:`BCF`
    file.

    Parameters
    ----------
    mode : string
        *mode* should be ``r`` for reading or ``w`` for writing. The default is
        text mode (:term:`VCF`).  For binary (:term:`BCF`) I/O you should append
        ``b`` for compressed or ``u`` for uncompressed :term:`BCF` output.

        If ``b`` is present, it must immediately follow ``r`` or ``w``.  Valid
        modes are ``r``, ``w``, ``wh``, ``rb``, ``wb``, ``wbu`` and ``wb0``.
        For instance, to open a :term:`BCF` formatted file for reading, type::

            f = pysam.VariantFile('ex1.bcf','r')

        If mode is not specified, we will try to auto-detect the file type.  All
        of the following should work::

            f1 = pysam.VariantFile('ex1.bcf')
            f2 = pysam.VariantFile('ex1.vcf')
            f3 = pysam.VariantFile('ex1.vcf.gz')

    index_filename : string
        Explicit path to an index file.

    header : VariantHeader
        :class:`VariantHeader` object required for writing.

    drop_samples: bool
        Ignore sample information when reading.

    duplicate_filehandle: bool
        By default, file handles passed either directly or through
        File-like objects will be duplicated before passing them to
        htslib. The duplication prevents issues where the same stream
        will be closed by htslib and through destruction of the
        high-level python object. Set to False to turn off
        duplication.

    ignore_truncation: bool
        Issue a warning, instead of raising an error if the current file
        appears to be truncated due to a missing EOF marker.  Only applies
        to bgzipped formats. (Default=False)

    threads: integer
        Number of threads to use for compressing/decompressing VCF/BCF files.
        Setting threads to > 1 cannot be combined with `ignore_truncation`.
        (Default=1)

    """
    def __cinit__(self, *args, **kwargs):
        self.htsfile = NULL

    def __init__(self, *args, **kwargs):
        self.header         = None
        self.index          = None
        self.filename       = None
        self.mode           = None
        self.threads        = 1
        self.index_filename = None
        self.is_stream      = False
        self.is_remote      = False
        self.is_reading     = False
        self.drop_samples   = False
        self.header_written = False
        self.start_offset   = -1

        self.open(*args, **kwargs)

    def __dealloc__(self):
        if not self.htsfile or not self.header:
            return

        # Write header if no records were written
        if self.htsfile.is_write and not self.header_written:
            with nogil:
                bcf_hdr_write(self.htsfile, self.header.ptr)

        cdef int ret = hts_close(self.htsfile)
        self.htsfile = NULL
        self.header = self.index = None

        if ret < 0:
            global errno
            if errno == EPIPE:
                errno = 0
            else:
                raise IOError(errno, force_str(strerror(errno)))

    def close(self):
        """closes the :class:`pysam.VariantFile`."""
        if not self.htsfile:
            return

        # Write header if no records were written
        if self.htsfile.is_write and not self.header_written:
            with nogil:
                bcf_hdr_write(self.htsfile, self.header.ptr)

        cdef int ret = hts_close(self.htsfile)
        self.htsfile = NULL
        self.header = self.index = None

        if ret < 0:
            global errno
            if errno == EPIPE:
                errno = 0
            else:
                raise IOError(errno, force_str(strerror(errno)))

    def __iter__(self):
        if not self.is_open:
            raise ValueError('I/O operation on closed file')

        if self.htsfile.is_write:
            raise ValueError('cannot iterate over Variantfile opened for writing')

        self.is_reading = 1
        return self

    def __next__(self):
        cdef int ret
        cdef int errcode
        cdef bcf1_t *record = bcf_init1()

        if not record:
            raise MemoryError('unable to allocate BCF record')

        record.pos = -1
        if self.drop_samples:
            record.max_unpack = BCF_UN_SHR

        with nogil:
            ret = bcf_read1(self.htsfile, self.header.ptr, record)

        if ret < 0:
            errcode = record.errcode
            bcf_destroy1(record)
            if errcode:
                raise IOError('unable to parse next record')
            if ret == -1:
                raise StopIteration
            elif ret == -2:
                raise IOError('truncated file')
            elif errno:
                raise IOError(errno, strerror(errno))
            else:
                raise IOError('unable to fetch next record')

        return makeVariantRecord(self.header, record)

    def copy(self):
        if not self.is_open:
            raise ValueError

        cdef VariantFile vars = VariantFile.__new__(VariantFile)
        cdef bcf_hdr_t *hdr

        # FIXME: re-open using fd or else header and index could be invalid
        vars.htsfile = self._open_htsfile()

        if not vars.htsfile:
            raise ValueError('Cannot re-open htsfile')

        # minimize overhead by re-using header and index.  This approach is
        # currently risky, but see above for how this can be mitigated.
        vars.header         = self.header
        vars.index          = self.index

        vars.filename       = self.filename
        vars.mode           = self.mode
        vars.threads        = self.threads
        vars.index_filename = self.index_filename
        vars.drop_samples   = self.drop_samples
        vars.is_stream      = self.is_stream
        vars.is_remote      = self.is_remote
        vars.is_reading     = self.is_reading
        vars.start_offset   = self.start_offset
        vars.header_written = self.header_written

        if self.htsfile.is_bin:
            vars.seek(self.tell())
        else:
            with nogil:
                hdr = bcf_hdr_read(vars.htsfile)
            makeVariantHeader(hdr)

        return vars

    def open(self, filename, mode='r',
             index_filename=None,
             VariantHeader header=None,
             drop_samples=False,
             duplicate_filehandle=True,
             ignore_truncation=False,
             threads=1):
        """open a vcf/bcf file.

        If open is called on an existing VariantFile, the current file will be
        closed and a new file will be opened.
        """
        cdef bcf_hdr_t *hdr
        cdef BGZF *bgzfp
        cdef hts_idx_t *idx
        cdef tbx_t *tidx
        cdef char *cfilename
        cdef char *cindex_filename = NULL
        cdef char *cmode

        if threads > 1 and ignore_truncation:
            # This won't raise errors if reaching a truncated alignment,
            # because bgzf_mt_reader in htslib does not deal with
            # bgzf_mt_read_block returning non-zero values, contrary
            # to bgzf_read (https://github.com/samtools/htslib/blob/1.7/bgzf.c#L888)
            # Better to avoid this (for now) than to produce seemingly correct results.
            raise ValueError('Cannot add extra threads when "ignore_truncation" is True')
        self.threads = threads

        # close a previously opened file
        if self.is_open:
            self.close()

        if not mode or mode[0] not in 'rwa':
            raise ValueError('mode must begin with r, w or a')

        self.duplicate_filehandle = duplicate_filehandle

        format_modes = [m for m in mode[1:] if m in 'bcguz']
        if len(format_modes) > 1:
            raise ValueError('mode contains conflicting format specifiers: {}'.format(''.join(format_modes)))

        invalid_modes = [m for m in mode[1:] if m not in 'bcguz0123456789ex']
        if invalid_modes:
            raise ValueError('invalid mode options: {}'.format(''.join(invalid_modes)))

        # Autodetect mode from filename
        if mode == 'w' and isinstance(filename, str):
            if filename.endswith('.gz'):
                mode = 'wz'
            elif filename.endswith('.bcf'):
                mode = 'wb'

        # for htslib, wbu seems to not work
        if mode == 'wbu':
            mode = 'wb0'

        self.mode = mode = force_bytes(mode)
        try:
            filename = encode_filename(filename)
            self.is_remote = hisremote(filename)
            self.is_stream = filename == b'-'
        except TypeError:
            filename = filename
            self.is_remote = False
            self.is_stream = True

        self.filename = filename

        if index_filename is not None:
            self.index_filename = index_filename = encode_filename(index_filename)
        else:
            self.index_filename = None

        self.drop_samples = bool(drop_samples)
        self.header = None

        self.header_written = False

        if mode.startswith(b'w'):
            # open file for writing
            if index_filename is not None:
                raise ValueError('Cannot specify an index filename when writing a VCF/BCF file')

            # header structure (used for writing)
            if header:
                self.header = header.copy()
            else:
                self.header = VariantHeader()
                #raise ValueError('a VariantHeader must be specified')

            # Header is not written until the first write or on close
            self.htsfile = self._open_htsfile()

            if not self.htsfile:
                raise ValueError("could not open file `{}` (mode='{}')".format(filename, mode))

        elif mode.startswith(b'r'):
            # open file for reading
            self.htsfile = self._open_htsfile()

            if not self.htsfile:
                if errno:
                    raise IOError(errno, 'could not open variant file `{}`: {}'.format(filename, force_str(strerror(errno))))
                else:
                    raise ValueError('could not open variant file `{}`'.format(filename))

            if self.htsfile.format.format not in (bcf, vcf):
                raise ValueError('invalid file `{}` (mode=`{}`) - is it VCF/BCF format?'.format(filename, mode))

            self.check_truncation(ignore_truncation)

            with nogil:
                hdr = bcf_hdr_read(self.htsfile)

            try:
                self.header = makeVariantHeader(hdr)
            except ValueError:
                raise ValueError('file `{}` does not have valid header (mode=`{}`) - is it VCF/BCF format?'.format(filename, mode))

            if isinstance(self.filename, bytes):
                cfilename = self.filename
            else:
                cfilename = NULL

            # check for index and open if present
            if self.htsfile.format.format == bcf and cfilename:
                if index_filename is not None:
                    cindex_filename = index_filename
                with nogil:
                    idx = bcf_index_load2(cfilename, cindex_filename)
                self.index = makeBCFIndex(self.header, idx)

            elif self.htsfile.format.compression == bgzf and cfilename:
                if index_filename is not None:
                    cindex_filename = index_filename
                with nogil:
                    tidx = tbx_index_load2(cfilename, cindex_filename)
                self.index = makeTabixIndex(tidx)

            if not self.is_stream:
                self.start_offset = self.tell()
        else:
            raise ValueError('unknown mode {}'.format(mode))

    def reset(self):
        """reset file position to beginning of file just after the header."""
        return self.seek(self.start_offset)

    def is_valid_tid(self, tid):
        """
        return True if the numerical :term:`tid` is valid; False otherwise.

        returns -1 if reference is not known.
        """
        if not self.is_open:
            raise ValueError('I/O operation on closed file')

        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef int rid = tid
        return 0 <= rid < hdr.n[BCF_DT_CTG]

    def get_tid(self, reference):
        """
        return the numerical :term:`tid` corresponding to
        :term:`reference`

        returns -1 if reference is not known.
        """
        if not self.is_open:
            raise ValueError('I/O operation on closed file')

        cdef vdict_t *d = <vdict_t*>self.header.ptr.dict[BCF_DT_CTG]
        reference = force_bytes(reference)
        cdef khint_t k = kh_get_vdict(d, reference)
        return kh_val_vdict(d, k).id if k != kh_end(d) else -1

    def get_reference_name(self, tid):
        """
        return :term:`reference` name corresponding to numerical :term:`tid`
        """
        if not self.is_open:
            raise ValueError('I/O operation on closed file')

        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef int rid = tid
        if rid < 0 or rid >= hdr.n[BCF_DT_CTG]:
            raise ValueError('Invalid tid')
        return bcf_str_cache_get_charptr(bcf_hdr_id2name(hdr, rid))

    def fetch(self, contig=None, start=None, stop=None, region=None, reopen=False, end=None, reference=None):
        """fetch records in a :term:`region`, specified either by
        :term:`contig`, *start*, and *end* (which are 0-based, half-open);
        or alternatively by a samtools :term:`region` string (which is
        1-based inclusive).

        Without *contig* or *region* all mapped records will be fetched.  The
        records will be returned ordered by contig, which will not necessarily
        be the order within the file.

        Set *reopen* to true if you will be using multiple iterators on the
        same file at the same time.  The iterator returned will receive its
        own copy of a filehandle to the file effectively re-opening the
        file.  Re-opening a file incurrs some overhead, so use with care.

        If only *contig* is set, all records on *contig* will be fetched.
        If both *region* and *contig* are given, an exception is raised.

        Note that a bgzipped :term:`VCF`.gz file without a tabix/CSI index
        (.tbi/.csi) or a :term:`BCF` file without a CSI index can only be
        read sequentially.
        """
        if not self.is_open:
            raise ValueError('I/O operation on closed file')

        if self.htsfile.is_write:
            raise ValueError('cannot fetch from Variantfile opened for writing')

        if contig is None and region is None:
            self.is_reading = 1
            bcf = self.copy() if reopen else self
            bcf.seek(self.start_offset)
            return iter(bcf)

        if self.index is None:
            raise ValueError('fetch requires an index')

        _, tid, start, stop = self.parse_region(contig, start, stop, region,
                                                None, end=end, reference=reference)

        if contig is None:
            contig = self.get_reference_name(tid)

        self.is_reading = 1
        return self.index.fetch(self, contig, start, stop, reopen)

    def new_record(self, *args, **kwargs):
        """Create a new empty :class:`VariantRecord`.

        See :meth:`VariantHeader.new_record`
        """
        return self.header.new_record(*args, **kwargs)

    cpdef int write(self, VariantRecord record) except -1:
        """
        write a single :class:`pysam.VariantRecord` to disk.

        returns the number of bytes written.
        """
        if record is None:
            raise ValueError('record must not be None')

        if not self.is_open:
            return ValueError('I/O operation on closed file')

        if not self.htsfile.is_write:
            raise ValueError('cannot write to a Variantfile opened for reading')

        if not self.header_written:
            self.header_written = True
            with nogil:
                bcf_hdr_write(self.htsfile, self.header.ptr)

        #if record.header is not self.header:
        #    record.translate(self.header)
        #    raise ValueError('Writing records from a different VariantFile is not yet supported')

        if record.ptr.n_sample != bcf_hdr_nsamples(self.header.ptr):
            msg = 'Invalid VariantRecord.  Number of samples does not match header ({} vs {})'
            raise ValueError(msg.format(record.ptr.n_sample, bcf_hdr_nsamples(self.header.ptr)))

        # Sync END annotation before writing
        bcf_sync_end(record)

        cdef int ret

        with nogil:
            ret = bcf_write1(self.htsfile, self.header.ptr, record.ptr)

        if ret < 0:
            raise IOError(errno, strerror(errno))

        return ret

    def subset_samples(self, include_samples):
        """
        Read only a subset of samples to reduce processing time and memory.
        Must be called prior to retrieving records.
        """
        if not self.is_open:
            raise ValueError('I/O operation on closed file')

        if self.htsfile.is_write:
            raise ValueError('cannot subset samples from Variantfile opened for writing')

        if self.is_reading:
            raise ValueError('cannot subset samples after fetching records')

        self.header._subset_samples(include_samples)

        # potentially unnecessary optimization that also sets max_unpack
        if not include_samples:
            self.drop_samples = True

