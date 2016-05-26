# cython: embedsignature=True
# cython: profile=True
###############################################################################
###############################################################################
## Cython wrapper for htslib VCF/BCF reader/writer
###############################################################################
#
# NOTICE: This code is incomplete and preliminary.  It offers a nearly
#         complete Pythonic interface to VCF/BCF metadata and data with
#         reading and writing capability.  It has limited capability to to
#         mutate the resulting data.  Documentation and a unit test suite
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
# Here is a quick tour through the API::
#
#     VariantFile(filename, mode=None, header=None, drop_samples=False)
#
#         Attributes / Properties
#
#             htsfile:      htsFile*                             [private]
#             start_offset: BGZF offset of first record          [private]
#             filename:     filename                             [read only]
#             mode:         mode                                 [read only]
#             header:       VariantHeader object                 [read only]
#             index:        TabixIndex, BCFIndex or None         [read only]
#             drop_samples: sample information is to be ignored  [read only]
#
#             is_stream:    file is stdin/stdout                 [read only]
#             is_remote:    file is not on the local filesystem  [read only]
#             is_reading:   file has begun reading records       [read only]
#             category:     file format general category         [read only]
#             format:       file format                          [read only]
#             version:      tuple of (major, minor) format version [read only]
#             compression:  file compression
#             description:  vaguely human readable description of  [read only]
#                           file format.
#
#         Methods:
#             copy()
#             close()
#             open(filename, mode=None, header=None, drop_samples=False)
#             reset()
#             seek(offset)
#             tell()
#             fetch(contig=None, start=None, stop=None, region=None, reopen=False)
#             subset_samples(include_samples)
#
#     VariantHeader()
#
#         version:      VCF version
#         samples:      sequence-like access to samples
#         records:      sequence-like access to partially parsed headers
#         contigs:      mapping-like object for contig name -> VariantContig
#
#         filters:      mapping-like object for filter name -> VariantMetadata
#         info:         mapping-like object for info name -> VariantMetadata
#         formats:      mapping-like object for formats name -> VariantMetadata
#
#     VariantRecord(...)
#
#         header:       VariantHeader object
#         rid:          reference id (i.e. tid)
#         chrom:        chromosome/contig string
#         contig:       synonym for chrom
#         pos:          1-based start position (inclusive)
#         start:        0-based start position (inclusive)
#         stop:         0-based stop position (exclusive)
#         rlen:         reference length (stop - start)
#         id:           record identifier
#         ref:          reference allele
#         alleles:      alleles (ref followed by alts)
#         alts:         alt alleles
#         qual:         quality (float)
#         filter:       mapping-like object for filter name -> type info
#         info:         mapping-like object for info name -> value
#         format:       mapping-like object for format name -> type info
#         samples:      mapping-like object of sample genotypes & attrs
#
#     VariantRecordSample(...)
#
#         name:           sample name
#         index:          sample index
#         allele_indices: tuple of allele indices (ref=0, alt=1..len(alts), missing=-1)
#         alleles:        tuple of alleles (missing=None)
#
#         VariantRecordSample is also a mapping object from formats to values
#
#     VariantContig(...)
#
#         id:           reference id (i.e. tid)
#         name:         chromosome/contig string
#         length:       contig length if provided, else None
#         header:       defining VariantHeaderRecord
#
#     VariantMetadata(...) # for FILTER, INFO and FORMAT metadata
#
#         id:           internal id
#         name:         metadata name
#         type:         value data type
#         number:       number of values
#         header:       defining VariantHeaderRecord
#
#    VariantHeaderRecord(...)  # replace with single tuple of key/value pairs?
#
#        type:          record type
#        key:           first record key
#        value:         first record value
#        attrs:         remaining key/value pairs
#
###############################################################################
#
# TODO list for next major sprint:
#
#   * more genotype methods
#   * unit test suite (perhaps py.test based)
#   * documentation
#   * htslib 1.2 format info
#
# For later sprints:
#
#   * ability to create indices
#   * mutable header and record data
#   * pickle support
#   * Python 3 support
#   * left/right locus normalization
#   * parallel iteration (like synced_bcf_reader)
#   * fix reopen to re-use fd
#
###############################################################################
#
# The MIT License
#
# Copyright (c) 2015 Kevin Jacobs (jacobs@bioinformed.com)
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

from libc.string cimport strcmp, strpbrk
from libc.stdint cimport INT8_MAX, INT16_MAX, INT32_MAX

cimport cython

from cpython.object  cimport PyObject
from cpython.ref     cimport Py_INCREF
from cpython.dict    cimport PyDict_GetItemString, PyDict_SetItemString
from cpython.tuple   cimport PyTuple_New, PyTuple_SET_ITEM
from cpython.bytes   cimport PyBytes_FromStringAndSize
from cpython.unicode cimport PyUnicode_DecodeASCII
from cpython.version cimport PY_MAJOR_VERSION

from pysam.chtslib   cimport hisremote


from warnings         import warn


__all__ = ['VariantFile',
           'VariantHeader',
           'VariantHeaderRecord',
           'VariantRecord']

########################################################################
########################################################################
## Constants
########################################################################

cdef int   MAX_POS = 2 << 29
cdef tuple VALUE_TYPES = ('Flag', 'Integer', 'Float', 'String')
cdef tuple METADATA_TYPES = ('FILTER', 'INFO', 'FORMAT', 'CONTIG', 'STRUCTURED', 'GENERIC')
cdef tuple METADATA_LENGTHS = ('FIXED', 'VARIABLE', 'A', 'G', 'R')

cdef tuple FORMAT_CATEGORIES = ('UNKNOWN', 'ALIGNMENTS', 'VARIANTS', 'INDEX', 'REGIONS')
cdef tuple FORMATS = ('UNKNOWN', 'BINARY_FORMAT', 'TEXT_FORMAT', 'SAM', 'BAM', 'BAI', 'CRAM', 'CRAI',
                      'VCF', 'BCF', 'CSI', 'GZI', 'TBI', 'BED')
cdef tuple COMPRESSION = ('NONE', 'GZIP', 'BGZF', 'CUSTOM')

########################################################################
########################################################################
## Python 3 compatibility functions
########################################################################

from pysam.cutils cimport force_bytes, force_str, charptr_to_str, charptr_to_str_w_len
from pysam.cutils cimport encode_filename, from_string_and_size


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

    if PY_MAJOR_VERSION < 3:
        val = s
    else:
        val = PyUnicode_DecodeASCII(s, strlen(s), NULL)

    PyDict_SetItemString(bcf_str_cache, s, val)

    return val


########################################################################
########################################################################
## Low level type conversion helpers
########################################################################


cdef inline int is_gt_fmt(bcf_hdr_t *hdr, int fmt_id):
    return strcmp(bcf_hdr_int2id(hdr, BCF_DT_ID, fmt_id), "GT") == 0


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

    if not data or n <= 0:
        return None

    if type == BCF_BT_CHAR:
        datac = <char *>data
        while n and datac[n-1] == bcf_str_vector_end:
            n -= 1
        value = charptr_to_str_w_len(datac, n) if datac[0] != bcf_str_missing else None
        # FIXME: Need to know length?  Report errors?  Pad with missing values?  Not clear what to do.

        value = tuple(v or None for v in value.split(',')) if value else ()
        # FIXME: Need to know length?  Report errors?  Pad with missing values?  Not clear what to do.
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

    assert(value_count <= n)

    if bt_type == BCF_BT_CHAR:
        if not isinstance(values, (str, bytes)):
            values = b','.join(force_bytes(v) if v is not None else b'' for v in values)
            value_count = len(values)
        assert(value_count <= n)
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


cdef bcf_copy_expand_array(void *src_data, int src_type, ssize_t src_values,
                           void *dst_data, int dst_type, ssize_t dst_values,
                           int vlen):
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
        memcpy(src_datac, dst_datac, src_values)
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


cdef bcf_get_value_count(VariantRecord record, int hl_type, int id, ssize_t *count, int *scalar):
    cdef bcf_hdr_t *hdr = record.header.ptr
    cdef bcf1_t *r = record.ptr
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
        count[0] = r.n_allele * (r.n_allele + 1) // 2
    elif length == BCF_VL_VAR:
        count[0] = -1
    else:
        raise ValueError('Unknown format length')


cdef object bcf_info_get_value(VariantRecord record, const bcf_info_t *z):
    cdef bcf_hdr_t *hdr = record.header.ptr

    cdef char *s
    cdef ssize_t count
    cdef int scalar

    bcf_get_value_count(record, BCF_HL_INFO, z.key, &count, &scalar)

    if z.len == 0:
        if  bcf_hdr_id2type(hdr, BCF_HL_INFO, z.key) == BCF_HT_FLAG:
            value = True
        elif scalar:
            value = None
        else:
            value = ()
    elif z.len == 1:
        if z.type == BCF_BT_INT8:
            if z.v1.i == bcf_int8_missing:
                value = None
            elif z.v1.i == bcf_int8_vector_end:
                value = ()
            else:
                value = z.v1.i
        elif z.type == BCF_BT_INT16:
            if z.v1.i == bcf_int16_missing:
                value = None
            elif z.v1.i == bcf_int16_vector_end:
                value = ()
            else:
                value = z.v1.i
        elif z.type == BCF_BT_INT32:
            if z.v1.i == bcf_int32_missing:
                value = None
            elif z.v1.i == bcf_int32_vector_end:
                value = ()
            else:
                value = z.v1.i
        elif z.type == BCF_BT_FLOAT:
            if bcf_float_is_missing(z.v1.f):
                value = None
            elif bcf_float_is_vector_end(z.v1.f):
                value = ()
            else:
                value = z.v1.f
        elif z.type == BCF_BT_CHAR:
            value = force_str(chr(z.v1.i))
        else:
            raise TypeError('unsupported info type code')

        if not scalar and value != ():
            value = (value,)
    else:
        value = bcf_array_to_object(z.vptr, z.type, z.len, count, scalar)

    return value


cdef object bcf_check_values(VariantRecord record, value, int hl_type, int ht_type,
                             int id, int bt_type, ssize_t bt_len, ssize_t *value_count,
                             int *scalar, int *realloc):

    bcf_get_value_count(record, hl_type, id, value_count, scalar)

    # Validate values now that we know the type and size
    values = (value,) if not isinstance(value, tuple) else value

    # Validate values now that we know the type and size
    if ht_type == BCF_HT_FLAG:
        value_count[0] = 1

    if value_count[0] != -1 and value_count[0] != len(values):
        if scalar[0]:
            raise TypeError('value expected to be scalar'.format(value_count[0]))
        else:
            raise TypeError('values expected to be {:d}-tuple'.format(value_count[0]))

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
    cdef bcf1_t *r = record.ptr
    cdef int32_t nalleles = r.n_allele
    cdef list gt_values = []
    cdef char *s
    cdef int i

    if not values:
        return ()

    if not isinstance(values, (list, tuple)):
        values = (values,)

    for value in values:
        if value is None:
            gt_values.append(None)
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
    cdef bcf_hdr_t *hdr = record.header.ptr
    cdef bcf1_t *r = record.ptr
    cdef vdict_t *d
    cdef khiter_t k
    cdef int info_id, info_type, scalar, dst_type, realloc, vlen = 0
    cdef ssize_t i, value_count, alloc_len, alloc_size, dst_size

    if bcf_unpack(r, BCF_UN_INFO) < 0:
        raise ValueError('Error unpacking VariantRecord')

    bkey = force_bytes(key)
    cdef bcf_info_t *info = bcf_get_info(hdr, r, bkey)

    if info:
        info_id = info.key
    else:
        d = <vdict_t *>hdr.dict[BCF_DT_ID]
        k = kh_get_vdict(d, bkey)

        if k == kh_end(d) or kh_val_vdict(d, k).info[BCF_HL_INFO] & 0xF == 0xF:
            raise KeyError('unknown INFO')

        info_id = kh_val_vdict(d, k).id

    info_type = bcf_hdr_id2type(hdr, BCF_HL_INFO, info_id)
    values = bcf_check_values(record, value, BCF_HL_INFO, info_type, info_id,
                              info.type if info else -1, info.len if info else -1,
                              &value_count, &scalar, &realloc)

    if info_type == BCF_HT_FLAG:
        if bcf_update_info(hdr, r, bkey, NULL, bool(values[0]), info_type) < 0:
            raise ValueError('Unable to update INFO values')
        return

    vlen = value_count < 0
    value_count = len(values)

    # If we can, write updated values to existing allocated storage
    if info and not realloc:
        r.d.shared_dirty |= BCF1_DIRTY_INF

        if value_count == 0:
            info.len = 0
            # FIXME: Check if need to free vptr if info.len > 0?
        elif value_count == 1:
            # FIXME: Check if need to free vptr if info.len > 0?
            if info.type == BCF_BT_INT8 or info.type == BCF_BT_INT16 or info.type == BCF_BT_INT32:
                bcf_object_to_array(values, &info.v1.i, BCF_BT_INT32, 1, vlen)
            elif info.type == BCF_BT_FLOAT:
                bcf_object_to_array(values, &info.v1.f, BCF_BT_FLOAT, 1, vlen)
            else:
                raise TypeError('unsupported info type code')
            info.len = 1
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
    cdef bcf_hdr_t *hdr = record.header.ptr
    cdef bcf1_t *r = record.ptr
    cdef ssize_t value_count
    cdef int scalar

    if bcf_unpack(r, BCF_UN_INFO) < 0:
        raise ValueError('Error unpacking VariantRecord')

    bkey = force_bytes(key)
    cdef bcf_info_t *info = bcf_get_info(hdr, r, bkey)

    if not info:
        raise KeyError(key)

    bcf_get_value_count(record, BCF_HL_INFO, info.key, &value_count, &scalar)

    if value_count <= 0:
        null_value = ()
    elif scalar:
        null_value = None
    else:
        null_value = (None,)*value_count

    bcf_info_set_value(record, bkey, null_value)


cdef bcf_format_get_value(VariantRecordSample sample, key):
    cdef bcf_hdr_t *hdr = sample.record.header.ptr
    cdef bcf1_t *r = sample.record.ptr
    cdef ssize_t count
    cdef int scalar

    if bcf_unpack(r, BCF_UN_ALL) < 0:
        raise ValueError('Error unpacking VariantRecord')

    bkey = force_bytes(key)
    cdef bcf_fmt_t *fmt = bcf_get_fmt(hdr, r, bkey)

    if not fmt or not fmt.p:
        raise KeyError('invalid FORMAT')

    if is_gt_fmt(hdr, fmt.id):
        return bcf_format_get_allele_indices(sample)

    bcf_get_value_count(sample.record, BCF_HL_FMT, fmt.id, &count, &scalar)

    if fmt.p and fmt.n and fmt.size:
        return bcf_array_to_object(fmt.p + sample.index * fmt.size, fmt.type, fmt.n, count, scalar)
    elif scalar:
        return None
    elif count <= 0:
        return ()
    else:
        return (None,)*count


cdef bcf_format_set_value(VariantRecordSample sample, key, value):
    cdef bcf_hdr_t *hdr = sample.record.header.ptr
    cdef bcf1_t *r = sample.record.ptr
    cdef int fmt_id
    cdef vdict_t *d
    cdef khiter_t k
    cdef int fmt_type, scalar, realloc, dst_type, vlen = 0
    cdef ssize_t i, n, value_count, alloc_size, alloc_len, dst_size

    if bcf_unpack(r, BCF_UN_ALL) < 0:
        raise ValueError('Error unpacking VariantRecord')

    bkey = force_bytes(key)
    cdef bcf_fmt_t *fmt = bcf_get_fmt(hdr, r, bkey)

    if fmt:
        fmt_id = fmt.id
    else:
        d = <vdict_t *>hdr.dict[BCF_DT_ID]
        k = kh_get_vdict(d, bkey)

        if k == kh_end(d) or kh_val_vdict(d, k).info[BCF_HL_FMT] & 0xF == 0xF:
            raise KeyError('unknown format')

        fmt_id = kh_val_vdict(d, k).id

    fmt_type = bcf_hdr_id2type(hdr, BCF_HL_FMT, fmt_id)

    if fmt_type == BCF_HT_FLAG:
        raise ValueError('Flag types are not allowed on FORMATs')

    if is_gt_fmt(hdr, fmt_id):
        value = bcf_encode_alleles(sample.record, value)

    values = bcf_check_values(sample.record, value, BCF_HL_FMT, fmt_type, fmt_id,
                              fmt.type if fmt else -1, fmt.n if fmt else -1,
                              &value_count, &scalar, &realloc)

    vlen = value_count < 0
    value_count = len(values)

    # If we can, write updated values to existing allocated storage
    if fmt and not realloc:
        r.d.indiv_dirty = 1
        bcf_object_to_array(values, fmt.p + sample.index * fmt.size, fmt.type, fmt.n, vlen)
        return

    alloc_len = max(1, value_count)
    if fmt and fmt.n > alloc_len:
        alloc_len = fmt.n

    n = bcf_hdr_nsamples(hdr)
    new_values = bcf_empty_array(fmt_type, n*alloc_len, vlen)
    cdef char *valp = <char *>new_values

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

    if fmt and n > 1:
        for i in range(n):
            bcf_copy_expand_array(fmt.p + i*fmt.size, fmt.type, fmt.n,
                                  valp  + i*dst_size, dst_type, alloc_len,
                                  vlen)

    bcf_object_to_array(values, valp + sample.index*dst_size, dst_type, alloc_len, vlen)

    if bcf_update_format(hdr, r, bkey, valp, <int>(n*alloc_len), fmt_type) < 0:
        raise ValueError('Unable to update format values')


cdef bcf_format_del_value(VariantRecordSample sample, key):
    cdef bcf_hdr_t *hdr = sample.record.header.ptr
    cdef bcf1_t *r = sample.record.ptr
    cdef ssize_t value_count
    cdef int scalar

    if bcf_unpack(r, BCF_UN_ALL) < 0:
        raise ValueError('Error unpacking VariantRecord')

    bkey = force_bytes(key)
    cdef bcf_fmt_t *fmt = bcf_get_fmt(hdr, r, bkey)

    if not fmt or not fmt.p:
        raise KeyError(key)

    bcf_get_value_count(sample.record, BCF_HL_FMT, fmt.id, &value_count, &scalar)

    if value_count <= 0:
        null_value = ()
    elif scalar:
        null_value = None
    else:
        null_value = (None,)*value_count

    bcf_format_set_value(sample, bkey, null_value)


cdef bcf_format_get_allele_indices(VariantRecordSample sample):
    cdef bcf_hdr_t *hdr = sample.record.header.ptr
    cdef bcf1_t *r = sample.record.ptr
    cdef int32_t n = bcf_hdr_nsamples(hdr)

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
            elif data8[i] == bcf_int8_missing:
                a = -1
            else:
                a = bcf_gt_allele(data8[i])
            alleles.append(a if 0 <= a < nalleles else None)
    elif fmt0.type == BCF_BT_INT16:
        data16 = <int16_t *>(fmt0.p + sample.index * fmt0.size)
        for i in range(fmt0.n):
            if data16[i] == bcf_int16_vector_end:
                break
            elif data16[i] == bcf_int16_missing:
                a = -1
            else:
                a = bcf_gt_allele(data16[i])
            alleles.append(a if 0 <= a < nalleles else None)
    elif fmt0.type == BCF_BT_INT32:
        data32 = <int32_t *>(fmt0.p + sample.index * fmt0.size)
        for i in range(fmt0.n):
            if data32[i] == bcf_int32_vector_end:
                break
            elif data32[i] == bcf_int32_missing:
                a = -1
            else:
                a = bcf_gt_allele(data32[i])
            alleles.append(a if 0 <= a < nalleles else None)

    return tuple(alleles)


cdef bcf_format_get_alleles(VariantRecordSample sample):
    cdef bcf_hdr_t *hdr = sample.record.header.ptr
    cdef bcf1_t *r = sample.record.ptr
    cdef int32_t nsamples = bcf_hdr_nsamples(hdr)

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
    cdef bcf_hdr_t *hdr = sample.record.header.ptr
    cdef bcf1_t *r = sample.record.ptr
    cdef int32_t n = bcf_hdr_nsamples(hdr)

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
    cdef bcf_hdr_t *hdr = sample.record.header.ptr
    cdef bcf1_t *r = sample.record.ptr
    cdef int32_t n = bcf_hdr_nsamples(hdr)

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


########################################################################
########################################################################
## Variant Header objects
########################################################################

#FIXME: implement a full mapping interface
#FIXME: passing bcf_hrec_t*  may not be the safest approach once mutating
#       operations are allowed.
cdef class VariantHeaderRecord(object):
    """header record from a :class:`VariantHeader` object"""

    property type:
        """header type: FILTER, INFO, FORMAT, CONTIG, STRUCTURED, or GENERIC"""
        def __get__(self):
            cdef bcf_hrec_t *r = self.ptr
            return METADATA_TYPES[r.type]

    property key:
        """header key (the part before '=', in FILTER/INFO/FORMAT/contig/fileformat etc.)"""
        def __get__(self):
            cdef bcf_hrec_t *r = self.ptr
            return bcf_str_cache_get_charptr(r.key) if r.key else None

    property value:
        """header value.  Set only for generic lines, None for FILTER/INFO, etc."""
        def __get__(self):
            cdef bcf_hrec_t *r = self.ptr
            return charptr_to_str(r.value) if r.value else None

    property attrs:
        """sequence of additional header attributes"""
        def __get__(self):
            cdef bcf_hrec_t *r = self.ptr
            cdef int i
            return tuple((bcf_str_cache_get_charptr(r.keys[i]) if r.keys[i] else None,
                          charptr_to_str(r.vals[i]) if r.vals[i] else None)
                         for i in range(r.nkeys))

    def __len__(self):
        cdef bcf_hrec_t *r = self.ptr
        return r.nkeys

    def __bool__(self):
        cdef bcf_hrec_t *r = self.ptr
        return r.nkeys != 0

    def __getitem__(self, key):
        """get attribute value"""
        cdef bcf_hrec_t *r = self.ptr
        cdef int i
        bkey = force_bytes(key)
        for i in range(r.nkeys):
            if r.keys[i] and r.keys[i] == bkey:
                return charptr_to_str(r.vals[i]) if r.vals[i] else None
        raise KeyError('cannot find metadata key')

    def __iter__(self):
        cdef bcf_hrec_t *r = self.ptr
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
        cdef int i
        for i in range(r.nkeys):
            if r.keys[i]:
                yield charptr_to_str(r.vals[i]) if r.vals[i] else None

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        cdef bcf_hrec_t *r = self.ptr
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

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None

    #TODO: implement __richcmp__

    def __str__(self):
        cdef bcf_hrec_t *r = self.ptr
        if r.type == BCF_HL_GEN:
            return '##{}={}'.format(self.key, self.value)
        else:
            attrs = ','.join('{}={}'.format(k, v) for k,v in self.attrs if k != 'IDX')
            return '##{}=<{}>'.format(self.key or self.type, attrs)


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
    """filter, info or format metadata record from a :class:`VariantHeader`
    object"""

    property name:
        """metadata name"""
        def __get__(self):
            cdef bcf_hdr_t *hdr = self.header.ptr
            return bcf_str_cache_get_charptr(hdr.id[BCF_DT_ID][self.id].key)

    # Q: Should this be exposed?
    property id:
        """metadata internal header id number"""
        def __get__(self):
            return self.id

    property number:
        """metadata number (i.e. cardinality)"""
        def __get__(self):
            cdef bcf_hdr_t *hdr = self.header.ptr
            if not bcf_hdr_idinfo_exists(hdr, self.type, self.id) or self.type == BCF_HL_FLT:
                return None
            cdef int l = bcf_hdr_id2length(hdr, self.type, self.id)
            if l == BCF_VL_FIXED:
                return bcf_hdr_id2number(hdr, self.type, self.id)
            elif l == BCF_VL_VAR:
                return '.'
            else:
                return METADATA_LENGTHS[l]

    property type:
        """metadata value type"""
        def __get__(self):
            cdef bcf_hdr_t *hdr = self.header.ptr
            if not bcf_hdr_idinfo_exists(hdr, self.type, self.id) or \
               self.type == BCF_HL_FLT:
                return None
            return VALUE_TYPES[bcf_hdr_id2type(hdr, self.type, self.id)]

    property description:
        """metadata description (or None if not set)"""
        def __get__(self):
            descr = self.record.get('Description')
            if descr:
                descr = descr.strip('"')
            return force_str(descr)

    property record:
        """:class:`VariantHeaderRecord` associated with this
        :class:`VariantMetadata` object"""
        def __get__(self):
            cdef bcf_hdr_t *hdr = self.header.ptr
            if not bcf_hdr_idinfo_exists(hdr, self.type, self.id):
                return None
            cdef bcf_hrec_t *hrec = hdr.id[BCF_DT_ID][self.id].val.hrec[self.type]
            if not hrec:
                return None
            return makeVariantHeaderRecord(self.header, hrec)


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

    def add(self, id, number, type, description, **kwargs):
        """Add a new filter, info or format record"""
        if id in self:
            raise ValueError('Header already exists for id={}'.format(id))

        if self.type == BCF_HL_FLT:
            if number is not None:
                raise ValueError('Number must be None when adding a filter')
            if type is not None:
                raise ValueError('Type must be None when adding a filter')

            items = [('ID', id), ('Description', description)]
        else:
            if type not in VALUE_TYPES:
                raise ValueError('unknown type specified: {}'.format(type))
            if number is None:
                number = '.'

            items = [('ID', id),
                     ('Number', number),
                     ('Type', type),
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

        bkey = force_bytes(key)
        cdef khiter_t k = kh_get_vdict(d, bkey)

        if k == kh_end(d) or kh_val_vdict(d, k).info[self.type] & 0xF == 0xF:
            raise KeyError('invalid filter')

        return makeVariantMetadata(self.header, self.type, kh_val_vdict(d, k).id)

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

    property name:
        """contig name"""
        def __get__(self):
            cdef bcf_hdr_t *hdr = self.header.ptr
            return bcf_str_cache_get_charptr(hdr.id[BCF_DT_CTG][self.id].key)

    property id:
        """contig internal id number"""
        def __get__(self):
            return self.id

    property length:
        """contig length or None if not available"""
        def __get__(self):
            cdef bcf_hdr_t *hdr = self.header.ptr
            cdef uint32_t length = hdr.id[BCF_DT_CTG][self.id].val.info[0]
            return length if length else None

    property header:
        """:class:`VariantHeaderRecord` associated with this :class:`VariantContig` object"""
        def __get__(self):
            cdef bcf_hdr_t *hdr = self.header.ptr
            cdef bcf_hrec_t *hrec = hdr.id[BCF_DT_CTG][self.id].val.hrec[0]
            return makeVariantHeaderRecord(self.header, hrec)


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
        bkey = force_bytes(key)
        cdef khiter_t k = kh_get_vdict(d, bkey)

        if k == kh_end(d):
            raise KeyError('invalid contig')

        cdef int id = kh_val_vdict(d, k).id

        return makeVariantContig(self.header, id)

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

    def add(self, id, **kwargs):
        """Add a new contig record"""
        if id in self:
            raise ValueError('Header already exists for contig {}'.format(id))

        items = [('ID', id)] + kwargs.items()
        self.header.add_meta('contig', items=items)


cdef VariantHeaderContigs makeVariantHeaderContigs(VariantHeader header):
    if not header:
        raise ValueError('invalid VariantHeader')

    cdef VariantHeaderContigs contigs = VariantHeaderContigs.__new__(VariantHeaderContigs)
    contigs.header = header

    return contigs


cdef class VariantHeaderSamples(object):
    """sequence of sample names from a :class:`VariantHeader` object"""

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
        bkey = force_bytes(key)
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
        # self.ptr == NULL should be impossible
        return self.ptr != NULL

    def copy(self):
        return makeVariantHeader(bcf_hdr_dup(self.ptr))

    property version:
        """VCF version"""
        def __get__(self):
            return force_str(bcf_hdr_get_version(self.ptr))

    property samples:
        """samples (:class:`VariantHeaderSamples`)"""
        def __get__(self):
            return makeVariantHeaderSamples(self)

    property records:
        """header records (:class:`VariantHeaderRecords`)"""
        def __get__(self):
            return makeVariantHeaderRecords(self)

    property contigs:
        """contig information (:class:`VariantHeaderContigs`)"""
        def __get__(self):
            return makeVariantHeaderContigs(self)

    property filters:
        """filter metadata (:class:`VariantHeaderMetadata`)"""
        def __get__(self):
            return makeVariantHeaderMetadata(self, BCF_HL_FLT)

    property info:
        """info metadata (:class:`VariantHeaderMetadata`)"""
        def __get__(self):
            return makeVariantHeaderMetadata(self, BCF_HL_INFO)

    property formats:
        """format metadata (:class:`VariantHeaderMetadata`)"""
        def __get__(self):
            return makeVariantHeaderMetadata(self, BCF_HL_FMT)

    property alts:
        """alt metadata (:class:`dict` ID->record).

        The data returned just a snapshot of alt records, is created
        every time the property is requested, and modifications will
        not be reflected in the header metadata and vice versa.

        i.e. it is just a dict that reflects the state of alt records
        at the time it is created.

        """
        def __get__(self):
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
        cdef char *hstr = bcf_hdr_fmt_text(self.ptr, 0, &hlen)

        try:
            return charptr_to_str_w_len(hstr, hlen)
        finally:
            free(hstr)

    def add_record(self, VariantHeaderRecord record):
        """Add an existing :class:`VariantHeaderRecord` to this header"""
        cdef bcf_hrec_t *r = record.ptr

        if r.type == BCF_HL_GEN:
            self.add_meta(r.key, r.value)
        else:
            items = [(k,v) for k,v in record.attrs if k != 'IDX']
            self.add_meta(r.key, items=items)

    def add_line(self, line):
        """Add a metadata line to this header"""
        bline = force_bytes(line)
        if bcf_hdr_append(self.ptr, bline) < 0:
            raise ValueError('invalid header line')

        if self.ptr.dirty:
            bcf_hdr_sync(self.ptr)

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
                    key = force_bytes(key)
                    bcf_hrec_add_key(hrec, key, <int>len(key))

                    value = force_bytes(str(value))
                    quoted = strpbrk(value, ' ;,"\t<>') != NULL
                    bcf_hrec_set_val(hrec, hrec.nkeys-1, value, <int>len(value), quoted)
        except:
            bcf_hrec_destroy(hrec)
            raise

        bcf_hdr_add_hrec(self.ptr, hrec)

        if self.ptr.dirty:
            bcf_hdr_sync(self.ptr)

    def add_sample(self, name):
        """Add a new sample to this header"""
        bname = force_bytes(name)
        if bcf_hdr_add_sample(self.ptr, bname) < 0:
            raise ValueError('Duplicated sample name: {}'.format(name))
        if self.ptr.dirty:
            bcf_hdr_sync(self.ptr)


cdef VariantHeader makeVariantHeader(bcf_hdr_t *hdr):
    if not hdr:
        raise ValueError('cannot create VariantHeader')

    cdef VariantHeader header = VariantHeader.__new__(VariantHeader)
    header.ptr = hdr

    return header


########################################################################
########################################################################
## Variant Record objects
########################################################################

cdef class VariantRecordFilter(object):
    """Filters set on a :class:`VariantRecord` object, presented as a mapping from
       filter index or name to :class:`VariantMetadata` object"""

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

            if not bcf_hdr_idinfo_exists(hdr, BCF_HL_FLT, id) \
               or not bcf_has_filter(hdr, self.record.ptr, bkey):
                raise KeyError('Invalid filter')

        return makeVariantMetadata(self.record.header, BCF_HL_FLT, id)

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

            if not bcf_hdr_idinfo_exists(hdr, BCF_HL_FLT, id) \
               or not bcf_has_filter(hdr, self.record.ptr, bkey):
                raise KeyError('Invalid filter')

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
        bkey = force_bytes(key)
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

        bkey = force_bytes(key)
        cdef bcf_fmt_t *fmt = bcf_get_fmt(hdr, r, bkey)

        if not fmt or not fmt.p:
            raise KeyError('unknown format')

        return makeVariantMetadata(self.record.header, BCF_HL_FMT, fmt.id)

    def __delitem__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr

        bkey = force_bytes(key)
        cdef bcf_fmt_t *fmt = bcf_get_fmt(hdr, r, bkey)

        if not fmt or not fmt.p:
            raise KeyError('unknown format')

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
        bkey = force_bytes(key)
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

    cdef VariantRecordFormat format = VariantRecordFormat.__new__(
        VariantRecordFormat)
    format.record = record

    return format


#TODO: Add a getmeta method to return the corresponding VariantMetadata?
cdef class VariantRecordInfo(object):
    """Info data stored in a :class:`VariantRecord` object, presented as a
       mapping from info metadata name to value."""

    def __len__(self):
        return self.record.ptr.n_info

    def __bool__(self):
        return self.record.ptr.n_info != 0

    def __getitem__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef vdict_t *d
        cdef khiter_t k
        cdef info_id

        if bcf_unpack(r, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        bkey = force_bytes(key)
        cdef bcf_info_t *info = bcf_get_info(hdr, r, bkey)

        if not info:
            d = <vdict_t *>hdr.dict[BCF_DT_ID]
            k = kh_get_vdict(d, bkey)

            if k == kh_end(d) or kh_val_vdict(d, k).info[BCF_HL_INFO] & 0xF == 0xF:
                raise KeyError('Unknown INFO field: {}'.format(key))

            info_id = kh_val_vdict(d, k).id
        else:
            info_id = info.key

        if bcf_hdr_id2type(hdr, BCF_HL_INFO, info_id) == BCF_HT_FLAG:
            return info != NULL and info.vptr != NULL

        if not info or not info.vptr:
            raise KeyError('Invalid INFO field: {}'.format(key))

        return bcf_info_get_value(self.record, info)

    def __setitem__(self, key, value):
        bcf_info_set_value(self.record, key, value)

    def __delitem__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr

        if bcf_unpack(r, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        bkey = force_bytes(key)
        cdef bcf_info_t *info = bcf_get_info(hdr, r, bkey)

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
                if bcf_update_info(hdr, r, key, NULL, 0, info.type) < 0:
                    raise ValueError('Unable to delete INFO')

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info
        cdef int i

        for i in range(r.n_info):
            info = &r.d.info[i]
            if info and info.vptr:
                yield bcf_str_cache_get_charptr(bcf_hdr_int2id(hdr, BCF_DT_ID, info.key))

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr

        if bcf_unpack(r, BCF_UN_INFO) < 0:
            raise ValueError('Error unpacking VariantRecord')

        bkey = force_bytes(key)
        cdef bcf_info_t *info = bcf_get_info(hdr, r, bkey)

        return info != NULL

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

    def itervalues(self):
        """D.itervalues() -> an iterator over the values of D"""
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info
        cdef int i

        for i in range(r.n_info):
            info = &r.d.info[i]
            if info and info.vptr:
                yield bcf_info_get_value(self.record, info)

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info
        cdef int i

        for i in range(r.n_info):
            info = &r.d.info[i]
            if info and info.vptr:
                key = bcf_hdr_int2id(hdr, BCF_DT_ID, info.key)
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

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None

    #TODO: implement __richcmp__


cdef VariantRecordInfo makeVariantRecordInfo(VariantRecord record):
    if not record:
        raise ValueError('invalid VariantRecord')

    cdef VariantRecordInfo info = VariantRecordInfo.__new__(VariantRecordInfo)
    info.record = record

    return info


cdef class VariantRecordSamples(object):
    """mapping from sample index or name to :class:`VariantRecordSample` object."""

    def __len__(self):
        return bcf_hdr_nsamples(self.record.header.ptr)

    def __bool__(self):
        return bcf_hdr_nsamples(self.record.header.ptr) != 0

    def __getitem__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int n = bcf_hdr_nsamples(hdr)
        cdef int sample_index
        cdef vdict_t *d
        cdef khiter_t k

        if isinstance(key, int):
            sample_index = key
        else:
            bkey = force_bytes(key)
            sample_index = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, bkey)
            if sample_index < 0:
                raise KeyError('invalid sample name')

        if sample_index < 0 or sample_index >= n:
            raise IndexError('invalid sample index')

        return makeVariantRecordSample(self.record, sample_index)

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int32_t i, n = bcf_hdr_nsamples(hdr)

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
        cdef int n = bcf_hdr_nsamples(hdr)
        cdef int sample_index
        cdef vdict_t *d
        cdef khiter_t k

        if isinstance(key, int):
            sample_index = key
        else:
            bkey = force_bytes(key)
            sample_index = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, bkey)
            if sample_index < 0:
                raise KeyError('invalid sample name')

        return 0 <= sample_index < n

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

    def itervalues(self):
        """D.itervalues() -> an iterator over the values of D"""
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int32_t i, n = bcf_hdr_nsamples(hdr)

        for i in range(n):
            yield makeVariantRecordSample(self.record, i)

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int32_t i, n = bcf_hdr_nsamples(hdr)

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

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None

    #TODO: implement __richcmp__


cdef VariantRecordSamples makeVariantRecordSamples(VariantRecord record):
    if not record:
        raise ValueError('invalid VariantRecord')

    cdef VariantRecordSamples samples = VariantRecordSamples.__new__(
        VariantRecordSamples)
    samples.record = record

    return samples


cdef class VariantRecord(object):
    """Variant record"""

    def __dealloc__(self):
        if self.ptr:
            bcf_destroy1(self.ptr)
            self.ptr = NULL

    property rid:
        """internal reference id number"""
        def __get__(self):
            return self.ptr.rid
        def __set__(self, rid):
            cdef bcf_hdr_t *hdr = self.header.ptr
            cdef int r = rid
            if rid < 0 or r >= hdr.n[BCF_DT_CTG] or not hdr.id[BCF_DT_CTG][r].val:
                raise ValueError('invalid reference id')
            self.ptr.rid = r

    property chrom:
        """chromosome/contig name"""
        def __get__(self):
            return bcf_str_cache_get_charptr(bcf_hdr_id2name(self.header.ptr, self.ptr.rid))
        def __set__(self, chrom):
            cdef vdict_t *d = <vdict_t*>self.header.ptr.dict[BCF_DT_CTG]
            bchrom = force_bytes(chrom)
            cdef khint_t k = kh_get_vdict(d, bchrom)
            if k == kh_end(d):
                raise ValueError('Invalid chromosome/contig')
            self.ptr.rid = kh_val_vdict(d, k).id

    property contig:
        """chromosome/contig name"""
        def __get__(self):
            return bcf_str_cache_get_charptr(bcf_hdr_id2name(self.header.ptr, self.ptr.rid))
        def __set__(self, chrom):
            cdef vdict_t *d = <vdict_t*>self.header.ptr.dict[BCF_DT_CTG]
            bchrom = force_bytes(chrom)
            cdef khint_t k = kh_get_vdict(d, bchrom)
            if k == kh_end(d):
                raise ValueError('Invalid chromosome/contig')
            self.ptr.rid = kh_val_vdict(d, k).id

    property pos:
        """record start position on chrom/contig (1-based inclusive)"""
        def __get__(self):
            return self.ptr.pos + 1
        def __set__(self, pos):
            if pos < 1:
                raise ValueError('Position must be positive')
            # FIXME: check start <= stop?
            #   KBJ: Can't or else certain mutating operations will become
            #        difficult or impossible.  e.g.  having to delete
            #        info['END'] before being able to reset pos is going to
            #        create subtle bugs.  Better to check this when writing
            #        records.
            self.ptr.pos = pos - 1

    property start:
        """record start position on chrom/contig (0-based inclusive)"""
        def __get__(self):
            return self.ptr.pos
        def __set__(self, start):
            if start < 0:
                raise ValueError('Start coordinate must be non-negative')
            # FIXME: check start <= stop?
            #   KBJ: See above.
            self.ptr.pos = start

    property stop:
        """record stop position on chrom/contig (0-based exclusive)"""
        def __get__(self):
            return self.ptr.pos + self.ptr.rlen
        def __set__(self, stop):
            if stop < self.ptr.pos:
                raise ValueError('Stop coordinate must be greater than or equal to start')
            self.ptr.rlen = stop - self.ptr.pos

    property rlen:
        """record length on chrom/contig (typically rec.stop - rec.start unless END info is supplied)"""
        def __get__(self):
            return self.ptr.rlen
        def __set__(self, rlen):
            if rlen < 0:
                raise ValueError('Reference length must be non-negative')
            self.ptr.rlen = rlen

    property qual:
        """phred scaled quality score or None if not available"""
        def __get__(self):
            return self.ptr.qual if not bcf_float_is_missing(self.ptr.qual) else None
        def __set__(self, qual):
            if qual is not None:
                self.ptr.qual = qual
            else:
                bcf_float_set(&self.ptr.qual, bcf_float_missing)

#   property n_allele:
#       def __get__(self):
#           return self.ptr.n_allele

#   property n_sample:
#       def __get__(self):
#           return self.ptr.n_sample

    property id:
        """record identifier or None if not available"""
        def __get__(self):
            cdef bcf1_t *r = self.ptr
            if bcf_unpack(r, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking VariantRecord')
            return bcf_str_cache_get_charptr(r.d.id) if r.d.id != b'.' else None
        def __set__(self, id):
            cdef bcf1_t *r = self.ptr
            if bcf_unpack(r, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking VariantRecord')
            cdef char *idstr = NULL
            if id is not None:
                bid = force_bytes(id)
                idstr = bid
            if bcf_update_id(self.header.ptr, self.ptr, idstr) < 0:
                raise ValueError('Error updating id')

    property ref:
        """reference allele"""
        def __get__(self):
            cdef bcf1_t *r = self.ptr
            if bcf_unpack(r, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking VariantRecord')
            return charptr_to_str(r.d.allele[0]) if r.d.allele else None
        def __set__(self, ref):
            cdef bcf1_t *r = self.ptr
            if bcf_unpack(r, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking VariantRecord')
            #FIXME: Set alleles directly -- this is stupid
            if not ref:
                raise ValueError('ref allele cannot be null')
            ref = force_bytes(ref)
            if r.d.allele and r.n_allele:
                alleles = [r.d.allele[i] for i in range(r.n_allele)]
                alleles[0] = ref
            else:
                alleles = [ref]
            self.alleles = alleles

    property alleles:
        """tuple of reference allele followed by alt alleles"""
        def __get__(self):
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
        def __set__(self, values):
            cdef bcf1_t *r = self.ptr
            if bcf_unpack(r, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking VariantRecord')
            values = [force_bytes(v) for v in values]
            if b'' in values:
                raise ValueError('cannot set null allele')
            values = b','.join(values)
            if bcf_update_alleles_str(self.header.ptr, r, values) < 0:
                raise ValueError('Error updating alleles')

    property alts:
        """tuple of alt alleles"""
        def __get__(self):
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
        def __set__(self, values):
            #FIXME: Set alleles directly -- this is stupid
            cdef bcf1_t *r = self.ptr
            if bcf_unpack(r, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking VariantRecord')
            values = [force_bytes(v) for v in values]
            if b'' in values:
                raise ValueError('cannot set null alt allele')
            ref  = [r.d.allele[0] if r.d.allele and r.n_allele else b'.']
            self.alleles = ref + values

    property filter:
        """filter information (see :class:`VariantRecordFilter`)"""
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_FLT) < 0:
                raise ValueError('Error unpacking VariantRecord')
            return makeVariantRecordFilter(self)

    property info:
        """info data (see :class:`VariantRecordInfo`)"""
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_INFO) < 0:
                raise ValueError('Error unpacking VariantRecord')
            return makeVariantRecordInfo(self)

    property format:
        """sample format metadata (see :class:`VariantRecordFormat`)"""
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_FMT) < 0:
                raise ValueError('Error unpacking VariantRecord')
            return makeVariantRecordFormat(self)

    property samples:
        """sample data (see :class:`VariantRecordSamples`)"""
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_ALL) < 0:
                raise ValueError('Error unpacking VariantRecord')
            return makeVariantRecordSamples(self)

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

    property name:
        """sample name"""
        def __get__(self):
            cdef bcf_hdr_t *hdr = self.record.header.ptr
            cdef bcf1_t *r = self.record.ptr
            cdef int32_t n = bcf_hdr_nsamples(hdr)

            if self.index < 0 or self.index >= n:
                raise ValueError('invalid sample index')

            return charptr_to_str(hdr.samples[self.index])

    property allele_indices:
        """allele indices for called genotype, if present.  Otherwise None"""
        def __get__(self):
            return bcf_format_get_allele_indices(self)
        def __set__(self, values):
            self['GT'] = values
        def __del__(self):
            self['GT'] = ()

    property alleles:
        """alleles for called genotype, if present.  Otherwise None"""
        def __get__(self):
            return bcf_format_get_alleles(self)
        def __set__(self, values):
            self['GT'] = values
        def __del__(self):
            self['GT'] = ()

    property phased:
        """False if genotype is missing or any allele is unphased.  Otherwise True."""
        def __get__(self):
            return bcf_sample_get_phased(self)
        def __set__(self, value):
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
        bkey = force_bytes(key)
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

        if not refs:
            raise ValueError('Cannot retrieve reference sequence names')

        self.refs = char_array_to_tuple(refs, n, free_after=1)
        self.refmap = { r:i for i,r in enumerate(self.refs) }

    def __dealloc__(self):
        if self.ptr:
            hts_idx_destroy(self.ptr)
            self.ptr = NULL

    def fetch(self, bcf, contig, start, stop, region, reopen):
        return BCFIterator(bcf, contig, start, stop, region, reopen)


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

        if not refs:
            raise ValueError('Cannot retrieve reference sequence names')

        self.refs = char_array_to_tuple(refs, n, free_after=1)
        self.refmap = { r:i for i,r in enumerate(self.refs) }

    def __dealloc__(self):
        if self.ptr:
            tbx_destroy(self.ptr)
            self.ptr = NULL

    def fetch(self, bcf, contig, start, stop, region, reopen):
        return TabixIterator(bcf, contig, start, stop, region, reopen)


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


# Interal function to clean up after iteration stop or failure.
# This would be a nested function if it weren't a cdef function.
cdef void _stop_BCFIterator(BCFIterator self, bcf1_t *record):
    bcf_destroy1(record)

    # destroy iter so future calls to __next__ raise StopIteration
    bcf_itr_destroy(self.iter)
    self.iter = NULL


cdef class BCFIterator(BaseIterator):
    def __init__(self, VariantFile bcf, contig=None, start=None, stop=None, region=None, reopen=True):

        if not isinstance(bcf.index, BCFIndex):
            raise ValueError('bcf index required')

        cdef BCFIndex index = bcf.index
        cdef int rid, cstart, cstop
        cdef char *cregion

        if not index:
            raise ValueError('bcf index required')

        if reopen:
            bcf = bcf.copy()

        if region is not None:
            if contig is not None or start is not None or stop is not None:
                raise ValueError  # FIXME

            bregion = force_bytes(region)
            cregion = bregion
            with nogil:
                self.iter = bcf_itr_querys(index.ptr, bcf.header.ptr, cregion)
        else:
            if contig is None:
                raise ValueError  # FIXME

            try:
                rid = index.refmap[contig]
            except KeyError:
                raise('Unknown contig specified')

            if start is None:
                start = 0
            if stop is None:
                stop = MAX_POS

            cstart, cstop = start, stop

            with nogil:
                self.iter = bcf_itr_queryi(index.ptr, rid, cstart, cstop)

        # Do not fail on self.iter == NULL, since it signifies a null query.

        self.bcf = bcf
        self.index = index

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
            else:
                raise ValueError('error reading BCF file')

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

    def __init__(self, VariantFile bcf, contig=None, start=None, stop=None, region=None, reopen=True):
        if not isinstance(bcf.index, TabixIndex):
            raise ValueError('tabix index required')

        cdef TabixIndex index = bcf.index

        if not index:
            raise ValueError('bcf index required')

        if reopen:
            bcf = bcf.copy()

        if region is not None:
            if contig is not None or start is not None or stop is not None:
                raise ValueError  # FIXME

            self.iter = tbx_itr_querys(index.ptr, region)
        else:
            if contig is None:
                raise ValueError  # FIXME

            rid = index.refmap.get(contig, -1)

            if start is None:
                start = 0
            if stop is None:
                stop = MAX_POS

            self.iter = tbx_itr_queryi(index.ptr, rid, start, stop)

        # Do not fail on self.iter == NULL, since it signifies a null query.

        self.bcf = bcf
        self.index = index

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
            else:
                raise ValueError('error reading indexed VCF file')

        cdef bcf1_t *record = bcf_init1()

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


cdef class VariantFile(object):
    """*(filename, mode=None, index_filename=None, header=None, drop_samples=False)*

    A :term:`VCF`/:term:`BCF` formatted file. The file is automatically
    opened.

    *mode* should be ``r`` for reading or ``w`` for writing. The default is
    text mode (:term:`VCF`).  For binary (:term:`BCF`) I/O you should append
    ``b`` for compressed or ``u`` for uncompressed :term:`BCF` output.

    If ``b`` is present, it must immediately follow ``r`` or ``w``.  Valid
    modes are ``r``, ``w``, ``wh``, ``rb``, ``wb``, ``wbu`` and ``wb0``.
    For instance, to open a :term:`BCF` formatted file for reading, type::

        f = pysam.VariantFile('ex1.bcf','rb')

    If mode is not specified, we will try to auto-detect in the order 'rb',
    'r', thus both the following should work::

        f1 = pysam.VariantFile('ex1.bcf')
        f2 = pysam.VariantFile('ex1.vcf')

    If an index for a variant file exists (.csi or .tbi), it will be opened
    automatically.  Without an index random access to records via
    :meth:`fetch` is disabled.

    For writing, a :class:`VariantHeader` object must be provided, typically
    obtained from another :term:`VCF` file/:term:`BCF` file.
    """
    def __cinit__(self, *args, **kwargs):
        self.htsfile = NULL

    def __init__(self, *args, **kwargs):
        self.header         = None
        self.index          = None
        self.filename       = None
        self.mode           = None
        self.index_filename = None
        self.is_stream      = False
        self.is_remote      = False
        self.is_reading     = False
        self.drop_samples   = False
        self.start_offset   = -1

        self.open(*args, **kwargs)

    def __dealloc__(self):
        if self.htsfile:
            hts_close(self.htsfile)
            self.htsfile = NULL

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    property category:
        """General file format category.  One of UNKNOWN, ALIGNMENTS,
        VARIANTS, INDEX, REGIONS"""
        def __get__(self):
            if not self.htsfile:
                raise ValueError('metadata not available on closed file')
            return FORMAT_CATEGORIES[self.htsfile.format.category]

    property format:
        """File format.

        One of UNKNOWN, BINARY_FORMAT, TEXT_FORMAT, SAM, BAM,
        BAI, CRAM, CRAI, VCF, BCF, CSI, GZI, TBI, BED.
        """
        def __get__(self):
            if not self.htsfile:
                raise ValueError('metadata not available on closed file')
            return FORMATS[self.htsfile.format.format]

    property version:
        """Tuple of file format version numbers (major, minor)"""
        def __get__(self):
            if not self.htsfile:
                raise ValueError('metadata not available on closed file')
            return (self.htsfile.format.version.major,
                    self.htsfile.format.version.minor)

    property compression:
        """File compression.

        One of NONE, GZIP, BGZF, CUSTOM."""
        def __get__(self):
            if not self.htsfile:
                raise ValueError('metadata not available on closed file')
            return COMPRESSION[self.htsfile.format.compression]

    property description:
        """Vaguely human readable description of the file format"""
        def __get__(self):
            if not self.htsfile:
                raise ValueError('metadata not available on closed file')
            cdef char *desc = hts_format_description(&self.htsfile.format)
            try:
                return charptr_to_str(desc)
            finally:
                free(desc)

    def close(self):
        """closes the :class:`pysam.VariantFile`."""
        if self.htsfile:
            hts_close(self.htsfile)
            self.htsfile = NULL
        self.header = self.index = None

    property is_open:
        def __get__(self):
            """return True if VariantFile is open and in a valid state."""
            return self.htsfile != NULL

    def __iter__(self):
        if not self.is_open:
            raise ValueError('I/O operation on closed file')

        if not self.mode.startswith(b'r'):
            raise ValueError(
                'cannot iterate over Variantfile opened for writing')

        self.is_reading = 1
        return self

    def __next__(self):
        cdef int ret
        cdef bcf1_t *record = bcf_init1()

        record.pos = -1
        if self.drop_samples:
            record.max_unpack = BCF_UN_SHR

        with nogil:
            ret = bcf_read1(self.htsfile, self.header.ptr, record)

        if ret < 0:
            bcf_destroy1(record)
            if ret == -1:
                raise StopIteration
            elif ret == -2:
                raise IOError('truncated file')
            else:
                raise ValueError('Variant read failed')

        return makeVariantRecord(self.header, record)

    def copy(self):
        if not self.is_open:
            raise ValueError

        cdef VariantFile vars = VariantFile.__new__(VariantFile)
        cdef bcf_hdr_t *hdr
        cdef char *cfilename
        cdef char *cmode

        # FIXME: re-open using fd or else header and index could be invalid
        cfilename, cmode = self.filename, self.mode
        with nogil:
            vars.htsfile = hts_open(cfilename, cmode)

        if not vars.htsfile:
            raise ValueError('Cannot re-open htsfile')

        # minimize overhead by re-using header and index.  This approach is
        # currently risky, but see above for how this can be mitigated.
        vars.header         = self.header
        vars.index          = self.index

        vars.filename       = self.filename
        vars.mode           = self.mode
        vars.index_filename = self.index_filename
        vars.drop_samples   = self.drop_samples
        vars.is_stream      = self.is_stream
        vars.is_remote      = self.is_remote
        vars.is_reading     = self.is_reading
        vars.start_offset   = self.start_offset

        if self.htsfile.is_bin:
            vars.seek(self.tell())
        else:
            with nogil:
                hdr = bcf_hdr_read(vars.htsfile)
            makeVariantHeader(hdr)

        return vars

    def open(self, filename, mode='rb',
             index_filename=None,
             VariantHeader header=None,
             drop_samples=False):
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

        # close a previously opened file
        if self.is_open:
            self.close()

        if mode not in ('r','w','rb','wb', 'wh', 'wbu', 'rU', 'wb0'):
            raise ValueError('invalid file opening mode `{}`'.format(mode))

        # for htslib, wbu seems to not work
        if mode == 'wbu':
            mode = 'wb0'

        self.mode = mode = force_bytes(mode)
        self.filename = filename = encode_filename(filename)
        if index_filename is not None:
            self.index_filename = index_filename = encode_filename(index_filename)
        else:
            self.index_filename = None
        self.drop_samples = bool(drop_samples)
        self.header = None

        self.is_remote = hisremote(filename)
        self.is_stream = filename == b'-'

        if mode.startswith(b'w'):
            # open file for writing
            if index_filename is not None:
                raise ValueError('Cannot specify an index filename when writing a VCF/BCF file')

            # header structure (used for writing)
            if header:
                self.header = header.copy()
            else:
                raise ValueError('a VariantHeader must be specified')

            # open file. Header gets written to file at the same time
            # for bam files and sam files (in the latter case, the
            # mode needs to be wh)
            cfilename, cmode = filename, mode
            with nogil:
                self.htsfile = hts_open(cfilename, cmode)

            if not self.htsfile:
                raise ValueError("could not open file `{}` (mode='{}')".format((filename, mode)))

            with nogil:
                bcf_hdr_write(self.htsfile, self.header.ptr)

        elif mode.startswith(b'r'):
            # open file for reading
            if filename != b'-' and not self.is_remote and not os.path.exists(filename):
                raise IOError('file `{}` not found'.format(filename))

            cfilename, cmode = filename, mode
            with nogil:
                self.htsfile = hts_open(cfilename, cmode)

            if not self.htsfile:
                raise ValueError("could not open file `{}` (mode='{}') - is it VCF/BCF format?".format(filename, mode))

            if self.htsfile.format.format not in (bcf, vcf):
                raise ValueError("invalid file `{}` (mode='{}') - is it VCF/BCF format?".format(filename, mode))

            if self.htsfile.format.compression == bgzf:
                bgzfp = hts_get_bgzfp(self.htsfile)
                if bgzfp and bgzf_check_EOF(bgzfp) == 0:
                    warn('[%s] Warning: no BGZF EOF marker; file may be truncated'.format(filename))

            with nogil:
                hdr = bcf_hdr_read(self.htsfile)

            try:
                self.header = makeVariantHeader(hdr)
            except ValueError:
                raise ValueError("file `{}` does not have valid header (mode='{}') - is it VCF/BCF format?".format(filename, mode))

            # check for index and open if present
            if self.htsfile.format.format == bcf:
                if index_filename is not None:
                    cindex_filename = index_filename
                with nogil:
                    idx = bcf_index_load2(cfilename, cindex_filename)
                self.index = makeBCFIndex(self.header, idx)

            elif self.htsfile.format.compression == bgzf:
                if index_filename is not None:
                    cindex_filename = index_filename
                with nogil:
                    tidx = tbx_index_load2(cfilename, cindex_filename)
                self.index = makeTabixIndex(tidx)

            if not self.is_stream:
                self.start_offset = self.tell()
        else:
            raise ValueError("unknown mode {}".format(mode))

    def reset(self):
        """reset file position to beginning of file just after the header."""
        return self.seek(self.start_offset, 0)

    def seek(self, uint64_t offset):
        """move file pointer to position *offset*, see
        :meth:`pysam.VariantFile.tell`."""
        if not self.is_open:
            raise ValueError('I/O operation on closed file')
        if self.is_stream:
            raise OSError('seek not available in streams')

        cdef int64_t ret
        if self.htsfile.format.compression != no_compression:
            with nogil:
                ret = bgzf_seek(hts_get_bgzfp(self.htsfile), offset, SEEK_SET)
        else:
            with nogil:
                ret = hts_useek(self.htsfile, <int>offset, SEEK_SET)
        return ret

    def tell(self):
        """return current file position, see :meth:`pysam.VariantFile.seek`."""
        if not self.is_open:
            raise ValueError('I/O operation on closed file')
        if self.is_stream:
            raise OSError('tell not available in streams')

        cdef int64_t ret
        if self.htsfile.format.compression != no_compression:
            with nogil:
                ret = bgzf_tell(hts_get_bgzfp(self.htsfile))
        else:
            with nogil:
                ret = hts_utell(self.htsfile)
        return ret

    def fetch(self, contig=None, start=None, stop=None, region=None, reopen=False):
        """fetch records in a :term:`region` using 0-based indexing. The
        region is specified by :term:`contig`, *start* and *end*.
        Alternatively, a samtools :term:`region` string can be supplied.

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

        if not self.mode.startswith(b'r'):
            raise ValueError('cannot fetch from Variantfile opened '
                             'for writing')

        if contig is None and region is None:
            self.is_reading = 1
            bcf = self.copy() if reopen else self
            bcf.seek(self.start_offset)
            return iter(bcf)

        if not self.index:
            raise ValueError('fetch requires an index')

        self.is_reading = 1
        return self.index.fetch(self, contig, start, stop, region, reopen)

    cpdef int write(self, VariantRecord record) except -1:
        """
        write a single :class:`pysam.VariantRecord` to disk.

        returns the number of bytes written.
        """
        if not self.is_open:
            return ValueError('I/O operation on closed file')

        if not self.mode.startswith(b'w'):
            raise ValueError('cannot write to a Variantfile opened for reading')

        #if record.header is not self.header:
        #    raise ValueError('Writing records from a different VariantFile is not yet supported')

        cdef int ret

        with nogil:
            ret = bcf_write1(self.htsfile, self.header.ptr, record.ptr)

        if ret < 0:
            raise ValueError('write failed')

        return ret

    def subset_samples(self, include_samples):
        """
        Read only a subset of samples to reduce processing time and memory.
        Must be called prior to retrieving records.
        """
        if not self.is_open:
            raise ValueError('I/O operation on closed file')

        if not self.mode.startswith(b'r'):
            raise ValueError('cannot subset samples from Variantfile '
                             'opened for writing')

        if self.is_reading:
            raise ValueError('cannot subset samples after fetching records')

        self.header._subset_samples(include_samples)

        # potentially unnecessary optimization that also sets max_unpack
        if not include_samples:
            self.drop_samples = True
