# cython: embedsignature=True
# cython: profile=True
###############################################################################
###############################################################################
## Cython wrapper for htslib VCF/BCF reader/writer
###############################################################################
#
# NOTICE: This code is incomplete and preliminary.  It is nearly complete as
#         an immutable interface, but has no capability (yet) to mutate the
#         resulting data (beyond dropping all samples).  Documentation still
#         needs to be written and a unit test suite is in the works.  The
#         code is also specific to Python 2 and will require a bit of work
#         to properly adapt to Python 3.
#
# Here is a minimal example of how to use the API:
#
#     $ cat bcfview.py
#     import sys
#     from pysam.cbcf import BCFFile
#
#     bcf_in = BCFFile(sys.argv[1]) # auto-detect input format
#     bcf_out = BCFFile('-', 'w', header=bcf_in.header)
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
#     real	0m58.336s
#     user	1m6.599s
#     sys	0m3.595s
#
#     $ time bcftools view ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf |wc -l
#      1103800  # bcftools adds an extra header
#
#     real	0m55.126s
#     user	1m3.502s
#     sys	0m3.459s
#
# Here is a quick tour through the objects:
#
#     BCFFile(filename, mode=None, header=None, drop_samples=False)
#
#         htsfile:      htsFile*                             [private]
#         start_offset: BGZF offset of first record          [private]
#         filename:     filename                             [read only]
#         mode:         mode                                 [read only]
#         header:       BCFHeader object                     [read only]
#         index:        TabixIndex, BCFIndex or None         [read only]
#         drop_samples: sample information is to be ignored  [read only]
#
#         # temporary until htslib 1.2's better format metadata
#         is_bcf:       file is a bcf file                   [read only]
#         is_stream:    file is stdin/stdout                 [read only]
#         is_remote:    file is not on the local filesystem  [read only]
#
#     BCFHeader(mode)   # mode='r' for reading, mode='w' for writing
#
#         version:      VCF version
#         samples:      sequence-like access to samples
#         records:      sequence-like access to partially parsed headers
#         contigs:      mapping-like object for contig name -> id
#
#         # type info value objects are not yet implemented
#         filters:      mapping-like object for filter name -> type info
#         info:         mapping-like object for info name -> type info
#         formats:      mapping-like object for formats name -> type info
#
#     BCFRecord(...)
#
#         header:       BCFHeader object
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
#         genos:        sequence-like object of sample genotypes & attrs
#
#     BCFGeno(...)
#
#         sample:         sample name
#         sample_index:   sample index
#         allele_indices: tuple of allele indices (ref=0, alt=1..len(alts), missing=-1)
#         alleles:        tuple of alleles (missing=None)
#
#         BCFGeno is also a mapping object from formats to values
#
###############################################################################
#
# TODO list for next major sprint:
#
#   * Finish header metadata types
#   * more genotype methods
#   * unit test suite (perhaps py.test based)
#   * documentation
#
# For later sprints:
#
#   * ability to create indices
#   * mutable header and record data
#   * pickle support
#   * Python 3 support
#   * left/right locus normalization
#   * parallel iteration (like synced_bcf_reader)
#   * htslib 1.2 format info
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

from libc.string cimport strcmp

cimport cython

from cpython cimport PyBytes_Check, PyUnicode_Check
from cpython.version cimport PY_MAJOR_VERSION


__all__ = ['BCFFile', 'BCFHeader']


########################################################################
########################################################################
## Constants
########################################################################

cdef int MAX_POS = 2 << 29


########################################################################
########################################################################
## Python 3 compatibility functions
########################################################################

IS_PYTHON3 = PY_MAJOR_VERSION >= 3


cdef inline from_string_and_size(char* s, size_t length):
    if PY_MAJOR_VERSION < 3:
        return s[:length]
    else:
        return s[:length].decode('ascii')


# filename encoding (copied from lxml.etree.pyx)
cdef str FILENAME_ENCODING
FILENAME_ENCODING = sys.getfilesystemencoding()
if FILENAME_ENCODING is None:
    FILENAME_ENCODING = sys.getdefaultencoding()
if FILENAME_ENCODING is None:
    FILENAME_ENCODING = 'ascii'


cdef bytes encode_filename(object filename):
    """Make sure a filename is 8-bit encoded (or None)."""
    if filename is None:
        return None
    elif PyBytes_Check(filename):
        return filename
    elif PyUnicode_Check(filename):
        return filename.encode(FILENAME_ENCODING)
    else:
        raise TypeError('Argument must be string or unicode.')


cdef force_str(object s):
    """Return s converted to str type of current Python (bytes in Py2, unicode in Py3)"""
    if s is None:
        return None
    if PY_MAJOR_VERSION < 3:
        return s
    elif PyBytes_Check(s):
        return s.decode('ascii')
    else:
        # assume unicode
        return s


cdef bytes force_bytes(object s):
    """convert string or unicode object to bytes, assuming ascii encoding."""
    if PY_MAJOR_VERSION < 3:
        return s
    elif s is None:
        return None
    elif PyBytes_Check(s):
        return s
    elif PyUnicode_Check(s):
        return s.encode('ascii')
    else:
        raise TypeError('Argument must be string, bytes or unicode.')


cdef charptr_to_str(const char* s):
    if PY_MAJOR_VERSION < 3:
        return s
    else:
        return s.decode('ascii')


########################################################################
########################################################################
## Low level type conversion helpers
########################################################################


cdef tuple char_array_to_tuple(const char **a, int n, int free_after=0):
    if not a:
        return None
    try:
         return tuple( charptr_to_str(a[i]) for i in range(n) )
    finally:
        if free_after and a:
            free(a)


cdef bcf_array_to_object(void *data, int type, int n, int scalar=0):
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
        value = datac[:n] if datac[0] != bcf_str_missing else None
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

        if not value:
            value = None
        elif scalar and len(value) == 1:
            value = value[0]
        else:
            value = tuple(value)

    return value


cdef object bcf_info_value(const bcf_info_t *z):
    cdef char *s

    if not z:
        return None
    elif z.len == 0:
        value = True
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
            s = <char *>&z.v1.i
            value = s if not s or s[0] != bcf_str_missing else None
        else:
            raise TypeError('unsupported info type code')
    else:
        value = bcf_array_to_object(z.vptr, z.type, z.len)

    return value


cdef inline int is_gt_fmt(bcf_hdr_t *h, bcf_fmt_t *fmt):
    return strcmp(bcf_hdr_int2id(h, BCF_DT_ID, fmt.id), "GT") == 0


########################################################################
########################################################################
## BCF Header objects
########################################################################

#FIXME: passing hrec index may not be the safest approach once mutating
#       operations are allowed.
cdef class BCFHeaderRecord(object):
    property type:
        def __get__(self):
            cdef bcf_hdr_t *h = self.header.ptr
            cdef bcf_hrec_t *r = h.hrec[self.header_index]

            if r.type == BCF_HL_CTG:
                return 'CONTIG'
            elif r.type == BCF_HL_INFO:
                return 'INFO'
            elif r.type == BCF_HL_FLT:
                return 'FILTER'
            elif r.type == BCF_HL_FMT:
                return 'FORMAT'
            elif r.type == BCF_HL_STR:
                return 'STRUCTURED'
            elif r.type == BCF_HL_GEN:
                return 'GENERIC'
            else:
                return 'UNKNOWN'

    property key:
        def __get__(self):
            cdef bcf_hdr_t *h = self.header.ptr
            cdef bcf_hrec_t *r = h.hrec[self.header_index]
            return r.key if r.key else None

    property value:
        def __get__(self):
            cdef bcf_hdr_t *h = self.header.ptr
            cdef bcf_hrec_t *r = h.hrec[self.header_index]
            return r.value if r.value else None

    property attrs:
        def __get__(self):
            cdef bcf_hdr_t *h = self.header.ptr
            cdef bcf_hrec_t *r = h.hrec[self.header_index]
            cdef int i
            return tuple( (r.keys[i] if r.keys[i] else None,
                           r.vals[i] if r.vals[i] else None) for i in range(r.nkeys) )


cdef BCFHeaderRecord makeBCFHeaderRecord(BCFHeader header, int i):
    if not header:
        raise ValueError('invalid BCFHeader')

    cdef BCFHeaderRecord record = BCFHeaderRecord.__new__(BCFHeaderRecord)
    record.header = header
    record.header_index = i
    return record


cdef class BCFHeaderRecords(object):
    def __len__(self):
        return self.header.ptr.nhrec

    def __bool__(self):
        return self.header.ptr.nhrec != 0

    def __getitem__(self, index):
        cdef int32_t i = index
        if i < 0 or i >= self.header.ptr.nhrec:
            raise IndexError('invalid header record index')
        return makeBCFHeaderRecord(self.header, i)

    def __iter__(self):
        cdef int32_t i
        for i in range(self.header.ptr.nhrec):
            yield makeBCFHeaderRecord(self.header, i)

    __hash__ = None


cdef BCFHeaderRecords makeBCFHeaderRecords(BCFHeader header):
    if not header:
        raise ValueError('invalid BCFHeader')

    cdef BCFHeaderRecords records = BCFHeaderRecords.__new__(BCFHeaderRecords)
    records.header = header
    return records


cdef class BCFHeaderMetadata(object):
    def __len__(self):
        cdef bcf_hdr_t *h = self.header.ptr
        cdef bcf_idpair_t *idpair
        cdef int32_t i, n = 0

        for i in range(h.n[BCF_DT_ID]):
            idpair = h.id[BCF_DT_ID] + i
            if idpair.key and idpair.val and idpair.val.info[self.type] & 0xF != 0xF:
                n += 1

        return n

    def __bool__(self):
        cdef bcf_hdr_t *h = self.header.ptr
        cdef bcf_idpair_t *idpair
        cdef int32_t i

        for i in range(h.n[BCF_DT_ID]):
            idpair = h.id[BCF_DT_ID] + i
            if idpair.key and idpair.val and idpair.val.info[self.type] & 0xF != 0xF:
                return True

        return False

    def __getitem__(self, key):
        cdef bcf_hdr_t *h = self.header.ptr
        cdef vdict_t *d = <vdict_t *>h.dict[BCF_DT_ID]
        cdef khiter_t k = kh_get_vdict(d, key)

        if k == kh_end(d) or kh_val_vdict(d, k).info[self.type] & 0xF == 0xF:
            raise KeyError('invalid filter')

        return None

    def __iter__(self):
        cdef bcf_hdr_t *h = self.header.ptr
        cdef bcf_idpair_t *idpair
        cdef int32_t i

        for i in range(h.n[BCF_DT_ID]):
            idpair = h.id[BCF_DT_ID] + i
            if idpair.key and idpair.val and idpair.val.info[self.type] & 0xF != 0xF:
                yield idpair.key

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


cdef BCFHeaderMetadata makeBCFHeaderMetadata(BCFHeader header, int32_t type):
    if not header:
        raise ValueError('invalid BCFHeader')

    cdef BCFHeaderMetadata meta = BCFHeaderMetadata.__new__(BCFHeaderMetadata)
    meta.header = header
    meta.type = type

    return meta


cdef class BCFHeaderContigs(object):
    def __len__(self):
        cdef bcf_hdr_t *h = self.header.ptr
        assert kh_size(<vdict_t *>h.dict[BCF_DT_CTG]) == h.n[BCF_DT_CTG]
        return h.n[BCF_DT_CTG]

    def __bool__(self):
        cdef bcf_hdr_t *h = self.header.ptr
        assert kh_size(<vdict_t *>h.dict[BCF_DT_CTG]) == h.n[BCF_DT_CTG]
        return h.n[BCF_DT_CTG] != 0

    def __getitem__(self, key):
        cdef bcf_hdr_t *h = self.header.ptr
        cdef vdict_t *d = <vdict_t *>h.dict[BCF_DT_CTG]
        cdef khiter_t k = kh_get_vdict(d, key)

        if k == kh_end(d):
            raise KeyError('invalid filter')

        # FIXME: Insert correct value type
        return kh_val_vdict(d, k).id

    def __iter__(self):
        cdef bcf_hdr_t *h = self.header.ptr
        cdef vdict_t *d = <vdict_t *>h.dict[BCF_DT_CTG]
        cdef uint32_t n = kh_size(<vdict_t *>h.dict[BCF_DT_CTG])

        assert n == h.n[BCF_DT_CTG]

        for i in range(n):
            yield bcf_hdr_id2name(h, i)

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


cdef BCFHeaderContigs makeBCFHeaderContigs(BCFHeader header):
    if not header:
        raise ValueError('invalid BCFHeader')

    cdef BCFHeaderContigs contigs = BCFHeaderContigs.__new__(BCFHeaderContigs)
    contigs.header = header

    return contigs


cdef class BCFHeaderSamples(object):
    def __len__(self):
        return bcf_hdr_nsamples(self.header.ptr)

    def __bool__(self):
        return bcf_hdr_nsamples(self.header.ptr) != 0

    def __getitem__(self, index):
        cdef bcf_hdr_t *h = self.header.ptr
        cdef int32_t n = bcf_hdr_nsamples(h)
        cdef int32_t i = index

        if i < 0 or i >= n:
            raise IndexError('invalid sample index')

        return h.samples[i]

    def __iter__(self):
        cdef bcf_hdr_t *h = self.header.ptr
        cdef int32_t n = bcf_hdr_nsamples(h)
        cdef int32_t i

        for i in range(n):
            yield h.samples[i]

    def __contains__(self, key):
        cdef bcf_hdr_t *h = self.header.ptr
        cdef vdict_t *d = <vdict_t *>h.dict[BCF_DT_SAMPLE]
        cdef khiter_t k = kh_get_vdict(d, key)

        return k != kh_end(d)

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None

    #TODO: implement __richcmp__


cdef BCFHeaderSamples makeBCFHeaderSamples(BCFHeader header):
    if not header:
        raise ValueError('invalid BCFHeader')

    cdef BCFHeaderSamples samples = BCFHeaderSamples.__new__(BCFHeaderSamples)
    samples.header = header

    return samples


cdef class BCFHeader(object):
    #FIXME: Add structured proxy
    #FIXME: Add generic proxy
    #FIXME: Add immutable flag
    #FIXME: Add mutable methods

    # See makeBCFHeader for C constructor
    def __cinit__(self, mode):
        self.ptr = NULL

    # Python constructor
    def __init__(self, mode):
        if mode not in 'rw':
            raise ValueError("invalid header mode specified '{}'".format(mode))

        mode = force_bytes(mode)
        self.ptr = bcf_hdr_init(mode)

        if not self.ptr:
            raise ValueError('cannot create BCFHeader')

    def __dealloc__(self):
        if self.ptr:
            bcf_hdr_destroy(self.ptr)
            self.ptr = NULL

    def __bool__(self):
        # self.ptr == NULL should be impossible
        return self.ptr != NULL

    def copy(self):
        return makeBCFHeader(bcf_hdr_dup(self.ptr))

    property version:
        def __get__(self):
            return bcf_hdr_get_version(self.ptr)

    property samples:
        def __get__(self):
            return makeBCFHeaderSamples(self)

    property records:
        def __get__(self):
            return makeBCFHeaderRecords(self)

    property contigs:
        def __get__(self):
            return makeBCFHeaderContigs(self)

    property filters:
        def __get__(self):
            return makeBCFHeaderMetadata(self, BCF_HL_FLT)

    property info:
        def __get__(self):
            return makeBCFHeaderMetadata(self, BCF_HL_INFO)

    property formats:
        def __get__(self):
            return makeBCFHeaderMetadata(self, BCF_HL_FMT)


cdef BCFHeader makeBCFHeader(bcf_hdr_t *h):
    if not h:
        raise ValueError('cannot create BCFHeader')

    cdef BCFHeader header = BCFHeader.__new__(BCFHeader, None)
    header.ptr = h

    return header


########################################################################
########################################################################
## BCF Record objects
########################################################################

#TODO: Implement BCFFormat value class and finish Mapping interface
cdef class BCFRecordFilter(object):
    def __len__(self):
        return self.record.ptr.d.n_flt

    def __bool__(self):
        return self.record.ptr.d.n_flt != 0

    def __getitem__(self, index):
        # FIXME: Switch to value class
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i = index
        cdef int n = r.d.n_flt

        if i < 0 or i >= n:
            raise IndexError('invalid filter index')

        return bcf_hdr_int2id(h, BCF_DT_ID, r.d.flt[i])

    def __iter__(self):
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i, n = r.d.n_flt

        for i in range(n):
            yield bcf_hdr_int2id(h, BCF_DT_ID, r.d.flt[i])

#   def get(self, key, default=None):
#       """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
#       try:
#           return self[key]
#       except KeyError:
#           return default

    def __contains__(self, key):
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        return bcf_has_filter(h, r, key) == 1

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

#   def itervalues(self):
#       """D.itervalues() -> an iterator over the values of D"""
#       cdef bcf1_t *r = self.record.ptr
#       cdef bcf_info_t *info
#       cdef int i, n = r.n_info
#
#       for i in range(n):
#           info = &r.d.info[i]
#           yield bcf_info_value(info)

#   def iteritems(self):
#       """D.iteritems() -> an iterator over the (key, value) items of D"""
#       cdef bcf_hdr_t *h = self.record.header.ptr
#       cdef bcf1_t *r = self.record.ptr
#       cdef bcf_info_t *info
#       cdef int i, n = r.n_info
#
#       for i in range(n):
#           info = &r.d.info[i]
#           key = bcf_hdr_int2id(h, BCF_DT_ID, info.key)
#           value = bcf_info_value(info)
#           yield key, value

    def keys(self):
        """D.keys() -> list of D's keys"""
        return list(self)

#   def items(self):
#       """D.items() -> list of D's (key, value) pairs, as 2-tuples"""
#       return list(self.iteritems())

#   def values(self):
#       """D.values() -> list of D's values"""
#       return list(self.itervalues())

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None

    #TODO: implement __richcmp__


cdef BCFRecordFilter makeBCFRecordFilter(BCFRecord record):
    if not record:
        raise ValueError('invalid BCFRecord')

    cdef BCFRecordFilter filter = BCFRecordFilter.__new__(BCFRecordFilter)
    filter.record = record

    return filter


#TODO: Implement VCFFormat value class and finish Mapping interface
cdef class BCFRecordFormat(object):
    def __len__(self):
        return self.record.ptr.n_fmt

    def __bool__(self):
        return self.record.ptr.n_fmt != 0

    def __getitem__(self, index):
        # FIXME: Switch to value class
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i = index
        cdef int n = r.n_fmt

        if i < 0 or i >= n:
            raise IndexError('invalid format index')

        return bcf_hdr_int2id(h, BCF_DT_ID, r.d.fmt[i].id)

    def __iter__(self):
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i, n = r.n_fmt

        for i in range(n):
            yield bcf_hdr_int2id(h, BCF_DT_ID, r.d.fmt[i].id)

#   def get(self, key, default=None):
#       """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
#       try:
#           return self[key]
#       except KeyError:
#           return default

    def __contains__(self, key):
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_fmt_t *fmt = bcf_get_fmt(h, r, key)
        return fmt != NULL

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

#   def itervalues(self):
#       """D.itervalues() -> an iterator over the values of D"""
#       cdef bcf1_t *r = self.record.ptr
#       cdef bcf_info_t *info
#       cdef int i, n = r.n_info
#
#       for i in range(n):
#           info = &r.d.info[i]
#           yield bcf_info_value(info)

#   def iteritems(self):
#       """D.iteritems() -> an iterator over the (key, value) items of D"""
#       cdef bcf_hdr_t *h = self.record.header.ptr
#       cdef bcf1_t *r = self.record.ptr
#       cdef bcf_info_t *info
#       cdef int i, n = r.n_info
#
#       for i in range(n):
#           info = &r.d.info[i]
#           key = bcf_hdr_int2id(h, BCF_DT_ID, info.key)
#           value = bcf_info_value(info)
#           yield key, value

    def keys(self):
        """D.keys() -> list of D's keys"""
        return list(self)

#   def items(self):
#       """D.items() -> list of D's (key, value) pairs, as 2-tuples"""
#       return list(self.iteritems())

#   def values(self):
#       """D.values() -> list of D's values"""
#       return list(self.itervalues())

    # Mappings are not hashable by default, but subclasses can change this
    __hash__ = None

    #TODO: implement __richcmp__


cdef BCFRecordFormat makeBCFRecordFormat(BCFRecord record):
    if not record:
        raise ValueError('invalid BCFRecord')

    cdef BCFRecordFormat format = BCFRecordFormat.__new__(BCFRecordFormat)
    format.record = record

    return format


cdef class BCFRecordInfo(object):
    def __len__(self):
        return self.record.ptr.n_info

    def __bool__(self):
        return self.record.ptr.n_info != 0

    def __getitem__(self, key):
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info = bcf_get_info(h, r, key)

        if not info:
            raise KeyError('Unknown INFO field: {}'.format(key))

        return bcf_info_value(info)

    def __iter__(self):
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i, n = r.n_info

        for i in range(n):
            yield bcf_hdr_int2id(h, BCF_DT_ID, r.d.info[i].key)

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info = bcf_get_info(h, r, key)

        return info != NULL

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

    def itervalues(self):
        """D.itervalues() -> an iterator over the values of D"""
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info
        cdef int i, n = r.n_info

        for i in range(n):
            info = &r.d.info[i]
            yield bcf_info_value(info)

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info
        cdef int i, n = r.n_info

        for i in range(n):
            info = &r.d.info[i]
            key = bcf_hdr_int2id(h, BCF_DT_ID, info.key)
            value = bcf_info_value(info)
            yield key, value

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


cdef BCFRecordInfo makeBCFRecordInfo(BCFRecord record):
    if not record:
        raise ValueError('invalid BCFRecord')

    cdef BCFRecordInfo info = BCFRecordInfo.__new__(BCFRecordInfo)
    info.record = record

    return info


cdef class BCFRecordGenos(object):
    def __len__(self):
        return bcf_hdr_nsamples(self.record.header.ptr)

    def __bool__(self):
        return bcf_hdr_nsamples(self.record.header.ptr) != 0

    def __getitem__(self, key):
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int n = bcf_hdr_nsamples(h)
        cdef int sample_index
        cdef vdict_t *d
        cdef khiter_t k

        if isinstance(key, int):
            sample_index = key
        else:
            sample_index = bcf_hdr_id2int(h, BCF_DT_SAMPLE, key)
            if sample_index < 0:
                raise KeyError('invalid sample name')

        if sample_index < 0 or sample_index >= n:
            raise IndexError('invalid sample index')

        return makeBCFGeno(self.record, sample_index)

    def __iter__(self):
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int32_t i, n = bcf_hdr_nsamples(h)

        for i in range(n):
            yield h.samples[i]

#   def get(self, key, default=None):
#       """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
#       try:
#           return self[key]
#       except KeyError:
#           return default

    def __contains__(self, key):
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int n = bcf_hdr_nsamples(h)
        cdef int sample_index
        cdef vdict_t *d
        cdef khiter_t k

        if isinstance(key, int):
            sample_index = key
        else:
            sample_index = bcf_hdr_id2int(h, BCF_DT_SAMPLE, key)
            if sample_index < 0:
                raise KeyError('invalid sample name')

        return 0 <= sample_index < n

    def iterkeys(self):
        """D.iterkeys() -> an iterator over the keys of D"""
        return iter(self)

    def itervalues(self):
        """D.itervalues() -> an iterator over the values of D"""
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int32_t i, n = bcf_hdr_nsamples(h)

        for i in range(n):
            yield makeBCFGeno(self.record, i)

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int32_t i, n = bcf_hdr_nsamples(h)

        for i in range(n):
            yield h.samples[i], makeBCFGeno(self.record, i)

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


cdef BCFRecordGenos makeBCFRecordGenos(BCFRecord record):
    if not record:
        raise ValueError('invalid BCFRecord')

    cdef BCFRecordGenos genos = BCFRecordGenos.__new__(BCFRecordGenos)
    genos.record = record

    return genos


cdef class BCFRecord(object):
    def __dealloc__(self):
        if self.ptr:
            bcf_destroy1(self.ptr)
            self.ptr = NULL

    property rid:
        def __get__(self):
            return self.ptr.rid

    property chrom:
        def __get__(self):
            return bcf_hdr_id2name(self.header.ptr, self.ptr.rid)

    property contig:
        def __get__(self):
            return bcf_hdr_id2name(self.header.ptr, self.ptr.rid)

    property pos:
        def __get__(self):
            return self.ptr.pos + 1

    property start:
        def __get__(self):
            return self.ptr.pos

    property stop:
        def __get__(self):
            return self.ptr.pos + self.ptr.rlen

    property rlen:
        def __get__(self):
            return self.ptr.rlen

    property qual:
        def __get__(self):
            return self.ptr.qual if not bcf_float_is_missing(self.ptr.qual) else None

#   property n_info:
#       def __get__(self):
#            if bcf_unpack(self.ptr, BCF_UN_INFO) < 0:
#               raise ValueError('Error unpacking BCFRecord')
#           return self.ptr.n_info

#   property n_allele:
#       def __get__(self):
#           return self.ptr.n_allele

#   property n_fmt:
#       def __get__(self):
#           return self.ptr.n_fmt

#   property n_sample:
#       def __get__(self):
#           return self.ptr.n_sample

#   property shared:
#       def __get__(self):
#           return self.ptr.shared.s

#   property indiv:
#       def __get__(self):
#           return self.ptr.indiv.s

#   property n_flt:
#       def __get__(self):
#           if bcf_unpack(self.ptr, BCF_UN_FLT) < 0:
#               raise ValueError('Error unpacking BCFRecord')
#           return self.ptr.d.n_flt

    property id:
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking BCFRecord')
            id = self.ptr.d.id
            return id if id != b'.' else None

    property ref:
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking BCFRecord')
            return self.ptr.d.allele[0] if self.ptr.d.allele else None

    property alleles:
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking BCFRecord')
            if not self.ptr.d.allele:
                return None
            return tuple(self.ptr.d.allele[i] for i in range(self.ptr.n_allele))

    property alts:
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking BCFRecord')
            if self.ptr.n_allele < 2 or not self.ptr.d.allele:
                return None
            return tuple(self.ptr.d.allele[i] for i in range(1,self.ptr.n_allele))

    property filter:
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_FLT) < 0:
                raise ValueError('Error unpacking BCFRecord')
            return makeBCFRecordFilter(self)

    property info:
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_INFO) < 0:
                raise ValueError('Error unpacking BCFRecord')
            return makeBCFRecordInfo(self)

    property format:
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_FMT) < 0:
                raise ValueError('Error unpacking BCFRecord')
            return makeBCFRecordFormat(self)

    property genos:
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_IND) < 0:
                raise ValueError('Error unpacking BCFRecord')
            return makeBCFRecordGenos(self)


cdef BCFRecord makeBCFRecord(BCFHeader header, bcf1_t *r):
    if not header:
        raise ValueError('invalid BCFHeader')

    if not r:
        raise ValueError('cannot create BCFRecord')

    cdef BCFRecord record = BCFRecord.__new__(BCFRecord)
    record.header = header
    record.ptr = r

    return record


########################################################################
########################################################################
## BCF Genotype object
########################################################################


cdef class BCFGeno(object):
    property sample:
        def __get__(self):
            cdef bcf_hdr_t *h = self.record.header.ptr
            cdef bcf1_t *r = self.record.ptr
            cdef int32_t n = bcf_hdr_nsamples(h)

            if self.sample_index < 0 or self.sample_index >= n:
                raise ValueError('invalid sample index')

            return h.samples[self.sample_index]

    property allele_indices:
        def __get__(self):
            cdef bcf_hdr_t *h = self.record.header.ptr
            cdef bcf1_t *r = self.record.ptr
            cdef int32_t n = bcf_hdr_nsamples(h)

            if self.sample_index < 0 or self.sample_index >= n or not r.n_fmt:
                return None

            cdef bcf_fmt_t *fmt0 = r.d.fmt
            cdef int gt0 = strcmp(bcf_hdr_int2id(h, BCF_DT_ID, fmt0.id), "GT") == 0

            if not gt0 or not fmt0.n:
                return None

            cdef int8_t  *data8
            cdef int16_t *data16
            cdef int32_t *data32
            alleles = []

            if fmt0.type == BCF_BT_INT8:
                data8 = <int8_t *>(fmt0.p + self.sample_index * fmt0.size)
                for i in range(fmt0.n):
                    if data8[i] == bcf_int8_vector_end:
                        break
                    alleles.append( (data8[i] >> 1) - 1 )
            elif fmt0.type == BCF_BT_INT16:
                data16 = <int16_t *>(fmt0.p + self.sample_index * fmt0.size)
                for i in range(fmt0.n):
                    if data16[i] == bcf_int16_vector_end:
                        break
                    alleles.append( (data16[i] >> 1) - 1 )
            elif fmt0.type == BCF_BT_INT32:
                data32 = <int32_t *>(fmt0.p + self.sample_index * fmt0.size)
                for i in range(fmt0.n):
                    if data32[i] == bcf_int32_vector_end:
                        break
                    alleles.append( (data32[i] >> 1) - 1 )

            return tuple(alleles)

    property alleles:
        def __get__(self):
            cdef bcf_hdr_t *h = self.record.header.ptr
            cdef bcf1_t *r = self.record.ptr
            cdef int32_t nsamples = bcf_hdr_nsamples(h)
            cdef int32_t nalleles = r.n_allele

            if self.sample_index < 0 or self.sample_index >= nsamples or not r.n_fmt:
                return None

            cdef bcf_fmt_t *fmt0 = r.d.fmt
            cdef int gt0 = strcmp(bcf_hdr_int2id(h, BCF_DT_ID, fmt0.id), "GT") == 0

            if not gt0 or not fmt0.n:
                return None

            cdef int32_t  a
            cdef int8_t  *data8
            cdef int16_t *data16
            cdef int32_t *data32
            alleles = []

            if fmt0.type == BCF_BT_INT8:
                data8 = <int8_t *>(fmt0.p + self.sample_index * fmt0.size)
                for i in range(fmt0.n):
                    if data8[i] == bcf_int8_vector_end:
                        break
                    a = (data8[i] >> 1) - 1
                    alleles.append(r.d.allele[a] if 0 <= a < nalleles else None)
            elif fmt0.type == BCF_BT_INT16:
                data16 = <int16_t *>(fmt0.p + self.sample_index * fmt0.size)
                for i in range(fmt0.n):
                    if data16[i] == bcf_int16_vector_end:
                        break
                    a = (data16[i] >> 1) - 1
                    alleles.append(r.d.allele[a] if 0 <= a < nalleles else None)
            elif fmt0.type == BCF_BT_INT32:
                data32 = <int32_t *>(fmt0.p + self.sample_index * fmt0.size)
                for i in range(fmt0.n):
                    if data32[i] == bcf_int32_vector_end:
                        break
                    a = (data32[i] >> 1) - 1
                    alleles.append(r.d.allele[a] if 0 <= a < nalleles else None)

            return tuple(alleles)

    def __len__(self):
        return self.record.ptr.n_fmt

    def __bool__(self):
        return self.record.ptr.n_fmt != 0

    def __getitem__(self, key):
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_fmt_t *fmt
        cdef int index

        if isinstance(key, int):
            index = key
            if index < 0 or index >= r.n_fmt:
                raise IndexError('invalid format index')
            fmt = r.d.fmt + index
        else:
            fmt = bcf_get_fmt(h, r, key)

        if not fmt:
            raise KeyError('invalid format requested')

        if is_gt_fmt(h, fmt):
            return self.alleles
        elif fmt.p and fmt.n and fmt.size:
            return bcf_array_to_object(fmt.p + self.sample_index * fmt.size, fmt.type, fmt.n, scalar=1)
        else:
            return None

    def __iter__(self):
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i, n = r.n_fmt

        for i in range(n):
            yield bcf_hdr_int2id(h, BCF_DT_ID, r.d.fmt[i].id)

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        cdef bcf_hdr_t *h = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_fmt_t *fmt = bcf_get_fmt(h, r, key)
        return fmt != NULL

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


cdef BCFGeno makeBCFGeno(BCFRecord record, int32_t sample_index):
    if not record or sample_index < 0:
        raise ValueError('cannot create BCFGeno')

    cdef BCFGeno geno = BCFGeno.__new__(BCFGeno)
    geno.record = record
    geno.sample_index = sample_index

    return geno


########################################################################
########################################################################
## Index objects
########################################################################


cdef class BaseIndex(object):
    pass


cdef class BCFIndex(object):
    def __init__(self):
        self.refs = ()
        self.refmap = {}

        if not self.ptr:
            raise ValueError('Invalid index object')

        cdef const char **refs
        cdef int n

        refs = bcf_index_seqnames(self.ptr, self.header.ptr, &n)

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


cdef BCFIndex makeBCFIndex(BCFHeader header, hts_idx_t *idx):
    if not idx:
        return None

    if not header:
        raise ValueError('invalid BCFHeader')

    cdef BCFIndex index = BCFIndex.__new__(BCFIndex)
    index.header = header
    index.ptr = idx
    index.__init__()

    return index


cdef class TabixIndex(BaseIndex):
    def __init__(self):
        self.refs = ()
        self.refmap = {}

        if not self.ptr:
            raise ValueError('Invalid index object')

        cdef const char **refs
        cdef int n

        refs = tbx_seqnames(self.ptr, &n)

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


cdef class BCFIterator(BaseIterator):
    def __init__(self, BCFFile bcf, contig=None, start=None, stop=None, region=None, reopen=True):

        if not isinstance(bcf.index, BCFIndex):
            raise ValueError('bcf index required')

        cdef BCFIndex index = bcf.index

        if not index:
            raise ValueError('bcf index required')

        if reopen:
            bcf = bcf.copy()

        if region is not None:
            if contig is not None or start is not None or stop is not None:
                raise ValueError  # FIXME

            self.iter = bcf_itr_querys(index.ptr, bcf.header.ptr, region)
        else:
            if contig is None:
                raise ValueError  # FIXME

            rid = index.refmap.get(contig, -1)

            if rid is not None and (start is not None or stop is not None):
                raise ValueError  # FIXME

            if start is None:
                start = 0
            if stop is None:
                stop = MAX_POS

            self.iter = bcf_itr_queryi(index.ptr, rid, start, stop)

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

        cdef int ret = bcf_itr_next(self.bcf.htsfile, self.iter, <void *>record)

        if ret < 0:
            bcf_destroy1(record)
            bcf_itr_destroy(self.iter)
            self.iter = NULL
            raise StopIteration

        return makeBCFRecord(self.bcf.header, record)


cdef class TabixIterator(BaseIterator):
    def __cinit__(self, *args, **kwargs):
        self.line_buffer.l = 0
        self.line_buffer.m = 0
        self.line_buffer.s = NULL

    def __init__(self, BCFFile bcf, contig=None, start=None, stop=None, region=None, reopen=True):
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

            if rid is not None and (start is not None or stop is not None):
                raise ValueError  # FIXME

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

        cdef int ret = tbx_itr_next(self.bcf.htsfile, self.index.ptr, self.iter, <void *>&self.line_buffer)

        if ret < 0:
            tbx_itr_destroy(self.iter)
            self.iter = NULL
            raise StopIteration

        cdef bcf1_t *record = bcf_init1()

        record.pos = -1
        if self.bcf.drop_samples:
            record.max_unpack = BCF_UN_SHR

        ret = vcf_parse1(&self.line_buffer, self.bcf.header.ptr, record)

        if ret < 0:
            bcf_destroy1(record)
            raise ValueError('error in vcf_parse')

        return makeBCFRecord(self.bcf.header, record)


########################################################################
########################################################################
## BCF File
########################################################################


cdef class BCFFile(object):
    """*(filename, mode=None, header=None, drop_samples=False)*

    A :term:`VCF`/:term:`BCF` formatted file. The file is automatically
    opened.

    *mode* should be ``r`` for reading or ``w`` for writing. The default is
    text mode (:term:`VCF`).  For binary (:term:`BCF`) I/O you should append
    ``b`` for compressed or ``u`` for uncompressed :term:`BCF` output.  Use
    ``h`` to output header information in text (:term:`TAM`) mode.

    If ``b`` is present, it must immediately follow ``r`` or ``w``.  Valid
    modes are ``r``, ``w``, ``wh``, ``rb``, ``wb``, ``wbu`` and ``wb0``.
    For instance, to open a :term:`BCF` formatted file for reading, type::

        f = pysam.BCFFile('ex1.bcf','rb')

    If mode is not specified, we will try to auto-detect in the order 'rb',
    'r', thus both the following should work::

        f1 = pysam.BCFFile('ex1.bcf')
        f2 = pysam.BCFFile('ex1.vcf')

    If an index for a BCF file exists (.csi), it will be opened
    automatically.  Without an index random access to records via
    :meth:`fetch` is disabled.

    For writing, a :class:`BCFHeader` object of a :term:`VCF`
    file/:term:`BCF` file must be provided.
    """
    def __cinit__(self, *args, **kwargs):
        self.htsfile = NULL

    def __init__(self, *args, **kwargs):
        self.header       = None
        self.index        = None
        self.filename     = None
        self.mode         = None
        self.is_bcf       = False
        self.is_stream    = False
        self.is_remote    = False
        self.drop_samples = False
        self.start_offset = -1

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

    def __iter__(self):
        if not self.is_open:
            raise ValueError('I/O operation on closed file')
        return self

    def __next__(self):
        cdef int ret
        cdef bcf1_t *record = bcf_init1()

        record.pos = -1
        if self.drop_samples:
            record.max_unpack = BCF_UN_SHR

        ret = bcf_read1(self.htsfile, self.header.ptr, record)

        if ret == -1:
            bcf_destroy1(record)
            raise StopIteration
        elif ret == -2:
            bcf_destroy1(record)
            raise IOError('truncated file')

        return makeBCFRecord(self.header, record)

    def copy(self):
        if not self.is_open:
            raise ValueError

        cdef BCFFile bcf = BCFFile.__new__(BCFFile)

        # FIXME: re-open using fd or else header and index could be invalid
        #        adding an hts_reopen/hts_dup function will be easy once
        #        htslib 1.2 is a bit further along
        bcf.htsfile      = hts_open(self.filename, self.mode)

        if not bcf.htsfile:
            raise ValueError('Cannot re-open htsfile')

        # minimize overhead by re-using header and index.  This approach is
        # currently risky, but see above for how this can be mitigated.
        bcf.header       = self.header
        bcf.index        = self.index

        bcf.filename     = self.filename
        bcf.mode         = self.mode
        bcf.drop_samples = self.drop_samples
        bcf.is_bcf       = self.is_bcf
        bcf.is_stream    = self.is_stream
        bcf.is_remote    = self.is_remote
        bcf.start_offset = self.start_offset

        # FIXME: seeking here does not currently work for text files.
        # Falling back to re-reading the header for now.  There is no point
        # in debugging this now, but I will look into this more when htslib
        # 1.2 is more stable.

        if self.is_bcf:
            bcf.seek(self.tell())
        else:
            makeBCFHeader(bcf_hdr_read(bcf.htsfile))

        return bcf

    def close(self):
        """closes the :class:`pysam.BCFFile`."""
        if self.htsfile:
            hts_close(self.htsfile)
            self.htsfile = NULL
        self.header = self.index = None

    property is_open:
        def __get__(self):
            """return True if BCFFile is open and in a valid state."""
            return self.htsfile != NULL

    def open(self, filename, mode=None, BCFHeader header=None, drop_samples=False):
        """open a vcf/bcf file.

        If open is called on an existing BCFFile, the current file will be
        closed and a new file will be opened.
        """
        # close a previously opened file
        if self.is_open:
            self.close()

        # read mode autodetection
        if mode is None:
            try:
                self.open(filename, 'rb', header=header)
                return
            except ValueError, msg:
                pass

            self.open(filename, 'r', header=header)
            return

        if mode not in ('r','w','rb','wb', 'wh', 'wbu', 'rU', 'wb0'):
            raise ValueError('invalid file opening mode `{}`'.format(mode))

        mode = mode.encode('ascii')

        # for htslib, wbu seems to not work
        if mode == b'wbu':
            mode = b'wb0'

        self.mode = mode
        self.filename = filename = encode_filename(filename)
        self.drop_samples = bool(drop_samples)

        # FIXME: Use htsFormat when it is available
        self.is_bcf = filename.endswith('.bcf')
        self.is_remote = filename.startswith(b'http:') or filename.startswith(b'ftp:')
        self.is_stream = filename == b'-'

        if mode[0] == b'w':
            # open file for writing

            # header structure (used for writing)
            if header:
                self.header = header.copy()
            else:
                raise ValueError('a BCFHeader must be specified')

            # open file. Header gets written to file at the same time for bam files
            # and sam files (in the latter case, the mode needs to be wh)
            self.htsfile = hts_open(filename, mode)

            if not self.htsfile:
                raise ValueError("could not open file `{}` (mode='{}')".format((filename, mode)))

            bcf_hdr_write(self.htsfile, self.header.ptr)

        elif mode[0] == b'r':
            # open file for reading
            if filename != b'-' and not self.is_remote and not os.path.exists(filename):
                raise IOError('file `{}` not found'.format(filename))

            self.htsfile = hts_open(filename, mode)

            if not self.htsfile:
                raise ValueError("could not open file `{}` (mode='{}') - is it VCF/BCF format?".format((filename, mode)))

            self.header = makeBCFHeader(bcf_hdr_read(self.htsfile))

            if not self.header:
                raise ValueError("file `{}` does not have valid header (mode='{}') - is it BCF format?".format((filename, mode)))

            # check for index and open if present
            if self.is_bcf:
                self.index = makeBCFIndex(self.header, bcf_index_load(filename))
            else:
                self.index = makeTabixIndex(tbx_index_load(filename + '.tbi'))

            if not self.is_stream:
                self.start_offset = self.tell()

    def reset(self):
        """reset file position to beginning of file just after the header."""
        return self.seek(self.start_offset, 0)

    def seek(self, uint64_t offset):
        """move file pointer to position *offset*, see :meth:`pysam.BCFFile.tell`."""
        if not self.is_open:
            raise ValueError('I/O operation on closed file')
        if self.is_stream:
            raise OSError('seek not available in streams')

        return bgzf_seek(hts_get_bgzfp(self.htsfile), offset, SEEK_SET)

    def tell(self):
        """
        return current file position, see :meth:`pysam.BCFFile.seek`.
        """
        if not self.is_open:
            raise ValueError('I/O operation on closed file')
        if self.is_stream:
            raise OSError('tell not available in streams')

        return bgzf_tell(hts_get_bgzfp(self.htsfile))

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

        Note that a :term:`VCF` file without a tabix index (.tbi) or a
        :term:`BCF` file without a CSI index can only be read sequentially.
        """
        if not self.is_open:
            raise ValueError('I/O operation on closed file')

        if contig is None and region is None:
            bcf = self.copy() if reopen else self
            bcf.seek(self.start_offset)
            return bcf

        if not self.index:
            raise ValueError('fetch requires an index')

        return self.index.fetch(self, contig, start, stop, region, reopen)

    cpdef int write(self, BCFRecord record) except -1:
        """
        write a single :class:`pysam.BCFRecord` to disk.

        returns the number of bytes written.
        """
        if not self.is_open:
            return 0

        cdef int ret = bcf_write1(self.htsfile, self.header.ptr, record.ptr)

        if ret < 0:
            raise ValueError('write failed')

        return ret
