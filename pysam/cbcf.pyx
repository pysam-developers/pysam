# cython: embedsignature=True
# cython: profile=True
###############################################################################
###############################################################################
## Cython wrapper for htslib VCF/BCF reader/writer
###############################################################################
#
# NOTICE: This code is incomplete and preliminary.  It does offer a nearly
#         complete immutable Pythonic interface to VCF/BCF metadata and data
#         with reading and writing capability, but has no capability (yet)
#         to mutate the resulting data (beyond dropping all samples).
#         Documentation still needs to be written and a unit test suite is
#         in the works.  The code is also superficially specific to Python 2
#         and will require a bit of work to properly adapt to Python 3.
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

cimport cython

from cpython cimport PyBytes_Check, PyUnicode_Check
from cpython.version cimport PY_MAJOR_VERSION

__all__ = ['VariantFile', 'VariantHeader']

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

from pysam.cutils cimport force_bytes, force_str, charptr_to_str
from pysam.cutils cimport encode_filename, from_string_and_size


########################################################################
########################################################################
## Low level type conversion helpers
########################################################################


cdef tuple char_array_to_tuple(const char **a, int n, int free_after=0):
    if not a:
        return None
    try:
         return tuple(charptr_to_str(a[i]) for i in range(n))
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


cdef inline int is_gt_fmt(bcf_hdr_t *hdr, bcf_fmt_t *fmt):
    return strcmp(bcf_hdr_int2id(hdr, BCF_DT_ID, fmt.id), "GT") == 0


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
            return r.key if r.key else None

    property value:
        """header value.  Set only for generic lines, None for FILTER/INFO, etc."""
        def __get__(self):
            cdef bcf_hrec_t *r = self.ptr
            return r.value if r.value else None

    property attrs:
        """sequence of additional header attributes"""
        def __get__(self):
            cdef bcf_hrec_t *r = self.ptr
            cdef int i
            return tuple( (r.keys[i] if r.keys[i] else None,
                           r.vals[i] if r.vals[i] else None) for i in range(r.nkeys) )

    def __len__(self):
        cdef bcf_hrec_t *r = self.ptr
        return r.nkeys

    def __bool__(self):
        cdef bcf_hrec_t *r = self.ptr
        cdef int i
        for i in range(r.nkeys):
            yield r.keys[i]

    def __getitem__(self, key):
        """get attribute value"""
        cdef bcf_hrec_t *r = self.ptr
        cdef int i
        for i in range(r.nkeys):
            if r.keys[i] and r.keys[i] == key:
                return r.vals[i] if r.vals[i] else None
        raise KeyError('cannot find metadata key')

    def __iter__(self):
        cdef bcf_hrec_t *r = self.ptr
        cdef int i
        for i in range(r.nkeys):
            if r.keys[i]:
                yield r.keys[i]

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
                yield r.vals[i] if r.vals[i] else None

    def iteritems(self):
        """D.iteritems() -> an iterator over the (key, value) items of D"""
        cdef bcf_hrec_t *r = self.ptr
        cdef int i
        for i in range(r.nkeys):
            if r.keys[i]:
                yield r.keys[i], r.vals[i] if r.vals[i] else None

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
    """filter, info or format metadata record from a :class:`VariantHeader` object"""
    property name:
        """metadata name"""
        def __get__(self):
            cdef bcf_hdr_t *hdr = self.header.ptr
            return hdr.id[BCF_DT_ID][self.id].key

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
            if not bcf_hdr_idinfo_exists(hdr, self.type, self.id) or self.type == BCF_HL_FLT:
                return None
            return VALUE_TYPES[bcf_hdr_id2type(hdr, self.type, self.id)]

    property description:
        """metadata description (or None if not set)"""
        def __get__(self):
            descr = self.record.get('Description')
            if descr:
                descr = descr.strip('"')
            return descr

    property record:
        """:class:`VariantHeaderRecord` associated with this :class:`VariantMetadata` object"""
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

            items = [('ID', id), ('Number', number), ('Type', type), ('Description', description)]

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
        cdef khiter_t k = kh_get_vdict(d, key)

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
            return hdr.id[BCF_DT_CTG][self.id].key

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
        cdef khiter_t k = kh_get_vdict(d, key)

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
            yield bcf_hdr_id2name(hdr, i)

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

        return hdr.samples[i]

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef int32_t i, n = bcf_hdr_nsamples(hdr)

        for i in range(n):
            yield hdr.samples[i]

    def __contains__(self, key):
        cdef bcf_hdr_t *hdr = self.header.ptr
        cdef vdict_t *d = <vdict_t *>hdr.dict[BCF_DT_SAMPLE]
        cdef khiter_t k = kh_get_vdict(d, key)

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
            return bcf_hdr_get_version(self.ptr)

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
        """
        alt metadata (:class:`dict` ID->record).  The data returned just a snapshot of alt records,
        is created every time the property is requested, and modifications will not be reflected
        in the header metadata and vice versa.

        i.e. it is just a dict that reflects the state of alt records at the time it is created.
        """
        def __get__(self):
            return { record['ID']:record for record in self.records if record.key.upper() == 'ALT' }


    # only safe to do when opening an htsfile
    cdef _subset_samples(self, include_samples):
        keep_samples    = set(self.samples)
        include_samples = set(include_samples)
        missing_samples = include_samples - keep_samples
        keep_samples   &= include_samples

        if missing_samples:
            # FIXME: add specialized exception with payload
            raise ValueError('missing {:d} requested samples'.format(len(missing_samples)))

        keep_samples = ','.join(keep_samples)
        cdef char *keep = <char *>keep_samples if keep_samples else NULL
        cdef ret = bcf_hdr_set_samples(self.ptr, keep, 0)

        if ret != 0:
            raise ValueError('bcf_hdr_set_samples failed: ret = {}'.format(ret))

    def __str__(self):
        cdef int hlen
        cdef char *hstr = bcf_hdr_fmt_text(self.ptr, 0, &hlen)

        ret = hstr[:hlen]
        free(hstr)
        return force_str(hstr)

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
        if bcf_hdr_append(self.ptr, line) < 0:
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
            hrec.key = strdup(key)

            if value is not None:
                hrec.value = strdup(value)
            else:
                for key, value in items:
                    bcf_hrec_add_key(hrec, key, len(key))

                    value = str(value)
                    quoted = strpbrk(value, ' ;,"\t<>') != NULL
                    bcf_hrec_set_val(hrec, hrec.nkeys-1, value, len(value), quoted)
        except:
            bcf_hrec_destroy(hrec)
            raise

        bcf_hdr_add_hrec(self.ptr, hrec)

        if self.ptr.dirty:
            bcf_hdr_sync(self.ptr)

    def add_sample(self, name):
        """Add a new sample to this header"""
        if bcf_hdr_add_sample(self.ptr, name) < 0:
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
    """mapping from filter index or name to :class:`VariantMetadata` object for filters set on a :class:`VariantRecord` object."""

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

            id = bcf_hdr_id2int(hdr, BCF_DT_ID, key)

            if not bcf_hdr_idinfo_exists(hdr, BCF_HL_FLT, id) or not bcf_has_filter(hdr, self.record.ptr, key):
                raise KeyError('Invalid filter')

        return makeVariantMetadata(self.record.header, BCF_HL_FLT, id)

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i, n = r.d.n_flt

        for i in range(n):
            yield bcf_hdr_int2id(hdr, BCF_DT_ID, r.d.flt[i])

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        return bcf_has_filter(hdr, r, key) == 1

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
    """mapping from format name or index to :class:`VariantMetadata` object for formats present in a :class:`VariantRecord` object."""

    def __len__(self):
        return self.record.ptr.n_fmt

    def __bool__(self):
        return self.record.ptr.n_fmt != 0

    def __getitem__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_fmt_t *fmt
        cdef int index
        cdef int n = r.n_fmt

        if isinstance(key, int):
            index = key
            if index < 0 or index >= n:
                raise IndexError('invalid format index')
            fmt = &r.d.fmt[index]
        else:
            fmt = bcf_get_fmt(hdr, r, key)
            if not fmt:
                raise KeyError('unknown format')

        return makeVariantMetadata(self.record.header, BCF_HL_FMT, fmt.id)

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i, n = r.n_fmt

        for i in range(n):
            yield bcf_hdr_int2id(hdr, BCF_DT_ID, r.d.fmt[i].id)

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_fmt_t *fmt = bcf_get_fmt(hdr, r, key)
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


cdef VariantRecordFormat makeVariantRecordFormat(VariantRecord record):
    if not record:
        raise ValueError('invalid VariantRecord')

    cdef VariantRecordFormat format = VariantRecordFormat.__new__(VariantRecordFormat)
    format.record = record

    return format


#TODO: Add a getmeta method to return the corresponding VariantMetadata?
cdef class VariantRecordInfo(object):
    """mapping from info metadata name to value for info data present in a :class:`VariantRecord` object."""

    def __len__(self):
        return self.record.ptr.n_info

    def __bool__(self):
        return self.record.ptr.n_info != 0

    def __getitem__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info = bcf_get_info(hdr, r, key)

        if not info:
            raise KeyError('Unknown INFO field: {}'.format(key))

        return bcf_info_value(info)

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i, n = r.n_info

        for i in range(n):
            yield bcf_hdr_int2id(hdr, BCF_DT_ID, r.d.info[i].key)

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info = bcf_get_info(hdr, r, key)

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
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_info_t *info
        cdef int i, n = r.n_info

        for i in range(n):
            info = &r.d.info[i]
            key = bcf_hdr_int2id(hdr, BCF_DT_ID, info.key)
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
            sample_index = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, key)
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
            yield hdr.samples[i]

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
            sample_index = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, key)
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
            yield hdr.samples[i], makeVariantRecordSample(self.record, i)

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

    cdef VariantRecordSamples samples = VariantRecordSamples.__new__(VariantRecordSamples)
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
            return bcf_hdr_id2name(self.header.ptr, self.ptr.rid)
        def __set__(self, chrom):
            cdef vdict_t *d = <vdict_t*>self.header.ptr.dict[BCF_DT_CTG]
            cdef khint_t k = kh_get_vdict(d, chrom)
            if k == kh_end(d):
                raise ValueError('Invalid chromosome/contig')
            self.ptr.rid = kh_val_vdict(d, k).id

    property contig:
        """chromosome/contig name"""
        def __get__(self):
            return bcf_hdr_id2name(self.header.ptr, self.ptr.rid)
        def __set__(self, chrom):
            cdef vdict_t *d = <vdict_t*>self.header.ptr.dict[BCF_DT_CTG]
            cdef khint_t k = kh_get_vdict(d, chrom)
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
            self.ptr.pos = pos - 1

    property start:
        """record start position on chrom/contig (0-based inclusive)"""
        def __get__(self):
            return self.ptr.pos
        def __set__(self, start):
            if start < 0:
                raise ValueError('Start coordinate must be non-negative')
            # FIXME: check start <= stop?
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
                memcpy(&self.ptr.qual, &bcf_float_missing, 4)

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
#               raise ValueError('Error unpacking VariantRecord')
#           return self.ptr.d.n_flt

    property id:
        """record identifier or None if not available"""
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking VariantRecord')
            id = self.ptr.d.id
            return id if id != b'.' else None
        def __set__(self, id):
            cdef char *idstr = NULL
            if id is not None:
                idstr = id
            if bcf_update_id(self.header.ptr, self.ptr, idstr) < 0:
                raise ValueError('Error updating id')

    property ref:
        """reference allele"""
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking VariantRecord')
            return self.ptr.d.allele[0] if self.ptr.d.allele else None
        def __set__(self, ref):
            alleles = list(self.alleles)
            alleles[0] = ref
            self.alleles = alleles

    property alleles:
        """tuple of reference allele followed by alt alleles"""
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking VariantRecord')
            if not self.ptr.d.allele:
                return None
            return tuple(self.ptr.d.allele[i] for i in range(self.ptr.n_allele))
        def __set__(self, values):
            if bcf_unpack(self.ptr, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking VariantRecord')
            values = ','.join(values)
            if bcf_update_alleles_str(self.header.ptr, self.ptr, values) < 0:
                raise ValueError('Error updating alleles')

    property alts:
        """tuple of alt alleles"""
        def __get__(self):
            if bcf_unpack(self.ptr, BCF_UN_STR) < 0:
                raise ValueError('Error unpacking VariantRecord')
            if self.ptr.n_allele < 2 or not self.ptr.d.allele:
                return None
            return tuple(self.ptr.d.allele[i] for i in range(1,self.ptr.n_allele))
        def __set__(self, alts):
            alleles = [self.ref]
            alleles.extend(alts)
            self.alleles = alleles

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
            if bcf_unpack(self.ptr, BCF_UN_IND) < 0:
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

        ret = line.s[:line.l]
        ret = force_str(ret)

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
       Provides data accessors for genotypes and a mapping interface from format name to values.
    """

    property name:
        """sample name"""
        def __get__(self):
            cdef bcf_hdr_t *hdr = self.record.header.ptr
            cdef bcf1_t *r = self.record.ptr
            cdef int32_t n = bcf_hdr_nsamples(hdr)

            if self.index < 0 or self.index >= n:
                raise ValueError('invalid sample index')

            return hdr.samples[self.index]

    property allele_indices:
        """allele indices for called genotype, if present.  Otherwise None"""
        def __get__(self):
            cdef bcf_hdr_t *hdr = self.record.header.ptr
            cdef bcf1_t *r = self.record.ptr
            cdef int32_t n = bcf_hdr_nsamples(hdr)

            if self.index < 0 or self.index >= n or not r.n_fmt:
                return None

            cdef bcf_fmt_t *fmt0 = r.d.fmt
            cdef int gt0 = is_gt_fmt(hdr, fmt0)

            if not gt0 or not fmt0.n:
                return None

            cdef int8_t  *data8
            cdef int16_t *data16
            cdef int32_t *data32
            alleles = []

            if fmt0.type == BCF_BT_INT8:
                data8 = <int8_t *>(fmt0.p + self.index * fmt0.size)
                for i in range(fmt0.n):
                    if data8[i] == bcf_int8_vector_end:
                        break
                    alleles.append( (data8[i] >> 1) - 1 )
            elif fmt0.type == BCF_BT_INT16:
                data16 = <int16_t *>(fmt0.p + self.index * fmt0.size)
                for i in range(fmt0.n):
                    if data16[i] == bcf_int16_vector_end:
                        break
                    alleles.append( (data16[i] >> 1) - 1 )
            elif fmt0.type == BCF_BT_INT32:
                data32 = <int32_t *>(fmt0.p + self.index * fmt0.size)
                for i in range(fmt0.n):
                    if data32[i] == bcf_int32_vector_end:
                        break
                    alleles.append( (data32[i] >> 1) - 1 )

            return tuple(alleles)

    property alleles:
        """alleles for called genotype, if present.  Otherwise None"""
        def __get__(self):
            cdef bcf_hdr_t *hdr = self.record.header.ptr
            cdef bcf1_t *r = self.record.ptr
            cdef int32_t nsamples = bcf_hdr_nsamples(hdr)
            cdef int32_t nalleles = r.n_allele

            if self.index < 0 or self.index >= nsamples or not r.n_fmt:
                return None

            cdef bcf_fmt_t *fmt0 = r.d.fmt
            cdef int gt0 = is_gt_fmt(hdr, fmt0)

            if not gt0 or not fmt0.n:
                return None

            cdef int32_t  a
            cdef int8_t  *data8
            cdef int16_t *data16
            cdef int32_t *data32
            alleles = []

            if fmt0.type == BCF_BT_INT8:
                data8 = <int8_t *>(fmt0.p + self.index * fmt0.size)
                for i in range(fmt0.n):
                    if data8[i] == bcf_int8_vector_end:
                        break
                    a = (data8[i] >> 1) - 1
                    alleles.append(r.d.allele[a] if 0 <= a < nalleles else None)
            elif fmt0.type == BCF_BT_INT16:
                data16 = <int16_t *>(fmt0.p + self.index * fmt0.size)
                for i in range(fmt0.n):
                    if data16[i] == bcf_int16_vector_end:
                        break
                    a = (data16[i] >> 1) - 1
                    alleles.append(r.d.allele[a] if 0 <= a < nalleles else None)
            elif fmt0.type == BCF_BT_INT32:
                data32 = <int32_t *>(fmt0.p + self.index * fmt0.size)
                for i in range(fmt0.n):
                    if data32[i] == bcf_int32_vector_end:
                        break
                    a = (data32[i] >> 1) - 1
                    alleles.append(r.d.allele[a] if 0 <= a < nalleles else None)

            return tuple(alleles)

    property phased:
        """False if genotype is missing or any allele is unphased.  Otherwise True."""
        def __get__(self):
            cdef bcf_hdr_t *hdr = self.record.header.ptr
            cdef bcf1_t *r = self.record.ptr
            cdef int32_t n = bcf_hdr_nsamples(hdr)

            if self.index < 0 or self.index >= n or not r.n_fmt:
                return False

            cdef bcf_fmt_t *fmt0 = r.d.fmt
            cdef int gt0 = is_gt_fmt(hdr, fmt0)

            if not gt0 or not fmt0.n:
                return False

            cdef int8_t  *data8
            cdef int16_t *data16
            cdef int32_t *data32

            phased = False

            if fmt0.type == BCF_BT_INT8:
                data8 = <int8_t *>(fmt0.p + self.index * fmt0.size)
                for i in range(fmt0.n):
                    if data8[i] == bcf_int8_vector_end:
                        break
                    if i and data8[i] & 1 == 0:
                        return False
                    phased = True
            elif fmt0.type == BCF_BT_INT16:
                data16 = <int16_t *>(fmt0.p + self.index * fmt0.size)
                for i in range(fmt0.n):
                    if data16[i] == bcf_int16_vector_end:
                        break
                    if i and data16[i] & 1 == 0:
                        return False
                    phased = True
            elif fmt0.type == BCF_BT_INT32:
                data32 = <int32_t *>(fmt0.p + self.index * fmt0.size)
                for i in range(fmt0.n):
                    if data32[i] == bcf_int32_vector_end:
                        break
                    if i and data32[i] & 1 == 0:
                        return False
                    phased = True

            return phased

    def __len__(self):
        return self.record.ptr.n_fmt

    def __bool__(self):
        return self.record.ptr.n_fmt != 0

    def __getitem__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_fmt_t *fmt
        cdef int index

        if isinstance(key, int):
            index = key
            if index < 0 or index >= r.n_fmt:
                raise IndexError('invalid format index')
            fmt = r.d.fmt + index
        else:
            fmt = bcf_get_fmt(hdr, r, key)

        if not fmt:
            raise KeyError('invalid format requested')

        if is_gt_fmt(hdr, fmt):
            return self.alleles
        elif fmt.p and fmt.n and fmt.size:
            return bcf_array_to_object(fmt.p + self.index * fmt.size, fmt.type, fmt.n, scalar=1)
        else:
            return None

    def __iter__(self):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef int i, n = r.n_fmt

        for i in range(n):
            yield bcf_hdr_int2id(hdr, BCF_DT_ID, r.d.fmt[i].id)

    def get(self, key, default=None):
        """D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        try:
            return self[key]
        except KeyError:
            return default

    def __contains__(self, key):
        cdef bcf_hdr_t *hdr = self.record.header.ptr
        cdef bcf1_t *r = self.record.ptr
        cdef bcf_fmt_t *fmt = bcf_get_fmt(hdr, r, key)
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

            cregion = region
            with nogil:
                self.iter = bcf_itr_querys(index.ptr, bcf.header.ptr, cregion)
        else:
            if contig is None:
                raise ValueError  # FIXME

            rid = index.refmap.get(contig, -1)

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
    """*(filename, mode=None, header=None, drop_samples=False)*

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
        self.header       = None
        self.index        = None
        self.filename     = None
        self.mode         = None
        self.is_stream    = False
        self.is_remote    = False
        self.is_reading   = False
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

    property category:
        """General file format category.  One of UNKNOWN, ALIGNMENTS, VARIANTS, INDEX, REGIONS"""
        def __get__(self):
            if not self.htsfile:
                raise ValueError('metadata not available on closed file')
            return FORMAT_CATEGORIES[self.htsfile.format.category]

    property format:
        """File format.
           One of UNKNOWN, BINARY_FORMAT, TEXT_FORMAT, SAM, BAM, BAI, CRAM, CRAI, VCF, BCF, CSI, GZI, TBI, BED.
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
            return self.htsfile.format.version.major, self.htsfile.format.version.minor

    property compression:
        """File compression.  One of NONE, GZIP, BGZF, CUSTOM."""
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
                return force_str(desc)
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

        if self.mode[0] != b'r':
            raise ValueError('cannot iterate over Variantfile opened for writing')

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
        cdef char *cfilename, *cmode

        # FIXME: re-open using fd or else header and index could be invalid
        cfilename, cmode = self.filename, self.mode
        with nogil:
            vars.htsfile = hts_open(cfilename, cmode)

        if not vars.htsfile:
            raise ValueError('Cannot re-open htsfile')

        # minimize overhead by re-using header and index.  This approach is
        # currently risky, but see above for how this can be mitigated.
        vars.header       = self.header
        vars.index        = self.index

        vars.filename     = self.filename
        vars.mode         = self.mode
        vars.drop_samples = self.drop_samples
        vars.is_stream    = self.is_stream
        vars.is_remote    = self.is_remote
        vars.is_reading   = self.is_reading
        vars.start_offset = self.start_offset

        if self.htsfile.is_bin:
            vars.seek(self.tell())
        else:
            with nogil:
                hdr = bcf_hdr_read(vars.htsfile)
            makeVariantHeader(hdr)

        return vars

    def open(self, filename, mode=None, VariantHeader header=None, drop_samples=False):
        """open a vcf/bcf file.

        If open is called on an existing VariantFile, the current file will be
        closed and a new file will be opened.
        """
        cdef bcf_hdr_t *hdr
        cdef hts_idx_t *idx
        cdef tbx_t *tidx
        cdef char *cfilename, *cmode

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
        self.is_remote = filename.startswith(b'http:') or filename.startswith(b'ftp:')
        self.is_stream = filename == b'-'

        if mode[0] == b'w':
            # open file for writing

            # header structure (used for writing)
            if header:
                self.header = header.copy()
            else:
                raise ValueError('a VariantHeader must be specified')

            # open file. Header gets written to file at the same time for bam files
            # and sam files (in the latter case, the mode needs to be wh)
            cfilename, cmode = filename, mode
            with nogil:
                self.htsfile = hts_open(cfilename, cmode)

            if not self.htsfile:
                raise ValueError("could not open file `{}` (mode='{}')".format((filename, mode)))

            with nogil:
                bcf_hdr_write(self.htsfile, self.header.ptr)

        elif mode[0] == b'r':
            # open file for reading
            if filename != b'-' and not self.is_remote and not os.path.exists(filename):
                raise IOError('file `{}` not found'.format(filename))

            cfilename, cmode = filename, mode
            with nogil:
                self.htsfile = hts_open(cfilename, cmode)

            if not self.htsfile:
                raise ValueError("could not open file `{}` (mode='{}') - is it VCF/BCF format?".format((filename, mode)))

            with nogil:
                hdr = bcf_hdr_read(self.htsfile)
            self.header = makeVariantHeader(hdr)

            if not self.header:
                raise ValueError("file `{}` does not have valid header (mode='{}') - is it BCF format?".format((filename, mode)))

            # check for index and open if present
            if self.htsfile.format.format == bcf:
                cfilename = filename
                with nogil:
                    idx = bcf_index_load(cfilename)
                self.index = makeBCFIndex(self.header, idx)
            else:
                tabix_filename = filename + '.tbi'
                cfilename = tabix_filename
                with nogil:
                    tidx = tbx_index_load(cfilename)
                self.index = makeTabixIndex(tidx)

            if not self.is_stream:
                self.start_offset = self.tell()

    def reset(self):
        """reset file position to beginning of file just after the header."""
        return self.seek(self.start_offset, 0)

    def seek(self, uint64_t offset):
        """move file pointer to position *offset*, see :meth:`pysam.VariantFile.tell`."""
        if not self.is_open:
            raise ValueError('I/O operation on closed file')
        if self.is_stream:
            raise OSError('seek not available in streams')

        cdef int ret
        if self.htsfile.format.compression != no_compression:
            with nogil:
                ret = bgzf_seek(hts_get_bgzfp(self.htsfile), offset, SEEK_SET)
        else:
            with nogil:
                ret = hts_useek(self.htsfile, offset, SEEK_SET)
        return ret


    def tell(self):
        """return current file position, see :meth:`pysam.VariantFile.seek`."""
        if not self.is_open:
            raise ValueError('I/O operation on closed file')
        if self.is_stream:
            raise OSError('tell not available in streams')

        cdef int ret
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

        Note that a :term:`VCF` file without a tabix index (.tbi) or a
        :term:`BCF` file without a CSI index can only be read sequentially.
        """
        if not self.is_open:
            raise ValueError('I/O operation on closed file')

        if self.mode[0] != b'r':
            raise ValueError('cannot fetch from Variantfile opened for writing')

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
            return 0

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

        if self.mode[0] != b'r':
            raise ValueError('cannot subset samples from Variantfile opened for writing')

        if self.is_reading:
            raise ValueError('cannot subset samples after fetching records')

        self.header._subset_samples(include_samples)

        # potentially unnecessary optimization that also sets max_unpack
        if not include_samples:
            self.drop_samples = True
