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

from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup

from chtslib     cimport *


cdef class VariantHeader(object):
    cdef bcf_hdr_t *ptr

    cdef _subset_samples(self, include_samples)


cdef class VariantHeaderRecord(object):
    cdef VariantHeader header
    cdef bcf_hrec_t *ptr


cdef class VariantHeaderRecords(object):
    cdef VariantHeader header


cdef class VariantHeaderContigs(object):
    cdef VariantHeader header


cdef class VariantHeaderSamples(object):
    cdef VariantHeader header


cdef class VariantContig(object):
    cdef VariantHeader header
    cdef int id


cdef class VariantMetadata(object):
    cdef VariantHeader header
    cdef int type
    cdef int id


cdef class VariantHeaderMetadata(object):
    cdef VariantHeader header
    cdef int32_t type


cdef class VariantRecord(object):
    cdef VariantHeader header
    cdef bcf1_t *ptr


cdef class VariantRecordFilter(object):
    cdef VariantRecord record


cdef class VariantRecordFormat(object):
    cdef VariantRecord record


cdef class VariantRecordInfo(object):
    cdef VariantRecord record


cdef class VariantRecordSamples(object):
    cdef VariantRecord record


cdef class VariantRecordSample(object):
    cdef VariantRecord record
    cdef readonly int32_t index


cdef class BaseIndex(object):
    cdef tuple refs
    cdef dict refmap


cdef class BCFIndex(BaseIndex):
    cdef VariantHeader header
    cdef hts_idx_t *ptr


cdef class TabixIndex(BaseIndex):
    cdef tbx_t *ptr


cdef class BaseIterator(object):
    cdef VariantFile bcf
    cdef hts_itr_t  *iter


cdef class BCFIterator(BaseIterator):
    cdef BCFIndex index


cdef class TabixIterator(BaseIterator):
    cdef TabixIndex index
    cdef kstring_t line_buffer


cdef class VariantFile(object):
    cdef htsFile *htsfile                  # pointer to htsFile structure
    cdef int64_t  start_offset             # BGZF offset of first record

    cdef readonly object     filename      # filename as supplied by user
    cdef readonly object     mode          # file opening mode

    cdef readonly VariantHeader  header
    cdef readonly BaseIndex  index

    cdef readonly bint       drop_samples  # true if sample information is to be ignored

    # FIXME: Temporary, use htsFormat when it is available
    cdef readonly bint       is_bcf        # true if file is a bcf file
    cdef readonly bint       is_stream     # true if not a seekable file but a stream
    cdef readonly bint       is_remote     # true if file is not on the local filesystem
    cdef readonly bint       is_reading    # true if file has begun reading records

    cpdef int write(self, VariantRecord record) except -1
