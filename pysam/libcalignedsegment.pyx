# cython: language_level=3
# cython: embedsignature=True
# cython: profile=True
###############################################################################
###############################################################################
# Cython wrapper for SAM/BAM/CRAM files based on htslib
###############################################################################
# The principal classes defined in this module are:
#
# class AlignedSegment  an aligned segment (read)
#
# class PileupColumn    a collection of segments (PileupRead) aligned to
#                       a particular genomic position.
#
# class PileupRead      an AlignedSegment aligned to a particular genomic
#                       position. Contains additional attributes with respect
#                       to this.
#
# Additionally this module defines numerous additional classes that are part
# of the internal API. These are:
#
# Various iterator classes to iterate over alignments in sequential (IteratorRow)
# or in a stacked fashion (IteratorColumn):
#
# class IteratorRow
# class IteratorRowRegion
# class IteratorRowHead
# class IteratorRowAll
# class IteratorRowAllRefs
# class IteratorRowSelection
#
###############################################################################
#
# The MIT License
#
# Copyright (c) 2015 Andreas Heger
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
import re
import array
import json
import string
import ctypes
import struct

cimport cython
from cpython cimport array as c_array
from cpython cimport PyBytes_FromStringAndSize
from libc.string cimport memset, strchr
from cpython cimport array as c_array
from libc.stdint cimport INT8_MIN, INT16_MIN, INT32_MIN, \
    INT8_MAX, INT16_MAX, INT32_MAX, \
    UINT8_MAX, UINT16_MAX, UINT32_MAX

from pysam.libchtslib cimport HTS_IDX_NOCOOR
from pysam.libcutils cimport force_bytes, force_str, \
    charptr_to_str, charptr_to_bytes
from pysam.libcutils cimport qualities_to_qualitystring, qualitystring_to_array, \
    array_to_qualitystring

# Constants for binary tag conversion
cdef char * htslib_types = 'cCsSiIf'
cdef char * parray_types = 'bBhHiIf'

# translation tables

# cigar code to character and vice versa
cdef char* CODE2CIGAR= "MIDNSHP=XB"
cdef int NCIGAR_CODES = 10

CIGAR2CODE = dict([y, x] for x, y in enumerate(CODE2CIGAR))
CIGAR_REGEX = re.compile("(\d+)([MIDNSHP=XB])")

# names for keys in dictionary representation of an AlignedSegment
KEY_NAMES = ["name", "flag", "ref_name", "ref_pos", "map_quality", "cigar",
             "next_ref_name", "next_ref_pos", "length", "seq", "qual", "tags"]

#####################################################################
# C multiplication with wrapping around
cdef inline uint32_t c_mul(uint32_t a, uint32_t b):
    return (a * b) & 0xffffffff


cdef inline uint8_t tolower(uint8_t ch):
    if ch >= 65 and ch <= 90:
        return ch + 32
    else:
        return ch


cdef inline uint8_t toupper(uint8_t ch):
    if ch >= 97 and ch <= 122:
        return ch - 32
    else:
        return ch


cdef inline uint8_t strand_mark_char(uint8_t ch, bam1_t *b):
    if ch == b'=':
        if bam_is_rev(b):
            return b','
        else:
            return b'.'
    else:
        if bam_is_rev(b):
            return tolower(ch)
        else:
            return toupper(ch)


cdef inline bint pileup_base_qual_skip(const bam_pileup1_t * p, uint32_t threshold):
    cdef uint32_t c
    if p.qpos < p.b.core.l_qseq:
        c = bam_get_qual(p.b)[p.qpos]
    else:
        c = 0
    if c < threshold:
        return True
    return False


cdef inline char map_typecode_htslib_to_python(uint8_t s):
    """map an htslib typecode to the corresponding python typecode
    to be used in the struct or array modules."""

    # map type from htslib to python array
    cdef char * f = strchr(htslib_types, s)

    if f == NULL:
        return 0
    return parray_types[f - htslib_types]


cdef inline uint8_t map_typecode_python_to_htslib(char s):
    """determine value type from type code of array"""
    cdef char * f = strchr(parray_types, s)
    if f == NULL:
        return 0
    return htslib_types[f - parray_types]


cdef inline void update_bin(bam1_t * src):
    if src.core.flag & BAM_FUNMAP:
        # treat alignment as length of 1 for unmapped reads
        src.core.bin = hts_reg2bin(
            src.core.pos,
            src.core.pos + 1,
            14,
            5)
    elif pysam_get_n_cigar(src):
        src.core.bin = hts_reg2bin(
            src.core.pos,
            bam_endpos(src),
            14,
            5)
    else:
        src.core.bin = hts_reg2bin(
            src.core.pos,
            src.core.pos + 1,
            14,
            5)


# optional tag data manipulation
cdef convert_binary_tag(uint8_t * tag):
    """return bytesize, number of values and array of values
    in aux_data memory location pointed to by tag."""
    cdef uint8_t auxtype
    cdef uint8_t byte_size
    cdef int32_t nvalues
    # get byte size
    auxtype = tag[0]
    byte_size = aux_type2size(auxtype)
    tag += 1
    # get number of values in array
    nvalues = (<int32_t*>tag)[0]
    tag += 4

    # define python array
    cdef c_array.array c_values = array.array(
        chr(map_typecode_htslib_to_python(auxtype)))
    c_array.resize(c_values, nvalues)

    # copy data
    memcpy(c_values.data.as_voidptr, <uint8_t*>tag, nvalues * byte_size)

    # no need to check for endian-ness as bam1_core_t fields
    # and aux_data are in host endian-ness. See sam.c and calls
    # to swap_data
    return byte_size, nvalues, c_values


cdef inline uint8_t get_tag_typecode(value, value_type=None):
    """guess type code for a *value*. If *value_type* is None, the type
    code will be inferred based on the Python type of *value*

    """
    # 0 is unknown typecode
    cdef char typecode = 0

    if value_type is None:
        if isinstance(value, int):
            if value < 0:
                if value >= INT8_MIN:
                    typecode = b'c'
                elif value >= INT16_MIN:
                    typecode = b's'
                elif value >= INT32_MIN:
                    typecode = b'i'
            # unsigned ints
            else:
                if value <= UINT8_MAX:
                    typecode = b'C'
                elif value <= UINT16_MAX:
                    typecode = b'S'
                elif value <= UINT32_MAX:
                    typecode = b'I'
        elif isinstance(value, float):
            typecode = b'f'
        elif isinstance(value, str):
            typecode = b'Z'
        elif isinstance(value, bytes):
            typecode = b'Z'
        elif isinstance(value, array.array) or \
                isinstance(value, list) or \
                isinstance(value, tuple):
            typecode = b'B'
    else:
        if value_type in 'aAsSIcCZidfH':
            typecode = force_bytes(value_type)[0]

    return typecode


cdef inline uint8_t get_btag_typecode(value, min_value=None, max_value=None):
    '''returns the value typecode of a value.

    If max is specified, the appropriate type is returned for a range
    where value is the minimum.

    Note that this method returns types from the extended BAM alphabet
    of types that includes tags that are not part of the SAM
    specification.
    '''


    cdef uint8_t typecode

    t = type(value)

    if t is float:
        typecode = b'f'
    elif t is int:
        if max_value is None:
            max_value = value
        if min_value is None:
            min_value = value
        # signed ints
        if min_value < 0:
            if min_value >= INT8_MIN and max_value <= INT8_MAX:
                typecode = b'c'
            elif min_value >= INT16_MIN and max_value <= INT16_MAX:
                typecode = b's'
            elif min_value >= INT32_MIN or max_value <= INT32_MAX:
                typecode = b'i'
            else:
                raise ValueError(
                    "at least one signed integer out of range of "
                    "BAM/SAM specification")
        # unsigned ints
        else:
            if max_value <= UINT8_MAX:
                typecode = b'C'
            elif max_value <= UINT16_MAX:
                typecode = b'S'
            elif max_value <= UINT32_MAX:
                typecode = b'I'
            else:
                raise ValueError(
                    "at least one integer out of range of BAM/SAM specification")
    else:
        # Note: hex strings (H) are not supported yet
        if t is not bytes:
            value = value.encode('ascii')
        if len(value) == 1:
            typecode = b'A'
        else:
            typecode = b'Z'

    return typecode


# mapping python array.array and htslib typecodes to struct typecodes
DATATYPE2FORMAT = {
    ord('c'): ('b', 1),
    ord('C'): ('B', 1),
    ord('s'): ('h', 2),
    ord('S'): ('H', 2),
    ord('i'): ('i', 4),
    ord('I'): ('I', 4),
    ord('f'): ('f', 4),
    ord('d'): ('d', 8),
    ord('A'): ('c', 1),
    ord('a'): ('c', 1)}


cdef inline pack_tags(tags):
    """pack a list of tags. Each tag is a tuple of (tag, tuple).

    Values are packed into the most space efficient data structure
    possible unless the tag contains a third field with the typecode.

    Returns a format string and the associated list of arguments to be
    used in a call to struct.pack_into.
    """
    fmts, args = ["<"], []

    # htslib typecode
    cdef uint8_t typecode
    for tag in tags:

        if len(tag) == 2:
            pytag, value = tag
            valuetype = None
        elif len(tag) == 3:
            pytag, value, valuetype = tag
        else:
            raise ValueError("malformatted tag: %s" % str(tag))

        if valuetype is None:
            typecode = 0
        else:
            # only first character in valuecode matters
            typecode = force_bytes(valuetype)[0]

        pytag = force_bytes(pytag)
        pytype = type(value)

        if pytype is tuple or pytype is list:
            # binary tags from tuples or lists
            if not typecode:
                # automatically determine value type - first value
                # determines type. If there is a mix of types, the
                # result is undefined.
                typecode = get_btag_typecode(min(value),
                                             min_value=min(value),
                                             max_value=max(value))

            if typecode not in DATATYPE2FORMAT:
                raise ValueError("invalid value type '{}'".format(chr(typecode)))

            datafmt = "2sBBI%i%s" % (len(value), DATATYPE2FORMAT[typecode][0])
            args.extend([pytag[:2],
                         ord("B"),
                         typecode,
                         len(value)] + list(value))

        elif isinstance(value, array.array):
            # binary tags from arrays
            if typecode == 0:
                typecode = map_typecode_python_to_htslib(ord(value.typecode))

                if typecode == 0:
                    raise ValueError("unsupported type code '{}'".format(value.typecode))

            if typecode not in DATATYPE2FORMAT:
                raise ValueError("invalid value type '{}' ({})".format(chr(typecode), array.typecode))

            # use array.tostring() to retrieve byte representation and
            # save as bytes
            datafmt = "2sBBI%is" % (len(value) * DATATYPE2FORMAT[typecode][1])
            args.extend([pytag[:2],
                         ord("B"),
                         typecode,
                         len(value),
                         value.tobytes()])

        else:
            if typecode == 0:
                typecode = get_tag_typecode(value)
                if typecode == 0:
                    raise ValueError("could not deduce typecode for value {}".format(value))

            if typecode == b'a' or typecode == b'A' or typecode == b'Z' or typecode == b'H':
                value = force_bytes(value)

            if typecode == b"a":
                typecode = b'A'

            if typecode == b'Z' or typecode == b'H':
                datafmt = "2sB%is" % (len(value)+1)
            else:
                datafmt = "2sB%s" % DATATYPE2FORMAT[typecode][0]

            args.extend([pytag[:2],
                         typecode,
                         value])

        fmts.append(datafmt)

    return "".join(fmts), args


cdef inline int32_t calculateQueryLengthWithoutHardClipping(bam1_t * src):
    """return query length computed from CIGAR alignment.

    Length ignores hard-clipped bases.

    Return 0 if there is no CIGAR alignment.
    """

    cdef uint32_t * cigar_p = pysam_bam_get_cigar(src)

    if cigar_p == NULL:
        return 0

    cdef uint32_t k, qpos
    cdef int op
    qpos = 0

    for k from 0 <= k < pysam_get_n_cigar(src):
        op = cigar_p[k] & BAM_CIGAR_MASK

        if op == BAM_CMATCH or \
           op == BAM_CINS or \
           op == BAM_CSOFT_CLIP or \
           op == BAM_CEQUAL or \
           op == BAM_CDIFF:
            qpos += cigar_p[k] >> BAM_CIGAR_SHIFT

    return qpos


cdef inline int32_t calculateQueryLengthWithHardClipping(bam1_t * src):
    """return query length computed from CIGAR alignment.

    Length includes hard-clipped bases.

    Return 0 if there is no CIGAR alignment.
    """

    cdef uint32_t * cigar_p = pysam_bam_get_cigar(src)

    if cigar_p == NULL:
        return 0

    cdef uint32_t k, qpos
    cdef int op
    qpos = 0

    for k from 0 <= k < pysam_get_n_cigar(src):
        op = cigar_p[k] & BAM_CIGAR_MASK

        if op == BAM_CMATCH or \
           op == BAM_CINS or \
           op == BAM_CSOFT_CLIP or \
           op == BAM_CHARD_CLIP or \
           op == BAM_CEQUAL or \
           op == BAM_CDIFF:
            qpos += cigar_p[k] >> BAM_CIGAR_SHIFT

    return qpos


cdef inline int32_t getQueryStart(bam1_t *src) except -1:
    cdef uint32_t * cigar_p
    cdef uint32_t start_offset = 0
    cdef uint32_t k, op

    cigar_p = pysam_bam_get_cigar(src);
    for k from 0 <= k < pysam_get_n_cigar(src):
        op = cigar_p[k] & BAM_CIGAR_MASK
        if op == BAM_CHARD_CLIP:
            if start_offset != 0 and start_offset != src.core.l_qseq:
                raise ValueError('Invalid clipping in CIGAR string')
        elif op == BAM_CSOFT_CLIP:
            start_offset += cigar_p[k] >> BAM_CIGAR_SHIFT
        else:
            break

    return start_offset


cdef inline int32_t getQueryEnd(bam1_t *src) except -1:
    cdef uint32_t * cigar_p = pysam_bam_get_cigar(src)
    cdef uint32_t end_offset = src.core.l_qseq
    cdef uint32_t k, op

    # if there is no sequence, compute length from cigar string
    if end_offset == 0:
        for k from 0 <= k < pysam_get_n_cigar(src):
            op = cigar_p[k] & BAM_CIGAR_MASK
            if op == BAM_CMATCH or \
               op == BAM_CINS or \
               op == BAM_CEQUAL or \
               op == BAM_CDIFF or \
              (op == BAM_CSOFT_CLIP and end_offset == 0):
                end_offset += cigar_p[k] >> BAM_CIGAR_SHIFT
    else:
        # walk backwards in cigar string
        for k from pysam_get_n_cigar(src) > k >= 1:
            op = cigar_p[k] & BAM_CIGAR_MASK
            if op == BAM_CHARD_CLIP:
                if end_offset != src.core.l_qseq:
                    raise ValueError('Invalid clipping in CIGAR string')
            elif op == BAM_CSOFT_CLIP:
                end_offset -= cigar_p[k] >> BAM_CIGAR_SHIFT
            else:
                break

    return end_offset


cdef inline bytes getSequenceInRange(bam1_t *src,
                                     uint32_t start,
                                     uint32_t end):
    """return python string of the sequence in a bam1_t object.
    """

    cdef uint8_t * p
    cdef uint32_t k
    cdef char * s

    if not src.core.l_qseq:
        return None

    seq = PyBytes_FromStringAndSize(NULL, end - start)
    s   = <char*>seq
    p   = pysam_bam_get_seq(src)

    for k from start <= k < end:
        # equivalent to seq_nt16_str[bam1_seqi(s, i)] (see bam.c)
        # note: do not use string literal as it will be a python string
        s[k-start] = seq_nt16_str[p[k//2] >> 4 * (1 - k%2) & 0xf]

    return charptr_to_bytes(seq)


cdef inline object getQualitiesInRange(bam1_t *src,
                                       uint32_t start,
                                       uint32_t end):
    """return python array of quality values from a bam1_t object"""

    cdef uint8_t * p
    cdef uint32_t k

    p = pysam_bam_get_qual(src)
    if p[0] == 0xff:
        return None

    # 'B': unsigned char
    cdef c_array.array result = array.array('B', [0])
    c_array.resize(result, end - start)

    # copy data
    memcpy(result.data.as_voidptr, <void*>&p[start], end - start)

    return result


#####################################################################
## factory methods for instantiating extension classes
cdef class AlignedSegment
cdef AlignedSegment makeAlignedSegment(bam1_t *src,
                                       AlignmentHeader header):
    '''return an AlignedSegment object constructed from `src`'''
    # note that the following does not call __init__
    cdef AlignedSegment dest = AlignedSegment.__new__(AlignedSegment)
    dest._delegate = bam_dup1(src)
    dest.header = header
    return dest


cdef class PileupColumn
cdef PileupColumn makePileupColumn(const bam_pileup1_t ** plp,
                      int tid,
                      int pos,
                      int n_pu,
                      uint32_t min_base_quality,
                      char * reference_sequence,
                      AlignmentHeader header):
    '''return a PileupColumn object constructed from pileup in `plp` and
    setting additional attributes.

    '''
    # note that the following does not call __init__
    cdef PileupColumn dest = PileupColumn.__new__(PileupColumn)
    dest.header = header
    dest.plp = plp
    dest.tid = tid
    dest.pos = pos
    dest.n_pu = n_pu
    dest.min_base_quality = min_base_quality
    dest.reference_sequence = reference_sequence
    dest.buf.l = dest.buf.m = 0
    dest.buf.s = NULL

    return dest


cdef class PileupRead
cdef PileupRead makePileupRead(const bam_pileup1_t *src,
                               AlignmentHeader header):
    '''return a PileupRead object construted from a bam_pileup1_t * object.'''
    # note that the following does not call __init__
    cdef PileupRead dest = PileupRead.__new__(PileupRead)
    dest._alignment = makeAlignedSegment(src.b, header)
    dest._qpos = src.qpos
    dest._indel = src.indel
    dest._level = src.level
    dest._is_del = src.is_del
    dest._is_head = src.is_head
    dest._is_tail = src.is_tail
    dest._is_refskip = src.is_refskip
    return dest


cdef inline uint32_t get_alignment_length(bam1_t *src):
    cdef uint32_t k = 0
    cdef uint32_t l = 0
    if src == NULL:
        return 0
    cdef uint32_t * cigar_p = bam_get_cigar(src)
    if cigar_p == NULL:
        return 0
    cdef int op
    cdef uint32_t n = pysam_get_n_cigar(src)
    for k from 0 <= k < n:
        op = cigar_p[k] & BAM_CIGAR_MASK
        if op == BAM_CSOFT_CLIP or op == BAM_CHARD_CLIP:
            continue
        l += cigar_p[k] >> BAM_CIGAR_SHIFT
    return l


cdef inline uint32_t get_md_reference_length(char * md_tag):
    cdef int l = 0
    cdef int md_idx = 0
    cdef int nmatches = 0

    while md_tag[md_idx] != 0:
        if md_tag[md_idx] >= 48 and md_tag[md_idx] <= 57:
            nmatches *= 10
            nmatches += md_tag[md_idx] - 48
            md_idx += 1
            continue
        else:
            l += nmatches
            nmatches = 0
            if md_tag[md_idx] == b'^':
                md_idx += 1
                while md_tag[md_idx] >= 65 and md_tag[md_idx] <= 90:
                    md_idx += 1
                    l += 1
            else:
                md_idx += 1
                l += 1

    l += nmatches
    return l

# TODO: avoid string copying for getSequenceInRange, reconstituneSequenceFromMD, ...
cdef inline bytes build_alignment_sequence(bam1_t * src):
    """return expanded sequence from MD tag.

    The sequence includes substitutions and both insertions in the
    reference as well as deletions to the reference sequence. Combine
    with the cigar string to reconstitute the query or the reference
    sequence.

    Positions corresponding to `N` (skipped region from the reference)
    or `P` (padding (silent deletion from padded reference)) in the CIGAR
    string will not appear in the returned sequence. The MD should
    correspondingly not contain these. Thus proper tags are::

       Deletion from the reference:    cigar=5M1D5M    MD=5^C5
       Skipped region from reference:  cigar=5M1N5M    MD=10
       Padded region in the reference: cigar=5M1P5M    MD=10

    Returns
    -------

    None, if no MD tag is present.

    """
    if src == NULL:
        return None

    cdef uint8_t * md_tag_ptr = bam_aux_get(src, "MD")
    if md_tag_ptr == NULL:
        return None

    cdef uint32_t start = getQueryStart(src)
    cdef uint32_t end = getQueryEnd(src)
    # get read sequence, taking into account soft-clipping
    r = getSequenceInRange(src, start, end)
    cdef char * read_sequence = r
    cdef uint32_t * cigar_p = pysam_bam_get_cigar(src)
    if cigar_p == NULL:
        return None

    cdef uint32_t r_idx = 0
    cdef int op
    cdef uint32_t k, i, l, x
    cdef int nmatches = 0
    cdef int s_idx = 0

    cdef uint32_t max_len = get_alignment_length(src)
    if max_len == 0:
        raise ValueError("could not determine alignment length")

    cdef char * s = <char*>calloc(max_len + 1, sizeof(char))
    if s == NULL:
        raise ValueError(
            "could not allocate sequence of length %i" % max_len)

    for k from 0 <= k < pysam_get_n_cigar(src):
        op = cigar_p[k] & BAM_CIGAR_MASK
        l = cigar_p[k] >> BAM_CIGAR_SHIFT
        if op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
            for i from 0 <= i < l:
                s[s_idx] = read_sequence[r_idx]
                r_idx += 1
                s_idx += 1
        elif op == BAM_CDEL:
            for i from 0 <= i < l:
                s[s_idx] = b'-'
                s_idx += 1
        elif op == BAM_CREF_SKIP:
            pass
        elif op == BAM_CINS or op == BAM_CPAD:
            for i from 0 <= i < l:
                # encode insertions into reference as lowercase
                s[s_idx] = read_sequence[r_idx] + 32
                r_idx += 1
                s_idx += 1
        elif op == BAM_CSOFT_CLIP:
            pass
        elif op == BAM_CHARD_CLIP:
            pass # advances neither

    cdef char * md_tag = <char*>bam_aux2Z(md_tag_ptr)
    cdef int md_idx = 0
    cdef char c
    s_idx = 0

    # Check if MD tag is valid by matching CIGAR length to MD tag defined length
    # Insertions would be in addition to what is described by MD, so we calculate
    # the number of insertions separately.
    cdef int insertions = 0

    while s[s_idx] != 0:
        if s[s_idx] >= b'a':
            insertions += 1
        s_idx += 1
    s_idx = 0

    cdef uint32_t md_len = get_md_reference_length(md_tag)
    if md_len + insertions > max_len:
        free(s)
        raise AssertionError(
            "Invalid MD tag: MD length {} mismatch with CIGAR length {} and {} insertions".format(
            md_len, max_len, insertions))

    while md_tag[md_idx] != 0:
        # c is numerical
        if md_tag[md_idx] >= 48 and md_tag[md_idx] <= 57:
            nmatches *= 10
            nmatches += md_tag[md_idx] - 48
            md_idx += 1
            continue
        else:
            # save matches up to this point, skipping insertions
            for x from 0 <= x < nmatches:
                while s[s_idx] >= b'a':
                    s_idx += 1
                s_idx += 1
            while s[s_idx] >= b'a':
                s_idx += 1

            r_idx += nmatches
            nmatches = 0
            if md_tag[md_idx] == b'^':
                md_idx += 1
                while md_tag[md_idx] >= 65 and md_tag[md_idx] <= 90:
                    # assert s[s_idx] == '-'
                    s[s_idx] = md_tag[md_idx]
                    s_idx += 1
                    md_idx += 1
            else:
                # save mismatch
                # enforce lower case
                c = md_tag[md_idx]
                if c <= 90:
                    c += 32
                s[s_idx] = c
                s_idx += 1
                r_idx += 1
                md_idx += 1

    # save matches up to this point, skipping insertions
    for x from 0 <= x < nmatches:
        while s[s_idx] >= b'a':
            s_idx += 1
        s_idx += 1
    while s[s_idx] >= b'a':
        s_idx += 1

    seq = PyBytes_FromStringAndSize(s, s_idx)
    free(s)

    return seq


cdef inline bytes build_reference_sequence(bam1_t * src):
    """return the reference sequence in the region that is covered by the
    alignment of the read to the reference.

    This method requires the MD tag to be set.

    """
    cdef uint32_t k, i, l
    cdef int op
    cdef int s_idx = 0
    ref_seq = build_alignment_sequence(src)
    if ref_seq is None:
        raise ValueError("MD tag not present")

    cdef char * s = <char*>calloc(len(ref_seq) + 1, sizeof(char))
    if s == NULL:
        raise ValueError(
            "could not allocate sequence of length %i" % len(ref_seq))

    cdef char * cref_seq = ref_seq
    cdef uint32_t * cigar_p = pysam_bam_get_cigar(src)
    cdef uint32_t r_idx = 0
    for k from 0 <= k < pysam_get_n_cigar(src):
        op = cigar_p[k] & BAM_CIGAR_MASK
        l = cigar_p[k] >> BAM_CIGAR_SHIFT
        if op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
            for i from 0 <= i < l:
                s[s_idx] = cref_seq[r_idx]
                r_idx += 1
                s_idx += 1
        elif op == BAM_CDEL:
            for i from 0 <= i < l:
                s[s_idx] = cref_seq[r_idx]
                r_idx += 1
                s_idx += 1
        elif op == BAM_CREF_SKIP:
            pass
        elif op == BAM_CINS or op == BAM_CPAD:
            r_idx += l
        elif op == BAM_CSOFT_CLIP:
            pass
        elif op == BAM_CHARD_CLIP:
            pass # advances neither

    seq = PyBytes_FromStringAndSize(s, s_idx)
    free(s)

    return seq


cdef class AlignedSegment:
    '''Class representing an aligned segment.

    This class stores a handle to the samtools C-structure representing
    an aligned read. Member read access is forwarded to the C-structure
    and converted into python objects. This implementation should be fast,
    as only the data needed is converted.

    For write access, the C-structure is updated in-place. This is
    not the most efficient way to build BAM entries, as the variable
    length data is concatenated and thus needs to be resized if
    a field is updated. Furthermore, the BAM entry might be
    in an inconsistent state.

    One issue to look out for is that the sequence should always
    be set *before* the quality scores. Setting the sequence will
    also erase any quality scores that were set previously.

    Parameters
    ----------

    header:
         :class:`~pysam.AlignmentHeader` object to map numerical
         identifiers to chromosome names. If not given, an empty
         header is created.
    '''

    # Now only called when instances are created from Python
    def __init__(self, AlignmentHeader header=None):
        # see bam_init1
        self._delegate = <bam1_t*>calloc(1, sizeof(bam1_t))
        if self._delegate == NULL:
            raise MemoryError("could not allocated memory of {} bytes".format(sizeof(bam1_t)))
        # allocate some memory. If size is 0, calloc does not return a
        # pointer that can be passed to free() so allocate 40 bytes
        # for a new read
        self._delegate.m_data = 40
        self._delegate.data = <uint8_t *>calloc(
            self._delegate.m_data, 1)
        if self._delegate.data == NULL:
            raise MemoryError("could not allocate memory of {} bytes".format(self._delegate.m_data))
        self._delegate.l_data = 0
        # set some data to make read approximately legit.
        # Note, SAM writing fails with q_name of length 0
        self._delegate.core.l_qname = 0
        self._delegate.core.tid = -1
        self._delegate.core.pos = -1
        self._delegate.core.mtid = -1
        self._delegate.core.mpos = -1

        # caching for selected fields
        self.cache_query_qualities = None
        self.cache_query_alignment_qualities = None
        self.cache_query_sequence = None
        self.cache_query_alignment_sequence = None

        self.header = header

    def __dealloc__(self):
        bam_destroy1(self._delegate)

    def __str__(self):
        """return string representation of alignment.

        The representation is an approximate :term:`SAM` format, because
        an aligned read might not be associated with a :term:`AlignmentFile`.
        As a result :term:`tid` is shown instead of the reference name.
        Similarly, the tags field is returned in its parsed state.

        To get a valid SAM record, use :meth:`to_string`.
        """
        # sam-parsing is done in sam.c/bam_format1_core which
        # requires a valid header.
        return "\t".join(map(str, (self.query_name,
                                   self.flag,
                                   "#%d" % self.reference_id if self.reference_id >= 0 else "*",
                                   self.reference_start + 1,
                                   self.mapping_quality,
                                   self.cigarstring,
                                   "#%d" % self.next_reference_id if self.next_reference_id >= 0 else "*",
                                   self.next_reference_start + 1,
                                   self.template_length,
                                   self.query_sequence,
                                   self.query_qualities,
                                   self.tags)))

    def __copy__(self):
        return makeAlignedSegment(self._delegate, self.header)

    def __deepcopy__(self, memo):
        return makeAlignedSegment(self._delegate, self.header)

    def compare(self, AlignedSegment other):
        '''return -1,0,1, if contents in this are binary
        <,=,> to *other*
        '''

        # avoid segfault when other equals None
        if other is None:
            return -1

        cdef int retval, x
        cdef bam1_t *t
        cdef bam1_t *o

        t = self._delegate
        o = other._delegate

        # uncomment for debugging purposes
        # cdef unsigned char * oo, * tt
        # tt = <unsigned char*>(&t.core)
        # oo = <unsigned char*>(&o.core)
        # for x from 0 <= x < sizeof( bam1_core_t): print x, tt[x], oo[x]
        # tt = <unsigned char*>(t.data)
        # oo = <unsigned char*>(o.data)
        # for x from 0 <= x < max(t.l_data, o.l_data): print x, tt[x], oo[x], chr(tt[x]), chr(oo[x])

        # Fast-path test for object identity
        if t == o:
            return 0

        cdef uint8_t *a = <uint8_t*>&t.core
        cdef uint8_t *b = <uint8_t*>&o.core

        retval = memcmp(&t.core, &o.core, sizeof(bam1_core_t))
        if retval:
            return retval

        # cmp(t.l_data, o.l_data)
        retval = (t.l_data > o.l_data) - (t.l_data < o.l_data)
        if retval:
            return retval
        return memcmp(t.data, o.data, t.l_data)

    def __richcmp__(self, AlignedSegment other, int op):
        if op == 2:  # == operator
            return self.compare(other) == 0
        elif op == 3:  # != operator
            return self.compare(other) != 0
        else:
            return NotImplemented

    def __hash__(self):
        cdef bam1_t * src = self._delegate
        cdef int x

        # see http://effbot.org/zone/python-hash.htm
        cdef uint8_t * c = <uint8_t *>&src.core
        cdef uint32_t hash_value = c[0]
        for x from 1 <= x < sizeof(bam1_core_t):
            hash_value = c_mul(hash_value, 1000003) ^ c[x]
        c = <uint8_t *>src.data
        for x from 0 <= x < src.l_data:
            hash_value = c_mul(hash_value, 1000003) ^ c[x]

        return hash_value

    cpdef to_string(self):
        """returns a string representation of the aligned segment.

        The output format is valid SAM format if a header is associated
        with the AlignedSegment.
        """
        cdef kstring_t line
        line.l = line.m = 0
        line.s = NULL

        if self.header:
            if sam_format1(self.header.ptr, self._delegate, &line) < 0:
                if line.m:
                    free(line.s)
                raise ValueError('sam_format failed')
        else:
            raise NotImplementedError("todo")

        ret = force_str(line.s[:line.l])

        if line.m:
            free(line.s)

        return ret

    @classmethod
    def fromstring(cls, sam, AlignmentHeader header):
        """parses a string representation of the aligned segment.

        The input format should be valid SAM format.

        Parameters
        ----------
        sam:
            :term:`SAM` formatted string

        """
        cdef AlignedSegment dest = cls.__new__(cls)
        dest._delegate = <bam1_t*>calloc(1, sizeof(bam1_t))
        dest.header = header

        cdef kstring_t line
        line.l = line.m = len(sam)
        _sam = force_bytes(sam)
        line.s = _sam

        sam_parse1(&line, dest.header.ptr, dest._delegate)

        return dest

    cpdef tostring(self, htsfile=None):
        """deprecated, use :meth:`to_string()` instead.

        Parameters
        ----------

        htsfile:
            (deprecated) AlignmentFile object to map numerical
            identifiers to chromosome names. This parameter is present
            for backwards compatibility and ignored.
        """

        return self.to_string()

    def to_dict(self):
        """returns a json representation of the aligned segment.

        Field names are abbreviated versions of the class attributes.
        """
        # let htslib do the string conversions, but treat optional field properly as list
        vals = self.to_string().split("\t")
        n = len(KEY_NAMES) - 1
        return dict(list(zip(KEY_NAMES[:-1], vals[:n])) + [(KEY_NAMES[-1], vals[n:])])

    @classmethod
    def from_dict(cls, sam_dict, AlignmentHeader header):
        """parses a dictionary representation of the aligned segment.

        Parameters
        ----------
        sam_dict:
            dictionary of alignment values, keys corresponding to output from
            :meth:`todict()`.

        """
        # let htslib do the parsing
        # the tags field can be missing
        return cls.fromstring(
            "\t".join((sam_dict[x] for x in KEY_NAMES[:-1])) +
            "\t" +
            "\t".join(sam_dict.get(KEY_NAMES[-1], [])), header)

    ########################################################
    ## Basic attributes in order of appearance in SAM format
    property query_name:
        """the query template name (None if not present)"""
        def __get__(self):

            cdef bam1_t * src = self._delegate
            if src.core.l_qname == 0:
                return None

            return charptr_to_str(<char *>pysam_bam_get_qname(src))

        def __set__(self, qname):

            if qname is None or len(qname) == 0:
                return

            if len(qname) > 254:
                raise ValueError("query length out of range {} > 254".format(
                    len(qname)))

            qname = force_bytes(qname)
            cdef bam1_t * src = self._delegate
            # the qname is \0 terminated
            cdef uint8_t l = len(qname) + 1

            cdef char * p = pysam_bam_get_qname(src)
            cdef uint8_t l_extranul = 0

            if l % 4 != 0:
                l_extranul = 4 - l % 4

            cdef bam1_t * retval = pysam_bam_update(src,
                                                    src.core.l_qname,
                                                    l + l_extranul,
                                                    <uint8_t*>p)
            if retval == NULL:
                raise MemoryError("could not allocate memory")

            src.core.l_extranul = l_extranul
            src.core.l_qname = l + l_extranul

            # re-acquire pointer to location in memory
            # as it might have moved
            p = pysam_bam_get_qname(src)

            strncpy(p, qname, l)
            # x might be > 255
            cdef uint16_t x = 0

            for x from l <= x < l + l_extranul:
                p[x] = b'\0'

    property flag:
        """properties flag"""
        def __get__(self):
            return self._delegate.core.flag
        def __set__(self, flag):
            self._delegate.core.flag = flag

    property reference_name:
        """:term:`reference` name"""
        def __get__(self):
            if self._delegate.core.tid == -1:
                return None
            if self.header:
                return self.header.get_reference_name(self._delegate.core.tid)
            else:
                raise ValueError("reference_name unknown if no header associated with record")
        def __set__(self, reference):
            cdef int tid
            if reference is None or reference == "*":
                self._delegate.core.tid = -1
            elif self.header:
                tid = self.header.get_tid(reference)
                if tid < 0:
                    raise ValueError("reference {} does not exist in header".format(
                        reference))
                self._delegate.core.tid = tid
            else:
                raise ValueError("reference_name can not be set if no header associated with record")

    property reference_id:
        """:term:`reference` ID

        .. note::

            This field contains the index of the reference sequence in
            the sequence dictionary. To obtain the name of the
            reference sequence, use :meth:`get_reference_name()`

        """
        def __get__(self):
            return self._delegate.core.tid
        def __set__(self, tid):
            if tid != -1 and self.header and not self.header.is_valid_tid(tid):
                raise ValueError("reference id {} does not exist in header".format(
                    tid))
            self._delegate.core.tid = tid

    property reference_start:
        """0-based leftmost coordinate"""
        def __get__(self):
            return self._delegate.core.pos
        def __set__(self, pos):
            ## setting the position requires updating the "bin" attribute
            cdef bam1_t * src
            src = self._delegate
            src.core.pos = pos
            update_bin(src)

    property mapping_quality:
        """mapping quality"""
        def __get__(self):
            return pysam_get_qual(self._delegate)
        def __set__(self, qual):
            pysam_set_qual(self._delegate, qual)

    property cigarstring:
        '''the :term:`cigar` alignment as a string.

        The cigar string is a string of alternating integers
        and characters denoting the length and the type of
        an operation.

        .. note::
            The order length,operation is specified in the
            SAM format. It is different from the order of
            the :attr:`cigar` property.

        Returns None if not present.

        To unset the cigarstring, assign None or the
        empty string.
        '''
        def __get__(self):
            c = self.cigartuples
            if c is None:
                return None
            # reverse order
            else:
                return "".join([ "%i%c" % (y,CODE2CIGAR[x]) for x,y in c])

        def __set__(self, cigar):
            if cigar is None or len(cigar) == 0:
                self.cigartuples = []
            else:
                parts = CIGAR_REGEX.findall(cigar)
                # reverse order
                self.cigartuples = [(CIGAR2CODE[ord(y)], int(x)) for x,y in parts]

    # TODO
    # property cigar:
    #     """the cigar alignment"""

    property next_reference_id:
        """the :term:`reference` id of the mate/next read."""
        def __get__(self):
            return self._delegate.core.mtid
        def __set__(self, mtid):
            if mtid != -1 and self.header and not self.header.is_valid_tid(mtid):
                raise ValueError("reference id {} does not exist in header".format(
                    mtid))
            self._delegate.core.mtid = mtid

    property next_reference_name:
        """:term:`reference` name of the mate/next read (None if no
        AlignmentFile is associated)"""
        def __get__(self):
            if self._delegate.core.mtid == -1:
                return None
            if self.header:
                return self.header.get_reference_name(self._delegate.core.mtid)
            else:
                raise ValueError("next_reference_name unknown if no header associated with record")

        def __set__(self, reference):
            cdef int mtid
            if reference is None or reference == "*":
                self._delegate.core.mtid = -1
            elif reference == "=":
                self._delegate.core.mtid = self._delegate.core.tid
            elif self.header:
                mtid = self.header.get_tid(reference)
                if mtid < 0:
                    raise ValueError("reference {} does not exist in header".format(
                        reference))
                self._delegate.core.mtid = mtid
            else:
                raise ValueError("next_reference_name can not be set if no header associated with record")

    property next_reference_start:
        """the position of the mate/next read."""
        def __get__(self):
            return self._delegate.core.mpos
        def __set__(self, mpos):
            self._delegate.core.mpos = mpos

    property query_length:
        """the length of the query/read.

        This value corresponds to the length of the sequence supplied
        in the BAM/SAM file. The length of a query is 0 if there is no
        sequence in the BAM/SAM file. In those cases, the read length
        can be inferred from the CIGAR alignment, see
        :meth:`pysam.AlignedSegment.infer_query_length`.

        The length includes soft-clipped bases and is equal to
        ``len(query_sequence)``.

        This property is read-only but is updated when a new query
        sequence is assigned to this AlignedSegment.

        Returns 0 if not available.

        """
        def __get__(self):
            return self._delegate.core.l_qseq

    property template_length:
        """the observed query template length"""
        def __get__(self):
            return self._delegate.core.isize
        def __set__(self, isize):
            self._delegate.core.isize = isize

    property query_sequence:
        """read sequence bases, including :term:`soft clipped` bases
        (None if not present).

        Assigning to this attribute will invalidate any quality scores.
        Thus, to in-place edit the sequence and quality scores, copies of
        the quality scores need to be taken. Consider trimming for example::

           q = read.query_qualities
           read.query_sequence = read.query_sequence[5:10]
           read.query_qualities = q[5:10]

        The sequence is returned as it is stored in the BAM file. (This will
        be the reverse complement of the original read sequence if the mapper
        has aligned the read to the reverse strand.)
        """
        def __get__(self):
            if self.cache_query_sequence:
                return self.cache_query_sequence

            cdef bam1_t * src
            cdef char * s
            src = self._delegate

            if src.core.l_qseq == 0:
                return None

            self.cache_query_sequence = force_str(getSequenceInRange(
                src, 0, src.core.l_qseq))
            return self.cache_query_sequence

        def __set__(self, seq):
            # samtools manages sequence and quality length memory together
            # if no quality information is present, the first byte says 0xff.
            cdef bam1_t * src
            cdef uint8_t * p
            cdef char * s
            cdef int l, k
            cdef Py_ssize_t nbytes_new, nbytes_old

            if seq == None:
                l = 0
            else:
                l = len(seq)
                seq = force_bytes(seq)

            src = self._delegate

            # as the sequence is stored in half-bytes, the total length (sequence
            # plus quality scores) is (l+1)/2 + l
            nbytes_new = (l + 1) // 2 + l
            nbytes_old = (src.core.l_qseq + 1) // 2 + src.core.l_qseq

            # acquire pointer to location in memory
            p = pysam_bam_get_seq(src)
            src.core.l_qseq = l

            # change length of data field
            cdef bam1_t * retval = pysam_bam_update(src,
                                                    nbytes_old,
                                                    nbytes_new,
                                                    p)

            if retval == NULL:
                raise MemoryError("could not allocate memory")

            if l > 0:
                # re-acquire pointer to location in memory
                # as it might have moved
                p = pysam_bam_get_seq(src)
                for k from 0 <= k < nbytes_new:
                    p[k] = 0
                # convert to C string
                s = seq
                for k from 0 <= k < l:
                    p[k // 2] |= seq_nt16_table[<unsigned char>s[k]] << 4 * (1 - k % 2)

                # erase qualities
                p = pysam_bam_get_qual(src)
                memset(p, 0xff, l)

            self.cache_query_sequence = force_str(seq)

            # clear cached values for quality values
            self.cache_query_qualities = None
            self.cache_query_alignment_qualities = None

    property query_qualities:
        """read sequence base qualities, including :term:`soft clipped` bases 
        (None if not present).

        Quality scores are returned as a python array of unsigned
        chars. Note that this is not the ASCII-encoded value typically
        seen in FASTQ or SAM formatted files. Thus, no offset of 33
        needs to be subtracted.

        Note that to set quality scores the sequence has to be set
        beforehand as this will determine the expected length of the
        quality score array.

        This method raises a ValueError if the length of the
        quality scores and the sequence are not the same.

        """
        def __get__(self):

            if self.cache_query_qualities:
                return self.cache_query_qualities

            cdef bam1_t * src
            cdef char * q

            src = self._delegate

            if src.core.l_qseq == 0:
                return None

            self.cache_query_qualities = getQualitiesInRange(src, 0, src.core.l_qseq)
            return self.cache_query_qualities

        def __set__(self, qual):

            # note that memory is already allocated via setting the sequence
            # hence length match of sequence and quality needs is checked.
            cdef bam1_t * src
            cdef uint8_t * p
            cdef int l

            src = self._delegate
            p = pysam_bam_get_qual(src)
            if qual is None or len(qual) == 0:
                # if absent and there is a sequence: set to 0xff
                memset(p, 0xff, src.core.l_qseq)
                return

            # check for length match
            l = len(qual)
            if src.core.l_qseq != l:
                raise ValueError(
                    "quality and sequence mismatch: %i != %i" %
                    (l, src.core.l_qseq))

            # create a python array object filling it
            # with the quality scores

            # NB: should avoid this copying if qual is
            # already of the correct type.
            cdef c_array.array result = c_array.array('B', qual)

            # copy data
            memcpy(p, result.data.as_voidptr, l)

            # save in cache
            self.cache_query_qualities = qual

    property bin:
        """properties bin"""
        def __get__(self):
            return self._delegate.core.bin
        def __set__(self, bin):
            self._delegate.core.bin = bin


    ##########################################################
    # Derived simple attributes. These are simple attributes of
    # AlignedSegment getting and setting values.
    ##########################################################
    # 1. Flags
    ##########################################################
    property is_paired:
        """true if read is paired in sequencing"""
        def __get__(self):
            return (self.flag & BAM_FPAIRED) != 0
        def __set__(self,val):
            pysam_update_flag(self._delegate, val, BAM_FPAIRED)

    property is_proper_pair:
        """true if read is mapped in a proper pair"""
        def __get__(self):
            return (self.flag & BAM_FPROPER_PAIR) != 0
        def __set__(self,val):
            pysam_update_flag(self._delegate, val, BAM_FPROPER_PAIR)

    property is_unmapped:
        """true if read itself is unmapped"""
        def __get__(self):
            return (self.flag & BAM_FUNMAP) != 0
        def __set__(self, val):
            pysam_update_flag(self._delegate, val, BAM_FUNMAP)
            # setting the unmapped flag requires recalculation of
            # bin as alignment length is now implicitly 1
            update_bin(self._delegate)

    property is_mapped:
        """true if read itself is mapped
        (implemented in terms of :attr:`is_unmapped`)"""
        def __get__(self):
            return (self.flag & BAM_FUNMAP) == 0
        def __set__(self, val):
            pysam_update_flag(self._delegate, not val, BAM_FUNMAP)
            update_bin(self._delegate)

    property mate_is_unmapped:
        """true if the mate is unmapped"""
        def __get__(self):
            return (self.flag & BAM_FMUNMAP) != 0
        def __set__(self,val):
            pysam_update_flag(self._delegate, val, BAM_FMUNMAP)

    property mate_is_mapped:
        """true if the mate is mapped
        (implemented in terms of :attr:`mate_is_unmapped`)"""
        def __get__(self):
            return (self.flag & BAM_FMUNMAP) == 0
        def __set__(self,val):
            pysam_update_flag(self._delegate, not val, BAM_FMUNMAP)

    property is_reverse:
        """true if read is mapped to reverse strand"""
        def __get__(self):
            return (self.flag & BAM_FREVERSE) != 0
        def __set__(self,val):
            pysam_update_flag(self._delegate, val, BAM_FREVERSE)

    property is_forward:
        """true if read is mapped to forward strand
        (implemented in terms of :attr:`is_reverse`)"""
        def __get__(self):
            return (self.flag & BAM_FREVERSE) == 0
        def __set__(self,val):
            pysam_update_flag(self._delegate, not val, BAM_FREVERSE)

    property mate_is_reverse:
        """true if the mate is mapped to reverse strand"""
        def __get__(self):
            return (self.flag & BAM_FMREVERSE) != 0
        def __set__(self,val):
            pysam_update_flag(self._delegate, val, BAM_FMREVERSE)

    property mate_is_forward:
        """true if the mate is mapped to forward strand
        (implemented in terms of :attr:`mate_is_reverse`)"""
        def __get__(self):
            return (self.flag & BAM_FMREVERSE) == 0
        def __set__(self,val):
            pysam_update_flag(self._delegate, not val, BAM_FMREVERSE)

    property is_read1:
        """true if this is read1"""
        def __get__(self):
            return (self.flag & BAM_FREAD1) != 0
        def __set__(self,val):
            pysam_update_flag(self._delegate, val, BAM_FREAD1)
    property is_read2:
        """true if this is read2"""
        def __get__(self):
            return (self.flag & BAM_FREAD2) != 0
        def __set__(self, val):
            pysam_update_flag(self._delegate, val, BAM_FREAD2)
    property is_secondary:
        """true if not primary alignment"""
        def __get__(self):
            return (self.flag & BAM_FSECONDARY) != 0
        def __set__(self, val):
            pysam_update_flag(self._delegate, val, BAM_FSECONDARY)
    property is_qcfail:
        """true if QC failure"""
        def __get__(self):
            return (self.flag & BAM_FQCFAIL) != 0
        def __set__(self, val):
            pysam_update_flag(self._delegate, val, BAM_FQCFAIL)
    property is_duplicate:
        """true if optical or PCR duplicate"""
        def __get__(self):
            return (self.flag & BAM_FDUP) != 0
        def __set__(self, val):
            pysam_update_flag(self._delegate, val, BAM_FDUP)
    property is_supplementary:
        """true if this is a supplementary alignment"""
        def __get__(self):
            return (self.flag & BAM_FSUPPLEMENTARY) != 0
        def __set__(self, val):
            pysam_update_flag(self._delegate, val, BAM_FSUPPLEMENTARY)

    # 2. Coordinates and lengths
    property reference_end:
        '''aligned reference position of the read on the reference genome.

        reference_end points to one past the last aligned residue.
        Returns None if not available (read is unmapped or no cigar
        alignment present).

        '''
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            if (self.flag & BAM_FUNMAP) or pysam_get_n_cigar(src) == 0:
                return None
            return bam_endpos(src)

    property reference_length:
        '''aligned length of the read on the reference genome.

        This is equal to `reference_end - reference_start`. 
        Returns None if not available.
        '''
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            if (self.flag & BAM_FUNMAP) or pysam_get_n_cigar(src) == 0:
                return None
            return bam_endpos(src) - \
                self._delegate.core.pos

    property query_alignment_sequence:
        """aligned portion of the read.

        This is a substring of :attr:`query_sequence` that excludes flanking
        bases that were :term:`soft clipped` (None if not present). It
        is equal to ``query_sequence[query_alignment_start:query_alignment_end]``.

        SAM/BAM files may include extra flanking bases that are not
        part of the alignment.  These bases may be the result of the
        Smith-Waterman or other algorithms, which may not require
        alignments that begin at the first residue or end at the last.
        In addition, extra sequencing adapters, multiplex identifiers,
        and low-quality bases that were not considered for alignment
        may have been retained.

        """

        def __get__(self):
            if self.cache_query_alignment_sequence:
                return self.cache_query_alignment_sequence

            cdef bam1_t * src
            cdef uint32_t start, end

            src = self._delegate

            if src.core.l_qseq == 0:
                return None

            start = getQueryStart(src)
            end   = getQueryEnd(src)

            self.cache_query_alignment_sequence = force_str(
                getSequenceInRange(src, start, end))
            return self.cache_query_alignment_sequence

    property query_alignment_qualities:
        """aligned query sequence quality values (None if not present). These
        are the quality values that correspond to 
        :attr:`query_alignment_sequence`, that is, they exclude qualities of 
        :term:`soft clipped` bases. This is equal to 
        ``query_qualities[query_alignment_start:query_alignment_end]``.

        Quality scores are returned as a python array of unsigned
        chars. Note that this is not the ASCII-encoded value typically
        seen in FASTQ or SAM formatted files. Thus, no offset of 33
        needs to be subtracted.

        This property is read-only.

        """
        def __get__(self):

            if self.cache_query_alignment_qualities:
                return self.cache_query_alignment_qualities

            cdef bam1_t * src
            cdef uint32_t start, end

            src = self._delegate

            if src.core.l_qseq == 0:
                return None

            start = getQueryStart(src)
            end   = getQueryEnd(src)
            self.cache_query_alignment_qualities = \
                getQualitiesInRange(src, start, end)
            return self.cache_query_alignment_qualities

    property query_alignment_start:
        """start index of the aligned query portion of the sequence (0-based,
        inclusive).

        This the index of the first base in :attr:`query_sequence` 
        that is not soft-clipped.
        """
        def __get__(self):
            return getQueryStart(self._delegate)

    property query_alignment_end:
        """end index of the aligned query portion of the sequence (0-based,
        exclusive)

        This the index just past the last base in :attr:`query_sequence` 
        that is not soft-clipped.
        """
        def __get__(self):
            return getQueryEnd(self._delegate)

    property modified_bases:
        """Modified bases annotations from Ml/Mm tags. The output is
        Dict[(canonical base, strand, modification)] -> [ (pos,qual), ...]
        with qual being (256*probability), or -1 if unknown.
        Strand==0 for forward and 1 for reverse strand modification
        """
        def __get__(self):
            cdef bam1_t * src
            cdef hts_base_mod_state *m = hts_base_mod_state_alloc()
            cdef hts_base_mod mods[5]
            cdef int pos

            ret = {}
            src = self._delegate
            
            if bam_parse_basemod(src, m) < 0:        
                return None
            
            n = bam_next_basemod(src, m, mods, 5, &pos)

            while n>0:
                for i in range(n):
                    mod_code = chr(mods[i].modified_base) if mods[i].modified_base>0 else -mods[i].modified_base
                    mod_strand = mods[i].strand
                    if self.is_reverse:
                        mod_strand = 1 - mod_strand
                    key = (chr(mods[i].canonical_base), 
                            mod_strand,
                            mod_code )
                    ret.setdefault(key,[]).append((pos,mods[i].qual))
                    
                n = bam_next_basemod(src, m, mods, 5, &pos)

            if n<0:
                return None

            hts_base_mod_state_free(m)
            return ret

    property modified_bases_forward:
        """Modified bases annotations from Ml/Mm tags. The output is
        Dict[(canonical base, strand, modification)] -> [ (pos,qual), ...]
        with qual being (256*probability), or -1 if unknown.
        Strand==0 for forward and 1 for reverse strand modification.
        The positions are with respect to the original sequence from get_forward_sequence()
        """
        def __get__(self):
            pmods = self.modified_bases
            if pmods and self.is_reverse:                
                rmod = {}

                # Try to find the length of the original sequence
                rlen = self.infer_read_length()
                if rlen is None and self.query_sequence is None:
                    return rmod
                else:
                    rlen = len(self.query_sequence)
                    
                for k,mods in pmods.items():
                    nk = k[0],1 - k[1],k[2]
                    for i in range(len(mods)):
                        
                        mods[i] = (rlen - 1 -mods[i][0], mods[i][1])
                    rmod[nk] = mods
                return rmod
            
            return pmods

 
    property query_alignment_length:
        """length of the aligned query sequence.

        This is equal to :attr:`query_alignment_end` - 
        :attr:`query_alignment_start`
        """
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            return getQueryEnd(src) - getQueryStart(src)

    #####################################################
    # Computed properties

    def get_reference_positions(self, full_length=False):
        """a list of reference positions that this read aligns to.

        By default, this method only returns positions in the
        reference that are within the alignment. If *full_length* is
        set, None values will be included for any soft-clipped or
        unaligned positions within the read. The returned list will
        thus be of the same length as the read.

        """
        cdef uint32_t k, i, l, pos
        cdef int op
        cdef uint32_t * cigar_p
        cdef bam1_t * src
        cdef bint _full = full_length

        src = self._delegate
        if pysam_get_n_cigar(src) == 0:
            return []

        result = []
        pos = src.core.pos
        cigar_p = pysam_bam_get_cigar(src)

        for k from 0 <= k < pysam_get_n_cigar(src):
            op = cigar_p[k] & BAM_CIGAR_MASK
            l = cigar_p[k] >> BAM_CIGAR_SHIFT

            if op == BAM_CSOFT_CLIP or op == BAM_CINS:
                if _full:
                    for i from 0 <= i < l:
                        result.append(None)
            elif op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
                for i from pos <= i < pos + l:
                    result.append(i)
                pos += l
            elif op == BAM_CDEL or op == BAM_CREF_SKIP:
                pos += l

        return result

    def infer_query_length(self, always=False):
        """infer query length from CIGAR alignment.

        This method deduces the query length from the CIGAR alignment
        but does not include hard-clipped bases.

        Returns None if CIGAR alignment is not present.

        If *always* is set to True, `infer_read_length` is used instead.
        This is deprecated and only present for backward compatibility.
        """
        if always is True:
            return self.infer_read_length()
        cdef int32_t l = calculateQueryLengthWithoutHardClipping(self._delegate)
        if l > 0:
            return l
        else:
            return None

    def infer_read_length(self):
        """infer read length from CIGAR alignment.

        This method deduces the read length from the CIGAR alignment
        including hard-clipped bases.

        Returns None if CIGAR alignment is not present.
        """
        cdef int32_t l = calculateQueryLengthWithHardClipping(self._delegate)
        if l > 0:
            return l
        else:
            return None

    def get_reference_sequence(self):
        """return the reference sequence in the region that is covered by the
        alignment of the read to the reference.

        This method requires the MD tag to be set.

        """
        return force_str(build_reference_sequence(self._delegate))

    def get_forward_sequence(self):
        """return the original read sequence.

        Reads mapped to the reverse strand are stored reverse complemented in
        the BAM file. This method returns such reads reverse complemented back
        to their original orientation.

        Returns None if the record has no query sequence.
        """
        if self.query_sequence is None:
            return None
        s = force_str(self.query_sequence)
        if self.is_reverse:
            s = s.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))[::-1]
        return s

    def get_forward_qualities(self):
        """return the original base qualities of the read sequence,
        in the same format as the :attr:`query_qualities` property.

        Reads mapped to the reverse strand have their base qualities stored
        reversed in the BAM file. This method returns such reads' base qualities
        reversed back to their original orientation.
        """
        if self.is_reverse:
            return self.query_qualities[::-1]
        else:
            return self.query_qualities


    def get_aligned_pairs(self, matches_only=False, with_seq=False):
        """a list of aligned read (query) and reference positions.

        For inserts, deletions, skipping either query or reference
        position may be None.

        For padding in the reference, the reference position will
        always be None.

        Parameters
        ----------

        matches_only : bool
          If True, only matched bases are returned - no None on either
          side.
        with_seq : bool
          If True, return a third element in the tuple containing the
          reference sequence. For CIGAR 'P' (padding in the reference)
          operations, the third tuple element will be None. Substitutions
          are lower-case. This option requires an MD tag to be present.

        Returns
        -------

        aligned_pairs : list of tuples

        """
        cdef uint32_t k, i, pos, qpos, r_idx, l
        cdef int op
        cdef uint32_t * cigar_p
        cdef bam1_t * src = self._delegate
        cdef bint _matches_only = bool(matches_only)
        cdef bint _with_seq = bool(with_seq)

        # TODO: this method performs no checking and assumes that
        # read sequence, cigar and MD tag are consistent.

        if _with_seq:
            # force_str required for py2/py3 compatibility
            ref_seq = force_str(build_reference_sequence(src))
            if ref_seq is None:
                raise ValueError("MD tag not present")

        r_idx = 0

        if pysam_get_n_cigar(src) == 0:
            return []

        result = []
        pos = src.core.pos
        qpos = 0
        cigar_p = pysam_bam_get_cigar(src)
        for k from 0 <= k < pysam_get_n_cigar(src):
            op = cigar_p[k] & BAM_CIGAR_MASK
            l = cigar_p[k] >> BAM_CIGAR_SHIFT

            if op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
                if _with_seq:
                    for i from pos <= i < pos + l:
                        result.append((qpos, i, ref_seq[r_idx]))
                        r_idx += 1
                        qpos += 1
                else:
                    for i from pos <= i < pos + l:
                        result.append((qpos, i))
                        qpos += 1
                pos += l

            elif op == BAM_CINS or op == BAM_CSOFT_CLIP or op == BAM_CPAD:
                if not _matches_only:
                    if _with_seq:
                        for i from pos <= i < pos + l:
                            result.append((qpos, None, None))
                            qpos += 1
                    else:
                        for i from pos <= i < pos + l:
                            result.append((qpos, None))
                            qpos += 1
                else:
                    qpos += l

            elif op == BAM_CDEL:
                if not _matches_only:
                    if _with_seq:
                        for i from pos <= i < pos + l:
                            result.append((None, i, ref_seq[r_idx]))
                            r_idx += 1
                    else:
                        for i from pos <= i < pos + l:
                            result.append((None, i))
                else:
                    r_idx += l
                pos += l

            elif op == BAM_CHARD_CLIP:
                pass # advances neither

            elif op == BAM_CREF_SKIP:
                if not _matches_only:
                    if _with_seq:
                        for i from pos <= i < pos + l:
                            result.append((None, i, None))
                    else:
                        for i from pos <= i < pos + l:
                            result.append((None, i))

                pos += l

        return result

    def get_blocks(self):
        """ a list of start and end positions of
        aligned gapless blocks.

        The start and end positions are in genomic
        coordinates.

        Blocks are not normalized, i.e. two blocks
        might be directly adjacent. This happens if
        the two blocks are separated by an insertion
        in the read.
        """

        cdef uint32_t k, pos, l
        cdef int op
        cdef uint32_t * cigar_p
        cdef bam1_t * src

        src = self._delegate
        if pysam_get_n_cigar(src) == 0:
            return []

        result = []
        pos = src.core.pos
        cigar_p = pysam_bam_get_cigar(src)
        l = 0

        for k from 0 <= k < pysam_get_n_cigar(src):
            op = cigar_p[k] & BAM_CIGAR_MASK
            l = cigar_p[k] >> BAM_CIGAR_SHIFT
            if op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
                result.append((pos, pos + l))
                pos += l
            elif op == BAM_CDEL or op == BAM_CREF_SKIP:
                pos += l

        return result

    def get_overlap(self, uint32_t start, uint32_t end):
        """return number of aligned bases of read overlapping the interval
        *start* and *end* on the reference sequence.

        Return None if cigar alignment is not available.
        """
        cdef uint32_t k, i, pos, overlap
        cdef int op, o
        cdef uint32_t * cigar_p
        cdef bam1_t * src

        overlap = 0

        src = self._delegate
        if pysam_get_n_cigar(src) == 0:
            return None
        pos = src.core.pos
        o = 0

        cigar_p = pysam_bam_get_cigar(src)
        for k from 0 <= k < pysam_get_n_cigar(src):
            op = cigar_p[k] & BAM_CIGAR_MASK
            l = cigar_p[k] >> BAM_CIGAR_SHIFT

            if op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
                o = min( pos + l, end) - max( pos, start )
                if o > 0: overlap += o

            if op == BAM_CMATCH or op == BAM_CDEL or op == BAM_CREF_SKIP or op == BAM_CEQUAL or op == BAM_CDIFF:
                pos += l

        return overlap

    def get_cigar_stats(self):
        """summary of operations in cigar string.

        The output order in the array is "MIDNSHP=X" followed by a
        field for the NM tag. If the NM tag is not present, this
        field will always be 0.

        +-----+--------------+-----+
        |M    |BAM_CMATCH    |0    |
        +-----+--------------+-----+
        |I    |BAM_CINS      |1    |
        +-----+--------------+-----+
        |D    |BAM_CDEL      |2    |
        +-----+--------------+-----+
        |N    |BAM_CREF_SKIP |3    |
        +-----+--------------+-----+
        |S    |BAM_CSOFT_CLIP|4    |
        +-----+--------------+-----+
        |H    |BAM_CHARD_CLIP|5    |
        +-----+--------------+-----+
        |P    |BAM_CPAD      |6    |
        +-----+--------------+-----+
        |=    |BAM_CEQUAL    |7    |
        +-----+--------------+-----+
        |X    |BAM_CDIFF     |8    |
        +-----+--------------+-----+
        |B    |BAM_CBACK     |9    |
        +-----+--------------+-----+
        |NM   |NM tag        |10   |
        +-----+--------------+-----+

        If no cigar string is present, empty arrays will be returned.

        Returns:
            arrays :
                two arrays. The first contains the nucleotide counts within
                each cigar operation, the second contains the number of blocks
                for each cigar operation.

        """

        cdef int nfields = NCIGAR_CODES + 1

        cdef c_array.array base_counts = array.array(
            "I",
            [0] * nfields)
        cdef uint32_t [:] base_view = base_counts
        cdef c_array.array block_counts = array.array(
            "I",
            [0] * nfields)
        cdef uint32_t [:] block_view = block_counts

        cdef bam1_t * src = self._delegate
        cdef int op
        cdef uint32_t l
        cdef int32_t k
        cdef uint32_t * cigar_p = pysam_bam_get_cigar(src)

        if cigar_p == NULL:
            return None

        for k from 0 <= k < pysam_get_n_cigar(src):
            op = cigar_p[k] & BAM_CIGAR_MASK
            l = cigar_p[k] >> BAM_CIGAR_SHIFT
            base_view[op] += l
            block_view[op] += 1

        cdef uint8_t * v = bam_aux_get(src, 'NM')
        if v != NULL:
            base_view[nfields - 1] = <int32_t>bam_aux2i(v)

        return base_counts, block_counts

    #####################################################
    ## Unsorted as yet
    # TODO: capture in CIGAR object
    property cigartuples:
        """the :term:`cigar` alignment. The alignment
        is returned as a list of tuples of (operation, length).

        If the alignment is not present, None is returned.

        The operations are:

        +-----+--------------+-----+
        |M    |BAM_CMATCH    |0    |
        +-----+--------------+-----+
        |I    |BAM_CINS      |1    |
        +-----+--------------+-----+
        |D    |BAM_CDEL      |2    |
        +-----+--------------+-----+
        |N    |BAM_CREF_SKIP |3    |
        +-----+--------------+-----+
        |S    |BAM_CSOFT_CLIP|4    |
        +-----+--------------+-----+
        |H    |BAM_CHARD_CLIP|5    |
        +-----+--------------+-----+
        |P    |BAM_CPAD      |6    |
        +-----+--------------+-----+
        |=    |BAM_CEQUAL    |7    |
        +-----+--------------+-----+
        |X    |BAM_CDIFF     |8    |
        +-----+--------------+-----+
        |B    |BAM_CBACK     |9    |
        +-----+--------------+-----+

        .. note::
            The output is a list of (operation, length) tuples, such as
            ``[(0, 30)]``.
            This is different from the SAM specification and
            the :attr:`cigarstring` property, which uses a
            (length, operation) order, for example: ``30M``.

        To unset the cigar property, assign an empty list
        or None.
        """
        def __get__(self):
            cdef uint32_t * cigar_p
            cdef bam1_t * src
            cdef uint32_t op, l
            cdef uint32_t k

            src = self._delegate
            if pysam_get_n_cigar(src) == 0:
                return None

            cigar = []

            cigar_p = pysam_bam_get_cigar(src);
            for k from 0 <= k < pysam_get_n_cigar(src):
                op = cigar_p[k] & BAM_CIGAR_MASK
                l = cigar_p[k] >> BAM_CIGAR_SHIFT
                cigar.append((op, l))
            return cigar

        def __set__(self, values):
            cdef uint32_t * p
            cdef bam1_t * src
            cdef op, l
            cdef int k

            k = 0

            src = self._delegate

            # get location of cigar string
            p = pysam_bam_get_cigar(src)

            # empty values for cigar string
            if values is None:
                values = []

            cdef uint32_t ncigar = len(values)

            cdef bam1_t * retval = pysam_bam_update(src,
                                                    pysam_get_n_cigar(src) * 4,
                                                    ncigar * 4,
                                                    <uint8_t*>p)

            if retval == NULL:
                raise MemoryError("could not allocate memory")

            # length is number of cigar operations, not bytes
            pysam_set_n_cigar(src, ncigar)

            # re-acquire pointer to location in memory
            # as it might have moved
            p = pysam_bam_get_cigar(src)

            # insert cigar operations
            for op, l in values:
                p[k] = l << BAM_CIGAR_SHIFT | op
                k += 1

            ## setting the cigar string requires updating the bin
            update_bin(src)

    cpdef set_tag(self,
                  tag,
                  value,
                  value_type=None,
                  replace=True):
        """sets a particular field *tag* to *value* in the optional alignment
        section.

        *value_type* describes the type of *value* that is to entered
        into the alignment record. It can be set explicitly to one of
        the valid one-letter type codes. If unset, an appropriate type
        will be chosen automatically based on the python type of
        *value*.

        An existing value of the same *tag* will be overwritten unless
        *replace* is set to False. This is usually not recommended as a
        tag may only appear once in the optional alignment section.

        If *value* is `None`, the tag will be deleted.

        This method accepts valid SAM specification value types, which
        are::

           A: printable char
           i: signed int
           f: float
           Z: printable string
           H: Byte array in hex format
           B: Integer or numeric array

        Additionally, it will accept the integer BAM types ('cCsSI')

        For htslib compatibility, 'a' is synonymous with 'A' and the
        method accepts a 'd' type code for a double precision float.

        When deducing the type code by the python type of *value*, the
        following mapping is applied::

            i: python int
            f: python float
            Z: python str or bytes
            B: python array.array, list or tuple

        Note that a single character string will be output as 'Z' and
        not 'A' as the former is the more general type.
        """

        cdef int value_size
        cdef uint8_t tc
        cdef uint8_t * value_ptr
        cdef uint8_t *existing_ptr
        cdef float float_value
        cdef double double_value
        cdef int32_t int32_t_value
        cdef uint32_t uint32_t_value
        cdef int16_t int16_t_value
        cdef uint16_t uint16_t_value
        cdef int8_t int8_t_value
        cdef uint8_t uint8_t_value
        cdef bam1_t * src = self._delegate
        cdef char * _value_type
        cdef c_array.array array_value
        cdef object buffer

        if len(tag) != 2:
            raise ValueError('Invalid tag: %s' % tag)

        tag = force_bytes(tag)
        if replace:
            existing_ptr = bam_aux_get(src, tag)
            if existing_ptr:
                bam_aux_del(src, existing_ptr)

        # setting value to None deletes a tag
        if value is None:
            return

        cdef uint8_t typecode = get_tag_typecode(value, value_type)
        if typecode == 0:
            raise ValueError("can't guess type or invalid type code specified: {} {}".format(
                value, value_type))

        # sam_format1 for typecasting
        if typecode == b'Z':
            value = force_bytes(value)
            value_ptr = <uint8_t*><char*>value
            value_size = len(value)+1
        elif typecode == b'H':
            # Note that hex tags are stored the very same
            # way as Z string.s
            value = force_bytes(value)
            value_ptr = <uint8_t*><char*>value
            value_size = len(value)+1
        elif typecode == b'A' or typecode == b'a':
            value = force_bytes(value)
            value_ptr = <uint8_t*><char*>value
            value_size = sizeof(char)
            typecode = b'A'
        elif typecode == b'i':
            int32_t_value = value
            value_ptr = <uint8_t*>&int32_t_value
            value_size = sizeof(int32_t)
        elif typecode == b'I':
            uint32_t_value = value
            value_ptr = <uint8_t*>&uint32_t_value
            value_size = sizeof(uint32_t)
        elif typecode == b's':
            int16_t_value = value
            value_ptr = <uint8_t*>&int16_t_value
            value_size = sizeof(int16_t)
        elif typecode == b'S':
            uint16_t_value = value
            value_ptr = <uint8_t*>&uint16_t_value
            value_size = sizeof(uint16_t)
        elif typecode == b'c':
            int8_t_value = value
            value_ptr = <uint8_t*>&int8_t_value
            value_size = sizeof(int8_t)
        elif typecode == b'C':
            uint8_t_value = value
            value_ptr = <uint8_t*>&uint8_t_value
            value_size = sizeof(uint8_t)
        elif typecode == b'd':
            double_value = value
            value_ptr = <uint8_t*>&double_value
            value_size = sizeof(double)
        elif typecode == b'f':
            float_value  = value
            value_ptr = <uint8_t*>&float_value
            value_size = sizeof(float)
        elif typecode == b'B':
            # the following goes through python, needs to be cleaned up
            # pack array using struct
            fmt, args = pack_tags([(tag, value, value_type)])

            # remove tag and type code as set by bam_aux_append
            # first four chars of format (<2sB)
            fmt = '<' + fmt[4:]
            # first two values to pack
            args = args[2:]
            value_size = struct.calcsize(fmt)
            # buffer will be freed when object goes out of scope
            buffer = ctypes.create_string_buffer(value_size)
            struct.pack_into(fmt, buffer, 0, *args)
            # bam_aux_append copies data from value_ptr
            bam_aux_append(src,
                           tag,
                           typecode,
                           value_size,
                           <uint8_t*>buffer.raw)
            return
        else:
            raise ValueError('unsupported value_type {} in set_option'.format(typecode))

        bam_aux_append(src,
                       tag,
                       typecode,
                       value_size,
                       value_ptr)

    cpdef has_tag(self, tag):
        """returns true if the optional alignment section
        contains a given *tag*."""
        cdef uint8_t * v
        cdef int nvalues
        btag = force_bytes(tag)
        v = bam_aux_get(self._delegate, btag)
        return v != NULL

    cpdef get_tag(self, tag, with_value_type=False):
        """
        retrieves data from the optional alignment section
        given a two-letter *tag* denoting the field.

        The returned value is cast into an appropriate python type.

        This method is the fastest way to access the optional
        alignment section if only few tags need to be retrieved.

        Possible value types are "AcCsSiIfZHB" (see BAM format
        specification) as well as additional value type 'd' as
        implemented in htslib.

        Parameters:

            tag :
                data tag.

            with_value_type : Optional[bool]
                if set to True, the return value is a tuple of (tag value, type
                code). (default False)

        Returns:

            A python object with the value of the `tag`. The type of the
            object depends on the data type in the data record.

        Raises:

            KeyError
                If `tag` is not present, a KeyError is raised.

        """
        cdef uint8_t * v
        cdef int nvalues
        btag = force_bytes(tag)
        v = bam_aux_get(self._delegate, btag)
        if v == NULL:
            raise KeyError("tag '%s' not present" % tag)
        if chr(v[0]) == "B":
            auxtype = chr(v[0]) + chr(v[1])
        else:
            auxtype = chr(v[0])

        if auxtype in "iIcCsS":
            value = bam_aux2i(v)
        elif auxtype == 'f' or auxtype == 'F':
            value = bam_aux2f(v)
        elif auxtype == 'd' or auxtype == 'D':
            value = bam_aux2f(v)
        elif auxtype == 'A' or auxtype == 'a':
            # force A to a
            v[0] = b'A'
            # there might a more efficient way
            # to convert a char into a string
            value = '%c' % <char>bam_aux2A(v)
        elif auxtype == 'Z' or auxtype == 'H':
            # Z and H are treated equally as strings in htslib
            value = charptr_to_str(<char*>bam_aux2Z(v))
        elif auxtype[0] == 'B':
            bytesize, nvalues, values = convert_binary_tag(v + 1)
            value = values
        else:
            raise ValueError("unknown auxiliary type '%s'" % auxtype)

        if with_value_type:
            return (value, auxtype)
        else:
            return value

    def get_tags(self, with_value_type=False):
        """the fields in the optional alignment section.

        Returns a list of all fields in the optional
        alignment section. Values are converted to appropriate python
        values. For example: ``[(NM, 2), (RG, "GJP00TM04")]``

        If *with_value_type* is set, the value type as encode in
        the AlignedSegment record will be returned as well:

        [(NM, 2, "i"), (RG, "GJP00TM04", "Z")]

        This method will convert all values in the optional alignment
        section. When getting only one or few tags, please see
        :meth:`get_tag` for a quicker way to achieve this.

        """

        cdef char * ctag
        cdef bam1_t * src
        cdef uint8_t * s
        cdef char auxtag[3]
        cdef char auxtype
        cdef uint8_t byte_size
        cdef int32_t nvalues

        src = self._delegate
        if src.l_data == 0:
            return []
        s = pysam_bam_get_aux(src)
        result = []
        auxtag[2] = 0
        while s < (src.data + src.l_data):
            # get tag
            auxtag[0] = s[0]
            auxtag[1] = s[1]
            s += 2
            auxtype = s[0]
            if auxtype in (b'c', b'C'):
                value = <int>bam_aux2i(s)
                s += 1
            elif auxtype in (b's', b'S'):
                value = <int>bam_aux2i(s)
                s += 2
            elif auxtype in (b'i', b'I'):
                value = <int32_t>bam_aux2i(s)
                s += 4
            elif auxtype == b'f':
                value = <float>bam_aux2f(s)
                s += 4
            elif auxtype == b'd':
                value = <double>bam_aux2f(s)
                s += 8
            elif auxtype in (b'A', b'a'):
                value = "%c" % <char>bam_aux2A(s)
                s += 1
            elif auxtype in (b'Z', b'H'):
                value = charptr_to_str(<char*>bam_aux2Z(s))
                # +1 for NULL terminated string
                s += len(value) + 1
            elif auxtype == b'B':
                s += 1
                byte_size, nvalues, value = convert_binary_tag(s)
                # 5 for 1 char and 1 int
                s += 5 + (nvalues * byte_size) - 1
            else:
                raise KeyError("unknown type '%s'" % auxtype)

            s += 1

            if with_value_type:
                result.append((charptr_to_str(auxtag), value, chr(auxtype)))
            else:
                result.append((charptr_to_str(auxtag), value))

        return result

    def set_tags(self, tags):
        """sets the fields in the optional alignment section with
        a list of (tag, value) tuples.

        The value type of the values is determined from the
        python type. Optionally, a type may be given explicitly as
        a third value in the tuple, For example:

        x.set_tags([(NM, 2, "i"), (RG, "GJP00TM04", "Z")]

        This method will not enforce the rule that the same tag may appear
        only once in the optional alignment section.
        """

        cdef bam1_t * src
        cdef uint8_t * s
        cdef char * temp
        cdef int new_size = 0
        cdef int old_size
        src = self._delegate

        # convert and pack the data
        if tags is not None and len(tags) > 0:
            fmt, args = pack_tags(tags)
            new_size = struct.calcsize(fmt)
            buffer = ctypes.create_string_buffer(new_size)
            struct.pack_into(fmt,
                             buffer,
                             0,
                             *args)


        # delete the old data and allocate new space.
        # If total_size == 0, the aux field will be
        # empty
        old_size = pysam_bam_get_l_aux(src)
        cdef bam1_t * retval = pysam_bam_update(src,
                                                old_size,
                                                new_size,
                                                pysam_bam_get_aux(src))
        if retval == NULL:
            raise MemoryError("could not allocate memory")

        # copy data only if there is any
        if new_size > 0:

            # get location of new data
            s = pysam_bam_get_aux(src)

            # check if there is direct path from buffer.raw to tmp
            p = buffer.raw
            # create handle to make sure buffer stays alive long
            # enough for memcpy, see issue 129
            temp = p
            memcpy(s, temp, new_size)


    ########################################################
    # Compatibility Accessors
    # Functions, properties for compatibility with pysam < 0.8
    #
    # Several options
    #     change the factory functions according to API
    #         * requires code changes throughout, incl passing
    #           handles to factory functions
    #     subclass functions and add attributes at runtime
    #         e.g.: AlignedSegments.qname = AlignedSegments.query_name
    #         * will slow down the default interface
    #     explicit declaration of getters/setters
    ########################################################
    property qname:
        """deprecated, use :attr:`query_name` instead."""
        def __get__(self): return self.query_name
        def __set__(self, v): self.query_name = v
    property tid:
        """deprecated, use :attr:`reference_id` instead."""
        def __get__(self): return self.reference_id
        def __set__(self, v): self.reference_id = v
    property pos:
        """deprecated, use :attr:`reference_start` instead."""
        def __get__(self): return self.reference_start
        def __set__(self, v): self.reference_start = v
    property mapq:
        """deprecated, use :attr:`mapping_quality` instead."""
        def __get__(self): return self.mapping_quality
        def __set__(self, v): self.mapping_quality = v
    property rnext:
        """deprecated, use :attr:`next_reference_id` instead."""
        def __get__(self): return self.next_reference_id
        def __set__(self, v): self.next_reference_id = v
    property pnext:
        """deprecated, use :attr:`next_reference_start` instead."""
        def __get__(self):
            return self.next_reference_start
        def __set__(self, v):
            self.next_reference_start = v
    property cigar:
        """deprecated, use :attr:`cigarstring` or :attr:`cigartuples` instead."""
        def __get__(self):
            r = self.cigartuples
            if r is None:
                r = []
            return r
        def __set__(self, v): self.cigartuples = v
    property tlen:
        """deprecated, use :attr:`template_length` instead."""
        def __get__(self):
            return self.template_length
        def __set__(self, v):
            self.template_length = v
    property seq:
        """deprecated, use :attr:`query_sequence` instead."""
        def __get__(self):
            return self.query_sequence
        def __set__(self, v):
            self.query_sequence = v
    property qual:
        """deprecated, use :attr:`query_qualities` instead."""
        def __get__(self):
            return array_to_qualitystring(self.query_qualities)
        def __set__(self, v):
            self.query_qualities = qualitystring_to_array(v)
    property alen:
        """deprecated, use :attr:`reference_length` instead."""
        def __get__(self):
            return self.reference_length
        def __set__(self, v):
            self.reference_length = v
    property aend:
        """deprecated, use :attr:`reference_end` instead."""
        def __get__(self):
            return self.reference_end
        def __set__(self, v):
            self.reference_end = v
    property rlen:
        """deprecated, use :attr:`query_length` instead."""
        def __get__(self):
            return self.query_length
        def __set__(self, v):
            self.query_length = v
    property query:
        """deprecated, use :attr:`query_alignment_sequence` 
        instead."""
        def __get__(self):
            return self.query_alignment_sequence
        def __set__(self, v):
            self.query_alignment_sequence = v
    property qqual:
        """deprecated, use :attr:`query_alignment_qualities` 
        instead."""
        def __get__(self):
            return array_to_qualitystring(self.query_alignment_qualities)
        def __set__(self, v):
            self.query_alignment_qualities = qualitystring_to_array(v)
    property qstart:
        """deprecated, use :attr:`query_alignment_start` instead."""
        def __get__(self):
            return self.query_alignment_start
        def __set__(self, v):
            self.query_alignment_start = v
    property qend:
        """deprecated, use :attr:`query_alignment_end` instead."""
        def __get__(self):
            return self.query_alignment_end
        def __set__(self, v):
            self.query_alignment_end = v
    property qlen:
        """deprecated, use :attr:`query_alignment_length` 
        instead."""
        def __get__(self):
            return self.query_alignment_length
        def __set__(self, v):
            self.query_alignment_length = v
    property mrnm:
        """deprecated, use :attr:`next_reference_id` instead."""
        def __get__(self):
            return self.next_reference_id
        def __set__(self, v):
            self.next_reference_id = v
    property mpos:
        """deprecated, use :attr:`next_reference_start` 
        instead."""
        def __get__(self):
            return self.next_reference_start
        def __set__(self, v):
            self.next_reference_start = v
    property rname:
        """deprecated, use :attr:`reference_id` instead."""
        def __get__(self):
            return self.reference_id
        def __set__(self, v):
            self.reference_id = v
    property isize:
        """deprecated, use :attr:`template_length` instead."""
        def __get__(self):
            return self.template_length
        def __set__(self, v):
            self.template_length = v
    property blocks:
        """deprecated, use :meth:`get_blocks()` instead."""
        def __get__(self):
            return self.get_blocks()
    property aligned_pairs:
        """deprecated, use :meth:`get_aligned_pairs()` instead."""
        def __get__(self):
            return self.get_aligned_pairs()
    property inferred_length:
        """deprecated, use :meth:`infer_query_length()` instead."""
        def __get__(self):
            return self.infer_query_length()
    property positions:
        """deprecated, use :meth:`get_reference_positions()` instead."""
        def __get__(self):
            return self.get_reference_positions()
    property tags:
        """deprecated, use :meth:`get_tags()` instead."""
        def __get__(self):
            return self.get_tags()
        def __set__(self, tags):
            self.set_tags(tags)
    def overlap(self):
        """deprecated, use :meth:`get_overlap()` instead."""
        return self.get_overlap()
    def opt(self, tag):
        """deprecated, use :meth:`get_tag()` instead."""
        return self.get_tag(tag)
    def setTag(self, tag, value, value_type=None, replace=True):
        """deprecated, use :meth:`set_tag()` instead."""
        return self.set_tag(tag, value, value_type, replace)


cdef class PileupColumn:
    '''A pileup of reads at a particular reference sequence position
    (:term:`column`). A pileup column contains all the reads that map
    to a certain target base.

    This class is a proxy for results returned by the samtools pileup
    engine.  If the underlying engine iterator advances, the results
    of this column will change.
    '''
    def __init__(self):
        raise TypeError("this class cannot be instantiated from Python")

    def __str__(self):
        return "\t".join(map(str,
                              (self.reference_id,
                               self.reference_pos,
                               self.nsegments))) +\
            "\n" +\
            "\n".join(map(str, self.pileups))

    def __dealloc__(self):
        free(self.buf.s)

    def set_min_base_quality(self, min_base_quality):
        """set the minimum base quality for this pileup column.
        """
        self.min_base_quality = min_base_quality

    def __len__(self):
        """return number of reads aligned to this column.

        see :meth:`get_num_aligned`
        """
        return self.get_num_aligned()

    property reference_id:
        '''the reference sequence number as defined in the header'''
        def __get__(self):
            return self.tid

    property reference_name:
        """:term:`reference` name (None if no AlignmentFile is associated)"""
        def __get__(self):
            if self.header is not None:
                return self.header.get_reference_name(self.tid)
            return None

    property nsegments:
        '''number of reads mapping to this column.

        Note that this number ignores the base quality filter.'''
        def __get__(self):
            return self.n_pu
        def __set__(self, n):
            self.n_pu = n

    property reference_pos:
        '''the position in the reference sequence (0-based).'''
        def __get__(self):
            return self.pos

    property pileups:
        '''list of reads (:class:`pysam.PileupRead`) aligned to this column'''
        def __get__(self):
            if self.plp == NULL or self.plp[0] == NULL:
                raise ValueError("PileupColumn accessed after iterator finished")

            cdef int x
            cdef const bam_pileup1_t * p = NULL
            pileups = []

            # warning: there could be problems if self.n and self.buf are
            # out of sync.
            for x from 0 <= x < self.n_pu:
                p = &(self.plp[0][x])
                if p == NULL:
                    raise ValueError(
                        "pileup buffer out of sync - most likely use of iterator "
                        "outside loop")
                if pileup_base_qual_skip(p, self.min_base_quality):
                    continue
                pileups.append(makePileupRead(p, self.header))
            return pileups

    ########################################################
    # Compatibility Accessors
    # Functions, properties for compatibility with pysam < 0.8
    ########################################################
    property pos:
        """deprecated, use :attr:`reference_pos` instead."""
        def __get__(self):
            return self.reference_pos
        def __set__(self, v):
            self.reference_pos = v

    property tid:
        """deprecated, use :attr:`reference_id` instead."""
        def __get__(self):
            return self.reference_id
        def __set__(self, v):
            self.reference_id = v

    property n:
        """deprecated, use :attr:`nsegments` instead."""
        def __get__(self):
            return self.nsegments
        def __set__(self, v):
            self.nsegments = v

    def get_num_aligned(self):
        """return number of aligned bases at pileup column position.

        This method applies a base quality filter and the number is
        equal to the size of :meth:`get_query_sequences`,
        :meth:`get_mapping_qualities`, etc.

        """
        cdef uint32_t x = 0
        cdef uint32_t c = 0
        cdef uint32_t cnt = 0
        cdef const bam_pileup1_t * p = NULL
        if self.plp == NULL or self.plp[0] == NULL:
            raise ValueError("PileupColumn accessed after iterator finished")

        for x from 0 <= x < self.n_pu:
            p = &(self.plp[0][x])
            if p == NULL:
                raise ValueError(
                    "pileup buffer out of sync - most likely use of iterator "
                    "outside loop")
            if pileup_base_qual_skip(p, self.min_base_quality):
                continue
            cnt += 1
        return cnt

    def get_query_sequences(self, bint mark_matches=False, bint mark_ends=False, bint add_indels=False):
        """query bases/sequences at pileup column position.

        Optionally, the bases/sequences can be annotated according to the samtools
        mpileup format. This is the format description from the samtools mpileup tool::

           Information on match, mismatch, indel, strand, mapping
           quality and start and end of a read are all encoded at the
           read base column. At this column, a dot stands for a match
           to the reference base on the forward strand, a comma for a
           match on the reverse strand, a '>' or '<' for a reference
           skip, `ACGTN' for a mismatch on the forward strand and
           `acgtn' for a mismatch on the reverse strand. A pattern
           `\\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion
           between this reference position and the next reference
           position. The length of the insertion is given by the
           integer in the pattern, followed by the inserted
           sequence. Similarly, a pattern `-[0-9]+[ACGTNacgtn]+'
           represents a deletion from the reference. The deleted bases
           will be presented as `*' in the following lines. Also at
           the read base column, a symbol `^' marks the start of a
           read. The ASCII of the character following `^' minus 33
           gives the mapping quality. A symbol `$' marks the end of a
           read segment

        To reproduce samtools mpileup format, set all of mark_matches,
        mark_ends and add_indels to True.

        Parameters
        ----------

        mark_matches: bool

          If True, output bases matching the reference as "." or ","
          for forward and reverse strand, respectively. This mark
          requires the reference sequence. If no reference is
          present, this option is ignored.

        mark_ends : bool

          If True, add markers "^" and "$" for read start and end, respectively.

        add_indels : bool

          If True, add bases for bases inserted into or skipped from the
          reference. The latter requires a reference sequence file to have
          been given, e.g. via `pileup(fastafile = ...)`. If no reference
          sequence is available, skipped bases are represented as 'N's.

        Returns
        -------

        a list of bases/sequences per read at pileup column position. : list

        """
        cdef uint32_t x = 0
        cdef uint32_t j = 0
        cdef uint32_t c = 0
        cdef uint8_t cc = 0
        cdef uint8_t rb = 0
        cdef kstring_t * buf = &self.buf
        cdef const bam_pileup1_t * p = NULL

        if self.plp == NULL or self.plp[0] == NULL:
            raise ValueError("PileupColumn accessed after iterator finished")

        buf.l = 0

        # todo: reference sequence to count matches/mismatches
        # todo: convert assertions to exceptions
        for x from 0 <= x < self.n_pu:
            p = &(self.plp[0][x])
            if p == NULL:
                raise ValueError(
                    "pileup buffer out of sync - most likely use of iterator "
                    "outside loop")
            if pileup_base_qual_skip(p, self.min_base_quality):
                continue
            # see samtools pileup_seq
            if mark_ends and p.is_head:
                kputc(b'^', buf)

                if p.b.core.qual > 93:
                    kputc(126, buf)
                else:
                    kputc(p.b.core.qual + 33, buf)
            if not p.is_del:
                if p.qpos < p.b.core.l_qseq:
                    cc = <uint8_t>seq_nt16_str[bam_seqi(bam_get_seq(p.b), p.qpos)]
                else:
                    cc = b'N'

                if mark_matches and self.reference_sequence != NULL:
                    rb = self.reference_sequence[self.reference_pos]
                    if seq_nt16_table[cc] == seq_nt16_table[rb]:
                        cc = b'='
                kputc(strand_mark_char(cc, p.b), buf)
            elif add_indels:
                if p.is_refskip:
                    if bam_is_rev(p.b):
                        kputc(b'<', buf)
                    else:
                        kputc(b'>', buf)
                else:
                    kputc(b'*', buf)
            if add_indels:
                if p.indel > 0:
                    kputc(b'+', buf)
                    kputw(p.indel, buf)
                    for j from 1 <= j <= p.indel:
                        cc = seq_nt16_str[bam_seqi(bam_get_seq(p.b), p.qpos + j)]
                        kputc(strand_mark_char(cc, p.b), buf)
                elif p.indel < 0:
                    kputc(b'-', buf)
                    kputw(-p.indel, buf)
                    for j from 1 <= j <= -p.indel:
                        # TODO: out-of-range check here?
                        if self.reference_sequence == NULL:
                            cc = b'N'
                        else:
                            cc = self.reference_sequence[self.reference_pos + j]
                        kputc(strand_mark_char(cc, p.b), buf)
            if mark_ends and p.is_tail:
                kputc(b'$', buf)

            kputc(b':', buf)

        if buf.l == 0:
            # could be zero if all qualities are too low
            return ""
        else:
            # quicker to ensemble all and split than to encode all separately.
            # ignore last ":"
            return force_str(PyBytes_FromStringAndSize(buf.s, buf.l-1)).split(":")

    def get_query_qualities(self):
        """query base quality scores at pileup column position.

        Returns
        -------

        a list of quality scores : list
        """
        cdef uint32_t x = 0
        cdef const bam_pileup1_t * p = NULL
        cdef uint32_t c = 0
        result = []
        for x from 0 <= x < self.n_pu:
            p = &(self.plp[0][x])
            if p == NULL:
                raise ValueError(
                    "pileup buffer out of sync - most likely use of iterator "
                    "outside loop")

            if p.qpos < p.b.core.l_qseq:
                c = bam_get_qual(p.b)[p.qpos]
            else:
                c = 0
            if c < self.min_base_quality:
                continue
            result.append(c)
        return result

    def get_mapping_qualities(self):
        """query mapping quality scores at pileup column position.

        Returns
        -------

        a list of quality scores : list
        """
        if self.plp == NULL or self.plp[0] == NULL:
            raise ValueError("PileupColumn accessed after iterator finished")

        cdef uint32_t x = 0
        cdef const bam_pileup1_t * p = NULL
        result = []
        for x from 0 <= x < self.n_pu:
            p = &(self.plp[0][x])
            if p == NULL:
                raise ValueError(
                    "pileup buffer out of sync - most likely use of iterator "
                    "outside loop")

            if pileup_base_qual_skip(p, self.min_base_quality):
                continue
            result.append(p.b.core.qual)
        return result

    def get_query_positions(self):
        """positions in read at pileup column position.

        Returns
        -------

        a list of read positions : list
        """
        if self.plp == NULL or self.plp[0] == NULL:
            raise ValueError("PileupColumn accessed after iterator finished")

        cdef uint32_t x = 0
        cdef const bam_pileup1_t * p = NULL
        result = []
        for x from 0 <= x < self.n_pu:
            p = &(self.plp[0][x])
            if p == NULL:
                raise ValueError(
                    "pileup buffer out of sync - most likely use of iterator "
                    "outside loop")

            if pileup_base_qual_skip(p, self.min_base_quality):
                continue
            result.append(p.qpos)
        return result

    def get_query_names(self):
        """query/read names aligned at pileup column position.

        Returns
        -------

        a list of query names at pileup column position. : list
        """
        if self.plp == NULL or self.plp[0] == NULL:
            raise ValueError("PileupColumn accessed after iterator finished")

        cdef uint32_t x = 0
        cdef const bam_pileup1_t * p = NULL
        result = []
        for x from 0 <= x < self.n_pu:
            p = &(self.plp[0][x])
            if p == NULL:
                raise ValueError(
                    "pileup buffer out of sync - most likely use of iterator "
                    "outside loop")

            if pileup_base_qual_skip(p, self.min_base_quality):
                continue
            result.append(charptr_to_str(pysam_bam_get_qname(p.b)))
        return result


cdef class PileupRead:
    '''Representation of a read aligned to a particular position in the
    reference sequence.

    '''

    def __init__(self):
        raise TypeError(
            "this class cannot be instantiated from Python")

    def __str__(self):
        return "\t".join(
            map(str,
                (self.alignment, self.query_position,
                 self.indel, self.level,
                 self.is_del, self.is_head,
                 self.is_tail, self.is_refskip)))

    property alignment:
        """a :class:`pysam.AlignedSegment` object of the aligned read"""
        def __get__(self):
            return self._alignment

    property query_position:
        """position of the read base at the pileup site, 0-based.
        None if :attr:`is_del` or :attr:`is_refskip` is set.

        """
        def __get__(self):
            if self.is_del or self.is_refskip:
                return None
            else:
                return self._qpos

    property query_position_or_next:
        """position of the read base at the pileup site, 0-based.

        If the current position is a deletion, returns the next
        aligned base.

        """
        def __get__(self):
            return self._qpos

    property indel:
        """indel length for the position following the current pileup site.

        This quantity peeks ahead to the next cigar operation in this
        alignment. If the next operation is an insertion, indel will
        be positive. If the next operation is a deletion, it will be
        negation. 0 if the next operation is not an indel.

        """
        def __get__(self):
            return self._indel

    property level:
        """the level of the read in the "viewer" mode. Note that this value
        is currently not computed."""
        def __get__(self):
            return self._level

    property is_del:
        """1 iff the base on the padded read is a deletion"""
        def __get__(self):
            return self._is_del

    property is_head:
        """1 iff the base on the padded read is the left-most base."""
        def __get__(self):
            return self._is_head

    property is_tail:
        """1 iff the base on the padded read is the right-most base."""
        def __get__(self):
            return self._is_tail

    property is_refskip:
        """1 iff the base on the padded read is part of CIGAR N op."""
        def __get__(self):
            return self._is_refskip



cpdef enum CIGAR_OPS:
    CMATCH = 0
    CINS = 1
    CDEL = 2
    CREF_SKIP = 3
    CSOFT_CLIP = 4
    CHARD_CLIP = 5
    CPAD = 6
    CEQUAL = 7
    CDIFF = 8
    CBACK = 9


cpdef enum SAM_FLAGS:
    # the read is paired in sequencing, no matter whether it is mapped in a pair
    FPAIRED = 1
    # the read is mapped in a proper pair
    FPROPER_PAIR = 2
    # the read itself is unmapped; conflictive with FPROPER_PAIR
    FUNMAP = 4
    # the mate is unmapped
    FMUNMAP = 8
    # the read is mapped to the reverse strand
    FREVERSE = 16
    # the mate is mapped to the reverse strand
    FMREVERSE = 32
    # this is read1
    FREAD1 = 64
    # this is read2
    FREAD2 = 128
    # not primary alignment
    FSECONDARY = 256
    # QC failure
    FQCFAIL = 512
    # optical or PCR duplicate
    FDUP = 1024
    # supplementary alignment
    FSUPPLEMENTARY = 2048


__all__ = [
    "AlignedSegment",
    "PileupColumn",
    "PileupRead",
    "CMATCH",
    "CINS",
    "CDEL",
    "CREF_SKIP",
    "CSOFT_CLIP",
    "CHARD_CLIP",
    "CPAD",
    "CEQUAL",
    "CDIFF",
    "CBACK",
    "FPAIRED",
    "FPROPER_PAIR",
    "FUNMAP",
    "FMUNMAP",
    "FREVERSE",
    "FMREVERSE",
    "FREAD1",
    "FREAD2",
    "FSECONDARY",
    "FQCFAIL",
    "FDUP",
    "FSUPPLEMENTARY",
    "KEY_NAMES"]
