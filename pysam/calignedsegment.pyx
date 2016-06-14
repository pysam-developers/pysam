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
import ctypes
import struct

cimport cython
from cpython cimport array as c_array
from cpython.version cimport PY_MAJOR_VERSION
from cpython cimport PyErr_SetString, PyBytes_FromStringAndSize
from libc.string cimport strchr
from cpython cimport array as c_array

from pysam.cutils cimport force_bytes, force_str, \
    charptr_to_str, charptr_to_bytes
from pysam.cutils cimport qualities_to_qualitystring, qualitystring_to_array, \
    array_to_qualitystring

# Constants for binary tag conversion
cdef char * htslib_types = 'cCsSiIf'
cdef char * parray_types = 'bBhHiIf'

# translation tables

# cigar code to character and vice versa
cdef char* CODE2CIGAR= "MIDNSHP=XB"
cdef int NCIGAR_CODES = 10

if PY_MAJOR_VERSION >= 3:
    CIGAR2CODE = dict([y, x] for x, y in enumerate(CODE2CIGAR))
else:
    CIGAR2CODE = dict([ord(y), x] for x, y in enumerate(CODE2CIGAR))

CIGAR_REGEX = re.compile("(\d+)([MIDNSHP=XB])")

#####################################################################
# typecode guessing
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


cdef inline uint8_t get_value_code(value, value_type=None):
    '''guess type code for a *value*. If *value_type* is None,
    the type code will be inferred based on the Python type of
    *value*'''
    cdef uint8_t  typecode
    cdef char * _char_type

    if value_type is None:
        if isinstance(value, int):
            typecode = 'i'
        elif isinstance(value, float):
            typecode = 'd'
        elif isinstance(value, str):
            typecode = 'Z'
        elif isinstance(value, bytes):
            typecode = 'Z'
        elif isinstance(value, array.array) or \
                isinstance(value, list) or \
                isinstance(value, tuple):
            typecode = 'B'
        else:
            return 0
    else:
        if value_type not in 'Zidf':
            return 0
        value_type = force_bytes(value_type)
        _char_type = value_type
        typecode = (<uint8_t*>_char_type)[0]

    return typecode


cdef inline bytes getTypecode(value, maximum_value=None):
    '''returns the value typecode of a value.

    If max is specified, the approprite type is
    returned for a range where value is the minimum.
    '''

    if maximum_value is None:
        maximum_value = value

    cdef bytes valuetype

    t = type(value)

    if t is float:
        valuetype = b'f'
    elif t is int:
        # signed ints
        if value < 0:
            if value >= -128 and maximum_value < 128:
                valuetype = b'c'
            elif value >= -32768 and maximum_value < 32768:
                valuetype = b's'
            elif value < -2147483648 or maximum_value >= 2147483648:
                raise ValueError(
                    "at least one signed integer out of range of "
                    "BAM/SAM specification")
            else:
                valuetype = b'i'
        # unsigned ints
        else:
            if maximum_value < 256:
                valuetype = b'C'
            elif maximum_value < 65536:
                valuetype = b'S'
            elif maximum_value >= 4294967296:
                raise ValueError(
                    "at least one integer out of range of BAM/SAM specification")
            else:
                valuetype = b'I'
    else:
        # Note: hex strings (H) are not supported yet
        if t is not bytes:
            value = value.encode('ascii')
        if len(value) == 1:
            valuetype = b'A'
        else:
            valuetype = b'Z'

    return valuetype


cdef inline packTags(tags):
    """pack a list of tags. Each tag is a tuple of (tag, tuple).

    Values are packed into the most space efficient data structure
    possible unless the tag contains a third field with the typecode.

    Returns a format string and the associated list of arguments
    to be used in a call to struct.pack_into.
    """
    fmts, args = ["<"], []
    
    cdef char array_typecode

    datatype2format = {
        b'c': ('b', 1),
        b'C': ('B', 1),
        b's': ('h', 2),
        b'S': ('H', 2),
        b'i': ('i', 4),
        b'I': ('I', 4),
        b'f': ('f', 4),
        b'A': ('c', 1)}

    for tag in tags:

        if len(tag) == 2:
            pytag, value = tag
            valuetype = None
        elif len(tag) == 3:
            pytag, value, valuetype = tag
        else:
            raise ValueError("malformatted tag: %s" % str(tag))

        pytag = force_bytes(pytag)
        valuetype = force_bytes(valuetype)
        t = type(value)

        if t is tuple or t is list:
            # binary tags from tuples or lists
            if valuetype is None:
                # automatically determine value type - first value
                # determines type. If there is a mix of types, the
                # result is undefined.
                valuetype = getTypecode(min(value), max(value))

            if valuetype not in datatype2format:
                raise ValueError("invalid value type '%s'" % valuetype)

            datafmt = "2sccI%i%s" % (len(value), datatype2format[valuetype][0])
            args.extend([pytag[:2],
                         b"B",
                         valuetype,
                         len(value)] + list(value))

        elif isinstance(value, array.array):
            # binary tags from arrays
            if valuetype is None:
                array_typecode = map_typecode_python_to_htslib(ord(value.typecode))

                if array_typecode == 0:
                    raise ValueError("unsupported type code '{}'"
                                     .format(value.typecode))

                valuetype = force_bytes(chr(array_typecode))
                    
            if valuetype not in datatype2format:
                raise ValueError("invalid value type '%s' (%s)" %
                                 (valuetype, type(valuetype)))

            # use array.tostring() to retrieve byte representation and
            # save as bytes
            datafmt = "2sccI%is" % (len(value) * datatype2format[valuetype][1])
            args.extend([pytag[:2],
                         b"B",
                         valuetype,
                         len(value),
                         force_bytes(value.tostring())])

        else:
            if valuetype is None:
                valuetype = getTypecode(value)

            if valuetype in b"AZ":
                value = force_bytes(value)

            if valuetype == b"Z":
                datafmt = "2sc%is" % (len(value)+1)
            else:
                datafmt = "2sc%s" % datatype2format[valuetype][0]

            args.extend([pytag[:2],
                         valuetype,
                         value])

        fmts.append(datafmt)

    return "".join(fmts), args


cdef inline int32_t calculateQueryLength(bam1_t * src):
    """return query length computed from CIGAR alignment.

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

        if op == BAM_CMATCH or op == BAM_CINS or \
           op == BAM_CSOFT_CLIP or \
           op == BAM_CEQUAL or op == BAM_CDIFF:
            qpos += cigar_p[k] >> BAM_CIGAR_SHIFT

    return qpos


cdef inline int32_t getQueryStart(bam1_t *src) except -1:
    cdef uint32_t * cigar_p
    cdef uint32_t k, op
    cdef uint32_t start_offset = 0

    if pysam_get_n_cigar(src):
        cigar_p = pysam_bam_get_cigar(src);
        for k from 0 <= k < pysam_get_n_cigar(src):
            op = cigar_p[k] & BAM_CIGAR_MASK
            if op == BAM_CHARD_CLIP:
                if start_offset != 0 and start_offset != src.core.l_qseq:
                    PyErr_SetString(ValueError, 'Invalid clipping in CIGAR string')
                    return -1
            elif op == BAM_CSOFT_CLIP:
                start_offset += cigar_p[k] >> BAM_CIGAR_SHIFT
            else:
                break

    return start_offset


cdef inline int32_t getQueryEnd(bam1_t *src) except -1:
    cdef uint32_t * cigar_p
    cdef uint32_t k, op
    cdef uint32_t end_offset = src.core.l_qseq

    # if there is no sequence, compute length from cigar string
    if end_offset == 0:
        end_offset = calculateQueryLength(src)

    # walk backwards in cigar string
    if pysam_get_n_cigar(src) > 1:
        cigar_p = pysam_bam_get_cigar(src);
        for k from pysam_get_n_cigar(src) > k >= 1:
            op = cigar_p[k] & BAM_CIGAR_MASK
            if op == BAM_CHARD_CLIP:
                if end_offset != 0 and end_offset != src.core.l_qseq:
                    PyErr_SetString(ValueError,
                                    'Invalid clipping in CIGAR string')
                    return -1
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
        s[k-start] = seq_nt16_str[p[k/2] >> 4 * (1 - k%2) & 0xf]

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
## private factory methods
cdef class AlignedSegment
cdef makeAlignedSegment(bam1_t * src, AlignmentFile alignment_file):
    '''return an AlignedSegment object constructed from `src`'''
    # note that the following does not call __init__
    cdef AlignedSegment dest = AlignedSegment.__new__(AlignedSegment)
    dest._delegate = bam_dup1(src)
    dest._alignment_file = alignment_file
    return dest


cdef class PileupColumn
cdef makePileupColumn(bam_pileup1_t ** plp, int tid, int pos,
                      int n_pu, AlignmentFile alignment_file):
    '''return a PileupColumn object constructed from pileup in `plp` and
    setting additional attributes.

    '''
    # note that the following does not call __init__
    cdef PileupColumn dest = PileupColumn.__new__(PileupColumn)
    dest._alignment_file = alignment_file
    dest.plp = plp
    dest.tid = tid
    dest.pos = pos
    dest.n_pu = n_pu
    return dest

cdef class PileupRead
cdef inline makePileupRead(bam_pileup1_t * src, AlignmentFile alignment_file):
    '''return a PileupRead object construted from a bam_pileup1_t * object.'''
    cdef PileupRead dest = PileupRead.__new__(PileupRead)
    dest._alignment = makeAlignedSegment(src.b, alignment_file)
    dest._qpos = src.qpos
    dest._indel = src.indel
    dest._level = src.level
    dest._is_del = src.is_del
    dest._is_head = src.is_head
    dest._is_tail = src.is_tail
    dest._is_refskip = src.is_refskip
    return dest


cdef inline uint32_t get_alignment_length(bam1_t * src):
    cdef int k = 0
    cdef uint32_t l = 0
    if src == NULL:
        return 0
    cdef uint32_t * cigar_p = bam_get_cigar(src)
    if cigar_p == NULL:
        return 0
    cdef int op
    cdef int n = pysam_get_n_cigar(src)
    for k from 0 <= k < n:
        op = cigar_p[k] & BAM_CIGAR_MASK
        if op == BAM_CSOFT_CLIP or op == BAM_CHARD_CLIP:
            continue
        l += cigar_p[k] >> BAM_CIGAR_SHIFT
    return l


# TODO: avoid string copying for getSequenceInRange, reconstituneSequenceFromMD, ...
cdef inline bytes build_alignment_sequence(bam1_t * src):
    """return expanded sequence from MD tag.

    The sequence includes substitutions and both insertions in the
    reference as well as deletions to the reference sequence. Combine
    with the cigar string to reconstitute the query or the reference
    sequence.

    Positions corresponding to `N` (skipped region from the reference)
    in the CIGAR string will not appear in the returned sequence. The
    MD should correspondingly not contain these. Thus proper tags are::
    
       Deletion from the reference:   cigar=5M1D5M    MD=5^C5
       Skipped region from reference: cigar=5M1N5M    MD=10

    Returns
    -------

    None, if no MD tag is present.

    """
    if src == NULL:
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
            "could not allocated sequence of length %i" % max_len)

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
                s[s_idx] = '-'
                s_idx += 1
        elif op == BAM_CREF_SKIP:
            pass
        elif op == BAM_CINS:
            for i from 0 <= i < l:
                # encode insertions into reference as lowercase
                s[s_idx] = read_sequence[r_idx] + 32
                r_idx += 1
                s_idx += 1
        elif op == BAM_CSOFT_CLIP:
            pass
        elif op == BAM_CHARD_CLIP:
            pass # advances neither
        elif op == BAM_CPAD:
            raise NotImplementedError(
                "Padding (BAM_CPAD, 6) is currently not supported. "
                "Please implement. Sorry about that.")

    cdef uint8_t * md_tag_ptr = bam_aux_get(src, "MD")    
    if md_tag_ptr == NULL:
        seq = PyBytes_FromStringAndSize(s, s_idx)
        free(s)
        return seq

    cdef char * md_tag = <char*>bam_aux2Z(md_tag_ptr)
    cdef int md_idx = 0
    s_idx = 0

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
                while s[s_idx] >= 'a':
                    s_idx += 1
                s_idx += 1
            while s[s_idx] >= 'a':
                s_idx += 1

            r_idx += nmatches
            nmatches = 0
            if md_tag[md_idx] == '^':
                md_idx += 1
                while md_tag[md_idx] >= 65 and md_tag[md_idx] <= 90:
                    assert s[s_idx] == '-'
                    s[s_idx] = md_tag[md_idx]
                    s_idx += 1
                    md_idx += 1
            else:
                # save mismatch and change to lower case
                s[s_idx] = md_tag[md_idx] + 32
                s_idx += 1
                r_idx += 1
                md_idx += 1

    # save matches up to this point, skipping insertions
    for x from 0 <= x < nmatches:
        while s[s_idx] >= 'a':
            s_idx += 1
        s_idx += 1
    while s[s_idx] >= 'a':
        s_idx += 1

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
    '''

    # Now only called when instances are created from Python
    def __init__(self):
        # see bam_init1
        self._delegate = <bam1_t*>calloc(1, sizeof(bam1_t))
        # allocate some memory. If size is 0, calloc does not return a
        # pointer that can be passed to free() so allocate 40 bytes
        # for a new read
        self._delegate.m_data = 40
        self._delegate.data = <uint8_t *>calloc(
            self._delegate.m_data, 1)
        self._delegate.l_data = 0

        # caching for selected fields
        self.cache_query_qualities = None
        self.cache_query_alignment_qualities = None
        self.cache_query_sequence = None
        self.cache_query_alignment_sequence = None

    def __dealloc__(self):
        bam_destroy1(self._delegate)

    def __str__(self):
        """return string representation of alignment.

        The representation is an approximate :term:`SAM` format, because
        an aligned read might not be associated with a :term:`AlignmentFile`.
        As a result :term:`tid` is shown instead of the reference name.
        Similarly, the tags field is returned in its parsed state.

        To get a valid SAM record, use :meth:`tostring`.
        """
        # sam-parsing is done in sam.c/bam_format1_core which
        # requires a valid header.
        return "\t".join(map(str, (self.query_name,
                                   self.flag,
                                   self.reference_id,
                                   self.reference_start,
                                   self.mapping_quality,
                                   self.cigarstring,
                                   self.next_reference_id,
                                   self.next_reference_start,
                                   self.query_alignment_length,
                                   self.query_sequence,
                                   self.query_qualities,
                                   self.tags)))

    def __copy__(self):
        return makeAlignedSegment(self._delegate, self._alignment_file)

    def __deepcopy__(self, memo):
        return makeAlignedSegment(self._delegate, self._alignment_file)

    def compare(self, AlignedSegment other):
        '''return -1,0,1, if contents in this are binary
        <,=,> to *other*

        '''

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
        cdef bam1_t * src
        src = self._delegate
        # shift and xor values in the core structure
        # make sure tid and mtid are shifted by different amounts
        # should variable length data be included?
        cdef uint32_t hash_value = src.core.tid << 24 ^ \
            src.core.pos << 16 ^ \
            src.core.qual << 8 ^ \
            src.core.flag ^ \
            src.core.isize << 24 ^ \
            src.core.mtid << 16 ^ \
            src.core.mpos << 8

        return hash_value

    cpdef tostring(self, AlignmentFile_t htsfile):
        """returns a string representation of the aligned segment.

        The output format is valid SAM format.

        Parameters
        ----------

        htsfile -- AlignmentFile object to map numerical
                   identifers to chromosome names.
        """

        cdef kstring_t line
        line.l = line.m = 0
        line.s = NULL

        if sam_format1(htsfile.header, self._delegate, &line) < 0:
            if line.m:
                free(line.s)
            raise ValueError('sam_format failed')

        ret = force_str(line.s[:line.l])
        
        if line.m:
            free(line.s)

        return ret

    ########################################################
    ## Basic attributes in order of appearance in SAM format
    property query_name:
        """the query template name (None if not present)"""
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            if pysam_get_l_qname(src) == 0:
                return None
            return charptr_to_str(<char *>pysam_bam_get_qname(src))

        def __set__(self, qname):
            if qname is None or len(qname) == 0:
                return
            qname = force_bytes(qname)
            cdef bam1_t * src
            cdef int l
            cdef char * p

            src = self._delegate
            p = pysam_bam_get_qname(src)

            # the qname is \0 terminated
            l = len(qname) + 1
            pysam_bam_update(src,
                             pysam_get_l_qname(src),
                             l,
                             <uint8_t*>p)


            pysam_set_l_qname(src, l)

            # re-acquire pointer to location in memory
            # as it might have moved
            p = pysam_bam_get_qname(src)

            strncpy(p, qname, l)

    property flag:
        """properties flag"""
        def __get__(self):
            return pysam_get_flag(self._delegate)
        def __set__(self, flag):
            pysam_set_flag(self._delegate, flag)

    property reference_name:
        """:term:`reference` name (None if no AlignmentFile is associated)"""
        def __get__(self):
            if self._alignment_file is not None:
                return self._alignment_file.getrname(self._delegate.core.tid)
            return None

    property reference_id:
        """:term:`reference` ID

        .. note::

            This field contains the index of the reference sequence in
            the sequence dictionary. To obtain the name of the
            reference sequence, use
            :meth:`pysam.AlignmentFile.getrname()`

        """
        def __get__(self): return self._delegate.core.tid
        def __set__(self, tid): self._delegate.core.tid = tid

    property reference_start:
        """0-based leftmost coordinate"""
        def __get__(self): return self._delegate.core.pos
        def __set__(self, pos):
            ## setting the position requires updating the "bin" attribute
            cdef bam1_t * src
            src = self._delegate
            src.core.pos = pos
            if pysam_get_n_cigar(src):
                pysam_set_bin(src,
                              hts_reg2bin(
                                  src.core.pos,
                                  bam_endpos(src),
                                  14,
                                  5))
            else:
                pysam_set_bin(src,
                              hts_reg2bin(
                                  src.core.pos,
                                  src.core.pos + 1,
                                  14,
                                  5))

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
        def __get__(self): return self._delegate.core.mtid
        def __set__(self, mtid):
            self._delegate.core.mtid = mtid

    property next_reference_name:
        """:term:`reference` name of the mate/next read (None if no
        AlignmentFile is associated)"""
        def __get__(self):
            if self._alignment_file is not None:
                return self._alignment_file.getrname(self._delegate.core.mtid)
            return None

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
        :meth:`pysam.AlignmentFile.infer_query_length.`.

        The length includes soft-clipped bases and is equal to
        ``len(query_sequence)``.

        This property is read-only but can be set by providing a
        sequence.

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

        Note that assigning to seq will invalidate any quality scores.
        Thus, to in-place edit the sequence and quality scores, copies of
        the quality scores need to be taken. Consider trimming for example::

           q = read.query_qualities
           read.query_squence = read.query_sequence[5:10]
           read.query_qualities = q[5:10]

        The sequence is returned as it is stored in the BAM file. Some mappers
        might have stored a reverse complement of the original read
        sequence.
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
            nbytes_new = (l + 1) / 2 + l
            nbytes_old = (src.core.l_qseq + 1) / 2 + src.core.l_qseq

            # acquire pointer to location in memory
            p = pysam_bam_get_seq(src)
            src.core.l_qseq = l

            # change length of data field
            pysam_bam_update(src,
                             nbytes_old,
                             nbytes_new,
                             p)

            if l > 0:
                # re-acquire pointer to location in memory
                # as it might have moved
                p = pysam_bam_get_seq(src)
                for k from 0 <= k < nbytes_new:
                    p[k] = 0
                # convert to C string
                s = seq
                for k from 0 <= k < l:
                    p[k/2] |= seq_nt16_table[<unsigned char>s[k]] << 4 * (1 - k % 2)

                # erase qualities
                p = pysam_bam_get_qual(src)
                p[0] = 0xff

            self.cache_query_sequence = force_str(seq)

            # clear cached values for quality values
            self.cache_query_qualities = None
            self.cache_query_alignment_qualities = None

    property query_qualities:
        """read sequence base qualities, including :term:`soft
        clipped` bases (None if not present).

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
                if src.core.l_qseq != 0:
                    p[0] = 0xff
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
            return pysam_get_bin(self._delegate)
        def __set__(self, bin):
            pysam_set_bin(self._delegate, bin)


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
    property mate_is_unmapped:
        """true if the mate is unmapped"""
        def __get__(self):
            return (self.flag & BAM_FMUNMAP) != 0
        def __set__(self,val):
            pysam_update_flag(self._delegate, val, BAM_FMUNMAP)
    property is_reverse:
        """true if read is mapped to reverse strand"""
        def __get__(self):
            return (self.flag & BAM_FREVERSE) != 0
        def __set__(self,val):
            pysam_update_flag(self._delegate, val, BAM_FREVERSE)
    property mate_is_reverse:
        """true is read is mapped to reverse strand"""
        def __get__(self):
            return (self.flag & BAM_FMREVERSE) != 0
        def __set__(self,val):
            pysam_update_flag(self._delegate, val, BAM_FMREVERSE)
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

        This is equal to `aend - pos`. Returns None if not available.'''
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            if (self.flag & BAM_FUNMAP) or pysam_get_n_cigar(src) == 0:
                return None
            return bam_endpos(src) - \
                self._delegate.core.pos

    property query_alignment_sequence:
        """aligned portion of the read.

        This is a substring of :attr:`seq` that excludes flanking
        bases that were :term:`soft clipped` (None if not present). It
        is equal to ``seq[qstart:qend]``.

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
        are the quality values that correspond to :attr:`query`, that
        is, they exclude qualities of :term:`soft clipped` bases. This
        is equal to ``qual[qstart:qend]``.

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

        This the index of the first base in :attr:`seq` that is not
        soft-clipped.

        """
        def __get__(self):
            return getQueryStart(self._delegate)

    property query_alignment_end:
        """end index of the aligned query portion of the sequence (0-based,
        exclusive)"""
        def __get__(self):
            return getQueryEnd(self._delegate)

    property query_alignment_length:
        """length of the aligned query sequence.

        This is equal to :attr:`qend` - :attr:`qstart`"""
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
        cdef uint32_t k, i, pos
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
            elif op == BAM_CMATCH:
                for i from pos <= i < pos + l:
                    result.append(i)
                pos += l
            elif op == BAM_CDEL or op == BAM_CREF_SKIP:
                pos += l

        return result

    def infer_query_length(self, always=True):
        """inferred read length from CIGAR string.

        If *always* is set to True, the read length
        will be always inferred. If set to False, the length
        of the read sequence will be returned if it is
        available.

        Returns None if CIGAR string is not present.
        """

        cdef uint32_t * cigar_p
        cdef bam1_t * src

        src = self._delegate

        if not always and src.core.l_qseq:
            return src.core.l_qseq

        return calculateQueryLength(src)

    def get_reference_sequence(self):
        """return the reference sequence.

        This method requires the MD tag to be set.
        """
        cdef uint32_t k, i
        cdef int op
        cdef bam1_t * src = self._delegate
        ref_seq = force_str(build_alignment_sequence(src))
        if ref_seq is None:
            raise ValueError("MD tag not present")

        cdef uint32_t * cigar_p = pysam_bam_get_cigar(src)
        cdef uint32_t r_idx = 0
        result = []
        for k from 0 <= k < pysam_get_n_cigar(src):
            op = cigar_p[k] & BAM_CIGAR_MASK
            l = cigar_p[k] >> BAM_CIGAR_SHIFT
            if op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
                for i from 0 <= i < l:
                    result.append(ref_seq[r_idx])
                    r_idx += 1
            elif op == BAM_CDEL:
                for i from 0 <= i < l:
                    result.append(ref_seq[r_idx])
                    r_idx += 1
            elif op == BAM_CREF_SKIP:
                pass
            elif op == BAM_CINS:
                r_idx += l
            elif op == BAM_CSOFT_CLIP:
                pass
            elif op == BAM_CHARD_CLIP:
                pass # advances neither
            elif op == BAM_CPAD:
                raise NotImplementedError(
                    "Padding (BAM_CPAD, 6) is currently not supported. "
                    "Please implement. Sorry about that.")

        return "".join(result)

    def get_aligned_pairs(self, matches_only=False, with_seq=False):
        """a list of aligned read (query) and reference positions.

        For inserts, deletions, skipping either query or reference
        position may be None.

        Padding is currently not supported and leads to an exception.

        Parameters
        ----------

        matches_only : bool
          If True, only matched bases are returned - no None on either
          side.
        with_seq : bool
          If True, return a third element in the tuple containing the
          reference sequence. Substitutions are lower-case. This option
          requires an MD tag to be present.

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
            ref_seq = force_str(self.get_reference_sequence())
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

            elif op == BAM_CINS or op == BAM_CSOFT_CLIP:
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

            elif op == BAM_CPAD:
                raise NotImplementedError(
                    "Padding (BAM_CPAD, 6) is currently not supported. "
                    "Please implement. Sorry about that.")

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
            if op == BAM_CMATCH:
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

            if op == BAM_CMATCH:
                o = min( pos + l, end) - max( pos, start )
                if o > 0: overlap += o

            if op == BAM_CMATCH or op == BAM_CDEL or op == BAM_CREF_SKIP:
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
        |NM   |NM tag        |9    |
        +-----+--------------+-----+

        If no cigar string is present, empty arrays will be returned.

        Parameters
        ----------

        Returns
        -------

        arrays : two arrays. The first contains the nucleotide counts within
           each cigar operation, the second contains the number of blocks for
           each cigar operation.

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
            cdef int k

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
            cdef int k, ncigar

            k = 0

            src = self._delegate

            # get location of cigar string
            p = pysam_bam_get_cigar(src)

            # empty values for cigar string
            if values is None:
                values = []

            ncigar = len(values)
            # create space for cigar data within src.data
            pysam_bam_update(src,
                             pysam_get_n_cigar(src) * 4,
                             ncigar * 4,
                             <uint8_t*>p)

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
            pysam_set_bin(src,
                          hts_reg2bin(
                              src.core.pos,
                              bam_endpos(src),
                              14,
                              5))


    cpdef set_tag(self,
                  tag,
                  value,
                  value_type=None,
                  replace=True):
        """sets a particular field *tag* to *value* in the optional alignment
        section.

        *value_type* describes the type of *value* that is to entered
        into the alignment record.. It can be set explicitly to one
        of the valid one-letter type codes. If unset, an appropriate
        type will be chosen automatically.

        An existing value of the same *tag* will be overwritten unless
        replace is set to False. This is usually not recommened as a
        tag may only appear once in the optional alignment section.

        If *value* is None, the tag will be deleted.
        """

        cdef int value_size
        cdef uint8_t * value_ptr
        cdef uint8_t *existing_ptr
        cdef uint8_t typecode
        cdef float float_value
        cdef double double_value
        cdef int32_t int_value
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

        typecode = get_value_code(value, value_type)
        if typecode == 0:
            raise ValueError("can't guess type or invalid type code specified")

        # Not Endian-safe, but then again neither is samtools!
        if typecode == 'Z':
            value = force_bytes(value)
            value_ptr = <uint8_t*><char*>value
            value_size = len(value)+1
        elif typecode == 'i':
            int_value = value
            value_ptr = <uint8_t*>&int_value
            value_size = sizeof(int32_t)
        elif typecode == 'd':
            double_value = value
            value_ptr = <uint8_t*>&double_value
            value_size = sizeof(double)
        elif typecode == 'f':
            float_value  = value
            value_ptr = <uint8_t*>&float_value
            value_size = sizeof(float)
        elif typecode == 'B':
            # the following goes through python, needs to be cleaned up
            # pack array using struct
            if value_type is None:
                fmt, args = packTags([(tag, value)])
            else:
                fmt, args = packTags([(tag, value, value_type)])

            # remove tag and type code as set by bam_aux_append
            # first four chars of format (<2sc)
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
            raise ValueError('unsupported value_type in set_option')

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

        Parameters
        ----------

        tag :
            data tag.

        with_value_type : Optional[bool]
            if set to True, the return value is a tuple of (tag value, type code).
            (default False)

        Returns
        -------

        A python object with the value of the `tag`. The type of the
        object depends on the data type in the data record.

        Raises
        ------

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

        if auxtype == 'c' or auxtype == 'C' or auxtype == 's' or auxtype == 'S':
            value = <int>bam_aux2i(v)
        elif auxtype == 'i' or auxtype == 'I':
            value = <int32_t>bam_aux2i(v)
        elif auxtype == 'f' or auxtype == 'F':
            value = <float>bam_aux2f(v)
        elif auxtype == 'd' or auxtype == 'D':
            value = <double>bam_aux2f(v)
        elif auxtype == 'A':
            # there might a more efficient way
            # to convert a char into a string
            value = '%c' % <char>bam_aux2A(v)
        elif auxtype == 'Z':
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
        """the fields in the optional aligment section.

        Returns a list of all fields in the optional
        alignment section. Values are converted to appropriate python
        values. For example:

        [(NM, 2), (RG, "GJP00TM04")]

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
            if auxtype in ('c', 'C'):
                value = <int>bam_aux2i(s)
                s += 1
            elif auxtype in ('s', 'S'):
                value = <int>bam_aux2i(s)
                s += 2
            elif auxtype in ('i', 'I'):
                value = <int32_t>bam_aux2i(s)
                s += 4
            elif auxtype == 'f':
                value = <float>bam_aux2f(s)
                s += 4
            elif auxtype == 'd':
                value = <double>bam_aux2f(s)
                s += 8
            elif auxtype == 'A':
                value = "%c" % <char>bam_aux2A(s)
                s += 1
            elif auxtype in ('Z', 'H'):
                value = charptr_to_str(<char*>bam_aux2Z(s))
                # +1 for NULL terminated string
                s += len(value) + 1
            elif auxtype == 'B':
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
        """sets the fields in the optional alignmest section with
        a list of (tag, value) tuples.

        The :term:`value type` of the values is determined from the
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
            fmt, args = packTags(tags)
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
        pysam_bam_update(src,
                         old_size,
                         new_size,
                         pysam_bam_get_aux(src))

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
        """deprecated, use query_name instead"""
        def __get__(self): return self.query_name
        def __set__(self, v): self.query_name = v
    property tid:
        """deprecated, use reference_id instead"""
        def __get__(self): return self.reference_id
        def __set__(self, v): self.reference_id = v
    property pos:
        """deprecated, use reference_start instead"""
        def __get__(self): return self.reference_start
        def __set__(self, v): self.reference_start = v
    property mapq:
        """deprecated, use mapping_quality instead"""
        def __get__(self): return self.mapping_quality
        def __set__(self, v): self.mapping_quality = v
    property rnext:
        """deprecated, use next_reference_id instead"""
        def __get__(self): return self.next_reference_id
        def __set__(self, v): self.next_reference_id = v
    property pnext:
        """deprecated, use next_reference_start instead"""
        def __get__(self):
            return self.next_reference_start
        def __set__(self, v):
            self.next_reference_start = v
    property cigar:
        """deprecated, use cigartuples instead"""
        def __get__(self):
            r = self.cigartuples
            if r is None:
                r = []
            return r
        def __set__(self, v): self.cigartuples = v
    property tlen:
        """deprecated, use template_length instead"""
        def __get__(self):
            return self.template_length
        def __set__(self, v):
            self.template_length = v
    property seq:
        """deprecated, use query_sequence instead"""
        def __get__(self):
            return self.query_sequence
        def __set__(self, v):
            self.query_sequence = v
    property qual:
        """deprecated, query_qualities instead"""
        def __get__(self):
            return array_to_qualitystring(self.query_qualities)
        def __set__(self, v):
            self.query_qualities = qualitystring_to_array(v)
    property alen:
        """deprecated, reference_length instead"""
        def __get__(self):
            return self.reference_length
        def __set__(self, v):
            self.reference_length = v
    property aend:
        """deprecated, reference_end instead"""
        def __get__(self):
            return self.reference_end
        def __set__(self, v):
            self.reference_end = v
    property rlen:
        """deprecated, query_length instead"""
        def __get__(self):
            return self.query_length
        def __set__(self, v):
            self.query_length = v
    property query:
        """deprecated, query_alignment_sequence instead"""
        def __get__(self):
            return self.query_alignment_sequence
        def __set__(self, v):
            self.query_alignment_sequence = v
    property qqual:
        """deprecated, query_alignment_qualities instead"""
        def __get__(self):
            return array_to_qualitystring(self.query_alignment_qualities)
        def __set__(self, v):
            self.query_alignment_qualities = qualitystring_to_array(v)
    property qstart:
        """deprecated, use query_alignment_start instead"""
        def __get__(self):
            return self.query_alignment_start
        def __set__(self, v):
            self.query_alignment_start = v
    property qend:
        """deprecated, use query_alignment_end instead"""
        def __get__(self):
            return self.query_alignment_end
        def __set__(self, v):
            self.query_alignment_end = v
    property qlen:
        """deprecated, use query_alignment_length instead"""
        def __get__(self):
            return self.query_alignment_length
        def __set__(self, v):
            self.query_alignment_length = v
    property mrnm:
        """deprecated, use next_reference_id instead"""
        def __get__(self):
            return self.next_reference_id
        def __set__(self, v):
            self.next_reference_id = v
    property mpos:
        """deprecated, use next_reference_start instead"""
        def __get__(self):
            return self.next_reference_start
        def __set__(self, v):
            self.next_reference_start = v
    property rname:
        """deprecated, use reference_id instead"""
        def __get__(self):
            return self.reference_id
        def __set__(self, v):
            self.reference_id = v
    property isize:
        """deprecated, use template_length instead"""
        def __get__(self):
            return self.template_length
        def __set__(self, v):
            self.template_length = v
    property blocks:
        """deprecated, use get_blocks() instead"""
        def __get__(self):
            return self.get_blocks()
    property aligned_pairs:
        """deprecated, use get_aligned_pairs() instead"""
        def __get__(self):
            return self.get_aligned_pairs()
    property inferred_length:
        """deprecated, use infer_query_length() instead"""
        def __get__(self):
            return self.infer_query_length()
    property positions:
        """deprecated, use get_reference_positions() instead"""
        def __get__(self):
            return self.get_reference_positions()
    property tags:
        """deprecated, use get_tags() instead"""
        def __get__(self):
            return self.get_tags()
        def __set__(self, tags):
            self.set_tags(tags)
    def overlap(self):
        """deprecated, use get_overlap() instead"""
        return self.get_overlap()
    def opt(self, tag):
        """deprecated, use get_tag() instead"""
        return self.get_tag(tag)
    def setTag(self, tag, value, value_type=None, replace=True):
        """deprecated, use set_tag() instead"""
        return self.set_tag(tag, value, value_type, replace)


cdef class PileupColumn:
    '''A pileup of reads at a particular reference sequence postion
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

    property reference_id:
        '''the reference sequence number as defined in the header'''
        def __get__(self):
            return self.tid

    property reference_name:
        """:term:`reference` name (None if no AlignmentFile is associated)"""
        def __get__(self):
            if self._alignment_file is not None:
                return self._alignment_file.getrname(self.tid)
            return None

    property nsegments:
        '''number of reads mapping to this column.'''
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
            cdef int x
            pileups = []

            if self.plp == NULL or self.plp[0] == NULL:
                raise ValueError("PileupColumn accessed after iterator finished")

            # warning: there could be problems if self.n and self.buf are
            # out of sync.
            for x from 0 <= x < self.n_pu:
                pileups.append(makePileupRead(&(self.plp[0][x]),
                                              self._alignment_file))
            return pileups

    ########################################################
    # Compatibility Accessors
    # Functions, properties for compatibility with pysam < 0.8
    ########################################################
    property pos:
        def __get__(self):
            return self.reference_pos
        def __set__(self, v):
            self.reference_pos = v

    property tid:
        def __get__(self):
            return self.reference_id
        def __set__(self, v):
            self.reference_id = v

    property n:
        def __get__(self):
            return self.nsegments
        def __set__(self, v):
            self.nsegments = v


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
        None if is_del or is_refskip is set.

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
        """indel length for the position follwing the current pileup site.

        This quantity peeks ahead to the next cigar operation in this
        alignment. If the next operation is and insertion, indel will
        be positve. If the next operation is a deletion, it will be
        negation. 0 if the next operation is not an indel.

        """
        def __get__(self):
            return self._indel

    property level:
        """the level of the read in the "viewer" mode"""
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
        def __get__(self):
            return self._is_refskip

__all__ = [
    "AlignedSegment",
    "PileupColumn",
    "PileupRead"]
