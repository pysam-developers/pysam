# cython: embedsignature=True
# cython: profile=True
# adds doc-strings for sphinx
import tempfile
import os
import sys
import types
import itertools
import struct
import ctypes
import collections
import re
import platform
import warnings
import array

from cpython cimport PyErr_SetString, \
    PyBytes_Check, \
    PyUnicode_Check, \
    PyBytes_FromStringAndSize

from cpython cimport array

from cpython.version cimport PY_MAJOR_VERSION

cimport cython

########################################################################
########################################################################
########################################################################
## Python 3 compatibility functions
########################################################################
IS_PYTHON3 = PY_MAJOR_VERSION >= 3
cdef from_string_and_size(char* s, size_t length):
    if PY_MAJOR_VERSION < 3:
        return s[:length]
    else:
        return s[:length].decode("ascii")

# filename encoding (copied from lxml.etree.pyx)
cdef str _FILENAME_ENCODING
_FILENAME_ENCODING = sys.getfilesystemencoding()
if _FILENAME_ENCODING is None:
    _FILENAME_ENCODING = sys.getdefaultencoding()
if _FILENAME_ENCODING is None:
    _FILENAME_ENCODING = 'ascii'

#cdef char* _C_FILENAME_ENCODING
#_C_FILENAME_ENCODING = <char*>_FILENAME_ENCODING

cdef bytes _encodeFilename(object filename):
    """Make sure a filename is 8-bit encoded (or None)."""
    if filename is None:
        return None
    elif PyBytes_Check(filename):
        return filename
    elif PyUnicode_Check(filename):
        return filename.encode(_FILENAME_ENCODING)
    else:
        raise TypeError, u"Argument must be string or unicode."

cdef bytes _forceBytes(object s):
    u"""convert string or unicode object to bytes, assuming
    ascii encoding.
    """
    if PY_MAJOR_VERSION < 3:
        return s
    elif s is None:
        return None
    elif PyBytes_Check(s):
        return s
    elif PyUnicode_Check(s):
        return s.encode('ascii')
    else:
        raise TypeError, u"Argument must be string, bytes or unicode."

cdef inline bytes _forceCmdlineBytes(object s):
    return _forceBytes(s)

cdef _charptr_to_str(char* s):
    if PY_MAJOR_VERSION < 3:
        return s
    else:
        return s.decode("ascii")

cdef _forceStr(object s):
    """Return s converted to str type of current Python
    (bytes in Py2, unicode in Py3)"""
    if s is None:
        return None
    if PY_MAJOR_VERSION < 3:
        return s
    elif PyBytes_Check(s):
        return s.decode('ascii')
    else:
        # assume unicode
        return s

########################################################################
########################################################################
########################################################################
## Constants and global variables
########################################################################

# defines imported from samtools
DEF SEEK_SET = 0
DEF SEEK_CUR = 1
DEF SEEK_END = 2

cdef char* CODE2CIGAR= "MIDNSHP=X"
if IS_PYTHON3:
    CIGAR2CODE = dict( [y,x] for x,y in enumerate( CODE2CIGAR) )
else:
    CIGAR2CODE = dict( [ord(y),x] for x,y in enumerate( CODE2CIGAR) )
CIGAR_REGEX = re.compile( "(\d+)([MIDNSHP=X])" )

#####################################################################
# hard-coded constants
cdef int max_pos = 2 << 29

#####################################################################
#####################################################################
#####################################################################
## private factory methods
#####################################################################
cdef class AlignedSegment
cdef object makeAlignedSegment(bam1_t * src):
    '''enter src into AlignedSegment.'''
    # note that the following does not call __init__
    cdef AlignedSegment dest = AlignedSegment.__new__(AlignedSegment)
    dest._delegate = bam_dup1(src)
    return dest


cdef class PileupColumn
cdef makePileupColumn(bam_pileup1_t ** plp, int tid, int pos, int n_pu):
    # note that the following does not call __init__
     cdef PileupColumn dest = PileupColumn.__new__(PileupColumn)
     dest.plp = plp
     dest.tid = tid
     dest.pos = pos
     dest.n_pu = n_pu
     return dest

cdef class PileupRead
cdef makePileupRead(bam_pileup1_t * src):
    '''fill a  PileupRead object from a bam_pileup1_t * object.'''
    cdef PileupRead dest = PileupRead.__new__(PileupRead)
    dest._alignment = makeAlignedSegment(src.b)
    dest._qpos = src.qpos
    dest._indel = src.indel
    dest._level = src.level
    dest._is_del = src.is_del
    dest._is_head = src.is_head
    dest._is_tail = src.is_tail
    dest._is_refskip = src.is_refskip
    return dest

cdef convertBinaryTagToList( uint8_t * s ):
    """return bytesize, number of values list of values in s."""
    cdef char auxtype
    cdef uint8_t byte_size
    cdef int32_t nvalues

    # get byte size
    auxtype = s[0]
    byte_size = aux_type2size( auxtype )
    s += 1
    # get number of values in array
    nvalues = (<int32_t*>s)[0]
    s += 4
    # get values
    values = []
    if auxtype == 'c':
        for x from 0 <= x < nvalues:
            values.append((<int8_t*>s)[0])
            s += 1
    elif auxtype == 'C':
        for x from 0 <= x < nvalues:
            values.append((<uint8_t*>s)[0])
            s += 1
    elif auxtype == 's':
        for x from 0 <= x < nvalues:
            values.append((<int16_t*>s)[0])
            s += 2
    elif auxtype == 'S':
        for x from 0 <= x < nvalues:
            values.append((<uint16_t*>s)[0])
            s += 2
    elif auxtype == 'i':
        for x from 0 <= x < nvalues:
            values.append((<int32_t*>s)[0])
            s += 4
    elif auxtype == 'I':
        for x from 0 <= x < nvalues:
            values.append((<uint32_t*>s)[0])
            s += 4
    elif auxtype == 'f':
        for x from 0 <= x < nvalues:
            values.append((<float*>s)[0])
            s += 4

    return byte_size, nvalues, values


# valid types for SAM headers
VALID_HEADER_TYPES = {"HD" : dict,
                      "SQ" : list,
                      "RG" : list,
                      "PG" : list,
                      "CO" : list}

# order of records within SAM headers
VALID_HEADERS = ("HD", "SQ", "RG", "PG", "CO")

# default type conversions within SAM header records
VALID_HEADER_FIELDS = {"HD" : {"VN" : str, "SO" : str, "GO" : str},
                       "SQ" : {"SN" : str, "LN" : int, "AS" : str, 
                               "M5" : str, "SP" : str, "UR" : str,},
                       "RG" : {"ID" : str, "CN" : str, "DS" : str,
                               "DT" : str, "FO" : str, "KS" : str,
                               "LB" : str, "PG" : str, "PI" : str,
                               "PL" : str, "PU" : str, "SM" : str,},
                       "PG" : {"ID" : str, "PN" : str, "CL" : str, 
                               "PP" : str, "DS" : str, "VN" : str,},}

# output order of fields within records
VALID_HEADER_ORDER = {"HD" : ("VN", "SO", "GO"),
                      "SQ" : ("SN", "LN", "AS", "M5",
                               "UR", "SP"),
                      "RG" : ("ID", "SM", "LB", "DS", 
                              "PU", "PI", "CN", "DT",
                              "PL", "FO", "KS", "PG"),
                      "PG" : ("PN", "ID", "VN", "CL", 
                              "PP"),}


cdef class AlignmentFile:
    '''*(filename, mode=None, template = None,
         reference_names=None, reference_lengths = None,
         text=NULL, header=None,
         add_sq_text=False, check_header=True,
         check_sq=True)*

    A :term:`SAM`/:term:`BAM` formatted file. The file is
    automatically opened.

    *mode* should be ``r`` for reading or ``w`` for writing. The
    default is text mode (:term:`SAM`). For binary (:term:`BAM`) I/O
    you should append ``b`` for compressed or ``u`` for uncompressed
    :term:`BAM` output.  Use ``h`` to output header information in
    text (:term:`TAM`) mode.

    If ``b`` is present, it must immediately follow ``r`` or ``w``.
    Valid modes are ``r``, ``w``, ``wh``, ``rb``, ``wb``, ``wbu`` and
    ``wb0``. For instance, to open a :term:`BAM` formatted file for
    reading, type::

        f = pysam.AlignmentFile('ex1.bam','rb')

    If mode is not specified, we will try to auto-detect in the order
    'rb', 'r', thus both the following should work::

        f1 = pysam.AlignmentFile('ex1.bam')
        f2 = pysam.AlignmentFile('ex1.sam')

    If an index for a BAM file exists (.bai), it will be opened
    automatically. Without an index random access to reads via
    :meth:`fetch` and :meth:`pileup` is disabled.

    For writing, the header of a :term:`SAM` file/:term:`BAM` file can
    be constituted from several sources (see also the samtools format
    specification):

        1. If *template* is given, the header is copied from a another
           *AlignmentFile* (*template* must be of type *AlignmentFile*).

        2. If *header* is given, the header is built from a
           multi-level dictionary. The first level are the four types
           ('HD', 'SQ', ...). The second level are a list of lines,
           with each line being a list of tag-value pairs. The header
           is constructed first from all the defined fields, followed
           by user tags in alphabetical order.

        3. If *text* is given, new header text is copied from raw
           text.

        4. The names (*reference_names*) and lengths
           (*reference_lengths*) are supplied directly as lists.  By
           default, 'SQ' and 'LN' tags will be added to the header
           text. This option can be changed by unsetting the flag
           *add_sq_text*.

    For writing a CRAM file, the filename of the reference can be 
    added through a fasta formatted file (*reference_filename*)

    By default, if a file is opened in mode 'r', it is checked
    for a valid header (*check_header* = True) and a definition of
    chromosome names (*check_sq* = True).

    '''

    def __cinit__(self, *args, **kwargs ):
        self.htsfile = NULL
        self._filename = None
        self.is_bam = False
        self.is_stream = False
        self.is_cram = False
        self.is_remote = False
        
        self._open(*args, **kwargs)

        # allocate memory for iterator
        self.b = <bam1_t*>calloc(1, sizeof(bam1_t))

    def _isOpen(self):
        '''return true if htsfile has been opened.'''
        return self.htsfile != NULL

    def _hasIndex(self):
        '''return true if htsfile has an existing (and opened) index.'''
        return self.index != NULL

    def _open(self,
              filename,
              mode=None,
              AlignmentFile template=None,
              reference_names=None,
              reference_lengths=None,
              reference_filename=None,
              text=None,
              header=None,
              port=None,
              add_sq_text=True,
              check_header=True,
              check_sq=True,
              referencenames=None,
              referencelengths=None):
        '''open a sam, bam or cram formatted file.

        If _open is called on an existing file, the current file
        will be closed and a new file will be opened.
        '''
        # for backwards compatibility:
        if referencenames is not None:
            reference_names = referencenames
        if referencelengths is not None:
            reference_lengths = referencelengths

        # read mode autodetection
        if mode is None:
            try:
                self._open(filename, 'rb',
                           template=template,
                           reference_names=reference_names,
                           reference_lengths=reference_lengths,
                           reference_filename=reference_filename,
                           text=text,
                           header=header,
                           port=port,
                           check_header=check_header,
                           check_sq=check_sq)
                return
            except ValueError, msg:
                pass

            self._open(filename, 'r',
                       template=template,
                       reference_names=reference_names,
                       reference_lengths=reference_lengths,
                       reference_filename=reference_filename,
                       text=text,
                       header=header,
                       port=port,
                       check_header=check_header,
                       check_sq=check_sq)
            return

        assert mode in ("r","w","rb","wb", "wh",
                        "wbu", "rU", "wb0",
                        "rc", "wc"), \
            "invalid file opening mode `%s`" % mode

        # close a previously opened file
        if self.htsfile != NULL:
            self.close()

        # for htslib, wbu seems to not work
        if mode == "wbu":
            mode = "wb0"

        cdef bytes bmode = mode.encode('ascii')
        self._filename = filename = _encodeFilename(filename)

        # FIXME: Use htsFormat when it is available
        self.is_bam = len(mode) > 1 and mode[1] == 'b'
        self.is_cram = len(mode) > 1 and mode[1] == 'c'
        self.is_stream = filename == b"-"
        self.is_remote = filename.startswith(b"http:") or \
                         filename.startswith(b"ftp:")

        cdef char * ctext
        ctext = NULL

        if mode[0] == 'w':
            # open file for writing

            # header structure (used for writing)
            if template:
                self.header = bam_hdr_dup(template.header)
            elif header:
                self.header = self._buildHeader(header)
            else:
                # build header from a target names and lengths
                assert reference_names and reference_lengths, \
                    ("either supply options `template`, `header` "
                     "or  both `reference_names` and `reference_lengths` "
                     "for writing")
                assert len(reference_names) == len(reference_lengths), \
                    "unequal names and lengths of reference sequences"

                # allocate and fill header
                reference_names = [_forceBytes(ref) for ref in reference_names]
                self.header = bam_hdr_init()
                self.header.n_targets = len(reference_names)
                n = 0
                for x in reference_names:
                    n += len(x) + 1
                self.header.target_name = <char**>calloc(
                    n, sizeof(char*))
                self.header.target_len = <uint32_t*>calloc(
                    n, sizeof(uint32_t))
                for x from 0 <= x < self.header.n_targets:
                    self.header.target_len[x] = reference_lengths[x]
                    name = reference_names[x]
                    self.header.target_name[x] = <char*>calloc(
                        len(name) + 1, sizeof(char))
                    strncpy(self.header.target_name[x], name, len(name))

                # Optionally, if there is no text, add a SAM
                # compatible header to output file.
                if text is None and add_sq_text:
                    text = []
                    for x from 0 <= x < self.header.n_targets:
                        text.append("@SQ\tSN:%s\tLN:%s\n" % \
                                    (_forceStr(reference_names[x]), 
                                     reference_lengths[x]))
                    text = ''.join(text)

                if text is not None:
                    # copy without \0
                    text = _forceBytes(text)
                    ctext = text
                    self.header.l_text = strlen(ctext)
                    self.header.text = <char*>calloc(
                        strlen(ctext), sizeof(char))
                    memcpy(self.header.text, ctext, strlen(ctext))

            # open file (hts_open is synonym with sam_open)
            self.htsfile = hts_open(filename, bmode)

            # set filename with reference sequences. If no filename
            # is given, the CRAM reference arrays will be built from
            # the @SQ header in the header
            if self.is_cram and reference_filename:
                # note that fn_aux takes ownership, so create
                # a copy
                fn = _encodeFilename(reference_filename)
                self.htsfile.fn_aux = strdup(fn)

            # write header to htsfile
            if self.is_bam or self.is_cram or "h" in mode:
                sam_hdr_write(self.htsfile, self.header)

        elif mode[0] == "r":
            # open file for reading
            if (filename != b"-"
                and not self.is_remote
                and not os.path.exists(filename)):
                raise IOError("file `%s` not found" % filename)

            # open file (hts_open is synonym with sam_open)
            self.htsfile = hts_open(filename, bmode)
            if self.htsfile == NULL:
                raise ValueError(
                    "could not open file (mode='%s') - "
                    "is it SAM/BAM format?" % mode)

            # bam files require a valid header
            if self.is_bam or self.is_cram:
                self.header = sam_hdr_read(self.htsfile)
                if self.header == NULL:
                    raise ValueError(
                        "file does not have valid header (mode='%s') "
                        "- is it BAM format?" % mode )
            else:
                # in sam files it is optional (htsfile full of
                # unmapped reads)
                if check_header:
                    self.header = sam_hdr_read(self.htsfile)
                    if self.header == NULL:
                        raise ValueError(
                            "file does not have valid header (mode='%s') "
                            "- is it SAM format?" % mode )
                    # self.header.ignore_sam_err = True

            # disabled for autodetection to work needs to be disabled
            # so that reading from sam-files without headers works
            if check_sq and self.header.n_targets == 0:
                raise ValueError(
                    ("file header is empty (mode='%s') - "
                     "is it SAM/BAM format?") % mode)

        if self.htsfile == NULL:
            raise IOError("could not open file `%s`" % filename )

        # check for index and open if present
        cdef int format_index = -1
        if self.is_bam:
            format_index = HTS_FMT_BAI
        elif self.is_cram:
            format_index = HTS_FMT_CRAI

        if mode[0] == "r" and (self.is_bam or self.is_cram):

            # open index for remote files
            if self.is_remote:
                self.index = hts_idx_load(filename, format_index)
                if self.index == NULL:
                    warnings.warn(
                        "unable to open remote index for '%s'" % filename)
            else:
                if self.is_bam \
                   and not os.path.exists(filename + b".bai") \
                   and not os.path.exists(filename[:-4] + b".bai"):
                    self.index = NULL
                elif self.is_cram \
                     and not os.path.exists(filename + b".crai") \
                     and not os.path.exists(filename[:-4] + b".crai"):
                    self.index = NULL
                else:
                    # returns NULL if there is no index or index could
                    # not be opened
                    self.index = sam_index_load(self.htsfile,
                                                filename)
                    if self.index == NULL:
                        raise IOError(
                            "error while opening index for '%s'" %
                            filename)

            # save start of data section
            if not self.is_stream:
                self.start_offset = self.tell()

    def gettid(self, reference):
        '''
        convert :term:`reference` name into numerical :term:`tid`

        returns -1 if reference is not known.
        '''
        if not self._isOpen():
            raise ValueError("I/O operation on closed file")
        reference = _forceBytes(reference)
        return bam_name2id(self.header, reference)

    def getrname(self, tid):
        '''
        convert numerical :term:`tid` into :term:`reference` name.'''
        if not self._isOpen():
            raise ValueError("I/O operation on closed file")
        if not 0 <= tid < self.header.n_targets:
            raise ValueError("reference_id %i out of range 0<=tid<%i" % 
                             (tid, self.header.n_targets))
        return _charptr_to_str(self.header.target_name[tid])

    cdef char * _getrname(self, int tid):   # TODO unused
        '''
        convert numerical :term:`tid` into :term:`reference` name.'''
        if not self._isOpen():
            raise ValueError("I/O operation on closed file")

        if not 0 <= tid < self.header.n_targets:
            raise ValueError("tid %i out of range 0<=tid<%i" %
                             (tid, self.header.n_targets ))
        return self.header.target_name[tid]

    def _parseRegion(self,
                     reference=None,
                     start=None,
                     end=None,
                     region=None,
                     tid=None):
        '''parse region information.

        Raises ValueError for invalid regions.

        Returns a tuple of a flag, :term:`tid`, start and end. The
        flag indicates whether some coordinates were supplied.

        Note that region strings are 1-based, while *start* and *end* denote
        an interval in python coordinates.

        '''
        cdef int rtid
        cdef long long rstart
        cdef long long rend

        rtid = -1
        rstart = 0
        rend = max_pos
        if start != None:
            try:
                rstart = start
            except OverflowError:
                raise ValueError('start out of range (%i)' % start)

        if end != None:
            try:
                rend = end
            except OverflowError:
                raise ValueError('end out of range (%i)' % end)

        if region:
            region = _forceStr(region)
            parts = re.split("[:-]", region)
            reference = parts[0]
            if len(parts) >= 2:
                rstart = int(parts[1]) - 1
            if len(parts) >= 3:
                rend = int(parts[2])

        if not reference:
            return 0, 0, 0, 0

        if tid is not None:
            rtid = tid
        else:
            rtid = self.gettid(reference)

        if rtid < 0:
            raise ValueError(
                "invalid reference `%s`" % reference)
        if rstart > rend:
            raise ValueError(
                'invalid coordinates: start (%i) > end (%i)' % (rstart, rend))
        if not 0 <= rstart < max_pos:
            raise ValueError('start out of range (%i)' % rstart)
        if not 0 <= rend <= max_pos:
            raise ValueError('end out of range (%i)' % rend)

        return 1, rtid, rstart, rend

    def reset(self):
        '''reset file position to beginning of file just after
        the header.'''
        return self.seek(self.start_offset, 0)

    def seek(self, uint64_t offset, int where = 0):
        '''move file pointer to position *offset*, see
        :meth:`pysam.AlignmentFile.tell`.
        '''

        if not self._isOpen():
            raise ValueError("I/O operation on closed file")
        if not self.is_bam:
            raise NotImplementedError(
                "seek only available in bam files")
        if self.is_stream:
            raise OSError("seek no available in streams")

        return bgzf_seek(hts_get_bgzfp(self.htsfile), offset, where)

    def tell(self):
        '''
        return current file position.
        '''
        if not self._isOpen():
            raise ValueError("I/O operation on closed file")
        if not (self.is_bam or self.is_cram):
            raise NotImplementedError(
                "seek only available in bam files")

        return bgzf_tell(hts_get_bgzfp(self.htsfile))

    def fetch(self,
              reference=None,
              start=None,
              end=None,
              region=None,
              tid=None,
              callback=None,
              until_eof=False,
              multiple_iterators=False):
        '''fetch aligned, i.e. mapped, reads in a :term:`region`
        using 0-based
        indexing. The region is specified by :term:`reference`,
        *start* and *end*. Alternatively, a samtools :term:`region`
        string can be supplied.

        Without *reference* or *region* all mapped reads will be
        fetched. The reads will be returned ordered by reference
        sequence, which will not necessarily be the order within the
        file. 

        If *until_eof* is given, all reads from the current file
        position will be returned in order as they are within the
        file. Using this option will also fetch unmapped reads.

        Set *multiple_iterators* to true if you will be using multiple
        iterators on the same file at the same time. The iterator
        returned will receive its own copy of a filehandle to the file
        effectively re-opening the file. Re-opening a file creates
        some overhead, so beware.

        If only *reference* is set, all reads aligned to *reference*
        will be fetched.

        Note that a :term:`SAM` file does not allow random access. If
        *region* or *reference* are given, an exception is raised.

        '''
        cdef int rtid, rstart, rend, has_coord

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        has_coord, rtid, rstart, rend = self._parseRegion(reference,
                                                          start,
                                                          end,
                                                          region,
                                                          tid)

        # Turn of re-opening if htsfile is a stream
        if self.is_stream:
            multiple_iterators = False

        if self.is_bam or self.is_cram:
            if not until_eof and not self.is_remote:
                if not self._hasIndex():
                    raise ValueError(
                        "fetch called on bamfile without index")

            if has_coord:
                return IteratorRowRegion(
                    self, rtid, rstart, rend, 
                    multiple_iterators=multiple_iterators)
            else:
                if until_eof:
                    return IteratorRowAll(
                        self,
                        multiple_iterators=multiple_iterators)
                else:
                    # AH: check - reason why no multiple_iterators for
                    # AllRefs?
                    return IteratorRowAllRefs(
                        self,
                        multiple_iterators=multiple_iterators)
        else:
            if has_coord:
                raise ValueError(
                    "fetching by region is not available for sam files")

            if callback:
                raise NotImplementedError(
                    "callback not implemented yet")

            if self.header == NULL:
                raise ValueError(
                    "fetch called for htsfile without header")

            # check if targets are defined
            # give warning, sam_read1 segfaults
            if self.header.n_targets == 0:
                warnings.warn("fetch called for htsfile without header")
                
            return IteratorRowAll(self,
                                  multiple_iterators=multiple_iterators)

    def head(self, n, multiple_iterators=True):
        '''return iterator over the first n alignments. 

        This is useful for inspecting the bam-file.

        *multiple_iterators* is set to True by default in order to
        avoid changing the current file position.
        '''
        return IteratorRowHead(self, n,
                               multiple_iterators=multiple_iterators)

    def mate(self,
             AlignedSegment read):
        '''return the mate of :class:`AlignedSegment` *read*.

        Throws a ValueError if read is unpaired or the mate
        is unmapped.

        .. note::

            Calling this method will change the file position.
            This might interfere with any iterators that have
            not re-opened the file.

        .. note::
  
           This method is too slow for high-throughput processing.
           If a read needs to be processed with its mate, work
           from a read name sorted file or, better, cache reads.

        '''
        cdef uint32_t flag = read._delegate.core.flag

        if flag & BAM_FPAIRED == 0:
            raise ValueError("read %s: is unpaired" %
                             (read.query_name))
        if flag & BAM_FMUNMAP != 0:
            raise ValueError("mate %s: is unmapped" %
                             (read.query_name))

        # xor flags to get the other mate
        cdef int x = BAM_FREAD1 + BAM_FREAD2
        flag = (flag ^ x) & x

        # Make sure to use a separate file to jump around
        # to mate as otherwise the original file position
        # will be lost
        # The following code is not using the C API and
        # could thus be made much quicker, for example
        # by using tell and seek.
        for mate in self.fetch(
                read._delegate.core.mpos,
                read._delegate.core.mpos + 1,
                tid=read._delegate.core.mtid,
                multiple_iterators=True):
            if mate.flag & flag != 0 and \
               mate.query_name == read.query_name:
                break
        else:
            raise ValueError("mate not found")
        
        return mate

    def count(self,
              reference=None,
              start=None,
              end=None,
              region=None,
              until_eof=False):
        '''*(reference = None, start = None, end = None,
        region = None, callback = None, until_eof = False)*

        count reads :term:`region` using 0-based indexing. The region
        is specified by :term:`reference`, *start* and
        *end*. Alternatively, a samtools :term:`region` string can be
        supplied.

        Note that a :term:`SAM` file does not allow random access. If
        *region* or *reference* are given, an exception is raised.
        '''
        cdef AlignedSegment read
        cdef long counter = 0

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )
            
        for read in self.fetch(reference=reference,
                               start=start,
                               end=end,
                               region=region,
                               until_eof=until_eof):
            counter += 1

        return counter

    def pileup( self,
                reference = None,
                start = None,
                end = None,
                region = None,
                **kwargs ):
        '''perform a :term:`pileup` within a :term:`region`. The region is
        specified by :term:`reference`, *start* and *end* (using
        0-based indexing).  Alternatively, a samtools *region* string
        can be supplied.

        Without *reference* or *region* all reads will be used for the
        pileup. The reads will be returned ordered by
        :term:`reference` sequence, which will not necessarily be the
        order within the file.

        The method returns an iterator of type
        :class:`pysam.IteratorColumn` unless a *callback is
        provided. If a *callback* is given, the callback will be
        executed for each column within the :term:`region`.

        Note that :term:`SAM` formatted files do not allow random
        access.  In these files, if a *region* or *reference* are
        given an exception is raised.

        Optional *kwargs* to the iterator:

        stepper
           The stepper controlls how the iterator advances.
           Possible options for the stepper are

           ``all``
              skip reads in which any of the following flags are set:
              BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP

           ``nofilter``
              uses every single read
              

           ``samtools``
              same filter and read processing as in :term:`csamtools`
              pileup. This requires a *fastafile* to be given.


        fastafile
           A :class:`~pysam.FastaFile` object. This is required for
           some of the steppers.

        mask
           Skip all reads with bits set in mask if mask=True.

        max_depth
           Maximum read depth permitted. The default limit is *8000*.

        truncate

           By default, the samtools pileup engine outputs all reads
           overlapping a region (see note below).  If truncate is True
           and a region is given, only output columns in the exact
           region specificied.

        .. note::

            *all* reads which overlap the region are returned. The
            first base returned will be the first base of the first
            read *not* necessarily the first base of the region used
            in the query.

        '''
        cdef int rtid, rstart, rend, has_coord

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        has_coord, rtid, rstart, rend = self._parseRegion(
            reference, start, end, region )

        if self.is_bam or self.is_cram:
            if not self._hasIndex():
                raise ValueError("no index available for pileup")

            if has_coord:
                return IteratorColumnRegion(self,
                                            tid = rtid,
                                            start = rstart,
                                            end = rend,
                                            **kwargs )
            else:
                return IteratorColumnAllRefs(self, **kwargs )

        else:
            raise NotImplementedError( "pileup of samfiles not implemented yet" )

    @cython.boundscheck(False)  # we do manual bounds checking
    def count_coverage(self, chr, start, stop, quality_threshold = 15,
                       read_callback = 'all'):
        """Count ACGT in a part of a AlignmentFile. 
        Return 4 array.arrays of length = stop - start,
        in order A C G T.
        
        @quality_threshold is the minimum quality score (in phred) a
        base has to reach to be counted.  Possible @read_callback
        values are

        ``all``
`            skip reads in which any of the following
             flags are set: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL,
             BAM_FDUP

         ``nofilter``
             uses every single read

        Alternatively, @read_callback can be a function ```check_read(read)``1
        that should return True only for those reads that shall be included in
        the counting.

        """
        
        cdef int _start = start
        cdef int _stop = stop
        cdef int length = _stop - _start
        cdef array.array int_array_template = array.array('L', [])
        cdef array.array count_a
        cdef array.array count_c
        cdef array.array count_g
        cdef array.array count_t
        count_a = array.clone(int_array_template, length, zero=True)
        count_c = array.clone(int_array_template, length, zero=True)
        count_g = array.clone(int_array_template, length, zero=True)
        count_t = array.clone(int_array_template, length, zero=True)

        cdef char * seq
        cdef array.array quality
        cdef int qpos
        cdef int refpos
        cdef int c = 0
        cdef int _threshold = quality_threshold
        for read in self.fetch(chr, start, stop):
            if read_callback == 'all':
                if (read.flag & (0x4 | 0x100 | 0x200 | 0x400)):
                    continue
            elif read_callback == 'nofilter':
                pass
            else:
                if not read_callback(read):
                    continue
            seq = read.seq
            quality = read.query_qualities
            for qpos, refpos in read.get_aligned_pairs(True):
                if qpos is not None and refpos is not None and _start <= refpos < _stop:
                    if quality[qpos] > quality_threshold:
                        if seq[qpos] == 'A':
                            count_a.data.as_ulongs[refpos - _start] += 1
                        if seq[qpos] == 'C':
                            count_c.data.as_ulongs[refpos - _start] += 1
                        if seq[qpos] == 'G':
                            count_g.data.as_ulongs[refpos - _start] += 1
                        if seq[qpos] == 'T':
                            count_t.data.as_ulongs[refpos - _start] += 1
        return count_a, count_c, count_g, count_t

    def close(self):
        '''
        closes the :class:`pysam.AlignmentFile`.'''
        if self.htsfile != NULL:
            hts_close(self.htsfile)
            hts_idx_destroy(self.index);
            self.htsfile = NULL

    def __dealloc__(self):
        # remember: dealloc cannot call other methods
        # note: no doc string
        # note: __del__ is not called.

        # FIXME[kbj]: isn't self.close a method?  I've been duplicating
        # close within __dealloc__ (see BCFFile.__dealloc__).  Not a pretty
        # solution and perhaps unnecessary given that calling self.close has
        # been working for years.

        self.close()
        bam_destroy1(self.b)
        if self.header != NULL:
            bam_hdr_destroy(self.header)
            
    cpdef int write(self, AlignedSegment read) except -1:
        '''
        write a single :class:`pysam.AlignedSegment` to disk.

        returns the number of bytes written.
        '''
        if not self._isOpen():
            return 0

        cdef int ret = sam_write1(self.htsfile,
                                  self.header,
                                  read._delegate)

        # kbj: Still need to raise an exception with except -1. Otherwise
        #      when ret == -1 we get a "SystemError: error return without
        #      exception set".
        if ret < 0:
            raise ValueError('sam write failed')

        return ret

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    ###############################################################
    ###############################################################
    ###############################################################
    ## properties
    ###############################################################
    property filename:
        '''filename associated with this object.'''
        def __get__(self):
            return self._filename

    property nreferences:
        '''number of :term:`reference` sequences in the file.'''
        def __get__(self):
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
            return self.header.n_targets

    property references:
        """tuple with the names of :term:`reference` sequences."""
        def __get__(self):
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
            t = []
            for x from 0 <= x < self.header.n_targets:
                t.append(_charptr_to_str(self.header.target_name[x]))
            return tuple(t)

    property lengths:
        """tuple of the lengths of the :term:`reference` sequences. The
        lengths are in the same order as
        :attr:`pysam.AlignmentFile.references`

        """
        def __get__(self):
            if not self._isOpen():
                raise ValueError("I/O operation on closed file")
            t = []
            for x from 0 <= x < self.header.n_targets:
                t.append(self.header.target_len[x])
            return tuple(t)

    property mapped:
        """total number of mapped alignments according
        to the statistics recorded in the index.
        """
        def __get__(self):
            self._checkIndex()
            cdef int tid
            cdef uint64_t total = 0
            cdef uint64_t mapped, unmapped
            for tid from 0 <= tid < self.header.n_targets:
                hts_idx_get_stat(self.index, tid, &mapped, &unmapped)
                total += mapped
            return total

    def _checkIndex(self):
        '''check if index is present. Otherwise raise
        an error.'''
        if not self._isOpen():
            raise ValueError("I/O operation on closed file")
        if not self.is_bam and not self.is_cram:
            raise AttributeError(
                "AlignmentFile.mapped only available in bam files")
        if self.index == NULL:
            raise ValueError(
                "mapping information not recorded in index "
                "or index not available")


    property unmapped:
        """total number of unmapped reads according
        to the statistics recorded in the index.
        """
        def __get__(self):
            self._checkIndex()
            cdef int tid
            cdef uint64_t total = 0
            cdef uint64_t mapped, unmapped
            for tid from 0 <= tid < self.header.n_targets:
                hts_idx_get_stat(self.index, tid, &mapped, &unmapped)
                total += unmapped
            return total

    property nocoordinate:
        """total number of reads without coordinates according
        to the statistics recorded in the index.
        """
        def __get__(self):
            self._checkIndex()
            return hts_idx_get_n_no_coor(self.index)

    property text:
        '''full contents of the :term:`sam file` header as a string
    
        See :attr:`pysam.AlignmentFile.header` to get a parsed
        representation of the header.

        '''
        def __get__(self):
            if not self._isOpen():
                raise ValueError( "I/O operation on closed file" )
            return from_string_and_size(self.header.text, self.header.l_text)

    property header:
        '''header information within the :term:`sam file`. The records and
        fields are returned as a two-level dictionary. 

        The first level contains the record (``HD``, ``SQ``, etc) and
        the second level contains the fields (``VN``, ``LN``, etc).
        
        The parser is validating and will raise an AssertionError if
        if encounters any record or field tags that are not part of
        the SAM specification. Use the
        :attr:`pysam.AlignmentFile.text` attribute to get the unparsed
        header.

        The parsing follows the SAM format specification with the
        exception of the ``CL`` field. This option will consume the
        rest of a header line irrespective of any additional fields.
        This behaviour has been added to accommodate command line
        options that contain characters that are not valid field
        separators.

        '''
        def __get__(self):
            if not self._isOpen():
                raise ValueError( "I/O operation on closed file" )

            result = {}
            
            if self.header.text != NULL:
                # convert to python string (note: call self.text to
                # create 0-terminated string)
                t = self.text
                for line in t.split("\n"):
                    if not line.strip(): continue
                    assert line.startswith("@"), \
                        "header line without '@': '%s'" % line
                    fields = line[1:].split("\t")
                    record = fields[0]
                    assert record in VALID_HEADER_TYPES, \
                        "header line with invalid type '%s': '%s'" % (record, line)

                    # treat comments
                    if record == "CO":
                        if record not in result:
                            result[record] = []
                        result[record].append("\t".join( fields[1:]))
                        continue
                    # the following is clumsy as generators do not work?
                    x = {}

                    for idx, field in enumerate(fields[1:]):
                        if ":" not in field: 
                            raise ValueError("malformatted header: no ':' in field" )
                        key, value = field.split(":", 1)
                        if key in ("CL",):
                            # special treatment for command line
                            # statements (CL). These might contain
                            # characters that are non-conformant with
                            # the valid field separators in the SAM
                            # header. Thus, in contravention to the
                            # SAM API, consume the rest of the line.
                            key, value = "\t".join(fields[idx+1:]).split(":", 1)
                            x[key] = VALID_HEADER_FIELDS[record][key](value)
                            break

                        # uppercase keys must be valid
                        if key in VALID_HEADER_FIELDS[record]:
                            x[key] = VALID_HEADER_FIELDS[record][key](value)
                        # lowercase are permitted for user fields
                        elif not key.isupper():
                            x[key] = value
                        else:
                            raise ValueError(
                                "unknown field code '%s' in record '%s'" %
                                (key, record))

                    if VALID_HEADER_TYPES[record] == dict:
                        if record in result:
                            raise ValueError(
                                "multiple '%s' lines are not permitted" % record)

                        result[record] = x
                    elif VALID_HEADER_TYPES[record] == list:
                        if record not in result: result[record] = []
                        result[record].append(x)

                # if there are no SQ lines in the header, add the reference names
                # from the information in the bam file.
                #
                # Background: c-samtools keeps the textual part of the
                # header separate from the list of reference names and
                # lengths. Thus, if a header contains only SQ lines,
                # the SQ information is not part of the textual header
                # and thus are missing from the output. See issue 84.
                if "SQ" not in result:
                    sq = []
                    for ref, length in zip(self.references, self.lengths):
                        sq.append({'LN': length, 'SN': ref })
                    result["SQ"] = sq

            return result

    def _buildLine(self, fields, record):
        '''build a header line from *fields* dictionary for *record*'''

        # TODO: add checking for field and sort order
        line = ["@%s" % record]
        # comment
        if record == "CO":
            line.append(fields)
        # user tags
        elif record.islower():
            for key in sorted(fields):
                line.append("%s:%s" % (key, str(fields[key])))
        # defined tags
        else:
            # write fields of the specification
            for key in VALID_HEADER_ORDER[record]:
                if key in fields:
                    line.append("%s:%s" % (key, str(fields[key])))
            # write user fields
            for key in fields:
                if not key.isupper():
                    line.append("%s:%s" % (key, str(fields[key])))

        return "\t".join(line)

    cdef bam_hdr_t * _buildHeader(self, new_header):
        '''return a new header built from a dictionary in *new_header*.

        This method inserts the text field, target_name and target_len.
        '''

        lines = []

        # check if hash exists

        # create new header and copy old data
        cdef bam_hdr_t * dest

        dest = bam_hdr_init()

        # first: defined tags
        for record in VALID_HEADERS:
            if record in new_header:
                ttype = VALID_HEADER_TYPES[record]
                data = new_header[record]
                if type(data) != type(ttype()):
                    raise ValueError(
                        "invalid type for record %s: %s, expected %s" %
                        (record, type(data), type(ttype())))
                if type(data) is dict:
                    lines.append(self._buildLine(data, record))
                else:
                    for fields in new_header[record]:
                        lines.append(self._buildLine(fields, record))

        # then: user tags (lower case), sorted alphabetically
        for record, data in sorted(new_header.items()):
            if record in VALID_HEADERS: continue
            if type( data ) is dict:
                lines.append( self._buildLine( data, record ) )
            else:
                for fields in new_header[record]:
                    lines.append( self._buildLine( fields, record ) )

        text = "\n".join(lines) + "\n"
        if dest.text != NULL: free( dest.text )
        dest.text = <char*>calloc( len(text), sizeof(char))
        dest.l_text = len(text)
        cdef bytes btext = text.encode('ascii')
        strncpy( dest.text, btext, dest.l_text )

        cdef bytes bseqname
        # collect targets
        if "SQ" in new_header:
            seqs = []
            for fields in new_header["SQ"]:
                try:
                    seqs.append( (fields["SN"], fields["LN"] ) )
                except KeyError:
                    raise KeyError( "incomplete sequence information in '%s'" % str(fields))

            dest.n_targets = len(seqs)
            dest.target_name = <char**>calloc(dest.n_targets, sizeof(char*))
            dest.target_len = <uint32_t*>calloc(dest.n_targets, sizeof(uint32_t))

            for x from 0 <= x < dest.n_targets:
                seqname, seqlen = seqs[x]
                dest.target_name[x] = <char*>calloc(
                    len(seqname) + 1, sizeof(char))
                bseqname = seqname.encode('ascii')
                strncpy(dest.target_name[x], bseqname,
                        len(seqname) + 1)
                dest.target_len[x] = seqlen

        return dest

    ###############################################################
    ###############################################################
    ###############################################################
    ## file-object like iterator access
    ## note: concurrent access will cause errors (see IteratorRow
    ## and multiple_iterators)
    ## Possible solutions: deprecate or open new file handle
    ###############################################################
    def __iter__(self):
        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        if not self.is_bam and self.header.n_targets == 0:
            raise NotImplementedError(
                "can not iterate over samfile without header")
        return self

    cdef bam1_t * getCurrent( self ):
        return self.b

    cdef int cnext(self):
        '''
        cversion of iterator. Used by :class:`pysam.AlignmentFile.IteratorColumn`.
        '''
        return sam_read1(self.htsfile,
                         self.header,
                         self.b)

    def __next__(self):
        """
        python version of next().
        """
        cdef int ret = self.cnext()
        if (ret >= 0):
            return makeAlignedSegment(self.b)
        elif ret == -2:
            raise IOError('truncated file')
        else:
            raise StopIteration


cdef class IteratorRow:
    '''abstract base class for iterators over mapped reads.

    Various iterators implement different behaviours for wrapping around
    contig boundaries. Examples include:

    :class:`pysam.IteratorRowRegion`
        iterate within a single contig and a defined region.

    :class:`pysam.IteratorRowAll`
        iterate until EOF. This iterator will also include unmapped reads.

    :class:`pysam.IteratorRowAllRefs`
        iterate over all reads in all reference sequences.

    The method :meth:`AlignmentFile.fetch` returns an IteratorRow.

    .. note::

        It is usually not necessary to create an object of this class
        explicitely. It is returned as a result of call to a
        :meth:`AlignmentFile.fetch`.

    '''

    def __init__(self, AlignmentFile samfile, int multiple_iterators=False):
        
        if not samfile._isOpen():
            raise ValueError("I/O operation on closed file")

        # makes sure that samfile stays alive as long as the
        # iterator is alive
        self.samfile = samfile

        # reopen the file - note that this makes the iterator
        # slow and causes pileup to slow down significantly.
        if multiple_iterators:
            self.htsfile = hts_open(samfile._filename, 'r')
            assert self.htsfile != NULL
            # read header - required for accurate positioning
            # could a tell/seek work?
            self.header = sam_hdr_read(self.htsfile)
            assert self.header != NULL
            self.owns_samfile = True
        else:
            self.htsfile = self.samfile.htsfile
            self.owns_samfile = False
            self.header = self.samfile.header

        self.retval = 0

        self.b = bam_init1()

    def __dealloc__(self):
        bam_destroy1(self.b)
        if self.owns_samfile:
            hts_close(self.htsfile)
            bam_hdr_destroy(self.header)


cdef class IteratorRowRegion(IteratorRow):
    """*(AlignmentFile samfile, int tid, int beg, int end,
    int multiple_iterators=False)*

    iterate over mapped reads in a region.

    .. note::

        It is usually not necessary to create an object of this class
        explicitely. It is returned as a result of call to a
        :meth:`AlignmentFile.fetch`.

    """

    def __init__(self, AlignmentFile samfile,
                 int tid, int beg, int end,
                 int multiple_iterators=False):

        IteratorRow.__init__(self, samfile,
                             multiple_iterators=multiple_iterators)

        if not samfile._hasIndex():
            raise ValueError("no index available for iteration")

        self.iter = sam_itr_queryi(
            self.samfile.index,
            tid,
            beg,
            end)
    
    def __iter__(self):
        return self

    cdef bam1_t * getCurrent(self):
        return self.b

    cdef int cnext(self):
        '''cversion of iterator. Used by IteratorColumn'''
        self.retval = hts_itr_next(hts_get_bgzfp(self.htsfile),
                                   self.iter,
                                   self.b,
                                   self.htsfile)

    def __next__(self):
        """python version of next().
        """
        self.cnext()
        if self.retval >= 0:
            return makeAlignedSegment(self.b)
        elif self.retval == -2:
            # Note: it is currently not the case that hts_iter_next
            # returns -2 for a truncated file.
            # See https://github.com/pysam-developers/pysam/pull/50#issuecomment-64928625
            raise IOError('truncated file')
        else:
            raise StopIteration

    def __dealloc__(self):
        hts_itr_destroy(self.iter)


cdef class IteratorRowHead(IteratorRow):
    """*(AlignmentFile samfile, n, int multiple_iterators=False)*

    iterate over first n reads in *samfile*

    .. note::
        It is usually not necessary to create an object of this class
        explicitly. It is returned as a result of call to a
        :meth:`AlignmentFile.head`.

    """

    def __init__(self, AlignmentFile samfile, int n,
                 int multiple_iterators=False):

        IteratorRow.__init__(self, samfile,
                             multiple_iterators=multiple_iterators)

        self.max_rows = n
        self.current_row = 0

    def __iter__(self):
        return self

    cdef bam1_t * getCurrent( self ):
        return self.b

    cdef int cnext(self):
        '''cversion of iterator. Used by IteratorColumn'''
        return sam_read1(self.htsfile,
                         self.samfile.header,
                         self.b)

    def __next__(self):
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        if self.current_row >= self.max_rows:
            raise StopIteration

        cdef int ret = self.cnext()
        if (ret >= 0):
            self.current_row += 1
            return makeAlignedSegment( self.b )
        elif (ret == -2):
            raise IOError('truncated file')
        else:
            raise StopIteration


cdef class IteratorRowAll(IteratorRow):
    """*(AlignmentFile samfile, int multiple_iterators=False)*

    iterate over all reads in *samfile*

    .. note::

        It is usually not necessary to create an object of this class
        explicitly. It is returned as a result of call to a
        :meth:`AlignmentFile.fetch`.

    """

    def __init__(self, AlignmentFile samfile,
                 int multiple_iterators=False):

        IteratorRow.__init__(self, samfile,
                             multiple_iterators=multiple_iterators)

    def __iter__(self):
        return self

    cdef bam1_t * getCurrent( self ):
        return self.b

    cdef int cnext(self):
        '''cversion of iterator. Used by IteratorColumn'''
        return sam_read1(self.htsfile,
                         self.samfile.header,
                         self.b)

    def __next__(self):
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        cdef int ret = self.cnext()
        if (ret >= 0):
            return makeAlignedSegment(self.b)
        elif (ret == -2):
            raise IOError('truncated file')
        else:
            raise StopIteration


cdef class IteratorRowAllRefs(IteratorRow):
    """iterates over all mapped reads by chaining iterators over each
    reference

    .. note::
        It is usually not necessary to create an object of this class
        explicitely. It is returned as a result of call to a
        :meth:`AlignmentFile.fetch`.

    """

    def __init__(self, AlignmentFile samfile,
                 multiple_iterators=False):

        IteratorRow.__init__(self, samfile,
                             multiple_iterators=multiple_iterators)

        if not samfile._hasIndex():
            raise ValueError("no index available for fetch")

        self.tid = -1

    def nextiter(self):
        # get a new iterator for a chromosome. The file
        # will not be re-opened.
        self.rowiter = IteratorRowRegion(self.samfile,
                                         self.tid,
                                         0,
                                         1<<29)
        # set htsfile and header of the rowiter
        # to the values in this iterator to reflect multiple_iterators
        self.rowiter.htsfile = self.htsfile
        self.rowiter.header = self.header

        # make sure the iterator understand that IteratorRowAllRefs
        # has ownership
        self.rowiter.owns_samfile = False

    def __iter__(self):
        return self

    def __next__(self):
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        # Create an initial iterator
        if self.tid == -1:
            if not self.samfile.nreferences:
                raise StopIteration
            self.tid = 0
            self.nextiter()

        while 1:
            self.rowiter.cnext()

            # If current iterator is not exhausted, return aligned read
            if self.rowiter.retval > 0:
                return makeAlignedSegment(self.rowiter.b)

            self.tid += 1

            # Otherwise, proceed to next reference or stop
            if self.tid < self.samfile.nreferences:
                self.nextiter()
            else:
                raise StopIteration


cdef class IteratorRowSelection(IteratorRow):
    """*(AlignmentFile samfile)*

    iterate over reads in *samfile* at a given list of file positions.

    .. note::
        It is usually not necessary to create an object of this class
        explicitely. It is returned as a result of call to a :meth:`AlignmentFile.fetch`.
    """

    def __init__(self, AlignmentFile samfile, positions, int multiple_iterators=True):

        IteratorRow.__init__(self, samfile, multiple_iterators=multiple_iterators)

        self.positions = positions
        self.current_pos = 0

    def __iter__(self):
        return self

    cdef bam1_t * getCurrent(self):
        return self.b

    cdef int cnext(self):
        '''cversion of iterator'''

        # end iteration if out of positions
        if self.current_pos >= len(self.positions): return -1

        bgzf_seek(hts_get_bgzfp(self.htsfile),
                  self.positions[self.current_pos],
                  0)
        self.current_pos += 1
        return sam_read1(self.htsfile,
                         self.samfile.header,
                         self.b)

    def __next__(self):
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        cdef int ret = self.cnext()
        if (ret >= 0):
            return makeAlignedSegment(self.b)
        elif (ret == -2):
            raise IOError('truncated file')
        else:
            raise StopIteration


cdef int __advance_nofilter(void *data, bam1_t *b):
    '''advance without any read filtering.
    '''
    cdef __iterdata * d
    d = <__iterdata*>data
    return sam_itr_next(d.htsfile, d.iter, b)


cdef int __advance_all(void *data, bam1_t *b):
    '''only use reads for pileup passing basic
    filters:

    BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
    '''

    cdef __iterdata * d
    cdef mask = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP
    d = <__iterdata*>data
    cdef int ret = sam_itr_next(d.htsfile, d.iter, b)
    while ret >= 0 and b.core.flag & mask:
        ret = sam_itr_next(d.htsfile, d.iter, b)

    return ret


cdef int __advance_snpcalls(void * data, bam1_t * b):
    '''advance using same filter and read processing as in
    the samtools pileup.
    '''

    # Note that this method requries acces to some 
    # functions in the samtools code base and is thus
    # not htslib only.
    # The functions accessed in samtools are:
    # 1. bam_prob_realn
    # 2. bam_cap_mapQ
    cdef __iterdata * d
    d = <__iterdata*>data

    cdef int ret = sam_itr_next(d.htsfile, d.iter, b)
    cdef int skip = 0
    cdef int q
    cdef int is_cns = 1
    cdef int is_nobaq = 0
    cdef int capQ_thres = 0

    # reload sequence
    if d.fastafile != NULL and b.core.tid != d.tid:
        if d.seq != NULL:
            free(d.seq)
        d.tid = b.core.tid
        d.seq = faidx_fetch_seq(
            d.fastafile,
            d.header.target_name[d.tid],
            0, max_pos,
            &d.seq_len)

        if d.seq == NULL:
            raise ValueError(
                "reference sequence for '%s' (tid=%i) not found" % \
                (d.header.target_name[d.tid],
                 d.tid))

    while ret >= 0:
        skip = 0

        # realign read - changes base qualities
        if d.seq != NULL and is_cns and not is_nobaq: 
            bam_prob_realn(b, d.seq)

        if d.seq != NULL and capQ_thres > 10:
            q = bam_cap_mapQ(b, d.seq, capQ_thres)
            if q < 0:
                skip = 1
            elif b.core.qual > q:
                b.core.qual = q
        if b.core.flag & BAM_FUNMAP:
            skip = 1
        elif b.core.flag & 1 and not b.core.flag & 2:
            skip = 1

        if not skip:
            break
        # additional filters

        ret = sam_itr_next(d.htsfile, d.iter, b)

    return ret

cdef class IteratorColumn:
    '''abstract base class for iterators over columns.

    IteratorColumn objects wrap the pileup functionality of samtools.

    For reasons of efficiency, the iterator points to the current
    pileup buffer. The pileup buffer is updated at every iteration.
    This might cause some unexpected behavious. For example,
    consider the conversion to a list::

       f = AlignmentFile("file.bam", "rb")
       result = list( f.pileup() )

    Here, ``result`` will contain ``n`` objects of type
    :class:`PileupColumn` for ``n`` columns, but each object in
    ``result`` will contain the same information.

    The desired behaviour can be achieved by list comprehension::

       result = [ x.pileups() for x in f.pileup() ]

    ``result`` will be a list of ``n`` lists of objects of type
    :class:`PileupRead`.

    If the iterator is associated with a :class:`Fastafile` using the
    :meth:`addReference` method, then the iterator will export the
    current sequence via the methods :meth:`getSequence` and
    :meth:`seq_len`.

    Optional kwargs to the iterator:

    stepper
       The stepper controls how the iterator advances.

       Valid values are None, "all" (default), "nofilter" or "samtools".

       See AlignmentFile.pileup for description.
    
    fastafile
       A :class:`FastaFile` object

    max_depth
       maximum read depth. The default is 8000.

    '''

    def __cinit__( self, AlignmentFile samfile, **kwargs ):
        self.samfile = samfile
        # TODO
        # self.mask = kwargs.get("mask", BAM_DEF_MASK )
        self.fastafile = kwargs.get("fastafile", None)
        self.stepper = kwargs.get("stepper", None)
        self.max_depth = kwargs.get("max_depth", 8000)
        self.iterdata.seq = NULL
        self.tid = 0
        self.pos = 0
        self.n_plp = 0
        self.plp = NULL
        self.pileup_iter = <bam_plp_t>NULL

    def __iter__(self):
        return self

    cdef int cnext(self):
        '''perform next iteration.
        '''
        self.plp = bam_plp_auto( self.pileup_iter,
                                 &self.tid,
                                 &self.pos,
                                 &self.n_plp )

    cdef char * getSequence( self ):
        '''return current reference sequence underlying the iterator.
        '''
        return self.iterdata.seq

    property seq_len:
        '''current sequence length.'''
        def __get__(self): return self.iterdata.seq_len

    def addReference(self, Fastafile fastafile):
       '''
       add reference sequences in *fastafile* to iterator.'''
       self.fastafile = fastafile
       if self.iterdata.seq != NULL: free(self.iterdata.seq)
       self.iterdata.tid = -1
       self.iterdata.fastafile = self.fastafile.fastafile

    def hasReference(self):
        '''
        return true if iterator is associated with a reference'''
        return self.fastafile

    cdef setMask(self, mask):
        '''set masking flag in iterator.

        reads with bits set in *mask* will be skipped.
        '''
        raise NotImplementedError()
        # self.mask = mask
        # bam_plp_set_mask( self.pileup_iter, self.mask )

    cdef setupIteratorData( self,
                            int tid,
                            int start,
                            int end,
                            int multiple_iterators = 0 ):
        '''setup the iterator structure'''

        self.iter = IteratorRowRegion(self.samfile, tid, start, end, multiple_iterators)
        self.iterdata.htsfile = self.samfile.htsfile
        self.iterdata.iter = self.iter.iter
        self.iterdata.seq = NULL
        self.iterdata.tid = -1
        self.iterdata.header = self.samfile.header

        if self.fastafile is not None:
            self.iterdata.fastafile = self.fastafile.fastafile
        else:
            self.iterdata.fastafile = NULL

        # Free any previously allocated memory before reassigning
        # pileup_iter
        self._free_pileup_iter()

        if self.stepper is None or self.stepper == "all":
            self.pileup_iter = bam_plp_init(
                <bam_plp_auto_f>&__advance_all,
                &self.iterdata)
        elif self.stepper == "nofilter":
            self.pileup_iter = bam_plp_init(
                <bam_plp_auto_f>&__advance_nofilter,
                &self.iterdata)
        elif self.stepper == "samtools":
            self.pileup_iter = bam_plp_init(
                <bam_plp_auto_f>&__advance_snpcalls,
                &self.iterdata)
        else:
            raise ValueError(
                "unknown stepper option `%s` in IteratorColumn" % self.stepper)

        if self.max_depth:
            bam_plp_set_maxcnt(self.pileup_iter, self.max_depth)

        # bam_plp_set_mask( self.pileup_iter, self.mask )

    cdef reset( self, tid, start, end ):
        '''reset iterator position.

        This permits using the iterator multiple times without
        having to incur the full set-up costs.
        '''
        self.iter = IteratorRowRegion( self.samfile, tid, start, end, multiple_iterators = 0 )
        self.iterdata.iter = self.iter.iter

        # invalidate sequence if different tid
        if self.tid != tid:
            if self.iterdata.seq != NULL: free( self.iterdata.seq )
            self.iterdata.seq = NULL
            self.iterdata.tid = -1

        # self.pileup_iter = bam_plp_init( &__advancepileup, &self.iterdata )
        bam_plp_reset(self.pileup_iter)

    cdef _free_pileup_iter(self):
        '''free the memory alloc'd by bam_plp_init.

        This is needed before setupIteratorData allocates
        another pileup_iter, or else memory will be lost.
        '''
        if self.pileup_iter != <bam_plp_t>NULL:
            bam_plp_reset(self.pileup_iter)
            bam_plp_destroy(self.pileup_iter)
            self.pileup_iter = <bam_plp_t>NULL

    def __dealloc__(self):
        # reset in order to avoid memory leak messages for iterators
        # that have not been fully consumed
        self._free_pileup_iter()
        self.plp = <bam_pileup1_t*>NULL

        if self.iterdata.seq != NULL:
            free(self.iterdata.seq)
            self.iterdata.seq = NULL


cdef class IteratorColumnRegion(IteratorColumn):
    '''iterates over a region only.
    '''
    def __cinit__(self, AlignmentFile samfile,
                  int tid = 0,
                  int start = 0,
                  int end = max_pos,
                  int truncate = False,
                  **kwargs ):

        # initialize iterator
        self.setupIteratorData( tid, start, end, 1 )
        self.start = start
        self.end = end
        self.truncate = truncate

    def __next__(self):
        """python version of next().
        """

        while 1:
            self.cnext()
            if self.n_plp < 0:
                raise ValueError("error during iteration" )

            if self.plp == NULL:
                raise StopIteration
            
            if self.truncate:
                if self.start > self.pos: continue
                if self.pos >= self.end: raise StopIteration

            return makePileupColumn(&self.plp,
                                   self.tid,
                                   self.pos,
                                   self.n_plp)


cdef class IteratorColumnAllRefs(IteratorColumn):
    """iterates over all columns by chaining iterators over each reference
    """

    def __cinit__(self,
                  AlignmentFile samfile,
                  **kwargs):

        # no iteration over empty files
        if not samfile.nreferences:
            raise StopIteration

        # initialize iterator
        self.setupIteratorData(self.tid, 0, max_pos, 1)

    def __next__(self):
        """python version of next().
        """

        while 1:
            self.cnext()

            if self.n_plp < 0:
                raise ValueError("error during iteration" )

            # return result, if within same reference
            if self.plp != NULL:
                return makePileupColumn(&self.plp,
                                       self.tid,
                                       self.pos,
                                       self.n_plp)
                
            # otherwise, proceed to next reference or stop
            self.tid += 1
            if self.tid < self.samfile.nreferences:
                self.setupIteratorData(self.tid, 0, max_pos, 0)
            else:
                raise StopIteration


cdef inline int32_t _getQueryStart(bam1_t *src) except -1:
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


cdef inline int32_t _getQueryEnd(bam1_t *src) except -1:
    cdef uint32_t * cigar_p
    cdef uint32_t k, op
    cdef uint32_t end_offset = src.core.l_qseq

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

    if end_offset == 0:
        end_offset = src.core.l_qseq

    return end_offset


cdef inline object _getSequenceRange(bam1_t *src,
                                     uint32_t start, uint32_t end):
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

    return _charptr_to_str(seq)


cdef inline object _getQualitiesRange(bam1_t *src,
                                      uint32_t start, 
                                      uint32_t end):
    '''return an array of quality values.'''

    cdef uint8_t * p
    cdef uint32_t k

    p = pysam_bam_get_qual(src)
    if p[0] == 0xff:
        return None

    # 'B': unsigned char
    cdef array.array result = array.array('B', [0])
    array.resize(result, end - start)

    # copy data
    memcpy(result.data.as_voidptr, <void*>&p[start], end - start) 

    return result


def toQualityString(qualities):
    '''convert a list of quality score to the string
    representation used in the SAM format.'''
    if qualities is None:
        return None
    return "".join([chr(x+33) for x in qualities])
    

def fromQualityString(quality_string):
    '''return a list of quality scores from the
    stringn representation of quality scores used
    in the SAM format.'''
    if quality_string is None:
        return None
    return array.array('B', [ord(x)-33 for x in quality_string])


cdef inline uint8_t _get_value_code(value, value_type=None):
    '''guess type code for a *value*. If *value_type* is None,
    the type code will be inferred based on the Python type of
    *value*'''
    cdef uint8_t  type_code    
    cdef char * _char_type

    if value_type is None:
        if isinstance(value, int):
            type_code = 'i'
        elif isinstance(value, float):
            type_code = 'd'
        elif isinstance(value, str):
            type_code = 'Z'
        elif isinstance(value, bytes):
            type_code = 'Z'
        else:
            return 0
    else:
        if value_type not in 'Zidf':
            return 0
        value_type = _forceBytes(value_type)
        _char_type = value_type
        type_code = (<uint8_t*>_char_type)[0]

    return type_code


cdef inline _get_value_type(value, maximum_value=None):
    '''returns the value type of a value.

    If max is specified, the approprite type is
    returned for a range where value is the minimum.
    '''
    
    if maximum_value is None:
        maximum_value = value

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
            valuetype = b"A"
        else:
            valuetype = b'Z'

    return valuetype


cdef inline _pack_tags(tags):
    """pack a list of tags. Each tag is a tuple of (tag, tuple).
    
    Values are packed into the most space efficient data structure
    possible unless the tag contains a third field with the type code.

    Returns a fmt string and the associated list of arguments
    to used in a call to struct.pack_into.
    """
    fmts, args = ["<"], []

    for tag in tags:

        if len(tag) == 2:
            pytag, value = tag
            valuetype = None
        elif len(tag) == 3:
            pytag, value, valuetype = tag
        else:
            raise ValueError("malformatted tag: %s" % str(tag))

        if not type(pytag) is bytes:
            pytag = pytag.encode('ascii')

        datatype2format = {'c': 'b',
                           's': 'h',
                           'i': 'i',
                           'C': 'B',
                           'S': 'H',
                           'I': 'I',
                           'f': 'f',
                           'A': 'c',}

        t = type(value)
        if t is tuple or t is list:
            # binary tags are treated separately
            if valuetype is None:
                # automatically determine value type - first value
                # determines type. If there is a mix of types, the
                # result is undefined.
                valuetype = _get_value_type(min(value), max(value))
                            
            if valuetype not in datatype2format:
                raise ValueError("invalid value type '%s'" % valuetype)
            datafmt = "2sccI%i%s" % (len(value), datatype2format[valuetype])

            args.extend([pytag[:2], 
                         b"B",
                         valuetype,
                         len(value)] + list(value))
            fmts.append(datafmt)

        else:
            
            if valuetype is None:
                valuetype = _get_value_type(value)

            if valuetype == b"Z":
                fmt = "2sc%is" % (len(value)+1)
            else:
                fmt = "2sc%s" % datatype2format[valuetype]

            args.extend([pytag[:2],
                         valuetype,
                         value])

            fmts.append(fmt)

    return "".join(fmts), args
    

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

    def __dealloc__(self):
        bam_destroy1(self._delegate)

    def __str__(self):
        """return string representation of alignment.

        The representation is an approximate :term:`sam` format.

        An aligned read might not be associated with a :term:`AlignmentFile`.
        As a result :term:`tid` is shown instead of the reference name.

        Similarly, the tags field is returned in its parsed state.
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

    # Disabled so long as __cmp__ is a special method
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

    ########################################################
    ## Basic attributes in order of appearance in SAM format
    property query_name:
        """the query template name (None if not present)"""
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            if pysam_get_l_qname(src) == 0:
                return None
            return _charptr_to_str(<char *>pysam_bam_get_qname(src))

        def __set__(self, qname):
            if qname is None or len(qname) == 0:
                return
            qname = _forceBytes(qname)
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

           q = read.qual
           read.seq = read.seq[5:10]
           read.qual = q[5:10]

        The sequence is returned as it is stored in the BAM file. Some mappers
        might have stored a reverse complement of the original read 
        sequence.
        """
        def __get__(self):
            cdef bam1_t * src
            cdef char * s
            src = self._delegate

            if src.core.l_qseq == 0: return None

            return _getSequenceRange(src, 0, src.core.l_qseq)

        def __set__(self, seq):
            # samtools manages sequence and quality length memory together
            # if no quality information is present, the first byte says 0xff.
            cdef bam1_t * src
            cdef uint8_t * p
            cdef char * s
            cdef int l, k, nbytes_new, nbytes_old

            if seq == None:
                l = 0
            else:
                l = len(seq)                
                seq = _forceBytes(seq)

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

            cdef bam1_t * src
            cdef char * q

            src = self._delegate

            if src.core.l_qseq == 0:
                return None

            return _getQualitiesRange(src, 0, src.core.l_qseq)

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
            cdef array.array result = array.array('B', qual)

            # copy data
            memcpy(p, result.data.as_voidptr, l)
    

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
            cdef bam1_t * src
            cdef uint32_t start, end

            src = self._delegate

            if src.core.l_qseq == 0:
                return None

            start = _getQueryStart(src)
            end   = _getQueryEnd(src)

            return _getSequenceRange(src, start, end)

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
            cdef bam1_t * src
            cdef uint32_t start, end

            src = self._delegate

            if src.core.l_qseq == 0:
                return None

            start = _getQueryStart(src)
            end   = _getQueryEnd(src)

            return _getQualitiesRange(src, start, end)

    property query_alignment_start:
        """start index of the aligned query portion of the sequence (0-based,
        inclusive).

        This the index of the first base in :attr:`seq` that is not
        soft-clipped.

        """
        def __get__(self):
            return _getQueryStart(self._delegate)

    property query_alignment_end:
        """end index of the aligned query portion of the sequence (0-based,
        exclusive)"""
        def __get__(self):
            return _getQueryEnd(self._delegate)

    property query_alignment_length:
        """length of the aligned query sequence.

        This is equal to :attr:`qend` - :attr:`qstart`"""
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            return _getQueryEnd(src) - _getQueryStart(src)

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
        cdef uint32_t k, qpos
        cdef int op
        cdef uint32_t * cigar_p
        cdef bam1_t * src 

        src = self._delegate

        if not always and src.core.l_qseq:
            return src.core.l_qseq

        if pysam_get_n_cigar(src) == 0:
            return None

        qpos = 0
        cigar_p = pysam_bam_get_cigar(src)

        for k from 0 <= k < pysam_get_n_cigar(src):
            op = cigar_p[k] & BAM_CIGAR_MASK

            if op == BAM_CMATCH or op == BAM_CINS or \
               op == BAM_CSOFT_CLIP or \
               op == BAM_CEQUAL or op == BAM_CDIFF:
                qpos += cigar_p[k] >> BAM_CIGAR_SHIFT

        return qpos
            
    def get_aligned_pairs(self, matches_only = False):
        """a list of aligned read (query) and reference positions.
        For inserts, deletions, skipping either query or reference position may be None.

        If @matches_only is True, only matched bases are returned - no None on either side.

        Padding is currently not supported and leads to an exception
        
        """
        cdef uint32_t k, i, pos, qpos
        cdef int op
        cdef uint32_t * cigar_p
        cdef bam1_t * src 
        cdef int _matches_only

        _matches_only = bool(matches_only)

        src = self._delegate
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
                for i from pos <= i < pos + l:
                    result.append((qpos, i))
                    qpos += 1
                pos += l

            elif op == BAM_CINS or op == BAM_CSOFT_CLIP:
                if not _matches_only:
                    for i from pos <= i < pos + l:
                        result.append((qpos, None))
                        qpos += 1
                else:
                    qpos += l

            elif op == BAM_CDEL or op == BAM_CREF_SKIP:
                if not _matches_only:
                    for i from pos <= i < pos + l:
                        result.append((None, i))
                pos += l

            elif op == BAM_CHARD_CLIP:
                pass # advances neither

            elif op == BAM_CPAD:
                raise NotImplementedError("Padding (BAM_CPAD, 6) is currently not supported. Please implement. Sorry about that.")

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
        into the alignment record.. It can be set explicitely to one
        of the valid one-letter type codes. If unset, an appropriate
        type will be chosen automatically.

        An existing value of the same *tag* will be overwritten unless
        replace is set to False. This is usually not recommened as a
        tag may only appear once in the optional alignment section.

        If *value* is None, the tag will be deleted.
        """

        cdef int      value_size
        cdef uint8_t * value_ptr
        cdef uint8_t *existing_ptr
        cdef uint8_t  type_code
        cdef float    float_value
        cdef double   double_value
        cdef int32_t  int_value
        cdef bam1_t * src = self._delegate
        cdef char * _value_type
        
        if len(tag) != 2:
            raise ValueError('Invalid tag: %s' % tag)

        tag = _forceBytes(tag)
        if replace:
            existing_ptr = bam_aux_get(src, tag)
            if existing_ptr:
                bam_aux_del(src, existing_ptr)

        # setting value to None deletes a tag
        if value is None:
            return
        
        type_code = _get_value_code(value, value_type)
        if type_code == 0:
            raise ValueError("can't guess type or invalid type code specified")

        # Not Endian-safe, but then again neither is samtools!
        if type_code == 'Z':
            value = _forceBytes(value)
            value_ptr    = <uint8_t*><char*>value
            value_size   = len(value)+1
        elif type_code == 'i':
            int_value    = value
            value_ptr    = <uint8_t*>&int_value
            value_size   = sizeof(int32_t)
        elif type_code == 'd':
            double_value = value
            value_ptr    = <uint8_t*>&double_value
            value_size   = sizeof(double)
        elif type_code == 'f':
            float_value  = value
            value_ptr    = <uint8_t*>&float_value
            value_size   = sizeof(float)
        else:
            raise ValueError('Unsupported value_type in set_option')


        bam_aux_append(src,
                       tag,
                       type_code, 
                       value_size,
                       value_ptr)

    cpdef has_tag(self, tag):
        """returns true if the optional alignment section
        contains a given *tag*."""
        cdef uint8_t * v
        cdef int nvalues
        btag = _forceBytes(tag)
        v = bam_aux_get(self._delegate, btag)
        return v != NULL

    cpdef get_tag(self, tag):
        """retrieves data from the optional alignment section
        given a two-letter *tag* denoting the field.

        If *tag* is not present, a KeyError is raised.

        The returned value is cast into an appropriate python type.

        This method is the fastest way to access the optional
        alignment section if only few tags need to be retrieved.
        """
        cdef uint8_t * v
        cdef int nvalues
        btag = _forceBytes(tag)
        v = bam_aux_get(self._delegate, btag)
        if v == NULL:
            raise KeyError("tag '%s' not present" % tag)
        auxtype = chr(v[0])
        if auxtype == 'c' or auxtype == 'C' or auxtype == 's' or auxtype == 'S':
            return <int>bam_aux2i(v)
        elif auxtype == 'i' or auxtype == 'I':
            return <int32_t>bam_aux2i(v)
        elif auxtype == 'f' or auxtype == 'F':
            return <float>bam_aux2f(v)
        elif auxtype == 'd' or auxtype == 'D':
            return <double>bam_aux2f(v)
        elif auxtype == 'A':
            # there might a more efficient way
            # to convert a char into a string
            return '%c' % <char>bam_aux2A(v)
        elif auxtype == 'Z':
            return _charptr_to_str(<char*>bam_aux2Z(v))
        elif auxtype == 'B':
            bytesize, nvalues, values = convertBinaryTagToList(v + 1)
            return values
        else:
            raise ValueError("unknown auxilliary type '%s'" % auxtype)

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
                value = _charptr_to_str(<char*>bam_aux2Z(s))
                # +1 for NULL terminated string
                s += len(value) + 1
            elif auxtype == 'B':
                s += 1
                byte_size, nvalues, value = convertBinaryTagToList(s)
                # 5 for 1 char and 1 int
                s += 5 + (nvalues * byte_size) - 1
            else:
                raise KeyError("unknown type '%s'" % auxtype)

            s += 1

            result.append((_charptr_to_str(auxtag), value))

        return result

    def set_tags(self, tags):
        """sets the fields in the optional alignmest section with
        a list of (tag, value) tuples.

        The :term:`value type` of the values is determined from the
        python type. Optionally, a type may be given explicitely as
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
            fmt, args =_pack_tags(tags)
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
        def __get__(self): return self.query_name
        def __set__(self, v): self.query_name = v
    property tid:
        def __get__(self): return self.reference_id
        def __set__(self, v): self.reference_id = v
    property pos:
        def __get__(self): return self.reference_start
        def __set__(self, v): self.reference_start = v
    property mapq:
        def __get__(self): return self.mapping_quality
        def __set__(self, v): self.mapping_quality = v
    property rnext:
        def __get__(self): return self.next_reference_id
        def __set__(self, v): self.next_reference_id = v
    property pnext:
        def __get__(self):
            return self.next_reference_start
        def __set__(self, v):
            self.next_reference_start = v
    property cigar:
        def __get__(self):
            r = self.cigartuples
            if r is None:
                r = []
            return r
        def __set__(self, v): self.cigartuples = v
    property tlen:
        def __get__(self):
            return self.template_length
        def __set__(self, v):
            self.template_length = v
    property seq:
        def __get__(self): return self.query_sequence
        def __set__(self, v): self.query_sequence = v
    property qual:
        def __get__(self):
            return toQualityString(self.query_qualities)
        def __set__(self, v):
            self.query_qualities = fromQualityString(v)
    property alen:
        def __get__(self):
            return self.reference_length
        def __set__(self, v):
            self.reference_length = v
    property aend:
        def __get__(self):
            return self.reference_end
        def __set__(self, v):
            self.reference_end = v
    property rlen:
        def __get__(self):
            return self.query_length
        def __set__(self, v):
            self.query_length = v
    property query:
        def __get__(self):
            return self.query_alignment_sequence
        def __set__(self, v):
            self.query_alignment_sequence = v
    property qqual:
        def __get__(self):
            return toQualityString(self.query_alignment_qualities)
        def __set__(self, v):
            self.query_alignment_qualities = fromQualityString(v)
    property qstart:
        def __get__(self):
            return self.query_alignment_start
        def __set__(self, v):
            self.query_alignment_start = v
    property qend:
        def __get__(self):
            return self.query_alignment_end
        def __set__(self, v):
            self.query_alignment_end = v
    property qlen:
        def __get__(self):
            return self.query_alignment_length
        def __set__(self, v):
            self.query_alignment_length = v
    property mrnm:
        def __get__(self):
            return self.next_reference_id
        def __set__(self, v):
            self.next_reference_id = v
    property mpos:
        def __get__(self):
            return self.next_reference_start
        def __set__(self, v):
            self.next_reference_start = v
    property rname:
        def __get__(self):
            return self.reference_id
        def __set__(self, v):
            self.reference_id = v
    property isize:
        def __get__(self):
            return self.template_length
        def __set__(self, v):
            self.template_length = v
    property blocks:
        def __get__(self):
            return self.get_blocks()
    property aligned_pairs:
        def __get__(self):
            return self.get_aligned_pairs()
    property inferred_length:
        def __get__(self):
            return self.infer_query_length()
    property positions:
        def __get__(self):
            return self.get_reference_positions()
    property tags:
        def __get__(self):
            return self.get_tags()
        def __set__(self, tags):
            self.set_tags(tags)
    def overlap(self):
        return self.get_overlap()
    def opt(self, tag):
        return self.get_tag(tag)
    def setTag(self, tag, value, value_type=None, replace=True):
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
                pileups.append(makePileupRead(&(self.plp[0][x])))
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

    property indel:
        """indel length; 0 for no indel, positive for ins and negative            for del"""
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
        def __get__(self):
            return self._is_head

    property is_tail:
        def __get__(self):
            return self._is_tail

    property is_refskip:
        def __get__(self):
            return self._is_refskip


cdef class SNPCall:
    '''the results of a SNP call.'''
    cdef int _tid
    cdef int _pos
    cdef char _reference_base
    cdef char _genotype
    cdef int _consensus_quality
    cdef int _snp_quality
    cdef int _rms_mapping_quality
    cdef int _coverage

    property tid:
        '''the chromosome ID as is defined in the header'''
        def __get__(self):
            return self._tid

    property pos:
       '''nucleotide position of SNP.'''
       def __get__(self): return self._pos

    property reference_base:
       '''reference base at pos. ``N`` if no reference sequence supplied.'''
       def __get__(self): return from_string_and_size( &self._reference_base, 1 )

    property genotype:
       '''the genotype called.'''
       def __get__(self): return from_string_and_size( &self._genotype, 1 )

    property consensus_quality:
       '''the genotype quality (Phred-scaled).'''
       def __get__(self): return self._consensus_quality

    property snp_quality:
       '''the snp quality (Phred scaled) - probability of consensus being
       identical to reference sequence.'''
       def __get__(self): return self._snp_quality

    property mapping_quality:
       '''the root mean square (rms) of the mapping quality of all reads
       involved in the call.'''
       def __get__(self): return self._rms_mapping_quality

    property coverage:
       '''coverage or read depth - the number of reads involved in the call.'''
       def __get__(self): return self._coverage

    def __str__(self):

        return "\t".join( map(str, (
                    self.tid,
                    self.pos,
                    self.reference_base,
                    self.genotype,
                    self.consensus_quality,
                    self.snp_quality,
                    self.mapping_quality,
                    self.coverage ) ) )


cdef class IndexedReads:
    """index a Sam/BAM-file by query name.

    The index is kept in memory and can be substantial.

    By default, the file is re-openend to avoid conflicts if multiple
    operators work on the same file. Set *multiple_iterators* = False
    to not re-open *samfile*.
    """

    def __init__(self, AlignmentFile samfile, int multiple_iterators=True):

        # makes sure that samfile stays alive as long as this
        # object is alive.
        self.samfile = samfile

        assert samfile.is_bam, "can only IndexReads on bam files"

        # multiple_iterators the file - note that this makes the iterator
        # slow and causes pileup to slow down significantly.
        if multiple_iterators:
            self.htsfile = hts_open(samfile._filename, 'r')
            assert self.htsfile != NULL
            # read header - required for accurate positioning
            self.header = sam_hdr_read(self.htsfile)
            self.owns_samfile = True
        else:
            self.htsfile = self.samfile.htsfile
            self.header = self.samfile.header
            self.owns_samfile = False

    def build(self):
        '''build index.'''

        self.index = collections.defaultdict(list)

        # this method will start indexing from the current file
        # position if you decide
        cdef int ret = 1
        cdef bam1_t * b = <bam1_t*>calloc(1, sizeof( bam1_t))

        cdef uint64_t pos

        while ret > 0:
            pos = bgzf_tell(hts_get_bgzfp(self.htsfile))
            ret = sam_read1(self.htsfile,
                            self.samfile.header,
                            b)
            if ret > 0:
                qname = _charptr_to_str(pysam_bam_get_qname(b))
                self.index[qname].append(pos)

        bam_destroy1(b)

    def find(self, query_name):
        '''find *query_name* in index.

        Returns an iterator over all reads with query_name.

        Raise a KeyError if the *query_name* is not in the index.
        '''
        if query_name in self.index:
            return IteratorRowSelection(
                self.samfile,
                self.index[query_name],
                multiple_iterators = False)
        else:
            raise KeyError("read %s not found" % query_name)

    def __dealloc__(self):
        if self.owns_samfile:
            hts_close(self.htsfile)
            bam_hdr_destroy(self.header)

cpdef set_verbosity(int verbosity):
    u"""Set htslib's hts_verbose global variable to the specified value.
    """
    return hts_set_verbosity(verbosity)

cpdef get_verbosity():
    u"""Return the value of htslib's hts_verbose global variable.
    """
    return hts_get_verbosity()

__all__ = ["AlignmentFile",
           "IteratorRow",
           "IteratorColumn",
           "AlignedSegment",
           "PileupColumn",
           "PileupRead",
           "IndexedReads",
           "toQualityString",
           "fromQualityString",
           "get_verbosity",
           "set_verbosity"]
           # "IteratorSNPCalls",
           # "SNPCaller",
           # "IndelCaller",
           # "IteratorIndelCalls",



