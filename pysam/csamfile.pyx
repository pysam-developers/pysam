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
from cpython cimport PyErr_SetString, \
    PyBytes_Check, \
    PyUnicode_Check, \
    PyBytes_FromStringAndSize

from cpython.version cimport PY_MAJOR_VERSION

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
    u"""convert string or unicode object to bytes, assuming ascii encoding.
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
cdef class AlignedRead
cdef makeAlignedRead(bam1_t * src):
    '''enter src into AlignedRead.'''
    cdef AlignedRead dest = AlignedRead.__new__(AlignedRead)
    dest._delegate = bam_dup1(src)
    return dest

cdef class PileupProxy
cdef makePileupProxy(bam_pileup1_t ** plp, int tid, int pos, int n):
     cdef PileupProxy dest = PileupProxy.__new__(PileupProxy)
     dest.plp = plp
     dest.tid = tid
     dest.pos = pos
     dest.n = n
     return dest

cdef class PileupRead
cdef makePileupRead( bam_pileup1_t * src ):
    '''fill a  PileupRead object from a bam_pileup1_t * object.'''
    cdef PileupRead dest = PileupRead.__new__(PileupRead)
    dest._alignment = makeAlignedRead( src.b )
    dest._qpos = src.qpos
    dest._indel = src.indel
    dest._level = src.level
    dest._is_del = src.is_del
    dest._is_head = src.is_head
    dest._is_tail = src.is_tail
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


class PileupColumn(object):
    '''A pileup column. A pileup column contains
    all the reads that map to a certain target base.

    tid
        chromosome ID as is defined in the header
    pos
        the target base coordinate (0-based)
    n
        number of reads mapping to this column
    pileups
        list of reads (:class:`pysam.PileupRead`) aligned to this column
    '''
    def __str__(self):
        return "\t".join( map(str, (self.tid, self.pos, self.n))) +\
            "\n" + "\n".join( map(str, self.pileups) )


# valid types for sam headers
VALID_HEADER_TYPES = {"HD" : dict,
                      "SQ" : list,
                      "RG" : list,
                      "PG" : list,
                      "CO" : list}

# order of records within sam headers
VALID_HEADERS = ("HD", "SQ", "RG", "PG", "CO")

# type conversions within sam header records
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


cdef class Samfile:
    '''*(filename, mode=None, template = None,
         referencenames = None, referencelengths = None,
         text = NULL, header = None,
         add_sq_text = False, check_header = True,
         check_sq = True )*

    A :term:`SAM`/:term:`BAM` formatted file. The file is
    automatically opened.

    *mode* should be ``r`` for reading or ``w`` for writing. The
    default is text mode (:term:`SAM`). For binary (:term:`BAM`) I/O
    you should append ``b`` for compressed or ``u`` for uncompressed
    :term:`BAM` output.  Use ``h`` to output header information in
    text (:term:`TAM`) mode.

    If ``b`` is present, it must immediately follow ``r`` or ``w``.
    Valid modes are ``r``, ``w``, ``wh``, ``rb``, ``wb`` and
    ``wbu``. For instance, to open a :term:`BAM` formatted file for
    reading, type::

        f = pysam.Samfile('ex1.bam','rb')

    If mode is not specified, we will try to auto-detect in the order
    'rb', 'r', thus both the following should work::

        f1 = pysam.Samfile('ex1.bam' )
        f2 = pysam.Samfile('ex1.sam' )

    If an index for a BAM file exists (.bai), it will be opened
    automatically. Without an index random access to reads via
    :meth:`fetch` and :meth:`pileup` is disabled.

    For writing, the header of a :term:`SAM` file/:term:`BAM` file can
    be constituted from several sources (see also the samtools format
    specification):

        1. If *template* is given, the header is copied from a another
           *Samfile* (*template* must be of type *Samfile*).

        2. If *header* is given, the header is built from a
           multi-level dictionary. The first level are the four types
           ('HD', 'SQ', ...). The second level are a list of lines,
           with each line being a list of tag-value pairs. The header
           is constructed first from all the defined fields, followed
           by user tags in alphabetical order.

        3. If *text* is given, new header text is copied from raw
           text.

        4. The names (*referencenames*) and lengths
           (*referencelengths*) are supplied directly as lists.  By
           default, 'SQ' and 'LN' tags will be added to the header
           text. This option can be changed by unsetting the flag
           *add_sq_text*.

    By default, if a file is opened in mode 'r', it is checked
    for a valid header (*check_header* = True) and a definition of
    chromosome names (*check_sq* = True).

    '''

    def __cinit__(self, *args, **kwargs ):
        self.htsfile = NULL
        self._filename = None
        self.isbam = False
        self.isstream = False
        self._open(*args, **kwargs)

        # allocate memory for iterator
        self.b = <bam1_t*>calloc(1, sizeof(bam1_t))

    def _isOpen( self ):
        '''return true if htsfile has been opened.'''
        return self.htsfile != NULL

    def _hasIndex( self ):
        '''return true if htsfile has an existing (and opened) index.'''
        return self.index != NULL

    def _open(self,
              filename,
              mode=None,
              Samfile template=None,
              referencenames=None,
              referencelengths=None,
              text=None,
              header=None,
              port=None,
              add_sq_text=True,
              check_header=True,
              check_sq=True):
        '''open a sam/bam file.

        If _open is called on an existing bamfile, the current file will be
        closed and a new file will be opened.
        '''

        # read mode autodetection
        if mode is None:
            try:
                self._open(filename, 'rb',
                           template=template,
                           referencenames=referencenames,
                           referencelengths=referencelengths,
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
                       referencenames=referencenames,
                       referencelengths=referencelengths,
                       text=text,
                       header=header,
                       port=port,
                       check_header=check_header,
                       check_sq=check_sq)
            return

        assert mode in ( "r","w","rb","wb", "wh", "wbu", "rU" ), \
            "invalid file opening mode `%s`" % mode

        # close a previously opened file
        if self.htsfile != NULL:
            self.close()

        cdef bytes bmode = mode.encode('ascii')
        self._filename = filename = _encodeFilename(filename)
        self.isstream = filename == b"-"

        self.isbam = len(mode) > 1 and mode[1] == 'b'

        self.isremote = filename.startswith(b"http:") or \
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
                assert referencenames and referencelengths, \
                    ("either supply options `template`, `header` "
                     "or  both `referencenames` and `referencelengths` "
                     "for writing")
                assert len(referencenames) == len(referencelengths), \
                    "unequal names and lengths of reference sequences"

                # allocate and fill header
                referencenames = [_forceBytes(ref) for ref in referencenames]
                self.header = bam_hdr_init()
                self.header.n_targets = len(referencenames)
                n = 0
                for x in referencenames:
                    n += len(x) + 1
                self.header.target_name = <char**>calloc(n, sizeof(char*))
                self.header.target_len = <uint32_t*>calloc(n, sizeof(uint32_t))
                for x from 0 <= x < self.header.n_targets:
                    self.header.target_len[x] = referencelengths[x]
                    name = referencenames[x]
                    self.header.target_name[x] = <char*>calloc(
                        len(name) + 1, sizeof(char))
                    strncpy(self.header.target_name[x], name, len(name))

                # Optionally, if there is no text, add a SAM
                # compatible header to output file.
                if text is None and add_sq_text:
                    text = []
                    for x from 0 <= x < self.header.n_targets:
                        text.append("@SQ\tSN:%s\tLN:%s\n" % \
                                    (_forceStr(referencenames[x]), 
                                     referencelengths[x]))
                    text = ''.join(text)

                if text != None:
                    # copy without \0
                    text = _forceBytes(text)
                    ctext = text
                    self.header.l_text = strlen(ctext)
                    self.header.text = <char*>calloc(
                        strlen(ctext), sizeof(char))
                    memcpy(self.header.text, ctext, strlen(ctext))

            # open file. Header gets written to file at the same time for bam files
            # and sam files (in the latter case, the mode needs to be wh)
            self.htsfile = hts_open(filename, bmode)
            
            # for compatibility - "w" writes sam file without header
            if self.isbam or "h" in mode:
                # write header to htsfile
                sam_hdr_write(self.htsfile, self.header)
                
        elif mode[0] == "r":
            # open file for reading
            if (filename != b"-"
                and not self.isremote
                and not os.path.exists(filename)):
                raise IOError("file `%s` not found" % filename)

            # try to detect errors
            self.htsfile = hts_open(filename, bmode)
            if self.htsfile == NULL:
                raise ValueError(
                    "could not open file (mode='%s') - "
                    "is it SAM/BAM format?" % mode)

            # get file pointer
            # TODO: this is specific to BAM files
            #       refactor to make generalizable
            self.fp = self.htsfile.fp.bgzf

            # bam files require a valid header
            if self.isbam:
                self.header = sam_hdr_read(self.htsfile)
                if self.header == NULL:
                    raise ValueError(
                        "file does not have valid header (mode='%s') "
                        "- is it BAM format?" % mode )
            else:
                # in sam files it is optional (htsfile full of unmapped reads)
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
        if mode[0] == "r" and self.isbam:

            if not self.isremote:
                if not os.path.exists(filename + b".bai") \
                        and not os.path.exists( filename[:-4] + b".bai"):
                    self.index = NULL
                else:
                    # returns NULL if there is no index or index could not be opened
                    self.index = hts_idx_load(filename, HTS_FMT_BAI)
                    if self.index == NULL:
                        raise IOError("error while opening index `%s` " % filename )
            else:
                self.index = hts_idx_load(filename, HTS_FMT_BAI)
                if self.index == NULL:
                    warnings.warn("unable to open index for `%s` " % filename)

            if not self.isstream:
                self.start_offset = bgzf_tell(self.fp)

    def gettid( self, reference ):
        '''
        convert :term:`reference` name into numerical :term:`tid`

        returns -1 if reference is not known.
        '''
        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )
        reference = _forceBytes(reference)
        return bam_name2id(self.header, reference)

    def getrname( self, tid ):
        '''
        convert numerical :term:`tid` into :term:`reference` name.'''
        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )
        if not 0 <= tid < self.header.n_targets:
            raise ValueError("tid %i out of range 0<=tid<%i" % 
                             (tid, self.header.n_targets ) )
        return _charptr_to_str(self.header.target_name[tid])

    cdef char * _getrname( self, int tid ): # TODO unused
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
        '''
        parse region information.

        raise ValueError for for invalid regions.

        returns a tuple of flag, tid, start and end. Flag indicates
        whether some coordinates were supplied.

        Note that regions are 1-based, while start,end are python coordinates.
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
                raise ValueError( 'start out of range (%i)' % start )

        if end != None:
            try:
                rend = end
            except OverflowError:
                raise ValueError( 'end out of range (%i)' % end )

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
        '''reset file position to beginning of read section.'''
        return self.seek(self.start_offset, 0)

    def seek(self, uint64_t offset, int where = 0):
        '''
        move file pointer to position *offset*, see :meth:`pysam.Samfile.tell`.
        '''

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )
        if not self.isbam:
            raise NotImplementedError("seek only available in bam files")
        if self.isstream:
            raise OSError("seek no available in streams")

        return bgzf_seek(self.fp, offset, where)

    def tell(self):
        '''
        return current file position
        '''
        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )
        if not self.isbam:
            raise NotImplementedError("seek only available in bam files")

        return bgzf_tell(self.fp)

    def fetch(self,
              reference=None,
              start=None,
              end=None,
              region=None,
              tid=None,
              callback=None,
              until_eof=False,
              reopen=True):
        '''fetch aligned reads in a :term:`region` using 0-based indexing. The
        region is specified by :term:`reference`, *start* and
        *end*. Alternatively, a samtools :term:`region` string can be
        supplied.

        Without *reference* or *region* all mapped reads will be
        fetched. The reads will be returned ordered by reference
        sequence, which will not necessarily be the order within the
        file.

        If *until_eof* is given, all reads from the current file
        position will be returned in order as they are within the
        file. Using this option will also fetch unmapped reads.

        If *reopen* is set to true, the iterator returned will receive
        its own filehandle to the htsfile effectively opening its own
        copy of the file. The default behaviour is to re-open in order
        to safely work with multiple concurrent iterators on the same
        file. Re-opening a htsfile creates some overhead, so when
        using many calls to fetch() *reopen* can be set to False to
        gain some speed. Also, the tell() method will only work if
        *reopen* is set to False.

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
        if self.isstream:
            reopen = False

        if self.isbam:
            if not until_eof and not self._hasIndex() and not self.isremote:
                raise ValueError( "fetch called on bamfile without index" )

            if has_coord:
                return IteratorRowRegion(self, rtid, rstart, rend, 
                                         reopen=reopen)
            else:
                if until_eof:
                    return IteratorRowAll(self, reopen=reopen)
                else:
                    # AH: check - reason why no reopen for AllRefs?
                    return IteratorRowAllRefs(self) # , reopen=reopen )
        else:
            if has_coord:
                raise ValueError ("fetching by region is not available for sam files")

            if callback:
                raise NotImplementedError("callback not implemented yet")

            if self.header == NULL:
                raise ValueError("fetch called for htsfile without header")

            # check if targets are defined
            # give warning, sam_read1 segfaults
            if self.header.n_targets == 0:
                warnings.warn("fetch called for htsfile without header")
                
            return IteratorRowAll(self, reopen=reopen)

    def head(self, n):
        '''
        return iterator over the first n alignments. 

        This is useful for inspecting the bam-file.
        '''
        return IteratorRowHead(self, n)

    def mate(self,
             AlignedRead read):
        '''return the mate of :class:`AlignedRead` *read*.

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
            raise ValueError("read %s: is unpaired" % (read.qname))
        if flag & BAM_FMUNMAP != 0:
            raise ValueError("mate %s: is unmapped" % (read.qname))

        # xor flags to get the other mate
        cdef int x = BAM_FREAD1 + BAM_FREAD2
        flag = (flag ^ x) & x

        # the following code is not using the C API and
        # could thus be made much quicker
        for mate in self.fetch(
                read._delegate.core.mpos,
                read._delegate.core.mpos + 1,
                tid=read._delegate.core.mtid):
            if mate.flag & flag != 0 and \
               mate.qname == read.qname:
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

        Note that a :term:`TAM` file does not allow random access. If
        *region* or *reference* are given, an exception is raised.
        '''
        cdef AlignedRead read
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
        '''
        perform a :term:`pileup` within a :term:`region`. The region is specified by
        :term:`reference`, *start* and *end* (using 0-based indexing).
        Alternatively, a samtools *region* string can be supplied.

        Without *reference* or *region* all reads will be used for the pileup. The reads will be returned
        ordered by :term:`reference` sequence, which will not necessarily be the order within the file.

        The method returns an iterator of type :class:`pysam.IteratorColumn` unless
        a *callback is provided. If a *callback* is given, the callback will be executed
        for each column within the :term:`region`.

        Note that :term:`SAM` formatted files do not allow random access.
        In these files, if a *region* or *reference* are given an exception is raised.

        Optional *kwargs* to the iterator:

        stepper
           The stepper controlls how the iterator advances.
           Possible options for the stepper are

           ``all``
              use all reads for pileup.
           ``samtools``
              same filter and read processing as in :term:`csamtools` pileup

        fastafile
           A :class:`FastaFile` object

         mask
           Skip all reads with bits set in mask if mask=True.

         max_depth
           Maximum read depth permitted. The default limit is *8000*.

         truncate
           By default, the samtools pileup engine outputs all reads overlapping a region (see note below).
           If truncate is True and a region is given, only output columns in the exact region
           specificied.

        .. note::

            *all* reads which overlap the region are returned. The first base returned will be the
            first base of the first read *not* necessarily the first base of the region used in the query.

        '''
        cdef int rtid, rstart, rend, has_coord

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        has_coord, rtid, rstart, rend = self._parseRegion(
            reference, start, end, region )

        if self.isbam:
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

    def close( self ):
        '''
        closes the :class:`pysam.Samfile`.'''
        if self.htsfile != NULL:
            hts_close(self.htsfile)
            hts_idx_destroy(self.index);
            self.htsfile = NULL

    def __dealloc__( self ):
        # remember: dealloc cannot call other methods
        # note: no doc string
        # note: __del__ is not called.
        self.close()
        bam_destroy1(self.b)
        if self.header != NULL:
            bam_hdr_destroy(self.header)
            
    cpdef int write( self, AlignedRead read ) except -1:
        '''
        write a single :class:`pysam.AlignedRead` to disk.

        returns the number of bytes written.
        '''
        if not self._isOpen():
            return 0

        x = sam_write1(self.htsfile,
                       self.header,
                       read._delegate)

        return x

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
        '''number of :term:`filename` associated with this object.'''
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
                t.append( _charptr_to_str(self.header.target_name[x]) )
            return tuple(t)

    property lengths:
        """tuple of the lengths of the :term:`reference` sequences. The lengths are in the same order as
        :attr:`pysam.Samfile.references`
        """
        def __get__(self):
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
            t = []
            for x from 0 <= x < self.header.n_targets:
                t.append( self.header.target_len[x] )
            return tuple(t)

    property mapped:
        """total number of mapped alignments in file.
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
        if not self.isbam:
            raise AttributeError("Samfile.mapped only available in bam files")
        if self.index == NULL:
            raise ValueError("mapping information not recorded in index "
                                 "or index not available")


    property unmapped:
        """total number of unmapped reads in file.
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
        """total number of reads without coordinates
        """
        def __get__(self):
            self._checkIndex()
            return hts_idx_get_n_no_coor(self.index)

    property text:
        '''full contents of the :term:`sam file` header as a string.'''
        def __get__(self):
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
            return from_string_and_size(self.header.text, self.header.l_text)

    property header:
        '''header information within the :term:`sam file`. The records and fields are returned as
        a two-level dictionary.
        '''
        def __get__(self):
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )

            result = {}
            
            if self.header.text != NULL:
                # convert to python string (note: call self.text to create 0-terminated string)
                t = self.text
                for line in t.split("\n"):
                    if not line.strip(): continue
                    assert line.startswith("@"), "header line without '@': '%s'" % line
                    fields = line[1:].split("\t")
                    record = fields[0]
                    assert record in VALID_HEADER_TYPES, "header line with invalid type '%s': '%s'" % (record, line)

                    # treat comments
                    if record == "CO":
                        if record not in result: result[record] = []
                        result[record].append( "\t".join( fields[1:] ) )
                        continue
                    # the following is clumsy as generators do not work?
                    x = {}
                    for field in fields[1:]:
                        if ":" not in field: 
                            raise ValueError("malformatted header: no ':' in field" )
                        key, value = field.split(":",1)
                        # uppercase keys must be valid
                        # lowercase are permitted for user fields
                        if key in VALID_HEADER_FIELDS[record]:
                            x[key] = VALID_HEADER_FIELDS[record][key](value)
                        elif not key.isupper():
                            x[key] = value
                        else:
                            raise ValueError( "unknown field code '%s' in record '%s'" % (key, record) )

                    if VALID_HEADER_TYPES[record] == dict:
                        if record in result:
                            raise ValueError( "multiple '%s' lines are not permitted" % record )
                        result[record] = x
                    elif VALID_HEADER_TYPES[record] == list:
                        if record not in result: result[record] = []
                        result[record].append( x )

                # if there are no SQ lines in the header, add the reference names
                # from the information in the bam file.
                # Background: c-samtools keeps the textual part of the header separate from
                # the list of reference names and lengths. Thus, if a header contains only 
                # SQ lines, the SQ information is not part of the textual header and thus
                # are missing from the output. See issue 84.
                if "SQ" not in result:
                    sq = []
                    for ref, length in zip( self.references, self.lengths ):
                        sq.append( {'LN': length, 'SN': ref } )
                    result["SQ"] = sq

            return result

    def _buildLine( self, fields, record ):
        '''build a header line from *fields* dictionary for *record*'''

        # TODO: add checking for field and sort order
        line = ["@%s" % record ]
        # comment
        if record == "CO":
            line.append( fields )
        # user tags
        elif record.islower():
            for key in sorted(fields):
                line.append( "%s:%s" % (key, str(fields[key])))
        # defined tags
        else:
            # write fields of the specification
            for key in VALID_HEADER_ORDER[record]:
                if key in fields:
                    line.append( "%s:%s" % (key, str(fields[key])))
            # write user fields
            for key in fields:
                if not key.isupper():
                    line.append( "%s:%s" % (key, str(fields[key])))

        return "\t".join( line )

    cdef bam_hdr_t * _buildHeader( self, new_header ):
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
            dest.target_name = <char**>calloc( dest.n_targets, sizeof(char*) )
            dest.target_len = <uint32_t*>calloc( dest.n_targets, sizeof(uint32_t) )

            for x from 0 <= x < dest.n_targets:
                seqname, seqlen = seqs[x]
                dest.target_name[x] = <char*>calloc( len( seqname ) + 1, sizeof(char) )
                bseqname = seqname.encode('ascii')
                strncpy( dest.target_name[x], bseqname, len(seqname) + 1 )
                dest.target_len[x] = seqlen

        return dest

    ###############################################################
    ###############################################################
    ###############################################################
    ## file-object like iterator access
    ## note: concurrent access will cause errors (see IteratorRow
    ## and reopen)
    ## Possible solutions: deprecate or open new file handle
    ###############################################################
    def __iter__(self):
        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        if not self.isbam and self.header.n_targets == 0:
            raise NotImplementedError(
                "can not iterate over samfile without header")
        return self

    cdef bam1_t * getCurrent( self ):
        return self.b

    cdef int cnext(self):
        '''
        cversion of iterator. Used by :class:`pysam.Samfile.IteratorColumn`.
        '''
        cdef int ret
        return sam_read1(self.htsfile,
                         self.header,
                         self.b)

    def __next__(self):
        """
        python version of next().
        """
        cdef int ret
        ret = sam_read1(self.htsfile, self.header, self.b)
        if (ret >= 0):
            return makeAlignedRead(self.b)
        else:
            raise StopIteration

##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
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

    The method :meth:`Samfile.fetch` returns an IteratorRow.

    .. note::
        It is usually not necessary to create an object of this class
        explicitely. It is returned as a result of call to a :meth:`Samfile.fetch`.

    '''

    def __init__(self, Samfile samfile, int reopen=True):
        
        if not samfile._isOpen():
            raise ValueError( "I/O operation on closed file" )

        # makes sure that samfile stays alive as long as the
        # iterator is alive
        self.samfile = samfile

        # reopen the file - note that this makes the iterator
        # slow and causes pileup to slow down significantly.
        if reopen:
            self.htsfile = hts_open(samfile._filename, 'r')
            assert self.htsfile != NULL
            # read header - required for accurate positioning
            sam_hdr_read(self.htsfile)
            self.owns_samfile = True
        else:
            self.htsfile = self.samfile.htsfile
            self.owns_samfile = False

        self.retval = 0

        self.b = bam_init1()

    def __dealloc__(self):
        bam_destroy1(self.b)
        if self.owns_samfile:
            hts_close(self.htsfile)

cdef class IteratorRowRegion(IteratorRow):
    """*(Samfile samfile, int tid, int beg, int end, int reopen = True )*

    iterate over mapped reads in a region.

    By default, the file is re-openend to avoid conflicts between
    multiple iterators working on the same file. Set *reopen* = False
    to not re-open *samfile*.

    The samtools iterators assume that the file
    position between iterations do not change.
    As a consequence, no two iterators can work
    on the same file. To permit this, each iterator
    creates its own file handle by re-opening the
    file.

    Note that the index will be shared between
    samfile and the iterator.

    .. note::
        It is usually not necessary to create an object of this class
        explicitely. It is returned as a result of call to a :meth:`Samfile.fetch`.

    """

    def __init__(self, Samfile samfile,
                 int tid, int beg, int end,
                 int reopen=True):

        IteratorRow.__init__(self, samfile, reopen=reopen)

        if not samfile._hasIndex():
            raise ValueError( "no index available for iteration" )

        self.iter = sam_itr_queryi(
            self.samfile.index,
            tid,
            beg,
            end)
    
    def __iter__(self):
        return self

    cdef bam1_t * getCurrent( self ):
        return self.b

    cdef int cnext(self):
        '''cversion of iterator. Used by IteratorColumn'''
        self.retval = hts_itr_next(self.htsfile.fp.bgzf,
                                   self.iter,
                                   self.b,
                                   NULL)

    def __next__(self):
        """python version of next().
        """
        self.cnext()
        if self.retval < 0:
            raise StopIteration
        return makeAlignedRead(self.b)

    def __dealloc__(self):
        hts_itr_destroy(self.iter)

cdef class IteratorRowHead(IteratorRow):
    """*(Samfile samfile, n, int reopen = True)*

    iterate over first n reads in *samfile*

    By default, the file is re-openend to avoid conflicts between
    multiple iterators working on the same file. Set *reopen* = False
    to not re-open *samfile*.

    .. note::
        It is usually not necessary to create an object of this class
        explicitely. It is returned as a result of call to a :meth:`Samfile.head`.
        

    """

    def __init__(self, Samfile samfile, int n, int reopen=True):

        IteratorRow.__init__(self, samfile, reopen=reopen)

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

        cdef int ret
        ret = sam_read1(self.htsfile,
                        self.samfile.header, self.b)
        if (ret >= 0):
            self.current_row += 1
            return makeAlignedRead( self.b )
        else:
            raise StopIteration


cdef class IteratorRowAll(IteratorRow):
    """*(Samfile samfile, int reopen = True)*

    iterate over all reads in *samfile*

    By default, the file is re-openend to avoid conflicts between
    multiple iterators working on the same file. Set *reopen* = False
    to not re-open *samfile*.

    .. note::
        It is usually not necessary to create an object of this class
        explicitely. It is returned as a result of call to a :meth:`Samfile.fetch`.
        

    """

    def __init__(self, Samfile samfile, int reopen = True):

        IteratorRow.__init__(self, samfile, reopen=reopen)

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
        cdef int ret
        ret = sam_read1(self.htsfile,
                        self.samfile.header,
                        self.b)
        if (ret >= 0):
            return makeAlignedRead(self.b)
        else:
            raise StopIteration


cdef class IteratorRowAllRefs(IteratorRow):
    """iterates over all mapped reads by chaining iterators over each reference

    .. note::
        It is usually not necessary to create an object of this class
        explicitely. It is returned as a result of call to a :meth:`Samfile.fetch`.
    """

    def __init__(self, Samfile samfile, reopen=True):

        IteratorRow.__init__(self, samfile, reopen=reopen)

        if not samfile._hasIndex():
            raise ValueError("no index available for fetch")

        self.tid = -1

    def nextiter(self):
        self.rowiter = IteratorRowRegion(self.samfile,
                                         self.tid,
                                         0,
                                         1<<29)

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
            if self.rowiter.retval>0:
                return makeAlignedRead(self.rowiter.b)

            self.tid += 1

            # Otherwise, proceed to next reference or stop
            if self.tid < self.samfile.nreferences:
                self.nextiter()
            else:
                raise StopIteration


cdef class IteratorRowSelection(IteratorRow):
    """*(Samfile samfile)*

    iterate over reads in *samfile* at a given list of file positions.

    .. note::
        It is usually not necessary to create an object of this class
        explicitely. It is returned as a result of call to a :meth:`Samfile.fetch`.
    """

    def __init__(self, Samfile samfile, positions, int reopen=True):

        IteratorRow.__init__(self, samfile, reopen=reopen)

        self.positions = positions
        self.current_pos = 0

        self.fp = self.htsfile.fp.bgzf

    def __iter__(self):
        return self

    cdef bam1_t * getCurrent( self ):
        return self.b

    cdef int cnext(self):
        '''cversion of iterator'''

        # end iteration if out of positions
        if self.current_pos >= len(self.positions): return -1

        bgzf_seek(self.fp,
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
            return makeAlignedRead(self.b)
        else:
            raise StopIteration


cdef int __advance_all(void *data, bam1_t *b):
    '''advance without any read filtering.
    '''
    cdef __iterdata * d
    d = <__iterdata*>data
    return sam_itr_next(d.htsfile, d.iter, b)


# cdef int __advance_snpcalls( void * data, bam1_t * b ):
#     '''advance using same filter and read processing as in
#     the samtools pileup.
#     '''
#     cdef __iterdata * d
#     d = <__iterdata*>data

#     cdef int ret = bam_iter_read( d.samfile.x.bam, d.iter, b )
#     cdef int skip = 0
#     cdef int q
#     cdef int is_cns = 1
#     cdef int is_nobaq = 0
#     cdef int capQ_thres = 0

#     # reload sequence
#     if d.fastafile != NULL and b.core.tid != d.tid:
#         if d.seq != NULL: free(d.seq)
#         d.tid = b.core.tid
#         d.seq = faidx_fetch_seq(d.fastafile,
#                                 d.samfile.header.target_name[d.tid],
#                                 0, max_pos,
#                                 &d.seq_len)
#         if d.seq == NULL:
#             raise ValueError( "reference sequence for '%s' (tid=%i) not found" % \
#                                   (d.samfile.header.target_name[d.tid],
#                                    d.tid))


#     while ret >= 0:

#         skip = 0

#         # realign read - changes base qualities
#         if d.seq != NULL and is_cns and not is_nobaq: 
#             bam_prob_realn( b, d.seq )

#         if d.seq != NULL and capQ_thres > 10:
#             q = bam_cap_mapQ(b, d.seq, capQ_thres)
#             if q < 0: skip = 1
#             elif b.core.qual > q: b.core.qual = q
#         if b.core.flag & BAM_FUNMAP: skip = 1
#         elif b.core.flag & 1 and not b.core.flag & 2: skip = 1

#         if not skip: break
#         # additional filters

#         ret = bam_iter_read( d.samfile.x.bam, d.iter, b )

#     return ret

cdef class IteratorColumn:
    '''abstract base class for iterators over columns.

    IteratorColumn objects wrap the pileup functionality of samtools.

    For reasons of efficiency, the iterator points to the current
    pileup buffer. The pileup buffer is updated at every iteration.
    This might cause some unexpected behavious. For example,
    consider the conversion to a list::

       f = Samfile("file.bam", "rb")
       result = list( f.pileup() )

    Here, ``result`` will contain ``n`` objects of type :class:`PileupProxy` for ``n`` columns,
    but each object in ``result`` will contain the same information.

    The desired behaviour can be achieved by list comprehension::

       result = [ x.pileups() for x in f.pileup() ]

    ``result`` will be a list of ``n`` lists of objects of type :class:`PileupRead`.

    If the iterator is associated with a :class:`Fastafile` using the :meth:`addReference`
    method, then the iterator will export the current sequence via the methods :meth:`getSequence`
    and :meth:`seq_len`.

    Optional kwargs to the iterator

    stepper
       The stepper controls how the iterator advances.
       Possible options for the stepper are

       all
           use all reads for pileup.
       samtools
           same filter and read processing as in :term:`csamtools` pileup

       The default is to use "all" if no stepper is given.

    fastafile
       A :class:`FastaFile` object
    mask
       Skip all reads with bits set in mask.
    max_depth
       maximum read depth. The default is 8000.
    '''

    def __cinit__( self, Samfile samfile, **kwargs ):
        self.samfile = samfile
        # TODO
        # self.mask = kwargs.get("mask", BAM_DEF_MASK )
        self.fastafile = kwargs.get( "fastafile", None )
        self.stepper = kwargs.get( "stepper", None )
        self.max_depth = kwargs.get( "max_depth", 8000 )
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

    def addReference( self, Fastafile fastafile ):
       '''
       add reference sequences in *fastafile* to iterator.'''
       self.fastafile = fastafile
       if self.iterdata.seq != NULL: free(self.iterdata.seq)
       self.iterdata.tid = -1
       self.iterdata.fastafile = self.fastafile.fastafile

    def hasReference( self ):
        '''
        return true if iterator is associated with a reference'''
        return self.fastafile

    cdef setMask( self, mask ):
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
                            int reopen = 0 ):
        '''setup the iterator structure'''

        self.iter = IteratorRowRegion(self.samfile, tid, start, end, reopen)
        self.iterdata.htsfile = self.samfile.htsfile
        self.iterdata.iter = self.iter.iter
        self.iterdata.seq = NULL
        self.iterdata.tid = -1

        if self.fastafile != None:
            self.iterdata.fastafile = self.fastafile.fastafile
        else:
            self.iterdata.fastafile = NULL

        if self.stepper == None or self.stepper == "all":
            self.pileup_iter = bam_plp_init(
                <bam_plp_auto_f>&__advance_all,
                &self.iterdata)
        # elif self.stepper == "samtools":
        #     self.pileup_iter = bam_plp_init(&__advance_snpcalls, &self.iterdata)
        else:
            raise ValueError(
                "unknown stepper option `%s` in IteratorColumn" % self.stepper)

        if self.max_depth:
            bam_plp_set_maxcnt( self.pileup_iter, self.max_depth )

        # bam_plp_set_mask( self.pileup_iter, self.mask )

    cdef reset( self, tid, start, end ):
        '''reset iterator position.

        This permits using the iterator multiple times without
        having to incur the full set-up costs.
        '''
        self.iter = IteratorRowRegion( self.samfile, tid, start, end, reopen = 0 )
        self.iterdata.iter = self.iter.iter

        # invalidate sequence if different tid
        if self.tid != tid:
            if self.iterdata.seq != NULL: free( self.iterdata.seq )
            self.iterdata.seq = NULL
            self.iterdata.tid = -1

        # self.pileup_iter = bam_plp_init( &__advancepileup, &self.iterdata )
        bam_plp_reset(self.pileup_iter)

    def __dealloc__(self):
        # reset in order to avoid memory leak messages for iterators 
        # that have not been fully consumed
        if self.pileup_iter != <bam_plp_t>NULL:
            bam_plp_reset(self.pileup_iter)
            bam_plp_destroy(self.pileup_iter)
            self.pileup_iter = <bam_plp_t>NULL
            self.plp = <bam_pileup1_t*>NULL

        if self.iterdata.seq != NULL:
            free(self.iterdata.seq)
            self.iterdata.seq = NULL

cdef class IteratorColumnRegion(IteratorColumn):
    '''iterates over a region only.
    '''
    def __cinit__(self, Samfile samfile,
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

            return makePileupProxy(&self.plp,
                                   self.tid,
                                   self.pos,
                                   self.n_plp)

cdef class IteratorColumnAllRefs(IteratorColumn):
    """iterates over all columns by chaining iterators over each reference
    """

    def __cinit__(self,
                  Samfile samfile,
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
                return makePileupProxy(&self.plp,
                                       self.tid,
                                       self.pos,
                                       self.n_plp)
                
            # otherwise, proceed to next reference or stop
            self.tid += 1
            if self.tid < self.samfile.nreferences:
                self.setupIteratorData(self.tid, 0, max_pos, 0)
            else:
                raise StopIteration

##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
cdef inline int32_t query_start(bam1_t *src) except -1:
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

##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
cdef inline int32_t query_end(bam1_t *src) except -1:
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


cdef inline object get_seq_range(bam1_t *src, uint32_t start, uint32_t end):
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

    return seq


cdef inline object get_qual_range(bam1_t *src, uint32_t start, uint32_t end):
    cdef uint8_t * p
    cdef uint32_t k
    cdef char * q

    p = pysam_bam_get_qual(src)
    if p[0] == 0xff:
        return None

    qual = PyBytes_FromStringAndSize(NULL, end - start)
    q    = <char*>qual

    for k from start <= k < end:
        ## equivalent to t[i] + 33 (see bam.c)
        q[k-start] = p[k] + 33

    return qual

cdef inline uint8_t get_type_code(value, value_type = None):
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
        value_type = _forceBytes( value_type )
        _char_type = value_type
        type_code = (<uint8_t*>_char_type)[0]

    return type_code

cdef inline convert_python_tag(pytag, value, fmts, args):
    
    if not type(pytag) is bytes:
        pytag = pytag.encode('ascii')
    t = type(value)

    if t is tuple or t is list:
        # binary tags - treat separately
        pytype = 'B'
        # get data type - first value determines type
        if type(value[0]) is float:
            datafmt, datatype = "f", "f"
        else:
            mi, ma = min(value), max(value)
            absmax = max( abs(mi), abs(ma) )
            # signed ints
            if mi < 0: 
                if mi >= -127: datafmt, datatype = "b", 'c'
                elif mi >= -32767: datafmt, datatype = "h", 's'
                elif absmax < -2147483648: raise ValueError( "integer %i out of range of BAM/SAM specification" % value )
                else: datafmt, datatype = "i", 'i'

            # unsigned ints
            else:
                if absmax <= 255: datafmt, datatype = "B", 'C'
                elif absmax <= 65535: datafmt, datatype = "H", 'S'
                elif absmax > 4294967295: raise ValueError( "integer %i out of range of BAM/SAM specification" % value )
                else: datafmt, datatype = "I", 'I'

        datafmt = "2sccI%i%s" % (len(value), datafmt)
        args.extend( [pytag[:2], 
                      pytype.encode('ascii'),
                      datatype.encode('ascii'),
                      len(value)] + list(value) )
        fmts.append( datafmt )
        return

    if t is float:
        fmt, pytype = "2scf", 'f'
    elif t is int:
        # negative values
        if value < 0:
            if value >= -127: fmt, pytype = "2scb", 'c'
            elif value >= -32767: fmt, pytype = "2sch", 's'
            elif value < -2147483648: raise ValueError( "integer %i out of range of BAM/SAM specification" % value )
            else: fmt, pytype = "2sci", 'i'
        # positive values
        else:
            if value <= 255: fmt, pytype = "2scB", 'C'
            elif value <= 65535: fmt, pytype = "2scH", 'S'
            elif value > 4294967295: raise ValueError( "integer %i out of range of BAM/SAM specification" % value )
            else: fmt, pytype = "2scI", 'I'
    else:
        # Note: hex strings (H) are not supported yet
        if t is not bytes:
            value = value.encode('ascii')
        if len(value) == 1:
            fmt, pytype = "2scc", 'A'
        else:
            fmt, pytype = "2sc%is" % (len(value)+1), 'Z'

    args.extend( [pytag[:2],
                  pytype.encode('ascii'),
                  value ] )

    fmts.append( fmt )
    
###########################################################
###########################################################
###########################################################
cdef class AlignedRead:
    '''
    Class representing an aligned read. See the SAM format specification for
    the meaning of fields (http://samtools.sourceforge.net/).

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

    In Python 3, the fields containing sequence and quality
    (seq, query, qual and qqual) data are of type bytes. Other
    string data, such as the qname field and strings in the
    tags tuple, is represented as unicode strings. On assignment,
    both bytes and unicode objects are allowed, but unicode strings
    must contain only ASCII characters.
    '''

    # Now only called when instances are created from Python
    def __init__(self):
        # see bam_init1
        self._delegate = <bam1_t*>calloc( 1, sizeof( bam1_t) )
        # allocate some memory
        # If size is 0, calloc does not return a pointer that can be passed to free()
        # so allocate 40 bytes for a new read
        self._delegate.m_data = 40
        self._delegate.data = <uint8_t *>calloc( self._delegate.m_data, 1 )
        self._delegate.l_data = 0

    def __dealloc__(self):
        bam_destroy1(self._delegate)

    def __str__(self):
        """return string representation of alignment.

        The representation is an approximate :term:`sam` format.

        An aligned read might not be associated with a :term:`Samfile`.
        As a result :term:`tid` is shown instead of the reference name.

        Similarly, the tags field is returned in its parsed state.
        """
        # sam-parsing is done in sam.c/bam_format1_core which
        # requires a valid header.
        if sys.version_info[0] < 3:
            seq = self.seq
            qual = self.qual
        else:
            seq = self.seq.decode('ascii')
            qual = self.qual.decode('ascii')
        return "\t".join(map(str, (self.qname,
                                   self.flag,
                                   self.rname,
                                   self.pos,
                                   self.mapq,
                                   self.cigar,
                                   self.mrnm,
                                   self.mpos,
                                   self.rlen,
                                   seq,
                                   qual,
                                   self.tags )))

    def compare(self, AlignedRead other):
        '''return -1,0,1, if contents in this are binary <,=,> to *other*'''

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
        if t==o:
            return 0

        retval = memcmp(&t.core, &o.core, sizeof(bam1_core_t))

        if retval: return retval
        retval = (t.l_data > o.l_data) - (t.l_data < o.l_data) # cmp(t.l_data, o.l_data)
        if retval: return retval
        return memcmp(t.data, o.data, t.l_data)

    # Disabled so long as __cmp__ is a special method
    def __hash__(self):
        return _Py_HashPointer(<void *>self)

    def _convert_python_tag(self, pytag, value, fmts, args):

        if not type(pytag) is bytes:
            pytag = pytag.encode('ascii')
        t = type(value)

        if t is tuple or t is list:
            # binary tags - treat separately
            pytype = 'B'
            # get data type - first value determines type
            if type(value[0]) is float:
                datafmt, datatype = "f", "f"
            else:
                mi, ma = min(value), max(value)
                absmax = max( abs(mi), abs(ma) )
                # signed ints
                if mi < 0: 
                    if mi >= -127: datafmt, datatype = "b", 'c'
                    elif mi >= -32767: datafmt, datatype = "h", 's'
                    elif absmax < -2147483648: raise ValueError( "integer %i out of range of BAM/SAM specification" % value )
                    else: datafmt, datatype = "i", 'i'

                # unsigned ints
                else:
                    if absmax <= 255: datafmt, datatype = "B", 'C'
                    elif absmax <= 65535: datafmt, datatype = "H", 'S'
                    elif absmax > 4294967295: raise ValueError( "integer %i out of range of BAM/SAM specification" % value )
                    else: datafmt, datatype = "I", 'I'
                    
            datafmt = "2sccI%i%s" % (len(value), datafmt)
            args.extend( [pytag[:2], 
                          pytype.encode('ascii'),
                          datatype.encode('ascii'),
                          len(value)] + list(value) )
            fmts.append( datafmt )
            return

        if t is float:
            fmt, pytype = "2scf", 'f'
        elif t is int:
            # negative values
            if value < 0:
                if value >= -127: fmt, pytype = "2scb", 'c'
                elif value >= -32767: fmt, pytype = "2sch", 's'
                elif value < -2147483648: raise ValueError( "integer %i out of range of BAM/SAM specification" % value )
                else: fmt, pytype = "2sci", 'i'
            # positive values
            else:
                if value <= 255: fmt, pytype = "2scB", 'C'
                elif value <= 65535: fmt, pytype = "2scH", 'S'
                elif value > 4294967295: raise ValueError( "integer %i out of range of BAM/SAM specification" % value )
                else: fmt, pytype = "2scI", 'I'
        else:
            # Note: hex strings (H) are not supported yet
            if t is not bytes:
                value = value.encode('ascii')
            if len(value) == 1:
                fmt, pytype = "2scc", 'A'
            else:
                fmt, pytype = "2sc%is" % (len(value)+1), 'Z'

        args.extend( [pytag[:2],
                      pytype.encode('ascii'),
                      value ] )
        
        fmts.append( fmt )


    #######################################################################
    #######################################################################
    ## Basic properties
    #######################################################################
    property qname:
        """the query name (None if not present)"""
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            if pysam_get_l_qname(src) == 0:
                return None
            return _charptr_to_str(<char *>pysam_bam_get_qname(src))

        def __set__(self, qname ):
            if qname == None or len(qname) == 0: return
            qname = _forceBytes(qname)
            cdef bam1_t * src
            cdef int l
            cdef char * p

            src = self._delegate
            p = pysam_bam_get_qname( src )

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

    property cigar:
        """the :term:`cigar` alignment. The alignment
        is returned as a list of tuples of (operation, length). 

        If the alignment is not present, an empty list is
        returned.

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
            cdef op, l, cigar
            cdef int k
            cigar = []

            src = self._delegate
            if pysam_get_n_cigar(src) == 0:
                return cigar

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

    property cigarstring:
        '''the :term:`cigar` alignment as a string.
        
        The cigar string is a string of alternating integers
        and characters denoting the length and the type of
        an operation.

        .. note::
            The order length,operation is specified in the
            SAM format. It is different from the order of
            the :attr:`cigar` property.

        Returns the empty string if not present.

        To unset the cigarstring, assign None or the
        empty string.
        '''
        def __get__(self):
            c = self.cigar
            if c == None: return ""
            # reverse order
            else: return "".join([ "%i%c" % (y,CODE2CIGAR[x]) for x,y in c])
            
        def __set__(self, cigar):
            if cigar is None or len(cigar) == 0:
                self.cigar = []
            else:
                parts = CIGAR_REGEX.findall(cigar)
                # reverse order
                self.cigar = [(CIGAR2CODE[ord(y)], int(x)) for x,y in parts]

    property seq:
        """read sequence bases, including :term:`soft clipped` bases 
        (None if not present).

        In Python 3, this property is of type bytes and assigning a
        unicode string to it consisting of ASCII characters only will
        work, but is inefficient.

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

            return get_seq_range(src, 0, src.core.l_qseq)

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

    property qual:
        """read sequence base qualities, including :term:`soft
        clipped` bases (None if not present).

        In Python 3, this property is of type bytes and assigning a
        unicode string to it consisting of ASCII characters only will
        work, but is inefficient.

        Note that to set quality scores the sequence has to be set
        previously as this will determine the permitted length of
        the quality score array.

        This method raises a ValueError if the length of the 
        quality scores and the sequence are not the same.
        """
        def __get__(self):

            cdef bam1_t * src
            cdef char * q

            src = self._delegate

            if src.core.l_qseq == 0: return None

            return get_qual_range(src, 0, src.core.l_qseq)

        def __set__(self,qual):
            # note that space is already allocated via the sequences
            cdef bam1_t * src
            cdef uint8_t * p
            cdef char * q
            cdef int k

            src = self._delegate
            p = pysam_bam_get_qual(src)
            if qual == None or len(qual) == 0:
                # if absent - set to 0xff
                p[0] = 0xff
                return
            qual = _forceBytes(qual)
            cdef int l
            # convert to C string
            q = qual
            l = len(qual)
            if src.core.l_qseq != l:
                raise ValueError("quality and sequence mismatch: %i != %i" % (l, src.core.l_qseq))
            assert src.core.l_qseq == l
            for k from 0 <= k < l:
                p[k] = <uint8_t>q[k] - 33

    property query:
        """aligned portion of the read.

        This is a substring of :attr:`seq` that excludes flanking bases that were
        :term:`soft clipped` (None if not present). It is equal to ``seq[qstart:qend]``.

        In Python 3, this property is of type bytes. Assigning a
        unicode string to it consisting of ASCII characters only will
        work, but is inefficient.

        SAM/BAM files may include extra flanking bases that are
        not part of the alignment.  These bases may be the result of the
        Smith-Waterman or other algorithms, which may not require alignments
        that begin at the first residue or end at the last.  In addition,
        extra sequencing adapters, multiplex identifiers, and low-quality bases that
        were not considered for alignment may have been retained."""

        def __get__(self):
            cdef bam1_t * src
            cdef uint32_t start, end
            cdef char * s

            src = self._delegate

            if src.core.l_qseq == 0: return None

            start = query_start(src)
            end   = query_end(src)

            return get_seq_range(src, start, end)

    property qqual:
        """aligned query sequence quality values (None if not
        present). These are the quality values that correspond to :attr:`query`, that is,
        they exclude qualities of :term:`soft clipped` bases. This is equal to
        ``qual[qstart:qend]``.

        This property is read-only.

        In Python 3, this property is of type bytes."""
        def __get__(self):
            cdef bam1_t * src
            cdef uint32_t start, end

            src = self._delegate

            if src.core.l_qseq == 0: return None

            start = query_start(src)
            end   = query_end(src)

            return get_qual_range(src, start, end)

    property qstart:
        """start index of the aligned query portion of the sequence (0-based, inclusive).

        This the index of the first base in :attr:`seq` that is not soft-clipped.
        """
        def __get__(self):
            return query_start(self._delegate)

    property qend:
        """end index of the aligned query portion of the sequence (0-based, exclusive)"""
        def __get__(self):
            return query_end(self._delegate)

    property qlen:
        """length of the aligned query sequence.

        This is equal to :attr:`qend` - :attr:`qstart`"""
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            return query_end(src)-query_start(src)

    property tags:
        """the tags in the AUX field.

        This property permits convenience access to
        the tags. Changes it the returned list will
        not update the tags automatically. Instead,
        the following is required for adding a
        new tag::

            read.tags = read.tags + [("RG",0)]

        This method will happily write the same tag
        multiple times.
        """
        def __get__(self):
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
                    byte_size, nvalues, value = convertBinaryTagToList( s )
                    # 5 for 1 char and 1 int
                    s += 5 + ( nvalues * byte_size) - 1
                else:
                    raise KeyError("unknown type '%s'" % auxtype)

                s += 1

                result.append((_charptr_to_str(auxtag), value))

            return result

        def __set__(self, tags):
            cdef bam1_t * src
            cdef uint8_t * s
            cdef char * temp
            cdef int new_size = 0
            cdef int old_size
            src = self._delegate
            fmts, args = ["<"], []
            
            if tags is not None and len(tags) > 0:
                for pytag, value in tags:
                    convert_python_tag(pytag, value, fmts, args)
                fmt = "".join(fmts)
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

    cpdef setTag(self, tag, value, 
                 value_type = None, 
                 replace = True):
        '''
        Set optional field of alignment *tag* to *value*.  *value_type* may be specified,
        but if not the type will be inferred based on the Python type of *value*

        An existing value of the same tag will be overwritten unless
        *replace* is set to False.
        '''

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
        
        type_code = get_type_code(value, value_type)
        if type_code == 0:
            raise ValueError("can't guess type or invalid type code specified")

        # Not Endian-safe, but then again neither is samtools!
        if type_code == 'Z':
            value = _forceBytes( value )
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

        tag = _forceBytes( tag )
        if replace:
            existing_ptr = bam_aux_get(src, tag)
            if existing_ptr:
                bam_aux_del(src, existing_ptr)

        bam_aux_append(src,
                       tag,
                       type_code, 
                       value_size,
                       value_ptr)

    property flag:
        """properties flag"""
        def __get__(self):
            return pysam_get_flag(self._delegate)
        def __set__(self, flag):
            pysam_set_flag(self._delegate, flag)

    property rname:
        """
        :term:`target` ID

        DEPRECATED from pysam-0.4 - use tid in the future.
        The rname field caused a lot of confusion as it returns
        the :term:`target` ID instead of the reference sequence
        name.

        .. note::

            This field contains the index of the reference sequence
            in the sequence dictionary. To obtain the name
            of the reference sequence, use :meth:`pysam.Samfile.getrname()`

        """
        def __get__(self): return self._delegate.core.tid
        def __set__(self, tid): self._delegate.core.tid = tid

    property tid:
        """
        :term:`target` ID

        .. note::

            This field contains the index of the reference sequence
            in the sequence dictionary. To obtain the name
            of the reference sequence, use :meth:`pysam.Samfile.getrname()`

        """
        def __get__(self): return self._delegate.core.tid
        def __set__(self, tid): self._delegate.core.tid = tid

    property pos:
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

    property bin:
        """properties bin"""
        def __get__(self):
            return pysam_get_bin(self._delegate)
        def __set__(self, bin):
            pysam_set_bin(self._delegate, bin)

    property rlen:
        """length of the read. This includes soft-clipped bases
        and is equal to ``len(seq)``.

        This property is read-only.

        Returns 0 if not available."""
        def __get__(self):
            return self._delegate.core.l_qseq

    property aend:
        '''aligned reference position of the read on the reference genome.  
        
        aend points to one past the last aligned residue.
        Returns None if not available.'''
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            if (self.flag & BAM_FUNMAP) or pysam_get_n_cigar(src) == 0:
                return None
            return bam_endpos(src)

    property alen:
        '''aligned length of the read on the reference genome.

        This is equal to `aend - pos`. Returns None if not available.'''
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            if (self.flag & BAM_FUNMAP) or pysam_get_n_cigar(src) == 0:
                return None
            return bam_endpos(src) - \
                self._delegate.core.pos

    property mapq:
        """mapping quality"""
        def __get__(self):
            return pysam_get_qual(self._delegate)
        def __set__(self, qual):
            pysam_set_qual(self._delegate, qual)

    property mrnm:
        """the :term:`reference` id of the mate
        deprecated, use RNEXT instead.
        """
        def __get__(self):
            return self._delegate.core.mtid
        def __set__(self, mtid):
            self._delegate.core.mtid = mtid

    property rnext:
        """the :term:`reference` id of the mate """
        def __get__(self): return self._delegate.core.mtid
        def __set__(self, mtid): self._delegate.core.mtid = mtid

    property mpos:
        """the position of the mate
        deprecated, use PNEXT instead."""
        def __get__(self):
            return self._delegate.core.mpos
        def __set__(self, mpos):
            self._delegate.core.mpos = mpos

    property pnext:
        """the position of the mate"""
        def __get__(self):
            return self._delegate.core.mpos
        def __set__(self, mpos):
            self._delegate.core.mpos = mpos

    #######################################################################
    #######################################################################
    ## Flags
    #######################################################################
    property isize:
        """the insert size
        deprecated: use tlen instead"""
        def __get__(self):
            return self._delegate.core.isize
        def __set__(self, isize):
            self._delegate.core.isize = isize
    property tlen:
        """the template length"""
        def __get__(self):
            return self._delegate.core.isize
        def __set__(self, isize):
            self._delegate.core.isize = isize
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

    #######################################################################
    #######################################################################
    ## Derived properties
    #######################################################################
    property positions:
        """a list of reference positions that this read aligns to."""
        def __get__(self):
            cdef uint32_t k, i, pos
            cdef int op
            cdef uint32_t * cigar_p
            cdef bam1_t * src

            src = self._delegate
            if pysam_get_n_cigar(src) == 0:
                return []

            result = []
            pos = src.core.pos
            cigar_p = pysam_bam_get_cigar(src)

            for k from 0 <= k < pysam_get_n_cigar(src):
                op = cigar_p[k] & BAM_CIGAR_MASK
                l = cigar_p[k] >> BAM_CIGAR_SHIFT
                if op == BAM_CMATCH:
                    for i from pos <= i < pos + l:
                        result.append( i )

                if op == BAM_CMATCH or op == BAM_CDEL or op == BAM_CREF_SKIP:
                    pos += l

            return result

    property inferred_length:
        """inferred read length from CIGAR string.

        Returns 0 if CIGAR string is not present.
        """
        def __get__(self):
           cdef uint32_t k, qpos
           cdef int op
           cdef uint32_t * cigar_p
           cdef bam1_t * src 

           src = self._delegate
           if pysam_get_n_cigar(src) == 0: return 0

           qpos = 0
           cigar_p = pysam_bam_get_cigar(src)

           for k from 0 <= k < pysam_get_n_cigar(src):
               op = cigar_p[k] & BAM_CIGAR_MASK

               if op == BAM_CMATCH or op == BAM_CINS or op == BAM_CSOFT_CLIP:
                   qpos += cigar_p[k] >> BAM_CIGAR_SHIFT

           return qpos
            
    property aligned_pairs:
        """a list of aligned read and reference positions.

        Unaligned position are marked by None.
        """
        def __get__(self):
            cdef uint32_t k, i, pos, qpos
            cdef int op
            cdef uint32_t * cigar_p
            cdef bam1_t * src 

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

                if op == BAM_CMATCH:
                    for i from pos <= i < pos + l:
                        result.append( (qpos, i) )
                        qpos += 1
                    pos += l

                elif op == BAM_CINS:
                    for i from pos <= i < pos + l:
                        result.append( (qpos, None) )
                        qpos += 1

                elif op == BAM_CDEL or op == BAM_CREF_SKIP:
                    for i from pos <= i < pos + l:
                        result.append( (None, i) )
                    pos += l
                       
            return result

    property blocks:
        """ a list of start and end positions of
        aligned gapless blocks.

        The start and end positions are in genomic 
        coordinates. 
      
        Blocks are not normalized, i.e. two blocks 
        might be directly adjacent. This happens if
        the two blocks are separated by an insertion 
        in the read.
        """

        def __get__(self):
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

    #######################################################################
    #######################################################################
    ## 
    #######################################################################
    def overlap( self, uint32_t start, uint32_t end ):
        """return number of aligned bases of read overlapping the interval *start* and *end*
        on the reference sequence.
        """
        cdef uint32_t k, i, pos, overlap
        cdef int op, o
        cdef uint32_t * cigar_p
        cdef bam1_t * src

        overlap = 0

        src = self._delegate
        if pysam_get_n_cigar(src) == 0:
            return 0
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

    def opt(self, tag):
        """retrieves optional data given a two-letter *tag*"""
        #see bam_aux.c: bam_aux_get() and bam_aux2i() etc
        cdef uint8_t * v
        cdef int nvalues
        btag = _forceBytes(tag)
        v = bam_aux_get(self._delegate, btag)
        if v == NULL: raise KeyError( "tag '%s' not present" % tag )
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
            bytesize, nvalues, values = convertBinaryTagToList( v + 1 )
            return values
        else:
            raise ValueError("unknown auxilliary type '%s'" % auxtype)


    def fancy_str (self):
        """returns list of fieldnames/values in pretty format for debugging
        """
        ret_string = []
        field_names = {
           "tid":           "Contig index",
           "pos":           "Mapped position on contig",
           "mtid":          "Contig index for mate pair",
           "mpos":          "Position of mate pair",
           "isize":         "Insert size",
           "flag":          "Binary flag",
           "n_cigar":       "Count of cigar entries",
           "cigar":         "Cigar entries",
           "qual":          "Mapping quality",
           "bin":           "Bam index bin number",
           "l_qname":       "Length of query name",
           "qname":         "Query name",
           "l_qseq":        "Length of query sequence",
           "qseq":          "Query sequence",
           "bqual":         "Quality scores",
           "l_data":         "Length of auxilary data",
           "m_data":        "Maximum data length",
           }
        fields_names_in_order = ["tid", "pos", "mtid", "mpos", "isize", "flag",
                                 "n_cigar", "cigar", "qual", "bin", "l_qname", "qname",
                                 "l_qseq", "qseq", "bqual", "l_data", "m_data"]

        for f in fields_names_in_order:
            if not f in self.__dict__:
                continue
            ret_string.append("%-30s %-10s= %s" % (field_names[f], "(" + f + ")", self.__getattribute__(f)))

        for f in self.__dict__:
            if not f in field_names:
                ret_string.append("%-30s %-10s= %s" % (f, "", self.__getattribute__(f)))
        return ret_string

cdef class PileupProxy:
    '''A pileup column. A pileup column contains
    all the reads that map to a certain target base.

    tid
        chromosome ID as is defined in the header
    pos
        the target base coordinate (0-based)
    n
        number of reads mapping to this column
    pileups
        list of reads (:class:`pysam.PileupRead`) aligned to this column

    This class is a proxy for results returned by the samtools pileup engine.
    If the underlying engine iterator advances, the results of this column
    will change.
    '''
    def __init__(self):
        raise TypeError("This class cannot be instantiated from Python")

    def __str__(self):
        return "\t".join( map(str, (self.tid, self.pos, self.n))) +\
            "\n" +\
            "\n".join( map(str, self.pileups) )

    property tid:
        '''the chromosome ID as is defined in the header'''
        def __get__(self): return self.tid

    property n:
        '''number of reads mapping to this column.'''
        def __get__(self): return self.n_pu
        def __set__(self, n): self.n_pu = n

    property pos:
        def __get__(self): return self.pos

    property pileups:
        '''list of reads (:class:`pysam.PileupRead`) aligned to this column'''
        def __get__(self):
            cdef int x
            pileups = []

            if self.plp == NULL or self.plp[0] == NULL:
                raise ValueError("PileupProxy accessed after iterator finished")

            # warning: there could be problems if self.n and self.buf are
            # out of sync.
            for x from 0 <= x < self.n_pu:
                pileups.append(makePileupRead(&(self.plp[0][x])))
            return pileups

cdef class PileupRead:
    '''A read aligned to a column.
    '''

    def __init__(self):
        raise TypeError("This class cannot be instantiated from Python")

    def __str__(self):
        return "\t".join( map(str, (self.alignment, self.qpos, self.indel, self.level, self.is_del, self.is_head, self.is_tail ) ) )

    property alignment:
        """a :class:`pysam.AlignedRead` object of the aligned read"""
        def __get__(self):
            return self._alignment
    property qpos:
        """position of the read base at the pileup site, 0-based"""
        def __get__(self):
            return self._qpos
    property indel:
        """indel length; 0 for no indel, positive for ins and negative for del"""
        def __get__(self):
            return self._indel
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
    property level:
        def __get__(self):
            return self._level


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
       '''the snp quality (Phred scaled) - probability of consensus being identical to reference sequence.'''
       def __get__(self): return self._snp_quality

    property mapping_quality:
       '''the root mean square (rms) of the mapping quality of all reads involved in the call.'''
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
    """index a bamfile by read.

    The index is kept in memory.

    By default, the file is re-openend to avoid conflicts if
    multiple operators work on the same file. Set *reopen* = False
    to not re-open *samfile*.
    """

    def __init__(self, Samfile samfile, int reopen=True):

        # makes sure that samfile stays alive as long as this
        # object is alive.
        self.samfile = samfile

        assert samfile.isbam, "can only IndexReads on bam files"

        # reopen the file - note that this makes the iterator
        # slow and causes pileup to slow down significantly.
        if reopen:
            self.htsfile = hts_open(samfile._filename, 'r')
            assert self.htsfile != NULL
            # read header - required for accurate positioning
            sam_hdr_read(self.htsfile)
            self.owns_samfile = True
        else:
            self.htsfile = self.samfile.htsfile
            self.owns_samfile = False

        # TODO: BAM file specific
        self.fp = self.htsfile.fp.bgzf

    def build(self):
        '''build index.'''

        self.index = collections.defaultdict(list)

        # this method will start indexing from the current file position
        # if you decide
        cdef int ret = 1
        cdef bam1_t * b = <bam1_t*>calloc(1, sizeof( bam1_t))

        cdef uint64_t pos

        while ret > 0:
            pos = bgzf_tell(self.fp)
            ret = sam_read1(self.htsfile,
                            self.samfile.header,
                            b)
            if ret > 0:
                qname = _charptr_to_str(pysam_bam_get_qname(b))
                self.index[qname].append(pos)

        bam_destroy1(b)

    def find(self, qname):
        if qname in self.index:
            return IteratorRowSelection(
                self.samfile,
                self.index[qname],
                reopen = False)
        else:
            raise KeyError("read %s not found" % qname)

    def __dealloc__(self):
        if self.owns_samfile:
            hts_close(self.htsfile)

__all__ = ["Samfile",
           "IteratorRow",
           "IteratorColumn",
           "AlignedRead",
           "PileupColumn",
           "PileupProxy",
           "PileupRead",
           "IndexedReads" ]
           # "IteratorSNPCalls",
           # "SNPCaller",
           # "IndelCaller",
           # "IteratorIndelCalls",



