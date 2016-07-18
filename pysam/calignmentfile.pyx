# cython: embedsignature=True
# cython: profile=True
########################################################
########################################################
# Cython wrapper for SAM/BAM/CRAM files based on htslib
########################################################
# The principal classes defined in this module are:
#
# class AlignmentFile   read/write access to SAM/BAM/CRAM formatted files
# 
# class IndexedReads    index a SAM/BAM/CRAM file by query name while keeping
#                       the original sort order intact
# 
# Additionally this module defines numerous additional classes that
# are part of the internal API. These are:
# 
# Various iterator classes to iterate over alignments in sequential
# (IteratorRow) or in a stacked fashion (IteratorColumn):
# 
# class IteratorRow
# class IteratorRowRegion
# class IteratorRowHead
# class IteratorRowAll
# class IteratorRowAllRefs
# class IteratorRowSelection
# class IteratorColumn
# class IteratorColumnRegion
# class IteratorColumnAllRefs
#
########################################################
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
########################################################
import os
import collections
import re
import warnings
import array

from cpython cimport array as c_array
from cpython.version cimport PY_MAJOR_VERSION

from pysam.cutils cimport force_bytes, force_str, charptr_to_str
from pysam.cutils cimport encode_filename, from_string_and_size
from pysam.calignedsegment cimport makeAlignedSegment, makePileupColumn
from pysam.chtslib cimport hisremote

if PY_MAJOR_VERSION >= 3:
    from io import StringIO
else:
    from StringIO import StringIO

cimport cython

########################################################
## Constants and global variables

# defines imported from samtools
DEF SEEK_SET = 0
DEF SEEK_CUR = 1
DEF SEEK_END = 2

# maximum genomic coordinace
cdef int MAX_POS = 2 << 29

# valid types for SAM headers
VALID_HEADER_TYPES = {"HD" : dict,
                      "SQ" : list,
                      "RG" : list,
                      "PG" : list,
                      "CO" : list}

# order of records within SAM headers
VALID_HEADERS = ("HD", "SQ", "RG", "PG", "CO")

# default type conversions within SAM header records
KNOWN_HEADER_FIELDS = {"HD" : {"VN" : str, "SO" : str, "GO" : str},
                       "SQ" : {"SN" : str, "LN" : int, "AS" : str, 
                               "M5" : str, "SP" : str, "UR" : str,},
                       "RG" : {"ID" : str, "CN" : str, "DS" : str,
                               "DT" : str, "FO" : str, "KS" : str,
                               "LB" : str, "PG" : str, "PI" : str,
                               "PL" : str, "PM" : str, "PU" : str,
                               "SM" : str,},
                       "PG" : {"ID" : str, "PN" : str, "CL" : str, 
                               "PP" : str, "DS" : str, "VN" : str,},}

# output order of fields within records. Ensure that CL is at
# the end as parsing a CL will ignore any subsequent records.
VALID_HEADER_ORDER = {"HD" : ("VN", "SO", "GO"),
                      "SQ" : ("SN", "LN", "AS", "M5",
                               "UR", "SP"),
                      "RG" : ("ID", "CN", "SM", "LB",
                              "PU", "PI", "DT", "DS",
                              "PL", "FO", "KS", "PG",
                              "PM"),
                      "PG" : ("PN", "ID", "VN", "PP",
                              "DS", "CL"),}


def build_header_line(fields, record):
    '''build a header line from `fields` dictionary for `record`'''

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

cdef bam_hdr_t * build_header(new_header):
    '''return a new header built from a dictionary in `new_header`.

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
                lines.append(build_header_line(data, record))
            else:
                for fields in new_header[record]:
                    lines.append(build_header_line(fields, record))

    # then: user tags (lower case), sorted alphabetically
    for record, data in sorted(new_header.items()):
        if record in VALID_HEADERS: continue
        if type(data) is dict:
            lines.append(build_header_line(data, record))
        else:
            for fields in new_header[record]:
                lines.append(build_header_line(fields, record))

    text = "\n".join(lines) + "\n"
    if dest.text != NULL: free( dest.text )
    dest.text = <char*>calloc(len(text), sizeof(char))
    dest.l_text = len(text)
    cdef bytes btext = text.encode('ascii')
    strncpy(dest.text, btext, dest.l_text)

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


cdef class AlignmentFile:
    """AlignmentFile(filepath_or_object, mode=None, template=None,
    reference_names=None, reference_lengths=None, text=NULL,
    header=None, add_sq_text=False, check_header=True, check_sq=True,
    reference_filename=None, filename=None)

    A :term:`SAM`/:term:`BAM` formatted file. 

    If `filepath_or_object` is a string, the file is automatically
    opened. If `filepath_or_object` is a python File object, the
    already opened file will be used.

    If the file is opened for reading an index for a BAM file exists
    (.bai), it will be opened automatically. Without an index random
    access via :meth:`~pysam.AlignmentFile.fetch` and
    :meth:`~pysam.AlignmentFile.pileup` is disabled.

    For writing, the header of a :term:`SAM` file/:term:`BAM` file can
    be constituted from several sources (see also the samtools format
    specification):

        1. If `template` is given, the header is copied from a another
           `AlignmentFile` (`template` must be a
           :class:`~pysam.AlignmentFile`).

        2. If `header` is given, the header is built from a
           multi-level dictionary. 

        3. If `text` is given, new header text is copied from raw
           text.

        4. The names (`reference_names`) and lengths
           (`reference_lengths`) are supplied directly as lists.

    When reading or writing a CRAM file, the filename of a FASTA-formatted
    reference can be specified with `reference_filename`.

    By default, if a file is opened in mode 'r', it is checked
    for a valid header (`check_header` = True) and a definition of
    chromosome names (`check_sq` = True).

    Parameters
    ----------
    mode : string
        `mode` should be ``r`` for reading or ``w`` for writing. The
        default is text mode (:term:`SAM`). For binary (:term:`BAM`)
        I/O you should append ``b`` for compressed or ``u`` for
        uncompressed :term:`BAM` output.  Use ``h`` to output header
        information in text (:term:`TAM`) mode. Use ``c`` for
        :term:`CRAM` formatted files.

        If ``b`` is present, it must immediately follow ``r`` or
        ``w``.  Valid modes are ``r``, ``w``, ``wh``, ``rb``, ``wb``,
        ``wbu``, ``wb0``, ``rc`` and ``wc``. For instance, to open a
        :term:`BAM` formatted file for reading, type::

           f = pysam.AlignmentFile('ex1.bam','rb')

        If mode is not specified, the method will try to auto-detect
        in the order 'rb', 'r', thus both the following should work::

            f1 = pysam.AlignmentFile('ex1.bam')
            f2 = pysam.AlignmentFile('ex1.sam')

    template : AlignmentFile
        when writing, copy  header frem `template`.

    header :  dict
        when writing, build header from a multi-level dictionary. The
        first level are the four types ('HD', 'SQ', ...). The second
        level are a list of lines, with each line being a list of
        tag-value pairs. The header is constructed first from all the
        defined fields, followed by user tags in alphabetical order.

    text : string
        when writing, use the string provided as the header

    reference_names : list
        see referece_lengths

    reference_lengths : list
        when writing, build header from list of chromosome names and
        lengths.  By default, 'SQ' and 'LN' tags will be added to the
        header text. This option can be changed by unsetting the flag
        `add_sq_text`.

    add_sq_text : bool
        do not add 'SQ' and 'LN' tags to header. This option permits
        construction :term:`SAM` formatted files without a header.

    check_header : bool
        when reading, check if header is present (default=True)

    check_sq : bool
        when reading, check if SQ entries are present in header
        (default=True)

    reference_filename : string
        Path to a FASTA-formatted reference file. Valid only for CRAM files.
        When reading a CRAM file, this overrides both ``$REF_PATH`` and the URL
        specified in the header (``UR`` tag), which are normally used to find
        the reference.

    filename : string
        Alternative to filepath_or_object. Filename of the file
        to be opened.

    """

    def __cinit__(self, *args, **kwargs):

        self.htsfile = NULL
        self._filename = None
        self.is_bam = False
        self.is_stream = False
        self.is_cram = False
        self.is_remote = False

        if "filename" in kwargs:
            args = [kwargs["filename"]]
            del kwargs["filename"]

        self._open(*args, **kwargs)

        # allocate memory for iterator
        self.b = <bam1_t*>calloc(1, sizeof(bam1_t))

    def is_open(self):
        '''return true if htsfile has been opened.'''
        return self.htsfile != NULL

    def has_index(self):
        """return true if htsfile has an existing (and opened) index.
        """
        return self.index != NULL

    def check_index(self):
        """return True if index is present.

        Raises
        ------

        AttributeError
            if htsfile is :term:`SAM` formatted and thus has no index.

        ValueError
            if htsfile is closed or index could not be opened.
        """

        if not self.is_open():
            raise ValueError("I/O operation on closed file")
        if not self.is_bam and not self.is_cram:
            raise AttributeError(
                "AlignmentFile.mapped only available in bam files")
        if self.index == NULL:
            raise ValueError(
                "mapping information not recorded in index "
                "or index not available")
        return True

    def _open(self,
              filepath_or_object,
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
              filepath_index=None,
              referencenames=None,
              referencelengths=None):
        '''open a sam, bam or cram formatted file.

        If _open is called on an existing file, the current file
        will be closed and a new file will be opened.
        '''
        cdef char *cfilename
        cdef char *creference_filename
        cdef char *cindexname
        cdef char *cmode

        # for backwards compatibility:
        if referencenames is not None:
            reference_names = referencenames
        if referencelengths is not None:
            reference_lengths = referencelengths

        # autodetection for read
        if mode is None:
            mode = "r"

        assert mode in ("r", "w", "rb", "wb", "wh",
                        "wbu", "rU", "wb0",
                        "rc", "wc"), \
            "invalid file opening mode `%s`" % mode

        # close a previously opened file
        if self.htsfile != NULL:
            self.close()

        # StringIO not supported
        if isinstance(filepath_or_object, StringIO):
            filename = "stringio"
            raise NotImplementedError(
                "access from StringIO objects not supported")
            if filepath_or_object.closed:
                raise ValueError('I/O operation on closed StringIO object')
        # check if we are working with a File object
        elif hasattr(filepath_or_object, "fileno"):
            filename = filepath_or_object.name
            if filepath_or_object.closed:
                raise ValueError('I/O operation on closed file')
        else:
            filename = filepath_or_object

        # for htslib, wbu seems to not work
        if mode == "wbu":
            mode = "wb0"

        cdef bytes bmode = mode.encode('ascii')
        self._filename = filename = encode_filename(filename)
        self._reference_filename = reference_filename = encode_filename(
            reference_filename)

        # FIXME: Use htsFormat when it is available
        self.is_stream = filename == b"-"
        self.is_remote = hisremote(filename)

        cdef char * ctext
        cdef hFILE * fp
        ctext = NULL

        if mode[0] == 'w':
            # open file for writing

            # header structure (used for writing)
            if template:
                self.header = bam_hdr_dup(template.header)
            elif header:
                self.header = build_header(header)
            else:
                # build header from a target names and lengths
                assert reference_names and reference_lengths, \
                    ("either supply options `template`, `header` "
                     "or  both `reference_names` and `reference_lengths` "
                     "for writing")
                assert len(reference_names) == len(reference_lengths), \
                    "unequal names and lengths of reference sequences"

                # allocate and fill header
                reference_names = [force_bytes(ref) for ref in reference_names]
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
                                    (force_str(reference_names[x]), 
                                     reference_lengths[x]))
                    text = ''.join(text)

                if text is not None:
                    # copy without \0
                    text = force_bytes(text)
                    ctext = text
                    self.header.l_text = strlen(ctext)
                    self.header.text = <char*>calloc(
                        strlen(ctext), sizeof(char))
                    memcpy(self.header.text, ctext, strlen(ctext))

            # open file (hts_open is synonym with sam_open)
            cfilename, cmode = filename, bmode
            if hasattr(filepath_or_object, "fileno"):
                fp = hdopen(filepath_or_object.fileno(), cmode)
                with nogil:
                    self.htsfile = hts_hopen(fp, cfilename, cmode)
            else:
                with nogil:
                    self.htsfile = hts_open(cfilename, cmode)

            # htsfile.format does not get set until writing, so use
            # the format specifier explicitely given by the user.
            self.is_bam = "b" in mode
            self.is_cram = "c" in mode

            # set filename with reference sequences. If no filename
            # is given, the CRAM reference arrays will be built from
            # the @SQ header in the header
            if self.is_cram and reference_filename:
                # note that fn_aux takes ownership, so create a copy
                self.htsfile.fn_aux = strdup(self._reference_filename)

            # write header to htsfile
            if self.is_bam or self.is_cram or "h" in mode:
                with nogil:
                    sam_hdr_write(self.htsfile, self.header)

        elif mode[0] == "r":
            # open file for reading
            if (filename != b"-"
                and not self.is_remote
                and not os.path.exists(filename)):
                raise IOError("file `%s` not found" % filename)

            # open file (hts_open is synonym with sam_open)
            cfilename, cmode = filename, bmode
            if hasattr(filepath_or_object, "fileno"):
                fp = hdopen(filepath_or_object.fileno(), cmode)
                with nogil:
                    self.htsfile = hts_hopen(fp, cfilename, cmode)
            else:
                with nogil:
                    self.htsfile = hts_open(cfilename, cmode)

            if self.htsfile == NULL:
                raise ValueError(
                    "could not open file (mode='%s') - "
                    "is it SAM/BAM format?" % mode)

            self.is_bam = self.htsfile.format.format == bam
            self.is_cram = self.htsfile.format.format == cram

            # bam files require a valid header
            if self.is_bam or self.is_cram:
                with nogil:
                    self.header = sam_hdr_read(self.htsfile)
                if self.header == NULL:
                    raise ValueError(
                        "file does not have valid header (mode='%s') "
                        "- is it BAM format?" % mode )
            else:
                # in sam files it is optional (htsfile full of
                # unmapped reads)
                if check_header:
                    with nogil:
                        self.header = sam_hdr_read(self.htsfile)
                    if self.header == NULL:
                        raise ValueError(
                            "file does not have valid header (mode='%s') "
                            "- is it SAM format?" % mode )
                    # self.header.ignore_sam_err = True

            # set filename with reference sequences
            if self.is_cram and reference_filename:
                creference_filename = self._reference_filename
                hts_set_opt(self.htsfile,
                            CRAM_OPT_REFERENCE,
                            creference_filename)

            if check_sq and self.header.n_targets == 0:
                raise ValueError(
                    ("file has no sequences defined (mode='%s') - "
                     "is it SAM/BAM format? Consider opening with "
                     "check_sq=False") % mode)

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
            if self.is_remote and not filepath_index:
                cfilename = filename

                with nogil:
                    self.index = hts_idx_load(cfilename, format_index)
                if self.index == NULL:
                    warnings.warn(
                        "unable to open remote index for '%s'" % cfilename)
            else:
                has_index = True
                cfilename = filename
                if filepath_index:
                    if not os.path.exists(filepath_index):
                        warnings.warn(
                            "unable to open index at %s" % cfilename)
                        self.index = NULL
                        has_index = False
                else:
                    if self.is_bam \
                            and not os.path.exists(filename + b".bai") \
                            and not os.path.exists(filename[:-4] + b".bai"):
                        self.index = NULL
                        has_index = False
                    elif self.is_cram \
                            and not os.path.exists(filename + b".crai") \
                            and not os.path.exists(filename[:-5] + b".crai"):
                        self.index = NULL
                        has_index = False

                if has_index:
                    # returns NULL if there is no index or index could
                    # not be opened
                    if filepath_index:
                        cindexname = filepath_index = encode_filename(filepath_index)
                        with nogil:
                            self.index = sam_index_load2(self.htsfile,
                                                         cfilename,
                                                         cindexname)

                    else:
                        with nogil:
                            self.index = sam_index_load(self.htsfile,
                                                        cfilename)
                    if self.index == NULL:
                        raise IOError(
                            "error while opening index for '%s'" %
                            filename)

            # save start of data section
            if not self.is_stream:
                self.start_offset = self.tell()

    def get_tid(self, reference):
        """
        return the numerical :term:`tid` corresponding to
        :term:`reference`

        returns -1 if reference is not known.
        """
        if not self.is_open():
            raise ValueError("I/O operation on closed file")
        reference = force_bytes(reference)
        return bam_name2id(self.header, reference)

    def get_reference_name(self, tid):
        """
        return :term:`reference` name corresponding to numerical :term:`tid`
        """
        if not self.is_open():
            raise ValueError("I/O operation on closed file")
        if not 0 <= tid < self.header.n_targets:
            raise ValueError("reference_id %i out of range 0<=tid<%i" % 
                             (tid, self.header.n_targets))
        return charptr_to_str(self.header.target_name[tid])

    def reset(self):
        """reset file position to beginning of file just after
        the header.

        Returns
        -------

        The file position after moving the file pointer.

        """
        return self.seek(self.start_offset, 0)

    def seek(self, uint64_t offset, int where=0):
        """move file pointer to position `offset`, see
        :meth:`pysam.AlignmentFile.tell`.

        Parameters
        ----------
        
        offset : int

        position of the read/write pointer within the file.

        where : int
    
        optional and defaults to 0 which means absolute file
        positioning, other values are 1 which means seek relative to
        the current position and 2 means seek relative to the file's
        end.
        
        Returns
        -------
        
        the file position after moving the file pointer

        """

        if not self.is_open():
            raise ValueError("I/O operation on closed file")
        if not self.is_bam:
            raise NotImplementedError(
                "seek only available in bam files")
        if self.is_stream:
            raise OSError("seek no available in streams")

        cdef uint64_t pos
        with nogil:
            pos = bgzf_seek(hts_get_bgzfp(self.htsfile), offset, where)
        return pos

    def tell(self):
        """
        return current file position.
        """
        if not self.is_open():
            raise ValueError("I/O operation on closed file")
        if not (self.is_bam or self.is_cram):
            raise NotImplementedError(
                "seek only available in bam files")

        cdef uint64_t pos
        with nogil:
            pos = bgzf_tell(hts_get_bgzfp(self.htsfile))
        return pos

    def parse_region(self,
                     reference=None,
                     start=None,
                     end=None,
                     region=None,
                     tid=None):
        """parse alternative ways to specify a genomic region. A region can
        either be specified by :term:`reference`, `start` and
        `end`. `start` and `end` denote 0-based, half-open
        intervals.

        Alternatively, a samtools :term:`region` string can be
        supplied.
        
        If any of the coordinates are missing they will be replaced by the
        minimum (`start`) or maximum (`end`) coordinate.

        Note that region strings are 1-based, while `start` and `end` denote
        an interval in python coordinates.

        Returns
        -------
        
        tuple :  a tuple of `flag`, :term:`tid`, `start` and `end`. The
        flag indicates whether no coordinates were supplied and the
        genomic region is the complete genomic space.

        Raises
        ------
        
        ValueError
           for invalid or out of bounds regions.

        """
        cdef int rtid
        cdef long long rstart
        cdef long long rend

        rtid = -1
        rstart = 0
        rend = MAX_POS
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
            region = force_str(region)
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
        if not 0 <= rstart < MAX_POS:
            raise ValueError('start out of range (%i)' % rstart)
        if not 0 <= rend <= MAX_POS:
            raise ValueError('end out of range (%i)' % rend)

        return 1, rtid, rstart, rend

    def fetch(self,
              reference=None,
              start=None,
              end=None,
              region=None,
              tid=None,
              until_eof=False,
              multiple_iterators=False):
        """fetch reads aligned in a :term:`region`. 

        See :meth:`AlignmentFile.parse_region` for more information
        on genomic regions.

        Without a `reference` or `region` all mapped reads in the file
        will be fetched. The reads will be returned ordered by reference
        sequence, which will not necessarily be the order within the
        file. This mode of iteration still requires an index. If there is
        no index, use `until_eof=True`.

        If only `reference` is set, all reads aligned to `reference`
        will be fetched.

        A :term:`SAM` file does not allow random access. If `region`
        or `reference` are given, an exception is raised.

        :class:`~pysam.FastaFile`
        :class:`~pysam.IteratorRow`
        :class:`~pysam.IteratorRow`
        :class:`~IteratorRow`
        :class:`IteratorRow`

        Parameters
        ----------
        
        until_eof : bool

           If `until_eof` is True, all reads from the current file
           position will be returned in order as they are within the
           file. Using this option will also fetch unmapped reads.

        multiple_iterators : bool
           
           If `multiple_iterators` is True, multiple
           iterators on the same file can be used at the same time. The
           iterator returned will receive its own copy of a filehandle to
           the file effectively re-opening the file. Re-opening a file
           creates some overhead, so beware.

        Returns
        -------

        An iterator over a collection of reads.

        Raises
        ------

        ValueError
            if the genomic coordinates are out of range or invalid or the
            file does not permit random access to genomic coordinates.

        """
        cdef int rtid, rstart, rend, has_coord

        if not self.is_open():
            raise ValueError( "I/O operation on closed file" )

        has_coord, rtid, rstart, rend = self.parse_region(
            reference,
            start,
            end,
            region,
            tid)

        # Turn of re-opening if htsfile is a stream
        if self.is_stream:
            multiple_iterators = False

        if self.is_bam or self.is_cram:
            if not until_eof and not self.is_remote:
                if not self.has_index():
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
        '''return an iterator over the first n alignments. 

        This iterator is is useful for inspecting the bam-file.

        Parameters
        ----------

        multiple_iterators : bool
        
            is set to True by default in order to
            avoid changing the current file position.
        
        Returns
        -------
        
        an iterator over a collection of reads
        
        '''
        return IteratorRowHead(self, n,
                               multiple_iterators=multiple_iterators)

    def mate(self, AlignedSegment read):
        '''return the mate of :class:`~pysam.AlignedSegment` `read`.

        .. note::

            Calling this method will change the file position.
            This might interfere with any iterators that have
            not re-opened the file.

        .. note::
  
           This method is too slow for high-throughput processing.
           If a read needs to be processed with its mate, work
           from a read name sorted file or, better, cache reads.

        Returns
        -------
        
        :class:`~pysam.AlignedSegment` : the mate

        Raises
        ------

        ValueError
            if the read is unpaired or the mate is unmapped

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

    def pileup(self,
               reference=None,
               start=None,
               end=None,
               region=None,
               **kwargs):
        """perform a :term:`pileup` within a :term:`region`. The region is
        specified by :term:`reference`, 'start' and 'end' (using
        0-based indexing).  Alternatively, a samtools 'region' string
        can be supplied.

        Without 'reference' or 'region' all reads will be used for the
        pileup. The reads will be returned ordered by
        :term:`reference` sequence, which will not necessarily be the
        order within the file.

        Note that :term:`SAM` formatted files do not allow random
        access.  In these files, if a 'region' or 'reference' are
        given an exception is raised.

        .. note::

            'all' reads which overlap the region are returned. The
            first base returned will be the first base of the first
            read 'not' necessarily the first base of the region used
            in the query.

        Parameters
        ----------

        stepper : string
           The stepper controlls how the iterator advances.
           Possible options for the stepper are

           ``all``
              skip reads in which any of the following flags are set:
              BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP

           ``nofilter``
              uses every single read

           ``samtools``
              same filter and read processing as in :term:`csamtools`
              pileup. This requires a 'fastafile' to be given.


        fastafile : :class:`~pysam.FastaFile` object.

           This is required for some of the steppers.

        max_depth : int
           Maximum read depth permitted. The default limit is '8000'.

        truncate : bool

           By default, the samtools pileup engine outputs all reads
           overlapping a region. If truncate is True and a region is
           given, only columns in the exact region specificied are
           returned.

        Returns
        -------

        an iterator over genomic positions.

        """
        cdef int rtid, rstart, rend, has_coord

        if not self.is_open():
            raise ValueError("I/O operation on closed file")

        has_coord, rtid, rstart, rend = self.parse_region(
            reference, start, end, region)

        if self.is_bam or self.is_cram:
            if not self.has_index():
                raise ValueError("no index available for pileup")

            if has_coord:
                return IteratorColumnRegion(self,
                                            tid=rtid,
                                            start=rstart,
                                            end=rend,
                                            **kwargs )
            else:
                return IteratorColumnAllRefs(self, **kwargs )

        else:
            raise NotImplementedError(
                "pileup of samfiles not implemented yet")

    def count(self,
              reference=None,
              start=None,
              end=None,
              region=None,
              until_eof=False,
              read_callback="nofilter"):
        '''count the number of reads in :term:`region`

        The region is specified by :term:`reference`, `start` and
        `end`. Alternatively, a :term:`samtools` :term:`region` string
        can be supplied.

        A :term:`SAM` file does not allow random access and if
        `region` or `reference` are given, an exception is raised.

        Parameters
        ----------
        
        reference : string
            reference_name of the genomic region (chromosome)

        start : int
            start of the genomic region

        end : int
            end of the genomic region
        
        region : string
            a region string in samtools format.

        until_eof : bool
            count until the end of the file, possibly including 
            unmapped reads as well.

        read_callback: string or function

            select a call-back to ignore reads when counting. It can
            be either a string with the following values:

            ``all``
                skip reads in which any of the following
                flags are set: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL,
                BAM_FDUP

            ``nofilter``
                uses every single read

            Alternatively, `read_callback` can be a function
            ``check_read(read)`` that should return True only for
            those reads that shall be included in the counting.

        Raises
        ------

        ValueError
            if the genomic coordinates are out of range or invalid.

        '''
        cdef AlignedSegment read
        cdef long counter = 0

        if not self.is_open():
            raise ValueError( "I/O operation on closed file" )

        cdef int filter_method = 0
        if read_callback == "all":
            filter_method = 1
        elif read_callback == "nofilter":
            filter_method = 2

        for read in self.fetch(reference=reference,
                               start=start,
                               end=end,
                               region=region,
                               until_eof=until_eof):
            # apply filter
            if filter_method == 1:
                # filter = "all"
                if (read.flag & (0x4 | 0x100 | 0x200 | 0x400)):
                    continue
            elif filter_method == 2:
                # filter = "nofilter"
                pass
            else:
                if not read_callback(read):
                    continue
            counter += 1

        return counter

    @cython.boundscheck(False)  # we do manual bounds checking
    def count_coverage(self, 
                       reference=None,
                       start=None,
                       end=None,
                       region=None,
                       quality_threshold=15,
                       read_callback='all'):
        """count the coverage of genomic positions by reads in :term:`region`.

        The region is specified by :term:`reference`, `start` and
        `end`. Alternatively, a :term:`samtools` :term:`region` string
        can be supplied. The coverage is computed per-base [ACGT].

        Parameters
        ----------
        
        reference : string
            reference_name of the genomic region (chromosome)

        start : int
            start of the genomic region

        end : int
            end of the genomic region

        region : int
            a region string.

        quality_threshold : int
            quality_threshold is the minimum quality score (in phred) a
            base has to reach to be counted. 

        read_callback: string or function

            select a call-back to ignore reads when counting. It can
            be either a string with the following values:

            ``all``
                skip reads in which any of the following
                flags are set: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL,
                BAM_FDUP

            ``nofilter``
                uses every single read

            Alternatively, `read_callback` can be a function
            ``check_read(read)`` that should return True only for
            those reads that shall be included in the counting.

        Raises
        ------

        ValueError
            if the genomic coordinates are out of range or invalid.

        Returns
        -------

        four array.arrays of the same length in order A C G T : tuple

        """
        
        cdef int _start = start
        cdef int _stop = end
        cdef int length = _stop - _start
        cdef c_array.array int_array_template = array.array('L', [])
        cdef c_array.array count_a
        cdef c_array.array count_c
        cdef c_array.array count_g
        cdef c_array.array count_t
        count_a = c_array.clone(int_array_template, length, zero=True)
        count_c = c_array.clone(int_array_template, length, zero=True)
        count_g = c_array.clone(int_array_template, length, zero=True)
        count_t = c_array.clone(int_array_template, length, zero=True)

        cdef AlignedSegment read
        cdef cython.str seq
        cdef c_array.array quality
        cdef int qpos
        cdef int refpos
        cdef int c = 0
        cdef int filter_method = 0
        if read_callback == "all":
            filter_method = 1
        elif read_callback == "nofilter":
            filter_method = 2
    
        cdef int _threshold = quality_threshold
        for read in self.fetch(reference=reference,
                               start=start,
                               end=end,
                               region=region):
            # apply filter
            if filter_method == 1:
                # filter = "all"
                if (read.flag & (0x4 | 0x100 | 0x200 | 0x400)):
                    continue
            elif filter_method == 2:
                # filter = "nofilter"
                pass
            else:
                if not read_callback(read):
                    continue

            # count
            seq = read.seq
            quality = read.query_qualities
            for qpos, refpos in read.get_aligned_pairs(True):
                if qpos is not None and refpos is not None and \
                   _start <= refpos < _stop:
                    if quality[qpos] >= quality_threshold:
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
        # AH: I have removed the call to close. Even though it is working,
        # it seems to be dangerous according to the documentation as the
        # object be partially deconstructed already.
        if self.htsfile != NULL:
            hts_close(self.htsfile)
            hts_idx_destroy(self.index);
            self.htsfile = NULL

        bam_destroy1(self.b)
        if self.header != NULL:
            bam_hdr_destroy(self.header)
            
    cpdef int write(self, AlignedSegment read) except -1:
        '''
        write a single :class:`pysam.AlignedSegment` to disk.

        Raises
        ------
        ValueError
            if the writing failed

        Returns
        -------
            
        int : the number of bytes written. If the file is closed,
              this will be 0.
        '''
        if not self.is_open():
            return 0

        cdef int ret

        with nogil:
            ret = sam_write1(self.htsfile,
                             self.header,
                             read._delegate)

        # kbj: Still need to raise an exception with except -1. Otherwise
        #      when ret == -1 we get a "SystemError: error return without
        #      exception set".
        if ret < 0:
            raise ValueError('sam write failed')

        return ret

    # context manager interface
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
    property closed:
        """bool indicating the current state of the file object. 
        This is a read-only attribute; the close() method changes the value. 
        """
        def __get__(self):
            return not self.is_open()

    property filename:
        """filename associated with this object. This is a read-only attribute."""
        def __get__(self):
            return self._filename

    property nreferences:
        """"int with the number of :term:`reference` sequences in the file.
        This is a read-only attribute."""
        def __get__(self):
            if not self.is_open():
                raise ValueError("I/O operation on closed file")
            return self.header.n_targets

    property references:
        """tuple with the names of :term:`reference` sequences. This is a 
        read-only attribute"""
        def __get__(self):
            if not self.is_open(): raise ValueError( "I/O operation on closed file" )
            t = []
            for x from 0 <= x < self.header.n_targets:
                t.append(charptr_to_str(self.header.target_name[x]))
            return tuple(t)

    property lengths:
        """tuple of the lengths of the :term:`reference` sequences. This is a
        read-only attribute. The lengths are in the same order as
        :attr:`pysam.AlignmentFile.references`

        """
        def __get__(self):
            if not self.is_open():
                raise ValueError("I/O operation on closed file")
            t = []
            for x from 0 <= x < self.header.n_targets:
                t.append(self.header.target_len[x])
            return tuple(t)

    property mapped:
        """int with total number of mapped alignments according to the
        statistics recorded in the index. This is a read-only
        attribute.
        """
        def __get__(self):
            self.check_index()
            cdef int tid
            cdef uint64_t total = 0
            cdef uint64_t mapped, unmapped
            for tid from 0 <= tid < self.header.n_targets:
                with nogil:
                    hts_idx_get_stat(self.index, tid, &mapped, &unmapped)
                total += mapped
            return total

    property unmapped:
        """int with total number of unmapped reads according to the statistics
        recorded in the index. This number of reads includes the number of reads
        without coordinates. This is a read-only attribute.
        """
        def __get__(self):
            self.check_index()
            cdef int tid
            cdef uint64_t total = hts_idx_get_n_no_coor(self.index)
            cdef uint64_t mapped, unmapped
            for tid from 0 <= tid < self.header.n_targets:
                with nogil:
                    hts_idx_get_stat(self.index, tid, &mapped, &unmapped)
                total += unmapped
            return total

    property nocoordinate:
        """int with total number of reads without coordinates according to the
        statistics recorded in the index. This is a read-only attribute.
        """
        def __get__(self):
            self.check_index()
            cdef uint64_t n
            with nogil:
                n = hts_idx_get_n_no_coor(self.index)
            return n

    property format:
        '''string describing the file format'''
        def __get__(self):
            if not self.is_open():
                raise ValueError( "I/O operation on closed file" )
            return hts_format_description(&self.htsfile.format)

    property text:
        '''string with the full contents of the :term:`sam file` header as a
        string. 

        This is a read-only attribute.
        
        See :attr:`pysam.AlignmentFile.header` to get a parsed
        representation of the header.
        '''
        def __get__(self):
            if not self.is_open():
                raise ValueError( "I/O operation on closed file" )
            return from_string_and_size(self.header.text, self.header.l_text)

    property header:
        """two-level dictionay with header information from the file. 
        
        This is a read-only attribute.

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

        """
        def __get__(self):
            if not self.is_open():
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
                            x[key] = KNOWN_HEADER_FIELDS[record][key](value)
                            break

                        # interpret type of known header record tags, default to str
                        x[key] = KNOWN_HEADER_FIELDS[record].get(key, str)(value)

                    if VALID_HEADER_TYPES[record] == dict:
                        if record in result:
                            raise ValueError(
                                "multiple '%s' lines are not permitted" % record)

                        result[record] = x
                    elif VALID_HEADER_TYPES[record] == list:
                        if record not in result: result[record] = []
                        result[record].append(x)

                # if there are no SQ lines in the header, add the
                # reference names from the information in the bam
                # file.
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

    ###############################################################
    ## file-object like iterator access
    ## note: concurrent access will cause errors (see IteratorRow
    ## and multiple_iterators)
    ## Possible solutions: deprecate or open new file handle
    def __iter__(self):
        if not self.is_open():
            raise ValueError("I/O operation on closed file")

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
        cdef int ret
        with nogil:
            ret = sam_read1(self.htsfile,
                            self.header,
                            self.b)
        return ret

    def __next__(self):
        cdef int ret = self.cnext()
        if (ret >= 0):
            return makeAlignedSegment(self.b, self)
        elif ret == -2:
            raise IOError('truncated file')
        else:
            raise StopIteration
            
    # Compatibility functions for pysam < 0.8.3
    def gettid(self, reference):
        """deprecated, use get_tid() instead"""
        return self.get_tid(reference)
        
    def getrname(self, tid):
        """deprecated, use get_reference_name() instead"""
        return self.get_reference_name(tid)


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
        explicitly. It is returned as a result of call to a
        :meth:`AlignmentFile.fetch`.

    '''

    def __init__(self, AlignmentFile samfile, int multiple_iterators=False):
        cdef char *cfilename
        cdef char *creference_filename
        
        if not samfile.is_open():
            raise ValueError("I/O operation on closed file")

        # makes sure that samfile stays alive as long as the
        # iterator is alive
        self.samfile = samfile

        # reopen the file - note that this makes the iterator
        # slow and causes pileup to slow down significantly.
        if multiple_iterators:
            cfilename = samfile._filename
            with nogil:
                self.htsfile = hts_open(cfilename, 'r')
            assert self.htsfile != NULL
            # read header - required for accurate positioning
            # could a tell/seek work?
            with nogil:
                self.header = sam_hdr_read(self.htsfile)
            assert self.header != NULL
            self.owns_samfile = True
            # options specific to CRAM files
            if samfile.is_cram and samfile._reference_filename:
                creference_filename = samfile._reference_filename
                hts_set_opt(self.htsfile,
                            CRAM_OPT_REFERENCE,
                            creference_filename)

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
        explicitly. It is returned as a result of call to a
        :meth:`AlignmentFile.fetch`.

    """

    def __init__(self, AlignmentFile samfile,
                 int tid, int beg, int end,
                 int multiple_iterators=False):

        IteratorRow.__init__(self, samfile,
                             multiple_iterators=multiple_iterators)

        if not samfile.has_index():
            raise ValueError("no index available for iteration")

        with nogil:
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
        with nogil:
            self.retval = hts_itr_next(hts_get_bgzfp(self.htsfile),
                                       self.iter,
                                       self.b,
                                       self.htsfile)

    def __next__(self):
        self.cnext()
        if self.retval >= 0:
            return makeAlignedSegment(self.b, self.samfile)
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

    iterate over first n reads in `samfile`

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
        cdef int ret
        with nogil:
            ret = sam_read1(self.htsfile,
                            self.samfile.header,
                            self.b)
        return ret

    def __next__(self):
        if self.current_row >= self.max_rows:
            raise StopIteration

        cdef int ret = self.cnext()
        if ret >= 0:
            self.current_row += 1
            return makeAlignedSegment(self.b, self.samfile)
        elif ret == -2:
            raise IOError('truncated file')
        else:
            raise StopIteration


cdef class IteratorRowAll(IteratorRow):
    """*(AlignmentFile samfile, int multiple_iterators=False)*

    iterate over all reads in `samfile`

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
        cdef int ret
        with nogil:
            ret = sam_read1(self.htsfile,
                            self.samfile.header,
                            self.b)
        return ret

    def __next__(self):
        cdef int ret = self.cnext()
        if ret >= 0:
            return makeAlignedSegment(self.b, self.samfile)
        elif ret == -2:
            raise IOError('truncated file')
        else:
            raise StopIteration


cdef class IteratorRowAllRefs(IteratorRow):
    """iterates over all mapped reads by chaining iterators over each
    reference

    .. note::
        It is usually not necessary to create an object of this class
        explicitly. It is returned as a result of call to a
        :meth:`AlignmentFile.fetch`.

    """

    def __init__(self, AlignmentFile samfile,
                 multiple_iterators=False):

        IteratorRow.__init__(self, samfile,
                             multiple_iterators=multiple_iterators)

        if not samfile.has_index():
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
                return makeAlignedSegment(self.rowiter.b, self.samfile)

            self.tid += 1

            # Otherwise, proceed to next reference or stop
            if self.tid < self.samfile.nreferences:
                self.nextiter()
            else:
                raise StopIteration


cdef class IteratorRowSelection(IteratorRow):
    """*(AlignmentFile samfile)*

    iterate over reads in `samfile` at a given list of file positions.

    .. note::
        It is usually not necessary to create an object of this class
        explicitly. It is returned as a result of call to a :meth:`AlignmentFile.fetch`.
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

        cdef uint64_t pos = self.positions[self.current_pos]
        with nogil:
            bgzf_seek(hts_get_bgzfp(self.htsfile),
                      pos,
                      0)
        self.current_pos += 1

        cdef int ret
        with nogil:
            ret = sam_read1(self.htsfile,
                            self.samfile.header,
                            self.b)
        return ret

    def __next__(self):
        cdef int ret = self.cnext()
        if (ret >= 0):
            return makeAlignedSegment(self.b, self.samfile)
        elif (ret == -2):
            raise IOError('truncated file')
        else:
            raise StopIteration


cdef int __advance_nofilter(void *data, bam1_t *b):
    '''advance without any read filtering.
    '''
    cdef __iterdata * d
    d = <__iterdata*>data
    cdef int ret
    with nogil:
        ret = sam_itr_next(d.htsfile, d.iter, b)
    return ret


cdef int __advance_all(void *data, bam1_t *b):
    '''only use reads for pileup passing basic
    filters:

    BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
    '''

    cdef __iterdata * d
    cdef mask = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP
    d = <__iterdata*>data
    cdef int ret
    with nogil:
        ret = sam_itr_next(d.htsfile, d.iter, b)
    while ret >= 0 and b.core.flag & mask:
        with nogil:
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

    cdef int ret
    cdef int skip = 0
    cdef int q
    cdef int is_cns = 1
    cdef int is_nobaq = 0
    cdef int capQ_thres = 0

    with nogil:
        ret = sam_itr_next(d.htsfile, d.iter, b)

    # reload sequence
    if d.fastafile != NULL and b.core.tid != d.tid:
        if d.seq != NULL:
            free(d.seq)
        d.tid = b.core.tid
        with nogil:
            d.seq = faidx_fetch_seq(
                d.fastafile,
                d.header.target_name[d.tid],
                0, MAX_POS,
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

        with nogil:
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
    :class:`~pysam.PileupColumn` for ``n`` columns, but each object in
    ``result`` will contain the same information.

    The desired behaviour can be achieved by list comprehension::

       result = [ x.pileups() for x in f.pileup() ]

    ``result`` will be a list of ``n`` lists of objects of type
    :class:`~pysam.PileupRead`.

    If the iterator is associated with a :class:`~pysam.Fastafile` using the
    :meth:`addReference` method, then the iterator will export the
    current sequence via the methods :meth:`getSequence` and
    :meth:`seq_len`.

    Optional kwargs to the iterator:

    stepper
       The stepper controls how the iterator advances.

       Valid values are None, "all" (default), "nofilter" or "samtools".

       See AlignmentFile.pileup for description.
    
    fastafile
       A :class:`~pysam.FastaFile` object

    max_depth
       maximum read depth. The default is 8000.

    '''

    def __cinit__( self, AlignmentFile samfile, **kwargs ):
        self.samfile = samfile
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
        # do not release gil here because of call-backs
        self.plp = bam_plp_auto(self.pileup_iter,
                                &self.tid,
                                &self.pos,
                                &self.n_plp)

    cdef char * getSequence(self):
        '''return current reference sequence underlying the iterator.
        '''
        return self.iterdata.seq

    property seq_len:
        '''current sequence length.'''
        def __get__(self):
            return self.iterdata.seq_len

    def addReference(self, Fastafile fastafile):
       '''
       add reference sequences in `fastafile` to iterator.'''
       self.fastafile = fastafile
       if self.iterdata.seq != NULL:
           free(self.iterdata.seq)
       self.iterdata.tid = -1
       self.iterdata.fastafile = self.fastafile.fastafile

    def hasReference(self):
        '''
        return true if iterator is associated with a reference'''
        return self.fastafile

    cdef setMask(self, mask):
        '''set masking flag in iterator.

        reads with bits set in `mask` will be skipped.
        '''
        raise NotImplementedError()
        # self.mask = mask
        # bam_plp_set_mask( self.pileup_iter, self.mask )

    cdef setupIteratorData( self,
                            int tid,
                            int start,
                            int end,
                            int multiple_iterators=0 ):
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
            with nogil:
                self.pileup_iter = bam_plp_init(
                    <bam_plp_auto_f>&__advance_all,
                    &self.iterdata)
        elif self.stepper == "nofilter":
            with nogil:
                self.pileup_iter = bam_plp_init(
                    <bam_plp_auto_f>&__advance_nofilter,
                    &self.iterdata)
        elif self.stepper == "samtools":
            with nogil:
                self.pileup_iter = bam_plp_init(
                    <bam_plp_auto_f>&__advance_snpcalls,
                    &self.iterdata)
        else:
            raise ValueError(
                "unknown stepper option `%s` in IteratorColumn" % self.stepper)

        if self.max_depth:
            with nogil:
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
            if self.iterdata.seq != NULL:
                free(self.iterdata.seq)
            self.iterdata.seq = NULL
            self.iterdata.tid = -1

        # self.pileup_iter = bam_plp_init( &__advancepileup, &self.iterdata )
        with nogil:
            bam_plp_reset(self.pileup_iter)

    cdef _free_pileup_iter(self):
        '''free the memory alloc'd by bam_plp_init.

        This is needed before setupIteratorData allocates
        another pileup_iter, or else memory will be lost.
        '''
        if self.pileup_iter != <bam_plp_t>NULL:
            with nogil:
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
                  int end = MAX_POS,
                  int truncate = False,
                  **kwargs ):

        # initialize iterator
        self.setupIteratorData(tid, start, end, 1)
        self.start = start
        self.end = end
        self.truncate = truncate

    def __next__(self):

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
                                   self.n_plp,
                                   self.samfile)


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
        self.setupIteratorData(self.tid, 0, MAX_POS, 1)

    def __next__(self):

        while 1:
            self.cnext()

            if self.n_plp < 0:
                raise ValueError("error during iteration" )

            # return result, if within same reference
            if self.plp != NULL:
                return makePileupColumn(&self.plp,
                                        self.tid,
                                        self.pos,
                                        self.n_plp,
                                        self.samfile)
                
            # otherwise, proceed to next reference or stop
            self.tid += 1
            if self.tid < self.samfile.nreferences:
                self.setupIteratorData(self.tid, 0, MAX_POS, 0)
            else:
                raise StopIteration


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
    """*(AlignmentFile samfile, multiple_iterators=True)

    Index a Sam/BAM-file by query name while keeping the
    original sort order intact.

    The index is kept in memory and can be substantial.

    By default, the file is re-openend to avoid conflicts if multiple
    operators work on the same file. Set `multiple_iterators` = False
    to not re-open `samfile`.

    Parameters
    ----------

    samfile : AlignmentFile
        File to be indexed.

    multiple_iterators : bool
        Flag indicating whether the file should be reopened. Reopening prevents
        existing iterators being affected by the indexing.

    """

    def __init__(self, AlignmentFile samfile, int multiple_iterators=True):
        cdef char *cfilename

        # makes sure that samfile stays alive as long as this
        # object is alive.
        self.samfile = samfile

        assert samfile.is_bam, "can only IndexReads on bam files"

        # multiple_iterators the file - note that this makes the iterator
        # slow and causes pileup to slow down significantly.
        if multiple_iterators:
            cfilename = samfile._filename
            with nogil:
                self.htsfile = hts_open(cfilename, 'r')
            assert self.htsfile != NULL
            # read header - required for accurate positioning
            with nogil:
                self.header = sam_hdr_read(self.htsfile)
            self.owns_samfile = True
        else:
            self.htsfile = self.samfile.htsfile
            self.header = self.samfile.header
            self.owns_samfile = False

    def build(self):
        '''build the index.'''

        self.index = collections.defaultdict(list)

        # this method will start indexing from the current file
        # position if you decide
        cdef int ret = 1
        cdef bam1_t * b = <bam1_t*>calloc(1, sizeof( bam1_t))

        cdef uint64_t pos

        while ret > 0:
            with nogil:
                pos = bgzf_tell(hts_get_bgzfp(self.htsfile))
                ret = sam_read1(self.htsfile,
                                self.samfile.header,
                                b)
            if ret > 0:
                qname = charptr_to_str(pysam_bam_get_qname(b))
                self.index[qname].append(pos)

        bam_destroy1(b)

    def find(self, query_name):
        '''find `query_name` in index.

        Returns
        -------

        IteratorRowSelection
            Returns an iterator over all reads with query_name.

        Raises
        ------
        
        KeyError
            if the `query_name` is not in the index.

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

__all__ = [
    "AlignmentFile",
    "IteratorRow",
    "IteratorColumn",
    "IndexedReads"]
