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
# class AlignmentHeader manage SAM/BAM/CRAM header data
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
# class IteratorColumnAll
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
try:
    from collections.abc import Sequence, Mapping  # noqa
except ImportError:
    from collections import Sequence, Mapping  # noqa
import re
import warnings
import array
from libc.errno  cimport errno, EPIPE
from libc.string cimport strcmp, strpbrk, strerror
from libc.stdint cimport INT32_MAX

from cpython cimport array as c_array

from pysam.libcutils cimport force_bytes, force_str, charptr_to_str
from pysam.libcutils cimport encode_filename, from_string_and_size
from pysam.libcalignedsegment cimport makeAlignedSegment, makePileupColumn
from pysam.libchtslib cimport HTSFile, hisremote

from io import StringIO

cimport cython


__all__ = [
    "AlignmentFile",
    "AlignmentHeader",
    "IteratorRow",
    "IteratorColumn",
    "IndexedReads"]

IndexStats = collections.namedtuple("IndexStats",
                                    ("contig",
                                     "mapped",
                                     "unmapped",
                                     "total"))

########################################################
## global variables
# maximum genomic coordinace
# for some reason, using 'int' causes overflow
cdef int MAX_POS = (1 << 31) - 1

# valid types for SAM headers
VALID_HEADER_TYPES = {"HD" : Mapping,
                      "SQ" : Sequence,
                      "RG" : Sequence,
                      "PG" : Sequence,
                      "CO" : Sequence}

# order of records within SAM headers
VALID_HEADERS = ("HD", "SQ", "RG", "PG", "CO")

# default type conversions within SAM header records
KNOWN_HEADER_FIELDS = {"HD" : {"VN" : str, "SO" : str, "GO" : str},
                       "SQ" : {"SN" : str, "LN" : int, "AS" : str,
                               "M5" : str, "SP" : str, "UR" : str,
                               "AH" : str,},
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
                               "UR", "SP", "AH"),
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


cdef AlignmentHeader makeAlignmentHeader(bam_hdr_t *hdr):
    if not hdr:
        raise ValueError('cannot create AlignmentHeader, received NULL pointer')

    # check: is AlignmetHeader.__cinit__ called?
    cdef AlignmentHeader header = AlignmentHeader.__new__(AlignmentHeader)
    header.ptr = hdr

    return header

def read_failure_reason(code):
    if code == -2:
        return 'truncated file'
    else:
        return "error {} while reading file".format(code)


# the following should be class-method for VariantHeader, but cdef @classmethods
# are not implemented in cython.
cdef int fill_AlignmentHeader_from_list(bam_hdr_t *dest,
                                        reference_names,
                                        reference_lengths,
                                        add_sq_text=True,
                                        text=None) except -1:
    """build header from list of reference names and lengths.
    """

cdef class AlignmentHeader(object):
    """header information for a :class:`AlignmentFile` object

    Parameters
    ----------
    header_dict : dict
        build header from a multi-level dictionary. The
        first level are the four types ('HD', 'SQ', ...). The second
        level are a list of lines, with each line being a list of
        tag-value pairs. The header is constructed first from all the
        defined fields, followed by user tags in alphabetical
        order. Alternatively, an :class:`~pysam.AlignmentHeader`
        object can be passed directly.

    text : string
        use the string provided as the header

    reference_names : list
        see reference_lengths

    reference_lengths : list
        build header from list of chromosome names and lengths.  By
        default, 'SQ' and 'LN' tags will be added to the header
        text. This option can be changed by unsetting the flag
        `add_sq_text`.

    add_sq_text : bool
        do not add 'SQ' and 'LN' tags to header. This option permits
        construction :term:`SAM` formatted files without a header.

    """

    # See makeVariantHeader for C constructor
    def __cinit__(self):
        self.ptr = NULL

    # Python constructor
    def __init__(self):
        self.ptr = bam_hdr_init()
        if self.ptr is NULL:
            raise MemoryError("could not create header")

    @classmethod
    def _from_text_and_lengths(cls, text, reference_names, reference_lengths):

        cdef AlignmentHeader self = AlignmentHeader()
        cdef char *ctext
        cdef int l_text
        cdef int n, x
        if text is not None:
            btext = force_bytes(text)
            ctext = btext
            l_text = len(btext)
            self.ptr.text = <char*>calloc(l_text + 1, sizeof(char))
            if self.ptr.text == NULL:
                raise MemoryError("could not allocate {} bytes".format(l_text + 1), sizeof(char))
            self.ptr.l_text = l_text
            memcpy(self.ptr.text, ctext, l_text + 1)

        if reference_names and reference_lengths:
            reference_names = [force_bytes(ref) for ref in reference_names]

            self.ptr.n_targets = len(reference_names)

            n = sum([len(reference_names) + 1])
            self.ptr.target_name = <char**>calloc(n, sizeof(char*))
            if self.ptr.target_name == NULL:
                raise MemoryError("could not allocate {} bytes".format(n, sizeof(char *)))

            self.ptr.target_len = <uint32_t*>calloc(n, sizeof(uint32_t))
            if self.ptr.target_len == NULL:
                raise MemoryError("could not allocate {} bytes".format(n, sizeof(uint32_t)))

            for x from 0 <= x < self.ptr.n_targets:
                self.ptr.target_len[x] = reference_lengths[x]
                name = reference_names[x]
                self.ptr.target_name[x] = <char*>calloc(len(name) + 1, sizeof(char))
                if self.ptr.target_name[x] == NULL:
                    raise MemoryError("could not allocate {} bytes".format(len(name) + 1, sizeof(char)))
                strncpy(self.ptr.target_name[x], name, len(name))

        return self

    @classmethod
    def from_text(cls, text):

        reference_names, reference_lengths = [], []
        for line in text.splitlines():
            if line.startswith("@SQ"):
                fields = dict([x.split(":", 1) for x in line.split("\t")[1:]])
                try:
                    reference_names.append(fields["SN"])
                    reference_lengths.append(int(fields["LN"]))
                except KeyError:
                    raise KeyError("incomplete sequence information in '%s'" % str(fields))
                except ValueError:
                    raise ValueError("wrong sequence information in '%s'" % str(fields))

        return cls._from_text_and_lengths(text, reference_names, reference_lengths)

    @classmethod
    def from_dict(cls, header_dict):

        cdef list lines = []
        # first: defined tags
        for record in VALID_HEADERS:
            if record in header_dict:
                data = header_dict[record]
                if not isinstance(data, VALID_HEADER_TYPES[record]):
                    raise ValueError(
                        "invalid type for record %s: %s, expected %s".format(
                            record, type(data), VALID_HEADER_TYPES[record]))
                if isinstance(data, Mapping):
                    lines.append(build_header_line(data, record))
                else:
                    for fields in header_dict[record]:
                        lines.append(build_header_line(fields, record))

        # then: user tags (lower case), sorted alphabetically
        for record, data in sorted(header_dict.items()):
            if record in VALID_HEADERS:
                continue
            if isinstance(data, Mapping):
                lines.append(build_header_line(data, record))
            else:
                for fields in header_dict[record]:
                    lines.append(build_header_line(fields, record))

        text = "\n".join(lines) + "\n"

        reference_names, reference_lengths = [], []
        if "SQ" in header_dict:
            for fields in header_dict["SQ"]:
                try:
                    reference_names.append(fields["SN"])
                    reference_lengths.append(fields["LN"])
                except KeyError:
                    raise KeyError("incomplete sequence information in '%s'" % str(fields))

        return cls._from_text_and_lengths(text, reference_names, reference_lengths)

    @classmethod
    def from_references(cls, reference_names, reference_lengths, text=None, add_sq_text=True):

        if len(reference_names) != len(reference_lengths):
            raise ValueError("number of reference names and lengths do not match")

        # optionally, if there is no text, add a SAM compatible header to output file.
        if text is None and add_sq_text:
            text = "".join(["@SQ\tSN:{}\tLN:{}\n".format(x, y) for x, y in zip(
                reference_names, reference_lengths)])

        return cls._from_text_and_lengths(text, reference_names, reference_lengths)

    def __dealloc__(self):
        bam_hdr_destroy(self.ptr)
        self.ptr = NULL

    def __bool__(self):
        return self.ptr != NULL

    def copy(self):
        return makeAlignmentHeader(bam_hdr_dup(self.ptr))

    property nreferences:
        """int with the number of :term:`reference` sequences in the file.

        This is a read-only attribute."""
        def __get__(self):
            return self.ptr.n_targets

    property references:
        """tuple with the names of :term:`reference` sequences. This is a
        read-only attribute"""
        def __get__(self):
            t = []
            cdef int x
            for x in range(self.ptr.n_targets):
                t.append(charptr_to_str(self.ptr.target_name[x]))
            return tuple(t)

    property lengths:
        """tuple of the lengths of the :term:`reference` sequences. This is a
        read-only attribute. The lengths are in the same order as
        :attr:`pysam.AlignmentFile.references`
        """
        def __get__(self):
            t = []
            cdef int x
            for x in range(self.ptr.n_targets):
                t.append(self.ptr.target_len[x])
            return tuple(t)

    def _build_sequence_section(self):
        """return sequence section of header.

        The sequence section is built from the list of reference names and
        lengths stored in the BAM-file and not from any @SQ entries that
        are part of the header's text section.
        """

        cdef int x
        text = []
        for x in range(self.ptr.n_targets):
            text.append("@SQ\tSN:{}\tLN:{}\n".format(
                force_str(self.ptr.target_name[x]),
                self.ptr.target_len[x]))
        return "".join(text)

    def to_dict(self):
        """return two-level dictionary with header information from the file.

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

        If no @SQ entries are within the text section of the header,
        this will be automatically added from the reference names and
        lengths stored in the binary part of the header.
        """
        result = collections.OrderedDict()

        # convert to python string
        t = self.__str__()
        for line in t.split("\n"):
            line = line.strip(' \0')
            if not line:
                continue
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

            if VALID_HEADER_TYPES[record] == Mapping:
                if record in result:
                    raise ValueError(
                        "multiple '%s' lines are not permitted" % record)

                result[record] = x
            elif VALID_HEADER_TYPES[record] == Sequence:
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

    def as_dict(self):
        """deprecated, use :meth:`to_dict()` instead"""
        return self.to_dict()

    def get_reference_name(self, tid):
        if tid == -1:
            return None
        if not 0 <= tid < self.ptr.n_targets:
            raise ValueError("reference_id %i out of range 0<=tid<%i" %
                             (tid, self.ptr.n_targets))
        return charptr_to_str(self.ptr.target_name[tid])

    def get_reference_length(self, reference):
        cdef int tid = self.get_tid(reference)
        if tid < 0:
            raise KeyError("unknown reference {}".format(reference))
        else:
            return self.ptr.target_len[tid]

    def is_valid_tid(self, int tid):
        """
        return True if the numerical :term:`tid` is valid; False otherwise.

        Note that the unmapped tid code (-1) counts as an invalid.
        """
        return 0 <= tid < self.ptr.n_targets

    def get_tid(self, reference):
        """
        return the numerical :term:`tid` corresponding to
        :term:`reference`

        returns -1 if reference is not known.
        """
        reference = force_bytes(reference)
        tid = bam_name2id(self.ptr, reference)
        if tid < -1:
            raise ValueError('could not parse header')
        return tid

    def __str__(self):
        '''string with the full contents of the :term:`sam file` header as a
        string.

        If no @SQ entries are within the text section of the header,
        this will be automatically added from the reference names and
        lengths stored in the binary part of the header.

        See :attr:`pysam.AlignmentFile.header.to_dict()` to get a parsed
        representation of the header.
        '''
        text = from_string_and_size(self.ptr.text, self.ptr.l_text)
        if "@SQ" not in text:
            text += "\n" + self._build_sequence_section()
        return text

    # dictionary access methods, for backwards compatibility.
    def __setitem__(self, key, value):
        raise TypeError("AlignmentHeader does not support item assignment (use header.to_dict()")

    def __getitem__(self, key):
        return self.to_dict().__getitem__(key)

    def items(self):
        return self.to_dict().items()

    # PY2 compatibility
    def iteritems(self):
        return self.to_dict().items()

    def keys(self):
        return self.to_dict().keys()

    def values(self):
        return self.to_dict().values()

    def get(self, *args):
        return self.to_dict().get(*args)

    def __len__(self):
        return self.to_dict().__len__()

    def __contains__(self, key):
        return self.to_dict().__contains__(key)


cdef class AlignmentFile(HTSFile):
    """AlignmentFile(filepath_or_object, mode=None, template=None,
    reference_names=None, reference_lengths=None, text=NULL,
    header=None, add_sq_text=False, check_header=True, check_sq=True,
    reference_filename=None, filename=None, index_filename=None,
    filepath_index=None, require_index=False, duplicate_filehandle=True,
    ignore_truncation=False, threads=1)

    A :term:`SAM`/:term:`BAM`/:term:`CRAM` formatted file.

    If `filepath_or_object` is a string, the file is automatically
    opened. If `filepath_or_object` is a python File object, the
    already opened file will be used.

    If the file is opened for reading and an index exists (if file is BAM, a
    .bai file or if CRAM a .crai file), it will be opened automatically.
    `index_filename` may be specified explicitly. If the index is not named
    in the standard manner, not located in the same directory as the
    BAM/CRAM file, or is remote.  Without an index, random access via
    :meth:`~pysam.AlignmentFile.fetch` and :meth:`~pysam.AlignmentFile.pileup`
    is disabled.

    For writing, the header of a :term:`SAM` file/:term:`BAM` file can
    be constituted from several sources (see also the samtools format
    specification):

        1. If `template` is given, the header is copied from another
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
        when writing, copy header from file `template`.

    header :  dict or AlignmentHeader
        when writing, build header from a multi-level dictionary. The
        first level are the four types ('HD', 'SQ', ...). The second
        level are a list of lines, with each line being a list of
        tag-value pairs. The header is constructed first from all the
        defined fields, followed by user tags in alphabetical
        order. Alternatively, an :class:`~pysam.AlignmentHeader`
        object can be passed directly.

    text : string
        when writing, use the string provided as the header

    reference_names : list
        see reference_lengths

    reference_lengths : list
        when writing or opening a SAM file without header build header
        from list of chromosome names and lengths.  By default, 'SQ'
        and 'LN' tags will be added to the header text. This option
        can be changed by unsetting the flag `add_sq_text`.

    add_sq_text : bool
        do not add 'SQ' and 'LN' tags to header. This option permits
        construction :term:`SAM` formatted files without a header.

    add_sam_header : bool
        when outputting SAM the default is to output a header. This is
        equivalent to opening the file in 'wh' mode. If this option is
        set to False, no header will be output. To read such a file,
        set `check_header=False`.

    check_header : bool
        obsolete: when reading a SAM file, check if header is present
        (default=True)

    check_sq : bool
        when reading, check if SQ entries are present in header
        (default=True)

    reference_filename : string
        Path to a FASTA-formatted reference file. Valid only for CRAM files.
        When reading a CRAM file, this overrides both ``$REF_PATH`` and the URL
        specified in the header (``UR`` tag), which are normally used to find
        the reference.

    index_filename : string
        Explicit path to the index file.  Only needed if the index is not
        named in the standard manner, not located in the same directory as
        the BAM/CRAM file, or is remote.  An IOError is raised if the index
        cannot be found or is invalid.

    filepath_index : string
        Alias for `index_filename`.

    require_index : bool
        When reading, require that an index file is present and is valid or
        raise an IOError.  (default=False)

    filename : string
        Alternative to filepath_or_object. Filename of the file
        to be opened.

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

    format_options: list
        A list of key=value strings, as accepted by --input-fmt-option and
        --output-fmt-option in samtools.
    threads: integer
        Number of threads to use for compressing/decompressing BAM/CRAM files.
        Setting threads to > 1 cannot be combined with `ignore_truncation`.
        (Default=1)
    """

    def __cinit__(self, *args, **kwargs):
        self.htsfile = NULL
        self.filename = None
        self.mode = None
        self.threads = 1
        self.is_stream = False
        self.is_remote = False
        self.index = NULL

        if "filename" in kwargs:
            args = [kwargs["filename"]]
            del kwargs["filename"]

        self._open(*args, **kwargs)

        # allocate memory for iterator
        self.b = <bam1_t*>calloc(1, sizeof(bam1_t))
        if self.b == NULL:
            raise MemoryError("could not allocate memory of size {}".format(sizeof(bam1_t)))

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

        if not self.is_open:
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
              add_sam_header=True,
              check_header=True,
              check_sq=True,
              index_filename=None,
              filepath_index=None,
              require_index=False,
              referencenames=None,
              referencelengths=None,
              duplicate_filehandle=True,
              ignore_truncation=False,
              format_options=None,
              threads=1):
        '''open a sam, bam or cram formatted file.

        If _open is called on an existing file, the current file
        will be closed and a new file will be opened.

        '''
        cdef char *cfilename = NULL
        cdef char *creference_filename = NULL
        cdef char *cindexname = NULL
        cdef char *cmode = NULL
        cdef bam_hdr_t * hdr = NULL

        if threads > 1 and ignore_truncation:
           # This won't raise errors if reaching a truncated alignment,
           # because bgzf_mt_reader in htslib does not deal with
           # bgzf_mt_read_block returning non-zero values, contrary
           # to bgzf_read (https://github.com/samtools/htslib/blob/1.7/bgzf.c#L888)
           # Better to avoid this (for now) than to produce seemingly correct results.
           raise ValueError('Cannot add extra threads when "ignore_truncation" is True')
        self.threads = threads

        # for backwards compatibility:
        if referencenames is not None:
            reference_names = referencenames
        if referencelengths is not None:
            reference_lengths = referencelengths

        # close a previously opened file
        if self.is_open:
            self.close()

        # autodetection for read
        if mode is None:
            mode = "r"

        if add_sam_header and mode == "w":
            mode = "wh"

        assert mode in ("r", "w", "rb", "wb", "wh",
                        "wbu", "rU", "wb0",
                        "rc", "wc"), \
            "invalid file opening mode `%s`" % mode

        self.duplicate_filehandle = duplicate_filehandle

        # StringIO not supported
        if isinstance(filepath_or_object, StringIO):
            raise NotImplementedError(
                "access from StringIO objects not supported")
        # reading from a file descriptor
        elif isinstance(filepath_or_object, int):
            self.filename = filepath_or_object
            filename = None
            self.is_remote = False
            self.is_stream = True
        # reading from a File object or other object with fileno
        elif hasattr(filepath_or_object, "fileno"):
            if filepath_or_object.closed:
                raise ValueError('I/O operation on closed file')
            self.filename = filepath_or_object
            # .name can be TextIOWrapper
            try:
                filename = encode_filename(str(filepath_or_object.name))
                cfilename = filename
            except AttributeError:
                filename = None
            self.is_remote = False
            self.is_stream = True
        # what remains is a filename
        else:
            self.filename = filename = encode_filename(filepath_or_object)
            cfilename = filename
            self.is_remote = hisremote(cfilename)
            self.is_stream = self.filename == b'-'

        # for htslib, wbu seems to not work
        if mode == "wbu":
            mode = "wb0"

        self.mode = force_bytes(mode)
        self.reference_filename = reference_filename = encode_filename(
            reference_filename)

        if mode[0] == 'w':
            # open file for writing

            if not (template or header or text or (reference_names and reference_lengths)):
                raise ValueError(
                    "either supply options `template`, `header`, `text` or  both `reference_names` "
                    "and `reference_lengths` for writing")

            if template:
                # header is copied, though at the moment not strictly
                # necessary as AlignmentHeader is immutable.
                self.header = template.header.copy()
            elif isinstance(header, AlignmentHeader):
                self.header = header.copy()
            elif isinstance(header, Mapping):
                self.header = AlignmentHeader.from_dict(header)
            elif reference_names and reference_lengths:
                self.header = AlignmentHeader.from_references(
                    reference_names,
                    reference_lengths,
                    add_sq_text=add_sq_text,
                    text=text)
            elif text:
                self.header = AlignmentHeader.from_text(text)
            else:
                raise ValueError("not enough information to construct header. Please provide template, "
                                 "header, text or reference_names/reference_lengths")
            self.htsfile = self._open_htsfile()

            if self.htsfile == NULL:
                if errno:
                    raise IOError(errno, "could not open alignment file `{}`: {}".format(
                        force_str(filename),
                        force_str(strerror(errno))))
                else:
                    raise ValueError("could not open alignment file `{}`".format(force_str(filename)))
            if format_options and len(format_options):
                self.add_hts_options(format_options)
            # set filename with reference sequences. If no filename
            # is given, the CRAM reference arrays will be built from
            # the @SQ header in the header
            if "c" in mode and reference_filename:
                if (hts_set_fai_filename(self.htsfile, self.reference_filename) != 0):
                    raise ValueError("failure when setting reference filename")

            # write header to htsfile
            if "b" in mode or "c" in mode or "h" in mode:
                hdr = self.header.ptr
                with nogil:
                    sam_hdr_write(self.htsfile, hdr)

        elif mode[0] == "r":
            # open file for reading
            self.htsfile = self._open_htsfile()

            if self.htsfile == NULL:
                if errno:
                    raise IOError(errno, "could not open alignment file `{}`: {}".format(force_str(filename),
                                  force_str(strerror(errno))))
                else:
                    raise ValueError("could not open alignment file `{}`".format(force_str(filename)))

            if self.htsfile.format.category != sequence_data:
                raise ValueError("file does not contain alignment data")

            if format_options and len(format_options):
                self.add_hts_options(format_options)

            self.check_truncation(ignore_truncation)

            # bam/cram files require a valid header
            if self.is_bam or self.is_cram:
                with nogil:
                    hdr = sam_hdr_read(self.htsfile)
                if hdr == NULL:
                    raise ValueError(
                        "file does not have a valid header (mode='%s') "
                        "- is it BAM/CRAM format?" % mode)
                self.header = makeAlignmentHeader(hdr)
            else:
                # in sam files a header is optional. If not given,
                # user may provide reference names and lengths to built
                # an on-the-fly header.
                if reference_names and reference_lengths:
                    # build header from a target names and lengths
                    self.header = AlignmentHeader.from_references(
                        reference_names=reference_names,
                        reference_lengths=reference_lengths,
                        add_sq_text=add_sq_text,
                        text=text)
                else:
                    with nogil:
                        hdr = sam_hdr_read(self.htsfile)
                    if hdr == NULL:
                        raise ValueError(
                            "SAM? file does not have a valid header (mode='%s'), "
                            "please provide reference_names and reference_lengths")
                    self.header = makeAlignmentHeader(hdr)

            # set filename with reference sequences
            if self.is_cram and reference_filename:
                creference_filename = self.reference_filename
                hts_set_opt(self.htsfile,
                            CRAM_OPT_REFERENCE,
                            creference_filename)

            if check_sq and self.header.nreferences == 0:
                raise ValueError(
                    ("file has no sequences defined (mode='%s') - "
                     "is it SAM/BAM format? Consider opening with "
                     "check_sq=False") % mode)

            if self.is_bam or self.is_cram:
                self.index_filename = index_filename or filepath_index
                if self.index_filename:
                    cindexname = bfile_name = encode_filename(self.index_filename)

                if cfilename or cindexname:
                    with nogil:
                        self.index = sam_index_load2(self.htsfile, cfilename, cindexname)

                    if not self.index and (cindexname or require_index):
                        if errno:
                            raise IOError(errno, force_str(strerror(errno)))
                        else:
                            raise IOError('unable to open index file `%s`' % self.index_filename)

                elif require_index:
                    raise IOError('unable to open index file')

                # save start of data section
                if not self.is_stream:
                    self.start_offset = self.tell()

    def fetch(self,
              contig=None,
              start=None,
              stop=None,
              region=None,
              tid=None,
              until_eof=False,
              multiple_iterators=False,
              reference=None,
              end=None):
        """fetch reads aligned in a :term:`region`.

        See :meth:`~pysam.HTSFile.parse_region` for more information
        on how genomic regions can be specified. :term:`reference` and
        `end` are also accepted for backward compatibility as synonyms
        for :term:`contig` and `stop`, respectively.

        Without a `contig` or `region` all mapped reads in the file
        will be fetched. The reads will be returned ordered by reference
        sequence, which will not necessarily be the order within the
        file. This mode of iteration still requires an index. If there is
        no index, use `until_eof=True`.

        If only `contig` is set, all reads aligned to `contig`
        will be fetched.

        A :term:`SAM` file does not allow random access. If `region`
        or `contig` are given, an exception is raised.

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

		An iterator over a collection of reads. : IteratorRow

        Raises
        ------

        ValueError
            if the genomic coordinates are out of range or invalid or the
            file does not permit random access to genomic coordinates.

        """
        cdef int rtid, rstart, rstop, has_coord

        if not self.is_open:
            raise ValueError( "I/O operation on closed file" )

        has_coord, rtid, rstart, rstop = self.parse_region(
            contig, start, stop, region, tid,
            end=end, reference=reference)

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
                    self, rtid, rstart, rstop,
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
                    "fetching by region is not available for SAM files")

            if multiple_iterators == True:
                raise ValueError(
                    "multiple iterators not implemented for SAM files")

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

		an iterator over a collection of reads : IteratorRowHead

        '''
        return IteratorRowHead(self, n,
                               multiple_iterators=multiple_iterators)

    def mate(self, AlignedSegment read):
        '''return the mate of :class:`pysam.AlignedSegment` `read`.

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

        the mate : AlignedSegment

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
               contig=None,
               start=None,
               stop=None,
               region=None,
               reference=None,
               end=None,
               **kwargs):
        """perform a :term:`pileup` within a :term:`region`. The region is
        specified by :term:`contig`, `start` and `stop` (using
        0-based indexing).  :term:`reference` and `end` are also accepted for
        backward compatibility as synonyms for :term:`contig` and `stop`,
        respectively.  Alternatively, a samtools 'region' string
        can be supplied.

        Without 'contig' or 'region' all reads will be used for the
        pileup. The reads will be returned ordered by
        :term:`contig` sequence, which will not necessarily be the
        order within the file.

        Note that :term:`SAM` formatted files do not allow random
        access.  In these files, if a 'region' or 'contig' are
        given an exception is raised.

        .. note::

            'all' reads which overlap the region are returned. The
            first base returned will be the first base of the first
            read 'not' necessarily the first base of the region used
            in the query.

        Parameters
        ----------

        truncate : bool

           By default, the samtools pileup engine outputs all reads
           overlapping a region. If truncate is True and a region is
           given, only columns in the exact region specified are
           returned.

        max_depth : int
           Maximum read depth permitted. The default limit is '8000'.

        stepper : string
           The stepper controls how the iterator advances.
           Possible options for the stepper are

           ``all``
              skip reads in which any of the following flags are set:
              BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP

           ``nofilter``
              uses every single read turning off any filtering.

           ``samtools``
              same filter and read processing as in samtools
              pileup. For full compatibility, this requires a
              'fastafile' to be given. The following options all pertain
              to filtering of the ``samtools`` stepper.

        fastafile : :class:`~pysam.FastaFile` object.

           This is required for some of the steppers.

        ignore_overlaps: bool

           If set to True, detect if read pairs overlap and only take
           the higher quality base. This is the default.

        flag_filter : int

           ignore reads where any of the bits in the flag are set. The default is
           BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP.

        flag_require : int

           only use reads where certain flags are set. The default is 0.

        ignore_orphans: bool

            ignore orphans (paired reads that are not in a proper pair).
            The default is to ignore orphans.

        min_base_quality: int

           Minimum base quality. Bases below the minimum quality will
           not be output. The default is 13.

        adjust_capq_threshold: int

           adjust mapping quality. The default is 0 for no
           adjustment. The recommended value for adjustment is 50.

        min_mapping_quality : int

           only use reads above a minimum mapping quality. The default is 0.

        compute_baq: bool

           re-alignment computing per-Base Alignment Qualities (BAQ). The
           default is to do re-alignment. Realignment requires a reference
           sequence. If none is present, no realignment will be performed.

        redo_baq: bool

           recompute per-Base Alignment Quality on the fly ignoring
           existing base qualities. The default is False (use existing
           base qualities).

        Returns
        -------

        an iterator over genomic positions. : IteratorColumn

        """
        cdef int rtid, has_coord
        cdef int32_t rstart, rstop

        if not self.is_open:
            raise ValueError("I/O operation on closed file")

        has_coord, rtid, rstart, rstop = self.parse_region(
            contig, start, stop, region, reference=reference, end=end)

        if has_coord:
            if not self.has_index():
                raise ValueError("no index available for pileup")

            return IteratorColumnRegion(self,
                                        tid=rtid,
                                        start=rstart,
                                        stop=rstop,
                                        **kwargs)
        else:
            if self.has_index():
                return IteratorColumnAllRefs(self, **kwargs)
            else:
                return IteratorColumnAll(self, **kwargs)

    def count(self,
              contig=None,
              start=None,
              stop=None,
              region=None,
              until_eof=False,
              read_callback="nofilter",
              reference=None,
              end=None):
        '''count the number of reads in :term:`region`

        The region is specified by :term:`contig`, `start` and `stop`.
        :term:`reference` and `end` are also accepted for backward
        compatibility as synonyms for :term:`contig` and `stop`,
        respectively.  Alternatively, a `samtools`_ :term:`region`
        string can be supplied.

        A :term:`SAM` file does not allow random access and if
        `region` or `contig` are given, an exception is raised.

        Parameters
        ----------

        contig : string
            reference_name of the genomic region (chromosome)

        start : int
            start of the genomic region (0-based inclusive)

        stop : int
            end of the genomic region (0-based exclusive)

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

        reference : string
            backward compatible synonym for `contig`

        end : int
            backward compatible synonym for `stop`

        Raises
        ------

        ValueError
            if the genomic coordinates are out of range or invalid.

        '''
        cdef AlignedSegment read
        cdef long counter = 0

        if not self.is_open:
            raise ValueError("I/O operation on closed file")

        cdef int filter_method = 0
        if read_callback == "all":
            filter_method = 1
        elif read_callback == "nofilter":
            filter_method = 2

        for read in self.fetch(contig=contig,
                               start=start,
                               stop=stop,
                               reference=reference,
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
                       contig,
                       start=None,
                       stop=None,
                       region=None,
                       quality_threshold=15,
                       read_callback='all',
                       reference=None,
                       end=None):
        """count the coverage of genomic positions by reads in :term:`region`.

        The region is specified by :term:`contig`, `start` and `stop`.
        :term:`reference` and `end` are also accepted for backward
        compatibility as synonyms for :term:`contig` and `stop`,
        respectively.  Alternatively, a `samtools`_ :term:`region`
        string can be supplied.  The coverage is computed per-base [ACGT].

        Parameters
        ----------

        contig : string
            reference_name of the genomic region (chromosome)

        start : int
            start of the genomic region (0-based inclusive). If not
            given, count from the start of the chromosome.

        stop : int
            end of the genomic region (0-based exclusive). If not given,
            count to the end of the chromosome.

        region : string
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

        reference : string
            backward compatible synonym for `contig`

        end : int
            backward compatible synonym for `stop`

        Raises
        ------

        ValueError
            if the genomic coordinates are out of range or invalid.

        Returns
        -------

        four array.arrays of the same length in order A C G T : tuple

        """

        cdef uint32_t contig_length = self.get_reference_length(contig)
        cdef int _start = start if start is not None else 0
        cdef int _stop = stop if stop is not None else contig_length
        _stop = _stop if _stop < contig_length else contig_length

        if _stop == _start:
            raise ValueError("interval of size 0")
        if _stop < _start:
            raise ValueError("interval of size less than 0")

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

        cdef int _threshold = quality_threshold or 0
        for read in self.fetch(contig=contig,
                               reference=reference,
                               start=start,
                               stop=stop,
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
            if seq is None:
                continue
            quality = read.query_qualities

            for qpos, refpos in read.get_aligned_pairs(True):
                if qpos is not None and refpos is not None and \
                   _start <= refpos < _stop:

                    # only check base quality if _threshold > 0
                    if (_threshold and quality and quality[qpos] >= _threshold) or not _threshold:
                        if seq[qpos] == 'A':
                            count_a.data.as_ulongs[refpos - _start] += 1
                        if seq[qpos] == 'C':
                            count_c.data.as_ulongs[refpos - _start] += 1
                        if seq[qpos] == 'G':
                            count_g.data.as_ulongs[refpos - _start] += 1
                        if seq[qpos] == 'T':
                            count_t.data.as_ulongs[refpos - _start] += 1

        return count_a, count_c, count_g, count_t

    def find_introns_slow(self, read_iterator):
        """Return a dictionary {(start, stop): count}
        Listing the intronic sites in the reads (identified by 'N' in the cigar strings),
        and their support ( = number of reads ).

        read_iterator can be the result of a .fetch(...) call.
        Or it can be a generator filtering such reads. Example
        samfile.find_introns((read for read in samfile.fetch(...) if read.is_reverse)
        """
        res = collections.Counter()
        for r in read_iterator:
            if 'N' in r.cigarstring:
                last_read_pos = False
                for read_loc, genome_loc in r.get_aligned_pairs():
                    if read_loc is None and last_read_pos:
                        start = genome_loc
                    elif read_loc and last_read_pos is None:
                        stop = genome_loc  # we are right exclusive ,so this is correct
                        res[(start, stop)] += 1
                        del start
                        del stop
                    last_read_pos = read_loc
        return res

    def find_introns(self, read_iterator):
        """Return a dictionary {(start, stop): count}
        Listing the intronic sites in the reads (identified by 'N' in the cigar strings),
        and their support ( = number of reads ).

        read_iterator can be the result of a .fetch(...) call.
        Or it can be a generator filtering such reads. Example
        samfile.find_introns((read for read in samfile.fetch(...) if read.is_reverse)
        """
        cdef:
            uint32_t base_position, junc_start, nt
            int op
            AlignedSegment r
            int BAM_CREF_SKIP = 3 #BAM_CREF_SKIP

        res = collections.Counter()

        match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
        for r in read_iterator:
            base_position = r.pos
            cigar = r.cigartuples
            if cigar is None:
                continue

            for op, nt in cigar:
                if op in match_or_deletion:
                    base_position += nt
                elif op == BAM_CREF_SKIP:
                    junc_start = base_position
                    base_position += nt
                    res[(junc_start, base_position)] += 1
        return res


    def close(self):
        '''closes the :class:`pysam.AlignmentFile`.'''

        if self.htsfile == NULL:
            return

        if self.index != NULL:
            hts_idx_destroy(self.index)
            self.index = NULL

        cdef int ret = hts_close(self.htsfile)
        self.htsfile = NULL

        self.header = None

        if ret < 0:
            global errno
            if errno == EPIPE:
                errno = 0
            else:
                raise IOError(errno, force_str(strerror(errno)))

    def __dealloc__(self):
        cdef int ret = 0

        if self.index != NULL:
            hts_idx_destroy(self.index)
            self.index = NULL

        if self.htsfile != NULL:
            ret = hts_close(self.htsfile)
            self.htsfile = NULL

        self.header = None

        if self.b:
            bam_destroy1(self.b)
            self.b = NULL

        if ret < 0:
            global errno
            if errno == EPIPE:
                errno = 0
            else:
                raise IOError(errno, force_str(strerror(errno)))

    cpdef int write(self, AlignedSegment read) except -1:
        '''
        write a single :class:`pysam.AlignedSegment` to disk.

        Raises:
            ValueError
                if the writing failed

        Returns:
            int :
                the number of bytes written. If the file is closed,
                this will be 0.
        '''
        if not self.is_open:
            return 0

        if self.header.ptr.n_targets <= read._delegate.core.tid:
            raise ValueError(
                "AlignedSegment refers to reference number {} that "
                "is larger than the number of references ({}) in the header".format(
                    read._delegate.core.tid, self.header.ptr.n_targets))

        cdef int ret
        with nogil:
            ret = sam_write1(self.htsfile,
                             self.header.ptr,
                             read._delegate)

        # kbj: Still need to raise an exception with except -1. Otherwise
        #      when ret == -1 we get a "SystemError: error return without
        #      exception set".
        if ret < 0:
            raise IOError(
            "sam_write1 failed with error code {}".format(ret))

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
    property mapped:
        """int with total number of mapped alignments according to the
        statistics recorded in the index. This is a read-only
        attribute.
        (This will be 0 for a CRAM file indexed by a .crai index, as that
        index format does not record these statistics.)
        """
        def __get__(self):
            self.check_index()
            cdef int tid
            cdef uint64_t total = 0
            cdef uint64_t mapped, unmapped
            for tid from 0 <= tid < self.header.nreferences:
                with nogil:
                    hts_idx_get_stat(self.index, tid, &mapped, &unmapped)
                total += mapped
            return total

    property unmapped:
        """int with total number of unmapped reads according to the statistics
        recorded in the index. This number of reads includes the number of reads
        without coordinates. This is a read-only attribute.
        (This will be 0 for a CRAM file indexed by a .crai index, as that
        index format does not record these statistics.)
        """
        def __get__(self):
            self.check_index()
            cdef int tid
            cdef uint64_t total = hts_idx_get_n_no_coor(self.index)
            cdef uint64_t mapped, unmapped
            for tid from 0 <= tid < self.header.nreferences:
                with nogil:
                    hts_idx_get_stat(self.index, tid, &mapped, &unmapped)
                total += unmapped
            return total

    property nocoordinate:
        """int with total number of reads without coordinates according to the
        statistics recorded in the index, i.e., the statistic printed for "*"
        by the ``samtools idxstats`` command. This is a read-only attribute.
        (This will be 0 for a CRAM file indexed by a .crai index, as that
        index format does not record these statistics.)
        """
        def __get__(self):
            self.check_index()
            cdef uint64_t n
            with nogil:
                n = hts_idx_get_n_no_coor(self.index)
            return n

    def get_index_statistics(self):
        """return statistics about mapped/unmapped reads per chromosome as
        they are stored in the index, similarly to the statistics printed
        by the ``samtools idxstats`` command.

        CRAI indexes do not record these statistics, so for a CRAM file
        with a .crai index the returned statistics will all be 0.

        Returns:
            list :
                a list of records for each chromosome. Each record has the
                attributes 'contig', 'mapped', 'unmapped' and 'total'.
        """

        self.check_index()
        cdef int tid
        cdef uint64_t mapped, unmapped
        results = []
        # TODO: use header
        for tid from 0 <= tid < self.nreferences:
            with nogil:
                hts_idx_get_stat(self.index, tid, &mapped, &unmapped)
            results.append(
                IndexStats._make((
                    self.get_reference_name(tid),
                    mapped,
                    unmapped,
                    mapped + unmapped)))

        return results

    ###############################################################
    ## file-object like iterator access
    ## note: concurrent access will cause errors (see IteratorRow
    ## and multiple_iterators)
    ## Possible solutions: deprecate or open new file handle
    def __iter__(self):
        if not self.is_open:
            raise ValueError("I/O operation on closed file")

        if not self.is_bam and self.header.nreferences == 0:
            raise NotImplementedError(
                "can not iterate over samfile without header")
        return self

    cdef bam1_t * getCurrent(self):
        return self.b

    cdef int cnext(self):
        '''
        cversion of iterator. Used by :class:`pysam.AlignmentFile.IteratorColumn`.
        '''
        cdef int ret
        cdef bam_hdr_t * hdr = self.header.ptr
        with nogil:
            ret = sam_read1(self.htsfile,
                            hdr,
                            self.b)
        return ret

    def __next__(self):
        cdef int ret = self.cnext()
        if ret >= 0:
            return makeAlignedSegment(self.b, self.header)
        elif ret == -1:
            raise StopIteration
        else:
            raise IOError(read_failure_reason(ret))

    ###########################################
    # methods/properties referencing the header
    def is_valid_tid(self, int tid):
        """
        return True if the numerical :term:`tid` is valid; False otherwise.

        Note that the unmapped tid code (-1) counts as an invalid.
        """
        if self.header is None:
            raise ValueError("header not available in closed files")
        return self.header.is_valid_tid(tid)

    def get_tid(self, reference):
        """
        return the numerical :term:`tid` corresponding to
        :term:`reference`

        returns -1 if reference is not known.
        """
        if self.header is None:
            raise ValueError("header not available in closed files")
        return self.header.get_tid(reference)

    def get_reference_name(self, tid):
        """
        return :term:`reference` name corresponding to numerical :term:`tid`
        """
        if self.header is None:
            raise ValueError("header not available in closed files")
        return self.header.get_reference_name(tid)

    def get_reference_length(self, reference):
        """
        return :term:`reference` length corresponding to numerical :term:`tid`
        """
        if self.header is None:
            raise ValueError("header not available in closed files")
        return self.header.get_reference_length(reference)

    property nreferences:
        """int with the number of :term:`reference` sequences in the file.
        This is a read-only attribute."""
        def __get__(self):
            if self.header:
                return self.header.nreferences
            else:
                raise ValueError("header not available in closed files")

    property references:
        """tuple with the names of :term:`reference` sequences. This is a
        read-only attribute"""
        def __get__(self):
            if self.header:
                return self.header.references
            else:
                raise ValueError("header not available in closed files")

    property lengths:
        """tuple of the lengths of the :term:`reference` sequences. This is a
        read-only attribute. The lengths are in the same order as
        :attr:`pysam.AlignmentFile.references`

        """
        def __get__(self):
            if self.header:
                return self.header.lengths
            else:
                raise ValueError("header not available in closed files")

    # Compatibility functions for pysam < 0.14
    property text:
        """deprecated, use :attr:`references` and :attr:`lengths` instead"""
        def __get__(self):
            if self.header:
                return self.header.__str__()
            else:
                raise ValueError("header not available in closed files")

    # Compatibility functions for pysam < 0.8.3
    def gettid(self, reference):
        """deprecated, use :meth:`get_tid` instead"""
        return self.get_tid(reference)

    def getrname(self, tid):
        """deprecated, use :meth:`get_reference_name` instead"""
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
        cdef char *cindexname = NULL

        if not samfile.is_open:
            raise ValueError("I/O operation on closed file")

        # makes sure that samfile stays alive as long as the
        # iterator is alive
        self.samfile = samfile

        # reopen the file - note that this makes the iterator
        # slow and causes pileup to slow down significantly.
        if multiple_iterators:

            cfilename = samfile.filename
            with nogil:
                self.htsfile = hts_open(cfilename, 'r')
            assert self.htsfile != NULL

            if samfile.has_index():
                if samfile.index_filename:
                    cindexname = bindex_filename = encode_filename(samfile.index_filename)
                with nogil:
                    self.index = sam_index_load2(self.htsfile, cfilename, cindexname)
            else:
                self.index = NULL

            # need to advance in newly opened file to position after header
            # better: use seek/tell?
            with nogil:
                hdr = sam_hdr_read(self.htsfile)
            if hdr is NULL:
                raise IOError("unable to read header information")
            self.header = makeAlignmentHeader(hdr)

            self.owns_samfile = True

            # options specific to CRAM files
            if samfile.is_cram and samfile.reference_filename:
                creference_filename = samfile.reference_filename
                hts_set_opt(self.htsfile,
                            CRAM_OPT_REFERENCE,
                            creference_filename)

        else:
            self.htsfile = samfile.htsfile
            self.index = samfile.index
            self.owns_samfile = False
            self.header = samfile.header

        self.retval = 0

        self.b = bam_init1()

    def __dealloc__(self):
        bam_destroy1(self.b)
        if self.owns_samfile:
            hts_idx_destroy(self.index)
            hts_close(self.htsfile)


cdef class IteratorRowRegion(IteratorRow):
    """*(AlignmentFile samfile, int tid, int beg, int stop,
    int multiple_iterators=False)*

    iterate over mapped reads in a region.

    .. note::

        It is usually not necessary to create an object of this class
        explicitly. It is returned as a result of call to a
        :meth:`AlignmentFile.fetch`.

    """

    def __init__(self, AlignmentFile samfile,
                 int tid, int beg, int stop,
                 int multiple_iterators=False):

        if not samfile.has_index():
            raise ValueError("no index available for iteration")

        IteratorRow.__init__(self, samfile,
                             multiple_iterators=multiple_iterators)
        with nogil:
            self.iter = sam_itr_queryi(
                self.index,
                tid,
                beg,
                stop)

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
            return makeAlignedSegment(self.b, self.header)
        elif self.retval == -1:
            raise StopIteration
        elif self.retval == -2:
            # Note: it is currently not the case that hts_iter_next
            # returns -2 for a truncated file.
            # See https://github.com/pysam-developers/pysam/pull/50#issuecomment-64928625
            raise IOError('truncated file')
        else:
            raise IOError("error while reading file {}: {}".format(self.samfile.filename, self.retval))

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

    def __init__(self,
                 AlignmentFile samfile,
                 int n,
                 int multiple_iterators=False):

        IteratorRow.__init__(self, samfile,
                             multiple_iterators=multiple_iterators)

        self.max_rows = n
        self.current_row = 0

    def __iter__(self):
        return self

    cdef bam1_t * getCurrent(self):
        return self.b

    cdef int cnext(self):
        '''cversion of iterator. Used by IteratorColumn'''
        cdef int ret
        cdef bam_hdr_t * hdr = self.header.ptr
        with nogil:
            ret = sam_read1(self.htsfile,
                            hdr,
                            self.b)
        return ret

    def __next__(self):
        if self.current_row >= self.max_rows:
            raise StopIteration

        cdef int ret = self.cnext()
        if ret >= 0:
            self.current_row += 1
            return makeAlignedSegment(self.b, self.header)
        elif ret == -1:
            raise StopIteration
        else:
            raise IOError(read_failure_reason(ret))


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

    cdef bam1_t * getCurrent(self):
        return self.b

    cdef int cnext(self):
        '''cversion of iterator. Used by IteratorColumn'''
        cdef int ret
        cdef bam_hdr_t * hdr = self.header.ptr
        with nogil:
            ret = sam_read1(self.htsfile,
                            hdr,
                            self.b)
        return ret

    def __next__(self):
        cdef int ret = self.cnext()
        if ret >= 0:
            return makeAlignedSegment(self.b, self.header)
        elif ret == -1:
            raise StopIteration
        else:
            raise IOError(read_failure_reason(ret))


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
                                         MAX_POS)
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
                return makeAlignedSegment(self.rowiter.b, self.header)

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
        cdef bam_hdr_t * hdr = self.header.ptr
        with nogil:
            ret = sam_read1(self.htsfile,
                            hdr,
                            self.b)
        return ret

    def __next__(self):
        cdef int ret = self.cnext()
        if ret >= 0:
            return makeAlignedSegment(self.b, self.header)
        elif ret == -1:
            raise StopIteration
        else:
            raise IOError(read_failure_reason(ret))


cdef int __advance_nofilter(void *data, bam1_t *b):
    '''advance without any read filtering.
    '''
    cdef __iterdata * d = <__iterdata*>data
    cdef int ret
    with nogil:
        ret = sam_itr_next(d.htsfile, d.iter, b)
    return ret


cdef int __advance_raw_nofilter(void *data, bam1_t *b):
    '''advance (without iterator) without any read filtering.
    '''
    cdef __iterdata * d = <__iterdata*>data
    cdef int ret
    with nogil:
        ret = sam_read1(d.htsfile, d.header, b)
    return ret


cdef int __advance_all(void *data, bam1_t *b):
    '''only use reads for pileup passing basic filters such as

    BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
    '''

    cdef __iterdata * d = <__iterdata*>data
    cdef mask = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP
    cdef int ret
    while 1:
        with nogil:
            ret = sam_itr_next(d.htsfile, d.iter, b)
        if ret < 0:
            break
        if b.core.flag & d.flag_filter:
            continue
        break
    return ret


cdef int __advance_raw_all(void *data, bam1_t *b):
    '''only use reads for pileup passing basic filters such as

    BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
    '''

    cdef __iterdata * d = <__iterdata*>data
    cdef int ret
    while 1:
        with nogil:
            ret = sam_read1(d.htsfile, d.header, b)
        if ret < 0:
            break
        if b.core.flag & d.flag_filter:
            continue
        break
    return ret


cdef int __advance_samtools(void * data, bam1_t * b):
    '''advance using same filter and read processing as in
    the samtools pileup.
    '''
    cdef __iterdata * d = <__iterdata*>data
    cdef int ret
    cdef int q

    while 1:
        with nogil:
            ret = sam_itr_next(d.htsfile, d.iter, b) if d.iter else sam_read1(d.htsfile, d.header, b)
        if ret < 0:
            break
        if b.core.flag & d.flag_filter:
            continue
        if d.flag_require and not (b.core.flag & d.flag_require):
            continue

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
                    "reference sequence for '{}' (tid={}) not found".format(
                        d.header.target_name[d.tid], d.tid))

        # realign read - changes base qualities
        if d.seq != NULL and d.compute_baq:
            # 4th option to realign is flag:
            # apply_baq = flag&1, extend_baq = flag&2, redo_baq = flag&4
            if d.redo_baq:
                sam_prob_realn(b, d.seq, d.seq_len, 7)
            else:
                sam_prob_realn(b, d.seq, d.seq_len, 3)

        if d.seq != NULL and d.adjust_capq_threshold > 10:
            q = sam_cap_mapq(b, d.seq, d.seq_len, d.adjust_capq_threshold)
            if q < 0:
                continue
            elif b.core.qual > q:
                b.core.qual = q

        if b.core.qual < d.min_mapping_quality:
            continue
        if d.ignore_orphans and b.core.flag & BAM_FPAIRED and not (b.core.flag & BAM_FPROPER_PAIR):
            continue

        break

    return ret


cdef class IteratorColumn:
    '''abstract base class for iterators over columns.

    IteratorColumn objects wrap the pileup functionality of samtools.

    For reasons of efficiency, the iterator points to the current
    pileup buffer. The pileup buffer is updated at every iteration.
    This might cause some unexpected behaviour. For example,
    consider the conversion to a list::

       f = AlignmentFile("file.bam", "rb")
       result = list(f.pileup())

    Here, ``result`` will contain ``n`` objects of type
    :class:`~pysam.PileupColumn` for ``n`` columns, but each object in
    ``result`` will contain the same information.

    The desired behaviour can be achieved by list comprehension::

       result = [x.pileups() for x in f.pileup()]

    ``result`` will be a list of ``n`` lists of objects of type
    :class:`~pysam.PileupRead`.

    If the iterator is associated with a :class:`~pysam.Fastafile`
    using the :meth:`add_reference` method, then the iterator will
    export the current sequence via the methods :meth:`get_sequence`
    and :meth:`seq_len`.

    See :class:`~AlignmentFile.pileup` for kwargs to the iterator.
    '''

    def __cinit__( self, AlignmentFile samfile, **kwargs):
        self.samfile = samfile
        self.fastafile = kwargs.get("fastafile", None)
        self.stepper = kwargs.get("stepper", "samtools")
        self.max_depth = kwargs.get("max_depth", 8000)
        self.ignore_overlaps = kwargs.get("ignore_overlaps", True)
        self.min_base_quality = kwargs.get("min_base_quality", 13)
        self.iterdata.seq = NULL
        self.iterdata.min_mapping_quality = kwargs.get("min_mapping_quality", 0)
        self.iterdata.flag_require = kwargs.get("flag_require", 0)
        self.iterdata.flag_filter = kwargs.get("flag_filter", BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)
        self.iterdata.adjust_capq_threshold = kwargs.get("adjust_capq_threshold", 0)
        self.iterdata.compute_baq = kwargs.get("compute_baq", True)
        self.iterdata.redo_baq = kwargs.get("redo_baq", False)
        self.iterdata.ignore_orphans = kwargs.get("ignore_orphans", True)

        self.tid = 0
        self.pos = 0
        self.n_plp = 0
        self.plp = NULL
        self.pileup_iter = <bam_mplp_t>NULL

    def __iter__(self):
        return self

    cdef int cnext(self):
        '''perform next iteration.
        '''
        # do not release gil here because of call-backs
        cdef int ret = bam_mplp_auto(self.pileup_iter,
                                     &self.tid,
                                     &self.pos,
                                     &self.n_plp,
                                     &self.plp)
        return ret

    cdef char * get_sequence(self):
        '''return current reference sequence underlying the iterator.
        '''
        return self.iterdata.seq

    property seq_len:
        '''current sequence length.'''
        def __get__(self):
            return self.iterdata.seq_len

    def add_reference(self, FastaFile fastafile):
       '''
       add reference sequences in `fastafile` to iterator.'''
       self.fastafile = fastafile
       if self.iterdata.seq != NULL:
           free(self.iterdata.seq)
       self.iterdata.tid = -1
       self.iterdata.fastafile = self.fastafile.fastafile

    def has_reference(self):
        '''
        return true if iterator is associated with a reference'''
        return self.fastafile

    cdef _setup_iterator(self,
                         int tid,
                         int start,
                         int stop,
                         int multiple_iterators=0):
        '''setup the iterator structure'''

        self.iter = IteratorRowRegion(self.samfile, tid, start, stop, multiple_iterators)
        self.iterdata.htsfile = self.samfile.htsfile
        self.iterdata.iter = self.iter.iter
        self.iterdata.seq = NULL
        self.iterdata.tid = -1
        self.iterdata.header = self.samfile.header.ptr

        if self.fastafile is not None:
            self.iterdata.fastafile = self.fastafile.fastafile
        else:
            self.iterdata.fastafile = NULL

        # Free any previously allocated memory before reassigning
        # pileup_iter
        self._free_pileup_iter()

        cdef void * data[1]
        data[0] = <void*>&self.iterdata

        if self.stepper is None or self.stepper == "all":
            with nogil:
                self.pileup_iter = bam_mplp_init(1,
                                                 <bam_plp_auto_f>&__advance_all,
                                                 data)
        elif self.stepper == "nofilter":
            with nogil:
                self.pileup_iter = bam_mplp_init(1,
                                                 <bam_plp_auto_f>&__advance_nofilter,
                                                 data)
        elif self.stepper == "samtools":
            with nogil:
                self.pileup_iter = bam_mplp_init(1,
                                                 <bam_plp_auto_f>&__advance_samtools,
                                                 data)
        else:
            raise ValueError(
                "unknown stepper option `%s` in IteratorColumn" % self.stepper)

        if self.max_depth:
            with nogil:
                bam_mplp_set_maxcnt(self.pileup_iter, self.max_depth)

        if self.ignore_overlaps:
            with nogil:
                bam_mplp_init_overlaps(self.pileup_iter)

    cdef _setup_raw_rest_iterator(self):
        '''set up an "iterator" that just uses sam_read1(), similar to HTS_IDX_REST'''

        self.iter = None
        self.iterdata.iter = NULL
        self.iterdata.htsfile = self.samfile.htsfile
        self.iterdata.seq = NULL
        self.iterdata.tid = -1
        self.iterdata.header = self.samfile.header.ptr

        if self.fastafile is not None:
            self.iterdata.fastafile = self.fastafile.fastafile
        else:
            self.iterdata.fastafile = NULL

        # Free any previously allocated memory before reassigning
        # pileup_iter
        self._free_pileup_iter()

        cdef void * data[1]
        data[0] = <void*>&self.iterdata

        if self.stepper is None or self.stepper == "all":
            with nogil:
                self.pileup_iter = bam_mplp_init(1,
                                                 <bam_plp_auto_f>&__advance_raw_all,
                                                 data)
        elif self.stepper == "nofilter":
            with nogil:
                self.pileup_iter = bam_mplp_init(1,
                                                 <bam_plp_auto_f>&__advance_raw_nofilter,
                                                 data)
        elif self.stepper == "samtools":
            with nogil:
                self.pileup_iter = bam_mplp_init(1,
                                                 <bam_plp_auto_f>&__advance_samtools,
                                                 data)
        else:
            raise ValueError(
                "unknown stepper option `%s` in IteratorColumn" % self.stepper)

        if self.max_depth:
            with nogil:
                bam_mplp_set_maxcnt(self.pileup_iter, self.max_depth)

        if self.ignore_overlaps:
            with nogil:
                bam_mplp_init_overlaps(self.pileup_iter)

    cdef reset(self, tid, start, stop):
        '''reset iterator position.

        This permits using the iterator multiple times without
        having to incur the full set-up costs.
        '''
        if self.iter is None:
            raise TypeError("Raw iterator set up without region cannot be reset")

        self.iter = IteratorRowRegion(self.samfile, tid, start, stop, multiple_iterators=0)
        self.iterdata.iter = self.iter.iter

        # invalidate sequence if different tid
        if self.tid != tid:
            if self.iterdata.seq != NULL:
                free(self.iterdata.seq)
            self.iterdata.seq = NULL
            self.iterdata.tid = -1

        # self.pileup_iter = bam_mplp_init(1
        #                                  &__advancepileup,
        #                                  &self.iterdata)
        with nogil:
            bam_mplp_reset(self.pileup_iter)

    cdef _free_pileup_iter(self):
        '''free the memory alloc'd by bam_plp_init.

        This is needed before setup_iterator allocates another
        pileup_iter, or else memory will be lost.  '''
        if self.pileup_iter != <bam_mplp_t>NULL:
            with nogil:
                bam_mplp_reset(self.pileup_iter)
                bam_mplp_destroy(self.pileup_iter)
                self.pileup_iter = <bam_mplp_t>NULL

    def __dealloc__(self):
        # reset in order to avoid memory leak messages for iterators
        # that have not been fully consumed
        self._free_pileup_iter()
        self.plp = <const bam_pileup1_t*>NULL

        if self.iterdata.seq != NULL:
            free(self.iterdata.seq)
            self.iterdata.seq = NULL

    # backwards compatibility

    def hasReference(self):
        return self.has_reference()
    cdef char * getSequence(self):
        return self.get_sequence()
    def addReference(self, FastaFile fastafile):
        return self.add_reference(fastafile)


cdef class IteratorColumnRegion(IteratorColumn):
    '''iterates over a region only.
    '''
    def __cinit__(self,
                  AlignmentFile samfile,
                  int tid = 0,
                  int start = 0,
                  int stop = MAX_POS,
                  int truncate = False,
                  int multiple_iterators = True,
                  **kwargs ):

        # initialize iterator. Multiple iterators not available
        # for CRAM.
        if multiple_iterators and samfile.is_cram:
            warnings.warn("multiple_iterators not implemented for CRAM")
            multiple_iterators = False

        self._setup_iterator(tid, start, stop, multiple_iterators)
        self.start = start
        self.stop = stop
        self.truncate = truncate

    def __next__(self):

        cdef int n

        while 1:
            n = self.cnext()
            if n < 0:
                raise ValueError("error during iteration" )

            if n == 0:
                raise StopIteration

            if self.truncate:
                if self.start > self.pos:
                    continue
                if self.pos >= self.stop:
                    raise StopIteration

            return makePileupColumn(&self.plp,
                                    self.tid,
                                    self.pos,
                                    self.n_plp,
                                    self.min_base_quality,
                                    self.iterdata.seq,
                                    self.samfile.header)


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
        self._setup_iterator(self.tid, 0, MAX_POS, 1)

    def __next__(self):

        cdef int n
        while 1:
            n = self.cnext()
            if n < 0:
                raise ValueError("error during iteration")

            # proceed to next reference or stop
            if n == 0:
                self.tid += 1
                if self.tid < self.samfile.nreferences:
                    self._setup_iterator(self.tid, 0, MAX_POS, 0)
                else:
                    raise StopIteration
                continue

            # return result, if within same reference
            return makePileupColumn(&self.plp,
                                    self.tid,
                                    self.pos,
                                    self.n_plp,
                                    self.min_base_quality,
                                    self.iterdata.seq,
                                    self.samfile.header)


cdef class IteratorColumnAll(IteratorColumn):
    """iterates over all columns, without using an index
    """

    def __cinit__(self,
                  AlignmentFile samfile,
                  **kwargs):

        self._setup_raw_rest_iterator()

    def __next__(self):

        cdef int n
        n = self.cnext()
        if n < 0:
            raise ValueError("error during iteration")

        if n == 0:
            raise StopIteration

        return makePileupColumn(&self.plp,
                                self.tid,
                                self.pos,
                                self.n_plp,
                                self.min_base_quality,
                                self.iterdata.seq,
                                self.samfile.header)


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
    """Index a Sam/BAM-file by query name while keeping the
    original sort order intact.

    The index is kept in memory and can be substantial.

    By default, the file is re-opened to avoid conflicts if multiple
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
        cdef bam_hdr_t * hdr = NULL
        assert samfile.is_bam, "can only apply IndexReads on bam files"

        # multiple_iterators the file - note that this makes the iterator
        # slow and causes pileup to slow down significantly.
        if multiple_iterators:
            cfilename = samfile.filename
            with nogil:
                self.htsfile = hts_open(cfilename, 'r')
            if self.htsfile == NULL:
                raise OSError("unable to reopen htsfile")

            # need to advance in newly opened file to position after header
            # better: use seek/tell?
            with nogil:
                hdr = sam_hdr_read(self.htsfile)
            if hdr == NULL:
                raise OSError("unable to read header information")
            self.header = makeAlignmentHeader(hdr)
            self.owns_samfile = True
        else:
            self.htsfile = self.samfile.htsfile
            self.header = samfile.header
            self.owns_samfile = False

    def build(self):
        '''build the index.'''

        self.index = collections.defaultdict(list)

        # this method will start indexing from the current file position
        cdef int ret = 1
        cdef bam1_t * b = <bam1_t*>calloc(1, sizeof( bam1_t))
        if b == NULL:
            raise MemoryError("could not allocate {} bytes".format(sizeof(bam1_t)))

        cdef uint64_t pos
        cdef bam_hdr_t * hdr = self.header.ptr

        while ret > 0:
            with nogil:
                pos = bgzf_tell(hts_get_bgzfp(self.htsfile))
                ret = sam_read1(self.htsfile,
                                hdr,
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
