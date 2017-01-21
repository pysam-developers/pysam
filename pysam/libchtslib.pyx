# cython: embedsignature=True
# cython: profile=True
# adds doc-strings for sphinx
import os

from posix.unistd cimport dup

from pysam.libchtslib cimport *

from pysam.libcutils cimport force_bytes, force_str, charptr_to_str, charptr_to_str_w_len
from pysam.libcutils cimport encode_filename, from_string_and_size


__all__ = ["get_verbosity", "set_verbosity"]


########################################################################
########################################################################
## Constants
########################################################################

cdef int   MAX_POS = 2 << 29
cdef tuple FORMAT_CATEGORIES = ('UNKNOWN', 'ALIGNMENTS', 'VARIANTS', 'INDEX', 'REGIONS')
cdef tuple FORMATS = ('UNKNOWN', 'BINARY_FORMAT', 'TEXT_FORMAT', 'SAM', 'BAM', 'BAI', 'CRAM', 'CRAI',
                      'VCF', 'BCF', 'CSI', 'GZI', 'TBI', 'BED')
cdef tuple COMPRESSION = ('NONE', 'GZIP', 'BGZF', 'CUSTOM')


cpdef set_verbosity(int verbosity):
    """Set htslib's hts_verbose global variable to the specified value."""
    return hts_set_verbosity(verbosity)

cpdef get_verbosity():
    """Return the value of htslib's hts_verbose global variable."""
    return hts_get_verbosity()


class CallableValue(object):
    def __init__(self, value):
        self.value = value
    def __call__(self):
        return self.value
    def __bool__(self):
        return self.value
    def __nonzero__(self):
        return self.value
    def __eq__(self, other):
        return self.value == other
    def __ne__(self, other):
        return self.value != other


CTrue = CallableValue(True)
CFalse = CallableValue(False)


cdef class HTSFile(object):
    """
    Base class for HTS file types
    """
    def __cinit__(self, *args, **kwargs):
        self.htsfile = NULL
        self.duplicate_filehandle = True

    def __dealloc__(self):
        if self.htsfile:
            hts_close(self.htsfile)
            self.htsfile = NULL

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    @property
    def category(self):
        """General file format category.  One of UNKNOWN, ALIGNMENTS,
        VARIANTS, INDEX, REGIONS"""
        if not self.htsfile:
            raise ValueError('metadata not available on closed file')
        return FORMAT_CATEGORIES[self.htsfile.format.category]

    @property
    def format(self):
        """File format.

        One of UNKNOWN, BINARY_FORMAT, TEXT_FORMAT, SAM, BAM,
        BAI, CRAM, CRAI, VCF, BCF, CSI, GZI, TBI, BED.
        """
        if not self.htsfile:
            raise ValueError('metadata not available on closed file')
        return FORMATS[self.htsfile.format.format]

    @property
    def version(self):
        """Tuple of file format version numbers (major, minor)"""
        if not self.htsfile:
            raise ValueError('metadata not available on closed file')
        return self.htsfile.format.version.major, self.htsfile.format.version.minor

    @property
    def compression(self):
        """File compression.

        One of NONE, GZIP, BGZF, CUSTOM."""
        if not self.htsfile:
            raise ValueError('metadata not available on closed file')
        return COMPRESSION[self.htsfile.format.compression]

    @property
    def description(self):
        """Vaguely human readable description of the file format"""
        if not self.htsfile:
            raise ValueError('metadata not available on closed file')
        cdef char *desc = hts_format_description(&self.htsfile.format)
        try:
            return charptr_to_str(desc)
        finally:
            free(desc)

    @property
    def is_open(self):
        """return True if HTSFile is open and in a valid state."""
        return CTrue if self.htsfile != NULL else CFalse

    @property
    def is_closed(self):
        """return True if HTSFile is closed."""
        return self.htsfile == NULL

    @property
    def closed(self):
        """return True if HTSFile is closed."""
        return self.htsfile == NULL

    @property
    def is_write(self):
        """return True if HTSFile is open for writing"""
        return self.htsfile != NULL and self.htsfile.is_write != 0

    @property
    def is_read(self):
        """return True if HTSFile is open for reading"""
        return self.htsfile != NULL and self.htsfile.is_write == 0

    @property
    def is_sam(self):
        """return True if HTSFile is reading or writing a SAM alignment file"""
        return self.htsfile != NULL and self.htsfile.format.format == sam

    @property
    def is_bam(self):
        """return True if HTSFile is reading or writing a BAM alignment file"""
        return self.htsfile != NULL and self.htsfile.format.format == bam

    @property
    def is_cram(self):
        """return True if HTSFile is reading or writing a BAM alignment file"""
        return self.htsfile != NULL and self.htsfile.format.format == cram

    @property
    def is_vcf(self):
        """return True if HTSFile is reading or writing a VCF variant file"""
        return self.htsfile != NULL and self.htsfile.format.format == vcf

    @property
    def is_bcf(self):
        """return True if HTSFile is reading or writing a BCF variant file"""
        return self.htsfile != NULL and self.htsfile.format.format == bcf

    def reset(self):
        """reset file position to beginning of file just after the header.

        Returns
        -------

        The file position after moving the file pointer.

        """
        return self.seek(self.start_offset)

    def seek(self, uint64_t offset):
        """move file pointer to position *offset*, see :meth:`pysam.HTSFile.tell`."""
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
        """return current file position, see :meth:`pysam.HTSFile.seek`."""
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

    cdef htsFile *_open_htsfile(self) except? NULL:
        cdef char *cfilename
        cdef char *cmode = self.mode
        cdef int fd, dup_fd

        if isinstance(self.filename, bytes):
            cfilename = self.filename
            with nogil:
                return hts_open(cfilename, cmode)
        else:
            if isinstance(self.filename, int):
                fd = self.filename
            else:
                fd = self.filename.fileno()
               
            if self.duplicate_filehandle:
                dup_fd = dup(fd)
            else:
                dup_fd = fd

            # Replicate mode normalization done in hts_open_format
            smode = self.mode.replace(b'b', b'').replace(b'c', b'')
            if b'b' in self.mode:
                smode += b'b'
            elif b'c' in self.mode:
                smode += b'c'
            cmode = smode

            hfile = hdopen(dup_fd, cmode)
            if hfile == NULL:
                raise IOError('Cannot create hfile')

            try:
                # filename.name can be an int
                filename = str(self.filename.name)
            except AttributeError:
                filename = '<fd:{}>'.format(fd)

            filename = encode_filename(filename)
            cfilename = filename
            with nogil:
                return hts_hopen(hfile, cfilename, cmode)

    def _exists(self):
        """return False iff file is local, a file and exists.
        """
        return (not isinstance(self.filename, (str, bytes)) or
                self.filename == b'-' or
                self.is_remote or
                os.path.exists(self.filename))
