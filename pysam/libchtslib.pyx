# cython: embedsignature=True
# cython: profile=True
# adds doc-strings for sphinx
import os
import io

from posix.unistd cimport dup
from libc.errno  cimport errno
from cpython cimport PyBytes_FromStringAndSize

from pysam.libchtslib cimport *

from pysam.libcutils cimport force_bytes, force_str, charptr_to_str, charptr_to_str_w_len
from pysam.libcutils cimport encode_filename, from_string_and_size


from warnings         import warn


__all__ = ['get_verbosity', 'set_verbosity', 'HFile', 'HTSFile']

# defines imported from samtools
DEF SEEK_SET = 0
DEF SEEK_CUR = 1
DEF SEEK_END = 2

########################################################################
########################################################################
## Constants
########################################################################

# maximum genomic coordinace
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


cdef class HFile(object):
    cdef hFILE *fp
    cdef readonly object name, mode

    def __init__(self, name, mode='r', closedf=True):
        self._open(name, mode, closefd=True)

    def __dealloc__(self):
        self.close()

    @property
    def closed(self):
        return self.fp == NULL

    cdef _open(self, name, mode, closefd=True):
        self.name = name
        self.mode = mode

        mode = force_bytes(mode)

        if isinstance(name, int):
            if self.fp != NULL:
                name = dup(name)
            self.fp = hdopen(name, mode)
        else:
            name = encode_filename(name)
            self.fp = hopen(name, mode)

        if not self.fp:
            raise OSError(errno, 'failed to open HFile', self.name)

    def close(self):
        if self.fp == NULL:
            return

        cdef hFILE *fp = self.fp
        self.fp = NULL

        if hclose(fp) != 0:
            raise OSError(herrno(self.fp), 'failed to close HFile', self.name)

    def fileno(self):
        if self.fp == NULL:
            raise OSError('operation on closed HFile')
        if isinstance(self.name, int):
            return self.name
        else:
            raise AttributeError('fileno not available')

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        self.close()

    def __iter__(self):
        return self

    def __next__(self):
        line = self.readline()
        if not line:
            raise StopIteration()
        return line

    def flush(self):
        if self.fp == NULL:
            raise OSError('operation on closed HFile')
        if hflush(self.fp) != 0:
            raise OSError(herrno(self.fp), 'failed to flush HFile', self.name)

    def isatty(self):
        if self.fp == NULL:
            raise OSError('operation on closed HFile')
        return False

    def readable(self):
        return self.fp != NULL and 'r' in self.mode

    def read(self, Py_ssize_t size=-1):
        if self.fp == NULL:
            raise OSError('operation on closed HFile')

        if size == 0:
            return b''

        cdef list parts = []
        cdef bytes part
        cdef Py_ssize_t chunk_size, ret, bytes_read = 0
        cdef char *cpart

        while size == -1 or bytes_read < size:
            chunk_size = 4096
            if size != -1:
                chunk_size = min(chunk_size, size - bytes_read)

            part = PyBytes_FromStringAndSize(NULL, chunk_size)
            cpart = <char *>part
            ret = hread(self.fp, <void *>cpart, chunk_size)

            if ret < 0:
                OSError(herrno(self.fp), 'failed to read HFile', self.name)
            elif not ret:
                break

            bytes_read += ret

            if ret < chunk_size:
                part = cpart[:ret]

            parts.append(part)

        return b''.join(parts)

    def readall(self):
        return self.read()

    def readinto(self, buf):
        if self.fp == NULL:
            raise OSError('operation on closed HFile')

        size = len(buf)

        if size == 0:
            return size

        mv = memoryview(buf)
        ret = hread(self.fp, <void *>mv, size)

        if ret < 0:
            OSError(herrno(self.fp), 'failed to read HFile', self.name)

        return ret

    def readline(self, Py_ssize_t size=-1):
        if self.fp == NULL:
            raise OSError('operation on closed HFile')

        if size == 0:
            return b''

        cdef list parts = []
        cdef bytes part
        cdef Py_ssize_t chunk_size, ret, bytes_read = 0
        cdef char *cpart

        while size == -1 or bytes_read < size:
            chunk_size = 4096
            if size != -1:
                chunk_size = min(chunk_size, size - bytes_read)

            part = PyBytes_FromStringAndSize(NULL, chunk_size)
            cpart = <char *>part

            # Python bytes objects allocate an extra byte for a null terminator
            ret = hgetln(cpart, chunk_size+1, self.fp)

            if ret < 0:
                OSError(herrno(self.fp), 'failed to read HFile', self.name)
            elif not ret:
                break

            bytes_read += ret

            if ret < chunk_size:
                part = cpart[:ret]
                cpart = <char *>part

            parts.append(part)

            if cpart[ret-1] == b'\n':
               break

        return b''.join(parts)

    def readlines(self):
        return list(self)

    def seek(self, Py_ssize_t offset, int whence=SEEK_SET):
        if self.fp == NULL:
            raise OSError('operation on closed HFile')

        cdef Py_ssize_t off = hseek(self.fp, offset, whence)

        if off < 0:
            raise OSError(herrno(self.fp), 'seek failed on HFile', self.name)

        return off

    def tell(self):
        if self.fp == NULL:
            raise OSError('operation on closed HFile')

        ret = htell(self.fp)

        if ret < 0:
            raise OSError(herrno(self.fp), 'tell failed on HFile', self.name)

        return ret

    def seekable(self):
        return self.fp != NULL

    def truncate(self, size=None):
        raise NotImplementedError()

    def writable(self):
        return self.fp != NULL and 'w' in self.mode

    def write(self, bytes b):
        if self.fp == NULL:
            raise OSError('operation on closed HFile')

        got = hwrite(self.fp, <void *>b, len(b))

        if got < 0:
            raise OSError(herrno(self.fp), 'write failed on HFile', self.name)

        return got

    def writelines(self, lines):
        for line in lines:
            self.write(line)


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

    def close(self):
        if self.htsfile:
            hts_close(self.htsfile)
            self.htsfile = NULL

    def __dealloc__(self):
        if self.htsfile:
            hts_close(self.htsfile)
            self.htsfile = NULL

    def check_truncation(self, ignore_truncation=False):
        """Check if file is truncated."""
        if not self.htsfile:
            return

        if self.htsfile.format.compression != bgzf:
            return

        cdef BGZF *bgzfp = hts_get_bgzfp(self.htsfile)
        if not bgzfp:
            return

        cdef int ret = bgzf_check_EOF(bgzfp)
        if ret < 0:
            raise OSError(errno, 'error checking for EOF marker')
        elif ret == 0:
            msg = 'no BGZF EOF marker; file may be truncated'.format(self.filename)
            if ignore_truncation:
                warn(msg)
            else:
                raise OSError(msg)

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
        if self.htsfile.format.compression == bgzf:
            with nogil:
                ret = bgzf_seek(hts_get_bgzfp(self.htsfile), offset, SEEK_SET)
        elif self.htsfile.format.compression == no_compression:
            with nogil:
                ret = hts_useek(self.htsfile, <int>offset, SEEK_SET)
        else:
            raise NotImplementedError("seek not implemented in files compressed by method {}".format(
                self.htsfile.format.compression))
        return ret

    def tell(self):
        """return current file position, see :meth:`pysam.HTSFile.seek`."""
        if not self.is_open:
            raise ValueError('I/O operation on closed file')
        if self.is_stream:
            raise OSError('tell not available in streams')

        cdef int64_t ret
        if self.htsfile.format.compression == bgzf:
            with nogil:
                ret = bgzf_tell(hts_get_bgzfp(self.htsfile))
        elif self.htsfile.format.compression == no_compression:
            with nogil:
                ret = hts_utell(self.htsfile)
        elif self.htsfile.format.format == cram:
            with nogil:
                ret = htell(cram_fd_get_fp(self.htsfile.fp.cram))
        else:
            raise NotImplementedError("seek not implemented in files compressed by method {}".format(
                self.htsfile.format.compression))

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
