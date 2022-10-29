# cython: language_level=3
"""Functions that read and write block gzipped files.

The user of the file doesn't have to worry about the compression
and random access is allowed if an index file is present."""

# based on Python 3.5's gzip module

import io

from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free

from cpython.object cimport PyObject
from cpython.bytes  cimport PyBytes_FromStringAndSize, _PyBytes_Resize

from pysam.libcutils   cimport force_bytes, encode_filename
from pysam.libchtslib  cimport bgzf_open, bgzf_index_build_init, bgzf_write, bgzf_read, \
                               bgzf_flush, bgzf_index_dump, bgzf_close, bgzf_seek, \
                               bgzf_tell, bgzf_getline, kstring_t, SEEK_SET, BGZF

__all__ = ["BGZFile"]


BUFFER_SIZE = io.DEFAULT_BUFFER_SIZE


cdef class BGZFile(object):
    """The BGZFile class simulates most of the methods of a file object with
    the exception of the truncate() method.

    This class only supports opening files in binary mode. If you need to open a
    compressed file in text mode, use the gzip.open() function.
    """
    cdef BGZF* bgzf
    cdef readonly object name, index

    def __init__(self, filename, mode=None, index=None):
        """Constructor for the BGZFile class.

        The mode argument can be any of 'r', 'rb', 'a', 'ab', 'w', 'wb', 'x', or
        'xb' depending on whether the file will be read or written.  The default
        is the mode of fileobj if discernible; otherwise, the default is 'rb'.
        A mode of 'r' is equivalent to one of 'rb', and similarly for 'w' and
        'wb', 'a' and 'ab', and 'x' and 'xb'.
        """
        if mode and ('t' in mode or 'U' in mode):
            raise ValueError("Invalid mode: {!r}".format(mode))
        if not mode:
            mode = 'rb'
        elif mode and 'b' not in mode:
            mode += 'b'

        mode = force_bytes(mode)

        self.name = encode_filename(filename)
        self.index = encode_filename(index) if index is not None else None

        self.bgzf = bgzf_open(self.name, mode)

        if self.bgzf.is_write and index is not None and bgzf_index_build_init(self.bgzf) < 0:
            raise IOError('Error building bgzf index')

    def __dealloc__(self):
        self.close()

    def write(self, data):
        if not self.bgzf:
            raise ValueError("write() on closed BGZFile object")

        if not self.bgzf.is_write:
            import errno
            raise IOError(errno.EBADF, "write() on read-only BGZFile object")

        if isinstance(data, bytes):
            length = len(data)
        else:
            # accept any data that supports the buffer protocol
            data = memoryview(data)
            length = data.nbytes

        if length > 0 and bgzf_write(self.bgzf, <char *>data, length) < 0:
            raise IOError('BGZFile write failed')

        return length

    def read(self, size=-1):
        cdef ssize_t read_size

        if not self.bgzf:
            raise ValueError("read() on closed BGZFile object")

        if self.bgzf.is_write:
            import errno
            raise IOError(errno.EBADF, "read() on write-only BGZFile object")

        if size < 0:
            chunks = []
            while 1:
                chunk = PyBytes_FromStringAndSize(NULL, BUFFER_SIZE)
                cdata = <bytes>chunk
                read_size = bgzf_read(self.bgzf, <char *>chunk, BUFFER_SIZE)
                if read_size < 0:
                    raise IOError('Error reading from BGZFile')
                elif not read_size:
                    break
                elif read_size < BUFFER_SIZE:
                    chunk = chunk[:read_size]
                chunks.append(chunk)
            return b''.join(chunks)

        elif size > 0:
            chunk = PyBytes_FromStringAndSize(NULL, size)
            read_size = bgzf_read(self.bgzf, <char *>chunk, size)
            if read_size < 0:
                raise IOError('Error reading from BGZFile')
            elif read_size < size:
                chunk = chunk[:read_size]
            return chunk
        else:
            return b''

    @property
    def closed(self):
        return self.bgzf == NULL

    def close(self):
        if not self.bgzf:
            return

        if self.bgzf.is_write and bgzf_flush(self.bgzf) < 0:
            raise IOError('Error flushing BGZFile object')

        if self.index and bgzf_index_dump(self.bgzf, self.index, NULL) < 0:
            raise IOError('Cannot write index')

        cdef ret = bgzf_close(self.bgzf)
        self.bgzf = NULL

        if ret < 0:
            raise IOError('Error closing BGZFile object')

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        self.close()

    def flush(self):
        if not self.bgzf:
            return

        if self.bgzf.is_write and bgzf_flush(self.bgzf) < 0:
            raise IOError('Error flushing BGZFile object')

    def fileno(self):
        """Invoke the underlying file object's fileno() method.

        This will raise AttributeError if the underlying file object
        doesn't support fileno().
        """
        raise AttributeError('fileno')

    def rewind(self):
        '''Return the uncompressed stream file position indicator to the
        beginning of the file'''
        if not self.bgzf:
            raise ValueError("rewind() on closed BGZFile object")
        if not self.bgzf.is_write:
            raise IOError("Can't rewind in write mode")
        if bgzf_seek(self.bgzf, 0, SEEK_SET) < 0:
            raise IOError('Error seeking BGZFFile object')

    def readable(self):
        if not self.bgzf:
            raise ValueError("readable() on closed BGZFile object")
        return self.bgzf != NULL and not self.bgzf.is_write

    def writable(self):
        return self.bgzf != NULL and self.bgzf.is_write

    def seekable(self):
        return True

    def tell(self):
        if not self.bgzf:
            raise ValueError("seek() on closed BGZFile object")
        cdef int64_t off = bgzf_tell(self.bgzf)
        if off < 0:
            raise IOError('Error in tell on BGZFFile object')

        return off

    def seek(self, offset, whence=io.SEEK_SET):
        if not self.bgzf:
            raise ValueError("seek() on closed BGZFile object")
        if whence is not io.SEEK_SET:
            raise ValueError('Seek from end not supported')

        cdef int64_t off = bgzf_seek(self.bgzf, offset, SEEK_SET)
        if off < 0:
            raise IOError('Error seeking BGZFFile object')

        return off

    def readline(self, size=-1):
        if not self.bgzf:
            raise ValueError("readline() on closed BGZFile object")

        cdef kstring_t line
        cdef char c

        line.l = line.m = 0
        line.s = NULL

        cdef int ret = bgzf_getline(self.bgzf, b'\n', &line)
        if ret == -1:
            s = b''
        elif ret == -2:
            if line.m:
                free(line.s)
            raise IOError('Error reading line in BGZFFile object')
        else:
            s = line.s[:line.l]

        if line.m:
            free(line.s)

        return s

    def __iter__(self):
        return self

    def __next__(self):
        line = self.readline()
        if not line:
            raise StopIteration()
        return line
