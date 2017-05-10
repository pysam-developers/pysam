#cdef extern from "Python.h":
#    ctypedef struct FILE

from libc.stdint cimport uint8_t, int32_t, uint32_t, int64_t, uint64_t

cdef class TupleProxy:

    cdef:
        char * data
        char ** fields
        int nfields
        int index
        int nbytes
        int offset
        bint is_modified

    cdef encoding

    cpdef int getMaxFields(self)
    cpdef int getMinFields(self)
#    cdef char * _getindex(self, int idx)

    cdef take(self, char * buffer, size_t nbytes)
    cdef present(self, char * buffer, size_t nbytes)
    cdef copy(self, char * buffer, size_t nbytes, bint reset=*)
    cdef update(self, char * buffer, size_t nbytes)


cdef class NamedTupleProxy(TupleProxy):
    pass


cdef class GTFProxy(NamedTupleProxy):
    cdef object attribute_dict
    cpdef int getMaxFields(self)
    cpdef int getMinFields(self)


cdef class GFF3Proxy(GTFProxy):
    pass


cdef class BedProxy(NamedTupleProxy):

    cdef:
        char * contig
        uint32_t start
        uint32_t end
        int bedfields

    cpdef int getMaxFields(self)
    cpdef int getMinFields(self)
    cdef update(self, char * buffer, size_t nbytes)

cdef class VCFProxy(NamedTupleProxy) :

    cdef:
        char * contig
        uint32_t pos

    cdef update(self, char * buffer, size_t nbytes)
