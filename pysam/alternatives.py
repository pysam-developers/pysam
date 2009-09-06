# This file contains an alternative implementation
# to call samtools functions. These are direct calls.
# Plus: less overhead
# Minus: more trouble in maintaining

#in csamtools.pxd
    # samtools toolkit functions
    ctypedef int (*pysam_samtools_f)(int argc, char *argv[])

    int bam_taf2baf(int argc, char *argv[])
    int bam_pileup(int argc, char *argv[])
    int bam_merge(int argc, char *argv[])
    int bam_index(int argc, char *argv[])
    int bam_sort(int argc, char *argv[])
    int bam_tview_main(int argc, char *argv[])
    int bam_mating(int argc, char *argv[])
    int bam_rmdup(int argc, char *argv[])
    int bam_rmdupse(int argc, char *argv[])
    int bam_flagstat(int argc, char *argv[])
    int bam_fillmd(int argc, char *argv[])
    int main_samview(int argc, char *argv[])
    int main_import(int argc, char *argv[])
    int faidx_main(int argc, char *argv[])
    int glf3_view_main(int argc, char *argv[])


## Alternative code in csamtools.pyx
cdef class SamtoolsWrapper:
    '''generic wrapper around samtools functions'''
    cdef pysam_samtools_f f 

    def __init__(self): self.f = NULL

    def call(self, *args ):

        if self.f == NULL: raise NotImplementedError("invalid call to base class" )

        cdef char ** cargs
        cdef int i, n, retval
        n = len(args)
        # allocate one more for first (dummy) argument (contains command)
        cargs = <char**>calloc( n+1, sizeof( char *) )
        cargs[0] = "method"
        for i from 0 <= i < n:
            cargs[i+1] = args[i]
        for i from 0 <= i < n+1:
            print cargs[i]
        retval = self.f(n+1, cargs)
        free( cargs )
        return retval

cdef class SamtoolsWrapperImport( SamtoolsWrapper ):
    def __init__(self): self.f = main_import
cdef class SamtoolsWrapperPileup( SamtoolsWrapper ):
    def __init__(self): self.f = bam_pileup
cdef class SamtoolsWrapperMerge( SamtoolsWrapper ):
    def __init__(self): self.f = bam_merge
cdef class SamtoolsWrapperSort( SamtoolsWrapper ):
    def __init__(self): self.f = bam_sort
cdef class SamtoolsWrapperIndex( SamtoolsWrapper ):
    def __init__(self): self.f = bam_index
cdef class SamtoolsWrapperFaidx( SamtoolsWrapper ):
    def __init__(self): self.f = faidx_main
cdef class SamtoolsWrapperFixMate( SamtoolsWrapper ):
    def __init__(self): self.f = bam_mating
cdef class SamtoolsWrapperRmDup( SamtoolsWrapper ):
    def __init__(self): self.f = bam_rmdup
cdef class SamtoolsWrapperRmDupSe( SamtoolsWrapper ):
    def __init__(self): self.f = bam_rmdupse
cdef class SamtoolsWrapperGlf3Viwen( SamtoolsWrapper ):
    def __init__(self): self.f = glf3_view_main
cdef class SamtoolsWrapperFlagStat( SamtoolsWrapper ):
    def __init__(self): self.f = bam_flagstat
cdef class SamtoolsWrapperFillMd( SamtoolsWrapper ):
    def __init__(self): self.f = bam_fillmd
cdef class SamtoolsWrapperCalMd( SamtoolsWrapper ):
    def __init__(self): self.f = bam_fillmd

automatic creation of these functions does not work 
due to pyrex/cython

def sort( *args, **kwargs ): return SamtoolsWrapperSort().call(*args, **kwargs)
