# cython: embedsignature=True
# adds doc-strings for sphinx

ctypedef struct tmpstruct_t:
    int beg
    int end
    samfile_t * input   

cdef int fetch_func( bam1_t *b, void *data):
   cdef bam_plbuf_t * buf
   buf = <bam_plbuf_t*>data
   bam_plbuf_push(b, buf)
   return 0

cdef int pileup_func( uint32_t tid, uint32_t pos, int n, bam_pileup1_t *pl, void *data):
   cdef tmpstruct_t * tmp 
   tmp = <tmpstruct_t*>data
   if pos >= tmp.beg and pos < tmp.end:
       print "%s\t%d\t%d" % ( tmp.input.header.target_name[tid], pos + 1, n )
   return 0

cdef samfile_t * openSam( filename, mode ):
    """open a samfile *filename*."""
    cdef samfile_t * fh
    
    fh = samopen( filename, mode, NULL)

    if fh == NULL:
        raise IOError("could not open file `%s`" % filename )

    return fh
  
cdef bam_index_t * openIndex( filename ):
    """open index for *filename*."""
    cdef bam_index_t *idx
    idx = bam_index_load(filename)

    if idx == NULL:
        raise IOError("could not open index `%s` " % filename )

    return idx

def test_pileup( filename ):
    """pileup test"""
    cdef tmpstruct_t tmp

    tmp.beg = 0
    tmp.end = 0x7fffffff
    tmp.input = openSam( filename, "rb")

    sampileup( tmp.input, -1, pileup_func, &tmp);  
    
    samclose( tmp.input )

def test_fetch( filename, region ):

    cdef tmpstruct_t tmp

    tmp.beg = 0
    tmp.end = 0x7fffffff
    tmp.input = openSam( filename, "rb")
    
    cdef int ref
    cdef bam_index_t *idx
    cdef bam_plbuf_t *buf

    # load bam index
    idx = openIndex( filename )

    # parse the region
    bam_parse_region(tmp.input.header, region, &ref, &tmp.beg, &tmp.end)
    if ref < 0:
        raise ValueError( "Invalid region %s" % region )

    # initialize pileup
    buf = bam_plbuf_init(pileup_func, &tmp)
    bam_fetch(tmp.input.x.bam, idx, ref, tmp.beg, tmp.end, buf, fetch_func)
    # finalize pileup
    bam_plbuf_push( NULL, buf)
    bam_index_destroy(idx);
    bam_plbuf_destroy(buf);

    samclose( tmp.input )
