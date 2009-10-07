# cython: embedsignature=True
# adds doc-strings for sphinx

# defines imported from samtools

DEF SEEK_SET = 0
DEF SEEK_CUR = 1
DEF SEEK_END = 2

## These are bits set in the flag.
## have to put these definitions here, in csamtools.pxd they got ignored
## @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
DEF BAM_FPAIRED       =1
## @abstract the read is mapped in a proper pair */
DEF BAM_FPROPER_PAIR  =2
## @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
DEF BAM_FUNMAP        =4
## @abstract the mate is unmapped */
DEF BAM_FMUNMAP       =8
## @abstract the read is mapped to the reverse strand */
DEF BAM_FREVERSE      =16
## @abstract the mate is mapped to the reverse strand */
DEF BAM_FMREVERSE     =32
## @abstract this is read1 */
DEF BAM_FREAD1        =64
## @abstract this is read2 */
DEF BAM_FREAD2       =128
## @abstract not primary alignment */
DEF BAM_FSECONDARY   =256
## @abstract QC failure */
DEF BAM_FQCFAIL      =512
## @abstract optical or PCR duplicate */
DEF BAM_FDUP        =1024

#####################################################################
#####################################################################
#####################################################################
## private methods
#####################################################################

cdef class AlignedRead
cdef samfile_t * openSam( char * filename, char * mode ):
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

## there might a better way to build an Aligment object,
## but I could not get the following to work, usually
## getting a "can not convert to Python object"
##   1. struct to dict conversion 
##   2. using "cdef class AlignedRead" with __cinit__(self, bam1_t * )
DEF BAM_CIGAR_SHIFT=4
DEF BAM_CIGAR_MASK=((1 << BAM_CIGAR_SHIFT) - 1)

cdef samtoolsToAlignedRead(AlignedRead dest, bam1_t * src):
    '''enter src into AlignedRead.'''
    dest._delegate = bam_dup1(src)
    return dest



cdef bam1_t * alignedReadToSamtools( bam1_t * dest, src ):
    '''convert to samtools data structure'''
    dest.core.tid = src.tid
    dest.core.pos = src.pos
    dest.core.bin = src.bin
    dest.core.qual = src.qual
    dest.core.l_qname = src.l_qname
    dest.core.flag = src.flag
    dest.core.n_cigar = src.n_cigar
    dest.core.l_qseq = src.l_qseq
    dest.data_len = src.data_len
    dest.l_aux = src.l_aux
    dest.m_data = src.m_data
    dest.core.mtid = src.mtid
    dest.core.mpos = src.mpos
    dest.core.isize = src.isize
    return dest

cdef samtoolsToPileupRead( dest, bam_pileup1_t * src ):
    '''fill a  PileupRead object from a bam_pileup1_t * object.'''
    dest.alignment = samtoolsToAlignedRead( AlignedRead(), src.b )
    dest.qpol = src.qpos
    dest.indel = src.indel
    dest.level = src.level
    dest.is_del = src.is_del
    dest.is_head = src.is_head
    dest.is_tail = src.is_tail
    return dest

#####################################################################
#####################################################################
#####################################################################
## Generic callbacks for inserting python callbacks.
#####################################################################
cdef int fetch_callback( bam1_t *alignment, void *f):
    '''callback for bam_fetch. 
    
    calls function in *f* with a new :class:`AlignedRead` object as parameter.
    '''
    a = samtoolsToAlignedRead( AlignedRead(), alignment )
    (<object>f)(a)

cdef int pileup_callback( uint32_t tid, uint32_t pos, int n, bam_pileup1_t *pl, void *f):
    '''callback for pileup.

    calls function in *f* with a new :class:`Pileup` object as parameter.

    tid
        chromosome ID as is defined in the header
    pos
        start coordinate of the alignment, 0-based
    n
        number of elements in pl array
    pl
        array of alignments
    data
        user provided data
    '''

    p = PileupColumn()
    p.tid = tid
    p.pos = pos
    p.n = n
    pileups = []

    for x from 0 <= x < n:
        pileups.append( samtoolsToPileupRead( PileupRead(), &(pl[x]) ) )
    p.pileups = pileups
        
    (<object>f)(p)

cdef int pileup_fetch_callback( bam1_t *b, void *data):
    '''callback for bam_fetch. 

    Fetches reads and submits them to pileup.
    '''
    cdef bam_plbuf_t * buf
    buf = <bam_plbuf_t*>data
    bam_plbuf_push(b, buf)
    return 0

######################################################################
######################################################################
######################################################################
## Public methods
######################################################################
cdef class Samfile:
    '''A bam/sam file
    '''

    cdef samfile_t * samfile
    cdef bam_index_t *index
    cdef char * filename
    cdef int isbam

    def __cinit__(self, *args, **kwargs ):
       self.samfile = NULL
       self.isbam = False
       self.open( *args, **kwargs )

    def isOpen( self ):
        '''return true if samfile has been opened.'''
        return self.samfile != NULL

    def hasIndex( self ):
        '''return true if samfile has an existing (and opened) index.'''
        return self.index != NULL

    def open( self, 
              char * filename, 
              mode,
              Samfile template = None,
              targetnames = None,
              targetlengths = None,
              ):
        '''open a sam/bam file. If an index exists, it will be opened as well.

        For writing, the names (*targetnames*) and lengths (*targetlengths*)
        of targets need to be supplied or the header is taken from a *template*.

        The *mode* corresponds to the syntax in samtools. Valid modes are
        ``/[rw](b?)(u?)(h?)/``: 
        
        r for reading
        w for writing, 
        b for BAM I/O, i.e., input/output is in binary format 
        u for uncompressed BAM output 
        h for outputing header in SAM 
         
        If ``b`` is present, it must immediately follow ``r`` or ``w``. 
        Currently valid modes are ``r``, ``w``, ``wh``, ``rb``, ``wb`` and ``wbu``.
        '''
        
        assert mode in ("r","w","rb","wb"), "invalid file opening mode `%s`" % mode

        ## TODO: implement automatic indexing
        cdef bam_header_t * header
        header = NULL

        if self.samfile != NULL: self.close()
        self.filename = filename

        self.isbam = len(mode) > 1 and mode[1] == 'b'

        if mode[0] == 'w':
            # open file for writing
            
            # header structure (used for writing)
            if self.isbam:

                if template:
                    # copy header from another file
                    header = template.samfile.header
                else:
                    # build header from a target names and lengths
                    assert targetnames and targetlengths, "supply names and lengths of targets for writing"
                    assert len(targetnames) == len(targetlengths), "unequal names and lengths of targets"

                    # allocate and fill header
                    header = bam_header_init()
                    header.n_targets = len(targetnames)
                    n = 0
                    for x in targetnames: n += len(x) + 1
                    header.target_name = <char**>calloc(n, sizeof(char*))
                    header.target_len = <uint32_t*>calloc(n, sizeof(uint32_t))
                    for x from 0 <= x < header.n_targets:
                        header.target_len[x] = targetlengths[x]
                        name = targetnames[x]
                        header.target_name[x] = <char*>calloc(len(name)+1, sizeof(char))
                        strncpy( header.target_name[x], name, len(name) )

                    header.l_text = 0
                    header.text = NULL
                    header.hash = NULL
                    header.rg2lib = NULL
                    
            # open file. Header gets written to file at the same time
            self.samfile = samopen( filename, mode, header )

            # bam_header_destroy takes care of cleaning up of all the members
            if not template and header != NULL:
                bam_header_destroy( header )

        elif mode[0] == "r":
            # open file for reading
            self.samfile = samopen( filename, mode, NULL )

        if self.samfile == NULL:
            raise IOError("could not open file `%s`" % filename )

        if mode[0] == "r" and self.isbam:
            self.index = bam_index_load(filename)
            if self.index == NULL:
                raise IOError("could not open index `%s` " % filename )

    def getTarget( self, tid ):
        '''convert numerical :term:`tid` into target name.'''
        if not 0 <= tid < self.samfile.header.n_targets:
            raise ValueError( "tid out of range 0<=tid<%i" % self.samfile.header.n_targets )
        return self.samfile.header.target_name[tid]

    def getNumTargets( self ):
        '''return the number of targets in file.'''
        return self.samfile.header.n_targets

    def getTargets( self ):
        '''return tuple of all targets.'''
        l = []
        for x from 0 <= x < self.samfile.header.n_targets:
            l.append( self.samfile.header.target_name[x] )
        return tuple(l)

    def getTargetLengths( self ):
        '''return tuple of all target lengths.'''
        l = []
        for x from 0 <= x < self.samfile.header.n_targets:
            l.append( self.samfile.header.target_len[x] )
        return tuple(l)

    def fetch(self, region, callback ):
        '''fetch a :term:`region` from the samfile.

        This method will execute a callback on each aligned read in
        a region.
        '''
        assert self.isbam, "fetch is only available for bam files"
        assert self.hasIndex(), "no index available for fetch"

        cdef int tid
        cdef int start
        cdef int end

        # parse the region
        bam_parse_region( self.samfile.header, region, &tid, &start, &end)
        if tid < 0:
            raise ValueError( "invalid region `%s`" % region )
        
        bam_fetch(self.samfile.x.bam, self.index, tid, start, end, <void*>callback, fetch_callback )

    def pileup( self, region, callback ):
        '''fetch a :term:`region` from the samfile and run a callback
        on each postition
        '''
        cdef int tid
        cdef int start
        cdef int end
        cdef bam_plbuf_t *buf

        if self.isbam:
            assert self.hasIndex(), "no index available for pileup"

            # parse the region
            bam_parse_region( self.samfile.header, region, &tid, &start, &end)
            if tid < 0: raise ValueError( "invalid region `%s`" % region )

            buf = bam_plbuf_init(pileup_callback, <void*>callback )

            bam_fetch(self.samfile.x.bam, self.index, tid, start, end, buf, pileup_fetch_callback )

            # finalize pileup
            bam_plbuf_push( NULL, buf)
            bam_plbuf_destroy(buf);

        else:
            raise NotImplementedError( "pileup of samfiles not wrapped" )

    def close( self ):
        '''close file.'''
        
        if self.samfile != NULL:
            samclose( self.samfile )
            bam_index_destroy(self.index);
            self.samfile = NULL

    def __dealloc__( self ):
        '''clean up.'''
        # Note that __del__ is not called.
        self.close()

    def write( self, read ):
        '''write read to file.

        return the number of bytes written.
        '''
        cdef bam1_t dest
        alignedReadToSamtools( &dest, read )
        return bam_write1( self.samfile.x.bam, &dest)

## turning callbacks elegantly into iterators is an unsolved problem, see the following threads:
## http://groups.google.com/group/comp.lang.python/browse_frm/thread/0ce55373f128aa4e/1d27a78ca6408134?hl=en&pli=1
## http://www.velocityreviews.com/forums/t359277-turning-a-callback-function-into-a-generator.html
## Thus I chose to rewrite the functions requiring callbacks. The downside is that if the samtools C-API or code
## changes, the changes have to be manually entered.

cdef class IteratorRow:
    """iterate over mapped reads in a region.
    """
    
    cdef bam_fetch_iterator_t*  bam_iter # iterator state object
    cdef bam1_t *               b
    cdef                        error_msg
    cdef int                    error_state
    def __cinit__(self, Samfile samfile, region ):
        self.bam_iter = NULL

        assert samfile.isOpen()
        assert samfile.hasIndex()

        # parse the region
        cdef int      tid
        cdef int      beg
        cdef int      end
        self.error_state = 0
        self.error_msg = None
        bam_parse_region( samfile.samfile.header, region, &tid, &beg, &end)
        if tid < 0: 
            self.error_state = 1
            self.error_msg = "invalid region `%s`" % region
            return

        cdef bamFile  fp
        fp = samfile.samfile.x.bam
        self.bam_iter = bam_init_fetch_iterator(fp, samfile.index, tid, beg, end)

    def __iter__(self):
        return self 

    cdef bam1_t * getCurrent( self ):
        return self.b

    def __next__(self): 
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        if self.error_state:
            raise ValueError( self.error_msg)
        
        self.b = bam_fetch_iterate(self.bam_iter)
        if self.b <> NULL:
            dest = samtoolsToAlignedRead( AlignedRead(), self.b )
            return dest
        else:
            raise StopIteration

    def __dealloc__(self):
        '''remember: dealloc cannot call other methods!'''
        if self.bam_iter:
            bam_cleanup_fetch_iterator(self.bam_iter)
        
cdef class IteratorRowAll:
    """iterate over all mapped reads
    """

    cdef bam1_t *               b
    cdef samfile_t *            fp

    def __cinit__(self, Samfile samfile):

        assert samfile.isOpen()

        # parse the region
        self.fp = samfile.samfile

        # allocate memory for alignment
        self.b = <bam1_t*>calloc(1, sizeof(bam1_t))

        # seek past the header not necessary for bam files

    def __iter__(self):
        return self 

    cdef bam1_t * getCurrent( self ):
        return self.b

    def __next__(self): 
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        cdef int ret
        ret = samread(self.fp, self.b)
        if (ret > 0):
            return samtoolsToAlignedRead( AlignedRead(), self.b )
        else:
            raise StopIteration

    def __dealloc__(self):
        '''remember: dealloc cannot call other methods!'''
        bam_destroy1(self.b);
        

cdef class IteratorColumn:
    '''iterate over columns.

    This iterator wraps the pileup functionality of the samtools
    function bam_plbuf_push.
    '''
    cdef bam_plbuf_t *buf

    # check if first iteration
    cdef int notfirst
    cdef int n_pu
    cdef int eof 
    cdef IteratorRow iter

    def __cinit__(self, Samfile samfile, region ):

        self.iter = IteratorRow( samfile, region )
        self.buf = bam_plbuf_init(NULL, NULL )
        self.n_pu = 0
        self.eof = 0

    def __iter__(self):
        return self 

    cdef int cnext(self):

        cdef int retval1, retval2

        # check if previous plbuf was incomplete. If so, continue within
        # the loop and yield if necessary
        if self.n_pu > 0:
            self.n_pu = pysam_bam_plbuf_push( self.iter.getCurrent(), self.buf, 1)
            if self.n_pu > 0: return 1

        if self.eof: return 0

        # get next alignments and submit until plbuf indicates that
        # an new column has finished
        while self.n_pu == 0:
            retval1 = self.iter.cnext()
            # wrap up if no more input
            if retval1 == 0: 
                self.n_pu = pysam_bam_plbuf_push( NULL, self.buf, 0)            
                self.eof = 1
                return 1

            # submit to plbuf
            self.n_pu = pysam_bam_plbuf_push( self.iter.getCurrent(), self.buf, 0)            
            if self.n_pu < 0: raise ValueError( "error while iterating" )

        # plbuf has yielded
        return 1

    def __next__(self): 
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        cdef int ret
        ret = self.cnext()
        cdef bam_pileup1_t * pl

        if ret > 0 :
            p = PileupColumn()
            p.tid = pysam_get_tid( self.buf )
            p.pos = pysam_get_pos( self.buf )
            p.n = self.n_pu
            pl = pysam_get_pileup( self.buf )

            pileups = []
            
            for x from 0 <= x < p.n:
                pileups.append( samtoolsToPileupRead( PileupRead(), &pl[x]) )
            p.pileups = pileups

            return p
        else:
            raise StopIteration

    def __dealloc__(self):
        bam_plbuf_destroy(self.buf);

cdef class AlignedRead:
    '''
    todo:update
    Python wrapper around an aligned read. The following
    fields are exported:

    tid
        chromosome/target ID
    pos
        0-based leftmost coordinate
    bin
        bin calculated by bam_reg2bin()
    qual
        mapping quality
    l_qname
        length of the query name
    flag
        bitwise flag
    n_cigar
        number of CIGAR operations
    l_qseq
        length of the query sequence (read)
    qname
        the query name (None if not present)
    cigar
        the :term:`cigar` alignment string (None if not present)
    qseq
        the query sequence (None if not present)
    bqual
        the base quality (None if not present)
    mtid
        the :term:`target` id of the mate 
    mpos  
        the position of the mate
    izise
        the insert size
    '''

    def __dealloc__(self):
        """todo is this enough or do we need to free() each string? eg 'qual' etc"""
        bam_destroy1(self._delegate)
            
    cdef:
         bam1_t * _delegate 
         
    
    def __str__(self):
        """todo"""
        return "\t".join(map(str, (self.qname,
                                    self.rname,
                                    self.pos,
                                    self.qual,
                                    self.flag,
                                    self.seq,
                                    self.mapq , self.opt('MF'), self.is_paired)))
                                 
 
                                    
    #For all the following we follow the samfile specification for field names, changing uppercase to  lowercase only
    
    property cigar:
        def __get__(self):
            cdef uint32_t * cigar_p
            cdef bam1_t * src 
            cdef op, l, cigar
            src = self._delegate
            if src.core.n_cigar > 0:
                cigar = []
                cigar_p = < uint32_t *> (src.data + src.core.l_qname)
                for k from 0 <= k < src.core.n_cigar:
                    op = cigar_p[k] & BAM_CIGAR_MASK
                    l = cigar_p[k] >> BAM_CIGAR_SHIFT
                    cigar.append((op, l))
                return cigar
            return None
    property qname:
        def __get__(self):
            cdef bam1_t * src 
            src = self._delegate
            ## parse qname (bam1_qname)
            return < char *> src.data
    property seq:
        def __get__(self):
            cdef bam1_t * src
            cdef uint8_t * seq_p 
            cdef char * s
            src = self._delegate
            bam_nt16_rev_table = "=ACMGRSVTWYHKDBN"
            ## parse qseq (bam1_seq)
            if src.core.l_qseq: 
                s = < char *> calloc(src.core.l_qseq + 1 , sizeof(char))
                seq_p = < uint8_t *> (src.data + src.core.n_cigar * 4 + src.core.l_qname)
                for k from 0 <= k < src.core.l_qseq:
                ## equivalent to bam_nt16_rev_table[bam1_seqi(s, i)] (see bam.c)
                    s[k] = "=ACMGRSVTWYHKDBN"[((seq_p)[(k) / 2] >> 4 * (1 - (k) % 2) & 0xf)]
                return s
            return None
    property qual:
        def __get__(self):
            #todo is this still needed
            cdef bam1_t * src 
            cdef uint8_t * qual_p
            cdef char * q
            src = self._delegate
            qual_p = < uint8_t *> (src.data + src.core.n_cigar * 4 + src.core.l_qname + (src.core.l_qseq + 1) / 2)        
            if qual_p[0] != 0xff:
                q = < char *> calloc(src.core.l_qseq + 1 , sizeof(char))
                for k from 0 <= k < src.core.l_qseq:
                ## equivalent to t[i] + 33 (see bam.c)
                    q[k] = qual_p[k] + 33
                return q
            return None


    #todo check which of these are needed (see specification) do we need l_qseq, n_cigar etc? also check tid=rname
    property flag: 
        def __get__(self): 
            return self._delegate.core.flag
    property rname: 
        def __get__(self): 
            return self._delegate.core.tid
    property pos: 
        def __get__(self): 
            return self._delegate.core.pos
    property mapq: 
        def __get__(self): 
            return self._delegate.core.qual
    property mrnm: 
        def __get__(self): 
            return self._delegate.core.mtid
    property mpos: 
        def __get__(self): 
            return self._delegate.core.mpos
    property isize: 
        def __get__(self): 
            return self._delegate.core.isize
    
    
    #todo: finish these off
    property is_paired: 
        def __get__(self): 
            """true if read is paired in sequencing"""
            return (self._delegate.core.flag & BAM_FPAIRED) != 0
    property is_proper_pair:
        """true if read is mapped in a proper pair"""
        def __get__(self): return (self.flag & BAM_FPROPER_PAIR) != 0
    #def _getfunmap( self ): return (self.flag & BAM_FUNMAP) != 0
    #is_unmapped = property( _getfunmap, doc="true if read itself is unmapped" )
    #def _getfmunmap( self ): return (self.flag & BAM_FMUNMAP) != 0
    #mate_is_unmapped = property( _getfmunmap, doc="true if the mate is unmapped" )
    #def _getreverse( self ): return (self.flag & BAM_FREVERSE) != 0
    #is_reverse = property( _getreverse, doc = "true if read is mapped to reverse strand" )
    #def _getmreverse( self ): return (self.flag & BAM_FMREVERSE) != 0
    #mate_is_reverse = property( _getmreverse, doc= "true is read is mapped to reverse strand" )
    #def _getread1( self ): return (self.flag & BAM_FREAD1) != 0
    #is_read1 = property( _getread1, doc = "true if this is read1" )
    #def _getread2( self ): return (self.flag & BAM_FREAD2) != 0
    #is_read2 = property( _getread2, doc = "true if this is read2" )
    #def _getfsecondary( self ): return (self.flag & BAM_FSECONDARY) != 0
    #is_secondary = property( _getfsecondary, doc="true if not primary alignment" )
    #def _getfqcfail( self ): return (self.flag & BAM_FSECONDARY) != 0
    #is_qcfail = property( _getfqcfail, doc="true if QC failure" )
    #def _getfdup( self ): return (self.flag & BAM_FDUP) != 0
    #is_duplicate = property( _getfdup, doc="true if optical or PCR duplicate" )
    
    def opt(self, tag):
        """retrieves optional data given a two-letter *tag*"""
        #see bam_aux_get.c: bam_aux_get() and bam_aux2i() etc 
        cdef uint8_t * v
        v = bam_aux_get(self._delegate, tag)
        type = chr(v[0])
        if type == 'c' or type == 'C' or type == 's' or type == 'S' or type == 'i':
            return bam_aux2i(v)            
        elif type == 'f':
            return bam_aux2f(v)
        elif type == 'd':
            return bam_aux2d(v)
        elif type == 'A':
            return bam_aux2A(v)
        elif type == 'Z':
            return bam_aux2Z(v)
    
    def fancy_str (self):
        """
        returns list of fieldnames/values in pretty format for debugging
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


                        "l_aux":         "Length of auxilary data",
                        "m_data":        "Maximum data length",
                        "data_len":      "Current data length",
                        }
        fields_names_in_order = ["tid", "pos", "mtid", "mpos", "isize", "flag", 
                                 "n_cigar", "cigar", "qual", "bin", "l_qname", "qname", 
                                 "l_qseq", "qseq", "bqual", "l_aux", "m_data", "data_len"]

        for f in fields_names_in_order:
            if not f in self.__dict__:
                continue
            ret_string.append("%-30s %-10s= %s" % (field_names[f], "(" + f + ")", self.__getattribute__(f)))

        for f in self.__dict__:
            if not f in field_names:
                ret_string.append("%-30s %-10s= %s" % (f, "", self.__getattribute__(f)))
        return ret_string

        

   


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
            "\n" +\
            "\n".join( map(str, self.pileups) )

class PileupRead(object):
    '''A read aligned to a column.

    alignment
       a :class:`pysam.AlignedRead` object of the aligned read
    qpos
       position of the read base at the pileup site, 0-based
    indel
       indel length; 0 for no indel, positive for ins and negative for del
    is_del
       1 iff the base on the padded read is a deletion
    level
       the level of the read in the "viewer" mode 
    '''

    def __init__(self, **kwargs ):
        
        ## there must be an easier way to construct this class, but it won't accept
        ## a typed parameter.
        self.alignment = kwargs.get( "alignment", AlignedRead() )
        self.qpos      = kwargs.get("qpos", 0)
        self.indel     = kwargs.get("indel", 0)
        self.level     = kwargs.get("level", 0)
        self.is_del    = kwargs.get("is_del", 0)
        self.is_head   = kwargs.get("is_head", 0)
        self.is_tail   = kwargs.get("is_tail", 0)

    def __str__(self):
        return "\t".join( map(str, (self.alignment, self.qpos, self.indel, self.level, self.is_del, self.is_head, self.is_tail ) ) )

import tempfile, os, sys

class Outs:
    '''http://mail.python.org/pipermail/python-list/2000-June/038406.html'''
    def __init__(self, id = 1):
        self.streams = []
        self.id = id

    def setdevice(self, filename):
        '''open an existing file, like "/dev/null"'''
        fd = os.open(filename, os.O_WRONLY)
        self.setfd(fd)

    def setfile(self, filename):
        '''open a new file.'''
        fd = os.open(filename, os.O_WRONLY|os.O_CREAT, 0660);
        self.setfd(fd)

    def setfd(self, fd):
        ofd = os.dup(self.id)      #  Save old stream on new unit.
        self.streams.append(ofd)
        sys.stdout.flush()    #  Buffered data goes to old stream.
        os.dup2(fd, self.id)        #  Open unit 1 on new stream.
        os.close(fd)          #  Close other unit (look out, caller.)
            
    def restore(self):
        '''restore previous output stream'''
        if self.streams:
            # the following was not sufficient, hence flush both stderr and stdout
            # os.fsync( self.id )
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(self.streams[-1], self.id)
            os.close(self.streams[-1])
            del self.streams[-1]

def _samtools_dispatch( method, args = () ):
    '''call ``method`` in samtools providing arguments in args.
    
    .. note:: 
       This method redirects stdout and stderr to capture it 
       from samtools. If for some reason stdout/stderr disappears
       the reason might be in this method.

    .. note::
       The current implementation might only work on linux.
       
    .. note:: 
       This method captures stdout and stderr using temporary files, 
       which are then read into memory in their entirety. This method
       is slow and might cause large memory overhead. 

    See http://bytes.com/topic/c/answers/487231-how-capture-stdout-temporarily
    on the topic of redirecting stderr/stdout.
    '''

    # note that debugging this module can be a problem
    # as stdout/stderr will not appear

    # redirect stderr and stdout to file

    # open files and redirect into it
    stderr_h, stderr_f = tempfile.mkstemp()
    stdout_h, stdout_f = tempfile.mkstemp()

    stdout_save = Outs( sys.stdout.fileno() )
    stdout_save.setfd( stdout_h )
    stderr_save = Outs( sys.stderr.fileno() )
    stderr_save.setfd( stderr_h )

    # do the function call to samtools
    cdef char ** cargs
    cdef int i, n, retval

    n = len(args)
    # allocate two more for first (dummy) argument (contains command)
    cargs = <char**>calloc( n+2, sizeof( char *) )
    cargs[0] = "samtools"
    cargs[1] = method
    for i from 0 <= i < n: cargs[i+2] = args[i]
    retval = pysam_dispatch(n+2, cargs)
    free( cargs )

    # restore stdout/stderr. This will also flush, so
    # needs to be before reading back the file contents
    stdout_save.restore()
    stderr_save.restore()

    # capture stderr/stdout.
    out_stderr = open( stderr_f, "r").readlines()
    out_stdout = open( stdout_f, "r").readlines()

    # clean up files
    os.remove( stderr_f )
    os.remove( stdout_f )

    return retval, out_stderr, out_stdout

__all__ = ["Samfile", "IteratorRow", "IteratorRowAll", "IteratorColumn", "AlignedRead", "PileupColumn", "PileupRead" ]

               

