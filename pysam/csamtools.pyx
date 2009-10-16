# cython: embedsignature=True
# adds doc-strings for sphinx

import tempfile, os, sys, types, itertools

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

DEF BAM_CIGAR_SHIFT=4
DEF BAM_CIGAR_MASK=((1 << BAM_CIGAR_SHIFT) - 1)

#####################################################################
#####################################################################
#####################################################################
## private factory methods
#####################################################################
cdef class AlignedRead
cdef makeAlignedRead( bam1_t * src):
    '''enter src into AlignedRead.'''
    cdef AlignedRead dest
    dest = AlignedRead()
    dest._delegate = bam_dup1(src)
    return dest

cdef class PileupRead
cdef makePileupRead( bam_pileup1_t * src ):
    '''fill a  PileupRead object from a bam_pileup1_t * object.'''
    cdef PileupRead dest
    dest = PileupRead()
    dest._alignment = makeAlignedRead( src.b )
    dest._qpos = src.qpos
    dest._indel = src.indel
    dest._level = src.level
    dest._is_del = src.is_del
    dest._is_head = src.is_head
    dest._is_tail = src.is_tail
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
    a = makeAlignedRead( alignment )
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
        pileups.append( makePileupRead( &(pl[x]) ) )
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
# valid types for sam headers
VALID_HEADER_TYPES = { "HD" : dict, 
                       "SQ" : list, 
                       "RG" : list, 
                       "PG" : dict, 
                       "CO" : list }

# order of records within sam headers
VALID_HEADERS = ("HD", "SQ", "RG", "PG", "CO" )

# type conversions within sam header records
VALID_HEADER_FIELDS = { "HD" : { "VN" : str, "SO" : str, "GO" : str },
                        "SQ" : { "SN" : str, "LN" : int, "AS" : str, "M5" : str, "UR" : str, "SP" : str },
                        "RG" : { "ID" : str, "SM" : str, "LB" : str, "DS" : str, "PU" : str, "PI" : str, "CN" : str, "DT" : str, "PL" : str, },
                        "PG" : { "ID" : str, "VN" : str, "CL" : str }, }

# output order of fields within records
VALID_HEADER_ORDER = { "HD" : ( "VN", "SO", "GO" ),
                       "SQ" : ( "SN", "LN", "AS", "M5" , "UR" , "SP" ),
                       "RG" : ( "ID", "SM", "LB", "DS" , "PU" , "PI" , "CN" , "DT", "PL" ),
                       "PG" : ( "ID", "VN", "CL" ), }

######################################################################
######################################################################
######################################################################
## Public methods
######################################################################
cdef class Samfile:
    '''*(filename, mode='r', template = None, referencenames = None, referencelengths = None, text = NULL, header = None)*
              
    A *SAM* file. The file is automatically opened.

    
    *mode* should be ``r`` for reading or ``w`` for writing. The default is text mode so for binary (:term:`BAM`) I/O you should append 
    ``b`` for compressed or ``u`` for uncompressed :term:`BAM` output. Use ``h`` to output header information  in text (:term:`TAM`)  mode.

    If ``b`` is present, it must immediately follow ``r`` or ``w``. 
    Currently valid modes are ``r``, ``w``, ``wh``, ``rb``, ``wb`` and ``wbu``.
    
    so to open a :term:`BAM` file for reading::

        f=Samfile('ex1.bam','rb')
    
    

    For writing, the header of a :term:`TAM` file/:term:`BAM` file can be constituted from several
    sources:

        1. If *template* is given, the header is copied from a another *Samfile* (*template* must be of type *Samfile*).

        2. If *header* is given, the header is build from a multi-level dictionary. The first level are the four types ('HD', 'SQ', ...). The second level is then a list of lines, with each line being a list of tag-value pairs.

        3. If *text* is given, new header text is copied from raw text.

        4. The names (*referencenames*) and lengths (*referencelengths*) are supplied directly as lists. 

    
    

    '''

    cdef char * filename
    # pointer to samfile
    cdef samfile_t * samfile
    # pointer to index
    cdef bam_index_t *index
    # true if file is a bam file
    cdef int isbam

    def __cinit__(self, *args, **kwargs ):
        self.samfile = NULL
        self.isbam = False
        self._open( *args, **kwargs )

    def _isOpen( self ):
        '''return true if samfile has been opened.'''
        return self.samfile != NULL

    def _hasIndex( self ):
        '''return true if samfile has an existing (and opened) index.'''
        return self.index != NULL

    def _open( self, 
              char * filename, 
              mode ='r',
              Samfile template = None,
              referencenames = None,
              referencelengths = None,
              char * text = NULL,
              header = None,
              ):
        '''open a sam/bam file.

        If _open is called on an existing bamfile, the current file will be
        closed and a new file will be opened.
        '''

        assert mode in ( "r","w","rb","wb", "wh", "wbu" ), "invalid file opening mode `%s`" % mode

        # close a previously opened file
        if self.samfile != NULL: self.close()

        cdef bam_header_t * header_to_write
        header_to_write = NULL

        if self.samfile != NULL: self.close()
        self.filename = filename

        self.isbam = len(mode) > 1 and mode[1] == 'b'

        if mode[0] == 'w':
            # open file for writing
            
            # header structure (used for writing)
            if template:
                # copy header from another file
                header_to_write = template.samfile.header

            elif header:
                header_to_write = self._buildHeader( header )
            else:
                # build header from a target names and lengths
                assert referencenames and referencelengths, "supply names and lengths of reference sequences for writing"
                assert len(referencenames) == len(referencelengths), "unequal names and lengths of reference sequences"

                # allocate and fill header
                header_to_write = bam_header_init()
                header_to_write.n_targets = len(referencenames)
                n = 0
                for x in referencenames: n += len(x) + 1
                header_to_write.target_name = <char**>calloc(n, sizeof(char*))
                header_to_write.target_len = <uint32_t*>calloc(n, sizeof(uint32_t))
                for x from 0 <= x < header_to_write.n_targets:
                    header_to_write.target_len[x] = referencelengths[x]
                    name = referencenames[x]
                    header_to_write.target_name[x] = <char*>calloc(len(name)+1, sizeof(char))
                    strncpy( header_to_write.target_name[x], name, len(name) )

                if text != NULL:
                    # copy without \0
                    header_to_write.l_text = strlen(text)
                    header_to_write.text = <char*>calloc( strlen(text), sizeof(char) )
                    memcpy( header_to_write.text, text, strlen(text) )

                header_to_write.hash = NULL
                header_to_write.rg2lib = NULL
                    
            # open file. Header gets written to file at the same time for bam files
            # and sam files (in the latter case, the mode needs to be wh)
            self.samfile = samopen( filename, mode, header_to_write )

            # bam_header_destroy takes care of cleaning up of all the members
            if not template and header_to_write != NULL:
                bam_header_destroy( header_to_write )

        elif mode[0] == "r":
            # open file for reading
            self.samfile = samopen( filename, mode, NULL )

        if self.samfile == NULL:
            raise IOError("could not open file `%s`" % filename )

        if mode[0] == "r" and self.isbam:
            self.index = bam_index_load(filename)
            if self.index == NULL:
                raise IOError("could not open index `%s` " % filename )

    def _getTarget( self, tid ):
        '''(tid )
        convert numerical :term:`tid` into target name.'''
        if not 0 <= tid < self.samfile.header.n_targets:
            raise ValueError( "tid out of range 0<=tid<%i" % self.samfile.header.n_targets )
        return self.samfile.header.target_name[tid]

    def getNumReferences( self ):
        """return the number of :term:`reference` sequences."""
        return self.samfile.header.n_targets

    def _parseRegion( self, reference = None, start = None, end = None, 
                       region = None ):
        '''parse region information.

        raise Value for for invalid regions.

        returns a tuple of region, tid, start and end. Region
        is a valid samtools :term:`region` or None if the region
        extends over the whole file.
        '''
        
        cdef int rtid
        cdef int rstart
        cdef int rend

        rtid = rstart = rend = 0

        # translate to a region
        if reference:
            if start != None and end != None:
                region = "%s:%i-%i" % (reference, start, end)
            else:
                region = reference

        if region:
            bam_parse_region( self.samfile.header, region, &rtid, &rstart, &rend)        
            if rtid < 0: raise ValueError( "invalid region `%s`" % region )

            if rstart > rend: raise ValueError( 'invalid region: start (%i) > end (%i)' % (rstart, rend) )
            if rstart < 0: raise ValueError( 'negative start coordinate (%i)' % rstart )
            if rend < 0: raise ValueError( 'negative end coordinate (%i)' % rend )
            
        return region, rtid, rstart, rend

    def fetch( self, 
               reference = None, start = None, end = None, 
               region = None, 
               callback = None,
               until_eof = False ):
        '''*(reference = None, start = None, end = None, region = None, callback = None, until_eof = False)*
               
        fetch :meth:`AlignedRead` objects in a :term:`region` using 0-based indexing. The region is specified by
        :term:`reference`, *start* and *end*. Alternatively, a samtools :term:`region` string can be supplied.

        Without *reference* or *region* all reads will be fetched. The reads will be returned
        ordered by reference sequence, which will not necessarily be the order within the file.
        If *until_eof* is given, all reads from the current file position will be returned
        *as they are sorted within the file*.  
        
        If only *reference* is set, all reads matching on *reference* will be fetched.

        The method returns an iterator of type :class:`pysam.IteratorRow` unless
        a *callback is provided. If *callback* is given, the callback will be executed 
        for each position within the :term:`region`. Note that callbacks currently work
        only, if *region* or *reference* is given.

        Note that a :term:`TAM` file does not allow random access. If *region* or *reference* are given,
        an exception is raised.
        '''
        cdef int rtid
        cdef int rstart
        cdef int rend

        region, rtid, rstart, rend = self._parseRegion( reference, start, end, region )

        if self.isbam:
            assert self._hasIndex(), "no index available for fetch"            
            if callback:
                if not region:
                    raise ValueError( "callback functionality requires a region/reference" )
                return bam_fetch(self.samfile.x.bam, self.index, rtid, rstart, rend, <void*>callback, fetch_callback )
            else:
                if region:
                    return IteratorRow( self, region )
                else:
                    if until_eof:
                        return IteratorRowAll( self )
                    else:
                        # return all targets by chaining the individual targets together.
                        i = []
                        for x in self.references: i.append( IteratorRow( self, x))
                        return itertools.chain( *i )
        else:                    
            if region != None:
                raise ValueError ("fetch for a region is not available for sam files" )
            if callback:
                raise NotImplementedError( "callback not implemented yet" )
            else:
                return IteratorRowAll( self )

    def pileup( self, reference = None, start = None, end = None, region = None, callback = None ):
        '''run a pileup within a :term:`region` using 0-based indexing. The region is specified by
        :term:`reference`, *start* and *end*. Alternatively, a samtools *region* string can be supplied.

        Without *reference* or *region* all reads will be fetched. The reads will be returned
        ordered by :term:`reference` sequence, which will not necessarily be the order within the file.

        The method returns an iterator of type :class:`pysam.IteratorColumn` unless
        a *callback is provided. If *callback* is given, the callback will be executed 
        for each position within the :term:`region`. 

        Note that samfiles do not allow random access. If *region* or *reference* are given,
        an exception is raised.
        
        .. Note::
            Note that *all* reads which overlap the region are returned. So the first base returned will be the first base of the first read *not* the first base of the region.
        '''
        cdef int rtid
        cdef int rstart
        cdef int rend
        cdef bam_plbuf_t *buf

        region, rtid, rstart, rend = self._parseRegion( reference, start, end, region )
        
        if self.isbam:
            assert self._hasIndex(), "no index available for pileup"

            if callback:
                if not region:
                    raise ValueError( "callback functionality requires a region/reference" )

                buf = bam_plbuf_init(pileup_callback, <void*>callback )
                bam_fetch(self.samfile.x.bam, 
                          self.index, rtid, rstart, rend, 
                          buf, pileup_fetch_callback )
                
                # finalize pileup
                bam_plbuf_push( NULL, buf)
                bam_plbuf_destroy(buf)
            else:
                if region:
                    return IteratorColumn( self, region )
                else:
                    # return all targets by chaining the individual targets together.
                    i = []
                    for x in self.references: i.append( IteratorColumn( self, x))
                    return itertools.chain( *i )

        else:
            raise NotImplementedError( "pileup of samfiles not implemented yet" )

    def close( self ):
        '''closes file.'''
        
        if self.samfile != NULL:
            samclose( self.samfile )
            bam_index_destroy(self.index);
            self.samfile = NULL

    def __dealloc__( self ):
        '''clean up.'''
        # Note that __del__ is not called.
        self.close()

    def write( self, AlignedRead read ):
        '''(AlignedRead read )
        write a single :class:`pysam.AlignedRead`..

        return the number of bytes written.
        '''
        return samwrite( self.samfile, read._delegate )

    property references:
        """tuple with the names of :term:`reference` sequences."""
        def __get__(self): 
            t = []
            for x from 0 <= x < self.samfile.header.n_targets:
                t.append( self.samfile.header.target_name[x] )
            return tuple(t)

    property lengths:
        """tuple of the lengths of the :term:`reference` sequences. The lengths are in the same order as :attr:`pysam.Samfile.reference`
        """
        def __get__(self): 
            t = []
            for x from 0 <= x < self.samfile.header.n_targets:
                t.append( self.samfile.header.target_len[x] )
            return tuple(t)

    property text:
        '''full contents of the :term:`sam file` header as a string.'''
        def __get__(self):
            # create a temporary 0-terminated copy
            cdef char * t
            t = <char*>calloc( self.samfile.header.l_text + 1, sizeof(char) )
            memcpy( t, self.samfile.header.text, self.samfile.header.l_text )
            result = t
            free(t)
            return result

    property header:
        '''header information within the :term:`sam file`. The records and fields are returned as 
        a two-level dictionary.
        '''
        def __get__(self):
            result = {}

            if self.samfile.header.text != NULL:
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
                        key,value = field.split(":")
                        if key not in VALID_HEADER_FIELDS[record]:
                            raise ValueError( "unknown field code '%s' in record '%s'" % (key, record) )
                        x[key] = VALID_HEADER_FIELDS[record][key](value)

                    if VALID_HEADER_TYPES[record] == dict:
                        if record in result:
                            raise ValueError( "multiple '%s' lines are not permitted" % record )
                        result[record] = x
                    elif VALID_HEADER_TYPES[record] == list:
                        if record not in result: result[record] = []
                        result[record].append( x )

            return result

    def _buildLine( self, fields, record ):
        '''build a header line from *fields* dictionary for *record*'''

        # TODO: add checking for field and sort order
        line = ["@%s" % record ]
        if record == "CO":
            line.append( fields )
        else:
            for key in VALID_HEADER_ORDER[record]:
                if key in fields:
                    line.append( "%s:%s" % (key, str(fields[key])))
        return "\t".join( line ) 

    cdef bam_header_t * _buildHeader( self, new_header ):
        '''return a new header built from a dictionary in *new_header*.

        This method inserts the text field, target_name and target_len.
        '''

        lines = []

        # check if hash exists

        # create new header and copy old data
        cdef bam_header_t * dest

        dest = bam_header_init()
                
        for record in VALID_HEADERS:
            if record in new_header:
                ttype = VALID_HEADER_TYPES[record]
                data = new_header[record]
                if type( data ) != type( ttype() ):
                    raise ValueError( "invalid type for record %s: %s, expected %s" % (record, type(data), type(ttype()) ) )
                if type( data ) == types.DictType:
                    lines.append( self._buildLine( data, record ) )
                else:
                    for fields in new_header[record]:
                        lines.append( self._buildLine( fields, record ) )

        text = "\n".join(lines) + "\n"
        if dest.text != NULL: free( dest.text )
        dest.text = <char*>calloc( len(text), sizeof(char))
        dest.l_text = len(text)
        strncpy( dest.text, text, dest.l_text )

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
                strncpy( dest.target_name[x], seqname, len(seqname) + 1 )
                dest.target_len[x] = seqlen

        return dest

## turning callbacks elegantly into iterators is an unsolved problem, see the following threads:
## http://groups.google.com/group/comp.lang.python/browse_frm/thread/0ce55373f128aa4e/1d27a78ca6408134?hl=en&pli=1
## http://www.velocityreviews.com/forums/t359277-turning-a-callback-function-into-a-generator.html
## Thus I chose to rewrite the functions requiring callbacks. The downside is that if the samtools C-API or code
## changes, the changes have to be manually entered.

cdef class IteratorRow:
    """iterates over mapped reads in a region.
    """
    
    cdef bam_fetch_iterator_t*  bam_iter # iterator state object
    cdef bam1_t *               b
    cdef                        error_msg
    cdef int                    error_state
    def __cinit__(self, Samfile samfile, region ):
        self.bam_iter = NULL

        assert samfile._isOpen()
        assert samfile._hasIndex()

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

    cdef int cnext(self):
        '''cversion of iterator. Used by IteratorColumn'''
        self.b = bam_fetch_iterate(self.bam_iter)
        if self.b == NULL: return 0
        return 1

    def __next__(self): 
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        if self.error_state:
            raise ValueError( self.error_msg)
        
        self.b = bam_fetch_iterate(self.bam_iter)
        if self.b != NULL:
            return makeAlignedRead( self.b )
        else:
            raise StopIteration

    def __dealloc__(self):
        '''remember: dealloc cannot call other methods!'''
        if self.bam_iter:
            bam_cleanup_fetch_iterator(self.bam_iter)
        
cdef class IteratorRowAll:
    """iterates over all mapped reads
    """

    cdef bam1_t * b
    cdef samfile_t * fp

    def __cinit__(self, Samfile samfile):

        assert samfile._isOpen()

        self.fp = samfile.samfile

        # allocate memory for alignment
        self.b = <bam1_t*>calloc(1, sizeof(bam1_t))

    def __iter__(self):
        return self 

    cdef bam1_t * getCurrent( self ):
        return self.b

    cdef int cnext(self):
        '''cversion of iterator. Used by IteratorColumn'''
        cdef int ret
        return samread(self.fp, self.b)

    def __next__(self): 
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        cdef int ret
        ret = samread(self.fp, self.b)
        if (ret > 0):
            return makeAlignedRead( self.b )
        else:
            raise StopIteration

    def __dealloc__(self):
        '''remember: dealloc cannot call other methods!'''
        bam_destroy1(self.b)
        
cdef class IteratorColumn:
    '''iterates over columns.

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
                pileups.append( makePileupRead( &pl[x]) )
            p.pileups = pileups

            return p
        else:
            raise StopIteration

    def __dealloc__(self):
        bam_plbuf_destroy(self.buf);

cdef class AlignedRead:
    '''
    Class representing an aligned read. see SAM format specification for meaning of fields (http://samtools.sourceforge.net/).
        
    '''
    cdef:
         bam1_t * _delegate 

    def __cinit__( self ):
        self._delegate = <bam1_t*>calloc( sizeof( bam1_t), 1 )

    def __dealloc__(self):
        """todo is this enough or do we need to free() each string? eg 'qual' etc"""
        bam_destroy1(self._delegate)
    
    def __str__(self):
        """todo"""
        return "\t".join(map(str, (self.qname,
                                   self.rname,
                                   self.pos,
                                   self.qual,
                                   self.flag,
                                   self.seq,
                                   self.mapq )))
    
    property cigar:
        """the :term:`cigar` alignment string (None if not present)"""
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
        """the query name (None if not present)"""
        def __get__(self):
            cdef bam1_t * src 
            src = self._delegate
            ## parse qname (bam1_qname)
            return < char *> src.data
    property seq:
        """the query sequence (None if not present)"""
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
                retval=s
                free(s)
                return retval
            return None
    property qual:
        """the base quality (None if not present)"""
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
                retval=q
                free(q)
                return retval
            return None

    property flag: 
        """properties flag"""
        def __get__(self): 
            return self._delegate.core.flag
    property rname: 
        """chromosome/target ID"""
        def __get__(self): 
            return self._delegate.core.tid
    property pos: 
        """0-based leftmost coordinate"""
        def __get__(self): 
            return self._delegate.core.pos
    property mapq: 
        """mapping quality"""
        def __get__(self): 
            return self._delegate.core.qual
    property mrnm:
        """the :term:`reference` id of the mate """     
        def __get__(self): 
            return self._delegate.core.mtid
    property mpos: 
        """the position of the mate"""
        def __get__(self): 
            return self._delegate.core.mpos
    property isize: 
        """the insert size"""
        def __get__(self): 
            return self._delegate.core.isize
    property is_paired: 
        """true if read is paired in sequencing"""
        def __get__(self): return (self._delegate.core.flag & BAM_FPAIRED) != 0
    property is_proper_pair:
        """true if read is mapped in a proper pair"""
        def __get__(self): return (self.flag & BAM_FPROPER_PAIR) != 0
    property is_unmapped:
        """true if read itself is unmapped"""
        def __get__(self): return (self.flag & BAM_FUNMAP) != 0
    property mate_is_unmapped: 
        """true if the mate is unmapped""" 
        def __get__(self): return (self.flag & BAM_FMUNMAP) != 0
    property is_reverse:
        """true if read is mapped to reverse strand"""
        def __get__(self):return (self.flag & BAM_FREVERSE) != 0
    property mate_is_reverse:
        """true is read is mapped to reverse strand"""
        def __get__(self): return (self.flag & BAM_FMREVERSE) != 0
    property is_read1: 
        """true if this is read1"""
        def __get__(self): return (self.flag & BAM_FREAD1) != 0
    property is_read2:
        """true if this is read2"""
        def __get__(self): return (self.flag & BAM_FREAD2) != 0
    property is_secondary:
        """true if not primary alignment"""
        def __get__(self): return (self.flag & BAM_FSECONDARY) != 0
    property is_qcfail:
        """true if QC failure"""
        def __get__(self): return (self.flag & BAM_FSECONDARY) != 0
    property is_duplicate:
        """ true if optical or PCR duplicate"""
        def __get__(self): return (self.flag & BAM_FDUP) != 0
    
    def opt(self, tag):
        """retrieves optional data given a two-letter *tag*"""
        #see bam_aux.c: bam_aux_get() and bam_aux2i() etc 
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
    
    #    def fancy_str (self):
    #        """
    #        returns list of fieldnames/values in pretty format for debugging
    #        """
    #        ret_string = []
    #        field_names = {
    #                        "tid":           "Contig index",
    #                        "pos":           "Mapped position on contig",
    #
    #                        "mtid":          "Contig index for mate pair",
    #                        "mpos":          "Position of mate pair",
    #                        "isize":         "Insert size",
    #
    #                        "flag":          "Binary flag",
    #                        "n_cigar":       "Count of cigar entries",
    #                        "cigar":         "Cigar entries",
    #                        "qual":          "Mapping quality",
    #
    #                        "bin":           "Bam index bin number",
    #
    #                        "l_qname":       "Length of query name",
    #                        "qname":         "Query name",
    #
    #                        "l_qseq":        "Length of query sequence",
    #                        "qseq":          "Query sequence",
    #                        "bqual":         "Quality scores",
    #
    #
    #                        "l_aux":         "Length of auxilary data",
    #                        "m_data":        "Maximum data length",
    #                        "data_len":      "Current data length",
    #                        }
    #        fields_names_in_order = ["tid", "pos", "mtid", "mpos", "isize", "flag", 
    #                                 "n_cigar", "cigar", "qual", "bin", "l_qname", "qname", 
    #                                 "l_qseq", "qseq", "bqual", "l_aux", "m_data", "data_len"]
    #
    #        for f in fields_names_in_order:
    #            if not f in self.__dict__:
    #                continue
    #            ret_string.append("%-30s %-10s= %s" % (field_names[f], "(" + f + ")", self.__getattribute__(f)))
    #
    #        for f in self.__dict__:
    #            if not f in field_names:
    #                ret_string.append("%-30s %-10s= %s" % (f, "", self.__getattribute__(f)))
#        return ret_string

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

cdef class PileupRead:
    '''A read aligned to a column.
    '''

    cdef:
         AlignedRead _alignment
         int32_t  _qpos
         int _indel
         int _level
         uint32_t _is_del
         uint32_t _is_head
         uint32_t _is_tail

    def __cinit__( self ):
        pass

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
    property head:
        def __get__(self):
            return self._is_head
    property tail:
        def __get__(self):
            return self._is_tail
    
        
        
    
            
            
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
        sys.stdout.flush()          #  Buffered data goes to old stream.
        os.dup2(fd, self.id)        #  Open unit 1 on new stream.
        os.close(fd)                #  Close other unit (look out, caller.)
            
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

    # patch for `samtools view`
    # samtools `view` closes stdout, from which I can not
    # recover. Thus redirect output to file with -o option.
    if method == "view":
        if "-o" in args: raise ValueError("option -o is forbidden in samtools view")
        args = ( "-o", stdout_f ) + args

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

               

