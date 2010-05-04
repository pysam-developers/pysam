# cython: embedsignature=True
# adds doc-strings for sphinx

import tempfile, os, sys, types, itertools, struct, ctypes

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
    # destroy dummy delegate created in constructor
    # to prevent memory leak.
    pysam_bam_destroy1(dest._delegate)
    dest._delegate = bam_dup1(src)
    return dest

cdef class PileupProxy
cdef makePileupProxy( bam_plbuf_t * buf, int n ):
     cdef PileupProxy dest
     dest = PileupProxy()
     dest.buf = buf
     dest.n = n
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
            "\n" + "\n".join( map(str, self.pileups) )

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

class StderrStore():
    '''
    stderr is captured. 
    '''
    def __init__(self):
        self.stderr_h, self.stderr_f = tempfile.mkstemp()
        self.stderr_save = Outs( sys.stderr.fileno() )
        self.stderr_save.setfd( self.stderr_h )
        
    def release(self):
        self.stderr_save.restore()
        if os.path.exists(self.stderr_f):
            os.remove( self.stderr_f )

    def __del__(self):
        self.release()

######################################################################
######################################################################
######################################################################
# valid types for sam headers
VALID_HEADER_TYPES = { "HD" : dict, 
                       "SQ" : list, 
                       "RG" : list, 
                       "PG" : list, 
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
    
    *mode* should be ``r`` for reading or ``w`` for writing. The default is text mode so for binary 
    (:term:`BAM`) I/O you should append ``b`` for compressed or ``u`` for uncompressed :term:`BAM` output. 
    Use ``h`` to output header information  in text (:term:`TAM`)  mode.

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

    If an index for a BAM file exists (.bai), it will be opened automatically. Without an index random
    access to reads via :meth:`fetch` and :meth:`pileup` is disabled.
    '''

    cdef char * filename
    # pointer to samfile
    cdef samfile_t * samfile
    # pointer to index
    cdef bam_index_t *index
    # true if file is a bam file
    cdef int isbam

    # current read within iteration
    cdef bam1_t * b

    def __cinit__(self, *args, **kwargs ):
        self.samfile = NULL
        self.isbam = False
        self._open( *args, **kwargs )

        # allocate memory for iterator
        self.b = <bam1_t*>calloc(1, sizeof(bam1_t))

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
        self.samfile = NULL

        cdef bam_header_t * header_to_write
        header_to_write = NULL

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
                assert referencenames and referencelengths, "either supply options `template`, `header` or  both `refernencenames` and `referencelengths` for writing"
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
            store = StderrStore()
            self.samfile = samopen( filename, mode, header_to_write )
            store.release()

            # bam_header_destroy takes care of cleaning up of all the members
            if not template and header_to_write != NULL:
                bam_header_destroy( header_to_write )

        elif mode[0] == "r":
            # open file for reading
            if strncmp( filename, "-", 1) != 0 and not os.path.exists( filename ):
                raise IOError( "file `%s` not found" % filename)

            store = StderrStore()
            self.samfile = samopen( filename, mode, NULL )
            store.release()

        if self.samfile == NULL:
            raise IOError("could not open file `%s`" % filename )

        if mode[0] == "r" and self.isbam:
            if not os.path.exists(filename + ".bai"):
                self.index = NULL
            else:
                # returns NULL if there is no index or index could not be opened
                self.index = bam_index_load(filename)
                if self.index == NULL:
                    raise IOError("error while opening index `%s` " % filename )

    def getrname( self, tid ):
        '''(tid )
        convert numerical :term:`tid` into :ref:`reference` name.'''
        if not 0 <= tid < self.samfile.header.n_targets:
            raise ValueError( "tid out of range 0<=tid<%i" % self.samfile.header.n_targets )
        return self.samfile.header.target_name[tid]

    def _parseRegion( self, 
                      reference = None, 
                      start = None, 
                      end = None, 
                      region = None ):
        '''parse region information.

        raise Value for for invalid regions.

        returns a tuple of region, tid, start and end. Region
        is a valid samtools :term:`region` or None if the region
        extends over the whole file.

        Note that regions are 1-based, while start,end are python coordinates.
        '''
        
        cdef int rtid
        cdef int rstart
        cdef int rend
        cdef int max_pos
        max_pos = 2 << 29

        rtid = rstart = rend = 0

        # translate to a region
        if reference:
            if start != None and end != None:
                region = "%s:%i-%i" % (reference, start+1, end)
            else:
                region = reference

        if region:
            store = StderrStore()
            bam_parse_region( self.samfile.header, region, &rtid, &rstart, &rend)        
            store.release()
            if rtid < 0: raise ValueError( "invalid region `%s`" % region )
            if rstart > rend: raise ValueError( 'invalid region: start (%i) > end (%i)' % (rstart, rend) )
            if not 0 <= rstart < max_pos: raise ValueError( 'start out of range (%i)' % rstart )
            if not 0 <= rend < max_pos: raise ValueError( 'end out of range (%i)' % rend )

        return region, rtid, rstart, rend

    def fetch( self, 
               reference = None, 
               start = None, 
               end = None, 
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

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        region, rtid, rstart, rend = self._parseRegion( reference, start, end, region )

        if self.isbam:
            if callback:
                if not region:
                    raise ValueError( "callback functionality requires a region/reference" )
                if not self._hasIndex(): raise ValueError( "no index available for fetch" )
                return bam_fetch(self.samfile.x.bam, 
                                 self.index, rtid, rstart, rend, <void*>callback, fetch_callback )
            else:
                if region:
                    return IteratorRow( self, rtid, rstart, rend )
                else:
                    if until_eof:
                        return IteratorRowAll( self )
                    else:
                        # return all targets by chaining the individual targets together.
                        if not self._hasIndex(): raise ValueError( "no index available for fetch" )
                        i = []
                        rstart = 0
                        rend = 1<<29
                        for rtid from 0 <= rtid < self.nreferences: 
                            i.append( IteratorRow( self, rtid, rstart, rend))
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

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        region, rtid, rstart, rend = self._parseRegion( reference, start, end, region )
        
        if self.isbam:
            if not self._hasIndex(): raise ValueError( "no index available for pileup" )

            if callback:
                if not region:
                    raise ValueError( "callback functionality requires a region/reference" )

                buf = bam_plbuf_init( <bam_pileup_f>pileup_callback, <void*>callback )
                bam_fetch(self.samfile.x.bam, 
                          self.index, rtid, rstart, rend, 
                          buf, pileup_fetch_callback )
                
                # finalize pileup
                bam_plbuf_push( NULL, buf)
                bam_plbuf_destroy(buf)
            else:
                if region:
                    return IteratorColumn( self, rtid, rstart, rend )
                else:
                    # return all targets by chaining the individual targets together.
                    i = []
                    rstart = 0
                    rend = 1<<29
                    for rtid from 0 <= rtid < self.nreferences: 
                        i.append( IteratorColumn( self, rtid, rstart, rend))
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
        # remember: dealloc cannot call other methods
        # Note that __del__ is not called.
        self.close()
        pysam_bam_destroy1(self.b)

    def write( self, AlignedRead read ):
        '''(AlignedRead read )
        write a single :class:`pysam.AlignedRead`..

        return the number of bytes written.
        '''
        return samwrite( self.samfile, read._delegate )

    property nreferences:
        '''number of :term:`reference` sequences in the file.'''
        def __get__(self):
            return self.samfile.header.n_targets

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
                        key, value = field.split(":",1)
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

    def __iter__(self):
        return self 

    cdef bam1_t * getCurrent( self ):
        return self.b

    cdef int cnext(self):
        '''cversion of iterator. Used by IteratorColumn'''
        cdef int ret
        return samread(self.samfile, self.b)

    def __next__(self): 
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        cdef int ret
        ret = samread(self.samfile, self.b)
        if (ret > 0):
            return makeAlignedRead( self.b )
        else:
            raise StopIteration

cdef class Fastafile:
    '''*(filename)*
              
    A *FASTA* file. The file is automatically opened.

    The file expects an indexed fasta file.

    TODO: 
        add automatic indexing.
        add function to get sequence names.
    '''

    cdef char * filename
    # pointer to fastafile
    cdef faidx_t * fastafile

    def __cinit__(self, *args, **kwargs ):
        self.fastafile = NULL
        self._open( *args, **kwargs )

    def _isOpen( self ):
        '''return true if samfile has been opened.'''
        return self.fastafile != NULL

    def _open( self, 
               char * filename ):
        '''open an indexed fasta file.

        This method expects an indexed fasta file.
        '''

        # close a previously opened file
        if self.fastafile != NULL: self.close()
        self.filename = filename
        self.fastafile = fai_load( filename )

        if self.fastafile == NULL:
            raise IOError("could not open file `%s`" % filename )

    def close( self ):
        if self.fastafile != NULL:
            fai_destroy( self.fastafile )
            self.fastafile = NULL

    def fetch( self, 
               reference = None, 
               start = None, 
               end = None,
               region = None):
               
        '''*(reference = None, start = None, end = None, region = None)*
               
        fetch :meth:`AlignedRead` objects in a :term:`region` using 0-based indexing. The region is specified by
        :term:`reference`, *start* and *end*. Alternatively, a samtools :term:`region` string can be supplied.
        '''
        
        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        cdef int len, max_pos
        cdef char * seq
        max_pos = 2 << 29

        if not region:
            if reference == None: raise ValueError( 'no sequence/region supplied.' )
            if start == None and end == None:
                region = "%s" % str(reference)
            elif start == None or end == None:
                raise ValueError( 'only start or only end of region supplied' )
            else:
                if start > end: raise ValueError( 'invalid region: start (%i) > end (%i)' % (start, end) )
		# valid ranges are from 0 to 2^29-1
                if not 0 <= start < max_pos: raise ValueError( 'start out of range (%i)' % start )
                if not 0 <= end < max_pos: raise ValueError( 'end out of range (%i)' % end )
                region = "%s:%i-%i" % (reference, start+1, end )

        # samtools adds a '\0' at the end
        seq = fai_fetch( self.fastafile, region, &len )
        # copy to python
        result = seq
        # clean up
        free(seq)
        
        return result

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
    cdef Samfile                samfile
    def __cinit__(self, Samfile samfile, int tid, int beg, int end ):
        self.bam_iter = NULL

        assert samfile._isOpen()
        assert samfile._hasIndex()
        
        # makes sure that samfile stays alive as long as the
        # iterator is alive.
        self.samfile = samfile

        # parse the region
        self.error_state = 0
        self.error_msg = None

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
        pysam_bam_destroy1(self.b)
        
cdef class IteratorColumn:
    '''iterates over columns.

    This iterator wraps the pileup functionality of samtools.
    
    For reasons of efficiency, the iterator returns the current 
    pileup buffer. As this buffer is updated at every iteration, 
    the contents of this iterator will change accordingly. Hence the conversion to
    a list will not produce the expected result::
    
       f = Samfile("file.bam", "rb")
       result = list( f.pileup() )

    Here, result will contain ``n`` objects of type :class:`PileupProxy` for ``n`` columns, 
    but each object will contain the same information.
    
    If the results of several columns are required at the same time, the results
    need to be stored explicitely::

       result = [ x.pileups() for x in f.pileup() ]

    Here, result will be a list of ``n`` lists of objects of type :class:`PileupRead`.

    '''
    cdef bam_plbuf_t *buf

    # check if first iteration
    cdef int notfirst
    # result of the last plbuf_push
    cdef int n_pu
    cdef int eof 
    cdef IteratorRow iter

    def __cinit__(self, Samfile samfile, int tid, int start, int end ):

        self.iter = IteratorRow( samfile, tid, start, end )
        self.buf = bam_plbuf_init(NULL, NULL )
        self.n_pu = 0
        self.eof = 0

    def __iter__(self):
        return self 

    cdef int cnext(self):
        '''perform next iteration.
        
        return 1 if there is a buffer to emit. Return 0 for end of iteration.
        '''

        cdef int retval1, retval2

        # pysam bam_plbuf_push returns:
        # 1: if buf is full and can be emitted
        # 0: if b has been added
        # -1: if there was an error

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
                return self.n_pu

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
            return makePileupProxy( self.buf, self.n_pu )
        else:
            raise StopIteration

    def __dealloc__(self):
        bam_plbuf_destroy(self.buf);

cdef class AlignedRead:
    '''
    Class representing an aligned read. see SAM format specification for meaning of fields (http://samtools.sourceforge.net/).

    This class stores a handle to the samtools C-structure representing
    an aligned read. Member read access is forwarded to the C-structure
    and converted into python objects. This implementation should be fast,
    as only the data needed is converted.

    For write access, the C-structure is updated in-place. This is
    not the most efficient way to build BAM entries, as the variable
    length data is concatenated and thus needs to resized if
    a field is updated. Furthermore, the BAM entry might be
    in an inconsistent state. The :meth:`~validate` method can
    be used to check if an entry is consistent.

    One issue to look out for is that the sequence should always
    be set *before* the quality scores. Setting the sequence will
    also erase any quality scores that were set previously.
    '''
    cdef:
         bam1_t * _delegate 

    def __cinit__( self ):
        # see bam_init1
        self._delegate = <bam1_t*>calloc( 1, sizeof( bam1_t) )
        # allocate some memory 
        # If size is 0, calloc does not return a pointer that can be passed to free()
        # so allocate 40 bytes for a new read
        self._delegate.m_data = 40
        self._delegate.data = <uint8_t *>calloc( self._delegate.m_data, 1 )
        self._delegate.data_len = 0

    def __dealloc__(self):
        '''clear up memory.'''
        pysam_bam_destroy1(self._delegate)
    
    def __str__(self):
        """todo"""
        return "\t".join(map(str, (self.qname,
                                   self.rname,
                                   self.pos,
                                   self.cigar,
                                   self.qual,
                                   self.flag,
                                   self.seq,
                                   self.mapq,
                                   self.tags)))
    
       
    def __cmp__(self, AlignedRead other):
        '''return true, if contents in this are binary equal to ``other``.'''
        cdef int retval, x
        cdef bam1_t *t, *o
        t = self._delegate
        o = other._delegate

        # uncomment for debugging purposes
        # cdef unsigned char * oo, * tt
        # tt = <unsigned char*>(&t.core)
        # oo = <unsigned char*>(&o.core)
        # for x from 0 <= x < sizeof( bam1_core_t): print x, tt[x], oo[x]
        # tt = <unsigned char*>(t.data)
        # oo = <unsigned char*>(o.data)
        # for x from 0 <= x < max(t.data_len, o.data_len): print x, tt[x], oo[x], chr(tt[x]), chr(oo[x])

        retval = memcmp( &t.core, 
                          &o.core, 
                          sizeof( bam1_core_t ))

        if retval: return retval
        retval = cmp( t.data_len, o.data_len)
        if retval: return retval
        return memcmp( t.data, 
                       o.data, 
                       sizeof( t.data_len ))

    property qname:
        """the query name (None if not present)"""
        def __get__(self):
            cdef bam1_t * src 
            src = self._delegate
            if src.core.l_qname == 0: return None
            return <char *>pysam_bam1_qname( src )

        def __set__(self, qname ):
            if qname == None or len(qname) == 0: return
            cdef bam1_t * src 
            cdef int l 
            cdef char * p

            src = self._delegate            
            p = pysam_bam1_qname( src )

            # the qname is \0 terminated
            l = len(qname) + 1
            pysam_bam_update( src, 
                              src.core.l_qname, 
                              l, 
                              <uint8_t*>p )

            src.core.l_qname = l

            # re-acquire pointer to location in memory
            # as it might have moved
            p = pysam_bam1_qname(src)

            strncpy( p, qname, l )
            
    property cigar:
        """the :term:`cigar` alignment (None if not present).
        """
        def __get__(self):
            cdef uint32_t * cigar_p
            cdef bam1_t * src 
            cdef op, l, cigar
            src = self._delegate
            if src.core.n_cigar == 0: return None
            
            cigar = []
            cigar_p = pysam_bam1_cigar(src);
            for k from 0 <= k < src.core.n_cigar:
                op = cigar_p[k] & BAM_CIGAR_MASK
                l = cigar_p[k] >> BAM_CIGAR_SHIFT
                cigar.append((op, l))
            return cigar

        def __set__(self, values ):
            if values == None or len(values) == 0: return
            cdef uint32_t * p
            cdef bam1_t * src 
            cdef op, l
            cdef int k

            k = 0

            src = self._delegate

            # get location of cigar string
            p = pysam_bam1_cigar(src)

            # create space for cigar data within src.data
            pysam_bam_update( src, 
                              src.core.n_cigar * 4,
                              len(values) * 4, 
                              p )
            
            # length is number of cigar operations, not bytes
            src.core.n_cigar = len(values)

            # re-acquire pointer to location in memory
            # as it might have moved
            p = pysam_bam1_cigar(src)

            # insert cigar operations
            for op, l in values:
                p[k] = l << BAM_CIGAR_SHIFT | op
                k += 1

            ## setting the cigar string also updates the "bin" attribute
            src.core.bin = bam_reg2bin( src.core.pos, bam_calend( &src.core, p))

    property seq:
        """the query sequence (None if not present)"""
        def __get__(self):
            cdef bam1_t * src
            cdef uint8_t * p 
            cdef char * s
            src = self._delegate
            bam_nt16_rev_table = "=ACMGRSVTWYHKDBN"
            ## parse qseq (bam1_seq)
            if src.core.l_qseq == 0: return None

            s = < char *> calloc(src.core.l_qseq + 1 , sizeof(char))
            p = pysam_bam1_seq( src )
            for k from 0 <= k < src.core.l_qseq:
            ## equivalent to bam_nt16_rev_table[bam1_seqi(s, i)] (see bam.c)
                s[k] = "=ACMGRSVTWYHKDBN"[((p)[(k) / 2] >> 4 * (1 - (k) % 2) & 0xf)]
            retval=s
            free(s)
            return retval

        def __set__(self,seq):
            # samtools manages sequence and quality length memory together
            # if no quality information is present, the first byte says 0xff.
            
            if seq == None or len(seq) == 0: return
            cdef bam1_t * src
            cdef uint8_t * p 
            cdef char * s
            src = self._delegate
            cdef int l, k, nbytes_new, nbytes_old

            l = len(seq)
            
            # as the sequence is stored in half-bytes, the total length (sequence
            # plus quality scores) is (l+1)/2 + l
            nbytes_new = (l+1)/2 + l
            nbytes_old = (src.core.l_qseq+1)/2 + src.core.l_qseq
            # acquire pointer to location in memory
            p = pysam_bam1_seq( src )
            src.core.l_qseq = l

            pysam_bam_update( src, 
                              nbytes_old,
                              nbytes_new,
                              p)
            # re-acquire pointer to location in memory
            # as it might have moved
            p = pysam_bam1_seq( src )
            for k from 0 <= k < nbytes_new: p[k] = 0
            # convert to C string
            s = seq
            for k from 0 <= k < l:
                p[k/2] |= pysam_translate_sequence(s[k]) << 4 * (1 - k % 2)

            # erase qualities
            p = pysam_bam1_qual( src )
            p[0] = 0xff

    property qual:
        """the base quality (None if not present)"""
        def __get__(self):
            cdef bam1_t * src 
            cdef uint8_t * p
            cdef char * q
            src = self._delegate
            if src.core.l_qseq == 0: return None

            p = pysam_bam1_qual( src )
            if p[0] == 0xff: return None

            q = < char *>calloc(src.core.l_qseq + 1 , sizeof(char))
            for k from 0 <= k < src.core.l_qseq:
            ## equivalent to t[i] + 33 (see bam.c)
                q[k] = p[k] + 33
            # convert to python string
            retval=q
            # clean up
            free(q)
            return retval

        def __set__(self,qual):
            # note that space is already allocated via the sequences
            cdef bam1_t * src
            cdef uint8_t * p
            cdef char * q 
            src = self._delegate
            p = pysam_bam1_qual( src )
            if qual == None or len(qual) == 0:
                # if absent - set to 0xff
                p[0] = 0xff
                return
            cdef int l
            # convert to C string
            q = qual
            l = len(qual)
            if src.core.l_qseq != l:
                raise ValueError("quality and sequence mismatch: %i != %i" % (l, src.core.l_qseq))
            assert src.core.l_qseq == l
            for k from 0 <= k < l:
                p[k] = <uint8_t>q[k] - 33

    property tags:
        """the tags in the AUX field."""
        def __get__(self):
            cdef char * ctag
            cdef bam1_t * src
            cdef uint8_t * s
            cdef char tpe
            
            src = self._delegate
            if src.l_aux == 0: return None
            
            s = pysam_bam1_aux( src )
            result = []
            ctag = <char*>calloc( 3, sizeof(char) )
            cdef int x
            while s < (src.data + src.data_len):
                # get tag
                ctag[0] = s[0]
                ctag[1] = s[1]
                pytag = ctag

                s += 2

                # convert type - is there a better way?
                ctag[0] = s[0]
                ctag[1] = 0
                pytype = ctag

                # get type and value 
                # how do I do char literal comparison in cython?
                # the code below works (i.e, is C comparison)
                tpe = toupper(s[0])
                if tpe == 'S'[0]:
                    value = <int>bam_aux2i(s)            
                    s += 2
                elif tpe == 'I'[0]:
                    value = <int>bam_aux2i(s)            
                    s += 4
                elif tpe == 'F'[0]:
                    value = <float>bam_aux2f(s)
                    s += 4
                elif tpe == 'D'[0]:
                    value = <double>bam_aux2d(s)
                    s += 8
                elif tpe == 'C'[0]:
                    value = <int>bam_aux2i(s)
                    s += 1
                elif tpe == 'A'[0]:
                    # there might a more efficient way
                    # to convert a char into a string
                    value = "%c" % <char>bam_aux2A(s)
                    s += 1
                elif tpe == 'Z'[0]:
                    value = <char*>bam_aux2Z(s)
                    # +1 for NULL terminated string
                    s += len(value) + 1

                # skip over type
                s += 1

                # ignore pytype
                result.append( (pytag, value) )

            free( ctag )
            return result

        def __set__(self, tags):
            cdef char * ctag
            cdef bam1_t * src
            cdef uint8_t * s
            cdef uint8_t * new_data
            cdef int guessed_size, control_size
            src = self._delegate
            cdef int max_size, size
            max_size = 4000

            # map samtools code to python.struct code and byte size
            buffer = ctypes.create_string_buffer(max_size)

            offset = 0
            for pytag, value in tags:
                t = type(value)
                if t == types.FloatType:
                    fmt = "<cccf"
                elif t == types.IntType:
                    if value < 0:
                        if value >= -127: fmt, pytype = "<cccb", 'c'
                        elif value >= -32767: fmt, pytype = "<ccch", 's'
                        elif value < -2147483648: raise ValueError( "integer %i out of range of BAM/SAM specification" % value )
                        else: fmt, ctype = "<ccci", 'i'[0]
                    else:
                        if value <= 255: fmt, pytype = "<cccB", 'C'
                        elif value <= 65535: fmt, pytype = "<cccH", 'S'
                        elif value > 4294967295: raise ValueError( "integer %i out of range of BAM/SAM specification" % value )
                        else: fmt, pytype = "<cccI", 'I'
                else:
                    # Note: hex strings (H) are not supported yet
                    if len(value) == 1:
                        fmt, pytype = "<cccc", 'A'
                    else:
                        fmt, pytype = "<ccc%is" % (len(value)+1), 'Z'

                size = struct.calcsize(fmt)
                if offset + size > max_size:
                    raise NotImplementedError("tags field too large")

                struct.pack_into( fmt,
                                  buffer,
                                  offset,
                                  pytag[0],
                                  pytag[1],
                                  pytype,
                                  value )
                offset += size
            
            # delete the old data and allocate new
            pysam_bam_update( src, 
                              src.l_aux,
                              offset,
                              pysam_bam1_aux( src ) )
            
            src.l_aux = offset

            if offset == 0: return

            # get location of new data
            s = pysam_bam1_aux( src )            
            
            # check if there is direct path from buffer.raw to tmp
            cdef char * temp 
            temp = buffer.raw
            memcpy( s, temp, offset )            

    property flag: 
        """properties flag"""
        def __get__(self): return self._delegate.core.flag
        def __set__(self, flag): self._delegate.core.flag = flag
    property rname: 
        """:term:`reference` ID"""
        def __get__(self): return self._delegate.core.tid
        def __set__(self, tid): self._delegate.core.tid = tid
    property pos: 
        """0-based leftmost coordinate"""
        def __get__(self): return self._delegate.core.pos
        def __set__(self, pos): 
            ## setting the cigar string also updates the "bin" attribute
            cdef bam1_t * src
            src = self._delegate
            if src.core.n_cigar:
                src.core.bin = bam_reg2bin( src.core.pos, bam_calend( &src.core, pysam_bam1_cigar(src)) )
            else:
                src.core.bin = bam_reg2bin( src.core.pos, src.core.pos + 1)
            self._delegate.core.pos = pos
    property bin: 
        """properties bin"""
        def __get__(self): return self._delegate.core.bin
        def __set__(self, bin): self._delegate.core.bin = bin
    property rlen:
        '''length of the read (read only). Returns 0 if not given.'''
        def __get__(self): return self._delegate.core.l_qseq
    property mapq: 
        """mapping quality"""
        def __get__(self): return self._delegate.core.qual
        def __set__(self, qual): self._delegate.core.qual = qual
    property mrnm:
        """the :term:`reference` id of the mate """     
        def __get__(self): return self._delegate.core.mtid
        def __set__(self, mtid): self._delegate.core.mtid = mtid
    property mpos: 
        """the position of the mate"""
        def __get__(self): return self._delegate.core.mpos
        def __set__(self, mpos): self._delegate.core.mpos = mpos
    property isize: 
        """the insert size"""
        def __get__(self): return self._delegate.core.isize
        def __set__(self, isize): self._delegate.core.isize = isize
    property is_paired: 
        """true if read is paired in sequencing"""
        def __get__(self): return (self._delegate.core.flag & BAM_FPAIRED) != 0
        def __set__(self,val): 
            if val: self._delegate.core.flag |= BAM_FPAIRED
            else: self._delegate.core.flag &= ~BAM_FPAIRED
    property is_proper_pair:
        """true if read is mapped in a proper pair"""
        def __get__(self): return (self.flag & BAM_FPROPER_PAIR) != 0
        def __set__(self,val): 
            if val: self._delegate.core.flag |= BAM_FPROPER_PAIR
            else: self._delegate.core.flag &= ~BAM_FPROPER_PAIR
    property is_unmapped:
        """true if read itself is unmapped"""
        def __get__(self): return (self.flag & BAM_FUNMAP) != 0
        def __set__(self,val): 
            if val: self._delegate.core.flag |= BAM_FUNMAP
            else: self._delegate.core.flag &= ~BAM_FUNMAP
    property mate_is_unmapped: 
        """true if the mate is unmapped""" 
        def __get__(self): return (self.flag & BAM_FMUNMAP) != 0
        def __set__(self,val): 
            if val: self._delegate.core.flag |= BAM_FMUNMAP
            else: self._delegate.core.flag &= ~BAM_FMUNMAP
    property is_reverse:
        """true if read is mapped to reverse strand"""
        def __get__(self):return (self.flag & BAM_FREVERSE) != 0
        def __set__(self,val): 
            if val: self._delegate.core.flag |= BAM_FREVERSE
            else: self._delegate.core.flag &= ~BAM_FREVERSE
    property mate_is_reverse:
        """true is read is mapped to reverse strand"""
        def __get__(self): return (self.flag & BAM_FMREVERSE) != 0
        def __set__(self,val): 
            if val: self._delegate.core.flag |= BAM_FMREVERSE
            else: self._delegate.core.flag &= ~BAM_FMREVERSE
    property is_read1: 
        """true if this is read1"""
        def __get__(self): return (self.flag & BAM_FREAD1) != 0
        def __set__(self,val): 
            if val: self._delegate.core.flag |= BAM_FREAD1
            else: self._delegate.core.flag &= ~BAM_FREAD1
    property is_read2:
        """true if this is read2"""
        def __get__(self): return (self.flag & BAM_FREAD2) != 0
        def __set__(self,val): 
            if val: self._delegate.core.flag |= BAM_FREAD2
            else: self._delegate.core.flag &= ~BAM_FREAD2
    property is_secondary:
        """true if not primary alignment"""
        def __get__(self): return (self.flag & BAM_FSECONDARY) != 0
        def __set__(self,val): 
            if val: self._delegate.core.flag |= BAM_FSECONDARY
            else: self._delegate.core.flag &= ~BAM_FSECONDARY
    property is_qcfail:
        """true if QC failure"""
        def __get__(self): return (self.flag & BAM_FQCFAIL) != 0
        def __set__(self,val): 
            if val: self._delegate.core.flag |= BAM_FQCFAIL
            else: self._delegate.core.flag &= ~BAM_FQCFAIL
    property is_duplicate:
        """ true if optical or PCR duplicate"""
        def __get__(self): return (self.flag & BAM_FDUP) != 0
        def __set__(self,val): 
            if val: self._delegate.core.flag |= BAM_FDUP
            else: self._delegate.core.flag &= ~BAM_FDUP
    
    def opt(self, tag):
        """retrieves optional data given a two-letter *tag*"""
        #see bam_aux.c: bam_aux_get() and bam_aux2i() etc 
        cdef uint8_t * v
        v = bam_aux_get(self._delegate, tag)
        if v == NULL: raise KeyError( "tag '%s' not present" % tag )
        type = chr(v[0])
        if type == 'c' or type == 'C' or type == 's' or type == 'S' or type == 'i':
            return <int>bam_aux2i(v)            
        elif type == 'f':
            return <float>bam_aux2f(v)
        elif type == 'd':
            return <double>bam_aux2d(v)
        elif type == 'A':
            # there might a more efficient way
            # to convert a char into a string
            return '%c' % <char>bam_aux2A(v)
        elif type == 'Z':
            return <char*>bam_aux2Z(v)
    
    def fancy_str (self):
        """returns list of fieldnames/values in pretty format for debugging
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

cdef class PileupProxy:
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

    This class is a proxy for results returned by the samtools pileup engine.
    If the underlying engine iterator advances, the results of this column
    will change.
    '''
    cdef bam_plbuf_t * buf
    cdef int n_pu

    def __cinit__(self ):
        pass

    def __str__(self):
        return "\t".join( map(str, (self.tid, self.pos, self.n))) +\
            "\n" +\
            "\n".join( map(str, self.pileups) )

    property tid:
        '''the chromosome ID as is defined in the header'''
        def __get__(self): return pysam_get_tid( self.buf )

    property n:
        '''number of reads mapping to this column.'''
        def __get__(self): return self.n_pu
        def __set__(self, n): self.n_pu = n

    property pos:
        def __get__(self): return pysam_get_pos( self.buf )

    property pileups:
        '''list of reads (:class:`pysam.PileupRead`) aligned to this column'''
        def __get__(self):
            cdef bam_pileup1_t * pl
            pl = pysam_get_pileup( self.buf )
            pileups = []
            # warning: there could be problems if self.n and self.buf are
            # out of sync.
            for x from 0 <= x < self.n_pu:
                pileups.append( makePileupRead( &pl[x]) )
            return pileups

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
    property is_head:
        def __get__(self):
            return self._is_head
    property is_tail:
        def __get__(self):
            return self._is_tail
    property level:
        def __get__(self):
            return self._level

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

__all__ = ["Samfile", 
           "Fastafile",
           "IteratorRow", 
           "IteratorRowAll", 
           "IteratorColumn", 
           "AlignedRead", 
           "PileupColumn", 
           "PileupProxy", 
           "PileupRead" ]

               

