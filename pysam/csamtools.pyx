# cython: embedsignature=True
# cython: profile=True
# adds doc-strings for sphinx

import tempfile, os, sys, types, struct, ctypes, collections, re

from python_string cimport PyString_FromStringAndSize, PyString_AS_STRING
from python_exc    cimport PyErr_SetString

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

DEF BAM_CMATCH     = 0
DEF BAM_CINS       = 1
DEF BAM_CDEL       = 2
DEF BAM_CREF_SKIP  = 3
DEF BAM_CSOFT_CLIP = 4
DEF BAM_CHARD_CLIP = 5
DEF BAM_CPAD       = 6

#####################################################################
## hard-coded constants
cdef char * bam_nt16_rev_table = "=ACMGRSVTWYHKDBN"
cdef int max_pos = 2 << 29

#####################################################################
#####################################################################
#####################################################################
## private factory methods
#####################################################################
cdef class AlignedRead
cdef makeAlignedRead(bam1_t * src):
    '''enter src into AlignedRead.'''
    cdef AlignedRead dest = AlignedRead.__new__(AlignedRead)
    dest._delegate = bam_dup1(src)
    return dest

cdef class PileupProxy
cdef makePileupProxy( bam_pileup1_t * plp, int tid, int pos, int n ):
     cdef PileupProxy dest = PileupProxy.__new__(PileupProxy)
     dest.plp = plp
     dest.tid = tid
     dest.pos = pos
     dest.n = n
     return dest

cdef class PileupRead
cdef makePileupRead( bam_pileup1_t * src ):
    '''fill a  PileupRead object from a bam_pileup1_t * object.'''
    cdef PileupRead dest = PileupRead.__new__(PileupRead)
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

    cdef int x
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
        
    def readAndRelease( self ):
        self.stderr_save.restore()
        lines = []
        if os.path.exists(self.stderr_f):
            lines = open( self.stderr_f, "r" ).readlines()
            os.remove( self.stderr_f )
        return lines

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
    # true if file is not on the local filesystem
    cdef int isremote
    # current read within iteration
    cdef bam1_t * b
    # file opening mode
    cdef char * mode

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
               mode = 'r',
               Samfile template = None,
               referencenames = None,
               referencelengths = None,
               text = None,
               header = None,
               port = None,
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

        self.isremote = strncmp(filename,"http:",5) == 0 or \
            strncmp(filename,"ftp:",4) == 0 

        cdef char * ctext
        ctext = NULL

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
                assert referencenames and referencelengths, "either supply options `template`, `header` or  both `referencenames` and `referencelengths` for writing"
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

                if text != None:
                    # copy without \0
                    ctext = text
                    header_to_write.l_text = strlen(ctext)
                    header_to_write.text = <char*>calloc( strlen(ctext), sizeof(char) )
                    memcpy( header_to_write.text, ctext, strlen(ctext) )

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
            if strncmp( filename, "-", 1) != 0 and \
                    not self.isremote and \
                    not os.path.exists( filename ):
                raise IOError( "file `%s` not found" % filename)

            store = StderrStore()
            self.samfile = samopen( filename, mode, NULL )
            result = store.readAndRelease()
            # test for specific messages as open also outputs status messages
            # that can be ignored.
            if "[bam_header_read] invalid BAM binary header (this is not a BAM file).\n" in result:
                raise ValueError( "invalid BAM binary header (is this a BAM file?)" )
            elif '[samopen] no @SQ lines in the header.\n' in result:
                raise ValueError( "no @SQ lines in the header (is this a SAM file?)")

        if self.samfile == NULL:
            raise IOError("could not open file `%s`" % filename )

        # check for index and open if present
        if mode[0] == "r" and self.isbam:

            if not self.isremote:
                if not os.path.exists(filename +".bai"): 
                    self.index = NULL
                else:
                    # returns NULL if there is no index or index could not be opened
                    self.index = bam_index_load(filename)
                    if self.index == NULL:
                        raise IOError("error while opening index `%s` " % filename )
            else:
                self.index = bam_index_load(filename)
                if self.index == NULL:
                    raise IOError("error while opening index `%s` " % filename )
                                    
    def getrname( self, tid ):
        '''(tid )
        convert numerical :term:`tid` into :ref:`reference` name.'''
        if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
        if not 0 <= tid < self.samfile.header.n_targets:
            raise ValueError( "tid out of range 0<=tid<%i" % self.samfile.header.n_targets )
        return self.samfile.header.target_name[tid]

    def gettid( self, reference ):
        '''(reference)
        convert :ref:`reference` name into numerical :term:`tid`

        returns -1 if reference is not known.
        '''
        if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
        return pysam_reference2tid( self.samfile.header, reference )

    def _parseRegion( self, 
                      reference = None, 
                      start = None, 
                      end = None,
                      region = None ):
        '''parse region information.

        raise ValueError for for invalid regions.

        returns a tuple of flag, tid, start and end. Flag indicates
        whether some coordinates were supplied.

        Note that regions are 1-based, while start,end are python coordinates.
        '''
        # This method's main objective is to translate from a reference to a tid. 
        # For now, it calls bam_parse_region, which is clumsy. Might be worth
        # implementing it all in pysam (makes use of khash).
        
        cdef int rtid
        cdef int rstart
        cdef int rend

        rtid = -1
        rstart = 0
        rend = max_pos
        if start != None: rstart = start
        if end != None: rend = end

        if region:
            parts = re.split( "[:-]", region )
            reference = parts[0]
            if len(parts) >= 2: rstart = int(parts[1]) - 1
            if len(parts) >= 3: rend = int(parts[2])

        if not reference: return 0, 0, 0, 0

        rtid = self.gettid( reference )
        if rtid < 0: raise ValueError( "invalid reference `%s`" % reference )
        if rstart > rend: raise ValueError( 'invalid coordinates: start (%i) > end (%i)' % (rstart, rend) )
        if not 0 <= rstart < max_pos: raise ValueError( 'start out of range (%i)' % rstart )
        if not 0 <= rend <= max_pos: raise ValueError( 'end out of range (%i)' % rend )

        return 1, rtid, rstart, rend
    
    def seek( self, uint64_t offset, int where = 0):
        '''move to current file to position *offset*'''

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )
        if not self.isbam:
            raise NotImplementedError("seek only available in bam files")
        return bam_seek( self.samfile.x.bam, offset, where )

    def tell( self ):
        '''return current file position'''
        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )
        if not self.isbam:
            raise NotImplementedError("seek only available in bam files")

        return bam_tell( self.samfile.x.bam )

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
        cdef int rtid, rstart, rend, has_coord

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        has_coord, rtid, rstart, rend = self._parseRegion( reference, start, end, region )
        
        if self.isbam:
            if not until_eof and not self._hasIndex() and not self.isremote: 
                raise ValueError( "fetch called on bamfile without index" )

            if callback:
                if not has_coord: raise ValueError( "callback functionality requires a region/reference" )
                if not self._hasIndex(): raise ValueError( "no index available for fetch" )
                return bam_fetch(self.samfile.x.bam, 
                                 self.index, 
                                 rtid, 
                                 rstart, 
                                 rend, 
                                 <void*>callback, 
                                 fetch_callback )
            else:
                if has_coord:
                    return IteratorRowRegion( self, rtid, rstart, rend )
                else:
                    if until_eof:
                        return IteratorRowAll( self )
                    else:
                        return IteratorRowAllRefs(self)
        else:   
            # check if header is present - otherwise sam_read1 aborts
            # this happens if a bamfile is opened with mode 'r'
            if self.samfile.header.n_targets == 0:
                raise ValueError( "fetch called for samfile without header")
                  
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

            *all* reads which overlap the region are returned. The first base returned will be the 
            first base of the first read *not* necessarily the first base of the region used in the query.
        '''
        cdef int rtid, rstart, rend, has_coord
        cdef bam_plbuf_t *buf

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        has_coord, rtid, rstart, rend = self._parseRegion( reference, start, end, region )
        
        if self.isbam:
            if not self._hasIndex(): raise ValueError( "no index available for pileup" )

            if callback:
                if not has_coord: raise ValueError( "callback functionality requires a region/reference" )

                buf = bam_plbuf_init( <bam_pileup_f>pileup_callback, <void*>callback )
                bam_fetch(self.samfile.x.bam, 
                          self.index, rtid, rstart, rend, 
                          buf, pileup_fetch_callback )
                
                # finalize pileup
                bam_plbuf_push( NULL, buf)
                bam_plbuf_destroy(buf)
            else:
                if has_coord:
                    return IteratorColumnRegion( self, rtid, rstart, rend )
                else:
                    return IteratorColumnAllRefs(self)

        else:
            raise NotImplementedError( "pileup of samfiles not implemented yet" )

    def close( self ):
        '''closes file.'''
        if self.samfile != NULL:
            samclose( self.samfile )
            bam_index_destroy(self.index);
            self.samfile = NULL

    def __dealloc__( self ):
        # remember: dealloc cannot call other methods
        # note: no doc string
        # note: __del__ is not called.
        self.close()
        bam_destroy1(self.b)

    def write( self, AlignedRead read ):
        '''(AlignedRead read )
        write a single :class:`pysam.AlignedRead`..

        return the number of bytes written.
        '''
        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        return samwrite( self.samfile, read._delegate )

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    property nreferences:
        '''number of :term:`reference` sequences in the file.'''
        def __get__(self):
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
            return self.samfile.header.n_targets

    property references:
        """tuple with the names of :term:`reference` sequences."""
        def __get__(self): 
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
            t = []
            for x from 0 <= x < self.samfile.header.n_targets:
                t.append( self.samfile.header.target_name[x] )
            return tuple(t)

    property lengths:
        """tuple of the lengths of the :term:`reference` sequences. The lengths are in the same order as :attr:`pysam.Samfile.reference`
        """
        def __get__(self): 
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
            t = []
            for x from 0 <= x < self.samfile.header.n_targets:
                t.append( self.samfile.header.target_len[x] )
            return tuple(t)

    property text:
        '''full contents of the :term:`sam file` header as a string.'''
        def __get__(self):
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
            return PyString_FromStringAndSize(self.samfile.header.text, self.samfile.header.l_text)

    property header:
        '''header information within the :term:`sam file`. The records and fields are returned as 
        a two-level dictionary.
        '''
        def __get__(self):
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )

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
        if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
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

    def __len__(self):
        if self.fastafile == NULL:
            raise ValueError( "calling len() on closed file" )

        return faidx_fetch_nseq(self.fastafile)

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
               
        fetch :meth:`AlignedRead` objects in a :term:`region` using 0-based indexing. 
        
        The region is specified by :term:`reference`, *start* and *end*. 
        
        fetch returns an empty string if the region is out of range or addresses an unknown *reference*.

        If *reference* is given and *start* is None, the sequence from the 
        first base is returned. Similarly, if *end* is None, the sequence 
        until the last base is returned.
        
        Alternatively, a samtools :term:`region` string can be supplied.
        '''
        
        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        cdef int length
        cdef char * seq

        if not region:
            if reference is None: raise ValueError( 'no sequence/region supplied.' )
            if start is None: start = 0
            if end is None: end = max_pos -1

            if start > end: raise ValueError( 'invalid region: start (%i) > end (%i)' % (start, end) )
            if start == end: return ""
            # valid ranges are from 0 to 2^29-1
            if not 0 <= start < max_pos: raise ValueError( 'start out of range (%i)' % start )
            if not 0 <= end < max_pos: raise ValueError( 'end out of range (%i)' % end )
            # note: faidx_fetch_seq has a bug such that out-of-range access
            # always returns the last residue. Hence do not use faidx_fetch_seq,
            # but use fai_fetch instead 
            # seq = faidx_fetch_seq(self.fastafile, 
            #                       reference, 
            #                       start,
            #                       end-1, 
            #                       &length)
            region = "%s:%i-%i" % (reference, start+1, end)
            seq = fai_fetch( self.fastafile, 
                             region,
                             &length )
        else:
            # samtools adds a '\0' at the end
            seq = fai_fetch( self.fastafile, region, &length )

        # copy to python
        if seq == NULL:
            return ""
        else:
            try:
                py_seq = PyString_FromStringAndSize(seq, length)
            finally:
                free(seq)

        return py_seq

    cdef char * _fetch( self, char * reference, int start, int end, int * length ):
        '''fetch sequence for reference, start and end'''
        
        return faidx_fetch_seq(self.fastafile, 
                               reference, 
                               start,
                               end-1, 
                               length )

###########################################################################
###########################################################################
###########################################################################
## turning callbacks elegantly into iterators is an unsolved problem, see the following threads:
## http://groups.google.com/group/comp.lang.python/browse_frm/thread/0ce55373f128aa4e/1d27a78ca6408134?hl=en&pli=1
## http://www.velocityreviews.com/forums/t359277-turning-a-callback-function-into-a-generator.html
## Thus I chose to rewrite the functions requiring callbacks. The downside is that if the samtools C-API or code
## changes, the changes have to be manually entered.
cdef class IteratorRow:
    '''abstract base class for iterators over mapped reads.'''

    pass

cdef class IteratorRowRegion(IteratorRow):
    """iterates over mapped reads in a region.

    The samtools iterators assume that the file
    position between iterations do not change.
    As a consequence, no two iterators can work
    on the same file. To permit this, each iterator
    creates its own file handle by re-opening the
    file.

    Note that the index will be shared between 
    samfile and the iterator.
    """
    
    cdef bam_iter_t             iter # iterator state object
    cdef bam1_t *               b
    cdef int                    retval
    cdef Samfile                samfile
    cdef samfile_t              * fp

    def __cinit__(self, Samfile samfile, int tid, int beg, int end ):

        if not samfile._isOpen():
            raise ValueError( "I/O operation on closed file" )
        
        if not samfile._hasIndex():
            raise ValueError( "no index available for iteration" )
        
        # makes sure that samfile stays alive as long as the
        # iterator is alive
        self.samfile = samfile

        if samfile.isbam: mode = "rb"
        else: mode = "r"

        # reopen the file - note that this makes the iterator
        # slow and causes pileup to slow down significantly.
        store = StderrStore()
        self.fp = samopen( samfile.filename, mode, NULL )
        store.release()

        self.retval = 0

        self.iter = bam_iter_query(self.samfile.index, 
                                   tid, 
                                   beg, 
                                   end)
        self.b = bam_init1()

    def __iter__(self):
        return self 

    cdef bam1_t * getCurrent( self ):
        return self.b

    cdef int cnext(self):
        '''cversion of iterator. Used by IteratorColumn'''
        self.retval = bam_iter_read( self.fp.x.bam, 
                                     self.iter, 
                                     self.b)
        
    def __next__(self): 
        """python version of next().
        """
        self.cnext()
        if self.retval < 0: raise StopIteration
        return makeAlignedRead( self.b )

    def __dealloc__(self):
        bam_destroy1(self.b)
        samclose( self.fp )

cdef class IteratorRowAll(IteratorRow):
    """iterates over all reads
    """

    cdef bam1_t * b
    cdef samfile_t * fp

    def __cinit__(self, Samfile samfile):

        if not samfile._isOpen():
            raise ValueError( "I/O operation on closed file" )

        if samfile.isbam: mode = "rb"
        else: mode = "r"

        # reopen the file to avoid iterator conflict
        store = StderrStore()
        self.fp = samopen( samfile.filename, mode, NULL )
        store.release()

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
        bam_destroy1(self.b)
        samclose( self.fp )


cdef class IteratorRowAllRefs(IteratorRow):
    """iterates over all mapped reads by chaining iterators over each reference
    """
    cdef Samfile     samfile
    cdef int         tid
    cdef IteratorRowRegion rowiter

    def __cinit__(self, Samfile samfile):
        assert samfile._isOpen()
        if not samfile._hasIndex(): raise ValueError("no index available for fetch")
        self.samfile = samfile
        self.tid = -1

    def nextiter(self):
        self.rowiter = IteratorRowRegion(self.samfile, self.tid, 0, 1<<29)

    def __iter__(self):
        return self

    def __next__(self):
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        # Create an initial iterator
        if self.tid==-1:
            if not self.samfile.nreferences:
                raise StopIteration
            self.tid = 0
            self.nextiter()

        while 1:
            self.rowiter.cnext()

            # If current iterator is not exhausted, return aligned read
            if self.rowiter.retval>0:
                return makeAlignedRead(self.rowiter.b)

            self.tid += 1

            # Otherwise, proceed to next reference or stop
            if self.tid<self.samfile.nreferences:
                self.nextiter()
            else:
                raise StopIteration

ctypedef struct __iterdata:
    bamFile fp
    bam_iter_t iter

cdef int __advance( void * data, bam1_t * b ):
    cdef __iterdata * d
    d = <__iterdata*>data
    return bam_iter_read( d.fp, d.iter, b )

cdef class IteratorColumn:
    '''abstract base class for iterators over columns.

    IteratorColumn objects wrap the pileup functionality of samtools.
    
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

    # result of the last plbuf_push
    cdef IteratorRowRegion iter
    cdef int tid
    cdef int pos
    cdef int n_plp
    cdef const_bam_pileup1_t_ptr plp
    cdef bam_plp_t pileup_iter
    cdef __iterdata iterdata 
    cdef Samfile samfile

    def __iter__(self):
        return self 

    cdef int cnext(self):
        '''perform next iteration.
        '''
        self.plp = bam_plp_auto( self.pileup_iter, 
                                 &self.tid,
                                 &self.pos,
                                 &self.n_plp )

cdef class IteratorColumnRegion(IteratorColumn):
    '''iterates over a region only.
    '''
    def __cinit__(self, Samfile samfile, int tid, int start, int end ):

        self.n_plp = 0
        self.tid = 0
        self.pos = 0
        self.plp = NULL
        self.samfile = samfile
        
        # initialize iterator
        self.iter = IteratorRowRegion( samfile, tid, start, end )
        self.iterdata.fp = samfile.samfile.x.bam
        self.iterdata.iter = self.iter.iter
        self.pileup_iter = bam_plp_init( &__advance, &self.iterdata )

    def __next__(self): 
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """

        while 1:
            self.cnext()
            if self.n_plp < 0:
                raise ValueError("error during iteration" )
        
            if self.plp == NULL:
                raise StopIteration
            
            return makePileupProxy( <bam_pileup1_t*>self.plp, self.tid, self.pos, self.n_plp )

    def __dealloc__(self):
        # reset in order to avoid memory leak messages for iterators that have
        # not been fully consumed
        bam_plp_reset(self.pileup_iter)
        bam_plp_destroy(self.pileup_iter)

cdef class IteratorColumnAllRefs(IteratorColumn):
    """iterates over all columns by chaining iterators over each reference
    """

    def __cinit__(self, Samfile samfile):

        self.n_plp = 0
        self.tid = 0
        self.pos = 0
        self.plp = NULL
        self.samfile = samfile

        # no iteration over empty files
        if not samfile.nreferences: 
            raise StopIteration

        # initialize iterator
        self.iter = IteratorRowRegion( samfile, self.tid, 0, max_pos )
        self.iterdata.fp = samfile.samfile.x.bam
        self.iterdata.iter = self.iter.iter
        self.pileup_iter = bam_plp_init( &__advance, &self.iterdata )

    def __next__(self):
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        
        while 1:
            self.cnext()

            if self.n_plp < 0:
                raise ValueError("error during iteration" )
        
            # return result, if within same reference
            if self.plp != NULL:
                return makePileupProxy( <bam_pileup1_t*>self.plp, self.tid, self.pos, self.n_plp )

            # Otherwise, proceed to next reference or stop
            self.tid += 1
            if self.tid < self.samfile.nreferences:
                self.iter = IteratorRowRegion( self.samfile, self.tid, 0, max_pos )
                self.iterdata.fp = self.samfile.samfile.x.bam
                self.iterdata.iter = self.iter.iter
                self.pileup_iter = bam_plp_init( &__advance, &self.iterdata )
            else:
                raise StopIteration

    def __dealloc__(self):
        # reset in order to avoid memory leak messages for iterators that have
        # not been fully consumed
        bam_plp_reset(self.pileup_iter)
        bam_plp_destroy(self.pileup_iter)

cdef inline int32_t query_start(bam1_t *src) except -1:
    cdef uint32_t * cigar_p, op
    cdef uint32_t k
    cdef uint32_t start_offset = 0

    if src.core.n_cigar:
        cigar_p = bam1_cigar(src);
        for k from 0 <= k < src.core.n_cigar:
            op = cigar_p[k] & BAM_CIGAR_MASK
            if op==BAM_CHARD_CLIP:
                if start_offset!=0 and start_offset!=src.core.l_qseq:
                    PyErr_SetString(ValueError, 'Invalid clipping in CIGAR string')
                    return -1
            elif op==BAM_CSOFT_CLIP:
                start_offset += cigar_p[k] >> BAM_CIGAR_SHIFT
            else:
                break

    return start_offset


cdef inline int32_t query_end(bam1_t *src) except -1:
    cdef uint32_t * cigar_p, op
    cdef uint32_t k
    cdef uint32_t end_offset = src.core.l_qseq

    if src.core.n_cigar>1:
        cigar_p = bam1_cigar(src);
        for k from src.core.n_cigar > k >= 1:
            op = cigar_p[k] & BAM_CIGAR_MASK
            if op==BAM_CHARD_CLIP:
                if end_offset!=0 and end_offset!=src.core.l_qseq:
                    PyErr_SetString(ValueError, 'Invalid clipping in CIGAR string')
                    return -1
            elif op==BAM_CSOFT_CLIP:
                end_offset -= cigar_p[k] >> BAM_CIGAR_SHIFT
            else:
                break

    if end_offset==0:
        end_offset = src.core.l_qseq

    return end_offset


cdef inline object get_seq_range(bam1_t *src, uint32_t start, uint32_t end):
    cdef uint8_t * p
    cdef uint32_t k
    cdef char * s

    if not src.core.l_qseq:
        return None

    seq = PyString_FromStringAndSize(NULL, end-start)
    s   = PyString_AS_STRING(seq)
    p   = bam1_seq(src)

    for k from start <= k < end:
        # equivalent to bam_nt16_rev_table[bam1_seqi(s, i)] (see bam.c)
        # note: do not use string literal as it will be a python string
        s[k-start] = bam_nt16_rev_table[p[k/2] >> 4 * (1 - k%2) & 0xf]

    return seq


cdef inline object get_qual_range(bam1_t *src, uint32_t start, uint32_t end):
    cdef uint8_t * p
    cdef uint32_t k
    cdef char * q

    p = bam1_qual(src)
    if p[0] == 0xff:
        return None

    qual = PyString_FromStringAndSize(NULL, end-start)
    q    = PyString_AS_STRING(qual)

    for k from start <= k < end:
        ## equivalent to t[i] + 33 (see bam.c)
        q[k-start] = p[k] + 33

    return qual

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

    # Now only called when instances are created from Python
    def __init__(self):
        # see bam_init1
        self._delegate = <bam1_t*>calloc( 1, sizeof( bam1_t) )
        # allocate some memory 
        # If size is 0, calloc does not return a pointer that can be passed to free()
        # so allocate 40 bytes for a new read
        self._delegate.m_data = 40
        self._delegate.data = <uint8_t *>calloc( self._delegate.m_data, 1 )
        self._delegate.data_len = 0

    def __dealloc__(self):
        bam_destroy1(self._delegate)
    
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
    
       
    def compare(self, AlignedRead other):
        '''return -1,0,1, if contents in this are binary <,=,> to *other*'''

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

        # Fast-path test for object identity
        if t==o:
            return 0

        retval = memcmp(&t.core, &o.core, sizeof(bam1_core_t))

        if retval: return retval
        retval = cmp(t.data_len, o.data_len)
        if retval: return retval
        return memcmp(t.data, o.data, t.data_len)

    # Disabled so long as __cmp__ is a special method
    def __hash__(self):
        return _Py_HashPointer(<void *>self)

    property qname:
        """the query name (None if not present)"""
        def __get__(self):
            cdef bam1_t * src 
            src = self._delegate
            if src.core.l_qname == 0: return None
            return <char *>bam1_qname( src )

        def __set__(self, qname ):
            if qname == None or len(qname) == 0: return
            cdef bam1_t * src 
            cdef int l 
            cdef char * p

            src = self._delegate            
            p = bam1_qname( src )

            # the qname is \0 terminated
            l = len(qname) + 1
            pysam_bam_update( src, 
                              src.core.l_qname, 
                              l, 
                              <uint8_t*>p )

            src.core.l_qname = l

            # re-acquire pointer to location in memory
            # as it might have moved
            p = bam1_qname(src)

            strncpy( p, qname, l )
            
    property cigar:
        """the :term:`cigar` alignment (None if not present).
        """
        def __get__(self):
            cdef uint32_t * cigar_p
            cdef bam1_t * src 
            cdef op, l, cigar
            cdef int k

            src = self._delegate
            if src.core.n_cigar == 0: return None
            
            cigar = []
            cigar_p = bam1_cigar(src);
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
            p = bam1_cigar(src)

            # create space for cigar data within src.data
            pysam_bam_update( src, 
                              src.core.n_cigar * 4,
                              len(values) * 4, 
                              <uint8_t*>p )
            
            # length is number of cigar operations, not bytes
            src.core.n_cigar = len(values)

            # re-acquire pointer to location in memory
            # as it might have moved
            p = bam1_cigar(src)

            # insert cigar operations
            for op, l in values:
                p[k] = l << BAM_CIGAR_SHIFT | op
                k += 1

            ## setting the cigar string also updates the "bin" attribute
            src.core.bin = bam_reg2bin( src.core.pos, bam_calend( &src.core, p))

    property seq:
        """read sequence bases, including :term:`soft clipped` bases (None if not present)"""
        def __get__(self):
            cdef bam1_t * src
            cdef char * s
            src = self._delegate

            if src.core.l_qseq == 0: return None

            return get_seq_range(src, 0, src.core.l_qseq)

        def __set__(self,seq):
            # samtools manages sequence and quality length memory together
            # if no quality information is present, the first byte says 0xff.
            
            if seq == None or len(seq) == 0: return
            cdef bam1_t * src
            cdef uint8_t * p 
            cdef char * s
            cdef int l, k, nbytes_new, nbytes_old

            src = self._delegate

            l = len(seq)
            
            # as the sequence is stored in half-bytes, the total length (sequence
            # plus quality scores) is (l+1)/2 + l
            nbytes_new = (l+1)/2 + l
            nbytes_old = (src.core.l_qseq+1)/2 + src.core.l_qseq
            # acquire pointer to location in memory
            p = bam1_seq( src )
            src.core.l_qseq = l

            pysam_bam_update( src, 
                              nbytes_old,
                              nbytes_new,
                              p)
            # re-acquire pointer to location in memory
            # as it might have moved
            p = bam1_seq( src )
            for k from 0 <= k < nbytes_new: p[k] = 0
            # convert to C string
            s = seq
            for k from 0 <= k < l:
                p[k/2] |= pysam_translate_sequence(s[k]) << 4 * (1 - k % 2)

            # erase qualities
            p = bam1_qual( src )
            p[0] = 0xff


    property qual:
        """read sequence base qualities, including :term:`soft clipped` bases (None if not present)"""
        def __get__(self):

            cdef bam1_t * src
            cdef char * q

            src = self._delegate

            if src.core.l_qseq == 0: return None

            return get_qual_range(src, 0, src.core.l_qseq)

        def __set__(self,qual):
            # note that space is already allocated via the sequences
            cdef bam1_t * src
            cdef uint8_t * p
            cdef char * q 
            cdef int k

            src = self._delegate
            p = bam1_qual( src )
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

    property query:
        """aligned portion of the read and excludes any flanking bases that were :term:`soft clipped` (None if not present)

        SAM/BAM files may included extra flanking bases sequences that were
        not part of the alignment.  These bases may be the result of the
        Smith-Waterman or other algorithms, which may not require alignments
        that begin at the first residue or end at the last.  In addition,
        extra sequencing adapters, multiplex identifiers, and low-quality bases that
        were not considered for alignment may have been retained."""

        def __get__(self):
            cdef bam1_t * src
            cdef uint32_t start, end
            cdef char * s

            src = self._delegate

            if src.core.l_qseq == 0: return None

            start = query_start(src)
            end   = query_end(src)

            return get_seq_range(src, start, end)

    property qqual:
        """aligned query sequence quality values (None if not present)"""
        def __get__(self):
            cdef bam1_t * src
            cdef uint32_t start, end
            cdef char * q

            src = self._delegate

            if src.core.l_qseq == 0: return None

            start = query_start(src)
            end   = query_end(src)

            return get_qual_range(src, start, end)

    property qstart:
        """start index of the aligned query portion of the sequence (0-based, inclusive)"""
        def __get__(self):
            return query_start(self._delegate)

    property qend:
        """end index of the aligned query portion of the sequence (0-based, exclusive)"""
        def __get__(self):
            return query_end(self._delegate)

    property qlen:
        """Length of the aligned query sequence"""
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            return query_end(src)-query_start(src)

    property tags:
        """the tags in the AUX field.

        This property permits convenience access to 
        the tags. Changes it the returned list will
        not update the tags automatically. Instead,
        the following is required for adding a 
        new tag::

            read.tags = read.tags + [("RG",0)]

        """
        def __get__(self):
            cdef char * ctag
            cdef bam1_t * src
            cdef uint8_t * s
            cdef char auxtag[3]
            cdef char auxtype
            
            src = self._delegate
            if src.l_aux == 0: return None
            
            s = bam1_aux( src )
            result = []
            auxtag[2] = 0
            while s < (src.data + src.data_len):
                # get tag
                auxtag[0] = s[0]
                auxtag[1] = s[1]
                s += 2
                auxtype = s[0]

                if auxtype in ('c', 'C'):
                    value = <int>bam_aux2i(s)            
                    s += 1
                elif auxtype in ('s', 'S'):
                    value = <int>bam_aux2i(s)            
                    s += 2
                elif auxtype in ('i', 'I'):
                    value = <float>bam_aux2i(s)
                    s += 4
                elif auxtype == 'f':
                    value = <float>bam_aux2f(s)
                    s += 4
                elif auxtype == 'd':
                    value = <double>bam_aux2d(s)
                    s += 8
                elif auxtype == 'A':
                    value = "%c" % <char>bam_aux2A(s)
                    s += 1
                elif auxtype in ('Z', 'H'):
                    value = <char*>bam_aux2Z(s)
                    # +1 for NULL terminated string
                    s += len(value) + 1
                 # 
                s += 1
  
                result.append( (auxtag, value) )

            return result

        def __set__(self, tags):
            cdef char * ctag
            cdef bam1_t * src
            cdef uint8_t * s
            cdef uint8_t * new_data
            cdef char * temp 
            cdef int guessed_size, control_size
            cdef int max_size, size, offset

            src = self._delegate
            max_size = 4000
            offset = 0

            if tags != None: 

                # map samtools code to python.struct code and byte size
                buffer = ctypes.create_string_buffer(max_size)

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
            # if offset == 0, the aux field will be 
            # empty
            pysam_bam_update( src, 
                              src.l_aux,
                              offset,
                              bam1_aux( src ) )
            
            src.l_aux = offset

            # copy data only if there is any
            if offset != 0:

                # get location of new data
                s = bam1_aux( src )            
            
                # check if there is direct path from buffer.raw to tmp
                temp = buffer.raw
                memcpy( s, temp, offset )            

    property flag: 
        """properties flag"""
        def __get__(self): return self._delegate.core.flag
        def __set__(self, flag): self._delegate.core.flag = flag
    property rname: 
        """
        :term:`target` ID

        .. note::

            This field contains the index of the reference sequence 
            in the sequence dictionary. To obtain the name
            of the reference sequence, use :meth:`pysam.Samfile.getrname()`

        """
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
                src.core.bin = bam_reg2bin( src.core.pos, bam_calend( &src.core, bam1_cigar(src)) )
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
    property aend:
        '''aligned end position of the read (read only).  Returns
        None if not available.'''
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            if (self.flag & BAM_FUNMAP) or src.core.n_cigar == 0:
                return None
            return bam_calend(&src.core, bam1_cigar(src))
    property alen:
        '''aligned length of the read (read only).  Returns None if
        not available.'''
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            if (self.flag & BAM_FUNMAP) or src.core.n_cigar == 0:
                return None
            return bam_calend(&src.core, 
                               bam1_cigar(src)) - \
                               self._delegate.core.pos

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
    cdef bam_pileup1_t * plp
    cdef int tid
    cdef int pos
    cdef int n_pu
    
    def __init__(self):
        raise TypeError("This class cannot be instantiated from Python")

    def __str__(self):
        return "\t".join( map(str, (self.tid, self.pos, self.n))) +\
            "\n" +\
            "\n".join( map(str, self.pileups) )

    property tid:
        '''the chromosome ID as is defined in the header'''
        def __get__(self): return self.tid

    property n:
        '''number of reads mapping to this column.'''
        def __get__(self): return self.n_pu
        def __set__(self, n): self.n_pu = n

    property pos:
        def __get__(self): return self.pos

    property pileups:
        '''list of reads (:class:`pysam.PileupRead`) aligned to this column'''
        def __get__(self):
            cdef int x
            pileups = []
            # warning: there could be problems if self.n and self.buf are
            # out of sync.
            for x from 0 <= x < self.n_pu:
                pileups.append( makePileupRead( &(self.plp[x])) )
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

    def __init__(self):
        raise TypeError("This class cannot be instantiated from Python")

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

cdef class SNPCall:
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
       def __get__(self): return PyString_FromStringAndSize( &self._reference_base, 1 )

    property genotype:
       '''the genotype called.'''
       def __get__(self): return PyString_FromStringAndSize( &self._genotype, 1 )

    property consensus_quality:
       '''the genotype quality (Phred-scaled).'''
       def __get__(self): return self._consensus_quality

    property snp_quality:
       '''the snp quality (Phred scaled) - probability of consensus being identical to reference sequence.'''
       def __get__(self): return self._snp_quality

    property mapping_quality:
       '''the root mean square (rms) of the mapping quality of all reads involved in the call.'''
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

cdef class SNPCaller:
    '''Proxy for samtools SNP caller.

    This caller is fast for calling few SNPs in selected regions.

    It is slow, if called over large genomic regions.
     '''

    cdef Fastafile fasta
    cdef bam_maqcns_t * c
    cdef char * seq
    cdef int seq_tid 
    cdef int seq_length
    cdef Samfile samfile

    def __cinit__(self, 
                  Samfile samfile,
                  Fastafile fasta ):

        self.samfile = samfile
        self.c =  bam_maqcns_init()
        self.c.is_soap = 1
        self.fasta = fasta
        self.seq = NULL
        self.seq_tid = -1
        self.seq_length = 0
        bam_maqcns_prepare( self.c )

    def call(self, reference, int pos ): 
        """return a :class:`SNPCall` object.
        """

        cdef int rb, tid

        tid = self.samfile.gettid( reference )

        # update sequence
        if tid != self.seq_tid:
            if self.seq != NULL: free( self.seq )
            self.seq = NULL
            self.seq_tid = tid
     
        # reload sequence if necessary
        if self.seq == NULL:
            rname = self.samfile.getrname( self.seq_tid )
            self.seq = self.fasta._fetch( rname,
                                          0, max_pos, 
                                          &self.seq_length )
            if self.seq == NULL:
                raise ValueError( "reference sequence for '%s' (tid=%i) not found" % \
                                      (rname,
                                       self.seq_tid))

        # reference base
        if pos >= self.seq_length:
            raise ValueError( "position %i out of bounds on reference sequence (len=%i)" % (pos, self.seq_length) )

        rb = self.seq[pos]

        # initialize pileup engine
        cdef IteratorColumn itr = IteratorColumnRegion( self.samfile, tid, pos, pos + 1 )

        while 1:
            itr.cnext()
            
            if itr.n_plp < 0:
                raise ValueError("error during iteration" )

            if itr.plp == NULL:
                raise ValueError( "no reads in region - no call" )
             
            if itr.pos == pos: break

        cdef uint32_t cns = bam_maqcns_call( itr.n_plp, 
                                             itr.plp, 
                                             self.c )

        cdef int ref_q, rb4 
        rb4 = bam_nt16_table[rb]
        ref_q = 0
        if rb4 != 15 and cns>>28 != 15 and cns>>28 != rb4:
            # a SNP
            if cns >> 24 & 0xf == rb4:
                ref_q = cns >> 8 & 0xff 
            else: 
                ref_q = (cns >> 8 & 0xff) + (cns & 0xff)
                
            if ref_q > 255: ref_q = 255

        cdef SNPCall call

        call = SNPCall()
        call._tid = itr.tid
        call._pos = itr.pos
        call._reference_base = rb
        call._genotype = bam_nt16_rev_table[cns>>28]
        call._consensus_quality = cns >> 8 & 0xff
        call._snp_quality = ref_q
        call._rms_mapping_quality = cns >> 16 & 0xff
        call._coverage = itr.n_plp

        return call 

    def __dealloc__(self):
        if self.seq != NULL: free( self.seq )
        bam_maqcns_destroy( self.c )

cdef class IteratorSNPCalls:
    """call SNPs within a region.

    This caller is fast if SNPs are called over large continuous
    regions. It is slow, if instantiated frequently.

    .. todo:

        * allow options similar to the samtools pileup command
        * check for efficiency
    """

    cdef IteratorColumn iter
    cdef Fastafile fasta
    cdef bam_maqcns_t * c
    cdef char * seq
    cdef int seq_tid 
    cdef int seq_length
    def __cinit__(self, 
                  IteratorColumn iterator_column,
                  Fastafile fasta ):

        self.iter = iterator_column
        self.c =  bam_maqcns_init()
        self.c.is_soap = 1
        self.fasta = fasta
        self.seq = NULL
        self.seq_tid = -1
        self.seq_length = 0
        bam_maqcns_prepare( self.c )

    def __iter__(self):
        return self 

    def __next__(self): 
        """python version of next().
        """

        # the following code was adapted from bam_plcmd.c:pileup_func()
        self.iter.cnext()

        if self.iter.n_plp < 0:
            raise ValueError("error during iteration" )

        if self.iter.plp == NULL:
            raise StopIteration

        cdef int rb

        # update sequence
        if self.iter.tid != self.seq_tid:
            if self.seq != NULL: free( self.seq )
            self.seq = NULL
            self.seq_tid = self.iter.tid

        # reload sequence if necessary
        if self.seq == NULL:
            # TODO: look for a  shortcut 
            self.seq = self.fasta._fetch( self.iter.iter.samfile.samfile.header.target_name[self.seq_tid], 
                                          0, max_pos, 
                                          &self.seq_length )
            if self.seq == NULL:
                raise ValueError( "reference sequence for '%s' (tid=%i) not found" % \
                                       (self.iter.iter.samfile.samfile.header.target_name[self.seq_tid], 
                                        self.seq_tid))

        # reference base
        if self.iter.pos >= self.seq_length:
            raise ValueError( "position %i out of bounds on reference sequence (len=%i)" % (self.iter.pos, self.seq_length) )

        rb = self.seq[self.iter.pos]

        cdef uint32_t cns 
        cns = bam_maqcns_call( self.iter.n_plp, 
                               self.iter.plp, 
                               self.c )

        cdef int ref_q, rb4 
        rb4 = bam_nt16_table[rb]
        ref_q = 0
        if rb4 != 15 and cns>>28 != 15 and cns>>28 != rb4:
            # a SNP
            if cns >> 24 & 0xf == rb4:
                ref_q = cns >> 8 & 0xff 
            else: 
                ref_q = (cns >> 8 & 0xff) + (cns & 0xff)

            if ref_q > 255: ref_q = 255

        cdef SNPCall call

        call = SNPCall()
        call._tid = self.iter.tid
        call._pos = self.iter.pos
        call._reference_base = rb
        call._genotype = bam_nt16_rev_table[cns>>28]
        call._consensus_quality = cns >> 8 & 0xff
        call._snp_quality = ref_q
        call._rms_mapping_quality = cns >> 16&0xff
        call._coverage = self.iter.n_plp

        return call 

    def __dealloc__(self):
        if self.seq != NULL: free( self.seq )
        bam_maqcns_destroy( self.c )

         
__all__ = ["Samfile", 
           "Fastafile",
           "IteratorRow", 
           "IteratorColumn", 
           "AlignedRead", 
           "PileupColumn", 
           "PileupProxy", 
           "PileupRead",
           "IteratorSNPCalls",
           "SNPCaller" ]

               

