# cython: embedsignature=True
# cython: profile=True
# adds doc-strings for sphinx
import tempfile
import os
import sys
import types
import itertools
import struct
import ctypes
import collections
import re
import platform
from cpython cimport PyString_FromStringAndSize, PyString_AS_STRING
from cpython cimport PyErr_SetString

#from cpython.string cimport PyString_FromStringAndSize, PyString_AS_STRING
#from cpython.exc    cimport PyErr_SetString, PyErr_NoMemory

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
# hard-coded constants
cdef char * bam_nt16_rev_table = "=ACMGRSVTWYHKDBN"
cdef int max_pos = 2 << 29

# redirect stderr to 0
_logfile = open(os.path.devnull, "w")
pysam_set_stderr( PyFile_AsFile( _logfile ) )

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
        return
        self.stderr_h, self.stderr_f = tempfile.mkstemp()
        self.stderr_save = Outs( sys.stderr.fileno() )
        self.stderr_save.setfd( self.stderr_h )
        
    def readAndRelease( self ):
        return []
        self.stderr_save.restore()
        lines = []
        if os.path.exists(self.stderr_f):
            lines = open( self.stderr_f, "r" ).readlines()
            os.remove( self.stderr_f )
        return lines

    def release(self):
        return
        self.stderr_save.restore()
        if os.path.exists(self.stderr_f):
            os.remove( self.stderr_f )

    def __del__(self):
        self.release()

class StderrStoreWindows():
    '''does nothing. stderr can't be redirected on windows'''
    def __init__(self): pass
    def readAndRelease(self): return []
    def release(self): pass

if platform.system()=='Windows':
    del StderrStore
    StderrStore = StderrStoreWindows


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
                        "PG" : { "PN" : str, "ID" : str, "VN" : str, "CL" : str }, }

# output order of fields within records
VALID_HEADER_ORDER = { "HD" : ( "VN", "SO", "GO" ),
                       "SQ" : ( "SN", "LN", "AS", "M5" , "UR" , "SP" ),
                       "RG" : ( "ID", "SM", "LB", "DS" , "PU" , "PI" , "CN" , "DT", "PL" ),
                       "PG" : ( "PN", "ID", "VN", "CL" ), }


######################################################################
######################################################################
######################################################################
## Public methods
######################################################################
cdef class Fastafile:
    '''*(filename)*
              
    A *FASTA* file. The file is automatically opened.

    The file expects an indexed fasta file.

    TODO: 
        add automatic indexing.
        add function to get sequence names.
    '''

    cdef char * _filename
    # pointer to fastafile
    cdef faidx_t * fastafile

    def __cinit__(self, *args, **kwargs ):
        self.fastafile = NULL
        self._filename = NULL
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
        if self._filename != NULL: free(self._filename)
        self._filename = strdup(filename)
        self.fastafile = fai_load( filename )

        if self.fastafile == NULL:
            raise IOError("could not open file `%s`" % filename )

    def close( self ):
        if self.fastafile != NULL:
            fai_destroy( self.fastafile )
            self.fastafile = NULL

    def __dealloc__(self):
        self.close()
        if self._filename != NULL: free(self._filename)

    property filename:
        '''number of :term:`filename` associated with this object.'''
        def __get__(self):
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
            return self._filename

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

#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
cdef int count_callback( bam1_t *alignment, void *f):
     '''callback for bam_fetch - count number of reads.
     '''
     cdef int* counter = (<int*>f)
     counter[0] += 1;

ctypedef struct MateData:
     char * name
     bam1_t * mate
     uint32_t flag

#------------------------------------------------------------------------
#------------------------------------------------------------------------
#------------------------------------------------------------------------
cdef int mate_callback( bam1_t *alignment, void *f):
     '''callback for bam_fetch = filter mate
     '''
     cdef MateData * d = (<MateData*>f)
     # printf("mate = %p, name1 = %s, name2=%s\t%i\t%i\t%i\n", 
     #        d.mate, d.name, bam1_qname(alignment),
     #        d.flag, alignment.core.flag, alignment.core.flag & d.flag)

     if d.mate == NULL: 
         # could be sped up by comparing the lengths of query strings first
         # using l_qname
         #
         # also, make sure that we get the other read by comparing 
         # the flags
         if alignment.core.flag & d.flag != 0 and \
                 strcmp( bam1_qname( alignment ), d.name ) == 0:
             d.mate = bam_dup1( alignment )


cdef class Samfile:
    '''*(filename, mode=None, template = None, referencenames = None, referencelengths = None, text = NULL, header = None,
         add_sq_text = False )*
              
    A :term:`SAM`/:term:`BAM` formatted file. The file is automatically opened.
    
    *mode* should be ``r`` for reading or ``w`` for writing. The default is text mode (:term:`SAM`). For binary 
    (:term:`BAM`) I/O you should append ``b`` for compressed or ``u`` for uncompressed :term:`BAM` output. 
    Use ``h`` to output header information in text (:term:`TAM`)  mode.

    If ``b`` is present, it must immediately follow ``r`` or ``w``. 
    Valid modes are ``r``, ``w``, ``wh``, ``rb``, ``wb`` and ``wbu``. For instance, to open 
    a :term:`BAM` formatted file for reading, type::

        f = pysam.Samfile('ex1.bam','rb')

    If mode is not specified, we will try to auto-detect in the order 'rb', 'r', thus both the following
    should work::

        f1 = pysam.Samfile('ex1.bam' )
        f2 = pysam.Samfile('ex1.sam' )

    If an index for a BAM file exists (.bai), it will be opened automatically. Without an index random
    access to reads via :meth:`fetch` and :meth:`pileup` is disabled.

    For writing, the header of a :term:`SAM` file/:term:`BAM` file can be constituted from several
    sources (see also the samtools format specification):

        1. If *template* is given, the header is copied from a another *Samfile* 
           (*template* must be of type *Samfile*).

        2. If *header* is given, the header is built from a multi-level dictionary. The first level 
           are the four types ('HD', 'SQ', ...). The second level are a list of lines, with each line 
           being a list of tag-value pairs.

        3. If *text* is given, new header text is copied from raw text.

        4. The names (*referencenames*) and lengths (*referencelengths*) are supplied directly as lists. 
           By default, 'SQ' and 'LN' tags will be added to the header text. This option can be
           changed by unsetting the flag *add_sq_text*. 

    '''

    def __cinit__(self, *args, **kwargs ):
        self.samfile = NULL
        self._filename = NULL
        self.isbam = False
        self.isstream = False
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
               mode = None,
               Samfile template = None,
               referencenames = None,
               referencelengths = None,
               text = None,
               header = None,
               port = None,
               add_sq_text = True,
              ):
        '''open a sam/bam file.

        If _open is called on an existing bamfile, the current file will be
        closed and a new file will be opened.
        '''

        # read mode autodetection
        if mode is None:
            try:
                self._open(filename, 'rb', template=template,
                           referencenames=referencenames,
                           referencelengths=referencelengths,
                           text=text, header=header, port=port)
                return
            except ValueError, msg:
                pass
            
            self._open(filename, 'r', template=template,
                       referencenames=referencenames,
                       referencelengths=referencelengths,
                       text=text, header=header, port=port)
            return

        assert mode in ( "r","w","rb","wb", "wh", "wbu", "rU" ), "invalid file opening mode `%s`" % mode
        assert filename != NULL

        # close a previously opened file
        if self.samfile != NULL: self.close()
        self.samfile = NULL

        cdef bam_header_t * header_to_write
        header_to_write = NULL
        
        if self._filename != NULL: free(self._filename )
        self._filename = strdup( filename )
        self.isstream = strcmp( filename, "-" ) == 0

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

                # Optionally, if there is no text, add a SAM compatible header to output
                # file.
                if text is None and add_sq_text:
                    text = ''
                    for x from 0 <= x < header_to_write.n_targets:
                        text += "@SQ\tSN:%s\tLN:%s\n" % (referencenames[x], referencelengths[x] )

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

            # try to detect errors
            self.samfile = samopen( filename, mode, NULL )
            if self.samfile == NULL:
                raise ValueError( "could not open file (mode='%s') - is it SAM/BAM format?" % mode)

            if self.samfile.header == NULL:
                raise ValueError( "file does not have valid header (mode='%s') - is it SAM/BAM format?" % mode )
            
            #disabled for autodetection to work
            # needs to be disabled so that reading from sam-files without headers works
            #if self.samfile.header.n_targets == 0:
            #    raise ValueError( "file header is empty (mode='%s') - is it SAM/BAM format?" % mode)

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

            if not self.isstream:
                self.start_offset = bam_tell( self.samfile.x.bam )

    def gettid( self, reference ):
        '''
        convert :term:`reference` name into numerical :term:`tid`

        returns -1 if reference is not known.
        '''
        if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
        return pysam_reference2tid( self.samfile.header, reference )

    def getrname( self, tid ):
        '''
        convert numerical :term:`tid` into :term:`reference` name.'''
        if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
        if not 0 <= tid < self.samfile.header.n_targets:
            raise ValueError( "tid %i out of range 0<=tid<%i" % (tid, self.samfile.header.n_targets ) )
        return self.samfile.header.target_name[tid]

    cdef char * _getrname( self, int tid ):
        '''
        convert numerical :term:`tid` into :term:`reference` name.'''
        if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
        if not 0 <= tid < self.samfile.header.n_targets:
            raise ValueError( "tid %i out of range 0<=tid<%i" % (tid, self.samfile.header.n_targets ) )
        return self.samfile.header.target_name[tid]

    def _parseRegion( self, 
                      reference = None, 
                      start = None, 
                      end = None,
                      region = None ):
        '''
        parse region information.

        raise ValueError for for invalid regions.

        returns a tuple of flag, tid, start and end. Flag indicates
        whether some coordinates were supplied.

        Note that regions are 1-based, while start,end are python coordinates.
        '''
        # This method's main objective is to translate from a reference to a tid. 
        # For now, it calls bam_parse_region, which is clumsy. Might be worth
        # implementing it all in pysam (makes use of khash).
        
        cdef int rtid
        cdef long long rstart
        cdef long long rend

        rtid = -1
        rstart = 0
        rend = max_pos
        if start != None: 
            try:
                rstart = start
            except OverflowError:
                raise ValueError( 'start out of range (%i)' % start )
            
        if end != None: 
            try:
                rend = end
            except OverflowError:
                raise ValueError( 'end out of range (%i)' % end )

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
    
    def reset( self ):
        '''reset file position to beginning of read section.'''
        return self.seek( self.start_offset, 0 )

    def seek( self, uint64_t offset, int where = 0):
        '''
        move file pointer to position *offset*, see :meth:`pysam.Samfile.tell`.
        '''

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )
        if not self.isbam:
            raise NotImplementedError("seek only available in bam files")
        if self.isstream:
            raise OSError("seek no available in streams")
        
        return bam_seek( self.samfile.x.bam, offset, where )

    def tell( self ):
        '''
        return current file position
        '''
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
        '''
        fetch aligned reads in a :term:`region` using 0-based indexing. The region is specified by
        :term:`reference`, *start* and *end*. Alternatively, a samtools :term:`region` string can 
        be supplied.

        Without *reference* or *region* all mapped reads will be fetched. The reads will be returned
        ordered by reference sequence, which will not necessarily be the order within the file.

        If *until_eof* is given, all reads from the current file position will be returned
        in order as they are within the file. Using this option will also fetch unmapped reads. 
        
        If only *reference* is set, all reads aligned to *reference* will be fetched.

        The method returns an iterator of type :class:`pysam.IteratorRow` unless
        a *callback is provided. If *callback* is given, the callback will be executed 
        for each position within the :term:`region`. Note that callbacks currently work
        only, if *region* or *reference* is given.

        Note that a :term:`SAM` file does not allow random access. If *region* or *reference* are given,
        an exception is raised.
        '''
        cdef int rtid, rstart, rend, has_coord

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        has_coord, rtid, rstart, rend = self._parseRegion( reference, start, end, region )
        
        if self.isstream: reopen = False
        else: reopen = True

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
                    return IteratorRowRegion( self, rtid, rstart, rend, reopen=reopen )
                else:
                    if until_eof:
                        return IteratorRowAll( self, reopen=reopen )
                    else:
                        # AH: check - reason why no reopen for AllRefs?
                        return IteratorRowAllRefs(self ) # , reopen=reopen )
        else:   
            # check if header is present - otherwise sam_read1 aborts
            # this happens if a bamfile is opened with mode 'r'
            if has_coord:
                raise ValueError ("fetching by region is not available for sam files" )

            if self.samfile.header.n_targets == 0:
                raise ValueError( "fetch called for samfile without header")

            if callback:
                raise NotImplementedError( "callback not implemented yet" )
            else:
                return IteratorRowAll( self, reopen=reopen )

    def mate( self, 
              AlignedRead read ):
        '''return the mate of :class:`AlignedRead` *read*.

        Throws a ValueError if read is unpaired or the mate
        is unmapped.

        .. note::
            Calling this method will change the file position.
            This might interfere with any iterators that have
            not re-opened the file.

        '''
        cdef uint32_t flag = read._delegate.core.flag

        if flag & BAM_FPAIRED == 0:
            raise ValueError( "read %s: is unpaired" % (read.qname))
        if flag & BAM_FMUNMAP != 0:
            raise ValueError( "mate %s: is unmapped" % (read.qname))
        
        cdef MateData mate_data

        mate_data.name = <char *>bam1_qname(read._delegate)
        mate_data.mate = NULL
        # xor flags to get the other mate
        cdef int x = BAM_FREAD1 + BAM_FREAD2
        mate_data.flag = ( flag ^ x) & x

        bam_fetch(self.samfile.x.bam, 
                  self.index, 
                  read._delegate.core.mtid, 
                  read._delegate.core.mpos,
                  read._delegate.core.mpos + 1,
                  <void*>&mate_data, 
                  mate_callback )

        if mate_data.mate == NULL:
            raise ValueError( "mate not found" )

        cdef AlignedRead dest = AlignedRead.__new__(AlignedRead)
        dest._delegate = mate_data.mate
        return dest

    def count( self, 
               reference = None, 
               start = None, 
               end = None, 
               region = None, 
               until_eof = False ):
        '''*(reference = None, start = None, end = None, region = None, callback = None, until_eof = False)*
               
        count  reads :term:`region` using 0-based indexing. The region is specified by
        :term:`reference`, *start* and *end*. Alternatively, a samtools :term:`region` string can be supplied.

        Note that a :term:`TAM` file does not allow random access. If *region* or *reference* are given,
        an exception is raised.
        '''
        cdef int rtid
        cdef int rstart
        cdef int rend

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )
        
        region, rtid, rstart, rend = self._parseRegion( reference, start, end, region )

        cdef int counter
        counter = 0;

        if self.isbam:
            if not until_eof and not self._hasIndex() and not self.isremote: 
                raise ValueError( "fetch called on bamfile without index" )

            if not region:
                raise ValueError( "counting functionality requires a region/reference" )
            if not self._hasIndex(): raise ValueError( "no index available for fetch" )
            bam_fetch(self.samfile.x.bam, 
                             self.index, 
                             rtid, 
                             rstart, 
                             rend, 
                             <void*>&counter, 
                             count_callback )
            return counter
        else:   
            raise ValueError ("count for a region is not available for sam files" )

    def pileup( self, 
                reference = None, 
                start = None, 
                end = None, 
                region = None, 
                callback = None,
                **kwargs ):
        '''
        perform a :term:`pileup` within a :term:`region`. The region is specified by
        :term:`reference`, *start* and *end* (using 0-based indexing). 
        Alternatively, a samtools *region* string can be supplied.

        Without *reference* or *region* all reads will be used for the pileup. The reads will be returned
        ordered by :term:`reference` sequence, which will not necessarily be the order within the file.

        The method returns an iterator of type :class:`pysam.IteratorColumn` unless
        a *callback is provided. If a *callback* is given, the callback will be executed 
        for each column within the :term:`region`. 

        Note that :term:`SAM` formatted files do not allow random access. 
        In these files, if a *region* or *reference* are given an exception is raised.
        
        Optional *kwargs* to the iterator:

        stepper
           The stepper controlls how the iterator advances. 
           Possible options for the stepper are 
       
           ``all``
              use all reads for pileup.
           ``samtools``
              same filter and read processing as in :term:`csamtools` pileup

        fastafile
           A :class:`FastaFile` object

         mask
           Skip all reads with bits set in mask.

         max_depth
           Maximum read depth permitted. The default limit is *8000*.

        .. note::

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
                    return IteratorColumnRegion( self, 
                                                 tid = rtid, 
                                                 start = rstart, 
                                                 end = rend, 
                                                 **kwargs )
                else:
                    return IteratorColumnAllRefs(self, **kwargs )

        else:
            raise NotImplementedError( "pileup of samfiles not implemented yet" )

    def close( self ):
        '''
        closes the :class:`pysam.Samfile`.'''
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
        if self._filename != NULL: free( self._filename )

    cpdef int write( self, AlignedRead read ) except -1:
        '''
        write a single :class:`pysam.AlignedRead` to disk.

        returns the number of bytes written.
        '''
        if not self._isOpen():
            return 0

        return samwrite( self.samfile, read._delegate )

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    ###############################################################
    ###############################################################
    ###############################################################
    ## properties
    ###############################################################
    property filename:
        '''number of :term:`filename` associated with this object.'''
        def __get__(self):
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
            return self._filename
        
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
        """tuple of the lengths of the :term:`reference` sequences. The lengths are in the same order as 
        :attr:`pysam.Samfile.references`
        """
        def __get__(self): 
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
            t = []
            for x from 0 <= x < self.samfile.header.n_targets:
                t.append( self.samfile.header.target_len[x] )
            return tuple(t)

    property mapped:
        """total number of mapped reads in file.
        """
        def __get__(self):
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
            if not self.isbam: raise AttributeError( "Samfile.mapped only available in bam files" )
            
            cdef int tid
            cdef uint32_t total = 0
            for tid from 0 <= tid < self.samfile.header.n_targets:
                total += pysam_get_mapped( self.index, tid )
            return total

    property unmapped:
        """total number of unmapped reads in file.
        """
        def __get__(self):
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
            if not self.isbam: raise AttributeError( "Samfile.unmapped only available in bam files" )
            cdef int tid
            cdef uint32_t total = 0
            for tid from 0 <= tid < self.samfile.header.n_targets:
                total += pysam_get_unmapped( self.index, tid )
            # get unmapped reads without coordinates
            total += pysam_get_unmapped( self.index, -1 )
            return total

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
                        # uppercase keys must be valid
                        # lowercase are permitted for user fields
                        if key in VALID_HEADER_FIELDS[record]:
                            x[key] = VALID_HEADER_FIELDS[record][key](value)
                        elif not key.isupper():
                            x[key] = value
                        else:
                            raise ValueError( "unknown field code '%s' in record '%s'" % (key, record) )

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
            # write fields of the specification
            for key in VALID_HEADER_ORDER[record]:
                if key in fields:
                    line.append( "%s:%s" % (key, str(fields[key])))
            # write user fields
            for key in fields:
                if not key.isupper():
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

    ###############################################################
    ###############################################################
    ###############################################################
    ## file-object like iterator access
    ## note: concurrent access will cause errors (see IteratorRow
    ## and reopen)
    ## Possible solutions: deprecate or open new file handle
    ###############################################################
    def __iter__(self):
        if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
        if not self.isbam and self.samfile.header.n_targets == 0:
                raise NotImplementedError( "can not iterate over samfile without header")
        return self 

    cdef bam1_t * getCurrent( self ):
        return self.b

    cdef int cnext(self):
        '''
        cversion of iterator. Used by :class:`pysam.Samfile.IteratorColumn`.
        '''
        cdef int ret
        return samread(self.samfile, self.b)

    def __next__(self): 
        """
        python version of next().
        """
        cdef int ret
        ret = samread(self.samfile, self.b)
        if (ret > 0):
            return makeAlignedRead( self.b )
        else:
            raise StopIteration

##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
cdef class IteratorRow:
    '''abstract base class for iterators over mapped reads.

    Various iterators implement different behaviours for wrapping around
    contig boundaries. Examples include:

    :class:`pysam.IteratorRowRegion`
        iterate within a single contig and a defined region.

    :class:`pysam.IteratorRowAll`
        iterate until EOF. This iterator will also include unmapped reads.

    :class:`pysam.IteratorRowAllRefs`
        iterate over all reads in all reference sequences.
        
    The method :meth:`Samfile.fetch` returns an IteratorRow.
    '''
    pass

cdef class IteratorRowRegion(IteratorRow):
    """*(Samfile samfile, int tid, int beg, int end, int reopen = True )*

    iterate over mapped reads in a region.

    By default, the file is re-openend to avoid conflicts between
    multiple iterators working on the same file. Set *reopen* = False
    to not re-open *samfile*.

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
    # true if samfile belongs to this object
    cdef int owns_samfile

    def __cinit__(self, Samfile samfile, int tid, int beg, int end, int reopen = True ):

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
        if reopen:
            store = StderrStore()            
            self.fp = samopen( samfile._filename, mode, NULL )
            store.release()
            assert self.fp != NULL
            self.owns_samfile = True
        else:
            self.fp = self.samfile.samfile
            self.owns_samfile = False

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
        if self.owns_samfile: samclose( self.fp )

cdef class IteratorRowAll(IteratorRow):
    """*(Samfile samfile, int reopen = True)*

    iterate over all reads in *samfile*

    By default, the file is re-openend to avoid conflicts between
    multiple iterators working on the same file. Set *reopen* = False
    to not re-open *samfile*.
    """

    # cdef bam1_t * b
    # cdef samfile_t * fp
    # # true if samfile belongs to this object
    # cdef int owns_samfile

    def __cinit__(self, Samfile samfile, int reopen = True ):

        if not samfile._isOpen():
            raise ValueError( "I/O operation on closed file" )

        if samfile.isbam: mode = "rb"
        else: mode = "r"

        # reopen the file to avoid iterator conflict
        if reopen:
            store = StderrStore()
            self.fp = samopen( samfile._filename, mode, NULL )
            store.release()
            assert self.fp != NULL
            self.owns_samfile = True
        else:
            self.fp = samfile.samfile
            self.owns_samfile = False

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
        if self.owns_samfile: samclose( self.fp )

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

cdef class IteratorRowSelection(IteratorRow):
    """*(Samfile samfile)*

    iterate over reads in *samfile* at a given list of file positions.
    """

    cdef bam1_t * b
    cdef int current_pos 
    cdef samfile_t * fp
    cdef positions
    # true if samfile belongs to this object
    cdef int owns_samfile

    def __cinit__(self, Samfile samfile, positions, int reopen = True ):

        if not samfile._isOpen():
            raise ValueError( "I/O operation on closed file" )

        if not samfile._isOpen():
            raise ValueError( "I/O operation on closed file" )

        assert samfile.isbam, "can only use this iterator on bam files"
        mode = "rb"

        # reopen the file to avoid iterator conflict
        if reopen:
            store = StderrStore()
            self.fp = samopen( samfile._filename, mode, NULL )
            store.release()
            assert self.fp != NULL
            self.owns_samfile = True
        else:
            self.fp = samfile.samfile
            self.owns_samfile = False

        # allocate memory for alignment
        self.b = <bam1_t*>calloc(1, sizeof(bam1_t))

        self.positions = positions
        self.current_pos = 0

    def __iter__(self):
        return self 

    cdef bam1_t * getCurrent( self ):
        return self.b

    cdef int cnext(self):
        '''cversion of iterator'''

        # end iteration if out of positions
        if self.current_pos >= len(self.positions): return -1

        bam_seek( self.fp.x.bam, self.positions[self.current_pos], 0 ) 
        self.current_pos += 1
        return samread(self.fp, self.b)

    def __next__(self): 
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """

        cdef int ret = self.cnext()
        if (ret > 0):
            return makeAlignedRead( self.b )
        else:
            raise StopIteration

    def __dealloc__(self):
        bam_destroy1(self.b)
        if self.owns_samfile: samclose( self.fp )

##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
ctypedef struct __iterdata:
    samfile_t * samfile
    bam_iter_t iter
    faidx_t * fastafile
    int tid
    char * seq
    int seq_len

cdef int __advance_all( void * data, bam1_t * b ):
    '''advance without any read filtering.
    '''
    cdef __iterdata * d
    d = <__iterdata*>data
    return bam_iter_read( d.samfile.x.bam, d.iter, b )

cdef int __advance_snpcalls( void * data, bam1_t * b ):
    '''advance using same filter and read processing as in 
    the samtools pileup.
    '''
    cdef __iterdata * d
    d = <__iterdata*>data

    cdef int ret = bam_iter_read( d.samfile.x.bam, d.iter, b )
    cdef int skip = 0
    cdef int q
    cdef int is_cns = 1
    cdef int is_nobaq = 0
    cdef int capQ_thres = 0

    # reload sequence
    if d.fastafile != NULL and b.core.tid != d.tid:
        if d.seq != NULL: free(d.seq)
        d.tid = b.core.tid
        d.seq = faidx_fetch_seq(d.fastafile, 
                                d.samfile.header.target_name[d.tid],
                                0, max_pos, 
                                &d.seq_len)
        if d.seq == NULL:
            raise ValueError( "reference sequence for '%s' (tid=%i) not found" % \
                                  (d.samfile.header.target_name[d.tid], 
                                   d.tid))


    while ret >= 0:

        skip = 0

        # realign read - changes base qualities
        if d.seq != NULL and is_cns and not is_nobaq: bam_prob_realn( b, d.seq )

        if d.seq != NULL and capQ_thres > 10:
            q = bam_cap_mapQ(b, d.seq, capQ_thres)
            if q < 0: skip = 1
            elif b.core.qual > q: b.core.qual = q
        if b.core.flag & BAM_FUNMAP: skip = 1
        elif b.core.flag & 1 and not b.core.flag & 2: skip = 1

        if not skip: break
        # additional filters

        ret = bam_iter_read( d.samfile.x.bam, d.iter, b )

    return ret

cdef class IteratorColumn:
    '''abstract base class for iterators over columns.

    IteratorColumn objects wrap the pileup functionality of samtools.
    
    For reasons of efficiency, the iterator points to the current 
    pileup buffer. The pileup buffer is updated at every iteration.
    This might cause some unexpected behavious. For example,
    consider the conversion to a list::
    
       f = Samfile("file.bam", "rb")
       result = list( f.pileup() )

    Here, ``result`` will contain ``n`` objects of type :class:`PileupProxy` for ``n`` columns, 
    but each object in ``result`` will contain the same information.
    
    The desired behaviour can be achieved by list comprehension::

       result = [ x.pileups() for x in f.pileup() ]

    ``result`` will be a list of ``n`` lists of objects of type :class:`PileupRead`.

    If the iterator is associated with a :class:`Fastafile` using the :meth:`addReference`
    method, then the iterator will export the current sequence via the methods :meth:`getSequence`
    and :meth:`seq_len`.

    Optional kwargs to the iterator

    stepper
       The stepper controlls how the iterator advances. 
       Possible options for the stepper are 
       
       all
           use all reads for pileup.
       samtools
           same filter and read processing as in :term:`csamtools` pileup
    fastafile
       A :class:`FastaFile` object
    mask
       Skip all reads with bits set in mask.
    max_depth
       maximum read depth. The default is 8000.
    '''

    # result of the last plbuf_push
    cdef IteratorRowRegion iter
    cdef int tid
    cdef int pos
    cdef int n_plp
    cdef int mask
    cdef const_bam_pileup1_t_ptr plp
    cdef bam_plp_t pileup_iter
    cdef __iterdata iterdata 
    cdef Samfile samfile
    cdef Fastafile fastafile
    cdef stepper
    cdef int max_depth

    def __cinit__( self, Samfile samfile, **kwargs ):
        self.samfile = samfile
        self.mask = kwargs.get("mask", BAM_DEF_MASK )
        self.fastafile = kwargs.get( "fastafile", None )
        self.stepper = kwargs.get( "stepper", None )
        self.max_depth = kwargs.get( "max_depth", 8000 )
        self.iterdata.seq = NULL
        self.tid = 0
        self.pos = 0
        self.n_plp = 0
        self.plp = NULL
        self.pileup_iter = <bam_plp_t>NULL


    def __iter__(self):
        return self 

    cdef int cnext(self):
        '''perform next iteration.

        This method is analogous to the samtools bam_plp_auto method.
        It has been re-implemented to permit for filtering.
        '''
        self.plp = bam_plp_auto( self.pileup_iter, 
                                 &self.tid,
                                 &self.pos,
                                 &self.n_plp )

    cdef char * getSequence( self ):
        '''return current reference sequence underlying the iterator.
        '''
        return self.iterdata.seq

    property seq_len:
        '''current sequence length.'''
        def __get__(self): return self.iterdata.seq_len

    def addReference( self, Fastafile fastafile ):
       '''
       add reference sequences in *fastafile* to iterator.'''
       self.fastafile = fastafile
       if self.iterdata.seq != NULL: free(self.iterdata.seq)
       self.iterdata.tid = -1
       self.iterdata.fastafile = self.fastafile.fastafile

    def hasReference( self ):
        '''
        return true if iterator is associated with a reference'''
        return self.fastafile

    cdef setMask( self, mask ):
        '''set masking flag in iterator.

        reads with bits set in *mask* will be skipped.
        '''
        self.mask = mask
        bam_plp_set_mask( self.pileup_iter, self.mask )

    cdef setupIteratorData( self, 
                            int tid, 
                            int start, 
                            int end, 
                            int reopen = 0 ):
        '''setup the iterator structure'''

        self.iter = IteratorRowRegion( self.samfile, tid, start, end, reopen )
        self.iterdata.samfile = self.samfile.samfile
        self.iterdata.iter = self.iter.iter
        self.iterdata.seq = NULL
        self.iterdata.tid = -1

        if self.fastafile != None:
            self.iterdata.fastafile = self.fastafile.fastafile
        else:
            self.iterdata.fastafile = NULL

        if self.stepper == None or self.stepper == "all":
            self.pileup_iter = bam_plp_init( &__advance_all, &self.iterdata )
        elif self.stepper == "samtools":        
            self.pileup_iter = bam_plp_init( &__advance_snpcalls, &self.iterdata )
        else:
            raise ValueError( "unknown stepper option `%s` in IteratorColumn" % self.stepper)

        if self.max_depth:
            bam_plp_set_maxcnt( self.pileup_iter, self.max_depth )

        bam_plp_set_mask( self.pileup_iter, self.mask )

    cdef reset( self, tid, start, end ):
        '''reset iterator position.

        This permits using the iterator multiple times without
        having to incur the full set-up costs.
        '''
        self.iter = IteratorRowRegion( self.samfile, tid, start, end, reopen = 0 )
        self.iterdata.iter = self.iter.iter
       
        # invalidate sequence if different tid
        if self.tid != tid:
            if self.iterdata.seq != NULL: free( self.iterdata.seq )
            self.iterdata.seq = NULL            
            self.iterdata.tid = -1
            
        # self.pileup_iter = bam_plp_init( &__advancepileup, &self.iterdata )
        bam_plp_reset(self.pileup_iter)

    def __dealloc__(self):
        # reset in order to avoid memory leak messages for iterators that have
        # not been fully consumed
        if self.pileup_iter != <bam_plp_t>NULL:
            bam_plp_reset(self.pileup_iter)
            bam_plp_destroy(self.pileup_iter)
            self.pileup_iter = <bam_plp_t>NULL

        if self.iterdata.seq != NULL: 
            free(self.iterdata.seq)
            self.iterdata.seq = NULL

cdef class IteratorColumnRegion(IteratorColumn):
    '''iterates over a region only.
    '''
    def __cinit__(self, Samfile samfile, 
                  int tid = 0, 
                  int start = 0, 
                  int end = max_pos,
                  **kwargs ):

        # initialize iterator
        self.setupIteratorData( tid, start, end, 1 )

    def __next__(self): 
        """python version of next().
        """

        while 1:
            self.cnext()
            if self.n_plp < 0:
                raise ValueError("error during iteration" )
        
            if self.plp == NULL:
                raise StopIteration
            
            return makePileupProxy( <bam_pileup1_t*>self.plp, 
                                     self.tid, 
                                     self.pos, 
                                     self.n_plp )

cdef class IteratorColumnAllRefs(IteratorColumn):
    """iterates over all columns by chaining iterators over each reference
    """

    def __cinit__(self, 
                  Samfile samfile,
                  **kwargs ):

        # no iteration over empty files
        if not samfile.nreferences: raise StopIteration

        # initialize iterator
        self.setupIteratorData( self.tid, 0, max_pos, 1 )

    def __next__(self):
        """python version of next().
        """
        
        while 1:
            self.cnext()

            if self.n_plp < 0:
                raise ValueError("error during iteration" )
        
            # return result, if within same reference
            if self.plp != NULL:
                return makePileupProxy( <bam_pileup1_t*>self.plp, 
                                         self.tid, 
                                         self.pos, 
                                         self.n_plp )

            # otherwise, proceed to next reference or stop
            self.tid += 1
            if self.tid < self.samfile.nreferences:
                self.setupIteratorData( self.tid, 0, max_pos, 0 )
            else:
                raise StopIteration

##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
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

##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
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
    Class representing an aligned read. see SAM format specification for 
    the meaning of fields (http://samtools.sourceforge.net/).

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
        """return string representation of alignment.

        The representation is an approximate :term:`sam` format.

        An aligned read might not be associated with a :term:`Samfile`.
        As a result :term:`tid` is shown instead of the reference name.

        Similarly, the tags field is returned in its parsed state.
        """
        # sam-parsing is done in sam.c/bam_format1_core which
        # requires a valid header.
        return "\t".join(map(str, (self.qname,
                                   self.flag,
                                   self.rname,
                                   self.pos,
                                   self.mapq,
                                   self.cigar,
                                   self.mrnm,
                                   self.mpos,
                                   self.rlen,
                                   self.seq,
                                   self.qual,
                                   self.tags )))
       
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


        This method will happily write the same tag
        multiple times.
        """
        def __get__(self):
            cdef char * ctag
            cdef bam1_t * src
            cdef uint8_t * s
            cdef char auxtag[3]
            cdef char auxtype
            
            src = self._delegate
            if src.l_aux == 0: return []
            
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
                        fmt, pytype = "<cccf", 'f'
                    elif t == types.IntType:
                        if value < 0:
                            if value >= -127: fmt, pytype = "<cccb", 'c'
                            elif value >= -32767: fmt, pytype = "<ccch", 's'
                            elif value < -2147483648: raise ValueError( "integer %i out of range of BAM/SAM specification" % value )
                            else: fmt, pytype = "<ccci", 'i'[0]
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

        DEPRECATED from pysam-0.4 - use tid in the future.
        The rname field caused a lot of confusion as it returns
        the :term:`target` ID instead of the reference sequence
        name.

        .. note::

            This field contains the index of the reference sequence 
            in the sequence dictionary. To obtain the name
            of the reference sequence, use :meth:`pysam.Samfile.getrname()`
            
        """
        def __get__(self): return self._delegate.core.tid
        def __set__(self, tid): self._delegate.core.tid = tid

    property tid: 
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
        '''aligned end position of the read on the reference genome.  Returns
        None if not available.'''
        def __get__(self):
            cdef bam1_t * src
            src = self._delegate
            if (self.flag & BAM_FUNMAP) or src.core.n_cigar == 0:
                return None
            return bam_calend(&src.core, bam1_cigar(src))
    property alen:
        '''aligned length of the read on the reference genome.  Returns None if
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
        """the :term:`reference` id of the mate 
        deprecated, use RNEXT instead.
        """     
        def __get__(self): return self._delegate.core.mtid
        def __set__(self, mtid): self._delegate.core.mtid = mtid
    property rnext:
        """the :term:`reference` id of the mate """     
        def __get__(self): return self._delegate.core.mtid
        def __set__(self, mtid): self._delegate.core.mtid = mtid
    property mpos: 
        """the position of the mate
        deprecated, use PNEXT instead."""
        def __get__(self): return self._delegate.core.mpos
        def __set__(self, mpos): self._delegate.core.mpos = mpos
    property pnext: 
        """the position of the mate"""
        def __get__(self): return self._delegate.core.mpos
        def __set__(self, mpos): self._delegate.core.mpos = mpos
    property isize: 
        """the insert size
        deprecated: use tlen instead"""
        def __get__(self): return self._delegate.core.isize
        def __set__(self, isize): self._delegate.core.isize = isize
    property tlen: 
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
        def __get__(self): return (self.flag & BAM_FREVERSE) != 0
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
        """true if optical or PCR duplicate"""
        def __get__(self): return (self.flag & BAM_FDUP) != 0
        def __set__(self,val): 
            if val: self._delegate.core.flag |= BAM_FDUP
            else: self._delegate.core.flag &= ~BAM_FDUP
    property positions:
        """a list of reference positions that this read aligns to."""
        def __get__(self):
            cdef uint32_t k, i, pos
            cdef int op
            cdef uint32_t * cigar_p
            cdef bam1_t * src 

            result = []
            src = self._delegate
            if src.core.n_cigar == 0: return []

            pos = src.core.pos

            cigar_p = bam1_cigar(src)
            for k from 0 <= k < src.core.n_cigar:
                op = cigar_p[k] & BAM_CIGAR_MASK
                l = cigar_p[k] >> BAM_CIGAR_SHIFT
                if op == BAM_CMATCH:
                    for i from pos <= i < pos + l:
                        result.append( i )

                if op == BAM_CMATCH or op == BAM_CDEL or op == BAM_CREF_SKIP:
                    pos += l

            return result

    def overlap( self, uint32_t start, uint32_t end ):
        """return number of aligned bases of read overlapping the interval *start* and *end*
        on the reference sequence.
        """
        cdef uint32_t k, i, pos, overlap
        cdef int op, o
        cdef uint32_t * cigar_p
        cdef bam1_t * src 

        overlap = 0

        src = self._delegate
        if src.core.n_cigar == 0: return 0
        pos = src.core.pos
        o = 0

        cigar_p = bam1_cigar(src)
        for k from 0 <= k < src.core.n_cigar:
            op = cigar_p[k] & BAM_CIGAR_MASK
            l = cigar_p[k] >> BAM_CIGAR_SHIFT

            if op == BAM_CMATCH:
                o = min( pos + l, end) - max( pos, start )
                if o > 0: overlap += o

            if op == BAM_CMATCH or op == BAM_CDEL or op == BAM_CREF_SKIP:
                pos += l

        return overlap

    def opt(self, tag):
        """retrieves optional data given a two-letter *tag*"""
        #see bam_aux.c: bam_aux_get() and bam_aux2i() etc 
        cdef uint8_t * v
        v = bam_aux_get(self._delegate, tag)
        if v == NULL: raise KeyError( "tag '%s' not present" % tag )
        type = chr(v[0])
        if type == 'c' or type == 'C' or type == 's' or type == 'S':
            return <int>bam_aux2i(v)            
        elif type == 'i' or type == 'I':
            return <int32_t>bam_aux2i(v)            
        elif type == 'f' or type == 'F':
            return <float>bam_aux2f(v)
        elif type == 'd' or type == 'D':
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
        sys.stderr.flush()          #  Buffered data goes to old stream.
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

def _samtools_dispatch( method, 
                        args = (),
                        catch_stdout = True,
                        catch_stderr = False,
                        ):
    '''call ``method`` in samtools providing arguments in args.
    
    .. note:: 
       This method redirects stdout and (optionally) stderr to capture it 
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
    # as stdout/stderr will not appear on the terminal
    
    # some special cases
    if method == "index":
        if not os.path.exists( args[0] ):
            raise IOError( "No such file or directory: '%s'" % args[0] )

    # redirect stderr and stdout to file
    if catch_stderr:
        stderr_h, stderr_f = tempfile.mkstemp()
        stderr_save = Outs( sys.stderr.fileno() )
        stderr_save.setfd( stderr_h )

    if catch_stdout:
        stdout_h, stdout_f = tempfile.mkstemp()
        stdout_save = Outs( sys.stdout.fileno() )
        stdout_save.setfd( stdout_h )

        # patch for `samtools view`
        # samtools `view` closes stdout, from which I can not
        # recover. Thus redirect output to file with -o option.
        if method == "view":
            if "-o" in args: raise ValueError("option -o is forbidden in samtools view")
            args = ( "-o", stdout_f ) + args

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
    if catch_stdout:
        stdout_save.restore()
        out_stdout = open( stdout_f, "r").readlines()
        os.remove( stdout_f )
    else:
        out_stdout = []

    if catch_stderr:
        stderr_save.restore()
        out_stderr = open( stderr_f, "r").readlines()
        os.remove( stderr_f )
    else:
        out_stderr = []
    
    return retval, out_stderr, out_stdout

cdef class SNPCall:
    '''the results of a SNP call.'''
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


# cdef class SNPCallerBase:
#     '''Base class for SNP callers.

#     *min_baseQ*
#        minimum base quality (possibly capped by BAQ)
#     *capQ_threshold*
#        coefficient for adjusting mapQ of poor mappings
#     *theta*
#        theta in maq consensus calling model
#     *n_haplotypes*
#        number of haplotypes in the sample
#     *het_rate*
#        prior of a difference between two haplotypes
#     '''

#     cdef bam_maqcns_t * c
#     cdef IteratorColumn iter

#     def __cinit__(self, 
#                   IteratorColumn iterator_column, 
#                   **kwargs ):

#         self.iter = iterator_column
#         self.c =  bam_maqcns_init()

#         # set the default parameterization according to
#         # samtools

#         # new default mode for samtools >0.1.10
#         self.c.errmod = kwargs.get( "errmod", BAM_ERRMOD_MAQ2 )

#         self.c.min_baseQ = kwargs.get( "min_baseQ", 13 )
#         # self.c.capQ_thres = kwargs.get( "capQ_threshold", 60 )
#         self.c.n_hap = kwargs.get( "n_haplotypes", 2 )
#         self.c.het_rate = kwargs.get( "het_rate", 0.001 )
#         self.c.theta = kwargs.get( "theta", 0.83 )

#         if self.c.errmod != BAM_ERRMOD_MAQ2:
#             self.c.theta += 0.02

#         # call prepare AFTER setting parameters
#         bam_maqcns_prepare( self.c )

#     def __dealloc__(self):
#         bam_maqcns_destroy( self.c )

    # cdef __dump( self, glf1_t * g, uint32_t cns, int rb ):
    #     '''debugging output.'''

    #     pysam_dump_glf( g, self.c );
    #     print ""
    #     for x in range(self.iter.n_plp):
    #         print "--> read %i %s %i" % (x, 
    #                                      bam1_qname(self.iter.plp[x].b),
    #                                      self.iter.plp[x].qpos,
    #                                      )

    #     print "pos=%i, cns=%i, q_r = %f, depth=%i, n=%i, rb=%i, cns-cq=%i %i %i %i" \
    #         % (self.iter.pos, 
    #            cns, 
    #            self.c.q_r,
    #            self.iter.n_plp,
    #            self.iter.n_plp,
    #            rb,
    #            cns >> 8 & 0xff,
    #            cns >> 16 & 0xff,
    #            cns & 0xff,
    #            cns >> 28,
    #            )

    #     printf("-------------------------------------\n");
    #     sys.stdout.flush()

# cdef class IteratorSNPCalls( SNPCallerBase ):
#     """*(IteratorColumn iterator)*

#     call SNPs within a region.

#     *iterator* is a pileup iterator. SNPs will be called
#     on all positions returned by this iterator.

#     This caller is fast if SNPs are called over large continuous
#     regions. It is slow, if instantiated frequently and in random
#     order as the sequence will have to be reloaded.

#     """
    
#     def __cinit__(self, 
#                   IteratorColumn iterator_column,
#                   **kwargs ):

#         assert self.iter.hasReference(), "IteratorSNPCalls requires an pileup iterator with reference sequence"

#     def __iter__(self):
#         return self 

#     def __next__(self): 
#         """python version of next().
#         """

#         # the following code was adapted from bam_plcmd.c:pileup_func()
#         self.iter.cnext()

#         if self.iter.n_plp < 0:
#             raise ValueError("error during iteration" )

#         if self.iter.plp == NULL:
#            raise StopIteration

#         cdef char * seq = self.iter.getSequence()
#         cdef int seq_len = self.iter.seq_len

#         assert seq != NULL

#         # reference base
#         if self.iter.pos >= seq_len:
#             raise ValueError( "position %i out of bounds on reference sequence (len=%i)" % (self.iter.pos, seq_len) )

#         cdef int rb = seq[self.iter.pos]
#         cdef uint32_t cns 
#        cdef glf1_t * g

#        g = bam_maqcns_glfgen( self.iter.n_plp,
#                               self.iter.plp,
#                               bam_nt16_table[rb],
#                               self.c )

#        if pysam_glf_depth( g ) == 0:
#            cns = 0xfu << 28 | 0xf << 24
#        else:
#            cns = glf2cns(g, <int>(self.c.q_r + .499))
           
#        free(g)
            
#         cdef SNPCall call

#         call = SNPCall()
#         call._tid = self.iter.tid
#         call._pos = self.iter.pos
#         call._reference_base = rb
#         call._genotype = bam_nt16_rev_table[cns>>28]
#         call._consensus_quality = cns >> 8 & 0xff
#         call._snp_quality = cns & 0xff
#         call._rms_mapping_quality = cns >> 16&0xff
#         call._coverage = self.iter.n_plp

#         return call 

# cdef class SNPCaller( SNPCallerBase ):
#     '''*(IteratorColumn iterator_column )*

#     The samtools SNP caller.

#     This object will call SNPs in *samfile* against the reference
#     sequence in *fasta*.

#     This caller is fast for calling few SNPs in selected regions.

#     It is slow, if called over large genomic regions.
#     '''


#     def __cinit__(self, 
#                   IteratorColumn iterator_column, 
#                   **kwargs ):

#         pass

#     def call(self, reference, int pos ): 
#         """call a snp on chromosome *reference*
#         and position *pos*.

#         returns a :class:`SNPCall` object.
#         """

#         cdef int tid = self.iter.samfile.gettid( reference )

#         self.iter.reset( tid, pos, pos + 1 )

#         while 1:
#             self.iter.cnext()
            
#             if self.iter.n_plp < 0:
#                 raise ValueError("error during iteration" )

#             if self.iter.plp == NULL:
#                 raise ValueError( "no reads in region - no call" )
             
#             if self.iter.pos == pos: break

#         cdef char * seq = self.iter.getSequence()
#         cdef int seq_len = self.iter.seq_len

#         assert seq != NULL

#         # reference base
#         if self.iter.pos >= seq_len:
#             raise ValueError( "position %i out of bounds on reference sequence (len=%i)" % (self.iter.pos, seq_len) )

#         cdef int rb = seq[self.iter.pos]
#         cdef uint32_t cns 
# #        cdef glf1_t * g
# #
# #        g = bam_maqcns_glfgen( self.iter.n_plp,
# #                               self.iter.plp,
# #                               bam_nt16_table[rb],
# #                               self.c )
# ##
# #
# #        if pysam_glf_depth( g ) == 0:
# #            cns = 0xfu << 28 | 0xf << 24
# #        else:
# #            cns = glf2cns(g, <int>(self.c.q_r + .499))
# #
# #        free(g)
            
#         cdef SNPCall call

#         call = SNPCall()
#         call._tid = self.iter.tid
#         call._pos = self.iter.pos
#         call._reference_base = rb
#         call._genotype = bam_nt16_rev_table[cns>>28]
#         call._consensus_quality = cns >> 8 & 0xff
#         call._snp_quality = cns & 0xff
#         call._rms_mapping_quality = cns >> 16&0xff
#         call._coverage = self.iter.n_plp

#         return call 

# cdef class IndelCall:
#     '''the results of an indel call.'''
#     cdef int _tid
#     cdef int _pos
#     cdef int _coverage
#     cdef int _rms_mapping_quality
#     cdef bam_maqindel_ret_t * _r 

#     def __cinit__(self):
#         #assert r != NULL
#         #self._r = r
#         pass

#     property tid:
#         '''the chromosome ID as is defined in the header'''
#         def __get__(self): 
#             return self._tid
    
#     property pos:
#        '''nucleotide position of SNP.'''
#        def __get__(self): return self._pos

#     property genotype:
#        '''the genotype called.'''
#        def __get__(self): 
#            if self._r.gt == 0:
#                s = PyString_FromStringAndSize( self._r.s[0], self._r.indel1 + 1)
#                return "%s/%s" % (s,s)
#            elif self._r.gt == 1:
#                s = PyString_FromStringAndSize( self._r.s[1], self._r.indel2 + 1)
#                return "%s/%s" % (s,s)
#            else:
#                return "%s/%s" % (self.first_allele, self.second_allele )

#     property consensus_quality:
#        '''the genotype quality (Phred-scaled).'''
#        def __get__(self): return self._r.q_cns

#     property snp_quality:
#        '''the snp quality (Phred scaled) - probability of consensus being identical to reference sequence.'''
#        def __get__(self): return self._r.q_ref

#     property mapping_quality:
#        '''the root mean square (rms) of the mapping quality of all reads involved in the call.'''
#        def __get__(self): return self._rms_mapping_quality

#     property coverage:
#        '''coverage or read depth - the number of reads involved in the call.'''
#        def __get__(self): return self._coverage

#     property first_allele:
#        '''sequence of first allele.'''
#        def __get__(self): return PyString_FromStringAndSize( self._r.s[0], self._r.indel1 + 1)

#     property second_allele:
#        '''sequence of second allele.'''
#        def __get__(self): return PyString_FromStringAndSize( self._r.s[1], self._r.indel2 + 1)

#     property reads_first:
#        '''reads supporting first allele.'''
#        def __get__(self): return self._r.cnt1

#     property reads_second:
#        '''reads supporting first allele.'''
#        def __get__(self): return self._r.cnt2

#     property reads_diff:
#        '''reads supporting first allele.'''
#        def __get__(self): return self._r.cnt_anti

#     def __str__(self):

#         return "\t".join( map(str, (
#                     self.tid,
#                     self.pos,
#                     self.genotype,
#                     self.consensus_quality,
#                     self.snp_quality,
#                     self.mapping_quality,
#                     self.coverage,
#                     self.first_allele,
#                     self.second_allele,
#                     self.reads_first,
#                     self.reads_second,
#                     self.reads_diff ) ) )

#     def __dealloc__(self ):
#         bam_maqindel_ret_destroy(self._r)

# cdef class IndelCallerBase:
#     '''Base class for SNP callers.

#     *min_baseQ*
#        minimum base quality (possibly capped by BAQ)
#     *capQ_threshold*
#        coefficient for adjusting mapQ of poor mappings
#     *theta*
#        theta in maq consensus calling model
#     *n_haplotypes*
#        number of haplotypes in the sample
#     *het_rate*
#        prior of a difference between two haplotypes
#     '''

#     cdef bam_maqindel_opt_t * options
#     cdef IteratorColumn iter
#     cdef int cap_mapQ
#     cdef int max_depth

#     def __cinit__(self, 
#                   IteratorColumn iterator_column, 
#                   **kwargs ):


#         self.iter = iterator_column

#         assert iterator_column.hasReference(), "IndelCallerBase requires an pileup iterator with reference sequence"

#         self.options = bam_maqindel_opt_init()

#         # set the default parameterization according to
#         # samtools

#         self.options.r_indel = kwargs.get( "r_indel", 0.00015 )
#         self.options.q_indel = kwargs.get( "q_indel", 40 )
#         self.cap_mapQ = kwargs.get( "cap_mapQ", 60 )
#         self.max_depth = kwargs.get( "max_depth", 1024 )

#     def __dealloc__(self):
#         free( self.options )

#     def _call( self ):

#         cdef char * seq = self.iter.getSequence()
#         cdef int seq_len = self.iter.seq_len

#         assert seq != NULL

#         # reference base
#         if self.iter.pos >= seq_len:
#             raise ValueError( "position %i out of bounds on reference sequence (len=%i)" % (self.iter.pos, seq_len) )

#         cdef bam_maqindel_ret_t * r 
        
#         cdef int m = min( self.max_depth, self.iter.n_plp )

#         # printf("pysam: m=%i, q_indel=%i, r_indel=%f, r_snp=%i, mm_penalty=%i, indel_err=%i, ambi_thres=%i\n",
#         #        m, self.options.q_indel, self.options.r_indel, self.options.r_snp, self.options.mm_penalty,
#         #        self.options.indel_err, self.options.ambi_thres );

#         r = bam_maqindel(m, 
#                          self.iter.pos, 
#                          self.options,
#                          self.iter.plp, 
#                          seq,
#                          0, 
#                          NULL)
        
#         if r == NULL: return None

#         cdef IndelCall call
#         call = IndelCall()
#         call._r = r
#         call._tid = self.iter.tid
#         call._pos = self.iter.pos
#         call._coverage = self.iter.n_plp

#         cdef uint64_t rms_aux = 0
#         cdef int i = 0
#         cdef bam_pileup1_t * p
#         cdef int tmp

#         for i from 0 <= i < self.iter.n_plp:
#             p = self.iter.plp + i
#             if p.b.core.qual < self.cap_mapQ:
#                 tmp = p.b.core.qual 
#             else:
#                 tmp = self.cap_mapQ
#             rms_aux += tmp * tmp

#         call._rms_mapping_quality = <uint64_t>(sqrt(<double>rms_aux / self.iter.n_plp) + .499)

#         return call 

# cdef class IndelCaller( IndelCallerBase ):
#     '''*(IteratorColumn iterator_column )*

#     The samtools SNP caller.

#     This object will call SNPs in *samfile* against the reference
#     sequence in *fasta*.

#     This caller is fast for calling few SNPs in selected regions.

#     It is slow, if called over large genomic regions.
#     '''

#     def __cinit__(self, 
#                   IteratorColumn iterator_column, 
#                   **kwargs ):

#         pass

#     def call(self, reference, int pos ): 
#         """call a snp on chromosome *reference*
#         and position *pos*.

#         returns a :class:`SNPCall` object or None, if no indel call could be made.
#         """

#         cdef int tid = self.iter.samfile.gettid( reference )

#         self.iter.reset( tid, pos, pos + 1 )

#         while 1:
#             self.iter.cnext()
            
#             if self.iter.n_plp < 0:
#                 raise ValueError("error during iteration" )

#             if self.iter.plp == NULL:
#                 raise ValueError( "no reads in region - no call" )
             
#             if self.iter.pos == pos: break

#         return self._call()

# cdef class IteratorIndelCalls( IndelCallerBase ):
#     """*(IteratorColumn iterator)*

#     call indels within a region.

#     *iterator* is a pileup iterator. SNPs will be called
#     on all positions returned by this iterator.

#     This caller is fast if SNPs are called over large continuous
#     regions. It is slow, if instantiated frequently and in random
#     order as the sequence will have to be reloaded.

#     """
    
#     def __cinit__(self, 
#                   IteratorColumn iterator_column,
#                   **kwargs ):
#         pass


#     def __iter__(self):
#         return self 

#     def __next__(self): 
#         """python version of next().
#         """

#         # the following code was adapted from bam_plcmd.c:pileup_func()
#         self.iter.cnext()

#         if self.iter.n_plp < 0:
#             raise ValueError("error during iteration" )

#         if self.iter.plp == NULL:
#            raise StopIteration

#         return self._call()



cdef class IndexedReads:
    """index a bamfile by read.

    The index is kept in memory.

    By default, the file is re-openend to avoid conflicts if
    multiple operators work on the same file. Set *reopen* = False
    to not re-open *samfile*.
    """

    cdef Samfile samfile
    cdef samfile_t * fp
    cdef index
    # true if samfile belongs to this object
    cdef int owns_samfile

    def __init__(self, Samfile samfile, int reopen = True ):
        self.samfile = samfile

        if samfile.isbam: mode = "rb"
        else: mode = "r"

        # reopen the file - note that this makes the iterator
        # slow and causes pileup to slow down significantly.
        if reopen:
            store = StderrStore()            
            self.fp = samopen( samfile._filename, mode, NULL )
            store.release()
            assert self.fp != NULL
            self.owns_samfile = True
        else:
            self.fp = samfile.samfile
            self.owns_samfile = False

        assert samfile.isbam, "can only IndexReads on bam files"

    def build( self ):
        '''build index.'''
        
        self.index = collections.defaultdict( list )

        # this method will start indexing from the current file position
        # if you decide
        cdef int ret = 1
        cdef bam1_t * b = <bam1_t*> calloc(1, sizeof( bam1_t) )
        
        cdef uint64_t pos

        while ret > 0:
            pos = bam_tell( self.fp.x.bam ) 
            ret = samread( self.fp, b)
            if ret > 0:
                qname = bam1_qname( b )
                self.index[qname].append( pos )                
            
        bam_destroy1( b )

    def find( self, qname ):
        if qname in self.index:
            return IteratorRowSelection( self.samfile, self.index[qname], reopen = False )
        else:
            raise KeyError( "read %s not found" % qname )

    def __dealloc__(self):
        if self.owns_samfile: samclose( self.fp )

__all__ = ["Samfile", 
           "Fastafile",
           "IteratorRow", 
           "IteratorColumn", 
           "AlignedRead", 
           "PileupColumn", 
           "PileupProxy", 
           "PileupRead",
           # "IteratorSNPCalls",
           # "SNPCaller",
           # "IndelCaller",
           # "IteratorIndelCalls", 
           "IndexedReads" ]

               

