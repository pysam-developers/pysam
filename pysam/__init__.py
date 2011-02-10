from csamtools import *
from ctabix import *
import csamtools
import ctabix
import Pileup
import sys
import os

class SamtoolsError( Exception ):
    '''exception raised in case of an error incurred in the samtools library.'''

    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class SamtoolsDispatcher(object):
    '''samtools dispatcher. 

    Emulates the samtools command line as module calls.
    
    Captures stdout and stderr. 

    Raises a :class:`pysam.SamtoolsError` exception in case
    samtools exits with an error code other than 0.

    Some command line options are associated with parsers.
    For example, the samtools command "pileup -c" creates
    a tab-separated table on standard output. In order to 
    associate parsers with options, an optional list of 
    parsers can be supplied. The list will be processed
    in order checking for the presence of each option.

    If no parser is given or no appropriate parser is found, 
    the stdout output of samtools commands will be returned.
    '''
    dispatch=None
    parsers=None

    def __init__(self,dispatch, parsers): 
        self.dispatch = dispatch
        self.parsers = parsers
        self.stderr = []

    def __call__(self, *args, **kwargs):
        '''execute a samtools command
        '''
        retval, stderr, stdout = csamtools._samtools_dispatch( self.dispatch, args )
        if retval: raise SamtoolsError( "\n".join( stderr ) )
        self.stderr = stderr
        # samtools commands do not propagate the return code correctly.
        # I have thus added this patch to throw if there is output on stderr.
        # Note that there is sometimes output on stderr that is not an error,
        # for example: [sam_header_read2] 2 sequences loaded.
        # Ignore messages like these
        stderr = [ x for x in stderr \
                       if not (x.startswith( "[sam_header_read2]" ) or \
                                   x.startswith("[bam_index_load]") or \
                                   x.startswith("[bam_sort_core]") or \
                                   x.startswith("[samopen] SAM header is present") )
                   ]
        if stderr: raise SamtoolsError( "\n".join( stderr ) )

        # call parser for stdout:
        if not kwargs.get("raw") and stdout and self.parsers:
            for options, parser in self.parsers:
                for option in options: 
                    if option not in args: break
                else:
                    return parser(stdout)

        return stdout

    def getMessages( self ):
        return self.stderr

    def usage(self):
        '''return the samtools usage information for this command'''
        retval, stderr, stdout = csamtools._samtools_dispatch( self.dispatch )
        return "".join(stderr)

#
# samtools command line options to export in python
#
# import is a python reserved word.
SAMTOOLS_DISPATCH = { 
    "view" : ( "view", None ),
    "sort" : ( "sort", None),
    "samimport": ( "import", None),
    "pileup" : ( "pileup", ( (("-c",), Pileup.iterate ), ), ),
    "faidx" : ("faidx", None),
    "tview" : ("tview", None),
    "index" : ("index", None),
    "fixmate" : ("fixmate", None),
    "glfview" : ("glfview", None),
    "flagstat" : ("flagstat", None),
    "calmd" : ("calmd", None),
    "merge" : ("merge", None),  
    "rmdup" : ("rmdup", None) }

# instantiate samtools commands as python functions
for key, options in SAMTOOLS_DISPATCH.iteritems():
    cmd, parser = options
    globals()[key] = SamtoolsDispatcher(cmd, parser)

# hack to export all the symbols from csamtools
__all__ = csamtools.__all__ + \
    ctabix.__all__ + \
    [ "SamtoolsError", "SamtoolsDispatcher" ] + list(SAMTOOLS_DISPATCH) +\
    ["Pileup",] 

from version import __version__, __samtools_version__
