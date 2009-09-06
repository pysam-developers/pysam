from csamtools import *

class SamtoolsError( Exception ):
    '''exception raised in case of an error incurred in the samtools library.'''

    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

# list of samtools command line options to export in python
SAMTOOLS_DISPATCH = ( "view",
                      "sort",
                      "pileup",
                      "faidx",
                      "tview",
                      "index",
                      "fixmate",
                      "glfview",
                      "flagstat",
                      "calmd",
                      "merge",  
                      "rmdup" )

class SamtoolsDispatcher(object):
    '''samtools dispatcher. 

    Emulates the samtools command line.
    '''
    dispatch=None

    def __init__(self,dispatch): 
        self.dispatch = dispatch

    def __call__(self,*args,**kwargs):
        retval, stderr = csamtools._samtools_dispatch( self.dispatch, args )
        if retval: raise SamtoolsError( "\n".join( stderr ) )
        return retval

    def usage(self):
        '''return the samtools usage information for this command'''
        retval, stderr = csamtools._samtools_dispatch( self.dispatch )
        return "".join(stderr)

# instantiate samtools commands as python functions
for key in SAMTOOLS_DISPATCH:
    globals()[key] = SamtoolsDispatcher(key)

# hack to export all the symbols from csamtools
__all__ = csamtools.__all__ + [ "SamtoolsError", "SamtoolsDispatcher",] + list(SAMTOOLS_DISPATCH)

