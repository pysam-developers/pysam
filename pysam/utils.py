from typing import (
    Callable,
    List,
    Tuple,
    Iterable,
    Union,
)

from pysam.libcutils import _pysam_dispatch


class SamtoolsError(Exception):
    '''exception raised in case of an error incurred in the samtools
    library.'''

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class PysamDispatcher(object):
    '''The dispatcher emulates the samtools/bctools command line.

    Captures stdout and stderr.

    Raises a :class:`pysam.SamtoolsError` exception in case samtools
    exits with an error code other than 0.

    Some command line options are associated with parsers.  For
    example, the samtools command "pileup -c" creates a tab-separated
    table on standard output. In order to associate parsers with
    options, an optional list of parsers can be supplied. The list
    will be processed in order checking for the presence of each
    option.

    If no parser is given or no appropriate parser is found, the
    stdout output of samtools/bcftools commands will be returned.

    '''

    dispatch = None
    parsers = None
    collection = None

    def __init__(
        self,
        collection: str,
        dispatch: str,
        parsers: Iterable[Tuple[str, Callable[[Union[str, List[str]]], Union[str, List[str]]]]],
    ):
        self.collection = collection
        self.dispatch = dispatch
        self.parsers = parsers
        self.stderr = []

    def __call__(self, *args: str, **kwargs) -> Union[str, List[str]]:
        '''
        execute a samtools command.

        Keyword arguments:
        catch_stdout -- redirect stdout from the samtools command and
            return as variable (default True)
        save_stdout -- redirect stdout to a filename.
        raw -- ignore any parsers associated with this samtools command.
        split_lines -- return stdout (if catch_stdout is True and stderr
                       as a list of strings.
        '''
        retval, stderr, stdout = _pysam_dispatch(
            self.collection,
            self.dispatch,
            args,
            catch_stdout=kwargs.get("catch_stdout", True),
            save_stdout=kwargs.get("save_stdout", None))

        if kwargs.get("split_lines", False):
            stdout = stdout.splitlines()
            if stderr:
                stderr = stderr.splitlines()

        if retval:
            raise SamtoolsError(
                "%s returned with error %i: "
                "stdout=%s, stderr=%s" %
                (self.collection,
                 retval,
                 stdout,
                 stderr))

        self.stderr = stderr

        # call parser for stdout:
        if not kwargs.get("raw") and stdout and self.parsers:
            for options, parser in self.parsers:
                for option in options:
                    if option not in args:
                        break
                else:
                    return parser(stdout)

        return stdout

    def get_messages(self):
        return self.stderr

    def usage(self):
        '''return the samtools usage information for this command'''
        retval, stderr, stdout = _pysam_dispatch(
            self.collection,
            self.dispatch,
            is_usage=True,
            catch_stdout=True)
        # some tools write usage to stderr, such as mpileup
        if stderr:
            return stderr
        else:
            return stdout


class unquoted_str(str):
    '''Tag a value as an unquoted string. Meta-information in the VCF
    header takes the form of key=value pairs. By default, pysam will
    enclose the value in quotation marks. Tagging that value with
    unquoted_str will prevent this quoting.'''
