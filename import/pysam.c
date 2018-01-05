#include <ctype.h>
#include <assert.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "@pysam@.pysam.h"

FILE * @pysam@_stderr = NULL;
FILE * @pysam@_stdout = NULL;
const char * @pysam@_stdout_fn = NULL;
int @pysam@_stdout_fileno = STDOUT_FILENO;


FILE * @pysam@_set_stderr(int fd)
{
  if (@pysam@_stderr != NULL)
    fclose(@pysam@_stderr);
  @pysam@_stderr = fdopen(fd, "w");
  return @pysam@_stderr;
}

void @pysam@_unset_stderr(void)
{
  if (@pysam@_stderr != NULL)
    fclose(@pysam@_stderr);
  @pysam@_stderr = fopen("/dev/null", "w");
}

FILE * @pysam@_set_stdout(int fd)
{
  if (@pysam@_stdout != NULL)
    fclose(@pysam@_stdout);
  @pysam@_stdout = fdopen(fd, "w");
  if (@pysam@_stdout == NULL)
    {
      fprintf(@pysam@_stderr, "could not set stdout to fd %i", fd);
    }
  @pysam@_stdout_fileno = fd;
  return @pysam@_stdout;
}

void @pysam@_set_stdout_fn(const char *fn)
{
  @pysam@_stdout_fn = fn;
}

void @pysam@_unset_stdout(void)
{
  if (@pysam@_stdout != NULL)
    fclose(@pysam@_stdout);
  @pysam@_stdout = fopen("/dev/null", "w");
  @pysam@_stdout_fileno = STDOUT_FILENO;
}

void @pysam@_set_optind(int val)
{
  // setting this in cython via 
  // "from posix.unistd cimport optind"
  // did not work.
  //
  // setting to 0 forces a complete re-initialization
  optind = val;
}



