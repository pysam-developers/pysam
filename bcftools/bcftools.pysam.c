#include <ctype.h>
#include <assert.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bcftools.pysam.h"

FILE * bcftools_stderr = NULL;
FILE * bcftools_stdout = NULL;
const char * bcftools_stdout_fn = NULL;
int bcftools_stdout_fileno = STDOUT_FILENO;


FILE * bcftools_set_stderr(int fd)
{
  if (bcftools_stderr != NULL)
    fclose(bcftools_stderr);
  bcftools_stderr = fdopen(fd, "w");
  return bcftools_stderr;
}

void bcftools_unset_stderr(void)
{
  if (bcftools_stderr != NULL)
    fclose(bcftools_stderr);
  bcftools_stderr = fopen("/dev/null", "w");
}

FILE * bcftools_set_stdout(int fd)
{
  if (bcftools_stdout != NULL)
    fclose(bcftools_stdout);
  bcftools_stdout = fdopen(fd, "w");
  if (bcftools_stdout == NULL)
    {
      fprintf(bcftools_stderr, "could not set stdout to fd %i", fd);
    }
  bcftools_stdout_fileno = fd;
  return bcftools_stdout;
}

void bcftools_set_stdout_fn(const char *fn)
{
  bcftools_stdout_fn = fn;
}

void bcftools_unset_stdout(void)
{
  if (bcftools_stdout != NULL)
    fclose(bcftools_stdout);
  bcftools_stdout = fopen("/dev/null", "w");
  bcftools_stdout_fileno = STDOUT_FILENO;
}

int bcftools_puts(const char *s)
{
  if (fputs(s, bcftools_stdout) == EOF) return EOF;
  return putc('\n', bcftools_stdout);
}

void bcftools_set_optind(int val)
{
  // setting this in cython via 
  // "from posix.unistd cimport optind"
  // did not work.
  //
  // setting to 0 forces a complete re-initialization
  optind = val;
}



