#include <ctype.h>
#include <assert.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "samtools.pysam.h"

FILE * samtools_stderr = NULL;
FILE * samtools_stdout = NULL;
const char * samtools_stdout_fn = NULL;
int samtools_stdout_fileno = STDOUT_FILENO;


FILE * samtools_set_stderr(int fd)
{
  if (samtools_stderr != NULL)
    fclose(samtools_stderr);
  samtools_stderr = fdopen(fd, "w");
  return samtools_stderr;
}

void samtools_unset_stderr(void)
{
  if (samtools_stderr != NULL)
    fclose(samtools_stderr);
  samtools_stderr = fopen("/dev/null", "w");
}

FILE * samtools_set_stdout(int fd)
{
  if (samtools_stdout != NULL)
    fclose(samtools_stdout);
  samtools_stdout = fdopen(fd, "w");
  if (samtools_stdout == NULL)
    {
      fprintf(samtools_stderr, "could not set stdout to fd %i", fd);
    }
  samtools_stdout_fileno = fd;
  return samtools_stdout;
}

void samtools_set_stdout_fn(const char *fn)
{
  samtools_stdout_fn = fn;
}

void samtools_unset_stdout(void)
{
  if (samtools_stdout != NULL)
    fclose(samtools_stdout);
  samtools_stdout = fopen("/dev/null", "w");
  samtools_stdout_fileno = STDOUT_FILENO;
}

int samtools_puts(const char *s)
{
  if (fputs(s, samtools_stdout) == EOF) return EOF;
  return putc('\n', samtools_stdout);
}

void samtools_set_optind(int val)
{
  // setting this in cython via 
  // "from posix.unistd cimport optind"
  // did not work.
  //
  // setting to 0 forces a complete re-initialization
  optind = val;
}



