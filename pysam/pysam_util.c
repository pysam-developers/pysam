#include <ctype.h>
#include <assert.h>
#include <unistd.h>
#include <stdio.h>
#include "bam.h"
#include "bam_endian.h"
#include "htslib/khash.h"
#include "htslib/ksort.h"
#include "htslib/knetfile.h"
#include "pysam_util.h"


FILE * pysam_stderr = NULL;
FILE * pysam_stdout = NULL;
const char * pysam_stdout_fn = NULL;
int PYSAM_STDOUT_FILENO = STDOUT_FILENO;


FILE * pysam_set_stderr(int fd)
{
  if (pysam_stderr != NULL)
    fclose(pysam_stderr);
  pysam_stderr = fdopen(fd, "w");
  return pysam_stderr;
}

void pysam_unset_stderr(void)
{
  if (pysam_stderr != NULL)
    fclose(pysam_stderr);
  pysam_stderr = fopen("/dev/null", "w");
}

FILE * pysam_set_stdout(int fd)
{
  if (pysam_stdout != NULL)
    fclose(pysam_stdout);
  pysam_stdout = fdopen(fd, "w");
  if (pysam_stdout == NULL)
    {
      fprintf(pysam_stderr, "could not set stdout to fd %i", fd);
    }
  PYSAM_STDOUT_FILENO = fd;
  return pysam_stdout;
}

void pysam_set_stdout_fn(const char *fn)
{
  pysam_stdout_fn = fn;
}

void pysam_unset_stdout(void)
{
  if (pysam_stdout != NULL)
    fclose(pysam_stdout);
  pysam_stdout = fopen("/dev/null", "w");
  PYSAM_STDOUT_FILENO = STDOUT_FILENO;
}

void set_optind(int val)
{
  // setting this in cython via 
  // "from posix.unistd cimport optind"
  // did not work.
  //
  // setting to 0 forces a complete re-initialization
  optind = val;
}



