#include <getopt.h>
#include <unistd.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bcftools.pysam.h"

FILE * bcftools_stderr = NULL;
FILE * bcftools_stdout = NULL;
const char * bcftools_stdout_fn = NULL;


FILE * bcftools_set_stderr(int fd)
{
  if (bcftools_stderr != NULL)
    fclose(bcftools_stderr);
  bcftools_stderr = fdopen(fd, "w");
  return bcftools_stderr;
}

void bcftools_close_stderr(void)
{
  fclose(bcftools_stderr);
  bcftools_stderr = NULL;
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
  return bcftools_stdout;
}

void bcftools_set_stdout_fn(const char *fn)
{
  bcftools_stdout_fn = fn;
}

void bcftools_close_stdout(void)
{
  fclose(bcftools_stdout);
  bcftools_stdout = NULL;
}

int bcftools_puts(const char *s)
{
  if (fputs(s, bcftools_stdout) == EOF) return EOF;
  return putc('\n', bcftools_stdout);
}


static jmp_buf bcftools_jmpbuf;
static int bcftools_status = 0;

int bcftools_dispatch(int argc, char *argv[])
{
  /* Reset getopt()/getopt_long() processing. */
#if defined __GLIBC__
  optind = 0;
#elif defined _OPTRESET || defined _OPTRESET_DECLARED
  optreset = optind = 1;
#else
  optind = 1;
#endif

  if (setjmp(bcftools_jmpbuf) == 0)
    return bcftools_main(argc, argv);
  else
    return bcftools_status;
}

void bcftools_exit(int status)
{
  bcftools_status = status;
  longjmp(bcftools_jmpbuf, 1);
}
