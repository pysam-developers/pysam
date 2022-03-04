#include <getopt.h>
#include <unistd.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "@pysam@.pysam.h"

FILE * @pysam@_stderr = NULL;
FILE * @pysam@_stdout = NULL;
const char * @pysam@_stdout_fn = NULL;


FILE * @pysam@_set_stderr(int fd)
{
  if (@pysam@_stderr != NULL)
    fclose(@pysam@_stderr);
  @pysam@_stderr = fdopen(fd, "w");
  return @pysam@_stderr;
}

void @pysam@_close_stderr(void)
{
  fclose(@pysam@_stderr);
  @pysam@_stderr = NULL;
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
  return @pysam@_stdout;
}

void @pysam@_set_stdout_fn(const char *fn)
{
  @pysam@_stdout_fn = fn;
}

void @pysam@_close_stdout(void)
{
  fclose(@pysam@_stdout);
  @pysam@_stdout = NULL;
}

int @pysam@_puts(const char *s)
{
  if (fputs(s, @pysam@_stdout) == EOF) return EOF;
  return putc('\n', @pysam@_stdout);
}


static jmp_buf @pysam@_jmpbuf;
static int @pysam@_status = 0;

int @pysam@_dispatch(int argc, char *argv[])
{
  /* Reset getopt()/getopt_long() processing. */
#if defined __GLIBC__
  optind = 0;
#elif defined _OPTRESET || defined _OPTRESET_DECLARED
  optreset = optind = 1;
#else
  optind = 1;
#endif

  if (setjmp(@pysam@_jmpbuf) == 0)
    return @pysam@_main(argc, argv);
  else
    return @pysam@_status;
}

void @pysam@_exit(int status)
{
  @pysam@_status = status;
  longjmp(@pysam@_jmpbuf, 1);
}
