#ifndef PYSAM_H
#define PYSAM_H

#include "stdio.h"

extern FILE * samtools_stderr;

extern FILE * samtools_stdout;

extern const char * samtools_stdout_fn;

/*! set pysam standard error to point to file descriptor

  Setting the stderr will close the previous stderr.
 */
FILE * samtools_set_stderr(int fd);

/*! set pysam standard output to point to file descriptor

  Setting the stderr will close the previous stdout.
 */
FILE * samtools_set_stdout(int fd);

/*! set pysam standard output to point to filename

 */
void samtools_set_stdout_fn(const char * fn);

/*! set pysam standard error to /dev/null.
  
  Unsetting the stderr will close the previous stderr.
 */
void samtools_unset_stderr(void);

/*! set pysam standard error to /dev/null.
  
  Unsetting the stderr will close the previous stderr.
 */
void samtools_unset_stdout(void);

int samtools_dispatch(int argc, char *argv[]);

void samtools_set_optind(int);

extern int samtools_main(int argc, char *argv[]);
  
#endif
