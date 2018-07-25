#ifndef PYSAM_H
#define PYSAM_H

#include "stdio.h"

extern FILE * bcftools_stderr;

extern FILE * bcftools_stdout;

extern const char * bcftools_stdout_fn;

/*! set pysam standard error to point to file descriptor

  Setting the stderr will close the previous stderr.
 */
FILE * bcftools_set_stderr(int fd);

/*! set pysam standard output to point to file descriptor

  Setting the stderr will close the previous stdout.
 */
FILE * bcftools_set_stdout(int fd);

/*! set pysam standard output to point to filename

 */
void bcftools_set_stdout_fn(const char * fn);

/*! set pysam standard error to /dev/null.
  
  Unsetting the stderr will close the previous stderr.
 */
void bcftools_unset_stderr(void);

/*! set pysam standard error to /dev/null.
  
  Unsetting the stderr will close the previous stderr.
 */
void bcftools_unset_stdout(void);

int bcftools_puts(const char *s);

int bcftools_dispatch(int argc, char *argv[]);

void bcftools_set_optind(int);

extern int bcftools_main(int argc, char *argv[]);
  
#endif
