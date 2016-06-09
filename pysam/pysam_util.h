#ifndef PYSAM_UTIL_H
#define PYSAM_UTIL_H

/*! set pysam standard error to point to file descriptor

  Setting the stderr will close the previous stderr.
 */
FILE * pysam_set_stderr(int fd);

/*! set pysam standard output to point to file descriptor

  Setting the stderr will close the previous stdout.
 */
FILE * pysam_set_stdout(int fd);

/*! set pysam standard output to point to filename

 */
void pysam_set_stdout_fn(const char * fn);

/*! set pysam standard error to /dev/null.
  
  Unsetting the stderr will close the previous stderr.
 */
void pysam_unset_stderr(void);

/*! set pysam standard error to /dev/null.
  
  Unsetting the stderr will close the previous stderr.
 */
void pysam_unset_stdout(void);

int pysam_dispatch(int argc, char *argv[]);

void set_optind(int);

#endif
