#ifndef samtools_PYSAM_H
#define samtools_PYSAM_H

#include <stdio.h>

#ifndef __has_attribute
#define __has_attribute(attribute) 0
#endif
#ifndef PYSAM_NORETURN
#if __has_attribute(__noreturn__) || __GNUC__ >= 3
#define PYSAM_NORETURN __attribute__((__noreturn__))
#else
#define PYSAM_NORETURN
#endif
#endif

extern FILE * samtools_stderr;

extern FILE * samtools_stdout;

extern const char * samtools_stdout_fn;

/*! set pysam standard error to point to file descriptor

  Setting the stderr will close the previous stderr.
 */
FILE * samtools_set_stderr(int fd);

/*! set pysam standard output to point to file descriptor

  Setting the stdout will close the previous stdout.
 */
FILE * samtools_set_stdout(int fd);

/*! set pysam standard output to point to filename

 */
void samtools_set_stdout_fn(const char * fn);

/*! close pysam standard error and set to NULL
  
 */
void samtools_close_stderr(void);

/*! close pysam standard output and set to NULL
  
 */
void samtools_close_stdout(void);

int samtools_puts(const char *s);

int samtools_dispatch(int argc, char *argv[]);

void PYSAM_NORETURN samtools_exit(int status);

extern int samtools_main(int argc, char *argv[]);

/* Define these only in samtools/bcftools C source, not Cython code. */
#if !(defined CYTHON_ABI || defined CYTHON_HEX_VERSION)

/*! Several non-static function names are used in both samtools and bcftools.
    Both libcsamtools.so and libcbcftools.so are loaded simultaneously, leading
    to collisions and wrong functions being called. #define these names so the
    actual symbol names include distinct prefixes to avoid collisions.
 */
#define main_consensus samtools_main_consensus
#define main_reheader samtools_main_reheader
#define bam_smpl_init samtools_bam_smpl_init
#define bam_smpl_destroy samtools_bam_smpl_destroy
#define read_file_list samtools_read_file_list

#endif

#endif
