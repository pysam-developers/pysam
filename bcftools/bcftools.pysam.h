#ifndef bcftools_PYSAM_H
#define bcftools_PYSAM_H

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

extern FILE * bcftools_stderr;

extern FILE * bcftools_stdout;

extern const char * bcftools_stdout_fn;

/*! set pysam standard error to point to file descriptor

  Setting the stderr will close the previous stderr.
 */
FILE * bcftools_set_stderr(int fd);

/*! set pysam standard output to point to file descriptor

  Setting the stdout will close the previous stdout.
 */
FILE * bcftools_set_stdout(int fd);

/*! set pysam standard output to point to filename

 */
void bcftools_set_stdout_fn(const char * fn);

/*! close pysam standard error and set to NULL
  
 */
void bcftools_close_stderr(void);

/*! close pysam standard output and set to NULL
  
 */
void bcftools_close_stdout(void);

int bcftools_puts(const char *s);

int bcftools_dispatch(int argc, char *argv[]);

void PYSAM_NORETURN bcftools_exit(int status);

extern int bcftools_main(int argc, char *argv[]);

/* Define these only in samtools/bcftools C source, not Cython code. */
#if !(defined CYTHON_ABI || defined CYTHON_HEX_VERSION)

/*! Several non-static function names are used in both samtools and bcftools.
    Both libcsamtools.so and libcbcftools.so are loaded simultaneously, leading
    to collisions and wrong functions being called. #define these names so the
    actual symbol names include distinct prefixes to avoid collisions.
 */
#define main_consensus bcftools_main_consensus
#define main_reheader bcftools_main_reheader
#define bam_smpl_init bcftools_bam_smpl_init
#define bam_smpl_destroy bcftools_bam_smpl_destroy
#define read_file_list bcftools_read_file_list

#endif

#endif
