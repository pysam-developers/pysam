#include <ctype.h>
#include <assert.h>
#include "bam.h"
#include "bam_endian.h"
#include "htslib/khash.h"
#include "htslib/ksort.h"
#include "htslib/knetfile.h"
#include "pysam_util.h"

// Definition of pysamerr
#include "stdio.h"
FILE * pysamerr = NULL;

FILE * pysam_set_stderr(int fd)
{
  if (pysamerr != NULL)
    fclose(pysamerr);
  pysamerr = fdopen(fd, "w");
  return pysamerr;
}

void pysam_unset_stderr(void)
{
  if (pysamerr != NULL)
    fclose(pysamerr);
  pysamerr = fopen("/dev/null", "w");
}


// dummy function - required for samtools integration
void print_error(const char *format, ...)
{
}

// dummy function - required for samtools integration
void print_error_errno(const char *format, ...)
{
}

const char *samtools_version()
{
}


// pysam dispatch function to emulate the samtools
// command line within python.
// taken from the main function in bamtk.c
// added code to reset getopt
int bam_taf2baf(int argc, char *argv[]);
int bam_mpileup(int argc, char *argv[]);
int bam_merge(int argc, char *argv[]);
int bam_index(int argc, char *argv[]);
int bam_sort(int argc, char *argv[]);
int bam_tview_main(int argc, char *argv[]);
int bam_mating(int argc, char *argv[]);
int bam_rmdup(int argc, char *argv[]);
int bam_flagstat(int argc, char *argv[]);
int bam_fillmd(int argc, char *argv[]);
int bam_idxstats(int argc, char *argv[]);
int main_samview(int argc, char *argv[]);
int main_import(int argc, char *argv[]);
int main_reheader(int argc, char *argv[]);
int main_cut_target(int argc, char *argv[]);
int main_phase(int argc, char *argv[]);
int main_cat(int argc, char *argv[]);
int main_depth(int argc, char *argv[]);
int main_bam2fq(int argc, char *argv[]);
int main_pad2unpad(int argc, char *argv[]);
int main_bedcov(int argc, char *argv[]);
int main_bamshuf(int argc, char *argv[]);

int faidx_main(int argc, char *argv[]);

int pysam_dispatch(int argc, char *argv[] )
{
  extern int optind;
#ifdef _WIN32
  setmode(fileno(stdout), O_BINARY);
  setmode(fileno(stdin),  O_BINARY);
#ifdef _USE_KNETFILE
  knet_win32_init();
#endif
#endif

  // reset getopt
  optind = 1;

  if (argc < 2) return 1;
  int retval = 0;
  
  if (strcmp(argv[1], "view") == 0) retval = main_samview(argc-1, argv+1);
  else if (strcmp(argv[1], "import") == 0) retval = main_import(argc-1, argv+1);
  else if (strcmp(argv[1], "mpileup") == 0) retval = bam_mpileup(argc-1, argv+1);
  else if (strcmp(argv[1], "merge") == 0) retval = bam_merge(argc-1, argv+1);
  else if (strcmp(argv[1], "sort") == 0) retval = bam_sort(argc-1, argv+1);
  else if (strcmp(argv[1], "index") == 0) retval = bam_index(argc-1, argv+1);
  else if (strcmp(argv[1], "faidx") == 0) retval = faidx_main(argc-1, argv+1);
  else if (strcmp(argv[1], "idxstats") == 0) retval = bam_idxstats(argc-1, argv+1);
  else if (strcmp(argv[1], "fixmate") == 0) retval = bam_mating(argc-1, argv+1);
  else if (strcmp(argv[1], "rmdup") == 0) retval = bam_rmdup(argc-1, argv+1);
  else if (strcmp(argv[1], "flagstat") == 0) retval = bam_flagstat(argc-1, argv+1);
  else if (strcmp(argv[1], "calmd") == 0) retval = bam_fillmd(argc-1, argv+1);
  else if (strcmp(argv[1], "fillmd") == 0) retval = bam_fillmd(argc-1, argv+1);
  else if (strcmp(argv[1], "reheader") == 0) retval = main_reheader(argc-1, argv+1);
  else if (strcmp(argv[1], "cat") == 0) retval = main_cat(argc-1, argv+1);
  else if (strcmp(argv[1], "targetcut") == 0) retval = main_cut_target(argc-1, argv+1);
  else if (strcmp(argv[1], "phase") == 0) retval = main_phase(argc-1, argv+1);
  else if (strcmp(argv[1], "depth") == 0) retval = main_depth(argc-1, argv+1);
  else if (strcmp(argv[1], "bam2fq") == 0) retval = main_bam2fq(argc-1, argv+1);
  else if (strcmp(argv[1], "pad2unpad") == 0) retval = main_pad2unpad(argc-1, argv+1);
  else if (strcmp(argv[1], "depad") == 0) retval = main_pad2unpad(argc-1, argv+1);
  else if (strcmp(argv[1], "bedcov") == 0) retval = main_bedcov(argc-1, argv+1);
  else if (strcmp(argv[1], "bamshuf") == 0) retval = main_bamshuf(argc-1, argv+1);
  
#if _CURSES_LIB != 0
  else if (strcmp(argv[1], "tview") == 0) retval = bam_tview_main(argc-1, argv+1);
#endif
  else 
    {
      fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
      return 1;
    }
  fflush(stdout);

  return retval;
}





