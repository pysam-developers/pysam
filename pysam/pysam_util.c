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

int samtools_main(int argc, char *argv[]);
  
int pysam_dispatch(int argc, char *argv[])
{
  // reset getopt
  optind = 1;

  return samtools_main(argc, argv);
}





