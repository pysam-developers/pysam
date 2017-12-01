#include <ctype.h>
#include <assert.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "htslib/khash.h"
#include "htslib/ksort.h"
#include "htslib/knetfile.h"

#if !(_POSIX_C_SOURCE >= 200809L || _XOPEN_SOURCE >= 700)
/*
 * A rudimentary emulation of getline() for systems that dont support it
 * natively.  Since this is used for PPD file reading, it assumes (possibly
 * falsely) that BUFSIZ is big enough.
 */
ssize_t
getline(char **line, size_t *linelen, FILE *fp)
{
  if (*linelen == 0) 
    {
      *linelen = BUFSIZ;
      *line = malloc(*linelen);
    }

  memset(*line, 0, *linelen);
  fgets(*line, *linelen, fp);

  return (strlen(*line));

}
#endif



