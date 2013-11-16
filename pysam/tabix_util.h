/* See issue 122
   On some MACOSX systems getline is not defined.
 */
#include <zlib.h>
#include "kseq.h"

#if !(_POSIX_C_SOURCE >= 200809L || _XOPEN_SOURCE >= 700)
#include "unistd.h"
ssize_t getline(char **line, size_t *linelen, FILE *fp);
#endif

KSTREAM_INIT( gzFile, gzread, 16384)

