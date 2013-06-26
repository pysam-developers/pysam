/* See issue 122
   On some MACOSX systems getline is not defined.
 */
#if !(_POSIX_C_SOURCE >= 200809L || _XOPEN_SOURCE >= 700)
ssize_t getline(char **line, size_t *linelen, FILE *fp);
#endif
