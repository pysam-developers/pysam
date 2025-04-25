#include <string.h>
#include <unistd.h>

// C++-style comment; for-decl; optind

int main(int argc, char **argv) {
    int sum = 0;
    for (int i = optind; i < argc; i++)
        sum += strlen(argv[i]);
    return sum;
}
