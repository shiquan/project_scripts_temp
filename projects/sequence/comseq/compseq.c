#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include "sequence.h"
int main(int argc, char **argv)
{
    if ( argc == 1 && !isatty(fileno(stdin))) {
        fprintf(stderr, "Usage : compseq ATCGAG");
        return 1;
    }
    int length = strlen(argv[1]);
    if (length == 1)
        return 1;
    char *rev = rev_seqs(argv[1], length);
    if ( rev ) {
        printf("%d\n", length);
        printf("%s\n", rev);
        return 0;
    }
    return 1;    
}
