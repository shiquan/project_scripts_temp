#include "utils.h"
#include <string.h>

int usage()
{
    fprintf(stderr, "bwtFind  <ref> <sequence|stdin>\n");
    return 1;
}

struct args {
    const char *ref_file;
    char *sequence;
    char *buffer;
    char *bwt;
    int *offsets;
    int size;
    int max_size;
    int start;
    int end;
} args = {
    .ref_file = NULL,
    .sequence = NULL,
    .buffer = NULL,
    .bwt = NULL,
    .offsets = NULL,
    .size = 0,
    .max_size = 0,
    .start = -1,
    .end = -1,
};

int parse_args(int argc, char **argv)
{
    if ( argc != 3 )
        return usage();

    args.ref_file = argv[1];
    args.sequence = (char*)strdup(argv[2]);
    int l = strlen(args.sequence);

    if ( l < 4 || l > 100)
        error("Query sequence is must be longer than 10 bases and shorted than 100 bases.");
    
    FILE *fp = fopen(args.ref_file, "r");
    if ( fp == NULL )
        error("%s : %s.", args.ref_file, strerror(errno));
    
    for ( ;; ) {
        int a = fgetc(fp);
        if ( a == EOF )
            break;
        if ( args.size == args.max_size ) {
            if ( args.max_size == 0 ) {
                args.max_size = 1024;
                args.buffer = (char*)malloc(1024*sizeof(char));
            }
            else {
                args.max_size = args.max_size << 1;
                args.buffer = realloc(args.buffer, sizeof(char)*args.max_size);
            }
            if ( args.buffer == NULL ) {
                error_print("Insufficient memeory.");
                //free(args.buffer);
                return 1;
            }
            //args.buffer = (char*)buffer;
        }
        args.buffer[args.size++] = a;
    }
    
    fclose(fp);
    return 0;
}

int find_occ(char *bwt, char c, int pos)
{
    if ( pos < 0 )
        return 0;
    int i;
    int count = 0;
    for ( i = 0; i < pos; ++i )
        if ( bwt[i] == c )
            ++count;
    return count;
}
int bwt_cmp_buff(const void *va, const void *vb)
{
    int a = *((const int *)va);
    int b = *((const int *)vb);
    return strcmp(args.buffer+a, args.buffer+b);
}

int bwt_builder()
{
    int i;
    args.offsets = (int*)malloc((args.size+1)*sizeof(int));
    for ( i = 0; i < args.size+1; i++ )  
        args.offsets[i] = i;
    mergesort(args.offsets, args.size, sizeof(int), bwt_cmp_buff);
    args.bwt = (char*)malloc((args.size+1) *sizeof(char));
    for ( i = 0; i < args.size + 1; i++)
        args.bwt[i] = args.offsets[i] == 0 ? '\0' : args.buffer[args.offsets[i]-1];
    int l = strlen(args.sequence);
    for ( i = 0; i < args.size; i++) {
        if ( args.buffer[args.offsets[i]] == args.sequence[l-1]) {
            args.start = i;
            break;
        }
    }
    for ( ; i < args.size; i++) {
        if ( args.buffer[args.offsets[i]] != args.sequence[l-1]) {
            args.end = i;
            break;
        }
    }

    if (args.end == -1 )
        args.end = args.size -1;
        
    //free(offsets);
    return 0;
}

int bwt_finder_core(char *bwt, int size, char *str)
{
    if ( size == 0 )
        return 0;

    int *occ = (int*)calloc(size+1, sizeof(int));

    int i, counts[256], offsets[256];
    for ( i = 0; i < 256; ++i)
        counts[i] = 0;

    for ( i = 0; i < size; ++i ) {
        char c = bwt[i];
        occ[i] = counts[c];
        counts[c] += 1;
    }
    int total = 0;
    for ( i = 0; i < 256; ++i ) {
        offsets[i] = total;
        total += counts[i];
    }
    int l = strlen(str);
    i = l - 2;
    //int start = 0, end = size -1;
    
    while ( args.start <= args.end && i >= 0) {
        char c = str[i];
        int startOcc = find_occ(bwt, c, args.start);
        int endOcc = find_occ(bwt, c, args.end);
        
        args.start = offsets[c] + startOcc;
        args.end = offsets[c] + endOcc;

        printf("start : %d, end: %d, offsets : %d, startOcc : %d, endOcc : %d\n", args.start, args.end, offsets[c], startOcc, endOcc);                

        if ( startOcc == endOcc)
            break;
        
        i = i - 1;
    }

    //if ( i == 0 )
    //  return 0;
    
    printf("mapped region : start=%u, end=%u\n", args.start, args.end);
    int j = 0, k;
    for ( i = args.start; i <= args.end; i++, j++) {
        printf ("[%d] pos : %d\n", j, args.offsets[i]);        
        int b = args.offsets[i] -2 ;        
        int e = args.offsets[i] + l + 2;
        putc(' ', stdout);
        putc(' ', stdout);
        puts(args.sequence);
        for ( k = b; k < e; ++k)
            putc(args.buffer[k], stdout);
        putc('\n',stdout);
    }
    
    return args.end - args.start + 1;
}

void memory_release()
{
    free(args.buffer);
    free(args.bwt);
    free(args.offsets);
}

void bwt_finder()
{
    bwt_builder();
    bwt_finder_core(args.bwt, args.size, args.sequence);
    memory_release();
}
                
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    bwt_finder();
    
    return 0;
}
