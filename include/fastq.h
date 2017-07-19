#ifndef BARCODE_HEADER
#define BARCODE_HEADER

#include "utils.h"
#include "htslib/bgzf.h"

struct name {
    char *barcode;
    char *name;
    BGZF *fp1;
    BGZF *fp2;
};

struct barcode {
    int n, m;
    struct name *names;
};



extern int check_match(char *s1, const char *s2, int m, int l);
extern int check_match2(char *s1, const char *s2, int m, int l);

extern int load_barcode_file(const char*fn, struct barcode*b);

extern int clean_barcode_struct(struct barcode *bc);

// return -1 for unknown or error, 0 for fastq, 1 for fasta.
extern int check_file_is_fastq(const char *fn);

extern int check_acgt(const char *seq, int l);
#endif
