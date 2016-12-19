#ifndef CNV_BED_HEADER
#define CNV_BED_HEADER
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <htslib/kstring.h>
//   bits 
//   0000 ref-ref
//   0001 ref-del
//   0011 del-del
//   0100 ref-dup
//   1100 dup-dup
//   0101 del-dup
//   1111 malformed

#define check_dup(flag) (flag & (1<<2))
#define check_del(flag) (flag & 1)
#define check_malformed(flag) (((flag&0xf) ^ 0xf) == 0)
struct cnv_bed {
    struct cnv_bed *next;
    int id;
    int start;
    int end;
    int flag;
};

struct cnv_spec {
    const char *fname; // point to input fname
    gzFile fp;
    void *ks;
    kstring_t string;
    // chromosome name to id hash
    void *chrhash;
    // sample name to id hash
    void *samhash;
    // cutoff value to filter CNV regions
    int min_length;
    int max_length;

    int n_samples;
    int m_samples;
    // sample names array
    char **samples;

    int n_chrom;
    int m_chrom;
    // chromosome names array
    char **chrom;    
};


#define KSTRING_INIT {0, 0, 0}

#define CNV_REF       0x0
#define CNV_DEL_HET   0x1
#define CNV_DEL_HOM   0x2
#define CNV_DUP_HET   0x4
#define CNV_DEL_DUP   0x5
#define CNV_DUP_HOM   0x8
#define CNV_MASK      0xf

extern struct cnv_spec * cnv_spec_init();
extern void cnv_spec_destroy(struct cnv_spec *spec);
extern int flag_inconsis(int flag);
extern int combine_flag(int flag1, int flag2);

extern int parse_type(char *type);
const char *explain_type(int flag);       
extern int cnv_load_fname(struct cnv_spec *spec, const char *fname);
extern int cnv_read(struct cnv_spec *spec, struct cnv_bed *line);

#endif
