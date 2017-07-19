#include "utils.h"
#include "fastq.h"
#include "zlib.h"
#include "string.h"
#include "htslib/bgzf.h"
#include "htslib/kseq.h"


KSEQ_INIT(gzFile, gzread)

int check_match(char *seq1, const char *seq2, int mismatch, int length) {
    int i;
    int m = 0;
    for ( i = 0; i < length; ++i ) {
        if ( seq1[i] != seq2[i] ) {
            m++;
            if ( m > mismatch )
                return -1;
        }
    }
    return mismatch;
}
int check_match2(char *seq1, const char *seq2, int mismatch, int length) {
    int i;
    int m = 0;
    for ( i = 0; i < length; ++i ) {
        if ( seq2[i] == 'n' || seq2[i] == 'N')
            continue;
        if ( seq1[i] != seq2[i] ) {            
            m++;
            if ( m > mismatch )
                return -1;
        }
    }
    return mismatch;
}

int load_barcode_file(const char *fn, struct barcode *bc)
{
    BGZF *fp = bgzf_open(fn, "r");
    if ( fp == NULL )
        error("%s : %s", fn, strerror(errno));

    kstring_t string = {0, 0, 0};

    do {
        if ( bgzf_getline(fp, '\n', &string) < 0 )
            break;
        int i, j;
        for ( i = 0; i < string.l; ++i ) {
            if ( string.s[i] == '\t') {
                string.s[i] = '\0';
                break;
            }
        }
        if ( i == string.l )
            error("Failed to parse barcode file: %s", fn);
        
        for ( j = i+1; j < string.l; ++j ) {
            switch (string.s[j]) {
                case 'a':
                case 'A':
                case 'c':
                case 'C':
                case 'g':
                case 'G':
                case 't':
                case 'T':
                case 'n':
                case 'N':
                    continue;
                default:
                    error("Unsupport barcode sequence. %s", string.s);
            }
        }
        if ( bc->n == bc->m ) {
            bc->m = bc->n+8;
            bc->names = (struct name*)realloc(bc->names, bc->m *sizeof(struct name));
        }
        struct name *name = &bc->names[bc->n];
        int k;
        for ( k = 0; k < bc->n; ++k) {
            if ( strncmp(bc->names[k].barcode, string.s+i+1, j-i) == 0 ||
                 strncmp(bc->names[k].name, string.s, i) == 0 ) {
                error_print("Duplicated barcode lines. %s, %s", string.s, string.s+i+1);
                break;
            }
        }
        if ( k == bc->n ) {            
            name->barcode = strndup(string.s+i+1, string.l-i);
            name->name = strndup(string.s, i);
            bc->n++;
        }        
    } while (1);
    
    if ( bc->n == 0 )
        return 1;

    free(string.s);
    bgzf_close(fp);    
    return 0;
}


int clean_barcode_struct(struct barcode *bc)
{
    int i;
    for ( i = 0; i < bc->n; ++i ) {
        struct name *name = &bc->names[i];
        free(name->name);
        free(name->barcode);
        if ( name->fp1 )
            bgzf_close(name->fp1);
        if ( name->fp2 )
            bgzf_close(name->fp2);
    }
    free(bc->names);
    return 0;
}
int check_file_is_fastq(const char *fn)
{
    gzFile fp;
    kseq_t *seq = NULL;
    int l;
    int ret = 0;
    fp = gzopen(fn, "r");
    if ( fp == NULL ) {
        error_print("%s : %s", fn, strerror(errno));
        return -1;
    }
    seq = kseq_init(fp);
    l = kseq_read(seq);
    if ( l >= 0 ) {
        ret =  seq->qual.l > 0 ? 0 : 1;
    } else {   
        ret = -1;
    }

    gzclose(fp);
    kseq_destroy(seq);
    return ret;
}

int check_acgt(const char *seq, int length)
{
    int i;
    for ( i = 0; i < length; ++i ) {
        switch( seq[i] ) {
            case 'a':
            case 'A':
            case 'c':
            case 'C':
            case 'g':
            case 'G':
            case 't':
            case 'T':
                continue;
            default:
                return 1;
        }
    }
    return 0;
}


    
