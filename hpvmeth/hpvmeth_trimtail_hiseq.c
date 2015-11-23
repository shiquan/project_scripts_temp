/* This program is used to trim adaptor sequnence and randam barcode sequences for hpv methylation projects
 * shiquan.cn@gmail.com   2015/11/23
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <klib/kseq.h>
#include <klib/kstring.h>
#include <klib/khash.h>
#include <errno.h>
#include <assert.h>
#include <zlib.h>

struct pair {
    gzFile fq1;
    gzFile fq2;
};
struct {
    char *fastq1;
    char *fastq2;
    char *samplcode;
    char *outdir;
    int n_samples;
    char **sam;
    struct pair fin;
    struct pair *fout;
    struct pair error;
    struct {
        uint64_t all;
        uint64_t bias;
        uint64_t *samples;
    } s;

} args = {
    NULL,
    NULL,
    NULL,
    NULL,
    0,
    NULL,
    0,
    NULL,
    0,
    { 0, 0, NULL }
};

KSEQ_INIT(gzFile, gzread)

struct adaptor {
    char * name;
    char *barcode;
};

struct adaptor adaptors[] = {
    {"AZ1", "ATATTATTA"},
    {"AZ2", "TTATTATTA"},
    {"AZ3", "TTAATATTA"},
    {"AZ4", "TTAATATAA"},
    {"AZ5", "TTAATATAT"},
};

void release_memory()
{
    free(args.fastq1);
    free(args.fastq2);
    if (args.outdir) free(args.outdir);
    int i;
    for (i=0; i < args.n_samples; ++i) free(args.sam[i]);
    if (args.n_samples) {
        free(args.sam);
        free(args.s.samples);
    }
}
void prase_argv(char ***argv, int *argc)
{
    const char *name = (*argv)[0];
    (*argv)++; (*argc)--;
    if (*argc == 0) {
	    fprintf(stderr, "Usage : %s -f1 fastq1.fq.gz -f2 fastq2.fq.gz -s barcode.txt -o out_dir\n", name);
	    exit(1);
    }
    const char *cmd = (*argv)[0];
    while (*argc > 0) {
        if (!strcmp(cmd, "-h") || !strcmp(cmd, "--help")) {
            fprintf(stderr, "Usage : %s -f1 fastq1.fq.gz -f2 fastq2.fq.gz -s barcode.txt -o out_dir\n", name);
            exit(1);
        }
        if (!strcmp(cmd, "-f1") || !strcmp(cmd, "--fastq1")) {
            (*argv)++;
            (*argc)--;
            args.fastq1 = strdup((*argv)[0]);
        } else if (!strcmp(cmd, "-f2") || !strcmp(cmd, "--fastq2")) {
            (*argv)++;
            (*argc)--;
            args.fastq2 = strdup((*argv)[0]);
        /* } else if (!strcmp(cmd, "-s") || !strcmp(cmd, "--barcode")) { */
        /*     (*argv)++; */
        /*     (*argc)--; */
        /*     args.samplcode = strdup((*argv)[0]); */
        } else if (!strcmp(cmd, "-o") || !strcmp(cmd, "--out-dir")) {
            (*argv)++;
            (*argc)--;
            args.outdir = strdup((*argv)[0]);
        } else {
            fprintf(stderr, "error parameter: %s\n", cmd);
            exit(1);
        }
        (*argv)++;
        (*argc)--;
    }
    if ( args.fastq1 == NULL || args.fastq2 == NULL ) {
        fprintf(stderr, "Set --fastq1 and --fastq2 parameters please. See more information with -h.\n");
        exit(1);
    }
    /* if ( args.samplcode == NULL ) { */
    /*     fprintf(stderr, "Specify barcode and sample list by --barcode. [sample\tbarcode]\n"); */
    /*     exit(1); */
    /* } */
}
char comp(char c)
{
    if (c == 'A') return 'T';
    if (c == 'T') return 'A';
    if (c == 'C') return 'G';
    if (c == 'G') return 'C';
    return c;
}
int check_sample(char *seq1, char *seq2)
{
    char seq[9];
    if (!strncmp(seq1+23, "ATT", 3) && !strncmp(seq1+35,"TGA", 3)) {
        memcpy(seq, seq1+26, 9);
    }
    else if (!strncmp(seq2+23, "ATT", 3) && !strncmp(seq2+35,"TGA", 3)) {
        memcpy(seq, seq2+26, 9);
    }
    else if (!strncmp(seq1+23, "TCC", 3) && !strncmp(seq2+35,"ACT", 3)) {
        memcpy(seq, seq1+26, 9);
        int i;
        for (i=0; i< 9; ++i) seq[i] = comp(seq[i]);
    }
    else if (!strncmp(seq2+23, "TCC", 3) && !strncmp(seq2+35,"ACT", 3)) {
        memcpy(seq, seq2+26, 9);
    }
    else
    {
        return -1;
    }

    int i;
    for (i = 0; i < args.n_samples; ++i)
        if (!strcmp(seq, adaptors[i].barcode)) return i;

    return i;
}
void load_fastq()
{
    args.n_samples = sizeof(adaptors)/ sizeof(adaptors[0]);
    assert(args.n_samples > 0);
    args.s.samples = calloc(args.n_samples, sizeof(uint64_t));
    args.fout = calloc(args.n_samples, sizeof(struct pair));

#define BRANCH(__file, __fp, __mode) do {                               \
        __fp = gzopen(__file, __mode);                                  \
        if ( __fp == NULL ) {                                           \
            fprintf(stderr, "[load_fastq] %s : %s\n", __file, strerror(errno)); \
            exit(1);                                                    \
        }                                                               \
    } while(0)

    BRANCH(args.fastq1, args.fin.fq1, "r");
    BRANCH(args.fastq2, args.fin.fq2, "r");

    kstring_t str = { 0, 0, 0 };
    int i;

    for (i=0; i<args.n_samples; ++i) {
        str.m = 0;
        if (args.outdir != NULL) {
            kputs(args.outdir, &str);
            if (str.s[str.l-1] != '/') kputc('/', &str);
        }
        kputs(adaptors[i].name, &str);
        kputs("_1.fq.gz", &str);
        BRANCH(str.s,args.fout[i].fq1,"w");
        str.l -= 8;
        BRANCH(str.s,args.fout[i].fq2,"w");
    }
    str.m = 0;
    if (args.outdir == NULL) {
        kputs(args.outdir, &str);
        if (str.s[str.l-1] != '/') kputc('/', &str);
    }
    kputs("error_1.fq.gz", &str);
    BRANCH(str.s, args.error.fq1, "w");
    str.l -= 13;
    BRANCH(str.s, args.error.fq2, "w");
    if (str.m) free(str.s);
#undef BRANCH

    kseq_t *seq1;
    kseq_t *seq2;

    seq1 = kseq_init(args.fin.fq1);
    seq2 = kseq_init(args.fin.fq2);

    while ( kseq_read(seq1)>0 && kseq_read(seq2)>0 ) {
        int i;
        str.l = 0;
        i = check_sample(seq1->seq.s, seq2->seq.s);
        if (i == -1) {
            kputc('@',&str);
            kputs(seq1->name.s, &str);      kputc('\n', &str);
            kputs(seq1->seq.s, &str);       kputc('\n', &str);
            kputs(seq1->qual.s, &str);      kputc('\n', &str);
            gzwrite(args.error.fq1, str.s, str.l);
            str.l = 0;
            kputc('@',&str);
            kputs(seq2->name.s, &str);      kputc('\n', &str);
            kputs(seq2->seq.s, &str);       kputc('\n', &str);
            kputs(seq2->qual.s, &str);      kputc('\n', &str);
            gzwrite(args.error.fq2, str.s, str.l);
            str.l = 0;
            args.s.bias ++;
        } else {
            seq1->seq.s += 37;
            seq1->qual.s += 37;
            seq2->seq.s += 37;
            seq2->qual.s += 37;
            kputc('@',&str);
            kputs(seq1->name.s, &str);      kputc('\n', &str);
            kputs(seq1->seq.s, &str);       kputc('\n', &str);
            kputs(seq1->qual.s, &str);      kputc('\n', &str);
            gzwrite(args.error.fq1, str.s, str.l);
            str.l = 0;
            kputc('@',&str);
            kputs(seq2->name.s, &str);      kputc('\n', &str);
            kputs(seq2->seq.s, &str);       kputc('\n', &str);
            kputs(seq2->qual.s, &str);      kputc('\n', &str);
            gzwrite(args.error.fq2, str.s, str.l);
            args.s.samples[i]++;
        }

    }
    gzclose(args.fin.fq1);
    gzclose(args.fin.fq2);
    kseq_destroy(seq1);
    kseq_destroy(seq2);
    for (i=0; i<args.n_samples; ++i) { gzclose(args.fout[i].fq1); gzclose(args.fout[i].fq2); }    
}
int main(int argc, char **argv)
{
    prase_argv(&argv, &argc);
    load_fastq();
    release_memory();
    return 0;
}
