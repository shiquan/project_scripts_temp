/* Author:    Shi Quan            */
/* Email:     shiquan@genomcis.cn */

#include "utils.h"
#include "number.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include <zlib.h>

static char * program_name =  "dyncutadaptor";
static char * Version = "v0.1.4";

unsigned char bam_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

char *bam_nt16_rev_table = "=ACMGRSVTWYHKDBN";

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };

// Illumina
//static uint8_t illumina[19] = { 2, 8, 4, 8, 2, 8, 2, 8, 8, 1, 8, 1, 2, 1, 2, 1, 8, 2, 8};
// TTCAGCCT
//static uint8_t bgiseq_ad153[] = {};

struct args {
    int seed;
    int slave_mode;
    int minimum;
    int mismatch;
    int tail;
    uint8_t *adaptor;
    const char *fastq1;
    const char *fastq2;
    const char *out1;
    const char *out2;
} args = {
    .seed = 5,
    .slave_mode = 1,
    .minimum = 0,
    .mismatch = 1,
    .tail = 3,
    .adaptor = 0,
    .fastq1 = 0,
    .fastq2 = 0,
    .out1 = 0,
    .out2 = 0,
};

void seq_comp(uint8_t * seq, int length)
{
    int i;
    for ( i = 0; i < length>>1; ++i ) {
        int8_t t = seq_comp_table[seq[length - 1 - i]];
        seq[length - 1 - i] = seq_comp_table[seq[i]];
        seq[i] = t;
    }
    if ( length & 1 )
        seq[i] = seq_comp_table[seq[i]];
}

KSEQ_INIT(gzFile, gzread)

uint8_t * seq2code(char * str, int n) 
{
    int i;
    uint8_t *a;
    a = calloc(n, sizeof(uint8_t));
    for(i = 0; i < n; ++i) {
        a[i] = bam_nt16_table[str[i]];
    }
    return a;
}

static int * BMprep(const uint8_t * pat, int m)
{
    int i, *suff, *prep, *bmGs, *bmBc;
    prep = (int*) calloc(m + 15, sizeof(int));
    bmGs = prep; bmBc = prep + m;
    for (i = 0; i < 15; ++i)
        bmBc[i] = m;
    for (i = 0; i < m - 1; ++i)
        bmBc[pat[i]] = m - i - 1;
    suff = (int*) calloc(m, sizeof(int));
    int f = 0, g;
    suff[m - 1] = m;
    g = m - 1;
    for(i = m - 2; i >= 0; --i) {
        if (i > g && suff[i + m - 1 - f] < i - g) {
            suff[i] = suff[i + m - 1 - f];
        } else {
            if (i < g) g = i;
            f = i;
            while (g >= 0 && pat[g] == pat[g + m - 1 - f])
                --g;
            suff[i] = f - g;
        }
    }
    int j = 0;
    for (i = 0; i < m; ++i)
        bmGs[i] = m;
    for (i = m -1; i >= 0; --i) {
        if (suff[i] == i + 1) {
            for (; j < m - 1 - i; ++j) {
                bmGs[j] = m - 1 - i;
            }
        }
    }
    for (i = 0; i <= m - 2; ++i) {
        bmGs[m - 1 - suff[i]] = m - 1 -i;
    }
    free(suff);
    return prep;
}

int location(const uint8_t * str, int n, const uint8_t * pat, int m, int *prep)
{
    int i, j, *bmGs, *bmBc;
    //prep = BMprep(pat, m);
    bmGs = prep;
    bmBc = prep + m;
    j  = 0;
    while ( j <= n - m) {
        int b = 0, k;
        for ( i = m - 1; i >= 0; --i) {
            if (!(pat[i] & str[i+j])) {
                if (b) {
                    i = k;
                    break;
                } else {
                    k = i;
                    b++;
                }
            }
        }
        if (i >= 0) {
            int max = bmBc[str[i+j]] - m + 1 + i;
            if (max < bmGs[i])
                max = bmGs[i];
            j += max;
        } else {
            //free(prep);
            return j;
        }
    }
    j = n - m;
    while (j < n - args.tail) {
        for (i = 0 ; i < n - j && pat[i] == str[j + i]; ++i);
        if ( i < n - j) {
            ++j;
        } else {
            //free(prep);
            return j;
        }
    }
    //free(prep);
    return 0;
}

int loc_adaptor(uint8_t * str, int n, uint8_t * pat, int m, int * prep, int *k)
{
    *k = location(str, n, pat, m, prep);
    if ( *k ) return 1;
    seq_comp(str, n);
    *k = location(str, n, pat, m, prep);
    if ( *k ) return -1;
    return 0;		
}

int check_loc(uint8_t * seq, int l, uint8_t * a, int m, int loc) 
{
    int i, j, k = 0;
    j = l - loc > m ? m + loc: l;
    uint8_t * t = a;
    for (i = loc; i < j; ++i) {
        if(seq[i] != *t)
            k ++;
        t++;
    }
    if ( k > j -2 || k > args.mismatch)
        return 0;
    return 1;
}

void seq_reverse(char * seq, int len) 
{
    int i;
    for ( i = 0; i < len>>1; ++i) {
        unsigned char t = seq[i];
        seq[i] = seq[len - i -1];
        seq[len - i - 1] = t;
    }
}

int seq_reloc(char * seq, unsigned long * len, unsigned long loc) 
{
    seq_reverse(seq, *len);
    seq[loc] = 0;
    seq_reverse(seq, loc);
    *len = loc;
    return 1;
}

long long filter_reads = 0;

int cut_adaptor(kseq_t * seq1, kseq_t *seq2, uint8_t * pat, int len, int * prep) 
{
    int m = 0, n = 0;
    int l, loc = 0;		
    uint8_t * s, * p;
    s = seq2code(seq1->seq.s, seq1->seq.l);
    p = seq2code(seq2->seq.s, seq2->seq.l);
    int length = seq1->seq.l;
    if ( args.slave_mode ) {
        m = location(s, seq1->seq.l, pat, len, prep);
        if ( check_loc(s, seq1->seq.l, pat, len, m) ) {
            loc = m;
        } else {
            n = location(p, seq2->seq.l, pat, len, prep);
            if ( check_loc(p, seq2->seq.l, pat, len, n) ) {
                loc = n;
            }
        }
        if (loc) {
            if ( args.minimum && loc < args.minimum) {
                goto MINI;
            }
            seq1->seq.s[loc] = 0; seq1->seq.l = l;
            seq1->qual.s[loc] = 0; seq1->qual.l = l;
            seq2->seq.s[loc] = 0; seq2->seq.l = l;
            seq2->qual.s[loc] = 0; seq2->qual.l = l;
        } else {
            seq_comp(s, seq1->seq.l);
            seq_comp(p, seq2->seq.l);
            
            m = location(s, seq1->seq.l, pat, len, prep);
            if (check_loc(s, seq1->seq.l, pat, len, m)) {
                loc = m;
                
            } else {
                n = location(p, seq2->seq.l, pat, len, prep);
                if (check_loc(p, seq2->seq.l, pat, len, n))
                    loc = n;
            }
            if ( loc ) {
                if (args.minimum && loc < args.minimum) {
                    goto MINI;
                }
                seq_reloc(seq1->seq.s, &seq1->seq.l, loc);
                seq_reloc(seq2->seq.s, &seq2->seq.l, loc);
                seq_reloc(seq1->qual.s, &seq1->qual.l, loc);
                seq2->qual.s[loc] = 0; seq2->qual.l = l;
            } else {
                free(s); free(p);
                return 0;
            }
        }
    } else {
        m = location(s, seq1->seq.l, pat, len, prep);
        n = location(p, seq2->seq.l, pat, len, prep);
        if (m || n) {
            loc = m > n ? n > 0 ? n : m : m > 0 ? m : n;
            if (args.minimum && loc < args.minimum) {
                goto MINI;
            }
            seq1->seq.s[loc] = 0; seq1->seq.l = l;
            seq1->qual.s[loc] = 0; seq1->qual.l = l;
            seq2->seq.s[loc] = 0; seq2->seq.l = l;
            seq2->qual.s[loc] = 0; seq2->qual.l = l;
        } else { 
            seq_comp(s, seq1->seq.l);
            seq_comp(p, seq2->seq.l);
            m = location(s, seq1->seq.l, pat, len, prep);
            n = location(p, seq2->seq.l, pat, len, prep);
            if (m || n) {
                loc = m > n ? n > 0 ? n : m : m > 0 ? m : n;
                if ( args.minimum && loc < args.minimum) {
                    goto MINI;
                }
                seq_reloc(seq1->seq.s, &seq1->seq.l, loc);
                seq_reloc(seq2->seq.s, &seq2->seq.l, loc);
                seq_reloc(seq1->qual.s, &seq1->qual.l, loc);
                seq2->qual.s[loc] = 0; seq2->qual.l = l;
            } else {
                free(s); free(p);
                return 0;
            }
        }
    }
    free(s); free(p);
    return 1;
MINI:
    filter_reads++;
    free(s); free(p);
    return -1;
}

int check_name(kstring_t name1, kstring_t name2)
{
    char * s1 = name1.s;
    char * s2 = name2.s;
    size_t n;
    for(n = 0; n < name1.l - 1; ++n, ++s1, ++s2) {
        if (*s1 != *s2) {
            return 0;
        }
    }
    return 1;
}

int loadfastq_pe(const char * pe1, const char * pe2, uint8_t * a)
{
    gzFile fq1, fq2;
    kseq_t * seq1, * seq2;
    int l, m;
    long long all_reads = 0;
    long long cutted_reads = 0;
    fq1 = gzopen(pe1, "r");
    fq2 = gzopen(pe2, "r");
    if (fq1 == NULL) {
        fprintf(stderr, "[loadfastq_pe] %s : %s\n", pe1, strerror(errno));
        return 0;
    }
    if (fq2 == NULL) {
        fprintf(stderr, "[loadfastq_pe] %s : %s\n", pe2, strerror(errno));
        return 0;
    }
    seq1 = kseq_init(fq1);
    seq2 = kseq_init(fq2);
    FILE * fp1, * fp2;
    fp1 = fopen(args.out1, "w");
    fp2 = fopen(args.out2, "w");
    int *prep;
    prep = BMprep(a, args.seed);
    while (kseq_read(seq1) > 0 && kseq_read(seq2) > 0) {
        all_reads++;
        if (check_name(seq1->name, seq2->name) < 1) {
            fprintf(stderr, "[loadfastq_pe] The fastq files not paired:\n %s\n%s\n",
                    seq1->name.s, seq2->name.s);
            goto LOADFQ_ERROR;
        }
        int n;
        n = cut_adaptor(seq1, seq2, a, args.seed, prep);
        if ( n )
            cutted_reads++;
        if ( n < 0 )
            continue;
        fputc('@', fp1);
        fputs(seq1->name.s, fp1);
        fputc('\n', fp1); 
        fputs(seq1->seq.s, fp1);
        fputs("\n", fp1);
        fputs("+\n",fp1);
        fputs(seq1->qual.s, fp1);
        fputs("\n", fp1);
        fputc('@', fp2);
        fputs(seq2->name.s, fp2);
        fputs("\n", fp2); 
        fputs(seq2->seq.s, fp2);
        fputs("\n", fp2);
        fputs("+\n",fp2);
        fputs(seq2->qual.s, fp2);
        fputs("\n", fp2);
    }
    free(prep);
    kseq_destroy(seq1);
    kseq_destroy(seq2);
    gzclose(fq1);
    gzclose(fq2);
    fclose(fp1);
    fclose(fp2);
    fprintf(stdout, "dealed with %lld reads.\nTotal reads is %lld\nfilter %lld reads\n", cutted_reads, all_reads, filter_reads);
    return 1;
    
  LOADFQ_ERROR:
    kseq_destroy(seq1);
    kseq_destroy(seq2);
    gzclose(fq1);
    gzclose(fq2);
    return 0;
}


int usage(int s)
{
    if ( s ) {
        fprintf(stderr,"\
Usage: %s --fastq1 [in1.fq.gz] --fastq2 [in2.fq.gz] \n\
======================================================\n\
-fastq1, -f         Fastq file of read1.\n\
-fastq2, -r         Fastq file of read2.\n\
-outfq1, -o         New fastq file of read1.\n\
-outfq2, -p         New fastq file of read2.\n\
-seed,   -s         Initial length of adaptor.[5]\n\
-slave,  -d         Cut all sequence like adaptor.\n\
-mis,    -i         Tolerate mismatchs.[1]\n\
-tail,   -t         Don't cut tail in the last n bp.[3]\n\
-min,    -m         Don't keep seqences shorter.[0]\n\
-help,   -h         See this information.\n\
-adaptor            Adaptor sequence.\n\
======================================================\n\
Author: Shi Quan (shiquan@genomics.cn)\n\
Pages:\n\
https://github.com/shiquan/small_projects_collections\n\
Version: %s\n\
", program_name, Version);
    } else {
        fprintf(stderr,"\
Usage: %s --fastq1 [in1.fq.gz] --fastq2 [in2.fq.gz]\n\
use -help for more information.\n\
",program_name);
    }
    return 1;
}

int parse_args(int ac, char **av)
{
    const char *minimum = 0;
    const char *seed = 0;
    const char *tail = 0;
    const char *mismatch = 0;
    const char *adaptor = 0;
    if ( ac == 0 )
        return usage(0);
    
    int i;
    for ( i = 0; i < ac; ) {
        const char *a = av[i++];
        if ( strcmp(a, "-h") == 0 || strcmp(a, "-help") == 0 )
            return usage(1);

        const char **arg_var = 0;
        if ( (strcmp(a, "-fastq1") == 0 || strcmp(a, "-f") == 0)  && args.fastq1 == 0 )
            arg_var = &args.fastq1;
        else if ( (strcmp(a, "-fastq2") == 0 || strcmp(a, "-r") == 0 ) && args.fastq2 == 0 )
            arg_var = &args.fastq2;
        else if ( (strcmp(a, "-outfq1") == 0 || strcmp(a, "-o") == 0 ) && args.out1 == 0 )
            arg_var = &args.out1;
        else if ( (strcmp(a, "-outfq2") == 0 || strcmp(a, "-p") == 0 ) && args.out2 == 0 )
            arg_var = &args.out2;
        else if ( (strcmp(a, "-min") == 0 || strcmp(a, "-m") == 0 ) && minimum == 0 )
            arg_var = &minimum;
        else if ( (strcmp(a, "-seed") == 0 || strcmp(a, "-s") == 0 ) && seed == 0 )
            arg_var = &seed;
        else if ( (strcmp(a, "-mis") == 0 || strcmp(a, "-i") == 0 ) && mismatch == 0 )
            arg_var = &mismatch;
        else if ( strcmp(a, "-adaptor") == 0 && adaptor == 0 )
            arg_var = &adaptor;
        
        if ( arg_var != 0 ) {
            if ( i == ac ) {
                error("Missing argument after %s.", a);
            }
            *arg_var = av[i++];
            continue;
        }

        if ( strcmp(a, "-slave") == 0 ) {
            args.slave_mode = 1;
            continue;
        }
        error("Unknown argument %s.", a);
    }

    if ( args.fastq1 == NULL || args.fastq2 == NULL ) {
        error("No fastq1 and/or fastq2 specified.");
    }

    if ( seed ) {
        args.seed = str2int((char*)seed);
    }

    if ( minimum ) {
        args.minimum = str2int((char*)minimum);
    }

    if ( mismatch ) {
        args.mismatch = str2int((char*)mismatch);
    }

    if ( tail ) {
        args.tail = str2int((char*)tail);
    }

    if ( adaptor ) {
        int length = strlen(adaptor);
        args.adaptor = seq2code((char*)adaptor, length);
        if (args.seed > length ) {
            args.seed = length;
            warnings("Force set seed to %d.", length);
        }
        if ( args.adaptor == NULL )
            error("Failed to recongnise adaptor sequence. %s.", adaptor);
    } else {
        // Default is BGISEQ AD153 adaptor       
        args.adaptor = seq2code("TTCAGCCT",8);
        if (args.seed > 8) {
            args.seed = 8;
            warnings("Force seed to 8.");
        }
    }

    if ( args.seed > 19 || args.seed < 5) {
        warnings("Seed must set in [5, 19]. Force to 5..");
        args.seed =5;
    }
    if ( args.mismatch > 5 ) {
        warnings("No more than 5 mismatches allowed. Force to 1.");
        args.mismatch = 1;
    }

    if ( args.tail > 5 ) {
        warnings("Do NOT support 5nt longer tails. Force to 5.");
        args.tail = 5;
    }

    if ( args.out1 == NULL )
        args.out1 = "reads1.fq";

    if ( args.out2 == NULL )
        args.out2 = "reads2.fq";
    
    return 0;
}

int main(int argc, char *argv[])
{
    if ( parse_args(--argc, ++argv) )
        return 1;
    loadfastq_pe(args.fastq1, args.fastq2, args.adaptor);
    return 1;
}
