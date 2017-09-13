// Check the genepred data and transcripts sequences is consistant.
#include "utils.h"
#include "genepred.h"
//#include "faidx_def.h"
#include "htslib/faidx.h"
#include "htslib/tbx.h"
#include "htslib/kstring.h"
#include "ksw.h"
#include "number.h"
#include "sequence.h"
#include "htslib/kseq.h"
#include <pthread.h>
#include <zlib.h>

KSTREAM_INIT(gzFile, gzread, 16384)

int usage()
{
    fprintf(stderr,
            "realign_trans - realgin transcripts to reference genome and generate gap information.\n"
            "options:\n"
            " -data     Gene prediction format database.\n"
            " -rna      Transcripts sequence database included in FASTA format, indexed by samtools faidx.\n"
            " -ref      Genome reference database in FASTA format, indexed by samtools faidx.\n"
            " -format   Format of gene prediction database, refgene is default. [genepred,refflat,refgene]\n"
            //" -p        Set threads [1].\n"
	    " -gapo     Penalty score for gap open.\n"
	    " -gape     Penalty score for gap extension.\n"
	    " -sa       Penalty score for match.\n"
	    " -sb       Penalty score for mismatch.\n"	   
            "\nwebsite: https://github.com/shiquan/small_projects\n");
    return 1;
}
struct buffer {
    int n, m;
    kstring_t *buffer;
};
struct args {
    const char *genepred_fname;
    const char *refrna_fname;
    const char *reference_fname;
    faidx_t *rna_fai;
    faidx_t *ref_fai;
    const char *format;
    int gapo;
    int gape;
    int sa;
    int sb;
    int8_t mat[25];
    int threads;
    int chunk_size;
    kstream_t *ks;
    struct buffer *pbuf;
} args = {
    .genepred_fname = NULL,
    .refrna_fname = NULL,
    .reference_fname = NULL,
    .rna_fai = NULL,
    .ref_fai = NULL,
    .format = "refgene",
    .gapo = 5,
    .gape = 2,
    .sa = 4,
    .sb = 4,
    .threads = 1,
    .chunk_size = 10000000,
    //.fp = NULL,
    .ks = NULL,
    .pbuf = NULL,
};

int parse_args(int ac, char **av)
{
    if ( ac == 1 )
        return usage();

    int i;
    const char *gapo = NULL;
    const char *gape = NULL;
    const char *sa = NULL;
    const char *sb = NULL;
    const char *threads = NULL;
    for ( i = 1; i < ac; ) {
        const char *a = av[i++];
        const char **var = 0;
        if ( strcmp(a, "-h") == 0 )
            return usage();

        if ( strcmp(a, "-data") == 0 && args.genepred_fname == NULL )
            var = &args.genepred_fname;
        else if ( strcmp(a, "-rna") == 0 && args.refrna_fname == NULL )
            var = &args.refrna_fname;
        else if ( strcmp(a, "-ref") == 0 && args.reference_fname == NULL )
            var = &args.reference_fname;
        else if ( strcmp(a, "-format") == 0 )
            var = &args.format;
        else if ( strcmp(a, "-gape") == 0 )
            var = &gape;
        else if ( strcmp(a, "-gapo") == 0 )
            var = &gapo;
        else if ( strcmp(a, "-sa") == 0 )
            var = &sa;
        else if ( strcmp(a, "-sb") == 0 )
            var = &sb;        
        else if ( strcmp(a, "-p") == 0 )
            var = &threads;
        
        if ( var != 0 ) {
            if ( i == ac )
                error("Missing an argument after %s.", a);
            *var = av[i++];
            continue;
        }
        error("Unknown argument: %s.", a);
        return 1;               
    }

    if ( args.genepred_fname == NULL )
        error("Specify genepred database with -data.");

    if ( args.refrna_fname == NULL )
        error("Specify RNA sequence database with -rna.");

    if ( args.reference_fname == NULL )
        error("Specify reference genome database with -ref.");


    args.rna_fai = fai_load(args.refrna_fname);
    if ( args.rna_fai == NULL )
        error("Failed to load fai index of %s.", args.refrna_fname);

    args.ref_fai = fai_load(args.reference_fname);
    if ( args.ref_fai == NULL )
        error("Failed to load fai index of %s.", args.reference_fname);

    if ( args.format ) {
        if ( strcmp(args.format, "genepred") == 0 )
            set_format_genepred();
        else if ( strcmp(args.format, "refgene") == 0 )
            set_format_refgene();
        else if ( strcmp(args.format, "refflat") == 0 )
            set_format_refflat();
        else
            error("Unknown format %s.", args.format);
    }

    if ( sa != NULL ) {
        int temp = str2int((char*)sa);
        args.sa = temp > 0 ?  temp : args.sa;
    }

    if ( sb != NULL ) {
        int temp = str2int((char*)sb);
        args.sb = temp > 0 ?  temp : args.sb;
    }

    if ( gape != NULL ) {
        int temp = str2int((char*)gape);
        args.gape = temp > 0 ?  temp : args.gape;
    }

    if ( gapo != NULL ) {
        int temp = str2int((char*)gapo);
        args.gapo = temp > 0 ?  temp : args.gapo;
    }

    /* if ( threads != NULL ) { */
    /*     args.threads = str2int((char*)threads); */
    /*     if ( args.threads < 1 ) */
    /*         args.threads = 1; */
    /* } */
    // Since BGZF is NOT thread safe, here only accept thread == 1.
    // init memory pool for threads
    args.pbuf = (struct buffer*)calloc(args.threads, sizeof(struct buffer));
    for ( i = 0; i < args.threads; ++i ) {
        args.pbuf[i].n = args.pbuf[i].m = 0;
        args.pbuf[i].buffer = NULL;
    }

    int j, k;
    for ( i = k = 0; i < 4; ++i ) {
        for ( j = 0; j < 4; ++j )
            args.mat[k++] = i == j ? args.sa : -args.sb;
        args.mat[k++] = 0; 
    }
    for ( j = 0; j < 5; ++j ) args.mat[k++] = 0;

    return 0;
}

struct data {
    int n_buffers;
    faidx_t *rna_fai;
    faidx_t *ref_fai;
    struct args *args;
    kstring_t *buffers;
};

void kt_for(int n_threads, void (*func)(void*, long, int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

kstring_t *read_buffer(kstream_t *ks, int chunk_size, int *_n)
{
    kstring_t *buffers = NULL;
    kstring_t str = {0, 0, 0};
    int ret;
    int n, m;
    m = n = 0;
    int size = 0;
    int dret;
    for ( ;; ) {
        //ret = bgzf_getline(fp, '\n', &str);
        ret = ks_getuntil2(ks, KS_SEP_LINE, &str, &dret, 0);
        if ( ret < 0 )
            break;
        if ( n == m ) {
            m = m ? m<<1 : 256;
            buffers = (kstring_t*)realloc(buffers, m*sizeof(kstring_t));
        }
        kstring_t *s = &buffers[n];
        s->l = s->m = str.l;
        s->s = strndup(str.s, str.l);
        size += buffers[n++].l;
        if ( size > chunk_size )
            break;
    }
    *_n = n;
    if ( str.m )
        free(str.s);
    
    return buffers;
}

int process(struct args *args, struct data *data, kstring_t *str)
{
    // extern function from faidx_def.c
    // retrieve transcript version number from transcript reference database
    // the format of transcript title in FASTA should be format as >TRANSCRIPT[space]VERSION 
    extern int trans_retrieve_version(void *_fai, const char *trans);
    
    static const int sntab[256] = {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    };
   
    struct genepred_line *line = genepred_line_create();
    parse_line(str, line);
    int len;
    int ver = trans_retrieve_version((void*)data->rna_fai, line->name1);
    char *trans = faidx_fetch_seq((const faidx_t*)data->rna_fai, line->name1, 0, 10000, &len);

    if ( trans == 0 || len == 0 ) {
        genepred2line(line, str);
        ksprintf(str, "\t%d\t*", ver);
        return 1;
    }
    
    int i;
    kstring_t tmp = { 0, 0, 0, };
    for ( i = 0; i < line->exon_count; ++i ) {
        int l = 0;

        char *seq = faidx_fetch_seq((const faidx_t*)data->ref_fai, line->chrom, line->exons[BLOCK_START][i]-1, line->exons[BLOCK_END][i]-1, &l);
        if ( seq == NULL || l <= 0 ) {
            tmp.l = 0;
            break;
        }
            
        if ( l != line->exons[BLOCK_END][i] - line->exons[BLOCK_START][i]+1)
            error("Inconsistant block length.");
        if ( seq == NULL ) {
            warnings("Failed to fetch sequence from reference. %s:%d-%d.",line->chrom, line->exons[BLOCK_START][i], line->exons[BLOCK_END][i]);
            break;
        }
        kputs(seq, &tmp);
        free(seq);        
    }

    if ( tmp.l == 0 ) {
        genepred2line(line, str);
        ksprintf(str, "\t%d\t*", ver);
        return 1;
    }
    //debug_print("%s\n%s\n%s\n", line->name1, trans, tmp.s);
    for ( i = 0; i < tmp.l; ++i )
        tmp.s[i] = (uint8_t)sntab[tmp.s[i]];
        
    if ( line->strand == '-') {
        for (i = 0; i < tmp.l>>1; ++i ) {
            int t = tmp.s[i] == 4 ? 4 : 3 - tmp.s[i]; tmp.s[i] = 3 - tmp.s[tmp.l-i-1]; tmp.s[tmp.l-i-1] = t;
        }
        if ( tmp.l % 2 )
            tmp.s[tmp.l/2+1] = 3 - tmp.s[tmp.l/2+1];
    }
    
    for ( i = 0; i < len; ++i )
        trans[i] = (uint8_t)sntab[trans[i]];
    
    uint32_t *cigar = NULL;
    int n_cigar = 0;
    int score = ksw_global(len, (uint8_t*)trans, tmp.l, (uint8_t*)tmp.s, 5, args->mat, args->gapo, args->gape, 100, &n_cigar, &cigar);
    int ret;
    genepred2line(line, str);
    kputc('\t', str); kputw(ver, str); kputc('\t', str);
    if ( n_cigar && score > 0 ) {	
        for ( i = 0; i < n_cigar; ++i ) { 
            int c = cigar[i]&0xf;
            ksprintf(str, "%d%c", cigar[i]>>4, "MIDSH"[c]);
        }
        ret = 0;
    } else {
        kputs("*", str);
	warnings("%s not properly checked.", line->name1);
        ret = 1;
    }
    
    free(cigar);
    free(trans);
    free(tmp.s);
    return ret;
}

static void worker_for(void *_data, long i, int tid)
{
    struct data *data = (struct data*)_data;
    process(data->args, data, &data->buffers[i]);
}

static void *worker_pipeline(void *shared, int step, void *_data)
{
    int i;
    struct args *args = (struct args *)shared;
    if ( step == 0 ) {
        struct data *data;
        data = (struct data*)malloc(sizeof(struct data));
        data->buffers = read_buffer(args->ks, args->chunk_size, &data->n_buffers);
        data->args = args;
        data->rna_fai = fai_load(args->refrna_fname);
        data->ref_fai = fai_load(args->reference_fname);
        if ( data->n_buffers )
            return data;
        else
            free(data);
    } else if ( step == 1 ) {
        struct data *data = (struct data*)_data;
        kt_for(args->threads, worker_for, data, data->n_buffers);
        return data;
    } else if ( step == 2 ) {
        struct data *data = (struct data*)_data;
        for ( i = 0; i < data->n_buffers; ++i ) {
            puts(data->buffers[i].s);
            free(data->buffers[i].s);
        }
        data->n_buffers = 0;
        fai_destroy(data->rna_fai);
        fai_destroy(data->ref_fai);
        free(data->buffers);
        free(data);
    }
    return 0;
}

int realign_trans()
{    
    gzFile fp = gzopen(args.genepred_fname, "r");
    if ( fp == NULL )
        error("%s : %s.", args.genepred_fname, strerror(errno));

    args.ks = ks_init(fp);
    kt_pipeline(2, worker_pipeline, &args, 3);
    
    ks_destroy(args.ks);
    gzclose(fp);
    fai_destroy(args.rna_fai);
    fai_destroy(args.ref_fai);
    return 0;
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    if ( realign_trans() )
        return 1;
    // release_memory();

    return 0;
}
