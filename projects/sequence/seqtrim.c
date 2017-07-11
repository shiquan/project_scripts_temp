#include "utils.h"
#include <string.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <zlib.h>
#include "sequence.h"
#include "pkg_version.h"

KSEQ_INIT(gzFile, gzread)

struct args {
    const char *input_fname;
    int trim_start; // the start location of the sequences
    int trim_end; // the end location of the sequences
    int print_title;
    int compl;
} args = {
    .input_fname = 0,
    .trim_start = 0,
    .trim_end = 0,
    .print_title = 1,
    .compl = 0,
};

int usage()
{
    fprintf(stderr, "Trim FASTA and FASTQ file with predefined range.\n");
    fprintf(stderr, "Usage :\nseqtrim [options] in.fasta.gz|in.fastq.gz\n");
    fprintf(stderr, "    -start loc    // start location of the sequences, default is the first of the sequences.\n");
    fprintf(stderr, "    -end loc      // end location of the sequences, default is the end of the sequences.\n");
    fprintf(stderr, "    -seq          // only export the sequences. no titles.\n");
    fprintf(stderr, "    -comp         // export complement strand of sequences\n");
    fprintf(stderr, "\nVersion: %s\n", PROJECTS_VERSION);
    fprintf(stderr, "Homepage: https://github.com/shiquan/small_projects\n");
    return 1;
};
int parse_args(int argc, char **argv)
{
    if ( argc == 1)
        return usage();
    int i;
    const char *start = 0;
    const char *end = 0;
    for (i = 1; i < argc; )  {
        const char *a = argv[i++];
        if ( strcmp(a, "-h") == 0)
            return usage();
        const char **var = 0;
        if ( strcmp(a, "-start") == 0 ) 
            var = &start;
        else if ( strcmp(a, "-end") == 0 )
            var = &end;
        else if ( strcmp(a, "-seq") == 0 ) {            
            args.print_title = 0;
            continue;
        }
        else if ( strcmp(a, "-comp") == 0 ) {
            args.compl = 1;
            continue;
        }

        if ( var != 0) {
            if ( i == argc )
                error("Missing argument after %s", a);
            *var = argv[i++];
            continue;
        }

        if ( args.input_fname == 0 ) {
            args.input_fname = a;
            continue;
        }
        error("Unknown argument, %s.", a);
    }
    if ( args.input_fname == 0 && (!isatty(fileno(stdin))) )
        args.input_fname = "-";
    if ( args.input_fname == 0 )
        error("No input file.");
    if ( start != NULL) {
        args.trim_start = atoi(start);
        if (args.trim_start < 0)
            error("Bad start value, should > 0, %d", args.trim_start);
    }

    if ( end != NULL) {
        args.trim_end = atoi(end);
        if ( args.trim_end <= 0 )
            error("Bad end value, should > 0, %d", args.trim_end);
    }
    if (args.trim_end && args.trim_start >= args.trim_end)
        error("Should set location start smaller than end.");
    
    return 0;
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    gzFile fp;
    kseq_t *seq;
    fp = strcmp(args.input_fname, "-") == 0 ? gzdopen(fileno(stdin), "r") : gzopen(args.input_fname, "r");
    if ( fp == 0 )
        error("%s : %s", args.input_fname, strerror(errno));
    seq = kseq_init(fp);

    int l;
    while ((l = kseq_read(seq)) >= 0) {
        if ( args.trim_start && args.trim_start > seq->seq.l)
            continue;
        
        if ( args.print_title ) {
            putchar(seq->qual.l ? '@' : '>');
            fputs(seq->name.s, stdout);
            putchar('\n');
        }

        if ( args.trim_end && args.trim_end < seq->seq.l ) {
            seq->seq.s[args.trim_end] = '\0';
            seq->seq.l = args.trim_end;
            if ( seq->qual.l ) {
                seq->qual.s[args.trim_end] = '\0';
                seq->qual.l = args.trim_end;
            }
        }

        if ( args.trim_start ) {
            if ( args.compl == 1 ) {
                int i;
                for ( i = seq->seq.l -1 ; i >= args.trim_start-1; --i )
                    putchar("TGCAN"[seq2code4(seq->seq.s[i])]);
            } else {
                fputs(seq->seq.s + args.trim_start -1, stdout);
            }
            putchar('\n');
            if ( args.print_title && seq->qual.l) {
                putchar('+'); putchar('\n');
                if ( args.compl == 1 ) {
                    int i;
                    for ( i = seq->qual.l - 1; i >= args.trim_start-1; --i)
                        putchar(seq->qual.s[i]);
                } else {
                    fputs(seq->qual.s + args.trim_start -1, stdout);
                }
                putchar('\n');
            }
        } else {
            if ( args.compl == 1 ) {
                int i;
                for ( i = seq->seq.l -1 ; i >= args.trim_start-1; --i )
                    putchar("TGCAN"[seq2code4(seq->seq.s[i])]);            
            } else {
                fputs(seq->seq.s, stdout);
            }
            putchar('\n');
            if ( args.print_title && seq->qual.l) {
                putchar('+'); putchar('\n');
                fputs(seq->qual.s, stdout);
                putchar('\n');
            }            
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    return 0;
}
