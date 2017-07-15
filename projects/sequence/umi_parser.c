#include "utils.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include "sequence.h"
#include "number.h"
#include <zlib.h>
#include <string.h>
#include "pkg_version.h"

KSEQ_INIT(gzFile, gzread)
struct pair {
    int id;
    int start;
    int end;
};
struct args {
    const char *input1_fname;
    const char *input2_fname;
    const char *trim_reg;
    const char *umi_reg;
    const char *output1_fname;
    const char *output2_fname;
    kstring_t str1;
    kstring_t str2;
    struct pair trim;
    struct pair umi;
} args = {
    .input1_fname = NULL,
    .input2_fname = NULL,
    .trim_reg = NULL,
    .umi_reg = NULL,
    .output1_fname = NULL,
    .output2_fname = NULL,
    .str1 = {0, 0, 0},
    .str2 = {0, 0, 0},
    .trim = {0, 0, 0},
    .umi = {0, 0, 0},
};

int usage()
{
    fprintf(stderr,
            "Trim primer sequences and barcode sequences in fastq files and rename read names with UMI id.\n"
            "Usage : \n"
            "UMI_parser -umi 2:101-113 read1.fq.gz [read2.fq.gz]\n"
            "   -umi   2:101-106           // UMI barcode region in read sequence, format is [1|2]:start-end\n"
            //"   -trim  1:1-6               // primer region, this region will be trimmed\n"
            "   -t     INT                 // number of threads\n"
            "   -out1  FILE                // output file for read1\n"
            "   -out2  FILE                // output file for read2\n"
            "Version : %s\n"
            "Homepage : https://github.com/shiquan/small_projects\n",
            PROJECTS_VERSION
        );
    return 1;
}

int parse_reg(const char *str, int l, struct pair *pair)
{
    if ( str[0] == '1'  && str[1] == ':') {
        pair->id = 1;
    } else if ( str[0] == '2' && str[1] == ':') {
        pair->id = 2;
    } else {
        error_print("%s does not like a region.", str);
        return 1;
    }
    
    int i;
    for ( i = 2; i < l; ++i )
        if ( str[i] == '-')
            break;
    pair->start = str2int_l((char*)str+2, i-2);
    pair->end = str2int_l((char*)str+i+1, l-i);
    if ( pair->start < 1 || pair->end < 1 ) {
        error_print("Unsupport region, %d - %d.", pair->start, pair->end);
        return 1;
    }
        
    return 0;
}
    
int parse_args(int ac, char **av)
{
    if ( ac == 1 )
        return usage();
    
    int i;
    const char *threads = NULL;
    for ( i = 1; i < ac; ) {
        const char *a = av[i++];
        const char **var = 0;
        if ( strcmp(a, "-h") == 0 )
            return usage();

        if ( strcmp(a, "-umi") == 0 && args.umi_reg == NULL )
            var = &args.umi_reg;
        else if ( strcmp(a, "-trim") == 0 && args.trim_reg == NULL )
            var = &args.trim_reg;
        else if ( strcmp(a, "-t") == 0 )
            var = &threads;
        else if ( strcmp(a, "-out1") == 0 && args.output1_fname == NULL )
            var = &args.output1_fname;
        else if ( strcmp(a, "-out2") == 0 && args.output2_fname == NULL )
            var = &args.output2_fname;

        if ( var != 0 ) {
            if ( i == ac ) {
                error("Miss an argument after %s.", a);
                return 1;
            }
            *var = av[i++];
            continue;
        }

        if ( args.input1_fname == NULL )
            args.input1_fname = a;
        else if ( args.input2_fname == NULL )
            args.input2_fname = a;
        else
            error("Unknown argument : %s.", a);
    }

    if ( args.input1_fname == NULL )
        error("No sequence file specified.");

    if ( args.umi_reg == NULL )
        error("No UMI region specified.");

    if ( args.output1_fname == NULL ) {
        ksprintf(&args.str1, "UMI_%s", args.input1_fname);
        args.output1_fname = (const char*)args.str1.s;
    }

    if ( args.output2_fname == NULL && args.input2_fname != NULL ) {
        ksprintf(&args.str2, "UMI_%s", args.input2_fname);
        args.output2_fname = (const char*)args.str2.s;
    }

    if ( args.trim_reg && parse_reg(args.trim_reg, strlen(args.trim_reg), &args.trim) == 1 )
        return 1;
    if ( args.umi_reg && parse_reg(args.umi_reg, strlen(args.umi_reg), &args.umi) == 1)
        return 1;
    
    return 0;
}

int parse_UMI()
{
    // check UMI regions
    gzFile fp1, fp2;
    kseq_t *seq1 = NULL;
    kseq_t *seq2 = NULL;
    int l1, l2;
    int i;
    BGZF *out1 = NULL;
    BGZF *out2 = NULL;
    
    fp1 = gzopen(args.input1_fname, "r");
    if ( fp1 == NULL )
        error("%s : %s.", args.input1_fname, strerror(errno));
    out1 = bgzf_open(args.output1_fname, "w");
    if ( out1 == NULL )
        error("%s : %s.", args.output1_fname, strerror(errno));
    
    seq1 = kseq_init(fp1);

    if ( args.input2_fname != NULL ) {
        fp2 = gzopen(args.input2_fname, "r");
        if (fp2 == NULL )
            error("%s : %s.", args.input2_fname, strerror(errno));
        seq2 = kseq_init(fp2);
        out2 = bgzf_open(args.output2_fname, "w");
        if ( out2 == NULL )
            error("%s : %s.", args.output2_fname, strerror(errno));
    }

    int length = args.umi.end - args.umi.start +1;
    
    if ( seq2 == NULL ) {
        if ( args.umi.id == 2)
            error("Inconsistant UMI region. %s", args.umi_reg);
        
        if ( args.trim.id == 2 )
            error("Inconsistant trim region. %s", args.trim_reg);
        
        kstring_t string = {0, 0, 0};
        
        while ( (l1 = kseq_read(seq1)) >= 0 ) {
            string.l = 0;
            if ( seq1->qual.l == 0)
                error("Only support FASTQ for now.");
            
            kputc('@', &string);
            kputs(seq1->name.s, &string);
            if ( string.s[string.l-2] == '/' && string.s[string.l-1] == '1')
                string.l -= 2;
            kputsn(seq1->seq.s+args.umi.start-1, length, &string);
            kputc('\n', &string);
            if ( args.umi.start == 1 ) 
                kputs(seq1->seq.s+args.umi.end, &string);
            else
                kputsn(seq1->seq.s, args.umi.start-1, &string);
            kputs("\n+\n", &string);
            if ( args.umi.start == 1 ) 
                kputs(seq1->qual.s+args.umi.end, &string);
            else
                kputsn(seq1->qual.s, args.umi.start-1, &string);
            kputc('\n', &string);            
        
            if ( bgzf_write(out1, string.s, string.l) != string.l )
                error("Write error : %d", out1->errcode);
        }

        if ( string.m )
            free(string.s);
        
    } else {
        kstring_t str1 = { 0, 0, 0};
        kstring_t str2 = { 0, 0, 0};
        
        do {
            l1 = kseq_read(seq1);
            l2 = kseq_read(seq2);
            if ( l1 < 0 || l2 < 0 )
                break;
            
            for ( i = 0; i < seq1->name.l-1; ++i )
                if ( seq1->name.s[i] != seq2->name.s[i])
                    error("Inconsistant read name. %s vs %s.", seq1->name.s, seq2->name.s);
            
            str1.l = 0;
            str2.l = 0;

            if ( args.umi.id == 1 ) {
                kputc('@', &str1);
                kputs(seq1->name.s, &str1);
                if ( str1.s[str1.l-2] == '/' && str1.s[str1.l-1] == '1')
                    str1.l -= 2;
                kputs("_UID:",&str1);
                kputsn(seq1->seq.s+args.umi.start-1, length, &str1);
                kputc('\n', &str1);
                if ( args.umi.start == 1 ) 
                    kputs(seq1->seq.s+args.umi.end, &str1);
                else
                    kputsn(seq1->seq.s, args.umi.start-1, &str1);
                kputs("\n+\n", &str1);
                if ( args.umi.start == 1 ) 
                    kputs(seq1->qual.s+args.umi.end, &str1);
                else
                    kputsn(seq1->qual.s, args.umi.start-1, &str1);
                kputc('\n', &str1);            

                kputc('@', &str2);
                kputs(seq2->name.s, &str2);
                if ( str2.s[str2.l-2] == '/' && str2.s[str2.l-1] == '2')
                    str2.l -= 2;
                kputsn(seq1->seq.s+args.umi.start-1, length, &str2);
                kputc('\n', &str2);
                ksprintf(&str2, "%s\n+\n%s\n", seq2->seq.s, seq2->qual.s);
                
            } else {
                kputc('@', &str1);
                kputs(seq1->name.s, &str1);
                if ( str1.s[str1.l-2] == '/' && str1.s[str1.l-1] == '1')
                    str1.l -= 2;
                kputs("_UID:",&str1);
                kputsn(seq2->seq.s+args.umi.start-1, length, &str1);
                kputc('\n', &str1);
                ksprintf(&str1, "%s\n+\n%s\n", seq1->seq.s, seq1->qual.s);

                kputc('@', &str2);
                kputs(seq2->name.s, &str2);
                if ( str2.s[str2.l-2] == '/' && str2.s[str2.l-1] == '2')
                    str2.l -= 2;
                kputs("_UID:",&str2);
                kputsn(seq2->seq.s+args.umi.start-1, length, &str2);
                kputc('\n', &str2);
                if ( args.umi.start == 1 ) 
                    kputs(seq2->seq.s+args.umi.end, &str2);
                else
                    kputsn(seq2->seq.s, args.umi.start-1, &str2);
                kputs("\n+\n", &str2);
                if ( args.umi.start == 1 ) 
                    kputs(seq2->qual.s+args.umi.end, &str2);
                else
                    kputsn(seq2->qual.s, args.umi.start-1, &str2);
                kputc('\n', &str2);
            }
            if ( bgzf_write(out1, str1.s, str1.l) != str1.l )
                error("Write error : %d", out1->errcode);

            if ( bgzf_write(out2, str2.s, str2.l) != str2.l )
                error("Write error : %d", out2->errcode);

        } while(1);
        if ( str1.m )
            free(str1.s);
        if ( str2.m)
            free(str2.s);
    }

    bgzf_close(out1);
    bgzf_close(out2);
    if ( args.str1.m)
        free(args.str1.s);
    if ( args.str2.m)
        free(args.str2.s);

    kseq_destroy(seq1);
    if ( seq2 )
        kseq_destroy(seq2);
    return 0;
}
int main (int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    if ( parse_UMI() )
        return 1;

    return 0;
}
