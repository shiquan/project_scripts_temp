#include "utils.h"
#include "number.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "htslib/bgzf.h"
#include "fastq.h"
#include "pkg_version.h"

KSEQ_INIT(gzFile, gzread)

int usage()
{
    fprintf(stderr,
            "Trim adaptor pollution sequence and rename read id with UID.\n"
            "dyncut_adaptor [options] read1.fq [read2.fq]\n"
            "\n"
            "Options :\n"
            "    -adaptor      adaptor pollution sequences for read 1, should be complementary with R2 sequencing primer in BGISEQ500 platform\n"
            "    -barcode      barcode file, consist of the barcode name and the sequence columns\n"
            "    -out          output directory\n"
            "    -mis_ada      mismatch allowed in the adaptor sequence alignment\n"
            "    -mis_bar      mismatch allowed in the barcode sequence alignment [0-3]\n"
            "    -min_length   minimual sequence length for trimmed reads\n"
            "    -report       specify report file instead of print to stdout\n"
            "    -uid          rename read id with extra barcode tag\n"
            "    -dropr2       drop R2 fastq file adaptor pollution is detected in R1\n"
            "    -trim         trim ends even if partly adaptor detected, conflict with -uid and -barcode\n"
            "\n"
            "Version : %s"
            "Homepage : https://github.com/shiquan/small_projects\n",
            PROJECTS_VERSION            
        );
    return 1;
}

struct report { int n, m; uint64_t *a; } report = { 0, 0, 0, };

struct args {
    const char *adaptor;
    const char *barcode_fname;
    const char *output_dir;
    int adaptor_length;
    int mis_ada;
    int mis_bar;
    int minimual_length;
    const char *report_fname;
    int rename_uid_flag;
    int trim_tails_flag;
    int drop_read2;
    const char *read1_file;
    const char *read2_file;
    BGZF *failed_1;
    BGZF *failed_2;
    const char *report;
    struct barcode barcode;    
} args = {
    .adaptor = NULL,
    .barcode_fname = NULL,
    .output_dir = NULL,
    .adaptor_length = 0,
    .mis_ada = 0,
    .mis_bar = 0,
    .report_fname = NULL,
    .rename_uid_flag = 0,
    .trim_tails_flag = 0,
    .read1_file = NULL,
    .read2_file = NULL,
    .drop_read2 = 0,
    .failed_1 = NULL,
    .failed_2 = NULL,
    .report = NULL,
    .barcode = { 0, 0, 0,},
};

int parse_args(int argc, char **argv)
{
    if ( argc == 1 )
        return usage();

    const char *mis_ada = NULL;
    const char *mis_bar = NULL;
    const char *minimual_length = NULL;
    int i;
    for ( i = 1; i < argc; ) {
        const char *a = argv[i++];
        if ( strcmp(a, "-h") == 0 )
            return usage();

        const char **var = 0;
        if ( strcmp(a, "-adaptor") == 0 && args.adaptor == NULL )
            var = &args.adaptor;
        else if ( strcmp(a, "-barcode") == 0 && args.barcode_fname == NULL )
            var = &args.barcode_fname;
        else if ( (strcmp(a, "-output") == 0 || strcmp (a, "-out") == 0 )&& args.output_dir == NULL )
            var = &args.output_dir;
        else if ( strcmp(a, "-mis_ada") == 0 && mis_ada == NULL )           
            var = &mis_ada;
        else if ( strcmp(a, "-mis_bar") == 0 && mis_bar == NULL )
            var = &mis_bar;
        else if ( strcmp(a, "-min_length") == 0 && minimual_length == NULL )
            var = &minimual_length;
        else if ( strcmp(a, "-report") == 0 && args.report == NULL )
            var = &args.report;

        if ( var != 0 ) {
            if ( i == argc ) {
                error_print("Miss an argument after %s.", a);
                return 1;
            }
            *var = argv[i++];
            continue;
        }
        if ( strcmp(a, "-uid") == 0 ) {
            args.rename_uid_flag = 1;
            continue;
        } else if ( strcmp(a, "-dropr2") == 0 ) {
            args.drop_read2 = 1;
            continue;
        }
        
        if ( args.read1_file == NULL )
            args.read1_file = a;
        else if ( args.read2_file == NULL )
            args.read2_file = a;
        else
            error("Unknown argument, %s.", a);        
    }

    if ( args.read1_file == NULL )
        error("No sequence file specified.");

    if ( args.adaptor == NULL )
        error("Please specify adaptor sequence with -adaptor");

    args.adaptor_length = strlen(args.adaptor);
    if ( check_acgt(args.adaptor, args.adaptor_length) )
        error("%s looks not like an adaptor sequence.", args.adaptor);
    
    if ( mis_ada ) {
        args.mis_ada = str2int((char*)mis_ada);
        if ( args.mis_ada < 0 )
            args.mis_ada = 0;
    }

    if ( mis_bar ) {
        args.mis_bar = str2int((char*)mis_bar);
        if ( args.mis_bar < 0 )
            args.mis_bar = 0;
    }

    if ( minimual_length )
        args.minimual_length = str2int((char*)minimual_length);
    
    if ( load_barcode_file(args.barcode_fname, &args.barcode) ) {
        error_print("Failed to load barcode file.");
        return 1;
    }    
    return 0;
}

#define STR_INIT {0, 0, 0}

// trim adaptor pollution sequences and export read1 fastq file into output directory
int trim_adaptor()
{    
    int l1, l2;
    int i;
    gzFile fp1, fp2;
    kseq_t *seq1 = NULL;
    kseq_t *seq2 = NULL;
    int file_is_fastq;
    file_is_fastq = check_file_is_fastq(args.read1_file);
    if ( file_is_fastq == -1)
        return 1;

    fp1 = gzopen(args.read1_file, "r");
    if ( fp1 == 0 )
        error ("%s : %s", args.read1_file, strerror(errno));
    seq1 = kseq_init(fp1);
    
    if ( args.read2_file != NULL) {
        int check_file = check_file_is_fastq(args.read2_file);
        if ( check_file != file_is_fastq )
            error("Inconsistant file type. %d vs %d.", file_is_fastq, check_file );        
            
        fp2 = gzopen(args.read2_file, "r");
        if ( fp2 == NULL)
            error("%s : %s.", args.read2_file, strerror(errno));
        
        seq2 = kseq_init(fp2);    
    }

    kstring_t string = STR_INIT;
    for ( i = 0; i < args.barcode.n; ++i ) {
        struct name *name = &args.barcode.names[i];
        string.l = 0;
        if ( args.output_dir )
            ksprintf(&string,"%s/%s_1.%s.gz", args.output_dir, name->name, file_is_fastq == 0 ? "fq" : "fa");
        else
            ksprintf(&string,"%s_1.%s.gz",name->name, file_is_fastq == 0 ? "fq" : "fa");

        name->fp1 = bgzf_open(string.s, "w");
        if ( name->fp1 == NULL )
            error("%s : %s.", string.s, strerror(errno));
        name->fp2 = NULL;
        if ( args.read2_file != NULL && args.drop_read2 == 0) {
            string.l = 0;
            if ( args.output_dir )
                ksprintf(&string,"%s/%s_2.%s.gz", args.output_dir, name->name, file_is_fastq == 0 ? "fq" : "fa");
            else
                ksprintf(&string,"%s_2.%s.gz",name->name, file_is_fastq == 0 ? "fq" : "fa");
            name->fp2 = bgzf_open(string.s, "w");
            if ( name->fp2 == NULL )
                error("%s : %s.", string.s, strerror(errno));            
        }
    }
    do {
            string.l = 0;
            if ( args.output_dir )
                ksprintf(&string, "%s/untrimmed_1.%s.gz", args.output_dir, file_is_fastq == 0 ? "fq" : "fa");
            else
                ksprintf(&string, "untrimmed_1.%s.gz", file_is_fastq == 0 ? "fq" : "fa");
            args.failed_1 = bgzf_open(string.s, "w");
            if ( args.failed_1 == NULL )
                error("%s : %s.", string.s, strerror(errno));
            
            if ( args.read2_file == NULL || args.drop_read2 == 0)
                break;
            string.l = 0;
            if ( args.output_dir )
                ksprintf(&string, "%s/untrimmed_2.%s.gz", args.output_dir, file_is_fastq == 0 ? "fq" : "fa");
            else
                ksprintf(&string, "untrimmed_2.%s.gz", file_is_fastq == 0 ? "fq" : "fa");
            args.failed_2 = bgzf_open(string.s, "w");
            if ( args.failed_2 == NULL )
                error("%s : %s.", string.s, strerror(errno));
    } while ( 0 );
    
    if ( args.read2_file == NULL ) {
        // TODO : improve the efficient of this function.
        // complexity of this function is O(n*k*(m-l)), n is reads count, k is the barcode count, m is the read length, and  l is the adaptor length
        BGZF *fp = NULL;
        do {
            l1 = kseq_read(seq1);
            if ( l1 < 0 )
                break;
            if ( l1 == 0 )
                continue;
            int check_length = l1 - args.adaptor_length;
            string.l = 0;
            for ( i = 0; i <  check_length; ++i ) {
                
                if ( check_match2(seq1->seq.s+i, args.adaptor, args.mis_ada, args.adaptor_length) < 0 ) { 
                    continue;
                } else { // check the barcode, failed barcode reads will also export to failed fastqs.
                    if ( args.minimual_length &&  i < args.minimual_length )
                        goto skip_record;

                    int j;
                    for ( j = 0; j < args.barcode.n; ++j ) {
                        struct name *name = &args.barcode.names[j];
                        int l = strlen(name->barcode);
                        if ( l1 < i + l + args.adaptor_length ) { // export to failed fastqs
                            continue;
                        }                             
                        int m = check_match2(seq1->seq.s+i+args.adaptor_length, name->barcode, args.mis_bar, l);
                        if ( m == -1 )
                            continue;
                        
                        // export trimmed fastqs
                        if ( args.rename_uid_flag == 1 ) {
                            if ( file_is_fastq == 0 ) {
                                kputc('@', &string);
                                kputs(seq1->name.s, &string);
                                if ( string.s[string.l-2] == '/' && string.s[string.l-1] == '1')
                                    string.l -= 2;
                                kputs("_UID:",&string);
                                kputsn(seq1->seq.s+i+args.adaptor_length, l, &string);
                                kputc('\n', &string);
                                kputsn(seq1->seq.s, i, &string);                            
                                kputs("\n+\n", &string);
                                kputsn(seq1->qual.s, i, &string);
                                kputc('\n', &string);
                            } else {
                                kputc('>', &string);
                                kputs(seq1->name.s, &string);
                                if ( string.s[string.l-2] == '/' && string.s[string.l-1] == '1')
                                    string.l -= 2;
                                kputs("_UID:",&string);
                                kputsn(seq1->seq.s+i+args.adaptor_length, l, &string);
                                kputc('\n', &string);
                                kputsn(seq1->seq.s, i, &string);                            
                                kputc('\n', &string);                                
                            }
                        } else {
                            if ( file_is_fastq == 0 ) {
                                kputc('@', &string);
                                kputs(seq1->name.s, &string);
                                kputc('\n', &string);
                                kputsn(seq1->seq.s, i, &string);
                                kputs("\n+\n", &string);
                                kputsn(seq1->qual.s, i, &string);
                                kputc('\n', &string);
                            } else {
                                kputc('>', &string);
                                kputs(seq1->name.s, &string);
                                kputc('\n', &string);
                                kputsn(seq1->seq.s, i, &string);
                                kputc('\n', &string);
                            }
                        } // end parse uid
                        fp = name->fp1;
                        break;                     
                    }
                    if ( j == args.barcode.n ) {
                        fp = args.failed_1;
                        break;
                    }
                }
                break;
            }
            if ( i == check_length ) {
                if ( file_is_fastq == 0 ) {
                    ksprintf(&string, "@%s\n%s\n+\n%s\n", seq1->name.s, seq1->seq.s, seq1->qual.s);                
                } else {
                    ksprintf(&string, ">%s\n%s\n", seq1->name.s, seq1->seq.s);
                }
                fp = args.failed_1;
            }
            if ( fp == NULL )
                error("Error file handler. Please report this bug to developer.");
            
            if ( bgzf_write(fp, string.s, string.l ) != string.l )
                error ("Write error : %d", fp->errcode);

          skip_record:
            fp = NULL;
            
        } while (1);
        
    } else {
        kstring_t str1 = STR_INIT;
        kstring_t str2 = STR_INIT;
        BGZF *fp1 = NULL, *fp2 = NULL;
        do {
            l1 = kseq_read(seq1);
            l2 = kseq_read(seq2);
            if ( l1 < 0 || l2 < 0 )
                break;
            str1.l = 0;
            str2.l = 0;

            for ( i = 0; i < seq1->name.l -1; ++i )
                if ( seq1->name.s[i] != seq2->name.s[i])
                    error("Inconsistant read name. %s vs %s.", seq1->name.s, seq2->name.s);

            // only check adaptor pollution in the read 1, if success trim read 1 and read 2
            int check_length = l1 - args.adaptor_length;            
            for ( i = 0; i < check_length; ++i ) {
                if ( check_match2(seq1->seq.s+i, args.adaptor, args.mis_ada, args.adaptor_length) < 0 ) {
                    continue;
                } else {
                    if ( args.minimual_length &&  i < args.minimual_length )
                        goto skip_record2;

                    int j;
                    for ( j = 0; j < args.barcode.n; ++j ) {
                        struct name *name = &args.barcode.names[j];
                        int l = strlen(name->barcode);
                        if ( l1 < i + l )
                            continue;
                    
                        int m;
                        m = check_match2(seq1->seq.s+i+args.adaptor_length, name->barcode, args.mis_bar, l);
                        if ( m == -1 )
                            continue;

                        // export trimed fastqs
                        if ( args.rename_uid_flag == 1 ) {
                            if ( file_is_fastq == 0 ) {
                                kputc('@', &str1);
                                kputs(seq1->name.s, &str1);
                                if ( str1.s[str1.l-2] == '/' && str1.s[str1.l-1] == '1')
                                    str1.l -= 2;
                                kputs("_UID:",&str1);
                                kputsn(seq1->seq.s+i+args.adaptor_length, l, &str1);
                                kputc('\n', &str1);
                                kputsn(seq1->seq.s, i, &str1);                            
                                kputs("\n+\n", &str1);
                                kputsn(seq1->qual.s, i, &str1);
                                kputc('\n', &str1);
                                if ( args.drop_read2 == 0 ) {
                                    kputc('@', &str2);
                                    kputs(seq1->name.s, &str2);
                                    if ( str2.s[str2.l-2] == '/' && str2.s[str2.l-1] == '2')
                                        str2.l -= 2;
                                    kputs("_UID:",&str2);
                                    kputsn(seq1->seq.s+i+args.adaptor_length, l, &str2);
                                    kputc('\n', &str2);
                                    kputsn(seq1->seq.s, i, &str2);                            
                                    kputs("\n+\n", &str2);
                                    kputsn(seq1->qual.s, i, &str2);
                                    kputc('\n', &str2);
                                }
                            } else {
                                kputc('>', &str1);
                                kputs(seq1->name.s, &str1);
                                if ( str1.s[str1.l-2] == '/' && str1.s[str1.l-1] == '1')
                                    str1.l -= 2;
                                kputs("_UID:",&str1);
                                kputsn(seq1->seq.s+i+args.adaptor_length, l, &str1);
                                kputc('\n', &str1);
                                kputsn(seq1->seq.s, i, &str1);                            
                                kputc('\n', &str1);

                                // read 2
                                if ( args.drop_read2 == 0 ) {
                                    kputc('>', &str2);
                                    kputs(seq1->name.s, &str2);
                                    if ( str2.s[str2.l-2] == '/' && str2.s[str2.l-1] == '2')
                                        str2.l -= 2;
                                    kputs("_UID:",&str2);
                                    kputsn(seq1->seq.s+i+args.adaptor_length, l, &str2);
                                    kputc('\n', &str2);
                                    kputsn(seq1->seq.s, i, &str2);                            
                                    kputc('\n', &str2);
                                }
                            }
                        } else {
                            if ( file_is_fastq == 0 ) {
                                kputc('@', &str1);
                                kputs(seq1->name.s, &str1);
                                kputc('\n', &str1);
                                kputsn(seq1->seq.s, i, &str1);
                                kputs("\n+\n", &str1);
                                kputsn(seq1->qual.s, i, &str1);
                                kputc('\n', &str1);

                                if ( args.drop_read2 == 0 ) {
                                    // read 2
                                    kputc('@', &str2);
                                    kputs(seq1->name.s, &str2);
                                    kputc('\n', &str2);
                                    kputsn(seq1->seq.s, i, &str2);
                                    kputs("\n+\n", &str2);
                                    kputsn(seq1->qual.s, i, &str2);
                                    kputc('\n', &str2);
                                }
                            } else {
                                kputc('>', &str1);
                                kputs(seq1->name.s, &str1);
                                kputc('\n', &str1);
                                kputsn(seq1->seq.s, i, &str1);
                                kputc('\n', &str1);

                                // read 2
                                if ( args.drop_read2 == 0 ) {
                                    kputc('>', &str2);
                                    kputs(seq1->name.s, &str2);
                                    kputc('\n', &str2);
                                    kputsn(seq1->seq.s, i, &str2);
                                    kputc('\n', &str2);
                                }
                            }                        
                        } // end uid parse
                        fp1 = name->fp1;

                        if ( args.drop_read2  == 0 )
                            fp2 = name->fp2;
                        break;                    
                    } // end barcode loop
                    if ( j == args.barcode.n ) { // if no barcode supported
                        i = check_length;
                        break;
                    }
                } // end match
                // i = check_length;
                break;
            }
            if ( i == check_length ) {
                if ( file_is_fastq == 0 ) {
                    ksprintf(&str1, "@%s\n%s\n+\n%s\n", seq1->name.s, seq1->seq.s, seq1->qual.s);
                    ksprintf(&str2, "@%s\n%s\n+\n%s\n", seq2->name.s, seq2->seq.s, seq2->qual.s);                
                } else {
                    ksprintf(&str1, ">%s\n%s\n", seq1->name.s, seq1->seq.s);
                    ksprintf(&str2, ">%s\n%s\n", seq2->name.s, seq2->seq.s);
                }
                fp1 = args.failed_1;
                fp2 = args.failed_2;
            }
            if ( fp1 == NULL)
                error("Error file handler. Please report this bug to developer. %s", seq1->name.s);
            
            if ( bgzf_write(fp1, str1.s, str1.l ) != str1.l )
                error ("Write error : %d", fp1->errcode);            
            if ( fp2 == args.failed_2 && bgzf_write(fp2, str2.s, str2.l ) != str2.l )                
                error ("Write error : %d", fp1->errcode);
          skip_record2:
            fp1 = NULL;
            fp2 = NULL;
        } while (1);
        if ( l1 != l2 )
            error("Inconsistant read records. %d vs %d", l1, l2);

        if ( str1.m )
            free(str1.s);
        if ( str2.m )
            free(str2.s);
    }

    if ( string.m)
        free(string.s);

    kseq_destroy(seq1);
    if ( seq2 )
        kseq_destroy(seq2);
    clean_barcode_struct(&args.barcode);
    if ( args.failed_1 != NULL )
        bgzf_close(args.failed_1);
    if ( args.failed_2 != NULL )
        bgzf_close(args.failed_2);    
    return 0;
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    if ( trim_adaptor() )
        return 1;

    return 0;
}
