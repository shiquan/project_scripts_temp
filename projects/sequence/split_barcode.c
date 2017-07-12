#include "utils.h"
#include <string.h>
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include <zlib.h>
#include "number.h"
#include "pkg_version.h"

KSEQ_INIT(gzFile, gzread)

int usage()
{
    fprintf(stderr,
            "A lightweight program to split barcode reads from FASTQ/FASTA file.\n"
            "Usage :\n"
            "split_barcode -reg 2:101-113 -barcode barcode.txt reads1.fq.gz [reads2.fq.gz]\n"
            "    -reg       // barcode region in read sequence, format is <read 1|2> : <start> - <end>\n"
            "    -comp      // use the complement of barcode sequences\n"
            "    -mismatch  // maximum mismatch tolerant, [1-3]\n"
            "    -barcode   // barcode file, this parameter is mandontory\n"
            "    -out       // output directory.\n"
            "\nAbout the barcode file, it should consist of barcode name and barcode sequences columns, and\n"
            "seperated by tab.\n"
            "Homepage: \n"
            "https://github.com/shiquan/small_projects\n"
        );
    return 1;
}

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
struct args {
    const char *barcode_file;
    const char *output_dir;
    const char *barcode_region;
    const char *read1_file;
    const char *read2_file;
    int mismatch;
    int compl_flag;    
    int start;
    int end;
    struct barcode barcode;
    BGZF *failed_1;
    BGZF *failed_2;
} args = {
    .barcode_file = NULL,
    .output_dir = NULL,
    .barcode_region = NULL,
    .read1_file = NULL,
    .read2_file = NULL,
    .mismatch = 0,
    .compl_flag = 0,
    .barcode = {0, 0, 0},
    .failed_1 = NULL,
    .failed_2 = NULL,
};

static int parse_args(int ac, char **av)
{
    if ( ac == 1 )
        return usage();
    
    int i;
    const char *mismatch = NULL;
    for ( i = 1; i < ac; ) {
        const char *a = av[i++];
        if ( strcmp(a, "-h") == 0 )
            return usage();

        const char **var = 0;
        if ( strcmp(a, "-barcode") == 0 && args.barcode_file == NULL )
            var = &args.barcode_file;
        else if ( strcmp(a, "-out") == 0 && args.output_dir == NULL )
            var = &args.output_dir;
        else if ( strcmp(a, "-reg") == 0 && args.barcode_region == NULL )
            var = &args.barcode_region;
        else if ( strcmp(a, "-mismatch") == 0 && mismatch == NULL )
            var = &mismatch;

        if ( var != 0 ) {
            if ( i == ac ) {
                error("Miss an argument after %s.", a);
                return 1;
            }
            *var = av[i++];
            continue;
        }
        if ( strcmp(a, "-comp") == 0 ) {
            args.compl_flag = 1;
            continue;
        }

        if ( args.read1_file == NULL )
            args.read1_file = a;
        else if ( args.read2_file == NULL )
            args.read2_file = a;
        else
            error("Unknown argument, %s.",a);        
    }

    if ( args.read1_file == NULL )
        error("No sequence file specified.");

    if ( args.barcode_file == NULL )
        error("Barcode index file should be specified by -barcode.");

    if ( args.barcode_region == NULL )
        error("Barcode region should be specified by -reg.");

    if ( mismatch ) {
        args.mismatch = str2int((char*)mismatch);
        if (args.mismatch > 3 && args.mismatch < 1 )
            error("-mismatch should be defined between [1,3].");
    }

    int l;
    l = strlen(args.barcode_region);
    if ( l > 4 && (args.barcode_region[0] == 1 || args.barcode_region[0] == 2) && args.barcode_region[1] == ':') {
        for ( i = 2; i < l; ++i )
            if ( args.barcode_region[i] == '-')
                break;
        args.start = str2int_l((char*)args.barcode_region+2, i - 2);
        args.end = str2int_l((char*)args.barcode_region+i, l -i);
        if ( args.start < 1 || args.end < 1 )
            error("Barcode region unsupport. %d-%d", args.start, args.end);
    } else {
        error("-reg does not look like the supported format. %s", args.barcode_region);
    }
    
    return 0;
}
static int check_match(char *seq1, char *seq2, int mismatch, int length) {
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
static int load_barcode_file(const char *fn, struct barcode *bc)
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
            error("Failed to parse barcode file: %s", string.s);
        
        for ( j = i; j < string.l; ++j ) {
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
            if ( strncmp(bc->names[k].barcode, string.s+i, j-i) == 0 ||
                 strncmp(bc->names[k].name, string.s, i) == 0 ) {
                error_print("Duplicated barcode lines. %s, %s", string.s, string.s+i);
                break;
            }
        }
        if ( k < bc->n ) {            
            name->barcode = strndup(string.s+i, string.l-i);
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

int split_barcode()
{
    if ( load_barcode_file(args.barcode_file, &args.barcode) ) {
        error_print("Failed to load barcode file.");
        return 1;
    }    

    int i;
    int l1, l2;
    int file_is_fastq;

    gzFile fp1, fp2;
    kseq_t *seq1 = NULL;
    kseq_t *seq2 = NULL;
    fp1 = gzopen(args.read1_file, "r");
    if ( fp1 == 0 )
        error("%s : %s", args.read1_file, strerror(errno));

    // check file type
    seq1 = kseq_init(fp1);
    l1 = kseq_read(seq1);
    if ( l1 >= 0 ) {
        file_is_fastq = seq1->qual.l > 0 ? 1 : 0;
        kseq_destroy(seq1);
        gzclose(fp1);
        fp1 = gzopen(args.read1_file, "r");
        seq1 = kseq_init(fp1);
    } else {
        error("%s is empty.", args.read1_file);
    }

    // check read 2
    if ( args.read2_file != NULL ) {
        fp2 = gzopen(args.read2_file, "r");
        if ( fp2 == 0 )
            error("%s : %s", args.read2_file, strerror(errno));
        seq2 = kseq_init(fp2);
        l2 = kseq_read(seq2);
        if ( l2 >= 0 ) {
            if ( file_is_fastq != (seq2->qual.l != 0) )
                error("Inconsistant read type, read1 and read2 must be in same format.");
            kseq_destroy(seq2);
            gzclose(fp2);
            fp2 = gzopen(args.read2_file, "r");
            seq2 = kseq_init(fp2);
        } else {
            error("%s is empty.", args.read2_file);
        }
    }

    // open file handlers
    kstring_t temp = { 0, 0, 0};
    for ( i = 0; i < args.barcode.n; ++i ) {
        struct name *name = &args.barcode.names[i];
        temp.l = 0;
        if ( args.output_dir )
            ksprintf(&temp,"%s/%s_1.%s.gz", args.output_dir, name->name, file_is_fastq ? "fq" : "fa");
        else
            ksprintf(&temp,"%s_1.%s.gz",name->name, file_is_fastq ? "fq" : "fa");

        name->fp1 = bgzf_open(temp.s, "w");
        if ( name->fp1 == NULL )
            error("%s : %s.", temp.s, strerror(errno));
        name->fp2 = NULL;
        if ( seq2 ) {
            temp.l = 0;
            if ( args.output_dir )
                ksprintf(&temp,"%s/%s_2.%s.gz", args.output_dir, name->name, file_is_fastq ? "fq" : "fa");
            else
                ksprintf(&temp,"%s_2.%s.gz",name->name, file_is_fastq ? "fq" : "fa");
            name->fp2 = bgzf_open(temp.s, "w");
            if ( name->fp2 == NULL )
                error("%s : %s.", temp.s, strerror(errno));            
        }
    }
    temp.l = 0;
    if ( args.output_dir )
        ksprintf(&temp, "%s/failed_1.%s.gz", args.output_dir, file_is_fastq ? "fq" : "fa");
    else
        ksprintf(&temp,"failed_1.%s.gz", file_is_fastq ? "fq" : "fa");

    args.failed_1 = bgzf_open(temp.s, "w");
    args.failed_2 = NULL;
    if ( args.failed_1 == NULL )
        error("%s : %s.", temp.s, strerror(errno));

    if ( seq2 ) {
        temp.l = 0;
        if ( args.output_dir )
            ksprintf(&temp, "%s/failed_2.%s.gz", args.output_dir, file_is_fastq ? "fq" : "fa");
        else
            ksprintf(&temp,"failed_2.%s.gz", file_is_fastq ? "fq" : "fa");

        args.failed_2 = bgzf_open(temp.s, "w");
        if ( args.failed_2 == NULL )
            error("%s : %s.", temp.s, strerror(errno));
    }
    free(temp.s);

    // check barcodes
    int read_flag = args.barcode_region[0] == 1 ? 1 : 2;    
    int length;
    length = args.end - args.start;
    assert(length > 0 );
    // SE
    if ( seq2 == NULL ) {        
        if ( read_flag == 2 )
            error("Inconsistant region specified. %s", args.barcode_region);
        kstring_t string = {0, 0, 0};
        while ( (l1 = kseq_read(seq1)) >= 0 ) {
            string.l = 0;
            if (seq1->qual.l) {
                ksprintf(&string, "@%s\n%s\n+\n%s\n", seq1->name.s, seq1->seq.s, seq1->qual.s);
            } else {
                ksprintf(&string, ">%s\n%s\n", seq1->name.s, seq1->seq.s);
            }
            for ( i = 0; i < args.barcode.n; ++i ) {
                struct name *name = &args.barcode.names[i];
                int mismatch = check_match(seq1->seq.s+args.start,name->barcode, args.mismatch, length);
                if ( mismatch == -1 )
                    continue;
                if ( bgzf_write(name->fp1, string.s, string.l) != string.l)
                    error("Write error : %d", name->fp1->errcode);
            }
            if ( i == args.barcode.n)
                if ( bgzf_write(args.failed_1, string.s, string.l) != string.l)
                    error("Write error : %d", args.failed_1->errcode);            
        }
        if ( string.l )
            free(string.s);
        
    } else { // PE
        kstring_t string1 = {0, 0, 0};
        kstring_t string2 = {0, 0, 0};
        while ( (l1 = kseq_read(seq1)) >=0 && (l2 = kseq_read(seq2)) >= 0) {
            // check the read name
            for ( i = 0; i < seq1->name.l; ++i )
                if ( seq1->name.s[i] != seq2->name.s[i])
                    error("Inconsistant read name. %s vs %s.", seq1->name.s, seq2->name.s);
            
            string1.l = 0;
            string2.l = 0;
            
            if ( file_is_fastq ) {
                ksprintf(&string1, "@%s\n%s\n+\n%s\n", seq1->name.s, seq1->seq.s, seq1->qual.s);
                ksprintf(&string2, "@%s\n%s\n+\n%s\n", seq2->name.s, seq2->seq.s, seq2->qual.s);
            } else {
                ksprintf(&string1, ">%s\n%s\n", seq1->name.s, seq1->seq.s);
                ksprintf(&string2, ">%s\n%s\n", seq2->name.s, seq2->seq.s);
            }
            
            for ( i = 0; i < args.barcode.n; ++i ) {
                struct name *name = &args.barcode.names[i];
                int mis;
                if ( read_flag == 1 )
                    mis = check_match(seq1->seq.s+args.start, name->barcode, args.mismatch, length);
                else
                    mis = check_match(seq2->seq.s+args.start, name->barcode, args.mismatch, length);
                if ( mis == -1 )
                    continue;
                if ( bgzf_write(name->fp1, string1.s, string1.l) != string1.l)
                    error("Write error : %d", name->fp1->errcode);
                if ( bgzf_write(name->fp2, string2.s, string2.l) != string2.l)
                    error("Write error : %d", name->fp2->errcode);
            }
            if ( i == args.barcode.n) {
                if ( bgzf_write(args.failed_1, string1.s, string1.l) != string1.l)
                    error("Write error : %d", args.failed_1->errcode);            
                if ( bgzf_write(args.failed_2, string2.s, string2.l) != string2.l)
                    error("Write error : %d", args.failed_2->errcode);                            
            }

        }
        if (l1 < 0 || l2 < 0 )
            error("Inconsistant read records. %d vs %d",l1, l2);
    }

    for ( i = 0; i < args.barcode.n; ++i ) {
        struct name *name = &args.barcode.names[i];
        free(name->name);
        free(name->barcode);
        bgzf_close(name->fp1);
        if ( name->fp2 )
            bgzf_close(name->fp2);
    }
    bgzf_close(args.failed_1);    
    if ( args.failed_2)
        bgzf_close(args.failed_2);
    
    free(args.barcode.names);

    return 0;
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    if ( split_barcode() )
        return 1;

    return 0;
}
