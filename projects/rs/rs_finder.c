// rs_finder.c - Convert RS id to genome coordinate
// 
#include "utils.h"
#include "number.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include <ctype.h>

KHASH_MAP_INIT_INT64(dat, uint64_t)

typedef kh_dat_t hash_t;
                      
struct args {
    const char * data_fname;
    const char * bad_fname;
    const char * input_fname;
    hash_t *hash;
} args = {
    .data_fname = NULL,
    .bad_fname = NULL,
    .input_fname = NULL,
    .hash = NULL,
};

int usage()
{
    fprintf(stderr, "rs_finder -data RS_id.tsv -bad bad_output.tsv input.tsv\n"
            "The format of RS_id.tsv should be (seperated by tab):\n"
            "  chrom, pos, ref, alt, rs_id\n"
            "The format of input.tsv should be (seperated by tab):\n"
            "  rs_id, extend information ...\n"
        );
    
    return 1;
}

int parse_args(int argc, char **argv)
{
    if ( argc == 1)
        return usage();

    int i;
    for ( i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if ( strcmp(a, "-data") == 0 )
            var = &args.data_fname;
        else if ( strcmp(a, "-bad") == 0 )
            var = &args.bad_fname;

        if ( var != 0 ) {
            if ( argc == i )
                error("Missing argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if ( args.input_fname == 0 ) {
            args.input_fname = a;
            continue;
        }

        error("Unknown argument, %s.", a);                
    }
    return 0;
}

int load_dataset()
{
    FILE *fp = fopen(args.data_fname, "r");
    if ( fp == NULL )
        error("%s : %s", args.data_fname, strerror(errno));

    args.hash = kh_init(dat);
    kstring_t str = KSTRING_INIT;
    uint64_t last_offset = 0;

    // rs_flag set to check '\t', 'r', 's', start parse the downstream characters
    // into number if and only if rs_flag == 3
    int rs_flag = 0;
    int new_line = 1;
    
    for ( ;; ) {
        char c = fgetc(fp);
        if ( feof(fp) )
            break;

        // skip comments
        if ( new_line == 1 && c == '#' ) {
            for ( ; !feof(fp) && c != '\n'; ) c = fgetc(fp);
            last_offset = ftell(fp);
        }

        // unset new line flag
        if ( new_line == '\n' )
            new_line = 1;
        else
            new_line = 0;
        
        if ( c == '\t' || c == 'r' || c == 's' ) {
            rs_flag ++;
            continue;
        }

        // put rs id to hash as a key, offset is value
        if ( c == '\n' && rs_flag == 3 && str.l != 0 ) {
            int id = str2int(str.s);
            int ret;
            khiter_t k = kh_get(dat, args.hash, id);
            if ( k == kh_end(args.hash) )
                k = kh_put(dat, args.hash, id, &ret);
            else
                warnings("Duplicated record rs%d.", id);

            kh_val(args.hash, k) = last_offset;
            last_offset = ftell(fp);
            str.l = 0;

            // set flag to mark new line
            new_line = 1;

            continue;
        }
        
        if ( rs_flag == 3 ) {            

            if ( isdigit((int)c) )
                kputc(c, &str);

            continue;
        }
        
        rs_flag = 0;
    }
    
    fclose(fp);
    if ( str.m )
        free(str.s);
    
    return 0;
}

int parse_line()
{
    
    kstring_t str = KSTRING_INIT;
    int rs_flag = 0;
    int new_line = 1;
    
    FILE *fp = fopen(args.input_fname, "r");
    if ( fp == NULL )
        error("%s : %s.", args.input_fname, strerror(errno));

    FILE *data = fopen(args.data_fname, "r");
    FILE *bad = args.bad_fname == NULL ? stderr : fopen(args.bad_fname, "r");
    if ( bad == NULL )
        error("%s : %s.", args.bad_fname, strerror(errno));
    
    for ( ;; ) {

        char c = fgetc(fp);
        if ( feof(fp) )
            break;

        // skip comments
        if ( new_line == 1 && c == '#' )
            for ( ; !feof(fp) || c != '\n'; ) c = fgetc(fp);

        // reset new line flag
        if ( c == '\n' || c == 'r' || c == 's' )
            new_line = 1;
        else
            new_line = 0;

        if ( c == 'r' || c == 's' ) {
            rs_flag ++;
            continue;
        }

        // all records init with rs_id column
        if ( rs_flag == 2 ) {
            for ( ;; ) {
                if ( !isdigit(c) )
                    break;
                kputc(c, &str);
                c = fgetc(fp);
            }
            int rs_id = str2int(str.s);
            khiter_t k = kh_get(dat, args.hash, rs_id);
            if ( k == kh_end(args.hash) ) {
                fprintf(bad, "rs%d", rs_id);
                for ( ;; ) {
                    fputc(c, bad);
                    if ( c == '\n' || feof(bad) )
                        break;
                    c = fgetc(fp);
                }
            }
            else {
                uint64_t offset = kh_val(args.hash, k);
                fseek(data, offset, SEEK_SET);
                for ( ;; ) {
                    int c = fgetc(data);
                    // emssion '\n'
                    if ( c == '\n' || feof(data) )
                        break;
                    putc(c, stdout);
                }

                for ( ;; ) {
                    putc(c, stdout);
                    if ( c == '\n' || feof(fp) )
                        break;
                    c = fgetc ( fp);
                }
            }
            str.l = 0;
        }

        if ( c == '\n' ) {
            rs_flag = 0;
            continue;
        }
        error("%s : Unrecognised format. %c", args.input_fname, c);
    }
    fclose(fp);
    fclose(data);
    fclose(bad);
    return 0;
}

void memory_release()
{
    if ( args.hash ) {
        kh_destroy(dat, args.hash);
    }
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    if ( load_dataset() )
        error("Failed to load dataset. %s", args.data_fname);

    parse_line();
    memory_release();

    return 0;
}
