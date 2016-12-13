// CNV_regions_format_per_sample.c - format CNV regions
// Some CNV called use splited reads and paired reads to detect the SV for single allele, but they
// failed to consider the diploid nature of the genome to process the allele check. So we may found
// the overlaped region between two nearby heterozygote copy number variantions. This program will
// check the discordances and split the CNV regions into severl pieces and rebuild the genotypes for
// the overlaped regions.

#include "utils.h"
#include "number.h"
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/khash.h>
#include <zlib.h>
#define KSTRING_INIT {0, 0, 0}
KHASH_MAP_INIT_STR(name, int)
typedef kh_name_t namehash_type;
KSTREAM_INIT(gzFile, gzread, 8192)
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
struct bed_cnv {
    struct bed_cnv *next;
    int id;
    int start;
    int end;
    int flag;
};

struct args {
    const char * input_fname;
    const char * output_fname;
    FILE *fp_out;
    int min_length;
    int max_length;
    int report;
    int n_names;
    int m_names;
    char **names;
    namehash_type *hash;
    char *sample_name;
    struct bed_cnv *node;
} args = {
    .input_fname = 0,
    .output_fname = 0,
    .fp_out = NULL,
    .min_length = 0,
    .max_length = 0,
    .report = 0,
    .n_names = 0,
    .m_names = 0,
    .hash = NULL,
    .names = NULL,
    .sample_name = NULL,
    .node = NULL,
};


static int flag_inconsis(int flag)
{
    int i, j;
    int c = 0;
    for (i = 0, j = 0; i < 4; ++i ) {
        c += (flag>>i) & 1;
    }
    return c < 3 ? 0 : 1;
}
int parse_type(char *type)
{
    if ( !strcmp(type, "N/N") || !strcmp(type, "<REF>/<REF>")  || !strcmp(type, "N") || !strcmp(type, "<REF>") )
        return 0x0;
    if ( !strcmp(type, "N/<DEL>") || !strcmp(type, "<DEL>/N") )
        return 0x1;
    if ( !strcmp(type, "N/<DUP>") || !strcmp(type, "<DUP>/N") )
        return 0x4;
    if ( !strcmp(type, "<DUP>/<DUP>"))
        return 0xc;
    if ( !strcmp(type, "<DEL>/<DEL>"))
        return 0x3;
    return 0x0;
}
// accept bed format, Chrom\tStart(1-based)\tEnd\tType\t[Sample]
int parse_line(kstring_t *string, struct bed_cnv *line)
{
    if ( args.hash == NULL )
        args.hash = kh_init(name);
    
    if ( string->l == 0 )
        return 1;
    if (string->s[0] == '#')
        return 1;
    
    char *raw_line = strdup(string->s);
    int nfields = 0;
    int *splits = ksplit(string, 0, &nfields);
    if (splits == NULL )
        return 1;
    if ( nfields < 4)
        goto skip_line;
    
    char *chrom = string->s+splits[0];
    char *start = string->s+splits[1];
    char *end = string->s+splits[2];
    char *type = string->s+splits[3];
    char *sample = nfields > 4 ? string->s+splits[4] : NULL;

    // Only check DUP/DEL regions, skip other types
    if (*start == '.' || *end == '.')
        goto skip_line;

    line->next = NULL;
    khint_t k;
    k = kh_get(name, args.hash, chrom);
    if ( k == kh_end(args.hash) ) {
        if ( args.n_names == args.m_names) {
            args.m_names += 4;
            args.names = (char**)realloc(args.names, args.m_names*sizeof(char*));
        }
        int ret;
        args.names[args.n_names] = strdup(chrom);
        k = kh_put(name, args.hash, args.names[args.n_names], &ret);
        kh_val(args.hash, k) = args.n_names;
        args.n_names++;
    }
    // init chromosome
    line->id = kh_val(args.hash, k);
    
    if ( check_num_likely(start) ) {
        // Convert the 1 based start position to 0 based
        line->start = str2int(start) -1;
    } else {
        fprintf(stderr, "Failed to parse line %s.\n", raw_line);
        goto skip_line;
    }

    if ( check_num_likely(end) ) {
        line->end = str2int(end);
    } else {
        fprintf(stderr, "Failed to parse line %s.\n", raw_line);
        goto skip_line;
    }

    if ( args.min_length && line->end - line->start < args.min_length)
        goto skip_line;
    if ( args.max_length && line->end - line->start > args.max_length)
        goto skip_line;
    line->flag = parse_type(type);
    if ( check_malformed(line->flag) || flag_inconsis(line->flag) ) {
        fprintf(stderr,"Failed to parse type %s, %d\n", type, line->flag);
        goto skip_line;
    }

    // If sample column is set, check the sample names then.
    if ( sample ) {
        if ( args.sample_name == NULL) {
            args.sample_name = strdup(sample);
        } else if ( strcmp(args.sample_name, sample) ) {
            error("Inconsistance sample names. %s vs. %s.", args.sample_name, sample);
        }
    }        
    free(raw_line);
    return 0;
    
  skip_line:    
    free(raw_line);
    free(splits);
    return 1;
}
int combine_flag(int flag1, int flag2)
{
    if ( (flag1 & 0x3) && (flag2 & 0x3) )
        return 0x3;
    if ( (flag1 & 0x4) && (flag2 & 0x4) )
        return 0xc;
    return flag1 | flag2;            
}
const char * const types[] = {
    "<REF>/<REF>",
    "<REF>/<DEL>",
    "<DEL>/<DEL>",
    "<REF>/<DUP>",
    "<DUP>/<DUP>",
    "<DEL>/<DUP>",
    "<UNKNOWN>",
};
const char *explain_type(int flag)
{
    flag &= 0xf;
    if ( !flag )
        return types[0];
    if ( flag & 0x1 ) {
        if ( flag & 0x2 )
            return types[2];                
        if ( flag & 0xc )
            return types[5];
        return types[1];
    } 
    if ( flag & 0x4) {
        if (flag & 0x8)
            return types[4];
        return types[3];
    }
    return types[6];
}
int print_node(struct bed_cnv *node)
{
    if (node == NULL)
        return 1;
    fprintf(args.fp_out,"%s\t%d\t%d\t%s\n", args.names[node->id], node->start, node->end, explain_type(node->flag));
    return 0;
}
int push_node(kstring_t *string)
{
    struct bed_cnv *node = (struct bed_cnv*)malloc(sizeof(struct bed_cnv));
    struct bed_cnv *temp = args.node;
    if ( parse_line(string, node) ) {
        free(node);
        return 0;
    }
        
    if (args.node == NULL) {
        goto update_line;
    }
    
    if ( node->id != temp->id || temp->end <= node->start) {
        print_node(temp);
        free(temp);
        goto update_line;
    }
    
    if (node->start < temp->start){                
        fprintf(stderr, "Regions are not properly sorted. %s : %d %d; %d %d\n", args.names[node->id], node->start+1, node->end, temp->start, temp->end);
        // int end = node->end;
        // node->end = temp->start;
        // node->flag = 0xf;
        // print_node(node);
        node->start = temp->start;
        // node->end = end;
        if ( node->end < temp->end) {
            node->flag = combine_flag(node->flag, temp->flag);
            print_node(node);
            temp->start = node->end;
            free(node);
            node = temp;
        } else if ( node->end == temp->end) {
            node->flag = combine_flag(node->flag, temp->flag);
            print_node(node);
            free(node);
            free(temp);
            node = NULL;
        } else {
            temp->flag = combine_flag(node->flag, temp->flag);
            print_node(temp);
            node->start = temp->end;
            free(temp);
        }
        // return 0;
    } else if ( node->start == temp->start) {            
        if ( node->end < temp->end ) {
            fprintf(stderr, "Regions are not properly sorted. %s : %d %d; %d %d\n", args.names[node->id], node->start+1, node->end, temp->start, temp->end);
            node->flag = combine_flag(node->flag, temp->flag);
            print_node(node);
            temp->start = node->end;
            free(node);
            node = temp;
            // return 1;
        } else if ( node->end == temp->end )  {                
            temp->flag = combine_flag(temp->flag, node->flag);
            print_node(temp);
            free(temp);
            free(node);
            node = NULL;
        } else {
            temp->flag = combine_flag(temp->flag, node->flag);
            print_node(temp);
            node->start = temp->end;
            free(temp);
        }
    } else {
        int end = temp->end;
        temp->end = node->start;
        print_node(temp);
        temp->end = end;
        temp->start = node->start;
        if ( temp->end > node->end) {
            node->flag = combine_flag(temp->flag, node->flag);
            print_node(node);
            temp->start = node->end;            
            free(node);
            node = temp;
        } else if ( temp->end == node->end) {
            node->flag = combine_flag(temp->flag, node->flag);
            print_node(node);
            free(temp);
            free(node);
            node = NULL;
        } else {
            temp->flag = combine_flag(temp->flag, node->flag);
            print_node(temp);
            node->start = temp->end;
            free(temp);                    
        }
    }
  update_line:
    args.node = node;
    return 0;
}
int usage(char *name)
{
    fprintf(stderr,
            "%s [options] input_cnv.bed\n"
            "-o <file>         output.bed\n"
            "-report           with report.\n"
            "-min <length>     minimal length to skip.\n"
            "-max <length>     maximal length capped to.\n"
            "\nHomepage:\n"
            "https://github.com/shiquan/small_projects_collections/tree/master/projects/CNV_database"
            , name); 
    return 1;
}

int parse_args(int ac, char **av)
{
    if ( ac == 0 )
        return usage(av[0]);
    
    int i;
    const char *minimal = NULL;
    const char *maximal = NULL;
    for ( i = 1; i < ac; ) {
        const char *a = av[i++];
        const char **var =0;
        if ( strcmp(a, "-o") == 0 && args.output_fname == 0 )
            var = &args.output_fname;
        else if ( strcmp(a, "-report") == 0)
            args.report = 1;
        else if ( strcmp(a, "-min") == 0 && minimal == 0 )
            var = &minimal;
        else if ( strcmp(a, "-max") == 0 && maximal == 0 )
            var = &maximal;
        if ( var != 0 ) {
            if ( i == ac ) 
                error("Missing an argument after %s.",a);
            *var = av[i++];
            continue;
        }
        if ( args.input_fname == 0 ) {
            args.input_fname = a;
            continue;
        }        
        error("Unknown argument %s", a);
    }
    if ( args.input_fname == NULL) {
        if ( !isatty(fileno(stdin)) )
            args.input_fname = "-";
        else 
            return usage(av[0]);
    }

    if ( minimal && check_num_likely(minimal) )
        args.min_length = str2int((char*)minimal);
    if ( maximal && check_num_likely(maximal) )
        args.max_length = str2int((char*)maximal);
    if ( args.output_fname ) {
        args.fp_out = fopen(args.output_fname, "w");
        if ( args.fp_out == NULL) {
            fprintf(stderr, "%s : %s\n", args.output_fname, strerror(errno));
            args.fp_out = stdout;
        }
    } else {
        args.fp_out = stdout;
    }
    return 0;
}
int rebuild_regions()
{
    gzFile fp = strcmp(args.input_fname, "-") == 0 ? gzdopen(fileno(stdin), "r") : gzopen(args.input_fname, "r");
    if ( fp == NULL )
        error("%s : %s\n", args.input_fname, strerror(errno));

    kstream_t *ks;
    int dret = 0;
    ks = ks_init(fp);
    khiter_t k;
    kstring_t string = KSTRING_INIT;

    while( ks_getuntil(ks, 2, &string, &dret) >= 0 ) {
        if ( string.l == 0 )
            continue;
        if ( string.s[0] == '#')
            continue;
        if ( push_node(&string) )
            return 1;
        string.l = 0;
    }
    if ( args.node )
        print_node(args.node);
    
    if ( string.m )
        free(string.s);
    gzclose(fp);
    return 0;
}
int export_report()
{
    if ( args.report == 0)
        return 1;

    return 0;
}
int release_memory()
{
    fclose(args.fp_out);
    if ( args.hash )
        kh_destroy(name, args.hash);
    int i;
    for ( i = 0; i < args.n_names; ++i )
        free(args.names[i]);
    free(args.names);
    return 0;
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    if ( rebuild_regions() )
        return 1;
    
    export_report();
    
    release_memory();
    return 0;
}
