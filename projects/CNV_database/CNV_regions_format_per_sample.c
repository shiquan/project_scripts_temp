// CNV_regions_format_per_sample.c - format CNV regions
// Some CNV called use splited reads and paired reads to detect the SV for single allele, but they
// failed to consider the diploid nature of the genome to process the allele check. So we may found
// the overlaped region between two nearby heterozygote copy number variantions. This program will
// check the discordances and split the CNV regions into severl pieces and rebuild the genotypes for
// the overlaped regions.

#include "utils.h"
#include "number.h"
#include "cnv_bed.h"
#include <htslib/kstring.h>

struct args {
    const char * input_fname;
    const char * output_fname;
    FILE *fp_out;
    int report;
    struct cnv_bed *node;
    struct cnv_spec *spec;
} args = {
    .input_fname = 0,
    .output_fname = 0,
    .fp_out = NULL,
    .report = 0,
    .node = NULL,
    .spec = NULL,
};

int print_node(struct cnv_bed *node)
{
    if (node == NULL)
        return 1;
    fprintf(args.fp_out,"%s\t%d\t%d\t%s", args.spec->chrom[node->id], node->start, node->end, explain_type(node->flag));
    if ( args.spec->n_samples == 1 )
        fprintf(args.fp_out, "\t%s", args.spec->samples[0]);
    fputc('\n', args.fp_out);
    return 0;
}
int push_node(struct cnv_bed *node)
{
    struct cnv_bed *temp = args.node;
    struct cnv_bed **pp = &args.node;
    if (args.node == NULL) {
        goto update_line;
    }   
    // LOG_print("temp: %s\t%d\t%d",args.spec->chrom[temp->id], temp->start+1, temp->end);
    // LOG_print("node: %s\t%d\t%d",args.spec->chrom[node->id], node->start+1, node->end);

    if ( node->id != temp->id || temp->end <= node->start) {
        print_node(temp);
        free(temp);
        goto update_line;
    }
    
    if ( node->start < temp->start ){                
        fprintf(stderr, "Regions are not properly sorted. %s : %d %d; %d %d\n", args.spec->chrom[node->id], node->start+1, node->end, temp->start, temp->end);
        node->start = temp->start;
        goto check_end;
    } else if ( node->start == temp->start) {            
        //if ( node->end < temp->end ) {
        //    fprintf(stderr, "Regions are not properly sorted. %s : %d %d; %d %d\n", args.spec->chrom[node->id], node->start+1, node->end, temp->start, temp->end);
        //}
        goto check_end;
    } else {
        int end = temp->end;
        temp->end = node->start;
        print_node(temp);
        temp->end = end;
        temp->start = node->start;
        goto check_end;
    }
    
  check_end:
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

  update_line:
    *pp = node;
    return 0;
}
int usage(char *name)
{
    fprintf(stderr,
            "* Format CNV outputs by split all overlap and overhang regions.\n"
            "Usage:"
            "%s [options] input_cnv.bed\n"
            "-o <file>         output.bed\n"
            "-report           with report.\n"
            "-min <length>     minimal length to skip.\n"
            "-max <length>     maximal length capped to.\n"
            "\nHomepage:\n"
            "https://github.com/shiquan/small_projects_collections/tree/master/projects/CNV_database\n"
            , name); 
    return 1;
}

int parse_args(int ac, char **av)
{
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
        else if ( strcmp(a, "-h") == 0 )
            return usage(av[0]);    

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
    args.spec = cnv_spec_init();
    
    if ( minimal && check_num_likely(minimal) )
        args.spec->min_length = str2int((char*)minimal);
    if ( maximal && check_num_likely(maximal) )
        args.spec->max_length = str2int((char*)maximal);
    cnv_load_fname(args.spec, args.input_fname);
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
    for (;;) {
        struct cnv_bed *node = (struct cnv_bed*)malloc(sizeof(struct cnv_bed));
        if ( cnv_read(args.spec, node) == -1 )
            break;
        // convert 1 based position to 0 based start
        node->start--; 
        push_node(node);
    }
    if ( args.node )
        print_node(args.node);
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
    cnv_spec_destroy(args.spec);
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
