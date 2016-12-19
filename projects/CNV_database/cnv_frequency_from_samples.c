#include "utils.h"
#include "cnv_bed.h"
#include "number.h"
#include "sort_list.h"
#include <errno.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/khash.h>
#include <zlib.h>

struct args {
    const char * input_fname;
    const char * output_fname;
    const char * sample_list;
    FILE *fp_out;    
    struct cnv_bed *node;
    struct cnv_spec *spec;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .sample_list = NULL,
    .fp_out = NULL,
    .node = NULL,
};

int usage(const char *name)
{
    fprintf(stderr, "%s -o output.tsv cnv_merged.bed [-min <min_length> -max <max_length>]\n", name);    
    return 1; 
}
int parse_args(int ac, char **av)
{
    if ( ac == 1)
        return usage(av[0]);
    int i;
    const char *minimal = 0;
    const char *maximal = 0;
    for (i = 1; i < ac; ) {
        const char *a = av[i++];
        const char **var = 0;
        if ( strcmp(a, "-o") == 0 && args.output_fname == 0)
            var = &args.output_fname;
        else if ( strcmp(minimal, "-min") == 0 && minimal == 0 )
            var = &minimal;
        else if ( strcmp(maximal, "-max") == 0 && maximal == 0 )
            var = &maximal;
        
        if ( var != 0 ) {
            if ( i == ac )
                error("Missing an argument after %s.", a);
            *var = av[i++];
            continue;
        }
        if (args.input_fname == 0) {
            args.input_fname = a;
            continue;
        }
        error("Unknown argument %s.", a);
    }

    if ( args.input_fname == NULL ) {
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

    if ( args.output_fname ) {
        args.fp_out = fopen(args.output_fname, "w");
        if ( args.fp_out == NULL ) {
            fprintf(stderr, "%s : %s\n", args.output_fname, strerror(errno));
            args.fp_out = stdout;
        }
    } else {
        args.fp_out = stdout;
    }
    return 0;    
}
static void print_node(struct cnv_bed *node, struct cnv_bed *end)
{
    int last_id = -1;
    int last_start = -1;
    int last_end = -1;
    int del_count = 0;
    int dup_count = 0;
    for ( ; node && node != end; node = node->next ) {
        if ( last_id == -1 ) {
            last_id = node->id;
            last_start = node->start;
            last_end = node->end;
        } else if ( last_id == node->id && last_start == node->start && last_end == node->end) {
            if ( node->flag & CNV_DEL_HET ) {                
                 del_count ++;
                 if ( node->flag & CNV_DEL_HOM ) {
                     del_count++;
                 } else if ( node->flag & CNV_DUP_HET ) {
                     dup_count ++;
                 }
            } else if ( node->flag & CNV_DUP_HET ) {
                dup_count++;
                if ( node->flag & CNV_DUP_HOM ) {
                    dup_count++;
                }
            }
        } else {
            fprintf(args.fp_out, "%s\t%d\t%d<DUP>,<DEL>\t%d,%d\n",
                    args.spec->chrom[last_id], last_start, last_end, dup_count, del_count);
            last_id = node->id;
            last_start = node->start;
            last_end = node->end;
        }
    }
}

struct cnv_splitors {
    int n;
    int m;
    int *splitors;
} splitors = {
    .n = 0,
    .m = 0,
    .splitors = 0,
};

void push_splitor(int pos)
{
    int i;
    if ( splitors.n == splitors.m ) {
        splitors.m += 10;
        splitors.splitors = realloc( splitors.splitors, splitors.m * sizeof(int));
    }
    for ( i = 0; i < splitors.n; ++i ) {
        if ( splitors.splitors[i] == pos ) {
            break;
        } else if ( splitors.splitors[i] > pos) {
            memmove(splitors.splitors + i+1, splitors.splitors + i, (splitors.n - i) * sizeof(int));
            splitors.splitors[i] = pos;
            break;
        }        
    }
}
void set_cutoff_splitor(int pos)
{
    int i;
    for ( i = 0; i < splitors.n; ++i ) {
        if ( splitors.splitors[i] >= pos ) {
            splitors.splitors[i] = pos;
            splitors.n = i+1;
            break;
        }
    }
    if ( i == splitors.n )
        push_splitor(pos);
}
int bed_comp(const void *a, const void *b)
{
    struct cnv_bed *bed1 = (struct cnv_bed *)a;
    struct cnv_bed *bed2 = (struct cnv_bed *)b;
    if ( bed1->start == bed2->start )
        return bed1->end - bed2->end;
    return bed1->start - bed2->start;
}
struct cnv_bed *split_and_count_allele(int pos)
{
    struct cnv_bed * temp = args.node;
    struct cnv_bed **pp = &args.node;

    // generate the splitors    
    while ( 1 )  {
        if ( temp->start == pos )
            break;
        if ( temp->start > pos )
            error("The cached list is not properly sorted.");
        push_splitor(temp->start);
        push_splitor(temp->end);
        temp = temp->next;
    }
    if ( pos > 0)
        set_cutoff_splitor(pos);

    // generate blocks
    int i;
    // int last_pos;
    for ( temp = *pp; temp && temp->start < pos; ) {
        // last_pos = -1;
        for (i = 0; i < splitors.n; ++i ) {
            if ( splitors.splitors[i] > temp->start && splitors.splitors[i] <= temp->end) {
                struct cnv_bed * node = malloc(sizeof(struct cnv_bed));
                node->id = temp->id;
                node->start = splitors.splitors[i];
                node->end = temp->start;
                node->next = temp->next;
                temp->next = node;
                temp->end = node->start;
                temp = node;
            } else {
                break;
            }
        }
        temp = temp->next;
    }    
    splitors.n = 0;
    sort_list(pp, bed_comp);
    for ( temp = *pp; temp && temp->start < pos; temp = temp->next);
    print_node(*pp, temp);
    for ( ; *pp && *pp != temp; ) {
        struct cnv_bed *node = *pp;
        *pp = (*pp)->next;
        free(node);
    }
    return *pp;
}
int finish_node()
{
    split_and_count_allele(0);
    return 0;
}
static int push_node(struct cnv_bed *node)
{
    struct cnv_bed *temp= args.node;
    struct cnv_bed **pp = &args.node;
    // output format "chrom\tstart\tend\t<DUP>,<DEL>\tDEL allele count,Dup allele count\n"
    if ( node->id != temp->id ) {
        finish_node();
        *pp = node;
        return 0;
    }
    while ( 1 ) {
        // The node->start come from any upstream position than cached nodes, suggest the input regions are not porperly sorted.
        if (node->start < temp->start) 
            error("Position not properly sorted. %s: %d vs %d.", args.spec->chrom[node->id], node->start, temp->start);
        
        // Assume the cached regions are sorted, if some regions are located in the upstream of this node, which means the end positions
        // of these are smaller than node->start, all the regions should be printed out. And the next region closed with the printed regions
        // will be updated to the header of the list.

        // the new node will be added at the tail of the list.
        if ( temp->next == NULL) {
            temp->next = node;
            break;
        }
        temp = temp->next;
    }
    *pp = split_and_count_allele(node->start);
    
    return 0;
}
int freq_calculate()
{    
    for (;;) {
        struct cnv_bed *node = malloc(sizeof(struct cnv_bed));
        if ( cnv_read(args.spec, node) == -1)
            break;
        push_node(node);
    }
    if ( args.node )
        print_node(args.node, NULL);
    return 0;
}

int release_memory()
{
    fclose(args.fp_out);
    cnv_spec_destroy(args.spec);
    return 0;
}

int main ( int argc, char ** argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    freq_calculate();

    release_memory();
    
    return 0;
}
