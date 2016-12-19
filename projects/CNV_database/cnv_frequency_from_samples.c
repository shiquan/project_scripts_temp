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
    int end;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .sample_list = NULL,
    .fp_out = NULL,
    .node = NULL,
    .end = 0,
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
    LOG_print("%d", ac);
    for (i = 1; i < ac; ) {
        const char *a = av[i++];
        const char **var = 0;
        if ( strcmp(a, "-o") == 0 && args.output_fname == 0)
            var = &args.output_fname;
        else if ( strcmp(a, "-min") == 0 && minimal == 0 )
            var = &minimal;
        else if ( strcmp(a, "-max") == 0 && maximal == 0 )
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
    LOG_print("input : %s", args.input_fname);
    args.spec = cnv_spec_init();
    cnv_load_fname(args.spec, args.input_fname);
    LOG_print("load success.");
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
    LOG_print("Parse arguments success.");
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
        } else if ( last_id != node->id || last_start != node->start || last_end != node->end) {
            fprintf(args.fp_out, "%s\t%d\t%d\t<DUP>,<DEL>\t%d,%d\n",
                    args.spec->chrom[last_id], last_start, last_end, dup_count, del_count);
            last_id = node->id;
            last_start = node->start;
            last_end = node->end;
            del_count = 0;
            dup_count = 0;
        }
        if ( node->flag & CNV_DEL_HET ) {                
            del_count ++;
        }
        if ( node->flag & CNV_DEL_HOM ) {
            del_count+=2;
        }
        if ( node->flag & CNV_DUP_HET ) {
            dup_count++;
        }
        if ( node->flag & CNV_DUP_HOM ) {
            dup_count+=2;
        }
    }
    if ( last_id != -1) {
        fprintf(args.fp_out, "%s\t%d\t%d\t<DUP>,<DEL>\t%d,%d\n",
                args.spec->chrom[last_id], last_start, last_end, dup_count, del_count);
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
            splitors.n++;
            break;
        }        
    }
    if ( i == splitors.n )
        splitors.splitors[splitors.n++] = pos;
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
    struct cnv_bed **bed1 = (struct cnv_bed **)a;
    struct cnv_bed **bed2 = (struct cnv_bed **)b;
    // if ( bed1->start == bed2->start )
    // return bed2->end - bed1->end;
    return (*bed1)->start - (*bed2)->start;
}
struct cnv_bed *split_blocks(struct cnv_bed *node)
{
    if ( splitors.n == 0 )
        return node->next;
    int i;
    for (i = 0; i < splitors.n; i++) {
        int splitor = splitors.splitors[i];

        if ( splitor >= node->end)
            break;

        if ( splitor <= node->start)
            continue;

        struct cnv_bed *temp = malloc(sizeof(struct cnv_bed));
        temp->id = node->id;
        temp->flag = node->flag;
        temp->start = splitor;
        temp->end = node->end;
        temp->next = node->next;
        node->end = splitor;
        node->next = temp;
        node = temp;
    }
    return node->next;
}
void split_and_count_allele()
{
    struct cnv_bed * temp = args.node;
    struct cnv_bed **pp = &args.node;
    if ( temp == NULL)
        return;
    // clean splitors buffer
    splitors.n = 0;
    // generate the splitors    
    while ( temp )  {
        // if ( temp->start == pos )
        // break;
        // if ( pos > 0 && temp->start > pos )
        //    error("The cached list is not properly sorted. %s %d, %d.",  args.spec->chrom[temp->id], temp->start, pos);
        push_splitor(temp->start);
        push_splitor(temp->end);
        temp = temp->next;
    }
    //if ( pos > 0)
    //   set_cutoff_splitor(pos);
    int i;
    
    // generate blocks
    // int last_pos;
    temp = *pp;
    for ( ;; ) {
        // last_pos = -1;
        temp = split_blocks(temp);
        if ( temp == NULL)
            break;
    }    
    // splitors.n = 0;


    sort_list(pp, bed_comp);
    print_node(*pp, NULL);

    for ( ; *pp; ) {
        struct cnv_bed *node = *pp;
        *pp = (*pp)->next;
        free(node);
    }
    args.end = 0;
    // return *pp;
}
int finish_node()
{
    split_and_count_allele();
    return 0;
}
void make_header(struct cnv_bed *node)
{
    if (args.node != NULL)
        error("Try to realloc a new header for a nonfree header point.");
    args.node = node;
    args.end = node->end;
}
static int push_node(struct cnv_bed *node)
{
    struct cnv_bed *temp= args.node;
    struct cnv_bed **pp = &args.node;    
    // output format "chrom\tstart\tend\t<DUP>,<DEL>\tDEL allele count,Dup allele count\n"

    if ( temp== NULL || node->id != temp->id ) {
        finish_node();
        make_header(node);
        return 0;
    }
    for ( ;; ) {       
        // The node->start come from any upstream position than cached nodes, suggest the input regions are not porperly sorted.
        if (node->start < temp->start) 
            error("Position not properly sorted. %s: %d vs %d.", args.spec->chrom[node->id], node->start, temp->start);
        
        // the new node will be added at the tail of the list.
        if ( temp->next == NULL) {            

            if ( args.end && args.end <= node->start) {
                split_and_count_allele();
                make_header(node);                
            } else {
                temp->next = node;
                if ( args.end < node->end)
                    args.end = node->end;
            }
            break;
        }
        temp = temp->next;
    }
    // for ( temp = *pp; temp && temp != node; temp = temp->next ) {
        // Assume the cached regions are sorted, if some regions are located in the upstream of this node, which means the end positions
        // of these are smaller than node->start, all the regions should be printed out. And the next region closed with the printed regions
        // will be updated to the header of the list.
    // if ( temp->start >= node->start) {
    // *pp = split_and_count_allele(node->start);
            // break;
// }
//    }
    
    return 0;
}
int freq_calculate()
{
    LOG_print("Start calculate frequencies.");
    while ( 1 ) {
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

int main(int argc, char ** argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    freq_calculate();

    release_memory();
    
    return 0;
}
