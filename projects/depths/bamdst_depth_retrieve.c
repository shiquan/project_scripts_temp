#include "utils.h"
#include "number.h"
#include "htslib/kstring.h"
#include "htslib/tbx.h"

int usage()
{
    fprintf(stderr,
            "bamdst_depth_retrieve [options] depth.tsv.gz\n"
            " -cutoff depth1,depth2     Depth cutoff values to stat coverage, seperated with \",\".\n"
            " -out output.tsv           Output average,median,[coverage].\n"
            " -sum summary.txt          Summary file, export all bases, coverages.\n"
            //" -select <raw|rmdup|cov>   Select depth column in depth.tsv.gz.\n"
            " -col <INT>                Select column to calculate depth.\n"
            " -reg target.bed           Target region in BED format.\n"
        );
    return 1;
}

enum depth_col {
    raw_col = 1,
    rmdup_col,
    coverage_col,
};

struct bed {
    char *chrom;
    int start;
    int end;
    uint64_t total;
    uint64_t uncover;
};

struct args {

    const char    *input_fname;
    const char    *output_fname;
    const char    *summary_fname;
    const char    *data_fname;
    FILE          *fp_input;
    FILE          *fp_output;
    FILE          *fp_summary;
    htsFile       *fp_data;
    
    // enum depth_col col;
    int            col;
    tbx_t         *idx;
    int            n_depth;
    // cache for cutoff values
    int           *depths;

    //uint64_t      *cov_bases;
    uint64_t       total_base;
    uint64_t       total_length;
    uint64_t       n_lines;

    // count depths for each cutoff
    uint64_t      *depths_cutoff;
    uint64_t      *depths_cutoff_per_reg;
    
    // buffers
    struct bed     bed;        
} args = {
    .input_fname   = NULL,
    .output_fname  = NULL,
    .summary_fname = NULL,
    .fp_input      = NULL,
    .fp_output     = NULL,
    .fp_summary    = NULL,
    .col           = 3,
    .idx           = NULL,
    .n_depth       = 0,
    .depths        = NULL,
    //.cov_bases     = NULL,
    .total_base    = 0,
    .total_length  = 0,
    .n_lines       = 0,

    .depths_cutoff = NULL,
    .depths_cutoff_per_reg = NULL,
    .bed          = { 0, 0, 0, 0, 0},
};
int comp(const void *a, const void *b)
{
    int l = *(const int*)a;
    int r = *(const int*)b;
    return l - r;
}
static int *str2intArray(const char *_s, int *n_arr)
{
    char *ss = (char*)_s;
    char *se = (char*)_s;
    int i;
    int l;
    int n = 0;
    int m = 0;
    int *d;

    if ( ss == NULL ) {
        *n_arr = 1;
        d = (int*)malloc(sizeof(int));
        d[0] = 0;
        return d;
    }
    
    l = strlen(ss);
    int itr = -1;
    int start = -1;
    for ( i = 0; i < l; ++i) {        
        if ( ss[i] == ',') {
            // emit first comma
            if ( itr == -1)
                error("Unrecognize format %s", _s);
            
            int num = str2int_l(ss+itr, i - itr);
            
            if ( num < 0 )
                error("Unrecognize format %s", _s);

            if ( start >= 0 ) {
                if ( num < start )
                    error("Unsorted range, %d-%d.", start, num);
                int range = num - start + 1;
                if ( m < n + range) {
                    m = n+ range;
                    d = (int*)realloc(d, m*sizeof(int));
                }
                int j;
                for ( j = start; j <= num; ++j) {                    
                    d[n] = j;
                    n++;
                }
                start = -1;
            }
            else {
                if (m == n ) {
                    m+=2;
                    d = (int*)realloc(d, m*sizeof(int));
                }
                d[n++] = num;
            }
            itr = i+1;
        }
        else if ( ss[i] == '-' ) {
            if ( start >= 0 )
                error("Unrecognize format %s", _s);

            start = str2int_l(ss+itr, i-itr);
            if ( start < 0 )
                error("Unrecognize format %s", _s);

            itr = i + 1;
        }
        else if ( i == l -1 ) {
            int num = str2int_l(ss+itr, i - itr+1);
            
            if ( num < 0 )
                error("Unrecognize format %s", _s);

            if ( start >= 0 ) {
                if ( num < start )
                    error("Unsorted range, %d-%d.", start, num);
                int range = num - start + 1;
                if ( m < n + range) {
                    m = n+ range;
                    d = (int*)realloc(d, m*sizeof(int));
                }
                int j;
                for ( j = start; j <= num; ++j) {                    
                    d[n] = j;
                    n++;
                }
                start = -1;
            }
            else {
                if (m == n ) {
                    m+=2;
                    d = (int*)realloc(d, m*sizeof(int));
                }
                d[n++] = num;
            }         
        }
        else {
            if ( itr == -1 )
                itr = i;
        }
    }

    qsort((void*)d, n, sizeof(int), comp);

    *n_arr = n;
    
    if ( d[0] != 0 ) {
        d = realloc(d, (n+1)*sizeof(int));
        memmove(d, d+1, n*sizeof(int));
        d[0] = 0;
    }
    
    return d;
}

int parse_args(int argc, char **argv)
{
    if ( argc == 1 )
        return usage();

    int i;
    const char *cutoff_string = 0;
    const char *col_str = 0;
    for ( i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;

        if ( strcmp(a, "-reg") == 0 )
            var = &args.input_fname;
        else if ( strcmp(a, "-sum") == 0 )
            var = &args.summary_fname;
        else if ( strcmp(a, "-out") == 0 )
            var = &args.output_fname;
        else if ( strcmp(a, "-cutoff") == 0 )
            var = &cutoff_string;
        else if ( strcmp(a, "-col") == 0 )
            var = &col_str;

        if ( var != 0 ) {
            if ( argc == i )
                error("Missing an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if ( args.data_fname == 0 ) {
            args.data_fname = a;
            continue;
        }
            
        error("Unknown argument, %s.", a);
    }

    if ( args.input_fname == NULL )
        error("Region must be set with -reg.");

    if ( args.data_fname == NULL )
        error("depth.tsv.gz must be specified.");
    
    args.idx = tbx_index_load(args.data_fname);
    if ( args.idx == NULL )
        error("Failed to load index of %s.", args.data_fname);
    
    args.fp_data = hts_open(args.data_fname, "r");
    if ( args.fp_data == NULL )
        error("%s : %s.", args.data_fname, strerror(errno));

    args.fp_input = fopen(args.input_fname, "r");
    if ( args.fp_input == NULL )
        error("%s : %s.", args.input_fname, strerror(errno));

    if ( args.output_fname ) {
        args.fp_output = fopen(args.output_fname, "w");
        if ( args.fp_output == NULL )
            error("%s : %s.", args.output_fname, strerror(errno));
    }

    args.fp_summary = args.summary_fname == NULL ? stdout : fopen(args.summary_fname, "w");
    if ( args.fp_summary == NULL )
        error("%s : %s.", args.summary_fname, strerror(errno));

    if ( col_str )
        args.col = str2int(col_str);

    args.depths = str2intArray(cutoff_string, &args.n_depth);
    args.depths_cutoff = (uint64_t*)calloc(args.n_depth, sizeof(uint64_t));
    args.depths_cutoff_per_reg = (uint64_t*)calloc(args.n_depth, sizeof(uint64_t));
    memset(args.depths_cutoff, 0, args.n_depth*sizeof(uint64_t));
    memset(args.depths_cutoff_per_reg, 0, args.n_depth*sizeof(uint64_t));

    if ( args.fp_output ) {
        fprintf(args.fp_output, "#Chrom\tstart\tend\tave");
        for ( i = 0; i < args.n_depth; ++i )
            fprintf(args.fp_output, "\tcov %dx", args.depths[i]);
        fputc('\n', args.fp_output);
    }
    return 0;
}

static int parse_depthData(char *str, int l)
{
    int i;
    int col = 1;
    int start = 0;

    for ( i = 0; i < l; ++i) {
        if ( col == args.col ) {
            if ( str[i] == '\t' || i == l -1 ) {
                return str2int_l(str+start, i-start);
            }
        }
        else {
            if ( str[i] == '\t' ) {
                col++;
                start = i+1;
            }
        }
    }
    if ( col < args.col)
        error("Unsufficient column in line %s.", str);
    
    return 0;
}
void clean_bed(struct bed *bed)
{
    if ( bed->chrom)
        free(bed->chrom);
    memset(bed, 0, sizeof(struct bed));
}

int read_bed()
{
    if ( feof(args.fp_input ))
        return 1;

    struct bed *bed = &args.bed;
    
    clean_bed(bed);
    
    kstring_t str = { 0, 0, 0};
    int       col = 0;
    
    for ( ;; ) {
        // end of line
        if ( feof(args.fp_input) )
            return 1;
        
        char c = fgetc(args.fp_input);

        // emit comments
        if ( c == '#' ) {
            for ( ; c != '\n' && !feof(args.fp_input); c = fgetc(args.fp_input) );
            continue;
        }
        
#define BRANCH(_key, _func, _len) do {                                       \
            if ( c == '\t' || c == '\n' || feof(args.fp_input) ) {   \
                if ( _len == 0 ) break;\
                _key = _func(str.s, _len);                                   \
                str.l = 0;\
                col++;\
            }\
            else {\
                kputc(c, &str);\
            }\
        } while(0)\
            
        if ( col == 0 ) {
            BRANCH(bed->chrom, strndup, str.l);
        }
        else if ( col == 1 ) {
            BRANCH(bed->start, str2int_l, str.l);
        }
        else if ( col == 2 ) {
            BRANCH(bed->end, str2int_l, str.l);
        }
        else { 
            for ( ; c != '\n' && !feof(args.fp_input); )
                c = fgetc(args.fp_input);
            str.l = 0;                
            break;
        }

#undef BRANCH
    }
    // if truncated return -1
    if ( bed->end == 0 ) {
        bed->end = bed->start;
        bed->start--;
        if ( bed->start < 0 ) {
            warnings("Truncated line. %lld ", args.n_lines);
            return -1;
        }
    }

    args.n_lines ++;
    // clean memory
    if ( str.m ) free(str.s);
    
    return 0;
}

void memory_release()
{
    clean_bed(&args.bed);
    fclose ( args.fp_input );
    if ( args.fp_output )
        fclose ( args.fp_output );
    if ( args.fp_summary )
        fclose ( args.fp_summary);
    if ( args.n_depth )
        free ( args.depths );
    hts_close( args.fp_data );
    
    free ( args.depths_cutoff );
    free ( args.depths_cutoff_per_reg );
    // free ( args.cov_bases );
}
void depths_retrieve()
{
    int id;
    kstring_t str = {0, 0, 0};
    kstring_t reg = {0, 0, 0};
    struct bed *bed = &args.bed;
    for ( ;; ) {
        
        if ( read_bed() )
            break;

        memset(args.depths_cutoff_per_reg, 0, args.n_depth*sizeof(uint64_t));
        
        id = tbx_name2id(args.idx, bed->chrom);
        if ( id == -1 ) {
            warnings("Chromosome %s not found.", bed->chrom);
            continue;
        }

        hts_itr_t *itr = tbx_itr_queryi(args.idx, id, bed->start, bed->end);
        if ( itr == NULL ) {
            warnings("Region %s\t%d\t%d not found.", bed->chrom, bed->start, bed->end);
            continue;
        }
        int depth = 0; 
        int i;
        int total_bases_reg = 0;
        int length  = bed->end - bed->start;
        int uncover = bed->end - bed->start;
        assert(length > 0);
        args.total_length += length;
        
        //for ( i = 0; i < args.n_depth; ++i )
        //  args.depths_cutoff_per_reg[i] = 0;
        
        ksprintf(&reg, "%s\t%d\t%d\t", bed->chrom, bed->start, bed->end);
        
        while ( tbx_itr_next(args.fp_data, args.idx, itr, &str) >= 0 ) {
            depth = parse_depthData(str.s, str.l);
            args.total_base += depth;
            total_bases_reg += depth;
            uncover --;            
            for ( i = 0; i < args.n_depth; ++i) {
                if ( depth >= args.depths[i] ) {
                    args.depths_cutoff[i]++;
                    args.depths_cutoff_per_reg[i]++;
                }
                else {
                    break;
                }
            }
            str.l = 0;
        }
        double average_depth_reg = (double)total_bases_reg/length;
        ksprintf(&reg, "%.4f", average_depth_reg);
        for ( i = 0; i < args.n_depth; ++i ) {
            double cov = (double)args.depths_cutoff_per_reg[i]/length;
            ksprintf(&reg, "\t%.4f", cov);
        }

        if ( args.fp_output )
            fprintf(args.fp_output, "%s\n", reg.s);
        reg.l = 0;
        tbx_itr_destroy(itr);            
    }

    free(reg.s);
    
}

void summary_output()
{
    fprintf(args.fp_summary, "Total bases : %llu\n", args.total_base);
    int i;
    for ( i = 0; i < args.n_depth; ++i )
        fprintf(args.fp_summary, "Coverage above %d fold: %.4f\n", args.depths[i], (float)args.depths_cutoff[i]/args.total_length);    
}

int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    depths_retrieve();
    summary_output();
    memory_release();
    
    return 0;
}
    
