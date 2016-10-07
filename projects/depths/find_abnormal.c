// find abnormal points from tabix indexed depth files
// find_abnormal -target target.bed -mode right|left|both -cutoff [3sigma] depth.txt.gz -o out.txt
#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>
#include <zlib.h>
#include "utils.h"
#include "options.h"
#include "bed_utils.h"

enum mode {
    both_mode,
    right_mode,
    left_mode,
};

struct args {
    const char *input_fname;
    const char *target_fname;
    const char *output_fname;
    float cutoff;
    enum mode mode;
    struct bedaux *target_aux;
};

struct args args = {
    .input_fname = 0,
    .target_fname = 0,
    .output_fname = 0,
    .cutoff = 3.0,
    .mode = both,
    .target_aux = NULL,
};

const struct option_spec options[] = {
    { "-target", argument_path, argument_optional, "target regions not set, stat all positions"},          
    { "-mode", argument_string, argument_optional, "no mode set, use both sides mode for test"},
    { "-cutoff", argument_float, argument_optional, "cutoff value for abnormal value, default is 3 sigma"},
    { "-out", argument_path, argument_optional, "export stdout in default"},
    { "-h", argument_bool, argument_optional, NULL},
    { "-help", argument_bool, argument_optional, NULL},
    { NULL, _unknown, _unknown, NULL},
};

static int target_is_set = 0;

int parse_args(int argc, char **args)
{
    if ( init_args(options, argc, args) )
        return 1;

    if ( argument_exist("-h") ||  argument_exist("-help") )
        return usage();
    
    if ( argc == 1 & CHECK_STDIN ) {
        args.input_fname = "-";
    } 

    if ( argc >= 2 ) {
        args.input_fname = argv[1];
        if ( argc > 2) 
            warnings("Only accept one input file.");
    }

    if ( args.input_fname == NULL )
        error("No input file found. Use -h for more information.");

    if ( argument_exist("-cutoff") ) 
        args.cutoff = get_argument_float("-cutoff");
    
    if ( argument_exist("-mode") ) {
        char *mode = get_argument_string("-mode");
        if ( strcmp(mode, "left") == 0 )
            args.mode = left_mode;
        else if ( strcmp(mode, "right") == 0 )
            args.mode = right_mode;
        else if ( strcmp(mode, "both") == 0 )
            args.mode = both_mode;
        else
            error("Unknown mode type %s; [right | left | both ]");
    }

    if ( argument_exist("-out") ) 
        args.out = get_argument_path("-out");

    if ( argument_exist("-target") ) {
        args.target = get_argument_path("-target");
        args.target_aux = bedaux_init();
        bed_read(args.target_aux, args.target);
        target_is_set = 1;
    }
    
    return 0;    
}

int cutoff_check()
{
    if ( target_is_set ) {
        tbx_t *idx = NULL;
        idx = tbx_index_load(args.input_fname);
        if ( idx == NULL )
            error("Failed to load tabix index of %s.", args.input_fname);

        struct bed_line line = BED_LINE_INIT;
        
        return 0;
    }

    
    return 0;
}

int export_abnormal()
{
    return 0;
}

void memory_release()
{
    if ( target_is_set ) 
        bed_destroy(args.target_aux);
    
}
int main(int argc, char **argv)
{
    // parse arguments
    parse_args(argc, argv);
    // calculate the cutoff depths for low or high boundaries
    cutoff_check();
    // export the abnormal positions
    export_abnormal();
    
    return 0;
}
