#include <stdio.h>
#include <stdlib.h>
#include <htslib/kstring.h>
#include <zlib.h>

struct args {
    const char *input_fname;
    const char *output_fname;
    const char *name;
    const char *color;
    const char *visibility;    
};
int parse_args(int argc, char **argv)
{
    if (argc == 1)
	error("Usage: depth2wig depth.tsv.gz -o out.wig");

    const char *a = 0;
    int i;
    for (i = 1; i < argc; ) {
    }
}
int main(int argc, char **argv)
{
    parse_args(argc, args);

    export_wig();
    return 0;
}
