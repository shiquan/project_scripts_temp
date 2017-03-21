#include "utils.h"
#include "hgvs.h"
#include "hgvs_vcf.h"
#include "genepred.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include <string.h>

struct args {

    const char *fname_input;
    const char *data;
    const char *fasta;
    
    htsFile *fp_input;
    htsFile *fp_output;
    bcf_hdr_t *hdr_in;
    bcf_hdr_t *hdr_out;
    
} args = {
    .fname_input = NULL,
    .data = NULL,
    .fasta = NULL,
    
    .fp_input = NULL,
    .fp_output = NULL,
    .hdr_in = NULL,
    .hdr_out = NULL,    
};

static const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF )
	return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF )
	return "wb";      // compressed BCF
    if ( file_type & FT_GZ )
	return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}

int usage()
{
    fprintf(stderr, "vcf_hgvs\n"
            " -data        < genepred data >\n"
            " -fasta       < transcript rna fasta file >\n"
            " -gene_list   [ genes list file]\n"
            " -trans_list  [ transcripts list file]\n"
            "  input.vcf.gz[bcf/vcf]\n"
            "\n"
            "shiquan@genomics.cn\n"
        );
    
    return 1;
}

int parse_args(int ac, char **av)
{
    if ( ac == 0 )
        return usage();

    const char *fname_output = NULL;
    const char *out_type_string = NULL;
    const char *genes_fname = NULL;
    const char *transcripts_fname = NULL;
    
    int i;

    for ( i = 0; i < ac; ) {
        const char *a = av[i++];
        const char **arg_var = 0;
        
        if ( strcmp(a, "-h") == 0 )
            return usage();

        if ( strcmp(a, "-data") == 0 )
            arg_var = &args.data;
        else if ( strcmp(a, "-fasta") == 0 )
            arg_var = &args.fasta;
        else if ( strcmp(a, "-o") == 0 )
            arg_var = &fname_output;
        else if ( strcmp(a, "-O") == 0 )
            arg_var = &out_type_string;
        else if ( strcmp(a, "-gene_list") == 0 )
            arg_var = &genes_fname;
        else if ( strcmp(a, "-trans_list") == 0 )
            arg_var = &transcripts_fname;
        
        if ( arg_var != 0 ) {
            if ( i == ac )
                error("Miss an argument after %s.", a);
            *arg_var = av[i++];
            continue;
        }

        if ( args.fname_input == 0 ) {
            args.fname_input = a;
            continue;
        }

        error("Unknown argument %s.", a);            
    }

    if ( args.data == NULL )
        error("-data required genepred database.");

    if ( args.fasta == NULL )
        error("-fasta required transcripts fasta file.");

    // Load input file.
    if ( args.fname_input == NULL && (!isatty(fileno(stdin)) ) ) {
        args.fname_input = "-";        
    }
    args.fp_input = hts_open(args.fname_input, "r");

    htsFormat type = *hts_get_format(args.fp_input);
    if ( type.format != vcf && type.format != bcf )
        error("Unsupported input format %s.", args.fname_input);

    int out_type = FT_VCF;
    if (out_type_string != 0) {
	switch (out_type_string[0]) {
	    case 'b':
		out_type = FT_BCF_GZ; break;
	    case 'u':
		out_type = FT_BCF; break;
	    case 'z':
		out_type = FT_VCF_GZ; break;
	    case 'v':
		out_type = FT_VCF; break;
	    default :
		error("The output type %s not recognised\n", out_type_string);
	};
    }

    args.hdr_in = bcf_hdr_read(args.fp_input);

    args.fp_output = fname_output == 0 ? hts_open("-", hts_bcf_wmode(out_type)) :
        hts_open(fname_output, hts_bcf_wmode(out_type));

    args.hdr_out = bcf_hdr_dup(args.hdr_in);

    if ( init_hgvs_anno(args.data, args.fasta, args.hdr_out) )
        return 1;

    if ( genes_fname != NULL ) {
        if ( set_genes_list ( genes_fname ) ) {
            warnings("Failed to load gene list : %s", genes_fname);
        }
    }

    if ( transcripts_fname != NULL ) {
        if ( set_transcripts_list ( transcripts_fname ) ) {
            warnings("Failed to load transcripts list : %s", transcripts_fname);
        }
    }
    
    return 0;
}

int memory_release()
{
    hts_close(args.fp_input);
    hts_close(args.fp_output);
    bcf_hdr_destroy(args.hdr_in);
    bcf_hdr_destroy(args.hdr_out);
    return 1;
}

int annotate_hgvs()
{

    bcf_hdr_write(args.fp_output, args.hdr_out);
    bcf1_t *line = bcf_init();
    for ( ;; ) {
        if ( bcf_read(args.fp_input, args.hdr_in, line) != 0 )
            break;
        setter_hgvs_vcf(args.hdr_out, line);
        bcf_write1(args.fp_output, args.hdr_out, line);
    }
    bcf_destroy(line);
    return 0;
}

int main(int argc, char **argv)
{
    if ( parse_args(--argc, ++argv) )
        return 1;

    annotate_hgvs();

    memory_release();

    return 0;
}
