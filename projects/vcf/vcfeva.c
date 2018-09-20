// gcc -I. -lhts -o vcfeva vcfeva.c
// 
//This program is used to evalue the differnece between the two samples, which merged in a VCF file.

// Reference
// https://genomebiology.biomedcentral.com/articles/10.1186/gb-2009-10-3-r32


#include "utils.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/bgzf.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "pkg_version.h"

enum bool {
    _promote_to_int = -1, // according to C07, enum should promote to int
    true = 0,
    false,
};

struct args {
    const char *input_fname;
    const char *output_fname;
    const char *report_fname;
    enum bool dotasref;
    enum bool export_uncov;
    int output_type;
} args = {
    .input_fname = 0,
    .output_fname = 0,
    .report_fname = 0,
    .dotasref = false,
    .export_uncov = false,
    .output_type = FT_VCF,
};


static int bias = 0;

int usage()
{
    fprintf(stderr, "Evaluate consistant of two samples in vcf file.\n");
    fprintf(stderr, "Version : %s\n", PROJECTS_VERSION);
    fprintf(stderr, "Usage: vcfeva [option] in.vcf.gz\n");
    fprintf(stderr, "      -dotasref      treat the dot in the vcf file as a ref allele\n");
    fprintf(stderr, "      -export_uncov  export the uncover positions in the outfile\n");
    fprintf(stderr, "      -report [FILE] export the report to this file instead of stdout\n");
    fprintf(stderr, "      -h, -help      see help information\n");
    fprintf(stderr, "      -o  [FILE]     export the inconsistent positions in this file\n");
    fprintf(stderr, "      -O  <u|v|z|b>  output type. default is u\n");
    return 1;
}


htsFile *read_vcf_file(char * fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp )
	error ("Could not read file %s : %s", fname, strerror(errno));
    return fp;
}

static const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF ) return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF ) return "wb";      // compressed BCF
    if ( file_type & FT_GZ ) return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}

int parse_args(int argc, char **argv)
{
    int i;
    const char *output_type = NULL;
    for (i = 1; i < argc; ) {
	const char *a = argv[i++];
	if ( strcmp(a, "-h") == 0 || strcmp(a, "-help") == 0)
	    return usage();
	if ( strcmp(a, "-dotasref") == 0 ) {
	    args.dotasref = true;
	    continue;
	}
	if ( strcmp(a, "-export_uncov") == 0 ) {
	    args.export_uncov = true;
	    continue;
	}
	const char **arg_var = 0;
	if ( strcmp(a, "-o") == 0 )
	    arg_var = &args.output_fname;
	else if ( strcmp(a, "-report") == 0 )
	    arg_var = &args.report_fname;
	else if ( strcmp(a, "-O") == 0 )
            arg_var = &output_type;

        if ( arg_var != 0 ) {
            if ( i == argc ) {
                error_print("Missing an argument after %s", a);
                return 1;
            }
            *arg_var = argv[i++];
            continue;
        }

        if ( args.input_fname == NULL ) {
            args.input_fname = a;
            continue;
        }
        error("Unknown argument %s", a);
    }


    if ( args.input_fname == NULL && (!isatty(fileno(stdin))) )
        args.input_fname = "-";

    if ( args.input_fname == NULL )
        error("No input file.");

    args.output_type = FT_VCF;
    if ( output_type != 0 ) {
	switch (output_type[0]) {
	    case 'b':
		args.output_type = FT_BCF_GZ; break;
	    case 'u':
		args.output_type = FT_BCF; break;
	    case 'z':
		args.output_type = FT_VCF_GZ; break;
	    case 'v':
		args.output_type = FT_VCF; break;
	    default :
		error("The output type \"%d\" not recognised\n", args.output_type);
	};
    }
    return 0;
}

#define ALT_IS_REF 0
#define ALT_IS_HET 1
#define ALT_IS_HOM 2
#define ALT_IS_UNC 3

void write_report(uint32_t *m, bcf_hdr_t *hdr)
{
    int i, j;
    const char *title[] = { "REF", "HET", "HOM", "UNCOVER" };
    FILE *fp;
    if ( args.report_fname ) {
	fp = fopen(args.report_fname, "w");
	if ( fp == NULL )
            error("%s : %s", args.report_fname, strerror(errno));
    } else {
	fp = stdout;
    }
    fprintf(fp,"ref sample: %s\n", hdr->samples[0]);
    fprintf(fp,"test sample: %s\n", hdr->samples[1]);
    fprintf(fp,"=============================================================================\n");
    for ( i = -1; i < 4; ++i )
    {
	for ( j=-1; j<4; ++j )
	{
	    if ( i == -1 ) {
		if ( j == -1 ) fprintf(fp, "%10s", "/");
		else fprintf(fp, "\t%10s", title[j]);
	    } else {
		if ( j == -1 ) fprintf(fp, "%10s", title[i]);
		else if ( j == ALT_IS_HET && i == ALT_IS_HET ) fprintf(fp, "\t%8u(%d)", m[5], bias);
		else fprintf(fp, "\t%10u",m[j*4+i]);
	    }
	}
	fprintf(fp, "\n");
    }
    fprintf(fp,"=============================================================================\n");
    float flsp = 0;
    float flsn = 0;
    float des_rate = 0;
    float cov_rate = 0;
    flsn = (float)(m[4] + m[8])/(m[4] + m[5] + m[6] + m[8] + m[9] + m[10] + bias);
    flsp = (float)(m[1] + m[2])/(m[1] + m[2] + m[5] + m[6] + m[9] + m[10] + bias);
    des_rate = (float)(m[6] + m[9])/(m[5] + m[6] + m[9] + m[10]);
    cov_rate = 1.0 - (float)(m[3] + m[7] + m[11])/(m[0] + m[1] + m[2] + m[3] + m[4] + m[5] + m[6] + m[7] + m[8] + m[9] + m[10] + m[11] + bias);
    fprintf(fp, "False Positive rate : %.5f\n", flsp);
    fprintf(fp, "False Negative rate : %.5f\n", flsn);
    fprintf(fp, "Variant discrepancy rate : %.5f\n", des_rate);
    fprintf(fp, "Coverage rate : %.5f\n", cov_rate);
    fclose(fp);
}

typedef struct {
    int alt;
    int min;
    int max;
} corr_t;

int evaluate_vcfs()
{
    htsFile *fp_input = hts_open(args.input_fname, "r");
    if ( fp_input == NULL )
        error("%s : %s.", args.input_fname, strerror(errno));
    enum htsExactFormat fmt = hts_get_format(fp_input)->format;

    if ( fmt != vcf && fmt != bcf )
        error("This is not a VCF/BCF file : %s.", args.input_fname);

    htsFile *fp_output = NULL;
    if ( args.output_fname != NULL ) {
        fp_output = hts_open(args.output_fname, hts_bcf_wmode(args.output_type));

        if ( fp_output == NULL )
            error("%s : %s.", args.output_fname, strerror(errno));
    }

    bcf_hdr_t *hdr = bcf_hdr_read(fp_input);
    int n_samples = bcf_hdr_nsamples(hdr);
    if ( n_samples != 2 )
        error("Input file must contain two samples. %d", n_samples);

    LOG_print("Using sample %s as reference.", hdr->samples[0]);
    LOG_print("Using sample %s as test.", hdr->samples[1]);

    uint32_t matrix[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
    bcf1_t * v = bcf_init1();
    kstring_t str = { 0, 0, 0 };

    if ( fp_output != NULL )
        bcf_hdr_write(fp_output, hdr);

    int i;
    uint64_t line = 0;
    while ( bcf_read1(fp_input, hdr, v) >= 0 ) {
	bcf_unpack(v, BCF_UN_STR|BCF_UN_FMT);
	int k;
	str.l = 0;
	int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
	if ( !bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, tag_id) )
            warnings("There is no 'GT' in the header!");
	for ( i = 0; i < v->n_fmt; ++i )
	    if ( v->d.fmt[i].id == tag_id ) break;
	if ( i == v->n_fmt ) {
	    vcf_format1(hdr, v, &str);
	    warnings("There is no tag GT in this line : %s", str.s);
	    continue;
	}
	corr_t xy[2] = { {-1, -2, -2}, {-1, -2, -2} };
	bcf_fmt_t * fmt = &v->d.fmt[i];

	for ( i = 0; i < 2; ++i ) {
	    int corr = i;
	    if ( fmt == NULL ) {
		if ( args.dotasref == true )
                    xy[corr].alt = ALT_IS_REF;
		else
                    xy[corr].alt = ALT_IS_UNC;
		continue;
	    }
	    int last = -2;
	    uint8_t *d = (uint8_t*)((char*)fmt->p + fmt->size*i);
	    for ( k = 0; k < fmt->n && d[k] != (uint8_t)bcf_int8_vector_end; ++k ) {
		int curr = d[k]>>1;
		if ( last != curr ) {
		    if ( curr ) {
			if ( last == -2 ) xy[corr].alt = curr > 1 ? ALT_IS_HOM : ALT_IS_REF;
			else xy[corr].alt =  ALT_IS_HET;
		    } else {
			xy[corr].alt =  args.dotasref == true ? ALT_IS_REF : ALT_IS_UNC;
		    }
		} else {
		    if ( curr ) {
			xy[corr].alt = curr > 1 ? ALT_IS_HOM : ALT_IS_REF;
		    } else {
			xy[corr].alt = args.dotasref == true ? ALT_IS_REF : ALT_IS_UNC;
		    }
		}
		if (last == -2 ) {
		    xy[corr].min = xy[corr].max = curr;
		} else {
		    if ( curr < xy[corr].min ) xy[corr].min = curr;
		    else if ( curr > xy[corr].max ) xy[corr].max = curr;
		}
		last = curr;
	    }
	}
	matrix[xy[0].alt][xy[1].alt]++;
	if ( xy[0].alt != xy[1].alt && fp_output != NULL) {
	    if ( xy[0].alt == ALT_IS_UNC || xy[1].alt == ALT_IS_UNC ) {
		if ( args.export_uncov == true) {
		    str.l = 0;
		    vcf_format1(hdr, v, &str);
		    vcf_write(fp_output, hdr, v);
		}
	    } else {
		str.l = 0;
		vcf_format1(hdr, v, &str);
		vcf_write(fp_output, hdr, v);
	    }
	}
	if ( xy[0].alt == ALT_IS_HET && xy[1].alt == ALT_IS_HET && (xy[0].min != xy[1].min || xy[0].max != xy[1].max ) ) {
	    bias++;
	    matrix[ALT_IS_HET][ALT_IS_HET]--;
	    if ( fp_output != NULL ) {
		str.l = 0;
		vcf_format1(hdr, v, &str);
		vcf_write(fp_output, hdr, v);
	    }
	}
	line++;
    }
    
    write_report(matrix, hdr);

    bcf_hdr_destroy(hdr);
    hts_close(fp_input);
    if ( fp_output )
        hts_close(fp_output);

    if ( str.m )
        free(str.s);
    
    return 0;
}

int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;

    if ( evaluate_vcfs() )
        return 1;

    return 0;
}


