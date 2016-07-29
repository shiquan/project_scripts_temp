// gcc -I. -lhts -o vcfeva common.c vcfeva.c
//
//This program is used to evalue the differnece between the two samples, which merged in a VCF file.

#include "commons.h"
#include "htslib/vcf.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/kseq.h"

static char * Version = "0.1";
static char * outfile = NULL;
static char * report = NULL;
static bool dotasref = FALSE;
static bool export_uncover = FALSE;

static int bias = 0;

int usage()
{
    fprintf(stderr, "Usage: mergedvcfevalue [option] in.vcf.gz\n");
    fprintf(stderr, "                 --out  [FILE]   export the inconsistent positions in this file\n");
    fprintf(stderr, "                 --dotasref      treat the dot in the vcf file as a ref allele\n");
    fprintf(stderr, "                 --export_uncov  export the uncover positions in the outfile\n");
    fprintf(stderr, "                 --report [FILE] export the report to this file instead of stdout\n");
    fprintf(stderr, "                 --help          see this message\n");
    fprintf(stderr, "                 --version       show version\n");
    return 1;
}

int show_version()
{
    fprintf(stderr, "%s\n", Version);
    return 1;
}

htsFile *read_vcf_file(char * fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) errabort("Could not read file : %s", fname);
    return fp;
}

#define ALT_IS_REF 0
#define ALT_IS_HET 1
#define ALT_IS_HOM 2
#define ALT_IS_UNC 3

static char *mode = "zw";

int check_filename(char *fn)
{
    char *ss = fn;
    int length = strlen(ss);
    if ( length < 5 ) {
	warnings("unrecongnized file format name, only support vcf and vcf.gz till now. %s", fn);
	return 1;
    }
    if ( length > 4 && !strcmp(ss+length - 4, ".vcf") ) {
	mode = "w";
	return 0;
    } else if ( length > 7 && !strcmp(ss+length - 7, ".vcf.gz") ) {
	mode = "zw";
	return 0;
    } else {
	fn = (char*)realloc(fn, (length+8) *sizeof(char));
	memcpy(fn + length, ".vcf.gz", 7*sizeof(char));
	fn[length+7] = '\0';
	return 0;
    }
    return 1;
}

void write_report(uint32_t *m, bcf_hdr_t *hdr)
{
    int i, j;
    const char *title[] = { "REF", "HET", "HOM", "UNCOVER" };
    FILE *fp;
    if ( report ) {
	fp = fopen(report, "w");
	if ( fp == NULL ) errabort("%s : %s", report, strerror(errno));
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

int main(int argc, char **argv)
{
    int i, n;
    static struct option const long_opts[] =
    {
	{"out", required_argument, NULL, 1},
	{"report", required_argument, NULL, 2},
	{"dotasref", no_argument, NULL, 3},
	{"help", no_argument, NULL, 0},
	{"version", no_argument, NULL, 4},
	{"export_uncov", no_argument, NULL, 5}
    };
    bool help = FALSE;
    bool report_version = FALSE;
    while ((n = getopt_long(argc, argv, "1:2:304", long_opts, NULL)) >= 0)
    {
	switch (n)
	{
	case 1 : outfile = strdup(optarg); break;
	case 2 : report = strdup(optarg); break;
	case 3 : dotasref = TRUE; break;
	case 0 : help = TRUE; break;
	case 4 : report_version = TRUE; break;
	case 5 : export_uncover = TRUE; break;
	default : return 1;
	}
	if ( help ) return usage();
	if ( report_version ) return show_version();
    }
    n = argc - optind;
    if ( n > 1 ) errabort("only accept one input vcf");
    if ( export_uncover == TRUE && outfile == FALSE) {
	warnings("export uncove region only used with option --out");
	export_uncover = FALSE;
    }
    char * input;
    if ( n == 0 ) input = strdup("-");
    else input = strdup(argv[optind]);
    htsFile * fp = read_vcf_file(input);
    enum htsExactFormat fmt = hts_get_format(fp)->format;
    if ( fmt != vcf && fmt != bcf ) errabort("This is not a VCF/BCF file : %s", input);
    bcf_hdr_t * hdr = bcf_hdr_read(fp);
    int n_samples = bcf_hdr_nsamples(hdr);
    if ( n_samples != 2 ) errabort("the input VCF/BCF file must contain only two samples! %d", n_samples);
    LOG("Using sample %s as ref ...", hdr->samples[0]);
    LOG("Using sample %s as test ...", hdr->samples[1]);
    uint32_t matrix[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };
    bcf1_t * v = bcf_init1();
    kstring_t str = { 0, 0, 0 };
    uint32_t line = 0;
    htsFile *out = NULL;
    if ( outfile && !check_filename(outfile) ) out = hts_open(outfile, mode);
    if ( out != NULL ) bcf_hdr_write(out, hdr);
    while ( bcf_read1(fp, hdr, v) >= 0 )
    {
	bcf_unpack(v, BCF_UN_STR|BCF_UN_FMT);
	int k;
	str.l = 0;
	int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
	if ( !bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, tag_id) ) warnings("There is no 'GT' in the header!");
	for ( i = 0; i < v->n_fmt; ++i )
	    if ( v->d.fmt[i].id == tag_id ) break;
	if ( i == v->n_fmt ) {
	    vcf_format1(hdr, v, &str);
	    LOG("There is no tag GT in this line : %s", str.s);
	    continue;
	}
	corr_t xy[2] = { {-1, -2, -2}, {-1, -2, -2} };
	bcf_fmt_t * fmt = &v->d.fmt[i];

	for ( i = 0; i < 2; ++i )
	{
	    int corr = i;
	    if ( fmt == NULL ) {
		if ( dotasref == TRUE ) xy[corr].alt = ALT_IS_REF;
		else xy[corr].alt = ALT_IS_UNC;
		continue;
	    }
	    int last = -2;
	    uint8_t *d = (uint8_t*)((char*)fmt->p + fmt->size*i);
	    for ( k = 0; k < fmt->n && d[k] != (uint8_t)bcf_int8_vector_end; ++k )
	    {
		int curr = d[k]>>1;
		if ( last != curr ) {
		    if ( curr ) {
			if ( last == -2 ) xy[corr].alt = curr > 1 ? ALT_IS_HOM : ALT_IS_REF;
			else xy[corr].alt =  ALT_IS_HET;
		    } else {
			xy[corr].alt =  dotasref == TRUE ? ALT_IS_REF : ALT_IS_UNC;
		    }
		} else {
		    if ( curr ) {
			xy[corr].alt = curr > 1 ? ALT_IS_HOM : ALT_IS_REF;
		    } else {
			xy[corr].alt = dotasref == TRUE ? ALT_IS_REF : ALT_IS_UNC;
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
	if ( xy[0].alt != xy[1].alt && out != NULL) {
	    if ( xy[0].alt == ALT_IS_UNC || xy[1].alt == ALT_IS_UNC ) {
		if ( export_uncover == TRUE ) {
		    str.l = 0;
		    vcf_format1(hdr, v, &str);
		    vcf_write(out, hdr, v);
		}
	    } else {
		str.l = 0;
		vcf_format1(hdr, v, &str);
		vcf_write(out, hdr, v);
	    }
	}
	if ( xy[0].alt == ALT_IS_HET && xy[1].alt == ALT_IS_HET && (xy[0].min != xy[1].min || xy[0].max != xy[1].max ) ) {
	    bias++;
	    matrix[ALT_IS_HET][ALT_IS_HET]--;
	    if ( out != NULL ) {
		str.l = 0;
		vcf_format1(hdr, v, &str);
		vcf_write(out, hdr, v);
	    }
	}
	line++;
    }
    if ( out ) hts_close(out);
    if ( str.m ) free(str.s);
    write_report(matrix, hdr);
    bcf_hdr_destroy(hdr);
    free(input);
    bcf_destroy1(v);
    if ( outfile ) free(outfile);
    if ( report ) free(report);
    if ( hts_close(fp) ) warnings("hts_close returned non-zero status: %s", input);
    return 0;
}
