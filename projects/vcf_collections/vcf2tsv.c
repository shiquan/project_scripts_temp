// retrieve the genotype information in TSV format from a vcf/bcf file
// !! This is just a demo program to show newbies How To Use Htslib. To Export Vcf File into a tsv file, I perfer to use `bcftools query`
// shiquan@genomics.cn

#include "commons.h"
#include <sys/stat.h>
//#include "htslib/tbx.h"
#include "htslib/vcf.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/kseq.h"
//#include "htslib/kstring.h"

static int help = 0;

int usage()
{
    fprintf(stderr, "Usage: vcf2tsv in.vcf.gz > out.tsv\n");
    return 1;
}

htsFile *read_vcf_file(char * fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) errabort("Could not read file : %s", fname);
    return fp;
}

int main(int argc, char **argv)
{
    int i, n;
    if ( argc == optind ) return usage();
    char *fname = argv[optind];
    htsFile *fp = read_vcf_file(fname);
  
    enum htsExactFormat format = hts_get_format(fp)->format;
    if ( format != vcf && format != bcf ) errabort("only accept vcf/bcf file : %s", fname);

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    int n_samples = bcf_hdr_nsamples(hdr);
    kstring_t title = {0, 0, 0};
    if ( n_samples ) {
	ksprintf(&title,"#CHROM\tPOS\tREF");
	for (i = 0; i < n_samples; ++i)
	    ksprintf(&title, "\t%s", hdr->samples[i]);
	printf("%s\n", title.s);
	freemem(title.s);
    }
    bcf1_t *v = bcf_init1();
    while ( bcf_read1(fp, hdr, v) >= 0 )
    {
	int i;
	kstring_t str = {0, 0, 0};
	bcf_unpack(v, BCF_UN_STR|BCF_UN_FMT);
	ksprintf(&str, "%s\t%d\t%s",hdr->id[BCF_DT_CTG][v->rid].key ,v->pos+1, v->d.allele[0]);
	int j, k;
	int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
	if ( !bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, tag_id) ) errabort("There is not tag GT in the vcf header!");
	bcf_unpack(v, BCF_UN_ALL);

	for (i = 0; i < v->n_fmt; ++i)
	    if ( v->d.fmt[i].id == tag_id ) break;
	if (i == v->n_fmt ) {
	    for ( j = 0; j < n_samples; ++j ) kputs("\t.", &str); // no tag in this line
	} else {
	    bcf_fmt_t *fmt = &v->d.fmt[i];
	    for ( j = 0; j < n_samples; ++j )
	    {
		if (fmt == NULL) {
		    kputc('.', &str);
		    continue;
		}
		uint8_t *d = (uint8_t*)((char*)fmt->p + fmt->size * j);
		kputc('\t', &str);
		int last = -2;
		for ( k = 0; k < fmt->n && d[k] != (uint8_t)bcf_int8_vector_end; ++k )
		{
		    if ( last != d[k]>>1 ) {
			if ( k ) kputc("/|"[d[k]&1], &str); // phased or not
			if ( d[k] >> 1 ) {
			    kputs(v->d.allele[(d[k]>>1) -1], &str);
			}
			else
			    kputc('.', &str);
		    }
		    last = d[k]>>1;
		}
		if ( k == 0 ) kputc('.', &str);
	    }
	    printf("%s\n", str.s);
	}
	free(str.s);
    }
    
    bcf_hdr_destroy(hdr);
    bcf_destroy1(v);
    if ( hts_close(fp) ) errabort("hts_close returned non-zero status: %s", fname);
  
    return 0;
}
