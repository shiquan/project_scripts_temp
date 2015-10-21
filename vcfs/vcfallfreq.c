/* vcfallfreq.c - count the allele freq of VCF/BCF
 * SHI Quan (shiquan@genomics.cn)
 */

#include "commons.h"
#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
int usage()
{
  fprintf(stderr, "Usage: vcfallfreq in.vcf.gz\n");
  return 1;
}
htsFile *read_vcf_file(char *fn)
{
  htsFile *fp = hts_open(fn, "r");
  if ( !fp ) errabort("Could not read file : %s", fn);
  assert(fp->format.format == vcf || fp->format.format == bcf);
  return fp;
}

struct _cnt_t {
    int id;
    uint64_t n;
};

int main(int argc, char ** argv)
{
  int n;
  n = argc - optind;
  if ( n>1 ) errabort(" only accept one input vcf file");
  char *input;
  if ( n == 0 ) input = strdup("-");
  else input = strdup(argv[n]);
  htsFile *fp = read_vcf_file(input);
  free(input);
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  int n_samples = bcf_hdr_nsamples(hdr);
  assert(n_samples>1);
  LOG("Total Samples : %d\n", n_samples);
  bcf1_t *v = bcf_init();
  fprintf(stdout, "#Chr\tPos\tRef\tAlt\tN_samples\tN_alleles\tHom_ref\tHet\tHom_alt\tminor_allele\tMinor_allele_count\tMAF\tmajor_alele\tmajor_allele_count\tmajorAF\n");
  int i, k;
  kstring_t str = {0, 0, 0};
  while ( bcf_read1(fp, hdr, v) >= 0 )
    {
      // output : chr pos ref alt N_sam N_alleles ref het hom minorA minorF majorA majorF
      bcf_unpack(v, BCF_UN_STR|BCF_UN_FMT);
      str.l =0;
      int gt = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
      if ( !bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gt)) warnings("There is no 'GT' tag in the header!");
      for ( i=0; i<v->n_fmt; ++i ) if ( v->d.fmt[i].id==gt ) break;
      if ( i==v->n_fmt )
	{
	  vcf_format1(hdr, v, &str);
	  LOG("There is no tag 'GT' in this line : %s", str.s);
	  continue;
	}
      bcf_fmt_t *fmt = &v->d.fmt[i];
      int nal = v->n_allele;
      int ref = 0, het = 0, hom = 0;
      int n_sams = 0;
      uint64_t cnt[nal];
      for ( i=0; i<nal; ++i ) cnt[i]=0;
      for ( i=0; i<n_samples; ++i )
	{
	  uint8_t *d = (uint8_t*)((char*)fmt->p+fmt->size*i);
	  for ( n=0; n<fmt->n && d[n] != (uint8_t)bcf_int8_vector_end; ++n);
	  if ( n==0 ) continue;
	  n_sams++;
	  if ( n==1 ) {
	    int g = (d[0]>>1)-1;
	      if ( g==0 ) ref++; 
	      else hom++;
	      cnt[g] += 2; // 2 copy number
	  }
	  else {
	    int last = -2;
	    for ( k=0; k<n; ++k)
	      {
		int g = (d[k]>>1)-1;
		//fprintf(stderr, "g: %d\n", g);
		if ( last != g ) {
		  if ( last == -2 ) last=g;
		  else het++;
		} 
		else {
		  if (g==0) ref++;
		  else hom++;
		}
		cnt[g]++;
	      }
	  }
	}
      kputs(hdr->id[BCF_DT_CTG][v->rid].key, &str); // CHROM
      kputc('\t', &str); kputw(v->pos+1, &str); // POS
      kputc('\t', &str); // REF
      if (v->n_allele > 0) kputs(v->d.allele[0], &str);
      else kputc('.', &str);
      kputc('\t', &str); // ALT
      if (v->n_allele > 1) {
	for ( k=1; k<v->n_allele; ++k) {
	  if ( k>1 ) kputc(',', &str);
	  kputs(v->d.allele[k], &str);
	}
      }
      else kputc('.', &str);
      kputc('\t', &str); kputw(n_sams, &str); // N_sams
      struct _cnt_t minor = {0, 0};
      struct _cnt_t major = {0, 0};
      uint64_t all=0;
      for ( k=0; k<nal; ++k ) 
	{
	  all+=cnt[k];
	  if (k==0) {
	      minor.n = cnt[0];
	      minor.id = 0;
	  }
	  if (minor.n > cnt[k]) {
	      minor.n=cnt[k];
	      minor.id=k;
	  }
	  if (major.n < cnt[k]) {
	      major.n=cnt[k];
	      major.id=k;
	  }
	}
      if (minor.id == major.id) {
	  minor.id = -1;
	  minor.n = 0;
      }
      double minF=0, majF=0;
      if (all) {
	minF = (double)minor.n/all;
	majF = (double)major.n/all;
      }
      kputc('\t', &str); kputw(all, &str); // all alleles
      kputc('\t', &str); kputw(ref, &str); // ref
      kputc('\t', &str); kputw(het, &str); // het
      kputc('\t', &str); kputw(hom, &str); // hom
      kputc('\t', &str);
      if (minor.id == -1) kputc('.', &str); // minor allele
      else kputs(v->d.allele[minor.id], &str); // minor allele
      kputc('\t', &str); kputw(minor.n, &str); // minorN
      kputc('\t', &str); ksprintf(&str, "%g", minF); // minorF
      kputc('\t', &str); kputs(v->d.allele[major.id], &str); // major allele
      kputc('\t', &str); kputw(major.n, &str); // majorN
      kputc('\t', &str); ksprintf(&str, "%g", majF); // majorF
      fprintf(stdout, "%s\n", str.s);
    }
  bcf_hdr_destroy(hdr);
  bcf_destroy(v);
  if (str.m) free(str.s);
  if ( hts_close(fp) ) warnings("hts_close return non-zero status");
  return 0;
}
