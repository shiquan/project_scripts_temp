/* Author:    Shi Quan            */
/* Email:     shiquan@genomcis.cn */

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <getopt.h>
// #include <pthread.h>
#include <zlib.h>

#include "kseq.h"
#include "kstring.h"

static char * program_name =  "dyncutadaptor";
static char * Version = "v0.1.2";

unsigned char bam_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

char *bam_nt16_rev_table = "=ACMGRSVTWYHKDBN";

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };

// Illumina
static uint8_t illumina[19] = { 2, 8, 4, 8, 2, 8, 2, 8, 8, 1, 8, 1, 2, 1, 2, 1, 8, 2, 8};
//static uint8_t illumina[10] = { 2, 8, 4, 8, 2, 8, 2, 8, 8, 1 };
static int length_adaptor = 10;
static int slave = 1;
static int minimum = 0;
static int bias = 1;
static int cest = 3;

void seq_comp(uint8_t * seq, int length)
{
	int i;
	for (i = 0; i < length>>1; ++i) {
		int8_t t = seq_comp_table[seq[length - 1 - i]];
		seq[length - 1 - i] = seq_comp_table[seq[i]];
		seq[i] = t;
	}
	if (length & 1) seq[i] = seq_comp_table[seq[i]];
}

KSEQ_INIT(gzFile, gzread)

static char * out1 = 0;
static char * out2 = 0;

uint8_t * seq2code(char * str, int n) 
{
	int i;
	uint8_t *a;
	a = calloc(n, sizeof(uint8_t));
	for(i = 0; i < n; ++i) {
		a[i] = bam_nt16_table[str[i]];
	}
	return a;
}

static int * BMprep(const uint8_t * pat, int m)
{
	int i, *suff, *prep, *bmGs, *bmBc;
	prep = (int*) calloc(m + 15, sizeof(int));
	bmGs = prep; bmBc = prep + m;
	for (i = 0; i < 15; ++i) bmBc[i] = m;
	for (i = 0; i < m - 1; ++i) bmBc[pat[i]] = m - i - 1;
	suff = (int*) calloc(m, sizeof(int));
	int f = 0, g;
	suff[m - 1] = m;
	g = m - 1;
	for(i = m - 2; i >= 0; --i) {
		if (i > g && suff[i + m - 1 - f] < i - g) 
			suff[i] = suff[i + m - 1 - f];
		else {
			if (i < g) g = i;
			f = i;
			while (g >= 0 && pat[g] == pat[g + m - 1 - f]) --g;
			suff[i] = f - g;
		}
	}
	int j = 0;
	for (i = 0; i < m; ++i) bmGs[i] = m;
	for (i = m -1; i >= 0; --i)
		if (suff[i] == i + 1)
			for (; j < m - 1 - i; ++j)
				bmGs[j] = m - 1 - i;
	for (i = 0; i <= m - 2; ++i)
		bmGs[m - 1 - suff[i]] = m - 1 -i;
	free(suff);
	return prep;
}

int location(const uint8_t * str, int n, const uint8_t * pat, int m, int *prep)
{
	int i, j, *bmGs, *bmBc;
	//prep = BMprep(pat, m);
	bmGs = prep;
	bmBc = prep + m;
	j  = 0;
	while ( j <= n - m) {
		int b = 0, k;
		for ( i = m - 1; i >= 0; --i) {
			if (!(pat[i] & str[i+j])) {
				if (b) {
					i = k;
					break;
				} else {
					k = i;
					b++;
				}
			}
		}
		if (i >= 0) {
			int max = bmBc[str[i+j]] - m + 1 + i;
			if (max < bmGs[i]) max = bmGs[i];
			j += max;
		} else {
			//free(prep);
			return j;
		}
	}
	j = n - m;
	while (j < n - cest) {
		for (i = 0 ; i < n - j && pat[i] == str[j + i]; ++i);
		if ( i < n - j) ++j;
		else {
			//free(prep);
			return j;
		}
	}
	//free(prep);
	return 0;
}

int loc_adaptor(uint8_t * str, int n, uint8_t * pat, int m, int * prep, int *k)
{
	*k = location(str, n, pat, m, prep);
	if ( *k) return 1;
	seq_comp(str, n);
	*k = location(str, n, pat, m, prep);
	if ( *k) return -1;
	return 0;		
}

int check_loc(uint8_t * seq, int l, uint8_t * a, int m, int loc) 
{
	int i, j, k = 0;
	j = l - loc > m ? m + loc: l;
	uint8_t * t = a;
	for (i = loc; i < j; ++i) {
		if(seq[i] != *t) k ++;
		t++;
	}
	if ( k > j -2 || k > bias) return 0;
	return 1;
}

void seq_reverse(char * seq, int len) 
{
	int i;
	for ( i = 0; i < len>>1; ++i) {
		unsigned char t = seq[i];
		seq[i] = seq[len - i -1];
		seq[len - i - 1] = t;
	}
}

int seq_reloc(char * seq, unsigned long * len, unsigned long loc) 
{
	seq_reverse(seq, *len);
	seq[loc] = 0;
	seq_reverse(seq, loc);
	*len = loc;
	return 1;
}

long long filter_reads = 0;

int cut_adaptor(kseq_t * seq1, kseq_t *seq2, uint8_t * pat, int len, int * prep) 
{
	int m = 0, n = 0;
	int l, loc = 0;		
	uint8_t * s, * p;
	s = seq2code(seq1->seq.s, seq1->seq.l);
	p = seq2code(seq2->seq.s, seq2->seq.l);
	int length = seq1->seq.l;
	if (slave) {
		m = location(s, seq1->seq.l, pat, len, prep);
		if (check_loc(s, seq1->seq.l, pat, len, m)) 
			loc = m;
		else {
			n = location(p, seq2->seq.l, pat, len, prep);
			if(check_loc(p, seq2->seq.l, pat, len, n)) loc = n;
		}
		if (loc) {
			if (minimum && loc < minimum) goto MINI;
			seq1->seq.s[loc] = 0; seq1->seq.l = l;
			seq1->qual.s[loc] = 0; seq1->qual.l = l;
			seq2->seq.s[loc] = 0; seq2->seq.l = l;
			seq2->qual.s[loc] = 0; seq2->qual.l = l;
		} else {
			seq_comp(s, seq1->seq.l); seq_comp(p, seq2->seq.l);
			m = location(s, seq1->seq.l, pat, len, prep);
			if (check_loc(s, seq1->seq.l, pat, len, m)) loc = m;
			else {
				n = location(p, seq2->seq.l, pat, len, prep);
				if(check_loc(p, seq2->seq.l, pat, len, n)) loc = n;
			}
			if (loc) {
				if (minimum && loc < minimum) goto MINI;
				seq_reloc(seq1->seq.s, &seq1->seq.l, loc);
				seq_reloc(seq2->seq.s, &seq2->seq.l, loc);
				seq_reloc(seq1->qual.s, &seq1->qual.l, loc);
				seq2->qual.s[loc] = 0; seq2->qual.l = l;
			} else {
				free(s); free(p);
				return 0;
			}
		}
	} else {
		m = location(s, seq1->seq.l, pat, len, prep);
		n = location(p, seq2->seq.l, pat, len, prep);
		if (m || n) {
			loc = m > n ? n > 0 ? n : m : m > 0 ? m : n;
			if (minimum && loc < minimum) goto MINI;
			seq1->seq.s[loc] = 0; seq1->seq.l = l;
			seq1->qual.s[loc] = 0; seq1->qual.l = l;
			seq2->seq.s[loc] = 0; seq2->seq.l = l;
			seq2->qual.s[loc] = 0; seq2->qual.l = l;
		} else { 
			seq_comp(s, seq1->seq.l); seq_comp(p, seq2->seq.l);
			m = location(s, seq1->seq.l, pat, len, prep);
			n = location(p, seq2->seq.l, pat, len, prep);
			if (m || n) {
				loc = m > n ? n > 0 ? n : m : m > 0 ? m : n;
				if (minimum && loc < minimum) goto MINI;
				seq_reloc(seq1->seq.s, &seq1->seq.l, loc);
				seq_reloc(seq2->seq.s, &seq2->seq.l, loc);
				seq_reloc(seq1->qual.s, &seq1->qual.l, loc);
				seq2->qual.s[loc] = 0; seq2->qual.l = l;
			} else {
				free(s); free(p);
				return 0;
			}
		}
	}
	free(s); free(p);
	return 1;
MINI:
	filter_reads++;
	free(s); free(p);
	return -1;
}

int check_name(kstring_t name1, kstring_t name2)
{
	char * s1 = name1.s;
	char * s2 = name2.s;
	size_t n;
	for(n = 0; n < name1.l - 1; ++n, ++s1, ++s2) if (*s1 != *s2) return 0;
	return 1;
}

int loadfastq_pe(const char * pe1, const char * pe2, uint8_t * a)
{
	gzFile fq1, fq2;
	kseq_t * seq1, * seq2;
	int l, m;
	long long all_reads = 0;
	long long cutted_reads = 0;
	fq1 = gzopen(pe1, "r"); fq2 = gzopen(pe2, "r");
	if (fq1 == NULL) {
		fprintf(stderr, "[loadfastq_pe] %s : %s\n", pe1, strerror(errno));
		return 0;
	}
	if (fq2 == NULL) {
		fprintf(stderr, "[loadfastq_pe] %s : %s\n", pe2, strerror(errno));
		return 0;
	}
	seq1 = kseq_init(fq1); seq2 = kseq_init(fq2);
	FILE * fp1, * fp2;
	fp1 = fopen(out1, "w"); fp2 = fopen(out2, "w");
	int *prep;
	prep = BMprep(a, length_adaptor);
	while (kseq_read(seq1) > 0 && kseq_read(seq2) > 0) {
		all_reads++;
		if (check_name(seq1->name, seq2->name) < 1) {
			fprintf(stderr, "[loadfastq_pe] The fastq files not paired:\n %s\n%s\n",
				seq1->name.s, seq2->name.s);
			goto LOADFQ_ERROR;
		}
		int n;
		n = cut_adaptor(seq1, seq2, a, length_adaptor, prep);
		if (n) cutted_reads++;
		if (n < 0) continue;
		fputc('@', fp1); fputs(seq1->name.s, fp1); fputc('\n', fp1); 
		fputs(seq1->seq.s, fp1); fputs("\n", fp1); fputs("+\n",fp1);
		fputs(seq1->qual.s, fp1); fputs("\n", fp1);
		fputc('@', fp2); fputs(seq2->name.s, fp2); fputs("\n", fp2); 
		fputs(seq2->seq.s, fp2); fputs("\n", fp2); fputs("+\n",fp2);
		fputs(seq2->qual.s, fp2); fputs("\n", fp2);
	}
	free(prep);
	kseq_destroy(seq1); kseq_destroy(seq2);
	gzclose(fq1); gzclose(fq2);
	fclose(fp1); fclose(fp2);
	fprintf(stdout, "dealed with %lld reads.\nTotal reads is %lld\nfilter %lld reads\n", cutted_reads, all_reads, filter_reads);
	return 1;

LOADFQ_ERROR:
	kseq_destroy(seq1);
	kseq_destroy(seq2);
	gzclose(fq1);
	gzclose(fq2);
	return 0;
}


int usage(int s)
{
	if (s) {
		fprintf(stderr,"\
Usage: %s --fastq1 [in1.fq.gz] --fastq2 [in2.fq.gz] \n\
======================================================\n\
--fastq1, -f         Fastq file of read1.\n\
--fastq2, -r         Fastq file of read2.\n\
--outfq1, -o         New fastq file of read1.\n\
--outfq2, -p         New fastq file of read2.\n\
--seed,   -s         Initial length of adaptor.[10]\n\
--slave,  -d         Cut all sequence like adaptor.\n\
--mis,    -i         Tolerate mismatchs.[1]\n\
--tail,   -t         Don't cut tail in the last n bp.[3]\n\
--min,    -m         Don't keep seqences shorter.[0]\n\
--help,   -h         See this information.\n\
======================================================\n\
Author: SHI Quan (shiquan@genomics.cn)\n\
Version: %s\n\
", program_name, Version);
	} else {
		fprintf(stderr,"\
Usage: %s --fastq1 [in1.fq.gz] --fastq2 [in2.fq.gz]\n\
use --help for more information.\n\
",program_name);
	}
	return 0;
}

static struct option const long_opts[] = {
	{"fastq1", required_argument, NULL, 'f'},
	{"fastq2", required_argument, NULL, 'r'},
	{"outfq1", required_argument, NULL, 'o'},
	{"outfq2", required_argument, NULL, 'p'},
	{"seed",   required_argument, NULL, 's'},
	{"mis",    required_argument, NULL, 'i'},
	{"tail",   required_argument, NULL, 't'},
	{"slave",  no_argument, NULL, 'd'},
	{"min",    required_argument, NULL, 'm'},
	{"help", no_argument, NULL, 'h'}
};

int main(int argc, char *argv[])
{
	char *fastq1 = 0;
	char *fastq2 = 0;
	int help = 0;
	int n, i;
	while ((n = getopt_long(argc, argv, "f:r:o:p:s:dhi:t:", long_opts, NULL)) >= 0) {
		switch (n) {
		case 'o': out1 = optarg; break;
		case 'p': out2 = optarg; break;
		case 'f': fastq1  = optarg;  break;
		case 'r': fastq2  = optarg;  break;
		case 'h': help = 1; break;
		case 't': cest = atoi(optarg); break;
		case 'i': bias =atoi(optarg); break;
		case 's': length_adaptor = atoi(optarg); break;
		case 'd': slave = 0; break;
		case 'm': minimum = atoi(optarg); break;
		default:  return usage(0);
		}
		if (help) return usage(1);
	}
	if (length_adaptor > 19 || length_adaptor < 5) {
		fprintf(stderr, "seed must be in [5, 19]\n");
		return 0;
	}
	if (minimum > 40 || minimum < 0) {
		fprintf(stderr, "min must be in [1,40]\n");
		return 0;
	}
	if (bias > 5 || bias < 0) {
		fprintf(stderr, "mismatch must below 5!\n");
		return 0;
	}
	if (cest > 5 || cest < 0) {
		fprintf(stderr, "not support keep tail longer than 5 bp.\n");
		return 0;
	}
	if (!fastq1 || !fastq2 ) return usage(0);
	if (!out1) out1 = "reads1.fq";
	if (!out2) out2 = "reads2.fq";
	loadfastq_pe(fastq1, fastq2, illumina);
	return 1;
}
