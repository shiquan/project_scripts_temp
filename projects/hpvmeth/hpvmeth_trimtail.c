/* 
 *  This program used to cut tails of torrent raw sequence for HPV-meth project 
 *  Usage : hpvmeth_trimtail in.bam out.fq.gz
 *  Author: shiquan@genomics.cn 2015/10/21
 *  
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>

#define SKIP_LENGTH 30
#define BARCODE_GAP 12 // 9+ 3
#define MAXTAIL 25

/*
  struct of adaptors:
  5' GTATGGTAGGTTGATGTGATTTAGTxxxxxxxxxCGA
  3' TGTGTTGTTTGATTGTGTTGAGTTGxxxxxxxxxGTT
 */
//const unsigned char tail5l[] = {2,3,0,3,2,2,3,0,2,2,3,3,2,0,3,2,3,2,0,3,3,3,0,2,3,}
const uint8_t tail5l[] = {1, 8, 8};
const uint8_t tail5r[] = {8, 4, 1};
//const char tail3l[] = {3,2,3,2,3,3,2,3,3,3,2,0,3,3,2,3,2,3,3,2,0,2,3,3,2,};
const uint8_t tail3l[] = {8, 4, 1};
const uint8_t tail3r[] = {4, 8, 8};

uint8_t seq_nt16_rev_table[] = { 0, 8, 4, 0, 2,  0, 0, 0, 1 };
int usage()
{
    fprintf(stderr, "Usage: \n");
    fprintf(stderr, "hpvmeth_trimtail in.bam out.fq\n");
    fprintf(stderr, "This program was designed to cut library prepare tails in the HPV-methylation project.\n");
    fprintf(stderr, "shiquan@genomics.cn\n");
    return 1;
}
void reverse(uint8_t *seq, int qlen)
{
    int i;
    for (i=0; i<qlen;++i)
    {
	uint8_t t = seq[i];
	seq[i] = seq[qlen-i-1];
	seq[qlen-i-1] = t;
    }
}
int trim_tail(uint8_t *buf, uint8_t *qual, int qlen, kstring_t *string)
{
    int i, j1, j2;
    int t1, t2;
    int l = 5;
    int b = 0;
    for (i=0, j1=0; i<MAXTAIL;)
    {

	if ( ((buf[i] ==tail5l[j1]) && (i+BARCODE_GAP>qlen-1 ? 1 : buf[i+BARCODE_GAP]==tail5r[j1]))
	    || ((buf[i] ==seq_nt16_rev_table[tail5l[j1]]) && (i+BARCODE_GAP>qlen-1 ? 1 : buf[i+BARCODE_GAP]==seq_nt16_rev_table[tail5r[j1]]))) {
	    j1++;
	    if (j1==3) break;

	}
	else {
	    i = j1 > 0? i--: i;
	    j1=0;
	}
	i++;
    }
    for (i=0, j2=0; j1 != 3 && i<MAXTAIL;)
    {
	if ( ((buf[i] ==tail3l[j1]) && (i+BARCODE_GAP>qlen-1 ? 1 : buf[i+BARCODE_GAP]==tail3r[j2]))
	    || ((buf[i] ==seq_nt16_rev_table[tail5l[j1]]) && (i+BARCODE_GAP>qlen-1 ? 1 : buf[i+BARCODE_GAP]==seq_nt16_rev_table[tail3r[j2]]))) {
	    j2++;
	    if (j2==3) break;
	}
	else {
	    i = j2 > 0? i-- : i;
	    j2=0;
	}
	++i;
    }

    if (j1>j2) {
	l = 5;
	t1 = j1;
    }
    else {
	l = 3;
	t1 = j2;
    }
	
    t1 = t1==0 ? 0 : i-t1+BARCODE_GAP;

    reverse(buf, qlen);
    for (i=0, j1=0; i<MAXTAIL;)
    {
	if ( l==5 && (((buf[i] ==tail5l[j1]) && (i+BARCODE_GAP>qlen-1 ? 1 : buf[i+BARCODE_GAP]==tail5r[j1]))
		      || ((buf[i] ==seq_nt16_rev_table[tail5l[j1]]) && (i+BARCODE_GAP>qlen-1 ? 1 : buf[i+BARCODE_GAP]==seq_nt16_rev_table[tail5r[j1]])))) {
	    j1++;
	    if (j1==3) break;
	}
	else {
	    i = j1>0 ? i-- : i;
	    j1=0;
	}
	i++;
    }
    for (i=0, j2=0; j1 != 3 && i<MAXTAIL;)
    {
	if ( l==3 && (((buf[i] ==tail3l[j1]) && (i+BARCODE_GAP>qlen-1 ? 1 : buf[i+BARCODE_GAP]==tail3r[j2]))
		      || ((buf[i] ==seq_nt16_rev_table[tail5l[j1]]) && (i+BARCODE_GAP>qlen-1 ? 1 : buf[i+BARCODE_GAP]==seq_nt16_rev_table[tail3r[j2]])))) {
	    j2++;
	    if (j2==3) break;
	}
	else {
	    i = j2 > 0 ? i-- : i;
	    j2=0;
	}
	i++;
    }
    t2 = j1>j2 ? j1 :j2;
    t2 = t2==0 ? 0 : qlen -i+t2-BARCODE_GAP;
    reverse(buf, qlen);
    int length = t2 > t1 ? t2-t1 : qlen - t1;
    if (length<SKIP_LENGTH) {
	return -1;
    }

    if (t1==0) {
	    for (i=0; i<length; ++i) buf[i] = seq_nt16_str[buf[i]];
	    kputsn((char*)buf, length-1, string);kputs("\n+\n",string);
	    for (i=0; i<length; ++i) buf[i] = qual[i]+33;
	    kputsn((char*)buf, length-1, string);kputc('\n',string);
    }
    else {
	memmove(buf, buf+t1, length);
	buf[length]='\0';
	for (i=0; i<length; ++i) buf[i] = seq_nt16_str[buf[i]];
	kputsn((char*)buf, length-1, string);kputs("\n+\n",string);
	for (i=t1; i<length; ++i) buf[i] = qual[i]+33;
	kputsn((char*)buf, length-1, string);kputc('\n',string);
    }
    return length==qlen ? 0 : 1;
}
int main(int argc, char **argv)
{
    if (argc != 3) return usage();
    samFile *fp = NULL;
    FILE *fq = NULL;
    bam_hdr_t *h;
    bam1_t *b;
    fp = sam_open(argv[1], "r");
    if ( fp == NULL ) {
	fprintf(stderr, "%s : %s\n", argv[1], strerror(errno));
	exit(1);
    }
    if (!strcmp(argv[2],"-")) fq = stdout;
    else fq = fopen(argv[2], "w");
    if ( fq == NULL )  {
	fprintf(stderr, "%s : %s\n", argv[2], strerror(errno));
	exit(1);
    }
    h = sam_hdr_read(fp);
    b = bam_init1();
    kstring_t string = {0, 0, NULL};
    uint64_t n_reads = 0;
    uint64_t n_skip = 0;
    uint64_t n_trim = 0;
    int i;
    uint8_t *buf=NULL;
    uint32_t max_buf=0;
    while ( sam_read1(fp, h, b)>=0 )
    {
	n_reads++;
	int32_t qlen = b->core.l_qseq;
	assert(qlen>0);
	if ( qlen < SKIP_LENGTH ) {
	    n_skip++;
	    continue;
	}
	uint8_t *seq;
	uint8_t *qual = bam_get_qual(b);
	if (max_buf<qlen+1) {
	    max_buf=qlen+1;
	    kroundup32(max_buf);
	    buf = realloc(buf, max_buf);
	    if ( buf==NULL) {
		fprintf(stderr, "Out of memory\n");
		return 1;
	    }
	}
	buf[qlen] = '\0';
	
	string.l = 0;
	kputc('@', &string);
	kputs(bam_get_qname(b), &string);
	kputc('\n', &string);
	seq = bam_get_seq(b);
	for (i=0; i<qlen; i++) buf[i] = bam_seqi(seq, i);
	int ret = trim_tail(buf,qual,qlen,&string);
	if (ret == 1) n_trim++;
	else if ( ret == -1) {
	    n_skip++;
	    continue;
	}
	fputs(string.s,fq);
    }
    free(string.s);
    free(buf);
    bam_destroy1(b);
    bam_hdr_destroy(h);
    sam_close(fp);
    fprintf(stderr,"n_reads : %lld\nn_skip : %lld\nn_trim : %lld\n",n_reads, n_skip, n_trim);
    return 0;
}
