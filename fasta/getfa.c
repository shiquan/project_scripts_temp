/* getfa.c -- this is a fork of samtools faidx. The script is demo for htslib learner.
 * shiquan@link.cuhk.edu.hk
 * getfa in.fa regs
 */

#include<stdio.h>
#include<stdlib.h>
#include<htslib/faidx.h>
int main(int argc, char **argv)
{
    if (argc != 3) {
	fprintf(stderr, "Usage: %s <in.fa> <regions>\n", argv[0]);
	return 1;
    }
    faidx_t *fai = fai_load(argv[1]);\
    if (fai == NULL) {
	fprintf(stderr, "Failed to load fasta file\n");
	return 1;
    }
    int len;
    char *s = fai_fetch(fai, argv[2], &len);
    fprintf(stdout, "%s\n%s\n", argv[2], s);
    fai_destroy(fai);
    return 0;
}
