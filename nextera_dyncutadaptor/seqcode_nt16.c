#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

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

int main(int argc, char * argv[]) 
{
	if (argc < 2) {
		fprintf(stderr, "Usage: %s <sequence>\n", argv[0]);
		return 0;
	}
	char * seq;
	seq = strdup(argv[1]);
	int len, i;
	len = strlen(seq);
	uint8_t * compl;
	compl = calloc(len, sizeof(uint8_t));
	fputs("seqen: ", stdout); fputs(seq, stdout); putchar('\n'); fputs("code: ", stdout);
	for (i = 0; i < len; ++i) {
		compl[i] = bam_nt16_table[seq[i]];
		fprintf(stdout, "%d, ",compl[i]);
	}
	putchar('\n');
	for (i = 0; i < len; ++i) {
		compl[i] = seq_comp_table[compl[i]];
		seq[i] = bam_nt16_rev_table[compl[i]];
	}
	fputs("compl: ", stdout); fputs(seq, stdout); putchar('\n'); fputs("code: ", stdout);
	for (i = 0; i < len; ++i) fprintf(stdout, "%d, ",compl[i]);
	for (i = 0; i < len>>1; ++ i) {
		uint8_t t = compl[i] ;
		compl[i] = compl[len - i -1];
		compl[len - i - 1] = t;
		seq[i] = bam_nt16_rev_table[compl[i]];
		seq[len - i - 1] = bam_nt16_rev_table[compl[len - i -1]];
	}
	putchar('\n');
	fputs("trans: ", stdout); fputs(seq, stdout); putchar('\n'); fputs("code: ", stdout);
	for (i = 0; i < len; ++i) {
		fprintf(stdout, "%d, ",compl[i]);
		compl[i] = seq_comp_table[compl[i]];
		seq[i]  = bam_nt16_rev_table[compl[i]];
	}
	putchar('\n');
	fputs("reves: ", stdout); fputs(seq, stdout); putchar('\n'); fputs("code: ", stdout);
	for (i = 0; i < len; ++i) fprintf(stdout, "%d, ",compl[i]);
	putchar('\n');
	free(seq);
	free(compl);
	return 1;
}

