/* pick predefined columns from BAM/SAM/CRAM file
 * shiquan@link.cuhk.edu.hk
 * Usage: bam_pick -f format_string in.bam > retrieved.tsv
 * format_string is a predefined tags string, like CHROM,POS,SEQ,QUAL,AS,...
 */

#include<stdio.h>
#include<stdlib.h>
#include<htslib/sam.h>
#include<htslib/kstring.h>
struct col1 {
    const char *key; // CHROM, POS, SEQ, ...
    int id;
    int (*setter)(bam_hdr_t *, bam1_t *, struct col*, kstring_t *str);
};
struct cols {
    int m, l;
    struct col1 *a;
};
void cols_clear(struct cols *cols)
{
    int i;
    for (i = 0; i < cols->l; ++i) {
	free(cols->a[i].key);
    }
    free(cols->a);    
}
struct args {
    const char *file_in;
    const char *format_string;
    htsFile *fp;
    bam_hdr_t *h;
    struct cols *cols;
};
struct args args = {
    .file_in = NULL,
    .format_string = NULL,
    .fp = NULL,
    .h = NULL,
    .cols = NULL,
};
void args_destroy()
{
    hts_close(args.fp);
    bam_hdr_destroy(args.h);
    cols_clear(args.cols);
}
int setter_chrom(bam_hdr_t *, bam1_t *, struct col*, kstring_t *str);
int setter_pos(bam_hdr_t *, bam1_t *, struct col*, kstring_t *str);
int setter_mapq(bam_hdr_t *, bam1_t *, struct col*, kstring_t *str);
int setter_queryname(bam_hdr_t *, bam1_t *, struct col*, kstring_t *str);
int setter_flag(bam_hdr_t *, bam1_t *, struct col*, kstring_t *str);
int setter_cigar(bam_hdr_t *, bam1_t *, struct col*, kstring_t *str);
int setter_mchrom(bam_hdr_t *, bam1_t *, struct col*, kstring_t *str);
int setter_mpos(bam_hdr_t *, bam1_t *, struct col*, kstring_t *str);
int setter_seqs(bam_hdr_t *, bam1_t *, struct col*, kstring_t *str);
int setter_quals(bam_hdr_t *, bam1_t *, struct col*, kstring_t *str);
int setter_tag(bam_hdr_t *, bam1_t *, struct col*, kstring_t *str);

void format_string_init(const char *_string, struct cols *cols)
{
    char *string = (char*)strdup(_string);
    char *ss = string, *se;

    // CHROM,POS,QUAL,SEQ,...,TAG
    while(*ss) {
	while (*se && *se != ',') se++;
	char *key =(char*)strndup(ss, se-ss);
	if (cols->m == cols->l ) {
	    cols->m = cols->m == 0 ? 2 : cols->m << 1;
	    cols->a = (struct col1*)realloc(cols->a, cols->m*sizeof(struct col1));
	}
	struct col1 *c = &cols->a[cols->l++];
	c->key = key;
	c->id = -1;
	if ( strcmp(key, "CHROM") == 0 ) {
	    c->setter = setter_chrom;
	} else if ( strcmp(key, "POS") == 0 ) {
	} else if ( strcmp(key, "MAPQ") == 0 ) {
	} else if ( strcmp(key, "QNAME") == 0) {
	}
	ss = se + 1;
    }
    
}
void args_praser(int ac, char**av)
{
}
// bam_pick_core convert the bam structure to predefined cols into str
void bam_pick_core(bam_hdr_t *h, bam1_t *line, int n_cols, struct col *cols, kstring_t *str)
{
    int i;
    for (i = 0; i < n_cols; ++i)
	cols[i].setter(h, line, &cols[i], str);
}
int main(inta argc, char **argv)
{
    args_praser(argc, argv);
    bam1_t *b = bam_init1();
    int r;
    while ( (r = sam_read1(args.fp, args.h, b) ) >= 0) {
	bam_pick_core(args.h, b, );
    }

    bam_destroy1(b);
    args_destroy();
    return 0;
}
