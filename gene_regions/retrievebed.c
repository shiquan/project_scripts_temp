/* retrievebed.c -- retrieve bed region by a gene or transcript list from dbref databases
  INPUT: 1) transcript.list/gene.list 2)dbref database
  OUT: bed like tab seperated file/ stdout
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <assert.h>
#include <errno.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/khash.h>

#define DBREF_PATH  "/Users/shiquan/Documents/02.codes/project_scripts_temp/ncbi_anno_rel104.dbref.gz"

typedef char* string;

KHASH_MAP_INIT_STR(list, string)

typedef kh_list_t listHash_t;

listHash_t *lh = NULL;
char **fn_read = NULL;
int lines = 0;
char *dbref_fn;

KSTREAM_INIT(gzFile, gzread, 8192)

enum line_type {
    is_header,
    is_body,
    is_empty,
    is_comment,
    is_unknown,
};

void init_list(char *fn)
{
    assert(lh == NULL);
    if (fn == NULL)
	return;
    int i;
    khiter_t k;
    int ret;

    fn_read = hts_readlines((const char *)fn, &lines);
    if (lines == 0)
	return;
    lh = kh_init(list);
    for (i=0; i<lines; ++i) {
	char *name = fn_read[i];
	k = kh_get(list, lh, name);
	if (k == kh_end(lh)) {
	    k = kh_put(list, lh, name, &ret);
	}
    }
}
void destroy_list()
{
    int i;
    if (lh)
	kh_destroy(list,lh);
    if (fn_read) {
	for (i=0; i<lines; ++i) 
	    free(fn_read[i]);
	free(fn_read);
    }
}
enum line_type check_line_type(char *line)
{
    if (line == NULL)
	return is_unknown;
    char *s = (char*)line;
    int length;
    length = strlen(s);
    if (length == 0 || length == 1)
	return is_empty;
    if (s[0] == '#')
	return is_comment;
    if (s[0] == 'N' || s[0] == 'X')
	return is_header;
    return is_body;
}


int usage()
{
    fprintf(stderr, "retrievebed nm|gene list_file [dbref database]\n");
    return 1;
}

enum list_type {
    is_trans,
    is_genes,
    is_error,
};
int retrieve_from_dbref(enum list_type list_type);

int main(int argc, char **argv)
{
    if (argc == optind || argc < 2)
	return usage();

    init_list(argv[2]);    
    dbref_fn = argc == 3 ? argv[3] : NULL;
    int ret = 1;
    enum list_type list_type = is_error;
    if (!strcmp(argv[1], "nm"))
	list_type = is_trans;
    else if (!strcmp(argv[1], "gene"))
	list_type = is_genes;
    else
	ret = usage();
    ret = retrieve_from_dbref(list_type);
    //free(dbref_fn);
    destroy_list();    
    return ret;
}

int retrieve_from_dbref(enum list_type list_type)
{
    if (list_type == is_error)
	return 1;
    if (dbref_fn == NULL) dbref_fn = DBREF_PATH;

    gzFile fp = gzopen(dbref_fn, "r");
    if (fp == NULL) {
	fprintf(stderr, "%s: %s\n", dbref_fn, strerror(errno));
	return 1;
    }

    kstream_t *ks;
    int dret;
    ks = ks_init(fp);
    kstring_t *str;
    str = (kstring_t*) malloc(sizeof(kstring_t));
    int body_included = 0;
    while (ks_getuntil(ks, 2, str, &dret) >= 0) {
	enum line_type line_type = check_line_type(str->s);
	if (line_type == is_comment || line_type == is_unknown || line_type == is_empty)
	    continue;
      header:
	body_included = 0;
	if (line_type == is_header) {
	    int nfields = 0;
	    int *splits = ksplit(str, 0, &nfields);
	    assert(nfields >  9);
	    khiter_t k =0;

	    char *name = list_type == is_trans ? str->s + splits[0] : str->s + splits[8];
	    k = kh_get(list, lh, name);
	    if (k != kh_end(lh)) body_included = 1;
	    if (body_included == 1) {
		char *gene_name = strdup(str->s +splits[8]);
		char *trans_name = strdup(str->s + splits[0]);
		char *chrom = strdup(str->s + splits[1]);
		str->l = 0;
		while (ks_getuntil(ks, 2, str, &dret) >= 0) {
		    line_type = check_line_type(str->s);
		    if (line_type == is_body) {
			int *ss = ksplit(str, 0, &nfields);
			char *func = str->s + ss[0];
			int32_t start = atoi(str->s + ss[1]);
			int32_t stop = atoi(str->s+ ss[2]);
			char *exIn = str->s + ss[4];
			fprintf(stdout,"%s\t%d\t%d\t%s\t%s\t%s\t%s\n",chrom,start,stop, gene_name, trans_name, func,exIn);
			free(ss);
		    } else if ( line_type == is_header) {
			goto header;
		    } 
		    str->l = 0;
		}
	    }
	    free(splits);   
	}
	
	str->l = 0;
    }
    free(str);
    gzclose(fp);
    return 0;
}
