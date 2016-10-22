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
#include <htslib/tbx.h>
#include "genepred.h"

typedef char * string;

KHASH_MAP_INIT_STR(list, string)

typedef kh_list_t list_hash_t;

KSTREAM_INIT(gzFile, gzread, 8192)

struct list {
    list_hash_t *hash;
    int lines;
    char **reads;
};

#define LIST_INIT { 0, 0, 0}

struct args {
    const char *refgene_fname;
    const char *format;
    int noheader;
    const char *fast;
    struct list *genes;
    struct list *transcripts;    
} args = {
    .refgene_fname = NULL,
    .format = NULL,
    .fast = NULL,
    .noheader = 0,
    .genes = NULL,
    .transcripts = NULL,
};

struct list *init_list(const char *fn)
{
    if (fn == NULL)
	return NULL;
    int i;
    khiter_t k;
    int ret;
    struct list *list = (struct list*)malloc(sizeof(struct list));
        
    list->reads = hts_readlines(fn, &list->lines);
    if ( list->reads == NULL) {
	fprintf(stderr, "%s : %s\n", fn, strerror(errno));
        free(list);
        return NULL;
    }
    
    if (list->lines == 0)
        goto empty_list;

    list->hash = kh_init(list);
    for ( i=0; i< list->lines; ++i ) {
	char *name = list->reads[i];
        if (name == NULL)
            continue;
        if (*name == '#' || *name == '/')
            continue;

	k = kh_get(list, list->hash, name);
	if (k == kh_end(list->hash))
	    k = kh_put(list, list->hash, name, &ret);	
    }
    
    return list;
    
  empty_list:
    free(list->reads);
    free(list);
    return NULL;
}
void destroy_list(struct list* list)
{
    if (list == NULL)
        return;
    if ( list->reads ) {
        int i;
        for (i=0; i<list->lines; ++i) 
	    free(list->reads[i]);
	free(list->reads);
    }
   
    if ( list->hash )
	kh_destroy(list, list->hash);
    free(list);
}
void destroy_args()
{
    destroy_list(args.genes);
    destroy_list(args.transcripts);
}
int usage()
{
    fprintf(stderr, "retrievebed\n");
    fprintf(stderr, "    -nm transcripts.txt\n");
    fprintf(stderr, "    -gene genes.txt\n");
    fprintf(stderr, "    -fast < one gene or transcript name>\n");
    fprintf(stderr, "    -format [ genepred | refgene | refflat ]\n");
    fprintf(stderr, "    -noheader\n");
    fprintf(stderr, "   [genepred.tsv.gz]\n");
    return 1;
}
int parse_args(int argc, char **argv)
{
    if ( argc == 1)
        return usage();
    --argc, ++argv;
    int i;
    const char *transcripts = 0;
    const char *genes = 0;
    for ( i = 0; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if ( strcmp(a, "-nm") == 0 && transcripts == 0 )
            var = &transcripts;
        else if ( strcmp(a, "-gene") == 0 && genes == 0 )
            var = &genes;
        else if ( strcmp(a, "-format") == 0 && args.format == 0)
            var = &args.format;
        else if ( strcmp(a, "-fast") == 0 && args.fast == 0)
            var = &args.fast;
        if (var != 0 ) {
            if ( i == argc ) {
                fprintf(stderr,"Missing an argument after %s", a);
                return 1;
            }
            *var = argv[i++];
            continue;
        }
        if ( strcmp(a, "-noheader") == 0 ) {
            args.noheader = 1;
            continue;
        }
        if ( args.refgene_fname == 0) {
            args.refgene_fname = a;
            continue;
        }
        fprintf(stderr,"Unknown argument : %s", a);
        return 1;
    }
    if ( args.refgene_fname == 0) {
        args.refgene_fname = getenv("REFGENE");
        if (args.refgene_fname == 0) {
            fprintf(stderr, "No genepred or refgene databases is set.");
            return 1;
        }            
    }

    if ( genes == 0 && transcripts == 0 && args.fast == 0)
        return usage();
    // if -fast specified, ignore gene or transcript list
    if ( args.fast == 0) {
        args.genes = init_list(genes);
        args.transcripts = init_list(transcripts);
    }
    
    if ( args.format == 0 )
        args.format = "genepred";

    if ( strcmp(args.format, "genepred") == 0 ) {
        set_format_genepred();
    } else if ( strcmp(args.format, "refgene") == 0 ) {
        set_format_refgene();
    } else if ( strcmp(args.format, "refflat") == 0 ) {
        set_format_refflat();
    } else {
        fprintf(stderr, "Unknown format, %s", args.format);
        return 1;
    }
    
    return 0;
}

int retrieve_from_dbref()
{
    gzFile fp = gzopen(args.refgene_fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "%s: %s\n", args.refgene_fname, strerror(errno));
	return 1;
    }

    kstream_t *ks;
    int dret = 0;
    ks = ks_init(fp);
    khiter_t k;
    kstring_t string = {0, 0, 0};
    struct genepred line;
    if (args.noheader == 0 ) 
        fprintf(stdout, "#Chrom\tStart(0based)\tEnd(1based)\tStrand\tGene\tTranscript\tExon\tStart loc\tEnd loc\n");
    while (ks_getuntil(ks, 2, &string, &dret) >= 0) {
        if ( string.l == 0 )
            continue;
        if ( string.s[0] == '#' || string.s[0] == '/')
            continue;
        genepred_parser(&string, &line);
        if ( args.fast ) {
            if ( strcmp(line.name1, args.fast) == 0 || strcmp(line.name2, args.fast) == 0)
                generate_dbref_database(&line);
        } else if ( args.genes ) {
            list_hash_t *hash = args.genes->hash;
            k = kh_get(list, hash, line.name2);
            if ( k == kh_end(hash) ) {
                if ( args.transcripts ) {
                    hash = args.transcripts->hash;
                    k = kh_get(list, hash, line.name1);
                    if ( k != kh_end(hash) ) {
                        generate_dbref_database(&line);
                    }
                }
            } else {
                generate_dbref_database(&line);
            }
        } else if ( args.transcripts ) {
            list_hash_t *hash = args.transcripts->hash;
            k = kh_get(list, hash, line.name1);
            if ( k != kh_end(hash) ) {
                generate_dbref_database(&line);
            }            
        }
	string.l = 0;
    }
    ks_destroy(ks);
    free(string.s);
    gzclose(fp);
    return 0;
}
int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) )
        return 1;
    int ret;
    ret = retrieve_from_dbref();
    destroy_args();    
    return ret;
}

