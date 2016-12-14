#include <stdio.h>
#include <stdlib.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <zlib.h>
#include <unistd.h>
#include "utils.h"
KSEQ_INIT(gzFile, gzread)

#define KSTRING_INIT { 0, 0, 0}
struct args {
    gzFile fp;
    FILE *out;
    const char *output_fname;
    kseq_t *seq;
};

struct args args = {
    .fp = NULL,
    .out = NULL,
    .output_fname = 0,
    .seq = NULL,
};
struct buffer {
    kstring_t name;
    kstring_t seq;
};
struct kstring_buffer {
    int m, n;
    struct buffer *buffer;
} buffer = {
    .m = 0,
    .n = 0,
    .buffer = 0,
};

void args_destroy()
{
    fclose(args.out);
    kseq_destroy(args.seq);
}

int generate_entry(kstring_t *name1, kstring_t *name2, kstring_t *string1, kstring_t *string2)
{
    assert(string1->l == string2->l);
    kstring_t string = KSTRING_INIT;
    int i;
    for ( i = 0; i < string1->l; ++i ) {
	char base = 'N';
	char A1 = string1->s[i];
	char A2 = string2->s[i];
	if (A1 == 'A' || A1 == 'a') {
	    switch (A2) {
		case 'C':
		case 'c':
		    base = 'M';
		    break;
		case 'G':
		case 'g':
		    base = 'R';
		    break;
		case 'T':
		case 't':
		    base = 'W';
		    break;
		default:
		    base = 'A';
		    break;
	    }
	} else if (A1 == 'C' || A1 == 'c') {
	    switch (A2) {
		case 'A':
		case 'a':
		    base = 'M';
		    break;
		case 'G':
		case 'g':
		    base = 'S';
		    break;
		case 'T':
		case 't':
		    base = 'Y';
		    break;
		default:
		    base = 'C';
		    break;
	    }
	} else if (A1 == 'G' || A1 == 'g') {
	    switch (A2) {
		case 'A':
		case 'a':
		    base = 'R';
		    break;
		case 'C':
		case 'c':
		    base = 'S';
		    break;
		case 'T':
		case 't':
		    base = 'K';
		    break;
		default:
		    base = 'G';
		    break;
	    }
	} else if (A1 == 'T' || A1 == 't') {
	    switch (A2) {
		case 'A':
		case 'a':
		    base = 'W';
		    break;
		case 'C':
		case 'c':
		    base = 'Y';
		    break;
		case 'G':
		case 'g':
		    base = 'K';
		    break;
		default:
		    base = 'T';
		    break;
	    }
	} else {
	    switch (A2) {
		case 'A':
		case 'a':
		    base = 'A';
		    break;
		case 'C':
		case 'c':
		    base = 'C';
			break;
		case 'G':
		case 'g':
		    base = 'G';
		    break;
		case 'T':
		case 't':
		    base = 'T';
		    break;
		default:
			base = 'N';
			break;
	    }
	}
	if ( base != 'N') kputc(base, &string);
    }
    kstring_t name = KSTRING_INIT;
    kputs(name1->s, &name);
    kputs(",", &name );
    kputs(name2->s, &name);
    fprintf(args.out, ">%s\n%s\n", name.s, string.s);
    free(name.s);
    free(string.s);
    return 0;
}
int buffer_push(kstring_t *name, kstring_t *seq)
{
    if ( buffer.m == buffer.n ) {
	buffer.m = buffer.m == 0 ? 4 : buffer.m<<1;
	buffer.buffer = (struct buffer*)realloc(buffer.buffer, buffer.m * sizeof(struct buffer));
    }
    memset(&buffer.buffer[buffer.n].name, 0, sizeof(kstring_t));
    memset(&buffer.buffer[buffer.n].seq, 0, sizeof(kstring_t));
    kputs(name->s, &buffer.buffer[buffer.n].name);
    kputs(seq->s, &buffer.buffer[buffer.n].seq);
    buffer.n++;
    return 0;
}
void buffer_clear()
{
    int i;
    for (i = 0; i < buffer.n; ++i ) {
	free(buffer.buffer[i].name.s);
	free(buffer.buffer[i].seq.s);	
    }
    free(buffer.buffer);
    buffer.m = buffer.n = 0;
    buffer.buffer = NULL;	
}
int export_database()
{
    int l;
    int i;
    while ( (l = kseq_read(args.seq)) >= 0) {
	for (i = 0; i < buffer.n; ++i )
	    generate_entry(&args.seq->name, &buffer.buffer[i].name, &args.seq->seq, &buffer.buffer[i].seq);
	buffer_push(&args.seq->name, &args.seq->seq);
    }
    buffer_clear();
    return 0;
}
int parse_args(int argc, char **argv)
{
    const char *input_fname = 0;
    int i;
    for (i = 0; i < argc;) {
	const char *a = argv[i++];
	if ( strcmp(a, "-h") == 0 )
	    error("database_construct -o out.fa in.fa");
	const char ** var = 0;
	if ( strcmp(a, "-o") == 0 && args.output_fname == 0) 
	    var = &args.output_fname;

	if (var != 0) {
	    if ( i == argc )
		error("Missing an argument after %s", a);
	    *var = argv[i++];
	    continue;
	}
	if (input_fname == 0) {
	    input_fname = a;
	    continue;
	}
	error("Unknow argument : %s.", a);
    }
    if ( input_fname == 0) {
        if ( !isatty(fileno(stdin)) ) {
            args.fp = gzdopen(fileno(stdin), "r");
        } else {
            error("database_construct -o out.fa in.fa");
        }
    } else {
        args.fp = gzopen(input_fname, "r");
    }

    if (args.fp == NULL)
        error("Failed to open %s : %s.", input_fname == 0 ? "-" : input_fname, strerror(errno));
    args.seq = kseq_init(args.fp);
    if (args.output_fname == 0)
        args.out = stdout;
    else
        args.out = fopen(args.output_fname, "w");
    return 0;
}
int main (int argc, char **argv)
{
    parse_args(argc, argv);
    export_database();
    return 0;
}
