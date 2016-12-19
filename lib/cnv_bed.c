#include "utils.h"
#include "number.h"
#include "cnv_bed.h"
#include <htslib/khash.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <zlib.h>

KSTREAM_INIT(gzFile, gzread, 8192)

KHASH_MAP_INIT_STR(name, int)
typedef kh_name_t namehash_type;

int flag_inconsis(int flag)
{
    int i, j;
    int c = 0;
    for (i = 0, j = 0; i < 4; ++i ) {
        c += (flag>>i) & 1;
    }
    return c < 3 ? 0 : 1;
}
int combine_flag(int flag1, int flag2)
{
    if ( (flag1 & 0x3) && (flag2 & 0x3) )
        return CNV_DEL_HOM;
    if ( (flag1 & 0x4) && (flag2 & 0x4) )
        return CNV_DUP_HOM;
    return flag1 | flag2;            
}

struct cnv_spec *cnv_spec_init()
{
    struct cnv_spec *spec = (struct cnv_spec*)malloc(sizeof(struct cnv_spec));
    spec->fp = NULL;
    spec->ks = NULL;
    spec->string.l = spec->string.m = 0;
    spec->string.s = NULL;
    spec->chrhash = kh_init(name);
    spec->samhash = kh_init(name);
    spec->n_samples = spec->m_samples = 0;
    spec->samples = NULL;
    spec->n_chrom = spec->m_chrom = 0;
    spec->chrom = NULL;
    spec->min_length = spec->max_length = 0;
    return spec;
}

void cnv_spec_destroy(struct cnv_spec *spec)
{
    int i;

    gzclose(spec->fp);
    ks_destroy(spec->ks);
    if ( spec->string.m )
        free(spec->string.s);
    
    kh_destroy(name, (namehash_type*)spec->chrhash);
    kh_destroy(name, (namehash_type*)spec->samhash);
    for ( i = 0; i < spec->n_samples; ++i )
        free(spec->samples[i]);
    if ( spec->samples )
        free(spec->samples);
    for ( i = 0; i < spec->n_chrom; ++i )
        free(spec->chrom[i]);
    if ( spec->chrom )
        free(spec->chrom);
    free(spec);
    spec = NULL;
}

int parse_type(char *type)
{
    if ( !strcmp(type, "N/N") || !strcmp(type, "<REF>/<REF>")  || !strcmp(type, "N") || !strcmp(type, "<REF>") )
        return CNV_REF;
    if ( !strcmp(type, "N/<DEL>") || !strcmp(type, "<DEL>/N") || !strcmp(type, "<REF>/<DEL>") )
        return CNV_DEL_HET;
    if ( !strcmp(type, "N/<DUP>") || !strcmp(type, "<DUP>/N") || !strcmp(type, "<REF>/<DUP>") )
        return CNV_DUP_HET;
    if ( !strcmp(type, "<DEL>/<DUP>") || !strcmp(type, "<DUP>/<DEL>"))
        return CNV_DEL_DUP;
    if ( !strcmp(type, "<DUP>/<DUP>"))
        return CNV_DUP_HOM;
    if ( !strcmp(type, "<DEL>/<DEL>"))
        return CNV_DEL_HOM;
    return CNV_REF;
}

const char *explain_type(int flag)
{
    static const char * const cnv_types[] = {
        "<REF>/<REF>",
        "<REF>/<DEL>",
        "<DEL>/<DEL>",
        "<REF>/<DUP>",
        "<DUP>/<DUP>",
        "<DEL>/<DUP>",
        "<UNKNOWN>",
    };        
    flag &= CNV_MASK;
    if ( !flag )
        return cnv_types[0];
    if ( flag == CNV_DEL_HET ) 
        return cnv_types[1];
    if ( flag == CNV_DEL_HOM )
        return cnv_types[2];                
    if ( flag == CNV_DEL_DUP )
        return cnv_types[5];
    if ( flag == CNV_DUP_HET)
        return cnv_types[3];
    if (flag == CNV_DUP_HOM)
        return cnv_types[4];

    return cnv_types[6];
}

struct cnv_bed_format {
    int chrom_col;
    int start_col;
    int end_col;
    int type_col;
    int sample_col;
} format = {
    .chrom_col = 0,
    .start_col = 1,
    .end_col = 2,
    .type_col = 3,
    .sample_col = 4,
};

// accept bed format, Chrom\tStart(1-based)\tEnd\tType\t[Sample]
static int parse_line(struct cnv_spec *spec, kstring_t *string, struct cnv_bed *line)
{
    if ( string->l == 0 )
        return 1;
    if (string->s[0] == '#')
        return 1;
    
    char *raw_line = strdup(string->s);
    int nfields = 0;
    int *splits = ksplit(string, 0, &nfields);
    if (splits == NULL )
        return 1;
    if ( nfields < 4)
        goto skip_line;
    
    char *chrom = string->s+splits[format.chrom_col];
    char *start = string->s+splits[format.start_col];
    char *end = string->s+splits[format.end_col];
    char *type = string->s+splits[format.type_col];
    char *sample = nfields > format.sample_col ? string->s+splits[format.sample_col] : NULL;

    // Only check DUP/DEL regions, skip other types
    if (*start == '.' || *end == '.')
        goto skip_line;

    line->next = NULL;
    khint_t k;
    k = kh_get(name, (namehash_type*)spec->chrhash, chrom);
    if ( k == kh_end((namehash_type*)spec->chrhash) ) {
        if ( spec->n_chrom == spec->m_chrom) {
            spec->m_chrom += 4;
            spec->chrom = (char**)realloc(spec->chrom, spec->m_chrom*sizeof(char*));
        }
        int ret;
        spec->chrom[spec->n_chrom] = strdup(chrom);
        k = kh_put(name, (namehash_type*)spec->chrhash, spec->chrom[spec->n_chrom], &ret);
        kh_val((namehash_type*)spec->chrhash, k) = spec->n_chrom;
        spec->n_chrom++;
    }
    // init chromosome
    line->id = kh_val((namehash_type*)spec->chrhash, k);
    
    if ( check_num_likely(start) ) {
        line->start = str2int(start);
    } else {
        fprintf(stderr, "Failed to parse line %s.\n", raw_line);
        goto skip_line;
    }

    if ( check_num_likely(end) ) {
        line->end = str2int(end);
    } else {
        fprintf(stderr, "Failed to parse line %s.\n", raw_line);
        goto skip_line;
    }

    if ( spec->min_length && line->end - line->start < spec->min_length)
        goto skip_line;
    if ( spec->max_length && line->end - line->start > spec->max_length)
        goto skip_line;
    line->flag = parse_type(type);
    if ( check_malformed(line->flag) || flag_inconsis(line->flag) ) {
        fprintf(stderr,"Failed to parse type %s, %d\n", type, line->flag);
        goto skip_line;
    }

    // If sample column is set, check the sample names then.
    if ( sample ) {
        if ( spec->n_samples == spec->m_samples ) {
            spec->m_samples = spec->m_samples == 0 ? 2 : spec->m_samples + 10;
            spec->samples = (char**)realloc(spec->samples, spec->m_samples*sizeof(char*));
        }
        khint_t k = kh_get(name, (namehash_type*)spec->samhash, sample);        
        if ( k == kh_end((namehash_type*)spec->samhash) ) {
            int ret;
            spec->samples[spec->n_samples] = strdup(sample);
            k = kh_put(name, (namehash_type*)spec->samhash, spec->samples[spec->n_samples], &ret);
            kh_val((namehash_type*)spec->samhash, k) = spec->n_samples;
            spec->n_samples++;
        }
    }
    free(raw_line);
    return 0;
    
  skip_line:    
    free(raw_line);
    free(splits);
    return 1;
}

int cnv_load_fname(struct cnv_spec *spec, const char *fname)
{
    spec->fp = strcmp(fname, "-") == 0 ? gzdopen(fileno(stdin), "r") : gzopen(fname, "r");
    if (spec->fp == NULL)
        error("%s : %s.", fname, strerror(errno));
    spec->fname = fname;
    spec->ks = ks_init(spec->fp);    
    return 0;    
}

// -1 for end, 0 for read success
int cnv_read(struct cnv_spec *spec, struct cnv_bed *line)
{
    if (spec->fp == NULL)
        return -1;
    int ret;
    while ( ks_getuntil(spec->ks, 2, &spec->string, &ret) >= 0) {
        if ( spec->string.l == 0 )
            continue;
        if ( spec->string.s[0] == '#')
            continue;
        if ( parse_line(spec, &spec->string, line) )
            continue;
        return 0;
    }    
    return -1;
}
