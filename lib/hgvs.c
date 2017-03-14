#include "utils.h"
#include "hgvs.h"
#include "genepred.h"
#include "number.h"
#include "sort_list.h"
#include "htslib/kstring.h"
#include <string.h>
#include <regex.h>

struct hgvs_spec {
    struct genepred_spec *data;
    struct hgvs_des des;
    // To parse hgvs nomen.
    regex_t *exp;
} spec;

void hgvs_core_clear(struct hgvs_core *core)
{
    if ( core->name.name1 != NULL)
        free(core->name.name1);
    if ( core->name.name2 != NULL)
        free(core->name.name2);
    memset(core, 0, sizeof(struct hgvs_core));
}
void hgvs_des_clear(struct hgvs_des *des)
{
    int i;
    for ( i = 0; i < des->l; ++i) {
        hgvs_core_clear(&des->a[i]);
    }
    free(des->a);
    if ( des->chrom != NULL )
        free(des->chrom);
    if ( des->ref != NULL)
        free(des->ref);
    if ( des->alt != NULL)
        free(des->alt);
    memset(des, 0, sizeof(struct hgvs_des));
}
void hgvs_spec_destroy()
{
    genepred_spec_destroy(spec.data);
    hgvs_des_clear(&spec.des);
    free(spec.exp);
}
int init_hgvs_spec(const char *fname, const char *fasta)
{
    memset(&spec, 0, sizeof(struct hgvs_spec));
    spec.data = genepred_spec_init();
    if ( genepred_load_data(spec.data, fname) == NULL )
        return 1;
    if ( genepred_load_fasta(spec.data, fasta) == NULL )
        return 1;
    spec.exp = (regex_t*)malloc(sizeof(regex_t));
    int rv;
    rv = regcomp(spec.exp, "(.*):[gcn][.](.*)", REG_EXTENDED);
    if ( rv != 0 ) {
        error("regcomp failed with %d.", rv);
    }    
    return 0;
}
enum nametype {
    name_is_unknown = -1,
    name_is_chromosome,
    name_is_coding_rna,
    name_is_noncoding_rna,
};
static int convert_loc2position(struct genepred_line *line, int location, int offset)
{
    int i;
    if ( line->strand == '+') {
        for ( i = 0; i < line->exon_count; ++i ) {
            if ( location >= line->loc[BLOCK_START][i] && location <= line->loc[BLOCK_END][i])
                break;            
        }
        
        if ( location > line->loc[BLOCK_END][i])
            return 0;
        if ( offset != 0 ) {
            if ( line->exons[BLOCK_START][i] != location && line->exons[BLOCK_END][i] != location ) {
                error("Inconsist between offset and exon edge. exon location : %d, offset: %d.", location, offset);
            }
        }
        return line->exons[BLOCK_START][i] + location - line->loc[BLOCK_START][i] + offset;
    } else {
        for ( i = 0; i < line->exon_count; ++i ) {
            if ( location >= line->loc[BLOCK_END][i] && location <= line->loc[BLOCK_START][i])
                break;            
        }
        if ( offset != 0 ) {
            if ( line->exons[BLOCK_START][i] != location && line->exons[BLOCK_END][i] != location ) {
                error("Inconsist between offset and exon edge. exon location : %d, offset: %d.", location, offset);
            }
        }
        return line->exons[BLOCK_START][i] + line->loc[BLOCK_START][i] - location + offset;                              
    }
    return 0;
}
// ss point the start of string, se point the end of the convert string.
static int parse_position(char *ss, char *se, struct genepred_line *line)
{
    // Parse the position string. The possible position may be look like 123, *123, -123, 123+12, 123-12, etc.        
    int position = 0;
    int offset = 0;
    // Assume all the position located in the cds region first. For intron region, offset should NOT be 0, and
    // pos_type is the type of nearby region.        
    enum func_region_type pos_type = func_region_cds;
    char *s1 = ss;
    char *s2;
    if ( *s1 == '*' ) {
        pos_type = func_region_utr3;
        ++s1;
    } else if ( *s1 == '-') {
        pos_type = func_region_utr5;
        ++s1;
    }
    // For 
    int i;
    for ( s2 = s1, i = 0; *s2 >= 48 && *s2 <= 57; ++s2, ++i );
    position = str2int_l(s1, i);        
    if ( s2 != se ) {
        for ( i = 0, s1 = s2; s2 != se; ++s2, ++i);
        if ( check_num_likely_l(s1, i) ) {
            offset = str2int_l(s1, i);
        } else {
            error("Failed to parse offset. %s.", s1);
        }
    }

    // Convert functional location to gene location.
    if ( pos_type == func_region_cds ) {
        position += line->forward_length;
        if ( line->strand == '-' )
            position = line->reference_length - position + 1;
    } else if ( pos_type == func_region_utr5 ) {
        if ( line->strand == '+' ) {
            position = line->forward_length - position + 1;
        } else {
            position = line->reference_length - ( line->backward_length - position + 1 ) + 1;
        }
    } else if ( pos_type == func_region_utr3 ) {
        if ( line->strand == '+' ) {                    
            position = line->reference_length - (line->backward_length - position + 1 ) + 1;
        } else {
            position = line->reference_length - (line->forward_length - position + 1 ) + 1;                    
        }
    } else {
        error("Impossible situation.");
    }
    return convert_loc2position(line, position, offset);
}
// Standard HGVS name should be NM_0001.2:c.123A>G; tolerant format could be NM_0001:c.123A>G (no version number);
// Check string format could be parsed and convert cds position to genome position.
// This code is fragile, improve me.
int check_hgvs_name(const char *name)
{
    // Only allow match once.
    regmatch_t matches[1]; 
    if ( regexec(spec.exp, name, 1, matches, 0) ) {
        int i, l;
        l = strlen(name);
        fprintf(stderr, "%s\n", name);
        for (i = 0; i < matches[0].rm_so; ++i ) {
            fprintf(stderr, ANSI_COLOR_RED "%c" ANSI_COLOR_RESET, '^');
        }
        for (i = matches[0].rm_so; i < matches[0].rm_eo; ++i ) {
            fprintf(stderr, ANSI_COLOR_GREEN "%c" ANSI_COLOR_RESET, '~');
        }
        for (i = matches[0].rm_eo; i < l; ++i) {
            fprintf(stderr, ANSI_COLOR_RED "%c" ANSI_COLOR_RESET, '^');
        }
        fprintf(stderr, "\n");
        return 1;
    }
    return 0;
}
int parse_hgvs_name(const char *name)
{
    if ( check_hgvs_name(name) )
        return 1;

    // If format test success.    
    struct hgvs_des *des = &spec.des;
    struct genepred_spec *data = spec.data;
    // Check the chromosome or transcript name.
    enum nametype type = name_is_unknown;
    char *ss = (char*)name;
    char *safe_lock = ss;
    char *se = ss;
    char *se1 = ss;
    int i;
    while ( safe_lock && *safe_lock )
        safe_lock++;
    
    for ( ; se != safe_lock && *se != ':'; ++se );    
    for ( i = 0; *se1 != '(' && se1 != se; ++se1, ++i);
    
    kstring_t string = { 0, 0, 0};
    kputsn(ss, i, &string);   
    se++;
    if ( *se == 'g' ) {
        // Assume genome locations.
        type = name_is_chromosome;
    } else if ( *se == 'c' ) {
        type = name_is_coding_rna;
    } else if ( *se == 'n' ) {
        type = name_is_noncoding_rna;
    } else {
        error("Undefined variant type. %c.", se[1]);
    }
    // Check the start position.
    for ( ; se != safe_lock && *se != '.'; ++se);
    se1 = ++se;

    // The possible character of position should be numberic and :
    // -  alias UTR5 region
    // *  alias UTR3 or intron region
    // +  alias intron
    for ( i = 0; se != safe_lock && (se[0] == '*' || se[0] == '-' || se[0] == '+' || (se[0] >= 48 && se[0] <= 57)); ++se, ++i);
    if ( se == se1 )
        error("No position found.");

    if ( type == name_is_chromosome ) {
        des->chrom = strdup(string.s);
        // For genome locations, only number is accept.
        if ( check_num_likely_l(se1, i)) {
            des->start = str2int_l(se1, i);
        } else {
            error("Genome position not readable. %s, %d.", se1, i);
        }
        if ( se[0] == '_') {
            se1 = ++se;
            for ( i = 0; se[0] == '?' || se[0] == '*' || se[0] == '-' || se[0] == '+' || (se[0] >= 48 && se[0] <= 57); ++se, ++i);
            if ( se == se1 )
                error("No position found.");            
            if ( check_num_likely_l(se1, i) == 0 ) {
                des->end = str2int_l(se1, i);
            } else {
                error("Genome position not readable. %s.", se1);
            }            
        } else {
            des->end = des->start;
        }
        
    } else {
        // Here only retrieve one record for each transcript. For some reasons, like alternative locus, one transcript
        // may align to several genome regions, that usually mislead researchers one gene will be expressed by different
        // genome locus in one cell which it is not true. For these multi records, I will give a warning message.
        struct genepred_line *line = genepred_retrieve_trans(data, string.s);
        if ( line == NULL )
            error("No transcript found in data. %s.", string.s);        
        parse_line_locs(line);        
        des->start = parse_position(se1, se, line);
        des->chrom = strdup(line->chrom);
        // Check the end position.    
        if ( se[0] == '_') {
            se1 = ++se;
            for ( i = 0; se[0] == '?' || se[0] == '*' || se[0] == '-' || se[0] == '+' || (se[0] >= 48 && se[0] <= 57); ++se, ++i);
            if ( se == se1 )
                error("No position found.");
            des->end = parse_position(se1, se, line);
            list_lite_delete(&line, genepred_line_destroy);
        } else {
            des->end = des->start;
        }

        if ( line->next != NULL ) {
            warnings("Multiple transcript hits. Only random pick one record. %s.", line->name1);
        }
        list_lite_delete(&line, genepred_line_destroy);
    }

    // Check the variant type.
    for ( se1 = se, i = 0; *se == 'A' || *se == 'a' || *se == 'C' || *se == 'c' || *se == 'G' || *se == 'g' || *se == 'T' || *se == 't'; ++se, ++i);
    
    // Deletion, insertion
    if ( se == se1 ) {        
        if ( se[0] == 'd' && se[1] == 'e' && se[2] == 'l' ) {
            des->alt_length = 0;
            des->alt = 0;
            des->ref_length = des->end - des->start + 1;
            // Ref sequence will be allocated.
            se += 3;
            for ( se1 = se, i = 0; *se == 'A' || *se == 'a' || *se == 'C' || *se == 'c' || *se == 'G' || *se == 'g' || *se == 'T' || *se == 't'; ++se, ++i);
            if ( i > 0 ) {
                if ( i != des->ref_length) {
                    error("Inconsistand deletion length. %s.", ss);
                }
                des->ref = strndup(se1, i);
            }
        } else if ( se[0] == 'i' && se[1] == 'n' && se[2] == 's' ) {
            des->ref_length = 0;
            des->ref = 0;
            se += 3;
            for ( se1 = se, i = 0; *se == 'A' || *se == 'a' || *se == 'C' || *se == 'c' || *se == 'G' || *se == 'g' || *se == 'T' || *se == 't'; ++se, ++i);
            if ( i > 0 ) {
                des->alt_length = i;
                des->alt= strndup(se1, i);
            } else {
                des->alt_length = 0;
                des->alt = 0;
            }            
        } else {
            error("Failed to parse the variant type. %s.", se);
        }
    } else {
        des->ref_length = i;
        des->ref = strndup(se1, i);
        if ( des->end != 0 && des->start + i -1 != des->end ) {
            error("Inconsistant position and sequences. %s.", ss);
        } else {
            des->end = des->start + i - 1;
        }
        // Standard variants.
        if ( *se == '>' ) {
            for ( se1 = ++se, i = 0; *se == 'A' || *se == 'a' || *se == 'C' || *se == 'c' || *se == 'G' || *se == 'g' || *se == 'T' || *se == 't'; ++se, ++i);
            des->alt_length = i;
            des->alt = strndup(se1, i);
        } else {
            error("Failed to parse the variant type. %s.", se);
        }
    }
    return 0;
}
// return the exon id, for intergenic return -1.
int find_the_block(struct genepred_line *line, int *blk_start, int *blk_end, int pos)
{
    *blk_start = -1;
    *blk_end = -1;
    int i;
    for ( i = 0; i < line->exon_count; ++i ) {
        int start = line->exons[BLOCK_START][i];
        int end = line->exons[BLOCK_END][i];
        if ( pos < start) {
            *blk_end = i;
            break;
        } else {
            *blk_start = i;
            if ( pos <= end ) {
                *blk_end = i;
                break;
            }           
        }
    }
    // Always return 0 ? if out of range return 1??
    return 0;
}
int find_locate(struct genepred_line *line, int *pos, int *offset, int start)
{
    int block1 = 0;
    int block2 = 0;
    if ( find_the_block(line, &block1, &block2, start ) )
        return 1;
    if ( block1 == block2 ) {
        if ( line->strand == '+' ) {
            *pos = line->loc[BLOCK_START][block1] + start - line->exons[BLOCK_START][block1];
        } else {
            *pos = line->loc[BLOCK_START][block1] - (start - line->exons[BLOCK_START][block1]);
        }
        *offset = 0;
    } else {
        int upstream =  start - line->exons[BLOCK_END][block1];
        int downstream = line->exons[BLOCK_START][block2];
        if ( upstream > downstream ) {
            *pos = line->loc[BLOCK_START][block2];
            *offset = line->strand == '+' ? -downstream : downstream;
        } else {
            *pos = line->loc[BLOCK_END][block1];
            *offset = line->strand == '+' ? upstream : -upstream;
        }        
    }
    return 0;
}
// return 0 on success, 1 on out of range.
int generate_hgvs_core(struct genepred_line *line, struct hgvs_core *core, int start, int end)
{
    int i;
    int blk_start = 0, blk_end = 0;
    struct hgvs_name *name = &core->name;
    if ( line->loc_parsed == 0 )
        parse_line_locs(line);

    if ( find_locate(line, &name->pos, &name->offset, start) )
        return 1;
    // Locate end. For most cases variants are snps, start == end.
    if ( end != start ) {
        find_locate(line, &name->end_pos, &name->end_offset, end);
    }
    debug_print("reference length : %d, forward length : %d, backward length : %d.", line->reference_length, line->forward_length, line->backward_length);
    if ( line->cdsstart == line->cdsend ) {
        core->type.func = func_region_noncoding;
    } else {
        if ( line->strand == '+' ) {
            if ( name->pos <= line->forward_length ) {
                core->type.func = func_region_utr5;
                name->pos = line->forward_length - name->pos + 1;
            } else if ( line->reference_length - name->pos <= line->backward_length ) {
                core->type.func = func_region_utr3;
                name->pos = line->backward_length - (line->reference_length - name->pos);
            } else {
                core->type.func = func_region_cds;
                name->pos = name->pos - line->backward_length;
            }
        } else {
            if ( name->pos <= line->backward_length ) {
                core->type.func = func_region_utr5;
                name->pos = line->forward_length - name->pos + 1;                
            } else if ( line->reference_length - name->pos <= line->forward_length ) {
                core->type.func = func_region_utr3;
                name->pos = line->backward_length - (line->reference_length - name->pos);                
            } else {
                core->type.func = func_region_cds;
                name->pos = name->pos - line->backward_length;
            }
        }
    }
    name->name1 = strdup(line->name1);
    name->name2 = strdup(line->name2);
    return 0;
}
// Fill all possible HGVS names for this variant.
int fill_hgvs_name()
{
    // Check the position inited.
    struct hgvs_des *des = &spec.des;
    if ( des->start == 0 )
        error("Variant position is not inited.");

    struct genepred_line *line = genepred_retrieve_region(spec.data, des->chrom, des->start-1, des->end);

    for ( ;; ) {
        if ( line == NULL )
            break;
        
        if ( des->l == des->m ) {
            des->m += 2;
            des->a = (struct hgvs_core*)realloc(des->a, des->m*sizeof(struct hgvs_core));
            memset(&des->a[des->l], 0, sizeof(struct genepred_line));
        }

        hgvs_core_clear(&des->a[des->l]);
        if ( generate_hgvs_core(line, &des->a[des->l], des->start, des->end) == 0 )
            des->l++;
        generate_dbref_database(line);
        struct genepred_line * temp = line;
        line = line->next;
        genepred_line_destroy(temp);
    }
    return 0;
}

int print_hgvs_core()
{
    int i;
    struct hgvs_des *des = &spec.des;
    
    for ( i = 0; i < des->l; ++i ) {
        struct hgvs_core *core = &des->a[i];
        if ( core->name.name1 == NULL )
            error("No transcript name found.");
        printf("%s:", core->name.name1);
        if ( core->type.func == func_region_noncoding ) {
            printf("n.");
        } else if ( core->type.func == func_region_cds ) {
            printf("c.");            
        } else if ( core->type.func == func_region_utr5 ) {
            printf("c.-");
        } else if ( core->type.func == func_region_utr3 ) {
            printf("c.*");
        }
        printf("%d", core->name.pos);
        if ( core->name.offset  > 0 ) {
            printf("+%d", core->name.offset);
        } else if ( core->name.offset < 0 ) {
            printf("%d", core->name.offset);
        }
        printf("\n");
    }
    return 0;
}

#ifdef HGVS_MAIN
int usage()
{
    fprintf(stderr,
            "Usage: hgvs_converter NM001:c.123A>T\n"
        );
    return 1;
}
int parse_args(int ac, char **av)
{
    if ( ac == 0 )
        return usage();
    const char *data_fname = 0;
    const char *fasta = 0;
    const char *name = 0;
    int i;
    for ( i =0; i < ac; ) {
        const char *a = av[i++];
        const char **var = 0;
        if ( strcmp(a, "-data") == 0  && data_fname == 0 )
            var = &data_fname;
        else if ( strcmp(a, "-fasta") == 0 && fasta == 0 )
            var = &fasta;

        if ( var != 0) {
            if (i == ac) {
                fprintf(stderr, "Missing an argument after %s.", a);
                return 1;
            }
            *var = av[i++];
            continue;
        }
        if ( name == 0 ) {
            name = a;
            continue;
        }
        fprintf(stderr, "Unknown argument : %s", a);
        return 1;
    }

    if ( data_fname == NULL )
        error("-data genepred databases is required.");

    if ( fasta == NULL )
        error("-fasta transcripts fasta file is required.");
    if ( init_hgvs_spec(data_fname, fasta) )
        return 1;
    if ( parse_hgvs_name(name) ) {
        error("Failed to parse name string.");
        return 1;
    }
    
    
    return 0;
}

void convert_hgvs()
{
    fill_hgvs_name();
    print_hgvs_core();
}
void release_memory()
{
}
int main(int argc, char **argv)
{
    if ( parse_args(--argc, ++argv) )
        return 1;
    convert_hgvs();
    release_memory();
    return 0;
}
#endif
