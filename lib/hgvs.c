#include "utils.h"
#include "genepred.h"

struct hgvs_spec {
    struct genepred_spec *data;
    struct hgvs_des des;
} spec;

void hgvs_core_clear(struct hgvs_core *core)
{
    if ( core->name.name1 )
        free(core->name.name1);
    if ( core->name.name2 )
        free(core->name.name2);
    memset(core, 0, sizeof(struct hgvs_core));
}
void hgvs_des_clear(struct hgvs_des *des)
{
    int i;
    for ( i = 0; i < des->i; ++i) {
        hgvs_core_clear(&des->a[i]);
    }
    free(des->a);
    if ( des->chrom )
        free(des->chrom);
    if ( des->ref_length )
        free(des->ref);
    if ( des->alt_length )
        free(des->alt);
    memset(des, 0, sizeof(struct hgvs_des));
}

int init_hgvs_spec(const char *fname, const char *fasta)
{
    memset(&spec, 0, sizeof(struct hgvs_spec));
    spec.data = genepred_spec_init();
    if ( genepred_load_data(spec.data, fname) == NULL )
        return 1;
    if ( genepred_load_fasta(spec.data, fasta) == NULL )
        return 1;
    
    return 0;
}
// Standard HGVS name should be NM_0001.2:c.123A>G; tolerant format could be NM_0001:c.123A>G (no version number);
// 
int check_hgvs_name(char *name)
{
}
int parse_hgvs_name(char *name)
{
    if ( check_hgvs_name(name) )
        return 1;

    return 0;
}
int find_the_block(struct genepred_line *line, int *blk_start, int *blk_end, int pos)
{
}
// return 0 on success, 1 on out of range.
int generate_hgvs_core(struct genepred_line *line, struct hgvs_core *name, int start, int end)
{
    int i;
    int blk_start = 0, blk_end = 0;
    // Locate start.
    for (i = 0; i < line->exoncount; ++i) {
        if ( find_the_block(line, &blk_start, &blk_end, start ) )
            return 1;
        
    }
    // Locate end. For most cases variants are snps, start == end.
    if ( end != start ) {
        for ( i = 0; i < line->exoncount; ++i ) {
         if ( find_the_block(line, &blk_start, &blk_end, end ) )
             return 1;

         
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

    struct genepred_line *line = genepred_retrieve_region(spec->data, des->name, des->start, des->end);

    for ( ;; ) {
        if ( line == NULL )
            break;
        
        if ( des->l == des->m ) {
            des->m += 2;
            des->a = (struct hgvs_core*)realloc(des->a, des->m*sizeof(struct hgvs_core));
        }

        hgvs_core_clear(&des->a[des->l]);
        if ( generate_hgvs_core(line, &des->a[des->l]) == 0 )
            des->l ++;
        struct genepred_line * temp = line;
        line = line->next;
        genepred_line_destroy(temp);
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
}

void convert_hgvs()
{
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
