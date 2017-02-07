#include "utils.h"
#include "genepred.h"

struct hgvs_spec {
    struct genepred_spec *data;
    struct hgvs_des des;
} spec;

void hgvs_core_clear(struct hgvs_core *core)
{
    if ( core->name.name )
        free(core->name.name);
    if ( core->name.name1 )
        free(core->name.name1);
    memset(core, 0, sizeof(struct hgvs_core));
}
void hgvs_des_clear(struct hgvs_des *des)
{
    int i;
    for ( i = 0; i < des->i; ++i) {
        hgvs_core_clear(&des->a[i]);
    }
    free(des->a);
    if ( des->ref_length )
        free(des->ref);
    if ( des->alt_length )
        free(des->alt);
    memset(des, 0, sizeof(struct hgvs_des));
}

int init_hgvs_spec(const char *fname)
{
    memset(&spec, 0, sizeof(struct hgvs_spec));

    return 0;
}

int parse_hgvs_name(char *name)
{
}

int fill_hgvs_name()
{
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

int main(int argc, char **argv)
{
    if ( parse_args(--argc, ++argv) )
        return 1;
    convert_hgvs();
    release_memory();
    return 0;
}
#endif
