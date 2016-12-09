#include "val_list.h"

struct vlist_spec *vlist_init()
{
    struct vlist_spec *list = malloc(sizeof(struct vlist_spec));
    memset(list, 0, sizeof (struct vlist_spec) );
    return list;
}
int vlist_add(struct vlist_spec *spec, union value val)
{
    int i;
    for ( i = 0; i < spec->n_nodes; ++i ) {
        struct 
        struct val_node *node = spec
        if (
    }
    
}
