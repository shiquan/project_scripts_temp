// val_list.h is designed for caching uncontinue values, the values will be kept in
// one way list, sorted in default

#ifndef VAL_LIST_HEADER
#define VAL_LIST_HEADER

union value {
    int i;
    float f;
};

struct val_node {
    uint32_t count;
    union value val;
};

struct vlist_node {
    // point to next node
    struct vlist_node *next;
    struct val_node node;
};

enum val_type {
    val_is_integer,
    val_is_float,
};

struct vlist_spec {
    // length of list
    int n_nodes;
    // total counts, sum of count value of each node
    int n_counts;
    // head node
    struct vlist_node *head;
    struct vlist_node *tail;
};

extern struct vlist_spec *vlist_init();

extern int vlist_add(struct vlist_spec *spec, union value val);

extern float vlist_median(struct vlist_spec *spec);

extern float vlist_average(struct vlist_spec *spec);

extern float vlist_sd(struct vlist_spec *spec);

extern void vlist_destroy(struct vlist_spec *spec); 

// convert list structure to array structure
extern struct val_node *convert_vlist_array(struct vlist_spec *spec);

#endif
