#ifndef SORT_LIST_HEADER
#define SORT_LIST_HEADER

struct list_lite {
    struct list_lite *next;
};

extern int count_list(const void *list);

typedef int comp_func(const void *elem1, const void *elem2);

extern int sort_list(void *list, comp_func *func);

#endif


