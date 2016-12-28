#ifndef SORT_LIST_HEADER
#define SORT_LIST_HEADER

struct list_lite {
    struct list_lite *next;
};

extern int count_list(const void *list);

typedef int comp_func(const void *elem1, const void *elem2);

typedef void delete_func(void *list);

extern int sort_list(void *list, comp_func *func);
extern void list_lite_delete(void *list, delete_func delete);
#endif


