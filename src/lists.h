
#include <string.h>

#define LIST_INIT_DEFAULT_NMAX 10

typedef struct list {
	void  *items;        /* pointer to raw bytes */
	size_t elem_size;    /* sizeof(element type) */
	unsigned int N;
	unsigned int Nmax;
} s_list;


s_list list_initialize(size_t elem_size, size_t Nmax);
int list_ensure_capacity(s_list *list, size_t need);
int list_push(s_list *list, const void *elem);
void *list_get_ptr(s_list *list, size_t idx);
int list_pop(s_list *list, void *out_elem);
void free_list(s_list *list);


