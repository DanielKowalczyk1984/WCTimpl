#ifndef _HEAP_H
#define _HEAP_H

#ifdef __cplusplus
extern "C" {
#endif


#include <assert.h>
#define CCutil_MAXDBL (1e30)
#define CCutil_MAXINT (2147483647)
typedef struct heapelm {
    int key;
    void *obj;
} heapelm;

typedef struct Heap_t {
    int     end;
    int     size;
    int    *perm;
    int    *iperm;

    heapelm *elms;
} pmcheap;


int pmcheap_init(pmcheap **heap, int size),
    pmcheap_free(pmcheap  *heap),
    pmcheap_free_all(pmcheap *heap),
    pmcheap_insert(pmcheap *heap, int key, void *obj),
    pmcheap_remove(pmcheap *heap, int href),
    pmcheap_get_key(const pmcheap *heap, int href),
    pmcheap_size(const pmcheap *heap),
    pmcheap_decrease_key(pmcheap *heap, int href, int new_key),
    pmcheap_relabel(pmcheap *heap, int href, int new_key);
void *pmcheap_get_obj(const pmcheap *heap, int href);

void *pmcheap_min(pmcheap *heap),
     pmcheap_reset(pmcheap *heap),
     pmcheap_reset_free(pmcheap *heap);
#endif

#ifdef __cplusplus
}
#endif