////////////////////////////////////////////////////////////////
//                                                            //
//  heap.c                                                    //
//  PMC                                                       //
//                                                            //
//  Created by Daniel on 21/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#include <assert.h>
#include "defs.h"
#include "util.h"

/* #define HEAP_INTEGRITY_CHECKS 1 */


static int pmcheap_empty(pmcheap *heap) {
    assert(heap->end >= 0);
    return !(heap->end);
}

MAYBE_UNUSED
static int pmcheap_integrity(pmcheap *heap) {
    int val = 0;
    int i;
    int *perm = heap->perm;
    int *iperm = heap->iperm;

    for (i = 0 ; i < heap->end; ++i) {
        val = !(iperm[perm[i]] == i);
        CCcheck_val(val, "Failed: iperm[perm[i]] == i");
        val = !(perm[iperm[i]] == i);
        CCcheck_val(val, "Failed: perm[iperm[i]] == i");

        if (i > 0) {
            int parent = i >> 1;
            val = (heap->elms[perm[parent]].key < heap->elms[perm[i]].key);
            CCcheck_val(val, "Failed: heap order");
        }
    }

    return val;
}

#ifdef HEAP_INTEGRITY_CHECKS
#define HEAP_INTEGRITY(rval,heap,msg) {              \
        rval = pmcheap_integrity(heap);           \
        COLORcheck_rval(rval,msg);                     \
    }
#else
#define HEAP_INTEGRITY(rval,heap,msg)
#endif


int pmcheap_init(pmcheap **heap, int size) {
    int val = 0;

    if (size == 0) {
        val = 1;
        goto CLEAN;
    }

    *heap = (pmcheap *) CC_SAFE_MALLOC(1, pmcheap);
    CCcheck_NULL_2(*heap, "Failed to allocate heap");
    (*heap)->perm  = (int *) NULL;
    (*heap)->iperm = (int *) NULL;
    (*heap)->elms  = (heapelm *) NULL;
    (*heap)->end = 0;
    size += 2;
    (*heap)->size = size;
    (*heap)->perm = (int *) CC_SAFE_MALLOC(size, int);
    CCcheck_NULL_2((*heap)->perm, "Failed to allocate (*heap)->perm");
    (*heap)->iperm = (int *) CC_SAFE_MALLOC(size, int);
    CCcheck_NULL_2((*heap)->iperm, "Failed to allocate (*heap)->iperm");
    (*heap)->elms = (heapelm *) CC_SAFE_MALLOC(size, heapelm);
    CCcheck_NULL_2((*heap)->elms, "Failed to allocate (*heap)->elms");
    /* Use sentenials at beginning and end.*/
    (*heap)->elms[0].key       =   -CCutil_MAXINT;
    (*heap)->perm[0] = (*heap)->iperm[0] = 0;
    (*heap)->elms[size - 1].key  =  CCutil_MAXINT;
    (*heap)->perm[size - 1] = (*heap)->iperm[size - 1] = size - 1;
    pmcheap_reset(*heap);
    HEAP_INTEGRITY(val, *heap,
                   "pmcheap_integrity failed in pmcheap_relabel.");
CLEAN:

    if (val) {
        pmcheap_free(*heap);
        *heap = (pmcheap *) NULL;
    }

    return val;
}

int pmcheap_free(pmcheap  *heap) {
    if (heap) {
        if (heap->perm) {
            free(heap->perm);
        }

        if (heap->iperm) {
            free(heap->iperm);
        }

        if (heap->elms) {
            free(heap->elms);
        }

        free(heap);
    }

    return 0;
}

int pmcheap_free_all(pmcheap *heap) {
    if (heap) {
        for (int i = 1; i <= heap->end; ++i) {
            if (heap->elms[heap->perm[i]].obj) {
                free(heap->elms[heap->perm[i]].obj);
            }
        }

        CC_IFFREE(heap->elms, heapelm)
        CC_IFFREE(heap->perm, int)
        CC_IFFREE(heap->iperm, int)
        free(heap);
    }

    return 0;
}

void pmcheap_reset_free(pmcheap *heap) {
    int i;

    for (i = 1; i  <= heap->size; ++i) {
        if (heap->elms[heap->perm[i]].obj) {
            free(heap->elms[heap->perm[i]].obj);
        }
    }

    for (i = 1; i <= heap->size; i++) {
        heap->elms[i].obj = NULL;
        heap->elms[i].key = CCutil_MAXINT;
        heap->perm[i] = heap->iperm[i] = i;
    }

    heap->end = 0;
}

void pmcheap_reset(pmcheap *heap) {
    int i;
    heap->end = 0;

    for (i = 1; i + 1 < heap->size; ++i) {
        heap->elms[i].obj = NULL;
        heap->elms[i].key = CCutil_MAXINT;
        heap->perm[i] = heap->iperm[i] = i;
    }

#ifdef HEAP_INTEGRITY_CHECKS
    assert(!pmcheap_integrity(heap));
#endif
}

static int pmcheap_liftup(pmcheap *heap,
                          int           pos) {
    int swaps   = 0;
    int href     = heap->perm[pos];
    int *perm   = heap->perm;
    int *iperm  = heap->iperm;
    int  parent = pos >> 1;
    int key = heap->elms[href].key;

    /* The sentinel at index 0 will stop the loop.*/
    while (heap->elms[perm[parent]].key > key) {
        /* Move the parent down .*/
        perm[pos] = perm[parent];
        iperm[perm[pos]] = pos;
        pos     = parent;
        parent  >>= 1;
        ++swaps;
    }

    /* If elm at href was lifted up (to pos), update perm arrays.*/
    if (href != heap->perm[pos]) {
        perm[pos]         = href;
        iperm[href]        = pos;
    }

    return swaps;
}

static int pmcheap_siftdown(pmcheap *heap, int pos) {
    int swaps    = 0;
    int end_half = heap->end / 2;
    int *perm  = heap->perm;
    int *iperm = heap->iperm;
    heapelm *elms = heap->elms;
    int  minc, rightc;
    int  ref = perm[pos];
    int  key = heap->elms[ref].key;

    while (pos <= end_half) {
        minc  = pos << 1;  /* j = k*2 */
        rightc = minc + 1;

        /* set minc to minimum of left and right child */
        if (elms[perm[minc]].key  >  elms[perm[rightc]].key) {
            minc = rightc;
        }

        if (key <= elms[perm[minc]].key) {
            break;
        }

        /* move element 'minc' up in the heap */
        perm[pos] = perm[minc];
        iperm[perm[pos]] = pos;
        ++swaps;
        pos = minc;
    }

    perm[pos] = ref;
    iperm[ref] = pos;
    return swaps;
}


int pmcheap_insert(pmcheap *heap, int key, void *obj) {
    int val = 0;
    (heap->end)++;

    if (heap->end  >= heap->size) {
        int i;
        heap->size = heap->size * 2 + 2;
        /* realloc memory */
        heap->perm = (int *) CCutil_reallocrus(heap->perm,
                                               heap->size * sizeof(int));
        CCcheck_NULL_2(heap->perm, "Failed to reallocate heap->perm");
        heap->iperm = (int *) CCutil_reallocrus(heap->iperm,
                                                heap->size * sizeof(int));
        CCcheck_NULL_2(heap->iperm, "Failed to reallocate heap->iperm");
        heap->elms = (heapelm *) CCutil_reallocrus(heap->elms,
                     heap->size * sizeof(heapelm));
        CCcheck_NULL_2(heap->elms, "Failed to reallocate heap->elms");

        for (i = heap->end; i < heap->size; ++i) {
            heap->elms[i].obj = NULL;
            heap->elms[i].key = CCutil_MAXINT;
            heap->perm[i] = heap->iperm[i] = i;
        }

#ifdef HEAP_INTEGRITY_CHECKS
        assert(!pmcheap_integrity(heap));
#endif
    }

    heap->elms[heap->perm[heap->end]].obj  = obj;
    heap->elms[heap->perm[heap->end]].key  = key;
    pmcheap_liftup(heap, heap->end);
    HEAP_INTEGRITY(rval, heap, "pmcheap_integrity failed in pmcheap_insert.");
CLEAN:
    return val;
}

int pmcheap_remove(pmcheap *heap, int href) {
    int rval = 0;
    int heap_pos = heap->iperm[href];

    if (heap_pos > heap->end) {
        fprintf(stderr, "pmcheap_remove error: href does not exist!\n");
        rval = 1;
        goto CLEANUP;
    }

    heap->elms[href].key = CCutil_MAXINT;

    /* A successive uplift of the minimum child might be faster but for
     the time beeing I'm too lazy to implement this.  Instead I lift
     the last element to the hole and sift it down.
     */
    if (heap_pos < heap->end) {
        /* swap last elm to heap_pos  */
        heap->perm[heap_pos] = heap->perm[heap->end];
        heap->perm[heap->end] = href;
        /* swap iperm */
        heap->iperm[heap->perm[heap_pos]] = heap_pos;
        heap->iperm[heap->perm[heap->end]] = heap->end;
        /* Move down elm at index heap_pos.
         It cannot travel to heap->end, as that element has
         key of CCutil_MAXINT now.
         */
        rval = pmcheap_relabel(heap, heap->perm[heap_pos],
                               heap->elms[heap->perm[heap_pos]].key);
    }

    (heap->end)--;
    HEAP_INTEGRITY(rval, heap, "pmcheap_integrity failed in pmcheap_remove.");
CLEANUP:
    return rval;
}

void *pmcheap_min(pmcheap *heap) {
    int href;
    void *obj;
#ifdef HEAP_INTEGRITY_CHECKS
    assert(!pmcheap_integrity(heap));
#endif

    if (pmcheap_empty(heap)) {
        return (void *) NULL;
    }

    href = heap->perm[1];
    obj = heap->elms[href].obj;
    /* swap last elm to index 1 */
    heap->perm[1] = heap->perm[heap->end];
    heap->perm[heap->end] = href;
    /* swap iperm */
    heap->iperm[heap->perm[1]] = 1;
    heap->iperm[heap->perm[heap->end]] = heap->end;
    heap->elms[heap->perm[heap->end]].key = CCutil_MAXINT;
    heap->elms[heap->perm[heap->end]].obj = NULL;
    (heap->end)--;
    /* Move down elm at index 1. */
    pmcheap_siftdown(heap, 1);
#ifdef HEAP_INTEGRITY_CHECKS
    assert(!pmcheap_integrity(heap));
#endif
    return obj;
}

int pmcheap_get_key(const pmcheap *heap, int href) {
    assert(href < heap->size);
    return heap->elms[href].key;
}

int pmcheap_size(const pmcheap *heap) {
    return heap->end;
}

void *pmcheap_get_obj(const pmcheap *heap, int href) {
    assert(href < heap->size);
    return heap->elms[href].obj;
}




int pmcheap_decrease_key(pmcheap *heap, int href, int new_key) {
    int rval = 0;
    int heap_pos = heap->iperm[href];

    if (heap->elms[href].key < new_key) {
        fprintf(stderr,
                "pmcheap_decrease_key error: new_key is greater than old key!\n");
        rval = 1;
        goto CLEANUP;
    }

    heap->elms[href].key = new_key;
    pmcheap_liftup(heap, heap_pos);
    HEAP_INTEGRITY(rval, heap,
                   "pmcheap_integrity failed in pmcheap_decrease_key.\n");
CLEANUP:
    return rval;
}

int pmcheap_relabel(pmcheap *heap, int href, int new_key) {
    int rval = 0;
    int heap_pos = heap->iperm[href];
    heap->elms[href].key = new_key;

    if (!pmcheap_liftup(heap, heap_pos)) {
        pmcheap_siftdown(heap, heap_pos);
    }

    HEAP_INTEGRITY(rval, heap, "pmcheap_integrity failed in pmcheap_relabel.");
#ifdef HEAP_INTEGRITY_CHECKS
CLEANUP:
#endif
    return rval;
}
