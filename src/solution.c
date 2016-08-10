#include "datastructsol.h"
#include "util.h"
#include <string.h>

gint comparefunc(const void *a, const void *b, void *data);
gint compare_func(const void *a, const void *b);
gint order_weight(gconstpointer a, gconstpointer b, void *data);


void solution_init(solution *sol) {
    if (sol) {
        sol->part     = (partlist *)NULL;
        sol->vlist    = (joblist *)NULL;
        sol->perm     = (int *)NULL;
        sol->nmachines   = 0;
        sol->njobs   = 0;
        sol->totalweightcomptime = 0;
        sol->dist     = 0;
        sol->iter     = 0;
        sol->fitness  = .0;
    }
}

void solution_free(solution *sol) {
    int i;

    if (sol) {
        for (i = 0; i < sol->nmachines; ++i) {
            partlist_free(sol->part + i);
        }

        CC_IFFREE(sol->part, partlist);
        CC_IFFREE(sol->vlist, joblist);
        CC_IFFREE(sol->perm, int);
        sol->nmachines   = 0;
        sol->totalweightcomptime = 0;
        sol->njobs   = 0;
        sol->dist     = 0;
        sol->iter     = 0;
        sol->fitness  = .0;
    }
}

int solution_alloc(solution *sol, int nmachines, int njobs) {
    int val = 0;
    sol->nmachines  = nmachines;
    sol->njobs = njobs;
    int i;

    if (sol->vlist != NULL || sol->part != NULL) {
        fprintf(stderr, "Error we already allocated memory to part or vlist\n");
        val = 1;
        goto CLEAN;
    }

    sol->part = CC_SAFE_MALLOC(nmachines, partlist);
    CCcheck_NULL_2(sol->part, "Failed to allocate memory to part");

    for (i = 0; i < nmachines; ++i) {
        partlist_init(sol->part + i);
        sol->part[i].sumtimes = CC_SAFE_MALLOC(njobs, int);
        CCcheck_NULL_2(sol->part[i].sumtimes, "Failed to allocate memory");
        fill_int(sol->part[i].sumtimes, njobs, 0);
        sol->part[i].sumweights = CC_SAFE_MALLOC(njobs, int);
        CCcheck_NULL_2(sol->part[i].sumweights, "Failed to allocate memory");
        fill_int(sol->part[i].sumweights, njobs, 0);
        (sol->part + i)->key = i;
    }

    sol->vlist = CC_SAFE_MALLOC(njobs, joblist);
    CCcheck_NULL_2(sol->vlist, "Failed to allocate memory to vlist");
    sol->perm = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(sol->perm, "Failed to allocate memory to perm");

    for (i = 0; i < njobs; ++i) {
        joblist_init(sol->vlist + i);
        sol->perm[i] = -1;
    }

CLEAN:

    if (val) {
        for (i = 0; i < nmachines; i++) {
            partlist_free(sol->part + i);
        }

        CC_IFFREE(sol->part, partlist);
        CC_IFFREE(sol->vlist, joblist);
        CC_IFFREE(sol->perm, int);
    }

    return val;
}

gint comparefunc(const void *a, const void *b, void *data) {
    (void) data;
    const int *v = &(((const Job *)a)->job);
    const int *w = &(((const Job *)b)->job);
    return *v - *w;
}

gint compare_func(const void *a, const void *b) {
    const int *v = &(((const partlist *)a)->completiontime);
    const int *w = &(((const partlist *)b)->completiontime);

    if (*v != *w) {
        return *v - *w;
    } else {
        if (*v == 0 || *w == 0) {
            return *v - *w;
        }

        const int *vv = &(((Job *)((const partlist *)a)->list->head->data)->job);
        const int *ww = &(((Job *)((const partlist *)b)->list->head->data)->job);
        return *vv - *ww;
    }
}

gint order_weight(gconstpointer a, gconstpointer b, void *data) {
    (void) data;
    const int *v = &(((const Job *)a)->weight);
    const int *w = &(((const Job *)b)->weight);
    return -(*v - *w);
}

void solution_max(solution *sol) {
    int max = 0;

    for (int i = 0; i < sol->nmachines; ++i) {
        if (sol->part[i].completiontime > max) {
            max = sol->part[i].completiontime;
        }
    }

    sol->totalweightcomptime = max;
}

void solution_unique(solution *sol) {
    int i;
    int nmachines = sol->nmachines;
    GList *it = (GList *) NULL;
    partlist *temp_partlist = (partlist *) NULL;
    int counter = 0;
    qsort(sol->part, sol->nmachines, sizeof(partlist), compare_func);

    /** Compute permutation */
    if (sol->perm == NULL) {
        sol->perm = CC_SAFE_MALLOC(sol->njobs, int);
    }

    for (i = 0; i < nmachines; i++) {
        temp_partlist = sol->part + i;
        temp_partlist->key = i;

        for (it = temp_partlist->list->head; it; it = it->next) {
            sol->perm[counter] = ((Job *)it->data)->job;
            sol->vlist[((Job *)it->data)->job].part = sol->part + i;
            counter++;
        }
    }
}

void solution_print(solution *sol) {
    for (int i = 0; i < sol->nmachines; ++i) {
        printf("Machine %d: ", sol->part[i].key);

        for (GList *it = sol->part[i].list->head; it; it = it->next) {
            printf("%d ", ((Job *)it->data)->job);
        }

        printf("with C =  %d, wC = %d and %d jobs\n", sol->part[i].completiontime,
               sol->part[i].totcompweight
               , g_queue_get_length(sol->part[i].list));
    }

    printf("with total weighted completion time %d\n", sol->totalweightcomptime);
}

void test_SOLUTION(solution *sol) {
    int i;

    for (i = 0; i <  sol->njobs; i++) {
        printf("%d ", i);

        for (GList *it = sol->vlist[i].part->list->head; it; it = it->next) {
            printf("%d ", ((Job *)it->data)->job);
        }

        printf("\n");
    }
}

int solution_copy(solution *dest, solution src) {
    int val = 0;
    int counter = 0;
    solution_init(dest);
    dest->dist     = src.dist;
    dest->fitness  = src.fitness;
    dest->iter     = src.iter;
    dest->totalweightcomptime = src.totalweightcomptime;
    dest->nmachines   = src.nmachines;
    dest->njobs   = src.njobs;
    val = solution_alloc(dest, dest->nmachines, dest->njobs);
    CCcheck_val_2(val, "Failed in  solution_alloc");

    for (int i = 0; i < dest->nmachines; i++) {
        dest->part[i].key = src.part[i].key;
        dest->part[i].totcompweight = src.part[i].totcompweight;
        g_queue_free(dest->part[i].list);
        dest->part[i].list = (GQueue *) NULL;
        dest->part[i].list = g_queue_copy(src.part[i].list);
        memcpy(dest->part[i].sumtimes, src.part[i].sumtimes, sizeof(int)*dest->njobs);
        memcpy(dest->part[i].sumweights, src.part[i].sumweights,
               sizeof(int)*dest->njobs);
        dest->part[i].completiontime = src.part[i].completiontime;

        for (GList *it = dest->part[i].list->head; it; it = it->next) {
            dest->perm[counter] = ((Job *)it->data)->job;
            dest->vlist[((Job *)it->data)->job].part = dest->part + i;
            counter++;
        }
    }

CLEAN:

    if (val) {
        solution_free(dest);
        CC_IFFREE(dest, solution);
    }

    return val;
}

int solution_update(solution *dest, solution src) {
    int val = 0;
    int counter = 0;
    dest->dist     = src.dist;
    dest->fitness  = src.fitness;
    dest->iter     = src.iter;
    dest->totalweightcomptime = src.totalweightcomptime;
    dest->nmachines   = src.nmachines;
    dest->njobs   = src.njobs;

    for (int i = 0; i < dest->nmachines; i++) {
        g_queue_free(dest->part[i].list);
        dest->part[i].totcompweight = src.part[i].totcompweight;
        dest->part[i].list = g_queue_copy(src.part[i].list);
        dest->part[i].completiontime = src.part[i].completiontime;
        memcpy(dest->part[i].sumtimes, src.part[i].sumtimes, sizeof(int)*dest->njobs);
        memcpy(dest->part[i].sumweights, src.part[i].sumweights,
               sizeof(int)*dest->njobs);

        for (GList *it = dest->part[i].list->head; it; it = it->next) {
            dest->perm[counter] = ((Job *)it->data)->job;
            dest->vlist[((Job *)it->data)->job].part = dest->part + i;
            counter++;
        }
    }

    return val;
}

void solution_calc(solution *sol, Job *jobarray) {
    int i, j;
    int njobs = sol->njobs;
    int nmachines = sol->nmachines;
    GList *it = (GList *) NULL;
    Job *temp_job = (Job *) NULL;
    partlist *temp_partlist = (partlist *) NULL;
    sol->totalweightcomptime = 0;

    /** Order in WSPT order and compute objective value of this solution */
    for (i = 0; i < nmachines; ++i) {
        temp_partlist = sol->part + i;
        temp_partlist->completiontime = 0;
        g_queue_sort(temp_partlist->list, (GCompareDataFunc)comparefunc, NULL);

        for (it = temp_partlist->list->head; it; it = g_list_next(it)) {
            temp_job = ((Job *)it->data);
            temp_partlist->completiontime += temp_job->processingime ;
            temp_partlist->totcompweight += temp_job->weight *
                                            temp_partlist->completiontime;
            sol->totalweightcomptime += temp_partlist->completiontime * temp_job->weight;
        }
    }

    /** Initialize sum weights and sum times */
    for (i = 0; i < nmachines; ++i) {
        temp_partlist = sol->part + i;

        if (sol->vlist[0].part == temp_partlist) {
            temp_partlist->sumtimes[0] = jobarray[0].processingime;
        } else  {
            temp_partlist->sumtimes[0] = 0;
        }

        temp_partlist->sumweights[sol->njobs - 1] = 0;
    }

    /** Compute sum weights and sum times */
    for (j = 1; j < njobs; j++) {
        for (i = 0; i < nmachines; i++) {
            temp_partlist = sol->part + i;

            if (sol->vlist[j].part == temp_partlist) {
                temp_partlist->sumtimes[j] = temp_partlist->sumtimes[j - 1] +
                                             jobarray[j].processingime;
            } else {
                temp_partlist->sumtimes[j] = temp_partlist->sumtimes[j - 1];
            }
        }
    }

    for (j = njobs - 1; j > 0; --j) {
        for (i = 0; i < nmachines; ++i) {
            temp_partlist = sol->part + i;

            if (sol->vlist[j].part == temp_partlist) {
                temp_partlist->sumweights[j - 1] = temp_partlist->sumweights[j] +
                                                   jobarray[j].weight;
            } else {
                temp_partlist->sumweights[j - 1] = temp_partlist->sumweights[j];
            }
        }
    }
}

void solution_wct(solution *sol) {
    int i;
    int nmachines = sol->nmachines;
    GList *it = (GList *) NULL;
    Job *temp_job = (Job *) NULL;
    partlist *temp_partlist = (partlist *) NULL;
    sol->totalweightcomptime = 0;

    /** Order in WSPT order and compute objective value of this solution */
    for (i = 0; i < nmachines; ++i) {
        temp_partlist = sol->part + i;
        temp_partlist->completiontime = 0;
        g_queue_sort(temp_partlist->list, (GCompareDataFunc)comparefunc, NULL);

        for (it = temp_partlist->list->head; it; it = g_list_next(it)) {
            temp_job = ((Job *)it->data);
            temp_partlist->completiontime += temp_job->processingime;
            sol->totalweightcomptime += temp_partlist->completiontime * temp_job->weight;
        }
    }
}

void partlist_permquicksort(int *perm, partlist *part, int nbpart,
                            int (*functionPtr)(partlist *, partlist *)) {
    int i, j, temp;
    partlist t;

    if (nbpart <= 1) {
        return;
    }

    CC_SWAP(perm[0], perm[(nbpart - 1) / 2], temp);
    i = 0;
    j = nbpart;
    memcpy(&t, &(part[perm[0]]), sizeof(partlist));

    while (1) {
        do {
            i++;
        } while (i < nbpart && (*functionPtr)(&(part[perm[i]]), &t));

        do {
            j--;
        } while ((*functionPtr)(&t, &(part[perm[j]])));

        if (j < i) {
            break;
        }

        CC_SWAP(perm[i], perm[j], temp);
    }

    CC_SWAP(perm[0], perm[j], temp);
    partlist_permquicksort(perm, part, j, (*functionPtr));
    partlist_permquicksort(perm + i, part, nbpart - i, (*functionPtr));
}



