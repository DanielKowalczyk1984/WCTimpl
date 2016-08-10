#include <stdio.h>
#include "datastructsol.h"
#include "util.h"

int partition_order(const void *a, const void *b, void *data) {
    (void) data;
    const int *v1 = (const int *) & (((const Job *)a)->job);
    const int *w1 = (const int *) & (((const Job *)b)->job);
    return (*v1 - *w1);
}

void partlist_free(partlist *part) {
    if (part) {
        if (part->list != (GQueue *) NULL) {
            g_queue_free(part->list);
        }

        part->list = (GQueue *) NULL;
        CC_IFFREE(part->sumweights, int);
        CC_IFFREE(part->sumtimes, int);
    }
}

void partlist_init(partlist *part) {
    if (part) {
        part->list = g_queue_new();
        part->sumtimes = (int *) NULL;
        part->sumweights = (int *) NULL;
        part->completiontime = 0;
        part->totcompweight = 0;
    }
}

void joblist_init(joblist *jlist) {
    if (jlist) {
        jlist->part = (partlist *) NULL;
    }
}

void partition_init(partlist *part, joblist *jlist, int nbpart, int njobs) {
    int i;

    for (i = 0; i < nbpart; i++) {
        partlist_init(&part[i]);
        part[i].sumtimes = CC_SAFE_MALLOC(njobs, int);
        part[i].sumweights = CC_SAFE_MALLOC(njobs, int);
        part[i].key = i;
    }

    for (i = 0; i < njobs; i++) {
        joblist_init(&jlist[i]);
    }
}

int partlist_insert_order(partlist *part, joblist *jlist, Job *job, int njobs) {
    int val = 0;
    int i;

    if (jlist[job->job].part != NULL) {
        fprintf(stderr, "Error: double insertion !!!!\n");
        val = 1;
        goto CLEAN;
    }

    jlist[job->job].part = part;
    g_queue_insert_sorted(part->list, job, (GCompareDataFunc)partition_order, NULL);

    for (i = job->job; i < njobs; ++i) {
        part->sumtimes[i] += job->processingime;
    }

    for (i = 0; i < job->job; ++i) {
        part->sumweights[i] += job->weight;
    }

    part->totcompweight += job->weight * part->sumtimes[job->job] +
                           job->processingime * part->sumweights[job->job];
    part->completiontime += job->processingime;
CLEAN:
    return val;
}

int partlist_insert(partlist *part, joblist *jlist, Job *job) {
    int val = 0;

    if (jlist[job->job].part != NULL) {
        fprintf(stderr, "Error: double insertion\n");
        val = 1;
        goto CLEAN;
    }

    jlist[job->job].part = part;
    g_queue_push_tail(part->list, job);
    part->completiontime += job->processingime;
    part->totcompweight += part->completiontime * job->weight;
CLEAN:
    return val;
}

int partlist_delete_custom(joblist *jlist, Job *job, int njobs) {
    int i, val = 0;
    partlist *p = (partlist *)NULL;

    if (jlist[job->job].part == (partlist *) NULL) {
        fprintf(stderr, "Error deleting a job that is not assigned\n");
        val = 1;
        goto CLEAN;
    }

    p = (jlist + job->job)->part;

    if (g_queue_remove(p->list, job)) {
        jlist[job->job].part = (partlist *)NULL;
        p->totcompweight -= job->weight * p->sumtimes[job->job] + job->processingime *
                            p->sumweights[job->job];
        p->completiontime -= job->processingime;

        for (i = job->job; i < njobs; ++i) {
            p->sumtimes[i] -= job->processingime;
        }

        for (i = 0; i < job->job; ++i) {
            p->sumweights[i] -= job->weight;
        }
    } else {
        printf("We didn't find the job\n");
    }

CLEAN:
    return val;
}

int partlist_delete(joblist *jlist, Job *job) {
    int val = 0;
    partlist *p = NULL;

    if (jlist[job->job].part == NULL) {
        fprintf(stderr, "Error deleting a job that is not assigned\n");
        val = 1;
        goto CLEAN;
    }

    p = jlist[job->job].part;
    g_queue_remove(p->list, job);
    p->completiontime -= job->processingime;
    jlist[job->job].part = NULL;
CLEAN:
    return val;
}

void partlist_move(partlist *part, joblist *jlist, Job *job) {
    if (jlist[job->job].part != NULL) {
        partlist_delete(jlist, job);
        partlist_insert(part, jlist, job);
    } else {
        partlist_insert(part, jlist, job);
    }
}

void partlist_move_order(partlist *part, joblist *jlist, Job *job, int njobs) {
    if (jlist[job->job].part != NULL) {
        partlist_delete_custom(jlist, job, njobs);
        partlist_insert_order(part, jlist, job, njobs);
    } else {
        partlist_insert_order(part, jlist, job, njobs);
    }
}
