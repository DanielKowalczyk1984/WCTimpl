#include <stdio.h>
#include "datastructsol.h"



double partition_order(const void *a, const void *b, void *data){
    (void) data;

    const double *v1 = (double *)&(((const Job*)a)->weight);
    const double *w1 = (double *)&(((const Job*)b)->weight);

    const double *v2 = (double *)&(((const Job*)a)->processingime);
    const double *w2 = (double *)&(((const Job*)b)->processingime);

    double v,w;
    v = *v1/(*v2);
    w = *w1/(*w2); 

    return -(v - w);
}


void partlist_free(partlist *part){
    if(part){
        if(part->list != (GQueue *) NULL) g_queue_free(part->list);
        part->list = (GQueue*) NULL;
    }
}

void partlist_init(partlist *part){
    if(part){
        part->list = g_queue_new();
    }
}

void joblist_init(joblist *jlist){
    if(jlist){
        jlist->part = (partlist *) NULL;
    }
}

void partition_init(partlist *part, joblist *jlist, int nbpart, int jcount){
    int i;

    for(i = 0; i < nbpart; i++){
        partlist_init(&part[i]);
        part[i].key = i;
    }

    for(i = 0; i < jcount;i++){
        joblist_init(&jlist[i]);
    }
}


int partlist_insert_order(partlist *part,joblist *jlist,Job *job){
    int val = 0;

    if(jlist[job->job].part != NULL){
        fprintf(stderr, "Error: double insertion !!!!\n");
        val = 1;
        goto CLEAN;
    }

    jlist[job->job].part = part;

    g_queue_insert_sorted(part->list, job, (GCompareDataFunc)partition_order,NULL);

    CLEAN:
    return val;
}

int partlist_insert(partlist *part, joblist *jlist,Job *job){
    int val = 0;

    if(jlist[job->job].part != NULL){
        fprintf(stderr, "Error: double insertion\n");
        val = 1;
        goto CLEAN;
    }

    jlist[job->job].part = part;

    g_queue_push_head(part->list, job);

    CLEAN:
    return val;
}

int partlist_delete_custom(joblist *jlist, Job *job){
    int val = 0;
    partlist *p = (partlist *)NULL;

    if (jlist[job->job].part == (partlist*) NULL){
        fprintf(stderr, "Error deleting a job that is not assigned\n");
        val = 1;
        goto CLEAN;
    }

    p = (jlist + job->job)->part;
    if(g_queue_remove(p->list, job)){
        jlist[job->job].part = (partlist*)NULL;
    } else {
        printf("We didn't find the job\n");
    }

    CLEAN:
    return val;
}

int partlist_delete(joblist *jlist,Job *job){
    int val = 0;

    partlist *p = NULL;

    if(jlist[job->job].part == NULL){
        fprintf(stderr, "Error deleting a job that is not assigned\n");
        val = 1;
        goto CLEAN;
    }

    p = jlist[job->job].part;
    g_queue_remove(p->list, job);
    jlist[job->job].part = NULL;

    CLEAN:
    return val;
}

void partlist_move(partlist *part, joblist *jlist,Job *job){
    if(jlist[job->job].part !=NULL){
        partlist_delete(jlist,job);
        partlist_insert(part, jlist,job);
    } else {
        partlist_insert(part,jlist,job);
    }
}

void partlist_move_order(partlist *part, joblist *jlist,Job *job){
    if(jlist[job->job].part != NULL){
        partlist_delete(jlist,job);
        partlist_insert_order(part, jlist,job);
    } else {
        partlist_insert_order(part,jlist,job);
    }
}