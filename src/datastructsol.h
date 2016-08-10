#ifndef DATASTRUCTSOL_H
#define DATASTRUCTSOL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <glib.h>

typedef struct _Job {
    int job;
    int weight;
    int processingime;
    int releasetime;
    int duetime;
    GList *ready_list;
    GList *due_list;
} Job;


////////////////////////////////////////
// Definition of adjNode and adjGraph //
////////////////////////////////////////


typedef struct adjNode {
    int *adj;
    int node;
    int weight;
    int residual_weight;
    int active;
    int degree;
    int color;
    int profit;
} adjNode;

typedef struct adjGraph {
    adjNode *nodelist;
    int *adjspace;
    int **adjMatrix;
    int vcount;
    int ecount;

    int *perm;
    int *weightorder;
    int *makespan;
    int makespancolor;
    int ncolor;
    int totweight;
    int flag;

} adjGraph;

///////////////////////////
//Definition of partlist //
///////////////////////////

typedef struct _partlist {
    GQueue *list;
    int *sumtimes;
    int *sumweights;
    int completiontime;
    int totcompweight;
    int key;
} partlist;

typedef struct _joblist {
    partlist *part;
} joblist;


/////////////////////////////
// Definition of solution  //
/////////////////////////////

typedef struct _solution {
    partlist *part;
    joblist *vlist;
    int *perm;
    int totalweightcomptime;
    int njobs;
    int nmachines;
    int dist;
    int iter;
    double fitness;
} solution;


//////////////////////////////////////////
// Definition of functions of partlist  //
//////////////////////////////////////////

void partlist_free(partlist *part);
void partlist_init(partlist *part);
void joblist_init(joblist *vlist);
void partition_init(partlist *part, joblist *vlist, int nbpart, int jcount);
int partlist_insert_order(partlist *part, joblist *vlist, struct _Job *job,
                          int njobs);
int partlist_insert(partlist *part, joblist *vlist, struct _Job *job);
int partlist_delete(joblist *vlist, struct _Job *job);
int partlist_delete_custom(joblist *vlist, struct _Job *job, int njobs);
void partlist_move(partlist *part, joblist *vlist, struct _Job *job);
void partlist_move_order(partlist *part, joblist *vlist, struct _Job *job,
                         int jobs);
int partition_order(const void *a, const void *b, void *data);
int find_vertex(const void *a, const void *b);

void partlist_permquicksort(int *perm, partlist *part, int nbpart,
                            int (*functionPtr)(partlist *, partlist *));
int partlist_more_totweight(partlist *c1, partlist *c2);


////////////////////////////////////////
// Definition of functions of solution//
////////////////////////////////////////

void solution_init(solution *sol);
void solution_free(solution *sol);
int solution_alloc(solution *sol, int nbpart, int jcount);

void solution_unique(solution *sol);
void solution_calc(solution *sol, struct _Job *joblist);
void solution_max(solution *sol);
void solution_print(solution *sol);
void test_SOLUTION(solution *sol);
void solution_wct(solution *sol);
int solution_copy(solution *dest, solution src);
int solution_update(solution *dest, solution src);
int solution_check(partlist *part, int jcount);

/**
 * scheduleset.c
 */

typedef struct Scheduleset {
    int count;
    int size;
    int age;
    int totweight;
    int totwct;
    GHashTable *table;
    int *members;
    int *C;
    int id;
} Scheduleset;

void Scheduleset_SWAP(Scheduleset *c1, Scheduleset *c2, Scheduleset *t);

/*Initialization and free memory for the colorset*/
void Scheduleset_init(Scheduleset *set);
void Scheduleset_free(Scheduleset *set);
void Schedulesets_free(Scheduleset **set, int *nsets);
int  COLORcopy_sets(Scheduleset **dsts, int *nsets, Scheduleset *src_s,
                    int src_nsets);

/*Check if the coloring is feasible*/
int Scheduleset_check_set(Scheduleset *set, int vcount);
int Scheduleset_check(Scheduleset *set, int ccount, int vcount);

/*Transformation of covers*/
int transform_into_coloring(int vcount, int *ncolors,
                            Scheduleset **colorclasses);

/*Sorting Schedulesets*/
void Scheduleset_quicksort(Scheduleset *cclasses, int ccount,
                           int (*functionPtr)(Scheduleset *, Scheduleset *));
void Scheduleset_permquicksort(int *perm, Scheduleset *cclasses, int ccount,
                               int (*functionPtr)(Scheduleset *, Scheduleset *));
int Scheduleset_less(Scheduleset *c1, Scheduleset *c2);
int Scheduleset_more(Scheduleset *c1, Scheduleset *c2);
int Scheduleset_less_wct(Scheduleset *c1, Scheduleset *c2);
int print_schedule(Scheduleset *cclasses, int ccount);
int Scheduleset_max(Scheduleset *cclasses, int ccount);
int update_Schedulesets(Scheduleset **dst, int *ndst, Scheduleset *src,
                        int nsrc);
int add_Schedulesets(Scheduleset **dst, int *ndst, Scheduleset *src, int nsrc);
int Scheduleset_less_totweight(Scheduleset *c1, Scheduleset *c2);
int Scheduleset_more_totweight(Scheduleset *c1, Scheduleset *c2);
int partlist_to_Scheduleset(partlist *part, int nbpart, int njobs,
                            Scheduleset **classes, int *ccount);


///////////////////////////////////////////////////////
// Definition of functions of adjGraph and adjNode   //
///////////////////////////////////////////////////////


/* Build and Free adjGraph and adjNode */
int adjGraph_build(adjGraph *G, int vcount, int ecount, const int elist[],
                   const int weights[]);
int adjGraph_buildquick(adjGraph *G, int vcount, int ecount, int *elist);
int adjGraph_copy(adjGraph *Gdst, const adjGraph *Gsrc);
void adjGraph_init(adjGraph *G);
void adjGraph_free(adjGraph *G);
void adjGraph_freequick(adjGraph *G);
int adjGraph_simplify(adjGraph *G);
int adjGraph_simplifyquick(adjGraph *G);
int adjGraph_complement(adjGraph *Gc, const adjGraph *G);
int adjGraph_get_elist(int *ecount, int *elist[], const adjGraph *G);
int read_adjlist(char *f, int *pvcount, int *pecount, int **pelist,
                 int **weightlist);

/* Some usefull functions*/
int adjGraph_delete_unweighted(adjGraph *G, int **new_nweights,
                               const int nweights[]);
int adjGraph_connectedness(const adjGraph *G);
void graph_print(int ecount, const int elist[]);
int adjGraph_reset_schedule(adjGraph *G);
int adjGraph_adjust_schedule(adjGraph *Gdst, const adjGraph *Gsrc);
void reset_color(adjGraph *G);

/* Some Sorting functions on adjGraph*/
void adjGraph_quicksort_perm(adjNode *nodelist, int *perm, int vcount,
                             int (*compareFunction)(adjNode *, adjNode *));
void adjGraph_quicksort(adjNode *nodelist, int vcount,
                        int (*compareFunction)(adjNode *, adjNode *));
void adjGraph_sort_adjlists_by_id(adjGraph *G);
int adjGraph_degree(adjNode *n1, adjNode *n2);
int adjGraph_weight(adjNode *n1, adjNode *n2);
int adjGraph_invdegree(adjNode *n1, adjNode *n2);
int  adjGraph_print(int ecount, const int elist[]);
#ifdef __cplusplus
}
#endif

#endif // DATASTRUCTSOL_H
