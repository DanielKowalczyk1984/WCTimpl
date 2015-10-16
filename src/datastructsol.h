#include <glib.h>

typedef struct _Job{
    GList *before;
    GList *after;
    int job;
    int weight;
    int processingime;
    int releasetime;
    int duetime;
} Job;


///////////////////////////
//Definition of partlist //
///////////////////////////

typedef struct _partlist{
    GQueue *list;
    int completiontime;
    int key;
} partlist;

typedef struct _joblist{
    partlist *part;
} joblist;


/////////////////////////////
// Definition of solution  //
/////////////////////////////

typedef struct _solution{
    partlist *part;
    joblist *vlist;
    int *perm;
    int **sumtimes;
    int **sumweights;
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
void vertexlist_init(joblist *vlist);
void partition_init(partlist *part, joblist*vlist,int nbpart, int jcount);
int partlist_insert_order(partlist *part,joblist *vlist,Job *job);
int partlist_insert(partlist *part, joblist *vlist,Job *job);
int partlist_delete(joblist *vlist, Job *job);
int partlist_delete_custom(joblist *vlist,Job *vertex);
void partlist_move(partlist *part, joblist *vlist,Job *vertex);
void partlist_move_order(partlist *part, joblist *vlist,Job *vertex);
double partition_order(const void *a,const void *b,void *data);
int find_vertex(const void *a,const void *b);


////////////////////////////////////////
// Definition of functions of solution//
////////////////////////////////////////

void solution_init(solution *sol);
void solution_free(solution *sol);
int solution_alloc(solution *sol,int nbpart,int jcount);

void solution_unique(solution *sol);
void solution_max(solution *sol);
void solution_print(solution *sol);
void test_SOLUTION(solution *sol);
int solution_copy(solution *dest,solution src);
int solution_update(solution *dest,solution src);
int solution_check(partlist *part,int jcount);