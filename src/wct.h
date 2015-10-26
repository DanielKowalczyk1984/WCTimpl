#ifndef _WCT_H
#define _WCT_H


#include "wctprivate.h"

/**
 * greedy.c
 */

typedef struct _pair_job_machine
{
    int job;
    int machine;
} pair_job_machine;


int random_rcl_assignment(Job *jobarray, int njobs, int nmachines, solution* new_sol, GRand *rand_);
int construct_wspt(Job* jobarray, int njobs, int  nmachines, solution*  new_sol);
int random_assignment( Job *jobarray, int njobs, int nmachines, solution *new_sol, GRand *rand_ );

/*
compare_functions
 */
gint compare_func1(const void *a, const void *b, void *user_data);


/**
 * lowerbound.c
 */

typedef struct _MACHINE {
    double totweight;
    double totcompletion;
} MACHINE;


int lowerbound_cw(Job *array, int njobs, int nmachines);
int lowerbound_cp(Job *array, int njobs, int nmachines);
int lowerbound_eei(Job* array, int njobs, int nmachines);

/**
 * compare_functions
 */

int compare_cw(BinomialHeapValue a, BinomialHeapValue b);


/**
 * scatter.c
 */

typedef struct min_max{
    int n;
    int min;
    solution *new_sol;
} min_max;



/**
 * Init and free functions
 */
void REFSET_init(REFSET *rs);
void REFSET_free(REFSET *rs);
void P_init(P *p);
void P_free(P *p);
void SS_init(SS *problem,int Q,int D,double timelimit);
void SS_free(SS *problem);
void free_list2(REFSET *rs);

/**
 * Compare functions
 */

int order_refset(const void *a,const void *b);
int order_makespan(const void *a, const void *b,void *data);
int order_distance(const void *a,const void *b,void *data);
int order_makespan_list(const void *a, const void *b);


/**
 * For each functions
 */

void distance_min_max(void *data,void *user_data);
void max_dist(void *data,void *user_data);
void free_sol(void *data,void *user_data);
void assign_iter(void *data, void *user_data);
void print_sol(void *data,void *user_data);
void print_makespan(void *data, void *user_data);
void refset_dist(void *data,void *user_data);
void for_each_comp_fitness(void *data,void *user_data);

/**
 * Print functions
 */

void print_pool(SS *scatter_search);
void print_refset(SS *scatter_search);
void print_distance(SS *scatter_search);
void print_pool_n(SS *scatter_search,int n);
void print_pool_makespan(SS *scatter_search);
void print_refset_makespan(SS *scatter_search);
void print_list1(SS *scatter_search);

/**
 * Scatter search functions
 */

int SSproblem_definition(
    SS *problem,
    int b1,
    int b2,
    double timelimit,
    int combine_method,
    int njobs,
    int nmachine,
    Job *joblist,
    int lowerbound);

void add_solution_pool(SS *scatter_search, solution *new_sol);
int solution_in_pool(SS *scatter_search,solution *new_sol);
int solution_in_refset(SS *scatter_search, solution *new_sol);

/**
 * Distance functions 
 */
int SSrefset_distance(SS *scatter_search,solution *new_sol);
void *maximum_distance(SS *scatter_search);
int update_distance(SS *scatter_search,solution *sol);


/**
 * Combination functions
 */
int combine_GPX(SS *scatter_search, GQueue *queue, int *subsetsol,int nbelements, solution *new_sol);
int combine_PM(SS *scatter_search,GQueue *queue,int* subsetsol,int nbelements,solution *new_sol);
 

int SSCreate_refset(SS *scatter_search);
int add_solution_refset(SS *scatter_search);

/**
 * Update functions
 */
int static_update(SS *scatter_search);
int dynamic_update(SS *scatter_search, GQueue *list, solution *new_sol);
int diversification_update(SS *scatter_search);

int SSrun_scatter_search(SS *scatter_search);

void print_pool_totalweightcomptime( SS *scatter_search );
void print_refset_totalweightcomptime( SS *scatter_search );

/**
 * Help functions
 */

int compute_fitness(SS *scatter_search);

/**
 * localsearch.c
 */

/* local search methods */
typedef struct LS_data
{
    int key;
    void *obj;
}LS_data;


int compare_max(const void *, const void *);
int compare_min(const void *, const void *);
LS_data *construct_data(void *data,int key);
void free_LS(void *data);

/**
 * localsearch functions
 */

int k_l_move_general(Job** K_jobs, Job** L_jobs, partlist*m_k,partlist *m_l, solution *sol,int k, int l);
int local_search_machine_general_best( solution *sol,int lowerbound, int k, int l);
int local_search_machine_general_first( solution *sol,int lowerbound, int k, int l);

/**
 * compare_functions
 */
int sort_jobs(gconstpointer a, gconstpointer b);

/**
 * localsearch_wrappers
 */

void localsearch(solution *sol,int lowerbound);
void localsearch_SS(solution *sol,int lowerbound);
void localsearch_random(solution *sol,int lowerbound);
void localsearch_random_k(solution *sol,int lowerbound,int nb);

/**
 * wct.c
 */

/*Initialization and free memory for the problem*/
void wctproblem_init(wctproblem *problem);
void wctproblem_free(wctproblem *problem);

int compute_schedule(wctproblem *problem);

/*Computation functions*/
int compute_BPPC(wctproblem *problem);
int compute_lower_bound_BPPC(wctproblem *problem,wctdata *pd,int compute_lb);
int sequential_branching_BPPC(wctproblem *problem);
int create_branches(wctdata* pd,wctproblem* problem);
int check_integrality(wctdata *pd,int nmachine,int *result);
int build_lp_BPPC(wctdata* pd);
void make_pi_feasible(wctdata *pd);
int heur_colors_with_stable_sets( wctdata *pd );

int compute_PMCD(wctproblem *problem);
int compute_lower_bound_PMCD(wctproblem *problem, wctdata *pd);
int sequential_branching_PMCD(wctproblem *problem);
int create_branches_PMCD(wctdata *pd,wctproblem *problem);
int check_integrality_PMCD(wctdata *pd,int nmachines);
int build_lp_PMCD(wctdata *pd,int nmachines);
void make_pi_feasible_PMCD(wctdata* pd);
int compute_objective(wctdata *pd);


/*Help functions for branching*/
int insert_into_branching_heap(wctdata* pd,wctproblem* problem);
int skip_wctdata(wctdata* pd, wctproblem* problem);
int branching_msg(wctdata* pd, wctproblem* problem);
int recover_elist(wctdata *pd);
int collect_same_children(wctdata* pd);
int collect_diff_children(wctdata* pd);
void temporary_data_free(wctdata *pd);
int new_eindex(int v, int v1, int v2);
int prune_duplicated_sets (wctdata* pd);
int double2int(int *kpc_pi,int *scalef,const double *pi, int vcount);
double safe_lower_dbl(int numerator,int denominator);
int reset_root_pd(wctproblem *problem,wctdata *pd);
int reset_root_pd1(wctproblem *problem,wctdata *pd);
int add_newsets(wctdata *pd);
int add_to_init(wctproblem *problem,Scheduleset *sets,int nbsets);
int delete_to_bigcclasses(wctdata *pd, int capacity);
int add_some_maximal_stable_sets(wctdata *pd,int number);

int remove_finished_subtreebis(wctdata* child);

void make_pi_feasible(wctdata *pd);

/*Backup*/
int backup_wctdata(wctdata* pd, wctproblem* problem);

/**
 * wct.c
 */

int read_problem(char *f, int *njobs, int **duration, int **weights);

/**
 * heurdiving.c
 */

int heur_exec(wctproblem *problem,wctdata * pd,int *result);
void heur_init(wctdata * pd);
void heur_free(wctdata * pd);

#endif

