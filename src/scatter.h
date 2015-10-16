#ifndef SCATTER_H
#define SCATTER_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glib.h>
#include "datastructsol.h"

typedef struct _REFSET{ 
    int newsol;
    GQueue *list1;
    GQueue *list2;
} REFSET;

typedef struct _P{
    int PSize;
    GList *list;
    int lenght;
} P;

typedef struct _SS{
    REFSET *rs;
    P *p;
    int b1;
    int b2;
    
    Job **joblist;
    int njobs;
    int nmachines;
    int lowerbound; 
    int upperbound;

    double timelimit;
    GRand *random;
    int iter;
    int combine_method;
    
    enum {
        init   = 0,
        add    = 1,
        update = 3,
        opt    = 4
    } status;
} SS;

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

/**
 * Help functions
 */

int compute_fitness(SS *scatter_search);


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
void reinsert_fibheap(void *data, void *user_data);

int local_search_machine_move( solution* sol , int lowerbound);
int local_search_machine_exchange( solution *sol , int lowerbound);
int local_search_machine_2_0( solution *sol , int lowerbound);
int local_search_machine_2_1( solution *sol , int lowerbound);
int local_search_machine_2_2( solution *sol , int lowerbound);
int local_search_machine_general( solution *sol , int lowerbound   , int k , int l);
int local_search_machine_general_order( solution *sol , int lowerbound   , int k , int l);
void localsearch(solution *sol,int lowerbound);
void localsearch_SS(solution *sol,int lowerbound);
void localsearch_random(solution *sol,int lowerbound);
void localsearch_random_k(solution *sol,int lowerbound,int nb);

#endif

