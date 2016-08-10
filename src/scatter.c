#include <math.h>
#include <assert.h>
#include "wct.h"

static int select_parrent(GRand *Rand, int counter, int *sum, int totsum);
static int minimum_partition(solution *sol, int lowerbound);

int sort_vertices(const void *a, const void *b, void *data);

static int nodepair_ref_key(int v1, int v2);
static void inodepair_ref_key(int *v1, int *v2, int index);

void print_totalweightcomptime(void *data, void *user_data);
int order_totalweightcomptime(const void *a, const void *b);
int order_totalweightcomptime_list(const void *a, const void *b);

int move_scatter(Job *j, partlist *m_i);


/** compute row-index v1 and column-index v2 from array-index.*/
static void inodepair_ref_key(int *v1, int *v2, int index) {
    *v2 = (int) floor(sqrt(2 * ((double)index) + 0.25) - 0.5);
    *v1 = index - (*v2 * (*v2 + 1) / 2);
}

static int nodepair_ref_key(int v1, int v2) {
    /* We store only the elements of the upper right triangle within the
     njobs x njobs matrix. */
    assert(v1 <= v2);
    return v2 * (v2 + 1) / 2 + v1;
}

static int minimum_partition(solution *sol, int lowerbound) {
    int val = 0;

    for (int i = 1; i < sol->nmachines; ++i) {
        if (sol->part[i].completiontime - lowerbound <= 0) {
            if ((CC_OURABS(sol->part[i].completiontime - lowerbound) <
                    CC_OURABS(sol->part[val].completiontime - lowerbound))
                    || ((CC_OURABS(sol->part[i].completiontime - lowerbound) ==
                         CC_OURABS(sol->part[val].completiontime - lowerbound) &&
                         g_queue_get_length(sol->part[i].list) > g_queue_get_length(
                             sol->part[val].list)))) {
                val = i;
            }
        }
    }

    return val;
}
/*
For each functions
 */

void distance_min_max(void *data, void *user_data) {
    int k, l;
    min_max *temp = (min_max *)user_data;
    solution *new_sol = temp->new_sol;
    int n = temp->n;
    solution *sol = (solution *)data;
    int curval = 0;

    for (l = 0; l < n - 1; l++) {
        for (k = l + 1; k < n; k++) {
            if ((new_sol->vlist[l].part == new_sol->vlist[k].part
                    && sol->vlist[l].part != sol->vlist[k].part) ||
                    (new_sol->vlist[l].part != new_sol->vlist[k].part
                     && sol->vlist[l].part == sol->vlist[k].part)) {
                curval++;
            }

            if (curval > temp->min) {
                break;
            }
        }

        if (curval > temp->min) {
            break;
        }
    }

    if (curval <= temp->min) {
        temp->min = curval;
    }
}

void max_dist(void *data, void *user_data) {
    solution *sol = (solution *) data;
    solution *max = (solution *) user_data;

    if ((sol->dist > max->dist) || (sol->dist == max->dist
                                    && sol->totalweightcomptime < max->totalweightcomptime)) {
        max = (solution *) data;
    }
}

void free_sol(void *data, void *user_data) {
    solution *sol = (solution *)data;
    solution_free(sol);
    CC_IFFREE(sol, solution);
    sol = (solution *)user_data;
}

void assign_iter(void *data, void *user_data) {
    solution *sol = (solution *)data;
    int *iter = (int *)user_data;
    sol->iter = *iter;
}
void print_sol(void *data, void *user_data) {
    solution *sol = (solution *)data;
    int *i = (int *)user_data;
    printf("solution %d with wct %d \n", (*i)++, sol->totalweightcomptime);
    solution_print(sol);
}

void print_totalweightcomptime(void *data, void *user_data) {
    solution *sol  = (solution *)data;
    int *i = (int *)user_data;
    printf("solution %d with wct %d, distance %d \n", (*i)++,
           sol->totalweightcomptime,
           sol->dist);
}

void refset_dist(void *data, void *user_data) {
    SSrefset_distance((SS *)user_data, (solution *)data);
}

void for_each_comp_fitness(void *data, void *user_data) {
    solution *sol = (solution *) data;
    int *lowerbound = ((int *)user_data);
    sol->fitness = (float)(*lowerbound) / (sol->totalweightcomptime - *lowerbound);
}

/*
Compare functions
 */

int order_totalweightcomptime(const void *a, const void *b) {
    const int *aa = &(((const solution *)((const GPtrArray *)
                                          a)->pdata)->totalweightcomptime);
    const int *bb = &(((const solution *)((const GPtrArray *)
                                          b)->pdata)->totalweightcomptime);
    return *aa - *bb;
}
int order_distance(const void *a, const void *b) {
    const int *aa = &(((const solution *)((const GPtrArray *)a)->pdata)->dist);
    const int *bb = &(((const solution *)((const GPtrArray *)b)->pdata)->dist);
    return (*aa - *bb);
}

int order_totalweightcomptime_list(const void *a, const void *b) {
    const int *aa = &(((const solution *)a)->totalweightcomptime);
    const int *bb = &(((const solution *)b)->totalweightcomptime);
    return *aa - *bb;
}

int sort_vertices(const void *a, const void *b, void *data) {
    (void) data;
    const int *aa = &(((const Job *)a)->weight);
    const int *bb = &(((const Job *)b)->weight);
    return -(*aa - *bb);
}

/*
Init and free functions
 */
void SS_init(SS *problem, int b1, int b2, double timelimit) {
    problem->p         = (P *)NULL;
    problem->rs        = (REFSET *)NULL;
    problem->jobarray = (Job **)NULL;
    problem->b1        = b1;
    problem->b2        = b2;
    problem->timelimit = timelimit;
    problem->nmachines = 0;
    problem->njobs     = 0;
    problem->status    = init;
    problem->random    = (GRand *) NULL;
    problem->iter = 0;
    problem->combine_method = 1;
}

void SS_free(SS *problem) {
    REFSET_free(problem->rs);
    P_free(problem->p);
    CC_IFFREE(problem->rs, REFSET);
    CC_IFFREE(problem->p, P);
    CC_IFFREE(problem->jobarray, Job *);
    problem->b2        = 0;
    problem->b1        = 0;
    problem->timelimit = .0;
    problem->njobs        = 0;

    if (problem->random != NULL) {
        g_rand_free(problem->random);
    }
}

void REFSET_init(REFSET *rs) {
    if (rs) {
        rs->newsol        = 1;
        rs->list1         = g_ptr_array_new();
        rs->list2         = g_ptr_array_new();
    }
}

void REFSET_free(REFSET *rs) {
    if (rs) {
        rs->newsol = 0;
        g_ptr_array_foreach(rs->list1, free_sol, NULL);
        g_ptr_array_foreach(rs->list2, free_sol, NULL);
        g_ptr_array_free(rs->list1, TRUE);
        g_ptr_array_free(rs->list2 , TRUE);
        rs->newsol = 1;
        rs->list1 = (GPtrArray *) NULL;
        rs->list2 = (GPtrArray *) NULL;
    }
}

void P_init(P *p) {
    if (p) {
        p->PSize = 0;
        p->list = (GList *)NULL;
    }
}

void P_free(P *p) {
    if (p) {
        for (GList *it = p->list; it; it = it->next) {
            solution_free((solution *)it->data);
            CC_IFFREE(it->data, solution);
        }

        p->PSize = 0;
        g_list_free(p->list);
        P_init(p);
    }
}

void free_list2(REFSET *rs) {
    if (rs) {
        g_ptr_array_foreach(rs->list2, free_sol, NULL);
        g_ptr_array_free(rs->list2, TRUE);
        rs->list2 = (GPtrArray *)NULL;
    }
}

/*
Scatter Search functions
 */

int SSproblem_definition(
    SS *problem,
    int b1,
    int b2,
    double timelimit,
    int combine_method,
    int njobs,
    int nmachines,
    Job *jobarray,
    int lowerbound) {
    int i, val         = 0;
    REFSET *temp_rs = NULL;
    P *temp_p       = NULL;
    /*initialize scatter search data structure */
    SS_init(problem, b1, b2, timelimit);
    problem->combine_method = combine_method;
    problem->nmachines = nmachines;
    problem->njobs = njobs;
    problem->lowerbound = lowerbound;
    problem->upperbound = INT_MAX;
    /* Initialize pool */
    problem->p = CC_SAFE_MALLOC(1, P);
    CCcheck_NULL_2(problem->p, "Failed to allocate memory to problem->p");
    temp_p  = problem->p;
    P_init(temp_p);
    /* Initialize refset */
    problem->rs         = CC_SAFE_MALLOC(1, REFSET);
    CCcheck_NULL_2(problem->rs, "Failed to allocate memory to problem->rs");
    temp_rs             = problem->rs;
    REFSET_init(temp_rs);
    /* Initialize Jobarray of scatter search data structure */
    problem->jobarray = CC_SAFE_MALLOC(njobs, Job *);
    CCcheck_NULL_2(problem->jobarray, "Failed to allocate memory");

    for (i = 0; i < njobs; i++) {
        problem->jobarray[i] = jobarray + i ;
    }

    problem->random = g_rand_new_with_seed(48654642);
    CCcheck_NULL_2(problem->random, "Failed in g_rand_new_with_seed");
CLEAN:

    if (val) {
        SS_free(problem);
    }

    return val;
}


void add_solution_pool(SS *scatter_search, solution *new_sol) {
    P *p          = scatter_search->p;
    p->list = g_list_append(p->list, new_sol);
}

int add_solution_refset(SS *scatter_search) {
    int i;
    int val        = 0;
    REFSET *refset = scatter_search->rs;
    P *pool        = scatter_search->p;

    switch (scatter_search->status) {
    case init:
        pool->list = g_list_sort(pool->list, order_totalweightcomptime_list);

        for (i = 0; i < scatter_search->b1; ++i) {
            void *data = pool->list->data;
            pool->list = g_list_remove(pool->list, data);
            g_ptr_array_add(refset->list1, data);
        }

        scatter_search->status = add;
        break;

    case add:
        for (i = 0; i < scatter_search->b2; ++i) {
            void *data = maximum_distance(scatter_search);
            CCcheck_NULL_2(data, "Failed in maximum_distance");
            pool->list = g_list_remove(pool->list, data);
            update_distance(scatter_search, (solution *)data);
            g_ptr_array_add(refset->list2, data);
        }

        g_ptr_array_sort(refset->list2, order_distance);
        scatter_search->status = update;

        for (GList *it = pool->list; it; it = it->next) {
            solution_free((solution *)it->data);
            CC_IFFREE(it->data, solution);
        }

        g_list_free(pool->list);
        pool->list = (GList *) NULL;
        break;

    case update:
        diversification_update(scatter_search);
        break;

    case opt:
        break;
    }

CLEAN:
    return val;
}

int solution_in_pool(SS *scatter_search, solution *new_sol) {
    int val = 1;
    int njobs = new_sol->njobs;
    GList *it = scatter_search->p->list;

    if (it == (GList *) NULL || g_list_length(it) == 0) {
        val = 0;
        return val;
    }

    do {
        int i = 0;

        for (i = 0; i < njobs; ++i) {
            if (((solution *)it->data)->perm[i] != new_sol->perm[i]) {
                break;
            }
        }

        if (i == njobs) {
            return val;
        }

        it = it->next;
    } while (it != (GList *) NULL);

    val = 0;
    return val;
}

int solution_in_refset(SS *scatter_search, solution *new_sol) {
    int val    = 1;
    GPtrArray *list1 =  scatter_search->rs->list1;
    int njobs = new_sol->njobs;
    int i = 0;
    guint j = 0;

    if (list1->len == 0) {
        val = 0;
        return val;
    }

    do {
        for (i = 0; i < njobs; ++i) {
            solution *sol = (solution *) g_ptr_array_index(list1, j);

            if (sol->perm[i] != new_sol->perm[i]) {
                break;
            }
        }

        if (i == njobs) {
            return val;
        }

        j++;
    } while (j < list1->len);

    GPtrArray *list2 = scatter_search->rs->list2;

    if (list2->len == 0) {
        val = 0;
        return val;
    }

    j = 0;

    do {
        for (i = 0; i < njobs; ++i) {
            solution *sol = (solution *) g_ptr_array_index(list2, j);

            if ((sol)->perm[i] != new_sol->perm[i]) {
                break;
            }
        }

        if (i == njobs) {
            return val;
        }

        j++;
    } while (j < list2->len);

    val = 0;
    return val;
}

void print_pool(SS *scatter_search) {
    P *pool = scatter_search->p;
    int i = 0;
    g_list_foreach(pool->list, print_sol, &i);
}

void print_pool_totalweightcomptime(SS *scatter_search) {
    P *pool = scatter_search->p;
    int i = 0;
    g_list_foreach(pool->list, print_totalweightcomptime, &i);
}


void print_refset(SS *scatter_search) {
    REFSET *rs = scatter_search->rs;
    int i = 0;
    g_ptr_array_foreach(rs->list1, print_sol, &i);
    i = 0;
    g_ptr_array_foreach(rs->list2, print_sol, &i);
}

void print_refset_totalweightcomptime(SS *scatter_search) {
    REFSET *rs = scatter_search->rs;
    int i = 0;
    g_ptr_array_foreach(rs->list1, print_totalweightcomptime, &i);
    g_ptr_array_foreach(rs->list2, print_totalweightcomptime, &i);
}

void print_pool_n(SS *scatter_search, int n) {
    P *pool = scatter_search->p;
    solution *sol = (solution *)NULL;

    if ((guint)n > g_list_length(pool->list)) {
        printf("n is too big\n");
    } else {
        for (int i = 0; i < n; ++i) {
            sol = ((solution *)g_list_nth(pool->list, i)->data);
            printf("solution %d with wct %d\n", i, sol->totalweightcomptime);
            solution_print(sol);
        }
    }
}

void print_list1(SS *scatter_search) {
    REFSET *refset = scatter_search->rs;
    int i = 0;
    g_ptr_array_foreach(refset->list1, print_sol, &i);
}

void print_distance(SS *scatter_search) {
    P *pool = scatter_search->p;
    GList *it = (GList *)NULL;
    int i = 0;

    for (it = pool->list; it; it = it->next) {
        printf("Distance to refset = %d and solution %d\n",
               ((solution *)it->data)->dist, i++);
    }
}


int SSrefset_distance(SS *scatter_search, solution *new_sol) {
    int val = 0;
    REFSET *refset = scatter_search->rs;
    min_max temp = {scatter_search->njobs, INT_MAX, new_sol};

    if (refset->list1->len ==  0) {
        printf("We can't compute the distance between pool and refset\n");
        val = 1;
        return val;
    }

    g_ptr_array_foreach(refset->list1, distance_min_max, &temp);
    g_ptr_array_foreach(refset->list2, distance_min_max, &temp);
    new_sol->dist = temp.min;
    return val;
}



void *maximum_distance(SS *scatter_search) {
    solution *val = (solution *)NULL;
    P *pool = scatter_search->p;
    val = (solution *)pool->list->data;
    g_list_foreach(pool->list, max_dist, val);
    return val;
}

int update_distance(SS *scatter_search, solution *sol) {
    int l, k, val = 0;
    GList *it = (GList *)NULL;
    P *pool = scatter_search->p;

    if (pool->list == (GList *)NULL) {
        printf("We can't update the distances. The pool is empty!!!!\n");
        val = 1;
        goto CLEAN;
    }

    for (it = pool->list; it; it = g_list_next(it)) {
        int curval = 0;

        for (l = 0; l < scatter_search->njobs - 1; l++) {
            for (k = l + 1; k < scatter_search->njobs; k++) {
                if ((((solution *)it->data)->vlist[l].part == ((solution *)
                        it->data)->vlist[k].part && sol->vlist[l].part != sol->vlist[k].part) ||
                        (((solution *)it->data)->vlist[l].part != ((solution *)
                                it->data)->vlist[k].part && sol->vlist[l].part == sol->vlist[k].part)) {
                    curval++;
                }

                if (curval > ((solution *)it->data)->dist) {
                    break;
                }
            }

            if (curval > ((solution *)it->data)->dist) {
                break;
            }
        }

        if (curval <= ((solution *)it->data)->dist) {
            ((solution *)it->data)->dist = curval;
        }
    }

CLEAN:
    return val;
}



int SSCreate_refset(SS *scatter_search) {
    int val   = 0;
    P *pool   = scatter_search->p;
    add_solution_refset(scatter_search);
    g_list_foreach(pool->list, refset_dist, scatter_search);
    add_solution_refset(scatter_search);
    return val;
}

int compute_fitness(SS *scatter_search) {
    int val = 0;
    REFSET *refset = scatter_search->rs;
    int lowerbound = scatter_search->lowerbound;
    g_ptr_array_foreach(refset->list1, for_each_comp_fitness, &lowerbound);
    g_ptr_array_foreach(refset->list2, for_each_comp_fitness, &lowerbound);
    return val;
}

int SSrun_scatter_search(SS *scatter_search, CCutil_timer *timer) {
    guint i;
    REFSET *refset = scatter_search->rs;
    int flag, val = 0;
    int nbnew_sol = 0;
    int nbjobs = scatter_search->njobs;
    int nmachines = scatter_search->nmachines;
    int nb_noimprovements = 0;
    int *totsubset = (int *) NULL;
    double limit = scatter_search->timelimit;
    totsubset = CC_SAFE_MALLOC(scatter_search->b1 + 1, int);
    CCcheck_NULL_2(totsubset, "Failed tot allocate memory to tot subset");

    if (refset->list1->len == 0) {
        printf("We can't run scatter search, refset is empty\n");
        val = 1;
        goto CLEAN;
    }

    printf("PMCombination method\n");
    CCutil_suspend_timer(timer);
    CCutil_resume_timer(timer);

    while (refset->newsol && scatter_search->status != opt &&
            timer->cum_zeit < limit) {
        GPtrArray *list = g_ptr_array_new();

        for (i = 0; i < refset->list1->len; ++i) {
            g_ptr_array_add(list, g_ptr_array_index(refset->list1, i));
        }

        for (i = 0; i < refset->list2->len; ++i) {
            g_ptr_array_add(list, g_ptr_array_index(refset->list2, i));
        }

        compute_fitness(scatter_search);
        refset->newsol = 0;

        if (scatter_search->iter > 0) {
            int best = scatter_search->upperbound;
            CCutil_suspend_timer(timer);
            CCutil_resume_timer(timer);
            double rel_error = ((double)best - (double)scatter_search->lowerbound) /
                               (double)scatter_search->lowerbound;
            printf("iteration %d with best wct %d and rel error %f, number of new solutions %d, time = %f\n",
                   scatter_search->iter, best, rel_error, nbnew_sol, timer->cum_zeit);
        }

        nbnew_sol = 0;
        k_subset_init(list->len, 2, totsubset, &flag);

        while (flag && scatter_search->status != opt) {
            solution *new_sol = CC_SAFE_MALLOC(1, solution);
            solution_init(new_sol);
            solution_alloc(new_sol, nmachines, nbjobs);
            int rval = 1;
            rval = combine_PM(scatter_search, list, totsubset, 2, new_sol);

            if (!rval) {
                solution_calc(new_sol, *(scatter_search->jobarray));
                new_sol->iter = scatter_search->iter + 1;
                localsearch_random_k(new_sol, scatter_search->lowerbound, 2);
                solution_unique(new_sol);

                if (!dynamic_update(scatter_search, list, new_sol)) {
                    if (new_sol->totalweightcomptime < scatter_search->upperbound) {
                        scatter_search->upperbound = new_sol->totalweightcomptime;
                    }

                    if (new_sol->totalweightcomptime == scatter_search->lowerbound) {
                        scatter_search->status = opt;
                        printf("Found optimal with SS with subset generation 1\n");
                    }

                    nbnew_sol++;
                } else {
                    solution_free(new_sol);
                    CC_IFFREE(new_sol, solution);
                }
            } else {
                solution_free(new_sol);
                CC_IFFREE(new_sol, solution);
            }

            k_subset_lex_successor(list->len, 2, totsubset, &flag);
        }

        k_subset_init(list->len, 3, totsubset, &flag);
        solution *sol = (solution *)NULL;
        sol = (solution *) g_ptr_array_index(list, 0);

        while (flag && scatter_search->status != opt) {
            solution *new_sol = CC_SAFE_MALLOC(1, solution);
            solution_init(new_sol);
            solution_alloc(new_sol, nmachines, nbjobs);
            int rval = 1;
            rval = combine_PM(scatter_search, list, totsubset, 3,
                              new_sol);  //Dell'Amico et al.

            if (!rval) {
                solution_calc(new_sol, *(scatter_search->jobarray));
                new_sol->iter = scatter_search->iter + 1;
                localsearch_random_k(new_sol, scatter_search->lowerbound, 2);
                solution_unique(new_sol);

                if (!dynamic_update(scatter_search, list, new_sol)) {
                    if (new_sol->totalweightcomptime < scatter_search->upperbound) {
                        scatter_search->upperbound = new_sol->totalweightcomptime;
                    }

                    if (new_sol->totalweightcomptime == scatter_search->lowerbound) {
                        scatter_search->status = opt;
                        printf("Found optimal with SS\n");
                    }

                    nbnew_sol++;
                } else {
                    solution_free(new_sol);
                    CC_IFFREE(new_sol, solution);
                }
            } else {
                solution_free(new_sol);
                CC_IFFREE(new_sol, solution);
            }

            k_subset_lex_successor(list->len, 3, totsubset, &flag);
            sol = (solution *) g_ptr_array_index(list, totsubset[1] - 1);
        }

        k_subset_init(list->len, 4, totsubset, &flag);
        sol = (solution *) g_ptr_array_index(list, 0);
        solution *temp_sol = (solution *) g_ptr_array_index(list, 1);

        while (flag && scatter_search->status != opt
                && sol->totalweightcomptime <= scatter_search->upperbound
                && temp_sol->totalweightcomptime <= scatter_search->upperbound) {
            solution *new_sol = CC_SAFE_MALLOC(1, solution);
            solution_init(new_sol);
            solution_alloc(new_sol, nmachines, nbjobs);
            int rval = 1;
            rval = combine_PM(scatter_search, list, totsubset, 4, new_sol);

            if (!rval) {
                solution_calc(new_sol, *(scatter_search->jobarray));
                new_sol->iter = scatter_search->iter + 1;
                localsearch_random_k(new_sol, scatter_search->lowerbound, 2);
                solution_unique(new_sol);

                if (!dynamic_update(scatter_search, list, new_sol)) {
                    if (new_sol->totalweightcomptime < scatter_search->upperbound) {
                        scatter_search->upperbound = new_sol->totalweightcomptime;
                    }

                    if (new_sol->totalweightcomptime == scatter_search->lowerbound) {
                        scatter_search->status = opt;
                        printf("Found optimal with SS\n");
                    }

                    nbnew_sol++;
                } else {
                    solution_free(new_sol);
                    CC_IFFREE(new_sol, solution);
                }
            } else {
                solution_free(new_sol);
                CC_IFFREE(new_sol, solution);
            }

            k_subset_lex_successor(list->len, 4, totsubset, &flag);
            sol = (solution *) g_ptr_array_index(list, totsubset[1] - 1);
            temp_sol = (solution *)g_ptr_array_index(list, totsubset[2] - 1);
        }

        for (i = 0; i < refset->list2->len + 1; ++i) {
            totsubset[i] = i;
        }

        for (i = 5; i < refset->list2->len + 1
                && scatter_search->status != opt; i++) {
            solution *new_sol = CC_SAFE_MALLOC(1, solution);
            solution_init(new_sol);
            solution_alloc(new_sol, nmachines, nbjobs);
            int rval = 1;
            rval = combine_PM(scatter_search, list, totsubset, i, new_sol);

            if (!rval) {
                solution_calc(new_sol, *(scatter_search->jobarray));
                new_sol->iter = scatter_search->iter + 1;
                localsearch_random_k(new_sol, scatter_search->lowerbound, 2);
                solution_unique(new_sol);

                if (!dynamic_update(scatter_search, list, new_sol)) {
                    if (new_sol->totalweightcomptime < scatter_search->upperbound) {
                        scatter_search->upperbound = new_sol->totalweightcomptime;
                    }

                    if (new_sol->totalweightcomptime == scatter_search->lowerbound) {
                        scatter_search->status = opt;
                        printf("Found optimal with SS\n");
                    }

                    nbnew_sol++;
                } else {
                    solution_free(new_sol);
                    CC_IFFREE(new_sol, solution);
                }
            } else {
                solution_free(new_sol);
                CC_IFFREE(new_sol, solution);
            }
        }

        scatter_search->iter++;

        if (nbnew_sol == 0) {
            nb_noimprovements++;
        }

        if (refset->newsol == 0 && scatter_search->iter < 10
                && nb_noimprovements < 5) {
            add_solution_refset(scatter_search);
            refset->newsol = 1;
        }

        g_ptr_array_free(list, TRUE);
    }

CLEAN:
    CC_IFFREE(totsubset, int);
    return val;
}

int moveSS(Job *j, partlist *m_j, partlist *m_i) {
    int nb_job = j->job;
    return j->processingime * (m_j->sumweights[nb_job] - m_i->sumweights[nb_job])
           + j->weight * (m_j->sumtimes[nb_job] - m_i->sumtimes[nb_job] -
                          j->processingime);
}

int combinePathRelinking(SS *scatter_search, GPtrArray *array, int *subsetsol,
                         int *found) {
    int i, val = 0;
    int njobs = scatter_search->njobs;
    int nmachines = scatter_search->nmachines;
    int dist = 0;
    int counter = 0;
    Job **jobarray = scatter_search->jobarray;
    REFSET *refset = scatter_search->rs;
    GList *list = (GList *) NULL;
    GList *it = (GList *) NULL;
    GList *pool = (GList *) NULL;
    solution *prev = (solution *) NULL;
    solution sol_g;
    solution sol_c;
    solution *sol1 = (solution *)g_ptr_array_index(array, subsetsol[1] - 1);
    solution *sol2 = (solution *)g_ptr_array_index(array, subsetsol[2] - 1);

    if (sol1->totalweightcomptime < sol2->totalweightcomptime) {
        val = solution_copy(&sol_g, *sol1);
        val = solution_copy(&sol_c, *sol2);
    } else {
        val = solution_copy(&sol_c, *sol1);
        val = solution_copy(&sol_g, *sol2);
    }

    for (i = 0; i < njobs; ++i) {
        partlist *temp1 = sol_g.vlist[i].part;
        partlist *temp2 = sol_c.vlist[i].part;

        if (temp1->key != temp2->key) {
            dist++;
            list = g_list_append(list, jobarray[i]);
        }
    }

    prev = &(sol_c);

    while (dist > 0) {
        Job *minjob = (Job *) NULL;
        partlist *minmach = (partlist *) NULL;
        solution *sol = CC_SAFE_MALLOC(1, solution);
        solution_init(sol);
        solution_alloc(sol, nmachines, njobs);
        solution_update(sol, *prev);
        int max = INT_MIN;
        int temp;

        for (it = list; it; it = it->next) {
            Job *j = (Job *)it->data;
            partlist *mach_c = sol->vlist[j->job].part;
            partlist *mach_g = sol->part + sol_g.vlist[j->job].part->key;
            temp = moveSS(j, mach_c, mach_g);

            if (temp > max) {
                minjob = j;
                minmach = mach_g;
                max = temp;
            }
        }

        partlist_move_order(minmach, sol->vlist, minjob, njobs);
        sol->totalweightcomptime -= max;
        pool = g_list_append(pool, sol);
        list = g_list_remove(list, minjob);
        dist--;
        prev = sol;
    }

    counter = 0;
    GList *l = pool;

    while (l != NULL) {
        GList *next = l->next;

        if (counter % 10 == 0) {
            solution *sol = (solution *)l->data;
            localsearch_wrap(sol, scatter_search->lowerbound, 0);
            solution_unique(sol);
        } else {
            solution *sol = (solution *)l->data;
            solution_free(sol);
            CC_IFFREE(sol, solution);
            pool = g_list_delete_link(pool, l);
        }

        counter++;
        l = next;
    }

    l = pool;

    while (l != NULL) {
        solution *last1 = (solution *) g_ptr_array_index(refset->list1,
                          refset->list1->len - 1);
        solution *last2 = (solution *) g_ptr_array_index(refset->list2,  0);
        GList *next = l->next;
        solution *new_sol = (solution *)l->data;
        int not_in_refset = !solution_in_refset(scatter_search, new_sol);

        if (new_sol->totalweightcomptime <  last1->totalweightcomptime &&
                not_in_refset) {
            new_sol->dist = 0;
            refset->newsol = 1;
            new_sol->iter = scatter_search->iter + 1;
            g_ptr_array_remove(refset->list1, last1);
            g_ptr_array_remove(array, last1);
            g_ptr_array_add(refset->list1, new_sol);
            g_ptr_array_add(array, new_sol);
            g_ptr_array_sort(refset->list1, order_totalweightcomptime);
            g_ptr_array_sort(array, order_totalweightcomptime);
            solution_free(last1);
            CC_IFFREE(last1, solution);
            *found = 1;
        } else if (not_in_refset) {
            SSrefset_distance(scatter_search, new_sol);

            if (new_sol->dist > last2->dist) {
                refset->newsol = 1;
                new_sol->iter = scatter_search->iter + 1;
                g_ptr_array_remove(refset->list2, last2);
                g_ptr_array_remove(array, last2);
                g_ptr_array_add(refset->list2, new_sol);
                g_ptr_array_add(array, new_sol);
                g_ptr_array_sort(refset->list2, order_distance);
                g_ptr_array_sort(array, order_distance);
                solution_free(last2);
                CC_IFFREE(last2, solution);
                *found = 1;
            } else {
                solution_free(new_sol);
                CC_IFFREE(new_sol, solution);
            }
        } else {
            solution_free(new_sol);
            CC_IFFREE(new_sol, solution);
        }

        pool = g_list_delete_link(pool, l);
        l = next;
    }

    g_list_free(pool);
    solution_free(&sol_g);
    solution_free(&sol_c);
    return val;
}

int SSrun_scatter_searchPR(SS *scatter_search, CCutil_timer *timer) {
    guint i;
    int flag, val         = 0;
    REFSET *refset    = scatter_search->rs;
    int nbnew_sol = 0;
    int nb_noimprovements = 0;
    double limit = scatter_search->timelimit;
    int *totsubset = (int *) NULL;
    totsubset = CC_SAFE_MALLOC(3, int);
    CCcheck_NULL_2(totsubset, "Failed tot allocate memory to tot subset");

    if (refset->list1->len == 0) {
        printf("We can't run scatter search, refset is empty\n");
        val = 1;
        goto CLEAN;
    }

    CCutil_suspend_timer(timer);
    CCutil_resume_timer(timer);

    while (refset->newsol && scatter_search->status != opt &&
            timer->cum_zeit < limit) {
        refset->newsol = 0;
        GPtrArray *array = g_ptr_array_new();

        for (i = 0; i < refset->list1->len; ++i) {
            g_ptr_array_add(array, g_ptr_array_index(refset->list1, i));
        }

        for (i = 0; i < refset->list2->len; ++i) {
            g_ptr_array_add(array, g_ptr_array_index(refset->list2, i));
        }

        refset->newsol = 0;

        if (scatter_search->iter > 0) {
            int best = scatter_search->upperbound;
            CCutil_suspend_timer(timer);
            CCutil_resume_timer(timer);
            double rel_error = ((double)best - (double)scatter_search->lowerbound) /
                               (double)scatter_search->lowerbound;
            printf("iteration %d with best wct %d and rel error %f, number of new solutions %d, time %f\n",
                   scatter_search->iter, best, rel_error, nbnew_sol, timer->cum_zeit);
        }

        nbnew_sol = 0;
        k_subset_init(array->len, 2, totsubset, &flag);

        while (flag && scatter_search->status != opt && timer->cum_zeit < limit) {
            int rval = 0;
            solution *sol1 = (solution *)g_ptr_array_index(array, totsubset[1] - 1);
            solution *sol2 = (solution *)g_ptr_array_index(array, totsubset[2] - 1);

            if (sol1->iter >= scatter_search->iter || sol2->iter >= scatter_search->iter) {
                combinePathRelinking(scatter_search, array, totsubset, &rval);
            }

            if (rval) {
                solution *temp = (solution *) g_ptr_array_index(array, 0);

                if (temp->totalweightcomptime < scatter_search->upperbound) {
                    scatter_search->upperbound = temp->totalweightcomptime;
                }

                if (temp->totalweightcomptime == scatter_search->lowerbound) {
                    scatter_search->status = opt;
                }

                nbnew_sol++;
            }

            k_subset_lex_successor(array->len, 2, totsubset, &flag);
            CCutil_suspend_timer(timer);
            CCutil_resume_timer(timer);
        }

        if (nbnew_sol == 0) {
            nb_noimprovements++;
        }

        scatter_search->iter++;
        CCutil_suspend_timer(timer);
        CCutil_resume_timer(timer);

        if (refset->newsol == 0 && scatter_search->iter < scatter_search->njobs &&
                nb_noimprovements < 8  && timer->cum_zeit < limit) {
            add_solution_refset(scatter_search);
            refset->newsol = 1;
        }

        g_ptr_array_free(array, TRUE);
        CCutil_suspend_timer(timer);
        CCutil_resume_timer(timer);
    }

CLEAN:
    CC_IFFREE(totsubset, int);
    return val;
}

/* functions for constructing new solutions*/

static int select_parrent(GRand *Rand, int nbelements, int *sum, int totsum) {
    int val = -1;
    int counter = 0;

    for (int i = 0; i < nbelements; ++i) {
        if (sum[i] > 0) {
            counter++;
        }
    }

    int temp = counter - 1;
    int *parents = (int *)NULL;

    if (counter > 0) {
        parents = CC_SAFE_MALLOC(counter, int);
    } else {
        goto CLEAN;
    }

    if (totsum <= 0) {
        goto CLEAN;
    }

    for (int i = 0; i < nbelements; ++i) {
        if (sum[i] > 0) {
            parents[temp--] = i;
        }
    }

    val = parents[g_rand_int_range(Rand, 0, counter)];
CLEAN:
    CC_IFFREE(parents, int);
    return val;
}

int combine_GPX(SS *scatter_search, GPtrArray *queue, int *subsetsol,
                int nbelements, solution *new_sol) {
    int k, val    = 1;
    int nmachines = ((solution *)g_ptr_array_index(scatter_search->rs->list1,
                     0))->nmachines;
    int njobs = ((solution *)    g_ptr_array_index(scatter_search->rs->list2,
                 0))->njobs;
    int part, i, j;
    int totsum = 0;
    int *sum = (int *) NULL;
    solution *sol = (solution *)NULL;
    solution *temp_sol = (solution *) NULL;
    sum = CC_SAFE_MALLOC(nbelements, int);
    CCcheck_NULL_2(sum, "Failed to allocate memory");

    for (i = 0; i < nbelements; i++) {
        sum[i] = njobs;
        totsum += njobs;
    }

    for (k = 1;  k <= nbelements && val; k++) {
        temp_sol  = (solution *) g_ptr_array_index(queue, subsetsol[k] - 1);

        if (temp_sol->iter >= scatter_search->iter) {
            val = 0;
        }
    }

    if (val) {
        goto CLEAN;
    }

    sol = CC_SAFE_MALLOC(nbelements, solution);
    CCcheck_NULL_2(sol, "Failed to allocate memory");

    for (i = 0; i  < nbelements; i++) {
        solution_copy(sol + i, *((solution *)g_ptr_array_index(queue,
                                 subsetsol[i + 1] - 1)));
    }

    val = solution_alloc(new_sol, nmachines,  njobs);
    CCcheck_val_2(val, "Failed in solution_alloc");
    new_sol->nmachines = 0;
    new_sol->njobs = 0;

    while (new_sol->nmachines < nmachines && totsum > 0) {
        i = select_parrent(scatter_search->random, nbelements, sum, totsum);

        if (i != -1) {
            part = minimum_partition(sol + i, scatter_search->lowerbound);

            for (GList *it = (sol + i)->part[part].list->head; it; it = it->next) {
                partlist_insert(&new_sol->part[new_sol->nmachines], new_sol->vlist,
                                (Job *)it->data);
            }

            for (GList *it = new_sol->part[new_sol->nmachines].list->head; it ;
                    it = it->next) {
                for (j = 0; j < nbelements; j++) {
                    //partlist_delete_custom( ( sol + j )->vlist, ( ( Job * )it->data ), njobs );
                }
            }

            for (j = 0; j < nbelements; j++) {
                sum[j] -= g_queue_get_length(new_sol->part[new_sol->nmachines].list);
            }

            totsum -= nbelements * g_queue_get_length(
                          new_sol->part[new_sol->nmachines].list);
            new_sol->njobs += g_queue_get_length(new_sol->part[new_sol->nmachines].list);
            new_sol->nmachines++;
        } else {
            solution_print(new_sol);

            for (j = 0; j < nbelements; j++) {
                printf("%d ", sum[j]);
            }

            printf("\n");
            printf("totsum = %d\n", totsum);
        }
    }

    if (new_sol->njobs != njobs) {
        pmcheap *heap = (pmcheap *)NULL;
        partlist *temp = (partlist *) NULL;
        GQueue *list = g_queue_new();

        for (i = 0; i < njobs; ++i) {
            if (sol->part[i].list != NULL) {
                for (GList *it = sol->part[i].list->head; it; it = it->next) {
                    g_queue_push_head(list, it->data);
                }
            }
        }

        g_queue_sort(list, sort_vertices, NULL);
        pmcheap_init(&heap, new_sol->nmachines);

        for (i = 0; i < new_sol->nmachines; ++i) {
            val = pmcheap_insert(heap, new_sol->part[i].completiontime, new_sol->part + i);

            if (val) {
                printf("Failed at pmcheap_insert in %s at %d\n", __FILE__, __LINE__);
                g_queue_free(list);
                pmcheap_free(heap);
                goto CLEAN;
            }
        }

        GList *it = list->head;

        while (it) {
            temp = (partlist *) pmcheap_min(heap);
            Job *job = (Job *)it->data;
            partlist_insert(temp, new_sol->vlist, job);
            new_sol->njobs++;
            pmcheap_insert(heap, temp->completiontime, temp);
            it = it->next;
        }

        pmcheap_free(heap);
        g_queue_free(list);
    }

    solution_calc(new_sol, *(scatter_search->jobarray));
CLEAN:

    if (val) {
        solution_free(new_sol);
    }

    CC_IFFREE(sum, int);

    if (sol) {
        for (i = 0; i < nbelements; i++) {
            solution_free(sol + i);
        }

        CC_IFFREE(sol, solution);
    }

    return val;
}

int move_scatter(Job *j, partlist *m_i) {
    int nb_job = j->job;
    return j->processingime * (- m_i->sumweights[nb_job])
           + j->weight * (- m_i->sumtimes[nb_job] - j->processingime);
}

int combine_PM(SS *scatter_search, GPtrArray *array, int *subsetsol,
               int nbelements, solution *new_sol) {
    int i, j, k, winner, val = 1;
    int nbsubsets = ((scatter_search->njobs  + 1) * scatter_search->njobs) / 2;
    int njobs = scatter_search->njobs, nmachines = scatter_search->nmachines, it,
        count, step, first;
    partlist *temp = (partlist *)NULL;
    solution *sol = (solution *)NULL;
    float *fitness = (float *)NULL;
    float *cumulfitness = (float *) NULL;
    pmcheap *heap = (pmcheap *) NULL;

    // if ( list->len >= ( guint )nbelements ) {
    //     solution *tempsol = ( solution * )g_ptr_array_index(list, 0);
    //     njobs = tempsol->njobs;
    //     nmachines = tempsol->nmachines;
    // } else {
    //     printf( "Error in combine solution.\n" );
    //     val = 1;
    //     goto CLEAN;
    // }

    for (k = 1;  k <= nbelements && val; k++) {
        sol  = (solution *) g_ptr_array_index(array, subsetsol[k] - 1);

        if (sol->iter >= scatter_search->iter) {
            val = 0;
        }
    }

    if (val) {
        goto CLEAN;
    }

    fitness = CC_SAFE_MALLOC(nbsubsets, float);
    CCcheck_NULL_2(fitness, "Failed to allocate memory");
    fill_float(fitness, nbsubsets, .0);
    cumulfitness = CC_SAFE_MALLOC(nbsubsets, float);
    CCcheck_NULL_2(cumulfitness, "Failed to allocate memory");
    fill_float(cumulfitness, nbsubsets, .0);
    new_sol->njobs = 0;
    val = pmcheap_init(&heap, nmachines);
    CCcheck_val_2(val, "Failed in pmcheap_init");

    for (i = 0; i < nmachines; ++i) {
        pmcheap_insert(heap, new_sol->part[i].completiontime, new_sol->part + i);
    }

    for (i = 0; i < scatter_search->njobs; i++) {
        for (j = i ; j < scatter_search->njobs; j++) {
            int key = nodepair_ref_key(i, j);
            fitness[key] = 0;

            if (i != j) {
                for (k = 1; k <= nbelements; k++) {
                    sol  = (solution *)g_ptr_array_index(array, subsetsol[k] - 1);

                    if (sol->vlist[i].part == sol->vlist[j].part) {
                        fitness[key] += sol->fitness;
                    }
                }
            }
        }
    }

    for (i = 1; i < nbsubsets; i++) {
        cumulfitness[i] = fitness[i] + cumulfitness[i - 1];
    }

    while (cumulfitness[nbsubsets - 1] > 0) {
        const float f = g_rand_double_range(scatter_search->random, .0,
                                            cumulfitness[nbsubsets - 1]);
        count = nbsubsets;
        first = 0;

        while (count > 0) {
            it = first;
            step = count / 2;
            it += step;

            if (cumulfitness[it] < f) {
                first = ++it;
                count -= step + 1;
            } else {
                count = step;
            }
        }

        winner = first;
        inodepair_ref_key(&i, &j, winner);
        Job *node1 = *(scatter_search->jobarray + i);
        Job *node2 = *(scatter_search->jobarray + j);
        /*int max = move_scatter(node1, new_sol->part) + move_scatter(node2, new_sol->part) - node1->processingime*node2->weight;
        temp = new_sol->part;
        for ( i = 1; i < nmachines; ++i)
        {
            int a = move_scatter(node1, new_sol->part + i) + move_scatter(node2, new_sol->part + i) - node1->processingime*node2->weight;
            if (a > max)
            {
                temp = new_sol->part + i;
                max = a;
            }
        } */
        temp = (partlist *) pmcheap_min(heap);
        val = partlist_insert(temp, new_sol->vlist, node1);
        CCcheck_val_2(val, "Failed partlist_insert");
        val = partlist_insert(temp, new_sol->vlist, node2);
        CCcheck_val_2(val, "Failed partlist_insert");
        new_sol->njobs += 2;
        pmcheap_insert(heap, temp->completiontime, temp);

        for (i = 0; i < scatter_search->njobs; i++) {
            if (i <= node1->job) {
                fitness[nodepair_ref_key(i, node1->job)] = .0;
            } else {
                fitness[nodepair_ref_key(node1->job, i)] = .0;
            }

            if (i <= node2->job) {
                fitness[nodepair_ref_key(i, node2->job)] = .0;
            } else {
                fitness[nodepair_ref_key(node2->job, i)] = .0;
            }
        }

        for (i = 1; i < nbsubsets; i++) {
            cumulfitness[i] = fitness[i] + cumulfitness[i - 1];
        }
    }

    if (new_sol->njobs != njobs) {
        for (i = 0; i < njobs; i++) {
            if (new_sol->vlist[scatter_search->jobarray[i]->job].part == NULL) {
                Job *node = scatter_search->jobarray[i];
                /*int max = move_scatter(node, new_sol->part);
                temp = new_sol->part;
                for ( j = 1; j < nmachines; ++j)
                {
                    int a = move_scatter(node, new_sol->part + j);
                    if (a < max)
                    {
                        temp = new_sol->part + j;
                        max = a;
                    }
                } */
                temp = (partlist *) pmcheap_min(heap);
                val = partlist_insert(temp, new_sol->vlist, node);
                CCcheck_val_2(val, "Failed in partlist order")
                new_sol->njobs++;
                pmcheap_insert(heap, temp->completiontime, temp);
            }
        }
    }

    new_sol->fitness = (double)scatter_search->lowerbound /
                       (new_sol->totalweightcomptime - scatter_search->lowerbound);
CLEAN:

    if (val) {
        solution_free(new_sol);
    }

    pmcheap_free(heap);
    CC_IFFREE(fitness, float);
    CC_IFFREE(cumulfitness, float);
    return val;
}



// int static_update( SS *scatter_search )
// {
//     int val = 0;
//     GList *it = ( GList * )NULL;
//     REFSET *refset = scatter_search->rs;
//     P *pool = scatter_search->p;

//     for ( it = refset->list1->head; it; it = it->next ) {
//         ( ( solution * )it->data )->iter = 1;
//     }

//     for ( it = refset->list2->head; it; it = it->next ) {
//         ( ( solution * )it->data )->iter = 1;
//     }

//     for ( it = pool->list; it; it = it->next ) {
//         SSrefset_distance( scatter_search, it->data );
//     }

//     refset->newsol = 0;
//     it = pool->list;

//     while ( it ) {
//         solution *sol = ( solution * )it->data;
//         solution *last = g_queue_peek_tail( refset->list1 );
//         solution *last1 = g_queue_peek_tail( refset->list2 );
//         int not_in_refset = !solution_in_refset( scatter_search, sol );

//         if ( sol->totalweightcomptime < last->totalweightcomptime && not_in_refset ) {
//             refset->newsol = 1;
//             void *data = g_queue_pop_tail(refset->list1);
//             solution_free( data );
//             CC_IFFREE( data, solution );
//             it = pool->list = g_list_remove( pool->list, sol );

//             if ( g_list_length( pool->list ) != 0 ) {
//                 val = update_distance( scatter_search, sol );
//                 CCcheck_val_2( val, "Failed in update_distance" );
//                 g_queue_insert_sorted( refset->list1, sol, order_totalweightcomptime, NULL );
//             }
//         } else if ( sol->dist > last1->dist ) {
//             refset->newsol = 1;
//             solution *data = g_queue_pop_tail( refset->list2 );
//             solution_free( data );
//             CC_IFFREE( data, solution );
//             it = pool->list = g_list_remove( pool->list, sol );

//             if ( g_list_length( pool->list ) != 0 ) {
//                 val = update_distance( scatter_search, sol );
//                 CCcheck_val_2( val, "Failed in update_distance" );
//                 g_queue_insert_sorted( refset->list2, sol, order_distance, NULL );
//             }
//         } else {
//             it = it->next;
//         }
//     }

//     printf( "update\n" );
// CLEAN:
//     return val;
// }

int dynamic_update(SS *scatter_search, GPtrArray *list, solution *new_sol) {
    int val = 1;
    REFSET *refset = scatter_search->rs;
    solution *last1 = (solution *) g_ptr_array_index(refset->list1,
                      refset->list1->len - 1);
    assert(last1 != NULL);
    solution *last2 = (solution *) g_ptr_array_index(refset->list2,  0);
    assert(last2 != NULL);
    int not_in_refset = !solution_in_refset(scatter_search, new_sol);

    if (new_sol->totalweightcomptime <  last1->totalweightcomptime &&
            not_in_refset) {
        refset->newsol = 1;
        new_sol->dist = 0;
        g_ptr_array_remove(refset->list1, last1);
        g_ptr_array_remove(list, last1);
        g_ptr_array_add(refset->list1, new_sol);
        g_ptr_array_add(list, new_sol);
        g_ptr_array_sort(refset->list1, order_totalweightcomptime);
        g_ptr_array_sort(list, order_totalweightcomptime);
        solution_free(last1);
        CC_IFFREE(last1, solution);
        val = 0;
        return val;
    } else if (not_in_refset) {
        SSrefset_distance(scatter_search, new_sol);

        if (new_sol->dist > last2->dist) {
            refset->newsol = 1;
            g_ptr_array_remove(refset->list2, last2);
            g_ptr_array_remove(list, last2);
            g_ptr_array_add(refset->list2, new_sol);
            g_ptr_array_add(list, new_sol);
            g_ptr_array_sort(refset->list2, order_distance);
            g_ptr_array_sort(list, order_distance);
            solution_free(last2);
            CC_IFFREE(last2, solution);
            val = 0;
            return val;
        }
    }

    return val;
}

int diversification_update(SS *scatter_search) {
    int  val = 0;
    int njobs = scatter_search->njobs;
    int nmachines = scatter_search->nmachines;
    Job *jobarray = *(scatter_search->jobarray);
    GRand *rand_ = scatter_search->random;
    REFSET *refset = scatter_search->rs;
    free_list2(refset);
    refset->list2 = (GPtrArray *) NULL;
    refset->list2 = g_ptr_array_new();

    while (refset->list2->len < (guint)scatter_search->b2
            && scatter_search->status != opt) {
        solution *new_sol = (solution *)NULL;
        new_sol = CC_SAFE_MALLOC(1, solution);
        CCcheck_NULL(new_sol, "Failed to allocate memory to new_sol");
        solution_init(new_sol);
        solution_alloc(new_sol, nmachines, njobs);

        if (g_rand_boolean(scatter_search->random)) {
            random_rcl_assignment(jobarray, njobs, nmachines, new_sol, rand_);
        } else {
            random_assignment(jobarray, njobs, nmachines, new_sol, rand_);
        }

        solution_calc(new_sol, jobarray);
        localsearch_random_k(new_sol, scatter_search->lowerbound, 3);
        solution_unique(new_sol);

        if (!solution_in_refset(scatter_search, new_sol)) {
            if (new_sol->totalweightcomptime < scatter_search->upperbound) {
                scatter_search->upperbound = new_sol->totalweightcomptime;
            }

            if (scatter_search->upperbound == scatter_search->lowerbound) {
                scatter_search->status = opt;
            }

            solution *data = (solution *) g_ptr_array_index(refset->list1,
                             refset->list1->len - 1);

            if (new_sol->totalweightcomptime < data->totalweightcomptime) {
                g_ptr_array_remove(refset->list1, data);
                SSrefset_distance(scatter_search, data);
                g_ptr_array_add(refset->list1, new_sol);
                g_ptr_array_add(refset->list2, data);
            } else {
                SSrefset_distance(scatter_search, new_sol);
                g_ptr_array_add(refset->list2, new_sol);
            }
        } else {
            solution_free(new_sol);
            CC_IFFREE(new_sol, solution);
        }
    }

    if (refset->list2->len == 0) {
        printf("NO new solutions\n");
        val = 1;
        goto CLEAN;
    }

    g_ptr_array_sort(refset->list2, order_distance);
    g_ptr_array_sort(refset->list1, order_totalweightcomptime);
    g_ptr_array_foreach(refset->list1, assign_iter, &scatter_search->iter);
    g_ptr_array_foreach(refset->list2, assign_iter, &scatter_search->iter);
CLEAN:
    return val;
}
