#include "wct.h"

static int add_feasible_solution(wctproblem *problem, solution *new_sol);

/**
 * comparefunctions
 */

int compare_func1(gconstpointer a, gconstpointer b, void *user_data) {
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

/**
 * greedy constructions
 */

int random_rcl_assignment(Job *jobarray, int njobs, int nmachines,
                          solution *new_sol, GRand *rand_) {
    int i;
    double max;
    double min;
    double lb, ub;
    double temp_dbl;
    partlist *temp = (partlist *) NULL;
    Job *temp_job = (Job *) NULL;
    //GList *it = (GList *) NULL;
    GQueue *to_do_list = (GQueue *) NULL;
    temp = new_sol->part;
    to_do_list = g_queue_new();

    for (i = 0; i < njobs; ++i) {
        g_queue_push_tail(to_do_list, jobarray + i);
    }

    while (!g_queue_is_empty(to_do_list)) {
        temp_job = (Job *)to_do_list->head->data;
        max = ((double)temp[0].completiontime + (double)temp_job->processingime);
        min = max;
        GArray *rcl = g_array_new(FALSE, FALSE, sizeof(pair_job_machine));

        /** Compute min and max */
        for (i = 1; i < nmachines; ++i) {
            //for (it = to_do_list->head; it; it = it->next)
            //{
            //temp_job = (Job*)it->data;
            temp_dbl = (temp[i].completiontime + temp_job->processingime);

            if (max < temp_dbl) {
                max = temp_dbl;
            }

            if (min > temp_dbl) {
                min = temp_dbl;
            }

            //}
        }

        /** Compute RCL */
        pair_job_machine temp_job_machine;
        lb = min;
        ub = min + 0.25 * (max - lb);

        for (i = 0; i < nmachines; ++i) {
            //for (it = to_do_list->head; it; it = g_list_next(it))
            //{
            //temp_job = ((Job*)it->data);
            double g = ((double)temp[i].completiontime + (double)temp_job->processingime);

            if (lb <= g && g <= ub) {
                temp_job_machine.job = temp_job->job;
                temp_job_machine.machine = i;
                g_array_append_val(rcl, temp_job_machine);
            }

            //}
        }

        /** Choose uniformaly an assignment of a job to a machine */
        int a = g_rand_int_range(rand_, 0, rcl->len);
        int job = g_array_index(rcl, pair_job_machine, a).job;
        int machine = g_array_index(rcl, pair_job_machine, a).machine;
        partlist_insert(temp + machine, new_sol->vlist, jobarray + job);
        g_queue_pop_nth(to_do_list, g_queue_index(to_do_list, jobarray + job));
        g_array_free(rcl, TRUE);
    }

    g_queue_free(to_do_list);
    return 0;
}

int random_assignment(Job *jobarray, int njobs, int nmachines,
                      solution *new_sol, GRand *rand_) {
    int i, val = 0;
    double n;
    partlist *temp = (partlist *) NULL;
    Job *j = (Job *) NULL;
    GQueue *queue = (GQueue *) NULL;
    queue = g_queue_new();

    for (i = 0; i < nmachines; ++i) {
        g_queue_push_head(queue, new_sol->part + i);
    }

    for (i = 0; i < njobs; ++i) {
        j = jobarray + i;
        n = g_rand_double_range(rand_, 0.0, 1.0);

        if (n < 0.8) {
            temp = (partlist *) g_queue_pop_head(queue);
        } else if (n >= 0.8 && n < 0.95) {
            temp = (partlist *) g_queue_pop_nth(queue, 1);
        } else {
            temp = (partlist *)g_queue_pop_nth(queue, 2);
        }

        val = partlist_insert(temp, new_sol->vlist, j);
        CCcheck_val_2(val, "Failed in partlist_insert_order");
        g_queue_insert_sorted(queue, temp, compare_func1, NULL);
    }

CLEAN:
    g_queue_free(queue);
    return val;
}

int construct_wspt(Job *jobarray, int njobs, int  nmachines,
                   solution  *new_sol) {
    int i, val = 0;
    pmcheap *heap = (pmcheap *) NULL;
    partlist *temp = (partlist *) NULL;
    Job *j = (Job *) NULL;
    pmcheap_init(&heap, nmachines);
    CCcheck_NULL_2(heap, "Failed to initialize heap");

    for (i = nmachines - 1; i >= 0; i--) {
        pmcheap_insert(heap, new_sol->part[i].completiontime, new_sol->part + i);
    }

    for (i = 0; i < njobs; i++) {
        j = jobarray + i;
        temp = (partlist *) pmcheap_min(heap);
        val = partlist_insert(temp, new_sol->vlist, j);
        CCcheck_val_2(val, "Failed in partlist_insert_order");
        pmcheap_insert(heap, temp->completiontime, temp);
    }

CLEAN:

    if (val) {
        solution_free(new_sol);
    }

    pmcheap_free(heap);
    return 0;
}

/** Construct feasible solutions */

void update_bestschedule(wctproblem *problem, solution *new_sol) {
    if (new_sol == NULL) {
        return;
    }

    if (new_sol->totalweightcomptime < problem->global_upper_bound) {
        problem->global_upper_bound = new_sol->totalweightcomptime;
        problem->rel_error = (double)(problem->global_upper_bound -
                                      problem->global_lower_bound) / (problem->global_lower_bound);
        partlist_to_Scheduleset(new_sol->part, new_sol->nmachines, new_sol->njobs,
                                &(problem->bestschedule), &(problem->nbestschedule));
    }

    if (problem->global_upper_bound == problem->global_lower_bound) {
        problem->status = optimal;
    }
}

static int add_feasible_solution(wctproblem *problem, solution *new_sol) {
    int val = 0;
    wctdata *root_pd = &(problem->root_pd);
    SS *scatter_search = &(problem->scatter_search);
    solution_calc(new_sol, root_pd->jobarray);
    localsearch_wrap(new_sol, problem->global_lower_bound, 0);
    solution_unique(new_sol);

    if (!solution_in_pool(scatter_search, new_sol)) {
        add_solution_pool(scatter_search, new_sol);
        update_bestschedule(problem, new_sol);
        scatter_search->p->PSize++;

        if (root_pd->ccount == 0 && problem->parms.construct != 0) {
            update_Schedulesets(&root_pd->cclasses, &root_pd->ccount, problem->bestschedule,
                                problem->nbestschedule);
            root_pd->gallocated = root_pd->ccount;
        } else if (problem->parms.construct != 0) {
            partlist_to_Scheduleset(new_sol->part, new_sol->nmachines, new_sol->njobs,
                                    &(root_pd->newsets), &(root_pd->nnewsets));
            add_newsets(root_pd);
        }
    } else {
        solution_free(new_sol);
        CC_IFFREE(new_sol, solution);
    }

    return val;
}

solution *new_sol_init(int nmachines, int vcount) {
    int val = 0;
    solution *sol = (solution *) NULL;
    sol = CC_SAFE_MALLOC(1, solution);
    CCcheck_NULL_2(sol, "Failed to allocate memory")
    solution_init(sol);
    val = solution_alloc(sol, nmachines, vcount);
    CCcheck_val_2(val, "Failed in solution_alloc");
CLEAN:

    if (val) {
        solution_free(sol);
        CC_IFFREE(sol, solution);
    }

    return sol;
}

int construct_feasible_solutions(wctproblem *problem) {
    int val = 0;
    int iterations = 0;
    wctdata *pd = &(problem->root_pd);
    wctparms *parms = &(problem->parms) ;
    SS *scatter_search = &(problem->scatter_search);
    CCutil_timer *timer = &(problem->tot_scatter_search);
    GRand *rand1 = g_rand_new_with_seed(1984);
    GRand *rand2 = g_rand_new_with_seed(1654651);
    CCutil_start_timer(timer);
    val = SSproblem_definition(scatter_search, 10, 8,
                               parms->scatter_search_cpu_limit, parms->combine_method,
                               pd->njobs, pd->nmachines, pd->jobarray, problem->global_lower_bound);
    CCcheck_val_2(val, "Failed in SSproblem_definition");

    while (scatter_search->p->PSize < parms->nb_feas_sol) {
        iterations++;
        solution *new_sol = (solution *) NULL;
        new_sol = new_sol_init(pd->nmachines, pd->njobs);
        CCcheck_NULL(new_sol, "Failed to allocate")

        if (problem->status == no_sol) {
            construct_wspt(pd->jobarray, pd->njobs, pd->nmachines, new_sol);
            problem->status = feasible;
        } else {
            if (g_rand_boolean(rand1)) {
                random_assignment(pd->jobarray, pd->njobs, pd->nmachines, new_sol, rand2);
            } else {
                random_rcl_assignment(pd->jobarray, pd->njobs, pd->nmachines, new_sol, rand2);
            }
        }

        val = add_feasible_solution(problem, new_sol);
        CCcheck_val(val, "Failed in add_feasible_solution");
    }

    CCutil_suspend_timer(timer);
    printf("We needed %f seconds to construct %d solutions in %d iterations\n",
           timer->cum_zeit, parms->nb_feas_sol, iterations);
    printf("upperbound = %d, lowerbound = %d\n", problem->global_upper_bound,
           problem->global_lower_bound);
    CCutil_resume_timer(timer);

    if (parms->scatter_search) {
        SSCreate_refset(scatter_search);
        scatter_search->upperbound = problem->global_upper_bound;

        if (parms->combine_method) {
            SSrun_scatter_search(scatter_search, &(problem->tot_scatter_search));
        } else {
            SSrun_scatter_searchPR(scatter_search, & (problem->tot_scatter_search));
        }

        for (unsigned i = 0; i < scatter_search->rs->list1->len ; i++) {
            solution *new_sol = (solution *)g_ptr_array_index(scatter_search->rs->list1 ,
                                i);
            partlist_to_Scheduleset(new_sol->part, pd->nmachines, pd->njobs, &(pd->newsets),
                                    &(pd->nnewsets));
            add_newsets(pd);
        }
    }

CLEAN:
    g_rand_free(rand1);
    g_rand_free(rand2);
    CCutil_stop_timer(&(problem->tot_scatter_search), 0);
    return val;
}
