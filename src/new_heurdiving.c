#include "new_heurdiving.h"

int execute_CG_heur(wctproblem *problem, wctdata *pd) {
    int val = 0;
    GHashTable *G_0 = g_hash_table_new(g_direct_hash, g_direct_equal);
    GHashTable *part_sol = g_hash_table_new(g_direct_hash, g_direct_equal);
    GHashTable *rounded_sol = g_hash_table_new(g_direct_hash, g_direct_equal);

    /** Initialize the columns */
    for (int i = 0; i < pd->ccount; ++i) {
        g_hash_table_insert(G_0, pd->cclasses + i, pd->cclasses + i);
    }

    val = pure_diving(problem, pd, G_0, pd->rhs, pd->njobs, (double) pd->nmachines,
                      part_sol, rounded_sol);

    g_hash_table_destroy(G_0);
    g_hash_table_destroy(part_sol);
    g_hash_table_destroy(rounded_sol);
    return val;
}

int pure_diving(wctproblem *problem, wctdata *pd, GHashTable *G, double *a,
                int njobs, double K, GHashTable *part_sol, GHashTable *rounded_sol) {
    int val = 0;
    LP_data_CG_heur *data = CC_SAFE_MALLOC(1, LP_data_CG_heur);
    CCcheck_NULL_2(data, "Failed to allocate memory");
    LP_data_CG_heur_init(data, pd);
    gpointer key, value;
    GHashTableIter iter;


    /**step 1: Initliaze columns, RHS, K */
    data->a = update_rhs(a, njobs, rounded_sol);
    CCcheck_NULL_2(data->a, "Failed to update rhs");
    data->K = K;
    g_hash_table_foreach(rounded_sol, update_K_iterator, &(data->K));
    g_hash_table_foreach(rounded_sol, update_partsol_iterator, data);
    g_hash_table_foreach(part_sol, update_partsol_iterator, data);
    /** step 2: create the proper column pool */
    //g_hash_table_foreach(G,  update_column_iterator, &user_data);
    g_hash_table_iter_init(&iter, G);

    while (g_hash_table_iter_next(&iter, &key, &value)) {
        int add = 1;
        gpointer key_test, value_test;
        GHashTableIter iter_test;
        Scheduleset *set = (Scheduleset *) key;
        g_hash_table_iter_init(&iter_test, rounded_sol);

        while (g_hash_table_iter_next(&iter_test, &key_test, &value_test) && add) {
            Scheduleset *set_test = (Scheduleset *) key_test;

            for (int i = 0; i < set_test->count && add; ++i) {
                if (g_hash_table_lookup(set->table, GINT_TO_POINTER(i))) {
                    add = 0;
                }
            }
        }

        if (add) {
            g_hash_table_insert(data->G, key, value);
        }
    }

    /** build LP */
    val = build_lp_CG(data);
    CCcheck_val_2(val, "Failed to build the LP");
    /** Compute LP, find lambda_t, check if feasible or LP lower bound is bigger than upper bound*/
    compute_lower_bound_CG_heur(problem, data);

    /** F = g_ptr_array of all columns with non-integral value + check if partsol with G min F is primal  solution and record this solution*/
    if (data->status != LP_bound_computed) {
        goto CLEAN;
    }

    /** rounded_sol is set to 0 If F is empty return */

    /** choose heuristacally set of columns and set the value to one */

    /** call pure diving with new_G, etc */


CLEAN:
    LP_data_CG_heur_free(data);
    CC_IFFREE(data, LP_data_CG_heur);
    return val;
}

double *update_rhs(double *a, int njobs, GHashTable *rounded_sol) {
    double *new_a = CC_SAFE_MALLOC(njobs, double);
    CCcheck_NULL_3(new_a, "Failed to allocate memory to a_new");
    memcpy(new_a, a, njobs * sizeof(double));
    g_hash_table_foreach(rounded_sol, update_rhs_iterator, new_a);
CLEAN:
    return new_a;
}

void update_rhs_iterator(gpointer key, gpointer value, gpointer user_data) {
    Scheduleset *column = (Scheduleset *) key;
    double *x = (double *) value;
    double *a = (double *) user_data;

    for (int i = 0; i < column->count; ++i) {
        a[column->members[i]] -= (*x);
    }
}

void update_K_iterator(gpointer key, gpointer value, gpointer user_data) {
    double *new_K = (double *) user_data;
    double *x = (double *) value;
    *new_K -= *x;
}

void update_partsol_iterator(gpointer key, gpointer value, gpointer user_data) {
    LP_data_CG_heur *data = (LP_data_CG_heur *) user_data;
    GHashTable *new_partsol = data->partial_sol;

    g_hash_table_insert(new_partsol, key, value);

}

void test_column_iterator(gpointer key, gpointer value, gpointer user_data) {
    int *add = &(((pair_test_column *)user_data)->add);
    Scheduleset *set = ((pair_test_column *) user_data)->set;
    Scheduleset *test_set = (Scheduleset *) key;

    for (int i = 0; i < test_set->count && *add; ++i) {
        if (g_hash_table_lookup(set->table, GINT_TO_POINTER(i))) {
            *add = 0;
        }
    }
}

void add_column_iterator(gpointer key, gpointer value, gpointer user_data) {
    int *val = ((pair_build_LP *) user_data)->val;
    int j;
    Scheduleset *set = (Scheduleset *) key;
    LP_data_CG_heur *data = ((pair_build_LP *) user_data)->data;
    int *covered = ((pair_build_LP *) user_data)->covered;
    int *counter = &(((pair_build_LP *) user_data)->counter);
    double rounded_sol = ((pair_build_LP *) user_data)->rounded_sol;

    wctlp *LP = data->LP;
    *val = wctlp_addcol(LP, set->count + 1, set->members, data->coef, set->totwct,
                        rounded_sol, GRB_INFINITY, wctlp_CONT, NULL);
    set->id = data->nb_cols++;

    if (*val) {
        wctlp_printerrorcode(*val);
    }

    for (j = 0; j < set->count && *counter < data->pd->njobs; j++) {
        if (!covered[set->members[j]]) {
            covered[set->members[j]] = 1;
            (*counter)++;
        }
    }
}

void complete_RM_iterator(gpointer key, gpointer value, gpointer user_data) {
    int *val = ((pair_build_LP *) user_data)->val;
    int j;
    Scheduleset *set = (Scheduleset *) key;
    LP_data_CG_heur *data = ((pair_build_LP *) user_data)->data;
    wctlp *LP = data->LP;
    int *covered = ((pair_build_LP *) user_data)->covered;
    int *counter = &(((pair_build_LP *) user_data)->counter);
    double rhs = 1.0;
    *val = wctlp_addrow(LP, 0, (int *) NULL, (double *) NULL, wctlp_EQUAL, rhs,
                        (char *) NULL);
    double coef = 1.0;
    *val = wctlp_addcol(LP, 1 , &data->nb_rows, &coef, set->totwct, 0.0,
                        GRB_INFINITY, wctlp_CONT, NULL);
    data->nb_rows++;
    set->id = data->nb_cols++;

    if (*val) {
        wctlp_printerrorcode(*val);
    }

    for (j = 0; j < set->count && *counter < data->pd->njobs; j++) {
        if (!covered[set->members[j]]) {
            covered[set->members[j]] = 1;
            (*counter)++;
        }
    }

}

void print_columns_iterator(gpointer key, gpointer value, gpointer user_data) {
    Scheduleset *set = (Scheduleset *) key;

    for (int i = 0; i < set->count; ++i) {
        printf("%d ", set->members[i]);
    }

    printf("number of jobs =%d, Completion time %d, Cost %d\n", set->count,
           set->totweight, set->totwct);
}

int build_lp_CG(LP_data_CG_heur *data) {
    int val = 0;
    int i, j;
    wctdata *pd = data->pd;
    int njobs = pd->njobs;
    double *a = data->a;
    double K = data->K;
    int  *covered = (int *)NULL;
    GHashTableIter iter;
    gpointer key ;
    gpointer value;
    covered = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(covered, "Failed to allocate memory to covered");
    fill_int(covered, njobs, 0);
    int counter = 0;
    val = wctlp_init(&(data->LP), NULL);
    CCcheck_val_2(val, "wctlp_init failed");
    data->coef = (double *)CCutil_reallocrus(data->coef,
                 (njobs + 1) * sizeof(double));
    CCcheck_NULL_2(data->coef, "out of memory for coef");
    fill_dbl(data->coef, njobs + 1, 1.0);
    GHashTable *G = data->G;

    /** add the set covering constraints */
    for (i = 0; i < njobs; i++) {
        val = wctlp_addrow(data->LP, 0, (int *)NULL, (double *)NULL,
                           wctlp_GREATER_EQUAL, a[i],
                           (char *)NULL);
        CCcheck_val_2(val, "Failed wctlp_addrow");
        data->nb_rows++;
    }

    /** add the number of machines constraints */
    wctlp_addrow(data->LP, 0  , (int *)NULL, (double *) NULL, wctlp_LESS_EQUAL, K ,
                 NULL);
    data->nb_rows++;

    /** Write the columns */
    g_hash_table_iter_init(&iter, G);

    while (g_hash_table_iter_next(&iter, &key, &value)) {
        Scheduleset *set = (Scheduleset *) key;
        wctlp *LP = data->LP;
        val = wctlp_addcol(LP, set->count + 1, set->members, data->coef, set->totwct,
                           0.0, GRB_INFINITY, wctlp_CONT, NULL);
        CCcheck_val_2(val, "failed to add the column");
        set->id = data->nb_cols++;
        g_hash_table_iter_replace(&iter, &(set->id));

        for (j = 0; j < set->count && counter < data->pd->njobs; j++) {
            if (!covered[set->members[j]]) {
                covered[set->members[j]] = 1;
                (counter)++;
            }
        }
    }

    /** add the columns that are already selected */
    // g_hash_table_foreach(partsol, complete_RM_iterator, &user_data);
    // CCcheck_val_2(val, "Failed to the partial solution");


    data->pi = (double *)CCutil_reallocrus(data->pi,
                                           (data->nb_rows) * sizeof(double));
    CCcheck_NULL_2(data->pi, "Failed to allocate memory to data->pi");
    data->pi_in = CC_SAFE_MALLOC(data->nb_rows , double);
    CCcheck_NULL_2(data->pi_in, "Failed to allocate memory");
    fill_dbl(data->pi_in, data->nb_rows , 0.0);
    data->eta_in = 0.0;
    data->pi_out = CC_SAFE_MALLOC(data->nb_rows , double);
    CCcheck_NULL_2(data->pi_out, "Failed to allocate memory");
    data->eta_out = 0.0;
    fill_dbl(data->pi_out, data->nb_rows, 0.0);
    data->pi_sep = CC_SAFE_MALLOC(data->nb_rows, double);
    CCcheck_NULL_2(data->pi_sep, "Failed to allocate memory");
    data->subgradient_in = CC_SAFE_MALLOC(data->nb_rows, double);
    CCcheck_NULL_2(data->subgradient_in, "Failed to allocate memory");
    data->subgradient = CC_SAFE_MALLOC(data->nb_rows, double);
    CCcheck_NULL_2(data->subgradient, "Failed to allocate memory");
    data->rhs = CC_SAFE_MALLOC(data->nb_rows , double);
    CCcheck_NULL_2(data->rhs, "Failed to allocate memory");
    val = wctlp_get_rhs(data->LP, data->rhs);
    CCcheck_val_2(val, "Failed to get RHS");
    wctlp_write(data->LP, "test_heur.lp");
CLEAN:

    if (val) {
        wctlp_free(&(data->LP));
        CC_IFFREE(data->coef, double);
        CC_IFFREE(data->pi, double);
        CC_IFFREE(data->pi_in, double)
        CC_IFFREE(data->pi_out, double)
        CC_IFFREE(data->pi_sep, double)
        CC_IFFREE(data->subgradient, double)
        CC_IFFREE(data->subgradient_in, double)
        CC_IFFREE(data->rhs, double)
    }

    CC_IFFREE(covered, int);
    return val;
}

void LP_data_CG_heur_init(LP_data_CG_heur *data, wctdata *pd) {
    /*Initialization  of the LP*/
    data->LP = (wctlp *)NULL;
    data->x = (double *)NULL;
    data->coef = (double *) NULL;
    data->pi = (double *) NULL;
    /**init stab data */
    data->pi_in = (double *) NULL;
    data->pi_out = (double *) NULL;
    data->pi_sep = (double *) NULL;
    data->subgradient = (double *) NULL;
    data->subgradient_in = (double *) NULL;
    data->rhs = (double *) NULL;
    data->alpha = 0.8;
    data->nb_rows = 0;
    data->nb_cols = 0;
    data->G = g_hash_table_new(g_direct_hash, g_direct_equal);
    data->rounded_sol = g_hash_table_new(g_direct_hash, g_direct_equal);
    data->partial_sol = g_hash_table_new(g_direct_hash, g_direct_equal);
    data->nb_rounded_sol = 0;
    data->nb_partial_sol = 0;
    /** init data to the node */
    data->pd = pd;
}

void LP_data_CG_heur_free(LP_data_CG_heur *data) {
    wctlp_free(&(data->LP));
    CC_IFFREE(data->x, double);
    CC_IFFREE(data->coef, double);
    CC_IFFREE(data->pi, double);
    CC_IFFREE(data->pi_in, double);
    CC_IFFREE(data->pi_out, double);
    CC_IFFREE(data->subgradient, double);
    CC_IFFREE(data->subgradient_in, double);
    CC_IFFREE(data->rhs, double);
    g_hash_table_destroy(data->G);
    g_hash_table_destroy(data->rounded_sol);
    g_hash_table_destroy(data->partial_sol);
    CC_IFFREE(data->a, double);
}

void pair_build_LP_init(pair_build_LP *user_data, int *covered,
                        LP_data_CG_heur *data, int *val) {
    user_data->counter = 0;
    user_data->covered = covered;
    user_data->data = data;
    user_data->val = val;
    user_data->rounded_sol = 0.0;
}

int compute_objective_CG_heur(LP_data_CG_heur *data) {
    int val = 0;
    int i;

    data->LP_lower_bound_dual = .0;

    /** compute lower bound with the dual variables */
    for (i = 0; i < data->nb_rows; i++) {
        data->LP_lower_bound_dual += (double) data->pi[i] * data->rhs[i] ;
    }

    /** Get the LP lower bound and compute the lower bound of WCT */
    val = wctlp_objval(data->LP, &(data->LP_lower_bound));
    CCcheck_val_2(val, "wctlp_objval failed");
    data->lower_bound = ((int) ceil(data->LP_lower_bound_dual) < (int) ceil(
                             data->LP_lower_bound)) ? (int) ceil(data->LP_lower_bound_dual) : (int) ceil(
                            data->LP_lower_bound) ;

CLEAN:
    return val;
}

int compute_lower_bound_CG_heur(wctproblem *problem, LP_data_CG_heur *data) {
    int  j, val = 0;
    int iterations = 0;
    int break_while_loop = 1;
    int    nnonimprovements     = 0;
    int status = GRB_LOADED;
    double cur_cputime;
    wctdata *pd = data->pd;
    wctparms *parms = &(problem->parms);

    if (pd->status == infeasible) {
        goto CLEAN;
    }

    if (dbg_lvl() > 1) {
        printf("Starting compute_lower_bound with lb %d and ub %d at depth %d(id = %d, opt_track = %d)\n",
               pd->lower_bound, pd->upper_bound, pd->depth, pd->id, pd->opt_track);
    }

    CCutil_start_resume_time(&(problem->tot_lb_lp));


    if (!data->LP) {
        val = build_lp_CG(data);
        CCcheck_val(val, "build_lp failed");
    }

    /** Init alpha */
    switch (parms->stab_technique) {
    case stab_wentgnes:
        data->alpha = 0.8;
        break;

    case stab_dynamic:
        data->alpha = 0.0;
        break;

    case no_stab:
        break;
    }

    /** Compute LP relaxation */
    cur_cputime = CCutil_zeit();
    val = wctlp_optimize(data->LP, &status);
    CCcheck_val_2(val, "wctlp_optimize failed");
    cur_cputime = CCutil_zeit() - cur_cputime;

    if (dbg_lvl() > 1) {
        printf("Simplex took %f seconds.\n", CCutil_zeit() - cur_cputime);
        fflush(stdout);
    }

    switch (status) {
    case GRB_OPTIMAL:
        /** get the dual variables and make them feasible */
        val = wctlp_pi(data->LP, data->pi);
        CCcheck_val_2(val, "wctlp_pi failed");
        /** Compute the objective function */
        val = compute_objective_CG_heur(data);
        CCcheck_val_2(val, "Failed in compute_objective");
        memcpy(data->pi_out, data->pi, sizeof(double) * (data->nb_rows));
        data->eta_out = data->LP_lower_bound_dual;
        break;

    case GRB_INFEASIBLE:
        /** get the dual variables and make them feasible */
        val = wctlp_pi_inf(data->LP, data->pi);

        for (int i = 0; i < data->nb_rows; ++i) {
            CCcheck_val_2(val, "wctlp_pi failed");
        }

        break;
    }

    break_while_loop = 0;
    CCutil_suspend_timer(&(problem->tot_cputime));
    CCutil_resume_timer(&(problem->tot_cputime));

    while ((iterations < pd->maxiterations)
            && !break_while_loop
            && problem->tot_cputime.cum_zeit <= problem->parms.branching_cpu_limit) {

        /** Solve the pricing problem*/
        CCutil_start_resume_time(&problem->tot_pricing);

        switch (status) {
        case GRB_OPTIMAL:
            iterations++;

            if (iterations < pd->maxiterations) {
                switch (parms->stab_technique) {
                case stab_wentgnes:
                    solve_stab_CG_heur(problem, data);

                case stab_dynamic:

                case no_stab:
                    solve_pricing_CG_heur(data, parms);
                    break;
                }
            }

            break;

        case GRB_INFEASIBLE:
            val = solve_farkas_dbl_CG_heur(data);
            CCcheck_val_2(val, "Failed in solving farkas");
            break;
        }

        CCutil_suspend_timer(&problem->tot_pricing);

        for (j = 0; j < pd->nnewsets; j++) {
            val = wctlp_addcol(data->LP, pd->newsets[j].count + 1, pd->newsets[j].members,
                               data->coef, pd->newsets[j].totwct, 0.0, GRB_INFINITY, wctlp_CONT, NULL);
            CCcheck_val_2(val, "wctlp_addcol failed");
        }

        switch (status) {
        case GRB_OPTIMAL:
            switch (parms->stab_technique) {
            case stab_wentgnes:
            case stab_dynamic:
                break_while_loop = (CC_OURABS(data->eta_out - data->eta_in) < 0.00001);
                break;

            case no_stab:
                break_while_loop = (pd->nnewsets == 0 || nnonimprovements > 5);
                break;
            }

            break;

        case GRB_INFEASIBLE:
            break_while_loop = (pd->nnewsets == 0);
            break;
        }

        add_newsets_CG_heur(data);
        /** Compute LP relaxation */
        cur_cputime = CCutil_zeit();
        val = wctlp_optimize(data->LP, &status);
        CCcheck_val_2(val, "wctlp_optimize failed");
        cur_cputime = CCutil_zeit() - cur_cputime;

        if (dbg_lvl() > 1) {
            printf("Simplex took %f seconds.\n", CCutil_zeit() - cur_cputime);
            fflush(stdout);
        }

        switch (status) {
        case GRB_OPTIMAL:

            /** get the dual variables and make them feasible */
            /** Compute the objective function */

            if (pd->update) {
                val = wctlp_pi(data->LP, data->pi);
                CCcheck_val_2(val, "wctlp_pi failed");
                val = compute_objective_CG_heur(data);
                CCcheck_val_2(val, "Failed in compute_objective");
                memcpy(data->pi_out, data->pi, sizeof(double) * (data->nb_rows));
                data->eta_out = pd->LP_lower_bound_dual;
            }

            break;

        case GRB_INFEASIBLE:
            /** get the dual variables and make them feasible */
            val = wctlp_pi_inf(data->LP, data->pi);
            CCcheck_val_2(val, "wctlp_pi failed");
            break;
        }

        CCutil_suspend_timer(&(problem->tot_cputime));
        CCutil_resume_timer(&(problem->tot_cputime));
    }

    if (iterations < pd->maxiterations &&
            problem->tot_cputime.cum_zeit <= problem->parms.branching_cpu_limit) {
        switch (status) {
        case GRB_OPTIMAL:

            /** change status of problem */
            if (problem->status == no_sol) {
                problem->status = lp_feasible;
            }

            if (dbg_lvl() > 1) {
                printf("Found lb = %d (%f) upper_bound = %d (id= %d, iterations = %d,opt_track = %d).\n",
                       pd->lower_bound, pd->LP_lower_bound, pd->upper_bound,
                       pd->id, iterations, pd->opt_track);
            }

            data->x = CC_SAFE_MALLOC(pd->ccount, double);
            CCcheck_NULL_2(data->x, "Failed to allocate memory to pd->x");
            val = wctlp_x(data->LP, data->x, 0);
            CCcheck_val_2(val, "Failed in wctlp_x");
            data->status = LP_bound_computed;
            break;

        case GRB_INFEASIBLE:
            data->status = infeasible;
        }
    } else  {
        data->status = LP_bound_estimated;
    }

    if (dbg_lvl() > 1) {
        printf("iterations = %d\n", iterations);
    }

    fflush(stdout);
    problem->nb_generated_col += iterations;
    CCutil_suspend_timer(&(problem->tot_lb_lp));
CLEAN:
    return val;
}


int add_newsets_CG_heur(LP_data_CG_heur *data) {
    int val = 0;
    Scheduleset *tmpsets = (Scheduleset *) NULL;
    wctdata *pd = data->pd;
    int i;
    GHashTable *G = data->G;

    if (pd->nnewsets == 0) {
        return val;
    }

    if (pd->ccount + pd->nnewsets > pd->gallocated) {
        g_hash_table_destroy(G);
        G = g_hash_table_new(g_direct_hash, g_direct_equal);
        pd->gallocated *= 2;
        tmpsets = CC_SAFE_MALLOC(pd->gallocated, Scheduleset);
        CCcheck_NULL_2(tmpsets, "Failed to allocate memory to tmpsets");
        memcpy(tmpsets, pd->cclasses, pd->ccount * sizeof(Scheduleset));
        free(pd->cclasses);
        pd->cclasses = tmpsets;
        data->nb_cols = 0;

        for (i = 0; i < pd->ccount; ++i) {
            g_hash_table_insert(G, pd->cclasses + i, pd->cclasses + i);
            pd->cclasses[i].id = data->nb_cols++;
        }

        tmpsets = NULL;
    }

    memcpy(pd->cclasses + pd->ccount, pd->newsets,
           pd->nnewsets * sizeof(Scheduleset));

    for (i = pd->ccount; i < pd->ccount + pd->nnewsets; ++i) {
        g_hash_table_insert(G, pd->cclasses + i, pd->cclasses + i);
        pd->cclasses[i].id = data->nb_cols++;
    }

    pd->ccount += pd->nnewsets;

    for (i = pd->ccount; i < pd->gallocated; i++) {
        Scheduleset_init(pd->cclasses + i);
    }

CLEAN:

    if (val) {
        CC_IFFREE(pd->cclasses, Scheduleset);
    }

    CC_IFFREE(pd->newsets, Scheduleset);
    pd->nnewsets = 0;
    return val;
}


