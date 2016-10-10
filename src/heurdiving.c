#include <assert.h>
#include <math.h>
#include "heurdiving.h"



int heur_exec(wctproblem *problem, wctdata *pd, int *result) {
    int i, status, val = 0;
    int searchbound;
    double objval;
    int cutoff;
    int lperror;
    GQueue *branchingvars = g_queue_new();
    int *selectedvars = (int *) NULL;
    int *tabulist = (int *) NULL;
    int *discrepancies = (int *) NULL;
    int nlpcands;
    int startnlpcands;
    //int maxdepth;
    int divedepth;
    int maxdivedepth;
    WCTbranchcand branchcand;
    WCTbranchcand_init(&branchcand);
    double *value = (double *)NULL;
    *result = DELEYAD;
    heur_diving *heur = &(pd->heur_data);
    CCcheck_NULL_2(heur, "heur is empty");
    discrepancies = CC_SAFE_MALLOC(heur->maxdiscdepth, int);
    CCcheck_NULL_2(discrepancies, "Failed to allocate memory")

    for (i = 0; i < heur->maxdiscdepth; ++i) {
        discrepancies[i] = 0;
    }

    tabulist = CC_SAFE_MALLOC(heur->maxdiscrepancy , int);
    CCcheck_NULL_2(tabulist, "Failed to allocate memory")

    for (i = 0; i < heur->maxdiscrepancy; ++i) {
        tabulist[i] = 0;
    }

    selectedvars = CC_SAFE_MALLOC(heur->maxdiscdepth, int);
    CCcheck_NULL_2(selectedvars, "Failed to allocate memory")
    val = branchcandlp(pd, &branchcand);
    CCcheck_val_2(val, "Failed in branchcandlp");

    if (branchcand.nlpcands == 0) {
        int success;
        val = constructsolution(pd, problem->parms.nmachines, &success);
        CCcheck_val_2(val, "Failed at constructsolution");
        goto CLEAN;
    }

    nlpcands = branchcand.nlpcands;
    *result = DIDNOTFIND;
    lperror = 0;
    cutoff = 0;
    divedepth = 0;
    startnlpcands = nlpcands;
    status = GRB_OPTIMAL;
    val =  wctlp_objval(pd->LP, &objval);
    CCcheck_val_2(val, "Failed in wctlp_objval");
    /** Adjusement here? */
    pd->LP_lower_bound_heur = pd->partial_sol + pd->LP_lower_bound_dual;
    objval = pd->LP_lower_bound_heur;
    searchbound = pd->lower_bound;
    maxdivedepth = problem->parms.nmachines + 1;

    if (dbg_lvl() > 0) {
        printf("executing primal heuristic: depth=%d, %d fractionals, searchbound=%d",
               divedepth, nlpcands, searchbound);
    }

    while (!lperror && !cutoff && status == GRB_OPTIMAL && nlpcands > 0
            && *result != FOUNDSOL
            && (divedepth < 10 || nlpcands <= startnlpcands - divedepth / 2
                || (divedepth < maxdivedepth))
            && (divedepth >= heur->maxdiscdepth
                || discrepancies[divedepth] <= heur->maxdiscrepancy)
            && (problem->tot_cputime.cum_zeit < problem->parms.branching_cpu_limit)
            && (pd->status != LP_bound_estimated)) {
        int bestcand;
        double bestcandsol;
        int bestcandmayround;
        int backtracked;
        int farkaspricing;
        divedepth++;
        bestcand = -1;
        bestcandmayround = 1;
        val = heur_divingselect_var(pd, tabulist, heur->maxdiscrepancy, &bestcand,
                                    &bestcandmayround, &branchcand);
        CCcheck_val_2(val, "Failed at heur_divingselect_var");

        if (bestcand == -1) {
            printf("No variable for diving could be selected, diving aborted\n");
            break;
        }

        assert(bestcand != -1);
        GRBgetdblattrelement(pd->LP->model, GRB_DBL_ATTR_X, bestcand, &bestcandsol);

        if (divedepth - 1 < heur->maxdiscdepth) {
            selectedvars[divedepth - 1] = bestcand;
        }

        g_queue_push_tail(branchingvars, GINT_TO_POINTER(bestcand));

        if (bestcandmayround) {
            int success;
            val = constructsolution(pd, problem->parms.nmachines, &success);
            CCcheck_val(val, "Failed in constructsolution");

            if (success) {
                if (dbg_lvl() > 0) {
                    printf(" -> solution was feasible and good enough\n");
                }

                *result = FOUNDSOL;
                break;
            }
        }

        if (dbg_lvl() == 0) {
            printf("  dive %d/%d,  var <%d>, sol=%g\n",
                   divedepth, maxdivedepth, bestcand, bestcandsol);
        }

        val = adjustLP_ceil(pd, bestcand, bestcandsol);
        CCcheck_val_2(val, "Failed to adjust LP");
        backtracked = 0;
        farkaspricing = 0;

        do {
            cutoff = (((status > GRB_OPTIMAL)
                       || (objval > (double)pd->upper_bound + lp_int_tolerance())
                       || (nlpcands == 0 && objval > (double)pd->upper_bound + lp_int_tolerance())));

            if (!cutoff || backtracked || farkaspricing) {
                /** adjustment here? */
                val = compute_lower_bound(problem, pd);
                wctlp_write(pd->LP, "testheur.lp");
                pd->LP_lower_bound_heur = pd->LP_lower_bound_dual + pd->partial_sol;

                if (val) {
                    printf("Failed at heur_compute_lower_bound\n");
                    lperror = 1;
                    break;
                }

                val = wctlp_status(pd->LP, &status);
                CCcheck_val(val, "Failed at wctlp_status");
                cutoff = ((status > GRB_OPTIMAL)  ? 1 : 0);
            }

            if (status == GRB_INFEASIBLE && !farkaspricing && !backtracked) {
                //if (dbg_lvl() > 0) {
                printf("*** infeasibility detected at level %d - perform Farkas pricing\n",
                       divedepth);
                //}
                farkaspricing = 1;
            } else {
                farkaspricing = 0;
            }

            CCutil_suspend_timer(&(problem->tot_cputime));
            CCutil_resume_timer(&(problem->tot_cputime));

            if (cutoff && !backtracked && heur->backtrack
                    && *result != FOUNDSOL &&
                    problem->tot_cputime.cum_zeit < problem->parms.branching_cpu_limit) {
                int discrepencie = 0;

                do {
                    int var = GPOINTER_TO_INT(g_queue_pop_tail(branchingvars));
                    adjustLP_floor(pd, var);
                    --divedepth;

                    if (divedepth < heur->maxdiscdepth) {
                        discrepencie = (discrepancies[divedepth] >= heur->maxdiscrepancy);
                    }
                } while (divedepth > 0 && (divedepth >= heur->maxdiscdepth
                                           || discrepencie));

                assert(divedepth < heur->maxdiscdepth);

                if (discrepancies[divedepth] < heur->maxdiscrepancy) {
                    tabulist[discrepancies[divedepth]] = selectedvars[divedepth];
                }

                ++discrepancies[divedepth];

                for (i = divedepth + 1; i < heur->maxdiscdepth; ++i) {
                    discrepancies[i] = discrepancies[divedepth];
                }

                backtracked = 1;
            } else {
                backtracked = 0;
            }
        } while (backtracked || farkaspricing);

        CCutil_suspend_timer(&(problem->tot_cputime));
        CCutil_resume_timer(&(problem->tot_cputime));

        if (!lperror && !cutoff && status == GRB_OPTIMAL &&
                problem->tot_cputime.cum_zeit < problem->parms.branching_cpu_limit) {
            val = wctlp_objval(pd->LP , &objval);
            objval = pd->LP_lower_bound_dual + pd->partial_sol;
            WCTbranchcand_free(&branchcand);
            val = branchcandlp(pd, &branchcand);
            CCcheck_val_2(val, "Failed at branchcandlp");
            nlpcands = branchcand.nlpcands;
        }

        //if (dbg_lvl() > 0) {
        printf(" -> status=%d, objval=%f/%d, nfrac=%d, depth = %d\n", status, objval,
               searchbound,
               nlpcands, divedepth);
        //}
        CCutil_suspend_timer(&(problem->tot_cputime));
        CCutil_resume_timer(&(problem->tot_cputime));
    }

    if (nlpcands == 0 && !lperror && !cutoff && status == GRB_OPTIMAL) {
        int success;
        val = constructsolution(pd, problem->parms.nmachines, &success);
        CCcheck_val_2(val, "Failed at constructsolution");

        if (success) {
            if (dbg_lvl() > 0) {
                printf(" -> solution was feasible and good enough\n");
            }

            *result = FOUNDSOL;
        }
    }

CLEAN:
    CC_IFFREE(value, double);
    CC_IFFREE(discrepancies, int)
    CC_IFFREE(tabulist, int)
    CC_IFFREE(selectedvars, int)
    WCTbranchcand_free(&branchcand);
    g_queue_free(branchingvars);
    return val;
}

int branchcandlp(wctdata *pd, WCTbranchcand *branchcand) {
    int i, status, val = 0;
    wctlp *LP = pd->LP;
    val = wctlp_status(LP, &status);
    CCcheck_val_2(val, "Failed in wctlp_status");

    if (!(status == GRB_OPTIMAL || status == GRB_UNBOUNDED)) {
        fprintf(stderr, "LP status is not GRB_OPTIMAL or GRB_UNBOUNDED\n");
        val = 1;
        goto CLEAN;
    }

    if (status == GRB_UNBOUNDED) {
        if (dbg_lvl() > 0) {
            printf("LP is unbounded -> no branching \n");
        }

        goto CLEAN;
    }

    if (dbg_lvl() > 0) {
        printf(" -> recalculating LP branching candidates\n");
    }

    int ncols;
    int insertpos;
    double primsol;
    double frac;
    double lb;
    val = GRBgetintattr(LP->model, GRB_INT_ATTR_NUMVARS, &ncols);
    CHECK_VAL_GRB2(val, "Failed to get the numbers of columns", LP->env);
    val = WCTbranchcand_alloc(branchcand, ncols);
    CCcheck_val_2(val, "Failed to allocate WCTbranchcand");

    for (i = 0; i < ncols; ++i) {
        val = GRBgetdblattrelement(LP->model, GRB_DBL_ATTR_X, i, &primsol);
        CHECK_VAL_GRB2(val, "Failed to get primsol", LP->env);
        val = GRBgetdblattrelement(LP->model, GRB_DBL_ATTR_LB, i, &lb);
        CHECK_VAL_GRB2(val, "Failed to get GRB_DBL_ATTR_LB", LP->env);

        if (lb >= 0.5) {
            continue;
        }

        frac = primsol - floor(primsol + EPSILON);

        if (frac <= EPSILON) {
            continue;
        }

        insertpos = branchcand->nlpcands;
        branchcand->lpcands[insertpos] = i;
        branchcand->lpcandssol[insertpos] = primsol;
        branchcand->lpcandsfrac[insertpos] = frac;
        branchcand->nlpcands++;
    }

CLEAN:

    if (val) {
        WCTbranchcand_free(branchcand);
    }

    return val;
}

void WCTbranchcand_free(WCTbranchcand *branchcand) {
    if (branchcand) {
        CC_IFFREE(branchcand->lpcands, int);
        CC_IFFREE(branchcand->lpcandssol, double);
        CC_IFFREE(branchcand->lpcandsfrac, double);
    }

    WCTbranchcand_init(branchcand);
}

void WCTbranchcand_init(WCTbranchcand *branchcand) {
    if (branchcand) {
        branchcand->lpcandssize = 0;
        branchcand->nlpcands = 0;
        branchcand->lpcands = (int *) NULL;
        branchcand->lpcandssol = (double *) NULL;
        branchcand->lpcandsfrac = (double *) NULL;
    }
}

int WCTbranchcand_alloc(WCTbranchcand *branchcand, int ncol) {
    int val = 0;

    if (ncol) {
        branchcand->lpcandssize = ncol;
        branchcand->lpcands = CC_SAFE_MALLOC(ncol, int);
        CCcheck_NULL_2(branchcand->lpcands, "Failed to allocate memory");
        branchcand->lpcandssol = CC_SAFE_MALLOC(ncol, double);
        CCcheck_NULL_2(branchcand->lpcandssol, "Failed to allocate memory");
        branchcand->lpcandsfrac = CC_SAFE_MALLOC(ncol, double);
        CCcheck_NULL_2(branchcand->lpcandsfrac, "Failed to allocate memory");
    } else {
        val = 1;
    }

CLEAN:
    return val;
}

int heur_divingselect_var(wctdata *pd, int *tabulist, int tabulistsize,
                          int *bestcand, int *bestcandmayround, WCTbranchcand *branchcand) {
    int val = 0;
    double bestobjgain;
    double bestfrac;
    int bestcandmayrounddown;
    int bestcandmayroundup;
    int i;
    assert(bestcand != NULL);
    assert(bestcand != NULL);
    assert(bestcandmayround != NULL);
    assert(branchcand->lpcands != NULL);
    assert(branchcand->lpcandsfrac != NULL);
    assert(branchcand->lpcandssol != NULL);
    bestcandmayrounddown = 1;
    bestcandmayroundup = 1;
    bestobjgain = DBL_MAX;
    bestfrac = DBL_MAX;

    for (i = 0; i < branchcand->nlpcands; ++i) {
        int var;
        double frac;
        double obj;
        double lb;
        double ub;
        int j;
        var = branchcand->lpcands[i];
        frac = branchcand->lpcandsfrac[i];
        val = GRBgetdblattrelement(pd->LP->model, GRB_DBL_ATTR_OBJ, var, &obj);
        CHECK_VAL_GRB2(val, "Failed at GRB_DBL_ATTR_OBJ", pd->LP->env);
        val = GRBgetdblattrelement(pd->LP->model, GRB_DBL_ATTR_LB, var, &lb);
        CHECK_VAL_GRB2(val, "Failed at GRB_DBL_ATTR_LB", pd->LP->env);
        val = GRBgetdblattrelement(pd->LP->model, GRB_DBL_ATTR_UB, var, &ub);
        CHECK_VAL_GRB2(val, "Failed at GRB_DBL_ATTR_UB", pd->LP->env);
        int mayrounddown = !(lb >= 1.0) ? 1 : 0;
        int mayroundup = !(ub <= .0) ? 1 : 0;

        for (j = 0; j < tabulistsize; ++j) {
            if (tabulist[j] == var) {
                break;
            }

            if (j < tabulistsize) {
                continue;
            }
        }

        if (mayrounddown || mayroundup) {
            if (bestcandmayrounddown || bestcandmayroundup) {
                double objgain;
                objgain = (1.0 - frac) * obj;

                if (CC_OURABS(1.0 - frac) < 0.01) {
                    objgain *= 1000;
                }

                if (objgain < bestobjgain || (objgain == bestobjgain && frac > bestfrac)) {
                    *bestcand = var;
                    bestobjgain = objgain;
                    bestfrac = frac;
                    //bestcandmayrounddown = mayrounddown;
                    //bestcandmayroundup = mayroundup;
                }
            }
        } else {
            if (CC_OURABS(1.0 - frac) < 0.01) {
                frac += 10.0;
            }

            if (bestcandmayrounddown || bestcandmayroundup || frac > bestfrac) {
                *bestcand = var;
                bestfrac = frac;
                bestcandmayrounddown = 0;
                bestcandmayroundup = 0;
            }

            assert(bestfrac < DBL_MAX);
        }
    }

    *bestcandmayround = bestcandmayrounddown || bestcandmayroundup;
CLEAN:
    WCTbranchcand_free(branchcand);
    return val;
}


int constructsolution(wctdata *pd, int nmachines, int *success) {
    int i, j, k, val = 0;
    int *perm = (int *) NULL;
    int *covered = (int *) NULL;
    double *x = (double *) NULL;
    double sum = .0;
    int count = 0;
    GList *list = (GList *) NULL;
    perm = CC_SAFE_MALLOC(pd->ccount, int);
    CCcheck_NULL_2(perm, "Failed to allocate memory");

    for (i = 0; i < pd->ccount; ++i) {
        perm[i] = i;
    }

    covered = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(covered, "Failed to allocate memory");
    fill_int(covered, pd->njobs, 0);
    x = CC_SAFE_MALLOC(pd->ccount, double);
    CCcheck_NULL_2(x, "Failed to allocate memory");
    val = wctlp_x(pd->LP, x , 0);
    CCcheck_val_2(val, "Failed to allocate x");
    Scheduleset_permquicksort(perm, pd->cclasses, pd->ccount, Scheduleset_more);

    for (i = 0; i < pd->ccount && count < pd->njobs && sum < nmachines; ++i) {
        int temp = perm[i];
        double upfrac;
        double downfrac;

        if (x[temp] > EPSILON) {
            if (x[temp] > 1.0) {
                x[temp] = 1.0;
            }

            upfrac = ceil(x[temp]) - x[temp];
            downfrac = x[temp] - floor(x[temp]);

            if (upfrac <= downfrac) {
                x[temp] = ceil(x[temp]);
            } else {
                x[temp] = floor(x[temp]);
            }

            if (x[temp] == 1.0) {
                sum += x[temp];
                list = g_list_append(list, GINT_TO_POINTER(temp));

                for (j = 0; j < pd->cclasses[temp].count; ++j) {
                    if (!covered[pd->cclasses[temp].members[j]]) {
                        covered[pd->cclasses[temp].members[j]] = 1;
                        count++;
                    }
                }
            }
        } else {
            x[temp] = .0;
        }
    }

    if (count == pd->njobs) {
        if (sum <= (double) nmachines) {
            GList *it;
            int counter = 0;
            fill_int(covered, pd->njobs, 0);
            Schedulesets_free(&(pd->bestcolors), &(pd->nbbest));
            pd->bestcolors = (Scheduleset *) realloc(pd->bestcolors,
                             nmachines * sizeof(Scheduleset));
            CCcheck_NULL_2(pd->bestcolors, "Failed to realloc pd->bestcolors");

            for (i = 0; i < nmachines; i++) {
                Scheduleset_init(pd->bestcolors + i);
            }

            pd->nbbest = 0;
            int tmp = 0;

            for (it = list; it && counter < pd->njobs; it = it->next) {
                i = GPOINTER_TO_INT(it->data);
                tmp = pd->nbbest;
                pd->bestcolors[tmp].members = CC_SAFE_MALLOC(pd->cclasses[i].count, int);

                for (k = 0; k < pd->cclasses[i].count; ++k) {
                    if (!covered[pd->cclasses[i].members[k]]) {
                        covered[pd->cclasses[i].members[k]] = 1;
                        pd->bestcolors[tmp].members[pd->bestcolors[tmp].count++] =
                            pd->cclasses[i].members[k];
                        pd->bestcolors[tmp].totweight += pd->duration[pd->cclasses[i].members[k]];
                        pd->bestcolors[tmp].totwct += pd->weights[pd->cclasses[i].members[k]] *
                                                      pd->bestcolors[tmp].totweight;
                        counter++;
                    }
                }

                pd->nbbest++;
            }

            pd->nbbest = nmachines;
            /** Check if solution is feasible */
        } else {
            *success = 0;
        }
    } else {
        *success = 0;
    }

CLEAN:
    g_list_free(list);
    CC_IFFREE(x, double);
    CC_IFFREE(perm, int);
    CC_IFFREE(covered, int);
    return val;
}

int heur_compute_lower_bound(wctproblem *problem, wctdata *pd) {
    int j, val = 0;
    int iterations = 0;
    int break_while_loop = 1;
    double cur_cputime;
    double lb_cputime;
    double start_time;
    double last_lower_bound = DBL_MAX;
    int    nnonimprovements     = 0;
    int status = 0;
    CCcheck_val_2(val, "Failed at adjGraph_build");

    if (dbg_lvl() > 1) {
        printf("Starting compute_lower_bound with lb %d and ub %d at depth %d(id = %d, opt_track = %d)\n",
               pd->lower_bound, pd->upper_bound, pd->depth, pd->id, pd->opt_track);
    }

    start_time = CCutil_zeit();
    assert(pd->ccount);
    assert(pd->gallocated >= pd->ccount);
    assert(pd->LP != NULL);
    CCcheck_val_2(val, "Failed in kpc_init_model");
    CC_IFFREE(pd->kpc_pi, int);
    pd->kpc_pi = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(pd->kpc_pi, "Failed to allocate memory to pd->kpc_pi");
    pd->retirementage = (int)sqrt((double)pd->njobs) + 30;
    /** Initialize pricing problem */
    CCutil_start_resume_time(&(problem->tot_lb_heur));

    do {
        iterations++;
        cur_cputime = CCutil_zeit();
        val = wctlp_optimize(pd->LP, &status);
        CCcheck_val(val, "wctlp_optimize failed");
        cur_cputime = CCutil_zeit() - cur_cputime;

        if (dbg_lvl() > 1) {
            printf("Simplex took %f seconds.\n", CCutil_zeit() - cur_cputime);
            fflush(stdout);
        }

        val = wctlp_pi(pd->LP, pd->pi);
        CCcheck_val_2(val, "wctlp_pi failed");
        make_pi_feasible(pd);
        double2int(pd->kpc_pi, &(pd->kpc_pi_scalef_heur), pd->pi, pd->njobs);
        val = compute_objective_heur(pd);
        CCcheck_val_2(val, "compute_objective failed");
        pd->dbl_est_lower_bound_heur = pd->LP_lower_bound_heur;

        if (iterations < pd->maxiterations) {
            if (fabs(last_lower_bound - pd->dbl_est_lower_bound) < 0.0001) {
                nnonimprovements++;
            } else {
                nnonimprovements = 0;
            }

            last_lower_bound = pd->LP_lower_bound_heur;
            CCutil_start_resume_time(&problem->tot_pricing);
            /** Solving Pricing Problem */
            CCutil_suspend_timer(&problem->tot_pricing);

            for (j = 0; j < pd->nnewsets; j++) {
                val = wctlp_addcol(pd->LP, pd->newsets[j].count,
                                   pd->newsets[j].members,
                                   pd->coef, 1.0, 0.0, 1.0,
                                   wctlp_CONT, NULL);
                CCcheck_val_2(val, "wctlp_addcol failed");
            }

            break_while_loop = ((pd->nnewsets == 0));
            add_newsets(pd);
        }

        lb_cputime = CCutil_zeit() - start_time;
    } while ((iterations < pd->maxiterations) && !break_while_loop
             && lb_cputime <= problem->parms.branching_cpu_limit);

    if (iterations < pd->maxiterations &&
            lb_cputime <= problem->parms.branching_cpu_limit) {
        val = wctlp_optimize(pd->LP, &status);
        CCcheck_val_2(val, "wctlp_optimize failed");
        val = compute_objective_heur(pd);
        CCcheck_val_2(val, "compute_objective failed");

        if (dbg_lvl() > 0) {
            printf("Found lb = %d (%f) upper_bound = %d (id= %d, iterations = %d,opt_track = %d).\n",
                   pd->lower_bound, pd->dbl_safe_lower_bound, pd->upper_bound,
                   pd->id, iterations, pd->opt_track);
        }

        if (pd->x != NULL) {
            CC_IFFREE(pd->x, double);
        }

        pd->x = CC_SAFE_MALLOC(pd->ccount, double);
        CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
        val = wctlp_x(pd->LP, pd->x, 0);
        CCcheck_val_2(val, "Failed in wctlp_x");
        pd->status = LP_bound_computed;
    } else {
        pd->status = LP_bound_estimated;
    }

    fflush(stdout);
    CCutil_suspend_timer(&(problem->tot_lb_heur));
CLEAN:
    return val;
}


void heur_init(wctdata *pd) {
    heur_diving *heur_data = &(pd->heur_data);
    heur_data->backtrack = 1;
    heur_data->usefarkasonly = 0;
    heur_data->maxdiscrepancy = 2;
    heur_data->maxdiscdepth = 3;
    heur_data->roundedsum = DBL_MAX;
    heur_data->ccount = pd->ccount;
    heur_data->roundedsol = (double *) NULL;
    heur_data->sol = (double *) NULL;
}

void heur_free(wctdata *pd) {
    heur_diving *heur_data = &(pd->heur_data);
    CC_IFFREE(heur_data->roundedsol, double);
    CC_IFFREE(heur_data->sol, double);
}

int adjustLP_ceil(wctdata *pd, int bestcand, double bestcandsol) {
    int i, val = 0;
    double *value = (double *) NULL;
    double *rhs = (double *) NULL;
    int *vind = (int *) NULL;
    pd->partial_sol += (double) pd->cclasses[bestcand].totwct;
    assert(bestcandsol <= 1.0 && bestcandsol >= -EPSILON);
    val = wctlp_setbound(pd->LP, bestcand, 'L', ceil(bestcandsol));
    CCcheck_val(val, "Failed to set bound");
    rhs = CC_SAFE_MALLOC(pd->njobs + 1, double);
    CCcheck_NULL_2(rhs, "Failed to allocate memory");
    val = GRBgetdblattrarray(pd->LP->model, GRB_DBL_ATTR_RHS, 0, pd->njobs + 1,
                             rhs);
    CHECK_VAL_GRB(val, "Failed in getting the RHS", pd->LP->env);

    for (i = 0; i < pd->cclasses[bestcand].count + 1 ; i++) {
        pd->rhs[pd->cclasses[bestcand].members[i]] -= 1.0;
        printf("test heur %d\n", pd->cclasses[bestcand].members[i]);
    }

    value = CC_SAFE_MALLOC(pd->cclasses[bestcand].count + 1, double);
    CCcheck_NULL_2(value, "Failed to allocate memory");
    fill_dbl(value, pd->cclasses[bestcand].count + 1, 0);
    vind = CC_SAFE_MALLOC(pd->cclasses[bestcand].count + 1, int);
    CCcheck_NULL_2(vind, "Failed to allocate memory");

    for (i = 0; i < pd->ccount; ++i) {
        //if (i != bestcandsol) {
        fill_int(vind, pd->cclasses[bestcand].count + 1, i);
        val = GRBchgcoeffs(pd->LP->model, pd->cclasses[bestcand].count,
                           pd->cclasses[bestcand].members, vind, value);
        CHECK_VAL_GRB(val, "Failed at GRBchgcoeffs", pd->LP->env);
        //}
    }

    val = GRBsetdblattrarray(pd->LP->model, GRB_DBL_ATTR_RHS, 0, pd->njobs ,
                             pd->rhs);
    CHECK_VAL_GRB(val, "Failed at GRB_DBL_ATTR_RHS", pd->LP->env);
    val = GRBupdatemodel(pd->LP->model);
    CHECK_VAL_GRB(val, "Failed to update model", pd->LP->env);
    wctlp_write(pd->LP, "testheur.lp");
    printf("test %d %f\n", bestcand, pd->rhs[pd->njobs]);
    getchar();
CLEAN:
    CC_IFFREE(value, double);
    CC_IFFREE(vind, int);
    CC_IFFREE(rhs, double);
    return val;
}

int adjustLP_floor(wctdata *pd, int var) {
    int val = 0;
    double *value = (double *) NULL;
    double *rhs = (double *) NULL;
    int *vind = (int *) NULL;
    wctlp *LP = pd->LP;
    pd->partial_sol -= (double) pd->cclasses[var].totwct;
    val = wctlp_setbound(LP, var, 'L', 0.0);
    CCcheck_val_2(val, "Failed in ");
    rhs = CC_SAFE_MALLOC(pd->njobs + 1, double);
    CCcheck_NULL_2(rhs, "Failed to allocate memory");
    GRBgetdblattrarray(pd->LP->model, GRB_DBL_ATTR_RHS, 0, pd->njobs + 1, rhs);
    value = CC_SAFE_MALLOC(pd->cclasses[var].count + 1, double);
    CCcheck_NULL_2(value, "Failed to allocate memory");
    fill_dbl(value, pd->cclasses[var].count + 1, 1.0);

    for (int i = 0; i < pd->cclasses[var].count + 1; i++) {
        pd->rhs[pd->cclasses[var].members[i]] += 1.0;
    }

    vind = CC_SAFE_MALLOC(pd->cclasses[var].count, int);
    CCcheck_NULL_2(vind, "Failed to allocate memory");
    fill_int(vind, pd->cclasses[var].count, var);
    val = GRBsetdblattrarray(pd->LP->model, GRB_DBL_ATTR_RHS, 0, pd->njobs,
                             pd->rhs);
    CHECK_VAL_GRB2(val, "Failed to update RHS", LP->env);
    val = GRBchgcoeffs(LP->model, pd->cclasses[var].count,
                       pd->cclasses[var].members, vind, value);
    val = GRBupdatemodel(pd->LP->model);
    CHECK_VAL_GRB(val, "Failed to update model", pd->LP->env);
CLEAN:
    CC_IFFREE(vind, int);
    CC_IFFREE(value, double);
    CC_IFFREE(rhs, double);
    return val;
}


int compute_objective_heur(wctdata *pd) {
    int val = 0;
    val = wctlp_objval(pd->LP, &(pd->LP_lower_bound_heur));
    CCcheck_val_2(val, "wctlp_objval failed");
    pd->lower_bound = (int) ceil(pd->lower_bound);

    if (dbg_lvl() > 0) {
        printf("Current primal LP objective: %19.16f (%lld / %lld) (LP-solver %19.16f).\n",
               pd->dbl_safe_lower_bound,
               (long long) pd->lower_scaled_bound,
               (long long) pd->kpc_pi_scalef, pd->LP_lower_bound_heur);
    }

CLEAN:
    return val;
}

