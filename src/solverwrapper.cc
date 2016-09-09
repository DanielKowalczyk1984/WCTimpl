#include "wct.h"
#include "new_heurdiving.h"
#include "wctparms.h"
#include "PricerSolver.hpp"
#include <iostream>
#include <vector>

template<typename T = double, bool reverse = false>
int construct_sol(Scheduleset **set, int *nnewsets, int *d,
                  Optimal_Solution<T> &sol, int nbjobs) {
    int val = 0;
    int nbset = 1;
    int C = 0;
    std::vector<int> *v = &(sol.jobs);
    Scheduleset *newset = CC_SAFE_MALLOC(1, Scheduleset);
    CCcheck_NULL_2(newset, "Failed to allocate memory newset");
    Scheduleset_init(newset);
    newset->members = CC_SAFE_MALLOC(sol.jobs.size() + 1, int);
    CCcheck_NULL_2(newset->members, "Failed to allocate memory members");
    newset->C = CC_SAFE_MALLOC(sol.jobs.size(), int);
    CCcheck_NULL_2(newset->C, "Failed to allocate memory");
    newset->table = g_hash_table_new(g_direct_hash, g_direct_equal);
    CCcheck_NULL_2(newset->table, "Failed to allocate memory");

    if (reverse) {
        std::copy(v->rbegin(), v->rend(), newset->members);
    } else {
        std::copy(v->begin(), v->end(), newset->members);
    }

    for (size_t i = 0; i < sol.jobs.size(); i++) {
        C += d[newset->members[i]];
        newset->C[i] = C;
        g_hash_table_insert(newset->table, GINT_TO_POINTER(newset->members[i]),
                            newset->C + i);
    }

    newset->totwct = sol.cost;
    newset->totweight = sol.C_max;
    newset->count = sol.jobs.size();
    newset->members[sol.jobs.size()] = nbjobs;
    *set = newset;
    *nnewsets = 1;
CLEAN:

    if (val) {
        Schedulesets_free(&(newset), &(nbset));
    }

    return val;
}

extern "C" {
    PricerSolver *newSolver(int *p, int *w, int *r, int *d, int nbjobs, int Hmin,
                            int Hmax) {
        return new PricerSolver(p, w, r, d, nbjobs, Hmin, Hmax);
    }

    PricerSolver *newSolverDP(int *p, int *w, int *r, int *d, int nbjobs, int Hmin,
                              int Hmax) {
        return new PricerSolver(p, w, r, d, nbjobs, Hmin, Hmax, false);
    }

    PricerSolver *copySolver(PricerSolver *src) {
        return new PricerSolver(*src);
    }

    void print_dot_file(PricerSolver *solver, char* name){
        solver->create_dot_zdd(name);
    }

    void freeSolver(PricerSolver *src) {
        delete src;
        src = (PricerSolver *) NULL;
    }

    int solvedblzdd(wctdata *pd) {
        int val = 0;
        Optimal_Solution<double> s = pd->solver->solve_duration_zdd_double(pd->pi);

        if (s.obj > 0.00001) {
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), pd->duration, s,
                                pd->njobs);
            CCcheck_val_2(val, "Failed in construct_sol_zdd");
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solvedblbdd(wctdata *pd) {
        int val = 0;
        Optimal_Solution<double> s = pd->solver->solve_duration_bdd_double(pd->pi);

        if (s.obj > 0.00001) {
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), pd->duration, s,
                                pd->njobs);
            CCcheck_val_2(val, "Failed to construct_sol_bdd");
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solve_dynamic_programming_ahv(wctdata *pd) {
        int val = 0;
        Optimal_Solution<double> s = pd->solver->dynamic_programming_ahv(pd->pi);

        if (s.obj < -0.00001) {
            val = construct_sol<double, true>(&(pd->newsets), &(pd->nnewsets), pd->duration,
                                              s, pd->njobs);
            CCcheck_val_2(val, "Failed in constructing sol");
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solve_dynamic_programming_ahv_CG_heur(LP_data_CG_heur *data) {
        int val = 0;
        wctdata *pd = data->pd;
        Optimal_Solution<double> s = pd->solver->dynamic_programming_ahv(data->pi);

        if (s.obj < -0.00001) {
            val = construct_sol<double, true>(&(pd->newsets), &(pd->nnewsets), pd->duration,
                                              s, pd->njobs);
            CCcheck_val_2(val, "Failed in constructing sol");
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solve_farkas_dbl(wctdata *pd) {
        int val = 0;
        Optimal_Solution<double> s = pd->solver->solve_farkas_double(pd->pi);

        if (s.obj < -0.00001) {
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), pd->duration, s,
                                pd->njobs);
            CCcheck_val_2(val, "Failed in constructing jobs");
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solve_farkas_dbl_DP(wctdata *pd) {
        int val = 0;
        Optimal_Solution<double> s = pd->solver->dynamic_programming_ahv_farkas(pd->pi);

        if (s.obj < -0.00001) {
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), pd->duration, s,
                                pd->njobs);
            CCcheck_val_2(val, "Failed in constructing jobs");
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solve_farkas_dbl_CG_heur(LP_data_CG_heur *data) {
        int val = 0;
        wctdata *pd = data->pd;
        Optimal_Solution<double> s = pd->solver->solve_farkas_double(data->pi);

        if (s.obj < -0.00001) {
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), pd->duration, s,
                                pd->njobs);
            CCcheck_val_2(val, "Failed in constructing jobs");
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solve_weight_dbl_bdd(wctdata *pd) {
        int val = 0;
        Optimal_Solution<double> s = pd->solver->solve_weight_bdd_double(pd->pi);

        if (s.obj > 0.00001) {
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), pd->duration, s,
                                pd->njobs);
            CCcheck_val_2(val, "Failed in construction")
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solve_weight_dbl_bdd_CG_heur(LP_data_CG_heur *data) {
        int val = 0;
        wctdata *pd = data->pd;
        Optimal_Solution<double> s = pd->solver->solve_weight_bdd_double(data->pi);

        if (s.obj > 0.00001) {
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), pd->duration, s,
                                pd->njobs);
            CCcheck_val_2(val, "Failed in construction")
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solve_weight_dbl_zdd(wctdata *pd) {
        int val = 0;
        Optimal_Solution<double> s = pd->solver->solve_weight_zdd_double(pd->pi);

        if (s.obj > 0.00001) {
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), pd->duration, s,
                                pd->njobs);
            CCcheck_val_2(val, "Failed in construction")
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solve_weight_dbl_zdd_CG_heur(LP_data_CG_heur *data) {
        int val = 0;
        wctdata *pd = data->pd;
        Optimal_Solution<double> s = pd->solver->solve_weight_zdd_double(data->pi);

        if (s.obj > 0.00001) {
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), pd->duration, s,
                                pd->njobs);
            CCcheck_val_2(val, "Failed in construction")
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    void deletePricerSolver(PricerSolver *solver) {
        if (solver) {
            delete solver;
            solver = (PricerSolver *) NULL;
        }
    }

    int calculate_table(PricerSolver *solver, wctparms *parms) {
        int val = 0;

        switch (parms->solver) {
        case bdd_solver:
            solver->init_bdd_table();
            break;

        case zdd_solver:
            solver->init_zdd_table();
            break;

        case DP_solver:
            break;
        }

        switch (parms->construct) {
        case yes_construct:
            break;

        case no_construct:
            solver->init_table_farkas();
            break;
        }

        return val = 0;
    }

    int add_conflict_constraints(PricerSolver *solver, wctparms *parms,
                                 int *elist_same, int ecount_same, int *elist_differ, int  ecount_differ) {
        int val = 0;

        switch (parms->solver) {
        case bdd_solver:
            solver->init_bdd_conflict_solver(elist_same, ecount_same, elist_differ,
                                             ecount_differ);
            break;

        case zdd_solver:
            solver->init_zdd_conflict_solver(elist_same, ecount_same, elist_differ,
                                             ecount_differ);
            break;

        case DP_solver:
            break;
        }

        return val;
    }

    void iterate_dd(PricerSolver *solver) {
        solver->iterate_dd();
    }

    void iterate_zdd(PricerSolver *solver) {
        solver->iterate_zdd();
    }

    int free_conflict_constraints(PricerSolver *solver, wctparms *parms,
                                  int ecount_same, int ecount_differ) {
        int val = 0;

        switch (parms->solver) {
        case bdd_solver:
            solver->free_bdd_solver(ecount_same, ecount_differ);
            break;

        case zdd_solver:
            solver->free_zdd_solver(ecount_same, ecount_differ);
            break;

        case DP_solver:
            break;
        }

        return val;
    }

    size_t get_datasize(PricerSolver *solver) {
        return solver->zdd->size();
    }

    size_t get_numberrows_zdd(PricerSolver *solver) {
        return solver->zdd->root().row();
    }

    size_t get_numberrows_bdd(PricerSolver *solver) {
        return solver->dd->root().row();
    }

    int add_one_conflict(PricerSolver *solver, wctparms *parms, int v1, int v2,
                         int same) {
        int val = 0;

        switch (parms->solver) {
        case bdd_solver:
            solver->init_bdd_one_conflict(v1, v2, same);
            break;

        case zdd_solver:
            solver->init_zdd_one_conflict(v1, v2, same);
            break;

        case DP_solver:
            break;
        }

        return val;
    }

    int init_tables(PricerSolver *solver) {
        int val = 0;
        solver->init_tables();
        return val;
    }

    void set_release_due_time(PricerSolver *solver, int *releasetime,
                              int *duetime) {
        solver->set_release_due_time(releasetime, duetime);
    }

    void compute_subgradient(Optimal_Solution<double> &sol, double *sub_gradient,
                             double *rhs, int nbjobs, int nbmachines) {
        fill_dbl(sub_gradient, nbjobs, 1.0);
        sub_gradient[nbjobs] = nbmachines;

        for (auto &v : sol.jobs) {
            sub_gradient[v] -= (double) nbmachines * 1.0;
        }
    }

    void adjust_alpha(double *pi_out, double *pi_in, double *subgradient,
                      int nbjobs, double &alpha) {
        double sum = 0.0;

        for (int i = 0; i <= nbjobs; ++i) {
            sum += subgradient[i] * (pi_out[i] - pi_in[i]);
        }

        if (sum > 0) {
            alpha = alpha + (1 - alpha) * 0.05;
        } else {
            alpha = CC_MAX(0, alpha - 0.05);
        }
    }

    void compute_pi_eta_sep(int vcount, double *pi_sep, double *eta_sep,
                            double alpha, double *pi_in, double *eta_in, double *pi_out, double *eta_out) {
        int i;
        double beta = 1.0 - alpha;

        for (i = 0; i <= vcount; ++i) {
            pi_sep[i] = alpha * pi_in[i] + beta * pi_out[i];
        }

        *eta_sep = alpha * (*eta_in) + beta * (*eta_out);
    }

    int solve_pricing(wctdata *pd, wctparms *parms) {
        int val = 0;

        switch (parms->solver) {
        case bdd_solver:
            val = solve_weight_dbl_bdd(pd);
            CCcheck_val_2(val, "Failed solve_weight_dbl_bdd");
            break;

        case zdd_solver:
            val = solve_weight_dbl_zdd(pd);
            CCcheck_val_2(val, "Failed solve_weight_dbl_zdd")
            break;

        case DP_solver:
            val = solve_dynamic_programming_ahv(pd);
            CCcheck_val_2(val, "Failed in solve_dynamic_programming_ahv")
            break;
        }

CLEAN:
        return val;
    }

    int solve_pricing_CG_heur(LP_data_CG_heur *data, wctparms *parms) {
        int val = 0;

        switch (parms->solver) {
        case bdd_solver:
            val = solve_weight_dbl_bdd_CG_heur(data);
            CCcheck_val_2(val, "Failed solve_weight_dbl_bdd");
            break;

        case zdd_solver:
            val = solve_weight_dbl_zdd_CG_heur(data);
            CCcheck_val_2(val, "Failed solve_weight_dbl_zdd")
            break;

        case DP_solver:
            val = solve_dynamic_programming_ahv_CG_heur(data);
            CCcheck_val_2(val, "Failed in solve_dynamic_programming_ahv")
            break;
        }

CLEAN:
        return val;
    }

    double compute_lagrange(Optimal_Solution<double> &sol, double *rhs, double *pi,
                            int nbjobs, int nbmachines) {
        double result = 0;
        double a = 0.0;
        int i;
        std::vector<int> *v = &(sol.jobs);

        for (std::vector<int>::iterator it = v->begin(); it != v->end(); ++it) {
            result += pi[*it];
        }

        for (i = 0; i < nbjobs; ++i) {
            a += rhs[i] * pi[i];
        }

        //a += rhs[nbjobs] * pi[nbjobs];
        result = CC_MIN(0.0, sol.cost - result);
        result = rhs[nbjobs] * result;
        result += a;
        return result;
    }

    void compute_lagrange_iterator(gpointer key, gpointer value,
                                   gpointer user_data) {
        double *a = (double *) user_data;
        Scheduleset *set = (Scheduleset *) key;
        *a += set->totwct;
    }

    double compute_lagrange_CG_heur(double *pi, int njobs, GHashTable *part_sol,
                                    double *a, double K, Optimal_Solution<double> &sol) {
        double result = pi[njobs];
        double temp = 0.0;
        int i;
        std::vector<int> *v = &(sol.jobs);

        for (std::vector<int>::iterator it = v->begin(); it != v->end(); ++it) {
            result += pi[*it];
        }

        for (i = 0; i < njobs; ++i) {
            temp += a[i] * pi[i];
        }

        temp += K * pi[njobs];
        g_hash_table_foreach(part_sol, compute_lagrange_iterator, &a);
        result = CC_MIN(0.0, sol.cost - result);
        result = K * result;
        result += temp;
        return result;
    }

    double compute_reduced_cost(Optimal_Solution<double> &sol, double *pi,
                                int nbjobs) {
        double result = pi[nbjobs];
        std::vector<int> *v = &(sol.jobs);

        for (std::vector<int>::iterator it = v->begin();  it != v->end(); it++) {
            result += pi[*it];
        }

        result = sol.cost - result;
        return result;
    }

    void compute_sol_stab(PricerSolver *solver, wctparms *parms, double *pi,
                          Optimal_Solution<double> *sol) {
        switch (parms->solver) {
        case bdd_solver:
            *sol = solver->solve_weight_bdd_double(pi);
            break;

        case zdd_solver:
            *sol = solver->solve_weight_zdd_double(pi);
            break;

        case DP_solver:
            *sol = solver->dynamic_programming_ahv(pi);
            break;
        }
    }

    int construct_sol_stab(wctdata *pd, wctparms *parms,
                           Optimal_Solution<double> &sol) {
        int val = 0;

        switch (parms->solver) {
        case bdd_solver:
        case zdd_solver:
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), pd->duration, sol,
                                pd->njobs);
            CCcheck_val_2(val, "Failed in construction of solution");
            break;

        case DP_solver:
            val = construct_sol<double, true>(&(pd->newsets), &(pd->nnewsets), pd->duration,
                                              sol, pd->njobs);
            CCcheck_val_2(val, "Failed in construction of solution");
            break;
        }

CLEAN:
        return val;
    }

    int solve_stab(wctdata *pd, wctparms *parms) {
        int val = 0;
        int heading_in = 0;
        PricerSolver *solver = pd->solver;
        Optimal_Solution<double> sol;
        heading_in = (pd->eta_in == 0.0) ? 1 : !(CC_OURABS((pd->eta_out -
                     pd->eta_in) / (pd->eta_in)) < 4.0);

        if (heading_in) {
            double result;
            compute_sol_stab(solver, parms, pd->pi, &sol);
            result = compute_lagrange(sol, pd->rhs, pd->pi, pd->njobs, pd->nmachines);

            if (result > pd->eta_in) {
                pd->eta_in = result;
                memcpy(pd->pi_in, pd->pi, sizeof(double) * (pd->njobs + 1));
            }

            val = construct_sol_stab(pd, parms, sol);
            CCcheck_val_2(val, "Failed in construct solution");
        } else {
            double k = 0.0;
            double alpha;
            bool mispricing = false;
            double  result_sep;
            double result_out;

            do {
                k += 1.0;
                alpha = CC_MAX(0, 1.0 - k * (1.0 - pd->alpha));
                compute_pi_eta_sep(pd->njobs, pd->pi_sep, &(pd->eta_sep), alpha, pd->pi_in,
                                   &(pd->eta_in), pd->pi_out, &(pd->eta_out));
                compute_sol_stab(solver, parms, pd->pi_sep, &sol);
                result_sep = compute_lagrange(sol, pd->rhs, pd->pi_sep, pd->njobs,
                                              pd->nmachines);

                if (result_sep < pd->eta_sep) {
                    val = construct_sol_stab(pd, parms, sol);
                    CCcheck_val_2(val, "Failed in construct_sol_stab");
                    pd->update = 1;
                    mispricing = false;
                } else {
                    result_out = compute_lagrange(sol, pd->rhs, pd->pi_out, pd->njobs,
                                                  pd->nmachines);

                    if (result_out < pd->eta_out) {
                        val = construct_sol_stab(pd, parms, sol);
                        CCcheck_val_2(val, "Failed in construct_sol_stab");
                        mispricing = false;
                        pd->update = 0;
                    } else {
                        mispricing = true;
                    }
                }
            } while (mispricing && alpha > 0); /** mispricing check */

            if (result_sep > pd->eta_in) {
                pd->eta_in = result_sep;
                memcpy(pd->pi_in, pd->pi_sep, sizeof(double) * (pd->njobs + 1));
            }
        }

        if (dbg_lvl() > 0) {
            printf("heading = %d, alpha = %f, result of primal bound and Lagragian bound: out =%f, in = %f\n",
                   heading_in, pd->alpha, pd->eta_out, pd->eta_in);
        }

CLEAN:
        return val;
    }

    int solve_stab_CG_heur(wctproblem *problem, LP_data_CG_heur *data) {
        int val = 0;
        int heading_in = 0;
        wctparms *parms = &(problem->parms);
        wctdata *pd = data->pd;
        PricerSolver *solver = pd->solver;
        Optimal_Solution<double> sol;
        heading_in = (pd->eta_in == 0.0) ? 1 : !(CC_OURABS((pd->eta_out -
                     pd->eta_in) / (pd->eta_in)) < 4.0);

        if (heading_in) {
            double result;
            compute_sol_stab(solver, parms, data->pi, &sol);
            result = compute_lagrange_CG_heur(data->pi, pd->njobs, data->partial_sol,
                                              data->a, data->K, sol);

            if (result > data->eta_in) {
                data->eta_in = result;
                memcpy(data->pi_in, data->pi, sizeof(double) * (pd->njobs + 1));
            }

            val = construct_sol_stab(pd, parms, sol);
            CCcheck_val_2(val, "Failed in construct solution");
        } else {
            double k = 0.0;
            double alpha;
            bool mispricing = false;
            double  result_sep;
            double result_out;

            do {
                k += 1.0;
                alpha = CC_MAX(0, 1.0 - k * (1.0 - data->alpha));
                compute_pi_eta_sep(pd->njobs, data->pi_sep, &(data->eta_sep), alpha,
                                   data->pi_in, &(data->eta_in), data->pi_out, &(data->eta_out));
                compute_sol_stab(solver, parms, data->pi_sep, &sol);
                result_sep = compute_lagrange_CG_heur(data->pi_sep, pd->njobs,
                                                      data->partial_sol, data->a, data->K, sol);

                if (result_sep < data->eta_sep) {
                    val = construct_sol_stab(pd, parms, sol);
                    CCcheck_val_2(val, "Failed in construct_sol_stab");
                    pd->update = 1;
                    mispricing = false;
                } else {
                    result_out = compute_lagrange_CG_heur(data->pi_out, pd->njobs,
                                                          data->partial_sol, data->a, data->K, sol);

                    if (result_out < data->eta_out) {
                        val = construct_sol_stab(pd, parms, sol);
                        CCcheck_val_2(val, "Failed in construct_sol_stab");
                        mispricing = false;
                        pd->update = 0;
                    } else {
                        mispricing = true;
                    }
                }
            } while (mispricing && alpha > 0); /** mispricing check */

            if (result_sep > data->eta_in) {
                data->eta_in = result_sep;
                memcpy(data->pi_in, data->pi_sep, sizeof(double) * (pd->njobs + 1));
            }
        }

        //if (dbg_lvl() > 0) {
        printf("heading = %d, alpha = %f, result of primal bound and Lagragian bound: out =%f, in = %f m = %f\n",
               heading_in, data->alpha, data->eta_out, data->eta_in, data->K);
        //}
CLEAN:
        return val;
    }

    int solve_stab_dynamic(wctdata *pd, wctparms *parms) {
        int val = 0;
        PricerSolver *solver = pd->solver;
        Optimal_Solution<double> sol;
        double k = 0.0;
        double alpha;
        double result_sep;
        double result_out;
        bool mispricing = false;

        do {
            k += 1.0;
            alpha = CC_MAX(0.0, 1.0 - k * (1 - pd->alpha));
            compute_pi_eta_sep(pd->njobs, pd->pi_sep, &(pd->eta_sep), alpha, pd->pi_in,
                               &(pd->eta_in), pd->pi_out, &(pd->eta_out));
            compute_sol_stab(solver, parms, pd->pi_sep, &(sol));
            result_sep = compute_lagrange(sol, pd->rhs, pd->pi_sep, pd->njobs,
                                          pd->nmachines);

            if (result_sep < pd->eta_sep) {
                val = construct_sol_stab(pd, parms, sol);
                CCcheck_val_2(val, "Failed in construct_sol_stab");
                pd->update = 1;
                compute_subgradient(sol, pd->subgradient, pd->rhs, pd->njobs, pd->nmachines);
                adjust_alpha(pd->pi_out, pd->pi_in, pd->subgradient, pd->njobs, alpha);
                pd->alpha = alpha;
                mispricing = false;
            } else {
                result_out = compute_lagrange(sol, pd->rhs, pd->pi_out, pd->njobs,
                                              pd->nmachines);

                if (result_out < pd->eta_out) {
                    val = construct_sol_stab(pd, parms, sol);
                    CCcheck_val_2(val, "Failed in construct_sol_stab");
                    mispricing = false;
                    pd->update = 0;
                } else {
                    mispricing = true;
                }
            }
        } while (mispricing && alpha > 0.0);

        if (result_sep > pd->eta_in) {
            pd->eta_in = result_sep;
            memcpy(pd->pi_in, pd->pi_sep, sizeof(double) * (pd->njobs + 1));
        }

        if (dbg_lvl() > 0) {
            printf(" alpha = %f, result of primal bound and Lagragian bound: out =%f, in = %f\n",
                   pd->alpha, pd->eta_out, pd->eta_in);
        }

CLEAN:
        return val;
    }
}


