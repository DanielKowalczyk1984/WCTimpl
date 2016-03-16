#include "wct.h"
#include "wctparms.h"
#include "PricerSolver.hpp"
#include <iostream>
#include <vector>

template<typename T = double, bool reverse = false>
int construct_sol(Scheduleset **set, int *nnewsets, Optimal_Solution<T> &sol, int nbjobs)
{
    int val = 0;
    int nbset = 1;
    std::vector<int> *v = &(sol.jobs);
    Scheduleset *newset = CC_SAFE_MALLOC(1, Scheduleset);
    CCcheck_NULL_2(newset, "Failed to allocate memory newset");
    Scheduleset_init(newset);
    newset->members = CC_SAFE_MALLOC(sol.jobs.size() + 1, int);
    CCcheck_NULL_2(newset->members, "Failed to allocate memory members");

    if (reverse) {
        std::copy(v->rbegin(), v->rend(), newset->members);
    } else {
        std::copy(v->begin(), v->end(), newset->members);
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
    PricerSolver *newSolver(int *p, int *w, int *r, int *d, int nbjobs, int Hmin, int Hmax)
    {
        return new PricerSolver(p, w, r, d, nbjobs, Hmin, Hmax, false, true);
    }

    int solvedblzdd(wctdata *pd)
    {
        int val = 0;
        Optimal_Solution<double> s = pd->solver->solve_duration_zdd_double(pd->pi);

        if (s.obj > 0.00001) {
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), s, pd->njobs);
            CCcheck_val_2(val, "Failed in construct_sol_zdd");
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solvedblbdd(wctdata *pd)
    {
        int val = 0;
        Optimal_Solution<double> s = pd->solver->solve_duration_bdd_double(pd->pi);

        if (s.obj > 0.00001) {
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), s, pd->njobs);
            CCcheck_val_2(val, "Failed to construct_sol_bdd");
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solve_dynamic_programming_ahv(wctdata *pd)
    {
        int val = 0;
        Optimal_Solution<double> s = pd->solver->dynamic_programming_ahv(pd->pi);
        printf("solve_dynamic_programming_ahv\n");

        if (s.obj < -0.000001) {
            val = construct_sol<double, true>(&(pd->newsets), &(pd->nnewsets), s, pd->njobs);
            CCcheck_val_2(val, "Failed in constructing sol");
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solve_farkas_dbl(wctdata *pd)
    {
        int val = 0;
        Optimal_Solution<double> s = pd->solver->solve_farkas_double(pd->pi);
        printf("solve_farkas_dbl\n");

        if (s.obj > 0.000001) {
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), s, pd->njobs);
            CCcheck_val_2(val, "Failed in constructing jobs");
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solve_weight_dbl_bdd(wctdata *pd)
    {
        int val = 0;
        Optimal_Solution<double> s = pd->solver->solve_weight_bdd_double(pd->pi);
        printf("solve_weight_dbl_bdd\n");

        if (s.obj > 0.00001) {
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), s, pd->njobs);
            CCcheck_val_2(val, "Failed in construction")
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }

    int solve_weight_dbl_zdd(wctdata *pd)
    {
        int val = 0;
        Optimal_Solution<double> s = pd->solver->solve_weight_zdd_double(pd->pi);
        printf("solve_weight_dbl_zdd\n" );

        if (s.obj > 0.00001) {
            val = construct_sol(&(pd->newsets), &(pd->nnewsets), s, pd->njobs);
            CCcheck_val_2(val, "Failed in construction")
        } else {
            pd->nnewsets = 0;
        }

CLEAN:
        return val;
    }


    void deletePricerSolver(PricerSolver *solver)
    {
        if (solver) {
            delete solver;
            solver = (PricerSolver *) NULL;
        }
    }

    void compute_subgradient(Optimal_Solution<double> &sol, double *sub_gradient, int nbjobs, int nbmachines)
    {
        fill_dbl(sub_gradient, nbjobs, 1.0);
        sub_gradient[nbjobs] = nbmachines;

        for (auto &v : sol.jobs) {
            sub_gradient[v] -= 1.0;
        }

        sub_gradient[nbjobs] -= 1.0;
    }

    void adjust_alpha(double *pi_out, double *pi_in, double *subgradient, int nbjobs, double &alpha)
    {
        double sum = 0.0;

        for (int i = 0; i <= nbjobs; ++i) {
            sum += subgradient[i] * (pi_out[i] - pi_in[i]);
        }

        if (sum > 0) {
            alpha = alpha + (1 - alpha) * 0.1;
        } else {
            alpha = CC_MAX(0, alpha - 0.01);
        }
    }

    void compute_pi_eta_sep(int vcount, double *pi_sep, double *eta_sep, double alpha, double *pi_in, double *eta_in, double *pi_out, double *eta_out)
    {
        int i;
        double beta = 1.0 - alpha;

        for (i = 0; i < vcount; ++i) {
            pi_sep[i] = alpha * pi_in[i] + beta * pi_out[i];
        }

        *eta_sep = alpha * (*eta_in) + beta * (*eta_out);
    }

    int solve_pricing(wctdata *pd, wctparms *parms)
    {
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

    double compute_lagrange(Optimal_Solution<double> &sol, double *pi, int nbjobs, int nbmachines)
    {
        double result = pi[nbjobs];
        double a = 0.0;
        int i;
        std::vector<int> *v = &(sol.jobs);

        for (std::vector<int>::iterator it = v->begin(); it != v->end(); ++it) {
            result += pi[*it];
        }

        for (i = 0; i < nbjobs; ++i) {
            a += pi[i];
        }

        a += nbmachines * pi[nbjobs];
        result = CC_MIN(0.0, sol.cost - result);
        result = (double)nbmachines * result;
        result += a;
        return result;
    }

    void compute_sol_stab(PricerSolver *solver, wctparms *parms, double *pi, Optimal_Solution<double> *sol)
    {
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

    int construct_sol_stab(wctdata *pd, wctparms *parms, Optimal_Solution<double> &sol)
    {
        int val = 0;

        switch (parms->solver) {
            case bdd_solver:
            case zdd_solver:
                if (sol.obj > 0.000001) {
                    val = construct_sol(&(pd->newsets), &(pd->nnewsets), sol, pd->njobs);
                    CCcheck_val_2(val, "Failed in construction of solution");
                } else {
                    pd->nnewsets = 0;
                }

                break;

            case DP_solver:
                if (sol.obj < -0.000001) {
                    val = construct_sol<double, true>(&(pd->newsets), &(pd->nnewsets), sol, pd->njobs);
                    CCcheck_val_2(val, "Failed in construction of solution");
                } else {
                    pd->nnewsets = 0;
                }

                break;
        }

CLEAN:
        return val;
    }

    int solve_stab(wctdata *pd, wctparms *parms)
    {
        int val = 0;
        int heading_in = 0;
        PricerSolver *solver = pd->solver;
        Optimal_Solution<double> sol;
        heading_in = (pd->eta_in == 0.0) ? 1 : !(CC_OURABS((pd->eta_out - pd->eta_in) / (pd->eta_in)) < 4.0);

        if (heading_in) {
            double result;
            compute_sol_stab(solver, parms, pd->pi, &sol);
            result = compute_lagrange(sol, pd->pi, pd->njobs, pd->nmachines);

            if (result > pd->eta_in) {
                pd->eta_in = result;
                memcpy(pd->pi_in, pd->pi, sizeof(double)*pd->njobs + 1);
            }

            val = construct_sol_stab(pd, parms, sol);
            CCcheck_val_2(val, "Failed in construct solution");
        } else {
            double k = 0.0;
            double alpha;

            do {
                k += 1.0;
                alpha = CC_MAX(0, 1.0 - k * (1.0 - pd->alpha));
                double  result;
                compute_pi_eta_sep(pd->njobs, pd->pi_sep, &(pd->eta_sep), alpha, pd->pi_in, &(pd->eta_in), pd->pi_out, &(pd->eta_out));
                compute_sol_stab(solver, parms, pd->pi_sep, &sol);
                result = compute_lagrange(sol, pd->pi_sep, pd->njobs, pd->nmachines);

                if (result < pd->eta_sep) {
                    val = construct_sol_stab(pd, parms, sol);
                    CCcheck_val_2(val, "Failed in construct_sol_stab");

                    if (result > pd->eta_in) {
                        pd->eta_in = result;
                        memcpy(pd->pi_in, pd->pi_sep, sizeof(double)*pd->njobs + 1);
                    }
                } else {
                    pd->eta_in = result;
                    memcpy(pd->pi_in, pd->pi_sep, sizeof(double)*pd->njobs + 1);
                }
            } while (pd->nnewsets == 0 && alpha > 0.0); /** mispricing check */
        }

        //if(dbg_lvl() > 1) {
        printf("heading = %d, alpha = %f, result of primal bound and Lagragian bound: out =%f, in = %f\n", heading_in, pd->alpha, pd->eta_out, pd->eta_in);
        //}
CLEAN:
        return val;
    }
}


