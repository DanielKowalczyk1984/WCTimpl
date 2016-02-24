#include "wct.h"
#include "PricerSolver.hpp"
#include <iostream>
#include <vector>

extern "C" {
    PricerSolver *newSolver(int *p, int *w, int *r, int *d, int nbjobs, int Hmin, int Hmax) {
        return new PricerSolver(p, w, r, d, nbjobs, Hmin, Hmax,false,true);
    }

    int construct_sol_bdd(Scheduleset **set, int *nnewsets, PricerInfoBDD<double> &sol, int nbjobs){
        int val = 0;
        int nbset = 1;
        std::vector<int> *v = &(sol.jobs);

        Scheduleset *newset = CC_SAFE_MALLOC(1, Scheduleset);
        CCcheck_NULL_2(newset, "Failed to allocate memory newset");

        Scheduleset_init(newset);
        newset->members = CC_SAFE_MALLOC(sol.jobs.size() + 1, int);
        CCcheck_NULL_2(newset->members, "Failed to allocate memory members");
        std::copy(v->begin(), v->end(), newset->members);
        newset->totwct = sol.cost;
        newset->count = v->size();
        newset->members[v->size()] = nbjobs;
        *set= newset;
        *nnewsets = 1;

        CLEAN:
        if(val) {
            Schedulesets_free(&(newset), &(nbset));
        }

        return val;
    }

    int construct_sol_zdd(Scheduleset **set, int *nnewsets, Optimal_ZDD<double> &sol, int nbjobs){
        int val = 0;
        int nbset = 1;

        Scheduleset *newset = CC_SAFE_MALLOC(1, Scheduleset);
        CCcheck_NULL_2(newset, "Failed to allocate memory newset");

        Scheduleset_init(newset);
        newset->members = CC_SAFE_MALLOC(sol.jobs.size() + 1, int);
        CCcheck_NULL_2(newset->members, "Failed to allocate memory members");
        std::copy(sol.jobs.begin(), sol.jobs.end(), newset->members);
        newset->totwct = sol.cost;
        newset->count = sol.jobs.size();
        newset->members[sol.jobs.size()] = nbjobs;
        *set= newset;
        *nnewsets = 1;

        CLEAN:
        if(val) {
            Schedulesets_free(&(newset), &(nbset));
        }
        return val;
    }

    void solveint(PricerSolver *solver, int *pi, int **sol, int *nsol, int *cost,int *newsol) {

    }

    void calc(PricerSolver *solver, double *pi){
        solver->calculate_reducedcost_iteration(pi);
    }

    int solvedblzdd(wctdata *pd) {
        int val = 0;
        Optimal_ZDD<double> s = pd->solver->solve_duration_zdd_double(pd->pi);

        if(s.obj > 0.00001) {
            val = construct_sol_zdd(&(pd->newsets), &(pd->nnewsets), s, pd->njobs);
            CCcheck_val_2(val, "Failed in construct_sol_zdd");
        } else {
            pd->nnewsets = 0;
        }

        CLEAN:
        return val;
    }

    int solvedblbdd(wctdata *pd){
        int val = 0;
        PricerInfoBDD<double> s = pd->solver->solve_duration_bdd_double(pd->pi);

        if(s.obj > 0.00001) {
            val = construct_sol_bdd(&(pd->newsets), &(pd->nnewsets), s, pd->njobs);
            CCcheck_val_2(val, "Failed to construct_sol_bdd");
        } else {
            pd->nnewsets = 0;
        }

        CLEAN:
        return val;
        
    }

    void solvefarkasint(PricerSolver *solver, int *pi, int **sol,int *nsol, int* cost, int *newsol){
        
    }

    void solvefarkasdbl(PricerSolver *solver, double *pi, int **sol,int *nsol, int* cost, int *newsol){
        
    }


    void deletePricerSolver(PricerSolver *solver) {
        if (solver)
        {
            delete solver;
            solver = (PricerSolver *) NULL;
        }
    }

    void compute_pi_eta_sep(int vcount, double *pi_sep, double *eta_sep, double alpha, double *pi_in, double *eta_in, double *pi_out, double *eta_out)
    {
        int i;
        double beta = 1.0 - alpha;

        for(i = 0; i < vcount; ++i) {
            pi_sep[i] = alpha * pi_in[i] + beta * pi_out[i];
        }

        *eta_sep = alpha * (*eta_in) + beta * (*eta_out);
    }

    double compute_lagrange_bdd(PricerInfoBDD<double> &sol, double *pi, int nbjobs, int nbmachines){
        double result = 0.0;
        double a = 0.0;
        int i;
        std::vector<int> *v = &(sol.jobs);

        for (std::vector<int>::iterator it = v->begin(); it != v->end(); ++it){
            result += pi[*it];
        }

        for (i = 0; i < nbjobs; ++i){
            a += pi[i];
        }

        result = CC_MIN(0.0, sol.cost - result);
        result = (double)nbmachines * result;
        result += a;

        return result;

    }

    int solve_stab_bdd(wctdata *pd)
    {
        int val = 0;
        int heading_in = 0;
        PricerSolver *solver = pd->solver;
        PricerInfoBDD<double> sol;
        //printf("test %f %f\n",pd->eta_in, (pd->eta_in == 0) ? 0.0 : CC_OURABS((pd->eta_out - pd->eta_in)/(pd->eta_in)) );
        if(  (pd->eta_in == 0.0) ? 1 : CC_OURABS((pd->eta_out - pd->eta_in)/(pd->eta_in)) < 0.0005 || CC_OURABS((pd->eta_out - pd->eta_in)/(pd->eta_in)) > 0.1) {
            heading_in = 1;
        }

        if(heading_in) {
            double result;
            sol = solver->solve_duration_bdd_double(pd->pi);
            result = compute_lagrange_bdd(sol, pd->pi, pd->njobs, pd->nmachines);
            printf("test result = %f\n", result);
            printf("test eta_in = %f\n", pd->eta_in);

            if(result > pd->eta_in) {
                pd->eta_in = result;
                memcpy(pd->pi_in, pd->pi, sizeof(double)*pd->njobs);
            }

            if(sol.obj > pd->kpc_pi_scalef) {
                val = construct_sol_bdd(&(pd->newsets), &(pd->nnewsets), sol, solver->nbjobs);
                CCcheck_val_2(val, "Failed in construction of solution");
            } else {
                pd->nnewsets = 0;
            }

        } else {
            double k = 0.0;
            double alpha;

            do {
                printf("Computining stab %f\n", k);
                k += 1.0;
                alpha = CC_MAX(0, 1.0 - k * (1 - pd->alpha));
                double  result;

                compute_pi_eta_sep(pd->njobs, pd->pi_sep, &(pd->eta_sep), alpha, pd->pi_in, &(pd->eta_in), pd->pi_out, &(pd->eta_out));
                sol = solver->solve_duration_bdd_double(pd->pi_sep);
                result = compute_lagrange_bdd(sol, pd->pi_sep, pd->njobs, pd->nmachines);

                if(result < pd->eta_sep) {
                    printf("tets test\n");
                    val = construct_sol_bdd(&(pd->newsets), &(pd->nnewsets), sol, solver->nbjobs);
                    CCcheck_val_2(val, "Failed in construct_sol_combo");

                    if(result > pd->eta_in) {
                        pd->eta_in = result;
                        memcpy(pd->pi_in, pd->pi_sep, sizeof(double)*pd->njobs);
                    }
                } else {
                    pd->eta_in = result;
                    memcpy(pd->pi_in, pd->pi_sep, sizeof(double)*pd->njobs);
                    result = compute_lagrange_bdd(sol, pd->pi_out, pd->njobs, pd->nmachines);

                    if(result < pd->eta_out) {
                        val = construct_sol_bdd(&(pd->newsets), &(pd->nnewsets), sol, solver->nbjobs);
                        CCcheck_val_2(val, "Failed in construct_sol_combo");
                    }

                }
            } while(pd->nnewsets == 0 && alpha > 0.0);
        }

        //if(dbg_lvl() > 1) {
            printf("result of primal bound and Lagragian bound: out =%f, in = %f\n", pd->eta_out, pd->eta_in);
        //}

    CLEAN:
        return val;
    }

}