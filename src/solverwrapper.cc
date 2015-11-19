#include "util.h"
#include "solverwrapper.h"
#include "PricerSolver.hpp"
#include  "datastructsol.h"
#include <iostream>
#include <vector>

extern "C" {
    PricerSolver *newSolver(int *p, int *w, int *r, int *d, int nbjobs, int Hmin, int Hmax) {
        return new PricerSolver(p, w, r, d, nbjobs, Hmin, Hmax);
    }

    void solveint(PricerSolver *solver, int *pi, int *sol, int *nsol, int *cost,int *newsol) {
        PricerInfo<int> s;
        s = solver->solveInt(pi);

        CC_IFFREE(sol, int);
        *nsol = 0;
        *cost = 0;
        *newsol = 0;

        if (s.jobs.size() > 0 && s.obj > 0) {
            *nsol = (int) s.jobs.size();
            sol = CC_SAFE_MALLOC(*nsol, int);
            std::copy(s.jobs.rbegin(), s.jobs.rend(), sol);
            *cost = s.cost;
            *newsol = 1;
        } 
    }

    void solvedbl(PricerSolver *solver, double *pi, int *sol, int *nsol, int *cost, int *newsol) {
        PricerInfo<double> s;
        s = solver->solveDbl(pi);

        CC_IFFREE(sol, int);
        *nsol = 0;
        *cost = 0;
        *newsol = 0;

        if (s.jobs.size() > 0 && s.obj > 0) {
            *nsol = (int) s.jobs.size();
            sol = CC_SAFE_MALLOC(*nsol, int);
            fill_int(sol, *nsol, 0);
            std::copy(s.jobs.rbegin(), s.jobs.rend(), sol);
            *cost = s.cost;
            *newsol = 1;
        }
    }

    void solvefarkasint(PricerSolver *solver, int *pi, int *sol,int *nsol, int* cost, int *newsol){
        PricerInfo<int> s;
        s = solver->solvefarkasInt(pi);

        CC_IFFREE(sol, int);
        *nsol = 0;
        *cost = 0;
        *newsol = 0;

        if (s.jobs.size() > 0 && s.obj > 0) {
            *nsol = (int) s.jobs.size();
            sol = CC_SAFE_MALLOC(*nsol, int);
            fill_int(sol, *nsol, 0);
            std::copy(s.jobs.rbegin(), s.jobs.rend(), sol);
            *cost = s.cost;
            *newsol = 1;
        }
    }

    void solvefarkasdbl(PricerSolver *solver, double *pi, int *sol,int *nsol, int* cost, int *newsol){
        PricerInfo<double> s;
        s = solver->solvefarkasDbl(pi);

        CC_IFFREE(sol, int);
        *nsol = 0;
        *cost = 0;
        *newsol = 0;

        if (s.jobs.size() > 0 && s.obj > 0) {
            *nsol = (int) s.jobs.size();
            sol = CC_SAFE_MALLOC(*nsol, int);
            fill_int(sol, *nsol, 0);
            std::copy(s.jobs.rbegin(), s.jobs.rend(), sol);
            *cost = s.cost;
            *newsol = 1;
        }
    }


    void deletePricerSolver(PricerSolver *solver) {
        if (solver)
        {
            delete solver;
            solver = (PricerSolver *) NULL;
        }
    }

}