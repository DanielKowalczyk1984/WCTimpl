#include "util.h"
#include "solverwrapper.h"
#include "PricerSolver.hpp"

extern "C" {
    PricerSolver *newSolver(int *p, int *w, int *r, int *d, int nbjobs, int Hmin, int Hmax){
        return new PricerSolver(p, w, r, d, nbjobs,Hmin,Hmax);
    }

    void solveint(PricerSolver *solver, int *pi, int *sol,int *nsol){
        PricerInfo<int> s;
        s = solver->solveInt(pi);
    }

    void solvedbl(PricerSolver *solver, double *pi, int *sol, int *nsol){
        PricerInfo<double> s;
        s = solver->solveDbl(pi);
    }

    void deletePricerSolver(PricerSolver *solver){
        if (solver)
        {
            delete solver;
            solver = (PricerSolver *) NULL;
        }
    }

}