#include "util.h"
#include "solverwrapper.h"
#include "PricerSolver.hpp"
#include <iostream>
#include <vector>

extern "C" {
    PricerSolver *newSolver(int *p, int *w, int *r, int *d, int nbjobs, int Hmin, int Hmax) {
        return new PricerSolver(p, w, r, d, nbjobs, Hmin, Hmax,false,true);
    }

    void solveint(PricerSolver *solver, int *pi, int **sol, int *nsol, int *cost,int *newsol) {
        PricerInfoBDD<int> s;
        //s = solver->solveInt(pi);
        CC_IFFREE(*sol, int);
        *nsol = 0;
        *cost = 0;
        *newsol = 0;

        if (s.jobs.size() > 0 && s.obj > 0) {
            int *temp = (int *) NULL;
            *nsol = (int) s.jobs.size();
            temp = CC_SAFE_MALLOC(*nsol + 1, int);
            fill_int(temp, *nsol + 1, 0);
            std::copy(s.jobs.rbegin(), s.jobs.rend(), temp);
            temp[*nsol] = solver->nbjobs;
            *cost = s.cost;
            *newsol = 1;
            *sol = temp;
        } 
    }

    void calc(PricerSolver *solver, double *pi){
        solver->calculate_reducedcost_iteration(pi);
    }

    void solvedblzdd(PricerSolver *solver, double *pi) {
        double s;
        s = solver->solve_duration_zdd_double(pi);
        

        printf("test %f\n", s);

        
    }

    void solvedblbdd(PricerSolver *solver, double *pi, int **sol, int *nsol,  int *cost, int *newsol){
        PricerInfoBDD<double> s;
        s = solver->solve_duration_bdd_double(pi);
        CC_IFFREE(*sol, int);
        *nsol = 0;
        *cost = 0;
        *newsol = 0;

        printf("test obj %f\n", s.obj );
        if (s.jobs.size() > 0 && s.obj > 0.0001) {
            int *temp = (int *) NULL;
            *nsol = (int) s.jobs.size();
            temp = CC_SAFE_MALLOC(*nsol + 1, int);
            fill_int(temp, *nsol + 1, 0);
            std::copy(s.jobs.rbegin(), s.jobs.rend(), temp);
            temp[*nsol] = solver->nbjobs;
            *cost = s.cost;
            *newsol = 1;
            *sol = temp;
        }
    }



    void solvefarkasint(PricerSolver *solver, int *pi, int **sol,int *nsol, int* cost, int *newsol){
        PricerInfoBDD<int> s;
        //s = solver->solvefarkasInt(pi);

        CC_IFFREE(*sol, int);
        *nsol = 0;
        *cost = 0;
        *newsol = 0;

        if (s.jobs.size() > 0 && s.obj > 0) {
            int *temp = (int *) NULL;
            *nsol = (int) s.jobs.size();
            temp = CC_SAFE_MALLOC(*nsol + 1, int);
            fill_int(temp, *nsol + 1, 0);
            std::copy(s.jobs.rbegin(), s.jobs.rend(), temp);
            temp[*nsol] = solver->nbjobs;
            *cost = s.cost;
            *newsol = 1;
            *sol = temp;
        }
    }

    void solvefarkasdbl(PricerSolver *solver, double *pi, int **sol,int *nsol, int* cost, int *newsol){
        PricerInfoBDD<double> s;
        //s = solver->solvefarkasDbl(pi);

        CC_IFFREE(*sol, int);
        *nsol = 0;
        *cost = 0;
        *newsol = 0;

        if (s.jobs.size() > 0 && s.obj > 0) {
            int *temp = (int *) NULL;
            *nsol = (int) s.jobs.size();
            temp = CC_SAFE_MALLOC(*nsol + 1, int);
            fill_int(temp, *nsol + 1, 0);
            std::copy(s.jobs.rbegin(), s.jobs.rend(), temp);
            temp[*nsol] = solver->nbjobs;
            *cost = s.cost;
            *newsol = 1;
            *sol = temp;
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