#ifndef SOLVERWRAPER_H
#define SOLVERWRAPER_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PricerSolver PricerSolver;

PricerSolver *newSolver(int *p, int *w, int *r, int *d, int nbjobs, int Hmin, int Hmax);

void solvedblbdd(PricerSolver *solver, double *pi, int **sol, int *nsol,  int *cost, int *newsol);

void solvedblzdd(PricerSolver *solver, double *pi);

void solvedbl(PricerSolver *solver, double *pi, int **sol, int *nsol, int* cost, int *newsol);

void solveint(PricerSolver *solver, int *pi, int **sol,int *nsol, int* cost, int *newsol);

void solvefarkasdbl(PricerSolver *solver, double *pi, int **sol, int *nsol, int* cost, int *newsol);

void solvefarkasint(PricerSolver *solver, int *pi, int **sol,int *nsol, int* cost, int *newsol);

PricerSolver *addDuedateConstraint(PricerSolver *solver);

PricerSolver *addReleasetimeConstraint(PricerSolver *solver);

void addRestriction(PricerSolver *solver, int *jobs);

void deletePricerSolver(PricerSolver *solver);

void calc(PricerSolver *solver, double *pi);

#ifdef __cplusplus
}
#endif

#endif // SOLVERWRAPER_H
