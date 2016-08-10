#ifndef __HEUR_DIVING_H
#define __HEUR_DIVING_H

#ifdef __cplusplus
extern "C" {
#endif


#include "wctparms.h"
#include "wct.h"
#include "lp.h"

#define EPSILON 0.00001

typedef struct BranchCand {
    int *lpcands;
    double *lpcandssol;
    double *lpcandsfrac;

    int lpcandssize;
    int nlpcands;
} WCTbranchcand;


enum resultdiving {
    DELEYAD = 0,
    DIDNOTRUN = 1,
    DIDNOTFIND = 2,
    FOUNDSOL = 3,
};

int heur_divingselect_var(wctdata *pd, int *tabulist, int tabulistsize,
                          int *bestcand, int *bestcanmayround, WCTbranchcand *branchcand);
int branchcandlp(wctdata *pd, WCTbranchcand *branchcand);
int constructsolution(wctdata *pd, int nmachines, int *success);
int heur_compute_lower_bound(wctproblem *, wctdata *);
int adjustLP_ceil(wctdata *pd, int bestcand, double bestcandsol);
int adjustLP_floor(wctdata *pd, int var);
int compute_objective_heur(wctdata *pd);


int WCTbranchcand_alloc(WCTbranchcand *branchcand, int ncol);
void WCTbranchcand_init(WCTbranchcand *branchcand);
void WCTbranchcand_free(WCTbranchcand *branchcand);



#ifdef __cplusplus
}
#endif
#endif // __HEUR_DIVING_H
