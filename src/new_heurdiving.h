#include "wctprivate.h"
#include "wct.h"

/** data struct CG heuristics */
typedef struct _pair_update_column {
    GHashTable *new_G;
    GHashTable *rounded_sol;
} pair_update_column;

typedef struct _pair_test_column {
    int add;
    Scheduleset *set;
} pair_test_column;


typedef struct _LP_data_CG_heur {
    /** Wentges smoothing technique */
    double *pi_in;
    double *pi_out;
    double *pi_sep;
    double *subgradient;
    double *subgradient_in;
    double *rhs;
    double eta_in;
    double eta_out;
    double eta_sep;
    double alpha;

    //The column generation lp information
    wctlp *LP;
    double *x;
    double *coef;
    double *pi;
    int nb_row;

    /** Pointer to the data of the node */
    wctdata *pd;

    /** set of all the columns */
    GHashTable *G;
    GHashTable *rounded_sol;
    GHashTable *partial_sol;
    int nb_rounded_sol;
    int nb_partial_sol;
    double *a;
    double K;

} LP_data_CG_heur;

typedef struct _pair_build_LP {
    int counter;
    int *covered;
    int *val;
    double rounded_sol;
    LP_data_CG_heur *data;
} pair_build_LP;

/** help functions for CG heuristics */
double *update_rhs(double *a, int njobs, GHashTable *rounded_sol);
void update_rhs_iterator(gpointer key, gpointer value, gpointer user_data);
void update_K_iterator(gpointer key, gpointer value, gpointer user_data);
void update_partsol_iterator(gpointer key, gpointer value, gpointer user_data);
void update_column_iterator(gpointer key, gpointer value, gpointer user_data);

void print_columns_iterator(gpointer key, gpointer value, gpointer user_data);
/** data constructers and destructors */
void LP_data_CG_heur_init(LP_data_CG_heur *data , wctdata *pd);
void LP_data_CG_heur_free(LP_data_CG_heur *data);
int build_lp_CG(LP_data_CG_heur *data);
int pure_diving(wctdata *pd, GHashTable *G, double *a, int njobs, double K, GHashTable *part_sol, GHashTable *rounded_sol);
int execute_CG_heur(wctproblem *problem, wctdata *pd);
