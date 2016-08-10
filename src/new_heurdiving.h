#include "wctprivate.h"
#include "wct.h"

#ifdef __cplusplus
extern "C" {
#endif

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
    int nb_rows;
    int nb_cols;

    /** Lower bound and dual bound value */
    int lower_bound;
    int upper_bound;
    double LP_lower_bound;
    double LP_lower_bound_dual;


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

    /** status of column genearation heuristic */
    int status;

} LP_data_CG_heur;

typedef struct _pair_build_LP {
    int counter;
    int *covered;
    int *val;
    double rounded_sol;
    LP_data_CG_heur *data;
} pair_build_LP;

void pair_build_LP_init(pair_build_LP *user_data, int *covered,
                        LP_data_CG_heur *data, int *val);
/** help functions for CG heuristics */
double *update_rhs(double *a, int njobs, GHashTable *rounded_sol);
void update_rhs_iterator(gpointer key, gpointer value, gpointer user_data);
void update_K_iterator(gpointer key, gpointer value, gpointer user_data);
void update_partsol_iterator(gpointer key, gpointer value, gpointer user_data);
void update_column_iterator(gpointer key, gpointer value, gpointer user_data);
void complete_RM_iterator(gpointer key, gpointer value, gpointer user_data);
int compute_objective_CG_heur(LP_data_CG_heur *data);

void print_columns_iterator(gpointer key, gpointer value, gpointer user_data);
/** data constructers and destructors */
void LP_data_CG_heur_init(LP_data_CG_heur *data , wctdata *pd);
void LP_data_CG_heur_free(LP_data_CG_heur *data);
int build_lp_CG(LP_data_CG_heur *data);
int pure_diving(wctproblem *problem, wctdata *pd, GHashTable *G, double *a,
                int njobs, double K, GHashTable *part_sol, GHashTable *rounded_sol);
int execute_CG_heur(wctproblem *problem, wctdata *pd);
int compute_lower_bound_CG_heur(wctproblem *problem, LP_data_CG_heur *data);
int solve_stab_CG_heur(wctproblem *problem, LP_data_CG_heur *data);
int add_newsets_CG_heur(LP_data_CG_heur *data);
/** Pricer Solvers */
int solve_pricing_CG_heur(LP_data_CG_heur *data, wctparms *parms);
int solve_dynamic_programming_ahv_CG_heur(LP_data_CG_heur *data);
int solve_farkas_dbl_CG_heur(LP_data_CG_heur *data);
int solve_weight_dbl_bdd_CG_heur(LP_data_CG_heur *data);
int solve_weight_dbl_zdd_CG_heur(LP_data_CG_heur *data);
#ifdef __cplusplus
}
#endif
