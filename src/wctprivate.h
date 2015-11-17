#ifndef WCT_PRIVATE_H
#define WCT_PRIVATE_H 


#include "util.h"
#include "datastructsol.h"
#include "wctparms.h"
#include "lp.h"
#include "solverwrapper.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * scatter search data types
 */

typedef struct _REFSET{ 
    int newsol;
    GPtrArray *list1;
    int nb1;
    GPtrArray *list2;
    int nb2;
} REFSET;

typedef struct _P{
    int PSize;
    GList *list;
    int lenght;
} P;

typedef struct _SS{
    REFSET *rs;
    P *p;
    int b1;
    int b2;
    
    Job **jobarray;
    int njobs;
    int nmachines;
    int lowerbound; 
    int upperbound;

    double timelimit;
    GRand *random;
    int iter;
    int combine_method;
    
    enum {
        init   = 0,
        add    = 1,
        update = 3,
        opt    = 4
    } status;
} SS;


typedef struct heur_diving {
    int backtrack;
    int usefarkasonly;
    
    int maxdiscdepth;
    int maxdiscrepancy;

    double *sol;
    double *roundedsol;
    int ccount;

    double objval;
    double roundedsum;
} heur_diving;

/**
 * wct data types nodes of branch and bound tree
 */

typedef struct wctdata wctdata;
struct wctdata {
    // The id and depth of the node in the B&B tree
    int id;
    int depth;
    
    enum{
        initialized             = 0,
        LP_bound_estimated      = 1,
        LP_bound_computed       = 2,
        submitted_for_branching = 3,
        finished                = 4,
    }status;
    
    // The job information
    int njobs;
    int nmachines;
    int *duetime;
    int *releasetime;
    int *duration;
    int *weights;
    int *orig_node_ids;
    //data for meta heuristic
    Job* jobarray;
    int H_max;
    int H_min;
    
    //The column generation lp information
    wctlp *LP;
    double *x;
    double *coef;
    double *pi;
    PricerSolver *solver;
    //Colorset(Assignments)
    int ccount;
    Scheduleset*cclasses;
    int dzcount;
    int gallocated;
    Scheduleset*newsets;
    int nnewsets;
    
    int kpc_pi_scalef;
    int kpc_pi_scalef_heur;
    int *kpc_pi;
    int *kpc_pi_heur;
    int lower_bound;
    int upper_bound;
    int lower_scaled_bound;
    double dbl_safe_lower_bound;
    double dbl_est_lower_bound;
    double dbl_est_lower_bound_heur;
    double LP_lower_bound;
    double LP_lower_bound_heur;
    
    
    //KPC instances
    
    
    // Best Solution
    Scheduleset*bestcolors;
    int besttotwct;
    int nbbest;
    
    const Scheduleset*debugcolors;
    int ndebugcolors;
    int opt_track;
    
    
    //maxiterations and retireage
    int maxiterations;
    int retirementage;
    
    //Branches
    int branch_job;
    int completiontime;
    
    wctdata* parent;
    wctdata* duetime_child;
    int nsame;
    wctdata* releasetime_child;
    int ndiff;

    //heur_diving
    heur_diving heur_data;
    
    char pname[MAX_PNAME_LEN];
};

/**
 * wct problem data type
 */


typedef struct wctproblem wctproblem;

struct wctproblem{
    wctparms parms;
    wctdata root_pd;
    
    int nwctdata;
    int global_upper_bound;
    int first_upper_bound;
    int global_lower_bound;
    int first_lower_bound;
    double first_rel_error;
    double rel_error;
    int maxdepth;

    enum  {
        no_sol = 0,     
        feasible   = 1,
        meta_heur  = 2,
        optimal    = 3
    } status;

    /* Scatter search*/
    SS scatter_search;
    
    /* All stable sets*/
    Scheduleset*initsets;
    int nbinitsets;
    int gallocated;
    /* Best Solution*/
    Scheduleset*bestschedule;
    int nbestschedule;
    /*heap variables*/
    pmcheap *br_heap;
    int mult_key;
    int found;
    /*Cpu time measurement*/
    CCutil_timer tot_cputime;
    CCutil_timer tot_branching_cputime;
    CCutil_timer tot_lb;
    CCutil_timer tot_lb_cpu_time;
    CCutil_timer tot_scatter_search;
    CCutil_timer tot_pricing;
    double real_time;
};

/*Initialize pmc data*/
void wctdata_init(wctdata* pd);
int set_id_and_name(wctdata* pd, int id, const char* fname);
int wctdata_init_unique(wctdata *pd, int id, const char *name);
int lp_build(wctdata *pd);

/*Free the wctdata*/
void lpwctdata_free(wctdata *pd);
void children_data_free(wctdata *pd);
void wctdata_free(wctdata *pd);
#ifdef __cplusplus
}
#endif

#endif
