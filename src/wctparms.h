////////////////////////////////////////////////////////////////
//                                                            //
//  wctparms.h                                                //
//  wct                                                       //
//                                                            //
//  Created by Daniel on 20/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#ifndef __WCTPARMS_H
#define __WCTPARMS_H

#ifdef __cplusplus
extern "C" {
#endif

enum BBSearchStrat {
    min_search_strategy    = 0,
    no_branching    = min_search_strategy,
    min_lb_strategy = 1,
    dfs_strategy    = 2,
    max_strategy    = 3,
};

enum CombineMethod {
    min_combine_method = 0,
    pm_combine_method  = min_combine_method,
    path_combine_method = 1,
};

enum BranchandBound {
    no = 0,
    yes = 1,
};

enum dualvariablesType {
    Dbl = 0,
    Int = 1,
};

enum stab_techniques {
    no_stab = 0,
    stab_wentgnes = 1,
    stab_dynamic = 2,
};

enum pricing_solver {
    min_solver = 0,
    bdd_solver = min_solver,
    zdd_solver = 1,
    DP_solver = 2,
};

enum construct_solutions {
    min_construct_solutions = 1,
    yes_construct = min_construct_solutions,
    no_construct = 0,
};

enum diving_heur {
    min_diving_heur = 0,
    select_diving_heur = 1,
};

enum test_ahv {
    min_use_test = 0,
    use_test = 1,
};

enum print {
    min_print_size = 0,
    use_print = 1,
};

enum BranchandBoundStrat {
    min_bb_strategy = 0,
    conflict_strategy = min_bb_strategy,
    ahv_strategy = 1,
    cbfs_conflict_strategy = 2,
    cbfs_ahv_strategy = 3,
};

enum Strong_Branching {
    min_strong_branching = 0,
    use_strong_branching = min_strong_branching,
    no_strong_branching = 1,
};



typedef struct wctparms {
    /**
     * General parameters
     */
    int init_upper_bound;
    int bb_search_strategy;
    int bb_branch_strategy;
    int strong_branching;
    int parallel_branching;
    int nb_feas_sol;
    double branching_cpu_limit;
    /**
     * scatter search
     */
    int combine_method;
    int scatter_search;
    double scatter_search_cpu_limit;
    /**
     * column generation
     */
    int branchandbound;
    int dual_var_type;
    int stab_technique;
    int solver;
    int construct;
    int diving_heuristic;
    int test_ahv;
    int print;

    int delete_elists;
    int delete_cclasses;

    char *jobfile;
    char *outfile;
    char *cclasses_infile;
    char *cclasses_outfile;
    char *color_infile;
    char *backupdir;

    int upper_bounds_only;
    int nmachines;
} wctparms;

/*Initialization and free memory*/
void wctparms_init(wctparms *parms);
void wctparms_free(wctparms *parms);

/*Functions for setting some parameters*/
int wctparms_set_init_upper_bound(wctparms *parms, int bound);
int wctparms_set_parallel(wctparms *parms, int parallel);
int wctparms_set_branching_cpu_limit(wctparms *parms, double limit);
int wctparms_set_search_strategy(wctparms *parms, int strategy);
int wctparms_set_branching_strategy(wctparms *parms, int strategy);
int wctparms_set_strong_branching(wctparms *parms, int strong);
int wctparms_set_nb_feas_sol(wctparms *parms, int nb_sol);

/**
 * column generation
 */
int wctparms_set_branchandbound(wctparms *parms, int bound);
int wctparms_set_dual_var_type(wctparms *parms, int);
int wctparms_set_stab_technique(wctparms *parms, int stab_technique);
int wctparms_set_solver(wctparms *parms, int solver);
int wctparms_set_construct(wctparms *parms, int construct);
int wctparms_set_diving_heuristic(wctparms *parms, int diving);
int wctparms_set_test_ahv(wctparms *parms, int use);
int wctparms_set_print(wctparms *parms, int print);

/**
 * scatter search
 */
int wctparms_set_scatter_search(wctparms *parms, int scatter);
int wctparms_set_combine_method(wctparms *parms, int combine_method);
int wctparms_set_scatter_search_cpu_limit(wctparms *parms, double limit);

/*Functions for defining the filesname*/
int wctparms_set_outfile(wctparms *parms, const char *fname);
int wctparms_set_file(wctparms *parms, const char *fname);
int wctparms_set_backupdir(wctparms *parms, const char *fname);
int wctparms_set_cclasses_infile(wctparms *parms, const char *fname);
int wctparms_set_cclasses_outfile(wctparms *parms, const char *fname);
int wctparms_set_color_infile(wctparms *parms, const char *fname);
int wctparms_set_nmachines(wctparms *parms, int nmachines);






#ifdef __cplusplus
}
#endif
#endif
