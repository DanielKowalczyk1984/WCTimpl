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

enum wctBranchStrat{
    min_strategy    = 0,
    no_branching    = min_strategy,
    min_lb_strategy = 1,
    dfs_strategy    = 2,
    max_strategy    = 3,
};

enum RoundingStrat{
    min_rounding      = 0,
    neighbor_rounding = min_rounding,
    uniform_rounding  = 1,
    no_rounding       = 2,
    max_rounding      = 3,
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


typedef struct wctparms{
    int init_upper_bound;
    int branching_strategy;
    int rounding_strategy;
    int parallel_branching;
    int nb_feas_sol;
    int combine_method;
    int scatter_search;
    int branchandbound;
    int symmetrybreaking;
    
    int delete_elists;
    int delete_cclasses;
    
    char *jobfile;
    char *outfile;
    char *cclasses_infile;
    char *cclasses_outfile;
    char *color_infile;
    char *backupdir;
    
    int upper_bounds_only;
    double branching_cpu_limit;
    double scatter_search_cpu_limit;
    
    int nmachines;
} wctparms;

/*Initialization and free memory*/
void wctparms_init(wctparms *parms);
void wctparms_free(wctparms *parms);

/*Functions for setting some parameters*/
int wctparms_set_init_upper_bound(wctparms* parms,int bound);
int wctparms_set_parallel(wctparms *parms,int parallel);
int wctparms_set_branching_cpu_limit(wctparms *parms, double limit);
int wctparms_set_branching_strategy(wctparms *parms, int strategy);
int wctparms_set_rounding_strategy(wctparms *parms, int strategy);
int wctparms_set_nb_feas_sol(wctparms *parms,int nb_sol);
int wctparms_set_combine_method(wctparms *parms,int combine_method);
int wctparms_set_scatter_search(wctparms* parms,int scatter);
int wctparms_set_branchandbound(wctparms *parms,int bound);
int wctparms_set_scatter_search_cpu_limit(wctparms *parms,double limit);

/*Functions for defining the filesname*/
int wctparms_set_outfile(wctparms *parms, const char *fname);
int wctparms_set_file(wctparms *parms, const char *fname);
int wctparms_set_backupdir(wctparms *parms, const char* fname);
int wctparms_set_cclasses_infile(wctparms *parms, const char *fname);
int wctparms_set_cclasses_outfile(wctparms *parms, const char *fname);
int wctparms_set_color_infile(wctparms *parms, const char *fname);
int wctparms_set_nmachines(wctparms *parms,int nmachines);





#endif

#ifdef __cplusplus
}
#endif
