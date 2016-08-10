#include "wctparms.h"
#include "util.h"
#include <string.h>
#include <limits.h>



void wctparms_init(wctparms *parms) {
    parms->init_upper_bound = INT_MAX;
    parms->parallel_branching = 0;
    parms->bb_branch_strategy = min_bb_strategy;
    parms->bb_search_strategy = min_search_strategy;
    parms->strong_branching = min_strong_branching;
    parms->nb_feas_sol = 10;
    parms->combine_method = min_combine_method;
    parms->scatter_search = 0;
    parms->branchandbound = no;
    parms->dual_var_type = Dbl;
    parms->stab_technique = no_stab;
    parms->solver = min_solver;
    parms->construct = min_construct_solutions;
    parms->diving_heuristic = min_diving_heur;
    parms->test_ahv = min_use_test;
    parms->print = min_print_size;
    parms->delete_elists = 1;
    parms->delete_cclasses = 0;
    parms->jobfile = (char *) NULL;
    parms->outfile = (char *) NULL;
    parms->cclasses_outfile = (char *) NULL;
    parms->cclasses_infile = (char *) NULL;
    parms->color_infile = (char *) NULL;
    parms->backupdir = (char *) NULL;
    parms->upper_bounds_only = 0;
    parms->branching_cpu_limit = 3600.0;
    parms->scatter_search_cpu_limit = 40.0;
}

void wctparms_free(wctparms *parms) {
    CC_IFFREE(parms->color_infile, char);
    CC_IFFREE(parms->jobfile, char);
    CC_IFFREE(parms->outfile, char);
    CC_IFFREE(parms->cclasses_outfile, char);
    CC_IFFREE(parms->cclasses_infile, char);
    CC_IFFREE(parms->backupdir, char);
}

static int copy_string(char **dst, const char *src) {
    int val = 0;
    int len;
    len = (int)strlen(src) + 1;
    CC_IFFREE(*dst, char);
    *dst = (char *) CC_SAFE_MALLOC(len, char);
    CCcheck_NULL_2(*dst, "Failed to allocate dst");
    strcpy(*dst, src);
CLEAN:
    return val;
}

/*Functions for setting the name of the files*/
int wctparms_set_outfile(wctparms *parms, const char *fname) {
    return copy_string(&(parms->outfile), fname);
}
int wctparms_set_file(wctparms *parms, const char *fname) {
    return copy_string(&(parms->jobfile), fname);
}
int wctparms_set_cclasses_infile(wctparms *parms, const char *fname) {
    return copy_string(&(parms->cclasses_infile), fname);
}
int wctparms_set_cclasses_outfile(wctparms *parms, const char *fname) {
    return copy_string(&(parms->cclasses_outfile), fname);
}
int wctparms_set_backupdir(wctparms *parms, const char *fname) {
    return copy_string(&(parms->backupdir), fname);
}
int wctparms_set_color_infile(wctparms *parms, const char *fname) {
    return copy_string(&(parms->color_infile), fname);
}


/*Functions for setting the branching parameters*/
int wctparms_set_init_upper_bound(wctparms *parms, int bound) {
    parms->init_upper_bound = bound;
    return 0;
}
int wctparms_set_parallel(wctparms *parms, int parallel) {
    parms->parallel_branching = parallel;
    return 0;
}
int wctparms_set_branching_cpu_limit(wctparms *parms, double limit) {
    parms->branching_cpu_limit = limit;
    return 0;
}
int wctparms_set_branching_strategy(wctparms *parms, int strategy) {
    parms->bb_branch_strategy = strategy;
    return 0;
}

int wctparms_set_search_strategy(wctparms *parms, int strategy) {
    parms->bb_search_strategy = strategy;
    return 0;
}

int wctparms_set_strong_branching(wctparms *parms, int strong) {
    parms->strong_branching = strong;
    return 0;
}

int wctparms_set_nmachines(wctparms *parms, int nmachines) {
    parms->nmachines = nmachines;
    return 0;
}

int wctparms_set_nb_feas_sol(wctparms *parms, int nb_sol) {
    parms->nb_feas_sol = nb_sol;
    return 0;
}

int wctparms_set_combine_method(wctparms *parms, int combine_method) {
    parms->combine_method = combine_method;
    return 0;
}

int wctparms_set_scatter_search(wctparms *parms, int scatter) {
    parms->scatter_search = scatter;
    return 0;
}

int wctparms_set_branchandbound(wctparms *parms, int bound) {
    parms->branchandbound = bound;
    return 0;
}

int wctparms_set_scatter_search_cpu_limit(wctparms *parms, double limit) {
    parms->scatter_search_cpu_limit = limit;
    return 0;
}

int wctparms_set_dual_var_type(wctparms *parms, int type) {
    parms->dual_var_type = type;
    return 0;
}

int wctparms_set_stab_technique(wctparms *parms, int stab_technique) {
    parms->stab_technique = stab_technique;
    return 0;
}

int wctparms_set_solver(wctparms *parms, int solver) {
    parms->solver = solver;
    return 0;
}

int wctparms_set_construct(wctparms *parms, int construct) {
    parms->construct = construct;
    return 0;
}

int wctparms_set_diving_heuristic(wctparms *parms, int diving) {
    parms->diving_heuristic = diving;
    return 0;
}

int wctparms_set_test_ahv(wctparms *parms, int use) {
    parms->test_ahv = use;
    return 0;
}

int wctparms_set_print(wctparms *parms, int print) {
    parms->print = print;
    return 0;
}

