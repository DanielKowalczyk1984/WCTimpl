////////////////////////////////////////////////////////////////
//                                                            //
//  main.c                                                    //
//  PMC                                                       //
//                                                            //
//  Created by Daniel on 20/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#include "wct.h"
#include "defs.h"

#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static void usage(char *f)
{
    fprintf(stderr, "Usage %s: [-see below-] adjlist_file NrMachines\n", f);
    fprintf(stderr, "   -d      turn on the debugging\n");
    fprintf(stderr, "   -f n    Number of feasible solutions that have to be constructed\n");
    fprintf(stderr, "   -s int  Node selection: 0 = none, 1= minimum lower bound(default), 2 = DFS\n");
    fprintf(stderr, "   -l dbl  Cpu time limit for branching\n");
    fprintf(stderr, "   -L dbl  Cpu time limit for scatter search\n");
    fprintf(stderr, "   -C int  Combine method scatter search: 0 = Pathrelinking, 1 = PM\n");
    fprintf(stderr, "   -r int  Scatter search use: 0 = no scatter search(default), 1 = scatter search\n");
    fprintf(stderr, "   -B int  Branch and Bound use: 0 = no branch and bound(default), 1 = use branch and bound\n");
    fprintf(stderr, "   -S int  Stabilization technique: 0 = no stabilization(default), 1 = stabilization wentgnes\n");
    fprintf(stderr, "   -z int  Pricing solver technique: 0 = BDD(default), 1 = ZDD, 2 = DP\n");
    fprintf(stderr, "   -c int  Construct heuristically solutions: 0 = yes(default), 1 = no\n");
}


static int parseargs(int ac, char **av, wctparms *parms)
{
    int c;
    int val = 0;
    int debug = dbg_lvl();

    while ((c = getopt(ac, av, "dr:f:s:l:L:C:B:z:S:c:")) != EOF) {
        switch (c) {
            case 'd':
                ++(debug);
                set_dbg_lvl(debug);
                break;

            case 'r':
                val = wctparms_set_scatter_search(parms, atoi(optarg));
                CCcheck_val(val, "Failed cclasses_infile");
                break;

            case 'f':
                val = wctparms_set_nb_feas_sol(parms, atoi(optarg));
                CCcheck_val(val, "Failed number feasible solutions");
                break;

            case 's':
                val = wctparms_set_branching_strategy(parms, atoi(optarg));
                CCcheck_val(val, "Failed set_branching_strategy");
                break;

            case 'l':
                val = wctparms_set_branching_cpu_limit(parms, atof(optarg));
                CCcheck_val(val, "Failed wctparms_set_branching_cpu_limit");
                break;

            case 'L':
                val = wctparms_set_scatter_search_cpu_limit(parms, atof(optarg));
                CCcheck_val(val, "Failed wctparms_set_scatter_search_cpu_limit");
                break;

            case 'C':
                val = wctparms_set_combine_method(parms, atoi(optarg));
                CCcheck_val(val, "Failed wctparms_set_combine_method");
                break;

            case 'B':
                val = wctparms_set_branchandbound(parms, atoi(optarg));
                CCcheck_val(val, "Failed wctparms_set_branchandbound");
                break;

            case 'S':
                val = wctparms_set_stab_technique(parms, atoi(optarg));
                CCcheck_val(val, "Failed in wctparms_set_stab_technique");
                break;

            case 'z':
                val = wctparms_set_solver(parms, atoi(optarg));
                CCcheck_val(val, "Failed in wctparms_set_solver");
                break;

            case 'c':
                val = wctparms_set_construct(parms, atoi(optarg));
                CCcheck_val(val, "Failed in construct sol");
                break;

            default:
                usage(av[0]);
                val = 1;
                goto CLEAN;
        }
    }

    if (ac <= optind) {
        val = 1;
        goto CLEAN;
    } else {
        val = wctparms_set_file(parms, av[optind++]);
        CCcheck_val(val, "Failed in wctparms_set_file");

        if (ac <= optind) {
            val = 1;
            goto CLEAN;
        }

        val = wctparms_set_nmachines(parms, atoi(av[optind++]));
        CCcheck_val(val, "Failed in wctparms_set_nmachines");
    }

CLEAN:

    if (val) {
        usage(av[0]);
    }

    return val;
}

static int get_problem_name(char *pname, const char *efname)
{
    int    rval = 0;
    int    len = 0;
    const char *fname = strrchr(efname, '/');
    const char *lastdot = strrchr(efname, '.');

    if (!fname) {
        fname = efname;
    } else {
        fname++;
    }

    if (lastdot) {
        len = lastdot - fname + 1;
    } else {
        len = strlen(fname);
    }

    if (snprintf(pname, len, "%s", fname) < 0) {
        rval = 1;
    }

    printf("Extracted problem name %s\n", pname);
    return rval;
}

MAYBE_UNUSED
static int print_to_csv(wctproblem *problem)
{
    int val = 0;
    wctdata *pd = &(problem->root_pd);
    FILE *file = (FILE *)NULL;
    char filenm[128];
    GDate date;
    g_date_set_time_t(&date, time(NULL));
    problem->real_time = getRealTime() - problem->real_time;
    CCutil_stop_timer(&(problem->tot_cputime), 0);
    file = fopen(filenm, "a+");

    if (file == NULL) {
        printf("We couldn't open %s in %s at line %d\n", filenm, __FILE__, __LINE__);
        val = 1;
        goto CLEAN;
    }

    fseek(file, 0, SEEK_END);
    int size = ftell(file);

    if (size == 0) {
        fprintf(file, "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s\n", "NameInstance", "tot_Real_time", "tot_cputime",
                "tot_lb", "tot_feasible_cputime",
                "tot_lb_cpu_time",
                "tot_branch_and_bound", "tot_scatter_search", "tot_kpc",
                "rel_error", "status", "global_lower_bound", "global_upper_bound",
                "first_lower_bound", "first_upper_bound", "first_rel_error",
                "maxdepth");
    }

    fprintf(file, "%s;%f;%f;%f;%f;%f;%f;%f;%f;%d;%d;%d;%d;%d;%f;%d\n", pd->pname,
            problem->real_time,
            problem->tot_cputime.cum_zeit,
            problem->tot_lb.cum_zeit,
            problem->tot_lb_lp.cum_zeit,
            problem->tot_branch_and_bound.cum_zeit,
            problem->tot_scatter_search.cum_zeit,
            problem->tot_pricing.cum_zeit,
            problem->rel_error,
            problem->status,
            problem->global_lower_bound,
            problem->global_upper_bound,
            problem->first_lower_bound,
            problem->first_upper_bound,
            problem->first_rel_error,
            problem->maxdepth);
    fclose(file);
CLEAN:
    return val;
}

MAYBE_UNUSED
static int print_to_screen(wctproblem *problem)
{
    int val = 0;

    switch (problem->status) {
        case no_sol:
            printf("We didn't decide if this instance is feasible or infeasible\n");
            break;

        case feasible:
        case lp_feasible:
        case meta_heur:
            printf("The suboptimal schedule with relative error %f is given by:\n", (double)(problem->global_upper_bound - problem->global_lower_bound) / (problem->global_lower_bound));
            print_schedule(problem->bestschedule, problem->nbestschedule);
            break;

        case optimal:
            printf("The optimal schedule is given by:\n");
            print_schedule(problem->bestschedule, problem->nbestschedule);
            break;
    }

    printf("Compute_schedule took %f seconds(tot_scatter_search %f, tot_branch_and_bound %f, tot_lb_lp_root %f, tot_lb_lp %f, tot_lb %f, tot_pricing %f, tot_build_dd %f) and %f seconds in real time\n",
           problem->tot_cputime.cum_zeit,
           problem->tot_scatter_search.cum_zeit,
           problem->tot_branch_and_bound.cum_zeit,
           problem->tot_lb_lp_root.cum_zeit,
           problem->tot_lb_lp.cum_zeit,
           problem->tot_lb.cum_zeit,
           problem->tot_pricing.cum_zeit,
           problem->tot_build_dd.cum_zeit,
           problem->real_time);
    return val;
}

int main(int ac, char **av)
{
    int val = 0;
    int i;
    wctproblem problem;
    wctproblem_init(&problem);
    CCcheck_val(val, "Failed in wctproblem_init");
    wctparms *parms = &(problem.parms);
    wctdata *pd = &(problem.root_pd);
    val = program_header(ac, av);
    CCcheck_val(val, "Failed in programheader");
    CCutil_start_timer(&(problem.tot_cputime));
    double start_time = CCutil_zeit();
    wctdata_init(pd);
    pd->id = 0;
    problem.nwctdata = 1;
    val = parseargs(ac, av, parms);

    if (val) {
        goto CLEAN;
    }

    get_problem_name(pd->pname, parms->jobfile);

    if (dbg_lvl() > 1) {
        printf("Debugging turned on\n");
    }

    fflush(stdout);
    /** Reading and preprocessing the data */
    val  = read_problem(parms->jobfile, &(pd->njobs), &(pd->duration), &(pd->weights));
    pd->nmachines = parms->nmachines;
    CCcheck_val(val, "read_adjlist failed");
    pd->orig_node_ids = (int *)CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(pd->orig_node_ids, "No memory to allocated orig_node_ids\n");

    for (i = 0; i < pd->njobs; i++) {
        pd->orig_node_ids[i] = i;
    }

    Preprocessdata(pd);
    printf("Reading and preprocessing of the data took %f\n", CCutil_zeit() - start_time);
    /** Computing initial lowerbound */
    CCutil_start_timer(&(problem.tot_lb));
    problem.global_lower_bound = lowerbound_eei(pd->jobarray, pd->njobs, pd->nmachines);
    problem.global_lower_bound = CC_MAX(problem.global_lower_bound, lowerbound_cp(pd->jobarray, pd->njobs, pd->nmachines));
    problem.global_lower_bound = CC_MAX(problem.global_lower_bound, lowerbound_cw(pd->jobarray, pd->njobs, pd->nmachines));
    CCutil_stop_timer(&(problem.tot_lb), 0);
    printf("Computing lowerbound EEI, CP and CW took %f\n", problem.tot_lb.cum_zeit);
    /** Construction Pricersolver */
    CCutil_start_resume_time(&(problem.tot_build_dd));
    problem.solver = newSolver(pd->duration, pd->weights, pd->releasetime, pd->duetime, pd->njobs, pd->H_min, pd->H_max);
    CCutil_suspend_timer(&(problem.tot_build_dd));

    /** Construct Feasible solutions */
    switch (parms->construct) {
        case yes_construct:
            construct_feasible_solutions(&problem);
            break;

        case no_construct:
            break;
    }

    /** Compute Schedule with Branch and Price */
    compute_schedule(&problem);
    /** Print all the information to screen and csv */
    //print_to_csv(&problem);
    print_to_screen(&problem);
CLEAN:
    wctproblem_free(&problem);
    return val;
}
