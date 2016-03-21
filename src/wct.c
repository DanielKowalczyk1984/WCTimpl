#include <string.h>
#include <math.h>
#include <assert.h>
#include "wct.h"
#include "wctparms.h"

int debug = 0;

static const double min_ndelrow_ratio = 0.5;

/*Information about debug*/
int dbg_lvl()
{
    return debug;
}
void set_dbg_lvl(int dbglvl)
{
    debug = dbglvl;
}

/**
 * Reading of the jobfile
 */
static int permute_nodes(int *invorder, int vcount, int *duration,
                         int *weight, int **durationlist, int **weightlist);
int read_problem(char *f, int *njobs, int **durationlist, int **weightlist)
{
    int i, val = 0;
    int nbjobs = 0,  prob = 0;
    int curduration, curweight, curjob = 0;
    int *duration = (int *) NULL;
    int *weight = (int *) NULL;
    double *ratio = (double *) NULL;
    char buf[256], *p;
    int bufsize;
    const char *delim = " \n";
    char *data = (char *) NULL;
    char *buf2 = (char *) NULL;
    FILE *in = (FILE *) NULL;
    int *perm = (int *) NULL;
    int *iperm = (int *)NULL;
    in = fopen(f, "r");

    if (!in) {
        fprintf(stderr, "Unable to open file %s\n", f);
        val = 1;
        goto CLEAN;
    }

    if (fgets(buf, 254, in) != NULL) {
        p = buf;

        if (p[0] == 'p') {
            if (prob) {
                fprintf(stderr, "ERROR: in this file we have to p lines\n");
                val = 1;
                goto CLEAN;
            }

            prob = 1;
            strtok(p, delim);
            data = strtok(NULL, delim);
            sscanf(data, "%d", &nbjobs);
            bufsize = 2 * nbjobs * (2 + (int) ceil(log((double)nbjobs + 10)));
            buf2 = (char *) CC_SAFE_MALLOC(bufsize, char);
            CCcheck_NULL_2(buf2, "Failed to allocate buf2");
            weight = CC_SAFE_MALLOC(nbjobs, int);
            CCcheck_NULL_2(weight, "out of memory for weight");
            duration = CC_SAFE_MALLOC(nbjobs, int);
            CCcheck_NULL_2(duration, "Failed to allocate memory");
            ratio = CC_SAFE_MALLOC(nbjobs, double);
            CCcheck_NULL_2(ratio, "Failed to allocate memory");
        } else {
            fprintf(stderr, "File has to give first the number vertices and edges.\n");
            val = 1;
            goto CLEAN;
        }
    } else {
        val = 1;
        goto CLEAN;
    }

    while (fgets(buf2, bufsize, in) != (char *)NULL) {
        p = buf2;

        if (p[0] == 'p') {
            if (prob) {
                fprintf(stderr, "ERROR: in this file we have to p lines\n");
                val = 1;
                goto CLEAN;
            }
        } else if (p[0] == 'n') {
            if (!prob) {
                fprintf(stderr, "ERROR n before p in file\n");
                val = 1;
                goto CLEAN;
            }

            strtok(p, delim);
            data = strtok(NULL, delim);
            sscanf(data, "%d", &curweight);
            data = strtok(NULL, delim);
            sscanf(data, "%d", &curduration);
            weight[curjob] = curweight;
            duration[curjob] = curduration;
            data = strtok(NULL, delim);
            curjob++;
        }
    }

    perm = CC_SAFE_MALLOC(nbjobs, int);
    CCcheck_NULL_2(perm, "Failed to allocate memory");

    for (i = 0; i < nbjobs; i++) {
        perm[i] = i;
        ratio[i] = (double) duration[i] / (double) weight[i];
    }

    CCutil_double_perm_quicksort(perm, ratio, nbjobs);
    iperm = CC_SAFE_MALLOC(nbjobs, int);
    CCcheck_NULL_2(iperm, "Failed to allocate memory");

    for (i = 0; i < nbjobs; ++i) {
        iperm[perm[i]] = i;
    }

    permute_nodes(iperm, nbjobs, duration, weight, durationlist, weightlist);
    *njobs = nbjobs;
CLEAN:

    if (val) {
        CC_IFFREE(*durationlist, int);
        CC_IFFREE(*weightlist, int);
    }

    CC_IFFREE(weight, int);
    CC_IFFREE(duration, int);
    CC_IFFREE(buf2, char);
    CC_IFFREE(perm, int);
    CC_IFFREE(iperm, int);
    CC_IFFREE(ratio, double);

    if (in) {
        fclose(in);
    }

    return val;
}

static int permute_nodes(int *invorder, int njobs,  int *duration,
                         int *weights, int **durationlist, int **weightlist)
{
    int i, val = 0;
    int *idurationlist = (int *) NULL, *iweightlist = (int *) NULL;
    iweightlist = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(iweightlist, "out of memory for iweights");
    idurationlist = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(iweightlist, "out of memory for iweights");

    for (i = 0; i < njobs; i++) {
        iweightlist[invorder[i]] = weights[i];
        idurationlist[invorder[i]] = duration[i];
    }

    *durationlist =  idurationlist;
    *weightlist = iweightlist;
CLEAN:

    if (val) {
        CC_IFFREE(iweightlist, int);
        CC_IFFREE(idurationlist, int);
    }

    return val;
}

/*Functions for initialization of the problem and freeing the problem*/
void wctproblem_init(wctproblem *problem)
{
    /*B&B info*/
    problem->nwctdata = 0;
    problem-> global_upper_bound = INT_MAX;
    problem->first_lower_bound = 0;
    problem->global_lower_bound = 0;
    problem->first_upper_bound = INT_MAX;
    problem->status = no_sol;
    problem->rel_error = 1.0;
    problem->nbestschedule = 0;
    problem->bestschedule = (Scheduleset *)NULL;
    problem->maxdepth = 0;
    problem->nbinitsets = 0;
    problem->gallocated = 0;
    problem->initsets = (Scheduleset *)NULL;
    /*data of the problem*/
    wctdata_init(&(problem->root_pd));
    /*parms of the problem*/
    wctparms_init(&(problem->parms));
    /*heap initialization*/
    problem->br_heap = (pmcheap *) NULL;
    pmcheap_init(&(problem->br_heap), 1000);
    problem->found = 0;
    /*CPU timer initialisation*/
    CCutil_init_timer(&(problem->tot_cputime), "tot_cputime");
    CCutil_init_timer(&(problem->tot_scatter_search), "tot_scatter_search");
    CCutil_init_timer(&(problem->tot_branch_and_bound), "tot_branch_and_bound");
    CCutil_init_timer(&(problem->tot_lb_lp_root), "tot_lb_lp_root");
    CCutil_init_timer(&(problem->tot_build_dd), "tot_build_dd");
    CCutil_init_timer(&(problem->tot_lb_lp), "tot_lb_lp");
    CCutil_init_timer(&(problem->tot_lb), "tot_lb");
    CCutil_init_timer(&(problem->tot_pricing), "tot_pricing");
    CCutil_init_timer(&(problem->tot_lb_heur), "tot_lb_heur");
    /* Scatter sear init*/
    SS_init(&problem->scatter_search, 15, 10, 0.1);
    /** Solver */
    problem->solver = (PricerSolver *) NULL;
}

void wctproblem_free(wctproblem *problem)
{
    /*free the parameters*/
    wctparms_free(&(problem->parms));
    CC_IFFREE(problem->root_pd.duration, int);
    CC_IFFREE(problem->root_pd.weights, int);
    wctdata_free(&(problem->root_pd));
    /*free the heap*/
    pmcheap_free(problem->br_heap);
    problem->br_heap = (pmcheap *) NULL;
    Schedulesets_free(&(problem->initsets), &(problem->gallocated));
    Schedulesets_free(&(problem->bestschedule), &(problem->nbestschedule));
    SS_free(&problem->scatter_search);
    deletePricerSolver(problem->solver);
}

/*Functions for initialization and free the data*/
void wctdata_init(wctdata *pd)
{
    /*Initialization B&B data*/
    pd->id = -1;
    pd->depth = 0;
    pd->status = initialized;
    sprintf(pd->pname, "temporary");
    /*Initialization graph data*/
    pd->njobs = 0;
    pd->orig_node_ids = (int *) NULL;
    pd->duration = (int *)NULL;
    pd->weights = (int *)NULL;
    pd->duetime = (int *) NULL;
    pd->releasetime = (int *) NULL;
    pd->jobarray = (Job *)NULL;
    pd->H_max = 0;
    pd->H_min = 0;
    pd->upper_bound = INT_MAX;
    pd->lower_bound = INT_MAX;
    pd->dbl_safe_lower_bound = 0.0;
    pd->dbl_est_lower_bound = 0.0;
    pd->lower_scaled_bound = 1;
    pd->kpc_pi_scalef = 1;
    pd->LP_lower_bound = 0.0;
    /*Initialization  of the LP*/
    pd->LP = (wctlp *)NULL;
    pd->x = (double *)NULL;
    pd->coef = (double *) NULL;
    pd->pi = (double *) NULL;
    pd->kpc_pi = (int *)NULL;
    /**init stab data */
    pd->pi_in = (double *) NULL;
    pd->pi_out = (double *) NULL;
    pd->pi_sep = (double *) NULL;
    pd->subgradient = (double *) NULL;
    pd->subgradient_in = (double *) NULL;
    pd->alpha = 0.8;
    /*Initialization pricing_problem*/
    pd->solver = (PricerSolver *) NULL;
    pd->nnonimprovements = 0;
    /*Initialization of Scheduleset*/
    pd->ccount = 0;
    pd->cclasses = (Scheduleset *)NULL;
    pd->dzcount = 0;
    pd->gallocated = 0;
    pd->newsets = (Scheduleset *)NULL;
    pd->nnewsets = 0;
    pd->bestcolors = (Scheduleset *) NULL;
    pd->nbbest = 0;
    pd->debugcolors = (Scheduleset *) NULL;
    pd->ndebugcolors = 0;
    pd->opt_track = 0;
    /*Initialization max and retirement age*/
    pd->maxiterations = 1000000;
    pd->retirementage = 1000000;
    /*initialization of branches*/
    pd->branch_job = -1;
    pd->parent = (wctdata *) NULL;
    pd->duetime_child = (wctdata *) NULL;
    pd->nsame = 0;
    pd->releasetime_child = (wctdata *) NULL;
    pd->ndiff = 0;
    heur_init(pd);
}

void lpwctdata_free(wctdata *pd)
{
    if (pd->LP) {
        wctlp_free(&(pd->LP));
    }

    if (pd->coef) {
        free(pd->coef);
        pd->coef = (double *)NULL;
    }

    if (pd->pi) {
        free(pd->pi);
        pd->pi = (double *) NULL;
    }

    if (pd->x) {
        free(pd->x);
        pd->x = (double *)NULL;
    }

    if (pd->kpc_pi) {
        free(pd->kpc_pi);
        pd->kpc_pi = (int *)NULL;
    }

    if (pd->solver) {
        deletePricerSolver(pd->solver);
    }

    CC_IFFREE(pd->pi_out, double);
    CC_IFFREE(pd->pi_in, double);
    CC_IFFREE(pd->pi_sep, double);
    CC_IFFREE(pd->subgradient, double);
    CC_IFFREE(pd->subgradient_in, double);
    heur_free(pd);
    Schedulesets_free(&(pd->newsets), &(pd->nnewsets));
    Schedulesets_free(&(pd->cclasses), &(pd->gallocated));
    pd->ccount = 0;
}

void children_data_free(wctdata *pd)
{
    int i;

    for (i = 0; i < pd->nsame; ++i) {
        wctdata_free(&(pd->duetime_child[i]));
    }

    CC_IFFREE(pd->duetime_child, wctdata);

    for (i = 0; i < pd->ndiff; ++i) {
        wctdata_free(&(pd->releasetime_child[i]));
    }

    CC_IFFREE(pd->releasetime_child, wctdata);
    pd->nsame = pd->ndiff = 0;
}

void temporary_data_free(wctdata *pd)
{
    children_data_free(pd);
    lpwctdata_free(pd);
}

void wctdata_free(wctdata *pd)
{
    temporary_data_free(pd);
    Schedulesets_free(&(pd->bestcolors), &(pd->nbbest));
    CC_IFFREE(pd->releasetime, int);
    CC_IFFREE(pd->duetime, int);
    CC_IFFREE(pd->orig_node_ids, int);
    CC_IFFREE(pd->jobarray, Job);
}

/** Preprocess data*/
gint compare_readytime(gconstpointer a, gconstpointer b)
{
    const int *x = &(((const Job *)a)->processingime);
    const int *y = &(((const Job *)b)->processingime);
    return *x - *y;
}

int calculate_ready_due_times(Job *jobarray, int njobs, int nmachines, int Hmin)
{
    int i, j, val = 0;
    int *sumleft = (int *) NULL;
    int *sumright = (int *) NULL;
    int temp_duration, temp_weight;
    sumleft = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(sumleft, "Failed to, allocate memory");
    sumright = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(sumright, "Failed to allocate memory");
    sumleft[0] = 0;

    for (i = 1; i < njobs; ++i) {
        sumleft[i] = sumleft[i - 1] + jobarray[i].processingime;
    }

    sumright[njobs - 1] = jobarray[njobs - 1].processingime;

    for (i = njobs - 2 ; i >= 0; --i) {
        sumright[i] = sumleft[i + 1] + jobarray[i].processingime;
    }

    for (i = nmachines  ; i < njobs; ++i) {
        temp_duration = jobarray[i].processingime;
        temp_weight = jobarray[i].weight;
        GList *temp_list = (GList *) NULL;

        for (j = 0; j < i; ++j) {
            if ((jobarray[j].processingime <= temp_duration && jobarray[i].weight >= temp_weight)
                    || (sumleft[j] <= Hmin - sumright[i])) {
                temp_list = g_list_append(temp_list, jobarray + j);
            }
        }

        if (g_list_length(temp_list) > (guint)nmachines - 1) {
            temp_list = g_list_sort(temp_list, compare_readytime);
            GList *it;
            int len = g_list_length(temp_list);
            int counter = 0;

            for (it = temp_list; it && counter < len - nmachines + 1; it = it->next) {
                jobarray[i].releasetime += ((Job *)it->data)->processingime;
                counter++;
            }

            jobarray[i].releasetime = (int) ceil((double) jobarray[i].releasetime / (double) nmachines);
        }

        g_list_free(temp_list);
    }

    for (i = 0; i < njobs - 1; ++i) {
        temp_duration = jobarray[i].processingime;
        temp_weight = jobarray[i].weight;
        int sum = jobarray[i].processingime;

        for (j = i + 1; j < njobs; ++j) {
            if ((jobarray[j].processingime >= temp_duration && jobarray[i].weight <= temp_weight)
                    || (sumleft[j] >= Hmin - sumright[i])) {
                sum += jobarray[j].processingime;
            }
        }

        int delta = (int)((double)jobarray[i].duetime - ceil((double)sum / (double) nmachines)) + jobarray[i].processingime;

        if (delta < jobarray[i].duetime) {
            jobarray[i].duetime = delta;
        }
    }

CLEAN:
    CC_IFFREE(sumleft, int);
    CC_IFFREE(sumright, int);
    return val;
}

int calculate_Hmax(int *durations, int nmachines, int njobs)
{
    int i, max = 0, val = 0;
    double temp;

    for (i = 0; i < njobs; ++i) {
        val += durations[i];

        if (max < durations[i]) {
            max = durations[i];
        }
    }

    val += (nmachines - 1) * max;
    temp = (double) val;
    temp = temp / (double)nmachines;
    val = (int) ceil(temp);
    return val;
}

int calculate_Hmin(int *durations, int nmachines, int njobs, int *perm)
{
    int i, val = 0;
    double temp;

    for (i = 0; i < njobs; ++i) {
        val += durations[i];
    }

    for (i = 0; i < nmachines - 1; ++i) {
        val -= durations[perm[i]];
    }

    temp = (double) val;
    temp = temp / (double)nmachines;
    val = (int) ceil(temp);
    return val;
}

int Preprocessdata(wctdata *pd)
{
    int i, val = 0;
    int njobs = pd->njobs;
    int nmachines = pd->nmachines;
    Job *jobarray = (Job *) NULL;
    int *perm = (int *) NULL;
    jobarray = CC_SAFE_MALLOC(pd->njobs, Job);
    CCcheck_NULL_2(jobarray, "Failed to allocate memory");
    perm = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(perm, "Failed to allocate memory");
    pd->releasetime = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(pd->releasetime, "Failed to allocate releasetime");
    pd->duetime = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(pd->duetime, "Failed to allocate duetime");

    /** Initialize jobarray of rootnode */
    for (i = 0; i < njobs; ++i) {
        jobarray[i].weight = pd->weights[i];
        jobarray[i].processingime = pd->duration[i];
        jobarray[i].releasetime = 0;
        jobarray[i].job = i;
        perm[i] = i;
    }

    CCutil_int_perm_quicksort_0(perm, pd->duration, njobs);
    /** Calculate H_max */
    pd->H_max = calculate_Hmax(pd->duration, pd->nmachines, njobs);

    for (i = 0; i < njobs; ++i) {
        jobarray[i].duetime = pd->H_max;
    }

    printf("H_max = %d\n", pd->H_max);
    /** Calculate H_min */
    pd->H_min = calculate_Hmin(pd->duration, pd->nmachines, njobs, perm);
    printf("H_min = %d\n", pd->H_min);
    /** Calculate ready times and due times */
    calculate_ready_due_times(jobarray, njobs, nmachines, pd->H_min);
    pd->jobarray = jobarray;

    for (i = 0; i < njobs; ++i) {
        pd->releasetime[i] = jobarray[i].releasetime;
        pd->duetime[i] = jobarray[i].duetime;
    }

CLEAN:

    if (val) {
        CC_IFFREE(jobarray, Job);
    }

    CC_IFFREE(perm, int);
    return val;
}

/** Help function for column generation */
static void print_ages(wctdata *pd)
{
    int i;
    printf("AGES:");

    for (i = 0; i < pd->ccount; ++i) {
        printf(" %4d", pd->cclasses[i].age);
    }

    printf("\n");
}

int compute_objective(wctdata *pd)
{
    int val = 0;
    int i;
    pd->LP_lower_bound_dual = .0;

    /** compute lower bound with the dual variables */
    for (i = 0; i < pd->njobs; i++) {
        pd->LP_lower_bound_dual += (double) pd->pi[i];
    }

    pd->LP_lower_bound_dual += pd->nmachines * pd->pi[pd->njobs];
    /** Get the LP lower bound and compute the lower bound of WCT */
    val = wctlp_objval(pd->LP, &(pd->LP_lower_bound));
    CCcheck_val_2(val, "pmclp_objval failed");

    if (pd->lower_bound > (int) ceil(pd->LP_lower_bound)) {
        pd->lower_bound = (int) ceil(pd->LP_lower_bound);
    }

    //if (dbg_lvl() > 0) {
        printf("Current primal LP objective: %19.16f  (LP_dual-bound %19.16f, lowerbound = %d).\n", pd->LP_lower_bound, pd->LP_lower_bound_dual, pd->lower_bound);
    //}

CLEAN:
    return val;
}

MAYBE_UNUSED
void make_pi_feasible(wctdata *pd)
{
    int c;

    for (c = 0; c < pd->ccount; ++c) {
        int i;
        double colsum = .0;
        double newcolsum = .0;

        for (i = 0; i < pd->cclasses[c].count; ++i) {
            if (signbit(pd->pi[pd->cclasses[c].members[i]])) {
                pd->pi[pd->cclasses[c].members[i]] = 0.0;
            }

            colsum += pd->pi[pd->cclasses[c].members[i]];
            colsum = nextafter(colsum, DBL_MAX);
        }

        if (!signbit(pd->pi[pd->cclasses[i].members[i]])) {
            pd->pi[pd->cclasses[c].members[i]] = 0;
        }

        colsum +=  pd->pi[pd->cclasses[c].members[i]];

        if (colsum > pd->cclasses[c].totwct) {
            for (i = 0; i < pd->cclasses[c].count; ++i) {
                pd->pi[pd->cclasses[c].members[i]] /= colsum;
                pd->pi[pd->cclasses[c].members[i]] *= pd->cclasses[c].totwct;
                newcolsum += pd->pi[pd->cclasses[c].members[i]];
            }

            pd->pi[pd->cclasses[c].members[i]] /= colsum;
            pd->pi[pd->cclasses[c].members[i]] *= pd->cclasses[c].totwct;
            newcolsum += pd->pi[pd->cclasses[c].members[i]];

            if (dbg_lvl() > 1) {
                printf("Decreased column sum of %5d from  %30.20f to  %30.20f\n", c, colsum,
                       newcolsum);
            }
        }
    }
}

MAYBE_UNUSED
void make_pi_feasible_farkas_pricing(wctdata *pd)
{
    int c;

    for (c = 0; c < pd->ccount; ++c) {
        int i;
        double colsum = .0;
        double newcolsum = .0;

        for (i = 0; i < pd->cclasses[c].count; ++i) {
            if (signbit(pd->pi[pd->cclasses[c].members[i]])) {
                pd->pi[pd->cclasses[c].members[i]] = 0.0;
            }

            colsum += pd->pi[pd->cclasses[c].members[i]];
            colsum = nextafter(colsum, DBL_MAX);
        }

        colsum += pd->pi[pd->cclasses[c].members[i]];

        if (colsum > pd->cclasses[c].totwct) {
            for (i = 0; i < pd->cclasses[c].count; ++i) {
                pd->pi[pd->cclasses[c].members[i]] /= colsum;
                pd->pi[pd->cclasses[c].members[i]] *= pd->cclasses[c].totwct;
                newcolsum += pd->pi[pd->cclasses[c].members[i]];
            }

            pd->pi[pd->cclasses[c].members[i]] /= colsum;
            pd->pi[pd->cclasses[c].members[i]] *= pd->cclasses[c].totwct;
            newcolsum += pd->pi[pd->cclasses[c].members[i]];

            if (dbg_lvl() > 1) {
                printf("Decreased column sum of %5d from  %30.20f to  %30.20f\n", c, colsum,
                       newcolsum);
            }
        }
    }
}

static void reset_ages(Scheduleset *cclasses, int cccount)
{
    int i;

    for (i = 0; i < cccount; i++) {
        cclasses[i].age = 0;
    }
}

int add_newsets(wctdata *pd)
{
    int val = 0;
    Scheduleset *tmpsets = (Scheduleset *) NULL;
    int i;

    if (pd->nnewsets == 0) {
        return val;
    }

    reset_ages(pd->newsets, pd->nnewsets);

    if (pd->ccount + pd->nnewsets > pd->gallocated) {
        pd->gallocated *= 2;
        tmpsets = CC_SAFE_MALLOC(pd->gallocated, Scheduleset);
        CCcheck_NULL_2(tmpsets, "Failed to allocate memory to tmpsets");
        memcpy(tmpsets, pd->cclasses, pd->ccount * sizeof(Scheduleset));
        free(pd->cclasses);
        pd->cclasses = tmpsets;
        tmpsets = NULL;
    }

    memcpy(pd->cclasses + pd->ccount, pd->newsets,
           pd->nnewsets * sizeof(Scheduleset));
    pd->ccount += pd->nnewsets;

    for (i = pd->ccount; i < pd->gallocated; i++) {
        Scheduleset_init(pd->cclasses + i);
    }

CLEAN:

    if (val) {
        CC_IFFREE(pd->cclasses, Scheduleset);
    }

    CC_IFFREE(pd->newsets, Scheduleset);
    pd->nnewsets = 0;
    return val;
}

MAYBE_UNUSED
int double2int(int *kpc_pi, int *scalef, const double *pi, int vcount)
{
    int    i;
    double max_dbl_nweight = -DBL_MAX;
    double max_prec_dbl = exp2(DBL_MANT_DIG - 1);
    static const double max_mwiswt   = (double) INT_MAX;
    double dbl_scalef = CC_MIN(max_prec_dbl, max_mwiswt);
    dbl_scalef /= (double) vcount;

    for (i = 0; i < vcount; ++i) {
        max_dbl_nweight =
            CC_MAX(max_dbl_nweight, pi[i]);
    }

    dbl_scalef /= CC_MAX(1.0, max_dbl_nweight);
    dbl_scalef  = floor(dbl_scalef);
    *scalef  = (int) dbl_scalef;

    for (i = 0; i < vcount; ++i) {
        double weight = pi[i] * dbl_scalef;
        assert(weight < (double) INT_MAX);
        kpc_pi[i] = (int) weight;
    }

    return 0;
}

/** Branch and Price Algorithm */

static int test_theorem_ahv(wctdata *pd, const double x[], GList *branchjobs, int *min_completiontime)
{
    int val = 0;
    int i, j;
    min_completiontime = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(min_completiontime, "failed to allocate memory");
    fill_int(min_completiontime, pd->njobs, INT_MAX);

    for (j = 0; j < pd->njobs; ++j) {
        int found = 0;
        int C = 0 ;

        for (i = 0; i < pd->ccount; ++i) {
            //C = GPOINTER_TO_INT(g_hash_table_lookup(pd->cclasses[i].completiontime, GINT_TO_POINTER(j)));
            if (x[i] <= 0.0 || !C) {
                continue;
            }

            if (C < min_completiontime[j]) {
                min_completiontime[j] = C;
            }

            if (C != min_completiontime[j] && !found) {
                branchjobs = g_list_append(branchjobs, GINT_TO_POINTER(j));
                found = 1;
            }
        }
    }

CLEAN:

    if (val) {
        CC_IFFREE(min_completiontime, int);
        g_list_free(branchjobs);
    }

    return val;
}

static int grab_integral_solution_ahv(wctdata *pd, int *completion_time)
{
    int val = 0;
    int incumbent;
    double testincumbent;
    int *perm = (int *) NULL;
    pmcheap *heap = (pmcheap *) NULL;
    solution *sol = (solution *) NULL;
    int i;
    perm = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(perm, "Failed to allocate memory");

    for (i = 0; i < pd->njobs; ++i) {
        perm[i] = i;
    }

    sol = CC_SAFE_MALLOC(1, solution);
    CCcheck_NULL_2(sol, "Failed to allocate memory")
    solution_init(sol);
    solution_alloc(sol, pd->nmachines, pd->njobs);
    val = pmcheap_init(&heap, pd->nmachines);

    for (i = 0; i < pd->nmachines; i++) {
        pmcheap_insert(heap, sol->part[i].completiontime, sol->part + i);
    }

    CCutil_int_perm_quicksort(perm, completion_time, pd->njobs);
    incumbent = 0;

    for (i = 0; i < pd->njobs; ++i) {
        Job *j = pd->jobarray + perm[i];
        partlist *part = (partlist *) pmcheap_min(heap);
        partlist_insert(part, sol->vlist, j);
        pmcheap_insert(heap, part->completiontime, part);
        incumbent += pd->jobarray[perm[i]].weight * (part->completiontime);
    }

    val = wctlp_objval(pd->LP, &testincumbent);
    CCcheck_val_2(val, "COLORlp_objval failed");
    Schedulesets_free(&(pd->bestcolors), &(pd->nbbest));
    pd->bestcolors = (Scheduleset *) realloc(pd->bestcolors,
                     pd->nmachines * sizeof(Scheduleset));
    CCcheck_NULL_2(pd->bestcolors, "Failed to realloc pd->bestcolors");
    /** Construct solution with the help of theorem of AHV */
    val = partlist_to_Scheduleset(sol->part, sol->nmachines, pd->njobs, &(pd->bestcolors), &(pd->nbbest));
    CCcheck_val_2(val, "Failed in conversion");
    printf("Intermediate solution:\n");
    /** Print Solution */

    if (pd->besttotwct < pd->upper_bound) {
        pd->upper_bound = pd->besttotwct;
    }

    if (pd->upper_bound == pd->lower_bound) {
        pd->status = finished;
    }

CLEAN:
    CC_IFFREE(perm, int);
    pmcheap_free(heap);
    solution_free(sol);
    CC_IFFREE(sol, solution);
    return val;
}

static void adapt_global_upper_bound(wctproblem *problem, int new_upper_bound)
{
    if (problem->global_upper_bound > new_upper_bound) {
        problem->global_upper_bound = new_upper_bound;
    }
}

static int collect_releasetime_child(wctdata *cd)
{
    int rval = 0;
    int c;

    for (c = 0; c < cd->ndiff; ++c) {
        if (cd->releasetime_child[c].nbbest && (!cd->nbbest || cd->releasetime_child[c].besttotwct < cd->upper_bound)) {
            if (cd->besttotwct) {
                Schedulesets_free(&(cd->bestcolors), &(cd->nbbest));
            }

            cd->upper_bound = cd->besttotwct = cd->releasetime_child[c].besttotwct;
            cd->releasetime_child[c].nbbest = 0;
            cd->bestcolors = cd->releasetime_child[c].bestcolors;
            cd->releasetime_child[c].bestcolors = (Scheduleset *) NULL;
            /** Check if the solution is feasible, i.e. every job is covered */
        }
    }

    return rval;
}

static int collect_duetime_child(wctdata *cd)
{
    int rval = 0;
    int c;

    for (c = 0; c < cd->nsame; ++c) {
        if (cd->duetime_child[c].nbbest && (!cd->nbbest || cd->duetime_child[c].besttotwct < cd->upper_bound)) {
            if (cd->besttotwct) {
                Schedulesets_free(&(cd->bestcolors), &(cd->nbbest));
            }

            cd->upper_bound = cd->besttotwct = cd->duetime_child[c].besttotwct;
            cd->duetime_child[c].nbbest = 0;
            cd->bestcolors = cd->duetime_child[c].bestcolors;
            cd->duetime_child[c].bestcolors = (Scheduleset *) NULL;
            /** Check if the solution is feasible, i.e. every job is covered */
        }
    }

    return rval;
}

static int remove_finished_subtree(wctdata *child)
{
    int rval = 0;
    int i;
    wctdata *cd = (wctdata *) child;
    int all_same_finished = 1;
    int all_diff_finished = 1;

    while (cd) {
        for (i = 0; i < cd->nsame; ++i) {
            if (cd->duetime_child[i].status < finished) {
                all_same_finished = 0;
                break;
            }
        }

        if (cd->nsame && all_same_finished) {
            rval = collect_duetime_child(cd);
            CCcheck_val_2(rval, "Failed in collect_same_children");

            for (i = 0;  i < cd->nsame; ++i) {
                wctdata_free(cd->duetime_child + i);
            }

            free(cd->duetime_child);
            cd->duetime_child = (wctdata *) NULL;
            cd->nsame = 0;
        }

        for (i = 0; i < cd->ndiff; ++i) {
            if (cd->releasetime_child[i].status < finished) {
                all_diff_finished = 0;
                break;
            }
        }

        if (cd->ndiff && all_diff_finished) {
            rval = collect_releasetime_child(cd);
            CCcheck_val_2(rval, "Failed in collect_diff_children");

            for (i = 0;  i < cd->ndiff; ++i) {
                wctdata_free(cd->releasetime_child + i);
            }

            free(cd->releasetime_child);
            cd->releasetime_child = (wctdata *) NULL;
            cd->ndiff = 0;
        }

        if (!cd->duetime_child && !cd->releasetime_child) {
            cd->status = finished;
            cd = cd->parent;
        } else {
            cd = (wctdata *) NULL;
        }
    }

CLEAN:
    return rval;
}

int set_id_and_name(wctdata *pd, int id, const char *fname)
{
    int val = 0;
    int sval = 0;
    CCcheck_NULL_2(pd, "np memory was allocated to pd");
    pd->id = id;
    sval = snprintf(pd->pname, MAX_PNAME_LEN, "%s", fname);

    if (sval < 0 || MAX_PNAME_LEN <= sval) {
        val = 1;
        CCcheck_val(val, "Failed to write pname")
    }

CLEAN:
    return val;
}

static int create_duetime(wctdata *parent_cd, int branch_job, int completiontime)
{
    int val = 0;
    wctdata    *pd = (wctdata *) NULL;
    pd = (wctdata *) CC_SAFE_MALLOC(1, wctdata);
    CCcheck_NULL_2(pd, "Failed to allocate pd");
    wctdata_init(pd);
    pd->depth = parent_cd->depth + 1;
    pd->duration = parent_cd->duration;
    pd->weights = parent_cd->weights;
    pd->jobarray = parent_cd->jobarray;
    parent_cd->nsame         = 1;
    parent_cd->duetime_child = pd;
    pd->branch_job = branch_job;
    pd->completiontime = completiontime;
    pd->upper_bound = parent_cd->upper_bound;
    pd->lower_bound = parent_cd->lower_bound;
    pd->dbl_safe_lower_bound = parent_cd->dbl_safe_lower_bound;
    pd->parent = parent_cd;
    /** adjusted release time and duetime */
    /** Construct PricerSolver */
    /** transfer feasible schedules */
    CCcheck_val_2(val, "Failed in transfer_same_cclasses");
CLEAN:

    if (val) {
        if (pd) {
            wctdata_free(pd);
            free(pd);
        }

        parent_cd->duetime_child = (wctdata *) NULL;
    }

    return val;
}

static int is_releasetime_child(wctdata *cd)
{
    int i;

    for (i = 0; cd->parent && i < cd->parent->ndiff; ++i) {
        if (cd == cd->parent->releasetime_child + i) {
            return 1;
        }
    }

    return 0;
}

static int recover_pricersolver(wctdata *cd)
{
    int val = 0;
    wctdata    **path  = (wctdata **) NULL;
    int            npath = 0;
    wctdata     *tmp_cd  = cd;
    wctdata     *root_cd = cd;
    int            ndiff   = 0;
    int            i;
    int           *new_orig_node_ids = (int *) NULL;

    while (tmp_cd) {
        npath++;
        tmp_cd = tmp_cd->parent;
    }

    path = CC_SAFE_MALLOC(npath, wctdata *);
    CCcheck_NULL_2(path, "Failed to allocate path.");
    tmp_cd = cd;
    i      = npath;

    while (tmp_cd) {
        i--;
        path[i] = tmp_cd;
        root_cd = tmp_cd;

        if (is_releasetime_child(tmp_cd)) {
            ndiff++;
        }

        tmp_cd = tmp_cd->parent;
    }

    assert(!path[0]->parent);

    if (!cd->orig_node_ids) {
        int v;
        new_orig_node_ids = CC_SAFE_MALLOC(root_cd->njobs, int);
        CCcheck_NULL_2(new_orig_node_ids, "Failed to allocate new_orig_node_ids.");

        for (v = 0; v < root_cd->njobs; ++v) {
            new_orig_node_ids [v] = v;
        }
    }

    for (i = 1; i < npath; ++i) {
        wctdata *cur_cd = path[i];

        if (is_releasetime_child(cur_cd)) {
        } else {
        }
    }

    if (!cd->orig_node_ids) {
        cd->orig_node_ids = new_orig_node_ids;
        new_orig_node_ids = (int *) NULL;
    }

CLEAN:
    CC_IFFREE(path, wctdata *);
    CC_IFFREE(new_orig_node_ids, int);
    return val;
}

static int create_releasetime(wctdata *parent_cd, int branch_job, int completiontime)
{
    int val = 0;
    wctdata    *pd = (wctdata *) NULL;
    pd = (wctdata *) CC_SAFE_MALLOC(1, wctdata);
    CCcheck_NULL_2(pd, "Failed to allocate pd");
    wctdata_init(pd);
    pd->depth = parent_cd->depth + 1;
    pd->duration = parent_cd->duration;
    pd->weights = parent_cd->weights;
    pd->jobarray = parent_cd->jobarray;
    parent_cd->nsame         = 1;
    parent_cd->releasetime_child = pd;
    pd->branch_job = branch_job;
    pd->completiontime = completiontime;
    pd->upper_bound = parent_cd->upper_bound;
    pd->lower_bound = parent_cd->lower_bound;
    pd->dbl_safe_lower_bound = parent_cd->dbl_safe_lower_bound;
    pd->parent = parent_cd;
    /** adjusted duetime */
    /** Construct PricerSolver */
    /** transfer feasible schedules */
    CCcheck_val_2(val, "Failed in transfer_same_cclasses");
CLEAN:

    if (val) {
        if (pd) {
            wctdata_free(pd);
            free(pd);
        }

        parent_cd->duetime_child = (wctdata *) NULL;
    }

    return val;
}

static int find_strongest_children(int *strongest_v1, wctdata *pd, wctproblem *problem, GList *branchnodes, int *completiontime)
{
    int    rval = 0;
    int    max_non_improving_branches  = 3; /* cd->njobs / 100 + 1; */
    int    remaining_branches          = max_non_improving_branches;
    double strongest_dbl_lb = 0.0;
    *strongest_v1 = -1;
    GList *it = branchnodes;

    while (it && (remaining_branches--)) {
        int v1 = GPOINTER_TO_INT(it->data);
        double dbl_child_lb;
        /* Create duetime and releasetime */
        rval = create_duetime(pd, v1, completiontime[v1]);
        CCcheck_val_2(rval, "Failed in create_duetime");
        rval = create_releasetime(pd, v1, completiontime[v1]);
        CCcheck_val_2(rval, "Failed in create_differ");
        pd->duetime_child->maxiterations = 5;
        pd->releasetime_child->maxiterations = 5;
        compute_lower_bound(problem, pd->duetime_child);
        compute_lower_bound(problem, pd->releasetime_child);
        dbl_child_lb = (pd->duetime_child->dbl_safe_lower_bound < pd->releasetime_child->dbl_safe_lower_bound) ?
                       pd->duetime_child->dbl_safe_lower_bound : pd->releasetime_child->dbl_safe_lower_bound;
        wctdata_free(pd->duetime_child);
        pd->nsame = 0;
        free(pd->duetime_child);
        pd->duetime_child = (wctdata *) NULL;
        wctdata_free(pd->releasetime_child);
        pd->ndiff = 0;
        free(pd->releasetime_child);
        pd->releasetime_child = (wctdata *) NULL;

        if (dbl_child_lb > strongest_dbl_lb) {
            strongest_dbl_lb = dbl_child_lb;
            *strongest_v1     = v1;
            remaining_branches = max_non_improving_branches;
        }
    }

CLEAN:
    return rval;
}

static void Scheduleset_unify(Scheduleset *cclasses, int *new_ccount, int ccount)
{
    int i;
    Scheduleset temp;
    Scheduleset_quicksort(cclasses, ccount, (*Scheduleset_less));
    *new_ccount = 0;
    i = 0;

    if (! ccount) {
        return;
    }

    /* Find first non-empty set */
    while (!cclasses[i].count) {
        Scheduleset_free(&(cclasses[i++]));
    }

    for (; i < ccount; ++i) {
        if (*new_ccount == 0
                || Scheduleset_less(&(cclasses[*new_ccount - 1]), &(cclasses[i]))) {
            (*new_ccount)++;

            if (*new_ccount  < i + 1) {
                Scheduleset_SWAP(&(cclasses[*new_ccount - 1]), & (cclasses[i]), &temp);
            }
        } else {
            Scheduleset_free(&(cclasses[i]));
        }
    }
}

int prune_duplicated_sets(wctdata *pd)
{
    int val = 0;
    int i, j;
    Scheduleset_unify(pd->cclasses, &(pd->ccount), pd->ccount);

    for (i = 0 ; i < pd->ccount; ++i) {
        if (dbg_lvl() > 1) {
            printf("TRANSSORT SET ");

            for (j = 0; j < pd->cclasses[i].count; ++j) {
                printf(" %d", pd->cclasses[i].members[j]);
            }

            printf("\n");
        }
    }

    return val;
}

int create_branches(wctdata *pd, wctproblem *problem)
{
    int val = 0;
    int result;
    int status;
    double *x = (double *)NULL;
    int strongest_v1 = -1;
    GList *branchjobs = (GList *) NULL;
    int *min_completiontime = (int *) NULL;

    if (!pd->LP) {
        val = build_lp(pd, 1);
        CCcheck_val_2(val, "Failed at build_lp");
    }

    if (!pd->ccount) {
        compute_lower_bound(problem, pd);
    }

    assert(pd->ccount != 0);
    x = CC_SAFE_MALLOC(pd->ccount, double);
    CCcheck_NULL_2(x, "Failed to allocate memory to x");
    val = wctlp_optimize(pd->LP, &status);
    CCcheck_val_2(val, "Failed at wctlp_optimize");
    val = wctlp_x(pd->LP, x, 0);
    CCcheck_val_2(val, "Failed at wctlp_x");
    CC_IFFREE(pd->x, double);
    pd->x = CC_SAFE_MALLOC(pd->ccount, double);
    CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
    memcpy(pd->x, x, pd->ccount * sizeof(double));
    val = test_theorem_ahv(pd, pd->x, branchjobs, min_completiontime);
    CCcheck_val_2(val, "Failed in test ahv");

    if (g_list_length(branchjobs) == 0) {
        printf("LP returned integral solution\n");
        val = grab_integral_solution_ahv(pd, min_completiontime);
        CCcheck_val_2(val, "Failed in grab_int_sol");
        assert(pd->status = finished);
        goto CLEAN;
    }

    val = heur_exec(problem, pd, &result);
    CCcheck_val_2(val, "Failed at heur_exec");
    val = find_strongest_children(&strongest_v1, pd, problem, branchjobs, min_completiontime);
    CCcheck_val_2(val, "Failed in find_strongest_children");
    val = create_duetime(pd, strongest_v1, min_completiontime[strongest_v1]);
    CCcheck_val(val, "Failed in create_same");
    val = set_id_and_name(pd->duetime_child, problem->nwctdata++, pd->pname);
    CCcheck_val_2(val, "Failed in set_id_and_name");
    val = compute_lower_bound(problem, pd->duetime_child);
    CCcheck_val_2(val, "Failed in compute_lower_bound");
    val = create_releasetime(pd, strongest_v1, min_completiontime[strongest_v1]);
    CCcheck_val_2(val, "Failed in create_differ");
    val = set_id_and_name(pd->releasetime_child, problem->nwctdata++, pd->pname);
    CCcheck_val_2(val, "Failed in set_id_and_name");
    val = compute_lower_bound(problem, pd->releasetime_child);
    CCcheck_val_2(val, "Failed in compute_lower_bound");
CLEAN:
    lpwctdata_free(pd);
    g_list_free(branchjobs);
    CC_IFFREE(min_completiontime, int);
    CC_IFFREE(x, double);
    return val;
}

double safe_lower_dbl(int numerator, int denominator)
{
    double result;
    double denom_mult;
    denom_mult = denominator;
    denom_mult = nextafter(denom_mult, DBL_MAX);
    denom_mult = 1 / denom_mult;
    denom_mult = nextafter(denom_mult, -DBL_MAX);
    result = (double) numerator * denom_mult;
    result = nextafter(result, -DBL_MAX);
    return result;
}

static int trigger_lb_changes(wctdata *child)
{
    int val = 0;
    int i;
    int new_lower_bound = child->lower_bound;
    wctdata *pd = (wctdata *) child->parent;

    while (pd) {
        for (i = 0;  i < pd->nsame; ++i) {
            if (pd->duetime_child[i].lower_bound < new_lower_bound) {
                new_lower_bound = pd->duetime_child[i].lower_bound;
            }
        }

        for (i = 0;  i < pd->ndiff; ++i) {
            if (pd->releasetime_child[i].lower_bound < new_lower_bound) {
                new_lower_bound = pd->releasetime_child[i].lower_bound;
            }
        }

        if (new_lower_bound > pd->lower_bound) {
            if (! pd->parent) {   /* i.e. pd == root_cd */
                time_t current_time;
                char   current_timestr[40] = "";
                (void) time(&current_time);
                strftime(current_timestr, 39, "%c", localtime(&current_time));
                printf("Lower bound increased from %d to %d (%s). \n",
                       pd->lower_bound, new_lower_bound, current_timestr);
            }

            pd->lower_bound = new_lower_bound;
            pd = pd->parent;
        } else {
            pd = (wctdata *) NULL;
        }
    }

    return val;
}

int skip_wctdata(wctdata *pd, wctproblem *problem)
{
    pmcheap *br_heap = problem->br_heap;

    if (dbg_lvl()) {
        printf("Skipping with lb %d and ub %d at depth %d (id = %d, "
               "opt_track = %d, unprocessed nodes = %d).\n",
               pd->lower_bound, pd->upper_bound,
               pd->depth,
               pd->id, pd->opt_track, pmcheap_size(br_heap));
    }

    pd->status = finished;
    return 0;
}

int insert_into_branching_heap(wctdata *pd, wctproblem *problem)
{
    int val = 0;
    int heap_key = 0;

    switch (problem->parms.branching_strategy) {
        case dfs_strategy:
            if (pd->parent) {
                heap_key = (int)(pd->dbl_est_lower_bound) - pd->depth * 1000000 ;
            }

            break;

        case min_lb_strategy:
        default:
            heap_key = (int)(pd->dbl_est_lower_bound * problem->mult_key) - pd->depth -
                       pd->id % 2;
    }

    int lb = (int) ceil(pd->dbl_est_lower_bound);

    if (lb < problem->global_upper_bound) {
        if (dbg_lvl()) {
            printf("Inserting into branching heap with lb %d and ub %d at depth %d (id = %d) heap_key = %d\n",
                   pd->lower_bound, pd->upper_bound, pd->depth, pd->id, heap_key);
        }

        val = pmcheap_insert(problem->br_heap, heap_key, (void *)pd);
        CCcheck_val(val, "Failed at pmcheap_insert");
    } else {
        skip_wctdata(pd, problem);
    }

    trigger_lb_changes(pd);
    CCcheck_val_2(val, "Failed in trigger_lb_changes");
CLEAN:
    return val;
}

int branching_msg(wctdata *pd, wctproblem *problem)
{
    pmcheap *br_heap = problem->br_heap;

    if (pd->lower_bound < pd->upper_bound) {
        CCutil_suspend_timer(&problem->tot_cputime);
        printf("Branching with lb %d (est. %f and %f) at depth %d (id = %d, "
               "time = %f, unprocessed nodes = %d, vcount = %d, upper bound = %d, lower bound = %d).\n",
               pd->lower_bound, pd->dbl_safe_lower_bound, pd->LP_lower_bound,
               pd->depth,
               pd->id, problem->tot_cputime.cum_zeit, pmcheap_size(br_heap), pd->njobs, problem->global_upper_bound, problem->global_lower_bound);
        CCutil_resume_timer(&problem->tot_cputime);
    }

    return 0;
}

int sequential_branching(wctproblem *problem)
{
    int val = 0;
    wctdata *pd;
    pmcheap *br_heap = problem->br_heap;
    wctparms *parms = &(problem->parms);
    printf("ENTERED SEQUANTIAL BRANCHING:\n");
    CCutil_suspend_timer(&problem->tot_branch_and_bound);

    while ((pd = (wctdata *) pmcheap_min(br_heap))
            && problem->tot_branch_and_bound.cum_zeit < parms->branching_cpu_limit) {
        CCutil_resume_timer(&problem->tot_branch_and_bound);
        int i;
        pd->upper_bound = problem->global_upper_bound;

        if (pd->lower_bound >= pd->upper_bound || pd->status == infeasible) {
            skip_wctdata(pd, problem);
            remove_finished_subtree(pd);
        } else {
            branching_msg(pd, problem);
            /** Construct PricerSolver */
            recover_pricersolver(pd);

            if (problem->maxdepth < pd->depth) {
                problem->maxdepth = pd->depth;
            }

            val = create_branches(pd, problem);
            CCcheck_val_2(val, "Failed at create_branches");

            for (i = 0; i < pd->nsame; i++) {
                val = insert_into_branching_heap(&(pd->duetime_child[i]), problem);
                CCcheck_val_2(val, "Failed in insert_into_branching_heap");
            }

            for (i = 0; i < pd->ndiff; i++) {
                val = insert_into_branching_heap(pd->releasetime_child + i, problem);
                CCcheck_val_2(val, "Faield at insert_into_branching_heap");
            }

            assert(pd->lower_bound <= pd->upper_bound);
            adapt_global_upper_bound(problem, pd->upper_bound);

            if (pd->upper_bound == pd->lower_bound) {
                remove_finished_subtree(pd);
            }

            /** Check for integer solutions */
        }

        CCutil_suspend_timer(&problem->tot_branch_and_bound);
    }

    CCutil_resume_timer(&problem->tot_branch_and_bound);

    if (pd) {
        printf("Branching timeout of %f second reached\n",
               parms->branching_cpu_limit);
    }

    children_data_free(&problem->root_pd);
CLEAN:
    return val;
}

static int grow_ages(wctdata *pd)
{
    int val = 0;
    int i;
    int *cstat;
    cstat = (int *)CC_SAFE_MALLOC(pd->ccount, int);
    CCcheck_NULL_2(cstat, "Failed to allocate cstat");
    val = wctlp_basis_cols(pd->LP, cstat, 0);
    CCcheck_val_2(val, "Failed in pmclp_basis_cols");
    pd->dzcount = 0;

    for (i = 0; i < pd->ccount; ++i) {
        if (cstat[i] == wctlp_LOWER || cstat[i] == wctlp_FREE) {
            pd->cclasses[i].age++;

            if (pd->cclasses[i].age > pd->retirementage) {
                pd->dzcount++;
            }
        } else {
            pd->cclasses[i].age = 0;
        }
    }

    /*    printf("%d out of %d are older than %d.\n", pd->dzcount, pd->ccount,  */
    /*           pd->retirementage); */
CLEAN:
    CC_IFFREE(cstat, int);
    return val;
}

static int delete_old_cclasses(wctdata *pd)
{
    int val   = 0;
    int i;
    int min_numdel = pd->njobs * min_ndelrow_ratio;
    int first_del = -1;
    int last_del  = -1;
    /** pd->dzcount can be deprecated! */
    pd->dzcount = 0;

    for (i = 0; i < pd->ccount; ++i) {
        if (pd->cclasses[i].age > pd->retirementage) {
            pd->dzcount++;
        }
    }

    if (pd->dzcount > min_numdel) {
        int       new_ccount = 0;
        Scheduleset *new_cclasses = (Scheduleset *) NULL;
        assert(pd->gallocated >= pd->ccount);
        new_cclasses = CC_SAFE_MALLOC(pd->gallocated, Scheduleset);
        CCcheck_NULL_2(new_cclasses, "Failed to allocate new_cclasses");

        for (i = 0; i < pd->gallocated; ++i) {
            Scheduleset_init(new_cclasses + i);
        }

        for (i = 0; i < pd->ccount; ++i) {
            if (pd->cclasses[i].age <= pd->retirementage) {
                if (first_del != -1) {
                    /** Delete recently found deletion range.*/
                    val = wctlp_deletecols(pd->LP, first_del, last_del);
                    CCcheck_val_2(val, "Failed in pmclp_deletecols");
                    first_del = last_del = -1;
                }

                memcpy(new_cclasses + new_ccount, pd->cclasses + i, sizeof(Scheduleset));
                new_ccount++;
            } else {
                Scheduleset_free(pd->cclasses + i);

                if (first_del == -1) {
                    first_del = new_ccount;
                    last_del  = first_del;
                } else {
                    last_del++;
                }
            }
        }

        if (first_del != -1) {
            /** Delete the final range. This can occur if the last
             element is to be deleted, e.g. when no further columns were
             added in a B&B branch.
             */
            wctlp_deletecols(pd->LP, first_del, last_del);
            CCcheck_val_2(val, "Failed in pmclp_deletecols");
        }

        assert(pd->dzcount == pd->ccount - new_ccount);
        CC_IFFREE(pd->cclasses, Scheduleset);
        pd->cclasses = new_cclasses;
        pd->ccount   = new_ccount;
        //if (dbg_lvl() > 0) {
        printf("Deleted %d out of %d columns with age > %d.\n",
               pd->dzcount, pd->dzcount + pd->ccount, pd->retirementage);
        //}
        getchar();
        pd->dzcount = 0;
    }

CLEAN:
    return val;
}

int build_lp(wctdata *pd, int construct)
{
    int val = 0;
    int i, j;
    int  *covered = (int *)NULL;
    int counter = 0;
    covered = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(covered, "Failed to allocate memory to covered");
    fill_int(covered, pd->njobs, 0);
    val = wctlp_init(&(pd->LP), NULL);
    CCcheck_val_2(val, "wctlp_init failed");

    for (i = 0; i < pd->njobs; i++) {
        val = wctlp_addrow(pd->LP, 0, (int *)NULL, (double *)NULL, wctlp_EQUAL, 1.0,
                           (char *)NULL);
        CCcheck_val_2(val, "Failed wctlp_addrow");
    }

    wctlp_addrow(pd->LP, 0  , (int *)NULL, (double *) NULL, wctlp_EQUAL, (double)pd->nmachines , NULL);
    pd->coef = (double *)CCutil_reallocrus(pd->coef, (pd->njobs + 1) * sizeof(double));
    CCcheck_NULL_2(pd->coef, "out of memory for coef");
    fill_dbl(pd->coef, pd->njobs + 1, 1.0);

    /** add constraint about number of machines */
    if (construct) {
        for (i = 0; i < pd->ccount; i++) {
            val = wctlp_addcol(pd->LP, pd->cclasses[i].count + 1,
                               pd->cclasses[i].members,
                               pd->coef, pd->cclasses[i].totwct, 0.0, 1.0, wctlp_CONT, NULL);

            if (val) {
                wctlp_printerrorcode(val);
            }

            CCcheck_val_2(val, "pmclp_addcol failed");

            for (j = 0; j < pd->cclasses[i].count && counter < pd->njobs; j++) {
                if (!covered[pd->cclasses[i].members[j]]) {
                    covered[pd->cclasses[i].members[j]] = 1;
                    counter++;
                }
            }
        }
    }

    if (counter < pd->njobs) {
        /** Farkas Pricing */
        for (i = 0; i < pd->njobs; i++) {
            if (!covered[i]) {
                pd->newsets = CC_SAFE_MALLOC(1, Scheduleset);
                Scheduleset_init(pd->newsets);
                pd->newsets[0].members = CC_SAFE_MALLOC(2, int);
                pd->nnewsets = 1;
                CCcheck_NULL_2(pd->newsets[0].members,
                               "Failed to allocate memory to pd->newsets->members");
                pd->newsets[0].count++;
                pd->newsets[0].members[0] = i;
                pd->newsets[0].members[1] = pd->njobs;
                pd->newsets[0].totwct = pd->weights[i] * pd->duration[i];
                pd->newsets[0].totweight = pd->duration[i];
                pd->newsets->age = 0;
                val = wctlp_addcol(pd->LP, 2, pd->newsets[0].members, pd->coef, pd->newsets[0].totwct, 0.0, 1.0,
                                   wctlp_CONT, NULL);
                CCcheck_val_2(val, "Failed in pmclp_addcol");

                if (pd->gallocated == 0 && pd->ccount == 0) {
                    pd->cclasses = CC_SAFE_MALLOC(pd->njobs, Scheduleset);

                    for (j = 0; j < pd->njobs; ++j) {
                        Scheduleset_init(pd->cclasses + i);
                    }

                    pd->gallocated = pd->njobs;
                }

                add_newsets(pd);
            }
        }
    }

    pd->pi = (double *)CCutil_reallocrus(pd->pi, (pd->njobs + 1) * sizeof(double));
    CCcheck_NULL_2(pd->pi, "Failed to allocate memory to pd->pi");
CLEAN:

    if (val) {
        wctlp_free(&(pd->LP));
        CC_IFFREE(pd->coef, double);
        CC_IFFREE(pd->pi, double);
    }

    CC_IFFREE(covered, int);
    return val;
}


int compute_lower_bound(wctproblem *problem, wctdata *pd)
{
    int  j, val = 0;
    int iterations = 0;
    int break_while_loop = 1;
    int    nnonimprovements     = 0;
    int status = GRB_LOADED;
    double cur_cputime;
    double last_lower_bound = DBL_MAX;
    wctparms *parms = &(problem->parms);
    /* Construction of solver*/
    CCutil_start_resume_time(&(problem->tot_build_dd));
    pd->solver = problem->solver;
    CCutil_suspend_timer(&(problem->tot_build_dd));
    /** Init table */

    if (dbg_lvl()) {
        printf("Starting compute_lower_bound with lb %d and ub %d at depth %d(id = %d, opt_track = %d)\n",
               pd->lower_bound, pd->upper_bound, pd->depth, pd->id, pd->opt_track);
    }

    CCutil_start_resume_time(&(problem->tot_lb_lp));

    /** Construct solutions if list of columns is empty */
    if (!pd->ccount && parms->construct) {
        solution *new_sol;
        new_sol = new_sol_init(pd->nmachines, pd->njobs);
        construct_wspt(pd->jobarray, pd->njobs, pd->nmachines, new_sol);
        partlist_to_Scheduleset(new_sol->part, pd->nmachines, pd->njobs, &(pd->cclasses), &(pd->ccount));
        solution_free(new_sol);
        assert(pd->gallocated >= pd->ccount);
    }

    if (!pd->LP) {
        val = build_lp(pd, parms->construct);
        CCcheck_val(val, "build_lp failed");
    }

    pd->pi_in = CC_SAFE_MALLOC(pd->njobs + 1, double);
    CCcheck_NULL_2(pd->pi_in, "Failed to allocate memory");
    fill_dbl(pd->pi_in, pd->njobs, 0.0);
    pd->eta_in = 0.0;
    pd->pi_out = CC_SAFE_MALLOC(pd->njobs + 1, double);
    CCcheck_NULL_2(pd->pi_out, "Failed to allocate memory");
    pd->pi_sep = CC_SAFE_MALLOC(pd->njobs + 1, double);
    CCcheck_NULL_2(pd->pi_sep, "Failed to allocate memory");
    pd->subgradient_in = CC_SAFE_MALLOC(pd->njobs + 1, double);
    CCcheck_NULL_2(pd->subgradient_in, "Failed to allocate memory");
    pd->subgradient = CC_SAFE_MALLOC(pd->njobs + 1, double);
    CCcheck_NULL_2(pd->subgradient, "Failed to allocate memory");

    /** Init alpha */
    switch (parms->stab_technique) {
        case stab_wentgnes:
            pd->alpha = 0.8;
            break;

        case stab_dynamic:
            pd->alpha = 0.0;
            break;

        case no_stab:
            break;
    }

    do {
        iterations++;

        /** delete old columns */
        if (pd->dzcount > pd->njobs * min_ndelrow_ratio && status == GRB_OPTIMAL) {
            val = delete_old_cclasses(pd);
        }

        /** Compute LP relaxation */
        cur_cputime = CCutil_zeit();
        val = wctlp_optimize(pd->LP, &status);
        CCcheck_val_2(val, "pmclp_optimize failed");
        cur_cputime = CCutil_zeit() - cur_cputime;

        if (dbg_lvl()) {
            printf("Simplex took %f seconds.\n", CCutil_zeit() - cur_cputime);
            fflush(stdout);
        }

        if (dbg_lvl() > 1) {
            print_ages(pd);
        }

        CCutil_start_resume_time(&problem->tot_pricing);

        switch (status) {
            case GRB_OPTIMAL:
                /** grow ages of the different columns */
                val = grow_ages(pd);
                CCcheck_val_2(val, "Failed in grow_ages");
                /** get the dual variables and make them feasible */
                val = wctlp_pi(pd->LP, pd->pi);
                CCcheck_val_2(val, "pmclp_pi failed");
                /** Compute the objective function */
                val = compute_objective(pd);
                CCcheck_val_2(val, "Failed in compute_objective");

                if (iterations < pd->maxiterations) {
                    /** nnonimprovements? */
                    last_lower_bound = pd->dbl_safe_lower_bound;

                    /** Solve the pricing problem */
                    switch (parms->stab_technique) {
                        case stab_wentgnes:
                            memcpy(pd->pi_out, pd->pi, sizeof(double)*pd->njobs + 1);
                            pd->eta_out = pd->LP_lower_bound;
                            val = solve_stab(pd, parms);
                            CCcheck_val_2(val, "Failed in solve_stab");
                            break;

                        case stab_dynamic:
                            memcpy(pd->pi_out, pd->pi, sizeof(double)*pd->njobs + 1);
                            pd->eta_out = pd->LP_lower_bound;
                            val = solve_stab_dynamic(pd, parms);
                            CCcheck_val_2(val, "Failed in solve_stab");
                            break;

                        case no_stab:
                            val = solve_pricing(pd, parms);
                            CCcheck_val_2(val, "Failed in solving pricing");
                            break;
                    }
                }

                break;

            case GRB_INFEASIBLE:
                /** get the dual variables and make them feasible */
                val = wctlp_pi(pd->LP, pd->pi);
                CCcheck_val_2(val, "pmclp_pi failed");

                if (iterations < pd->maxiterations) {
                    val = solve_farkas_dbl(pd);
                    CCcheck_val_2(val, "Failed in solving farkas");
                }

                break;
        }

        CCutil_suspend_timer(&problem->tot_pricing);

        for (j = 0; j < pd->nnewsets; j++) {
            val = wctlp_addcol(pd->LP, pd->newsets[j].count + 1, pd->newsets[j].members, pd->coef, pd->newsets[j].totwct, 0.0, 1.0, wctlp_CONT, NULL);
            CCcheck_val_2(val, "pmclp_addcol failed");
        }

        switch (status) {
            case GRB_OPTIMAL:
                switch (parms->stab_technique) {
                    case stab_wentgnes:
                    case stab_dynamic:
                        break_while_loop = (CC_OURABS(pd->eta_out - pd->eta_in) < 0.0001);
                        break;

                    case no_stab:
                        break_while_loop = (pd->nnewsets == 0 || nnonimprovements > 5);
                        break;
                }

                break;

            case GRB_INFEASIBLE:
                break_while_loop = (pd->nnewsets == 0);
                break;
        }
        printf("new col\n");
        getchar();
        add_newsets(pd);
        CCutil_suspend_timer(&(problem->tot_cputime));
        CCutil_resume_timer(&(problem->tot_cputime));
    } while ((iterations < pd->maxiterations)
             && !break_while_loop
             && problem->tot_cputime.cum_zeit <= problem->parms.branching_cpu_limit);

    if (iterations < pd->maxiterations && problem->tot_cputime.cum_zeit <= problem->parms.branching_cpu_limit) {
        switch (status) {
            case GRB_OPTIMAL:

                /** change status of problem */
                if (problem->status == no_sol) {
                    problem->status = lp_feasible;
                }

                val = wctlp_optimize(pd->LP, &status);
                CCcheck_val_2(val, "pmclp_optimize failed");
                val = wctlp_pi(pd->LP, pd->pi);
                CCcheck_val_2(val, "pmclp_pi failed");
                val = compute_objective(pd);
                CCcheck_val_2(val, "compute_objective failed");

                if (dbg_lvl() > 1) {
                    printf("Found lb = %d (%f) upper_bound = %d (id= %d, iterations = %d,opt_track = %d).\n",
                           pd->lower_bound, pd->LP_lower_bound, pd->upper_bound,
                           pd->id, iterations, pd->opt_track);
                }

                pd->x = CC_SAFE_MALLOC(pd->ccount, double);
                CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
                val = wctlp_x(pd->LP, pd->x, 0);
                CCcheck_val_2(val, "Failed in pmclp_x");
                pd->status = LP_bound_computed;
                break;

            case GRB_INFEASIBLE:
                pd->status = infeasible;
        }
    } else  {
        pd->status = LP_bound_estimated;
    }

    printf("iterations = %d\n", iterations);
    fflush(stdout);
    CCutil_suspend_timer(&(problem->tot_lb_lp));
CLEAN:
    pd->solver = (PricerSolver *) NULL;
    return val;
}

static int prefill_heap(wctdata *pd, wctproblem *problem)
{
    int val = 0;
    int insert_into_heap = 0;

    if (problem->nwctdata <= pd->id) {
        problem->nwctdata = pd->id + 1;
    }

    if (pd->status < LP_bound_computed) {
        printf("Found a node with LP not computed!\n");
        val = compute_lower_bound(problem, pd);
        CCcheck_val_2(val, "Failed at compute_lower_bound");
        insert_into_heap = 1;
    }

    if (pd->status < finished) {
        int i;

        if (!pd->nsame || !pd->ndiff) {
            insert_into_heap = 1;
        }

        for (i = 0; (!insert_into_heap) && i < pd->nsame; ++i) {
            if (pd->duetime_child[i].status < LP_bound_computed) {
                insert_into_heap = 1;
            }
        }

        for (i = 0; (!insert_into_heap) && i < pd->ndiff; ++i) {
            if (pd->releasetime_child[i].status < LP_bound_computed) {
                insert_into_heap = 1;
            }
        }
    }

    if (insert_into_heap) {
        val = insert_into_branching_heap(pd, problem);
        CCcheck_val_2(val, "Failed in insert_into_branching_heap");
        children_data_free(pd);
    } else {
        int i;

        for (i = 0; i < pd->nsame; ++i) {
            prefill_heap(pd->duetime_child + i, problem);
        }

        for (i = 0; i < pd->ndiff; ++i) {
            prefill_heap(pd->releasetime_child + i, problem);
        }
    }

CLEAN:
    return val;
}

int compute_schedule(wctproblem *problem)
{
    int val = 0;
    wctdata *root_pd = &(problem->root_pd);
    problem->mult_key = (double)(INT_MAX - 1) / root_pd->njobs;
    problem->first_upper_bound = problem->global_upper_bound;
    problem->first_lower_bound = problem->global_lower_bound;
    problem->first_rel_error = (double)(problem->global_upper_bound -
                                        problem->global_lower_bound) / ((double)problem->global_lower_bound);
    /** Transform columns into maximal schedule with  respect to the properties of optimal solutions */
    /** Add maximal schedule sets with respect to the properties of optimal solutions */
    prune_duplicated_sets(root_pd);

    if (root_pd->status >= LP_bound_computed) {
        val = prefill_heap(root_pd, problem);
        CCcheck_val(val, "Failed in prefill_heap");
    } else {
        CCutil_start_timer(&(problem->tot_lb_lp_root));
        val = compute_lower_bound(problem, root_pd);
        CCcheck_val_2(val, "Failed in compute_lower_bound");
        CCutil_stop_timer(&(problem->tot_lb_lp_root), 0);
        //val = insert_into_branching_heap( root_pd, problem );
        CCcheck_val_2(val, "insert_into_branching_heap failed");
    }

    CCutil_start_resume_time(&(problem->tot_branch_and_bound));
    //val = sequential_branching( problem );
    CCcheck_val(val, "Failed in sequential_branching");
    CCutil_suspend_timer(&(problem->tot_branch_and_bound));
CLEAN:
    return val;
}


