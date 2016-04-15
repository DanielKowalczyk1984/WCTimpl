#include <string.h>
#include <math.h>
#include <assert.h>
#include "wct.h"
#include "wctparms.h"
#include "heurdiving.h"

int recover_elist(wctdata *pd);

/** help functions for conflict branching */
static int create_same_conflict(wctproblem *problem, wctdata *parent_pd, wctdata **child, int v1, int v2);
static int create_differ_conflict(wctproblem *problem, wctdata *parent_pd, wctdata **child, int v1, int v2);
static int collect_same_child_conflict(wctdata *cd);
static int collect_diff_child_conflict(wctdata *cd);
static int remove_finished_subtree_conflict(wctdata *child);
static int trigger_lb_changes_conflict(wctdata *child);
/** help functions for wide branching */
static int create_same_wide(wctproblem *problem, wctdata *parent_pd, int *wide_v1, int *wide_v2, int nbwide);
static int create_diff_wide(wctproblem *problem, wctdata *parent_pd, int v1, int v2, int nbwide);
static int collect_diff_child_wide(wctdata *cd);
static int collect_same_child_wide(wctdata *cd);
static int remove_finished_subtree_wide(wctdata *child);
/** help functions for ahv branching */
static int create_releasetime_ahv(wctproblem *problem, wctdata *parent_pd, wctdata **child, int branch_job, int completion_time);
static int create_duetime_ahv(wctproblem *problem, wctdata *parent_pd, wctdata **child, int branch_job, int completiontime);
static int collect_duetime_child_ahv(wctdata *pd);
static int collect_releasetime_child_ahv(wctdata *pd);
static int remove_finished_subtree_ahv(wctdata *pd);
static int trigger_lb_changes_ahv(wctdata *child);

static int print_size_to_csv(wctproblem *problem, wctdata *pd);
/** small submipping heuristic */
static int submiping(wctdata *cd);

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
static int permute_nodes(int *invorder, int njobs, int *duration,
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

/** Printing sizes of ZDD */
MAYBE_UNUSED
static int print_size_to_csv(wctproblem *problem, wctdata *pd)
{
    int val = 0;
    int size;
    wctdata *root_node = &(problem->root_pd);
    char filenm[128];
    FILE *file = (FILE *)NULL;
    GDate date;
    g_date_set_time_t(&date, time(NULL));
    sprintf(filenm, "SizeZDD_%d_%d_%d_%d_%d.csv", pd->nmachines, pd->njobs, date.day, date.month, date.year);
    file = fopen(filenm, "a+");

    if (file == NULL) {
        printf("We couldn't open %s in %s at line %d\n", filenm, __FILE__, __LINE__);
        val = 1;
        goto CLEAN;
    }

    fseek(file, 0, SEEK_END);
    size = ftell(file);

    if (size == 0) {
        fprintf(file, "%s;%s;%s;%s;%s\n", "NameInstance", "depth", "size ZDD",
                "nb_same",
                "nb_diff");
    }

    fprintf(file, "%s;%d;%lu;%d;%d\n", root_node->pname, pd->depth, get_datasize(pd->solver), pd->ecount_same, pd->ecount_differ);
    fclose(file);
CLEAN:
    return val;
}

/*Functions for initialization of the problem and freeing the problem*/
void wctproblem_init(wctproblem *problem)
{
    /** Job data */
    problem->duration = (int *) NULL;
    problem->weight = (int *) NULL;
    problem->releasetime = (int *) NULL;
    problem->duetime = (int *) NULL;
    problem->jobarray = (Job *) NULL;
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
}

void wctproblem_free(wctproblem *problem)
{
    int n = problem->root_pd.njobs;
    /*free the parameters*/
    wctparms_free(&(problem->parms));
    wctdata_free(&(problem->root_pd));
    /*free the heap*/
    pmcheap_free(problem->br_heap);
    problem->br_heap = (pmcheap *) NULL;
    Schedulesets_free(&(problem->initsets), &(problem->gallocated));
    Schedulesets_free(&(problem->bestschedule), &(problem->nbestschedule));
    SS_free(&problem->scatter_search);
    CC_IFFREE(problem->duetime, int);
    CC_IFFREE(problem->releasetime, int);
    CC_IFFREE(problem->duration, int);
    CC_IFFREE(problem->weight, int);

    for (int i = 0; i < n; i++) {
        if (problem->jobarray[i].due_list) {
            g_list_free(problem->jobarray[i].due_list);
        }

        if (problem->jobarray[i].ready_list) {
            g_list_free(problem->jobarray[i].ready_list);
        }
    }

    CC_IFFREE(problem->jobarray, Job);
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
    pd->partial_sol = 0.0;
    pd->rhs = (double *) NULL;
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
    pd->update = 1;
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
    /** ahv branching */
    pd->duetime_child = (wctdata *) NULL;
    pd->nduetime = 0;
    pd->releasetime_child = (wctdata *) NULL;
    pd->nreleasetime = 0;
    pd->branch_job = -1;
    pd->completiontime = 0;
    /** conflict branching */
    pd->elist_same = (int *) NULL;
    pd->ecount_same = 0;
    pd->elist_differ = (int *) NULL;
    pd->ecount_differ = 0;
    pd->same_children = (wctdata *) NULL;
    pd->nsame = 0;
    pd->diff_children = (wctdata *) NULL;
    pd->ndiff = 0;
    pd->v1 = -1;
    pd->v2 = -1;
    /** Wide branching */
    pd->v1_wide = (int *) NULL;
    pd->v2_wide = (int *) NULL;
    pd->nb_wide = 0;
    pd->same_children_wide = (wctdata **) NULL;
    pd->diff_children_wide = (wctdata **) NULL;
    heur_init(pd);
}

static void free_elist(wctdata *cd, wctparms *parms)
{
    if (cd->parent && parms->delete_elists) {
        CC_IFFREE(cd->elist_same, int);
        CC_IFFREE(cd->elist_differ, int);
        CC_IFFREE(cd->v1_wide, int);
        CC_IFFREE(cd->v2_wide, int);
        CC_IFFREE(cd->orig_node_ids, int);
    }
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

    CC_IFFREE(pd->pi_out, double);
    CC_IFFREE(pd->pi_in, double);
    CC_IFFREE(pd->pi_sep, double);
    CC_IFFREE(pd->subgradient, double);
    CC_IFFREE(pd->subgradient_in, double);
    CC_IFFREE(pd->rhs, double);
    CC_IFFREE(pd->releasetime, int);
    CC_IFFREE(pd->duetime, int);
    heur_free(pd);

    if (pd->solver) {
        freeSolver(pd->solver);
        pd->solver = (PricerSolver *) NULL;
    }

    Schedulesets_free(&(pd->newsets), &(pd->nnewsets));
    Schedulesets_free(&(pd->cclasses), &(pd->gallocated));
    pd->ccount = 0;
}

void children_data_free(wctdata *pd)
{
    int i;

    for (i = 0; i < pd->nsame; ++i) {
        wctdata_free(&(pd->same_children[i]));
    }

    CC_IFFREE(pd->same_children, wctdata);

    for (i = 0; i < pd->ndiff; ++i) {
        wctdata_free(&(pd->diff_children[i]));
    }

    CC_IFFREE(pd->diff_children, wctdata);
    pd->nsame = pd->ndiff = 0;
}

void temporary_data_free(wctdata *pd)
{
    children_data_free(pd);
    lpwctdata_free(pd);
}

void wctdata_free(wctdata *pd)
{
    Schedulesets_free(&(pd->bestcolors), &(pd->nbbest));
    temporary_data_free(pd);
    CC_IFFREE(pd->elist_same, int);
    CC_IFFREE(pd->elist_differ, int);
    CC_IFFREE(pd->v1_wide, int);
    CC_IFFREE(pd->v2_wide, int);
    CC_IFFREE(pd->orig_node_ids, int);
    CC_IFFREE(pd->v1_wide, int);
    CC_IFFREE(pd->v2_wide, int);
}

/** help functions for heap srong branching */
static int nodepair_ref_key(int v1, int v2)
{
    /* We store only the elements of the upper right triangle within the
     vcount x vcount matrix. */
    assert(v1 <= v2);
    return v2 * (v2 + 1) / 2 + v1;
}

/** compute row-index v1 and column-index v2 from array-index.*/
static void inodepair_ref_key(int *v1, int *v2, int index)
{
    *v2 = (int) floor(sqrt(2 * ((double)index) + 0.25) - 0.5);
    *v1 = index - (*v2 * (*v2 + 1) / 2);
}

static int x_frac(const double x)
{
    double mean = 0.5;
    double frac = fabs(x - mean);
    assert(frac <= 1.0);
    return (int)(frac * (double) INT_MAX);
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

    for (i = 0  ; i < njobs; ++i) {
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

            for (it = temp_list; it && counter < len - nmachines; it = it->next) {
                jobarray[i].releasetime += ((Job *)it->data)->processingime;
                counter++;
            }

            jobarray[i].releasetime = (int) ceil((double) jobarray[i].releasetime / (double) nmachines);
        }

        jobarray[i].due_list = temp_list;
    }

    for (i = 0; i < njobs; ++i) {
        temp_duration = jobarray[i].processingime;
        temp_weight = jobarray[i].weight;
        int sum = jobarray[i].processingime;
        GList *templist = (GList *) NULL;

        for (j = i + 1; j < njobs; ++j) {
            if ((jobarray[j].processingime >= temp_duration && jobarray[i].weight <= temp_weight)
                    || (sumleft[j] >= Hmin - sumright[i])) {
                sum += jobarray[j].processingime;
                templist = g_list_append(templist, jobarray + j);
            }
        }

        int delta = (int)((double)jobarray[i].duetime - ceil((double)sum / (double) nmachines)) + jobarray[i].processingime;

        if (delta < jobarray[i].duetime) {
            jobarray[i].duetime = delta;
        }

        jobarray[i].ready_list = templist;
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

int Preprocessdata(wctproblem *problem, wctdata *pd)
{
    int i, val = 0;
    int njobs = pd->njobs;
    int nmachines = pd->nmachines;
    Job *jobarray = (Job *) NULL;
    int *perm = (int *) NULL;
    int *release_time = (int *) NULL;
    int *due_time = (int *) NULL;
    jobarray = CC_SAFE_MALLOC(pd->njobs, Job);
    CCcheck_NULL_2(jobarray, "Failed to allocate memory");
    perm = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(perm, "Failed to allocate memory");
    release_time = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(release_time, "Failed to allocate releasetime");
    due_time = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(due_time, "Failed to allocate duetime");

    /** Initialize jobarray of rootnode */
    for (i = 0; i < njobs; ++i) {
        jobarray[i].weight = problem->weight[i];
        jobarray[i].processingime = problem->duration[i];
        jobarray[i].releasetime = 0;
        jobarray[i].duetime = 0;
        jobarray[i].job = i;
        jobarray[i].due_list = (GList *) NULL;
        jobarray[i].ready_list = (GList *) NULL;
        perm[i] = i;
    }

    CCutil_int_perm_quicksort_0(perm, problem->duration, njobs);
    /** Calculate H_max */
    pd->H_max = calculate_Hmax(problem->duration, pd->nmachines, njobs);

    for (i = 0; i < njobs; ++i) {
        jobarray[i].duetime = pd->H_max;
    }

    printf("H_max = %d\n", pd->H_max);
    /** Calculate H_min */
    pd->H_min = calculate_Hmin(problem->duration, pd->nmachines, njobs, perm);
    printf("H_min = %d\n", pd->H_min);
    /** Calculate ready times and due times */
    calculate_ready_due_times(jobarray, njobs, nmachines, pd->H_min);
    problem->jobarray = jobarray;

    for (i = 0; i < njobs; ++i) {
        release_time[i] = jobarray[i].releasetime;
        due_time[i] = jobarray[i].duetime;
    }

    pd->jobarray = jobarray;
    pd->releasetime = release_time;
    pd->duetime = due_time;
    pd->duration = problem->duration;
    pd->weights = problem->weight;
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

int compute_objective(wctdata *pd, wctparms *parms)
{
    int val = 0;
    int i;
    pd->LP_lower_bound_dual = .0;

    /** compute lower bound with the dual variables */
    for (i = 0; i < pd->njobs + 1; i++) {
        pd->LP_lower_bound_dual += (double) pd->pi[i] * pd->rhs[i];
    }

    /** Get the LP lower bound and compute the lower bound of WCT */
    val = wctlp_objval(pd->LP, &(pd->LP_lower_bound));
    CCcheck_val_2(val, "wctlp_objval failed");
    pd->lower_bound = ((int) ceil(pd->LP_lower_bound_dual) < (int) ceil(pd->LP_lower_bound)) ? (int) ceil(pd->LP_lower_bound_dual) : (int) ceil(pd->LP_lower_bound) ;

    if (parms->stab_technique == stab_wentgnes || parms->stab_technique == stab_dynamic) {
        pd->lower_bound = (int) ceil(pd->eta_in);
    }

    if (dbg_lvl() > 0) {
        printf("Current primal LP objective: %19.16f  (LP_dual-bound %19.16f, lowerbound = %d).\n", pd->LP_lower_bound, pd->LP_lower_bound_dual, pd->lower_bound);
    }

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
int double2int(int *kpc_pi, int *scalef, const double *pi, int njobs)
{
    int    i;
    double max_dbl_nweight = -DBL_MAX;
    double max_prec_dbl = exp2(DBL_MANT_DIG - 1);
    static const double max_mwiswt   = (double) INT_MAX;
    double dbl_scalef = CC_MIN(max_prec_dbl, max_mwiswt);
    dbl_scalef /= (double) njobs;

    for (i = 0; i < njobs; ++i) {
        max_dbl_nweight =
            CC_MAX(max_dbl_nweight, pi[i]);
    }

    dbl_scalef /= CC_MAX(1.0, max_dbl_nweight);
    dbl_scalef  = floor(dbl_scalef);
    *scalef  = (int) dbl_scalef;

    for (i = 0; i < njobs; ++i) {
        double weight = pi[i] * dbl_scalef;
        assert(weight < (double) INT_MAX);
        kpc_pi[i] = (int) weight;
    }

    return 0;
}

/** Branch and Price Algorithm */

static int test_theorem_ahv(wctdata *pd, GList **branchjobs, int **min_completiontime)
{
    int val = 0;
    int i, j;
    int *temp_completiontime = (int *) NULL;
    GList *temp_branchjobs = (GList *) NULL;
    temp_completiontime = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(temp_completiontime, "Failed to allocate memory");
    fill_int(temp_completiontime, pd->njobs, INT_MAX);

    for (j = 0; j < pd->njobs; ++j) {
        int *C = (int *) NULL;

        for (i = 0; i < pd->ccount; ++i) {
            C = (int *) g_hash_table_lookup(pd->cclasses[i].table, GINT_TO_POINTER(j));

            if (pd->x[i] <= 0.0  || !C) {
                continue;
            }

            if (*C < temp_completiontime[j]) {
                temp_completiontime[j] = *C;
            }
        }
    }

    for (j = 0; j < pd->njobs; ++j) {
        int *C = (int *) NULL;
        int found = 0;

        for (i = 0; i < pd->ccount; ++i) {
            C = (int *) g_hash_table_lookup(pd->cclasses[i].table, GINT_TO_POINTER(j));

            if (pd->x[i] <= 0.0  || !C) {
                continue;
            }

            if (*C != temp_completiontime[j] && !found) {
                temp_branchjobs = g_list_append(temp_branchjobs, GINT_TO_POINTER(j));
                found = 1;
            }
        }
    }

CLEAN:

    if (val) {
        CC_IFFREE(temp_completiontime, int);
        g_list_free(*branchjobs);
    }

    *min_completiontime = temp_completiontime;
    *branchjobs = temp_branchjobs;
    return val;
}

static int grab_integral_solution_ahv(wctdata *pd, int *completion_time)
{
    int val = 0;
    int i;
    int totweight = 0;
    pmcheap *heap_sol = (pmcheap *) NULL;
    pmcheap *heap_jobs = (pmcheap *) NULL;
    solution *sol = (solution *) NULL;
    partlist *part = (partlist *) NULL;
    Job *job = (Job *) NULL;
    sol = CC_SAFE_MALLOC(1, solution);
    CCcheck_NULL_2(sol, "Failed to allocate memory")
    solution_init(sol);
    solution_alloc(sol, pd->nmachines, pd->njobs);
    val = pmcheap_init(&heap_sol, pd->nmachines);
    CCcheck_val_2(val, "Failed at pmcheap_init");
    val = pmcheap_init(&(heap_jobs), pd->njobs);
    CCcheck_val_2(val, "Failed at pmcheap_init");

    for (i = 0; i < pd->nmachines; i++) {
        pmcheap_insert(heap_sol, sol->part[i].completiontime, sol->part + i);
    }

    for (i = 0; i < pd->njobs; i++) {
        completion_time[i] -= pd->jobarray[i].processingime;
        pmcheap_insert(heap_jobs, completion_time[i], pd->jobarray + i);
    }

    while ((job = (Job *) pmcheap_min(heap_jobs))) {
        part = (partlist *) pmcheap_min(heap_sol);
        partlist_insert(part, sol->vlist, job);
        pmcheap_insert(heap_sol, part->completiontime, part);
        totweight += job->weight * part->completiontime;
    }

    Schedulesets_free(&(pd->bestcolors), &(pd->nbbest));
    pd->bestcolors = (Scheduleset *) realloc(pd->bestcolors,
                     pd->nmachines * sizeof(Scheduleset));
    CCcheck_NULL_2(pd->bestcolors, "Failed to realloc pd->bestcolors");
    /** Construct solution with the help of theorem of AHV */
    val = partlist_to_Scheduleset(sol->part, sol->nmachines, pd->njobs, &(pd->bestcolors), &(pd->nbbest));
    CCcheck_val_2(val, "Failed in conversion");
    val = Scheduleset_check(pd->bestcolors, pd->nbbest, pd->njobs);
    CCcheck_val_2(val, "Failed in Check scheduleset");
    /** Print Solution */
    printf("Intermediate solution:\n");
    print_schedule(pd->bestcolors, pd->nbbest);
    printf("with total weighted completion time = %d\n", totweight);

    if (totweight < pd->upper_bound) {
        pd->upper_bound = totweight;
        pd->besttotwct = totweight;
    }

    if (pd->upper_bound == pd->lower_bound) {
        pd->status = finished;
    }

CLEAN:
    pmcheap_free(heap_sol);
    pmcheap_free(heap_jobs);
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

static int collect_releasetime_child_ahv(wctdata *cd)
{
    int val = 0;
    int c;

    for (c = 0; c < cd->nreleasetime; ++c) {
        if (cd->releasetime_child[c].nbbest && (!cd->nbbest || cd->releasetime_child[c].upper_bound < cd->upper_bound)) {
            if (cd->nbbest) {
                Schedulesets_free(&(cd->bestcolors), &(cd->nbbest));
            }

            cd->upper_bound = cd->besttotwct = cd->releasetime_child[c].upper_bound;
            cd->nbbest = cd->releasetime_child[c].nbbest;
            cd->releasetime_child[c].nbbest = 0;
            cd->bestcolors = cd->releasetime_child[c].bestcolors;
            cd->releasetime_child[c].bestcolors = (Scheduleset *) NULL;
            /** Check if the solution is feasible, i.e. every job is covered */
        }
    }

    return val;
}

static int collect_duetime_child_ahv(wctdata *cd)
{
    int rval = 0;
    int c;

    for (c = 0; c < cd->nduetime; ++c) {
        if (cd->duetime_child[c].nbbest && (!cd->nbbest || cd->duetime_child[c].upper_bound < cd->upper_bound)) {
            if (cd->nbbest) {
                Schedulesets_free(&(cd->bestcolors), &(cd->nbbest));
            }

            cd->upper_bound = cd->duetime_child[c].upper_bound;
            cd->nbbest = cd->duetime_child[c].nbbest;
            cd->duetime_child[c].nbbest = 0;
            cd->bestcolors = cd->duetime_child[c].bestcolors;
            cd->duetime_child[c].bestcolors = (Scheduleset *) NULL;
            /** Check if the solution is feasible, i.e. every job is covered */
        }
    }

    return rval;
}

static int collect_same_child_wide(wctdata *cd)
{
    int rval = 0;
    int c;

    for (c = 0; c < cd->nsame; ++c) {
        if (cd->same_children_wide[c]->nbbest && (!cd->nbbest || cd->same_children_wide[c]->upper_bound < cd->upper_bound)) {
            if (cd->nbbest) {
                Schedulesets_free(&(cd->bestcolors), &(cd->nbbest));
            }

            cd->upper_bound =  cd->same_children_wide[c]->upper_bound;
            cd->nbbest = cd->same_children_wide[c]->nbbest;
            cd->same_children_wide[c]->nbbest = 0;
            cd->bestcolors = cd->same_children_wide[c]->bestcolors;
            cd->same_children_wide[c]->bestcolors = (Scheduleset *) NULL;
            /** Check if the solution is feasible, i.e. every job is covered */
        }
    }

    return rval;
}

static int collect_same_child_conflict(wctdata *cd)
{
    int rval = 0;
    int c;

    for (c = 0; c < cd->nsame; ++c) {
        if (cd->same_children[c].nbbest && (!cd->nbbest || cd->same_children[c].upper_bound < cd->upper_bound)) {
            if (cd->nbbest) {
                Schedulesets_free(&(cd->bestcolors), &(cd->nbbest));
            }

            cd->upper_bound = cd->same_children[c].upper_bound;
            cd->nbbest = cd->same_children[c].nbbest;
            cd->same_children[c].nbbest = 0;
            cd->bestcolors = cd->same_children[c].bestcolors;
            cd->same_children[c].bestcolors = (Scheduleset *) NULL;
            /** Check if the solution is feasible, i.e. every job is covered */
        }
    }

    return rval;
}

static int collect_diff_child_wide(wctdata *cd)
{
    int rval = 0;
    int c;

    for (c = 0; c < cd->ndiff; ++c) {
        if (cd->diff_children_wide[c]->nbbest && (!cd->nbbest || cd->diff_children_wide[c]->upper_bound < cd->upper_bound)) {
            if (cd->nbbest) {
                Schedulesets_free(&(cd->bestcolors), &(cd->nbbest));
            }

            cd->upper_bound = cd->besttotwct = cd->diff_children_wide[c]->besttotwct;
            cd->nbbest = cd->diff_children_wide[c]->nbbest;
            cd->diff_children_wide[c]->nbbest = 0;
            cd->bestcolors = cd->diff_children_wide[c]->bestcolors;
            cd->diff_children_wide[c]->bestcolors = (Scheduleset *) NULL;
            /** Check if the solution is feasible, i.e. every job is covered */
        }
    }

    return rval;
}

static int collect_diff_child_conflict(wctdata *cd)
{
    int rval = 0;
    int c;

    for (c = 0; c < cd->ndiff; ++c) {
        if (cd->diff_children[c].nbbest && (!cd->nbbest || cd->diff_children[c].upper_bound < cd->upper_bound)) {
            if (cd->nbbest) {
                Schedulesets_free(&(cd->bestcolors), &(cd->nbbest));
            }

            cd->upper_bound = cd->diff_children[c].upper_bound;
            cd->nbbest = cd->diff_children[c].nbbest;
            cd->diff_children[c].nbbest = 0;
            cd->bestcolors = cd->diff_children[c].bestcolors;
            cd->diff_children[c].bestcolors = (Scheduleset *) NULL;
            /** Check if the solution is feasible, i.e. every job is covered */
        }
    }

    return rval;
}

static int remove_finished_subtree_wide(wctdata *child)
{
    int val = 0;
    int i;
    wctdata *cd = (wctdata *) child;
    int all_same_finished = 1;
    int all_diff_finished = 1;

    while (cd) {
        for (i = 0; i < cd->nsame; ++i) {
            if (cd->same_children_wide[i]->status < finished) {
                all_same_finished = 0;
                break;
            }
        }

        if (cd->nsame && all_same_finished) {
            val = collect_same_child_wide(cd);
            CCcheck_val_2(val, "Failed in collect_same_children");

            for (i = 0;  i < cd->nsame; ++i) {
                wctdata_free(cd->same_children_wide[i]);
                CC_IFFREE(cd->same_children_wide[i], wctdata);
            }

            free(cd->same_children_wide);
            cd->same_children_wide = (wctdata **) NULL;
            cd->nsame = 0;
        }

        for (i = 0; i < cd->ndiff; ++i) {
            if (cd->diff_children_wide[i]->status < finished) {
                all_diff_finished = 0;
                break;
            }
        }

        if (cd->ndiff && all_diff_finished) {
            val = collect_diff_child_wide(cd);
            CCcheck_val_2(val, "Failed in collect_diff_children");

            for (i = 0;  i < cd->ndiff; ++i) {
                wctdata_free(cd->diff_children_wide[i]);
                CC_IFFREE(cd->diff_children_wide[i], wctdata);
            }

            CC_IFFREE(cd->diff_children_wide, wctdata *);
            cd->ndiff = 0;
        }

        if (!cd->same_children && !cd->diff_children) {
            cd->status = finished;
            CCcheck_val_2(val, "Failed to write_wctdata");
            cd = cd->parent;
        } else {
            cd = (wctdata *) NULL;
        }
    }

CLEAN:
    return val;
}


static int remove_finished_subtree_conflict(wctdata *child)
{
    int val = 0;
    int i;
    wctdata *cd = (wctdata *) child;
    int all_same_finished = 1;
    int all_diff_finished = 1;

    while (cd) {
        for (i = 0; i < cd->nsame; ++i) {
            if (cd->same_children[i].status < finished) {
                all_same_finished = 0;
                break;
            }
        }

        if (cd->nsame && all_same_finished) {
            val = collect_same_child_conflict(cd);
            CCcheck_val_2(val, "Failed in collect_same_children");

            for (i = 0;  i < cd->nsame; ++i) {
                wctdata_free(cd->same_children + i);
            }

            free(cd->same_children);
            cd->same_children = (wctdata *) NULL;
            cd->nsame = 0;
        }

        for (i = 0; i < cd->ndiff; ++i) {
            if (cd->diff_children[i].status < finished) {
                all_diff_finished = 0;
                break;
            }
        }

        if (cd->ndiff && all_diff_finished) {
            val = collect_diff_child_conflict(cd);
            CCcheck_val_2(val, "Failed in collect_diff_children");

            for (i = 0;  i < cd->ndiff; ++i) {
                wctdata_free(cd->diff_children + i);
            }

            free(cd->diff_children);
            cd->diff_children = (wctdata *) NULL;
            cd->ndiff = 0;
        }

        if (!cd->same_children && !cd->diff_children) {
            cd->status = finished;
            CCcheck_val_2(val, "Failed to write_wctdata");
            cd = cd->parent;
        } else {
            cd = (wctdata *) NULL;
        }
    }

CLEAN:
    return val;
}

MAYBE_UNUSED
static int remove_finished_subtree_ahv(wctdata *child)
{
    int rval = 0;
    int i;
    wctdata *cd = (wctdata *) child;
    int all_due_finished = 1;
    int all_ready_finished = 1;

    while (cd) {
        for (i = 0; i < cd->nduetime; ++i) {
            if (cd->duetime_child[i].status < infeasible) {
                all_due_finished = 0;
                break;
            }
        }

        if (cd->nduetime && all_due_finished) {
            rval = collect_duetime_child_ahv(cd);
            CCcheck_val_2(rval, "Failed in collect_same_children");

            for (i = 0;  i < cd->nduetime; ++i) {
                wctdata_free(cd->duetime_child + i);
            }

            free(cd->duetime_child);
            cd->duetime_child = (wctdata *) NULL;
            cd->nduetime = 0;
        }

        for (i = 0; i < cd->nreleasetime; ++i) {
            if (cd->releasetime_child[i].status < infeasible) {
                all_ready_finished = 0;
                break;
            }
        }

        if (cd->nreleasetime && all_ready_finished) {
            rval = collect_releasetime_child_ahv(cd);
            CCcheck_val_2(rval, "Failed in collect_diff_children");

            for (i = 0;  i < cd->nreleasetime; ++i) {
                wctdata_free(cd->releasetime_child + i);
            }

            free(cd->releasetime_child);
            cd->releasetime_child = (wctdata *) NULL;
            cd->nreleasetime = 0;
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

static int create_duetime_ahv(wctproblem *problem, wctdata *parent_cd, wctdata **child, int branch_job, int completiontime)
{
    int val = 0;
    wctparms *parms = &(problem->parms);
    wctdata    *pd = (wctdata *) NULL;
    pd = (wctdata *) CC_SAFE_MALLOC(1, wctdata);
    CCcheck_NULL_2(pd, "Failed to allocate pd");
    wctdata_init(pd);
    /** Init B&B data */
    pd->parent = parent_cd;
    pd->depth = parent_cd->depth + 1;
    // parent_cd->nduetime += 1;
    // parent_cd->duetime_child = pd;
    pd->branch_job = branch_job;
    pd->completiontime = completiontime;
    /** init jobs data */
    pd->njobs = parent_cd->njobs;
    pd->nmachines = parent_cd->nmachines;
    pd->duration = parent_cd->duration;
    pd->weights = parent_cd->weights;
    pd->jobarray = parent_cd->jobarray;
    pd->H_min = parent_cd->H_min;
    pd->H_max = parent_cd->H_max;
    /** init lower and upper bound data */
    pd->upper_bound = parent_cd->upper_bound;
    pd->lower_bound = parent_cd->lower_bound;
    pd->LP_lower_bound = parent_cd->LP_lower_bound;
    pd->LP_lower_bound_dual = parent_cd->LP_lower_bound_dual;
    pd->dbl_safe_lower_bound = parent_cd->dbl_safe_lower_bound;
    /** adjusted duetime */
    pd->releasetime = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(pd->releasetime, "Failed to allocate memory");
    pd->duetime = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(pd->duetime, "Failed to allocate memory");
    memcpy(pd->releasetime, parent_cd->releasetime, sizeof(int)*pd->njobs);
    memcpy(pd->duetime, parent_cd->duetime, sizeof(int)*pd->njobs);
    pd->duetime[branch_job] = completiontime;

    // for (GList *it = pd->jobarray[branch_job].due_list; it && (counter < len - pd->nmachines); it = g_list_next(it)) {
    //     int j = ((Job *)it->data)->job;
    //     pd->duetime[j] = CC_MIN(pd->duetime[j], pd->duetime[branch_job] - pd->duration[branch_job] + pd->duration[j]);
    //     counter++;
    // }

    /** Construct PricerSolver */
    if (pd->parent) {
        CCutil_start_timer(&(problem->tot_build_dd));
        pd->solver = newSolver(pd->duration, pd->weights, pd->releasetime, pd->duetime, pd->njobs, pd->H_min, pd->H_max);

        switch (parms->solver) {
            case bdd_solver:
                if ((size_t)pd->njobs != get_numberrows_bdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                } else {
                    init_tables(pd->solver);
                }

                break;

            case zdd_solver:
                if ((size_t)pd->njobs != get_numberrows_zdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                } else {
                    init_tables(pd->solver);
                }

                break;

            case DP_solver:
                break;
        }

        CCutil_suspend_timer(&(problem->tot_build_dd));
    }

    /** transfer feasible schedules */
    pd->gallocated = parent_cd->ccount;
    pd->cclasses = CC_SAFE_MALLOC(pd->gallocated, Scheduleset);
    CCcheck_NULL_2(pd->cclasses, "Failed to allocate memory");
    pd->ccount = 0;

    for (int i = 0; i < parent_cd->ccount; i++) {
        gboolean construct = 1;
        int *C = (int *) NULL;

        if ((C = (int *) g_hash_table_lookup(parent_cd->cclasses[i].table, GINT_TO_POINTER(branch_job)))) {
            if (pd->duetime[branch_job] < *C) {
                construct = 0;
            }
        }

        // for (GList *it = pd->jobarray[branch_job].due_list; it && construct && (counter < len - pd->nmachines); it = g_list_next(it)) {
        //     int job = ((Job *)it->data)->job;
        //     counter ++;

        //     if ((C = (int *) g_hash_table_lookup(parent_cd->cclasses[i].table, GINT_TO_POINTER(job)))) {
        //         if (pd->duetime[job] < *C) {
        //             construct = 0;
        //         }
        //     }
        // }

        if (construct) {
            int j;
            Scheduleset_init(pd->cclasses + pd->ccount);
            pd->cclasses[pd->ccount].members =  CC_SAFE_MALLOC(parent_cd->cclasses[i].count + 1, int);
            pd->cclasses[pd->ccount].C = CC_SAFE_MALLOC(parent_cd->cclasses[i].count, int);
            pd->cclasses[pd->ccount].table = g_hash_table_new(g_direct_hash, g_direct_equal);
            pd->cclasses[pd->ccount].count = 0;
            pd->cclasses[pd->ccount].age = 0;
            pd->cclasses[pd->ccount].totweight = 0;
            pd->cclasses[pd->ccount].totwct = 0;

            for (j = 0; j < parent_cd->cclasses[i].count && construct; ++j) {
                pd->cclasses[pd->ccount].members[(pd->cclasses[pd->ccount].count)] = parent_cd->cclasses[i].members[j];
                pd->cclasses[pd->ccount].totweight += pd->duration[parent_cd->cclasses[i].members[j]];
                pd->cclasses[pd->ccount].C[(pd->cclasses[pd->ccount].count)] = pd->cclasses[pd->ccount].totweight;
                g_hash_table_insert(pd->cclasses[pd->ccount].table, GINT_TO_POINTER(pd->cclasses[pd->ccount].members[(pd->cclasses[pd->ccount].count)]), pd->cclasses[pd->ccount].C + (pd->cclasses[pd->ccount].count));
                pd->cclasses[pd->ccount].totwct += pd->weights[parent_cd->cclasses[i].members[j]] * pd->cclasses[pd->ccount].totweight;
                pd->cclasses[pd->ccount].count++;
            }

            pd->cclasses[pd->ccount].members[pd->cclasses[pd->ccount].count] = pd->njobs;
            pd->ccount++;

            if (dbg_lvl() > 1) {
                printf("PARENT SET SAME ");

                for (j = 0; j < parent_cd->cclasses[i].count; ++j) {
                    printf(" %d", parent_cd->cclasses[i].members[j]);
                }

                printf("\n");
                printf("TRANS SET SAME");

                for (j = 0; j < pd->cclasses[pd->ccount - 1].count; ++j) {
                    printf(" %d", pd->cclasses[pd->ccount - 1].members[j]);
                }

                printf("\n");
            }
        }
    }

    for (int i = pd->ccount; i < pd->gallocated; i++) {
        Scheduleset_init(pd->cclasses + i);
    }

    val = prune_duplicated_sets(pd);
    CCcheck_val_2(val, "Failed in prune_duplicated_sets");
    *child = pd;
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

MAYBE_UNUSED
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

static int is_diff_child(wctdata *pd)
{
    int i;

    for (i = 0; pd->parent && i < pd->parent->ndiff; ++i) {
        if (pd == pd->parent->diff_children + i) {
            return 1;
        }
    }

    return 0;
}

static int is_same_child(wctdata *pd)
{
    int i;

    for (i = 0; pd->parent && i < pd->parent->nsame; ++i) {
        if (pd == pd->parent->same_children + i) {
            return 1;
        }
    }

    return 0;
}

static int create_differ_wide(wctproblem *problem, wctdata *parent_pd,
                              int v1, int v2)
{
    int val = 0;
    int i;
    wctdata    *pd =  CC_SAFE_MALLOC(1, wctdata);
    wctparms *parms = &(problem->parms);
    CCcheck_NULL_2(pd, "Failed to allocate pd");
    wctdata_init(pd);
    /** Init B&B data */
    pd->parent = parent_pd;
    pd->depth = parent_pd->depth + 1;
    parent_pd->diff_children_wide[parent_pd->ndiff++] = pd;
    pd->v1 = v1;
    pd->v2 = v2;
    /** Init jobs data */
    pd->njobs = parent_pd->njobs;
    pd->nmachines = parent_pd->nmachines;
    pd->duration = parent_pd->duration;
    pd->releasetime = parent_pd->releasetime;
    pd->duetime = parent_pd->duetime;
    pd->weights = parent_pd->weights;
    pd->jobarray = parent_pd->jobarray;
    /** Init lower bound and upper bound of node */
    pd->upper_bound = parent_pd->upper_bound;
    pd->lower_bound = parent_pd->lower_bound;
    pd->LP_lower_bound = parent_pd->LP_lower_bound;
    pd->LP_lower_bound_dual = parent_pd->LP_lower_bound_dual;
    pd->dbl_safe_lower_bound = parent_pd->dbl_safe_lower_bound;
    /* Create  graph with extra edge (v1,v2) */
    pd->ecount_differ = parent_pd->ecount_differ + 1;
    pd->ecount_same = parent_pd->ecount_same;
    // pd->elist_differ  = CC_SAFE_MALLOC(2 * pd->ecount_differ, int);
    // CCcheck_NULL_2(pd->elist_differ, "Failed to allocate pd->elist");
    // if (parent_pd->ecount_differ > 0) {
    //     memcpy(pd->elist_differ, parent_pd->elist_differ, 2 * parent_pd->ecount_differ * sizeof(int));
    // }
    // pd->elist_differ[ 2 * (pd->ecount_differ - 1)] = v1;
    // pd->elist_differ[ 2 * (pd->ecount_differ - 1) + 1] = v2;
    // pd->debugcolors = parent_pd->debugcolors;
    // pd->ndebugcolors = parent_pd->ndebugcolors;
    // /** Copy same list */
    // if (parent_pd->ecount_same > 0) {
    //     pd->elist_same = CC_SAFE_MALLOC(2 * pd->ecount_same, int);
    //     memcpy(pd->elist_same, parent_pd->elist_same, 2 * pd->ecount_same * sizeof(int));
    // }
    // /* END: Create  graph with extra edge (v1,v2) */
    // if (dbg_lvl() > 1) {
    //     printf("create_differ created following graph:\n");
    //     adjGraph_print(pd->ecount_differ, pd->elist_differ);
    // }

    /* Construction of solver*/
    if (pd->parent && (parms->solver == bdd_solver || parms->solver == zdd_solver)) {
        CCutil_start_resume_time(&(problem->tot_build_dd));
        pd->solver = copySolver(pd->parent->solver);
        add_one_conflict(pd->solver, parms, pd->v1, pd->v2, 0);

        switch (parms->solver) {
            case bdd_solver:
                if ((size_t)pd->njobs != get_numberrows_bdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                }

                break;

            case zdd_solver:
                if ((size_t)pd->njobs != get_numberrows_zdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                }

                break;

            case DP_solver:
                break;
        }

        init_tables(pd->solver);
        CCutil_suspend_timer(&(problem->tot_build_dd));
    }

    /* Transfer independent sets by removing v2 if both v1 and v2 are currently contained: */
    pd->gallocated  =  parent_pd->ccount;
    pd->cclasses =  CC_SAFE_MALLOC(pd->gallocated, Scheduleset);
    pd->ccount = 0;

    for (i = 0; i < parent_pd->ccount; ++i) {
        int j;
        int v1_found = 0;
        int construct = 1;

        for (j = 0; j < parent_pd->cclasses[i].count; ++j) {
            int current_elm = parent_pd->cclasses[i].members[j];

            if (current_elm ==  v1) {
                v1_found = 1;
            }

            if (current_elm ==  v2) {
                if (v1_found) {
                    construct = 0;
                }
            }
        }

        if (construct) {
            Scheduleset_init(pd->cclasses + pd->ccount);
            pd->cclasses[pd->ccount].members = CC_SAFE_MALLOC(parent_pd->cclasses[i].count + 1, int);
            CCcheck_NULL_2(pd->cclasses[pd->ccount].members, "Failed to allocate pd->cclasses[i].members");
            pd->cclasses[pd->ccount].C = CC_SAFE_MALLOC(parent_pd->cclasses[i].count, int);
            CCcheck_NULL_2(pd->cclasses[pd->ccount].members, "Failed to allocate memory");
            pd->cclasses[pd->ccount].table = g_hash_table_new(g_direct_hash, g_direct_equal);

            for (j = 0;  j < parent_pd->cclasses[i].count; j++) {
                pd->cclasses[pd->ccount].members[pd->cclasses[pd->ccount].count] = parent_pd->cclasses[i].members[j];
                pd->cclasses[pd->ccount].totweight += parent_pd->duration[parent_pd->cclasses[i].members[j]];
                pd->cclasses[pd->ccount].C[pd->cclasses[pd->ccount].count] = pd->cclasses[pd->ccount].totweight;
                g_hash_table_insert(pd->cclasses[pd->ccount].table, GINT_TO_POINTER(pd->cclasses[pd->ccount].members[pd->cclasses[pd->ccount].count]), pd->cclasses[pd->ccount].C + pd->cclasses[pd->ccount].count);
                pd->cclasses[pd->ccount].totwct += parent_pd->weights[parent_pd->cclasses[i].members[j]] * pd->cclasses[pd->ccount].totweight;
                (pd->cclasses[pd->ccount].count)++;
            }

            pd->cclasses[pd->ccount].members[pd->cclasses[pd->ccount].count] = pd->njobs;
            pd->ccount++;
        }

        if (dbg_lvl() > 1 && construct) {
            printf("PARENT SET DIFFER");

            for (j = 0; j < parent_pd->cclasses[i].count; ++j) {
                printf(" %d", parent_pd->cclasses[i].members[j]);
            }

            printf("\n");
            printf("TRANS SET DIFFER");

            for (j = 0; j < pd->cclasses[pd->ccount - 1].count; ++j) {
                printf(" %d", pd->cclasses[pd->ccount - 1].members[j]);
            }

            printf("\n");
        }

        CCcheck_val_2(val, "Illegal colorset created in create_differ\n!");
    }

    for (i = pd->ccount; i < pd->gallocated; i++) {
        Scheduleset_init(pd->cclasses + i);
    }

    val = prune_duplicated_sets(pd);
    CCcheck_val_2(val, "Failed in prune_duplicated_sets");
CLEAN:

    if (val) {
        if (pd) {
            wctdata_free(pd);
            free(pd);
        }

        parent_pd->diff_children = (wctdata *) NULL;
    }

    return val;
}


static int create_differ_conflict(wctproblem *problem, wctdata *parent_pd, wctdata **child, int v1, int v2)
{
    int val = 0;
    int i;
    wctdata    *pd =  CC_SAFE_MALLOC(1, wctdata);
    wctparms *parms = &(problem->parms);
    CCcheck_NULL_2(pd, "Failed to allocate pd");
    wctdata_init(pd);
    /** Init B&B data */
    pd->parent = parent_pd;
    pd->depth = parent_pd->depth + 1;
    // parent_pd->ndiff         += 1;
    // parent_pd->diff_children = pd;
    pd->v1 = v1;
    pd->v2 = v2;
    /** Init jobs data */
    pd->njobs = parent_pd->njobs;
    pd->nmachines = parent_pd->nmachines;
    pd->duration = parent_pd->duration;
    pd->weights = parent_pd->weights;
    pd->jobarray = parent_pd->jobarray;
    pd->releasetime = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(pd->releasetime, "Failed to allocate  memory");
    pd->duetime = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(pd->duetime, "Failed to allocate memory");
    memcpy(pd->releasetime, parent_pd->releasetime, sizeof(int)*pd->njobs);
    memcpy(pd->duetime, parent_pd->duetime, sizeof(int)*pd->njobs);
    /** Init lower bound and upper bound of node */
    pd->upper_bound = parent_pd->upper_bound;
    pd->lower_bound = parent_pd->lower_bound;
    pd->LP_lower_bound = parent_pd->LP_lower_bound;
    pd->LP_lower_bound_dual = parent_pd->LP_lower_bound_dual;
    pd->dbl_safe_lower_bound = parent_pd->dbl_safe_lower_bound;
    /* Create  graph with extra edge (v1,v2) */
    pd->ecount_differ = parent_pd->ecount_differ + 1;
    pd->ecount_same = parent_pd->ecount_same;

    /* Construction of solver*/
    if (pd->parent) {
        CCutil_start_resume_time(&(problem->tot_build_dd));
        pd->solver = copySolver(pd->parent->solver);
        add_one_conflict(pd->solver, parms, pd->v1, pd->v2, 0);
        set_release_due_time(pd->solver, pd->releasetime, pd->duetime);

        switch (parms->solver) {
            case bdd_solver:
                if ((size_t)pd->njobs != get_numberrows_bdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                }

                break;

            case zdd_solver:
                if ((size_t)pd->njobs != get_numberrows_zdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                }

                break;

            case DP_solver:
                break;
        }

        init_tables(pd->solver);
        CCutil_suspend_timer(&(problem->tot_build_dd));
    }

    /* Transfer independent sets by removing v2 if both v1 and v2 are currently contained: */
    pd->gallocated  =  parent_pd->ccount;
    pd->cclasses =  CC_SAFE_MALLOC(pd->gallocated, Scheduleset);
    pd->ccount = 0;

    for (i = 0; i < parent_pd->ccount; ++i) {
        int j;
        int construct = 1;
        gboolean v1_in;
        gboolean v2_in;
        v1_in = g_hash_table_contains(parent_pd->cclasses[i].table, GINT_TO_POINTER(v1));
        v2_in = g_hash_table_contains(parent_pd->cclasses[i].table, GINT_TO_POINTER(v2));
        construct = (v1_in && v2_in) ? 0 : 1;

        if (construct) {
            Scheduleset_init(pd->cclasses + pd->ccount);
            pd->cclasses[pd->ccount].members = CC_SAFE_MALLOC(parent_pd->cclasses[i].count + 1, int);
            CCcheck_NULL_2(pd->cclasses[pd->ccount].members, "Failed to allocate pd->cclasses[i].members");
            pd->cclasses[pd->ccount].C = CC_SAFE_MALLOC(parent_pd->cclasses[i].count, int);
            CCcheck_NULL_2(pd->cclasses[pd->ccount].members, "Failed to allocate memory");
            pd->cclasses[pd->ccount].table = g_hash_table_new(g_direct_hash, g_direct_equal);

            for (j = 0;  j < parent_pd->cclasses[i].count; j++) {
                pd->cclasses[pd->ccount].members[pd->cclasses[pd->ccount].count] = parent_pd->cclasses[i].members[j];
                pd->cclasses[pd->ccount].totweight += parent_pd->duration[parent_pd->cclasses[i].members[j]];
                pd->cclasses[pd->ccount].C[pd->cclasses[pd->ccount].count] = pd->cclasses[pd->ccount].totweight;
                g_hash_table_insert(pd->cclasses[pd->ccount].table, GINT_TO_POINTER(pd->cclasses[pd->ccount].members[pd->cclasses[pd->ccount].count]), pd->cclasses[pd->ccount].C + pd->cclasses[pd->ccount].count);
                pd->cclasses[pd->ccount].totwct += parent_pd->weights[parent_pd->cclasses[i].members[j]] * pd->cclasses[pd->ccount].totweight;
                (pd->cclasses[pd->ccount].count)++;
            }

            pd->cclasses[pd->ccount].members[pd->cclasses[pd->ccount].count] = pd->njobs;
            pd->ccount++;
        }

        if (dbg_lvl() > 1 && construct) {
            printf("PARENT SET DIFFER");

            for (j = 0; j < parent_pd->cclasses[i].count; ++j) {
                printf(" %d", parent_pd->cclasses[i].members[j]);
            }

            printf("\n");
            printf("TRANS SET DIFFER");

            for (j = 0; j < pd->cclasses[pd->ccount - 1].count; ++j) {
                printf(" %d", pd->cclasses[pd->ccount - 1].members[j]);
            }

            printf("\n");
        }

        CCcheck_val_2(val, "Illegal colorset created in create_differ\n!");
    }

    for (i = pd->ccount; i < pd->gallocated; i++) {
        Scheduleset_init(pd->cclasses + i);
    }

    val = prune_duplicated_sets(pd);
    CCcheck_val_2(val, "Failed in prune_duplicated_sets");
    *child = pd;
CLEAN:

    if (val) {
        if (pd) {
            wctdata_free(pd);
            free(pd);
        }
    }

    return val;
}

static int transfer_same_cclasses_wide(wctdata *pd,
                                       const Scheduleset *parent_cclasses,
                                       int   parent_ccount,
                                       int   *v1_wide,
                                       int   *v2_wide)
{
    int val = 0;
    int i;
    /* Transfer independent sets: */
    pd->gallocated = pd->ccount   =  parent_ccount ;
    pd->cclasses = CC_SAFE_MALLOC(pd->gallocated, Scheduleset);
    CCcheck_NULL_2(pd->cclasses, "Failed to allocate memory");
    pd->ccount = 0;

    for (i = 0; i < parent_ccount; ++i) {
        int j;
        int construct = 1;
        gboolean v1_in;
        gboolean v2_in;

        for (j = 0; j < pd->nb_wide; j++) {
            v1_in =  g_hash_table_contains(parent_cclasses[i].table, GINT_TO_POINTER(v1_wide[i]));
            v2_in =  g_hash_table_contains(parent_cclasses[i].table, GINT_TO_POINTER(v2_wide[i]));

            if ((v1_in == 1 && v2_in == 0) || (v1_in == 0 && v2_in == 1)) {
                construct = 0;
                break;
            }
        }

        if (construct) {
            Scheduleset_init(pd->cclasses + pd->ccount);
            pd->cclasses[pd->ccount].members =  CC_SAFE_MALLOC(parent_cclasses[i].count + 1, int);
            pd->cclasses[pd->ccount].C = CC_SAFE_MALLOC(parent_cclasses[i].count, int);
            pd->cclasses[pd->ccount].table = g_hash_table_new(g_direct_hash, g_direct_equal);
            pd->cclasses[pd->ccount].count = 0;
            pd->cclasses[pd->ccount].age = 0;
            pd->cclasses[pd->ccount].totweight = 0;
            pd->cclasses[pd->ccount].totwct = 0;
        } else {
            continue;
        }

        for (j = 0; j < parent_cclasses[i].count && construct; ++j) {
            pd->cclasses[pd->ccount].members[(pd->cclasses[pd->ccount].count)] = parent_cclasses[i].members[j];
            pd->cclasses[pd->ccount].totweight += pd->duration[parent_cclasses[i].members[j]];
            pd->cclasses[pd->ccount].C[(pd->cclasses[pd->ccount].count)] = pd->cclasses[pd->ccount].totweight;
            g_hash_table_insert(pd->cclasses[pd->ccount].table, GINT_TO_POINTER(pd->cclasses[pd->ccount].members[(pd->cclasses[pd->ccount].count)]), pd->cclasses[pd->ccount].C + (pd->cclasses[pd->ccount].count));
            pd->cclasses[pd->ccount].totwct += pd->weights[parent_cclasses[i].members[j]] * pd->cclasses[pd->ccount].totweight;
            pd->cclasses[pd->ccount].count++;
            /* else 'parent_cclasses[i].members[j] == v2' and we skip it*/
        }

        if (construct) {
            pd->cclasses[pd->ccount].members[pd->cclasses[pd->ccount].count] = pd->njobs;
            pd->ccount++;
        }

        if (dbg_lvl() > 1 && construct) {
            printf("PARENT SET SAME ");

            for (j = 0; j < parent_cclasses[i].count; ++j) {
                printf(" %d", parent_cclasses[i].members[j]);
            }

            printf("\n");
            printf("TRANS SET SAME");

            for (j = 0; j < pd->cclasses[pd->ccount - 1].count; ++j) {
                printf(" %d", pd->cclasses[pd->ccount - 1].members[j]);
            }

            printf("\n");
        }
    }

    for (i = pd->ccount; i < pd->gallocated; i++) {
        Scheduleset_init(pd->cclasses + i);
    }

    val = prune_duplicated_sets(pd);
    CCcheck_val_2(val, "Failed in prune_duplicated_sets");
    /* END Transfer independent sets: */
CLEAN:
    return val;
}

static int transfer_same_cclasses(wctdata *pd,
                                  const Scheduleset *parent_cclasses,
                                  int   parent_ccount,
                                  int   v1,
                                  int   v2)
{
    int val = 0;
    int i;
    /* Transfer independent sets: */
    pd->gallocated = pd->ccount   =  parent_ccount ;
    pd->cclasses = CC_SAFE_MALLOC(pd->gallocated, Scheduleset);
    CCcheck_NULL_2(pd->cclasses, "Failed to allocate memory");
    pd->ccount = 0;

    for (i = 0; i < parent_ccount; ++i) {
        int j;
        int construct = 1;
        gboolean v1_in;
        gboolean v2_in;
        v1_in = g_hash_table_contains(parent_cclasses[i].table, GINT_TO_POINTER(v1));
        v2_in = g_hash_table_contains(parent_cclasses[i].table, GINT_TO_POINTER(v2));

        if ((v1_in == 1 && v2_in == 0) || (v1_in == 0 && v2_in == 1)) {
            construct = 0;
        } else {
            Scheduleset_init(pd->cclasses + pd->ccount);
            pd->cclasses[pd->ccount].members =  CC_SAFE_MALLOC(parent_cclasses[i].count + 1, int);
            pd->cclasses[pd->ccount].C = CC_SAFE_MALLOC(parent_cclasses[i].count, int);
            pd->cclasses[pd->ccount].table = g_hash_table_new(g_direct_hash, g_direct_equal);
            pd->cclasses[pd->ccount].count = 0;
            pd->cclasses[pd->ccount].age = 0;
            pd->cclasses[pd->ccount].totweight = 0;
            pd->cclasses[pd->ccount].totwct = 0;
        }

        for (j = 0; j < parent_cclasses[i].count && construct; ++j) {
            pd->cclasses[pd->ccount].members[(pd->cclasses[pd->ccount].count)] = parent_cclasses[i].members[j];
            pd->cclasses[pd->ccount].totweight += pd->duration[parent_cclasses[i].members[j]];
            pd->cclasses[pd->ccount].C[(pd->cclasses[pd->ccount].count)] = pd->cclasses[pd->ccount].totweight;
            g_hash_table_insert(pd->cclasses[pd->ccount].table, GINT_TO_POINTER(pd->cclasses[pd->ccount].members[(pd->cclasses[pd->ccount].count)]), pd->cclasses[pd->ccount].C + (pd->cclasses[pd->ccount].count));
            pd->cclasses[pd->ccount].totwct += pd->weights[parent_cclasses[i].members[j]] * pd->cclasses[pd->ccount].totweight;
            pd->cclasses[pd->ccount].count++;
        }

        if (construct) {
            pd->cclasses[pd->ccount].members[pd->cclasses[pd->ccount].count] = pd->njobs;
            pd->ccount++;
        }

        if (dbg_lvl() > 1 && construct) {
            printf("PARENT SET SAME ");

            for (j = 0; j < parent_cclasses[i].count; ++j) {
                printf(" %d", parent_cclasses[i].members[j]);
            }

            printf("\n");
            printf("TRANS SET SAME");

            for (j = 0; j < pd->cclasses[pd->ccount - 1].count; ++j) {
                printf(" %d", pd->cclasses[pd->ccount - 1].members[j]);
            }

            printf("\n");
        }
    }

    for (i = pd->ccount; i < pd->gallocated; i++) {
        Scheduleset_init(pd->cclasses + i);
    }

    val = prune_duplicated_sets(pd);
    CCcheck_val_2(val, "Failed in prune_duplicated_sets");
    /* END Transfer independent sets: */
CLEAN:
    return val;
}

MAYBE_UNUSED
static int mark_neighborhood(int *neighbor_marker,
                             int vcount, int ecount, int elist[],
                             int v)

{
    int val = 0;
    int i;
    adjGraph G;
    adjGraph_init(&G);
    val = adjGraph_buildquick(&G, vcount, ecount, elist);
    CCcheck_val_2(val, "COLORadjgraph_build");

    for (i = 0; i < vcount; ++i) {
        neighbor_marker[i] = 0;
    }

    for (i = 0; i < G.nodelist[v].degree; ++i) {
        int v_i = G.nodelist[v].adj[i];
        neighbor_marker[v_i] = 1;
    }

    adjGraph_freequick(&G);
CLEAN:
    return val;
}


static int create_same_wide(wctproblem *problem, wctdata *parent_pd, int *v1_wide, int *v2_wide, int nb_wide)
{
    int val = 0;
    wctdata   *pd = (wctdata *) NULL;
    wctparms *parms = &(problem->parms);
    pd = (wctdata *) CC_SAFE_MALLOC(1, wctdata);
    CCcheck_NULL_2(pd, "Failed to allocate pd");
    wctdata_init(pd);
    /** Init B&B data */
    pd->depth = parent_pd->depth + 1;
    parent_pd->same_children_wide[parent_pd->nsame++] = pd;
    /** Init jobs data */
    pd->njobs = parent_pd->njobs;
    pd->nmachines = parent_pd->nmachines;
    pd->duration = parent_pd->duration;
    pd->releasetime = parent_pd->releasetime;
    pd->duetime = parent_pd->duetime;
    pd->weights = parent_pd->weights;
    pd->jobarray = parent_pd->jobarray;
    pd->ecount_same = parent_pd->ecount_same + 1;
    pd->ecount_differ = parent_pd->ecount_differ;
    /** Init lower bound and upper bound */
    pd->upper_bound = parent_pd->upper_bound;
    pd->lower_bound = parent_pd->lower_bound;
    pd->LP_lower_bound = parent_pd->LP_lower_bound;
    pd->LP_lower_bound_dual = parent_pd->LP_lower_bound_dual;
    pd->dbl_safe_lower_bound = parent_pd->dbl_safe_lower_bound;
    pd->parent = parent_pd;
    pd->debugcolors = parent_pd->debugcolors;
    pd->ndebugcolors = parent_pd->ndebugcolors;
    pd->nb_wide = nb_wide;
    pd->v1_wide = CC_SAFE_MALLOC(nb_wide, int);
    CCcheck_NULL_2(pd->v1_wide, "Failed to allocate memory");
    pd->v2_wide = CC_SAFE_MALLOC(nb_wide, int);
    CCcheck_NULL_2(pd->v2_wide, "Failed to allocate memory");
    memcpy(pd->v1_wide, v1_wide, sizeof(int)*nb_wide);
    memcpy(pd->v2_wide, v2_wide, sizeof(int)*nb_wide);
    pd->elist_same = CC_SAFE_MALLOC(2 * pd->nb_wide, int);
    CCcheck_NULL_2(pd->elist_same, "Failed to allocate memory");

    for (int i = 0; i < nb_wide; i++) {
        pd->elist_same[2 * i] = v1_wide[i];
        pd->elist_same[2 * i + 1] = v2_wide[i];
    }

    pd->ecount_same += nb_wide;

    /* Construction of solver*/
    if (pd->parent && (parms->solver == bdd_solver || parms->solver == zdd_solver)) {
        CCutil_start_resume_time(&(problem->tot_build_dd));
        pd->solver = copySolver(pd->parent->solver);
        add_conflict_constraints(pd->solver, &(problem->parms), pd->elist_same, pd->nb_wide, NULL, 0);

        switch (parms->solver) {
            case bdd_solver:
                if ((size_t)pd->njobs != get_numberrows_bdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                }

                break;

            case zdd_solver:
                if ((size_t)pd->njobs != get_numberrows_zdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                }

                break;

            case DP_solver:
                break;
        }

        init_tables(pd->solver);
        CCutil_suspend_timer(&(problem->tot_build_dd));
    }

    val = transfer_same_cclasses_wide(pd,
                                      parent_pd->cclasses,
                                      parent_pd->ccount,
                                      pd->v1_wide, pd->v2_wide);
    CCcheck_val_2(val, "Failed in transfer_same_cclasses");
CLEAN:

    if (val) {
        if (pd) {
            wctdata_free(pd);
            free(pd);
        }

        parent_pd->same_children = (wctdata *) NULL;
    }

    CC_IFFREE(pd->elist_same, int);
    pd->ecount_same = 0;
    return val;
}


static int create_same_conflict(wctproblem *problem, wctdata *parent_pd, wctdata **child, int v1, int v2)
{
    int val = 0;
    wctdata   *pd = (wctdata *) NULL;
    wctparms *parms = &(problem->parms);
    pd = (wctdata *) CC_SAFE_MALLOC(1, wctdata);
    CCcheck_NULL_2(pd, "Failed to allocate pd");
    wctdata_init(pd);
    /** Init B&B data */
    pd->parent = parent_pd;
    pd->depth = parent_pd->depth + 1;
    // parent_pd->nsame = 1;
    // parent_pd->same_children = pd;
    pd->v1 = v1;
    pd->v2 = v2;
    /** Init jobs data */
    pd->njobs = parent_pd->njobs;
    pd->nmachines = parent_pd->nmachines;
    pd->duration = parent_pd->duration;
    pd->weights = parent_pd->weights;
    pd->jobarray = parent_pd->jobarray;
    pd->ecount_same = parent_pd->ecount_same + 1;
    pd->ecount_differ = parent_pd->ecount_differ;
    pd->releasetime = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(pd->releasetime, "Failed to allocate memory");
    pd->duetime = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(pd->duetime, "Failed to allocate memory");
    memcpy(pd->releasetime, parent_pd->releasetime, sizeof(int)*pd->njobs);
    memcpy(pd->duetime, parent_pd->duetime, sizeof(int)*pd->njobs);
    /** Init lower bound and upper bound */
    pd->upper_bound = parent_pd->upper_bound;
    pd->lower_bound = parent_pd->lower_bound;
    pd->LP_lower_bound = parent_pd->LP_lower_bound;
    pd->LP_lower_bound_dual = parent_pd->LP_lower_bound_dual;
    pd->dbl_safe_lower_bound = parent_pd->dbl_safe_lower_bound;
    pd->parent = parent_pd;
    pd->debugcolors = parent_pd->debugcolors;
    pd->ndebugcolors = parent_pd->ndebugcolors;

    /* Construction of solver*/
    if (pd->parent) {
        CCutil_start_resume_time(&(problem->tot_build_dd));
        pd->solver = copySolver(pd->parent->solver);
        add_one_conflict(pd->solver, parms, pd->v1, pd->v2, 1);
        set_release_due_time(pd->solver, pd->releasetime, pd->duetime);

        switch (parms->solver) {
            case bdd_solver:
                if ((size_t)pd->njobs != get_numberrows_bdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                }

                break;

            case zdd_solver:
                if ((size_t)pd->njobs != get_numberrows_zdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                }

                break;

            case DP_solver:
                break;
        }

        init_tables(pd->solver);
        CCutil_suspend_timer(&(problem->tot_build_dd));
    }

    val = transfer_same_cclasses(pd, parent_pd->cclasses, parent_pd->ccount, v1, v2);
    CCcheck_val_2(val, "Failed in transfer_same_cclasses");
    *child = pd;
CLEAN:

    if (val) {
        if (pd) {
            wctdata_free(pd);
            free(pd);
        }
    }

    return val;
}

int recover_elist(wctdata *pd)
{
    int val = 0;
    wctdata    **path  = (wctdata **) NULL;
    int            npath = 0;
    wctdata     *tmp_cd  = pd;
    int           *elist_same   = (int *) NULL;
    int *elist_differ = (int *) NULL;
    int            ecount_differ = 0;
    int ecount_same = 0;
    int            ndiff   = 0;
    int nsame = 0;
    int            i;

    while (tmp_cd) {
        npath++;
        tmp_cd = tmp_cd->parent;
    }

    path = CC_SAFE_MALLOC(npath, wctdata *);
    CCcheck_NULL_2(path, "Failed to allocate path.");
    tmp_cd = pd;
    i      = npath;

    while (tmp_cd) {
        i--;
        path[i] = tmp_cd;

        if (is_diff_child(tmp_cd)) {
            ndiff++;
        } else if (is_same_child(tmp_cd)) {
            nsame++;
        }

        tmp_cd = tmp_cd->parent;
    }

    if (ndiff) {
        elist_differ = CC_SAFE_MALLOC(2 * ndiff, int);
        CCcheck_NULL(elist_differ, "Failed to allocate memory");
    }

    if (nsame) {
        elist_same = CC_SAFE_MALLOC(2 * nsame, int);
        CCcheck_NULL(elist_same, "Failed to allocate memory");
    }

    for (i = npath - 1; i >= 0; i--) {
        wctdata *cur_cd = path[i];

        if (is_diff_child(cur_cd)) {
            elist_differ[2 * ecount_differ] = cur_cd->v1;
            elist_differ[2 * ecount_differ + 1] = cur_cd->v2;
            ecount_differ++;
        } else  if (is_same_child(cur_cd)) {
            elist_same[2 * ecount_same] = cur_cd->v1;
            elist_same[2 * ecount_same + 1] = cur_cd->v2;
            ecount_same++;
        }
    }

    pd->elist_differ = elist_differ;
    pd->ecount_differ = ecount_differ;
    elist_differ     = (int *) NULL;
    pd->elist_same = elist_same;
    pd->ecount_same = ecount_same;
    elist_same = (int *) NULL;
CLEAN:
    CC_IFFREE(elist_differ, int);
    CC_IFFREE(elist_same, int);
    CC_IFFREE(path, wctdata *);
    return val;
}

static int create_releasetime_ahv(wctproblem *problem, wctdata *parent_cd, wctdata **child, int branch_job, int completiontime)
{
    int val = 0;
    wctparms *parms = &(problem->parms);
    wctdata    *pd = (wctdata *) NULL;
    pd = (wctdata *) CC_SAFE_MALLOC(1, wctdata);
    CCcheck_NULL_2(pd, "Failed to allocate pd");
    wctdata_init(pd);
    /** Init B&B data */
    pd->parent = parent_cd;
    pd->depth = parent_cd->depth + 1;
    // parent_cd->nreleasetime += 1;
    // parent_cd->releasetime_child = pd;
    pd->branch_job = branch_job;
    pd->completiontime = completiontime;
    /** Init jobs data */
    pd->njobs  = parent_cd->njobs;
    pd->nmachines = parent_cd->nmachines;
    pd->duration = parent_cd->duration;
    pd->weights = parent_cd->weights;
    pd->jobarray = parent_cd->jobarray;
    pd->H_min = parent_cd->H_min;
    pd->H_max = parent_cd->H_max;
    /** Init lower and upper bound at the node */
    pd->upper_bound = parent_cd->upper_bound;
    pd->lower_bound = parent_cd->lower_bound;
    pd->LP_lower_bound = parent_cd->LP_lower_bound;
    pd->LP_lower_bound_dual = parent_cd->LP_lower_bound_dual;
    pd->dbl_safe_lower_bound = parent_cd->dbl_safe_lower_bound;
    /** adjusted release_time */
    pd->releasetime = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(pd->releasetime, "Failed to allocate memory");
    pd->duetime = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(pd->duetime, "Failed to allocate memory");
    memcpy(pd->releasetime, parent_cd->releasetime, sizeof(int)*pd->njobs);
    memcpy(pd->duetime, parent_cd->duetime, sizeof(int)*pd->njobs);
    pd->releasetime[branch_job] = completiontime + 1 - pd->duration[branch_job];

    // for (GList *it = pd->jobarray[branch_job].ready_list; it; it = g_list_next(it)) {
    //     int j = ((Job *)it->data)->job;
    //     pd->releasetime[j] = CC_MAX(pd->releasetime[j], pd->releasetime[branch_job]);
    // }

    /** Construct PricerSolver */
    if (pd->parent) {
        CCutil_start_timer(&(problem->tot_build_dd));
        pd->solver = newSolver(pd->duration, pd->weights, pd->releasetime, pd->duetime, pd->njobs, pd->H_min, pd->H_max);

        switch (parms->solver) {
            case bdd_solver:
                if ((size_t)pd->njobs != get_numberrows_bdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                } else {
                    init_tables(pd->solver);
                }

                break;

            case zdd_solver:
                if ((size_t)pd->njobs != get_numberrows_zdd(pd->solver)) {
                    pd->status = infeasible;
                    CCutil_suspend_timer(&(problem->tot_build_dd));
                    goto CLEAN;
                } else {
                    init_tables(pd->solver);
                }

                break;

            case DP_solver:
                break;
        }

        CCutil_suspend_timer(&(problem->tot_build_dd));;
    }

    /** transfer feasible schedules */
    pd->gallocated = parent_cd->ccount;
    pd->cclasses = CC_SAFE_MALLOC(pd->gallocated, Scheduleset);
    CCcheck_NULL_2(pd->cclasses, "Failed to allocate memory");
    pd->ccount = 0;

    for (int i = 0; i < parent_cd->ccount; i++) {
        gboolean construct = 1;
        int *C = (int *) NULL;

        if ((C = (int *) g_hash_table_lookup(parent_cd->cclasses[i].table, GINT_TO_POINTER(branch_job)))) {
            if (pd->releasetime[branch_job] + pd->duration[branch_job] > *C) {
                construct = 0;
            }
        }

        // for (GList *it = pd->jobarray[branch_job].ready_list; it && construct; it = g_list_next(it)) {
        //     int job = ((Job *)it->data)->job;

        //     if ((C = (int *) g_hash_table_lookup(parent_cd->cclasses[i].table, GINT_TO_POINTER(job)))) {
        //         if (pd->releasetime[job] + pd->duration[job] > *C) {
        //             construct = 0;
        //         }
        //     }
        // }

        if (construct) {
            int j;
            Scheduleset_init(pd->cclasses + pd->ccount);
            pd->cclasses[pd->ccount].members =  CC_SAFE_MALLOC(parent_cd->cclasses[i].count + 1, int);
            pd->cclasses[pd->ccount].C = CC_SAFE_MALLOC(parent_cd->cclasses[i].count, int);
            pd->cclasses[pd->ccount].table = g_hash_table_new(g_direct_hash, g_direct_equal);
            pd->cclasses[pd->ccount].count = 0;
            pd->cclasses[pd->ccount].age = 0;
            pd->cclasses[pd->ccount].totweight = 0;
            pd->cclasses[pd->ccount].totwct = 0;

            for (j = 0; j < parent_cd->cclasses[i].count && construct; ++j) {
                pd->cclasses[pd->ccount].members[(pd->cclasses[pd->ccount].count)] = parent_cd->cclasses[i].members[j];
                pd->cclasses[pd->ccount].totweight += pd->duration[parent_cd->cclasses[i].members[j]];
                pd->cclasses[pd->ccount].C[(pd->cclasses[pd->ccount].count)] = pd->cclasses[pd->ccount].totweight;
                g_hash_table_insert(pd->cclasses[pd->ccount].table, GINT_TO_POINTER(pd->cclasses[pd->ccount].members[(pd->cclasses[pd->ccount].count)]), pd->cclasses[pd->ccount].C + (pd->cclasses[pd->ccount].count));
                pd->cclasses[pd->ccount].totwct += pd->weights[parent_cd->cclasses[i].members[j]] * pd->cclasses[pd->ccount].totweight;
                pd->cclasses[pd->ccount].count++;
            }

            pd->cclasses[pd->ccount].members[pd->cclasses[pd->ccount].count] = pd->njobs;
            pd->ccount++;

            if (dbg_lvl() > 1) {
                printf("PARENT SET SAME ");

                for (j = 0; j < parent_cd->cclasses[i].count; ++j) {
                    printf(" %d", parent_cd->cclasses[i].members[j]);
                }

                printf("\n");
                printf("TRANS SET SAME");

                for (j = 0; j < pd->cclasses[pd->ccount - 1].count; ++j) {
                    printf(" %d", pd->cclasses[pd->ccount - 1].members[j]);
                }

                printf("\n");
            }
        }
    }

    for (int i = pd->ccount; i < pd->gallocated; i++) {
        Scheduleset_init(pd->cclasses + i);
    }

    val = prune_duplicated_sets(pd);
    CCcheck_val_2(val, "Failed in prune_duplicated_sets");
    *child = pd;
CLEAN:

    if (val) {
        if (pd) {
            wctdata_free(pd);
            free(pd);
        }

        parent_cd->releasetime_child = (wctdata *) NULL;
    }

    return val;
}



static int find_strongest_children_conflict(int *strongest_v1,

        int *strongest_v2,
        wctdata    *pd,
        wctproblem *problem,
        pmcheap *cand_heap,
        int          *nodepair_refs,
        double       *nodepair_weights)
{
    int    val = 0;
    int    max_non_improving_branches  = 3; /* pd->njobs / 100 + 1; */
    int    remaining_branches          = max_non_improving_branches;
    double strongest_dbl_lb = 0.0;
    int   *min_nodepair;
    wctparms *parms = &(problem->parms);
    *strongest_v1 = -1;
    *strongest_v2 = -1;

    switch (parms->strong_branching) {
        case use_strong_branching:
            while ((min_nodepair = (int *) pmcheap_min(cand_heap))
                    && (remaining_branches--)) {
                int v1 = -1, v2 = -1;
                double dbl_child_lb;
                inodepair_ref_key(&v1, &v2, (int)(min_nodepair - nodepair_refs));
                assert(v1 < v2);
                wctdata *same_children = (wctdata *) NULL;
                wctdata *diff_children = (wctdata *) NULL;

                if (dbg_lvl() == 0) {
                    printf("Creating branches for v1 = %d, v2 = %d (node-pair weight %f)\n", v1,
                           v2,
                           nodepair_weights[(int)(min_nodepair - nodepair_refs)]);
                }

                /* Create DIFFER and SAME */
                val = create_same_conflict(problem, pd, &(same_children), v1, v2);
                CCcheck_val_2(val, "Failed in create_same");
                val = create_differ_conflict(problem, pd, &(diff_children), v1, v2);
                CCcheck_val_2(val, "Failed in create_differ");
                same_children->maxiterations = 5;
                diff_children->maxiterations = 5;
                compute_lower_bound(problem, same_children);
                compute_lower_bound(problem, diff_children);
                dbl_child_lb = (same_children->LP_lower_bound <
                                diff_children->LP_lower_bound) ?
                               same_children->LP_lower_bound :
                               diff_children->LP_lower_bound;

                if (dbl_child_lb > strongest_dbl_lb) {
                    strongest_dbl_lb = dbl_child_lb;
                    *strongest_v1     = v1;
                    *strongest_v2     = v2;
                    remaining_branches = max_non_improving_branches;

                    if (pd->same_children) {
                        wctdata_free(pd->same_children);
                        free(pd->same_children);
                        pd->nsame = 0;
                    }

                    same_children->maxiterations = 1000000;
                    pd->same_children = same_children;
                    pd->nsame = 1;

                    if (pd->diff_children) {
                        wctdata_free(pd->diff_children);
                        free(pd->diff_children);
                        pd->ndiff = 0;
                    }

                    diff_children->maxiterations = 1000000;
                    pd->diff_children = diff_children;
                    pd->ndiff = 1;
                } else {
                    wctdata_free(same_children);
                    free(same_children);
                    wctdata_free(diff_children);
                    free(diff_children);
                }

                if (dbg_lvl() > 1) {
                    printf("Found child bound of %f for v1 = %d, v2 = %d, nodepair_weight = %f .\n",
                           dbl_child_lb,
                           v1, v2,
                           nodepair_weights[(int)(min_nodepair - nodepair_refs)]);
                }
            }

            if (dbg_lvl() > 0) {
                int nodepair_ref = nodepair_ref_key(*strongest_v1, *strongest_v2);
                printf("Found strongest child bound of %f for v1 = %d, "
                       "v2 = %d, nodepair_weight = %f .\n",
                       strongest_dbl_lb, *strongest_v1,
                       *strongest_v2, nodepair_weights[nodepair_ref]);
            }

            break;

        case no_strong_branching:
            min_nodepair = (int *) pmcheap_min(cand_heap);
            int v1 = -1, v2 = -1;
            inodepair_ref_key(&v1, &v2, (int)(min_nodepair - nodepair_refs));
            assert(v1 < v2);

            if (dbg_lvl() == 0) {
                printf("Creating branches for v1 = %d, v2 = %d (node-pair weight %f)\n", v1,
                       v2,
                       nodepair_weights[(int)(min_nodepair - nodepair_refs)]);
            }

            val = create_same_conflict(problem, pd, &(pd->same_children), v1, v2);
            CCcheck_val_2(val, "Failed in create_same");
            pd->nsame = 1;
            val = create_differ_conflict(problem, pd, &(pd->diff_children), v1, v2);
            CCcheck_val_2(val, "Failed in create_differ");
            pd->ndiff = 1;
            break;
    }

CLEAN:
    return val;
}

static int find_strongest_children_ahv(int *strongest_v1, wctdata *pd, wctproblem *problem, GList *branchnodes, int *completiontime)
{
    int    rval = 0;
    int    max_non_improving_branches  = 3; /* cd->njobs / 100 + 1; */
    int    remaining_branches          = max_non_improving_branches;
    double strongest_dbl_lb = 0.0;
    wctparms *parms = &(problem->parms);
    *strongest_v1 = -1;
    GList *it = branchnodes;
    int v1;

    switch (parms->strong_branching) {
        case use_strong_branching:
            while (it && (remaining_branches--)) {
                v1 = GPOINTER_TO_INT(it->data);
                double dbl_child_lb;
                wctdata *duetime_child = (wctdata *) NULL;
                wctdata *releasetime_child = (wctdata *) NULL;

                if (dbg_lvl() == 0) {
                    printf("Creating branches for v1 = %d, C = %d.\n", v1, completiontime[v1]);
                }

                /* Create duetime and releasetime */
                rval = create_duetime_ahv(problem, pd, &duetime_child, v1, completiontime[v1]);
                CCcheck_val_2(rval, "Failed in create_duetime");
                rval = create_releasetime_ahv(problem, pd, &releasetime_child, v1, completiontime[v1]);
                CCcheck_val_2(rval, "Failed in create_differ");
                duetime_child->maxiterations = 5;
                releasetime_child->maxiterations = 5;
                compute_lower_bound(problem, duetime_child);
                compute_lower_bound(problem, releasetime_child);
                dbl_child_lb = (duetime_child->LP_lower_bound_dual < releasetime_child->LP_lower_bound_dual) ?
                               duetime_child->LP_lower_bound_dual : releasetime_child->LP_lower_bound_dual;

                if (dbl_child_lb > strongest_dbl_lb) {
                    strongest_dbl_lb = dbl_child_lb;
                    *strongest_v1     = v1;
                    remaining_branches = max_non_improving_branches;

                    if (pd->duetime_child) {
                        wctdata_free(pd->duetime_child);
                        free(pd->duetime_child);
                        pd->nduetime = 0;
                    }

                    pd->nduetime = 1;
                    pd->duetime_child = duetime_child;
                    duetime_child->maxiterations = 1000000;

                    if (pd->releasetime_child) {
                        wctdata_free(pd->releasetime_child);
                        free(pd->releasetime_child);
                        pd->nreleasetime = 0;
                    }

                    pd->nreleasetime = 1;
                    pd->releasetime_child = releasetime_child;
                    releasetime_child->maxiterations = 1000000;
                } else {
                    wctdata_free(duetime_child);
                    free(duetime_child);
                    wctdata_free(releasetime_child);
                    free(releasetime_child);
                }

                it = it->next;
            }

            if (dbg_lvl() > 0) {
                printf("Found strongest child bound of %f for v1 = %d.\n",
                       strongest_dbl_lb, *strongest_v1);
            }

            break;

        case no_strong_branching:
            v1 = GPOINTER_TO_INT(it->data);

            if (dbg_lvl() == 0) {
                printf("Creating branches for v1 = %d, C = %d.\n", v1, completiontime[v1]);
            }

            rval = create_duetime_ahv(problem, pd, &(pd->duetime_child), v1, completiontime[v1]);
            CCcheck_val_2(rval, "Failed in create_duetime");
            pd->nduetime = 1;
            rval = create_releasetime_ahv(problem, pd, &(pd->releasetime_child), v1, completiontime[v1]);
            CCcheck_val_2(rval, "Failed in create_differ");
            pd->nreleasetime = 1;
            break;
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
        if (*new_ccount == 0 || Scheduleset_less(&(cclasses[*new_ccount - 1]), &(cclasses[i]))) {
            (*new_ccount)++;

            if (*new_ccount  < i + 1) {
                Scheduleset_SWAP(&(cclasses[*new_ccount - 1]), &(cclasses[i]), &temp);
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

    if (dbg_lvl() > 1) {
        for (i = 0 ; i < pd->ccount; ++i) {
            printf("TRANSSORT SET ");

            for (j = 0; j < pd->cclasses[i].count; ++j) {
                printf(" %d", pd->cclasses[i].members[j]);
            }

            printf("\n");
        }
    }

    return val;
}

int create_branches_ahv(wctdata *pd, wctproblem *problem)
{
    int val = 0;
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
    val = test_theorem_ahv(pd, &branchjobs, &min_completiontime);
    CCcheck_val_2(val, "Failed in test ahv");

    if (g_list_length(branchjobs) == 0) {
        printf("LP returned integral solution\n");
        val = grab_integral_solution_ahv(pd, min_completiontime);
        CCcheck_val_2(val, "Failed in grab_int_sol");
        assert(pd->status = finished);
        goto CLEAN;
    }

    val = find_strongest_children_ahv(&strongest_v1, pd, problem, branchjobs, min_completiontime);
    CCcheck_val_2(val, "Failed in find_strongest_children");
    // val = create_duetime_ahv(problem, pd, strongest_v1, min_completiontime[strongest_v1]);
    // CCcheck_val(val, "Failed in create_same");
    val = set_id_and_name(pd->duetime_child, problem->nwctdata++, pd->pname);
    CCcheck_val_2(val, "Failed in set_id_and_name");
    val = compute_lower_bound(problem, pd->duetime_child);
    CCcheck_val_2(val, "Failed in compute_lower_bound");
    // val = create_releasetime_ahv(problem, pd, strongest_v1, min_completiontime[strongest_v1]);
    // CCcheck_val_2(val, "Failed in create_differ");
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

static int trigger_lb_changes_ahv(wctdata *child)
{
    int val = 0;
    int i;
    int new_lower_bound = child->lower_bound;
    wctdata *pd = (wctdata *) child->parent;

    while (pd) {
        for (i = 0;  i < pd->nduetime; ++i) {
            if (pd->duetime_child[i].lower_bound < new_lower_bound) {
                new_lower_bound = pd->duetime_child[i].lower_bound;
            }
        }

        for (i = 0;  i < pd->nreleasetime; ++i) {
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
        } else  {
            pd = (wctdata *) NULL;
        }
    }

    return val;
}

static int trigger_lb_changes_conflict(wctdata *child)
{
    int val = 0;
    int i;
    int new_lower_bound = child->lower_bound;
    wctdata *pd = (wctdata *) child->parent;

    while (pd) {
        for (i = 0;  i < pd->nsame; ++i) {
            if (pd->same_children[i].lower_bound < new_lower_bound) {
                new_lower_bound = pd->same_children[i].lower_bound;
            }
        }

        for (i = 0;  i < pd->ndiff; ++i) {
            if (pd->diff_children[i].lower_bound < new_lower_bound) {
                new_lower_bound = pd->diff_children[i].lower_bound;
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
        } else  {
            pd = (wctdata *) NULL;
        }
    }

    return val;
}

static int grab_int_sol(wctdata *pd, double *x, double tolerance)
{
    int val = 0;
    double test_incumbent = .0;
    double incumbent;
    int *colored = (int *) NULL;
    int i;
    int tot_weighted = 0;
    colored = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(colored, "Failed to allocate colored");
    fill_int(colored, pd->njobs, 0);
    val = wctlp_objval(pd->LP, &incumbent);
    CCcheck_val_2(val, "wctlp_objval failed");
    Schedulesets_free(&(pd->bestcolors), &(pd->nbbest));
    pd->bestcolors = CC_SAFE_MALLOC(pd->nmachines, Scheduleset);
    CCcheck_NULL_2(pd->bestcolors, "Failed to realloc pd->bestcolors");
    pd->nbbest = 0;

    for (i = 0; i < pd->ccount; ++i) {
        test_incumbent += x[i];

        if (x[i] >= 1.0  - tolerance) {
            int j = pd->nbbest;
            int k;
            Scheduleset_init(pd->bestcolors + j);
            pd->bestcolors[j].members = CC_SAFE_MALLOC(pd->cclasses[i].count, int);
            CCcheck_NULL_2(pd->bestcolors[j].members,
                           "Failed to realloc pd->bestcolors[j].members");

            for (k = 0; k < pd->cclasses[i].count; ++k) {
                if (!colored[pd->cclasses[i].members[k]]) {
                    colored[pd->cclasses[i].members[k]] = 1;
                    pd->bestcolors[j].members[pd->bestcolors[j].count++] = pd->cclasses[i].members[k];
                    pd->bestcolors[j].totweight += pd->duration[pd->cclasses[i].members[k]];
                    pd->bestcolors[j].totwct += pd->weights[pd->cclasses[i].members[k]] * pd->bestcolors[j].totweight;
                }
            }

            pd->nbbest++;
            tot_weighted += pd->bestcolors[j].totwct;

            if (pd->nbbest > pd->nmachines) {
                printf("ERROR: \"Integral\" solution turned out to be not integral!\n");
                fflush(stdout);
                val = 1;
                goto CLEAN;
            }
        }
    }

    val = Scheduleset_check(pd->bestcolors, pd->nbbest, pd->njobs);
    CCcheck_val_2(val, "ERROR: An incorrect coloring was created.");
    printf("Intermediate schedule:\n");
    print_schedule(pd->bestcolors, pd->nbbest);
    printf("with total weight %d\n", tot_weighted);
    assert(fabs((double)tot_weighted -  incumbent) <=
           integral_incumbent_tolerance);

    if (tot_weighted < pd->upper_bound) {
        pd->upper_bound = tot_weighted;
        pd->besttotwct = tot_weighted;
    }

    if (pd->upper_bound == pd->lower_bound) {
        pd->status = finished;
    }

CLEAN:
    CC_IFFREE(colored, int);
    return val;
}

static int insert_frac_pairs_into_heap(wctdata *pd, const double x[],
                                       int          *nodepair_refs,
                                       double       *nodepair_weights,
                                       int           npairs,
                                       pmcheap *heap)
{
    int val = 0;
    int i;
    int ref_key;

    for (i = 0; i < pd->ccount; ++i) {
        int j;

        if (x[i] <= 0.0  || x[i] >= 1.0 - lp_int_tolerance()) {
            continue;
        }

        for (j = 0; j < pd->cclasses[i].count; ++j) {
            int v1 = pd->cclasses[i].members[j];
            int k;
            ref_key  = nodepair_ref_key(v1, v1);
            nodepair_weights[ref_key] += x[i];

            for (k = j + 1 ; k < pd->cclasses[i].count; ++k) {
                assert(k != j);
                int v2 = pd->cclasses[i].members[k];
                assert(v1 < v2);
                ref_key  = nodepair_ref_key(v1, v2);
                nodepair_weights[ref_key] += x[i];
            }
        }
    }

    for (ref_key = 0; ref_key < npairs; ++ref_key) {
        int v1, v2;
        inodepair_ref_key(&v1, &v2, ref_key);

        if (v1 != v2 && nodepair_weights[ref_key] > 0.0) {
            int v1_key  = nodepair_ref_key(v1, v1);
            int v2_key  = nodepair_ref_key(v2, v2);
            double denom        = (nodepair_weights[v1_key] + nodepair_weights[v2_key]) /
                                  2;
            double dbl_heap_key = nodepair_weights[ref_key] / denom;
            int    int_heap_key =  x_frac(dbl_heap_key);
            val = pmcheap_insert(heap, int_heap_key + 1,
                                 (void *) & (nodepair_refs[ref_key]));
            CCcheck_val_2(val, "Failed in pmcheap_insert");
        }
    }

    if (dbg_lvl()) {
        printf("Size of frac heap is %d\n", pmcheap_size(heap));
    }

CLEAN:
    return val;
}

int create_branches_conflict(wctdata *pd, wctproblem *problem)
{
    int val = 0;
    int result = DELEYAD;
    int status;
    int i;
    double *x = (double *)NULL;
    wctparms *parms = &(problem->parms);
    int strongest_v1 = -1, strongest_v2 = -1;
    int *nodepair_refs = (int *) NULL;
    double *nodepair_weights = (double *) NULL;
    int npairs = pd->njobs * (pd->njobs + 1) / 2;
    int *mf_col = (int *)NULL;
    GList *branchjobs = (GList *) NULL;
    int *completion_time = (int *) NULL;
    pmcheap *heap = (pmcheap *)NULL;
    val = pmcheap_init(&heap, npairs);
    CCcheck_val_2(val, "Failed pmcheap_init");
    nodepair_refs = CC_SAFE_MALLOC(npairs, int);
    CCcheck_NULL_2(nodepair_refs, "Failed to allocate memory to nodepair_refs");
    nodepair_weights = CC_SAFE_MALLOC(npairs, double);
    CCcheck_NULL_2(nodepair_weights,
                   "Failed to allocate memory to nodepair_weights");

    for (i = 0; i < npairs; i++) {
        nodepair_refs[i] = -1;
        nodepair_weights[i] = .0;
    }

    mf_col = CC_SAFE_MALLOC(pd->njobs, int);
    CCcheck_NULL_2(mf_col, "Failed to allocate memory to mf_col");

    for (i = 0; i < pd->njobs; i++) {
        mf_col[i] = -1;
    }

    if (!pd->LP) {
        val = build_lp(pd, parms->construct);
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

    if (status == GRB_INFEASIBLE) {
        goto CLEAN;
    }

    val = wctlp_x(pd->LP, x, 0);
    CCcheck_val_2(val, "Failed at wctlp_x");
    CC_IFFREE(pd->x, double);
    pd->x = CC_SAFE_MALLOC(pd->ccount, double);
    CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
    memcpy(pd->x, x, pd->ccount * sizeof(double));

    if (parms->test_ahv) {
        val = test_theorem_ahv(pd, &branchjobs, &completion_time);
        CCcheck_val_2(val, "failed in test_theorem_ahv");

        if (g_list_length(branchjobs) == 0) {
            printf("Test ahv found a integral solution:\n");
            val = grab_integral_solution_ahv(pd, completion_time);
            CCcheck_val_2(val, "Failed in grab integral solution ahv");
            assert(pd->status == finished);
            goto CLEAN;
        }
    }

    val = insert_frac_pairs_into_heap(pd, pd->x, nodepair_refs, nodepair_weights, npairs, heap);
    CCcheck_val_2(val, "Failed in insert_frac_pairs_into_heap");

    if (pmcheap_size(heap) == 0) {
        printf("LP returned integral solution\n");
        val = grab_int_sol(pd, x, lp_int_tolerance());
        CCcheck_val_2(val, "Failed in grab_int_sol");
        assert(pd->status = finished);
        goto CLEAN;
    }

    if (pd->depth % 5 == 0 && parms->diving_heuristic) {
        parms->stab_technique = no_stab;
        val = heur_exec(problem, pd, &result);
        CCcheck_val_2(val, "Failed at heur_exec");

        if (result == FOUNDSOL) {
            printf("Heuristic found solution\n");
        }

        parms->stab_technique = stab_wentgnes;
    }

    if (dbg_lvl() > 1) {
        printf("Collected %d branching candidates.\n", pmcheap_size(heap));
    }

    val = find_strongest_children_conflict(&strongest_v1, &strongest_v2, pd, problem, heap, nodepair_refs, nodepair_weights);
    CCcheck_val_2(val, "Failed in find_strongest_children");
    // val = create_same_conflict(problem, pd, strongest_v1, strongest_v2);
    // CCcheck_val(val, "Failed in create_same");
    val = set_id_and_name(pd->same_children, problem->nwctdata++, pd->pname);
    CCcheck_val_2(val, "Failed in set_id_and_name");
    val = compute_lower_bound(problem, pd->same_children);
    CCcheck_val_2(val, "Failed in compute_lower_bound");
    // val = create_differ_conflict(problem, pd, strongest_v1, strongest_v2);
    // CCcheck_val_2(val, "Failed in create_differ");
    val = set_id_and_name(pd->diff_children, problem->nwctdata++, pd->pname);
    CCcheck_val_2(val, "Failed in set_id_and_name");
    val = compute_lower_bound(problem, pd->diff_children);
    CCcheck_val_2(val, "Failed in compute_lower_bound");
    free_elist(pd->same_children, &(problem->parms));
    free_elist(pd->diff_children, &(problem->parms));
CLEAN:
    lpwctdata_free(pd);
    free_elist(pd, &(problem->parms));

    if (heap) {
        pmcheap_free(heap);
        heap = (pmcheap *)NULL;
    }

    CC_IFFREE(x, double);
    CC_IFFREE(mf_col, int);
    CC_IFFREE(nodepair_refs, int);
    CC_IFFREE(nodepair_weights, double);
    CC_IFFREE(completion_time, int);
    g_list_free(branchjobs);
    return val;
}

int create_branches_wide(wctdata *pd, wctproblem *problem)
{
    int val = 0;
    int result = DELEYAD;
    int status;
    int i;
    double *x = (double *)NULL;
    int *v1_wide = (int *) NULL;
    int *v2_wide = (int *) NULL;
    int *min_nodepair;
    int nb_wide;
    wctparms *parms = &(problem->parms);
    int *nodepair_refs = (int *) NULL;
    double *nodepair_weights = (double *) NULL;
    int npairs = pd->njobs * (pd->njobs + 1) / 2;
    int *mf_col = (int *)NULL;
    pmcheap *heap = (pmcheap *)NULL;
    val = pmcheap_init(&heap, npairs);
    CCcheck_val_2(val, "Failed pmcheap_init");
    nodepair_refs = CC_SAFE_MALLOC(npairs, int);
    CCcheck_NULL_2(nodepair_refs, "Failed to allocate memory to nodepair_refs");
    nodepair_weights = CC_SAFE_MALLOC(npairs, double);
    CCcheck_NULL_2(nodepair_weights,
                   "Failed to allocate memory to nodepair_weights");

    for (i = 0; i < npairs; i++) {
        nodepair_refs[i] = -1;
        nodepair_weights[i] = .0;
    }

    if (!pd->LP) {
        val = build_lp(pd, parms->construct);
        CCcheck_val_2(val, "Failed at build_lp");
    }

    if (!pd->ccount) {
        compute_lower_bound(problem, pd);
    }

    assert(pd->ccount != 0);
    val = wctlp_optimize(pd->LP, &status);
    CCcheck_val_2(val, "Failed at wctlp_optimize");

    if (status == GRB_INFEASIBLE) {
        goto CLEAN;
    }

    x = CC_SAFE_MALLOC(pd->ccount, double);
    CCcheck_NULL_2(x, "Failed to allocate memory to x");
    val = wctlp_x(pd->LP, x, 0);
    CCcheck_val_2(val, "Failed at wctlp_x");
    CC_IFFREE(pd->x, double);
    pd->x = CC_SAFE_MALLOC(pd->ccount, double);
    CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
    memcpy(pd->x, x, pd->ccount * sizeof(double));
    val = insert_frac_pairs_into_heap(pd, pd->x, nodepair_refs, nodepair_weights, npairs, heap);
    CCcheck_val_2(val, "Failed in insert_frac_pairs_into_heap");

    if (pmcheap_size(heap) == 0) {
        printf("LP returned integral solution\n");
        val = grab_int_sol(pd, x, lp_int_tolerance());
        CCcheck_val_2(val, "Failed in grab_int_sol");
        assert(pd->status = finished);
        goto CLEAN;
    }

    if (pd->depth % 5 == 0 && parms->diving_heuristic) {
        val = heur_exec(problem, pd, &result);
        CCcheck_val_2(val, "Failed at heur_exec");

        if (result == FOUNDSOL) {
            printf("Heuristic found solution\n");
            assert(pd->status = finished);
            goto CLEAN;
        }
    }

    if (dbg_lvl() > 1) {
        printf("Collected %d branching candidates.\n", pmcheap_size(heap));
    }

    nb_wide = CC_MIN(pmcheap_size(heap), 5);
    pd->same_children_wide = CC_SAFE_MALLOC(1, wctdata *);
    CCcheck_NULL_2(pd->same_children_wide, "Failed to allocate memory");
    pd->diff_children_wide = CC_SAFE_MALLOC(nb_wide, wctdata *);
    CCcheck_NULL_2(pd->diff_children_wide, "Failed to allocate memory");
    v1_wide = CC_SAFE_MALLOC(nb_wide, int);
    CCcheck_NULL_2(v1_wide, "Failed to allocate");
    v2_wide = CC_SAFE_MALLOC(nb_wide, int);
    CCcheck_NULL_2(v2_wide, "Failed to allocate");

    for (i = 0; i < nb_wide; i++) {
        min_nodepair = (int *) pmcheap_min(heap);
        int v1 = -1, v2 = -1;
        inodepair_ref_key(&v1, &v2, (int)(min_nodepair - nodepair_refs));
        assert(v1 < v2);
        v1_wide[i] = v1;
        v2_wide[i] = v2;

        if (dbg_lvl() > 1) {
            printf("Inserted  v1 = %d and v2 = %d with nodepair_weight = %f .\n",
                   v1, v2,
                   nodepair_weights[(int)(min_nodepair - nodepair_refs)]);
        }
    }

    val = create_same_wide(problem, pd, v1_wide, v2_wide, nb_wide);
    CCcheck_val(val, "Failed in create_same");
    val = set_id_and_name(*(pd->same_children_wide), problem->nwctdata++, pd->pname);
    CCcheck_val_2(val, "Failed in set_id_and_name");
    val = compute_lower_bound(problem, *(pd->same_children_wide));
    CCcheck_val_2(val, "Failed in compute_lower_bound");

    for (i = 0; i < nb_wide; ++i) {
        create_differ_wide(problem, pd, v1_wide[i], v2_wide[i]);
        CCcheck_val_2(val, "Failed in create_differ");
        val = set_id_and_name(pd->diff_children_wide[i], problem->nwctdata++, pd->pname);
        CCcheck_val_2(val, "Failed in set_id_and_name");
        val = compute_lower_bound(problem, pd->diff_children_wide[i]);
        CCcheck_val_2(val, "Failed in compute_lower_bound");
        free_elist(pd->diff_children_wide[i], &(problem->parms));
    }

    free_elist(*(pd->same_children_wide), &(problem->parms));
CLEAN:
    lpwctdata_free(pd);
    free_elist(pd, &(problem->parms));

    if (heap) {
        pmcheap_free(heap);
        heap = (pmcheap *)NULL;
    }

    CC_IFFREE(v1_wide, int);
    CC_IFFREE(v2_wide, int);
    CC_IFFREE(x, double);
    CC_IFFREE(mf_col, int);
    CC_IFFREE(nodepair_refs, int);
    CC_IFFREE(nodepair_weights, double);
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
    wctparms *parms = &(problem->parms);
    problem->parms.bb_search_strategy = dfs_strategy;

    switch (parms->bb_search_strategy) {
        case dfs_strategy:
            if (pd->parent) {
                heap_key = (int)(pd->lower_bound) - pd->depth * 10000 - pd->id % 2 ;
            }

            break;

        case min_lb_strategy:
        default:
            heap_key = (int)(pd->LP_lower_bound * problem->mult_key) -
                       pd->id % 2;
    }

    int lb = pd->lower_bound;

    if (dbg_lvl()) {
        printf("Inserting into branching heap with lb %d and ub %d at depth %d (id = %d) heap_key = %d\n",
               pd->lower_bound, pd->upper_bound, pd->depth, pd->id, heap_key);
    }

    free_elist(pd, &(problem->parms));

    if (lb < pd->upper_bound) {
        val = pmcheap_insert(problem->br_heap, heap_key, (void *)pd);
        CCcheck_val(val, "Failed at pmcheap_insert");
    } else {
        skip_wctdata(pd, problem);
    }

    switch (parms->bb_branch_strategy) {
        case conflict_strategy:
            val = trigger_lb_changes_conflict(pd);
            CCcheck_val_2(val, "Failed in trigger_lb_changes_conflict");
            break;

        case ahv_strategy:
            val = trigger_lb_changes_ahv(pd);
            CCcheck_val_2(val, "Failed  in trigger_lb_changes_ahv")
            break;
    }

CLEAN:
    return val;
}

int branching_msg_wide(wctdata *pd, wctproblem *problem)
{
    pmcheap *br_heap = problem->br_heap;

    if (pd->lower_bound < pd->upper_bound) {
        CCutil_suspend_timer(&problem->tot_cputime);
        printf("Branching with lb %d (LP %f) at depth %d (id = %d, "
               "time = %f, unprocessed nodes = %d, nbjobs= %d, upper bound = %d, lower bound = %d, v1 = %d, v2 = %d, nbdiff = %d, nbsame = %d ).\n",
               pd->lower_bound, pd->LP_lower_bound,
               pd->depth,
               pd->id, problem->tot_cputime.cum_zeit, pmcheap_size(br_heap), pd->njobs, problem->global_upper_bound, problem->global_lower_bound, pd->v1, pd->v2, pd->ecount_differ, pd->ecount_same
              );
        CCutil_resume_timer(&problem->tot_cputime);
    }

    return 0;
}

int branching_msg(wctdata *pd, wctproblem *problem)
{
    pmcheap *br_heap = problem->br_heap;

    if (pd->lower_bound < pd->upper_bound) {
        CCutil_suspend_timer(&problem->tot_cputime);
        printf("Branching with lb %d (LP %f) at depth %d (id = %d, "
               "time = %f, unprocessed nodes = %d, nbjobs= %d, upper bound = %d, lower bound = %d, v1 = %d, v2 = %d, nbdiff = %d, nbsame = %d ).\n",
               pd->lower_bound, pd->LP_lower_bound,
               pd->depth,
               pd->id, problem->tot_cputime.cum_zeit, pmcheap_size(br_heap), pd->njobs, problem->global_upper_bound, problem->global_lower_bound, pd->v1, pd->v2, pd->ecount_differ, pd->ecount_same
              );
        CCutil_resume_timer(&problem->tot_cputime);
    }

    return 0;
}

int branching_msg_ahv(wctdata *pd, wctproblem *problem)
{
    pmcheap *br_heap = problem->br_heap;

    if (pd->lower_bound < pd->upper_bound) {
        CCutil_suspend_timer(&problem->tot_cputime);
        printf("Branching with lb %d (LP %f) at depth %d (id = %d, "
               "time = %f, unprocessed nodes = %d, nbjobs= %d, upper bound = %d, lower bound = %d, branch_job = %d, completion_time= %d).\n",
               pd->lower_bound, pd->LP_lower_bound,
               pd->depth,
               pd->id, problem->tot_cputime.cum_zeit, pmcheap_size(br_heap), pd->njobs, problem->global_upper_bound, problem->global_lower_bound, pd->branch_job, pd->completiontime);
        CCutil_resume_timer(&problem->tot_cputime);
    }

    return 0;
}

int sequential_branching_ahv(wctproblem *problem)
{
    int val = 0;
    wctdata *pd;
    pmcheap *br_heap = problem->br_heap;
    wctparms *parms = &(problem->parms);
    printf("ENTERED SEQUANTIAL BRANCHING AHV:\n");
    CCutil_suspend_timer(&problem->tot_branch_and_bound);

    while ((pd = (wctdata *) pmcheap_min(br_heap))
            && problem->tot_branch_and_bound.cum_zeit < parms->branching_cpu_limit) {
        CCutil_resume_timer(&problem->tot_branch_and_bound);
        int i;
        pd->upper_bound = problem->global_upper_bound;

        if (pd->lower_bound >= pd->upper_bound || pd->status == infeasible || pd->eta_in > pd->upper_bound - 1) {
            skip_wctdata(pd, problem);
            remove_finished_subtree_ahv(pd);
        } else {
            branching_msg_ahv(pd, problem);
            /** Construct PricerSolver */
            /*val = recover_elist(pd);
            CCcheck_val_2(val, "Failed in recover_elist");*/

            if (problem->maxdepth < pd->depth) {
                problem->maxdepth = pd->depth;
            }

            val = create_branches_ahv(pd, problem);
            CCcheck_val_2(val, "Failed at create_branches");

            for (i = 0; i < pd->nduetime; i++) {
                val = insert_into_branching_heap(pd->duetime_child + i, problem);
                CCcheck_val_2(val, "Faield at insert_into_branching_heap");
            }

            for (i = 0; i < pd->nreleasetime; i++) {
                val = insert_into_branching_heap(&(pd->releasetime_child[i]), problem);
                CCcheck_val_2(val, "Failed in insert_into_branching_heap");
            }

            assert(pd->lower_bound <= pd->upper_bound);
            adapt_global_upper_bound(problem, pd->upper_bound);

            if (pd->upper_bound == pd->lower_bound) {
                remove_finished_subtree_ahv(pd);
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

int sequential_branching_conflict(wctproblem *problem)
{
    int val = 0;
    wctdata *pd;
    pmcheap *br_heap = problem->br_heap;
    wctparms *parms = &(problem->parms);
    printf("ENTERED SEQUANTIAL BRANCHING CONFLICT:\n");
    CCutil_suspend_timer(&problem->tot_branch_and_bound);

    while ((pd = (wctdata *) pmcheap_min(br_heap))
            && problem->tot_branch_and_bound.cum_zeit < parms->branching_cpu_limit) {
        CCutil_resume_timer(&problem->tot_branch_and_bound);
        int i;
        pd->upper_bound = problem->global_upper_bound;

        if (pd->lower_bound >= pd->upper_bound || pd->status == infeasible || pd->eta_in > pd->upper_bound - 1) {
            skip_wctdata(pd, problem);
            remove_finished_subtree_conflict(pd);
        } else {
            branching_msg(pd, problem);
            /** Construct PricerSolver */
            /*val = recover_elist(pd);
            CCcheck_val_2(val, "Failed in recover_elist");*/

            if (problem->maxdepth < pd->depth) {
                problem->maxdepth = pd->depth;
            }

            val = create_branches_conflict(pd, problem);
            CCcheck_val_2(val, "Failed at create_branches");

            for (i = 0; i < pd->nsame; i++) {
                val = insert_into_branching_heap(&(pd->same_children[i]), problem);
                CCcheck_val_2(val, "Failed in insert_into_branching_heap");
            }

            for (i = 0; i < pd->ndiff; i++) {
                val = insert_into_branching_heap(pd->diff_children + i, problem);
                CCcheck_val_2(val, "Faield at insert_into_branching_heap");
            }

            assert(pd->lower_bound <= pd->upper_bound);
            adapt_global_upper_bound(problem, pd->upper_bound);

            if (pd->upper_bound == pd->lower_bound) {
                remove_finished_subtree_conflict(pd);
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

int sequential_branching_wide(wctproblem *problem)
{
    int val = 0;
    wctdata *pd;
    pmcheap *br_heap = problem->br_heap;
    wctparms *parms = &(problem->parms);
    printf("ENTERED SEQUANTIAL WIDE BRANCHING:\n");
    CCutil_suspend_timer(&problem->tot_branch_and_bound);

    while ((pd = (wctdata *) pmcheap_min(br_heap))
            && problem->tot_branch_and_bound.cum_zeit < parms->branching_cpu_limit) {
        CCutil_resume_timer(&problem->tot_branch_and_bound);
        int i;
        pd->upper_bound = problem->global_upper_bound;

        if (pd->lower_bound >= pd->upper_bound || pd->status == infeasible || pd->eta_in > pd->upper_bound - 1) {
            skip_wctdata(pd, problem);
            remove_finished_subtree_wide(pd);
        } else {
            branching_msg_wide(pd, problem);
            /** Construct PricerSolver */
            /*val = recover_elist(pd);
            CCcheck_val_2(val, "Failed in recover_elist");*/

            if (problem->maxdepth < pd->depth) {
                problem->maxdepth = pd->depth;
            }

            val = create_branches_wide(pd, problem);
            CCcheck_val_2(val, "Failed at create_branches");

            for (i = 0; i < pd->nsame; i++) {
                val = insert_into_branching_heap((pd->same_children_wide[i]), problem);
                CCcheck_val_2(val, "Failed in insert_into_branching_heap");
            }

            for (i = 0; i < pd->ndiff; i++) {
                val = insert_into_branching_heap(pd->diff_children_wide[i], problem);
                CCcheck_val_2(val, "Faield at insert_into_branching_heap");
            }

            assert(pd->lower_bound <= pd->upper_bound);
            adapt_global_upper_bound(problem, pd->upper_bound);

            if (pd->upper_bound == pd->lower_bound) {
                remove_finished_subtree_wide(pd);
            }
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

MAYBE_UNUSED
static int grow_ages(wctdata *pd)
{
    int val = 0;
    int i;
    int *cstat;
    cstat = (int *)CC_SAFE_MALLOC(pd->ccount, int);
    CCcheck_NULL_2(cstat, "Failed to allocate cstat");
    val = wctlp_basis_cols(pd->LP, cstat, 0);
    CCcheck_val_2(val, "Failed in wctlp_basis_cols");
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

MAYBE_UNUSED
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
                    CCcheck_val_2(val, "Failed in wctlp_deletecols");
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
            CCcheck_val_2(val, "Failed in wctlp_deletecols");
        }

        assert(pd->dzcount == pd->ccount - new_ccount);
        CC_IFFREE(pd->cclasses, Scheduleset);
        pd->cclasses = new_cclasses;
        pd->ccount   = new_ccount;

        if (dbg_lvl() > 1) {
            printf("Deleted %d out of %d columns with age > %d.\n",
                   pd->dzcount, pd->dzcount + pd->ccount, pd->retirementage);
        }

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
                pd->nnewsets = 1;
                pd->newsets = CC_SAFE_MALLOC(1, Scheduleset);
                Scheduleset_init(pd->newsets);
                pd->newsets[0].members = CC_SAFE_MALLOC(2, int);
                CCcheck_NULL_2(pd->newsets[0].members,
                               "Failed to allocate memory to pd->newsets->members");
                pd->newsets[0].C = CC_SAFE_MALLOC(1, int);
                CCcheck_NULL_2(pd->newsets[0].C, "Failed to allocate memory");
                pd->newsets[0].table = g_hash_table_new(g_direct_hash, g_direct_equal);
                CCcheck_NULL_2(pd->newsets[0].table, "Failed to allocate memory");
                pd->newsets[0].count++;
                pd->newsets[0].members[0] = i;
                pd->newsets[0].members[1] = pd->njobs;
                pd->newsets[0].totwct = pd->weights[i] * pd->duration[i];
                pd->newsets[0].totweight = pd->duration[i];
                pd->newsets[0].C[0] = pd->duration[i];
                g_hash_table_insert(pd->newsets[0].table, GINT_TO_POINTER(i), pd->newsets[0].C);
                pd->newsets->age = 0;
                val = wctlp_addcol(pd->LP, 2, pd->newsets[0].members, pd->coef, pd->newsets[0].totwct, 0.0, 1.0,
                                   wctlp_CONT, NULL);
                CCcheck_val_2(val, "Failed in wctlp_addcol");

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
    pd->pi_in = CC_SAFE_MALLOC(pd->njobs + 1, double);
    CCcheck_NULL_2(pd->pi_in, "Failed to allocate memory");
    fill_dbl(pd->pi_in, pd->njobs + 1, 0.0);
    pd->eta_in = 0.0;
    pd->pi_out = CC_SAFE_MALLOC(pd->njobs + 1, double);
    CCcheck_NULL_2(pd->pi_out, "Failed to allocate memory");
    pd->eta_out = 0.0;
    fill_dbl(pd->pi_out, pd->njobs + 1, 0.0);
    pd->pi_sep = CC_SAFE_MALLOC(pd->njobs + 1, double);
    CCcheck_NULL_2(pd->pi_sep, "Failed to allocate memory");
    pd->subgradient_in = CC_SAFE_MALLOC(pd->njobs + 1, double);
    CCcheck_NULL_2(pd->subgradient_in, "Failed to allocate memory");
    pd->subgradient = CC_SAFE_MALLOC(pd->njobs + 1, double);
    CCcheck_NULL_2(pd->subgradient, "Failed to allocate memory");
    pd->rhs = CC_SAFE_MALLOC(pd->njobs + 1, double);
    CCcheck_NULL_2(pd->rhs, "Failed to allocate memory");
    val = wctlp_get_rhs(pd->LP, pd->rhs);
    CCcheck_val_2(val, "Failed to get RHS");
CLEAN:

    if (val) {
        wctlp_free(&(pd->LP));
        CC_IFFREE(pd->coef, double);
        CC_IFFREE(pd->pi, double);
        CC_IFFREE(pd->pi_in, double)
        CC_IFFREE(pd->pi_out, double)
        CC_IFFREE(pd->pi_sep, double)
        CC_IFFREE(pd->subgradient, double)
        CC_IFFREE(pd->subgradient_in, double)
        CC_IFFREE(pd->rhs, double)
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
    wctparms *parms = &(problem->parms);

    if (parms->print) {
        print_size_to_csv(problem, pd);
    }

    if (pd->status == infeasible) {
        goto CLEAN;
    }

    if (dbg_lvl() > 1) {
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

    pd->retirementage = (int) sqrt(pd->njobs) + 30;

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

    /** Compute LP relaxation */
    cur_cputime = CCutil_zeit();
    val = wctlp_optimize(pd->LP, &status);
    CCcheck_val_2(val, "wctlp_optimize failed");
    cur_cputime = CCutil_zeit() - cur_cputime;

    if (dbg_lvl() > 1) {
        printf("Simplex took %f seconds.\n", CCutil_zeit() - cur_cputime);
        fflush(stdout);
    }

    if (dbg_lvl() > 1) {
        print_ages(pd);
    }

    switch (status) {
        case GRB_OPTIMAL:
            /** grow ages of the different columns */
            val = grow_ages(pd);
            CCcheck_val_2(val, "Failed in grow_ages");
            /** get the dual variables and make them feasible */
            val = wctlp_pi(pd->LP, pd->pi);
            CCcheck_val_2(val, "wctlp_pi failed");
            /** Compute the objective function */
            val = compute_objective(pd, parms);
            CCcheck_val_2(val, "Failed in compute_objective");
            memcpy(pd->pi_out, pd->pi, sizeof(double) * (pd->njobs + 1));
            pd->eta_out = pd->LP_lower_bound_dual;
            break;

        case GRB_INFEASIBLE:
            /** get the dual variables and make them feasible */
            val = wctlp_pi(pd->LP, pd->pi);
            CCcheck_val_2(val, "wctlp_pi failed");
            break;
    }

    break_while_loop = 0;
    CCutil_suspend_timer(&(problem->tot_cputime));
    CCutil_resume_timer(&(problem->tot_cputime));

    while ((iterations < pd->maxiterations)
            && !break_while_loop
            && problem->tot_cputime.cum_zeit <= problem->parms.branching_cpu_limit) {
        /** delete old columns */
        // if (pd->dzcount > pd->njobs * min_ndelrow_ratio && status == GRB_OPTIMAL) {
        //     val = delete_old_cclasses(pd);
        // }
        /** Solve the pricing problem*/
        CCutil_start_resume_time(&problem->tot_pricing);

        switch (status) {
            case GRB_OPTIMAL:
                iterations++;

                if (iterations < pd->maxiterations) {
                    switch (parms->stab_technique) {
                        case stab_wentgnes:
                            val = solve_stab(pd, parms);
                            CCcheck_val_2(val, "Failed in solve_stab");
                            break;

                        case stab_dynamic:
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
                val = solve_farkas_dbl(pd);
                CCcheck_val_2(val, "Failed in solving farkas");
                break;
        }

        CCutil_suspend_timer(&problem->tot_pricing);

        for (j = 0; j < pd->nnewsets; j++) {
            val = wctlp_addcol(pd->LP, pd->newsets[j].count + 1, pd->newsets[j].members, pd->coef, pd->newsets[j].totwct, 0.0, 1.0, wctlp_CONT, NULL);
            CCcheck_val_2(val, "wctlp_addcol failed");
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

        add_newsets(pd);
        /** Compute LP relaxation */
        cur_cputime = CCutil_zeit();
        val = wctlp_optimize(pd->LP, &status);
        CCcheck_val_2(val, "wctlp_optimize failed");
        cur_cputime = CCutil_zeit() - cur_cputime;

        if (dbg_lvl() > 1) {
            printf("Simplex took %f seconds.\n", CCutil_zeit() - cur_cputime);
            fflush(stdout);
        }

        if (dbg_lvl() > 1) {
            print_ages(pd);
        }

        switch (status) {
            case GRB_OPTIMAL:
                /** grow ages of the different columns */
                val = grow_ages(pd);
                CCcheck_val_2(val, "Failed in grow_ages");
                /** get the dual variables and make them feasible */
                val = wctlp_pi(pd->LP, pd->pi);
                CCcheck_val_2(val, "wctlp_pi failed");
                /** Compute the objective function */
                val = compute_objective(pd, parms);
                CCcheck_val_2(val, "Failed in compute_objective");

                if (pd->update) {
                    memcpy(pd->pi_out, pd->pi, sizeof(double) * (pd->njobs + 1));
                    pd->eta_out = pd->LP_lower_bound_dual;
                }

                break;

            case GRB_INFEASIBLE:
                /** get the dual variables and make them feasible */
                val = wctlp_pi(pd->LP, pd->pi);
                CCcheck_val_2(val, "wctlp_pi failed");
                break;
        }

        CCutil_suspend_timer(&(problem->tot_cputime));
        CCutil_resume_timer(&(problem->tot_cputime));
    }

    if (iterations < pd->maxiterations && problem->tot_cputime.cum_zeit <= problem->parms.branching_cpu_limit) {
        switch (status) {
            case GRB_OPTIMAL:

                /** change status of problem */
                if (problem->status == no_sol) {
                    problem->status = lp_feasible;
                }

                if (dbg_lvl() > 1) {
                    printf("Found lb = %d (%f) upper_bound = %d (id= %d, iterations = %d,opt_track = %d).\n",
                           pd->lower_bound, pd->LP_lower_bound, pd->upper_bound,
                           pd->id, iterations, pd->opt_track);
                }

                pd->x = CC_SAFE_MALLOC(pd->ccount, double);
                CCcheck_NULL_2(pd->x, "Failed to allocate memory to pd->x");
                val = wctlp_x(pd->LP, pd->x, 0);
                CCcheck_val_2(val, "Failed in wctlp_x");
                pd->status = LP_bound_computed;
                break;

            case GRB_INFEASIBLE:
                printf("LP is infeasible\n");
                pd->status = infeasible;
        }
    } else  {
        pd->status = LP_bound_estimated;
    }

    if (dbg_lvl() > 1) {
        printf("iterations = %d\n", iterations);
    }

    fflush(stdout);
    CCutil_suspend_timer(&(problem->tot_lb_lp));
CLEAN:
    return val;
}



MAYBE_UNUSED static int submiping(wctdata *cd)
{
    int val = 0;
    double incumbent;
    int status;
    double *colsol;
    int *colored = (int *) NULL;
    int i;
    colsol = (double *) CC_SAFE_MALLOC(cd->ccount, double);
    CCcheck_NULL_2(colsol, "Failed to allocate colsol");
    colored = (int *) CC_SAFE_MALLOC(cd->njobs, int);
    CCcheck_NULL_2(colored, "Failed to allocate colored");

    for (i = 0; i < cd->njobs; ++i) {
        colored[i] = 0;
    }

    val = wctlp_set_coltypes(cd->LP, wctlp_BIN);
    CCcheck_val_2(val,
                  "wctlp_set_all_coltypes "
                  "(this warning can be ignored if you are using "
                  "CPLEX or QSopt as an LP-solver)");
    /* wctlp_write (cd->LP, "lpheur.lp"); */
    val = wctlp_setnodelimit(cd->LP, 1);
    CCcheck_val_2(val, "wctlp_setnodelimit failed");
    val = wctlp_optimize(cd->LP, &status);
    CCcheck_val_2(val, "wctlp_optimize failed");
    val =  wctlp_x(cd->LP, colsol, 0);
    CCcheck_val_2(val, "wctlp_x failed");
    val = wctlp_objval(cd->LP, &incumbent);
    CCcheck_val_2(val, "wctlp_objval failed");
    grab_int_sol(cd, colsol, lp_int_tolerance());
    wctlp_set_coltypes(cd->LP, wctlp_CONT);
    CCcheck_val_2(val, "wctlp_set_all_coltypes");
    printf("Found lower bound of %lld and upper bound of %g.\n",
           (long long) cd->lower_bound, incumbent);
    CCcheck_val_2(val, "ERROR: An incorrect coloring was created.");
    print_schedule(cd->bestcolors, cd->nbbest);
    printf("status = %d\n", status);
    getchar();
CLEAN:
    CC_IFFREE(colsol, double)
    CC_IFFREE(colored, int)
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
    wctparms *parms = &(problem->parms);
    problem->mult_key = 1.0;
    problem->first_upper_bound = problem->global_upper_bound;
    problem->first_lower_bound = problem->global_lower_bound;
    problem->first_rel_error = (double)(problem->global_upper_bound -
                                        problem->global_lower_bound) / ((double)problem->global_lower_bound);
    prune_duplicated_sets(root_pd);

    if (root_pd->status >= LP_bound_computed) {
        val = prefill_heap(root_pd, problem);
        CCcheck_val(val, "Failed in prefill_heap");
    } else {
        CCutil_start_timer(&(problem->tot_lb_lp_root));
        val = compute_lower_bound(problem, root_pd);
        CCcheck_val_2(val, "Failed in compute_lower_bound");

        if (root_pd->lower_bound > problem->global_lower_bound) {
            problem->global_lower_bound = root_pd->lower_bound;
            problem->first_lower_bound = root_pd->lower_bound;
            problem->first_rel_error = (double)(problem->first_upper_bound -
                                                problem->first_lower_bound) / ((double)problem->first_lower_bound);
        }

        problem->parms.construct = 1;
        CCcheck_val_2(val, "Failed in compute_lower_bound");
        CCutil_stop_timer(&(problem->tot_lb_lp_root), 0);
        val = insert_into_branching_heap(root_pd, problem);
        CCcheck_val_2(val, "insert_into_branching_heap failed");
    }

    CCutil_start_resume_time(&(problem->tot_branch_and_bound));

    switch (parms->bb_branch_strategy) {
        case conflict_strategy:
            val = sequential_branching_conflict(problem);
            CCcheck_val(val, "Failed in sequential_branching_conflict");
            break;

        case ahv_strategy:
            val = sequential_branching_ahv(problem);
            CCcheck_val_2(val, "Failed in sequential_branching_ahv");
            break;
    }

    CCutil_stop_timer(&(problem->tot_branch_and_bound), 0);
    printf("Compute schedule finished with LB %d and UB %d\n", root_pd->lower_bound, problem->global_upper_bound);

    if (root_pd->lower_bound == problem->global_upper_bound) {
        problem->global_lower_bound = root_pd->lower_bound;
        problem->rel_error = (double)(problem->global_upper_bound -
                                      problem->global_lower_bound) / ((double)problem->global_lower_bound);
        problem->status = optimal;
        printf("The optimal schedule is given by:\n");
        print_schedule(root_pd->bestcolors, root_pd->nbbest);
        printf("with total weighted completion time %d\n", root_pd->upper_bound);
    } else {
        problem->global_lower_bound = root_pd->lower_bound;
        problem->rel_error = (double)(problem->global_upper_bound -
                                      problem->global_lower_bound) / ((double)problem->global_lower_bound);
        problem->status = meta_heur;
        problem->global_lower_bound = root_pd->lower_bound;
        printf("The suboptimal schedule is given by:\n");
        print_schedule(root_pd->bestcolors, root_pd->nbbest);
        printf("with total weighted completion time\n");
    }

CLEAN:
    return val;
}


