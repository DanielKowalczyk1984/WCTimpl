#include <glib.h>
#include <assert.h>
#include <math.h>
#include "wct.h"

int main(int argc, char const *argv[])
{
    int i;
    int val = 0;
    int njobs = 100 ;
    int nmachine = 5;
    int *weights = (int *) NULL;
    int *durations = (int *) NULL;
    int *perm = (int *) NULL;
    double *ratio = (double *) NULL;
    Job *jobarray = (Job *) NULL;
    solution *sol = CC_SAFE_MALLOC(1, solution);
    solution *sol2 = CC_SAFE_MALLOC(1, solution);
    GRand *rand_ = g_rand_new_with_seed(1984);

    jobarray = CC_SAFE_MALLOC(njobs, Job);
    CCcheck_NULL_2(jobarray, "Failed to allocate");
    durations = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(durations, "Failed to allocate")
    weights = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(weights, "Failed to allocate")
    ratio = CC_SAFE_MALLOC(njobs, double);
    CCcheck_NULL_2(ratio, "Failed to allocate")
    perm = CC_SAFE_MALLOC(njobs, int);
    CCcheck_NULL_2(perm, "Failed to allocate")


    for(i = 0; i < njobs; ++i) {
        durations[i]= g_random_int_range(1, 100);
        weights[i] = g_random_int_range(1, 20);
        ratio[i] = (double) durations[i]/(double) weights[i];
        perm[i] = i;
    }

    CCutil_double_perm_quicksort(perm, ratio, njobs);

        
    for (i = 0; i < njobs; ++i)
    {
      jobarray[i].weight = weights[perm[i]];
      jobarray[i].processingime = durations[perm[i]];
      jobarray[i].job = i;
    }

    int lowerbound = lowerbound_eei(jobarray, njobs, nmachine);
    
    solution_init(sol);
    solution_alloc(sol, nmachine, njobs);
    solution_init(sol2);
    solution_alloc(sol2, nmachine, njobs);
    construct_wspt(jobarray, njobs, nmachine, sol);
    construct_wspt(jobarray, njobs, nmachine, sol2);
    solution_calc(sol, jobarray);
    solution_calc(sol2, jobarray);
    printf("first solution1 = %d\n", sol->totalweightcomptime);
    printf("first solution2 = %d\n", sol2->totalweightcomptime);

    local_search_machine_general_first(jobarray, sol, 0, 1, 0);
    local_search_machine_general_first(jobarray, sol, 0, 1, 1);
    local_search_machine_general_first(jobarray, sol, 0, 2, 1);
    local_search_machine_general_first(jobarray, sol, 0, 2, 2);
    
    local_search_machine_general_best(jobarray, sol2, 0, 1, 0);
    local_search_machine_general_best(jobarray, sol2, 0, 1, 1);
    local_search_machine_general_best(jobarray, sol2, 0, 2, 1);
    local_search_machine_general_best(jobarray, sol2, 0, 2, 2);
    solution_wct(sol, jobarray);
    solution_wct(sol2, jobarray);
    printf("final solution1  = %d\n", sol->totalweightcomptime);
    printf("final solution2  = %d\n", sol2->totalweightcomptime);
    solution_free(sol);
    solution_free(sol2);

    CLEAN:
    CC_IFFREE(durations, int)
    CC_IFFREE(weights, int)
    CC_IFFREE(ratio, double)
    CC_IFFREE(perm, int)
    CC_IFFREE(jobarray, Job);
    CC_IFFREE(sol2, solution);
    CC_IFFREE(sol, solution);
    g_rand_free(rand_);
    return val;
}