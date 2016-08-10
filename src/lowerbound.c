#include <math.h>
#include "wct.h"

int gcd(int x, int y);
int gcd_weight(int n, Job *a);
int gcd_duration(int n, Job *a);

int gcd(int x, int y) {
    int wk;

    if (x < y) {
        wk = x;
        x = y;
        y = wk;
    }

    while (y) {
        wk = x % y;
        x = y;
        y = wk;
    }

    return x;
}


int gcd_weight(int n, Job *a) {
    if (n == 1) {
        return a[0].weight;
    }

    if (n == 2) {
        return gcd(a[0].weight, a[1].weight);
    }

    int h = n / 2;
    return gcd(gcd_weight(h, a + (h - 1)), gcd_weight(n - h, a + h));
}

int gcd_duration(int n, Job *a) {
    if (n == 1) {
        return a[0].processingime;
    }

    if (n == 2) {
        return gcd(a[0].processingime, a[1].processingime);
    }

    int h = n / 2;
    return gcd(gcd_duration(h, a + (h - 1)), gcd_duration(n - h, a + h));
}

int compare_cw(BinomialHeapValue a, BinomialHeapValue b) {
    double *aw = &(((MACHINE *)a)->totcompletion);
    double *bw = &(((MACHINE *)b)->totcompletion);

    if (*aw < *bw) {
        return -1;
    } else {
        return 1;
    }
}

int lowerbound_cw(Job *array, int njobs, int nmachines) {
    int val = 0;
    int i, j;
    double temp;
    double totweight = 0.0;
    int w = gcd_weight(njobs, array);
    int Cn = 0;
    double CN = 0.0;
    int *nbsubjobs = (int *) NULL;
    int *perm = (int *) NULL;
    double *durationsubjobs = (double *) NULL;
    double *ratio = (double *) NULL;
    MACHINE *machine = (MACHINE *) NULL;
    BinomialHeap *heapdbl;
    heapdbl = binomial_heap_new(BINOMIAL_HEAP_TYPE_MIN, compare_cw);
    nbsubjobs =  CC_SAFE_MALLOC(njobs, int);
    perm = CC_SAFE_MALLOC(njobs, int);
    durationsubjobs = CC_SAFE_MALLOC(njobs, double);
    ratio = CC_SAFE_MALLOC(njobs, double);
    machine  = CC_SAFE_MALLOC(nmachines, MACHINE);

    for (i = 0; i < nmachines; ++i) {
        machine[i].totcompletion = 0.0;
        machine[i].totweight = 0.0;
        binomial_heap_insert(heapdbl, machine + i);
    }

    for (i = 0; i < njobs; ++i) {
        nbsubjobs[i] = array[i].weight / w;
        durationsubjobs[i] = (double)array[i].processingime / (double)nbsubjobs[i];
        ratio[i] = durationsubjobs[i] / w;
        perm[i] = i;
    }

    CCutil_double_perm_quicksort(perm, ratio, njobs);

    for (i = 0; i < njobs; ++i) {
        temp = durationsubjobs[perm[i]];

        for (j = 0; j < nbsubjobs[perm[i]]; ++j) {
            MACHINE *temp_machine = binomial_heap_pop(heapdbl);
            totweight += w * ((*temp_machine).totcompletion + temp);
            (*temp_machine).totcompletion += temp;
            (*temp_machine).totweight += (double)w;
            binomial_heap_insert(heapdbl, temp_machine);
        }
    }

    for (i = 0; i < njobs; ++i) {
        Cn += array[i].processingime * array[i].weight;
        CN += (double)w * array[i].processingime;
    }

    double a = (double)Cn / 2.0 - CN / 2.0;
    val = (int) ceil(totweight + a) ;
    CC_IFFREE(ratio, double);
    CC_IFFREE(durationsubjobs, double);
    CC_IFFREE(machine, MACHINE);
    CC_IFFREE(perm, int);
    CC_IFFREE(nbsubjobs, int);
    binomial_heap_free(heapdbl);
    return val;
}

int lowerbound_cp(Job *array, int njobs, int nmachines) {
    int val = 0;
    int i, j;
    double temp;
    int p = gcd_duration(njobs, array);
    int *nbsubjobs = (int *) NULL;
    int *perm = (int *) NULL;
    double *weightsubjobs = (double *) NULL;
    double *ratio = (double *) NULL;
    double totweight = 0.0;
    double *completiontime = (double *) NULL;
    pmcheap *heap = (pmcheap *) NULL;
    int Cn = 0;
    double CN = 0.0;
    pmcheap_init(&heap, nmachines);
    nbsubjobs =  CC_SAFE_MALLOC(njobs, int);
    weightsubjobs = CC_SAFE_MALLOC(njobs, double);
    ratio = CC_SAFE_MALLOC(njobs, double);
    completiontime  = CC_SAFE_MALLOC(nmachines, double);
    perm = CC_SAFE_MALLOC(njobs, int);

    for (i = 0; i < nmachines; ++i) {
        completiontime[i] = 0;
        pmcheap_insert(heap, completiontime[i], completiontime + i);
    }

    for (i = 0; i < njobs; ++i) {
        nbsubjobs[i] = array[i].processingime / p;
        weightsubjobs[i] = (double)array[i].weight / (double)nbsubjobs[i];
        ratio[i] = (double)p / weightsubjobs[i];
        perm[i] = i;
    }

    CCutil_double_perm_quicksort(perm, ratio, njobs);

    for (i = 0; i < njobs; ++i) {
        temp = weightsubjobs[perm[i]];

        for (j = 0; j < nbsubjobs[perm[i]]; ++j) {
            int *temp_machine = pmcheap_min(heap);
            *temp_machine += p;
            totweight += temp * (*temp_machine);
            pmcheap_insert(heap, *temp_machine, temp_machine);
        }
    }

    for (i = 0; i < njobs; ++i) {
        Cn += array[i].processingime * array[i].weight;
        CN += array[i].weight * p;
    }

    double a = (double)Cn / 2.0 - CN / 2.0;
    val = (int) ceil(totweight + a);
    CC_IFFREE(ratio, double);
    CC_IFFREE(perm, int);
    CC_IFFREE(nbsubjobs, int);
    CC_IFFREE(weightsubjobs, double);
    CC_IFFREE(completiontime, double);
    pmcheap_free(heap);
    return val;
}

int lowerbound_eei(Job *array, int njobs, int nmachines) {
    int val = 0;
    int i;
    int C1 = 0;
    int Cn = 0;
    int C = 0;
    double lowerbound;

    for (i = 0; i < njobs; ++i) {
        Cn += array[i].processingime * array[i].weight;
    }

    for (i = 0; i < njobs; ++i) {
        C += array[i].processingime;
        C1 += array[i].weight * C;
    }

    lowerbound = (double) C1 / (double) nmachines + (double)(nmachines - 1) /
                 (2.0 * (double)nmachines) * (double)Cn;
    val = (int) ceil(lowerbound);
    return val;
}

int lowerbound_ak(Job *array, int njobs, int nmachines) {
    int i, val = 0;
    int *perm_p = (int *) NULL;
    int *perm_w = (int *) NULL;
    int *start = (int *) NULL;
    perm_p = CC_SAFE_MALLOC(njobs, int);
    perm_w = CC_SAFE_MALLOC(njobs, int);
    start = CC_SAFE_MALLOC(njobs, int);

    for (i = 0; i < njobs; ++i) {
        perm_p[i] = array[i].processingime;
        perm_w[i] = array[i].weight;
    }

    CCutil_int_array_quicksort(perm_p, njobs);
    CCutil_int_array_quicksort_0(perm_w, njobs);

    for (i = 0; i < njobs; ++i) {
        if (i < nmachines) {
            start[i] = 0;
        } else {
            start[i] = start[i - nmachines] + perm_p[i];
        }

        val += perm_w[i] * (start[i]) + array[i].processingime * array[i].weight;
    }

    CC_IFFREE(start, int);
    CC_IFFREE(perm_p, int);
    CC_IFFREE(perm_w, int);
    return val;
}
