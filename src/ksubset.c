#include "util.h"
#include "defs.h"
#include <sys/utsname.h>

int bin_coef(int n, int r) {
    int i, b;

    if ((r < 0) || (n < r)) {
        return (0);
    }

    if ((2 * r) > n) {
        r = n - r;
    }

    b = 1;

    if (r > 0) {
        for (i = 0; i <= (r - 1); i = i + 1) {
            b = (b * (n - i)) / (i + 1);
        }
    }

    return (b);
}

void k_subset_init(int n, int k, int *subset, int *flag) {
    if (k > n) {
        k = n;
    }

    int i;

    for (i = 1; i <= k; i++) {
        subset[i] = i;
    }

    *flag = 1;
}

int k_subset_lex_successor(int n, int k, int *subset, int *flag) {
    int val = 0, i, j;
    int *temp_set = CC_SAFE_MALLOC(k + 1, int);
    CCcheck_NULL_2(temp_set, "Failed to allocate memory to temp_set");
    *flag = 1;

    for (i = 1; i <= k; ++i) {
        temp_set[i] = subset[i];
    }

    i = k;

    while ((i >= 1) && (subset[i] == (n - k + i))) {
        i--;
    }

    if (i == 0) {
        (*flag) = 0;
    } else {
        for (j = i; j <= k; j++) {
            temp_set[j] = subset[i] + 1 + j - i;
        }

        for (j = 1; j <= k; j++) {
            subset[j] = temp_set[j];
        }
    }

CLEAN:
    CC_IFFREE(temp_set, int);
    return val;
}

void k_subset_lex_rank(int *subset, int k, int n, int *r) {
    int i, j, lo, hi;
    (*r) = 0;
    subset[0] = 0;

    for (i = 1; i <= k; i = i + 1) {
        lo = subset[i - 1] + 1;
        hi = subset[i] - 1;

        if (lo <=  hi) {
            for (j = lo; j <= hi; j = j + 1) {
                (*r) = (*r) + bin_coef(n - j, k - i);
            }
        }
    }
}

void k_subset_lex_unrank(int r, int *T, int n, int k) {
    int x, i, y;
    x = 1;

    for (i = 1; i <= k; i = i + 1) {
        y = bin_coef(n - x, k - i);

        while (y <= r) {
            r = r - y;
            x = x + 1;
            y = bin_coef(n - x, k - i);
        }

        T[i] = x;
        x = x + 1;
    }
}

void print_line() {
    printf("----------------------------------------------------------------------------------\n");
}

int bisearch(int *sorted, const void *target, int size, int
             (*compare)(const void *key1, const void *key2)) {
    int left, middle, right;
    left = 0;
    right = size - 1;

    while (left <= right) {
        middle = (left + right) / 2;

        switch (compare((sorted + (middle)), target)) {
        case -1:
            left = middle + 1;
            break;

        case 1:
            right = middle - 1;
            break;

        case 0:
            return middle;
        }
    }

    return -1;
}

int ksubset_init(int n, int k, ksubset_lex *set) {
    int val = 0;
    set->n = n;

    if (k > n) {
        k = set->n;
    }

    set->m = k;
    set->x = CC_SAFE_MALLOC(set->m + (set->m == 0), int);
    CCcheck_NULL_2(set->x, "Failed to allocate memory to data");
    set->j = 1;
    set->x[0] = 0;
CLEAN:

    if (val) {
        CC_IFFREE(set->x, int);
    }

    return val;
}

void ksubset_free(ksubset_lex *set) {
    set->m = 0;
    set->n = 0;
    set->j = 0;
    CC_IFFREE(set->x, int);
}

int ksubset_next(ksubset_lex *set) {
    int j1 = set->j - 1;
    int z1 = set->x[j1] + 1;

    if (z1 < set->n) {
        if (set->j < set->m) {
            set->x[set->j] = set->x[j1] + 1;
            ++set->j;
            return set->j;
        }

        set->x[j1] = z1;
        return set->j;
    } else {
        if (j1 == 0) {
            return 0;
        }

        --set->j;
        set->x[set->j - 1] += 1;
        return set->j;
    }
}

int *ksubset_data(ksubset_lex *set) {
    return set->x;
}

int ksubset_check(ksubset_lex *set) {
    if (set->x[set->j - 1] >= set->n) {
        return 0;
    }

    if (set->j > set->m) {
        return 0;
    }

    return 1;
}


int ksubset_rec_init(ksubset_rec *set, ulong n) {
    int val = 0;
    set->n = n;
    set->rv = CC_SAFE_MALLOC(set->n + 1, ulong);
    CCcheck_NULL_2(set->rv, "Failed to allocate memory to rv");
    ++set->rv;
    set->rv[-1] = -1UL;
CLEAN:

    if (val) {
        CC_IFFREE(set->rv, ulong);
    }

    return val;
}

void ksubset_rec_free(ksubset_rec *set) {
    --set->rv;
    CC_IFFREE(set->rv, ulong);
}

void ksubset_rec_generate(void *data, ksubset_rec *set, ulong kmin, ulong kmax,
                          ulong rq, ulong nq, void (*visit)(const void *, const void *, ulong)) {
    int temp;
    set->ct = 0;
    set->rct = 0;
    set->kmin = kmin;
    set->kmax = kmax;

    if (set->kmin > set->kmax) {
        CC_SWAP(set->kmin, set->kmax, temp);
    }

    if (set->kmax > set->n) {
        set->kmax = set->n;
    }

    if (set->kmin > set->n) {
        set->kmin = set->n;
    }

    set->visit = visit;
    set->rq = rq % 4;
    set->pq = (rq >> 2) % 4;
    set->nq = nq;
    ksubset_next_rec(data, set, 0);
}

void ksubset_next_rec(void *data, ksubset_rec *set, ulong d) {
    if (d > set->kmax) {
        return;
    }

    ++set->rct;
    ulong rv1 = set->rv[d - 1];
    int q = 1;

    switch (set->rq % 4) {
    case 0:
        q = 1;
        break;

    case 1:
        q = !(d & 1);
        break;

    case 2:
        q = rv1 & 1;

    case 3:
        q = (d ^ rv1) & 1;
        break;
    }

    if (set->nq) {
        q = !q;
    }

    ulong x0 = rv1 + 1;
    ulong rx = set->nq - (set->kmin - d);
    ulong x1 = CC_MIN(set->nq - 1, rx);
#define PCOND(x) if ( (set->pq==x) && (d>=set->kmin) )  { set->visit(data,set, d);  ++set->ct; }
    PCOND(0);

    if (q) {  // forward:
        PCOND(1);

        for (ulong x = x0; x <= x1; ++x)  {
            set->rv[d] = x;
            ksubset_next_rec(data, set, d + 1);
        }

        PCOND(2);
    } else { // backward:
        PCOND(2);

        for (ulong x = x1; x >= x0; --x)  {
            set->rv[d] = x;
            ksubset_next_rec(data, set, d + 1);
        }

        PCOND(1);
    }

    PCOND(3);
#undef PCOND
}

void dump_uname() {
    struct utsname  uts;
    uname(&uts);
    printf("sysname: %s\n", uts.sysname);
    printf("nodename: %s\n", uts.nodename);
    printf("release: %s\n", uts.release);
    printf("version: %s\n", uts.version);
    printf("machine: %s\n", uts.machine);
}

int program_header(int ac, char **av) {
    int val = 0;
    time_t starttime;
    int i;
    printf("########################################################################################################\n");
    printf("Running: ");

    for (i = 0; i < ac; ++i) {
        printf(" %s", av[i]);
    }

    printf("\n");
    dump_uname();
    (void) time(&starttime);
    printf("########################################################################################################\n");
    return val;
}
