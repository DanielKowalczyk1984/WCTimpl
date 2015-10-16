#ifndef _K_SUBSET_H
#define _K_SUBSET_H
#include "alloc.h"

/*Generation of k-subsets*/
typedef struct ksubset_lex{
    int n;
    int j;
    int m;
    int *x;
} ksubset_lex;

int ksubset_init(int n,int k,ksubset_lex *set);
void ksubset_free(ksubset_lex *set);
int ksubset_next(ksubset_lex *set);
int* ksubset_data(ksubset_lex *set);
int ksubset_check(ksubset_lex *set);

typedef  unsigned long ulong;

typedef struct ksubset_rec{
    ulong n;
    ulong kmin,kmax;
    ulong *rv;
    ulong ct;
    ulong rct;
    ulong rq;
    ulong pq;
    ulong nq;
    void (*visit)(const void *,const void*,ulong);
} ksubset_rec;

int ksubset_rec_init(ksubset_rec *set,ulong n);
void ksubset_rec_free(ksubset_rec *set);
void ksubset_rec_generate(void *data,ksubset_rec *set,ulong kmin,ulong kmax,ulong rq,ulong nq,void (*visit)(const void *,const void* ,ulong));
void ksubset_next_rec(void *data,ksubset_rec *set,ulong d);

/*Generation of k-subsets*/
int bin_coef(int n, int r);
void k_subset_init(int n,int k,int *subset,int *flag);
int k_subset_lex_successor(int n,int k,int *subset, int *flag);
void k_subset_lex_rank(int *subset,int k,int n,int * r);
void k_subset_lex_unrank(int r, int *T,int n,int k);

#endif



