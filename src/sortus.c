//
//  sortrus.c
//  PMC
//
//  Created by Daniel on 21/02/14.
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved.
//

/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                         SORTING ROUTINES                                 */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*   Written by:  Applegate, Bixby, Chvatal, and Cook                       */
/*   DATE:  February 24, 1994                                               */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  char *CCutil_linked_radixsort (char *data, char *datanext,              */
/*      char *dataval, int valsize)                                         */
/*    USAGE:                                                                */
/*      head = (bar *) CCutil_linked_radixsort ((char *) head,              */
/*         (char *) &(head->next), (char *) &(head->val), sizeof (int));    */
/*    Then head is the start of the linked list in increasing order of      */
/*    val, with next as the field that links the bars.                      */
/*    WARNING: DOES NOT HANDLE NEGATIVE NUMBERS PROPERLY.                   */
/*                                                                          */
/*  void CCutil_int_array_quicksort (int *len, int n)                       */
/*    len - the array to be sorted                                          */
/*    n - the number of elements in len                                     */
/*    Uses quicksort to put len in increasing order.                        */
/*                                                                          */
/*  void CCutil_int_perm_quicksort (int *perm, int *len, int n)             */
/*    perm - must be allocated and initialized by the calling routine,      */
/*           it will be arranged in increasing order of len.                */
/*    n - the number of elements in perm and len.                           */
/*                                                                          */
/*  void CCutil_double_perm_quicksort (int *perm, double *len, int n)       */
/*    perm - must be allocated and initialized by the calling routine,      */
/*           it will be arranged in increasing order of len.                */
/*    n - the number of elements in perm and len.                           */
/****************************************************************************/

#include "defs.h"
#include "util.h"


#define BITS_PER_PASS (8)

#define NBINS (1<<BITS_PER_PASS)




void CCutil_int_array_quicksort(int *len, int n)
{
    int i, j, temp, t;

    if (n <= 1) {
        return;
    }

    CC_SWAP(len[0], len[(n - 1) / 2], temp);
    i = 0;
    j = n;
    t = len[0];

    while (1) {
        do {
            i++;
        } while (i < n && len[i] < t);

        do {
            j--;
        } while (len[j] > t);

        if (j < i) {
            break;
        }

        CC_SWAP(len[i], len[j], temp);
    }

    CC_SWAP(len[0], len[j], temp);
    CCutil_int_array_quicksort(len, j);
    CCutil_int_array_quicksort(len + i, n - i);
}

void CCutil_int_array_quicksort_0(int *len, int n)
{
    int i, j, temp, t;

    if (n <= 1) {
        return;
    }

    CC_SWAP(len[0], len[(n - 1) / 2], temp);
    i = 0;
    j = n;
    t = len[0];

    while (1) {
        do {
            i++;
        } while (i < n && len[i] > t);

        do {
            j--;
        } while (len[j] < t);

        if (j < i) {
            break;
        }

        CC_SWAP(len[i], len[j], temp);
    }

    CC_SWAP(len[0], len[j], temp);
    CCutil_int_array_quicksort_0(len, j);
    CCutil_int_array_quicksort_0(len + i, n - i);
}

void CCutil_int_perm_quicksort(int *perm, int *len, int n)
{
    int i, j, temp, t;

    if (n <= 1) {
        return;
    }

    CC_SWAP(perm[0], perm[(n - 1) / 2], temp);
    i = 0;
    j = n;
    t = len[perm[0]];

    while (1) {
        do {
            i++;
        } while (i < n && len[perm[i]] < t);

        do {
            j--;
        } while (len[perm[j]] > t);

        if (j < i) {
            break;
        }

        CC_SWAP(perm[i], perm[j], temp);
    }

    CC_SWAP(perm[0], perm[j], temp);
    CCutil_int_perm_quicksort(perm, len, j);
    CCutil_int_perm_quicksort(perm + i, len, n - i);
}

void CCutil_int_perm_quicksort_0(int *perm, int *len, int n)
{
    int i, j, temp, t;

    if (n <= 1) {
        return;
    }

    CC_SWAP(perm[0], perm[(n - 1) / 2], temp);
    i = 0;
    j = n;
    t = len[perm[0]];

    while (1) {
        do {
            i++;
        } while (i < n && len[perm[i]] > t);

        do {
            j--;
        } while (len[perm[j]] < t);

        if (j < i) {
            break;
        }

        CC_SWAP(perm[i], perm[j], temp);
    }

    CC_SWAP(perm[0], perm[j], temp);
    CCutil_int_perm_quicksort_0(perm, len, j);
    CCutil_int_perm_quicksort_0(perm + i, len, n - i);
}


void CCutil_double_perm_quicksort(int *perm, double *len, int n)
{
    int i, j, temp;
    double t;

    if (n <= 1) {
        return;
    }

    CC_SWAP(perm[0], perm[(n - 1) / 2], temp);
    i = 0;
    j = n;
    t = len[perm[0]];

    while (1) {
        do {
            i++;
        } while (i < n && len[perm[i]] < t);

        do {
            j--;
        } while (len[perm[j]] > t);

        if (j < i) {
            break;
        }

        CC_SWAP(perm[i], perm[j], temp);
    }

    CC_SWAP(perm[0], perm[j], temp);
    CCutil_double_perm_quicksort(perm, len, j);
    CCutil_double_perm_quicksort(perm + i, len, n - i);
}


static int partition(int *input, int p, int r)
{
    int pivot = input[r];

    while (p < r) {
        while (input[p] < pivot) {
            p++;
        }

        while (input[r] > pivot) {
            r--;
        }

        if (input[p] == input[r]) {
            p++;
        } else if (p < r) {
            int tmp = input[p];
            input[p] = input[r];
            input[r] = tmp;
        }
    }

    return r;
}

int CCutil_quickselect(int *input, int p, int r, int k)
{
    if (p == r) {
        return input[p];
    }

    int j = partition(input, p, r);
    int length = j - p + 1;

    if (length == k) {
        return input[j];
    } else if (k < length) {
        return CCutil_quickselect(input, p, j - 1, k);
    } else {
        return CCutil_quickselect(input, j + 1, r, k - length);
    }
}

static int recursiveSelect(int *V, int *inds, int start, int end, int k)
{
    /*recursively partitions vector V to find kth element */
    int pivot, i, tmp;
    float tmpdbl;

    if (start == end - 1) { /* active set is only one entry    */
        return inds[start];
    }

    /*otherwise pick pivot and split */
    pivot = start;

    for (i = start + 1; i < end; i++) {
        /* put all entries greater than V[pivot] to the left of pivot*/
        if (V[i] > V[pivot]) {
            /* first move pivot value up one index */
            tmp = inds[pivot];
            tmpdbl = V[pivot];
            inds[pivot] = inds[pivot + 1];
            V[pivot] = V[pivot + 1];
            inds[pivot + 1] = tmp;
            V[pivot + 1] = tmpdbl;

            if (i > pivot + 1) {
                /* now we need to swap i th with pivot-1 th */
                tmp = inds[pivot];
                tmpdbl = V[pivot];
                inds[pivot] = inds[i];
                V[pivot] = V[i];
                inds[i] = tmp;
                V[i] = tmpdbl;
            }

            pivot++;
        }
    }

    if (pivot == k) {
        return inds[pivot];
    }
    /* now we decide if kth is to the right of the pivot or the left */
    else if (pivot > k) {
        return recursiveSelect(V, inds, start, pivot, k);
    } else {
        return recursiveSelect(V, inds, pivot + 1, end, k);
    }
}


int quickselect(int *V, int N, int k)
{
    /* returns the index of the kth greatest value in vector V */
    int *Vcopy;
    int i, *inds, ret;

    if (k >= N) {
        return -1;
    }

    if (k < 0) {
        return -1;
    }

    Vcopy = (int *)malloc(N * sizeof(float));
    inds = (int *)malloc(N * sizeof(int));

    for (i = 0; i < N; i++) {
        Vcopy[i] = V[i];
        inds[i] = i;
    }

    ret = recursiveSelect(Vcopy, inds, 0, N, k);
    free(Vcopy);
    free(inds);
    return ret;
}
