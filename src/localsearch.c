#include <assert.h>
#include "wct.h"

int move(Job *j, partlist *m_j, partlist *m_i);

int sort_jobs(gconstpointer a, gconstpointer b) {
    int aa = (((const Job *)a)->job);
    int bb = (((const Job *)b)->job);
    return aa - bb;
}

int move(Job *j, partlist *m_j, partlist *m_i) {
    int nb_job = j->job;
    return j->processingime * (m_j->sumweights[nb_job] - m_i->sumweights[nb_job])
           + j->weight * (m_j->sumtimes[nb_job] - m_i->sumtimes[nb_job] -
                          j->processingime);
}

int k_l_move_general(Job **K_jobs, Job **L_jobs, partlist *m_k, partlist *m_l,
                     solution *sol, int k, int l) {
    int i, val = 0;
    Job *ptr_job = (Job *) NULL;
    Job *ptr_job2 = (Job *) NULL;
    GList *list = (GList *) NULL;
    GList *it = (GList *) NULL;
    GList *it2 = (GList *) NULL;

    if (l > 0) {
        m_l = sol->vlist[L_jobs[0]->job].part;
    }

    for (i = 0;  i < k; i++) {
        list = g_list_insert_sorted(list, K_jobs[i], sort_jobs);
    }

    if (l > 0) {
        for (i = 0; i < l; i++) {
            list = g_list_insert_sorted(list, L_jobs[i], sort_jobs);
        }
    }

    for (it = list; it; it = it->next) {
        ptr_job = (Job *) it->data;

        if (l > 0) {
            if (sol->vlist[ptr_job->job].part == m_k) {
                val += move(ptr_job, m_k, m_l);
            } else {
                val += move(ptr_job, m_l, m_k);
            }

            it2 = list;

            while (it2 != it) {
                ptr_job2 = (Job *)it2->data;

                if (sol->vlist[ptr_job->job].part != sol->vlist[ptr_job2->job].part) {
                    val += 2 * ptr_job->weight * ptr_job2->processingime;
                } else {
                    val -= 2 * ptr_job->weight * ptr_job2->processingime;
                }

                it2 = it2->next;
            }
        } else {
            val += move(ptr_job, m_k, m_l);
            it2 = list;

            while (it2 != it) {
                ptr_job2 = (Job *)it2->data;

                if (sol->vlist[ptr_job->job].part != sol->vlist[ptr_job2->job].part) {
                    val += 2 * ptr_job->weight * ptr_job2->processingime;
                } else {
                    val -= 2 * ptr_job->weight * ptr_job2->processingime;
                }

                it2 = it2->next;
            }
        }
    }

    g_list_free(list);
    return val;
}

int local_search_machine_general_best(solution *sol, int lowerbound, int k,
                                      int l) {
    int i, j, n1, n2, it1, it2, K_flag, L_flag, moved, val = 0;
    int nbiter = 0;
    int njobs = sol->njobs;
    int nmachines = sol->nmachines;
    partlist *L_max = (partlist *) NULL;
    partlist *K_max = (partlist *) NULL;
    CCutil_timer time_move;
    CCutil_init_timer(&time_move, NULL);
    int max;
    int improvement = 0;
    Job **K_jobs = (Job **) NULL;
    Job **L_jobs = (Job **) NULL;
    Job **K_jobs_max = (Job **) NULL;
    Job **L_jobs_max = (Job **) NULL;
    int *K_subset = (int *) NULL;
    int *L_subset = (int *) NULL;
    partlist *temp_m1, *temp_m2;
    K_subset = CC_SAFE_MALLOC(k + 1, int);
    K_jobs = CC_SAFE_MALLOC(k, Job *);
    K_jobs_max = CC_SAFE_MALLOC(k, Job *);

    if (l > 0) {
        L_jobs = CC_SAFE_MALLOC(l, Job *);
        L_subset = CC_SAFE_MALLOC(l + 1, int);
        L_jobs_max = CC_SAFE_MALLOC(l, Job *);
    }

    if (dbg_lvl() > 0) {
        printf("Executing local search with k = %d and l = %d\n", k, l);
    }

    moved = 1;
    CCutil_start_timer(&time_move);

    while (moved && sol->totalweightcomptime != lowerbound) {
        nbiter++;
        moved = 0;
        max = 0;
        int temp;

        for (i = 0; i < nmachines - 1; ++i) {
            for (j = i + 1; j < nmachines; ++j) {
                temp_m1 = sol->part + i;
                temp_m2 = sol->part + j;
                n1 = g_queue_get_length(temp_m1->list);

                if (k > n1) {
                    continue;
                }

                k_subset_init(n1, k, K_subset, &K_flag);

                while (K_flag) {
                    for (it1 = 0; it1 < k; ++it1) {
                        K_jobs[it1] = ((Job *)g_queue_peek_nth(temp_m1->list, K_subset[it1 + 1] - 1));
                    }

                    if (l == 0) {
                        temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                        if (temp > max) {
                            for (it1 = 0; it1 < k; ++it1) {
                                K_jobs_max[it1] = K_jobs[it1];
                            }

                            L_max = temp_m2;
                            moved = 1;
                            max = temp;
                        }
                    } else if (l > 0) {
                        n2 = g_queue_get_length(temp_m2->list);

                        if (l > n2) {
                            continue;
                        }

                        k_subset_init(n2, l, L_subset, &L_flag);

                        while (L_flag) {
                            for (it2 = 0; it2 < l; ++it2) {
                                L_jobs[it2] = ((Job *)g_queue_peek_nth(temp_m2->list, L_subset[it2 + 1] - 1));
                            }

                            temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                            if (temp > max) {
                                for (it1 = 0; it1 < k; ++it1) {
                                    K_jobs_max[it1] = K_jobs[it1];
                                }

                                for (it1 = 0; it1 < l; ++it1) {
                                    L_jobs_max[it1] = L_jobs[it1];
                                }

                                L_max = temp_m2;
                                K_max = temp_m1;
                                moved = 1;
                                max = temp;
                            }

                            k_subset_lex_successor(n2, l, L_subset, &L_flag);
                        }
                    }

                    k_subset_lex_successor(n1, k, K_subset, &K_flag);
                }

                temp_m1 = sol->part + j;
                temp_m2 = sol->part + i;
                n1 = g_queue_get_length(temp_m1->list);

                if (k > n1) {
                    continue;
                }

                k_subset_init(n1, k, K_subset, &K_flag);

                while (K_flag) {
                    for (it1 = 0; it1 < k; ++it1) {
                        K_jobs[it1] = ((Job *)g_queue_peek_nth(temp_m1->list, K_subset[it1 + 1] - 1));
                    }

                    if (l == 0) {
                        temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                        if (temp > max) {
                            for (it1 = 0; it1 < k; ++it1) {
                                K_jobs_max[it1] = K_jobs[it1];
                            }

                            L_max = temp_m2;
                            moved = 1;
                            max = temp;
                        }
                    } else if (l > 0) {
                        n2 = g_queue_get_length(temp_m2->list);

                        if (l > n2) {
                            continue;
                        }

                        k_subset_init(n2, l, L_subset, &L_flag);

                        while (L_flag) {
                            for (it2 = 0; it2 < l; ++it2) {
                                L_jobs[it2] = ((Job *)g_queue_peek_nth(temp_m2->list, L_subset[it2 + 1] - 1));
                            }

                            temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                            if (temp > max) {
                                for (it1 = 0; it1 < k; ++it1) {
                                    K_jobs_max[it1] = K_jobs[it1];
                                }

                                for (it1 = 0; it1 < l; ++it1) {
                                    L_jobs_max[it1] = L_jobs[it1];
                                }

                                L_max = temp_m2;
                                K_max = temp_m1;
                                max = temp;
                                moved = 1;
                            }

                            k_subset_lex_successor(n2, l, L_subset, &L_flag);
                        }
                    }

                    k_subset_lex_successor(n1, k, K_subset, &K_flag);
                }
            }
        }

        if (moved && max > 0) {
            improvement += max;

            if (l == 0) {
                for (it1 = 0; it1 < k; ++it1) {
                    partlist_move_order(L_max, sol->vlist, K_jobs_max[it1], njobs);
                }
            } else if (l > 0) {
                for (it1 = 0; it1 < k; ++it1) {
                    partlist_move_order(L_max, sol->vlist, K_jobs_max[it1], njobs);
                }

                for (it1 = 0; it1 < l; ++it1) {
                    partlist_move_order(K_max, sol->vlist, L_jobs_max[it1], njobs);
                }
            }

            sol->totalweightcomptime -= max;
        }
    };

    CCutil_stop_timer(&time_move, 0);

    if (dbg_lvl()) {
        printf("local search %d  - %d -> number of iterations = %d, objective = %d, time = %f, time per iteration = %f, improvement = %d\n",
               k, l, nbiter, sol->totalweightcomptime, time_move.cum_zeit,
               time_move.cum_zeit / nbiter, improvement);
    }

    CC_IFFREE(L_jobs, Job *);
    CC_IFFREE(K_jobs, Job *);
    CC_IFFREE(K_jobs_max, Job *);
    CC_IFFREE(L_jobs_max, Job *);
    CC_IFFREE(L_subset, int);
    CC_IFFREE(K_subset, int);
    return val;
}

int local_search_machine_general_first(solution *sol, int lowerbound, int k,
                                       int l) {
    int i, j, n1, n2, it1, it2, K_flag, L_flag, moved, val = 0;
    int nbiter = 0;
    int njobs = sol->njobs;
    int nmachines = sol->nmachines;
    partlist *L_max = (partlist *) NULL;
    partlist *K_max = (partlist *) NULL;
    int max;
    CCutil_timer time_move;
    CCutil_init_timer(&time_move, NULL);
    int improvement = 0;
    Job **K_jobs = (Job **) NULL;
    Job **L_jobs = (Job **) NULL;
    Job **K_jobs_max = (Job **) NULL;
    Job **L_jobs_max = (Job **) NULL;
    int *K_subset = (int *) NULL;
    int *L_subset = (int *) NULL;
    partlist *temp_m1, *temp_m2;
    K_subset = CC_SAFE_MALLOC(k + 1, int);
    K_jobs = CC_SAFE_MALLOC(k, Job *);
    K_jobs_max = CC_SAFE_MALLOC(k, Job *);

    if (l > 0) {
        L_jobs = CC_SAFE_MALLOC(l, Job *);
        L_subset = CC_SAFE_MALLOC(l + 1, int);
        L_jobs_max = CC_SAFE_MALLOC(l, Job *);
    }

    if (dbg_lvl()) {
        printf("Executing local search with k = %d and l = %d\n", k, l);
    }

    moved = 1;
    CCutil_start_timer(&time_move);

    while (moved && sol->totalweightcomptime != lowerbound) {
        nbiter++;
        moved = 0;
        max = 0;
        int temp;

        for (i = 0; i < nmachines - 1 && !moved; ++i) {
            for (j = i + 1; j < nmachines && !moved; ++j) {
                temp_m1 = sol->part + i;
                temp_m2 = sol->part + j;
                n1 = g_queue_get_length(temp_m1->list);

                if (k > n1) {
                    continue;
                }

                k_subset_init(n1, k, K_subset, &K_flag);

                while (K_flag && !moved) {
                    for (it1 = 0; it1 < k; ++it1) {
                        K_jobs[it1] = ((Job *)g_queue_peek_nth(temp_m1->list, K_subset[it1 + 1] - 1));
                    }

                    if (l == 0) {
                        temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                        if (temp > max) {
                            for (it1 = 0; it1 < k; ++it1) {
                                K_jobs_max[it1] = K_jobs[it1];
                            }

                            L_max = temp_m2;
                            moved = 1;
                            max = temp;
                        }
                    } else if (l > 0) {
                        n2 = g_queue_get_length(temp_m2->list);
                        k_subset_init(n2, l, L_subset, &L_flag);

                        if (l > n2) {
                            K_flag = 0;
                            continue;
                        };

                        while (L_flag && !moved) {
                            for (it2 = 0; it2 < l; ++it2) {
                                L_jobs[it2] = ((Job *)g_queue_peek_nth(temp_m2->list, L_subset[it2 + 1] - 1));
                            }

                            temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                            if (temp > max) {
                                for (it1 = 0; it1 < k; ++it1) {
                                    K_jobs_max[it1] = K_jobs[it1];
                                }

                                for (it1 = 0; it1 < l; ++it1) {
                                    L_jobs_max[it1] = L_jobs[it1];
                                }

                                L_max = temp_m2;
                                K_max = temp_m1;
                                moved = 1;
                                max = temp;
                            }

                            k_subset_lex_successor(n2, l, L_subset, &L_flag);
                        }
                    }

                    k_subset_lex_successor(n1, k, K_subset, &K_flag);
                }

                if (k != l) {
                    temp_m1 = sol->part + j;
                    temp_m2 = sol->part + i;
                    n1 = g_queue_get_length(temp_m1->list);

                    if (k > n1) {
                        continue;
                    }

                    k_subset_init(n1, k, K_subset, &K_flag);

                    while (K_flag && !moved) {
                        for (it1 = 0; it1 < k; ++it1) {
                            K_jobs[it1] = ((Job *)g_queue_peek_nth(temp_m1->list, K_subset[it1 + 1] - 1));
                        }

                        if (l == 0) {
                            temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                            if (temp > max) {
                                for (it1 = 0; it1 < k; ++it1) {
                                    K_jobs_max[it1] = K_jobs[it1];
                                }

                                L_max = temp_m2;
                                moved = 1;
                                max = temp;
                            }
                        } else if (l > 0) {
                            n2 = g_queue_get_length(temp_m2->list);

                            if (l > n2) {
                                K_flag = 0;
                                continue;
                            };

                            k_subset_init(n2, l, L_subset, &L_flag);

                            while (L_flag && !moved) {
                                for (it2 = 0; it2 < l; ++it2) {
                                    L_jobs[it2] = ((Job *)g_queue_peek_nth(temp_m2->list, L_subset[it2 + 1] - 1));
                                }

                                temp = k_l_move_general(K_jobs, L_jobs, temp_m1, temp_m2, sol, k, l);

                                if (temp > max) {
                                    for (it1 = 0; it1 < k; ++it1) {
                                        K_jobs_max[it1] = K_jobs[it1];
                                    }

                                    for (it1 = 0; it1 < l; ++it1) {
                                        L_jobs_max[it1] = L_jobs[it1];
                                    }

                                    L_max = temp_m2;
                                    K_max = temp_m1;
                                    max = temp;
                                    moved = 1;
                                }

                                k_subset_lex_successor(n2, l, L_subset, &L_flag);
                            }
                        }

                        k_subset_lex_successor(n1, k, K_subset, &K_flag);
                    }
                }
            }
        }

        if (moved && max > 0) {
            improvement += max;

            for (it1 = 0; it1 < k; ++it1) {
                partlist_move_order(L_max, sol->vlist, K_jobs_max[it1], njobs);
            }

            if (l > 0) {
                for (it1 = 0; it1 < l; ++it1) {
                    partlist_move_order(K_max, sol->vlist, L_jobs_max[it1], njobs);
                }
            }

            sol->totalweightcomptime -= max;
        }
    };

    CCutil_stop_timer(&time_move, 0);

    if (dbg_lvl()) {
        printf("local search %d  - %d -> number of iterations = %d, objective = %d, time = %f, time per iteration = %f, improvement = %d\n",
               k, l, nbiter, sol->totalweightcomptime, time_move.cum_zeit,
               time_move.cum_zeit / nbiter, improvement);
    }

    CC_IFFREE(L_jobs, Job *);
    CC_IFFREE(K_jobs, Job *);
    CC_IFFREE(K_jobs_max, Job *);
    CC_IFFREE(L_jobs_max, Job *);
    CC_IFFREE(L_subset, int);
    CC_IFFREE(K_subset, int);
    return val;
}

void localsearch_wrap(solution *sol, int lowerbound, int best) {
    if (best) {
        local_search_machine_general_best(sol, lowerbound, 1, 0);
        local_search_machine_general_best(sol, lowerbound, 1, 1);
        local_search_machine_general_best(sol, lowerbound, 2, 0);
        local_search_machine_general_best(sol, lowerbound, 2, 1);
    } else {
        local_search_machine_general_first(sol, lowerbound, 1, 0);
        local_search_machine_general_first(sol, lowerbound, 1, 1);
        local_search_machine_general_first(sol, lowerbound, 2, 0);
        local_search_machine_general_first(sol, lowerbound, 2, 1);
    }
}

void localsearch_random_k(solution *sol, int lowerbound, int nb) {
    int i, j, l, k, tot;
    int **matrix = (int **) NULL;
    matrix = CC_SAFE_MALLOC(nb, int *);

    for (i = 0; i < nb; ++i) {
        matrix[i] = CC_SAFE_MALLOC(nb, int);

        for (j = 0; j < nb; ++j) {
            matrix[i][j] = 0;
        }
    }

    tot = ((2 + nb) * (nb - 1)) / 2 - 1;
    k = g_random_int_range(1, nb);
    l = g_random_int_range(0, k + 1);
    local_search_machine_general_first(sol, lowerbound, k, l);
    matrix[k][l] = 1;

    for (i = 0; i < tot && !(sol->totalweightcomptime == lowerbound); ++i) {
        while (matrix[k][l] != 0) {
            k = g_random_int_range(1, nb);
            l = g_random_int_range(0, k + 1);
        }

        local_search_machine_general_first(sol, lowerbound, k, l);
        matrix[k][l] = 1;
    }

    for (i = 0; i < nb; ++i) {
        CC_IFFREE(matrix[i], int);
    }

    CC_IFFREE(matrix, int *);
}
