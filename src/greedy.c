#include "wct.h"

/**
 * greedy constructions
 */

int random_rcl_assignment(Job *jobarray, int njobs, int nmachines, solution* new_sol, GRand *rand_) {
  int i;
  int machine_max;
  int job_max;
  int machine_min;
  int job_min;
  double max;
  double min;
  double lb, ub;
  double temp_dbl;

  partlist *temp = (partlist *) NULL;
  Job* temp_job = (Job *) NULL;
  GList *it = (GList *) NULL;
  GQueue *to_do_list = (GQueue *) NULL;

  temp = new_sol->part;
  to_do_list = g_queue_new();

  for (i = 0; i < njobs; ++i) {
    g_queue_push_head(to_do_list, jobarray + i);
  }

  while (!g_queue_is_empty(to_do_list)) {
    temp_job = (Job *)to_do_list->head->data;
    machine_max = 0;
    job_max = (temp_job)->job;
    max = ((double)temp[0].completiontime + (double)temp_job->processingime) / ((double) temp_job->weight);
    machine_min = 0;
    job_min = temp_job->job;
    min = max;
    GArray *rcl = g_array_new(FALSE, FALSE, sizeof(pair_job_machine));
    /** Compute min and max */
    for (i = 0; i < nmachines; ++i)
    {
      for (it = to_do_list->head; it; it = it->next)
      {
        temp_job = (Job*)it->data;
        temp_dbl = ((double)temp[i].completiontime + (double)temp_job->processingime) / ((double) temp_job->weight);
        if (max < temp_dbl)
        {
          machine_max = i;
          job_max = temp_job->job;
          max = temp_dbl;
        }
        if (min > temp_dbl)
        {
          machine_min = i;
          job_min = temp_job->job;
          min = temp_dbl;
        }
      }
    }

    /** Compute RCL */
    pair_job_machine temp_job_machine;
    lb = min;
    ub = min + 0.05 * (max - lb);
    for (i = 0; i < nmachines; ++i)
    {
      for (it = to_do_list->head; it; it = g_list_next(it))
      {
        temp_job = ((Job*)it->data);
        double g = (((double)temp[i].completiontime + (double)temp_job->processingime) / ((double) temp_job->weight));
        if (lb <= g && g <= ub) {
          temp_job_machine.job = temp_job->job;
          temp_job_machine.machine = i;
          g_array_append_val(rcl, temp_job_machine);
        }
      }
    }

    /** Choose uniformaly an assignment of a job to a machine */
    int a = g_rand_int_range(rand_, 0, rcl->len);
    int job = g_array_index(rcl, pair_job_machine, a).job;
    int machine = g_array_index(rcl, pair_job_machine, a).machine;
    partlist_insert(temp + machine, new_sol->vlist, jobarray + job);
    g_queue_pop_nth(to_do_list, g_queue_index(to_do_list, jobarray + job));
    g_array_free(rcl, TRUE);
  }

  g_queue_free(to_do_list);
  return 0;
}

int random_assignment( Job *jobarray, int njobs, int nmachines, solution *new_sol, GRand *rand_ )
{
  int i, val = 0;
  double n;
  GQueue *queue     = ( GQueue * ) NULL;
  partlist *temp = ( partlist * ) NULL;
  Job *j = (Job *) NULL;
  
  queue = g_queue_new();

  for ( i = 0; i < nmachines; i++ ) {
    g_queue_push_head( queue, new_sol->part + i );
  }

  for (i = 0; i < njobs; ++i) {
    j = jobarray + i;
    n = g_rand_double_range(rand_, 0.0, 1.0);
    if (n < 0.8) {
      temp = ( partlist * ) g_queue_pop_head(queue);
    } else if ( n >= 0.8 && n < 0.95) {
      temp = (partlist *) g_queue_pop_nth(queue, 1);
    } else {
      temp = (partlist *) g_queue_pop_nth(queue, 2);
    }
    val = partlist_insert( temp, new_sol->vlist, j );
    CCcheck_val_2( val, "Failed in partlist_insert_order" );
    g_queue_insert_sorted(queue, temp, compare_func1, NULL);
  }

CLEAN:
  g_queue_free( queue );
  return val;
}

int construct_wspt(Job* jobarray, int njobs, int  nmachines, solution*  new_sol) {
  int i, val = 0;
  pmcheap *heap = ( pmcheap * ) NULL;
  partlist *temp = (partlist *) NULL;
  Job *j = (Job *) NULL; 
  pmcheap_init( &heap, nmachines );
  CCcheck_NULL_2( heap, "Failed to initialize heap" );

  for ( i = nmachines - 1; i >= 0; i-- ) {
    pmcheap_insert( heap, new_sol->part[i].completiontime, new_sol->part + i );
  }

  for ( i = 0; i < njobs; i++ ) {
    j = jobarray + i;
    temp = ( partlist * ) pmcheap_min( heap );
    val = partlist_insert( temp, new_sol->vlist, j );
    CCcheck_val_2( val, "Failed in partlist_insert_order" );
    pmcheap_insert( heap, temp->completiontime, temp );
  }

CLEAN:
  if (val)
  {
    solution_free(new_sol);
  }
  pmcheap_free(heap);
  return 0;
}

/**
 * comparefunctions
 */

gint compare_func1(const void *a, const void *b, void *user_data) {
  const int *v = &(((const partlist*)a)->completiontime);
  const int *w = &(((const partlist*)b)->completiontime);

  user_data = NULL;

  if (*v != *w) {
    return *v - *w;
  } else {
    if (*v == 0 || *w == 0) {
      return *v - *w;
    }
    const int *vv = &(((Job*)((const partlist *)a)->list->head->data)->job);
    const int *ww = &(((Job*)((const partlist *)b)->list->head->data)->job);
    return *vv - *ww;
  }
}
