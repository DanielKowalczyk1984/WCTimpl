/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_init.c
 * @brief  initial primal heuristic for the weighted completion time problem based on scatter search paradigm
 * @author Daniel Kowalczyk
 * 
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "alloc.h"
#include "heur_init.h"
#include "pricer_wct.h"
#include "reader_wct.h"
#include "scip/cons_setppc.h"
#include "scatter.h"

#define HEUR_NAME             "initschedule"
#define HEUR_DESC             "initial scatter search for WCT"
#define HEUR_DISPCHAR         't'
#define HEUR_PRIORITY         1
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         0
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


/* default values for parameters for scatter search */
#define DEFAULT_USESCATTERSEARCH    TRUE
#define DEFAULT_COMBINEMETHOD    0
#define DEFAULT_B1   10
#define DEFAULT_B2   10
#define DEFAULT_ITER     50
#define DEFAULT_POOLSIZE   40
#define DEFAULT_TIMELIMIT 150.0



/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Bool usescattersearch; /* should the scatter search heuristic be used in order to improve the greedy-solution? */
   SS scatter_search;          /* Scatter search data */
   int combinationmethod;      /* which of the combination methods should be used for scatter search: 0 = GPX, 1 = DellAmico et al, 2 = Pathrelinking" */
   int b1;                     /* number of diverse solution in the refset */
   int b2;                     /* number of good solutions in the refset*/
   int poolsize;               /* number of initial solutions in the pool */
   SCIP_Real timelimit;        /* timelimit of the scatter search algorithm */
   int iter;                   /* number of iterations in the scatter search algorithm */
   int dispfreq;               /* frequency for displaying status information, only active with output verbosity level 2 */
};




/*
 * Local methods
 */

/* Creates a new solution */
solution *new_sol_init( int nmachines, int vcount )
{
    int val = 0;
    solution *sol = ( solution * ) NULL;
    sol = CC_SAFE_MALLOC( 1, solution );
    CCcheck_NULL_2( sol, "Failed to allocate memory" )
    solution_init( sol );
    val = solution_alloc( sol, nmachines, vcount );
    CCcheck_val_2( val, "Failed in solution_alloc" );
CLEAN:

    if ( val ) {
        solution_free( sol );
    }

    return sol;
}

/* update best schedule */
void update_bestschedule( SS* scatter_search, solution *new_sol )
{
    if ( new_sol == NULL ) {
        return;
    }

    if ( new_sol->totalweightcomptime < scatter_search->upperbound ) {
        scatter_search->upperbound = new_sol->totalweightcomptime;
       
    }

    if ( scatter_search->upperbound == scatter_search->lowerbound ) {
        scatter_search->status = opt;
    }
}

/* add solution to pool */
static int add_feasible_solution( SS *scatter_search, solution *new_sol)
{
    int val = 0;

    localsearch_random_k( new_sol, scatter_search->lowerbound, 3 );
    solution_unique( new_sol );

    if ( !solution_in_pool( scatter_search, new_sol ) ) {
        add_solution_pool( scatter_search, new_sol );
        update_bestschedule( scatter_search, new_sol );
        scatter_search->p->PSize++;
    } else {
        solution_free( new_sol );
        CC_IFFREE( new_sol, solution );
    }

    return val;
}

int construct_solution(int njobs,int  nmachines,Job* joblist,solution*  new_sol){
   

   return 0;
}

/** computes the initial schedules with a greedy method */
static
SCIP_RETCODE greedyInitialSchedule(
   SCIP*                 scip               /**< SCIP data structure */
   )
{
   

   return SCIP_OKAY;

}


/** runs scatter search heuristic
 */
static
SCIP_RETCODE runScatterSearch(
   SCIP_HEURDATA*        heurdata,           /**< data of the heuristic */
   SCIP_Bool*            success,             /**< pointer to store if something went wrong */
   int njobs,
   int nmachines,
   int lowerbound,
   Job* joblist,
   SCIP_CLOCK *timelimit
   )
{
   int val = 0;
   SS *scatter_search = &(heurdata->scatter_search);
   P *pool = scatter_search->p;
   REFSET *refset = scatter_search->rs;
   int indicator = 0;
   int nosolutions = 0;

   val = SSproblem_definition(scatter_search,
      heurdata->b1, heurdata->b2,
      heurdata->timelimit,heurdata->combinationmethod,
      njobs, nmachines,
      joblist ,lowerbound);



   while(pool->PSize < heurdata->poolsize
      && scatter_search->status != opt
      && nosolutions <= njobs){

      solution* new_sol = (solution*) NULL;
      new_sol = new_sol_init(nmachines, njobs);
      CCcheck_NULL_2(new_sol, "Failed to allocate memory to new_sol");

      /* construct new solution */
      val = construct_solution(njobs, nmachines, joblist, new_sol);
      CCcheck_val_2(val, "Failed in constructing solution");

      /* add solution to pool */
      add_feasible_solution(scatter_search, new_sol);

      if(scatter_search->status != opt && pool->PSize > 5) {
         /* CALCULATE STRONGER LOWERBOUNDS WITH COMBINATORIAL LOWERBOUND TECHNIQUES */
      }

      if ( !( indicator <  scatter_search->p->PSize ) ) {
          nosolutions++;
      } else {
          nosolutions = 0;
      }

   }

   printf("We have found %d solutions in %f seconds with timelimt = %f\n",pool->lenght, 40.0, scatter_search->timelimit );

   if(scatter_search->status != opt) {
      SSCreate_refset(scatter_search);
      SSrun_scatter_search(scatter_search);
   }

   CLEAN:
   if(val) {
      return SCIP_ERROR;
   }
   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyInit)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeInit)
{
   SCIP_HEURDATA* heurdata;

   /* free heuristic rule data */
   heurdata = SCIPheurGetData(heur);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

void WCTgetJobsPartlist(int *jobs, int *nscheduleset, partlist *list){
   GList *it = list->list->head;
   *nscheduleset = 0;

   for (; it != NULL; it = g_list_next(it))
   {
      jobs[(*nscheduleset)++] = ((Job*)it->data)->job;
   }
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecInit)
{
   int i;
   int j;
   int njobs;
   int nmachines;
   int njobsset;
   int *jobs;
   Job *joblist;
   int lowerbound;
   SCIP_SOL* sol;
   SCIP_Bool stored;
   SCIP_Bool success;
   int nschedulejobs;
   int setnumber;
   SCIP_VAR* var;
   SCIP_CONS** constraints;
   SCIP_HEURDATA* heurdata;
   SS* scatter_search;
   REFSET* refset;
   SCIP_CLOCK* timelimit;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);


   njobs = WCTprobGetNjobs(scip);
   nmachines = WCTprobGetNmachines(scip);
   joblist = WCTprobGetJoblist(scip);
   lowerbound = WCTprobGetLB(scip);

   SCIP_CALL(SCIPallocMemoryArray(scip,&jobs,njobs));

   SCIPcreateClock(scip, &timelimit);

   *result = SCIP_DIDNOTFIND;

   /* get the job-constraits */
   constraints = WCTprobGetConstraints(scip);
   assert(constraints != NULL);

   /* compute an initial schedules with a scatter_search algorithm */
   runScatterSearch(heurdata, &success, njobs, nmachines, lowerbound, joblist, timelimit);
   scatter_search = &(heurdata->scatter_search);
   refset = scatter_search->rs;


   /* create vars for the computed schedules */
   while(!g_queue_is_empty(refset->list1)){
      solution *sol = g_queue_pop_head(refset->list1);

      for (i = 0; i < sol->nmachines; ++i)
      {

         WCTgetJobsPartlist(jobs,&nschedulejobs,sol->part + i);
         SCIP_CALL( WCTprobAddNewScheduleSet(scip, jobs, nschedulejobs, &setnumber) );
         assert(setnumber != -1);

         /* create variable for the stable set and add it to SCIP */
         SCIP_CALL( SCIPcreateVar(scip, &var, NULL, 0.0, 1.0, sol->totalweightcomptime, SCIP_VARTYPE_BINARY,
               TRUE, TRUE, NULL, NULL, NULL, NULL, (SCIP_VARDATA*)(size_t)setnumber) ); /*lint !e571*/

         SCIP_CALL( WCTprobAddVarForScheduleSet(scip, setnumber, var) );
         SCIP_CALL( SCIPaddVar(scip, var) );
         SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );

         for( j = 0; j < nschedulejobs; j++ )
         {
            /* add variable to node constraints of nodes in the set */
            SCIP_CALL( SCIPaddCoefSetppc(scip, constraints[jobs[j]], var) );
         }
      }
      solution_free(sol);
   }
  
   SCIPfreeMemoryArray(scip, &jobs);
   
   /* create solution consisting of all yet created stable sets,
      that means all sets of the solution given by the solution file or created by the greedy and tabu search */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   assert(sol != NULL);
   for( i = 0; i < WCTprobGetNScheduleSets(scip); i++ )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, WCTprobGetVarForScheduleSet(scip, i), 1.0) );
   }
   SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, FALSE, FALSE, FALSE, &stored) );
   assert(stored);

   /* set maximal number of variables to be priced in each round */
   SCIP_CALL( SCIPsetIntParam(scip, "pricers/coloring/maxvarsround",
         WCTprobGetNScheduleSets(scip)*WCTprobGetNjobs(scip)/50) );

   *result = SCIP_FOUNDSOL;
   SS_free(scatter_search);
   return SCIP_OKAY;
}/*lint !e715*/

/*
 * primal heuristic specific interface methods
 */

/** creates the init primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurInit(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create init primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   heur = NULL;
   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecInit, heurdata) );
   assert(heur != NULL);

   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyInit) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeInit) );

   /* add parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/initschedule/usescattersearch",
         "should the scatter search heuristic be used in order to improve the greedy-solution?",
         &heurdata->usescattersearch, TRUE, DEFAULT_USESCATTERSEARCH, NULL, NULL) );


   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/initschedule/combinationmethod",
         "which of the combination methods should be used for scatter search: 0 = GPX, 1 = Dell'Amico et al., 2 = Pathrelinking\n",
         &heurdata->combinationmethod, TRUE, DEFAULT_COMBINEMETHOD, 0, 3, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/initschedule/b1",
         "number of diverse solutions in the refset",
         &heurdata->b1, TRUE, DEFAULT_B1, 0, 20, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/initschedule/b2",
         "number of good solutions in the refset",
         &heurdata->b2, TRUE, DEFAULT_B2, 0, 20, NULL, NULL) );


   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/initschedule/poolsize",
         "number of solutions in the initial poolsize",
         &heurdata->poolsize, TRUE, DEFAULT_POOLSIZE, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "/heuristics/initschedule/iter",
         "number of interations in the scatter search",
         &heurdata->iter, TRUE, DEFAULT_ITER, 0, 1000,NULL, NULL));

   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/initschedule/timelimit",
         "timelimit for the scatter search algorithm",
         &heurdata->timelimit, TRUE, DEFAULT_TIMELIMIT, 0.0, 150.0, NULL, NULL));


   return SCIP_OKAY;
}
