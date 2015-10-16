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

/**@file   probdata_coloring.c
 * @brief  problem data for the weighted completion time problem
 * @author Daniel Kowalczyk
 *
 * This file implements the problem data for the weighted completion time algorithm.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "probdata_wct.h"

#define EVENTHDLR_NAME         "probdatavardeleted"
#define EVENTHDLR_DESC         "event handler for variable deleted event"

struct SCIP_ProbData
{
   SCIP_CONS**      constraints;        /* array of added constraints */
   int*             releasetimes;       /* array of release times for every job */
   int*             duetimes;           /* array of due times for every job*/
   int              H_min;              /* minimum completion time of every machine */
   int              H_max;              /* maximum completion time of everey machine */      

   /* schedule set / variable - information*/
   int**            schedulesets;         /* array of schedule sets */
   int*             schedulesetlengths;   /* length of the array in schedulesets */
   int              maxschedulesets;      /* length of array schedulesets */
   int              nschedulesets;        /* number of schedulesets saved in array schedulesets */
   SCIP_VAR**       schedulesetvars;      /* variables belonging to schedule sets */

   /* general information */
   int              njobs;              /* number of jobs*/
   int              nmachines;           /* number of machines */
   int*             processingtimes;    /* processing times of all the jobs */
   int*             weights;            /* weights of all the jobs */
   int              lowerbound;          /* combinatorial lower bound of the current instance */
   Job              *jobs;               /* information for scatter search */
};

/*
 * Local methods
 */

/** 
 *  Preprocessing of 
 */

int max_processingtime(int* processingtimes, int njobs){
   int max = 0;

   for (int i = 0; i < njobs; ++i)
   {
      if(processingtimes[i] > max){
         max = processingtimes[i];
      }
   }

   return max;
}

void calculate_H_max_min(SCIP_PROBDATA *probdata){


   

}


static 
SCIP_RETCODE preprocessSchedule( 
   SCIP*                 scip               /**< SCIP data structure */
   )
{

   SCIP_PROBDATA *probdata;

   int i;
   int j;
   char opt;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   printf("\npreprocessing...\n");

   /* Calculate H_max and H_min */
   int* processingtimes = probdata->processingtimes;
   int* perm = (int*) NULL;
   int nmachine = probdata->nmachines;
   int njobs = probdata->njobs;

   SCIP_CALL(SCIPallocMemoryArray(scip,&(perm),njobs));
   
   double H_max = 0.0;
   double H_min = 0.0;

   
   for (i = 0; i < njobs; ++i)
   {
      H_max += (double) processingtimes[i];
      perm[i] = i;
   }

   H_min = H_max;

   H_max += (double) (nmachine - 1)*max_processingtime(processingtimes, njobs);
   H_max  = H_max/(double)nmachine;
   H_max = ceil(H_max);

   SCIPsortDownIntInt(processingtimes, perm, njobs);

   for (int i = 0; i < nmachine - 1; ++i)
   {
      H_min -= processingtimes[perm[i]];
   }

   H_min = H_min/(double)nmachine;

   probdata->H_max = (int) H_max;
   probdata->H_min = (int) H_min;


   SCIPfreeMemoryArray(scip, perm);
   return SCIP_OKAY;

}




/*
 * Callback methods of probdata
 */

/** transforms the problem */
static
SCIP_DECL_PROBTRANS(probtransWCT)
{
   int i,j;

   assert(scip != NULL);
   assert(sourcedata != NULL);
   assert(targetdata != NULL);

   /* allocate memory */
   SCIP_CALL(SCIPallocMemory(scip, targetdata));

   (*targetdata)->maxschedulesets = sourcedata->maxschedulesets;
   (*targetdata)->nschedulesets = sourcedata->nschedulesets;
   (*targetdata)->njobs = sourcedata->njobs;

   /* Allocate memory to sets and set of lengths */
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->schedulesetlengths), sourcedata->maxschedulesets));   
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->schedulesetvars), sourcedata->maxschedulesets));
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->schedulesets ), sourcedata->maxschedulesets));   

   for (i = 0; i < sourcedata->nschedulesets; ++i)
   {
      assert(sourcedata->schedulesets[i] != NULL);
      (*targetdata)->schedulesetlengths[i] = sourcedata->schedulesetlengths[i];
      SCIP_CALL( SCIPtransformVar(scip, sourcedata->schedulesetvars[i], &((*targetdata)->schedulesetvars[i])));
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->schedulesets[i]),sourcedata->schedulesetlengths[i]));
      for (j = 0; j < sourcedata->schedulesetlengths[i]; ++j)
      {
         (*targetdata)->schedulesets[i][j] = sourcedata->schedulesets[i][j];
      }
   }

   /* allocate memory to jobs dependent array's */
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->duetimes ), sourcedata->njobs));
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->releasetimes ), sourcedata->njobs));
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->processingtimes ), sourcedata->njobs));
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->weights ), sourcedata->njobs));

   for (int i = 0; i < sourcedata->njobs; ++i)
   {
      (*targetdata)->weights[i] = sourcedata->weights[i];
      (*targetdata)->duetimes[i] = sourcedata->duetimes[i];
      (*targetdata)->releasetimes[i] = sourcedata->releasetimes[i];
      (*targetdata)->processingtimes[i] = sourcedata->processingtimes[i];
   }

   /* create array for constraints */
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->constraints), sourcedata->njobs));
   /* transform constraints */
   SCIP_CALL( SCIPtransformConss(scip, sourcedata->njobs, sourcedata->constraints, (*targetdata)->constraints));

   return SCIP_OKAY;
}


/** deletes the transformed problem */
static
SCIP_DECL_PROBDELTRANS(probdeltransWCT)
{
   int i;

   assert(scip != NULL);
   assert(probdata != NULL);

   /*  release the constraints and free memory of the constraints */
   for(i = 0; i < WCTprobGetNjobs(scip); ++i) {
      SCIP_CALL(SCIPreleaseCons(scip, &((*probdata)->constraints[i])) );
   }
   SCIPfreeMemoryArray(scip, &((*probdata)->constraints));

   /* free the arrays of schedulesets and release the associated variables */
   for (i = (*probdata)->nschedulesets - 1; i >= 0; i--)
   {
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->schedulesets[i]), (*probdata)->schedulesetlengths[i]);
      SCIP_CALL(SCIPreleaseVar(scip, &((*probdata)->schedulesetvars[i])));
   }

   SCIPfreeMemoryArray(scip, &((*probdata)->releasetimes));
   SCIPfreeMemoryArray(scip, &((*probdata)->duetimes));
   SCIPfreeMemoryArray(scip, &((*probdata)->weights));
   SCIPfreeMemoryArray(scip, &((*probdata)->processingtimes));
   SCIPfreeBufferArray(scip, &((*probdata)->schedulesetlengths));


   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}

static
SCIP_DECL_PROBDELORIG(probdelorigWCT)
{
   int i;

   assert(scip != NULL);
   assert(probdata != NULL);

   /*  release the constraints and free memory of the constraints */
   for(i = 0; i < WCTprobGetNjobs(scip); ++i) {
      SCIP_CALL(SCIPreleaseCons(scip, &((*probdata)->constraints[i])) );
   }
   SCIPfreeMemoryArray(scip, &((*probdata)->constraints));

   /* free the arrays of schedulesets and release the associated variables */
   for (i = (*probdata)->nschedulesets - 1; i >= 0; i--)
   {
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->schedulesets[i]), (*probdata)->schedulesetlengths[i]);
      SCIP_CALL(SCIPreleaseVar(scip, &((*probdata)->schedulesetvars[i])));
   }

   SCIPfreeMemoryArray(scip, &((*probdata)->releasetimes));
   SCIPfreeMemoryArray(scip, &((*probdata)->duetimes));
   SCIPfreeMemoryArray(scip, &((*probdata)->weights));
   SCIPfreeMemoryArray(scip, &((*probdata)->processingtimes));
   SCIPfreeMemoryArray(scip, &((*probdata)->schedulesetlengths));
   
   /* jobs array */
   SCIPfreeMemoryArray(scip, &((*probdata)->jobs));


   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}


/*
 * Callback methods of event handler
 */

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecProbdatavardeleted)
{
   SCIP_VAR* var;
   SCIP_PROBDATA* probdata;
   int idx;

   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARDELETED);
   var = SCIPeventGetVar(event);
   probdata = (SCIP_PROBDATA*) eventdata;

   assert(probdata != NULL);
   assert(var != NULL);

   /* get index of variable in stablesets array */
   idx = (int)(size_t) SCIPvarGetData(var);

   SCIPdebugMessage("remove variable %s [%d] from list of stable sets\n", SCIPvarGetName(var), idx);


   /* remove variable from stablesets array and release it */


   /* move all subsequent variables to the front */




   return SCIP_OKAY;
}/*lint !e715*/



/*
 * probdata specific interface methods
 */

/** sets up the problem data */
SCIP_RETCODE SCIPcreateProbWCT(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< problem name */           
   int                   njobs,              /**< number of jobs */
   int*                  processingtimes,    /**< processing times for all the jobs */
   int*                  weights,            /**< weights for all the jobs */
   int*                  perm                /**< permutation with respect to the smith rule */
   )
{
   int i;
   SCIP_PROBDATA* probdata = NULL;

   assert(njobs > 0);  /* at least one node */



   printf("Creating problem: %s \n", name);
   
   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, &probdata) );

   /* create constraints */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->constraints), njobs) );
   
   /* original data of jobs and data for scatter search */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->jobs),njobs));
   for (int i = 0; i < njobs; ++i)
   {
      probdata->jobs[i].job = i;
      probdata->jobs[i].processingime = processingtimes[perm[i]];
      probdata->jobs[i].weight = weights[perm[i]];
      probdata->jobs[i].releasetime = 0;
      probdata->jobs[i].duetime = INT16_MAX;
      probdata->jobs[i].after = (GList*) NULL;
      probdata->jobs[i].before = (GList*) NULL;
   }

   /* at the beginning memory for 2 sets */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->schedulesetlengths), 2));
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->schedulesets), 2));
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->schedulesetvars), 2));

   probdata->maxschedulesets = 2;
   probdata->nschedulesets = 0;

   /* include variable deleted event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecProbdatavardeleted, NULL) );

   /* create problem in SCIP */
   SCIP_CALL( SCIPcreateProb(scip, name, probdelorigWCT, probtransWCT, probdeltransWCT, 
         NULL, NULL, NULL, probdata) );

   SCIP_CALL( preprocessSchedule(scip) );

   return SCIP_OKAY;
}


/* ----------------------------------- external methods -------------------------- */

/** returns the number of stable sets / variables */
int WCTprobGetNScheduleSets(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->nschedulesets;
}


/** prints all stable sets to standart output */
void WCTprobPrintScheduleSets(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   int i;
   int j;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   for ( i = 0; i < probdata->nschedulesets; i++ )
   {
      printf( "Set %d: ", i);
      for ( j = 0; j < probdata->schedulesetlengths[i]; j++ )
      {
         printf("%d, ", probdata->schedulesets[i][j]+1);
      }
      printf("ub = %f", SCIPvarGetUbLocal(probdata->schedulesetvars[i]));
      printf(", inLP = %u", SCIPvarIsInLP(probdata->schedulesetvars[i]));
      printf("\n");
   }
}


/** prints the requested stable set to standart output */
void WCTprobPrintScheduleSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setnumber           /**< the number of the requested set */
   )
{
   SCIP_PROBDATA* probdata;
   int i;
   int j;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   i = setnumber;
   printf( "Set %d: ", i);
   for ( j = 0; j < probdata->schedulesetlengths[i]; j++ )
   {
      printf("%d, ", probdata->schedulesets[i][j]+1);
   }
   if ( probdata->schedulesetvars[i] != NULL )
      printf("ub = %f", SCIPvarGetUbLocal(probdata->schedulesetvars[i]));
   printf("\n");
}




/** adds a variable that belongs to a given schedule set */
SCIP_RETCODE WCTprobAddVarForScheduleSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex,           /**< index of the schedule set */
   SCIP_VAR*             var                 /**< pointer to the variable */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   assert((setindex >= 0) && (setindex < probdata->nschedulesets));

   /* catch variable deleted event on the variable to update the schedulesetvars array in the problem data */
   SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_VARDELETED, SCIPfindEventhdlr(scip, EVENTHDLR_NAME),
         (SCIP_EVENTDATA*) probdata, NULL) );

   probdata->schedulesetvars[setindex] = var;

   return SCIP_OKAY;
}


/** gets the variable belonging to a given schedule set */
SCIP_VAR* WCTprobGetVarForScheduleSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex            /**< index of the stable set */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   assert ( (setindex >= 0) && (setindex < probdata->nschedulesets));

   return probdata->schedulesetvars[setindex];
}


/** checks whether a job is in a given schedule set, returns true iff it is */
SCIP_Bool WCTprobIsJobInScheduleSet( 
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex,           /**< index of the stable set */
   int                   node                /**< number of the node */
   )
{
   SCIP_PROBDATA* probdata;
   int l;
   int u;
   int m;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   l = 0;
   u = probdata->schedulesetlengths[setindex]-1;
   while ( l <= u )
   {
      m = (l+u)/2;
      if ( probdata->schedulesets[setindex][m] == node )
      {
         return TRUE;
      }
      if ( probdata->schedulesets[setindex][m] > node )
      {
         l = m+1;
      }
      if ( probdata->schedulesets[setindex][m] < node )
      {
         u = m-1;
      }
   }
   return FALSE;
}


/** checks whether the first set is equal to the second set, both sets have to be sorted in a decreasing way */
SCIP_Bool WCTprobScheduleSetsAreEqual(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  set1,               /**< array of jobs in the first set */ 
   int                   nset1jobs,         /**< number of jobs in the first set */
   int*                  set2,               /**< array of jobs in the second set */ 
   int                   nset2jobs          /**< number of jobs in the second set */
   )
{
   
   int i;

   assert(scip != NULL);
   assert(set1 != NULL && set2 != NULL);
   assert(nset1jobs > 0 && nset2jobs > 0);

   if ( nset1jobs != nset2jobs )
   {
      return FALSE;
   }
   for ( i = 0; i < nset1jobs; i++ )
   {
      if ( set1[i] != set2[i] )
      {
         return FALSE;
      }
   }
   return TRUE;

}


/** checks whether the given schedule set is new
    returns TRUE if the schedule is new, 
            FALSE if it is equal to an already existing schedule set */
SCIP_Bool WCTprobscheduleSetIsNew(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  schedulesetjobs,     /**< array of jobs in the schedule set */
   int                   nschedulesetjobs     /**< number of jobs in the schedule set */
   )
{
   SCIP_PROBDATA* probdata; 
   int i;
   
   assert(schedulesetjobs != NULL);
   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* sort the set */
   SCIPsortDownInt(schedulesetjobs, nschedulesetjobs);

   for ( i = 0; i < WCTprobGetNScheduleSets(scip); i++ )
   {
      if ( WCTprobScheduleSetsAreEqual(scip, schedulesetjobs, nschedulesetjobs, 
               probdata->schedulesets[i], 
               probdata->schedulesetlengths[i]) )
      {
         return FALSE;
      }
   }

   return TRUE;
}

/** adds a new stable set, the set must be sorted descendingly, 
 *  attention: you need to check whether it is new before adding it
 */
SCIP_RETCODE WCTprobAddNewScheduleSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  schedulesetjobs,     /**< array of jobs in the stable set */
   int                   nschedulesetjobs,    /**< number of jobs in the stable set */
   int*                  setindex            /**< return value: index of the stable set */
   )
{
   SCIP_PROBDATA* probdata; 
   int newsize;
   int i;
   
   assert(schedulesetjobs != NULL);
   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* the set should be sorted descendingly */
#ifndef NDEBUG
   for ( i = 0; i < nschedulesetjobs-2; i++ )
   {
      assert(schedulesetjobs[i]>schedulesetjobs[i+1]);
   }
#endif
 
   /* ensure that array is big enough */
   if ( (probdata->nschedulesets + 1) > probdata->maxschedulesets)
   {
      newsize = 2* probdata->maxschedulesets;
      assert(newsize >  probdata->nschedulesets + 1);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(probdata->schedulesets), newsize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(probdata->schedulesetlengths), newsize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(probdata->schedulesetvars), newsize) );
      probdata->maxschedulesets = newsize;
      SCIPdebugMessage("Set-array resized: %d --> %d\n", newsize/2, newsize);
   }

   /* alloc memory for the new stable set */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->schedulesets[probdata->nschedulesets]), nschedulesetjobs) ); /*lint !e866*/
   probdata->schedulesetlengths[probdata->nschedulesets] = nschedulesetjobs;
   probdata->schedulesetvars[probdata->nschedulesets] = NULL;
   for ( i = 0; i < nschedulesetjobs; i++ )
   {
      assert(schedulesetjobs[i] >= 0);
      probdata->schedulesets[probdata->nschedulesets][i] = schedulesetjobs[i];
   }
   *setindex = probdata->nschedulesets;

   probdata->nschedulesets++;

   return SCIP_OKAY;
}


/** returns the schedule set with the given index */
void WCTprobGetScheduleSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex,            /**< index of the stable set */
   int**                 scheduleset,          /**< return value: pointer to the stable set */
   int*                  nelements           /**< return value: number of elements in the stable set */
   )
{
   SCIP_PROBDATA* probdata; 

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   *scheduleset = probdata->schedulesets[setindex];
   *nelements = probdata->schedulesetlengths[setindex];
}


/** returns all stable sets */
void WCTprobGetScheduleSets(
   SCIP*                 scip,               /**< SCIP data structure */
   int***                schedulesets,         /**< return value: pointer to the schedule sets */
   int**                 nelements,          /**< return value: number of elements in the schedule sets */
   int*                  nschedulesets         /**< return value: number of sets */
   )
{
   SCIP_PROBDATA* probdata; 

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   *schedulesets = probdata->schedulesets;
   *nelements = probdata->schedulesetlengths;
   *nschedulesets = probdata->nschedulesets;
}


/** returns the number of jobs */
int WCTprobGetNjobs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->njobs;
}


/** returns the number of machines */
int WCTprobGetNmachines(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->nmachines;
}

/** returns the lowerbound */
int WCTprobGetLB(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->lowerbound;
}

Job* WCTprobGetJoblist(
   SCIP*                scip
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   assert(probdata->jobs != NULL);

   return probdata->jobs;
}         

/** returns all node-constraints */
SCIP_CONS** WCTprobGetConstraints(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->constraints;
}


/** returns the node-constraint belonging to a given node */
SCIP_CONS* WCTprobGetConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   jobs                /**< number of the jobs, for which this constraint assures coloring */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   assert(jobs >= 0 && jobs < probdata->njobs);

   return probdata->constraints[jobs];
}


/** creates all node-constraints and saves them in an array */
SCIP_RETCODE WCTprobSetUpArrayOfCons(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONS** constraints;
   int njobs;
   int i; 

   assert(scip != NULL);

   constraints = WCTprobGetConstraints(scip);
   assert(constraints != NULL);
   njobs = WCTprobGetNjobs(scip);
   for ( i = 0; i < njobs; i++ )
   {
      char consname[SCIP_MAXSTRLEN];
     
      /* create the constraint */
      sprintf(consname, "Node-Constraint%d", i+1);

      SCIP_CALL( SCIPcreateConsSetcover(scip, &constraints[i], consname, 0, NULL, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, constraints[i]) );
   }

   return SCIP_OKAY;
}


/** checks whether the given node is in the given array */
SCIP_Bool WCTprobIsNodeInArray(
   int                   job,               /**< the number of the job */
   int*                  arrayjobs,         /**< the jobs of the maximum scheduleset */
   int                   narrayjobs         /**< number of jobs in the maximum scheduleset */
   )
{
   int i;

   assert(arrayjobs != NULL);

   for ( i = 0; i < narrayjobs; i++ )
   {
      if ( arrayjobs[i] == job )
      {
         return TRUE;
      }
   }
   return FALSE;
}

/** checks whether the given two given arrays are equal */
SCIP_Bool WCTprobEqualSortedArrays(
   int*                  schedule1jobs,         /**< the jobs of the first set */
   int                   nschedule1jobs,        /**< number of jobs in the first set */
   int*                  schedule2jobs,         /**< the jobs of the second set */
   int                   nschedule2jobs         /**< number of jobs in the second set */
   )
{
   int i;

   assert(schedule1jobs != NULL);
   assert(nschedule1jobs > 0);
   assert(schedule2jobs != NULL);
   assert(nschedule2jobs > 0);

   if ( nschedule1jobs != nschedule2jobs )
   {
      return FALSE;
   }
   for ( i = 0; i < nschedule1jobs; i++ )
   {
      if ( schedule1jobs[i] != schedule2jobs[i] )
      {
         return FALSE;
      }
   }
   return TRUE;
}
