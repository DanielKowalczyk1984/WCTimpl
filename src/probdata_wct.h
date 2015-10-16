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

/**@file   probdata_wct.h
 * @brief  problem data for vertex coloring algorithm
 * @author Gerald Gamrath
 *
 * This file implements the problem data for the weighted completion time algorithm.
 *
 * The problem data contains
 *
 * The preprocessing 
 *
 *
 * Each variable has a pointer of type SCIP_VARDATA* that is used in this case to store an integer
 * representing the number of the schedule set. With the aid of this, the corresponding schedule set can
 * be found in the array returned by WCTprobGetscheduleSets().  This array contains all schedule sets
 * and is also used to check whether a schedule set found by the pricer is really new. This can be
 * done by calling COLORprobScheduleSetIsNew(). All sets are sorted decreasingly with respect to the
 * indices of the jobs. New candidates should also be sorted that way.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROBDATA_COLORING__
#define __SCIP_PROBDATA_COLORING__

#include <assert.h>

#include "scip/scip.h"
#include "scip/cons_setppc.h"
#include "reader_wct.h"
#include "scatter.h"

#ifdef __cplusplus
extern "C" {
#endif
/* schedule set / variable functions */

/** returns the number of schedule sets / variables */
extern
int WCTprobGetNScheduleSets(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the schedule set with the given index */
extern
void WCTprobGetScheduleSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex,           /**< index of the schedule set */
   int**                 scheduleset,          /**< return value: pointer to the schedule set */
   int*                  nelements           /**< return value: number of elements in the schedule set */
   );

/** returns all schedule sets */
extern
void WCTprobGetScheduleSets(
   SCIP*                 scip,               /**< SCIP data structure */
   int***                schedulesets,         /**< return value: pointer to the schedule sets */
   int**                 nelements,          /**< return value: number of elements in the schedule sets */
   int*                  nschedulesets         /**< return value: number of sets */
   );

/** adds a new schedule set, the set must be sorted descendingly, 
    attention: you need to check whether it is new before adding it*/
extern
SCIP_RETCODE WCTprobAddNewScheduleSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  schedulejobs,        /**< array of jobs in the schedule set */
   int                   nschedulejobs,       /**< number of jobs in the schedule set */
   int*                  setindex            /**< return value: index of the schedule set, -i-1 if set was not new 
                                              *   and is already saved as set i */
   );

/** adds a variable that belongs to a given schedule set */
extern
SCIP_RETCODE WCTprobAddVarForScheduleSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex,              /**< index of the schedule set */
   SCIP_VAR*             var                 /**< pointer to the variable */
   );

/** gets the variable belonging to a given schedule set */
extern
SCIP_VAR* WCTprobGetVarForScheduleSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex            /**< index of the schedule set */
   );

/** checks whether the given schedule set is new, returns TRUE if the stable is new and 
 *  FALSE if it is part of or equal to an already existing schedule set 
 */
extern
SCIP_Bool WCTprobScheduleSetIsNew(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  schedulesetjobs,     /**< array of jobs in the schedule set */
   int                   nschedulesetjobs     /**< number of jobs in the schedule set */
   );

/** checks whether the first set is equal to the second set, both sets have to be sorted in a decreasing way */
extern
SCIP_Bool WCTprobScheduleSetsAreEqual(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  set1,               /**< array of jobs in the first set */ 
   int                   nset1jobs,         /**< number of jobs in the first set */
   int*                  set2,               /**< array of jobs in the second set */ 
   int                   nset2jobs          /**< number of jobs in the second set */
   );

/** prints all schedule sets to standart output */
extern
void WCTprobPrintScheduleSets(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** prints the requested schedule set to standart output */
extern
void WCTprobPrintScheduleSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setnumber           /**< the number of the requested set */
   );

/** checks whether a Job is in a given schedule set, returns true iff it is */
extern
SCIP_Bool WCTprobIsJobInScheduleSet( 
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex,           /**< index of the schedule set */
   int                   job                /**< number of the Job */
   );

/* constraint functions */

/** creates all job-constraints and saves them in an array */
extern
SCIP_RETCODE WCTprobSetUpArrayOfCons(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns all job-constraints */
extern
SCIP_CONS** WCTprobGetConstraints(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the job-constraint belonging to a given Job */
extern
SCIP_CONS* WCTprobGetConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   job                /**< number of the Job, for which this constraint assures WCTing */
   );



/* preprocessing functions */





/** returns the number of jobs */
extern
int WCTprobGetNjobs(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
int WCTprobGetNmachines(
   SCIP*                scip
   );

extern
int WCTprobGetLB(
   SCIP*                scip
   );

extern
Job* WCTprobGetJoblist(
   SCIP*                scip
   ); 
/* miscellaneous functions */

/** checks whether the given Job is in the given array */ 
extern
SCIP_Bool WCTprobIsJobInArray(
   int                   job,               /**< the number of the Job */
   int*                  arrayjobs,         /**< the jobs of the maximum clique */
   int                   narrayjobs         /**< number of jobs in the maximum clique */
   );

/** checks whether the given two given arrays are equal */
extern
SCIP_Bool WCTprobEqualSortedArrays(
   int*                  array1jobs,         /**< the jobs of the first set */
   int                   narray1jobs,        /**< number of jobs in the first set */
   int*                  array2jobs,         /**< the jobs of the second set */
   int                   narray2jobs         /**< number of jobs in the second set */
   );


/* create probdate */

/** sets up the problem data */
extern
SCIP_RETCODE SCIPcreateProbWCT(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< problem name */           
   int                   njobs,              /**< number of jobs */
   int*                  processingtimes,    /**< processing times for all the jobs */
   int*                  weights,             /**< weights for all the jobs */
   int*                  perm                /**< permutation with respect to the smith rule */
   );

#ifdef __cplusplus
}
#endif

#endif
