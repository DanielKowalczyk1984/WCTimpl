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

/**@file   reader_col.c
 * @brief  file reader for vertex coloring instances
 * @author Gerald Gamrath
 *
 * This file implements the reader for vertex coloring problems in DIMACS standard format.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "reader_wct.h"



#define READER_NAME             "wctreader"
#define READER_DESC             "file reader for a WCT instances"
#define READER_EXTENSION        "txt"

#define COL_MAX_LINELEN 1024

/*
 * Local methods
 */

/** get next number from string s */
static
long getNextNumber(
  char**                s                   /**< pointer to the pointer of the current position in the string */
)
{
  long tmp;
  /* skip whitespaces */
  while ( isspace(**s) )
    ++(*s);
  /* read number */
  tmp = atol(*s);
  /* skip whitespaces */
  while ( (**s != 0) && (!isspace(**s)) )
    ++(*s);
  return tmp;
}

/** read LP in "COL File Format" */  
static
SCIP_RETCODE readWCT(
   SCIP*                 scip,               /**< SCIP data structure */   
   const char*           filename            /**< name of the input file */
   )
{
   SCIP_FILE* fp;               /* file-reader */
   char buf[COL_MAX_LINELEN];   /* maximal length of line */
   char* char_p;
   char* probname;
   int njobs;
   int nmachines;
   int* weights;
   int* processingtimes;
   int *perm;
   double *ratio;
   int i;
   int j;
   int weight;
   int processtime;

   
   assert(scip != NULL);
   assert(filename != NULL);
   
   if (NULL == (fp = SCIPfopen(filename, "r")))
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      perror(filename);
      return SCIP_NOFILE;
   }
   
   /* Get problem name from filename and save it */
   if( SCIPfgets(buf, (int) sizeof(buf), fp) == NULL)
      return SCIP_READERROR;

   i = 1;
   while ( (filename[i] != '/') && (filename[i] != '\0') )
   {
      i++;
   }
   if ( filename[i] != '/' )
   {
      j = i;
      i = -1;
   }
   else
   {
      j = i+1;
      while ( filename[i] == '/' && filename[j] != '\0' )
      {
         j = i+1;
         while ( filename[j] != '\0' )
         {
            j++;
            if ( filename[j] == '/' )
            {
               i = j;
               break;
            }
         }
      }
   }

   if( j-i-4 <= 0 )
      return SCIP_READERROR;

   SCIP_CALL( SCIPallocMemoryArray(scip, &probname, (j-i-4)) );
   strncpy(probname, &filename[i+1], (j-i-5)); /*lint !e732 !e776*/
   probname[j-i-5]= '\0';

   /* Read until information about graph starts */
   while( !SCIPfeof(fp) && (buf[0] != 'p') )
   {
      SCIPfgets(buf, (int) sizeof(buf), fp); /*lint !e534*/
   } 
   /* no graph information in file! */
   if ( SCIPfeof(fp) )
   {
      SCIPerrorMessage("Error! Could not find line starting with 'p'.\n");
      return SCIP_READERROR;
   }
   /* wrong format of the line containig number of nodes and edges */
   if ( buf[2] != 'j' || buf[3] != 'o' || buf[4] != 'b' || buf[5] != 's' )
   {
      SCIPerrorMessage("Line starting with 'p' must continue with 'jobs'!\n");
      return SCIP_READERROR;
   }
   char_p = &buf[6];

   /* read out number of jobs and machines, the pointer char_p will be changed */
   njobs = (int) getNextNumber(&char_p);
   nmachines = (int) getNextNumber(&char_p);
   if ( njobs <= 0 )
   {
      SCIPerrorMessage("Number of jobs must be positive!\n");
      return SCIP_READERROR;
   }
   if ( nmachines < 0 )
   {	  
      SCIPerrorMessage("Number of machines must be nonnegative!\n");
      return SCIP_READERROR;
   }
   /* create array for weights and processing times */
   SCIP_CALL( SCIPallocMemoryArray(scip, &weights, njobs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &processingtimes, njobs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &perm, njobs));
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ratio, njobs));      
   
   /* fill array for edges */
   i = 0;
   while ( !SCIPfeof(fp) )
   {
      SCIPfgets(buf, (int) sizeof(buf), fp); /*lint !e534*/
      if ( buf[0] == 'n')
      {
         char_p = &buf[2];
         
         weight = (int) getNextNumber(&char_p);
         processtime = (int) getNextNumber(&char_p);
         weights[i] = weight;
         processingtimes[i] = processtime;
         ratio[i] = (double) weight/ (double) processtime;
         perm[i] = i++;
      }
   }

   if (i != njobs)
   {
     SCIPerrorMessage("Number of jobs is not equal to the number of lines that begin with n");
     return SCIP_READERROR;
   }

   printf("Read WCT instance: %d njobs, %d nmachines\n", njobs, nmachines);

   SCIPsortDownRealInt(ratio, perm, njobs);

   /* create problem data */
   SCIP_CALL( SCIPcreateProbWCT(scip, probname, njobs, processingtimes, weights, perm) );

   /* create LP */
   SCIPdebugMessage("CreateLP...\n");
   SCIP_CALL( WCTprobSetUpArrayOfCons(scip) );

   
   /* activate the pricer */
   SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, "coloring")) );
   SCIP_CALL( SCIPsetObjIntegral(scip) );
   
   SCIPfreeMemoryArray(scip, &weights);
   SCIPfreeMemoryArray(scip, &processingtimes);
   SCIPfreeMemoryArray(scip, &perm);
   SCIPfreeMemoryArray(scip, &ratio);
   SCIPfreeMemoryArray(scip, &probname);
   SCIPfclose(fp);

   return SCIP_OKAY;
}




/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyTxt)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
 
   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadTxt)
{  /*lint --e{715}*/
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   
   SCIP_CALL( readWCT(scip, filename) );
   
   *result = SCIP_SUCCESS;
   
   return SCIP_OKAY;
}

/*
 * col file reader specific interface methods
 */

/** includes the col file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderCol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
  SCIP_READERDATA* readerdata;
  SCIP_READER* reader;

  /* create col reader data */
  readerdata = NULL;

  /* include col reader */
  SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION,
        readerdata) );

  SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyTxt) );
  SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadTxt) );


  return SCIP_OKAY;
}
