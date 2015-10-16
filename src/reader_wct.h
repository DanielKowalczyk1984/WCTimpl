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

/**@file   reader_col.h
 * @brief  file reader for Weighted completion time instances
 * @author Daniel Kowalczyk
 *
 * This file implements the reader for Weighted completion time instances.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_WCT_H__
#define __SCIP_READER_WCT_H__

#include "scip/scip.h"
#include "probdata_wct.h"
#include "scip/cons_setppc.h"
#include "scip/type_cons.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes the col file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderTxt(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
