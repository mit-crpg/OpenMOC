/**
 * @file MICQuery.h
 * @brief Routines to check machine for an Intel Many-Integrated Core (MIC)
 *        board.
 * @author June 18, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef MICQUERY_H_
#define MICQUERY_H_

#ifdef __cplusplus
#include <offload.h>
#include "../../log.h"
#endif

bool machineContainsMIC();
int machineContainsHowManyMICs();
bool amIRunningOnMIC();


#endif /* MICQUERY_H_ */
