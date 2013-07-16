/**
 * @file MICQuery.h
 * @brief Routines to check machine for an Intel Many-Integrated Core (MIC)
 *        board.
 * @author June 18, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef MICQUERY_H_
#define MICQUERY_H_


#ifdef SWIG
#define MIC_ATTRIBUTE
#else
#define MIC_ATTRIBUTE __attribute__((target(mic)))
#endif


#ifdef __cplusplus
#include <offload.h>
#include "../../log.h"
#endif

bool machineContainsMIC();
int machineContainsHowManyMICs();
MIC_ATTRIBUTE bool amIRunningOnMIC();


#endif /* MICQUERY_H_ */
