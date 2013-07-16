/**
 * @file DeviceQuery.h
 * @brief Routines to check machine for an Intel Xeon Phi (MIC) and print MIC 
 *        characteristics to the screen.
 * @author June 11, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef DEVICEQUERY_H_
#define DEVICEQUERY_H_

#ifdef __cplusplus
#include <offload.h>
#include "../host/log.h"
#endif

bool machineContainsMIC();
int getNumMICDevices();


#endif /* DEVICEQUERY_H_ */
