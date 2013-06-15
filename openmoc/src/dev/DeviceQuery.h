/**
 * @file DeviceQuery.h
 * @brief Routines to check machine for an NVIDIA GPU and print GPU 
 *        characteristics to the screen.
 * @author May 30, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef DEVICEQUERY_H_
#define DEVICEQUERY_H_

#ifdef __cplusplus
#include "../host/log.h"
#endif

bool machineContainsGPU();
void attachGPU(int id=0);
void printBasicDeviceInfo();
void printDetailedDeviceInfo();


#endif /* DEVICEQUERY_H_ */
