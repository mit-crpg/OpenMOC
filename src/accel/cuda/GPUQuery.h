/**
 * @file GPUQuery.h
 * @brief Routines to check machine for an NVIDIA GPU and print GPU 
 *        characteristics to the screen.
 * @author May 30, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef GPUQUERY_H_
#define GPUQUERY_H_

#ifdef __cplusplus
#include "../../log.h"
#endif

bool machineContainsGPU();
void attachGPU(int id=0);
void printBasicGPUInfo();
void printDetailedGPUInfo();
int getNumThreadsInWarp();


#endif /* GPUQUERY_H_ */
