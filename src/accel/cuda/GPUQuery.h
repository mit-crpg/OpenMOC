/**
 * @file GPUQuery.h
 * @brief Routines to check machine for an NVIDIA GPU and print GPU
 *        and CUDA hardware characteristics to the screen.
 * @author May 30, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef GPUQUERY_H_
#define GPUQUERY_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "../../log.h"
#endif

bool machine_contains_gpu();
void attach_gpu(int id=0);
void print_basic_gpu_info();
void print_detailed_gpu_info();
int get_num_threads_per_warp();

#endif /* GPUQUERY_H_ */
