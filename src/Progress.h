/**
 * @file Progress.h
 * @brief An object to track progress
 * @date January 11, 2016
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef PROGRESS_H_
#define PROGRESS_H_


#ifdef __cplusplus
#include <math.h>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iomanip>
#include <omp.h>
#include "log.h"
#endif

#ifdef MPIx
#include <mpi.h>
#endif

class Geometry;

class Progress {

private:

  std::string _name;
  int _counter;
  int _num_iterations;
  int _curr_interval;
  std::vector<int> _intervals;
  Geometry* _geometry;
  bool _mpi_comm;

public:
  Progress(int num_iterations, std::string name, double interval=0.1,
           Geometry* geometry=NULL, bool mpi_comm=false);
  virtual ~Progress();

  /* Worker functions */
  void incrementCounter();
  void reset();
};

#endif /* PROGRESS_H_ */
