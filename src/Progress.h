/**
 * @file Progress.h
 * @brief A progress object
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


class Progress {

private:

  /** A list of lists representing the vector */
  std::string _name;
  int _counter;
  int _num_iterations;
  double _interval;
  double _last_interval;

public:
  Progress(int num_iterations, std::string name, double interval=0.1);
  virtual ~Progress();

  /* Worker functions */
  void incrementCounter();
  void reset();
};

#endif /* PROGRESS_H_ */
