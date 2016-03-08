/**
 * @file Vector.h
 * @brief A vector object
 * @date May 5, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef VECTOR_H_
#define VECTOR_H_


#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include <math.h>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include "log.h"
#include "pairwise_sum.h"
#endif


class Vector {

private:

  /** An array representing the vector */
  FP_PRECISION* _array;
  int _num_rows;
  int _num_x;
  int _num_y;
  int _num_groups;

  /** OpenMP mutual exclusion locks for atomic cell updates */
  omp_lock_t* _cell_locks;

  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumGroups(int num_groups);

public:
  Vector(omp_lock_t* cell_locks, int num_x=1, int num_y=1, int num_groups=1);
  virtual ~Vector();

  /* Worker functions */
  void incrementValue(int cell, int group, FP_PRECISION val);
  void incrementValues(int cell, int group_start, int group_end,
      FP_PRECISION* vals);
  void clear();
  void scaleByValue(FP_PRECISION val);
  void printString();
  void copyTo(Vector* vector);

  /* Getter functions */
  FP_PRECISION getValue(int cell, int group);
  FP_PRECISION* getArray();
  int getNumX();
  int getNumY();
  int getNumGroups();
  int getNumRows();
  FP_PRECISION getSum();
  omp_lock_t* getCellLocks();

  /* Setter functions */
  void setValue(int cell, int group, FP_PRECISION val);
  void setValues(int cell, int group_start, int group_end, FP_PRECISION* vals);
  void setAll(FP_PRECISION val);
};

#endif /* VECTOR_H_ */
