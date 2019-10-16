/**
 * @file Vector.h
 * @brief A vector object
 * @date May 5, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef SRC_VECTOR_H_
#define SRC_VECTOR_H_


#ifdef __cplusplus
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "log.h"
#include "pairwise_sum.h"
#endif


class Vector {
 private:
  /** A list of lists representing the vector */
  CMFD_PRECISION* _array;
  int _num_rows;
  int _num_x;
  int _num_y;
  int _num_z;
  int _num_groups;

  /** OpenMP mutual exclusion locks for atomic cell updates */
  omp_lock_t* _cell_locks;

  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumZ(int num_z);
  void setNumGroups(int num_groups);

 public:
  Vector(omp_lock_t* cell_locks, int num_x = 1, int num_y = 1, int num_z = 1,
         int num_groups = 1);
  virtual ~Vector();

  /* Worker functions */
  void incrementValue(int cell, int group, CMFD_PRECISION val);
  void incrementValues(int cell, int group_start, int group_end,
                       CMFD_PRECISION* vals);
  void clear();
  void scaleByValue(CMFD_PRECISION val);
  void printString();
  void copyTo(Vector* vector);

  /* Getter functions */
  CMFD_PRECISION* getArray();
  int getNumX();
  int getNumY();
  int getNumZ();
  int getNumGroups();
  int getNumRows();
  double getSum();
  long getNumNegativeValues();
  omp_lock_t* getCellLocks();

  /* Setter functions */
  void setValue(int cell, int group, CMFD_PRECISION val);
  void setValues(int cell, int group_start, int group_end,
                 CMFD_PRECISION* vals);
  void setAll(CMFD_PRECISION val);

/**
 * @brief Get a value at location described by a given cell and group index.
 * @param cell The cell location index.
 * @param group The group location index.
 */
  inline CMFD_PRECISION getValue(int cell, int group) {
    return _array[cell*_num_groups + group];
  }

};

#endif  // SRC_VECTOR_H_
