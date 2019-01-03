/**
 * @file Matrix.h
 * @brief A matrix object
 * @date May 5, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef MATRIX_H_
#define MATRIX_H_


#ifdef __cplusplus
#include <math.h>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include "log.h"
#include "constants.h"
#endif


class Matrix {

private:

  /** A vector of maps representing the matrix */
  std::vector< std::map<int, CMFD_PRECISION> > _LIL;

  /** The CSR matrix variables */
  CMFD_PRECISION* _A;
  int* _IA;
  int* _JA;
  CMFD_PRECISION* _DIAG;

  bool _modified;
  int _num_x;
  int _num_y;
  int _num_z;
  int _num_groups;
  int _num_rows;

  /** OpenMP mutual exclusion locks for atomic cell updates */
  omp_lock_t* _cell_locks;

  void convertToCSR();
  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumZ(int num_z);
  void setNumGroups(int num_groups);

public:
  Matrix(omp_lock_t* cell_locks, int num_x=1, int num_y=1, int num_z=1,
         int num_groups=1);
  virtual ~Matrix();

  /* Worker functions */
  void incrementValue(int cell_from, int group_from, int cell_to, int group_to,
                      CMFD_PRECISION val);
  void clear();
  void printString();
  void transpose();

  /* Getter functions */
  CMFD_PRECISION getValue(int cell_from, int group_from, int cell_to,
                        int group_to);
  CMFD_PRECISION* getA();
  int* getIA();
  int* getJA();
  CMFD_PRECISION* getDiag();
  int getNumX();
  int getNumY();
  int getNumZ();
  int getNumGroups();
  int getNumRows();
  int getNNZ();
  omp_lock_t* getCellLocks();

  /* Setter functions */
  void setValue(int cell_from, int group_from, int cell_to, int group_to,
                CMFD_PRECISION val);
};

#endif /* MATRIX_H_ */
