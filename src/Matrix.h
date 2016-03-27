/**
 * @file Matrix.h
 * @brief A matrix object
 * @date May 5, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef MATRIX_H_
#define MATRIX_H_


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
#include <iomanip>
#include "log.h"
#endif


class Matrix {

private:

  /** A list of lists representing the matrix */
  std::vector< std::map<int, FP_PRECISION> > _LIL;

  /** The CSR matrix variables */
  FP_PRECISION* _A;
  FP_PRECISION* _LU;
  int* _IA;
  int* _JA;
  int* _ILU;
  int* _JLU;
  FP_PRECISION* _DIAG;

  bool _modified;
  int _num_x;
  int _num_y;
  int _num_groups;
  int _num_rows;
  int _NNZ;
  int _NNZLU;

  /** OpenMP mutual exclusion locks for atomic cell updates */
  omp_lock_t* _cell_locks;

  void convertToCSR();
  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumGroups(int num_groups);

public:
  Matrix(omp_lock_t* cell_locks, int num_x=1, int num_y=1, int num_groups=1);
  virtual ~Matrix();

  /* Worker functions */
  void incrementValue(int cell_from, int group_from, int cell_to, int group_to,
                      FP_PRECISION val);
  void clear();
  void printString();
  void transpose();

  /* Getter functions */
  FP_PRECISION getValue(int cell_from, int group_from, int cell_to,
                        int group_to);
  FP_PRECISION* getA();
  FP_PRECISION* getLU();
  int* getIA();
  int* getILU();
  int* getJA();
  int* getJLU();
  FP_PRECISION* getDiag();
  int getNumX();
  int getNumY();
  int getNumGroups();
  int getNumRows();
  int getNNZ();
  int getNNZLU();
  omp_lock_t* getCellLocks();

  /* Setter functions */
  void setValue(int cell_from, int group_from, int cell_to, int group_to,
                FP_PRECISION val);
};

#endif /* MATRIX_H_ */
