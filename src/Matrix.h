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
#endif


class Matrix {

private:

  /** A list of lists representing the matrix */
  std::vector< std::map<int, FP_PRECISION> > _LIL;

  /** The CSR matrix variables */
  FP_PRECISION* _A;
  int* _IA;
  int* _JA;
  FP_PRECISION* _DIAG;
  
  bool _modified;
  int _num_x;
  int _num_y;
  int _num_z;
  int _num_groups;
  int _num_rows;
  
public:
  Matrix(int num_x=1, int num_y=1, int num_z=1, int num_groups=1);
  virtual ~Matrix();

  /* Worker functions */
  void incrementValue(int row, int col, FP_PRECISION val);
  void incrementValueByCoords(int x_from, int y_from, int z_from, int g_from,
                              int x_to, int y_to, int z_to, int g_to, FP_PRECISION val);
  void incrementValueByCell(int cell_from, int g_from, int cell_to, int g_to, FP_PRECISION val);
  void setValue(int row, int col, FP_PRECISION val);
  void setValueByCoords(int x_from, int y_from, int z_from, int g_from,
                        int x_to, int y_to, int z_to, int g_to, FP_PRECISION val);
  void setValueByCell(int cell_from, int g_from, int cell_to, int g_to, FP_PRECISION val);
  void clear();
  void convertToCSR();  
  void printString();
  void random();
  
  /* Getter functions */
  FP_PRECISION getValue(int row, int col);
  FP_PRECISION getValueByCoords(int x_from, int y_from, int z_from, int g_from,
                          int x_to, int y_to, int z_to, int g_to);
  FP_PRECISION getValueByCell(int cell_from, int g_from, int cell_to, int g_to);
  FP_PRECISION* getA();
  int* getIA();
  int* getJA();
  FP_PRECISION* getDIAG();
  int getNumX();
  int getNumY();
  int getNumZ();
  int getNumGroups();
  int getNumRows();
  int getNNZ();
  
};

#endif /* MATRIX_H_ */
