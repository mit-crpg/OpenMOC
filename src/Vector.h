/**
 * @file Vector.h
 * @brief A vector object
 * @date May 5, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef VECTOR_H_
#define VECTOR_H_


#ifdef __cplusplus
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

  /** A list of lists representing the vector */
  FP_PRECISION* _array;
  int _num_rows;
  int _num_x;
  int _num_y;
  int _num_z;
  int _num_groups;
  
public:
  Vector(int num_x=1, int num_y=1, int num_z=1, int num_groups=1);
  virtual ~Vector();

  /* Worker functions */
  void incrementValue(int row, FP_PRECISION val);
  void incrementValueByCell(int cell, int g, FP_PRECISION val);
  void incrementValueByCoords(int x, int y, int z, int g, FP_PRECISION val);
  void setValue(int row, FP_PRECISION val);
  void setValueByCoords(int x, int y, int z, int g, FP_PRECISION val);
  void setValueByCell(int cell, int g, FP_PRECISION val);
  void setAll(FP_PRECISION val);
  void clear();
  void scaleByValue(FP_PRECISION val);  
  std::string toString();
  void printString();
  void copyTo(Vector* vector);
  void random();
  FP_PRECISION sum();
  
  /* Getter functions */
  FP_PRECISION getValue(int row);
  FP_PRECISION getValueByCoords(int x, int y, int z, int g);
  FP_PRECISION getValueByCell(int cell, int g=0);
  FP_PRECISION* getArray();
  int getNumX();
  int getNumY();
  int getNumZ();
  int getNumGroups();
  int getNumRows();

};

#endif /* VECTOR_H_ */
