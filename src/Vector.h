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
  double* _array;
  int _num_rows;
  int _num_x;
  int _num_y;
  int _num_z;
  int _num_groups;
  
public:
  Vector(int num_x=1, int num_y=1, int num_z=1, int num_groups=1);
  virtual ~Vector();

  /* Worker functions */
  void incrementValue(int row, double val);
  void incrementValueByCell(int cell, int g, double val);
  void incrementValueByCoords(int x, int y, int z, int g, double val);
  void setValue(int row, double val);
  void setValueByCoords(int x, int y, int z, int g, double val);
  void setValueByCell(int cell, int g, double val);
  void setAll(double val);
  void clear();
  void scaleByValue(double val);  
  std::string toString();
  void printString();
  void copyTo(Vector* vector);
  void random();
  double sum();
  
  /* Getter functions */
  double getValue(int row);
  double getValueByCoords(int x, int y, int z, int g);
  double getValueByCell(int cell, int g=0);
  double* getArray();
  int getNumX();
  int getNumY();
  int getNumZ();
  int getNumGroups();
  int getNumRows();

};

#endif /* VECTOR_H_ */
