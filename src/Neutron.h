/* 
 @file    Neutron.h
 @brief   contains Neutron class
 @author  Luke Eure
 @date    January 8 2016
*/

#ifndef NEUTRON_H
#define NEUTRON_H

#include <vector>
#include <cmath>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "Point.h"

class Neutron {
public:
  Neutron(int neutron_num);
  virtual ~Neutron();

  void kill();
  void move(double distance);
  void reflect(int axis);
  void setGroup(int new_group);
  void setPosition(int axis, double value);
  void sampleDirection();
  void getPositionVector(Point* &position);
  void setDirection(int axis, double magnitude);
  double arand();
  double getDirection(int axis);
  double getDistance(Point *coord);
  double getPosition(int axis);
  bool alive();
  int getGroup();
  int rand();
  int sampleEnergyGroup(std::vector <double> chi);
  int sampleScatteredGroup(std::vector <double> &scattering_matrix);

private:
  
  /** tells if the neutron is alive */
  bool _neutron_alive;

  /** energy group of the neutron */
  int _neutron_group;

  /** position of the neutron */
  Point _xyz;

  /** direction of travel of the neutron */
  std::vector <double> _neutron_direction;

  /** identification number */
  int _id;

  /** to be used in calling rand_r() */
  unsigned int _seed;
};

#endif
