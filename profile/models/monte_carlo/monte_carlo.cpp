/*
 @file    monte_carlo.cpp
 @brief   creates geometry and materials to run a Monte Carlo simulation
 @author  Luke Eure
 @date    January 9 2016
*/

#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "../../../src/Tally.h"
#include "../../../src/Neutron.h"
#include "../../../src/MCSolver.h"
#include "../../../src/Universe.h"
#include "../../../src/Material.h"
#include "../../../src/Cell.h"
#include "../../../src/Universe.h"
#include "../../../src/Geometry.h"

int main() {

/*  // number of energy groups
  int num_groups = 7;

  // create fuel
  Material* fuel = new Material(1, "UO2");
  fuel->setNumEnergyGroups(num_groups);
  
  fuel->setSigmaTByGroup(1.779490E-01, 1);
  fuel->setSigmaTByGroup(3.298050E-01, 2);
  fuel->setSigmaTByGroup(4.803880E-01, 3);
  fuel->setSigmaTByGroup(5.543670E-01, 4);
  fuel->setSigmaTByGroup(3.118010E-01, 5);
  fuel->setSigmaTByGroup(3.951680E-01, 6);
  fuel->setSigmaTByGroup(5.644060E-01, 7);

  fuel->setSigmaFByGroup(7.212060E-03, 1);
  fuel->setSigmaFByGroup(8.193010E-04, 2);
  fuel->setSigmaFByGroup(6.453200E-03, 3);
  fuel->setSigmaFByGroup(1.856480E-02, 4);
  fuel->setSigmaFByGroup(1.780840E-02, 5);
  fuel->setSigmaFByGroup(8.303480E-02, 6);
  fuel->setSigmaFByGroup(2.160040E-01, 7);

  fuel->setNuSigmaFByGroup(2.005998E-02, 1);
  fuel->setNuSigmaFByGroup(2.027303E-03, 2);
  fuel->setNuSigmaFByGroup(1.570599E-02, 3);
  fuel->setNuSigmaFByGroup(4.518301E-02, 4);
  fuel->setNuSigmaFByGroup(4.334208E-02, 5);
  fuel->setNuSigmaFByGroup(2.020901E-01, 6);
  fuel->setNuSigmaFByGroup(5.257105E-01, 7);

  fuel->setChiByGroup(5.87910E-01, 1);
  fuel->setChiByGroup(4.11760E-01, 2);
  fuel->setChiByGroup(3.39060E-04, 3);
  fuel->setChiByGroup(1.17610E-07, 4);
  fuel->setChiByGroup(0., 5);
  fuel->setChiByGroup(0., 6);
  fuel->setChiByGroup(0., 7);

  fuel->setSigmaSByGroup(1.275370E-01, 1, 1);
  fuel->setSigmaSByGroup(4.237800E-02, 1, 2);
  fuel->setSigmaSByGroup(9.437400E-06, 1, 3);
  fuel->setSigmaSByGroup(5.516300E-09, 1, 4);
  fuel->setSigmaSByGroup(0., 1, 5);
  fuel->setSigmaSByGroup(0., 1, 6);
  fuel->setSigmaSByGroup(0., 1, 7);
  fuel->setSigmaSByGroup(0., 2, 1);
  fuel->setSigmaSByGroup(3.244560E-01, 2, 2);
  fuel->setSigmaSByGroup(1.631400E-03, 2, 3);
  fuel->setSigmaSByGroup(3.142700E-09, 2, 4);
  fuel->setSigmaSByGroup(0., 2, 5);
  fuel->setSigmaSByGroup(0., 2, 6);
  fuel->setSigmaSByGroup(0., 2, 7);
  fuel->setSigmaSByGroup(0., 3, 1);
  fuel->setSigmaSByGroup(0., 3, 2);
  fuel->setSigmaSByGroup(4.509400E-01, 3, 3);
  fuel->setSigmaSByGroup(2.679200E-03, 3, 4);
  fuel->setSigmaSByGroup(0., 3, 5);
  fuel->setSigmaSByGroup(0., 3, 6);
  fuel->setSigmaSByGroup(0., 3, 7);
  fuel->setSigmaSByGroup(0., 4, 1);
  fuel->setSigmaSByGroup(0., 4, 2);
  fuel->setSigmaSByGroup(0., 4, 3);
  fuel->setSigmaSByGroup(4.525650E-01, 4, 4);
  fuel->setSigmaSByGroup(5.566400E-03, 4, 5);
  fuel->setSigmaSByGroup(0., 4, 6);
  fuel->setSigmaSByGroup(0., 4, 7);
  fuel->setSigmaSByGroup(0., 5, 1);
  fuel->setSigmaSByGroup(0., 5, 2);
  fuel->setSigmaSByGroup(0., 5, 3);
  fuel->setSigmaSByGroup(1.252500E-04, 5, 4);
  fuel->setSigmaSByGroup(2.714010E-01, 5, 5);
  fuel->setSigmaSByGroup(1.025500E-02, 5, 6);
  fuel->setSigmaSByGroup(1.002100E-08, 5, 7);
  fuel->setSigmaSByGroup(0., 6, 1);
  fuel->setSigmaSByGroup(0., 6, 2);
  fuel->setSigmaSByGroup(0., 6, 3);
  fuel->setSigmaSByGroup(0., 6, 4);
  fuel->setSigmaSByGroup(1.296800E-03, 6, 5);
  fuel->setSigmaSByGroup(2.658020E-01, 6, 6);
  fuel->setSigmaSByGroup(1.680900E-02, 6, 7);
  fuel->setSigmaSByGroup(0., 7, 1);
  fuel->setSigmaSByGroup(0., 7, 2);
  fuel->setSigmaSByGroup(0., 7, 3);
  fuel->setSigmaSByGroup(0., 7, 4);
  fuel->setSigmaSByGroup(0., 7, 5);
  fuel->setSigmaSByGroup(8.545800E-03, 7, 6);
  fuel->setSigmaSByGroup(2.730800E-01, 7, 7);

  // create moderator
  Material* moderator = new Material(2, "Water");
  moderator->setNumEnergyGroups(num_groups);

  moderator->setSigmaTByGroup(1.592060E-01, 1);
  moderator->setSigmaTByGroup(4.129700E-01, 2);
  moderator->setSigmaTByGroup(5.903100E-01, 3);
  moderator->setSigmaTByGroup(5.843500E-01, 4);
  moderator->setSigmaTByGroup(7.180000E-01, 5);
  moderator->setSigmaTByGroup(1.254450E+00, 6);
  moderator->setSigmaTByGroup(2.650380E+00, 7);
  
  moderator->setSigmaSByGroup(4.447770E-02, 1, 1);
  moderator->setSigmaSByGroup(1.134000E-01, 1, 2);
  moderator->setSigmaSByGroup(7.234700E-04, 1, 3);
  moderator->setSigmaSByGroup(3.749900E-06, 1, 4);
  moderator->setSigmaSByGroup(5.318400E-08, 1, 5);
  moderator->setSigmaSByGroup(0., 1, 6);
  moderator->setSigmaSByGroup(0., 1, 7);
  moderator->setSigmaSByGroup(0., 2, 1);
  moderator->setSigmaSByGroup(2.823340E-01, 2, 2);
  moderator->setSigmaSByGroup(1.299400E-01, 2, 3);
  moderator->setSigmaSByGroup(6.234000E-04, 2, 4);
  moderator->setSigmaSByGroup(4.800200E-05, 2, 5);
  moderator->setSigmaSByGroup(7.448600E-06, 2, 6);
  moderator->setSigmaSByGroup(1.045500E-06, 2, 7);
  moderator->setSigmaSByGroup(0., 3, 1);
  moderator->setSigmaSByGroup(0., 3, 2);
  moderator->setSigmaSByGroup(3.452560E-01, 3, 3);
  moderator->setSigmaSByGroup(2.245700E-01, 3, 4);
  moderator->setSigmaSByGroup(1.699900E-02, 3, 5);
  moderator->setSigmaSByGroup(2.644300E-03, 3, 6);
  moderator->setSigmaSByGroup(5.034400E-04, 3, 7);
  moderator->setSigmaSByGroup(0., 4, 1);
  moderator->setSigmaSByGroup(0., 4, 2);
  moderator->setSigmaSByGroup(0., 4, 3);
  moderator->setSigmaSByGroup(9.102840E-02, 4, 4);
  moderator->setSigmaSByGroup(4.155100E-01, 4, 5);
  moderator->setSigmaSByGroup(6.373200E-02, 4, 6);
  moderator->setSigmaSByGroup(1.213900E-02, 4, 7);
  moderator->setSigmaSByGroup(0., 5, 1);
  moderator->setSigmaSByGroup(0., 5, 2);
  moderator->setSigmaSByGroup(0., 5, 3);
  moderator->setSigmaSByGroup(7.143700E-05, 5, 4);
  moderator->setSigmaSByGroup(1.391380E-01, 5, 5);
  moderator->setSigmaSByGroup(5.118200E-01, 5, 6);
  moderator->setSigmaSByGroup(6.122900E-02, 5, 7);
  moderator->setSigmaSByGroup(0., 6, 1);
  moderator->setSigmaSByGroup(0., 6, 2);
  moderator->setSigmaSByGroup(0., 6, 3);
  moderator->setSigmaSByGroup(0., 6, 4);
  moderator->setSigmaSByGroup(2.215700E-03, 6, 5);
  moderator->setSigmaSByGroup(6.999130E-01, 6, 6);
  moderator->setSigmaSByGroup(5.373200E-01, 6, 7);
  moderator->setSigmaSByGroup(0., 7, 1);
  moderator->setSigmaSByGroup(0., 7, 2);
  moderator->setSigmaSByGroup(0., 7, 3);
  moderator->setSigmaSByGroup(0., 7, 4);
  moderator->setSigmaSByGroup(0., 7, 5);
  moderator->setSigmaSByGroup(1.324400E-01, 7, 6);
  moderator->setSigmaSByGroup(2.480700E+00, 7, 7);

  moderator->setNuSigmaFByGroup(0., 1);
  moderator->setNuSigmaFByGroup(0., 2);
  moderator->setNuSigmaFByGroup(0., 3);
  moderator->setNuSigmaFByGroup(0., 4);
  moderator->setNuSigmaFByGroup(0., 5);
  moderator->setNuSigmaFByGroup(0., 6);
  moderator->setNuSigmaFByGroup(0., 7);

  moderator->setSigmaFByGroup(0., 1);
  moderator->setSigmaFByGroup(0., 2);
  moderator->setSigmaFByGroup(0., 3);
  moderator->setSigmaFByGroup(0., 4);
  moderator->setSigmaFByGroup(0., 5);
  moderator->setSigmaFByGroup(0., 6);
  moderator->setSigmaFByGroup(0., 7);

  moderator->setChiByGroup(0., 1);
  moderator->setChiByGroup(0., 2);
  moderator->setChiByGroup(0., 3);
  moderator->setChiByGroup(0., 4);
  moderator->setChiByGroup(0., 5);
  moderator->setChiByGroup(0., 6);
  moderator->setChiByGroup(0., 7);

  moderator->setSigmaFByGroup(0.0, 1);
  moderator->setSigmaFByGroup(0.0, 2);
  moderator->setNuSigmaFByGroup(0.0, 1);
  moderator->setNuSigmaFByGroup(0.0, 2);
  moderator->setChiByGroup(0.0, 1);
  moderator->setChiByGroup(0.0, 2);
  */
  int num_groups = 1;

  Material* fuel = new Material(1, "UO2");
  fuel->setNumEnergyGroups(num_groups);
  fuel->setSigmaTByGroup(2, 1);
  fuel->setSigmaFByGroup(1.7, 1);
  fuel->setNuSigmaFByGroup(2, 1);
  fuel->setChiByGroup(1, 1);
  fuel->setSigmaSByGroup(0, 1, 1);

  // create openmoc surfaces and set their boundary types
//  ZCylinder* z_cylinder = new ZCylinder(0.0, 0.0, 1.0, 4, "pin");
  XPlane* x_min = new XPlane(-2.0, 0, "x_min");
  XPlane* x_max = new XPlane(2.0, 1, "x_max");
  YPlane* y_min = new YPlane(-2.0, 2, "y_min");
  YPlane* y_max = new YPlane(2.0, 3, "y_max");

  x_min->setBoundaryType(REFLECTIVE);
  x_max->setBoundaryType(REFLECTIVE);
  y_min->setBoundaryType(REFLECTIVE);
  y_max->setBoundaryType(REFLECTIVE);
  std::cout << "num groups " << num_groups << std::endl;

  /*
  // create cells
  Cell* fuel_cell = new Cell(2, "fuel");
  fuel_cell->setFill(fuel);
  fuel_cell->addSurface(-1, z_cylinder);

  Cell* moderator_cell = new Cell(1, "moderator");
  moderator_cell->setFill(moderator);
  moderator_cell->addSurface(1, x_min);
  moderator_cell->addSurface(-1, x_max);
  moderator_cell->addSurface(1, y_min);
  moderator_cell->addSurface(-1, y_max);
  moderator_cell->addSurface(1, z_cylinder);
*/
  Cell* fuel_cell = new Cell(2, "fuel");
  fuel_cell->setFill(fuel);
  fuel_cell->addSurface(1, x_min);
  fuel_cell->addSurface(-1, x_max);
  fuel_cell->addSurface(1, y_min);
  fuel_cell->addSurface(-1, y_max);

  // create universes
  Universe* root_universe = new Universe(0, "root universe");
//  root_universe->addCell(moderator_cell);
  root_universe->addCell(fuel_cell);

  // create geometry
  Geometry* geometry = new Geometry();
  geometry->setRootUniverse(root_universe);

  // simulate neutron histories
  int num_neutrons = 1000000;
  int num_batches = 100;

  std::cout << "makes solver\n";
  // create solver and run simulation
  MCSolver solver;
  solver.setGeometry(geometry);
  solver.initialize();

  std::cout << "computes eigenvalue\n";
  solver.computeEigenvalue(num_neutrons, num_batches, num_groups);

  std::cout << std::endl;

  delete geometry;
  delete x_min;
  delete x_max;
  delete y_min;
  delete y_max;

  return 0;
}
