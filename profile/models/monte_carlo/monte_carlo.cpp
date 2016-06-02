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

  // number of energy groups
  int num_groups = 2;

  // create fuel
  Material* fuel = new Material(1, "fuel");
  fuel->setNumEnergyGroups(num_groups);
  fuel->setSigmaTByGroup(2.0/9.0, 1);
  fuel->setSigmaTByGroup(5.0/6.0, 2);
  fuel->setSigmaFByGroup(1.0/480.0, 1);
  fuel->setSigmaFByGroup(1.0/16.0, 2);
  fuel->setNuSigmaFByGroup(2.4/480.0, 1);
  fuel->setNuSigmaFByGroup(2.4/16.0, 2);
  fuel->setSigmaSByGroup(71.0/360.0, 1, 1);
  fuel->setSigmaSByGroup(.02, 1, 2);
  fuel->setSigmaSByGroup(0.0, 2, 1);
  fuel->setSigmaSByGroup(11.0/15.0, 2, 2);
  fuel->setChiByGroup(1.0, 1);
  fuel->setChiByGroup(0.0, 2);

  // create moderator
  Material* moderator = new Material(2, "moderator");
  moderator->setNumEnergyGroups(num_groups);
  moderator->setSigmaTByGroup(2.0/9.0, 1);
  moderator->setSigmaTByGroup(5.0/3.0, 2);
  moderator->setSigmaFByGroup(0.0, 1);
  moderator->setSigmaFByGroup(0.0, 2);
  moderator->setNuSigmaFByGroup(0.0, 1);
  moderator->setNuSigmaFByGroup(0.0, 2);
  moderator->setSigmaSByGroup(71.0/360.0, 1, 1);
  moderator->setSigmaSByGroup(.025, 1, 2);
  moderator->setSigmaSByGroup(0.0, 2, 1);
  moderator->setSigmaSByGroup(47.0/30.0, 2, 2);
  moderator->setChiByGroup(0.0, 1);
  moderator->setChiByGroup(0.0, 2);

  // create openmoc surfaces and set their boundary types
  XPlane* x_min = new XPlane(-2.0, 0, "x_min");
  XPlane* x_max = new XPlane(2.0, 1, "x_max");
  YPlane* y_min = new YPlane(-2.0, 2, "y_min");
  YPlane* y_max = new YPlane(2.0, 3, "y_max");

  x_min->setBoundaryType(VACUUM);
  x_max->setBoundaryType(VACUUM);
  y_min->setBoundaryType(VACUUM);
  y_max->setBoundaryType(VACUUM);

  // create cells
  Cell* root_cell = new Cell(0, "root");
  root_cell->addSurface(1, x_min);
  root_cell->addSurface(-1, x_max);
  root_cell->addSurface(1, y_min);
  root_cell->addSurface(-1, y_max);

  Cell* moderator_cell = new Cell(1, "moderator");
  moderator_cell->setFill(moderator);

  Cell* fuel_cell = new Cell(2, "fuel");
  fuel_cell->setFill(fuel);

  // create universes
  Universe* root_universe = new Universe(0, "root universe");
  root_universe->addCell(root_cell);

  Universe* moderator_universe = new Universe(1, "moderator universe");
  moderator_universe->addCell(moderator_cell);

  Universe* fuel_universe = new Universe(2, "fuel universe");
  fuel_universe->addCell(fuel_cell);

  // create lattice
  const int numXLat = 9;
  const int numYLat = 9;
  const int numZLat = 1;
  Lattice* lattice = new Lattice();
  lattice->setWidth(4.0/9.0, 4.0/9.0);

  // fill latice with universes
  Universe* matrix[numXLat*numYLat];

  int mold[numXLat*numYLat] = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 1, 2, 2, 2, 1, 1, 1,
                                1, 1, 1, 2, 2, 2, 1, 1, 1,
                                1, 1, 1, 2, 2, 2, 1, 1, 1,
                                1, 1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, 1, 1, 1, 1, 1};

  std::map<int, Universe*> names = {{1, moderator_universe},
      {2, fuel_universe}};
  for (int n=0; n<numXLat*numYLat; n++)
      matrix[n] = names[mold[n]];

  lattice->setUniverses(numZLat, numYLat, numXLat, matrix);

  // fill root cell with lattice
  root_cell->setFill(lattice);

  // create geometry
  Geometry* geometry = new Geometry();
  geometry->setRootUniverse(root_universe);

  // simulate neutron histories
  int num_neutrons = 100000;
  int num_batches = 3;

  // create solver
  MCSolver solver;
  solver.setGeometry(geometry);
  solver.setRootCell(root_cell);

  // initialize solver function
  solver.initialize();

  // simulate neutrons
  solver.computeEigenValue(num_neutrons, num_batches, num_groups);

  std::cout << std::endl;

  delete geometry;
  delete x_min;
  delete x_max;
  delete y_min;
  delete y_max;

  return 0;
}
