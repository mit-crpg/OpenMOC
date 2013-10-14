/*
 * Mesh.h
 *
 *  Created on: May 6, 2012
 *      Author: samuelshaner
 */

#ifndef MESH_H_
#define MESH_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <utility>
#include <sstream>
#include "log.h"
#include "LocalCoords.h"
#include "Surface.h"
#include "Material.h"

class Mesh {
private:
  double _length_x;
  double _length_y;
  int _cells_x;
  int _cells_y;
  int _num_groups;
  int _num_fsrs;
  boundaryType* _boundaries;
  double* _volumes;
  FP_PRECISION* _currents;
  std::vector< std::vector<int> > _cell_fsrs;
  bool _acceleration;
  std::map<std::string, FP_PRECISION*> _fluxes;
  Material** _materials;
  int* _fsr_indices;
  double* _lengths_x;
  double* _lengths_y;
  double* _bounds_x;
  double* _bounds_y;


public:
  Mesh(bool acceleration=false);
  virtual ~Mesh();
  void initialize();
  double getLengthX();
  double getLengthY();
  void setLengthX(double length_x);
  void setLengthY(double length_y);
  void setCellLengthX(int cell_num, double length_x);
  void setCellLengthY(int cell_num, double length_y);
  int getCellsX();
  int getCellsY();
  void setCellsX(int cells_x);
  void setCellsY(int cells_y);
  void setFSRBounds();
  void setCellBounds();
  int findMeshCell(double x, double y);
  int findMeshSurface(int fsr_id, LocalCoords* coord);
  void printCurrents();
  void splitCorners();
  void setBoundary(int side, boundaryType boundary);
  int getNumCells();
  void setNumGroups(int num_groups);
  int getNumGroups();
  void setNumFSRs(int num_fsrs);
  int getNumFSRs();
  int getCellNext(int cell_num, int surface_id);
  boundaryType getBoundary(int side);
  FP_PRECISION getFlux(std::string flux_name, int cell_id, int group);
  int findCellId(LocalCoords* coord);
  void setSurfaceCurrents(FP_PRECISION* surface_currents);
  void setVolume(double volume, int cell_num);
  bool getAcceleration();
  std::vector<std::vector<int> >* getCellFSRs();
  Material** getMaterials();
  double* getVolumes();
  FP_PRECISION* getFluxes(std::string flux_name);
  double* getLengthsX();
  double* getLengthsY();
  FP_PRECISION* getCurrents();

  void initializeMaterials(std::map<int, Material*>* materials, int* fsrs_to_mats);
};

#endif /* MESH_H_ */
