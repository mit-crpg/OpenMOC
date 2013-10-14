/**
 * @file Cmfd.h
 * @brief The Cmfd class.
 * @date October 14, 2013
 * @author Sam Shaner, MIT, Course 22 (shaner@mit.edu)
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

  /* physical mesh size */
  double _length_x;
  double _length_y;

  /* number of cells in x and y directions */
  int _cells_x;
  int _cells_y;

  /* number of groups */
  int _num_groups;

  /* number of surface current values */
  int _num_currents;

  /* number of fsrs */
  int _num_fsrs;

  /* number of azim angles */
  int _num_azim;

  /* array of boundary enums */
  boundaryType* _boundaries;

  /* array of mesh cell volumes */
  double* _volumes;

  /* array of mesh surface currents */
  FP_PRECISION* _currents;

  /* vector of vectors of fsrs in each mesh cell */
  std::vector< std::vector<int> > _cell_fsrs;

  /* acceleration flag */
  bool _acceleration;

  /* map of fluxes mapped onto mesh */
  std::map<std::string, FP_PRECISION*> _fluxes;

  /* materials array */
  Material** _materials;

  /* array of fsr bounds */
  int* _fsr_indices;

  /* array of lenghts of each mesh cell in x and y directions */
  double* _lengths_x;
  double* _lengths_y;

  /* array of cell bounds in x and y direction */
  double* _bounds_x;
  double* _bounds_y;

public:
  Mesh(bool acceleration=false);
  virtual ~Mesh();
  void initialize();
  void setFSRBounds();
  void setCellBounds();

  /* get mesh parameters */
  double getLengthX();
  double getLengthY();
  int getCellsX();
  int getCellsY();
  int getNumCells();
  boundaryType getBoundary(int side);
  int getNumCurrents();
  FP_PRECISION getFlux(std::string flux_name, int cell_id, int group);
  std::vector<std::vector<int> >* getCellFSRs();
  Material** getMaterials();
  double* getVolumes();
  FP_PRECISION* getFluxes(std::string flux_name);
  double* getLengthsX();
  double* getLengthsY();
  FP_PRECISION* getCurrents();

  /* set mesh parameters */
  void setLengthX(double length_x);
  void setLengthY(double length_y);
  void setCellLengthX(int cell_num, double length_x);
  void setCellLengthY(int cell_num, double length_y);
  void setCellsX(int cells_x);
  void setCellsY(int cells_y);
  void setSurfaceCurrents(FP_PRECISION* surface_currents);
  void setVolume(double volume, int cell_num);

  /* set general problem specs */
  void setNumGroups(int num_groups);
  void setNumAzim(int num_azim);
  void setNumFSRs(int num_fsrs);
  bool getAcceleration();
  
  /* get generation problem specs */
  int getNumGroups();
  int getNumFSRs();

  /* worker functions */
  int findMeshCell(double x, double y);
  int findMeshSurface(int fsr_id, LocalCoords* coord, int angle);
  void printCurrents();
  void splitCorners();
  void setBoundary(int side, boundaryType boundary);
  int getCellNext(int cell_num, int surface_id);
  int findCellId(LocalCoords* coord);
  void initializeMaterials(std::map<int, Material*>* materials, int* fsrs_to_mats);
};

#endif /* MESH_H_ */
