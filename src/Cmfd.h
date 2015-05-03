/**
 * @file Cmfd.h
 * @brief The Cmfd class.
 * @date October 14, 2013
 * @author Sam Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef CMFD_H_
#define CMFD_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include "Python.h"
#include "log.h"
#include "Timer.h"
#include "Universe.h"
#include "Track.h"
#include "PolarQuad.h"
#include "linalg.h"
#include "pairwise_sum.h"
#include <utility>
#include <math.h>
#include <limits.h>
#include <string>
#include <sstream>
#include <queue>
#include <iostream>
#include <fstream>
#endif


/**
 * @class Cmfd Cmfd.h "src/Cmfd.h"
 * @brief A class for Coarse Mesh Finite Difference (CMFD) acceleration.
 */
class Cmfd {

private:

  /** Pointer to polar quadrature object */
  PolarQuad* _polar_quad;

  /** The keff eigenvalue */
  FP_PRECISION _k_eff;

  /** The A (destruction) matrix */
  FP_PRECISION** _A;

  /** The M (production) matrix */
  FP_PRECISION** _M;

  /** The old source vector */
  FP_PRECISION* _old_source;

  /** The new source vector */
  FP_PRECISION* _new_source;

  /** Vector representing the flux for each cmfd cell and cmfd enegy group at 
   * the end of a CMFD solve */
  FP_PRECISION* _new_flux;

  /** Vector representing the flux for each cmfd cell and cmfd enegy group at 
   * the beginning of a CMFD solve */
  FP_PRECISION* _old_flux;

  /** Vector representing the flux during the previous iteration of a 
   * cmfd solve */
  FP_PRECISION* _flux_temp;

  /** Gauss-Seidel SOR relaxation factor */
  FP_PRECISION _SOR_factor;

  /** cmfd source convergence threshold */
  FP_PRECISION _source_convergence_threshold;

  /** Number of cells in x direction */
  int _num_x;

  /** Number of cells in y direction */
  int _num_y;

  /** Number of energy groups */
  int _num_moc_groups;

  /** Number of polar angles */
  int _num_polar;

  /** Number of energy groups used in cmfd solver. Note that cmfd supports
   * energy condensation from the MOC */
  int _num_cmfd_groups;

  /** Coarse energy indices for fine energy groups */
  int* _group_indices;

  /** Map of MOC groups to CMFD groups */
  int* _group_indices_map;

  /** Number of FSRs */
  int _num_FSRs;

  /** The volumes (areas) for each FSR */
  FP_PRECISION* _FSR_volumes;

  /** Pointers to Materials for each FSR */
  Material** _FSR_materials;

  /** The FSR scalar flux in each energy group */
  FP_PRECISION* _FSR_fluxes;

  /** Array of CMFD cell volumes */
  FP_PRECISION* _volumes;

  /** Array of material pointers for CMFD cell materials */
  Material** _materials;

  /** Physical dimensions of the geometry and each CMFD cell */
  double _width;
  double _height;
  double _cell_width;
  double _cell_height;

  /** Array of geometry boundaries */
	//  int* _boundaries;
  boundaryType* _boundaries;

  /** Array of surface currents for each CMFD cell */
  FP_PRECISION* _surface_currents;

  /** Vector of vectors of FSRs containing in each cell */
  std::vector< std::vector<int> > _cell_fsrs;

  /** MOC flux update relaxation factor */
  FP_PRECISION _relax_factor;

  /** Flag indicating whether to use optically thick correction factor */
  bool _optically_thick;

  /** Pointer to Lattice object representing the CMFD mesh */
  Lattice* _lattice;

  /** Flag indicating whether to update the MOC flux */
  bool _flux_update_on;

public:

  Cmfd();
  virtual ~Cmfd();

  /* Worker functions */
  void constructMatrices();
  void computeDs(int moc_iteration);
  void computeXS();
  void updateMOCFlux();
  FP_PRECISION computeDiffCorrect(FP_PRECISION d, FP_PRECISION h);
  FP_PRECISION computeKeff(int moc_iteration);
  void initializeCellMap();
  void initializeGroupMap();
  void initializeFlux();
  void initializeMaterials();
  void rescaleFlux();
  void linearSolve(FP_PRECISION** mat, FP_PRECISION* vec_x, FP_PRECISION* vec_b,
                   FP_PRECISION conv, int max_iter=10000);
  void splitCorners();
  int getCellNext(int cell_num, int surface_id);
  int findCmfdCell(LocalCoords* coords);
  int findCmfdSurface(int cell, LocalCoords* coords);
  void addFSRToCell(int cmfd_cell, int fsr_id);
  void updateBoundaryFlux(Track** tracks, FP_PRECISION* boundary_flux, 
                          int num_tracks);

  /* Get parameters */
  int getNumCmfdGroups();
  int getNumMOCGroups();
  int getNumCells();
  int getCmfdGroup(int group);
  bool isOpticallyThick();
  FP_PRECISION getMOCRelaxationFactor();
  int getBoundary(int side);
  Lattice* getLattice();
  int getNumX();
  int getNumY();
  int convertFSRIdToCmfdCell(int fsr_id);
  std::vector< std::vector<int> > getCellFSRs();
  bool isFluxUpdateOn();
  FP_PRECISION getFluxRatio(int cmfd_cell, int moc_group);

  /* Set parameters */
  void setSORRelaxationFactor(FP_PRECISION SOR_factor);
  void setWidth(double width);
  void setHeight(double height);
  void setNumX(int num_x);
  void setNumY(int num_y);
  void setSurfaceCurrents(FP_PRECISION* surface_currents);
  void setNumFSRs(int num_fsrs);
  void setNumMOCGroups(int num_moc_groups);
  void setOpticallyThick(bool thick);
  void setMOCRelaxationFactor(FP_PRECISION relax_factor);
  void setBoundary(int side, boundaryType boundary);
  void setLattice(Lattice* lattice);
  void setLatticeStructure(int num_x, int num_y);
  void setFluxUpdateOn(bool flux_update_on);
  void setGroupStructure(int* group_indices, int length_group_indices);
  void setSourceConvergenceThreshold(FP_PRECISION source_thresh);
  void setPolarQuadrature(PolarQuad* polar_quad);
  
  /* Set FSR parameters */
  void setFSRMaterials(Material** FSR_materials);
  void setFSRVolumes(FP_PRECISION* FSR_volumes);
  void setFSRFluxes(FP_PRECISION* scalar_flux);
  void setCellFSRs(std::vector< std::vector<int> > cell_fsrs);
};

#endif /* CMFD_H_ */
