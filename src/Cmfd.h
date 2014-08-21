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
#include <utility>
#include <math.h>
#include <limits.h>
#include <string>
#include <sstream>
#include <queue>
#include <iostream>
#include <fstream>
#include "Quadrature.h"
#include "log.h"
#include "Material.h"
#include "Surface.h"
#include "Timer.h"
#include "Universe.h"
#endif


/**
 * @class Cmfd Cmfd.h "src/Cmfd.h"
 * @brief A class for Coarse Mesh Finite Difference (CMFD) acceleration.
 */
class Cmfd {

private:

  /** Pointer to polar Quadrature object */
  Quadrature* _quad;

  /** The keff eigenvalue */
  double _k_eff;

  /** The A (destruction) matrix */
  double** _A;

  /** The M (production) matrix */
  double** _M;

  /** The old source vector */
  double* _old_source;

  /** The new source vector */
  double* _new_source;

  /** flux vectors */
  double* _new_flux;
  double* _old_flux;
  double* _flux_temp;

  /** Gauss-Seidel SOR relaxation factor */
  double _SOR_factor;

  /** L2 norm */
  double _l2_norm;

  /** Keff convergence criteria */
  double _conv_criteria;

  /** Number of cells in x direction */
  int _cx;

  /** Number of cells in y direction */
  int _cy;

  /** Number of energy groups */
  int _num_groups;

  /** Number of coarse energy groups */
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
  double* _volumes;

  /** Array of material pointers for CMFD cell materials */
  Material** _materials;

  /** Physical dimensions of the geometry and each CMFD cell */
  double _width;
  double _height;
  double _cell_width;
  double _cell_height;

  /** Array of geometry boundaries */
  boundaryType* _boundaries;

  /** Array of surface currents for each CMFD cell */
  double* _currents;

  /** Vector of vectors of FSRs containing in each cell */
  std::vector< std::vector<int> > _cell_fsrs;

  /** Flag indicating whether CMFD mesh should be overlaid on geometry */
  bool _overlay_mesh;

  /** MOC flux update relaxation factor */
  double _relax_factor;

  /** Flag indicating whether to use optically thick correction factor */
  bool _optically_thick;

  /** Pointer to Lattice object representing the CMFD mesh */
  Lattice* _lattice;

  /** Flag indicating whether to update the MOC flux */
  bool _update_flux;

public:

  Cmfd(double criteria=1e-8);
  virtual ~Cmfd();

  /* Worker functions */
  void constructMatrices();
  void computeDs();
  void computeXS();
  void updateMOCFlux();
  double computeDiffCorrect(double d, double h);
  double computeKeff();
  void initializeCellMap();
  void initializeGroupMap();
  void initializeFlux();
  void initializeMaterials();
  void rescaleFlux();
  void linearSolve(double** mat, double* vec_x, double* vec_b,
                   double conv, int max_iter=10000);
  void splitCorners();
  int getCellNext(int cell_num, int surface_id);
  int findCmfdCell(LocalCoords* coords);
  int findCmfdSurface(int cell, LocalCoords* coords);
  void addFSRToCell(int cmfd_cell, int fsr_id);

  /* Matrix and Vector functions */
  void matZero(double** mat, int width);
  void vecCopy(double* vec_from, double* vec_to);
  double vecSum(double* vec);
  void matMultM(double** mat, double* vec_x, double* vec_y);
  void vecSet(double* vec, double val);
  void vecScale(double* vec, double scale_val);

  /* Get parameters */
  double getKeff();
  int getNumCmfdGroups();
  int getNumGroups();
  int getNumCells();
  int getCmfdGroup(int group);
  bool getOverlayMesh();
  bool getOpticallyThick();
  double getMOCRelaxationFactor();
  boundaryType getBoundary(int side);
  Lattice* getLattice();
  int getNumX();
  int getNumY();
  int convertFSRIdToCmfdCell(int fsr_id);
  std::vector< std::vector<int> > getCellFSRs();
  bool getUpdateFlux();

  /* Set parameters */
  void setSORRelaxationFactor(double SOR_factor);
  void setWidth(double width);
  void setHeight(double height);
  void setNumX(int cx);
  void setNumY(int cy);
  void setSurfaceCurrents(double* surface_currents);
  void setNumFSRs(int num_fsrs);
  void setNumGroups(int num_groups);
  void setOverlayMesh(bool overlay_mesh);
  void setOpticallyThick(bool thick);
  void setMOCRelaxationFactor(double relax_factor);
  void setBoundary(int side, boundaryType boundary);
  void setLattice(Lattice* lattice);
  void setLatticeStructure(int num_x, int num_y);
  void setUpdateFlux(bool update_flux);
  void setGroupStructure(int* group_indices, int ncg);

  /* Set FSR parameters */
  void setFSRMaterials(Material** FSR_materials);
  void setFSRVolumes(FP_PRECISION* FSR_volumes);
  void setFSRFluxes(FP_PRECISION* scalar_flux);
  void setCellFSRs(std::vector< std::vector<int> > cell_fsrs);
};

#endif /* CMFD_H_ */
