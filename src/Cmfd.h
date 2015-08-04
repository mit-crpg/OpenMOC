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
#ifdef SWIG
#include "Python.h"
#endif
#include "log.h"
#include "Timer.h"
#include "Universe.h"
#include "Track.h"
#include "PolarQuad.h"
#include "linalg.h"
#include "Geometry.h"
#include <utility>
#include <math.h>
#include <limits.h>
#include <string>
#include <sstream>
#include <queue>
#include <iostream>
#include <fstream>
#endif

/** Forward declaration of Geometry class */
class Geometry;

/** Comparitor for sorting k-nearest stencil std::pair objects */
inline bool stencilCompare(const std::pair<int, FP_PRECISION>& firstElem,
                           const std::pair<int, FP_PRECISION>& secondElem) {
  return firstElem.second < secondElem.second;
}

#undef track_flux

/** Indexing macro for the angular fluxes for each polar angle and energy
 *  group for either the forward or reverse direction for a given Track */ 
#define track_flux(p,e) (track_flux[(p)*_num_moc_groups + (e)])

/**
 * @class Cmfd Cmfd.h "src/Cmfd.h"
 * @brief A class for Coarse Mesh Finite Difference (CMFD) acceleration.
 */
class Cmfd {

private:

  /** Pointer to polar quadrature object */
  PolarQuad* _polar_quad;

  /** Pointer to geometry object */
  Geometry* _geometry;

  /** The keff eigenvalue */
  FP_PRECISION _k_eff;

  /** The A (destruction) matrix */
  Matrix* _A;

  /** The M (production) matrix */
  Matrix* _M;

  /** The old source vector */
  Vector* _old_source;

  /** The new source vector */
  Vector* _new_source;

  /** Vector representing the flux for each cmfd cell and cmfd enegy group at 
   * the end of a CMFD solve */
  Vector* _new_flux;

  /** Vector representing the flux for each cmfd cell and cmfd enegy group at 
   * the beginning of a CMFD solve */
  Vector* _old_flux;

  /** Vector representing the ratio of the new to old CMFD flux */
  Vector* _flux_ratio;

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

  /** Vector of CMFD cell volumes */
  Vector* _volumes;

  /** Array of material pointers for CMFD cell materials */
  Material** _materials;

  /** Physical dimensions of the geometry and each CMFD cell */
  double _width;
  double _height;
  double _cell_width;
  double _cell_height;

  /** Array of geometry boundaries */
  boundaryType* _boundaries;

  /** Vector of surface currents for each CMFD cell */
  Vector* _surface_currents;

  /** Vector of corner currents for each CMFD cell */
  Vector* _corner_currents;

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

  /** Flag indicating whether to use centroid updating (default true) */
  bool _centroid_update_on;

  /** Number of cells used in updating MOC flux (default 3) */
  int _k_nearest;

  /** Map storing the k-nearest stencil for each fsr */
  std::map<int, std::vector< std::pair<int, FP_PRECISION> > >
    _k_nearest_stencils;

  /* Private worker functions */
  FP_PRECISION computeDiffCorrect(FP_PRECISION d, FP_PRECISION h);
  void constructMatrices();
  void computeDs(int moc_iteration);
  void computeXS();
  void updateMOCFlux();
  void rescaleFlux();
  void splitCorners();
  void initializeMaterials();
  void initializeCurrents();
  void generateKNearestStencils();

  /* Private getter functions */
  int getCellNext(int cell_num, int surface_id);
  FP_PRECISION getUpdateRatio(int cmfd_cell, int moc_group, int fsr);
  FP_PRECISION getDistanceToCentroid(Point* centroid, int cell,
                                     int stencil_index);

public:

  Cmfd();
  virtual ~Cmfd();

  /* Worker functions */
  FP_PRECISION computeKeff(int moc_iteration);
  void initialize();
  void initializeCellMap();
  void initializeGroupMap();
  int findCmfdCell(LocalCoords* coords);
  int findCmfdSurface(int cell, LocalCoords* coords);
  int findCmfdCorner(int cell, LocalCoords* coords);
  void addFSRToCell(int cmfd_cell, int fsr_id);
  void zeroCurrents();
  void tallyCurrent(segment* curr_segment, FP_PRECISION* track_flux,
                    FP_PRECISION* polar_weights, bool fwd);
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
  std::vector< std::vector<int> >* getCellFSRs();
  bool isFluxUpdateOn();
  bool isCentroidUpdateOn();

  /* Set parameters */
  void setSORRelaxationFactor(FP_PRECISION SOR_factor);
  void setGeometry(Geometry* geometry);
  void setWidth(double width);
  void setHeight(double height);
  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumFSRs(int num_fsrs);
  void setNumMOCGroups(int num_moc_groups);
  void setOpticallyThick(bool thick);
  void setMOCRelaxationFactor(FP_PRECISION relax_factor);
  void setBoundary(int side, boundaryType boundary);
  void setLattice(Lattice* lattice);
  void setLatticeStructure(int num_x, int num_y);
  void setFluxUpdateOn(bool flux_update_on);
  void setCentroidUpdateOn(bool centroid_update_on);
  void setGroupStructure(int* group_indices, int length_group_indices);
  void setSourceConvergenceThreshold(FP_PRECISION source_thresh);
  void setPolarQuadrature(PolarQuad* polar_quad);
  void setKNearest(int k_nearest);

  /* Set FSR parameters */
  void setFSRMaterials(Material** FSR_materials);
  void setFSRVolumes(FP_PRECISION* FSR_volumes);
  void setFSRFluxes(FP_PRECISION* scalar_flux);
  void setCellFSRs(std::vector< std::vector<int> >* cell_fsrs);
};

#endif /* CMFD_H_ */
