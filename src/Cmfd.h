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
#include "Universe.h"
#include "Track2D.h"
#include "Track3D.h"
#include "Quadrature.h"
#include "linalg.h"
#include "Geometry.h"
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
#define track_flux(p,e) (track_flux[(p)*_num_moc_groups + (e)]

/**
 * @class Cmfd Cmfd.h "src/Cmfd.h"
 * @brief A class for Coarse Mesh Finite Difference (CMFD) acceleration.
 */
class Cmfd {

private:

  /** Pointer to polar quadrature object */
  Quadrature* _quadrature;

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

  /** Number of cells in z direction */
  int _num_z;

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

  /** If the user specified fine-to-coarse group indices */
  bool _user_group_indices;

  /** Number of FSRs */
  int _num_FSRs;

  /** The volumes (areas) for each FSR */
  FP_PRECISION* _FSR_volumes;

  /** Pointers to Materials for each FSR */
  Material** _FSR_materials;

  /** The FSR scalar flux in each energy group */
  FP_PRECISION* _FSR_fluxes;

  /** Array of CMFD cell volumes */
  Vector* _volumes;

  /** Array of material pointers for CMFD cell materials */
  Material** _materials;

  /** Physical dimensions of the geometry and each CMFD cell */
  FP_PRECISION _width_x;
  FP_PRECISION _width_y;
  FP_PRECISION _width_z;
  FP_PRECISION _cell_width_x;
  FP_PRECISION _cell_width_y;
  FP_PRECISION _cell_width_z;

  /** Array of geometry boundaries */
  boundaryType* _boundaries;

  /** Array of surface currents for each CMFD cell */
  Vector* _surface_currents;

  /** Vector of vectors of FSRs containing in each cell */
  std::vector< std::vector<int> > _cell_fsrs;

  /** Pointer to Lattice object representing the CMFD mesh */
  Lattice* _lattice;

  /** Flag indicating whether to update the MOC flux */
  bool _flux_update_on;

  /** Flag indicating whether to us centroid updating */
  bool _centroid_update_on;

  /** Number of cells to used in updating MOC flux */
  int _k_nearest;

  /** Map storing the k-nearest stencil for each fsr */
  std::map<int, std::vector< std::pair<int, FP_PRECISION> > >
    _k_nearest_stencils;

  /** OpenMP mutual exclusion locks for atomic CMFD cell operations */
  omp_lock_t* _cell_locks;

  /** Flag indicating whether the problem is 2D or 3D */
  bool _solve_3D;

  /** Array of azimuthal track spacings */
  FP_PRECISION* _azim_spacings;

  /** 2D array of polar track spacings */
  FP_PRECISION** _polar_spacings;

  /* Private worker functions */
  FP_PRECISION computeLarsensEDCFactor(FP_PRECISION dif_coef,
                                       FP_PRECISION delta);
  void constructMatrices(int moc_iteration);
  void collapseXS();
  void updateMOCFlux();
  void rescaleFlux();
  void splitVertexCurrents();
  void splitEdgeCurrents();
  void getVertexSplitSurfaces(int cell, int vertex, std::vector<int>* surfaces);
  void getEdgeSplitSurfaces(int cell, int edge, std::vector<int>* surfaces);
  void initializeMaterials();
  void initializeCurrents();
  void generateKNearestStencils();

  /* Private getter functions */
  int getCellNext(int cell_id, int surface_id);
  int getCellByStencil(int cell_id, int stencil_id);
  FP_PRECISION getUpdateRatio(int cell_id, int moc_group, int fsr);
  FP_PRECISION getDistanceToCentroid(Point* centroid, int cell_id,
                                     int stencil_index);
  FP_PRECISION getSurfaceDiffusionCoefficient(int cmfd_cell, int surface,
                                              int group, int moc_iteration,
                                              bool correction);
  FP_PRECISION getDiffusionCoefficient(int cmfd_cell, int group);
  FP_PRECISION getSurfaceWidth(int surface);
  FP_PRECISION getPerpendicularSurfaceWidth(int surface);
  int getSense(int surface);


public:

  Cmfd();
  virtual ~Cmfd();

  /* Worker functions */
  FP_PRECISION computeKeff(int moc_iteration);
  void initialize();
  void initializeCellMap();
  void initializeGroupMap();
  void initializeLattice(Point* offset);
  int findCmfdCell(LocalCoords* coords);
  int findCmfdSurface(int cell_id, LocalCoords* coords);
  int findCmfdSurfaceOTF(int cell_id, double z, int surface_2D);
  void addFSRToCell(int cell_id, int fsr_id);
  void zeroCurrents();
  void tallyCurrent(segment* curr_segment, FP_PRECISION* track_flux,
                    int azim_index, int polar_index, bool fwd);

  /* Get parameters */
  int getNumCmfdGroups();
  int getNumMOCGroups();
  int getNumCells();
  int getCmfdGroup(int group);
  int getBoundary(int side);
  Lattice* getLattice();
  int getNumX();
  int getNumY();
  int getNumZ();
  int convertFSRIdToCmfdCell(int fsr_id);
  std::vector< std::vector<int> >* getCellFSRs();
  bool isFluxUpdateOn();
  bool isCentroidUpdateOn();

  /* Set parameters */
  void setSORRelaxationFactor(FP_PRECISION SOR_factor);
  void setGeometry(Geometry* geometry);
  void setWidthX(double width);
  void setWidthY(double width);
  void setWidthZ(double width);
  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumZ(int num_z);
  void setNumFSRs(int num_fsrs);
  void setNumMOCGroups(int num_moc_groups);
  void setBoundary(int side, boundaryType boundary);
  void setLatticeStructure(int num_x, int num_y, int num_z=1);
  void setFluxUpdateOn(bool flux_update_on);
  void setCentroidUpdateOn(bool centroid_update_on);
  void setGroupStructure(int* group_indices, int length_group_indices);
  void setSourceConvergenceThreshold(FP_PRECISION source_thresh);
  void setQuadrature(Quadrature* quadrature);
  void setKNearest(int k_nearest);
  void setSolve3D(bool solve_3d);
  void setAzimSpacings(double* azim_spacings, int num_azim);
  void setPolarSpacings(double** polar_spacings, int num_azim,
                        int num_polar);

  /* Set FSR parameters */
  void setFSRMaterials(Material** FSR_materials);
  void setFSRVolumes(FP_PRECISION* FSR_volumes);
  void setFSRFluxes(FP_PRECISION* scalar_flux);
  void setCellFSRs(std::vector< std::vector<int> >* cell_fsrs);
};

#endif /* CMFD_H_ */
