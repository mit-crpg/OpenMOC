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
#include "constants.h"
#include "Universe.h"
#include "Track.h"
#include "Track3D.h"
#include "Quadrature.h"
#include "linalg.h"
#include "Geometry.h"
#include "Timer.h"
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
  double _k_eff;

  /** The A (destruction) matrix */
  Matrix* _A;

  /** The M (production) matrix */
  Matrix* _M;

  /** The old source vector */
  Vector* _old_source;

  /** The new source vector */
  Vector* _new_source;

  /* Domain boundary communication buffers */
  FP_PRECISION*** _boundary_volumes;
  FP_PRECISION*** _boundary_reaction;
  FP_PRECISION*** _boundary_diffusion;
  FP_PRECISION*** _old_boundary_flux;
  FP_PRECISION*** _boundary_surface_currents;

  FP_PRECISION*** _send_volumes;
  FP_PRECISION*** _send_reaction;
  FP_PRECISION*** _send_diffusion;
  FP_PRECISION*** _send_currents;

  FP_PRECISION* _send_split_current_data;
  FP_PRECISION* _receive_split_current_data;
  FP_PRECISION** _send_split_currents_array;
  FP_PRECISION** _receive_split_currents_array;
  FP_PRECISION*** _off_domain_split_currents;
  FP_PRECISION*** _received_split_currents;

  /** Vector representing the flux for each cmfd cell and cmfd enegy group at
   * the end of a CMFD solve */
  Vector* _new_flux;

  /** Vector representing the flux for each cmfd cell and cmfd enegy group at
   * the beginning of a CMFD solve */
  Vector* _old_flux;

  /** The corrected diffusion coefficients from the previous iteration */
  Vector* _old_dif_surf_corr;

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

  /** If a linear source approximation is used */
  bool _linear_source;

  /** If diffusion coefficients are limited by the flux */
  bool _flux_limiting;

  /** Number of FSRs */
  long _num_FSRs;

  /** The volumes (areas) for each FSR */
  FP_PRECISION* _FSR_volumes;

  /** Pointers to Materials for each FSR */
  Material** _FSR_materials;

  /** The FSR scalar flux in each energy group */
  NEW_FP_PRECISION* _FSR_fluxes;

  /** The source region flux moments (x, y, and z) for each energy group */
  FP_PRECISION* _flux_moments;

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

  /** Array of surface currents on all faces + edges and corners used in
      debugging */
  Vector* _full_surface_currents;

  /** Array of surface currents on edges and corners for each CMFD cell */
  std::map<int, FP_PRECISION> _edge_corner_currents;

  /** Vector of vectors of FSRs containing in each cell */
  std::vector< std::vector<long> > _cell_fsrs;

  /** Pointer to Lattice object representing the CMFD mesh */
  Lattice* _lattice;

  /** Flag indicating whether to update the MOC flux */
  bool _flux_update_on;

  /** Flag indicating whether to us centroid updating */
  bool _centroid_update_on;

  /** Flag indicating whether to check neutron balance on every CMFD solve */
  bool _check_neutron_balance;

  /** Number of cells to used in updating MOC flux */
  int _k_nearest;

  /** Relaxation factor to use for corrected diffusion coefficients */
  FP_PRECISION _relaxation_factor;

  /** Map storing the k-nearest stencil for each fsr */
  std::map<int, std::vector< std::pair<int, FP_PRECISION> > >
    _k_nearest_stencils;

  /** OpenMP mutual exclusion locks for atomic CMFD cell operations */
  omp_lock_t* _cell_locks;
  
  /** OpenMP mutual exclusion lock for edge/corner current tallies */
  omp_lock_t _edge_corner_lock;

  /** Flag indicating whether the problem is 2D or 3D */
  bool _solve_3D;

  /** Array of azimuthal track spacings */
  FP_PRECISION* _azim_spacings;

  /** 2D array of polar track spacings */
  FP_PRECISION** _polar_spacings;

  /** Whether to use axial interpolation for flux update ratios */
  bool _use_axial_interpolation;

  /** Axial interpolation constants */
  std::vector<double*> _axial_interpolants;

  //TODO: document
  ConvergenceData* _convergence_data;
  DomainCommunicator* _domain_communicator;
  FP_PRECISION* _inter_domain_data;
  FP_PRECISION* _send_domain_data;
  FP_PRECISION** _domain_data_by_surface;
  FP_PRECISION** _send_data_by_surface;
  std::vector<std::map<int, int> > _boundary_index_map;

  /* The number of on-domain cells in the x-direction */
  int _local_num_x;

  /* The number of on-domain cells in the y-direction */
  int _local_num_y;

  /* The number of on-domain cells in the z-direction */
  int _local_num_z;

  //TODO: document
  long _total_tally_size;
  FP_PRECISION* _tally_memory;
  FP_PRECISION** _reaction_tally;
  FP_PRECISION** _volume_tally;
  FP_PRECISION** _diffusion_tally;
  bool _tallies_allocated;

  /** A timer to record timing data for a simulation */
  Timer* _timer;

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
  int convertDirectionToSurface(int* direction);
  void convertSurfaceToDirection(int surface, int* direction);
  std::string getSurfaceNameFromDirection(int* direction);
  std::string getSurfaceNameFromSurface(int surface);

  /* Private getter functions */
  int getCellNext(int cell_id, int surface_id, bool global=true,
                  bool neighbor=false);
  int getCellByStencil(int cell_id, int stencil_id);
  FP_PRECISION getFluxRatio(int cell_id, int group, int fsr);
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
  int getLocalCMFDCell(int cmfd_cell); //TODO: optimize, document
  int getGlobalCMFDCell(int cmfd_cell); //TODO: optimize, document
    int getCellColor(int cmfd_cell); //TODO: optimize, document
  void packBuffers();
#ifdef MPIx
  void ghostCellExchange();
  void communicateSplits(bool faces);
#endif
  void unpackSplitCurrents(bool faces);
  void copyFullSurfaceCurrents();
  void checkNeutronBalance(bool pre_split=true);

public:

  Cmfd();
  virtual ~Cmfd();

  /* Worker functions */
  FP_PRECISION computeKeff(int moc_iteration);
  void initialize();
  void initializeCellMap();
  void initializeGroupMap();
  void allocateTallies();
  void initializeLattice(Point* offset);
  int findCmfdCell(LocalCoords* coords);
  int findCmfdSurface(int cell_id, LocalCoords* coords);
  int findCmfdSurfaceOTF(int cell_id, double z, int surface_2D);
  void addFSRToCell(int cell_id, long fsr_id);
  void zeroCurrents();
  void tallyCurrent(segment* curr_segment, float* track_flux,
                    int azim_index, int polar_index, bool fwd);
  void printTimerReport();
  void checkBalance();

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
  int convertFSRIdToCmfdCell(long fsr_id);
  int convertGlobalFSRIdToCmfdCell(long global_fsr_id);
  std::vector< std::vector<long> >* getCellFSRs();
  bool isFluxUpdateOn();
  bool isCentroidUpdateOn();

  /* Set parameters */
  void setSORRelaxationFactor(FP_PRECISION SOR_factor);
  void setCMFDRelaxationFactor(FP_PRECISION relaxation_factor);
  void setGeometry(Geometry* geometry);
  void setWidthX(double width);
  void setWidthY(double width);
  void setWidthZ(double width);
  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumZ(int num_z);
  void setNumFSRs(long num_fsrs);
  void setNumMOCGroups(int num_moc_groups);
  void setBoundary(int side, boundaryType boundary);
  void setLatticeStructure(int num_x, int num_y, int num_z=1);
  void setFluxUpdateOn(bool flux_update_on);
  void setCentroidUpdateOn(bool centroid_update_on);
  void setGroupStructure(std::vector< std::vector<int> > group_indices);
  void setSourceConvergenceThreshold(FP_PRECISION source_thresh);
  void setQuadrature(Quadrature* quadrature);
  void setKNearest(int k_nearest);
  void setSolve3D(bool solve_3d);
  void setAzimSpacings(FP_PRECISION* azim_spacings, int num_azim);
  void setPolarSpacings(FP_PRECISION** polar_spacings, int num_azim,
                        int num_polar);
  //TODO: clean, document
#ifdef MPIx
  void setNumDomains(int num_x, int num_y, int num_z);
  void setDomainIndexes(int idx_x, int idx_y, int idx_z);
#endif
  void setConvergenceData(ConvergenceData* convergence_data);
  void useAxialInterpolation(bool interpolate);
  void useFluxLimiting(bool flux_limiting);

  /* Set FSR parameters */
  void setFSRMaterials(Material** FSR_materials);
  void setFSRVolumes(FP_PRECISION* FSR_volumes);
  void setFSRFluxes(NEW_FP_PRECISION* scalar_flux);
  void setCellFSRs(std::vector< std::vector<long> >* cell_fsrs);
  void setFluxMoments(FP_PRECISION* flux_moments);
};


/**
 * @brief Get the CMFD group given an MOC group.
 * @param group the MOC energy group
 * @return the CMFD energy group
 */
inline int Cmfd::getCmfdGroup(int group) {
  return _group_indices_map[group];
}


/*
 * @brief Quickly finds a 3D CMFD surface given a cell, global coordinate, and
 *        2D CMFD surface. Intended for use in axial on-the-fly ray tracing.
 * @details If the coords is not on a surface, -1 is returned. If there is
 *          no 2D CMFD surface intersection, -1 should be input for the 2D CMFD
 *          surface.
 * @param cell_id The CMFD cell ID that the local coords is in.
 * @param z the axial height in the root universe of the point being evaluated.
 * @param surface_2D The ID of the 2D CMFD surface that the LocalCoords object
 *        intersects. If there is no 2D intersection, -1 should be input.
 */
inline int Cmfd::findCmfdSurfaceOTF(int cell_id, double z, int surface_2D) {
  int global_cell_id = getGlobalCMFDCell(cell_id);
  return _lattice->getLatticeSurfaceOTF(global_cell_id, z, surface_2D);
}


/**
 * @brief Tallies the current contribution from this segment across the
 *        the appropriate CMFD mesh cell surface.
 * @param curr_segment the current Track segment
 * @param track_flux the outgoing angular flux for this segment
 * @param polar_weights array of polar weights for some azimuthal angle
 * @param fwd boolean indicating direction of integration along segment
 */
inline void Cmfd::tallyCurrent(segment* curr_segment, float* track_flux,
                               int azim_index, int polar_index, bool fwd) {

  int surf_id, cell_id, cmfd_group;
  int ncg = _num_cmfd_groups;
  FP_PRECISION currents[_num_cmfd_groups];
  memset(currents, 0.0, sizeof(FP_PRECISION) * _num_cmfd_groups);
  std::map<int, FP_PRECISION>::iterator it;

  /* Check if the current needs to be tallied */
  bool tally_current = false;
  if (curr_segment->_cmfd_surface_fwd != -1 && fwd) {
    surf_id = curr_segment->_cmfd_surface_fwd % NUM_SURFACES;
    cell_id = curr_segment->_cmfd_surface_fwd / NUM_SURFACES;
    tally_current = true;
  }
  else if (curr_segment->_cmfd_surface_bwd != -1 && !fwd) {
    surf_id = curr_segment->_cmfd_surface_bwd % NUM_SURFACES;
    cell_id = curr_segment->_cmfd_surface_bwd / NUM_SURFACES;
    tally_current = true;
  }


  /* Tally current if necessary */
  if (tally_current) {

    int local_cell_id = getLocalCMFDCell(cell_id);

    if (_solve_3D) {
      FP_PRECISION wgt = _quadrature->getWeightInline(azim_index, polar_index);
      for (int e=0; e < _num_moc_groups; e++) {

        /* Get the CMFD group */
        cmfd_group = getCmfdGroup(e);

        /* Increment the surface group */
        currents[cmfd_group] += track_flux[e] * wgt;
      }

      /* Increment currents */
      if (surf_id < NUM_FACES) {
        _surface_currents->incrementValues
            (local_cell_id, surf_id*ncg, (surf_id+1)*ncg - 1, currents);
      }
      else {

        omp_set_lock(&_edge_corner_lock);

        int first_ind = (local_cell_id * NUM_SURFACES + surf_id) * ncg;
        it = _edge_corner_currents.find(first_ind);
        if (it == _edge_corner_currents.end())
          for (int g=0; g < ncg; g++)
            _edge_corner_currents[first_ind+g] = 0.0;

        for (int g=0; g < ncg; g++)
          _edge_corner_currents[first_ind+g] += currents[g];

        omp_unset_lock(&_edge_corner_lock);
      }
    }
    else {
      int pe = 0;
      for (int e=0; e < _num_moc_groups; e++) {

        /* Get the CMFD group */
        cmfd_group = getCmfdGroup(e);

        for (int p = 0; p < _num_polar/2; p++) {
          currents[cmfd_group] += track_flux[pe]
              * _quadrature->getWeightInline(azim_index, p);
          pe++;
        }
      }

      /* Increment currents */
      if (surf_id < NUM_FACES) {
        _surface_currents->incrementValues
            (local_cell_id, surf_id*ncg, (surf_id+1)*ncg - 1, currents);
      }
      else {
        omp_set_lock(&_edge_corner_lock);

        int first_ind = (local_cell_id * NUM_SURFACES + surf_id) * ncg;
        it = _edge_corner_currents.find(first_ind);
        if (it == _edge_corner_currents.end())
          for (int g=0; g < ncg; g++)
            _edge_corner_currents[first_ind+g] = 0.0;

        for (int g=0; g < ncg; g++)
          _edge_corner_currents[first_ind+g] += currents[g];
        
        omp_unset_lock(&_edge_corner_lock);

      }
    }
  }
}
#endif /* CMFD_H_ */
