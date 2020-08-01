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

/** Optimization macro for 3D calculations to avoid branch statements */
#ifdef THREED
#define _SOLVE_3D (true)
#endif

/** Forward declaration of Geometry class */
class Geometry;

/** Comparator for sorting k-nearest stencil std::pair objects */
inline bool stencilCompare(const std::pair<int, double>& firstElem,
                           const std::pair<int, double>& secondElem) {
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
  CMFD_PRECISION*** _boundary_volumes;
  CMFD_PRECISION*** _boundary_reaction;
  CMFD_PRECISION*** _boundary_diffusion;
  CMFD_PRECISION*** _old_boundary_flux;
  CMFD_PRECISION*** _boundary_surface_currents;

  CMFD_PRECISION*** _send_volumes;
  CMFD_PRECISION*** _send_reaction;
  CMFD_PRECISION*** _send_diffusion;
  CMFD_PRECISION*** _send_currents;

  CMFD_PRECISION* _send_split_current_data;
  CMFD_PRECISION* _receive_split_current_data;
  CMFD_PRECISION** _send_split_currents_array;
  CMFD_PRECISION** _receive_split_currents_array;
  CMFD_PRECISION*** _off_domain_split_currents;
  CMFD_PRECISION*** _received_split_currents;

  /** Vector representing the flux for each cmfd cell and cmfd energy group at
   * the end of a CMFD solve */
  Vector* _new_flux;

  /** Vector representing the flux for each cmfd cell and cmfd energy group at
   * the beginning of a CMFD solve */
  Vector* _old_flux;

  /** The corrected diffusion coefficients from the previous iteration */
  Vector* _old_dif_surf_corr;

  /** Whether the old diffusion coefficient has been set */
  bool _old_dif_surf_valid;

  /** Gauss-Seidel SOR relaxation factor */
  double _SOR_factor;

  /** cmfd source convergence threshold */
  double _source_convergence_threshold;

  /** Number of cells in x direction */
  int _num_x;

  /** Number of cells in y direction */
  int _num_y;

  /** Number of cells in z direction */
  int _num_z;

  /** Sweep number on MOC side */
  int _moc_iteration;

  /** Number of energy groups */
  int _num_moc_groups;

  /** Number of polar angles */
  int _num_polar;

  /** Number of azimuthal angles */
  int _num_azim;

  /** Number of energy groups used in cmfd solver. Note that cmfd supports
   * energy condensation from the MOC */
  int _num_cmfd_groups;

  /** Coarse energy indices for fine energy groups */
  int* _group_indices;

  /** Map of MOC groups to CMFD groups */
  int* _group_indices_map;

  /** Number of energy groups in the backup CMFD solver */
  int _num_backup_groups;

  /** Map of MOC groups to backup CMFD group structure */
  std::vector< std::vector<int> > _backup_group_structure;

  /** Map of CMFD groups to backup CMFD group structure */
  int* _cmfd_group_to_backup_group;

  /** If the user specified fine-to-coarse group indices */
  bool _user_group_indices;

  /** If a linear source approximation is used */
  bool _linear_source;

  /** If diffusion coefficients are limited by the flux */
  bool _flux_limiting;

  /** Whether to rebalance the computed sigma-t to be consistent with the MOC
   *  solution on every sweep */
  bool _balance_sigma_t;

  /** Number of FSRs */
  long _num_FSRs;

  /** The volumes (areas) for each FSR */
  FP_PRECISION* _FSR_volumes;

  /** Pointers to Materials for each FSR */
  Material** _FSR_materials;

  /** The FSR scalar flux in each energy group */
  FP_PRECISION* _FSR_fluxes;

  /** The FSR source in each energy group */
  FP_PRECISION* _FSR_sources;

  /** The source region flux moments (x, y, and z) for each energy group */
  FP_PRECISION* _flux_moments;

  /** Array of CMFD cell volumes */
  Vector* _volumes;

  /** Array of material pointers for CMFD cell materials */
  Material** _materials;

  /** Physical dimensions of the geometry and each CMFD cell */
  double _width_x;
  double _width_y;
  double _width_z;
  double _cell_width_x;
  double _cell_width_y;
  double _cell_width_z;

  /** Physical dimensions of non-uniform CMFD meshes (for whole geometry) */
  std::vector<double> _cell_widths_x;
  std::vector<double> _cell_widths_y;
  std::vector<double> _cell_widths_z;

  /** Distance of each mesh from the left-lower-bottom most point */
  std::vector<double> _accumulate_x;
  std::vector<double> _accumulate_y;
  std::vector<double> _accumulate_z;

  /** True if the cmfd meshes are non-uniform */
  bool _non_uniform;

  /** True if the cmfd mesh has been adjusted to fit the domain decomposition */
  bool _widths_adjusted_for_domains;

  /** Array of geometry boundaries */
  boundaryType* _boundaries;

  /** Array of surface currents for each CMFD cell */
  Vector* _surface_currents;

  /** Array of total current from starting boundary fluxes */
  Vector* _starting_currents;

  /** Array of net currents of all CMFD cells */
  Vector* _net_currents;

  /** Array of surface currents on all faces + edges and corners used in
      debugging */
  Vector* _full_surface_currents;

  /** Array of surface currents on edges and corners for each CMFD cell */
  std::map<int, CMFD_PRECISION> _edge_corner_currents;

  /** Vector of vectors of FSRs containing in each cell */
  std::vector< std::vector<long> > _cell_fsrs;

  /** Pointer to Lattice object representing the CMFD mesh */
  Lattice* _lattice;

  /** Flag indicating whether to update the MOC flux */
  bool _flux_update_on;

  /** Flag indicating whether to use centroid updating */
  bool _centroid_update_on;

  /** Flag indicating whether to check neutron balance on every CMFD solve */
  bool _check_neutron_balance;

  /** Flag indicating whether to print prolongation factors at every solve */
  bool _print_cmfd_prolongation_ratios;

  /** Whether to allow the CMFD solver to work with / return negative fluxes */
  bool _negative_fluxes_allowed;

  /** Number of MOC iterations before the CMFD update ratios are limited */
  int _num_unbounded_iterations;

  /** Number of cells to use in updating MOC flux */
  int _k_nearest;

  /** Relaxation factor to use for corrected diffusion coefficients */
  double _relaxation_factor;

  /** Map storing the k-nearest stencil for each fsr */
  std::map<long, std::vector< std::pair<int, double> > >
    _k_nearest_stencils;

  /** OpenMP mutual exclusion locks for atomic CMFD cell operations */
  omp_lock_t* _cell_locks;

  /** OpenMP mutual exclusion lock for edge/corner current tallies */
  omp_lock_t _edge_corner_lock;

#ifndef THREED
  /** Flag indicating whether the problem is 2D or 3D */
  bool _SOLVE_3D;
#endif

  /** Array of azimuthal track spacings */
  double* _azim_spacings;

  /** 2D array of polar track spacings */
  double** _polar_spacings;

  /** Whether to use axial interpolation for flux update ratios */
  int _use_axial_interpolation;

  /** Axial interpolation constants */
  std::vector<double*> _axial_interpolants;

  /* Structure to contain information about the convergence of the CMFD */
  ConvergenceData* _convergence_data;

  /* MPI communicator to transfer buffers, mainly currents at interfaces */
  DomainCommunicator* _domain_communicator;

  /* Buffer to contain received data */
  CMFD_PRECISION* _inter_domain_data;

  /* Buffer to contain sent data from domain */
  CMFD_PRECISION* _send_domain_data;

  /* For each face (1st dimension of the array), will contain data received */
  CMFD_PRECISION** _domain_data_by_surface;

  /* For each face (1st dimension of the array), will contain data to send */
  CMFD_PRECISION** _send_data_by_surface;

  /* Map of the indexes to each boundary in the tally arrays */
  std::vector<std::map<int, int> > _boundary_index_map;

  /* The number of on-domain cells in the x-direction */
  int _local_num_x;

  /* The number of on-domain cells in the y-direction */
  int _local_num_y;

  /* The number of on-domain cells in the z-direction */
  int _local_num_z;

  std::vector<int> _accumulate_lmx;
  std::vector<int> _accumulate_lmy;
  std::vector<int> _accumulate_lmz;

  /* Size of _tally_memory array */
  long _total_tally_size;

  /* 1D array that contains all tallies (diffusion, reaction and volume) */
  CMFD_PRECISION* _tally_memory;

  /* 2D array that contains reaction rates in each cell and group */
  CMFD_PRECISION** _reaction_tally;

  /* 2D array that contains volume tallies of each cell */
  CMFD_PRECISION** _volume_tally;

  /* 2D array that contains diffusion tallies for each cell and groups */
  CMFD_PRECISION** _diffusion_tally;

  /* Boolean to check if tallies are allocated */
  bool _tallies_allocated;

  /* Boolean to check if the domain communicator (for domain decomposed CMFD)
   * has been allocated */
  bool _domain_communicator_allocated;

  /** A timer to record timing data for a simulation */
  Timer* _timer;

  /** A one-group backup CMFD solver */
  Cmfd* _backup_cmfd;

  /* Private worker functions */
  CMFD_PRECISION computeLarsensEDCFactor(CMFD_PRECISION dif_coef,
                                         CMFD_PRECISION delta);
  void constructMatrices();
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
  CMFD_PRECISION getFluxRatio(int cell_id, int group, long fsr);
  CMFD_PRECISION getUpdateRatio(int cell_id, int moc_group, long fsr);
  double getDistanceToCentroid(Point* centroid, int cell_id, int local_cell_id,
                               int stencil_index);
  void getSurfaceDiffusionCoefficient(int cmfd_cell, int surface,
        int group, CMFD_PRECISION& dif_surf, CMFD_PRECISION& dif_surf_corr);
  CMFD_PRECISION getDiffusionCoefficient(int cmfd_cell, int group);
  CMFD_PRECISION getSurfaceWidth(int surface, int global_ind);
  CMFD_PRECISION getPerpendicularSurfaceWidth(int surface, int global_ind);
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
  void checkNeutronBalance(bool pre_split=true, bool old_source=false);
  void printProlongationFactors();

public:

  Cmfd();
  virtual ~Cmfd();

  /* Worker functions */
  double computeKeff(int moc_iteration);
  void initialize();
  void initializeCellMap();
  void initializeGroupMap();
  void allocateTallies();
  void initializeLattice(Point* offset, bool is_2D=false);
  void initializeBackupCmfdSolver();
  void copyCurrentsToBackup();
  int findCmfdCell(LocalCoords* coords);
  int findCmfdSurface(int cell_id, LocalCoords* coords);
  int findCmfdSurfaceOTF(int cell_id, double z, int surface_2D);
  void addFSRToCell(int cell_id, long fsr_id);
  void zeroCurrents();
  void tallyCurrent(segment* curr_segment, float* track_flux,
                    int azim_index, int polar_index, bool fwd);
  void tallyStartingCurrent(Point* point, double delta_x, double delta_y,
                            double delta_z, float* track_flux, double weight);
  void recordNetCurrents();

  /* Debug and information output */
  void printInputParamsSummary();
  void printTimerReport();
  void checkBalance();
  void printProlongation();

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
  int getLocalNumZ();
  Vector* getLocalCurrents();
  CMFD_PRECISION*** getBoundarySurfaceCurrents();
  int convertFSRIdToCmfdCell(long fsr_id);
  int convertGlobalFSRIdToCmfdCell(long global_fsr_id);
  std::vector< std::vector<long> >* getCellFSRs();
  bool isFluxUpdateOn();
  bool isCentroidUpdateOn();
  bool isSigmaTRebalanceOn();

  /* Set parameters */
  void setSORRelaxationFactor(double SOR_factor);
  void setCMFDRelaxationFactor(double relaxation_factor);
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
  void setSourceConvergenceThreshold(double source_thresh);
  void setQuadrature(Quadrature* quadrature);
  void setNumUnboundedIterations(int unbounded_iterations);
  void setKNearest(int k_nearest);
  void setSolve3D(bool solve_3d);
  void setAzimSpacings(const std::vector<double>& azim_spacings,
                       int num_azim);
  void setPolarSpacings(const std::vector< std::vector<double> >&
                        polar_spacings, int num_azim, int num_polar);
  void setKeff(double k_eff);
  void setBackupGroupStructure(std::vector< std::vector<int> > group_indices);

#ifdef MPIx
  void setNumDomains(int num_x, int num_y, int num_z);
  void setDomainIndexes(int idx_x, int idx_y, int idx_z);
#endif
  void setConvergenceData(ConvergenceData* convergence_data);
  void useAxialInterpolation(int interpolate);

  /* Methods to try to fix stability issues */
  void useFluxLimiting(bool flux_limiting);
  void enforceBalanceOnDiagonal(int cmfd_cell, int group);
  void rebalanceSigmaT(bool balance_sigma_t);

  /* Set FSR parameters */
  void setFSRMaterials(Material** FSR_materials);
  void setFSRVolumes(FP_PRECISION* FSR_volumes);
  void setFSRFluxes(FP_PRECISION* scalar_flux);
  void setFSRSources(FP_PRECISION* sources);
  void setCellFSRs(std::vector< std::vector<long> >* cell_fsrs);
  void setFluxMoments(FP_PRECISION* flux_moments);

  /* Set XYZ widths of non-uniform meshes */
  void setWidths(std::vector< std::vector<double> > widths);

  /* For debug use */
  void printCmfdCellSizes();

  /* For printing infomation about the CMFD object */
  std::string toString();
};


/**
 * @brief Get the CMFD group given an MOC group.
 * @param group the MOC energy group
 * @return the CMFD energy group
 */
inline int Cmfd::getCmfdGroup(int group) {
  return _group_indices_map[group];
}


/**
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
 * @brief Converts a local CMFD cell ID into its global ID
 * @param cmfd_cell The local CMFD cell ID
 * @return The global CMFD cell ID
 */
inline int Cmfd::getGlobalCMFDCell(int cmfd_cell) {

  int x_start = 0;
  int y_start = 0;
  int z_start = 0;
  if (_domain_communicator != NULL) {
    x_start = _accumulate_lmx[_domain_communicator->_domain_idx_x];
    y_start = _accumulate_lmy[_domain_communicator->_domain_idx_y];
    z_start = _accumulate_lmz[_domain_communicator->_domain_idx_z];
  }

  int ix = cmfd_cell % _local_num_x;
  int iy = (cmfd_cell % (_local_num_x * _local_num_y)) / _local_num_x;
  int iz = cmfd_cell / (_local_num_x * _local_num_y);

  return ((iz + z_start) * _num_y + iy + y_start) * _num_x
                + ix + x_start;
}


/**
 * @brief Tallies the current contribution from this segment across the
 *        the appropriate CMFD mesh cell surface.
 * @param curr_segment the current Track segment
 * @param track_flux the outgoing angular flux for this segment
 * @param azim_index azimuthal index of track angle
 * @param polar_index polar index of track angle
 * @param fwd boolean indicating direction of integration along segment
 */
inline void Cmfd::tallyCurrent(segment* curr_segment, float* track_flux,
                               int azim_index, int polar_index, bool fwd) {

  int surf_id, cell_id, cmfd_group;
  int ncg = _num_cmfd_groups;

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

    CMFD_PRECISION currents[_num_cmfd_groups]
         __attribute__ ((aligned(VEC_ALIGNMENT)));
    memset(currents, 0, _num_cmfd_groups * sizeof(CMFD_PRECISION));
    int local_cell_id = getLocalCMFDCell(cell_id);

    if (_SOLVE_3D) {
      double wgt = _quadrature->getWeightInline(azim_index, polar_index);
      for (int e=0; e < _num_moc_groups; e++) {

        /* Get the CMFD group */
        cmfd_group = getCmfdGroup(e);

        /* Increment the surface group current */
        currents[cmfd_group] += track_flux[e];
      }

#pragma omp simd aligned(currents)
      for (int g=0; g < ncg; g++)
        currents[g] *= wgt;

      /* Increment currents on faces */
      if (surf_id < NUM_FACES) {
        _surface_currents->incrementValues
            (local_cell_id, surf_id*ncg, (surf_id+1)*ncg - 1, currents);
      }
      /* Increment currents on corners and edges */
      else {

        int first_ind = (local_cell_id * NUM_SURFACES + surf_id) * ncg;
        omp_set_lock(&_edge_corner_lock);

#pragma omp simd aligned(currents)
        for (int g=0; g < ncg; g++)
          _edge_corner_currents[first_ind+g] += currents[g];

        omp_unset_lock(&_edge_corner_lock);
#ifdef INTEL
#pragma omp flush
#endif
      }
    }
    else {
      int pe = 0;
      for (int p=0; p < _num_polar/2; p++) {
        for (int e=0; e < _num_moc_groups; e++) {

          /* Get the CMFD group */
          cmfd_group = getCmfdGroup(e);

          currents[cmfd_group] += track_flux[pe]
              * _quadrature->getWeightInline(azim_index, p);
          pe++;
        }
      }

      /* Increment currents on face */
      if (surf_id < NUM_FACES) {
        _surface_currents->incrementValues
            (local_cell_id, surf_id*ncg, (surf_id+1)*ncg - 1, currents);
      }
      else {
        omp_set_lock(&_edge_corner_lock);

        int first_ind = (local_cell_id * NUM_SURFACES + surf_id) * ncg;

        /* Add contribution to corner current */
#pragma omp simd aligned(currents)
        for (int g=0; g < ncg; g++)
          _edge_corner_currents[first_ind+g] += currents[g];

        omp_unset_lock(&_edge_corner_lock);
#ifdef INTEL
#pragma omp flush
#endif
      }
    }
  }
}
#endif /* CMFD_H_ */
