/**
 * @file Solver.h
 * @brief The Solver class.
 * @date February 7, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#ifdef SWIG
#include "Python.h"
#endif
#include "constants.h"
#include "Timer.h"
#include "PolarQuad.h"
#include "TrackGenerator.h"
#include "Cmfd.h"
#include "ExpEvaluator.h"
#include <math.h>
#endif

/** Indexing macro for the scalar flux in each FSR and energy group */
#define _scalar_flux(r,e) (_scalar_flux[(r)*_num_groups + (e)])

/** Indexing macro for the old scalar flux in each FSR and energy group */
#define _old_scalar_flux(r,e) (_old_scalar_flux[(r)*_num_groups + (e)])

/** Indexing macro for the total source divided by the total cross-section
 *  (\f$ \frac{Q}{\Sigma_t} \f$) in each FSR and energy group */
#define _reduced_sources(r,e) (_reduced_sources[(r)*_num_groups + (e)])

/** Indexing macro for the polar quadrature weights multiplied by the
 *  azimuthal angle quadrature weights */
#define _polar_weights(i,p) (_polar_weights[(i)*_num_polar + (p)])

/** Indexing macro for the angular fluxes for each polar angle and energy
 *  group for the outgoing reflective track for both the forward and
 *  reverse direction for a given track */
#define _boundary_flux(i,j,p,e) (_boundary_flux[(i)*2*_polar_times_groups \
                                                + (j)*_polar_times_groups \
                                                + (p)*_num_groups + (e)])

/** Indexing macro for the angular fluxes for each polar angle and energy
 *  group for the outgoing reflective track for both the forward and
 *  reverse direction for a given track */
#define _boundary_flux_copy(i,j,p,e) (_boundary_flux_copy[(i)*2*_polar_times_groups \
                                                          + (j)*_polar_times_groups \
                                                          + (p)*_num_groups + (e)])

/** Indexing scheme for fixed sources for each FSR and energy group */
#define _fixed_sources(r,e) (_fixed_sources[(r)*_num_groups + (e)])

/** Indexing scheme for the total fission source (\f$ \nu\Sigma_f\Phi \f$)
 *  for each FSR and energy group */
#define fission_sources(r,e) (fission_sources[(r)*_num_groups + (e)])

/** Indexing scheme for the total scatter source (\f$ Sigma_s\Phi \f$)
 *  for each FSR and energy group */
#define scatter_sources(r,e) (scatter_sources[(r)*_num_groups + (e)])


/**
 * @enum solverMode
 * @brief The solution mode used by the MOC solver.
*/
enum solverMode {

  /** The forward flux distribution */
  FORWARD,

  /** The adjoint flux distribution */
  ADJOINT,
};


/**
 * @enum residualType
 * @brief The type of residual used for the convergence criterion.
*/
enum residualType {

  /** A residual on the scalar flux distribution */
  SCALAR_FLUX,

  /** A residual on the fission source distribution */
  FISSION_SOURCE,

  /** A residual on the total source distribution */
  TOTAL_SOURCE,
};


/**
 * @class Solver Solver.h "src/Solver.h"
 * @brief This is an abstract base class which different Solver subclasses
 *        implement for different architectures or source iteration algorithms.
 */
class Solver {

protected:

  /** The number of azimuthal angles */
  int _num_azim;

  /** The number of energy groups */
  int _num_groups;

  /** The number of flat source regions */
  int _num_FSRs;

  /** The number of fissionable flat source regions */
  int _num_fissionable_FSRs;

  /** The FSR "volumes" (i.e., areas) indexed by FSR UID */
  FP_PRECISION* _FSR_volumes;

  /** The FSR Material pointers indexed by FSR UID */
  Material** _FSR_materials;

  /** The FSR Centroid pointers indexed by FSR UID */
  Point** _FSR_centroids;

  /** A pointer to a TrackGenerator which contains Tracks */
  TrackGenerator* _track_generator;

  /** A pointer to a Geometry with initialized FSR offset maps */
  Geometry* _geometry;

  /** The number of Materials */
  int _num_materials;

  /** A pointer to a polar quadrature */
  PolarQuad* _polar_quad;

  /** A boolean indicating if a user-defined PolarQuad was assigned */
  bool _user_polar_quad;

  /** The number of polar angles */
  int _num_polar;

  /** The number of polar angles times energy groups */
  int _polar_times_groups;

  /** A pointer to the 2D ragged array of Tracks */
  Track** _tracks;

  /** The total number of Tracks */
  int _tot_num_tracks;

  /** The weights for each polar angle in the polar angle quadrature */
  FP_PRECISION* _polar_weights;

  /** The angular fluxes for each Track for all energy groups, polar angles,
   *  and azimuthal angles. This array stores the boundary fluxes for a
   *  a Track along both "forward" and "reverse" directions. */
  FP_PRECISION* _boundary_flux;
  FP_PRECISION* _boundary_flux_copy;

  /** The scalar flux for each energy group in each FSR */
  FP_PRECISION* _scalar_flux;

  /** The old scalar flux for each energy group in each FSR */
  FP_PRECISION* _old_scalar_flux;

  /** Optional user-specified fixed sources in each FSR and energy group */
  FP_PRECISION* _fixed_sources;

  /** A mapping of fixed sources keyed by the pair (FSR ID, energy group) */
  std::map< std::pair<int, int>, FP_PRECISION > _fix_src_FSR_map;

  /** A mapping of fixed sources keyed by the pair (Cell*, energy group) */
  std::map< std::pair<Cell*, int>, FP_PRECISION > _fix_src_cell_map;

  /** A mapping of fixed sources keyed by the pair (Material*, energy group) */
  std::map< std::pair<Material*, int>, FP_PRECISION > _fix_src_material_map;

  /** Ratios of source to total cross-section for each FSR and energy group */
  FP_PRECISION* _reduced_sources;

  /** The current iteration's approximation to k-effective */
  FP_PRECISION _k_eff;

  /** The number of source iterations needed to reach convergence */
  int _num_iterations;

  /** The tolerance for converging the source/flux */
  FP_PRECISION _converge_thresh;

  /** An ExpEvaluator to compute exponentials in the transport equation */
  ExpEvaluator* _exp_evaluator;

  /** Indicator of whether the flux array is defined by the user */
  bool _user_fluxes;

  /** A timer to record timing data for a simulation */
  Timer* _timer;

  /** A pointer to a Coarse Mesh Finite Difference (CMFD) acceleration object */
  Cmfd* _cmfd;

  /** The number of groups of tracks that can be looped over in parallel
   *  without data races between threads */
  int _num_parallel_track_groups;

  void clearTimerSplits();

public:
  Solver(TrackGenerator* track_generator=NULL);
  virtual ~Solver();

  virtual void setGeometry(Geometry* geometry);

  Geometry* getGeometry();
  TrackGenerator* getTrackGenerator();
  PolarQuad* getPolarQuad();
  FP_PRECISION getFSRVolume(int fsr_id);
  int getNumPolarAngles();
  int getNumIterations();
  double getTotalTime();
  FP_PRECISION getKeff();
  FP_PRECISION getConvergenceThreshold();
  FP_PRECISION getMaxOpticalLength();
  bool isUsingDoublePrecision();
  bool isUsingExponentialInterpolation();

  virtual FP_PRECISION getFSRSource(int fsr_id, int group);
  virtual FP_PRECISION getFlux(int fsr_id, int group);
  virtual void getFluxes(FP_PRECISION* out_fluxes, int num_fluxes) = 0;

  virtual void setTrackGenerator(TrackGenerator* track_generator);
  virtual void setPolarQuadrature(PolarQuad* polar_quad);
  virtual void setConvergenceThreshold(FP_PRECISION threshold);
  virtual void setFluxes(FP_PRECISION* in_fluxes, int num_fluxes) = 0;
  void setFixedSourceByFSR(int fsr_id, int group, FP_PRECISION source);
  void setFixedSourceByCell(Cell* cell, int group, FP_PRECISION source);
  void setFixedSourceByMaterial(Material* material, int group,
                                FP_PRECISION source);
  void setMaxOpticalLength(FP_PRECISION max_optical_length);
  void setExpPrecision(FP_PRECISION precision);
  void useExponentialInterpolation();
  void useExponentialIntrinsic();

  virtual void initializePolarQuadrature();
  virtual void initializeExpEvaluator();
  virtual void initializeMaterials(solverMode mode=FORWARD);
  virtual void initializeFSRs();
  virtual void countFissionableFSRs();
  virtual void initializeFixedSources();
  virtual void initializeCmfd();

  virtual void resetMaterials(solverMode mode=FORWARD);
  virtual void fissionTransportSweep();
  virtual void scatterTransportSweep();

  /**
   * @brief Initializes Track boundary angular and FSR scalar flux arrays.
   */
  virtual void initializeFluxArrays() = 0;

  /**
   * @brief Allocates memory for FSR source arrays.
   */
  virtual void initializeSourceArrays() = 0;

  /**
   * @brief Zero each Track's boundary fluxes for each energy group and polar
   *        angle in the "forward" and "reverse" directions.
   */
  virtual void zeroTrackFluxes() = 0;

  /**
   * @brief Set the scalar flux for each FSR and energy group to some value.
   * @param value the value to assign to each FSR scalar flux
   */
  virtual void flattenFSRFluxes(FP_PRECISION value) = 0;

  /**
   * @brief Stores the current scalar fluxes in the old scalar flux array.
   */
  virtual void storeFSRFluxes() = 0;

  /**
   * @brief Normalizes all FSR scalar fluxes and Track boundary angular
   *        fluxes to the total fission source (times \f$ \nu \f$).
   */
  virtual void normalizeFluxes() = 0;

  /**
   * @brief Computes the total source (fission, scattering, fixed) for
   *        each FSR and energy group.
   */
  virtual void computeFSRSources() = 0;

  /**
   * @brief Computes the total fission source for each FSR and energy group.
   */
  virtual void computeFSRFissionSources() = 0;

  /**
   * @brief Computes the total scattering source for each FSR and energy group.
   */
  virtual void computeFSRScatterSources() = 0;

  /**
   * @brief Computes the residual between successive flux/source iterations.
   * @param res_type the type of residual (FLUX, FISSION_SOURCE, TOTAL_SOURCE)
   * @return the total residual summed over FSRs and energy groups
   */
  virtual double computeResidual(residualType res_type) = 0;

  /**
   * @brief Compute \f$ k_{eff} \f$ from total fission and absorption rates
   *        in each FSR and energy group.
   */
  virtual void computeKeff() = 0;

  /**
   * @brief Add the source term contribution in the transport equation to
   *        the FSR scalar flux.
   */
  virtual void addSourceToScalarFlux() = 0;

  /**
   * @brief This method performs one transport swep.
   */
  virtual void transportSweep() = 0;

  void computeFlux(int max_iters=1000, solverMode mode=FORWARD,
                   bool only_fixed_source=true);
  void computeSource(int max_iters=1000, solverMode mode=FORWARD,
                     double k_eff=1.0, residualType res_type=TOTAL_SOURCE);
  void computeEigenvalue(int max_iters=1000, solverMode mode=FORWARD,
                         residualType res_type=FISSION_SOURCE);

 /**
  * @brief Computes the volume-weighted, energy integrated fission rate in
  *        each FSR and stores them in an array indexed by FSR ID.
  * @details This is a helper method for SWIG to allow users to retrieve
  *          FSR fission rates as a NumPy array. An example of how this method
  *          can be called from Python is as follows:
  *
  * @code
  *          num_FSRs = geometry.getNumFSRs()
  *          fission_rates = solver.computeFSRFissionRates(num_FSRs)
  * @endcode
  *
  * @param fission_rates an array to store the fission rates (implicitly passed
  *                      in as a NumPy array from Python)
  * @param num_FSRs the number of FSRs passed in from Python
  */
  virtual void computeFSRFissionRates(double* fission_rates, int num_FSRs) = 0;

  void printTimerReport();
};


#endif /* SOLVER_H_ */
