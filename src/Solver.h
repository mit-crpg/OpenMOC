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
#include "Python.h"
#include "Timer.h"
#include "PolarQuad.h"
#include "TrackGenerator.h"
#include "Cmfd.h"
#include <math.h>
#endif

/** Indexing macro for the scalar flux in each FSR and energy group */
#define _scalar_flux(r,e) (_scalar_flux[(r)*_num_groups + (e)])

/** Indexing macro for the surface currents for each CMFD Mesh surface and
 *  each energy group */
#define _surface_currents(r,e) (_surface_currents[(r)*_cmfd->getNumCmfdGroups() \
                                                  + _cmfd->getCmfdGroup((e))])

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

/** Indexing macro for the leakage for each polar angle and energy group
 *  for both the forward and reverse direction for each track */
#define _boundary_leakage(i,pe2) (_boundary_leakage[2*(i)*_polar_times_groups \
                                                    +(pe2)])

/** Indexing scheme for the total fission source (\f$ \nu\Sigma_f\Phi \f$)
 *  for each FSR and energy group */
#define _fission_sources(r,e) (_fission_sources[(r)*_num_groups + (e)])

/** Indexing scheme for the total in-scatter source (\f$ \Sigma_s\Phi \f$)
 *  for each FSR and energy group */
#define _scatter_sources(r,e) (_scatter_sources[(r)*_num_groups + (e)])

/** The value of 4pi: \f$ 4\pi \f$ */
#define FOUR_PI 12.5663706143

/** The values of 1 divided by 4pi: \f$ \frac{1}{4\pi} \f$ */
#define ONE_OVER_FOUR_PI 0.0795774715


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

  /** The number of mesh cells */
  int _num_mesh_cells;

  /** The FSR "volumes" (i.e., areas) indexed by FSR UID */
  FP_PRECISION* _FSR_volumes;

  /** The FSR Material pointers indexed by FSR UID */
  Material** _FSR_materials;

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

  /** Twice the number of polar angles */
  int _two_times_num_polar;

  /** The number of polar angles times energy groups */
  int _polar_times_groups;

  /** A pointer to the 2D ragged array of Tracks */
  Track** _tracks;

  /** A pointer to an array with the number of Tracks per azimuthal angle */
  int* _num_tracks;

  /** The total number of Tracks */
  int _tot_num_tracks;

  /** The weights for each azimuthal angle */
  FP_PRECISION* _azim_weights;

  /** The weights for each polar angle in the polar angle quadrature */
  FP_PRECISION* _polar_weights;

  /** The angular fluxes for each Track for all energy groups, polar angles,
   *  and azimuthal angles. This array stores the boundary fluxes for a
   *  a Track along both "forward" and "reverse" directions. */
  FP_PRECISION* _boundary_flux;

  /** The angular leakages for each Track for all energy groups, polar angles,
   *  and azimuthal angles. This array stores the weighted outgoing fluxes
   *  for a Track along both "forward" and "reverse" directions. */
  FP_PRECISION* _boundary_leakage;

  /** The scalar flux for each energy group in each FSR */
  FP_PRECISION* _scalar_flux;

  /** The CMFD Mesh surface currents in each energy group */
  FP_PRECISION* _surface_currents;

  /** The fission source in each FSR and energy group */
  FP_PRECISION* _fission_sources;

  /** The in-scatter source in each FSR and energy group */
  FP_PRECISION* _scatter_sources;

  /** The old fission source in each FSR from the previous iteration */
  FP_PRECISION* _old_fission_sources;

  /** Ratios of source to total cross-section for each FSR and energy group */
  FP_PRECISION* _reduced_sources;

  /** An array of the residuals between the old source and the new source
   *  on each iteration in each FSR and energy group */
  FP_PRECISION* _source_residuals;

  /** The current iteration's approximation to k-effective */
  FP_PRECISION _k_eff;

  /** An array of k-effective at each iteration */
  std::vector<FP_PRECISION> _residual_vector;

  /** The total leakage across vacuum boundaries */
  FP_PRECISION _leakage;

  /** The number of source iterations needed to reach convergence */
  int _num_iterations;

  /** Whether or not the Solver has converged the source */
  bool _converged_source;

  /** The tolerance for converging the source */
  FP_PRECISION _source_convergence_thresh;

  /** A boolean indicating whether or not to use linear interpolation
   *  to comptue the exponential in the transport equation */
  bool _interpolate_exponential;

  /** The exponential linear interpolation table */
  FP_PRECISION* _exp_table;

  /** The size of the exponential linear interpolation table */
  int _exp_table_size;

  /** The maximum index of the exponential linear interpolation table */
  int _exp_table_max_index;

  /** The spacing for the exponential linear interpolation table */
  FP_PRECISION _exp_table_spacing;

  /** The inverse spacing for the exponential linear interpolation table */
  FP_PRECISION _inverse_exp_table_spacing;

  /** A timer to record timing data for a simulation */
  Timer* _timer;

  /** A pointer to a Coarse Mesh Finite Difference (CMFD) acceleration object */
  Cmfd* _cmfd;

  int round_to_int(float x);
  int round_to_int(double x);

  /**
   * @brief Initializes Track boundary angular flux and leakage and
   *        FSR scalar flux arrays.
   */
  virtual void initializeFluxArrays() =0;

  /**
   * @brief Allocates memory for FSR source arrays.
   */
  virtual void initializeSourceArrays() =0;

  /**
   * @brief Builds the exponential linear interpolation table.
   */
  virtual void buildExpInterpTable() =0;

  virtual void initializePolarQuadrature();
  virtual void initializeFSRs();
  virtual void initializeCmfd();
  virtual void checkTrackSpacing();

  /**
   * @brief Zero each Track's boundary fluxes for each energy group and polar
   *        angle in the "forward" and "reverse" directions.
   */
  virtual void zeroTrackFluxes() =0;

  /**
   * @brief Set the scalar flux for each FSR and energy group to some value.
   * @param value the value to assign to each FSR scalar flux
   */
  virtual void flattenFSRFluxes(FP_PRECISION value) =0;

  /**
   * @brief Set the source for each FSR and energy group to some value.
   * @param value the value to assign to each FSR source
   */
  virtual void flattenFSRSources(FP_PRECISION value) =0;

  /**
   * @brief Normalizes all FSR scalar fluxes and Track boundary angular
   *        fluxes to the total fission source (times \f$ \nu \f$).
   */
  virtual void normalizeFluxes() =0;

  /**
   * @brief Computes the total source (fission and scattering) for each FSR
   *        and energy group.
   * @return the residual between this source and the previous source
   */
  virtual FP_PRECISION computeFSRSources() =0;

  /**
   * @brief Compute \f$ k_{eff} \f$ from total fission and absorption rates
   *        in each FSR and energy group.
   */
  virtual void computeKeff() =0;

  /**
   * @brief Add the source term contribution in the transport equation to
   *        the FSR scalar flux.
   */
  virtual void addSourceToScalarFlux() =0;

  /**
   * @brief This method performs one transport sweep of all azimuthal angles,
   *        Tracks, segments, polar angles and energy groups.
   */
  virtual void transportSweep() =0;

  void clearTimerSplits();


public:
  Solver(Geometry* geom=NULL, TrackGenerator* track_generator=NULL);
  virtual ~Solver();

  Geometry* getGeometry();
  FP_PRECISION getFSRVolume(int fsr_id);
  TrackGenerator* getTrackGenerator();
  int getNumPolarAngles();
  int getNumIterations();
  double getTotalTime();
  FP_PRECISION getKeff();
  FP_PRECISION getSourceConvergenceThreshold();

  bool isUsingSinglePrecision();
  bool isUsingDoublePrecision();
  bool isUsingExponentialInterpolation();
  bool isUsingExponentialIntrinsic();

  /**
   * @brief Returns the scalar flux for a FSR and energy group.
   * @param fsr_id the ID for the FSR of interest
   * @param energy_group the energy group of interest
   * @return the FSR scalar flux
   */
  virtual FP_PRECISION getFSRScalarFlux(int fsr_id, int energy_group) =0;

  /**
   * @brief Returns an array of the scalar flux in each FSR and energy group.
   * @return an array of FSR scalar fluxes
   */
  virtual FP_PRECISION* getFSRScalarFluxes() =0;

  /**
   * @brief Returns the source for a FSR and energy group.
   * @param fsr_id the ID for the FSR of interest
   * @param energy_group the energy group of interest
   * @return the FSR source
   */
  virtual FP_PRECISION getFSRSource(int fsr_id, int energy_group) =0;

  virtual void setGeometry(Geometry* geometry);
  virtual void setTrackGenerator(TrackGenerator* track_generator);
  virtual void setPolarQuadrature(PolarQuad* polar_quad);
  virtual void setSourceConvergenceThreshold(FP_PRECISION source_thresh);

  void useExponentialInterpolation();
  void useExponentialIntrinsic();

  virtual FP_PRECISION convergeSource(int max_iterations);

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
  virtual void computeFSRFissionRates(double* fission_rates, int num_FSRs) =0;

  void printTimerReport();
};


/**
 * @brief Rounds a single precision floating point value to an integer.
 * @param x a float precision floating point value
 * @brief the rounded integer value
 */
inline int Solver::round_to_int(float x) {
  return lrintf(x);
}


/**
 * @brief Rounds a double precision floating point value to an integer.
 * @param x a double precision floating point value
 * @brief the rounded integer value
 */
inline int Solver::round_to_int(double x) {
  return lrint(x);
}


#endif /* SOLVER_H_ */
