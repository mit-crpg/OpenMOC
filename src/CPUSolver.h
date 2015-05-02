/**
 * @file CPUSolver.h
 * @brief The CPUSolver class.
 * @date May 28, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef CPUSOLVER_H_
#define CPUSOLVER_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include "Solver.h"
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#endif

/** Indexing macro for the angular fluxes for each polar angle and energy
 *  group for either the forward or reverse direction for a given Track */ 
#define track_flux(p,e) (track_flux[(p)*_num_groups + (e)])

/** Indexing macro for the angular fluxes for each polar angle and energy
 *  group for the outgoing reflective track from a given Track */
#define track_out_flux(p,e) (track_out_flux[(p)*_num_groups + (e)])

/** Indexing macro for the leakage for each polar angle and energy group
 *  for either the forward or reverse direction for a given Track */
#define track_leakage(p,e) (track_leakage[(p)*_num_groups + (e)])


/**
 * @class CPUSolver CPUSolver.h "src/CPUSolver.h"
 * @brief This a subclass of the Solver class for multi-core CPUs using
 *        OpenMP multi-threading.
 */
class CPUSolver : public Solver {

protected:

  /** The number of shared memory OpenMP threads */
  int _num_threads;

  /** OpenMP mutual exclusion locks for atomic FSR scalar flux updates */
  omp_lock_t* _FSR_locks;

  /** OpenMP mutual exclusion locks for atomic surface current updates */
  omp_lock_t* _cmfd_surface_locks;

  void initializeFSRs();
  void initializeFluxArrays();
  void initializeSourceArrays();
  void buildExpInterpTable();
  void initializeCmfd();

  void zeroTrackFluxes();
  void flattenFSRFluxes(FP_PRECISION value);
  void zeroSurfaceCurrents();
  void flattenFSRSources(FP_PRECISION value);
  void normalizeFluxes();
  FP_PRECISION computeFSRSources();

  /**
   * @brief Computes the contribution to the FSR flux from a Track segment.
   * @param curr_segment a pointer to the Track segment of interest
   * @param azim_index a pointer to the azimuthal angle index for this segment
   * @param track_flux a pointer to the Track's angular flux
   * @param fsr_flux a pointer to the temporary FSR scalar flux buffer
   * @param fwd
   */
  virtual void scalarFluxTally(segment* curr_segment, int azim_index,
                               FP_PRECISION* track_flux, FP_PRECISION* fsr_flux,
                               bool fwd);

  /**
   * @brief Updates the boundary flux for a Track given boundary conditions.
   * @param track_id the ID number for the Track of interest
   * @param azim_index a pointer to the azimuthal angle index for this segment
   * @param direction the Track direction (forward - true, reverse - false)
   * @param track_flux a pointer to the Track's outgoing angular flux
   */
  virtual void transferBoundaryFlux(int track_id, int azim_index,
                                    bool direction,
                                    FP_PRECISION* track_flux);
  void addSourceToScalarFlux();
  void computeKeff();
  void transportSweep();

  /**
   * @brief Computes the exponential term in the transport equation for a
   *        track segment.
   * @details This method uses either a linear interpolation table (default)
   *          or the exponential intrinsic exp(...) function if requested by
   *          the user through a call to the Solver::useExponentialIntrinsic()
   *          routine.
   * @param sigma_t the total group cross-section at this energy
   * @param length the length of the Track segment projected in the xy-plane
   * @param p the polar angle index
   * @return the evaluated exponential
   */
  virtual FP_PRECISION computeExponential(FP_PRECISION sigma_t,
                                          FP_PRECISION length, int p);

public:
  CPUSolver(Geometry* geometry=NULL, TrackGenerator* track_generator=NULL);
  virtual ~CPUSolver();

  int getNumThreads();
  FP_PRECISION getFSRScalarFlux(int fsr_id, int energy_group);
  FP_PRECISION* getFSRScalarFluxes();
  FP_PRECISION getFSRSource(int fsr_id, int energy_group);
  FP_PRECISION* getSurfaceCurrents();

  void setNumThreads(int num_threads);

  void computeFSRFissionRates(double* fission_rates, int num_FSRs);

};


#endif /* CPUSOLVER_H_ */
