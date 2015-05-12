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


#undef track_flux

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

  void initializeFluxArrays();
  void initializeSourceArrays();
  void initializeFSRs();

  void zeroTrackFluxes();
  void flattenFSRFluxes(FP_PRECISION value);
  void flattenFSRSources(FP_PRECISION value);
  void normalizeFluxes();
  FP_PRECISION computeFSRSources();

  /**
   * @brief Computes the contribution to the FSR flux from a Track segment.
   * @param curr_segment a pointer to the Track segment of interest
   * @param azim_index a pointer to the azimuthal angle index for this segment
   * @param track_flux a pointer to the Track's angular flux
   * @param fsr_flux a pointer to the temporary FSR scalar flux buffer
   */
  virtual void tallyScalarFlux(segment* curr_segment, int azim_index,
                               FP_PRECISION* track_flux, FP_PRECISION* fsr_flux);

  /**
   * @brief Computes the contribution to surface current from a Track segment.
   * @param curr_segment a pointer to the Track segment of interest
   * @param azim_index a pointer to the azimuthal angle index for this segment
   * @param track_flux a pointer to the Track's angular flux
   * @param fwd the direction of integration along the segment
   */
  virtual void tallySurfaceCurrent(segment* curr_segment, int azim_index,
                                   FP_PRECISION* track_flux, bool fwd);

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

public:
  CPUSolver(TrackGenerator* track_generator=NULL);
  virtual ~CPUSolver();

  int getNumThreads();
  FP_PRECISION getFSRScalarFlux(int fsr_id, int energy_group);
  FP_PRECISION* getFSRScalarFluxes();
  FP_PRECISION getFSRSource(int fsr_id, int energy_group);

  void setNumThreads(int num_threads);

  void computeFSRFissionRates(double* fission_rates, int num_FSRs);
};


#endif /* CPUSOLVER_H_ */
