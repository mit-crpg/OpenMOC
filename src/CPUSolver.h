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
#include "TrackTraversingAlgorithms.h"
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#endif

#undef track_flux

/** Indexing macro for the angular fluxes for each polar angle and energy
 *  group for either the forward or reverse direction for a given Track */
#define track_flux(pe) (track_flux[(pe)])

/** Indexing macro for the angular fluxes for each polar angle and energy
 *  group for the outgoing reflective track from a given Track */
#define track_out_flux(pe) (track_out_flux[(pe)])

/** Indexing macro for the leakage for each polar angle and energy group
 *  for either the forward or reverse direction for a given Track */
#define track_leakage(pe) (track_leakage[(pe)])


//FIXME
struct sendInfo {
  long track_id;
  int domain;
  bool fwd;
};


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

  //FIXME
#ifdef MPIx
  int _track_message_size;
  std::vector<FP_PRECISION*> _send_buffers;
  std::vector<FP_PRECISION*> _receive_buffers;
  std::vector<int> _neighbor_domains;
  MPI_Request* _MPI_requests;
  bool* _MPI_sends;
  bool* _MPI_receives;
#endif


  void initializeFluxArrays();
  void initializeSourceArrays();
  void initializeFSRs();

  void zeroTrackFluxes();
  void copyBoundaryFluxes();
#ifdef MPIx
  void setupMPIBuffers();
  void deleteMPIBuffers();
  void packBuffers(std::vector<long> &packing_indexes,
                   std::vector<int> &buffer_indexes);
  void transferAllInterfaceFluxesNew();
  void transferAllInterfaceFluxes();
  void printCycle(long track_start, int domain_start, int length);
#endif
  void flattenFSRFluxes(FP_PRECISION value);
  void storeFSRFluxes();
  void normalizeFluxes();
  void computeFSRSources();
  void transportSweep();
  void addSourceToScalarFlux();
  void computeKeff();
  double computeResidual(residualType res_type);

public:
  CPUSolver(TrackGenerator* track_generator=NULL);
  virtual ~CPUSolver();

  int getNumThreads();
  void setNumThreads(int num_threads);
  virtual void setFixedSourceByFSR(int fsr_id, int group, FP_PRECISION source);
  void computeFSRFissionRates(double* fission_rates, int num_FSRs);

  /**
   * @brief Computes the contribution to the FSR flux from a Track segment.
   * @param curr_segment a pointer to the Track segment of interest
   * @param azim_index a pointer to the azimuthal angle index for this segment
   * @param track_flux a pointer to the Track's angular flux
   * @param fsr_flux a pointer to the temporary FSR scalar flux buffer
   */
  virtual void tallyScalarFlux(segment* curr_segment, int azim_index,
                               int polar_index, FP_PRECISION* track_flux,
                               FP_PRECISION* fsr_flux);

  /**
   * @brief Computes the contribution to surface current from a Track segment.
   * @param curr_segment a pointer to the Track segment of interest
   * @param azim_index a pointer to the azimuthal angle index for this segment
   * @param track_flux a pointer to the Track's angular flux
   * @param fwd the direction of integration along the segment
   */
  virtual void tallyCurrent(segment* curr_segment, int azim_index,
                            int polar_index, FP_PRECISION* track_flux,
                            bool fwd);

  /**
   * @brief Updates the boundary flux for a Track given boundary conditions.
   * @param track_id the ID number for the Track of interest FIXME
   * @param azim_index a pointer to the azimuthal angle index for this segment
   * @param direction the Track direction (forward - true, reverse - false)
   * @param track_flux a pointer to the Track's outgoing angular flux
   */
  virtual void transferBoundaryFlux(Track* track, int azim_index,
                                    int polar_index, bool direction,
                                    FP_PRECISION* track_flux);

  virtual void getFluxes(FP_PRECISION* out_fluxes, int num_fluxes);
};


#endif /* CPUSOLVER_H_ */
