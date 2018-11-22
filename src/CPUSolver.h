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


/* Structure containing the info to send about a track (used in printCycle) */
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

#ifdef MPIx
  /* Message size when communicating track angular fluxes at interfaces */
  int _track_message_size;

  /* Buffer to send track angular fluxes and associated information */
  std::vector<float*> _send_buffers;

  /* Buffer to receive track angular fluxes and associated information */
  std::vector<float*> _receive_buffers;

  /* Vector of vectors containing boundary track ids and direction */
  std::vector<std::vector<long> > _boundary_tracks;

  /* Vector of vectors containing the connecting track id and direction */
  std::vector<std::vector<long> > _track_connections;

  /* Rank of domains neighboring local domain */
  std::vector<int> _neighbor_domains;

  /* Array to check whether MPI communications are finished */
  MPI_Request* _MPI_requests;

  /* Arrays of booleans to know whether a send/receive call was made */
  bool* _MPI_sends;
  bool* _MPI_receives;
#endif


  virtual void initializeFluxArrays();
  virtual void initializeSourceArrays();
  virtual void initializeFSRs();


  void zeroTrackFluxes();
  void copyBoundaryFluxes();
  void tallyStartingCurrents();
#ifdef MPIx
  void setupMPIBuffers();
  void deleteMPIBuffers();
  void packBuffers(std::vector<long> &packing_indexes);
  void transferAllInterfaceFluxes();
  void printCycle(long track_start, int domain_start, int length);
  void boundaryFluxChecker();
#endif
  virtual void flattenFSRFluxes(FP_PRECISION value);
  void flattenFSRFluxesChiSpectrum();
  void storeFSRFluxes();
  virtual double normalizeFluxes();
  virtual void computeFSRSources(int iteration);
  void transportSweep();
  virtual void computeStabilizingFlux();
  virtual void stabilizeFlux();
  virtual void addSourceToScalarFlux();
  void computeKeff();
  double computeResidual(residualType res_type);

public:
  CPUSolver(TrackGenerator* track_generator=NULL);
  virtual ~CPUSolver();

  int getNumThreads();
  void setNumThreads(int num_threads);
  void setFixedSourceByFSR(long fsr_id, int group, FP_PRECISION source);
  void computeFSRFissionRates(double* fission_rates, long num_FSRs);
  void printInputParamsSummary();

  virtual void tallyScalarFlux(segment* curr_segment, int azim_index,
                               int polar_index, float* track_flux);

  void tallyCurrent(segment* curr_segment, int azim_index,
                            int polar_index, float* track_flux,
                            bool fwd);

  void transferBoundaryFlux(Track* track, int azim_index,
                                    int polar_index, bool direction,
                                    float* track_flux);

  void getFluxes(FP_PRECISION* out_fluxes, int num_fluxes);

  void initializeFixedSources();

  void printFSRFluxes(std::vector<double> dim1,
                      std::vector<double> dim2, double offset,
                      const char* plane);
  void printFluxesTemp();
  void printNegativeSources(int iteration, int num_x, int num_y, int num_z);
};


#endif /* CPUSOLVER_H_ */
