/**
 * @file GPUSolver.h
 * @brief The GPUSolver class and CUDA physics kernels.
 * @date August 5, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef GPUSOLVER_H_
#define GPUSOLVER_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "../../Solver.h"
#endif

#define PySys_WriteStdout printf

#include <thrust/copy.h>
#include <iostream>
#include <vector>

#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/reduce.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include "clone.h"
#include "dev_exponential.h"
#include "GPUQuery.h"

/** If number of groups is known at compile time */
#ifdef NGROUPS
#define NUM_GROUPS (NGROUPS)
#define _NUM_GROUPS (NGROUPS)
#else
#define _NUM_GROUPS (_num_groups)
#endif

/** Indexing macro for the scalar flux in each FSR and energy group */
#define scalar_flux(tid,e) (scalar_flux[(tid)*NUM_GROUPS + (e)])

/** Indexing macro for the old scalar flux in each FSR and energy group */
#define old_scalar_flux(tid,e) (old_scalar_flux[(tid)*NUM_GROUPS + (e)])

/** Indexing macro for the total source divided by the total cross-section,
 *  \f$ \frac{Q}{\Sigma_t} \f$, in each FSR and energy group */
#define reduced_sources(tid,e) (reduced_sources[(tid)*NUM_GROUPS + (e)])

/** Indexing scheme for fixed sources for each FSR and energy group */
#define fixed_sources(r,e) (fixed_sources[(r)*NUM_GROUPS + (e)])

/** Indexing macro for the azimuthal and polar weights */
#define weights(i,p) (weights[(i)*num_polar_2 + (p)])

/** Indexing macro for the angular fluxes for each polar angle and energy
 *  group for a given Track */
#define boundary_flux(t,pe2) (boundary_flux[2*(t)*polar_times_groups+(pe2)])

/** Indexing macro for the starting angular fluxes for each polar angle and
 *  energy group for a given Track. These are copied to the boundary fluxes
 *  array at the beginning of each transport sweep */
#define start_flux(t,pe2) (start_flux[2*(t)*polar_times_groups+(pe2)])

/**
 * @class GPUSolver GPUSolver.h "openmoc/src/dev/gpu/GPUSolver.h"
 * @brief This a subclass of the Solver class for NVIDIA Graphics
 *        Processing Units (GPUs).
 * @details The source code for this class includes C++ coupled with
 *          compute intensive CUDA kernels for execution on the GPU.
 */
class GPUSolver : public Solver {

private:

  /** The number of thread blocks */
  int _B;

  /** The number of threads per thread block */
  int _T;

  /** The FSR Material pointers index by FSR ID */
  int* _FSR_materials;

  /** A pointer to an array of the Materials on the device */
  dev_material* _materials;

  /** Pointer to chi spectrum material on the device */
  dev_material* _dev_chi_spectrum_material;

  /** A pointer to the array of Tracks on the device */
  dev_track* _dev_tracks;

  /** Thrust vector of angular fluxes for each track */
  thrust::device_vector<FP_PRECISION> _boundary_flux;

  /** Thrust vector of starting angular fluxes for each track */
  thrust::device_vector<FP_PRECISION> _start_flux;

  /** Thrust vector of FSR scalar fluxes */
  thrust::device_vector<FP_PRECISION> _scalar_flux;

  /** Thrust vector of old FSR scalar fluxes */
  thrust::device_vector<FP_PRECISION> _old_scalar_flux;

  /** Thrust vector of stabilizing flux */
  thrust::device_vector<FP_PRECISION> _stabilizing_flux;

  /** Thrust vector of fixed sources in each FSR */
  thrust::device_vector<FP_PRECISION> _fixed_sources;

  /** Thrust vector of source / sigma_t in each FSR */
  thrust::device_vector<FP_PRECISION> _reduced_sources;

  /** Map of Material IDs to indices in _materials array */
  std::map<int, int> _material_IDs_to_indices;

  void copyQuadrature();

public:

  GPUSolver(TrackGenerator* track_generator=NULL);
  virtual ~GPUSolver();

  int getNumThreadBlocks();

  /**
   * @brief Returns the number of threads per block to execute on the GPU.
   * @return the number of threads per block
   */
  int getNumThreadsPerBlock();
  double getFSRSource(long fsr_id, int group) override;
  double getFlux(long fsr_id, int group) override;
  void getFluxes(FP_PRECISION* out_fluxes, int num_fluxes) override;

  void setNumThreadBlocks(int num_blocks);
  void setNumThreadsPerBlock(int num_threads);
  void setGeometry(Geometry* geometry) override;
  void setTrackGenerator(TrackGenerator* track_generator) override;
  void setFluxes(FP_PRECISION* in_fluxes, int num_fluxes) override;

  void initializeExpEvaluators() override;
  void initializeMaterials(solverMode mode) override;
  void initializeFSRs() override;
  void initializeTracks();
  void initializeFluxArrays() override;
  void initializeSourceArrays() override;
  void initializeFixedSources() override;
  void initializeCmfd() override;

  void zeroTrackFluxes();
  void flattenFSRFluxes(FP_PRECISION value);
  void flattenFSRFluxesChiSpectrum();
  void storeFSRFluxes();
  void computeStabilizingFlux();
  void stabilizeFlux();
  void computeFSRSources(int iteration);
  void computeFSRFissionSources();
  void computeFSRScatterSources();
  void transportSweep();
  void addSourceToScalarFlux();
  void computeKeff();
  double normalizeFluxes();
  double computeResidual(residualType res_type);

  void computeFSRFissionRates(double* fission_rates, long num_FSRs, bool nu = false);
};


#endif /* GPUSOLVER_H_ */
