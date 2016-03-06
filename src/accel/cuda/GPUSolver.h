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
#include "../../constants.h"
#include "../../Solver.h"
#endif

#define PySys_WriteStdout printf

#include <thrust/copy.h>
#include <iostream>

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
#include <sm_20_atomic_functions.h>
#include "clone.h"
#include "GPUExpEvaluator.h"

/** Indexing macro for the scalar flux in each SR and energy group */
#define scalar_flux(tid,e) (scalar_flux[(tid)*(*num_groups) + (e)])

/** Indexing macro for the old scalar flux in each SR and energy group */
#define old_scalar_flux(tid,e) (old_scalar_flux[(tid)*(*num_groups) + (e)])

/** Indexing macro for the total source divided by the total cross-section,
 *  \f$ \frac{Q}{\Sigma_t} \f$, in each SR and energy group */
#define reduced_sources(tid,e) (reduced_sources[(tid)*(*num_groups) + (e)])

/** Indexing scheme for fixed sources for each SR and energy group */
#define fixed_sources(r,e) (fixed_sources[(r)*(*num_groups) + (e)])

/** Indexing macro for the azimuthal and polar weights */
#define polar_weights(i,p) (polar_weights[(i)*(*num_polar) + (p)])

/** Indexing macro for the angular fluxes for each polar angle and energy
 *  group for a given Track */
#define boundary_flux(t,pe2) (boundary_flux[2*(t)*(*polar_times_groups)+(pe2)])


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

  /** Twice the number of polar angles */
  int _two_times_num_polar;

  /** The SR Material pointers index by SR ID */
  int* _SR_materials;

  /** A pointer to an array of the Materials on the device */
  dev_material* _materials;

  /** A pointer to the array of Tracks on the device */
  dev_track* _dev_tracks;

  /** Thrust vector of angular fluxes for each track */
  thrust::device_vector<FP_PRECISION> _boundary_flux;

  /** Thrust vector of SR scalar fluxes */
  thrust::device_vector<FP_PRECISION> _scalar_flux;

  /** Thrust vector of old SR scalar fluxes */
  thrust::device_vector<FP_PRECISION> _old_scalar_flux;

  /** Thrust vector of fixed sources in each SR */
  thrust::device_vector<FP_PRECISION> _fixed_sources;

  /** Thrust vector of source / sigma_t in each SR */
  thrust::device_vector<FP_PRECISION> _reduced_sources;

  /** Map of Material IDs to indices in _materials array */
  std::map<int, int> _material_IDs_to_indices;

public:

  GPUSolver(TrackGenerator* track_generator=NULL);
  virtual ~GPUSolver();

  int getNumThreadBlocks();

  /**
   * @brief Returns the number of threads per block to execute on the GPU.
   * @return the number of threads per block
   */
  int getNumThreadsPerBlock();
  FP_PRECISION getSRSource(int sr_id, int group);
  FP_PRECISION getFlux(int sr_id, int group);
  void getFluxes(FP_PRECISION* out_fluxes, int num_fluxes);

  void setNumThreadBlocks(int num_blocks);
  void setNumThreadsPerBlock(int num_threads);
  void setGeometry(Geometry* geometry);
  void setTrackGenerator(TrackGenerator* track_generator);
  void setFluxes(FP_PRECISION* in_fluxes, int num_fluxes);

  void initializePolarQuadrature();
  void initializeExpEvaluator();
  void initializeMaterials(solverMode mode=ADJOINT);
  void initializeSRs();
  void initializeTracks();
  void initializeFluxArrays();
  void initializeSourceArrays();
  void initializeFixedSources();

  void zeroTrackFluxes();
  void flattenSRFluxes(FP_PRECISION value);
  void storeSRFluxes();
  void normalizeFluxes();
  void computeSRSources();
  void computeSRFissionSources();
  void computeSRScatterSources();
  void transportSweep();
  void addSourceToScalarFlux();
  void computeKeff();
  double computeResidual(residualType res_type);

  void computeSRFissionRates(double* fission_rates, int num_SRs);
};


#endif /* GPUSOLVER_H_ */
