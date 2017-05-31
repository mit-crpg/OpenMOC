/**
 * @file GPUSolver.h
 * @brief The GPUSolver class and CUDA physics kernels.
 * @date August 5, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef GPUSOLVER_H_
#define GPUSOLVER_H_

#ifdef __cplusplus
#include "Python.h"
#include "../../constants.h"
#include "../../Solver.h"
#endif

#define PySys_WriteStdout printf

#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/reduce.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <thrust/iterator/constant_iterator.h>
#include <sm_20_atomic_functions.h>
#include "clone.h"
#include "GPUExpEvaluator.h"

/** Indexing macro for the scalar flux in each FSR and energy group */
#define scalar_flux(tid,e) (scalar_flux[(tid)*(*num_groups) + (e)])

/** Indexing macro for the old scalar flux in each FSR and energy group */
#define old_scalar_flux(tid,e) (old_scalar_flux[(tid)*(*num_groups) + (e)])

/** Indexing macro for the total source divided by the total cross-section,
 *  \f$ \frac{Q}{\Sigma_t} \f$, in each FSR and energy group */
#define reduced_sources(tid,e) (reduced_sources[(tid)*(*num_groups) + (e)])

/** Indexing scheme for fixed sources for each FSR and energy group */
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

  /** The FSR Material pointers index by FSR ID */
  int* _FSR_materials;

  /** A pointer to an array of the Materials on the device */
  dev_material* _materials;

  /** A pointer to the array of Tracks on the device */
  dev_track* _dev_tracks;

  /** Thrust vector of angular fluxes for each track */
  thrust::device_vector<NEW_PRECISION> _boundary_flux;

  /** Thrust vector of leakages for each track */
  thrust::device_vector<NEW_PRECISION> _boundary_leakage;

  /** Thrust vector of FSR scalar fluxes */
  thrust::device_vector<NEW_PRECISION> _scalar_flux;

  /** Thrust vector of old FSR scalar fluxes */
  thrust::device_vector<NEW_PRECISION> _old_scalar_flux;

  /** Thrust vector of fixed sources in each FSR */
  thrust::device_vector<NEW_PRECISION> _fixed_sources;

  /** Thrust vector of source / sigma_t in each FSR */
  thrust::device_vector<NEW_PRECISION> _reduced_sources;

  /** Map of Material IDs to indices in _materials array */
  std::map<int, int> _material_IDs_to_indices;

  int computeScalarTrackIndex(int i, int j);

  void initializePolarQuadrature();
  void initializeExpEvaluator();
  void initializeFSRs();
  void initializeMaterials();
  void initializeTracks();
  void initializeFluxArrays();
  void initializeSourceArrays();

  void zeroTrackFluxes();
  void flattenFSRFluxes(NEW_PRECISION value);
  void storeFSRFluxes();
  void normalizeFluxes();
  void computeFSRSources();
  void transportSweep();
  void addSourceToScalarFlux();
  void computeKeff();
  double computeResidual(residualType res_type);

public:

  GPUSolver(TrackGenerator* track_generator=NULL);
  virtual ~GPUSolver();

  int getNumThreadBlocks();

  /**
   * @brief Returns the number of threads per block to execute on the GPU.
   * @return the number of threads per block
   */
  int getNumThreadsPerBlock();
  NEW_PRECISION getFSRScalarFlux(int fsr_id, int group);
  NEW_PRECISION getFSRSource(int fsr_id, int group);

  void setNumThreadBlocks(int num_blocks);
  void setNumThreadsPerBlock(int num_threads);
  void setFixedSourceByFSR(int fsr_id, int group, 
                           NEW_PRECISION source);
  void setGeometry(Geometry* geometry);
  void setTrackGenerator(TrackGenerator* track_generator);

  void computeFSRFissionRates(double* fission_rates, int num_FSRs);
};


#endif /* GPUSOLVER_H_ */
