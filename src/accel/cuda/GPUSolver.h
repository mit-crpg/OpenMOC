/**
 * @file GPUSolver.h
 * @brief The GPUSolver class and CUDA physics kernels.
 * @date August 5, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef GPUSOLVER_H_
#define GPUSOLVER_H_

#ifdef __cplusplus
#include "../../constants.h"
#include "../../Solver.h"
#endif

#include <thrust/reduce.h>
#include <thrust/device_vector.h>
#include <sm_13_double_functions.h>
#include <sm_20_atomic_functions.h>
#include "clone.h"
#include "GPUExpEvaluator.h"


/** Indexing macro for the scalar flux in each FSR and energy group */
#define scalar_flux(tid,e) (scalar_flux[(tid)*(*num_groups) + (e)])

/** Indexing macro for the total source divided by the total cross-section,
 *  \f$ \frac{Q}{\Sigma_t} \f$, in each FSR and energy group */
#define reduced_sources(tid,e) (reduced_sources[(tid)*(*num_groups) + (e)])

/** Indexing macro for the azimuthal and polar weights */
#define polar_weights(i,p) (polar_weights[(i)*(*num_polar) + (p)])

/** Indexing macro for the angular fluxes for each polar angle and energy
 *  group for a given Track */
#define boundary_flux(t,pe2) (boundary_flux[2*(t)*(*polar_times_groups)+(pe2)])

/** The maximum number of polar angles to reserve constant memory on GPU */
#define MAX_POLAR_ANGLES 10

/** The maximum number of azimuthal angles to reserve constant memory on GPU */
#define MAX_AZIM_ANGLES 256


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

  /** A pointer to the array of Tracks on the device */
  dev_track* _dev_tracks;

  /** An array of the cumulative number of Tracks for each azimuthal angle */
  int* _track_index_offsets;

  /** A pointer to the Thrust vector of fission rates in each FSR */
  FP_PRECISION* _fission_sources;

  /** A pointer to the Thrust vector of scatter rates in each FSR */
  FP_PRECISION* _scatter_sources;

  /** Thrust vector of fission sources in each FSR */
  thrust::device_vector<FP_PRECISION> _fission_sources_vec;

  /** Thrust vector of total reaction rates in each FSR */
  thrust::device_vector<FP_PRECISION> _total_vec;

  /** Thrust vector of fission rates in each FSR */
  thrust::device_vector<FP_PRECISION> _fission_vec;

  /** Thrust vector of scatter rates in each FSR */
  thrust::device_vector<FP_PRECISION> _scatter_vec;

  /** Thrust vector of source residuals in each FSR */
  thrust::device_vector<FP_PRECISION> _source_residuals_vec;

  /** Thrust vector of leakages for each track */
  thrust::device_vector<FP_PRECISION> _boundary_leakage_vec;

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
  void initializeThrustVectors();

  void zeroTrackFluxes();
  void flattenFSRFluxes(FP_PRECISION value);
  void flattenFSRSources(FP_PRECISION value);
  void normalizeFluxes();
  FP_PRECISION computeFSRSources();
  void addSourceToScalarFlux();
  void computeKeff();
  void transportSweep();

public:

  GPUSolver(TrackGenerator* track_generator=NULL);
  virtual ~GPUSolver();

  int getNumThreadBlocks();
  int getNumThreadsPerBlock();
  FP_PRECISION getFSRScalarFlux(int fsr_id, int energy_group);
  FP_PRECISION* getFSRScalarFluxes();
  FP_PRECISION getFSRSource(int fsr_id, int energy_group);

  void setNumThreadBlocks(int num_blocks);
  void setNumThreadsPerBlock(int num_threads);
  void setGeometry(Geometry* geometry);
  void setTrackGenerator(TrackGenerator* track_generator);

  void computeFSRFissionRates(double* fission_rates, int num_FSRs);
};


#endif /* GPUSOLVER_H_ */
