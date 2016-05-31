#include "GPUSolver.h"

/** The number of azimuthal angles */
__constant__ int num_azim[1];

/** The number of energy groups */
__constant__ int num_groups[1];

/** The number of FSRs */
__constant__ int num_FSRs[1];

/** The number of polar angles */
__constant__ int num_polar[1];

/** Half the number of polar angles */
__constant__ int num_polar_2[1];

/** The number of polar angles times energy groups */
__constant__ int polar_times_groups[1];

/** An array for the sines of the polar angle in the Quadrature set */
__constant__ FP_PRECISION sin_thetas[MAX_POLAR_ANGLES_GPU];

/** An array of the weights from the Quadrature set */
__constant__ FP_PRECISION weights[MAX_POLAR_ANGLES_GPU*MAX_AZIM_ANGLES_GPU];

/** The total number of Tracks */
__constant__ int tot_num_tracks[1];

/** An GPUExpEvaluator object to compute exponentials */
__constant__ GPUExpEvaluator exp_evaluator;


/**
 * @brief A struct used to check if a value on the GPU is equal to INF.
 * @details This is used as a predicate in Thrust routines.
 */
struct isinf_test {
  /**
   * @brief Checks if a double precision value is INF.
   * @param a the value to check
   * @return true if equal to INF, false otherwise
   */
  __host__ __device__ bool operator()(double a) {
    return isinf(a);
  }

  /**
   * @brief Checks if a single precision value is INF.
   * @param a the value to check
   * @return true if equal to INF, false otherwise
   */
  __host__ __device__ bool operator()(float a) {
    return isinf(a);
  }
};


/**
 * @brief A struct used to check if a value on the GPU is equal to NaN.
 * @details This is used as a predicate in Thrust routines.
 */
struct isnan_test {
  /**
   * @brief Checks if a double precision value is NaN.
   * @param a the value to check
   * @return true if equal to NaN, false otherwise
   */
  __host__ __device__ bool operator()(double a) {
    return isnan(a);
  }

  /**
   * @brief Checks if a single precision value is NaN.
   * @param a the value to check
   * @return true if equal to NaN, false otherwise
   */
  __host__ __device__ bool operator()(float a) {
    return isnan(a);
  }
};


/**
 * @brief A functor to multiply all elements in a Thrust vector by a constant.
 * @param constant the constant to multiply the vector
 */
template< typename T >
struct multiplyByConstant {

public:
  /* The constant to multiply by */
  const T constant;

  /**
   * @brief Constructor for the functor.
   * @param constant to multiply each element in a Thrust vector
   */
  multiplyByConstant(T constant) : constant(constant) {}

  /**
   * @brief Multiply an element in a Thrust vector.
   * @param VecElem the element to multiply
   */
  __host__ __device__ void operator()(T& VecElem) const {
    VecElem = VecElem * constant;
  }
};


/**
 * @class This provides a templated interface for a strided iterator over
 *        a Thrust device_vector on a GPU.
 * @details This code is taken from the Thrust examples site on 1/20/2015:
 *           https://github.com/thrust/thrust/blob/master/examples/strided_range.cu
 */
template <typename Iterator>
class strided_range {

public:

  typedef typename thrust::iterator_difference<Iterator>::type difference_type;

  struct stride_functor : public thrust::unary_function<difference_type,difference_type> {

    difference_type stride;

    stride_functor(difference_type stride) : stride(stride) { }

    __host__ __device__ difference_type operator()(const difference_type& i) const {
      return stride * i;
    }
  };

  typedef typename thrust::counting_iterator<difference_type> CountingIterator;
  typedef typename thrust::transform_iterator<stride_functor, CountingIterator>
    TransformIterator;
  typedef typename thrust::permutation_iterator<Iterator,TransformIterator>
    PermutationIterator;
  typedef PermutationIterator iterator;

  /**
   * @brief The strided iterator constructor.
   */
  strided_range(Iterator first, Iterator last, difference_type stride)
    : first(first), last(last), stride(stride) { }

  /**
   * @brief Get the first element in the iterator.
   * @return the first element in the iterator
   */
  iterator begin(void) const {
    return PermutationIterator(first,
      TransformIterator(CountingIterator(0), stride_functor(stride)));
  }

  /**
   * @brief Get the last element in the iterator.
   * @return the last element in the iterator
   */
  iterator end(void) const {
    return begin() + ((last - first) + (stride - 1)) / stride;
  }

protected:

  /** The first element in the underlying device_vector as set by the constructor */
  Iterator first;

  /** The last element in the underlying device_vector as set by the constructor */
  Iterator last;

  /** The stride to use when iterating over the underlying device_vector */
  difference_type stride;

};


/**
 * @brief Compute the total fission source from all FSRs.
 * @param FSR_volumes an array of FSR volumes
 * @param FSR_materials an array of FSR Material indices
 * @param materials an array of dev_materials on the device
 * @param scalar_flux the scalar flux in each FSR and energy group
 * @param fission_sources array of fission sources in each FSR and energy group
 */
__global__ void computeFissionSourcesOnDevice(FP_PRECISION* FSR_volumes,
                                              int* FSR_materials,
                                              dev_material* materials,
                                              FP_PRECISION* scalar_flux,
                                              FP_PRECISION* fission_sources) {

  /* Use a shared memory buffer for each thread's fission source */
  extern __shared__ FP_PRECISION shared_fission_source[];

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  dev_material* curr_material;
  FP_PRECISION* nu_sigma_f;
  FP_PRECISION volume, source;

  /* Initialize fission source to zero */
  shared_fission_source[threadIdx.x] = 0;

  /* Iterate over all FSRs */
  while (tid < *num_FSRs) {

    curr_material = &materials[FSR_materials[tid]];
    nu_sigma_f = curr_material->_nu_sigma_f;
    volume = FSR_volumes[tid];

    /* Iterate over energy groups and update fission source for
     * this thread block */
    for (int e=0; e < *num_groups; e++) {
      source = nu_sigma_f[e] * scalar_flux(tid,e) * volume;
      shared_fission_source[threadIdx.x] += source;
    }

    /* Increment thread id */
    tid += blockDim.x * gridDim.x;
  }

  /* Copy this thread's fission source to global memory */
  tid = threadIdx.x + blockIdx.x * blockDim.x;
  fission_sources[tid] = shared_fission_source[threadIdx.x];
}


/**
 * @brief Computes the total source (fission, scattering, fixed) in each FSR.
 * @details This method computes the total source in each region based on
 *          this iteration's current approximation to the scalar flux.
 * @param FSR_materials an array of FSR Material indices
 * @param materials an array of dev_material pointers
 * @param scalar_flux an array of FSR scalar fluxes
 * @param fixed_sources an array of fixed (user-defined) sources
 * @param reduced_sources an array of FSR sources / total xs
 * @param inverse_k_eff the inverse of keff
 */
__global__ void computeFSRSourcesOnDevice(int* FSR_materials,
                                          dev_material* materials,
                                          FP_PRECISION* scalar_flux,
                                          FP_PRECISION* fixed_sources,
                                          FP_PRECISION* reduced_sources,
                                          FP_PRECISION inverse_k_eff) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  FP_PRECISION fission_source;
  FP_PRECISION scatter_source;

  dev_material* curr_material;
  FP_PRECISION* sigma_t;
  FP_PRECISION* sigma_s;
  FP_PRECISION* fiss_mat;

  /* Iterate over all FSRs */
  while (tid < *num_FSRs) {

    curr_material = &materials[FSR_materials[tid]];

    sigma_t = curr_material->_sigma_t;
    sigma_s = curr_material->_sigma_s;
    fiss_mat = curr_material->_fiss_matrix;

    /* Compute scatter + fission source for group g */
    for (int g=0; g < *num_groups; g++) {
      scatter_source = 0;
      fission_source = 0;

      for (int g_prime=0; g_prime < *num_groups; g_prime++) {
        scatter_source += sigma_s[g*(*num_groups)+g_prime] * scalar_flux(tid,g_prime);
        fission_source += fiss_mat[g*(*num_groups)+g_prime] * scalar_flux(tid,g_prime);
      }

      fission_source *= inverse_k_eff;

      /* Compute total (scatter+fission+fixed) reduced source */
      reduced_sources(tid,g) = fixed_sources(tid,g);
      reduced_sources(tid,g) += scatter_source + fission_source;
      reduced_sources(tid,g) *= ONE_OVER_FOUR_PI;
      reduced_sources(tid,g) = __fdividef(reduced_sources(tid,g), sigma_t[g]);
    }

    /* Increment the thread id */
    tid += blockDim.x * gridDim.x;
  }
}


/**
 * @brief Computes the total fission source in each FSR in each energy group
 * @details This method is a helper routine for the openmoc.krylov submodule.
 *          This routine computes the total fission source in each FSR. If the
 *          divide_sigma_t parameter is true then the fission source will be
 *          divided by the total cross-section in each FSR.
 * @param FSR_materials an array of FSR Material indices
 * @param materials an array of dev_material pointers
 * @param divide_sigma_t a boolean indicating whether to divide by the total xs
 * @param scalar_flux an array of FSR scalar fluxes
 * @param reduced_sources an array of FSR fission sources
 */
__global__ void computeFSRFissionSourcesOnDevice(int* FSR_materials,
                                                 dev_material* materials,
						 bool divide_sigma_t,
                                                 FP_PRECISION* scalar_flux,
                                                 FP_PRECISION* reduced_sources) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  FP_PRECISION fission_source;

  dev_material* curr_material;
  FP_PRECISION* sigma_t;
  FP_PRECISION* fiss_mat;

  /* Iterate over all FSRs */
  while (tid < *num_FSRs) {

    curr_material = &materials[FSR_materials[tid]];

    sigma_t = curr_material->_sigma_t;
    fiss_mat = curr_material->_fiss_matrix;

    /* Compute fission source for group g */
    for (int g=0; g < *num_groups; g++) {
      fission_source = 0;

      for (int g_prime=0; g_prime < *num_groups; g_prime++)
        fission_source += fiss_mat[g*(*num_groups)+g_prime] * scalar_flux(tid,g_prime);

      /* Set the reduced fission source for FSR tid in group g */
      reduced_sources(tid,g) = fission_source;
      reduced_sources(tid,g) *= ONE_OVER_FOUR_PI;
      if (divide_sigma_t)
        reduced_sources(tid,g) = __fdividef(reduced_sources(tid,g), sigma_t[g]);
    }

    /* Increment the thread id */
    tid += blockDim.x * gridDim.x;
  }
}


/**
 * @brief Computes the total scattering source in each FSR and energy group.
 * @details This method is a helper routine for the openmoc.krylov submodule.
 *          This routine computes the total scatter source in each FSR. If the
 *          divide_sigma_t parameter is true then the scatter source will be
 *          divided by the total cross-section in each FSR.
 * @param FSR_materials an array of FSR Material indices
 * @param materials an array of dev_material pointers
 * @param divide_sigma_t a boolean indicating whether to divide by the total xs
 * @param scalar_flux an array of FSR scalar fluxes
 * @param reduced_sources an array of FSR scatter sources
 */
__global__ void computeFSRScatterSourcesOnDevice(int* FSR_materials,
                                                 dev_material* materials,
						 bool divide_sigma_t,
                                                 FP_PRECISION* scalar_flux,
                                                 FP_PRECISION* reduced_sources) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  FP_PRECISION scatter_source;

  dev_material* curr_material;
  FP_PRECISION* sigma_s;
  FP_PRECISION* sigma_t;

  /* Iterate over all FSRs */
  while (tid < *num_FSRs) {

    curr_material = &materials[FSR_materials[tid]];

    sigma_s = curr_material->_sigma_s;
    sigma_t = curr_material->_sigma_t;

    /* Compute total scattering source for this FSR in group g */
    for (int g=0; g < *num_groups; g++) {
      scatter_source = 0;

      for (int g_prime=0; g_prime < *num_groups; g_prime++)
        scatter_source += sigma_s[g*(*num_groups)+g_prime] * scalar_flux(tid,g_prime);

      /* Set the reduced scatter source for FSR tid in group g */
      reduced_sources(tid,g) = scatter_source;
      reduced_sources(tid,g) *= ONE_OVER_FOUR_PI;
      if (divide_sigma_t)
        reduced_sources(tid,g) = __fdividef(reduced_sources(tid,g), sigma_t[g]);
    }

    /* Increment the thread id */
    tid += blockDim.x * gridDim.x;
  }
}


/**
 * @brief Perform an atomic addition in double precision to an array address.
 * @details This method is straight out of CUDA C Developers Guide (cc 2013).
 * @param address the array memory address
 * @param val the value to add to the array
 * @return the atomically added array value and input value
 */
__device__ double atomicAdd(double* address, double val) {

  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;

  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val +
                    __longlong_as_double(assumed)));
  } while (assumed != old);

  return __longlong_as_double(old);
}


/**
 * @brief Computes the contribution to the FSR scalar flux from a Track
 *        segment in a single energy group.
 * @details This method integrates the angular flux for a Track segment across
 *        energy groups and polar angles, and tallies it into the FSR scalar
 *        flux, and updates the Track's angular flux.
 * @param curr_segment a pointer to the Track segment of interest
 * @param azim_index a pointer to the azimuthal angle index for this segment
 * @param energy_group the energy group of interest
 * @param materials the array of dev_material pointers
 * @param track_flux a pointer to the Track's angular flux
 * @param reduced_sources the array of FSR sources / total xs
 * @param scalar_flux the array of FSR scalar fluxes
 */
__device__ void tallyScalarFlux(dev_segment* curr_segment,
                                int azim_index,
                                int energy_group,
                                dev_material* materials,
                                FP_PRECISION* track_flux,
                                FP_PRECISION* reduced_sources,
                                FP_PRECISION* scalar_flux) {

  int fsr_id = curr_segment->_region_uid;
  FP_PRECISION length = curr_segment->_length;
  dev_material* curr_material = &materials[curr_segment->_material_index];
  FP_PRECISION* sigma_t = curr_material->_sigma_t;

  /* The change in angular flux long this Track segment in this FSR */
  FP_PRECISION delta_psi;
  FP_PRECISION exponential;

  /* Zero the FSR scalar flux contribution from this segment and energy group */
  FP_PRECISION fsr_flux = 0.0;

  /* Loop over polar angles */
  for (int p=0; p < *num_polar_2; p++) {
    exponential =
      exp_evaluator.computeExponential(sigma_t[energy_group] * length, p);
    delta_psi = (track_flux[p] - reduced_sources(fsr_id,energy_group));
    delta_psi *= exponential;
    fsr_flux += delta_psi * weights(azim_index,p);
    track_flux[p] -= delta_psi;
  }

  /* Atomically increment the scalar flux for this FSR */
  atomicAdd(&scalar_flux(fsr_id,energy_group), fsr_flux);
}


/**
 * @brief Updates the boundary flux for a Track given boundary conditions.
 * @details For reflective and periodic boundary conditions, the outgoing
 *          boundary flux for the Track is given to the corresponding reflecting
 *          or periodic Track. For vacuum boundary conditions, the outgoing flux
 *          is tallied as leakage. Note: Only one energy group is transferred
 *          by this routine.
 * @param curr_track a pointer to the Track of interest
 * @param azim_index a pointer to the azimuthal angle index for this segment
 * @param track_flux an array of the outgoing Track flux
 * @param boundary_flux an array of all angular fluxes
 * @param weights an array of Quadrature weights
 * @param energy_angle_index the energy group index
 * @param direction the Track direction (forward - true, reverse - false)
 */
__device__ void transferBoundaryFlux(dev_track* curr_track,
                                     int azim_index,
                                     FP_PRECISION* track_flux,
                                     FP_PRECISION* boundary_flux,
                                     int energy_angle_index,
                                     bool direction) {

  int start = energy_angle_index;
  bool transfer_flux;
  int track_out_id;

  /* For the "forward" direction */
  if (direction) {
    transfer_flux = curr_track->_transfer_flux_out;
    track_out_id = curr_track->_track_out;
    start += curr_track->_next_out * (*polar_times_groups);
  }

  /* For the "reverse" direction */
  else {
    transfer_flux = curr_track->_transfer_flux_in;
    track_out_id = curr_track->_track_in;
    start += curr_track->_next_in * (*polar_times_groups);
  }

  FP_PRECISION* track_out_flux = &boundary_flux(track_out_id,start);

  /* Put Track's flux in the shared memory temporary flux array */
  for (int p=0; p < *num_polar_2; p++)
    track_out_flux[p] = track_flux[p] * transfer_flux;
}


/**
 * @brief This method performs one transport sweep of one halfspace of all
 *        azimuthal angles, tracks, segments, polar angles and energy groups.
 * @details The method integrates the flux along each track and updates the
 *          boundary fluxes for the corresponding output Track, while updating
 *          the scalar flux in each FSR.
 * @param scalar_flux an array of FSR scalar fluxes
 * @param boundary_flux an array of Track boundary fluxes
 * @param reduced_sources an array of FSR sources / total xs
 * @param materials an array of dev_material pointers
 * @param tracks an array of Tracks
 * @param tid_offset the Track offset for azimuthal angle halfspace
 * @param tid_max the upper bound on the Track IDs for this azimuthal
 *                angle halfspace
 */
__global__ void transportSweepOnDevice(FP_PRECISION* scalar_flux,
                                       FP_PRECISION* boundary_flux,
                                       FP_PRECISION* reduced_sources,
                                       dev_material* materials,
                                       dev_track* tracks,
                                       int tid_offset,
                                       int tid_max) {

  /* Shared memory buffer for each thread's angular flux */
  extern __shared__ FP_PRECISION temp_flux[];
  FP_PRECISION* track_flux;

  int tid = tid_offset + threadIdx.x + blockIdx.x * blockDim.x;
  int track_id = tid / *num_groups;
  int track_flux_index = threadIdx.x * (*num_polar);
  int energy_group = tid % (*num_groups);
  int energy_angle_index = energy_group * (*num_polar_2);

  dev_track* curr_track;
  int azim_index;
  int num_segments;
  dev_segment* curr_segment;

  /* Iterate over Track with azimuthal angles in (0, pi/2) */
  while (track_id < tid_max) {

    /* Initialize local registers with important data */
    curr_track = &tracks[track_id];
    azim_index = curr_track->_azim_angle_index;
    num_segments = curr_track->_num_segments;

    /* Retrieve pointer to thread's shared memory buffer for angular flux */
    track_flux = &temp_flux[track_flux_index];

    /* Put Track's flux in the shared memory temporary flux array */
    for (int p=0; p < *num_polar_2; p++) {

      /* Forward flux along this Track */
      track_flux[p] = boundary_flux(track_id,p+energy_angle_index);

      /* Reverse flux along this Track */
      track_flux[(*num_polar_2) + p] =
            boundary_flux(track_id,p+energy_angle_index+(*polar_times_groups));
    }

    /* Loop over each Track segment in forward direction */
    for (int i=0; i < num_segments; i++) {
      curr_segment = &curr_track->_segments[i];
      tallyScalarFlux(curr_segment, azim_index, energy_group, materials,
                      track_flux, reduced_sources, scalar_flux);
    }

    /* Transfer boundary angular flux to outgoing Track */
    transferBoundaryFlux(curr_track, azim_index, track_flux, boundary_flux,
                         energy_angle_index, true);

    /* Loop over each Track segment in reverse direction */
    track_flux = &temp_flux[track_flux_index + (*num_polar_2)];

    for (int i=num_segments-1; i > -1; i--) {
      curr_segment = &curr_track->_segments[i];
      tallyScalarFlux(curr_segment, azim_index, energy_group, materials,
                      track_flux, reduced_sources, scalar_flux);
    }

    /* Transfer boundary angular flux to outgoing Track */
    transferBoundaryFlux(curr_track, azim_index, track_flux, boundary_flux,
                         energy_angle_index, false);

    /* Update the indices for this thread to the next Track, energy group */
    tid += blockDim.x * gridDim.x;
    track_id = tid / *num_groups;
    energy_group = tid % (*num_groups);
    energy_angle_index = energy_group * (*num_polar_2);
  }
}


/**
 * @brief Add the source term contribution in the transport equation to
 *        the FSR scalar flux on the GPU.
 * @param scalar_flux an array of FSR scalar fluxes
 * @param reduced_sources an array of FSR sources / total xs
 * @param FSR_volumes an array of FSR volumes
 * @param FSR_materials an array of FSR material indices
 * @param materials an array of dev_material pointers
 */
__global__ void addSourceToScalarFluxOnDevice(FP_PRECISION* scalar_flux,
                                              FP_PRECISION* reduced_sources,
                                              FP_PRECISION* FSR_volumes,
                                              int* FSR_materials,
                                              dev_material* materials) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  FP_PRECISION volume;

  dev_material* curr_material;
  FP_PRECISION* sigma_t;

  /* Iterate over all FSRs */
  while (tid < *num_FSRs) {

    curr_material = &materials[FSR_materials[tid]];
    volume = FSR_volumes[tid];
    sigma_t = curr_material->_sigma_t;

    /* Iterate over all energy groups */
    for (int i=0; i < *num_groups; i++) {
      scalar_flux(tid,i) = __fdividef(scalar_flux(tid,i),
                                     (sigma_t[i] * volume));
      scalar_flux(tid,i) += FOUR_PI * reduced_sources(tid,i);
    }

    /* Increment thread id */
    tid += blockDim.x * gridDim.x;
  }
}


/**
 * @brief Compute the total volume-intergrated fission source from
 *        all FSRs and energy groups.
 * @param FSR_volumes an array of the FSR volumes
 * @param FSR_materials an array of the FSR Material indices
 * @param materials an array of the dev_material pointers
 * @param scalar_flux an array of FSR scalar fluxes
 * @param fission an array of FSR nu-fission rates
 */
__global__ void computeFSRFissionRatesOnDevice(FP_PRECISION* FSR_volumes,
                                               int* FSR_materials,
                                               dev_material* materials,
                                               FP_PRECISION* scalar_flux,
                                               FP_PRECISION* fission) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  dev_material* curr_material;
  FP_PRECISION* nu_sigma_f;
  FP_PRECISION volume;

  FP_PRECISION fiss = 0.;

  /* Iterate over all FSRs */
  while (tid < *num_FSRs) {

    curr_material = &materials[FSR_materials[tid]];
    nu_sigma_f = curr_material->_nu_sigma_f;
    volume = FSR_volumes[tid];

    FP_PRECISION curr_fiss = 0.;

    /* Compute fission rates rates for this thread block */
    for (int e=0; e < *num_groups; e++)
      curr_fiss += nu_sigma_f[e] * scalar_flux(tid,e);

    fiss += curr_fiss * volume;

    /* Increment thread id */
    tid += blockDim.x * gridDim.x;
  }

  /* Copy this thread's fission to global memory */
  tid = threadIdx.x + blockIdx.x * blockDim.x;
  fission[tid] = fiss;
}


/**
 * @brief Constructor initializes arrays for dev_tracks and dev_materials..
 * @details The constructor initalizes the number of CUDA threads and thread
 *          blocks each to a default of 64.
 * @param track_generator an optional pointer to the TrackjGenerator
 */
GPUSolver::GPUSolver(TrackGenerator* track_generator) :

  Solver(track_generator) {

  /* The default number of thread blocks and threads per thread block */
  _B = 64;
  _T = 64;

  _materials = NULL;
  _dev_tracks = NULL;
  _FSR_materials = NULL;

  if (track_generator != NULL)
    setTrackGenerator(track_generator);
}


/**
 * @brief Solver destructor frees all memory on the device, including arrays
 *        for the FSR scalar fluxes and sources and Track boundary fluxes.
 */
GPUSolver::~GPUSolver() {

  if (_FSR_volumes != NULL) {
    cudaFree(_FSR_volumes);
    _FSR_volumes = NULL;
  }

  if (_FSR_materials != NULL) {
    cudaFree(_FSR_materials);
    _FSR_materials = NULL;
  }

  if (_materials != NULL) {
    cudaFree(_materials);
    _materials = NULL;
  }

  if (_dev_tracks != NULL) {
    cudaFree(_dev_tracks);
    _dev_tracks = NULL;
  }

  /* Clear Thrust vectors's memory on the device */
  _boundary_flux.clear();
  _scalar_flux.clear();
  _old_scalar_flux.clear();
  _fixed_sources.clear();
  _reduced_sources.clear();
}


/**
 * @brief Returns the number of thread blocks to execute on the GPU.
 * @return the number of thread blocks
 */
int GPUSolver::getNumThreadBlocks() {
  return _B;
}


/**
 * @brief Returns the number of threads per block to execute on the GPU.
 * @return the number of threads per block
 */
int GPUSolver::getNumThreadsPerBlock() {
  return _T;
}


/**
 * @brief Returns the source for some energy group for a flat source region
 * @details This is a helper routine used by the openmoc.process module.
 * @param fsr_id the ID for the FSR of interest
 * @param group the energy group of interest
 * @return the flat source region source
 */
FP_PRECISION GPUSolver::getFSRSource(int fsr_id, int group) {

  if (fsr_id >= _num_FSRs)
    log_printf(ERROR, "Unable to return a source for FSR ID = %d "
               "since the max FSR ID = %d", fsr_id, _num_FSRs-1);

  else if (fsr_id < 0)
    log_printf(ERROR, "Unable to return a source for FSR ID = %d "
               "since FSRs do not have negative IDs", fsr_id);

  else if (group-1 >= _num_groups)
    log_printf(ERROR, "Unable to return a source in group %d "
               "since there are only %d groups", group, _num_groups);

  else if (group <= 0)
    log_printf(ERROR, "Unable to return a source in group %d "
               "since groups must be greater or equal to 1", group);

  else if (_scalar_flux.size() == 0)
    log_printf(ERROR, "Unable to return a source "
               "since it has not yet been computed");

  /* Get host material */
  Material* host_material = _geometry->findFSRMaterial(fsr_id);

  /* Get cross sections and scalar flux */
  FP_PRECISION* sigma_s = host_material->getSigmaS();
  FP_PRECISION* fiss_mat = host_material->getFissionMatrix();

  FP_PRECISION* fsr_scalar_fluxes = new FP_PRECISION[_num_groups];
  FP_PRECISION* scalar_flux =
       thrust::raw_pointer_cast(&_scalar_flux[0]);
  cudaMemcpy((void*)fsr_scalar_fluxes, (void*)&scalar_flux[fsr_id*_num_groups],
             _num_groups * sizeof(FP_PRECISION),
             cudaMemcpyDeviceToHost);

  FP_PRECISION fission_source = 0.0;
  FP_PRECISION scatter_source = 0.0;
  FP_PRECISION total_source;

  /* Compute total scattering and fission sources for this FSR */
  for (int g=0; g < _num_groups; g++) {
    scatter_source += sigma_s[(group-1)*(_num_groups)+g]
                      * fsr_scalar_fluxes[g];
    fission_source += fiss_mat[(group-1)*(_num_groups)+g]
                      * fsr_scalar_fluxes[g];
  }

  fission_source /= _k_eff;

  /* Compute the total source */
  total_source = fission_source + scatter_source;

  /* Add in fixed source (if specified by user) */
  total_source += _fixed_sources(fsr_id,group-1);

  /* Normalize to solid angle for isotropic approximation */
  total_source *= ONE_OVER_FOUR_PI;

  delete [] fsr_scalar_fluxes;

  return total_source;
}


/**
 * @brief Returns the scalar flux for some FSR and energy group.
 * @param fsr_id the ID for the FSR of interest
 * @param group the energy group of interest
 * @return the FSR scalar flux
 */
FP_PRECISION GPUSolver::getFlux(int fsr_id, int group) {

  if (fsr_id >= _num_FSRs)
    log_printf(ERROR, "Unable to return a scalar flux for FSR ID = %d "
               "since the max FSR ID = %d", fsr_id, _num_FSRs-1);

  else if (fsr_id < 0)
    log_printf(ERROR, "Unable to return a scalar flux for FSR ID = %d "
               "since FSRs do not have negative IDs", fsr_id);

  else if (group-1 >= _num_groups)
    log_printf(ERROR, "Unable to return a scalar flux in group %d "
               "since there are only %d groups", group, _num_groups);

  else if (group <= 0)
    log_printf(ERROR, "Unable to return a scalar flux in group %d "
               "since groups must be greater or equal to 1", group);

  if (_scalar_flux.size() == 0)
    log_printf(ERROR, "Unable to return a scalar flux "
               "since it has not yet been computed");

  return _scalar_flux(fsr_id,group-1);
}


/**
 * @brief Fills an array with the scalar fluxes on the GPU.
 * @details This class method is a helper routine called by the OpenMOC
 *          Python "openmoc.krylov" module for Krylov subspace methods.
 *          Although this method appears to require two arguments, in
 *          reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_fluxes = num_groups * num_FSRs
 *          fluxes = solver.getFluxes(num_fluxes)
 * @endcode
 *
 * @param fluxes an array of FSR scalar fluxes in each energy group
 * @param num_fluxes the total number of FSR flux values
 */
void GPUSolver::getFluxes(FP_PRECISION* out_fluxes, int num_fluxes) {

  if (num_fluxes != _num_groups * _num_FSRs)
    log_printf(ERROR, "Unable to get FSR scalar fluxes since there are "
               "%d groups and %d FSRs which does not match the requested "
               "%d flux values", _num_groups, _num_FSRs, num_fluxes);

  else if (_scalar_flux.size() == 0)
    log_printf(ERROR, "Unable to get FSR scalar fluxes since they "
               "have not yet been allocated on the device");

  FP_PRECISION* scalar_flux =
       thrust::raw_pointer_cast(&_scalar_flux[0]);

  /* Copy the fluxes from the GPU to the input array */
  cudaMemcpy((void*)out_fluxes, (void*)scalar_flux,
            num_fluxes * sizeof(FP_PRECISION), cudaMemcpyDeviceToHost);
}


/**
 * @brief Sets the number of thread blocks (>0) for CUDA kernels.
 * @param num_blocks the number of thread blocks
 */
void GPUSolver::setNumThreadBlocks(int num_blocks) {

  if (num_blocks < 0)
    log_printf(ERROR, "Unable to set the number of CUDA thread blocks "
               "to %d since it is a negative number", num_blocks);

  _B = num_blocks;
}


/**
 * @brief Sets the number of threads per block (>0) for CUDA kernels.
 * @param num_threads the number of threads per block
 */
void GPUSolver::setNumThreadsPerBlock(int num_threads) {

  if (num_threads < 0)
    log_printf(ERROR, "Unable to set the number of CUDA threads per block "
               "to %d since it is a negative number", num_threads);

  _T = num_threads;
}


/**
 * @brief Sets the Geometry for the Solver.
 * @details This is a private setter method for the Solver and is not
 *          intended to be called by the user.
 * @param geometry a pointer to a Geometry object
 */
void GPUSolver::setGeometry(Geometry* geometry) {

  Solver::setGeometry(geometry);

  std::map<int, Material*> host_materials=_geometry->getAllMaterials();
  std::map<int, Material*>::iterator iter;
  int material_index = 0;

  /* Iterate through all Materials and clone them as dev_material structs
   * on the device */
  for (iter=host_materials.begin(); iter != host_materials.end(); ++iter) {
    _material_IDs_to_indices[iter->second->getId()] = material_index;
    material_index++;
  }
}


/**
 * @brief Sets the Solver's TrackGenerator with characteristic Tracks.
 * @details The TrackGenerator must already have generated Tracks and have
 *          used ray tracing to segmentize them across the Geometry. This
 *          should be initated in Python prior to assigning the TrackGenerator
 *          to the Solver:
 *
 * @code
 *          track_generator.generateTracks()
 *          solver.setTrackGenerator(track_generator)
 * @endcode
 *
 * @param track_generator a pointer to a TrackGenerator object
 */
void GPUSolver::setTrackGenerator(TrackGenerator* track_generator) {
  Solver::setTrackGenerator(track_generator);
  initializeTracks();
  copyQuadrature();
}


/**
 * @brief Set the flux array for use in transport sweep source calculations.
 * @detail This is a helper method for the checkpoint restart capabilities,
 *         as well as the IRAMSolver in the openmoc.krylov submodule. This
 *         routine may be used as follows from within Python:
 *
 * @code
 *          num_FSRs = solver.getGeometry.getNumFSRs()
 *          num_groups = solver.getGeometry.getNumEnergyGroups()
 *          fluxes = numpy.random.rand(num_FSRs * num_groups, dtype=np.float)
 *          solver.setFluxes(fluxes)
 * @endcode
 *
 *          NOTE: This routine stores a pointer to the fluxes for the Solver
 *          to use during transport sweeps and other calculations. Hence, the
 *          flux array pointer is shared between NumPy and the Solver.
 *
 * @param in_fluxes an array with the fluxes to use
 * @param num_fluxes the number of flux values (# groups x # FSRs)
 */
void GPUSolver::setFluxes(FP_PRECISION* in_fluxes, int num_fluxes) {
  if (num_fluxes != _num_groups * _num_FSRs)
    log_printf(ERROR, "Unable to set an array with %d flux values for %d "
               " groups and %d FSRs", num_fluxes, _num_groups, _num_FSRs);

  /* Allocate array if flux arrays have not yet been initialized */
  if (_scalar_flux.size() == 0)
    initializeFluxArrays();

  FP_PRECISION* scalar_flux =
       thrust::raw_pointer_cast(&_scalar_flux[0]);

  /* Copy the input fluxes onto the GPU */
  cudaMemcpy((void*)scalar_flux, (void*)in_fluxes,
             num_fluxes * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);
  _user_fluxes = true;
}


/**
 * @brief Creates a polar quadrature object for the GPUSolver on the GPU.
 */
void GPUSolver::copyQuadrature() {

  log_printf(INFO, "Copying quadrature on the GPU...");

  if (_num_polar_2 > MAX_POLAR_ANGLES_GPU)
    log_printf(ERROR, "Unable to initialize a polar quadrature with %d "
               "angles for the GPUSolver which is limited to %d polar "
               "angles. Update the MAX_POLAR_ANGLES_GPU macro in constants.h "
               "and recompile.", _num_polar_2, MAX_POLAR_ANGLES_GPU);

  /* Copy half the number of polar angles to constant memory on the GPU */
  cudaMemcpyToSymbol(num_polar_2, (void*)&_num_polar_2, sizeof(int), 0,
                     cudaMemcpyHostToDevice);

  /* Copy the number of polar angles to constant memory on the GPU */
  int n_polar = 2 * _num_polar_2;
  cudaMemcpyToSymbol(num_polar, (void*)&n_polar,
                     sizeof(int), 0, cudaMemcpyHostToDevice);

  /* Copy the weights to constant memory on the GPU */
  int num_azim_2 = _quadrature->getNumAzimAngles() / 2;
  FP_PRECISION total_weights[num_azim_2 * _num_polar_2];
  for (int a=0; a < num_azim_2; a++)
    for (int p=0; p < _num_polar_2; p++)
      total_weights[a*_num_polar_2 + p] = _quadrature->getWeight(a, p);
  cudaMemcpyToSymbol(weights, (void*)total_weights,
      _num_polar_2 * num_azim_2 * sizeof(FP_PRECISION), 0, cudaMemcpyHostToDevice);

  /* Copy the sines of the polar angles which is needed if the user
   * requested the use of the exp intrinsic to evaluate exponentials */
  cudaMemcpyToSymbol(sin_thetas, (void*)_quadrature->getSinThetas(),
                     _num_polar_2 * sizeof(FP_PRECISION), 0,
                     cudaMemcpyHostToDevice);
}


/**
 * @brief Initializes new GPUExpEvaluator object to compute exponentials.
 */
void GPUSolver::initializeExpEvaluator() {

  Solver::initializeExpEvaluator();

  log_printf(INFO, "Initializing the exponential evaluator on the GPU...");

  /* Allocate memory for a GPUExpEvaluator on the device */
  GPUExpEvaluator* dev_exp_evaluator;
  cudaMalloc((void**)&dev_exp_evaluator, sizeof(GPUExpEvaluator));

  /* Clone ExpEvaluator from the host into GPUExpEvaluator on the device */
  clone_exp_evaluator(_exp_evaluator, dev_exp_evaluator);

  /* Copy the GPUExpEvaluator into constant memory on the GPU */
  cudaMemcpyToSymbol(exp_evaluator, (void*)dev_exp_evaluator,
                     sizeof(GPUExpEvaluator), 0, cudaMemcpyDeviceToDevice);
}


/**
 * @brief Initializes the FSR volumes and dev_materials array on the GPU.
 * @details This method assigns each FSR a unique, monotonically increasing
 *          ID, sets the Material for each FSR, and assigns a volume based on
 *          the cumulative length of all of the segments inside the FSR.
 */
void GPUSolver::initializeFSRs() {

  log_printf(NORMAL, "Initializing FSRs on the GPU...");

  /* Delete old FSRs array if it exists */
  if (_FSR_volumes != NULL) {
    cudaFree(_FSR_volumes);
    _FSR_volumes = NULL;
  }

  if (_FSR_materials != NULL) {
    cudaFree(_FSR_materials);
    _FSR_materials = NULL;
  }

  Solver::initializeFSRs();

  /* Allocate memory for all FSR volumes and dev_materials on the device */
  try{

    /* Store pointers to arrays of FSR data created on the host by the
     * the parent class Solver::initializeFSRs() routine */
    FP_PRECISION* host_FSR_volumes = _FSR_volumes;
    int* host_FSR_materials = _FSR_materials;

    /* Allocate memory on device for FSR volumes and Material indices */
    cudaMalloc((void**)&_FSR_volumes, _num_FSRs * sizeof(FP_PRECISION));
    cudaMalloc((void**)&_FSR_materials, _num_FSRs * sizeof(int));

    /* Create a temporary FSR to material indices array */
    int* FSRs_to_material_indices = new int[_num_FSRs];

    /* Populate FSR Material indices array */
    for (int i = 0; i < _num_FSRs; i++)
      FSRs_to_material_indices[i] = _material_IDs_to_indices[_geometry->
        findFSRMaterial(i)->getId()];

    /* Copy the arrays of FSR data to the device */
    cudaMemcpy((void*)_FSR_volumes, (void*)host_FSR_volumes,
      _num_FSRs * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)_FSR_materials, (void*)FSRs_to_material_indices,
      _num_FSRs * sizeof(int), cudaMemcpyHostToDevice);

    /* Copy the number of FSRs into constant memory on the GPU */
    cudaMemcpyToSymbol(num_FSRs, (void*)&_num_FSRs, sizeof(int), 0,
      cudaMemcpyHostToDevice);

    /* Free the array of FSRs data allocated by the Solver parent class */
    free(host_FSR_materials);

    /* Free the temporary array of FSRs to material indices on the host */
    free(FSRs_to_material_indices);
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for FSRs on GPU");
  }
}


/**
 * @brief Allocates all Materials data on the GPU.
 * @details This method loops over the materials in the host_materials map.
 *          Since CUDA does not support std::map data types on the device,
 *          the materials map must be converted to an array and a map created
 *          that maps a material ID to an indice in the new materials array. In
 *          initializeTracks, this map is used to convert the Material ID
 *          associated with every segment to an index in the materials array.
 * @param mode the solution type (FORWARD or ADJOINT)
 */
void GPUSolver::initializeMaterials(solverMode mode) {

  Solver::initializeMaterials(mode);

  log_printf(INFO, "Initializing materials on the GPU...");

  /* Copy the number of energy groups to constant memory on the GPU */
  cudaMemcpyToSymbol(num_groups, (void*)&_num_groups, sizeof(int), 0,
                     cudaMemcpyHostToDevice);

  /* Copy the number of polar angles times energy groups to constant memory
   * on the GPU */
  cudaMemcpyToSymbol(polar_times_groups, (void*)&_polar_times_groups,
                     sizeof(int), 0, cudaMemcpyHostToDevice);

  /* Delete old materials array if it exists */
  if (_materials != NULL)
    cudaFree(_materials);

  /* Allocate memory for all dev_materials on the device */
  try{

    std::map<int, Material*> host_materials=_geometry->getAllMaterials();
    std::map<int, Material*>::iterator iter;
    int material_index = 0;

    /* Iterate through all Materials and clone them as dev_material structs
     * on the device */
    cudaMalloc((void**)&_materials, _num_materials * sizeof(dev_material));
    for (iter=host_materials.begin(); iter != host_materials.end(); ++iter) {
      clone_material(iter->second, &_materials[material_index]);
      material_index++;
    }
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for Materials on GPU");
  }
}


/**
 * @brief Allocates memory for all Tracks on the GPU
 */
void GPUSolver::initializeTracks() {

  log_printf(INFO, "Initializing tracks on the GPU...");

  /* Delete old Tracks array if it exists */
  if (_dev_tracks != NULL)
    cudaFree(_dev_tracks);

  /* Allocate memory for all Tracks and Track offset indices on the device */
  try{

    /* Allocate array of dev_tracks */
    cudaMalloc((void**)&_dev_tracks, _tot_num_tracks * sizeof(dev_track));

    /* Iterate through all Tracks and clone them as dev_tracks on the device */
    int index;

    for (int i=0; i < _tot_num_tracks; i++) {

      clone_track(_tracks[i], &_dev_tracks[i], _material_IDs_to_indices);

      /* Get indices to next tracks along "forward" and "reverse" directions */
      index = _tracks[i]->getTrackIn()->getUid();
      cudaMemcpy((void*)&_dev_tracks[i]._track_in,
                 (void*)&index, sizeof(int), cudaMemcpyHostToDevice);

      index = _tracks[i]->getTrackOut()->getUid();
      cudaMemcpy((void*)&_dev_tracks[i]._track_out,
                 (void*)&index, sizeof(int), cudaMemcpyHostToDevice);
    }

    /* Copy the total number of Tracks into constant memory on GPU */
    cudaMemcpyToSymbol(tot_num_tracks, (void*)&_tot_num_tracks,
                       sizeof(int), 0, cudaMemcpyHostToDevice);
  }

  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for Tracks on GPU");
  }
}


/**
 * @brief Allocates memory for Track boundary angular and FSR scalar fluxes.
 * @details Deletes memory for old flux vectors if they were allocated for a
 *          previous simulation.
 */
void GPUSolver::initializeFluxArrays() {

  log_printf(INFO, "Initializing flux vectors on the GPU...");

  /* Clear Thrust vectors' memory if previously allocated */
  _boundary_flux.clear();
  _scalar_flux.clear();
  _old_scalar_flux.clear();

  /* Allocate memory for all flux arrays on the device */
  try {
    int size = 2 * _tot_num_tracks * _polar_times_groups;
    _boundary_flux.resize(size);

    size = _num_FSRs * _num_groups;
    _scalar_flux.resize(size);
    _old_scalar_flux.resize(size);
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for fluxes on GPU");
  }
}


/**
 * @brief Allocates memory for FSR source vectors on the GPU.
 * @details Deletes memory for old source vectors if they were allocated
 *          for a previous simulation.
 */
void GPUSolver::initializeSourceArrays() {

  log_printf(INFO, "Initializing source vectors on the GPU...");

  /* Clear Thrust vectors' memory if previously allocated */
  _reduced_sources.clear();
  _fixed_sources.clear();

  int size = _num_FSRs * _num_groups;

  /* Allocate memory for all source arrays on the device */
  try{
    _reduced_sources.resize(size);
    _fixed_sources.resize(size);
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for sources on GPU");
  }

  /* Initialize fixed sources to zero */
  thrust::fill(_fixed_sources.begin(), _fixed_sources.end(), 0.0);

  /* Fill fixed sources with those assigned by Cell, Material or FSR */
  initializeFixedSources();
}


/**
 * @brief Populates array of fixed sources assigned by FSR.
 */
void GPUSolver::initializeFixedSources() {

  Solver::initializeFixedSources();

  int fsr_id, group;
  std::pair<int, int> fsr_group_key;
  std::map< std::pair<int, int>, FP_PRECISION >::iterator fsr_iter;

  /* Populate fixed source array with any user-defined sources */
  for (fsr_iter = _fix_src_FSR_map.begin();
       fsr_iter != _fix_src_FSR_map.end(); ++fsr_iter) {

    /* Get the FSR with an assigned fixed source */
    fsr_group_key = fsr_iter->first;
    fsr_id = fsr_group_key.first;
    group = fsr_group_key.second;

    if (group <= 0 || group > _num_groups)
      log_printf(ERROR,"Unable to use fixed source for group %d in "
                 "a %d energy group problem", group, _num_groups);

    if (fsr_id < 0 || fsr_id >= _num_FSRs)
      log_printf(ERROR,"Unable to use fixed source for FSR %d with only "
                 "%d FSRs in the geometry", fsr_id, _num_FSRs);

    _fixed_sources(fsr_id, group-1) = _fix_src_FSR_map[fsr_group_key];
  }
}


/**
 * @brief Zero each Track's boundary fluxes for each energy group and polar
 *        angle in the "forward" and "reverse" directions.
 */
void GPUSolver::zeroTrackFluxes() {
  thrust::fill(_boundary_flux.begin(), _boundary_flux.end(), 0.0);
}


/**
 * @brief Set the scalar flux for each FSR and energy group to some value.
 * @param value the value to assign to each FSR scalar flux
 */
void GPUSolver::flattenFSRFluxes(FP_PRECISION value) {
  thrust::fill(_scalar_flux.begin(), _scalar_flux.end(), value);
}


/**
 * @brief Stores the FSR scalar fluxes in the old scalar flux array.
 */
void GPUSolver::storeFSRFluxes() {
  thrust::copy(_scalar_flux.begin(), _scalar_flux.end(),
               _old_scalar_flux.begin());
}


/**
 * @brief Normalizes all FSR scalar fluxes and Track boundary angular
 *        fluxes to the total fission source (times \f$ \nu \f$).
 */
void GPUSolver::normalizeFluxes() {

  /** Create Thrust vector of fission sources in each FSR */
  thrust::device_vector<FP_PRECISION> fission_sources_vec;
  fission_sources_vec.resize(_B * _T);
  FP_PRECISION* fission_sources =
       thrust::raw_pointer_cast(&fission_sources_vec[0]);

  FP_PRECISION* scalar_flux =
       thrust::raw_pointer_cast(&_scalar_flux[0]);

  int shared_mem = sizeof(FP_PRECISION) * _T;

  computeFissionSourcesOnDevice<<<_B, _T, shared_mem>>>(_FSR_volumes,
                                                        _FSR_materials,
                                                        _materials,
                                                        scalar_flux,
                                                        fission_sources);

  FP_PRECISION norm_factor = 1.0 / thrust::reduce(fission_sources_vec.begin(),
                                                  fission_sources_vec.end());

  /* Multiply all scalar and angular fluxes by the normalization constant */
  thrust::transform(_scalar_flux.begin(), _scalar_flux.end(),
                    thrust::constant_iterator<FP_PRECISION>(norm_factor),
                    _scalar_flux.begin(), thrust::multiplies<FP_PRECISION>());
  thrust::transform(_old_scalar_flux.begin(), _old_scalar_flux.end(),
                    thrust::constant_iterator<FP_PRECISION>(norm_factor),
                    _old_scalar_flux.begin(),
                    thrust::multiplies<FP_PRECISION>());
  thrust::transform(_boundary_flux.begin(), _boundary_flux.end(),
                    thrust::constant_iterator<FP_PRECISION>(norm_factor),
                    _boundary_flux.begin(), thrust::multiplies<FP_PRECISION>());

  /* Clear Thrust vector of FSR fission sources */
  fission_sources_vec.clear();
}


/**
 * @brief Computes the total source (fission, scattering, fixed) in each FSR.
 * @details This method computes the total source in each FSR based on
 *          this iteration's current approximation to the scalar flux.
 */
void GPUSolver::computeFSRSources() {

  FP_PRECISION* scalar_flux =
       thrust::raw_pointer_cast(&_scalar_flux[0]);
  FP_PRECISION* fixed_sources =
       thrust::raw_pointer_cast(&_fixed_sources[0]);
  FP_PRECISION* reduced_sources =
       thrust::raw_pointer_cast(&_reduced_sources[0]);

  computeFSRSourcesOnDevice<<<_B, _T>>>(_FSR_materials, _materials,
                                        scalar_flux, fixed_sources,
                                        reduced_sources, 1.0 / _k_eff);
}


/**
 * @brief Computes the fission source in each FSR.
 * @details This method computes the fission source in each FSR based on
 *          this iteration's current approximation to the scalar flux.
 */
void GPUSolver::computeFSRFissionSources() {

  FP_PRECISION* scalar_flux =
       thrust::raw_pointer_cast(&_scalar_flux[0]);
  FP_PRECISION* reduced_sources =
       thrust::raw_pointer_cast(&_reduced_sources[0]);

  computeFSRFissionSourcesOnDevice<<<_B, _T>>>(_FSR_materials, _materials, true,
                                               scalar_flux, reduced_sources);
}


/**
 * @brief Computes the scatter source in each FSR.
 * @details This method computes the scatter source in each FSR based on
 *          this iteration's current approximation to the scalar flux.
 */
void GPUSolver::computeFSRScatterSources() {

  FP_PRECISION* scalar_flux =
       thrust::raw_pointer_cast(&_scalar_flux[0]);
  FP_PRECISION* reduced_sources =
       thrust::raw_pointer_cast(&_reduced_sources[0]);

  computeFSRScatterSourcesOnDevice<<<_B, _T>>>(_FSR_materials, _materials, true,
                                               scalar_flux, reduced_sources);
}


/**
 * @brief This method performs one transport sweep of all azimuthal angles,
 *        Tracks, Track segments, polar angles and energy groups.
 * @details The method integrates the flux along each Track and updates the
 *          boundary fluxes for the corresponding output Track, while updating
 *          the scalar flux in each flat source region.
 */
void GPUSolver::transportSweep() {

  int shared_mem = _T * _num_polar_2 * 2 * sizeof(FP_PRECISION);
  int tid_offset = 0;
  int tid_max = 0;

  log_printf(DEBUG, "Transport sweep on device with %d blocks and %d threads",
             _B, _T);

  /* Get device pointer to the Thrust vectors */
  FP_PRECISION* scalar_flux =
       thrust::raw_pointer_cast(&_scalar_flux[0]);
  FP_PRECISION* boundary_flux =
       thrust::raw_pointer_cast(&_boundary_flux[0]);
  FP_PRECISION* reduced_sources =
       thrust::raw_pointer_cast(&_reduced_sources[0]);

  /* Initialize flux in each FSR to zero */
  flattenFSRFluxes(0.0);

  /* Loop over the parallel track groups and perform transport sweep on tracks
   * in that group */
  for (int g=0; g < _num_parallel_track_groups; g++) {
    tid_offset = tid_max * _num_groups;
    tid_max += _track_generator->getNumTracksByParallelGroup(g);

    transportSweepOnDevice<<<_B, _T, shared_mem>>>(scalar_flux, boundary_flux,
                                                   reduced_sources,
                                                   _materials, _dev_tracks,
                                                   tid_offset, tid_max);

    cudaDeviceSynchronize();
  }
}


/**
 * @brief Add the source term contribution in the transport equation to
 *        the FSR scalar flux.
 */
void GPUSolver::addSourceToScalarFlux() {
  FP_PRECISION* scalar_flux =
       thrust::raw_pointer_cast(&_scalar_flux[0]);
  FP_PRECISION* reduced_sources =
       thrust::raw_pointer_cast(&_reduced_sources[0]);

  addSourceToScalarFluxOnDevice<<<_B,_T>>>(scalar_flux, reduced_sources,
                                           _FSR_volumes, _FSR_materials,
                                           _materials);
}


/**
 * @brief Compute \f$ k_{eff} \f$ from successive fission sources.
 * @details This method computes the current approximation to the
 *          multiplication factor on this iteration as follows:
 *          \f$ k_{eff} = \frac{\displaystyle\sum_{i \in I}
 *                        \displaystyle\sum_{g \in G} \nu \Sigma^F_g \Phi V_{i}}
 *                        {\displaystyle\sum_{i \in I}
 *                        \displaystyle\sum_{g \in G} (\Sigma^T_g \Phi V_{i} -
 *                        \Sigma^S_g \Phi V_{i} - L_{i,g})} \f$
 */
void GPUSolver::computeKeff() {

  FP_PRECISION fission;

  thrust::device_vector<FP_PRECISION> fission_vec;
  fission_vec.resize(_B * _T);

  FP_PRECISION* fiss_ptr = thrust::raw_pointer_cast(&fission_vec[0]);
  FP_PRECISION* flux = thrust::raw_pointer_cast(&_scalar_flux[0]);

  /* Compute the total, fission and scattering reaction rates on device.
   * This kernel stores partial rates in a Thrust vector with as many
   * entries as CUDAthreads executed by the kernel */
  computeFSRFissionRatesOnDevice<<<_B, _T>>>(_FSR_volumes, _FSR_materials,
                                             _materials, flux, fiss_ptr);

  /* Compute the total fission source */
  fission = thrust::reduce(fission_vec.begin(), fission_vec.end());

  _k_eff *= fission;

  fission_vec.clear();
}


/**
 * @brief Computes the residual between source/flux iterations.
 * @param res_type the type of residuals to compute
 *        (SCALAR_FLUX, FISSION_SOURCE, TOTAL_SOURCE)
 * @return the average residual in each flat source region
 */
double GPUSolver::computeResidual(residualType res_type) {

  int norm;
  double residual;
  isinf_test inf_test;
  isnan_test nan_test;

  /* Allocate Thrust vector for residuals in each FSR */
  thrust::device_vector<double> residuals(_num_FSRs);

  if (res_type == SCALAR_FLUX) {

    norm = _num_FSRs;

    /* Allocate Thrust vector for residuals */
    thrust::device_vector<FP_PRECISION> fp_residuals(_num_FSRs * _num_groups);
    thrust::device_vector<FP_PRECISION> FSR_fp_residuals(_num_FSRs);

    /* Compute the relative flux change in each FSR and group */
    thrust::transform(_scalar_flux.begin(), _scalar_flux.end(),
                      _old_scalar_flux.begin(), fp_residuals.begin(),
                      thrust::minus<FP_PRECISION>());
    thrust::transform(fp_residuals.begin(), fp_residuals.end(),
                      _old_scalar_flux.begin(), fp_residuals.begin(),
                      thrust::divides<FP_PRECISION>());

    /* Replace INF and NaN values (from divide by zero) with 0. */
    thrust::replace_if(fp_residuals.begin(), fp_residuals.end(), inf_test, 0);
    thrust::replace_if(fp_residuals.begin(), fp_residuals.end(), nan_test, 0);

    /* Square the residuals */
    thrust::transform(fp_residuals.begin(), fp_residuals.end(),
                      fp_residuals.begin(), fp_residuals.begin(),
                      thrust::multiplies<FP_PRECISION>());

    typedef thrust::device_vector<FP_PRECISION>::iterator Iterator;

    /* Reduce flux residuals across energy groups within each FSR */
    for (int e=0; e < _num_groups; e++) {
      strided_range<Iterator> strider(fp_residuals.begin() + e,
                                      fp_residuals.end(), _num_groups);
      thrust::transform(FSR_fp_residuals.begin(), FSR_fp_residuals.end(),
                        strider.begin(), FSR_fp_residuals.begin(),
                        thrust::plus<FP_PRECISION>());
    }

    /* Copy the FP_PRECISION residual to the double precision residual */
    thrust::copy(FSR_fp_residuals.begin(),
                 FSR_fp_residuals.end(), residuals.begin());

    /* Sum up the residuals */
    residual = thrust::reduce(residuals.begin(), residuals.end());

    /* Deallocate memory for Thrust vectors */
    fp_residuals.clear();
    FSR_fp_residuals.clear();
    residuals.clear();

    /* Normalize the residual */
    residual = sqrt(residual / norm);

    return residual;
  }

  else if (res_type == FISSION_SOURCE) {

    if (_num_fissionable_FSRs == 0)
      log_printf(ERROR, "The Solver is unable to compute a "
                 "FISSION_SOURCE residual without fissionable FSRs");

    norm = _num_fissionable_FSRs;

    /* Allocate Thrust vectors for fission sources in each FSR, group */
    thrust::device_vector<FP_PRECISION> new_fission_sources_vec(_num_FSRs * _num_groups);
    thrust::device_vector<FP_PRECISION> old_fission_sources_vec(_num_FSRs * _num_groups);

    /* Allocate Thrust vectors for energy-integrated fission sources in each FSR */
    thrust::device_vector<FP_PRECISION> FSR_old_fiss_src(_num_FSRs);
    thrust::device_vector<FP_PRECISION> FSR_new_fiss_src(_num_FSRs);

    /* Cast Thrust vectors as array pointers */
    FP_PRECISION* old_fission_sources =
         thrust::raw_pointer_cast(&old_fission_sources_vec[0]);
    FP_PRECISION* new_fission_sources =
         thrust::raw_pointer_cast(&new_fission_sources_vec[0]);
    FP_PRECISION* scalar_flux =
         thrust::raw_pointer_cast(&_scalar_flux[0]);
    FP_PRECISION* old_scalar_flux =
         thrust::raw_pointer_cast(&_old_scalar_flux[0]);

    /* Compute the old and new nu-fission sources in each FSR, group */
    computeFSRFissionSourcesOnDevice<<<_B, _T>>>(_FSR_materials, _materials, false,
                                                 old_scalar_flux, old_fission_sources);
    computeFSRFissionSourcesOnDevice<<<_B, _T>>>(_FSR_materials, _materials, false,
                                                 scalar_flux, new_fission_sources);

    typedef thrust::device_vector<FP_PRECISION>::iterator Iterator;

    /* Reduce nu-fission sources across energy groups within each FSR */
    for (int e=0; e < _num_groups; e++) {
      strided_range<Iterator> old_strider(old_fission_sources_vec.begin() + e,
                                          old_fission_sources_vec.end(), _num_groups);
      strided_range<Iterator> new_strider(new_fission_sources_vec.begin() + e,
                                          new_fission_sources_vec.end(), _num_groups);
      thrust::transform(FSR_old_fiss_src.begin(), FSR_old_fiss_src.end(),
                        old_strider.begin(), FSR_old_fiss_src.begin(),
                        thrust::plus<FP_PRECISION>());
      thrust::transform(FSR_new_fiss_src.begin(), FSR_new_fiss_src.end(),
                        new_strider.begin(), FSR_new_fiss_src.begin(),
                        thrust::plus<FP_PRECISION>());
    }

    /* Compute the relative nu-fission source change in each FSR */
    thrust::transform(FSR_new_fiss_src.begin(), FSR_new_fiss_src.end(),
                      FSR_old_fiss_src.begin(), residuals.begin(),
                      thrust::minus<FP_PRECISION>());
    thrust::transform(residuals.begin(), residuals.end(),
                      FSR_old_fiss_src.begin(), residuals.begin(),
                      thrust::divides<FP_PRECISION>());

    /* Deallocate memory for Thrust vectors */
    old_fission_sources_vec.clear();
    new_fission_sources_vec.clear();
    FSR_old_fiss_src.clear();
    FSR_new_fiss_src.clear();
  }

  else if (res_type == TOTAL_SOURCE) {

    norm = _num_FSRs;

    /* Allocate Thrust vectors for fission/scatter sources in each FSR, group */
    thrust::device_vector<FP_PRECISION> new_sources_vec(_num_FSRs * _num_groups);
    thrust::device_vector<FP_PRECISION> old_sources_vec(_num_FSRs * _num_groups);
    thrust::fill(new_sources_vec.begin(), new_sources_vec.end(), 0.0);
    thrust::fill(old_sources_vec.begin(), old_sources_vec.end(), 0.0);

    /* Allocate Thrust vectors for energy-integrated fission/scatter sources in each FSR */
    thrust::device_vector<FP_PRECISION> FSR_old_src(_num_FSRs);
    thrust::device_vector<FP_PRECISION> FSR_new_src(_num_FSRs);
    thrust::fill(FSR_old_src.begin(), FSR_old_src.end(), 0.);
    thrust::fill(FSR_new_src.begin(), FSR_new_src.end(), 0.);

    /* Cast Thrust vectors as array pointers */
    FP_PRECISION* old_sources =
         thrust::raw_pointer_cast(&old_sources_vec[0]);
    FP_PRECISION* new_sources =
         thrust::raw_pointer_cast(&new_sources_vec[0]);
    FP_PRECISION* scalar_flux =
         thrust::raw_pointer_cast(&_scalar_flux[0]);
    FP_PRECISION* old_scalar_flux =
         thrust::raw_pointer_cast(&_old_scalar_flux[0]);

    /* Compute nu-fission source */

    /* Compute the old and new nu-fission sources in each FSR, group */
    computeFSRFissionSourcesOnDevice<<<_B, _T>>>(_FSR_materials, _materials, false,
                                                 old_scalar_flux, old_sources);
    computeFSRFissionSourcesOnDevice<<<_B, _T>>>(_FSR_materials, _materials, false,
                                                 scalar_flux, new_sources);

    typedef thrust::device_vector<FP_PRECISION>::iterator Iterator;

    /* Reduce nu-fission sources across energy groups within each FSR */
    for (int e=0; e < _num_groups; e++) {
      strided_range<Iterator> old_strider(old_sources_vec.begin() + e,
                                          old_sources_vec.end(), _num_groups);
      strided_range<Iterator> new_strider(new_sources_vec.begin() + e,
                                          new_sources_vec.end(), _num_groups);
      thrust::transform(FSR_old_src.begin(), FSR_old_src.end(),
                        old_strider.begin(), FSR_old_src.begin(),
                        thrust::plus<FP_PRECISION>());
      thrust::transform(FSR_new_src.begin(), FSR_new_src.end(),
                        new_strider.begin(), FSR_new_src.begin(),
                        thrust::plus<FP_PRECISION>());
    }

    /* Multiply fission sources by inverse keff */
    thrust::for_each(FSR_new_src.begin(), FSR_new_src.end(),
                     multiplyByConstant<FP_PRECISION>(1. / _k_eff));
    thrust::for_each(FSR_old_src.begin(), FSR_old_src.end(),
                     multiplyByConstant<FP_PRECISION>(1. / _k_eff));

    /* Compute scatter source */

    /* Reset sources Thrust vectors to zero */
    thrust::fill(new_sources_vec.begin(), new_sources_vec.end(), 0.0);
    thrust::fill(old_sources_vec.begin(), old_sources_vec.end(), 0.0);

    /* Compute the old and new scattering sources in each FSR, group */
    computeFSRScatterSourcesOnDevice<<<_B, _T>>>(_FSR_materials, _materials, false,
                                                 old_scalar_flux, old_sources);
    computeFSRScatterSourcesOnDevice<<<_B, _T>>>(_FSR_materials, _materials, false,
                                                 scalar_flux, new_sources);

    /* Reduce scatter sources across energy groups within each FSR */
    for (int e=0; e < _num_groups; e++) {
      strided_range<Iterator> old_strider(old_sources_vec.begin() + e,
                                          old_sources_vec.end(), _num_groups);
      strided_range<Iterator> new_strider(new_sources_vec.begin() + e,
                                          new_sources_vec.end(), _num_groups);
      thrust::transform(FSR_old_src.begin(), FSR_old_src.end(),
                        old_strider.begin(), FSR_old_src.begin(),
                        thrust::plus<FP_PRECISION>());
      thrust::transform(FSR_new_src.begin(), FSR_new_src.end(),
                        new_strider.begin(), FSR_new_src.begin(),
                        thrust::plus<FP_PRECISION>());
    }

    /* Compute the relative total source change in each FSR */
    thrust::transform(FSR_new_src.begin(), FSR_new_src.end(),
                      FSR_old_src.begin(), residuals.begin(),
                      thrust::minus<FP_PRECISION>());
    thrust::transform(residuals.begin(), residuals.end(),
                      FSR_old_src.begin(), residuals.begin(),
                      thrust::divides<FP_PRECISION>());

    /* Deallocate memory for Thrust vectors */
    old_sources_vec.clear();
    new_sources_vec.clear();
    FSR_old_src.clear();
    FSR_new_src.clear();
  }

  /* Replace INF and NaN values (from divide by zero) with 0. */
  thrust::replace_if(residuals.begin(), residuals.end(), inf_test, 0);
  thrust::replace_if(residuals.begin(), residuals.end(), nan_test, 0);

  /* Square the residuals */
  thrust::transform(residuals.begin(), residuals.end(),
                    residuals.begin(), residuals.begin(),
                    thrust::multiplies<double>());

  /* Sum up the residuals */
  residual = thrust::reduce(residuals.begin(), residuals.end());

  /* Deallocate memory for residuals vector */
  residuals.clear();

  /* Normalize the residual */
  residual = sqrt(residual / norm);

  return residual;
}


/**
 * @brief Computes the volume-averaged, energy-integrated nu-fission rate in
 *        each FSR and stores them in an array indexed by FSR ID.
 * @details This is a helper method for SWIG to allow users to retrieve
 *          FSR nu-fission rates as a NumPy array. An example of how this method
 *          can be called from Python is as follows:
 *
 * @code
 *          num_FSRs = geometry.getNumFSRs()
 *          fission_rates = solver.computeFSRFissionRates(num_FSRs)
 * @endcode
 *
 * @param fission_rates an array to store the nu-fission rates (implicitly
 *                      passed in as a NumPy array from Python)
 * @param num_FSRs the number of FSRs passed in from Python
 */
void GPUSolver::computeFSRFissionRates(double* fission_rates, int num_FSRs) {

  log_printf(INFO, "Computing FSR fission rates...");

  /* Allocate memory for the FSR nu-fission rates on the device and host */
  FP_PRECISION* dev_fission_rates;
  cudaMalloc((void**)&dev_fission_rates, _num_FSRs * sizeof(FP_PRECISION));
  FP_PRECISION* host_fission_rates = new FP_PRECISION[_num_FSRs];

  FP_PRECISION* scalar_flux =
       thrust::raw_pointer_cast(&_scalar_flux[0]);

  /* Compute the FSR nu-fission rates on the device */
  computeFSRFissionRatesOnDevice<<<_B, _T>>>(_FSR_volumes, _FSR_materials,
                                             _materials, scalar_flux,
                                             dev_fission_rates);

  /* Copy the nu-fission rate array from the device to the host */
  cudaMemcpy((void*)host_fission_rates, (void*)dev_fission_rates,
             _num_FSRs * sizeof(FP_PRECISION), cudaMemcpyDeviceToHost);

  /* Populate the double precision NumPy array for the output */
  for (int i=0; i < _num_FSRs; i++)
    fission_rates[i] = host_fission_rates[i];

  /* Deallocate the memory assigned to store the fission rates on the device */
  cudaFree(dev_fission_rates);
  delete [] host_fission_rates;
}
