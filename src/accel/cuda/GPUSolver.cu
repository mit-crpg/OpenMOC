#include "GPUSolver.h"

/** The number of azimuthal angles */
__constant__ int num_azim[1];

/** The number of energy groups */
__constant__ int num_groups[1];

/** The number of FSRs */
__constant__ int num_FSRs[1];

/** The number of polar angles */
__constant__ int num_polar[1];

/** Twice the number of polar angles */
__constant__ int two_times_num_polar[1];

/** The number of polar angles times energy groups */
__constant__ int polar_times_groups[1];

/** An array for the sines of the polar angle in the polar Quadrature set */
__constant__ FP_PRECISION sinthetas[MAX_POLAR_ANGLES];

/** An array of the weights for the polar angles from the Quadrature set */
__constant__ FP_PRECISION polar_weights[MAX_POLAR_ANGLES*MAX_AZIM_ANGLES];

/** A pointer to an array with the number of tracks per azimuthal angle */
__constant__ int num_tracks[MAX_AZIM_ANGLES/2];

/** The total number of Tracks */
__constant__ int tot_num_tracks[1];

/** A boolean indicating whether or not to use linear interpolation
 *  to comptue the exponential in the transport equation */
__constant__ bool interpolate_exponential[1];

/** The maximum index of the exponential linear interpolation table */
__constant__ int exp_table_max_index[1];

/** The spacing for the exponential linear interpolation table */
__constant__ FP_PRECISION exp_table_spacing[1];

/** The inverse spacing for the exponential linear interpolation table */
__constant__ FP_PRECISION inverse_exp_table_spacing[1];



/**
 * @brief Fast method to round a single precision floating point value
 *        to an integer on the GPU.
 * @param x float floating point value to round
 * @return the rounded down integer value
 */
__device__ int round_to_int(float x) {
  return __float2int_rd(x);
}


/**
 * @brief Fast method to round a double precision floating point value
 *        to an integer on the GPU.
 * @param x double floating point value to round
 * @return the rounded down integer value
 */
__device__ int round_to_int(double x) {
  return __double2int_rd(x);
}


/**
 * @brief Compute the total fission source from all FSRs on the GPU.
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
  FP_PRECISION volume;
  FP_PRECISION source;

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

  return;
}


/**
 * @brief Normalizes all FSR scalar fluxes and Track boundary angular
 *        fluxes to the total fission source (times \f$ \nu \f$).
 * @param scalar_flux an array of the FSR scalar fluxes
 * @param boundary_flux an array of the Track boundary fluxes
 * @param norm_factor the normalization factor
 */
__global__ void normalizeFluxesOnDevice(FP_PRECISION* scalar_flux,
                                        FP_PRECISION* boundary_flux,
                                        FP_PRECISION norm_factor) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  /* Normalize scalar fluxes for each FSR */
  while(tid < *num_FSRs) {

    for (int e=0; e < *num_groups; e++)
      scalar_flux(tid,e) *= norm_factor;

    tid += blockDim.x * gridDim.x;
  }

  tid = threadIdx.x + blockIdx.x * blockDim.x;

  /* Normalize angular boundary fluxes for each Track */
  while(tid < *tot_num_tracks) {

    for (int pe2=0; pe2 < 2*(*polar_times_groups); pe2++)
      boundary_flux(tid,pe2) *= norm_factor;

    tid += blockDim.x * gridDim.x;
  }

  return;
}


/**
 * @brief Computes the total source (fission and scattering) in each FSR
 *        on the GPU.
 * @details This method computes the total source in each region based on
 *          this iteration's current approximation to the scalar flux. A
 *          residual for the source with respect to the source compute on
 *          the previous iteration is computed and returned. The residual
 *          is determined as follows:
 *          \f$ res = \sqrt{\frac{\displaystyle\sum \displaystyle\sum
 *                    \left(\frac{Q^i - Q^{i-1}{Q^i}\right)^2}{# FSRs}}} $\f
 *
 * @param FSR_materials an array of FSR Material indices
 * @param materials an array of dev_material pointers
 * @param scalar_flux an array of FSR scalar fluxes
 * @param old_fission_sources an array of current FSR sources from previous iteration
 * @param reduced_sources an array of FSR sources / total xs
 * @param inverse_k_eff the inverse of keff
 * @param source_residuals an array of the FSR source residuals
 * @return the residual between this source and the previous source
 */
__global__ void computeFSRSourcesOnDevice(int* FSR_materials,
                                          dev_material* materials,
                                          FP_PRECISION* scalar_flux,
                                          FP_PRECISION* old_fission_sources,
                                          FP_PRECISION* reduced_sources,
                                          FP_PRECISION inverse_k_eff,
                                          FP_PRECISION* source_residuals) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  /* Reset the residual for the old and new fission sources to zero */
  source_residuals[threadIdx.x + blockIdx.x * blockDim.x] = 0.0;

  FP_PRECISION fission_source;
  FP_PRECISION scatter_source;
  FP_PRECISION fsr_fission_source;

  dev_material* curr_material;

  FP_PRECISION* nu_sigma_f;
  FP_PRECISION* sigma_s;
  FP_PRECISION* sigma_t;
  FP_PRECISION* chi;

  /* Iterate over all FSRs */
  while (tid < *num_FSRs) {

    curr_material = &materials[FSR_materials[tid]];

    nu_sigma_f = curr_material->_nu_sigma_f;
    sigma_s = curr_material->_sigma_s;
    sigma_t = curr_material->_sigma_t;
    chi = curr_material->_chi;

    /* Initialize the fission source to zero for this FSR */
    fission_source = 0;
    fsr_fission_source = 0.0;

    /* Compute total fission source for current region */
    for (int e=0; e < *num_groups; e++)
      fission_source += scalar_flux(tid,e) * nu_sigma_f[e];

    /* Compute total scattering source for this FSR in group G */
    for (int G=0; G < *num_groups; G++) {

      scatter_source = 0;

      for (int g=0; g < *num_groups; g++)
        scatter_source += sigma_s[G*(*num_groups)+g] * scalar_flux(tid,g);

      /* Set the fission source for FSR r in group G */
      fsr_fission_source += fission_source * chi[G] * inverse_k_eff;

      reduced_sources(tid,G) = __fdividef((inverse_k_eff * fission_source 
                               * chi[G] + scatter_source) * ONE_OVER_FOUR_PI, 
                               sigma_t[G]);
    }

    /* Compute the norm of residuals of the sources for convergence */
    if (fabs(fsr_fission_source) > 1E-10)
      source_residuals[threadIdx.x + blockIdx.x * blockDim.x] +=
                      pow((fsr_fission_source - old_fission_sources[tid]) 
                      / fsr_fission_source, 2);


    /* Update the old fission source */
    old_fission_sources[tid] = fsr_fission_source;

    /* Increment the thread id */
    tid += blockDim.x * gridDim.x;
  }

  return;
}


/**
 * @brief Compute the total fission source from all FSRs and energy groups
 *        on the GPU.
 * @param FSR_volumes an array of the FSR volumes
 * @param FSR_materials an array of the FSR Material indices
 * @param materials an array of the dev_material pointers
 * @param scalar_flux an array of FSR scalar fluxes
 * @param total array of FSR total reaction rates
 * @param fission an array of FSR fission rates
 * @param scatter an array of FSR scattering rates
 */
__global__ void computeKeffReactionRates(FP_PRECISION* FSR_volumes,
                                         int* FSR_materials,
                                         dev_material* materials,
                                         FP_PRECISION* scalar_flux,
                                         FP_PRECISION* total,
                                         FP_PRECISION* fission,
                                         FP_PRECISION* scatter) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  dev_material* curr_material;
  FP_PRECISION* sigma_t;
  FP_PRECISION* nu_sigma_f;
  FP_PRECISION* sigma_s;
  FP_PRECISION volume;

  FP_PRECISION tot = 0.;
  FP_PRECISION fiss = 0.;
  FP_PRECISION scatt = 0.;

  /* Iterate over all FSRs */
  while (tid < *num_FSRs) {

    curr_material = &materials[FSR_materials[tid]];
    sigma_t = curr_material->_sigma_t;
    nu_sigma_f = curr_material->_nu_sigma_f;
    sigma_s = curr_material->_sigma_s;
    volume = FSR_volumes[tid];

    FP_PRECISION curr_tot = 0.;
    FP_PRECISION curr_fiss = 0.;
    FP_PRECISION curr_scatt = 0.;

    /* Iterate over all energy groups and update total and fission
     * rates for this thread block */
    for (int e=0; e < *num_groups; e++) {
      curr_tot += sigma_t[e] * scalar_flux(tid,e);
      curr_fiss += nu_sigma_f[e] * scalar_flux(tid,e);
    }

    tot += curr_tot * volume;
    fiss += curr_fiss * volume;

    /* Iterate over all energy groups and update scattering
     * rates for this thread block */
    for (int G=0; G < *num_groups; G++) {
      for (int g=0; g < *num_groups; g++)
        curr_scatt += sigma_s[G*(*num_groups)+g] * scalar_flux(tid,g);
    }

    scatt += curr_scatt * volume;

    /* Increment thread id */
    tid += blockDim.x * gridDim.x;
  }

  /* Copy this thread's total and scatter rates to global memory */
  tid = threadIdx.x + blockIdx.x * blockDim.x;
  total[tid] = tot;
  fission[tid] = fiss;
  scatter[tid] = scatt;

  return;
}


/**
 * @brief Perform an atomic addition in double precision to an array address
 *        on the GPU.
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
 * @brief Computes the exponential term in the transport equation for a
 *        Track segment on the GPU.
 * @details This method computes \f$ 1 - exp(-l\Sigma^T_g/sin(\theta_p)) \f$
 *          for a segment with total group cross-section and for
 *          some polar angle.
 * @param sigma_t the total group cross-section at this energy
 * @param length the length of the line segment projected in the xy-plane
 * @param _exp_table the exponential linear interpolation table
 * @param p the polar angle index
 * @return the evaluated exponential
 */
__device__ FP_PRECISION computeExponential(FP_PRECISION sigma_t,
                                           FP_PRECISION length,
                                           FP_PRECISION* _exp_table,
                                           int p) {

  FP_PRECISION exponential;
  FP_PRECISION tau = sigma_t * length;

  /* Evaluate the exponential using the linear interpolation table */
  if (*interpolate_exponential) {
    int index;

    index = round_to_int(tau * (*inverse_exp_table_spacing));
    index *= (*two_times_num_polar);
    exponential = (1. - (_exp_table[index+2 * p] * tau +
                         _exp_table[index + 2 * p +1]));
  }

  /* Evalute the exponential using the intrinsic exp(...) function */
  else {
    FP_PRECISION sintheta = sinthetas[p];
    #ifdef SINGLE
    exponential = 1.0 - __expf(- tau / sintheta);
    #else
    exponential = 1.0 - exp(- tau / sintheta);
    #endif
  }

  return exponential;
}


/**
 * @brief Computes the contribution to the FSR scalar flux from a Track segment
 *        in a single energy group on the GPU.
 * @details This method integrates the angular flux for a Track segment across
 *        energy groups and polar angles, and tallies it into the FSR scalar
 *        flux, and updates the Track's angular flux.
 * @param curr_segment a pointer to the Track segment of interest
 * @param azim_index a pointer to the azimuthal angle index for this segment
 * @param energy_group the energy group of interest
 * @param materials the array of dev_material pointers
 * @param track_flux a pointer to the Track's angular flux
 * @param reduced_sources the array of FSR sources / total xs
 * @param polar_weights the array of polar Quadrature weights
 * @param _exp_table the exponential interpolation table
 * @param scalar_flux the array of FSR scalar fluxes
 */
__device__ void scalarFluxTally(dev_segment* curr_segment,
                                int azim_index,
                                int energy_group,
                                dev_material* materials,
                                FP_PRECISION* track_flux,
                                FP_PRECISION* reduced_sources,
                                FP_PRECISION* polar_weights,
                                FP_PRECISION* _exp_table,
                                FP_PRECISION* scalar_flux) {

  int fsr_id = curr_segment->_region_uid;
  FP_PRECISION length = curr_segment->_length;
  dev_material* curr_material = &materials[curr_segment->_material_index];
  FP_PRECISION *sigma_t = curr_material->_sigma_t;

  /* The change in angular flux long this Track segment in this FSR */
  FP_PRECISION delta_psi;
  FP_PRECISION exponential;

  /* Zero the FSR scalar flux contribution from this segment and energy group */
  FP_PRECISION fsr_flux = 0.0;

  /* Compute the exponential interpolation table index */

  /* Loop over polar angles */
  for (int p=0; p < *num_polar; p++) {
    exponential = computeExponential(sigma_t[energy_group],
                                     length, _exp_table, p);
    delta_psi = (track_flux[p] - reduced_sources(fsr_id,energy_group)) *
               exponential;
    fsr_flux += delta_psi * polar_weights(azim_index,p);
    track_flux[p] -= delta_psi;
  }

  /* Atomically increment the scalar flux for this FSR */
  atomicAdd(&scalar_flux(fsr_id,energy_group), fsr_flux);
}


/**
 * @brief Updates the boundary flux for a Track given boundary conditions
 *        on the GPU.
 * @details For reflective boundary conditions, the outgoing boundary flux
 *          for the Track is given to the reflecting track. For vacuum
 *          boundary conditions, the outgoing flux tallied as leakage.
 *          Note: Only one energy group is transferred by this routine.
 * @param curr_track a pointer to the Track of interest
 * @param azim_index a pointer to the azimuthal angle index for this segment
 * @param track_flux an array of the outgoing Track flux
 * @param boundary_flux an array of all angular fluxes
 * @param leakage an array of leakages for each CUDA thread
 * @param polar_weights an array of polar Quadrature weights
 * @param energy_angle_index the energy group index
 * @param direction the Track direction (forward - true, reverse - false)
 */
__device__ void transferBoundaryFlux(dev_track* curr_track,
                                     int azim_index,
                                     FP_PRECISION* track_flux,
                                     FP_PRECISION* boundary_flux,
                                     FP_PRECISION* leakage,
                                     FP_PRECISION* polar_weights,
                                     int energy_angle_index,
                                     bool direction) {

  int start = energy_angle_index;
  bool bc;
  int track_out_id;

  /* Extract boundary conditions for this Track and the pointer to the
   * outgoing reflective Track, and index into the leakage array */

  /* For the "forward" direction */
  if (direction) {
    bc = curr_track->_bc_out;
    track_out_id = curr_track->_track_out;
    start += curr_track->_refl_out * (*polar_times_groups);
  }

  /* For the "reverse" direction */
  else {
    bc = curr_track->_bc_in;
    track_out_id = curr_track->_track_in;
    start += curr_track->_refl_in * (*polar_times_groups);
  }

  FP_PRECISION* track_out_flux = &boundary_flux(track_out_id,start);

  /* Put Track's flux in the shared memory temporary flux array */
  for (int p=0; p < *num_polar; p++) {
    track_out_flux[p] = track_flux[p] * bc;
    leakage[0] += track_flux[p] * polar_weights(azim_index,p) * (!bc);
  }
}


/**
 * @brief This method performs one transport sweep of one halfspace of all
 *        azimuthal angles, tracks, segments, polar angles and energy groups
 *        on the GPU.
 * @details The method integrates the flux along each track and updates the
 *          boundary fluxes for the corresponding output Track, while updating
 *          the scalar flux in each FSR.
 * @param scalar_flux an array of FSR scalar fluxes
 * @param boundary_flux an array of Track boundary fluxes
 * @param reduced_sources an array of FSR sources / total xs
 * @param leakage an array of angular flux leakaages
 * @param materials an array of dev_material pointers
 * @param tracks an array of Tracks
 * @param _exp_table an array for the exponential interpolation table
 * @param tid_offset the Track offset for azimuthal angle halfspace
 * @param tid_max the upper bound on the Track IDs for this azimuthal
 *                angle halfspace
 */
__global__ void transportSweepOnDevice(FP_PRECISION* scalar_flux,
                                       FP_PRECISION* boundary_flux,
                                       FP_PRECISION* reduced_sources,
                                       FP_PRECISION* leakage,
                                       dev_material* materials,
                                       dev_track* tracks,
                                       FP_PRECISION* _exp_table,
                                       int tid_offset,
                                       int tid_max) {

  /* Shared memory buffer for each thread's angular flux */
  extern __shared__ FP_PRECISION temp_flux[];
  FP_PRECISION* track_flux;

  int tid = tid_offset + threadIdx.x + blockIdx.x * blockDim.x;
  int track_id = tid / *num_groups;
  int track_flux_index = threadIdx.x * (*two_times_num_polar);
  int energy_group = tid % (*num_groups);
  int energy_angle_index = energy_group * (*num_polar);

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
    for (int p=0; p < *num_polar; p++) {

      /* Forward flux along this Track */
      track_flux[p] = boundary_flux(track_id,p+energy_angle_index);

      /* Reverse flux along this Track */
      track_flux[(*num_polar) + p] =
            boundary_flux(track_id,p+energy_angle_index+(*polar_times_groups));
    }

    /* Loop over each Track segment in forward direction */
    for (int i=0; i < num_segments; i++) {
      curr_segment = &curr_track->_segments[i];
      scalarFluxTally(curr_segment, azim_index, energy_group, materials,
                      track_flux, reduced_sources, polar_weights,
                      _exp_table, scalar_flux);
    }

    /* Transfer boundary angular flux to outgoing Track */
    transferBoundaryFlux(curr_track, azim_index, track_flux, boundary_flux,
                         &leakage[threadIdx.x + blockIdx.x * blockDim.x],
                         polar_weights, energy_angle_index, true);

    /* Loop over each Track segment in reverse direction */
    track_flux = &temp_flux[track_flux_index + (*num_polar)];

    for (int i=num_segments-1; i > -1; i--) {
      curr_segment = &curr_track->_segments[i];
      scalarFluxTally(curr_segment, azim_index, energy_group, materials,
                      track_flux, reduced_sources, polar_weights,
                      _exp_table, scalar_flux);
  }

    /* Transfer boundary angular flux to outgoing Track */
    transferBoundaryFlux(curr_track, azim_index, track_flux, boundary_flux,
                         &leakage[threadIdx.x + blockIdx.x * blockDim.x],
                         polar_weights, energy_angle_index, false);

    /* Update the indices for this thread to the next Track, energy group */
    tid += blockDim.x * gridDim.x;
    track_id = tid / *num_groups;
    energy_group = tid % (*num_groups);
    energy_angle_index = energy_group * (*num_polar);
  }

  return;
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
      scalar_flux(tid,i) *= 0.5;
      scalar_flux(tid,i) = FOUR_PI * reduced_sources(tid,i) +
                           __fdividef(scalar_flux(tid,i), 
                           (sigma_t[i] * volume));
    }

    /* Increment thread id */
    tid += blockDim.x * gridDim.x;
  }

    return;
}


/**
 * @brief Computes the volume-weighted, energy integrated fission rate in
 *        each FSR and stores them in an array indexed by FSR ID on the GPU.
 * @details This is a helper method for the
 *          GPUSolver::computeFSRFissionRates(...) method.
 * @param fission_rates an array to store the fission rates
 * @param fission_rates an array in which to store the FSR fission rates
 * @param FSR_materials an array of FSR material indices
 * @param materials an array of dev_material pointers
 * @param scalar_flux an array of FSR scalar fluxes
 */
__global__ void computeFSRFissionRatesOnDevice(double* fission_rates,
                                               int* FSR_materials,
                                               dev_material* materials,
                                               FP_PRECISION* scalar_flux) {

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  dev_material* curr_material;
  FP_PRECISION* sigma_f;

  /* Loop over all FSRs and compute the volume-weighted fission rate */
  while (tid < *num_FSRs) {

    curr_material = &materials[FSR_materials[tid]];
    sigma_f = curr_material->_sigma_f;

    /* Initialize the fission rate for this FSR to zero */
    fission_rates[tid] = 0.0;

    for (int i=0; i < *num_groups; i++)
      fission_rates[tid] += sigma_f[i] * scalar_flux(tid,i);

    /* Increment thread id */
    tid += blockDim.x * gridDim.x;
  }

  return;
}


/**
 * @brief Constructor initializes arrays for dev_tracks and dev_materials..
 * @details The constructor retrieves the number of energy groups and FSRs
 *          and azimuthal angles from the Geometry and TrackGenerator if
 *          passed in as parameters by the user. The constructor initalizes
 *          the number of CUDA threads and thread blocks each to a default
 *          of 64.
 * @param geometry an optional pointer to the Geometry
 * @param track_generator an optional pointer to the TrackjGenerator
 */
GPUSolver::GPUSolver(Geometry* geometry, TrackGenerator* track_generator) :

  Solver(geometry, track_generator) {

  /* The default number of thread blocks and threads per thread block */
  _B = 64;
  _T = 64;

  _materials = NULL;
  _dev_tracks = NULL;
  _FSR_materials = NULL;

  _total = NULL;
  _fission = NULL;
  _scatter = NULL;
  _leakage = NULL;

  if (geometry != NULL)
    setGeometry(geometry);

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

  if (_boundary_flux != NULL) {
    cudaFree(_boundary_flux);
    _boundary_flux = NULL;
  }

  if (_scalar_flux != NULL) {
    cudaFree(_scalar_flux);
    _scalar_flux = NULL;
  }

  if (_old_fission_sources != NULL) {
    cudaFree(_old_fission_sources);
    _old_fission_sources = NULL;
  }

  if (_reduced_sources != NULL) {
    cudaFree(_reduced_sources);
    _reduced_sources = NULL;
  }

  if (_fission_sources != NULL) {
    _fission_sources_vec.clear();
    _fission_sources = NULL;
  }

  if (_total != NULL) {
    _total_vec.clear();
    _total = NULL;
  }

  if (_fission != NULL) {
    _fission_vec.clear();
    _fission = NULL;
  }

  if (_scatter != NULL) {
    _scatter_vec.clear();
    _scatter = NULL;
  }

  if (_source_residuals != NULL) {
    _source_residuals_vec.clear();
    _source_residuals = NULL;
  }

  if (_leakage != NULL) {
    _leakage_vec.clear();
    _leakage = NULL;
  }

  if (_exp_table != NULL) {
    cudaFree(_exp_table);
    _exp_table = NULL;
  }
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
 * @brief Returns the FSR scalar flux for some energy group.
 * @param fsr_id the ID for the FSR of interest
 * @param energy_group the energy group of interest
 */
FP_PRECISION GPUSolver::getFSRScalarFlux(int fsr_id, int energy_group) {

  /* Error checking */
  if (fsr_id >= _num_FSRs)
    log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d in energy"
               " group %d since the solver only contains FSR with IDs greater "
               "than or equal to %d", fsr_id, energy_group, _num_FSRs-1);

  if (fsr_id < 0)
    log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d in energy"
               " group %d since FSRs do not have negative IDs",
               fsr_id, energy_group);

  if (energy_group-1 >= _num_groups)
    log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d in energy"
               " group %d since the solver only has %d energy groups",
               fsr_id, energy_group, _num_groups);

  if (energy_group <= 0)
    log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d in energy"
               " group %d since energy groups are greater than 1",
               fsr_id, energy_group);

  /* Copy the scalar flux for this FSR and energy group from the device */
  FP_PRECISION fsr_scalar_flux;
  int flux_index = fsr_id * _num_groups + (energy_group - 1);
  cudaMemcpy((void*)&fsr_scalar_flux, (void*)&_scalar_flux[flux_index],
             sizeof(FP_PRECISION), cudaMemcpyDeviceToHost);

  return fsr_scalar_flux;
}


/**
 * @brief Return the scalar flux array indexed by FSR IDs and energy groups.
 *        which contains the corresponding fluxes for each flat source region.
 * @return an array of FSR scalar fluxes
 */
FP_PRECISION* GPUSolver::getFSRScalarFluxes() {

  if (_scalar_flux == NULL)
    log_printf(ERROR, "Unable to returns the GPUSolver's scalar flux "
               "array since it has not yet been allocated in memory");

  /* Copy the scalar flux for all FSRs from the device to the host */
  FP_PRECISION* fsr_scalar_fluxes = new FP_PRECISION[_num_FSRs * _num_groups];
  cudaMemcpy((void*)fsr_scalar_fluxes, (void*)_scalar_flux,
             _num_FSRs * _num_groups * sizeof(FP_PRECISION),
             cudaMemcpyDeviceToHost);

  return fsr_scalar_fluxes;
}


/**
 * @brief Returns the FSR source for some energy group.
 * @param fsr_id the ID for the FSR of interest
 * @param energy_group the energy group of interest
 */
FP_PRECISION GPUSolver::getFSRSource(int fsr_id, int energy_group) {

  /* Error checking */
  if (fsr_id >= _num_FSRs)
    log_printf(ERROR, "Unable to return a source for FSR id = %d in energy"
      " group %d since the solver only contains FSR with IDs greater than "
      "or equal to %d", fsr_id, energy_group, _num_FSRs-1);

  if (fsr_id < 0)
    log_printf(ERROR, "Unable to return a source for FSR id = %d in energy"
               " group %d since FSRs do not have negative IDs",
               fsr_id, energy_group);

  if (energy_group-1 >= _num_groups)
    log_printf(ERROR, "Unable to return a source for FSR id = %d in energy"
              " group %d since the solver only has %d energy groups",
               fsr_id, energy_group, _num_groups);

  if (energy_group <= 0)
    log_printf(ERROR, "Unable to return a source for FSR id = %d in energy"
               " group %d since energy groups are greater than 1",
               fsr_id, energy_group);

  /* Get host material */
  Material* host_material = _geometry->findFSRMaterial(fsr_id);

  /* Get cross sections and scalar flux */
  FP_PRECISION* nu_sigma_f = host_material->getNuSigmaF();
  FP_PRECISION* sigma_s = host_material->getSigmaS();
  FP_PRECISION* chi = host_material->getChi();
  FP_PRECISION* fsr_scalar_fluxes = new FP_PRECISION[_num_groups];
  cudaMemcpy((void*)fsr_scalar_fluxes, (void*)&_scalar_flux[fsr_id*_num_groups],
             _num_groups * sizeof(FP_PRECISION),
             cudaMemcpyDeviceToHost);

  /* Initialize variables */
  FP_PRECISION fission_source = 0.0;
  FP_PRECISION scatter_source = 0.0;
  FP_PRECISION total_source;

  /* Compute total fission source for current region */
  for (int e=0; e < _num_groups; e++){
    fission_source += fsr_scalar_fluxes[e] * nu_sigma_f[e];
  }

  fission_source /= _k_eff;

  /* Compute total scattering source for this FSR */
  for (int g=0; g < _num_groups; g++){
    scatter_source += sigma_s[(energy_group-1)*(_num_groups)+g] 
                    * fsr_scalar_fluxes[g];
  }

  /* Compute the total source */
  total_source = (fission_source * chi[energy_group-1] + scatter_source) *
      ONE_OVER_FOUR_PI;

  delete [] fsr_scalar_fluxes;

  return total_source;
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
 * @brief Sets the Geometry pointer for the GPUSolver.
 * @details The Geometry must already have initialized FSR offset maps
 *          and segmentized the TrackGenerator's tracks. Each of these
 *          should be initiated in Python prior to assigning a Geometry
 *          to the GPUSolver:
 *
 * @code
 *          geometry.initializeFlatSourceRegions()
 *          track_generator.generateTracks()
 * @endcode
 *
 * @param geometry a pointer to a Geometry
 */
void GPUSolver::setGeometry(Geometry* geometry) {

  Solver::setGeometry(geometry);

  initializeMaterials();

  /* Copy the number of energy groups to constant memory on the GPU */
  cudaMemcpyToSymbol(num_groups, (void*)&_num_groups, sizeof(int), 0,
                     cudaMemcpyHostToDevice);
}


/**
 * @brief Sets the TrackGenerator with characteristic tracks for the GPUSolver.
 * @details The TrackGenerator must already have generated Tracks and have
 *          used ray tracing to segmentize them across the Geometry. This
 *          should be initated in Python prior to assigning the TrackGenerator
 *          to the GPUSolver:
 *
 * @code
 *          track_generator.generateTracks()
 * @endcode
 *
 * @param track_generator a pointer to a TrackGenerator
 */
void GPUSolver::setTrackGenerator(TrackGenerator* track_generator) {
  Solver::setTrackGenerator(track_generator);
  initializeTracks();
}


/**
 * @brief Creates a polar Quadrature object for the GPUSolver on the GPU.
 */
void GPUSolver::initializePolarQuadrature() {

  log_printf(INFO, "Initializing polar quadrature on the GPU...");

  Solver::initializePolarQuadrature();

  if (_num_polar > MAX_POLAR_ANGLES)
    log_printf(ERROR, "Unable to initialize a polar quadrature with %d "
               "angles for the GPUSolver which is limited to %d polar "
               "angles. Update the MAX_POLAR_ANGLES macro in GPUSolver.h "
               "and recompile.", _num_polar, MAX_POLAR_ANGLES);

  /* Copy the number of polar angles to constant memory on the GPU */
  cudaMemcpyToSymbol(num_polar, (void*)&_num_polar, sizeof(int), 0,
                     cudaMemcpyHostToDevice);

  /* Copy twice the number of polar angles to constant memory on the GPU */
  cudaMemcpyToSymbol(two_times_num_polar, (void*)&_two_times_num_polar,
                     sizeof(int), 0, cudaMemcpyHostToDevice);

  /* Copy the number of polar angles times energy groups to constant memory
   * on the GPU */
  cudaMemcpyToSymbol(polar_times_groups, (void*)&_polar_times_groups,
                     sizeof(int), 0, cudaMemcpyHostToDevice);

  /* Copy the polar weights to constant memory on the GPU */
  cudaMemcpyToSymbol(polar_weights, (void*)_polar_weights,
      _num_polar * _num_azim * sizeof(FP_PRECISION), 0, cudaMemcpyHostToDevice);
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
    free(host_FSR_volumes);
    free(host_FSR_materials);

    /* Free the temporary array of FSRs to material indices on the host */
    free(FSRs_to_material_indices);
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the GPUSolver's FSRs "
               "on the device. Backtrace:%s", e.what());
  }

  initializeThrustVectors();
}


/**
 * @brief Allocates data on the GPU for all Materials data.
 * @details This method loops over the materials in the host_materials map.
 *          Since cuda does not support std::map data types on the device (GPU), 
 *          the materials map must be converted to an array and a map created
 *          that maps a material ID to an indice in the new materials array. In
 *          initializeTracks, this map is used to convert the Material ID
 *          associated with every segment to an index in the materials array.
 */
void GPUSolver::initializeMaterials() {

  log_printf(INFO, "Initializing materials on the GPU...");

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
    for (iter=host_materials.begin(); iter != host_materials.end(); ++iter){
      clone_material_on_gpu(iter->second, &_materials[material_index]);
      _material_IDs_to_indices[iter->second->getId()] = material_index;
      material_index++;
    }
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the GPUSolver's "
               "dev_materials. Backtrace:%s", e.what());
  }
}


/**
 * @brief Allocates memory on the GPU for all Tracks in the simulation.
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

      clone_track_on_gpu(_tracks[i], &_dev_tracks[i], _material_IDs_to_indices);

      /* Make Track reflective */
      index = computeScalarTrackIndex(_tracks[i]->getTrackInI(),
      _tracks[i]->getTrackInJ());
      cudaMemcpy((void*)&_dev_tracks[i]._track_in,
                 (void*)&index, sizeof(int), cudaMemcpyHostToDevice);

      index = computeScalarTrackIndex(_tracks[i]->getTrackOutI(),
      _tracks[i]->getTrackOutJ());
      cudaMemcpy((void*)&_dev_tracks[i]._track_out,
                 (void*)&index, sizeof(int), cudaMemcpyHostToDevice);
    }

    /* Copy the array of number of Tracks for each azimuthal angle into
     * constant memory on GPU */
    cudaMemcpyToSymbol(num_tracks, (void*)_num_tracks,
                       _num_azim * sizeof(int), 0, cudaMemcpyHostToDevice);

    /* Copy the total number of Tracks into constant memory on GPU */
    cudaMemcpyToSymbol(tot_num_tracks, (void*)&_tot_num_tracks,
                       sizeof(int), 0, cudaMemcpyHostToDevice);

    /* Copy the number of azimuthal angles into constant memory on GPU */
    cudaMemcpyToSymbol(num_azim, (void*)&_num_azim, sizeof(int), 0,
                       cudaMemcpyHostToDevice);

    /* Copy the array of number of Tracks for each azimuthal angles into
     * constant memory on GPU */
    cudaMemcpyToSymbol(num_tracks, (void*)_num_tracks,
                       _num_azim * sizeof(int), 0, cudaMemcpyHostToDevice);

    /* Copy the total number of Tracks into constant memory on GPU */
    cudaMemcpyToSymbol(tot_num_tracks, (void*)&_tot_num_tracks,
                       sizeof(int), 0, cudaMemcpyHostToDevice);
  }

  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the GPUSolver's "
               "dev_tracks on the device. Backtrace:%s", e.what());
  }
}


/**
 * @brief Allocates memory for Track boundary angular fluxes and leakages
 *        and FSR scalar fluxes on the GPU.
 * @details Deletes memory for old flux arrays if they were allocated for a
 *          previous simulation.
 */
void GPUSolver::initializeFluxArrays() {

  log_printf(INFO, "Initializing flux arrays on the GPU...");

  /* Delete old flux arrays if they exist */
  if (_boundary_flux != NULL)
    cudaFree(_boundary_flux);

  if (_scalar_flux != NULL)
    cudaFree(_scalar_flux);

  /* Allocate memory for all flux arrays on the device */
  try{
    cudaMalloc((void**)&_boundary_flux,
               2*_tot_num_tracks * _polar_times_groups*sizeof(FP_PRECISION));
    cudaMalloc((void**)&_scalar_flux,
               _num_FSRs * _num_groups * sizeof(FP_PRECISION));
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the GPUSolver's fluxes "
               "on the device. Backtrace:%s", e.what());
  }
}


/**
 * @brief Allocates memory for FSR source arrays on the GPU.
 * @details Deletes memory for old source arrays if they were allocated for a
 *          previous simulation.
 */
void GPUSolver::initializeSourceArrays() {

  log_printf(INFO, "Initializing source arrays on the GPU...");

  /* Delete old sources arrays if they exist */
  if (_old_fission_sources != NULL)
    cudaFree(_old_fission_sources);

  if (_reduced_sources != NULL)
    cudaFree(_reduced_sources);

  /* Allocate memory for all source arrays on the device */
  try{

    cudaMalloc((void**)&_old_fission_sources,
               _num_FSRs * sizeof(FP_PRECISION));

    cudaMalloc((void**)&_reduced_sources,
               _num_FSRs * _num_groups * sizeof(FP_PRECISION));
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the GPUSolver's FSR "
               "sources array on the device. Backtrace:%s", e.what());
  }
}


/**
 * @brief Initialize Thrust vectors for the fission and absorption rates,
 *        source residuals, leakage and fission sources.
 */
void GPUSolver::initializeThrustVectors() {

  log_printf(INFO, "Initializing Thrust vectors on the GPU...");

  /* Delete old vectors if they exist */
  if (_fission_sources != NULL) {
    _fission_sources = NULL;
    _fission_sources_vec.clear();
  }

  if (_total != NULL) {
    _total = NULL;
    _total_vec.clear();
  }

  if (_fission != NULL) {
    _fission = NULL;
    _fission_vec.clear();
  }

  if (_scatter != NULL) {
    _scatter = NULL;
    _scatter_vec.clear();
  }

  if (_source_residuals != NULL) {
    _source_residuals = NULL;
    _source_residuals_vec.clear();
  }

  if (_leakage != NULL) {
    _leakage = NULL;
    _leakage_vec.clear();
  }

  /* Allocate memory for fission, absorption and source vectors on device */
  try{
    /* Allocate fission source array on device */
    _fission_sources_vec.resize(_B * _T);
    _fission_sources = thrust::raw_pointer_cast(&_fission_sources_vec[0]);

    /* Allocate total reaction rate array on device */
    _total_vec.resize(_B * _T);
    _total = thrust::raw_pointer_cast(&_total_vec[0]);

    /* Allocate fission reaction rate array on device */
    _fission_vec.resize(_B * _T);
    _fission = thrust::raw_pointer_cast(&_fission_vec[0]);

    /* Allocate scattering reaction rate array on device */
    _scatter_vec.resize(_B * _T);
    _scatter = thrust::raw_pointer_cast(&_scatter_vec[0]);

    /* Allocate source residual array on device */
    _source_residuals_vec.resize(_B * _T);
    _source_residuals = thrust::raw_pointer_cast(&_source_residuals_vec[0]);

    /* Allocate leakage array on device */
    _leakage_vec.resize(_B * _T);
    _leakage = thrust::raw_pointer_cast(&_leakage_vec[0]);
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the GPUSolver's "
               "Thrust vectors. Backtrace:%s", e.what());
  }
}


/**
 * @brief This method computes the index for the Track j at azimuthal angle i.
 * @details This method is necessary since the array of dev_tracks on the device
 *          is a 1D array which needs a one-to-one mapping from the 2D jagged
 *          array of Tracks on the host.
 * @param i azimuthal angle number
 * @param j the jth track at angle i
 * @return an index into the device track array
 */
int GPUSolver::computeScalarTrackIndex(int i, int j) {

  int index =0;
  int p = 0;

  /* Iterate over each azimuthal angle and increment index by the number of
   * Tracks at each angle */
  while (p < i) {
    index += _num_tracks[p];
    p++;
  }

  /* Update index for this Track since it is the jth Track at angle i */
  index += j;

  return index;
}


/**
 * @brief Builds a linear interpolation table to compute exponentials for
 *        each segment of each Track for each polar angle on the GPU.
 */
void GPUSolver::buildExpInterpTable(){

  log_printf(INFO, "Building exponential interpolation table on device...");

  /* Copy a boolean indicating whether or not to use the linear interpolation
   * table or the exp intrinsic function */
  cudaMemcpyToSymbol(interpolate_exponential,(void*)&_interpolate_exponential,
                     sizeof(bool), 0, cudaMemcpyHostToDevice);

  /* Copy the sines of the polar angles which is needed if the user
   * requested the use of the exp intrinsic to evaluate exponentials */
  cudaMemcpyToSymbol(sinthetas, (void*)_polar_quad->getSinThetas(),
                     _num_polar * sizeof(FP_PRECISION), 0,
                     cudaMemcpyHostToDevice);

  /* Find largest optical path length track segment */
  FP_PRECISION tau = _track_generator->getMaxOpticalLength();

  /* Expand tau slightly to accomodate track segments which have a
   * length very nearly equal to the maximum value */
  tau *= 1.01;

  /* Set size of interpolation table */
  int num_array_values =
          tau * sqrt(1. / (8. * _source_convergence_thresh * 1e-2));
  _exp_table_spacing = tau / num_array_values;
  _inverse_exp_table_spacing = 1.0 / _exp_table_spacing;
  _exp_table_size = _two_times_num_polar * num_array_values;
  _exp_table_max_index = _exp_table_size - _two_times_num_polar - 1;

  /* Allocate array for the table */
  if (_exp_table != NULL)
    delete [] _exp_table;

  FP_PRECISION* exp_table = new FP_PRECISION[_exp_table_size];

  FP_PRECISION expon;
  FP_PRECISION intercept;
  FP_PRECISION slope;

  /* Create exponential interpolation table */
  for (int i = 0; i < num_array_values; i ++){
    for (int p = 0; p < _num_polar; p++){
      expon = exp(- (i * _exp_table_spacing) / _polar_quad->getSinTheta(p));
      slope = - expon / _polar_quad->getSinTheta(p);
      intercept = expon * (1 + (i * _exp_table_spacing) / 
                  _polar_quad->getSinTheta(p));
      exp_table[_two_times_num_polar * i + 2 * p] = slope;
      exp_table[_two_times_num_polar * i + 2 * p + 1] = intercept;
    }
  }

  /* Allocate memory for the interpolation table on the device */
  cudaMalloc((void**)&_exp_table, _exp_table_size * sizeof(FP_PRECISION));

  /* Copy exponential interpolation table to the device */
  cudaMemcpy((void*)_exp_table, (void*)exp_table,
             _exp_table_size * sizeof(FP_PRECISION),
             cudaMemcpyHostToDevice);

  /* Copy table size and spacing to constant memory on the device */
  cudaMemcpyToSymbol(exp_table_spacing, (void*)&_exp_table_spacing,
                     sizeof(FP_PRECISION), 0, cudaMemcpyHostToDevice);

  cudaMemcpyToSymbol(inverse_exp_table_spacing,
                     (void*)&_inverse_exp_table_spacing,
                     sizeof(FP_PRECISION), 0, cudaMemcpyHostToDevice);

  cudaMemcpyToSymbol(exp_table_max_index, (void*)&_exp_table_max_index,
                    sizeof(int), 0, cudaMemcpyHostToDevice);

  free(exp_table);

  return;
}


/**
 * @brief Zero each Track's boundary fluxes for each energy group and polar
 *        angle in the "forward" and "reverse" directions.
 */
void GPUSolver::zeroTrackFluxes() {
  int size = 2 * _tot_num_tracks * _num_polar * _num_groups;
  size *= sizeof(FP_PRECISION);
  cudaMemset(_boundary_flux, 0.0, size);
  return;
}


/**
 * @brief Set the FSR scalar flux for each energy group to some value.
 * @param value the value to assign to each FSR scalar flux
 */
void GPUSolver::flattenFSRFluxes(FP_PRECISION value) {

  int size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);

  cudaMemset(_scalar_flux, value, size);

  return;
}


/**
 * @brief Set the FSR source for each energy group to some value.
 * @param value the value to assign to each FSR source
 */
void GPUSolver::flattenFSRSources(FP_PRECISION value) {

  int size = _num_FSRs * sizeof(FP_PRECISION);

  cudaMemset(_old_fission_sources, value, size);

  return;
}


/**
 * @brief Normalizes all FSR scalar fluxes and Track boundary angular
 *        fluxes to the total fission source (times \f$ \nu \f$).
 */
void GPUSolver::normalizeFluxes() {

  int shared_mem = sizeof(FP_PRECISION) * _T;

  computeFissionSourcesOnDevice<<<_B, _T, shared_mem>>>(_FSR_volumes,
                                                        _FSR_materials,
                                                        _materials,
                                                        _scalar_flux,
                                                        _fission_sources);

  FP_PRECISION norm_factor = 1.0 / thrust::reduce(_fission_sources_vec.begin(),
                                                  _fission_sources_vec.end());
  normalizeFluxesOnDevice<<<_B, _T>>>(_scalar_flux, _boundary_flux,norm_factor);
}


/**
 * @brief Computes the total source (fission and scattering) in each FSR.
 * @details This method computes the total source in each FSR based on
 *          this iteration's current approximation to the scalar flux. A
 *          residual for the source with respect to the source compute on
 *          the previous iteration is computed and returned. The residual
 *          is determined as follows:
 *          /f$ res = \sqrt{\frac{\displaystyle\sum \displaystyle\sum
 *                    \left(\frac{Q^i - Q^{i-1}{Q^i}\right)^2}{\# FSRs}}} /f$
 *
 * @return the residual between this source and the previous source
 */
FP_PRECISION GPUSolver::computeFSRSources() {

  computeFSRSourcesOnDevice<<<_B, _T>>>(_FSR_materials, _materials,
                                        _scalar_flux, _old_fission_sources,
                                        _reduced_sources, 1.0 / _k_eff,
                                        _source_residuals);

  FP_PRECISION residual = thrust::reduce(_source_residuals_vec.begin(),
                                         _source_residuals_vec.end());
  residual = sqrt(residual / (_num_groups * _num_fissionable_FSRs));

  return residual;
}



/**
 * @brief This method performs one transport sweep of all azimuthal angles,
 *        Tracks, Track segments, polar angles and energy groups.
 * @details The method integrates the flux along each Track and updates the
 *          boundary fluxes for the corresponding output Track, while updating
 *          the scalar flux in each flat source region.
 */
void GPUSolver::transportSweep() {

  int shared_mem = _T * _two_times_num_polar * sizeof(FP_PRECISION);
  int tid_offset, tid_max;

  log_printf(DEBUG, "Transport sweep on device with %d blocks and %d threads",
             _B, _T);

  /* Initialize leakage to zero */
  thrust::fill(_leakage_vec.begin(), _leakage_vec.end(), 0.0);

  /* Initialize flux in each FSR to zero */
  flattenFSRFluxes(0.0);

  /* Sweep the first halfspace of azimuthal angle space */
  tid_offset = 0;
  tid_max = (_tot_num_tracks / 2);

  transportSweepOnDevice<<<_B, _T, shared_mem>>>(_scalar_flux, _boundary_flux,
                                                 _reduced_sources, _leakage,
                                                 _materials, _dev_tracks,
                                                 _exp_table,
                                                 tid_offset, tid_max);

  /* Sweep the second halfspace of azimuthal angle space */
  tid_offset = tid_max * _num_groups;
  tid_max = _tot_num_tracks;

  transportSweepOnDevice<<<_B, _T, shared_mem>>>(_scalar_flux, _boundary_flux,
                                                 _reduced_sources, _leakage,
                                                 _materials, _dev_tracks,
                                                 _exp_table,
                                                 tid_offset, tid_max);
}


/**
 * @brief Add the source term contribution in the transport equation to
 *        the FSR scalar flux.
 */
void GPUSolver::addSourceToScalarFlux() {

  addSourceToScalarFluxOnDevice<<<_B,_T>>>(_scalar_flux, _reduced_sources,
                                           _FSR_volumes, _FSR_materials,
                                           _materials);
}


/**
 * @brief Compute \f$ k_{eff} \f$ from the total, fission and scattering
 *        reaction rates and leakage.
 * @details This method computes the current approximation to the
 *          multiplication factor on this iteration as follows:
 *          \f$ k_{eff} = \frac{\displaystyle\sum_{i \in I}
 *                        \displaystyle\sum_{g \in G} \nu \Sigma^F_g \Phi V_{i}}
 *                        {\displaystyle\sum_{i \in I}
 *                        \displaystyle\sum_{g \in G} (\Sigma^T_g \Phi V_{i} -
 *                        \Sigma^S_g \Phi V_{i} - L_{i,g})} \f$
 */
void GPUSolver::computeKeff() {

  FP_PRECISION total;
  FP_PRECISION fission;
  FP_PRECISION scatter;
  FP_PRECISION leakage;

  /* Compute the total, fission and scattering reaction rates on device.
   * This kernel stores partial rates in a Thrust vector with as many
   * entries as CUDAthreads executed by the kernel */
  computeKeffReactionRates<<<_B, _T>>>(_FSR_volumes, _FSR_materials,
                                       _materials, _scalar_flux,
                                       _total, _fission, _scatter);

  cudaDeviceSynchronize();

  /* Compute the total reaction rate by reducing the partial total
   * rates compiled in the Thrust vector */
  total = thrust::reduce(_total_vec.begin(), _total_vec.end());

  /* Compute the fission rate by reducing the partial fission
   * rates compiled in the Thrust vector */
  fission = thrust::reduce(_fission_vec.begin(), _fission_vec.end());

  cudaMemcpy((void*)&fission, (void*)_fission,
             _B * _T * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);

  /* Compute the scattering rate by reducing the partial fission
   * rates compiled in the Thrust vector */
  scatter = thrust::reduce(_scatter_vec.begin(), _scatter_vec.end());

  /* Compute the total leakage by reducing the partial leakage
   * rates compiled in the Thrust vector */
  leakage = 0.5 * thrust::reduce(_leakage_vec.begin(), _leakage_vec.end());


  /* Compute the new keff from the total, fission, scatter and leakage */
  _k_eff = fission / (total - scatter + leakage);

  log_printf(DEBUG, "tot = %f, fiss = %f, scatt = %f, leak = %f,"
             " keff = %f", total, fission, scatter, leakage, _k_eff);
}


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
void GPUSolver::computeFSRFissionRates(double* fission_rates, int num_FSRs) {

  log_printf(INFO, "Computing FSR fission rates...");

  /* Allocate memory for the FSR fission rates on the device */
  double* dev_fission_rates;
  cudaMalloc((void**)&dev_fission_rates, _num_FSRs * sizeof(double));

  /* Compute the FSR fission rates on the device */
  computeFSRFissionRatesOnDevice<<<_B,_T>>>(dev_fission_rates,
                                           _FSR_materials,
                                           _materials,
                                           _scalar_flux);

  /* Copy the fission rate array from the device to the host */
  cudaMemcpy((void*)fission_rates, (void*)dev_fission_rates,
             _num_FSRs * sizeof(double), cudaMemcpyDeviceToHost);

  /* Deallocate the memory assigned to store the fission rates on the device */
  cudaFree(dev_fission_rates);

  return;
}
