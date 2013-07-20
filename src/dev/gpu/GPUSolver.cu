#include "GPUSolver.h"


/** The number of azimuthal angles */
__constant__ int num_azim[1];

/** The number of energy groups */
__constant__ int num_groups[1];

/** The number of flat source regions */
__constant__ int num_FSRs[1];

/** The number of polar angles */
__constant__ int num_polar[1];

/** Twice the number of polar angles */
__constant__ int two_times_num_polar[1];

/** The number of polar angles times energy groups */
__constant__ int polar_times_groups[1];

/** An array for the sines of the polar angle in the polar quadrature set */
__constant__ FP_PRECISION sinthetas[MAX_POLAR_ANGLES];

/** An array of the weights for the polar angles from the quadrature set */
__constant__ FP_PRECISION polar_weights[MAX_POLAR_ANGLES*MAX_AZIM_ANGLES];

/** A pointer to an array with the number of tracks per azimuthal angle */
__constant__ int num_tracks[MAX_AZIM_ANGLES/2];

/** The total number of tracks */
__constant__ int tot_num_tracks[1];

__constant__ bool interpolate_exponential[1];

/** The maximum index of the exponential prefactor array */
__constant__ int prefactor_max_index[1];

/** The spacing for the exponential prefactor array */
__constant__ FP_PRECISION prefactor_spacing[1];

/** The inverse spacing for the exponential prefactor array */
__constant__ FP_PRECISION inverse_prefactor_spacing[1];



/**
 * @brief Compute the total fission source from all flat source regions.
 * @param FSR_volumes an array of flat source region volumes
 * @param FSR_materials an array of flat source region materials
 * @param materials an array of materials on the device
 * @param scalar_flux the scalar flux in each flat source region
 * @param fission_sources array of fission sources in each flat source region
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
    double* nu_sigma_f;
    FP_PRECISION volume;
    FP_PRECISION source;

    /* Initialize fission source to zero */
    shared_fission_source[threadIdx.x] = 0;

    /* Iterate over all FSRs */
    while (tid < *num_FSRs) {

	curr_material = &materials[FSR_materials[tid]];
	nu_sigma_f = curr_material->_nu_sigma_f;
	volume = FSR_volumes[tid];

	/* Iterate over all energy groups and update
	 * fission source for this block */
	for (int e=0; e < *num_groups; e++) {
	    source = nu_sigma_f[e] * scalar_flux(tid,e) * volume;
	    shared_fission_source[threadIdx.x] += source;
	}
		 
	/* Increment thread id */
	tid += blockDim.x * gridDim.x;
    }

    /* Copy this threads fission source to global memory */
    tid = threadIdx.x + blockIdx.x * blockDim.x;
    fission_sources[tid] = shared_fission_source[threadIdx.x];
    
    return;
}


/**
 * @brief Normalizes all flatsourceregion scalar fluxes and track boundary
 *        angular fluxes to the total fission source (times \f$ \nu \f$).
 * @param scalar_flux an array of the flat source region scalar fluxes
 * @param boundary_flux an array of the boundary fluxes
 * @param norm_factor the normalization factor
 */
__global__ void normalizeFluxesOnDevice(FP_PRECISION* scalar_flux, 
					FP_PRECISION* boundary_flux, 
					FP_PRECISION norm_factor) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    
    /* Normalize scalar fluxes for each flat source region */
    while(tid < *num_FSRs) {

        for (int e=0; e < *num_groups; e++)
	    scalar_flux(tid,e) *= norm_factor;

	tid += blockDim.x * gridDim.x;
    }

    tid = threadIdx.x + blockIdx.x * blockDim.x;

    /* Normalize angular boundary fluxes for each track */
    while(tid < *tot_num_tracks) {
        for (int pe2=0; pe2 < 2*(*polar_times_groups); pe2++)
	    boundary_flux(tid,pe2) *= norm_factor;

	tid += blockDim.x * gridDim.x;
    }

    return;
}


/**
 * @brief Computes the total source (fission and scattering) in each flat 
 *        source region.
 * @details This method computes the total source in each region based on
 *          this iteration's current approximation to the scalar flux. A
 *          residual for the source with respect to the source compute on
 *          the previous iteration is computed and returned. The residual
 *          is determined as follows:
 *          /f$ res = \sqrt{\frac{\displaystyle\sum \displaystyle\sum 
 *                    \left(\frac{Q^i - Q^{i-1}{Q^i}\right)^2}{# FSRs}}} \f$
 *
 * @param FSR_materials an array of flat source region material UIDs
 * @param materials an array of material pointers
 * @param scalar_flux an array of flat source region scalar fluxes
 * @param source an array of flat source region sources
 * @param reduced_source an array of flat source region sources / total xs
 * @param inverse_k_eff the inverse of keff
 * @param an array of the source residuals 
 * @return the residual between this source and the previous source
 */
__global__ void computeFSRSourcesOnDevice(int* FSR_materials,
					  dev_material* materials,
					  FP_PRECISION* scalar_flux,
					  FP_PRECISION* source,
					  FP_PRECISION* old_source,
					  FP_PRECISION* reduced_source,
					  FP_PRECISION inverse_k_eff,
					  FP_PRECISION* source_residuals) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    /* Reset the residual for the old and new fission sources to zero */
    source_residuals[threadIdx.x + blockIdx.x * blockDim.x] = 0.0;

    FP_PRECISION fission_source;
    FP_PRECISION scatter_source;

    dev_material* curr_material;

    double* nu_sigma_f;
    double* sigma_s;
    double* sigma_t;
    double* chi;

    /* Iterate over all FSRs */
    while (tid < *num_FSRs) {

	curr_material = &materials[FSR_materials[tid]];

	nu_sigma_f = curr_material->_nu_sigma_f;
	sigma_s = curr_material->_sigma_s;
	sigma_t = curr_material->_sigma_t;
	chi = curr_material->_chi;

	/* Initialize the fission source to zero for this region */
	fission_source = 0;
	
	/* Compute total fission source for current region */
	for (int e=0; e < *num_groups; e++)
	    fission_source += scalar_flux(tid,e) * nu_sigma_f[e];
      
	/* Compute total scattering source for region for group G */
	for (int G=0; G < *num_groups; G++) {
	    scatter_source = 0;
	
	    for (int g=0; g < *num_groups; g++)
	        scatter_source += 
		    sigma_s[G*(*num_groups)+g] * scalar_flux(tid,g);
	
	    /* Set the total source for this region in this group */
	    source(tid,G) = (inverse_k_eff * fission_source * chi[G] +
			     scatter_source) * ONE_OVER_FOUR_PI;

	    reduced_source(tid,G) = __fdividef(source(tid,G), sigma_t[G]);

	    /* Compute the norm of residuals of the sources for convergence */
	    if (fabs(source(tid,G)) > 1E-10)
	        source_residuals[threadIdx.x + blockIdx.x * blockDim.x] +=
		    pow((source(tid,G) - old_source(tid,G)) /
		         source(tid,G), 2);

	    /* Update the old source */	
	    old_source(tid,G) = source(tid,G);
	}
	
	/* Increment the thread id */
	tid += blockDim.x * gridDim.x;
    }

    return;
}





/**
 * @brief Compute the total fission source from all flat source regions.
 * @param FSR_volumes an array of the flat source region volumes
 * @param FSR_materials an array of the flat source region material UIDs
 * @param materials an array of the material pointers
 * @param scalar_flux an array of flat source region scalar fluxes
 * @param tot_absorption an array of flat source region absorption rates
 * @param tot_fission an array of flat source region fission rates
 */
__global__ void computeFissionAndAbsorption(FP_PRECISION* FSR_volumes,
					    int* FSR_materials,
					    dev_material* materials,
					    FP_PRECISION* scalar_flux,
					    FP_PRECISION* tot_absorption,
					    FP_PRECISION* tot_fission) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    dev_material* curr_material;
    double* nu_sigma_f;
    double* sigma_a;
    FP_PRECISION volume;

    double absorption = 0.;
    double fission = 0.;

    /* Iterate over all FSRs */
    while (tid < *num_FSRs) {
        
	curr_material = &materials[FSR_materials[tid]];
	nu_sigma_f = curr_material->_nu_sigma_f;
	sigma_a = curr_material->_sigma_a;
	volume = FSR_volumes[tid];

	double curr_abs = 0.;
	double curr_fission = 0.;

	/* Iterate over all energy groups and update
	 * fission and absorption rates for this block */
	for (int e=0; e < *num_groups; e++) {
	    curr_abs += sigma_a[e] * scalar_flux(tid,e);
	    curr_fission += nu_sigma_f[e] * scalar_flux(tid,e);
	}

	absorption += curr_abs * volume;
	fission += curr_fission * volume;

	/* Increment thread id */
	tid += blockDim.x * gridDim.x;
    }

    /* Copy this thread's fission and absorption rates to global memory */
    tid = threadIdx.x + blockIdx.x * blockDim.x;
    tot_absorption[tid] = absorption;
    tot_fission[tid] = fission;

    return;
}


/**
 * @brief Compute the index into the exponential prefactor hashtable.
 * @details This method computes the index into the exponential prefactor
 *          hashtable for a segment length multiplied by the total 
 *          cross-section of the material the segment resides in.
 * @param sigm_t_l the cross-section multiplied by segment length
 * @return the hasthable index
 */ 
__device__ int computePrefactorIndex(FP_PRECISION tau) {
    int index = tau * *inverse_prefactor_spacing;
    index *= *two_times_num_polar;
    return index;
}


/**
 * @brief Perform an atomic addition in double precision to an array address.
 * @details This method is straight out of CUDA C Developers Guide (cc 2013)
 * @param address the array memory address
 * @param value the value to add to the array 
 * @return the atomically added array value and input value
 *
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
 *        track segment.
 * @details This method computes \f$ 1 - exp(-l\Sigma^T_g/sin(\theta_p)) \f$ 
 *          for a segment with total group cross-section and for
 *          some polar angle.
 * @brief sigma_t the total group cross-section at this energy
 * @brief length the length of the line segment projected in the xy-plane
 * @param _prefactor_array the exponential prefactor interpolation table
 * @brief p the polar angle index
 * @return the evaluated exponential
 */
__device__ FP_PRECISION computeExponential(FP_PRECISION sigma_t, 
					   FP_PRECISION length,
					   FP_PRECISION* _prefactor_array,
					   int p) {

    FP_PRECISION exponential;
    FP_PRECISION tau = sigma_t * length;

    /* Evaluate the exponential using the lookup table - linear interpolation */
    if (*interpolate_exponential) {
        int index;

	index = int(tau * (*inverse_prefactor_spacing)) * (*two_times_num_polar);
	exponential = (1. - (_prefactor_array[index+2 * p] * tau + 
			  _prefactor_array[index + 2 * p +1]));
    }

    /* Evalute the exponential using the intrinsic exp function */
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
 * @brief Computes the contribution to the flat source region scalar flux
 *        from a single track segment and a single energy group.
 * @details This method integrates the angular flux for a track segment across
 *        energy groups and polar angles, and tallies it into the flat
 *        source region scalar flux, and updates the track's angular flux.
 * @param curr_segment a pointer to the segment of interest
 * @param energy_group the energy group of interest
 * @param materials the array of materials
 * @param track_flux a pointer to the track's angular flux
 * @param reduced_source the array of flat source region sources / total xs
 * @param polar_weights the array of polar quadrature weights
 * @param _prefactor_array the exponential prefactor interpolation table
 * @param scalar_flux the array of flat source region scalar fluxes
 */
__device__ void scalarFluxTally(dev_segment* curr_segment, 
				int energy_group,
				dev_material* materials,
				FP_PRECISION* track_flux,
				FP_PRECISION* reduced_source,
				FP_PRECISION* polar_weights,
				FP_PRECISION* _prefactor_array,
				FP_PRECISION* scalar_flux) {

    int fsr_id = curr_segment->_region_uid;
    FP_PRECISION length = curr_segment->_length;
    dev_material* curr_material = &materials[curr_segment->_material_uid];
    double *sigma_t = curr_material->_sigma_t;

    /* The average flux long this segment in this flat source region */
    FP_PRECISION psibar;
    FP_PRECISION exponential;

    /* Zero the FSR scalar flux contribution from this segment 
     * and energy group */
    FP_PRECISION fsr_flux = 0.0;
    
    /* Compute the exponential prefactor hashtable index */
    
    /* Loop over polar angles */
    for (int p=0; p < *num_polar; p++) {
        exponential = computeExponential(sigma_t[energy_group], 
					 length, _prefactor_array, p);
        psibar = (track_flux[p] - reduced_source(fsr_id,energy_group)) * 
	         exponential;
	fsr_flux += psibar * polar_weights[p];
	track_flux[p] -= psibar;
    }

    /* Atomically increment the scalar flux for this flat source region */
    atomicAdd(&scalar_flux(fsr_id,energy_group), fsr_flux);
}


/**
 * @brief Updates the boundary flux for a track given boundary conditions.
 * @details For reflective boundary conditions, the outgoing boundary flux
 *          for the track is given to the reflecting track. For vacuum
 *          boundary conditions, the outgoing flux tallied as leakage.
 *          NOTE: Only one energy group is transferred by this routine.
 * @param curr_track a pointer to the track of interest
 * @param track_flux an array of the outgoing track flux
 * @param boundary_flux an array of all angular fluxes
 * @param leakage an array of leakages for each CUDA thread
 * @param polar_weights an array of polar quadrature weights
 * @param energy_angle_index the energy group index
 * @param direction the track direction (forward - true, reverse - false)
 */
__device__ void transferBoundaryFlux(dev_track* curr_track,
				     FP_PRECISION* track_flux, 
				     FP_PRECISION* boundary_flux, 
				     FP_PRECISION* leakage,
				     FP_PRECISION* polar_weights,
				     int energy_angle_index,
				     bool direction) {

    int start = energy_angle_index;
    bool bc;
    int track_out_id;

    /* Extract boundary conditions for this track and the pointer to the
     * outgoing reflective track, and index into the leakage array */

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

    /* Put track's flux in the shared memory temporary flux array */
    for (int p=0; p < *num_polar; p++) {
	track_out_flux[p] = track_flux[p] * bc;
	leakage[0] += track_flux[p] * polar_weights[p] * (!bc);
    }
}


/**
 * @brief This method performs one transport sweep of one halfspace of all 
 *        azimuthal angles, tracks, segments, polar angles and energy groups.
 * @details The method integrates the flux along each track and updates the 
 *          boundary fluxes for the corresponding output track, while updating 
 *          the scalar flux in each flat source region
 * @param scalar_flux an array of flat source region scalar fluxes
 * @param boundary_flux an array of boundary fluxes
 * @param reduced_source an array of flat source region sources / total xs
 * @param leakage an array of angular flux leakaages
 * @param materials an array of material pointers
 * @param tracks an array of tracks
 * @param _prefactor_array an array for the exponential prefactor table
 * @param tid_offset the track offset for azimuthal angle halfspace
 * @param tid_max the upper bound on the track IDs for this azimuthal 
 *                angle halfspace
 */
__global__ void transportSweepOnDevice(FP_PRECISION* scalar_flux,
				       FP_PRECISION* boundary_flux,
				       FP_PRECISION* reduced_source,
				       FP_PRECISION* leakage,
				       dev_material* materials,
				       dev_track* tracks,
				       FP_PRECISION* _prefactor_array,
				       int tid_offset,
				       int tid_max) {

    /* Shared memory buffer for each thread's angular flux */
    extern __shared__ FP_PRECISION temp_flux[];
    FP_PRECISION* track_flux;

    int tid = tid_offset + threadIdx.x + blockIdx.x * blockDim.x;

    int track_id = int(tid / *num_groups);
    int track_flux_index = threadIdx.x * (*two_times_num_polar);
    int energy_group = tid % (*num_groups);
    int energy_angle_index = energy_group * (*num_polar);

    dev_track* curr_track;
    int num_segments;
    dev_segment* curr_segment;

    /* Iterate over track with azimuthal angles in (0, pi/2) */
    while (track_id < tid_max) {

        /* Initialize local registers with important data */
        curr_track = &tracks[track_id];
      	num_segments = curr_track->_num_segments;

	/* Retrieve a pointer to this thread's shared memory buffer for angular flux */
	track_flux = &temp_flux[track_flux_index];
      
	/* Put track's flux in the shared memory temporary flux array */
      	for (int p=0; p < *num_polar; p++) {
	
	    /* Forward flux along this track */
	    track_flux[p] = boundary_flux(track_id,p+energy_angle_index);

	    /* Reverse flux along this track */
	    track_flux[(*num_polar) + p] = boundary_flux(track_id,p+energy_angle_index+(*polar_times_groups));
      	}
      
	/* Loop over each segment in forward direction */
	for (int i=0; i < num_segments; i++) {
	    curr_segment = &curr_track->_segments[i];
	    scalarFluxTally(curr_segment, energy_group, materials,
			    track_flux, reduced_source, polar_weights,
			    _prefactor_array, scalar_flux);
	}

	/* Transfer flux to outgoing track */
	transferBoundaryFlux(curr_track, track_flux, boundary_flux, 
			     &leakage[threadIdx.x + blockIdx.x * blockDim.x], 
			     polar_weights, energy_angle_index, true);

	/* Loop over each segment in reverse direction */
	track_flux = &temp_flux[track_flux_index + (*num_polar)];

	for (int i=num_segments-1; i > -1; i--) {
	    curr_segment = &curr_track->_segments[i];
	    scalarFluxTally(curr_segment, energy_group, materials,
			    track_flux, reduced_source, polar_weights,
			    _prefactor_array, scalar_flux);
	}

	/* Transfer flux to outgoing track */
	transferBoundaryFlux(curr_track, track_flux, boundary_flux, 
			     &leakage[threadIdx.x + blockIdx.x * blockDim.x], 
			     polar_weights, energy_angle_index, false);

        /* Update the indices for this thread to the next track, energy group */
	tid += blockDim.x * gridDim.x;
        track_id = int(tid / *num_groups);
	energy_group = tid % (*num_groups);
	energy_angle_index = energy_group * (*num_polar);
    }

    return;
}


/**
 * @brief Add the source term contribution in the transport equation to 
 *        the flat source region scalar flux
 * @param scalar_flux an array of flat source region scalar fluxes
 * @param reduced_source an array of flat source region sources / total xs
 * @param FSR_volumes an array of flat source region volumes
 * @param FSR_materials an array of flat source region material UIDs
 * @param materials an array of material pointers
 */
__global__ void addSourceToScalarFluxOnDevice(FP_PRECISION* scalar_flux,
					      FP_PRECISION* reduced_source,
					      FP_PRECISION* FSR_volumes,
					      int* FSR_materials,
					      dev_material* materials) {
  
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    FP_PRECISION volume;
    
    dev_material* curr_material;
    double* sigma_t;

    /* Iterate over all FSRs */
    while (tid < *num_FSRs) {

	curr_material = &materials[FSR_materials[tid]];
	volume = FSR_volumes[tid];
	sigma_t = curr_material->_sigma_t;
	
	/* Iterate over all energy groups */
	for (int i=0; i < *num_groups; i++) {
	    scalar_flux(tid,i) *= 0.5;
	    scalar_flux(tid,i) = FOUR_PI * reduced_source(tid,i) + 
	      __fdividef(scalar_flux(tid,i), (sigma_t[i] * volume));
	}

	/* Increment thread id */
	tid += blockDim.x * gridDim.x;
    }

    return;
}




/**
 * @brief Constructor initializes array pointers for tracks and materials.
 * @details The constructor retrieves the number of energy groups and flat
 *          source regions and azimuthal angles from the geometry and track
 *          generator.
 * @param geometry an optional pointer to the geometry
 * @param track_generator an optional pointer to the track generator
 */
GPUSolver::GPUSolver(Geometry* geometry, TrackGenerator* track_generator) :

    Solver(geometry, track_generator) {

    /**************************************************************************/
    /*                        Host data initialization                        */
    /**************************************************************************/

    /* The default number of threadblocks and threads per threadblock */
     _B = 64;
    _T = 64;


    /**************************************************************************/
    /*                       Device data initialization                       */
    /**************************************************************************/

    _materials = NULL;
    _dev_tracks = NULL;

    _tot_absorption = NULL;
    _tot_fission = NULL;
    _leakage = NULL;

    if (track_generator != NULL)
        setTrackGenerator(track_generator);

    if (geometry != NULL)
        setGeometry(geometry);
}



/**
 * @brief Solver destructor frees all memory on the device, including arrays
 *        for the fluxes and sources.
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

    if (_source != NULL) {
        cudaFree(_source);
	_source = NULL;
    }

    if (_old_source != NULL) {
        cudaFree(_old_source);
	_old_source = NULL;
    }

    if (_reduced_source != NULL) {
        cudaFree(_reduced_source);
	_reduced_source = NULL;
    }

    if (_FSRs_to_powers != NULL) {
        cudaFree(_FSRs_to_powers);
        _FSRs_to_powers = NULL;
    }

    if (_FSRs_to_pin_powers != NULL) {
        cudaFree(_FSRs_to_pin_powers);
    	_FSRs_to_pin_powers = NULL;
    }

    if (_fission_sources != NULL) {
        _fission_sources_vec.clear();
	_fission_sources = NULL;
    }

    if (_tot_absorption != NULL) {
        _tot_absorption_vec.clear();
	_tot_absorption = NULL;
    }

    if (_tot_fission != NULL) {
        _tot_fission_vec.clear();
	_tot_fission = NULL;
    }

    if (_source_residuals != NULL) {
        _source_residuals_vec.clear();
	_source_residuals = NULL;
    }

    if (_leakage != NULL) {
        _leakage_vec.clear();
	_leakage = NULL;
    }

    if (_prefactor_array != NULL) {
        cudaFree(_prefactor_array);
	_prefactor_array = NULL;
    }
}


/**
 * @brief Returns the scalar flux for some energy group for a flat source region
 * @param fsr_id the ID for the FSR of interest
 * @param energy_group the energy group of interest
 */
FP_PRECISION GPUSolver::getFSRScalarFlux(int fsr_id, int energy_group) {

    /* Error checking */
    if (fsr_id >= _num_FSRs)
        log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d "
		 "in enery group %d since the solver only contains FSR with "
		   "IDs greater than or equal to %d", 
		   fsr_id, energy_group, _num_FSRs-1);

    if (fsr_id < 0)
        log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d "
		  "in energy group %d since FSRs do not have negative IDs", 
		  fsr_id, energy_group);

    if (energy_group-1 >= _num_groups)
        log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d "
		   "in energy group %d since the solver only has %d energy "
		   "groups", fsr_id, energy_group, _num_groups);

    if (energy_group <= 0)
        log_printf(ERROR, "Unable to return a scalar flux for FSR id = %d "
		 "in energy group %d since energy groups are greater than 1",
		 fsr_id, energy_group);

    /* Copy the scalar flux for this FSR and energy group from the device */
    FP_PRECISION fsr_scalar_flux;
    int flux_index = fsr_id * _num_groups + energy_group - 1;
    cudaMemcpy((void*)&fsr_scalar_flux, (void*)&_scalar_flux[flux_index], 
	       sizeof(FP_PRECISION), cudaMemcpyDeviceToHost);

    return fsr_scalar_flux;
}


/**
 * @brief Return an array indexed by flat source region IDs and energy groups 
 *        which contains the corresponding fluxes for each flat source region.
 * @return an array of flat source region scalar fluxes
 */
FP_PRECISION* GPUSolver::getFSRScalarFluxes() {

    if (_scalar_flux == NULL)
        log_printf(ERROR, "Unable to returns the device solver's scalar flux "
		   "array since it has not yet been allocated in memory");

    /* Copy the scalar flux for all FSRs from the device */
    FP_PRECISION* fsr_scalar_fluxes = new FP_PRECISION[_num_FSRs * _num_groups];
    cudaMemcpy((void*)fsr_scalar_fluxes, (void*)_scalar_flux,
	       _num_FSRs * _num_groups * sizeof(FP_PRECISION),
	       cudaMemcpyDeviceToHost);

    return fsr_scalar_fluxes;
}


/**
 * @brief Return an array indexed by flat source region IDs with the
 *        corresponding flat source region power.
 * @return an array of flat source region powers
 */
FP_PRECISION* GPUSolver::getFSRPowers() {
    if (_FSRs_to_powers == NULL)
        log_printf(ERROR, "Unable to returns the device solver's FSR power "
		   "array since it has not yet been allocated in memory");

    return _FSRs_to_powers;
}


/**
 * @brief Return an array indexed by flat source region IDs with the
 *        corresponding pin cell power.
 * @return an array of flat source region pin powers
 */
FP_PRECISION* GPUSolver::getFSRPinPowers() {
    if (_FSRs_to_pin_powers == NULL)
        log_printf(ERROR, "Unable to returns the device solver's FSR pin power "
		   "array since it has not yet been allocated in memory");

    return _FSRs_to_pin_powers;
}


/**
 * @brief Sets the number of threadblocks (>0) for device kernels
 * @param num_blocks the number of threadblocks
 */
void GPUSolver::setNumThreadBlocks(int num_blocks) {

    if (num_blocks < 0)
        log_printf(ERROR, "Unable to set the number of threadblocks to %d since "
		   "it is a negative number", num_blocks);

    _B = num_blocks;
}


/**
 * @brief Sets the number of threads per block (>0) for device kernels
 * @param num_threads the number of threads per block
 */
void GPUSolver::setNumThreadsPerBlock(int num_threads) {

    if (num_threads < 0)
        log_printf(ERROR, "Unable to set the number of threads per block to %d "
		   "since it is a negative number", num_threads);

    _T = num_threads;
}


/**
 * @brief Sets the geometry for the solver.
 * @details The geometry must already have initialized flat source region maps
 *          and segmentized the trackgenerator's tracks.
 * @param geometry a pointer to a geometry
 */
void GPUSolver::setGeometry(Geometry* geometry) {
    Solver::setGeometry(geometry);
    initializeMaterials();

    /* Copy the number of energy groups to constant memory on the GPU */
    cudaMemcpyToSymbol(num_groups, (void*)&_num_groups, sizeof(int), 0,
		       cudaMemcpyHostToDevice);
}


/**
 * @brief Sets the track generator with characteristic tracks for the solver.
 * @details The track generator must already have generated tracks and have
 *          segmentized them across the geometry.
 * @param track_generator a pointer to a trackgenerator
 */
void GPUSolver::setTrackGenerator(TrackGenerator* track_generator) {
    Solver::setTrackGenerator(track_generator);
    initializeTracks();
}


/**
 * @brief Creates a polar quadrature object for the solver.
 */
void GPUSolver::initializePolarQuadrature() {

    log_printf(INFO, "Initializing polar quadrature on the GPU...");

    /* Deletes the old quadrature if one existed */
    if (_quad != NULL)
        delete _quad;

    _quad = new Quadrature(_quadrature_type, _num_polar);
    _polar_times_groups = _num_groups * _num_polar;

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

    /* Compute polar times azimuthal angle weights */
    if (_polar_weights != NULL)
        delete [] _polar_weights;

    _polar_weights =
        (FP_PRECISION*)malloc(_num_polar * _num_azim * sizeof(FP_PRECISION));

    FP_PRECISION* multiples = _quad->getMultiples();
    double* azim_weights = _track_generator->getAzimWeights();

    for (int i=0; i < _num_azim; i++) {
        for (int j=0; j < _num_polar; j++)
	    _polar_weights[i*_num_polar+j] = azim_weights[i]*multiples[j]*FOUR_PI;
    }

    /* Copy the polar weights to constant memory on the GPU */
    cudaMemcpyToSymbol(polar_weights, (void*)_polar_weights,
		       _num_polar * _num_azim * sizeof(FP_PRECISION),
		       0, cudaMemcpyHostToDevice);
}


/**
 * @brief Allocates memory for flat source region power arrays.
 * @details Deletes memory for power arrays if they were allocated from
 *          previous simulation.
 */
void GPUSolver::initializePowerArrays() {

     log_printf(INFO, "Initializing FSR power arrays on the GPU...");

    /* Delete old power arrays if they exist */
    if (_FSRs_to_powers != NULL)
        delete [] _FSRs_to_powers;

    if (_FSRs_to_pin_powers != NULL)
        delete [] _FSRs_to_pin_powers;

    /* Allocate memory for FSR power and pin power arrays */
    try{
	_FSRs_to_powers = new FP_PRECISION[_num_FSRs];
	_FSRs_to_pin_powers = new FP_PRECISION[_num_FSRs];
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the device solver's "
		   "FSR power arrays. Backtrace:%s", e.what());
    }
}


/**
 * @brief Initializes each of the flat source region objects inside the solver's
 *        array of flatsourceregions. 
 * @details This method assigns each flat source region a unique, monotonically
 *          increasing ID, sets the material for each flat source region, and 
 *          assigns a volume based on the cumulative length of all of the 
 *          segments inside the flat source region.
 */
void GPUSolver::initializeFSRs() {

    log_printf(INFO, "Initializing FSRs on the GPU...");

    /* Delete old FSRs array if it exists */
    if (_FSR_volumes != NULL)
        cudaFree(_FSR_volumes);

    if (_FSR_materials != NULL)
        cudaFree(_FSR_materials);

    /* Allocate memory for all tracks and track offset indices on the device */
    try{

        /* Allocate memory on device for FSR volumes and material uids */
        cudaMalloc((void**)&_FSR_volumes, _num_FSRs * sizeof(FP_PRECISION));
        cudaMalloc((void**)&_FSR_materials, _num_FSRs * sizeof(int));
	
	/* Create a temporary FSR array to populate and then copy to device */
	FP_PRECISION* temp_FSR_volumes = new FP_PRECISION[_num_FSRs];

	/* Get the array indexed by FSR IDs with material ID values */
	int* FSRs_to_materials = _geometry->getFSRtoMaterialMap();

	/* Initialize each FSRs volume to 0 to avoid NaNs */
	memset(temp_FSR_volumes, FP_PRECISION(0.), 
	       _num_FSRs*sizeof(FP_PRECISION));

	Track* track;
	int num_segments;
	segment* curr_segment;
	segment* segments;
	FP_PRECISION volume;

	double* azim_weights = _track_generator->getAzimWeights();

	/* Set each FSR's volume by accumulating the total length of all
	   tracks inside the FSR. Iterate over azimuthal angle, track, segment*/
	for (int i=0; i < _num_azim; i++) {
	    for (int j=0; j < _num_tracks[i]; j++) {

	        track = &_track_generator->getTracks()[i][j];
		num_segments = track->getNumSegments();
		segments = track->getSegments();

		/* Iterate over the track's segments to update FSR volumes */
		for (int s = 0; s < num_segments; s++) {
		    curr_segment = &segments[s];
		    volume = curr_segment->_length * azim_weights[i];
		    temp_FSR_volumes[curr_segment->_region_id] += volume;
		}
	    }
	}

	/* Copy the temporary array of FSRs to the device */
	cudaMemcpy((void*)_FSR_volumes, (void*)temp_FSR_volumes,
		   _num_FSRs * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);
	cudaMemcpy((void*)_FSR_materials, (void*)FSRs_to_materials,
		   _num_FSRs * sizeof(int), cudaMemcpyHostToDevice);

	/* Copy the number of flat source regions into constant memory on 
	 * the GPU */
	cudaMemcpyToSymbol(num_FSRs, (void*)&_num_FSRs, sizeof(int), 0,
		       cudaMemcpyHostToDevice);

	/* Free the temporary array of FSRs on the host */
	free(temp_FSR_volumes);
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's flat "
		   "source regions on the device. Backtrace:%s", e.what());
    }

    initializeThrustVectors();
}


/**
 * @brief Allocates data on the GPU for all materials data.
 */
void GPUSolver::initializeMaterials() {

    log_printf(INFO, "Initializing materials on the GPU...");

    /* Delete old materials array if it exists */
    if (_materials != NULL)
        cudaFree(_materials);

    /* Allocate memory for all tracks and track offset indices on the device */
    try{

	std::map<short int, Material*> host_materials=_geometry->getMaterials();
	std::map<short int, Material*>::iterator iter;

        /* Iterate through all materials and clone them on the device */
        cudaMalloc((void**)&_materials, _num_materials * sizeof(dev_material));
	for (iter=host_materials.begin(); iter != host_materials.end(); ++iter)
	    cloneMaterialOnGPU(iter->second, &_materials[iter->second->getUid()]);
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the device solver's "
		   "materials. Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory on the GPU for all tracks in the simulation.
 */
void GPUSolver::initializeTracks() {

    log_printf(INFO, "Initializing tracks on the GPU...");

    /* Delete old tracks array if it exists */
    if (_dev_tracks != NULL)
        cudaFree(_dev_tracks);

    /* Allocate memory for all tracks and track offset indices on the device */
    try{

        /* Allocate array of tracks */
    	cudaMalloc((void**)&_dev_tracks, _tot_num_tracks * sizeof(dev_track));

        /* Iterate through all tracks and clone them on the device */
	int index;

	for (int i=0; i < _tot_num_tracks; i++) {

	    cloneTrackOnGPU(_tracks[i], &_dev_tracks[i]);

	    /* Make track reflective */
	    index = computeScalarTrackIndex(_tracks[i]->getTrackInI(),
		        		       _tracks[i]->getTrackInJ());
	    cudaMemcpy((void*)&_dev_tracks[i]._track_in,
		   (void*)&index, sizeof(int), cudaMemcpyHostToDevice);

	    index = computeScalarTrackIndex(_tracks[i]->getTrackOutI(),
		        		       _tracks[i]->getTrackOutJ());
	    cudaMemcpy((void*)&_dev_tracks[i]._track_out,
		   (void*)&index, sizeof(int), cudaMemcpyHostToDevice);
	}

	/* Copy the array of number of tracks for each azimuthal angles into 
	 * constant memory on GPU */
	cudaMemcpyToSymbol(num_tracks, (void*)_num_tracks, 
			   _num_azim * sizeof(int), 0, cudaMemcpyHostToDevice);
    
	/* Copy the total number of tracks into constant memory on GPU */
	cudaMemcpyToSymbol(tot_num_tracks, (void*)&_tot_num_tracks,
			   sizeof(int), 0, cudaMemcpyHostToDevice);

	/* Copy the number of azimuthal angles into constant memory on GPU */
	cudaMemcpyToSymbol(num_azim, (void*)&_num_azim, sizeof(int), 0, 
			   cudaMemcpyHostToDevice);
	
	/* Copy the array of number of tracks for each azimuthal angles into 
	 * constant memory on GPU */
	cudaMemcpyToSymbol(num_tracks, (void*)_num_tracks, 
			   _num_azim * sizeof(int), 0, cudaMemcpyHostToDevice);
	
	/* Copy the total number of tracks into constant memory on GPU */
	cudaMemcpyToSymbol(tot_num_tracks, (void*)&_tot_num_tracks,
			   sizeof(int), 0, cudaMemcpyHostToDevice);
    }

    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's tracks "
		   "on the device. Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory for track boundary angular fluxes and leakages
 *        flat source region scalar fluxes.
 * @details Deletes memory for old flux arrays if they were allocated from
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
        log_printf(ERROR, "Could not allocate memory for the solver's fluxes "
		   "on the device. Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory for flat source region source arrays.
 * @details Deletes memory for old source arrays if they were allocated from
 *          previous simulation.
 */
void GPUSolver::initializeSourceArrays() {

    log_printf(INFO, "Initializing source arrays on the GPU...");

    /* Delete old sources arrays if they exist */
    if (_source != NULL)
        cudaFree(_source);

    if (_old_source != NULL)
        cudaFree(_old_source);

    if (_reduced_source != NULL)
        cudaFree(_reduced_source);

    /* Allocate memory for all source arrays on the device */
    try{

        cudaMalloc((void**)&_source, 
		   _num_FSRs * _num_groups * sizeof(FP_PRECISION));

	cudaMalloc((void**)&_old_source,
		   _num_FSRs * _num_groups * sizeof(FP_PRECISION));

	cudaMalloc((void**)&_reduced_source,
		   _num_FSRs * _num_groups * sizeof(FP_PRECISION));
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's flat "
		   "source region sources array on the device. "
		   "Backtrace:%s", e.what());
    }
}


/**
 * @brief Initialize Thrust vectors for the fission and absorption rates,
 *        source residuals, leakage and fission sources.
 */
void GPUSolver::initializeThrustVectors() {

    log_printf(INFO, "Initializing thrust vectors on the GPU...");

    /* Delete old vectors if they exist */
    if (_fission_sources != NULL) {
        _fission_sources = NULL;
        _fission_sources_vec.clear();
    }

    if (_tot_absorption != NULL) {
        _tot_absorption = NULL;
        _tot_absorption_vec.clear();
    }

    if (_tot_fission != NULL) {
        _tot_fission = NULL;
        _tot_fission_vec.clear();
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
      
	/* Allocate total absorption reaction rate array on device */
	_tot_absorption_vec.resize(_B * _T);
	_tot_absorption = thrust::raw_pointer_cast(&_tot_absorption_vec[0]);

	/* Allocate fission reaction rate array on device */
	_tot_fission_vec.resize(_B * _T);
	_tot_fission = thrust::raw_pointer_cast(&_tot_fission_vec[0]);

	/* Allocate source residual array on device */
	_source_residuals_vec.resize(_B * _T);
	_source_residuals = thrust::raw_pointer_cast(&_source_residuals_vec[0]);

	/* Allocate leakage array on device */
	_leakage_vec.resize(_B * _T);
	_leakage = thrust::raw_pointer_cast(&_leakage_vec[0]);
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's "
		   "Thrust vectors.Backtrace:%s", e.what());
    }
}


/**
 * @brief This method computes the index for the jth track at azimuthal angle i.
 * @details This method is necessary since the array of tracks on the device 
 *          is a 1D array which needs a one-to-one mapping from the 2D jagged 
 *          array of tracks on the host.
 * @param i azimuthal angle number
 * @param j the jth track at angle i
 * @return an index into the device track array
 */
int GPUSolver::computeScalarTrackIndex(int i, int j) {

    int index =0;
    int p = 0;

    /* Iterate over each azimuthal angle and increment index by the number of
       tracks at each angle */
    while (p < i) {
        index += _num_tracks[p];
	p++;
    }

    /* Update index for this track since it is the jth track at angle i */
    index += j;
    
    return index;
}


/**
 * @brief Pre-computes exponential pre-factors for each segment of each track 
 *        for each polar angle. 
 * @details This method will generate a hashmap which contains values of the 
 *          prefactor for specific segment lengths (the keys into the hashmap).
 */
void GPUSolver::precomputePrefactors(){

    log_printf(INFO, "Building exponential prefactor hashtable on device...");

    /* Copy a boolean indicating whether or not to use the linear interpolation
     * table or the exp intrinsic function */
    cudaMemcpyToSymbol(interpolate_exponential, 
		       (void*)&_interpolate_exponential, 
		       sizeof(bool), 0, cudaMemcpyHostToDevice);

    /* Copy the sines of the polar angles which is needed if the user
     * requested the use of the exp intrinsic to evaluate exponentials */
    cudaMemcpyToSymbol(sinthetas, (void*)_quad->getSinThetas(),
		       _num_polar * sizeof(FP_PRECISION), 0, 
		       cudaMemcpyHostToDevice);

    /* Set size of prefactor array */
    int num_array_values = 10 * sqrt(1. / (8. * _source_convergence_thresh));
    _prefactor_spacing = 10. / num_array_values;
    _inverse_prefactor_spacing = 1.0 / _prefactor_spacing;
    _prefactor_array_size = _two_times_num_polar * num_array_values;
    _prefactor_max_index = _prefactor_array_size - _two_times_num_polar - 1;
    
    /* allocate arrays */
    FP_PRECISION* prefactor_array = new FP_PRECISION[_prefactor_array_size];
    
    FP_PRECISION expon;
    FP_PRECISION intercept;
    FP_PRECISION slope;

    /* Create prefactor array */
    for (int i = 0; i < num_array_values; i ++){
        for (int p = 0; p < _num_polar; p++){
	    expon = exp(- (i * _prefactor_spacing) / _quad->getSinTheta(p));
	    slope = - expon / _quad->getSinTheta(p);
	    intercept = expon * (1 + (i * _prefactor_spacing) /
				 _quad->getSinTheta(p));
	    prefactor_array[_two_times_num_polar * i + 2 * p] = slope;
	    prefactor_array[_two_times_num_polar * i + 2 * p + 1] = intercept;
	}
    }

    /* Allocate memory for the prefactor array on the device */
    cudaMalloc((void**)&_prefactor_array, 
	       _prefactor_array_size * sizeof(FP_PRECISION));

    /* Copy prefactor array to the device */
    cudaMemcpy((void*)_prefactor_array, (void*)prefactor_array, 
	       _prefactor_array_size * sizeof(FP_PRECISION),
	       cudaMemcpyHostToDevice);

    /* Copy prefactor array size and spacing to constant memory on the device */
    cudaMemcpyToSymbol(prefactor_spacing, (void*)&_prefactor_spacing, 
		       sizeof(FP_PRECISION), 0, cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(inverse_prefactor_spacing, 
		       (void*)&_inverse_prefactor_spacing, 
		       sizeof(FP_PRECISION), 0, cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(prefactor_max_index, (void*)&_prefactor_max_index,
		       sizeof(int), 0, cudaMemcpyHostToDevice);

    free(prefactor_array);

    return;
}


/**
 * @brief Zero each track's boundary fluxes for each energy group and polar
 *        angle in the "forward" and "reverse" directions.
 */
void GPUSolver::zeroTrackFluxes() {
    int size = 2 * _tot_num_tracks * _num_polar * _num_groups;
    size *= sizeof(FP_PRECISION);
    cudaMemset(_boundary_flux, 0.0, size);
    return;
}


/**
 * @brief Set the scalar flux for each energy group inside each flat source 
 *        region to a constant value.
 * @param value the value to assign to each flat source region flux
 */
void GPUSolver::flattenFSRFluxes(FP_PRECISION value) {
    int size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);
    cudaMemset(_scalar_flux, value, size);
    return;
}


/**
 * @brief Set the source for each energy group inside each flat source region
 *        to a constant value.
 * @param value the value to assign to each flat source region source
 */
void GPUSolver::flattenFSRSources(FP_PRECISION value) {
    int size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);
    cudaMemset(_source, value, size);
    cudaMemset(_old_source, value, size);
    return;
}


/**
 * @brief Normalizes all flat source region scalar fluxes and track boundary
 *        angular fluxes to the total fission source (times \f$ \nu \f$).
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

    normalizeFluxesOnDevice<<<_B, _T>>>(_scalar_flux, _boundary_flux, 
					norm_factor);
}


FP_PRECISION GPUSolver::computeFSRSources() {
    
    computeFSRSourcesOnDevice<<<_B, _T>>>(_FSR_materials, _materials, 
					  _scalar_flux, _source, 
					  _old_source, _reduced_source,
					  1.0 / _k_eff, _source_residuals);

    FP_PRECISION residual = thrust::reduce(_source_residuals_vec.begin(), 
					   _source_residuals_vec.end());
    residual = sqrt(residual / _num_FSRs);

    return residual;
}



void GPUSolver::transportSweep() {

    int shared_mem = _T * _two_times_num_polar * sizeof(FP_PRECISION);
    int tid_offset, tid_max;

    log_printf(DEBUG, "Transport sweep on device with %d blocks" 
                      " and %d threads", _B, _T);

    /* Initialize leakage to zero */
    thrust::fill(_leakage_vec.begin(), _leakage_vec.end(), 0.0);
    
    /* Initialize flux in each region to zero */
    flattenFSRFluxes(0.0);
    
    /* Sweep the first halfspace of azimuthal angle space */
    tid_offset = 0;
    tid_max = (_tot_num_tracks / 2);

    transportSweepOnDevice<<<_B, _T, shared_mem>>>(_scalar_flux, 
						   _boundary_flux,
						   _reduced_source, _leakage,
						   _materials, _dev_tracks,
						   _prefactor_array, 
						   tid_offset, tid_max);

    /* Sweep the second halfspace of azimuthal angle space */
    tid_offset = tid_max * _num_groups;
    tid_max = _tot_num_tracks;

    transportSweepOnDevice<<<_B, _T, shared_mem>>>(_scalar_flux,
						   _boundary_flux,
						   _reduced_source, _leakage,
						   _materials, _dev_tracks,
						   _prefactor_array,
						   tid_offset, tid_max);
}


/**
 * @brief Add the source term contribution in the transport equation to 
 *        the flat source region scalar flux
 */
void GPUSolver::addSourceToScalarFlux() {

    addSourceToScalarFluxOnDevice<<<_B,_T>>>(_scalar_flux, _reduced_source,
					     _FSR_volumes, _FSR_materials,
					     _materials);
}


/**
 * @brief Compute \f$ k_{eff} \f$ from the total fission and absorption rates.
 * @details This method computes the current approximation to the 
 *          multiplication factor on this iteration as follows:
 *          \f$ k_{eff} = \frac{\displaystyle\sum \displaystyle\sum \nu
 *                        \Sigma_f \Phi V}{\displaystyle\sum 
 *                        \displaystyle\sum \Sigma_a \Phi V} \f$
 */
void GPUSolver::computeKeff() {

    FP_PRECISION tot_absorption;
    FP_PRECISION tot_fission;
    FP_PRECISION tot_leakage;

    /* Compute the total fission and absorption rates on the device.
     * This kernel stores partial rates in a Thrust vector with as many
     * entries as GPU threads executed by the kernel */
    computeFissionAndAbsorption<<<_B, _T>>>(_FSR_volumes, _FSR_materials,
					    _materials, _scalar_flux,
					    _tot_absorption, _tot_fission);

    cudaDeviceSynchronize();

    /* Compute the total absorption rate by reducing the partial absorption
     * rates compiled in the Thrust vector */
    tot_absorption = thrust::reduce(_tot_absorption_vec.begin(),
				    _tot_absorption_vec.end());

    /* Compute the total fission rate by reducing the partial fission
     * rates compiled in the Thrust vector */
    tot_fission = thrust::reduce(_tot_fission_vec.begin(),
				 _tot_fission_vec.end());

    cudaMemcpy((void*)&tot_fission, (void*)_tot_fission, 
	       _B * _T * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);


    /* Compute the total leakage by reducing the partial leakage
     * rates compiled in the Thrust vector */
    tot_leakage = 0.5 * thrust::reduce(_leakage_vec.begin(),
				       _leakage_vec.end());


    /* Compute the new keff from the fission and absorption rates */
    _k_eff = tot_fission / (tot_absorption + tot_leakage);

    log_printf(DEBUG, "abs = %f, fiss = %f, leak = %f, keff = %f", 
	       tot_absorption, tot_fission, tot_leakage, _k_eff);
}


void GPUSolver::computePinPowers() {
    log_printf(ERROR, "Pin power computation on the GPU is not implemented!");
}
