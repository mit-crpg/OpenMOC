#include "DeviceSolver.h"


/** The number of azimuthal angles */
__constant__ int _num_azim_devc[1];

/** The number of energy groups */
__constant__ int _num_groups_devc[1];

/** The number of flat source regions */
__constant__ int _num_FSRs_devc[1];

/** The number of polar angles */
__constant__ int _num_polar_devc[1];

/** Twice the number of polar angles */
__constant__ int _two_times_num_polar_devc[1];

/** The number of polar angles times energy groups */
__constant__ int _polar_times_groups_devc[1];

/** An array of the weights for the polar angles from the quadrature set */
__constant__ FP_PRECISION _polar_weights_devc[MAX_POLAR_ANGLES*MAX_AZIM_ANGLES];

/** A pointer to an array with the number of tracks per azimuthal angle */
__constant__ int _num_tracks_devc[MAX_AZIM_ANGLES/2];

/** The total number of tracks */
__constant__ int _tot_num_tracks_devc[1];

/** An array of the cumulative number of tracks for each azimuthal angle */
__constant__ int _track_index_offsets_devc[MAX_AZIM_ANGLES/2];

/** The maximum index of the exponential prefactor array */
__constant__ int _prefactor_max_index_devc[1];

/** The spacing for the exponential prefactor array */
__constant__ FP_PRECISION _prefactor_spacing_devc[1];

/** The inverse spacing for the exponential prefactor array */
__constant__ FP_PRECISION _inverse_prefactor_spacing_devc[1];


/**
 * @brief Set the scalar flux for each energy group inside each 
 *        dev_flatsourceregion to a constant value.
 * @param value the value to assign to each flat source region flux
 */
__global__ void flattenFSRFluxes(FP_PRECISION* scalar_flux, 
					       FP_PRECISION* old_scalar_flux,
					       FP_PRECISION value) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    /* Loop over all FSRs and energy groups */
    while (tid < *_num_FSRs_devc) {
        for (int e=0; e < *_num_groups_devc; e++) {
            scalar_flux(tid,e) = value;
  	    old_scalar_flux(tid,e) = value;
         }

	tid += blockDim.x * gridDim.x;
     }

    return;
}


/**
 * @brief Set the source for each energy group inside each dev_flatsourceregion
 *        to a constant value.
 * @param value the value to assign to each flat source region source
 */
__global__ void flattenFSRSources(FP_PRECISION* source, FP_PRECISION* old_source,
				  FP_PRECISION value) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    while (tid < *_num_FSRs_devc) {
        for (int e=0; e < *_num_groups_devc; e++) {
	    source(tid,e) = value;
	    old_source(tid,e) = value;
	}

	tid += blockDim.x * gridDim.x;
    }

    return;
}


/**
 * @brief Zero each track's boundary fluxes for each energy group and polar
 *        angle in the "forward" and "reverse" directions.
 * @param boundary_flux array of angular fluxes for each track and energy group
 */
__global__ void zeroTrackFluxes(FP_PRECISION* boundary_flux) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    /* Loop over all tracks and energy groups and set each track's 
     * incoming and outgoing flux to zero */
    while(tid < *_tot_num_tracks_devc) {
        for (int pe2=0; pe2 < 2*(*_polar_times_groups_devc); pe2++)
    	    boundary_flux(tid,pe2) = 0.0;

	tid += blockDim.x * gridDim.x;
    }

    return;
}


/**
* Compute the total fission source from all flat source regions
* @param FSRs pointer to the flat source region array on the device
* @param num_FSRs pointer to an int of the number of flat source regions
* @param materials pointer an array of materials on the device
* @param fission_source pointer to the value for the total fission source
*/
__global__ void computeFissionSource(dev_flatsourceregion* FSRs,
				     dev_material* materials,
				     FP_PRECISION* scalar_flux,
				     FP_PRECISION* fission_source) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    extern __shared__ FP_PRECISION shared_fission_source[];
    dev_flatsourceregion* curr_FSR;
    dev_material* curr_material;
    double* nu_sigma_f;
    FP_PRECISION volume;

    /* Initialize fission source to zero */
    shared_fission_source[threadIdx.x] = 0;

    /* Iterate over all FSRs */
    while (tid < *_num_FSRs_devc) {

        curr_FSR = &FSRs[tid];
	curr_material = &materials[curr_FSR->_material_uid];
	nu_sigma_f = curr_material->_nu_sigma_f;
	volume = curr_FSR->_volume;

	/* Iterate over all energy groups and update
	 * fission source for this block */
	for (int e=0; e < *_num_groups_devc; e++)
	    shared_fission_source[threadIdx.x] += nu_sigma_f[e] * scalar_flux(tid,e) * volume;

	/* Increment thread id */
	tid += blockDim.x * gridDim.x;
    }

    /* Copy this threads fission source to global memory */
    tid = threadIdx.x + blockIdx.x * blockDim.x;
    fission_source[tid] = shared_fission_source[threadIdx.x];
    
    return;
}


/**
 * @brief Normalizes all flatsourceregion scalar fluxes and track boundary
 *        angular fluxes to the total fission source (times nu).
 */
__global__ void normalizeFluxes(FP_PRECISION* scalar_flux, 
				FP_PRECISION* boundary_flux, 
				FP_PRECISION norm_factor) {

    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    
    /* Normalize scalar fluxes for each flat source region */
    
    while(tid < *_num_FSRs_devc) {
        for (int e=0; e < *_num_groups_devc; e++)
	  scalar_flux(tid,e) *= norm_factor;

	tid += blockDim.x * gridDim.x;
    }

    tid = threadIdx.x + blockIdx.x * blockDim.x;

    /* Normalize angular boundary fluxes for each track */
    while(tid < *_tot_num_tracks_devc) {
        for (int pe2=0; pe2 < 2*(*_polar_times_groups_devc); pe2++)
	    boundary_flux(tid,pe2) *= norm_factor;

	tid += blockDim.x * gridDim.x;
    }

    return;
}


/**
 * DeviceSolver constructor
 * @param geom pointer to the geometry
 * @param track_generator pointer to the TrackGenerator on the CPU
 */
DeviceSolver::DeviceSolver(Geometry* geometry, TrackGenerator* track_generator) {

    /**************************************************************************/
    /*                        Host data initialization                        */
    /**************************************************************************/

    /* The default number of threadblocks and threads per threadblock */
    _num_blocks = 64;
    _num_threads = 64;

    if (geometry != NULL)
        setGeometry(geometry);
    else
        _geometry = NULL;

    if (track_generator != NULL)
        setTrackGenerator(track_generator);
    else {
        _track_generator = NULL;
	_host_tracks = NULL;
	_num_tracks = NULL;
    }

    /* Default polar quadrature */
    _quad = NULL;
    _quadrature_type = TABUCHI;
    _num_polar = 3;
    _two_times_num_polar = 2 * _num_polar;

    _leakage = 0.0;
    _num_iterations = 0;
    _converged_source = false;
    _source_convergence_thresh = 1E-3;
    _flux_convergence_thresh = 1E-5;


    /**************************************************************************/
    /*                       Device data initialization                       */
    /**************************************************************************/

    _FSRs = NULL;
    _materials = NULL;
    _dev_tracks = NULL;
    _track_index_offsets = NULL;

    _boundary_flux = NULL;
    _scalar_flux = NULL;
    _old_scalar_flux = NULL;
    _source = NULL;
    _old_source = NULL;
    _ratios = NULL;

    _fission_source = NULL;
    _tot_abs = NULL;
    _tot_fission = NULL;
    _source_residual = NULL;

    _FSRs_to_powers = NULL;
    _FSRs_to_pin_powers = NULL;

    _prefactor_array = NULL;
    _k_eff_pinned = NULL;
    _norm_factor_pinned = NULL;
}



/**
 * Solver destructor frees all memory on the device
 */
//TODO: Fix this!!!
DeviceSolver::~DeviceSolver() {

    log_printf(NORMAL, "Cleaning up memory on the device...");

    /* Free FSRs, materials and tracks on device */
    if (_FSRs != NULL)

        cudaFree(_FSRs);
    if (_materials != NULL)
        cudaFree(_materials);

    if (_dev_tracks != NULL)
        cudaFree(_dev_tracks);

    if (_track_index_offsets != NULL)
        cudaFree(_track_index_offsets);

    if (_boundary_flux != NULL)
        cudaFree(_boundary_flux);

    if (_scalar_flux != NULL)
        cudaFree(_scalar_flux);

    if (_old_scalar_flux != NULL)
        cudaFree(_old_scalar_flux);

    if (_source != NULL)
        cudaFree(_source);

    if (_old_source != NULL)
        cudaFree(_old_source);

    if (_ratios != NULL)
        cudaFree(_ratios);

    if (_FSRs_to_powers != NULL)
        cudaFree(_FSRs_to_powers);

    if (_FSRs_to_pin_powers != NULL)
        cudaFree(_FSRs_to_pin_powers);

    if (_fission_source != NULL)
        _fission_source_vec.clear();

    if (_tot_abs != NULL)
        _tot_abs_vec.clear();

    if (_tot_fission != NULL)
        _tot_fission_vec.clear();

    if (_source_residual != NULL)
        _source_residual_vec.clear();

    if (_prefactor_array != NULL)
        cudaFree(_prefactor_array);

    if (_k_eff_pinned != NULL)
        cudaFreeHost(_k_eff_pinned);

    if (_norm_factor_pinned != NULL)
        cudaFreeHost(_norm_factor_pinned);
}


/**
 * @brief Returns a pointer to the geometry for this solver.
 * @return a pointer to the geometry
 */
Geometry* DeviceSolver::getGeometry() {

    if (_geometry == NULL)
        log_printf(ERROR, "Unable to return the device solver's geometry since "
		   "it has not yet been set");

    return _geometry;
}


/**
 * @brief Returns a pointer to the geometry for this solver.
 * @return a pointer to the geometry
 */
TrackGenerator* DeviceSolver::getTrackGenerator() {

    if (_track_generator == NULL)
        log_printf(ERROR, "Unable to return the device solver's track "
		   "genetrator since it has not yet been set");

    return _track_generator;
}


/**
 * @brief Returns the number of angles used for the polar quadrature.
 * @return the number of polar angles
 */
int DeviceSolver::getNumPolarAngles() {
    return _num_polar;
}


/**
 * @brief Returns the type of polar quadrature in use (TABUCHI or LEONARD).
 * @return the type of polar quadrature
 */
quadratureType DeviceSolver::getPolarQuadratureType() {
    return _quadrature_type;
}


/**
 * @brief Returns the number of transport sweeps to converge the source.
 * @return the number of iterations
 */
int DeviceSolver::getNumIterations() {
    return _num_iterations;
}


/**
 * @brief Returns the threshold for source convergence.
 * @return the threshold for source convergence
 */
FP_PRECISION DeviceSolver::getSourceConvergenceThreshold() {
    return _source_convergence_thresh;
}


/**
 * @brief Returns the threshold for flux convergence in fixed source iteration
 *        after the source has converged.
 * @return the threshold for flux convergence
 */
FP_PRECISION DeviceSolver::getFluxConvergenceThreshold() {
    return _flux_convergence_thresh;
}


/**
 * @brief
 * @details
 */
FP_PRECISION DeviceSolver::getFSRScalarFlux(int fsr_id, int energy_group) {

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
 * @brief Return a 2D array indexed by flatsourceregion IDs and energy groups 
 *        which contains the corresponding fluxes for each flatsourceregion.
 * @return a 2D array of dev_flatsourceregion scalar fluxes
 */
FP_PRECISION* DeviceSolver::getFSRScalarFluxes() {

    if (_scalar_flux == NULL)
        log_printf(ERROR, "Unable to returns the device solver's scalar flux "
		   "array since it has not yet been allocated in memory");

    /* Copy the scalar flux for all FSRs from the device */
    FP_PRECISION* fsr_scalar_fluxes = new FP_PRECISION[_num_FSRs * _num_groups];
    cudaMemcpy((void*)fsr_scalar_fluxes, (void*)_scalar_flux,
	       _num_FSRs * _num_groups * sizeof(FP_PRECISION),
	       cudaMemcpyDeviceToHost);

    return _scalar_flux;
}


/**
 * @brief Return an array indexed by flatsourceregion IDs with the
 *        corresponding flatsourceregion power.
 * @return an array of flatsourceregion powers
 */
FP_PRECISION* DeviceSolver::getFSRPowers() {
    if (_FSRs_to_powers == NULL)
        log_printf(ERROR, "Unable to returns the device solver's FSR power "
		   "array since it has not yet been allocated in memory");

    return _FSRs_to_powers;
}


/**
 * @brief Return an array indexed by flatsourceregion IDs with the
 *        corresponding pin cell power.
 * @return an array of flatsourceregion pin powers
 */
FP_PRECISION* DeviceSolver::getFSRPinPowers() {
    if (_FSRs_to_pin_powers == NULL)
        log_printf(ERROR, "Unable to returns the device solver's FSR pin power "
		   "array since it has not yet been allocated in memory");

    return _FSRs_to_pin_powers;
}


/**
 * @brief Sets the geometry for the solver.
 * @details The geometry must already have initialized flat source region maps
 *          and segmentized the trackgenerator's tracks.
 * @param geometry a pointer to a geometry
 */
void DeviceSolver::setGeometry(Geometry* geometry) {

    if (geometry->getNumFSRs() == 0)
        log_printf(ERROR, "Unable to set the Geometry for the Solver "
		 "since the Geometry has not yet initialized flat "
		 "source regions");

    if (geometry->getNumEnergyGroups() == 0)
        log_printf(ERROR, "Unable to set the Geometry for the Solver "
		 "since the Geometry does not contain any energy groups");

    if (geometry->getNumMaterials() == 0)
        log_printf(ERROR, "Unable to set the Geometry for the Solver "
		 "since the Geometry does not contain any materials");

    _geometry = geometry;
    _num_FSRs = _geometry->getNumFSRs();
    _num_groups = _geometry->getNumEnergyGroups();
    _polar_times_groups = _num_polar * _num_groups;
    _num_materials = _geometry->getNumMaterials();
}


/**
 * @brief Sets the trackgenerator with characteristic tracks for the solver.
 * @details The trackgenerator must already have generated tracks and have
 *          segmentized them using the geometry.
 * @param track_generator a pointer to a trackgenerator
 */
void DeviceSolver::setTrackGenerator(TrackGenerator* track_generator) {

    if (!track_generator->containsTracks())
        log_printf(ERROR, "Unable to set the TrackGenerator for the Solver "
		 "since the TrackGenerator has not yet generated tracks");

    _track_generator = track_generator;
    _num_azim = _track_generator->getNumAzim() / 2;
    _host_tracks = _track_generator->getTracks();
    _num_tracks = _track_generator->getNumTracksArray();
    _tot_num_tracks = _track_generator->getNumTracks();
}


/**
 * @brief Sets the type of polar angle quadrature set to use (ie, TABUCHI 
 *        or LEONARD).
 * @param type the polar angle quadrature type
 */
void DeviceSolver::setPolarQuadratureType(quadratureType quadrature_type) {
    _quadrature_type = quadrature_type;
}


/**
 * @brief Sets the number of polar angles to use (only 1, 2, or 3 currently
 *        supported).
 * @param num_polar the number of polar angles
 */
void DeviceSolver::setNumPolarAngles(int num_polar) {

    if (num_polar <= 0)
        log_printf(ERROR, "Unable to set the Solver's number of polar angles "
		   "to %d since this is a negative number", num_polar);

    if (num_polar > 3)
        log_printf(ERROR, "Unable to set the DeviceSolver's number of polar "
		   "angles to %d since this is not a supported value (only 1, "
		   "2 or 3 are currently supported)", num_polar);

    _num_polar = num_polar;
    _two_times_num_polar = 2 * _num_polar;
    _polar_times_groups = _num_polar * _num_groups;
}


/**
 * @brief Sets the threshold for source convergence (>0)
 * @param source_thresh the threshold for source convergence
 */
void DeviceSolver::setSourceConvergenceThreshold(FP_PRECISION source_thresh) {

    if (source_thresh <= 0.0)
        log_printf(ERROR, "Unable to set the source convergence threshold to "
		   "%f since the threshold must be a positive number", 
		   source_thresh);

    _source_convergence_thresh = source_thresh;
}


/**
 * @brief Sets the threshold for flux convergence (>0) in fixed source
 *        iteration after the source has converged.
 * @param source_thresh the threshold for flux convergence
 */
void DeviceSolver::setFluxConvergenceThreshold(FP_PRECISION flux_thresh) {

    if (flux_thresh <= 0.0)
        log_printf(ERROR, "Unable to set the flux convergence threshold to "
	       "%f since the threshold must be a positive number",
	       flux_thresh);

    _flux_convergence_thresh = flux_thresh;
}


/**
 * @brief Sets the number of threadblocks (>0) for device kernels
 * @param num_blocks the number of threadblocks
 */
void DeviceSolver::setNumThreadBlocks(int num_blocks) {

    if (num_blocks < 0)
        log_printf(ERROR, "Unable to set the number of threadblocks to %d since "
		   "it is a negative number", num_blocks);

    _num_blocks = num_blocks;
}


/**
 * @brief Sets the number of threads per block (>0) for device kernels
 * @param num_threads the number of threads per block
 */
void DeviceSolver::setNumThreadsPerBlock(int num_threads) {

    if (num_threads < 0)
        log_printf(ERROR, "Unable to set the number of threads per block to %d "
		   "since it is a negative number", num_threads);

    _num_threads = num_threads;
}


/**
 * @brief Allocates and initializes all memory on the device.
 * @details Memory allocated includes data necessary for transport sweeping,
 *          including tracks, segments, flat source regions, materials, 
 *          and the polar quadrature.
 */
void DeviceSolver::allocateDeviceData() {

    log_printf(NORMAL, "Allocating memory for the device solver...");


    /**************************************************************************/
    /*                             Error checking                             */
    /**************************************************************************/

    if (_track_generator == NULL)
        log_printf(ERROR, "Unable to allocate memory on the device since "
		   "the device solver does not have a pointer to the "
		   "track generator");

    if (_geometry == NULL)
        log_printf(ERROR, "Unable to allocate memory on the device since "
		   "the device solver does not have a pointer to the geometry");


    /**************************************************************************/
    /*                     Initiailze each type of memory                     */
    /**************************************************************************/

    initializeHostMemory();
    initializeGlobalMemory();
    initializeConstantMemory();
    initializePinnedMemory();

    //FIXME: Need to compute prefactor hashtable on the device
    precomputePrefactors();

    return;
}


/**
 * @brief Allocates memory for the solver on the host.
 * @details Memory allocation includes the polar quadrature, 
 */
void DeviceSolver::initializeHostMemory() {

    log_printf(INFO, "Initializing host memory for the device solver...");

    /* Initialize the a polar quadrature object on the host */
    initializePolarQuadrature();

    /* Initialize arrays of FSR powers and pin powers */
    initializePowerArrays();
}


/**
 * @brief Creates a polar quadrature object for the solver.
 */
void DeviceSolver::initializePolarQuadrature() {

    /* Deletes the old quadrature if one existed */
    if (_quad != NULL)
        delete _quad;

    _quad = new Quadrature(_quadrature_type, _num_polar);
    _polar_times_groups = _num_groups * _num_polar;
}


/**
 * @brief Allocates memory for flatsourceregion power arrays.
 * @details Deletes memory for power arrays if they were allocated from
 *          previous simulation.
 */
void DeviceSolver::initializePowerArrays() {

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
 * @brief
 * @details
 */
void DeviceSolver::initializeGlobalMemory() {

    log_printf(INFO, "Initializing global memory on the device...");

    initializeFSRs();
    initializeMaterials();
    initializeTracks();
    initializeFluxArrays();
    initializeSourceArrays();
    initializePowerArrays();
    initializeThrustVectors();
}


/**
 * This is a helper method for the allocateDeviceMemory method. It
 * initializes an array of dev_flatsourceregion structs on the host
 * with the appropriate values (volume, region uid, and material uid)
 * so that allocateDeviceMemory can copy the array in full to the device.
 */
void DeviceSolver::initializeFSRs() {

    log_printf(INFO, "Initializing FSRs on the device...");

    /* Delete old FSRs array if it exists */
    if (_FSRs != NULL)
        cudaFree(_FSRs);

    /* Allocate memory for all tracks and track offset indices on the device */
    try{

        /* Allocate memory on device for FSRs */
        cudaMalloc((void**)&_FSRs, _num_FSRs * sizeof(dev_flatsourceregion));

	/* Create a temporary FSR array to populate and then copy to device */
	dev_flatsourceregion* temp_FSRs = new dev_flatsourceregion[_num_FSRs];

	/* Get the array indexed by FSR IDs with material ID values */
	int* FSRs_to_materials = _geometry->getFSRtoMaterialMap();

	/* Iterate over all FSRs and set the UIDs and material IDs */
	for (int r=0; r < _num_FSRs; r++) {
	    temp_FSRs[r]._uid = r;
	    temp_FSRs[r]._material_uid = FSRs_to_materials[r];
	}

	/* Initialize each FSRs volume to 0 to avoid NaNs */
	for (int r=0; r < _num_FSRs; r++)
	    temp_FSRs[r]._volume = 0.0;

	Track* track;
	segment* seg;
	dev_flatsourceregion* fsr;

	double* azim_weights = _track_generator->getAzimWeights();


	/* Set each FSR's volume by accumulating the total length of all
	   tracks inside the FSR. Iterate over azimuthal angle, track, segment */
	for (int i=0; i < _num_azim; i++) {
	    for (int j=0; j < _num_tracks[i]; j++) {
	        track = &_track_generator->getTracks()[i][j];

		/* Iterate over the track's segments to update FSR volumes */
		for (int s = 0; s < track->getNumSegments(); s++) {
		    seg = track->getSegment(s);
		    fsr = &temp_FSRs[seg->_region_id];
		    fsr->_volume += seg->_length * azim_weights[i];
		}
	    }
	}

	/* Copy the temporary array of FSRs to the device */
	cudaMemcpy((void*)_FSRs, (void*)temp_FSRs, 
		   _num_FSRs * sizeof(dev_flatsourceregion), 
		   cudaMemcpyHostToDevice);

	/* Free the temporary array of FSRs on the host */
	free(temp_FSRs);
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's flat "
		   "source regions on the device. Backtrace:%s", e.what());
    }
}


/**
 * @brief
 * @details
 */
void DeviceSolver::initializeMaterials() {

    log_printf(INFO, "Initializing materials on the device...");

    /* Delete old materials array if it exists */
    if (_materials != NULL)
        cudaFree(_materials);

    /* Allocate memory for all tracks and track offset indices on the device */
    try{

	std::map<short int, Material*> host_materials=_geometry->getMaterials();
	std::map<short int, Material*>::iterator iter;

        /* Iterate through all materials and clone them on the device */
        cudaMalloc((void**)&_materials, _num_materials * sizeof(dev_material));
	for (iter = host_materials.begin(); iter != host_materials.end(); ++iter) {
	    cloneOnDevice(iter->second, &_materials[iter->second->getUid()]);
	}

    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the device solver's "
		   "materials. Backtrace:%s", e.what());
    }
}


/**
 * @brief
 * @details
 */
void DeviceSolver::initializeTracks() {

    log_printf(INFO, "Initializing tracks on the device...");

    /* Delete old tracks array if it exists */
    if (_dev_tracks != NULL)
        cudaFree(_dev_tracks);

    /* Delete old track index offsets array if it exists */
    if (_track_index_offsets != NULL)
        delete [] _track_index_offsets;

    /* Allocate array of tracks */
    cudaMalloc((void**)&_dev_tracks, _tot_num_tracks * sizeof(dev_track));

    /* An array of the cumulative number of tracks for each azimuthal angle */
    _track_index_offsets = new int[_num_azim];

    /* Allocate memory for all tracks and track offset indices on the device */
    try{

        /* Iterate through all tracks and clone them on the device */
        int counter = 0;
	int index;
	for (int i=0; i < _num_azim; i++) {

            _track_index_offsets[i] = counter;
  
	    for (int j=0; j < _num_tracks[i]; j++) {

	        /* Clone this track on the device */
	        cloneTrack(&_host_tracks[i][j], &_dev_tracks[counter]);

		/* Make track reflective */
		index = computeScalarTrackIndex(_host_tracks[i][j].getTrackInI(),
					       _host_tracks[i][j].getTrackInJ());
		cudaMemcpy((void*)&_dev_tracks[counter]._track_in,
			   (void*)&index, sizeof(int), cudaMemcpyHostToDevice);

		index = computeScalarTrackIndex(_host_tracks[i][j].getTrackOutI(),					_host_tracks[i][j].getTrackOutJ());
		cudaMemcpy((void*)&_dev_tracks[counter]._track_out, 
			   (void*)&index, sizeof(int), cudaMemcpyHostToDevice);

		counter++;
	    }
	}

	/* Copy the cumulative track index offsets for each azimuthal angle */
	cudaMemcpy((void*)&_track_index_offsets[_num_azim], (void*)&counter, 
		   sizeof(int), cudaMemcpyHostToDevice);

    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's tracks "
		   "on the device. Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory for track boundary angular fluxes and 
 *        flatsourceregion scalar fluxes on the device.
 * @details Deletes memory for old flux arrays if they were allocated from
 *          previous simulation.
 */
void DeviceSolver::initializeFluxArrays() {

    log_printf(INFO, "Initializing flux arrays on the device...");

    /* Delete old flux arrays if they exist */
    if (_boundary_flux != NULL)
        cudaFree(_boundary_flux);
    if (_scalar_flux != NULL)
        cudaFree(_scalar_flux);
    if (_old_scalar_flux != NULL)
        cudaFree(_old_scalar_flux);

    /* Allocate memory for all flux arrays on the device */
    try{
        cudaMalloc((void**)&_boundary_flux,
		   2*_tot_num_tracks * _polar_times_groups*sizeof(FP_PRECISION));
        cudaMalloc((void**)&_scalar_flux, 
		   _num_FSRs * _num_groups * sizeof(FP_PRECISION));
        cudaMalloc((void**)&_old_scalar_flux, 
		   _num_FSRs * _num_groups * sizeof(FP_PRECISION));
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's fluxes "
		   "on the device. Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory for flatsourceregion source arrays on the device.
 * @details Deletes memory for old source arrays if they were allocated from
 *          previous simulation.
 */
void DeviceSolver::initializeSourceArrays() {

    log_printf(INFO, "Initializing source arrays on the device...");

    /* Delete old sources arrays if they exist */
    if (_source != NULL)
        cudaFree(_source);
    if (_old_source != NULL)
        cudaFree(_old_source);
    if (_ratios != NULL)
        cudaFree(_ratios);

    /* Allocate memory for all source arrays on the device */
    try{

        cudaMalloc((void**)&_source, 
		   _num_FSRs * _num_groups * sizeof(FP_PRECISION));
	cudaMalloc((void**)&_old_source,
		   _num_FSRs * _num_groups * sizeof(FP_PRECISION));
	cudaMalloc((void**)&_ratios,
		   _num_FSRs * _num_groups * sizeof(FP_PRECISION));
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's flat "
		   "source region sources array on the device. "
		   "Backtrace:%s", e.what());
    }
}


/**
 * @brief
 * @details
 */
void DeviceSolver::initializeThrustVectors() {

    log_printf(INFO, "Initializing thrust vectors on the device...");

    /* Delete old vectors if they exist */
    if (_fission_source != NULL) {
        _fission_source = NULL;
        _fission_source_vec.clear();
    }
    if (_tot_abs != NULL) {
        _tot_abs = NULL;
        _tot_abs_vec.clear();
    }
    if (_tot_fission != NULL) {
        _tot_fission = NULL;
        _tot_fission_vec.clear();
    }
    if (_source_residual != NULL) {
        _source_residual = NULL;
        _source_residual_vec.clear();
    }


    /* Allocate memory for fission, absorption and source vectors on device */
    try{
        /* Allocate fission source array on device */
        _fission_source_vec.resize(_num_blocks * _num_threads);
	_fission_source = thrust::raw_pointer_cast(&_fission_source_vec[0]);
      
	/* Allocate total absorption reaction rate array on device */
	_tot_abs_vec.resize(_num_blocks * _num_threads);
	_tot_abs = thrust::raw_pointer_cast(&_tot_abs_vec[0]);

	/* Allocate fission reaction rate array on device */
	_tot_fission_vec.resize(_num_blocks * _num_threads);
	_tot_fission = thrust::raw_pointer_cast(&_tot_fission_vec[0]);

	/* Allocate source residual array on device */
	_source_residual_vec.resize(_num_blocks * _num_threads);
	_source_residual = thrust::raw_pointer_cast(&_source_residual_vec[0]);
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's "
		   "Thrust vectors.Backtrace:%s", e.what());
    }
}


/**
 * @brief Initializes data in constant memory on the device.
 * @details
 */
void DeviceSolver::initializeConstantMemory() {

    log_printf(INFO, "Initializing constant memory on the device...");

    /* Number of azimuthal angles */
    cudaMemcpyToSymbol(_num_azim_devc, (void*)&_num_azim, sizeof(int), 0, 
		       cudaMemcpyHostToDevice);

    /* Number of energy groups */
    cudaMemcpyToSymbol(_num_groups_devc, (void*)&_num_groups, sizeof(int), 0,
		       cudaMemcpyHostToDevice);

    /* Number of flat source regions */
    cudaMemcpyToSymbol(_num_FSRs_devc, (void*)&_num_FSRs, sizeof(int), 0,
		       cudaMemcpyHostToDevice);

    /* Number of polar angles */
    cudaMemcpyToSymbol(_num_polar_devc, (void*)&_num_polar, sizeof(int), 0,
		       cudaMemcpyHostToDevice);

    /* Twice the number of polar angles */
    cudaMemcpyToSymbol(_two_times_num_polar_devc, (void*)&_two_times_num_polar, 
		       sizeof(int), 0, cudaMemcpyHostToDevice);

    /* Number of polar angles times energy groups */
    cudaMemcpyToSymbol(_polar_times_groups_devc, (void*)&_polar_times_groups, 
		       sizeof(int), 0, cudaMemcpyHostToDevice);

    /* Compute polar times azimuthal angle weights */
    FP_PRECISION* polar_weights =
        (FP_PRECISION*)malloc(_num_polar * _num_azim * sizeof(FP_PRECISION));
    FP_PRECISION* multiples = _quad->getMultiples();
    double* azim_weights = _track_generator->getAzimWeights();

    for (int i=0; i < _num_azim; i++) {
        for (int j=0; j < _num_polar; j++)
	    polar_weights[i*_num_polar+j] = azim_weights[i]*multiples[j]*FOUR_PI;
    }

    cudaMemcpyToSymbol(_polar_weights_devc, (void*)polar_weights,
		       _num_polar * _num_azim * sizeof(FP_PRECISION),
		       0, cudaMemcpyHostToDevice);
    free(polar_weights);

    /* Array of number of tracks for each azimuthal angles */
    cudaMemcpyToSymbol(_num_tracks_devc, (void*)&_num_tracks, 
		       _num_azim * sizeof(int), 0, cudaMemcpyHostToDevice);
    
    /* Total number of tracks */
    cudaMemcpyToSymbol(_tot_num_tracks_devc, (void*)&_tot_num_tracks,
		       sizeof(int), 0, cudaMemcpyHostToDevice);

    /* Copy the cumulative index offset for the current azimuthal angle */
    cudaMemcpyToSymbol(_track_index_offsets_devc, 
		       (void*)&_track_index_offsets, 
		       _num_azim * sizeof(int), 0, cudaMemcpyHostToDevice);

}


/**
 * @brief
 * @details
 */
void DeviceSolver::initializePinnedMemory() {

    log_printf(INFO, "Initializing pinned memory on the device...");

    /* Pinned host memory for keff */
    unsigned int flags = cudaHostAllocWriteCombined;
    cudaHostAlloc((void**)&_k_eff_pinned, sizeof(FP_PRECISION), flags);
    cudaHostAlloc((void**)&_norm_factor_pinned, sizeof(FP_PRECISION), flags);
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
int DeviceSolver::computeScalarTrackIndex(int i, int j) {

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
 *        for each polar angle and copies the table to the device. 
 * @details This method will generate a hashmap which contains values of the 
 *          pre-factor for specific segment lengths (the keys into the hashmap).
 */
void DeviceSolver::precomputePrefactors(){

    log_printf(INFO, "Building exponential prefactor hashtable on device...");

    /* Set size of prefactor array */
    int num_array_values = 10 * sqrt(1. / (8. * _source_convergence_thresh));
    FP_PRECISION prefactor_spacing = 10. / num_array_values;
    FP_PRECISION inverse_prefactor_spacing = 1.0 / prefactor_spacing;
    int prefactor_array_size = _two_times_num_polar * num_array_values;
    int prefactor_max_index = prefactor_array_size - _two_times_num_polar - 1;
    
    /* allocate arrays */
    FP_PRECISION* prefactor_array = new FP_PRECISION[prefactor_array_size];
    
    FP_PRECISION expon;
    FP_PRECISION intercept;
    FP_PRECISION slope;


    /* Create prefactor array */
    for (int i = 0; i < num_array_values; i ++){
        for (int p = 0; p < _num_polar; p++){
	    expon = exp(- (i * prefactor_spacing) / _quad->getSinTheta(p));
	    slope = - expon / _quad->getSinTheta(p);
	    intercept = expon * (1 + (i * prefactor_spacing) /
				 _quad->getSinTheta(p));
	    prefactor_array[_two_times_num_polar * i + 2 * p] = slope;
	    prefactor_array[_two_times_num_polar * i + 2 * p + 1] = intercept;
	}
    }

    /* Allocate memory for the prefactor array on the device */
    cudaMalloc((void**)&_prefactor_array, 
	       prefactor_array_size * sizeof(FP_PRECISION));

    /* Copy prefactor array to the device */
    cudaMemcpy((void*)_prefactor_array, (void*)prefactor_array, 
	       prefactor_array_size * sizeof(FP_PRECISION),
	       cudaMemcpyHostToDevice);

    /* Copy prefactor array size and spacing to constant memory on the device */
    cudaMemcpyToSymbol(_prefactor_spacing_devc, (void*)&prefactor_spacing, 
		       sizeof(int), 0, cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(_inverse_prefactor_spacing_devc, 
		       (void*)&inverse_prefactor_spacing, 
		       sizeof(int), 0, cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(_prefactor_max_index_devc, (void*)&prefactor_max_index,
		       sizeof(int), 0, cudaMemcpyHostToDevice);

    free(prefactor_array);

    return;
}


/**
 * @brief Checks that each flat source region has at least one segment within 
 *        it and if not, throw an exception and prints an error message.
 */
void DeviceSolver::checkTrackSpacing() {

    log_printf(INFO, "Checking track spacing...");

    int* FSR_segment_tallies = new int[_num_FSRs];
    std::vector<segment*> segments;
    std::vector<segment*>::iterator iter;
    Cell* cell;

    /* Set each tally to zero to begin with */
    for (int r=0; r < _num_FSRs; r++)
        FSR_segment_tallies[r] = 0;

    /* Iterate over all azimuthal angles, all tracks, and all segments
     * and tally each segment in the corresponding FSR */
    for (int i=0; i < _num_azim; i++) {
        for (int j=0; j < _num_tracks[i]; j++) {
	    segments = _host_tracks[i][j].getSegments();

            for (iter=segments.begin(); iter != segments.end(); ++iter)
	        FSR_segment_tallies[(*iter)->_region_id]++;
	}
    }

    /* Loop over all FSRs and if one FSR does not have tracks in it, print
     * error message to the screen and exit program */
    for (int r=0; r < _num_FSRs; r++) {
        if (FSR_segment_tallies[r] == 0) {
	    cell = _geometry->findCellContainingFSR(r);
	    log_printf(ERROR, "No tracks were tallied inside FSR id = %d which "
		       "is cell id = %d. Please reduce your track spacing,"
		       " increase the number of azimuthal angles, or increase "
		       "the size of the flat source regions", r, cell->getId());
	}
    }

    delete [] FSR_segment_tallies;
}


FP_PRECISION DeviceSolver::convergeSource(int max_iterations, int B, int T){
  
    /* Error checking */
    if (_geometry == NULL)
        log_printf(ERROR, "The DeviceSolver is unable to converge the source "
		   "since it does not contain a Geometry");
    if (_track_generator == NULL)
        log_printf(ERROR, "The DeviceSolver is unable to converge the source "
		   "since it does not contain a TrackGenerator");

    /** A reference to _k_eff_pinned in the device's address space */
    FP_PRECISION* k_eff_pointer;
    FP_PRECISION* norm_factor_pointer;

    /* Get keff pointer from the device */
    cudaHostGetDevicePointer((void**)&k_eff_pointer, (void*)_k_eff_pinned, 0);
    cudaHostGetDevicePointer((void**)&norm_factor_pointer, 
			     (void*)_norm_factor_pinned, 0);

    FP_PRECISION residual = 0.0;

    setNumThreadBlocks(B);
    setNumThreadsPerBlock(T);

    /* Initialize data structures on the device */
    allocateDeviceData();

    /* Counter for the number of iterations to converge the source */
    _num_iterations = 0;

    /* An initial guess for the eigenvalue */
    *_k_eff_pinned = 1.0;

    /* Check that each FSR has at least one segment crossing it */
    checkTrackSpacing();

    /* Set scalar flux to unity for each region */
    flattenFSRFluxes<<<_num_blocks, _num_threads>>>(_scalar_flux, 
						    _old_scalar_flux, 1.0);
    flattenFSRSources<<<_num_blocks, _num_threads>>>(_source, _old_source, 1.0);
    zeroTrackFluxes<<<_num_blocks, _num_threads>>>(_boundary_flux);

    log_printf(NORMAL, "Converging the source on the device...");

    /* Source iteration loop */
    for (int i=0; i < max_iterations; i++) {

        log_printf(NORMAL, "Iteration %d on host: \tk_eff = %1.6f"
		 "\tres = %1.3E", i, *_k_eff_pinned, residual);

	computeFissionSource<<<_num_blocks, _num_threads, 
	  sizeof(FP_PRECISION)*_num_threads>>>(_FSRs, _materials, 
					       _scalar_flux, _fission_source);

	*_norm_factor_pinned = 1.0 / thrust::reduce(_fission_source_vec.begin(),
						    _fission_source_vec.end());

        log_printf(DEBUG, "norm_factor = %f", *_norm_factor_pinned);

	normalizeFluxes<<<_num_blocks, _num_threads>>>(_scalar_flux, 
						       _boundary_flux, 
						       *_norm_factor_pinned);

	//	residual = computeFSRSources();
	//	transportSweep(1);
	//	computeKeff();
	_num_iterations++;

	if (i > 1 && residual < _source_convergence_thresh){
	  //	    transportSweep(1000);
	    return *_k_eff_pinned;
	}
    }

    log_printf(WARNING, "Unable to converge the source after %d iterations",
	       max_iterations);

    return *_k_eff_pinned;
}
