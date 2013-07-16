#include "MICSolver.h"


MICSolver::MICSolver(Geometry* geom, TrackGenerator* track_generator) :
  Solver(geom, track_generator) {

    setNumThreads(1);

    _FSR_materials = NULL;
    
    _materials = NULL;
    _dev_tracks = NULL;
    _track_index_offsets = NULL;
    track_index_offsets = 0;

    materials = 0;
    dev_tracks = 0;

    if (track_generator != NULL)
        setTrackGenerator(track_generator);

    if (geom != NULL)
        setGeometry(geom);

    FSR_volumes = 0;
    FSR_materials = 0;

    boundary_flux = 0;
    scalar_flux = 0;
    old_scalar_flux = 0;
    fission_source = 0;
    source = 0;
    old_source = 0;
    ratios = 0;

    polar_weights = 0;
    prefactor_array = 0;
}



MICSolver::~MICSolver() { 

    if (_materials != NULL) {
        delete [] _materials;
	_materials = NULL;

	/*
	size_t materials = this->materials;
        #pragma offload target(mic) in(materials)
	{
	    dev_material* ptr = (dev_material*)(materials);
	    free(ptr);
	}
	*/
    }

    if (_dev_tracks != NULL) {
        delete [] _dev_tracks;
	_dev_tracks = NULL;

	/*	size_t dev_tracks = this->dev_tracks;
        #pragma offload target(mic) in(dev_tracks)
	{
	    dev_track* ptr = (dev_track*)(dev_tracks);
	    free(ptr);
	}
	*/
    }

    if (_track_index_offsets != NULL) {
        delete [] _track_index_offsets;
	_track_index_offsets = NULL;
	/*
	size_t track_index_offsets = this->track_index_offsets;
        #pragma offload target(mic) in(track_index_offsets)
	{
	    int* ptr = (int*)track_index_offsets;
	    free(ptr);
	}
	*/
    }
    /*
    if (FSR_materials != 0) {
	size_t FSR_materials = this->FSR_materials;
        #pragma offload target(mic) in(FSR_materials)
	{
	    free((int*)FSR_materials);
	}
    }

    if (FSR_volumes != 0) {
        size_t FSR_volumes = this->FSR_volumes;
        #pragma offload target(mic) in(FSR_volumes)
	{
	    free((FP_PRECISION*)FSR_volumes);
	}      
    }

    if (prefactor_array != 0) {
        size_t prefactor_array = this->prefactor_array;
        #pragma offload target(mic) in(prefactor_array)
	{
	  free((FP_PRECISION*)prefactor_array);
	}
    }

    if (polar_weights != 0) {
        size_t polar_weights = this->polar_weights;
        #pragma offload target(mic) in(polar_weights)
	{
	  free((FP_PRECISION*)polar_weights);
	}
    }


    size_t boundary_flux = this->boundary_flux;
    size_t scalar_flux = this->scalar_flux;
    size_t old_scalar_flux = this->old_scalar_flux;
    size_t fission_source = this->fission_source;
    size_t source = this->source;
    size_t old_source = this->old_source;
    size_t ratios = this->ratios;

    #pragma offload target(mic) in(boundary_flux) \
      in(scalar_flux) in(old_scalar_flux) in(fission_source) \
      in(source) in(old_source) in(ratios)
    {
        FP_PRECISION* ptr1 = (FP_PRECISION*)(boundary_flux);
	FP_PRECISION* ptr2 = (FP_PRECISION*)(scalar_flux);
	FP_PRECISION* ptr3 = (FP_PRECISION*)(old_scalar_flux);
	FP_PRECISION* ptr4 = (FP_PRECISION*)(fission_source);
	FP_PRECISION* ptr5 = (FP_PRECISION*)(source);
	FP_PRECISION* ptr6 = (FP_PRECISION*)(old_source);
	FP_PRECISION* ptr7 = (FP_PRECISION*)(ratios);

	free(ptr1);
	free(ptr2);
	free(ptr3);
	free(ptr4);
	free(ptr5);
	free(ptr6);
	free(ptr7);
    }
    */
}


/**
 * @brief Returns the number of shared memory OpenMP threads in use.
 * @return the number of threads
 */
int MICSolver::getNumThreads() {
    return _num_threads;
}


FP_PRECISION MICSolver::getFSRScalarFlux(int fsr_id, int energy_group) {

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

    return _scalar_flux(fsr_id,energy_group);
}


/**
 * @brief Return a 2D array indexed by flatsourceregion IDs and energy groups 
 *        which contains the corresponding fluxes for each flatsourceregion.
 * @return a 2D array of flatsourceregion scalar fluxes
 */
FP_PRECISION* MICSolver::getFSRScalarFluxes() {
    if (_scalar_flux == NULL)
        log_printf(ERROR, "Unable to returns the Solver's scalar flux array "
		 "since it has not yet been allocated in memory");

    return _scalar_flux;
}


/**
 * @brief Return an array indexed by flatsourceregion IDs with the
 *        corresponding flatsourceregion power.
 * @return an array of flatsourceregion powers
 */
FP_PRECISION* MICSolver::getFSRPowers() {
    if (_FSRs_to_powers == NULL)
        log_printf(ERROR, "Unable to returns the Solver's FSR power array "
		 "since it has not yet been allocated in memory");

    return _FSRs_to_powers;
}


/**
 * @brief Return an array indexed by flatsourceregion IDs with the
 *        corresponding pin cell power.
 * @return an array of flatsourceregion pin powers
 */
FP_PRECISION* MICSolver::getFSRPinPowers() {
    if (_FSRs_to_pin_powers == NULL)
        log_printf(ERROR, "Unable to returns the Solver's FSR pin power array "
		 "since it has not yet been allocated in memory");

    return _FSRs_to_pin_powers;
}


/**
 * @brief Sets the number of shared memory OpenMP threads to use (>0).
 * @details This method sets the number of threads to be no larger than the
 *          number of complementary angle pairs - reflecting angles which form
 *          cyclical tracks - since that is how OpenMOC parallizes loops.
 * @param num_threads the number of threads
 */
void MICSolver::setNumThreads(int num_threads) {
    if (num_threads <= 0)
        log_printf(ERROR, "Unable to set the number of threads for the Solver "
		   "to %d since it is less than or equal to 0", num_threads);

    _num_threads = num_threads;

    /* Set the number of threads for OpenMP */
    omp_set_num_threads(_num_threads);

    /* Set the number of threads on the target device */
    #pragma offload target(mic)
    {
        omp_set_num_threads(_num_threads);
    }
}


/**
 * @brief Sets the geometry for the solver.
 * @details The geometry must already have initialized flat source region maps
 *          and segmentized the trackgenerator's tracks.
 * @param geometry a pointer to a geometry
 */
void MICSolver::setGeometry(Geometry* geometry) {
    Solver::setGeometry(geometry);
    initializeMaterials();
}


/**
 * @brief Sets the trackgenerator with characteristic tracks for the solver.
 * @details The trackgenerator must already have generated tracks and have
 *          segmentized them using the geometry.
 * @param track_generator a pointer to a trackgenerator
 */
void MICSolver::setTrackGenerator(TrackGenerator* track_generator) {
    Solver::setTrackGenerator(track_generator);
    initializeTracks();
}


/**
 * @brief
 * @details
 */
void MICSolver::initializeMaterials() {

    double* sigma_t;
    double* sigma_a;
    double* sigma_f;
    double* nu_sigma_f;
    double* sigma_s;
    double* chi;
    int id;
    int uid;
  

    log_printf(INFO, "Initializing materials on the MIC...");

    /* Delete old materials array if it exists */
    if (_materials != NULL)
        delete[] _materials;

    /* Allocate memory for all tracks and track offset indices on the device */
    try{

	std::map<short int, Material*> host_materials=_geometry->getMaterials();
	std::map<short int, Material*>::iterator iter;

	/* Allocate memory for the array of materials */
        _materials = new dev_material[_num_materials];

	/* Allocate persistent array of materials on the MIC */
	size_t materials;
        #pragma offload target(mic) out(materials)
	{
	    dev_material* ptr = new dev_material[_num_materials];
	    materials = size_t(ptr);
	}

	this->materials = materials;

        /* Iterate through all materials and clone them into structures */
	for (iter=host_materials.begin(); iter != host_materials.end(); ++iter) {

	    uid = iter->second->getUid();
	    id = iter->second->getId();

	    cloneMaterialOnMIC(iter->second, &_materials[uid]);
	    
	    sigma_t = _materials[uid]._sigma_t;
	    sigma_a = _materials[uid]._sigma_a;
	    sigma_f = _materials[uid]._sigma_f;
	    nu_sigma_f = _materials[uid]._nu_sigma_f;
	    sigma_s = _materials[uid]._sigma_s;
	    chi = _materials[uid]._chi;

            #pragma offload target(mic) in(materials) in(id) in (uid) \
	      in(sigma_t : length(_num_groups) alloc_if(1) free_if(0)) \
	      in(sigma_a : length(_num_groups) alloc_if(1) free_if(0)) \
	      in(sigma_f : length(_num_groups) alloc_if(1) free_if(0)) \
	      in(nu_sigma_f : length(_num_groups) alloc_if(1) free_if(0)) \
	      in(sigma_s : length(_num_groups*_num_groups) alloc_if(1) free_if(0)) \
	      in(chi : length(_num_groups) alloc_if(1) free_if(0))
	    {
	        dev_material* ptr = (dev_material*)(materials);
		ptr[uid]._uid = uid;
		ptr[uid]._id = id;
		ptr[uid]._sigma_t = sigma_t;
		ptr[uid]._sigma_a = sigma_a;
		ptr[uid]._sigma_f = sigma_f;
		ptr[uid]._nu_sigma_f = nu_sigma_f;
		ptr[uid]._sigma_s = sigma_s;
		ptr[uid]._chi = chi;
	    }
	}
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the MIC solver's "
		   "materials. Backtrace:%s", e.what());
    }
}


/**
 * @brief
 * @details
 */
void MICSolver::initializeTracks() {

    log_printf(INFO, "Initializing tracks on the MIC...");

    /* Delete old tracks array if it exists */
    if (_dev_tracks != NULL)
        delete [] _dev_tracks;

    /* Delete old track index offsets array if it exists */
    if (_track_index_offsets != NULL)
        delete [] _track_index_offsets;

    /* Allocate array of tracks and memory for track index offsets */
    _dev_tracks = new dev_track[_tot_num_tracks];
    _track_index_offsets = new int[_num_azim+1];

    size_t dev_tracks;
    #pragma offload target(mic) out(dev_tracks)
    {
        dev_track* ptr = new dev_track[_tot_num_tracks];
        dev_tracks = size_t(ptr);
    }

    this->dev_tracks = dev_tracks;

    /* Allocate memory for all tracks and track offset indices on the device */
    try{

        /* Iterate through all tracks and clone them on the device */
        int counter = 0;
	int index;
	for (int i=0; i < _num_azim; i++) {

            _track_index_offsets[i] = counter;

	    for (int j=0; j < _num_tracks[i]; j++) {

	        /* Clone this track on the device */
	        cloneTrackOnMIC(&_tracks[i][j], &_dev_tracks[counter]);

		/* Make track reflective */
		index = computeScalarTrackIndex(_tracks[i][j].getTrackInI(),
					       _tracks[i][j].getTrackInJ());
		_dev_tracks[counter]._track_in = index;

		index = computeScalarTrackIndex(_tracks[i][j].getTrackOutI(), 
						_tracks[i][j].getTrackOutJ());
		_dev_tracks[counter]._track_out = index;

		counter++;



		size_t dev_tracks = this->dev_tracks;

		int uid = _tracks[i][j].getUid();
		int azim_angle_index = _tracks[i][j].getAzimAngleIndex();
		bool refl_in = _tracks[i][j].isReflIn();
		bool refl_out = _tracks[i][j].isReflOut();
		bool bc_in = _tracks[i][j].getBCIn();
		bool bc_out = _tracks[i][j].getBCOut();

		int track_in = computeScalarTrackIndex(_tracks[i][j].getTrackInI(),
					       _tracks[i][j].getTrackInJ());
		int track_out = computeScalarTrackIndex(
						     _tracks[i][j].getTrackOutI(),
						     _tracks[i][j].getTrackOutJ());

		int num_segments = _tracks[i][j].getNumSegments();

                #pragma offload target(mic) in(uid) in(azim_angle_index) \
		  in(num_segments) in(refl_in) in(refl_out) in(bc_in) \
		  in(bc_out) in(track_in) in(track_out) in(dev_tracks)
		{
		    dev_track* ptr = (dev_track*)(dev_tracks);
		    ptr[uid]._uid = uid;
		    ptr[uid]._azim_angle_index = azim_angle_index;
		    ptr[uid]._refl_in = refl_in;
		    ptr[uid]._refl_out = refl_out;
		    ptr[uid]._bc_in = bc_in;
		    ptr[uid]._bc_out = bc_out;
		    ptr[uid]._track_in = track_in;
		    ptr[uid]._track_out = track_out;
		    ptr[uid]._num_segments = num_segments;
		    ptr[uid]._segments = new dev_segment[num_segments];
		}

		segment* curr_segment;
		FP_PRECISION length;
		int region_uid;
		int material_uid;

		for (int s=0; s < _tracks[i][j].getNumSegments(); s++) {
		    curr_segment = _tracks[i][j].getSegment(s);
		    length = curr_segment->_length;
		    region_uid = curr_segment->_region_id;
		    material_uid = curr_segment->_material->getUid();

                    #pragma offload target(mic) in(dev_tracks) in(uid) \
		      in(s) in(length) in(region_uid) in(material_uid)
		    {
		        dev_track* ptr = (dev_track*)(dev_tracks);
		        ptr[uid]._segments[s]._length = length;
		        ptr[uid]._segments[s]._region_uid = region_uid;
			ptr[uid]._segments[s]._material_uid = material_uid;
		    }
		}
	    }
	}

	_track_index_offsets[_num_azim] = counter;

	size_t track_index_offsets;

	int* offsets = _track_index_offsets;

        #pragma offload target(mic) out(track_index_offsets) \
	  in(offsets : length(_num_azim+1))
	{
	    int* ptr = new int[_num_azim+1];
	    memcpy(ptr, offsets, (_num_azim+1)*sizeof(int));
	    track_index_offsets = size_t(ptr);
	}

	this->track_index_offsets = track_index_offsets;
    }

    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's tracks "
		   "on the MIC. Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory for track boundary angular fluxes and 
 *        flatsourceregion scalar fluxes.
 * @details Deletes memory for old flux arrays if they were allocated from
 *          previous simulation.
 */
void MICSolver::initializeFluxArrays() {

    /* Delete old flux arrays if they exist */
    if (_boundary_flux != NULL)
        delete [] _boundary_flux;
    if (_scalar_flux != NULL)
        delete [] _scalar_flux;
    if (_old_scalar_flux != NULL)
        delete [] _old_scalar_flux;

    /* Allocate memory for all flux arrays */
    try{
        _boundary_flux = new FP_PRECISION[2*_tot_num_tracks*_polar_times_groups];
	_scalar_flux = new FP_PRECISION[_num_FSRs*_num_groups];
	_old_scalar_flux = new FP_PRECISION[_num_FSRs*_num_groups];

	size_t boundary_flux;
	size_t scalar_flux;
	size_t old_scalar_flux;

        #pragma offload target(mic) out(boundary_flux) \
          out(scalar_flux) out(old_scalar_flux)
	{
	  FP_PRECISION* ptr1 = new FP_PRECISION[2*_tot_num_tracks*_polar_times_groups];
	  FP_PRECISION* ptr2 = new FP_PRECISION[_num_FSRs*_num_groups];
	  FP_PRECISION* ptr3 = new FP_PRECISION[_num_FSRs*_num_groups];

	  boundary_flux = size_t(ptr1);
	  scalar_flux = size_t(ptr2);
	  old_scalar_flux = size_t(ptr3);
	}

	this->boundary_flux = boundary_flux;
	this->scalar_flux = scalar_flux;
	this->old_scalar_flux = old_scalar_flux;
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's fluxes. "
		   "Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory for flatsourceregion source arrays.
 * @details Deletes memory for old source arrays if they were allocated from
 *          previous simulation.
 */
void MICSolver::initializeSourceArrays() {

    /* Delete old sources arrays if they exist */
    if (_fission_source != NULL)
        delete [] _fission_source;
    if (_source != NULL)
        delete [] _source;
    if (_old_source != NULL)
        delete [] _old_source;
    if (_ratios != NULL)
        delete [] _ratios;

    /* Allocate memory for all source arrays */
    try{
	_fission_source = new FP_PRECISION[_num_FSRs*_num_groups];
	_source = new FP_PRECISION[_num_FSRs*_num_groups];
	_old_source = new FP_PRECISION[_num_FSRs*_num_groups];
	_ratios = new FP_PRECISION[_num_FSRs*_num_groups];


	size_t fission_source;
	size_t source;
	size_t old_source;
	size_t ratios;

        #pragma offload target(mic) out(fission_source) \
	  out(source) out(old_source) out(ratios)
	{
	    FP_PRECISION* ptr1 = new FP_PRECISION[_num_FSRs*_num_groups];
	    FP_PRECISION* ptr2 = new FP_PRECISION[_num_FSRs*_num_groups];
	    FP_PRECISION* ptr3 = new FP_PRECISION[_num_FSRs*_num_groups];
	    FP_PRECISION* ptr4 = new FP_PRECISION[_num_FSRs*_num_groups];

	    fission_source = size_t(ptr1);
	    source = size_t(ptr2);
	    old_source = size_t(ptr3);
	    ratios = size_t(ptr4);
	}

	this->fission_source = fission_source;
	this->source = source;
	this->old_source = old_source;
	this->ratios = ratios;
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's flat "
		   "source region sources array. Backtrace:%s", e.what());
    }
}


/**
 * @brief Allocates memory for flatsourceregion power arrays.
 * @details Deletes memory for power arrays if they were allocated from
 *          previous simulation.
 */
void MICSolver::initializePowerArrays() {

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
        log_printf(ERROR, "Could not allocate memory for the solver's FSR "
		   "power arrays. Backtrace:%s", e.what());
    }
}


/**
 * @brief Creates a polar quadrature object for the solver.
 */
void MICSolver::initializePolarQuadrature() {
    /* Deletes the old quadrature if one existed */
    if (_quad != NULL)
        delete _quad;

    _quad = new Quadrature(_quadrature_type, _num_polar);
    _polar_times_groups = _num_groups * _num_polar;
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
#pragma offload_attribute(push, target(mic))
int MICSolver::computeScalarTrackIndex(int i, int j) {

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
#pragma offload_attribute(pop)



 
/**
 * @brief Pre-computes exponential pre-factors for each segment of each track 
 *        for each polar angle. 
 * @details This method will generate a hashmap which contains values of the 
 *          pre-factor for specific segment lengths (the keys into the hashmap).
 */
void MICSolver::precomputePrefactors() {

    /* Build exponential prefactor array based on table look up with linear 
     * interpolation */
    log_printf(INFO, "Building exponential prefactor hashtable...");

    FP_PRECISION azim_weight;

    _polar_weights = new FP_PRECISION[_num_azim*_num_polar];

    /* Precompute the total azimuthal weight for tracks at each polar angle */
    #pragma omp parallel for private (azim_weight)
    for (int i=0; i < _num_azim; i++) {
        azim_weight = _azim_weights[i];

        for (int p=0; p < _num_polar; p++)
	    _polar_weights(i,p) = azim_weight*_quad->getMultiple(p)*FOUR_PI;
    }

    /* Set size of prefactor array */
    int num_array_values = 10 * sqrt(1. / (8. * _source_convergence_thresh));
    _prefactor_spacing = 10. / num_array_values;
    _prefactor_array_size = _two_times_num_polar * num_array_values;
    _prefactor_max_index = _prefactor_array_size - _two_times_num_polar - 1.;
    
    log_printf(DEBUG, "Prefactor array size: %i, max index: %i",
	       _prefactor_array_size, _prefactor_max_index);

    /* allocate arrays */
    _prefactor_array = new FP_PRECISION[_prefactor_array_size];

    FP_PRECISION expon;
    FP_PRECISION intercept;
    FP_PRECISION slope;

    /* Create prefactor array */
    for (int i=0; i < num_array_values; i ++){
        for (int p=0; p < _num_polar; p++){
	    expon = exp(- (i * _prefactor_spacing) / _quad->getSinTheta(p));
	    slope = - expon / _quad->getSinTheta(p);
	    intercept = expon * (1 + (i * _prefactor_spacing) /
				 _quad->getSinTheta(p));
	    _prefactor_array[_two_times_num_polar * i + 2 * p] = slope;
	    _prefactor_array[_two_times_num_polar * i + 2 * p + 1] = intercept;
	}
    }

    /* Compute the reciprocal of the prefactor spacing */
    _inverse_prefactor_spacing = 1.0 / _prefactor_spacing;


    FP_PRECISION* _polar_weights = this->_polar_weights;
    FP_PRECISION* _prefactor_array = this->_prefactor_array;

    size_t polar_weights;
    size_t prefactor_array;

    #pragma offload target(mic) in(_polar_weights : length(_num_azim*_num_polar) \
				   alloc_if(1) free_if(0)) \
      in(_prefactor_array : length(_prefactor_array_size) alloc_if(1) free_if(0)) \
      out(polar_weights) out(prefactor_array)
    {
      polar_weights = size_t(_polar_weights);
      prefactor_array = size_t(_prefactor_array);
    }
      
    this->polar_weights = polar_weights;
    this->prefactor_array = prefactor_array;

    return;
}


/**
 * @brief Initializes each of the flatsourceregion objects inside the solver's
 *        array of flatsourceregions. 
 * @details This method assigns each flatsourceregion a unique, monotonically
 *          increasing ID, sets the material for each flatsourceregion, and 
 *          assigns a volume based on the cumulative length of all of the 
 *          segments inside the flatsourceregion.
 */
void MICSolver::initializeFSRs() {

    log_printf(INFO, "Initializing flat source regions...");

    /* Delete old FSRs array if it exists */
    if (_FSR_volumes != NULL)
        delete [] _FSR_volumes;

    if (_FSR_materials != NULL)
        delete [] _FSR_materials;

    _FSR_volumes = new FP_PRECISION[_num_FSRs];
    _FSR_materials = new int[_num_FSRs];

    int num_segments;
    segment* curr_segment;
    //    std::vector<segment*> segments;
    //    std::vector<segment*>::iterator iter;
    FP_PRECISION volume;
    CellBasic* cell;
    Material* material;
    Universe* univ_zero = _geometry->getUniverse(0);

    /* Initialize the FSR volumes to zero */
    memset(_FSR_volumes, FP_PRECISION(0.), _num_FSRs*sizeof(FP_PRECISION));

    /* Set each FSR's volume by accumulating the total length of all tracks
     * inside the FSR. Loop over azimuthal angle, track and segment. 
     * Note: this code region cannot be parallelized without a mutex lock
     * on FSR volume due to race conditions. */
    for (int i=0; i < _num_azim; i++) {
        for (int j=0; j < _num_tracks[i]; j++) {
	    num_segments = _tracks[i][j].getNumSegments();
	    //  	    segments = _tracks[i][j].getSegments();

	    //            for (iter=segments.begin(); iter != segments.end(); ++iter) {
	    for (int s=0; s < num_segments; s++) {
	        curr_segment = _tracks[i][j].getSegment(s);
	        volume = curr_segment->_length * _azim_weights[i];
		_FSR_volumes[curr_segment->_region_id] += volume;
	    }

	    //	        volume = (*iter)->_length * _azim_weights[i];
	    //		_FSR_volumes[(*iter)->_region_id] += volume;
	    //	    }
	}
    }

    /* Loop over all FSRs */
    #pragma omp parallel for private(cell, material)
    for (int r=0; r < _num_FSRs; r++) {

        /* Get the cell corresponding to this FSR from the geometry */
        cell = static_cast<CellBasic*>(_geometry->findCell(univ_zero, r));

	/* Get the cell's material and assign it to the FSR */
	material = _geometry->getMaterial(cell->getMaterial());
	_FSR_materials[r] = material->getUid();

	log_printf(DEBUG, "FSR id = %d has cell id = %d and material id = %d "
                  "and volume = %f", r, cell->getId(), 
                   _materials[_FSR_materials[r]]._uid, _FSR_volumes[r]);
    }

    size_t FSR_volumes;
    size_t FSR_materials;

    FP_PRECISION* _FSR_volumes = this->_FSR_volumes;
    int* _FSR_materials = this->_FSR_materials;

    #pragma offload target(mic) in(_FSR_materials : length(_num_FSRs) \
				   alloc_if(1) free_if(0)) \
      in(_FSR_volumes : length(_num_FSRs) alloc_if(1) free_if(0)) \
      out(FSR_volumes) out(FSR_materials)
    {
      FSR_volumes = size_t(_FSR_volumes);
      FSR_materials = size_t(_FSR_materials);
    }

    this->FSR_volumes = FSR_volumes;
    this->FSR_materials = FSR_materials;

    return;
}


/**
 * @brief Zero each track's boundary fluxes for each energy group and polar
 *        angle in the "forward" and "reverse" directions.
 */
MIC_ATTRIBUTE void MICSolver::zeroTrackFluxes() {

    /* Loop over azimuthal angle, track, polar angle, energy group
     * and set each track's incoming and outgoing flux to zero */
    #pragma offload target(mic) in(boundary_flux) 
    {
        FP_PRECISION* _boundary_flux = (FP_PRECISION*)(boundary_flux);

        #pragma omp parallel for
        for (int i=0; i < _tot_num_tracks; i++) {
	    for (int pe2=0; pe2 < 2*_polar_times_groups; pe2++)
	        _boundary_flux(i,pe2) = 0.0;
	}
    }

    return;
}


/**
 * @brief Set the scalar flux for each energy group inside each 
 *        flatsourceregion to a constant value.
 * @param value the value to assign to each flat source region flux
 */
MIC_ATTRIBUTE void MICSolver::flattenFSRFluxes(FP_PRECISION value) {

    /* Loop over all FSRs and energy groups */
    #pragma offload target(mic) in(scalar_flux, old_scalar_flux)
    {
        FP_PRECISION* _scalar_flux = (FP_PRECISION*)(scalar_flux);
	FP_PRECISION* _old_scalar_flux = (FP_PRECISION*)(old_scalar_flux);

        #pragma omp parallel for
        for (int r=0; r < _num_FSRs; r++) {
            for (int e=0; e < _num_groups; e++) {
	        _scalar_flux(r,e) = value;
		_old_scalar_flux(r,e) = value;
	    }
	}
    }

    return;
}


/**
 * @brief Set the source for each energy group inside each flatsourceregion
 *        to a constant value.
 * @param value the value to assign to each flat source region source
 */
void MICSolver::flattenFSRSources(FP_PRECISION value) {

    /* Loop over all FSRs and energy groups */
    #pragma offload target(mic) in(source, old_source)
    {

        FP_PRECISION* _source = (FP_PRECISION*)(source);
	FP_PRECISION* _old_source = (FP_PRECISION*)(old_source);

        #pragma omp parallel for
        for (int r=0; r < _num_FSRs; r++) {
	    for (int e=0; e < _num_groups; e++) {
	   
	        _source(r,e) = value;
		_old_source(r,e) = value;
	    }
	}
    }

    return;
}


/**
 * @brief Normalizes all flatsourceregion scalar fluxes and track boundary
 *        angular fluxes to the total fission source (times nu).
 */
void MICSolver::normalizeFluxes() {

    FP_PRECISION tot_fission_source;
    FP_PRECISION norm_factor;

    size_t materials = this->materials;

    /* Compute total fission source for each region, energy group */
    /* Loop over all FSRs and energy groups */
    #pragma offload target(mic) in(fission_source) in(scalar_flux) \
      in(FSR_volumes) in(FSR_materials) in(materials) inout(tot_fission_source)
    {
        FP_PRECISION volume;
	double* nu_sigma_f;

	FP_PRECISION* _scalar_flux = (FP_PRECISION*)(scalar_flux);
	FP_PRECISION* _fission_source = (FP_PRECISION*)(fission_source);
	FP_PRECISION* _FSR_volumes = (FP_PRECISION*)(FSR_volumes);
	int* _FSR_materials = (int*)(FSR_materials);
	dev_material* _materials = (dev_material*)(materials);

        #pragma omp parallel for private(volume, nu_sigma_f)
        for (int r=0; r < _num_FSRs; r++) {

	    /* Get pointers to important data structures */
	    nu_sigma_f = _materials[_FSR_materials[r]]._nu_sigma_f;
	    volume = _FSR_volumes[r];

	    for (int e=0; e < _num_groups; e++)
		_fission_source(r,e) = nu_sigma_f[e] * _scalar_flux(r,e) * volume;
	}

	/* Compute the total fission source */
        tot_fission_source = pairwise_sum<FP_PRECISION>(_fission_source, 
							    _num_FSRs*_num_groups);
    }

    /* Normalize scalar fluxes in each region */
    norm_factor = 1.0 / tot_fission_source;

    log_printf(DEBUG, "Normalization factor = %f", norm_factor);

    #pragma offload target(mic) in(norm_factor) in(scalar_flux) 
    {
        FP_PRECISION* _scalar_flux = (FP_PRECISION*)(scalar_flux);

        #pragma omp parallel for
        for (int r=0; r < _num_FSRs; r++) {
	    for (int e=0; e < _num_groups; e++) {
	        _scalar_flux(r,e) *= norm_factor;
	    }
	}
    }

    /* Normalize angular boundary fluxes for each track */
    #pragma offload target(mic) in(norm_factor) in(boundary_flux)
    {
        FP_PRECISION* _boundary_flux = (FP_PRECISION*)(boundary_flux);

        #pragma omp parallel for
        for (int i=0; i < _tot_num_tracks; i++) {
	    for (int pe2=0; pe2 < 2*_polar_times_groups; pe2++) {
	        _boundary_flux(i,pe2) *= norm_factor;
	    }
	}
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
 *                    \left(\frac{Q^i - Q^{i-1}{Q^i}\right)^2}{# FSRs}} \f$
 *
 * @return the residual between this source and the previous source
 */
FP_PRECISION MICSolver::computeFSRSources() {

    FP_PRECISION source_residual;

    #pragma offload target(mic) in(boundary_flux) in(FSR_materials) \
      in(materials) in(scalar_flux) in(source) in(old_source) in(ratios) \
      inout(source_residual)
    {
        FP_PRECISION* _boundary_flux = (FP_PRECISION*)(boundary_flux);
	FP_PRECISION* _scalar_flux = (FP_PRECISION*)(scalar_flux);
	FP_PRECISION* _source = (FP_PRECISION*)(source);
	FP_PRECISION* _old_source = (FP_PRECISION*)(old_source);
	FP_PRECISION* _ratios = (FP_PRECISION*)(ratios);
	int* _FSR_materials = (int*)(FSR_materials);
	dev_material* _materials = (dev_material*)(materials);

	FP_PRECISION* source_residuals = new FP_PRECISION[_num_groups*_num_FSRs];

	FP_PRECISION scatter_source;
	FP_PRECISION fission_source;
	double* nu_sigma_f;
	double* sigma_s;
	double* sigma_t;
	double* chi;
	dev_material* material;

        #pragma omp parallel for private(material, nu_sigma_f, \
	  chi, sigma_s, sigma_t, fission_source, scatter_source)
        for (int r=0; r < _num_FSRs; r++) {

	    FP_PRECISION* scatter_sources = new FP_PRECISION[_num_groups];
	    FP_PRECISION* fission_sources = new FP_PRECISION[_num_groups];

	    material = &_materials[_FSR_materials[r]];
	    nu_sigma_f = material->_nu_sigma_f;
	    chi = material->_chi;
	    sigma_s = material->_sigma_s;
	    sigma_t = material->_sigma_t;

	    /* Compute fission source for each group */
	    for (int e=0; e < _num_groups; e++)
	        fission_sources[e] = _scalar_flux(r,e) * nu_sigma_f[e];		

	    fission_source = pairwise_sum<FP_PRECISION>(fission_sources, 
							    _num_groups);
	
	    /* Compute total scattering source for group G */
	    for (int G=0; G < _num_groups; G++) {
	        scatter_source = 0;

		for (int g=0; g < _num_groups; g++)
		    scatter_sources[g] =sigma_s[G*_num_groups+g]*_scalar_flux(r,g);

		scatter_source = pairwise_sum<FP_PRECISION>(scatter_sources, 
                                                        _num_groups);

		/* Set the total source for region r in group G */
		_source(r,G) = ((1.0 / _k_eff) * fission_source *
				chi[G] + scatter_source) * ONE_OVER_FOUR_PI;

		_ratios(r,G) = _source(r,G) / sigma_t[G];

		/* Compute the norm of residual of the source in region, group */
		if (fabs(_source(r,G)) > 1E-10)
		    source_residuals(r,G) = pow((_source(r,G) - _old_source(r,G)) 
						/ _source(r,G), 2);
	    
		/* Update the old source */
		_old_source(r,G) = _source(r,G);
	    }

	    delete [] scatter_sources;
	    delete [] fission_sources;
	}

	source_residual = pairwise_sum<FP_PRECISION>(source_residuals, 
							 _num_FSRs*_num_groups);
    
	delete [] source_residuals;
    }

    /* Sum up the residuals from each group and in each region */
    source_residual = sqrt(source_residual / _num_FSRs);

    return source_residual;
}


/**
 * @brief Compute \f$ k_{eff} \f$ from the total fission and absorption rates.
 * @details This method computes the current approximation to the 
 *          multiplication factor on this iteration as follows:
 *          \f$ k_{eff} = \frac{\displaystyle\sum \displaystyle\sum \nu
 *                        \Sigma_f \Phi V}{\displaystyle\sum 
 *                        \displaystyle\sum \Sigma_a \Phi V} \f$
 */
void MICSolver::computeKeff() {

    double tot_abs = 0.0;
    double tot_fission = 0.0;

    #pragma offload target(mic) in(FSR_volumes) in(FSR_materials) \
      in(materials) in(scalar_flux) inout(tot_abs, tot_fission)
    {
        FP_PRECISION* _FSR_volumes = (FP_PRECISION*)(FSR_volumes);
	int* _FSR_materials = (int*)(FSR_materials);
	dev_material* _materials = (dev_material*)(materials);
	FP_PRECISION* _scalar_flux = (FP_PRECISION*)(scalar_flux);

	dev_material* material;
	double* sigma_a;
	double* nu_sigma_f;
	FP_PRECISION volume;
	
	FP_PRECISION* absorption_rates = new FP_PRECISION[_num_FSRs*_num_groups];
	FP_PRECISION* fission_rates = new FP_PRECISION[_num_FSRs*_num_groups];

        #pragma omp parallel for private(volume, material, sigma_a, nu_sigma_f)
	for (int r=0; r < _num_FSRs; r++) {

	    volume = _FSR_volumes[r];
	    material = &_materials[_FSR_materials[r]];
	    sigma_a = material->_sigma_a;
	    nu_sigma_f = material->_nu_sigma_f;
	    
	    for (int e=0; e < _num_groups; e++) {
	        absorption_rates[r*_num_groups+e] = sigma_a[e] * _scalar_flux(r,e) * volume;
		fission_rates[r*_num_groups+e] = nu_sigma_f[e] * _scalar_flux(r,e) * volume;
	    }
	}

	tot_abs = pairwise_sum<FP_PRECISION>(absorption_rates, _num_FSRs*_num_groups);
	tot_fission = pairwise_sum<FP_PRECISION>(fission_rates, _num_FSRs*_num_groups);
    
	delete [] absorption_rates;
	delete [] fission_rates;
    }

    _k_eff = tot_fission / (tot_abs + _leakage);

    printf("abs = %1.15f, fiss = %1.15f, leak = %1.15f, keff = %1.15f\n", 
	   tot_abs, tot_fission, _leakage, _k_eff);

    log_printf(DEBUG, "tot_abs = %f, tot_fission = %f, leakage = %f, "
	       "k_eff = %f", tot_abs, tot_fission, _leakage, _k_eff);

    return;
}


/**
 * @brief Checks if scalar flux has converged within the threshold.
 * @return true if converged, false otherwise
 */
#ifdef MIC
#pragma offload_attribute(push, target(mic))
#endif
bool MICSolver::isScalarFluxConverged() {

    bool converged = true;

    #pragma offload target(mic) in(scalar_flux) in(old_scalar_flux) \
      inout(converged)
    {
        FP_PRECISION* _scalar_flux = (FP_PRECISION*)(scalar_flux);
        FP_PRECISION* _old_scalar_flux = (FP_PRECISION*)(old_scalar_flux);

        #pragma omp parallel for
        for (int r=0; r < _num_FSRs; r++) {
	    for (int e=0; e < _num_groups; e++) {

	        if (fabs((_scalar_flux(r,e) - _old_scalar_flux(r,e)) /
			 _old_scalar_flux(r,e)) > _flux_convergence_thresh)

		  converged = false;

		/* Update old scalar flux */
		_old_scalar_flux(r,e) = _scalar_flux(r,e);
	    }
	}
    }

    return converged;
}
#ifdef MIC
#pragma offload_attribute(pop)
#endif


/**
 * This method performs on or more fixed source iterations by integrating
 * the flux along each track and updating the boundary fluxes for the
 * corresponding output track, while updating the scalar flux in each
 * flat source region
 * @param max_iterations the maximum number of iterations allowed
 */
void MICSolver::transportSweep(int max_iterations) {

    bool converged = false;

    log_printf(DEBUG, "Transport sweep with max_iterations = %d and "
	       "# threads = %d", max_iterations, _num_threads);

    FP_PRECISION leakage;

    /* Loop for until converged or max_iterations is reached */
    for (int i=0; i < max_iterations; i++) {

        /* Initialize the global total leakage tally to zero */
        leakage = 0.0;

        /* Initialize flux in each region to zero */
        flattenFSRFluxes(0.0);

        #pragma offload target(mic) in(dev_tracks) in(materials) \
	  in(scalar_flux) in(boundary_flux) in(ratios) in(prefactor_array) \
	  in(polar_weights) in(track_index_offsets) in(FSR_volumes) \
	  in(FSR_materials) inout(leakage)
        {

	    FP_PRECISION volume;

	    double* sigma_t;

	    /* Allocate memory for thread FSR scalar fluxes and leakages */
	    FP_PRECISION* thread_flux = new FP_PRECISION[_num_threads * 
							 _num_FSRs * 
							 _num_groups];
	    
	    FP_PRECISION* thread_leakage = new FP_PRECISION[2*_tot_num_tracks];

	    /* Initialize thread fluxes and leakages to zero */
	    memset(thread_flux, FP_PRECISION(0.), 
		   _num_threads*_num_FSRs*_num_groups*sizeof(FP_PRECISION));
	    memset(thread_leakage, FP_PRECISION(0.), 
		   2*_tot_num_tracks*sizeof(FP_PRECISION));

	    dev_track* _dev_tracks = (dev_track*)(dev_tracks);
	    dev_material* _materials = (dev_material*)(materials);
	    FP_PRECISION* _scalar_flux = (FP_PRECISION*)(scalar_flux);
	    FP_PRECISION* _boundary_flux = (FP_PRECISION*)(boundary_flux);
	    FP_PRECISION* _ratios = (FP_PRECISION*)(ratios);
	    FP_PRECISION* _prefactor_array = (FP_PRECISION*)(prefactor_array);
	    FP_PRECISION* _polar_weights = (FP_PRECISION*)(polar_weights);
	    int* _track_index_offsets = (int*)(track_index_offsets);
	    FP_PRECISION* _FSR_volumes = (FP_PRECISION*)(FSR_volumes);
	    int* _FSR_materials = (int*)(FSR_materials);

	    /* Loop over each thread and azimuthal angle.
	     * If we are using more than 1 thread then we create 
	     * separate threads for each pair of complementary  
	     * azimuthal angles - angles which wrap into cycles */
	    for (int t=0; t < 2; t++) {

	        int thread_id;

		int fsr_id;
		int track_out_id;
		
		FP_PRECISION fsr_flux;
		
		FP_PRECISION sigma_t_l;
		int index;
		
		bool bc;
		int start;
		int pe;
	    
		dev_track* curr_track;
		dev_segment* curr_segment;
		dev_material* curr_material;
		int num_segments;
		FP_PRECISION delta;

		int tid_max[2] = { _track_index_offsets[_num_azim / 2], 
				   _track_index_offsets[_num_azim]    };
		int tid_min[2] = {0, _track_index_offsets[_num_azim / 2] };

                #pragma omp parallel for private(track_out_id, bc, fsr_id, \
		  curr_track, curr_segment, num_segments, curr_material, \
	          sigma_t, fsr_flux, delta, pe, sigma_t_l, index, thread_id, start)
		for (int track_id=tid_min[t]; track_id < tid_max[t]; track_id++){
	        
		    thread_id = omp_get_thread_num();

		    /* Initialize local pointers to important data structures */
		    curr_track = &_dev_tracks[track_id];
		    num_segments = curr_track->_num_segments;

		    /* Loop over each segment in forward direction */
		    for (int s=0; s < num_segments; s++) {

		        curr_segment = &curr_track->_segments[s];
			fsr_id = curr_segment->_region_uid;
			curr_material = &_materials[curr_segment->_material_uid];
			sigma_t = curr_material->_sigma_t;

			/* Initialize polar angle and energy group counter */
			pe = 0;
			
			/* Loop over energy groups */
			for (int e=0; e < _num_groups; e++) {

			    fsr_flux = 0.;
			    sigma_t_l = sigma_t[e] * curr_segment->_length;
			    index = int(sigma_t_l * _inverse_prefactor_spacing)*_two_times_num_polar;

			    /* Loop over polar angles */
			    for (int p=0; p < _num_polar; p++){
			        delta = (_boundary_flux(track_id,pe) - 
					 _ratios(fsr_id,e)) * 
				        prefactor(index,p,sigma_t_l);
				fsr_flux += delta * _polar_weights[p];
				_boundary_flux(track_id,pe) -= delta;
				pe++;
			    }
		    
			    /* Increment the scalar flux for this thread' copy
			     * of this flat source region */
			    thread_flux(thread_id,fsr_id,e) += fsr_flux;
			}
		    }
		
		    /* Transfer flux to outgoing track */
		    track_out_id = curr_track->_track_out;
		    bc = curr_track->_bc_out;
		    start = curr_track->_refl_out * _polar_times_groups;
		    
		    for (pe=0; pe < _polar_times_groups; pe++) {
		      _boundary_flux(track_out_id,start+pe) = 
			_boundary_flux(track_id,pe) * bc;		  
		      thread_leakage[2*track_id] += _boundary_flux(track_id,pe) 
			* _polar_weights[pe%_num_polar] * (!bc);
		    }		    

		    /* Loop over each segment in reverse direction */
		    for (int s=num_segments-1; s > -1; s--) {
		      
		        curr_segment = &curr_track->_segments[s];
			fsr_id = curr_segment->_region_uid;
			curr_material = &_materials[curr_segment->_material_uid];
			sigma_t = curr_material->_sigma_t;		    
			
			/* Initialize polar angle and energy group counter */
			pe = _polar_times_groups;
		    
			/* Loop over energy groups */
			for (int e=0; e < _num_groups; e++) {
			    fsr_flux = 0.;
			    sigma_t_l = sigma_t[e] * curr_segment->_length;

			    index = int(sigma_t_l*_inverse_prefactor_spacing)*_two_times_num_polar;			
			    /* Loop over polar angles */
			    for (int p=0; p < _num_polar; p++){
			        delta = (_boundary_flux(track_id,pe) - 
					_ratios(fsr_id,e)) * 
				        prefactor(index,p,sigma_t_l);
				fsr_flux += delta * _polar_weights[p];
				_boundary_flux(track_id,pe) -= delta;
				pe++;
			    }
			
			    /* Increment the scalar flux for this thread' copy
			     * of this flat source region */
			    thread_flux(thread_id,fsr_id,e) += fsr_flux;
			}
		    }

		    /* Transfer flux to outgoing track */
		    track_out_id = curr_track->_track_in;
		    bc = curr_track->_bc_in;
		    start = curr_track->_refl_in * _polar_times_groups;
		    
		    for (pe=0; pe < _polar_times_groups; pe++) {
		        _boundary_flux(track_out_id,start+pe) = 
			  _boundary_flux(track_id,_polar_times_groups+pe) * bc;
		      
			thread_leakage[2*track_id+1] += 
			  _boundary_flux(track_id,_polar_times_groups+pe) 
			  * _polar_weights[pe%_num_polar] * (!bc);
		    }
		}       
	    }

	    /** Reduce leakage across threads */
	    leakage = pairwise_sum<FP_PRECISION>(thread_leakage, 
						 2*_tot_num_tracks);

	    /* Reduce scalar fluxes across threads from transport sweep and
	     * add in source term and normalize flux to volume for each region */
	    /* Loop over flat source regions, energy groups */
            #pragma omp parallel for private(volume, sigma_t)
	    for (int r=0; r < _num_FSRs; r++) {

	        volume = _FSR_volumes[r];
		sigma_t = _materials[_FSR_materials[r]]._sigma_t;
	      
		for (int e=0; e < _num_groups; e++) {
		  
		  /* Reduce flux across threads from transport sweep */
		  for (int t=0; t < _num_threads; t++)
		      _scalar_flux(r,e) += thread_flux(t,r,e);
	      
		  _scalar_flux(r,e) *= 0.5;
		  _scalar_flux(r,e) = FOUR_PI * _ratios(r,e) + 
		    (_scalar_flux(r,e) / (sigma_t[e] * volume));

		}
	    }

	    delete [] thread_flux;
	    delete [] thread_leakage;
	}
	      
        _leakage = leakage * 0.5;
	    
	    /* Check for convergence if max_iterations > 1 */
	if (max_iterations == 1 || isScalarFluxConverged())
	    return;
    }
	
	log_printf(WARNING, "Scalar flux did not converge after %d iterations",
		   max_iterations);
	
    return;
}





/**
 * Computes keff on the by performing a series of fixed source
 * iterations and updating the fission and scattering sources in each
 * flat source region of the geometry
 * @param max_iterations the maximum number of iterations allowed
 * @return the value of keff computed
 */
FP_PRECISION MICSolver::convergeSource(int max_iterations) {

    /* Error checking */
    if (_geometry == NULL)
        log_printf(ERROR, "The Solver is unable to converge the source "
		   "since it does not contain a Geometry");
    if (_track_generator == NULL)
        log_printf(ERROR, "The Solver is unable to converge the source "
		   "since it does not contain a TrackGenerator");

    log_printf(NORMAL, "Converging the source...");

    /* Counter for the number of iterations to converge the source */
    _num_iterations = 0;

    /* An initial guess for the eigenvalue */
    _k_eff = 1.0;

    /* The residual on the source */
    FP_PRECISION residual = 0.0;

    /* Initialize data structures */
    initializePolarQuadrature();
    initializeFluxArrays();
    initializeSourceArrays();
    initializePowerArrays();
    precomputePrefactors();
    initializeFSRs();

    /* Check that each FSR has at least one segment crossing it */
    checkTrackSpacing();

    /* Set scalar flux to unity for each region */
    flattenFSRFluxes(1.0);
    flattenFSRSources(1.0);
    zeroTrackFluxes();

    /* Source iteration loop */
    for (int i=0; i < max_iterations; i++) {

        log_printf(NORMAL, "Iteration %d: \tk_eff = %1.6f"
		 "\tres = %1.3E", i, _k_eff, residual);

	normalizeFluxes();
	residual = computeFSRSources();
	transportSweep(1);
	computeKeff();
	_num_iterations++;

	if (i > 1 && residual < _source_convergence_thresh){
//	    transportSweep(1000);
	    return _k_eff;
	}
    }

    log_printf(WARNING, "Unable to converge the source after %d iterations",
	       max_iterations);

    return _k_eff;
}



void MICSolver::computePinPowers() {
}


/**
 * @brief Compute the index into the exponential prefactor hashtable.
 * @details This method computes the index into the exponential prefactor
 *          hashtable for a segment length multiplied by the total 
 *          cross-section of the material the segment resides in.
 * @param sigm_t_l the cross-section multiplied by segment length
 * @return the hasthable index
 */ 
#pragma offload_attribute(push, target(mic))
int MICSolver::computePrefactorIndex(FP_PRECISION sigma_t_l) {
    return int(sigma_t_l * _inverse_prefactor_spacing) * _two_times_num_polar;
}
#pragma offload_attribute(pop)
