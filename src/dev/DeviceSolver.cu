/*
 * DeviceSolver.cu
 *
 *  Created on: Aug 5, 2012
 *      Author: will
 */


#include "DeviceSolver.h"


__constant__ int _num_azim[1];
__constant__ int _num_fsrs[1];
__constant__ int _tot_num_tracks[1];
__constant__ FP_PRECISION _sinthetas[NUM_POLAR_ANGLES];
__constant__ FP_PRECISION _weights[NUM_POLAR_ANGLES];
__constant__ FP_PRECISION _multiples[NUM_POLAR_ANGLES];

/* NOTE: The constant memory for polar weights will only
 * work for up to 256 azimuthal angles! This important since
 * constant memory cannot be allocated at runtime! */
__constant__ FP_PRECISION _polar_weights[NUM_POLAR_ANGLES*128];

#if !STORE_PREFACTORS
__constant__ int _prefactor_max_index[1];
__constant__ FP_PRECISION _prefactor_spacing[1];
#endif


/**
 * DeviceSolver constructor
 * @param geom pointer to the geometry
 * @param track_generator pointer to the TrackGenerator on the CPU
 */
DeviceSolver::DeviceSolver(Geometry* geom, TrackGenerator* track_generator,
										Plotter* plotter, Options* options) {

	_geom = geom;
	_track_generator = track_generator;
	_host_tracks = track_generator->getTracks();
	_plotter = plotter;

	_B = options->getNumThreadBlocks();
	_T = options->getNumThreadsPerBlock();

	_quad = new Quadrature(TABUCHI);

	_FSRs_to_powers = new FP_PRECISION[_geom->getNumFSRs()];
	_FSRs_to_pin_powers = new FP_PRECISION[_geom->getNumFSRs()];

	for (int e = 0; e <= NUM_ENERGY_GROUPS; e++)
		_FSRs_to_fluxes[e] = new FP_PRECISION[_geom->getNumFSRs()];


	/* Allocate memory on the device */
	allocateDeviceMemory();
}



/**
 * Solver destructor frees all memory on the device
 */
DeviceSolver::~DeviceSolver() {

    log_printf(NORMAL, "Cleaning up memory on the device...");

    delete [] _FSRs_to_powers;
	delete [] _FSRs_to_pin_powers;

	for (int e = 0; e <= NUM_ENERGY_GROUPS; e++)
		delete [] _FSRs_to_fluxes[e];

	/* Free FSRs and materials on device */
	CUDA_SAFE_CALL(cudaFree(_FSRs));
	CUDA_SAFE_CALL(cudaFree(_materials));
	CUDA_SAFE_CALL(cudaFree(_dev_tracks));
	CUDA_SAFE_CALL(cudaFree(_track_index_offsets));

#if !STORE_PREFACTORS
	/* Free the prefactor array on the device */
	CUDA_SAFE_CALL(cudaFree(_prefactor_array));
#endif

	CUDA_SAFE_CALL(cudaFree(_num_tracks));

	_fission_source_vec.clear();
	CUDA_SAFE_CALL(cudaFreeHost(_renorm_factor));
	CUDA_SAFE_CALL(cudaFreeHost(_k_eff));

	_tot_abs_vec.clear();
	_tot_fission_vec.clear();
	_source_residual_norm_vec.clear();

	return;
}


/**
 * Allocates and initializes all memory on the device necessary for fixed
 * source iteration, including tracks, segments, flat source regions
 * materials, and the polar quadrature
 */
void DeviceSolver::allocateDeviceMemory() {

    log_printf(NORMAL, "Allocating memory on the device.");

	/* Copy variables for number of angles, tracks, and FSRs */
	int num_fsrs = _geom->getNumFSRs();
	int num_azim = _track_generator->getNumAzim();

    int tot_num_tracks = 0;
    int tot_num_segments = 0;
	for (int i=0; i < num_azim; i++)
		tot_num_tracks += _track_generator->getNumTracks()[i];

	/* Pinned host memory */
	unsigned int flags = cudaHostAllocWriteCombined;
	CUDA_SAFE_CALL(cudaHostAlloc((void**)&_renorm_factor,
								sizeof(FP_PRECISION), flags));
	CUDA_SAFE_CALL(cudaHostAlloc((void**)&_k_eff,
										sizeof(FP_PRECISION), flags));

	/* Allocate Flat source regions array */
	CUDA_SAFE_CALL(cudaMalloc((void**)&_FSRs,
						_geom->getNumFSRs() * sizeof(dev_flatsourceregion)));

	/* Allocate materials array */
	CUDA_SAFE_CALL(cudaMalloc((void**)&_materials,
						_geom->getNumMaterials() * sizeof(dev_material)));

	/* Allocate tracks array */
	CUDA_SAFE_CALL(cudaMalloc((void**)&_dev_tracks,
				  tot_num_tracks*sizeof(dev_track)));

	/* Allocate variables for the number of angles, tracks and FSRs */
	CUDA_SAFE_CALL(cudaMalloc((void**)&_num_tracks, num_azim*sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(_num_azim, (void*)&num_azim,
								sizeof(int), 0, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(_tot_num_tracks, &tot_num_tracks,
								sizeof(int), 0, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMalloc((void**)&_track_index_offsets,
												(num_azim+1)*sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(_num_fsrs, &num_fsrs,
								sizeof(int), 0, cudaMemcpyHostToDevice));


	/* Allocate fission source array on device */
	_fission_source_vec.resize(_B*_T);
	_fission_source = thrust::raw_pointer_cast(&_fission_source_vec[0]);

	/* Allocate residual norm array to determine convergence of sources */
	_source_residual_norm_vec.resize(_B*_T);
	_source_residual_norm =
			thrust::raw_pointer_cast(&_source_residual_norm_vec[0]);


	/* Allocate array for total absorption and fission */
	_tot_abs_vec.resize(_B*_T);
	_tot_fission_vec.resize(_B*_T);
	_tot_abs = thrust::raw_pointer_cast(&_tot_abs_vec[0]);
	_tot_fission = thrust::raw_pointer_cast(&_tot_fission_vec[0]);

	/* Allocate memory on device for FSRs */
    dev_flatsourceregion* temp_fsrs = (dev_flatsourceregion*)
    		malloc(num_fsrs * sizeof(dev_flatsourceregion));
	initializeFSRs(temp_fsrs);
	CUDA_SAFE_CALL(cudaMemcpy((void*)_FSRs, (void*)temp_fsrs,
							num_fsrs*sizeof(dev_flatsourceregion),
							cudaMemcpyHostToDevice));
	free(temp_fsrs);

	/* Iterate through all materials and clone them on the device */
	int index = 0;
	std::map<short int, Material*> host_materials = _geom->getMaterials();
	std::map<short int, Material*>::iterator iter;
	for (iter = host_materials.begin(); iter != host_materials.end(); ++iter) {
		cloneOnDevice(iter->second, &_materials[(*iter).second->getUid()]);
		index++;
	}

	/* Iterate through all tracks and clone them on the device */
	int counter = 0;
	for (int i=0; i < num_azim; i++) {

		/* Copy the cumulative index offset for the current azimuthal angle */
		CUDA_SAFE_CALL(cudaMemcpy((void*)&_track_index_offsets[i],
					(void*)&counter, sizeof(int), cudaMemcpyHostToDevice));

		for (int j=0; j < _track_generator->getNumTracks()[i]; j++) {
			cloneTrack(&_host_tracks[i][j], &_dev_tracks[counter]);

			/* Make tracks reflective */
			index = computeScalarTrackIndex(_host_tracks[i][j].getTrackInI(),
											_host_tracks[i][j].getTrackInJ());
			CUDA_SAFE_CALL(cudaMemcpy((void*)&_dev_tracks[counter]._track_in,
													(void*)&index, sizeof(int),
													cudaMemcpyHostToDevice));

			index = computeScalarTrackIndex(_host_tracks[i][j].getTrackOutI(),
	    									_host_tracks[i][j].getTrackOutJ());
			CUDA_SAFE_CALL(cudaMemcpy((void*)&_dev_tracks[counter]._track_out,
													(void*)&index, sizeof(int),
													cudaMemcpyHostToDevice));

			tot_num_segments += _host_tracks[i][j].getNumSegments();
			counter++;
		}
	}

	log_printf(NORMAL, "Total number of tracks = %d, segments = %d",
									tot_num_tracks, tot_num_segments);


	/* Compute polar/azim angle weights */
	FP_PRECISION* polar_weights =
		(FP_PRECISION*)malloc(NUM_POLAR_ANGLES*num_azim*sizeof(FP_PRECISION));
	FP_PRECISION* azim_weights = _track_generator->getAzimWeights();

	for (int i=0; i < num_azim; i++) {
		for (int j=0; j < NUM_POLAR_ANGLES; j++)
			polar_weights[i*NUM_POLAR_ANGLES + j] = azim_weights[i] *
								_quad->getMultiples()[j] * FOUR_PI;
	}

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(_polar_weights, (void*)polar_weights,
							NUM_POLAR_ANGLES*num_azim*sizeof(FP_PRECISION),
							0, cudaMemcpyHostToDevice));
	free(polar_weights);

   /* Copy the cumulative index offset for the current azimuthal angle */
	CUDA_SAFE_CALL(cudaMemcpy((void*)&_track_index_offsets[num_azim],
										(void*)&counter, sizeof(int),
										cudaMemcpyHostToDevice));

	CUDA_SAFE_CALL(cudaMemcpy((void*)_num_tracks,
					(void*)_track_generator->getNumTracks(),
					num_azim*sizeof(int), cudaMemcpyHostToDevice));

	/* Copy polar quadrature to the device */
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(_sinthetas, (void*)_quad->getSinThetas(),
			NUM_POLAR_ANGLES*sizeof(FP_PRECISION), 0, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(_weights, (void*)_quad->getWeights(),
			NUM_POLAR_ANGLES*sizeof(FP_PRECISION), 0, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(_multiples, (void*)_quad->getMultiples(),
			NUM_POLAR_ANGLES*sizeof(FP_PRECISION), 0, cudaMemcpyHostToDevice));

	computePrefactors();

	return;
}


/**
 * Copies memory from the device and validates it with the host version of
 * each object, including tracks, segments, materials, quadrature, and flat
 * source regions. If any attribute is of any object is not identical or
 * within a threshold of the corresponding value on the device, this method
 * prints an error message and exits the program
 */
void DeviceSolver::validateDeviceDataIntegrity() {

	log_printf(NORMAL, "Validating device data...");

	/* Validate materials */
	int num_materials = _geom->getMaterials().size();
	dev_material dev_material;
	Material* host_material = NULL;

	for (int i=0; i < num_materials; i++) {

		/* Copy material back from the device */
		CUDA_SAFE_CALL(cudaMemcpy((void*)&dev_material, (void*)&_materials[i],
								sizeof(dev_material), cudaMemcpyDeviceToHost));

		/* Find the host material with the same uid */
		std::map<short int, Material*>::iterator iter;
		for (iter = _geom->getMaterials().begin();
								iter != _geom->getMaterials().end(); ++iter) {

			if (dev_material._uid == (*iter).second->getUid()) {
				host_material = (*iter).second;
				break;
			}

		}

		if (host_material == NULL)
		  log_printf(ERROR, "Could not find a material with uid = %d on the"
				  	  	  	  	  	  	  	  	  "host", dev_material._uid);

		/* Verify all cross sections */
		for (int p=0; p < NUM_ENERGY_GROUPS; p++) {

			/* Validate sigma_t */
			if (dev_material._sigma_t[p] != host_material->getSigmaT()[p])
				log_printf(ERROR, "Could not validate sigma_t on device for "
						"material with uid = %d and energy group %d: device "
						"sigma_t = %f, host sigma_t = %f",
						dev_material._uid, p, dev_material._sigma_t[p],
		       	   	   	host_material->getSigmaT()[p]);

			/* Validate sigma_a */
			if (dev_material._sigma_a[p] != host_material->getSigmaA()[p])
				log_printf(ERROR, "Could not validate sigma_a on device for "
						"material with uid = %d and energy group %d: device "
						"sigma_a = %f, host sigma_a = %f",
						dev_material._uid, p, dev_material._sigma_a[p],
						host_material->getSigmaA()[p]);

			/* validate sigma_f */
			if (dev_material._sigma_f[p] != host_material->getSigmaF()[p])
				log_printf(ERROR, "Could not validate sigma_f on device for "
						"material with uid = %d and energy group %d: device "
						"sigma_f = %f, host sigma_f = %f",
						dev_material._uid, p, dev_material._sigma_f[p],
						host_material->getSigmaF()[p]);

			/* Validate nu_sigma_f */
			if (dev_material._nu_sigma_f[p] != host_material->getNuSigmaF()[p])
				log_printf(ERROR, "Could not validate sigma_f on device for "
						"material with uid = %d and energy group %d: device "
						"nu_sigma_t = %f, host nu_sigma_f = %f",
						dev_material._uid, p, dev_material._nu_sigma_f[p],
						host_material->getNuSigmaF()[p]);

			/* Validate chi */
			if (dev_material._chi[p] != host_material->getChi()[p])
				log_printf(ERROR, "Could not validate chi on device for "
						"material with uid = %d and energy group %d: device "
						"chi = %f, host chi = %f",
						dev_material._uid, p, dev_material._chi[p],
						host_material->getChi()[p]);

			/* Validate sigma_s */
			for (int q=0; q < NUM_ENERGY_GROUPS; q++) {
				if (dev_material._sigma_s[p][q] !=
							host_material->getSigmaS()[p*NUM_ENERGY_GROUPS + q])
				  	  log_printf(ERROR, "Could not validate sigma_s on device"
			  			  " for material with uid = %d and energy groups "
			  			  " %d, %d: device sigma_s = %f, host sigma_s = %f",
						  dev_material._uid, p, q, dev_material._sigma_s[p][q],
						  host_material->getSigmaS()[p*NUM_ENERGY_GROUPS + q]);
			}
		}

		/* Reset host material to NULL for next iteration */
		host_material = NULL;
	}


	/* Validate tracks and segments */
    dev_track test_track;
    dev_segment* test_segments;
    int counter = 0;
    int track_in, track_out;
    int azim_angle_index;
    int num_azim = _track_generator->getNumAzim();

    FP_PRECISION* polar_weights =
   		(FP_PRECISION*)malloc(num_azim*NUM_POLAR_ANGLES*sizeof(FP_PRECISION));

    CUDA_SAFE_CALL(cudaMemcpyFromSymbol((void*)polar_weights, _polar_weights,
    				num_azim*NUM_POLAR_ANGLES*sizeof(FP_PRECISION),
    										0, cudaMemcpyDeviceToHost));

    /* Iterate over all tracks */
    for (int i=0; i < num_azim; i++) {
    	for (int j=0; j < _track_generator->getNumTracks()[i]; j++) {

    		/* Copy track back from the device */
    		CUDA_SAFE_CALL(cudaMemcpy((void*)&test_track,
			  	  	  (void*)&_dev_tracks[counter],
			  	  	  sizeof(dev_track), cudaMemcpyDeviceToHost));

    		azim_angle_index = test_track._azim_angle_index;

    		/* Validate azimuthal angle index */
    		if (i != azim_angle_index)
    			log_printf(WARNING, "Could not validate the azimuthal angle "
    					"index for device track: host index = %d, device index"
    					"= %d", i, azim_angle_index);

    		/* Validate polar fluxes */
    		for (int p=0; p < 2*GRP_TIMES_ANG; p++) {
    			if (fabs(test_track._polar_fluxes[p] -
    						_host_tracks[i][j].getPolarFluxes()[p]) > 1E-5)
    				log_printf(WARNING, "Could not validate polar flux for "
					  "device track i = %d j = %d for angle %d: device polar "
					  "flux = %f, host polar flux = %f", i, j, p,
					  test_track._polar_fluxes[p],
					  _host_tracks[i][j].getPolarFluxes()[p]);
    		}


    		/* Validate the number of segments */
    		if (test_track._num_segments != _host_tracks[i][j].getNumSegments())
    			log_printf(WARNING, "Could not validate the number of segments"
    				" for device track i = %d, j = %d: device num segments "
    				" = %d, host num segments = %d", i, j,
    				test_track._num_segments,
    				_host_tracks[i][j].getNumSegments());

    		/* Validate the incoming and outgoing tracks for boundary conditions */
    		track_in = computeScalarTrackIndex(_host_tracks[i][j].getTrackInI(),
											  _host_tracks[i][j].getTrackInJ());
    		track_out = computeScalarTrackIndex(_host_tracks[i][j].getTrackOutI(),
											 _host_tracks[i][j].getTrackOutJ());

    		if (test_track._track_in != track_in)
    			log_printf(WARNING, "Could not validate track in for device "
    					"track i = %d, j = %d: device track in = %d, host "
    					"track in = %d", i, j, test_track._track_in, track_in);

    		if (test_track._track_out != track_out)
    			log_printf(WARNING, "Could not validate track out for device "
    					"track i = %d, j = %d: device track out = %d, host "
    					"track out = %d", i, j, test_track._track_out,
    																track_out);

    		if (test_track._refl_in  != _host_tracks[i][j].isReflIn())
    			log_printf(WARNING, "Could not validate refl_in for device "
    					"track i = %d, j = %d: device track refl_in = %d, "
    					"host track refl_in = %d", i, j, test_track._refl_in,
    					_host_tracks[i][j].isReflIn());

    		if (test_track._refl_out  !=  _host_tracks[i][j].isReflOut())
    			log_printf(WARNING, "Could not validate refl_out for device "
    					"track i = %d, j = %d: device track refl_out = %d, "
    					"host track refl_out = %d", i, j, test_track._refl_out,
    					_host_tracks[i][j].isReflOut());

    		test_segments = new dev_segment[test_track._num_segments];

    		/* Copy back the segments from the device */
    		CUDA_SAFE_CALL(cudaMemcpy((void*)test_segments,
    									(void*)test_track._segments,
    		test_track._num_segments*sizeof(dev_segment),
    										cudaMemcpyDeviceToHost));

    		/* Validate all of the track's segments */
    		test_segments = new dev_segment[test_track._num_segments];
			/* Copy back the segments from the device */
			CUDA_SAFE_CALL(cudaMemcpy((void*)test_segments,
						(void*)test_track._segments,
						test_track._num_segments*sizeof(dev_segment),
						cudaMemcpyDeviceToHost));

    		for (int p=0; p < test_track._num_segments; p++) {
    			/* Validate segment length */
    			if (test_segments[p]._length !=
    								_host_tracks[i][j].getSegment(p)->_length)

    				log_printf(WARNING, "Could not validate the segment "
    						"length for device track i = %d, j = %d, segment "
    						"# = %d: device seg length = %f, host seg length "
    						"= %f", i, j, p, test_segments[p]._length,
    						_host_tracks[i][j].getSegment(p)->_length);

    			/* Validate segment material uid */
    			if (test_segments[p]._material_uid !=
						_host_tracks[i][j].getSegment(p)->_material->getUid())

    				log_printf(WARNING, "Could not validate the segment "
						"material uid for device track i = %d, j = %d, "
						"segment # = %d: device material uid = %d, host "
						"seg material uid = %d", i, j, p,
						test_segments[p]._material_uid,
   						_host_tracks[i][j].getSegment(p)->_material->getUid());

    			/* Validate segment flat source region uid */
    			if (test_segments[p]._region_uid !=
   								_host_tracks[i][j].getSegment(p)->_region_id)

				   log_printf(WARNING, "Could not validate the segment "
						   "region uid for device track i = %d, j = %d, "
						   "segment # = %d: device seg region uid = %d, "
						   "host seg region uid = %d", i, j, p,
						   test_segments[p]._region_uid,
						   _host_tracks[i][j].getSegment(p)->_region_id);
    		}

    		free(test_segments);

    		/* increment the counter as an index for the device tracks */
    		counter++;
    	}
    }


	/* Validate variables for the number of angles , tracks, and FSRs */
	int* num_tracks;
	int tot_num_tracks;
	int num_fsrs;

	CUDA_SAFE_CALL(cudaMemcpyFromSymbol((void*)&num_azim, _num_azim,
							  sizeof(int), 0, cudaMemcpyDeviceToHost));

	num_tracks = (int*)malloc(num_azim*sizeof(int));


	CUDA_SAFE_CALL(cudaMemcpy((void*)num_tracks, (void*)_num_tracks,
						  num_azim*sizeof(int), cudaMemcpyDeviceToHost));

	CUDA_SAFE_CALL(cudaMemcpyFromSymbol((void*)&tot_num_tracks,
					_tot_num_tracks, sizeof(int), 0, cudaMemcpyDeviceToHost));

	CUDA_SAFE_CALL(cudaMemcpyFromSymbol((void*)&num_fsrs, _num_fsrs,
							sizeof(int), 0, cudaMemcpyDeviceToHost));

	/* Validate number of azimuthal angles */
	if (num_azim != _track_generator->getNumAzim())
		log_printf(ERROR, "Could not validate num_azim on device: device "
			"num_azim = %d, host num_azim = %d", num_azim,
			_track_generator->getNumAzim());

	/* Validate the number of tracks per azimuthal angle */
	for (int i=0; i < num_azim; i++) {

		if (num_tracks[i] != _track_generator->getNumTracks()[i])
			log_printf(ERROR, "Could not validate num_tracks for azimuthal "
					"angle %d: device num_tracks = %d, host num_tracks = %d",
					i, num_tracks[i], _track_generator->getNumTracks()[i]);
	}


	/* Validate the polar quadrature on device */
	FP_PRECISION sinthetas[NUM_POLAR_ANGLES];
	FP_PRECISION multiples[NUM_POLAR_ANGLES];
	FP_PRECISION weights[NUM_POLAR_ANGLES];

	CUDA_SAFE_CALL(cudaMemcpyFromSymbol((void*)sinthetas, _sinthetas,
			NUM_POLAR_ANGLES*sizeof(FP_PRECISION), 0, cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpyFromSymbol((void*)weights, _weights,
			NUM_POLAR_ANGLES*sizeof(FP_PRECISION), 0, cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpyFromSymbol((void*)multiples, _multiples,
			NUM_POLAR_ANGLES*sizeof(FP_PRECISION), 0, cudaMemcpyDeviceToHost));

	for (int i=0; i < NUM_POLAR_ANGLES; i++) {

		if (sinthetas[i] != _quad->getSinThetas()[i])
			  log_printf(ERROR, "Could not validate polar sin theta on device "
				  "for angle i = %d: device sin theta = %f, host sin "
				  "theta = %f", i, sinthetas[i], _quad->getSinThetas()[i]);

		if (multiples[i] != _quad->getMultiples()[i])

			log_printf(ERROR, "Could not validate polar multiple on device "
					  "for angle i = %d: device multiple = %f, host multiple "
					  "= %f", i, multiples[i], _quad->getMultiples()[i]);

		if (weights[i] != _quad->getWeights()[i])

			log_printf(ERROR, "Could not validate polar weight on device "
					  "for angle i = %d: device weight = %f, host weight "
					  "= %f", i, weights[i], _quad->getWeights()[i]);
	}

	return;
}


/**
 * This is a helper method for the allocateDeviceMemory method. It
 * initializes an array of dev_flatsourceregion structs on the host
 * with the appropriate values (volume, region uid, and material uid)
 * so that allocateDeviceMemory can copy the array in full to the device.
 * @param dev_FSRs an array of dev_flatsourceregions on the device
 */
void DeviceSolver::initializeFSRs(dev_flatsourceregion* dev_FSRs) {

	log_printf(NORMAL, "Initializing FSRs...");

	CellBasic* cell;
	Material* material;
	Universe* univ_zero = _geom->getUniverse(0);
	Track* track;
	segment* seg;
	dev_flatsourceregion* fsr;

	/* Initialize each FSRs volume to 0 to avoid NaNs */
	for (int r = 0; r < _geom->getNumFSRs(); r++)
	  dev_FSRs[r]._volume = 0.0;


	/* Set each FSR's volume by accumulating the total length of all
	   tracks inside the FSR. Iterate over azimuthal angle, track and segment */
	for (int i = 0; i < _track_generator->getNumAzim(); i++) {
		for (int j = 0; j < _track_generator->getNumTracks()[i]; j++) {

			track = &_track_generator->getTracks()[i][j];

			/* Iterate over all of the track's segments to update FSR volumes */
			for (int s = 0; s < track->getNumSegments(); s++) {
				seg = track->getSegment(s);
				fsr = &dev_FSRs[seg->_region_id];
				fsr->_volume += seg->_length * track->getAzimuthalWeight();
			}
		}
	}

	/* Iterate over all FSRs */
	for (int r = 0; r < _geom->getNumFSRs(); r++) {

		/* Set the uid for this FSR */
		dev_FSRs[r]._uid = r;

		/* Get the cell corresponding to this FSR from the geometry */
		cell = static_cast<CellBasic*>(_geom->findCell(univ_zero, r));

		/* Get the cell's material and assign it to the FSR */
		material = _geom->getMaterial(cell->getMaterial());
		dev_FSRs[r]._material_uid = material->getUid();
	}

	return;
}


/**
 * This method computes the index for the jth track at azimuthal angle i.
 * This method is necessary since the array of tracks on the device is a 1D
 * array which needs a one-to-one mapping from the 2D jagged array of tracks
 * on the host
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
		index += _track_generator->getNumTracks()[p];
		p++;
	}

	/* Update index for this track since it is the jth track at angle i */
	index += j;

	return index;
}



/**
 * This method computes the exponential prefactors from the transport
 * equation and stores them as part of each track segment or in a
 * table depending on STORE_PREFACTORS is set to true or false,
 * respectively, in configurations.h
 */
void DeviceSolver::computePrefactors(){

#if STORE_PREFACTORS
	computePrefactorsOnDevice<<<_B, _T>>>(_dev_tracks, _materials);
#else

	/* set size of prefactor array */
	int num_array_values = 10 * sqrt(1 / (8 * SOURCE_CONVERG_THRESH));
	FP_PRECISION prefactor_spacing = 10.0 / num_array_values;
	int prefactor_array_size = 2 * NUM_POLAR_ANGLES * num_array_values;
	int prefactor_max_index = prefactor_array_size - 2*NUM_POLAR_ANGLES - 1;

	/* allocate arrays */
	FP_PRECISION* prefactor_array = new FP_PRECISION[prefactor_array_size];

	FP_PRECISION expon;
	FP_PRECISION intercept;
	FP_PRECISION slope;


	/* Create prefactor array */
	for (int i = 0; i < num_array_values; i ++){
		for (int j = 0; j < NUM_POLAR_ANGLES; j++){
			expon = exp(- (i * prefactor_spacing) / _quad->getSinTheta(j));
			slope = - expon / _quad->getSinTheta(j);
			intercept = expon * (1 + (i * prefactor_spacing) /
													_quad->getSinTheta(j));
			prefactor_array[2 * NUM_POLAR_ANGLES*i + 2*j] = slope;
			prefactor_array[2 * NUM_POLAR_ANGLES*i + 2*j + 1] = intercept;
		}
	}

	CUDA_SAFE_CALL(cudaMalloc((void**)&_prefactor_array,
					prefactor_array_size*sizeof(FP_PRECISION)));

	CUDA_SAFE_CALL(cudaMemcpy((void*)_prefactor_array, (void*)prefactor_array,
									prefactor_array_size*sizeof(FP_PRECISION),
													cudaMemcpyHostToDevice));

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(_prefactor_spacing,
								(void*)&prefactor_spacing,
								sizeof(int), 0, cudaMemcpyHostToDevice));

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(_prefactor_max_index,
								(void*)&prefactor_max_index,
								sizeof(int), 0, cudaMemcpyHostToDevice));

	free(prefactor_array);

#endif
	return;
}


/**
 * This method performs on or more fixed source iterations on the
 * device by integrating the flux along each track and updating the
 * boundary fluxes for the corresponding output track, while updating
 * the scalar flux in each flat source region
 * @param max_iterations the maximum number of iterations allowed
 */
void DeviceSolver::fixedSourceIteration(int max_iterations) {

	int shared_mem_size = sizeof(FP_PRECISION)*_T*(2*NUM_POLAR_ANGLES + 1);

	for (int i=0; i < max_iterations; i++) {

		/* Initialize the flux in each region to zero */
		zeroFSRFluxesOnDevice<<<_B, _T>>>(_FSRs);

#if STORE_PREFACTORS
		/* First angle for each pair of tracks */
		propagateTrackFluxForwardOnDevice<<<_B, _T, shared_mem_size >>>
									(_dev_tracks, _track_index_offsets, _FSRs);

		/* First angle for each pair of tracks */
		propagateTrackFluxReverseOnDevice<<<_B, _T, shared_mem_size>>>
									(_dev_tracks, _track_index_offsets, _FSRs);

#else
		/* First angle for each pair of tracks */
		propagateTrackFluxForwardOnDevice<<<_B, _T, shared_mem_size>>>
									(_dev_tracks, _track_index_offsets,
										_materials, _FSRs, _prefactor_array);

		/* First angle for each pair of tracks */
		propagateTrackFluxReverseOnDevice<<<_B, _T, shared_mem_size>>>
									(_dev_tracks, _track_index_offsets,
										_materials, _FSRs, _prefactor_array);
#endif

		/* Add in source term, normalize fluxes to volume and save old flux */
		normalizeFluxToVolumeOnDevice<<<_B, _T>>>(_FSRs, _materials);
	}
}


/**
 * Computes keff on the device by performing a series of fixed source
 * iterations and updating the fission and scattering sources in each
 * flat source region of the geometry
 * @param max_iterations the maximum number of iterations allowed
 * @return the value of keff computed
 */
FP_PRECISION DeviceSolver::computeKeff(int max_iterations) {


	FP_PRECISION tot_abs, tot_fission;
	FP_PRECISION tot_source_residual = 0.0;
	FP_PRECISION* renorm_factor;
	FP_PRECISION* k_eff;

	FILE *keff_file;
	keff_file = fopen("keff_file.txt", "w");

	CUDA_SAFE_CALL(cudaHostGetDevicePointer((void**)&renorm_factor,
											(void*)_renorm_factor, 0));
	CUDA_SAFE_CALL(cudaHostGetDevicePointer((void**)&k_eff,
												(void*)_k_eff, 0));

	log_printf(NORMAL, "Computing k_eff on device...");

	/* Make initial guess for keff */
	*_k_eff = 1.0;

	//NOTE: Assume here that checkTrackSpacing was already called on host
	//to check that each FSR has at least one track in it

	/* Set scalar flux and source to unity for each region (initial guess) */
	initDeviceFSRsOnDevice<<<_B, _T>>>(_FSRs);

	/* Set the incoming and outgoing fluxes for each track to zero */
	zeroTrackFluxesOnDevice<<<_B, _T>>>(_dev_tracks);

	/* Source iteration loop */
	for (int i=0; i < max_iterations; i++) {

		log_printf(NORMAL, "Iteration %d on device: \t\tk_eff = %1.6f"
							"\t\tres = %1.3E", i, *_k_eff, tot_source_residual);

		/*********************************************************************
		 * ReDEBUGize scalar and boundary fluxes
		 *********************************************************************/

		computeTotalFissionSourceOnDevice<<<_B, _T, sizeof(FP_PRECISION)*_T>>>
										(_FSRs, _materials, _fission_source);

		*_renorm_factor = 1.0 / thrust::reduce(_fission_source_vec.begin(),
													_fission_source_vec.end());

		log_printf(DEBUG, "renorm_factor = %f", *_renorm_factor);

		normalizeTrackFluxesOnDevice<<<_B, _T>>>(_dev_tracks, renorm_factor);
		normalizeFSRFluxesOnDevice<<<_B, _T>>>(_FSRs, renorm_factor);

		/*********************************************************************
		 * Compute the source and ratios for each region
		 *********************************************************************/

		computeFSRSourcesOnDevice<<<_B, _T, sizeof(FP_PRECISION)*_T*2>>>(_FSRs,
									_materials, k_eff, _source_residual_norm);

		/*********************************************************************
		 * Update flux and check for convergence
		 *********************************************************************/

		/* Iteration the flux with the new source */
		fixedSourceIteration(1);

		/* Update k_eff */
		computeAbsAndFissionOnDevice<<<_B, _T, sizeof(FP_PRECISION)*_T*2>>>
								(_FSRs, _materials, _tot_abs, _tot_fission);

		tot_abs = thrust::reduce(_tot_abs_vec.begin(), _tot_abs_vec.end());
		tot_fission = thrust::reduce(_tot_fission_vec.begin(),
													_tot_fission_vec.end());

		/* Compute the residual of the norm of the new and old sources */
		tot_source_residual = thrust::reduce(_source_residual_norm_vec.begin(),
											_source_residual_norm_vec.end());;
		tot_source_residual = sqrt(tot_source_residual / _geom->getNumFSRs());

		log_printf(DEBUG, "tot_abs = %f, tot_fission = %f",
														tot_abs, tot_fission);

		*_k_eff = tot_fission / tot_abs;

		fprintf(keff_file, "%1.6f\n", *_k_eff);

		/* If keff converged, return keff */
		if (i > 1 && tot_source_residual < SOURCE_CONVERG_THRESH) {

			/* Plot the fluxes if the user requested it */
			if (_plotter->plotFlux() == true){

				/* Copy fluxes back from FSRs on the device */
				int num_fsrs = _geom->getNumFSRs();

				/* Allocate memory on device for FSRs */
			    dev_flatsourceregion* temp_fsrs = (dev_flatsourceregion*)
							malloc(num_fsrs * sizeof(dev_flatsourceregion));

				CUDA_SAFE_CALL(cudaMemcpy((void*)temp_fsrs, (void*)_FSRs,
										num_fsrs * sizeof(dev_flatsourceregion),
										cudaMemcpyDeviceToHost));

				/* Load fluxes into FSR to flux map */
				for (int r=0; r < _geom->getNumFSRs(); r++) {

					FP_PRECISION* fluxes = temp_fsrs[r]._flux;

					for (int e=0; e < NUM_ENERGY_GROUPS; e++){

						_FSRs_to_fluxes[e][r] = fluxes[e];
						_FSRs_to_fluxes[NUM_ENERGY_GROUPS][r] =
							_FSRs_to_fluxes[NUM_ENERGY_GROUPS][r] + fluxes[e];
					}
				}

				plotFluxes();

				free(temp_fsrs);
			}

			fclose(keff_file);
			return *_k_eff;
		}

		cudaDeviceSynchronize();

	}

	log_printf(WARNING, "Unable to converge the source after %d iterations",
															max_iterations);
	return 0;
}


/**
 * Plots the fluxes for each flat source region using the
 * Plotter class and quickplot utility wrapper functions for
 * ImageMagick
 */
void DeviceSolver::plotFluxes(){

	/* create BitMaps for plotting */
	BitMap<int>* bitMapFSR = new BitMap<int>;
	BitMap<FP_PRECISION>* bitMap = new BitMap<FP_PRECISION>;
	bitMapFSR->pixel_x = _plotter->getBitLengthX();
	bitMapFSR->pixel_y = _plotter->getBitLengthX();
	bitMap->pixel_x = _plotter->getBitLengthX();
	bitMap->pixel_y = _plotter->getBitLengthY();
	initialize(bitMapFSR);
	initialize(bitMap);
	bitMap->geom_x = _geom->getWidth();
	bitMap->geom_y = _geom->getHeight();
	bitMapFSR->color_type = RANDOM;
	bitMap->color_type = SCALED;

	/* make FSR BitMap */
	_plotter->makeFSRMap(bitMapFSR->pixels);

	for (int i = 0; i < NUM_ENERGY_GROUPS; i++){

		std::stringstream string;
		string << "flux" << i + 1 << "group";
		std::string title_str = string.str();

		log_printf(NORMAL, "Plotting group %d flux...", (i+1));
		_plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels,
												_FSRs_to_fluxes[i]);
		plot(bitMap, title_str, _plotter->getExtension());
	}

	log_printf(NORMAL, "Plotting total flux...");
	_plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels,
							_FSRs_to_fluxes[NUM_ENERGY_GROUPS]);
	plot(bitMap, "flux_total", _plotter->getExtension());

	/* delete bitMaps */
	deleteBitMap(bitMapFSR);
	deleteBitMap(bitMap);

	return;
}


/**
 * Compute the fission rates in each FSR and save them in a map of
 * FSR ids to fission rates. Plots the pin powers using the Plotter
 * class and quickplot wrapper utility functions for ImageMagick
 */
void DeviceSolver::computePinPowers() {

	log_printf(NORMAL, "Computing pin powers...");

	dev_flatsourceregion* fsr;
	FP_PRECISION tot_pin_power = 0;
	FP_PRECISION avg_pin_power = 0;
	FP_PRECISION num_nonzero_pins = 0;
	FP_PRECISION curr_pin_power = 0;
	FP_PRECISION prev_pin_power = 0;

	/* create BitMaps for plotting */
	BitMap<int>* bitMapFSR = new BitMap<int>;
	BitMap<FP_PRECISION>* bitMap = new BitMap<FP_PRECISION>;
	bitMapFSR->pixel_x = _plotter->getBitLengthX();
	bitMapFSR->pixel_y = _plotter->getBitLengthY();
	bitMap->pixel_x = _plotter->getBitLengthX();
	bitMap->pixel_y = _plotter->getBitLengthY();
	initialize(bitMapFSR);
	initialize(bitMap);
	bitMap->geom_x = _geom->getWidth();
	bitMap->geom_y = _geom->getHeight();
	bitMapFSR->color_type = RANDOM;
	bitMap->color_type = SCALED;

	/* Copy fluxes back from FSRs on the device */
	int num_fsrs = _geom->getNumFSRs();

	/* Allocate memory on device for FSRs */
    dev_flatsourceregion* temp_fsrs = (dev_flatsourceregion*)
				malloc(num_fsrs * sizeof(dev_flatsourceregion));

	CUDA_SAFE_CALL(cudaMemcpy((void*)temp_fsrs, (void*)_FSRs,
							num_fsrs * sizeof(dev_flatsourceregion),
							cudaMemcpyDeviceToHost));

	/* make FSR BitMap */
	_plotter->makeFSRMap(bitMapFSR->pixels);

	/* Loop over all FSRs and compute the fission rate*/
	for (int i=0; i < _geom->getNumFSRs(); i++) {
		fsr = &temp_fsrs[i];
		int material_id = _geom->getFSRtoMaterialMap()[fsr->_uid];
		Material* material = _geom->getMaterial(material_id);
		_FSRs_to_powers[i] = computeFissionRate(fsr, material);
	}

	/* Compute the pin powers by adding up the powers of FSRs in each
	 * lattice cell, saving lattice cell powers to files, and saving the
	 * pin power corresponding to each FSR id in FSR_to_pin_powers */
	_geom->computePinPowers(_FSRs_to_powers, _FSRs_to_pin_powers);


	/* Compute the total power based by accumulating the power of each unique
	 * pin with a nonzero power */
	for (int i=0; i < _geom->getNumFSRs(); i++) {
		curr_pin_power = _FSRs_to_pin_powers[i];

		/* If this pin power is unique and nozero (doesn't match the previous
		 * pin's power), then tally it
		 */
		if (curr_pin_power > 0 && curr_pin_power != prev_pin_power) {
			tot_pin_power += curr_pin_power;
			num_nonzero_pins++;	int

			prev_pin_power = curr_pin_power;
		}
	}

	/* Compute the average pin power */
	avg_pin_power = tot_pin_power / num_nonzero_pins;

	/* Normalize each pin power to the average non-zero pin power */
	for (int i=0; i < _geom->getNumFSRs(); i++) {
		_FSRs_to_pin_powers[i] /= avg_pin_power;
	}


	log_printf(NORMAL, "Plotting pin powers...");
	_plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels,
											_FSRs_to_pin_powers);
	plot(bitMap, "pin_powers", _plotter->getExtension());

	/* delete bitMaps */
	deleteBitMap(bitMapFSR);
	deleteBitMap(bitMap);

	free(temp_fsrs);

	return;
}


/**
 * Computes exponential prefactors for each segment on the device
 * @param tracks pointer to the array of device tracks
 * @param tot_num_tracks pointer to an array of the number of tracks per angle
 * @param materials pointer to an array of device materials
 * @param sinthetas pointer to an array of sintheta values
 */
#if STORE_PREFACTORS
__global__ void computePrefactorsOnDevice(dev_track* tracks,
										dev_material* materials) {

	  int tid = threadIdx.x + blockIdx.x*blockDim.x;

	  dev_track* curr_track;
	  dev_segment* curr_segment;
	  int num_segments;

	  /* Iterate over all tracks */
	  while (tid < *_tot_num_tracks) {

		curr_track = &tracks[tid];
		num_segments = curr_track->_num_segments;

		/* Iterate over all segments */
		for (int i=0; i < num_segments; i++) {

			curr_segment = &curr_track->_segments[i];

			/* Iterate over all energy groups and polar angles */
			for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
				for (int p = 0; p < NUM_POLAR_ANGLES; p++) {

					curr_segment->_prefactors[e][p] = computePrefactorOnDevice(
							materials[curr_segment->_material_uid]._sigma_t[e],
							curr_segment->_length, _sinthetas[p]);
				}
			}
		}

		/* Update thread id */
		tid += blockDim.x*gridDim.x;
	  }

	  return;
}
#endif


/**
 * Compute the the prefactor for a single segment, energy group, and
 * polar angle
 * @param sigma_t the total cross-section
 * @param length the lenght of the segment
 * @param sintheta the sin of the polar angle
 */
__device__ FP_PRECISION computePrefactorOnDevice(FP_PRECISION sigma_t,
						FP_PRECISION length, FP_PRECISION sintheta) {

	FP_PRECISION prefactor = 1.0 - exp (-sigma_t * length / sintheta);
	return prefactor;
}


/**
 * Sets the boundary fluxes for each track to zero
 * @param tracks pointer to an array of tracks on the device
 * @param tot_num_tracks pointer to an int for the total number of tracks
 */
__global__ void zeroTrackFluxesOnDevice(dev_track* tracks) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	dev_track* curr_track;
	FP_PRECISION* polar_fluxes;

	/* Iterate over all tracks */
	while (tid < *_tot_num_tracks) {

		curr_track = &tracks[tid];
		polar_fluxes = curr_track->_polar_fluxes;

		/* Iterate over all polar angles and energy groups */
		for (int i = 0; i < GRP_TIMES_ANG * 2; i++)
			polar_fluxes[i] = 0.0;

		/* Increment thread id */
		tid += blockDim.x * gridDim.x;
	}

	return;
}



/**
 * Normalize the boundary fluxes for each track
 * @param tracks pointer to an array of tracks on the device
 * @param tot_num_tracks pointer to an int for the total number of tracks
 * @param renorm_factor factor to normalize each flux by
 */
__global__ void normalizeTrackFluxesOnDevice(dev_track* tracks,
										FP_PRECISION* renorm_factor) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	FP_PRECISION factor = *renorm_factor;
	dev_track* curr_track;
	FP_PRECISION* polar_fluxes;

	/* Iterate over all tracks */
	while (tid < *_tot_num_tracks) {

		curr_track = &tracks[tid];
		polar_fluxes = curr_track->_polar_fluxes;

		for (int i=0; i < GRP_TIMES_ANG * 2; i++)
			polar_fluxes[i] *= factor;

		/* Increment thread id */
		tid += blockDim.x * gridDim.x;
	}

	return;
}


/**
 * Initialize the flat source regions on the device by setting their
 * fluxes and old sources (from previous iteration) to 1.0
 * @param FSRs pointer to an array of
 * @param num_FSRs the number of flat source regions on the device
 */
__global__ void initDeviceFSRsOnDevice(dev_flatsourceregion* FSRs) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	dev_flatsourceregion* curr_FSR;

	/* Iterate over all flat source regions */
	while (tid < *_num_fsrs) {
		curr_FSR = &FSRs[tid];

		for (int i=0; i < NUM_ENERGY_GROUPS; i++) {
			curr_FSR->_flux[i] = 1.0;
			curr_FSR->_new_source[i] = 1.0;
			curr_FSR->_old_source[i] = 1.0;
		}

		/* Increment the thread id */
		tid +=  blockDim.x * gridDim.x;
	}

	return;
}


/**
 * Set the FSR fluxes to zero for each flat source region
 * @param FSRs pointer to array of flat source regions on the device
 * @param nuM_FSRs pointer to an int of the number of flat source regions
 */
__global__ void zeroFSRFluxesOnDevice(dev_flatsourceregion* FSRs) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	dev_flatsourceregion* curr_FSR;

	/* Iterate over all flat source regions */
	while (tid < *_num_fsrs) {
		curr_FSR = &FSRs[tid];

		/* Iterate over all energy groups */
		for (int i=0; i < NUM_ENERGY_GROUPS; i++)
			curr_FSR->_flux[i] = 0.0;

		/* Increment the thread id */
		tid +=  blockDim.x * gridDim.x;
	}

	return;
}


/**
 * Normalize the flux in each flat source region
 * @param FSRs pointer to flat source regions array on the device
 * @param num_FSRs pointer to an int of the number of flat source regions
 * @param renorm_factor the factor to normalize each flux by
 */
__global__ void normalizeFSRFluxesOnDevice(dev_flatsourceregion* FSRs,
											FP_PRECISION* renorm_factor) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	FP_PRECISION factor = *renorm_factor;
	dev_flatsourceregion* curr_FSR;

	/* Iterate over all flat source regions */
	while (tid < *_num_fsrs) {
		curr_FSR = &FSRs[tid];

		/* Iterate over each energy group */
		for (int i=0; i < NUM_ENERGY_GROUPS; i++)
			  curr_FSR->_flux[i] *= factor;

		 /* Increment thread id */
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
__global__ void computeTotalFissionSourceOnDevice(dev_flatsourceregion* FSRs,
											dev_material* materials,
											FP_PRECISION* fission_source) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	extern __shared__ FP_PRECISION shared_fission_source[];
	dev_flatsourceregion* curr_FSR;
	dev_material* curr_material;
	FP_PRECISION* nu_sigma_f;
	FP_PRECISION* scalar_flux;
	FP_PRECISION volume;

	/* Initialize fission source to zero */
	shared_fission_source[threadIdx.x] = 0;

	/* Iterate over all FSRs */
	while (tid < *_num_fsrs) {

		curr_FSR = &FSRs[tid];
		curr_material = &materials[curr_FSR->_material_uid];
		nu_sigma_f = curr_material->_nu_sigma_f;
		scalar_flux = curr_FSR->_flux;
		volume = curr_FSR->_volume;

		/* Iterate over all energy groups and update
		 * fission source for this block */
		for (int i=0; i < NUM_ENERGY_GROUPS; i++)
			shared_fission_source[threadIdx.x] +=
								nu_sigma_f[i] * scalar_flux[i] * volume;

		/* Increment thread id */
		tid += blockDim.x * gridDim.x;
   }

	/* Copy this threads fission source to global memory */
	tid = threadIdx.x + blockIdx.x * blockDim.x;
	fission_source[tid] = shared_fission_source[threadIdx.x];

	return;
}


/**
 * @param FSRs pointer to the flat source region array on the device
 * @param num_FSRs pointer to an int of the number of flat source regions
 * @param materials pointer an array of materials on the device
 */
 __global__ void computeFSRSourcesOnDevice(dev_flatsourceregion* FSRs,
										dev_material* materials,
										FP_PRECISION* curr_keff,
										FP_PRECISION* source_residual_norm) {

	 extern __shared__ FP_PRECISION shared_tot_abs_and_fission[];
	 int tid = threadIdx.x + blockIdx.x * blockDim.x;
	 shared_tot_abs_and_fission[threadIdx.x] = 0.0;
	 shared_tot_abs_and_fission[threadIdx.x+blockDim.x] = 0.0;

	 /* Reset the residual for the old and new fission sources to zero */
	 source_residual_norm[threadIdx.x + blockIdx.x * blockDim.x] = 0.0;

	 FP_PRECISION abs, fission, scatter;

	 FP_PRECISION volume;

	 dev_flatsourceregion* curr_FSR;
	 dev_material* curr_material;
	 FP_PRECISION* flux;
	 FP_PRECISION* new_source;
	 FP_PRECISION* old_source;
	 FP_PRECISION* sigma_a;
	 FP_PRECISION* nu_sigma_f;
	 FP_PRECISION* chi;
	 FP_PRECISION* sigma_s;
	 FP_PRECISION* sigma_t;
	 FP_PRECISION* ratios;

	 /* Iterate over all FSRs */
	 while (tid < *_num_fsrs) {

		 abs = 0;
		 fission = 0;
		 scatter = 0;

		 curr_FSR = &FSRs[tid];
		 flux = curr_FSR->_flux;
		 new_source = curr_FSR->_new_source;
		 old_source = curr_FSR->_old_source;
		 volume = curr_FSR->_volume;

		 curr_material = &materials[curr_FSR->_material_uid];
		 nu_sigma_f = curr_material->_nu_sigma_f;
		 sigma_a = curr_material->_sigma_a;
		 sigma_s = *curr_material->_sigma_s;
		 sigma_t = materials[curr_FSR->_material_uid]._sigma_t;
		 chi = curr_material->_chi;
		 ratios = curr_FSR->_ratios;

		 /* Compute total fission and absorption for region */
		 for (int i=0; i < NUM_ENERGY_GROUPS; i++) {
			 abs += flux[i] * sigma_a[i];
			 fission += flux[i] * nu_sigma_f[i];
		 }

		 /* Compute total scatter source for region */
		 for (int i=0; i < NUM_ENERGY_GROUPS; i++) {

			 scatter = 0;

			 /* Iterate over all energy groups */
			 for (int j=0; j < NUM_ENERGY_GROUPS; j++)
				 scatter += sigma_s[i*NUM_ENERGY_GROUPS + j] * flux[j];

			 /* Set the total source for this region in this group */
			 new_source[i] = (__fdividef(1.0, *curr_keff) * fission * chi[i] +
								scatter) * ONE_OVER_FOUR_PI;

			 //NOTE: You must index int a pointer to a thrust
			 //vector using the threadIdx.x, blockIdx.x
			 //and blockDim.x as shown below
			 if (fabs(new_source[i]) > 1E-10)
				 source_residual_norm[threadIdx.x + blockIdx.x * blockDim.x] +=
								 pow((new_source[i] - old_source[i]) /
										 	 	 	 	 	 new_source[i], 2);

			 old_source[i] = new_source[i];

			 ratios[i] = __fdividef(new_source[i], sigma_t[i]);

		 }

		 abs *= volume;
		 fission *= volume;

		 /* Update the temporary array for this block */
		 shared_tot_abs_and_fission[threadIdx.x] += abs;
		 shared_tot_abs_and_fission[threadIdx.x+blockDim.x] += fission;

		 /* Increment the thread id */
		 tid += blockDim.x * gridDim.x;
	 }

	 return;
}



/**
 * Computes the total fission and absorption rates from all flat source regions
 * and stores them in a device vector in global memory on the device
 * @param FSRs pointer to the flat source region array on the device
 * @param num_FSRs pointer to an int of the number of flat source regions
 * @param materials pointer an array of materials on the device
 * @param tot_abs device_vector of cumulative absorption rates from each block
 * @param tot_fission device_vector of cumulative fission rates froms each block
 */
__global__ void computeAbsAndFissionOnDevice(dev_flatsourceregion* FSRs,
							dev_material* materials,
							FP_PRECISION* tot_abs, FP_PRECISION* tot_fission) {

	extern __shared__ FP_PRECISION shared_tot_abs_and_fission[];
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	dev_flatsourceregion* curr_fsr;
	dev_material* curr_material;
	shared_tot_abs_and_fission[threadIdx.x] = 0.0;
	shared_tot_abs_and_fission[threadIdx.x+blockDim.x] = 0.0;
	FP_PRECISION abs = 0;
	FP_PRECISION fission;
	FP_PRECISION* sigma_a;
	FP_PRECISION* nu_sigma_f;
	FP_PRECISION* flux;
	FP_PRECISION volume;

	/* Iterate over all FSRs */
	while (tid < *_num_fsrs) {

		abs = 0;
		fission = 0;
		curr_fsr = &FSRs[tid];
		flux = curr_fsr->_flux;
		volume = curr_fsr->_volume;

		curr_material = &materials[curr_fsr->_material_uid];
		sigma_a = curr_material->_sigma_a;
		nu_sigma_f = curr_material->_nu_sigma_f;

		/* Iterate over all energy groups */
		for (int i=0; i < NUM_ENERGY_GROUPS; i++) {
			abs += sigma_a[i] * flux[i] * volume;
			fission += nu_sigma_f[i] * flux[i] * volume;
		}

		/* Update the temporary array for this block */
		shared_tot_abs_and_fission[threadIdx.x] += abs;
		shared_tot_abs_and_fission[threadIdx.x+blockDim.x] += fission;

		/* Increment the thread id */
		tid += blockDim.x * gridDim.x;
	}

	/* Store the shared memory total absorption and fission to global memory */
	//NOTE: You must index int a pointer to a thrust
	//vector using the threadIdx.x, blockIdx.x
	//and blockDim.x as shown below
	tot_abs[threadIdx.x + blockIdx.x * blockDim.x] =
							shared_tot_abs_and_fission[threadIdx.x];
	tot_fission[threadIdx.x + blockIdx.x * blockDim.x] =
							shared_tot_abs_and_fission[threadIdx.x+blockDim.x];

	return;
}


/**
 * Normalizes the flux to the volume of each FSR and adds in the source term
 * computed and stored in the ratios attribute for each FSR
 * @param FSRs pointer to the flat source region array on the device
 * @param num_FSRs pointer to an int of the number of flat source regions
 * @param materials pointer an array of materials on the device
 */
__global__ void normalizeFluxToVolumeOnDevice(dev_flatsourceregion* FSRs,
									dev_material* materials) {

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	dev_flatsourceregion* curr_FSR;
	FP_PRECISION* ratios;
	FP_PRECISION volume;

	dev_material* curr_material;
	FP_PRECISION* sigma_t;

	/* Iterate over all FSRs */
	while (tid < *_num_fsrs) {

		curr_FSR = &FSRs[tid];
		curr_material = &materials[curr_FSR->_material_uid];
		ratios = curr_FSR->_ratios;
		volume = curr_FSR->_volume;
		sigma_t = curr_material->_sigma_t;

		/* Iterate over all energy groups */
		for (int i=0; i < NUM_ENERGY_GROUPS; i++) {
			curr_FSR->_flux[i] = __fdividef(curr_FSR->_flux[i], 2.0);
			curr_FSR->_flux[i] = FOUR_PI * ratios[i] +
					__fdividef(curr_FSR->_flux[i], (sigma_t[i] * volume));
		}

		/* Increment thread id */
		tid += blockDim.x * gridDim.x;
	}

   return;
}



#if STORE_PREFACTORS
/**
 * This kernel integrates the neutron transport equation in the "forward"
 * direction along each track in the geometry using exponential prefactors
 * which are precomputed and stored for each segment, energy group and polar
 * angle
 * @param track a pointer to the array of tracks on the device
 * @param num_azim a pointer the number of azimuthal angles on the device
 * @param track_index_offsets a pointer to an array of cumulative track indices
 * @param materials a pointer to the array of materials on the device
 * @param FSRs a pointer to the flat source regions on the device
 */
__global__ void propagateTrackFluxForwardOnDevice(dev_track* tracks,
										int* track_index_offsets,
										dev_flatsourceregion* FSRs) {
#else
/**
 * This kernel integrates the neutron transport equation in the "forward"
 * direction along each track in the geometry using exponential prefactors
 * which are precomputed and stored in a hash table for O(1) lookup and
 * interpolation
 * @param track a pointer to the array of tracks on the device
 * @param num_azim a pointer the number of azimuthal angles on the device
 * @param track_index_offsets a pointer to an array of cumulative track indices
 * @param materials a pointer to the array of materials on the device
 * @param FSRs a pointer to the flat source regions on the device
 * @param prefactor_array a pointer to the hash table on the device
 * @param prefactor_spacing a pointer to the interval between prefactors
 * @param prefactor_max_index a pointer to the number of entries in the table
 */
__global__ void propagateTrackFluxForwardOnDevice(dev_track* tracks,
										int* track_index_offsets,
										dev_material* materials,
										dev_flatsourceregion* FSRs,
										FP_PRECISION* prefactor_array) {
#endif

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int index_offset = threadIdx.x * (2*NUM_POLAR_ANGLES + 1);
	int energy_group = tid % NUM_ENERGY_GROUPS;
	int energy_angle_index = energy_group * NUM_POLAR_ANGLES;
	int fsr_flux_index = index_offset + 2*NUM_POLAR_ANGLES;
	int track_flux_index;
	int azim_angle_index;

	dev_track* curr_track;
	dev_segment* curr_segment;
	dev_flatsourceregion* curr_FSR;
	int num_segments;
	FP_PRECISION delta;

	/* temporary flux for track and fsr */
	extern __shared__ FP_PRECISION temp_flux[];
	int index;

#if !STORE_PREFACTORS
	dev_material* curr_material;
	FP_PRECISION sigma_t_l;
#endif

	/* Iterate over half of the tracks */
	while (int(tid / NUM_ENERGY_GROUPS) < track_index_offsets[*_num_azim / 2]) {

		curr_track = &tracks[int(tid / NUM_ENERGY_GROUPS)];
		azim_angle_index = curr_track->_azim_angle_index;
		num_segments = curr_track->_num_segments;

		/* Put tracks flux in the shared memory temporary flux array */
		for (int p=0; p < NUM_POLAR_ANGLES; p++) {

			temp_flux[index_offset + p] =
							curr_track->_polar_fluxes[energy_angle_index + p];
//			if (curr_track->_refl_out)
//			if (curr_track->_refl_in)
//				temp_flux[index_offset + p] *= curr_track->_bc_out;
//			else
//				temp_flux[index_offset + p] *= curr_track->_bc_in;

			temp_flux[index_offset + NUM_POLAR_ANGLES + p] =
				curr_track->_polar_fluxes[GRP_TIMES_ANG+energy_angle_index+p];
//			if (curr_track->_refl_out)
//			if (curr_track->_refl_in)
//				temp_flux[index_offset + NUM_POLAR_ANGLES + p] *= curr_track->_bc_out;
////				temp_flux[index_offset + NUM_POLAR_ANGLES + p] *= curr_track->_bc_in;
//			else
//				temp_flux[index_offset + NUM_POLAR_ANGLES + p] *= curr_track->_bc_in;
		}

		/* Iterate over each segment in forward direction */
		track_flux_index = index_offset;

		for (int i=0; i < num_segments; i++) {
			curr_segment = &curr_track->_segments[i];
			curr_FSR = &FSRs[curr_segment->_region_uid];


			/* Iterate over all energy groups */
			temp_flux[fsr_flux_index] = 0.0;

#if STORE_PREFACTORS

				/* Iterate over all polar angles for this segment */
				for (int p=0; p < NUM_POLAR_ANGLES; p++) {
					delta = (temp_flux[track_flux_index + p] -
							curr_FSR->_ratios[energy_group]) *
							curr_segment->_prefactors[energy_group][p];
					temp_flux[fsr_flux_index] += delta *
						_polar_weights[azim_angle_index*NUM_POLAR_ANGLES + p];
					temp_flux[track_flux_index + p] -= delta;
				}
#else
				curr_material = &materials[curr_segment->_material_uid];

				sigma_t_l = min(curr_material->_sigma_t[energy_group] *
												curr_segment->_length, 10.0);
				index = min(int(__fdividef(sigma_t_l, *_prefactor_spacing)) *
										2.0 * NUM_POLAR_ANGLES,
										(FP_PRECISION)*_prefactor_max_index);

				for (int p=0; p < NUM_POLAR_ANGLES; p++) {
					delta = (temp_flux[track_flux_index + p] -
							curr_FSR->_ratios[energy_group]) *
							(1.0 - (prefactor_array[index + 2*p] * sigma_t_l +
							prefactor_array[index + 2*p + 1]));
					temp_flux[fsr_flux_index] += delta *
						_polar_weights[azim_angle_index*NUM_POLAR_ANGLES + p];
					temp_flux[track_flux_index + p] -= delta;
				}

#endif

			/* Update the scalar flux for this fsr */
			atomicAdd(&curr_FSR->_flux[energy_group], temp_flux[fsr_flux_index]);
		}

		/* Transfer flux to outgoing track */
		index = curr_track->_refl_out * GRP_TIMES_ANG;

		for (int p=0; p < NUM_POLAR_ANGLES; p++) {
			if (curr_track->_refl_out)
				temp_flux[track_flux_index + p] *= curr_track->_bc_in;
			else
				temp_flux[track_flux_index + p] *= curr_track->_bc_out;
		}

		memcpy((void*)&tracks[curr_track->_track_out]._polar_fluxes[index + energy_angle_index],
								(void*)&temp_flux[track_flux_index],
								sizeof(FP_PRECISION)*NUM_POLAR_ANGLES);

		/* Loop over each segment in reverse direction */
		track_flux_index = index_offset + NUM_POLAR_ANGLES;

		for (int i=num_segments-1; i > -1; i--) {
			curr_segment = &curr_track->_segments[i];
			curr_FSR = &FSRs[curr_segment->_region_uid];

			temp_flux[fsr_flux_index] = 0.0;

#if STORE_PREFACTORS

				/* Iterate over all polar angles for this segment */
				for (int p=0; p < NUM_POLAR_ANGLES; p++) {
					delta = (temp_flux[track_flux_index + p] -
							curr_FSR->_ratios[energy_group]) *
							curr_segment->_prefactors[energy_group][p];
					temp_flux[fsr_flux_index] += delta *
						_polar_weights[azim_angle_index*NUM_POLAR_ANGLES + p];
					temp_flux[track_flux_index + p] -= delta;
				}
#else
				curr_material = &materials[curr_segment->_material_uid];

				sigma_t_l = min(curr_material->_sigma_t[energy_group] *
										curr_segment->_length, 10.0);
				index = min(int(__fdividef(sigma_t_l, *_prefactor_spacing)) *
										2.0 * NUM_POLAR_ANGLES,
										(FP_PRECISION)*_prefactor_max_index);

				for (int p=0; p < NUM_POLAR_ANGLES; p++) {
					delta = (temp_flux[track_flux_index + p] -
							curr_FSR->_ratios[energy_group]) *
							(1.0 - (prefactor_array[index + 2*p] * sigma_t_l +
							prefactor_array[index + 2*p + 1]));
					temp_flux[fsr_flux_index] += delta *
						_polar_weights[azim_angle_index*NUM_POLAR_ANGLES + p];
					temp_flux[track_flux_index + p] -= delta;
				}


#endif

			/* Update the scalar flux for this fsr */
			atomicAdd(&curr_FSR->_flux[energy_group], temp_flux[fsr_flux_index]);
		}


		/* Transfer flux to outgoing track */
		index = curr_track->_refl_in * GRP_TIMES_ANG;

		for (int p=0; p < NUM_POLAR_ANGLES; p++) {
			if (curr_track->_refl_in)
				temp_flux[track_flux_index + p] *= curr_track->_bc_in;
			else
				temp_flux[track_flux_index + p] *= curr_track->_bc_out;
		}

		memcpy((void*)&tracks[curr_track->_track_in]._polar_fluxes[index + energy_angle_index],
									(void*)&temp_flux[track_flux_index],
									sizeof(FP_PRECISION)*NUM_POLAR_ANGLES);


		tid += blockDim.x * gridDim.x;
		energy_group = tid % NUM_ENERGY_GROUPS;
		energy_angle_index = energy_group * NUM_POLAR_ANGLES;
	}

	return;
}

#if STORE_PREFACTORS
/**
 * This kernel integrates the neutron transport equation in the "forward"
 * direction along each track in the geometry using exponential prefactors
 * which are precomputed and stored for each segment, energy group and polar
 * angle
 * @param track a pointer to the array of tracks on the device
 * @param num_azim a pointer the number of azimuthal angles on the device
 * @param track_index_offsets a pointer to an array of cumulative track indices
 * @param materials a pointer to the array of materials on the device
 * @param FSRs a pointer to the flat source regions on the device
 */
__global__ void propagateTrackFluxReverseOnDevice(dev_track* tracks,
										int* track_index_offsets,
										dev_flatsourceregion* FSRs) {
#else
/**
 * This kernel integrates the neutron transport equation in the "forward"
 * direction along each track in the geometry using exponential prefactors
 * which are precomputed and stored in a hash table for O(1) lookup and
 * interpolation
 * @param track a pointer to the array of tracks on the device
 * @param num_azim a pointer the number of azimuthal angles on the device
 * @param track_index_offsets a pointer to an array of cumulative track indices
 * @param materials a pointer to the array of materials on the device
 * @param FSRs a pointer to the flat source regions on the device
 * @param prefactor_array a pointer to the hash table on the device
 * @param prefactor_spacing a pointer to the interval between prefactors
 * @param prefactor_max_index a pointer to the number of entries in the table
 */
__global__ void propagateTrackFluxReverseOnDevice(dev_track* tracks,
									int* track_index_offsets,
									dev_material* materials,
									dev_flatsourceregion* FSRs,
									FP_PRECISION* prefactor_array) {
#endif

	int tid = track_index_offsets[*_num_azim / 2]*NUM_ENERGY_GROUPS +
								threadIdx.x + blockIdx.x * blockDim.x;
	int index_offset = threadIdx.x * (2*NUM_POLAR_ANGLES + 1);
	int energy_group = tid % NUM_ENERGY_GROUPS;
	int energy_angle_index = energy_group * NUM_POLAR_ANGLES;
	int fsr_flux_index = index_offset + 2*NUM_POLAR_ANGLES;
	int track_flux_index;
	int azim_angle_index;

	dev_track* curr_track;
	dev_segment* curr_segment;
	dev_flatsourceregion* curr_FSR;
	int num_segments;
	FP_PRECISION delta;

	/* temporary flux for track and fsr */
	extern __shared__ FP_PRECISION temp_flux[];
	int index;

#if !STORE_PREFACTORS
	dev_material* curr_material;
	FP_PRECISION sigma_t_l;
#endif

	while (int(tid / NUM_ENERGY_GROUPS) < track_index_offsets[*_num_azim]) {

		curr_track = &tracks[int(tid / NUM_ENERGY_GROUPS)];
		azim_angle_index = curr_track->_azim_angle_index;
		num_segments = curr_track->_num_segments;

		/* Put tracks flux in the shared memory temporary flux array
		 * (forward flux) */
		for (int p=0; p < NUM_POLAR_ANGLES; p++) {

			temp_flux[index_offset + p] =
					curr_track->_polar_fluxes[energy_angle_index + p];
//			if (curr_track->_refl_in)
//			if (curr_track->_refl_out)
//				temp_flux[index_offset + p] *= curr_track->_bc_out;
//			else
//				temp_flux[index_offset + p] *= curr_track->_bc_in;

			temp_flux[index_offset + NUM_POLAR_ANGLES + p] =
					curr_track->_polar_fluxes[energy_angle_index +
					                          GRP_TIMES_ANG + p];
//			if (curr_track->_refl_out)
//			if (curr_track->_refl_in)
//				temp_flux[index_offset + NUM_POLAR_ANGLES + p] *= curr_track->_bc_out;
//				temp_flux[index_offset + NUM_POLAR_ANGLES + p] *= curr_track->_bc_in;
//			else
//				temp_flux[index_offset + NUM_POLAR_ANGLES + p] *= curr_track->_bc_in;
//				temp_flux[index_offset + NUM_POLAR_ANGLES + p] *= curr_track->_bc_out;
		}

		/* Iterate over each segment in forward direction */
		track_flux_index = index_offset;
		for (int i=0; i < num_segments; i++) {
			curr_segment = &curr_track->_segments[i];
			curr_FSR = &FSRs[curr_segment->_region_uid];

			/* Iterate over all energy groups */
			temp_flux[fsr_flux_index] = 0.0;

#if STORE_PREFACTORS

				/* Iterate over all polar angles for this segment */
				for (int p=0; p < NUM_POLAR_ANGLES; p++) {
					delta = (temp_flux[track_flux_index + p] -
							curr_FSR->_ratios[energy_group]) *
							curr_segment->_prefactors[energy_group][p];
					temp_flux[fsr_flux_index] += delta *
						_polar_weights[azim_angle_index*NUM_POLAR_ANGLES + p];
					temp_flux[track_flux_index + p] -= delta;
				}
#else
				curr_material = &materials[curr_segment->_material_uid];

				sigma_t_l = min(curr_material->_sigma_t[energy_group] *
								curr_segment->_length, 10.0);
				index = min(int(__fdividef(sigma_t_l, *_prefactor_spacing)) *
								2.0 * NUM_POLAR_ANGLES,
								(FP_PRECISION)*_prefactor_max_index);

				for (int p=0; p < NUM_POLAR_ANGLES; p++) {
					delta = (temp_flux[track_flux_index + p] -
							curr_FSR->_ratios[energy_group]) *
							(1.0 - (prefactor_array[index + 2*p] * sigma_t_l +
							prefactor_array[index + 2*p + 1]));
					temp_flux[fsr_flux_index] += delta *
						_polar_weights[azim_angle_index*NUM_POLAR_ANGLES + p];
					temp_flux[track_flux_index + p] -= delta;
				}

#endif

			/* Update the scalar flux for this fsr */
			atomicAdd(&curr_FSR->_flux[energy_group],
							temp_flux[fsr_flux_index]);
		}

		/* Transfer flux to outgoing track */
		index = curr_track->_refl_out * GRP_TIMES_ANG;

		for (int p=0; p < NUM_POLAR_ANGLES; p++) {
			if (curr_track->_refl_out)
				temp_flux[track_flux_index + p] *= curr_track->_bc_in;
			else
				temp_flux[track_flux_index + p] *= curr_track->_bc_out;
		}

		memcpy((void*)&tracks[curr_track->_track_out]._polar_fluxes[index + energy_angle_index],
									(void*)&temp_flux[track_flux_index],
									sizeof(FP_PRECISION)*NUM_POLAR_ANGLES);


		/* Loop over each segment in reverse direction */
		track_flux_index = index_offset + NUM_POLAR_ANGLES;
		for (int i=num_segments-1; i > -1; i--) {
			curr_segment = &curr_track->_segments[i];
			curr_FSR = &FSRs[curr_segment->_region_uid];

			temp_flux[fsr_flux_index] = 0.0;

#if STORE_PREFACTORS

				/* Iterate over all polar angles for this segment */
				for (int p=0; p < NUM_POLAR_ANGLES; p++) {
					delta = (temp_flux[track_flux_index + p] -
							curr_FSR->_ratios[energy_group]) *
									curr_segment->_prefactors[energy_group][p];
					temp_flux[fsr_flux_index] += delta *
						_polar_weights[azim_angle_index*NUM_POLAR_ANGLES + p];
					temp_flux[track_flux_index + p] -= delta;
				}
#else
				curr_material = &materials[curr_segment->_material_uid];

				sigma_t_l = min(curr_material->_sigma_t[energy_group] *
											curr_segment->_length, 10.0);
				index = min(int(__fdividef(sigma_t_l, *_prefactor_spacing)) *
										2.0 * NUM_POLAR_ANGLES,
										(FP_PRECISION)*_prefactor_max_index);

				for (int p=0; p < NUM_POLAR_ANGLES; p++) {
					delta = (temp_flux[track_flux_index + p] -
							curr_FSR->_ratios[energy_group]) *
							(1.0 - (prefactor_array[index + 2*p] * sigma_t_l +
							prefactor_array[index + 2*p + 1]));
					temp_flux[fsr_flux_index] += delta *
						_polar_weights[azim_angle_index*NUM_POLAR_ANGLES + p];
					temp_flux[track_flux_index + p] -= delta;
				}


#endif

			/* Update the scalar flux for this fsr */
			atomicAdd(&curr_FSR->_flux[energy_group], temp_flux[fsr_flux_index]);
		}


		/* Transfer flux to outgoing track */
		index = curr_track->_refl_in * GRP_TIMES_ANG;

		for (int p=0; p < NUM_POLAR_ANGLES; p++) {
			if (curr_track->_refl_in)
				temp_flux[track_flux_index + p] *= curr_track->_bc_in;
			else
				temp_flux[track_flux_index + p] *= curr_track->_bc_out;
		}

		memcpy((void*)&tracks[curr_track->_track_in]._polar_fluxes[index + energy_angle_index],
								(void*)&temp_flux[track_flux_index],
								sizeof(FP_PRECISION)*NUM_POLAR_ANGLES);


		tid += blockDim.x * gridDim.x;
		energy_group = tid % NUM_ENERGY_GROUPS;
		energy_angle_index = energy_group * NUM_POLAR_ANGLES;
	}

	return;
}


__device__ double atomicAdd(double* address, double val) {

	unsigned long long int* address_as_ull =
						(unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;

	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed,
		__double_as_longlong(val +
		__longlong_as_double(assumed)));
	} while (assumed != old);

	return __longlong_as_double(old);
}
