/*
 * Solver.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Solver.h"


/**
 * Solver constructor
 * @param geom pointer to the geometry
 * @param track_generator pointer to the trackgenerator
 */
Solver::Solver(Geometry* geom, TrackGenerator* track_generator,
											Plotter* plotter) {

	_geom = geom;
	_quad = new Quadrature(TABUCHI);
	_num_FSRs = geom->getNumFSRs();
	_tracks = track_generator->getTracks();
	_num_tracks = track_generator->getNumTracks();
	_num_azim = track_generator->getNumAzim();
	_k_eff = 1.0;

	/* Count the total number of tracks */
	_tot_num_tracks = 0;
	for (int i=0; i < _num_azim; i++)
		_tot_num_tracks += _num_tracks[i];

	_plotter = plotter;

	try{
		_FSRs = new FlatSourceRegion[_num_FSRs];
		_FSRs_to_powers = new FP_PRECISION[_num_FSRs];
		_FSRs_to_pin_powers = new FP_PRECISION[_num_FSRs];

		for (int e = 0; e <= NUM_ENERGY_GROUPS; e++)
			_FSRs_to_fluxes[e] = new FP_PRECISION[_num_FSRs];
	}
	catch(std::exception &e) {
		log_printf(ERROR, "Could not allocate memory for the solver's flat "
					"source region array. Backtrace:%s", e.what());
	}

	/* Pre-compute exponential pre-factors */
	precomputePrefactors();
	initializeFSRs();
}


/**
 * Solver destructor deletes flat source regions array
 */
Solver::~Solver() {

	delete [] _FSRs;
	delete [] _FSRs_to_powers;
	delete [] _FSRs_to_pin_powers;
	delete _quad;

	for (int e = 0; e <= NUM_ENERGY_GROUPS; e++)
		delete [] _FSRs_to_fluxes[e];

#if !STORE_PREFACTORS
	delete [] _prefactor_array;
#endif

}

 
/**
 * Pre-computes exponential pre-factors for each segment of each track for
 * each polar angle. This method will store each pre-factor in an array inside
 * each segment if STORE_PREFACTORS is set to true inside the configurations.h
 * file. If it is not set to true then a hashmap will be generated which will
 * contain values of the pre-factor at for specific segment lengths (the keys
 * into the hashmap).
 */
void Solver::precomputePrefactors() {

	log_printf(INFO, "Pre-computing exponential pre-factors...");

	Track* curr_track;
	FP_PRECISION azim_weight;

	/* Precompute the total azimuthal weight for tracks at each polar angle */
	for (int i = 0; i < _num_azim; i++) {
		for (int j = 0; j < _num_tracks[i]; j++) {
			curr_track = &_tracks[i][j];
			azim_weight = curr_track->getAzimuthalWeight();

			for (int p = 0; p < NUM_POLAR_ANGLES; p++)
				curr_track->setPolarWeight(p,
									azim_weight*_quad->getMultiple(p)*FOUR_PI);
		}
	}


/*Store pre-factors inside each segment */
#if STORE_PREFACTORS

	log_printf(INFO, "Pre-factors will be stored inside each segment...");

	segment* curr_seg;

	/* Loop over azimuthal angle, track, segment, polar angle, energy group */
	for (int i = 0; i < _num_azim; i++) {
		for (int j = 0; j < _num_tracks[i]; j++) {
			curr_track = &_tracks[i][j];

			for (int s = 0; s < curr_track->getNumSegments(); s++) {
				curr_seg = curr_track->getSegment(s);

				for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
					for (int p = 0; p < NUM_POLAR_ANGLES; p++)
					  curr_seg->_prefactors[e][p] =
						  computePrefactor(curr_seg->_material->getSigmaT()[e],
						  	  	  curr_seg->_length, _quad->getSinTheta(p));
				}
			}
		}
	}



/* Use hash map */
#else

	/* make pre factor array based on table look up with linear interpolation */

	log_printf(NORMAL, "Making Prefactor array...");

	/* set size of prefactor array */
	int num_array_values = 10 * sqrt(1 / (8 * SOURCE_CONVERG_THRESH));
	_prefactor_spacing = 10.0 / num_array_values;
	_prefactor_array_size = 2 * NUM_POLAR_ANGLES * num_array_values;
	_prefactor_max_index = _prefactor_array_size - 2*NUM_POLAR_ANGLES - 1;

	log_printf(DEBUG, "prefactor array size: %i, max index: %i",
					_prefactor_array_size, _prefactor_max_index);

	/* allocate arrays */
	_prefactor_array = new FP_PRECISION[_prefactor_array_size];

	FP_PRECISION expon;
	FP_PRECISION intercept;
	FP_PRECISION slope;

	/* Create prefactor array */
	for (int i = 0; i < num_array_values; i ++){
		for (int j = 0; j < NUM_POLAR_ANGLES; j++){
			expon = exp(- (i * _prefactor_spacing) / _quad->getSinTheta(j));
			slope = - expon / _quad->getSinTheta(j);
			intercept = expon * (1 + (i * _prefactor_spacing) /
													_quad->getSinTheta(j));
			_prefactor_array[2 * NUM_POLAR_ANGLES*i + 2*j] = slope;
			_prefactor_array[2 * NUM_POLAR_ANGLES*i + 2*j + 1] = intercept;
		}
	}

#endif

	return;
}


/**
 * Function to compute the exponential prefactor for the transport equation for
 * a given segment
 * @param seg pointer to a segment
 * @param energy energy group index
 * @param angle polar angle index
 * @return the pre-factor
 */
FP_PRECISION Solver::computePrefactor(FP_PRECISION sigma_t,
								FP_PRECISION length, FP_PRECISION sintheta) {
        log_printf(INFO, "computing host prefactor...");
	FP_PRECISION prefactor = 1.0 - exp (-sigma_t * length / sintheta);
	return prefactor;
}


/**
 * Compute the ratio of source / sigma_t for each energy group in each flat
 * source region for efficient fixed source iteration
 */
void Solver::computeRatios() {

	#pragma omp parallel for
	for (int i = 0; i < _num_FSRs; i++)
		_FSRs[i].computeRatios();

	return;
}


void Solver::updateFSRFluxes() {


  FP_PRECISION* scalar_flux;
  FP_PRECISION* old_scalar_flux;
  FlatSourceRegion* fsr;
	
	#pragma omp parallel for private(fsr, scalar_flux, old_scalar_flux)
	for (int r = 0; r < _num_FSRs; r++) {
		fsr = &_FSRs[r];
		scalar_flux = fsr->getFlux();
		old_scalar_flux = fsr->getOldFlux();

		/* Update old scalar flux */
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
			old_scalar_flux[e] = scalar_flux[e];
	}

  return;
}


/**
 * Initializes each of the FlatSourceRegion objects inside the solver's
 * array of FSRs. This includes assigning each one a unique, monotonically
 * increasing id, setting the material for each FSR, and assigning a volume
 * based on the cumulative length of all of the segments inside the FSR.
 */
void Solver::initializeFSRs() {

	log_printf(NORMAL, "Initializing FSRs...");

	CellBasic* cell;
	Material* material;
	Universe* univ_zero = _geom->getUniverse(0);
	Track* track;
	segment* seg;
	FlatSourceRegion* fsr;

	/* Set each FSR's volume by accumulating the total length of all
	   tracks inside the FSR. Loop over azimuthal angle, track and segment */
	for (int i = 0; i < _num_azim; i++) {
		for (int j = 0; j < _num_tracks[i]; j++) {
			track = &_tracks[i][j];

			for (int s = 0; s < track->getNumSegments(); s++) {
				seg = track->getSegment(s);
				fsr = &_FSRs[seg->_region_id];
				fsr->incrementVolume(seg->_length*track->getAzimuthalWeight());
			}
		}
	}

	/* Loop over all FSRs */
	#pragma omp parallel for private(cell, material)
	for (int r = 0; r < _num_FSRs; r++) {

		/* Get the cell corresponding to this FSR from the geometry */
		cell = static_cast<CellBasic*>(_geom->findCell(univ_zero, r));

		/* Get the cell's material and assign it to the FSR */
		material = _geom->getMaterial(cell->getMaterial());
		_FSRs[r].setMaterial(material);

		log_printf(INFO, "FSR id = %d has cell id = %d and material id = %d "
				"and volume = %f", r, cell->getId(), material->getId(),
				_FSRs[r].getVolume());
	}

	return;
}


/**
 * Zero each track's incoming and outgoing polar fluxes
 */
void Solver::zeroTrackFluxes() {

	log_printf(INFO, "Setting all track polar fluxes to zero...");

	FP_PRECISION* polar_fluxes;

	/* Loop over azimuthal angle, track, polar angle, energy group
	 * and set each track's incoming and outgoing flux to zero */
	#pragma omp parallel for private(polar_fluxes)
	for (int i = 0; i < _num_azim; i++) {
		for (int j = 0; j < _num_tracks[i]; j++) {
			polar_fluxes = _tracks[i][j].getPolarFluxes();

			for (int i = 0; i < GRP_TIMES_ANG * 2; i++)
				polar_fluxes[i] = 0.0;
		}
	}
}


/**
 * Set the scalar flux for each energy group inside each FSR
 * to unity
 */
void Solver::initFSRs() {

	log_printf(INFO, "Setting all FSR scalar fluxes to unity...");
	FlatSourceRegion* fsr;

	/* Loop over all FSRs and energy groups */
	#pragma omp parallel for private(fsr)
	for (int r = 0; r < _num_FSRs; r++) {
		fsr = &_FSRs[r];
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
			fsr->setFlux(e, 1.0);
			fsr->setOldSource(e, 1.0);
		}
	}

	return;
}


/**
 * Set the scalar flux for each energy group inside each FSR to zero
 */
void Solver::zeroFSRFluxes() {

	log_printf(INFO, "Setting all FSR scalar fluxes to zero...");
	FlatSourceRegion* fsr;

	/* Loop over all FSRs and energy groups */
	for (int r = 0; r < _num_FSRs; r++) {
		fsr = &_FSRs[r];
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
			fsr->setFlux(e, 0.0);
	}

	return;
}


void Solver::renormalizeFluxes() {

	FP_PRECISION fission_source;
	FP_PRECISION renorm_factor, volume;
	FP_PRECISION* nu_sigma_f;
	FP_PRECISION* scalar_flux;
	FlatSourceRegion* fsr;
	Material* material;

	/* Initialize fission source to zero */
	fission_source = 0;

	/* Compute total fission source to zero for this region */
	for (int r = 0; r < _num_FSRs; r++) {

		/* Get pointers to important data structures */
		fsr = &_FSRs[r];
		material = fsr->getMaterial();
		nu_sigma_f = material->getNuSigmaF();
		scalar_flux = fsr->getFlux();
		volume = fsr->getVolume();

		for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
			fission_source += nu_sigma_f[e] * scalar_flux[e] * volume;
	}

	/* Renormalize scalar fluxes in each region */
	renorm_factor = 1.0 / fission_source;

	log_printf(DEBUG, "renorm_factor = %f", renorm_factor);

	#pragma omp parallel for
	for (int r = 0; r < _num_FSRs; r++)
		_FSRs[r].normalizeFluxes(renorm_factor);

	/* Renormalization angular boundary fluxes for each track */
	#pragma omp parallel for
	for (int j = 0; j < _num_azim; j++) {
		for (int k = 0; k < _num_tracks[j]; k++)
			_tracks[j][k].normalizeFluxes(renorm_factor);
	}

	return;
}


FP_PRECISION Solver::computeFSRSources() {

  FP_PRECISION scatter_source, fission_source;
  FP_PRECISION* nu_sigma_f;
  FP_PRECISION* sigma_s;
  FP_PRECISION* chi;
  FP_PRECISION* scalar_flux;
  FP_PRECISION* source;
  FP_PRECISION* old_source;
  FlatSourceRegion* fsr;
  Material* material;
  FP_PRECISION source_residual_norm = 0.0;

  /* For all regions, find the source */
  	  for (int r = 0; r < _num_FSRs; r++) {

  		  fsr = &_FSRs[r];
  		  material = fsr->getMaterial();

		/* Initialize the fission source to zero for this region */
		fission_source = 0;
		scalar_flux = fsr->getFlux();
		source = fsr->getSource();
		old_source = fsr->getOldSource();
		material = fsr->getMaterial();
		nu_sigma_f = material->getNuSigmaF();
		chi = material->getChi();
		sigma_s = material->getSigmaS();

		/* Compute total fission source for current region */
		for (int e = 0; e < NUM_ENERGY_GROUPS; e++)
		  fission_source += scalar_flux[e] * nu_sigma_f[e];

		/* Compute total scattering source for group G */
		for (int G = 0; G < NUM_ENERGY_GROUPS; G++) {
			scatter_source = 0;

			/* Compute total fission source for current region */
			for (int g = 0; g < NUM_ENERGY_GROUPS; g++)
				scatter_source += sigma_s[G*NUM_ENERGY_GROUPS + g] * scalar_flux[g];

			/* Set the total source for region r in group G */
			source[G] = ((1.0 / _k_eff) * fission_source *
									chi[G] + scatter_source) * ONE_OVER_FOUR_PI;

			/* Compute the norm of the residuals of the sources for convergence */
			if (fabs(source[G]) > 1E-10)
				source_residual_norm += pow((source[G]-old_source[G])/source[G], 2);

			/* Update the old source */
			old_source[G] = source[G];
		}
  	 }

  	 source_residual_norm = sqrt(source_residual_norm / _geom->getNumFSRs());

  	 return source_residual_norm;
}


/**
 * Compute k_eff from the new and old source and the value of k_eff from
 * the previous iteration
 */
void Solver::updateKeff() {

	FP_PRECISION tot_abs = 0.0;
	FP_PRECISION tot_fission = 0.0;
	FP_PRECISION abs = 0;
	FP_PRECISION fission;
	FP_PRECISION* sigma_a;
	FP_PRECISION* nu_sigma_f;
	FP_PRECISION* flux;
	Material* material;
	FlatSourceRegion* fsr;

	for (int r = 0; r < _num_FSRs; r++) {
		abs = 0;
		fission = 0;
		fsr = &_FSRs[r];
		material = fsr->getMaterial();
		sigma_a = material->getSigmaA();
		nu_sigma_f = material->getNuSigmaF();
		flux = fsr->getFlux();

		for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
			abs += sigma_a[e] * flux[e] * fsr->getVolume();
			fission += nu_sigma_f[e] * flux[e] * fsr->getVolume();
		}

			tot_abs += abs;

			tot_fission += fission;
	}

	_k_eff = tot_fission/tot_abs;
	log_printf(DEBUG, "tot_abs = %f, tot_fission = %f", tot_abs, tot_fission);

	return;
}


/**
 * Return an array indexed by FSR ids which contains the corresponding
 * fluxes for each FSR
 * @return an array map of FSR to fluxes
 */
FP_PRECISION** Solver::getFSRtoFluxMap() {
	return _FSRs_to_fluxes;
}


/**
 * Checks that each flat source region has at least one segment within it
 * and if not, exits the program with an error message
 */
void Solver::checkTrackSpacing() {

	int* FSR_segment_tallies = new int[_num_FSRs];
	Track* track;
	std::vector<segment*> segments;
	segment* segment;
	int num_segments;
	Cell* cell;

	/* Set each tally to zero to begin with */
	for (int i=0; i < _num_FSRs; i++)
		FSR_segment_tallies[i] = 0;


	/* Iterate over all azimuthal angles, all tracks, and all segments
	 * and tally each segment in the corresponding FSR */
	for (int i = 0; i < _num_azim; i++) {
		for (int j = 0; j < _num_tracks[i]; j++) {
			track = &_tracks[i][j];
			segments = track->getSegments();
			num_segments = track->getNumSegments();

			for (int s = 0; s < num_segments; s++) {
				segment = segments.at(s);
				FSR_segment_tallies[segment->_region_id]++;
			}
		}
	}


	/* Loop over all FSRs and if one FSR does not have tracks in it, print
	 * error message to the screen and exit program */
	for (int i=0; i < _num_FSRs; i++) {
		if (FSR_segment_tallies[i] == 0) {
			cell = _geom->findCell(i);
			log_printf(ERROR, "No tracks were tallied inside FSR id = %d which "
					"is cell id = %d. Please reduce your track spacing,"
					" increase the number of azimuthal angles, or increase the"
					" size of the flat source regions", i, cell->getId());
		}
	}

	delete [] FSR_segment_tallies;
}


/**
 * Compute the fission rates in each FSR and save them in a map of
 * FSR ids to fission rates
 */
void Solver::computePinPowers() {

	log_printf(NORMAL, "Computing pin powers...");

	FlatSourceRegion* fsr;
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

	/* make FSR BitMap */
	_plotter->makeFSRMap(bitMapFSR->pixels);

	/* Loop over all FSRs and compute the fission rate*/
	for (int i=0; i < _num_FSRs; i++) {
		fsr = &_FSRs[i];
		_FSRs_to_powers[i] = fsr->computeFissionRate();
	}

	/* Compute the pin powers by adding up the powers of FSRs in each
	 * lattice cell, saving lattice cell powers to files, and saving the
	 * pin power corresponding to each FSR id in FSR_to_pin_powers */
	_geom->computePinPowers(_FSRs_to_powers, _FSRs_to_pin_powers);


	/* Compute the total power based by accumulating the power of each unique
	 * pin with a nonzero power */
	for (int i=0; i < _num_FSRs; i++) {
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
	for (int i=0; i < _num_FSRs; i++) {
		_FSRs_to_pin_powers[i] /= avg_pin_power;
	}


	log_printf(NORMAL, "Plotting pin powers...");
	_plotter->makeRegionMap(bitMapFSR->pixels, bitMap->pixels,
											_FSRs_to_pin_powers);
	plot(bitMap, "pin_powers", _plotter->getExtension());

	/* delete bitMaps */
	deleteBitMap(bitMapFSR);
	deleteBitMap(bitMap);

	return;
}


/**
 * This method performs on or more fixed source iterations by integrating
 * the flux along each track and updating the boundary fluxes for the
 * corresponding output track, while updating the scalar flux in each
 * flat source region
 * @param max_iterations the maximum number of iterations allowed
 */
void Solver::fixedSourceIteration(int max_iterations) {

	Track* track;
	int num_segments;
	std::vector<segment*> segments;
	FP_PRECISION* weights;
	segment* segment;
	FP_PRECISION* polar_fluxes;
	FP_PRECISION* scalar_flux;
	FP_PRECISION* old_scalar_flux;
	FP_PRECISION* sigma_t;
	FlatSourceRegion* fsr;
	FP_PRECISION fsr_flux[NUM_ENERGY_GROUPS];
	FP_PRECISION* ratios;
	FP_PRECISION delta;
	FP_PRECISION volume;
	int t, j, k, s, p, e, pe;
	int max_num_threads = _num_azim / 2;

#if !STORE_PREFACTORS
	FP_PRECISION sigma_t_l;
	int index;
#endif

	log_printf(INFO, "Fixed source iteration with max_iterations = %d and "
			"# threads = %d", max_iterations, max_num_threads);

	/* Loop for until converged or max_iterations is reached */
	for (int i = 0; i < max_iterations; i++) {

		/* Initialize flux in each region to zero */
		zeroFSRFluxes();

		/* Loop over azimuthal each thread and azimuthal angle*
		 * If we are using OpenMP then we create a separate thread
		 * for each pair of reflecting azimuthal angles - angles which
		 * wrap into cycles on each other */
		/* Loop over each thread */
		#if STORE_PREFACTORS
		#pragma omp parallel for num_threads(std::min(omp_get_max_threads(), 16)) \
				private(t, k, j, s, p, e, pe, track, segments, \
						num_segments, weights, polar_fluxes, \
						segment, fsr, ratios, delta, fsr_flux)
		#else
		#pragma omp parallel for num_threads(std::min(omp_get_max_threads(), 16)) \
				private(t, k, j, s, p, e, pe, track, segments, \
						num_segments, weights, polar_fluxes,\
						segment, fsr, ratios, delta, fsr_flux,\
						sigma_t, sigma_t_l, index)
		#endif
		for (t=0; t < max_num_threads; t++) {

			/* Loop over the pair of azimuthal angles for this thread */
			j = t;
			while (j < _num_azim) {

			/* Loop over all tracks for this azimuthal angles */
			for (k = 0; k < _num_tracks[j]; k++) {

				/* Initialize local pointers to important data structures */
				track = &_tracks[j][k];
				segments = track->getSegments();
				num_segments = track->getNumSegments();
				weights = track->getPolarWeights();
				polar_fluxes = track->getPolarFluxes();

				/* Loop over each segment in forward direction */
				for (s = 0; s < num_segments; s++) {
					segment = segments.at(s);
					fsr = &_FSRs[segment->_region_id];

					ratios = fsr->getRatios();

					/* Zero out temporary FSR flux array */
					for (e = 0; e < NUM_ENERGY_GROUPS; e++)
						fsr_flux[e] = 0.0f;

					/* Initialize the polar angle and energy group counter */
					pe = 0;

#if !STORE_PREFACTORS
					sigma_t = segment->_material->getSigmaT();

					for (e = 0; e < NUM_ENERGY_GROUPS; e++) {

						fsr_flux[e] = 0.0f;
						sigma_t_l = sigma_t[e] * segment->_length;
						sigma_t_l = std::min(sigma_t_l, FP_PRECISION(10.0));
						index = sigma_t_l / _prefactor_spacing;
						index = std::min(index * 2 * NUM_POLAR_ANGLES,
												_prefactor_max_index);

						for (p = 0; p < NUM_POLAR_ANGLES; p++){
							delta = (polar_fluxes[pe] - ratios[e]) *
							(1 - (_prefactor_array[index + 2 * p] * sigma_t_l
							+ _prefactor_array[index + 2 * p + 1]));
							fsr_flux[e] += delta * weights[p];
							polar_fluxes[pe] -= delta;
							pe++;
						}
					}

#else
					/* Loop over all polar angles and energy groups */
					for (e = 0; e < NUM_ENERGY_GROUPS; e++) {
						for (p = 0; p < NUM_POLAR_ANGLES; p++) {
							delta = (polar_fluxes[pe] -ratios[e]) *
													segment->_prefactors[e][p];
							fsr_flux[e] += delta * weights[p];
							polar_fluxes[pe] -= delta;
							pe++;
						}
					}

#endif

					/* Increment the scalar flux for this FSR */
					fsr->incrementFlux(fsr_flux);
				}


				/* Transfer flux to outgoing track */
				track->getTrackOut()->setPolarFluxes(track->isReflOut(),
															0, polar_fluxes);

				/* Loop over each segment in reverse direction */
				for (s = num_segments-1; s > -1; s--) {
					segment = segments.at(s);
					fsr = &_FSRs[segment->_region_id];
					ratios = fsr->getRatios();

					/* Zero out temporary FSR flux array */
					for (e = 0; e < NUM_ENERGY_GROUPS; e++)
						fsr_flux[e] = 0.0;

					/* Initialize the polar angle and energy group counter */
					pe = GRP_TIMES_ANG;

#if !STORE_PREFACTORS
					sigma_t = segment->_material->getSigmaT();

					for (e = 0; e < NUM_ENERGY_GROUPS; e++) {

						fsr_flux[e] = 0;
						sigma_t_l = sigma_t[e] * segment->_length;
						sigma_t_l = std::min(sigma_t_l, FP_PRECISION(10.0));
						index = sigma_t_l / _prefactor_spacing;
						index = std::min(index * 2 * NUM_POLAR_ANGLES,
												_prefactor_max_index);

						for (p = 0; p < NUM_POLAR_ANGLES; p++){
							delta = (polar_fluxes[pe] - ratios[e]) *
							(1 - (_prefactor_array[index + 2 * p] * sigma_t_l
							+ _prefactor_array[index + 2 * p + 1]));
							fsr_flux[e] += delta * weights[p];
							polar_fluxes[pe] -= delta;
							pe++;
						}
					}

#else
					/* Loop over all polar angles and energy groups */
					for (e = 0; e < NUM_ENERGY_GROUPS; e++) {
						for (p = 0; p < NUM_POLAR_ANGLES; p++) {
							delta = (polar_fluxes[pe] - ratios[e]) *
											segment->_prefactors[e][p];
							fsr_flux[e] += delta * weights[p];
							polar_fluxes[pe] -= delta;
							pe++;
						}
					}
#endif

					/* Increment the scalar flux for this FSR */
					fsr->incrementFlux(fsr_flux);
				}

				/* Transfer flux to incoming track */
				track->getTrackIn()->setPolarFluxes(track->isReflIn(),
											GRP_TIMES_ANG, polar_fluxes);

			}

			/* Update the azimuthal angle index for this thread
			 * such that the next azimuthal angle is the one that reflects
			 * out of the current one. If instead this is the 2nd (final)
			 * angle to be used by this thread, break loop */
			if (j < max_num_threads)
				j = _num_azim - j - 1;
			else
				break;

			}
			
		}


		/* Add in source term and normalize flux to volume for each region */
		/* Loop over flat source regions, energy groups */
		#pragma omp parallel for private(fsr, scalar_flux, ratios, \
													sigma_t, volume)
		for (int r = 0; r < _num_FSRs; r++) {
			fsr = &_FSRs[r];
			scalar_flux = fsr->getFlux();
			ratios = fsr->getRatios();
			sigma_t = fsr->getMaterial()->getSigmaT();
			volume = fsr->getVolume();

			for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
				fsr->setFlux(e, scalar_flux[e] / 2.0);
				fsr->setFlux(e, FOUR_PI * ratios[e] + (scalar_flux[e] /
												(sigma_t[e] * volume)));
			}
		}


		/* Check for convergence if max_iterations > 1 */
		if (max_iterations > 1) {
			bool converged = true;

			#pragma omp parallel for private(fsr, scalar_flux, \
								old_scalar_flux) shared(converged)
			for (int r = 0; r < _num_FSRs; r++) {
				fsr = &_FSRs[r];
				scalar_flux = fsr->getFlux();
				old_scalar_flux = fsr->getOldFlux();

				for (int e = 0; e < NUM_ENERGY_GROUPS; e++) {
					if (fabs((scalar_flux[e] - old_scalar_flux[e]) /
							old_scalar_flux[e]) > FLUX_CONVERGENCE_THRESH)
						converged = false;

					/* Update old scalar flux */
					old_scalar_flux[e] = scalar_flux[e];
				}
			}

			if (converged)
				return;
		}

		/* Update the old scalar flux for each region, energy group */
		updateFSRFluxes();
	}

	if (max_iterations > 1)
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
FP_PRECISION Solver::computeKeff(int max_iterations) {

	FP_PRECISION* source;
	FlatSourceRegion* fsr;
	FP_PRECISION residual = 0.0;

	log_printf(NORMAL, "Computing k_eff on the host...");

	/* Check that each FSR has at least one segment crossing it */
	checkTrackSpacing();

	/* Set scalar flux to unity for each region */
    initFSRs();
	zeroTrackFluxes();

	/* Source iteration loop */
	for (int i = 0; i < max_iterations; i++) {

		log_printf(NORMAL, "Iteration %d on host: \t\tk_eff = %1.6f"
					"\t\tres = %1.3E", i, _k_eff, residual);

		/*********************************************************************
		 * Renormalize scalar and boundary fluxes
		 *********************************************************************/

		renormalizeFluxes();

		/*********************************************************************
		 * Compute the source for each region
		 *********************************************************************/

		residual = computeFSRSources();

		/*********************************************************************
		 * Update flux and check for convergence
		 *********************************************************************/

		/* Update pre-computed source / sigma_t ratios */
		computeRatios();

		/* Iteration the flux with the new source */
		fixedSourceIteration(1);

		/* Update k_eff */
		updateKeff();


		/* If k_eff converged, return k_eff */
		if (i > 1 && residual < SOURCE_CONVERG_THRESH){

			/* Converge the scalar flux spatially within geometry to plot */
			fixedSourceIteration(1000);

			if (_plotter->plotFlux() == true){
				/* Load fluxes into FSR to flux map */
				for (int r=0; r < _num_FSRs; r++) {
					FP_PRECISION* fluxes = _FSRs[r].getFlux();
					for (int e=0; e < NUM_ENERGY_GROUPS; e++){
						_FSRs_to_fluxes[e][r] = fluxes[e];
						_FSRs_to_fluxes[NUM_ENERGY_GROUPS][r] =
							_FSRs_to_fluxes[NUM_ENERGY_GROUPS][r] + fluxes[e];
					}
				}
				plotFluxes();
			}


			return _k_eff;
		}
	}

	log_printf(WARNING, "Unable to converge the source after %d iterations",
															max_iterations);

	return _k_eff;
}


// only plots flux
void Solver::plotFluxes(){

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
}
