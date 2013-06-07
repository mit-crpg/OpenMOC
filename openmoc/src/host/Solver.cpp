#include "Solver.h"


/**
 * @brief Constructor initializes empty arrays for source, flux, etc.
 * @details The construcor retrieves the number of energy groups and flat
 *          source regions and azimuthal angles from the geometry and track
 *          generator, and uses this to initialie empty arrays for the 
 *          flat source regions, boundary angular fluxes, scalar flatsourcergion
 *          fluxes, flatsourceregion sources and flatsourceregion powers. The 
 *          constructor initalizes the number of threads to a default of 1.
 * @param geometry an optional pointer to the geometry
 * @param track_generator an optional pointer to the trackgenerator
 * @param num_threads an optional number of threads
 */
Solver::Solver(Geometry* geometry, TrackGenerator* track_generator, 
	       int num_threads) {

    /* Default values */
    _num_groups = 0;
    _num_azim = 0;
    _num_FSRs = 0;
    _polar_times_groups = 0;

    _quad = NULL;
    _track_generator = NULL;
    _geometry = NULL;

    _tracks = NULL;
    _azim_weights = NULL;
    _polar_weights = NULL;
    _boundary_flux = NULL;

    _FSRs = NULL;
    _scalar_flux = NULL;
    _old_scalar_flux = NULL;
    _source = NULL;
    _old_source = NULL;
    _ratios = NULL;
    _FSR_to_power = NULL;
    _FSR_to_pin_power = NULL;

    _prefactor_array = NULL;

    if (geometry != NULL)
        setGeometry(geometry);
    if (track_generator != NULL)
        setTrackGenerator(track_generator);

    setNumThreads(num_threads);

    /* Default polar quadrature */
    _quadrature_type = TABUCHI;
    _num_polar = 3;
    _two_times_num_polar = 2 * _num_polar;

    _num_iterations = 0;
    _source_convergence_thresh = 1E-3;
    _flux_convergence_thresh = 1E-5;
    _converged_source = false;
}


/**
 * @brief Destructor deletes arrays of boundary angular flux for all tracks,
 *        scalar flux and source for each flatsourceregion.
 */
Solver::~Solver() {

    if (_polar_weights != NULL)
        delete [] _polar_weights;
    if (_boundary_flux != NULL)
        delete [] _boundary_flux;
    if (_FSRs != NULL)
        delete [] _FSRs;
    if (_scalar_flux != NULL)
        delete [] _scalar_flux;
    if (_old_scalar_flux != NULL)
        delete [] _old_scalar_flux;
    if (_source != NULL)
        delete [] _source;
    if (_old_source != NULL)
        delete [] _old_source;
    if (_ratios != NULL)
        delete [] _ratios;
    if (_FSR_to_power != NULL)
        delete [] _FSR_to_power;
    if (_FSR_to_pin_power != NULL)
        delete [] _FSR_to_pin_power;
    if (_prefactor_array != NULL)
        delete [] _prefactor_array;
    if (_quad != NULL)
        delete _quad;
}


/**
 * @brief Returns a pointer to the geometry for this solver.
 * @return a pointer to the geometry
 */
Geometry* Solver::getGeometry() {
    if (_geometry == NULL)
        log_printf(ERROR, "Unable to return the solver's geometry since it "
		 "has not yet been set");

    return _geometry;
}


/**
 * @brief Returns a pointer to the geometry for this solver.
 * @return a pointer to the geometry
 */
TrackGenerator* Solver::getTrackGenerator() {
    if (_track_generator == NULL)
        log_printf(ERROR, "Unable to return the solver's track genetrator "
		   "since it has not yet been set");

    return _track_generator;
}


/**
 * @brief Returns the number of shared memory OpenMP threads in use.
 * @return the number of threads
 */
int Solver::getNumThreads() {
    return _num_threads;
}


/**
 * @brief Returns the number of angles used for the polar quadrature.
 * @return the number of polar angles
 */
int Solver::getNumPolarAngles() {
    return _num_polar;
}


/**
 * @brief Returns the type of polar quadrature in use (TABUCHI or LEONARD).
 * @return the type of polar quadrature
 */
quadratureType Solver::getPolarQuadratureType() {
    return _quadrature_type;
}


/**
 * @brief Returns the number of transport sweeps to converge the source.
 * @return the number of iterations
 */
int Solver::getNumIterations() {
    return _num_iterations;
}


/**
 * @brief Returns the threshold for source convergence.
 * @return the threshold for source convergence
 */
FP_PRECISION Solver::getSourceConvergenceThreshold() {
    return _source_convergence_thresh;
}


/**
 * @brief Returns the threshold for flux convergence in fixed source iteration
 *        after the source has converged.
 * @return the threshold for flux convergence
 */
FP_PRECISION Solver::getFluxConvergenceThreshold() {
    return _flux_convergence_thresh;
}


FP_PRECISION Solver::getFSRScalarFlux(int fsr_id, int energy_group) {

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
FP_PRECISION* Solver::getFSRScalarFluxes() {
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
FP_PRECISION* Solver::getFSRPowers() {
    if (_FSR_to_power == NULL)
        log_printf(ERROR, "Unable to returns the Solver's FSR power array "
		 "since it has not yet been allocated in memory");

    return _FSR_to_power;
}


/**
 * @brief Return an array indexed by flatsourceregion IDs with the
 *        corresponding pin cell power.
 * @return an array of flatsourceregion pin powers
 */
FP_PRECISION* Solver::getFSRPinPowers() {
    if (_FSR_to_pin_power == NULL)
        log_printf(ERROR, "Unable to returns the Solver's FSR pin power array "
		 "since it has not yet been allocated in memory");

    return _FSR_to_pin_power;
}


/**
 * @brief Sets the number of shared memory OpenMP threads to use (>0).
 * @details This method sets the number of threads to be no larger than the
 *          number of complementary angle pairs - reflecting angles which form
 *          cyclical tracks - since that is how OpenMOC parallizes loops.
 * @param num_threads the number of threads
 */
void Solver::setNumThreads(int num_threads) {
    if (num_threads <= 0)
        log_printf(ERROR, "Unable to set the number of threads for the Solver "
		   "to %d since it is less than or equal to 0", num_threads);

    /* Set the number of threads to be no larger than half the number of 
     * azimuthal angles if the trackgenerator has been set */
    if (_num_azim == 0 || num_threads < _num_azim / 2)
      _num_threads = num_threads;
    else
        _num_threads = _num_azim / 2;

    /* Set the number of threads for OpenMP */
    omp_set_num_threads(_num_threads);
}


/**
 * @brief Sets the geometry for the solver.
 * @details The geometry must already have initialized flat source region maps
 *          and segmentized the trackgenerator's tracks.
 * @param geometry a pointer to a geometry
 */
void Solver::setGeometry(Geometry* geometry) {
    if (geometry->getNumFSRs() == 0)
        log_printf(ERROR, "Unable to set the Geometry for the Solver "
		 "since the Geometry has not yet initialized flat "
		 "source regions");
    if (geometry->getNumEnergyGroups() == 0)
        log_printf(ERROR, "Unable to set the Geometry for the Solver "
		 "since the Geometry does noet contain any materials");

    _geometry = geometry;
    _num_FSRs = _geometry->getNumFSRs();
    _num_groups = _geometry->getNumEnergyGroups();
    _polar_times_groups = _num_groups * _num_polar;
}


/**
 * @brief Sets the trackgenerator with characteristic tracks for the solver.
 * @details The trackgenerator must already have generated tracks and have
 *          segmentized them using the geometry.
 * @param track_generator a pointer to a trackgenerator
 */
void Solver::setTrackGenerator(TrackGenerator* track_generator) {
    if (!track_generator->containsTracks())
        log_printf(ERROR, "Unable to set the TrackGenerator for the Solver "
		 "since the TrackGenerator has not yet generated tracks");

    _track_generator = track_generator;
    _num_azim = _track_generator->getNumAzim() / 2;
    _tracks = _track_generator->getTracks();
    _num_tracks = _track_generator->getNumTracksArray();
    _azim_weights = _track_generator->getAzimWeights();
}


/**
 * @brief Sets the type of polar angle quadrature set to use (ie, TABUCHI 
 *        or LEONARD).
 * @param type the polar angle quadrature type
 */
void Solver::setPolarQuadratureType(quadratureType quadrature_type) {
    _quadrature_type = quadrature_type;
}


/**
 * @brief Sets the number of polar angles to use (only 1, 2, or 3 currently
 *        supported).
 * @param num_polar the number of polar angles
 */
void Solver::setNumPolarAngles(int num_polar) {

    if (num_polar <= 0)
        log_printf(ERROR, "Unable to set the Solver's number of polar angles "
		 "to %d since this is a negative number", num_polar);

    if (num_polar > 3)
        log_printf(ERROR, "Unable to set the Solver's number of polar angles "
		   "to %d since this is not a supported value (only 1, 2 or 3 "
		   " are currently supported)", num_polar);

    _num_polar = num_polar;
    _two_times_num_polar = 2 * _num_polar;
    _polar_times_groups = _num_groups * _num_polar;
}


/**
 * @brief Sets the threshold for source convergence (>0)
 * @param source_thresh the threshold for source convergence
 */
void Solver::setSourceConvergenceThreshold(FP_PRECISION source_thresh) {
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
void Solver::setFluxConvergenceThreshold(FP_PRECISION flux_thresh) {
    if (flux_thresh <= 0.0)
        log_printf(ERROR, "Unable to set the flux convergence threshold to "
	       "%f since the threshold must be a positive number",
	       flux_thresh);

    _flux_convergence_thresh = flux_thresh;
}


/**
 * @brief Allocates memory for track boundary angular fluxes and 
 *        flatsourceregion scalar fluxes.
 * @details Deletes memory for old flux arrays if they were allocated from
 *          previous simulation.
 */
void Solver::initializeFluxArrays() {

    /* Delete old flux arrays if they exist */
    if (_boundary_flux != NULL)
        delete [] _boundary_flux;
    if (_scalar_flux != NULL)
        delete [] _scalar_flux;
    if (_old_scalar_flux != NULL)
        delete [] _old_scalar_flux;

    int tot_num_tracks = _track_generator->getNumTracks();

    /* Allocate memory for all flux arrays */
    try{
        _boundary_flux = new FP_PRECISION[2*tot_num_tracks*_polar_times_groups];
	_scalar_flux = new FP_PRECISION[_num_FSRs*_num_groups];
	_old_scalar_flux = new FP_PRECISION[_num_FSRs*_num_groups];
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
void Solver::initializeSourceArrays() {

    /* Delete old sources arrays if they exist */
    if (_source != NULL)
        delete [] _source;
    if (_old_source != NULL)
        delete [] _old_source;
    if (_ratios != NULL)
        delete [] _ratios;

    /* Allocate memory for all source arrays */
    try{
	_source = new FP_PRECISION[_num_FSRs*_num_groups];
	_old_source = new FP_PRECISION[_num_FSRs*_num_groups];
	_ratios = new FP_PRECISION[_num_FSRs*_num_groups];
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
void Solver::initializePowerArrays() {

    /* Delete old power arrays if they exist */
    if (_FSR_to_power != NULL)
        delete [] _FSR_to_power;
    if (_FSR_to_pin_power != NULL)
        delete [] _FSR_to_pin_power;

    /* Allocate memory for FSR power and pin power arrays */
    try{
	_FSR_to_power = new FP_PRECISION[_num_FSRs];
	_FSR_to_pin_power = new FP_PRECISION[_num_FSRs];
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's FSR "
		   "power arrays. Backtrace:%s", e.what());
    }
}


/**
 * @brief Creates a polar quadrature object for the solver.
 */
void Solver::initializePolarQuadrature() {
    /* Deletes the old quadrature if one existed */
    if (_quad != NULL)
        delete _quad;

    _quad = new Quadrature(_quadrature_type, _num_polar);
    _polar_times_groups = _num_groups * _num_polar;
}

 
/**
 * @brief Pre-computes exponential pre-factors for each segment of each track 
 *        for each polar angle. 
 * @details This method will generate a hashmap which contains values of the 
 *          pre-factor for specific segment lengths (the keys into the hashmap).
 */
void Solver::precomputePrefactors() {

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
void Solver::initializeFSRs() {

    log_printf(INFO, "Initializing flat source regions...");

    /* Delete old FSRs array if it exists */
    if (_FSRs != NULL)
        delete [] _FSRs;

    _FSRs = new FlatSourceRegion[_num_FSRs];

    std::vector<segment*> segments;
    std::vector<segment*>::iterator iter;
    FP_PRECISION volume;
    CellBasic* cell;
    Material* material;
    Universe* univ_zero = _geometry->getUniverse(0);

    /* Set each FSR's volume by accumulating the total length of all tracks
     * inside the FSR. Loop over azimuthal angle, track and segment. 
     * Note: this code region cannot be parallelized without a mutex lock
     * on FSR volume due to race conditions. */
    for (int i=0; i < _num_azim; i++) {
        for (int j=0; j < _num_tracks[i]; j++) {
  	    segments = _tracks[i][j].getSegments();

            for (iter=segments.begin(); iter != segments.end(); ++iter) {
	        volume = (*iter)->_length * _azim_weights[i];
		_FSRs[(*iter)->_region_id].incrementVolume(volume);
	    }
	}
    }

    /* Loop over all FSRs */
    #pragma omp parallel for private(cell, material)
    for (int r=0; r < _num_FSRs; r++) {

        /* Get the cell corresponding to this FSR from the geometry */
        cell = static_cast<CellBasic*>(_geometry->findCell(univ_zero, r));

	/* Get the cell's material and assign it to the FSR */
	material = _geometry->getMaterial(cell->getMaterial());
	_FSRs[r].setMaterial(material);

	log_printf(DEBUG, "FSR id = %d has cell id = %d and material id = %d "
		   "and volume = %f", r, cell->getId(), material->getId(),
		   _FSRs[r].getVolume());
    }

    return;
}


/**
 * @brief Checks that each flat source region has at least one segment within 
 *        it and if not, throw an exception and prints an error message.
 */
void Solver::checkTrackSpacing() {

    int* FSR_segment_tallies = new int[_num_FSRs];
    std::vector<segment*> segments;
    std::vector<segment*>::iterator iter;
    Cell* cell;

    /* Set each tally to zero to begin with */
    #pragma omp parallel for
    for (int r=0; r < _num_FSRs; r++)
        FSR_segment_tallies[r] = 0;

    /* Iterate over all azimuthal angles, all tracks, and all segments
     * and tally each segment in the corresponding FSR */
    #pragma omp parallel for private (segments, iter)
    for (int i=0; i < _num_azim; i++) {
        for (int j=0; j < _num_tracks[i]; j++) {
	    segments = _tracks[i][j].getSegments();

            for (iter=segments.begin(); iter != segments.end(); ++iter)
	        FSR_segment_tallies[(*iter)->_region_id]++;
	}
    }

    /* Loop over all FSRs and if one FSR does not have tracks in it, print
     * error message to the screen and exit program */
    #pragma omp parallel for private (cell)
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


/**
 * @brief Zero each track's boundary fluxes for each energy group and polar
 *        angle in the "forward" and "reverse" directions.
 */
void Solver::zeroTrackFluxes() {

    /* Loop over azimuthal angle, track, polar angle, energy group
     * and set each track's incoming and outgoing flux to zero */
    #pragma omp parallel for
    for (int i=0; i < _track_generator->getNumTracks(); i++) {
        for (int pe2=0; pe2 < 2*_polar_times_groups; pe2++)
    	    _boundary_flux(i,pe2) = 0.0;
    }

    return;
}


/**
 * @brief Set the scalar flux for each energy group inside each 
 *        flatsourceregion to a constant value.
 * @param value the value to assign to each flat source region flux
 */
void Solver::flattenFSRFluxes(FP_PRECISION value) {

    /* Loop over all FSRs and energy groups */
    #pragma omp parallel for
    for (int r=0; r < _num_FSRs; r++) {
        for (int e=0; e < _num_groups; e++) {
            _scalar_flux(r,e) = value;
  	    _old_scalar_flux(r,e) = value;
         }
     }

    return;
}


/**
 * @brief Set the source for each energy group inside each flatsourceregion
 *        to a constant value.
 * @param value the value to assign to each flat source region source
 */
void Solver::flattenFSRSources(FP_PRECISION value) {

    /* Loop over all FSRs and energy groups */
    #pragma omp parallel for
    for (int r=0; r < _num_FSRs; r++) {
        for (int e=0; e < _num_groups; e++) {
            _source(r,e) = value;
    	    _old_source(r,e) = value;
        }
    }

    return;
}


/**
 * @brief Normalizes all flatsourceregion scalar fluxes and track boundary
 *        angular fluxes to the total fission source (times nu).
 */
void Solver::normalizeFluxes() {

    double* nu_sigma_f;
    FP_PRECISION volume;
    FP_PRECISION fission_source;
    FP_PRECISION norm_factor;

    /* Initialize fission source to zero */
    fission_source = 0;

    /* Compute total fission source for this region */
    #pragma omp parallel for private(volume, nu_sigma_f) \
      reduction(+: fission_source)
    for (int r=0; r < _num_FSRs; r++) {

        /* Get pointers to important data structures */
	nu_sigma_f = _FSRs[r].getMaterial()->getNuSigmaF();
	volume = _FSRs[r].getVolume();

	for (int e=0; e < _num_groups; e++)
	    fission_source += nu_sigma_f[e] * _scalar_flux(r,e) * volume;
    }

    /* Normalize scalar fluxes in each region */
    norm_factor = 1.0 / fission_source;

    log_printf(DEBUG, "Normalization factor = %f", norm_factor);

    #pragma omp parallel for
    for (int r=0; r < _num_FSRs; r++) {
        for (int e=0; e < _num_groups; e++)
	  _scalar_flux(r,e) *= norm_factor;
    }

    /* Normalize angular boundary fluxes for each track */
    #pragma omp parallel for
    for (int i=0; i < _track_generator->getNumTracks(); i++) {
        for (int pe2=0; pe2 < 2*_polar_times_groups; pe2++)
	    _boundary_flux(i,pe2) *= norm_factor;
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
FP_PRECISION Solver::computeFSRSources() {

    FP_PRECISION scatter_source;
    FP_PRECISION fission_source;
    double* nu_sigma_f;
    double* sigma_s;
    double* sigma_t;
    double* chi;
    Material* material;
    FP_PRECISION source_residual = 0.0;

    /* For all regions, find the source */
    //TODO: This can be parallelized! Need to privatize some variable and
    //      reduce the residual at the end
    #pragma omp parallel for private(material, nu_sigma_f, chi, sigma_s, \
      sigma_t, fission_source, scatter_source) \
      reduction(+:source_residual)
    for (int r=0; r < _num_FSRs; r++) {

        material = _FSRs[r].getMaterial();
	nu_sigma_f = material->getNuSigmaF();
	chi = material->getChi();
	sigma_s = material->getSigmaS();
        sigma_t = material->getSigmaT();

	/* Initialize the fission source to zero for this region */
	fission_source = 0;

	/* Compute total fission source for current region */
	for (int e=0; e < _num_groups; e++)
  	    fission_source += _scalar_flux(r,e) * nu_sigma_f[e];

	/* Compute total scattering source for group G */
        for (int G=0; G < _num_groups; G++) {
            scatter_source = 0;

	    for (int g=0; g < _num_groups; g++)
	        scatter_source += sigma_s[G*_num_groups+g] * _scalar_flux(r,g);

	    /* Set the total source for region r in group G */
	    _source(r,G) = ((1.0 / _k_eff) * fission_source *
                           chi[G] + scatter_source) * ONE_OVER_FOUR_PI;
	
	    _ratios(r,G) = _source(r,G) / sigma_t[G];


	    /* Compute the norm of residuals of the sources for convergence */
	    if (fabs(_source(r,G)) > 1E-10)
	        source_residual += pow((_source(r,G) - _old_source(r,G)) 
                                       / _source(r,G), 2);
	    
	    /* Update the old source */
	    _old_source(r,G) = _source(r,G);
        }
    }

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
void Solver::computeKeff() {

    Material* material;
    double* sigma_a;
    double* nu_sigma_f;
    FP_PRECISION volume;

    FP_PRECISION tot_abs = 0.0;
    FP_PRECISION tot_fission = 0.0;

    #pragma omp parallel for private(volume, material, sigma_a, nu_sigma_f) \
      reduction(+:tot_abs) reduction(+:tot_fission)
    for (int r=0; r < _num_FSRs; r++) {

        volume = _FSRs[r].getVolume();
	material = _FSRs[r].getMaterial();
	sigma_a = material->getSigmaA();
	nu_sigma_f = material->getNuSigmaF();

	for (int e=0; e < _num_groups; e++) {
            tot_abs += sigma_a[e] * _scalar_flux(r,e) * volume;
	    tot_fission += nu_sigma_f[e] * _scalar_flux(r,e) * volume;
        }
    }

    log_printf(NORMAL, "abs = %f, fiss = %f, leak = %f", tot_abs, tot_fission, _leakage);

    _k_eff = tot_fission / (tot_abs + _leakage);

    log_printf(DEBUG, "tot_abs = %f, tot_fission = %f, leakage = %f, "
	       "k_eff = %f", tot_abs, tot_fission, _leakage, _k_eff);

    return;
}


/**
 * @brief Checks if scalar flux has converged within the threshold.
 * @return true if converged, false otherwise
 */
bool Solver::isScalarFluxConverged() {

    bool converged = true;

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

    return converged;
}


/**
 * This method performs on or more fixed source iterations by integrating
 * the flux along each track and updating the boundary fluxes for the
 * corresponding output track, while updating the scalar flux in each
 * flat source region
 * @param max_iterations the maximum number of iterations allowed
 */
void Solver::transportSweep(int max_iterations) {

    bool converged = false;
    int thread_id;
    int track_id;
    int track_out_id;
    bool bc;
    int fsr_id;
    std::vector<segment*> segments;
    std::vector<segment*>::iterator iter;
    std::vector<segment*>::reverse_iterator riter;
    double* sigma_t;
    FP_PRECISION fsr_flux;
    FP_PRECISION delta;
    FP_PRECISION volume;
    FP_PRECISION sigma_t_l;
    int index;
    int pe;
    int max_num_threads = _num_azim / 2;

    /* Allocate memory for each thread's FSR scalar fluxes and leakages */
    FP_PRECISION* thread_flux = new FP_PRECISION[_num_threads * _num_FSRs * 
                                                                 _num_groups];
    FP_PRECISION* thread_leakage = new FP_PRECISION[_num_threads];

    /* Initialize thread fluxes and leakages to zero */
    memset(thread_flux, FP_PRECISION(0.), 
	   _num_threads*_num_FSRs*_num_groups*sizeof(FP_PRECISION));
    memset(thread_leakage, FP_PRECISION(0.), _num_threads*sizeof(FP_PRECISION));

    int start;

    log_printf(DEBUG, "Transport sweep with max_iterations = %d and "
	       "# threads = %d", max_iterations, _num_threads);

    /* Loop for until converged or max_iterations is reached */
    for (int i=0; i < max_iterations; i++) {

        /* Initialize the global total leakage tally to zero */
        _leakage = 0.0;

        /* Initialize flux in each region to zero */
        flattenFSRFluxes(0.0);

	/* Loop over each thread and azimuthal angle.
	 * If we are using more than 1 thread then we create 
	 * separate threads for each pair of complementary  
	 * azimuthal angles - angles which wrap into cycles */
	/* Loop over each thread */
        #pragma omp parallel for private(track_id, track_out_id, bc, fsr_id, \
	          segments, iter, riter, sigma_t, fsr_flux, delta, volume, \
	          pe, sigma_t_l, index, thread_id, start)
	for (int t=0; t < max_num_threads; t++) {

            thread_id = omp_get_thread_num();

            /* Loop over the pair of azimuthal angles for this thread */
	    int j = t;
	    while (j < _num_azim) {

	        /* Loop over all tracks for this azimuthal angles */
	        for (int k=0; k < _num_tracks[j]; k++) {

		    /* Initialize local pointers to important data structures */
  		    track_id = _tracks[j][k].getUid();
		    segments = _tracks[j][k].getSegments();

		    /* Loop over each segment in forward direction */
                    for (iter=segments.begin(); iter!=segments.end(); ++iter) {
                        fsr_id = (*iter)->_region_id;

			/* Initialize polar angle and energy group counter */
			pe = 0;
			sigma_t = (*iter)->_material->getSigmaT();

			/* Loop over energy groups */
			for (int e=0; e < _num_groups; e++) {
			    fsr_flux = 0.;
			    sigma_t_l = sigma_t[e] * (*iter)->_length;
			    index = computePrefactorIndex(sigma_t_l);
			    
			    /* Loop over polar angles */
			    for (int p=0; p < _num_polar; p++){
			        delta = (_boundary_flux(track_id,pe) - 
					_ratios(fsr_id,e)) * 
				        prefactor(index,p,sigma_t_l);
				fsr_flux += delta * _polar_weights[p];
				_boundary_flux(track_id,pe) -= delta;

				if (track_id == 1 && p == 0 && e == 0) {
				    printf("delta =%f, fsr_flux = %f, boundary flux = %f\n",
					     delta, fsr_flux, _boundary_flux(track_id,pe));
				    printf("polar weights = %f, ratios = %f, prefactor = %f\n", 
					   _polar_weights[p], _ratios(fsr_id,e),
					   prefactor(index,p,sigma_t_l));
				}

				pe++;
			    }

			    /* Increment the scalar flux for this thread' copy
			     * of this flat source region */
  		            thread_flux(thread_id,fsr_id,e) += fsr_flux;
			}			    
		    }

		    /* Transfer flux to outgoing track */
		    track_out_id = _tracks[j][k].getTrackOut()->getUid();
		    bc = _tracks[j][k].getBCOut();
		    start = _tracks[j][k].isReflOut() * _polar_times_groups;
		    for (pe=0; pe < _polar_times_groups; pe++) {
		        _boundary_flux(track_out_id,start+pe) = 
			    _boundary_flux(track_id,pe) * bc;
		        thread_leakage[thread_id] += 
			    _boundary_flux(track_id,pe) 
			    * _polar_weights[pe%_num_polar] * (!bc);
		    }


		    /* Loop over each segment in reverse direction */
		    for (riter=segments.rbegin(); riter != segments.rend(); 
                                                                      ++riter){

                        fsr_id = (*riter)->_region_id;

			/* Initialize polar angle and energy group counter */
			pe = _polar_times_groups;
			sigma_t = (*riter)->_material->getSigmaT();
			
			/* Loop over energy groups */
			for (int e=0; e < _num_groups; e++) {
			    fsr_flux = 0.;
			    sigma_t_l = sigma_t[e] * (*riter)->_length;
			    index = computePrefactorIndex(sigma_t_l);

			    /* Loop over polar angles */
			    for (int p=0; p < _num_polar; p++){
			        delta = (_boundary_flux(track_id,pe) - 
					 _ratios(fsr_id,e)) * 
				         prefactor(index,p,sigma_t_l);
				fsr_flux += delta * _polar_weights[p];
				_boundary_flux(track_id,pe) -= delta;

				if (track_id == 1 && p == 0 && e == 0) {
				    printf("delta =%f, fsr_flux = %f, boundary flux = %f\n",
					     delta, fsr_flux, _boundary_flux(track_id,pe));
				    printf("polar weights = %f, ratio = %f, prefactor = %f\n", 
					   _polar_weights[p], _ratios(fsr_id,e), 
					   prefactor(index,p,sigma_t_l));
				}

				pe++;
			    }

			    /* Increment the scalar flux for this thread' copy
			     * of this flat source region */
  		            thread_flux(thread_id,fsr_id,e) += fsr_flux;
			}
		    }

		    /* Transfer flux to outgoing track */
		    track_out_id = _tracks[j][k].getTrackIn()->getUid();
		    bc = _tracks[j][k].getBCIn();
		    start = _tracks[j][k].isReflIn() * _polar_times_groups;
		    for (pe=0; pe < _polar_times_groups; pe++) {
		        _boundary_flux(track_out_id,start+pe) = 
			   _boundary_flux(track_id,_polar_times_groups+pe) * bc;
		        thread_leakage[thread_id] += 
			    _boundary_flux(track_id,_polar_times_groups+pe) 
			    * _polar_weights[pe%_num_polar] * (!bc);
		    }
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

        /** Reduce leakage across threads */
        for (int t=0; t < _num_threads; t++)
	    _leakage += thread_leakage[t] * 0.5;

	/* Reduce scalar fluxes across threads from transport sweep and
	 * add in source term and normalize flux to volume for each region */
	/* Loop over flat source regions, energy groups */
        #pragma omp parallel for private(volume, sigma_t)
	for (int r=0; r < _num_FSRs; r++) {
	    volume = _FSRs[r].getVolume();
	    sigma_t = _FSRs[r].getMaterial()->getSigmaT();

	    for (int e=0; e < _num_groups; e++) {
	      
	        /* Reduce flux across threads from transport sweep */
                for (int t=0; t < _num_threads; t++)
		    _scalar_flux(r,e) += thread_flux(t,r,e);
	            _scalar_flux(r,e) *= 0.5;
		    _scalar_flux(r,e) = FOUR_PI * _ratios(r,e) + 
		      (_scalar_flux(r,e) / (sigma_t[e] * volume));
	    }
	}

	/* Check for convergence if max_iterations > 1 */
	if (max_iterations == 1 || isScalarFluxConverged()) {
            delete [] thread_flux;
            delete [] thread_leakage;
            return;
	}
    }

    log_printf(WARNING, "Scalar flux did not converge after %d iterations",
	                                                   max_iterations);

    delete [] thread_flux;
    delete [] thread_leakage;
    return;
}


/**
 * Computes keff on the by performing a series of fixed source
 * iterations and updating the fission and scattering sources in each
 * flat source region of the geometry
 * @param max_iterations the maximum number of iterations allowed
 * @return the value of keff computed
 */
FP_PRECISION Solver::convergeSource(int max_iterations) {

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

        log_printf(NORMAL, "Iteration %d on host: \tk_eff = %1.6f"
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


/**
 * @brief Compute the fission rates in each flatsourceregion and stores them 
 *        in an array indexed by flatsourceregion ID.
 */
void Solver::computePinPowers() {

    log_printf(INFO, "Computing FSR pin powers...");

    double* sigma_f;
    FP_PRECISION tot_pin_power = 0.;
    FP_PRECISION avg_pin_power = 0.;
    FP_PRECISION num_nonzero_pins = 0.;
    FP_PRECISION curr_pin_power = 0.;
    FP_PRECISION prev_pin_power = 0.;

    /* Loop over all FSRs and compute the fission rate*/
    #pragma omp parallel for private (sigma_f)
    for (int r=0; r < _num_FSRs; r++) {
        sigma_f = _FSRs[r].getMaterial()->getSigmaF();

        for (int e=0; e < _num_groups; e++)
	    _FSR_to_power[r] += sigma_f[e] * _scalar_flux(r,e);
    }

    /* Compute the pin powers by adding up the powers of FSRs in each
     * lattice cell, saving lattice cell powers to files, and saving the
     * pin power corresponding to each FSR id in FS_to_pin_power */
    _geometry->computePinPowers(_FSR_to_power, _FSR_to_pin_power);


    /* Compute the total power based by accumulating the power of each unique
     * pin with a nonzero power */
    for (int r=0; r < _num_FSRs; r++) {
        curr_pin_power = _FSR_to_pin_power[r];

	/* If this pin power is unique and nozero (doesn't match the previous
	 * pin's power), then tally it */
	if (curr_pin_power > 0. && curr_pin_power != prev_pin_power) {
	    tot_pin_power += curr_pin_power;
	    num_nonzero_pins++;
	    prev_pin_power = curr_pin_power;
	}
    }

    /* Compute the average pin power */
    avg_pin_power = tot_pin_power / num_nonzero_pins;

    /* Normalize each pin power to the average non-zero pin power */
    #pragma omp parallel for
    for (int r=0; r < _num_FSRs; r++)
        _FSR_to_pin_power[r] /= avg_pin_power;

    return;
}
