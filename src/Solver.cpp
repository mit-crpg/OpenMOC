#include "Solver.h"


/**
 * @brief Constructor initializes array pointers for tracks and materials.
 * @details The construcor retrieves the number of energy groups and flat
 *          source regions and azimuthal angles from the geometry and track
 *          generator. The constructor initalizes the number of threads to a 
 *          default of 1.
 * @param geometry an optional pointer to the geometry
 * @param track_generator an optional pointer to the trackgenerator
 */
Solver::Solver(Geometry* geometry, TrackGenerator* track_generator) {

    /* Default values */
    _num_materials = 0;
    _num_groups = 0;
    _num_azim = 0;
    _polar_times_groups = 0;

    _num_FSRs = 0;
    _FSR_volumes = NULL;
    _FSR_materials = NULL;

    _quad = NULL;
    _track_generator = NULL;
    _geometry = NULL;

    _tracks = NULL;
    _azim_weights = NULL;
    _polar_weights = NULL;
    _boundary_flux = NULL;
    _boundary_leakage = NULL;

    _scalar_flux = NULL;
    _fission_sources = NULL;
    _scatter_sources = NULL;
    _source = NULL;
    _old_source = NULL;
    _reduced_source = NULL;
    _source_residuals = NULL;

    _FSRs_to_powers = NULL;
    _FSRs_to_pin_powers = NULL;

    _prefactor_array = NULL;

    if (geometry != NULL)
        setGeometry(geometry);

    if (track_generator != NULL)
        setTrackGenerator(track_generator);

    /* Default polar quadrature */
    _quadrature_type = TABUCHI;
    _num_polar = 3;
    _two_times_num_polar = 2 * _num_polar;

    _num_iterations = 0;
    _source_convergence_thresh = 1E-3;
    _converged_source = false;

    _timer = new Timer();
}


/**
 * @brief Destructor deletes arrays of boundary angular fluxes,
 *        scalar fluxes and sources for each flat source region.
 * @details Deallocates memory for all arrays allocated for the Solver,
 *          including fluxes, sources, quadrature weights, and exponential
 *          prefactor interpolation table.
 */
Solver::~Solver() {

    if (_FSR_volumes != NULL)
        delete [] _FSR_volumes;

    if (_FSR_materials != NULL)
        delete [] _FSR_materials;

    if (_polar_weights != NULL)
        delete [] _polar_weights;

    if (_boundary_flux != NULL)
        delete [] _boundary_flux;

    if (_scalar_flux != NULL)
        delete [] _scalar_flux;

    if (_fission_sources != NULL)
        delete [] _fission_sources;

    if (_scatter_sources != NULL)
        delete [] _scatter_sources;

    if (_source != NULL)
        delete [] _source;

    if (_old_source != NULL)
        delete [] _old_source;

    if (_reduced_source != NULL)
        delete [] _reduced_source;

    if (_source_residuals != NULL)
        delete [] _source_residuals;

    if (_FSRs_to_powers != NULL)
        delete [] _FSRs_to_powers;

    if (_FSRs_to_pin_powers != NULL)
        delete [] _FSRs_to_pin_powers;

    if (_prefactor_array != NULL)
        delete [] _prefactor_array;

    if (_quad != NULL)
        delete _quad;
}


/**
 * @brief Returns a pointer to the geometry.
 * @return a pointer to the geometry
 */
Geometry* Solver::getGeometry() {

    if (_geometry == NULL)
        log_printf(ERROR, "Unable to return the solver's geometry since it "
		 "has not yet been set");

    return _geometry;
}


/**
 * @brief Returns a pointer to the track generator.
 * @return a pointer to the geometry
 */
TrackGenerator* Solver::getTrackGenerator() {

    if (_track_generator == NULL)
        log_printf(ERROR, "Unable to return the solver's track genetrator "
		   "since it has not yet been set");

    return _track_generator;
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
 * @brief Sets the geometry for the solver.
 * @details The geometry must already have initialized flat source region maps
 *          and segmentized the track generator's tracks.
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
    _num_materials = _geometry->getNumMaterials();
}


/**
 * @brief Sets the track generator with characteristic tracks for the solver.
 * @details The track generator must already have generated tracks and have
 *          segmentized them across the geometry.
 * @param track_generator a pointer to a trackgenerator
 */
void Solver::setTrackGenerator(TrackGenerator* track_generator) {

    if (!track_generator->containsTracks())
        log_printf(ERROR, "Unable to set the TrackGenerator for the Solver "
		 "since the TrackGenerator has not yet generated tracks");

    _track_generator = track_generator;
    _num_azim = _track_generator->getNumAzim() / 2;
    _num_tracks = _track_generator->getNumTracksArray();
    _tot_num_tracks = _track_generator->getNumTracks();
    _azim_weights = _track_generator->getAzimWeights();
    _tracks = new Track*[_tot_num_tracks];

    /* Initialize the tracks array */
    int counter = 0;

    for (int i=0; i < _num_azim; i++) {
        for (int j=0; j < _num_tracks[i]; j++) {
	  _tracks[counter] = &_track_generator->getTracks()[i][j];
	  counter++;
	}
    }
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
 * @brief Checks that each flat source region has at least one segment within 
 *        it and if not, throws an exception and prints an error message.
 */
void Solver::checkTrackSpacing() {

    int* FSR_segment_tallies = new int[_num_FSRs];
    int num_segments;
    segment* curr_segment; 
    segment* segments;
    Cell* cell;

    /* Set each tally to zero to begin with */
    #pragma omp parallel for
    for (int r=0; r < _num_FSRs; r++)
        FSR_segment_tallies[r] = 0;

    /* Iterate over all azimuthal angles, all tracks, and all segments
     * and tally each segment in the corresponding FSR */
    #pragma omp parallel for private (num_segments, curr_segment)
    for (int i=0; i < _tot_num_tracks; i++) {
     
        num_segments = _tracks[i]->getNumSegments();
	segments = _tracks[i]->getSegments();

	for (int s=0; s < num_segments; s++) {
	    curr_segment = &segments[s];
	    FSR_segment_tallies[curr_segment->_region_id]++;
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
 * @brief Computes keff by performing a series of transport sweep and 
 *        source updates.
 * @details This is the main method exposed to the user through the Python
 *          interface to run a simulation. The method makes an initial guess
 *          for the scalar and boundary fluxes and peforms transport sweeps
 *          and source updates until convergence.
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

    /* Clear all timing data from a previous simulation run */
    clearTimerSplits();

    /* Start the timer to record the total time to converge the source */
    _timer->startTimer();

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
	transportSweep();	
	addSourceToScalarFlux();
	computeKeff();

	_num_iterations++;

	/* Check for convergence of the fission source distribution */
	if (i > 1 && residual < _source_convergence_thresh) {
	    _timer->stopTimer();
	    _timer->recordSplit("Total time to converge the source");
	    return _k_eff;
	}
    }

    _timer->stopTimer();
    _timer->recordSplit("Total time to converge the source");

    log_printf(WARNING, "Unable to converge the source after %d iterations",
	       max_iterations);

    return _k_eff;
}


/**
 * @brief Deletes the Timer's timing entries for each timed code section
 *        code in the source convergence loop
 */
void Solver::clearTimerSplits() {
    _timer->clearSplit("Total time to converge the source");
}


/**
 * @brief Prints a report of the timing statistics to the console
 */
void Solver::printTimerReport() {

    std::string msg_string;
    
    log_printf(TITLE, "TIMING REPORT");

    /* Get the total runtime */
    double tot_time = _timer->getSplit("Total time to converge the source");
    msg_string = "Total time to solution";
    msg_string.resize(53, '.');
    log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), tot_time);

    /* Time per unknown */
    double time_per_unknown = tot_time / (_num_FSRs * _num_groups);
    msg_string = "Solution time per unknown";
    msg_string.resize(53, '.');
    log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), time_per_unknown);

    /* Time per iteration */
    double time_per_iter = tot_time / _num_iterations;
    msg_string = "Solution time per iteration";
    msg_string.resize(53, '.');
    log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), time_per_iter);

    /* Time per segment */
    int num_segments = _track_generator->getNumSegments();
    double time_per_segment = (time_per_iter / num_segments);
    msg_string = "Integration time per track segment";
    msg_string.resize(53, '.');
    log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), time_per_segment);

    log_printf(SEPARATOR, "*");
}
