#include "Solver.h"
/**
 * @brief Constructor initializes an empty Solver class with array pointers
 *        set to NULL.
 * @details If the constructor receives Geometry and TrackGenerator 
 *          objects it will retrieve the number of energy groups
 *          and FSRs from the Geometry and azimuthal angles from the
 *          TrackGenerator.
 * @param geometry an optional pointer to a Geometry object
 * @param track_generator an optional pointer to a TrackGenerator object
 */
Solver::Solver(Geometry* geometry, TrackGenerator* track_generator) {

  /* Default values */
  _num_materials = 0;
  _num_groups = 0;
  _num_azim = 0;

  _num_FSRs = 0;
  _num_fissionable_FSRs = 0;
  _num_mesh_cells = 0;
  _FSR_volumes = NULL;
  _FSR_materials = NULL;
  _surface_currents = NULL;

  _track_generator = NULL;
  _geometry = NULL;
  _cmfd = NULL;
  _exp_evaluator = new ExpEvaluator();
  _solve_3D = false;
  
  _tracks = NULL;
  _azim_spacings = NULL;
  _polar_spacings = NULL;
  _boundary_flux = NULL;
  _boundary_leakage = NULL;

  _scalar_flux = NULL;
  _fission_sources = NULL;
  _scatter_sources = NULL;
  _old_fission_sources = NULL;
  _reduced_sources = NULL;
  _source_residuals = NULL;

  /* Default polar quadrature */
  _fluxes_per_track = 0;
  
  if (track_generator != NULL)
    setTrackGenerator(track_generator);

  if (geometry != NULL){
    _cmfd = geometry->getCmfd();
    setGeometry(geometry);
  }
  
  _num_iterations = 0;
  _source_convergence_thresh = 1E-3;
  _converged_source = false;

  _timer = new Timer();

}


/**
 * @brief Destructor deletes arrays of boundary angular fluxes,
 *        scalar fluxes and sources for each FSR and energy group.
 * @details Deallocates memory for all arrays allocated for the Solver,
 *          including fluxes, sources, quadrature weights, and exponential
 *          linear interpolation table.
 */
Solver::~Solver() {

  if (_FSR_volumes != NULL)
    delete [] _FSR_volumes;

  if (_FSR_materials != NULL)
    delete [] _FSR_materials;

  if (_boundary_flux != NULL)
    delete [] _boundary_flux;

  if (_scalar_flux != NULL)
    delete [] _scalar_flux;

  if (_fission_sources != NULL)
    delete [] _fission_sources;

  if (_scatter_sources != NULL)
    delete [] _scatter_sources;

  if (_old_fission_sources != NULL)
    delete [] _old_fission_sources;

  if (_reduced_sources != NULL)
    delete [] _reduced_sources;

  if (_source_residuals != NULL)
    delete [] _source_residuals;

  if (_exp_evaluator != NULL)
    delete _exp_evaluator;

  if (_quad != NULL)
    delete _quad;
}


/**
 * @brief Returns a pointer to the Geometry.
 * @return a pointer to the Geometry
 */
Geometry* Solver::getGeometry() {

  if (_geometry == NULL)
    log_printf(ERROR, "Unable to return the Solver's Geometry since it "
               "has not yet been set");

  return _geometry;
}


/**
 * @brief Returns the calculated volume for a flat source region.
 * @param fsr_id the flat source region ID of interest
 * @return the flat source region volume
 */
FP_PRECISION Solver::getFSRVolume(int fsr_id) {

  if (fsr_id < 0 || fsr_id > _num_FSRs)
    log_printf(ERROR, "Unable to get the volume for FSR %d since the FSR "
               "IDs lie in the range (0, %d)", fsr_id, _num_FSRs);

  else if (_FSR_volumes == NULL)
    log_printf(ERROR, "Unable to get the volume for FSR %d since the FSR "
               "volumes have not yet been computed", fsr_id);

  return _FSR_volumes[fsr_id];
}


/**
 * @brief Returns a pointer to the TrackGenerator.
 * @return a pointer to the TrackGenerator
 */
TrackGenerator* Solver::getTrackGenerator() {

  if (_track_generator == NULL)
    log_printf(ERROR, "Unable to return the Solver's TrackGenetrator "
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
 * @brief Returns the number of source iterations to converge the source.
 * @return the number of iterations
 */
int Solver::getNumIterations() {
  return _num_iterations;
}


/**
 * @brief Returns the total time to converge the source (seconds).
 * @return the time to converge the source (seconds)
 */
double Solver::getTotalTime() {
  return _timer->getSplit("Total time to converge the source");
}


/**
 * @brief Returns the converged eigenvalue \f$ k_{eff} \f$.
 * @return the converged eigenvalue \f$ k_{eff} \f$
 */
FP_PRECISION Solver::getKeff() {
  return _k_eff;
}


/**
 * @brief Returns the threshold for source convergence.
 * @return the threshold for source convergence
 */
FP_PRECISION Solver::getSourceConvergenceThreshold() {
  return _source_convergence_thresh;
}


/**
 * @brief Returns whether the solver is using double floating point precision.
 * @return true if so, false otherwise
 */
bool Solver::isUsingDoublePrecision() {
#ifdef DOUBLE
  return true;
#else
  return false;
#endif
}


/**
 * @brief Returns whether the Solver uses linear interpolation to
 *        compute exponentials.
 * @return true if so, false otherwise
 */
bool Solver::isUsingExponentialInterpolation() {
  return _exp_evaluator->isUsingInterpolation();
}


/**
 * @brief Sets the Geometry for the Solver.
 * @details The Geometry must already have initialized FSR offset maps
 *          and segmentized the TrackGenerator's tracks. Each of these
 *          should be initiated in Python prior to assigning a Geometry
 *          to the Solver:
 *
 * @code
 *          geometry.initializeFlatSourceRegions()
 *          track_generator.generateTracks()
 * @endcode
 *
 * @param geometry a pointer to a Geometry object
 */
void Solver::setGeometry(Geometry* geometry) {

  if (geometry->getNumFSRs() == 0)
    log_printf(ERROR, "Unable to set the Geometry for the Solver since the "
               "Geometry has not yet initialized FSRs");

  _geometry = geometry;
  _num_FSRs = _geometry->getNumFSRs();
  _num_groups = _geometry->getNumEnergyGroups();
  _num_materials = _geometry->getNumMaterials();

  if (_solve_3D)
    _fluxes_per_track = _num_groups;
  else
    _fluxes_per_track = _num_groups * _num_polar/2;
  
  if (_cmfd != NULL)
    _num_mesh_cells = _cmfd->getNumCells();
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
void Solver::setTrackGenerator(TrackGenerator* track_generator) {

  if ((!track_generator->contains2DTracks() && track_generator->isSolve2D()) ||
      (!track_generator->contains3DTracks() && track_generator->isSolve3D()))
    log_printf(ERROR, "Unable to set the TrackGenerator for the Solver "
               "since the TrackGenerator has not yet generated tracks");

  _track_generator = track_generator;
  _solve_3D = _track_generator->isSolve3D();  
  _num_azim = _track_generator->getNumAzim();
  _tracks_per_cycle = _track_generator->getTracksPerCycle();
  _cycles_per_azim = _track_generator->getCyclesPerAzim();
  _tracks_per_plane = _track_generator->getTracksPerPlane();
  _azim_spacings = _track_generator->getAzimSpacings();
  _polar_spacings = _track_generator->getPolarSpacings();
  _quad = _track_generator->getQuadrature();
  _num_polar = _quad->getNumPolarAngles();
  
  /* Initialize the tracks array */
  int counter = 0;
  _num_tracks = new int*[_num_azim/4];
  
  if (!_solve_3D){
    _fluxes_per_track = _num_groups * _num_polar/2;
    _tot_num_tracks = _track_generator->getNum2DTracks();
    _tracks = new Track*[_tot_num_tracks];

    for (int a=0; a < _num_azim/4; a++){
      _num_tracks[a] = new int[_cycles_per_azim[a]];
      for (int c=0; c < _cycles_per_azim[a]; c++){
        _num_tracks[a][c] = counter;
        for (int t=0; t < _tracks_per_cycle[a]; t++){
          _tracks[counter] = &_track_generator->get2DTracks()[a][c][t];
          counter++;
        }
      }
    }

  }
  else{
    _fluxes_per_track = _num_groups;
    _tot_num_tracks = _track_generator->getNum3DTracks();
    _tracks = new Track*[_tot_num_tracks];
    
    for (int a=0; a < _num_azim/4; a++){
      _num_tracks[a] = new int[_cycles_per_azim[a]];
      for (int c=0; c < _cycles_per_azim[a]; c++){
        _num_tracks[a][c] = counter;
        for (int p=0; p < _num_polar; p++){
          for (int i=0; i < _track_generator->getNumZ(a,p)
                 + _track_generator->getNumL(a,p); i++){
            for (int t=0; t < _tracks_per_plane[a][c][p][i]; t++){
              _tracks[counter] = &_track_generator->get3DTracks()[a][c][p][i][t];
              counter++;
            }
          }
        }
      }
    }
  }
}


/**
 * @brief Sets the threshold for source convergence (>0)
 * @param source_thresh the threshold for source convergence
 */
void Solver::setSourceConvergenceThreshold(FP_PRECISION source_thresh) {

  if (source_thresh <= 0.0)
    log_printf(ERROR, "Unable to set the source convergence threshold to %f"
               "since the threshold must be a positive number", source_thresh);

  _source_convergence_thresh = source_thresh;
}


/**
 * @brief Informs the Solver to use linear interpolation to compute the
 *        exponential in the transport equation.
 */
void Solver::useExponentialInterpolation() {
  _exp_evaluator->useInterpolation();
}


/**
 * @brief Informs the Solver to use the exponential intrinsic exp(...) function
 *        to compute the exponential in the transport equation
 */
void Solver::useExponentialIntrinsic() {
  _exp_evaluator->useIntrinsic();
}


/**
 * @brief Initializes new ExpEvaluator object to compute exponentials.
 */
void Solver::initializeExpEvaluator() {
  double max_tau = _track_generator->getMaxOpticalLength();
  double tolerance = _source_convergence_thresh;

  _exp_evaluator->setSolve3D(_solve_3D);
  _exp_evaluator->setQuadrature(_quad);
  _exp_evaluator->initialize(max_tau, tolerance);
}


/**
 * @brief Initializes the FSR volumes and Materials array.
 * @details This method assigns each FSR a unique, monotonically increasing
 *          ID, sets the Material for each FSR, and assigns a volume based on
 *          the cumulative length of all of the segments inside the FSR.
 */
void Solver::initializeFSRs() {

  log_printf(INFO, "Initializing flat source regions...");

  /* Delete old FSR arrays if they exist */
  if (_FSR_volumes != NULL)
    delete [] _FSR_volumes;

  if (_FSR_materials != NULL)
    delete [] _FSR_materials;

  /* Get an array of volumes indexed by FSR  */
  if (_solve_3D)
    _FSR_volumes = _track_generator->get3DFSRVolumes();
  else
    _FSR_volumes = _track_generator->get2DFSRVolumes();
  
  /* Allocate an array of Material pointers indexed by FSR */
  _FSR_materials = new Material*[_num_FSRs];

  /* Compute the number of fissionable Materials */
  _num_fissionable_FSRs = 0;

  /* Loop over all FSRs to extract FSR material pointers */
  for (int r=0; r < _num_FSRs; r++) {

    /* Assign the Material corresponding to this FSR */
    _FSR_materials[r] = _geometry->findFSRMaterial(r);

    /* Increment number of fissionable FSRs */
    if (_FSR_materials[r]->isFissionable())
      _num_fissionable_FSRs++;

    log_printf(DEBUG, "FSR ID = %d has Material ID = %d and volume = %f ",
               r, _FSR_materials[r]->getId(), _FSR_volumes[r]);
  }

  return;
}


/**
 * @brief Initializes a Cmfd object for acceleratiion prior to source iteration.
 * @details Instantiates a dummy Cmfd object if one was not assigned to
 *          the Solver by the user and initializes FSRs, materials, fluxes
 *          and the Mesh object. This method is for internal use only and is
 *          called by the Solver::convergeSource() method and should not be
 *          called directly by the user.
 */
void Solver::initializeCmfd(){

  log_printf(INFO, "Initializing CMFD...");

  /* Give CMFD number of FSRs and FSR property arrays */
  _cmfd->setNumFSRs(_num_FSRs);
  _cmfd->setFSRVolumes(_FSR_volumes);
  _cmfd->setFSRMaterials(_FSR_materials);
  _cmfd->setFSRFluxes(_scalar_flux);
  _cmfd->setQuadrature(_quad);
}



/**
 * @brief Checks that each FSR has at least one Track segment crossing it
 *        and if not, throws an exception and prints an error message.
 * @details This method is for internal use only and is called by the
 *          Solver::convergeSource() method and should not be called
 *          directly by the user.
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

  /* Iterate over all azimuthal angles, all tracks, and all Track segments
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
      log_printf(ERROR, "No tracks were tallied inside FSR id = %d. Please "
                 "reduce your track spacing, increase the number of azimuthal"
                 "angles, or increase the size of the FSRs", r);
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
 *          and source updates until convergence. The method may be called
 *          by the user from Python as follows:
 *
 * @code
 *          max_iters = 1000
 *          solver.convergeSource(max_iters)
 * @endcode
 *
 * @param max_iterations the maximum number of source iterations to allow
 * @return the value of the computed eigenvalue \f$ k_{eff} \f$
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

  /* The old residual and k_eff */
  FP_PRECISION residual_old = 1.0;
  FP_PRECISION keff_old = 1.0;

  /* Initialize data structures */
  initializeExpEvaluator();
  initializeFluxArrays();
  initializeSourceArrays();
  initializeFSRs();

  if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
    initializeCmfd();

  /* Check that each FSR has at least one segment crossing it */
  checkTrackSpacing();
  
  /* Set scalar flux to unity for each region */
  flattenFSRSources(1.0);
  flattenFSRFluxes(1.0);
  zeroTrackFluxes();

  /* Source iteration loop */
  for (int i=0; i < max_iterations; i++) {

    log_printf(NORMAL, "Iteration %d: \tk_eff = %1.6f"
               "\tres = %1.3E", i, _k_eff, residual);

    normalizeFluxes();
    residual = computeFSRSources();
    transportSweep();
    addSourceToScalarFlux();

    /* Solve CMFD diffusion problem and update MOC flux */
    if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
      _k_eff = _cmfd->computeKeff(i);
    else
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
 *        code in the source convergence loop.
 */
void Solver::clearTimerSplits() {
  _timer->clearSplit("Total time to converge the source");
}


/**
 * @brief Prints a report of the timing statistics to the console.
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
  int num_segments = 0;
  if (_solve_3D)
    num_segments = _track_generator->getNum3DSegments();
  else
    num_segments = _track_generator->getNum2DSegments();

  int num_integrations = _fluxes_per_track * num_segments;
  double time_per_integration = (time_per_iter / num_integrations);
  msg_string = "Integration time per segment integration";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), time_per_integration);

  set_separator_character('-');
  log_printf(SEPARATOR, "-");

  msg_string = "           # tracks          # segments          # FSRs";
  log_printf(RESULT, "%s", msg_string.c_str());
  log_printf(SEPARATOR, "-");

  int num_digits = (int) log10((double) _tot_num_tracks);
  num_digits += (int) log10((double) num_segments);
  num_digits += (int) log10((double) _num_FSRs);

  num_digits = 66 - num_digits;
  num_digits /= 4;

  std::stringstream msg;

  for (int i=0; i < 4; i++) {
    for (int j=0; j < num_digits; j++)
      msg << " ";

    if (i == 0)
      msg << _tot_num_tracks;
    else if (i == 1)
      msg << num_segments;
    else if (i == 2)
      msg << _num_FSRs;
  }

  log_printf(RESULT, "%s", msg.str().c_str());
  log_printf(SEPARATOR, "-");
}
