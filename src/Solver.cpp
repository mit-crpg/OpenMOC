#include "Solver.h"
/**
 * @brief Constructor initializes an empty Solver class with array pointers
 *        set to NULL.
 * @param track_generator an optional pointer to a TrackGenerator object
 */
Solver::Solver(TrackGenerator* track_generator) {

  /* Default values */
  _num_materials = 0;
  _num_groups = 0;
  _num_azim = 0;

  _num_FSRs = 0;
  _num_fissionable_FSRs = 0;
  _FSR_volumes = NULL;
  _FSR_materials = NULL;

  _track_generator = NULL;
  _geometry = NULL;
  _cmfd = NULL;
  _exp_evaluator = new ExpEvaluator();
  _max_optical_length = 10;

  _tracks = NULL;
  _polar_weights = NULL;
  _boundary_flux = NULL;
  _boundary_leakage = NULL;

  _scalar_flux = NULL;
  _fission_sources = NULL;
  _scatter_sources = NULL;
  _fixed_sources = NULL;
  _old_fission_sources = NULL;
  _reduced_sources = NULL;
  _source_residuals = NULL;

  if (track_generator != NULL)
    setTrackGenerator(track_generator);

  /* Default polar quadrature */
  _user_polar_quad = false;
  _polar_quad = new TYPolarQuad();
  _num_polar = 3;
  _two_times_num_polar = 2 * _num_polar;
  _polar_times_groups = 0;

  _num_iterations = 0;
  _converge_thresh = 1E-5;
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

  if (_fixed_sources != NULL)
    delete [] _fixed_sources;

  if (_old_fission_sources != NULL)
    delete [] _old_fission_sources;

  if (_reduced_sources != NULL)
    delete [] _reduced_sources;

  if (_source_residuals != NULL)
    delete [] _source_residuals;

  if (_exp_evaluator != NULL)
    delete _exp_evaluator;

  if (_polar_quad != NULL && !_user_polar_quad)
    delete _polar_quad;
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
  return _converge_thresh;
}


/**
 * @brief Get the maximum allowable optical length for a track segment
 * @return The max optical length
 */
FP_PRECISION Solver::getMaxOpticalLength() {
  return _max_optical_length;
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
 * @details This is a private setter method for the Solver and is not
 *          intended to be called by the user.
 * @param geometry a pointer to a Geometry object
 */
void Solver::setGeometry(Geometry* geometry) {

  if (geometry->getNumFSRs() == 0)
    log_printf(ERROR, "Unable to set the Geometry for the Solver since the "
               "Geometry has not yet initialized FSRs");

  _geometry = geometry;
  _cmfd = geometry->getCmfd();
  _num_FSRs = _geometry->getNumFSRs();
  _num_groups = _geometry->getNumEnergyGroups();
  _polar_times_groups = _num_groups * _num_polar;
  _num_materials = _geometry->getNumMaterials();
}


/**
 * @brief Set the maximum allowable optical length for a track segment
 * @param max_optical_length The max optical length
 */
void Solver::setMaxOpticalLength(FP_PRECISION max_optical_length) {

  if (max_optical_length <= 0)
    log_printf(ERROR, "Cannot set max optical length to %f because it "
               "must be positive.", max_optical_length); 
        
  _max_optical_length = max_optical_length;
}


/**
 * @brief Sets the Solver's TrackGenerator with characteristic Tracks.
 * @details The TrackGenerator must already have generated Tracks and have
 *          used ray tracing to segmentize them across the Geometry. This
 *          should be initated in Python prior to assigning the TrackGenerator
 *          to the Solver:
 *
 * @code
 *          geometry.initializeFlatSourceRegions()
 *          track_generator.generateTracks()
 *          solver.setTrackGenerator(track_generator)
 * @endcode
 *
 * @param track_generator a pointer to a TrackGenerator object
 */
void Solver::setTrackGenerator(TrackGenerator* track_generator) {

  if (!track_generator->containsTracks())
    log_printf(ERROR, "Unable to set the TrackGenerator for the Solver "
               "since the TrackGenerator has not yet generated tracks");

  _track_generator = track_generator;
  _num_azim = _track_generator->getNumAzim() / 2;
  int* num_tracks = _track_generator->getNumTracksArray();
  _tot_num_tracks = _track_generator->getNumTracks();
  _tracks = new Track*[_tot_num_tracks];

  /* Initialize the tracks array */
  int counter = 0;

  for (int i=0; i < _num_azim; i++) {
    for (int j=0; j < num_tracks[i]; j++) {
      _tracks[counter] = &_track_generator->getTracks()[i][j];
      counter++;
    }
  }

  /* Retrieve and store the Geometry from the TrackGenerator */  
  setGeometry(_track_generator->getGeometry());
}


/**
 * @brief Assign a PolarQuad object to the Solver.
 * @details This routine allows use of a PolarQuad with any polar angle
 *          quadrature. Alternatively, this routine may take in any subclass
 *          of the PolarQuad parent class, including TYPolarQuad (default),
 *          LeonardPolarQuad, GLPolarQuad, etc.
 *
 *          Users may assign a PolarQuad object to the Solver from 
 *          Python script as follows:
 *
 * @code
 *          polar_quad = openmoc.LeonardPolarQuad()
 *          polar_quad.setNumPolarAngles(2)
 *          solver.setPolarQuadrature(polar_quad)
 * @endcode
 *
 * @param polar_quad a pointer to a PolarQuad object
 */
void Solver::setPolarQuadrature(PolarQuad* polar_quad) {
  _user_polar_quad = true;
  _polar_quad = polar_quad;
  _num_polar = _polar_quad->getNumPolarAngles();
  _two_times_num_polar = 2 * _num_polar;
  _polar_times_groups = _num_groups * _num_polar;
}


/**
 * @brief Sets the threshold for source/flux convergence (>0)
 * @brief The default threshold for convergence is 1E-5.
 * @param source_thresh the threshold for source/flux convergence
 */
void Solver::setSourceConvergenceThreshold(FP_PRECISION source_thresh) {

  if (source_thresh <= 0.0)
    log_printf(ERROR, "Unable to set the convergence threshold to %f "
               "since it is not a positive number", source_thresh);

  _converge_thresh = source_thresh;
}


/**
 * @brief Assign a fixed source for a flat source region and energy group.
 * @details This is a helper routine to perform error checking for the
 *          subclasses which store the source in the appropriate array.
 * @param fsr_id the flat source region ID
 * @param group the energy group
 * @param source the volume-averaged source in this group
 */
void Solver::setFixedSourceByFSR(int fsr_id, int group, FP_PRECISION source) {
  
  if (group <= 0 || group > _num_groups)
    log_printf(ERROR,"Unable to set fixed source for group %d in "
               "in a %d group problem", group, _num_groups);

  if (fsr_id < 0 || fsr_id >= _num_FSRs)
    log_printf(ERROR,"Unable to set fixed source for FSR %d with only "
               "%d FSRs in the geometry", fsr_id, _num_FSRs);
}


/**
 * @brief Assign a fixed source for a Cell and energy group.
 * @details This routine will add the fixed source to all instances of the
 *          Cell in the geometry (e.g., all FSRs for this Cell).
 * @param fsr_id the Cell of interest
 * @param group the energy group
 * @param source the volume-averaged source in this group
 */
void Solver::setFixedSourceByCell(Cell* cell, int group, FP_PRECISION source) {

  /* If the Cell is filled by a Universe, recursively
   * add the source to all Cells within it */
  if (cell->getType() == FILL) {
    std::map<int, Cell*> cells = cell->getAllCells();
    std::map<int, Cell*>::iterator iter;
    for (iter = cells.begin(); iter != cells.end(); ++iter)
      setFixedSourceByCell(iter->second, group, source);
  }

  /* If the Cell is filled by a Material, add the 
   * source to all FSRs for this Cell */
  else {
    CellBasic* fsr_cell;
    
    for (int r=0; r < _num_FSRs; r++) {
      fsr_cell = _geometry->findCellContainingFSR(r);
      if (cell->getId() == fsr_cell->getId())
        setFixedSourceByFSR(r, group, source);
    }
  }
}


/**
 * @brief Assign a fixed source for a Material and energy group.
 * @details This routine will add the fixed source to all instances of the
 *          Material in the geometry (e.g., all FSRs with this Material).
 * @param fsr_id the Material of interest
 * @param group the energy group
 * @param source the volume-averaged source in this group
 */
void Solver::setFixedSourceByMaterial(Material* material, int group, 
                                      FP_PRECISION source) {

  Material* fsr_material;

  /* Add the source to all FSRs for this Material */
  for (int r=0; r < _num_FSRs; r++) {
    fsr_material = _geometry->findFSRMaterial(r);
    if (material->getId() == fsr_material->getId())
      setFixedSourceByFSR(group, r, source);
  }
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
 * @brief Initializes new PolarQuad object.
 * @details Deletes memory old PolarQuad if one was previously allocated.
 */
void Solver::initializePolarQuadrature() {

  FP_PRECISION azim_weight;
  FP_PRECISION* azim_weights = _track_generator->getAzimWeights();

  /* Create Tabuchi-Yamamoto polar quadrature if a
   * PolarQuad was not assigned by the user */
  if (_polar_quad == NULL)
    _polar_quad = new TYPolarQuad();

  /* Initialize the PolarQuad object */
  _polar_quad->setNumPolarAngles(_num_polar);
  _polar_quad->initialize();
  _polar_times_groups = _num_groups * _num_polar;

  /* Deallocate polar weights if previously assigned */
  if (_polar_weights != NULL)
    delete [] _polar_weights;

  _polar_weights = new FP_PRECISION[_num_azim*_num_polar];

  /* Compute the total azimuthal weight for tracks at each polar angle */
  #pragma omp parallel for private(azim_weight) schedule(guided)
  for (int i=0; i < _num_azim; i++) {
    azim_weight = azim_weights[i];

    for (int p=0; p < _num_polar; p++)
      _polar_weights(i,p) = 
           azim_weight * _polar_quad->getMultiple(p) * FOUR_PI;
  }
}


/**
 * @brief Initializes new ExpEvaluator object to compute exponentials.
 */
void Solver::initializeExpEvaluator() {
  _exp_evaluator->setPolarQuadrature(_polar_quad);
  _exp_evaluator->initialize(_max_optical_length, _converge_thresh);

  if (_exp_evaluator->isUsingInterpolation())
    _track_generator->splitSegments(_max_optical_length);
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
  _FSR_volumes = _track_generator->getFSRVolumes();

  /* Allocate an array of Material pointers indexed by FSR */
  _FSR_materials = new Material*[_num_FSRs];

  /* Compute the number of fissionable Materials */
  std::map<int, Material*> all_materials = _geometry->getAllMaterials();
  _num_fissionable_FSRs = 0;

  /* Loop over all FSRs to extract FSR material pointers */
  for (int r=0; r < _num_FSRs; r++) {

    /* Assign the Material corresponding to this FSR */
    _FSR_materials[r] =  _geometry->findFSRMaterial(r);

    /* Increment number of fissionable FSRs */
    if (_FSR_materials[r]->isFissionable())
      _num_fissionable_FSRs++;

    log_printf(DEBUG, "FSR ID = %d has Material ID = %d and volume = %f ",
               r, _FSR_materials[r]->getId(), _FSR_volumes[r]);
  }
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
  _cmfd->setPolarQuadrature(_polar_quad);
  _cmfd->initializeSurfaceCurrents();
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
                 "reduce your track spacing, increase the number of "
                 "azimuthal angles, or increase the size of the FSRs", r);
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
 */
void Solver::convergeSource(int max_iterations) {

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

  /* The new/old residuals on the fission source */
  FP_PRECISION residual = 0.0;

  /* Initialize data structures */
  initializePolarQuadrature();
  initializeExpEvaluator();
  initializeFluxArrays();
  initializeSourceArrays();
  initializeFSRs();

  if (_num_fissionable_FSRs == 0.0)
    log_printf(ERROR, "The Solver is unable to converge the "
               "source without fissionable flat source regions");

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
    if (_cmfd != NULL && _cmfd->isFluxUpdateOn()){
      _k_eff = _cmfd->computeKeff(i);
      _cmfd->updateBoundaryFlux(_tracks, _boundary_flux, _tot_num_tracks);
    }
    else
      computeKeff();

    _num_iterations++;

    /* Check for convergence of the fission source distribution */
    if (i > 1 && residual < _converge_thresh) {
      _timer->stopTimer();
      _timer->recordSplit("Total time to converge the source");
      return;
    }
  }

  log_printf(WARNING, "Unable to converge the source");

  _timer->stopTimer();
  _timer->recordSplit("Total time to converge the source");
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
  int num_segments = _track_generator->getNumSegments();
  int num_integrations = 2 * _num_polar * _num_groups * num_segments;
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


/**
 * @brief Solves fixed source problem
 * @code
 *          max_iters = 1000
 *          solver.convergeFlux(max_iters)
 * @endcode
 *
 * @param max_iterations the maximum number of iterations to allow
 */
void Solver::convergeFlux(int max_iterations) {

  if (_track_generator == NULL)
    log_printf(ERROR, "The Solver is unable to converge the flux "
               "since it does not contain a TrackGenerator");

  log_printf(NORMAL, "Converging the flux...");

  /* Clear all timing data from a previous simulation run */
  clearTimerSplits();

  /* Start the timer to record the total time to converge the flux */
  _timer->startTimer();

  /* Counter for the number of iterations to converge the flux */
  _num_iterations = 0;

  //FIXME
  /* An initial guess for the eigenvalue */
  //  _k_eff = 1.0;

  /* The residual on the source */
  FP_PRECISION residual = 0.0;

  /* Initialize data structures */
  initializePolarQuadrature();
  initializeExpEvaluator();
  initializeFluxArrays();
  initializeSourceArrays();
  initializeFSRs();

  //FIXME
  /*
  if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
    initializeCmfd();
  */

  /* Check that each FSR has at least one segment crossing it */
  checkTrackSpacing();

  /* Set scalar flux to unity for each region */
  flattenFSRSources(1.0);
  flattenFSRFluxes(1.0);
  zeroTrackFluxes();

  /* Source iteration loop */
  for (int i=0; i < max_iterations; i++) {

    log_printf(NORMAL, "Iteration %d: \tres = %1.3E", i, residual);

    //FIXME
    //    residual = computeFSRSourcesForFixedSource();
    residual = computeFSRSources();
    transportSweep();
    addSourceToScalarFlux();

    //FIXME
    /* Solve CMFD diffusion problem and update MOC flux */
    /*
    if (_cmfd != NULL && _cmfd->isFluxUpdateOn()){
      _k_eff = _cmfd->computeKeff(i);
      _cmfd->updateBoundaryFlux(_tracks, _boundary_flux, _tot_num_tracks);
    }
    else
      computeKeff();
    */

    _num_iterations++;

    /* Check for convergence */
    if (i > 1 && residual < _converge_thresh) {
      _timer->stopTimer();
      _timer->recordSplit("Total time to converge the flux");
      return;
    }
  }

  _timer->stopTimer();
  _timer->recordSplit("Total time to converge the flux");

  log_printf(WARNING, "Unable to converge the flux");
}
