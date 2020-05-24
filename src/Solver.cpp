#include "Solver.h"
#include <unordered_map>
#include <fstream>
#include <sys/stat.h>
#ifdef BGQ
#include <spi/include/kernel/memory.h>
#endif

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
  _chi_spectrum_material = NULL;

  _k_eff = 1.;
  _keff_from_fission_rates = true;

  _track_generator = NULL;
  _geometry = NULL;
  _cmfd = NULL;
  _num_exp_evaluators_azim = 1;
  _num_exp_evaluators_polar = 1;
  _exp_evaluators = new ExpEvaluator**[_num_exp_evaluators_azim];
  _exp_evaluators[0] = new ExpEvaluator*[_num_exp_evaluators_polar];
  _exp_evaluators[0][0] = new ExpEvaluator();
#ifndef THREED
  _SOLVE_3D = false;
#endif
  _segment_formation = EXPLICIT_2D;

  /* Initialize pointers to NULL */
  _tracks = NULL;
  _boundary_flux = NULL;
  _start_flux = NULL;
  _boundary_leakage = NULL;

  _scalar_flux = NULL;
  _old_scalar_flux = NULL;
  _reference_flux = NULL;
  _stabilizing_flux = NULL;
  _fixed_sources = NULL;
  _reduced_sources = NULL;
  _source_type = "None";

  _regionwise_scratch = NULL;

  _fluxes_per_track = 0;

  if (track_generator != NULL)
    setTrackGenerator(track_generator);

  _num_iterations = 0;
  _converge_thresh = 1E-5;

  _timer = new Timer();

  /* Default settings */
  _solver_mode = FORWARD;
  _is_restart = false;
  _user_fluxes = false;
  _fixed_sources_on = false;
  _fixed_sources_initialized = false;
  _correct_xs = false;
  _stabilize_transport = false;
  _verbose = false;
  _calculate_initial_spectrum = false;
  _initial_spectrum_thresh = 1.0;
  _load_initial_FSR_fluxes = false;
  _calculate_residuals_by_reference = false;
  _negative_fluxes_allowed = false;
  _print_negative_sources = false;
  _OTF_transport = false;
  _xs_log_level = ERROR;
  _gpu_solver = false;

  //FIXME Parameters for xs modification, should be deleted
  _reset_iteration = -1;
  _limit_xs = false;
}

/**
 * @brief Destructor deletes arrays of boundary angular fluxes,
 *        scalar fluxes and sources for each FSR and energy group.
 * @details Deallocates memory for all arrays allocated for the Solver,
 *          including fluxes, sources, quadrature weights, and exponential
 *          linear interpolation table.
 */
Solver::~Solver() {

  if (_FSR_materials != NULL)
    delete [] _FSR_materials;

  if (_boundary_flux != NULL)
    delete [] _boundary_flux;

  if (_start_flux != NULL)
    delete [] _start_flux;

  if (_scalar_flux != NULL && !_user_fluxes)
    delete [] _scalar_flux;

  if (_old_scalar_flux != NULL)
    delete [] _old_scalar_flux;

  if (_reference_flux != NULL)
    delete [] _reference_flux;

  if (_stabilizing_flux != NULL)
    delete [] _stabilizing_flux;

  if (_fixed_sources != NULL)
    delete [] _fixed_sources;

  if (_reduced_sources != NULL)
    delete [] _reduced_sources;

  if (_boundary_leakage != NULL)
    delete [] _boundary_leakage;

  if (_regionwise_scratch != NULL)
    delete [] _regionwise_scratch;

  for (int i=0; i < _groupwise_scratch.size(); i++)
    delete [] _groupwise_scratch.at(i);
  _groupwise_scratch.clear();

  /** Delete exponential evaluators */
  if (_exp_evaluators != NULL) {
    for (int a=0; a < _num_exp_evaluators_azim; a++) {
      for (int p=0; p < _num_exp_evaluators_polar; p++) {
        delete _exp_evaluators[a][p];
      }
    }
    for (int a=0; a < _num_azim/2; a++) {
      /* Handle edge case when exp_evaluators are not initialized */
      if (a < 1 or _num_exp_evaluators_azim > 1)
        delete [] _exp_evaluators[a];
    }
    delete [] _exp_evaluators;
  }

  delete _timer;
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
FP_PRECISION Solver::getFSRVolume(long fsr_id) {

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
  return _timer->getSplit("Total time");
}


/**
 * @brief Returns the converged eigenvalue \f$ k_{eff} \f$.
 * @return the converged eigenvalue \f$ k_{eff} \f$
 */
double Solver::getKeff() {
  return _k_eff;
}


/**
 * @brief Returns the threshold for source/flux convergence.
 * @return the threshold for source/flux convergence
 */
double Solver::getConvergenceThreshold() {
  return _converge_thresh;
}


/**
 * @brief Get the maximum allowable optical length for a track segment
 * @return The max optical length
 */
FP_PRECISION Solver::getMaxOpticalLength() {
  return _exp_evaluators[0][0]->getMaxOpticalLength();
}


/**
 * @brief Returns whether the solver is using double floating point precision.
 * @return true if using double precision float point arithmetic
 */
bool Solver::isUsingDoublePrecision() {
#ifdef SINGLE
  return false;
#else
  return true;
#endif
}


/**
 * @brief Returns whether the Solver uses linear interpolation to
 *        compute exponentials.
 * @return true if using linear interpolation to compute exponentials
 */
bool Solver::isUsingExponentialInterpolation() {
  return _exp_evaluators[0][0]->isUsingInterpolation();
}


/**
 * @brief Returns whether the Solver is tackling a 3D problem.
 * @return true if the solver is set up with a 3D track generator
 */
bool Solver::is3D() {
  return _SOLVE_3D;
}


/**
 * @brief Returns the scalar flux for some FSR and energy group.
 * @param fsr_id the ID for the FSR of interest
 * @param group the energy group of interest
 * @return the FSR scalar flux
 */
double Solver::getFlux(long fsr_id, int group) {

  if (fsr_id >= _num_FSRs)
    log_printf(ERROR, "Unable to return a scalar flux for FSR ID = %d "
               "since the max FSR ID = %d", fsr_id, _num_FSRs-1);

  else if (fsr_id < 0)
    log_printf(ERROR, "Unable to return a scalar flux for FSR ID = %d "
               "since FSRs do not have negative IDs", fsr_id);

  else if (group-1 >= _num_groups)
    log_printf(ERROR, "Unable to return a scalar flux in group %d "
               "since there are only %d groups", group, _num_groups);

  else if (group <= 0)
    log_printf(ERROR, "Unable to return a scalar flux in group %d "
               "since groups must be greater or equal to 1", group);

  else if (_scalar_flux == NULL)
    log_printf(ERROR, "Unable to return a scalar flux "
             "since it has not yet been computed");

  return _scalar_flux(fsr_id,group-1);
}


/**
 * @brief Returns the source for some energy group for a flat source region
 * @details This is a helper routine used by the openmoc.process module.
 * @param fsr_id the ID for the FSR of interest
 * @param group the energy group of interest
 * @return the flat source region source
 */
double Solver::getFSRSource(long fsr_id, int group) {

  if (fsr_id >= _num_FSRs)
    log_printf(ERROR, "Unable to return a source for FSR ID = %d "
               "since the max FSR ID = %d", fsr_id, _num_FSRs-1);

  else if (fsr_id < 0)
    log_printf(ERROR, "Unable to return a source for FSR ID = %d "
               "since FSRs do not have negative IDs", fsr_id);

  else if (group-1 >= _num_groups)
    log_printf(ERROR, "Unable to return a source in group %d "
               "since there are only %d groups", group, _num_groups);

  else if (group <= 0)
    log_printf(ERROR, "Unable to return a source in group %d "
               "since groups must be greater or equal to 1", group);

  else if (_scalar_flux == NULL)
    log_printf(ERROR, "Unable to return a source "
               "since it has not yet been computed");

  Material* material = _FSR_materials[fsr_id];
  FP_PRECISION* nu_sigma_f = material->getNuSigmaF();
  FP_PRECISION* chi = material->getChi();
  double source = 0.;

  /* Compute fission source */
  if (material->isFissionable()) {
    for (int e=0; e < _num_groups; e++)
      source += _scalar_flux(fsr_id,e) * nu_sigma_f[e];
    source /= _k_eff * chi[group-1];
  }

  /* Compute scatter source */
  for (int g=0; g < _num_groups; g++)
    source += material->getSigmaSByGroup(g+1,group)
              * _scalar_flux(fsr_id,g);

  /* Add in fixed source (if specified by user) */
  if (_fixed_sources_on)
    source += _fixed_sources(fsr_id,group-1);

  /* Normalize to solid angle for isotropic approximation */
  source *= ONE_OVER_FOUR_PI;

  return source;
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

  /* Retrieve the segmentation type and get the 3D track generator */
  segmentationType segment_formation = track_generator->getSegmentFormation();
  TrackGenerator3D* track_generator_3D =
    dynamic_cast<TrackGenerator3D*>(track_generator);

  /* Check to make sure proper segments have been generated */
  if (!track_generator->containsSegments())
    log_printf(ERROR, "Unable to set the TrackGenerator for the Solver "
               "since the TrackGenerator has not yet generated tracks");

  _track_generator = track_generator;
  _segment_formation = _track_generator->getSegmentFormation();
  _num_azim = _track_generator->getNumAzim();
  _quad = _track_generator->getQuadrature();
  _num_polar = _quad->getNumPolarAngles();
  _tracks = _track_generator->getTracksArray();

  /* Set the number of tracks and fluxes per track */
  if (track_generator_3D != NULL) {
    _fluxes_per_track = _num_groups;
    _tot_num_tracks = track_generator_3D->getNum3DTracks();
    _tracks_per_stack = track_generator_3D->getTracksPerStack();
#ifndef THREED
    _SOLVE_3D = true;
#endif
  }
  else {
    _fluxes_per_track = _num_groups * _num_polar/2;
    _tot_num_tracks = _track_generator->getNum2DTracks();
#ifdef THREED
    log_printf(ERROR, "OpenMOC has been compiled for 3D cases only, please "
               "recompile without the -DTHREED optimization flag.");
#else
    _SOLVE_3D = false;
#endif
  }

  /* Retrieve and store the Geometry from the TrackGenerator */
  setGeometry(_track_generator->getGeometry());
}


/**
 * @brief Sets the threshold for source/flux convergence.
 * @brief The default threshold for convergence is 1E-5.
 * @param threshold the threshold for source/flux convergence
 */
void Solver::setConvergenceThreshold(double threshold) {

  if (threshold <= 0.0)
    log_printf(ERROR, "Unable to set the convergence threshold to %f "
               "since it is not a positive number", threshold);

  _converge_thresh = threshold;
}


/**
 * @brief Assign a fixed source for a flat source region and energy group.
 * @details This is a helper routine to perform error checking for the
 *          subclasses which store the source in the appropriate array.
 * @param fsr_id the flat source region ID
 * @param group the energy group
 * @param source the volume-averaged source in this group
 */
void Solver::setFixedSourceByFSR(long fsr_id, int group, double source) {

  /* Initialize part of the solver to be able to set FSR fixed sources */
  if (_num_groups == 0)
    _num_groups = _geometry->getNumEnergyGroups();

  if (_num_FSRs == 0)
    _num_FSRs = _geometry->getNumFSRs();

  if (group <= 0 || group > _num_groups)
    log_printf(ERROR,"Unable to set fixed source for group %d in "
               "in a %d energy group problem", group, _num_groups);

  if (fsr_id < 0 || fsr_id >= _num_FSRs)
    log_printf(ERROR,"Unable to set fixed source for FSR %d with only "
               "%d FSRs in the geometry", fsr_id, _num_FSRs);

  _fixed_sources_on = true;
  _fixed_sources_initialized = false;
  _fix_src_FSR_map[std::pair<int, int>(fsr_id, group)] = source;
}


/**
 * @brief Assign a fixed source for a Cell and energy group.
 * @details This routine will add the fixed source to all instances of the
 *          Cell in the geometry (e.g., all FSRs for this Cell).
 * @param cell a pointer to the Cell of interest
 * @param group the energy group
 * @param source the volume-averaged source in this group
 */
void Solver::setFixedSourceByCell(Cell* cell, int group, double source) {

  _fixed_sources_on = true;
  _fixed_sources_initialized = false;

  /* Recursively add the source to all Cells within a FILL type Cell */
  if (cell->getType() == FILL) {
    std::map<int, Cell*> cells = cell->getAllCells();
    std::map<int, Cell*>::iterator iter;
    for (iter = cells.begin(); iter != cells.end(); ++iter)
      setFixedSourceByCell(iter->second, group, source);
  }

  /* Add the source to all FSRs for this MATERIAL type Cell */
  else {
    _fix_src_cell_map[std::pair<Cell*, int>(cell, group)] = source;
  }
}


/**
 * @brief Assign a fixed source for a Material and energy group.
 * @details This routine will add the fixed source to all instances of the
 *          Material in the geometry (e.g., all FSRs with this Material).
 * @param material a pointer to the Material of interest
 * @param group the energy group
 * @param source the volume-averaged source in this group
 */
void Solver::setFixedSourceByMaterial(Material* material, int group,
                                      double source) {
  _fixed_sources_on = true;
  _fixed_sources_initialized = false;
  _fix_src_material_map[std::pair<Material*, int>(material, group)] = source;
}


/**
 * @brief Set the maximum allowable optical length for a track segment
 * @param max_optical_length The max optical length
 */
void Solver::setMaxOpticalLength(FP_PRECISION max_optical_length) {
  for (int a=0; a < _num_exp_evaluators_azim; a++)
    for (int p=0; p < _num_exp_evaluators_polar; p++)
      _exp_evaluators[a][p]->setMaxOpticalLength(max_optical_length);
}


/**
 * @brief Set the precision, or maximum allowable approximation error, of the
 *        the exponential interpolation table.
 * @details By default, the precision is 1E-5 based on the analysis in
 *          Yamamoto's 2003 paper.
 * @param precision the precision of the exponential interpolation table,
 */
void Solver::setExpPrecision(double precision) {
  for (int a=0; a < _num_exp_evaluators_azim; a++)
    for (int p=0; p < _num_exp_evaluators_polar; p++)
      _exp_evaluators[a][p]->setExpPrecision(precision);
}


/**
 * @brief Informs the Solver to use linear interpolation to compute the
 *        exponential in the transport equation.
 */
void Solver::useExponentialInterpolation() {
  for (int a=0; a < _num_exp_evaluators_azim; a++)
    for (int p=0; p < _num_exp_evaluators_polar; p++)
      _exp_evaluators[a][p]->useInterpolation();
}


/**
 * @brief Informs the Solver to use the exponential intrinsic exp(...)
 *        function to compute the exponential in the transport equation.
 */
void Solver::useExponentialIntrinsic() {
  for (int a=0; a < _num_exp_evaluators_azim; a++)
    for (int p=0; p < _num_exp_evaluators_polar; p++)
      _exp_evaluators[a][p]->useIntrinsic();
}


/**
 * @brief Choose between direct and adjoint mode.
 * @param solver_mode openmoc.FORWARD or .ADJOINT
 */
void Solver::setSolverMode(solverMode solver_mode) {

  _solver_mode = solver_mode;
  if (solver_mode == ADJOINT)
    log_printf(NORMAL, "Solver set to perform an adjoint calculation");
}


/**
 * @brief Informs the Solver that this is a 'restart' calculation and therefore
 *        k_eff, track angular and region scalar fluxes should not be reset.
 * @param is_restart whether solver run is a restart run
 */
void Solver::setRestartStatus(bool is_restart) {

  if (is_restart)
    log_printf(NORMAL, "Solver is in restart mode, no fluxes will be reset.");
  _is_restart = is_restart;
}


/**
 * @brief Informs the Solver that this calculation may involve negative fluxes
 *        for computing higher eigenmodes for example.
 * @param negative_fluxes_on whether to allow negative fluxes
 */
void Solver::allowNegativeFluxes(bool negative_fluxes_on) {
  _negative_fluxes_allowed = negative_fluxes_on;
}


/**
 * @brief Set whether to print negative sources at each iteration.
 * @param print_negative_sources whether to print negative sources
 */
void Solver::printAllNegativeSources(bool print_negative_sources) {
  _print_negative_sources = print_negative_sources;
}


/**
 * @brief Directs OpenMOC to correct unphysical cross-sections.
 * @details If a material is found with greater total scattering cross-section
 *          than total cross-section, the total cross-section is set to the
 *          scattering cross-section.
 */
void Solver::correctXS() {
  _correct_xs = true;
}


/**
 * @brief Directs OpenMOC to use the diagonal stabilizing correction to
 *        the source iteration transport sweep.
 * @details The source iteration process which MOC uses can be unstable
 *          if negative cross-sections arise from transport correction. This
 *          instability causes issues in convergence. The stabilizing
 *          correction fixes this by adding a diagonal matrix to both sides
 *          of the discretized transport equation which introduces no bias
 *          but transforms the iteration matrix into one that is stable.
 *
 *          Three stabilization options exist: DIAGONAL, YAMAMOTO, and GLOBAL.
 *
 *          DIAGONAL: The stabilization is only applied to fluxes where the
 *                    associated in-scatter cross-section is negative. The
 *                    added stabilizing flux is equal to the magnitude of the
 *                    in-scatter cross-section divided by the total
 *                    cross-section and scaled by the stabilization factor.
 *          YAMAMOTO: This is the same as DIAGONAL except that the largest
 *                    stabilization is applied to all regions, not just those
 *                    containing negative in-scatter cross-sections.
 *          GLOBAL: This method applies a global stabilization factor to all
 *                  fluxes defined by the user. In addition, the stabilization
 *                  factor in this option refers to a damping factor, not
 *                  the magnitude of the stabilizing correction.
 *
 * @param stabilization_factor The factor applied to the stabilizing correction
 * @param stabilization_type The type of stabilization to use
 */
void Solver::stabilizeTransport(double stabilization_factor,
                                stabilizationType stabilization_type) {
  _stabilize_transport = true;
  _stabilization_factor = stabilization_factor;
  _stabilization_type = stabilization_type;
}


/**
 * @brief Instructs OpenMOC to perform an initial spectrum calculation
 * @param threshold The convergence threshold of the spectrum calculation
 */
void Solver::setInitialSpectrumCalculation(double threshold) {
  _calculate_initial_spectrum = true;
  _initial_spectrum_thresh = threshold;
}


/**
 * @brief Determines which log level to set cross-section warnings
 * @details The default log level is ERROR
 * @param log_level The log level for outputing cross-section inconsistencies
 */
void Solver::setCheckXSLogLevel(logLevel log_level) {
  _xs_log_level = log_level;
}


/**
 * @brief Sets the chi spectrum for use as an inital flux guess
 * @param material The material used for obtaining the chi spectrum
 */
void Solver::setChiSpectrumMaterial(Material* material) {
  _chi_spectrum_material = material;
}


/**
 * @brief Initializes new ExpEvaluator object to compute exponentials.
 */
void Solver::initializeExpEvaluators() {

  /* Compute the first exponential evaluator */
  ExpEvaluator* first_evaluator = _exp_evaluators[0][0];
  first_evaluator->setQuadrature(_quad);

  /* Find minimum of optional user-specified and actual max taus */
  FP_PRECISION max_tau_a = _track_generator->getMaxOpticalLength();
  FP_PRECISION max_tau_b = first_evaluator->getMaxOpticalLength();

  /* Give track generator a max optical length for segments */
  _track_generator->setMaxOpticalLength(max_tau_b);

  /* Split Track segments so that none has a greater optical length */
  if (max_tau_a > max_tau_b) {

    log_printf(NODAL, "Splitting segments since the maximum optical length in "
               "the domain is %f and the maximum supported by the exponential "
               "evaluator is %f.", max_tau_a, max_tau_b);

    if (_segment_formation == EXPLICIT_3D || _segment_formation == EXPLICIT_2D)
      _track_generator->splitSegments(max_tau_b);
    /* For non explicit ray tracing, segments are split on-the-fly */
    else
      _track_generator->countSegments();
  }

  /* Delete old exponential evaluators */
  for (int a=0; a < _num_exp_evaluators_azim; a++) {
    for (int p=0; p < _num_exp_evaluators_polar; p++)
      if (_exp_evaluators[a][p] != first_evaluator)
        delete _exp_evaluators[a][p];
    delete [] _exp_evaluators[a];
  }
  delete [] _exp_evaluators;

  /* Determine number of exponential evaluators */
  _num_exp_evaluators_azim = _num_azim / 4;
  if (_SOLVE_3D)
    _num_exp_evaluators_polar = _num_polar / 2;
  else
    _num_exp_evaluators_polar = 1;

  /* Allocate new exponential evaluators */
  _exp_evaluators = new ExpEvaluator**[_num_azim/2];
  for (int a=0; a < _num_azim/2; a++)
    _exp_evaluators[a] = new ExpEvaluator*[_num_polar];
  for (int a=0; a < _num_exp_evaluators_azim; a++) {
    for (int p=0; p < _num_exp_evaluators_polar; p++) {

      /* Create a new exponential evaluator if necessary */
      if (a == 0 && p == 0)
        _exp_evaluators[a][p] = first_evaluator;
      else
        _exp_evaluators[a][p] = first_evaluator->deepCopy();

      /* Copy evaluators to supplementary positions */
      int sup_azim = _num_azim / 2 - a - 1;
      int sup_polar = _num_polar - p - 1;
      _exp_evaluators[sup_azim][p] = _exp_evaluators[a][p];
      _exp_evaluators[a][sup_polar] = _exp_evaluators[a][p];
      _exp_evaluators[sup_azim][sup_polar] = _exp_evaluators[a][p];
    }
  }

  /* Initialize exponential interpolation table */
  for (int a=0; a < _num_exp_evaluators_azim; a++)
    for (int p=0; p < _num_exp_evaluators_polar; p++)
      _exp_evaluators[a][p]->initialize(a, p, _SOLVE_3D);
}


/**
 * @brief Initializes the Material's production matrices.
 * @details In an adjoint calculation, this routine will transpose the
 *          scattering and fission matrices in each material.
 * @param mode the solution type (FORWARD or ADJOINT)
 */
void Solver::initializeMaterials(solverMode mode) {

  log_printf(INFO, "Initializing materials...");
  _solver_mode = mode;

  std::map<int, Material*> materials = _geometry->getAllMaterials();
  std::map<int, Material*>::iterator m_iter;

  for (m_iter = materials.begin(); m_iter != materials.end(); ++m_iter) {
    m_iter->second->buildFissionMatrix();

    if (mode == ADJOINT)
      m_iter->second->transposeProductionMatrices();
  }

  /* GPU solver needs this */
  _num_materials = _geometry->getNumMaterials();
}


/**
 * @brief Initializes the FSR volumes and Materials array.
 * @details This method assigns each FSR a unique, monotonically increasing
 *          ID, sets the Material for each FSR, and assigns a volume based on
 *          the cumulative length of all of the segments inside the FSR.
 */
void Solver::initializeFSRs() {

  log_printf(NORMAL, "Initializing solver arrays...");

  /* Delete old FSR arrays if they exist */
  if (_FSR_materials != NULL)
    delete [] _FSR_materials;

  /* Get an array of volumes indexed by FSR  */
  _track_generator->initializeFSRVolumesBuffer();
  _FSR_volumes = _track_generator->getFSRVolumes();

#ifdef NGROUPS
  if (_geometry->getNumEnergyGroups() != NGROUPS)
    log_printf(ERROR, "OpenMOC has been compiled for %d groups, and the "
               "current case is in %d groups, please re-compile with the right "
               "number of groups for the -DNGROUPS flag or without that flag.",
               NGROUPS, _geometry->getNumEnergyGroups());
#endif

  /* Retrieve simulation parameters from the Geometry */
  _num_FSRs = _geometry->getNumFSRs();
  _num_groups = _geometry->getNumEnergyGroups();
  _num_materials = _geometry->getNumMaterials();

  if (_SOLVE_3D) {
    _fluxes_per_track = _num_groups;
  }
  else {
    _fluxes_per_track = _num_groups * _num_polar/2;
  }

  /* Allocate scratch memory */
  for (int i=0; i < _groupwise_scratch.size(); i++)
    delete [] _groupwise_scratch.at(i);
  if (_regionwise_scratch != NULL)
    delete [] _regionwise_scratch;

  int num_threads = omp_get_max_threads();
  _groupwise_scratch.resize(num_threads);
  for (int i=0; i < num_threads; i++)
    _groupwise_scratch.at(i) = new FP_PRECISION[_num_groups];
  _regionwise_scratch = new double[_num_FSRs];

  /* Generate the FSR centroids */
  _track_generator->generateFSRCentroids(_FSR_volumes);

  /* Allocate an array of Material pointers indexed by FSR */
  _FSR_materials = new Material*[_num_FSRs];

  /* Loop over all FSRs to extract FSR material pointers */
  for (long r=0; r < _num_FSRs; r++)
    _FSR_materials[r] = _geometry->findFSRMaterial(r);
}


/**
 * @brief Counts the number of fissionable flat source regions.
 * @details This routine is used by the Solver::computeEigenvalue(...)
 *          routine which uses the number of fissionable FSRs to normalize
 *          the residual on the fission source distribution.
 */
void Solver::countFissionableFSRs() {

  log_printf(INFO_ONCE, "Counting fissionable FSRs...");

  /* Count the number of fissionable FSRs */
  _num_fissionable_FSRs = 0;
  for (long r=0; r < _num_FSRs; r++) {
    if (_FSR_materials[r]->isFissionable())
      _num_fissionable_FSRs++;
  }
}


/**
 * @brief Checks to see if limited XS should be reset.
 * @param iteration The MOC iteration number
 */
void Solver::checkLimitXS(int iteration) {

  if (iteration == _reset_iteration)
    log_printf(NORMAL, "Re-setting material cross-sections");
  else
    return;

  /* Create a set of material pointers */
  std::map<int, Material*> materials_set;

  /* Check all unique materials */
  for (std::map<int, Material*>::iterator it = _limit_materials.begin();
          it != _limit_materials.end(); ++it) {

    /* Get the material */
    int id = it->first;
    Material* material = it->second;
    Material* original_material = _original_materials[id];

    /* Extract cross-sections */
    FP_PRECISION* sigma_t = material->getSigmaT();
    FP_PRECISION* sigma_s = material->getSigmaS();

    /* Extract cross-sections */
    FP_PRECISION* original_sigma_t = original_material->getSigmaT();
    FP_PRECISION* original_sigma_s = original_material->getSigmaS();

    /* Loop over all energy groups */
    for (int e=0; e < _num_groups; e++) {
      sigma_t[e] = original_sigma_t[e];
      sigma_s[e*_num_groups+e] = original_sigma_s[e*_num_groups+e];
    }
  }
  log_printf(NORMAL, "Material re-set complete");
}


/**
 * @brief Instructs MOC to limit negative cross-sections for early iterations.
 * @param material_ids The material IDs of the cross-sections to limit
 * @param reset_iteration The iteration to reset cross-sections to their
 *        defaults
 */
void Solver::setLimitingXSMaterials(std::vector<int> material_ids,
                                    int reset_iteration) {
  _limit_xs_materials = material_ids;
  _reset_iteration = reset_iteration;
  _limit_xs = true;
}


/**
 * @brief Limits cross-sections so that there are no negative cross-sections.
 * @details A copy of the original cross-section is saved
 */
void Solver::limitXS() {

  log_printf(NORMAL, "Limiting negative cross-sections in %d materials",
             _limit_xs_materials.size());
  std::map<int, Material*> all_materials = _geometry->getAllMaterials();
  for (int i=0; i < _limit_xs_materials.size(); i++) {
    int mat_id = _limit_xs_materials.at(i);
    Material* material = all_materials[mat_id];
    Material* material_copy = material->clone();
    _original_materials[mat_id] = material_copy;
    _limit_materials[mat_id] = material;
    FP_PRECISION* scattering_matrix = material->getSigmaS();
    FP_PRECISION* sigma_t = material->getSigmaT();
    for (int e=0; e < _num_groups; e++) {
      double scattering_value = scattering_matrix[e*_num_groups+e];
      if (scattering_value < 0.0) {
        scattering_matrix[e*_num_groups+e] = 0.0;
        sigma_t[e] -= scattering_value;
      }
    }
  }
  log_printf(NORMAL, "Cross-section adjustment complete");
}


/**
 * @brief All material cross-sections in the geometry are checked for
 *        consistency.
 * @details Each cross-section is checked to ensure that the total
 *          cross-section is greater than or equal to the scattering
 *          cross-section for each energy group and that all cross-sections
 *          are positive.
 */
void Solver::checkXS() {

  log_printf(NORMAL, "Checking material cross-sections");

  /* Create a set of material pointers */
  std::set<Material*> materials_set;

  /* Get a set of the materials over all FSR */
  logLevel level = _xs_log_level;
#pragma omp parallel for
  for (long r=0; r < _num_FSRs; r++) {

    /* Get the material */
    Material* material = _FSR_materials[r];

    /* Check to see that this material hasn't been checked yet */
    if (materials_set.find(material) == materials_set.end()) {
#pragma omp critical
      {
        if (materials_set.find(material) == materials_set.end())
          materials_set.insert(material);
      }
    }
  }

  /* Check all unique materials */
  for (std::set<Material*>::iterator it = materials_set.begin();
          it != materials_set.end(); ++it) {

    /* Get the material */
    Material* material = *it;

    /* Extract cross-sections */
    char* name = material->getName();
    FP_PRECISION* sigma_t = material->getSigmaT();
    FP_PRECISION* sigma_f = material->getSigmaF();
    FP_PRECISION* nu_sigma_f = material->getNuSigmaF();
    FP_PRECISION* scattering_matrix = material->getSigmaS();
    FP_PRECISION* chi = material->getChi();

    /* Loop over all energy groups */
    for (int e=0; e < _num_groups; e++) {

      /* Check that the total cross-section is greater than or equal to the
         scattering cross-section */
      FP_PRECISION sigma_s = 0.0;
      for (int g=0; g < _num_groups; g++) {
        sigma_s += scattering_matrix[g*_num_groups+e];
        if (scattering_matrix[g*_num_groups+e] < 0)
          log_printf(level, "Negative scattering cross-section encountered "
                     "in material ID %d", material->getId());
      }
      if (sigma_s > sigma_t[e]) {
        if (_correct_xs) {
          log_printf(WARNING, "Invalid cross-sections encountered. The "
                     "scattering cross-section has value %6.4f which is "
                     "greater than the total cross-section of value %6.4f in"
                     " material ID %d for group %d", sigma_s, sigma_t[e],
                     material->getId(), e);
          sigma_t[e] = sigma_s;
          log_printf(WARNING, "The total cross-section has been corrected to "
                     " %6.4f in material ID %d for group %d", sigma_s,
                     material->getId(), e);
        }
        else {
          log_printf(level, "Invalid cross-sections encountered. The "
                     "scattering cross-section has value %6.4f which is "
                     "greater than the total cross-section of value %6.4f in"
                     " material ID %d for group %d", sigma_s, sigma_t[e],
                     material->getId(), e);
        }
      }

      /* Check for negative cross-section values */
      if (sigma_t[e] < 0 || sigma_f[e] < 0 || nu_sigma_f[e] < 0 || chi[e] < 0)
        log_printf(level, "Negative cross-section encountered in material "
                   "ID %d", material->getId());
    }
  }
  log_printf(NORMAL, "Material cross-section checks complete");
}


/**
 * @brief Initializes most components of Solver. Mostly needed from the
 *        Python side.
 */
void Solver::initializeSolver(solverMode solver_mode) {

  initializeFSRs();
  initializeSourceArrays();
  initializeExpEvaluators();
  initializeMaterials(solver_mode);
  initializeFluxArrays();
  countFissionableFSRs();
  zeroTrackFluxes();
}


/**
 * @brief Assigns fixed sources assigned by Cell, Material to FSRs.
 */
void Solver::initializeFixedSources() {

  log_printf(INFO, "Transferring fixed sources from cells/materials to FSRs");

  Cell* fsr_cell;
  Material* fsr_material;
  int group;
  FP_PRECISION source;
  std::pair<Cell*, int> cell_group_key;
  std::pair<Material*, int> mat_group_key;
  std::map< std::pair<Cell*, int>, FP_PRECISION >::iterator cell_iter;
  std::map< std::pair<Material*, int>, FP_PRECISION >::iterator mat_iter;

  /* Fixed sources assigned by Cell */
  for (cell_iter = _fix_src_cell_map.begin();
       cell_iter != _fix_src_cell_map.end(); ++cell_iter) {

    /* Get the Cell with an assigned fixed source */
    cell_group_key = cell_iter->first;
    group = cell_group_key.second;
    source = _fix_src_cell_map[cell_group_key];

    /* Search for this Cell in all FSRs */
    for (long r=0; r < _num_FSRs; r++) {
      fsr_cell = _geometry->findCellContainingFSR(r);
      if (cell_group_key.first->getId() == fsr_cell->getId())
        setFixedSourceByFSR(r, group, source);
    }
  }

  /* Fixed sources assigned by Material */
  for (mat_iter = _fix_src_material_map.begin();
       mat_iter != _fix_src_material_map.end(); ++mat_iter) {

    /* Get the Material with an assigned fixed source */
    mat_group_key = mat_iter->first;
    group = mat_group_key.second;
    source = _fix_src_material_map[mat_group_key];

    /* Search for this Material in all FSRs */
    for (long r=0; r < _num_FSRs; r++) {
      fsr_material = _geometry->findFSRMaterial(r);
      if (mat_group_key.first->getId() == fsr_material->getId())
        setFixedSourceByFSR(r, group, source);
    }
  }
}


/**
 * @brief Initializes a Cmfd object for acceleration prior to source iteration.
 * @details Instantiates a dummy Cmfd object if one was not assigned to
 *          the Solver by the user and initializes FSRs, materials, fluxes
 *          and the Mesh object. This method is for internal use only
 *          and should not be called directly by the user.
 */
void Solver::initializeCmfd() {

  log_printf(INFO_ONCE, "Initializing CMFD...");

  /* Retrieve CMFD from the Geometry */
  _cmfd = _geometry->getCmfd();

  /* If the user did not initialize Cmfd, simply return */
  if (_cmfd == NULL)
    return;
  else if (!_cmfd->isFluxUpdateOn())
    return;

  /* Initialize the CMFD energy group structure */
  _cmfd->setSourceConvergenceThreshold(_converge_thresh*1.e-1); //FIXME
  _cmfd->setNumMOCGroups(_num_groups);
  _cmfd->initializeGroupMap();

  /* Give CMFD number of FSRs and FSR property arrays */
  _cmfd->setSolve3D(_SOLVE_3D);
  _cmfd->setNumFSRs(_num_FSRs);
  _cmfd->setFSRVolumes(_FSR_volumes);
  _cmfd->setFSRMaterials(_FSR_materials);
  _cmfd->setFSRFluxes(_scalar_flux);
  _cmfd->setFSRSources(_reduced_sources);
  _cmfd->setQuadrature(_quad);
  _cmfd->setGeometry(_geometry);
  _cmfd->setAzimSpacings(_quad->getAzimSpacings(), _num_azim);
  if (!_is_restart)
    _cmfd->initialize();


  TrackGenerator3D* track_generator_3D =
    dynamic_cast<TrackGenerator3D*>(_track_generator);
  if (track_generator_3D != NULL)
    _cmfd->setPolarSpacings(_quad->getPolarSpacings(), _num_azim, _num_polar);
}


/**
 * @brief Performs a spectrum calculation to update the scalar fluxes.
 * @details This function is meant to be used before transport sweeps in an
 *          eigenvalue calculation in order to gain a better initial guess
 *          on the flux shape. It is equivalent to performing a CMFD update
 *          with no current tallies (pure diffusion solve) across a coarse
 *          mesh with one mesh cell per domain.
 * @param threshold The convergence threshold of the calculation
 */
void Solver::calculateInitialSpectrum(double threshold) {

  log_printf(NORMAL, "Calculating initial spectrum with threshold %3.2e",
             threshold);

  /* Setup the spectrum calclator as a CMFD solver in MOC group structure */
  Cmfd spectrum_calculator;
  std::vector<std::vector<int> > group_structure;
  group_structure.resize(_num_groups);
  for (int g=0; g < _num_groups; g++)
    group_structure.at(g).push_back(g+1);
  spectrum_calculator.setGroupStructure(group_structure);

  /* Set CMFD settings for the spectrum calculator */
  spectrum_calculator.setSORRelaxationFactor(1.6);
  spectrum_calculator.useFluxLimiting(true);
  spectrum_calculator.setKNearest(1);
  _geometry->initializeSpectrumCalculator(&spectrum_calculator);

  /* If 2D Solve, set z-direction mesh size to 1 and depth to 1.0 */
  if (!_SOLVE_3D) {
    spectrum_calculator.setNumZ(1);
    spectrum_calculator.setBoundary(SURFACE_Z_MIN, REFLECTIVE);
    spectrum_calculator.setBoundary(SURFACE_Z_MAX, REFLECTIVE);
  }

  /* Initialize the energy group structure */
  spectrum_calculator.setSourceConvergenceThreshold(threshold);
  spectrum_calculator.setNumMOCGroups(_num_groups);
  spectrum_calculator.initializeGroupMap();

  /* Give the spectrum calculator the number of FSRs and FSR property arrays */
  spectrum_calculator.setSolve3D(_SOLVE_3D);
  spectrum_calculator.setNumFSRs(_num_FSRs);
  spectrum_calculator.setFSRVolumes(_FSR_volumes);
  spectrum_calculator.setFSRMaterials(_FSR_materials);
  spectrum_calculator.setFSRFluxes(_scalar_flux);
  spectrum_calculator.setFSRSources(_reduced_sources);
  spectrum_calculator.setQuadrature(_quad);
  spectrum_calculator.setAzimSpacings(_quad->getAzimSpacings(), _num_azim);
  spectrum_calculator.initialize();

  TrackGenerator3D* track_generator_3D =
    dynamic_cast<TrackGenerator3D*>(_track_generator);
  if (track_generator_3D != NULL)
    spectrum_calculator.setPolarSpacings(_quad->getPolarSpacings(), _num_azim,
         _num_polar);

  /* Solve the system */
  log_printf(NORMAL, "Computing K-eff");
  _k_eff = spectrum_calculator.computeKeff(0);
  log_printf(NORMAL, "Normalizing Fluxes");
  normalizeFluxes();

  /* Copy k-eff to CMFD solver if applicable */
  if (_cmfd != NULL)
    _cmfd->setKeff(_k_eff);

  log_printf(NORMAL, "Calculated initial spectrum with k-eff = %6.6f", _k_eff);
}


/**
 * @brief Returns the Material data to its original state.
 * @details In an adjoint calculation, the scattering and fission matrices
 *          in each material are transposed during initialization. This
 *          routine returns both matrices to their original (FORWARD)
 *          state at the end of a calculation.
 * @param mode the solution type (FORWARD or ADJOINT)
 */
void Solver::resetMaterials(solverMode mode) {

  if (mode == FORWARD)
    return;

  log_printf(INFO, "Resetting materials...");

  std::map<int, Material*> materials = _geometry->getAllMaterials();
  std::map<int, Material*>::iterator m_iter;

  for (m_iter = materials.begin(); m_iter != materials.end(); ++m_iter)
    m_iter->second->transposeProductionMatrices();
}


/**
 * @brief This method performs one transport sweep using the fission source.
 * @details This is a helper routine used for Krylov subspace methods.
 */
void Solver::fissionTransportSweep() {
  computeFSRFissionSources();
  transportSweep();
  addSourceToScalarFlux();
}


/**
 * @brief This method performs one transport sweep using the scatter source.
 * @details This is a helper routine used for Krylov subspace methods.
 */
void Solver::scatterTransportSweep() {
  computeFSRScatterSources();
  transportSweep();
  addSourceToScalarFlux();
}


/**
 * @brief Computes the scalar flux distribution by performing a series of
 *        transport sweeps.
 * @details This is the main method exposed to the user through the Python
 *          interface to compute the scalar flux distribution, e.g., for a
 *          fixed source calculation. This routine makes an initial guess for
 *          scalar and boundary fluxes and performs transport sweep until
 *          convergence.
 *
 *          By default, this method will perform a maximum of 1000 transport
 *          sweeps with a 1E-5 threshold on the average FSR scalar flux. These
 *          values may be freely modified by the user at runtime.
 *
 *          The only_fixed_source runtime parameter may be used to control
 *          the type of source distribution used in the calculation. By
 *          default, this paramter is true and only the fixed sources specified
 *          by the user will be considered. Alternatively, when the parameter
 *          is false, the source will be computed as the scattering and fission
 *          sources resulting from a previously computed flux distribution
 *          (e.g., an eigenvalue calculation) in addition to any user-defined
 *          fixed sources.
 *
 *          This method may be called by the user to compute the scalar flux
 *          for a fixed source distribution from Python as follows:
 *
 * @code
 *          // Assign fixed sources
 *          // ...
 *
 *          // Find the flux distribution resulting from the fixed sources
 *          solver.computeFlux(max_iters=100)
 * @endcode
 *
 *          Alternatively, as described above, this method may be called by
 *          the user in Python to compute the flux from a superposition of
 *          fixed and / or eigenvalue sources as follows:
 *
 * @code
 *          // Solve for sources and scalar flux distribution
 *          solver.computeEigenvalue(max_iters=1000)
 *
 *          // Add fixed source(s)
 *          // ...
 *
 *          // Find fluxes from superposition of eigenvalue and fixed sources
 *          solver.computeFlux(max_iters=100, only_fixed_source=False)
 * @endcode
 *
 *
 * @param max_iters the maximum number of source iterations to allow
 * @param only_fixed_source use only fixed sources (true by default)
 */
void Solver::computeFlux(int max_iters, bool only_fixed_source) {

  if (_track_generator == NULL)
    log_printf(ERROR, "The Solver is unable to compute the flux "
               "since it does not contain a TrackGenerator");

  log_printf(NORMAL, "Computing the flux...");

  /* Clear all timing data from a previous simulation run */
  clearTimerSplits();

  /* Start the timers to record the total solve and initialization times */
  _timer->startTimer();
  _timer->startTimer();

  /* Initialize keff to 1 for FSR source calculations */
  _k_eff = 1.;

  double residual = 0.;

  /* Initialize data structures */
  initializeMaterials(_solver_mode);
  initializeFSRs();
  countFissionableFSRs();
  initializeSourceArrays();
  initializeExpEvaluators();

  /* Initialize new flux arrays if a) the user requested the use of
   * only fixed sources or b) no previous simulation was performed which
   * initialized and computed the flux (e.g., an eigenvalue calculation) */
  if (only_fixed_source || _num_iterations == 0) {
    if (!_is_restart)
      initializeFluxArrays();
    flattenFSRFluxes(0.);
    storeFSRFluxes();
  }

  /* Compute the sum of fixed, total and scattering sources */
  computeFSRSources(0);

  /* Stop timer for solver initialization */
  _timer->stopTimer();
  _timer->recordSplit("Solver initialization");

  /* Source iteration loop */
  for (int i=0; i < max_iters; i++) {

    transportSweep();
    addSourceToScalarFlux();
    residual = computeResidual(SCALAR_FLUX);
    storeFSRFluxes();

    log_printf(NORMAL, "Iteration %d:\tres = %1.3E", i, residual);

    /* Check for convergence of the fission source distribution */
    if (i > 1 && residual < _converge_thresh) {
      _num_iterations = i;
      _timer->stopTimer();
      _timer->recordSplit("Total time");
      return;
    }
  }

  log_printf(WARNING, "Unable to converge the flux");

  _num_iterations = max_iters;
  _timer->stopTimer();
  _timer->recordSplit("Total time");
}


/**
 * @brief Computes the total source distribution by performing a series of
 *        transport sweep and source updates.
 * @details This is the main method exposed to the user through the Python
 *          interface to compute the source distribution, e.g., for a fixed
 *          and/or external source calculation. This routine makes an initial
 *          guess for the scalar and boundary fluxes and performs transport
 *          sweeps and source updates until convergence.
 *
 *          By default, this method will perform a maximum of 1000 transport
 *          sweeps with a 1E-5 threshold on the integrated FSR total source.
 *          These values may be freely modified by the user at runtime.
 *
 *          The k_eff parameter may be used for fixed source calculations
 *          with fissionable material (e.g., start-up in a reactor from
 *          a fixed external source). In this case, the user must "guess"
 *          the critical eigenvalue to be be used to scale the fission source.
 *
 *          The res_type parameter may be used to control the convergence
 *          criterion - SCALAR_FLUX, TOTAL_SOURCE (default) and FISSION_SOURCE
 *          are all supported options in OpenMOC at this time.
 *
 *          This method may be called by the user from Python as follows:
 *
 * @code
 *          // Assign fixed sources
 *          // ...
 *
 *          // Find the source distribution resulting from the fixed sources
 *          solver.computeSource(max_iters=100, k_eff=0.981)
 * @endcode
 *
 * @param max_iters the maximum number of source iterations to allow
 * @param k_eff the sub/super-critical eigenvalue (default 1.0)
 * @param res_type the type of residual used for the convergence criterion
 */
void Solver::computeSource(int max_iters, double k_eff, residualType res_type) {

  if (_track_generator == NULL)
    log_printf(ERROR, "The Solver is unable to compute the source "
               "since it does not contain a TrackGenerator");

  else if (k_eff <= 0.)
    log_printf(ERROR, "The Solver is unable to compute the source with "
               "keff = %f since it is not a positive value", k_eff);

  log_printf(NORMAL, "Computing the source...");

  /* Clear all timing data from a previous simulation run */
  clearTimerSplits();

  _k_eff = k_eff;
  double residual = 0.;

  /* Initialize data structures */
  initializeMaterials(_solver_mode);
  initializeFSRs();
  initializeExpEvaluators();
  if (!_is_restart)
    initializeFluxArrays();
  initializeSourceArrays();

  /* Compute a starting guess for the fluxes */
  computeInitialFluxGuess(true);

  /* Start the timer to record the total time to converge the flux */
  _timer->startTimer();

  /* Source iteration loop */
  for (int i=0; i < max_iters; i++) {

    computeFSRSources(i);
    transportSweep();
    addSourceToScalarFlux();
    residual = computeResidual(res_type);
    storeFSRFluxes();

    log_printf(NORMAL, "Iteration %d:\tres = %1.3E", i, residual);

    /* Check for convergence of the fission source distribution */
    if (i > 1 && residual < _converge_thresh) {
      _num_iterations = i;
      _timer->stopTimer();
      _timer->recordSplit("Total time");
      return;
    }
  }

  log_printf(WARNING, "Unable to converge the source");

  _num_iterations = max_iters;
  _timer->stopTimer();
  _timer->recordSplit("Total time");
}


/**
 * @brief Computes keff by performing a series of transport sweep and
 *        source updates.
 * @details This is the main method exposed to the user through the Python
 *          interface to perform an eigenvalue calculation. The method makes
 *          an initial guess for the scalar and boundary fluxes and performs
 *          transport sweeps and source updates until convergence.
 *
 *          By default, this method will perform a maximum of 1000 transport
 *          sweeps with a 1E-5 threshold on the integrated FSR fission source.
 *          These values may be freely modified by the user at runtime.
 *
 *          The res_type parameter may be used to control the convergence
 *          criterion - SCALAR_FLUX, TOTAL_SOURCE and FISSION_SOURCE (default)
 *          are all supported options in OpenMOC at this time.
 *
 * @code
 *          solver.computeEigenvalue(max_iters=100, res_type=FISSION_SOURCE)
 * @endcode
 *
 * @param max_iters the maximum number of source iterations to allow
 * @param res_type the type of residual used for the convergence criterion
 */
void Solver::computeEigenvalue(int max_iters, residualType res_type) {

  if (_track_generator == NULL)
    log_printf(ERROR, "The Solver is unable to compute the eigenvalue "
               "since it does not contain a TrackGenerator");

  log_printf(NORMAL, "Initializing MOC eigenvalue solver...");

  /* Clear all timing data from a previous simulation run */
  clearTimerSplits();

  /* Reset number of iterations, start at 1 if restarting */
  _num_iterations = _is_restart;

  /* Start the timers to record the total solve and initialization times */
  _timer->startTimer();
  _timer->startTimer();

  /* Clear convergence data from a previous simulation run */
  double previous_residual = 1.0;
  double residual = 0.;
  if (!_is_restart)
    _k_eff = 1.;

  /* Initialize data structures */
  initializeMaterials(_solver_mode);
  initializeFSRs();
  countFissionableFSRs();
  initializeExpEvaluators();
  if (!_is_restart)
    initializeFluxArrays();
  initializeSourceArrays();
  initializeCmfd();
  _geometry->fixFSRMaps();
#ifdef MPIx
  if (_geometry->isDomainDecomposed())
    MPI_Barrier(_geometry->getMPICart());
#endif
  printInputParamsSummary();

  /* Print memory report */
#ifdef BGQ
  printBGQMemory();
#endif

  /* Compute a starting guess for the fluxes */
  computeInitialFluxGuess();

#ifdef MPIx
  if (_geometry->isDomainDecomposed())
    MPI_Barrier(_geometry->getMPICart());
#endif

  /* Create object to track convergence data if requested */
  ConvergenceData convergence_data;
  if (_verbose && _cmfd != NULL) {
    _cmfd->setConvergenceData(&convergence_data);
    log_printf(NORMAL, "iter   k-eff   eps-k  eps-MOC   D.R.   "
               "eps-FS1   eps-FSN   #FS  eps-flux1 eps-fluxN"
               "  #FX1 #FXN  MAX P.F.");
  }

  /* Stop timer for solver initialization */
  _timer->stopTimer();
  _timer->recordSplit("Solver initialization");

  log_printf(NORMAL, "Computing the eigenvalue...");

  /* Record the starting eigenvalue guess */
  double k_prev = _k_eff;

  /* Source iteration loop */
  for (int i=0; i < max_iters; i++) {

    /* Compute the stabilizing flux if necessary */
    if (i > 0 && _stabilize_transport) {
      computeStabilizingFlux();
    }

    /* Perform the source iteration */
    computeFSRSources(i);
    transportSweep();
    addSourceToScalarFlux();

    /* Solve CMFD diffusion problem and update MOC flux */
    if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
      _k_eff = _cmfd->computeKeff(_num_iterations);
    else
      computeKeff();

    /* Apply the flux adjustment if transport stabilization is on */
    if (i > 0 && _stabilize_transport) {
      stabilizeFlux();
    }

    /* Normalize the flux and compute residuals */
    normalizeFluxes();
    residual = computeResidual(res_type);

    /* Compute difference in k and apparent dominance ratio */
    double dr = residual / previous_residual;
    int dk = 1e5 * (_k_eff - k_prev);
    previous_residual = residual;
    k_prev = _k_eff;

    /* Ouptut iteration report */
    if (_verbose && _cmfd != NULL) {

      /* Unpack convergence data */
      double pf = convergence_data.pf;
      double cmfd_res_1 = convergence_data.cmfd_res_1;
      double cmfd_res_end = convergence_data.cmfd_res_end;
      double linear_res_1 = convergence_data.linear_res_1;
      double linear_res_end = convergence_data.linear_res_end;
      int cmfd_iters = convergence_data.cmfd_iters;
      int linear_iters_1 = convergence_data.linear_iters_1;
      int linear_iters_end = convergence_data.linear_iters_end;
      log_printf(NORMAL, "%3d  %1.6f  %5d  %1.6f  %1.3f  %1.6f  %1.6f"
                 "  %3d  %1.6f  %1.6f  %3d  %3d    %.4e", i, _k_eff,
                 dk, residual, dr, cmfd_res_1, cmfd_res_end,
                 cmfd_iters, linear_res_1, linear_res_end,
                 linear_iters_1, linear_iters_end, pf);
    }
    else {
      log_printf(NORMAL, "Iteration %d:  k_eff = %1.6f   "
                 "res = %1.3E  delta-k (pcm) = %d D.R. = %1.4f", i, _k_eff,
                 residual, dk, dr);
    }

    if (_cmfd != NULL) {
      if (residual <= 0)
        residual = 1e-6;
      _cmfd->setSourceConvergenceThreshold(0.01*residual);
    }
    storeFSRFluxes();
    _num_iterations++;

    /* Check for convergence of the fission source distribution */
    if (residual < _converge_thresh && std::abs(dk) < 1)
      break;
  }

  if (_num_iterations == max_iters)
    log_printf(WARNING, "Unable to converge the source distribution");

  _timer->stopTimer();
  _timer->recordSplit("Total time");
}


/**
 * @brief Load or compute a starting guess for scalar fluxes.
 * @details By default, OpenMOC assumes an initial flux guess flat in space and
 *          energy. Other options are available:
 *          - reference fluxes can be loaded as an initial guess. They will
 *            also be used to compute the residuals
 *          - a chi-spectrum from a fissile material can be assumed. Only
 *            fissile region will have a non-zero initial flux guess. Note
 *            that this may induce 0 reaction rates in CMFD cells, and startup
 *            unaccelerated iterations may be required
 *          - scalar fluxes can be loaded from a flux file, to perform a
 *            restart calculation. The angular fluxes need to converged
 *            again since they are not saved in the flux file
 *          - a 1-cell 1-group CMFD calculation can be run to compute an
 *            initial spectrum guess
 * @param is_source_computation whether the solver is computing an eigenvalue
 *        or a source distribution
 */
void Solver::computeInitialFluxGuess(bool is_source_computation) {

  /* Load reference solution if necessary */
  if (_calculate_residuals_by_reference) {
    loadFSRFluxes(_reference_file, false);
    long size = _num_FSRs * _num_groups;
    _reference_flux = new FP_PRECISION[size];
    memcpy(_reference_flux, _scalar_flux, size * sizeof(FP_PRECISION));
  }

  /* Guess flat spatial scalar flux for each region */
  if (!_is_restart) {
    if (_chi_spectrum_material == NULL)
      flattenFSRFluxes(1.0);
    else
      flattenFSRFluxesChiSpectrum();

    /* Normalize flux guess for eigenvalue computations */
    if (!is_source_computation)
      normalizeFluxes();
  }
  storeFSRFluxes();

  /* Load initial FSR fluxes from file if requested */
  if (_load_initial_FSR_fluxes) {
    loadFSRFluxes(_initial_FSR_fluxes_file, true);
    normalizeFluxes();
    storeFSRFluxes();

#ifdef MPIx
    if (_geometry->isDomainDecomposed())
      MPI_Barrier(_geometry->getMPICart());
#endif

    /* Perform startup sweeps to converge angular fluxes */
    int startup_iterations = 30;
    computeFSRSources(0);
    for (int i=0; i < startup_iterations; i++) {
      log_printf(NORMAL, "Startup sweep %d / %d", i, startup_iterations);
      transportSweep();
    }
    addSourceToScalarFlux();
  }

  /* Perform initial spectrum calculation if requested */
  if (_calculate_initial_spectrum)
    calculateInitialSpectrum(_initial_spectrum_thresh);
}


/**
 * @brief Deletes the Timer's timing entries for each timed code section
 *        code in the source convergence loop.
 */
void Solver::clearTimerSplits() {

  /* Solver and sweep timers */
  _timer->clearSplit("Total time");
  _timer->clearSplit("Solver initialization");
  _timer->clearSplit("Transport Sweep");
#ifdef MPIx
  _timer->clearSplit("Total transfer time");
  _timer->clearSplit("Packing time");
  _timer->clearSplit("Communication time");
  _timer->clearSplit("Unpacking time");
  _timer->clearSplit("Idle time");
#endif

  /* CMFD timers */
  _timer->clearSplit("Total CMFD time");
  _timer->clearSplit("Total collapse time");
  _timer->clearSplit("Matrix construction time");
#ifdef MPIx
  _timer->clearSplit("CMFD MPI communication time");
#endif
  _timer->clearSplit("Total solver time");
  _timer->clearSplit("Total MOC flux update time");
}


/**
 * @brief Sets the solver to print extra information for each iteration
 */
void Solver::setVerboseIterationReport() {
  set_line_length(120);
  _verbose = true;
}


/**
 * @brief Prints a report of the timing statistics to the console.
 */
void Solver::printTimerReport() {

  std::string msg_string;

  /* Collapse timer to average values in domain decomposition */
#ifdef MPIx
  if (_geometry->isDomainDecomposed())
    _timer->reduceTimer(_geometry->getMPICart());
#endif

  log_printf(TITLE, "TIMING REPORT");

  /* Print track generation time */
  _track_generator->printTimerReport(false);

  /* Get the total runtime */
  double tot_time = _timer->getSplit("Total time");
  msg_string = "Total time to solution";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), tot_time);

  /* Solver initialization */
  double initialization = _timer->getSplit("Solver initialization");
  msg_string = "  Solver initialization";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), initialization);

  /* Transport sweep */
  double transport_sweep = _timer->getSplit("Transport Sweep");
  msg_string = "  Transport Sweep";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), transport_sweep);

#ifdef MPIx
  /* Boundary track angular fluxes transfer */
  double transfer_time = _timer->getSplit("Total transfer time");
  msg_string = "  Angular Flux Transfer";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), transfer_time);

  /* Boundary track angular fluxes packing into buffers */
  double pack_time = _timer->getSplit("Packing time");
  msg_string = "    Angular Flux Packing Time";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), pack_time);

  /* Communication of track angular fluxes */
  double comm_time = _timer->getSplit("Communication time");
  msg_string = "    Angular Flux Communication Time";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), comm_time);

  /* Boundary track angular fluxes packing into buffers */
  double unpack_time = _timer->getSplit("Unpacking time");
  msg_string = "    Angular Flux Unpacking Time";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), unpack_time);

  /* Idle time between transport sweep and angular fluxes transfer */
  double idle_time = _timer->getSplit("Idle time");
  msg_string = "  Total Idle Time Between Sweep and Transfer";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), idle_time);
#endif

  /* CMFD acceleration time */
  if (_cmfd != NULL)
    _cmfd->printTimerReport();

  /* Time per segment */
  long num_segments = 0;
  TrackGenerator3D* track_generator_3D =
    dynamic_cast<TrackGenerator3D*>(_track_generator);
  if (track_generator_3D != NULL)
    num_segments = track_generator_3D->getNum3DSegments();
  else
    num_segments = _track_generator->getNum2DSegments();

  /* Reduce number of Tracks and segments if necessary */
  long total_num_segments = num_segments;
  long total_num_tracks = _tot_num_tracks;
#ifdef MPIx
  if (_geometry->isDomainDecomposed()) {
    MPI_Comm MPI_cart = _geometry->getMPICart();
    MPI_Allreduce(&num_segments, &total_num_segments, 1, MPI_LONG,
                  MPI_SUM, MPI_cart);
    MPI_Allreduce(&_tot_num_tracks, &total_num_tracks, 1, MPI_LONG,
                  MPI_SUM, MPI_cart);
  }
#endif

  int num_ranks = 1;
#ifdef MPIx
  if (_geometry->isDomainDecomposed()) {
    MPI_Comm MPI_cart = _geometry->getMPICart();
    MPI_Comm_size(MPI_cart, &num_ranks);
  }
#endif

  long num_integrations = 2 * _fluxes_per_track * total_num_segments *
      _num_iterations;
  double time_per_integration = transport_sweep / num_integrations *
                                (omp_get_max_threads() * num_ranks);
  double time_per_integ_total = tot_time / num_integrations *
                                (omp_get_max_threads() * num_ranks);

  if (!_gpu_solver) {
    msg_string = "Integration time by segment-group-thread (sweep)";
    msg_string.resize(53, '.');
    log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(),
               time_per_integration);

    msg_string = "Integration time by segment-group-thread (total)";
    msg_string.resize(53, '.');
    log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(),
               time_per_integ_total);
  }
  else {
    msg_string = "Integration time by segment-group-device (sweep)";
    msg_string.resize(53, '.');
    log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), time_per_integration
               / omp_get_max_threads());

    msg_string = "Integration time by segment-group-device (total)";
    msg_string.resize(53, '.');
    log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), time_per_integ_total
               / omp_get_max_threads());
  }

  /* Print footer with number of tracks, segments and fsrs */
  set_separator_character('-');
  log_printf(SEPARATOR, "-");

  msg_string = "           # tracks          # segments          # FSRs";
  log_printf(RESULT, "%s", msg_string.c_str());
  log_printf(SEPARATOR, "-");

  int num_digits = (int) log10((double) total_num_tracks);
  num_digits += (int) log10((double) total_num_segments);
  num_digits += (int) log10((double) _geometry->getNumTotalFSRs());

  num_digits = 66 - num_digits;
  num_digits /= 4;

  std::stringstream msg;

  for (int i=0; i < 4; i++) {
    for (int j=0; j < num_digits; j++)
      msg << " ";

    if (i == 0)
      msg << total_num_tracks;
    else if (i == 1)
      msg << total_num_segments;
    else if (i == 2)
      msg << _geometry->getNumTotalFSRs();
  }

  log_printf(RESULT, "%s", msg.str().c_str());
  log_printf(SEPARATOR, "-");
}


/**
 * @brief Prints fission rates to a binary file.
 * @param fname the name of the file to dump the fission rates to
 * @param nx number of mesh cells in the x-direction
 * @param ny number of mesh cells in the y-direction
 * @param nz number of mesh cells in the z-direction
 */
void Solver::printFissionRates(std::string fname, int nx, int ny, int nz) {

  Universe* root_universe = _geometry->getRootUniverse();
  double x_min = root_universe->getMinX();
  double x_max = root_universe->getMaxX();
  double y_min = root_universe->getMinY();
  double y_max = root_universe->getMaxY();
  double z_min = root_universe->getMinZ();
  double z_max = root_universe->getMaxZ();

  double* fission_rates = new double[nx*ny*nz];
  for (int i=0; i < nx*ny*nz; i++)
    fission_rates[i] = 0;

  int num_fsrs = _geometry->getNumTotalFSRs();
  double* fsr_fission_rates = new double[num_fsrs];
  computeFSRFissionRates(fsr_fission_rates, num_fsrs);

  for (long r=0; r < num_fsrs; r++) {

    std::vector<double> pt = _geometry->getGlobalFSRCentroidData(r);

    int x_ind = nx * (pt.at(0) - x_min) / (x_max - x_min);
    int y_ind = ny * (pt.at(1) - y_min) / (y_max - y_min);
    int z_ind = nz * (pt.at(2) - z_min) / (z_max - z_min);

    int ind = z_ind * nx * ny + y_ind * nx + x_ind;

    fission_rates[ind] += fsr_fission_rates[r];
  }

  int rank = 0;
#ifdef MPIx
  if (_geometry->isDomainDecomposed()) {
    MPI_Comm comm = _geometry->getMPICart();
    MPI_Comm_rank(comm, &rank);
  }
#endif

  if (rank == 0) {
    std::ofstream out(fname.c_str());
    out << "Fission rates for (" << nx << ", " << ny << ", " << nz <<
      ")" << std::endl;
    for (int i=0; i < nx; i++) {
      for (int j=0; j < ny; j++) {
        for (int k=0; k < nz; k++) {
          int ind = k * nx * ny + j * nx + i;
          out << "Region " << i << ", " << j << ", " << k << " at point " <<
            "(" << x_min + (i+0.5) * (x_max - x_min) << ", " <<
            y_min + (j+0.5) * (y_max - y_min) << ", " <<
            z_min + (k+0.5) * (z_max - z_min) << ") -> " <<
            fission_rates[ind] << std::endl;
        }
      }
    }
  }
  delete [] fission_rates;
  delete [] fsr_fission_rates;
}


/**
 * @brief A function that returns the array of scalar fluxes.
 * @return The scalar fluxes
 */
FP_PRECISION* Solver::getFluxesArray() {
  return _scalar_flux;
}


/**
 * @brief Sets computation method of k-eff from fission, absorption, and leakage
 *        rates rather than from fission rates.
 *        keff = fission/(absorption + leakage)
 */
void Solver::setKeffFromNeutronBalance() {
  _keff_from_fission_rates = false;
}


/**
 * @brief Sets residuals to be computed a error relative to a reference.
 * @param fname The file containing the flux solution of the reference
 */
void Solver::setResidualByReference(std::string fname) {
  _calculate_residuals_by_reference = true;
  _reference_file = fname;
}


/**
 * @brief Prints scalar fluxes to a binary file.
 * @param fname the name of the file to dump the fluxes to
 */
void Solver::dumpFSRFluxes(std::string fname) {

  /* Determine the FSR fluxes file name */
  std::string filename = fname;
  if (_geometry->isDomainDecomposed()) {
    int indexes[3];
    if (_geometry->isRootDomain())
      mkdir(filename.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#ifdef MPIx
    MPI_Barrier(_geometry->getMPICart());
#endif
    filename += "/node";
    _geometry->getDomainIndexes(indexes);
    for (int i=0; i < 3; i++) {
      long long int num = indexes[i];
      std::string str = std::to_string(num);
      filename += "_" + str;
    }
  }

  /* Write the FSR fluxes file */
  FILE* out;
  out = fopen(filename.c_str(), "w");
  if (out == NULL)
    log_printf(ERROR, "Fluxes file %s could not be written.", &filename[0]);

  /* Write k-eff */
  fwrite(&_k_eff, sizeof(double), 1, out);

  /* Write number of energy groups */
  fwrite(&_num_groups, sizeof(int), 1, out);

  /* Write number of energy groups */
  fwrite(&_num_FSRs, sizeof(long), 1, out);

  /* Write centroid and flux data */
  for (long r=0; r < _num_FSRs; r++) {
    Point* centroid = _geometry->getFSRCentroid(r);
    double x = centroid->getX();
    double y = centroid->getY();
    double z = centroid->getZ();
    fwrite(&x, sizeof(double), 1, out);
    fwrite(&y, sizeof(double), 1, out);
    fwrite(&z, sizeof(double), 1, out);
    fwrite(&_scalar_flux[r*_num_groups], sizeof(double), _num_groups, out);
  }
  fclose(out);
}


/**
 * @brief Load the initial scalar flux distribution from a binary file.
 * @param fname The file containing the scalar fluxes
 */
void Solver::loadInitialFSRFluxes(std::string fname) {
  _initial_FSR_fluxes_file = fname;
  _load_initial_FSR_fluxes = true;
}


/**
 * @brief Load scalar fluxes from a binary file.
 * @details The matching source regions between the current calculation and
 *          those in the loaded file are determined by comparing centroids
 * @param fname The file containing the scalar fluxes
 * @param assign_k_eff Whether to set k-eff to the one loaded in the binary file
 * @param tolerance The width of the region in which to search for the matching
 *        centroid
 */
void Solver::loadFSRFluxes(std::string fname, bool assign_k_eff,
                           double tolerance) {

  /* Determine the FSR fluxes file name */
  std::string filename = fname;
  if (_geometry->isDomainDecomposed()) {
    int indexes[3];
    filename += "/node";
    _geometry->getDomainIndexes(indexes);
    for (int i=0; i < 3; i++) {
      long long int num = indexes[i];
      std::string str = std::to_string(num);
      filename += "_" + str;
    }
  }

  /* Load the FSR fluxes file */
  FILE* in;
  in = fopen(filename.c_str(), "r");
  if (in == NULL)
    log_printf(ERROR, "Failed to find file %s", filename.c_str());
  log_printf(NORMAL, "Reading fluxes from %s", fname.c_str());

  /* Read number of energy groups */
  double k_eff;
  int ret = fread(&k_eff, sizeof(double), 1, in);
  if (assign_k_eff) {
    log_printf(NORMAL, "Loaded k-eff %6.6f", k_eff);
    _k_eff = k_eff;
  }

  /* Read number of energy groups */
  int num_groups;
  ret = fread(&num_groups, sizeof(int), 1, in);

  /* Read number of energy groups */
  long num_FSRs;
  ret = fread(&num_FSRs, sizeof(long), 1, in);

  /* Check that the number of FSRs and the number of groups match */
  if (num_FSRs != _num_FSRs)
    log_printf(ERROR, "The number of FSRs in the current Geometry do not match"
               " the number of FSRs in the binary flux data file");
  if (num_groups != _num_groups)
    log_printf(ERROR, "The number of energy groups in the current Geometry do "
               "not match the number of energy groups in the binary flux data "
               "file");

  /* Setup array structures */
  double* x_coord = new double[num_FSRs];
  double* y_coord = new double[num_FSRs];
  double* z_coord = new double[num_FSRs];
  double* fluxes = new double[num_FSRs * num_groups];

  /* Load data into structures */
  for (long r=0; r < num_FSRs; r++) {
    ret = fread(&x_coord[r], sizeof(double), 1, in);
    ret = fread(&y_coord[r], sizeof(double), 1, in);
    ret = fread(&z_coord[r], sizeof(double), 1, in);
    ret = fread(&fluxes[r*num_groups], sizeof(double), num_groups, in);
  }
  fclose(in);

  /* Setup cell index mapping */
  int* cell_x = new int[num_FSRs];
  int* cell_y = new int[num_FSRs];
  int* cell_z = new int[num_FSRs];
#pragma omp parallel for
  for (long r=0; r < num_FSRs; r++) {
    cell_x[r] = x_coord[r] / tolerance;
    cell_y[r] = y_coord[r] / tolerance;
    cell_z[r] = z_coord[r] / tolerance;
  }
  int* cell_indexes[3] = {cell_x, cell_y, cell_z};

  /* Find min/max indexes in order to make integer-index mapping */
  int min_ind[3] = {cell_x[0], cell_y[0], cell_z[0]};
  int max_ind[3] = {cell_x[0], cell_y[0], cell_z[0]};
  for (long r=0; r < num_FSRs; r++) {
    for (int i=0; i < 3; i++) {
      if (cell_indexes[i][r] < min_ind[i])
        min_ind[i] = cell_indexes[i][r];
      if (cell_indexes[i][r] > max_ind[i])
        max_ind[i] = cell_indexes[i][r];
    }
  }

  /* Create mapping of FSRs to cell indexes */
  std::unordered_map<long, std::vector<long> > hashed_lookup;
  long nx = max_ind[0] - min_ind[0] + 1;
  long ny = max_ind[1] - min_ind[1] + 1;
  long nz = max_ind[2] - min_ind[2] + 1;
  for (long r=0; r < num_FSRs; r++) {
    long index = (cell_z[r] - min_ind[2]) * nx * ny +
                (cell_y[r] - min_ind[1]) * nx + cell_x[r] - min_ind[0];
    if (hashed_lookup.find(index) == hashed_lookup.end())
      hashed_lookup.insert(std::make_pair(index, std::vector<long>()));
    hashed_lookup[index].push_back(r);
  }

  /* Generate centroids if they have not been generated yet */
  double max_centroid_error = 0.0;
  if (!_geometry->containsFSRCentroids())
    _track_generator->generateFSRCentroids(_FSR_volumes);

  /* Assign starting fluxes to the scalar fluxes array */
#pragma omp parallel for
  for (long r=0; r < _num_FSRs; r++) {

    /* Get the cell coordinates */
    Point* centroid = _geometry->getFSRCentroid(r);
    double* centroid_xyz = centroid->getXYZ();
    int cell_xyz[3];
    for (int i=0; i < 3; i++)
      cell_xyz[i] = centroid_xyz[i] / tolerance;

    /* Look at all cell combinations */
    double min_dist = std::numeric_limits<double>::max();
    long load_fsr = -1;
    for (int dx=-1; dx <= 1; dx++) {
      for (int dy=-1; dy <= 1; dy++) {
        for (int dz=-1; dz <= 1; dz++) {

          int new_cell_xyz[3];
          int d[3] = {dx, dy, dz};
          for (int i=0; i < 3; i++)
            new_cell_xyz[i] = cell_xyz[i] + d[i];

          /* Make sure index is within origin min/max bounds */
          for (int i=0; i < 3; i++) {
            if (new_cell_xyz[i] > max_ind[i])
              new_cell_xyz[i] = max_ind[i];
            if (cell_xyz[i] < min_ind[i])
              new_cell_xyz[i] = min_ind[i];
          }

          /* Calculate index */
          long index = (new_cell_xyz[2] - min_ind[2]) * nx * ny +
                       (new_cell_xyz[1] - min_ind[1]) * nx +
                       new_cell_xyz[0] - min_ind[0];

          /* Lookup all FSRs in the cell and check for distance to centroid */
          if (hashed_lookup.find(index) != hashed_lookup.end()) {
            for (int j =0; j < hashed_lookup[index].size(); j++) {
              long fsr_id = hashed_lookup[index].at(j);
              long dist = centroid->distance(x_coord[fsr_id], y_coord[fsr_id],
                                             z_coord[fsr_id]);
              if (dist < min_dist) {
                min_dist = dist;
                load_fsr = fsr_id;
              }
            }
          }
        }
      }
    }

    /* Check against maximum centroid mismatch */
    if (min_dist > max_centroid_error)
      max_centroid_error = min_dist;

    /* Check to ensure the loaded FSR is positive */
    if (load_fsr < 0)
      log_printf(ERROR, "Loaded FSR %d with location (%3.2f, %3.2f, %3.2f) "
                 "and cell (%d, %d, %d)", load_fsr, centroid_xyz[0],
                 centroid_xyz[1], centroid_xyz[2], cell_xyz[0], cell_xyz[1],
                 cell_xyz[2]);

    /* Assign scalar fluxes */
    for (int e=0; e < _num_groups; e++)
      _scalar_flux[r*num_groups+e] = fluxes[load_fsr * num_groups + e];
  }

  /* Delete auxilary data structures */
  delete [] x_coord;
  delete [] y_coord;
  delete [] z_coord;
  delete [] fluxes;
  delete [] cell_x;
  delete [] cell_y;
  delete [] cell_z;

  log_printf(NORMAL, "FSR fluxes successfully loaded with maximum centroid "
             "error %6.4e", max_centroid_error);
}


/**
 * @brief A function that prints a summary of the input parameters.
 */
void Solver::printInputParamsSummary() {

  /* Print track laydown parameters */
  log_printf(NORMAL, "Number of azimuthal angles = %d",
             _quad->getNumAzimAngles());
  log_printf(NORMAL, "Azimuthal ray spacing = %f",
             _track_generator->getDesiredAzimSpacing());
  log_printf(NORMAL, "Number of polar angles = %d",
             _quad->getNumPolarAngles());
  if (_SOLVE_3D) {
    TrackGenerator3D* track_generator_3D =
      static_cast<TrackGenerator3D*>(_track_generator);
    log_printf(NORMAL, "Z-spacing = %f",
               track_generator_3D->getDesiredZSpacing());
  }

  /* Print source type */
  log_printf(NORMAL, "Source type = %s", _source_type.c_str());

  /* Print MOC stabilization */
  if (_stabilize_transport) {

    std::string stabilization_str;

    if (_stabilization_type == DIAGONAL)
      stabilization_str = "DIAGONAL";
    else if (_stabilization_type == YAMAMOTO)
      stabilization_str = "TY";
    else if (_stabilization_type == GLOBAL)
      stabilization_str = "GLOBAL";

    log_printf(NORMAL, "MOC Damping = %s (%3.2f)", stabilization_str.c_str(),
               _stabilization_factor);
  }
  else {
    log_printf(NORMAL, "MOC transport undamped");
  }

  /* Print CMFD parameters */
  if (_cmfd != NULL)
    _cmfd->printInputParamsSummary();
  else
    log_printf(NORMAL, "CMFD acceleration: OFF");
}


/**
 * @brief Prints the memory report for the BG/Q architecture
 */
#ifdef BGQ
void Solver::printBGQMemory() {
  uint64_t shared, persist, heapavail, stackavail, stack, heap, guard, mmap;
  Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);

  log_printf(NORMAL, "Allocated heap: %.2f MB, avail. heap: %.2f MB",
             (double)heap/(1024*1024),(double)heapavail/(1024*1024));
  log_printf(NORMAL, "Allocated stack: %.2f MB, avail. stack: %.2f MB",
             (double)stack/(1024*1024), (double)stackavail/(1024*1024));
  log_printf(NORMAL, "Memory: shared: %.2f MB, persist: %.2f MB, guard: %.2f "
             "MB, mmap: %.2f MB\n", (double)shared/(1024*1024),
             (double)persist/(1024*1024), (double)guard/(1024*1024),
             (double)mmap/(1024*1024));
}
#endif
