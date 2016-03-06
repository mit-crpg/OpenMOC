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
  _num_parallel_track_groups = 0;

  _num_SRs = 0;
  _num_fissionable_SRs = 0;
  _SR_volumes = NULL;
  _SR_materials = NULL;

  _track_generator = NULL;
  _geometry = NULL;
  _cmfd = NULL;
  _exp_evaluator = new ExpEvaluator();

  _tracks = NULL;
  _polar_weights = NULL;
  _boundary_flux = NULL;

  _scalar_flux = NULL;
  _old_scalar_flux = NULL;
  _fixed_sources = NULL;
  _reduced_sources = NULL;

  if (track_generator != NULL)
    setTrackGenerator(track_generator);

  /* Default polar quadrature */
  _user_polar_quad = false;
  _polar_quad = new TYPolarQuad();
  _num_polar = 3;
  _polar_times_groups = 0;

  _num_iterations = 0;
  setConvergenceThreshold(1E-5);
  _user_fluxes = false;

  _timer = new Timer();
}


/**
 * @brief Destructor deletes arrays of boundary angular fluxes,
 *        scalar fluxes and sources for each SR and energy group.
 * @details Deallocates memory for all arrays allocated for the Solver,
 *          including fluxes, sources, quadrature weights, and exponential
 *          linear interpolation table.
 */
Solver::~Solver() {

  if (_SR_volumes != NULL)
    delete [] _SR_volumes;

  if (_SR_materials != NULL)
    delete [] _SR_materials;

  if (_polar_weights != NULL)
    delete [] _polar_weights;

  if (_boundary_flux != NULL)
    delete [] _boundary_flux;

  if (_scalar_flux != NULL && !_user_fluxes)
    delete [] _scalar_flux;

  if (_old_scalar_flux != NULL)
    delete [] _old_scalar_flux;

  if (_fixed_sources != NULL)
    delete [] _fixed_sources;

  if (_reduced_sources != NULL)
    delete [] _reduced_sources;

  if (_exp_evaluator != NULL)
    delete _exp_evaluator;

  if (_timer != NULL)
    delete _timer;

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
 * @brief Returns a pointer to the PolarQuad.
 * @return a pointer to the PolarQuad
 */
PolarQuad* Solver::getPolarQuad() {

  if (_polar_quad == NULL)
    log_printf(ERROR, "Unable to return the Solver's PolarQuad "
               "since it has not yet been set");

  return _polar_quad;
}


/**
 * @brief Returns the calculated volume for a flat source region.
 * @param sr_id the flat source region ID of interest
 * @return the flat source region volume
 */
FP_PRECISION Solver::getSRVolume(int sr_id) {

  if (sr_id < 0 || sr_id > _num_SRs)
    log_printf(ERROR, "Unable to get the volume for SR %d since the SR "
               "IDs lie in the range (0, %d)", sr_id, _num_SRs);

  else if (_SR_volumes == NULL)
    log_printf(ERROR, "Unable to get the volume for SR %d since the SR "
               "volumes have not yet been computed", sr_id);

  return _SR_volumes[sr_id];
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
FP_PRECISION Solver::getKeff() {
  return _k_eff;
}


/**
 * @brief Returns the threshold for source/flux convergence.
 * @return the threshold for source/flux convergence
 */
FP_PRECISION Solver::getConvergenceThreshold() {
  return _converge_thresh;
}


/**
 * @brief Get the maximum allowable optical length for a track segment
 * @return The max optical length
 */
FP_PRECISION Solver::getMaxOpticalLength() {
  return _exp_evaluator->getMaxOpticalLength();
}


/**
 * @brief Returns whether the solver is using double floating point precision.
 * @return true if using double precision float point arithmetic
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
 * @return true if using linear interpolation to compute exponentials
 */
bool Solver::isUsingExponentialInterpolation() {
  return _exp_evaluator->isUsingInterpolation();
}


/**
 * @brief Returns the source for some energy group for a flat source region
 * @details This is a helper routine used by the openmoc.process module.
 * @param sr_id the ID for the SR of interest
 * @param group the energy group of interest
 * @return the flat source region source
 */
FP_PRECISION Solver::getSRSource(int sr_id, int group) {

  if (sr_id >= _num_SRs)
    log_printf(ERROR, "Unable to return a source for SR ID = %d "
               "since the max SR ID = %d", sr_id, _num_SRs-1);

  else if (sr_id < 0)
    log_printf(ERROR, "Unable to return a source for SR ID = %d "
               "since SRs do not have negative IDs", sr_id);

  else if (group-1 >= _num_groups)
    log_printf(ERROR, "Unable to return a source in group %d "
               "since there are only %d groups", group, _num_groups);

  else if (group <= 0)
    log_printf(ERROR, "Unable to return a source in group %d "
               "since groups must be greater or equal to 1", group);

  else if (_scalar_flux == NULL)
    log_printf(ERROR, "Unable to return a source "
               "since it has not yet been computed");

  /* Get Material and cross-sections */
  Material* material = _SR_materials[sr_id];
  FP_PRECISION* sigma_s = material->getSigmaS();
  FP_PRECISION* fiss_mat = material->getFissionMatrix();

  FP_PRECISION fission_source = 0.0;
  FP_PRECISION scatter_source = 0.0;
  FP_PRECISION total_source;

  /* Compute total scattering and fission sources for this SR */
  for (int g=0; g < _num_groups; g++) {
    scatter_source += sigma_s[(group-1)*(_num_groups)+g]
                      * _scalar_flux(sr_id-1,g);
    fission_source += fiss_mat[(group-1)*(_num_groups)+g]
                      * _scalar_flux(sr_id-1,g);
  }

  fission_source /= _k_eff;

  /* Compute the total source */
  total_source = fission_source + scatter_source;

  /* Add in fixed source (if specified by user) */
  total_source += _fixed_sources(sr_id,group-1);

  /* Normalize to solid angle for isotropic approximation */
  total_source *= ONE_OVER_FOUR_PI;

  return total_source;
}


/**
 * @brief Returns the scalar flux for some SR and energy group.
 * @param sr_id the ID for the SR of interest
 * @param group the energy group of interest
 * @return the SR scalar flux
 */
FP_PRECISION Solver::getFlux(int sr_id, int group) {

  if (sr_id >= _num_SRs)
    log_printf(ERROR, "Unable to return a scalar flux for SR ID = %d "
               "since the max SR ID = %d", sr_id, _num_SRs-1);

  else if (sr_id < 0)
    log_printf(ERROR, "Unable to return a scalar flux for SR ID = %d "
               "since SRs do not have negative IDs", sr_id);

  else if (group-1 >= _num_groups)
    log_printf(ERROR, "Unable to return a scalar flux in group %d "
               "since there are only %d groups", group, _num_groups);

  else if (group <= 0)
    log_printf(ERROR, "Unable to return a scalar flux in group %d "
               "since groups must be greater or equal to 1", group);

  else if (_scalar_flux == NULL)
    log_printf(ERROR, "Unable to return a scalar flux "
             "since it has not yet been computed");

  return _scalar_flux(sr_id,group-1);
}


/**
 * @brief Sets the Geometry for the Solver.
 * @details This is a private setter method for the Solver and is not
 *          intended to be called by the user.
 * @param geometry a pointer to a Geometry object
 */
void Solver::setGeometry(Geometry* geometry) {

  if (geometry->getNumSRs() == 0)
    log_printf(ERROR, "Unable to set the Geometry for the Solver since the "
               "Geometry has not yet initialized SRs");

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
  _num_parallel_track_groups = _track_generator->getNumParallelTrackGroups();
  _tot_num_tracks = _track_generator->getNumTracks();
  _tracks = _track_generator->getTracksByParallelGroup();

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
  _polar_times_groups = _num_groups * _num_polar;
}


/**
 * @brief Sets the threshold for source/flux convergence.
 * @brief The default threshold for convergence is 1E-5.
 * @param source_thresh the threshold for source/flux convergence
 */
void Solver::setConvergenceThreshold(FP_PRECISION threshold) {

  if (threshold <= 0.0)
    log_printf(ERROR, "Unable to set the convergence threshold to %f "
               "since it is not a positive number", threshold);

  _converge_thresh = threshold;
}


/**
 * @brief Assign a fixed source for a flat source region and energy group.
 * @param sr_id the flat source region ID
 * @param group the energy group
 * @param source the volume-averaged source in this group
 */
void Solver::setFixedSourceBySR(int sr_id, int group, FP_PRECISION source) {
  _fix_src_SR_map[std::pair<int, int>(sr_id, group)] = source;
}


/**
 * @brief Assign a fixed source for a Cell and energy group.
 * @param cell the Cell of interest
 * @param group the energy group
 * @param source the volume-averaged source in this group
 */
void Solver::setFixedSourceByCell(Cell* cell, int group, FP_PRECISION source) {

  /* Recursively add the source to all Cells within a FILL type Cell */
  if (cell->getType() == FILL) {
    std::map<int, Cell*> cells = cell->getAllCells();
    std::map<int, Cell*>::iterator iter;
    for (iter = cells.begin(); iter != cells.end(); ++iter)
      setFixedSourceByCell(iter->second, group, source);
  }
  else
    _fix_src_cell_map[std::pair<Cell*, int>(cell, group)] = source;
}


/**
 * @brief Assign a fixed source for a Material and energy group.
 * @param material the Material of interest
 * @param group the energy group
 * @param source the volume-averaged source in this group
 */
void Solver::setFixedSourceByMaterial(Material* material, int group,
                                      FP_PRECISION source) {
  _fix_src_material_map[std::pair<Material*, int>(material, group)] = source;
}


/**
 * @brief Set the maximum allowable optical length for a track segment
 * @param max_optical_length The max optical length
 */
void Solver::setMaxOpticalLength(FP_PRECISION max_optical_length) {
  _exp_evaluator->setMaxOpticalLength(max_optical_length);
}


/**
 * @brief Set the precision, or maximum allowable approximation error, of the
 *        the exponential interpolation table.
 * @details By default, the precision is 1E-5 based on the analysis in
 *          Yamamoto's 2003 paper.
 * @param precision the precision of the exponential interpolation table,
 */
void Solver::setExpPrecision(FP_PRECISION precision) {
  _exp_evaluator->setExpPrecision(precision);
}


/**
 * @brief Informs the Solver to use linear interpolation to compute the
 *        exponential in the transport equation.
 */
void Solver::useExponentialInterpolation() {
  _exp_evaluator->useInterpolation();
}


/**
 * @brief Informs the Solver to use the exponential intrinsic exp(...)
 *        function to compute the exponential in the transport equation.
 */
void Solver::useExponentialIntrinsic() {
  _exp_evaluator->useIntrinsic();
}


/**
 * @brief Initializes a new PolarQuad object.
 * @details Deletes memory old PolarQuad if one was previously allocated.
 */
void Solver::initializePolarQuadrature() {

  FP_PRECISION* azim_weights = _track_generator->getAzimWeights();

  /* Initialize the PolarQuad object */
  _polar_quad->setNumPolarAngles(_num_polar);
  _polar_quad->initialize();
  _polar_times_groups = _num_groups * _num_polar;

  /* Deallocate polar weights if previously assigned */
  if (_polar_weights != NULL)
    delete [] _polar_weights;

  _polar_weights = new FP_PRECISION[_num_azim*_num_polar];

  /* Compute the total azimuthal weight for tracks at each polar angle */
#pragma omp parallel for schedule(guided)
  for (int i=0; i < _num_azim; i++) {
    for (int p=0; p < _num_polar; p++)
      _polar_weights(i,p) =
           azim_weights[i] * _polar_quad->getMultiple(p) * FOUR_PI;
  }
}


/**
 * @brief Initializes new ExpEvaluator object to compute exponentials.
 */
void Solver::initializeExpEvaluator() {

  _exp_evaluator->setPolarQuadrature(_polar_quad);

  if (_exp_evaluator->isUsingInterpolation()) {

    /* Find minimum of optional user-specified and actual max taus */
    FP_PRECISION max_tau_a = _track_generator->getMaxOpticalLength();
    FP_PRECISION max_tau_b = _exp_evaluator->getMaxOpticalLength();
    FP_PRECISION max_tau = std::min(max_tau_a, max_tau_b) + TAU_NUDGE;

    /* Split Track segments so that none has a greater optical length */
    _track_generator->splitSegments(max_tau);

    /* Initialize exponential interpolation table */
    _exp_evaluator->setMaxOpticalLength(max_tau);
    _exp_evaluator->initialize();
  }
}


/**
 * @brief Initializes the Material fission matrices.
 * @details In an adjoint calculation, this routine will transpose the
 *          scattering and fission matrices in each material.
 * @param mode the solution type (FORWARD or ADJOINT)
 */
void Solver::initializeMaterials(solverMode mode) {

  log_printf(INFO, "Initializing materials...");

  std::map<int, Material*> materials = _geometry->getAllMaterials();
  std::map<int, Material*>::iterator m_iter;

  for (m_iter = materials.begin(); m_iter != materials.end(); ++m_iter) {
    m_iter->second->buildFissionMatrix();

    if (mode == ADJOINT)
      m_iter->second->transposeProductionMatrices();
  }
}


/**
 * @brief Initializes the SR volumes and Materials array.
 * @details This method assigns each SR a unique, monotonically increasing
 *          ID, sets the Material for each SR, and assigns a volume based on
 *          the cumulative length of all of the segments inside the SR.
 */
void Solver::initializeSRs() {

  log_printf(INFO, "Initializing flat source regions...");

  /* Delete old SR arrays if they exist */
  if (_SR_volumes != NULL)
    delete [] _SR_volumes;

  if (_SR_materials != NULL)
    delete [] _SR_materials;

  /* Retrieve simulation parameters from the Geometry */
  _num_SRs = _geometry->getNumSRs();
  _num_groups = _geometry->getNumEnergyGroups();
  _polar_times_groups = _num_groups * _num_polar;
  _num_materials = _geometry->getNumMaterials();

  /* Get an array of volumes indexed by SR  */
  _SR_volumes = _track_generator->getSRVolumes();

  /* Generate the SR centroids */
  _track_generator->generateSRCentroids();

  /* Attach the correct materials to each track segment */
  _track_generator->initializeSegments();

  /* Allocate an array of Material pointers indexed by SR */
  _SR_materials = new Material*[_num_SRs];

  /* Loop over all SRs to extract SR material pointers */
  for (int r=0; r < _num_SRs; r++) {
    _SR_materials[r] = _geometry->findSRMaterial(r);
    log_printf(INFO, "SR ID = %d has Material ID = %d and volume = %f ",
               r, _SR_materials[r]->getId(), _SR_volumes[r]);
  }
}


/**
 * @brief Counts the number of fissionable flat source regions.
 * @details This routine is used by the Solver::computeEigenvalue(...)
 *          routine which uses the number of fissionable SRs to normalize
 *          the residual on the fission source distribution.
 */
void Solver::countFissionableSRs() {

  log_printf(INFO, "Counting fissionable SRs...");

  /* Count the number of fissionable SRs */
  _num_fissionable_SRs = 0;
  for (int r=0; r < _num_SRs; r++) {
    if (_SR_materials[r]->isFissionable())
      _num_fissionable_SRs++;
  }
}


/**
 * @brief Assigns fixed sources assigned by Cell, Material to SRs.
 */
void Solver::initializeFixedSources() {

  Cell* sr_cell;
  Material* sr_material;
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

    /* Search for this Cell in all SRs */
    for (int r=0; r < _num_SRs; r++) {
      sr_cell = _geometry->findCellContainingSR(r);
      if (cell_group_key.first->getId() == sr_cell->getId())
        setFixedSourceBySR(r, group, source);
    }
  }

  /** Fixed sources assigned by Material */
  for (mat_iter = _fix_src_material_map.begin();
       mat_iter != _fix_src_material_map.end(); ++mat_iter) {

    /* Get the Material with an assigned fixed source */
    mat_group_key = mat_iter->first;
    group = mat_group_key.second;
    source = _fix_src_material_map[mat_group_key];

    for (int r=0; r < _num_SRs; r++) {
      sr_material = _geometry->findSRMaterial(r);
      if (mat_group_key.first->getId() == sr_material->getId())
        setFixedSourceBySR(r, group, source);
    }
  }
}


/**
 * @brief Initializes a Cmfd object for acceleratiion prior to source iteration.
 * @details Instantiates a dummy Cmfd object if one was not assigned to
 *          the Solver by the user and initializes SRs, materials, fluxes
 *          and the Mesh object. This method is for internal use only
 *          and should not be called directly by the user.
 */
void Solver::initializeCmfd() {

  log_printf(INFO, "Initializing CMFD...");

  /* Retrieve CMFD from the Geometry */
  _cmfd = _geometry->getCmfd();

  /* If the user did not initialize Cmfd, simply return */
  if (_cmfd == NULL)
    return;
  else if (!_cmfd->isFluxUpdateOn())
    return;

  /* Intialize the CMFD energy group structure */
  _cmfd->setSourceConvergenceThreshold(_converge_thresh*1.e-1);
  _cmfd->setNumMOCGroups(_num_groups);
  _cmfd->initializeGroupMap();

  /* Give CMFD number of SRs and SR property arrays */
  _cmfd->setNumSRs(_num_SRs);
  _cmfd->setSRVolumes(_SR_volumes);
  _cmfd->setSRMaterials(_SR_materials);
  _cmfd->setSRFluxes(_scalar_flux);
  _cmfd->setPolarQuadrature(_polar_quad);
  _cmfd->setGeometry(_geometry);
  _cmfd->initialize();
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
  computeSRFissionSources();
  transportSweep();
  addSourceToScalarFlux();
}


/**
 * @brief This method performs one transport sweep using the scatter source.
 * @details This is a helper routine used for Krylov subspace methods.
 */
void Solver::scatterTransportSweep() {
  computeSRScatterSources();
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
 *          sweeps with a 1E-5 threshold on the average SR scalar flux. These
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
 * @param mode the solution type (FORWARD or ADJOINT)
 * @param only_fixed_source use only fixed sources (true by default)
 */
void Solver::computeFlux(int max_iters, solverMode mode,
                         bool only_fixed_source) {

  if (_track_generator == NULL)
    log_printf(ERROR, "The Solver is unable to compute the flux "
               "since it does not contain a TrackGenerator");

  log_printf(NORMAL, "Computing the flux...");

  /* Clear all timing data from a previous simulation run */
  clearTimerSplits();

  /* Start the timer to record the total time to converge the flux */
  _timer->startTimer();

  /* Initialize keff to 1 for SR source calculations */
  _k_eff = 1.;

  _num_iterations = 0;
  FP_PRECISION residual = 0.;

  /* Initialize data structures */
  initializeSRs();
  initializeMaterials(mode);
  countFissionableSRs();
  initializePolarQuadrature();
  initializeExpEvaluator();

  /* Initialize new flux arrays if a) the user requested the use of
   * only fixed sources or b) no previous simulation was performed which
   * initialized and computed the flux (e.g., an eigenvalue calculation) */
  if (only_fixed_source || _num_iterations == 0) {
    initializeFluxArrays();
    flattenSRFluxes(0.0);
    storeSRFluxes();
  }

  initializeSourceArrays();
  zeroTrackFluxes();

  /* Compute the sum of fixed, total and scattering sources */
  computeSRSources();

  /* Source iteration loop */
  for (int i=0; i < max_iters; i++) {

    transportSweep();
    addSourceToScalarFlux();
    residual = computeResidual(SCALAR_FLUX);
    storeSRFluxes();
    _num_iterations++;

    log_printf(NORMAL, "Iteration %d:\tres = %1.3E", i, residual);

    /* Check for convergence */
    if (i > 1 && residual < _converge_thresh)
      break;
  }

  if (_num_iterations == max_iters-1)
    log_printf(WARNING, "Unable to converge the flux");

  resetMaterials(mode);

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
 *          sweeps with a 1E-5 threshold on the integrated SR total source.
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
 *          // Find the flux distribution resulting from the fixed sources
 *          solver.computeSource(max_iters=100, k_eff=0.981)
 * @endcode
 *
 * @param max_iters the maximum number of source iterations to allow
 * @param mode the solution type (FORWARD or ADJOINT)
 * @param k_eff the sub/super-critical eigenvalue (default 1.0)
 * @param res_type the type of residual used for the convergence criterion
 */
void Solver::computeSource(int max_iters, solverMode mode,
                           double k_eff, residualType res_type) {

  if (_track_generator == NULL)
    log_printf(ERROR, "The Solver is unable to compute the source "
               "since it does not contain a TrackGenerator");

  else if (k_eff <= 0.)
    log_printf(ERROR, "The Solver is unable to compute the source with "
               "keff = %f since it is not a positive value", k_eff);

  log_printf(NORMAL, "Computing the source...");

  /* Clear all timing data from a previous simulation run */
  clearTimerSplits();

  /* Start the timer to record the total time to converge the flux */
  _timer->startTimer();

  /* Set the eigenvalue to the user-specified value */
  _k_eff = k_eff;

  _num_iterations = 0;
  FP_PRECISION residual = 0.;

  /* Initialize data structures */
  initializeSRs();
  initializeMaterials(mode);
  initializePolarQuadrature();
  initializeExpEvaluator();
  initializeFluxArrays();
  initializeSourceArrays();

  /* Guess unity scalar flux for each region */
  flattenSRFluxes(1.0);
  storeSRFluxes();
  zeroTrackFluxes();

  /* Source iteration loop */
  for (int i=0; i < max_iters; i++) {

    computeSRSources();
    transportSweep();
    addSourceToScalarFlux();
    residual = computeResidual(res_type);
    storeSRFluxes();
    _num_iterations++;

    log_printf(NORMAL, "Iteration %d:\tres = %1.3E", i, residual);

    /* Check for convergence */
    if (i > 1 && residual < _converge_thresh)
      break;
  }

  if (_num_iterations == max_iters-1)
    log_printf(WARNING, "Unable to converge the source distribution");

  resetMaterials(mode);

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
 *          sweeps with a 1E-5 threshold on the integrated SR fission source.
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
 * @param mode the solution type (FORWARD or ADJOINT)
 * @param res_type the type of residual used for the convergence criterion
 */
void Solver::computeEigenvalue(int max_iters, solverMode mode,
                               residualType res_type) {

  if (_track_generator == NULL)
    log_printf(ERROR, "The Solver is unable to compute the eigenvalue "
               "since it does not contain a TrackGenerator");

  log_printf(NORMAL, "Computing the eigenvalue...");

  /* Clear all timing data from a previous simulation run */
  clearTimerSplits();

  /* Start the timer to record the total time to converge the source */
  _timer->startTimer();

  _num_iterations = 0;
  FP_PRECISION residual = 0.;

  /* An initial guess for the eigenvalue */
  _k_eff = 1.0;

  /* Initialize data structures */
  initializeSRs();
  initializeMaterials(mode);
  countFissionableSRs();
  initializePolarQuadrature();
  initializeExpEvaluator();
  initializeFluxArrays();
  initializeSourceArrays();
  initializeCmfd();

  /* Set scalar flux to unity for each region */
  flattenSRFluxes(1.0);
  storeSRFluxes();
  zeroTrackFluxes();

  /* Source iteration loop */
  for (int i=0; i < max_iters; i++) {
    normalizeFluxes();
    computeSRSources();
    transportSweep();
    addSourceToScalarFlux();

    /* Solve CMFD diffusion problem and update MOC flux */
    if (_cmfd != NULL && _cmfd->isFluxUpdateOn()) {
      _k_eff = _cmfd->computeKeff(i);
      _cmfd->updateBoundaryFlux(_tracks, _boundary_flux, _tot_num_tracks);
    }
    else
      computeKeff();

    log_printf(NORMAL, "Iteration %d:\tk_eff = %1.6f"
               "\tres = %1.3E", i, _k_eff, residual);

    residual = computeResidual(res_type);
    storeSRFluxes();
    _num_iterations++;

    /* Check for convergence */
    if (i > 1 && residual < _converge_thresh)
      break;
  }

  if (_num_iterations == max_iters-1)
    log_printf(WARNING, "Unable to converge the source distribution");

  resetMaterials(mode);

  _timer->stopTimer();
  _timer->recordSplit("Total time");
}


/**
 * @brief Deletes the Timer's timing entries for each timed code section
 *        code in the source convergence loop.
 */
void Solver::clearTimerSplits() {
  _timer->clearSplit("Total time");
}


/**
 * @brief Prints a report of the timing statistics to the console.
 */
void Solver::printTimerReport() {

  std::string msg_string;

  log_printf(TITLE, "TIMING REPORT");

  /* Get the total runtime */
  double tot_time = _timer->getSplit("Total time");
  msg_string = "Total time to solution";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), tot_time);

  /* Time per iteration */
  double time_per_iter = tot_time / _num_iterations;
  msg_string = "Solution time per iteration";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), time_per_iter);

  /* Time per segment */
  int num_segments = _track_generator->getNumSegments();
  int num_integrations = 2 * _num_polar * _num_groups * num_segments;
  double time_per_integration = (time_per_iter / num_integrations);
  msg_string = "Time per segment integration";
  msg_string.resize(53, '.');
  log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), time_per_integration);

  set_separator_character('-');
  log_printf(SEPARATOR, "-");

  msg_string = "           # tracks          # segments          # SRs";
  log_printf(RESULT, "%s", msg_string.c_str());
  log_printf(SEPARATOR, "-");

  int num_digits = (int) log10((double) _tot_num_tracks);
  num_digits += (int) log10((double) num_segments);
  num_digits += (int) log10((double) _num_SRs);

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
      msg << _num_SRs;
  }

  log_printf(RESULT, "%s", msg.str().c_str());
  log_printf(SEPARATOR, "-");
}
