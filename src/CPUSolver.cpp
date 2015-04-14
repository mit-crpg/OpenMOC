#include "CPUSolver.h"


/**
 * @brief Constructor initializes array pointers for Tracks and Materials.
 * @details The constructor retrieves the number of energy groups and FSRs
 *          and azimuthal angles from the Geometry and TrackGenerator if
 *          passed in as parameters by the user. The constructor initalizes
 *          the number of OpenMP threads to a default of 1.
 * @param track_generator an optional pointer to the TrackGenerator
 */
CPUSolver::CPUSolver(TrackGenerator* track_generator)
  : Solver(track_generator) {

  setNumThreads(1);
  _FSR_locks = NULL;
}


/**
 * @brief Destructor deletes array for OpenMP mutual exclusion locks for
 *        FSR scalar flux updates, and calls Solver parent class destructor
 *        to deletes arrays for fluxes and sources.
 */
CPUSolver::~CPUSolver() {
  if (_FSR_locks != NULL)
    delete [] _FSR_locks;
}


/**
 * @brief Returns the number of shared memory OpenMP threads in use.
 * @return the number of threads
 */
int CPUSolver::getNumThreads() {
  return _num_threads;
}


/**
 * @brief Returns the scalar flux for some FSR and energy group.
 * @param fsr_id the ID for the FSR of interest
 * @param group the energy group of interest
 * @return the FSR scalar flux
 */
FP_PRECISION CPUSolver::getFSRScalarFlux(int fsr_id, int group) {

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
    log_printf(ERROR, "Unable to return a scalar flux"
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
FP_PRECISION CPUSolver::getFSRSource(int fsr_id, int group) {

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
  FP_PRECISION source = 0.;

  /* Compute fission source */
  if (material->isFissionable()) {
    for (int e=0; e < _num_groups; e++)
      source += _scalar_flux(fsr_id,e) * nu_sigma_f[e];
    source /= _k_eff * chi[group-1];
  }

  /* Compute scatter source */
  for (int g=0; g < _num_groups; g++)
    source += material->getSigmaSByGroupInline(g,group-1)
              * _scalar_flux(fsr_id,g);

  /* Compute fixed source (if specified by user) */
  source += _fixed_sources(fsr_id,group-1);

  /* Normalize to solid angle for isotropic approximation */
  source *= ONE_OVER_FOUR_PI;

  return source;
}


/**
 * @brief Return a scalar flux array indexed by FSR IDs and energy groups.
 * @details The energy groups are the innermost index and the FSR ID is
 *          the outermost index.
 * @return an array of flat source region scalar fluxes
 */
FP_PRECISION* CPUSolver::getFSRScalarFluxes() {

  if (_scalar_flux == NULL)
    log_printf(ERROR, "Unable to return FSR scalar flux array "
               "since it has not yet been allocated in memory");

  return _scalar_flux;
}


/**
 * @brief Sets the number of shared memory OpenMP threads to use (>0).
 * @param num_threads the number of threads
 */
void CPUSolver::setNumThreads(int num_threads) {

  if (num_threads <= 0)
    log_printf(ERROR, "Unable to set the number of threads to %d "
               "since it is less than or equal to 0", num_threads);

  /* Set the number of threads for OpenMP */
  _num_threads = num_threads;
  omp_set_num_threads(_num_threads);
}


/**
 * @brief Assign a fixed source for a flat source region and energy group.
 * @details Fixed sources should be scaled to reflect the fact that OpenMOC 
 *          normalizes the scalar flux such that the total energy- and 
 *          volume-integrated production rate sums to 1.0.
 * @param fsr_id the flat source region ID
 * @param group the energy group
 * @param source the volume-averaged source in this group
 */
void CPUSolver::setFixedSourceByFSR(int fsr_id, int group, 
                                    FP_PRECISION source) {

  Solver::setFixedSourceByFSR(fsr_id, group, source);

  /* Allocate the fixed sources array if not yet allocated */
  if (_fixed_sources == NULL) {
    int size = _num_FSRs * _num_groups;
    _fixed_sources = new FP_PRECISION[size];
    memset(_fixed_sources, 0.0, sizeof(FP_PRECISION) * size);
  }

  /* Store the fixed source for this FSR and energy group */
  _fixed_sources(fsr_id,group-1) = source;  
}


/**
 * @brief Initializes the FSR volumes and Materials array.
 * @details This method allocates and initializes an array of OpenMP
 *          mutual exclusion locks for each FSR for use in the
 *          transport sweep algorithm.
 */
void CPUSolver::initializeFSRs() {

  Solver::initializeFSRs();

  /* Allocate array of mutex locks for each FSR */
  _FSR_locks = new omp_lock_t[_num_FSRs];

  /* Loop over all FSRs to initialize OpenMP locks */
  #pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_FSRs; r++)
    omp_init_lock(&_FSR_locks[r]);
}


/**
 * @brief Allocates memory for Track boundary angular flux and leakage
 *        and FSR scalar flux arrays.
 * @details Deletes memory for old flux arrays if they were allocated
 *          for a previous simulation.
 */
void CPUSolver::initializeFluxArrays() {

  /* Delete old flux arrays if they exist */
  if (_boundary_flux != NULL)
    delete [] _boundary_flux;

  if (_boundary_leakage != NULL)
    delete [] _boundary_leakage;

  if (_scalar_flux != NULL)
    delete [] _scalar_flux;

  /* Allocate memory for the Track boundary flux and leakage arrays */
  try{
    int size = 2 * _tot_num_tracks * _polar_times_groups;
    _boundary_flux = new FP_PRECISION[size];
    _boundary_leakage = new FP_PRECISION[size];

    /* Allocate an array for the FSR scalar flux */
    size = _num_FSRs * _num_groups;
    _scalar_flux = new FP_PRECISION[size];
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the fluxes");
  }
}


/**
 * @brief Allocates memory for FSR source arrays.
 * @details Deletes memory for old source arrays if they were allocated for a
 *          previous simulation.
 */
void CPUSolver::initializeSourceArrays() {

  /* Delete old sources arrays if they exist */
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

  /* Allocate memory for all source arrays */
  try{
    int size = _num_FSRs * _num_groups;
    _fission_sources = new FP_PRECISION[size];
    _reduced_sources = new FP_PRECISION[size];

    /* If no fixed sources were assigned, use a zeroes array */
    if (_fixed_sources == NULL) {
      _fixed_sources = new FP_PRECISION[size];
      memset(_fixed_sources, 0.0, sizeof(FP_PRECISION) * size);
    }

    size = _num_threads * _num_groups;
    _scatter_sources = new FP_PRECISION[size];

    size = _num_FSRs;
    _old_fission_sources = new FP_PRECISION[size];
    _source_residuals = new FP_PRECISION[size];

  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for FSR sources");
  }
}


/**
 * @brief Zero each Track's boundary fluxes for each energy group
 *        and polar angle in the "forward" and "reverse" directions.
 */
void CPUSolver::zeroTrackFluxes() {

  #pragma omp parallel for schedule(guided)
  for (int t=0; t < _tot_num_tracks; t++) {
    for (int d=0; d < 2; d++) {
      for (int p=0; p < _num_polar; p++) {
        for (int e=0; e < _num_groups; e++) {
          _boundary_flux(t,d,p,e) = 0.0;
        }
      }
    }
  }
}


/**
 * @brief Set the scalar flux for each FSR and energy group to some value.
 * @param value the value to assign to each FSR scalar flux
 */
void CPUSolver::flattenFSRFluxes(FP_PRECISION value) {

  #pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _num_groups; e++)
      _scalar_flux(r,e) = value;
  }
}


/**
 * @brief Set the source for each FSR and energy group to some value.
 * @param value the value to assign to each FSR source
 */
void CPUSolver::flattenFSRSources(FP_PRECISION value) {

  #pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_FSRs; r++) 
    _old_fission_sources[r] = value;
}


/**
 * @brief Normalizes all FSR scalar fluxes and Track boundary angular
 *        fluxes to the total fission source (times \f$ \nu \f$).
 */
void CPUSolver::normalizeFluxes() {

  FP_PRECISION* nu_sigma_f;
  FP_PRECISION volume;
  FP_PRECISION tot_fission_source;
  FP_PRECISION norm_factor;

  /* Compute total fission source for each FSR, energy group */
  #pragma omp parallel for private(volume, nu_sigma_f) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {

    /* Get pointers to important data structures */
    nu_sigma_f = _FSR_materials[r]->getNuSigmaF();
    volume = _FSR_volumes[r];

    for (int e=0; e < _num_groups; e++)
      _fission_sources(r,e) = nu_sigma_f[e] * _scalar_flux(r,e) * volume;
  }

  /* Compute the total fission source */
  tot_fission_source = pairwise_sum<FP_PRECISION>(_fission_sources,
                                                  _num_FSRs*_num_groups);

  /* Normalize scalar fluxes in each FSR */
  norm_factor = 1.0 / tot_fission_source;

  log_printf(DEBUG, "Tot. Fiss. Src. = %f, Norm. factor = %f",
             tot_fission_source, norm_factor);

  #pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _num_groups; e++)
      _scalar_flux(r,e) *= norm_factor;
  }

  /* Normalize angular boundary fluxes for each Track */
  #pragma omp parallel for schedule(guided)
  for (int t=0; t < _tot_num_tracks; t++) {
    for (int d=0; d < 2; d++) {
      for (int p=0; p < _num_polar; p++) {
        for (int e=0; e < _num_groups; e++) {
          _boundary_flux(t,d,p,e) *= norm_factor;
        }
      }
    }
  }
}


/**
 * @brief Computes the total source (fission and scattering) in each FSR.
 * @details This method computes the total source in each FSR based on
 *          this iteration's current approximation to the scalar flux. A
 *          residual for the source with respect to the source from the
 *          previous iteration is computed and returned. The residual
 *          is determined as follows:
 *          \f$ res = \sqrt{\frac{\displaystyle\sum \displaystyle\sum
 *                    \left(\frac{Q^i - Q^{i-1}}{Q^i}\right)^2}
 *                    {\# FSRs \times # groups}} \f$
 *
 * @return the residual between this source and the previous source
 */
FP_PRECISION CPUSolver::computeFSRSources() {

  int tid;
  FP_PRECISION scatter_source;
  FP_PRECISION fission_source;
  FP_PRECISION fsr_fission_source;
  FP_PRECISION* nu_sigma_f;
  FP_PRECISION* sigma_s;
  FP_PRECISION* sigma_t;
  FP_PRECISION* chi;
  Material* material;

  FP_PRECISION source_residual = 0.0;
  FP_PRECISION inverse_k_eff = 1.0 / _k_eff;

  /* For all FSRs, find the source */
  #pragma omp parallel for private(tid, material, nu_sigma_f, chi, \
    sigma_s, sigma_t, fission_source, scatter_source, fsr_fission_source) \
    schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {

    tid = omp_get_thread_num();
    material = _FSR_materials[r];
    nu_sigma_f = material->getNuSigmaF();
    chi = material->getChi();
    sigma_s = material->getSigmaS();
    sigma_t = material->getSigmaT();

    /* Initialize the source residual to zero */
    _source_residuals[r] = 0.;

    /* Initialize the fission sources to zero */
    fission_source = 0.0;
    fsr_fission_source = 0.0;

    /* Compute fission source for each group */
    if (material->isFissionable()) {
      for (int e=0; e < _num_groups; e++)
        _fission_sources(r,e) = _scalar_flux(r,e) * nu_sigma_f[e];

      fission_source = pairwise_sum<FP_PRECISION>(&_fission_sources(r,0),
                                                  _num_groups);
      fission_source *= inverse_k_eff;
    }

    /* Compute total scattering source for group G */
    for (int G=0; G < _num_groups; G++) {
      for (int g=0; g < _num_groups; g++)
        _scatter_sources(tid,g) = material->getSigmaSByGroupInline(g,G)
                                  * _scalar_flux(r,g);
      scatter_source=pairwise_sum<FP_PRECISION>(&_scatter_sources(tid,0),
                                                _num_groups);

      fsr_fission_source += fission_source * chi[G];
      _reduced_sources(r,G) = fission_source * chi[G];
      _reduced_sources(r,G) += scatter_source + _fixed_sources(r,G);
      _reduced_sources(r,G) *= ONE_OVER_FOUR_PI / sigma_t[G];
    }

    /* Compute the norm of residual of the source in the FSR */
    if (fsr_fission_source > 0.0)
      _source_residuals[r] = 
           pow((fsr_fission_source - _old_fission_sources[r]) / 
                fsr_fission_source, 2);

    /* Update the old source */
    _old_fission_sources[r] = fsr_fission_source;
  }

  /* Sum up the residuals from each FSR */
  source_residual = pairwise_sum<FP_PRECISION>(_source_residuals, 
                                               _num_FSRs);
  source_residual = sqrt(source_residual \
                         / (_num_fissionable_FSRs * _num_groups));

  return source_residual;
}


/**
 * @brief Compute \f$ k_{eff} \f$ from the total, fission and scattering
 *        reaction rates and leakage.
 * @details This method computes the current approximation to the
 *          multiplication factor on this iteration as follows:
 *          \f$ k_{eff} = \frac{\displaystyle\sum_{i \in I}
 *                        \displaystyle\sum_{g \in G} \nu \Sigma^F_g \Phi V_{i}}
 *                        {\displaystyle\sum_{i \in I}
 *                        \displaystyle\sum_{g \in G} (\Sigma^T_g \Phi V_{i} -
 *                        \Sigma^S_g \Phi V_{i} - L_{i,g})} \f$
 */
void CPUSolver::computeKeff() {

  int tid;
  Material* material;
  FP_PRECISION* sigma;
  FP_PRECISION volume;

  FP_PRECISION total = 0., fission = 0., scatter = 0., leakage = 0.;
  FP_PRECISION* FSR_rates = new FP_PRECISION[_num_FSRs];
  FP_PRECISION* group_rates = new FP_PRECISION[_num_threads * _num_groups];

  /* Loop over all FSRs and compute the volume-integrated total rates */
  #pragma omp parallel for private(tid, volume, \
    material, sigma) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {

    tid = omp_get_thread_num() * _num_groups;
    volume = _FSR_volumes[r];
    material = _FSR_materials[r];
    sigma = material->getSigmaT();

    for (int e=0; e < _num_groups; e++)
      group_rates[tid+e] = sigma[e] * _scalar_flux(r,e);

    FSR_rates[r]=pairwise_sum<FP_PRECISION>(&group_rates[tid], _num_groups);
    FSR_rates[r] *= volume;
  }

  /* Reduce total rates across FSRs */
  total = pairwise_sum<FP_PRECISION>(FSR_rates, _num_FSRs);

  /* Loop over all FSRs and compute the volume-integrated fission rates */
  #pragma omp parallel for private(tid, volume, \
    material, sigma) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {

    tid = omp_get_thread_num() * _num_groups;
    volume = _FSR_volumes[r];
    material = _FSR_materials[r];
    sigma = material->getNuSigmaF();

    for (int e=0; e < _num_groups; e++)
      group_rates[tid+e] = sigma[e] * _scalar_flux(r,e);

    FSR_rates[r]=pairwise_sum<FP_PRECISION>(&group_rates[tid], _num_groups);
    FSR_rates[r] *= volume;
  }

  /* Reduce fission rates across FSRs */
  fission = pairwise_sum<FP_PRECISION>(FSR_rates, _num_FSRs);

  /* Loop over all FSRs and compute the volume-integrated scattering rates */
  #pragma omp parallel for private(tid, volume, \
    material) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {

    tid = omp_get_thread_num() * _num_groups;
    volume = _FSR_volumes[r];
    material = _FSR_materials[r];

    FSR_rates[r] = 0.;

    for (int G=0; G < _num_groups; G++) {
      for (int g=0; g < _num_groups; g++)
        group_rates[tid+g] = material->getSigmaSByGroupInline(g,G)
                             * _scalar_flux(r,g);

      FSR_rates[r]+=pairwise_sum<FP_PRECISION>(&group_rates[tid], _num_groups);
    }

    FSR_rates[r] *= volume;
  }

  /* Reduce scattering rates across FSRs */
  scatter = pairwise_sum<FP_PRECISION>(FSR_rates, _num_FSRs);

  /* Reduce leakage array across Tracks, energy groups, polar angles */
  int size = 2 * _tot_num_tracks * _polar_times_groups;
  leakage = pairwise_sum<FP_PRECISION>(_boundary_leakage, size) * 0.5;

  _k_eff = fission / (total - scatter + leakage);

  log_printf(DEBUG, "tot = %f, fiss = %f, scatt = %f, leak = %f,"
             "k_eff = %f", total, fission, scatter, leakage, _k_eff);

  delete [] FSR_rates;
  delete [] group_rates;
}


/**
 * @brief This method performs one transport sweep of all azimuthal angles,
 *        Tracks, Track segments, polar angles and energy groups.
 * @details The method integrates the flux along each Track and updates the
 *          boundary fluxes for the corresponding output Track, while updating
 *          the scalar flux in each flat source region.
 */
void CPUSolver::transportSweep() {

  int tid;
  int min_track, max_track;
  int azim_index, num_segments;
  Track* curr_track;
  segment* curr_segment;
  segment* segments;
  FP_PRECISION* track_flux;

  log_printf(DEBUG, "Transport sweep with %d OpenMP threads", _num_threads);

  /* Initialize flux in each FSr to zero */
  flattenFSRFluxes(0.0);

  if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
    _cmfd->zeroSurfaceCurrents();

  /* Loop over azimuthal angle halfspaces */
  for (int i=0; i < 2; i++) {

    /* Compute the minimum and maximum Track IDs corresponding to
     * this azimuthal angular halfspace */
    min_track = i * (_tot_num_tracks / 2);
    max_track = (i + 1) * (_tot_num_tracks / 2);

    /* Loop over each thread within this azimuthal angle halfspace */
    #pragma omp parallel for private(curr_track, azim_index, num_segments, \
      curr_segment, segments, track_flux, tid) schedule(guided)
    for (int track_id=min_track; track_id < max_track; track_id++) {

      tid = omp_get_thread_num();

      /* Use local array accumulator to prevent false sharing*/
      FP_PRECISION* thread_fsr_flux;
      thread_fsr_flux = new FP_PRECISION[_num_groups];

      /* Initialize local pointers to important data structures */
      curr_track = _tracks[track_id];
      azim_index = curr_track->getAzimAngleIndex();
      num_segments = curr_track->getNumSegments();
      segments = curr_track->getSegments();
      track_flux = &_boundary_flux(track_id,0,0,0);

      /* Loop over each Track segment in forward direction */
      for (int s=0; s < num_segments; s++) {
        curr_segment = &segments[s];
        tallyScalarFlux(curr_segment, azim_index, track_flux, thread_fsr_flux);
        tallySurfaceCurrent(curr_segment, azim_index, track_flux, true);
      }

      /* Transfer boundary angular flux to outgoing Track */
      transferBoundaryFlux(track_id, azim_index, true, track_flux);

      /* Loop over each Track segment in reverse direction */
      track_flux += _polar_times_groups;

      for (int s=num_segments-1; s > -1; s--) {
        curr_segment = &segments[s];
        tallyScalarFlux(curr_segment, azim_index, track_flux, thread_fsr_flux);
        tallySurfaceCurrent(curr_segment, azim_index, track_flux, false);
      }
      delete thread_fsr_flux;

      /* Transfer boundary angular flux to outgoing Track */
      transferBoundaryFlux(track_id, azim_index, false, track_flux);
    }
  }

  return;
}


/**
 * @brief Computes the contribution to the FSR scalar flux from a Track segment.
 * @details This method integrates the angular flux for a Track segment across
 *          energy groups and polar angles, and tallies it into the FSR
 *          scalar flux, and updates the Track's angular flux.
 * @param curr_segment a pointer to the Track segment of interest
 * @param azim_index a pointer to the azimuthal angle index for this segment
 * @param track_flux a pointer to the Track's angular flux
 * @param fsr_flux a pointer to the temporary FSR flux buffer
 * @param fwd
 */
void CPUSolver::tallyScalarFlux(segment* curr_segment, int azim_index,
                                FP_PRECISION* track_flux,
                                FP_PRECISION* fsr_flux){

  int fsr_id = curr_segment->_region_id;
  FP_PRECISION length = curr_segment->_length;
  FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();
  FP_PRECISION delta_psi, exponential;

  /* Set the FSR scalar flux buffer to zero */
  memset(fsr_flux, 0.0, _num_groups * sizeof(FP_PRECISION));

  /* Compute change in angular flux along segment in this FSR */
  for (int e=0; e < _num_groups; e++) {
    for (int p=0; p < _num_polar; p++){
      exponential = _exp_evaluator->computeExponential(sigma_t[e] * length, p);
      delta_psi = (track_flux(p,e)-_reduced_sources(fsr_id,e)) * exponential;
      fsr_flux[e] += delta_psi * _polar_weights(azim_index,p);
      track_flux(p,e) -= delta_psi;
    }
  }

  /* Atomically increment the FSR scalar flux from the temporary array */
  omp_set_lock(&_FSR_locks[fsr_id]);
  {
    for (int e=0; e < _num_groups; e++)
      _scalar_flux(fsr_id,e) += fsr_flux[e];
  }
  omp_unset_lock(&_FSR_locks[fsr_id]);
}


/**
 * @brief Tallies the current contribution from this segment across the
 *        the appropriate CMFD mesh cell surface.
 * @param curr_segment a pointer to the Track segment of interest
 * @param azim_index a pointer to the azimuthal angle index for this segment
 * @param track_flux a pointer to the Track's angular flux
 * @param fwd
 */
void CPUSolver::tallySurfaceCurrent(segment* curr_segment,
                                    int azim_index,
                                    FP_PRECISION* track_flux,
                                    bool fwd){

  FP_PRECISION surf_current;

  /* Tally surface currents if CMFD is in use */
  if (_cmfd != NULL && _cmfd->isFluxUpdateOn()){
    for (int e=0; e < _num_groups; e++) {
      surf_current = 0.;

      for (int p=0; p < _num_polar; p++)
        surf_current += track_flux(p,e) * _polar_weights(azim_index,p);

      surf_current /= 2.;
      _cmfd->tallySurfaceCurrent(curr_segment, surf_current, fwd, e);
    }
  }
}


/**
 * @brief Updates the boundary flux for a Track given boundary conditions.
 * @details For reflective boundary conditions, the outgoing boundary flux
 *          for the Track is given to the reflecting Track. For vacuum
 *          boundary conditions, the outgoing flux tallied as leakage.
 * @param track_id the ID number for the Track of interest
 * @param azim_index a pointer to the azimuthal angle index for this segment
 * @param direction the Track direction (forward - true, reverse - false)
 * @param track_flux a pointer to the Track's outgoing angular flux
 */
void CPUSolver::transferBoundaryFlux(int track_id,
                                     int azim_index,
                                     bool direction,
                                     FP_PRECISION* track_flux) {
  int start;
  int bc;
  FP_PRECISION* track_leakage;
  int track_out_id;

  /* Extract boundary conditions for this Track and the pointer to the
   * outgoing reflective Track, and index into the leakage array */

  /* For the "forward" direction */
  if (direction) {
    start = _tracks[track_id]->isReflOut() * _polar_times_groups;
    bc = (int)_tracks[track_id]->getBCOut();
    track_leakage = &_boundary_leakage(track_id,0);
    track_out_id = _tracks[track_id]->getTrackOut()->getUid();
  }

  /* For the "reverse" direction */
  else {
    start = _tracks[track_id]->isReflIn() * _polar_times_groups;
    bc = (int)_tracks[track_id]->getBCIn();
    track_leakage = &_boundary_leakage(track_id,_polar_times_groups);
    track_out_id = _tracks[track_id]->getTrackIn()->getUid();
  }

  FP_PRECISION* track_out_flux = &_boundary_flux(track_out_id,0,0,start);

  /* Loop over polar angles and energy groups */
  for (int e=0; e < _num_groups; e++) {
    for (int p=0; p < _num_polar; p++) {
      track_out_flux(p,e) = track_flux(p,e) * bc;
      track_leakage(p,e) = track_flux(p,e) *
                           _polar_weights(azim_index,p) * (!bc);
    }
  }
}


/**
 * @brief Add the source term contribution in the transport equation to
 *        the FSR scalar flux.
 */
void CPUSolver::addSourceToScalarFlux() {

  FP_PRECISION volume;
  FP_PRECISION* sigma_t;

  /* Add in source term and normalize flux to volume for each FSR */
  /* Loop over FSRs, energy groups */
  #pragma omp parallel for private(volume, sigma_t) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    volume = _FSR_volumes[r];
    sigma_t = _FSR_materials[r]->getSigmaT();

    for (int e=0; e < _num_groups; e++) {
      _scalar_flux(r,e) *= 0.5;
      _scalar_flux(r,e) /= (sigma_t[e] * volume);
      _scalar_flux(r,e) += (FOUR_PI * _reduced_sources(r,e));
    }
  }
}


/**
 * @brief Computes the volume-averaged, energy-integrated fission rate in
 *        each FSR and stores them in an array indexed by FSR ID.
 * @details This is a helper method for SWIG to allow users to retrieve
 *          FSR fission rates as a NumPy array. An example of how this method 
 *          can be called from Python is as follows:
 *
 * @code
 *          num_FSRs = geometry.getNumFSRs()
 *          fission_rates = solver.computeFSRFissionRates(num_FSRs)
 * @endcode
 *
 * @param fission_rates an array to store the fission rates (implicitly passed
 *                      in as a NumPy array from Python)
 * @param num_FSRs the number of FSRs passed in from Python
 */
void CPUSolver::computeFSRFissionRates(double* fission_rates, int num_FSRs) {

  log_printf(INFO, "Computing FSR fission rates...");

  FP_PRECISION* sigma_f;
  FP_PRECISION* scalar_flux = getFSRScalarFluxes();

  /* Initialize fission rates to zero */
  for (int r=0; r < _num_FSRs; r++)
    fission_rates[r] = 0.0;

  /* Loop over all FSRs and compute the volume-averaged fission rate */
  #pragma omp parallel for private (sigma_f) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    sigma_f = _FSR_materials[r]->getSigmaF();

    for (int e=0; e < _num_groups; e++)
      fission_rates[r] += sigma_f[e] * _scalar_flux(r,e);
  }
}
