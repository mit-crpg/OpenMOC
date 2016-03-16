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
CPUSolver::~CPUSolver() {}


/**
 * @brief Returns the number of shared memory OpenMP threads in use.
 * @return the number of threads
 */
int CPUSolver::getNumThreads() {
  return _num_threads;
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

  /* Warn the user if a fixed source has already been assigned to this FSR */
  if (_fixed_sources(fsr_id,group-1) != 0.)
    log_printf(WARNING, "Overriding fixed source %f in FSR ID=%d with %f",
               _fixed_sources(fsr_id,group-1), fsr_id, source);

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

  /* Get FSR locks from TrackGenerator */
  _FSR_locks = _track_generator->getFSRLocks();
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

  if (_start_flux != NULL)
    delete [] _start_flux;

  if (_boundary_leakage != NULL)
    delete [] _boundary_leakage;

  if (_scalar_flux != NULL)
    delete [] _scalar_flux;

  if (_old_scalar_flux != NULL)
    delete [] _old_scalar_flux;

  int size;

  /* Allocate memory for the Track boundary flux and leakage arrays */
  try {
    size = 2 * _tot_num_tracks * _fluxes_per_track;

    _boundary_flux = new FP_PRECISION[size];
    _start_flux = new FP_PRECISION[size];
    _boundary_leakage = new FP_PRECISION[size];

    /* Allocate an array for the FSR scalar flux */
    size = _num_FSRs * _num_groups;
    _scalar_flux = new FP_PRECISION[size];
    _old_scalar_flux = new FP_PRECISION[size];
    memset(_scalar_flux, 0., size * sizeof(FP_PRECISION));
    memset(_old_scalar_flux, 0., size * sizeof(FP_PRECISION));
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
  if (_reduced_sources != NULL)
    delete [] _reduced_sources;

  /* Allocate memory for all source arrays */
  try{
    int size = _num_FSRs * _num_groups;
    _reduced_sources = new FP_PRECISION[size];

    /* If no fixed sources were assigned, use a zeroes array */
    if (_fixed_sources == NULL) {
      _fixed_sources = new FP_PRECISION[size];
      memset(_fixed_sources, 0.0, sizeof(FP_PRECISION) * size);
    }
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
      for (int pe=0; pe < _fluxes_per_track; pe++) {
        _boundary_flux(t, d, pe) = 0.0;
        _start_flux(t, d, pe) = 0.0;
      }
    }
  }
}


/**
 * @brief Copies values from the start flux into the boundary flux array
 *        for both the "forward" and "reverse" directions.
 */
void CPUSolver::copyBoundaryFluxes() {

  #pragma omp parallel for schedule(guided)
  for (int t=0; t < _tot_num_tracks; t++) {
    for (int d=0; d < 2; d++) {
      for (int pe=0; pe < _fluxes_per_track; pe++)
        _boundary_flux(t,d,pe) = _start_flux(t, d, pe);
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
 * @brief Stores the FSR scalar fluxes in the old scalar flux array.
 */
void CPUSolver::storeFSRFluxes() {

  #pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _num_groups; e++)
      _old_scalar_flux(r,e) = _scalar_flux(r,e);
  }
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

  int size = _num_FSRs * _num_groups;
  FP_PRECISION* fission_sources = new FP_PRECISION[_num_FSRs * _num_groups];

  /* Compute total fission source for each FSR, energy group */
  #pragma omp parallel for private(volume, nu_sigma_f) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {

    /* Get pointers to important data structures */
    nu_sigma_f = _FSR_materials[r]->getNuSigmaF();
    volume = _FSR_volumes[r];

    for (int e=0; e < _num_groups; e++)
      fission_sources(r,e) = nu_sigma_f[e] * _scalar_flux(r,e) * volume;
  }

  /* Compute the total fission source */

  tot_fission_source = pairwise_sum<FP_PRECISION>(fission_sources,size);

  /* Deallocate memory for fission source array */
  delete [] fission_sources;

  /* Normalize scalar fluxes in each FSR */
  norm_factor = 1.0 / tot_fission_source;

  log_printf(DEBUG, "Tot. Fiss. Src. = %f, Norm. factor = %f",
             tot_fission_source, norm_factor);

  #pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _num_groups; e++)
      _scalar_flux(r, e) *= norm_factor;
  }

  /* Normalize angular boundary fluxes for each Track */
  #pragma omp parallel for schedule(guided)
  for (int i=0; i < _tot_num_tracks; i++) {
    for (int j=0; j < 2; j++) {
      for (int pe=0; pe < _fluxes_per_track; pe++) {
        _start_flux(i, j, pe) *= norm_factor;
        _boundary_flux(i, j, pe) *= norm_factor;
      }
    }
  }
}


/**
 * @brief Computes the total source (fission, scattering, fixed) in each FSR.
 * @details This method computes the total source in each FSR based on
 *          this iteration's current approximation to the scalar flux.
 */
void CPUSolver::computeFSRSources() {

  int tid;
  FP_PRECISION scatter_source, fission_source;
  FP_PRECISION* nu_sigma_f;
  FP_PRECISION* sigma_t;
  FP_PRECISION* chi;
  Material* material;

  int size = _num_FSRs * _num_groups;
  FP_PRECISION* fission_sources = new FP_PRECISION[size];
  size = _num_threads * _num_groups;
  FP_PRECISION* scatter_sources = new FP_PRECISION[size];

  /* For all FSRs, find the source */
  #pragma omp parallel for private(tid, material, nu_sigma_f, chi, \
    sigma_t, fission_source, scatter_source) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {

    tid = omp_get_thread_num();
    material = _FSR_materials[r];
    nu_sigma_f = material->getNuSigmaF();
    chi = material->getChi();
    sigma_t = material->getSigmaT();

    /* Initialize the fission sources to zero */
    fission_source = 0.0;

    /* Compute fission source for each group */
    if (material->isFissionable()) {
      for (int e=0; e < _num_groups; e++)
        fission_sources(r,e) = _scalar_flux(r,e) * nu_sigma_f[e];

      fission_source = pairwise_sum<FP_PRECISION>(&fission_sources(r,0),
                                                  _num_groups);
      fission_source /= _k_eff;
    }

    /* Compute total (fission+scatter+fixed) source for group G */
    for (int G=0; G < _num_groups; G++) {
      for (int g=0; g < _num_groups; g++)
        scatter_sources(tid,g) = material->getSigmaSByGroupInline(g,G)
                                  * _scalar_flux(r,g);
      scatter_source = pairwise_sum<FP_PRECISION>(&scatter_sources(tid,0),
                                                _num_groups);

      _reduced_sources(r,G) = fission_source * chi[G];
      _reduced_sources(r,G) += scatter_source + _fixed_sources(r,G);
      _reduced_sources(r,G) *= ONE_OVER_FOUR_PI / sigma_t[G];
    }
  }

  delete [] fission_sources;
  delete [] scatter_sources;
}


/**
 * @brief Computes the residual between source/flux iterations.
 * @param res_type the type of residuals to compute
 *        (SCALAR_FLUX, FISSION_SOURCE, TOTAL_SOURCE)
 * @return the average residual in each FSR
 */
double CPUSolver::computeResidual(residualType res_type) {

  int norm;
  double residual;
  double* residuals = new double[_num_FSRs];
  memset(residuals, 0., _num_FSRs * sizeof(double));

  if (res_type == SCALAR_FLUX) {

    norm = _num_FSRs;

    for (int r=0; r < _num_FSRs; r++) {
      for (int e=0; e < _num_groups; e++)
        if (_old_scalar_flux(r,e) > 0.) {
          residuals[r] += pow((_scalar_flux(r,e) - _old_scalar_flux(r,e)) /
                              _old_scalar_flux(r,e), 2);
      }
    }
  }

  else if (res_type == FISSION_SOURCE) {

    if (_num_fissionable_FSRs == 0)
      log_printf(ERROR, "The Solver is unable to compute a "
                 "FISSION_SOURCE residual without fissionable FSRs");

    norm = _num_fissionable_FSRs;

    double new_fission_source, old_fission_source;
    FP_PRECISION* nu_sigma_f;
    Material* material;

    for (int r=0; r < _num_FSRs; r++) {
      new_fission_source = 0.;
      old_fission_source = 0.;
      material = _FSR_materials[r];

      if (material->isFissionable()) {
        nu_sigma_f = material->getNuSigmaF();

        for (int e=0; e < _num_groups; e++) {
          new_fission_source += _scalar_flux(r,e) * nu_sigma_f[e];
          old_fission_source += _old_scalar_flux(r,e) * nu_sigma_f[e];
        }

        if (old_fission_source > 0.)
          residuals[r] = pow((new_fission_source -  old_fission_source) /
                              old_fission_source, 2);
      }
    }
  }

  else if (res_type == TOTAL_SOURCE) {

    norm = _num_FSRs;

    double new_total_source, old_total_source;
    FP_PRECISION inverse_k_eff = 1.0 / _k_eff;
    FP_PRECISION* nu_sigma_f;
    Material* material;

    for (int r=0; r < _num_FSRs; r++) {
      new_total_source = 0.;
      old_total_source = 0.;
      material = _FSR_materials[r];

      if (material->isFissionable()) {
        nu_sigma_f = material->getNuSigmaF();

        for (int e=0; e < _num_groups; e++) {
          new_total_source += _scalar_flux(r,e) * nu_sigma_f[e];
          old_total_source += _old_scalar_flux(r,e) * nu_sigma_f[e];
        }

        new_total_source *= inverse_k_eff;
        old_total_source *= inverse_k_eff;
      }

      /* Compute total scattering source for group G */
      for (int G=0; G < _num_groups; G++) {
        for (int g=0; g < _num_groups; g++) {
          new_total_source += material->getSigmaSByGroupInline(g,G)
                              * _scalar_flux(r,g);
          old_total_source += material->getSigmaSByGroupInline(g,G)
                              * _old_scalar_flux(r,g);
        }
      }

      if (old_total_source > 0.)
        residuals[r] = pow((new_total_source -  old_total_source) /
                            old_total_source, 2);
    }
  }

  /* Sum up the residuals from each FSR and normalize */
  residual = pairwise_sum<double>(residuals, _num_FSRs);
  residual = sqrt(residual / norm);

  /* Deallocate memory for residuals array */
  delete [] residuals;

  return residual;
}


/**
 * @brief Compute \f$ k_{eff} \f$ from successive fission sources.
 */
void CPUSolver::computeKeff() {

  int tid;
  Material* material;
  FP_PRECISION* sigma;
  FP_PRECISION volume;

  FP_PRECISION fission;
  FP_PRECISION* FSR_rates = new FP_PRECISION[_num_FSRs];
  FP_PRECISION* group_rates = new FP_PRECISION[_num_threads * _num_groups];

  /* Loop over all FSRs and compute the volume-integrated total rates */
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

  /* Reduce new fission rates across FSRs */
  fission = pairwise_sum<FP_PRECISION>(FSR_rates, _num_FSRs);

  /* The old_source is normalized to sum to _k_eff; therefore, the new
   * _k_eff is the old _k_eff * sum(new_source). Implicitly, we are just doing
   * _k_eff = sum(new_source) / sum(old_source). */
  _k_eff *= fission;

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

  log_printf(DEBUG, "Transport sweep with %d OpenMP threads",
      _num_threads);

  if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
    _cmfd->zeroCurrents();

  /* Initialize flux in each FSR to zero */
  flattenFSRFluxes(0.0);

  /* Copy starting flux to current flux */
  copyBoundaryFluxes();

  /* Tracks are traversed and the MOC equations from this CPUSolver are applied
     to all Tracks and corresponding segments */
  TransportSweep sweep_tracks(_track_generator);
  sweep_tracks.setCPUSolver(this);
  sweep_tracks.execute();
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
 */
void CPUSolver::tallyScalarFlux(segment* curr_segment,
                                int azim_index, int polar_index,
                                FP_PRECISION* track_flux,
                                FP_PRECISION* fsr_flux) {

  int fsr_id = curr_segment->_region_id;
  FP_PRECISION length = curr_segment->_length;
  FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();
  int p, pe;
  int a = azim_index;

  /* The change in angular flux along this Track segment in the FSR */
  FP_PRECISION delta_psi, exponential, tau;

  /* Set the FSR scalar flux buffer to zero */
  memset(fsr_flux, 0.0, _num_groups * sizeof(FP_PRECISION));

  if (_solve_3D) {

    for (int e=0; e < _num_groups; e++) {
      exponential = _exp_evaluator->computeExponential
        (sigma_t[e] * length, a, p);
      delta_psi = (track_flux(e)-_reduced_sources(fsr_id, e)) * exponential;
      fsr_flux[e] += delta_psi * _quad->getWeight(a, p);
      track_flux(e) -= delta_psi;
    }
  }
  else {

    pe = 0;

    /* Loop over energy groups */
    for (int e=0; e < _num_groups; e++) {

      tau = sigma_t[e] * length;

      /* Loop over polar angles */
      for (p=0; p < _num_polar/2; p++) {
        exponential = _exp_evaluator->computeExponential(tau, a, p);
        delta_psi = (track_flux(pe)-_reduced_sources(fsr_id,e)) * exponential;
        fsr_flux[e] += delta_psi * _quad->getWeight(a, p);
        track_flux(pe) -= delta_psi;
        pe++;
      }
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
 * @param azim_index the azimuthal index for this segmenbt
 * @param track_flux a pointer to the Track's angular flux
 * @param fwd boolean indicating direction of integration along segment
 */
void CPUSolver::tallyCurrent(segment* curr_segment, int azim_index,
                             int polar_index, FP_PRECISION* track_flux,
                             bool fwd) {

  /* Tally surface currents if CMFD is in use */
  if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
    _cmfd->tallyCurrent(curr_segment, track_flux, azim_index, polar_index, fwd);
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
                                     int azim_index, int polar_index,
                                     bool direction,
                                     FP_PRECISION* track_flux) {
  int start;
  boundaryType bc;
  FP_PRECISION* track_leakage;
  int track_out_id;
  int a = azim_index;

  /* Extract boundary conditions for this Track and the pointer to the
   * outgoing reflective Track, and index into the leakage array */

  /* For the "forward" direction */
  if (direction) {
    bc = _tracks[track_id]->getBCFwd();
    track_leakage = &_boundary_leakage(track_id, 0);
    if (bc == PERIODIC) {
      start = 0;
      track_out_id = _tracks[track_id]->getTrackPrdcFwd()->getUid();
    }
    else {
      start = _fluxes_per_track * (!_tracks[track_id]->getReflFwdFwd());
      track_out_id = _tracks[track_id]->getTrackReflFwd()->getUid();
    }
  }

  /* For the "reverse" direction */
  else {
    bc = _tracks[track_id]->getBCBwd();
    track_leakage = &_boundary_leakage(track_id,_fluxes_per_track);
    if (bc == PERIODIC) {
      start = _fluxes_per_track;
      track_out_id = _tracks[track_id]->getTrackPrdcBwd()->getUid();
    }
    else {
      start = _fluxes_per_track * (!_tracks[track_id]->getReflBwdFwd());
      track_out_id = _tracks[track_id]->getTrackReflBwd()->getUid();
    }
  }

  FP_PRECISION* track_out_flux = &_start_flux(track_out_id, 0, start);

  /* Set bc to 1 if bc is PERIODIC (bc == 2) */
  if (bc == PERIODIC)
    bc = REFLECTIVE;

  if (_solve_3D) {
    int p = _quad->getFirstOctantPolar(polar_index);
    for (int e=0; e < _num_groups; e++) {
      track_out_flux(e) = track_flux(e) * bc;
      track_leakage(e) = track_flux(e) * (!bc) * _quad->getWeight(a, p);
    }
  }
  else {

    int pe = 0;

    /* Loop over polar angles and energy groups */
    for (int e=0; e < _num_groups; e++) {
      for (int p=0; p < _num_polar/2; p++) {
        track_out_flux(pe) = track_flux(pe) * bc;
        track_leakage(pe) = track_flux(pe) * (!bc) * _quad->getWeight(a, p);
        pe++;
      }
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
      _scalar_flux(r, e) *= 0.5;
      _scalar_flux(r, e) /= (sigma_t[e] * volume);
      _scalar_flux(r, e) += FOUR_PI * _reduced_sources(r, e);
    }
  }

  return;
}


/**
 * @brief Computes the volume-averaged, energy-integrated nu-fission rate in
 *        each FSR and stores them in an array indexed by FSR ID.
 * @details This is a helper method for SWIG to allow users to retrieve
 *          FSR nu-fission rates as a NumPy array. An example of how this
 *          method can be called from Python is as follows:
 *
 * @code
 *          num_FSRs = geometry.getNumFSRs()
 *          fission_rates = solver.computeFSRFissionRates(num_FSRs)
 * @endcode
 *
 * @param fission_rates an array to store the nu-fission rates (implicitly
 *                      passed in as a NumPy array from Python)
 * @param num_FSRs the number of FSRs passed in from Python
 */
void CPUSolver::computeFSRFissionRates(double* fission_rates, int num_FSRs) {

  if (_scalar_flux == NULL)
    log_printf(ERROR, "Unable to compute FSR fission rates since the "
               "source distribution has not been calculated");

  log_printf(INFO, "Computing FSR fission rates...");

  FP_PRECISION* nu_sigma_f;

  /* Initialize fission rates to zero */
  for (int r=0; r < _num_FSRs; r++)
    fission_rates[r] = 0.0;

  /* Loop over all FSRs and compute the volume-weighted fission rate */
  #pragma omp parallel for private (nu_sigma_f) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    nu_sigma_f = _FSR_materials[r]->getNuSigmaF();

    for (int e=0; e < _num_groups; e++)
      fission_rates[r] += nu_sigma_f[e] * _scalar_flux(r,e);
  }
}
