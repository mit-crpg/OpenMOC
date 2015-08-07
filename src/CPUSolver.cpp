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
 * @brief Fills an array with the scalar fluxes.
 * @details This class method is a helper routine called by the OpenMOC
 *          Python "openmoc.krylov" module for Krylov subspace methods. 
 *          Although this method appears to require two arguments, in
 *          reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_fluxes = num_groups * num_FSRs
 *          fluxes = solver.getFluxes(num_fluxes)
 * @endcode
 *
 * @param fluxes an array of FSR scalar fluxes in each energy group
 * @param num_fluxes the total number of FSR flux values
 */
void CPUSolver::getFluxes(FP_PRECISION* out_fluxes, int num_fluxes) {

  if (num_fluxes != _num_groups * _num_FSRs)
    log_printf(ERROR, "Unable to get FSR scalar fluxes since there are "
               "%d groups and %d FSRs which does not match the requested "
               "%d flux values", _num_groups, _num_FSRs, num_fluxes);

  else if (_scalar_flux == NULL)
    log_printf(ERROR, "Unable to get FSR scalar fluxes since they "
               "have not yet been allocated");

  /* If the user called setFluxes(...) they already have the flux */
  if (_user_fluxes && _scalar_flux == out_fluxes)
    return;

  /* Otherwise, copy the fluxes into the input array */
  else {
    #pragma omp parallel for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {
      for (int e=0; e < _num_groups; e++)
        out_fluxes[r*_num_groups+e] = _scalar_flux(r,e);
    }
  }
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
    log_printf(WARNING, "Over-riding fixed source %f in FSR ID=%d with %f",
               _fixed_sources(fsr_id,group-1), fsr_id, source);

  /* Store the fixed source for this FSR and energy group */
  _fixed_sources(fsr_id,group-1) = source;  
}


/**
 * @brief Set the flux array for use in transport sweep source calculations.
 * @detail This is a helper method for the checkpoint restart capabilities,
 *         as well as the IRAMSolver in the openmoc.krylov submodule. This
 *         routine may be used as follows from within Python:
 *
 * @code
 *          fluxes = numpy.random.rand(num_FSRs * num_groups, dtype=np.float)
 *          solver.setFluxes(fluxes)
 * @endcode
 *
 *          NOTE: This routine stores a pointer to the fluxes for the Solver
 *          to use during transport sweeps and other calculations. Hence, the 
 *          flux array pointer is shared between NumPy and the Solver.
 *
 * @param in_fluxes an array with the fluxes to use
 * @param num_fluxes the number of flux values (# groups x # FSRs)
 */
void CPUSolver::setFluxes(FP_PRECISION* in_fluxes, int num_fluxes) {
  if (num_fluxes != _num_groups * _num_FSRs)
    log_printf(ERROR, "Unable to set an array with %d flux values for %d "
               " groups and %d FSRs", num_fluxes, _num_groups, _num_FSRs);

  /* Allocate array if flux arrays have not yet been initialized */
  if (_scalar_flux == NULL)
    initializeFluxArrays();

  /* Set the scalar flux array pointer to the array passed in from NumPy */
  _scalar_flux = in_fluxes;
  _user_fluxes = true;
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
 * @brief Allocates memory for Track boundary angular and FSR scalar fluxes.
 * @details Deletes memory for old flux arrays if they were allocated
 *          for a previous simulation.
 */
void CPUSolver::initializeFluxArrays() {

  /* Delete old flux arrays if they exist */
  if (_boundary_flux != NULL)
    delete [] _boundary_flux;

  if (_scalar_flux != NULL)
    delete [] _scalar_flux;

  if (_old_scalar_flux != NULL)
    delete [] _old_scalar_flux;

  /* Allocate memory for the Track boundary flux arrays */
  try{
    int size = 2 * _tot_num_tracks * _polar_times_groups;
    _boundary_flux = new FP_PRECISION[size];

    /* Allocate an array for the FSR scalar flux */
    size = _num_FSRs * _num_groups;
    _scalar_flux = new FP_PRECISION[size];
    _old_scalar_flux = new FP_PRECISION[size];
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
    for (int e=0; e < _num_groups; e++) {
      _scalar_flux(r,e) *= norm_factor;
      _old_scalar_flux(r,e) *= norm_factor;
    }
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
 * @brief Computes the total source (fission, scattering, fixed) in each FSR.
 * @details This method computes the total source in each FSR based on
 *          this iteration's current approximation to the scalar flux.
 */
void CPUSolver::computeFSRSources() {

  #pragma omp parallel default(none)
  {
    int tid;
    Material* material;
    FP_PRECISION* sigma_t;
    FP_PRECISION sigma_s, fiss_mat;
    FP_PRECISION scatter_source, fission_source;
    FP_PRECISION* fission_sources = new FP_PRECISION[_num_groups];
    FP_PRECISION* scatter_sources = new FP_PRECISION[_num_groups];

    /* Compute the total source for each FSR */
    #pragma omp for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

      tid = omp_get_thread_num();
      material = _FSR_materials[r];
      sigma_t = material->getSigmaT();

      /* Compute scatter + fission source for group g */
      for (int g=0; g < _num_groups; g++) {
        for (int g_prime=0; g_prime < _num_groups; g_prime++) {
          sigma_s = material->getSigmaSByGroup(g_prime+1,g+1);
          fiss_mat = material->getFissionMatrixByGroup(g_prime+1,g+1);
          scatter_sources[g_prime] = sigma_s * _scalar_flux(r,g_prime);
          fission_sources[g_prime] = fiss_mat * _scalar_flux(r,g_prime);
        }

        scatter_source = pairwise_sum<FP_PRECISION>(scatter_sources, 
                                                    _num_groups);
        fission_source = pairwise_sum<FP_PRECISION>(fission_sources,
                                                    _num_groups);
        fission_source /= _k_eff;

        /* Compute total (scatter+fission+fixed) reduced source */
        _reduced_sources(r,g) = _fixed_sources(r,g);
        _reduced_sources(r,g) += scatter_source + fission_source;
        _reduced_sources(r,g) *= ONE_OVER_FOUR_PI / sigma_t[g];
      }
    }

    delete [] fission_sources;
    delete [] scatter_sources;
  }
}

/**
 * @brief Computes the total fission source in each FSR.
 * @details This method is a helper routine for the openmoc.krylov submodule.
 */
void CPUSolver::computeFSRFissionSources() {

  #pragma omp parallel default(none)
  {
    int tid;
    Material* material;
    FP_PRECISION* sigma_t;
    FP_PRECISION fiss_mat;
    FP_PRECISION fission_source;
    FP_PRECISION* fission_sources = new FP_PRECISION[_num_groups];

    /* Compute the total source for each FSR */
    #pragma omp for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

      tid = omp_get_thread_num();
      material = _FSR_materials[r];
      sigma_t = material->getSigmaT();

      /* Compute scatter + fission source for group g */
      for (int g=0; g < _num_groups; g++) {
        for (int g_prime=0; g_prime < _num_groups; g_prime++) {
          fiss_mat = material->getFissionMatrixByGroup(g_prime+1,g+1);
          fission_sources[g_prime] = fiss_mat * _scalar_flux(r,g_prime);
        }
        
        fission_source = pairwise_sum<FP_PRECISION>(fission_sources,
                                                    _num_groups);

        /* Compute total (fission) reduced source */
        _reduced_sources(r,g) = fission_source;
        _reduced_sources(r,g) *= ONE_OVER_FOUR_PI / sigma_t[g];
      }
    }

    delete [] fission_sources;
  }
}

/**
 * @brief Computes the total scattering source in each FSR.
 * @details This method is a helper routine for the openmoc.krylov submodule.
 */
void CPUSolver::computeFSRScatterSources() {
  
  #pragma omp parallel default(none)
  {
    int tid;
    Material* material;
    FP_PRECISION* sigma_t;
    FP_PRECISION sigma_s;
    FP_PRECISION scatter_source;
    FP_PRECISION* scatter_sources = new FP_PRECISION[_num_groups];

    /* Compute the total source for each FSR */
    #pragma omp for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

      tid = omp_get_thread_num();
      material = _FSR_materials[r];
      sigma_t = material->getSigmaT();

      /* Compute scatter + fission source for group g */
      for (int g=0; g < _num_groups; g++) {
        for (int g_prime=0; g_prime < _num_groups; g_prime++) {
          sigma_s = material->getSigmaSByGroup(g_prime+1,g+1);
          scatter_sources[g_prime] = sigma_s * _scalar_flux(r,g_prime);
        }

        scatter_source = pairwise_sum<FP_PRECISION>(scatter_sources, 
                                                    _num_groups);

        /* Compute total (scatter) reduced source */
        _reduced_sources(r,g) = scatter_source;
        _reduced_sources(r,g) *= ONE_OVER_FOUR_PI / sigma_t[g];
      }
    }

    delete [] scatter_sources;
  }
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
          new_total_source += material->getSigmaSByGroup(g+1,G+1)
                              * _scalar_flux(r,g);
          old_total_source += material->getSigmaSByGroup(g+1,G+1)
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

  FP_PRECISION old_fission, new_fission;
  FP_PRECISION* FSR_rates = new FP_PRECISION[_num_FSRs];
  FP_PRECISION* group_rates = new FP_PRECISION[_num_threads * _num_groups];

  /* Compute the old nu-fission rates in each FSR */
  #pragma omp parallel for private(tid, volume, \
    material, sigma) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {

    tid = omp_get_thread_num() * _num_groups;
    volume = _FSR_volumes[r];
    material = _FSR_materials[r];
    sigma = material->getNuSigmaF();

    for (int e=0; e < _num_groups; e++)
      group_rates[tid+e] = sigma[e] * _old_scalar_flux(r,e);

    FSR_rates[r]=pairwise_sum<FP_PRECISION>(&group_rates[tid], _num_groups);
    FSR_rates[r] *= volume;
  }

  /* Reduce old fission rates across FSRs */
  old_fission = pairwise_sum<FP_PRECISION>(FSR_rates, _num_FSRs);

  /* Compute the new nu-fission rates in each FSR */
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
  new_fission = pairwise_sum<FP_PRECISION>(FSR_rates, _num_FSRs);

  _k_eff *= new_fission / old_fission;

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
    _cmfd->zeroCurrents();

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
        tallyCurrent(curr_segment, azim_index, track_flux, true);
      }

      /* Transfer boundary angular flux to outgoing Track */
      transferBoundaryFlux(track_id, azim_index, true, track_flux);

      /* Loop over each Track segment in reverse direction */
      track_flux += _polar_times_groups;

      for (int s=num_segments-1; s > -1; s--) {
        curr_segment = &segments[s];
        tallyScalarFlux(curr_segment, azim_index, track_flux, thread_fsr_flux);
        tallyCurrent(curr_segment, azim_index, track_flux, false);
      }
      delete [] thread_fsr_flux;

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
                                FP_PRECISION* fsr_flux) {

  int fsr_id = curr_segment->_region_id;
  FP_PRECISION length = curr_segment->_length;
  FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();
  FP_PRECISION delta_psi, exponential;

  /* Set the FSR scalar flux buffer to zero */
  memset(fsr_flux, 0.0, _num_groups * sizeof(FP_PRECISION));

  /* Compute change in angular flux along segment in this FSR */
  for (int e=0; e < _num_groups; e++) {
    for (int p=0; p < _num_polar; p++) {
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
 *        the appropriate CMFD mesh cell surface or corner.
 * @param curr_segment a pointer to the Track segment of interest
 * @param azim_index the azimuthal index for this segmenbt
 * @param track_flux a pointer to the Track's angular flux
 * @param fwd boolean indicating direction of integration along segment
 */
void CPUSolver::tallyCurrent(segment* curr_segment, int azim_index,
                             FP_PRECISION* track_flux, bool fwd) {

  /* Tally surface or corner currents if CMFD is in use */
  if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
    _cmfd->tallyCurrent(curr_segment, track_flux,
                        &_polar_weights(azim_index,0), fwd);
}


/**
 * @brief Updates the boundary flux for a Track given boundary conditions.
 * @details For reflective boundary conditions, the outgoing boundary flux
 *          for the Track is given to the reflecting Track.
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
  int track_out_id;

  /* Extract boundary conditions for this Track */

  /* For the "forward" direction */
  if (direction) {
    start = _tracks[track_id]->isReflOut() * _polar_times_groups;
    bc = (int)_tracks[track_id]->getBCOut();
    track_out_id = _tracks[track_id]->getTrackOut()->getUid();
  }

  /* For the "reverse" direction */
  else {
    start = _tracks[track_id]->isReflIn() * _polar_times_groups;
    bc = (int)_tracks[track_id]->getBCIn();
    track_out_id = _tracks[track_id]->getTrackIn()->getUid();
  }

  FP_PRECISION* track_out_flux = &_boundary_flux(track_out_id,0,0,start);

  /* Loop over polar angles and energy groups */
  for (int e=0; e < _num_groups; e++) {
    for (int p=0; p < _num_polar; p++)
      track_out_flux(p,e) = track_flux(p,e) * bc;
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

  /* Loop over all FSRs and compute the volume-averaged nu-fission rate */
  #pragma omp parallel for private (nu_sigma_f) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    nu_sigma_f = _FSR_materials[r]->getNuSigmaF();

    for (int e=0; e < _num_groups; e++)
      fission_rates[r] += nu_sigma_f[e] * _scalar_flux(r,e);
  }
}
