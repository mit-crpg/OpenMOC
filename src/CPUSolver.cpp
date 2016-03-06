#include "CPUSolver.h"


/**
 * @brief Constructor initializes array pointers for Tracks and Materials.
 * @details The constructor retrieves the number of energy groups and SRs
 *          and azimuthal angles from the Geometry and TrackGenerator if
 *          passed in as parameters by the user. The constructor initalizes
 *          the number of OpenMP threads to a default of 1.
 * @param track_generator an optional pointer to the TrackGenerator
 */
CPUSolver::CPUSolver(TrackGenerator* track_generator)
  : Solver(track_generator) {

  setNumThreads(1);
  _SR_locks = NULL;
}


/**
 * @brief Destructor deletes array for OpenMP mutual exclusion locks for
 *        SR scalar flux updates, and calls Solver parent class destructor
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
 * @brief Fills an array with the scalar fluxes.
 * @details This class method is a helper routine called by the OpenMOC
 *          Python "openmoc.krylov" module for Krylov subspace methods.
 *          Although this method appears to require two arguments, in
 *          reality it only requires one due to SWIG and would be called
 *          from within Python as follows:
 *
 * @code
 *          num_fluxes = num_groups * num_SRs
 *          fluxes = solver.getFluxes(num_fluxes)
 * @endcode
 *
 * @param fluxes an array of SR scalar fluxes in each energy group
 * @param num_fluxes the total number of SR flux values
 */
void CPUSolver::getFluxes(FP_PRECISION* out_fluxes, int num_fluxes) {

  if (num_fluxes != _num_groups * _num_SRs)
    log_printf(ERROR, "Unable to get SR scalar fluxes since there are "
               "%d groups and %d SRs which does not match the requested "
               "%d flux values", _num_groups, _num_SRs, num_fluxes);

  else if (_scalar_flux == NULL)
    log_printf(ERROR, "Unable to get SR scalar fluxes since they "
               "have not yet been allocated");

  /* If the user called setFluxes(...) they already have the flux */
  if (_user_fluxes && _scalar_flux == out_fluxes)
    return;

  /* Otherwise, copy the fluxes into the input array */
  else {
#pragma omp parallel for schedule(guided)
    for (int r=0; r < _num_SRs; r++) {
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
 * @brief Set the flux array for use in transport sweep source calculations.
 * @detail This is a helper method for the checkpoint restart capabilities,
 *         as well as the IRAMSolver in the openmoc.krylov submodule. This
 *         routine may be used as follows from within Python:
 *
 * @code
 *          fluxes = numpy.random.rand(num_SRs * num_groups, dtype=np.float)
 *          solver.setFluxes(fluxes)
 * @endcode
 *
 *          NOTE: This routine stores a pointer to the fluxes for the Solver
 *          to use during transport sweeps and other calculations. Hence, the
 *          flux array pointer is shared between NumPy and the Solver.
 *
 * @param in_fluxes an array with the fluxes to use
 * @param num_fluxes the number of flux values (# groups x # SRs)
 */
void CPUSolver::setFluxes(FP_PRECISION* in_fluxes, int num_fluxes) {
  if (num_fluxes != _num_groups * _num_SRs)
    log_printf(ERROR, "Unable to set an array with %d flux values for %d "
               " groups and %d SRs", num_fluxes, _num_groups, _num_SRs);

  /* Allocate array if flux arrays have not yet been initialized */
  if (_scalar_flux == NULL)
    initializeFluxArrays();

  /* Set the scalar flux array pointer to the array passed in from NumPy */
  _scalar_flux = in_fluxes;
  _user_fluxes = true;
}


/**
 * @brief Initializes the SR volumes and Materials array.
 * @details This method gets an array of OpenMP mutual exclusion locks
 *          for each SR for use in the transport sweep algorithm.
 */
void CPUSolver::initializeSRs() {
  Solver::initializeSRs();
  _SR_locks = _track_generator->getSRLocks();
}


/**
 * @brief Allocates memory for Track boundary angular and SR scalar fluxes.
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

    /* Allocate an array for the SR scalar flux */
    size = _num_SRs * _num_groups;
    _scalar_flux = new FP_PRECISION[size];
    _old_scalar_flux = new FP_PRECISION[size];
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the fluxes");
  }
}


/**
 * @brief Allocates memory for SR source arrays.
 * @details Deletes memory for old source arrays if they were allocated for a
 *          previous simulation.
 */
void CPUSolver::initializeSourceArrays() {

  /* Delete old sources arrays if they exist */
  if (_reduced_sources != NULL)
    delete [] _reduced_sources;
  if (_fixed_sources != NULL)
    delete [] _fixed_sources;

  int size = _num_SRs * _num_groups;

  /* Allocate memory for all source arrays */
  try{
    _reduced_sources = new FP_PRECISION[size];
    _fixed_sources = new FP_PRECISION[size];
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for SR sources");
  }

  /* Initialize fixed sources to zero */
  memset(_fixed_sources, 0.0, sizeof(FP_PRECISION) * size);

  /* Populate fixed source array with any user-defined sources */
  initializeFixedSources();
}


/**
 * @brief Populates array of fixed sources assigned by SR.
 */
void CPUSolver::initializeFixedSources() {

  Solver::initializeFixedSources();

  int sr_id, group;
  std::pair<int, int> sr_group_key;
  std::map< std::pair<int, int>, FP_PRECISION >::iterator sr_iter;

  /* Populate fixed source array with any user-defined sources */
  for (sr_iter = _fix_src_SR_map.begin();
       sr_iter != _fix_src_SR_map.end(); ++sr_iter) {

    /* Get the SR with an assigned fixed source */
    sr_group_key = sr_iter->first;
    sr_id = sr_group_key.first;
    group = sr_group_key.second;

    if (group <= 0 || group > _num_groups)
      log_printf(ERROR,"Unable to use fixed source for group %d in "
                 "a %d energy group problem", group, _num_groups);

    if (sr_id < 0 || sr_id >= _num_SRs)
      log_printf(ERROR,"Unable to use fixed source for SR %d with only "
                 "%d SRs in the geometry", sr_id, _num_SRs);

    _fixed_sources(sr_id, group-1) = _fix_src_SR_map[sr_group_key];
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
 * @brief Set the scalar flux for each SR and energy group to some value.
 * @param value the value to assign to each SR scalar flux
 */
void CPUSolver::flattenSRFluxes(FP_PRECISION value) {

#pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_SRs; r++) {
    for (int e=0; e < _num_groups; e++)
      _scalar_flux(r,e) = value;
  }
}


/**
 * @brief Stores the SR scalar fluxes in the old scalar flux array.
 */
void CPUSolver::storeSRFluxes() {

#pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_SRs; r++) {
    for (int e=0; e < _num_groups; e++)
      _old_scalar_flux(r,e) = _scalar_flux(r,e);
  }
}


/**
 * @brief Normalizes all SR scalar fluxes and Track boundary angular
 *        fluxes to the total fission source (times \f$ \nu \f$).
 */
void CPUSolver::normalizeFluxes() {

  FP_PRECISION* nu_sigma_f;
  FP_PRECISION volume;
  FP_PRECISION tot_fission_source;
  FP_PRECISION norm_factor;

  int size = _num_SRs * _num_groups;
  FP_PRECISION* fission_sources = new FP_PRECISION[_num_SRs * _num_groups];

  /* Compute total fission source for each SR, energy group */
#pragma omp parallel for private(volume, nu_sigma_f) schedule(guided)
  for (int r=0; r < _num_SRs; r++) {

    /* Get pointers to important data structures */
    nu_sigma_f = _SR_materials[r]->getNuSigmaF();
    volume = _SR_volumes[r];

    for (int e=0; e < _num_groups; e++)
      fission_sources(r,e) = nu_sigma_f[e] * _scalar_flux(r,e) * volume;
  }

  /* Compute the total fission source */
  tot_fission_source = pairwise_sum<FP_PRECISION>(fission_sources,size);

  /* Deallocate memory for fission source array */
  delete [] fission_sources;

  /* Normalize scalar fluxes in each SR */
  norm_factor = 1.0 / tot_fission_source;

  log_printf(DEBUG, "Tot. Fiss. Src. = %f, Norm. factor = %f",
             tot_fission_source, norm_factor);

#pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_SRs; r++) {
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
 * @brief Computes the total source (fission, scattering, fixed) in each SR.
 * @details This method computes the total source in each SR based on
 *          this iteration's current approximation to the scalar flux.
 */
void CPUSolver::computeSRSources() {

#pragma omp parallel default(none)
  {
    int tid;
    Material* material;
    FP_PRECISION* sigma_t;
    FP_PRECISION sigma_s, fiss_mat;
    FP_PRECISION scatter_source, fission_source;
    FP_PRECISION* fission_sources = new FP_PRECISION[_num_groups];
    FP_PRECISION* scatter_sources = new FP_PRECISION[_num_groups];

    /* Compute the total source for each SR */
#pragma omp for schedule(guided)
    for (int r=0; r < _num_SRs; r++) {

      tid = omp_get_thread_num();
      material = _SR_materials[r];
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
 * @brief Computes the total fission source in each SR.
 * @details This method is a helper routine for the openmoc.krylov submodule.
 */
void CPUSolver::computeSRFissionSources() {

#pragma omp parallel default(none)
  {
    int tid;
    Material* material;
    FP_PRECISION* sigma_t;
    FP_PRECISION fiss_mat;
    FP_PRECISION fission_source;
    FP_PRECISION* fission_sources = new FP_PRECISION[_num_groups];

    /* Compute the total source for each SR */
#pragma omp for schedule(guided)
    for (int r=0; r < _num_SRs; r++) {

      tid = omp_get_thread_num();
      material = _SR_materials[r];
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
 * @brief Computes the total scattering source in each SR.
 * @details This method is a helper routine for the openmoc.krylov submodule.
 */
void CPUSolver::computeSRScatterSources() {

#pragma omp parallel default(none)
  {
    int tid;
    Material* material;
    FP_PRECISION* sigma_t;
    FP_PRECISION sigma_s;
    FP_PRECISION scatter_source;
    FP_PRECISION* scatter_sources = new FP_PRECISION[_num_groups];

    /* Compute the total source for each SR */
#pragma omp for schedule(guided)
    for (int r=0; r < _num_SRs; r++) {

      tid = omp_get_thread_num();
      material = _SR_materials[r];
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
 * @return the average residual in each SR
 */
double CPUSolver::computeResidual(residualType res_type) {

  int norm;
  double residual;
  double* residuals = new double[_num_SRs];
  memset(residuals, 0., _num_SRs * sizeof(double));

  if (res_type == SCALAR_FLUX) {

    norm = _num_SRs;

#pragma omp parallel for schedule(guided)
    for (int r=0; r < _num_SRs; r++) {
      for (int e=0; e < _num_groups; e++)
        if (_old_scalar_flux(r,e) > 0.) {
          residuals[r] += pow((_scalar_flux(r,e) - _old_scalar_flux(r,e)) /
                               _old_scalar_flux(r,e), 2);
      }
    }
  }

  else if (res_type == FISSION_SOURCE) {

    if (_num_fissionable_SRs == 0)
      log_printf(ERROR, "The Solver is unable to compute a "
                 "FISSION_SOURCE residual without fissionable SRs");

    norm = _num_fissionable_SRs;

#pragma omp parallel
    {

      double new_fission_source, old_fission_source;
      FP_PRECISION* nu_sigma_f;
      Material* material;

#pragma omp for schedule(guided)
      for (int r=0; r < _num_SRs; r++) {
        new_fission_source = 0.;
        old_fission_source = 0.;
        material = _SR_materials[r];

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
  }

  else if (res_type == TOTAL_SOURCE) {

    norm = _num_SRs;

#pragma omp parallel
    {

      double new_total_source, old_total_source;
      FP_PRECISION inverse_k_eff = 1.0 / _k_eff;
      FP_PRECISION* nu_sigma_f;
      Material* material;
      FP_PRECISION* sigma_s;

#pragma omp for schedule(guided)
      for (int r=0; r < _num_SRs; r++) {
        new_total_source = 0.;
        old_total_source = 0.;
        material = _SR_materials[r];
        sigma_s = material->getSigmaS();

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
            new_total_source += sigma_s[G*_num_groups+g]
                * _scalar_flux(r,g);
            old_total_source += sigma_s[G*_num_groups+g]
                * _old_scalar_flux(r,g);
          }
        }

        if (old_total_source > 0.)
          residuals[r] = pow((new_total_source -  old_total_source) /
                             old_total_source, 2);
      }
    }
  }

  /* Sum up the residuals from each SR and normalize */
  residual = pairwise_sum<double>(residuals, _num_SRs);
  residual = sqrt(residual / norm);

  /* Deallocate memory for residuals array */
  delete [] residuals;

  return residual;
}


/**
 * @brief Compute \f$ k_{eff} \f$ from successive fission sources.
 */
void CPUSolver::computeKeff() {

  FP_PRECISION fission;
  FP_PRECISION* SR_rates = new FP_PRECISION[_num_SRs];
  FP_PRECISION* group_rates = new FP_PRECISION[_num_threads * _num_groups];

  /* Compute the old nu-fission rates in each SR */
#pragma omp parallel
  {

    int tid = omp_get_thread_num() * _num_groups;
    Material* material;
    FP_PRECISION* sigma;
    FP_PRECISION volume;

#pragma omp for schedule(guided)
    for (int r=0; r < _num_SRs; r++) {

      volume = _SR_volumes[r];
      material = _SR_materials[r];
      sigma = material->getNuSigmaF();

      for (int e=0; e < _num_groups; e++)
        group_rates[tid+e] = sigma[e] * _scalar_flux(r,e);

      SR_rates[r]=pairwise_sum<FP_PRECISION>(&group_rates[tid], _num_groups);
      SR_rates[r] *= volume;
    }
  }

  /* Reduce new fission rates across SRs */
  fission = pairwise_sum<FP_PRECISION>(SR_rates, _num_SRs);

  _k_eff *= fission;

  delete [] SR_rates;
  delete [] group_rates;
}


/**
 * @brief This method performs one transport sweep of all azimuthal angles,
 *        Tracks, Track segments, polar angles and energy groups.
 * @details The method integrates the flux along each Track and updates the
 *          boundary fluxes for the corresponding output Track, while updating
 *          the scalar flux in each source region.
 */
void CPUSolver::transportSweep() {

  log_printf(DEBUG, "Transport sweep with %d OpenMP threads", _num_threads);

  int min_track = 0;
  int max_track = 0;

  /* Initialize flux in each Sr to zero */
  flattenSRFluxes(0.0);

  if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
    _cmfd->zeroCurrents();

  /* Loop over the parallel track groups */
  for (int i=0; i < _num_parallel_track_groups; i++) {

    /* Compute the minimum and maximum Track IDs corresponding to
     * this parallel track group */
    min_track = max_track;
    max_track += _track_generator->getNumTracksByParallelGroup(i);

#pragma omp parallel
    {

      int tid = omp_get_thread_num();
      int azim_index, num_segments;
      Track* curr_track;
      segment* curr_segment;
      segment* segments;
      FP_PRECISION* track_flux;

      /* Use local array accumulator to prevent false sharing */
      FP_PRECISION thread_sr_flux[_num_groups];

      /* Loop over each thread within this azimuthal angle halfspace */
#pragma omp for schedule(guided)
      for (int track_id=min_track; track_id < max_track; track_id++) {

        /* Initialize local pointers to important data structures */
        curr_track = _tracks[track_id];
        azim_index = curr_track->getAzimAngleIndex();
        num_segments = curr_track->getNumSegments();
        segments = curr_track->getSegments();
        track_flux = &_boundary_flux(track_id,0,0,0);

        /* Loop over each Track segment in forward direction */
        for (int s=0; s < num_segments; s++) {
          curr_segment = &segments[s];
          tallyScalarFlux(curr_segment, azim_index, track_flux,
                          thread_sr_flux);
          tallyCurrent(curr_segment, azim_index, track_flux, true);
        }

        /* Transfer boundary angular flux to outgoing Track */
        transferBoundaryFlux(track_id, azim_index, true, track_flux);

        /* Loop over each Track segment in reverse direction */
        track_flux += _polar_times_groups;

        for (int s=num_segments-1; s > -1; s--) {
          curr_segment = &segments[s];
          tallyScalarFlux(curr_segment, azim_index, track_flux,
                          thread_sr_flux);
          tallyCurrent(curr_segment, azim_index, track_flux, false);
        }

        /* Transfer boundary angular flux to outgoing Track */
        transferBoundaryFlux(track_id, azim_index, false, track_flux);
      }
    }
  }

  return;
}


/**
 * @brief Computes the contribution to the SR scalar flux from a Track segment.
 * @details This method integrates the angular flux for a Track segment across
 *          energy groups and polar angles, and tallies it into the SR
 *          scalar flux, and updates the Track's angular flux.
 * @param curr_segment a pointer to the Track segment of interest
 * @param azim_index a pointer to the azimuthal angle index for this segment
 * @param track_flux a pointer to the Track's angular flux
 * @param sr_flux a pointer to the temporary SR flux buffer
 * @param fwd
 */
void CPUSolver::tallyScalarFlux(segment* curr_segment, int azim_index,
                                FP_PRECISION* track_flux,
                                FP_PRECISION* sr_flux) {

  int sr_id = curr_segment->_region_id;
  FP_PRECISION length = curr_segment->_length;
  FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();
  FP_PRECISION delta_psi, exponential;

  /* Set the SR scalar flux buffer to zero */
  memset(sr_flux, 0.0, _num_groups * sizeof(FP_PRECISION));

  /* Compute change in angular flux along segment in this SR */
  for (int e=0; e < _num_groups; e++) {
    for (int p=0; p < _num_polar; p++) {
      exponential = _exp_evaluator->computeExponential(sigma_t[e] * length, p);
      delta_psi = (track_flux(p,e)-_reduced_sources(sr_id,e)) * exponential;
      sr_flux[e] += delta_psi * _polar_weights(azim_index,p);
      track_flux(p,e) -= delta_psi;
    }
  }

  /* Atomically increment the SR scalar flux from the temporary array */
  omp_set_lock(&_SR_locks[sr_id]);
  {
    for (int e=0; e < _num_groups; e++)
      _scalar_flux(sr_id,e) += sr_flux[e];
  }
  omp_unset_lock(&_SR_locks[sr_id]);
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
                             FP_PRECISION* track_flux, bool fwd) {

  /* Tally surface currents if CMFD is in use */
  if (_cmfd != NULL && _cmfd->isFluxUpdateOn())
    _cmfd->tallyCurrent(curr_segment, track_flux,
                        &_polar_weights(azim_index,0), fwd);
}


/**
 * @brief Updates the boundary flux for a Track given boundary conditions.
 * @details For reflective and periodic boundary conditions, the outgoing
 *          boundary flux for the Track is given to the corresponding reflecting
 *          or periodic Track. For vacuum boundary conditions, the outgoing flux
 *          is tallied as leakage.
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
  bool transfer_flux;
  int track_out_id;

  /* For the "forward" direction */
  if (direction) {
    start = _tracks[track_id]->isNextOut() * _polar_times_groups;
    transfer_flux = _tracks[track_id]->getTransferFluxOut();
    track_out_id = _tracks[track_id]->getTrackOut()->getUid();
  }

  /* For the "reverse" direction */
  else {
    start = _tracks[track_id]->isNextIn() * _polar_times_groups;
    transfer_flux = _tracks[track_id]->getTransferFluxIn();
    track_out_id = _tracks[track_id]->getTrackIn()->getUid();
  }

  FP_PRECISION* track_out_flux = &_boundary_flux(track_out_id,0,0,start);

  /* Loop over polar angles and energy groups */
  for (int e=0; e < _num_groups; e++) {
    for (int p=0; p < _num_polar; p++)
      track_out_flux(p,e) = track_flux(p,e) * transfer_flux;
  }
}


/**
 * @brief Add the source term contribution in the transport equation to
 *        the SR scalar flux.
 */
void CPUSolver::addSourceToScalarFlux() {

  FP_PRECISION volume;
  FP_PRECISION* sigma_t;

  /* Add in source term and normalize flux to volume for each SR */
  /* Loop over SRs, energy groups */
#pragma omp parallel for private(volume, sigma_t) schedule(guided)
  for (int r=0; r < _num_SRs; r++) {
    volume = _SR_volumes[r];
    sigma_t = _SR_materials[r]->getSigmaT();

    for (int e=0; e < _num_groups; e++) {
      _scalar_flux(r,e) *= 0.5;
      _scalar_flux(r,e) /= (sigma_t[e] * volume);
      _scalar_flux(r,e) += (FOUR_PI * _reduced_sources(r,e));
    }
  }
}


/**
 * @brief Computes the volume-integrated, energy-integrated nu-fission rate in
 *        each SR and stores them in an array indexed by SR ID.
 * @details This is a helper method for SWIG to allow users to retrieve
 *          SR nu-fission rates as a NumPy array. An example of how this
 *          method can be called from Python is as follows:
 *
 * @code
 *          num_SRs = geometry.getNumSRs()
 *          fission_rates = solver.computeSRFissionRates(num_SRs)
 * @endcode
 *
 * @param fission_rates an array to store the nu-fission rates (implicitly
 *                      passed in as a NumPy array from Python)
 * @param num_SRs the number of SRs passed in from Python
 */
void CPUSolver::computeSRFissionRates(double* fission_rates, int num_SRs) {

  if (_scalar_flux == NULL)
    log_printf(ERROR, "Unable to compute SR fission rates since the "
               "source distribution has not been calculated");

  log_printf(INFO, "Computing SR fission rates...");

  FP_PRECISION* nu_sigma_f;
  FP_PRECISION volume;

  /* Initialize fission rates to zero */
  for (int r=0; r < _num_SRs; r++)
    fission_rates[r] = 0.0;

  /* Loop over all SRs and compute the volume-averaged nu-fission rate */
#pragma omp parallel for private (nu_sigma_f, volume) schedule(guided)
  for (int r=0; r < _num_SRs; r++) {
    nu_sigma_f = _SR_materials[r]->getNuSigmaF();
    volume = _SR_volumes[r];

    for (int e=0; e < _num_groups; e++)
      fission_rates[r] += nu_sigma_f[e] * _scalar_flux(r,e) * volume;
  }
}
