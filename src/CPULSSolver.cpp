#include "CPULSSolver.h"


/**
 * @brief Constructor initializes array pointers for Tracks and Materials.
 * @details The constructor retrieves the number of energy groups and FSRs
 *          and azimuthal angles from the Geometry and TrackGenerator if
 *          passed in as parameters by the user. The constructor initalizes
 *          the number of OpenMP threads to a default of 1.
 * @param track_generator an optional pointer to the TrackGenerator
 */
CPULSSolver::CPULSSolver(TrackGenerator* track_generator)
    : CPUSolver(track_generator) {
  _FSR_Cs = NULL;
  _FSR_Ms = NULL;

  _scalar_flux_x = NULL;
  _scalar_flux_y = NULL;

  _reduced_sources_x = NULL;
  _reduced_sources_y = NULL;

  _sin_phi = NULL;
  _cos_phi = NULL;
}


/**
 * @brief Destructor deletes array for OpenMP mutual exclusion locks for
 *        FSR scalar flux updates, and calls Solver parent class destructor
 *        to deletes arrays for fluxes and sources.
 */
CPULSSolver::~CPULSSolver() {

  if (_FSR_Cs != NULL)
    delete [] _FSR_Cs;

  if (_FSR_Ms != NULL)
    delete [] _FSR_Ms;

  if (_scalar_flux_x != NULL)
    delete [] _scalar_flux_x;

  if (_scalar_flux_y != NULL)
    delete [] _scalar_flux_y;

  if (_reduced_sources_x != NULL)
    delete [] _reduced_sources_x;

  if (_reduced_sources_y != NULL)
    delete [] _reduced_sources_y;

  if (_sin_phi != NULL)
    delete [] _sin_phi;

  if (_cos_phi != NULL)
    delete [] _cos_phi;
}


/**
 * @brief Allocates memory for Track boundary angular and FSR scalar fluxes.
 * @details Deletes memory for old flux arrays if they were allocated
 *          for a previous simulation.
 */
void CPULSSolver::initializeFluxArrays() {
  CPUSolver::initializeFluxArrays();

  /* Delete old flux arrays if they exist */
  if (_scalar_flux_x != NULL)
    delete [] _scalar_flux_x;

  if (_scalar_flux_y != NULL)
    delete [] _scalar_flux_y;

  try{
    /* Allocate an array for the FSR scalar flux */
    int size = _num_FSRs * _num_groups;
    _scalar_flux_x = new FP_PRECISION[size];
    _scalar_flux_y = new FP_PRECISION[size];
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the scalar flux moments");
  }
}


/**
 * @brief Allocates memory for FSR source arrays.
 * @details Deletes memory for old source arrays if they were allocated for a
 *          previous simulation.
 */
void CPULSSolver::initializeSourceArrays() {
  CPUSolver::initializeSourceArrays();

  /* Delete old sources arrays if they exist */
  if (_reduced_sources_x != NULL)
    delete [] _reduced_sources_x;

  if (_reduced_sources_y != NULL)
    delete [] _reduced_sources_y;

  int size = _num_FSRs * _num_groups;

  /* Allocate memory for all source arrays */
  try{
    _reduced_sources_x = new FP_PRECISION[size];
    _reduced_sources_y = new FP_PRECISION[size];
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for FSR source moments");
  }

  /* Initialize fixed sources to zero */
  memset(_reduced_sources_x, 0.0, sizeof(FP_PRECISION) * size);
  memset(_reduced_sources_y, 0.0, sizeof(FP_PRECISION) * size);

  if (_sin_phi != NULL)
    delete [] _sin_phi;

  if (_cos_phi != NULL)
    delete [] _cos_phi;

  int num_azim = _track_generator->getNumAzim();
  _sin_phi = new double[num_azim];
  _cos_phi = new double[num_azim];

  for (int i=0; i < num_azim; i++) {
    _sin_phi[i] = sin(_track_generator->getPhi(i));
    _cos_phi[i] = cos(_track_generator->getPhi(i));
  }

  if (_FSR_Cs != NULL)
    delete [] _FSR_Cs;

  if (_FSR_Ms != NULL)
    delete [] _FSR_Ms;

  _FSR_Ms = _track_generator->getFSRMs(_polar_quad);
  _FSR_Cs = _track_generator->getFSRCs(_polar_quad, _exp_evaluator);
}


/**
 * @brief Set the scalar flux for each FSR and energy group to some value.
 * @param value the value to assign to each FSR scalar flux
 */
void CPULSSolver::flattenFSRFluxes(FP_PRECISION value) {
  CPUSolver::flattenFSRFluxes(value);

#pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _num_groups; e++) {
      _scalar_flux_x(r,e) = 0.0;
      _scalar_flux_y(r,e) = 0.0;
    }
  }
}


/**
 * @brief Normalizes all FSR scalar fluxes and Track boundary angular
 *        fluxes to the total fission source (times \f$ \nu \f$).
 */
void CPULSSolver::normalizeFluxes() {

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

  for (int r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _num_groups; e++) {
      _scalar_flux(r,e) *= norm_factor;
      _scalar_flux_x(r,e) *= norm_factor;
      _scalar_flux_y(r,e) *= norm_factor;
      //log_printf(NORMAL, "(%d, %d): (%f, %f, %f)", r, e, _scalar_flux(r,e),
      //          _scalar_flux_x(r,e), _scalar_flux_y(r,e));
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
void CPULSSolver::computeFSRSources() {
  CPUSolver::computeFSRSources();

#pragma omp parallel
  {
    Material* material;
    FP_PRECISION* sigma_t;
    FP_PRECISION sigma_s, fiss_mat;
    FP_PRECISION scatter_source_x, fission_source_x;
    FP_PRECISION scatter_source_y, fission_source_y;
    FP_PRECISION* fission_sources_x = new FP_PRECISION[_num_groups];
    FP_PRECISION* scatter_sources_x = new FP_PRECISION[_num_groups];
    FP_PRECISION* fission_sources_y = new FP_PRECISION[_num_groups];
    FP_PRECISION* scatter_sources_y = new FP_PRECISION[_num_groups];
    FP_PRECISION det, src_x, src_y;

    /* Compute the total source for each FSR */
#pragma omp for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

      material = _FSR_materials[r];
      sigma_t = material->getSigmaT();
      det = 1.0 / (_FSR_Ms[r*3  ] * _FSR_Ms[r*3+1] -
                   _FSR_Ms[r*3+2] * _FSR_Ms[r*3+2]);

      /* Compute scatter + fission source for group g */
      for (int g=0; g < _num_groups; g++) {
        for (int g_prime=0; g_prime < _num_groups; g_prime++) {
          sigma_s = material->getSigmaSByGroup(g_prime+1,g+1);
          fiss_mat = material->getFissionMatrixByGroup(g_prime+1,g+1);
          scatter_sources_x[g_prime] = sigma_s * _scalar_flux_x(r,g_prime);
          fission_sources_x[g_prime] = fiss_mat * _scalar_flux_x(r,g_prime);
          scatter_sources_y[g_prime] = sigma_s * _scalar_flux_y(r,g_prime);
          fission_sources_y[g_prime] = fiss_mat * _scalar_flux_y(r,g_prime);
        }

        scatter_source_x = pairwise_sum<FP_PRECISION>(scatter_sources_x,
                                                      _num_groups);
        scatter_source_y = pairwise_sum<FP_PRECISION>(scatter_sources_y,
                                                      _num_groups);
        fission_source_x = pairwise_sum<FP_PRECISION>(fission_sources_x,
                                                      _num_groups);
        fission_source_y = pairwise_sum<FP_PRECISION>(fission_sources_y,
                                                      _num_groups);
        fission_source_x /= _k_eff;
        fission_source_y /= _k_eff;

        src_x = scatter_source_x + fission_source_x;
        src_y = scatter_source_y + fission_source_y;

        /* Compute total (scatter+fission+fixed) reduced source */
        _reduced_sources_x(r,g) = det * (_FSR_Ms[r*3+1] * src_x - _FSR_Ms[r*3+2] * src_y)
            * ONE_OVER_FOUR_PI / sigma_t[g];
        _reduced_sources_y(r,g) = det * (_FSR_Ms[r*3  ] * src_y - _FSR_Ms[r*3+2] * src_x)
            * ONE_OVER_FOUR_PI / sigma_t[g];
      }
    }

    delete [] fission_sources_x;
    delete [] scatter_sources_x;
    delete [] fission_sources_y;
    delete [] scatter_sources_y;
  }
}


/**
 * @brief This method performs one transport sweep of all azimuthal angles,
 *        Tracks, Track segments, polar angles and energy groups.
 * @details The method integrates the flux along each Track and updates the
 *          boundary fluxes for the corresponding output Track, while updating
 *          the scalar flux in each flat source region.
 */
void CPULSSolver::transportSweep() {

  log_printf(DEBUG, "Transport sweep with %d OpenMP threads", _num_threads);

  int min_track = 0;
  int max_track = 0;

  /* Initialize flux in each FSr to zero */
  flattenFSRFluxes(0.0);

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
      FP_PRECISION s_aki;
      double X, Y, x, y;
      Point* centroid;

      /* Use local array accumulator to prevent false sharing */
      FP_PRECISION thread_fsr_flux[_num_groups*3];

      /* Loop over each thread within this azimuthal angle halfspace */
#pragma omp for schedule(guided)
      for (int track_id=min_track; track_id < max_track; track_id++) {

        /* Initialize local pointers to important data structures */
        curr_track = _tracks[track_id];
        azim_index = curr_track->getAzimAngleIndex();
        num_segments = curr_track->getNumSegments();
        segments = curr_track->getSegments();
        track_flux = &_boundary_flux(track_id,0,0,0);
        X = curr_track->getStart()->getX();
        Y = curr_track->getStart()->getY();

        /* Loop over each Track segment in forward direction */
        for (int s=0; s < num_segments; s++) {
          curr_segment = &segments[s];
          s_aki = curr_segment->_length;
          centroid = _geometry->getFSRCentroid(curr_segment->_region_id);
          x = X - centroid->getX();
          y = Y - centroid->getY();

          tallyLSScalarFlux(curr_segment, azim_index, track_flux,
                            thread_fsr_flux, x, y, true);
          tallyCurrent(curr_segment, azim_index, track_flux, true);

          X += s_aki * _cos_phi[azim_index];
          Y += s_aki * _sin_phi[azim_index];
        }

        /* Transfer boundary angular flux to outgoing Track */
        transferBoundaryFlux(track_id, azim_index, true, track_flux);

        /* Loop over each Track segment in reverse direction */
        track_flux += _polar_times_groups;

        for (int s=num_segments-1; s > -1; s--) {
          curr_segment = &segments[s];
          s_aki = curr_segment->_length;
          centroid = _geometry->getFSRCentroid(curr_segment->_region_id);
          x = X - centroid->getX();
          y = Y - centroid->getY();

          tallyLSScalarFlux(curr_segment, azim_index, track_flux,
                            thread_fsr_flux, x, y, false);
          tallyCurrent(curr_segment, azim_index, track_flux, false);

          X -= s_aki * _cos_phi[azim_index];
          Y -= s_aki * _sin_phi[azim_index];
        }

        /* Transfer boundary angular flux to outgoing Track */
        transferBoundaryFlux(track_id, azim_index, false, track_flux);
      }
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
void CPULSSolver::tallyLSScalarFlux(segment* curr_segment, int azim_index,
                                    FP_PRECISION* track_flux,
                                    FP_PRECISION* fsr_flux,
                                    double x, double y, bool fwd) {

  int fsr_id = curr_segment->_region_id;
  FP_PRECISION s_aki = curr_segment->_length;
  FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();
  FP_PRECISION delta_psi, exp_F1, exp_F2, exp_G1, tau_aki, tau_mki;
  FP_PRECISION src_constant, src_moment;
  FP_PRECISION ax, ay;
  FP_PRECISION* sin_thetas = _polar_quad->getSinThetas();
  double cos_phi, sin_phi;

  if (fwd) {
    cos_phi = _cos_phi[azim_index];
    sin_phi = _sin_phi[azim_index];
  }
  else {
    cos_phi = - _cos_phi[azim_index];
    sin_phi = - _sin_phi[azim_index];
  }

  double xc_aki = x + s_aki * cos_phi / 2.0;
  double yc_aki = y + s_aki * sin_phi / 2.0;

  /* Set the FSR scalar flux buffer to zero */
  memset(fsr_flux, 0.0, _num_groups * 3 * sizeof(FP_PRECISION));

  /* Compute change in angular flux along segment in this FSR */
  for (int e=0; e < _num_groups; e++) {

    tau_aki = sigma_t[e] * s_aki;
    src_constant = _reduced_sources(fsr_id, e) +
        _reduced_sources_x(fsr_id, e) * xc_aki +
        _reduced_sources_y(fsr_id, e) * yc_aki;

    for (int p=0; p < _num_polar; p++) {
      exp_F1 = _exp_evaluator->computeExponential  (tau_aki, p);
      exp_F2 = _exp_evaluator->computeExponentialF2(tau_aki, p);
      exp_G1 = _exp_evaluator->computeExponentialG1(tau_aki, p);
      ax = cos_phi * sin_thetas[p];
      ay = sin_phi * sin_thetas[p];
      tau_mki = tau_aki / sin_thetas[p];

      src_moment = (ax * _reduced_sources_x(fsr_id, e) +
                    ay * _reduced_sources_y(fsr_id, e)) /
          (2.0 * sigma_t[e]);

      delta_psi = (track_flux(p,e) - src_constant) * exp_F1 -
          src_moment * exp_F2;

      fsr_flux[e*3    ] += _polar_weights(azim_index,p) * delta_psi;
      fsr_flux[e*3 + 1] += _polar_weights(azim_index,p) *
          (cos_phi * s_aki * track_flux(p,e) *
           (tau_mki / 2.0 - exp_G1) + x * delta_psi);
      fsr_flux[e*3 + 2] += _polar_weights(azim_index,p) *
          (sin_phi * s_aki * track_flux(p,e) *
           (tau_mki / 2.0 - exp_G1) + y * delta_psi);

      track_flux(p,e) -= delta_psi;

      /*
      if (track_flux(p,e) < 0.0) {
        log_printf(NORMAL, "encountered track with negative flux, (%f, %f)",
                   _geometry->getFSRCentroid(fsr_id)->getX(),
                   _geometry->getFSRCentroid(fsr_id)->getY());
      }
      */

      track_flux(p,e) = std::max(track_flux(p,e), 0.0);
    }
  }

  /* Atomically increment the FSR scalar flux from the temporary array */
  omp_set_lock(&_FSR_locks[fsr_id]);

  for (int e=0; e < _num_groups; e++) {
    _scalar_flux  (fsr_id,e) += fsr_flux[e*3    ];
    _scalar_flux_x(fsr_id,e) += fsr_flux[e*3 + 1];
    _scalar_flux_y(fsr_id,e) += fsr_flux[e*3 + 2];
  }

  omp_unset_lock(&_FSR_locks[fsr_id]);
}


/**
 * @brief Add the source term contribution in the transport equation to
 *        the FSR scalar flux.
 */
void CPULSSolver::addSourceToScalarFlux() {

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

      if (_scalar_flux(r,e) < 0.0)
        log_printf(NORMAL, "encountered negative scalar flux (%f, %f)",
                   _geometry->getFSRCentroid(r)->getX(),
                   _geometry->getFSRCentroid(r)->getY());

      _scalar_flux_x(r,e) *= 0.5;
      _scalar_flux_x(r,e) /= (sigma_t[e] * volume);
      _scalar_flux_x(r,e) += (FOUR_PI * _reduced_sources_x(r,e) * _FSR_Cs[r*_num_groups*3 + 3*e]);
      _scalar_flux_x(r,e) += (FOUR_PI * _reduced_sources_y(r,e) * _FSR_Cs[r*_num_groups*3 + 3*e + 2]);

      _scalar_flux_y(r,e) *= 0.5;
      _scalar_flux_y(r,e) /= (sigma_t[e] * volume);
      _scalar_flux_y(r,e) += (FOUR_PI * _reduced_sources_y(r,e) * _FSR_Cs[r*_num_groups*3 + 3*e + 1]);
      _scalar_flux_y(r,e) += (FOUR_PI * _reduced_sources_x(r,e) * _FSR_Cs[r*_num_groups*3 + 3*e + 2]);
    }
  }
}


/**
 * @brief Initializes new ExpEvaluator object to compute exponentials.
 */
void CPULSSolver::initializeExpEvaluator() {
  _exp_evaluator->useIntrinsic();
  Solver::initializeExpEvaluator();
}
