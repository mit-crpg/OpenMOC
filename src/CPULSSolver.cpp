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

  _FSR_source_constants = NULL;
  _FSR_lin_exp_matrix = NULL;
  _scalar_flux_xy = NULL;
  _reduced_sources_xy = NULL;
  _sin_phi = NULL;
  _cos_phi = NULL;
}


/**
 * @brief Destructor deletes array for OpenMP mutual exclusion locks for
 *        FSR scalar flux updates, and calls Solver parent class destructor
 *        to deletes arrays for fluxes and sources.
 */
CPULSSolver::~CPULSSolver() {

  if (_FSR_source_constants != NULL)
    delete [] _FSR_source_constants;

  if (_FSR_lin_exp_matrix != NULL)
    delete [] _FSR_lin_exp_matrix;

  if (_scalar_flux_xy != NULL)
    delete [] _scalar_flux_xy;

  if (_reduced_sources_xy != NULL)
    delete [] _reduced_sources_xy;

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

  /* Delete old flux moment arrays if they exist */
  if (_scalar_flux_xy != NULL)
    delete [] _scalar_flux_xy;

  try{
    /* Allocate an array for the FSR scalar flux */
    int size = _num_FSRs * _num_groups * 2;
    _scalar_flux_xy = new FP_PRECISION[size];
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

  /* Delete old sources moment arrays if they exist */
  if (_reduced_sources_xy != NULL)
    delete [] _reduced_sources_xy;

  int size = _num_FSRs * _num_groups * 2;

  /* Allocate memory for all source arrays */
  try {
    _reduced_sources_xy = new FP_PRECISION[size];
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for FSR source moments");
  }

  /* Initialize source moments to zero */
  memset(_reduced_sources_xy, 0.0, sizeof(FP_PRECISION) * size);

  /* Delete old sin(phi) and cos(phi) arrays if they exist */
  if (_sin_phi != NULL)
    delete [] _sin_phi;

  if (_cos_phi != NULL)
    delete [] _cos_phi;

  /* Get the sin theta terms from the polar quad */
  _sin_thetas = _polar_quad->getSinThetas();
  _inv_sin_thetas = _polar_quad->getInverseSinThetas();

  /* Generate the sin(phi) and cos(phi) arrays */
  int num_azim = _track_generator->getNumAzim();
  _sin_phi = new FP_PRECISION[num_azim];
  _cos_phi = new FP_PRECISION[num_azim];

  for (int i=0; i < num_azim; i++) {
    _sin_phi[i] = FP_PRECISION(sin(_track_generator->getPhi(i)));
    _cos_phi[i] = FP_PRECISION(cos(_track_generator->getPhi(i)));
  }

  /* Delete old FSR source constant and lin exp matrix arrays if they exist */
  if (_FSR_source_constants != NULL)
    delete [] _FSR_source_constants;

  if (_FSR_lin_exp_matrix != NULL)
    delete [] _FSR_lin_exp_matrix;

  /* Get the FSR lin exp matrix and source constants arrays */
  _FSR_lin_exp_matrix =
      _track_generator->getFSRLinearExpansionCoeffs(_polar_quad);
  _FSR_source_constants =
      _track_generator->getFSRSourceConstants(_polar_quad, _exp_evaluator);
}


/**
 * @brief Set the scalar flux constants for each FSR and energy group to some
 *        value and the scalar flux moments to zero.
 * @param value the value to assign to each FSR scalar flux
 */
void CPULSSolver::flattenFSRFluxes(FP_PRECISION value) {
  CPUSolver::flattenFSRFluxes(value);

#pragma omp parallel for schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _num_groups; e++) {
      _scalar_flux_xy(r,e,0) = 0.0;
      _scalar_flux_xy(r,e,1) = 0.0;
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
      _scalar_flux_xy(r,e,0) *= norm_factor;
      _scalar_flux_xy(r,e,1) *= norm_factor;
      _old_scalar_flux(r,e) *= norm_factor;
    }
  }

  /* Normalize angular boundary fluxes for each Track */
#pragma omp parallel for schedule(guided)
  for (int t=0; t < _tot_num_tracks; t++) {
    for (int d=0; d < 2; d++) {
      for (int p=0; p < _num_polar; p++) {
        for (int e=0; e < _num_groups; e++)
          _boundary_flux(t,d,p,e) *= norm_factor;
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
    FP_PRECISION scatter_source_x, scatter_source_y;
    FP_PRECISION fission_source_x, fission_source_y;
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

      /* Compute the determinat of the linear expansion coefficient matrix
       * for this FSR */
      det = _FSR_lin_exp_matrix[r*3  ] * _FSR_lin_exp_matrix[r*3+1] -
          _FSR_lin_exp_matrix[r*3+2] * _FSR_lin_exp_matrix[r*3+2];

      /* Compute scatter + fission source for group g */
      for (int g=0; g < _num_groups; g++) {
        for (int g_prime=0; g_prime < _num_groups; g_prime++) {
          sigma_s = material->getSigmaSByGroup(g_prime+1,g+1);
          fiss_mat = material->getFissionMatrixByGroup(g_prime+1,g+1);
          scatter_sources_x[g_prime] = sigma_s * _scalar_flux_xy(r,g_prime,0);
          fission_sources_x[g_prime] = fiss_mat * _scalar_flux_xy(r,g_prime,0);
          scatter_sources_y[g_prime] = sigma_s * _scalar_flux_xy(r,g_prime,1);
          fission_sources_y[g_prime] = fiss_mat * _scalar_flux_xy(r,g_prime,1);
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

        /* Compute total (scatter+fission+fixed) reduced source moments */
        _reduced_sources_xy(r,g,0) = ONE_OVER_FOUR_PI / (det * sigma_t[g]) *
            (_FSR_lin_exp_matrix[r*3+1] * src_x -
             _FSR_lin_exp_matrix[r*3+2] * src_y);
        _reduced_sources_xy(r,g,1) = ONE_OVER_FOUR_PI / (det * sigma_t[g]) *
            (_FSR_lin_exp_matrix[r*3  ] * src_y -
             _FSR_lin_exp_matrix[r*3+2] * src_x);
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

  /* Initialize flux in each FSR to zero */
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

      int azim_index, num_segments;
      Track* curr_track;
      segment* curr_segment;
      segment* segments;
      FP_PRECISION* track_flux;
      FP_PRECISION length;
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

        /* Get the starting point for the first sesgment in the global
         * coordinate system */
        X = curr_track->getStart()->getX();
        Y = curr_track->getStart()->getY();

        /* Loop over each Track segment in forward direction */
        for (int s=0; s < num_segments; s++) {
          curr_segment = &segments[s];
          length = curr_segment->_length;
          centroid = _FSR_centroids[curr_segment->_region_id];

          /* Get the starting point of the segment in local coordinates */
          x = X - centroid->getX();
          y = Y - centroid->getY();

          tallyLSScalarFlux(curr_segment, azim_index, track_flux,
                            thread_fsr_flux, x, y, 1);
          tallyCurrent(curr_segment, azim_index, track_flux, true);

          /* Increment the segment starting point to the next segment */
          X += length * _cos_phi[azim_index];
          Y += length * _sin_phi[azim_index];
        }

        /* Transfer boundary angular flux to outgoing Track */
        transferBoundaryFlux(track_id, azim_index, true, track_flux);

        /* Loop over each Track segment in reverse direction */
        track_flux += _polar_times_groups;

        for (int s=num_segments-1; s > -1; s--) {
          curr_segment = &segments[s];
          length = curr_segment->_length;
          centroid = _FSR_centroids[curr_segment->_region_id];

          /* Get the starting point of the segment in local coordinates */
          x = X - centroid->getX();
          y = Y - centroid->getY();

          tallyLSScalarFlux(curr_segment, azim_index, track_flux,
                            thread_fsr_flux, x, y, -1);
          tallyCurrent(curr_segment, azim_index, track_flux, false);

          X -= length * _cos_phi[azim_index];
          Y -= length * _sin_phi[azim_index];
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
 * @param x the x-coord of the segment starting point
 * @param y the y-coord of the segment starting point
 * @param fwd int indicating whether the segment is pointing forward (1) or
 *            backwards (-1)
 */
void CPULSSolver::tallyLSScalarFlux(segment* curr_segment, int azim_index,
                                    FP_PRECISION* track_flux,
                                    FP_PRECISION* fsr_flux,
                                    double x, double y, int fwd) {

  int fsr_id = curr_segment->_region_id;
  FP_PRECISION length = curr_segment->_length;
  FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();
  FP_PRECISION delta_psi, exp_F1, exp_F2, exp_H, tau;
  FP_PRECISION src_flat, src_linear, dt, dt2, ax, ay;
  FP_PRECISION polar_wgt_d_psi;
  FP_PRECISION polar_wgt_exp_h;
  int exp_index;
  int azim_times_polar = azim_index * _num_polar;

  /* Compute the segment midpoint */
  double cos_phi = fwd * _cos_phi[azim_index];
  double sin_phi = fwd * _sin_phi[azim_index];
  double xc = x + 0.5 * length * cos_phi;
  double yc = y + 0.5 * length * sin_phi;

  /* Set the FSR scalar flux buffer to zero */
  memset(fsr_flux, 0.0, _num_groups * 3 * sizeof(FP_PRECISION));

  /* Compute change in angular flux along segment in this LSR */
  for (int e=0; e < _num_groups; e++) {

    /* Compute the optical length */
    tau = sigma_t[e] * length;

    /* Get the location of the optical length in the exp look-up table */
    exp_index = floor(tau * _inv_spacing);

    /* Compute the distance and distance squared from the optical length to the
     * nearest tau in the table */
    dt = tau - exp_index * _spacing;
    dt2 = dt * dt;

    /* Increment the exp index to account for the staggered values in the
     * exp table */
    exp_index *= 9 * _num_polar;

    polar_wgt_d_psi = 0.0;
    polar_wgt_exp_h = 0.0;

    /* Compute the flat component of the reduced source */
    src_flat = _reduced_sources(fsr_id, e) +
        _reduced_sources_xy(fsr_id, e, 0) * xc +
        _reduced_sources_xy(fsr_id, e, 1) * yc;

    for (int p=0; p < _num_polar; p++) {

      /* Compute the exponential terms */
      exp_F1 = _exp_evaluator->computeExponentialInline(exp_index    , p, dt, dt2);
      exp_F2 = _exp_evaluator->computeExponentialInline(exp_index + 3, p, dt, dt2);
      exp_H  = _exp_evaluator->computeExponentialInline(exp_index + 6, p, dt, dt2) * length * track_flux(p,e);

      /* Increment the exp index for the next polar angle */
      exp_index += 9;

      ax = cos_phi * _sin_thetas[p];
      ay = sin_phi * _sin_thetas[p];

      /* Compute the moment component of the source */
      src_linear = (ax * _reduced_sources_xy(fsr_id, e, 0) +
                    ay * _reduced_sources_xy(fsr_id, e, 1)) * exp_F2;

      /* Compute the change in flux across the segment */
      delta_psi = (track_flux(p,e) - src_flat) * exp_F1 - src_linear * exp_F2;

      polar_wgt_d_psi += _polar_weights[azim_times_polar + p] * delta_psi;
      polar_wgt_exp_h += _polar_weights[azim_times_polar + p] * exp_H;

      /* Decrement the track flux */
      track_flux(p,e) -= delta_psi;
    }

    /* Increment the fsr scalar flux and scalar flux moments */
    fsr_flux[e*3    ] += polar_wgt_d_psi;
    fsr_flux[e*3 + 1] += polar_wgt_exp_h * cos_phi + polar_wgt_d_psi * x;
    fsr_flux[e*3 + 2] += polar_wgt_exp_h * sin_phi + polar_wgt_d_psi * y;
  }

  /* Atomically increment the FSR scalar flux from the temporary array */
  omp_set_lock(&_FSR_locks[fsr_id]);

  for (int e=0; e < _num_groups; e++) {
    _scalar_flux   (fsr_id,e)   += fsr_flux[e*3    ];
    _scalar_flux_xy(fsr_id,e,0) += fsr_flux[e*3 + 1];
    _scalar_flux_xy(fsr_id,e,1) += fsr_flux[e*3 + 2];
  }

  omp_unset_lock(&_FSR_locks[fsr_id]);
}


/**
 * @brief Add the source term contribution in the transport equation to
 *        the FSR scalar flux.
 */
void CPULSSolver::addSourceToScalarFlux() {

#pragma omp parallel
  {
    FP_PRECISION volume;
    FP_PRECISION* sigma_t;

    /* Add in source term and normalize flux to volume for each FSR */
    /* Loop over FSRs, energy groups */
#pragma omp for
    for (int r=0; r < _num_FSRs; r++) {
      volume = _FSR_volumes[r];
      sigma_t = _FSR_materials[r]->getSigmaT();

      for (int e=0; e < _num_groups; e++) {

        _scalar_flux(r,e) *= 0.5;
        _scalar_flux(r,e) /= (sigma_t[e] * volume);
        _scalar_flux(r,e) += (FOUR_PI * _reduced_sources(r,e));

        _scalar_flux_xy(r,e,0) *= 0.5;
        _scalar_flux_xy(r,e,0) /= (sigma_t[e] * volume);
        _scalar_flux_xy(r,e,0) += FOUR_PI * _reduced_sources_xy(r,e,0)
            * _FSR_source_constants[r*_num_groups*3 + 3*e];
        _scalar_flux_xy(r,e,0) += FOUR_PI * _reduced_sources_xy(r,e,1)
            * _FSR_source_constants[r*_num_groups*3 + 3*e + 2];

        _scalar_flux_xy(r,e,1) *= 0.5;
        _scalar_flux_xy(r,e,1) /= (sigma_t[e] * volume);
        _scalar_flux_xy(r,e,1) += FOUR_PI * _reduced_sources_xy(r,e,1)
            * _FSR_source_constants[r*_num_groups*3 + 3*e + 1];
        _scalar_flux_xy(r,e,1) += FOUR_PI * _reduced_sources_xy(r,e,0)
            * _FSR_source_constants[r*_num_groups*3 + 3*e + 2];
      }
    }
  }
}


/**
 * @brief Get the flux at a specific point in the geometry.
 * @param coords The coords of the point to get the flux at
 * @param group the energy group
 */
FP_PRECISION CPULSSolver::getFluxByCoords(LocalCoords* coords, int group) {

  double x, y, xc, yc, det;

  coords->setUniverse(_geometry->getRootUniverse());
  Cell* cell = _geometry->findCellContainingCoords(coords);
  int fsr_id = _geometry->getFSRId(coords);
  Point* centroid = _geometry->getFSRCentroid(fsr_id);
  x = coords->getX();
  y = coords->getY();
  xc = centroid->getX();
  yc = centroid->getY();

  det = _FSR_lin_exp_matrix[fsr_id*3] * _FSR_lin_exp_matrix[fsr_id*3+1] -
      _FSR_lin_exp_matrix[fsr_id*3+2] * _FSR_lin_exp_matrix[fsr_id*3+2];

  FP_PRECISION flux = _scalar_flux(fsr_id, group);
  flux += 1.0 / det *
      (_FSR_lin_exp_matrix[fsr_id*3+1] * _scalar_flux_xy(fsr_id, group, 0) -
       _FSR_lin_exp_matrix[fsr_id*3+2] * _scalar_flux_xy(fsr_id, group, 1))
      * (x - xc);
  flux += 1.0 / det *
      (_FSR_lin_exp_matrix[fsr_id*3  ] * _scalar_flux_xy(fsr_id, group, 1) -
       _FSR_lin_exp_matrix[fsr_id*3+2] * _scalar_flux_xy(fsr_id, group, 0))
      * (y - yc);

  return flux;
}


void CPULSSolver::initializeCmfd() {
  Solver::initializeCmfd();
  if (_cmfd != NULL) {
    _cmfd->setLinearSourceOn(true);
    _cmfd->setLSRFluxMoments(_scalar_flux_xy);
  }
}


void CPULSSolver::initializeExpEvaluator() {
  _exp_evaluator->useLinearSource();
  Solver::initializeExpEvaluator();

  _inv_spacing = _exp_evaluator->getInverseTableSpacing();
  _spacing = _exp_evaluator->getTableSpacing();
}
