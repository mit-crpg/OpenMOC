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
  _scalar_flux_xyz = NULL;
  _reduced_sources_xyz = NULL;
  _stabalizing_flux_xyz = NULL;
  _stabalize_moments = true;
  _source_type = "Linear";
}


/**
 * @brief Destructor deletes array for OpenMP mutual exclusion locks for
 *        FSR scalar flux updates, and calls Solver parent class destructor
 *        to deletes arrays for fluxes and sources.
 */
CPULSSolver::~CPULSSolver() {

  if (_scalar_flux_xyz != NULL)
    delete [] _scalar_flux_xyz;

  if (_reduced_sources_xyz != NULL)
    delete [] _reduced_sources_xyz;
  
  if (_stabalizing_flux_xyz != NULL)
    delete [] _stabalizing_flux_xyz;

  if (_FSR_lin_exp_matrix != NULL)
    delete [] _FSR_lin_exp_matrix;

  if (_FSR_source_constants != NULL)
    delete [] _FSR_source_constants;
}


/**
 * @brief Allocates memory for Track boundary angular and FSR scalar fluxes.
 * @details Deletes memory for old flux arrays if they were allocated
 *          for a previous simulation.
 */
void CPULSSolver::initializeFluxArrays() {
  CPUSolver::initializeFluxArrays();

  /* Delete old flux moment arrays if they exist */
  if (_scalar_flux_xyz != NULL)
    delete [] _scalar_flux_xyz;

  try {
    /* Allocate an array for the FSR scalar flux */
    long size = _num_FSRs * _num_groups * 3;
    long max_size = size;
#ifdef MPIX
    if (_geometry->isDomainDecomposed())
      MPI_Allreduce(&size, &max_size, 1, MPI_LONG, MPI_MAX,
                    _geometry->getMPICart());
#endif
    double max_size_mb = (double) (max_size * sizeof(FP_PRECISION))
        / (double) (1e6);

    if (_stabalize_transport && _stabalize_moments)
      max_size_mb *= 2;

    log_printf(NORMAL, "Max linear flux storage per domain = %6.2f MB",
               max_size_mb);

    _scalar_flux_xyz = new FP_PRECISION[size];
    memset(_scalar_flux_xyz, 0., size * sizeof(FP_PRECISION));
    
    if (_stabalize_transport && _stabalize_moments) {
      _stabalizing_flux_xyz = new FP_PRECISION[size];
      memset(_stabalizing_flux_xyz, 0., size * sizeof(FP_PRECISION));
    }
  }
  catch (std::exception &e) {
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
  if (_reduced_sources_xyz != NULL)
    delete [] _reduced_sources_xyz;

  long size = _num_FSRs * _num_groups * 3;

  /* Allocate memory for all source arrays */
  try {
    long max_size = size;
#ifdef MPIX
    if (_geometry->isDomainDecomposed())
      MPI_Allreduce(&size, &max_size, 1, MPI_LONG, MPI_MAX,
                    _geometry->getMPICart());
#endif
    double max_size_mb = (double) (max_size * sizeof(FP_PRECISION))
        / (double) (1e6);
    log_printf(NORMAL, "Max linear source storage per domain = %6.2f MB",
               max_size_mb);

    _reduced_sources_xyz = new FP_PRECISION[size];
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for FSR source moments");
  }

  /* Initialize source moments to zero */
  memset(_reduced_sources_xyz, 0.0, sizeof(FP_PRECISION) * size);
}


/**
 * @brief Initializes the FSR volumes and Materials array.
 * @details This method allocates and initializes an array of OpenMP
 *          mutual exclusion locks for each FSR for use in the
 *          transport sweep algorithm.
 */
void CPULSSolver::initializeFSRs() {

  CPUSolver::initializeFSRs();

  /* Generate linear source coefficients */
  log_printf(NORMAL, "Generating linear expansion coefficients");
  LinearExpansionGenerator lin_src_coeffs(this);
  lin_src_coeffs.execute();
  log_printf(NORMAL, "Linear expansion coefficient generation complete");
}


/**
 * @brief Set the scalar flux constants for each FSR and energy group to some
 *        value and the scalar flux moments to zero.
 * @param value the value to assign to each FSR scalar flux
 */
void CPULSSolver::flattenFSRFluxes(FP_PRECISION value) {
  CPUSolver::flattenFSRFluxes(value);

#pragma omp parallel for schedule(guided)
  for (long r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _num_groups; e++) {
      _scalar_flux_xyz(r,e,0) = 0.0;
      _scalar_flux_xyz(r,e,1) = 0.0;
      _scalar_flux_xyz(r,e,2) = 0.0;
    }
  }
}


/**
 * @brief Normalizes all FSR scalar fluxes and Track boundary angular
 *        fluxes to the total fission source (times \f$ \nu \f$).
 */
double CPULSSolver::normalizeFluxes() {

  /* Normalize scalar fluxes in each FSR */
  double norm_factor = CPUSolver::normalizeFluxes();

#pragma omp parallel for schedule(guided)
  for (long r=0; r < _num_FSRs; r++) {
    for (int e=0; e < _num_groups; e++) {
      _scalar_flux_xyz(r,e,0) *= norm_factor;
      _scalar_flux_xyz(r,e,1) *= norm_factor;
      _scalar_flux_xyz(r,e,2) *= norm_factor;
    }
  }

  return norm_factor;
}


//FIXME
void CPULSSolver::checkLimitXS(int iteration) {

  Solver::checkLimitXS(iteration);

  if (iteration != _reset_iteration)
    return;

  /* Generate linear source coefficients */
  log_printf(NORMAL, "Generating linear expansion coefficients");
  LinearExpansionGenerator lin_src_coeffs(this);
  lin_src_coeffs.execute();
  log_printf(NORMAL, "Linear expansion coefficient generation complete");
}


/**
 * @brief Computes the total source (fission, scattering, fixed) in each FSR.
 * @details This method computes the total source in each FSR based on
 *          this iteration's current approximation to the scalar flux.
 */
void CPULSSolver::computeFSRSources(int iteration) {
  CPUSolver::computeFSRSources(iteration);

  int num_coeffs = 3;
  if (_solve_3D)
    num_coeffs = 6;

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    Material* material;
    FP_PRECISION* sigma_t;
    FP_PRECISION* nu_sigma_f;
    FP_PRECISION* chi;
    FP_PRECISION* sigma_s;
    FP_PRECISION src_x, src_y, src_z;

    /* Compute the total source for each FSR */
#pragma omp for schedule(guided)
    for (long r=0; r < _num_FSRs; r++) {

      material = _FSR_materials[r];
      nu_sigma_f = material->getNuSigmaF();
      chi = material->getChi();
      sigma_t = material->getSigmaT();
      sigma_s = material->getSigmaS();

      /* Compute scatter + fission source for group g */
      for (int g=0; g < _num_groups; g++) {

        int first_scattering_index = g * _num_groups;

        /* Compute scatter sources */
        FP_PRECISION* scatter_sources_x = _groupwise_scratch.at(tid);
        for (int g_prime=0; g_prime < _num_groups; g_prime++) {
          int idx = first_scattering_index + g_prime;
          if (sigma_s[idx] > 0) {
            scatter_sources_x[g_prime] = sigma_s[idx] *
              _scalar_flux_xyz(r,g_prime,0);
          }
          else {
            scatter_sources_x[g_prime] = 1.0e-20;
          }
        }
        double scatter_source_x =
            pairwise_sum<FP_PRECISION>(scatter_sources_x, _num_groups);

        FP_PRECISION* scatter_sources_y = _groupwise_scratch.at(tid);
        for (int g_prime=0; g_prime < _num_groups; g_prime++) {
          int idx = first_scattering_index + g_prime;
          if (sigma_s[idx] > 0) {
            scatter_sources_y[g_prime] = sigma_s[idx]
              * _scalar_flux_xyz(r,g_prime,1);
          }
          else {
            scatter_sources_y[g_prime] = 1.0e-20;
          }
        }
        double scatter_source_y =
            pairwise_sum<FP_PRECISION>(scatter_sources_y, _num_groups);

        FP_PRECISION* scatter_sources_z = _groupwise_scratch.at(tid);
        for (int g_prime=0; g_prime < _num_groups; g_prime++) {
          int idx = first_scattering_index + g_prime;
          if (sigma_s[idx] > 0.0) {
            scatter_sources_z[g_prime] = sigma_s[idx]
              * _scalar_flux_xyz(r,g_prime,2);
          }
          else {
            scatter_sources_y[g_prime] = 1.0e-20;
          }
        }
        double scatter_source_z =
            pairwise_sum<FP_PRECISION>(scatter_sources_z, _num_groups);


        /* Compute fission sources */
        FP_PRECISION* fission_sources_x = _groupwise_scratch.at(tid);
        for (int g_prime=0; g_prime < _num_groups; g_prime++)
          fission_sources_x[g_prime] =
              nu_sigma_f[g_prime] * _scalar_flux_xyz(r,g_prime,0);
        double fission_source_x =
            pairwise_sum<FP_PRECISION>(fission_sources_x, _num_groups);

        FP_PRECISION* fission_sources_y = _groupwise_scratch.at(tid);
        for (int g_prime=0; g_prime < _num_groups; g_prime++)
          fission_sources_y[g_prime] =
              nu_sigma_f[g_prime] * _scalar_flux_xyz(r,g_prime,1);
        double fission_source_y =
            pairwise_sum<FP_PRECISION>(fission_sources_y, _num_groups);

        FP_PRECISION* fission_sources_z = _groupwise_scratch.at(tid);
        for (int g_prime=0; g_prime < _num_groups; g_prime++)
          fission_sources_z[g_prime] =
              nu_sigma_f[g_prime] * _scalar_flux_xyz(r,g_prime,2);
        double fission_source_z =
            pairwise_sum<FP_PRECISION>(fission_sources_z, _num_groups);

        fission_source_x /= _k_eff;
        fission_source_y /= _k_eff;
        fission_source_z /= _k_eff;

        src_x = scatter_source_x + chi[g] * fission_source_x;
        src_y = scatter_source_y + chi[g] * fission_source_y;
        src_z = scatter_source_z + chi[g] * fission_source_z;

        /* Compute total (scatter+fission+fixed) reduced source moments */
        if (_solve_3D) {
          _reduced_sources_xyz(r,g,0) = ONE_OVER_FOUR_PI / 2 *
              (_FSR_lin_exp_matrix[r*num_coeffs  ] * src_x +
               _FSR_lin_exp_matrix[r*num_coeffs+2] * src_y +
               _FSR_lin_exp_matrix[r*num_coeffs+3] * src_z);
          _reduced_sources_xyz(r,g,1) = ONE_OVER_FOUR_PI / 2 *
              (_FSR_lin_exp_matrix[r*num_coeffs+2] * src_x +
               _FSR_lin_exp_matrix[r*num_coeffs+1] * src_y +
               _FSR_lin_exp_matrix[r*num_coeffs+4] * src_z);
          _reduced_sources_xyz(r,g,2) = ONE_OVER_FOUR_PI / 2 *
              (_FSR_lin_exp_matrix[r*num_coeffs+3] * src_x +
               _FSR_lin_exp_matrix[r*num_coeffs+4] * src_y +
               _FSR_lin_exp_matrix[r*num_coeffs+5] * src_z);
        }
        else {
          _reduced_sources_xyz(r,g,0) = ONE_OVER_FOUR_PI / 2 *
              (_FSR_lin_exp_matrix[r*num_coeffs  ] * src_x +
               _FSR_lin_exp_matrix[r*num_coeffs+2] * src_y);
          _reduced_sources_xyz(r,g,1) = ONE_OVER_FOUR_PI / 2 *
              (_FSR_lin_exp_matrix[r*num_coeffs+2] * src_x +
               _FSR_lin_exp_matrix[r*num_coeffs+1] * src_y);
        }
      }
    }
  }
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
 * @param z the z-coord of the segment starting point
 * @param fwd int indicating whether the segment is pointing forward (1) or
 *            backwards (-1)
//FIXME: NOT THE CORRECT PARAMS
 */
void CPULSSolver::tallyLSScalarFlux(segment* curr_segment, int azim_index,
                                    int polar_index,
                                    float* track_flux,
                                    FP_PRECISION* fsr_flux,
                                    FP_PRECISION* scratch_pad,
                                    FP_PRECISION direction[3]) {

  long fsr_id = curr_segment->_region_id;
  FP_PRECISION length = curr_segment->_length;
  FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();
  FP_PRECISION* position = curr_segment->_starting_position;
  ExpEvaluator* exp_evaluator = _exp_evaluators[azim_index][polar_index];

  long start_scalar_idx = fsr_id * _num_groups;
  long start_linear_idx = start_scalar_idx * 3;

  /* Set the FSR scalar flux buffer to zero */
  memset(fsr_flux, 0.0, _num_groups * 4 * sizeof(FP_PRECISION));

  if (_solve_3D) {

    /* Compute the segment midpoint */
    double center_x2[3];
    for (int i=0; i<3; i++)
      center_x2[i] = 2 * (position[i] + 0.5 * length * direction[i]);

    FP_PRECISION wgt = _quad->getWeightInline(azim_index, polar_index);
    FP_PRECISION length_2D = exp_evaluator->convertDistance3Dto2D(length);

    // Compute the exponential terms
    FP_PRECISION* exponentials = &scratch_pad[0];
    for (int e=0; e < _num_groups; e++) {
      int idx = 3*e;
      FP_PRECISION tau = sigma_t[e] * length_2D;
      exp_evaluator->retrieveExponentialComponents(tau, 0, &scratch_pad[idx],
                                                   &scratch_pad[idx+1],
                                                   &scratch_pad[idx+2]);
    }

    // Compute the sources
    FP_PRECISION* sources = &scratch_pad[3*_num_groups];
    for (int e=0; e < _num_groups; e++) {

      // Compute indexes upfront for performance
      long scalar_idx = start_scalar_idx + e;
      long first_linear_idx = start_linear_idx + 3*e;
      int flat_idx = 2*e;

      // Compute the flat and linear components of the source
      sources[flat_idx] = 0.0;
      sources[flat_idx+1] = 0.0;
      for (int i=0; i<3; i++) {
        int linear_idx = first_linear_idx + i;
        sources[flat_idx] += _reduced_sources_xyz[linear_idx] * center_x2[i];
        sources[flat_idx+1] += _reduced_sources_xyz[linear_idx] * direction[i];
      }
      sources[flat_idx] += _reduced_sources[scalar_idx];
    }

    for (int e=0; e < _num_groups; e++) {

      // Load sources
      int first_idx = 2*e;
      FP_PRECISION src_flat = sources[first_idx];
      FP_PRECISION src_linear = sources[first_idx+1];

      // Load exponential terms
      first_idx += e; // equivalent to 3*e
      FP_PRECISION exp_F1 = exponentials[first_idx];
      FP_PRECISION exp_F2 = exponentials[first_idx+1];
      FP_PRECISION exp_H = exponentials[first_idx+2];
      FP_PRECISION tau = sigma_t[e] * length_2D;

      // Compute the change in flux across the segment
      exp_H *= length * track_flux[e] * length_2D * wgt;
      FP_PRECISION delta_psi = (tau * track_flux[e] - length_2D * src_flat) *
          exp_F1 - src_linear * length_2D * length_2D * exp_F2;

      // Increment the fsr scalar flux and scalar flux moments
      first_idx += e; // equivalent to 4*e
      fsr_flux[first_idx] += wgt * delta_psi;
      first_idx++;
      FP_PRECISION reduced_delta = wgt * delta_psi / sigma_t[e];
      for (int i=0; i<3; i++)
        fsr_flux[first_idx + i] += exp_H * direction[i] + reduced_delta * 
                                    position[i];

      // Decrement the track flux
      track_flux[e] -= delta_psi;
    }
  }
  else {

    int pe = 0;

    /* Compute the segment midpoint */
    double center[2];
    for (int i=0; i<2; i++)
      center[i] = position[i] + 0.5 * length * direction[i];

    // Loop over energy groups
    int num_polar_2 = _num_polar / 2;
    for (int e=0; e < _num_groups; e++) {

      FP_PRECISION tau = sigma_t[e] * length;
      int exp_index = exp_evaluator->getExponentialIndex(tau);
      FP_PRECISION dt = exp_evaluator->getDifference(exp_index, tau);
      FP_PRECISION dt2 = dt * dt;

      FP_PRECISION polar_wgt_d_psi = 0.0;
      FP_PRECISION polar_wgt_exp_h = 0.0;

      // Compute the flat component of the reduced source
      FP_PRECISION src_flat = 0.0;
      for (int i=0; i<2; i++)
        src_flat += _reduced_sources_xyz(fsr_id, e, i) * center[i];

      src_flat *= 2;
      src_flat += _reduced_sources(fsr_id, e);

      // Loop over polar angles
      for (int p=0; p < num_polar_2; p++) {

        /* Get the sine of the polar angle */
        FP_PRECISION sin_theta = _quad->getSinTheta(azim_index, p);

        // Compute the exponential terms
        FP_PRECISION exp_F1 = exp_evaluator->computeExponentialF1(exp_index, p,
                                                                  dt, dt2);
        FP_PRECISION exp_F2 = exp_evaluator->computeExponentialF2(exp_index, p,
                                                                  dt, dt2);
        FP_PRECISION exp_H = exp_evaluator->computeExponentialH(exp_index, p,
                                                                 dt, dt2)
            * tau * length * track_flux[pe];

        // Compute the moment component of the source
        FP_PRECISION src_linear = 0.0;
        for (int i=0; i<2; i++)
          src_linear += direction[i] * sin_theta *
              _reduced_sources_xyz(fsr_id, e, i);

        // Compute the change in flux across the segment
        FP_PRECISION delta_psi = (tau * track_flux[pe] - length * src_flat)
            * exp_F1 - length * length * src_linear * exp_F2;

        FP_PRECISION track_weight = _quad->getWeightInline(azim_index, p);
        polar_wgt_d_psi += track_weight * delta_psi;
        polar_wgt_exp_h += track_weight * exp_H;

        track_flux[pe] -= delta_psi;
        pe++;
      }

      // Increment the fsr scalar flux and scalar flux moments
      fsr_flux[e*4] += polar_wgt_d_psi;
      for (int i=0; i<2; i++)
        fsr_flux[e*4 + i + 1] += polar_wgt_exp_h * direction[i]
              + polar_wgt_d_psi * position[i];
    }
  }

  // Atomically increment the FSR scalar flux from the temporary array
  omp_set_lock(&_FSR_locks[fsr_id]);

  for (int e=0; e < _num_groups; e++) {

    // Compute indexes upfront for performance
    long scalar_idx = start_scalar_idx + e;
    long first_linear_idx = start_linear_idx + 3*e;
    int first_idx = e*4;

    // Add to global scalar flux vector
    _scalar_flux[scalar_idx] += fsr_flux[first_idx];
    first_idx++;
    for (int i=0; i<3; i++)
      _scalar_flux_xyz[first_linear_idx+i] += fsr_flux[first_idx + i];
  }

  omp_unset_lock(&_FSR_locks[fsr_id]);

  for (int i=0; i < 3; i++)
    position[i] += direction[i] * length;
}


/**
 * @brief Add the source term contribution in the transport equation to
 *        the FSR scalar flux.
 */
void CPULSSolver::addSourceToScalarFlux() {

  int nc = 3;
  if (_solve_3D)
    nc = 6;

#pragma omp parallel
  {
    FP_PRECISION volume;
    FP_PRECISION flux_const;
    FP_PRECISION* sigma_t;

    /* Add in source term and normalize flux to volume for each FSR */
    /* Loop over FSRs, energy groups */
#pragma omp for
    for (long r=0; r < _num_FSRs; r++) {
      volume = _FSR_volumes[r];
      sigma_t = _FSR_materials[r]->getSigmaT();

      for (int e=0; e < _num_groups; e++) {

        flux_const = FOUR_PI * 2 / sigma_t[e];

        _scalar_flux(r,e) /= volume;
        _scalar_flux(r,e) += (FOUR_PI * _reduced_sources(r,e));
        _scalar_flux(r,e) /= sigma_t[e];

        _scalar_flux_xyz(r,e,0) /= volume;
        _scalar_flux_xyz(r,e,0) += flux_const * _reduced_sources_xyz(r,e,0)
            * _FSR_source_constants[r*_num_groups*nc + nc*e    ];
        _scalar_flux_xyz(r,e,0) += flux_const * _reduced_sources_xyz(r,e,1)
            * _FSR_source_constants[r*_num_groups*nc + nc*e + 2];

        _scalar_flux_xyz(r,e,1) /= volume;
        _scalar_flux_xyz(r,e,1) += flux_const * _reduced_sources_xyz(r,e,0)
            * _FSR_source_constants[r*_num_groups*nc + nc*e + 2];
        _scalar_flux_xyz(r,e,1) += flux_const * _reduced_sources_xyz(r,e,1)
            * _FSR_source_constants[r*_num_groups*nc + nc*e + 1];

        if (_solve_3D) {
          _scalar_flux_xyz(r,e,0) += flux_const * _reduced_sources_xyz(r,e,2)
              * _FSR_source_constants[r*_num_groups*nc + nc*e + 3];
          _scalar_flux_xyz(r,e,1) += flux_const * _reduced_sources_xyz(r,e,2)
              * _FSR_source_constants[r*_num_groups*nc + nc*e + 4];

          _scalar_flux_xyz(r,e,2) /= volume;
          _scalar_flux_xyz(r,e,2) += flux_const * _reduced_sources_xyz(r,e,0)
              * _FSR_source_constants[r*_num_groups*nc + nc*e + 3];
          _scalar_flux_xyz(r,e,2) += flux_const * _reduced_sources_xyz(r,e,1)
              * _FSR_source_constants[r*_num_groups*nc + nc*e + 4];
          _scalar_flux_xyz(r,e,2) += flux_const * _reduced_sources_xyz(r,e,2)
              * _FSR_source_constants[r*_num_groups*nc + nc*e + 5];
        }
      }
    }
  }
}


/**
 * @brief Computes the stabalizing flux for transport stabalization
 */
void CPULSSolver::computeStabalizingFlux() {

  /* Compute flat stabalizing flux */
  CPUSolver::computeStabalizingFlux();

  /* Check whether moment stabalization is requested */ 
  if (!_stabalize_moments)
    return;

  if (_stabalization_type == DIAGONAL) {
    /* Loop over all flat source regions, compute stabalizing flux moments */
#pragma omp parallel for
    for (long r=0; r < _num_FSRs; r++) {

      /* Extract the scattering matrix */
      FP_PRECISION* scattering_matrix = _FSR_materials[r]->getSigmaS();
    
      /* Extract total cross-sections */
      FP_PRECISION* sigma_t = _FSR_materials[r]->getSigmaT();

      for (int e=0; e < _num_groups; e++) {
      
        /* Extract the in-scattering (diagonal) element */
        FP_PRECISION sigma_s = scattering_matrix[e*_num_groups+e];
      
        /* For negative cross-sections, add the absolute value of the 
           in-scattering rate to the stabalizing flux */
        if (sigma_s < 0.0) {
          for (int i=0; i < 3; i++) {
            _stabalizing_flux_xyz(r, e, i) = -_scalar_flux_xyz(r,e,i) * 
                _stabalization_factor * sigma_s / sigma_t[e];
          }
        }
      }
    }
  }
  else if (_stabalization_type == YAMAMOTO) {

    /* Treat each group */
#pragma omp parallel for
    for (int e=0; e < _num_groups; e++) {

      /* Look for largest absolute scattering ratio */
      FP_PRECISION max_ratio = 0.0;
      for (long r=0; r < _num_FSRs; r++) {
        
        /* Extract the scattering value */
        FP_PRECISION scat = _FSR_materials[r]->getSigmaSByGroup(e+1, e+1);
    
        /* Extract total cross-sections */
        FP_PRECISION total = _FSR_materials[r]->getSigmaTByGroup(e+1);

        /* Determine scattering ratio */
        FP_PRECISION ratio = std::abs(scat / total);
        if (ratio > max_ratio)
          ratio = max_ratio;
      }
      max_ratio *= _stabalization_factor;
      for (long r=0; r < _num_FSRs; r++) {
        for (int i=0; i < 3; i++) {
          _stabalizing_flux_xyz(r, e, i) = _scalar_flux_xyz(r,e,i) * max_ratio;
        }
      }
    }
  }
  else if (_stabalization_type == GLOBAL) {
    
    /* Get the multiplicative factor */
    FP_PRECISION mult_factor = 1.0 / _stabalization_factor - 1.0;
   
    /* Apply the global muliplicative factor */ 
#pragma omp parallel for
    for (long r=0; r < _num_FSRs; r++)
      for (int e=0; e < _num_groups; e++)
        for (int i=0; i <3; i++)
          _stabalizing_flux_xyz(r, e, i) = _scalar_flux_xyz(r, e, i)
             * mult_factor;
  }
}
      

/**
 * @brief Adjusts the scalar flux for transport stabalization
 */
void CPULSSolver::stabalizeFlux() {

  /* Stabalize the flat scalar flux */
  CPUSolver::stabalizeFlux();

  /* Check whether moment stabalization is requested */ 
  if (!_stabalize_moments)
    return;

  if (_stabalization_type == DIAGONAL) {
    /* Loop over all flat source regions, apply stabalizing flux moments */
#pragma omp parallel for
    for (long r=0; r < _num_FSRs; r++) {

      /* Extract the scattering matrix */
      FP_PRECISION* scattering_matrix = _FSR_materials[r]->getSigmaS();
    
      /* Extract total cross-sections */
      FP_PRECISION* sigma_t = _FSR_materials[r]->getSigmaT();
    
      for (int e=0; e < _num_groups; e++) {
      
        /* Extract the in-scattering (diagonal) element */
        FP_PRECISION sigma_s = scattering_matrix[e*_num_groups+e];
      
        /* For negative cross-sections, add the stabalizing flux
           and divide by the diagonal matrix element used to form it so that
           no bias is introduced but the source iteration is stabalized */
        if (sigma_s < 0.0) {
          for (int i=0; i < 3; i++) {
            _scalar_flux_xyz(r, e, i) += _stabalizing_flux_xyz(r, e, i);
            _scalar_flux_xyz(r, e, i) /= (1.0 - _stabalization_factor * sigma_s /
                                         sigma_t[e]);
          }
        }
      }
    }
  }
  else if (_stabalization_type == YAMAMOTO) {

    /* Treat each group */
#pragma omp parallel for
    for (int e=0; e < _num_groups; e++) {

      /* Look for largest absolute scattering ratio */
      FP_PRECISION max_ratio = 0.0;
      for (long r=0; r < _num_FSRs; r++) {
        
        /* Extract the scattering value */
        FP_PRECISION scat = _FSR_materials[r]->getSigmaSByGroup(e+1, e+1);
    
        /* Extract total cross-sections */
        FP_PRECISION total = _FSR_materials[r]->getSigmaTByGroup(e+1);

        /* Determine scattering ratio */
        FP_PRECISION ratio = std::abs(scat / total);
        if (ratio > max_ratio)
          ratio = max_ratio;
      }
      max_ratio *= _stabalization_factor;
      for (long r=0; r < _num_FSRs; r++) {
        for (int i=0; i < 3; i++) {
          _scalar_flux_xyz(r, e, i) += _stabalizing_flux_xyz(r, e, i);
          _scalar_flux_xyz(r, e, i) /= (1 + max_ratio);
        }
      }
    }
  }
  else if (_stabalization_type == GLOBAL) {

    /* Apply the damping factor */    
#pragma omp parallel for
    for (long r=0; r < _num_FSRs; r++) {
      for (int e=0; e < _num_groups; e++) {
        for (int i=0; i <3; i++) {
          _scalar_flux_xyz(r, e, i) += _stabalizing_flux_xyz(r, e, i);
          _scalar_flux_xyz(r, e, i) *= _stabalization_factor;
        }
      }
    }
  }
}


/**
 * @brief Get the flux at a specific point in the geometry.
 * @param coords The coords of the point to get the flux at
 * @param group the energy group
 */
double CPULSSolver::getFluxByCoords(LocalCoords* coords, int group) {

  double x, y, z, xc, yc, zc;

  coords->setUniverse(_geometry->getRootUniverse());
  Cell* cell = _geometry->findCellContainingCoords(coords);
  long fsr = _geometry->getFSRId(coords);
  Point* centroid = _geometry->getFSRCentroid(fsr);
  x = coords->getX();
  y = coords->getY();
  z = coords->getZ();
  xc = centroid->getX();
  yc = centroid->getY();
  zc = centroid->getZ();

  double flux = _scalar_flux(fsr, group);
  double flux_x = 0.0;
  double flux_y = 0.0;
  double flux_z = 0.0;


  if (_solve_3D) {
    flux_x = (x - xc) *
        (_FSR_lin_exp_matrix[fsr*6  ] * _scalar_flux_xyz(fsr, group, 0) +
         _FSR_lin_exp_matrix[fsr*6+2] * _scalar_flux_xyz(fsr, group, 1) +
         _FSR_lin_exp_matrix[fsr*6+3] * _scalar_flux_xyz(fsr, group, 2));
    flux_y = (y - yc) *
        (_FSR_lin_exp_matrix[fsr*6+2] * _scalar_flux_xyz(fsr, group, 0) +
         _FSR_lin_exp_matrix[fsr*6+1] * _scalar_flux_xyz(fsr, group, 1) +
         _FSR_lin_exp_matrix[fsr*6+4] * _scalar_flux_xyz(fsr, group, 2));
    flux_z = (z - zc) *
        (_FSR_lin_exp_matrix[fsr*6+3] * _scalar_flux_xyz(fsr, group, 0) +
         _FSR_lin_exp_matrix[fsr*6+4] * _scalar_flux_xyz(fsr, group, 1) +
         _FSR_lin_exp_matrix[fsr*6+5] * _scalar_flux_xyz(fsr, group, 2));
  }
  else {
    flux_x = (x - xc) *
        (_FSR_lin_exp_matrix[fsr*3  ] * _scalar_flux_xyz(fsr, group, 0) +
         _FSR_lin_exp_matrix[fsr*3+2] * _scalar_flux_xyz(fsr, group, 1));
    flux_y = (y - yc) *
        (_FSR_lin_exp_matrix[fsr*3+2] * _scalar_flux_xyz(fsr, group, 0) +
         _FSR_lin_exp_matrix[fsr*3+1] * _scalar_flux_xyz(fsr, group, 1));
  }

  flux += flux_x + flux_y + flux_z;
  return flux;
}


//FIXME
void CPULSSolver::initializeCmfd() {
  Solver::initializeCmfd();
  if (_cmfd != NULL)
    _cmfd->setFluxMoments(_scalar_flux_xyz);
}


//FIXME
void CPULSSolver::initializeExpEvaluators() {
  for (int a=0; a < _num_exp_evaluators_azim; a++)
    for (int p=0; p < _num_exp_evaluators_polar; p++)
      _exp_evaluators[a][p]->useLinearSource();
  Solver::initializeExpEvaluators();
}


/**
 FIXME
 */
double* CPULSSolver::getLinearExpansionCoeffsBuffer() {
#pragma omp critical
  {
    if (_FSR_lin_exp_matrix == NULL) {
      long size = _geometry->getNumFSRs() * 3;
      if (_solve_3D)
        size *= 2;
      _FSR_lin_exp_matrix = new double[size];
      memset(_FSR_lin_exp_matrix, 0., size * sizeof(double));
    }
  }

  return _FSR_lin_exp_matrix;
}


/**
FIXME
 */
FP_PRECISION* CPULSSolver::getSourceConstantsBuffer() {
#pragma omp critical
  {
    if (_FSR_source_constants == NULL) {
      long size = 3 * _geometry->getNumEnergyGroups() * _geometry->getNumFSRs();
      if (_solve_3D)
        size *= 2;

      long max_size = size;
#ifdef MPIX
      if (_geometry->isDomainDecomposed())
        MPI_Allreduce(&size, &max_size, 1, MPI_LONG, MPI_MAX,
                      _geometry->getMPICart());
#endif
      double max_size_mb = (double) (max_size * sizeof(FP_PRECISION))
          / (double) (1e6);
      log_printf(NORMAL, "Max linear constant storage per domain = %6.2f MB",
                 max_size_mb);

      _FSR_source_constants = new FP_PRECISION[size];
      memset(_FSR_source_constants, 0., size * sizeof(FP_PRECISION));
    }
  }

  return _FSR_source_constants;
}
