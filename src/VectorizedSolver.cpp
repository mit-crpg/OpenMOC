#include "VectorizedSolver.h"


/**
 * @brief Constructor initializes NULL arrays for source, flux, etc.
 * @param track_generator an optional pointer to a TrackGenerator object
 */
VectorizedSolver::VectorizedSolver(TrackGenerator* track_generator) :
  CPUSolver(track_generator) {

  if (_cmfd != NULL)
    log_printf(ERROR, "The VectorizedSolver is not yet configured for CMFD");

  _delta_psi = NULL;
  _thread_taus = NULL;
  _thread_exponentials = NULL;

  if (track_generator != NULL)
    setTrackGenerator(track_generator);

  vmlSetMode(VML_EP);
}


/**
 * @brief Destructor deletes Track boundary angular flux and
 *        and FSR scalar flux and source arrays.
 */
VectorizedSolver::~VectorizedSolver() {

  if (_boundary_flux != NULL) {
    MM_FREE(_boundary_flux);
    _boundary_flux = NULL;
  }

  if (_scalar_flux != NULL && !_user_fluxes) {
    MM_FREE(_scalar_flux);
    _scalar_flux = NULL;
  }

  if (_old_scalar_flux != NULL) {
    MM_FREE(_old_scalar_flux);
    _old_scalar_flux = NULL;
  }

  if (_reduced_sources != NULL) {
    MM_FREE(_reduced_sources);
    _reduced_sources = NULL;
  }

  if (_fixed_sources != NULL) {
    MM_FREE(_fixed_sources);
    _fixed_sources = NULL;
  }

  if (_delta_psi != NULL) {
    MM_FREE(_delta_psi);
    _delta_psi = NULL;
  }

  if (_thread_taus != NULL) {
    MM_FREE(_thread_taus);
    _thread_taus = NULL;
  }

  if (_thread_exponentials != NULL) {
    MM_FREE(_thread_exponentials);
    _thread_exponentials = NULL;
  }
}


/**
 * @brief Returns the number of vector lengths required to fit the number
 *        of energy groups.
 * @details If the number of energy groups is 35 and the vector width is
 *          4, this method will return 9 since 9*4 = 36 is the nearest
 *          integer greater than or equal to 35.
 * @return The number of vector widths
 */
int VectorizedSolver::getNumVectorWidths() {
  return _num_vector_lengths;
}


/**
 * @brief Sets the Geometry for the Solver.
 * @param geometry a pointer to the Geometry
 */
void VectorizedSolver::setGeometry(Geometry* geometry) {

  CPUSolver::setGeometry(geometry);

  /* Compute the number of SIMD vector widths needed to fit energy groups */
  _num_vector_lengths = (_num_groups / VEC_LENGTH) + 1;

  /* Reset the number of energy groups by rounding up for the number
   * of vector widths needed to accomodate the energy groups */
  _num_groups = _num_vector_lengths * VEC_LENGTH;
  _polar_times_groups = _num_groups * _num_polar;
}


/**
 * @brief Allocates memory for the exponential linear interpolation table.
 */
void VectorizedSolver::initializeExpEvaluator() {

  CPUSolver::initializeExpEvaluator();

  /* Deallocates memory for the exponentials if it was allocated for a
   * previous simulation */
  if (_thread_exponentials != NULL)
    MM_FREE(_thread_exponentials);

  /* Allocates memory for an array of exponential values for each thread
   * - this is not used by default, but can be to allow for vectorized
   * evaluation of the exponentials. Unfortunately this does not appear
   * to give any performance boost. */
  int size = _num_threads * _polar_times_groups * sizeof(FP_PRECISION);
  _thread_exponentials = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
}


/**
 * @brief Aligns all Material cross-section data for SIMD vector instructions.
 * @param mode the solution type (FORWARD or ADJOINT)
 */
void VectorizedSolver::initializeMaterials(solverMode mode) {
  /* Build fission matrices */
  Solver::initializeMaterials(mode);

  std::map<int, Material*> materials = _geometry->getAllMaterials();
  std::map<int, Material*>::iterator m_iter;

  /* Iterate over each Material and replace its cross-section with a new one
   * array that is a multiple of VEC_LENGTH long */
  for (m_iter = materials.begin(); m_iter != materials.end(); ++m_iter)
    m_iter->second->alignData();
}


/**
 * @brief Allocates memory for Track boundary angular and FSR scalar fluxes.
 * @details Deletes memory for old flux arrays if they were allocated for a
 *          previous simulation.
 */
void VectorizedSolver::initializeFluxArrays() {

  /* Delete old flux arrays if they exist */
  if (_boundary_flux != NULL)
    MM_FREE(_boundary_flux);

  if (_scalar_flux != NULL && !_user_fluxes)
    MM_FREE(_scalar_flux);

  if (_old_scalar_flux != NULL)
    MM_FREE(_old_scalar_flux);

  if (_delta_psi != NULL)
    MM_FREE(_delta_psi);

  if (_thread_taus != NULL)
    MM_FREE(_thread_taus);

  int size;

  /* Allocate aligned memory for all flux arrays */
  try{

    size = 2 * _tot_num_tracks * _num_groups * _num_polar;
    size *= sizeof(FP_PRECISION);
    _boundary_flux = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);

    size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);
    _scalar_flux = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
    _old_scalar_flux = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);

    size = _num_threads * _num_groups * sizeof(FP_PRECISION);
    _delta_psi = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);

    size = _num_threads * _polar_times_groups * sizeof(FP_PRECISION);
    _thread_taus = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
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
void VectorizedSolver::initializeSourceArrays() {

  /* Delete old sources arrays if they exist */
  if (_reduced_sources != NULL)
    MM_FREE(_reduced_sources);
  if (_fixed_sources != NULL)
    MM_FREE(_fixed_sources);

  int size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);

  /* Allocate aligned memory for all source arrays */
  try{
    _reduced_sources = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
    _fixed_sources = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for FSR sources");
  }

  /* Initialize fixed sources to zero */
  memset(_fixed_sources, 0.0, size);

  /* Populate fixed source array with any user-defined sources */
  initializeFixedSources();
}


/**
 * @brief Populates array of fixed sources assigned by FSR.
 */
void VectorizedSolver::initializeFixedSources() {

  Solver::initializeFixedSources();

  int fsr_id, group;
  std::pair<int, int> fsr_group_key;
  std::map< std::pair<int, int>, FP_PRECISION >::iterator fsr_iter;

  /* Populate fixed source array with any user-defined sources */
  for (fsr_iter = _fix_src_FSR_map.begin();
       fsr_iter != _fix_src_FSR_map.end(); ++fsr_iter) {

    /* Get the FSR with an assigned fixed source */
    fsr_group_key = fsr_iter->first;
    fsr_id = fsr_group_key.first;
    group = fsr_group_key.second;

    if (group <= 0 || group > _num_groups)
      log_printf(ERROR,"Unable to use fixed source for group %d in "
                 "a %d energy group problem", group, _num_groups);

    if (fsr_id < 0 || fsr_id >= _num_FSRs)
      log_printf(ERROR,"Unable to use fixed source for FSR %d with only "
                 "%d FSRs in the geometry", fsr_id, _num_FSRs);

    _fixed_sources(fsr_id, group-1) = _fix_src_FSR_map[fsr_group_key];
  }
}


/**
 * @brief Initializes the FSR volumes and Materials array.
 */
void VectorizedSolver::initializeFSRs() {

  CPUSolver::initializeFSRs();

  /* Compute the number of SIMD vector widths needed to fit energy groups */
  _num_vector_lengths = (_num_groups / VEC_LENGTH) + 1;

  /* Reset the number of energy groups by rounding up for the number
   * of vector widths needed to accomodate the energy groups */
  _num_groups = _num_vector_lengths * VEC_LENGTH;
  _polar_times_groups = _num_groups * _num_polar;
}



/**
 * @brief Normalizes all FSR scalar fluxes and Track boundary angular
 *        fluxes to the total fission source (times \f$ \nu \f$).
 */
void VectorizedSolver::normalizeFluxes() {

  FP_PRECISION* nu_sigma_f;
  FP_PRECISION volume;
  FP_PRECISION tot_fission_source;
  FP_PRECISION norm_factor;

  int size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);
  FP_PRECISION* fission_sources = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);

  /* Compute total fission source for each FSR, energy group */
#pragma omp parallel for private(volume, nu_sigma_f) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {

    /* Get pointers to important data structures */
    nu_sigma_f = _FSR_materials[r]->getNuSigmaF();
    volume = _FSR_volumes[r];

    /* Loop over energy group vector lengths */
    for (int v=0; v < _num_vector_lengths; v++) {

      /* Loop over each energy group within this vector */
#pragma simd vectorlength(VEC_LENGTH)
      for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
        fission_sources(r,e) = nu_sigma_f[e] * _scalar_flux(r,e);

      /* Loop over each energy group within this vector */
#pragma simd vectorlength(VEC_LENGTH)
      for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
        fission_sources(r,e) *= volume;
    }
  }

  /* Compute the total fission source */
  size = _num_FSRs * _num_groups;
#ifdef SINGLE
  tot_fission_source = cblas_sasum(size, fission_sources, 1);
#else
  tot_fission_source = cblas_dasum(size, fission_sources, 1);
#endif

  /* Deallocate memory for fission source array */
  MM_FREE(fission_sources);

  /* Compute the normalization factor */
  norm_factor = 1.0 / tot_fission_source;

  log_printf(DEBUG, "Tot. Fiss. Src. = %f, Normalization factor = %f",
             tot_fission_source, norm_factor);

  /* Normalize the FSR scalar fluxes */
#ifdef SINGLE
  cblas_sscal(size, norm_factor, _scalar_flux, 1);
  cblas_sscal(size, norm_factor, _old_scalar_flux, 1);
#else
  cblas_dscal(size, norm_factor, _scalar_flux, 1);
  cblas_dscal(size, norm_factor, _old_scalar_flux, 1);
#endif

  /* Normalize the Track angular boundary fluxes */
  size = 2 * _tot_num_tracks * _num_polar * _num_groups;

#ifdef SINGLE
  cblas_sscal(size, norm_factor, _boundary_flux, 1);
#else
  cblas_dscal(size, norm_factor, _boundary_flux, 1);
#endif

  return;
}


/**
 * @brief Computes the total source (fission, scattering, fixed) in each FSR.
 * @details This method computes the total source in each FSR based on
 *          this iteration's current approximation to the scalar flux.
 */
void VectorizedSolver::computeFSRSources() {

#pragma omp parallel default(none)
  {
    int tid;
    Material* material;
    FP_PRECISION* sigma_t;
    FP_PRECISION* sigma_s;
    FP_PRECISION* fiss_mat;
    FP_PRECISION scatter_source, fission_source;

    int size = _num_groups * sizeof(FP_PRECISION);
    FP_PRECISION* fission_sources =
      (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);
    FP_PRECISION* scatter_sources =
      (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);

    /* For all FSRs, find the source */
#pragma omp for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

      tid = omp_get_thread_num();
      material = _FSR_materials[r];
      sigma_t = material->getSigmaT();
      sigma_s = material->getSigmaS();
      fiss_mat = material->getFissionMatrix();

      /* Compute scatter + fission source for group G */
      for (int G=0; G < _num_groups; G++) {
        for (int v=0; v < _num_vector_lengths; v++) {

#pragma simd vectorlength(VEC_LENGTH)
          for (int g=v*VEC_LENGTH; g < (v+1)*VEC_LENGTH; g++) {
            scatter_sources[g] = sigma_s[G*_num_groups+g] * _scalar_flux(r,g);
            fission_sources[g] = fiss_mat[G*_num_groups+g] * _scalar_flux(r,g);
          }
        }

#ifdef SINGLE
        scatter_source=cblas_sasum(_num_groups, scatter_sources, 1);
        fission_source=cblas_sasum(_num_groups, fission_sources, 1);
#else
        scatter_source=cblas_dasum(_num_groups, scatter_sources, 1);
        fission_source=cblas_dasum(_num_groups, fission_sources, 1);
#endif

        fission_source /= _k_eff;

        /* Compute total (scatter+fission+fixed) reduced source */
        _reduced_sources(r,G) = _fixed_sources(r,G);
        _reduced_sources(r,G) += scatter_source + fission_source;
        _reduced_sources(r,G) *= ONE_OVER_FOUR_PI / sigma_t[G];
      }
    }

    MM_FREE(fission_sources);
    MM_FREE(scatter_sources);
  }
}


/**
 * @brief Add the source term contribution in the transport equation to
 *        the FSR scalar flux
 */
void VectorizedSolver::addSourceToScalarFlux() {

  FP_PRECISION volume;
  FP_PRECISION* sigma_t;

  /* Add in source term and normalize flux to volume for each FSR */
  /* Loop over FSRs, energy groups */
#pragma omp parallel for private(volume, sigma_t) schedule(guided)
  for (int r=0; r < _num_FSRs; r++) {

    volume = _FSR_volumes[r];
    sigma_t = _FSR_materials[r]->getSigmaT();

    /* Loop over each energy group vector length */
    for (int v=0; v < _num_vector_lengths; v++) {

      /* Loop over energy groups within this vector */
#pragma simd vectorlength(VEC_LENGTH)
      for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
        _scalar_flux(r,e) *= 0.5;

      /* Loop over energy groups within this vector */
#pragma simd vectorlength(VEC_LENGTH)
      for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
        _scalar_flux(r,e) = _scalar_flux(r,e) / (sigma_t[e] * volume);

      /* Loop over energy groups within this vector */
#pragma simd vectorlength(VEC_LENGTH)
      for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
        _scalar_flux(r,e) += FOUR_PI * _reduced_sources(r,e);
    }
  }

  return;
}


/**
 * @brief Compute \f$ k_{eff} \f$ from successive fission sources.
 */
void VectorizedSolver::computeKeff() {

  FP_PRECISION fission;

  int size = _num_FSRs * sizeof(FP_PRECISION);
  FP_PRECISION* FSR_rates = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);

  size = _num_threads * _num_groups * sizeof(FP_PRECISION);
  FP_PRECISION* group_rates = (FP_PRECISION*)MM_MALLOC(size, VEC_ALIGNMENT);

#pragma omp parallel
  {

    int tid = omp_get_thread_num() * _num_groups;
    Material* material;
    FP_PRECISION* sigma;
    FP_PRECISION volume;

    /* Compute the new nu-fission rates in each FSR */
#pragma omp for schedule(guided)
    for (int r=0; r < _num_FSRs; r++) {

      volume = _FSR_volumes[r];
      material = _FSR_materials[r];
      sigma = material->getNuSigmaF();

      /* Loop over each energy group vector length */
      for (int v=0; v < _num_vector_lengths; v++) {

        /* Loop over energy groups within this vector */
#pragma simd vectorlength(VEC_LENGTH)
        for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
          group_rates[tid+e] = sigma[e] * _scalar_flux(r,e);
      }

#ifdef SINGLE
      FSR_rates[r] = cblas_sasum(_num_groups, &group_rates[tid], 1) * volume;
#else
      FSR_rates[r] = cblas_dasum(_num_groups, &group_rates[tid], 1) * volume;
#endif
    }
  }

  /* Reduce new fission rates across FSRs */
#ifdef SINGLE
  fission = cblas_sasum(_num_FSRs, FSR_rates, 1);
#else
  fission = cblas_dasum(_num_FSRs, FSR_rates, 1);
#endif

  _k_eff *= fission;

  MM_FREE(FSR_rates);
  MM_FREE(group_rates);
}



/**
 * @brief Computes the contribution to the FSR scalar flux from a segment.
 * @details This method integrates the angular flux for a Track segment across
 *        energy groups and polar angles, and tallies it into the FSR scalar
 *        flux, and updates the Track's angular flux.
 * @param curr_segment a pointer to the Track segment of interest
 * @param azim_index a pointer to the azimuthal angle index for this segment
 * @param track_flux a pointer to the Track's angular flux
 * @param fsr_flux a pointer to the temporary FSR flux buffer
 */
void VectorizedSolver::tallyScalarFlux(segment* curr_segment,
                                       int azim_index,
                                       FP_PRECISION* track_flux,
                                       FP_PRECISION* fsr_flux) {

  int tid = omp_get_thread_num();
  int fsr_id = curr_segment->_region_id;
  FP_PRECISION* delta_psi = &_delta_psi[tid*_num_groups];
  FP_PRECISION* exponentials = &_thread_exponentials[tid*_polar_times_groups];

  computeExponentials(curr_segment, exponentials);

  /* Set the FSR scalar flux buffer to zero */
  memset(fsr_flux, 0.0, _num_groups * sizeof(FP_PRECISION));

  /* Tally the flux contribution from segment to FSR's scalar flux */
  /* Loop over polar angles */
  for (int p=0; p < _num_polar; p++) {

    /* Loop over each energy group vector length */
    for (int v=0; v < _num_vector_lengths; v++) {

      /* Loop over energy groups within this vector */
#pragma simd vectorlength(VEC_LENGTH)
      for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
        delta_psi[e] = track_flux(p,e) - _reduced_sources(fsr_id,e);

      /* Loop over energy groups within this vector */
#pragma simd vectorlength(VEC_LENGTH)
      for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
        delta_psi[e] *= exponentials(p,e);

      /* Loop over energy groups within this vector */
#pragma simd vectorlength(VEC_LENGTH)
      for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
        fsr_flux[e] += delta_psi[e] * _polar_weights(azim_index,p);

      /* Loop over energy groups within this vector */
#pragma simd vectorlength(VEC_LENGTH)
      for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
        track_flux(p,e) -= delta_psi[e];
    }
  }

  /* Atomically increment the FSR scalar flux from the temporary array */
  omp_set_lock(&_FSR_locks[fsr_id]);
  {
#ifdef SINGLE
    vsAdd(_num_groups, &_scalar_flux(fsr_id,0), fsr_flux,
          &_scalar_flux(fsr_id,0));
#else
    vdAdd(_num_groups, &_scalar_flux(fsr_id,0), fsr_flux,
          &_scalar_flux(fsr_id,0));
#endif
  }
  omp_unset_lock(&_FSR_locks[fsr_id]);
}


/**
 * @brief Computes an array of the exponentials in the transport equation,
 *        \f$ exp(-\frac{\Sigma_t * l}{sin(\theta)}) \f$, for each energy group
 *        and polar angle for a given Track segment.
 * @param curr_segment pointer to the Track segment of interest
 * @param exponentials the array to store the exponential values
 */
void VectorizedSolver::computeExponentials(segment* curr_segment,
                                           FP_PRECISION* exponentials) {

  FP_PRECISION length = curr_segment->_length;
  FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();

  /* Evaluate the exponentials using the linear interpolation table */
  if (_exp_evaluator->isUsingInterpolation()) {
    FP_PRECISION tau;

    for (int e=0; e < _num_groups; e++) {
      tau = length * sigma_t[e];
      for (int p=0; p < _num_polar; p++)
        exponentials(p,e) = _exp_evaluator->computeExponential(tau, p);
    }
  }

  /* Evalute the exponentials using the intrinsic exp(...) function */
  else {

    int tid = omp_get_thread_num();
    FP_PRECISION* sin_thetas = _polar_quad->getSinThetas();
    FP_PRECISION* taus = &_thread_taus[tid*_polar_times_groups];

    /* Initialize the tau argument for the exponentials */
    for (int p=0; p < _num_polar; p++) {

      for (int v=0; v < _num_vector_lengths; v++) {

#pragma simd vectorlength(VEC_LENGTH)
        for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
          taus(p,e) = -sigma_t[e] * length;

#pragma simd vectorlength(VEC_LENGTH)
        for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
          taus(p,e) /= sin_thetas[p];
      }
    }

    /* Evaluate the negative of the exponentials using Intel's MKL */
#ifdef SINGLE
    vsExp(_polar_times_groups, taus, exponentials);
#else
    vdExp(_polar_times_groups, taus, exponentials);
#endif

    /* Compute one minus the exponentials */
    for (int p=0; p < _num_polar; p++) {

      for (int v=0; v < _num_vector_lengths; v++) {

#pragma simd vectorlength(VEC_LENGTH)
        for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
          exponentials(p,e) = 1.0 - exponentials(p,e);
      }
    }
  }
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
void VectorizedSolver::transferBoundaryFlux(int track_id, int azim_index,
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
  for (int p=0; p < _num_polar; p++) {

    /* Loop over each energy group vector length */
    for (int v=0; v < _num_vector_lengths; v++) {

      /* Loop over energy groups within this vector */
#pragma simd vectorlength(VEC_LENGTH)
      for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++)
        track_out_flux(p,e) = track_flux(p,e) * transfer_flux;
    }
  }
}
