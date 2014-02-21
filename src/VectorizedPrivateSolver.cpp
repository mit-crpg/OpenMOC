#include "VectorizedPrivateSolver.h"


/**
 * @brief Constructor initializes empty arrays for source, flux, etc.
 * @details The construcor retrieves the number of energy groups and FSRs
 *          and azimuthal angles from the Geometry and TrackGenerator if
 *          they were provided by the user, and uses this to initialize
 *          empty arrays for the FSRs, boundary angular fluxes, FSR scalar
 *          fluxes, FSR sources and FSR fission rates. The constructor
 *          initalizes the number of threads to a default of 1.
 * @param geometry an optional pointer to the Geometry
 * @param track_generator an optional pointer to the TrackGenerator
 * @param cmfd an optional pointer to a Cmfd object object
 */
VectorizedPrivateSolver::VectorizedPrivateSolver(Geometry* geometry,
                                                TrackGenerator* track_generator,
                                                Cmfd* cmfd) :

  VectorizedSolver(geometry, track_generator, cmfd) {

  _thread_flux = NULL;
}


/**
 * @brief Destructor deletes arrays of Track boundary angular flux and
 *        FSR scalar flux and source arrays.
 */
VectorizedPrivateSolver::~VectorizedPrivateSolver() {

  if (_thread_flux != NULL) {

    for (int t=0; t < _num_threads; t++)
      _mm_free(_thread_flux[t]);

    delete [] _thread_flux;
    _thread_flux = NULL;
  }

}


/**
 * @brief Set the scalar flux for each energy group inside each FSR
 *        to some value.
 * @details This method also flattens the thread private FSR scalar flux array.
 * @param value the value to assign to each FSR flux
 */
void VectorizedPrivateSolver::flattenFSRFluxes(FP_PRECISION value) {

  CPUSolver::flattenFSRFluxes(value);

  /* Flatten the thread private FSR scalar flux array */
  #pragma omp parallel for schedule(guided)
  for (int tid=0; tid < _num_threads; tid++) {
    for (int r=0; r < _num_FSRs; r++) {
      for (int e=0; e < _num_groups; e++) {
        _thread_flux(tid,r,e) = 0.0;
      }
    }
  }

  return;
}



/**
 * @brief Allocates memory for Track boundary angular fluxes and
 *        FSR scalar fluxes and leakages.
 * @details Deletes memory for old flux arrays if they were allocated for a
 *          previous simulation.
 */
void VectorizedPrivateSolver::initializeFluxArrays() {

  VectorizedSolver::initializeFluxArrays();

  /* Delete old flux arrays if they exist */
  if (_thread_flux != NULL) {

    for (int t=0; t < _num_threads; t++)
      _mm_free(_thread_flux[t]);

    delete [] _thread_flux;
  }

  int size;

  /* Allocate aligned memory for all flux arrays */
  try{
    _thread_flux = new FP_PRECISION*[_num_threads];
    size = _num_FSRs * _num_groups * sizeof(FP_PRECISION);

    for (int t=0; t < _num_threads; t++)
      _thread_flux[t] = (FP_PRECISION*)_mm_malloc(size, VEC_ALIGNMENT);
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the "
               "VectorizedPrivateSolver's fluxes. Backtrace:%s", e.what());
    }
}


/**
 * @brief Computes the contribution to the FSR scalar flux from a Track segment.
 * @details This method integrates the angular flux for a Track segment across
 *          energy groups and polar angles, and tallies it into the FSR scalar
 *          flux, and updates the Track's angular flux.
 * @param curr_segment a pointer to the Track segment of interest
 * @param azim_index a pointer to the azimuthal angle index for this segment
 * @param track_flux a pointer to the Track's angular flux
 * @param fsr_flux a pointer to the temporary FSR scalar flux buffer
 */
void VectorizedPrivateSolver::scalarFluxTally(segment* curr_segment,
                                              int azim_index,
                                              FP_PRECISION* track_flux,
                                              FP_PRECISION* fsr_flux){

  int tid = omp_get_thread_num();
  int fsr_id = curr_segment->_region_id;
  FP_PRECISION length = curr_segment->_length;
  FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();

  /* The change in angular flux along this Track segment in the FSR */
  FP_PRECISION delta_psi;
  FP_PRECISION* exponentials = &_thread_exponentials[tid*_polar_times_groups];

  computeExponentials(curr_segment, exponentials);

  /* Tally the flux contribution from segment to FSR's scalar flux */
  /* Loop over polar angles */
  for (int p=0; p < _num_polar; p++){

    /* Loop over each energy group vector length */
    for (int v=0; v < _num_vector_lengths; v++) {

      /* Loop over energy groups within this vector */
      #pragma simd vectorlength(VEC_LENGTH) private(delta_psi)
      for (int e=v*VEC_LENGTH; e < (v+1)*VEC_LENGTH; e++) {
        delta_psi = (track_flux(p,e) - _reduced_source(fsr_id,e)) *
                   exponentials(p,e);
        fsr_flux[e] += delta_psi * _polar_weights(azim_index,p);
        track_flux(p,e) -= delta_psi;
      }
    }
  }

  return;
}


/**
 * @brief This method performs one transport sweep of all azimuthal angles,
 *        Tracks, Track segments, polar angles and energy groups.
 * @details The method integrates the flux along each track and updates the
 *          Track boundary fluxes for the corresponding output track, while
 *          updating the scalar flux in each FSR.
 */
void VectorizedPrivateSolver::transportSweep() {

  int tid;
  int fsr_id;
  Track* curr_track;
  int azim_index;
  int num_segments;
  segment* curr_segment;
  segment* segments;
  FP_PRECISION* track_flux;

  log_printf(DEBUG, "Transport sweep with %d OpenMP threads", _num_threads);

  /* Initialize flux in each FSR to zero */
  flattenFSRFluxes(0.0);

  /* Loop over azimuthal angle halfspaces */
  for (int i=0; i < 2; i++) {

    /* Compute the minimum and maximum Track IDs corresponding to
     * this azimuthal angular halfspace */
    int min = i * (_tot_num_tracks / 2);
    int max = (i + 1) * (_tot_num_tracks / 2);

    /* Loop over each thread within this azimuthal angle halfspace */
    #pragma omp parallel for private(tid, fsr_id, curr_track, azim_index, \
      num_segments, segments, curr_segment, track_flux) schedule(guided)
    for (int track_id=min; track_id < max; track_id++) {

      tid = omp_get_thread_num();

      /* Initialize local pointers to important data structures */
      curr_track = _tracks[track_id];
      azim_index = curr_track->getAzimAngleIndex();
      num_segments = curr_track->getNumSegments();
      segments = curr_track->getSegments();
      track_flux = &_boundary_flux(track_id,0,0,0);

      /* Loop over each Track segment in forward direction */
      for (int s=0; s < num_segments; s++) {
        curr_segment = &segments[s];
        fsr_id = curr_segment->_region_id;
        scalarFluxTally(curr_segment, azim_index, track_flux,
                        &_thread_flux(tid,fsr_id,0));
      }

      /* Transfer flux to outgoing Track */
      transferBoundaryFlux(track_id, azim_index, true, track_flux);

      /* Loop over each Track segment in reverse direction */
      track_flux += _polar_times_groups;

      for (int s=num_segments-1; s > -1; s--) {
        curr_segment = &segments[s];
        fsr_id = curr_segment->_region_id;
        scalarFluxTally(curr_segment, azim_index, track_flux,
                        &_thread_flux(tid,fsr_id,0));
      }

      /* Transfer flux to outgoing Track */
      transferBoundaryFlux(track_id, azim_index, false, track_flux);
    }
  }

  reduceThreadScalarFluxes();

  return;
}


/**
 * @brief Reduces the FSR scalar fluxes from private thread
 *        array to a global array.
 */
void VectorizedPrivateSolver::reduceThreadScalarFluxes() {

  for (int tid=0; tid < _num_threads; tid++) {
    for (int r=0; r < _num_FSRs; r++) {
      for (int e=0; e < _num_groups; e++) {
        _scalar_flux(r,e) += _thread_flux(tid,r,e);
      }
    }
  }

  return;
}
