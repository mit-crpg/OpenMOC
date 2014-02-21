#include "ThreadPrivateSolver.h"


/**
 * @brief Constructor initializes array pointers for Tracks and Materials.
 * @details The constructor retrieves the number of energy groups and FSRs
 *          and azimuthal angles from the Geometry and TrackGenerator if
 *          passed in as parameters by the user. The constructor initalizes
 *          the number of OpenMP threads to a default of 1.
 * @param geometry an optional pointer to the Geometry
 * @param track_generator an optional pointer to the TrackgGenerator
 * @param cmfd an optional pointer to a Cmfd object object
 */
ThreadPrivateSolver::ThreadPrivateSolver(Geometry* geometry,
                                         TrackGenerator* track_generator,
                                         Cmfd* cmfd) :
  CPUSolver(geometry, track_generator, cmfd) {

  _thread_flux = NULL;
  _thread_currents = NULL;
}


/**
 * @brief Destructor calls Solver subclass destructor to deletes arrays
 *        for fluxes and sources.
 */
ThreadPrivateSolver::~ThreadPrivateSolver() {

  if (_thread_flux != NULL) {

    for (int t=0; t < _num_threads; t++)
      delete [] _thread_flux[t];

    delete [] _thread_flux;
    _thread_flux = NULL;
  }

  if (_thread_currents != NULL) {
    delete [] _thread_currents;
    _thread_currents = NULL;
  }

}


/**
 * @brief Allocates memory for Track boundary angular flux and leakage and
 *        FSR scalar flux arrays.
 * @details Deletes memory for old flux arrays if they were allocated for a
 *          previous simulation.
 */
void ThreadPrivateSolver::initializeFluxArrays() {

  CPUSolver::initializeFluxArrays();

  /* Delete old flux arrays if they exist */
  if (_thread_flux != NULL) {
    for (int t=0; t < _num_threads; t++)
      delete [] _thread_flux[t];

    delete [] _thread_flux;
  }

  int size;

  /* Allocate memory for the flux and leakage arrays */
  try{

    /* Allocate a thread local array of FSR scalar fluxes */
    _thread_flux = new FP_PRECISION*[_num_threads];
    for (int t=0; t < _num_threads; t++)
      _thread_flux[t] = new FP_PRECISION[_num_FSRs * _num_groups];
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the Solver's fluxes. "
               "Backtrace:%s", e.what());
  }
}


/**
 * @brief Initializes Cmfd object for acceleration prior to source iteration.
 * @details Instantiates a dummy Cmfd object if one was not assigned to
 *          the Solver by the user and initializes FSRs, Materials, fluxes
 *          and the Mesh. This method intializes thread private arrays
 *          for the Cmfd Mesh surface currents.
 */
void ThreadPrivateSolver::initializeCmfd() {

  /* Call parent class method */
  CPUSolver::initializeCmfd();

  /* Delete old thread private Cmfd Mesh surface currents array it it exists */
  if (_thread_currents != NULL)
    delete [] _thread_currents;

  int size;

  /* Allocate memory for the thread private Cmfd Mesh surface currents array */
  try{
    /* Allocate a thread local array of Cmfd Mesh cell surface currents */
    if (_cmfd->getMesh()->getCmfdOn()){
      size = _num_threads * _num_mesh_cells * 8 * _cmfd->getNumCmfdGroups();
      _thread_currents = new FP_PRECISION[size];
    }
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not allocate memory for the Solver's Cmfd Mesh"
               " surface currents. Backtrace:%s", e.what());
  }

  return;
}



/**
 * @brief Set the FSR scalar flux for each energy group to some value.
 * @details This method also flattens the thread private FSR scalar flux array.
 * @param value the value to assign to each FSR scalar flux
 */
void ThreadPrivateSolver::flattenFSRFluxes(FP_PRECISION value) {

  CPUSolver::flattenFSRFluxes(value);

  /* Flatten the thread private FSR scalar flux array */
  #pragma omp parallel for schedule(guided)
  for (int tid=0; tid < _num_threads; tid++) {
    for (int r=0; r < _num_FSRs; r++) {
      for (int e=0; e < _num_groups; e++)
        _thread_flux(tid,r,e) = 0.0;
    }
  }

  return;
}


/**
 * @brief Set the surface currents for each energy group inside each Cmfd
 *        Mesh cell to zero.
 */
void ThreadPrivateSolver::zeroSurfaceCurrents() {

  CPUSolver::zeroSurfaceCurrents();

  #pragma omp parallel for schedule(guided)
  for (int tid=0; tid < _num_threads; tid++){
    for (int r=0; r < _num_mesh_cells; r++) {
      for (int s=0; s < 8; s++) {
        for (int e=0; e < _num_groups; e++)
          _thread_currents(tid,r*8+s,e) = 0.0;
      }
    }
  }

  return;
}


/**
 * @brief This method performs one transport sweep of all azimuthal angles,
 *        Tracks, Track segments, polar angles and energy groups.
 * @details The method integrates the flux along each track and updates the
 *          boundary fluxes for the corresponding output Track, while updating
 *          the scalar flux in each flat source region.
 */
void ThreadPrivateSolver::transportSweep() {

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

  if (_cmfd->getMesh()->getCmfdOn())
    zeroSurfaceCurrents();

  /* Loop over azimuthal angle halfspaces */
  for (int i=0; i < 2; i++) {

    /* Compute the minimum and maximum Track IDs corresponding to this
     * this azimuthal angular halfspace */
    int min = i * (_tot_num_tracks / 2);
    int max = (i + 1) * (_tot_num_tracks / 2);

    /* Loop over each thread within this azimuthal angle halfspace */
    #pragma omp parallel for private(tid, fsr_id, curr_track, azim_index, \
      num_segments, segments, curr_segment,  track_flux) schedule(guided)
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
                        &_thread_flux(tid,fsr_id,0),true);
      }

      /* Transfer boundary angular flux to outgoing track */
      transferBoundaryFlux(track_id, azim_index, true, track_flux);

     /* Loop over each Track segment in reverse direction */
      track_flux += _polar_times_groups;

      for (int s=num_segments-1; s > -1; s--) {
        curr_segment = &segments[s];
        fsr_id = curr_segment->_region_id;
        scalarFluxTally(curr_segment, azim_index, track_flux,
                        &_thread_flux(tid,fsr_id,0),false);
      }

      /* Transfer boundary angular flux to outgoing Track */
      transferBoundaryFlux(track_id, azim_index, false, track_flux);
    }
  }

  reduceThreadScalarFluxes();

  if (_cmfd->getMesh()->getCmfdOn())
    reduceThreadSurfaceCurrents();

  return;
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
 * @param fwd
 */
void ThreadPrivateSolver::scalarFluxTally(segment* curr_segment,
                                          int azim_index,
                                          FP_PRECISION* track_flux,
                                          FP_PRECISION* fsr_flux,
                                          bool fwd){

  int tid = omp_get_thread_num();
  int fsr_id = curr_segment->_region_id;
  FP_PRECISION length = curr_segment->_length;
  FP_PRECISION* sigma_t = curr_segment->_material->getSigmaT();

  /* The change in angular flux along this Track segment in the FSR */
  FP_PRECISION delta_psi;
  FP_PRECISION exponential;

  /* Loop over energy groups */
  for (int e=0; e < _num_groups; e++) {

    /* Loop over polar angles */
    for (int p=0; p < _num_polar; p++){
      exponential = computeExponential(sigma_t[e], length, p);
      delta_psi = (track_flux(p,e)-_reduced_source(fsr_id,e))*exponential;
      fsr_flux[e] += delta_psi * _polar_weights(azim_index,p);
      track_flux(p,e) -= delta_psi;
    }
  }

  if (_cmfd->getMesh()->getCmfdOn()){
    if (curr_segment->_mesh_surface_fwd != -1 && fwd){

      int pe = 0;

      /* Loop over energy groups */
      for (int e = 0; e < _num_groups; e++) {

        /* Loop over polar angles */
        for (int p = 0; p < _num_polar; p++){

          /* Increment current (polar and azimuthal weighted flux, group)*/
          _thread_currents(tid,curr_segment->_mesh_surface_fwd,e) +=
                             track_flux(p,e)*_polar_weights(azim_index, p)/2.0;
          pe++;
        }
      }
    }
    else if (curr_segment->_mesh_surface_bwd != -1 && !fwd){

      /* Set polar angle * energy group to 0 */
      int pe = 0;

      /* Loop over energy groups */
      for (int e = 0; e < _num_groups; e++) {

        /* Loop over polar angles */
        for (int p = 0; p < _num_polar; p++){

          /* increment current (polar and azimuthal weighted flux, group)*/
          _thread_currents(tid,curr_segment->_mesh_surface_bwd,e) +=
                             track_flux(p,e)*_polar_weights(azim_index, p)/2.0;

          pe++;
        }
      }
    }
  }

  return;
}


/**
 * @brief Reduces the FSR scalar fluxes from private thread private arrays to a
 *        global array FSR scalar flux array.
 */
void ThreadPrivateSolver::reduceThreadScalarFluxes() {

  for (int tid=0; tid < _num_threads; tid++) {
    for (int r=0; r < _num_FSRs; r++) {
      for (int e=0; e < _num_groups; e++)
        _scalar_flux(r,e) += _thread_flux(tid,r,e);
    }
  }

  return;
}


/**
 * @brief Reduces the Cmfd Mesh surface currents from private thread arrays to
 *        a global Mesh surface current array.
 */
void ThreadPrivateSolver::reduceThreadSurfaceCurrents() {

  for (int tid=0; tid < _num_threads; tid++){
    for (int r=0; r < _num_mesh_cells; r++) {
      for (int s=0; s < 8; s++) {
        for (int e=0; e < _cmfd->getNumCmfdGroups(); e++){
          _surface_currents[(r*8+s)*_cmfd->getNumCmfdGroups() + e] +=
                             _thread_currents[(tid)*_num_mesh_cells*8*
                             _cmfd->getNumCmfdGroups() + (r*8+s)*
                             _cmfd->getNumCmfdGroups() + e];
        }
      }
    }
  }

  return;
}
