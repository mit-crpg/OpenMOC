#include "TrackTraversingAlgorithms.h"
#include "CPUSolver.h"
#include "Quadrature.h"


/**
 * @brief Constructor for VolumeCalculator calls the TraverseTracks
 *        constructor
 * @param track_generator The TrackGenerator to pull tracking information
 */
VolumeCalculator::VolumeCalculator(TrackGenerator* track_generator)
                                  : TraverseTracks(track_generator) {
}


/**
 * @brief FSR volumes are calculated and saved in the TrackGenerator's FSR
 *        volumes buffer
 * @details VolumeKernels are created and used to loop over all segments and
 *          tally each segments contribution to FSR volumes.
 */
void VolumeCalculator::execute() {
#pragma omp parallel
  {
    VolumeKernel kernel(_track_generator);
    loopOverTracks(&kernel);
  }
}


/**
 * @brief Constructor for TransportSweep calls the TraverseTracks
 *        constructor and initializes the associated CPUSolver to NULL
 * @param track_generator The TrackGenerator to pull tracking information from
 */
TransportSweepFS::TransportSweepFS(TrackGenerator* track_generator)
    : TraverseTracks(track_generator) {
  _cpu_solver = NULL;
}


/**
 * @brief MOC equations are applied to every segment in the TrackGenerator
 * @details onTrack(...) applies the MOC equations to each segment and
 *          transfers boundary fluxes for the corresponding Track.
 */
void TransportSweepFS::execute() {
#pragma omp parallel
  {
    loopOverTracks(NULL);
  }
}


/**
 * @brief Sets the CPUSolver so that TransportSweep can apply MOC equations
 * @details This allows TransportSweep to transfer boundary fluxes from the
 *          CPUSolver and tally scalar fluxes
 * @param cpu_solver The CPUSolver which applies the MOC equations
 */
void TransportSweepFS::setCPUFSSolver(CPUFSSolver* cpu_solver) {
  _cpu_solver = cpu_solver;
}


/**
 * @brief Applies the MOC equations the Track and segments
 * @details The MOC equations are applied to each segment, attenuating the
 *          Track's angular flux and tallying FSR contributions. Finally,
 *          Track boundary fluxes are transferred.
 * @param track The Track for which the angular flux is attenuated and
 *        transferred
 * @param segments The segments over which the MOC equations are applied
 */
void TransportSweepFS::onTrack(Track* track, segment* segments) {

  /* Allocate temporary FSR flux locally */
  int num_groups = _track_generator->getGeometry()->getNumEnergyGroups();
  FP_PRECISION thread_fsr_flux[num_groups];

  /* Extract Track information */
  int track_id = track->getUid();
  int azim_index = track->getAzimAngleIndex();
  int num_segments = track->getNumSegments();
  FP_PRECISION* track_flux;

  /* Correct azimuthal index to first octant */
  Quadrature* quad = _track_generator->getQuadrature();

  /* Get the forward track flux */
  track_flux = _cpu_solver->getBoundaryFlux(track_id, true);

  /* Loop over each Track segment in forward direction */
  for (int s=0; s < num_segments; s++) {
    segment* curr_segment = &segments[s];
    _cpu_solver->tallyScalarFlux(curr_segment, azim_index, track_flux,
                                 thread_fsr_flux);
    _cpu_solver->tallyCurrent(curr_segment, azim_index, track_flux, true);
  }

  /* Transfer boundary angular flux to outgoing Track */
  _cpu_solver->transferBoundaryFlux(track_id, azim_index, true, track_flux);

  /* Get the backward track flux */
  track_flux = _cpu_solver->getBoundaryFlux(track_id, false);

  /* Loop over each Track segment in reverse direction */
  for (int s=num_segments-1; s >= 0; s--) {
    segment* curr_segment = &segments[s];
    _cpu_solver->tallyScalarFlux(curr_segment, azim_index, track_flux,
                                 thread_fsr_flux);
    _cpu_solver->tallyCurrent(curr_segment, azim_index, track_flux, false);
  }

  /* Transfer boundary angular flux to outgoing Track */
  _cpu_solver->transferBoundaryFlux(track_id, azim_index, false, track_flux);
}


/**
 * @brief Constructor for TransportSweep calls the TraverseTracks
 *        constructor and initializes the associated CPUSolver to NULL
 * @param track_generator The TrackGenerator to pull tracking information from
 */
TransportSweepLS::TransportSweepLS(TrackGenerator* track_generator)
    : TraverseTracks(track_generator) {
  _cpu_solver = NULL;
}


/**
 * @brief MOC equations are applied to every segment in the TrackGenerator
 * @details onTrack(...) applies the MOC equations to each segment and
 *          transfers boundary fluxes for the corresponding Track.
 */
void TransportSweepLS::execute() {
#pragma omp parallel
  {
    loopOverTracks(NULL);
  }
}


/**
 * @brief Sets the CPUSolver so that TransportSweep can apply MOC equations
 * @details This allows TransportSweep to transfer boundary fluxes from the
 *          CPUSolver and tally scalar fluxes
 * @param cpu_solver The CPUSolver which applies the MOC equations
 */
void TransportSweepLS::setCPULSSolver(CPULSSolver* cpu_solver) {
  _cpu_solver = cpu_solver;
}


/**
 * @brief Applies the MOC equations the Track and segments
 * @details The MOC equations are applied to each segment, attenuating the
 *          Track's angular flux and tallying FSR contributions. Finally,
 *          Track boundary fluxes are transferred.
 * @param track The Track for which the angular flux is attenuated and
 *        transferred
 * @param segments The segments over which the MOC equations are applied
 */
void TransportSweepLS::onTrack(Track* track, segment* segments) {

  segment* curr_segment;
  FP_PRECISION length;
  double x, y;
  Point* centroid;
  Point** FSR_centroids = _cpu_solver->getFSRCentroids();

  /* Allocate temporary FSR flux locally */
  int num_groups = _track_generator->getGeometry()->getNumEnergyGroups();
  FP_PRECISION thread_fsr_flux[num_groups*3];

  /* Extract Track information */
  int track_id = track->getUid();
  int azim_index = track->getAzimAngleIndex();
  int num_segments = track->getNumSegments();

  /* Correct azimuthal index to first octant */
  Quadrature* quad = _track_generator->getQuadrature();
  double cos_phi = cos(quad->getPhi(azim_index));
  double sin_phi = sin(quad->getPhi(azim_index));

  /* Get the forward track flux */
  FP_PRECISION* track_flux = _cpu_solver->getBoundaryFlux(track_id, true);

  /* Get the starting point for the first sesgment in the global
   * coordinate system */
  double X = track->getStart()->getX();
  double Y = track->getStart()->getY();

  /* Loop over each Track segment in forward direction */
  for (int s=0; s < num_segments; s++) {
    curr_segment = &segments[s];

    length = curr_segment->_length;
    centroid = FSR_centroids[curr_segment->_region_id];

    /* Get the starting point of the segment in local coordinates */
    x = X - centroid->getX();
    y = Y - centroid->getY();

    _cpu_solver->tallyLSScalarFlux(curr_segment, azim_index, track_flux,
                                     thread_fsr_flux, x, y, 1);
    _cpu_solver->tallyCurrent(curr_segment, azim_index, track_flux, true);

    /* Increment the segment starting point to the next segment */
    X += length * cos_phi;
    Y += length * sin_phi;
  }

  /* Transfer boundary angular flux to outgoing Track */
  _cpu_solver->transferBoundaryFlux(track_id, azim_index, true, track_flux);

  /* Get the backward track flux */
  track_flux = _cpu_solver->getBoundaryFlux(track_id, false);

  /* Loop over each Track segment in reverse direction */
  for (int s=num_segments-1; s >= 0; s--) {
    segment* curr_segment = &segments[s];

    length = curr_segment->_length;
    centroid = FSR_centroids[curr_segment->_region_id];

    /* Get the starting point of the segment in local coordinates */
    x = X - centroid->getX();
    y = Y - centroid->getY();

    _cpu_solver->tallyLSScalarFlux(curr_segment, azim_index, track_flux,
                                     thread_fsr_flux, x, y, -1);
    _cpu_solver->tallyCurrent(curr_segment, azim_index, track_flux, false);

    X -= length * cos_phi;
    Y -= length * sin_phi;
  }

  /* Transfer boundary angular flux to outgoing Track */
  _cpu_solver->transferBoundaryFlux(track_id, azim_index, false, track_flux);
}
