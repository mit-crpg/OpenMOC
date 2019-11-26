#include "MOCKernel.h"
#include "TrackGenerator3D.h"
#include "CPUSolver.h"

/**
 * @brief Constructor for the MOCKernel assigns default values
 * @param track_generator the TrackGenerator used to pull relevant tracking
 *        data from
 */
MOCKernel::MOCKernel(TrackGenerator* track_generator) {
  _count = 0;
  _max_tau = track_generator->retrieveMaxOpticalLength();
#ifndef NGROUPS
  _num_groups = track_generator->getGeometry()->getNumEnergyGroups();
#endif
}


/**
 * @brief Constructor for the VolumeKernel assigns default values, calls
 *        the MOCKernel constructor, and pulls references to FSR locks and FSR
 *        volumes from the provided TrackGenerator.
 * @param track_generator the TrackGenerator used to pull relevant tracking
 *        data from
 */
VolumeKernel::VolumeKernel(TrackGenerator* track_generator) :
                           MOCKernel(track_generator) {

  _FSR_locks = track_generator->getFSRLocks();

  if (_FSR_locks == NULL)
    log_printf(ERROR, "Unable to create a VolumeKernel without first creating "
               "FSR locks");

  _FSR_volumes = track_generator->getFSRVolumesBuffer();
  _quadrature = track_generator->getQuadrature();
  _weight = 0;
}


/**
 * @brief Constructor for the SegmentationKernel assigns default values, calls
 *        the MOCKernel constructor, and pulls a reference to temporary segment
 *        data from the provided TrackGenerator.
 * @param track_generator the TrackGenerator used to pull relevant tracking
 *        data from
 */
SegmentationKernel::SegmentationKernel(TrackGenerator* track_generator)
                                       : MOCKernel(track_generator) {

  int thread_id = omp_get_thread_num();
  TrackGenerator3D* track_generator_3D =
    dynamic_cast<TrackGenerator3D*>(track_generator);
  if (track_generator_3D != NULL)
    _segments = track_generator_3D->getTemporarySegments(thread_id);
  else
    _segments = NULL;
}


/**
 * @brief Constructor for the CounterKernel assigns default values and calls
 *        the MOCKernel constructor.
 * @param track_generator the TrackGenerator used to pull relevant tracking
 *        data from
 */
CounterKernel::CounterKernel(TrackGenerator* track_generator) :
                           MOCKernel(track_generator) {}


/**
 * @brief Prepares an MOCKernel for a new Track.
 * @details Resets the segment count
 * @param track The new Track the MOCKernel prepares to handle
 */
void MOCKernel::newTrack(Track* track) {
  _count = 0;
}


/**
 * @brief Prepares a VolumeKernel for a new Track.
 * @details Resets the segment count and updates the weight for the new Track
 * @param track The new Track the MOCKernel prepares to handle
 */
void VolumeKernel::newTrack(Track* track) {

  /* Compute the Track azimuthal weight */
  int azim_index = track->getAzimIndex();
  _weight = _quadrature->getAzimSpacing(azim_index)
      * _quadrature->getAzimWeight(azim_index);

  /* Multiply by polar weight if 3D */
  Track3D* track_3D = dynamic_cast<Track3D*>(track);
  if (track_3D != NULL) {
    int polar_index = track_3D->getPolarIndex();
    _weight *= _quadrature->getPolarSpacing(azim_index, polar_index)
        * _quadrature->getPolarWeight(azim_index, polar_index);
  }

  /* Reset the count //NOTE Segment cuts arent counted. */
  _count = 0;
}


/**
 * @brief Destructor for MOCKernel.
 */
MOCKernel::~MOCKernel() {};


/**
 * @brief Reads and returns the current count.
 * @details MOC kernels count how many times they are accessed. This value
 *          returns the value of the counter (number of execute accesses)
 *          since kernel creation or last call to newTrack.
 * @return _count the counter value
 */
int MOCKernel::getCount() {
  return _count;
}


/**
 * @brief Resets the maximum optcal path length for a segment.
 * @details MOC kernels ensure that there are no segments with an optical path
 *          length greater than the maximum optical path length by splitting
 *          them when they get too large.
 * @param max_tau the maximum optical path length for a segment
 */
void MOCKernel::setMaxOpticalLength(FP_PRECISION max_tau) {
  _max_tau = max_tau;
}


/**
 * @brief Adds segment contribution to the FSR volume.
 * @details The VolumeKernel execute function adds the product of the
 *          track length and track weight to the buffer array at index
 *          id, referring to the array of FSR volumes.
 * @param length segment length
 * @param mat Material associated with the segment
 * @param fsr_id the FSR ID of the FSR associated with the segment
 */
void VolumeKernel::execute(FP_PRECISION length, Material* mat, long fsr_id,
                           int track_idx, int cmfd_surface_fwd,
                           int cmfd_surface_bwd, FP_PRECISION x_start,
                           FP_PRECISION y_start, FP_PRECISION z_start,
                           FP_PRECISION phi, FP_PRECISION theta) {

  /* Set omp lock for FSRs */
  omp_set_lock(&_FSR_locks[fsr_id]);

  /* Add value to buffer */
  _FSR_volumes[fsr_id] += _weight * length;

  /* Unset lock */
  omp_unset_lock(&_FSR_locks[fsr_id]);

  /* Increment the count of the number of times the kernel is executed */
  _count++;
}


/**
 * @brief Increments the counter for the number of segments on the track.
 * @details The CounterKernel execute function counts the number of segments
 *          in a track by incrementing the counter variable upon execution. Due
 *          to restrictions on maximum optical path length, the counter may be
 *          incremented by more than one to account for splitting of the
 *          segment into segments of allowed optical path length.
 * @param length segment length
 * @param mat Material associated with the segment
 * @param fsr_id the FSR ID of the FSR associated with the segment
 */
void CounterKernel::execute(FP_PRECISION length, Material* mat, long fsr_id,
                            int track_idx, int cmfd_surface_fwd,
                            int cmfd_surface_bwd, FP_PRECISION x_start,
                            FP_PRECISION y_start, FP_PRECISION z_start,
                            FP_PRECISION phi, FP_PRECISION theta) {

  /* Determine the number of cuts on the segment */
  FP_PRECISION sin_theta = sin(theta);
  FP_PRECISION max_sigma_t = mat->getMaxSigmaT();

  int num_cuts = 1;
  if (length * max_sigma_t * sin_theta > _max_tau)
    num_cuts = length * max_sigma_t * sin_theta / _max_tau + 1;

  /* Increment segment count */
  _count += num_cuts;
}


/**
 * @brief Writes segment information to the segmentation data array.
 * @details The SegmentationKernel execute function writes segment information
 *          to the segmentation data referenced by _segments. Due to
 *          restrictions on maximum optical path length, the counter may be
 *          incremented by more than one to account for splitting of the
 *          segment into segments of allowed optical path length.
 * @param length segment length
 * @param mat Material associated with the segment
 * @param fsr_id the FSR ID of the FSR associated with the segment
 * @param track_idx the track index in stack
 * @param cmfd_surface_fwd CMFD surface at the end of the segment in the forward
 *        direction
 * @param cmfd_surface_bwd CMFD surface at the end of the segment in the
 *        backward direction
 * @param x_start x coordinate of the start of the segment
 * @param y_start y coordinate of the start of the segment
 * @param z_start z coordinate of the start of the sement
 * @param phi azimuthal angle of this segment
 * @param theta polar angle of this segment
 */
void SegmentationKernel::execute(FP_PRECISION length, Material* mat, long fsr_id,
                                int track_idx, int cmfd_surface_fwd,
                                int cmfd_surface_bwd, FP_PRECISION x_start,
                                FP_PRECISION y_start, FP_PRECISION z_start,
                                FP_PRECISION phi, FP_PRECISION theta) {

  /* Check if segments have not been set, if so return */
  if (_segments == NULL)
    return;

  /* Determine the number of cuts on the segment */
  FP_PRECISION sin_theta = sin(theta);
  FP_PRECISION max_sigma_t = mat->getMaxSigmaT();

  int num_cuts = 1;
  if (length * max_sigma_t * sin_theta > _max_tau)
    num_cuts = length * max_sigma_t * sin_theta / _max_tau + 1;
  FP_PRECISION temp_length = _max_tau / (max_sigma_t * sin_theta);

  /* Add segment information */
  for (int i=0; i < num_cuts-1; i++) {
    _segments[_count]._length = temp_length;
    _segments[_count]._material = mat;
    _segments[_count]._region_id = fsr_id;
    _segments[_count]._track_idx = track_idx;
    _segments[_count]._starting_position[0] = x_start;
    _segments[_count]._starting_position[1] = y_start;
    _segments[_count]._starting_position[2] = z_start;
    _segments[_count]._cmfd_surface_fwd = -1;
    if (i == 0)
      _segments[_count]._cmfd_surface_bwd = cmfd_surface_bwd;
    else
      _segments[_count]._cmfd_surface_bwd = -1;
    length -= temp_length;
    x_start += temp_length * sin_theta * cos(phi);
    y_start += temp_length * sin_theta * sin(phi);
    z_start += temp_length * cos(theta);
    _count++;
  }
  _segments[_count]._length = length;
  _segments[_count]._material = mat;
  _segments[_count]._region_id = fsr_id;
  _segments[_count]._track_idx = track_idx;
  _segments[_count]._starting_position[0] = x_start;
  _segments[_count]._starting_position[1] = y_start;
  _segments[_count]._starting_position[2] = z_start;
  _segments[_count]._cmfd_surface_fwd = cmfd_surface_fwd;
  if (num_cuts > 1)
    _segments[_count]._cmfd_surface_bwd = -1;
  else
    _segments[_count]._cmfd_surface_bwd = cmfd_surface_bwd;
  _count++;
}


/**
 * @brief Constructor for the TransportKernel.
 * @param track_generator the TrackGenerator used to pull relevant tracking
 *        data from
 */
TransportKernel::TransportKernel(TrackGenerator* track_generator)
                                 : MOCKernel(track_generator) {
  _track_generator = track_generator;
  _direction = true;
  _min_track_idx = 0;
  _max_track_idx = 0;
  _azim_index = 0;
  _xy_index = 0;
  _polar_index = 0;
  _track_id = 0;
}


/**
 * @brief Destructor for the TransportKernel.
 */
TransportKernel::~TransportKernel() {
}


/**
 * @brief Sets a pointer to the CPUSolver to enable use of transport functions.
 * @param cpu_solver pointer to the CPUSolver
 */
void TransportKernel::setCPUSolver(CPUSolver* cpu_solver) {
  _cpu_solver = cpu_solver;
}


/**
 * @brief Create a new track3D from an existing one.
 * @param track track to create the new track from
 */
void TransportKernel::newTrack(Track* track) {
  Track3D* track_3D = dynamic_cast<Track3D*>(track);
  _azim_index = track_3D->getAzimIndex();
  _xy_index = track_3D->getXYIndex();
  _polar_index = track_3D->getPolarIndex();
  _track_id = track_3D->getUid();
  _count = 0;
}


/**
 * @brief Sets the direction of the current track.
 * @param direction the direction of the track: true = Forward, false = Backward
 */
void TransportKernel::setDirection(bool direction) {
  _direction = direction;
}


/**
 * @brief Returns the direction of the current track.
 * @return _direction the direction of the track: true = Forward, false = Backward
 */
bool TransportKernel::getDirection() {
  return _direction;
}


/**
 * @brief Apply MOC equations, tally fluxes and tally CMFD currents.
 * @param length length of the segment
 * @param mat pointer to the Material around the segment
 * @param fsr_id id of the FSR the segment lies in
 * @param track_idx id of the track
 * @param cmfd_surface_fwd CMFD surface at the end of the segment in the forward
 *        direction
 * @param cmfd_surface_bwd CMFD surface at the end of the segment in the
 *        backward direction
 * @param x_start x coordinate of the start of the segment
 * @param y_start y coordinate of the start of the segment
 * @param z_start z coordinate of the start of the sement
 * @param phi azimuthal angle of this segment
 * @param theta polar angle of this segment
 */
void TransportKernel::execute(FP_PRECISION length, Material* mat, long fsr_id,
                              int track_idx, int cmfd_surface_fwd,
                              int cmfd_surface_bwd, FP_PRECISION x_start,
                              FP_PRECISION y_start, FP_PRECISION z_start,
                              FP_PRECISION phi, FP_PRECISION theta) {

  /* Update lower and upper bounds for this track */
  if (track_idx < _min_track_idx)
    _min_track_idx = track_idx;
  else if (track_idx > _max_track_idx)
    _max_track_idx = track_idx;

  /* Determine the number of cuts on the segment */
  FP_PRECISION sin_theta = sin(theta);
  FP_PRECISION max_sigma_t = mat->getMaxSigmaT();

  int num_cuts = 1;
  if (length * max_sigma_t * sin_theta > _max_tau)
    num_cuts = length * max_sigma_t * sin_theta / _max_tau + 1;

  /* Handle the length of possible cuts of segment */
  FP_PRECISION remain_length = length;
  FP_PRECISION segment_length = std::min(_max_tau / (sin_theta * max_sigma_t),
                                         length);

  /* Get the track flux */
  long curr_track_id = _track_id + track_idx;
  float* track_flux = _cpu_solver->getBoundaryFlux(curr_track_id,
                                                   _direction);

  /* Allocate a buffer to store flux contribution of all cuts in this fsr */
  FP_PRECISION fsr_flux[_num_groups] __attribute__ ((aligned (VEC_ALIGNMENT)));
  memset(fsr_flux, 0, _num_groups * sizeof(FP_PRECISION));

  /* Apply MOC equations to segments */
  for (int i=0; i < num_cuts; i++) {

    /* Copy data into segment */
    segment curr_segment;
    curr_segment._length = std::min(segment_length, remain_length);
    curr_segment._material = mat;
    curr_segment._region_id = fsr_id;

    /* Apply MOC equations */
    _cpu_solver->tallyScalarFlux(&curr_segment, _azim_index, fsr_flux,
                                 track_flux);

    /* CMFD surfaces can only be on the segment ends, and the propagation is
       always in the forward direction */
    if (i == num_cuts - 1) {
      curr_segment._cmfd_surface_fwd = cmfd_surface_fwd;
      _cpu_solver->tallyCurrent(&curr_segment, _azim_index, _polar_index,
                                track_flux, true);
    }

    /* Shorten remaining 3D length */
    remain_length -= segment_length;
  }

  /* Accumulate contribution of all cuts to FSR scalar flux */
  FP_PRECISION wgt = _track_generator->getQuadrature()->getWeightInline(
       _azim_index, _polar_index);
  _cpu_solver->accumulateScalarFluxContribution(fsr_id, wgt, fsr_flux);
}


/**
 * @brief Obtain and transfer the boundary track angular fluxes.
 */  //TODO Optimize, create a postTwoWay
void TransportKernel::post() {
#ifndef ONLYVACUUMBC

  TrackGenerator3D* track_generator_3D =
    dynamic_cast<TrackGenerator3D*>(_track_generator);

  /* Get track */
  Track3D track;
  TrackStackIndexes tsi;
  tsi._azim = _azim_index;
  tsi._xy = _xy_index;
  tsi._polar = _polar_index;

  for (int i=_min_track_idx; i <= _max_track_idx; i++) {
    float* track_flux = _cpu_solver->getBoundaryFlux(_track_id+i,
                                                     _direction);

    tsi._z = i;
    track_generator_3D->getTrackOTF(&track, &tsi);
    _cpu_solver->transferBoundaryFlux(&track, _azim_index, _polar_index,
                                      _direction, track_flux);
  }
#endif

  /* Reset track indexes */
  _min_track_idx = 0;
  _max_track_idx = 0;
}
