#include "MOCKernel.h"
#include "TrackGenerator3D.h"
#include "CPUSolver.h"

/**
 * @brief Constructor for the MOCKernel assigns default values
 * @param track_generator the TrackGenerator used to pull relevant tracking
 *        data from
 * @param row_num the row index into the temporary segments matrix
 */
MOCKernel::MOCKernel(TrackGenerator* track_generator, int row_num) {
  _count = 0;
  _max_tau = track_generator->retrieveMaxOpticalLength();
}


/**
 * @brief Constructor for the VolumeKernel assigns default values, calls
 *        the MOCKernel constructor, and pulls refernces to FSR locks and FSR
 *        volumes from the provided TrackGenerator.
 * @param track_generator the TrackGenerator used to pull relevant tracking
 *        data from
 * @param row_num the row index into the temporary segments matrix
 */
VolumeKernel::VolumeKernel(TrackGenerator* track_generator, int row_num) :
                           MOCKernel(track_generator, row_num) {

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
 * @param row_num the row index into the temporary segments matrix
 */
SegmentationKernel::SegmentationKernel(TrackGenerator* track_generator,
                                       int row_num)
                                       : MOCKernel(track_generator, row_num) {

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
 *        the MOCKernel constructor
 * @param track_generator the TrackGenerator used to pull relevant tracking
 *        data from
 * @param row_num the row index into the temporary segments matrix
 */
CounterKernel::CounterKernel(TrackGenerator* track_generator, int row_num) :
                           MOCKernel(track_generator, row_num) {}


/**
 * @brief Prepares an MOCKernel for a new Track
 * @details Resets the segment count
 * @param track The new Track the MOCKernel prepares to handle
 */
void MOCKernel::newTrack(Track* track) {
  _count = 0;
}


/**
 * @brief Prepares a VolumeKernel for a new Track
 * @details Resets the segment count and updates the weight for the new Track
 * @param track The new Track the MOCKernel prepares to handle
 */
void VolumeKernel::newTrack(Track* track) {

  /* Compute the Track cross-sectional area */
  int azim_index = track->getAzimIndex();
  _weight = _quadrature->getAzimSpacing(azim_index)
      * _quadrature->getAzimWeight(azim_index);

  Track3D* track_3D = dynamic_cast<Track3D*>(track);
  if (track_3D != NULL) {
    int polar_index = track_3D->getPolarIndex();
    _weight *= _quadrature->getPolarSpacing(azim_index, polar_index)
        * _quadrature->getPolarWeight(azim_index, polar_index);
  }

  /* Reset the count */
  _count = 0;
}


/**
 * @brief Destructor for MOCKernel
 */
MOCKernel::~MOCKernel() {};


/*
 * @brief Reads and returns the current count
 * @details MOC kernels count how many times they are accessed. This value
 *          returns the value of the counter (number of execute accesses)
 *          since kernel creation or last call to newTrack.
 * @return _count the counter value
 */
int MOCKernel::getCount() {
  return _count;
}


/* @brief Resets the maximum optcal path length for a segment
 * @details MOC kernels ensure that there are no segments with an optical path
 *          length greater than the maximum optical path length by splitting
 *          then when they get too large.
 * @param the maximum optical path length for a segment
 */
void MOCKernel::setMaxOpticalLength(FP_PRECISION max_tau) {
  _max_tau = max_tau;
}


/*
 * @brief Adds segment contribution to the FSR volume
 * @details The VolumeKernel execute function adds the product of the
 *          track length and track weight to the buffer array at index
 *          id, referring to the array of FSR volumes.
 * @param length segment length
 * @param mat Material associated with the segment
 * @param id the FSR ID of the FSR associated with the segment
 */
void VolumeKernel::execute(FP_PRECISION length, Material* mat, int fsr_id,
                           int track_idx, int cmfd_surface_fwd,
                           int cmfd_surface_bwd) {

  /* Set omp lock for FSRs */
  omp_set_lock(&_FSR_locks[fsr_id]);

  /* Add value to buffer */
  _FSR_volumes[fsr_id] += _weight * length;

  /* Unset lock */
  omp_unset_lock(&_FSR_locks[fsr_id]);

  /* Determine the number of cuts on the segment */
  FP_PRECISION* sigma_t = mat->getSigmaT();
  double max_sigma_t = 0;
  for (int e=0; e < mat->getNumEnergyGroups(); e++)
    if (sigma_t[e] > max_sigma_t)
      max_sigma_t = sigma_t[e];

  int num_cuts = std::max((int) std::ceil(length * max_sigma_t / _max_tau), 1);

  /* Increment count */
  _count += num_cuts;
}


/*
 * @brief Increments the counter for the number of segments on the track
 * @details The CounterKernel execute function counts the number of segments
 *          in a track by incrementing the counter variable upon execution. Due
 *          to restrictions on maximum optical path length, the counter may be
 *          incremented by more than one to account for splitting of the
 *          segment into segments of allowed optical path length.
 * @param length segment length
 * @param mat Material associated with the segment
 * @param id the FSR ID of the FSR associated with the segment
 */
void CounterKernel::execute(FP_PRECISION length, Material* mat, int fsr_id,
                            int track_idx, int cmfd_surface_fwd,
                            int cmfd_surface_bwd) {

  /* Determine the number of cuts on the segment */
  FP_PRECISION* sigma_t = mat->getSigmaT();
  double max_sigma_t = 0;
  for (int e=0; e < mat->getNumEnergyGroups(); e++)
    if (sigma_t[e] > max_sigma_t)
      max_sigma_t = sigma_t[e];

  int num_cuts = std::max((int) std::ceil(length * max_sigma_t / _max_tau), 1);

  /* Increment count */
  _count += num_cuts;
}


/*
 * @brief Writes segment information to the segmentation data array
 * @details The SegmentationKernel execute function writes segment information
 *          to the segmentation data referenced by _segments. Due to
 *          restrictions on maximum optical path length, the counter may be
 *          incremented by more than one to account for splitting of the
 *          segment into segments of allowed optical path length.
 * @param length segment length
 * @param mat Material associated with the segment
 * @param id the FSR ID of the FSR associated with the segment
 */
void SegmentationKernel::execute(FP_PRECISION length, Material* mat, int fsr_id,
                                int track_idx, int cmfd_surface_fwd,
                                int cmfd_surface_bwd) {

  /* Check if segments have not been set, if so return */
  if (_segments == NULL)
    return;

  /* Determine the number of cuts on the segment */
  FP_PRECISION* sigma_t = mat->getSigmaT();
  double max_sigma_t = 0;
  for (int e=0; e < mat->getNumEnergyGroups(); e++)
    if (sigma_t[e] > max_sigma_t)
      max_sigma_t = sigma_t[e];

  int num_cuts = std::max((int) std::ceil(length * max_sigma_t / _max_tau), 1);

  /* Add segment information */
  for (int i=0; i < num_cuts-1; i++) {
    FP_PRECISION temp_length = _max_tau / max_sigma_t;
    _segments[_count]._length = temp_length;
    _segments[_count]._material = mat;
    _segments[_count]._region_id = fsr_id;
    _segments[_count]._track_idx = track_idx;
    _segments[_count]._cmfd_surface_fwd = -1;
    if (i == 0)
      _segments[_count]._cmfd_surface_bwd = cmfd_surface_bwd;
    else
      _segments[_count]._cmfd_surface_bwd = -1;
    length -= temp_length;
    _count++;
  }
  _segments[_count]._length = length;
  _segments[_count]._material = mat;
  _segments[_count]._region_id = fsr_id;
  _segments[_count]._track_idx = track_idx;
  _segments[_count]._cmfd_surface_fwd = cmfd_surface_fwd;
  if (num_cuts > 1)
    _segments[_count]._cmfd_surface_bwd = -1;
  else
    _segments[_count]._cmfd_surface_bwd = cmfd_surface_bwd;
  _count++;
}


//FIXME
TransportKernel::TransportKernel(TrackGenerator* track_generator, int row_num)
                                 : MOCKernel(track_generator, row_num) {
  _direction = true;
  _min_track_idx = 0;
  _max_track_idx = 0;
  _azim_index = 0;
  _polar_index = 0;
  _track_id = 0;
  int num_groups = track_generator->getGeometry()->getNumEnergyGroups();
  _thread_fsr_flux = new FP_PRECISION[num_groups];
}

TransportKernel::~TransportKernel() {
  delete [] _thread_fsr_flux;
}


/**
 * @brief Sets a pointer to the CPUSolver to enable use of transport functions
 * @param cpu_solver pointer to the CPUSolver
 */
//FIXME
void TransportKernel::setCPUSolver(CPUSolver* cpu_solver) {
  _cpu_solver = cpu_solver;
}


/**
 * @brief Sets the indexes of the current Track
 * @param axim_index the Track's azimuthal index
 * @param polar_index the Track's polar index
 */
//FIXME: delete?
void TransportKernel::newTrack(Track* track) {
  Track3D* track_3D = dynamic_cast<Track3D*>(track);
  _azim_index = track_3D->getAzimIndex();
  _polar_index = track_3D->getPolarIndex();
  _track_id = track_3D->getUid();
  _count = 0;
}


/**
 * @brief Sets the direction of the current track
 * @param direction the direction of the track: true = Forward, false = Backward
 */
//FIXME: delete?
void TransportKernel::setDirection(bool direction) {
  _direction = direction;
}


//FIXME document
void TransportKernel::execute(FP_PRECISION length, Material* mat, int fsr_id,
                              int track_idx, int cmfd_surface_fwd,
                              int cmfd_surface_bwd) {

  if (track_idx < _min_track_idx)
    _min_track_idx = track_idx;
  else if (track_idx > _max_track_idx)
    _max_track_idx = track_idx;

  /* Determine the number of cuts on the segment */
  FP_PRECISION* sigma_t = mat->getSigmaT();
  double max_sigma_t = 0;
  for (int e=0; e < mat->getNumEnergyGroups(); e++)
    if (sigma_t[e] > max_sigma_t)
      max_sigma_t = sigma_t[e];

  int num_cuts = std::max((int) std::ceil(length * max_sigma_t / _max_tau), 1);

  /* Determine common length */
  FP_PRECISION temp_length = std::min(_max_tau / max_sigma_t, length);

  /* Apply MOC equations to segments */
  for (int i=0; i < num_cuts; i++) {

    /* Copy data into segment */
    segment curr_segment;
    curr_segment._length = temp_length;
    curr_segment._material = mat;
    curr_segment._region_id = fsr_id;
    curr_segment._track_idx = track_idx; //FIXME

    /* Determine CMFD surfaces */
    if (i == 0)
      curr_segment._cmfd_surface_bwd = cmfd_surface_bwd;
    else
      curr_segment._cmfd_surface_bwd = -1;
    if (i == num_cuts - 1)
      curr_segment._cmfd_surface_fwd = cmfd_surface_fwd;
    else
      curr_segment._cmfd_surface_fwd = -1;

    /* Get the backward track flux */
    int curr_track_id = _track_id + track_idx;
    FP_PRECISION* track_flux = _cpu_solver->getBoundaryFlux(curr_track_id,
                                                            _direction);

    /* Apply MOC equations */
    _cpu_solver->tallyScalarFlux(&curr_segment, _azim_index, _polar_index,
                                 track_flux, _thread_fsr_flux);
    _cpu_solver->tallyCurrent(&curr_segment, _azim_index, _polar_index,
                                     track_flux, true);

    /* Shorten remaining 3D length */
    length -= temp_length;
  }
}


//FIXME
void TransportKernel::post() {
  for (int i=_min_track_idx; i <= _max_track_idx; i++) {
    FP_PRECISION* track_flux = _cpu_solver->getBoundaryFlux(_track_id+i,
                                                            _direction);
    _cpu_solver->transferBoundaryFlux(_track_id+i, _azim_index, _polar_index,
                                      _direction, track_flux);
  }
  _min_track_idx = 0;
  _max_track_idx = 0;
}
