#include "MOCKernel.h"
#include "TrackGenerator.h"

/**
 * @brief Constructor for the MOCKernel
 * @details The count for the number of segments traversed is set to zero and
 *          the maximum optical path length is retrieved from the
 *          TrackGenerator
 * @param track_generator the TrackGenerator used to pull relevant tracking
 *        data
 */
MOCKernel::MOCKernel(TrackGenerator* track_generator) {
  _count = 0;
  _max_tau = track_generator->retrieveMaxOpticalLength();
}


/**
 * @brief Constructor for the VolumeKernel assigns default values, calls
 *        the MOCKernel constructor, and pulls references to FSR locks and FSR
 *        volumes from the provided TrackGenerator.
 * @param track_generator the TrackGenerator used to pull relevant tracking
 *        data
 */
VolumeKernel::VolumeKernel(TrackGenerator* track_generator) :
                           MOCKernel(track_generator) {
  _FSR_locks = track_generator->getFSRLocks();

  if (_FSR_locks == NULL)
    log_printf(ERROR, "Unable to create a VolumeKernel without first creating "
               "FSR locks which are normally created during Track generation");

  _FSR_volumes = track_generator->getFSRVolumes();
  _quadrature = track_generator->getQuadrature();
  _weight = 0;
}


/**
 * @brief Constructor for the CounterKernel assigns default values and calls
 *        the MOCKernel constructor
 * @param track_generator the TrackGenerator used to pull relevant tracking
 *        data
 */
CounterKernel::CounterKernel(TrackGenerator* track_generator) :
                             MOCKernel(track_generator) {}


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
  int azim_index = track->getAzimAngleIndex();
  _weight = _quadrature->getAzimSpacing(azim_index)
      * _quadrature->getAzimWeight(azim_index);

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
 * @return the number of segments traversed
 */
int MOCKernel::getCount() {
  return _count;
}


/* @brief Resets the maximum optcal path length for a segment
 * @details MOC kernels ensure that there are no segments with an optical path
 *          length greater than the maximum optical path length by splitting
 *          them when they get too large.
 * @param max_tau the maximum optical path length for a segment
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
void VolumeKernel::execute(FP_PRECISION length, Material* mat, int id,
                           int cmfd_surface_fwd, int cmfd_surface_bwd) {

  /* Set omp lock for FSRs */
  omp_set_lock(&_FSR_locks[id]);

  /* Add value to buffer */
  _FSR_volumes[id] += _weight * length;

  /* Unset lock */
  omp_unset_lock(&_FSR_locks[id]);

  /* Increment count */
  _count++;
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
void CounterKernel::execute(FP_PRECISION length, Material* mat, int id,
                            int cmfd_surface_fwd, int cmfd_surface_bwd) {

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
