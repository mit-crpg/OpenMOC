#include "MOCKernel.h"

/**
 * @biief Constructor for the MOCKernel assigns default values
 */
MOCKernel::MOCKernel() {
  _count = 0;
  _max_tau = std::numeric_limits<FP_PRECISION>::max();
  _weight = 0;
}


/**
 * @biief Constructor for the MOCKernel assigns default values
 */
VolumeKernel::VolumeKernel(omp_lock_t* FSR_locks) : MOCKernel() {

  if (FSR_locks == NULL)
    log_printf(ERROR, "Unable to create a VolumeKernel without "
               "an array of FSR locks");

  _FSR_locks = FSR_locks;
}


/**
 * @brief Destructor for MOCKernel
 */
MOCKernel::~MOCKernel() {};


/**
 * @brief Destructor for MOCKernel
 */
VolumeKernel::~VolumeKernel() {};



/**
 * @brief Sets the location of segment data
 * @param segments pointer to segment data
 */
void MOCKernel::setSegments(segment* segments) {
  _segments = segments;
}


/**
 * @brief Resets the counter to zero
 */
void MOCKernel::resetCount() {
  _count = 0;
}


/**
 * @brief Sets the location of buffer data
 * @param buffer pointer to buffer data
 */
void MOCKernel::setBuffer(FP_PRECISION* buffer) {
  _buffer = buffer;
}


/**
 * @brief Sets the weight to apply to segment data
 * @param weight the associated track's weight
 */
void MOCKernel::setWeight(FP_PRECISION weight) {
  _weight = weight;
}

/**
 * @brief Sets the maximum allowed optical path length before segments are
 *        split
 * @param max_tau maximum optical path length
 */
void MOCKernel::setMaxVal(FP_PRECISION max_tau) {
  _max_tau = max_tau;
}


/*
 * @brief Reads and returns the current count
 * @details MOC kernels count how many times they are accessed. This value
 *          returns the value of the counter (number of execute accesses)
 *          since kernel creation or last reset.
 * @return _count the counter value
 */
int MOCKernel::getCount() {
  return _count;
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

  /* Set omp lock for corresponding block of 10 FSRs */
  omp_set_lock(&_FSR_locks[id]);

  /* Add value to buffer */
  _buffer[id] += _weight * length;

  /* Unset lock */
  omp_unset_lock(&_FSR_locks[id]);
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
void SegmentationKernel::execute(FP_PRECISION length, Material* mat, int id,
    int cmfd_surface_fwd, int cmfd_surface_bwd) {

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
    _segments[_count]._region_id = id;
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
  _segments[_count]._region_id = id;
  _segments[_count]._cmfd_surface_fwd = cmfd_surface_fwd;
  if (num_cuts > 1)
    _segments[_count]._cmfd_surface_bwd = -1;
  else
    _segments[_count]._cmfd_surface_bwd = cmfd_surface_bwd;
  _count++;
}
