#include "MOCKernel.h"

/**
 * @biief Constructor for the MOCKernel assigns default values
 */
MOCKernel::MOCKernel() {
  _count = 0;
  _max_tau = 1.79769e+308;
  _weight = 0;
}


/**
 * @brief Destructor for MOCKernel
 */
MOCKernel::~MOCKernel() {};


/**
 * @brief Sets the location of segment data
 * @details Sets pointer for writing segment data
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
 * @details Sets pointer for writing floating point data
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
 * @brief Sets the maximum allowed optical path length
 * @details Sets the maximum allowed optical path length before segments are
 *          split
 * @param max_tau maximum optical path length
 */
void MOCKernel::setMaxVal(FP_PRECISION max_tau) {
  _max_tau = max_tau;
}


/*
 * @brief Reads and returns the current count
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
void VolumeKernel::execute(FP_PRECISION length, Material* mat, int id) {
  
  /* Add value to buffer */
  _buffer[id] += _weight * length;
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
void CounterKernel::execute(FP_PRECISION length, Material* mat, int id) {
  
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
void SegmentationKernel::execute(FP_PRECISION length, Material* mat, int id) {
  
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
    length -= temp_length;
    _count++;
  }
  _segments[_count]._length = length;
  _segments[_count]._material = mat;
  _segments[_count]._region_id = id;
  _count++;
}
