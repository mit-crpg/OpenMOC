/**
 * @file DeviceTrack.h
 * @brief Structures for Tracks and Track segments on a GPU.
 * @date June 29, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef DEVICETRACK_H_
#define DEVICETRACK_H_


#ifdef __cplusplus
#include "../Track.h"
#endif

/**
 * @struct dev_segment
 * @brief A dev_segment represents a line segment within a single flat source
 *        region along a track.
 * @details The dev_segment is intended for use on the GPU.
 */
struct dev_segment {

  /** The length of the segment (cm) */
  FP_PRECISION _length;

  /** An index into the _materials array that contains Material pointers */
  int _material_index;

  /** The ID for flat source region in which this segment resides */
  int _region_uid;
};


/**
 * @struct dev_track
 * @brief A dev_track represents a characteristic line across the geometry.
 * @details A dev_track has particular starting and ending points on the
 *          boundaries of the geometry and an azimuthal angle. The dev_track
 *          is intended for use on the GPU.
 */
struct dev_track {

  /** A monotonically increasing unique ID for each Track created */
  int _uid;

  /** The azimuthal angle index into the global 2D ragged array of Tracks */
  int _azim_angle_index;

  /** A vector of segments making up this track */
  dev_segment* _segments;

  /** The number of segments making up this Track */
  int _num_segments;

  /** Index of the next Track when traveling along this Track in the "forward"
   * direction. */
  long _next_track_fwd;

  /** Index of the next Track when traveling along this Track in the "reverse"
   * direction. */
  long _next_track_bwd;

  /** A boolean to indicate whether to give the flux to the "forward" (true)
   *  or "backward" (false) direction of the next Track going in the "forward"
   *  direction. */
  bool _next_fwd_is_fwd;

  /** A boolean to indicate whether to give the flux to the "forward" (true)
   *  or "reverse" (false) direction of the next Track going in the "reverse"
   *  direction. */
  bool _next_bwd_is_fwd;

  /** A boolean to indicate whether the outgoing angular flux along this
   *  Track's "forward" direction should be transferred to the outgoing
   *  Track. */
  bool _transfer_flux_fwd;

  /** A boolean to indicate whether the outgoing angular flux along this
   *  Track's "reverse" direction should be transferred to the incoming
   *  Track. */
  bool _transfer_flux_bwd;
};


#endif /* DEVICETRACK_H_ */
