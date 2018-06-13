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

  /** Index of the Track which reflects out of this Track along its "forward"
   * direction for reflective boundary conditions. */
  int _track_in;

  /** Index of the Track which reflects out of this Track along its "reverse"
   * direction for reflective boundary conditions. */
  int _track_out;

  /** The first index into the global 2D ragged array of Tracks for the Track
   *  that reflects out of this Track along its "forward" direction for
   *  reflective boundary conditions. */
  bool _refl_in;

  /** A boolean to indicate whether to give the flux to the "forward"
   *  (false) or "reverse" (true) direction of the Track reflecting out of
   *  this one along its "forward" direction for reflective boundary
   *  conditions. */
  bool _refl_out;

  /** A boolean to indicate whether the outgoing angular flux along this
   *  Track's "forward" direction should be zeroed out for vacuum boundary
   *  conditions. */
  bool _bc_in;

  /** A boolean to indicate whether the outgoing angular flux along this
   *  Track's "reverse" direction should be zeroed out for vacuum boundary
   *  conditions. */
  bool _bc_out;
};


#endif /* DEVICETRACK_H_ */
