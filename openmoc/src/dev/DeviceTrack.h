/**
 * @file DeviceTrack.h
 * @brief Structures for tracks and segments on a GPU.
 * @date June 29, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef DEVICETRACK_H_
#define DEVICETRACK_H_


#ifdef __cplusplus
#include "../host/Track.h"
#include "../host/configurations.h"
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
    /** A pointer to the material in which this segment resides */
    int _material_uid;
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
    int _azim_angle_index;
    dev_segment* _segments;
    int _num_segments;
    int _track_in, _track_out;
    bool _refl_in, _refl_out;
    bool _bc_in, _bc_out;
};


void cloneTrack(Track* track_h, dev_track* track_d);


#endif /* DEVICETRACK_H_ */
