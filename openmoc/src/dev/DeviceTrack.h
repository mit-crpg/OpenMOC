/*
 * DeviceTrack.h
 *
 *  Created on: Jun 29, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */


#ifndef DEVICETRACK_H_
#define DEVICETRACK_H_


#include <cutil.h>
#include "../host/Track.h"
#include "../host/configurations.h"



/* Represent a segment along a given track on the device */
struct dev_segment {
	FP_PRECISION _length;
	int _material_uid;
	int _region_uid;
#if STORE_PREFACTORS
	FP_PRECISION _prefactors[NUM_ENERGY_GROUPS][NUM_POLAR_ANGLES];
#endif
};


/* Represent a track on the device */
struct dev_track {
	int _azim_angle_index;
	FP_PRECISION _polar_fluxes[2 * GRP_TIMES_ANG];
    dev_segment* _segments;
    int _num_segments;
    int _track_in, _track_out;
	bool _refl_in, _refl_out;
	bool _bc_in, _bc_out;
};


void cloneTrack(Track* track, dev_track* dev_track);



#endif /* DEVICETRACK_H_ */
