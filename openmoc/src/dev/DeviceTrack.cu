/*
 * DeviceTrack.cu
 *
 *  Created on: Jun 29, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */


#include "DeviceTrack.h"


/**
 * Given a pointer to a track on the host and a track on the device,
 * copy all of the properties and segments from the track on the host to
 * the device
 * @param htrack pointer to a track on the host
 * @param dtrack pointer to a track on the device
 */
void cloneTrack(Track* htrack, dev_track* dtrack) {

	dev_segment* dev_segments;
	dev_segment* host_segments = new dev_segment[htrack->getNumSegments()];
	dev_track new_track;

	new_track._refl_in = htrack->isReflIn();
	new_track._refl_out = htrack->isReflOut();
	new_track._bc_in = htrack->getBCIn();
	new_track._bc_out = htrack->getBCOut();
	new_track._num_segments = htrack->getNumSegments();
	new_track._azim_angle_index = htrack->getAzimAngleIndex();

	memcpy(new_track._polar_fluxes, htrack->getPolarFluxes(),
							2*GRP_TIMES_ANG*sizeof(FP_PRECISION));

	CUDA_SAFE_CALL(cudaMalloc((void**)&dev_segments,
					htrack->getNumSegments()*sizeof(dev_segment)));
	new_track._segments = dev_segments;

	/* Iterate over all segments and memcpy to device */
	for (int s = 0; s < htrack->getNumSegments(); s++) {
	  
		segment* curr = htrack->getSegment(s);
		host_segments[s]._length = curr->_length;
		host_segments[s]._region_uid = curr->_region_id;
		host_segments[s]._material_uid = curr->_material->getUid();

#if STORE_PREFACTORS
	  	for (int i=0; i < NUM_ENERGY_GROUPS; i++) {
	  		for (int j=0; j < NUM_POLAR_ANGLES; j++)
	  			host_segments[s]._prefactors[i][j] = curr->_prefactors[i][j];
	  	}
#endif
	}

	CUDA_SAFE_CALL(cudaMemcpy((void*)dev_segments, (void*)host_segments,
		  htrack->getNumSegments()*sizeof(dev_segment), cudaMemcpyHostToDevice));
 	CUDA_SAFE_CALL(cudaMemcpy((void*)dtrack, (void*)&new_track,
 		  	  	  	  	  	  	  	  sizeof(dev_track), cudaMemcpyHostToDevice));

 	return;
}
