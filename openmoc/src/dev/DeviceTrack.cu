#include "DeviceTrack.h"


/**
 * @brief Given a pointer to a track on the host and a track on the device,
 *        copy all of the properties and segments from the track on the host 
 *        to the device.
 * @param track_h pointer to a track on the host
 * @param track_d pointer to a track on the device
 */
void cloneTrack(Track* track_h, dev_track* track_d) {

    dev_segment* dev_segments;
    dev_segment* host_segments = new dev_segment[track_h->getNumSegments()];
    dev_track new_track;

    new_track._uid = track_h->getUid();
    new_track._num_segments = track_h->getNumSegments();
    new_track._azim_angle_index = track_h->getAzimAngleIndex();

    new_track._track_in_i = track_h->getTrackInI();
    new_track._track_in_j = track_h->getTrackInJ();
    new_track._track_out_i = track_h->getTrackOutI();
    new_track._track_out_j = track_h->getTrackOutJ();

    new_track._refl_in = track_h->isReflIn();
    new_track._refl_out = track_h->isReflOut();
    new_track._bc_in = track_h->getBCIn();
    new_track._bc_out = track_h->getBCOut();

    cudaMalloc((void**)&dev_segments,
	       track_h->getNumSegments() * sizeof(dev_segment));
    new_track._segments = dev_segments;

    /* Iterate over all segments and memcpy to device */
    for (int s=0; s < track_h->getNumSegments(); s++) {
        segment* curr = track_h->getSegment(s);
	host_segments[s]._length = curr->_length;
	host_segments[s]._region_uid = curr->_region_id;
	host_segments[s]._material_uid = curr->_material->getUid();
    }

    cudaMemcpy((void*)dev_segments, (void*)host_segments,
	       track_h->getNumSegments() * sizeof(dev_segment), 
	       cudaMemcpyHostToDevice);
    cudaMemcpy((void*)track_d, (void*)&new_track, sizeof(dev_track), 
	       cudaMemcpyHostToDevice);

    return;
}
