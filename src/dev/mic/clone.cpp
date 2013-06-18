#include "clone.h"


/**
 * @brief Given a pointer to a material on the host and a material on the 
 *        MIC, copy all of the properties from the material on the host 
 *        to the MIC.
 * @param material_h pointer to a material on the host
 * @param material_d pointer to a material on the MIC
 */
void cloneMaterialOnMIC(Material* material_h, dev_material* material_d) {

    /* Copy over the material's id and uid */
    material_d->_id = material_h->getId();
    material_d->_uid = material_h->getUid();
    int num_groups = material_h->getNumEnergyGroups();

    /* Allocate memory for each material data array */
    material_d->_sigma_t = (double*)malloc(num_groups * sizeof(double));
    material_d->_sigma_a = (double*)malloc(num_groups * sizeof(double));
    material_d->_sigma_s=(double*)malloc(num_groups*num_groups*sizeof(double));
    material_d->_sigma_f = (double*)malloc(num_groups * sizeof(double));
    material_d->_nu_sigma_f = (double*)malloc(num_groups * sizeof(double));
    material_d->_chi = (double*)malloc(num_groups * sizeof(double));

    /* Copy materials data from host materials to device materials */
    memcpy((void*)material_d->_sigma_t, (void*)material_h->getSigmaT(), 
	   num_groups * sizeof(double));
    memcpy((void*)material_d->_sigma_a, (void*)material_h->getSigmaA(), 
	   num_groups * sizeof(double));
    memcpy((void*)material_d->_sigma_s, (void*)material_h->getSigmaS(), 
	   num_groups * num_groups * sizeof(double));
    memcpy((void*)material_d->_sigma_f, (void*)material_h->getSigmaF(), 
	   num_groups * sizeof(double));
    memcpy((void*)material_d->_nu_sigma_f, (void*)material_h->getNuSigmaF(), 
	   num_groups * sizeof(double));
    memcpy((void*)material_d->_chi, (void*)material_h->getChi(), 
	   num_groups * sizeof(double));

    return;
}


/**
 * @brief Given a pointer to a track on the host and a track on the MIC,
 *        copy all of the properties and segments from the track on the host 
 *        to the MIC.
 * @param track_h pointer to a track on the host
 * @param track_d pointer to a track on the MIC
 */
void cloneTrackOnMIC(Track* track_h, dev_track* track_d) {

    track_d->_uid = track_h->getUid();
    track_d->_num_segments = track_h->getNumSegments();
    track_d->_azim_angle_index = track_h->getAzimAngleIndex();
    track_d->_refl_in = track_h->isReflIn();
    track_d->_refl_out = track_h->isReflOut();
    track_d->_bc_in = track_h->getBCIn();
    track_d->_bc_out = track_h->getBCOut();

    track_d->_segments = new dev_segment[track_h->getNumSegments()];

    for (int s=0; s < track_h->getNumSegments(); s++) {
        segment* curr = track_h->getSegment(s);
	    track_d->_segments[s]._length = curr->_length;
	    track_d->_segments[s]._region_uid = curr->_region_id;
	    track_d->_segments[s]._material_uid = curr->_material->getUid();
    }

    return;
}
