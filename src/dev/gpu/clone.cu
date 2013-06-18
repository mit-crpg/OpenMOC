#include "clone.h"


/**
 * @brief Given a pointer to a material on the host and a material on the 
 *        GPU, copy all of the properties from the material on the host 
 *        to the GPU.
 * @param material_h pointer to a material on the host
 * @param material_d pointer to a material on the GPU
 */
void cloneMaterialOnGPU(Material* material_h, dev_material* material_d) {

    /* Copy over the material's id and uid */
    int id = material_h->getId();
    int uid = material_h->getUid();
    int num_groups = material_h->getNumEnergyGroups();

    cudaMemcpy((void*)&material_d->_id, (void*)&id, sizeof(int), 
	       cudaMemcpyHostToDevice);
    cudaMemcpy((void*)&material_d->_uid, (void*)&uid, sizeof(int), 
	       cudaMemcpyHostToDevice);

    /* Allocate memory on the device for each material data array */
    double* sigma_t;
    double* sigma_a;
    double* sigma_s;
    double* sigma_f;
    double* nu_sigma_f;
    double* chi;

    /* Allocate memory on device for materials data arrays */
    cudaMalloc((void**)&sigma_t, num_groups * sizeof(double));
    cudaMalloc((void**)&sigma_a, num_groups * sizeof(double));
    cudaMalloc((void**)&sigma_s, num_groups * num_groups * sizeof(double));
    cudaMalloc((void**)&sigma_f, num_groups * sizeof(double));
    cudaMalloc((void**)&nu_sigma_f, num_groups * sizeof(double));
    cudaMalloc((void**)&chi, num_groups * sizeof(double));

    /* Copy materials data from host to arrays on the device */
    cudaMemcpy((void*)sigma_t, (void*)material_h->getSigmaT(), 
	       num_groups * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)sigma_a, (void*)material_h->getSigmaA(), 
	       num_groups * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)sigma_s, (void*)material_h->getSigmaS(), 
	       num_groups * num_groups * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)sigma_f, (void*)material_h->getSigmaF(), 
	       num_groups * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)nu_sigma_f, (void*)material_h->getNuSigmaF(), 
	       num_groups * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)chi, (void*)material_h->getChi(), 
	       num_groups * sizeof(double), cudaMemcpyHostToDevice);

    /* Copy materials data pointers to device material */
    cudaMemcpy((void*)&material_d->_sigma_t, (void*)&sigma_t, sizeof(double*), 
                cudaMemcpyHostToDevice);
    cudaMemcpy((void*)&material_d->_sigma_a, (void*)&sigma_a, sizeof(double*), 
                cudaMemcpyHostToDevice);
    cudaMemcpy((void*)&material_d->_sigma_s, (void*)&sigma_s, sizeof(double*), 
                cudaMemcpyHostToDevice);
    cudaMemcpy((void*)&material_d->_sigma_f, (void*)&sigma_f, sizeof(double*), 
                cudaMemcpyHostToDevice);
    cudaMemcpy((void*)&material_d->_nu_sigma_f, (void*)&nu_sigma_f, 
                sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)&material_d->_chi, (void*)&chi, sizeof(double*), 
                cudaMemcpyHostToDevice);

    return;
}


/**
 * @brief Given a pointer to a track on the host and a track on the GPU,
 *        copy all of the properties and segments from the track on the host 
 *        to the GPU.
 * @param track_h pointer to a track on the host
 * @param track_d pointer to a track on the GPU
 */
void cloneTrackOnGPU(Track* track_h, dev_track* track_d) {

    dev_segment* dev_segments;
    dev_segment* host_segments = new dev_segment[track_h->getNumSegments()];
    dev_track new_track;

    new_track._uid = track_h->getUid();
    new_track._num_segments = track_h->getNumSegments();
    new_track._azim_angle_index = track_h->getAzimAngleIndex();

    new_track._refl_in = track_h->isReflIn();
    new_track._refl_out = track_h->isReflOut();
    new_track._bc_in = track_h->getBCIn();
    new_track._bc_out = track_h->getBCOut();

    cudaMalloc((void**)&dev_segments,
	       track_h->getNumSegments() * sizeof(dev_segment));
    new_track._segments = dev_segments;

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

    delete [] host_segments;

    return;
}
