#include "clone.h"

/**
 * @brief Given a pointer to a Material on the host and a dev_material on the
 *        GPU, copy all of the properties from the Material object on the host
 *        struct to the GPU.
 * @details This routine is called by the GPUSolver::initializeMaterials()
 *          private class method and is not intended to be called directly.
 * @param material_h pointer to a Material on the host
 * @param material_d pointer to a dev_material on the GPU
 */
void clone_material(Material* material_h, dev_material* material_d) {

  /* Copy over the Material's ID */
  int id = material_h->getId();
  int num_groups = material_h->getNumEnergyGroups();

  cudaMemcpy(&material_d->_id, &id, sizeof(int),
             cudaMemcpyHostToDevice);
  getLastCudaError();

  /* Allocate memory on the device for each dev_material data array */
  FP_PRECISION* sigma_t;
  FP_PRECISION* sigma_s;
  FP_PRECISION* sigma_f;
  FP_PRECISION* nu_sigma_f;
  FP_PRECISION* chi;
  FP_PRECISION* fiss_matrix;

  /* Allocate memory on device for dev_material data arrays */
  cudaMalloc(&sigma_t, num_groups * sizeof(FP_PRECISION));
  cudaMalloc(&sigma_s, num_groups * num_groups * sizeof(FP_PRECISION));
  cudaMalloc(&sigma_f, num_groups * sizeof(FP_PRECISION));
  cudaMalloc(&nu_sigma_f, num_groups * sizeof(FP_PRECISION));
  cudaMalloc(&chi, num_groups * sizeof(FP_PRECISION));
  cudaMalloc(&fiss_matrix, num_groups * num_groups * sizeof(FP_PRECISION));
  getLastCudaError();

  /* Copy Material data from host to arrays on the device */
  cudaMemcpy(sigma_t, material_h->getSigmaT(),
             num_groups * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);
  cudaMemcpy(sigma_s, material_h->getSigmaS(),
             num_groups * num_groups * sizeof(FP_PRECISION),
             cudaMemcpyHostToDevice);
  cudaMemcpy(sigma_f, material_h->getSigmaF(),
             num_groups * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);
  cudaMemcpy(nu_sigma_f, material_h->getNuSigmaF(),
             num_groups * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);
  cudaMemcpy(chi, material_h->getChi(),
             num_groups * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);
  cudaMemcpy(fiss_matrix, material_h->getFissionMatrix(),
             num_groups * num_groups * sizeof(FP_PRECISION),
             cudaMemcpyHostToDevice);
  getLastCudaError();

  /* Copy Material data pointers to dev_material on GPU */
  cudaMemcpy(&material_d->_sigma_t, &sigma_t,
             sizeof(FP_PRECISION*), cudaMemcpyHostToDevice);
  cudaMemcpy(&material_d->_sigma_s, &sigma_s,
             sizeof(FP_PRECISION*), cudaMemcpyHostToDevice);
  cudaMemcpy(&material_d->_sigma_f, &sigma_f,
             sizeof(FP_PRECISION*), cudaMemcpyHostToDevice);
  cudaMemcpy(&material_d->_nu_sigma_f, &nu_sigma_f,
             sizeof(FP_PRECISION*), cudaMemcpyHostToDevice);
  cudaMemcpy(&material_d->_chi, &chi,
             sizeof(FP_PRECISION*), cudaMemcpyHostToDevice);
  cudaMemcpy(&material_d->_fiss_matrix, &fiss_matrix,
             sizeof(FP_PRECISION*), cudaMemcpyHostToDevice);
  getLastCudaError();
}


/**
 * @brief Given a pointer to a Track on the host, a dev_track on
 *        the GPU, and the map of material IDs to indices in the
 *        _materials array, copy all of the class attributes and
 *        segments from the Track object on the host to the GPU.
 * @details This routine is called by the GPUSolver::initializeTracks()
 *          private class method and is not intended to be called
 *          directly.
 * @param track_h pointer to a Track on the host
 * @param track_d pointer to a dev_track on the GPU
 * @param material_IDs_to_indices map of material IDs to indices
 *        in the _materials array.
 */
void clone_track(Track* track_h, dev_track* track_d,
     		 std::map<int, int> &material_IDs_to_indices) {

  dev_segment* dev_segments;
  dev_segment* host_segments = new dev_segment[track_h->getNumSegments()];
  dev_track new_track;

  new_track._uid = track_h->getUid();
  new_track._num_segments = track_h->getNumSegments();
  new_track._azim_angle_index = track_h->getAzimIndex();

  new_track._next_fwd_is_fwd = track_h->getNextFwdFwd();
  new_track._next_bwd_is_fwd = track_h->getNextBwdFwd();
  new_track._transfer_flux_fwd = track_h->getBCFwd();
  new_track._transfer_flux_bwd = track_h->getBCBwd();
  new_track._next_track_fwd = track_h->getTrackNextFwd();
  new_track._next_track_bwd = track_h->getTrackNextBwd();

  cudaMalloc(&dev_segments,
             track_h->getNumSegments() * sizeof(dev_segment));
  getLastCudaError();
  new_track._segments = dev_segments;

  for (int s=0; s < track_h->getNumSegments(); s++) {
    segment* curr = track_h->getSegment(s);
    host_segments[s]._length = curr->_length;
    host_segments[s]._region_uid = curr->_region_id;
    host_segments[s]._material_index =
      material_IDs_to_indices[curr->_material->getId()];
  }

  cudaMemcpy(dev_segments, host_segments,
             track_h->getNumSegments() * sizeof(dev_segment),
             cudaMemcpyHostToDevice);
  getLastCudaError();
  cudaMemcpy(track_d, &new_track, sizeof(dev_track),
             cudaMemcpyHostToDevice);
  getLastCudaError();

  delete [] host_segments;
}
