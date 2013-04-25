#include "DeviceMaterial.h"


/**
 * @brief Given a pointer to a material on the host and a material on the 
 *        device, copy all of the properties from the material on the host 
 *        to the device.
 * @param host_material pointer to a material on the host
 * @param dev_material pointer to a material on the device
 */
void cloneOnDevice(Material* host_material, dev_material* dev_material) {

    /* Copy over the material's id and uid */
    int id = host_material->getId();
    int uid = host_material->getUid();

    cudaMemcpy((void*)&dev_material->_id, (void*)&id, sizeof(int), 
	       cudaMemcpyHostToDevice);
    cudaMemcpy((void*)&dev_material->_uid, (void*)&uid, sizeof(int), 
	       cudaMemcpyHostToDevice);

    /* Allocate memory on the device for each material data array */
    int num_groups = host_material->getNumEnergyGroups();

    cudaMalloc((void**)dev_material->_sigma_t, 
	       num_groups * sizeof(FP_PRECISION));
    cudaMalloc((void**)dev_material->_sigma_a, 
	       num_groups * sizeof(FP_PRECISION));
    cudaMalloc((void**)dev_material->_sigma_s, 
	       num_groups * num_groups * sizeof(FP_PRECISION));
    cudaMalloc((void**)dev_material->_sigma_f, 
	       num_groups * sizeof(FP_PRECISION));
    cudaMalloc((void**)dev_material->_nu_sigma_f, 
	       num_groups * sizeof(FP_PRECISION));
    cudaMalloc((void**)dev_material->_chi, 
	       num_groups * sizeof(FP_PRECISION));

    /* Copy materials data in to device material arrays */
    cudaMemcpy((void*)dev_material->_sigma_t, (void*)host_material->getSigmaT(),
	   num_groups * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)dev_material->_sigma_a, (void*)host_material->getSigmaA(),
	   num_groups * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)dev_material->_sigma_s, 
	       (void*)host_material->getSigmaS(),
	       num_groups * num_groups * sizeof(FP_PRECISION),
	       cudaMemcpyHostToDevice);
    cudaMemcpy((void*)dev_material->_sigma_f, (void*)host_material->getSigmaF(),
	   num_groups * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)dev_material->_nu_sigma_f, 
	       (void*)host_material->getNuSigmaF(),
	       num_groups * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)dev_material->_chi, (void*)host_material->getChi(),
	   num_groups * sizeof(FP_PRECISION), cudaMemcpyHostToDevice);

    return;
}
