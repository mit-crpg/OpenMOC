#include "DeviceMaterial.h"


/**
 * @brief Given a pointer to a material on the host and a material on the 
 *        device, copy all of the properties from the material on the host 
 *        to the device.
 * @param material_h pointer to a material on the host
 * @param material_d pointer to a material on the device
 */
void cloneOnDevice(Material* material_h, dev_material* material_d) {

    /* Copy over the material's id and uid */
    int id = material_h->getId();
    int uid = material_h->getUid();

    cudaMemcpy((void*)&material_d->_id, (void*)&id, sizeof(int), 
	       cudaMemcpyHostToDevice);
    cudaMemcpy((void*)&material_d->_uid, (void*)&uid, sizeof(int), 
	       cudaMemcpyHostToDevice);

    /* Allocate memory on the device for each material data array */
    int num_groups = material_h->getNumEnergyGroups();

    cudaMalloc((void**)material_d->_sigma_t, 
	       num_groups * sizeof(double));
    cudaMalloc((void**)material_d->_sigma_a, 
	       num_groups * sizeof(FP_PRECISION));
    cudaMalloc((void**)material_d->_sigma_s, 
	       num_groups * num_groups * sizeof(double));
    cudaMalloc((void**)material_d->_sigma_f, 
	       num_groups * sizeof(double));
    cudaMalloc((void**)material_d->_nu_sigma_f, 
	       num_groups * sizeof(double));
    cudaMalloc((void**)material_d->_chi, 
	       num_groups * sizeof(double));

    /* Copy materials data in to device material arrays */
    cudaMemcpy((void*)material_d->_sigma_t, (void*)material_h->getSigmaT(),
	   num_groups * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)material_d->_sigma_a, (void*)material_h->getSigmaA(),
	   num_groups * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)material_d->_sigma_s, 
	       (void*)material_h->getSigmaS(),
	       num_groups * num_groups * sizeof(double),
	       cudaMemcpyHostToDevice);
    cudaMemcpy((void*)material_d->_sigma_f, (void*)material_h->getSigmaF(),
	   num_groups * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)material_d->_nu_sigma_f, 
	       (void*)material_h->getNuSigmaF(),
	       num_groups * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)material_d->_chi, (void*)material_h->getChi(),
	   num_groups * sizeof(double), cudaMemcpyHostToDevice);

    return;
}
