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
