/*
 * Material.cu
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */


#include "DeviceMaterial.h"

/**
 * Given a pointer to a material on the host and a material on the device,
 * copy all of the properties from the material on the host to the device
 * @param hmaterial pointer to a material on the host
 * @param dmaterial pointer to a material on the device
 */
void cloneOnDevice(Material* hmaterial, dev_material* dmaterial) {

	/* Create a temporary material struct on the host and populate it with
	 * the data from the Material objecton the host */
	dev_material temp;
	temp._uid = hmaterial->getUid();

	for (int i=0; i < NUM_ENERGY_GROUPS; i++) {
		temp._sigma_a[i] = hmaterial->getSigmaA()[i];
		temp._sigma_f[i] = hmaterial->getSigmaF()[i];
		temp._nu_sigma_f[i] = hmaterial->getNuSigmaF()[i];
		temp._sigma_t[i] = hmaterial->getSigmaT()[i];
		temp._chi[i] = hmaterial->getChi()[i];

		for (int j=0; j < NUM_ENERGY_GROUPS; j++)
			temp._sigma_s[i][j] = hmaterial->getSigmaS()[i*NUM_ENERGY_GROUPS+j];
	}

	CUDA_SAFE_CALL(cudaMemcpy((void*)dmaterial, (void*)&temp,
				sizeof(dev_material), cudaMemcpyHostToDevice));

	return;

}
