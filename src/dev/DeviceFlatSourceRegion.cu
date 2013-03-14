/*
 * DeviceFlatSourceRegion.cu
 *
 *  Created on: Feb 3, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */


#include "DeviceFlatSourceRegion.h"


/**
 * Given a pointer to a flatsourceregion on the host and a
 * flatsourceregion on the device, copy all of the properties
 * from the flatsourceregion on the host to the device
 * @param hfsr pointer to a flatsourceregion on the host
 * @param dfsr pointer to a flatsourceregion on the device
 */
void cloneOnDevice(FlatSourceRegion* hfsr, dev_flatsourceregion* dfsr) {


	/* Create a temporary flatsourceregion struct on the host and populate it
	 * with the data from the FlatSourceRegion object on the host */
    dev_flatsourceregion* temp =
    		(dev_flatsourceregion*)malloc(sizeof(dev_flatsourceregion));
	temp->_uid = hfsr->getUid();
	temp->_material_uid = hfsr->getMaterial()->getUid();
	temp->_volume = hfsr->getVolume();

	for (int i=0; i < NUM_ENERGY_GROUPS; i++) {
	  temp->_flux[i] = hfsr->getFlux()[i];
	  temp->_new_source[i] = hfsr->getSource()[i];
	  temp->_ratios[i] = hfsr->getRatios()[i];
	}

	CUDA_SAFE_CALL(cudaMemcpy((void*)dfsr, (void*)temp,
			sizeof(dev_flatsourceregion), cudaMemcpyHostToDevice));

       free(temp);

	return;
}


/**
 * Compute the volumetric fission rate in this flat source region by adding
 * up the fission rates in each energy group. This mehtod assumes that fixed
 * source iteration has already been run since it uses the flux stored in
 * this region
 */
FP_PRECISION computeFissionRate(dev_flatsourceregion* dev_fsr,
													Material* material) {

	FP_PRECISION power = 0.0;
	FP_PRECISION* sigma_f = material->getSigmaF();

	/* Add the fission rates from each energy group */
	for (int e=0; e < NUM_ENERGY_GROUPS; e++)
		power += sigma_f[e] * dev_fsr->_flux[e];

	/* Multiply by volume of FSR */
	power *= dev_fsr->_volume;

	return power;
}
