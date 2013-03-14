/*
 * DeviceMaterial.h
 *
 *  Created on: Jan 19, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef DEVICEMATERIAL_H_
#define DEVICEMATERIAL_H_

#include <cutil.h>
#include "../host/Material.h"
#include "../host/configurations.h"


typedef struct dev_material {
	int _uid;
	FP_PRECISION _sigma_t[NUM_ENERGY_GROUPS];
	FP_PRECISION _sigma_a[NUM_ENERGY_GROUPS];
	FP_PRECISION _sigma_f[NUM_ENERGY_GROUPS];
	FP_PRECISION _nu_sigma_f[NUM_ENERGY_GROUPS];
	FP_PRECISION _chi[NUM_ENERGY_GROUPS];

	/* first index is row number; second index is column number */
	FP_PRECISION _sigma_s[NUM_ENERGY_GROUPS][NUM_ENERGY_GROUPS];
} dev_material;

void cloneOnDevice(Material* hmaterial, dev_material* dmaterial);

#endif /* MATERIAL_H_ */
