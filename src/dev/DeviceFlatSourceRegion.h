/*
 * FlatSourceRegion.h
 *
 *  Created on: Feb 3, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef DEVICEFLATSOURCEREGION_H_
#define DEVICEFLATSOURCEREGION_H_


#include <cutil.h>
#include "../host/FlatSourceRegion.h"
#include "DeviceMaterial.h"
#include "DeviceTrack.h"


typedef struct dev_flatsourceregion {
	int _uid;
	int _material_uid;
    FP_PRECISION _volume;
	FP_PRECISION _flux[NUM_ENERGY_GROUPS];
	FP_PRECISION _new_source[NUM_ENERGY_GROUPS];
	FP_PRECISION _old_source[NUM_ENERGY_GROUPS];
    FP_PRECISION _ratios[NUM_ENERGY_GROUPS];
} dev_flatsourceregion;


void cloneOnDevice(FlatSourceRegion* hfsr, dev_flatsourceregion* dfsr);
FP_PRECISION computeFissionRate(dev_flatsourceregion* dev_fsr,
													Material* material);


#endif /* DEVICEFLATSOURCEREGION_H_ */
