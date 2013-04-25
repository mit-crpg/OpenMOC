/**
 * @file FlatSourceRegion.h
 * @brief
 * @date July 10, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef DEVICEFLATSOURCEREGION_H_
#define DEVICEFLATSOURCEREGION_H_


#ifdef __cplusplus
#include "../host/FlatSourceRegion.h"
#include "DeviceMaterial.h"
#include "DeviceTrack.h"
#endif

/**
 * @struct dev_flatsourceregion
 * @brief This struct represents a unique discretized region within the 
 *        geometry in which OpenMOC approximates the scalar flux and source
 *        as constant.
 * @details The dev_flatsourceregion struct is intended to be stored on a GPU.
 */
typedef struct dev_flatsourceregion {
    int _uid;
    int _material_uid;
    FP_PRECISION _volume;
} dev_flatsourceregion;


void cloneOnDevice(FlatSourceRegion* fsr_h, dev_flatsourceregion* fsr_d);


#endif /* DEVICEFLATSOURCEREGION_H_ */
