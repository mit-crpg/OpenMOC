/**
 * @file DeviceMaterial.h
 * @brief
 * @author July 19, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef DEVICEMATERIAL_H_
#define DEVICEMATERIAL_H_

#ifdef __cplusplus
//#include <cutil.h>
#include "../host/Material.h"
#include "../host/configurations.h"
#endif

typedef struct dev_material {
    int _uid;
    FP_PRECISION _sigma_t[7];
    FP_PRECISION _sigma_a[7];
    FP_PRECISION _sigma_f[7];
    FP_PRECISION _nu_sigma_f[7];
    FP_PRECISION _chi[7];

    /* first index is row number; second index is column number */
    FP_PRECISION _sigma_s[7][7];
} dev_material;


void cloneOnDevice(Material* hmaterial, dev_material* dmaterial);

#endif /* MATERIAL_H_ */
