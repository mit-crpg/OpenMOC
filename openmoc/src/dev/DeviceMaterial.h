/**
 * @file DeviceMaterial.h
 * @brief The struct of material's nuclear data to be stored on a GPU.
 * @author July 19, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef DEVICEMATERIAL_H_
#define DEVICEMATERIAL_H_

#ifdef __cplusplus
#include "../host/Material.h"
#include "../host/configurations.h"
#endif

/**
 * @struct dev_material
 * @brief A material's nuclear data to be stored on a GPU.
 */
typedef struct dev_material {
    int _id;
    int _uid;
    FP_PRECISION* _sigma_t;
    FP_PRECISION* _sigma_a;
    FP_PRECISION* _sigma_f;
    FP_PRECISION* _nu_sigma_f;
    FP_PRECISION* _chi;

    /* first index is row number; second index is column number */
    FP_PRECISION* _sigma_s;

    dev_material() {
      _sigma_t = NULL;
      _sigma_a = NULL;
      _sigma_s = NULL;
      _sigma_f = NULL;
      _nu_sigma_f = NULL;
      _chi = NULL;
    }

    ~dev_material() {
        if (_sigma_t != NULL)
	    delete [] _sigma_t;
	if (_sigma_a != NULL)
	    delete [] _sigma_a;
        if (_sigma_s != NULL)
	    delete [] _sigma_s;
        if (_sigma_f != NULL)
	    delete [] _sigma_f;
        if (_nu_sigma_f != NULL)
	    delete [] _nu_sigma_f;
        if (_chi != NULL)
	    delete [] _chi;
    }

} dev_material;


void cloneOnDevice(Material* host_material, dev_material* dev_material);

#endif /* MATERIAL_H_ */
