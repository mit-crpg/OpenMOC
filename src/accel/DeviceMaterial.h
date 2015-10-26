/**
 * @file DeviceMaterial.h
 * @brief The struct of Material's nuclear data to be stored on a GPU.
 * @author July 19, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef DEVICEMATERIAL_H_
#define DEVICEMATERIAL_H_

#ifdef __cplusplus
#include "../log.h"
#include "../Material.h"
#endif


/**
 * @struct dev_material
 * @brief A Material's nuclear data to be stored on a GPU.
 */
struct dev_material {

  /** A user-defined ID for each Material created */
  int _id;

  /** An array of the total cross-sections for each energy group */
  FP_PRECISION* _sigma_t;

  /** A 2D array of the scattering cross-section matrix. The first index is
   *  row number and second index is column number */
  FP_PRECISION* _sigma_f;

  /** An array of the fission cross-sections multiplied by nu \f$ \nu \f$
   *  for each energy group */
  FP_PRECISION* _nu_sigma_f;

  /** An array of the chi \f$ \chi \f$ values for each energy group */
  FP_PRECISION* _chi;

  /** A 2D array of the fission matrix from/into each group */
  FP_PRECISION* _fiss_matrix;

  /** A 2D array of the scattering cross-section matrix from/into each group */
  FP_PRECISION* _sigma_s;

  /**
   * @brief Constructor for a dev_material struct on a GPU.
   */
  dev_material() {
    _sigma_t = NULL;
    _sigma_s = NULL;
    _sigma_f = NULL;
    _nu_sigma_f = NULL;
    _chi = NULL;
    _fiss_matrix = NULL;
  }

  /**
   * @brief Destructor releases data for all Material's cross-sections on GPU.
   */
  ~dev_material() {
    if (_sigma_t != NULL)
      delete [] _sigma_t;
    if (_sigma_s != NULL)
      delete [] _sigma_s;
    if (_sigma_f != NULL)
      delete [] _sigma_f;
    if (_nu_sigma_f != NULL)
      delete [] _nu_sigma_f;
    if (_chi != NULL)
      delete [] _chi;
    if (_fiss_matrix != NULL)
      delete [] _fiss_matrix;
  }

};


#endif /* MATERIAL_H_ */
