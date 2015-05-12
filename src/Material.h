/**
 * @file Material.h
 * @brief
 * @date January 19, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef MATERIAL_H_
#define MATERIAL_H_

#ifdef __cplusplus
#include "Python.h"
#include "log.h"
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#endif

/** A negligible cross-section value to over-ride user-defined
 *  cross-sections very near zero (e.g., within (-1E-10, 1E-10)) */
#define ZERO_SIGMA_T 1E-10

#ifdef ICPC
/** Word-aligned memory allocation for Intel's compiler */
#define MM_FREE(array) _mm_free(array)

/** Word-aligned memory allocation for Intel's compiler */
#define MM_MALLOC(size,alignment) _mm_alloc(size, alignment)

#else

/** Word-aligned memory deallocation for GNU's compiler */
#define MM_FREE(array) free(array)

/** Word-aligned memory allocation for GNU's compiler */
#define MM_MALLOC(size,alignment) malloc(size)

#endif

/** Error threshold for determining how close the sum of \f$ \Sigma_a \f$
 *  and \f$ \Sigma_s \f$ must match that of \f$ \Sigma_t \f$ for each energy
 *  group
 */
#define SIGMA_T_THRESH 1E-3


int material_id();
void reset_material_id();


/**
 * @class Material Material.h "src/Material.h"
 * @brief The Material class represents a unique material and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class Material {

private:

  /** A user-defined ID for each Material created */
  int _id;

  /** A user-defined name for the Material */
  char* _name;

  /** The number of energy groups */
  int _num_groups;

  /** An array of the total cross-sections for each energy group */
  FP_PRECISION* _sigma_t;

  /** An array of the absorption cross-sections for each energy group */
  FP_PRECISION* _sigma_a;

  /** A 2D array of the scattering cross-section matrix. The first index is
   *  row number and second index is column number */
  FP_PRECISION* _sigma_s;

  /** An array of the fission cross-sections for each energy group */
  FP_PRECISION* _sigma_f;

  /** An array of the fission cross-sections multiplied by nu \f$ \nu \f$
   *  for each energy group */
  FP_PRECISION* _nu_sigma_f;

  /** An array of the chi \f$ \chi \f$ values for each energy group */
  FP_PRECISION* _chi;

  /** An array of the diffusion coefficients for each energy group */
  FP_PRECISION* _dif_coef;

  /** An array of the diffusion coefficients for each energy group */
  FP_PRECISION* _buckling;

  /** An array of the diffusion coefficient for each energy group
   *  for each surface */
  FP_PRECISION* _dif_hat;

  /** An array of the CMFD correction to the diffusion coefficient values
   *  for each energy group for each surface */
  FP_PRECISION* _dif_tilde;

  /** A boolean representing whether or not this Material contains a non-zero
   *  fission cross-section and is fissionable */
  bool _fissionable;

  /** A boolean to indicate whether or not the data has been
   * allocated to be vector aligned for SIMD instructions */
  bool _data_aligned;

  /** The number of vector widths needed to fit all energy groups */
  int _num_vector_groups;

public:
  Material(int id=0, const char* name="");
  virtual ~Material();

  int getId() const;
  char* getName() const;
  int getNumEnergyGroups() const;
  FP_PRECISION* getSigmaT();
  FP_PRECISION* getSigmaA();
  FP_PRECISION* getSigmaS();
  FP_PRECISION* getSigmaF();
  FP_PRECISION* getNuSigmaF();
  FP_PRECISION* getChi();
  FP_PRECISION* getDifCoef();
  FP_PRECISION* getBuckling();
  FP_PRECISION* getDifHat();
  FP_PRECISION* getDifTilde();
  FP_PRECISION getSigmaTByGroup(int group);
  FP_PRECISION getSigmaAByGroup(int group);
  FP_PRECISION getSigmaSByGroup(int origin, int destination);
  FP_PRECISION getSigmaSByGroupInline(int origin, int destination);
  FP_PRECISION getSigmaFByGroup(int group);
  FP_PRECISION getNuSigmaFByGroup(int group);
  FP_PRECISION getChiByGroup(int group);
  FP_PRECISION getDifCoefByGroup(int group);
  FP_PRECISION getBucklingByGroup(int group);
  FP_PRECISION getDifHatByGroup(int group, int surface);
  FP_PRECISION getDifTildeByGroup(int group);  
  bool isFissionable();
  bool isDataAligned();
  int getNumVectorGroups();

  void setName(const char* name);
  void setNumEnergyGroups(const int num_groups);

  void setSigmaT(double* xs, int num_groups);
  void setSigmaA(double* xs, int num_groups);
  void setSigmaS(double* xs, int num_groups);
  void setSigmaF(double* xs, int num_groups);
  void setNuSigmaF(double* xs, int num_groups);
  void setChi(double* xs, int num_groups);
  void setBuckling(double* xs, int num_groups);
  void setDifCoef(double* xs, int num_groups);
  void setDifHat(double* xs, int num_groups);
  void setDifTilde(double* xs, int num_groups);

  void setSigmaTByGroup(double xs, int group);
  void setSigmaAByGroup(double xs, int group);
  void setSigmaFByGroup(double xs, int group);
  void setNuSigmaFByGroup(double xs, int group);
  void setSigmaSByGroup(double xs, int origin, int destination);
  void setChiByGroup(double xs, int group);
  void setBucklingByGroup(double xs, int group);
  void setDifCoefByGroup(double xs, int group);
  void setDifHatByGroup(double xs, int group, int surface);
  void setDifTildeByGroup(double xs, int group, int surface);

  void checkSigmaT();
  std::string toString();
  void printString();

  void alignData();

  Material* clone();
};


/**
 * @brief inline function for efficient mapping for scattering, from
 *        1D as stored in memory to 2D matrix
 * @details Encapsulates the logic for indexing into the scattering
 *        matrix so it does not need to be repeated in other parts of 
 *        the code.  Note that this routine is 0-based, rather than 
 *        1-based indexing, as it is intended for use inside the code,
 *        not by users from Python.
 * @param origin the column index of the matrix element
 * @param destination the row index of the matrix element
 */
inline FP_PRECISION Material::getSigmaSByGroupInline(
          int origin, int destination) {
  return _sigma_s[destination*_num_groups + origin];
}

#endif /* MATERIAL_H_ */
