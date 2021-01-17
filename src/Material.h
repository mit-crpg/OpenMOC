/**
 * @file Material.h
 * @brief
 * @date January 19, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef MATERIAL_H_
#define MATERIAL_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "constants.h"
#include "log.h"
#include "linalg.h"
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#ifdef __linux__
#include <malloc.h>
#endif
#endif

#ifdef INTEL
/** Aligned memory allocation for Intel's compiler */
#define MM_MALLOC(size,alignment) _mm_malloc(size, alignment)

/** Aligned memory deallocation for Intel's compiler */
#define MM_FREE(array) _mm_free(array)

#else
#ifdef __linux__
/** Aligned memory allocation for GNU's compiler */
#define MM_MALLOC(size,alignment) memalign(alignment, size)
#else
/** macosx aligns on 16 byte boundaries */
#define MM_MALLOC(size,alignment) malloc(size)
#if VEC_ALIGNMENT>64
#error "VEC_ALIGNMENT should be set to 64 bits for macosx"
#endif
#endif

/** Aligned memory deallocation for GNU's compiler */
#define MM_FREE(array) free(array)
#endif


int material_id();
void reset_material_id();
void maximize_material_id(int material_id);


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

  /** The volume / area of the Material computed from overlapping segments */
  double _volume;

  /** The total number of instances of this Material in the Geometry */
  int _num_instances;

  /** The number of energy groups */
  int _num_groups;

  /** An array of the total cross-sections for each energy group */
  FP_PRECISION* _sigma_t;

  /** Max total cross section */
  FP_PRECISION _max_sigma_t;

  /** A 2D array of the scattering cross-section matrix from/into each group */
  FP_PRECISION* _sigma_s;

  /** An array of the absorption cross-sections for each energy group */
  FP_PRECISION* _sigma_a;

  /** An array of the fission cross-sections for each energy group */
  FP_PRECISION* _sigma_f;

  /** An array of the fission cross-sections multiplied by nu \f$ \nu \f$
   *  for each energy group */
  FP_PRECISION* _nu_sigma_f;

  /** An array of the chi \f$ \chi \f$ values for each energy group */
  FP_PRECISION* _chi;

  /** A 2D array of the fission matrix from/into each group */
  FP_PRECISION* _fiss_matrix;

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
  double getVolume();
  int getNumInstances();
  int getNumEnergyGroups() const;
  FP_PRECISION* getSigmaT();
  FP_PRECISION getMaxSigmaT();
  FP_PRECISION* getSigmaS();
  FP_PRECISION* getSigmaA();
  FP_PRECISION* getSigmaF();
  FP_PRECISION* getNuSigmaF();
  FP_PRECISION* getChi();
  FP_PRECISION* getFissionMatrix();
  FP_PRECISION getSigmaTByGroup(int group);
  FP_PRECISION getSigmaSByGroup(int origin, int destination);
  FP_PRECISION getSigmaAByGroup(int group);
  FP_PRECISION getSigmaFByGroup(int group);
  FP_PRECISION getNuSigmaFByGroup(int group);
  FP_PRECISION getChiByGroup(int group);
  FP_PRECISION getFissionMatrixByGroup(int origin, int destination);
  bool isFissionable();
  bool isDataAligned();
  int getNumVectorGroups();

  void setName(const char* name);
  void setVolume(double volume);
  void incrementVolume(double volume);
  void setNumInstances(int num_instances);
  void incrementNumInstances();
  void setNumEnergyGroups(const int num_groups);

  void setSigmaT(double* xs, int num_groups);
  void setSigmaS(double* xs, int num_groups);
  void setSigmaF(double* xs, int num_groups);
  void setNuSigmaF(double* xs, int num_groups);
  void setChi(double* xs, int num_groups);

  void setSigmaTByGroup(double xs, int group);
  void setSigmaFByGroup(double xs, int group);
  void setNuSigmaFByGroup(double xs, int group);
  void setSigmaSByGroup(double xs, int origin, int destination);
  void setChiByGroup(double xs, int group);
  void setSigmaAByGroup(double xs, int group);

  void buildFissionMatrix();
  void transposeProductionMatrices();
  void alignData();
  Material* clone();

  std::string toString();
  void printString();
};


#endif /* MATERIAL_H_ */
