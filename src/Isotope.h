/**
 * @file Isotope.h
 * @brief Container for isotopic microscopic cross sections
 * @date February 13, 2014
 * @author Nathan A. Gibson, MIT, Course 22 (ngibson@mit.edu)
 */

#ifndef ISOTOPE_H_
#define ISOTOPE_H_

#ifdef __cplusplus
#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <vector>
#include "log.h"
#endif


/**
 * @class Isotope Isotope.h "src/Isotope.h"
 * @brief The Isotope class represents an istope and its relevant
 *        nuclear data (i.e., multigroup cross-sections) for neutron transport.
 */
class Isotope {

protected:

  /** A static counter for the number of materials in a simulation */
  static int _n;

  /** A monotonically increasing unique ID for each Material created */
  int _uid;

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

  /** A matrix of group to group scattering cross-sections */
  FP_PRECISION* _sigma_s;

  /** An array of the fission cross-sections for each energy group */
  FP_PRECISION* _sigma_f;

  /** An array of the fission cross-sections multiplied by nu \f$ \nu \f$
   *  for each energy group */
  FP_PRECISION* _nu_sigma_f;

  /** An array of the chi \f$ \chi \f$ values for each energy group */
  FP_PRECISION* _chi;

  /** A boolean representing whether or not this Material contains a non-zero
   *  fission cross-section and is fissionable */
  bool _fissionable;

  
public:
   Isotope(int id);
  ~Isotope();
  
  int getUid() const;
  int getId() const;
  int getNumEnergyGroups() const;
  
  FP_PRECISION* getSigmaT();
  FP_PRECISION* getSigmaA(); 
  FP_PRECISION* getSigmaS();
  FP_PRECISION* getSigmaF();
  FP_PRECISION* getNuSigmaF();
  FP_PRECISION* getChi();
  
  FP_PRECISION getSigmaTByGroup(int group);
  FP_PRECISION getSigmaAByGroup(int group); 
  FP_PRECISION getSigmaSByGroup(int origin, int destination);
  FP_PRECISION getScatterSource(int group, FP_PRECISION* flux);
  FP_PRECISION getSigmaFByGroup(int group);
  FP_PRECISION getNuSigmaFByGroup(int group);
  FP_PRECISION getChiByGroup(int group);
  
  void setNumEnergyGroups(const int num_groups);

  void setSigmaT(double* xs, int num_groups);
  void setSigmaA(double* xs, int num_groups);
  void setSigmaF(double* xs, int num_groups);
  void setSigmaS(double* xs, int num_groups_squared);
  void setNuSigmaF(double* xs, int num_groups);
  void setChi(double* xs, int num_groups);
  
  void setSigmaTByGroup(double xs, int group);
  void setSigmaAByGroup(double xs, int group);
  void setSigmaSByGroup(double xs, int origin, int destination);
  void setSigmaFByGroup(double xs, int group);
  void setNuSigmaFByGroup(double xs, int group);
  void setChiByGroup(double xs, int group);

  
};



#endif

