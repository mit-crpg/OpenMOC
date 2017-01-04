/**
 * @file Quadrature.h
 * @brief The Quadrature abstract class and subclasses.
 * @details The Quadrature subclasses are defined with tabulated or functional
 *          quadrature sets given in the "Lattice Physics Computations",
 *          Handbook of Nuclear Engineering, Dave Knott, Akio Yamamoto, 2010.
 * @date April 8, 2016
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
*/


#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "constants.h"
#include "log.h"
#include <sstream>
#endif
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>


/**
 * @enum quadratureType
 * @brief The types of quadrature sets supported by OpenMOC.
 */
enum quadratureType {
  TABUCHI_YAMAMOTO,
  LEONARD,
  GAUSS_LEGENDRE,
  EQUAL_WEIGHT,
  EQUAL_ANGLE
};


/**
 * @class Quadrature Quadrature.h "src/Quadrature.h"
 * @brief The arbitrary quadrature parent class.
 */
class Quadrature {

protected:

  /** The quadrature type being used */
  quadratureType _quad_type;

  /** The number of azimuthal angles in (0, 2*PI) */
  int _num_azim;

  /** The number of polar angles in (0, PI) */
  int _num_polar;

  /** An array of the sines of quadrature polar angles */
  FP_PRECISION** _sin_thetas;

  /** An array of the quadrature polar angles */
  double** _thetas;

  /** An array of the quadrature azimuthal angles */
  double* _phis;

  /** The actual track azimuthal spacing (cm) by azimuthal angle */
  FP_PRECISION* _azim_spacings;

  /** An array of the quadrature azimuthal weights */
  FP_PRECISION* _azim_weights;

  /** An array of the quadrature polar weights */
  FP_PRECISION** _polar_weights;

  /** An array of the total weights for each azimuthal/polar angle pair */
  FP_PRECISION** _total_weights;


  /* Templates for setting the same values to complimentary and supplimentary
   * angles */
  template <typename T>
  void setPolarValues(T** array, int azim_index, int polar_index, T value) {
    array[azim_index][polar_index] = value;
    array[_num_azim/2 - azim_index - 1][polar_index] = value;
    array[azim_index][_num_polar - polar_index - 1] = value;
    array[_num_azim/2 - azim_index - 1][_num_polar - polar_index - 1] = value;
  }

  template <typename T>
  void setAzimuthalValues(T* array, int azim_index, T value) {
    array[azim_index] = value;
    array[_num_azim/2 - azim_index - 1] = value;
  }

  /* Functions to delete and reset allocated arrays */
  void deletePolarArrays();
  void deleteAllArrays();


public:

  Quadrature();
  virtual ~Quadrature();

  int getNumPolarAngles() const;
  int getNumAzimAngles() const;
  FP_PRECISION getSinTheta(int azim, int polar);
  double getTheta(int azim, int polar);
  double getPhi(int azim);
  FP_PRECISION getAzimWeight(int azim);
  FP_PRECISION getPolarWeight(int azim, int polar);
  FP_PRECISION getWeight(int azim, int polar);
  FP_PRECISION getWeightInline(int azim, int polar);
  FP_PRECISION** getSinThetas();
  double** getThetas();
  double* getPhis();
  FP_PRECISION* getAzimWeights();
  FP_PRECISION** getPolarWeights();
  quadratureType getQuadratureType();
  FP_PRECISION* getAzimSpacings();
  FP_PRECISION getAzimSpacing(int azim);

  virtual void setNumAzimAngles(const int num_azim);
  virtual void setNumPolarAngles(const int num_polar);
  void setThetas(double* thetas, int num_azim_times_polar);
  void setPolarWeights(FP_PRECISION* weights, int num_azim_times_polar);
  void setTheta(double theta, int azim, int polar);
  void setPhi(double phi, int azim);
  void setAzimSpacing(FP_PRECISION spacing, int azim);
  void setAzimWeight(double weight, int azim);
  void setPolarWeight(FP_PRECISION weight, int azim, int polar);

  virtual void initialize();
  virtual void precomputeWeights(bool solve_3D);

  std::string toString();
};



/**
 * @class TYPolarQuad Quadrature.h "src/Quadrature.h"
 * @brief Tabuchi-Yamamoto's polar quadrature.
 */
class TYPolarQuad: public Quadrature {

private:

public:
  TYPolarQuad();
  void setNumPolarAngles(const int num_polar);
  void initialize();
  void precomputeWeights(bool solve_3D);
};


/**
 * @class LeonardPolarQuad Quadrature.h "src/Quadrature.h"
 * @brief Leonard's polar quadrature.
 */
class LeonardPolarQuad: public Quadrature {

private:

public:
  LeonardPolarQuad();
  void setNumPolarAngles(const int num_polar);
  void initialize();
  void precomputeWeights(bool solve_3D);
};



/**
 * @class GLPolarQuad Quadrature.h "src/Quadrature.h"
 * @brief Gauss-Legendre's polar quadrature.
 */
class GLPolarQuad: public Quadrature {

private:

  /** the roots to the Legendre polynomial of degree _num_polar */
  std::vector <double> _roots;

public:
  GLPolarQuad();
  void setNumPolarAngles(const int num_polar);
  void initialize();
  void precomputeWeights(bool solve_3D);

  double legendrePolynomial(int n, double x);
  double logDerivLegendre(int n, double x);
  double secondLogDerivLegendre(int n, double x);
  std::vector <double> getLegendreRoots(int n);
  std::vector <double> getGLWeights(std::vector <double> roots, int n);
};



/**
 * @class EqualWeightPolarQuad Quadrature.h "src/Quadrature.h"
 * @brief Equal weight polar quadrature.
 */
class EqualWeightPolarQuad: public Quadrature {

private:

public:
  EqualWeightPolarQuad();
  void setNumPolarAngles(const int num_polar);
  void initialize();
  void precomputeWeights(bool solve_3D);
};



/**
 * @class EqualAnglePolarQuad Quadrature.h "src/Quadrature.h"
 * @brief Equal angle polar quadrature.
 */
class EqualAnglePolarQuad: public Quadrature {

private:

public:
  EqualAnglePolarQuad();
  void setNumPolarAngles(const int num_polar);
  void initialize();
  void precomputeWeights(bool solve_3D);
};


/**
 * @brief Returns the total weight for Tracks with the given azimuthal and
 *        polar indexes without error checking and inlined
 * @details Angular weights are multiplied by Track spcings
 * @param azim index of the azimuthal angle of interest
 * @param polar index of the polar angle of interest
 * @return the total weight of each Track with the given indexes
 */
inline FP_PRECISION Quadrature::getWeightInline(int azim, int polar) {
  return _total_weights[azim][polar];
}

#endif /* QUADRATURE_H_ */
