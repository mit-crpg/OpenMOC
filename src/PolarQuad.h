/**
 * @file PolarQuad.h
 * @brief The PolarQuad abstract class and subclasses.
 * @details The PolarQuad subclasses are defined with tabulated or functional
 *          quadrature sets given in the "Lattice Physics Computations",
 *          Handbook of Nuclear Engineering, Dave Knott, Akio Yamamoto, 2010.
 * @date April 6, 2015
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
*/


#ifndef POLARQUAD_H_
#define POLARQUAD_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "constants.h"
#include "log.h"
#include <sstream>
#endif


/**
 * @enum quadratureType
 * @brief The types of quadrature sets supported by OpenMOC.
 */
enum quadratureType {
  TABUCHI_YAMAMOTO,
  LEONARD,
  GAUSS_LEGENDRE,
  EQUAL_WEIGHTS,
  EQUAL_ANGLES,
  CUSTOM
};


/**
 * @class PolarQuad PolarQuad.h "src/PolarQuad.h"
 * @brief The arbitrary polar quadrature parent class.
 */
class PolarQuad {

protected:

  /** The quadrature type being used */
  quadratureType _quad_type;

  /** The number of polar angles */
  int _num_polar;

  /** An array of the sines of quadrature angles */
  FP_PRECISION* _sin_thetas;
  FP_PRECISION* _inverse_sin_thetas;

  /** An array of the quadrature weights */
  FP_PRECISION* _weights;

  /** An array of the sines multipled by the weights */
  FP_PRECISION* _multiples;

  void precomputeMultiples();

public:

  PolarQuad();
  virtual ~PolarQuad();

  int getNumPolarAngles() const;
  FP_PRECISION getSinTheta(const int n) const;
  FP_PRECISION getInverseSinTheta(const int n) const;
  FP_PRECISION getWeight(const int n) const;
  FP_PRECISION getMultiple(const int n) const;
  FP_PRECISION* getSinThetas();
  FP_PRECISION* getInverseSinThetas();
  FP_PRECISION* getWeights();
  FP_PRECISION* getMultiples();
  quadratureType getQuadratureType();

  virtual void setNumPolarAngles(const int num_polar);
  void setSinThetas(double* sin_thetas, int num_polar);
  void setWeights(double* weights, int num_polar);

  virtual void initialize();

  std::string toString();
};



/**
 * @class TYPolarQuad PolarQuad.h "src/PolarQuad.h"
 * @brief Tabuchi-Yamamoto's polar quadrature.
 */
class TYPolarQuad: public PolarQuad {

private:

public:
  TYPolarQuad();
  void setNumPolarAngles(const int num_polar);
  void initialize();
};


/**
 * @class LeonardPolarQuad PolarQuad.h "src/PolarQuad.h"
 * @brief Leonard's polar quadrature.
 */
class LeonardPolarQuad: public PolarQuad {

private:

public:
  LeonardPolarQuad();
  void setNumPolarAngles(const int num_polar);
  void initialize();
};



/**
 * @class GLPolarQuad PolarQuad.h "src/PolarQuad.h"
 * @brief Gauss-Legendre's polar quadrature.
 */
class GLPolarQuad: public PolarQuad {

private:

public:
  GLPolarQuad();
  void setNumPolarAngles(const int num_polar);
  void initialize();
};



/**
 * @class EqualWeightsPolarQuad PolarQuad.h "src/PolarQuad.h"
 * @brief Equal weights polar quadrature.
 */
class EqualWeightsPolarQuad: public PolarQuad {

private:

public:
  EqualWeightsPolarQuad();
  void setNumPolarAngles(const int num_polar);
  void initialize();
};



/**
 * @class EqualAnglesPolarQuad PolarQuad.h "src/PolarQuad.h"
 * @brief Equal angles polar quadrature.
 */
class EqualAnglesPolarQuad: public PolarQuad {

private:

public:
  EqualAnglesPolarQuad();
  void setNumPolarAngles(const int num_polar);
  void initialize();
};


#endif /* POLARQUAD_H_ */
