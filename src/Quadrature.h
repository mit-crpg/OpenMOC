/**
 * @file Quadrature.h
 * @brief The Quadrature abstract class and subclasses.
 * @details The Quadrature subclasses are defined with tabulated or functional
 *          quadrature sets given in the "Lattice Physics Computations", 
 *          Handbook of Nuclear Engineering, Dave Knott, Akio Yamamoto, 2010.
 * @date April 6, 2015
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
*/


#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#ifdef __cplusplus
#include "Python.h"
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
  
  /** An array of the quadrature azimuthal weights */
  FP_PRECISION* _azim_weights;

  /** An array of the quadrature polar weights */
  FP_PRECISION** _polar_weights;

  /** An array of the sines multipled by the weights */
  FP_PRECISION** _multiples;
  
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
  FP_PRECISION getMultiple(int azim, int polar);
  FP_PRECISION** getSinThetas();
  double** getThetas();
  double* getPhis();
  FP_PRECISION* getAzimWeights();
  FP_PRECISION** getPolarWeights();
  FP_PRECISION** getMultiples();
  int getFirstOctantPolar(int polar);
  int getFirstOctantAzim(int azim);
  quadratureType getQuadratureType();
  int getOrthant(int azim, int polar);
  
  virtual void setNumAzimAngles(const int num_azim);
  virtual void setNumPolarAngles(const int num_polar);
  void setThetas(double* thetas, int num_azim_times_polar);
  void setPolarWeights(double* weights, int num_azim_times_polar);
  void setTheta(double theta, int azim, int polar);
  void setPhi(double phi, int azim);
  void setAzimWeight(double weight, int azim);
  void setPolarWeight(double weight, int azim, int polar);
  
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

public:
  GLPolarQuad();
  void setNumPolarAngles(const int num_polar);
  void initialize();
  void precomputeWeights(bool solve_3D);
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


#endif /* QUADRATURE_H_ */
