/**
 * @file PolarQuad.h
 * @brief The PolarQuad abstract class and subclasses.
 * @date April 6, 2015
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
*/


#ifndef POLARQUAD_H_
#define POLARQUAD_H_

#ifdef __cplusplus
#include <sstream>
#include "log.h"
#endif


/**
 * @class PolarQuad PolarQuad.h "src/PolarQuad.h"
 * @brief The arbitrary polar quadrature parent class.
 */
class PolarQuad {

private:

  /** The number of polar angles */
  int _num_polar;

  /** An array of the sines of quadrature angles */
  FP_PRECISION* _sin_thetas;

  /** An array of the quadrature weights */
  FP_PRECISION* _weights;

  /** An array of the sines multipled by the weights */
  FP_PRECISION* _multiples;

public:

  PolarQuad();
  virtual ~PolarQuad();

  int getNumPolarAngles() const;
  FP_PRECISION getSinTheta(const int n) const;
  FP_PRECISION getWeight(const int n) const;
  FP_PRECISION getMultiple(const int n) const;
  FP_PRECISION* getSinThetas();
  FP_PRECISION* getWeights();
  FP_PRECISION* getMultiples();

  void setNumPolarAngles(const int num_polar);
  void setSinThetas(FP_PRECISION* sin_thetas, int num_polar);
  void setWeights(FP_PRECISION* weights, int num_polar);

  void initialize();

  std::string toString();
};



/**
 * @class TYPolarQuad PolarQuad.h "src/PolarQuad.h"
 * @brief Tabuchi-Yamamoto's polar quadrature.
 */
class TYPolarQuad : public PolarQuad {

private:

public:

  TYPolarQuad();
  virtual ~TYPolarQuad();

  void setNumPolarAngles();
  void initialize();
    
}


/**
 * @class LeonardPolarQuad PolarQuad.h "src/PolarQuad.h"
 * @brief Leonard's polar quadrature.
 */
class LeonardPolarQuad : public PolarQuad { 

private:

public:

  LeonardPolarQuad();
  virtual ~LeonardPolarQuad();
  void setNumPolarAngles();
  void initialize();
}



/**
 * @class GLPolarQuad PolarQuad.h "src/PolarQuad.h"
 * @brief Gauss-Legendre's polar quadrature.
 */
class GLPolarQuad : public PolarQuad { 

private:

public:
  GLPolarQuad();
  virtual ~GLPolarQuad();
  void setNumPolarAngles();
  void initialize();
}



/**
 * @class EqualAnglePolarQuad PolarQuad.h "src/PolarQuad.h"
 * @brief Equal angle polar quadrature.
 */
class EqualAnglePolarQuad : public PolarQuad { 

private:

public:

  EqualAnglePolarQuad(int num_polar);
  virtual ~EqualAnglePolarQuad();
  void initialize();
}



/**
 * @class EqualWeightPolarQuad PolarQuad.h "src/PolarQuad.h"
 * @brief Equal weight polar quadrature.
 */
class EqualWeightPolarQuad : public PolarQuad { 

private:

public:

  EqualWeightPolarQuad(int num_polar);
  virtual ~EqualWeightPolarQuad();
  void initialize();
}


#endif /* POLARQUAD_H_ */
