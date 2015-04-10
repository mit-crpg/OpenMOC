/**
 * @file ExpEvaluator.h
 * @brief The ExpEvaluator class.
 * @date April 9, 2015.
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef EXPEVALUATOR_H_
#define EXPEVALUATOR_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include "log.h"
#include "PolarQuad.h"
#endif

//FIXME
#ifdef NVCC
#define CUDA_CALLABLE __host__ __device__
#else
#define CUDA_CALLABLE
#endif


/**
 * @class ExpEvaluator ExpEvaluator.h "src/ExpEvaluator.h"
 * @brief This is a class for evaluating exponentials.
 * @details The ExpEvaluator includes different algorithms to evaluate
 *          exponentials with varying degrees of accuracy and speed.
 */
class ExpEvaluator {

private:

  /** A boolean indicating whether or not to use linear interpolation */
  bool _interpolate_exponential;

  /** The exponential linear interpolation table */
  FP_PRECISION* _exp_table;

  /** The spacing for the exponential linear interpolation table */
  FP_PRECISION _exp_table_spacing;

  /** The inverse spacing for the exponential linear interpolation table */
  FP_PRECISION _inverse_exp_table_spacing;

  /** The PolarQuad object of interest */
  PolarQuad* _polar_quad;

  /** Twice the number of polar angles */
  int _two_times_num_polar;
  
public:

  ExpEvaluator();
  virtual ~ExpEvaluator();

  void setPolarQuadrature(PolarQuad* polar_quad);
  void useExponentialInterpolation();
  void useExponentialIntrinsic();

  bool isUsingExponentialInterpolation();

  void initialize(double max_tau, double tolerance);
  FP_PRECISION computeExponential(FP_PRECISION tau, int polar);
};


/**
 * @brief Rounds a single precision floating point value to an integer.
 * @param x a float precision floating point value
 * @brief the rounded integer value
 */
CUDA_CALLABLE inline int round_to_int(float x) {
#ifdef NVCC
  return __float2int_rd(x);
#else
  return lrintf(x);
#endif
}


/**
 * @brief Rounds a double precision floating point value to an integer.
 * @param x a double precision floating point value
 * @brief the rounded integer value
 */
CUDA_CALLABLE inline int round_to_int(double x) {
#ifdef NVCC
  return __double2int_rd(x);
#else
  return lrint(x);
#endif
}


#endif /* EXPEVALUATOR_H_ */
