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


/**
 * @class ExpEvaluator ExpEvaluator.h "src/ExpEvaluator.h"
 * @brief This is a class for evaluating exponentials.
 * @details The ExpEvaluator includes different algorithms to evaluate
 *          exponentials with varying degrees of accuracy and speed. This
 *          is a helper class for the Solver and its subclasses and it not
 *          intended to be initialized as a standalone object.
 */
class ExpEvaluator {

private:

  /** A boolean indicating whether or not to use linear interpolation */
  bool _interpolate;

  /** The inverse spacing for the exponential linear interpolation table */
  FP_PRECISION _inverse_exp_table_spacing;

  /** The number of entries in the exponential linear interpolation table */
  int _table_size;

  /** The exponential linear interpolation table */
  FP_PRECISION* _exp_table;

  /** The PolarQuad object of interest */
  PolarQuad* _polar_quad;

  /** Twice the number of polar angles */
  int _two_times_num_polar;

  /** The maximum optical length a track is allowed to have */
  FP_PRECISION _max_optical_length;

  /** The maximum acceptable approximation error for exponentials */
  FP_PRECISION _exp_precision;
  
public:

  ExpEvaluator();
  virtual ~ExpEvaluator();

  void setPolarQuadrature(PolarQuad* polar_quad);
  void setMaxOpticalLength(FP_PRECISION max_optical_length);
  void setExpPrecision(FP_PRECISION exp_precision);
  void useInterpolation();
  void useIntrinsic();

  FP_PRECISION getMaxOpticalLength();
  FP_PRECISION getExpPrecision();
  bool isUsingInterpolation();
  FP_PRECISION getTableSpacing();
  int getTableSize();
  FP_PRECISION* getExpTable();

  void initialize(double tolerance);
  FP_PRECISION computeExponential(FP_PRECISION tau, int polar);
};


/**
 * @brief Rounds a single precision floating point value to an integer.
 * @param x a float precision floating point value
 * @brief the rounded integer value
 */
inline int round_to_int(float x) {
  return lrintf(x);
}


/**
 * @brief Rounds a double precision floating point value to an integer.
 * @param x a double precision floating point value
 * @brief the rounded integer value
 */
inline int round_to_int(double x) {
  return lrint(x);
}


 #endif /* EXPEVALUATOR_H_ */
