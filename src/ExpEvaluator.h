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
#include "log.h"
#include "Quadrature.h"
#include <math.h>
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
  Quadrature* _quadrature;

  /** The number of polar angles */
  int _num_polar;

  /** The maximum optical length a track is allowed to have */
  FP_PRECISION _max_optical_length;

  /** The maximum acceptable approximation error for exponentials */
  FP_PRECISION _exp_precision;

  bool _solve_3D;

public:

  ExpEvaluator();
  virtual ~ExpEvaluator();

  void setQuadrature(Quadrature* quadrature);
  void setMaxOpticalLength(FP_PRECISION max_optical_length);
  void setExpPrecision(FP_PRECISION exp_precision);
  void setSolve3D(bool solve_3D);
  void useInterpolation();
  void useIntrinsic();

  FP_PRECISION getMaxOpticalLength();
  FP_PRECISION getExpPrecision();
  bool isUsingInterpolation();
  FP_PRECISION getTableSpacing();
  int getTableSize();
  FP_PRECISION* getExpTable();
  bool isSolve3D();

  void initialize();
  FP_PRECISION computeExponential(FP_PRECISION tau, int azim, int polar);
};


/**
 * @brief Computes the exponential term for a optical length and polar angle.
 * @details This method computes \f$ 1 - exp(-\tau/sin(\theta_p)) \f$
 *          for some optical path length and polar angle. This method
 *          uses either a linear interpolation table (default) or the
 *          exponential intrinsic exp(...) function.
 * @param tau the optical path length (e.g., sigma_t times length)
 * @param polar the polar angle index
 * @return the evaluated exponential
 */
inline FP_PRECISION ExpEvaluator::computeExponential(FP_PRECISION tau, int azim,
                                              int polar) {

  FP_PRECISION exponential;

  /* Evaluate the exponential using the lookup table - linear interpolation */
  if (_interpolate) {
    tau = std::min(tau, (_max_optical_length));
    int index = floor(tau * _inverse_exp_table_spacing);
    if (_solve_3D) {
      exponential = (1. - (_exp_table[index * 2] * tau +
                           _exp_table[index * 2 + 1]));
    }
    else {
      index *= _num_polar;
      exponential = (1. - (_exp_table[index + 2 * polar] * tau +
                           _exp_table[index + 2 * polar + 1]));
    }
  }

  /* Evalute the exponential using the intrinsic exp(...) function */
  else {
    if (_solve_3D)
      exponential = 1.0 - exp(-tau);
    else {
      FP_PRECISION sintheta = _quadrature->getSinTheta(azim, polar);
      exponential = 1.0 - exp(-tau / sintheta);
    }
  }

  return exponential;
}

#endif /* EXPEVALUATOR_H_ */
