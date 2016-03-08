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
#include "PolarQuad.h"
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

  /** A boolean indicating whether or not linear source is being used */
  bool _linear_source;

  /** The inverse spacing for the exponential linear interpolation table */
  FP_PRECISION _inverse_exp_table_spacing;
  FP_PRECISION _exp_table_spacing;

  /** The number of entries in the exponential linear interpolation table */
  int _table_size;

  /** The exponential linear interpolation table */
  FP_PRECISION* _exp_table;
  FP_PRECISION* _sin_theta;
  FP_PRECISION* _inv_sin_theta;

  /** The PolarQuad object of interest */
  PolarQuad* _polar_quad;

  /** The number of polar angles */
  int _num_polar;

  /** Three times the number of exponentials stored contiguously in the table */
  int _three_times_num_exp;

  /** The maximum optical length a track is allowed to have */
  FP_PRECISION _max_optical_length;

  /** The maximum acceptable approximation error for exponentials */
  FP_PRECISION _exp_precision;

  /** Get the polar angle index given an index into the exponential table */
  int getPolar(int index);

  /** Get the optical length given an index into the exponential table */
  FP_PRECISION getTau(int index, FP_PRECISION dt);

public:

  ExpEvaluator();
  virtual ~ExpEvaluator();

  void setPolarQuadrature(PolarQuad* polar_quad);
  void setMaxOpticalLength(FP_PRECISION max_optical_length);
  void setExpPrecision(FP_PRECISION exp_precision);
  void useInterpolation();
  void useIntrinsic();
  void useLinearSource();
  void useFlatSource();

  FP_PRECISION getMaxOpticalLength();
  FP_PRECISION getExpPrecision();
  bool isUsingInterpolation();
  FP_PRECISION getTableSpacing();
  FP_PRECISION getInverseTableSpacing();
  int getTableSize();
  FP_PRECISION* getExpTable();

  void initialize();
  FP_PRECISION computeExponentialF1(int index, FP_PRECISION dt,
                                    FP_PRECISION dt2);
  FP_PRECISION computeExponentialF2(int index, FP_PRECISION dt,
                                    FP_PRECISION dt2);
  FP_PRECISION computeExponentialH(int index, FP_PRECISION dt,
                                   FP_PRECISION dt2);
  FP_PRECISION computeExponentialG2(FP_PRECISION tau, int polar);
};


inline int ExpEvaluator::getPolar(int index) {
  return (index % (_three_times_num_exp * _num_polar)) / (3 * _num_polar);
}


inline FP_PRECISION ExpEvaluator::getTau(int index, FP_PRECISION dt) {
  index /= (_three_times_num_exp * _num_polar);
  return index * _exp_table_spacing + dt;
}


/**
 * @brief Computes the F1 exponential term.
 * @details This method computes F1 exponential from Ferrer [1] given
 *          an index into the exponential look-up table, the distance (in units
 *          of optical length) from the corresponding table value and the
 *          requested tau, and that distance squared. This method uses either a
 *          linear interpolation table (default) or the exponential intrinsic
 *          exp(...) function.
 *
 *            [1] R. Ferrer and J. Rhodes III, "A Linear Source Approximation
 *                Scheme for the Method of Characteristics", Nuclear Science and
 *                Engineering, Volume 182, February 2016.
 *
 * @param index the index into the exponential look-up table
 * @param dt the distance to the corresponding look-up table bin
 * @return dt2 the distance to the corresponding look-up table bin squared
 */
inline FP_PRECISION ExpEvaluator::computeExponentialF1(int index,
                                                       FP_PRECISION dt,
                                                       FP_PRECISION dt2) {
  if (_interpolate)
    return _exp_table[index] + _exp_table[index + 1] * dt +
        _exp_table[index + 2] * dt2;
  else {
    int polar = getPolar(index);
    FP_PRECISION tau = getTau(index, dt);
    return 1.0 - exp(- tau * _inv_sin_theta[polar]);
  }
}



/**
 * @brief Computes the F2 exponential term.
 * @details This method computes F2 exponential from Ferrer [1] given
 *          an index into the exponential look-up table, the distance (in units
 *          of optical length) from the corresponding table value and the
 *          requested tau, and that distance squared. This method uses either a
 *          linear interpolation table (default) or the exponential intrinsic
 *          exp(...) function.
 *
 *            [1] R. Ferrer and J. Rhodes III, "A Linear Source Approximation
 *                Scheme for the Method of Characteristics", Nuclear Science and
 *                Engineering, Volume 182, February 2016.
 *
 * @param index the index into the exponential look-up table
 * @param dt the distance to the corresponding look-up table bin
 * @return dt2 the distance to the corresponding look-up table bin squared
 */
inline FP_PRECISION ExpEvaluator::computeExponentialF2(int index,
                                                       FP_PRECISION dt,
                                                       FP_PRECISION dt2) {
  if (_interpolate)
    return _exp_table[index + 3] + _exp_table[index + 4] * dt +
        _exp_table[index + 5] * dt2;
  else {
    int polar = getPolar(index);
    FP_PRECISION tau_m = getTau(index, dt) * _inv_sin_theta[polar];
    FP_PRECISION F1 = 1.0 - exp(- tau_m);
    return 2 * (tau_m - F1) - tau_m * F1;
  }
}


/**
 * @brief Computes the H exponential term.
 * @details This method computes H exponential from Ferrer [1] given
 *          an index into the exponential look-up table, the distance (in units
 *          of optical length) from the corresponding table value and the
 *          requested tau, and that distance squared. This method uses either a
 *          linear interpolation table (default) or the exponential intrinsic
 *          exp(...) function.
 *
 *            [1] R. Ferrer and J. Rhodes III, "A Linear Source Approximation
 *                Scheme for the Method of Characteristics", Nuclear Science and
 *                Engineering, Volume 182, February 2016.
 *
 * @param index the index into the exponential look-up table
 * @param dt the distance to the corresponding look-up table bin
 * @return dt2 the distance to the corresponding look-up table bin squared
 */
inline FP_PRECISION ExpEvaluator::computeExponentialH(int index,
                                                      FP_PRECISION dt,
                                                      FP_PRECISION dt2) {
  if (_interpolate)
    return _exp_table[index + 6] + _exp_table[index + 7] * dt +
        _exp_table[index + 8] * dt2;
  else {
    int polar = getPolar(index);
    FP_PRECISION tau_m = getTau(index, dt) * _inv_sin_theta[polar];
    FP_PRECISION F1 = 1.0 - exp(- tau_m);
    FP_PRECISION G1 = 1 + 0.5 * tau_m - (1 + 1.0 / tau_m) * F1;
    return 0.5 * tau_m - G1;
  }
}

#endif /* EXPEVALUATOR_H_ */
