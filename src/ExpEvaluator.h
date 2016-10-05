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

  /** A boolean indicating whether or not linear source is being used */
  bool _linear_source;

  /** The spacing for the exponential linear interpolation table */
  FP_PRECISION _exp_table_spacing;

  /** The inverse spacing for the exponential linear interpolation table */
  FP_PRECISION _inverse_exp_table_spacing;

  /** The sine of the base polar angle */
  FP_PRECISION _sin_theta_no_offset;

  /** The inverse of the sine of the base polar angle */
  FP_PRECISION _inverse_sin_theta_no_offset;

  /** The number of entries in the exponential linear interpolation table */
  int _table_size;

  /** The exponential linear interpolation table */
  FP_PRECISION* _exp_table;

  /** The PolarQuad object of interest */
  Quadrature* _quadrature;

  /** The maximum optical length a track is allowed to have */
  FP_PRECISION _max_optical_length;

  /** The maximum acceptable approximation error for exponentials */
  FP_PRECISION _exp_precision;

  /** The azimuthal angle index used by this exponential evaluator */
  int _azim_index;

  /** The base polar angle index used by this exponential evaluator */
  int _polar_index;

  /** The number of exponential terms per optical length in the table */
  int _num_exp_terms;

  /** The number of polar angles handled in the exponential table */
  int _num_polar_terms;


public:

  ExpEvaluator();
  virtual ~ExpEvaluator();

  void setQuadrature(Quadrature* quadrature);
  void setMaxOpticalLength(FP_PRECISION max_optical_length);
  void setExpPrecision(FP_PRECISION exp_precision);
  void useInterpolation();
  void useIntrinsic();
  void useLinearSource();

  FP_PRECISION getMaxOpticalLength();
  FP_PRECISION getExpPrecision();
  bool isUsingInterpolation();
  FP_PRECISION getTableSpacing();
  int getTableSize();
  FP_PRECISION* getExpTable();
  int getExponentialIndex(FP_PRECISION tau);
  FP_PRECISION getDifference(int index, FP_PRECISION tau);
  FP_PRECISION convertDistance3Dto2D(FP_PRECISION length);

  void initialize(int azim_index, int polar_index, bool solve_3D);
  FP_PRECISION computeExponential(FP_PRECISION tau, int polar_offset);
  FP_PRECISION computeExponentialF1(int index, int polar_offset,
                                    FP_PRECISION dt, FP_PRECISION dt2);
  FP_PRECISION computeExponentialF2(int index, int polar_offset,
                                    FP_PRECISION dt, FP_PRECISION dt2);
  FP_PRECISION computeExponentialH(int index, int polar_offset,
                                   FP_PRECISION dt, FP_PRECISION dt2);
  void retrieveExponentialComponents(FP_PRECISION tau, int polar_offset,
                                     FP_PRECISION* exp_F1,
                                     FP_PRECISION* exp_F2,
                                     FP_PRECISION* exp_H);

  FP_PRECISION computeExponentialG2(FP_PRECISION tau);
  ExpEvaluator* deepCopy();
};


//FIXME
inline int ExpEvaluator::getExponentialIndex(FP_PRECISION tau) {
  return floor(tau * _inverse_exp_table_spacing);
}


//FIXME
inline FP_PRECISION ExpEvaluator::getDifference(int index, FP_PRECISION tau) {
  return tau - index * _exp_table_spacing;
}


//FIXME
inline FP_PRECISION ExpEvaluator::convertDistance3Dto2D(FP_PRECISION length) {
  return length * _sin_theta_no_offset;
}


//FIXME
inline FP_PRECISION ExpEvaluator::computeExponential(FP_PRECISION tau,
                                                     int polar_offset) {

  /* Extract exponential indexes and differences */
  int exp_index = getExponentialIndex(tau);
  FP_PRECISION dt = getDifference(exp_index, tau);
  FP_PRECISION dt2 = dt * dt;

  /* Compute the exponential */
  return computeExponentialF1(exp_index, polar_offset, dt, dt2);
}


//FIXME
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
                                                       int polar_offset,
                                                       FP_PRECISION dt,
                                                       FP_PRECISION dt2) {

  /* Calculate full index */
  int full_index = (index * _num_polar_terms + polar_offset) * _num_exp_terms;

  if (_interpolate) {
    return _exp_table[full_index] + _exp_table[full_index + 1] * dt +
        _exp_table[full_index + 2] * dt2;
  }
  else {
    int polar_index = _polar_index + polar_offset;
    FP_PRECISION tau = index * _exp_table_spacing + dt;
    FP_PRECISION inv_sin_theta = 1.0 / _quadrature->getSinTheta(_azim_index,
                                                                polar_index);
    return 1.0 - exp(- tau * inv_sin_theta);
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
                                                       int polar_offset,
                                                       FP_PRECISION dt,
                                                       FP_PRECISION dt2) {
  /* Calculate full index */
  int full_index = (index * _num_polar_terms + polar_offset) * _num_exp_terms;

  if (_interpolate)
    return _exp_table[full_index + 3] + _exp_table[full_index + 4] * dt +
        _exp_table[full_index + 5] * dt2;
  else {
    int polar_index = _polar_index + polar_offset;
    FP_PRECISION tau = index * _exp_table_spacing + dt;
    FP_PRECISION inv_sin_theta = 1.0 / _quadrature->getSinTheta(_azim_index,
                                                                polar_index);
    FP_PRECISION tau_m = tau * inv_sin_theta;
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
                                                      int polar_offset,
                                                      FP_PRECISION dt,
                                                      FP_PRECISION dt2) {
  /* Calculate full index */
  int full_index = (index * _num_polar_terms + polar_offset) * _num_exp_terms;

  if (_interpolate)
    return _exp_table[full_index + 6] + _exp_table[full_index + 7] * dt +
        _exp_table[full_index + 8] * dt2;
  else {
    int polar_index = _polar_index + polar_offset;
    FP_PRECISION tau = index * _exp_table_spacing + dt;
    FP_PRECISION inv_sin_theta = 1.0 / _quadrature->getSinTheta(_azim_index,
                                                                polar_index);
    FP_PRECISION tau_m = tau * inv_sin_theta;
    FP_PRECISION F1 = 1.0 - exp(- tau_m);
    FP_PRECISION G1 = 1 + 0.5 * tau_m - (1 + 1.0 / tau_m) * F1;
    return 0.5 * tau_m - G1;
  }
}


//FIXME: description
inline void ExpEvaluator::retrieveExponentialComponents(FP_PRECISION tau,
                                                        int polar_offset,
                                                        FP_PRECISION* exp_F1,
                                                        FP_PRECISION* exp_F2,
                                                        FP_PRECISION* exp_H) {
  if (_interpolate) {

    int exp_index = getExponentialIndex(tau);
    FP_PRECISION dt = getDifference(exp_index, tau);
    FP_PRECISION dt2 = dt * dt;
    int full_index = (exp_index * _num_polar_terms + polar_offset)
      * _num_exp_terms;
    *exp_F1 = _exp_table[full_index] + _exp_table[full_index + 1] * dt +
        _exp_table[full_index + 2] * dt2;
    *exp_F2 = _exp_table[full_index + 3] + _exp_table[full_index + 4] * dt +
        _exp_table[full_index + 5] * dt2;
    *exp_H = _exp_table[full_index + 6] + _exp_table[full_index + 7] * dt +
        _exp_table[full_index + 8] * dt2;
  }
  else {
    int polar_index = _polar_index + polar_offset;
    FP_PRECISION inv_sin_theta = 1.0 / _quadrature->getSinTheta(_azim_index,
                                                                polar_index);
    FP_PRECISION tau_m = tau * inv_sin_theta;
    FP_PRECISION F1 = 1.0 - exp(- tau_m);
    FP_PRECISION G1 = 1 + 0.5 * tau_m - (1 + 1.0 / tau_m) * F1;
    *exp_F1 = F1;
    *exp_F2 = 2 * (tau_m - F1) - tau_m * F1;
    *exp_H = 0.5 * tau_m - G1;
  }
}


#endif /* EXPEVALUATOR_H_ */
