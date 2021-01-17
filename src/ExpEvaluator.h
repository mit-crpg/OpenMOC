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
#ifdef __linux__
#include <malloc.h>
#endif
#include <math.h>
#include "exponentials.h"
#endif


/**
 * @class ExpEvaluator ExpEvaluator.h "src/ExpEvaluator.h"
 * @brief This is a class for evaluating exponentials.
 * @details The ExpEvaluator includes different algorithms to evaluate
 *          exponentials with varying degrees of accuracy and speed. This is a
 *          helper class for the Solver and its subclasses and it is not
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
  FP_PRECISION getInverseSinTheta();

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


/**
 * @brief Get the index on the exponential interpolation grid of the value right
 *        beneath tau.
 * @param tau optical distance
 * @return the index on the exponential interpolation grid
 */
inline int ExpEvaluator::getExponentialIndex(FP_PRECISION tau) {
  return int(tau * _inverse_exp_table_spacing);
}


/**
 * @brief Compute the difference between an optical path and an indexed value in
 *        the exponential interpolation grid.
 * @param index index on the exponential interpolation grid
 * @param tau optical distance
 * @return the difference between tau and the value on the grid
 */
inline FP_PRECISION ExpEvaluator::getDifference(int index, FP_PRECISION tau) {
  return tau - index * _exp_table_spacing;
}


/**
 * @brief Convert a 3D distance to a 2D based on the evaluator's polar angle.
 * @param length the 3D distance
 * @return the 2D distance
 */
inline FP_PRECISION ExpEvaluator::convertDistance3Dto2D(FP_PRECISION length) {
  return length * _sin_theta_no_offset;
}


/**
 * @brief Get the inverse of sin theta from the ExpEvaluator
 * @return inverse sin theta for the first angle
 */
inline FP_PRECISION ExpEvaluator::getInverseSinTheta() {
  return _inverse_sin_theta_no_offset;
}


/**
 * @brief Computes the F1 exponential term.
 * @param tau the optical distance (2D)
 * @param polar_offset an offset to the index in the look-up table
 * @return the F1 exponential term
 */
inline FP_PRECISION ExpEvaluator::computeExponential(FP_PRECISION tau,
                                                     int polar_offset) {

#ifndef THREED
  FP_PRECISION inv_sin_theta = _quadrature->getInvSinThetaInline(_azim_index,
                                                  _polar_index + polar_offset);
#else
  FP_PRECISION inv_sin_theta = _inverse_sin_theta_no_offset;
#endif
  FP_PRECISION exp_F1;
  expF1_fractional(tau * inv_sin_theta, &exp_F1);

  return inv_sin_theta * exp_F1;
}


/**
 * @brief Computes the F1 exponential term.
 * @details This method computes F1 exponential from Ferrer [1] given
 *          an index into the exponential look-up table, the distance (in units
 *          of optical length) from the corresponding table value and the
 *          requested tau, and that distance squared. This method uses either a
 *          linear interpolation table (default) or the exponential intrinsic
 *          exp(...) function. //DEPRECATED
 *
 *            [1] R. Ferrer and J. Rhodes III, "A Linear Source Approximation
 *                Scheme for the Method of Characteristics", Nuclear Science and
 *                Engineering, Volume 182, February 2016.
 *
 * @param index the index into the exponential look-up table
 * @param polar_offset an offset to the index in the look-up table
 * @param dt the distance to the corresponding look-up table bin
 * @param dt2 the distance to the corresponding look-up table bin squared
 * @return the evaluated F1 exponential term
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
    return (1.0 - exp(- tau * inv_sin_theta)) / tau;
  }
}


/**
 * @brief Computes the F2 exponential term.
 * @details This method computes F2 exponential from Ferrer [1] given
 *          an index into the exponential look-up table, the distance (in units
 *          of optical length) from the corresponding table value and the
 *          requested tau, and that distance squared. This method uses either a
 *          linear interpolation table (default) or the exponential intrinsic
 *          exp(...) function. //DEPRECATED
 *
 *            [1] R. Ferrer and J. Rhodes III, "A Linear Source Approximation
 *                Scheme for the Method of Characteristics", Nuclear Science and
 *                Engineering, Volume 182, February 2016.
 *
 * @param index the index into the exponential look-up table
 * @param polar_offset an offset to the index in the look-up table
 * @param dt the distance to the corresponding look-up table bin
 * @param dt2 the distance to the corresponding look-up table bin squared
 * @return the evaluated F2 exponential term
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
    FP_PRECISION F1 = (1.0 - exp(- tau_m)) / tau;
    return 2.0 / tau * (inv_sin_theta - F1) - inv_sin_theta * F1;
  }
}


/**
 * @brief Computes the H exponential term.
 * @details This method computes H exponential from Ferrer [1] given
 *          an index into the exponential look-up table, the distance (in units
 *          of optical length) from the corresponding table value and the
 *          requested tau, and that distance squared. This method uses either a
 *          linear interpolation table (default) or the exponential intrinsic
 *          exp(...) function. //DEPRECATED
 *
 *            [1] R. Ferrer and J. Rhodes III, "A Linear Source Approximation
 *                Scheme for the Method of Characteristics", Nuclear Science and
 *                Engineering, Volume 182, February 2016.
 *
 * @param index the index into the exponential look-up table
 * @param polar_offset an offset to the index in the look-up table
 * @param dt the distance to the corresponding look-up table bin
 * @param dt2 the distance to the corresponding look-up table bin squared
 * @return the evaluated H exponential term
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
    FP_PRECISION F1 = (1.0 - exp(- tau_m)) / tau;
    FP_PRECISION G1 = 1.0 / tau + 0.5 * inv_sin_theta - (1.0 + 1.0 / tau_m) * F1;
    return 0.5 * inv_sin_theta - G1;
  }
}


/**
 * @brief Computes the G2 exponential term for a optical length and polar angle.
 * @details This method computes the G2 exponential term from Ferrer [1]
 *          for some optical path length and polar angle. This method uses a
 *          rational fraction approximation to compute the exponential term.
 *
 *            [1] R. Ferrer and J. Rhodes III, "A Linear Source Approximation
 *                Scheme for the Method of Characteristics", Nuclear Science and
 *                Engineering, Volume 182, February 2016.
 *
 * @param tau the optical path length (e.g., sigma_t times length)
 * @return the evaluated exponential
 */
inline FP_PRECISION ExpEvaluator::computeExponentialG2(FP_PRECISION tau) {

  FP_PRECISION exp_G2;
  expG2_fractional(tau, &exp_G2);

  return exp_G2;
}


/**
 * @brief Computes the F1, F2, H exponential term.
 * @details This method computes F1, F2, H exponential from Ferrer [1] given the
 *          requested tau. This method uses either a linear interpolation table
 *          (default) or the exponential intrinsicexp(...) function.
 *
 *            [1] R. Ferrer and J. Rhodes III, "A Linear Source Approximation
 *                Scheme for the Method of Characteristics", Nuclear Science and
 *                Engineering, Volume 182, February 2016.
 *
 * @param tau the requested tau
 * @param polar_offset an offset to the index in the look-up table
 * @param exp_F1 pointer to the F1 exponential term to be computed
 * @param exp_F2 pointer to the F2 exponential term to be computed
 * @param exp_H  pointer to the H exponential term to be computed
 */
inline void ExpEvaluator::retrieveExponentialComponents(FP_PRECISION tau,
                                                        int polar_offset,
#ifdef SWIG  //FIXME Find out how to use restrict with SWIG
                                                        FP_PRECISION* exp_F1,
                                                        FP_PRECISION* exp_F2,
                                                        FP_PRECISION* exp_H) {
#else
                                              FP_PRECISION* __restrict__ exp_F1,
                                              FP_PRECISION* __restrict__ exp_F2,
                                              FP_PRECISION* __restrict__ exp_H) {
#endif

#ifndef THREED
  FP_PRECISION inv_sin_theta = _quadrature->getInvSinThetaInline(_azim_index,
                                                  _polar_index + polar_offset);
#else
  FP_PRECISION inv_sin_theta = _inverse_sin_theta_no_offset;
#endif

  /* Limit range of tau to avoid numerical errors */
  tau = std::max(FP_PRECISION(1e-8), tau * inv_sin_theta);

  /* Compute exponentials from a common exponential */
   FP_PRECISION exp_G;
   expG_fractional(tau, &exp_G);
   *exp_F1 = 1.f - tau*exp_G;
   *exp_F1 *= inv_sin_theta;

   exp_G *= inv_sin_theta;
   *exp_F2 = 2.f*exp_G - *exp_F1;

   *exp_H = *exp_F1 - exp_G;

    /* Quadratic exponential interpolation tables */
//     __builtin_assume_aligned(_exp_table, VEC_ALIGNMENT);
//
//     tau /= inv_sin_theta;
//     int exp_index = getExponentialIndex(tau);
//     FP_PRECISION dt = getDifference(exp_index, tau);
//     FP_PRECISION dt2 = dt * dt;
//     int full_index = (exp_index * _num_polar_terms + polar_offset)
//       * _num_exp_terms;
//     *exp_F1 = _exp_table[full_index] + _exp_table[full_index + 1] * dt +
//         _exp_table[full_index + 2] * dt2;
//     *exp_F2 = _exp_table[full_index + 3] + _exp_table[full_index + 4] * dt +
//         _exp_table[full_index + 5] * dt2;
//     *exp_H = _exp_table[full_index + 6] + _exp_table[full_index + 7] * dt +
//         _exp_table[full_index + 8] * dt2;
}


#endif /* EXPEVALUATOR_H_ */
