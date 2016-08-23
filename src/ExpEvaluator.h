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

  /** The spacing for the exponential linear interpolation table */
  FP_PRECISION _exp_table_spacing;

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

  /** The azimuthal angle index used by this exponential evaluator */
  int _azim_index;

  /** The polar angle index used by this exponential evaluator */
  int _polar_index;

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
  FP_PRECISION computeExponentialG2(FP_PRECISION tau);
  EXPEvaluator* copy(); //FIXME
};


//FIXME
inline int ExpEvaluator::getExponentialIndex(FP_PRECISION tau) {
  return floor(tau * _inverse_exp_table_spacing);
}


//FIXME
inline FP_PRECISION ExpEvaluator::getDifference(int index, FP_PRECISION tau) {
  return tau - index * _exp_table_spacing;
}


inline FP_PRECISION ExpEvaluator::convertDistance3Dto2D(FP_PRECISION length) {
  return length * _inv_sin_theta_no_offset;
}


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

  /* Evaluate the exponential using the lookup table - linear interpolation */
  if (_interpolate) {
    tau = std::min(tau, (_max_optical_length));
    int index = floor(tau * _inverse_exp_table_spacing);
    if (_solve_3D) {
      return 1.0 - (_exp_table[index * 2] * tau + _exp_table[index * 2 + 1]);
    }
    else {
      index *= _num_polar;
      return 1.0 - (_exp_table[index + 2 * polar] * tau +
              _exp_table[index + 2 * polar + 1]);
    }
  }

  /* Evalute the exponential using the intrinsic exp(...) function */
  else {
    if (_solve_3D)
      return 1.0 - exp(-tau);
    else {
      FP_PRECISION sintheta = _quadrature->getSinTheta(azim, polar);
      return 1.0 - exp(-tau / sintheta);
    }
  }
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
  int full_index = index * _max_polar_offset + polar_offset;

  if (_interpolate)
    return _exp_table[full_index] + _exp_table[full_index + 1] * dt +
        _exp_table[full_index + 2] * dt2;
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
  int full_index = index * _max_polar_offset + polar_offset;

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
  int full_index = index * _max_polar_offset + polar_offset;

  if (_interpolate)
    return _exp_table[index + 6] + _exp_table[index + 7] * dt +
        _exp_table[index + 8] * dt2;
  else {
    int polar_index = _polar_index + polar_offset;
    FP_PRECISION tau = index * _exp_table_spacing + dt;
    FP_PRECISION inv_sin_theta = 1.0 / _quadrature->getSinTheta(_azim_index,
                                                                polar_index);
    FP_PRECISION tau_m = tau * _inv_sin_theta;
    FP_PRECISION F1 = 1.0 - exp(- tau_m);
    FP_PRECISION G1 = 1 + 0.5 * tau_m - (1 + 1.0 / tau_m) * F1;
    return 0.5 * tau_m - G1;
  }
}


#endif /* EXPEVALUATOR_H_ */
