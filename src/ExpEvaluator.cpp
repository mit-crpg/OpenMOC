#include "ExpEvaluator.h"


/**
 * @brief Constructor initializes array pointers to NULL.
 * @details The constructor sets the interpolation scheme as the default
 *          for computing exponentials.
 */
ExpEvaluator::ExpEvaluator() {
  _interpolate = true;
  _exp_table = NULL;
  _quadrature = NULL;
  _max_optical_length = MAX_OPTICAL_LENGTH;
  _exp_precision = EXP_PRECISION;
}


/**
 * @brief Destructor deletes table for linear interpolation of exponentials
 */
ExpEvaluator::~ExpEvaluator() {
  if (_exp_table != NULL)
    delete [] _exp_table;
}


/**
 * @brief Set the Quadrature to use when computing exponentials.
 * @param quadrature a Quadrature object pointer
 */
void ExpEvaluator::setQuadrature(Quadrature* quadrature) {
  _quadrature = quadrature;
  _num_polar = _quadrature->getNumPolarAngles();
}


/**
 * @brief Sets the maximum optical length covered in the exponential
 *        interpolation table.
 * @param max_optical_length the maximum optical length
 */
void ExpEvaluator::setMaxOpticalLength(FP_PRECISION max_optical_length) {

  if (max_optical_length <= 0)
    log_printf(ERROR, "Cannot set max optical length to %f because it "
               "must be positive.", max_optical_length);

  _max_optical_length = max_optical_length;
}


/**
 * @brief Sets the maximum acceptable approximation error for exponentials.
 * @details This routine only affects the construction of the linear
 *          interpolation table for exponentials, if in use. By default,
 *          a value of 1E-5 is used for the table, as recommended by the
 *          analysis of Yamamoto in his 2004 paper on the subject.
 * @param exp_precision the maximum exponential approximation error
 */
void ExpEvaluator::setExpPrecision(FP_PRECISION exp_precision) {

  if (exp_precision <= 0)
    log_printf(ERROR, "Cannot set exp precision to %f because it "
               "must be positive.", exp_precision);

  _exp_precision = exp_precision;
}


/**
 * @brief Use linear interpolation to compute exponentials.
 */
void ExpEvaluator::useInterpolation() {
  _interpolate = true;
}


/**
 * @brief Use the exponential intrinsic exp(...) to compute exponentials.
 */
void ExpEvaluator::useIntrinsic() {
  _interpolate = false;
}


/**
 * @brief Gets the maximum optical length covered with the exponential
 *        interpolation table.
 * @return max_optical_length the maximum optical length
 */
FP_PRECISION ExpEvaluator::getMaxOpticalLength() {
  return _max_optical_length;
}


/**
 * @brief Gets the maximum acceptable approximation error for exponentials.
 * @return the maximum exponential approximation error
 */
FP_PRECISION ExpEvaluator::getExpPrecision() {
  return _exp_precision;
}


/**
 * @brief Returns true if using linear interpolation to compute exponentials.
 * @return true if so, false otherwise
 */
bool ExpEvaluator::isUsingInterpolation() {
  return _interpolate;
}


/**
 * @brief Returns the exponential table spacing.
 * @return exponential table spacing
 */
FP_PRECISION ExpEvaluator::getTableSpacing() {

  if (_exp_table == NULL)
    log_printf(ERROR, "Unable to return the exponential table spacing "
               "since it has not yet been initialized");

  return 1.0 / _inverse_exp_table_spacing;
}


/**
 * @brief Get the number of entries in the exponential interpolation table.
 * @param entries in the interpolation table
 *
 */
int ExpEvaluator::getTableSize() {

  if (_exp_table == NULL)
    log_printf(ERROR, "Unable to return exponential table size "
               "since it has not yet been initialized");

  return _table_size;
}


/**
 * @brief Returns a pointer to the exponential interpolation table.
 * @return pointer to the exponential interpolation table
 */
FP_PRECISION* ExpEvaluator::getExpTable() {

  if (_exp_table == NULL)
    log_printf(ERROR, "Unable to return exponential table "
               "since it has not yet been initialized");

  return _exp_table;
}


/**
 * @brief If using linear interpolation, builds the table for each polar angle.
 * @param tolerance the minimum acceptable interpolation accuracy
 */
void ExpEvaluator::initialize() {

  /* If no exponential table is needed, return */
  if (!_interpolate)
    return;

  log_printf(INFO, "Initializing exponential interpolation table...");

  /* Set size of interpolation table */
  int num_polar = _quadrature->getNumPolarAngles();
  int num_array_values = _max_optical_length * sqrt(1. / (8. * _exp_precision));
  FP_PRECISION exp_table_spacing = _max_optical_length / num_array_values;

  /* Increment the number of vaues in the array to ensure that a tau equal to
   * max_optical_length resides as the final entry in the table */
  num_array_values += 1;

  /* Compute the reciprocal of the table entry spacing */
  _inverse_exp_table_spacing = 1.0 / exp_table_spacing;

  /* Allocate array for the table */
  if (_exp_table != NULL)
    delete [] _exp_table;

  _table_size = 2 * num_polar * num_array_values;
  _exp_table = new FP_PRECISION[_table_size];

  FP_PRECISION expon;
  FP_PRECISION intercept;
  FP_PRECISION slope;
  FP_PRECISION sin_theta;
  FP_PRECISION tau;

  /* Create exponential linear interpolation table */
  for (int i=0; i < num_array_values; i++) {
    for (int p=0; p < num_polar; p++) {
      sin_theta = _quadrature->getSinTheta(0, p);

      /* Use the optical length at the start of the interval for the first value
       * to avoid exponential values greater than one. */
      if (i == 0)
        tau = i * exp_table_spacing;

      /* Use the optical length at the interval mid-point to reduce error. */
      else
        tau = (i + 0.5) * exp_table_spacing;

      expon = exp(- tau / sin_theta);
      intercept = expon * (1 + tau / sin_theta);
      slope = - expon / sin_theta;
      _exp_table[_num_polar * i + 2 * p] = slope;
      _exp_table[_num_polar * i + 2 * p + 1] = intercept;
    }
  }
}
