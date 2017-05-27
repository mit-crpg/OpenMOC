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
  _sin_theta_no_offset = 0.0;
  _inverse_sin_theta_no_offset = 0.0;
  _linear_source = false;
  _num_exp_terms = 3;
  _azim_index = 0;
  _polar_index = 0;
  _num_polar_terms = 0;
  _inverse_exp_table_spacing = 1.0;
  _exp_table_spacing = 1.0;
}


/**
 * @brief Destructor deletes table for linear interpolation of exponentials
 */
ExpEvaluator::~ExpEvaluator() {
  if (_exp_table != NULL)
    delete [] _exp_table;
}


/**
 * @brief Set the PolarQuad to use when computing exponentials.
 * @param polar_quad a PolarQuad object pointer
 */
void ExpEvaluator::setQuadrature(Quadrature* quadrature) {
  _quadrature = quadrature;
}


/**
 * @brief Sets the maximum optical length covered in the exponential
 *        interpolation table.
 * @param max_optical_length the maximum optical length
 */
void ExpEvaluator::setMaxOpticalLength(NEW_PRECISION max_optical_length) {

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
void ExpEvaluator::setExpPrecision(NEW_PRECISION exp_precision) {

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
 * @brief Use linear source exponentials.
 */
void ExpEvaluator::useLinearSource() {
  _linear_source = true;
  _num_exp_terms = 9;
}


/**
 * @brief Gets the maximum optical length covered with the exponential
 *        interpolation table.
 * @return max_optical_length the maximum optical length
 */
NEW_PRECISION ExpEvaluator::getMaxOpticalLength() {
  return _max_optical_length;
}


/**
 * @brief Gets the maximum acceptable approximation error for exponentials.
 * @return the maximum exponential approximation error
 */
NEW_PRECISION ExpEvaluator::getExpPrecision() {
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
NEW_PRECISION ExpEvaluator::getTableSpacing() {

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
NEW_PRECISION* ExpEvaluator::getExpTable() {

  if (_exp_table == NULL)
    log_printf(ERROR, "Unable to return exponential table "
               "since it has not yet been initialized");

  return _exp_table;
}


/**
 * @brief If using linear interpolation, builds the table for each polar angle.
 */
void ExpEvaluator::initialize(int azim_index, int polar_index, bool solve_3D) {

  /* Check for a valid Quadrature */
  if (_quadrature == NULL)
    log_printf(ERROR, "A Quadrature must be set before an Exponential "
               "Evaluator can be initialized");

  /* Extract the number of azimuthal and polar angles */
  int num_azim = _quadrature->getNumAzimAngles();
  int num_polar = _quadrature->getNumPolarAngles();

  /* Check for a valid azimuthal angle index */
  if (azim_index < 0 || azim_index > num_azim / 4)
    log_printf(ERROR, "Invalid azimuthal angle index of %d. The index must be "
               "betwen 0 and %d", azim_index, num_azim/4);

  /* Check for a valid polar angle index */
  if (polar_index < 0 || polar_index > num_polar / 2)
    log_printf(ERROR, "Invalid polar angle index of %d. The index must be "
               "betwen 0 and %d", polar_index, num_polar/2);

  /* Record the azimuthal angle index */
  _azim_index = azim_index;

  /* Determine the base polar angle index and maximum offset */
  if (solve_3D) {
    _polar_index = polar_index;
    _num_polar_terms = 1;
  }
  else {
    _polar_index = 0;
    _num_polar_terms = num_polar / 2;
  }

  /* Record the inverse sine for the base polar angle */
  _sin_theta_no_offset = _quadrature->getSinTheta(azim_index,
                                                  polar_index);
  _inverse_sin_theta_no_offset = 1.0 / _sin_theta_no_offset;

  log_printf(INFO, "Initializing exponential interpolation table...");

  /* If no exponential table is needed, return */
  if (!_interpolate)
    return;

  /* Set size of interpolation table */
  int num_array_values = _max_optical_length * sqrt(1. / (8. * _exp_precision));

  if (num_array_values < 1000)
    num_array_values = 1000;

  if (num_array_values > 1e5) {
    log_printf(WARNING, "Reducing exponential table size from %d to %d",
               num_array_values, 1e5);
    num_array_values = 1e5;
  }

  _exp_table_spacing = _max_optical_length / num_array_values;

  /* Increment the number of vaues in the array to ensure that a tau equal to
   * max_optical_length resides as the final entry in the table */
  num_array_values += 1;

  /* Compute the reciprocal of the table entry spacing */
  _inverse_exp_table_spacing = 1.0 / _exp_table_spacing;

  /* Delete old table */
  if (_exp_table != NULL)
    delete [] _exp_table;

  /* Allocate array for the table */
  _table_size = num_array_values * _num_exp_terms * _num_polar_terms;
  _exp_table = new NEW_PRECISION[_table_size];

  /* Create exponential linear interpolation table */
  for (int i=0; i < num_array_values; i++) {
    for (int p=0; p < _num_polar_terms; p++) {

      int index = _num_exp_terms * (_num_polar_terms * i + p);

      int current_polar = _polar_index + p;
      NEW_PRECISION sin_theta = _quadrature->getSinTheta(azim_index,
                                                        current_polar);
      NEW_PRECISION inv_sin_theta = 1.0 / sin_theta;

      NEW_PRECISION tau_a = i * _exp_table_spacing;
      NEW_PRECISION tau_m = tau_a * inv_sin_theta;
      NEW_PRECISION exponential = exp(-tau_m);

      NEW_PRECISION inv_sin_theta_2 = inv_sin_theta * inv_sin_theta;
      NEW_PRECISION tau_a_2 = tau_a * tau_a;
      NEW_PRECISION sin_theta_2 = sin_theta * sin_theta;

      NEW_PRECISION exp_const_1;
      NEW_PRECISION exp_const_2;
      NEW_PRECISION exp_const_3;

      /* Compute F1 */
      if (tau_a < 0.01) {
        exp_const_1 = inv_sin_theta;
        exp_const_2 = -0.5 * inv_sin_theta * inv_sin_theta;
        exp_const_3 = inv_sin_theta * inv_sin_theta_2 / 6;
        exp_const_1 += exp_const_2 * tau_a + exp_const_3 * tau_a_2;
      }
      else {
        exp_const_1 = (1.0 - exponential) / tau_a;
        exp_const_2 = (exponential * (1.0 + tau_m) - 1) / (tau_a * tau_a);
        exp_const_3 = 0.5 / (tau_a * tau_a * tau_a) *
            (2.0 - exponential * (tau_m * tau_m + 2.0 * tau_m + 2.0));
      }
      
      _exp_table[index] = exp_const_1;
      _exp_table[index+1] = exp_const_2;
      _exp_table[index+2] = exp_const_3;

      if (_linear_source) {

        /* Compute F2 */
        if (tau_a < 0.01) {
          exp_const_2 = inv_sin_theta_2 * inv_sin_theta / 6;
          exp_const_3 = -inv_sin_theta_2 * inv_sin_theta_2 / 12;
          exp_const_1 = exp_const_2 * tau_a + exp_const_3 * tau_a_2;
        }
        else {
          exp_const_1 = (tau_m - 2.0 + exponential * (2.0 + tau_m)) / tau_a_2;
          exp_const_2 = -(tau_m - 4.0 + exponential * (tau_m * tau_m + 3 * tau_m
              + 4.0)) / (tau_a_2 * tau_a);
          exp_const_3 = 0.5 * (2.0 * tau_m - 12.0 + exponential * (12.0 +
              tau_m * (10.0 + tau_m * (4.0 + tau_m)))) / (tau_a_2 * tau_a_2);
        }
        
        _exp_table[index+3] = exp_const_1;
        _exp_table[index+4] = exp_const_2;
        _exp_table[index+5] = exp_const_3;

        /* Compute H */
        if (tau_a < 0.01) {
          exp_const_1 = 0.5 * inv_sin_theta;
          exp_const_2 = -1.0 * inv_sin_theta_2 / 3.0;
          exp_const_3 = inv_sin_theta_2 * inv_sin_theta / 8.0;
          exp_const_1 += exp_const_2 * tau_a + exp_const_3 * tau_a_2;
        }
        else {
          exp_const_1 = (1.0 - exponential * (1 + tau_m)) / (tau_a * tau_m);
          exp_const_2 = (exponential * (tau_m * (tau_m + 2.0) + 2.0) - 2.0) 
              / (tau_m * tau_a_2);
          exp_const_3 = 0.5 * (6.0 - exponential * (6.0 + tau_m * (6.0 + tau_m
              * (3 + tau_m)))) / (tau_m * tau_a_2 * tau_a);
        }
        _exp_table[index+6] = exp_const_1;
        _exp_table[index+7] = exp_const_2;
        _exp_table[index+8] = exp_const_3;
      }
    }
  }
}


/**
 * @brief Computes the G2 exponential term for a optical length and polar angle.
 * @details This method computes the H exponential term from Ferrer [1]
 *          for some optical path length and polar angle. This method
 *          uses either a linear interpolation table (default) or the
 *          exponential intrinsic exp(...) function.
 *
 *            [1] R. Ferrer and J. Rhodes III, "A Linear Source Approximation
 *                Scheme for the Method of Characteristics", Nuclear Science and
 *                Engineering, Volume 182, February 2016.
 *
 * @param tau the optical path length (e.g., sigma_t times length)
 * @param polar the polar angle index
 * @return the evaluated exponential
 */
NEW_PRECISION ExpEvaluator::computeExponentialG2(NEW_PRECISION tau) {

  if (tau == 0.0)
    return 0.0;

  if (tau < 0.01)
    return 7.0 * tau * tau / 120.0 - tau / 12.0;
  else
    return 2.0 / 3.0 - (1 + 2.0 / tau)
        * (1.0 / tau + 0.5 - (1.0 + 1.0 / tau) *
          (1.0 - exp(- tau)) / tau);
}


//FIXME
ExpEvaluator* ExpEvaluator::deepCopy() {

  ExpEvaluator* new_evaluator = new ExpEvaluator();

  if (_interpolate)
    new_evaluator->useInterpolation();
  else
    new_evaluator->useIntrinsic();

  if (_linear_source)
    new_evaluator->useLinearSource();

  new_evaluator->setQuadrature(_quadrature);
  new_evaluator->setMaxOpticalLength(_max_optical_length);
  new_evaluator->setExpPrecision(_exp_precision);

  return new_evaluator;
}
