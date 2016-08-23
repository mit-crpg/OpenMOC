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
  _solve_3D = false;
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


void ExpEvaluator::setSolve3D(bool solve_3D) {
  _solve_3D = solve_3D;
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


bool ExpEvaluator::isSolve3D() {
  return _solve_3D;
}


/**
 * @brief If using linear interpolation, builds the table for each polar angle.
 //FIXME
 */
void ExpEvaluator::initialize(int azim_index, int polar_index, bool solve_3D) {

  /* Check for a valid Quadrature */
  if (_quadrature == NULL)
    log_printf(ERROR, "A Quadrature must be set before an Exponential "
               "Evaluator can be initialized");

  /* Extract the number of azimuthal and polar angles */
  //FIXME

  //FIXME
  if (azim_index < 0 || azim_index > _quadrature->get) return;

  /* Record the azimuthal and polar angle indexes */
  _azim_index = azim_index;

  /* Record the inverse sine for the base polar angle */
  _inv_sin_theta_no_offset = 1.0 / _quadrature->getSinTheta(azim_index,
                                                            polar_index);
  /* If no exponential table is needed, return */
  if (!_interpolate)
    return;

  log_printf(INFO, "Initializing exponential interpolation table...");

  int num_polar = _quadrature->getNumPolarAngles();
  int _max_polar_offset;
  if (solve_3D) {
    _polar_index = polar_index;
    _max_polar_offset = 1;
  }
  else {
    _polar_index = 0;
    _max_polar_offset = num_polar / 2;
  }

  /* Set size of interpolation table */
  int num_array_values;
  if (_linear_source)
    num_array_values = _max_optical_length * sqrt(1. / (8. * _exp_precision));
  else
    num_array_values = _max_optical_length * pow(1. / (72. * sqrt(3.0)
                                                 * _exp_precision), 1.0/3.0);

  FP_PRECISION _exp_table_spacing = _max_optical_length / num_array_values;

  /* Increment the number of vaues in the array to ensure that a tau equal to
   * max_optical_length resides as the final entry in the table */
  num_array_values += 1;

  /* Compute the reciprocal of the table entry spacing */
  _inverse_exp_table_spacing = 1.0 / _exp_table_spacing;

  /* Allocate array for the table */
  if (_exp_table != NULL)
    delete [] _exp_table;

  if (_solve_3D) {
    _table_size = num_array_values * _num_exp_terms * _max_polar_offset;
    _exp_table = new FP_PRECISION[_table_size];

    /* Create exponential linear interpolation table */
    for (int i=0; i < num_array_values; i++) {
      for (int p=0; p < _max_polar_offset; p++) {

        int index = _num_exp_terms * (_max_polar_offset * i + p);

        int current_polar = _polar_index + p;
        FP_PRECISION sin_theta = _quadrature->getSinTheta(azim_index,
                                                          current_polar);
        FP_PRECISION inv_sin_theta = 1.0 / sin_theta;

        FP_PRECISION tau_a = i * _exp_table_spacing;
        FP_PRECISION tau_m = tau_a * inv_sin_theta;
        FP_PRECISION exponential = exp(-tau_m);

        FP_PRECISION inv_sin_theta_2 = inv_sin_theta * inv_sin_theta;
        FP_PRECISION tau_a_2 = tau_a * tau_a;
        FP_PRECISION sin_theta_2 = sin_theta * sin_theta;

        /* Compute F1 */
        FP_PRECISION exp_const_1 = 1.0 - expon;
        FP_PRECISION exp_const_2 = exponential * inv_sin_theta;
        FP_PRECISION exp_const_3 = -0.5 * exp_const_2 * inv_sin_theta;

        _exp_table[index] = exp_const_1;
        _exp_table[index+1] = exp_const_2;
        _exp_table[index+2] = exp_const_3;

        if (_linear_source) {

          /* Compute F2 */
          exp_const_1 = 2 * exponential - 2 + tau_m + tau_m * exponential;
          exp_const_2 = (-exponential * (tau_a + sin_theta) + sin_theta) *
              inv_sin_theta_2;
          exp_const_3 = 0.5 * tau_a * exponential * inv_sin_theta *
              inv_sin_theta_2;

          _exp_table[index+3] = exp_const_1;
          _exp_table[index+4] = exp_const_2;
          _exp_table[index+5] = exp_const_3;

          /* Compute H */
          if (tau_a == 0.0) {
            exp_const_1 = 0.0;
            exp_const_2 = 0.5 * inv_sin_theta;
            exp_const_3 = -1.0 * inv_sin_theta_2 / 3.0;
          }
          else {
            exp_const_1 = (-exponential * (tau_a + sin_theta) * sin_theta) /
                tau_a;
            exp_const_2 = (exponential * (tau_a_2 + tau_a * sin_theta +
                sin_theta_2) - sin_theta_2) / (tau_a_2 * sin_theta);
            exp_const_3 = 1.0 / (2 * tau_a_2 * tau_a * sin_theta_2) *
                (-exponential * tau_a_2 + tau_a_2 * sin_theta + 2 * tau_a *
                sin_theta_2 + 2 * sin_theta_2 * sin_theta) + 2 * sin_theta_2
                * sin_theta;
          }
          _exp_table[index+6] = exp_const_1;
          _exp_table[index+7] = exp_const_2;
          _exp_table[index+8] = exp_const_3;
        }
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
FP_PRECISION ExpEvaluator::computeExponentialG2(FP_PRECISION tau) {

  if (tau == 0.0)
    return 0.0;

  tau = std::max(tau, 1.e-5);

  return 2.0 * tau / 3.0 - (1 + 2.0 / tau)
      * (1.0 + tau / 2.0 - (1.0 + 1.0 / tau) *
         (1.0 - exp(- tau)));
}
