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


void ExpEvaluator::setSolve3D(bool solve_3D){
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


bool ExpEvaluator::isSolve3D(){
  return _solve_3D;
}


/**
 * @brief If using linear interpolation, builds the table for each polar angle.
 * @param max_tau the maximum optical path length in the input range
 * @param tolerance the minimum acceptable interpolation accuracy
 */
void ExpEvaluator::initialize(double max_tau, double tolerance) {

  /* If no exponential table is needed, return */
  if (!_interpolate)
    return;

  log_printf(INFO, "Initializing exponential interpolation table...");

  /* Expand max tau slightly to avoid roundoff error approximation */
  max_tau *= 1.01;

  /* Set size of interpolation table */
  int num_array_values = max_tau * sqrt(1. / (8.e-2 * tolerance));
  FP_PRECISION exp_table_spacing = max_tau / num_array_values;

  /* Compute the reciprocal of the table entry spacing */
  _inverse_exp_table_spacing = 1.0 / exp_table_spacing;

  /* Allocate array for the table */
  if (_exp_table != NULL)
    delete [] _exp_table;

  if (_solve_3D){
    _table_size = num_array_values * 2;
    _exp_table = new FP_PRECISION[_table_size];
    
    FP_PRECISION expon;
    FP_PRECISION intercept;
    FP_PRECISION slope;
    FP_PRECISION sin_theta;
    
    /* Create exponential linear interpolation table */
    for (int i=0; i < num_array_values; i++){
      expon = exp(- (i * exp_table_spacing));
      slope = - expon;
      intercept = expon * (1 + (i * exp_table_spacing));
      _exp_table[i * 2] = slope;
      _exp_table[i * 2 + 1] = intercept;
    }
  }
  else{
    _table_size = _num_polar * num_array_values;
    _exp_table = new FP_PRECISION[_table_size];
    
    FP_PRECISION expon;
    FP_PRECISION intercept;
    FP_PRECISION slope;
    FP_PRECISION sin_theta;
    
    /* Create exponential linear interpolation table */
    for (int i=0; i < num_array_values; i++){
      for (int p=0; p < _num_polar/2; p++){
        sin_theta = _quadrature->getSinTheta(0,p);
        expon = exp(- (i * exp_table_spacing) / sin_theta);
        slope = - expon / sin_theta;
        intercept = expon * (1 + (i * exp_table_spacing) / sin_theta);
        _exp_table[_num_polar * i + 2 * p] = slope;
        _exp_table[_num_polar * i + 2 * p + 1] = intercept;
      }
    }
  }
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
FP_PRECISION ExpEvaluator::computeExponential(FP_PRECISION tau, int azim, int polar) {

  FP_PRECISION exponential;

  /* Evaluate the exponential using the lookup table - linear interpolation */
  if (_interpolate) {
    if (_solve_3D){
      int index;
      index = round_to_int(tau * _inverse_exp_table_spacing) * 2;
      exponential = (1. - (_exp_table[index] * tau +
                           _exp_table[index + 1]));
    }
    else{
      int index;
      index = round_to_int(tau * _inverse_exp_table_spacing);
      index *= _num_polar;
      exponential = (1. - (_exp_table[index + 2 * polar] * tau +
                           _exp_table[index + 2 * polar + 1]));
    }
  }

  /* Evalute the exponential using the intrinsic exp(...) function */
  else {
    if (_solve_3D)
      exponential = 1.0 - exp(- tau);
    else{
      FP_PRECISION sintheta = _quadrature->getSinTheta(azim, polar);
      exponential = 1.0 - exp(- tau / sintheta);
    }
  }

  return exponential;
}
