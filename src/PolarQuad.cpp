#include "PolarQuad.h"


/**
 * @brief Dummy constructor sets the default number of angles to zero.
 */
PolarQuad::PolarQuad() {

  _num_polar = 0;
  _sin_thetas = NULL;
  _weights = NULL;
  _multiples = NULL;

  if (num_polar_angles != 1 && num_polar_angles != 2 && num_polar_angles != 3)
    log_printf(ERROR, "Unable to set the number of polar angles to %d. Only"
               "1, 2 and 3 polar angles are supported by OpenMOC.");

  _num_polar = num_polar;

  /** Initialize memory for arrays */
  _sinthetas = new FP_PRECISION[_num_polar];
  _weights = new FP_PRECISION[_num_polar];
  _multiples = new FP_PRECISION[_num_polar];

  /* If Quadrature type is TABUCHI */
  if (type == TABUCHI) {

    _type = TABUCHI;

    if (_num_polar == 1) {
      _sinthetas[0] = 0.798184;
      _weights[0] = 1.0;
      _multiples[0] = _sinthetas[0] * _weights[0];
    }

    else if (_num_polar == 2) {
      _sinthetas[0] = 0.363900;
      _sinthetas[1] = 0.899900;
      _weights[0] = 0.212854;
      _weights[1] = 0.787146;
      _multiples[0] = _sinthetas[0] * _weights[0];
      _multiples[1] = _sinthetas[1] * _weights[1];
    }

    else if (_num_polar == 3) {
      _sinthetas[0] = 0.166648;
      _sinthetas[1] = 0.537707;
      _sinthetas[2] = 0.932954;
      _weights[0] = 0.046233;
      _weights[1] = 0.283619;
      _weights[2] = 0.670148;
      _multiples[0] = _sinthetas[0] * _weights[0];
      _multiples[1] = _sinthetas[1] * _weights[1];
      _multiples[2] = _sinthetas[2] * _weights[2];
    }

    else
      log_printf(ERROR, "TABUCHI type polar Quadrature supports 1, 2, or "
                 "3 polar angles but %d are defined", _num_polar);
  }

  /* If quadrature type is LEONARD */
  else if (type == LEONARD) {

    _type = LEONARD;

    if (_num_polar == 2) {
      _sinthetas[0] = 0.273658;
      _sinthetas[1] = 0.865714;
      _weights[0] = 0.139473;
      _weights[1] = 0.860527;
      _multiples[0] = _sinthetas[0] * _weights[0];
      _multiples[1] = _sinthetas[1] * _weights[1];
    }

    else if (_num_polar == 3) {
      _sinthetas[0] = 0.099812;
      _sinthetas[1] = 0.395534;
      _sinthetas[2] = 0.891439;
      _weights[0] = 0.017620;
      _weights[1] = 0.188561;
      _weights[2] = 0.793819;
      _multiples[0] = _sinthetas[0] * _weights[0];
      _multiples[1] = _sinthetas[1] * _weights[1];
      _multiples[2] = _sinthetas[2] * _weights[2];
    }

    else {
      log_printf(ERROR, "LEONARD type polar Quadrature supports 2, or 3"
                 "polar angles but %d are defined", _num_polar);
    }

  }

  else
    log_printf(ERROR, "LEONARD and TABUCHI polar Quadrature types "
               "supported, but unknown type given");
}


/**
 * @brief Destructor deletes arrray of sines of the polar angles, the weights
 *        of the polar angles and the products of the sines and weights.
 */
PolarQuad::~PolarQuad() {

  if (_sin_thetas != NULL)
    delete [] _sinthetas;

  if (_weights != NULL)
    delete [] _weights;

  if (_multiples != NULL)
    delete [] _multiples;
}



/**
 * @brief Returns the number of polar angles.
 * @return the number of polar angles
 */
int PolarQuad::getNumPolarAngles() const {
  return _num_polar;
}


/**
 * @brief Returns the \f$ sin(\theta)\f$ value for a particular polar angle.
 * @param n index of the polar angle of interest
 * @return the value of \f$ \sin(\theta) \f$ for this polar angle
 */
FP_PRECISION PolarQuad::getSinTheta(const int n) const {

  if (n < 0 || n >= _num_polar)
    log_printf(ERROR, "Attempted to retrieve sin theta for polar angle = %d"
               " but only %d polar angles are defined", n, _num_polar);

  return _sin_thetas[n];
}


/**
 * @brief Returns the weight value for a particular polar angle.
 * @param n index of the polar angle of interest
 * @return the weight for a polar angle
 */
FP_PRECISION PolarQuad::getWeight(const int n) const {

  if (n < 0 || n >= _num_polar)
    log_printf(ERROR, "Attempted to retrieve the weight for polar angle = %d "
               "but only %d polar angles are defined", n, _num_polar);

  return _weights[n];
}


/**
 * @brief Returns the multiple value for a particular polar angle.
 * @details A multiple is the sine of a polar angle multiplied by its weight.
 * @param n index of the polar angle of interest
 * @return the value of the sine of the polar angle multiplied with its weight
 */
FP_PRECISION PolarQuad::getMultiple(const int n) const {

  if (n < 0 || n >= _num_polar)
    log_printf(ERROR, "Attempted to retrieve the multiple for polar angle = %d"
               "but only %d polar angles are defined", n, _num_polar);

  return _multiples[n];
}


/**
 * @brief Returns a pointer to the PolarQuad's array of \f$ sin\theta_{p} \f$.
 * @return a pointer to the array of \f$ sin\theta_{p} \f$
 */
FP_PRECISION* PolarQuad::getSinThetas() {
  return _sin_thetas;
}


/**
 * @brief Returns a pointer to the PolarQuad's array of polar weights.
 * @return a pointer to the polar weights array
 */
FP_PRECISION* PolarQuad::getWeights() {
  return _weights;
}


/**
 * @brief Returns a pointer to the PolarQuad's array of multiples.
 * @details A multiple is the sine of a polar angle multiplied by its weight.
 * @return a pointer to the multiples array
 */
FP_PRECISION* PolarQuad::getMultiples() {
  return _multiples;
}


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles
 */
void PolarQuad::setNumPolarAngles(const int num_polar) {

  if (num_polar <= 0)
    log_printf(ERROR, "Unable to set the number of polar angles to %d "
               "which is less than or equal to zero", num_polar);

  _num_polar = num_polar;
}


/**
 * @brief Set the PolarQuad's array of sines of each polar angle.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the PolarQuad's sin thetas in Python. A user must initialize a
 *          NumPy array of the correct size (i.e., a float64 array the length
 *          of the number of polar angles) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          PolarQuad's sin thetas. An example of how this function might be
 *          called in Python is as follows:
 *
 * @code
 *          sin_thetas= numpy.array([0.05, 0.1, 0.15, ... ])
 *          polar_quad = openmoc.PolarQuad()
 *          polar_quad.setNumPolarAngles(len(sin_thetas))
 *          polar_quad.setSinThetas(sin_thetas)
 * @endcode
 *
 * @param sin_thetas the array of sines of each polar angle
 * @param num_polar the number of polar angles
 */
void PolarQuad::setSinThetas(FP_PRECISION* sin_thetas, int num_polar) {

  if (_num_polar != num_polar)
    log_printf(ERROR, "Unable to set %d sin thetas for PolarQuad "
               "with %d polar angles", num_polar, _num_polar);

  for (int i=0; i < _num_polar; i++)
    _sin_thetas[i] = FP_PRECISION(sin_thetas[i]);

}


/**
 * @brief Set the PolarQuad's array of weights for each angle.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the PolarQuad's angular weights in Python. A user must initialize a
 *          NumPy array of the correct size (i.e., a float64 array the length
 *          of the number of polar angles) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          PolarQuad's weights. An example of how this function might be
 *          called in Python is as follows:
 *
 * @code
 *          weights = numpy.array([0.05, 0.1, 0.15, ... ])
 *          polar_quad = openmoc.PolarQuad()
 *          polar_quad.setNumPolarAngles(len(sin_thetas))
 *          polar_quad.setWeights(weights)
 * @endcode
 *
 * @param weights the array of weights for each polar angle
 * @param num_polar the number of polar angles
 */
void PolarQuad::setWeights(FP_PRECISION* weights, int num_polar) {

  if (_num_polar != num_polar)
    log_printf(ERROR, "Unable to set %d weights for PolarQuad "
               "with %d polar angles", num_polar, _num_polar);

  for (int i=0; i < _num_polar; i++)
    _weights[i] = FP_PRECISION(weights[i]);
}

void PolarQuad::initialize() {
}


/**
 * @brief Converts this Quadrature to a character array of its attributes.
 * @details The character array includes the number of polar angles, the
 *          the values of the sine and weight of each polar angle, and the 
 *          product of the sine and weight of each polar angle.
 * @return a character array of the PolarQuad's attributes
 */
std::string PolarQuad::toString() {

  std::stringstream string;

  string << "PolarQuad";
  string << "\n\t# angles = " << _num_polar;

  string << "\n\tsin thetas = ";
  if (_sin_thetas != NULL) {
    for (int p = 0; p < _num_polar; p++)
      string << _sinthetas[p] << ", ";
  }

  string << "\n\tweights = ";
  if (_weights != NULL) {
    for (int p = 0; p < _num_polar; p++)
      string << _weights[p] << ", ";
  }

  string << "\n\tmultiples= ";
  if (_multiples != NULL) {
    for (int p = 0; p < _num_polar; p++)
      string << _multiples[p] << ", ";
  }

  return string.str();
}
