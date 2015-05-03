#include "PolarQuad.h"


/**
 * @brief Dummy constructor sets the default number of angles to zero.
 */
PolarQuad::PolarQuad() {
  _num_polar = 0;
  _sin_thetas = NULL;
  _weights = NULL;
  _multiples = NULL;
}


/**
 * @brief Destructor deletes arrray of sines of the polar angles, the weights
 *        of the polar angles and the products of the sines and weights.
 */
PolarQuad::~PolarQuad() {

  if (_sin_thetas != NULL)
    delete [] _sin_thetas;

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
    log_printf(ERROR, "Attempted to retrieve sin theta for polar angle = "
               "%d but only %d polar angles are defined", n, _num_polar);

  else if (_sin_thetas == NULL)
    log_printf(ERROR, "Attempted to retrieve sin theta for polar angle = %d "
               "but the sin thetas have not been initialized", n);

  return _sin_thetas[n];
}


/**
 * @brief Returns the weight value for a particular polar angle.
 * @param n index of the polar angle of interest
 * @return the weight for a polar angle
 */
FP_PRECISION PolarQuad::getWeight(const int n) const {

  if (n < 0 || n >= _num_polar)
    log_printf(ERROR, "Attempted to retrieve the weight for polar angle = "
               "%d but only %d polar angles are defined", n, _num_polar);

  else if (_weights == NULL)
    log_printf(ERROR, "Attempted to retrieve weight for polar angle = %d "
               "but the weights have not been initialized", n);

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
    log_printf(ERROR, "Attempted to retrieve the multiple for polar angle = "
               "%d but only %d polar angles are defined", n, _num_polar);

  else if (_multiples == NULL)
    log_printf(ERROR, "Attempted to retrieve multiple for polar angle = %d "
               "but the multiples have not been initialized", n);

  return _multiples[n];
}


/**
 * @brief Returns a pointer to the PolarQuad's array of \f$ sin\theta_{p} \f$.
 * @return a pointer to the array of \f$ sin\theta_{p} \f$
 */
FP_PRECISION* PolarQuad::getSinThetas() {

  if (_sin_thetas == NULL)
    log_printf(ERROR, "Attempted to retrieve the sin thetas array "
               "but it has not been initialized");

  return _sin_thetas;
}


/**
 * @brief Returns a pointer to the PolarQuad's array of polar weights.
 * @return a pointer to the polar weights array
 */
FP_PRECISION* PolarQuad::getWeights() {

  if (_weights == NULL)
    log_printf(ERROR, "Attempted to retrieve the weights array "
               "but it has not been initialized");

  return _weights;
}


/**
 * @brief Returns a pointer to the PolarQuad's array of multiples.
 * @details A multiple is the sine of a polar angle multiplied by its weight.
 * @return a pointer to the multiples array
 */
FP_PRECISION* PolarQuad::getMultiples() {

  if (_multiples == NULL)
    log_printf(ERROR, "Attempted to retrieve the multiples array "
               "but it has not been initialized");

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
 *          NumPy array of the correct size (e.g., a float64 array the length
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
void PolarQuad::setSinThetas(double* sin_thetas, int num_polar) {

  if (_num_polar != num_polar)
    log_printf(ERROR, "Unable to set %d sin thetas for PolarQuad "
               "with %d polar angles", num_polar, _num_polar);

  /* Deallocate memory if it was allocated previously */
  if (_sin_thetas != NULL)
    delete [] _sin_thetas;

  /* Initialize memory for arrays */
  _sin_thetas = new FP_PRECISION[_num_polar];

  /* Extract sin thetas from user input */
  for (int p=0; p < _num_polar; p++) {
    if (sin_thetas[p] < 0. || sin_thetas[p] > 1.)
      log_printf(ERROR, "Unable to set sin theta to %f which is "
                 "not in the range [0,1]", sin_thetas[p]);
    _sin_thetas[p] = FP_PRECISION(sin_thetas[p]);
  }
}


/**
 * @brief Set the PolarQuad's array of weights for each angle.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the PolarQuad's angular weights in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 array the length
 *          of the number of polar angles) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          PolarQuad's weights. An example of how this function might be
 *          called in Python is as follows:
 *
 * @code
 *          weights = numpy.array([0.05, 0.1, 0.15, ... ])
 *          polar_quad = openmoc.PolarQuad()
 *          polar_quad.setNumPolarAngles(len(weights))
 *          polar_quad.setWeights(weights)
 * @endcode
 *
 * @param weights the array of weights for each polar angle
 * @param num_polar the number of polar angles
 */
void PolarQuad::setWeights(double* weights, int num_polar) {

  if (_num_polar != num_polar)
    log_printf(ERROR, "Unable to set %d weights for PolarQuad "
               "with %d polar angles", num_polar, _num_polar);

  /* Deallocate memory if it was allocated previously */
  if (_weights != NULL)
    delete [] _weights;

  /* Initialize memory for arrays */
  _weights = new FP_PRECISION[_num_polar];
  double tot_weight = 0.;

  /* Extract weights from user input */
  for (int p=0; p < _num_polar; p++) {
    _weights[p] = FP_PRECISION(weights[p]);
    tot_weight += _weights[p];
  }  

  /* Check that weights sum to 1 (half of the symmetric polar domain) */
  if (fabs(tot_weight - 1.0) > POLAR_WEIGHT_SUM_TOL)
    log_printf(ERROR, "The polar weights sum to %f instead of 1.0", tot_weight);
}


/**
 * @brief Dummy routine to initialize the polar quadrature.
 * @details The parent class routine simply checks that the number of polar
 *          angles has been set by the user and returns;
 */
void PolarQuad::initialize() {
  if (_num_polar == 0)
    log_printf(ERROR, "Unable to initialize PolarQuad with zero polar angles. "
               "Set the number of polar angles before initialization.");

  if (_sin_thetas != NULL && _weights != NULL)
    precomputeMultiples();
}


/**
 * @brief This private routine computes the produce of the sine thetas and
 *        weights for each angle in the polar quadrature.
 * @details Note that this routine must be called after populating the
 *          sine thetas and weights arrays.
 */
void PolarQuad::precomputeMultiples() {

  /* Deallocate memory if it was allocated previously */
  if (_multiples != NULL)
    delete [] _multiples;

  /* Initialize memory for arrays */
  _multiples = new FP_PRECISION[_num_polar];

  /* Compute multiples of sine thetas and weights */
  for (int p=0; p < _num_polar; p++)
    _multiples[p] = _sin_thetas[p] * _weights[p];
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
      string << _sin_thetas[p] << ", ";
  }

  string << "\n\tweights = ";
  if (_weights != NULL) {
    for (int p = 0; p < _num_polar; p++)
      string << _weights[p] << ", ";
  }

  string << "\n\tmultiples = ";
  if (_multiples != NULL) {
    for (int p = 0; p < _num_polar; p++)
      string << _multiples[p] << ", ";
  }

  return string.str();
}



/**
 * @brief Dummy constructor calls the parent constructor.
 */
TYPolarQuad::TYPolarQuad(): PolarQuad() { }


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles (maximum 3)
 */
void TYPolarQuad::setNumPolarAngles(const int num_polar) {

  if (num_polar > 3)
    log_printf(ERROR, "Unable to set the number of polar angles to %d "
               "for TYPolarQuad (max 3 angles)", num_polar);

  PolarQuad::setNumPolarAngles(num_polar);
}


/**
 * @brief Routine to initialize the polar quadrature.
 * @details This routine uses the tabulated values for the Tabuchi-Yamamoto
 *          polar angle quadrature, including the sine thetas and weights.
 */
void TYPolarQuad::initialize() {

  /* Call parent class initialize routine */
  PolarQuad::initialize();

  /* Allocate temporary arrays for tabulated quadrature values */
  double* sin_thetas = new double[_num_polar];
  double* weights = new double[_num_polar];

  /* Tabulated values for the sine thetas and weights for the
   * Tabuchi-Yamamoto polar angle quadrature */
  if (_num_polar == 1) {
    sin_thetas[0] = 0.798184;
    weights[0] = 1.0;
  }

  else if (_num_polar == 2) {
    sin_thetas[0] = 0.363900;
    sin_thetas[1] = 0.899900;
    weights[0] = 0.212854;
    weights[1] = 0.787146;
  }

  else if (_num_polar == 3) {
    sin_thetas[0] = 0.166648;
    sin_thetas[1] = 0.537707;
    sin_thetas[2] = 0.932954;
    weights[0] = 0.046233;
    weights[1] = 0.283619;
    weights[2] = 0.670148;
  }

  /* Set the arrays of sin thetas and weights */
  PolarQuad::setSinThetas(sin_thetas, _num_polar);
  PolarQuad::setWeights(weights, _num_polar);

  /* Deallocate temporary arrays */
  delete [] sin_thetas;
  delete [] weights;

  /* Compute the product of the sine thetas and weights */
  PolarQuad::precomputeMultiples();
}


/**
 * @brief Dummy constructor calls the parent constructor.
 */
LeonardPolarQuad::LeonardPolarQuad(): PolarQuad() { }


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles (2 or 3)
 */
void LeonardPolarQuad::setNumPolarAngles(const int num_polar) {

  if (num_polar < 2 || num_polar > 3)
    log_printf(ERROR, "Unable to set the number of polar angles to %d "
               "for LeonardPolarQuad (2 or 3 angles)", num_polar);

  PolarQuad::setNumPolarAngles(num_polar);
}


/**
 * @brief Routine to initialize the polar quadrature.
 * @details This routine uses the tabulated values for the Leonard
 *          polar angle quadrature, including the sine thetas and weights.
 */
void LeonardPolarQuad::initialize() {

  /* Call parent class initialize routine */
  PolarQuad::initialize();

  /* Allocate temporary arrays for tabulated quadrature values */
  double* sin_thetas = new double[_num_polar];
  double* weights = new double[_num_polar];

  /* Tabulated values for the sine thetas and weights for the
   * Leonard polar angle quadrature */
  if (_num_polar == 2) {
    sin_thetas[0] = 0.273658;
    sin_thetas[1] = 0.865714;
    weights[0] = 0.139473;
    weights[1] = 0.860527;
  }

  else if (_num_polar == 3) {
    sin_thetas[0] = 0.099812;
    sin_thetas[1] = 0.395534;
    sin_thetas[2] = 0.891439;
    weights[0] = 0.017620;
    weights[1] = 0.188561;
    weights[2] = 0.793819;
  }

  /* Set the arrays of sin thetas and weights */
  PolarQuad::setSinThetas(sin_thetas, _num_polar);
  PolarQuad::setWeights(weights, _num_polar);

  /* Deallocate temporary arrays */
  delete [] sin_thetas;
  delete [] weights;

  /* Compute the product of the sine thetas and weights */
  PolarQuad::precomputeMultiples();
}



/**
 * @brief Dummy constructor calls the parent constructor.
 */
GLPolarQuad::GLPolarQuad(): PolarQuad() { }


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles (maximum 6)
 */
void GLPolarQuad::setNumPolarAngles(const int num_polar) {

  if (num_polar > 6)
    log_printf(ERROR, "Unable to set the number of polar angles to %d "
               "for GLPolarQuad (max 6 angles)", num_polar);

  PolarQuad::setNumPolarAngles(num_polar);
}


/**
 * @brief Routine to initialize the polar quadrature.
 * @details This routine uses the tabulated values for the Gauss-Legendre
 *          polar angle quadrature, including the sine thetas and weights.
 */
void GLPolarQuad::initialize() {

  /* Call parent class initialize routine */
  PolarQuad::initialize();

  /* Allocate temporary arrays for tabulated quadrature values */
  double* sin_thetas = new double[_num_polar];
  double* weights = new double[_num_polar];

  /* Tabulated values for the sine thetas and weights for the
   * Leonard polar angle quadrature */
  if (_num_polar == 1) {
    sin_thetas[0] = sin(acos(0.5773502691));
    weights[0] = 1.0;
  }
  else if (_num_polar == 2) {
    sin_thetas[0] = sin(acos(0.3399810435));
    sin_thetas[1] = sin(acos(0.8611363115));
    weights[0] = 0.6521451549;
    weights[1] = 0.3478548451;
  }
  else if (_num_polar == 3) {
    sin_thetas[0] = sin(acos(0.2386191860));
    sin_thetas[1] = sin(acos(0.6612093864));
    sin_thetas[2] = sin(acos(0.9324695142));
    weights[0] = 0.4679139346;
    weights[1] = 0.3607615730;
    weights[2] = 0.1713244924;
  }
  else if (_num_polar == 4) {
    sin_thetas[0] = sin(acos(0.1834346424));
    sin_thetas[1] = sin(acos(0.5255324099));
    sin_thetas[2] = sin(acos(0.7966664774));
    sin_thetas[3] = sin(acos(0.9602898564));
    weights[0] = 0.3626837834;
    weights[1] = 0.3137066459;
    weights[2] = 0.2223810344;
    weights[3] = 0.1012285363;
  }
  else if (_num_polar == 5) {
    sin_thetas[0] = sin(acos(0.1488743387));
    sin_thetas[1] = sin(acos(0.4333953941));
    sin_thetas[2] = sin(acos(0.6794095682));
    sin_thetas[3] = sin(acos(0.8650633666));
    sin_thetas[4] = sin(acos(0.9739065285));
    weights[0] = 0.2955242247;
    weights[1] = 0.2692667193;
    weights[2] = 0.2190863625;
    weights[3] = 0.1494513492;
    weights[4] = 0.0666713443;
  }
  else if (_num_polar == 6) {
    sin_thetas[0] = sin(acos(0.1252334085));
    sin_thetas[1] = sin(acos(0.3678314989));
    sin_thetas[2] = sin(acos(0.5873179542));
    sin_thetas[3] = sin(acos(0.7699026741));
    sin_thetas[4] = sin(acos(0.9041172563));
    sin_thetas[5] = sin(acos(0.9815606342));
    weights[0] = 0.2491470458;
    weights[1] = 0.2334925365;
    weights[2] = 0.2031674267;
    weights[3] = 0.1600783286;
    weights[4] = 0.1069393260;
    weights[5] = 0.0471753364;
  }

  /* Set the arrays of sin thetas and weights */
  PolarQuad::setSinThetas(sin_thetas, _num_polar);
  PolarQuad::setWeights(weights, _num_polar);

  /* Deallocate temporary arrays */
  delete [] sin_thetas;
  delete [] weights;

  /* Compute the product of the sine thetas and weights */
  PolarQuad::precomputeMultiples();
}



/**
 * @brief Dummy constructor calls the parent constructor.
 */
EqualWeightPolarQuad::EqualWeightPolarQuad(): PolarQuad() { }


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles
 */
void EqualWeightPolarQuad::setNumPolarAngles(const int num_polar) {
  PolarQuad::setNumPolarAngles(num_polar);
}


/**
 * @brief Routine to initialize the polar quadrature.
 * @details This routine generates the sine thetas and weights.
 */
void EqualWeightPolarQuad::initialize() {

  /* Call parent class initialize routine */
  PolarQuad::initialize();

  /* Allocate temporary arrays for tabulated quadrature values */
  double* sin_thetas = new double[_num_polar];
  double* weights = new double[_num_polar];

  double cos_theta_a, cos_theta_b;
  cos_theta_a = 1.;

  /* Generate the sin thetas and weights using equations 420-422 of the
   * DOE Nucl. Eng. Handbook "Lattice Physics Computations" */
  for (int p=0; p < _num_polar; p++) {
    cos_theta_b = cos_theta_a - (1. / _num_polar);
    sin_thetas[p] = sin(acos(0.5 * (cos_theta_a + cos_theta_b)));
    weights[p] = cos_theta_a - cos_theta_b;
    cos_theta_a = cos_theta_b;
  }

  /* Set the arrays of sin thetas and weights */
  PolarQuad::setSinThetas(sin_thetas, _num_polar);
  PolarQuad::setWeights(weights, _num_polar);

  /* Deallocate temporary arrays */
  delete [] sin_thetas;
  delete [] weights;

  /* Compute the product of the sine thetas and weights */
  PolarQuad::precomputeMultiples();
}



/**
 * @brief Dummy constructor calls the parent constructor.
 */
EqualAnglePolarQuad::EqualAnglePolarQuad(): PolarQuad() { }


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles
 */
void EqualAnglePolarQuad::setNumPolarAngles(const int num_polar) {
  PolarQuad::setNumPolarAngles(num_polar);
}


/**
 * @brief Routine to initialize the polar quadrature.
 * @details This routine generates the sine thetas and weights.
 */
void EqualAnglePolarQuad::initialize() {

  /* Call parent class initialize routine */
  PolarQuad::initialize();

  /* Allocate temporary arrays for tabulated quadrature values */
  double* sin_thetas = new double[_num_polar];
  double* weights = new double[_num_polar];

  double cos_theta_a, cos_theta_b;
  double theta_a, theta_b;
  double delta_theta = M_PI / (2. * _num_polar);
  theta_a = 0.;

  /* Generate the sin thetas and weights using equations 420-422 of the
   * DOE Nucl. Eng. Handbook "Lattice Physics Computations" */
  for (int p=0; p < _num_polar; p++) {
    theta_b = theta_a + delta_theta;
    cos_theta_a = cos(theta_a);
    cos_theta_b = cos(theta_b);
    sin_thetas[p] = sin(acos((0.5 * (cos_theta_a + cos_theta_b))));
    weights[p] = cos_theta_a - cos_theta_b;
    theta_a = theta_b;
  }

  /* Set the arrays of sin thetas and weights */
  PolarQuad::setSinThetas(sin_thetas, _num_polar);
  PolarQuad::setWeights(weights, _num_polar);

  /* Deallocate temporary arrays */
  delete [] sin_thetas;
  delete [] weights;

  /* Compute the product of the sine thetas and weights */
  PolarQuad::precomputeMultiples();
}
