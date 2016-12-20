#include "Quadrature.h"


/**
 * @brief Dummy constructor sets the default number of angles to zero.
 */
Quadrature::Quadrature() {
  _num_azim = 0;
  _num_polar = 0;
  _sin_thetas = NULL;
  _thetas = NULL;
  _phis = NULL;
  _azim_weights = NULL;
  _polar_weights = NULL;
  _multiples = NULL;
}


/**
 * @brief Destructor deletes arrray of sines of the polar angles, the weights
 *        of the polar angles and the products of the sines and weights.
 */
Quadrature::~Quadrature() {

  if (_sin_thetas != NULL) {
    for (int i=0; i < _num_azim/4; i++)
      delete [] _sin_thetas[i];
    delete [] _sin_thetas;
  }

  if (_thetas != NULL) {
    for (int i=0; i < _num_azim/4; i++)
      delete [] _thetas[i];
    delete [] _thetas;
  }

  if (_phis != NULL)
    delete [] _phis;

  if (_azim_weights != NULL)
    delete [] _azim_weights;

  if (_polar_weights != NULL) {
    for (int i=0; i < _num_azim/4; i++)
      delete [] _polar_weights[i];
    delete [] _polar_weights;
  }

  if (_multiples != NULL) {
    for (int i=0; i < _num_azim/4; i++)
      delete [] _multiples[i];
    delete [] _multiples;
  }
}



/**
 * @brief Returns the number of polar angles.
 * @return the number of polar angles
 */
int Quadrature::getNumPolarAngles() const {
  return _num_polar;
}


/**
 * @brief Returns the number of azimuthal angles.
 * @return the number of azimuthal angles
 */
int Quadrature::getNumAzimAngles() const {
  return _num_azim;
}


/**
 * @brief Returns the \f$ sin(\theta)\f$ value for a particular polar angle.
 * @param azim index of the azimthal angle of interest
 * @param polar index of the polar angle of interest
 * @return the value of \f$ \sin(\theta) \f$ for this azimuthal and polar angle
 */
FP_PRECISION Quadrature::getSinTheta(int azim, int polar) {

  if (polar < 0 || polar >= _num_polar)
    log_printf(ERROR, "Attempted to retrieve sin theta for polar angle = "
               "%d and azim angle = %d when only %d polar angles are "
               "defined", polar, azim, _num_polar);

  else if (azim < 0 || azim >= _num_azim)
    log_printf(ERROR, "Attempted to retrieve sin theta for polar angle = "
               "%d and azim angle = %d when only %d azim angles are "
               "defined", polar, azim, _num_azim);

  else if (_sin_thetas == NULL)
    log_printf(ERROR, "Attempted to retrieve sin theta for polar angle = %d "
               "and azim angle = %d but the sin thetas have not been "
               "initialized", polar, azim);

  azim = getFirstOctantAzim(azim);
  polar = getFirstOctantPolar(polar);

  return _sin_thetas[azim][polar];
}


/**
 * @brief Returns the polar angle value for a given azimuthal and polar
 *        angle index.
 * @param azim index of the azimthal angle of interest
 * @param polar index of the polar angle of interest
 * @return the value of the polar angle for this azimuthal and polar angle index
 */
double Quadrature::getTheta(int azim, int polar) {

  if (polar < 0 || polar >= _num_polar)
    log_printf(ERROR, "Attempted to retrieve theta for polar angle = "
               "%d and azim angle = %d when only %d polar angles are "
               "defined", polar, azim, _num_polar);

  else if (azim < 0 || azim >= _num_azim)
    log_printf(ERROR, "Attempted to retrieve theta for polar angle = "
               "%d and azim angle = %d when only %d azim angles are "
               "defined", polar, azim, _num_azim);

  else if (_thetas == NULL)
    log_printf(ERROR, "Attempted to retrieve theta for polar angle = %d "
               "and azim angle = %d but the thetas have not been "
               "initialized", polar, azim);

  int first_polar = getFirstOctantPolar(polar);
  azim = getFirstOctantAzim(azim);

  if (polar < _num_polar/2)
    return _thetas[azim][first_polar];
  else
    return M_PI - _thetas[azim][first_polar];
}


/**
 * @brief Returns the azimuthal angle value.
 * @param azim index of the azimthal angle of interest
 * @return the value of the azimuthal angle
 */
double Quadrature::getPhi(int azim) {

  if (azim < 0 || azim >= _num_azim)
    log_printf(ERROR, "Attempted to retrieve theta for "
               "azim angle = %d when only %d azim angles are "
               "defined", azim, _num_azim);

  else if (_phis == NULL)
    log_printf(ERROR, "Attempted to retrieve phi for "
               "azim angle = %d but the phis have not been "
               "initialized", azim);

  int first_azim = getFirstOctantAzim(azim);

  if (azim < _num_azim/4)
    return _phis[first_azim];
  else if (azim < _num_azim/2)
    return M_PI - _phis[first_azim];
  else if (azim < 3*_num_azim/4)
    return M_PI + _phis[first_azim];
  else
    return 2.0 * M_PI - _phis[first_azim];
}


/**
 * @brief Returns the azimuthal angle weight value for a particular azimuthal angle.
 * @param azim index of the azimuthal angle of interest
 * @return the weight for an azimuthal angle
 */
FP_PRECISION Quadrature::getAzimWeight(int azim) {

  if (azim < 0 || azim >= _num_azim)
    log_printf(ERROR, "Attempted to retrieve the azimuthal weight for "
               "azimuthal angle = %d but only %d azimuthal angles "
               "are defined", azim, _num_azim);

  else if (_azim_weights == NULL)
    log_printf(ERROR, "Attempted to retrieve weight for azimuthal angle = %d "
               "but the azimuthal weights have not been initialized", azim);

  azim = getFirstOctantAzim(azim);
  return _azim_weights[azim];
}


/**
 * @brief Returns the polar weight for a particular azimuthal and polar angle.
 * @param azim index of the azimthal angle of interest
 * @param polar index of the polar angle of interest
 * @return the value of the polar weight for this azimuthal and polar angle
 */
FP_PRECISION Quadrature::getPolarWeight(int azim, int polar) {

  if (polar < 0 || polar >= _num_polar)
    log_printf(ERROR, "Attempted to retrieve polar weight for polar angle = "
               "%d and azim angle = %d when only %d polar angles are "
               "defined", polar, azim, _num_polar);

  else if (azim < 0 || azim >= _num_azim)
    log_printf(ERROR, "Attempted to retrieve polar weight for polar angle = "
               "%d and azim angle = %d when only %d azim angles are "
               "defined", polar, azim, _num_azim);

  else if (_polar_weights == NULL)
    log_printf(ERROR, "Attempted to retrieve polar weight for polar angle = %d "
               "and azim angle = %d but the thetas have not been "
               "initialized", polar, azim);

  azim = getFirstOctantAzim(azim);
  polar = getFirstOctantPolar(polar);

  return _polar_weights[azim][polar];
}


/**
 * @brief Returns the multiple value for a particular azimuthal and polar angle.
 * @details A multiple is the sine of a polar angle multiplied by its weight.
 * @param azim index of the azimuthal angle of interest
 * @param polar index of the polar angle of interest
 * @return the value of the sine of the polar angle multiplied with its weight
 */
FP_PRECISION Quadrature::getMultiple(int azim, int polar) {

  if (polar < 0 || polar >= _num_polar)
    log_printf(ERROR, "Attempted to retrieve the multiple for polar angle = "
               "%d and azimuthal angle = %d but only %d polar angles are "
               "defined", polar, azim, _num_polar);

  else if (azim < 0 || azim >= _num_azim)
    log_printf(ERROR, "Attempted to retrieve the multiple for polar angle = "
               "%d and azimuthal angle = %d but only %d azimuthal angles are "
               "defined", polar, azim, _num_azim);

  else if (_multiples == NULL)
    log_printf(ERROR, "Attempted to retrieve multiple for polar angle = %d "
               "and azimuthal angle = %d but the multiples have not been "
               "initialized", polar, azim);

  azim = getFirstOctantAzim(azim);
  polar = getFirstOctantPolar(polar);

  return _multiples[azim][polar];
}


/**
 * @brief Returns a pointer to the Quadrature's array of \f$ sin\theta_{p} \f$.
 * @return a pointer to the array of \f$ sin\theta_{p} \f$
 */
FP_PRECISION** Quadrature::getSinThetas() {

  if (_sin_thetas == NULL)
    log_printf(ERROR, "Attempted to retrieve the sin thetas array "
               "but it has not been initialized");

  return _sin_thetas;
}


/**
 * @brief Returns a pointer to the Quadrature's array of \f$ \theta_{p} \f$.
 * @return a pointer to the array of \f$ \theta_{p} \f$
 */
double** Quadrature::getThetas() {

  if (_thetas == NULL)
    log_printf(ERROR, "Attempted to retrieve the thetas array "
               "but it has not been initialized");

  return _thetas;
}


/**
 * @brief Returns a pointer to the Quadrature's array of \f$ \phi \f$.
 * @return a pointer to the array of \f$ \phi \f$
 */
double* Quadrature::getPhis() {

  if (_phis == NULL)
    log_printf(ERROR, "Attempted to retrieve the phis array "
               "but it has not been initialized");

  return _phis;
}


/**
 * @brief Returns a pointer to the Quadrature's array of azimuthal weights.
 * @return a pointer to the azimuthal weights array
 */
FP_PRECISION* Quadrature::getAzimWeights() {

  if (_azim_weights == NULL)
    log_printf(ERROR, "Attempted to retrieve the azimuthal weights array "
               "but it has not been initialized");

  return _azim_weights;
}


/**
 * @brief Returns a pointer to the Quadrature's array of polar weights.
 * @return a pointer to the polar weights array
 */
FP_PRECISION** Quadrature::getPolarWeights() {

  if (_polar_weights == NULL)
    log_printf(ERROR, "Attempted to retrieve the polar weights array "
               "but it has not been initialized");

  return _polar_weights;
}


/**
 * @brief Returns a pointer to the Quadrature's array of multiples.
 * @details A multiple is the sine of a polar angle multiplied by its weight.
 * @return a pointer to the multiples array
 */
FP_PRECISION** Quadrature::getMultiples() {

  if (_multiples == NULL)
    log_printf(ERROR, "Attempted to retrieve the multiples array "
               "but it has not been initialized");

  return _multiples;
}


int Quadrature::getFirstOctantAzim(int azim) {
  if (azim < _num_azim/4)
    return azim;
  else if (azim < _num_azim/2)
    return (_num_azim/2 - azim - 1);
  else if (azim < 3*_num_azim/4)
    return azim - _num_azim/2;
  else
    return (_num_azim - azim - 1);
}


int Quadrature::getFirstOctantPolar(int polar) {
  if (polar < _num_polar/2)
    return polar;
  else
    return (_num_polar - polar - 1);
}


int Quadrature::getOrthant(int azim, int polar) {

  int orthant;

  if (azim < _num_azim/4) {
    if (polar < _num_polar/2)
      orthant = 0;
    else
      orthant = 2;
  }
  else if (azim < _num_azim/2) {
    if (polar < _num_polar/2)
      orthant = 4;
    else
      orthant = 6;
  }
  else if (azim < 3*_num_azim/4) {
    if (polar < _num_polar/2)
      orthant = 3;
    else
      orthant = 1;
  }
  else{
    if (polar < _num_polar/2)
      orthant = 7;
    else
      orthant = 5;
  }

  return orthant;
}


/**
 * @brief Set the number of azimuthal angles to initialize.
 * @param num_azim the number of azimuthal angles
 */
void Quadrature::setNumAzimAngles(const int num_azim) {

  if (num_azim <= 0)
    log_printf(ERROR, "Unable to set the number of azimuthal angles to %d "
               "which is less than or equal to zero", num_azim);

  else if (num_azim % 4 != 0)
    log_printf(ERROR, "Unable to set the number of azimuthal angles to %d "
               "which is not divisible by 4", num_azim);

  _num_azim = num_azim;
}


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles
 */
void Quadrature::setNumPolarAngles(const int num_polar) {

  if (num_polar <= 0)
    log_printf(ERROR, "Unable to set the number of polar angles to %d "
               "which is less than or equal to zero", num_polar);

  else if (num_polar % 2 != 0)
    log_printf(ERROR, "Unable to set the number of polar angles to %d "
               "which is not divisible by 2", num_polar);

  _num_polar = num_polar;
}


/**
 * @brief Set the Quadrature's array of sines of each polar angle.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Quadrature's sin thetas in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 array the length
 *          of the number of polar angles) as input to this function. This
 *          function then fills the NumPy array with the data values for the
 *          Quadrature's sin thetas. An example of how this function might be
 *          called in Python is as follows:
 *
 * @code
 *          sin_thetas= numpy.array([0.05, 0.1, 0.15, ... ])
 *          polar_quad = openmoc.Quadrature()
 *          polar_quad.setNumPolarAngles(len(sin_thetas))
 *          polar_quad.setSinThetas(sin_thetas)
 * @endcode
 *
 * @param sin_thetas the array of sines of each polar angle
 * @param num_polar the number of polar angles
 */
void Quadrature::setThetas(double* thetas, int num_azim_times_polar) {

  if (_num_polar/2 * _num_azim/4 != num_azim_times_polar)
    log_printf(ERROR, "Unable to set %d thetas for Quadrature "
               "with %d polar angles and %d azimuthal angles"
               " in each octant",
               num_azim_times_polar, _num_polar/2, _num_azim/4);

  /* Deallocate memory if it was allocated previously */
  if (_thetas != NULL) {
    for (int i=0; i < _num_azim/4; i++)
      delete [] _thetas[i];
    delete [] _thetas;
  }

  /* Initialize memory for arrays */
  _thetas = new double*[_num_azim/4];
  for (int i=0; i < _num_azim/4; i++)
    _thetas[i] = new double[_num_polar/2];

  /* Extract sin thetas from user input */
  int ap=0;
  for (int a=0; a < _num_azim/4; a++) {
    for (int p=0; p < _num_polar/2; p++) {
      if (thetas[ap] < 0. || thetas[ap] > M_PI_2)
        log_printf(ERROR, "Unable to set theta to %f which is "
                   "not in the range [0,PI/2]", thetas[ap]);

      _thetas[a][p] = thetas[ap];
      ap++;
    }
  }
}


void Quadrature::setPolarWeights(double* weights, int num_azim_times_polar) {

  if (_num_polar/2 * _num_azim/4 != num_azim_times_polar)
    log_printf(ERROR, "Unable to set %d polar weights for Quadrature "
               "with %d polar angles and %d azimuthal angles"
               " in each octant",
               num_azim_times_polar, _num_polar/2, _num_azim/4);

  /* Deallocate memory if it was allocated previously */
  if (_polar_weights != NULL) {
    for (int i=0; i < _num_azim/4; i++)
      delete [] _polar_weights[i];
    delete [] _polar_weights;
  }

  /* Initialize memory for arrays */
  _polar_weights = new FP_PRECISION*[_num_azim/4];
  for (int i=0; i < _num_azim/4; i++)
    _polar_weights[i] = new FP_PRECISION[_num_polar/2];

  /* Extract polar weights from user input */
  int ap=0;
  for (int a=0; a < _num_azim/4; a++) {
    for (int p=0; p < _num_polar/2; p++) {
      if (weights[ap] < 0. || weights[ap] > M_PI_2)
        log_printf(ERROR, "Unable to polar weight to %f which is "
                   "not in the range [0,PI/2]", weights[ap]);

      _polar_weights[a][p] = FP_PRECISION(weights[ap]);
      ap++;
    }
  }
}


void Quadrature::setTheta(double theta, int azim, int polar) {

  if (theta <= 0.0 || theta >= M_PI_2)
    log_printf(ERROR, "Unable to set theta for azim = %d and polar = %d "
               "to %f which is not in the range (0.0, PI/2)", azim, polar,
               theta);

  else if (azim >= _num_azim/4)
    log_printf(ERROR, "Unable to set theta for azim = %d and polar = %d "
               "since azim is not in the range (0, _num_azim/4)", azim, polar);

  else if (polar >= _num_polar/2)
    log_printf(ERROR, "Unable to set theta for azim = %d and polar = %d "
               "since polar is not in the range (0, _num_polar/2)",
               azim, polar);


  if (_thetas == NULL) {
    _thetas = new double*[_num_azim/4];
    for (int i=0; i < _num_azim/4; i++)
      _thetas[i] = new double[_num_polar/2];
  }

  _thetas[azim][polar] = theta;
}


void Quadrature::setPhi(double phi, int azim) {

  if (phi <= 0.0 || phi >= M_PI_2)
    log_printf(ERROR, "Unable to set phi for azim = %d to %f which is not "
               "in the range (0.0, PI/2)", azim, phi);

  else if (azim >= _num_azim/4)
    log_printf(ERROR, "Unable to set phi for azim = %d since azim is not in"
               " the range (0, _num_azim/4)", azim);

  else if (_phis == NULL)
    _phis = new double[_num_azim/4];

  _phis[azim] = phi;
}


void Quadrature::setAzimWeight(double weight, int azim) {

  if (weight <= 0.0 || weight >= M_PI_2)
    log_printf(ERROR, "Unable to set azim weight for azim = %d to %f which is "
               "not in the range (0.0, PI/2)", azim, weight);

  else if (azim >= _num_azim/4)
    log_printf(ERROR, "Unable to set azim weight for azim = %d since azim is "
               "not in the range (0, _num_azim/4)", azim);

  else if (_azim_weights == NULL)
    _azim_weights = new FP_PRECISION[_num_azim/4];

  _azim_weights[azim] = FP_PRECISION(weight);
}


void Quadrature::setPolarWeight(double weight, int azim, int polar) {

  if (weight <= 0.0 || weight >= M_PI_2)
    log_printf(ERROR, "Unable to set polar weight for azim = %d and "
               "polar = %d to %f which is not in the range (0.0, PI/2)",
               azim, polar, weight);

  else if (azim >= _num_azim/4)
    log_printf(ERROR, "Unable to set polar weight for azim = %d and polar = %d "
               "since azim is not in the range (0, _num_azim/4)", azim, polar);

  else if (polar >= _num_polar/2)
    log_printf(ERROR, "Unable to set polar weight for azim = %d and polar = %d "
               "since polar is not in the range (0, _num_polar/2)", \
               azim, polar);

  else if (_polar_weights == NULL) {
    _polar_weights = new FP_PRECISION*[_num_azim/4];
    for (int i=0; i < _num_azim/4; i++)
      _polar_weights[i] = new FP_PRECISION[_num_polar/2];
  }

  _polar_weights[azim][polar] = FP_PRECISION(weight);
}


/**
 * @brief Dummy routine to initialize the polar quadrature.
 * @details The parent class routine simply checks that the number of polar
 *          angles has been set by the user and returns;
 */
void Quadrature::initialize() {

  if (_num_polar == 0)
    log_printf(ERROR, "Unable to initialize Quadrature with zero polar angles. "
               "Set the number of polar angles before initialization.");

  else if (_num_azim == 0)
    log_printf(ERROR, "Unable to initialize Quadrature with zero azimuthal angles. "
               "Set the number of azimuthal angles before initialization.");

  if (_phis == NULL)
    _phis = new double[_num_azim/4];

  /* Compute a desired azimuthal angles */
  for (int a = 0; a < _num_azim/4; a++)
    _phis[a] = 2.0 * M_PI / _num_azim * (0.5 + a);
}


/**
 * @brief This private routine computes the produce of the sine thetas and
 *        weights for each angle in the polar quadrature.
 * @details Note that this routine must be called after populating the
 *          sine thetas and weights arrays.
 */
void Quadrature::precomputeWeights(bool solve_3D) {

  double x1, x2;

  if (_azim_weights != NULL)
    delete [] _azim_weights;

  _azim_weights = new FP_PRECISION[_num_azim/4];

  /* Compute the azimuthal weights */
  for (int a = 0; a < _num_azim/4; a++) {

    /* The azimuthal weights (in radians) using equal weight quadrature */
    if (a < _num_azim/4 - 1)
      x1 = 0.5 * (_phis[a+1] - _phis[a]);
    else
      x1 = M_PI_2 - _phis[a];

    if (a >= 1)
      x2 = 0.5 * (_phis[a] - _phis[a-1]);
    else
      x2 = _phis[a];

    _azim_weights[a] = FP_PRECISION((x1 + x2) / M_PI);
  }

  /* Deallocate memory if it was allocated previously */
  if (_sin_thetas != NULL) {
    for (int a=0; a < _num_azim/4; a++)
      delete [] _sin_thetas[a];
    delete [] _sin_thetas;
  }

  /* Deallocate memory if it was allocated previously */
  if (_multiples != NULL) {
    for (int a=0; a < _num_azim/4; a++)
      delete [] _multiples[a];
    delete [] _multiples;
  }

  /* Initialize memory for arrays */
  _multiples = new FP_PRECISION*[_num_azim/4];
  _sin_thetas = new FP_PRECISION*[_num_azim/4];
  for (int a=0; a < _num_azim/4; a++) {
    _multiples[a] = new FP_PRECISION[_num_polar/2];
    _sin_thetas[a] = new FP_PRECISION[_num_polar/2];
  }


  /* Compute multiples of sine thetas and weights */
  for (int a=0; a < _num_azim/4; a++) {
    for (int p=0; p < _num_polar/2; p++) {
      _sin_thetas[a][p] = sin(_thetas[a][p]);
      _multiples[a][p] = _azim_weights[a] * _polar_weights[a][p];
      if (!solve_3D)
        _multiples[a][p] *= _sin_thetas[a][p];
    }
  }
}


/**
 * @brief Converts this Quadrature to a character array of its attributes.
 * @details The character array includes the number of polar angles, the
 *          the values of the sine and weight of each polar angle, and the
 *          product of the sine and weight of each polar angle.
 * @return a character array of the Quadrature's attributes
 */
std::string Quadrature::toString() {

  std::stringstream string;

  string << "Quadrature";
  string << "\n\t# azim angles  = " << _num_azim;
  string << "\n\t# polar angles = " << _num_polar;

  string << "\n\tphis = ";
  if (_phis != NULL) {
    for (int a = 0; a < _num_azim/4; a++)
      string << _phis[a] << ", ";
  }

  string << "\n\tazim weights = ";
  if (_azim_weights != NULL) {
    for (int a = 0; a < _num_azim/4; a++)
      string << _azim_weights[a] << ", ";
  }

  string << "\n\tthetas = ";
  if (_thetas != NULL) {
    for (int a = 0; a < _num_azim/4; a++) {
      for (int p = 0; p < _num_polar/2; p)
        string << " (" << a << "," << p << "): " << _thetas[a][p] << ", ";
    }
  }

  string << "\n\tpolar weights = ";
  if (_polar_weights != NULL) {
    for (int a = 0; a < _num_azim/4; a++) {
      for (int p = 0; p < _num_polar/2; p)
        string << " (" << a << "," << p << "): " << _polar_weights[a][p] << ", ";
    }
  }

  string << "\n\tsin thetas = ";
  if (_sin_thetas != NULL) {
    for (int a = 0; a < _num_azim/4; a++) {
      for (int p = 0; p < _num_polar/2; p)
        string << " (" << a << "," << p << "): " << _sin_thetas[a][p] << ", ";
    }
  }

  string << "\n\tmultiples = ";
  if (_multiples != NULL) {
    for (int a = 0; a < _num_azim/4; a++) {
      for (int p = 0; p < _num_polar/2; p)
        string << " (" << a << "," << p << "): " << _multiples[a][p] << ", ";
    }
  }

  return string.str();
}


quadratureType Quadrature::getQuadratureType() {
  return _quad_type;
}


/**
 * @brief Dummy constructor calls the parent constructor.
 */
TYPolarQuad::TYPolarQuad(): Quadrature() {
  _quad_type = TABUCHI_YAMAMOTO;
}


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles (maximum 3)
 */
void TYPolarQuad::setNumPolarAngles(const int num_polar) {

  if (num_polar > 6)
    log_printf(ERROR, "Unable to set the number of polar angles to %d "
               "for TYPolarQuad (max 6 angles)", num_polar);

  Quadrature::setNumPolarAngles(num_polar);
}


/**
 * @brief Routine to initialize the polar quadrature.
 * @details This routine uses the tabulated values for the Tabuchi-Yamamoto
 *          polar angle quadrature, including the sine thetas and weights.
 */
void TYPolarQuad::initialize() {

  /* Call parent class initialize routine */
  Quadrature::initialize();

  /* Allocate temporary arrays for tabulated quadrature values */
  double* thetas = new double[_num_polar/2*_num_azim/4];

  /* Tabulated values for the sine thetas and weights for the
   * Tabuchi-Yamamoto polar angle quadrature */
  if (_num_polar == 2) {
    for (int a=0; a < _num_azim/4; a++) {
      thetas[a] = asin(0.798184);
    }
  }

  else if (_num_polar == 4) {
    for (int a=0; a < _num_azim/4; a++) {
      thetas[a*(_num_polar/2)] = asin(0.363900);
      thetas[a*(_num_polar/2)+1] = asin(0.899900);
    }
  }

  else if (_num_polar == 6) {
    for (int a=0; a < _num_azim/4; a++) {
      thetas[a*(_num_polar/2)] = asin(0.166648);
      thetas[a*(_num_polar/2)+1] = asin(0.537707);
      thetas[a*(_num_polar/2)+2] = asin(0.932954);
    }
  }

  /* Set the arrays of thetas */
  Quadrature::setThetas(thetas, _num_polar/2*_num_azim/4);

  /* Deallocate temporary arrays */
  delete [] thetas;
}


void TYPolarQuad::precomputeWeights(bool solve_3D) {

  /* Allocate temporary arrays for tabulated quadrature values */
  double* weights = new double[_num_polar/2*_num_azim/4];

  /* Tabulated values for the sine thetas and weights for the
   * Tabuchi-Yamamoto polar angle quadrature */
  if (_num_polar == 2) {
    for (int a=0; a < _num_azim/4; a++) {
      weights[a] = 0.5;
    }
  }

  else if (_num_polar == 4) {
    for (int a=0; a < _num_azim/4; a++) {
      weights[a*(_num_polar/2)] = 0.212854 / 2.0;
      weights[a*(_num_polar/2)+1] = 0.787146 / 2.0;
    }
  }

  else if (_num_polar == 6) {
    for (int a=0; a < _num_azim/4; a++) {
      weights[a*(_num_polar/2)] = 0.046233 / 2.0;
      weights[a*(_num_polar/2)+1] = 0.283619 / 2.0;
      weights[a*(_num_polar/2)+2] = 0.670148 / 2.0;
    }
  }

  /* Set the arrays of sin thetas and weights */
  Quadrature::setPolarWeights(weights, _num_polar/2*_num_azim/4);

  /* Deallocate temporary arrays */
  delete [] weights;

  Quadrature::precomputeWeights(solve_3D);
}



/**
 * @brief Dummy constructor calls the parent constructor.
 */
LeonardPolarQuad::LeonardPolarQuad(): Quadrature() {
  _quad_type = LEONARD;
}


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles (2 or 3)
 */
void LeonardPolarQuad::setNumPolarAngles(const int num_polar) {

  if (num_polar != 4 && num_polar != 6)
    log_printf(ERROR, "Unable to set the number of polar angles to %d "
               "for LeonardPolarQuad (4 or 6 angles)", num_polar);

  Quadrature::setNumPolarAngles(num_polar);
}


/**
 * @brief Routine to initialize the polar quadrature.
 * @details This routine uses the tabulated values for the Leonard
 *          polar angle quadrature, including the sine thetas and weights.
 */
void LeonardPolarQuad::initialize() {

  /* Call parent class initialize routine */
  Quadrature::initialize();

  /* Allocate temporary arrays for tabulated quadrature values */
  double* thetas = new double[_num_polar/2*_num_azim/4];

  /* Tabulated values for the sine thetas and weights for the
   * Leonard polar angle quadrature */
  if (_num_polar == 4) {
    for (int a=0; a < _num_azim/4; a++) {
      thetas[a*(_num_polar/2)] = asin(0.273658);
      thetas[a*(_num_polar/2)+1] = asin(0.865714);
    }
  }

  else if (_num_polar == 6) {
    for (int a=0; a < _num_azim/4; a++) {
      thetas[a*(_num_polar/2)] = asin(0.099812);
      thetas[a*(_num_polar/2)+1] = asin(0.395534);
      thetas[a*(_num_polar/2)+2] = asin(0.891439);
    }
  }

  /* Set the arrays of thetas and weights */
  Quadrature::setThetas(thetas, _num_polar/2*_num_azim/4);

  /* Deallocate temporary arrays */
  delete [] thetas;
}


void LeonardPolarQuad::precomputeWeights(bool solve_3D) {

  /* Allocate temporary arrays for tabulated quadrature values */
  double* weights = new double[_num_polar/2*_num_azim/4];

  /* Tabulated values for the sine thetas and weights for the
   * Leonard polar angle quadrature */
  if (_num_polar == 4) {
    for (int a=0; a < _num_azim/4; a++) {
      weights[a*(_num_polar/2)] = 0.139473 / 2.0;
      weights[a*(_num_polar/2)+1] = 0.860527 / 2.0;
    }
  }

  else if (_num_polar == 6) {
    for (int a=0; a < _num_azim/4; a++) {
      weights[a*(_num_polar/2)] = 0.017620 / 2.0;
      weights[a*(_num_polar/2)+1] = 0.188561 / 2.0;
      weights[a*(_num_polar/2)+2] = 0.793819 / 2.0;
    }
  }

  /* Set the arrays of thetas and weights */
  Quadrature::setPolarWeights(weights, _num_polar/2*_num_azim/4);

  /* Deallocate temporary arrays */
  delete [] weights;

  Quadrature::precomputeWeights(solve_3D);
}



/**
 * @brief Dummy constructor calls the parent constructor.
 */
GLPolarQuad::GLPolarQuad(): Quadrature() {
  _quad_type = GAUSS_LEGENDRE;
}


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles (maximum 6)
 */
void GLPolarQuad::setNumPolarAngles(const int num_polar) {

  Quadrature::setNumPolarAngles(num_polar);
}


/**
  * @brief sets _use_adjusted_weights to true
  */
void GLPolarQuad::useAdjustedWeights() {
  _use_adjusted_weights = true;
}


/**
 * @brief Routine to initialize the polar quadrature.
 * @details This routine uses the tabulated values for the Gauss-Legendre
 *          polar angle quadrature, including the sine thetas and weights.
 */
void GLPolarQuad::initialize() {

  /* Call parent class initialize routine */
  Quadrature::initialize();

  /* By default don't adjust weights */
  _use_adjusted_weights = false;

  /* Allocate temporary arrays for tabulated quadrature values */
  double* thetas = new double[_num_polar/2*_num_azim/4];

  /* Tabulated values for the sine thetas and weights for the
   * Leonard polar angle quadrature */
  if (_num_polar == 2) {
    for (int a=0; a < _num_azim/4; a++) {
      thetas[a*(_num_polar/2)] = acos(0.5773502691);
    }
  }
  else if (_num_polar == 4) {
    for (int a=0; a < _num_azim/4; a++) {
      thetas[a*(_num_polar/2)] = acos(0.3399810435);
      thetas[a*(_num_polar/2)+1] = acos(0.8611363115);
    }
  }
  else if (_num_polar == 6) {
    for (int a=0; a < _num_azim/4; a++) {
      thetas[a*(_num_polar/2)] = acos(0.2386191860);
      thetas[a*(_num_polar/2)+1] = acos(0.6612093864);
      thetas[a*(_num_polar/2)+2] = acos(0.9324695142);
    }
  }
  else if (_num_polar == 8) {
    for (int a=0; a < _num_azim/4; a++) {
      thetas[a*(_num_polar/2)] = acos(0.1834346424);
      thetas[a*(_num_polar/2)+1] = acos(0.5255324099);
      thetas[a*(_num_polar/2)+2] = acos(0.7966664774);
      thetas[a*(_num_polar/2)+3] = acos(0.9602898564);
    }
  }
  else if (_num_polar == 10) {
    for (int a=0; a < _num_azim/4; a++) {
      thetas[a*(_num_polar/2)] = acos(0.1488743387);
      thetas[a*(_num_polar/2)+1] = acos(0.4333953941);
      thetas[a*(_num_polar/2)+2] = acos(0.6794095682);
      thetas[a*(_num_polar/2)+3] = acos(0.8650633666);
      thetas[a*(_num_polar/2)+4] = acos(0.9739065285);
    }
  }
  else if (_num_polar == 12) {
    for (int a=0; a < _num_azim/4; a++) {
      thetas[a*(_num_polar/2)] = acos(0.1252334085);
      thetas[a*(_num_polar/2)+1] = acos(0.3678314989);
      thetas[a*(_num_polar/2)+2] = acos(0.5873179542);
      thetas[a*(_num_polar/2)+3] = acos(0.7699026741);
      thetas[a*(_num_polar/2)+4] = acos(0.9041172563);
      thetas[a*(_num_polar/2)+5] = acos(0.9815606342);
    }
  }

  /* Set the arrays of sin thetas and weights */
  Quadrature::setThetas(thetas, _num_polar/2*_num_azim/4);

  /* Deallocate temporary arrays */
  delete [] thetas;
}


/**
 * @brief Calculates total weights for every azimuthal/polar combination based
 *        on the Gauss-Legendre polar quadrature.
 * @param solve_3D Boolean indicating whether this is a 3D quadrature
 */
void GLPolarQuad::precomputeWeights(bool solve_3D) {

  /* get weights */
  std::vector <double> weights_vec = getGLWeights(_roots, _num_polar);

  /* Allocate temporary arrays for tabulated quadrature values */
  FP_PRECISION weights[_num_polar/2*_num_azim/4];

  /* Tabulated values for the sine thetas and weights for the
   * Leonard polar angle quadrature */
  for (int a=0; a < _num_azim/4; a++) {

    /* Calculate simple weights */
    std::vector <double> nodes;
    for (int p=0; p<_num_polar/2; ++p) {
      nodes.push_back(cos(_thetas[a][p]));
    }

    std::vector <double> simple_weights = Quadrature::getSimpleWeights(nodes);

    for (int i=0; i<_num_polar/2; i++) {

      if (_use_adjusted_weights)
        
        /* Set weights based on adjusted polar angles */
        weights[a*(_num_polar/2)+i] = simple_weights[i] / 2.0;
      
      else 
        
        /* Set weights based on actual GL roots */
        weights[a*(_num_polar/2)+i] = weights_vec[i] / 2.0;
    }
  }

  /* Set the arrays of sin thetas and weights */
  Quadrature::setPolarWeights(weights, _num_polar/2*_num_azim/4);
  Quadrature::precomputeWeights(solve_3D);
}


/**
 * @brief Dummy constructor calls the parent constructor.
 */
EqualWeightPolarQuad::EqualWeightPolarQuad(): Quadrature() {
  _quad_type = EQUAL_WEIGHT;
}


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles
 */
void EqualWeightPolarQuad::setNumPolarAngles(const int num_polar) {
  Quadrature::setNumPolarAngles(num_polar);
}


/**
 * @brief Routine to initialize the polar quadrature.
 * @details This routine generates the sine thetas and weights.
 */
void EqualWeightPolarQuad::initialize() {

  /* Call parent class initialize routine */
  Quadrature::initialize();

  /* Allocate temporary arrays for tabulated quadrature values */
  double* thetas = new double[_num_polar/2*_num_azim/4];

  double cos_theta_a, cos_theta_b;

  /* Generate the sin thetas and weights using equations 420-422 of the
   * DOE Nucl. Eng. Handbook "Lattice Physics Computations" */
  int ap = 0;
  for (int a=0; a < _num_azim/4; a++) {
    cos_theta_a = 1.;
    for (int p=0; p < _num_polar/2; p++) {
      cos_theta_b = cos_theta_a - (1. / (_num_polar/2));
      thetas[ap] = acos(0.5 * (cos_theta_a + cos_theta_b));
      cos_theta_a = cos_theta_b;
      ap++;
    }
  }

  /* Set the arrays of sin thetas and weights */
  Quadrature::setThetas(thetas, _num_polar/2*_num_azim/4);

  /* Deallocate temporary arrays */
  delete [] thetas;
}


void EqualWeightPolarQuad::precomputeWeights(bool solve_3D) {

  /* Allocate temporary arrays for tabulated quadrature values */
  double* weights = new double[_num_polar/2*_num_azim/4];

  double y1, y2;

  /* Generate the sin thetas and weights using equations 420-422 of the
   * DOE Nucl. Eng. Handbook "Lattice Physics Computations" */
  int ap = 0;
  for (int a=0; a < _num_azim/4; a++) {
    for (int p=0; p < _num_polar/2; p++) {

      if (p < _num_polar/2 - 1)
        y1 = 0.5 * (cos(_thetas[a][p]) - cos(_thetas[a][p+1]));
      else
        y1 = cos(_thetas[a][p]);

      if (p >= 1)
        y2 = 0.5 * (cos(_thetas[a][p-1]) - cos(_thetas[a][p]));
      else
        y2 = 1.0 - cos(_thetas[a][p]);

      weights[ap] = (y1 + y2) / 2.0;
      ap++;
    }
  }

  /* Set the arrays of sin thetas and weights */
  Quadrature::setPolarWeights(weights, _num_polar/2*_num_azim/4);

  /* Deallocate temporary arrays */
  delete [] weights;

  /* Compute the product of the sine thetas and weights */
  Quadrature::precomputeWeights(solve_3D);
}



/**
 * @brief Dummy constructor calls the parent constructor.
 */
EqualAnglePolarQuad::EqualAnglePolarQuad(): Quadrature() {
  _quad_type = EQUAL_ANGLE;
}


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles
 */
void EqualAnglePolarQuad::setNumPolarAngles(const int num_polar) {
  Quadrature::setNumPolarAngles(num_polar);
}


/**
 * @brief Routine to initialize the polar quadrature.
 * @details This routine generates the sine thetas and weights.
 */
void EqualAnglePolarQuad::initialize() {

  /* Call parent class initialize routine */
  Quadrature::initialize();

  /* Allocate temporary arrays for tabulated quadrature values */
  double* thetas = new double[_num_polar/2*_num_azim/4];

  double cos_theta_a, cos_theta_b;
  double theta_a, theta_b;
  double delta_theta = M_PI / (_num_polar);


  /* Generate the sin thetas and weights using equations 420-422 of the
   * DOE Nucl. Eng. Handbook "Lattice Physics Computations" */
  int ap = 0;
  for (int a=0; a < _num_azim/4; a++) {
    theta_a = 0.;
    for (int p=0; p < _num_polar/2; p++) {
      theta_b = theta_a + delta_theta;
      cos_theta_a = cos(theta_a);
      cos_theta_b = cos(theta_b);
      thetas[ap] = acos((0.5 * (cos_theta_a + cos_theta_b)));
      theta_a = theta_b;
      ap++;
    }
  }

  /* Set the arrays of sin thetas and weights */
  Quadrature::setThetas(thetas, _num_polar/2*_num_azim/4);

  /* Deallocate temporary arrays */
  delete [] thetas;
}


void EqualAnglePolarQuad::precomputeWeights(bool solve_3D) {

  /* Allocate temporary arrays for tabulated quadrature values */
  double* weights = new double[_num_polar/2*_num_azim/4];

  double y1, y2;

  /* Generate the sin thetas and weights using equations 420-422 of the
   * DOE Nucl. Eng. Handbook "Lattice Physics Computations" */
  int ap = 0;
  for (int a=0; a < _num_azim/4; a++) {
    for (int p=0; p < _num_polar/2; p++) {

      if (p < _num_polar/2 - 1)
        y1 = 0.5 * (cos(_thetas[a][p]) - cos(_thetas[a][p+1]));
      else
        y1 = cos(_thetas[a][p]);

      if (p >= 1)
        y2 = 0.5 * (cos(_thetas[a][p-1]) - cos(_thetas[a][p]));
      else
        y2 = 1.0 - cos(_thetas[a][p]);

      weights[ap] = (y1 + y2) / 2.0;
      ap++;
    }
  }

  /* Set the arrays of sin thetas and weights */
  Quadrature::setPolarWeights(weights, _num_polar/2*_num_azim/4);

  /* Deallocate temporary arrays */
  delete [] weights;

  /* Compute the product of the sine thetas and weights */
  Quadrature::precomputeWeights(solve_3D);
}


/** 
 * @brief    the Legendre polynomial of degree n evaluated at x
 * @param    n an integer >=0: the order of the polynomial
 * @param    x in (-1,1), the point at which to evaluate the polynomial
 * @return   the value of the Legendre polynomial of degree n at x
 */
double Quadrature::legendrePolynomial(int n, double x) {
  if (n == 0)
    return 1;
  if (n == 1)
    return x;
  else {
    double c = 2.0*(n-2) + 3.0;
    double b = 1.0*(n-2) + 2.0;
    double a = 1.0*(n-2) + 1.0;
    double value =
      c/b * x * Quadrature::legendrePolynomial(n-1, x)
      - a/b * Quadrature::legendrePolynomial(n-2,x);
    return value;
  }
}
 

/**
 * @brief    the first logarithmic derivative of a Legendre polynomial
 * @param    m the order of the polynomial
 * @param    x point at which to evaluate the logarithmic derivative
 * @return   the value of the logarithmic derivative at x
 */
double GLPolarQuad::logDerivLegendre(int n, double x) {
  double num = n * x - n * Quadrature::legendrePolynomial(n-1,x)
    / Quadrature::legendrePolynomial(n,x);
  double denom = x*x - 1;
  return num/denom;
}


/**
 * @brief    the second logarithmic derivative of a Legendre polynomial
 * @param    m the order of the polynomial
 * @param    x point at which to evaluate the logarithmic derivative
 * @return   the value of the logarithmic derivative at x
 */
double GLPolarQuad::secondLogDerivLegendre(int n, double x) {
  double num =
    n*(n+1) + logDerivLegendre(n,x) * ((1-x*x)* logDerivLegendre(n,x) - 2 * x);
  double denom = x*x-1;
  return num/denom;
}
  

/**
 * @brief    finds the roots of Legendre polynomial of order n
 * @detail   guesses for positive roots are set at logarithmic intervals. 
 *           Positive roots are found simultaneously using an 
 *           Alberth-Householder-n method. Each guess is successively nudged
 *           towards a true root. Only the positive roots are calculated
 * @param    n the order of the polynomial
 * @return   a list of the roots of the polynomial
 */
std::vector <double> GLPolarQuad::getLegendreRoots(int n) {

  /* put these somewhere else */
  double E1 = 1e-10;
  double E2 = 1e-10;

  std::vector <double> roots;
  std::vector <bool> converged;
  std::vector <double> s1_tilde;
  std::vector <double> s2_tilde;

  /* set guesses with log scale*/
  for (int i=0; i < n/2; ++i) {
    roots.push_back(- pow(2, (-.5*(i+1))) +1);
    converged.push_back(false);
    s1_tilde.push_back(0);
    s2_tilde.push_back(0);

  }

  if (n%2 == 1) {
    roots.push_back(0);
    converged.push_back(true);
    s1_tilde.push_back(0);
    s2_tilde.push_back(0);
  }


  bool all_roots_converged = false;

  /* use the Alberth-Housholder_n method to nudge guesses towards roots */
  while (not all_roots_converged) {
    
    /* set S tildes */
    for (int i=0; i < (n+1)/2; ++i) {
      if (not converged[i]) {
        double sum1 = 0;
        double sum2 = 0;
        for (int j=0; j<= (n+1)/2; ++j) {
          if (j != i) {
            sum1 += 1/(roots[i] - roots[j]);
            sum2 += -1/((roots[i] - roots[j])*(roots[i] - roots[j]));
          }
        }
    
        s1_tilde[i] = logDerivLegendre(n, roots[i]) - sum1;
        s2_tilde[i] = secondLogDerivLegendre(n, roots[i]) - sum2;

        /* householder method 2  Halley     */
        double u_new =
          roots[i] - 2*s1_tilde[i] / (s1_tilde[i]*s1_tilde[i] - s2_tilde[i]);
        double u_old = roots[i];
        roots[i] = u_new;
       
        /* if this is the actual root */
        if (abs(u_new - u_old) < E1) {
          if (std::abs(Quadrature::legendrePolynomial(n, u_new)) < E2) {
            converged[i] = true;

            /* if this root equals another root or it is less than 0 */
            for (int j=0; j < (n+1)/2; ++j) { 
              if (j != i) {
                if (std::abs(roots[j] - roots[i]) < E1 or roots[i] <= 0) {

                  /* reset the root to its original guess */
                  roots[i] = - pow(2, (-.5*(i+1))) +1;
                  converged[i] = false;
                }
              }
            }
          }
        } /* if this is the actual root */
      } /* if not converged */
    } /* for each guess */

    for (int i=0; i<(n+1)/2; ++ i) {
      all_roots_converged = converged[i];
      if (not all_roots_converged)
        break;
    }
  } /* while not all roots converged */
 
  /* add negative roots */
  std::sort (roots.begin(), roots.end());
  return roots;
}


/**
 * @brief    calculates the weights to be used in Gauss-Legendre Quadrature
 * @param    roots a vector containing the roots of the Legendre polynomial
 * @param    n the order of the Legendre Polynomial
 * @return   a vector of weights matched by index to the vector of roots
 */
std::vector <double> GLPolarQuad::getGLWeights(std::vector <double> roots,
                                               int n) {
  std::vector <double> weights;
  for (int i; i<roots.size(); ++i){
    weights.push_back(
        - (2*roots[i]*roots[i] - 2)
        / (n*n*Quadrature::legendrePolynomial(n-1, roots[i])
          * Quadrature::legendrePolynomial(n-1, roots[i])));
  }

  return weights;
}


/**
 * @brief    calculates the weights to be used in Gauss-Legendre Quadrature
 * @param    root the root of the Legendre polynomial
 * @param    n the order of the Legendre Polynomial
 * @return   a vector of weights matched by index to the vector of roots
 */
double Quadrature::getSingleGLWeight(double root, int n) {
  double weight =
        - (2*root*root - 2) / (n*n*Quadrature::legendrePolynomial(n-1, root)
                                  * Quadrature::legendrePolynomial(n-1, root));

  return weight;
}


/**
 * @brief    calculates the weights to be used in numerical integration
 * @details  assumes the function will be integrated over (-1, 1)
 * @param    nodes a vector containing the x's that are evaluated
 * @return   a vector of weights matched by index to the vector of nodes
 */
std::vector <double> Quadrature::getSimpleWeights(std::vector <double> nodes) {

  int n = nodes.size();
  std::vector <double> weights;
  std::cout.precision(14);
  std::cout.setf(std::ios::fixed);
  
  //TEMPORARY
  for (int i=0; i<n; ++i)
    nodes.push_back(-nodes[i]);
  double a = -1;
  n = nodes.size();

  // declare an array to store the elements of the augmented-matrix
  double A[n][n+1];

  // the solution array
  double x[n];

  double b = 1;

  // populate A
  for (int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j) {
      A[i][j] = pow(nodes[j], i);
    }
    A[i][n] = (pow(b, i+1) - pow(a, i+1)) / (i+1);
  }

  // pivotisation 
  for (int i=0; i<n; i++) {
    for (int k=i+1; k<n; k++) {
      if (A[i][i] < A[k][i]) {
        for (int j=0; j<=n; j++) {
          double temp = A[i][j];
          A[i][j] = A[k][j];
          A[k][j] = temp;
        }
      }
    }
  }

  // perform gauss elimination
  for (int i=0; i<n-1; i++) {
    for (int k=i+1; k<n; k++) {
      double t = A[k][i] / A[i][i];

      // make elements below the pivot elements equal to zero or eliminate the
      // variables
      for (int j=0; j<=n; j++)
        A[k][j] = A[k][j] - t * A[i][j];
     }
  }

  // back-substitution
  for (int i=n-1; i>=0; --i) {
    double sub = 0;
    for (int j=n-1; j>i; --j) {
      sub += x[j]*A[i][j];
    }
    x[i] = (A[i][n] - sub) / A[i][i];
  }

  double sum = 0;
  for (int i=0; i<n; ++i){
    weights.push_back(x[i]);
    sum += x[i];
  }

  return weights;

}
