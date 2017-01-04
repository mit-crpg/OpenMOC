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
  _azim_spacings = NULL;
  _azim_weights = NULL;
  _polar_weights = NULL;
  _total_weights = NULL;
}


/**
 * @brief Destructor deletes arrray of sines of the polar angles, the weights
 *        of the polar angles and the products of the sines and weights.
 */
Quadrature::~Quadrature() {
  deleteAllArrays();
}


/**
 * @brief Deletes all arrays indexed by polar angle
 */
void Quadrature::deletePolarArrays() {

  if (_sin_thetas != NULL) {
    for (int a=0; a < _num_azim/2; a++)
      delete [] _sin_thetas[a];
    delete [] _sin_thetas;
    _sin_thetas = NULL;
  }

  if (_thetas != NULL) {
    for (int a=0; a < _num_azim/2; a++)
      delete [] _thetas[a];
    delete [] _thetas;
    _thetas = NULL;
  }

  if (_polar_weights != NULL) {
    for (int a=0; a < _num_azim/2; a++)
      delete [] _polar_weights[a];
    delete [] _polar_weights;
    _polar_weights = NULL;
  }

  if (_total_weights != NULL) {
    for (int a=0; a < _num_azim/2; a++)
      delete [] _total_weights[a];
    delete [] _total_weights;
    _total_weights = NULL;
  }
}


/**
 * @brief Deletes all arrays allocated by the Quadrature
 */
void Quadrature::deleteAllArrays() {

  if (_phis != NULL)
    delete [] _phis;
  _phis = NULL;

  if (_azim_spacings != NULL)
    delete [] _azim_spacings;
  _azim_spacings = NULL;

  if (_azim_weights != NULL)
    delete [] _azim_weights;
  _azim_weights = NULL;

  deletePolarArrays();
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

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

  return _sin_thetas[azim][polar];
}


/**
 * @brief Returns the polar angle in radians for a given azimuthal and polar
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

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

  return _thetas[azim][polar];
}


/**
 * @brief Returns the azimuthal angle value in radians.
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

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

  return _phis[azim];
}


/**
 * @brief Returns the azimuthal angle weight value for a particular azimuthal
 *        angle.
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

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

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

  if (azim < 0 || azim >= _num_azim)
    log_printf(ERROR, "Attempted to retrieve polar weight for polar angle = "
               "%d and azim angle = %d when only %d azim angles are "
               "defined", polar, azim, _num_azim);

  if (_polar_weights == NULL)
    log_printf(ERROR, "Attempted to retrieve polar weight for polar angle = %d "
               "and azim angle = %d but the thetas have not been "
               "initialized", polar, azim);

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

  return _polar_weights[azim][polar];
}


/**
 * @brief Returns the total weight for Tracks with the given azimuthal and
 *        polar indexes
 * @details Angular weights are multiplied by Track spcings
 * @param azim index of the azimuthal angle of interest
 * @param polar index of the polar angle of interest
 * @return the total weight of each Track with the given indexes
 */
FP_PRECISION Quadrature::getWeight(int azim, int polar) {

  if (polar < 0 || polar >= _num_polar)
    log_printf(ERROR, "Attempted to retrieve the weight for polar angle = "
               "%d and azimuthal angle = %d but only %d polar angles are "
               "defined", polar, azim, _num_polar);

  else if (azim < 0 || azim >= _num_azim)
    log_printf(ERROR, "Attempted to retrieve the weight for polar angle = "
               "%d and azimuthal angle = %d but only %d azimuthal angles are "
               "defined", polar, azim, _num_azim);

  else if (_total_weights == NULL)
    log_printf(ERROR, "Attempted to retrieve weight for polar angle = %d "
               "and azimuthal angle = %d but the multiples have not been "
               "initialized", polar, azim);

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

  return _total_weights[azim][polar];
}


/**
 * @brief Returns a pointer to the Quadrature's array of polar angle sines
          \f$ sin\theta_{p} \f$.
 * @return a pointer to the array of \f$ sin\theta_{p} \f$
 */
FP_PRECISION** Quadrature::getSinThetas() {

  if (_sin_thetas == NULL)
    log_printf(ERROR, "Attempted to retrieve the sin thetas array "
               "but it has not been initialized");

  return _sin_thetas;
}


/**
 * @brief Returns a pointer to the Quadrature's array of polar angles
          \f$ \theta_{p} \f$.
 * @return a pointer to the array of \f$ \theta_{p} \f$
 */
double** Quadrature::getThetas() {

  if (_thetas == NULL)
    log_printf(ERROR, "Attempted to retrieve the thetas array "
               "but it has not been initialized");

  return _thetas;
}


/**
 * @brief Returns a pointer to the Quadrature's array of azimuthal angles
          \f$ \phi \f$.
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

  if (num_azim != _num_azim) {

    /* Delete arrays with old settings */
    deleteAllArrays();
    _num_azim = num_azim;
  }
}


/**
 * @brief Returns an array of adjusted azimuthal spacings.
 * @details An array of azimuthal spacings after adjustment is returned,
 *          indexed by azimuthal angle
 * @return the array of azimuthal spacings
 */
FP_PRECISION* Quadrature::getAzimSpacings() {
  return _azim_spacings;
}


/**
 * @brief Returns the adjusted azimuthal spacing at the requested azimuthal
 *        angle index.
 * @details The aziumthal spacing depends on the azimuthal angle. This function
 *          returns the azimuthal spacing used at the desired azimuthal angle
 *          index.
 * @param azim the requested azimuthal angle index
 * @return the requested azimuthal spacing
 */
FP_PRECISION Quadrature::getAzimSpacing(int azim) {
  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;
  return _azim_spacings[azim];
}


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles
 */
void Quadrature::setNumPolarAngles(const int num_polar) {

  if (num_polar <= 0)
    log_printf(ERROR, "Unable to set the number of polar angles to %d "
               "which is less than or equal to zero", num_polar);

  if (num_polar % 2 != 0)
    log_printf(ERROR, "Unable to set the number of polar angles to %d "
               "which is not divisible by 2", num_polar);

  if (num_polar != _num_polar) {

    /* Delete arrays with old settings */
    deletePolarArrays();
    _num_polar = num_polar;
  }
}


/**
 * @brief Sets the Quadrature's array of polar angles.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Quadrature's polar angles in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 array the length
 *          of the number of azimuthal times polar angles) as input to this
 *          function. This function then fills the Quadrature's polar angles
 *          with the given values. An example of how this function might be
 *          called in Python is as follows:
 *
 * @code
 *          thetas = numpy.array([pi/6, pi/4, pi/3, ... ])
 *          polar_quad = openmoc.Quadrature()
 *          polar_quad.setNumAzimAngles(num_azim)
 *          polar_quad.setNumPolarAngles(len(thetas) / num_azim)
 *          polar_quad.setThetas(thetas)
 * @endcode
 *
 * @param thetas the array of polar angle for each azimuthal/polar angle
 *        combination
 * @param num_azim_times_polar the total number of angles (azimuthal x polar)
 */
void Quadrature::setThetas(double* thetas, int num_azim_times_polar) {

  if (_num_polar/2 * _num_azim/4 != num_azim_times_polar)
    log_printf(ERROR, "Unable to set %d thetas for Quadrature "
               "with %d polar angles and %d azimuthal angles"
               " in each octant",
               num_azim_times_polar, _num_polar/2, _num_azim/4);

  /* Initialize memory for arrays */
  if (_thetas == NULL) {
    _thetas = new double*[_num_azim/2];
    for (int i=0; i < _num_azim/2; i++)
      _thetas[i] = new double[_num_polar];
  }

  /* Extract sin thetas from user input */
  int ap=0;
  for (int a=0; a < _num_azim/4; a++) {
    for (int p=0; p < _num_polar/2; p++) {
      if (thetas[ap] < 0. || thetas[ap] > M_PI_2)
        log_printf(ERROR, "Unable to set theta to %f which is "
                   "not in the range [0,PI/2]", thetas[ap]);

      _thetas[a][p] = thetas[ap];
      _thetas[_num_azim/2 - a - 1][p] = thetas[ap];
      _thetas[a][_num_polar - p - 1] = M_PI - thetas[ap];
      _thetas[_num_azim/2 - a - 1][_num_polar - p - 1] = M_PI - thetas[ap];
      ap++;
    }
  }
}


/**
 * @brief Set the Quadrature's array of polar weights.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Quadrature's polar weights in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 array the length
 *          of the number of azimuthal times polar angles) as input to this
 *          function. This function then fills the Quadrature's polar weights
 *          with the given values. An example of how this function might be
 *          called in Python is as follows:
 *
 * @code
 *          polar_weights = numpy.array([0.1, 0.2, 0.05, ... ])
 *          polar_quad = openmoc.Quadrature()
 *          polar_quad.setNumAzimAngles(num_azim)
 *          polar_quad.setNumPolarAngles(len(polar_weights) / num_azim)
 *          polar_quad.setPolarWeights(polar_weights)
 * @endcode
 *
 * @param weights The polar weights
 * @param num_azim_times_polar the total number of angles in one octant
 *        (azimuthal x polar)
 */
void Quadrature::setPolarWeights(FP_PRECISION* weights,
                                 int num_azim_times_polar) {

  if (_num_polar/2 * _num_azim/4 != num_azim_times_polar)
    log_printf(ERROR, "Unable to set %d polar weights for Quadrature "
               "with %d polar angles and %d azimuthal angles"
               " in each octant",
               num_azim_times_polar, _num_polar/2, _num_azim/4);

  /* Initialize memory for arrays */
  if (_polar_weights == NULL) {
    _polar_weights = new FP_PRECISION*[_num_azim/2];
    for (int i=0; i < _num_azim/2; i++)
      _polar_weights[i] = new FP_PRECISION[_num_polar];
  }

  /* Extract polar weights from user input */
  int ap=0;
  for (int a=0; a < _num_azim/4; a++) {
    for (int p=0; p < _num_polar/2; p++) {
      if (weights[ap] < 0. || weights[ap] > M_PI_2)
        log_printf(ERROR, "Unable to set polar weight to %f which is "
                   "not in the range [0,PI/2]", weights[ap]);

      setPolarValues(_polar_weights, a, p, weights[ap]);
      ap++;
    }
  }
}


/**
 * @brief Sets the polar angle for the given indexes.
 * @param theta the value in radians of the polar angle to be set
 * @param azim the azimuthal index of the angle of interest
 * @param polar the polar index of the angle of interest
 */
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
    _thetas = new double*[_num_azim/2];
    for (int i=0; i < _num_azim/2; i++)
      _thetas[i] = new double[_num_polar];
  }

  _thetas[azim][polar] = theta;
  _thetas[_num_azim/2 - azim - 1][polar] = theta;
  _thetas[azim][_num_polar - polar - 1] = M_PI - theta;
  _thetas[_num_azim/2 - azim - 1][_num_polar - polar - 1] = M_PI - theta;
}


/**
 * @brief Sets the azimuthal angle for the given index.
 * @param phi the value in radians of the azimuthal angle to be set
 * @param azim the azimuthal index
 */
void Quadrature::setPhi(double phi, int azim) {

  if (phi <= 0.0 || phi >= M_PI_2)
    log_printf(ERROR, "Unable to set phi for azim = %d to %f which is not "
               "in the range (0.0, PI/2)", azim, phi);

  if (azim >= _num_azim/4)
    log_printf(ERROR, "Unable to set phi for azim = %d since azim is not in"
               " the range (0, _num_azim/4)", azim);

  if (_phis == NULL)
    _phis = new double[_num_azim/2];

  _phis[azim] = phi;
  _phis[_num_azim/2 - azim - 1] = M_PI - phi;
}


/**
 * @brief Sets the azimuthal spacing for the given index.
 * @param spacing the spacing (cm) in the azimuthal direction to be set
 * @param azim the azimuthal index
 */
void Quadrature::setAzimSpacing(FP_PRECISION spacing, int azim) {

  if (spacing <= 0.0)
    log_printf(ERROR, "Unable to set azimuthal spacing for azim = %d to %f "
                      "which is not strictly greater than zero", azim,
                      spacing);

  if (azim >= _num_azim/4)
    log_printf(ERROR, "Unable to set azimuthal spacing for azim = %d since "
                      " azim is not in the range (0, _num_azim/4)", azim);

  if (_azim_spacings == NULL)
    _azim_spacings = new FP_PRECISION[_num_azim/2];

  setAzimuthalValues(_azim_spacings, azim, spacing);
}


/**
 * @brief Sets the azimuthal weight for the given index.
 * @param weight the weight of the azimuthal angle
 * @param azim the azimuthal index
 */
void Quadrature::setAzimWeight(double weight, int azim) {

  if (weight <= 0.0 || weight >= M_PI_2)
    log_printf(ERROR, "Unable to set azim weight for azim = %d to %f which is "
               "not in the range (0.0, PI/2)", azim, weight);

  if (azim >= _num_azim/4)
    log_printf(ERROR, "Unable to set azim weight for azim = %d since azim is "
               "not in the range (0, _num_azim/4)", azim);

  if (_azim_weights == NULL)
    _azim_weights = new FP_PRECISION[_num_azim/2];

  setAzimuthalValues(_azim_weights, azim, FP_PRECISION(weight));
}


/**
 * @brief Sets the polar weight for the given indexes.
 * @param weight the weight of the polar angle
 * @param azim the azimuthal index corresponding to the angle
 * @param azim the polar index corresponding to the angle
 */
void Quadrature::setPolarWeight(FP_PRECISION weight, int azim, int polar) {

  if (weight <= 0.0 || weight >= M_PI_2)
    log_printf(ERROR, "Unable to set polar weight for azim = %d and "
               "polar = %d to %f which is not in the range (0.0, PI/2)",
               azim, polar, weight);

  if (azim >= _num_azim/4)
    log_printf(ERROR, "Unable to set polar weight for azim = %d and polar = %d "
               "since azim is not in the range (0, _num_azim/4)", azim, polar);

  if (polar >= _num_polar/2)
    log_printf(ERROR, "Unable to set polar weight for azim = %d and polar = %d "
               "since polar is not in the range (0, _num_polar/2)", \
               azim, polar);

  if (_polar_weights == NULL) {
    _polar_weights = new FP_PRECISION*[_num_azim/2];
    for (int a=0; a < _num_azim/2; a++)
      _polar_weights[a] = new FP_PRECISION[_num_polar];
  }

  setPolarValues(_polar_weights, azim, polar, FP_PRECISION(weight));
}


/**
 * @brief Initialize the polar quadrature azimuthal angles.
 * @details The parent class routine simply checks that number of polar and
 *          azimuthal angles have been set by the user and generates the
 *          azimuthal angles if not already generated.
 */
void Quadrature::initialize() {

  if (_num_polar == 0)
    log_printf(ERROR, "Unable to initialize Quadrature with zero polar angles. "
               "Set the number of polar angles before initialization.");

  if (_num_azim == 0)
    log_printf(ERROR, "Unable to initialize Quadrature with zero azimuthal angles. "
               "Set the number of azimuthal angles before initialization.");

  if (_phis == NULL)
    _phis = new double[_num_azim/2];

  /* Compute a desired azimuthal angles */
  for (int a = 0; a < _num_azim/2; a++)
    _phis[a] = 2.0 * M_PI / _num_azim * (0.5 + a);
}


/**
 * @brief This private routine computes the product of the sine thetas and
 *        weights for each angle in the polar quadrature.
 * @details Note that this routine must be called after populating the
 *          sine thetas and weights arrays.
 */
void Quadrature::precomputeWeights(bool solve_3D) {

  /* Check that track spacings have been set */
  if (_azim_spacings == NULL)
    log_printf(ERROR, "Unable to precompute weights since track spacings have "
                      "not yet been set");

  /* Check that polar angles have been set */
  if (_thetas == NULL)
    log_printf(ERROR, "Unable to precompute weights since polar angles have "
                      "not yet been set");

  /* Clear azimuthal weights */
  if (_azim_weights == NULL)
    _azim_weights = new FP_PRECISION[_num_azim/2];

  /* Create uncorrected weights if no angles have been set yet */
  if (_phis == NULL) {
    log_printf(NORMAL, "WARNING: Using uncorrected angles for weights");
    double phi = M_PI / _num_azim;
    for (int a = 0; a < _num_azim/4; a++) {
      setPhi(phi, a);
      phi += 2*M_PI / _num_azim;
    }
  }

  /* Compute the azimuthal weights */
  double x1, x2;
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

    setAzimuthalValues(_azim_weights, a, FP_PRECISION((x1 + x2) / M_PI));
  }

  /* Allocate memory if it was not allocated previously */
  if (_sin_thetas == NULL) {
    _sin_thetas = new FP_PRECISION*[_num_azim/2];
    for (int a=0; a < _num_azim/2; a++)
      _sin_thetas[a] = new FP_PRECISION[_num_polar];
  }

  /* Allocate memory if it was not allocated previously */
  if (_total_weights == NULL) {
    _total_weights = new FP_PRECISION*[_num_azim/2];
    for (int a=0; a < _num_azim/2; a++)
    _total_weights[a] = new FP_PRECISION[_num_polar];
  }

  /* Compute multiples of sine thetas and weights */
  for (int a=0; a < _num_azim/4; a++) {
    for (int p=0; p < _num_polar/2; p++) {
      _sin_thetas[a][p] = sin(_thetas[a][p]);
      FP_PRECISION weight = 2.0 * M_PI * _azim_weights[a]
                            * _azim_spacings[a] * _polar_weights[a][p];
      weight *= 2.0 * _sin_thetas[a][p];
      setPolarValues(_total_weights, a, p, weight);
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

  string << "\n\ttotal weights = ";
  if (_total_weights != NULL) {
    for (int a = 0; a < _num_azim/4; a++) {
      for (int p = 0; p < _num_polar/2; p)
        string << " (" << a << "," << p << "): " << _total_weights[a][p] << ", ";
    }
  }

  return string.str();
}


/**
 * @breif Returns the type of Quadrature created.
 * @return The quadrature type
 */
quadratureType Quadrature::getQuadratureType() {
  return _quad_type;
}


/**
 * @brief Dummy constructor calls the parent constructor.
 */
TYPolarQuad::TYPolarQuad(): Quadrature() {
  _quad_type = TABUCHI_YAMAMOTO;
  _num_polar = 6;
}


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles (maximum 6)
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
  double thetas[_num_polar/2*_num_azim/4];

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
}


/**
 * @brief Calculates total weights for every azimuthal/polar combination based
 *        on the TY quadrature.
 * @param solve_3D Boolean indicating whether this is a 3D quadrature
 */
void TYPolarQuad::precomputeWeights(bool solve_3D) {

  /* Allocate temporary arrays for tabulated quadrature values */
  FP_PRECISION weights[_num_polar/2*_num_azim/4];

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
  Quadrature::precomputeWeights(solve_3D);
}



/**
 * @brief Dummy constructor calls the parent constructor.
 */
LeonardPolarQuad::LeonardPolarQuad(): Quadrature() {
  _quad_type = LEONARD;
  _num_polar = 6;
}


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles (4 or 6)
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
  double thetas[_num_polar/2*_num_azim/4];

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
}


/**
 * @brief Calculates total weights for every azimuthal/polar combination based
 *        on the Leonard polar quadrature.
 * @param solve_3D Boolean indicating whether this is a 3D quadrature
 */
void LeonardPolarQuad::precomputeWeights(bool solve_3D) {

  /* Allocate temporary arrays for tabulated quadrature values */
  FP_PRECISION weights[_num_polar/2*_num_azim/4];

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
  Quadrature::precomputeWeights(solve_3D);
}



/**
 * @brief Dummy constructor calls the parent constructor.
 */
GLPolarQuad::GLPolarQuad(): Quadrature() {
  _quad_type = GAUSS_LEGENDRE;
  _num_polar = 6;
}


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles (maximum 12)
 */
void GLPolarQuad::setNumPolarAngles(const int num_polar) {

  Quadrature::setNumPolarAngles(num_polar);
}


/**
 * @brief Routine to initialize the polar quadrature.
 * @details This routine uses the tabulated values for the Gauss-Legendre
 *          polar angle quadrature, including the sine thetas and weights.
 */
void GLPolarQuad::initialize() {
  
  /* Call parent class initialize routine */
  Quadrature::initialize();

  /* Allocate temporary arrays for tabulated quadrature values */
  double thetas[_num_polar/2*_num_azim/4];

  /* Get roots of Legendre polynomial */
  _roots = getLegendreRoots(_num_polar);

  /* Set theta values for polar angles */
  for (int a=0; a < _num_azim/4; a++) {
    for (int i=0; i < _num_polar/2; i++) {
      thetas[a*(_num_polar/2)+i] = acos(_roots[i]);
    }
  }

  /* Set the arrays of sin thetas and weights */
  Quadrature::setThetas(thetas, _num_polar/2*_num_azim/4);
}


/**
 * @brief Calculates total weights for every azimuthal/polar combination based
 *        on the Gauss-Legendre polar quadrature.
 * @param solve_3D Boolean indicating whether this is a 3D quadrature
 */
void GLPolarQuad::precomputeWeights(bool solve_3D) {

  /* Get weights */
  std::vector <double> weights_vec = getGLWeights(_roots, _num_polar);

  /* Allocate temporary arrays for tabulated quadrature values */
  FP_PRECISION weights[_num_polar/2*_num_azim/4];

  /* Weights for Gauss Legendre polar angle quadrature */
  for (int a=0; a < _num_azim/4; a++) {
    for (int i=0; i < _num_polar/2; i++) {
      weights[a*(_num_polar/2)+i] = weights_vec[i] / 2.0;
    }
  }

  /* Set the arrays of sin thetas and weights */
  Quadrature::setPolarWeights(weights, _num_polar/2*_num_azim/4);
  Quadrature::precomputeWeights(solve_3D);
}


/** 
 * @brief    the Legendre polynomial of degree n evaluated at x
 * @param    n an integer >=0: the order of the polynomial
 * @param    x in (-1,1), the point at which to evaluate the polynomial
 * @return   the value of the Legendre polynomial of degree n at x
 */
double GLPolarQuad::legendrePolynomial(int n, double x) {
  if (n == 0)
    return 1;
  if (n == 1)
    return x;
  else {
    double c = 2.0*(n-2) + 3.0;
    double b = 1.0*(n-2) + 2.0;
    double a = 1.0*(n-2) + 1.0;
    double value =
      c/b * x * legendrePolynomial(n-1, x) - a/b * legendrePolynomial(n-2,x);
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
  double num = n * x - n * legendrePolynomial(n-1,x) / legendrePolynomial(n,x);
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

  double E1 = 1e-10;
  double E2 = 1e-10;

  std::vector <double> roots;
  std::vector <bool> converged;
  std::vector <double> s1_tilde;
  std::vector <double> s2_tilde;

  /* Set guesses with log scale */
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

  /* Placeholder element */
  roots.push_back(0);

  bool all_roots_converged = false;

  /* Use the Alberth-Housholder_n method to nudge guesses towards roots */
  while (not all_roots_converged) {
    
    /* Set S tildes */
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

        /* Householder method 2 Halley */
        double u_new =
          roots[i] - 2*s1_tilde[i] / (s1_tilde[i]*s1_tilde[i] - s2_tilde[i]);
        double u_old = roots[i];
        roots[i] = u_new;
        
        /* If this is the actual root */
        if (abs(u_new - u_old) < E1) {
          if (std::abs(legendrePolynomial(n, u_new)) < E2) {
            converged[i] = true;

            /* If this root equals another root or it is less than 0 */
            for (int j=0; j < (n+1)/2; ++j) { 
              if (j != i) {
                if (std::abs(roots[j] - roots[i]) < E1 or roots[i] <= 0) {

                  /* Reset the root to its original guess */
                  roots[i] = - pow(2, (-.5*(i+1))) +1;
                  converged[i] = false;
                }
              }
            }
          }
        } /* If this is the actual root */
      } /* If not converged */
    } /* For each guess */

    for (int i=0; i<(n+1)/2; ++ i) {
      all_roots_converged = converged[i];
      if (not all_roots_converged)
        break;
    }
  } /* While not all roots converged */
  
  /* Remove placeholder root */
  roots.pop_back();

  /* Put roots in order */
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
  for (int i; i<roots.size(); ++i)
    weights.push_back(
        - (2*roots[i]*roots[i] - 2) / 
        (n*n*legendrePolynomial(n-1, roots[i])
         *legendrePolynomial(n-1, roots[i])));
  
  return weights;
}


/**
 * @brief Dummy constructor calls the parent constructor.
 */
EqualWeightPolarQuad::EqualWeightPolarQuad(): Quadrature() {
  _quad_type = EQUAL_WEIGHT;
  _num_polar = 6;
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
  double thetas[_num_polar/2*_num_azim/4];

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
}


/**
 * @brief Calculates total weights for every azimuthal/polar combination based
 *        on the equal weight polar quadrature.
 * @param solve_3D Boolean indicating whether this is a 3D quadrature
 */
void EqualWeightPolarQuad::precomputeWeights(bool solve_3D) {

  /* Allocate temporary arrays for tabulated quadrature values */
  FP_PRECISION weights[_num_polar/2*_num_azim/4];

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

  /* Compute the product of the sine thetas and weights */
  Quadrature::precomputeWeights(solve_3D);
}


/**
 * @brief Dummy constructor calls the parent constructor.
 */
EqualAnglePolarQuad::EqualAnglePolarQuad(): Quadrature() {
  _quad_type = EQUAL_ANGLE;
  _num_polar = 6;
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
  double thetas[_num_polar/2*_num_azim/4];

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
}


/**
 * @brief Calculates total weights for every azimuthal/polar combination based
 *        on the equal angle polar quadrature.
 * @param solve_3D Boolean indicating whether this is a 3D quadrature
 */
void EqualAnglePolarQuad::precomputeWeights(bool solve_3D) {

  /* Allocate temporary arrays for tabulated quadrature values */
  FP_PRECISION weights[_num_polar/2*_num_azim/4];

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

  /* Compute the product of the sine thetas and weights */
  Quadrature::precomputeWeights(solve_3D);
}
