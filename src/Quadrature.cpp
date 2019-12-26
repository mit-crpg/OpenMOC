#include "Quadrature.h"

namespace {
  size_t MAX_LG_ITERS = 10000;
}


/**
 * @brief Dummy constructor sets the default number of angles to zero.
 */
Quadrature::Quadrature() :
  _num_azim(0),
  _num_polar(0) {
}


/**
 * @brief Returns the number of polar angles.
 * @return the number of polar angles
 */
size_t Quadrature::getNumPolarAngles() const {
  return _num_polar;
}


/**
 * @brief Returns the number of azimuthal angles.
 * @return the number of azimuthal angles
 */
size_t Quadrature::getNumAzimAngles() const {
  return _num_azim;
}


/**
 * @brief Returns the \f$ sin(\theta)\f$ value for a particular polar angle.
 * @param azim index of the azimthal angle of interest
 * @param polar index of the polar angle of interest
 * @return the value of \f$ \sin(\theta) \f$ for this azimuthal and polar angle
 */
double Quadrature::getSinTheta(size_t azim, size_t polar) const {

  if (polar < 0 || polar >= _num_polar)
    log_printf(ERROR, "Attempted to retrieve sin theta for polar angle = "
               "%d and azim angle = %d when only %d polar angles are "
               "defined", polar, azim, _num_polar);

  else if (azim < 0 || azim >= _num_azim)
    log_printf(ERROR, "Attempted to retrieve sin theta for polar angle = "
               "%d and azim angle = %d when only %d azim angles are "
               "defined", polar, azim, _num_azim);

  else if (_sin_thetas.size() == 0)
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
double Quadrature::getTheta(size_t azim, size_t polar) const {

  if (polar < 0 || polar >= _num_polar)
    log_printf(ERROR, "Attempted to retrieve theta for polar angle = "
               "%d and azim angle = %d when only %d polar angles are "
               "defined", polar, azim, _num_polar);

  else if (azim < 0 || azim >= _num_azim)
    log_printf(ERROR, "Attempted to retrieve theta for polar angle = "
               "%d and azim angle = %d when only %d azim angles are "
               "defined", polar, azim, _num_azim);

  else if (_thetas.size() == 0)
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
double Quadrature::getPhi(size_t azim) const {

  if (azim < 0 || azim >= _num_azim)
    log_printf(ERROR, "Attempted to retrieve theta for "
               "azim angle = %d when only %d azim angles are "
               "defined", azim, _num_azim);

  else if (_phis.size() == 0)
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
double Quadrature::getAzimWeight(size_t azim) const {

  if (azim < 0 || azim >= _num_azim)
    log_printf(ERROR, "Attempted to retrieve the azimuthal weight for "
               "azimuthal angle = %d but only %d azimuthal angles "
               "are defined", azim, _num_azim);

  else if (_azim_weights.size() == 0)
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
double Quadrature::getPolarWeight(size_t azim, size_t polar) const {

  if (polar < 0 || polar >= _num_polar)
    log_printf(ERROR, "Attempted to retrieve polar weight for polar angle = "
               "%d and azim angle = %d when only %d polar angles are "
               "defined", polar, azim, _num_polar);

  if (azim < 0 || azim >= _num_azim)
    log_printf(ERROR, "Attempted to retrieve polar weight for polar angle = "
               "%d and azim angle = %d when only %d azim angles are "
               "defined", polar, azim, _num_azim);

  if (_polar_weights.size() == 0)
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
 * @details Angular weights are multiplied by Track spacings
 * @param azim index of the azimuthal angle of interest
 * @param polar index of the polar angle of interest
 * @return the total weight of each Track with the given indexes
 */
double Quadrature::getWeight(size_t azim, size_t polar) const {

  if (polar < 0 || polar >= _num_polar)
    log_printf(ERROR, "Attempted to retrieve the weight for polar angle = "
               "%d and azimuthal angle = %d but only %d polar angles are "
               "defined", polar, azim, _num_polar);

  else if (azim < 0 || azim >= _num_azim)
    log_printf(ERROR, "Attempted to retrieve the weight for polar angle = "
               "%d and azimuthal angle = %d but only %d azimuthal angles are "
               "defined", polar, azim, _num_azim);

  else if (_total_weights.size() == 0)
    log_printf(ERROR, "Attempted to retrieve weight for polar angle = %d "
               "and azimuthal angle = %d but the multiples have not been "
               "initialized", polar, azim);

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

  return _total_weights[azim][polar];
}


/**
 * @brief Returns a pointer to the Quadrature's vector of polar angle sines
          \f$ sin\theta_{p} \f$.
 * @return a reference to the vector of \f$ sin\theta_{p} \f$
 */
const std::vector<DoubleVec>& Quadrature::getSinThetas() const {

  if (_sin_thetas.size() == 0)
    log_printf(ERROR, "Attempted to retrieve the sin thetas vector "
               "but it has not been initialized");

  return _sin_thetas;
}


/**
 * @brief Returns a reference to the Quadrature's vector of polar angles
          \f$ \theta_{p} \f$.
 * @return a reference to the vector of \f$ \theta_{p} \f$
 */
const std::vector<DoubleVec>& Quadrature::getThetas() const {

  if (_thetas.size() == 0)
    log_printf(ERROR, "Attempted to retrieve the thetas vector "
               "but it has not been initialized");

  return _thetas;
}


/**
 * @brief Returns a pointer to the Quadrature's vector of azimuthal angles
          \f$ \phi \f$.
 * @return a pointer to the vector of \f$ \phi \f$
 */
const DoubleVec& Quadrature::getPhis() const {

  if (_phis.size() == 0)
    log_printf(ERROR, "Attempted to retrieve the phis vector "
               "but it has not been initialized");

  return _phis;
}


/**
 * @brief Returns a pointer to the Quadrature's vector of azimuthal weights.
 * @return a pointer to the azimuthal weights vector
 */
const DoubleVec& Quadrature::getAzimWeights() const {

  if (_azim_weights.size() == 0)
    log_printf(ERROR, "Attempted to retrieve the azimuthal weights vector "
               "but it has not been initialized");

  return _azim_weights;
}


/**
 * @brief Returns a pointer to the Quadrature's vector of polar weights.
 * @return a pointer to the polar weights vector
 */
const std::vector<DoubleVec>& Quadrature::getPolarWeights() const {

  if (_polar_weights.size() == 0)
    log_printf(ERROR, "Attempted to retrieve the polar weights vector "
               "but it has not been initialized");

  return _polar_weights;
}


/**
 * @brief Set the number of azimuthal angles to initialize.
 * @param num_azim the number of azimuthal angles
 */
void Quadrature::setNumAzimAngles(size_t num_azim) {

  if (num_azim <= 0)
    log_printf(ERROR, "Unable to set the number of azimuthal angles to %d "
               "which is less than or equal to zero", num_azim);

  else if (num_azim % 4 != 0)
    log_printf(ERROR, "Unable to set the number of azimuthal angles to %d "
               "which is not divisible by 4", num_azim);

  if (num_azim != _num_azim) {

    /* Clear vectors with old settings */
    _thetas.clear();
    _sin_thetas.clear();
    _phis.clear();
    _azim_spacings.clear();
    _polar_spacings.clear();
    _azim_weights.clear();
    _polar_weights.clear();
    _total_weights.clear();

    _num_azim = num_azim;
  }
}


/**
 * @brief Returns an vector of adjusted azimuthal spacings.
 * @details An vector of azimuthal spacings after adjustment is returned,
 *          indexed by azimuthal angle
 * @return the vector of azimuthal spacings
 */
const DoubleVec& Quadrature::getAzimSpacings() const {
  return _azim_spacings;
}


/**
 * @brief Returns the adjusted azimuthal spacing at the requested azimuthal
 *        angle index.
 * @details The azimuthal spacing depends on the azimuthal angle. This function
 *          returns the azimuthal spacing used at the desired azimuthal angle
 *          index.
 * @param azim the requested azimuthal angle index
 * @return the requested azimuthal spacing
 */
double Quadrature::getAzimSpacing(size_t azim) const {
  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;
  return _azim_spacings[azim];
}

/**
 * @brief Returns a 2D vector of adjusted polar spacings.
 * @details An vector of polar spacings after adjustment is returned,
 *          indexed first by azimuthal angle and then by polar angle
 * @return the 2D vector of polar spacings
 */
const std::vector<DoubleVec>& Quadrature::getPolarSpacings() const {
  return _polar_spacings;
}


/**
 * @brief Returns the adjusted polar spacing at the requested azimuthal
 *        angle index and polar angle index.
 * @details The polar spacing depends on the azimuthal angle and the polar
 *          angle. This function returns the azimuthal spacing used at the
 *          desired azimuthal angle and polar angle indexes.
 * @param azim the requested azimuthal angle index
 * @param polar the requested polar angle index
 * @return the requested polar spacing
 */
double Quadrature::getPolarSpacing(size_t azim, size_t polar) const {
  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;
  return _polar_spacings[azim][polar];
}

/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles
 */
void Quadrature::setNumPolarAngles(size_t num_polar) {

  if (num_polar <= 0)
    log_printf(ERROR, "Unable to set the number of polar angles to %d "
               "which is less than or equal to zero", num_polar);

  if (num_polar % 2 != 0)
    log_printf(ERROR, "Unable to set the number of polar angles to %d "
               "which is not divisible by 2", num_polar);

  if (num_polar != _num_polar) {

    /* Clear vectors with old settings */
    _thetas.clear();
    _sin_thetas.clear();
    _polar_spacings.clear();
    _polar_weights.clear();
    _total_weights.clear();

    _num_polar = num_polar;
  }
}


/**
 * @brief Sets the Quadrature's vector of polar angles.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Quadrature's polar angles in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 vector the length
 *          of the number of azimuthal times polar angles) as input to this
 *          function. This function then fills the Quadrature's polar angles
 *          with the given values. An example of how this function might be
 *          called in Python is as follows:
 *
 * @code
 *          thetas = numpy.array([pi/6, pi/4, pi/3, ... ])
 *          quad = openmoc.Quadrature()
 *          quad.setNumAzimAngles(num_azim)
 *          quad.setNumPolarAngles(len(thetas) / num_azim)
 *          quad.setThetas(thetas)
 * @endcode
 *
 * @param thetas the vector of polar angle for each azimuthal/polar angle
 *        combination
 */
void Quadrature::setThetas(const DoubleVec& thetas) {

  if (_num_polar/2 * _num_azim/4 != thetas.size())
    log_printf(ERROR, "Unable to set %d thetas for Quadrature "
               "with %d polar angles and %d azimuthal angles"
               " in each octant", thetas.size(), _num_polar/2,
               _num_azim/4);

  resize2D(_thetas, _num_azim/2, _num_polar);

  /* Extract sin thetas from user input */
  size_t ap=0;
  for (size_t a=0; a < _num_azim/4; ++a) {
    for (size_t p=0; p < _num_polar/2; ++p) {
      if (thetas[ap] < 0. || thetas[ap] > M_PI_2)
        log_printf(ERROR, "Unable to set theta to %f which is "
                   "not in the range [0,PI/2]", thetas[ap]);

      _thetas[a][p] = thetas[ap];
      _thetas[_num_azim/2 - a - 1][p] = thetas[ap];
      _thetas[a][_num_polar - p - 1] = M_PI - thetas[ap];
      _thetas[_num_azim/2 - a - 1][_num_polar - p - 1] = M_PI - thetas[ap];
      ++ap;
    }
  }
}


/**
 * @brief Set the Quadrature's vector of polar weights.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Quadrature's polar weights in Python. A user must initialize a
 *          NumPy array of the correct size (e.g., a float64 vector the length
 *          of the number of azimuthal times polar angles) as input to this
 *          function. This function then fills the Quadrature's polar weights
 *          with the given values. An example of how this function might be
 *          called in Python is as follows:
 *
 * @code
 *          polar_weights = numpy.array([0.1, 0.2, 0.05, ... ])
 *          quad = openmoc.Quadrature()
 *          quad.setNumAzimAngles(num_azim)
 *          quad.setNumPolarAngles(len(polar_weights) / num_azim)
 *          quad.setPolarWeights(polar_weights)
 * @endcode
 *
 * @param weights The polar weights
 */
void Quadrature::setPolarWeights(const DoubleVec& weights) {

  if (_num_polar/2 * _num_azim/4 != weights.size())
    log_printf(ERROR, "Unable to set %d polar weights for Quadrature "
               "with %d polar angles and %d azimuthal angles"
               " in each octant",
               weights.size(), _num_polar/2, _num_azim/4);

  /* Initialize memory for vectors */
  resize2D(_polar_weights, _num_azim/2, _num_polar);

  /* Extract polar weights from user input */
  size_t ap=0;
  for (size_t a=0; a < _num_azim/4; ++a) {
    for (size_t p=0; p < _num_polar/2; ++p) {
      if (weights[ap] < 0. || weights[ap] > M_PI_2)
        log_printf(ERROR, "Unable to set polar weight to %f which is "
                   "not in the range [0,PI/2]", weights[ap]);

      setPolarValues(_polar_weights, a, p, weights[ap]);
      ++ap;
    }
  }
}


/**
 * @brief Sets the polar angle for the given indexes.
 * @param theta the value in radians of the polar angle to be set
 * @param azim the azimuthal index of the angle of interest
 * @param polar the polar index of the angle of interest
 */
void Quadrature::setTheta(double theta, size_t azim, size_t polar) {

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

  resize2D(_thetas, _num_azim/2, _num_polar);

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
void Quadrature::setPhi(double phi, size_t azim) {

  if (phi <= 0.0 || phi >= M_PI_2)
    log_printf(ERROR, "Unable to set phi for azim = %d to %f which is not "
               "in the range (0.0, PI/2)", azim, phi);

  if (azim >= _num_azim/4)
    log_printf(ERROR, "Unable to set phi for azim = %d since azim is not in"
               " the range (0, _num_azim/4)", azim);

  _phis.resize(_num_azim/2);

  _phis[azim] = phi;
  _phis[_num_azim/2 - azim - 1] = M_PI - phi;
}


/**
 * @brief Sets the azimuthal spacing for the given index.
 * @param spacing the spacing (cm) in the azimuthal direction to be set
 * @param azim the azimuthal index
 */
void Quadrature::setAzimSpacing(double spacing, size_t azim) {

  if (spacing <= 0.0)
    log_printf(ERROR, "Unable to set azimuthal spacing for azim = %d to %f "
                      "which is not strictly greater than zero", azim,
                      spacing);

  if (azim >= _num_azim/4)
    log_printf(ERROR, "Unable to set azimuthal spacing for azim = %d since "
                      " azim is not in the range (0, _num_azim/4)", azim);

  _azim_spacings.resize(_num_azim/2);

  setAzimuthalValues(_azim_spacings, azim, spacing);
}


/**
 * @brief Sets the polar spacing for the given indexes.
 * @param spacing the spacing in the polar direction to be set
 * @param azim the azimuthal index corresponding to the angle
 * @param polar the polar index corresponding to the angle
 */
void Quadrature::setPolarSpacing(double spacing, size_t azim, size_t polar) {

  if (spacing <= 0)
    log_printf(ERROR, "Unable to set polar spacing for azim = %d and polar = "
                      "%d to %f which is not strictly greater than zero", azim,
                      polar, spacing);

  if (azim >= _num_azim/4)
    log_printf(ERROR, "Unable to set polar spacing for azim = %d and polar = "
                      "%d since azim is not in the range (0, _num_azim/4)",
                      azim, polar);
  if (polar >= _num_polar/2)
    log_printf(ERROR, "Unable to set polar spacing for azim = %d and polar = "
                      "%d since polar is not in the range (0, _num_polar/2)",
                      azim, polar);

  resize2D(_polar_spacings, _num_azim/2, _num_polar);
  setPolarValues(_polar_spacings, azim, polar, spacing);
}


/**
 * @brief Sets the azimuthal weight for the given index.
 * @param weight the weight of the azimuthal angle
 * @param azim the azimuthal index
 */
void Quadrature::setAzimWeight(double weight, size_t azim) {

  if (weight <= 0.0 || weight >= M_PI_2)
    log_printf(ERROR, "Unable to set azim weight for azim = %d to %f which is "
               "not in the range (0.0, PI/2)", azim, weight);

  if (azim >= _num_azim/4)
    log_printf(ERROR, "Unable to set azim weight for azim = %d since azim is "
               "not in the range (0, _num_azim/4)", azim);

  _azim_weights.resize(_num_azim/2);

  setAzimuthalValues(_azim_weights, azim, double(weight));
}


/**
 * @brief Sets the polar weight for the given indexes.
 * @param weight the weight of the polar angle
 * @param azim the azimuthal index corresponding to the angle
 * @param polar the polar index corresponding to the angle
 */
void Quadrature::setPolarWeight(double weight, size_t azim, size_t polar) {

  if (weight <= 0.0 || weight >= M_PI_2) {
    log_printf(ERROR, "Unable to set polar weight for azim = %d and "
               "polar = %d to %f which is not in the range (0.0, PI/2)",
               azim, polar, weight);
  }

  if (azim >= _num_azim/4) {
    log_printf(ERROR, "Unable to set polar weight for azim = %d and polar = %d "
               "since azim is not in the range (0, _num_azim/4)", azim, polar);
  }

  if (polar >= _num_polar/2) {
    log_printf(ERROR, "Unable to set polar weight for azim = %d and polar = %d "
               "since polar is not in the range (0, _num_polar/2)", \
               azim, polar);
  }

  resize2D(_polar_weights, _num_azim/2, _num_polar);

  setPolarValues(_polar_weights, azim, polar, double(weight));
}


/**
 * @brief Initialize the polar quadrature azimuthal angles.
 * @details The parent class routine simply checks that number of polar and
 *          azimuthal angles have been set by the user and generates the
 *          azimuthal angles if not already generated.
 */
void Quadrature::initialize() {

  if (_num_polar == 0) {
    log_printf(ERROR, "Unable to initialize Quadrature with zero polar angles. "
               "Set the number of polar angles before initialization.");
  }

  if (_num_azim == 0) {
    log_printf(ERROR, "Unable to initialize Quadrature with zero azimuthal "
               "angles. Set the number of azimuthal angles before "
               "initialization.");
  }

  _phis.resize(_num_azim/2);

  /* Compute a desired set of azimuthal angles */
  for (size_t a = 0; a < _num_azim/2; ++a) {
    _phis[a] = 2.0 * M_PI / _num_azim * (0.5 + a);
  }
}


/**
 * @brief This private routine computes the product of the sine thetas and
 *        weights for each angle in the polar quadrature.
 * @details Note that this routine must be called after populating the
 *          sine thetas and weights vectors.
 */
void Quadrature::precomputeWeights(bool solve_3D) {

  /* Check that track spacings have been set */
  if (_azim_spacings.size() == 0)
    log_printf(ERROR, "Unable to precompute weights since track spacings have "
                      "not yet been set");

  /* Check that polar angles have been set */
  if (_thetas.size() == 0)
    log_printf(ERROR, "Unable to precompute weights since polar angles have "
                      "not yet been set");

  /* Check that polar spacings have been set in 3D cases */
  if (solve_3D && _polar_spacings.size() == 0)
    log_printf(ERROR, "Unable to precompute weights for the 3D quadrature since"
                      " polar spacings have not yet been set");

  /* Clear azimuthal weights */
  _azim_weights.resize(_num_azim/2);

  /* Create uncorrected weights if no angles have been set yet */
  if (_phis.size() == 0) {
    log_printf(NORMAL, "WARNING: Using uncorrected angles for weights");
    double phi = M_PI / _num_azim;
    for (size_t a = 0; a < _num_azim/4; ++a) {
      setPhi(phi, a);
      phi += 2*M_PI / _num_azim;
    }
  }

  /* Compute the azimuthal weights */
  double x1, x2;
  for (size_t a = 0; a < _num_azim/4; ++a) {

    /* The azimuthal weights (in radians) using equal weight quadrature */
    if (a < _num_azim/4 - 1)
      x1 = 0.5 * (_phis[a+1] - _phis[a]);
    else
      x1 = M_PI_2 - _phis[a];

    if (a >= 1)
      x2 = 0.5 * (_phis[a] - _phis[a-1]);
    else
      x2 = _phis[a];

    setAzimuthalValues(_azim_weights, a, double((x1 + x2) / M_PI));
  }

  /* Allocate memory if it was not allocated previously */
  resize2D(_sin_thetas, _num_azim/2, _num_polar);
  resize2D(_inv_sin_thetas, _num_azim/2, _num_polar);
  resize2D(_total_weights, _num_azim/2, _num_polar);

  /* Compute multiples of sine thetas and weights */
  for (size_t a=0; a < _num_azim/4; ++a) {
    for (size_t p=0; p < _num_polar/2; ++p) {
      double sin_theta = sin(_thetas[a][p]);
      double weight = 2.0 * M_PI * _azim_weights[a] * _azim_spacings[a]
          * _polar_weights[a][p];
      if (solve_3D)
        weight *= _polar_spacings[a][p];
      else
        weight *= 2.0 * sin_theta;
      setPolarValues(_sin_thetas, a, p, sin_theta);
      setPolarValues(_inv_sin_thetas, a, p, FP_PRECISION(1.f / sin_theta));
      setPolarValues(_total_weights, a, p, FP_PRECISION(weight));
    }
  }
}


/**
 * @brief Converts this Quadrature to a character vector of its attributes.
 * @details The character vector includes the number of polar angles, the
 *          the values of the sine and weight of each polar angle, and the
 *          product of the sine and weight of each polar angle.
 * @return a character vector of the Quadrature's attributes
 */
std::string Quadrature::toString() const {

  std::stringstream string;

  string << "Quadrature";
  string << "\n\t# azim angles  = " << _num_azim;
  string << "\n\t# polar angles = " << _num_polar;

  string << "\n\tphis = " << _phis;
  string << "\n\tazim weights = " << _azim_weights;
  string << "\n\tthetas = " << _thetas;

  string << "\n\tpolar weights = " << _polar_weights;
  string << "\n\tsin thetas = " << _sin_thetas;
  string << "\n\ttotal weights = " << _total_weights;

  return string.str();
}


/**
 * @brief Prints to the provided output stream
 * @details Allows printing the Quadrature using <<
 * @param os the provided stream to write to
 * @param quad the quadrature object which is printed
 * @return the provided stream
 */
std::ostream& operator<<(std::ostream& os, const Quadrature& quad) {
  os << quad.toString();
  return os;
}


/**
 * @brief Returns the type of Quadrature created.
 * @return The quadrature type
 */
QuadratureType Quadrature::getQuadratureType() const {
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
void TYPolarQuad::setNumPolarAngles(size_t num_polar) {

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

  /* Allocate temporary vectors for tabulated quadrature values */
  DoubleVec thetas(_num_polar/2 * _num_azim/4);

  /* Tabulated values for the sine thetas and weights for the
   * Tabuchi-Yamamoto polar angle quadrature */
  if (_num_polar == 2) {
    for (size_t a=0; a < _num_azim/4; ++a) {
      thetas[a] = asin(0.798184);
    }
  }

  else if (_num_polar == 4) {
    for (size_t a=0; a < _num_azim/4; ++a) {
      thetas[a*(_num_polar/2)] = asin(0.363900);
      thetas[a*(_num_polar/2)+1] = asin(0.899900);
    }
  }

  else if (_num_polar == 6) {
    for (size_t a=0; a < _num_azim/4; ++a) {
      thetas[a*(_num_polar/2)] = asin(0.166648);
      thetas[a*(_num_polar/2)+1] = asin(0.537707);
      thetas[a*(_num_polar/2)+2] = asin(0.932954);
    }
  }

  /* Set the vectors of thetas */
  Quadrature::setThetas(thetas);
}


/**
 * @brief Calculates total weights for every azimuthal/polar combination based
 *        on the TY quadrature.
 * @param solve_3D Boolean indicating whether this is a 3D quadrature
 */
void TYPolarQuad::precomputeWeights(bool solve_3D) {

  /* Allocate temporary vectors for tabulated quadrature values */
  DoubleVec weights(_num_polar/2 * _num_azim/4);

  /* Tabulated values for the sine thetas and weights for the
   * Tabuchi-Yamamoto polar angle quadrature */
  if (_num_polar == 2) {
    for (size_t a=0; a < _num_azim/4; ++a) {
      weights[a] = 0.5;
    }
  }

  else if (_num_polar == 4) {
    for (size_t a=0; a < _num_azim/4; ++a) {
      weights[a*(_num_polar/2)] = 0.212854 / 2.0;
      weights[a*(_num_polar/2)+1] = 0.787146 / 2.0;
    }
  }

  else if (_num_polar == 6) {
    for (size_t a=0; a < _num_azim/4; ++a) {
      weights[a*(_num_polar/2)] = 0.046233 / 2.0;
      weights[a*(_num_polar/2)+1] = 0.283619 / 2.0;
      weights[a*(_num_polar/2)+2] = 0.670148 / 2.0;
    }
  }

  /* Set the vectors of sin thetas and weights */
  Quadrature::setPolarWeights(weights);
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
void LeonardPolarQuad::setNumPolarAngles(size_t num_polar) {

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

  /* Allocate temporary vectors for tabulated quadrature values */
  DoubleVec thetas(_num_polar/2 * _num_azim/4);

  /* Tabulated values for the sine thetas and weights for the
   * Leonard polar angle quadrature */
  if (_num_polar == 4) {
    for (size_t a=0; a < _num_azim/4; ++a) {
      thetas[a*(_num_polar/2)] = asin(0.273658);
      thetas[a*(_num_polar/2)+1] = asin(0.865714);
    }
  }

  else if (_num_polar == 6) {
    for (size_t a=0; a < _num_azim/4; ++a) {
      thetas[a*(_num_polar/2)] = asin(0.099812);
      thetas[a*(_num_polar/2)+1] = asin(0.395534);
      thetas[a*(_num_polar/2)+2] = asin(0.891439);
    }
  }

  /* Set the vectors of thetas and weights */
  Quadrature::setThetas(thetas);
}


/**
 * @brief Calculates total weights for every azimuthal/polar combination based
 *        on the Leonard polar quadrature.
 * @param solve_3D Boolean indicating whether this is a 3D quadrature
 */
void LeonardPolarQuad::precomputeWeights(bool solve_3D) {

  /* Allocate temporary vectors for tabulated quadrature values */
  DoubleVec weights(_num_polar/2 * _num_azim/4);

  /* Tabulated values for the sine thetas and weights for the
   * Leonard polar angle quadrature */
  if (_num_polar == 4) {
    for (size_t a=0; a < _num_azim/4; ++a) {
      weights[a*(_num_polar/2)] = 0.139473 / 2.0;
      weights[a*(_num_polar/2)+1] = 0.860527 / 2.0;
    }
  }

  else if (_num_polar == 6) {
    for (size_t a=0; a < _num_azim/4; ++a) {
      weights[a*(_num_polar/2)] = 0.017620 / 2.0;
      weights[a*(_num_polar/2)+1] = 0.188561 / 2.0;
      weights[a*(_num_polar/2)+2] = 0.793819 / 2.0;
    }
  }

  /* Set the vectors of thetas and weights */
  Quadrature::setPolarWeights(weights);
  Quadrature::precomputeWeights(solve_3D);
}



/**
 * @brief Dummy constructor calls the parent constructor.
 */
GLPolarQuad::GLPolarQuad(): Quadrature() {
  _quad_type = GAUSS_LEGENDRE;
  _num_polar = 6;
  _correct_weights = false;
}


/**
 * @brief Set the number of polar angles to initialize.
 * @param num_polar the number of polar angles (maximum 20)
 */
void GLPolarQuad::setNumPolarAngles(size_t num_polar) {

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

  /* Allocate temporary vectors for tabulated quadrature values */
  DoubleVec thetas(_num_polar/2 * _num_azim/4);

  /* get roots of Legendre polynomial */
  _roots = getLegendreRoots(_num_polar);

  /* Set theta values for polar angles */
  for (size_t a=0; a < _num_azim/4; ++a) {
    for (size_t i=0; i < _num_polar/2; ++i) {
      thetas[a*(_num_polar/2)+i] = acos(_roots[i]);
    }
  }

  /* Set the vectors of sin thetas and weights */
  Quadrature::setThetas(thetas);
}


/**
 * @brief Indicates whether to correct weights based on altered polar angles
 * @param use_corrected_weights Whether to alter the weights
 */
void GLPolarQuad::useCorrectedWeights(bool use_corrected_weights) {
  _correct_weights = use_corrected_weights;
}


/**
 * @brief Calculates total weights for every azimuthal/polar combination based
 *        on the Gauss-Legendre polar quadrature.
 * @param solve_3D Boolean indicating whether this is a 3D quadrature
 */
void GLPolarQuad::precomputeWeights(bool solve_3D) {

  /* Get uncorrected weights */
  std::vector <double> weights_vec = getGLWeights(_roots, _num_polar);

  /* Allocate temporary vectors for tabulated quadrature values */
  DoubleVec weights(_num_polar/2*_num_azim/4);

  /* Determine gauss-legendre weights from computed values */
  for (size_t a=0; a < _num_azim/4; ++a) {

    if (_correct_weights)
      weights_vec = getCorrectedWeights(a);

    for (size_t p=0; p < _num_polar/2; ++p) {
      weights[a*(_num_polar/2)+p] = weights_vec[p] / 2.0;
    }
  }

  /* Set the vectors of sin thetas and weights */
  Quadrature::setPolarWeights(weights);
  Quadrature::precomputeWeights(solve_3D);
}


/**
 * @brief the Legendre polynomial of degree n evaluated at x
 * @param n an integer >=0: the order of the polynomial
 * @param x in (-1,1), the point at which to evaluate the polynomial
 * @return the value of the Legendre polynomial of degree n at x
 */
double GLPolarQuad::legendrePolynomial(size_t n, double x) {
  if (n == 0)
    return 1;
  if (n == 1)
    return x;
  else {
    double c = 2.0*(n-2) + 3.0;
    double b = 1.0*(n-2) + 2.0;
    double a = 1.0*(n-2) + 1.0;
    double value = c/b * x * legendrePolynomial(n-1, x)
        - a/b * legendrePolynomial(n-2,x);
    return value;
  }
}


/**
 * @brief The first logarithmic derivative of a Legendre polynomial
 * @param n the order of the polynomial
 * @param x point at which to evaluate the logarithmic derivative
 * @return the value of the logarithmic derivative at x
 */
double GLPolarQuad::logDerivLegendre(size_t n, double x) {
  double num = n * x - n * legendrePolynomial(n-1,x) / legendrePolynomial(n,x);
  double denom = x*x - 1;
  return num/denom;
}


/**
 * @brief The second logarithmic derivative of a Legendre polynomial
 * @param n the order of the polynomial
 * @param x point at which to evaluate the logarithmic derivative
 * @return the value of the logarithmic derivative at x
 */
double GLPolarQuad::secondLogDerivLegendre(size_t n, double x) {
  double num =
    n*(n+1) + logDerivLegendre(n,x) * ((1-x*x)* logDerivLegendre(n,x) - 2 * x);
  double denom = x*x-1;
  return num/denom;
}


/**
 * @brief Finds the roots of Legendre polynomial of order n.
 * @details Guesses for positive roots are set at logarithmic intervals.
 *          Positive roots are found simultaneously using an
 *          Alberth-Householder-n method. Each guess is successively nudged
 *          towards a true root. Only the positive roots are calculated
 * @param n the order of the polynomial
 * @return a list of the roots of the polynomial
 */
DoubleVec GLPolarQuad::getLegendreRoots(size_t n) {

  /* desired precision on roots */
  double E1 = 1e-8;
  double E2 = 1e-8;

  DoubleVec roots;
  DoubleVec s1_tilde;
  DoubleVec s2_tilde;
  std::vector<bool> converged;

  /* set guesses with log scale */
  for (size_t i=0; i < n/2; ++i) {
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

  /* use the Alberth-Householder_n method to nudge guesses towards roots */
  for (size_t iter=0; iter < MAX_LG_ITERS; ++iter) {

    /* set S tildes */
    for (size_t i=0; i < (n+1)/2; ++i) {
      if (!converged[i]) {
        double sum1 = 0;
        double sum2 = 0;
        for (size_t j=0; j <= (n+1)/2; ++j) {
          if (j != i) {
            double diff = (roots[i] - roots[j]);
            if (fabs(diff) > FLT_EPSILON) {
              sum1 += 1. / diff;
              sum2 += -1. / (diff*diff);
            }
          }
        }

        s1_tilde[i] = logDerivLegendre(n, roots[i]) - sum1;
        s2_tilde[i] = secondLogDerivLegendre(n, roots[i]) - sum2;

        /* Householder method 2 Halley */
        double denom = s1_tilde[i]*s1_tilde[i] - s2_tilde[i];
        double u_new = 0.0;
        if (fabs(denom) > FLT_EPSILON)
          u_new = roots[i] - 2*s1_tilde[i] / denom;

        double u_old = roots[i];
        roots[i] = u_new;

        /* if this is the actual root */
        if (std::abs(u_new - u_old) < E1) {
          if (std::abs(legendrePolynomial(n, u_new)) < E2) {
            converged[i] = true;

            /* if this root equals another root or it is less than 0 */
            for (size_t j=0; j < (n+1)/2; ++j) {
              if (j != i) {
                if (std::abs(roots[j] - roots[i]) < E1 || roots[i] <= 0) {

                  /* reset the root to its original guess */
                  roots[i] = - pow(2, (-.5*(i+1))) + 1;
                  converged[i] = false;
                }
              }
            }
          }
        } /* if this is the actual root */
      } /* if not converged */
    } /* for each guess */

    /* check for convergence */
    for (size_t i=0; i<(n+1)/2; ++i) {
      all_roots_converged = converged[i];
      if (!all_roots_converged)
        break;
    }

    if (all_roots_converged)
      break;
    else if (iter == MAX_LG_ITERS - 1) {
      log_printf(ERROR, "Failed to converge Gauss-Legendre roots for %d roots",
                         n);
    }

  } /* while not all roots converged */

  /* Remove placeholder root */
  roots.pop_back();

  /* Put roots in order */
  std::sort (roots.begin(), roots.end());
  return roots;
}


/**
 * @brief Calculates the weights to be used in Gauss-Legendre Quadrature.
 * @param roots a vector containing the roots of the Legendre polynomial
 * @param n the order of the Legendre Polynomial
 * @return a vector of weights matched by index to the vector of roots
 */
DoubleVec GLPolarQuad::getGLWeights(const DoubleVec& roots, size_t n) {

  DoubleVec weights;
  for (size_t i=0; i<roots.size(); ++i) {
    double value = - (2*roots[i]*roots[i] - 2) /
        (n*n*legendrePolynomial(n-1, roots[i])
         *legendrePolynomial(n-1, roots[i]));
    weights.push_back(value);
  }

  return weights;
}

/**
 * @brief Calculates the weights to be used in Gauss-Legendre Quadrature.
 * @param root the root of the Legendre polynomial
 * @param n the order of the Legendre Polynomial
 * @return a vector of weights matched by index to the vector of roots
 */
double GLPolarQuad::getSingleWeight(double root, size_t n) {
  double weight = - (2*root*root - 2) /
      (n*n*legendrePolynomial(n-1, root) * legendrePolynomial(n-1, root));

  return weight;
}

/**
 * @brief Calculates the weights to be used in numerical integration.
 * @details assumes the function will be integrated over (-1, 1)
 * @details azim the azimuthal angle index
 * @return the vector of weights
 */
DoubleVec GLPolarQuad::getCorrectedWeights(size_t azim) const {

  /* Calculate abscissa */
  std::vector <double> nodes;
  for (size_t p=0; p<_num_polar/2; ++p) {
    nodes.push_back(cos(_thetas[azim][p]));
  }
  for (size_t p=0; p<_num_polar/2; ++p) {
    nodes.push_back(-nodes[p]);
  }
  size_t n = nodes.size();

  std::vector<double> weights;

  // declare a vector to store the elements of the augmented-matrix
  std::vector< std::vector<long double> > A;
  resize2D(A, n, n);

  // the solution vector
  std::vector<long double> x(n);

  // index array, used to keep track of the order in which the nodes were passed
  std::vector<int> index(n);
  for (size_t i=0; i<n; ++i)
    index[i] = i;

  long double a = -1;
  long double b = 1;

  // populate A
  for (size_t i=0; i<n; ++i) {
    for (size_t j=0; j<n; ++j) {
      A[i][j] = pow(nodes[j], i);
    }
    A[i][n] = (pow(b, i+1) - pow(a, i+1)) / (i+1);
  }

  /* Select pivots */
  for (size_t i=0; i<n; ++i) {
    for (size_t k=i+1; k<n; ++k) {
      if (A[i][i] < A[k][i]) {

        // switch column k and column i
        for (size_t j=0; j<=n; ++j) {
          long double temp = A[i][j];
          A[i][j] = A[k][j];
          A[k][j] = temp;
        }

        // switch their corresponding indices
        long double temp_ind = index[i];
        index[i] = index[k];
        index[k] = temp_ind;
      }
    }
  }

  // perform gauss elimination
  for (size_t i=0; i<n-1; ++i) {
    for (size_t k=i+1; k<n; ++k) {
      long double t = A[k][i] / A[i][i];

      // make elements below the pivot elements equal to zero or eliminate the
      // variables
      for (size_t j=0; j<=n; ++j)
        A[k][j] = A[k][j] - t * A[i][j];
     }
  }

  // back-substitution
  for (size_t i=n-1; i>=0; --i) {
    long double sub = 0;
    for (size_t j=n-1; j>i; --j) {
      sub += x[j]*A[i][j];
    }
    x[i] = (A[i][n] - sub) / A[i][i];
  }

  for (size_t i=0; i<n; ++i){
    weights.push_back(double(x[i]));
  }

  // fill the vector of weights so that it is indexed to the original
  // vector of nodes
  for (size_t i=0; i<n; ++i)
    weights[index[i]] = double(x[i]);

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
void EqualWeightPolarQuad::setNumPolarAngles(size_t num_polar) {
  Quadrature::setNumPolarAngles(num_polar);
}


/**
 * @brief Routine to initialize the polar quadrature.
 * @details This routine generates the sine thetas and weights.
 */
void EqualWeightPolarQuad::initialize() {

  /* Call parent class initialize routine */
  Quadrature::initialize();

  /* Allocate temporary vectors for tabulated quadrature values */
  DoubleVec thetas(_num_polar/2 * _num_azim/4);

  double cos_theta_a, cos_theta_b;

  /* Generate the sin thetas and weights using equations 420-422 of the
   * DOE Nucl. Eng. Handbook "Lattice Physics Computations" */
  size_t ap = 0;
  for (size_t a=0; a < _num_azim/4; ++a) {
    cos_theta_a = 1.;
    for (size_t p=0; p < _num_polar/2; ++p) {
      cos_theta_b = cos_theta_a - (1. / (_num_polar/2));
      thetas[ap] = acos(0.5 * (cos_theta_a + cos_theta_b));
      cos_theta_a = cos_theta_b;
      ++ap;
    }
  }

  /* Set the vectors of sin thetas and weights */
  Quadrature::setThetas(thetas);
}


/**
 * @brief Calculates total weights for every azimuthal/polar combination based
 *        on the equal weight polar quadrature.
 * @param solve_3D Boolean indicating whether this is a 3D quadrature
 */
void EqualWeightPolarQuad::precomputeWeights(bool solve_3D) {

  /* Allocate temporary vectors for tabulated quadrature values */
  DoubleVec weights(_num_polar/2 * _num_azim/4);

  double y1, y2;

  /* Generate the sin thetas and weights using equations 420-422 of the
   * DOE Nucl. Eng. Handbook "Lattice Physics Computations" */
  size_t ap = 0;
  for (size_t a=0; a < _num_azim/4; ++a) {
    for (size_t p=0; p < _num_polar/2; ++p) {

      if (p < _num_polar/2 - 1)
        y1 = 0.5 * (cos(_thetas[a][p]) - cos(_thetas[a][p+1]));
      else
        y1 = cos(_thetas[a][p]);

      if (p >= 1)
        y2 = 0.5 * (cos(_thetas[a][p-1]) - cos(_thetas[a][p]));
      else
        y2 = 1.0 - cos(_thetas[a][p]);

      weights[ap] = (y1 + y2) / 2.0;
      ++ap;
    }
  }

  /* Set the vectors of sin thetas and weights */
  Quadrature::setPolarWeights(weights);

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
void EqualAnglePolarQuad::setNumPolarAngles(size_t num_polar) {
  Quadrature::setNumPolarAngles(num_polar);
}


/**
 * @brief Routine to initialize the polar quadrature.
 * @details This routine generates the sine thetas and weights.
 */
void EqualAnglePolarQuad::initialize() {

  /* Call parent class initialize routine */
  Quadrature::initialize();

  /* Allocate temporary vectors for tabulated quadrature values */
  DoubleVec thetas(_num_polar/2 * _num_azim/4);

  double cos_theta_a, cos_theta_b;
  double theta_a, theta_b;
  double delta_theta = M_PI / (_num_polar);

  /* Generate the sin thetas and weights using equations 420-422 of the
   * DOE Nucl. Eng. Handbook "Lattice Physics Computations" */
  size_t ap = 0;
  for (size_t a=0; a < _num_azim/4; ++a) {
    theta_a = 0.;
    for (size_t p=0; p < _num_polar/2; ++p) {
      theta_b = theta_a + delta_theta;
      cos_theta_a = cos(theta_a);
      cos_theta_b = cos(theta_b);
      thetas[ap] = acos((0.5 * (cos_theta_a + cos_theta_b)));
      theta_a = theta_b;
      ++ap;
    }
  }

  /* Set the vectors of sin thetas and weights */
  Quadrature::setThetas(thetas);
}


/**
 * @brief Calculates total weights for every azimuthal/polar combination based
 *        on the equal angle polar quadrature.
 * @param solve_3D Boolean indicating whether this is a 3D quadrature
 */
void EqualAnglePolarQuad::precomputeWeights(bool solve_3D) {

  /* Allocate temporary vectors for tabulated quadrature values */
  DoubleVec weights(_num_polar/2 * _num_azim/4);

  double y1, y2;

  /* Generate the sin thetas and weights using equations 420-422 of the
   * DOE Nucl. Eng. Handbook "Lattice Physics Computations" */
  size_t ap = 0;
  for (size_t a=0; a < _num_azim/4; ++a) {
    for (size_t p=0; p < _num_polar/2; ++p) {

      if (p < _num_polar/2 - 1)
        y1 = 0.5 * (cos(_thetas[a][p]) - cos(_thetas[a][p+1]));
      else
        y1 = cos(_thetas[a][p]);

      if (p >= 1)
        y2 = 0.5 * (cos(_thetas[a][p-1]) - cos(_thetas[a][p]));
      else
        y2 = 1.0 - cos(_thetas[a][p]);

      weights[ap] = (y1 + y2) / 2.0;
      ++ap;
    }
  }

  /* Set the vectors of sin thetas and weights */
  Quadrature::setPolarWeights(weights);

  /* Compute the product of the sine thetas and weights */
  Quadrature::precomputeWeights(solve_3D);
}
