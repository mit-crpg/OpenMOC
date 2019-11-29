/**
 * @file Quadrature.h
 * @brief The Quadrature abstract class and subclasses.
 * @details The Quadrature subclasses are defined with tabulated or functional
 *          quadrature sets given in the "Lattice Physics Computations",
 *          Handbook of Nuclear Engineering, Dave Knott, Akio Yamamoto, 2010.
 * @date April 8, 2016
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
*/


#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "constants.h"
#include "log.h"
#include <sstream>
#endif
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>


/**
 * @enum QuadratureType
 * @brief The types of quadrature sets supported by OpenMOC.
 */
enum QuadratureType {
  TABUCHI_YAMAMOTO,
  LEONARD,
  GAUSS_LEGENDRE,
  EQUAL_WEIGHT,
  EQUAL_ANGLE
};

/** Shorthand for vectors of floats */
typedef std::vector<double> DoubleVec;
typedef std::vector<FP_PRECISION> FloatVec;
typedef DoubleVec::const_iterator DVCI;

/** Template function for writing vectors to a stream */
template <typename T>
  std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
  if (vec.size() == 0) {
    os << "EMPTY";
  }

  else {

    using CI = typename std::vector<T>::const_iterator;

    os << "[";
    for (CI cit = vec.begin(); cit != vec.end(); ++cit) {
      if (cit != vec.begin()) {
        os << ", ";
      }
      os << *cit;
    }
    os << "]";
  }
  return os;
}



/**
 * @class Quadrature Quadrature.h "src/Quadrature.h"
 * @brief The arbitrary quadrature parent class.
 */
class Quadrature {

protected:

  /** The quadrature type being used */
  QuadratureType _quad_type;

  /** The number of azimuthal angles in (0, 2*PI) */
  size_t _num_azim;

  /** The number of polar angles in (0, PI) */
  size_t _num_polar;

  /** An array of the sines of quadrature polar angles */
  std::vector<DoubleVec> _sin_thetas;

  /** An array of the inverse sines of quadrature polar angles */
  std::vector<FloatVec> _inv_sin_thetas;

  /** An array of the quadrature polar angles */
  std::vector<DoubleVec> _thetas;

  /** An array of the quadrature azimuthal angles */
  DoubleVec _phis;

  /** The actual track azimuthal spacing (cm) by azimuthal angle */
  DoubleVec _azim_spacings;

  /** An array of the quadrature azimuthal weights */
  DoubleVec _azim_weights;

  /** The actual track polar spacing (cm) by (azim, polar) */
  std::vector<DoubleVec> _polar_spacings;

  /** An array of the quadrature polar weights */
  std::vector<DoubleVec> _polar_weights;

  /** An array of the total weights for each azimuthal/polar angle pair */
  std::vector<FloatVec> _total_weights;

  /* Templates for setting the same values to complimentary and supplementary
   * angles */
  template <typename T>
  void setPolarValues(std::vector< std::vector<T> >& vec, size_t azim_index,
                      size_t polar_index, T value) {
    size_t c_azim_index  = _num_azim/2 - azim_index - 1;
    size_t s_polar_index = _num_polar - polar_index - 1;
    vec.at(azim_index).at(polar_index)     = value;
    vec.at(c_azim_index).at(polar_index)   = value;
    vec.at(azim_index).at(s_polar_index)   = value;
    vec.at(c_azim_index).at(s_polar_index) = value;
  }

  template <typename T>
  void setAzimuthalValues(std::vector<T>& vec, size_t azim_index, T value) {
    vec.at(azim_index)                   = value;
    vec.at(_num_azim/2 - azim_index - 1) = value;
  }

  template <typename T>
  static void resize2D(std::vector< std::vector<T> >& vec, size_t dim1,
                       size_t dim2) {
    vec.resize(dim1);
    for (size_t i = 0; i < dim1; ++i)
    {
        vec.at(i).resize(dim2);
    }
  }

public:

  Quadrature();

  size_t getNumPolarAngles() const;
  size_t getNumAzimAngles() const;
  double getSinTheta(size_t azim, size_t polar) const;
  double getSinThetaInline(size_t azim, size_t polar) const;
  FP_PRECISION getInvSinThetaInline(size_t azim, size_t polar) const;
  double getTheta(size_t azim, size_t polar) const;
  double getPhi(size_t azim) const;
  double getAzimWeight(size_t azim) const;
  double getPolarWeight(size_t azim, size_t polar) const;
  double getWeight(size_t azim, size_t polar) const;
  FP_PRECISION getWeightInline(size_t azim, size_t polar) const;
  const std::vector<DoubleVec>& getSinThetas() const;
  const std::vector<DoubleVec>& getThetas() const;
  const DoubleVec& getPhis() const;
  const DoubleVec& getAzimWeights() const;
  const std::vector<DoubleVec>& getPolarWeights() const;
  QuadratureType getQuadratureType() const;
  const DoubleVec& getAzimSpacings() const;
  double getAzimSpacing(size_t azim) const;
  const std::vector<DoubleVec>& getPolarSpacings() const;
  double getPolarSpacing(size_t azim, size_t polar) const;

  virtual void setNumAzimAngles(size_t num_azim);
  virtual void setNumPolarAngles(size_t num_polar);
  void setThetas(const DoubleVec& thetas);
  void setPolarWeights(const DoubleVec& weights);
  void setTheta(double theta, size_t azim, size_t polar);
  void setPhi(double phi, size_t azim);
  void setAzimSpacing(double spacing, size_t azim);
  void setPolarSpacing(double spacing, size_t azim, size_t polar);
  void setAzimWeight(double weight, size_t azim);
  void setPolarWeight(double weight, size_t azim, size_t polar);

  virtual void initialize();
  virtual void precomputeWeights(bool solve_3D);

  std::string toString() const;

  friend std::ostream& operator<<(std::ostream& os, const Quadrature& quad);
};



/**
 * @class TYPolarQuad Quadrature.h "src/Quadrature.h"
 * @brief Tabuchi-Yamamoto's polar quadrature.
 */
class TYPolarQuad: public Quadrature {

private:

public:
  TYPolarQuad();
  void setNumPolarAngles(size_t num_polar);
  void initialize();
  void precomputeWeights(bool solve_3D);
};


/**
 * @class LeonardPolarQuad Quadrature.h "src/Quadrature.h"
 * @brief Leonard's polar quadrature.
 */
class LeonardPolarQuad: public Quadrature {

private:

public:
  LeonardPolarQuad();
  void setNumPolarAngles(size_t num_polar);
  void initialize();
  void precomputeWeights(bool solve_3D);
};



/**
 * @class GLPolarQuad Quadrature.h "src/Quadrature.h"
 * @brief Gauss-Legendre's polar quadrature.
 */
class GLPolarQuad: public Quadrature {

private:

  /** the roots to the Legendre polynomial of degree _num_polar */
  std::vector <double> _roots;

  /** the roots which have been adjusted to allow for periodic tracks */
  std::vector <double> _adjusted_roots;

  /** Whether to updated weights based on adjusted polar angles */
  bool _correct_weights;

public:
  GLPolarQuad();
  void setNumPolarAngles(size_t num_polar);
  void useCorrectedWeights(bool use_corrected_weights);
  void initialize();
  void precomputeWeights(bool solve_3D);
  DoubleVec getCorrectedWeights(size_t azim) const;

  static double legendrePolynomial(size_t n, double x);
  static double logDerivLegendre(size_t n, double x);
  static double secondLogDerivLegendre(size_t n, double x);
  static double getSingleWeight(double root, size_t n);
  static DoubleVec getLegendreRoots(size_t n);
  static DoubleVec getGLWeights(const DoubleVec& roots, size_t n);
};



/**
 * @class EqualWeightPolarQuad Quadrature.h "src/Quadrature.h"
 * @brief Equal weight polar quadrature.
 */
class EqualWeightPolarQuad: public Quadrature {

private:

public:
  EqualWeightPolarQuad();
  void setNumPolarAngles(size_t num_polar);
  void initialize();
  void precomputeWeights(bool solve_3D);
};



/**
 * @class EqualAnglePolarQuad Quadrature.h "src/Quadrature.h"
 * @brief Equal angle polar quadrature.
 */
class EqualAnglePolarQuad: public Quadrature {

private:

public:
  EqualAnglePolarQuad();
  void setNumPolarAngles(size_t num_polar);
  void initialize();
  void precomputeWeights(bool solve_3D);
};


/**
 * @brief Returns the \f$ sin(\theta)\f$ value for a particular polar angle.
 * @param azim index of the azimthal angle of interest
 * @param polar index of the polar angle of interest
 * @return the value of \f$ \sin(\theta) \f$ for this azimuthal and polar angle
 * @details azim must be between 0 and _num_azim / 2
 */
inline double Quadrature::getSinThetaInline(size_t azim, size_t polar) const {

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

  return _sin_thetas[azim][polar];
}


/**
 * @brief Returns the \f$ 1/sin(\theta)\f$ value for a particular polar angle.
 * @param azim index of the azimthal angle of interest
 * @param polar index of the polar angle of interest
 * @return the value of \f$ 1 / \sin(\theta) \f$ for this azimuthal and polar angle
 * @details azim must be between 0 and _num_azim / 2
 */
inline FP_PRECISION Quadrature::getInvSinThetaInline(size_t azim, size_t polar) const {

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

  return _inv_sin_thetas[azim][polar];
}


/**
 * @brief Returns the total weight for Tracks with the given azimuthal and
 *        polar indexes without error checking and inlined
 * @details Angular weights are multiplied by Track spacings
 * @param azim index of the azimuthal angle of interest
 * @param polar index of the polar angle of interest
 * @return the total weight of each Track with the given indexes
 */
inline FP_PRECISION Quadrature::getWeightInline(size_t azim, size_t polar) const {
  return _total_weights[azim][polar];
}

#endif /* QUADRATURE_H_ */
