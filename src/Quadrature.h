/**
 * @file Quadrature.h
 * @brief The Quadrature class.
 * @date January 20, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
*/


#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#ifdef __cplusplus
#include <sstream>
#include "log.h"
#endif

/**
 * @enum quadratureType
 * @brief The type of polar quadrature.
 */
enum quadratureType {
  /** The Leonard type polar quadrature */
  LEONARD,

  /** The Tabuchi-Yamamoto polar quadrature */
  TABUCHI
};


/**
 * @class Quadrature Quadrature.h "openmoc/src/host/Quadrature.h"
 * @brief Stores values for a variety of polar quadratures which may be used.
 */
class Quadrature {

private:
  /** The type of polar quadrature (LEONARD or TABUCHI) */
  quadratureType _type;

  /** The number of polar angles */
  int _num_polar_angles;

  /** An array of the sine of the polar angles from the quadrature set */
  FP_PRECISION* _sinthetas;

  /** An array of the weights for the polar angles from the quadrature set */
  FP_PRECISION* _weights;

  /** An array of the sine multipled by the weight for the polar angles from
   *  the quadrature set */
  FP_PRECISION* _multiples;

public:

  Quadrature(quadratureType type=TABUCHI, int num_polar_angles=3);
  virtual ~Quadrature();

  int getNumPolarAngles() const;
  quadratureType getType() const;
  FP_PRECISION getSinTheta(const int n) const;
  FP_PRECISION getWeight(const int n) const;
  FP_PRECISION getMultiple(const int n) const;
  FP_PRECISION* getSinThetas();
  FP_PRECISION* getWeights();
  FP_PRECISION* getMultiples();
  std::string toString();
};


#endif /* QUADRATURE_H_ */
