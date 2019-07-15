/**
 * @file Point.h
 * @brief The Point class.
 * @date January 18, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef POINT_H_
#define POINT_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "log.h"
#include <math.h>
#include <sstream>
#endif

/**
 * @class Point Point.h "src/Point.h"
 * @brief Class to represent a 2D/3D point in space.
 */
class Point {

private:

  /** The Point's xyz coordinates */
  double _xyz[3];

public:
  Point();
  virtual ~Point();
  void setCoords(const double x, const double y, const double z=0.0);
  double getX() const;
  double getY() const;
  double getZ() const;
  double* getXYZ();
  void setX(const double x);
  void setY(const double y);
  void setZ(const double z);
  double distance(const double x, const double y, const double z) const;
  double distanceToPoint(const Point* point);
  void copyCoords(Point* point);
  std::string toString();
  void setXYZ(double* xyz);
};


/**
 * @brief Initializes a Point with two-dimensional coordinates.
 * @param x x-coordinate
 * @param y y-coordinate
 * @param z z-coordinate
 */
inline void Point::setCoords(const double x, const double y, const double z) {
  _xyz[0] = x;
  _xyz[1] = y;
  _xyz[2] = z;
}


/**
 * @brief Returns this Point's x-coordinate.
 * @return the x-coordinate
 */
inline double Point::getX() const {
  return _xyz[0];
}


/**
 * @brief Returns this Point's y-coordinate.
 * @return the y-coordinate
 */
inline double Point::getY() const {
  return _xyz[1];
}


/**
 * @brief Returns this Point's z-coordinate.
 * @return the z-coordinate
 */
inline double Point::getZ() const {
  return _xyz[2];
}


/**
 * @brief Returns this Point's coordinates.
 * @return the xyz coordinates
 */
inline double* Point::getXYZ() {
  return _xyz;
}


/**
 * @brief Set the Point's x-coordinate.
 * @param x the new x-coordinate
 */
inline void Point::setX(const double x) {
  _xyz[0] = x;
}


/**
 * @brief Set the Point's y-coordinate.
 * @param y the new y-coordinate
 */
inline void Point::setY(const double y) {
  _xyz[1] = y;
}


/**
 * @brief Set the Point's z-coordinate.
 * @param z the new z-coordinate
 */
inline void Point::setZ(const double z) {
  _xyz[2] = z;
}


/**
 * @brief Set the Point's x, y and z-coordinates.
 * @param xyz array with the three coordinates
 */
inline void Point::setXYZ(double* xyz) {
  _xyz[0] = xyz[0];
  _xyz[1] = xyz[1];
  _xyz[2] = xyz[2];
}


/**
 * @brief Compute the distance from this Point to another Point of interest.
 * @param x the x-coordinate of the Point of interest
 * @param y the y-coordinate of the Point of interest
 * @param z the z-coordinate of the Point of interest
 * @return distance to the Point of interest
 */
inline double Point::distance(const double x, const double y, const double z) const {
  double deltax = _xyz[0] - x;
  double deltay = _xyz[1] - y;
  double deltaz = _xyz[2] - z;
  return sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);
}


/**
 * @brief Compute the distance from this Point to another Point of interest.
 * @param point a pointer to the Point of interest
 * @return distance to the Point of interest
 */
inline double Point::distanceToPoint(const Point* point) {
  double deltax = _xyz[0] - point->getX();
  double deltay = _xyz[1] - point->getY();
  double deltaz = _xyz[2] - point->getZ();
  return sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);
}


/**
 * @brief Copy the coordinates from another point.
 * @param point a pointer to the Point that has the coordinates of interest
 */
inline void Point::copyCoords(Point* point) {
  _xyz[0] = point->getX();
  _xyz[1] = point->getY();
  _xyz[2] = point->getZ();
}


#endif /* POINT_H_ */
