/**
 * @file Point.h
 * @brief The Point class.
 * @date January 18, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef POINT_H_
#define POINT_H_

#ifdef __cplusplus
#include "Python.h"
#include "log.h"
#include <math.h>
#include <sstream>
#endif

/**
 * @class Point Point.h "src/Point.h"
 * @brief Class to represent a 2D point in space.
 */
class Point {

private:

  /** The Point's x-coordinate */
  double _x;

  /** The Point's y-coordinate */
  double _y;

public:
  Point();
  virtual ~Point();
  void setCoords(const double x, const double y);
  double getX() const;
  double getY() const;
  void setX(const double x);
  void setY(const double y);
  double distance(const double x, const double y) const;
  double distanceToPoint(const Point* point);
  std::string toString();
};


/**
 * @brief Initializes a Point with two-dimensional coordinates.
 * @param x x-coordinate
 * @param y y-coordinate
 */
inline void Point::setCoords(const double x, const double y) {
  _x = x;
  _y = y;
}


/**
 * @brief Returns this Point's x-coordinate.
 * @return the x-coordinate
 */
inline double Point::getX() const {
  return _x;
}


/**
 * @brief Returns this Point's y-coordinate.
 * @return the y-coordinate
 */
inline double Point::getY() const {
  return _y;
}


/**
 * @brief Set the Point's x-coordinate.
 * @param x the new x-coordinate
 */
inline void Point::setX(const double x) {
  _x = x;
}


/**
 * @brief Set the Point's y-coordinate
 * @param y the new y-coordinate
 */
inline void Point::setY(const double y) {
  _y = y;
}


/**
 * @brief Compute the distance from this Point to another Point of interest.
 * @param x the x-coordinate of the Point of interest
 * @param y the y-coordinate of the Point of interest
 * @return distance to the Point of interest
 */
inline double Point::distance(const double x, const double y) const {
  double deltax = _x - x;
  double deltay = _y - y;
  return sqrt(deltax*deltax + deltay*deltay);
}


/**
 * @brief Compute the distance from this Point to another Point of interest.
 * @param point a pointer to the Point of interest
 * @return distance to the Point of interest
 */
inline double Point::distanceToPoint(const Point* point) {
  double deltax = _x - point->_x;
  double deltay = _y - point->_y;
  return sqrt(deltax*deltax + deltay*deltay);
}


#endif /* POINT_H_ */
