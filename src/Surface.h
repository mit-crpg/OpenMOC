/**
 * @file Surface.h
 * @details The Surface class and subclasses.
 * @date January 9, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef SURFACE_H_
#define SURFACE_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "constants.h"
#include "LocalCoords.h"
#include "boundary_type.h"
#include <limits>
#include <map>
#include <vector>
#include <algorithm>
#endif


/* Forward declarations to resolve circular dependencies */
class LocalCoords;
class Cell;


int surface_id();
void reset_surface_id();
void maximize_surface_id(int surface_id);


/**
 * @enum surfaceType
 * @brief The types of surfaces supported by OpenMOC.
 */
enum surfaceType {
  /** A general plane */
  PLANE,

  /** A cylinder with axis parallel to the z-axis */
  ZCYLINDER,

  /** A plane perpendicular to the x-axis */
  XPLANE,

  /** A plane perpendicular to the y-axis */
  YPLANE,

  /** A plane perpendicular to the z-axis */
  ZPLANE,

  /** A generalized quadratic surface */
  QUADRATIC
};



/**
 * @class Surface Surface.h "src/Surface.h"
 * @brief Represents a general Surface in 3D.
 * @details The Surface class and its subclasses are used to define the
 *          geometry for an OpenMOC simulation using a constructive solid
 *          geometry (CSG) formalism. Surfaces are used during ray tracing
 *          of charateristic tracks across the geometry.
 */
class Surface {

protected:

  /** A static counter for the number of Surfaces in a simulation */
  static int _n;

  /** A monotonically increasing unique ID for each Surface created */
  int _uid;

  /** A user-defined id for each Surface created */
  int _id;

  /** A user-defined name for the Surface */
  char* _name;

  /** The type of Surface (ie, XPLANE, ZCYLINDER, etc) */
  surfaceType _surface_type;

  /** The type of boundary condition to be used for this Surface
   *  (ie, VACUUM or REFLECTIVE) */
  boundaryType _boundary_type;

  /* Vector of neighboring Cells */
  std::map<int, std::vector<Cell*>* > _neighbors;

public:
  Surface(const int id=0, const char* name="");
  virtual ~Surface();

  int getUid() const;
  int getId() const;
  char* getName() const;
  surfaceType getSurfaceType();
  boundaryType getBoundaryType();


  /**
   * @brief Returns the minimum coordinate in the axis direction of the
   *        space defined by halfspace and this surface.
   * @param axis The axis of interest (0 = x, 1 = y, 2 = z)
   * @param halfspace the halfspace to consider
   * @return the minimum coordinate in the axis direction
   */
  double getMin(int axis, int halfspace);


  /**
   * @brief Returns the maximum coordinate in the axis direction of the
   *        space defined by halfspace and this surface.
   * @param axis The axis of interest (0 = x, 1 = y, 2 = z)
   * @param halfspace the halfspace to consider
   * @return the maximum coordinate in the axis direction
   */
  double getMax(int axis, int halfspace);

  /**
   * @brief Returns the minimum x value for one of this Surface's halfspaces.
   * @param halfspace the halfspace of the Surface to consider
   * @return the minimum x value
   */
  virtual double getMinX(int halfspace) = 0;

  /**
   * @brief Returns the maximum x value for one of this Surface's halfspaces.
   * @param halfspace the halfspace of the Surface to consider
   * @return the maximum x value
   */
  virtual double getMaxX(int halfspace) = 0;

  /**
   * @brief Returns the minimum y value for one of this Surface's halfspaces.
   * @param halfspace the halfspace of the Surface to consider
   * @return the minimum y value
   */
  virtual double getMinY(int halfspace) = 0;

  /**
   * @brief Returns the maximum y value for one of this Surface's halfspaces.
   * @param halfspace the halfspace of the Surface to consider
   * @return the maximum y value
   */
  virtual double getMaxY(int halfspace) = 0;

  /**
   * @brief Returns the minimum z value for one of this Surface's halfspaces.
   * @param halfspace the halfspace of the Surface to consider
   * @return the minimum z value
   */
  virtual double getMinZ(int halfspace) = 0;

  /**
   * @brief Returns the maximum z value for one of this Surface's halfspaces.
   * @param halfspace the halfspace of the Surface to consider
   * @return the maximum z value
   */
  virtual double getMaxZ(int halfspace) = 0;

  void setName(const char* name);
  void setBoundaryType(const boundaryType boundary_type);
  void addNeighborCell(int halfspace, Cell* cell);

  /**
   * @brief Evaluate a Point using the Surface's potential equation.
   * @details This method returns the values \f$ f(x,y) \f$ for the potential
   *          function \f$f\f$ representing this Surface.
   * @param point a pointer to the Point of interest
   * @return the value of Point in the Plane's potential equation.
   */
  virtual double evaluate(const Point* point) const = 0;

  /**
   * @brief Finds the intersection Point with this Surface from a given
   *        Point and trajectory defined by an angle.
   * @param point pointer to the Point of interest
   * @param azim the azimuthal angle (in radians)
   * @param polar the polar angle (in radians)
   * @param points array of Points to store the intersection locations
   * @return the number of intersection Points (0 or 1)
   */
  virtual int intersection(Point* point, double azim, double polar,
                           Point* points) = 0;

  bool isPointOnSurface(Point* point);
  bool isCoordOnSurface(LocalCoords* coord);
  double getMinDistance(Point* point, double azim, double polar);
  double getMinDistance(LocalCoords* coord);

  /**
   * @brief Converts this Surface's attributes to a character array.
   * @details The character array returned conatins the type of Surface (ie,
   *          PLANE) and the coefficients in the potential equation.
   * @return a character array of this Surface's attributes
   */
  virtual std::string toString() = 0;

  void printString();
};


/**
 * @class Plane Surface.h "src/Surface.h"
 * @brief Represents a Plane perpendicular to the xy-plane.
 */
class Plane: public Surface {

protected:

  /** The coefficient for the linear term in x */
  double _A;

  /** The coefficient for the linear term in y */
  double _B;

  /** The coefficient for the linear term in z */
  double _C;

  /** The constant offset */
  double _D;

  /** The Plane is a friend of class Surface */
  friend class Surface;

  /** The Plane is a friend of class Zcylinder */
  friend class ZCylinder;

public:

  Plane(const double A, const double B, const double C, const double D,
        const int id=0, const char* name="");

  double getMinX(int halfspace);
  double getMaxX(int halfspace);
  double getMinY(int halfspace);
  double getMaxY(int halfspace);
  double getMinZ(int halfspace);
  double getMaxZ(int halfspace);
  double getA();
  double getB();
  double getC();
  double getD();

  double evaluate(const Point* point) const;
  int intersection(Point* point, double azim, double polar, Point* points);

  std::string toString();
};


/**
 * @class XPlane Surface.h "src/Surface.h"
 * @brief Represents a Plane perpendicular to the x-axis.
 */
class XPlane: public Plane {

private:

  /** The location of the XPlane along the x-axis */
  double _x;

public:
  XPlane(const double x, const int id=0, const char* name="");

  void setX(const double x);

  double getX();
  double getMinX(int halfspace);
  double getMaxX(int halfspace);

  std::string toString();
};


/**
 * @class YPlane Surface.h "src/Surface.h"
 * @brief Represents a Plane perpendicular to the y-axis.
 */
class YPlane: public Plane {

private:

  /** The location of the YPlane along the y-axis */
  double _y;

public:
  YPlane(const double y, const int id=0, const char* name="");

  void setY(const double y);

  double getY();
  double getMinY(int halfspace);
  double getMaxY(int halfspace);

  std::string toString();
};


/**
 * @class ZPlane Surface.h "src/Surface.h"
 * @brief Represents a Plane perpendicular to the z-axis.
 */
class ZPlane: public Plane {

private:

  /** The location of the ZPlane along the z-axis */
  double _z;

public:
  ZPlane(const double z, const int id=0, const char* name="");

  void setZ(const double z);

  double getZ();
  double getMinZ(int halfspace);
  double getMaxZ(int halfspace);

  std::string toString();
};


/**
 * @class ZCylinder Surface.h "src/Surface.h"
 * @brief Represents a Cylinder with axis parallel to the z-axis.
 */
class ZCylinder: public Surface {

private:

  /** A point for the ZCylinder's center */
  Point _center;

  /** The ZCylinder's radius */
  double _radius;

  /** The coefficient of the x-squared term */
  double _A;

  /** The coefficient of the y-squared term */
  double _B;

  /** The coefficient of the linear term in x */
  double _C;

  /** The coefficient of the linear term in y */
  double _D;

  /** The constant offset */
  double _E;

  /** The ZCylinder is a friend of the Surface class */
  friend class Surface;

  /** The ZCylinder is a friend of the Plane class */
  friend class Plane;

public:
  ZCylinder(const double x, const double y, const double radius,
            const int id=0, const char* name="");

  double getX0();
  double getY0();
  double getRadius();
  double getMinX(int halfspace);
  double getMaxX(int halfspace);
  double getMinY(int halfspace);
  double getMaxY(int halfspace);
  double getMinZ(int halfspace);
  double getMaxZ(int halfspace);

  double evaluate(const Point* point) const;
  int intersection(Point* point, double azim, double polar, Point* points);

  std::string toString();
};


/**
 * @brief Finds the minimum distance to a Surface.
 * @details Finds the minimum distance to a Surface from a Point with a
 *          given trajectory defined by an azim/polar to this Surface. If the
 *          trajectory will not intersect the Surface, returns INFINITY.
 * @param point a pointer to the Point of interest
 * @param azim the azimuthal angle defining the trajectory in radians
 * @param polar the polar angle defining the trajectory in radians
 * @return the minimum distance to the Surface
 */
inline double Surface::getMinDistance(Point* point, double azim, double polar) {

  /* Point array for intersections with this Surface */
  Point intersections[2];

  /* Find the intersection Point(s) */
  int num_inters = intersection(point, azim, polar, intersections);
  double distance = INFINITY;

  /* If there is one intersection Point */
  if (num_inters == 1)
    distance = intersections[0].distanceToPoint(point);

  /* If there are two intersection Points */
  else if (num_inters == 2) {
    double dist1 = intersections[0].distanceToPoint(point);
    double dist2 = intersections[1].distanceToPoint(point);

    /* Determine which intersection Point is nearest */
    if (dist1 < dist2)
      distance = dist1;
    else
      distance = dist2;
  }

  return distance;
}


/**
 * @brief Evaluate a Point using the Plane's quadratic Surface equation.
 * @param point a pointer to the Point of interest
 * @return the value of Point in the Plane's quadratic equation
 */
inline double Plane::evaluate(const Point* point) const {
  double x = point->getX();
  double y = point->getY();
  double z = point->getZ();

  return (_A * x + _B * y + _C * z + _D);
}


/**
 * @brief Return the radius of the ZCylinder.
 * @return the radius of the ZCylinder
 */
inline double ZCylinder::getRadius() {
  return this->_radius;
}


/**
 * @brief Evaluate a Point using the ZCylinder's quadratic Surface equation.
 * @param point a pointer to the Point of interest
 * @return the value of Point in the equation
 */
inline double ZCylinder::evaluate(const Point* point) const {
  double x = point->getX();
  double y = point->getY();
  return (_A * x * x + _B * y * y + _C * x + _D * y + _E);
}


#endif /* SURFACE_H_ */
