/**
 * @file Region.h
 * @brief The Region class.
 * @date March 10, 2017
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#ifndef REGION_H_
#define REGION_H_

#ifdef __cplusplus
#ifdef SWIG
#include "Python.h"
#endif
#include "Surface.h"
#include "boundary_type.h"
#include <limits>
#endif


/* Forward declarations to resolve circular dependencies */
class Intersection;
class Union;
class Complement;
class Halfspace;


/**
 * @class Region Region.h "src/Region.h"
 * @brief A region of space that can be assigned to a Cell.
 */
class Region {

protected:

  // FIXME: Add doxygen comment
  std::vector<Region*> _nodes;

public:
  virtual void addNode(Region* node);
  virtual std::vector<Region*> getNodes();
  virtual std::map<int, Halfspace*> getAllSurfaces();

  virtual double getMinX();
  virtual double getMaxX();
  virtual double getMinY();
  virtual double getMaxY();
  virtual double getMinZ();
  virtual double getMaxZ();
  virtual boundaryType getMinXBoundaryType();
  virtual boundaryType getMaxXBoundaryType();
  virtual boundaryType getMinYBoundaryType();
  virtual boundaryType getMaxYBoundaryType();

  Intersection* getIntersection(Region* other);
  Union* getUnion(Region* other);
  Complement* getInversion();

  // FIXME: Use the notation used by the ray tracing code
  virtual bool containsPoint(Point* point) =0;
  virtual double minSurfaceDist(LocalCoords* coords) =0;
  virtual Region* clone() =0;
};


/**
 * @class Intersection Intersection.h "src/Region.h"
 * @brief An intersection of two or more Regions.
 */
class Intersection : public Region {

public:
  Intersection* getIntersection(Region* other);
  bool containsPoint(Point* point);
  double minSurfaceDist(LocalCoords* coords);
  Intersection* clone();
};


/**
 * @class Union Union.h "src/Region.h"
 * @brief A union of two or more Regions.
 */
class Union : public Region {

 public:
  Union* getUnion(Region* other);
  bool containsPoint(Point* point);
  double minSurfaceDist(LocalCoords* coords);    
  Union* clone();
};


/**
 * @class Complement Complement.h "src/Region.h"
 * @brief A complement of a Region.
 */
class Complement : public Region {

 public:
  bool containsPoint(Point* point);  
  double minSurfaceDist(LocalCoords* coords);
  Complement* clone();
};



/**
 * @class Halfspace Halfspace.h "src/Region.h"
 * @brief A positive or negative halfspace Region.
 */
class Halfspace : public Region {

private:

  /** A pointer to the Surface object */
  Surface* _surface;

  /** The halfspace associated with this surface */
  int _halfspace;
  
public:

  Halfspace(int halfspace, Surface* surface);
  virtual ~Halfspace();
  Halfspace* clone();

  Surface* getSurface();
  int getHalfspace();

  // FIXME: this may be bullshit
  void addNode(Region* node);
  std::vector<Region*> getNodes();
  std::map<int, Halfspace*> getAllSurfaces();

  double getMinX();
  double getMaxX();
  double getMinY();
  double getMaxY();
  double getMinZ();
  double getMaxZ();
  boundaryType getMinXBoundaryType();
  boundaryType getMaxXBoundaryType();
  boundaryType getMinYBoundaryType();
  boundaryType getMaxYBoundaryType();

  Intersection* getIntersection(Region* other);
  Union* getUnion(Region* other);
  Halfspace* getInversion();
  
  bool containsPoint(Point* point);  
  double minSurfaceDist(LocalCoords* coords);
};


/**
 * @class Halfspace Halfspace.h "src/Region.h"
 * @brief A positive or negative halfspace Region.
 */
class RectangularPrism : public Intersection {

public:
  RectangularPrism(double width_x, double width_y,
		   double origin_x=0., double origin_y=0.);
  void setBoundaryType(boundaryType boundary_type);
};

#endif /* REGION_H_ */
