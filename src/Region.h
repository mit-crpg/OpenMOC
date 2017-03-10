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
#include "Point.h"
#include <limits>
#include <string>
#endif


// FIXME: Add a Halfspace class
// FIXME: Add this to SWIG
// FIXME: Incorporate this into the Cell class


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

public:
  // FIXME: Are these needed for an abstract class???
  Region();
  virtual ~Region();

  virtual std::vector<Region*> getNodes() =0;

  Intersection* getIntersection(Region* other);
  Union* getUnion(Region* other);
  Complement* getInversion();

  // FIXME: Use the notation used by the ray tracing code
  virtual bool contains(Point* point) =0;
};


/**
 * @class Intersection Intersection.h "src/Region.h"
 * @brief An intersection of two or more Regions.
 */
class Intersection : public Region {

private:

  // FIXME: Add doxygen comment
  std::vector<Region*> _nodes;

public:
  Intersection();
  virtual ~Intersection();

  void addNode(Region* node);

  std::vector<Region*> getNodes();
  
  Intersection* getIntersection(Region* other);
  bool contains(Point* point);
};


/**
 * @class Union Union.h "src/Region.h"
 * @brief A union of two or more Regions.
 */
class Union : public Region {

private:

  // FIXME: Add doxygen comment
  std::vector<Region*> _nodes;

public:
  Union();
  virtual ~Union();

  void addNode(Region* node);

  std::vector<Region*> getNodes();
  
  Union* getUnion(Region* other);
  bool contains(Point* point);
    
};



/**
 * @class Complement Complement.h "src/Region.h"
 * @brief A complement of a Region.
 */
class Complement : public Region {

private:

  // FIXME: Add doxygen comment
  Region* _node;

public:
  Complement();
  virtual ~Complement();

  // FIXME: should this retain the original addNodes() syntax???
  void setNode(Region* node);

  // FIXME: should this retain the original getNodes() syntax???
  std::vector<Region*> getNodes();
  
  bool contains(Point* point);  
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
  Halfspace(Surface* surface, int halfspace);
  virtual ~Halfspace();

  // FIXME: this may be bullshit
  std::vector<Region*> getNodes();

  Intersection* getIntersection(Region* other);
  Union* getUnion(Region* other);
  Halfspace* getInversion();
  
  bool contains(Point* point);  
};


#endif /* REGION_H_ */
