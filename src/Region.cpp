#include "Region.h"

/**
 * @brief Destructor clears vector of the nodes within the Region.
 */
Region::~Region() {
  _nodes.clear();
}


/**
 * @brief Add a node to the Region.
 * @details NOTE: This method deep copies the Region and stores
 *          the copy. Any changes made to the Region will not be
 *          reflected in the Region copy stored by the Region.
 * @param node a Region node to add to this Region
 */
void Region::addNode(Region* node) {
  _nodes.push_back(node->clone());
}


/**
 * @brief Return a vector of all of the Region's nodes.
 * @details This method will return a list of the Region's nodes
 *          when called from Python.
 * @returns a vector of the Region's nodes
 */
std::vector<Region*> Region::getNodes() {
  return _nodes;
}


/**
 * @brief Extracts a map of all Halfspaces contained in the Region.
 * @details This method recurses through all of the Region's nodes
 *          and collects the Halfspaces into a map indexed by 
 *          Surface ID. This method will return a dictionary of the
 *          Halfspaces indexed by Surface ID called in Python.
 * @returns a map of all Halfspaces in the Region
 */
std::map<int, Halfspace*> Region::getAllSurfaces() {

  std::map<int, Halfspace*> all_surfaces;
  std::map<int, Halfspace*> node_surfaces;
  std::vector<Region*>::iterator iter;
  /* Recursively collect all Halfspace's from this Regions nodes */
  for (iter = _nodes.begin(); iter != _nodes.end(); iter++) {
    node_surfaces = (*iter)->getAllSurfaces();
    all_surfaces.insert(node_surfaces.begin(), node_surfaces.end());
  }
  return all_surfaces;
}


/**
 * @brief Return the type of Region (ie, UNION, INTERSECTION, etc).
 * @return the Region type
 */
regionType Region::getRegionType() {
  return _region_type;
}


/**
 * @brief Return the minimum reachable x-coordinate in the Region.
 * @return the minimum x-coordinate
 */
double Region::getMinX() {

  /* Set a default min-x */
  double min_x = -std::numeric_limits<double>::infinity();

  /* Loop over all nodes in the Region */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); ++iter)
    min_x = std::max(min_x, (*iter)->getMinX());
  return min_x;
}


/**
 * @brief Return the maximum reachable x-coordinate in the Region.
 * @return the maximum x-coordinate
 */
double Region::getMaxX() {

  /* Set a default max-x */
  double max_x = std::numeric_limits<double>::infinity();

  /* Loop over all nodes in the Region */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); ++iter)
    max_x = std::min(max_x, (*iter)->getMaxX());
  return max_x;
}
 /**
 * @brief Return the minimum reachable y-coordinate in the Region.
 * @return the minimum y-coordinate
 */
double Region::getMinY() {

  /* Set a default min-y */
  double min_y = -std::numeric_limits<double>::infinity();

  /* Loop over all nodes in the Region */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); ++iter)
    min_y = std::max(min_y, (*iter)->getMinY());
  return min_y;
}
 /**
 * @brief Return the maximum reachable y-coordinate in the Region.
 * @return the maximum y-coordinate
 */
double Region::getMaxY() {

  /* Set a default max-y */
  double max_y = std::numeric_limits<double>::infinity();

  /* Loop over all nodes in the Region */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); ++iter)
    max_y = std::min(max_y, (*iter)->getMaxY());
  return max_y;
}
 /**
 * @brief Return the minimum reachable z-coordinate in the Region.
 * @return the minimum z-coordinate
 */
double Region::getMinZ() {

  /* Set a default min-z */
  double min_z = -std::numeric_limits<double>::infinity();

  /* Loop over all nodes in the Region */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); ++iter)
    min_z = std::max(min_z, (*iter)->getMinZ());
  return min_z;
}


/**
 * @brief Return the maximum reachable z-coordinate in the Region.
 * @return the maximum z-coordinate
 */
double Region::getMaxZ() {

  /* Set a default max-z */
  double max_z = std::numeric_limits<double>::infinity();

  /* Loop over all nodes in the Region */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); ++iter)
    max_z = std::min(max_z, (*iter)->getMaxZ());
  return max_z;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE)
 *        at the minimum reachable x-coordinate in the Region.
 * @return the boundary condition at the minimum x-coordinate
 */
boundaryType Region::getMinXBoundaryType() {

  /* Set a default min-x and boundary type*/
  double min_x = -std::numeric_limits<double>::infinity();
  boundaryType bc = BOUNDARY_NONE;

  /* Loop over all nodes inside the Region */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); ++iter) {
    if (min_x < (*iter)->getMinX()) {
      min_x = (*iter)->getMinX();
      bc = (*iter)->getMinXBoundaryType();
    }
  }
  return bc;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE)
 *        at the maximum reachable x-coordinate in the Region.
 * @return the boundary condition at the maximum x-coordinate
 */
boundaryType Region::getMaxXBoundaryType() {

  /* Set a default max-x and boundary type*/
  double max_x = std::numeric_limits<double>::infinity();

  boundaryType bc = BOUNDARY_NONE;
  /* Loop over all nodes inside the Region */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); ++iter) {
    if (max_x > (*iter)->getMaxX()) {
      max_x = (*iter)->getMaxX();
      bc = (*iter)->getMaxXBoundaryType();
    }
  }
  return bc;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE)
 *        at the minimum reachable y-coordinate in the Region.
 * @return the boundary condition at the minimum y-coordinate
 */
boundaryType Region::getMinYBoundaryType() {

  /* Set a default min-y and boundary type*/
  double min_y = -std::numeric_limits<double>::infinity();

  boundaryType bc = BOUNDARY_NONE;
  /* Loop over all nodes inside the Region */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); ++iter) {
    if (min_y < (*iter)->getMinY()) {
      min_y = (*iter)->getMinY();
      bc = (*iter)->getMinYBoundaryType();
    }
  }
  return bc;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE)
 *        at the maximum reachable y-coordinate in the Region.
 * @return the boundary condition at the maximum y-coordinate
 */
boundaryType Region::getMaxYBoundaryType() {

  /* Set a default max-y and boundary type*/
  double max_y = std::numeric_limits<double>::infinity();
  boundaryType bc = BOUNDARY_NONE;

  /* Loop over all nodes inside the Region */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); ++iter) {
    if (max_y > (*iter)->getMaxY()) {
      max_y = (*iter)->getMaxY();
      bc = (*iter)->getMaxYBoundaryType();
    }
  }
  return bc;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE)
 *        at the minimum reachable z-coordinate in the Region.
 * @return the boundary condition at the minimum z-coordinate
 */
boundaryType Region::getMinZBoundaryType() {

  /* Set a default min-z and boundary type*/
  double min_z = -std::numeric_limits<double>::infinity();
  boundaryType bc = BOUNDARY_NONE;

  /* Loop over all nodes inside the Region */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); ++iter) {
    if (min_z < (*iter)->getMinZ()) {
      min_z = (*iter)->getMinZ();
      bc = (*iter)->getMinZBoundaryType();
    }
  }
  return bc;
}


 /**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE)
 *        at the maximum reachable z-coordinate in the Region.
 * @return the boundary condition at the maximum z-coordinate
 */
boundaryType Region::getMaxZBoundaryType() {

  /* Set a default max-z and boundary type*/
  double max_z = std::numeric_limits<double>::infinity();
  boundaryType bc = BOUNDARY_NONE;

  /* Loop over all nodes inside the Region */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); ++iter) {
    if (max_z > (*iter)->getMaxZ()) {
      max_z = (*iter)->getMaxZ();
      bc = (*iter)->getMaxZBoundaryType();
    }
  }
  return bc;
}


/**
 * @brief Computes the minimum distance to a Surface in the Region from
 *        a point with a given trajectory at a certain angle stored in a
 *        LocalCoords object.
 * @details If the trajectory will not intersect any of the Surfaces in the
 *          Region returns INFINITY.
 * @param coords a pointer to a localcoords
 */
double Region::minSurfaceDist(LocalCoords* coords) {

  double curr_dist;
  double min_dist = INFINITY;

  /* Find the minimum distance to one of the Region's nodes */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); ++iter) {
    curr_dist = (*iter)->minSurfaceDist(coords);

    /* If the distance to Cell is less than current min distance, update */
    if (curr_dist < min_dist)
      min_dist = curr_dist;
  }

  return min_dist;
}


/**
 * @brief Computes the minimum distance to a Surface in the Region from
 *        a point with a given trajectory at a certain angle stored in a
 *        LocalCoords object.
 * @details If the trajectory will not intersect any of the Surfaces in the
 *          Region returns INFINITY.
 * @param coords a pointer to a localcoords
 */
double Region::minSurfaceDist(Point* point, double azim, double polar) {

  double curr_dist;
  double min_dist = INFINITY;

  /* Find the minimum distance to one of the Region's nodes */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); ++iter) {
    curr_dist = (*iter)->minSurfaceDist(point, azim, polar);

    /* If the distance to Cell is less than current min distance, update */
    if (curr_dist < min_dist)
      min_dist = curr_dist;
  }
  return min_dist;
}


/**
 * @brief Create a duplicate of the Region.
 * @return a pointer to the clone
 */
Region* Region::clone() {

  /* Instantiate appropriate class type for the clone */
  Region* clone;
  if (dynamic_cast<Intersection*>(this))
    clone = new Intersection();
  else if (dynamic_cast<Union*>(this))
    clone = new Union();
  else if (dynamic_cast<Complement*>(this))
    clone = new Complement();

  /* Add this region's nodes to the cloned region */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); iter++)
    clone->addNode((*iter));
  return clone;
}


/**
 * @brief Constructor sets the type of Region (INTERSECTION).
 */
Intersection::Intersection() {
  _region_type = INTERSECTION;
}


/**
 * @brief Determines whether a Point is contained inside the Intersection.
 * @details Queries each of the Intersection's nodes to determine if the
 *          Point is within the Intersection. This point is only inside the
 *          Intersection if it is contained by each and every node.
 * @param point a pointer to a Point
 * @returns true if the Point is inside the Intersection; otherwise false
 */
bool Intersection::containsPoint(Point* point) {

  /* Query each of the Intersection's nodes */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); iter++) {
    if (!(*iter)->containsPoint(point))
      return false;
  }
  return true;
}


/**
 * @brief Constructor sets the type of Region (UNION).
 */
Union::Union() {
  _region_type = UNION;
}


/**
 * @brief Determines whether a Point is contained inside the Union.
 * @details Queries each of the Union's nodes to determine if the Point
 *          is within the Union. This point is only inside the
 *          Union if it is contained by at least one node.
 * @returns true if the Point is inside the Union; otherwise false
 */
bool Union::containsPoint(Point* point) {

  /* Query each of the Intersection's nodes */
  std::vector<Region*>::iterator iter;
  for (iter = _nodes.begin(); iter != _nodes.end(); iter++) {
    if ((*iter)->containsPoint(point))
      return true;
  }
  return false;
}


/**
 * @brief Constructor sets the type of Region (COMPLEMENT).
 */
Complement::Complement() {
  _region_type = COMPLEMENT;
}


/**
 * @brief Determines whether a Point is contained inside the Union.
 * @details Queries each of the Complement's nodes to determine if the Point
 *          is within the Complement. This point is only inside the
 *          Complement is not contained by the Complement's node.
 * @param point a pointer to a Point
 * @returns true if the Point is inside the Complement; otherwise false
 */
bool Complement::containsPoint(Point* point) {
  if (_nodes.size() == 0)
    return false;
  else
    return !_nodes[0]->containsPoint(point);
}


/**
 * @brief Constructor sets the type of Region (HALFSPACE).
 * @param halfspace the side of the Surface (+1 or -1)
 * @param surface a pointer to the Surface of interest
 */
Halfspace::Halfspace(int halfspace, Surface* surface) {
  if (halfspace != -1 && halfspace != +1)
    log_printf(ERROR, "Unable to create Halfspace from Surface %d since the "
	       "halfspace %d is not -1 or 1", surface->getId(), halfspace);
   _region_type = HALFSPACE;
  _surface = surface;
  _halfspace = halfspace;
}


/**
 * @brief Create a duplicate of the Halfspace.
 * @return a pointer to the clone
 */
Halfspace* Halfspace::clone() {
  Halfspace* clone = new Halfspace(_halfspace, _surface);
  return clone;
}


/**
 * @brief Return a pointer to the Halfspace's Surface.
 * @return a pointer to the Halfspace's Surface
 */
Surface* Halfspace::getSurface() {
  return _surface;
}


/**
 * @brief Return the side of the Surface for this Halfspace.
 * @param the side of the surface for this Halfspace (+1 or -1)
 */
int Halfspace::getHalfspace() {
  return _halfspace;
}


/**
 * @brief Extracts a map of this Halfspace indexed by its Surface ID.
 * @details This is a base case and a helper method for the parent class'
 *          Region::getAllSurfaces() method. This method will return a
 *          dictionary of the Halfspace indexed by Surface ID if called
 *          in Python.
 * @returns a map of this Halfspace indexed by its Surface ID
 */
std::map<int, Halfspace*> Halfspace::getAllSurfaces() {
  std::map<int, Halfspace*> all_surfaces;
  all_surfaces[_surface->getId()] = this;
  return all_surfaces;
}


/**
 * @brief Return the minimum reachable x-coordinate in the Halfspace.
 * @return the minimum x-coordinate
 */
double Halfspace::getMinX() {
  return _surface->getMinX(_halfspace);
}


/**
 * @brief Return the maximum reachable x-coordinate in the Halfspace.
 * @return the maximum x-coordinate
 */
double Halfspace::getMaxX() {
  return _surface->getMaxX(_halfspace);
}


/**
 * @brief Return the minimum reachable y-coordinate in the Halfspace.
 * @return the minimum y-coordinate
 */
double Halfspace::getMinY() {
  return _surface->getMinY(_halfspace);
}


/**
 * @brief Return the maximum reachable y-coordinate in the Halfspace.
 * @return the maximum y-coordinate
 */
double Halfspace::getMaxY() {
  return _surface->getMaxY(_halfspace);
}


/**
 * @brief Return the minimum reachable z-coordinate in the Halfspace.
 * @return the minimum z-coordinate
 */
double Halfspace::getMinZ() {
  return _surface->getMinZ(_halfspace);
}


/**
 * @brief Return the maximum reachable z-coordinate in the Halfspace.
 * @return the maximum z-coordinate
 */
double Halfspace::getMaxZ() {
  return _surface->getMaxZ(_halfspace);
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE)
 *        of the Halfspace's Surface.
 * @return the boundary condition
 */
boundaryType Halfspace::getMinXBoundaryType() {
  return _surface->getBoundaryType();
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE)
 *        of the Halfspace's Surface.
 * @return the boundary condition
 */
boundaryType Halfspace::getMaxXBoundaryType() {
  return _surface->getBoundaryType();
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE)
 *        of the Halfspace's Surface.
 * @return the boundary condition
 */
boundaryType Halfspace::getMinYBoundaryType() {
  return _surface->getBoundaryType();
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE)
 *        of the Halfspace's Surface.
 * @return the boundary condition
 */
boundaryType Halfspace::getMaxYBoundaryType() {
  return _surface->getBoundaryType();
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE)
 *        of the Halfspace's Surface.
 * @return the boundary condition
 */
boundaryType Halfspace::getMinZBoundaryType() {
  return _surface->getBoundaryType();
}


 /**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE)
 *        of the Halfspace's Surface.
 * @return the boundary condition
 */
boundaryType Halfspace::getMaxZBoundaryType() {
  return _surface->getBoundaryType();
}


/**
 * @brief Determines whether a Point is contained inside the Halfspace.
 * @details Queries whether the Point is on the same side of the
 *          Surface as the Halfspace.
 * @returns true if the Point is inside the Halfspace; otherwise false
 */
bool Halfspace::containsPoint(Point* point) {

  double evaluation = _surface->evaluate(point) * _halfspace;
  return (evaluation >= 0);
}


/**
 * @brief Computes the minimum distance to the Surface in the Halfspace from
 *        a point with a given trajectory at a certain angle stored in a
 *        LocalCoords object.
 * @details If the trajectory will not intersect the Surface in the
 *          Halfspace returns INFINITY.
 * @param coords a pointer to a localcoords
 */
double Halfspace::minSurfaceDist(Point* point, double azim, double polar) {
  return _surface->getMinDistance(point, azim, polar);
}


/**
 * @brief Computes the minimum distance to the Surface in the Halfspace from
 *        a point with a given trajectory at a certain angle stored in a
 *        LocalCoords object.
 * @details If the trajectory will not intersect the Surface in the
 *          Halfspace returns INFINITY.
 * @param coords a pointer to a localcoords
 */
double Halfspace::minSurfaceDist(LocalCoords* coords) {
  return _surface->getMinDistance(coords);
}


/**
 * @brief Constructor creates an Intersection of the Intersection of
 *        two XPlane and two YPlane objects.
 * @details This is a special subclass of the Intersection which
 *          represents the interior of a rectangular prism aligned
 *          with th z-axis.
 * @param width_x the width of the prism along the x-axis (in cm)
 * @param width_y the width of the prism along the y-axis (in cm)
 * @param origin_x the center of the prism along the x-axis (in cm)
 * @param origin_y the center of the prism along the y-axis (in cm)
 * @returns a pointer to an Intersection object
 */
RectangularPrism::RectangularPrism(double width_x, double width_y,
				   double origin_x, double origin_y):
  Intersection() {

  /* Instantiate the XPlane and YPlane objects bounding the prism */
  XPlane* min_x = new XPlane(origin_x-width_x/2.);
  XPlane* max_x = new XPlane(origin_x+width_x/2.);
  YPlane* min_y = new YPlane(origin_y-width_y/2.);
  YPlane* max_y = new YPlane(origin_y+width_y/2.);

  /* Instantiate Haflspace objects for each XPlane and YPlane */
  Halfspace* half_min_x = new Halfspace(+1, min_x);
  Halfspace* half_max_x = new Halfspace(-1, max_x);
  Halfspace* half_min_y = new Halfspace(+1, min_y);
  Halfspace* half_max_y = new Halfspace(-1, max_y);

  /* Add the Halfspace node to the Intersection */
  addNode(half_min_x);
  addNode(half_max_x);
  addNode(half_min_y);
  addNode(half_max_y);

  /* Deallocate memory for the Halfspace's since the Intersection
   * stores copies of each one */
  delete half_min_x;
  delete half_max_x;
  delete half_min_y;
  delete half_max_y;
}


/**
 * @brief Sets the boundary condition type (ie., VACUUM, REFLECTIVE, etc)
 *        to assign to each of the XPlanes and YPlanes bounding the prism.
 * @param boundary_type the boundary condition type for this Prism
 */
void RectangularPrism::setBoundaryType(boundaryType boundary_type) {

  std::map<int, Halfspace*> all_surfaces = getAllSurfaces();
  std::map<int, Halfspace*>::iterator iter;

  /* Assign the boundary to each of the bounding XPlanes and YPlanes */
  for (iter = all_surfaces.begin(); iter != all_surfaces.end(); iter++)
    iter->second->getSurface()->setBoundaryType(boundary_type);
}
