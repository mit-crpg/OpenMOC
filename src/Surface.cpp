#include "Surface.h"

int Surface::_n = 0;

static int auto_id = DEFAULT_INIT_ID;

/**
 * @brief Returns an auto-generated unique surface ID.
 * @details This method is intended as a utility mehtod for user's writing
 *          OpenMOC input files. The method makes use of a static surface
 *          ID which is incremented each time the method is called to enable
 *          unique generation of monotonically increasing IDs. The method's
 *          first ID begins at 10000. Hence, user-defined surface IDs greater
 *          than or equal to 10000 are prohibited.
 */
int surface_id() {
  int id = auto_id;
  auto_id++;
  return id;
}


/**
 * @brief Resets the auto-generated unique Surface ID counter to 10000.
 */
void reset_surface_id() {
  auto_id = DEFAULT_INIT_ID;
}


/**
 * @brief Maximize the auto-generated unique Surface ID counter.
 * @details This method updates the auto-generated unique Surface ID
 *          counter if the input parameter is greater than the present
 *          value. This is useful for the OpenCG compatibility module
 *          to ensure that the auto-generated Surface IDs do not
 *          collide with those created in OpenCG.
 * @param surface_id the id assigned to the auto-generated counter
 */
void maximize_surface_id(int surface_id) {
  if (surface_id > auto_id)
    auto_id = surface_id;
}


/**
 * @brief Constructor assigns unique ID and user-defined ID for a Surface.
 * @details Assigns a default boundary condition for this Surface to
 *          BOUNDARY_NONE.
 * @param id an optional user-defined Surface ID
 * @param name an optional user-defined Surface name
 */
Surface::Surface(const int id, const char* name) {

  /* If the user did not define an optional ID, create one */
  if (id == 0)
    _id = surface_id();

  /* Use the user-defined ID */
  else
    _id = id;

  _uid = _n;
  _n++;

  _name = NULL;
  setName(name);

  _boundary_type = BOUNDARY_NONE;

  /* Initialize empty vectors of neighbor Cells for each halfspace */
  _neighbors[-1] = new std::vector<Cell*>();
  _neighbors[+1] = new std::vector<Cell*>();
}


/**
 * @brief Destructor.
 */
Surface::~Surface() {
  if (_name != NULL)
    delete [] _name;

  if (!_neighbors.empty()) {
    _neighbors[-1]->clear();
    _neighbors[+1]->clear();
    delete _neighbors[-1];
    delete _neighbors[+1];
    _neighbors.clear();
  }
}


/**
 * @brief Return the Surface's unique ID.
 * @return the Surface's unique ID
 */
int Surface::getUid() const {
  return _uid;
}


/**
 * @brief Return the Surface's user-defined ID.
 * @return the Surface's user-defined ID
 */
int Surface::getId() const {
  return _id;
}


/**
 * @brief Return the user-defined name of the Surface
 * @return the Surface name
 */
char* Surface::getName() const {
  return _name;
}


/**
 * @brief Return the type of Surface (ie, XPLANE, ZCYLINDER, etc).
 * @return the Surface type
 */
surfaceType Surface::getSurfaceType() {
  return _surface_type;
}


/**
 * @brief Returns the type of boundary conditions for this Surface (REFLECTIVE,
 *        VACUUM or BOUNDARY_NONE)
 * @return the type of boundary condition type for this Surface
 */
boundaryType Surface::getBoundaryType() {
  return _boundary_type;
}


/**
 * @brief Sets the name of the Surface
 * @param name the Surface name string
 */
void Surface::setName(const char* name) {
  int length = strlen(name);

  if (_name != NULL)
    delete [] _name;

  /* Initialize a character array for the Surface's name */
  _name = new char[length+1];

  /* Copy the input character array Surface name to the class attribute name */
  for (int i=0; i <= length; i++)
    _name[i] = name[i];
}


/**
 * @brief Sets the boundary condition type (ie, VACUUM or REFLECTIVE) for this
 *        Surface.
 * @param boundary_type the boundary condition type for this Surface
 */
void Surface::setBoundaryType(boundaryType boundary_type) {
  _boundary_type = boundary_type;
}


/**
 * @brief Adds a neighbor Cell to this Surface's collection of neighbors.
 * @param halfspace the +/-1 halfspace for the neighboring Cell
 * @param cell a pointer to the neighboring Cell
 */
void Surface::addNeighborCell(int halfspace, Cell* cell) {

  if (halfspace != -1 && halfspace != +1)
    log_printf(ERROR, "Unable to add neighbor Cell %d to Surface %d since the "
               "halfspace %d is not -1 or 1", cell->getId(), _id, halfspace);

  /* Get pointer to vector of neighbor Cells for this halfspace */
  std::vector<Cell*>* neighbors = _neighbors[halfspace];

  /* Add the neighbor Cell if the collection does not already contain it*/
  if (std::find(neighbors->begin(), neighbors->end(), cell) == neighbors->end())
    neighbors->push_back(cell);

  /* Update Cells with the neighbor Cells on the opposite Surface halfspace */
  std::vector<Cell*>::iterator iter1;
  std::vector<Cell*>::iterator iter2;
  for (iter1 = _neighbors[-1]->begin();
       iter1 != _neighbors[-1]->end(); ++iter1) {
    for (iter2 = _neighbors[1]->begin();
         iter2 != _neighbors[1]->end(); ++iter2)
      (*iter1)->addNeighborCell(*iter2);
  }

  for (iter1 = _neighbors[1]->begin();
       iter1 != _neighbors[1]->end(); ++iter1) {
    for (iter2 = _neighbors[-1]->begin();
         iter2 != _neighbors[-1]->end(); ++iter2)
      (*iter1)->addNeighborCell(*iter2);
  }
}


/**
 * @brief Return true or false if a Point is on or off of a Surface.
 * @param point pointer to the Point of interest
 * @return on (true) or off (false) the Surface
 */
bool Surface::isPointOnSurface(Point* point) {

  /* Uses a threshold to determine whether the point is on the Surface */
  if (fabs(evaluate(point)) < ON_SURFACE_THRESH)
    return true;
  else
    return false;
}


/**
 * @brief Return true or false if a LocalCoord is on or off of a Surface.
 * @param coord pointer to the LocalCoord of interest
 * @return on (true) or off (false) the Surface
 */
bool Surface::isCoordOnSurface(LocalCoords* coord) {
  return isPointOnSurface(coord->getPoint());
}


/**
 * @brief Finds the minimum distance to a Surface.
 * @details Finds the miniumum distance to a Surface from a LocalCoords
 *          with a trajectory defined by an angle to this Surface. If the
 *          trajectory will not intersect the Surface, returns INFINITY.
 * @param coords a pointer to a localcoords object
 * @return the minimum distance to the Surface
 */
double Surface::getMinDistance(LocalCoords* coords) {

  Point* point = coords->getPoint();
  double phi = coords->getPhi();

  /* Point array for intersections with this Surface */
  Point intersections[2];

  /* Find the intersection Point(s) */
  int num_inters = this->intersection(point, phi, intersections);
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
 * @brief Prints a string representation of all of the Surface's objects to
 *        the console.
 */
void Surface::printString() {
  log_printf(RESULT, toString().c_str());
}


/**
 * @brief Constructor.
 * @param A the first coefficient in \f$ A * x + B * y + C * z + D = 0 \f$
 * @param B the second coefficient in \f$ A * x + B * y + C * z + D = 0 \f$
 * @param C the third coefficient in \f$ A * x + B * y + C * z + D = 0 \f$
 * @param D the fourth coefficient in \f$ A * x + B * y + C * z + D = 0 \f$
 * @param id the optional Surface ID
 * @param name the optional name of the Surface
 */
Plane::Plane(const double A, const double B,
             const double C, const double D, const int id, const char* name):
  Surface(id, name) {

  _surface_type = PLANE;
  _A = A;
  _B = B;
  _C = C;
  _D = D;
}


/**
 * @brief Returns the minimum x value of -INFINITY
 * @param halfspace the halfspace of the Surface to consider
 * @return the minimum x value of -INFINITY
 */
double Plane::getMinX(int halfspace) {
  return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum x value of INFINITY.
 * @param halfspace the halfspace of the Surface to consider
 * @return the maximum x value of INFINITY
 */
double Plane::getMaxX(int halfspace) {
  return std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the minimum y value of -INFINITY
 * @param halfspace the halfspace of the Surface to consider
 * @return the minimum y value of -INFINITY
 */
double Plane::getMinY(int halfspace) {
  return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum y value of INFINITY.
 * @param halfspace the halfspace of the Surface to consider
 * @return the maximum y value of INFINITY
 */
double Plane::getMaxY(int halfspace) {
  return std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the minimum z value of -INFINITY
 * @param halfspace the halfspace of the Surface to consider
 * @return the minimum z value of -INFINITY
 */
double Plane::getMinZ(int halfspace) {
  return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum z value of INFINITY.
 * @param halfspace the halfspace of the Surface to consider
 * @return the maximum z value of INFINITY
 */
double Plane::getMaxZ(int halfspace) {
  return std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the A coefficient multiplying x in the surface equation
 * @return the value for the A coefficient
 */
double Plane::getA() {
  return _A;
}


/**
 * @brief Returns the B coefficient multiplying y in the surface equation
 * @return the value for the B coefficient
 */
double Plane::getB() {
  return _B;
}


/**
 * @brief Returns the C coefficient multiplying z in the surface equation
 * @return the value for the C coefficient
 */
double Plane::getC() {
  return _C;
}


/**
 * @brief Returns the D constant coefficient
 * @return the value for the D coefficient
 */
double Plane::getD() {
  return _D;
}


/**
* @brief Finds the intersection Point with this Plane from a given Point and
*        trajectory defined by an angle.
* @param point pointer to the Point of interest
* @param angle the angle defining the trajectory in radians
* @param points pointer to a Point to store the intersection Point
* @return the number of intersection Points (0 or 1)
*/
inline int Plane::intersection(Point* point, double angle, Point* points) {

  double x0 = point->getX();
  double y0 = point->getY();
  double z0 = point->getZ();
  double l;

  int num = 0;                /* number of intersections */
  double xcurr, ycurr, zcurr; /* coordinates of current intersection point */
  double mx = cos(angle);
  double my = sin(angle);

  /* The track and plane are parallel */
  if ((fabs(mx) < 1.e-10 && fabs(_A) > 1.e-10) ||
      (fabs(my) < 1.e-10 && fabs(_B) > 1.e-10))
    return 0;

  /* The track is not parallel to the plane */
  else{

    l = - (_A*x0 + _B*y0 + _C*z0 + _D) /
      (_A * mx + _B * my);
    xcurr = x0 + l * mx;
    ycurr = y0 + l * my;
    zcurr = z0;

    if (l > 0.0) {
      points[num].setCoords(xcurr, ycurr, zcurr);
      num++;
    }
  }

  return num;
}


/**
 * @brief Converts this Plane's attributes to a character array.
 * @details The character array returned conatins the type of Plane (ie,
 *          PLANE) and the A, B, and C coefficients in the
 *          quadratic Surface equation.
 * @return a character array of this Plane's attributes
 */
std::string Plane::toString() {

  std::stringstream string;

  string << "Surface ID = " << _id
         << ", name = " << _name
         << ", type = PLANE "
         << ", A = " << _A << ", B = " << _B << ", C = " << _C
         << ", D = " << _D;

  return string.str();
}


/**
 * @brief Constructor for a Plane perpendicular to the x-axis.
 * @param x the location of the Plane along the x-axis
 * @param id the optional Surface id
 * @param name the optional name of the XPlane
 */
XPlane::XPlane(const double x, const int id, const char* name):
  Plane(1, 0, 0, -x, id, name) {

  _surface_type = XPLANE;
  _x = x;
}


/**
 * @brief Set the location of this XPlane on the x-axis.
 * @param x the location of the XPlane on the x-axis
 */
void XPlane::setX(const double x) {
  _x = x;
  _D = -x;
}


/**
 * @brief Returns the location of the XPlane on the x-axis.
 * @return the location of the XPlane on the x-axis
 */
double XPlane::getX() {
  return _x;
}


/**
 * @brief Returns the minimum x value for one of this XPlane's halfspaces.
 * @param halfspace the halfspace of the XPlane to consider
 * @return the minimum x value
 */
double XPlane::getMinX(int halfspace) {
  if (halfspace == +1)
    return _x;
  else
    return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum x value for one of this XPlane's halfspaces.
 * @param halfspace the halfspace of the XPlane to consider
 * @return the maximum x value
 */
double XPlane::getMaxX(int halfspace) {
  if (halfspace == -1)
    return _x;
  else
    return std::numeric_limits<double>::infinity();
}


/**
 * @brief Converts this XPlane's attributes to a character array.
 * @details The character array returned conatins the type of Plane (ie,
 *          XPLANE) and the A, B, C, and D coefficients in the
 *          quadratic Surface equation and the location of the Plane on
 *          the x-axis.
 * @return a character array of this XPlane's attributes
 */
std::string XPlane::toString() {

  std::stringstream string;

  string << "Surface ID = " << _id
         << ", name = " << _name
         << ", type = XPLANE "
         << ", A = " << _A << ", B = " << _B
         << ", C = " << _C << ", D = " << _D
         << ", x = " << _x;

  return string.str();
}


/**
 * @brief Constructor for a Plane perpendicular to the y-axis.
 * @param y the location of the Plane along the y-axis
 * @param id the optional Surface id
 * @param name the optional Surface name
 */
YPlane::YPlane(const double y, const int id, const char* name):
  Plane(0, 1, 0, -y, id, name) {

  _surface_type = YPLANE;
  _y = y;
}


/**
 * @brief Set the location of this YPlane on the y-axis.
 * @param y the location of the YPlane on the y-axis
 */
void YPlane::setY(const double y) {
  _y = y;
  _D = -y;
}


/**
 * @brief Returns the location of the YPlane on the y-axis.
 * @return the location of the YPlane on the y-axis
 */
double YPlane::getY() {
  return _y;
}


/**
 * @brief Returns the minimum y value for one of this YPlane's halfspaces.
 * @param halfspace the halfspace of the YPlane to consider
 * @return the minimum y value
 */
double YPlane::getMinY(int halfspace) {
  if (halfspace == +1)
    return _y;
  else
    return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum y value for one of this YPlane's halfspaces.
 * @param halfspace the halfspace of the YPlane to consider
 * @return the maximum y value
 */
double YPlane::getMaxY(int halfspace) {
  if (halfspace == -1)
    return _y;
  else
    return std::numeric_limits<double>::infinity();
}


/**
 * @brief Converts this yplane's attributes to a character array
 * @details The character array returned conatins the type of Plane (ie,
 *          YPLANE) and the A, B, C, and D coefficients in the quadratic
 *          Surface equation and the location of the Plane on the y-axis.
 * @return a character array of this YPlane's attributes
 */
std::string YPlane::toString() {

  std::stringstream string;

  string << "Surface ID = " << _id
         << ", name = " << _name
         << ", type = YPLANE "
         << ", A = " << _A << ", B = " << _B
         << ", C = " << _C << ", D = " << _D
         << ", y = " << _y;

  return string.str();
}


/**
 * @brief Constructor for a Plane perpendicular to the z-axis.
 * @param z the location of the Plane along the z-axis
 * @param id the optional Surface ID
 * @param name the optional Surface name
 */
ZPlane::ZPlane(const double z, const int id, const char* name):
  Plane(0, 0, 1, -z, id, name) {

  _surface_type = ZPLANE;
  _z = z;
}


/**
 * @brief Set the location of this ZPlane on the z-axis.
 * @param z the location of the ZPlane on the z-axis
 */
void ZPlane::setZ(const double z) {
  _z = z;
  _D = -z;
}


/**
 * @brief Returns the location of the ZPlane on the z-axis.
 * @return the location of the ZPlane on the z-axis
 */
double ZPlane::getZ() {
  return _z;
}


/**
 * @brief Returns the minimum z value for one of this ZPlane's halfspaces.
 * @param halfspace the halfspace of the ZPlane to consider
 * @return the minimum z value
 */
double ZPlane::getMinZ(int halfspace) {
  if (halfspace == +1)
    return _z;
  else
    return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum z value for one of this ZPlane's halfspaces.
 * @param halfspace the halfspace of the ZPlane to consider
 * @return the maximum z value
 */
double ZPlane::getMaxZ(int halfspace) {
  if (halfspace == -1)
    return _z;
  else
    return std::numeric_limits<double>::infinity();
}


/**
 * @brief Converts this ZPlane's attributes to a character array.
 * @details The character array returned conatins the type of Plane (ie,
 *          ZPLANE) and the A, B, C, and D coefficients in the
 *          quadratic Surface equation and the location of the Plane along
 *          the z-axis.
 * @return a character array of this ZPlane's attributes
 */
std::string ZPlane::toString() {

  std::stringstream string;

  string << "Surface ID = " << _id
         << ", name = " << _name
         << ", type = ZPLANE "
         << ", A = " << _A << ", B = " << _B
         << ", C = " << _C << ", D = " << _D
         << ", z = " << _z;

  return string.str();
}


/**
 * @brief constructor.
 * @param x the x-coordinte of the ZCylinder center
 * @param y the y-coordinate of the ZCylinder center
 * @param radius the radius of the ZCylinder
 * @param id the optional Surface ID
 * @param name the optional Surface name
 */
ZCylinder::ZCylinder(const double x, const double y,
                     const double radius, const int id, const char* name):
    Surface(id, name) {

  _surface_type = ZCYLINDER;
  _A = 1.;
  _B = 1.;
  _C = -2.*x;
  _D = -2.*y;
  _E = x*x + y*y - radius*radius;
  _radius = radius;
  _center.setX(x);
  _center.setY(y);
}


/**
 * @brief Return the x-coordinate of the ZCylinder's center Point.
 * @return the x-coordinate of the ZCylinder center
 */
double ZCylinder::getX0() {
  return _center.getX();
}


/**
 * @brief Return the y-coordinate of the ZCylinder's center Point.
 * @return the y-coordinate of the ZCylinder center
 */
double ZCylinder::getY0() {
  return _center.getY();
}


/**
 * @brief Returns the minimum x value for one of this ZCylinder's halfspaces.
 * @param halfspace the halfspace of the ZCylinder to consider
 * @return the minimum x value
 */
double ZCylinder::getMinX(int halfspace) {
  if (halfspace == -1)
    return _center.getX() - _radius;
  else
    return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum x value for one of this ZCylinder's halfspaces.
 * @param halfspace the halfspace of the ZCylinder to consider
 * @return the maximum x value
 */
double ZCylinder::getMaxX(int halfspace) {
  if (halfspace == -1)
    return _center.getX() + _radius;
  else
    return std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the minimum y value for one of this ZCylinder's halfspaces.
 * @param halfspace the halfspace of the ZCylinder to consider
 * @return the minimum y value
 */
double ZCylinder::getMinY(int halfspace) {
  if (halfspace == -1)
    return _center.getY() - _radius;
  else
    return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum y value for one of this ZCylinder's halfspaces.
 * @param halfspace the halfspace of the ZCylinder to consider
 * @return the maximum y value
 */
double ZCylinder::getMaxY(int halfspace) {
  if (halfspace == -1)
    return _center.getY() + _radius;
  else
    return std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the minimum z value of -INFINITY.
 * @param halfspace the halfspace of the ZCylinder to consider
 * @return the minimum z value of -INFINITY
 */
double ZCylinder::getMinZ(int halfspace) {
  return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum z value of INFINITY.
 * @param halfspace the halfspace of the ZCylinder to consider
 * @return the maximum z value of INFINITY
 */
double ZCylinder::getMaxZ(int halfspace) {
  return std::numeric_limits<double>::infinity();
}


/**
 * @brief Finds the intersection Point with this zcylinder from a given Point
 *        and trajectory defined by an angle (0, 1, or 2 points).
 * @param point pointer to the Point of interest
 * @param angle the angle defining the trajectory in radians
 * @param points pointer to a an array of Points to store intersection Points
 * @param polar the polar angle defining the trajectory in radians
 * @return the number of intersection Points (0 or 1)
 */
int ZCylinder::intersection(Point* point, double angle, Point* points) {

  double x0 = point->getX();
  double y0 = point->getY();
  double z0 = point->getZ();
  double xcurr, ycurr, zcurr;
  int num = 0;                        /* Number of intersection Points */
  double a, b, c, q, discr;
  /*std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n";
  std::cout << "            ZCylinder calculation \n";
  std::cout << "x0, y0, z0 " << x0 << " " << y0 << " " << z0 << std::endl;
  std::cout << "angle " << angle << std::endl;
  std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n";
*/
  /* If the track is vertical in y */
/*  if ((fabs(angle - M_PI_2)) < TINY_MOVE) {

    std::cout << "!!!!!!track is vertical " << std::endl;

    * Solve for where the line x = x0 and the Surface F(x,y) intersect
     * Find the y where F(x0, y) = 0
     * Substitute x0 into F(x,y) and rearrange to put in
     * the form of the quadratic formula: ay^2 + by + c = 0 
     * This is simplified for a z-cylinder with a vertical axis *
    a = 1.;
    b = _D;
    c = _A * x0 * x0 + _C * x0 + _E;

    discr = b*b - 4*c;

    * There are no intersections *
    if (discr < 0)
      return 0;

    * There is one intersection (ie on the Surface) *
    else if (discr == 0) {
      xcurr = x0;
      ycurr = -b / 2;
      zcurr = z0;
      points[num].setCoords(xcurr, ycurr, zcurr);

      * Check that point is in same direction as angle *
      if (angle < M_PI && ycurr > y0)
        num++;
      else if (angle > M_PI && ycurr < y0)
        num++;

      return num;
    }

    * There are two intersections *
    else {
      xcurr = x0;
      zcurr = z0;
      
      // solve quadratic equation for y in a way that reduces error
      if (b > 0)
        ycurr = 2*c / (-b - sqrt(discr));
      if (b < 0)
        ycurr = -b/2 + sqrt(discr) / 2;

      points[num].setCoords(xcurr, ycurr, zcurr);
      if (angle < M_PI && ycurr > y0)
        num++;
      else if (angle > M_PI && ycurr < y0)
        num++;

      xcurr = x0;
      zcurr = z0;

      // solve quadratic equation for y in a way that reduces error
      if (b > 0)
        ycurr = -b/2 - sqrt(discr) / 2;
      if (b < 0)
        ycurr = 2*c / (-b + sqrt(discr));

      points[num].setCoords(xcurr, ycurr, zcurr);
      if (angle < M_PI && ycurr > y0)
        num++;
      else if (angle > M_PI && ycurr < y0)
        num++;

      return num;
    }
  }
i*/
  // if track is near vertical
  if (abs(tan(angle)) >= 1) {
    //std::cout << "track is not vertical " << std::endl;

    /* Solve for where the line y-y0 = m*(x-x0) and the Surface F(x,y)
     * intersect. Find the (x,y) where F(x, y0 + m*(x-x0)) = 0
     * Substitute the point-slope formula for y into F(x,y) and
     * rearrange to put in the form of the quadratic formula:
     * ax^2 + bx + c = 0
     */
    double m = tan(angle);

    q = x0 - 1/m * y0;
    a = 1 + 1/m *1/ m;
    b = 2 * 1/m * q + _D + _C * m;
    c = q * q + _C * q + _E;

    discr = b*b - 4*a*c;

// check this    /* Boolean value describing whether the track is traveling to the right */
    bool right = angle < M_PI / 2. || angle > 3. * M_PI / 2.;

    /* There are no intersections */
    if (discr < 0)
      return 0;

    /* There is one intersection (ie on the Surface) */
    else if (discr == 0) {

      /* Determine the point of intersection */
      ycurr = -b / (2*a);
      xcurr = x0 + 1/m * points[num].getY() - 1/m*y0;
      zcurr = z0;
      points[num].setCoords(xcurr, ycurr, zcurr);

      /* Increase the number of intersections if the intersection is in the
       * direction of the track is heading */
      if (right && xcurr > x0)
        num++;
      else if (!right && xcurr < x0)
        num++;

      return num;
    }

    /* There are two intersections */
    else {

      /* Determine the point of intersection */
      if (b > 0)
        ycurr = 2*c / (-b - sqrt(discr));
      if (b < 0)
        ycurr = -b/(2*a) + sqrt(discr) / (2*a);

      // this could lead to lots of error when the track is near-vertical
      xcurr = x0 + 1/m*(ycurr -y0);

      zcurr = z0;
      points[num].setCoords(xcurr, ycurr, zcurr);

      /* Increase the number of intersections if the intersection is in the
       * direction of the track is heading */
      if (right && xcurr > x0)
        num++;
      else if (!right && xcurr < x0)
        num++;

      /* Determine the point of intersection */
      if (b > 0)
        ycurr = -b/(2*a) - sqrt(discr) / (2*a);
      if (b < 0)
        ycurr = 2*c / (-b + sqrt(discr));

      // this could lead to lots of error when track is near-vertical
      xcurr = x0 + 1/m*ycurr -1/m*y0;
      zcurr = z0;
      points[num].setCoords(xcurr, ycurr, zcurr);
      //std::cout << "intersection point " << xcurr << " " << ycurr << "\n";

      /* Increase the number of intersections if the intersection is in the
       * direction of the track is heading */
      if (right && xcurr > x0)
        num++;
      else if (!right && xcurr < x0)
        num++;

      return num;
    }
  }

  /* If the track is near horizontal */
  else {
    //std::cout << "track is not vertical " << std::endl;

    /* Solve for where the line y-y0 = m*(x-x0) and the Surface F(x,y)
     * intersect. Find the (x,y) where F(x, y0 + m*(x-x0)) = 0
     * Substitute the point-slope formula for y into F(x,y) and
     * rearrange to put in the form of the quadratic formula:
     * ax^2 + bx + c = 0
     */
    double m = tan(angle);

    q = y0 - m * x0;
    a = 1 + m * m;
    b = 2 * m * q + _C + _D * m;
    c = q * q + _D * q + _E;

    discr = b*b - 4*a*c;

    /* Boolean value describing whether the track is traveling to the right */
    bool right = angle < M_PI / 2. || angle > 3. * M_PI / 2.;

    /* There are no intersections */
    if (discr < 0)
      return 0;

    /* There is one intersection (ie on the Surface) */
    else if (discr == 0) {

      /* Determine the point of intersection */
      xcurr = -b / (2*a);
      ycurr = y0 + m * points[num].getX() - m*x0;
      zcurr = z0;
      points[num].setCoords(xcurr, ycurr, zcurr);

      /* Increase the number of intersections if the intersection is in the
       * direction of the track is heading */
      if (right && xcurr > x0)
        num++;
      else if (!right && xcurr < x0)
        num++;

      return num;
    }

    /* There are two intersections */
    else {

      /* Determine the point of intersection */
      if (b > 0)
        xcurr = 2*c / (-b - sqrt(discr));
      if (b < 0)
        xcurr = -b/(2*a) + sqrt(discr) / (2*a);

      // this could lead to lots of error when the track is near-vertical
      ycurr = y0 + m*(xcurr -x0);

      zcurr = z0;
      points[num].setCoords(xcurr, ycurr, zcurr);

      /* Increase the number of intersections if the intersection is in the
       * direction of the track is heading */
      if (right && xcurr > x0)
        num++;
      else if (!right && xcurr < x0)
        num++;

      /* Determine the point of intersection */
      if (b > 0)
        xcurr = -b/(2*a) - sqrt(discr) / (2*a);
      if (b < 0)
        xcurr = 2*c / (-b + sqrt(discr));

      // this could lead to lots of error when track is near-vertical
      ycurr = y0 + m*xcurr -m*x0;
/*      std::cout << "--r2: m*xcurr   = " << m*xcurr << std::endl;
      std::cout << "      m*x0      = " << m*x0 << std::endl;
      std::cout << " m*(xcurr - x0) = " << m*(xcurr - x0) << std::endl;
      std::cout << "      y0        = " << y0 << std::endl;
      std::cout << "    y0 + m*()   = " << ycurr << std::endl;
      std::cout << "   distribute m = " << y0 + m*xcurr - m*x0
        << std::endl << std::endl;
*/
      zcurr = z0;
      points[num].setCoords(xcurr, ycurr, zcurr);
      //std::cout << "intersection point " << xcurr << " " << ycurr << "\n";

      /* Increase the number of intersections if the intersection is in the
       * direction of the track is heading */
      if (right && xcurr > x0)
        num++;
      else if (!right && xcurr < x0)
        num++;

      return num;
    }
  }
}


/**
 * @brief Converts this ZCylinder's attributes to a character array.
 * @details The character array returned conatins the type of Plane (ie,
 *          ZCYLINDER) and the A, B, C, D and E coefficients in the
 *          quadratic Surface equation.
 * @return a character array of this ZCylinder's attributes
 */
std::string ZCylinder::toString() {

  std::stringstream string;

  string << "Surface ID = " << _id
         << ", name " << _name
         << ", type = ZCYLINDER "
         << ", A = " << _A << ", B = " << _B
         << ", C = " << _C << ", D = " << _D << ", E = " << _E
         << ", x0 = " << _center.getX()
         << ", y0 = " << _center.getY()
         << ", radius = " << _radius;

    return string.str();
}

double Surface::kahanSum(double i, double j, double k) {
  double sum = 0;
  double c = 0;
  
  double y = i - c;
  double t = sum + y;
  c = (t - sum) - y;
  sum = t;
  
  y = j - c;
  t = sum + y;
  c = (t - sum) - y;
  sum = t;

  y = k - c;
  t = sum + y;
  c = (t - sum) - y;
  sum = t;

  return sum;
}
