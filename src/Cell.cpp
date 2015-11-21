#include "Cell.h"


int Cell::_n = 0;

static int auto_id = 10000;


/**
 * @brief Returns an auto-generated unique Cell ID.
 * @details This method is intended as a utility method for users writing
 *          OpenMOC input files. The method makes use of a static Cell
 *          ID which is incremented each time the method is called to enable
 *          unique generation of monotonically increasing IDs. The method's
 *          first ID begins at 10000. Hence, user-defined Cell IDs greater
 *          than or equal to 10000 are prohibited.
 */
int cell_id() {
  int id = auto_id;
  auto_id++;
  return id;
}


/**
 * @brief Resets the auto-generated unique Cell ID counter to 10000.
 */
void reset_cell_id() {
  auto_id = 10000;
}


/**
 * @brief Constructor sets the unique and user-specifed IDs for this Cell.
 * @param id the user-specified optional Cell ID
 * @param name the user-specified optional Cell name
 */
Cell::Cell(int id, const char* name) {

  /* If the user did not define an optional ID, create one */
  if (id == 0)
    _id = cell_id();

  /* Use the user-defined ID */
  else
    _id = id;

  _uid = _n;
  _n++;

  _name = NULL;
  setName(name);

  _cell_type = UNFILLED;
  _fill = NULL;
  _rotated = false;

  _num_rings = 0;
  _num_sectors = 0;

  /* Set a default bounding box around the Cell */
  _min_x = -std::numeric_limits<double>::infinity();
  _max_x = std::numeric_limits<double>::infinity();
  _min_y = -std::numeric_limits<double>::infinity();
  _max_y = std::numeric_limits<double>::infinity();
  _min_z = -std::numeric_limits<double>::infinity();
  _max_z = std::numeric_limits<double>::infinity();

  /* Set the default boundaries to be REFLECTIVE */
  _min_x_bc = REFLECTIVE;
  _max_x_bc = REFLECTIVE;
  _min_y_bc = REFLECTIVE;
  _max_y_bc = REFLECTIVE;
}


/**
 * @brief Destructor clears vector of Surface pointers bounding the Cell.
 */
Cell::~Cell() {

  std::map<int, surface_halfspace*>::iterator iter;
  for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
    delete iter->second;
  _surfaces.clear();

  if (_name != NULL)
    delete [] _name;
}


/**
 * @brief Return the Cell's unique ID.
 * @return the Cell's unique ID
 */
int Cell::getUid() const {
  return _uid;
}


/**
 * @brief Return the Cell's user-specified ID.
 * @return the Cell's user-specified ID
 */
int Cell::getId() const {
  return _id;
}


/**
 * @brief Return the user-defined name of the Cell
 * @return the Cell name
 */
char* Cell::getName() const {
  return _name;
}


/**
 * @brief Return the Cell type (FILL or MATERIAL).
 * @return the Cell type
 */
cellType Cell::getType() const {
  return _cell_type;
}


/**
 * @brief Return a pointer to the Material filling this Cell.
 * @return the Material fill pointer
 */
Material* Cell::getFillMaterial() {
  if (_cell_type == FILL)
    log_printf(ERROR, "Unable to get Material fill from Cell ID=%d", _id);

  return (Material*)_fill;
}


/**
 * @brief Return a pointer to the Material filling this Cell.
 * @return the Material fill pointer
 */
Universe* Cell::getFillUniverse() {
  if (_cell_type == MATERIAL)
    log_printf(ERROR, "Unable to get Universe fill from Cell ID=%d", _id);

  return (Universe*)_fill;
}


/**
 * @brief Return a boolean indicating whether or not the Cell has been rotated.
 * @return whether the Cell has been rotated
 */
bool Cell::isRotated() {
  return _rotated;
}


/**
 * @brief Return the number of rings in the Cell.
 * @return the number of rings
 */
int Cell::getNumRings() {
  return _num_rings;
}


/**
 * @brief Return the number of sectors in the Cell.
 * @return the number of sectors
 */
int Cell::getNumSectors() {
  return _num_sectors;
}


/**
 * @brief Return the minimum reachable x-coordinate in the Cell.
 * @return the minimum x-coordinate
 */
double Cell::getMinX() {
  findBoundingBox();
  return _min_x;
}


/**
 * @brief Return the maximum reachable x-coordinate in the Cell.
 * @return the maximum x-coordinate
 */
double Cell::getMaxX() {
  findBoundingBox();
  return _max_x;
}


/**
 * @brief Return the minimum reachable y-coordinate in the Cell.
 * @return the minimum y-coordinate
 */
double Cell::getMinY() {
  findBoundingBox();
  return _min_y;
}


/**
 * @brief Return the maximum reachable y-coordinate in the Cell.
 * @return the maximum y-coordinate
 */
double Cell::getMaxY() {
  findBoundingBox();
  return _max_y;
}


/**
 * @brief Return the minimum reachable z-coordinate in the Cell.
 * @return the minimum z-coordinate
 */
double Cell::getMinZ() {
  findBoundingBox();
  return _min_z;
}


/**
 * @brief Return the maximum reachable z-coordinate in the Cell.
 * @return the maximum z-coordinate
 */
double Cell::getMaxZ() {
  findBoundingBox();
  return _max_z;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE) at
 *        the minimum reachable x-coordinate in the Cell.
 * @return the boundary condition at the minimum x-coordinate
 */
boundaryType Cell::getMinXBoundaryType() {
  findBoundingBox();
  return _min_x_bc;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE) at
 *        the maximum reachable x-coordinate in the Cell.
 * @return the boundary condition at the maximum x-coordinate
 */
boundaryType Cell::getMaxXBoundaryType() {
  findBoundingBox();
  return _max_x_bc;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE) at
 *        the minimum reachable y-coordinate in the Cell.
 * @return the boundary condition at the minimum y-coordinate
 */
boundaryType Cell::getMinYBoundaryType() {
  findBoundingBox();
  return _min_y_bc;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE) at
 *        the maximum reachable y-coordinate in the Cell.
 * @return the boundary condition at the maximum y-coordinate
 */
boundaryType Cell::getMaxYBoundaryType() {
  findBoundingBox();
  return _max_y_bc;
}


/**
 * @brief Return the number of Surfaces in the Cell.
 * @return the number of Surfaces
 */
int Cell::getNumSurfaces() const {
  return _surfaces.size();
}


/**
 * @brief Return the std::map of Surface pointers and halfspaces (+/-1) for all
 *        surfaces bounding the Cell.
 * @return std::map of Surface pointers and halfspaces
 */
std::map<int, surface_halfspace*> Cell::getSurfaces() const {
  return _surfaces;
}


/**
 * @brief Return the std::vector of neighbor Cells to this Cell.
 * @return std::vector of neighbor Cell pointers
 */
std::vector<Cell*> Cell::getNeighbors() const {
  return _neighbors;
}


/**
 * @brief Returns the std::map of Cell IDs and Cell pointers within any
 *        nested Universes filling this Cell.
 * @return std::map of Cell IDs and pointers
 */
std::map<int, Cell*> Cell::getAllCells() {

  std::map<int, Cell*> cells;

  if (_cell_type == FILL && _fill != NULL) {
    std::map<int, Cell*> nested_cells;
    Universe* univ_fill = static_cast<Universe*>(_fill);

    if (univ_fill->getType() == SIMPLE)
      nested_cells = univ_fill->getAllCells();
    else
      nested_cells = static_cast<Lattice*>(univ_fill)->getAllCells();

    cells.insert(nested_cells.begin(), nested_cells.end());
  }

  return cells;
}


/**
 * @brief Returns the std::map of all nested Universe IDs and Universe pointers
          filling this Cell.
 * @return std::map of Universe IDs and pointers
 */
std::map<int, Universe*> Cell::getAllUniverses() {

  std::map<int, Universe*> universes;

  if (_cell_type == FILL && _fill != NULL) {
    Universe* univ_fill = static_cast<Universe*>(_fill);
    universes[univ_fill->getId()] = univ_fill;

    std::map<int, Universe*> nested_universes;
    if (univ_fill->getType() == SIMPLE)
      nested_universes = static_cast<Universe*>(_fill)->getAllUniverses();
    else
      nested_universes = static_cast<Lattice*>(_fill)->getAllUniverses();
    universes.insert(nested_universes.begin(), nested_universes.end());
  }

  return universes;
}


/**
 * @brief Sets the name of the Cell
 * @param name the Cell name string
 */
void Cell::setName(const char* name) {
  int length = strlen(name);

  if (_name != NULL)
    delete [] _name;

  /* Initialize a character array for the Cell's name */
  _name = new char[length+1];

  /* Copy the input character array Cell name to the class attribute name */
  for (int i=0; i <= length; i++)
    _name[i] = name[i];
}


/**
 * @brief Sets the Material filling this Cell.
 * @param fill the Material filling this Cell
 */
void Cell::setFill(Material* fill) {
  _cell_type = MATERIAL;
  _fill = fill;
}


/**
 * @brief Sets the Universe filling this Cell.
 * @param fill the Universe filling this Cell
 */
void Cell::setFill(Universe* fill) {
  _cell_type = FILL;
  _fill = fill;
}


/**
 * @brief Set the Cell's rotation angles about the x, y and z axes.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Cell's rotation angles in Python. A user must initialize a
 *          length 3 NumPy array as input to this function. This function then
 *          stores the data values in the NumPy array in the Cell's rotation
 *          array. An example of how this function might be called in Python
 *          is as follows:
 *
 * @code
 *          rotation = numpy.array([90., 0., 0.])
 *          cell = openmoc.Cell()
 *          cell.setRotation(rotation)
 * @endcode
 *
 * @param rotation the array of rotation angles
 * @param num_axes the number of axes (this must always be 3)
 */
void Cell::setRotation(double* rotation, int num_axes) {

  if (num_axes != 3)
    log_printf(ERROR, "Unable to set rotation with %d axes for Cell %d. "
               "The rotation array should be length 3.", num_axes, _id);

  for (int i=0; i < 3; i++)
    _rotation[i] = rotation[i];

  /* Compute rotation angles in x,y,z directions */
  double phi = -_rotation[0] * M_PI / 180.;
  double theta = -_rotation[1] * M_PI / 180.;
  double psi = -_rotation[2] * M_PI / 180.;

  /* Calculate rotation matrix based on angles given */
  /* Indexed by (y,x) since the universe array is indexed by (z,y,z) */
  _rotation_matrix[0][0] = cos(theta) * cos(psi);
  _rotation_matrix[1][0] = cos(theta) * sin(psi);
  _rotation_matrix[2][0] = -sin(theta);
  _rotation_matrix[0][1] = -cos(phi) * sin(psi) + 
                           sin(phi) * sin(theta) * cos(psi);
  _rotation_matrix[1][1] = cos(phi) * cos(psi) +
                           sin(phi) * sin(theta) * sin(psi);
  _rotation_matrix[2][1] = sin(phi) * cos(theta);
  _rotation_matrix[0][2] = sin(phi) * sin(psi) + \
                           cos(phi) * sin(theta) * cos(psi);
  _rotation_matrix[1][2] = -sin(phi) * cos(psi) + 
                           cos(phi) * sin(theta) * sin(psi);
  _rotation_matrix[2][2] = cos(phi) * cos(theta);

  for (int i=0; i < 3; i++) {
    for (int j=0; j < 3; j++)
      printf("i = %d , j = %d, rot mat = %f\n", i, j, _rotation_matrix[i][j]);
  }

  _rotated = true;
}


/**
 * @brief Set the Cell's number of rings.
 * @param num_rings the number of rings in this Cell
 */
void Cell::setNumRings(int num_rings) {
  if (num_rings < 0)
    log_printf(ERROR, "Unable to give %d rings to Cell %d since this is "
               "a negative number", num_rings, _id);

  _num_rings = num_rings;
}


/**
 * @brief Set the Cell's number of sectors.
 * @param num_sectors the number of sectors in this Cell
 */
void Cell::setNumSectors(int num_sectors) {
  if (num_sectors < 0)
    log_printf(ERROR, "Unable to give %d sectors to Cell %d since this is "
               "a negative number", num_sectors, _id);

  /* By default, a ring is considered to have a single sector in [0, 2*pi] */
  if (num_sectors == 1)
    _num_sectors = 0;

  else
    _num_sectors = num_sectors;
}


/**
 * @brief Insert a Surface into this Cell's container of bounding Surfaces.
 * @param halfspace the Surface halfspace (+/-1)
 * @param surface a pointer to the Surface
 */
void Cell::addSurface(int halfspace, Surface* surface) {

  if (halfspace != -1 && halfspace != +1)
    log_printf(ERROR, "Unable to add surface %d to cell %d since the halfspace"
               " %d is not -1 or 1", surface->getId(), _id, halfspace);

  surface_halfspace* new_surf_half = new surface_halfspace;
  new_surf_half->_surface = surface;
  new_surf_half->_halfspace = halfspace;

  _surfaces[surface->getId()] = new_surf_half;
}


/**
 * @brief Removes a Surface from this Cell's container of bounding Surfaces.
 * @param surface a pointer to the Surface to remove
 */
void Cell::removeSurface(Surface* surface) {

  if (_surfaces.find(surface->getId()) != _surfaces.end()) {
    delete _surfaces[surface->getId()];
    _surfaces.erase(surface->getId());
  }
}


/**
 * @brief Add a neighboring Cell to this Cell's collection of neighbors.
 * @param cell a pointer to the neighboring Cell
 */
void Cell::addNeighborCell(Cell* cell) {

  /* Add the neighbor Cell if it is not already in the collection */
  if (std::find(_neighbors.begin(), _neighbors.end(), cell) == _neighbors.end())
    _neighbors.push_back(cell);
}


/**
 * @brief Finds and stores a bounding box for the entire geometry.
 */
void Cell::findBoundingBox() {

  /* Set a default bounding box around the Cell */
  _min_x = -std::numeric_limits<double>::infinity();
  _max_x = std::numeric_limits<double>::infinity();
  _min_y = -std::numeric_limits<double>::infinity();
  _max_y = std::numeric_limits<double>::infinity();
  _min_z = -std::numeric_limits<double>::infinity();
  _max_z = std::numeric_limits<double>::infinity();

  /* Loop over all Surfaces inside the Cell */
  std::map<int, surface_halfspace*>::iterator iter;
  Surface* surface;
  int halfspace;
  double min_x, max_x, min_y, max_y, min_z, max_z;

  for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {

    surface = iter->second->_surface;
    halfspace = iter->second->_halfspace;

    max_x = surface->getMaxX(halfspace);
    max_y = surface->getMaxY(halfspace);
    max_z = surface->getMaxZ(halfspace);

    min_x = surface->getMinX(halfspace);
    min_y = surface->getMinY(halfspace);
    min_z = surface->getMinZ(halfspace);

    if (max_x != std::numeric_limits<double>::infinity() && max_x < _max_x) {
      _max_x = max_x;
      _max_x_bc = surface->getBoundaryType();
    }
    if (max_y != std::numeric_limits<double>::infinity() && max_y < _max_y) {
      _max_y = max_y;
      _max_y_bc = surface->getBoundaryType();
    }
    if (max_z != std::numeric_limits<double>::infinity() && max_z < _max_z) {
      _max_z = max_z;
    }

    if (min_x != -std::numeric_limits<double>::infinity() && min_x > _min_x) {
      _min_x = min_x;
      _min_x_bc = surface->getBoundaryType();
    }
    if (min_y != -std::numeric_limits<double>::infinity() && min_y > _min_y) {
      _min_y = min_y;
      _min_y_bc = surface->getBoundaryType();
    }
    if (min_z != -std::numeric_limits<double>::infinity() && min_z > _min_z) {
      _min_z = min_z;
    }
  }

  /* If we could not find a bounds for any dimension, readjust
   * it to +/- infinity */
  if (_max_x == -std::numeric_limits<double>::infinity())
    _max_x = std::numeric_limits<double>::infinity();
  if (_max_y == -std::numeric_limits<double>::infinity())
    _max_y = std::numeric_limits<double>::infinity();
  if (_max_z == -std::numeric_limits<double>::infinity())
    _max_z = std::numeric_limits<double>::infinity();

  if (_min_x == std::numeric_limits<double>::infinity())
    _min_x = -std::numeric_limits<double>::infinity();
  if (_min_y == std::numeric_limits<double>::infinity())
    _min_y = -std::numeric_limits<double>::infinity();
  if (_min_z == std::numeric_limits<double>::infinity())
    _min_z = -std::numeric_limits<double>::infinity();
}


/**
 * @brief Determines whether a Point is contained inside a Cell.
 * @details Queries each Surface inside the Cell to determine if the Point
 *          is on the same side of the Surface. This point is only inside
 *          the Cell if it is on the same side of every Surface in the Cell.
 * @param point a pointer to a Point
 */
bool Cell::containsPoint(Point* point) {

  /* Loop over all Surfaces inside the Cell */
  std::map<int, surface_halfspace*>::iterator iter;

  for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {

    /* Return false if the Point is not in the correct Surface halfspace */
    if (iter->second->_surface->evaluate(point) * iter->second->_halfspace
        < -ON_SURFACE_THRESH)
      return false;
  }

  /* Return true if the Point is in the correct halfspace for each Surface */
  return true;
}


/**
 * @brief Determines whether a Point is contained inside a Cell.
 * @details Queries each Surface inside the Cell to determine if the Point
 *          is on the same side of the Surface. This Point is only inside
 *          the Cell if it is on the same side of every Surface in the Cell.
 * @param coords a pointer to a localcoord
 */
bool Cell::containsCoords(LocalCoords* coords) {
  return containsPoint(coords->getPoint());
}


/**
 * @brief Computes the minimum distance to a Surface from a Point with a given
 *        trajectory at a certain angle.
 * @details If the trajectory will not intersect any of the Surfaces in the
 *          Cell returns INFINITY.
 * @param point the Point of interest
 * @param angle the angle of the trajectory (in radians from \f$[0,2\pi]\f$)
 */
double Cell::minSurfaceDist(Point* point, double angle) {

  double curr_dist;
  double min_dist = INFINITY;

  std::map<int, surface_halfspace*>::iterator iter;

  /* Loop over all of the Cell's Surfaces */
  for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {

    /* Find the minimum distance from this surface to this Point */
    curr_dist = iter->second->_surface->getMinDistance(point, angle);

    /* If the distance to Cell is less than current min distance, update */
    if (curr_dist < min_dist)
      min_dist = curr_dist;
  }

  return min_dist;
}


/**
 * @brief Create a duplicate of the Cell.
 * @return a pointer to the clone
 */
Cell* Cell::clone() {

  /* Construct new Cell */
  Cell* new_cell = new Cell();
  new_cell->setName(_name);
  new_cell->setNumRings(_num_rings);
  new_cell->setNumSectors(_num_sectors);

  if (_cell_type == MATERIAL)
    new_cell->setFill((Material*)_fill);
  else
    new_cell->setFill((Universe*)_fill);

  /* Loop over all of this Cell's Surfaces and add them to the clone */
  std::map<int, surface_halfspace*>::iterator iter;

  for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
    new_cell->addSurface(iter->second->_halfspace, iter->second->_surface);

  return new_cell;
}


/**
 * @brief Subdivides the Cell into clones for fuel pin angular sectors.
 * @param subcells an empty vector to store all subcells
 */
void Cell::sectorize(std::vector<Cell*>* subcells) {

  /* If the user didn't request any sectors, don't make any */
  if (_num_sectors == 0)
    return;

  double azim_angle;
  double delta_azim = 2. * M_PI / _num_sectors;
  double A, B;

  /* A container for each of the bounding planes for the sector Cells */
  std::vector<Plane*> planes;

  log_printf(DEBUG, "Sectorizing Cell %d with %d sectors",_id, _num_sectors);

  /* Create each of the bounding planes for the sector Cells */
  for (int i=0; i < _num_sectors; i++) {

    /* Figure out the angle for this plane */
    azim_angle = i * delta_azim;

    /* Instantiate the plane */
    A = cos(azim_angle);
    B = sin(azim_angle);
    Plane* plane = new Plane(A, B, 0., 0.);
    planes.push_back(plane);

    log_printf(DEBUG, "Created sector Plane id = %d, angle = %f, A = %f, "
               "B = %f", i, azim_angle, A, B);
  }

  /* Create sectors using disjoint halfspaces of pairing Planes */
  for (int i=0; i < _num_sectors; i++) {

    /* Create new Cell clone for this sector Cell */
    Cell* sector = clone();

    sector->setNumSectors(0);
    sector->setNumRings(0);

    log_printf(DEBUG, "Creating a new sector Cell %d for Cell %d",
               sector->getId(), _id);

    /* Add new bounding planar Surfaces to the clone */
    sector->addSurface(+1, planes.at(i));

    if (_num_sectors != 2) {
      if (i+1 < _num_sectors)
        sector->addSurface(-1, planes.at(i+1));
      else
        sector->addSurface(-1, planes.at(0));
    }

    /* Store the clone in the parent Cell's container of sector Cells */
    subcells->push_back(sector);
  }
}


/**
 * @brief Subdivides the Cell into clones for fuel pin rings.
 * @param subcells an empty vector to store all subcells
 */
void Cell::ringify(std::vector<Cell*>* subcells) {

  /* If the user didn't request any rings, don't make any */
  if (_num_rings == 0)
        return;

  int num_zcylinders = 0;
  ZCylinder* zcylinder1 = NULL;
  ZCylinder* zcylinder2 = NULL;
  double radius1 = 0;
  double radius2 = 0;
  double x1 = 0.;
  double y1 = 0.;
  double x2 = 0.;
  double y2 = 0.;
  int halfspace1 = 0;
  int halfspace2 = 0;
  std::vector<ZCylinder*> zcylinders;
  std::vector<Cell*> rings;

  /* See if the Cell contains 1 or 2 ZCYLINDER Surfaces */
  std::map<int, surface_halfspace*>::iterator iter1;
  for (iter1=_surfaces.begin(); iter1 != _surfaces.end(); ++iter1) {

    /* Determine if any of the Surfaces is a ZCylinder */
    if (iter1->second->_surface->getSurfaceType() == ZCYLINDER) {
      int halfspace = iter1->second->_halfspace;
      ZCylinder* zcylinder = static_cast<ZCylinder*>(iter1->second->_surface);

      /* Outermost bounding ZCylinder */
      if (halfspace == -1) {
        halfspace1 = halfspace;
        zcylinder1 = zcylinder;
        radius1 = zcylinder1->getRadius();
        x1 = zcylinder1->getX0();
        y1 = zcylinder1->getY0();
      }

      /* Innermost bounding zcylinder */
      else if (halfspace == +1) {
        halfspace2 = halfspace;
        zcylinder2 = zcylinder;
        radius2 = zcylinder2->getRadius();
        x2 = zcylinder2->getX0();
        y2 = zcylinder2->getY0();
      }

      num_zcylinders++;
    }
  }

  /* Error checking */
  if (num_zcylinders == 0)
    log_printf(ERROR, "Unable to ringify Cell %d since it does not "
              "contain any ZCYLINDER type Surface(s)", _id);

  if (num_zcylinders > 2)
    log_printf(NORMAL, "Unable to ringify Cell %d since it "
               "contains more than 2 ZCYLINDER Surfaces", _id);

  if (x1 != x2 && num_zcylinders == 2)
    log_printf(ERROR, "Unable to ringify Cell %d since it contains "
               "ZCylinder %d centered at x=%f and ZCylinder %d at x=%f. "
               "Both ZCylinders must have the same center.",
               _id, zcylinder1->getId(), x1, zcylinder2->getId(), x2);

  if (y1 != y2 && num_zcylinders == 2)
    log_printf(ERROR, "Unable to ringify Cell %d since it contains "
               "ZCylinder %d centered at y=%f and ZCylinder %d at y=%f. "
               "Both ZCylinders must have the same center.",
               _id, zcylinder1->getId(), y1, zcylinder2->getId(), y2);

  if (zcylinder1 == NULL && zcylinder2 != NULL)
    log_printf(ERROR, "Unable to ringify Cell %d since it only contains "
               "the positive halfpsace of ZCylinder %d. Rings can only be "
               "created for Cells on the interior (negative halfspace) "
               "of a ZCYLINDER Surface.", _id, zcylinder2->getId());

  if (radius1 <= radius2)
    log_printf(ERROR, "Unable to ringify Cell %d since it contains 2 "
               "disjoint ZCYLINDER Surfaces: halfspace %d for ZCylinder %d "
               "and halfspace %d for ZCylinder %d. Switch the signs of "
               "the 2 halfspaces for each Surface.", _id, halfspace1,
               zcylinder1->getId(), halfspace2, zcylinder2->getId());

  /* Compute the area to fill with each equal volume ring */
  double area = M_PI * fabs(radius1*radius1 - radius2*radius2) / _num_rings;

  /* Generate successively smaller ZCylinder Surfaces */
  for (int i=0; i < _num_rings-1; i++) {
    radius2 = sqrt(radius1*radius1 - (area / M_PI));
    ZCylinder* zcylinder = new ZCylinder(x1, y1, radius1);
    zcylinders.push_back(zcylinder);
    radius1 = radius2;
  }

  /* Store smallest, innermost ZCylinder */
  ZCylinder* zcylinder = new ZCylinder(x1, y1, radius1);
  zcylinders.push_back(zcylinder);

  /* Loop over ZCylinders and create a new Cell clone for each ring */
  std::vector<ZCylinder*>::iterator iter2;
  std::vector<Cell*>::iterator iter3;

  for (iter2 = zcylinders.begin(); iter2 != zcylinders.end(); ++iter2) {

    /* Create ZCylinders for each of the sectorized Cells */
    if (subcells->size() != 0) {
      for (iter3 = subcells->begin(); iter3 != subcells->end(); ++iter3) {
        log_printf(DEBUG, "Creating a new ring in sector Cell %d",
                   (*iter3)->getId());

        /* Create a new Cell clone */
        Cell* ring = (*iter3)->clone();
        ring->setNumSectors(0);
        ring->setNumRings(0);

        /* Add new bounding ZCylinder surfaces to the clone */
        ring->addSurface(-1, (*iter2));

        /* Look ahead and check if we have an inner ZCylinder to add */
        if (iter2+1 == zcylinders.end()) {
          rings.push_back(ring);
          continue;
        }
        else
          ring->addSurface(+1, *(iter2+1));

        /* Store the clone in the parent Cell's container of ring Cells */
        rings.push_back(ring);
      }
    }

    /* Create ZCylinders for this un-sectorized Cell */
    else {
      log_printf(DEBUG, "Creating new ring in un-sectorized Cell %d",_id);

      /* Create a new Cell clone */
      Cell* ring = clone();
      ring->setNumSectors(0);
      ring->setNumRings(0);

      /* Add new bounding ZCylinder Surfaces to the clone */
      ring->addSurface(-1, (*iter2));

      /* Look ahead and check if we have an inner ZCylinder to add */
      if (iter2+1 == zcylinders.end()) {
        rings.push_back(ring);
        break;
      }
      else
        ring->addSurface(+1, *(iter2+1));

      /* Store the clone in the parent Cell's container of ring Cells */
      rings.push_back(ring);
    }
  }

  /* Store all of the rings in the parent Cell's subcells container */
  subcells->clear();
  subcells->insert(subcells->end(), rings.begin(), rings.end());
}


/**
 * @brief Subdivides a Cell into rings and sectors.
 * @details This method uses the Cell's clone method to produce a vector of
 *          this Cell's subdivided ring and sector Cells.
 * @return a vector of Cell pointers to the new subdivided Cells
 */
void Cell::subdivideCell() {

  /** A container of all Cell clones created for rings and sectors */
  std::vector<Cell*>* subcells = new std::vector<Cell*>();

  sectorize(subcells);
  ringify(subcells);

  /* Put any ring / sector subcells in a new Universe fill */
  if (subcells->size() != 0) {

    /* Create a new Universe to contain all of the subcells */
    Universe* new_fill = new Universe();

    /* Add each subcell to the new Universe fill */
    std::vector<Cell*>::iterator iter;
    for (iter = subcells->begin(); iter != subcells->end(); ++iter)
      new_fill->addCell(*iter);

    /* Set the new Universe as the fill for this Cell */
    setFill(new_fill);
  }
}


/**
 * @brief Build a collection of neighboring Cells for optimized ray tracing.
 */
void Cell::buildNeighbors() {

  Surface* surface;
  int halfspace;

  /* Add this Cell to all of the Surfaces in this Cell */
  std::map<int, surface_halfspace*>::iterator iter;
  for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {
    surface = iter->second->_surface;
    halfspace = iter->second->_halfspace;
    surface->addNeighborCell(halfspace, this);
  }

  /* Make recursive call to the Cell's fill Universe */
  if (_cell_type == FILL) {
    Universe* fill = static_cast<Universe*>(_fill);
    if (fill->getType() == SIMPLE)
      static_cast<Universe*>(fill)->buildNeighbors();
    else
      static_cast<Lattice*>(fill)->buildNeighbors();
  }
}


/**
 * @brief Convert this Cell's attributes to a string format.
 * @return a character array of this Cell's attributes
 */
std::string Cell::toString() {

  std::stringstream string;

  string << "Cell ID = " << _id
         << ", name = " << _name
         << ", # rings = " << _num_rings
         << ", # sectors = " << _num_sectors;

  if (_cell_type == FILL) {
    string << ", type = FILL, "
           << ", fill id = " << static_cast<Universe*>(_fill)->getId();
  }
  else if(_cell_type == MATERIAL) {
    string << ", type = MATERIAL"
           << ", fill id = " << static_cast<Material*>(_fill)->getId();
  }
  else
    string << ", type = UNFILLED";

  string << ", # surfaces = " << getNumSurfaces();

  /** Add string data for the Surfaces in this Cell */
  std::map<int, surface_halfspace*>::iterator iter;
  string << ", Surfaces: ";
  for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
    string <<  iter->second->_surface->toString() << ", ";

  return string.str();
}


/**
 * @brief Prints a string representation of all of the Cell's attributes
 *        to the console.
 */
void Cell::printString() {
  log_printf(NORMAL, toString().c_str());
}
