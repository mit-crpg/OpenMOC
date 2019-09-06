#include "Cell.h"
#include <set>
#include <cmath>


int Cell::_n = 0;

static int auto_id = DEFAULT_INIT_ID;
static std::set<int> used_ids;

/**
 * @brief Returns an auto-generated unique Cell ID.
 * @details This method is intended as a utility method for users writing
 *          OpenMOC input files. The method makes use of a static Cell
 *          ID which is incremented each time the method is called to enable
 *          unique generation of monotonically increasing IDs. The method's
 *          first ID begins at 1,000,000. Hence, user-defined Cell IDs greater
 *          than or equal to 1,000,000 are prohibited.
 */
int cell_id() {
  int id = auto_id;
  auto_id++;
  while (used_ids.find(id) != used_ids.end()) {
    id = auto_id;
    auto_id++;
  }
  return id;
}


/**
 * @brief Resets the auto-generated unique Cell ID counter to 1,000,000.
 */
void reset_cell_id() {
  auto_id = DEFAULT_INIT_ID;
}


/**
 * @brief Maximize the auto-generated unique Cell ID counter.
 * @details This method updates the auto-generated unique Cell ID
 *          counter if the input parameter is greater than the present
 *          value. This is useful for the OpenMC compatibility module
 *          to ensure that the auto-generated Cell IDs do not
 *          collide with those created in OpenMC.
 * @param cell_id the id assigned to the auto-generated counter
 */
void maximize_cell_id(int cell_id) {
  if (cell_id > auto_id)
    auto_id = cell_id;
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

  /* Add the ID to the used set */
  used_ids.insert(_id);

  _uid = _n;
  _n++;

  _name = NULL;
  setName(name);

  _cell_type = UNFILLED;
  _region = NULL;
  _current_region = NULL;
  _fill = NULL;
  _volume = 0.;
  _num_instances = 0;

  _rotated = false;
  memset(&_rotation, 0, 3*sizeof(double));
  memset(&_rotation_matrix, 0, 9*sizeof(double));

  _translated = false;
  memset(&_translation, 0, 3*sizeof(double));

  _num_rings = 0;
  _num_sectors = 0;
  _parent = NULL;
}


/**
 * @brief Destructor clears vector of Surface pointers bounding the Cell.
 */
Cell::~Cell() {

  if (_name != NULL)
    delete [] _name;
  if (_region != NULL)
    delete _region;
  /* Materials are deleted separately from cells, since multiple cells
      can share the same material */
  /* Universes are also deleted separately, since Universes can have been
     defined in your input scripts, rather than loaded from a Geometry file */
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
 * @brief Return the user-defined name of the Cell.
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
 * @brief Return the Cell's Region, its spatial domain.
 * @return the Cell's Region
 */
Region* Cell::getRegion() {
  return _region;
}


/**
 * @brief Return the aggregate volume/area of all instances of this Cell.
 * @details The volume/area of the Cell is computed from track segments which
 *          overlap this Cell during track generation.
 * @return the volume/area of the Cell
 */
double Cell::getVolume() {
  return _volume;
}


/**
 * @brief Return a boolean indicating whether the Cell has been rotated.
 * @return whether the Cell has been rotated
 */
bool Cell::isRotated() {
  return _rotated;
}


/**
 * @brief Return the number of instances of this Cell in the Geometry.
 * @details The number of instances of this Cell in the Geometry is
 *          determined during track generation.
 * @return the number of cell instances
 */
int Cell::getNumInstances() {
  return _num_instances;
}


/**
 * @brief Return a boolean indicating whether the Cell has been translated.
 * @return whether the Cell has been translated
 */
bool Cell::isTranslated() {
  return _translated;
}


/**
 * @brief Get the rotation angle about the x-axis in degrees.
 * @param units the angular units in "radians" or "degrees" (default)
 * @return the rotation angle about the x-axis
 */
double Cell::getPhi(std::string units) {
  std::string degrees("degrees");
  std::string radians("radians");

  /* Return phi in degrees or radians */
  if (degrees.compare(units) == 0)
    return _rotation[0] / M_PI * 180.;
  else if (radians.compare(units) == 0)
    return _rotation[0];

  log_printf(ERROR, "Unable to return phi in units %s", units.c_str());
  return -1;
}


/**
 * @brief Get the rotation angle about the y-axis in degrees.
 * @param units the angular units in "radians" or "degrees" (default)
 * @return the rotation angle about the y-axis
 */
double Cell::getTheta(std::string units) {
  std::string degrees("degrees");
  std::string radians("radians");

  /* Return theta in degrees or radians */
  if (degrees.compare(units) == 0)
    return _rotation[1] / M_PI * 180.;
  else if (radians.compare(units) == 0)
    return _rotation[1];

  log_printf(ERROR, "Unable to return theta in units %s", units.c_str());
  return -1;
}


/**
 * @brief Get the rotation angle about the z-axis in degrees.
 * @param units the angular units in "radians" or "degrees" (default)
 * @return the rotation angle about the z-axis
 */
double Cell::getPsi(std::string units) {
  std::string degrees("degrees");
  std::string radians("radians");

  /* Return psi in degrees or radians */
  if (degrees.compare(units) == 0)
    return _rotation[2] / M_PI * 180.;
  else if (radians.compare(units) == 0)
    return _rotation[2];

  log_printf(ERROR, "Unable to return psi in units %s", units.c_str());
  return -1;
}


/**
 * @brief Return pointer to array for the rotation matrix.
 * @return a pointer to an array of rotation angles
 */
double* Cell::getRotationMatrix() {
  return _rotation_matrix;
}


/**
 * @brief Fills an array with the rotation angles for x, y and z.
 * @details This class method is intended to be called by the OpenMOC
 *          Python OpenMC compatiblity module. Although this method appears to
 *          require two arguments, in reality it only requires one due to SWIG
 *          and would be called from within Python as follows:
 *
 * @code
 *          rotation = cell.getRotation(3)
 * @endcode
 *
 * @param rotations an array of rotation angles of length 3 for x, y and z
 * @param num_axes the number of axes (this must always be 3)
 * @param units the angular units in "radians" or "degrees" (default)
 */
void Cell::retrieveRotation(double* rotations, int num_axes,
                            std::string units) {
  if (num_axes != 3)
    log_printf(ERROR, "Unable to get rotation with %d axes for Cell %d. "
               "The rotation array should be length 3.", num_axes, _id);

  std::string degrees("degrees");
  std::string radians("radians");

  /* Return psi in degrees or radians */
  for (int i=0; i < 3; i++) {
    if (degrees.compare(units) == 0)
      rotations[i] = _rotation[i] * 180. / M_PI;
    else if (radians.compare(units) == 0)
      rotations[i] = _rotation[i];
    else
      log_printf(ERROR, "Unable to return rotation in units %s", units.c_str());
  }
}


/**
 * @brief Return pointer to array for the translations along x, y and z.
 * @return a pointer to an array of translations
 */
double* Cell::getTranslation() {
  return _translation;
}


/**
 * @brief Fills an array with the translations along x, y and z.
 * @details This class method is intended to be called by the OpenMOC
 *          Python OpenMC compatiblity module. Although this method appears to
 *          require two arguments, in reality it only requires one due to SWIG
 *          and would be called from within Python as follows:
 *
 * @code
 *          translation = cell.retrieveTranslation(3)
 * @endcode
 *
 * @param translations an array of translations of length 3 for x, y and z
 * @param num_axes the number of axes (this must always be 3)
 */
void Cell::retrieveTranslation(double* translations, int num_axes) {
  if (num_axes != 3)
    log_printf(ERROR, "Unable to get translation with %d axes for Cell %d. "
               "The translation array should be length 3.", num_axes, _id);

  for (int i=0; i < 3; i++)
    translations[i] = _translation[i];
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

  if (_rotated || _translated)
    log_printf(WARNING, "Cell rotations and translations are ignored, minX may "
               "be inaccurate");

  double min_x = -std::numeric_limits<double>::infinity();

  /* Look in region for minimum */
  if (_region != NULL)
    min_x = _region->getMinX();

  /* If region has an infinite min_x, it could be that some Halfspaces are only
     kept in the Parent's region */
  if (getParent() != NULL &&
      std::abs(min_x) == std::numeric_limits<double>::infinity())
    min_x = getParent()->getMinX();

  return min_x;
}


/**
 * @brief Return the maximum reachable x-coordinate in the Cell.
 * @return the maximum x-coordinate
 */
double Cell::getMaxX() {

  if (_rotated || _translated)
    log_printf(WARNING, "Cell rotations and translations are ignored, maxX may "
               "be inaccurate");

  double max_x = +std::numeric_limits<double>::infinity();

  /* Look in region for maximum */
  if (_region != NULL)
    max_x = _region->getMaxX();

  /* If region has an infinite max_x, it could be that some Halfspaces are only
     kept in the Parent's region */
  if (getParent() != NULL &&
      std::abs(max_x) == std::numeric_limits<double>::infinity())
    max_x = getParent()->getMaxX();

  return max_x;
}


/**
 * @brief Return the minimum reachable y-coordinate in the Cell.
 * @return the minimum y-coordinate
 */
double Cell::getMinY() {

  if (_rotated || _translated)
    log_printf(WARNING, "Cell rotations and translations are ignored, minY may "
               "be inaccurate");

  double min_y = -std::numeric_limits<double>::infinity();

  /* Look in region for minimum */
  if (_region != NULL)
    min_y = _region->getMinY();

  /* If region has an infinite min_y, it could be that some Halfspaces are only
     kept in the Parent's region */
  if (getParent() != NULL &&
      std::abs(min_y) == std::numeric_limits<double>::infinity())
    min_y = getParent()->getMinY();

  return min_y;
}


/**
 * @brief Return the maximum reachable y-coordinate in the Cell.
 * @return the maximum y-coordinate
 */
double Cell::getMaxY() {

  if (_rotated || _translated)
    log_printf(WARNING, "Cell rotations and translations are ignored, maxY may "
               "be inaccurate");

  double max_y = +std::numeric_limits<double>::infinity();

  /* Look in region for maximum */
  if (_region != NULL)
    max_y = _region->getMaxY();

  /* If region has an infinite max_y, it could be that some Halfspaces are only
     kept in the Parent's region */
  if (getParent() != NULL &&
      std::abs(max_y) == std::numeric_limits<double>::infinity())
    max_y = getParent()->getMaxY();

  return max_y;
}


/**
 * @brief Return the minimum reachable z-coordinate in the Cell.
 * @return the minimum z-coordinate
 */
double Cell::getMinZ() {

  if (_rotated || _translated)
    log_printf(DEBUG, "Cell rotations and translations are ignored, minZ may "
               "be inaccurate");

  double min_z = -std::numeric_limits<double>::infinity();

  /* Look in region for minimum */
  if (_region != NULL)
    min_z = _region->getMinZ();

  /* If region has an infinite min_z, it could be that some Halfspaces are only
     kept in the Parent's region */
  if (getParent() != NULL &&
      std::abs(min_z) == std::numeric_limits<double>::infinity())
    min_z = getParent()->getMinZ();

  return min_z;
}


/**
 * @brief Return the maximum reachable z-coordinate in the Cell.
 * @return the maximum z-coordinate
 */
double Cell::getMaxZ() {

  if (_rotated || _translated)
    log_printf(DEBUG, "Cell rotations and translations are ignored, maxZ may "
               "be inaccurate");

  double max_z = +std::numeric_limits<double>::infinity();

  /* Look in region for maximum */
  if (_region != NULL)
    max_z = _region->getMaxZ();

  /* If region has an infinite max_z, it could be that some Halfspaces are only
     kept in the Parent's region */
  if (getParent() != NULL &&
      std::abs(max_z) == std::numeric_limits<double>::infinity())
    max_z = getParent()->getMaxZ();

  return max_z;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE) at
 *        the minimum reachable x-coordinate in the Cell.
 * @return the boundary condition at the minimum x-coordinate
 */
boundaryType Cell::getMinXBoundaryType() {

  //FIXME Cell rotations are not taken into account
  boundaryType boundary = BOUNDARY_NONE;

  /* Look in region for min x boundary */
  if (_region != NULL)
    boundary = _region->getMinXBoundaryType();

  /* If no boundary was found in the region, it may be that the boundary
     is in the Parent cell's region */
  if (getParent() != NULL && boundary == BOUNDARY_NONE)
    boundary = getParent()->getMinXBoundaryType();

  return boundary;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE) at
 *        the maximum reachable x-coordinate in the Cell.
 * @return the boundary condition at the maximum x-coordinate
 */
boundaryType Cell::getMaxXBoundaryType() {

  //FIXME Cell rotations are not taken into account
  boundaryType boundary = BOUNDARY_NONE;

  /* Look in region for max x boundary */
  if (_region != NULL)
    boundary = _region->getMaxXBoundaryType();

  /* If no boundary was found in the region, it may be that the boundary
     is in the Parent cell's region */
  if (getParent() != NULL && boundary == BOUNDARY_NONE)
    boundary = getParent()->getMaxXBoundaryType();

  return boundary;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE) at
 *        the minimum reachable y-coordinate in the Cell.
 * @return the boundary condition at the minimum y-coordinate
 */
boundaryType Cell::getMinYBoundaryType() {

  //FIXME Cell rotations are not taken into account
  boundaryType boundary = BOUNDARY_NONE;

  /* Look in region for min y boundary */
  if (_region != NULL)
    boundary = _region->getMinYBoundaryType();

  /* If no boundary was found in the region, it may be that the boundary
     is in the Parent cell's region */
  if (getParent() != NULL && boundary == BOUNDARY_NONE)
    boundary = getParent()->getMinYBoundaryType();

  return boundary;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE) at
 *        the maximum reachable y-coordinate in the Cell.
 * @return the boundary condition at the maximum y-coordinate
 */
boundaryType Cell::getMaxYBoundaryType() {

  //FIXME Cell rotations are not taken into account
  boundaryType boundary = BOUNDARY_NONE;

  /* Look in region for max y boundary */
  if (_region != NULL)
    boundary = _region->getMaxYBoundaryType();

  /* If no boundary was found in the region, it may be that the boundary
     is in the Parent cell's region */
  if (getParent() != NULL && boundary == BOUNDARY_NONE)
    boundary = getParent()->getMaxYBoundaryType();

  return boundary;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE) at
 *        the minimum reachable z-coordinate in the Cell.
 * @return the boundary condition at the minimum z-coordinate
 */
boundaryType Cell::getMinZBoundaryType() {

  //FIXME Cell rotations are not taken into account
  boundaryType boundary = BOUNDARY_NONE;

  /* Look in region for min z boundary */
  if (_region != NULL)
    boundary = _region->getMinZBoundaryType();

  /* If no boundary was found in the region, it may be that the boundary
     is in the Parent cell's region */
  if (getParent() != NULL && boundary == BOUNDARY_NONE)
    boundary = getParent()->getMinZBoundaryType();

  return boundary;
}


/**
 * @brief Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE) at
 *        the maximum reachable z-coordinate in the Cell.
 * @return the boundary condition at the maximum z-coordinate
 */
boundaryType Cell::getMaxZBoundaryType() {

  //FIXME Cell rotations are not taken into account
  boundaryType boundary = BOUNDARY_NONE;

  /* Look in region for max z boundary */
  if (_region != NULL)
    boundary = _region->getMaxZBoundaryType();

  /* If no boundary was found in the region, it may be that the boundary
     is in the Parent cell's region */
  if (getParent() != NULL && boundary == BOUNDARY_NONE)
    boundary = getParent()->getMaxZBoundaryType();

  return boundary;
}


/**
 * @brief Return the number of Surfaces in the Cell.
 * @return the number of Surfaces
 */
int Cell::getNumSurfaces() const {
  std::map<int, Halfspace*> all_surfaces;
  if (_region != NULL)
    all_surfaces = _region->getAllSurfaces();
  return all_surfaces.size();
}


/**
 * @brief Return the std::map of Halfspace object pointers for all
 *        surfaces within the Region bounding the Cell.
 * @return std::map of Halfspace object pointers with surface ID as a key.
 */
std::map<int, Halfspace*> Cell::getSurfaces() const {
  std::map<int, Halfspace*> all_surfaces;
  if (_region != NULL)
    all_surfaces = _region->getAllSurfaces();
  return all_surfaces;
}


/**
 * @brief Return the std::vector of neighbor Cells to this Cell.
 * @return std::vector of neighbor Cell pointers
 */
std::vector<Cell*> Cell::getNeighbors() const {
  return _neighbors;
}


/**
 * @brief Return true if the Cell has a parent and false otherwise.
 * @return whether the Cell has a parent Cell
 */
bool Cell::hasParent() {
  if (_parent == NULL)
    return false;
  else
    return true;
}


/**
 * @brief Return this Cell's parent Cell.
 * @details If no parent Cell has been assigned from Cell cloning, then
 *          NULL is returned.
 * @return a pointer to the parent Cell
 */
Cell* Cell::getParent() {
  return _parent;
}


/**
 * @brief Get the oldest ancestor Cell for this Cell.
 * @details This method traverses the linked list of parent Cells to find
 *          the one at the root node. The oldest ancestor Cell is likely the
 *          one created by the user at runtime, while intermediate ancestors
 *          were created during radial and angular spatial discretization.
 * @return this Cell's oldest ancestor Cell
 */
Cell* Cell::getOldestAncestor() {

  /* If this Cell has no parent, return NULL */
  if (_parent == NULL)
    return _parent;

  /* Otherwise, navigate to the first parent Cell */
  else {
    /* Traverse linked list to the root node */
    Cell* curr = _parent;
    while (curr->hasParent())
      curr = curr->getParent();

    return curr;
  }
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
 * @brief Sets the Region this Cell lives in.
 * @details NOTE: This method deep copies the Region and stores
 *          the copy. Any changes made to the Region will not be
 *          reflected in the Region copy stored by the Cell.
 * @param region the Region bounding the Cell
 */
void Cell::setRegion(Region* region) {
  _region = region->clone();
}


/**
 * @brief Set the volume/area of the Cell.
 * @param volume the volume/area of the Cell
 */
void Cell::setVolume(double volume) {
  _volume = volume;
}


/**
 * @brief Increment the volume/area of the Cell by some amount.
 * @details This routine is called by the TrackGenerator during track
 *          generation and segmentation.
 * @param volume the amount to increment the current volume by
 */
void Cell::incrementVolume(double volume) {
  _volume += volume;
}


/**
 * @brief Set the number of instances of this Cell.
 * @param num_instances the number of instances of this Cell in the Geometry
 */
void Cell::setNumInstances(int num_instances) {
  _num_instances = num_instances;
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
 *          rotation = numpy.array([0., 0., 90.])
 *          cell = openmoc.Cell()
 *          cell.setRotation(rotation)
 * @endcode
 *
 * @param rotation the array of rotation angles
 * @param num_axes the number of axes (this must always be 3)
 * @param units the angular units in "radians" or "degrees" (default)
 */
void Cell::setRotation(double* rotation, int num_axes, std::string units) {

  if (num_axes != 3)
    log_printf(ERROR, "Unable to set rotation with %d axes for Cell %d. "
               "The rotation array should be length 3.", num_axes, _id);

  if (_cell_type == MATERIAL)
    log_printf(ERROR, "Rotations cannot be defined on material cells, "
               "create a rotated cell, containing a universe, containing "
               "the material cell.");

  std::string degrees("degrees");
  std::string radians("radians");

  /* Store rotation angles in radians */
  for (int i=0; i < 3; i++) {
    if (degrees.compare(units) == 0)
      _rotation[i] = rotation[i] * M_PI / 180.;
    else if (radians.compare(units) == 0)
      _rotation[i] = rotation[i];
    else
      log_printf(ERROR, "Unable to set rotation with units %s", units.c_str());
  }

  /* Use pitch-roll-yaw convention according to eqns 51-59 on Wolfram:
   * http://mathworld.wolfram.com/EulerAngles.html */
  double theta = _rotation[0];
  double psi = _rotation[1];
  double phi = _rotation[2];

  /* Calculate rotation matrix based on angles given */
  /* Indexed by (y,x) since the universe array is indexed by (z,y,z) */
  _rotation_matrix[0] = cos(theta) * cos(phi);
  _rotation_matrix[1] = cos(theta) * sin(phi);
  _rotation_matrix[2] = -sin(theta);
  _rotation_matrix[3] = sin(psi) * sin(theta) * cos(phi) -
                        cos(psi) * sin(phi);
  _rotation_matrix[4] = sin(psi) * sin(theta) * sin(phi) +
                        cos(psi) * cos(phi);
  _rotation_matrix[5] = cos(theta) * sin(psi);
  _rotation_matrix[6] = cos(psi) * sin(theta) * cos(phi) +
                        sin(psi) * sin(phi);
  _rotation_matrix[7] = cos(psi) * sin(theta) * sin(phi) -
                        sin(psi) * cos(phi);
  _rotation_matrix[8] = cos(theta) * cos(psi);

  _rotated = true;
}


/**
 * @brief Increment the number of instances of this Cell.
 * @details This routine is called by the TrackGenerator during track
 *          generation and segmentation.
 */
void Cell::incrementNumInstances() {
  _num_instances++;
}

/**
 * @brief Set the Cell's translation along the x, y and z axes.
 * @details This method is a helper function to allow OpenMOC users to assign
 *          the Cell's translations in Python. A user must initialize a
 *          length 3 NumPy array as input to this function. This function then
 *          stores the data values in the NumPy array in the Cell's translation
 *          array. An example of how this function might be called in Python
 *          is as follows:
 *
 * @code
 *          translation = numpy.array([0.25, 0.25, 0.])
 *          cell = openmoc.Cell()
 *          cell.setTranslation(translation)
 * @endcode
 *
 * @param translation the array of translations
 * @param num_axes the number of axes (this must always be 3)
 */
void Cell::setTranslation(double* translation, int num_axes) {

  if (num_axes != 3)
    log_printf(ERROR, "Unable to set translation for %d axes for Cell %d. "
               "The translation array should be length 3.", num_axes, _id);

  if (_cell_type == MATERIAL)
    log_printf(ERROR, "Translations cannot be defined on material cells, "
               "create a translated cell, containing a universe, containing "
               "the current material cell ID %d.", _id);

  for (int i=0; i < 3; i++)
    _translation[i] = translation[i];

  _translated = true;
}


/**
 * @brief Set the Cell's number of rings.
 * @param num_rings the number of rings in this Cell
 * @param inner_radius add (0,0,r) ZCylinder to ringify around (default none)
 */
void Cell::setNumRings(int num_rings, double inner_radius) {
  if (num_rings < 0)
    log_printf(ERROR, "Unable to give %d rings to Cell %d since this is "
               "a negative number", num_rings, _id);

  if (num_rings == 1)
    _num_rings = 0;
  else
    _num_rings = num_rings;

  if (inner_radius >= 0) {
    // Create an inner Z cylinder at coordinate (0, 0)
    ZCylinder* zcylinder = new ZCylinder(0, 0, inner_radius);
    zcylinder->setName("innermost ring cylinder");

    // The cell is defined as being outside this cylinder
    addSurface(+1, zcylinder);
  }
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
 * @brief Assign a parent Cell to this Cell.
 * @details This is used by Cell cloning when applied for radial and
 *          angular discretization.
 * @param parent a pointer to the parent Cell
 */
void Cell::setParent(Cell* parent) {
  _parent = parent;
}


/**
 * @brief Insert a Surface into this Cell's bounding Region, assuming that an
 *        intersection between the region and the halfspace is desired.
 * @param halfspace the Surface halfspace (+/-1)
 * @param surface a pointer to the Surface
 */
void Cell::addSurface(int halfspace, Surface* surface) {

  if (halfspace != -1 && halfspace != +1)
    log_printf(ERROR, "Unable to add surface %d to cell %d since the halfspace"
               " %d is not -1 or 1", surface->getId(), _id, halfspace);

  Halfspace* new_halfspace = new Halfspace(halfspace, surface);

  /* Assign the Halfspace as the Cell's Region if it has none */
  if (_region == NULL)
    _region = new_halfspace;
  else{
    if (dynamic_cast<Intersection*>(_region))
      _region->addNode(new_halfspace, false);
    else {
      Intersection* intersection = new Intersection();
      intersection->addNode(_region, false);
      intersection->addNode(new_halfspace, false);
      _region = intersection;
    }
  }
}


 /**
  * @brief Insert a Surface into this Cell's bounding Region.
  * @param halfspace the Surface halfspace (+/-1)
  * @param surface a pointer to the Surface
  */
 void Cell::addSurfaceInRegion(int halfspace, Surface* surface) {

   /* Create a new halfspace */
   Halfspace* new_halfspace = new Halfspace(halfspace, surface);

  /* Assign the Halfspace as the Cell's Region if it has none */
  if (_region == NULL)
    _region = new_halfspace;
  /* Else add Halfspace to region */
  else
    _current_region->addNode(new_halfspace, false);
}


/**
 * @brief Removes a Surface from this Cell's container of bounding Surfaces.
 * @param surface a pointer to the Surface to remove
 */
void Cell::removeSurface(Surface* surface) {

  //FIXME This map cannot be modified, it's a const
  std::map<int, Halfspace*> surfaces = getSurfaces();
  if (surface != NULL && surfaces.find(surface->getId()) != surfaces.end()) {
    delete surfaces[surface->getId()];
    surfaces.erase(surface->getId());
  }
}


 /**
  * @brief Insert a logical node (intersection or union) into the cell region.
  * @details This method creates a node in a tree of regions. The leaves, or the
  *          nodes at the very bottom of the tree, are halfspaces. The region
             is defined by that tree.
  * @param region_type the logical operation
  */
 void Cell::addLogicalNode(int region_type) {

   /* Create new region if void */
   if (_region == NULL) {
     if (region_type == INTERSECTION) {
       Intersection* intersection = new Intersection();
       _region = intersection;
     }
     else if (region_type == UNION) {
       Union* _union = new Union();
       _region = _union;
     }
     else if (region_type == COMPLEMENT) {
       Complement* complement = new Complement();
       _region = complement;
     }
     _current_region = _region;
  }
  /* Add node under current region, and move current region to said node */
  else {
    if (region_type == INTERSECTION) {
      Intersection* intersection = new Intersection();
      intersection->setParentRegion(_current_region);
      _current_region->addNode(intersection, false);
      _current_region = intersection;
    }
    else if (region_type == UNION) {
      Union* _union = new Union();
      _union->setParentRegion(_current_region);
      _current_region->addNode(_union, false);
      _current_region = _union;
    }
    else if (region_type == COMPLEMENT) {
      Complement* complement = new Complement();
      complement->setParentRegion(_current_region);
      _current_region->addNode(complement, false);
      _current_region = complement;
    }
  }
}


/**
 * @brief Climb up the logical tree of regions.
 */
void Cell::goUpOneRegionLogical() {
  _current_region = _current_region->getParentRegion();
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
 * @brief Determines whether a Point is contained inside a Cell.
 * @details Queries the Region bounding the Cell to determine if the Point
 *          is within the Region. This point is only inside the Cell if it
 *          is on the same side of every Surface bounding the Cell.
 * @param point a pointer to a Point
 * @return true if the Point is inside the Cell; otherwise false
 */
bool Cell::containsPoint(Point* point) {

  /* If a FILL Cell, query the filling Universe or Lattice */
  if (_region == NULL) {
    if (_cell_type == FILL) {
      Universe* univ = static_cast<Universe*>(_fill);
      if (univ->getType() == SIMPLE)
        return univ->containsPoint(point);
      else {
        Lattice* latt = static_cast<Lattice*>(_fill);
        return latt->containsPoint(point);
      }
    }
    else
      return true;
  }

  /* If not, query the Cell's bounding Region */
  else
    return _region->containsPoint(point);
}


/**
 * @brief Determines whether a Point is contained inside a Cell.
 * @details Queries each Surface inside the Cell to determine if the Point
 *          is on the same side of the Surface. This Point is only inside
 *          the Cell if it is on the same side of every Surface in the Cell.
 * @param coords a pointer to a localcoord
 * @return whether the cell contains the coords or not
 */
bool Cell::containsCoords(LocalCoords* coords) {
  return containsPoint(coords->getPoint());
}


/**
 * @brief Computes the minimum distance to a Surface in the Cell's Region from
 *        a point with a given trajectory at a certain angle stored in a
 *        LocalCoords object.
 * @details If the trajectory will not intersect any of the Surfaces in the
 *          Cell returns INFINITY.
 * @param coords a pointer to a localcoords
 * @return distance to nearest intersection with the cell's region boundaries
 */
double Cell::minSurfaceDist(LocalCoords* coords) {
  if (_region == NULL)
    return INFINITY;
  else
    return _region->minSurfaceDist(coords);
}


/**
 * @brief Computes the minimum distance to a Surface from a Point with a given
 *        trajectory at a certain angle.
 * @details If the trajectory will not intersect any of the Surfaces in the
 *          Cell returns INFINITY.
 * @param point the Point of interest
 * @param azim the azimuthal angle of the trajectory
 *        (in radians from \f$[0,2\pi]\f$)
 * @param polar the polar angle of the trajectory
 *        (in radians from \f$[0,\pi]\f$)
 * @return distance to nearest intersection with the cell's region boundaries
 */
double Cell::minSurfaceDist(Point* point, double azim, double polar) {
  if (_region == NULL)
    return INFINITY;
  else
    return _region->minSurfaceDist(point, azim, polar);
}


/**
 * @brief Returns true if this Cell is filled with a fissionable Material.
 * @details If the Cell is filled by a Material, this method will simply query
 *          the filling Material. If the Cell is filled by a Universe, this
 *          method will consider any Materials filling those Cells contained
 *          by the filling Universe. This method should not be called prior to
 *          the calling of the Geometry::computeFissionability() method.
 * @return true if contains a fissionable Material
 */
bool Cell::isFissionable() {
  if (_fill != NULL)
    if (_cell_type == FILL)
      return ((Universe*)_fill)->isFissionable();
    else
      return ((Material*)_fill)->isFissionable();
  else
    return false;
}


/**
 * @brief Create a duplicate of the Cell.
 * @return a pointer to the clone
 */
Cell* Cell::clone(bool clone_region) {

  /* Construct new Cell */
  Cell* new_cell = new Cell();
  new_cell->setName(_name);
  new_cell->setNumRings(_num_rings);
  new_cell->setNumSectors(_num_sectors);
  new_cell->setParent(this);

  if (_cell_type == MATERIAL)
    new_cell->setFill((Material*)_fill);
  else
    new_cell->setFill((Universe*)_fill);

  /* Clone the Cell's Region (cloning is done in setRegion) */
  if (_region != NULL && clone_region)
    new_cell->setRegion(_region);

  if (_rotated)
    new_cell->setRotation(_rotation, 3, "radians");
  if (_translated)
    new_cell->setTranslation(_translation, 3);

  return new_cell;
}


/**
 * @brief Subdivides the Cell into clones for fuel pin angular sectors.
 * @param subcells an empty vector to store all subcells
 */
void Cell::sectorize(std::vector<Cell*>& subcells) {

  /* If the user didn't request any sectors, don't make any */
  if (_num_sectors == 0)
    return;

  double azim_angle;
  double delta_azim = 2. * M_PI / _num_sectors;
  double A, B;

  /* A container for each of the bounding planes for the sector Cells */
  std::vector<Plane*> planes;

  log_printf(DEBUG, "Sectorizing Cell %d with %d sectors", _id, _num_sectors);

  /* Create each of the bounding planes for the sector Cells */
  for (int i=0; i < _num_sectors; i++) {

    /* Figure out the angle for this plane */
    azim_angle = i * delta_azim + M_PI / 4.0;

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
    Cell* sector = clone(false);

    sector->setNumSectors(0);
    sector->setNumRings(0);

    log_printf(DEBUG, "Creating a new sector Cell %d for Cell %d",
               sector->getId(), _id);

    if (_num_sectors != 2) {
      /* Add new bounding planar Surfaces to the clone */
      sector->addSurface(+1, planes.at(i));

      if (i+1 < _num_sectors)
        sector->addSurface(-1, planes.at(i+1));
      else
        sector->addSurface(-1, planes.at(0));
    }
    else {
      /* For _num_sectors==2, planes[0] and planes[1] are actually the same but
         opposite direction, so the two adjacent sectors will have the same
         Halfspace value, which will cause trouble when a point is on the plane.
         This is to avoid this trouble. */
      int halfspace = (i==0? +1 : -1);
      sector->addSurface(halfspace, planes.at(0));
    }

    /* Store the clone in the parent Cell's container of sector Cells */
    subcells.push_back(sector);
  }
}


/**
 * @brief Subdivides the Cell into clones for fuel pin rings.
 * @param subcells an empty vector to store all subcells
 * @param max_radius the maximum allowable radius used in the subdivisions
 */
void Cell::ringify(std::vector<Cell*>& subcells, double max_radius) {

  /* If the user didn't request any rings, don't make any */
  if (_num_rings == 0)
        return;

  int num_zcylinders = 0;
  ZCylinder* zcylinder1 = NULL;
  ZCylinder* zcylinder2 = NULL;
  double radius1 = max_radius;
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
  std::map<int, Halfspace*>::iterator iter1;
  std::map<int, Halfspace*> _surfaces = getSurfaces();
  for (iter1=_surfaces.begin(); iter1 != _surfaces.end(); ++iter1) {

    /* Determine if any of the Surfaces is a ZCylinder */
    if (iter1->second->getSurface()->getSurfaceType() == ZCYLINDER) {
      int halfspace = iter1->second->getHalfspace();
      ZCylinder* zcylinder = static_cast<ZCylinder*>(iter1->second->getSurface());

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
  if (num_zcylinders == 0) {
    log_printf(WARNING, "Unable to ringify Cell %d: %s since it does not "
              "contain any ZCYLINDER type Surface(s)", _id, _name);
    return;
  }

  if (num_zcylinders > 2)
    log_printf(ERROR, "Unable to ringify Cell %d since it "
               "contains more than 2 ZCYLINDER Surfaces", _id);

  if (fabs(x1 - x2) > FLT_EPSILON && num_zcylinders == 2)
    log_printf(ERROR, "Unable to ringify Cell %d since it contains "
               "ZCylinder %d centered at x=%f and ZCylinder %d at x=%f. "
               "Both ZCylinders must have the same center.",
               _id, zcylinder1->getId(), x1, zcylinder2->getId(), x2);

  if (fabs(y1 - y2) > FLT_EPSILON && num_zcylinders == 2)
    log_printf(ERROR, "Unable to ringify Cell %d since it contains "
               "ZCylinder %d centered at y=%f and ZCylinder %d at y=%f. "
               "Both ZCylinders must have the same center.",
               _id, zcylinder1->getId(), y1, zcylinder2->getId(), y2);

  if (radius1 <= radius2)
    log_printf(ERROR, "Unable to ringify Cell %d since it contains 2 "
               "disjoint ZCYLINDER Surfaces: halfspace %d for ZCylinder %d "
               "and halfspace %d for ZCylinder %d. Switch the signs of "
               "the 2 halfspaces for each Surface.", _id, halfspace1,
               zcylinder1->getId(), halfspace2, zcylinder2->getId());

  /* Loop over ZCylinders and create a new Cell clone for each ring */
  std::vector<ZCylinder*>::iterator iter2;
  std::vector<Cell*>::iterator iter3;

  /* Compute the increment, either by radius or area, to use to construct
   * the concentric rings */
  double increment;

  /* If there is no outer bounding surface, make the rings have the same
   * radius increment (e.g. moderator in a pin cell universe). */
  if (halfspace1 == 0){
    increment = fabs(radius1 - radius2) / _num_rings;

    /* Heuristic to improve area-balancing for low number of rings */
    if (fabs(radius1 - max_radius) < FLT_EPSILON && _num_rings < 3)
      increment = 1.5 * (radius1 - radius2) / _num_rings;
  }

  /* If there is an outer bounding surface, make the rings have the same
   * area (e.g. fuel in a pin cell universe).*/
  else
    increment = M_PI * fabs(radius1*radius1 - radius2*radius2) / _num_rings;

  /* Generate successively smaller ZCylinders */
  for (int i=0; i < _num_rings-1; i++) {

    /* Compute the outer radius of the next ring */
    if (halfspace1 == 0)
      radius2 = radius1 - increment;
    else
      radius2 = sqrt(radius1 * radius1 - (increment / M_PI));

    ZCylinder* zcylinder = new ZCylinder(x1, y1, radius1);
    zcylinders.push_back(zcylinder);
    radius1 = radius2;
  }

  /* Store smallest, innermost ZCylinder */
  ZCylinder* zcylinder = new ZCylinder(x1, y1, radius1);
  zcylinders.push_back(zcylinder);

  /* Create ring Cells with successively smaller ZCylinders */
  for (iter2 = zcylinders.begin(); iter2 != zcylinders.end(); ++iter2) {

    /* Create ZCylinders for each of the sectorized Cells */
    if (subcells.size() != 0) {
      for (iter3 = subcells.begin(); iter3 != subcells.end(); ++iter3) {
        log_printf(DEBUG, "Creating a new ring in sector Cell ID=%d",
                   (*iter3)->getId());

        /* Create a new Cell clone */
        Cell* ring = (*iter3)->clone();
        ring->setNumSectors(0);
        ring->setNumRings(0);

        /* Add ZCylinder only if this is not the outermost ring in an
         * unbounded Cell (i.e. the moderator in a fuel pin cell) */
        if ((*iter2)->getRadius() < max_radius)
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
      Cell* ring = clone(false);
      ring->setNumSectors(0);
      ring->setNumRings(0);

      /* Add ZCylinder only if this is not the outermost ring in an
       * unbounded Cell (i.e. the moderator in a fuel pin cell) */
      if ((*iter2)->getRadius() < max_radius)
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
  subcells.clear();
  subcells.insert(subcells.end(), rings.begin(), rings.end());
}


/**
 * @brief Subdivides a Cell into rings and sectors aligned with the z-axis.
 * @details This method uses the Cell's clone method to produce a vector of
 *          this Cell's subdivided ring and sector Cells.
 * @param max_radius the maximum allowable radius used in the subdivisions
 * @return a vector of Cell pointers to the new subdivided Cells
 */
void Cell::subdivideCell(double max_radius) {

  /** A container of all Cell clones created for rings and sectors */
  std::vector<Cell*> subcells;

  sectorize(subcells);
  ringify(subcells, max_radius);

  /* Put any ring / sector subcells in a new Universe fill */
  if (subcells.size() != 0) {

    /* Create a new Universe to contain all of the subcells */
    Universe* new_fill = new Universe();

    /* Add each subcell to the new Universe fill */
    std::vector<Cell*>::iterator iter;
    for (iter = subcells.begin(); iter != subcells.end(); ++iter)
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

  /* Add this Cell to the neighbor lists of all this cell's surfaces */
  std::map<int, Halfspace*>::iterator iter;
  std::map<int, Halfspace*> _surfaces = getSurfaces();
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
  else if (_cell_type == MATERIAL) {
    string << ", type = MATERIAL"
           << ", fill id = " << static_cast<Material*>(_fill)->getId();
  }
  else
    string << ", type = UNFILLED";

  if (_rotated) {
    string << ", (rotation = " << getPhi() << ", ";
    string << getTheta() << ", " << getPsi() << ")";
  }
  if (_translated) {
    string << ", (translation = " << _translation[0] << ", ";
    string << _translation[1] << ", " << _translation[2] << ")";
  }

  string << ", # surfaces = " << getNumSurfaces();

  /** Add string data for the Surfaces in this Cell */
  std::map<int, Halfspace*>::iterator iter;
  std::map<int, Halfspace*> _surfaces = getSurfaces();
  string << ", Surfaces: ";
  for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
    string << std::showpos << "\nhalfspace = " << iter->second->_halfspace
           << ", " << iter->second->_surface->toString();

  return string.str();
}


/**
 * @brief Prints a string representation of all of the Cell's attributes
 *        to the console.
 */
void Cell::printString() {
  log_printf(NORMAL, toString().c_str());
}


/**
 * @brief Obtain and return the number of ZCylinders in the cell's surfaces
 * @return the number of ZCylinders used to define this cell's region
 */
int Cell::getNumZCylinders() {

  std::map<int, Halfspace*>::iterator iter;
  std::map<int, Halfspace*> _surfaces = getSurfaces();
  int num = 0;
  for (iter=_surfaces.begin(); iter != _surfaces.end(); ++iter)
    if (iter->second->getSurface()->getSurfaceType() == ZCYLINDER)
      num++;

  return num;
}
