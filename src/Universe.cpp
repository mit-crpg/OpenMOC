#include "Universe.h"


int Universe::_n = 0;

static int auto_id = DEFAULT_INIT_ID;

/**
 * @brief Returns an auto-generated unique Universe ID.
 * @details This method is intended as a utility method for user's writing
 *          OpenMOC input files. The method makes use of a static Universe
 *          ID which is incremented each time the method is called to enable
 *          unique generation of monotonically increasing IDs. The method's
 *          first ID begins at 10000. Hence, user-defined Universe IDs greater
 *          than or equal to 10000 is prohibited.
 */
int universe_id() {
  int id = auto_id;
  auto_id++;
  return id;
}


/**
 * @brief Resets the auto-generated unique Universe ID counter to 10000.
 */
void reset_universe_id() {
  auto_id = DEFAULT_INIT_ID;
}


/**
 * @brief Maximize the auto-generated unique Universe ID counter.
 * @details This method updates the auto-generated unique Universe ID
 *          counter if the input parameter is greater than the present
 *          value. This is useful for the OpenMC compatibility module
 *          to ensure that the auto-generated Universe IDs do not
 *          collide with those created in OpenMC.
 * @param universe_id the id assigned to the auto-generated counter
 */
void maximize_universe_id(int universe_id) {
  if (universe_id > auto_id)
    auto_id = universe_id;
}


/**
 * @brief Constructor assigns a unique and user-specified ID for the Universe.
 * @param id the user-specified optional Universe ID
 * @param name the user-specified optional Universe ID
 */
Universe::Universe(const int id, const char* name) {

  /* If the user did not define an optional ID, create one */
  if (id == -1)
    _id = universe_id();

  /* Use the user-defined ID */
  else
    _id = id;

  _uid = _n;
  _n++;

  _name = NULL;
  setName(name);

  _type = SIMPLE;

  /* By default, the Universe's fissionability is unknown */
  _fissionable = false;
}


/**
 * @brief Destructor clears the Cell pointers container.
 */
Universe::~Universe() {

  if (_name != NULL)
    delete [] _name;

  /* Clear the map of Cells */
  _cells.clear();
}


/**
 * @brief Returns the Universe's unique ID.
 * @return the Universe's unique ID.
 */
int Universe::getUid() const {
  return _uid;
}


/**
 * Return the user-specified ID for this Universe.
 * @return the user-specified Universe ID
 */
int Universe::getId() const {
  return _id;
}


/**
 * @brief Return the user-defined name of the Universe
 * @return the Universe name
 */
char* Universe::getName() const {
  return _name;
}


/**
 * @brief Return the Universe type (SIMPLE or LATTICE).
 * @return the Universe type
 */
universeType Universe::getType() {
  return _type;
}

/**
 * @brief Return the number of Cells in this Universe.
 * @return the number of Cells
 */
int Universe::getNumCells() const {
  return _cells.size();
}


/**
 * @brief Aggregates a list (vector) of the IDs of all Materials within
 *        the MATERIAL type Cells filling this Universe.
 * @details Note that this method only searches the first level of Cells
 *          below this Universe within the nested Universe coordinate system.
 * @return a vector of Material IDs
 */
double Universe::getMinX() {

  double min_x = -std::numeric_limits<double>::infinity();
  std::map<int, Cell*>::iterator c_iter;
  std::map<int, surface_halfspace*>::iterator s_iter;
  Surface* surf;
  int halfspace;

  /* Check if the universe contains a cell with an x-min boundary */
  for (c_iter = _cells.begin(); c_iter != _cells.end(); ++c_iter) {
    std::map<int, surface_halfspace*> surfs = c_iter->second->getSurfaces();

    for (s_iter = surfs.begin(); s_iter != surfs.end(); ++s_iter) {
      surf = s_iter->second->_surface;
      halfspace = s_iter->second->_halfspace;

      if (surf->getSurfaceType() == XPLANE && halfspace == +1 &&
          surf->getBoundaryType() != BOUNDARY_NONE)
        return surf->getMinX(halfspace);
    }
  }

  /* If a x-min boundary was not found, get the x-min from the bounding boxes
   * of the cells */
  for (c_iter = _cells.begin(); c_iter != _cells.end(); ++c_iter)
    min_x = std::min(min_x, c_iter->second->getMinX());

  return min_x;
}


/**
 * @brief Returns the maximum reachable x-coordinate in the Universe.
 * @return the maximum reachable x-coordinate
 */
double Universe::getMaxX() {

  double max_x = std::numeric_limits<double>::infinity();
  std::map<int, Cell*>::iterator c_iter;
  std::map<int, surface_halfspace*>::iterator s_iter;
  Surface* surf;
  int halfspace;

  /* Check if the universe contains a cell with an x-max boundary */
  for (c_iter = _cells.begin(); c_iter != _cells.end(); ++c_iter) {
    std::map<int, surface_halfspace*> surfs = c_iter->second->getSurfaces();

    for (s_iter = surfs.begin(); s_iter != surfs.end(); ++s_iter) {
      surf = s_iter->second->_surface;
      halfspace = s_iter->second->_halfspace;
      if (surf->getSurfaceType() == XPLANE && halfspace == -1 &&
          surf->getBoundaryType() != BOUNDARY_NONE)
        return surf->getMaxX(halfspace);
    }
  }

  /* If a x-max boundary was not found, get the x-max from the bounding boxes
   * of the cells */
  for (c_iter = _cells.begin(); c_iter != _cells.end(); ++c_iter)
    max_x = std::max(max_x, c_iter->second->getMaxX());

  return max_x;
}


/**
 * @brief Returns the minimum reachable y-coordinate in the Universe.
 * @return the minimum reachable y-coordinate
 */
double Universe::getMinY() {

  double min_y = -std::numeric_limits<double>::infinity();
  std::map<int, Cell*>::iterator c_iter;
  std::map<int, surface_halfspace*>::iterator s_iter;
  Surface* surf;
  int halfspace;

  /* Check if the universe contains a cell with an y-min boundary */
  for (c_iter = _cells.begin(); c_iter != _cells.end(); ++c_iter) {
    std::map<int, surface_halfspace*> surfs = c_iter->second->getSurfaces();

    for (s_iter = surfs.begin(); s_iter != surfs.end(); ++s_iter) {
      surf = s_iter->second->_surface;
      halfspace = s_iter->second->_halfspace;

      if (surf->getSurfaceType() == YPLANE && halfspace == +1 &&
          surf->getBoundaryType() != BOUNDARY_NONE)
        return surf->getMinY(halfspace);
    }
  }

  /* If a y-min boundary was not found, get the y-min from the bounding boxes
   * of the cells */
  for (c_iter = _cells.begin(); c_iter != _cells.end(); ++c_iter)
    min_y = std::min(min_y, c_iter->second->getMinY());

  return min_y;
}


/**
 * @brief Returns the maximum reachable y-coordinate in the Universe.
 * @return the maximum reachable y-coordinate
 */
double Universe::getMaxY() {

  double max_y = std::numeric_limits<double>::infinity();
  std::map<int, Cell*>::iterator c_iter;
  std::map<int, surface_halfspace*>::iterator s_iter;
  Surface* surf;
  int halfspace;

  /* Check if the universe contains a cell with an y-max boundary */
  for (c_iter = _cells.begin(); c_iter != _cells.end(); ++c_iter) {
    std::map<int, surface_halfspace*>surfs = c_iter->second->getSurfaces();

    for (s_iter = surfs.begin(); s_iter != surfs.end(); ++s_iter) {
      surf = s_iter->second->_surface;
      halfspace = s_iter->second->_halfspace;
      if (surf->getSurfaceType() == YPLANE && halfspace == -1 &&
          surf->getBoundaryType() != BOUNDARY_NONE)
        return surf->getMaxY(halfspace);
    }
  }

  /* If a y-max boundary was not found, get the y-max from the bounding boxes
   * of the cells */
  for (c_iter = _cells.begin(); c_iter != _cells.end(); ++c_iter)
    max_y = std::max(max_y, c_iter->second->getMaxY());

  return max_y;
}


/**
 * @brief Returns the minimum reachable z-coordinate in the Universe.
 * @return the minimum reachable z-coordinate
 */
double Universe::getMinZ() {

  double min_z = -std::numeric_limits<double>::infinity();
  std::map<int, Cell*>::iterator iter;

  for (iter = _cells.begin(); iter != _cells.end(); ++iter)
    min_z = std::min(min_z, iter->second->getMinZ());

  return min_z;
}


/**
 * @brief Returns the maximum reachable z-coordinate in the Universe.
 * @return the maximum reachable z-coordinate
 */
double Universe::getMaxZ() {

  double max_z = std::numeric_limits<double>::infinity();
  std::map<int, Cell*>::iterator iter;

  for (iter = _cells.begin(); iter != _cells.end(); ++iter)
    max_z = std::max(max_z, iter->second->getMaxZ());

  return max_z;
}


/**
 * @brief Returns the boundary conditions (VACUUM or REFLECTIVE) at the minimum
 *         reachable x-coordinate in the Universe.
 * @return the boundary conditions at the minimum reachable x-coordinate
 */
boundaryType Universe::getMinXBoundaryType() {

  std::map<int, Cell*>::iterator c_iter;
  std::map<int, surface_halfspace*>::iterator s_iter;
  Surface* surf;
  int halfspace;
  boundaryType bc_x = BOUNDARY_NONE;

  /* Check if the universe contains a cell with an x-min boundary */
  for (c_iter = _cells.begin(); c_iter != _cells.end(); ++c_iter) {
    std::map<int, surface_halfspace*>surfs = c_iter->second->getSurfaces();

    for (s_iter = surfs.begin(); s_iter != surfs.end(); ++s_iter) {
      surf = s_iter->second->_surface;
      halfspace = s_iter->second->_halfspace;

      if (surf->getSurfaceType() == XPLANE && halfspace == +1 &&
          surf->getBoundaryType() != BOUNDARY_NONE)
        bc_x = surf->getBoundaryType();
    }
  }

  return bc_x;
}


/**
 * @brief Returns the boundary conditions (VACUUM or REFLECTIVE) at the maximum
 *         reachable x-coordinate in the Universe.
 * @return the boundary conditions at the maximum reachable x-coordinate
 */
boundaryType Universe::getMaxXBoundaryType() {

  std::map<int, Cell*>::iterator c_iter;
  std::map<int, surface_halfspace*>::iterator s_iter;
  Surface* surf;
  int halfspace;
  boundaryType bc_x = BOUNDARY_NONE;

  /* Check if the universe contains a cell with an x-min boundary */
  for (c_iter = _cells.begin(); c_iter != _cells.end(); ++c_iter) {
    std::map<int, surface_halfspace*>surfs = c_iter->second->getSurfaces();

    for (s_iter = surfs.begin(); s_iter != surfs.end(); ++s_iter) {
      surf = s_iter->second->_surface;
      halfspace = s_iter->second->_halfspace;

      if (surf->getSurfaceType() == XPLANE && halfspace == -1 &&
          surf->getBoundaryType() != BOUNDARY_NONE)
        bc_x = surf->getBoundaryType();
    }
  }

  return bc_x;
}


/**
 * @brief Returns the boundary conditions (VACUUM or REFLECTIVE) at the minimum
 *         reachable y-coordinate in the Universe.
 * @return the boundary conditions at the minimum reachable y-coordinate
 */
boundaryType Universe::getMinYBoundaryType() {

  std::map<int, Cell*>::iterator c_iter;
  std::map<int, surface_halfspace*>::iterator s_iter;
  Surface* surf;
  int halfspace;
  boundaryType bc_y = BOUNDARY_NONE;

  /* Check if the universe contains a cell with an x-min boundary */
  for (c_iter = _cells.begin(); c_iter != _cells.end(); ++c_iter) {
    std::map<int, surface_halfspace*>surfs = c_iter->second->getSurfaces();

    for (s_iter = surfs.begin(); s_iter != surfs.end(); ++s_iter) {
      surf = s_iter->second->_surface;
      halfspace = s_iter->second->_halfspace;

      if (surf->getSurfaceType() == YPLANE && halfspace == +1 &&
          surf->getBoundaryType() != BOUNDARY_NONE)
        bc_y = surf->getBoundaryType();
    }
  }

  return bc_y;
}


/**
 * @brief Returns the boundary conditions (VACUUM or REFLECTIVE) at the maximum
 *         reachable y-coordinate in the Universe.
 * @return the boundary conditions at the maximum reachable y-coordinate
 */
boundaryType Universe::getMaxYBoundaryType() {

  std::map<int, Cell*>::iterator c_iter;
  std::map<int, surface_halfspace*>::iterator s_iter;
  Surface* surf;
  int halfspace;
  boundaryType bc_y = BOUNDARY_NONE;

  /* Check if the universe contains a cell with an x-min boundary */
  for (c_iter = _cells.begin(); c_iter != _cells.end(); ++c_iter) {
    std::map<int, surface_halfspace*>surfs = c_iter->second->getSurfaces();

    for (s_iter = surfs.begin(); s_iter != surfs.end(); ++s_iter) {
      surf = s_iter->second->_surface;
      halfspace = s_iter->second->_halfspace;

      if (surf->getSurfaceType() == YPLANE && halfspace == -1 &&
          surf->getBoundaryType() != BOUNDARY_NONE)
        bc_y = surf->getBoundaryType();
    }
  }

  return bc_y;
}


/**
 * @brief Returns a Cell in this universe.
 * @param cell_id the integer the cell_id
 * @return Returns the cell pointer.
 */
Cell* Universe::getCell(int cell_id) {

  if (_cells.find(cell_id) == _cells.end())
    log_printf(ERROR, "Unable to return Cell with ID = %d from Universe with "
               "ID = %d since it does not contain this Cell", cell_id, _id);

    return _cells.at(cell_id);
}


/**
 * @brief Return the container of Cell IDs and Cell pointers in this Universe.
 * @return std::map of Cell IDs
 */
std::map<int, Cell*> Universe::getCells() const {
  return _cells;
}


/**
 * @brief Returns the std::map of Cell IDs and Cell pointers in this Universe
 *        at all nested Universe levels.
 * @return std::map of Cell IDs and pointers
 */
std::map<int, Cell*> Universe::getAllCells() {

  std::map<int, Cell*> cells;
  std::map<int, Cell*>::iterator iter;

  /* Add this Universe's Cells to the map */
  cells.insert(_cells.begin(), _cells.end());

  /* Append all Cells in each Cell in the Universe to the map */
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    std::map<int, Cell*> nested_cells = iter->second->getAllCells();
    cells.insert(nested_cells.begin(), nested_cells.end());
  }

  return cells;
}


/**
 * @brief Returns the std::map of all IDs and Material pointers filling
          this Universe.
 * @return std::map of Material IDs and pointers
 */
std::map<int, Material*> Universe::getAllMaterials() {

  std::map<int, Cell*> cells = getAllCells();
  std::map<int, Cell*>::iterator iter;
  std::map<int, Material*> materials;

  Cell* cell;
  Material* material;

  for (iter = cells.begin(); iter != cells.end(); ++iter) {
    cell = iter->second;

    if (cell->getType() == MATERIAL) {
      material = cell->getFillMaterial();
      materials[material->getId()] = material;
    }
  }

  return materials;
}


/**
 * @brief Returns the std::map of all nested Universe IDs and Universe pointers
          filling this Universe.
 * @return std::map of Universe IDs and pointers
 */
std::map<int, Universe*> Universe::getAllUniverses() {

  /* Get all Cells in this Universe */
  std::map<int, Cell*> cells = getAllCells();

  std::map<int, Universe*> universes;
  universes[_id] = this;
  std::map<int, Cell*>::iterator iter;
  Cell* cell;

  /* Append all Universes containing each Cell to the map */
  for (iter = cells.begin(); iter != cells.end(); ++iter) {
    cell = iter->second;
    std::map<int, Universe*> nested_universes = cell->getAllUniverses();
    universes.insert(nested_universes.begin(), nested_universes.end());
  }

  return universes;
}


/**
 * @brief Returns true if the Universe contains a Cell filled by a fissionable
 *        Material and false otherwise.
 * @details This method should not be called prior to the calling of the
 *          Geometry::computeFissionability() method.
 * @return true if contains a fissionable Material
 */
bool Universe::isFissionable() {
  return _fissionable;
}


/**
 * @brief Sets the name of the Universe
 * @param name the Universe name string
 */
void Universe::setName(const char* name) {
  int length = strlen(name);

  if (_name != NULL)
    delete [] _name;

  /* Initialize a character array for the Universe's name */
  _name = new char[length+1];

  /* Copy the input character array Universe name to the class attribute name */
  for (int i=0; i <= length; i++)
    _name[i] = name[i];
}


/**
 * @brief Sets the Universe type to SIMPLE or LATTICE.
 * @param type the Universe type
 */
void Universe::setType(universeType type) {
  _type = type;
}


/**
 * @brief Sets whether or not this Universe contains a fissionable Material
 *        with a non-zero fission cross-section.
 * @details This method is called by the Geometry::computeFissionability()
 *          class method.
 * @param fissionable true if the Universe contains a fissionable Material;
 *        false otherwise
 */
void Universe::setFissionability(bool fissionable) {
  _fissionable = fissionable;
}


/**
 * @brief Adds a Cell to this Universe.
 * @details Stores the user-specified Cell ID and Cell pointer in a std::map
 *          along with all of other Cells added to this Universe.
 * @param cell the Cell pointer
 */
void Universe::addCell(Cell* cell) {

  try {
    _cells.insert(std::pair<int, Cell*>(cell->getId(), cell));
    log_printf(DEBUG, "Added Cell with ID = %d to Universe with ID = %d",
               cell->getId(), _id);
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to add Cell with ID = %d to Universe with"
               " ID = %d. Backtrace:\n%s", cell, _id, e.what());
  }
}


/**
 * @brief Removes a Cell from this Universe's container of Cells.
 * @param cell a pointer to the Cell to remove
 */
void Universe::removeCell(Cell* cell) {
  if (_cells.find(cell->getId()) != _cells.end())
    _cells.erase(cell->getId());
}


/**
 * @brief Finds the Cell for which a LocalCoords object resides.
 * @details Finds the Cell that a LocalCoords object is located inside by
 *          checking each of this Universe's Cells. Returns NULL if the
 *          LocalCoords is not in any of the Cells.
 * @param coords a pointer to the LocalCoords of interest
 * @return a pointer the Cell where the LocalCoords is located
 */
Cell* Universe::findCell(LocalCoords* coords) {

  Cell* cell;
  std::vector<Cell*> cells;
  std::vector<Cell*>::iterator iter;

  /* Sets the LocalCoord type to UNIV at this level */
  coords->setType(UNIV);

  /* If the LocalCoords is populated with Universe/Cell already, we assume
   * that we are looking for the location in a neighboring Cell */
  if (coords->getCell() != NULL)
    cells = coords->getCell()->getNeighbors();

  /* Add all of Universe's Cells to the back of neighbor Cells vector */
  std::transform(_cells.begin(), _cells.end(),
                 std::back_inserter(cells), pair_second(_cells));

  /* Loop over all Cells */
  for (iter = cells.begin(); iter != cells.end(); ++iter) {
    cell = (*iter);

    if (cell->containsCoords(coords)) {

      /* Set the Cell on this level */
      coords->setCell(cell);

      /* MATERIAL type Cell - lowest level, terminate search for Cell */
      if (cell->getType() == MATERIAL)
        return cell;

      /* FILL type Cell - Cell contains a Universe at a lower level
       * Update coords to next level and continue search */
      else if (cell->getType() == FILL) {

        LocalCoords* next_coords =
            new LocalCoords(coords->getX(), coords->getY(), coords->getZ());
        next_coords->setPhi(coords->getPhi());

        /* Apply translation to position in the next coords */
	if (cell->isTranslated()){
	  double* translation = cell->getTranslation();
	  double new_x = coords->getX() + translation[0];
	  double new_y = coords->getY() + translation[1];
	  double new_z = coords->getZ() + translation[2];
	  next_coords->setX(new_x);
	  next_coords->setY(new_y);
	  next_coords->setZ(new_z);
	}

        /* Apply rotation to position and direction in the next coords */
	if (cell->isRotated()){
	  double x = coords->getX();
	  double y = coords->getY();
	  double z = coords->getZ();
	  double* matrix = cell->getRotationMatrix();
	  double new_x = matrix[0] * x + matrix[1] * y + matrix[2] * z;
	  double new_y = matrix[3] * x + matrix[4] * y + matrix[5] * z;
	  double new_z = matrix[6] * x + matrix[7] * y + matrix[8] * z;
	  next_coords->setX(new_x);
	  next_coords->setY(new_y);
	  next_coords->setZ(new_z);
          next_coords->incrementPhi(cell->getPsi() * M_PI / 180.);
	}

        Universe* univ = cell->getFillUniverse();
        next_coords->setUniverse(univ);

        coords->setNext(next_coords);
        next_coords->setPrev(coords);
        if (univ->getType() == SIMPLE)
          return univ->findCell(next_coords);
        else
          return static_cast<Lattice*>(univ)->findCell(next_coords);
      }
    }
  }

  return NULL;
}


/**
 * @brief Subdivides all of the Material-filled Cells within this Universe
 *        into rings and angular sectors aligned with the z-axis.
 * @param max_radius the maximum allowable radius used in the subdivisions
 */
void Universe::subdivideCells(double max_radius) {

  log_printf(DEBUG, "Subdividing Cells for Universe ID=%d "
	     "with max radius %f", _id, max_radius);

  std::map<int, Cell*>::iterator iter;

  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {

    /* Cells filled with Materials */
    if (iter->second->getType() == MATERIAL) {
      Cell* cell = iter->second;

      if (cell->getNumRings() > 0 || cell->getNumSectors() > 0)
        cell->subdivideCell(max_radius);
    }

    /* Cells filled with Universes */
    else {
      Universe* fill = iter->second->getFillUniverse();
      if (fill->getType() == SIMPLE)
        fill->subdivideCells(max_radius);
      else
        static_cast<Lattice*>(fill)->subdivideCells(max_radius);
    }
  }
}


/**
 * @brief Builds collections of neighboring Cells for all Cells in this
 *        Universe for optimized ray tracing.
 */
void Universe::buildNeighbors() {

  /* Loop over all of the Universe's Cells and make recursive call */
  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter)
    iter->second->buildNeighbors();
}


/**
 * @brief Convert the member attributes of this Universe to a character array.
 * @return a character array representing the Universe's attributes
 */
std::string Universe::toString() {

  std::stringstream string;
  std::map<int, Cell*>::iterator iter;

  string << "Universe ID = " << _id;
  string << ", name = " << _name;

  string << ", type = ";
  if (_type == SIMPLE)
    string << "SIMPLE";
  else
    string << "LATTICE";

  string << ", # cells = " << _cells.size() << ", cell IDs = ";

  for (iter = _cells.begin(); iter != _cells.end(); ++iter)
    string << iter->first << ", ";

  return string.str();
}


/**
 * @brief Prints a string representation of the Universe's attributes to
 *        the console.
 */
void Universe::printString() {
  log_printf(RESULT, toString().c_str());
}


/**
 * @brief Clones this Universe and copy cells map
 * @return a pointer to the Universe clone
 */
Universe* Universe::clone() {

  log_printf(DEBUG, "Cloning Universe %d", _id);

  /* Instantiate new Universe clone */
  Universe* clone = new Universe(universe_id(), _name);

  /* Loop over Cells in Universe and clone each one */
  std::map<int, Cell*>::iterator iter1;
  for (iter1 = _cells.begin(); iter1 != _cells.end(); ++iter1) {

    /* If the Cell is filled with a Material, clone it */
    if (iter1->second->getType() == MATERIAL) {

      /* Clone the Cell */
      Cell* parent = static_cast<Cell*>(iter1->second);
      Cell* cell_clone = parent->clone();

      /* Add Cell clone to the list */
      clone->addCell(cell_clone);
    }
    /* Throw error message if Cell is FILL type */
    else {
      log_printf(ERROR, "Unable to clone Universe %d since it contains Cell %d"
                 "which is filled with a Universe rather than a Material");
    }
  }

  return clone;
}

/**
 * @brief Constructor sets the user-specified and unique IDs for this Lattice.
 * @param id the user-specified optional Lattice (Universe) ID
 * @param name the user-specified optional Lattice (Universe) name
 */
Lattice::Lattice(const int id, const char* name): Universe(id, name) {

  _type = LATTICE;
  _offset.setCoords(0.0, 0.0, 0.0);

  /* Default width and number of Lattice cells along each dimension */
  _num_x = -1;
  _num_y = -1;
  _num_z = -1;
  _width_x = 0;
  _width_y = 0;
  _width_z = 0;
}


/**
 * @brief Destructor clears memory for all of Universes pointers.
 */
Lattice::~Lattice() {

  /* Clear the triple-nested vector of Universes */
  for (int k=0; k < _universes.size(); k++) {
    for (int j=0; j < _universes.at(k).size(); j++)
      _universes.at(k).at(j).clear();
    _universes.at(k).clear();
  }

  _universes.clear();
}


/**
 * @brief Set the offset in global coordinates for this Lattice.
 * @details A lattice is assumed to be a rectilinear grid with the center/origin
 *          of the grid located in the center of the Lattice's parent universe.
 *          The offset represents the offset of the lattice center/origin with
 *          respect to the center of the parent universe. Therefore an offset of
 *          (-1,2) would move the center/origin of the lattice to the left 1 cm
 *          and up 2 cm.
 * @param x the offset in the x direction
 * @param y the offset in the y direction
 * @param z the offset in the z direction
 */
void Lattice::setOffset(double x, double y, double z) {
  _offset.setX(x);
  _offset.setY(y);
  _offset.setZ(z);
}


/**
 * @brief Return a pointer to the offset for this Cell (in global coordinates).
 * @return the offset of the Cell
 */
Point* Lattice::getOffset() {
  return &_offset;
}


/**
 * @brief Return the number of Lattice cells along the x-axis.
 * @return the number of Lattice cells along x
 */
int Lattice::getNumX() const {
  return _num_x;
}


/**
 * @brief Return the number of Lattice cells along the y-axis
 * @return the number of Lattice cells along y
 */
int Lattice::getNumY() const {
  return _num_y;
}


/**
 * @brief Return the number of Lattice cells along the z-axis
 * @return the number of Lattice cells along z
 */
int Lattice::getNumZ() const {
  return _num_z;
}


/**
 * @brief Return the width of the Lattice along the x-axis.
 * @return the width of the Lattice cells along x
 */
double Lattice::getWidthX() const {
  return _width_x;
}


/**
 * @brief Return the width of the Lattice along the y-axis.
 * @return the width of the Lattice cells along y
 */
double Lattice::getWidthY() const {
  return _width_y;
}


/**
 * @brief Return the width of the Lattice along the z-axis.
 * @return the width of the Lattice cells along z
 */
double Lattice::getWidthZ() const {
  return _width_z;
}


/**
 * @brief Returns the minimum reachable x-coordinate in the Lattice.
 * @return the minimum reachable x-coordinate
 */
double Lattice::getMinX() {
  return _offset.getX() - (_num_x * _width_x / 2.);
}


/**
 * @brief Returns the maximum reachable x-coordinate in the Lattice.
 * @return the maximum reachable x-coordinate
 */
double Lattice::getMaxX() {
  return _offset.getX() + (_num_x * _width_x / 2.);
}


/**
 * @brief Returns the minimum reachable y-coordinate in the Lattice.
 * @return the minimum reachable y-coordinate
 */
double Lattice::getMinY() {
  return _offset.getY() - (_num_y * _width_y / 2.);
}


/**
 * @brief Returns the maximum reachable y-coordinate in the Lattice.
 * @return the maximum reachable y-coordinate
 */
double Lattice::getMaxY() {
  return _offset.getY() + (_num_y * _width_y / 2.);
}


/**
 * @brief Returns the minimum reachable z-coordinate in the Lattice.
 * @return the minimum reachable z-coordinate
 */
double Lattice::getMinZ() {
  return _offset.getZ() - (_num_z * _width_z / 2.);
}


/**
 * @brief Returns the maximum reachable z-coordinate in the Lattice.
 * @return the maximum reachable z-coordinate
 */
double Lattice::getMaxZ() {
  return _offset.getZ() + (_num_z * _width_z / 2.);
}


/**
 * @brief Returns a pointer to the Universe within a specific Lattice cell.
 * @param lat_x the x index to the Lattice cell
 * @param lat_y the y index to the Lattice cell
 * @param lat_z the z index to the Lattice cell
 * @return pointer to a Universe filling the Lattice cell
 */
Universe* Lattice::getUniverse(int lat_x, int lat_y, int lat_z) const {

  /* Checks that lattice indices are within the bounds of the lattice */
  if (lat_x >= _num_x || lat_y >= _num_y || lat_z >= _num_z ||
      lat_x < 0 || lat_y < 0 || lat_z < 0)
    log_printf(ERROR, "Cannot retrieve Universe from Lattice ID = %d: Index"
               "out of bounds: Tried to access Cell x = %d, y = %d, "
               "z = %d but bounds are x = %d, y = %d, z = %d",
               _id, lat_x, lat_y, lat_z, _num_x, _num_y, _num_z);

  return _universes.at(lat_z).at(lat_y).at(lat_x).second;
}


/**
 * @brief Return a 3D vector of the Universes in the Lattice.
 * @return 3D vector of Universes
 */
std::vector< std::vector< std::vector< std::pair<int, Universe*> > > >*
  Lattice::getUniverses() {
  return &_universes;
}


/**
 * @brief Aggregates a list (vector) of the IDs of all Universes within
 *        the FILL type Cells filling this Universe.
 * @details Note that this method only searches the first level of Cells
 *          below this Universe within the nested Universe coordinate system.
 * @return a vector of Universe IDs
 */
std::map<int, Universe*> Lattice::getUniqueUniverses() {

  std::map<int, Universe*> unique_universes;
  Universe* universe;
  int univ_id;

  for (int k = _num_z-1; k > -1; k--) {
    for (int j = _num_y-1; j > -1;  j--) {
      for (int i = 0; i < _num_x; i++) {
        universe = _universes.at(k).at(j).at(i).second;
        unique_universes[universe->getId()] = universe;
      }
    }
  }

  return unique_universes;
}


/**
 * @brief Returns the std::map of Cell IDs and Cell pointers in this Lattice
 *        at all nested Universe levels.
 * @return std::map of Cell IDs and pointers
 */
std::map<int, Cell*> Lattice::getAllCells() {

  std::map<int, Cell*> cells;
  std::map<int, Universe*> unique_universes = getUniqueUniverses();

  std::map<int, Universe*>::iterator iter;
  std::map<int, Cell*> nested_cells;

  for (iter = unique_universes.begin(); iter != unique_universes.end(); ++iter) {
    nested_cells = iter->second->getAllCells();
    cells.insert(nested_cells.begin(), nested_cells.end());
  }

  return cells;
}


/**
 * @brief Returns the std::map of all nested Universe IDs and Universe pointers
          filling this Lattice.
 * @return std::map of Universe IDs and pointers
 */
std::map<int, Universe*> Lattice::getAllUniverses() {

  /* Initialize a map of all Universes contained by the Lattice in each
   * nested Universe level */
  std::map<int, Universe*> all_universes;

  /* Get all unique Universes contained in each of the Lattice cells */
  std::map<int, Universe*> unique_universes = getUniqueUniverses();

  /* Add the unique Universes filling each Lattice cell */
  all_universes.insert(unique_universes.begin(), unique_universes.end());

  /* Append all Universes containing each Cell to the map */
  std::map<int, Universe*>::iterator iter;
  std::map<int, Universe*> nested_universes;

  for (iter = unique_universes.begin(); iter != unique_universes.end(); ++iter) {
    nested_universes = iter->second->getAllUniverses();
    all_universes.insert(nested_universes.begin(), nested_universes.end());
  }

  return all_universes;
}


/**
 * @brief Set the number of Lattice cells along the x-axis.
 * @param num_x the number of Lattice cells along x
 */
void Lattice::setNumX(int num_x) {
  _num_x = num_x;
}


/**
 * @brief Set the number of Lattice cells along the y-axis.
 * @param num_y the number of Lattice cells along y
 */
void Lattice::setNumY(int num_y) {
  _num_y = num_y;
}


/**
 * @brief Set the number of Lattice cells along the z-axis.
 * @param num_z the number of Lattice cells along z
 */
void Lattice::setNumZ(int num_z) {
  _num_z = num_z;
}


/**
 * @brief Set the width of each Lattice cell.
 * @param width_x the width along the x-axis in centimeters
 * @param width_y the width along the y-axis in centimeters
 * @param width_z the width along the z-axis in centimeters
 */
void Lattice::setWidth(double width_x, double width_y, double width_z) {

  if (width_x <= 0 || width_y <= 0 || width_z <= 0)
    log_printf(ERROR, "Unable to set the width of Lattice ID = %d "
               "for x = %f, y = %f, and z = %f since they are not positive "
               "values", _id, width_x, width_y, _width_z);

  _width_x = width_x;
  _width_y = width_y;
  _width_z = width_z;
}


/**
 * @brief Sets the array of Universe pointers filling each Lattice cell.
 * @details This is a helper method for SWIG to allow users to assign Universes
 *          to a Lattice using a 3D Python list (list of lists of lists). An
 *          example how this method can be called from Python is as follows:
 *
 * @code
 *          u1 = Universe(name='Universe 1')
 *          u2 = Universe(name='Universe 2')
 *          u3 = Universe(name='Universe 3')
 *          lattice.setUniverses([[[u1, u2, u1, u2],
 *                                 [u2, u3, u2, u3],
 *                                 [u1, u2, u1, u2],
 *                                 [u2, u3, u2, u3]],
 *                                [[u1, u2, u1, u2],
 *                                 [u2, u3, u2, u3],
 *                                 [u1, u2, u1, u2],
 *                                 [u2, u3, u2, u3]]])
 * @endcode
 *
 * @param num_z the number of Lattice cells along z
 * @param num_y the number of Lattice cells along y
 * @param num_x the number of Lattice cells along x
 * @param universes the array of Universes for each Lattice cell
 */
void Lattice::setUniverses(int num_z, int num_y, int num_x,
                           Universe** universes) {

  std::map<int, Universe*> unique_universes = getUniqueUniverses();
  std::map<int, Universe*>::iterator iter;

  /* Remove all Universes in the Lattice */
  for (iter = unique_universes.begin(); iter != unique_universes.end(); ++iter)
    removeUniverse(iter->second);

  /* Clear all Universe maps in the Lattice (from a previous run) */
  for (int i=0; i < _num_y; i++)
    _universes.at(i).clear();

  _universes.clear();

  /* Set the Lattice dimensions */
  setNumX(num_x);
  setNumY(num_y);
  setNumZ(num_z);

  Universe* universe;

  /* The Lattice cells are assumed input in row major order starting from the
   * upper left corner. This double loop reorders the Lattice cells from the
   * start from the lower left corner */
  for (int k = 0; k < _num_z; k++) {
    _universes.push_back(
        std::vector< std::vector< std::pair<int, Universe*> > >());
    for (int j = 0; j < _num_y; j++) {

      _universes.at(k).push_back(std::vector< std::pair<int, Universe*> >());

      for (int i = 0; i < _num_x; i++) {
        universe = universes
            [(_num_z-1-k)*_num_x*_num_y + (_num_y-1-j)*_num_x + i];
        _universes.at(k).at(j).push_back(std::pair<int, Universe*>(
            universe->getId(), universe));
      }
    }
  }
}


/**
 * @brief Update the Universe in a particular Lattice cell.
 * @details This method may only be used after an array of Universes
 *          has been assigned with the Lattice::setUniverses(...) method.
 * @param lat_x the Lattice cell index along x
 * @param lat_y the Lattice cell index along y
 * @param lat_z the Lattice cell index along z
 * @param universe the Universe to insert into the Lattice
 */
void Lattice::updateUniverse(int lat_x, int lat_y, int lat_z,
                             Universe* universe) {

  if (_num_x == -1 || _num_y == -1 || _num_z == -1)
    log_printf(ERROR, "Unable to set Universe %d in Lattice %d which "
	       "has not yet been assigned an array of Universes",
	       universe->getId(), _id);
  if (lat_x < 0 || lat_x >= _num_x)
    log_printf(ERROR, "Unable to set Universe %d in Lattice %d with "
	       "Lattice cell index lat_x=%d which is outside the "
	       "array of Universes", universe->getId(), _id, lat_x);
  if (lat_y < 0 || lat_y >= _num_y)
    log_printf(ERROR, "Unable to set Universe %d in Lattice %d with "
	       "Lattice cell index lat_z=%d which is outside the "
	       "array of Universes", universe->getId(), _id, lat_z);
  if (lat_z < 0 || lat_z >= _num_z)
    log_printf(ERROR, "Unable to set Universe %d in Lattice %d with "
	       "Lattice cell index lat_z=%d which is outside the "
	       "array of Universes", universe->getId(), _id, lat_z);

  /* Assign the Universe to the array */
  _universes.at(lat_z).at(lat_y).at(lat_x) =
    std::pair<int, Universe*>(universe->getId(), universe);
}


/**
 * @brief Removes all references to a Universe from the Lattice.
 * @param universe the Universe to remove
 */
void Lattice::removeUniverse(Universe* universe) {

  Universe* null = NULL;

  /* Clear any Universes in the Lattice (from a previous run) */
  for (int k=0; k < _num_z; k++) {
    for (int j=0; j < _num_y; j++) {
      for (int i=0; i < _num_x; i++) {
        if (universe->getId() == getUniverse(i,j,k)->getId())
          _universes.at(k).at(j).at(i) = std::pair<int,Universe*>(-1, null);
      }
    }
  }
}


/**
 * @brief Subdivides all of the Material-filled Cells within this Lattice
 *        into rings and angular sectors aligned with the z-axis.
 * @param max_radius the maximum allowable radius used in the subdivisions
 */
void Lattice::subdivideCells(double max_radius) {

  log_printf(DEBUG, "Subdividing Cells for Lattice ID=%d "
             "with max radius %f", _id, max_radius);

  std::map<int, Universe*>::iterator iter;
  std::map<int, Universe*> universes = getUniqueUniverses();

  /* Compute the max radius as the distance from the center to a corner
   * of a Lattice cell. */
  double radius = sqrt(_width_x*_width_x/4.0 + _width_y*_width_y/4.0);

  /* If the lattice pitch is smaller than max_radius parameter, over-ride it */
  if (radius < max_radius)
    max_radius = radius;

  /* Subdivide all Cells */
  for (iter = universes.begin(); iter != universes.end(); ++iter)
    iter->second->subdivideCells(max_radius);
}


/**
 * @brief Builds collections of neighboring Cells for all Cells in each
 *        Universe in the Lattice for optimized ray tracing.
 */
void Lattice::buildNeighbors() {

  /* Get list of unique Universes in this Lattice */
  std::map<int, Universe*> universes = getUniqueUniverses();

  /* Loop over each Universe and make recursive call */
  std::map<int, Universe*>::iterator iter;
  for (iter = universes.begin(); iter != universes.end(); ++iter)
    iter->second->buildNeighbors();
}



/**
 * @brief Checks if a Point is within the bounds of a Lattice.
 * @param point a pointer to the Point of interest
 * @return true if the Point is in the bounds, false if not
 */
bool Lattice::withinBounds(Point* point) {

  /* Computes the Lattice bounds */
  double bound_x_max = _offset.getX() + _num_x/2.0 * _width_x;
  double bound_x_min = _offset.getX() - _num_x/2.0 * _width_x;
  double bound_y_max = _offset.getY() + _num_y/2.0 * _width_y;
  double bound_y_min = _offset.getY() - _num_y/2.0 * _width_y;
  double bound_z_max = _offset.getZ() + _num_z/2.0 * _width_z;
  double bound_z_min = _offset.getZ() - _num_z/2.0 * _width_z;

  double x = point->getX();
  double y = point->getY();
  double z = point->getZ();

  /* If the Point is outside the x bounds */
  if (x > bound_x_max || x < bound_x_min)
    return false;

  /* If the Point is outside the y bounds */
  else if (y > bound_y_max || y < bound_y_min)
    return false;

  /* If the Point is outside the z bounds */
  else if (z > bound_z_max || z < bound_z_min)
    return false;

  /* If the Point is within the bounds */
  else
    return true;
}


/**
 * @brief Finds the Cell within this Lattice that a LocalCoords is in.
 * @details This method first find the Lattice cell, then searches the
 *          Universe inside that Lattice cell. If LocalCoords is outside
 *          the bounds of the Lattice, this method will return NULL.
 * @param coords the LocalCoords of interest
 * @return a pointer to the Cell this LocalCoord is in or NULL
 */
Cell* Lattice::findCell(LocalCoords* coords) {

  /* Set the LocalCoord to be a LAT type at this level */
  coords->setType(LAT);

  /* Compute the x, y, and z indices for the Lattice cell this coord is in */
  int lat_x = getLatX(coords->getPoint());
  int lat_y = getLatY(coords->getPoint());
  int lat_z = getLatZ(coords->getPoint());
  double next_x, next_y, next_z;

  /* If the indices are outside the bound of the Lattice */
  if (lat_x < 0 || lat_x >= _num_x ||
      lat_y < 0 || lat_y >= _num_y ||
      lat_z < 0 || lat_z >= _num_z) {
    return NULL;
  }

  /* Compute local position of Point in the next level Universe */
  next_x = coords->getX() -
    (-_width_x*_num_x/2.0 + _offset.getX() + (lat_x + 0.5) * _width_x) +
    getOffset()->getX();
  next_y = coords->getY() -
    (-_width_y*_num_y/2.0 + _offset.getY() + (lat_y + 0.5) * _width_y) +
    getOffset()->getY();

  if (_width_z == std::numeric_limits<double>::infinity())
    next_z = coords->getZ();
  else
    next_z = coords->getZ() -
      (-_width_z*_num_z/2.0 + _offset.getZ() + (lat_z + 0.5) * _width_z) +
      getOffset()->getZ();

  /* Create a new LocalCoords object for the next level Universe */
  LocalCoords* next_coords;

  if (coords->getNext() == NULL)
    next_coords = new LocalCoords(next_x, next_y, next_z);
  else
    next_coords = coords->getNext();

  Universe* univ = getUniverse(lat_x, lat_y, lat_z);
  next_coords->setUniverse(univ);
  next_coords->setPhi(coords->getPhi());

  /* Set Lattice indices */
  coords->setLattice(this);
  coords->setLatticeX(lat_x);
  coords->setLatticeY(lat_y);
  coords->setLatticeZ(lat_z);

  coords->setNext(next_coords);
  next_coords->setPrev(coords);

  /* Search the next lowest level Universe for the Cell */
  return univ->findCell(next_coords);
}


/**
 * @brief Finds the distance to the nearest surface.
 * @details Knowing that a Lattice must be cartesian, this function computes
 *          the distance to the nearest boundary between lattice cells
 *          in the direction of the track.
 * @param coords a pointer to a localcoords object
 * @return the distance to the nearest Lattice cell boundary
 */
double Lattice::minSurfaceDist(LocalCoords* coords) {

  /* Compute the x, y, and z indices for the Lattice cell this point is in */
  int lat_x = getLatX(coords->getPoint());
  int lat_y = getLatY(coords->getPoint());
  double phi = coords->getPhi();

  /* Create planes representing the boundaries of the lattice cell */
  XPlane xplane(0.0);
  YPlane yplane(0.0);

  /* Get the min distance for X PLANE  */
  if (phi < M_PI_2)
    xplane.setX(((lat_x+1) * _width_x - _width_x*_num_x/2.0 + _offset.getX()));
  else
    xplane.setX((lat_x * _width_x - _width_x*_num_x/2.0 + _offset.getX()));

  double dist_x = xplane.getMinDistance(coords);

  /* Get the min distance for Y PLANE */
  if (phi < M_PI)
    yplane.setY(((lat_y+1) * _width_y - _width_y*_num_y/2.0 + _offset.getY()));
  else
    yplane.setY((lat_y * _width_y - _width_y*_num_y/2.0 + _offset.getY()));

  double dist_y = yplane.getMinDistance(coords);

  /* return shortest distance to next lattice cell */
  return std::min(dist_x, dist_y);
}


/**
 * @brief Finds the Lattice cell x index that a point lies in.
 * @param point a pointer to a point being evaluated.
 * @return the Lattice cell x index.
 */
int Lattice::getLatX(Point* point) {

  /* Compute the x indice for the Lattice cell this point is in */
  int lat_x = (int)floor((point->getX() + _width_x*_num_x/2.0 -
                          _offset.getX()) / _width_x);

  /* get the distance to the left surface */
  double dist_to_left = point->getX() + _num_x*_width_x/2.0 - _offset.getX();

  /* Check if the Point is on the Lattice boundaries and if so adjust
   * x Lattice cell indice */
  if (fabs(dist_to_left) < ON_SURFACE_THRESH)
    lat_x = 0;
  else if (fabs(dist_to_left - _num_x*_width_x) < ON_SURFACE_THRESH)
    lat_x = _num_x - 1;
  else if (lat_x < 0 || lat_x > _num_x-1)
    log_printf(ERROR, "Trying to get lattice x index for point that is "
               "outside lattice bounds.");

  return lat_x;
}


/**
 * @brief Finds the Lattice cell y index that a point lies in.
 * @param point a pointer to a point being evaluated.
 * @return the Lattice cell y index.
 */
int Lattice::getLatY(Point* point) {

  /* Compute the y indice for the Lattice cell this point is in */
  int lat_y = (int)floor((point->getY() + _width_y*_num_y/2.0 -
                          _offset.getY()) / _width_y);

  /* get the distance to the bottom surface */
  double dist_to_bottom = point->getY() + _width_y*_num_y/2.0 - _offset.getY();

  /* Check if the Point is on the Lattice boundaries and if so adjust
   * y Lattice cell indice */
  if (fabs(dist_to_bottom) < ON_SURFACE_THRESH)
    lat_y = 0;
  else if (fabs(dist_to_bottom - _num_y*_width_y) < ON_SURFACE_THRESH)
    lat_y = _num_y - 1;
  else if (lat_y < 0 || lat_y > _num_y-1)
    log_printf(ERROR, "Trying to get lattice y index for point that is "
               "outside lattice bounds.");

  return lat_y;
}


/**
 * @brief Finds the Lattice cell z index that a point lies in.
 * @param point a pointer to a point being evaluated.
 * @return the Lattice cell z index.
 */
int Lattice::getLatZ(Point* point) {

  /* Check to see if lattice is infinite in z direction */
  if (_width_z == std::numeric_limits<double>::infinity())
    return 0;

  /* Compute the z indice for the Lattice cell this point is in */
  int lat_z = (int)floor((point->getZ() + _width_z*_num_z/2.0 -
                          _offset.getZ()) / _width_z);

  /* get the distance to the bottom surface */
  double dist_to_bottom = point->getZ() + _width_z*_num_z/2.0 - _offset.getZ();

  /* Check if the Point is on the Lattice boundaries and if so adjust
   * z Lattice cell indice */
  if (fabs(dist_to_bottom) < ON_SURFACE_THRESH)
    lat_z = 0;
  else if (fabs(dist_to_bottom - _num_z*_width_z) < ON_SURFACE_THRESH)
    lat_z = _num_z - 1;
  else if (lat_z < 0 || lat_z > _num_z-1)
    log_printf(ERROR, "Trying to get lattice z index for point that is "
               "outside lattice bounds: %f.", point->getZ());

  return lat_z;
}


/**
 * @brief Converts a Lattice's attributes to a character array representation.
 * @return character array of this Lattice's attributes
 */
std::string Lattice::toString() {

  std::stringstream string;

  string << "Lattice ID = " << _id
         << ", name = " << _name
         << ", # cells along x = " << _num_x
         << ", # cells along y = " << _num_y
         << ", # cells along z = " << _num_z
         << ", x width = " << _width_x
         << ", y width = " << _width_y
         << ", z width = " << _width_z;

  string << "\n\t\tUniverse IDs within this Lattice: ";

  for (int k = _num_z-1; k > -1;  k--) {
    for (int j = _num_y-1; j > -1;  j--) {
      for (int i = 0; i < _num_x; i++)
        string << _universes.at(k).at(j).at(i).first << ", ";
      string << "\n\t\t";
    }
  }

  return string.str();
}


/**
 * @brief Prints a string representation of all of the Lattice's attributes to
 *        the console.
 */
void Lattice::printString() {
  log_printf(RESULT, toString().c_str());
}


/**
 * @brief Finds the Lattice cell index that a point lies in.
 * @details Lattice cells are numbered starting with 0 x-min/y-min/z-min corner.
 *          Lattice cell IDs then increase monotonically from x-min to x-max,
 *          y-min to y-max, and z-min to z-max. Note that values increase first
 *          on the x-axis, followed by the y-axis, then on the z-axis.
 *          For example, the indices for a 4 x 4 x 1 lattice:
 *                  12  13  14  15
 *                  8    9  10  11
 *                  4    5   6   7
 *                  0    1   2   3
 * @param point a pointer to a point being evaluated.
 * @return the Lattice cell index.
 */
int Lattice::getLatticeCell(Point* point) {
  return getLatZ(point)*_num_x*_num_y + getLatY(point)*_num_x + getLatX(point);
}


/**
 * @brief Finds the distance from a point to a particular lattice cell surface.
 * @param cell the cell index that the point is in.
 * @param point a pointer to a point being evaluated.
 * @param surface a surface id to get the distance to.
 * @return the distance to the lattice cell surface of interest.
 */
double Lattice::getDistanceToSurface(int cell, Point* point, int surface) {

  /* Get coordinates of the point */
  double x = point->getX();
  double y = point->getY();

  /* Check that the surface is valid */
  if (surface < 0 || surface >= NUM_SURFACES)
    log_printf(ERROR, "Unable to get the distance from point (%f, %f) to "
               "lattice surface %d since there are only %d surfaces",
               x, y, surface, NUM_SURFACES);

  /* Get lattice indices of the cell */
  int lat_x = cell % _num_x;
  int lat_y = cell / _num_x;

  double dist;

  /* Return the distance to the corresponding lattice cell surface */
  if (surface == SURFACE_X_MIN) {
    dist = lat_x*_width_x - _width_x*_num_x/2.0 + _offset.getX();
    return fabs(x - dist);
  }
  else if (surface == SURFACE_X_MAX) {
    dist = (lat_x + 1)*_width_x - _width_x*_num_x/2.0 + _offset.getX();
    return fabs(x - dist);
  }
  else if (surface == SURFACE_Y_MIN) {
    dist = lat_y*_width_y - _width_y*_num_y/2.0 + _offset.getY();
    return fabs(y - dist);
  }
  else {
    dist = (lat_y + 1)*_width_y - _width_y*_num_y/2.0 + _offset.getY();
    return fabs(y - dist);
  }
}


/**
 * @brief Finds the Lattice cell surface that a point lies on.
 *        If the point is not on a surface, -1 is returned.
 * @details The surface indices are defined in constants.h as they
 *          need to be consistent with the surface constant definitions
 *          used in Cmfd. The index returned takes into account
 *         the cell index and returns NUM_SURFACES*cell_index + surface_index.
 * @param cell the cell index that the point is in.
 * @param point a pointer to a point being evaluated.
 * @return the Lattice surface index.
 */
int Lattice::getLatticeSurface(int cell, Point* point) {

  int surface = -1;

  /* Get coordinates of point and cell boundaries */
  double x = point->getX();
  double y = point->getY();
  double z = point->getZ();
  int lat_x = (cell % (_num_x*_num_y)) % _num_x;
  int lat_y = (cell % (_num_x*_num_y)) / _num_x;
  int lat_z = cell / (_num_x*_num_y);

  /* Create planes representing the boundaries of the lattice cell */
  XPlane xplane(0.0);
  YPlane yplane(0.0);

  /* Bools indicating if point is on each surface */
  bool on_min_x, on_max_x, on_min_y, on_max_y;

  /* Check if point is on X_MIN boundary */
  xplane.setX((lat_x*_width_x - _width_x*_num_x/2.0 + _offset.getX()));
  on_min_x = xplane.isPointOnSurface(point);

  /* Check if point is on X_MAX boundary */
  xplane.setX(((lat_x + 1)*_width_x - _width_x*_num_x/2.0 + _offset.getX()));
  on_max_x = xplane.isPointOnSurface(point);

  /* Check if point is on Y_MIN boundary */
  yplane.setY((lat_y*_width_y - _width_y*_num_y/2.0 + _offset.getY()));
  on_min_y = yplane.isPointOnSurface(point);

  /* Check if point is on Y_MAX boundary */
  yplane.setY(((lat_y + 1)*_width_y - _width_y*_num_y/2.0 + _offset.getY()));
  on_max_y = yplane.isPointOnSurface(point);

  if (on_min_x) {
    if (on_min_y)
      surface = SURFACE_X_MIN_Y_MIN;
    else if (on_max_y)
      surface = SURFACE_X_MIN_Y_MAX;
    else
      surface = SURFACE_X_MIN;
  }
  else if (on_max_x) {
    if (on_min_y)
      surface = SURFACE_X_MAX_Y_MIN;
    else if (on_max_y)
      surface = SURFACE_X_MAX_Y_MAX;
    else
      surface = SURFACE_X_MAX;
  }
  else if (on_min_y)
    surface = SURFACE_Y_MIN;
  else if (on_max_y)
    surface = SURFACE_Y_MAX;

  if (surface != -1)
    surface = NUM_SURFACES * cell + surface;

  return surface;
}
