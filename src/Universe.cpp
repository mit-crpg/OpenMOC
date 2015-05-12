#include "Universe.h"


int Universe::_n = 0;

static int auto_id = 10000;


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
  auto_id = 10000;
}


/**
 * @brief Constructor assigns a unique and user-specified ID for the Universe.
 * @param id the user-specified optional Universe ID
 * @param name the user-specified optional Universe ID
 */
Universe::Universe(const int id, const char* name) {

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

  _type = SIMPLE;

  /* By default, the Universe's fissionability is unknown */
  _fissionable = false;
}


/**
 * @brief Destructor clears the Cell pointers container.
 */
Universe::~Universe() {
  _cells.clear();

  if (_name != NULL)
    delete [] _name;
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

  double min_x = std::numeric_limits<double>::infinity();

  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {

    if (min_x > iter->second->getMinX())
      min_x = iter->second->getMinX();
  }

  return min_x;
}


/**
 * @brief Returns the maximum reachable x-coordinate in the Universe.
 * @return the maximum reachable x-coordinate
 */
double Universe::getMaxX() {

  double max_x = -std::numeric_limits<double>::infinity();

  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    if (max_x < iter->second->getMaxX())
      max_x = iter->second->getMaxX();
  }

  return max_x;
}


/**
 * @brief Returns the minimum reachable y-coordinate in the Universe.
 * @return the minimum reachable y-coordinate
 */
double Universe::getMinY() {

  double min_y = std::numeric_limits<double>::infinity();

  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    if (min_y > iter->second->getMinY())
      min_y = iter->second->getMinY();
  }

  return min_y;
}


/**
 * @brief Returns the maximum reachable y-coordinate in the Universe.
 * @return the maximum reachable y-coordinate
 */
double Universe::getMaxY(){

  double max_y = -std::numeric_limits<double>::infinity();

  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    if (max_y < iter->second->getMaxY())
      max_y = iter->second->getMaxY();
  }

  return max_y;
}


/**
 * @brief Returns the minimum reachable z-coordinate in the Universe.
 * @return the minimum reachable z-coordinate
 */
double Universe::getMinZ() {
  double min_z = std::numeric_limits<double>::infinity();

  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    if (min_z > iter->second->getMinZ())
      min_z = iter->second->getMinZ();
  }

  return min_z;
}


/**
 * @brief Returns the maximum reachable z-coordinate in the Universe.
 * @return the maximum reachable z-coordinate
 */
double Universe::getMaxZ() {

  double max_z = -std::numeric_limits<double>::infinity();

  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    if (max_z < iter->second->getMaxZ())
      max_z = iter->second->getMaxZ();
  }

  return max_z;
}


/**
 * @brief Returns the boundary conditions (VACUUM or REFLECTIVE) at the minimum
 *         reachable x-coordinate in the Universe.
 * @return the boundary conditions at the minimum reachable x-coordinate
 */
boundaryType Universe::getMinXBoundaryType() {

  double min_x = std::numeric_limits<double>::infinity();
  boundaryType bc_x = REFLECTIVE;

  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {

    if (min_x > iter->second->getMinX()) {
      min_x = iter->second->getMinX();
      bc_x = iter->second->getMinXBoundaryType();
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

  double max_x = -std::numeric_limits<double>::infinity();
  boundaryType bc_x = REFLECTIVE;

  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    if (max_x < iter->second->getMaxX()) {
      max_x = iter->second->getMaxX();
      bc_x = iter->second->getMaxXBoundaryType();
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

  double min_y = std::numeric_limits<double>::infinity();
  boundaryType bc_y = REFLECTIVE;

  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {

    if (min_y > iter->second->getMinY()) {
      min_y = iter->second->getMinY();
      bc_y = iter->second->getMinYBoundaryType();
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

  double max_y = -std::numeric_limits<double>::infinity();
  boundaryType bc_y = REFLECTIVE;

  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    if (max_y < iter->second->getMaxY()) {
      max_y = iter->second->getMaxY();
      bc_y = iter->second->getMaxYBoundaryType();
    }
  }

  return bc_y;
}


/**
 * @brief Returns the boundary conditions (VACUUM or REFLECTIVE) at the minimum
 *         reachable z-coordinate in the Universe.
 * @return the boundary conditions at the minimum reachable z-coordinate
 */
boundaryType Universe::getMinZBoundaryType() {

  double min_z = std::numeric_limits<double>::infinity();
  boundaryType bc_z = REFLECTIVE;

  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {

    if (min_z > iter->second->getMinZ()) {
      min_z = iter->second->getMinZ();
      bc_z = iter->second->getMinZBoundaryType();
    }
  }

  return bc_z;
}


/**
 * @brief Returns the boundary conditions (VACUUM or REFLECTIVE) at the maximum
 *         reachable z-coordinate in the Universe.
 * @return the boundary conditions at the maximum reachable z-coordinate
 */
boundaryType Universe::getMaxZBoundaryType() {

  double max_z = -std::numeric_limits<double>::infinity();
  boundaryType bc_z = REFLECTIVE;

  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    if (max_z < iter->second->getMaxZ()) {
      max_z = iter->second->getMaxZ();
      bc_z = iter->second->getMaxZBoundaryType();
    }
  }

  return bc_z;
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
 * @brief Returns a CellFill in this Universe.
 * @param cell_id the integer the cell_id
 * @return Returns the CellFill pointer.
 */
CellFill* Universe::getCellFill(int cell_id) {

  CellFill* cell = NULL;
  if (_cells.find(cell_id) == _cells.end())
    log_printf(ERROR, "Unable to return Cell with ID = %d from Universe with "
               "ID = %d since it does not contain this Cell", cell_id, _id);

  cell = static_cast<CellFill*>(_cells.at(cell_id));
  if (cell->getType() != FILL)
    log_printf(WARNING, "Retrieving Cell %d from Universe %d, but it "
               "is not a FILL type Cell", cell->getId(), _id);
  return cell;
}


/**
 * @brief Returns a CellBasic in this Universe.
 * @param cell_id the integer the cell_id
 * @return Returns the CellFill pointer.
 */
CellBasic* Universe::getCellBasic(int cell_id) {

  CellBasic* cell = NULL;
  if (_cells.find(cell_id) == _cells.end())
    log_printf(ERROR, "Unable to return Cell with ID = %d from Universe with "
               "ID = %d since the it does not contain this Cell", cell_id, _id);

  cell = static_cast<CellBasic*>(_cells.at(cell_id));
  if (cell->getType() != MATERIAL)
    log_printf(WARNING, "Retrieving Cell %d from Universe %d, but it "
               "is not a MATERIAL type Cell", cell->getId(), _id);

  return cell;
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
      material = static_cast<CellBasic*>(cell)->getMaterial();
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
    cell = (*iter).second;
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
    log_printf(INFO, "Added Cell with ID = %d to Universe with ID = %d",
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

  Cell* return_cell = NULL;
  Cell* cell;
  std::map<int, Cell*>::iterator iter;

  /* Sets the LocalCoord type to UNIV at this level */
  coords->setType(UNIV);

  /* Loop over all Cells in this Universe */
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    cell = iter->second;

    if (cell->cellContainsCoords(coords)) {

      /* Set the Cell on this level */
      coords->setCell(cell);

      /* MATERIAL type Cell - lowest level, terminate search for Cell */
      if (cell->getType() == MATERIAL) {
        coords->setCell(cell);
        return_cell = cell;
        return return_cell;
      }

      /* FILL type Cell - Cell contains a Universe at a lower level
       * Update coords to next level and continue search */
      else if (cell->getType() == FILL) {

        LocalCoords* next_coords;

        if (coords->getNext() == NULL)
          next_coords = new LocalCoords(coords->getX(), coords->getY());
        else
          next_coords = coords->getNext();

        CellFill* fill = static_cast<CellFill*>(cell);
        Universe* univ = fill->getFill();
        next_coords->setUniverse(univ);
        coords->setCell(cell);

        coords->setNext(next_coords);
        next_coords->setPrev(coords);
        if (univ->getType() == SIMPLE)
          return univ->findCell(next_coords);
        else
          return static_cast<Lattice*>(univ)->findCell(next_coords);
      }
    }
  }

  return return_cell;
}


/**
 * @brief Finds the distance to the nearest surface.
 * @details Loops over all the cells within the universe and computes
 *          the distance to each one following the direction of the track.
 *          Returns distance to nearest next cell's nearest surface.
 * @param point a pointer to a starting point
 * @param angle the azimuthal angle of the track
 * @return the distance to the nearest surface
 */
double Universe::minSurfaceDist(Point* point, double angle) {

  Point min_intersection;
  std::map<int, Cell*>::iterator iter;
  double dist;
  double min_dist = INFINITY;

  /* Loop over all Cells in this Universe */
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    dist = iter->second->minSurfaceDist(point, angle, &min_intersection);
    min_dist = std::min(dist, min_dist);
  }

  return min_dist;
}


/**
 * @brief Subdivides all of the Cells within this Universe into rings
 *        and angular sectors.
 */
void Universe::subdivideCells() {

  log_printf(DEBUG, "Subdividing Cells for Universe %d", _id);

  std::map<int, Cell*>::iterator iter1;

  while (iter1 != _cells.end()) {

    for (iter1 = _cells.begin(); iter1 != _cells.end(); ++iter1) {

      if (((*iter1).second)->getType() == MATERIAL) {
        CellBasic* cell = static_cast<CellBasic*>((*iter1).second);

        if (cell->getNumRings() > 0 || cell->getNumSectors() > 0) {
          std::vector<CellBasic*> newcells = cell->subdivideCell();

          log_printf(DEBUG, "Cell %d in Universe %d has %d subcells",
                     cell->getId(), _id, newcells.size());

          std::vector<CellBasic*>::iterator iter2;
          for (iter2=newcells.begin(); iter2!=newcells.end(); ++iter2)
            addCell((*iter2));

          _cells.erase(iter1);
          break;
        }
      }
    }
  }
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
    if ((*iter1).second->getType() == MATERIAL) {

      /* Clone the Cell */
      CellBasic* parent = static_cast<CellBasic*>((*iter1).second);
      CellBasic* cell_clone = parent->clone();

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
  _offset.setCoords(0.0, 0.0);

  /* Default width and number of Lattice cells along each dimension */
  _num_y = 0;
  _num_x = 0;
  _width_x = 0;
  _width_y = 0;
}


/**
 * @brief Destructor clears memory for all of Universes pointers.
 */
Lattice::~Lattice() {

  for (int i=0; i < _num_y; i++)
    _universes.at(i).clear();

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
 */
void Lattice::setOffset(double x, double y) {
  _offset.setX(x);
  _offset.setY(y);
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
 * @brief Returns the minimum reachable x-coordinate in the Lattice.
 * @return the minimum reachable x-coordinate
 */
double Lattice::getMinX() {
  return _offset.getX() - (_num_x + _width_x / 2.);
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
double Lattice::getMaxY(){
  return _offset.getY() + (_num_y * _width_y / 2.);
}


/**
 * @brief Returns the minimum reachable z-coordinate in the Lattice.
 * @return the minimum reachable z-coordinate
 */
double Lattice::getMinZ() {
  return -std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns the maximum reachable z-coordinate in the Lattice.
 * @return the maximum reachable z-coordinate
 */
double Lattice::getMaxZ() {
  return std::numeric_limits<double>::infinity();
}


/**
 * @brief Returns a pointer to the Universe within a specific Lattice cell.
 * @param lat_x the x index to the Lattice cell
 * @param lat_y the y index to the Lattice cell
 * @return pointer to a Universe filling the Lattice cell
 */
Universe* Lattice::getUniverse(int lat_x, int lat_y) const {

  /* Checks that lattice indices are within the bounds of the lattice */
  if (lat_x > _num_x || lat_y > _num_y)
    log_printf(ERROR, "Cannot retrieve Universe from Lattice ID = %d: Index"
               "out of bounds: Tried to access Cell x = %d, y = %d but bounds"
               "are x = %d, y = %d", _id, lat_x, lat_y, _num_x, _num_y);

  return _universes.at(lat_y).at(lat_x).second;
}


/**
 * @brief Return a 2D vector of the Universes in the Lattice.
 * @return 2D vector of Universes
 */
std::vector<std::vector<std::pair<int, Universe*>>>
  Lattice::getUniverses() const {

  return _universes;
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

  for (int i = _num_y-1; i > -1;  i--) {
    for (int j = 0; j < _num_x; j++) {
      universe = _universes.at(i).at(j).second;
      unique_universes[universe->getId()] = universe;
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

  for (iter = unique_universes.begin(); iter != unique_universes.end(); ++iter){
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

  for (iter = unique_universes.begin(); iter != unique_universes.end(); ++iter){
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
 * @brief Set the width of each Lattice cell.
 * @param width_x the width along the x-axis in centimeters
 * @param width_y the width along the y-axis in centimeters
 */
void Lattice::setWidth(double width_x, double width_y) {

  if (width_x <= 0 || width_y <= 0)
    log_printf(ERROR, "Unable to set the width of Lattice ID = %d "
               "for x = %f and y = %f since they are not positive values",
               _id, width_x, width_y);

  _width_x = width_x;
  _width_y = width_y;
}


/**
 * @brief Sets the array of Universe pointers filling each Lattice cell.
 * @details This is a helper method for SWIG to allow users to assign Universes
 *          to a Lattice using a 2D Python list (list of lists). An example
 *          how this method can be called from Python is as follows:
 *
 * @code
 *          u1 = Universe(name='Universe 1')
 *          u2 = Universe(name='Universe 2')
 *          u3 = Universe(name='Universe 3')
 *          lattice.setLatticeCells([[u1, u2, u1, u2],
 *                                   [u2, u3, u2, u3],
 *                                   [u1, u2, u1, u2],
 *                                   [u2, u3, u2, u3]])
 * @endcode
 *
 * @param num_y the number of Lattice cells along y
 * @param num_x the number of Lattice cells along x
 * @param universes the array of Universes for each Lattice cell
 */
void Lattice::setUniverses(int num_y, int num_x, Universe** universes) {

  /* Clear any Universes in the Lattice (from a previous run) */
  for (int i=0; i < _num_y; i++)
    _universes.at(i).clear();

  _universes.clear();

  /* Set the Lattice dimensions */
  setNumX(num_x);
  setNumY(num_y);

  Universe* universe;

  /* The Lattice cells are assumed input in row major order starting from the
   * upper left corner. This double loop reorders the Lattice cells from the
   * to start from the lower left corner */
  for (int j = 0; j < _num_y; j++) {

    _universes.push_back(std::vector< std::pair<int, Universe*> >());

    for (int i = 0; i < _num_x; i++){
      universe = universes[(_num_y-1-j)*_num_x + i];
      _universes.at(j).push_back(std::pair<int, Universe*>
                                 (universe->getId(), universe));
    }
  }
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

  double x = point->getX();
  double y = point->getY();

  /* If the Point is outside the x bounds */
  if (x > bound_x_max || x < bound_x_min)
    return false;

  /* If the Point is outside the y boounds */
  else if (y > bound_y_max || y < bound_y_min)
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

  /* Compute the x and y indices for the Lattice cell this coord is in */
  int lat_x = getLatX(coords->getPoint());
  int lat_y = getLatY(coords->getPoint());

  /* If the indices are outside the bound of the Lattice */
  if (lat_x < 0 || lat_x >= _num_x ||
      lat_y < 0 || lat_y >= _num_y) {
    return NULL;
  }

  /* Compute local position of Point in the next level Universe */
  double nextX = coords->getX()
      - (-_width_x*_num_x/2.0 + _offset.getX() + (lat_x + 0.5) * _width_x)
      + getOffset()->getX();
  double nextY = coords->getY()
      - (-_width_y*_num_y/2.0 + _offset.getY() + (lat_y + 0.5) * _width_y)
      + getOffset()->getY();

  /* Create a new LocalCoords object for the next level Universe */
  LocalCoords* next_coords;

  if (coords->getNext() == NULL)
    next_coords = new LocalCoords(nextX, nextY);
  else
    next_coords = coords->getNext();

  Universe* univ = getUniverse(lat_x, lat_y);
  next_coords->setUniverse(univ);

  /* Set Lattice indices */
  coords->setLattice(this);
  coords->setLatticeX(lat_x);
  coords->setLatticeY(lat_y);

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
 *          Returns distance to nearest Lattice cell boundary.
 * @param point a pointer to a starting point
 * @param angle the azimuthal angle of the track
 * @return the distance to the nearest Lattice cell boundary
 */
double Lattice::minSurfaceDist(Point* point, double angle) {

  double next_x, next_y;
  double dist_x, dist_y;
  double dist_row, dist_col;

  /* Compute the x and y indices for the Lattice cell this point is in */
  int lat_x = getLatX(point);
  int lat_y = getLatY(point);

  /* find distance to next x plane crossing */
  if (angle < M_PI / 2.0)
    next_x = (lat_x + 1) * _width_x - _width_x*_num_x/2.0 + _offset.getX();
  else
    next_x = lat_x * _width_x - _width_x*_num_x/2.0 + _offset.getX();

  /* get distance to the nearest cell in the current row */
  next_y = point->getY() + tan(angle) * (next_x - point->getX());
  dist_x = fabs(next_x - point->getX());
  dist_y = fabs(next_y - point->getY());
  dist_row = pow(pow(dist_x, 2) + pow(dist_y, 2), 0.5);

  /* find distance to next y plane crossing */
  next_y = (lat_y + 1) * _width_y - _width_y*_num_y/2.0 + _offset.getY();
  next_x = point->getX() + (next_y - point->getY()) / tan(angle);
  dist_x = fabs(next_x - point->getX());
  dist_y = fabs(next_y - point->getY());
  dist_col = pow(pow(dist_x, 2) + pow(dist_y, 2), 0.5);

  /* return shortest distance to next lattice cell */
  return std::min(dist_row, dist_col);
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
 * @brief Converts a Lattice's attributes to a character array representation.
 * @return character array of this Lattice's attributes
 */
std::string Lattice::toString() {

  std::stringstream string;

  string << "Lattice ID = " << _id
         << ", name = " << _name
         << ", # cells along x = " << _num_x
         << ", # cells along y = " << _num_y
         << ", x width = " << _width_x
         << ", y width = " << _width_y;

  string << "\n\t\tUniverse IDs within this Lattice: ";

  for (int i = _num_y-1; i > -1;  i--) {
    for (int j = 0; j < _num_x; j++)
      string << _universes.at(i).at(j).first << ", ";
    string << "\n\t\t";
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
 * @details Lattice cells are numbered starting with 0 in the lower left
 *          corner. Lattice cell IDs in all rows then increase monotonically 
 *          from left to right. For example, the indices for a 4 x 4 lattice:
 *                  12  13  14  15
 *                  8    9  10  11
 *                  4    5   6   7
 *                  0    1   2   3
 * @param point a pointer to a point being evaluated.
 * @return the Lattice cell index.
 */
int Lattice::getLatticeCell(Point* point){
  return (getLatY(point)*_num_x + getLatX(point));
}


/**
 * @brief Finds the Lattice cell surface that a point lies on.
 *        If the point is not on a surface, -1 is returned.
 * @details The surface indices for a lattice cell are 0 (left),
 *         1, (bottom), 2 (right), 3 (top), 4 (bottom-left corner),
 *         5 (bottom-right corner), 6 (top-right corner), and 
 *         7 (top-left corner). The index returned takes into account
 *         the cell index and returns 8*cell_index + surface_index.
 * @param cell the cell index that the point is in.
 * @param point a pointer to a point being evaluated.
 * @return the Lattice surface index.
 */
int Lattice::getLatticeSurface(int cell, Point* point) {

  /* Get coordinates of point and cell boundaries */
  double x = point->getX();
  double y = point->getY();
  int lat_x = cell % _num_x;
  int lat_y = cell / _num_x;
  double left = lat_x*_width_x - _width_x*_num_x/2.0 + _offset.getX();
  double right = (lat_x + 1)*_width_x - _width_x*_num_x/2.0 + _offset.getX();
  double bottom = lat_y*_width_y - _width_y*_num_y/2.0 + _offset.getY();
  double top = (lat_y+1)*_width_y - _width_y*_num_y/2.0 + _offset.getY();
  int surface = -1;

  /* Check if point is on left boundary */ 
  if (fabs(x - left) <= ON_SURFACE_THRESH){
    /* Check if point is on bottom boundary */ 
    if (fabs(y - bottom) <= ON_SURFACE_THRESH)
      surface = cell*8 + 4;
    /* Check if point is on top boundary */ 
    else if (fabs(y - top) <= ON_SURFACE_THRESH)
      surface = cell*8 + 7;
    else
      surface = cell*8;
  }
  /* Check if point is on right boundary */ 
  else if (fabs(x - right) <= ON_SURFACE_THRESH){
    /* Check if point is on bottom boundary */ 
    if (fabs(y - bottom) <= ON_SURFACE_THRESH)
      surface = cell*8 + 5;
    /* Check if point is on top boundary */ 
    else if (fabs(y - top) <= ON_SURFACE_THRESH)
      surface = cell*8 + 6;
    else
      surface = cell*8 + 2;
  }
  /* Check if point is on bottom boundary */ 
  else if (fabs(y - bottom) <= ON_SURFACE_THRESH)
    surface = cell*8 + 1;
  /* Check if point is on top boundary */ 
  else if (fabs(y - top) <= ON_SURFACE_THRESH)
    surface = cell*8 + 3;

  return surface;
}
