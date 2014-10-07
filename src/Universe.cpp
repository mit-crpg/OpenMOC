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
 * @brief Constructor assigns a unique and user-specified ID for the Universe.
 * @param id the user-specified optional Universe ID
 * @param name the user-specified optional Universe ID
 */
Universe::Universe(const int id, const char* name) {
  _uid = _n;
  _id = id;
  _n++;
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
 * @brief Return a pointer to the origin for this Cell (in global coordinates).
 * @return the origin of the Cell
 */
Point* Universe::getOrigin() {
  return &_origin;
}


/**
 * @brief Aggregates a list (vector) of the IDs of all Materials within
 *        the MATERIAL type Cells filling this Universe.
 * @details Note that this method only searches the first level of Cells
 *          below this Universe within the nested Universe coordinate system.
 * @return a vector of Material IDs
 */
std::vector<int> Universe::getMaterialIds() {

  std::vector<int> material_ids;

  /* Otherwise, we don't yet know whether or not this Universe contains a Cell
   * with a Material with a non-zero (fissionable) fission cross-section */
  std::map<int, Cell*>::iterator iter;

  /* Loop over all Cells in this Universe */
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {

    Cell* cell = iter->second;

    /* If this Cell is filled by a Material, add the Material's ID
     * to the list */
    if (cell->getType() == MATERIAL)
      material_ids.push_back(static_cast<CellBasic*>(cell)->getMaterial());
  }

  return material_ids;
}


/**
 * @brief Aggregates a list (vector) of the IDs of all Universes within
 *        the FILL type Cells filling this Universe.
 * @details Note that this method only searches the first level of Cells
 *          below this Universe within the nested Universe coordinate system.
 * @return a vector of Universe IDs
 */
std::vector<int> Universe::getNestedUniverseIds() {

  std::vector<int> universe_ids;

  /* Loop over all Cells in this Universe */
  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {

    Cell* cell = iter->second;

    /* If this Cell is filled by a Material, add the Material's ID
     * to the list */
    if (cell->getType() == FILL)
      universe_ids.push_back(static_cast<CellFill*>(cell)->getUniverseFillId());
    }

  return universe_ids;
}


/**
 * @brief Aggregates a list (vector) of the IDs of all Cells in the Universe.
 * @details Note that this method only searches the first level of Cells
 *          below this Universe within the nested Universe coordinate system.
 *          This is a helper method for SWIG to allow a user to pass in a single
 *          NumPy array as an argument. In particular, the method permits
 *          "Universe cloning" in Python, a technique which permits the easy
 *          creation of a large heterogeneous Contructive Solid Geometry with
 *          different Materials in each region.
 *
 *          An example of how this might be used in Python is given as follows:
 *
 * @code
 *          num_cells = universe.getNumCells()
 *          cell_ids = numpy.zeros(num_cells)
 *          universe.getCellIds(cell_ids)
 * @endcode
 *
 * @param cell_ids an array to populate with Cell IDs
 * @param num_cells the number of Cells in the Universe
 */
void Universe::getCellIds(int* cell_ids, int num_cells) {

  int counter = 0;

  /* Loop over all Cells in this Universe */
  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    Cell* cell = iter->second;
    cell_ids[counter] = cell->getId();
    counter++;
  }

  return;
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
 * @brief Returns the local ID for the FSR representing a Cell in this Universe.
 * @details This method is used when constructing an ID for an FSR and should
 *          only be called by the Geometry. Users should NOT call this method
 *          directly to retreive an FSR ID.
 * @param cell_id the ID of the cell of interest
 */
int Universe::getFSR(int cell_id) {

  if (_cells.find(cell_id) == _cells.end())
    log_printf(ERROR, "Tried to find FSR ID for Cell with ID = %d in "
               " Universe with ID = %d but no Cell exists", cell_id, _id);

  return _region_map.at(cell_id);
}


/**
 * @brief Return the container of Cell IDs and Cell pointers in this Universe.
 * @return std::map of Cell IDs
 */
std::map<int, Cell*> Universe::getCells() const {
  return _cells;
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
 * @brief Sets the name of the Universe
 * @param name the Universe name string
 */
void Universe::setName(const char* name) {
  int length = strlen(name);

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
 * @brief Set the origin in global coordinates for this Universe.
 * @param origin a pointer to the origin
 */
void Universe::setOrigin(Point* origin) {
  _origin.setX(origin->getX());
  _origin.setY(origin->getY());
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
 * @brief Finds the Cell for which a LocalCoords object resides.
 * @details Finds the Cell that a LocalCoords object is located inside by
 *          checking each of this Universe's Cells. Returns NULL if the
 *          LocalCoords is not in any of the Cells.
 * @param coords a pointer to the LocalCoords of interest
 * @param universes a container of all of the Universes passed in by Geometry
 * @return a pointer the Cell where the LocalCoords is located
 */
Cell* Universe::findCell(LocalCoords* coords,
                         std::map<int, Universe*> universes) {

  Cell* return_cell = NULL;
  std::map<int, Cell*>::iterator iter;

  /* Sets the LocalCoord type to UNIV at this level */
  coords->setType(UNIV);

  /* Loop over all Cells in this Universe */
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    Cell* cell = iter->second;

    if (cell->cellContainsCoords(coords)) {

      /* Set the Cell on this level */
      coords->setCell(cell->getId());

      /* MATERIAL type Cell - lowest level, terminate search for Cell */
      if (cell->getType() == MATERIAL) {
        coords->setCell(cell->getId());
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

        CellFill* cell_fill = static_cast<CellFill*>(cell);
        int universe_id = cell_fill->getUniverseFillId();
        next_coords->setUniverse(universe_id);
        Universe* univ = universes.at(universe_id);
        coords->setCell(cell->getId());

        coords->setNext(next_coords);
        next_coords->setPrev(coords);
        if (univ->getType() == SIMPLE)
          return univ->findCell(next_coords, universes);
        else
          return static_cast<Lattice*>(univ)->findCell(next_coords, universes);
      }
    }
  }

  return return_cell;
}


/**
 * @brief Convert the member attributes of this Universe to a character array.
 * @return a character array representing the Universe's attributes
 */
std::string Universe::toString() {

  std::stringstream string;
  std::map<int, Cell*>::iterator iter;

  string << "Universe ID = " << _id;
  string << ", name = " << _name << ", type = ";

  if (_type == SIMPLE)
    string << "SIMPLE";
  else
    string << "LATTICE";

  string << ", num cells = " << _cells.size() << ", cell IDs = ";

  for (iter = _cells.begin(); iter != _cells.end(); ++iter)
    string << iter->first << ", ";

  return string.str();
}


/**
 * @brief Compute the FSR offset maps for this Universe and return the number of
 *        FSRs inside the Universe.
 * @details The FSR map is simply a std::map of the local FSR IDs for each
 *          of the Cells making up this Universe. This is used when
 *          constructing a global FSR ID for a Lattice which may contain
 *          several repeating versions of this Universe.
 * @return the number of FSRs in the Universe
 */
int Universe::computeFSRMaps() {

  /* Initialize a counter */
  std::map<int, Cell*>::iterator iter;
  int count = 0;

  /* Iterate over Cells in the Universe to set the map and update count */
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {
    _region_map.insert(std::pair<int, int>(iter->first, count));
    count += iter->second->getNumFSRs();
  }

  return count;
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

      if ((*iter1).second->getType() == MATERIAL) {
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
 * @brief Prints a string representation of the Universe's attributes to
 *        the console.
 */
void Universe::printString() {
  log_printf(RESULT, toString().c_str());
}


/**
 * @brief Clones this Universe and all of the Cells within it and returns it
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
      cell_clone->setUniverse(clone->getId());
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
 * @param width_x the width of the Lattice cells along x
 * @param width_y the width of the Lattice cells along y
 * @param name the user-specified optional Lattice (Universe) name
 */
Lattice::Lattice(const int id, double width_x, double width_y,
                 const char* name):
  Universe(id, name) {

  _width_x = width_x;
  _width_y = width_y;
  _type = LATTICE;

  /* Default number of Lattice cells along each dimension */
  _num_y = 0;
  _num_x = 0;
}


/**
 * @brief Destructor clears memory for all of Universes pointers.
 */
Lattice::~Lattice() {

  for (int i=0; i < _num_x; i++)
    _universes.at(i).clear();

  _universes.clear();
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
 * @brief Return the origin of the Lattice.
 * @return a pointer to the origin of the Lattice
 */
Point* Lattice::getOrigin() {
  return &_origin;
}


/**
 * @brief Return a 2D vector of the Universes in the Lattice.
 * @return 2D vector of Universes
 */
std::vector< std::vector< std::pair<int, Universe*> > >
Lattice::getUniverses() const {
  return _universes;
}


/**
 * @brief Returns a pointer to the Universe within a specific Lattice cell.
 * @param lattice_x the x index to the Lattice cell
 * @param lattice_y the y index to the Lattice cell
 * @return pointer to a Universe filling the Lattice cell
 */
Universe* Lattice::getUniverse(int lattice_x, int lattice_y) const {

  /* Checks that lattice indices are within the bounds of the lattice */
  if (lattice_x > _num_x || lattice_y > _num_y)
    log_printf(ERROR, "Cannot retrieve Universe from Lattice ID = %d: Index"
               "out of bounds: Tried to access Cell x = %d, y = %d but bounds"
               "are x = %d, y = %d", _id, lattice_x, lattice_y, _num_x, _num_y);

  return _universes.at(lattice_y).at(lattice_x).second;
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
 * @brief Return the id of a flat source region base index (smallest FSR
 *        region id within a specific Lattice cell)
 * @details This method is used when constructing an ID for an FSR and should
 *          only be called by the Geometry. Users should NOT call this method
 *          directly to retreive an FSR ID.
 * @param lat_x the x index of the Lattice cell
 * @param lat_y the y index of the Lattice cell
 * @return the index of the base FSR
 */
int Lattice::getFSR(int lat_x, int lat_y) {

  /* Check if Lattice indices are out of bounds */
  if (lat_x > _num_x || lat_y > _num_y)
    log_printf(ERROR, "Tried to access FSR map of Lattice ID = %d, but indices"
               "lat_x = %d and lat_y = %d were out of bounds",
               _id, lat_x, lat_y);

  return _region_map[lat_y][lat_x].second;
}


/**
 * @brief Aggregates a list (vector) of the IDs of all Universes within
 *        the FILL type Cells filling this Universe.
 * @details Note that this method only searches the first level of Cells
 *          below this Universe within the nested Universe coordinate system.
 * @return a vector of Universe IDs
 */
std::vector<int> Lattice::getNestedUniverseIds() {

  std::vector<int> universe_ids;

  /* Loop over all Lattice cells */
  for (int i = 0; i < _num_y; i++) {
    for (int j = 0; j< _num_x; j++)
      universe_ids.push_back(_universes[i][j].second->getId());
  }

  return universe_ids;
}


/**
 * @brief Sets the pointer to a Universe filling one of this Lattice's
 *        Lattice cells.
 * @param universe pointer to the Universe
 */
void Lattice::setUniversePointer(Universe* universe) {

  /* Check that _surfaces contains this Surface ID and delete the ID
   *  otherwise throw an error */
  int universe_id = universe->getId();
  bool universe_not_found = true;

  /* Loop over all Universes in Lattice */
  for (int i = 0; i < _num_y; i++) {
    for (int j = 0; j< _num_x; j++) {

      /* Checks if the Universe ID matches what the Lattice's
       * Universes container expects */
      if (_universes.at(i).at(j).first == universe_id) {
        _universes[i][j].second = universe;
        universe_not_found = false;
      }
    }
  }

  if (universe_not_found)
    log_printf(WARNING, "Tried to set the Universe pointer for "
              "Lattice id = %d for Universe ID = %d but the Lattice "
               "does not contain the Universe", _id, universe_id);
  else
    log_printf(INFO, "Set the Universe pointer for Lattice "
               "ID = %d for Universe ID = %d", _id, universe_id);

  return;
}


/**
 * @brief Sets the arrary of Universe IDs filling each Lattice cell.
 * @details This is a helper method for SWIG to allow users to assign Universe
 *          IDs to a Lattice using a 2D Python list (list of lists). An example
 *          how this method can be called from Python is as follows:
 *
 * @code
 *          lattice.setLatticeCells([[1, 2, 1, 2],
 *                                   [2, 3, 2, 3],
 *                                   [1, 2, 1, 2],
 *                                   [2, 3, 2, 3]])
 * @endcode
 *
 * @param universes the array of Universe IDs for each Lattice cell
 * @param num_x the number of Lattice cells along x
 * @param num_y the number of Lattice cells along y
 */
void Lattice::setLatticeCells(int num_x, int num_y, int* universes) {

  /* Clear any Universes in the Lattice (from a previous run) */
  for (int i=0; i < _num_x; i++)
    _universes.at(i).clear();

  _universes.clear();

  /* Set the Lattice dimensions */
  _num_x = num_x;
  _num_y = num_y;

  /* Set origin to lower left corner with (x=0, y=0) in center of Lattice */
  _origin.setX(-_width_x*_num_x/2.0);
  _origin.setY(-_width_y*_num_y/2.0);

  Universe* empty_universe_pointer = NULL;

  /* The Lattice cells are assumed input in row major order starting from the
   * upper left corner. This double loop reorders the Lattice cells from the
   * to start from the lower left corner */
  for (int i = 0; i < _num_y; i++) {

    _universes.push_back(std::vector< std::pair<int, Universe*> >());

    for (int j = 0; j< _num_x; j++){
      _universes.at(i).push_back(std::pair<int, Universe*>
                  (universes[(_num_y-1-i)*_num_x + j], empty_universe_pointer));
    }
  }

  /* Intialize _region_map */
  for (int i = 0; i < _num_y; i++) {

    _region_map.push_back(std::vector< std::pair<int, int> >());

    for (int j = 0; j < _num_x; j++) {
      _region_map.at(i).push_back(std::pair<int, int>());
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
  double bound_x = _num_x/2.0 * _width_x;
  double bound_y = _num_y/2.0 * _width_y;

  double x = point->getX();
  double y = point->getY();

  /* If the Point is outside the x bounds */
  if (x > bound_x || x < -1*bound_x)
    return false;

  /* If the Point is outside the y boounds */
  else if (y > bound_y || y < -1*bound_y)
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
 * @param universes a std::map of all Universes passed in from the geometry
 * @return a pointer to the Cell this LocalCoord is in or NULL
 */
Cell* Lattice::findCell(LocalCoords* coords,
                        std::map<int, Universe*> universes) {

  /* Set the LocalCoord to be a LAT type at this level */
  coords->setType(LAT);

  /* Compute the x and y indices for the Lattice cell this coord is in */
  int lat_x = (int)floor((coords->getX() - _origin.getX()) / _width_x);
  int lat_y = (int)floor((coords->getY() - _origin.getY()) / _width_y);

  /* Check if the LocalCoord is on the Lattice boundaries and if so adjust
   * x or y Lattice cell indices */
  if (fabs(fabs(coords->getX()) - _num_x*_width_x*0.5) <
      ON_LATTICE_CELL_THRESH) {

    if (coords->getX() > 0)
      lat_x = _num_x - 1;
    else
      lat_x = 0;
  }
  if (fabs(fabs(coords->getY()) - _num_y*_width_y*0.5) <
      ON_LATTICE_CELL_THRESH) {
    if (coords->getY() > 0)
      lat_y = _num_y - 1;
    else
      lat_y = 0;
  }

  /* If the indices are outside the bound of the Lattice */
  if (lat_x < 0 || lat_x >= _num_x ||
      lat_y < 0 || lat_y >= _num_y) {
    return NULL;
  }

  /* Compute local position of Point in the next level Universe */
  double nextX = coords->getX() - (_origin.getX() + (lat_x + 0.5) * _width_x);
  double nextY = coords->getY() - (_origin.getY() + (lat_y + 0.5) * _width_y);

  /* Create a new LocalCoords object for the next level Universe */
  LocalCoords* next_coords;

  if (coords->getNext() == NULL)
    next_coords = new LocalCoords(nextX, nextY);
  else
    next_coords = coords->getNext();

  int universe_id = getUniverse(lat_x, lat_y)->getId();
  Universe* univ = universes.at(universe_id);
  next_coords->setUniverse(universe_id);

  /* Set Lattice indices */
  coords->setLattice(_id);
  coords->setLatticeX(lat_x);
  coords->setLatticeY(lat_y);

  coords->setNext(next_coords);
  next_coords->setPrev(coords);

  /* Search the next lowest level Universe for the Cell */
  return univ->findCell(next_coords, universes);
}


/**
 * @brief Finds the next Cell for a LocalCoords object along a trajectory
 *        defined by some angle (in radians from 0 to pi).
 * @details The method will update the LocalCoords passed in as an argument
 *          to be the one at the boundary of the next Cell crossed along the
 *          given trajectory. It will do this by recursively building a linked
 *          list of LocalCoords from the LocalCoords passed in as an argument
 *          down to the lowest level Cell found. In the process it will set
 *          the local coordinates for each LocalCoords in the linked list for
 *          the Lattice or Universe that it is in. If the LocalCoords is
 *          outside the bounds of the Lattice or on the boundaries this method
 *          will return NULL; otherwise it will return a pointer to the Cell
 *          that the LocalCoords will reach next along its trajectory.
 * @param coords pointer to a LocalCoords object
 * @param angle the angle of the trajectory
 * @param universes a map of all of the Universes passed in by the Geometry
 * @return a pointer to a Cell if found, NULL if no cell found
 */
Cell* Lattice::findNextLatticeCell(LocalCoords* coords, double angle,
                                   std::map<int, Universe*> universes) {

  /* Tests the upper, lower, left and right Lattice Cells adjacent to
   * the LocalCoord and uses the one with the shortest distance from
   * the current location of the LocalCoord */

  /* Initial distance is infinity */
  double distance = std::numeric_limits<double>::infinity();

  /* Current minimum distance */
  double d;

  /* Properties of the current LocalCoords */
  double x0 = coords->getX();
  double y0 = coords->getY();
  int lattice_x = coords->getLatticeX();
  int lattice_y = coords->getLatticeY();

  /* Slope of trajectory */
  double m = sin(angle) / cos(angle);

  /* Properties of the new location for LocalCoords */
  /* Current point of minimum distance */
  double x_curr, y_curr;

  /* x-coordinate on new Lattice cell */
  double x_new = x0;

  /* y-coordinate on new Lattice cell */
  double y_new = x0;

  /* New x Lattice cell index */
  int new_lattice_x;

  /* New y Lattice cell index */
  int new_lattice_y;

  /* Test Point for computing distance */
  Point test;

  /* Check lower Lattice cell */
  if (lattice_y >= 0 && angle >= M_PI) {
    y_curr = (lattice_y - _num_y/2.0) * _width_y;
    x_curr = x0 + (y_curr - y0) / m;
    test.setCoords(x_curr, y_curr);

    /* Check if the test Point is within the bounds of the Lattice */
    if (withinBounds(&test)) {
      d = test.distanceToPoint(coords->getPoint());

      /* Check if distance to test Point is current minimum */
      if (d < distance) {
        distance = d;
        x_new = x_curr;
        y_new = y_curr;
      }
    }
  }

  /* Upper Lattice cell */
  if (lattice_y <= _num_y-1 && angle <= M_PI) {
    y_curr = (lattice_y - _num_y/2.0 + 1) * _width_y;
    x_curr = x0 + (y_curr - y0) / m;
    test.setCoords(x_curr, y_curr);

    /* Check if the test Point is within the bounds of the Lattice */
    if (withinBounds(&test)) {
      d = test.distanceToPoint(coords->getPoint());

      /* Check if distance to test Point is current minimum */
      if (d < distance) {
        distance = d;
        x_new = x_curr;
        y_new = y_curr;
      }
    }
  }

  /* Left Lattice cell */
  if (lattice_x >= 0 && (angle >= M_PI/2 && angle <= 3*M_PI/2)) {
    x_curr = (lattice_x - _num_x/2.0) * _width_x;
    y_curr = y0 + m * (x_curr - x0);
    test.setCoords(x_curr, y_curr);

    /* Check if the test Point is within the bounds of the Lattice */
    if (withinBounds(&test)) {
      d = test.distanceToPoint(coords->getPoint());

      /* Check if distance to test Point is current minimum */
      if (d < distance) {
        distance = d;
        x_new = x_curr;
        y_new = y_curr;
      }
    }
  }

  /* Right Lattice cell */
  if (lattice_x <= _num_x-1 && (angle <= M_PI/2 || angle >= 3*M_PI/2)) {
    x_curr = (lattice_x - _num_x/2.0 + 1) * _width_x;
    y_curr = y0 + m * (x_curr - x0);
    test.setCoords(x_curr, y_curr);

    /* Check if the test Point is within the bounds of the Lattice */
    if (withinBounds(&test)) {
      d = test.distanceToPoint(coords->getPoint());

      /* Check if distance to test Point is current minimum */
      if (d < distance) {
        distance = d;
        x_new = x_curr;
        y_new = y_curr;
      }
    }
  }

  /* If no Point was found on the Lattice cell, then the LocalCoords was
   * already on the boundary of the Lattice */
  if (distance == INFINITY)
    return NULL;

  /* Otherwise a Point was found inside a new Lattice cell */
  else {
    /* Update the Localcoords location to the Point on the new Lattice cell
     * plus a small bit to ensure that its coordinates are inside cell */
    double delta_x = (x_new - coords->getX()) + cos(angle) * TINY_MOVE;
    double delta_y = (y_new - coords->getY()) + sin(angle) * TINY_MOVE;
    coords->adjustCoords(delta_x, delta_y);

    /* Compute the x and y indices for the new Lattice cell */
    new_lattice_x = (int)floor((coords->getX() - _origin.getX())/_width_x);
    new_lattice_y = (int)floor((coords->getY() - _origin.getY())/_width_y);

    /* Check if the LocalCoord is on the lattice boundaries and if so adjust
     * x or y Lattice cell indices i */
    if (fabs(fabs(coords->getX()) - _num_x*_width_x*0.5) <
        ON_LATTICE_CELL_THRESH) {

      if (coords->getX() > 0)
        new_lattice_x = _num_x - 1;
      else
        new_lattice_x = 0;
    }
    if (fabs(fabs(coords->getY()) - _num_y*_width_y*0.5) <
        ON_LATTICE_CELL_THRESH) {

      if (coords->getY() > 0)
        new_lattice_y = _num_y - 1;
      else
        new_lattice_y = 0;
    }

    /* Check if new Lattice cell indices are within the bounds, if not,
     * new LocalCoords is now on the boundary of the Lattice */
    if (new_lattice_x >= _num_x || new_lattice_x < 0)
      return NULL;
    else if (new_lattice_y >= _num_y || new_lattice_y < 0)
      return NULL;

    /* New LocalCoords is still within the interior of the Lattice */
    else {
      /* Update the LocalCoords Lattice cell indices */
      coords->setLatticeX(new_lattice_x);
      coords->setLatticeY(new_lattice_y);

      /* Move to next lowest level Universe */
      coords->prune();
      Universe* univ = _universes.at(new_lattice_y).at(new_lattice_x).second;
      LocalCoords* next_coords;

      /* Compute local position of Point in next level Universe */
      double nextX = coords->getX() - (_origin.getX()
                     + (new_lattice_x + 0.5) * _width_x);
      double nextY = coords->getY() - (_origin.getY()
                     + (new_lattice_y + 0.5) * _width_y);

      /* Set the coordinates at the next level LocalCoord */
      next_coords = new LocalCoords(nextX, nextY);
      next_coords->setPrev(coords);
      coords->setNext(next_coords);

      next_coords->setUniverse(univ->getId());

      /* Search lower level Universe */
      return findCell(coords, universes);
    }
  }
}


/**
 * @brief Computes the flat source region base indices for each of the
 *        Lattice cells within this Lattice (i.e., the minimum ID for the flat
 *        source regions within each Lattice cell). Returns the number of
 *        FSRs in the Lattice.
 * @return the number of FSRs
 */
int Lattice::computeFSRMaps() {

  /* Initialize a counter */
  int count = 0;

  /* loop over Universes in the Lattice to set the map and update count */
  for (int i = 0; i < _num_y; i++) {
    for (int j = 0; j < _num_x; j++) {
      Universe *u = _universes.at(i).at(j).second;
      _region_map[i][j].second = count;
      count += u->computeFSRMaps();
    }
  }

  return count;
}


/**
 * @brief Converts a Lattice's attributes to a character array representation.
 * @return character array of this Lattice's attributes
 */
std::string Lattice::toString() {

  std::stringstream string;

  string << "Lattice ID = " << _id << ", name = " << _name
         << ", num cells along x = " << _num_x
         << ", num cells along y = " << _num_y << ", x width = "
         << _width_x << ", y width = " << _width_y;

  string << "\n\t\tUniverse IDs within this Lattice:\n\t\t";

  for (int i = _num_y-1; i > -1;  i--) {
    for (int j = 0; j < _num_x; j++)
      string << _universes.at(i).at(j).first << "  ";
    string << "\n\t\t";
  }

  return string.str().c_str();
}


/**
 * @brief Prints a string representation of all of the Lattice's attributes to
 *        the console.
 */
void Lattice::printString() {
  log_printf(RESULT, toString().c_str());
}
