#include "Geometry.h"


/**
 * @brief Constructor initializes an empty Geometry.
 */
Geometry::Geometry() {

  /* Initializing the corners of the bounding box encapsulating
   * the Geometry to be infinite  */
  _x_min = std::numeric_limits<double>::max();
  _y_min = std::numeric_limits<double>::max();
  _x_max = -std::numeric_limits<double>::max();
  _y_max = -std::numeric_limits<double>::max();

  _max_seg_length = 0;
  _min_seg_length = std::numeric_limits<double>::infinity();

  /* Default boundary conditions are reflective */
  _top_bc    = REFLECTIVE;
  _bottom_bc = REFLECTIVE;
  _left_bc   = REFLECTIVE;
  _right_bc  = REFLECTIVE;

  _num_FSRs = 0;
  _num_groups = 0;

  /* Initialize CMFD object to NULL */
  _cmfd = NULL;

  /* initialize _num_FSRs lock */
  _num_FSRs_lock = new omp_lock_t;
  omp_init_lock(_num_FSRs_lock);
}


/**
 * @brief Destructor clears all references to Materials, Surfaces, Cells,
 *        Universes and Lattices.
 */
Geometry::~Geometry() {

  _materials.clear();
  _surfaces.clear();
  _cells.clear();
  _universes.clear();
  _lattices.clear();

  /* Free FSR  maps if they were initialized */
  if (_num_FSRs != 0) {
    _FSR_keys_map.clear();
    _FSRs_to_keys.clear();
    _FSRs_to_material_IDs.clear();
  }
}


/**
 * @brief initialize the CellFill objects with the Universes filling them.
 */
void Geometry::initializeCellFillPointers() {

  std::map<int, Cell*>::iterator iter;
  CellFill* cell;
  Universe* univ;

  /* Iterates over all Cells in the Geometry */
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {

    /* Checks if CellFill references this Universe and sets its pointer */
    if (iter->second->getType() == FILL) {
      cell = static_cast<CellFill*>(iter->second);
      univ = _universes.at(cell->getUniverseFillId());
      cell->setUniverseFillPointer(univ);
    }
  }

  return;
}


/**
 * @brief Returns the total height (y extent) of the Geometry in cm.
 * @return the total height of the Geometry (cm)
 */
double Geometry::getHeight() {
  return (_y_max - _y_min);
}


/**
 * @brief Returns the total width (x extent) of the Geometry in cm.
 * @return the total width of the Geometry (cm)
 */
double Geometry::getWidth() {
  return (_x_max - _x_min);
}


/**
 * @brief Return the minimum x-coordinate contained by the Geometry.
 * @return the minimum x-coordinate (cm)
 */
double Geometry::getXMin() {
  return _x_min;
}


/**
 * @brief Return the maximum x-coordinate contained by the Geometry.
 * @return the maximum x-coordinate (cm)
 */
double Geometry::getXMax() {
  return _x_max;
}


/**
 * @brief Return the minimum y-coordinate contained by the Geometry.
 * @return the minimum y-coordinate (cm)
 */
double Geometry::getYMin() {
  return _y_min;
}


/**
 * @brief Return the maximum y-coordinate contained by the Geometry.
 * @return the maximum y-coordinate (cm)
 */
double Geometry::getYMax() {
  return _y_max;
}


/**
 * @brief Returns the boundary condition for the top Surface of the Geometry.
 * @details The boundary conditions are vacuum (false) and reflective (false).
 * @return the boundary conditions for the top of the Geometry
 */
boundaryType Geometry::getBCTop() {
  return _top_bc;
}


/**
 * @brief Returns the boundary condition for the bottom Surface of the Geometry.
 * @details The boundary conditions are vacuum (false) and reflective (false).
 * @return the boundary conditions for the bottom of the Geometry
 */
boundaryType Geometry::getBCBottom() {
  return _bottom_bc;
}


/**
 * @brief Returns the boundary condition for the left Surface of the Geometry.
 * @details The boundary conditions are vacuum (false) and reflective (false).
 * @return the boundary conditions for the left surface of the Geometry
 */
boundaryType Geometry::getBCLeft() {
  return _left_bc;
}


/**
 * @brief Returns the boundary condition for the right Surface of the Geometry.
 * @details The boundary conditions are vacuum (false) and reflective (false).
 * @return the boundary conditions for the right surface of the Geometry
 */
boundaryType Geometry::getBCRight() {
  return _right_bc;
}


/**
 * @brief Returns the number of flat source regions (FSRs) in the Geometry.
 * @return number of FSRs
 */
int Geometry::getNumFSRs() {
  return _num_FSRs;
}


/**
 * @brief Sets the number of flat source regions (FSRs) in the Geometry.
 * @param num_fsrs number of FSRs
 */
void Geometry::setNumFSRs(int num_fsrs) {
  _num_FSRs = num_fsrs;
}


/**
 * @brief Returns the number of energy groups for each Material's nuclear data.
 * @return the number of energy groups
 */
int Geometry::getNumEnergyGroups() {
  if (getNumMaterials() == 0)
    log_printf(ERROR, "Unable to return the number of energy groups from "
               "the geometry since it does not contain any materials");

  return _num_groups;
}


/**
 * @brief Returns the number of Materials in the Geometry.
 * @return the number of Materials
 */
int Geometry::getNumMaterials() {
  return _materials.size();
}


/**
 * @brief Returns the number of Cells in the Geometry.
 * @return the number of Cells
 */
int Geometry::getNumCells() {
  return _cells.size();
}


/**
 * @brief Return the max Track segment length computed during segmentation (cm)
 * @return max Track segment length (cm)
 */
double Geometry::getMaxSegmentLength() {
  return _max_seg_length;
}


/**
 * @brief Return the min Track segment length computed during segmentation (cm)
 * @return min Track segment length (cm)
 */
double Geometry::getMinSegmentLength() {
  return _min_seg_length;
}


/**
 * @brief Return a std::map container of Material IDs (keys) with Materials
 *        pointers (values).
 * @return a std::map of Materials indexed by Material ID in the geometry
 */
std::map<int, Material*> Geometry::getMaterials() {
  return _materials;
}


/**
 * @brief Return a pointer to a Material object in the Geometry.
 * @param id the user-specified Material ID
 * @return a pointer to the Material object
 */
Material* Geometry::getMaterial(int id) {

  Material* material = NULL;

  try {
    material = _materials.at(id);
  }
  catch (std::exception & e) {
    log_printf(ERROR, "Attempted to retrieve Material with ID = %d which"
               " does not exist. Backtrace:\n%s", id, e.what());
  }

  return material;
}


/**
 * @brief Return a pointer to a Surface from the Geometry.
 * @param id the user-specified Surface ID
 * @return a pointer to the Surface object
 */
Surface* Geometry::getSurface(int id) {

  Surface* surface = NULL;

  try {
    surface = _surfaces.at(id);
  }
  catch (std::exception & e) {
    log_printf(ERROR, "Attempted to retrieve Surface with ID = %d which "
               "has not been declared. Backtrace:\n%s", id, e.what());
  }

  return surface;
}

/**
 * @brief Return a pointer to a Cell object in the Geometry.
 * @param id the user-specified Cell's ID
 * @return a pointer to the Cell object
 */
Cell* Geometry::getCell(int id) {

  Cell* cell = NULL;

  try {
    cell = static_cast<Cell*>(_cells.at(id));
  }
  catch (std::exception & e) {
    log_printf(ERROR, "Attempted to retrieve Cell with ID = %d which has "
               "not been declared. Backtrace:\n%s", id, e.what());
  }

  return cell;
}


/**
 * @brief Return a pointer to a Cell filled by a Material (CellBasic)
 *        from the Geometry.
 * @param id the user-specified Cell's ID
 * @return a pointer to the Cell object
 */
CellBasic* Geometry::getCellBasic(int id) {

  CellBasic* cell = NULL;

  try {
    cell = static_cast<CellBasic*>(_cells.at(id));
    if (cell->getType() != MATERIAL)
      log_printf(WARNING, "Retrieving Cell %d from the Geometry, but it "
                 "is not a MATERIAL type Cell", cell->getId());
  }
  catch (std::exception & e) {
    log_printf(ERROR, "Attempted to retrieve Cell with ID = %d which has "
               "not been declared. Backtrace:\n%s", id, e.what());
  }

  return cell;
}


/**
 * @brief Return a pointer to a Cell filled by a Universe (CellFill)
 *        from the Geometry.
 * @param id the user-specified Cell's ID
 * @return a pointer to the Cell object
 */
CellFill* Geometry::getCellFill(int id) {

  CellFill* cell = NULL;

  try {
    cell = static_cast<CellFill*>(_cells.at(id));
    if (cell->getType() != FILL)
      log_printf(WARNING, "Retrieving Cell %d from the Geometry, but it "
                 "is not a FILL type Cell", cell->getId());
  }
  catch (std::exception & e) {
    log_printf(ERROR, "Attempted to retrieve Cell with ID = %d which has "
               "not been declared. Backtrace:\n%s", id, e.what());
  }

  return cell;
}


/**
 * @brief Return a pointer to a Universe from the Geometry.
 * @param id the user-specified Universe ID
 * @return a pointer to the Universe object
 */
Universe* Geometry::getUniverse(int id) {

  Universe* universe = NULL;

  try {
    universe = _universes.at(id);
  }
  catch (std::exception & e) {
    log_printf(ERROR, "Attempted to retrieve Universe with ID = %d which "
               "has not been declared. Backtrace:\n%s", id, e.what());
  }

  return universe;
}


/**
 * @brief Return a pointer to a Lattice from the Geometry.
 * @param id the user-specified Lattice (Universe) ID
 * @return a pointer to the Lattice object
 */
Lattice* Geometry::getLattice(int id) {

  Lattice* lattice = NULL;
  try {
    lattice = _lattices.at(id);
  }
  catch (std::exception & e) {
    log_printf(ERROR, "Attempted to retrieve Lattice with ID = %d which "
               "has not been declared. Backtrace:\n%s", id, e.what());
  }

  return lattice;
}


/**
 * @brief Returns a pointer to the CMFD object.
 * @return A pointer to the CMFD object
 */
Cmfd* Geometry::getCmfd(){
  return _cmfd;
}


/**
 * @brief Sets the pointer to a CMFD object used for acceleration.
 * @param A pointer to the CMFD object
 */
void Geometry::setCmfd(Cmfd* cmfd){
  _cmfd = cmfd;
}


/**
 * @brief Add a Material to the Geometry.
 * @param material a pointer to a Material object
 */
void Geometry::addMaterial(Material* material) {

  /* Checks the number of energy groups */
  if (material->getNumEnergyGroups() == 0)
    log_printf(ERROR, "Unable to add Material %d since it does not "
               "contain any nuclear data", material->getId());

  if (_num_groups == 0)
    _num_groups = material->getNumEnergyGroups();

  else if (_num_groups != material->getNumEnergyGroups())
    log_printf(ERROR, "Unable to add Material %d with %d energy groups to the"
               " Geometry which contains Material(s) with %d energy groups",
               material->getId(), material->getNumEnergyGroups(), _num_groups);

  try {
    /* Check that the sum of the Material's absorption and scattering
     * cross-sections equals its total cross-section */
    material->checkSigmaT();
    _materials.insert(std::pair<int,Material*>(material->getId(), material));
    log_printf(INFO, "Added Material with ID = %d to Geometry",
               material->getId());
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to add Material with ID = %d to Geometry. "
               "Backtrace:\n%s", material->getId(), e.what());
  }
}


/**
 * @brief Add a Surface to the Geometry.
 * @details This method will update the boundary conditions and bounding box
 *          encapsulating the Geometry if the Surface has VACUUUM or REFLECTIVE
 *          boundary conditions.
 * @param surface a pointer to the Surface object
 */
void Geometry::addSurface(Surface* surface) {

  try {
    _surfaces.insert(std::pair<int, Surface*>(surface->getId(), surface));
    log_printf(INFO, "Added Surface with ID = %d to Geometry",surface->getId());
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to add Surface with ID = %d to Geometry. "
               "Backtrace:\n%s", surface->getId(), e.what());
  }


  /* Use new surface to update the boundaries of the Geometry */
  switch (surface->getBoundaryType()) {
  case REFLECTIVE:
    if (surface->getXMin() < _x_min &&
        surface->getXMin() !=-std::numeric_limits<double>::infinity()) {
      _x_min = surface->getXMin();
      _left_bc = REFLECTIVE;
     }

    if (surface->getXMax() > _x_max &&
        surface->getXMax() != std::numeric_limits<double>::infinity()) {
      _x_max = surface->getXMax();
      _right_bc = REFLECTIVE;
    }

    if (surface->getYMin() < _y_min &&
        surface->getYMin() !=-std::numeric_limits<double>::infinity()) {
      _y_min = surface->getYMin();
      _bottom_bc = REFLECTIVE;
    }
    if (surface->getYMax() > _y_max &&
        surface->getYMax() != std::numeric_limits<double>::infinity()) {
     _y_max = surface->getYMax();
     _top_bc = REFLECTIVE;
    }
    break;

  case VACUUM:
    if (surface->getXMin() < _x_min &&
        surface->getXMin() !=-std::numeric_limits<double>::infinity()) {
      _x_min = surface->getXMin();
      _left_bc = VACUUM;
    }

    if (surface->getXMax() > _x_max &&
        surface->getXMax() != std::numeric_limits<double>::infinity()) {
      _x_max = surface->getXMax();
      _right_bc = VACUUM;
    }

    if (surface->getYMin() < _y_min &&
        surface->getYMin() !=-std::numeric_limits<double>::infinity()) {
      _y_min = surface->getYMin();
      _bottom_bc = VACUUM;
    }

    if (surface->getYMax() > _y_max &&
        surface->getYMax() != std::numeric_limits<double>::infinity()) {
      _y_max = surface->getYMax();
      _top_bc = VACUUM;
    }
    break;

  case BOUNDARY_NONE:
    break;
  }

  return;
}


/**
 * @brief Add a Cell to the Geometry.
 * @details This method checks if the Universe the Cell is in already exists.
 *          If not, the Universe is created and added to the Geometry.
 * @param cell a pointer to the Cell object
 */
void Geometry::addCell(Cell* cell) {

  /* Prints error msg if the Cell is filled with a non-existent Material */
  if (cell->getType() == MATERIAL &&
           _materials.find(static_cast<CellBasic*>(cell)->getMaterial()) ==
           _materials.end()) {

    log_printf(ERROR, "Attempted to add Cell with Material with ID = %d,"
               " but the Geometry does not contain this Material",
    static_cast<CellBasic*>(cell)->getMaterial());
  }

  /* Add the Cell's surfaces to the Geometry's surfaces map */
  std::map<int, surface_halfspace> cells_surfaces = cell->getSurfaces();
  std::map<int, surface_halfspace>::iterator iter;
  for (iter = cells_surfaces.begin(); iter != cells_surfaces.end(); ++iter)
    addSurface(iter->second._surface);

  /* Insert the Cell into the Geometry's Cell container */
  try {
    _cells.insert(std::pair<int, Cell*>(cell->getId(), cell));
    log_printf(INFO, "Added Cell with ID = %d to Geometry", cell->getId());
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to add Cell with ID = %d to Geometry. "
               "Backtrace:\n%s", cell->getId(), e.what());
  }

  /* Checks if the Universe the Cell in exists; if not, creates Universe */
  if (_universes.find(cell->getUniverseId()) == _universes.end()) {
    try {
      Universe* univ = new Universe(cell->getUniverseId());
      addUniverse(univ);
      log_printf(INFO, "Created Universe with ID = %d", cell->getUniverseId());
    }
    catch (std::exception &e) {
      log_printf(ERROR, "Unable to create a new Universe with ID = %d "
                 "and add it to the Geometry. Backtrace:\n%s",
                 cell->getUniverseId(), e.what());
    }
  }

  /* Adds the cell to the appropriate universe */
  _universes.at(cell->getUniverseId())->addCell(cell);

  return;
}


/**
 * @brief Add a Universe to the Geometry.
 * @details This method sets the pointer to the Universe for each CellFill
 *          containing this Universe.
 * @param universe a pointer to the Universe object
 */
void Geometry::addUniverse(Universe* universe) {

  /* Add the Universe */
  try {
    _universes[universe->getId()] = universe;
    log_printf(INFO, "Added Universe with ID = %d to Geometry. ",
               universe->getId());
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to add Universe with ID = %d to Geometry. "
               "Backtrace:\n%s", universe->getId(), e.what());
  }

  /* Checks if any CellFill references this Universe and sets its pointer */
  std::map<int, Cell*>::iterator iter;
  for (iter = _cells.begin(); iter != _cells.end(); ++iter) {

    if (iter->second->getType() == FILL) {
      CellFill* cell = static_cast<CellFill*>(iter->second);

      if (cell->getUniverseFillId() == universe->getId()) {
        cell->setUniverseFillPointer(universe);
      }
    }
  }

  return;
}


/**
 * @brief Add a Lattice to the Geometry.
 * @details Adds the Lattice to both the Lattice and Universe containers. This
 *          method sets the pointers to the Universes contained by each
 *          Lattice cell.
 * @param lattice a pointer to the Lattice object
 */
void Geometry::addLattice(Lattice* lattice) {

  /* Sets the Universe pointers for the Lattice and checks if the Lattice
   * contains a Universe which does not exist */
  for (int i = 0; i < lattice->getNumY(); i++) {
    for (int j = 0; j < lattice->getNumX(); j++) {
      int universe_id = lattice->getUniverses().at(i).at(j).first;

      /* If the Universe does not exist */
      if (_universes.find(universe_id) == _universes.end())
        log_printf(ERROR, "Attempted to create Lattice containing Universe "
                   "with ID = %d, but the Geometry does not contain this "
                   "Universe", lattice->getUniverses().at(i).at(j).first);

      /* Set the Universe pointer */
      else
        lattice->setUniversePointer(_universes.at(universe_id));
    }
  }

  /* Add the Lattice to the Geometry's :attices container */
  try {
    _lattices.insert(std::pair<int, Lattice*>(lattice->getId(), lattice));
    log_printf(INFO, "Added Lattice with ID = %d to Geometry",
    lattice->getId());
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to add Lattice with ID = %d to Geometry. "
               ". Backtrace:\n%s", lattice->getId(), e.what());
  }

  /* Add the Lattice to the Universes container as well */
  addUniverse(lattice);
}


/**
 * @brief Removes a Material from the Geometry.
 * @details Note: this method does not remove the Cells filled by this Material
 *          from the Geometry.
 * @param id the Material ID
 */
void Geometry::removeMaterial(int id) {

  /* Checks if the Geometry contains this Material */
  if (_materials.find(id) == _materials.end())
    log_printf(WARNING, "Cannot remove a Material with ID = %d from the "
               "geometry since it doesn't contain that Material", id);

  try {
    std::map<int, Material*>::iterator iter;
    iter = _materials.find(id);
    _materials.erase(iter);
    log_printf(INFO, "Removed Material with ID = %d from geometry", id);
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to remove a Material with ID = %d from "
               "Geometry. Backtrace:\n%s", id, e.what());
  }
}


/**
 * @brief Removes a Cell from the Geometry.
 * @details Note: this method does not remove the Universe, Surface(s), or
 *          Material contained within the Cell.
 * @param id the Cell ID
 */
void Geometry::removeCell(int id) {

  /* Checks if the Geometry contains this Cell */
  if (_cells.find(id) == _cells.end())
    log_printf(WARNING, "Cannot remove a Cell with ID = %d from the "
               "Geometry since it doesn't contain that Cell", id);

  try {
    std::map<int, Cell*>::iterator iter;
    iter = _cells.find(id);
    _cells.erase(iter);
    log_printf(INFO, "Removed Cell with ID = %d from Geometry", id);
  }
  catch (std::exception &e) {
    log_printf(ERROR, "Unable to remove a Cell with ID = %d "
               "from Geometry. Backtrace:\n%s", id, e.what());
  }
}



/**
 * @brief Removes a Universe from the Geometry.
 * @details Note: this method does not remove the Cells contained within the
 *          Universe, or Cells filled by this Universe.
 * @param id the Universe ID
 */
void Geometry::removeUniverse(int id) {

  /* Checks if the Geometry contains this Universe */
  if (_universes.find(id) == _universes.end())
    log_printf(WARNING, "Cannot remove a Universe with ID = %d from the "
               "Geometry since it doesn't contain that Universe", id);

  /* Remove the Universe from the Geometry */
  else {
    try {
      std::map<int, Universe*>::iterator iter;
      iter = _universes.find(id);
      _universes.erase(iter);
      log_printf(INFO, "Removed Universe with ID = %d from Geometry", id);
    }
    catch (std::exception &e) {
      log_printf(ERROR, "Unable to remove Universe with ID = %d from "
                 "Geometry. Backtrace:\n%s", id, e.what());
    }
  }

  return;
}


/**
 * @brief Removes a Lattice from the Geometry.
 * @details Note: this method does not remove the Universes filling each
 *          Lattice cell in the Lattice.
 * @param id the Lattice ID
 */
void Geometry::removeLattice(int id) {

  /* Checks if the Geometry contains this Universe */
  if (_lattices.find(id) == _lattices.end())
    log_printf(WARNING, "Cannot remove a Lattice with ID = %d from the "
               "Geometry since it doesn't contain that Lattice", id);

  /* Remove the Universe from the Geometry */
  else {
    try {
      std::map<int, Lattice*>::iterator iter;
      iter = _lattices.find(id);
      _lattices.erase(iter);
      log_printf(INFO, "Removed Lattice with ID = %d from Geometry", id);
    }
    catch (std::exception &e) {
      log_printf(ERROR, "Unable to remove Lattice with ID = %d from "
                 "Geometry. Backtrace:\n%s", id, e.what());
    }
  }
}



/**
 * @brief Find the Cell that this LocalCoords object is in at the lowest level
 *        of the nested Universe hierarchy.
 * @details This method assumes that the LocalCoords has been initialized
 *          with coordinates and a Universe ID. The method will recursively
 *          find the Cell on the lowest level of the nested Universe hierarchy
 *          by building a linked list of LocalCoords from the LocalCoord
 *          passed in as an argument down to the lowest level Cell found. In
 *          the process it will set the coordinates at each level of the
 *          hierarchy for each LocalCoord in the linked list for the Lattice
 *          or Universe that it is in. If the LocalCoords is outside the bounds
 *          of the Geometry or on the boundaries this method will return NULL;
 *          otherwise it will return a pointer to the Cell that is found by the
 *          recursive Geometry::findCell(...) method.
 * @param coords pointer to a LocalCoords object
 * @return returns a pointer to a Cell if found, NULL if no Cell found
 */
CellBasic* Geometry::findCellContainingCoords(LocalCoords* coords) {

  int universe_id = coords->getUniverse();
  Universe* univ = _universes.at(universe_id);

  if (universe_id == 0){
    if (!withinBounds(coords))
      return NULL;
  }

  if (univ->getType() == SIMPLE)
    return static_cast<CellBasic*>(univ->findCell(coords, _universes));
  else
    return static_cast<CellBasic*>(static_cast<Lattice*>(univ)->findCell
                                   (coords, _universes));
}


/**
 * @brief Find the first Cell of a Track segment with a starting Point that is
 *        represented by the LocalCoords method parameter.
 * @details This method assumes that the LocalCoords has been initialized
 *          with coordinates and a Universe ID. This method will move the
 *          initial starting point by a small amount along the direction of
 *          the Track in order to ensure that the track starts inside of a
 *          distinct FSR rather than on the boundary between two of them.
 *          The method will recursively find the LocalCoords by building a
 *          linked list of LocalCoords from the LocalCoords passed in as an
 *          argument down to the Cell found in the lowest level of the nested
 *          Universe hierarchy. In the process, the method will set the
 *          coordinates at each level in the nested Universe hierarchy for
 *          each LocalCoord in the linked list for the Lattice or Universe
 *          that it is in.
 * @param coords pointer to a LocalCoords object
 * @param angle the angle for a trajectory projected from the LocalCoords
 * @return returns a pointer to a cell if found, NULL if no cell found
*/
CellBasic* Geometry::findFirstCell(LocalCoords* coords, double angle) {
  double delta_x = cos(angle) * TINY_MOVE;
  double delta_y = sin(angle) * TINY_MOVE;
  coords->adjustCoords(delta_x, delta_y);
  return findCellContainingCoords(coords);
}


/**
 * @brief Find the Material for a flat source region ID.
 * @details  This method finds the fsr_id within the 
 *           _FSR_to_material_IDs map and returns the corresponding
 *           pointer to the Material object.
 * @param fsr_id a FSR id
 * @return a pointer to the Material that this FSR is in
 */
Material* Geometry::findFSRMaterial(int fsr_id) {
  return getMaterial(_FSRs_to_material_IDs.at(fsr_id));
}


/**
 * @brief Finds the next Cell for a LocalCoords object along a trajectory
 *        defined by some angle (in radians from 0 to Pi).
 * @details The method will update the LocalCoords passed in as an argument
 *          to be the one at the boundary of the next Cell crossed along the
 *          given trajectory. It will do this by finding the minimum distance
 *          to the surfaces at all levels of the coords hierarchy.
 *          If the LocalCoords is outside the bounds of the Geometry or on 
 *          the boundaries this method will return NULL; otherwise it will 
 *          return a pointer to the Cell that the LocalCoords will reach 
 *          next along its trajectory.
 * @param coords pointer to a LocalCoords object
 * @param angle the angle of the trajectory
 * @return a pointer to a Cell if found, NULL if no Cell found
 */
CellBasic* Geometry::findNextCell(LocalCoords* coords, double angle) {

  Cell* cell = NULL;
  double dist;
  double min_dist = std::numeric_limits<double>::infinity();
  Point surf_intersection;

  /* Find the current Cell */
  cell = findCellContainingCoords(coords);

  /* Get lowest level coords */
  coords = coords->getLowestLevel();

  /* If the current coords is not in any Cell, return NULL */
  if (cell == NULL)
    return NULL;

  /* If the current coords is inside a Cell, look for next Cell */
  else {

    /* Ascend universes until at the highest level.
    * At each universe/lattice level get distance to next 
    * universe or lattice cell. Recheck min_dist. */
    while (coords != NULL) {

      /* If we reach a LocalCoord in a Lattice, find the distance to the
      * nearest lattice cell boundary */
      if (coords->getType() == LAT) {
        Lattice* lattice = _lattices.at(coords->getLattice());
        dist = lattice->minSurfaceDist(coords->getPoint(), angle);
      }
      /* If we reach a LocalCoord in a Universe, find the distance to the
      * nearest cell surface */
      else{
        cell = _cells.at(coords->getCell());
        dist = cell->minSurfaceDist(coords->getPoint(), angle, &surf_intersection);
      }

      /* Recheck min distance */
      min_dist = std::min(dist, min_dist);

      /* Ascend one level */
      if (coords->getUniverse() == 0)
        break;
      else{
        coords = coords->getPrev();
        coords->prune();
      }
    }

    /* Check for distance to nearest CMFD mesh cell boundary */
    if (_cmfd != NULL){
      Lattice* lattice = _cmfd->getLattice();
      dist = lattice->minSurfaceDist(coords->getPoint(), angle);
      min_dist = std::min(dist, min_dist);
    }

    /* Move point and get next cell */
    double delta_x = cos(angle) * (min_dist + TINY_MOVE);
    double delta_y = sin(angle) * (min_dist + TINY_MOVE);
    coords->adjustCoords(delta_x, delta_y);
    return findCellContainingCoords(coords);
  }
}


/**
 * @brief Find and return the ID of the flat source region that a given
 *        LocalCoords object resides within.
 * @param coords a LocalCoords object pointer
 * @return the FSR ID for a given LocalCoords object
 */
int Geometry::findFSRId(LocalCoords* coords) {

  int fsr_id = 0;
  LocalCoords* curr = coords;
  curr = coords->getLowestLevel();
  std::hash<std::string> key_hash_function;

  /* Generate unique FSR key */
  std::size_t fsr_key_hash = key_hash_function(getFSRKey(coords));

  /* If FSR has not been encountered, update FSR maps and vectors */
  if (_FSR_keys_map.find(fsr_key_hash) == _FSR_keys_map.end()){
      
    /* Get the cell that contains coords */
    CellBasic* cell = findCellContainingCoords(curr);
    
    /* Get the lock */
    omp_set_lock(_num_FSRs_lock);

    /* Recheck to see if FSR has been added to maps after getting the lock */
    if (_FSR_keys_map.find(fsr_key_hash) != _FSR_keys_map.end())
      fsr_id = _FSR_keys_map.at(fsr_key_hash)._fsr_id;
    else{

        /* Add FSR information to FSR key map and FSR_to vectors */
      fsr_id = _num_FSRs;
      fsr_data* fsr = new fsr_data;
      fsr->_fsr_id = fsr_id;
      Point* point = new Point();
      point->setCoords(coords->getHighestLevel()->getX(), 
                       coords->getHighestLevel()->getY());
      fsr->_point = point;
      _FSR_keys_map[fsr_key_hash] = *fsr;
      _FSRs_to_keys.push_back(fsr_key_hash);
      _FSRs_to_material_IDs.push_back(cell->getMaterial());

      /* If CMFD acceleration is on, add FSR to CMFD cell */
      if (_cmfd != NULL){
        int cmfd_cell = _cmfd->findCmfdCell(coords->getHighestLevel());
        _cmfd->addFSRToCell(cmfd_cell, fsr_id);
      }

      /* Increment FSR counter */
      _num_FSRs++;
    }

    /* Release lock */
    omp_unset_lock(_num_FSRs_lock);

  }
  /* If FSR has already been encountered, get the fsr id from map */
  else
    fsr_id = _FSR_keys_map.at(fsr_key_hash)._fsr_id;
  
  return fsr_id;
}


/**
 * @brief Return the ID of the flat source region that a given
 *        LocalCoords object resides within.
 * @param coords a LocalCoords object pointer
 * @return the FSR ID for a given LocalCoords object
 */
int Geometry::getFSRId(LocalCoords* coords) {
    
  int fsr_id = 0;
  std::string fsr_key;
  std::hash<std::string> key_hash_function;
  
  try{
    fsr_key = getFSRKey(coords);
    fsr_id = _FSR_keys_map.at(key_hash_function(fsr_key))._fsr_id;
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not find FSR ID with key: %s. Try creating "
               "geometry with finer track laydown. "
               "Backtrace:%s", fsr_key.c_str(), e.what());
  }

  return fsr_id;  
}


/**
 * @brief Return the characteristic point for a given FSR ID
 * @param fsr_id the FSR ID
 * @return the FSR's characteristic point
 */
Point* Geometry::getFSRPoint(int fsr_id) {

  Point* point;
    
  try{
    point = _FSR_keys_map.at(_FSRs_to_keys.at(fsr_id))._point;
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not find characteristic poitn in FSR: %i. "
               "Backtrace:%s", fsr_id, e.what());
  }

  return point;  
}


/**
 * @brief Generate a string FSR "key" that identifies an FSR by its
 *        unique hierarchical lattice/universe/cell structure.
 * @detail Since not all FSRs will reside on the absolute lowest universe
 *         level and Cells might overlap other cells, it is important to
 *         have a method for uniquely identifying FSRs. This method
 *         createds a unique FSR key by constructing a structured string
 *         that describes the hierarchy of lattices/universes/cells.
 * @param coords a LocalCoords object pointer
 * @return the FSR key
 */
std::string Geometry::getFSRKey(LocalCoords* coords) {

  std::stringstream key;
  LocalCoords* curr = coords->getHighestLevel();
  std::ostringstream curr_level_key;

  /* If CMFD is on, get CMFD latice cell and write to key */
  if (_cmfd != NULL){
      curr_level_key << _cmfd->getLattice()->getLatX(curr->getPoint());
      key << "CMFD = (" << curr_level_key.str() << ", ";
      curr_level_key.str(std::string());
      curr_level_key << _cmfd->getLattice()->getLatY(curr->getPoint());
      key << curr_level_key.str() << ") : ";
  }

  /* Descend the linked list hierarchy until the lowest level has
   * been reached */
  while(curr != NULL){
      
    /* Clear string stream */
    curr_level_key.str(std::string());
      
    if (curr->getType() == LAT) {

      /* Write lattice ID and lattice cell to key */
      curr_level_key << curr->getLattice();
      key << "LAT = " << curr_level_key.str() << " (";
      curr_level_key.str(std::string());
      curr_level_key << curr->getLatticeX();
      key << curr_level_key.str() << ", ";
      curr_level_key.str(std::string());
      curr_level_key << curr->getLatticeY();
      key << curr_level_key.str() << ") : ";
    }
    else{

      /* write universe ID to key */
      curr_level_key << curr->getUniverse();
      key << "UNIV = " << curr_level_key.str() << " : ";
    }

    /* If lowest coords reached break; otherwise get next coords */
    if (curr->getNext() == NULL)
      break;
    else
      curr = curr->getNext();
  }

  /* clear string stream */
  curr_level_key.str(std::string());

  /* write cell id to key */
  curr_level_key << curr->getCell();
  key << "CELL = " << curr_level_key.str();  

  return key.str();
}



/**
 * @brief Subidivides all Cells in the Geometry into rings and angular sectors.
 * @details This method is called by the Geometry::initializeFlatSourceRegions()
 *          method but may also be called by the user in Python if needed:
 *
 * @code
 *          geometry.subdivideCells()
 * @endcode
 */
void Geometry::subdivideCells() {

  std::map<int, Universe*>::iterator iter;
  Universe* curr;

  std::map<int, Cell*>::iterator iter1;
  std::map<int, Cell*> cells;

  /* Loop over all Universe in the Geometry and instruct each to inform
   * their Cells to subdivide into rings and sectors as specified by
   * the user during Cell instantiation */
  for (iter = _universes.begin(); iter != _universes.end(); ++iter) {
    curr = (*iter).second;
    curr->subdivideCells();
    cells = curr->getCells();
    for (iter1 = cells.begin(); iter1 != cells.end(); ++iter1)
      addCell((*iter1).second);
  }
}


/**
 * @brief Compute the number of flat source regions in the Geometry and
 *        initialize arrays for FSR IDs and maps.
 * @details This method is intended to be called by the user before initiating
 *          source iteration. This method first subdivides all Cells by calling
 *          the Geometry::subdivideCells() method. Then it computes the total
 *          number of FSRs in the Geometry and initializes integer arrays of
 *          maps to Cells and Materials UIDs/IDs indexed by FSR IDs.
 */
void Geometry::initializeFlatSourceRegions() {

  /* Initialize pointers from CellFills to Universes */
  initializeCellFillPointers();

  /* Subdivide Cells into sectors and rings */
  subdivideCells();

  /* Assign UIDs to materials */
  std::map<int, Material*>::iterator iter;
  int uid = 0;
  for (iter = _materials.begin(); iter != _materials.end(); ++iter){
    iter->second->setUid(uid);
    uid++;
  }
  
  /* Initialize CMFD */
  if (_cmfd != NULL)
    initializeCmfd();
}


/**
 * @brief This method performs ray tracing to create Track segments within each
 *        flat source region in the Geometry.
 * @details This method starts at the beginning of a Track and finds successive
 *          intersection points with FSRs as the Track crosses through the
 *          Geometry and creates segment structs and adds them to the Track.
 * @param track a pointer to a track to segmentize
 */
void Geometry::segmentize(Track* track) {

  /* Track starting Point coordinates and azimuthal angle */
  double x0 = track->getStart()->getX();
  double y0 = track->getStart()->getY();
  double phi = track->getPhi();

  /* Length of each segment */
  FP_PRECISION segment_length;
  Material* segment_material;
  int fsr_id;
  FP_PRECISION* sigma_t;
  int min_num_segments;
  int num_segments;

  /* Use a LocalCoords for the start and end of each segment */
  LocalCoords segment_start(x0, y0);
  LocalCoords segment_end(x0, y0);
  segment_start.setUniverse(0);
  segment_end.setUniverse(0);

  /* Find the Cell containing the Track starting Point */
  Cell* curr = findFirstCell(&segment_end, phi);
  Cell* prev;

  /* If starting Point was outside the bounds of the Geometry */
  if (curr == NULL)
    log_printf(ERROR, "Could not find a Cell containing the start Point "
               "of this Track: %s", track->toString().c_str());

  /* While the end of the segment's LocalCoords is still within the Geometry,
   * move it to the next Cell, create a new segment, and add it to the
   * Geometry */
  while (curr != NULL) {

    segment_end.copyCoords(&segment_start);

    /* Find the next Cell along the Track's trajectory */
    prev = curr;
    curr = findNextCell(&segment_end, phi);

    /* Checks to make sure that new Segment does not have the same start
     * and end Points */
    if (segment_start.getX() == segment_end.getX() &&
      segment_start.getY() == segment_end.getY()) {

      log_printf(ERROR, "Created a Track segment with the same start and end "
                 "point: x = %f, y = %f", segment_start.getX(),
                  segment_start.getY());
    }

    /* Find the segment length between the segment's start and end points */
    segment_length = FP_PRECISION(segment_end.getPoint()
                      ->distanceToPoint(segment_start.getPoint()));
    segment_material = _materials.at(static_cast<CellBasic*>(prev)
                       ->getMaterial());
    sigma_t = segment_material->getSigmaT();

    /* Find the ID of the FSR that contains the segment */
    fsr_id = findFSRId(&segment_start);

    /* Compute the number of Track segments to cut this segment into to ensure
     * that it's length is small enough for the exponential table */
    min_num_segments = 1;
    for (int e=0; e < _num_groups; e++) {
      num_segments = ceil(segment_length * sigma_t[e] / 10.0);
      if (num_segments > min_num_segments)
      min_num_segments = num_segments;
    }

    /* "Cut up" Track segment into sub-segments such that the length of each
     * does not exceed the size of the exponential table in the Solver */
    for (int i=0; i < min_num_segments; i++) {

      /* Create a new Track segment */
      segment* new_segment = new segment;
      new_segment->_material = segment_material;
      new_segment->_length = segment_length / FP_PRECISION(min_num_segments);

      /* Update the max and min segment lengths */
      if (segment_length > _max_seg_length)
        _max_seg_length = segment_length;
      if (segment_length < _min_seg_length)
        _min_seg_length = segment_length;

      log_printf(DEBUG, "segment start x = %f, y = %f, segment end "
                 "x = %f, y = %f", segment_start.getX(), segment_start.getY(),
                 segment_end.getX(), segment_end.getY());

      new_segment->_region_id = fsr_id;

      /* Save indicies of CMFD Mesh surfaces that the Track segment crosses */
      if (_cmfd != NULL){
          
        /* Find cmfd cell that segment lies in */
        int cmfd_cell = _cmfd->findCmfdCell(&segment_start);

        /* Reverse nudge from surface to determine whether segment start or end
         * points lie on a cmfd surface. */
        double delta_x = cos(phi) * TINY_MOVE;
        double delta_y = sin(phi) * TINY_MOVE;
        segment_start.adjustCoords(-delta_x, -delta_y);
        segment_end.adjustCoords(-delta_x, -delta_y);

        if (i == min_num_segments-1)
          new_segment->_cmfd_surface_fwd = 
              _cmfd->findCmfdSurface(cmfd_cell, &segment_end);
        else
          new_segment->_cmfd_surface_fwd = -1;

        if (i == 0)
          new_segment->_cmfd_surface_bwd = 
              _cmfd->findCmfdSurface(cmfd_cell, &segment_start);
        else
          new_segment->_cmfd_surface_bwd = -1;

        /* Re-nudge segments from surface. */
        segment_start.adjustCoords(delta_x, delta_y);
        segment_end.adjustCoords(delta_x, delta_y);

      }

      /* Add the segment to the Track */
      track->addSegment(new_segment);

    }
  }

  log_printf(DEBUG, "Created %d segments for Track: %s",
             track->getNumSegments(), track->toString().c_str());

  /* Truncate the linked list for the LocalCoords */
  segment_start.prune();
  segment_end.prune();

  log_printf(DEBUG, "Track %d max. segment length: %f",
             track->getUid(), _max_seg_length);
  log_printf(DEBUG, "Track %d min. segment length: %f",
             track->getUid(), _min_seg_length);

  return;
}


/**
 * @brief Determines the fissionability of each Universe within this Geometry.
 * @details A Universe is determined fissionable if it contains a CellBasic
 *          filled by a Material with a non-zero fission cross-section. Note
 *          that this method recurses through all Universes at each level in
 *          the nested Universe hierarchy. Users should only call this method
 *          without a parameter (the default) from Python as follows to ensure
 *          that the recursion starts from the uppermost Universe level:
 *
 * @code
 *          geometry.computeFissionability()
 * @endcode
 *
 * @param univ the Universe of interest (default is NULL)
 */
void Geometry::computeFissionability(Universe* univ) {

  bool fissionable = false;
  std::vector<int> material_ids;
  std::vector<int> universe_ids;

  /* If no Universe was passed in as an argument, then this is the first
   * recursive call from a user via Python, so get the base Universe */
  if (univ == NULL)
    univ = _universes.at(0);

  /* If a Universe was passed in as an argument, then this is a recursive
   * call with a Universe at a lower level in the nested Universe hierarchy */
  if (univ->getType() == SIMPLE) {
    material_ids = univ->getMaterialIds();
    universe_ids = univ->getNestedUniverseIds();
  }

  else
    universe_ids = static_cast<Lattice*>(univ)->getNestedUniverseIds();

  /* Loop over the nested Universes first to ensure that fissionability
   * is set at each nested Universe level */
  for (int i=0; i < universe_ids.size(); i++) {
    int universe_id = universe_ids[i];
    Universe* universe = _universes.at(universe_id);

    /* Recursively check whether this nested Universe is fissionable */
    computeFissionability(universe);

    if (universe->isFissionable())
      fissionable = true;
  }

  /* Loop over the Materials in this Universe at this level */
  for (int i=0; i < material_ids.size(); i++) {
    int material_id = material_ids[i];
    Material* material = _materials.at(material_id);

    /* Check whether this Material is fissionable or not */
    if (material->isFissionable())
      fissionable = true;
  }

  /* Set this Universe's fissionability based on the nested Universes
   * and Materials within it */
  univ->setFissionability(fissionable);
}


/**
 * @brief Converts this Geometry's attributes to a character array.
 * @details This method calls the toString() method for all Materials,
 *          Surfaces, Cell, Universes and Lattices contained by the Geometry.
 * @return a character array of this Geometry's class attributes
 */
std::string Geometry::toString() {

  std::stringstream string;
  std::map<int, Material*>::iterator iter1;
  std::map<int, Surface*>::iterator iter2;
  std::map<int, Cell*>::iterator iter3;
  std::map<int, Universe*>::iterator iter4;
  std::map<int, Lattice*>::iterator iter5;

  string << "Geometry: width = " << getWidth() << ", height = "
         << getHeight() << ", Bounding Box: ((" << _x_min << ", "
         << _y_min << "), (" << _x_max << ", " << _y_max << ")";

  string << "\n\tMaterials:\n\t\t";
  for (iter1 = _materials.begin(); iter1 != _materials.end(); ++iter1)
    string << iter1->second->toString() << "\n\n\t\t";

  string << "\n\tSurfaces:\n\t\t";
  for (iter2 = _surfaces.begin(); iter2 != _surfaces.end(); ++iter2)
    string << iter2->second->toString() << "\n\t\t";

  string << "\n\tCells:\n\t\t";
  for (iter3 = _cells.begin(); iter3 != _cells.end(); ++iter3)
    string << iter3->second->toString() << "\n\t\t";

  string << "\n\tUniverses:\n\t\t";
  for (iter4 = _universes.begin(); iter4 != _universes.end(); ++iter4)
    string << iter4->second->toString() << "\n\t\t";

  string << "\n\tLattices:\n\t\t";
  for (iter5 = _lattices.begin(); iter5 != _lattices.end(); ++iter5)
    string << iter5->second->toString()  << "\n\t\t";

  std::string formatted_string = string.str();
  formatted_string.erase(formatted_string.end()-3);

  return formatted_string;
}


/**
 * @brief Prints a string representation of all of the Geometry's attributes to
 *        the console.
 * @details This method calls the printString() method for all Materials,
 *          Surfaces, Cell, Universes and Lattices contained by the Geometry.
 */
void Geometry::printString() {
  log_printf(RESULT, toString().c_str());
}



/**
 * @brief This is a method that initializes the CMFD Lattice and sets
 *          CMFD parameters.
 */
void Geometry::initializeCmfd(){

  /* Get information about geometry and CMFD mesh */
  int num_x = _cmfd->getNumX();
  int num_y = _cmfd->getNumY();
  double height = getHeight();
  double width = getWidth();
  double cell_width = width / num_x;
  double cell_height = height / num_y;
  
  /* Create CMFD lattice and set properties */
  Lattice* lattice = new Lattice(0, cell_width, cell_height);
  lattice->setNumX(num_x);
  lattice->setNumY(num_y);
  lattice->setOffset(_x_min + getWidth()/2.0, _y_min + getHeight()/2.0);
  _cmfd->setLattice(lattice);


  /* Set CMFD mesh boundary conditions */
  _cmfd->setBoundary(0,getBCLeft());
  _cmfd->setBoundary(1,getBCBottom());
  _cmfd->setBoundary(2,getBCRight());
  _cmfd->setBoundary(3,getBCTop());

  /* Set CMFD mesh dimensions and number of groups */
  _cmfd->setWidth(width);
  _cmfd->setHeight(height);
  _cmfd->setNumMOCGroups(_num_groups);

  /* If user did not set CMFD group structure, create CMFD group
  * structure that is the same as the MOC group structure */
  if (_cmfd->getNumCmfdGroups() == 0)
    _cmfd->setGroupStructure(NULL, _num_groups+1);

  /* Intialize CMFD Maps */
  _cmfd->initializeCellMap();
  _cmfd->initializeGroupMap();
}


/**
 * @brief Returns the map that maps FSR keys to FSR IDs
 * @return _FSR_keys_map map of FSR keys to FSR IDs
 */
std::unordered_map<std::size_t, fsr_data> Geometry::getFSRKeysMap(){
  return _FSR_keys_map;
}


/**
 * @brief Returns the vector that maps FSR IDs to FSR key hashes
 * @return _FSR_keys_map map of FSR keys to FSR IDs
 */
std::vector<std::size_t> Geometry::getFSRsToKeys(){
  return _FSRs_to_keys;
}


/**
 * @brief Return a vector indexed by flat source region IDs which contain
 *        the corresponding Material IDs.
 * @return an integer vector of FSR-to-Material IDs indexed by FSR ID
 */
std::vector<int> Geometry::getFSRsToMaterialIDs() {
  if (_num_FSRs == 0)
    log_printf(ERROR, "Unable to return the FSR-to-Material map array since "
               "the Geometry has not initialized FSRs.");

  return _FSRs_to_material_IDs;
}


/**
 * @brief Sets the _FSR_keys_map map
 * @param FSR_keys_map map of FSR keys to FSR IDs
 */
void Geometry::setFSRKeysMap(std::unordered_map<std::size_t, fsr_data> FSR_keys_map){
  _FSR_keys_map = FSR_keys_map;
}


/**
 * @brief Sets the _FSRs_to_keys vector
 * @param FSRs_to_keys vector of FSR key hashes indexed by FSR IDs
 */
void Geometry::setFSRsToKeys(std::vector<std::size_t> FSRs_to_keys){
  _FSRs_to_keys = FSRs_to_keys;
}


/**
 * @brief Sets the _FSRs_to_material_IDs vector
 * @param FSRs_to_material_IDs vector mapping FSR IDs to cells
 */
void Geometry::setFSRsToMaterialIDs(std::vector<int> FSRs_to_material_IDs){
  _FSRs_to_material_IDs = FSRs_to_material_IDs;
}


bool Geometry::withinBounds(LocalCoords* coords){

  double x = coords->getX();
  double y = coords->getY();
  
  if (x < _x_min || x > _x_max || y < _y_min || y > _y_max)
    return false;
  else
    return true;
}
