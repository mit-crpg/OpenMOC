#include "Geometry.h"

/**
 * @brief Resets the auto-generated unique IDs for Materials, Surfaces,
 *        Cells and Universes/Lattices to 10000.
 */
void reset_auto_ids() {
  reset_material_id();
  reset_surf_id();
  reset_cell_id();
  reset_universe_id();
}


/**
 * @brief Constructor initializes an empty Geometry.
 */
Geometry::Geometry() {

  /* Initialize CMFD object to NULL */
  _cmfd = NULL;
}


/**
 * @brief Destructor clears FSR to Cells and Materials maps.
 */
Geometry::~Geometry() {

  /* Free FSR maps if they were initialized */
  if (_FSR_keys_map.size() != 0) {

    fsr_data **values = _FSR_keys_map.values();
    for (int i=0; i<_FSR_keys_map.size(); i++)
      delete values[i];
    delete[] values;

    ExtrudedFSR **extruded_fsrs = _extruded_FSR_keys_map.values();
    for (int i=0; i<_extruded_FSR_keys_map.size(); i++) {
      if (extruded_fsrs[i]->_mesh != NULL)
        delete[] extruded_fsrs[i]->_mesh;
      delete[] extruded_fsrs[i]->_fsr_ids;
      delete[] extruded_fsrs[i]->_materials;
      delete extruded_fsrs[i];
    }
    delete[] extruded_fsrs;

    _FSR_keys_map.clear();
    _extruded_FSR_keys_map.clear();
    _FSRs_to_keys.clear();
    _FSRs_to_material_IDs.clear();
  }
}


/**
 * @brief Returns the total width in the x-direction of the Geometry in cm.
 * @return the total width of the Geometry in the x-direction (cm)
 */
double Geometry::getWidthX() {
  return (getMaxX() - getMinX());
}


/**
 * @brief Returns the total width in the y-direction of the Geometry in cm.
 * @return the total width of the Geometry in the y-direction (cm)
 */
double Geometry::getWidthY() {
  return (getMaxY() - getMinY());
}


/**
 * @brief Returns the total width in the z-direction of the Geometry in cm.
 * @return the total width of the Geometry in the z-direction (cm)
 */
double Geometry::getWidthZ() {
  return (getMaxZ() - getMinZ());
}


/**
 * @brief Return the minimum x-coordinate contained by the Geometry.
 * @return the minimum x-coordinate (cm)
 */
double Geometry::getMinX() {
  return _root_universe->getMinX();
}


/**
 * @brief Return the maximum x-coordinate contained by the Geometry.
 * @return the maximum x-coordinate (cm)
 */
double Geometry::getMaxX() {
  return _root_universe->getMaxX();
}


/**
 * @brief Return the minimum y-coordinate contained by the Geometry.
 * @return the minimum y-coordinate (cm)
 */
double Geometry::getMinY() {
  return _root_universe->getMinY();
}


/**
 * @brief Return the maximum y-coordinate contained by the Geometry.
 * @return the maximum y-coordinate (cm)
 */
double Geometry::getMaxY() {
  return _root_universe->getMaxY();
}


/**
 * @brief Return the minimum z-coordinate contained by the Geometry.
 * @return the minimum z-coordinate (cm)
 */
double Geometry::getMinZ() {
  return _root_universe->getMinZ();
}


/**
 * @brief Return the maximum z-coordinate contained by the Geometry.
 * @return the maximum z-coordinate (cm)
 */
double Geometry::getMaxZ() {
  return _root_universe->getMaxZ();
}


/**
 * @brief Returns the boundary conditions (REFLECTIVE or VACUUM) at the
 *        minimum x-coordinate in the Geometry.
 * @return the boundary conditions for the minimum x-coordinate in the Geometry
 */
boundaryType Geometry::getMinXBoundaryType() {
  return _root_universe->getMinXBoundaryType();
}


/**
 * @brief Returns the boundary conditions (REFLECTIVE or VACUUM) at the
 *        maximum x-coordinate in the Geometry.
 * @return the boundary conditions for the maximum z-coordinate in the Geometry
 */
boundaryType Geometry::getMaxXBoundaryType() {
  return _root_universe->getMaxXBoundaryType();
}


/**
 * @brief Returns the boundary conditions (REFLECTIVE or VACUUM) at the
 *        minimum y-coordinate in the Geometry.
 * @return the boundary conditions for the minimum y-coordinate in the Geometry
 */
boundaryType Geometry::getMinYBoundaryType() {
  return _root_universe->getMinYBoundaryType();
}


/**
 * @brief Returns the boundary conditions (REFLECTIVE or VACUUM) at the
 *        maximum y-coordinate in the Geometry.
 * @return the boundary conditions for the maximum y-coordinate in the Geometry
 */
boundaryType Geometry::getMaxYBoundaryType() {
  return _root_universe->getMaxYBoundaryType();
}


/**
 * @brief Returns the boundary conditions (REFLECTIVE or VACUUM) at the
 *        minimum z-coordinate in the Geometry.
 * @return the boundary conditions for the minimum z-coordinate in the Geometry
 */
boundaryType Geometry::getMinZBoundaryType() {
  return _root_universe->getMinZBoundaryType();
}


/**
 * @brief Returns the boundary conditions (REFLECTIVE or VACUUM) at the
 *        maximum z-coordinate in the Geometry.
 * @return the boundary conditions for the maximum z-coordinate in the Geometry
 */
boundaryType Geometry::getMaxZBoundaryType() {
  return _root_universe->getMaxZBoundaryType();
}


/**
 * @brief Returns the number of flat source regions in the Geometry.
 * @return number of FSRs
 */
int Geometry::getNumFSRs() {
  return _FSRs_to_keys.size();
}

/**
 * @brief Returns the number of energy groups for each Material's nuclear data.
 * @return the number of energy groups
 */
int Geometry::getNumEnergyGroups() {

  std::map<int, Material*> materials = getAllMaterials();

  if (materials.size() == 0)
    log_printf(ERROR, "Unable to return the number of energy groups from "
               "the Geometry since it does not contain any Materials");

  int num_groups = materials.begin()->second->getNumEnergyGroups();
  std::map<int, Material*>::iterator iter;

  for (iter = materials.begin(); iter != materials.end(); ++iter) {
    if (iter->second->getNumEnergyGroups() != num_groups)
      log_printf(ERROR, "Unable to return the number of energy groups from "
                 "the Geometry since it contains different numbers of groups: "
                 "%d and %d", num_groups, iter->second->getNumEnergyGroups());
  }

  return num_groups;
}


/**
 * @brief Returns the number of Materials in the Geometry.
 * @return the number of Materials
 */
int Geometry::getNumMaterials() {

  std::map<int, Material*> all_materials;

  if (_all_materials.size() == 0)
    all_materials = getAllMaterials();
  else
    all_materials = _all_materials;

  int num_materials = all_materials.size();
  return num_materials;
}


/**
 * @brief Returns the number of Cells in the Geometry.
 * @return the number of Cells
 */
int Geometry::getNumCells() {

  int num_cells = 0;

  if (_root_universe != NULL) {
    std::map<int, Cell*> all_cells = _root_universe->getAllCells();
    num_cells = all_cells.size();
  }

  return num_cells;
}


/**
 * @brief Return a std::map container of Material IDs (keys) with Materials
 *        pointers (values).
 * @return a std::map of Materials indexed by Material ID in the geometry
 */
std::map<int, Material*> Geometry::getAllMaterials() {

  std::map<int, Material*> all_materials;
  Cell* cell;
  Material* material;

  if (_root_universe != NULL) {
    std::map<int, Cell*> all_cells = _root_universe->getAllCells();
    std::map<int, Cell*>::iterator iter;

    for (iter = all_cells.begin(); iter != all_cells.end(); ++iter) {
      cell = (*iter).second;

      if (cell->getType() == MATERIAL) {
        material = cell->getFillMaterial();

        if (material != NULL)
          all_materials[material->getId()] = material;
      }
    }
  }

  return all_materials;
}


/**
 * @brief Return a std::map container of Cell IDs (keys) with Cells
 *        pointers (values).
 * @return a std::map of Cells indexed by Cell ID in the geometry
 */
std::map<int, Cell*> Geometry::getAllMaterialCells() {

  std::map<int, Cell*> all_material_cells;
  Cell* cell;

  if (_root_universe != NULL) {
    std::map<int, Cell*> all_cells = _root_universe->getAllCells();
    std::map<int, Cell*>::iterator iter;

    for (iter = all_cells.begin(); iter != all_cells.end(); ++iter) {
      cell = (*iter).second;

      if (cell->getType() == MATERIAL)
        all_material_cells[cell->getId()] = cell;
    }
  }

  return all_material_cells;
}


/**
 * @brief Returns the Universe at the root node in the CSG tree.
 * @return the root Universe
 */
Universe* Geometry::getRootUniverse() {
  return _root_universe;
}


/**
 * @brief Returns a pointer to the CMFD object.
 * @return A pointer to the CMFD object
 */
Cmfd* Geometry::getCmfd() {
  return _cmfd;
}


/**
 * @brief Sets the root Universe for the CSG tree.
 * @param root_universe the root Universe of the CSG tree.
 */
void Geometry::setRootUniverse(Universe* root_universe) {
  _root_universe = root_universe;
}


/**
 * @brief Sets the pointer to a CMFD object used for acceleration.
 * @param cmfd a pointer to the CMFD object
 */
void Geometry::setCmfd(Cmfd* cmfd) {
  _cmfd = cmfd;
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
Cell* Geometry::findCellContainingCoords(LocalCoords* coords) {

  Universe* univ = coords->getUniverse();
  Cell* cell;

  if (univ->getId() == _root_universe->getId()) {
    if (!withinBounds(coords))
      return NULL;
  }

  if (univ->getType() == SIMPLE)
    cell = univ->findCell(coords);
  else
    cell = static_cast<Lattice*>(univ)->findCell(coords);

  return cell;
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
Cell* Geometry::findFirstCell(LocalCoords* coords, double azim, double polar) {
  double delta_x = cos(azim) * sin(polar) * TINY_MOVE;
  double delta_y = sin(azim) * sin(polar) * TINY_MOVE;
  double delta_z = cos(polar) * TINY_MOVE;
  coords->adjustCoords(delta_x, delta_y, delta_z);
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

  std::map<int, Material*> all_materials;

  if (_all_materials.size() == 0)
    all_materials = getAllMaterials();
  else
    all_materials = _all_materials;

  return all_materials[_FSRs_to_material_IDs.at(fsr_id)];
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
Cell* Geometry::findNextCell(LocalCoords* coords, double azim, double polar) {

  Cell* cell = NULL;
  double dist;
  double min_dist = std::numeric_limits<double>::infinity();

  /* Get lowest level coords */
  coords = coords->getLowestLevel();

  /* Get the current Cell */
  cell = coords->getCell();

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
        Lattice* lattice = coords->getLattice();
        dist = lattice->minSurfaceDist(coords->getPoint(), azim, polar);
      }
      /* If we reach a LocalCoord in a Universe, find the distance to the
       * nearest cell surface */
      else{
        Cell* cell = coords->getCell();
        dist = cell->minSurfaceDist(coords->getPoint(), azim, polar);
      }

      /* Recheck min distance */
      min_dist = std::min(dist, min_dist);

      /* Ascend one level */
      if (coords->getUniverse() == _root_universe)
        break;
      else {
        coords = coords->getPrev();
        coords->prune();
      }
    }

    /* Check for distance to nearest CMFD mesh cell boundary */
    if (_cmfd != NULL) {
      Lattice* lattice = _cmfd->getLattice();
      dist = lattice->minSurfaceDist(coords->getPoint(), azim, polar);
      min_dist = std::min(dist, min_dist);
    }

    /* Move point and get next cell */
    double delta_x = cos(azim) * sin(polar) * (min_dist + TINY_MOVE);
    double delta_y = sin(azim) * sin(polar) * (min_dist + TINY_MOVE);
    double delta_z = cos(polar) * (min_dist + TINY_MOVE);
    coords->adjustCoords(delta_x, delta_y, delta_z);

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

  int fsr_id;
  LocalCoords* curr = coords;
  curr = coords->getLowestLevel();
  std::hash<std::string> key_hash_function;

  /* Generate faulty FSR key */
  std::size_t fsr_key_hash = key_hash_function(getFSRKey(coords));

  /* If FSR has not been encountered, update FSR maps and vectors */
  if (!_FSR_keys_map.contains(fsr_key_hash)) {

    /* Try to get a clean copy of the fsr_id, adding the FSR data
       if necessary where -1 indicates the key was already added */
    fsr_id = _FSR_keys_map.insert_and_get_count(fsr_key_hash, NULL);
    if (fsr_id == -1) {
      fsr_data volatile* fsr;
      do {
        fsr = _FSR_keys_map.at(fsr_key_hash);
      } while (fsr == NULL);
      fsr_id = fsr->_fsr_id;
    }
    else {

      /* Add FSR information to FSR key map and FSR_to vectors */
      fsr_data* fsr = new fsr_data;
      fsr->_fsr_id = fsr_id;
      _FSR_keys_map.at(fsr_key_hash) = fsr;
      Point* point = new Point();
      point->setCoords(coords->getHighestLevel()->getX(),
                       coords->getHighestLevel()->getY(),
                       coords->getHighestLevel()->getZ());

      /* Get the cell that contains coords */
      Cell* cell = findCellContainingCoords(curr);
      fsr->_point = point;
      fsr->_mat_id = cell->getFillMaterial()->getId();

      /* If CMFD acceleration is on, add FSR CMFD cell to FSR data */
      if (_cmfd != NULL)
        fsr->_cmfd_cell = _cmfd->findCmfdCell(coords->getHighestLevel());
    }
  }

  /* If FSR has already been encountered, get the fsr id from map */
  else {
    fsr_data volatile* fsr;
    do {
      fsr = _FSR_keys_map.at(fsr_key_hash);
    } while (fsr == NULL);

    fsr_id = fsr->_fsr_id;
  }

  return fsr_id;
}


/**
 * @brief Finds and returns a pointer to the axially extruded flat source
 *        region that a given LocalCoords object resides within.
 * @param coords a LocalCoords object pointer
 * @return the ID of an extruded FSR for a given LocalCoords object
 */
int Geometry::findExtrudedFSR(LocalCoords* coords) {

  int fsr_id;
  LocalCoords* curr = coords;
  curr = coords->getLowestLevel();
  std::hash<std::string> key_hash_function;

  /* Generate unique FSR key */
  std::size_t fsr_key_hash = key_hash_function(getFSRKey(coords));

  /* If FSR has not been encountered, update FSR maps and vectors */
  if (!_extruded_FSR_keys_map.contains(fsr_key_hash)) {

    /* Try to get a clean copy of the fsr_id, adding the FSR data
       if necessary where -1 indicates the key was already added */
    fsr_id = _extruded_FSR_keys_map.insert_and_get_count(fsr_key_hash, NULL);
    if (fsr_id == -1) {
      ExtrudedFSR volatile* fsr;
      do {
        fsr = _extruded_FSR_keys_map.at(fsr_key_hash);
      } while (fsr == NULL);
      fsr_id = fsr->_fsr_id;
    }
    else {

      /* Add FSR information to FSR key map and FSR_to vectors */
      ExtrudedFSR* fsr = new ExtrudedFSR;
      fsr->_fsr_id = fsr_id;
      fsr->_num_fsrs = 0;
      fsr->_coords = new LocalCoords(0, 0, 0);
      _extruded_FSR_keys_map.at(fsr_key_hash) = fsr;
      coords->copyCoords(fsr->_coords);
    }
  }

  /* If FSR has already been encountered, get the fsr id from map */
  else {
    ExtrudedFSR volatile* fsr;
    do {
      fsr = _extruded_FSR_keys_map.at(fsr_key_hash);
    } while (fsr == NULL);

    fsr_id = fsr->_fsr_id;
  }

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
    fsr_id = _FSR_keys_map.at(key_hash_function(fsr_key))->_fsr_id;
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not find FSR ID with key: %s. Try creating "
               "geometry with finer track spacing", fsr_key.c_str());
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
    point = _FSR_keys_map.at(_FSRs_to_keys.at(fsr_id))->_point;
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not find characteristic point in FSR: %d", fsr_id);
  }

  return point;
}


/**
 * @brief Return the centroid for a given FSR ID
 * @param fsr_id the FSR ID
 * @return the FSR's centroid
 */
Point* Geometry::getFSRCentroid(int fsr_id) {

  Point* point;

  try{
    point = _FSR_keys_map.at(_FSRs_to_keys.at(fsr_id))->_centroid;
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not find centroid in FSR: %d.", fsr_id);
  }

  return point;
}


/**
 * @brief Return the CMFD cell for a given FSR ID
 * @param fsr_id the FSR ID
 * @return the CMFD cell
 */
int Geometry::getCmfdCell(int fsr_id) {

  int cmfd_cell;

  try{
    cmfd_cell = _FSR_keys_map.at(_FSRs_to_keys.at(fsr_id))->_cmfd_cell;
  }
  catch(std::exception &e) {
    log_printf(ERROR, "Could not find CMFD cell in FSR: %d", fsr_id);
  }

  return cmfd_cell;
}

/**
 * @brief Generate a string FSR "key" that identifies an FSR by its
 *        unique hierarchical lattice/universe/cell structure.
 * @details Since not all FSRs will reside on the absolute lowest universe
 *          level and Cells might overlap other cells, it is important to
 *          have a method for uniquely identifying FSRs. This method
 *          creates a unique FSR key by constructing a structured string
 *          that describes the hierarchy of lattices/universes/cells.
 * @param coords a LocalCoords object pointer
 * @return the FSR key
 */
std::string Geometry::getFSRKey(LocalCoords* coords) {

  std::stringstream key;
  LocalCoords* curr = coords->getHighestLevel();
  std::ostringstream curr_level_key;

  /* If CMFD is on, get CMFD latice cell and write to key */
  if (_cmfd != NULL) {
      curr_level_key << _cmfd->getLattice()->getLatX(curr->getPoint());
      key << "CMFD = (" << curr_level_key.str() << ", ";
      curr_level_key.str(std::string());
      curr_level_key << _cmfd->getLattice()->getLatY(curr->getPoint());
      key << curr_level_key.str() << ", ";
      curr_level_key.str(std::string());
      curr_level_key << _cmfd->getLattice()->getLatZ(curr->getPoint());
      key << curr_level_key.str() << ") : ";
  }

  /* Descend the linked list hierarchy until the lowest level has
   * been reached */
  while (curr != NULL) {

    /* Clear string stream */
    curr_level_key.str(std::string());

    if (curr->getType() == LAT) {

      /* Write lattice ID and lattice cell to key */
      curr_level_key << curr->getLattice()->getId();
      key << "LAT = " << curr_level_key.str() << " (";
      curr_level_key.str(std::string());
      curr_level_key << curr->getLatticeX();
      key << curr_level_key.str() << ", ";
      curr_level_key.str(std::string());
      curr_level_key << curr->getLatticeY();
      key << curr_level_key.str() << ", ";
      curr_level_key.str(std::string());
      curr_level_key << curr->getLatticeZ();
      key << curr_level_key.str() << ") : ";
    }
    else {
      /* write universe ID to key */
      curr_level_key << curr->getUniverse()->getId();
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
  curr_level_key << curr->getCell()->getId();
  key << "CELL = " << curr_level_key.str();

  return key.str();
}


/**
 * @brief Return a pointer to an ExtrudedFSR by its extruded FSR ID
 * @param extruded_fsr_id the extruded FSR ID
 * @return a pointer to the ExtrudedFSR
 */
ExtrudedFSR* Geometry::getExtrudedFSR(int extruded_fsr_id) {
  return _extruded_FSR_lookup[extruded_fsr_id];
}


/**
 * @brief Subidivides all Cells in the Geometry into rings and angular sectors
 *        aligned with the z-axis.
 * @details This method is called by the Geometry::initializeFlatSourceRegions()
 *          method but may also be called by the user in Python if needed:
 *
 * @code
 *          geometry.subdivideCells()
 * @endcode
 */
void Geometry::subdivideCells() {

  /* Compute equivalent radius with the same area as the Geometry */
  /* This is used as the maximum radius for all ringified Cells */
  double max_radius = sqrt(getWidthX() * getWidthY() / M_PI);

  /* Recursively subdivide Cells into rings and sectors */
  _root_universe->subdivideCells(max_radius);
}


/**
 * @brief Compute the number of flat source regions in the Geometry and
 *        initialize CMFD.
 * @details This method is intended to be called by the user before initiating
 *          source iteration. This method first subdivides all Cells by calling
 *          the Geometry::subdivideCells() method. Then it initializes the CMFD
 *          object.
 */
void Geometry::initializeFlatSourceRegions() {

  /* Subdivide Cells into sectors and rings */
  subdivideCells();

  /* Build collections of neighbor Cells for optimized ray tracing */
  _root_universe->buildNeighbors();

  /* Create map of Material IDs to Material pointers */
  _all_materials = getAllMaterials();

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
 * @param z_coord the axial height at which the 2D plane of the geometry is
 *        formed
 */
void Geometry::segmentize2D(Track* track, double z_coord) {

  /* Track starting Point coordinates and azimuthal angle */
  double x0 = track->getStart()->getX();
  double y0 = track->getStart()->getY();
  double z0 = z_coord;
  double phi = track->getPhi();
  double delta_x, delta_y, delta_z;

  /* Length of each segment */
  FP_PRECISION length;
  Material* material;
  int fsr_id;
  int num_segments;

  /* Use a LocalCoords for the start and end of each segment */
  LocalCoords start(x0, y0, z0);
  LocalCoords end(x0, y0, z0);
  start.setUniverse(_root_universe);
  end.setUniverse(_root_universe);

  /* Find the Cell containing the Track starting Point */
  Cell* curr = findFirstCell(&end, phi);
  Cell* prev;

  /* If starting Point was outside the bounds of the Geometry */
  if (curr == NULL)
    log_printf(ERROR, "Could not find a Cell containing the start Point "
               "of this Track: %s", track->toString().c_str());

  /* While the end of the segment's LocalCoords is still within the Geometry,
   * move it to the next Cell, create a new segment, and add it to the
   * Geometry */
  while (curr != NULL) {

    end.copyCoords(&start);

    /* Find the next Cell along the Track's trajectory */
    prev = curr;
    curr = findNextCell(&end, phi);

    /* Checks that segment does not have the same start and end Points */
    if (start.getX() == end.getX() && start.getY() == end.getY())
      log_printf(ERROR, "Created segment with same start and end "
                 "point: x = %f, y = %f", start.getX(), start.getY());

    /* Find the segment length, Material and FSR ID */
    length = FP_PRECISION(end.getPoint()->distanceToPoint(start.getPoint()));
    material = prev->getFillMaterial();
    fsr_id = findFSRId(&start);

    /* Create a new Track segment */
    segment* new_segment = new segment;
    new_segment->_material = material;
    new_segment->_length = length;
    new_segment->_region_id = fsr_id;

    log_printf(DEBUG, "segment start x = %f, y = %f; end x = %f, y = %f",
               start.getX(), start.getY(), end.getX(), end.getY());

    /* Save indicies of CMFD Mesh surfaces that the Track segment crosses */
    if (_cmfd != NULL) {

      /* Find cmfd cell that segment lies in */
      int cmfd_cell = _cmfd->findCmfdCell(&start);

      /* Reverse nudge from surface to determine whether segment start or end
       * points lie on a CMFD surface. */
      delta_x = cos(phi) * TINY_MOVE;
      delta_y = sin(phi) * TINY_MOVE;
      start.adjustCoords(-delta_x, -delta_y);
      end.adjustCoords(-delta_x, -delta_y);

      new_segment->_cmfd_surface_fwd =
        _cmfd->findCmfdSurface(cmfd_cell, &end);
      new_segment->_cmfd_surface_bwd =
        _cmfd->findCmfdSurface(cmfd_cell, &start);

      /* Re-nudge segments from surface. */
      start.adjustCoords(delta_x, delta_y);
      end.adjustCoords(delta_x, delta_y);
    }

    /* Add the segment to the Track */
    track->addSegment(new_segment);

  }

  log_printf(DEBUG, "Created %d segments for Track: %s",
             track->getNumSegments(), track->toString().c_str());

  /* Truncate the linked list for the LocalCoords */
  start.prune();
  end.prune();

  return;
}


/**
 * @brief This method performs ray tracing to create Track segments within each
 *        flat source region in the Geometry.
 * @details This method starts at the beginning of a Track and finds successive
 *          intersection points with FSRs as the Track crosses through the
 *          Geometry and creates segment structs and adds them to the Track.
 * @param track a pointer to a track to segmentize
 * @param max_optical_length the maximum optical length a segment is allowed to
 *          have
 */
void Geometry::segmentize3D(Track3D* track) {

  /* Track starting Point coordinates and azimuthal angle */
  double x0 = track->getStart()->getX();
  double y0 = track->getStart()->getY();
  double z0 = track->getStart()->getZ();
  double phi = track->getPhi();
  double theta = track->getTheta();
  double delta_x, delta_y, delta_z;

  /* Length of each segment */
  FP_PRECISION length;
  Material* material;
  int fsr_id;
  int num_segments;

  /* Use a LocalCoords for the start and end of each segment */
  LocalCoords start(x0, y0, z0);
  LocalCoords end(x0, y0, z0);
  start.setUniverse(_root_universe);
  end.setUniverse(_root_universe);

  /* Find the Cell containing the Track starting Point */
  Cell* curr = findFirstCell(&end, phi, theta);
  Cell* prev;

  /* If starting Point was outside the bounds of the Geometry */
  if (curr == NULL)
    log_printf(ERROR, "Could not find a Cell containing the start Point "
               "of this Track3D: %s", track->toString().c_str());

  /* While the end of the segment's LocalCoords is still within the Geometry,
   * move it to the next Cell, create a new segment, and add it to the
   * Geometry */
  while (curr != NULL) {

    end.copyCoords(&start);

    /* Find the next Cell along the Track's trajectory */
    prev = curr;
    curr = findNextCell(&end, phi, theta);

    /* Checks to make sure that new Segment does not have the same start
     * and end Points */
    if (start.getX() == end.getX() &&
        start.getY() == end.getY() &&
        start.getZ() == end.getZ()) {

      log_printf(ERROR, "Created a Track3D segment with the same start and end "
                 "point: x = %f, y = %f, z = %f", start.getX(),
                 start.getY(), start.getZ());
    }

    /* Find the segment length between the segment's start and end points */
    length = FP_PRECISION(end.getPoint()->distanceToPoint(start.getPoint()));
    material = prev->getFillMaterial();
    fsr_id = findFSRId(&start);

    /* Create a new Track segment */
    segment* new_segment = new segment;
    new_segment->_material = material;
    new_segment->_length = length;
    new_segment->_region_id = fsr_id;

    log_printf(DEBUG, "segment start x = %f, y = %f, z = %f; "
               "end x = %f, y = %f, z = %f",
               start.getX(), start.getY(), start.getZ(),
               end.getX(), end.getY(), end.getZ());

    /* Save indicies of CMFD Mesh surfaces that the Track segment crosses */
    if (_cmfd != NULL) {

      /* Find cmfd cell that segment lies in */
      int cmfd_cell = _cmfd->findCmfdCell(&start);

      /* Reverse nudge from surface to determine whether segment start or end
       * points lie on a cmfd surface. */
      double delta_x = cos(phi) * sin(theta) * TINY_MOVE;
      double delta_y = sin(phi) * sin(theta) * TINY_MOVE;
      double delta_z = cos(theta) * TINY_MOVE;
      start.adjustCoords(-delta_x, -delta_y, -delta_z);
      end.adjustCoords(-delta_x, -delta_y, -delta_z);

      new_segment->_cmfd_surface_fwd =
        _cmfd->findCmfdSurface(cmfd_cell, &end);
      new_segment->_cmfd_surface_bwd =
        _cmfd->findCmfdSurface(cmfd_cell, &start);

      /* Re-nudge segments from surface. */
      start.adjustCoords(delta_x, delta_y, delta_z);
      end.adjustCoords(delta_x, delta_y, delta_z);
    }

    /* Add the segment to the Track */
    track->addSegment(new_segment);
  }

  log_printf(DEBUG, "Created %d segments for Track3D: %s",
             track->getNumSegments(), track->toString().c_str());

  /* Truncate the linked list for the LocalCoords */
  start.prune();
  end.prune();
}


/**
 * @brief This method performs ray tracing to create extruded track segments
 *        within each flat source region in the implicit 2D superposition plane
 *        of the Geometry by 2D ray tracing across input heights that encompass
 *        all radial geometric detail.
 * @details This method starts at the beginning of an extruded track and finds
 *          successive intersection points with FSRs as the extruded track
 *          crosses radially through the Geometry at defined z-coords. The
 *          minimum distance to intersection of all z-coords is chosen leading
 *          to implicitly capturing all geometric radial detail at the defined
 *          z-heights, saving the lengths and region IDs to the extruded track
 *          and initializing ExtrudedFSR structs in the traversed FSRs.
 * @param flattend_track a pointer to a 2D track to segmentize into regions of
 *        extruded FSRs
 * @param z_coords a vector of axial heights in the root geometry at which
 *        the Geometry is segmentized radially
 */
void Geometry::segmentizeExtruded(Track* flattened_track,
    std::vector<double> z_coords) {

  /* Track starting Point coordinates and azimuthal angle */
  double x0 = flattened_track->getStart()->getX();
  double y0 = flattened_track->getStart()->getY();
  double z0 = z_coords[0];
  double phi = flattened_track->getPhi();
  double delta_x, delta_y, delta_z;

  /* Length of each segment */
  FP_PRECISION length;
  int min_z_ind = 0;
  int region_id;
  int num_segments;

  /* Use a LocalCoords for the start and end of each segment */
  LocalCoords start(x0, y0, z0);
  LocalCoords end(x0, y0, z0);
  start.setUniverse(_root_universe);
  end.setUniverse(_root_universe);

  /* Find the Cell containing the Track starting Point */
  Cell* curr = findFirstCell(&end, phi);

  /* If starting Point was outside the bounds of the Geometry */
  if (curr == NULL)
    log_printf(ERROR, "Could not find a Cell containing the start Point "
               "of this Track: %s", flattened_track->toString().c_str());

  /* While the end of the segment's LocalCoords is still within the Geometry,
   * move it to the next Cell, create a new segment, and add it to the
   * Geometry */
  while (curr != NULL) {

    /* Records the minimum length to a 2D intersection */
    FP_PRECISION min_length = std::numeric_limits<FP_PRECISION>::infinity();

    /* Copy end coordinates to start */
    end.copyCoords(&start);

    /* Loop over all z-heights to find shortest 2D intersection */
    for (int i=0; i < z_coords.size(); i++) {

      /* Change z-height and copy starting coordinates to end */
      start.setZ(z_coords[i]);
      start.copyCoords(&end);

      /* Find the next Cell along the Track's trajectory */
      curr = findNextCell(&end, phi);

      /* Checks that segment does not have the same start and end Points */
      if (start.getX() == end.getX() && start.getY() == end.getY())
        log_printf(ERROR, "Created segment with same start and end "
                   "point: x = %f, y = %f", start.getX(), start.getY());

      /* Find the segment length and extruded FSR */
      length = FP_PRECISION(end.getPoint()->distanceToPoint(start.getPoint()));

      /* Check if the segment length is the smallest found */
      if (length < min_length) {
        min_length = length;
        min_z_ind = i;
      }
    }

    /* Traverse across shortest segment */
    start.setZ(z_coords[min_z_ind]);
    start.copyCoords(&end);
    curr = findNextCell(&end, phi);

    /* Find FSR using starting coordinate */
    region_id = findExtrudedFSR(&start);

    log_printf(DEBUG, "segment start x = %f, y = %f; end x = %f, y = %f",
               start.getX(), start.getY(), end.getX(), end.getY());

    /* Create a new 2D Track segment with extruded region ID */
    segment* new_segment = new segment;
    new_segment->_length = min_length;
    new_segment->_region_id = region_id;

    /* Save indicies of CMFD Mesh surfaces that the Track segment crosses */
    if (_cmfd != NULL) {

      /* Find cmfd cell that segment lies in */
      int cmfd_cell = _cmfd->findCmfdCell(&start);

      /* Reverse nudge from surface to determine whether segment start or end
       * points lie on a CMFD surface. */
      delta_x = cos(phi) * TINY_MOVE;
      delta_y = sin(phi) * TINY_MOVE;
      start.adjustCoords(-delta_x, -delta_y, 0);
      end.adjustCoords(-delta_x, -delta_y, 0);

      new_segment->_cmfd_surface_fwd =
        _cmfd->findCmfdSurface(cmfd_cell, &end);
      new_segment->_cmfd_surface_bwd =
        _cmfd->findCmfdSurface(cmfd_cell, &start);

      /* Re-nudge segments from surface. */
      start.adjustCoords(delta_x, delta_y, 0);
      end.adjustCoords(delta_x, delta_y, 0);
    }

    /* Add the segment to the 2D track */
    flattened_track->addSegment(new_segment);
  }

  /* Truncate the linked list for the LocalCoords */
  start.prune();
  end.prune();

  return;
}


/**
 * @brief Rays are shot vertically through each ExtrudedFSR struct to calculate
 *        the axial mesh and initialize 3D FSRs
 * @details From a 2D point within each FSR, a temporary 3D track is created
 *          starting at the bottom of the geometry and extending vertically to
 *          the top of the geometry. These tracks are segmented using the
 *          segmentize3D routine to calculate the distances between axial
 *          intersections forming the axial meshes if necessary and
 *          initializing the 3D FSRs as new regions are traversed.
 * @param global_z_mesh A global z mesh used for ray tracing. If the vector's
 *        length is zero, z meshes are local and need to be created for every
 *        ExtrudedFSR.
 */
void Geometry::initializeAxialFSRs(std::vector<FP_PRECISION> global_z_mesh) {

  log_printf(NORMAL, "Initializing 3D FSRs in axially extruded regions...");

  /* Determine the extent of the axial geometry */
  FP_PRECISION min_z = getMinZ();
  FP_PRECISION max_z = getMaxZ();

  /* Extract list of extruded FSRs */
  ExtrudedFSR** extruded_FSRs = _extruded_FSR_keys_map.values();

  std::string msg = "initializing 3D FSRs";
  Progress progress(_extruded_FSR_keys_map.size(), msg);

  /* Loop over extruded FSRs */
  for (int i=0; i < _extruded_FSR_keys_map.size(); i++) {

    progress.incrementCounter();

    /* Extract coordinates of extruded FSR */
    ExtrudedFSR* extruded_FSR = extruded_FSRs[i];
    double x0 = extruded_FSR->_coords->getX();
    double y0 = extruded_FSR->_coords->getY();

    /* Determine if there is a global mesh or local meshes should be created */
    if (global_z_mesh.size() > 0) {

      /* Allocate materials in extruded FSR */
      size_t num_regions = global_z_mesh.size() - 1;
      extruded_FSR->_num_fsrs = num_regions;
      extruded_FSR->_materials = new Material*[num_regions];
      extruded_FSR->_fsr_ids = new int[num_regions];

      /* Loop over all regions in the global mesh */
      for (int n=0; n < num_regions; n++) {

        /* Set the axial coordinate at the midpoint of mesh boundaries */
        double midpt = (global_z_mesh[n] + global_z_mesh[n+1]) / 2;
        LocalCoords coord(x0, y0, midpt);
        coord.setUniverse(_root_universe);

        /* Get the FSR ID and material */
        Cell* cell = findCellContainingCoords(&coord);
        int fsr_id = findFSRId(&coord);
        Material* material = cell->getFillMaterial();

        /* Set the FSR ID and material */
        extruded_FSR->_fsr_ids[n] = fsr_id;
        extruded_FSR->_materials[n] = material;
      }
    }
    else {

      /* Create vertical track in the extruded FSR */
      Track3D track;
      track.setValues(x0, y0, min_z, x0, y0, max_z, 0, 0);

      /* Shoot vertical track through the geometry to initialize 3D FSRs */
      segmentize3D(&track);

      /* Extract segments from track */
      int num_segments = track.getNumSegments();
      segment* segments = track.getSegments();

      /* Allocate materials and mesh in extruded FSR */
      extruded_FSR->_num_fsrs = (size_t) num_segments;
      extruded_FSR->_materials = new Material*[num_segments];
      extruded_FSR->_fsr_ids = new int[num_segments];
      extruded_FSR->_mesh = new FP_PRECISION[num_segments+1];

      /* Initialize values in extruded FSR */
      for (int s=0; s < num_segments; s++) {
        extruded_FSR->_materials[s] = segments[s]._material;
        extruded_FSR->_fsr_ids[s] = segments[s]._region_id;
      }

      /* Initialize z mesh */
      FP_PRECISION level = min_z;
      extruded_FSR->_mesh[0] = level;
      for (int s=0; s < num_segments; s++) {
        level += segments[s]._length;
        if (std::abs(level - max_z) < 1e-12)
          level = max_z;
      extruded_FSR->_mesh[s+1] = level;
      }
    }
  }
  delete[] extruded_FSRs;
}


/**
 * @brief Initialize key and material ID vectors for lookup by FSR ID
 * @detail This function initializes and sets reverse lookup vectors by FSR ID.
 *      This is called after the FSRs have all been identified and allocated
 *      during segmentation. This function must be called after
 *      Geometry::segmentize() has completed. It should not be called if tracks
 *      are loaded from a file.
 */
void Geometry::initializeFSRVectors() {

  /* get keys and values from map */
  std::size_t *key_list = _FSR_keys_map.keys();
  fsr_data **value_list = _FSR_keys_map.values();

  /* allocate vectors */
  int num_FSRs = _FSR_keys_map.size();
  _FSRs_to_keys = std::vector<std::size_t>(num_FSRs);
  _FSRs_to_material_IDs = std::vector<int>(num_FSRs);

  /* fill vectors key and material ID information */
  #pragma omp parallel for
  for (int i=0; i < num_FSRs; i++)
  {
    std::size_t key = key_list[i];
    fsr_data* fsr = value_list[i];
    int fsr_id = fsr->_fsr_id;
    _FSRs_to_keys.at(fsr_id) = key;
    _FSRs_to_material_IDs.at(fsr_id) = fsr->_mat_id;
  }

  /* add cmfd information serially */
  if (_cmfd != NULL) {
    for (int i=0; i < num_FSRs; i++) {
      fsr_data* fsr = value_list[i];
      int fsr_id = fsr->_fsr_id;
      _cmfd->addFSRToCell(fsr->_cmfd_cell, fsr_id);
    }
  }

  /* Check if extruded FSRs are present */
  size_t num_extruded_FSRs = _extruded_FSR_keys_map.size();
  if (num_extruded_FSRs > 0) {

    /* Allocate extruded FSR lookup vector and fill with extruded FSRs by ID */
    _extruded_FSR_lookup = std::vector<ExtrudedFSR*>(num_extruded_FSRs);
    ExtrudedFSR **extruded_value_list = _extruded_FSR_keys_map.values();
    for (int i=0; i < num_extruded_FSRs; i++) {
      int fsr_id = extruded_value_list[i]->_fsr_id;
      _extruded_FSR_lookup[fsr_id] = extruded_value_list[i];
    }

    delete [] extruded_value_list;
  }

  /* Delete key and value lists */
  delete[] key_list;
  delete[] value_list;
}


/**
 * @brief Determines the fissionability of each Universe within this Geometry.
 * @details A Universe is determined fissionable if it contains a Cell
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

  Material* material;
  std::map<int, Material*> materials;
	std::map<int, Material*>::iterator mat_iter;

  Universe* universe;
  std::map<int, Universe*> universes;
	std::map<int, Universe*>::iterator univ_iter;

  /* If no Universe was passed in as an argument, then this is the first
   * recursive call from a user via Python, so get the base Universe */
  if (univ == NULL)
    univ = _root_universe;

  /* If a Universe was passed in as an argument, then this is a recursive
   * call with a Universe at a lower level in the nested Universe hierarchy */
  if (univ->getType() == SIMPLE) {
    materials = univ->getAllMaterials();
    universes = univ->getAllUniverses();
  }

  else
    universes = static_cast<Lattice*>(univ)->getAllUniverses();

  /* Loop over the nested Universes first to ensure that fissionability
   * is set at each nested Universe level */
  for (univ_iter=universes.begin(); univ_iter != universes.end(); ++univ_iter) {
    universe = univ_iter->second;

    /* Recursively check whether this nested Universe is fissionable */
    computeFissionability(universe);

    if (universe->isFissionable())
      fissionable = true;
  }

  /* Loop over the Materials in this Universe at this level */
  for (mat_iter=materials.begin(); mat_iter != materials.end(); ++mat_iter) {
    material = mat_iter->second;

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

  std::map<int, Cell*> all_cells = _root_universe->getAllCells();
  std::map<int, Universe*> all_universes = _root_universe->getAllUniverses();

  std::map<int, Cell*>::iterator cell_iter;
  std::map<int, Universe*>::iterator univ_iter;

  string << "\n\tCells:\n\t\t";
  for (cell_iter = all_cells.begin(); cell_iter != all_cells.end(); ++cell_iter)
    string << cell_iter->second->toString() << "\n\t\t";

  string << "\n\tUniverses:\n\t\t";
  for (univ_iter = all_universes.begin();
       univ_iter != all_universes.end(); ++univ_iter)
    string << univ_iter->second->toString() << "\n\t\t";

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
void Geometry::initializeCmfd() {

  /* Set CMFD mesh boundary conditions */
  _cmfd->setBoundary(SURFACE_X_MIN, getMinXBoundaryType());
  _cmfd->setBoundary(SURFACE_Y_MIN, getMinYBoundaryType());
  _cmfd->setBoundary(SURFACE_Z_MIN, getMinZBoundaryType());
  _cmfd->setBoundary(SURFACE_X_MAX, getMaxXBoundaryType());
  _cmfd->setBoundary(SURFACE_Y_MAX, getMaxYBoundaryType());
  _cmfd->setBoundary(SURFACE_Z_MAX, getMaxZBoundaryType());

  /* Set CMFD mesh dimensions and number of groups */
  _cmfd->setWidthX(getWidthX());
  _cmfd->setWidthY(getWidthY());

  /* Intialize CMFD Maps */
  _cmfd->initializeCellMap();

  /* Initialize the CMFD lattice */
  Point offset;
  offset.setX(getMinX() + getWidthX()/2.0);
  offset.setY(getMinY() + getWidthY()/2.0);

  /* If geometry is infinite in z, set Cmfd z-width to 1.0 and z-offset to 0 */
  if (getWidthZ() == std::numeric_limits<double>::infinity()) {
    _cmfd->setWidthZ(1.0);
    offset.setZ(0.0);
  }
  else {
    _cmfd->setWidthZ(getWidthZ());
    offset.setZ(getMinZ() + getWidthZ()/2.0);
  }

  _cmfd->initializeLattice(&offset);
}


/**
 * @brief Returns a pointer to the map that maps FSR keys to FSR IDs
 * @return pointer to _FSR_keys_map map of FSR keys to FSR IDs
 */
ParallelHashMap<std::size_t, fsr_data*>* Geometry::getFSRKeysMap() {
  return &_FSR_keys_map;
}


/**
 * @brief Returns the vector that maps FSR IDs to FSR key hashes
 * @return _FSR_keys_map map of FSR keys to FSR IDs
 */
std::vector<std::size_t>* Geometry::getFSRsToKeys() {
  return &_FSRs_to_keys;
}

/**
 * @brief Return a vector indexed by flat source region IDs which contain
 *        the corresponding Material IDs.
 * @return an integer vector of FSR-to-Material IDs indexed by FSR ID
 */
std::vector<int>* Geometry::getFSRsToMaterialIDs() {
  if (_FSR_keys_map.size() == 0)
    log_printf(ERROR, "Unable to return the FSR-to-Material map array since "
               "the Geometry has not initialized FSRs.");

  return &_FSRs_to_material_IDs;
}

/**
 * @brief Sets the _FSR_keys_map map
 * @details The _FSR_keys_map stores a hash of a std::string representing
 *          the Lattice/Cell/Universe hierarchy for a unique region
 *          and the associated FSR data. fsr_data is a struct that contains
 *          a unique FSR id and a Point located in the highest level Universe
 *          that is contained in the FSR. This method is used when the tracks
 *          are read from file to avoid unnecessary segmentation.
 * @param FSR_keys_map map of FSR keys to FSR data
 */
void Geometry::setFSRKeysMap(ParallelHashMap<std::size_t, fsr_data*>*
                             FSR_keys_map) {
  _FSR_keys_map = *FSR_keys_map;
}

/**
 * @brief Sets the _FSRs_to_keys vector
 * @param FSRs_to_keys vector of FSR key hashes indexed by FSR IDs
 */
void Geometry::setFSRsToKeys(std::vector<std::size_t>* FSRs_to_keys) {
  _FSRs_to_keys = *FSRs_to_keys;
}


/**
 * @brief Sets the _FSRs_to_material_IDs vector
 * @param FSRs_to_material_IDs vector mapping FSR IDs to cells
 */
void Geometry::setFSRsToMaterialIDs(std::vector<int>* FSRs_to_material_IDs) {
  _FSRs_to_material_IDs = *FSRs_to_material_IDs;
}


/**
 * @brief Determins whether a point is within the bounding box of the geometry.
 * @param coords a populated LocalCoords linked list
 * @return boolean indicating whether the coords is within the geometry
 */
bool Geometry::withinBounds(LocalCoords* coords) {

  double x = coords->getX();
  double y = coords->getY();
  double z = coords->getZ();

  if (x < getMinX() || x > getMaxX() || y < getMinY() || y > getMaxY()
      || z < getMinZ() || z > getMaxZ())
    return false;
  else
    return true;
}


/**
 * @brief Finds the Cell containing a given fsr ID.
 * @param fsr_id an FSR ID.
 */
Cell* Geometry::findCellContainingFSR(int fsr_id) {

  Point* point = _FSR_keys_map.at(_FSRs_to_keys[fsr_id])->_point;
  LocalCoords* coords = new LocalCoords(point->getX(), point->getY(),
                                        point->getZ());
  coords->setUniverse(_root_universe);
  Cell* cell = findCellContainingCoords(coords);

  delete coords;

  return cell;
}


/**
 * @brief Sets the centroid for an FSR
 * @details The _FSR_keys_map stores a hash of a std::string representing
 *          the Lattice/Cell/Universe hierarchy for a unique region
 *          and the associated FSR data. _centroid is a point that represents
 *          the numerical centroid of an FSR computed using all segments
 *          contained in the FSR. This method is used by the TrackGenerator
 *          to set the centroid after segments have been created. It is
 *          important to note that this method is a helper function for the
 *          TrackGenerator and should not be explicitly called by the user.
 * @param fsr a FSR id
 * @param centroid a Point representing the FSR centroid
 */
void Geometry::setFSRCentroid(int fsr, Point* centroid) {
  _FSR_keys_map.at(_FSRs_to_keys[fsr])->_centroid = centroid;
}


/*
 * @brief Returns a vector of z-coords defining a superposition of all axial
 *        boundaries in the Geometry.
 * @details The Geometry is traversed to retrieve all Z-planes and implicit
 *          z-boundaries, such as lattice boundaries. The levels of all these
 *          z-boundaries are rounded and added to a set containing no
 *          duplicates, creating a mesh.
 * @reutrn a vector of z-coords
 */
std::vector<FP_PRECISION> Geometry::getUniqueZHeights() {

  /* Get the bounds of the geometry */
  double min_z = getMinZ();
  double max_z = getMaxZ();

  /* Initialize set for axial mesh */
  std::set<double> unique_mesh;

  /* Initialize vector of unvisited universes and add the root universe */
  std::vector<Universe*> universes;
  universes.push_back(_root_universe);

  /* Initialize vector of offsets */
  std::vector<double> offsets;
  offsets.push_back(0.0);

  /* Cycle through known universes */
  while (!universes.empty()) {

    /* Get the last universe and explore it */
    Universe* curr_universe = universes.back();
    universes.pop_back();

    /* Get the z-offset of the universe */
    double z_offset = offsets.back();
    offsets.pop_back();

    /* Store a vector of the z_heights before rounding */
    std::vector<double> z_heights;

    /* Check if universe is actually a lattice */
    universeType type = curr_universe->getType();
    if (type == LATTICE) {

      /* Get lattice dimensions */
      Lattice* lattice = static_cast<Lattice*>(curr_universe);
      int nx = lattice->getNumX();
      int ny = lattice->getNumY();
      int nz = lattice->getNumZ();

      /* Get offset of the lattice */
      z_offset += lattice->getOffset()->getZ();

      /* Calculate z-intersections */
      double width = lattice->getWidthZ();
      double offset = z_offset - nz * width / 2;
      for (int k=0; k<nz+1; k++) {
        double z_height = k * width + offset;
        z_heights.push_back(z_height);
      }

      /* Add universes to unvisited universes vector */
      for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
          for (int k=0; k<nz; k++) {
            Universe* new_universe = lattice->getUniverse(i, j, k);
            universes.push_back(new_universe);
            offsets.push_back(z_offset);
          }
        }
      }
    }

    /* Otherwise check if universe is simple, contains cells */
    else if (type == SIMPLE) {

      /* Get all cells in the universe */
      std::map<int, Cell*> cells = curr_universe->getCells();

      /* Cycle through all cells */
      std::map<int, Cell*>::iterator cell_iter;
      for (cell_iter = cells.begin(); cell_iter != cells.end(); ++cell_iter) {

        /* Get surfaces bounding the cell */
        std::map<int, surface_halfspace*> surfaces =
          cell_iter->second->getSurfaces();

        /* Cycle through all surfaces and add them to the set */
        std::map<int, surface_halfspace*>::iterator surf_iter;
        for (surf_iter = surfaces.begin(); surf_iter != surfaces.end();
            ++surf_iter) {

          /* Extract surface type */
          Surface* surface = surf_iter->second->_surface;
          surfaceType surf_type = surface->getSurfaceType();

          /* Treat surface types */
          if (surf_type == PLANE) {

            /* Extract plane paramters */
            Plane* plane = static_cast<Plane*>(surface);
            double A = plane->getA();
            double B = plane->getB();
            double C = plane->getC();
            double D = plane->getD();

            /* Check if there is a z-component */
            if (C != 0) {

              /* Check if plane has a continuous varying slope */
              if (A != 0 || B != 0)
                log_printf(ERROR, "Continuous axial variation found in the "
                          "Geometry during axial on-the-fly ray tracing. "
                          "Axial on-the-fly ray tracing only supports "
                          "geometries that are capable of having an axially "
                          "extruded representation");

              /* Otherwise, surface is a z-plane */
              else
                z_heights.push_back(-D/C + z_offset);
            }
          }

          /* Treat explicit z-planes */
          else if (surf_type == ZPLANE) {
            ZPlane* zplane = static_cast<ZPlane*>(surface);
            z_heights.push_back(zplane->getZ() + z_offset);
          }
        }

        /* Add min and max z-height to the cell */
        double z_limits[2];
        z_limits[0] = cell_iter->second->getMinZ();
        z_limits[1] = cell_iter->second->getMaxZ();
        for (int i=0; i < 2; i++) {
          if (std::abs(z_limits[i]) != std::numeric_limits<double>::infinity())
            z_heights.push_back(z_limits[i] + z_offset);
        }

        /* See if cell is filled with universes or lattices */
        cellType cell_type = cell_iter->second->getType();
        if (cell_type == FILL) {
          Universe* new_universe = cell_iter->second->getFillUniverse();
          universes.push_back(new_universe);
          offsets.push_back(z_offset);
        }
      }
    }

    /* Add rounded z-heights to the set of unique z-heights */
    for (int i=0; i < z_heights.size(); i++) {

      /* Round z-height */
      z_heights[i] = floor(z_heights[i] / ON_SURFACE_THRESH)
        * ON_SURFACE_THRESH;

      /* Add the rounded z-height to the set */
      if (z_heights[i] >= min_z && z_heights[i] <= max_z)
        unique_mesh.insert(z_heights[i]);
    }
  }

  /* Add CMFD levels */
  if (_cmfd != NULL) {

    /* Cycle through CMFD mesh not included by boundaries */
    double cmfd_num_z = _cmfd->getNumZ();
    double width = (max_z - min_z) / cmfd_num_z;
    for (int i=1; i < cmfd_num_z; i++) {

      /* Calculate z-height */
      double z_height = min_z + width*i;

      /* Round z-height */
      int place = 8;
      z_height = floor(z_height * pow(10, place)) * pow(10, -place);

      /* Add height to set */
      unique_mesh.insert(z_height);
    }
  }

  /* Get a vector of the unique z-heights in the Geometry */
  std::vector<double> unique_heights;
  std::set<double>::iterator iter;
  for (iter = unique_mesh.begin(); iter != unique_mesh.end(); ++iter)
    unique_heights.push_back(static_cast<FP_PRECISION>(*iter));

  return unique_heights;
}


/**
 * @brief Returns a vector of z-coords defining potential unique radial planes
 *        in the Geometry
 * @details The Geometry is traversed to retrieve all Z-planes and implicit
 *          z-boundaries, such as lattice boundaries. The mid points of this
 *          mesh are then used to construcut a vector of all potential unique
 *          radial planes and returned to the user.
 * @reutrn a vector of z-coords
 */
std::vector<FP_PRECISION> Geometry::getUniqueZPlanes() {

  /* Get a vector of all unique z-heights in the Geometry */
  std::vector<FP_PRECISION> unique_heights = getUniqueZHeights();

  /* Use the midpoints to construct all possible unique radial planes */
  std::vector<FP_PRECISION> unique_z_planes;
  for (int i=1; i < unique_heights.size(); i++) {
    FP_PRECISION mid = (unique_heights[i-1] + unique_heights[i]) / 2;
    unique_z_planes.push_back(mid);
  }

  return unique_z_planes;
}
