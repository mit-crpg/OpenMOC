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

    _FSR_keys_map.clear();
    _FSRs_to_keys.clear();
    _FSRs_to_material_IDs.clear();
  }

  /* Remove all Materials in the Geometry */
  std::map<int, Material*> materials = getAllMaterials();
  std::map<int, Cell*> cells = getAllCells();
  std::map<int, Universe*> universes = getAllUniverses();
  std::map<int, Material*>::iterator m_iter;
  std::map<int, Cell*>::iterator c_iter;
  std::map<int, Universe*>::iterator u_iter;

  /* Remove all Materials in the Geometry */
  for (m_iter = materials.begin(); m_iter != materials.end(); ++m_iter)
    delete m_iter->second;

  /* Remove all Cells in the Geometry */
  for (c_iter = cells.begin(); c_iter != cells.end(); ++c_iter)
    delete c_iter->second;

  /* Remove all Universes in the Geometry */
  for (u_iter = universes.begin(); u_iter != universes.end(); ++u_iter)
    delete u_iter->second;
}


/**
 * @brief Returns the total height (y extent) of the Geometry in cm.
 * @return the total height of the Geometry (cm)
 */
double Geometry::getHeight() {
  return (getMaxY() - getMinY());
}


/**
 * @brief Returns the total width (x extent) of the Geometry in cm.
 * @return the total width of the Geometry (cm)
 */
double Geometry::getWidth() {
  return (getMaxX() - getMinX());
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
    std::map<int, Cell*> all_cells = getAllCells();
    num_cells = all_cells.size();
  }

  return num_cells;
}


/**
 * @brief Return a std::map container of Universe IDs (keys) with Unierses
 *        pointers (values).
 * @return a std::map of Universes indexed by Universe ID in the geometry
 */
std::map<int, Universe*> Geometry::getAllUniverses() {

  std::map<int, Universe*> all_universes;

  if (_root_universe != NULL)
    all_universes = _root_universe->getAllUniverses();

  return all_universes;
}


/**
 * @brief Return a std::map container of Cell IDs (keys) with Cells
 *        pointers (values).
 * @return a std::map of Cells indexed by Cell ID in the geometry
 */
std::map<int, Cell*> Geometry::getAllCells() {

  std::map<int, Cell*> all_cells;

  if (_root_universe != NULL)
    all_cells = _root_universe->getAllCells();

  return all_cells;
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
    std::map<int, Cell*> all_cells = getAllCells();
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
 * @brief Return a std::map container of Material IDs (keys) with Materials
 *        pointers (values).
 * @return a std::map of Materials indexed by Material ID in the geometry
 */
std::map<int, Material*> Geometry::getAllMaterials() {

  std::map<int, Material*> all_materials;
  Cell* cell;
  Material* material;

  if (_root_universe != NULL) {
    std::map<int, Cell*> all_cells = getAllMaterialCells();
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
Cell* Geometry::findFirstCell(LocalCoords* coords, double angle) {
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
Cell* Geometry::findNextCell(LocalCoords* coords, double angle) {

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
        dist = lattice->minSurfaceDist(coords->getPoint(), angle);
      }
      /* If we reach a LocalCoord in a Universe, find the distance to the
       * nearest cell surface */
      else {
        Cell* cell = coords->getCell();
        dist = cell->minSurfaceDist(coords->getPoint(), angle);
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

  int fsr_id;
  LocalCoords* curr = coords;
  curr = coords->getLowestLevel();

  /* Generate unique FSR key */
  std::string fsr_key = getFSRKey(coords);

  /* If FSR has not been encountered, update FSR maps and vectors */
  if (!_FSR_keys_map.contains(fsr_key)) {

    /* Try to get a clean copy of the fsr_id, adding the FSR data 
       if necessary where -1 indicates the key was already added */
    fsr_id = _FSR_keys_map.insert_and_get_count(fsr_key, NULL);
    if (fsr_id == -1)
    {
      fsr_data volatile* fsr;
      do {
        fsr = _FSR_keys_map.at(fsr_key);
      } while (fsr == NULL);
      fsr_id = fsr->_fsr_id;
    }
    else {
      
      /* Add FSR information to FSR key map and FSR_to vectors */
      fsr_data* fsr = new fsr_data;
      fsr->_fsr_id = fsr_id;
      _FSR_keys_map.at(fsr_key) = fsr;
      Point* point = new Point();
      point->setCoords(coords->getHighestLevel()->getX(), 
                       coords->getHighestLevel()->getY());
      
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
      fsr = _FSR_keys_map.at(fsr_key);
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

  try{
    fsr_key = getFSRKey(coords);
    fsr_id = _FSR_keys_map.at(fsr_key)->_fsr_id;
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
 * @brief Subidivides all Cells in the Geometry into rings and angular sectors.
 * @details This method is called by the Geometry::initializeFlatSourceRegions()
 *          method but may also be called by the user in Python if needed:
 *
 * @code
 *          geometry.subdivideCells()
 * @endcode
 */
void Geometry::subdivideCells() {

  std::map<int, Universe*> all_universes = _root_universe->getAllUniverses();
  std::map<int, Universe*>::iterator iter;

  /* Loop over all Universe in the Geometry and instruct each to inform
   * their Cells to subdivide into rings and sectors as specified by
   * the user during Cell instantiation */
  for (iter = all_universes.begin(); iter != all_universes.end(); ++iter)
    (*iter).second->subdivideCells();
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
  double delta_x, delta_y;

  /* Length of each segment */
  FP_PRECISION length;
  Material* material;
  int fsr_id;
  int num_segments;

  /* Use a LocalCoords for the start and end of each segment */
  LocalCoords start(x0, y0);
  LocalCoords end(x0, y0);
  start.setUniverse(_root_universe);
  end.setUniverse(_root_universe);

  /* Find the Cell containing the Track starting Point */
  Cell* curr = findFirstCell(&end, phi);
  Cell* prev;

  /* If starting Point was outside the bounds of the Geometry */
  if (curr == NULL)
    log_printf(ERROR, "Could not find a material-filled Cell containing the "
               "start Point of this Track: %s", track->toString().c_str());

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

      new_segment->_cmfd_surface_fwd = _cmfd->findCmfdSurface(cmfd_cell, &end);
      new_segment->_cmfd_surface_bwd =
        _cmfd->findCmfdSurface(cmfd_cell, &start);
      new_segment->_cmfd_corner_fwd = _cmfd->findCmfdCorner(cmfd_cell, &end);
      new_segment->_cmfd_corner_bwd = _cmfd->findCmfdCorner(cmfd_cell, &start);

      /* Re-nudge segments from surface */
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
  std::string *key_list = _FSR_keys_map.keys();
  fsr_data **value_list = _FSR_keys_map.values();

  /* allocate vectors */
  int num_FSRs = _FSR_keys_map.size();
  _FSRs_to_keys = std::vector<std::string>(num_FSRs);
  _FSRs_to_material_IDs = std::vector<int>(num_FSRs);

  /* fill vectors key and material ID information */
  #pragma omp parallel for
  for (int i=0; i < num_FSRs; i++)
  {
    std::string key = key_list[i];
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

  std::map<int, Cell*> all_cells = getAllCells();
  std::map<int, Universe*> all_universes = getAllUniverses();

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
  _cmfd->setBoundary(SURFACE_X_MAX, getMaxXBoundaryType());
  _cmfd->setBoundary(SURFACE_Y_MAX, getMaxYBoundaryType());

  /* Set CMFD mesh dimensions and number of groups */
  _cmfd->setWidth(getWidth());
  _cmfd->setHeight(getHeight());
  _cmfd->setNumMOCGroups(getNumEnergyGroups());

  /* If user did not set CMFD group structure, create CMFD group
  * structure that is the same as the MOC group structure */
  if (_cmfd->getNumCmfdGroups() == 0)
    _cmfd->setGroupStructure(NULL, getNumEnergyGroups()+1);

  /* Intialize CMFD Maps */
  _cmfd->initializeCellMap();
  _cmfd->initializeGroupMap();

  /* Initialize the CMFD lattice */
  Point offset;
  double offset_x = getMinX() + getWidth()/2.0;
  double offset_y = getMinY() + getHeight()/2.0;
  offset.setX(offset_x);
  offset.setY(offset_y);
  _cmfd->initializeLattice(&offset);
}


/**
 * @brief Returns a pointer to the map that maps FSR keys to FSR IDs
 * @return pointer to _FSR_keys_map map of FSR keys to FSR IDs
 */
ParallelHashMap<std::string, fsr_data*>* Geometry::getFSRKeysMap() {
  return &_FSR_keys_map;
}


/**
 * @brief Returns the vector that maps FSR IDs to FSR key hashes
 * @return _FSR_keys_map map of FSR keys to FSR IDs
 */
std::vector<std::string>* Geometry::getFSRsToKeys() {
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
void Geometry::setFSRKeysMap(ParallelHashMap<std::string, fsr_data*>* 
                             FSR_keys_map) {
  _FSR_keys_map = *FSR_keys_map;
}

/**
 * @brief Sets the _FSRs_to_keys vector
 * @param FSRs_to_keys vector of FSR key hashes indexed by FSR IDs
 */
void Geometry::setFSRsToKeys(std::vector<std::string>* FSRs_to_keys) {
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

  if (x < getMinX() || x > getMaxX() || y < getMinY() || y > getMaxY())
    return false;
  else
    return true;
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


Cell* Geometry::findCellContainingFSR(int fsr_id) {

  Point* point = _FSR_keys_map.at(_FSRs_to_keys[fsr_id])->_point;
  LocalCoords* coords = new LocalCoords(point->getX(), point->getY());
  coords->setUniverse(_root_universe);
  Cell* cell = findCellContainingCoords(coords);

  delete coords;

  return cell;
}
