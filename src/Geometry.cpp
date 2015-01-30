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

  _num_FSRs = 0;

  _max_seg_length = 0;
  _min_seg_length = std::numeric_limits<double>::infinity();

  /* Initialize CMFD object to NULL */
  _cmfd = NULL;

  /* initialize _num_FSRs lock */
  _num_FSRs_lock = new omp_lock_t;
  omp_init_lock(_num_FSRs_lock);
}


/**
 * @brief Destructor clears FSR to Cells and Materials maps.
 */
Geometry::~Geometry() {

  /* Free FSR  maps if they were initialized */
  if (_num_FSRs != 0) {
    _FSR_keys_map.clear();
    _FSRs_to_keys.clear();
    _FSRs_to_material_IDs.clear();
  }
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
        material = static_cast<CellBasic*>(cell)->getMaterial();
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
Cmfd* Geometry::getCmfd(){
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
void Geometry::setCmfd(Cmfd* cmfd){
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
CellBasic* Geometry::findCellContainingCoords(LocalCoords* coords) {

  Universe* univ = coords->getUniverse();
  Cell* cell;

  if (universe_id == 0){
    if (!withinBounds(coords))
      return NULL;
  }

  if (univ->getType() == SIMPLE)
    cell = univ->findCell(coords);
  else
    cell = static_cast<Lattice*>(univ)->findCell(coords);

  return static_cast<CellBasic*>(cell);
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
        Lattice* lattice = coords->getLattice();
        dist = lattice->minSurfaceDist(coords->getPoint(), angle);
      }
      /* If we reach a LocalCoord in a Universe, find the distance to the
       * nearest cell surface */
      else{
        Universe* universe = coords->getUniverse();
        dist = universe->minSurfaceDist(coords->getPoint(), angle);
      }

      /* Recheck min distance */
      min_dist = std::min(dist, min_dist);

      /* Ascend one level */
      if (coords->getUniverse() == _root_universe)
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
      _FSRs_to_material_IDs.push_back(cell->getMaterial()->getId());

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
    log_printf(ERROR, "Could not find characteristic point in FSR: %i. "
               "Backtrace:%s", fsr_id, e.what());
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
      curr_level_key << curr->getLattice()->getId();
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

  std::map<int, Cell*>::iterator iter1;
  std::map<int, Cell*> cells;

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
 * @param max_optical_length the maximum optical length a segment is allowed to
 *          have
 */
void Geometry::segmentize(Track* track, FP_PRECISION max_optical_length) {

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
  int num_groups;

  /* Use a LocalCoords for the start and end of each segment */
  LocalCoords segment_start(x0, y0);
  LocalCoords segment_end(x0, y0);
  segment_start.setUniverse(_root_universe);
  segment_end.setUniverse(_root_universe);

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
    segment_material = static_cast<CellBasic*>(prev)->getMaterial();
    sigma_t = segment_material->getSigmaT();

    /* Find the ID of the FSR that contains the segment */
    fsr_id = findFSRId(&segment_start);

    /* Compute the number of Track segments to cut this segment into to ensure
     * that it's length is small enough for the exponential table */
    min_num_segments = 1;
    num_groups = segment_material->getNumEnergyGroups();
    for (int g=0; g < num_groups; g++) {
      num_segments = ceil(segment_length * sigma_t[g] / max_optical_length);
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
void Geometry::initializeCmfd(){

  /* Get information about geometry and CMFD mesh */
  int num_x = _cmfd->getNumX();
  int num_y = _cmfd->getNumY();
  double height = getHeight();
  double width = getWidth();
  double cell_width = width / num_x;
  double cell_height = height / num_y;

  /* Create CMFD lattice and set properties */
  Lattice* lattice = new Lattice();
  lattice->setWidth(cell_width, cell_height);
  lattice->setNumX(num_x);
  lattice->setNumY(num_y);
  lattice->setOffset(getMinX() + getWidth()/2.0, 
                     getMinY() + getHeight()/2.0);
  _cmfd->setLattice(lattice);


  /* Set CMFD mesh boundary conditions */
  _cmfd->setBoundary(0, getMinXBoundaryType());
  _cmfd->setBoundary(1, getMinYBoundaryType());
  _cmfd->setBoundary(2, getMaxXBoundaryType());
  _cmfd->setBoundary(3, getMaxYBoundaryType());

  /* Set CMFD mesh dimensions and number of groups */
  _cmfd->setWidth(width);
  _cmfd->setHeight(height);
  _cmfd->setNumMOCGroups(getNumEnergyGroups());

  /* If user did not set CMFD group structure, create CMFD group
  * structure that is the same as the MOC group structure */
  if (_cmfd->getNumCmfdGroups() == 0)
    _cmfd->setGroupStructure(NULL, getNumEnergyGroups()+1);

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
 * @details The _FSR_keys_map stores a hash of a std::string representing
 *          the Lattice/Cell/Universe hierarchy for a unique region
 *          and the associated FSR data. fsr_data is a struct that contains
 *          a unique FSR id and a Point located in the highest level Universe
 *          that is contained in the FSR. This method is used when the tracks 
 *          are read from file to avoid unnecessary segmentation.  
 * @param FSR_keys_map map of FSR keys to FSR data
 */
void Geometry::setFSRKeysMap(std::unordered_map<std::size_t, fsr_data> 
                             FSR_keys_map){
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


/**
 * @brief Determins whether a point is within the bounding box of the geometry.
 * @param coords a populated LocalCoords linked list
 * @return boolean indicating whether the coords is within the geometry
 */
bool Geometry::withinBounds(LocalCoords* coords){

  double x = coords->getX();
  double y = coords->getY();

  if (x < getMinX() || x > getMaxX() || y < getMinY() || y > getMaxY())
    return false;
  else
    return true;
}
