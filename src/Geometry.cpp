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
Geometry::Geometry(Mesh* mesh) {

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

  if (mesh == NULL)
    _mesh = new Mesh();
  else
    _mesh = mesh;
}


/**
 * @brief Destructor clears FSR to Cells and Materials maps.
 */
Geometry::~Geometry() {

  /* Free FSR-to-Cells and Materials offset maps if they were initialized */
  if (_num_FSRs != 0) {
    delete [] _FSRs_to_cells;
    delete [] _FSRs_to_material_UIDs;
    delete [] _FSRs_to_material_IDs;
  }
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
  std::map<int, Material*> all_materials = getAllMaterials();
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
 * @brief Return an array indexed by flat source region IDs which contain
 *        the corresponding Cell IDs.
 * @return an integer array of FSR-to-Cell IDs indexed by FSR ID
 */
int* Geometry::getFSRtoCellMap() {
  if (_num_FSRs == 0)
    log_printf(ERROR, "Unable to return the FSR-to-Cell map array since "
               "the Geometry has not initialized FSRs.");

  return _FSRs_to_cells;
}


/**
 * @brief Return an array indexed by flat source region IDs which contain
 *        the corresponding Material IDs.
 * @return an integer array of FSR-to-Material UIDs indexed by FSR ID
 */
int* Geometry::getFSRtoMaterialMap() {
  if (_num_FSRs == 0)
    log_printf(ERROR, "Unable to return the FSR-to-Material map array since "
               "the Geometry has not initialized FSRs.");

  return _FSRs_to_material_UIDs;
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
 * @brief Returns the Universe at the root node in the CSG tree.
 * @return the root Universe
 */
Universe* Geometry::getRootUniverse() {
  return _root_universe;
}


/**
 * @brief Returns a pointer to the CMFD Mesh object.
 * @return A pointer to the CMFD Mesh object
 */
Mesh* Geometry::getMesh(){
    return _mesh;
}


/**
 * @brief Sets the root Universe for the CSG tree.
 * @param root_universe the root Universe of the CSG tree.
 */
void Geometry::setRootUniverse(Universe* root_universe) {
  _root_universe = root_universe;
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

  if (univ->getType() == SIMPLE)
    return univ->findCell(coords);
  else
    return static_cast<Lattice*>(univ)->findCell(coords);
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
 * @brief Find the Cell for a flat source region ID.
 * @details  This metthod calls the recursive Geometry::findCell(...)
 *           method with a pointer to the Universe with ID=0 which at the
 *           uppermost level in the nested Universe hierarchy.
 * @param fsr_id a FSr id
 * @return a pointer to the Cell that this FSR is in
 */
CellBasic* Geometry::findCellContainingFSR(int fsr_id) {
  return static_cast<CellBasic*>(findCell(_root_universe, fsr_id));
}


/**
 * @brief Find the Cell for an FSR ID at a certain level in the nested
 *        Universe hierarchy.
 * @details This is a recursive method which is intended to be called
 *          with the Universe with ID=0 which is at the uppermost level of the
 *          nested Universe hierarchy. The method will recursively call itself
 *          for each level in the hierarchy until it reaches the Cell which
 *          corresponds to the requested FSR ID.
 * @param univ a universe pointer for one level in the nested Universe hierarchy
 * @param fsr_id a FSr ID
 * @return a pointer to the Cell that this FSR is in
 */
Cell* Geometry::findCell(Universe* univ, int fsr_id) {

  Cell* cell = NULL;

  /* Check if the FSR ID is out of bounds */
  if (fsr_id < -1 || fsr_id > _num_FSRs)
    log_printf(ERROR, "Tried to find the Cell for FSR with ID = %d which "
               "does not exist", fsr_id);


  /* If the Universe is a SIMPLE type, then find the Cell the smallest FSR
   * offset map entry that is less than or equal to fsr_id */
  if (univ->getType() == SIMPLE) {
    std::map<int, Cell*>::iterator iter;
    std::map<int, Cell*> cells = univ->getCells();
    Cell* cell_min = NULL;
    int max_id = 0;
    int min_id = INT_MAX;
    int fsr_map_id;

    /* Loop over this Universe's Cells */
    for (iter = cells.begin(); iter != cells.end(); ++iter) {
      fsr_map_id = univ->getFSR(iter->first);
      if (fsr_map_id <= fsr_id && fsr_map_id >= max_id) {
        max_id = fsr_map_id;
        cell = iter->second;
      }
      if (fsr_map_id < min_id) {
        min_id = fsr_map_id;
        cell_min = iter->second;
      }
    }

    /* If the max_id is greater than the fsr_id, there has either been
     * an error or we are at Universe 0 and need to go down one level */
    if (max_id > fsr_id) {
      if (cell_min->getType() == MATERIAL)
        log_printf(ERROR, "Could not find Cell for FSR ID = %d: "
                   "max_id (%d) > fsr_id (%d)", fsr_id, max_id, fsr_id);
      else {
        CellFill* cellfill = static_cast<CellFill*>(cell_min);
        return findCell(cellfill->getFill(), fsr_id);
      }
    }
    /* Otherwise, decrement the fsr_id and make recursive call to next
     * Universe unless an error condition is met */
    else {
      fsr_id -= max_id;
      if (fsr_id == 0 && cell_min->getType() == MATERIAL)
        return cell;
      else if (fsr_id != 0 && cell_min->getType() == MATERIAL)
        log_printf(ERROR, "Could not find Cell for FSR ID = %d: "
                   "fsr_id = %d and Cell type = MATERIAL", fsr_id, fsr_id);
      else {
        CellFill* cellfill = static_cast<CellFill*>(cell_min);
        return findCell(cellfill->getFill(), fsr_id);
      }
    }
  }

  /* If the Universe is a Lattice then we find the Lattice cell with the
   * smallest FSR offset map entry that is less than or equal to fsr_id */
  else {
    Lattice* lat = static_cast<Lattice*>(univ);
    Universe* next_univ = NULL;
    int num_y = lat->getNumY();
    int num_x = lat->getNumX();
    int max_id = 0;
    int fsr_map_id;

    /* Loop over all Lattice cells */
    for (int i = 0; i < num_y; i++) {
      for (int j = 0; j < num_x; j++) {

        fsr_map_id = lat->getFSR(j, i);

        if (fsr_map_id <= fsr_id && fsr_map_id >= max_id) {
          max_id = fsr_map_id;
          next_univ = lat->getUniverse(j, i);
        }
      }
    }

    /* If the max_id is out of bounds, then query failed */
    if (max_id > fsr_id || next_univ == NULL)
      log_printf(ERROR, "No Lattice cell found for FSR ID = %d, max_id = "
                 "%d", fsr_id, max_id);

    /* Otherwise update fsr_id and make recursive call to next level */
    fsr_id -= max_id;
    return findCell(next_univ, fsr_id);
  }

  return cell;
}


/**
 * @brief Finds the next Cell for a LocalCoords object along a trajectory
 *        defined by some angle (in radians from 0 to Pi).
 * @details The method will update the LocalCoords passed in as an argument
 *          to be the one at the boundary of the next Cell crossed along the
 *          given trajectory. It will do this by recursively building a linked
 *          list of LocalCoords from the LocalCoords passed in as an argument
 *          down to the Cell in the lowest level of the Universe hiearchy. In
 *          the process, the method will set the coordinates for each
 *          LocalCoords in the linked list for the Lattice or Universe that it
 *          is in at that nested Universe level. If the LocalCoords is outside
 *          the bounds of the Geometry or on the boundaries this method will
 *          will return NULL; otherwise it will return a pointer to the Cell
 *          that the LocalCoords will reach next along its trajectory.
 * @param coords pointer to a LocalCoords object
 * @param angle the angle of the trajectory
 * @return a pointer to a Cell if found, NULL if no Cell found
 */
Cell* Geometry::findNextCell(LocalCoords* coords, double angle) {

  Cell* cell = NULL;
  double dist;

  /* Find the current Cell */
  cell = findCellContainingCoords(coords);

  /* If the current coords is not in any Cell, return NULL */
  if (cell == NULL)
    return NULL;

  /* If the current coords is inside a Cell, look for next Cell */
  else {

    /* Check the min distance to the next Surface in the current Cell */
    Point surf_intersection;
    LocalCoords* lowest_level = coords->getLowestLevel();
    dist = cell->minSurfaceDist(lowest_level->getPoint(), angle,
                                &surf_intersection);

    /* If the distance returned is not INFINITY, the trajectory will
     * intersect a Surface in the Cell */
    if (dist != std::numeric_limits<double>::infinity()) {
      LocalCoords test(0, 0);

      /* Move LocalCoords just to the next Surface in the Cell plus an
       * additional small bit into the next Cell */
      double delta_x = cos(angle) * TINY_MOVE;
      double delta_y = sin(angle) * TINY_MOVE;

      /* Copy coords to the test coords before moving it by delta and
       * finding the new Cell it is in - do this for testing purposes
       * in case the new Cell found is NULL or is in a new Lattice cell*/
      coords->copyCoords(&test);
      coords->updateMostLocal(&surf_intersection);
      coords->adjustCoords(delta_x, delta_y);

      /* Find new Cell for the perturbed coords */
      cell = findCellContainingCoords(coords);

     /* Check if the next Cell found is in the same lattice cell
      * as the previous cell */
      LocalCoords* test_curr = test.getLowestLevel();
      LocalCoords* coords_curr = coords->getLowestLevel();

      while (test_curr != NULL && test_curr->getUniverse() != NULL &&
             coords_curr != NULL && coords_curr->getUniverse() != NULL){

        if (coords_curr->getType() == LAT && test_curr->getType() == LAT) {

          if (coords_curr->getLatticeX() != test_curr->getLatticeX() ||
              coords_curr->getLatticeY() != test_curr->getLatticeY()) {
            dist = std::numeric_limits<double>::infinity();
            break;
          }
        }

        /* Update LocalCoords linked list to the next level up in the nested
         * Universe hierarchy */
        test_curr = test_curr->getPrev();
        coords_curr = coords_curr->getPrev();
      }

      /* Check if Cell is NULL - this means that intersection point
       * is outside the bounds of the Geometry and the old coords should
       * should be restored so that we can look for the next Lattice cell */
      if (cell == NULL)
        dist = std::numeric_limits<double>::infinity();

      /* If the distance is not INFINITY then the new Cell found is the
       * one to return */
      if (dist != std::numeric_limits<double>::infinity()) {
        test.prune();
        return cell;
      }

      /* If the distance is INFINITY then the new Cell found is not
       * the one to return and we should move to a new Lattice cell */
      else
        test.copyCoords(coords);

      test.prune();
    }

    /* If the distance returned is INFINITY, the trajectory will not
     * intersect a Surface in the Cell. We thus need to readjust to
     * the LocalCoord to the base Universe and check whether we need
     * to move to a new Lattice cell */
    if (dist == std::numeric_limits<double>::infinity()) {

      /* Get the lowest level LocalCoords in the linked list */
      LocalCoords* curr = coords->getLowestLevel();

      /* Retrace linked list from lowest level */
      while (curr != NULL && curr->getUniverse() != NULL) {
        curr = curr->getPrev();

        /* If we reach a LocalCoord in a lattice, delete all lower
         * level LocalCoords in linked list and break loop */
        if (curr->getType() == LAT) {
          curr->prune();
          curr = NULL;
        }
      }

      /* Get the lowest level Universe in linked list */
      curr = coords->getLowestLevel();

      /* Retrace through the Lattices in the LocalCoord and check for
       * Lattice cell crossings in each one. If we never find a crossing
       * and reach universe 0 then return NULL since this means we have
       * reached the edge of the Geometry */
      while (curr->getUniverse() != NULL) {

        /* If the lowest level LocalCoords is inside a Lattice, find
         * the next Lattice cell */
        if (curr->getType() == LAT) {

          Lattice* lattice = curr->getLattice();

          cell = lattice->findNextLatticeCell(curr, angle);

          /* If Cell returned is NULL, the LocalCoords are outside of current
          * Lattice, so move to a higher level Lattice if there is one */
          if (cell == NULL) {

            /* Delete current Lattice */
            curr->getPrev()->prune();

            /* Get the lowest level LocalCoords in linked list */
            curr = coords->getLowestLevel();

            /* Retrace linked list from lowest level */
            while (curr != NULL && curr->getUniverse() != 0) {

              curr = curr->getPrev();

              /* If we reach a LocalCoord in a Lattice, delete all lower level
               * LocalCoords in linked list and break loop. */
              if (curr->getType() == LAT) {
                curr->prune();
                curr = NULL;
              }
            }

            /* Get the lowest level Universe in linked list */
            curr = coords->getLowestLevel();
          }

          /* If lowest level Universe is not a Lattice, return current Cell */
          else
            return cell;
        }
      }
    }
  }

  /* If no Cell was found, return NULL */
  return NULL;
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

  /* Traverse linked list defined by the LocalCoords object and compute the FSR
   * ID using the offset maps at each level of the nested Universe hierarchy */
  while (curr != NULL) {

    /* If the current level is a Lattice, add an offset from the Lattice map
     * to the FSR ID */
    if (curr->getType() == LAT) {
      Lattice* lattice = curr->getLattice();
      fsr_id += lattice->getFSR(curr->getLatticeX(), curr->getLatticeY());
    }

    /* If the current level is a Lattice, add an offset from the Lattice map
     * to the FSR ID*/
    else if (curr->getType() == UNIV) {
      Universe* universe = curr->getUniverse();
      fsr_id += universe->getFSR(curr->getCell()->getId());
    }

    /* Get next LocalCoords node in the linked list */
    curr = curr->getNext();
  }

  return fsr_id;
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
  Universe* curr;

  /* Loop over all Universe in the Geometry and instruct each to inform
   * their Cells to subdivide into rings and sectors as specified by
   * the user during Cell instantiation */
  for (iter = all_universes.begin(); iter != all_universes.end(); ++iter) {
    curr = (*iter).second;

    if (curr->getType() == SIMPLE)
      curr->subdivideCells();
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

  /* Subdivide Cells into sectors and rings */
  subdivideCells();

  /* Generate flat source regions offset maps for each Universe and Lattice */
  _num_FSRs = _root_universe->computeFSRMaps();

  log_printf(NORMAL, "Number of flat source regions: %d", _num_FSRs);

  /* Allocate memory for maps between FSR IDs and Cell or Material IDs/UIDs */
  _FSRs_to_cells = new int[_num_FSRs];
  _FSRs_to_material_UIDs = new int[_num_FSRs];
  _FSRs_to_material_IDs = new int[_num_FSRs];

  /* Load maps with Cell and Material IDs/UIDs */
  for (int r=0; r < _num_FSRs; r++) {
    CellBasic* curr = static_cast<CellBasic*>(findCell(_root_universe, r));
    _FSRs_to_cells[r] = curr->getId();
    _FSRs_to_material_UIDs[r] = curr->getMaterial()->getUid();
    _FSRs_to_material_IDs[r] = curr->getMaterial()->getId();
  }

  if (_mesh->getCmfdOn())
    initializeMesh();
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

      /* Get pointer to CMFD Mesh surfaces that the Track segment crosses */
      if (_mesh->getCmfdOn()){

        new_segment->_mesh_surface_fwd =
                _mesh->findMeshSurface(new_segment->_region_id, &segment_end);
        new_segment->_mesh_surface_bwd =
                _mesh->findMeshSurface(new_segment->_region_id, &segment_start);
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
  _root_universe->computeFissionability();
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

  string << "Geometry: width = "
         << getWidth() << ", height = "
         << getHeight() << ", Bounding Box: (("
         << _x_min << ", "
         << _y_min << "), ("
         << _x_max << ", "
         << _y_max << ")";

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
 * @brief This is a recursive method which makes a mesh for solving the
 *        Course Mesh Finite Difference (CMFD) diffusion equations.
 * @details The CMFD Mesh must be a structured Cartesian mesh and defined
 *          as a certain nested Univeres level (counting from the top). This
 *          method defines the mesh dimensions (width and height) and mesh cell
 *          dimensions (width and height). This function adds new Mesh cell
 *          objects to the Mesh and defines the values in each Mesh cell.
 */
void Geometry::initializeMesh(){

  Universe* univ = _root_universe;

  int max_mesh_level = 0;
  int mesh_level = 0;

  /* Find the Mesh depth of the Geometry */
  max_mesh_level = findMeshDepth(univ, max_mesh_level);
  log_printf(DEBUG, "Max mesh depth is: %i level(s)", max_mesh_level);

  /* Set CMFD level to user specified value if possible */
  if (_mesh->getMeshLevel() == -1){
    mesh_level = max_mesh_level;
    _mesh->setMeshLevel(mesh_level);
  }

  else if (_mesh->getMeshLevel() >= 0 &&
           _mesh->getMeshLevel() <= max_mesh_level)
    mesh_level = _mesh->getMeshLevel();
  else{
    log_printf(WARNING, "User input CMFD Mesh level was outside the bounds of "
               "the Mesh level range (%d to %d)", 0, max_mesh_level);
    mesh_level = max_mesh_level;
    _mesh->setMeshLevel(max_mesh_level);
  }

  log_printf(INFO, "CMFD Mesh level: %i, max mesh level: %i",
             mesh_level, max_mesh_level);

  /* Find Cell width and height at Mesh nested Universe level */
  int width = 0;
  int height = 0;
  findMeshHeight(univ, &height, mesh_level);
  univ = _root_universe;
  findMeshWidth(univ, &width, mesh_level);

  /* Set Mesh boundary conditions */
  _mesh->setBoundary(0,getBCLeft());
  _mesh->setBoundary(1,getBCBottom());
  _mesh->setBoundary(2,getBCRight());
  _mesh->setBoundary(3,getBCTop());
  _mesh->setNumFSRs(_num_FSRs);

  /* Set the Mesh cell width and height */
  if (mesh_level > 0){
    _mesh->setCellsY(height);
    _mesh->setCellsX(width);
  }
  else{
    _mesh->setCellsY(1);
    _mesh->setCellsX(1);
  }

  /* Set Mesh dimensions and initialize Mesh variables */
  _mesh->setLengthY(getHeight());
  _mesh->setLengthX(getWidth());
  _mesh->initialize();
  log_printf(NORMAL, "Number of mesh cells: %i", height*width);
  log_printf(DEBUG, "mesh cell width: %i", _mesh->getCellsX());
  log_printf(DEBUG, "mesh cell height: %i", _mesh->getCellsY());

  /* Decide whether CMFD acceleration is needed for MOC acceleration */
  if (_num_FSRs <= 1000 && _mesh->getSolveType() == MOC){
    _mesh->setAcceleration(false);
    log_printf(INFO, "CMFD acceleration was turned off because there are "
               "<= 100 fsrs and CMFD is not needed for small geometries");
  }

  /* Decide whether optically thick correction factor is needed for MOC
   * acceleration */
  if (getHeight()*getWidth() / (height*width) >= 10.0 &&
      _mesh->getSolveType() == MOC){

    _mesh->setOpticallyThick(true);
    log_printf(INFO, "Optically thick correction factor turned on for CMFD "
               "acceleration because the average mesh cell size is >= 10 cm^2");
  }

  /* Make a vector of FSR IDs in each Mesh cell */
  int meshCellNum = 0;

  if (mesh_level > 0)
    defineMesh(_mesh, univ, mesh_level, &meshCellNum, 0, true, 0);

  else{
    _mesh->setCellLengthX(0, getWidth());
    _mesh->setCellLengthY(0, getHeight());

    for (int fsr = 0; fsr < _num_FSRs; fsr++)
      _mesh->getCellFSRs()->at(0).push_back(fsr);
  }

  /* set mesh fsr and cell bounds and initialize materials */
  _mesh->setFSRBounds();
  _mesh->setCellBounds();

  if (_mesh->getSolveType() == DIFFUSION) {
    std::map<int, Material*> all_materials = getAllMaterials();
    _mesh->initializeMaterialsDiffusion(&all_materials, _FSRs_to_material_IDs);
  }

  return;
}


/**
 * @brief This is a recursive method which stores the IDs of all FSRs located
 *        in a Mesh cell object in a std::vector owned by the Mesh cell.
 * @param univ a pointer to the Universe that contains the Mesh cell
 * @param cell_num the Mesh cell indice
 * @param fsr_id an integer pointer set to the first fsr_id in this Universe
 */
void Geometry::findFSRsInCell(Universe* univ, int cell_num, int* fsr_id){

  /* If the Universe is a SIMPLE type Universe */
  if (univ->getType() == SIMPLE) {

    std::map<int, Cell*> cells = univ->getCells();
    Cell* curr;

    /* For each of the Cells inside the Universe, check if it is a
     * MATERIAL or FILL type */
    std::map<int, Cell*>::iterator iter;
    for (iter = cells.begin(); iter != cells.end(); ++iter) {
      curr = iter->second;

      /* If the current Cell is a MATERIAL type cell, store its FSR ID */
      if (curr->getType() == MATERIAL) {

        log_printf(DEBUG, "Storing FSR ID = %i in CMFD Mesh", *fsr_id);
        _mesh->getCellFSRs()->at(cell_num).push_back(*fsr_id);
        log_printf(DEBUG, "cell num %i, fsr list: %i",
                   cell_num, _mesh->getCellFSRs()->at(cell_num).size());
        *fsr_id += 1;
      }

      /* If the current Cell is a FILL type cell recursively call findFSRsInCell() */
      else {
        CellFill* fill_cell = static_cast<CellFill*>(curr);
        Universe* fill = fill_cell->getFill();
        findFSRsInCell(fill, cell_num, fsr_id);
      }
    }
  }

  /* If the Universe is a LATTICE type Universe recursively call 
     findFSRsInCell() */
  else {

    Lattice* lattice = static_cast<Lattice*>(univ);
    Universe* curr;
    int num_x = lattice->getNumX();
    int num_y = lattice->getNumY();
    int baseFSR = *fsr_id;

    /* Loop over all Lattice cells in this Lattice */
    for (int i = num_y-1; i > -1; i--) {
      for (int j = 0; j < num_x; j++) {

        /* Get a pointer to the current Lattice cell */
        curr = lattice->getUniverse(j, i);
        log_printf(DEBUG, "Getting Lattice FSR offset= %i",
                   lattice->getFSR(j,i));
        *fsr_id = baseFSR + lattice->getFSR(j,i);

        /* Find all FSRs in this Lattice */
        findFSRsInCell(curr, cell_num, fsr_id);
      }
    }
  }
}


/**
 * @brief This is a recursive method which defines all the parameters of the
 *        the Mesh cell objects in a Mesh.
 * @details This method takes in the uppermost level Universe in the nested
 *          Universe hierarchy (with ID = 0) and recurses until it reaches the
 *          Universe level of the CMFD mesh. Then, the function loops over all
 *          the Cells in the Lattice and defines the corresponding Mesh cell
 *          object for each Lattice cell.
 * @param mesh the CMFD Mesh object to initialize
 * @param univ a pointer to a the base Universe (ID = 0)
 * @param depth the number of Lattices that must be descended to reach the
 *              CMFD Mesh level.
 * @param meshCellNum a pointer to an integer used to store the index of the
 *                    current Mesh cell object.
 * @param row the current row of the parent Lattice
 * @param base boolean indicating whether the current Lattice is the highest
 *             level Lattice
 * @param fsr_id a pointer to an integer that is set to first fsr_id in this
 *               Universe
 */
void Geometry::defineMesh(Mesh* mesh, Universe* univ, int depth,
                          int* meshCellNum, int row, bool base, int fsr_id){

  /* If the Universe is a SIMPLE type Universe */
  if (univ->getType() == SIMPLE){
    std::map<int, Cell*> cells = univ->getCells();
    Cell* curr;

    /* For each of the Cells inside the Lattice, check if it is
     * MATERIAL or FILL type */
    std::map<int, Cell*>::iterator iter;
    for (iter = cells.begin(); iter != cells.end(); ++iter) {
      curr = iter->second;CellFill* fill_cell = static_cast<CellFill*>(curr);
      Universe* fill = fill_cell->getFill();
      defineMesh(mesh, fill, depth, meshCellNum, row, base, fsr_id);
    }
  }

  /* If the Universe is a LATTICE type Universe */
  else {

    Lattice* lattice = static_cast<Lattice*>(univ);
    Universe* curr;
    int num_x = lattice->getNumX();
    int num_y = lattice->getNumY();
    log_printf(DEBUG, "num_x: %i num_y: %i", num_x, num_y);

    /* If the current LATTICE is the CMFD Mesh Lattice */
    if (depth == 1){

      /* If the current LATTICE is the base Lattice */
      if (base == true){
        for (int i = num_y-1; i > -1; i--) {
          for (int j = 0; j < num_x; j++) {

            curr = lattice->getUniverse(j,i);
            fsr_id = lattice->getFSR(j,i);
            log_printf(DEBUG, "Added FSR ID to counter -> fsr_id: %i", fsr_id);

            /* Store fsr_ids of FSRs in this LATTICE in a Mesh cell object */
            findFSRsInCell(curr, *meshCellNum, &fsr_id);
            mesh->setCellLengthX(*meshCellNum, lattice->getWidthX());
            mesh->setCellLengthY(*meshCellNum, lattice->getWidthY());

            log_printf(DEBUG, "Mesh cell: %i, width: %f, height: %f",
                      *meshCellNum, lattice->getWidthX(), lattice->getWidthY());

            *meshCellNum = *meshCellNum + 1;
          }
        }
      }

      /* If the current LATTICE is not the base Lattice */
      else{

        int baseFSR = fsr_id;

        for (int j = 0; j < num_x; j++) {

          curr = lattice->getUniverse(j,row);
          fsr_id = baseFSR + lattice->getFSR(j,row);

          log_printf(DEBUG, "Set FSr ID to: %i", fsr_id);

          /* Store fsr_ids of the FSRs in this LATTICE in a Mesh cell object */
          findFSRsInCell(curr, *meshCellNum, &fsr_id);
          mesh->setCellLengthX(*meshCellNum, lattice->getWidthX());
          mesh->setCellLengthY(*meshCellNum, lattice->getWidthY());

          log_printf(DEBUG, "Mesh cell num: %i, width: %f, height: %f",
                     *meshCellNum, lattice->getWidthX(), lattice->getWidthY());

          *meshCellNum = *meshCellNum + 1;
        }
      }
    }

    /* If the current LATTICE is not the CMFD Mesh Lattice */
    else {

      base = false;

      for (int i = num_y-1; i > -1; i--) {

        curr = lattice->getUniverse(0,i);
        int nextHeight = nextLatticeHeight(curr);

        log_printf(DEBUG, "next height: %i", nextHeight);

        for (int k = nextHeight-1; k > -1; k--) {
          for (int j = 0; j < num_x; j++) {

            curr = lattice->getUniverse(j,i);
            fsr_id = lattice->getFSR(j,i);

            /* Recursively call defineMesh(...) until LATTICE level of
             * CMFD mesh is reached */
            defineMesh(mesh, curr, depth - 1, meshCellNum, k, base, fsr_id);
          }
        }
      }
    }
  }

  return;
}


/**
 * @brief This is a recursive method that finds the Mesh cell height of the next
 *        lowest LATTICE in a given Universe.
 * @param univ a pointer to the Universe of interest
 * @return the level in the nested Universe hierarchy
 */
int Geometry::nextLatticeHeight(Universe* univ){

  int height = 1;

  /* If the Universe is a SIMPLE type Universe */
  if (univ->getType() == SIMPLE){

    std::map<int, Cell*> cells = univ->getCells();
    Cell* curr;

    std::map<int, Cell*>::iterator iter;
    iter = cells.begin();
    curr = iter->second;

    /* IF the Cell is FILL type recursively call nextLatticeHeight(...) */
    if (curr->getType() == FILL){
      CellFill* fill_cell = static_cast<CellFill*>(curr);
      Universe* fill = fill_cell->getFill();
      height = nextLatticeHeight(fill);
    }
  }

  /* If the Universe is a LATTICE type Universe return the height */
  else {
    Lattice* lattice = static_cast<Lattice*>(univ);
    height = lattice->getNumY();
  }

  return height;
}



/**
 * @brief This is a recursive method that finds the Mesh cell height of the
 *        LATTICE at the CMFD Mesh level.
 * @param univ a pointer to the base Universe (ID = 0)
 * @param height a pointer to the accumulator for the Mesh height
 * @param depth the number of nested Universe levels that must be descended to
 *              reach the CMFD Mesh level.
*/
void Geometry::findMeshHeight(Universe* univ, int* height, int depth){

  /* If the Universe is a SIMPLE type Universe */
  if (univ->getType() == SIMPLE){
    std::map<int, Cell*> cells = univ->getCells();
    Cell* curr;
    std::map<int, Cell*>::iterator iter;
    for (iter = cells.begin(); iter != cells.end(); ++iter) {
      curr = iter->second;

      /* IF the Cell is FILL type, recursively call findMeshHeight(...) */
      if (curr->getType() == FILL){
        CellFill* fill_cell = static_cast<CellFill*>(curr);
        Universe* fill = fill_cell->getFill();
        findMeshHeight(fill, height, depth);
      }
    }
  }

  /* If the Universe is a LATTICE type Universe */
  else {

    /* Check to see if CMFD Mesh Lattice depth has been reached */
    if (depth > 0){

      Lattice* lattice = static_cast<Lattice*>(univ);
      Universe* curr;
      int num_y = lattice->getNumY();

      /* If the current LATTICE is the CMFD Mesh Lattice add its numY to
       * height accumulator */
      if (depth == 1){
        *height = *height + num_y;
      }

      else{

        depth = depth - 1;

        for (int i = num_y-1; i > -1; i--) {
          curr = lattice->getUniverse(0, i);

          /* Find the height of the current Universe */
          findMeshHeight(curr, height, depth);
        }
      }
    }
  }

  return;
}


/**
 * @brief This is a recursive method that finds the Mesh cell width of the
 *        LATTICE at the CMFD Mesh level.
 * @param univ a pointer to the base Universe (ID = 0)
 * @param width a pointer to the accumulator for the Mesh width
 * @param depth the number of nested Universe levels that must be descended to
 *              reach the CMFD Mesh level.
 */
void Geometry::findMeshWidth(Universe* univ, int* width, int depth){

  /* If the Universe is a SIMPLE type Universe */
  if (univ->getType() == SIMPLE){
    std::map<int, Cell*> cells = univ->getCells();
    Cell* curr;
    std::map<int, Cell*>::iterator iter;
    for (iter = cells.begin(); iter != cells.end(); ++iter) {
      curr = iter->second;

      /* IF the Cell is FILL type, recursively call findMeshWidth(...) */
      if (curr->getType() == FILL){
        CellFill* fill_cell = static_cast<CellFill*>(curr);
        Universe* fill = fill_cell->getFill();
        findMeshWidth(fill, width, depth);
     }
    }
  }

  /* If the Universe is a LATTICE type Universe */
  else {

    /* check to see if CMFD Mesh Lattice depth has been reached */
    if (depth > 0){

      Lattice* lattice = static_cast<Lattice*>(univ);
      Universe* curr;
      int num_x = lattice->getNumX();
      int num_y = lattice->getNumY();
      int i = num_y-1;

      /* If the current LATTICE is the CMFD Mesh Lattice add its numX to width
      * accumulator */
      if (depth == 1)
        *width = *width + num_x;

      else{

        depth = depth - 1;

        for (int j = 0; j < num_x; j++) {
          curr = lattice->getUniverse(j, i);

          /* Find the width of the current Universe */
          findMeshWidth(curr, width, depth);
        }
      }
    }
  }

  return;
}


/**
 * @brief This is a recursive method that finds the depth of the Geometry Mesh.
 * @param univ a pointer to a the base Universe (ID = 0)
 * @param mesh_level a pointer to the accumulator for the CMFD Mesh level
 */
int Geometry::findMeshDepth(Universe* univ, int mesh_level){

  /* If the Universe is a SIMPLE type Universe */
  if (univ->getType() == SIMPLE){
    std::map<int, Cell*> cells = univ->getCells();

    Cell* curr;
    std::map<int, Cell*>::iterator iter;

    for (iter = cells.begin(); iter != cells.end(); ++iter) {
      curr = iter->second;

      /* IF the Cell is FILL type recursively call findMeshWidth(...) */
      if (curr->getType() == FILL){
        CellFill* fill_cell = static_cast<CellFill*>(curr);
        Universe* fill = fill_cell->getFill();
        mesh_level = findMeshDepth(fill, mesh_level);
      }
    }
  }

  /* If the Universe is a LATTICE type Universe */
  else {

    Lattice* lattice = static_cast<Lattice*>(univ);
    Universe* curr;
    int num_x = lattice->getNumX();
    int num_y = lattice->getNumY();
    mesh_level++;
    int i = num_y-1;
    int j = num_x-1;
    int min = 0;
    int levels = 0;

    /* Loop through Lattice cells and recursively call  */
    for (int x = 0; x <= j; x++){
      for (int y = 0; y <= i; y++){

        curr = lattice->getUniverse(x, y);

        /* Find the width of the current Universe */
        levels = findMeshDepth(curr, mesh_level);

        if (x == 0 && y == 0)
          min = levels;

        else{
          min = std::min(min, levels);
        }

      }
    }

    mesh_level = min;
  }

  return mesh_level;
}
