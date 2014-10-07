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
 * @brief Default constructor used in rings/sectors subdivision of Cells.
 */
Cell::Cell() {
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

  /* If the user-defined ID is in the prohibited range, return an error */
  else if (id >= 10000)
    log_printf(ERROR, "Unable to set the ID of a cell to %d since cell IDs "
               "greater than or equal to 10000 is probibited by OpenMOC.", id);

  /* Use the user-defined ID */
  else
    _id = id;

  _uid = _n;
  _n++;
  setName(name);

  /* By default, the Cell's fissionability is unknown */
  _fissionable = false;
  _fissionable = false;
}


/**
 * @brief Destructor clears vector of Surface pointers bounding the Cell.
 */
Cell::~Cell() {
  _surfaces.clear();
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
std::map<Surface*, int> Cell::getSurfaces() const {
  return _surfaces;
}


/**
 * @brief Returns true if the Cell contains a fissionable Material and
 *        false otherwise.
 * @return wehther or not the Cell contains a fissionable Material
 */
bool Cell::isFissionable() {
  return _fissionable;
}


/**
 * @brief Sets the name of the Cell
 * @param name the Cell name string
 */
void Cell::setName(const char* name) {
  int length = strlen(name);

  /* Initialize a character array for the Cell's name */
  _name = new char[length+1];

  /* Copy the input character array Cell name to the class attribute name */
  for (int i=0; i <= length; i++)
    _name[i] = name[i];
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

  _surfaces.insert(std::pair<Surface*, int>(surface, halfspace));
}


/**
 * @brief Removes a Surface from this Cell's container of bounding Surfaces.
 * @param surface a pointer to the Surface to remove
 */
void Cell::removeSurface(Surface* surface) {
  if (_surfaces.find(surface) != _surfaces.end())
    _surfaces.erase(surface);
}


/**
 * @brief Determines whether a Point is contained inside a Cell.
 * @details Queries each Surface inside the Cell to determine if the Point
 *          is on the same side of the Surface. This point is only inside
 *          the Cell if it is on the same side of every Surface in the Cell.
 * @param point a pointer to a Point
 */
bool Cell::cellContainsPoint(Point* point) {

  /* Loop over all Surfaces inside the Cell */
  std::map<Surface*, int>::iterator iter;
  for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {

    /* Return false if the Point is not in the correct Surface halfspace */
    if (iter->first->evaluate(point) * iter->second < -ON_SURFACE_THRESH)
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
bool Cell::cellContainsCoords(LocalCoords* coords) {
  return cellContainsPoint(coords->getPoint());
}


/**
 * @brief Computes the minimum distance to a Surface from a Point with a given
 *        trajectory at a certain angle.
 * @details If the trajectory will not intersect any of the Surfaces in the
 *          Cell returns INFINITY.
 * @param point the Point of interest
 * @param angle the angle of the trajectory (in radians from \f$[0,2\pi]\f$)
 * @param min_intersection a pointer to the intersection Point that is found
 */
double Cell::minSurfaceDist(Point* point, double angle,
                            Point* min_intersection) {

  double min_dist = INFINITY;
  double d;
  Point intersection;

  std::map<Surface*, int>::iterator iter;

  /* Loop over all of the Cell's Surfaces */
  for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter) {

    /* Find the minimum distance from this surface to this Point */
    d = iter->first->getMinDistance(point, angle, &intersection);

    /* If the distance to Cell is less than current min distance, update */
    if (d < min_dist) {
      min_dist = d;
      min_intersection->setX(intersection.getX());
      min_intersection->setY(intersection.getY());
    }
  }

  return min_dist;
}


/**
 * @brief Constructor sets the user-specified and unique IDs for this CellBasic.
 * @param id the user-specified optional Cell ID
 * @param name the user-specified optional Cell name
 * @param rings the number of equal volume rings to divide this Cell into
 *        (the default is zero)
 * @param sectors the number of angular sectors to divide this Cell into
 *        (the default is zero)
 */
CellBasic::CellBasic(int id, const char* name, int rings, int sectors):
  Cell(id, name) {

  _cell_type = MATERIAL;
  setNumRings(rings);
  setNumSectors(sectors);
}



/**
 * @brief Return the Material filling the CellBasic.
 * @return the Material's pointer
 */
Material* CellBasic::getMaterial() {
  return _material;
}


/**
 * @brief Return the number of rings in the Cell.
 * @return the number of rings
 */
int CellBasic::getNumRings() {
  return _num_rings;
}


/**
 * @brief Return the number of sectors in the Cell.
 * @return the number of sectors
 */
int CellBasic::getNumSectors() {
  return _num_sectors;
}


/**
 * @brief Return the number of flat source regions in this CellBasic.
 * @details This method is used when the Geometry recursively constructs flat
 *           source regions. By definition, CellBasics are the lowest level
 *           entity in the geometry and thus only have one flat source region
 *           within them, so this method always returns 1.
 * @return the number of FSRs in this Cell
 */
int CellBasic::getNumFSRs() {
  return 1;
}


/**
 * @brief Returns an empty std::map of Cell IDs and Cell pointers.
 * @return empty std::map of Cell IDs and pointers
 */
std::map<int, Cell*> CellBasic::getAllCells() {
  std::map<int, Cell*> cells;
  return cells;
}


/**
 * @brief Returns an empty std::map of nested Universe IDs and pointers
 * @return empty std::map of Universe IDs and pointers
 */
std::map<int, Universe*> CellBasic::getAllUniverses() {
  std::map<int, Universe*> universes;
  return universes;
}


/**
 * @brief Sets the Material filling this Cell.
 * @param material the Material filling this Cell
 */
void CellBasic::setMaterial(Material* material) {
  _material = material;
}


/**
 * @brief Set the Cell's number of rings.
 * @param num_rings the number of rings in this Cell
 */
void CellBasic::setNumRings(int num_rings) {
  if (num_rings < 0)
    log_printf(ERROR, "Unable to give %d rings to Cell %d since this is "
               "a negative number", num_rings, _id);

  _num_rings = num_rings;
}


/**
 * @brief Set the Cell's number of sectors.
 * @param num_sectors the number of sectors in this Cell
 */
void CellBasic::setNumSectors(int num_sectors) {
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
 * @brief Create a duplicate of the CellBasic.
 * @return a pointer to the clone
 */
CellBasic* CellBasic::clone() {

  /* Construct new Cell */
  CellBasic* new_cell = new CellBasic();
  new_cell->setNumRings(_num_rings);
  new_cell->setNumSectors(_num_sectors);

  /* Loop over all of this Cell's Surfaces and add them to the clone */
  std::map<Surface*, int>::iterator iter;

  for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
    new_cell->addSurface(iter->second, iter->first);

  return new_cell;
}


/**
 * @brief Subdivides the Cell into clones for fuel pin angular sectors.
 */
void CellBasic::sectorize() {

  /* If the user didn't request any sectors, don't make any */
  if (_num_sectors == 0)
    return;

  double azim_angle;
  double delta_azim = 2. * M_PI / _num_sectors;
  double A, B;

  /* A container for each of the bouding planes for the sector Cells */
  std::vector<Plane*> planes;

  log_printf(DEBUG, "Sectorizing Cell %d with %d sectors",_id, _num_sectors);

  /* Create each of the bounding planes for the sector Cells */
  for (int i=0; i < _num_sectors; i++) {

    /* Figure out the angle for this plane */
    azim_angle = i * delta_azim;

    /* Instantiate the plane */
    A = cos(azim_angle);
    B = sin(azim_angle);
    Plane* plane = new Plane(A, B, 0.);
    planes.push_back(plane);

    log_printf(DEBUG, "Created sector Plane id = %d, angle = %f, A = %f, "
               "B = %f", i, azim_angle, A, B);
  }

  /* Create sectors using disjoint halfspaces of pairing Planes */
  for (int i=0; i < _num_sectors; i++) {

    /* Create new CellBasic clone for this sector Cell */
    CellBasic* sector = clone();

    sector->setNumSectors(0);
    sector->setNumRings(0);

    log_printf(DEBUG, "Creating a new sector Cell with %d for Cell %d",
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
    _sectors.push_back(sector);
  }

  /* Store all of the sectors in the parent Cell's subcells container */
  _subcells.clear();
  _subcells.insert(_subcells.end(), _sectors.begin(), _sectors.end());
}


/**
 * @brief Subdivides the Cell into clones for fuel pin rings.
 */
void CellBasic::ringify() {

  /* If the user didn't request any rings, don't make any */
  if (_num_rings == 0)
        return;

  int num_circles = 0;
  Circle* circle1 = NULL;
  Circle* circle2 = NULL;
  double radius1 = 0;
  double radius2 = 0;
  double x1 = 0.;
  double y1 = 0.;
  double x2 = 0.;
  double y2 = 0.;
  int halfspace1 = 0;
  int halfspace2 = 0;
  std::vector<Circle*> circles;

  /* See if the Cell contains 1 or 2 CIRCLE Surfaces */
  std::map<Surface*, int>::iterator iter1;
  for (iter1=_surfaces.begin(); iter1 != _surfaces.end(); ++iter1) {

    /* Determine if any of the Surfaces is a Circle */
    if ((*iter1).first->getSurfaceType() == CIRCLE) {
      int halfspace = (*iter1).second;
      Circle* circle = static_cast<Circle*>((*iter1).first);

      /* Outermost bounding Circle */
      if (halfspace == -1) {
        halfspace1 = halfspace;
        circle1 = circle;
        radius1 = circle1->getRadius();
        x1 = circle1->getX0();
        y1 = circle1->getY0();
      }

      /* Innermost bounding circle */
      else if (halfspace == +1) {
        halfspace2 = halfspace;
        circle2 = circle;
        radius2 = circle2->getRadius();
        x2 = circle2->getX0();
        y2 = circle2->getY0();
      }

      num_circles++;
    }
  }

  /* Error checking */
  if (num_circles == 0)
    log_printf(ERROR, "Unable to ringify Cell %d since it does not "
              "contain any CIRCLE type Surface(s)", _id);

  if (num_circles > 2)
    log_printf(NORMAL, "Unable to ringify Cell %d since it "
               "contains more than 2 CIRCLE Surfaces", _id);

  if (x1 != x2 && num_circles == 2)
    log_printf(ERROR, "Unable to ringify Cell %d since it contains "
               "Circle %d centered at x=%f and Circle %d at x=%f. "
               "Both Circles must have the same center.",
               _id, circle1->getId(), x1, circle2->getId(), x2);

  if (y1 != y2 && num_circles == 2)
    log_printf(ERROR, "Unable to ringify Cell %d since it contains "
               "Circle %d centered at y=%f and Circle %d at y=%f. "
               "Both Circles must have the same center.",
               _id, circle1->getId(), y1, circle2->getId(), y2);

  if (circle1 == NULL && circle2 != NULL)
    log_printf(ERROR, "Unable to ringify Cell %d since it only contains "
               "the positive halfpsace of Circle %d. Rings can only be "
               "created for Cells on the interior (negative halfspace) "
               "of a CIRCLE Surface.", _id, circle2->getId());

  if (radius1 <= radius2)
    log_printf(ERROR, "Unable to ringify Cell %d since it contains 2 "
               "disjoint CIRCLE Surfaces: halfspace %d for Circle %d "
               "and halfspace %d for Circle %d. Switch the signs of "
               "the 2 halfspaces for each Surface.", _id, halfspace1,
               circle1->getId(), halfspace2, circle2->getId());

  /* Compute the area to fill with each equal volume ring */
  double area = M_PI * fabs(radius1*radius1 - radius2*radius2) / _num_rings;

  /* Generate successively smaller Circle Surfaces */
  for (int i=0; i < _num_rings-1; i++) {
    radius2 = sqrt(radius1*radius1 - (area / M_PI));
    Circle* circle = new Circle(x1, y1, radius1);
    circles.push_back(circle);
    radius1 = radius2;
  }

  /* Store smallest, innermost Circle */
  Circle* circle = new Circle(x1, y1, radius1);
  circles.push_back(circle);

  /* Loop over Circles and create a new Cell clone for each ring */
  std::vector<Circle*>::iterator iter2;
  std::vector<CellBasic*>::iterator iter3;

  for (iter2 = circles.begin(); iter2 != circles.end(); ++iter2) {

    /* Create Circles for each of the sectorized Cells */
    if (_sectors.size() != 0) {
      for (iter3 = _sectors.begin(); iter3 != _sectors.end(); ++iter3) {
        log_printf(DEBUG, "Creating a new ring in sector Cell %d",
                   (*iter3)->getId());

        /* Create a new CellBasic clone */
        CellBasic* ring = (*iter3)->clone();
        ring->setNumSectors(0);
        ring->setNumRings(0);

        /* Add new bounding Circle surfaces to the clone */
        ring->addSurface(-1, (*iter2));

        /* Look ahead and check if we have an inner Circle to add */
        if (iter2+1 == circles.end()) {
          _rings.push_back(ring);
          continue;
        }
        else
          ring->addSurface(+1, *(iter2+1));

        /* Store the clone in the parent Cell's container of ring Cells */
        _rings.push_back(ring);
      }
    }

    /* Create Circles for this un-sectorized Cell */
    else {
      log_printf(DEBUG, "Creating new ring in un-sectorized Cell %d",_id);

      /* Create a new CellBasic clone */
      CellBasic* ring = clone();
      ring->setNumSectors(0);
      ring->setNumRings(0);

      /* Add new bounding Circle Surfaces to the clone */
      ring->addSurface(-1, (*iter2));

      /* Look ahead and check if we have an inner Circle to add */
      if (iter2+1 == circles.end()) {
        _rings.push_back(ring);
        break;
      }
      else
        ring->addSurface(+1, *(iter2+1));

      /* Store the clone in the parent Cell's container of ring Cells */
      _rings.push_back(ring);
    }
  }

  /* Store all of the rings in the parent Cell's subcells container */
  _subcells.clear();
  _subcells.insert(_subcells.end(), _rings.begin(), _rings.end());
}


/**
 * @brief Subdivides a Cell into rings and sectors.
 * @details This method uses the Cell's clone method to produce a vector of
 *          this Cell's subdivided ring and sector Cells.
 * @return a vector of Cell pointers to the new subdivided Cells
 */
std::vector<CellBasic*> CellBasic::subdivideCell() {
  sectorize();
  ringify();
  return _subcells;
}


/**
 * @brief Computes whether or not the Material filling this Cell is fissionable.
 */
void CellBasic::computeFissionability() {
  _fissionable = _material->isFissionable();
}

/**
 * @brief Convert this CellBasic's attributes to a string format.
 * @return a character array of this CellBasic's attributes
 */
std::string CellBasic::toString() {

  std::stringstream string;

  string << "Cell ID = " << _id
         << ", name = " << _name
         << ", type = MATERIAL"
         << ", material id = " << _material->getId()
         << ", # surfaces = " << getNumSurfaces()
         << ", # rings = " << _num_rings
         << ", # sectors = " << _num_sectors;

  /* Append each of the surface ids to the string */
  std::map<Surface*, int>::iterator iter;
  string << ", surface ids = ";
  for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
    string << iter->first->getId() << ", ";

  return string.str();
}


/**
 * @brief Prints a string representation of all of the CellBasic's attributes
 *        to the console.
 */
void CellBasic::printString() {
  log_printf(NORMAL, toString().c_str());
}


/**
 * @brief CellFill constructor
 * @param id the user-specified optional Cell ID
 * @param name the user-specified optional Cell name
 */
CellFill::CellFill(int id, const char* name): Cell(id, name) {
  _cell_type = FILL;
}


/**
 * @brief Return a pointer to the Universe filling this Cell.
 * @return the Universe pointer
 */
Universe* CellFill::getFill() const {
  return _fill;
}


/**
 * @brief Return the number of flat source regions in this CellFill.
 * @details This method is used when the geometry recursively constructs flat
 *          source regions.
 * @return the number of FSRs in this CellFill
 */
int CellFill::getNumFSRs() {

  if (_fill != NULL && _fill->getType() == SIMPLE)
    return _fill->computeFSRMaps();

  else
    return static_cast<Lattice*>(_fill)->computeFSRMaps();
}


/**
 * @brief Returns the std::map of Cell IDs and Cell pointers within any
 *        nested Universes filling this Cell.
 * @return std::map of Cell IDs and pointers
 */
std::map<int, Cell*> CellFill::getAllCells() {

  std::map<int, Cell*> cells;

  if (_cell_type == FILL && _fill != NULL) {
    std::map<int, Cell*> nested_cells = _fill->getAllCells();
    cells.insert(nested_cells.begin(), nested_cells.end());
  }

  return cells;
}


/**
 * @brief Returns the std::map of all nested Universe IDs and Universe pointers
          filling this Cell.
 * @return std::map of Universe IDs and pointers
 */
std::map<int, Universe*> CellFill::getAllUniverses() {

  std::map<int, Universe*> universes;

  if (_cell_type == FILL && _fill != NULL) {
    universes[_fill->getId()] = _fill;
    std::map<int, Universe*> nested_universes = _fill->getAllUniverses();
    universes.insert(nested_universes.begin(), nested_universes.end());
  }

  return universes;
}


/**
 * @brief Set a pointer to the Universe filling this CellFill.
 * @param universe the Universe's pointer
 */
void CellFill::setFill(Universe* universe) {
  _fill = universe;
}


/**
 * @brief Computes whether or not the Universe filling this Cell is fissionable.
 * @return whether or not the Universe filling the Cell is fissionable
 */
void CellFill::computeFissionability() {
  _fill->computeFissionability();
  _fissionable = _fill->isFissionable();
}


/**
 * @brief Convert this CellFill's attributes to a string format.
 * @return a character array of this Cell's attributes
 */
std::string CellFill::toString() {

  std::stringstream string;

  string << "Cell ID = " << _id
         << ", name = " << _name
         << ", type = FILL, "
         << ", fill = " << _fill->getId()
         << ", # surfaces = " << getNumSurfaces();

  /** Add the IDs for the Surfaces in this Cell */
  std::map<Surface*, int>::iterator iter;
  string << ", surface ids = ";
  for (iter = _surfaces.begin(); iter != _surfaces.end(); ++iter)
    string << iter->first->getId() << ", ";

  return string.str();
}


/**
 * @brief Prints a string representation of all of the CellFill's attributes
 *        to the console.
 */
void CellFill::printString() {
  log_printf(NORMAL, toString().c_str());
}
